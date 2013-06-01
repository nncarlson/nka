
// Simple nonlinear PDE problem.
// This file shows how to solve the nonlinear problem
//
// -\Delta u - \lambda e^u = 0  in \Omega = (0,1) \times (0,1)
//                       u = 0  on \partial \Omega
//
// using NOX 
//
// derived from  packages/didasko/nox/ex2.cpp
//
// uses NOX and implements NKA or Broyden as a user supplied direction
// uses ML for preconditioning

// define either USENKA or USEBROYDEN
#define USENKA
// #define USEBROYDEN



#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Operator.h"
#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "EpetraExt_HDF5.h"
#include "NOX_Direction_Generic.H"
#include "NOX_Direction_Factory.H"
#include "NOX_Direction_UserDefinedFactory.H"
#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Utils.H"
#include "NOX_MeritFunction_Generic.H"
#include "NOX_GlobalData.H"
#include "NOX_MultiVector.H"
#include "NOX_Epetra_Vector.H"


#include "NKA.H"
#include "NKADirection.H"
#include "NKADirFactory.H"

#include "BroydenSMDirection.hpp"
#include "BroydenSMDirFactory.hpp"

#include "Teuchos_ENull.hpp"
#include <ml_MultiLevelPreconditioner.h>

// this is required to know the number of lower, upper, left and right
// node for each node of the Cartesian grid (composed by nx \timex ny 
// elements)

static void  get_neighbours( const int i, const int nx, const int ny,
			     int & left, int & right, 
			     int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}

// This function creates a CrsMatrix, whose elements corresponds
// to the discretization of a Laplacian over a Cartesian grid,
// with nx grid point along the x-axis and and ny grid points 
// along the y-axis. For the sake of simplicity, I suppose that
// all the nodes in the matrix are internal nodes (Dirichlet
// boundary nodes are supposed to have been already condensated)

Epetra_CrsMatrix * CreateLaplacian( const int nx, const int ny,
				    const Epetra_Comm * Comm)
{

  int NumGlobalElements = nx * ny;
    
  // create a map
  Epetra_Map * Map = new Epetra_Map(NumGlobalElements,0,*Comm);
  // local number of rows
  int NumMyElements = Map->NumMyElements();
  // get update list
  int * MyGlobalElements = Map->MyGlobalElements();

  double hx = 1.0/(nx-1);
  double hy = 1.0/(ny-1);
  double off_left  = -1.0/(hx*hx);
  double off_right = -1.0/(hx*hx);
  double off_lower = -1.0/(hy*hy);
  double off_upper = -1.0/(hy*hy);
  double diag      =  2.0/(hx*hx) + 2.0/(hy*hy);
  
  int left, right, lower, upper;
    
  // a bit overestimated the nonzero per row
  
  Epetra_CrsMatrix * A = new Epetra_CrsMatrix(Copy,*Map,5);
    
  // Add  rows one-at-a-time
    
  double *Values = new double[4];
  int *Indices = new int[4];
    
  for( int i=0 ; i<NumMyElements; ++i ) {
    int NumEntries=0;
    get_neighbours(  MyGlobalElements[i], nx, ny, 
		     left, right, lower, upper);
    if( left != -1 ) {
      Indices[NumEntries] = left;
      Values[NumEntries] = off_left;
      ++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = off_right;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = off_lower;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = off_upper;
      ++NumEntries;
    }
    // put the off-diagonal entries
    A->InsertGlobalValues(MyGlobalElements[i], NumEntries, Values, Indices);
    // Put in the diagonal entry
    A->InsertGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
  }

  // put matrix in local ordering
  A->FillComplete();

  delete [] Indices;
  delete [] Values;
  delete    Map;

  return A;
  
} /* createJacobian */

// ==========================================================================
// This class contians the main definition of the nonlinear problem at
// hand. A method is provided to compute F(x) for a given x, and another
// method to update the entries of the Jacobian matrix, for a given x.
// As the Jacobian matrix J can be written as
//    J = L - diag(lambda*exp(x[i])),
// where L corresponds to the discretization of a Laplacian, and diag
// is a diagonal matrix with lambda*exp(x[i]). Basically, to update
// the jacobian we simply update the diagonal entries. Similarly, to compute
// F(x), we reset J to be equal to L, then we multiply it by the
// (distributed) vector x, then we add the diagonal contribution
// ==========================================================================

class PDEProblem {

public:

  // constructor. Requires the number of nodes along the x-axis
  // and y-axis, the value of lambda, and the Epetra_Communicator
  // (to define a Map, which is a linear map in this case)
  PDEProblem(const int nx, const int ny, const double lambda,
	     const Epetra_Comm * Comm) :
    nx_(nx), ny_(ny), lambda_(lambda) 
  {
    hx_ = 1.0/(nx_-1);
    hy_ = 1.0/(ny_-1);
    Matrix_ = CreateLaplacian(nx_,ny_,Comm);
  }

  // destructor
  ~PDEProblem() 
  {
    delete Matrix_;
  }

  // compute F(x)
  void ComputeF(const Epetra_Vector & x, Epetra_Vector & f) 
  {
    // reset diagonal entries
    double diag      =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);
  
    int NumMyElements = Matrix_->Map().NumMyElements();
    // get update list
    int * MyGlobalElements = Matrix_->Map().MyGlobalElements( );
    
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry
      Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, &diag, MyGlobalElements+i);
    }
    // matrix-vector product (intra-processes communication occurs
    // in this call)
    Matrix_->Multiply(false,x,f);

    // add diagonal contributions
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry

      f[i] += lambda_*exp(x[i]);   
    }

  }

  // update the Jacobian matrix for a given x
  void UpdateJacobian(const Epetra_Vector & x) 
  {
    double diag      =  2.0/(hx_*hx_) + 2.0/(hy_*hy_);
  
    int NumMyElements = Matrix_->Map().NumMyElements();
    // get update list
    int * MyGlobalElements = Matrix_->Map().MyGlobalElements( );
  
    for( int i=0 ; i<NumMyElements; ++i ) {
      // Put in the diagonal entry

      double newdiag = diag +lambda_*exp(x[i]); 

      Matrix_->ReplaceGlobalValues(MyGlobalElements[i], 1, 
				   &newdiag, MyGlobalElements+i);
    }

  }

  // returns a pointer to the internally stored matrix
  Epetra_CrsMatrix * GetMatrix()
  {
    return Matrix_;
  }
  
private:
  
  int nx_, ny_;
  double hx_, hy_;
  Epetra_CrsMatrix * Matrix_;
  double lambda_;

}; /* class PDEProblem */
 
// ==========================================================================
// This is the main NOX class for this example. Here we define
// the interface between the nonlinear problem at hand, and NOX.
// The constructor accepts a PDEProblem object. Using a pointer
// to this object, we can update the Jacobian and compute F(x),
// using the definition of our problem. This interface is bit
// crude: For instance, no PrecMatrix nor Preconditioner is specified.
// ==========================================================================

class SimpleProblemInterface : public NOX::Epetra::Interface::Required,
                               public NOX::Epetra::Interface::Jacobian,
			       public NOX::Epetra::Interface::Preconditioner

{

public:
 
  //! Constructor
  SimpleProblemInterface( PDEProblem * Problem ) :
    Problem_(Problem) {};

  //! Destructor
  ~SimpleProblemInterface() 
  {
  }

  bool computeF(const Epetra_Vector & x, Epetra_Vector & f,
                NOX::Epetra::Interface::Required::FillType F )
  {
    Problem_->ComputeF(x,f);

    return true;
  }
  
  bool computeJacobian(const Epetra_Vector & x, Epetra_Operator & Jac)
  {
    Problem_->UpdateJacobian(x);
    return true;
  }

  bool computePrecMatrix(const Epetra_Vector & x, Epetra_RowMatrix & M) 
  {
    cout << "*ERR* SimpleProblem::preconditionVector()\n";
    cout << "*ERR* don't use explicit preconditioning" << endl;
    throw 1;
  }  
  
  bool computePreconditioner(const Epetra_Vector & x, Epetra_Operator & O,
		 Teuchos::ParameterList* pl)
  {
    return true;
  }  

private:
  
  PDEProblem * Problem_;
  
}; /* class SimpleProblemInterface */




// =========== //
// main driver //
// =========== //

int main( int argc, char **argv )
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // define the parameters of the nonlinear PDE problem
  int nx = 500;
  int ny = 500;  
  double lambda = 1.0;

  // create the problem
  PDEProblem Problem(nx,ny,lambda,&Comm);
 
  // starting solution, here a zero vector
  Epetra_Vector InitialGuess(Problem.GetMatrix()->Map());
  InitialGuess.PutScalar(0.0);

  // Set up the problem interface
  Teuchos::RCP<SimpleProblemInterface> interface = 
    Teuchos::rcp(new SimpleProblemInterface(&Problem) );
  
  // Create the top level parameter list
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr =
    Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList& nlParams = *(nlParamsPtr.get());

  // Set the nonlinear solver method
  nlParams.set("Nonlinear Solver", "Line Search Based");
  

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm.MyPID()); 
  printParams.set("Output Precision", 3);
  printParams.set("Output Processor", 0);
  printParams.set("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);

  // start definition of nonlinear solver parameters
  // Sublist for line search 
  Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  searchParams.set("Method","Full Step");

  // Sublist for direction
  // here we define our own direction method
  Teuchos::ParameterList &dirParams = nlParams.sublist("Direction");
  dirParams.set("Method", "User Defined");


#if defined(USENKA)

  // these are parameters that the NKA direction class
  // will look at
  Teuchos::ParameterList &mydirParams = nlParams.sublist("NKADirection");
  mydirParams.set("maxv", 10);
  mydirParams.set("vtol", 1e-8);
  mydirParams.set("beta", 1.0);

#elif defined(USEBROYDEN)

  // these are parameters that the BroydenSM direction class
  // will look at
  Teuchos::ParameterList &mydirParams = nlParams.sublist("BroydenSMDirection");
  mydirParams.set("nmax", 10);
  mydirParams.set("precondition", true);

#else
#error "Either USENKA or USEBROYDEN must be defined"
#endif


  // we need a GlobalData variable fo initialize the 
  // Direction factory
  Teuchos::RCP<NOX::GlobalData> gd = Teuchos::rcp(new NOX::GlobalData(nlParamsPtr));
  // Need a NOX::Epetra::Vector for constructors
  NOX::Epetra::Vector noxInitGuess(InitialGuess, NOX::DeepCopy);  


#if defined(USENKA)

  // create our direction factory
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> dir_factory 
    = Teuchos::rcp(new NKADirFactory(gd, nlParams, noxInitGuess));
  dirParams.set("User Defined Direction Factory", dir_factory);

#elif defined(USEBROYDEN)

  // create our direction factory
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> dir_factory 
    = Teuchos::rcp(new BroydenSMDirFactory(gd, nlParams, noxInitGuess));
  dirParams.set("User Defined Direction Factory", dir_factory);  

#else 
#error "Either USENKA or USEBROYDEN must be defined"
#endif

  // we need a linear solver to associate a preconditioner to it
  // we will then only use the preconditioner through the applyRightPreconditioning
  // method
  Teuchos::ParameterList& lsParams = mydirParams.sublist("Linear Solver");
  lsParams.set("Preconditioner","User Defined");
    
  Teuchos::RCP<NOX::Epetra::Interface::Required> iReq = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;

  // create an ML preconditioner
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPrec = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Problem.GetMatrix(), true) );

  // create the linear solver with its ML preconditioner
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
    Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
  					       iReq, iPrec, MLPrec, InitialGuess)); 


  
  Teuchos::RCP<NOX::Epetra::Group> grpPtr = 
    Teuchos::rcp(new NOX::Epetra::Group(printParams, 
					iReq, 
					noxInitGuess, 
  					linSys)); 

  // Set up the status tests
  Teuchos::RCP<NOX::StatusTest::NormF> testNormF = 
    Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-6));
  Teuchos::RCP<NOX::StatusTest::MaxIters> testMaxIters = 
    Teuchos::rcp(new NOX::StatusTest::MaxIters(100));
  // this will be the convergence test to be used
  Teuchos::RCP<NOX::StatusTest::Combo> combo = 
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR, 
					    testNormF, testMaxIters));

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = 
    NOX::Solver::buildSolver(grpPtr, combo, nlParamsPtr);

  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();

if( Comm.MyPID() == 0 )   
  if( NOX::StatusTest::Converged  == status )
    cout << "\n" << "-- NOX solver converged --" << "\n";
  else
    cout << "\n" << "-- NOX solver did not converge --" << "\n";
 
  // Print the answer
  if( Comm.MyPID() == 0 ) {
  cout << "\n" << "-- Parameter List From Solver --" << "\n";
  solver->getList().print(cout);
  }

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group & finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector & finalSolution = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();


  gd = Teuchos::null;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
} /* main */

