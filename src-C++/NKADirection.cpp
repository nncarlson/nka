#include "NKADirection.H"
#include "Teuchos_ENull.hpp"


NKADirection::NKADirection(const Teuchos::RCP<NOX::GlobalData>& gd, 
			   Teuchos::ParameterList& param, 
			   const NOX::Epetra::Vector &initvec)
{
  // get the parameters the parameter list
    
  Teuchos::ParameterList& p = param.sublist("NKADirection");
  
  int itmp = p.get("maxv", 10);
  int mvec  = itmp;
  
  double dtmp =  p.get("vtol", 1e-4);
  double vtol  = dtmp;
  
  state = new nka (mvec, vtol, initvec);
  
  reset(gd, param);
  
};
  

NKADirection::~NKADirection() 
{
  globalDataPtr = Teuchos::null;  

  delete state;
};


bool NKADirection::reset(const Teuchos::RCP<NOX::GlobalData> &gd, 
			 Teuchos::ParameterList &param)
{

  globalDataPtr = gd;
  utils = gd->getUtils();
  paramPtr = &param;
  
  state->nka_restart();
  
  return true;
  
};

  
bool NKADirection::compute (NOX::Abstract::Vector &dir, 
			    NOX::Abstract::Group &soln, 
			    const NOX::Solver::Generic  &solver)
{

  NOX::Abstract::Group::ReturnType status;
  
  // evaluate the nonlinear functional
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to compute F");

  NOX::Epetra::Vector f (dynamic_cast<const NOX::Epetra::Vector&>(soln.getF()), NOX::DeepCopy);  
  NOX::Epetra::Vector precond_f (f, NOX::ShapeCopy);

  status = soln.applyRightPreconditioning(false,paramPtr->sublist("NKADirection").sublist("Linear Solver"), f, precond_f);
  if (status != NOX::Abstract::Group::Ok) 
    throwError("compute", "Unable to apply preconditioner"); 
  
  NOX::Epetra::Vector d(precond_f);

  state->nka_correction(d, precond_f);

  d.scale(-1.0);

  dir = d;

  return true;
};
  
bool NKADirection::compute (NOX::Abstract::Vector  &dir, 
			    NOX::Abstract::Group  &soln, 
			    const NOX::Solver::LineSearchBased  &solver)
{
  
  return NOX::Direction::Generic::compute( dir, soln, solver );
};

  
void NKADirection::throwError(const string& functionName, 
			      const string& errorMsg)
{
  if (utils->isPrintType(NOX::Utils::Error))
    utils->err() << "MyDirection::" << functionName 
		 << " - " << errorMsg << endl;
  throw "NOX Error";
}
