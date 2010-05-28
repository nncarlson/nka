//   NONLINEAR_KRYLOV_ACCELERATOR - C++/Trilinos implementation

//   Neil N. Carlson <neil.n.carlson@gmail.com>
//   Markus Berndt <markus.berndt@gmail.com>

//   This code implements the nonlinear Krylov accelerator introduced in [1]
//   for inexact Newton's (IN) method, where the correction equation of
//   Newton's method is only approximately solved because the Jacobian matrix
//   is approximated and/or the linear system is not solved exactly.  Placed
//   in the iteration loop, this black-box accelerator listens to the sequence
//   of inexact corrections and replaces them with accelerated corrections;
//   the resulting method is a type of accelerated inexact Newton (AIN) method.
//   Note that an IN iteration is merely a standard fixed point iteration for
//   a preconditioned system, and so this accelerator is more generally
//   applicable to fixed point iterations.

//   This code is a straightforward translation of the original Fortran 95
//   implementation into C++/Trilinos.

//   [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
//       weighted moving finite element code I: in one dimension", SIAM J.
//       Sci. Comput;, 19 (1998), pp. 728-765.  See section 9.

//   ************************************************************************

//   Copyright (c) 2010  Neil N. Carlson, Markus Berndt

//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:

//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.

//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//   DEALINGS IN THE SOFTWARE.



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
