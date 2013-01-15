#ifndef __BroydenSMDIRECTION_H__
#define __BroydenSMDIRECTION_H__

#include "NOX_Common.H"

#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Direction_Generic.H"
#include "NOX_GlobalData.H"
#include "NOX_Utils.H"

class BroydenSMDirection : public NOX::Direction::Generic {

public:

  BroydenSMDirection(const Teuchos::RCP<NOX::GlobalData>&, 
	       Teuchos::ParameterList&, const NOX::Abstract::Vector&);

  ~BroydenSMDirection();

  bool reset (const Teuchos::RCP<NOX::GlobalData>&, 
	      Teuchos::ParameterList&);
  
  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::Generic&);

  bool compute (NOX::Abstract::Vector&, NOX::Abstract::Group&,
		const NOX::Solver::LineSearchBased&);

private:
  
  void throwError(const string&, const string&);

  // Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  // Global data pointer.  
  // Keep this so any stored parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  Teuchos::ParameterList *paramPtr; 

  Teuchos::RCP<NOX::Abstract::Vector> tmpVecPtr;

  int nmax;     // maximum number of Broyden iterations before a restart
  int n;        // number of currenr Broyden vectors 
  int itc;      // count the nonlinear iterations
  bool precond; // whether to use preconditioning

  std::list<Teuchos::RCP<NOX::Abstract::Vector> > sptr;
  std::list<double> snormsq;

  const NOX::Abstract::Vector& vec;
  
};

#endif
