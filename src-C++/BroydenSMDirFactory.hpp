#ifndef _BROYDENSMDIRFACTORY_H_
#define _BROYDENSMDIRFACTORY_H_


#include "NOX_Direction_UserDefinedFactory.H"
#include "BroydenSMDirection.hpp"



class BroydenSMDirFactory : public NOX::Direction::UserDefinedFactory {
public: 

  BroydenSMDirFactory(const Teuchos::RCP<NOX::GlobalData>&, Teuchos::ParameterList&, 
		const NOX::Abstract::Vector&);
  
  ~BroydenSMDirFactory(); 

  Teuchos::RCP<NOX::Direction::Generic> 
  buildDirection(const Teuchos::RCP<NOX::GlobalData>&, Teuchos::ParameterList&) const;

private:
  Teuchos::RCP<BroydenSMDirection>  my_dir;

};

#endif
