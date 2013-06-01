#ifndef __NKADIRFACTORY_H__
#define __NKADIRFACTORY_H__

#include "NOX_Direction_UserDefinedFactory.H"
#include "NKADirection.H"

class NKADirFactory : public NOX::Direction::UserDefinedFactory {
public: 

  NKADirFactory(const Teuchos::RCP<NOX::GlobalData>&, Teuchos::ParameterList&, 
		const NOX::Abstract::Vector&);
  
  ~NKADirFactory(); 

  Teuchos::RCP<NOX::Direction::Generic> 
  buildDirection(const Teuchos::RCP<NOX::GlobalData>&, Teuchos::ParameterList&) const;

private:
  Teuchos::RCP<NKADirection>  my_dir;

};

#endif
