
#include "BroydenSMDirFactory.hpp"


BroydenSMDirFactory::BroydenSMDirFactory(const Teuchos::RCP<NOX::GlobalData> &gd, 
			     Teuchos::ParameterList &params, 
			     const NOX::Abstract::Vector &initvec)

{
  my_dir = Teuchos::rcp(new BroydenSMDirection(gd, params, initvec));
};


BroydenSMDirFactory::~BroydenSMDirFactory() 
{
};


Teuchos::RCP<NOX::Direction::Generic> 
BroydenSMDirFactory::buildDirection(const Teuchos::RCP<NOX::GlobalData> &gd, 
			      Teuchos::ParameterList &params) const
{
  return my_dir;
};
