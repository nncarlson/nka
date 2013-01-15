#include "BroydenSMDirection.hpp"


BroydenSMDirection::BroydenSMDirection(const Teuchos::RCP<NOX::GlobalData>& gd, 
				       Teuchos::ParameterList& param, 
				       const NOX::Abstract::Vector& initvec):
  vec(initvec)

{
  reset(gd, param);  
};
  

BroydenSMDirection::~BroydenSMDirection() 
{

};


bool BroydenSMDirection::reset(const Teuchos::RCP<NOX::GlobalData> &gd, 
			       Teuchos::ParameterList &param)
{
  globalDataPtr = gd;
  utils = gd->getUtils();
  paramPtr = &param;
  
  nmax    = param.get<int>("nmax", 10);
  precond = param.get<bool>("precondition", true); 
  itc = -1;

  return true;
};

  
bool BroydenSMDirection::compute (NOX::Abstract::Vector &dir, 
			    NOX::Abstract::Group &soln, 
			    const NOX::Solver::Generic  &solver)
{
  NOX::Abstract::Group::ReturnType status;
  
  // evaluate the nonlinear functional
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok) {
    throwError("compute", "Unable to compute F");
  }

  Teuchos::RCP<NOX::Abstract::Vector> fptr 
    = soln.getF().clone(NOX::DeepCopy);

  Teuchos::RCP<NOX::Abstract::Vector> z 
    = (*fptr).clone(NOX::ShapeCopy);


  if (precond) 
    {
      status 
	= soln.applyRightPreconditioning(false,
					 paramPtr->sublist("Linear Solver"), 
					 *fptr, *z); 

      if (status != NOX::Abstract::Group::Ok) 
	{
	  throwError("compute", "Unable to apply preconditioner"); 
	}
    }
  else
    {
      *z = *fptr;
    }


  // compute the direction

  if (itc == -1) 
    {
      itc++;
      n = -1;
      sptr.push_back(z);
      sptr.front()->scale(-1.0);
      snormsq.push_back(pow(sptr.front()->norm(),2.0));
    } 
  else
    {
      itc++;
      n++;
      
      z->scale(-1.0);
      
      double sdotz;
      std::list<Teuchos::RCP<NOX::Abstract::Vector> >::iterator s;
      std::list<double>::iterator snsq = snormsq.begin();
    
      s = sptr.begin();
      for (int j=0; j<sptr.size()-1; j++)
	{
	  sdotz = (*s)->innerProduct(*z);
	  s++;
	  z->update(sdotz/(*snsq),*(*s),1.0);
	  snsq++;
	}
    
      sdotz = sptr.back()->innerProduct(*z);
      z->scale(1.0/(1.0-sdotz/snormsq.back()));
      sptr.push_back(z);
      snormsq.push_back(pow((*sptr.back()).norm(),2.0));



      if (sptr.size()>nmax) 
	{
	  sptr.erase(sptr.begin());
	  snormsq.erase(snormsq.begin());

	  n=nmax;
	}



    }

  dir = *sptr.back();

  return true;

};
  
bool BroydenSMDirection::compute (NOX::Abstract::Vector  &dir, 
			    NOX::Abstract::Group  &soln, 
			    const NOX::Solver::LineSearchBased  &solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
};

  
void BroydenSMDirection::throwError(const string& functionName, 
			      const string& errorMsg)
{
  if (utils->isPrintType(NOX::Utils::Error))
    utils->err() << "BroydenSMDirection::" << functionName 
		 << " - " << errorMsg << endl;
  throw "NOX Error";
}
