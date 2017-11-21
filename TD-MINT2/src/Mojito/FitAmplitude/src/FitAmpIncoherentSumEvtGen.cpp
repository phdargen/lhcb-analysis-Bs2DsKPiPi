// author: Philippe d'Argent
#include "Mint/FitAmpIncoherentSumEvtGen.h"

#include "Mint/FitAmplitude.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/FitAmplitude.h"
#include "Mint/IntegCalculator.h"

#include <iostream>
#include <complex>

using namespace std;
using namespace MINT;

std::string FitAmpIncoherentSumEvtGen::IncPrefix(){
  return "Inco_";
}

FitAmpIncoherentSumEvtGen::FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		     , const char* fname
		     , MinuitParameterSet* pset
		     , const std::string& prefix
		     , const std::string& lineshapePrefix
		     , const std::string& opt
		     )
  : FitAmpList(pat, fname, pset, FitAmpIncoherentSumEvtGen::IncPrefix()+ prefix, lineshapePrefix, opt)
  , _useAnalyticGradient("useAnalyticGradient",0)
{
  /*
    //Important! Ensures everything is initialised
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(1,pat);
    this->getVal(eventTest[0]);  
  */
}

FitAmpIncoherentSumEvtGen::FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		     , MinuitParameterSet* pset
		     , const std::string& prefix
		     , const std::string& lineshapePrefix
		     , const std::string& opt
		     )
  : FitAmpList(pat, pset, FitAmpIncoherentSumEvtGen::IncPrefix() + prefix, lineshapePrefix, opt)
  , _useAnalyticGradient("useAnalyticGradient",0)
{
  /*
    //Important! Ensures everything is initialised
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(1,pat);
    this->getVal(eventTest[0]);  
  */
}
FitAmpIncoherentSumEvtGen::FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		     , const std::string& prefix
		     , const std::string& lineshapePrefix
		     , const std::string& opt
		     )
  : FitAmpList(pat, FitAmpIncoherentSumEvtGen::IncPrefix() + prefix, lineshapePrefix, opt)
  , _useAnalyticGradient("useAnalyticGradient",0)
{
  /*
    //Important! Ensures everything is initialised
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(1,pat);
    this->getVal(eventTest[0]);  
  */
}

FitAmpIncoherentSumEvtGen::FitAmpIncoherentSumEvtGen(const FitAmpIncoherentSumEvtGen& other)
  : IReturnRealForEvent<IDalitzEvent>()
  , IFastAmplitudeIntegrable()
  , ILookLikeFitAmpSum()
  , FitAmpList(other), _useAnalyticGradient("useAnalyticGradient",0)
{
  /*
    //Important! Ensures everything is initialised
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(1,_pat);
    this->getVal(eventTest[0]);  
  */
}

FitAmpIncoherentSumEvtGen::FitAmpIncoherentSumEvtGen(const FitAmpList& other)
  : IReturnRealForEvent<IDalitzEvent>()
  , IFastAmplitudeIntegrable()
  , ILookLikeFitAmpSum()
  , FitAmpList(other), _useAnalyticGradient("useAnalyticGradient",0)
{
  /*
    //Important! Ensures everything is initialised
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(1,_pat);
    this->getVal(eventTest[0]);
  */
}
counted_ptr<FitAmpListBase> FitAmpIncoherentSumEvtGen::GetCloneSameFitParameters() const{ 
// need to reform these one day...
// ... it all relies on the copy-constructur/=operator in FitAmpitude
// not copying the fit parameters, but just their pointers
// which will need to be reviewed.
//
  bool dbThis=false;
  if(dbThis) cout << "FitAmpSum::GetCloneSameFitParameters()" << endl;
  /* 
     There'll be 'physical' copies of all Amplitudes, but the
     FitParameters remain the same (pointers to the same
     FitParameter Object).  This is useful for the CP-con coding
     as it is now, but perhaps a bit counter-intuitive.  needs to
     be reviewed at some point. This behaviour is defined in the
     copy constructor of the FitAmplitude class.
  */

  /*
  counted_ptr<FitAmpIncoherentSumEvtGen> 
    newList(new FitAmpIncoherentSumEvtGen((IDalitzEventList*) this->getEventRecord()
				    , _paraFName.c_str(), _minuitParaSet));
  */
  counted_ptr<FitAmpIncoherentSumEvtGen> newList(new FitAmpIncoherentSumEvtGen(*this));
  newList->deleteAll();

  newList->add(*this);
  if(dbThis) cout << "cloning FitAmpIncoherentSumEvtGen " << newList->size() << endl;
  return newList;
}

FitAmpIncoherentSumEvtGen& FitAmpIncoherentSumEvtGen::operator=(const FitAmpIncoherentSumEvtGen& other){
  if(&other == this) return *this;
  (FitAmpList)(*this) = (FitAmpList) (other);
  return *this;
}
FitAmpIncoherentSumEvtGen& FitAmpIncoherentSumEvtGen::operator=(const FitAmpList& other){
  if(&other == this) return *this;
  (FitAmpIncoherentSumEvtGen)(*this) = other;
  return *this;
}

double FitAmpIncoherentSumEvtGen::getVal(IDalitzEvent& evt){
  bool dbthis=false;

  double sum(0);

  if(0 == this->size()){
    createAllAmps(evt.eventPattern()
		  , FitAmpIncoherentSumEvtGen::IncPrefix());
  }

  for(unsigned int i=0; i< this->size(); i++){
    if(dbthis){
      cout << "FitAmpIncoherentSumEvtGen::getVal()"
	   << "\n     > for " << getAmpPtr(i)->theBareDecay().oneLiner()
	   << "\n     > I get " << getAmpPtr(i)->getVal(evt)
	   << endl;
    }
    evt.setPermutationIndex(0);
    sum += norm(this->getAmpPtr(i)->getNewOnePermutationsVal(evt));
    evt.setPermutationIndex(1);
    sum += norm(this->getAmpPtr(i)->getNewOnePermutationsVal(evt));
  }

  if(dbthis) cout << "FitAmpIncoherentSumEvtGen::getVal(evt):"
		  << " returning this: " << sum 
		  << endl;

  if(false && sum > 200){
    cout << "large FitAmpIncoherentSumEvtGen " << sum
	 << " the largest amplitude is: "
	 << endl;
    printLargestAmp();
  }

  return efficiency(evt)*sum;

}

void FitAmpIncoherentSumEvtGen::Gradient(IDalitzEvent& evt,vector<double>& grad,MinuitParameterSet* mps){
    
    for (unsigned int i=0; i<mps->size(); i++) {
        if(mps->getParPtr(i)->hidden())continue;

        string name_i= mps->getParPtr(i)->name();
        if(name_i.find("Inco")==std::string::npos)continue;
        if(name_i.find("_Re")!=std::string::npos){
	  if(mps->getParPtr(i)->iFixInit() && mps->getParPtr(i+1)->iFixInit()){
	    i++;
	    continue;
	  }
	  //name_i.replace(name_i.find("Inco_"),5,"");
	  name_i.replace(name_i.find("_Re"),3,"");
	  for(unsigned int j=0; j< this->size(); j++){
	    if(A_is_in_B(name_i, this->getAmpPtr(j)->name())){
	      if(i+1 >= grad.size()){
		cout << "WARNING in FitAmpIncoherentSumEvtGen::Gradient"
		     << " have to increase size of grad to avoid memory issues" << endl;
		grad.resize(i+2);
	      }
	      double tmp = 2.*std::norm(this->getAmpPtr(j)->getValWithoutFitParameters(evt));
	      grad[i]= tmp * this->getAmpPtr(j)->AmpPhase().real();
	      grad[i+1]= tmp * this->getAmpPtr(j)->AmpPhase().imag();
	      i++;
	      break;
	    }
	  }
	}
	else if(mps->getParPtr(i)->iFixInit())continue;
	else {
	  std::cout << "FitAmpIncoherentSumEvtGen::Gradient() called. Sorry, I don't know how to calculate the derivative with respect to the fit parameter " << mps->getParPtr(i)->name() << " ! Please implement me or set useAnalytic Gradient to 0 in your options file. I'll crash now. " << std::endl;
	  throw "crash";
        }
        
    }
    
} 

counted_ptr<IIntegrationCalculator> 
FitAmpIncoherentSumEvtGen::makeIntegrationCalculator(){
  return (counted_ptr<IIntegrationCalculator>) makeIntegCalculator();
}
counted_ptr<IntegCalculator> 
FitAmpIncoherentSumEvtGen::makeIntegCalculator(){
  counted_ptr<IntegCalculator> l(new IntegCalculator);
  for(unsigned int i=0; i < _fitAmps.size(); i++){
    if(_fitAmps[i]->canBeIgnored()) continue;
    l->addAmps( (_fitAmps[i]),  (_fitAmps[i]));
  }

  for(unsigned int i=0; i < _fitAmpLists.size(); i++){
    l->append(_fitAmpLists[i]->makeIntegCalculator());
  }

  cout << "FitAmpIncoherentSumEvtGen: setting efficiency POINTER "
       << " in integCalculator to " 
       << _efficiency.get();
  if(0 == _efficiency.get()){
    cout << " (0 means no pointer, 100% efficiency).";
  }
  cout << endl;

  l->setEfficiency(_efficiency);
  return l;
}

counted_ptr<FitAmpPairList> 
FitAmpIncoherentSumEvtGen::makeFitAmpPairList(){
  counted_ptr<FitAmpPairList> l(new FitAmpPairList);
  for(unsigned int i=0; i < _fitAmps.size(); i++){
    if(_fitAmps[i]->canBeIgnored()) continue;
    l->addAmps( (_fitAmps[i]),  (_fitAmps[i]));
  }

  for(unsigned int i=0; i < _fitAmpLists.size(); i++){
    l->append(_fitAmpLists[i]->makeFitAmpPairList());
  }

  cout << "FitAmpIncoherentSumEvtGen: setting efficiency POINTER "
       << " in integCalculator to " 
       << _efficiency.get();
  if(0 == _efficiency.get()){
    cout << " (0 means no pointer, 100% efficiency).";
  }
  cout << endl;

  l->setEfficiency(_efficiency);
  return l;
}

void FitAmpIncoherentSumEvtGen::print(std::ostream& os) const{
   os << "FitAmpIncoherentSumEvtGen::print\n====================";

  for(unsigned int i=0; i< this->size(); i++){
    os << "\n\t" << this->getAmpPtr(i)->theBareDecay().oneLiner()
       << endl;
  }
}
void FitAmpIncoherentSumEvtGen::printNonZero(std::ostream& os) const{
   os << "FitAmpSum::print\n====================";

  for(unsigned int i=0; i < this->size(); i++){
    if(this->getAmpPtr(i)->isZero()) continue;
    os << "\n\t" << this->getAmpPtr(i)->theBareDecay().oneLiner()
       << endl;
  }
}


FitAmpIncoherentSumEvtGen::~FitAmpIncoherentSumEvtGen(){
  deleteAll();
}

FitAmpIncoherentSumEvtGen& 
FitAmpIncoherentSumEvtGen::operator+=(const FitAmpIncoherentSumEvtGen& other){
  add(other);
  return *this;
}
FitAmpIncoherentSumEvtGen 
FitAmpIncoherentSumEvtGen::operator+(const FitAmpIncoherentSumEvtGen& rhs) const{
  FitAmpIncoherentSumEvtGen fas(*this);
  fas.add(rhs);
  return fas;
}


FitAmpIncoherentSumEvtGen& FitAmpIncoherentSumEvtGen::operator*=(double r){
  multiply(r);
  return *this;
}
FitAmpIncoherentSumEvtGen& FitAmpIncoherentSumEvtGen::operator*=(const complex<double>& z){
  multiply(z);
  return *this;
}
FitAmpIncoherentSumEvtGen& FitAmpIncoherentSumEvtGen::operator*=(const counted_ptr<IReturnComplex>& irc){
  multiply(irc);
  return *this;
}

FitAmpIncoherentSumEvtGen FitAmpIncoherentSumEvtGen::operator*(double r) const{
  FitAmpIncoherentSumEvtGen fas(*this);
  fas.multiply(r);
  return fas;
}
FitAmpIncoherentSumEvtGen FitAmpIncoherentSumEvtGen::operator*(const complex<double>& z) const{
  FitAmpIncoherentSumEvtGen fas(*this);
  fas.multiply(z);
  return fas;
}
FitAmpIncoherentSumEvtGen FitAmpIncoherentSumEvtGen::operator*(const counted_ptr<IReturnComplex>& irc) const{
  FitAmpIncoherentSumEvtGen fas(*this);
  fas.multiply(irc);
  return fas;
}


FitAmpIncoherentSumEvtGen operator*(double r, const FitAmpIncoherentSumEvtGen& rhs){
  FitAmpIncoherentSumEvtGen fas(rhs);
  fas.multiply(r);
  return fas;
}
FitAmpIncoherentSumEvtGen operator*(const complex<double>& z, const FitAmpIncoherentSumEvtGen& rhs){
  FitAmpIncoherentSumEvtGen fas(rhs);
  fas.multiply(z);
  return fas;
}
FitAmpIncoherentSumEvtGen operator*(const counted_ptr<IReturnComplex>& irc
		     , const FitAmpIncoherentSumEvtGen& rhs){
  FitAmpIncoherentSumEvtGen fas(rhs);
  fas.multiply(irc);
  return fas;
}



//
