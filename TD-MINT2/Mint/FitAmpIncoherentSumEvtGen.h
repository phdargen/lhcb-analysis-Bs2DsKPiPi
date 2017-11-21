#ifndef FitAmpIncoherentSumEvtGen_HH
#define FitAmpIncoherentSumEvtGen_HH
// author: Philippe d'Argent


#include <iostream>
#include "Mint/counted_ptr.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/IntegCalculator.h"
#include "Mint/IReturnRealForEvent.h"
#include "Mint/IReturnComplexForEvent.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/IIntegrationCalculator.h"
#include "Mint/ILookLikeFitAmpSum.h"
#include "Mint/FitAmpList.h"
#include <vector>

class FitAmpIncoherentSumEvtGen 
: virtual public MINT::IReturnRealForEvent<IDalitzEvent>
, virtual public IFastAmplitudeIntegrable
, virtual public ILookLikeFitAmpSum
, public FitAmpList
{
 protected:
  MINT::NamedParameter<int> _useAnalyticGradient;     
  static std::string IncPrefix();
 public:
  FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		      , const char* fname=0
		      , MINT::MinuitParameterSet* pset=0
		      , const std::string& prefix=""
		      , const std::string& lineshapePrefix=""
		      , const std::string& opt=""
		      );

  FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		      , MINT::MinuitParameterSet* pset
		      , const std::string& prefix=""
		      , const std::string& lineshapePrefix=""
		      , const std::string& opt=""
		      );
  FitAmpIncoherentSumEvtGen(const DalitzEventPattern& pat
		      , const std::string& prefix
		      , const std::string& lineshapePrefix=""
		      , const std::string& opt=""
		      );
  
  FitAmpIncoherentSumEvtGen(const FitAmpIncoherentSumEvtGen& other);
  FitAmpIncoherentSumEvtGen(const FitAmpList& other);
  /* 
     The copy constructor copies like this: There'll be 'physical'
     copies of all Amplitudes, but the FitParameters remain the
     same (pointers to the same FitParameter Object).  This is
     useful for the CP-conj coding as it is now, but perhaps a bit
     counter-intuitive.  needs to be reviewed at some point. This
     behaviour is defined in the copy constructor of the
     FitAmplitude class.
  */
  virtual MINT::counted_ptr<FitAmpListBase> GetCloneSameFitParameters() const;


  virtual DalitzBoxSet makeBoxes(const DalitzEventPattern& pat
				 , double nSigma=2){
    return FitAmpList::makeBoxes(pat, this, nSigma);}

  virtual DalitzBWBoxSet makeBWBoxes(const DalitzEventPattern& pat
				     , TRandom* rnd=gRandom){
    return FitAmpList::makeBWBoxes(pat, this, rnd);}


  double getVal(IDalitzEvent& evt);
  double getVal(IDalitzEvent* evtPtr){
    if(0 == evtPtr) return 0;
    return getVal(*evtPtr);
  }
  
  virtual void Gradient(IDalitzEvent& evt,std::vector<double>& grad, MINT::MinuitParameterSet* mps);
  virtual bool useAnalyticGradient() {return _useAnalyticGradient;}
  
  /*
  double getSmootherLargerVal();
  double getSmootherLargerVal(IDalitzEvent* evt);
  */

  virtual MINT::counted_ptr<IIntegrationCalculator> makeIntegrationCalculator();
  virtual MINT::counted_ptr<IntegCalculator> makeIntegCalculator();
  virtual MINT::counted_ptr<FitAmpPairList> makeFitAmpPairList();

  virtual double Prob(IDalitzEvent& evt){
    return getVal(evt);
  }

  /*
  virtual double SmootherLargerProb(){
    return getSmootherLargerVal();
  }
  */

  virtual double RealVal(IDalitzEvent& evt){
    return Prob(evt);
  }

  /*
  virtual double SmootherLargerRealVal(){
    return SmootherLargerProb();
  }
  */

  virtual MINT::counted_ptr<MINT::IUnweightedEventGenerator<IDalitzEvent> > 
    makeEventGenerator(const DalitzEventPattern& pat, TRandom* rnd=gRandom){
    MINT::counted_ptr<MINT::IUnweightedEventGenerator<IDalitzEvent> > 
      ptr(new DalitzBWBoxSet(makeBWBoxes(pat, rnd)));
    return ptr;
  }

  virtual ~FitAmpIncoherentSumEvtGen();

  void printLargestAmp(std::ostream& os = std::cout);
  virtual void printLargestAmp(IDalitzEvent& evt, std::ostream& os = std::cout){
    FitAmpList::printLargestAmp(evt, os);}

  virtual void print(std::ostream& os=std::cout) const;
  virtual void printNonZero(std::ostream& os=std::cout) const;

  friend class FitAmplitude;

  FitAmpIncoherentSumEvtGen& operator*=(double r);
  FitAmpIncoherentSumEvtGen& operator*=(const std::complex<double>& z);
  FitAmpIncoherentSumEvtGen& operator*=(const MINT::counted_ptr<MINT::IReturnComplex>& irc);

  FitAmpIncoherentSumEvtGen operator*(double r) const;
  FitAmpIncoherentSumEvtGen operator*(const std::complex<double>& z) const;
  FitAmpIncoherentSumEvtGen operator*(const MINT::counted_ptr<MINT::IReturnComplex>& irc) const;


  FitAmpIncoherentSumEvtGen& operator=(const FitAmpIncoherentSumEvtGen& other);
  FitAmpIncoherentSumEvtGen& operator=(const FitAmpList& other);
  FitAmpIncoherentSumEvtGen& operator+=(const FitAmpIncoherentSumEvtGen& other);
  FitAmpIncoherentSumEvtGen operator+(const FitAmpIncoherentSumEvtGen& other) const;

};

FitAmpIncoherentSumEvtGen operator*(double r, const FitAmpIncoherentSumEvtGen& rhs);
FitAmpIncoherentSumEvtGen operator*(const std::complex<double>& z, const FitAmpIncoherentSumEvtGen& rhs);
FitAmpIncoherentSumEvtGen operator*(const MINT::counted_ptr<MINT::IReturnComplex>& irc
		     , const FitAmpIncoherentSumEvtGen& rhs);


#endif
//
