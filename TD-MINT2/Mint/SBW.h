#ifndef SBW_LINESHAPE_HH
#define SBW_LINESHAPE_HH
// author: Philippe d'Argent (p.dargent@cern.ch)
// status: 19 March 2015 GMT

#include <complex>
#include <string>

#include "Mint/ILineshape.h"
#include "Mint/BW_BW.h"
#include "Mint/NamedParameter.h"
#include "Mint/FitParDependent.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/CLHEPPhysicalConstants.h"

class SBW : public BW_BW, virtual public ILineshape{
 public:
  
 SBW( const AssociatedDecayTree& tree, const std::string& namePrefix): BW_BW(tree, namePrefix){}

  virtual std::string name() const{
    return "SBW("+prefix()+_theDecay.oneLiner() +")";
  }

  virtual ~SBW(){}

 protected:

  virtual double GofM() {return mumsWidth();}
  virtual std::complex<double> BreitWigner(){
	  double mass = mumsMass();
  	  double width = mumsWidth();

          const double m2hh = mumsRecoMass2()/GeV/GeV;
	  //double p = twoBody_dgtPsq_in_MumsFrame(mumsRecoMass(), daughterPDGMass(0), daughterPDGMass(1));
	  //if(p <= 0) return 0.;

	  std::complex<double> invBW(mumsRecoMass()-mass, - width/2.);
  	  return 1.*GeV/invBW;//*pow(pABSq(),GetAlpha());//*(1.+pABSq()/GeV*c1()+pow(pABSq()/GeV,2.)*c2()+pow(pABSq()/GeV,3.)*c3());
  } 

//  double GetAlpha() const{
//	  return _RPL->get(mumsPID())->alpha();
//  }
  double c1() const {
        return _RPL->get(mumsPID())->c1();
  }
  double c2() const {
        return _RPL->get(mumsPID())->c2();
  }
  double c3() const {
        return _RPL->get(mumsPID())->c3();
  }
  double c4() const {
        return _RPL->get(mumsPID())->c4();
  }
  double c5() const {
        return _RPL->get(mumsPID())->c5();
  }

};

#endif
//
