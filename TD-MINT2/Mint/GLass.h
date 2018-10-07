#ifndef GLASS_LINESHAPE_HH
#define GLASS_LINESHAPE_HH
// author: Philippe d'Argent

#include "Mint/ILineshape.h"
#include "Mint/BW_BW.h"
#include "Mint/AssociatedDecayTree.h"
#include "Mint/IDalitzEvent.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/NamedParameter.h"
#include "Mint/ResonancePropertiesList.h"
#include "Mint/FitParDependent.h"
#include "Mint/ResonancePropertiesFitRef.h"

#include <complex>
#include <string>

/*
  Lass parameterised a la arXiv:1004.5053v3
*/

class GLass : public BW_BW, virtual public ILineshape{
 protected:

 virtual std::complex<double> BreitWigner();

 double a() const{return mumsFittableProperties().fitGLASS_a()/GeV ;}
 double r() const{return mumsFittableProperties().fitGLASS_r()/GeV ;}
 double R() const{return mumsFittableProperties().fitGLASS_R() ;}
 double phiR() const{return mumsFittableProperties().fitGLASS_phiR() ;}
 double F() const{return mumsFittableProperties().fitGLASS_F() ;}
 double phiF() const{return mumsFittableProperties().fitGLASS_phiF() ;}
 double alpha1() const{return mumsFittableProperties().fitGLASS_alpha1() ;}
 double alpha2() const{return mumsFittableProperties().fitGLASS_alpha2() ;}
 double alpha3() const{return mumsFittableProperties().fitGLASS_alpha3() ;}

 public:
  
    GLass( const AssociatedDecayTree& tree, const std::string& namePrefix=""): 
    BW_BW(tree, namePrefix)
    {}

  virtual std::string name() const{
    return "GLASS("+prefix()+_theDecay.oneLiner() +")";
  }
  virtual ~GLass(){}
};

#endif
//
