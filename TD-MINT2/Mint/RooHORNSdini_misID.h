#ifndef ROO_ROOHORNSDINI_MISID
#define ROO_ROOHORNSDINI_MISID

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TMath.h"


class RooRealVar;

class  RooHORNSdini_misID : public RooAbsPdf {


public:
  RooHORNSdini_misID(const char *name, const char *title,RooAbsReal& _m, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _csi, RooAbsReal& _m1, RooAbsReal& _s1, RooAbsReal& _m2, RooAbsReal& _s2, RooAbsReal& _m3, RooAbsReal& _s3, RooAbsReal& _m4, RooAbsReal& _s4, RooAbsReal& _f1, RooAbsReal& _f2,  RooAbsReal& _f3);

  RooHORNSdini_misID(const  RooHORNSdini_misID& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new  RooHORNSdini_misID(*this,newname); }
  inline virtual ~ RooHORNSdini_misID() {}


protected:
  RooRealProxy m;
  RooRealProxy a;
  RooRealProxy b;
  RooRealProxy csi;
  RooRealProxy m1;
  RooRealProxy s1;
  RooRealProxy m2;
  RooRealProxy s2;
  RooRealProxy m3;
  RooRealProxy s3;
  RooRealProxy m4;
  RooRealProxy s4;
  RooRealProxy f1;
  RooRealProxy f2;
  RooRealProxy f3;


  Double_t evaluate() const;


private:
//  ClassDef( RooHORNSdini_misID,0)
};

#endif
