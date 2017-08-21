/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   Katya Govorkova, Nikhef 2k17
 *                                                                           *
 *****************************************************************************/
#ifndef ROO_SPLINEPRODUCT
#define ROO_SPLINEPRODUCT

#include "RooRealProxy.h"
#include "RooListProxy.h"

#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooAbsGaussModelEfficiency.h"

class RooRealVar;
class RooArgList ;
class TH1;
class TGraph;
class TGraphErrors;

class RooSplineProduct : public  RooAbsGaussModelEfficiency {
public:
  RooSplineProduct() ;
  RooSplineProduct(const char* name, const char* title, RooRealVar& x,
    const RooCubicSplineFun& sp1, const RooCubicSplineFun& sp2);
  ~RooSplineProduct() ;

  RooSplineProduct(const RooSplineProduct& other, const char* name = 0);
  TObject* clone(const char* newname) const { return new RooSplineProduct(*this, newname); }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const;
  Double_t analyticalIntegral(Int_t code, const char* rangeName) const;

  // for use as RooAbsGaussModelEfficiency...
  std::complex<double> productAnalyticalIntegral(Double_t umin, Double_t umax
                                                ,Double_t scale, Double_t offset
                                                ,const std::complex<double>& z) const;

private:

  RooRealProxy _x;
  RooCubicSplineFun _sp1;
  RooCubicSplineFun _sp2;
  RooListProxy _coefList1;
  RooListProxy _coefList2;

  void init();

  Double_t evaluate() const;

  // for use in RooGaussEfficiencyModel...
  std::complex<double> gaussIntegralE(bool left, const RooGaussModelAcceptance::M_n<7U>& dM,
                           const RooGaussModelAcceptance::K_n& K,
                           double offset, double* sc) const ;

  //ClassDef(RooSplineProduct,1) // Product of two cubic splines polynomial PDF
};

#endif
