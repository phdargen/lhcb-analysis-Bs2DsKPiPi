/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 * @(#)root/roofit:$Id: RooSplineProduct.cxx
 * Authors:                                                                  *
 *   Katya Govorkova, Nikhef 2k17
 *                                                                           *
 *****************************************************************************/

// STD & STL
#include <cmath>
#include <complex>
#include <sstream>

// ROOT
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

// RooFit
#include "RooMsgService.h"
#include "RooMath.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgList.h"

#include "Mint/RooSplineProduct.h"

using namespace std;

//_____________________________________________________________________________
void RooSplineProduct::init() {
    // TODO: current implementation of product of two splines
    // only for splines with the same knot vector
    assert(_sp1.knots() == _sp2.knots());
    assert(_sp1.knotSize() == _sp2.knotSize());
}

//_____________________________________________________________________________
RooSplineProduct::RooSplineProduct()
{
}

//_____________________________________________________________________________
RooSplineProduct::RooSplineProduct(const char* name,
                                   const char* title, RooRealVar& x,
                                   const RooCubicSplineFun& sp1,
                                   const RooCubicSplineFun& sp2
                                 ) :
   RooAbsGaussModelEfficiency(name, title),
   _x("x", "Dependent", this, x),
   _sp1(sp1,"spline1"),
   _sp2(sp2,"spline2"),
   _coefList1("coefficients1","List of coefficients for first spline",this),
   _coefList2("coefficients2","List of coefficients for second spline",this)
{
  _coefList1.add(sp1.coefficients());
  _coefList2.add(sp2.coefficients());
  init();
}

//_____________________________________________________________________________
RooSplineProduct::RooSplineProduct(const RooSplineProduct& other, const char* name) :
  RooAbsGaussModelEfficiency(other, name),
  _x("x", this, other._x),
  _sp1(other._sp1,"spline1"),
  _sp2(other._sp2,"spline2"),
  _coefList1("coefList1", this, other._coefList1),
  _coefList2("coefList2", this, other._coefList2)
{
}

//_____________________________________________________________________________
RooSplineProduct::~RooSplineProduct()
{
}

//_____________________________________________________________________________
Double_t RooSplineProduct::evaluate() const
{
  return _sp1._aux.evaluate(_x, _coefList1) * _sp2._aux.evaluate(_x, _coefList2);
}

//_____________________________________________________________________________
Int_t RooSplineProduct::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const
{
  return 0;
}
//_____________________________________________________________________________
Double_t RooSplineProduct::analyticalIntegral(Int_t code, const char* /* rangeName */) const
{
  return 0;
}

//_____________________________________________________________________________
complex<double>  RooSplineProduct::gaussIntegralE(bool left, const RooGaussModelAcceptance::M_n<7U>& dM,
                                             const RooGaussModelAcceptance::K_n& K, double offset,
                                             double* sc) const
{
  RooCubicSplineKnot::S2_edge s2_jk( _sp1._aux.S2_jk_edge( left, _coefList1, _coefList2 ), offset );
  std::complex<double> sum(0,0);
  for (int j=0;j<3;++j) for (int k=0;k<3-j;++k) sum += dM(j)*s2_jk(j,k)*K(k)*sc[j+k];
  return sum;
}

//_____________________________________________________________________________
std::complex<double>
RooSplineProduct::productAnalyticalIntegral(Double_t umin, Double_t umax,
                                             Double_t scale, Double_t offset,
                                             const std::complex<double>& z) const
{
    RooGaussModelAcceptance::K_n K(z);
    assert(_sp1.knotSize()>1);
    double lo = scale*umin+offset;
    double hi = scale*umax+offset+1.e-7;
    typedef RooGaussModelAcceptance::M_n<7U> M_n;
    std::vector<M_n> M; M.reserve( _sp1.knotSize() );
    for (unsigned int i=0;i<_sp1.knotSize();++i) {
        double x = (_sp1.u(i)-offset)/scale ;
        if (lo>=_sp1.u(i)) x = umin ; // take M[i] if lo<=u(i) else M_n(lo)
        if (_sp1.u(i)>=hi) x = umax ; // take M[i+1] if u(i+1)<=hi else M_n(hi)
        M.push_back( M_n( x, z ) );
    }
    double sc[7]; for (int i=0;i<7;++i) sc[i] = pow(scale,i);
    std::complex<double> sum(0,0);
    if (lo<_sp1.u(0)) sum += gaussIntegralE(true,  M.front()-M_n( umin,z), K, offset, sc);
    for (unsigned i=0;i<_sp1.knotSize()-1 && _sp1.u(i)<hi ;++i) {
        M_n dM = M[i+1]-M[i];
        RooCubicSplineKnot::S2_jk s2_jk( _sp1._aux.S2_jk_sum(i, _coefList1, _coefList2), offset );
        for (int j=0;j<7;++j) for (int k=0;k<7-j;++k) sum += dM(j)*s2_jk(j,k)*K(k)*sc[j+k];
    }
    if (hi>_sp1.u(_sp1.knotSize()-1)) sum += gaussIntegralE(false, M_n(umax,z)-M.back(), K, offset, sc);
    return sum;
}
