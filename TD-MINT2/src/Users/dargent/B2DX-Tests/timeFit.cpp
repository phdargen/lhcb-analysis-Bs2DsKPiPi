// Performs simple toy fit with Mixing PDF 
// author: Philippe d'Argent
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/DalitzPdfBase.h"
#include "Mint/DalitzPdfBaseFastInteg.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/AmpRatios.h"
#include "Mint/IEventGenerator.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/SignalGenerator.h"
#include "Mint/FromFileGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/cexp.h"
#include "Mint/DalitzPdfNormChecker.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/DalitzPdfSaveInteg.h"
#include "Mint/Chi2Binning.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpList.h"
#include "Mint/DalitzPdfBaseMCInteg.h"
#include "RooRealConstant.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDecay.h"
#include "RooBDecay.h"
#include "RooPlot.h"
#include "RooEffProd.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using namespace RooFit ;
//using namespace RooStats;
using namespace MINT;

// DsK like time PDF with additional coherence factor 
class TimePdf : public MINT::PdfBase<IDalitzEvent>
{
protected:
    FitParameter& _C;
    FitParameter& _D;
    FitParameter& _D_bar;
    FitParameter& _S;
    FitParameter& _S_bar;
    FitParameter& _k;

    FitParameter& _tau;
    FitParameter& _dGamma;
    FitParameter& _dm;
    FitParameter& _eff_tag;
    FitParameter& _w;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
    }
    void endFit(){
    }

    inline double un_normalised(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        const double D = 1./2. * ((1.+f) * _D + (1.-f) * _D_bar);
        const double S = 1./2. * ((1.+f) * _S + (1.-f) * _S_bar);

        const double val =  exp(-fabs(t)/(double)_tau) *
        (
         (2.-fabs(q))*cosh((double)_dGamma/2.*t)
         +f*q*(1.-2.*w)* _C *cos((double)_dm*t)
         -(2.-fabs(q))*2.0*_k* D *sinh((double)_dGamma/2.*t)
         -f*2.0*q*(1.-2.*w)*_k* S *sin((double)_dm*t)
         )*e_eff;
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double D = 1./2. * ((1.+f) * _D + (1.-f) * _D_bar);

        const double Gamma = 1./((double) _tau);

        double val = un_normalised(evt);

        double norm = 2.* ((4.*Gamma/(4.*Gamma*Gamma-_dGamma*_dGamma)) - 2.0*_k* D * (2.*_dGamma/(4.*Gamma*Gamma-_dGamma*_dGamma)));
        
        return val/norm;
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }
    
    TimePdf(MINT::FitParameter& C,MINT::FitParameter& D, MINT::FitParameter& D_bar,MINT::FitParameter& S, MINT::FitParameter& S_bar,MINT::FitParameter& k, MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm, MINT::FitParameter& eff_tag,MINT::FitParameter& w ):
    _C(C),_D(D),_D_bar(D_bar),_S(S),_S_bar(S_bar),_k(k),_tau(tau),_dGamma(dGamma),_dm(dm),_eff_tag(eff_tag), _w(w) {;}
    
};


// Time PDF in term of r, gamma, delta instead of CP coefficients
class TimePdf_mod : public MINT::PdfBase<IDalitzEvent>
{
protected:
    FitParameter& _r;
    FitParameter& _delta;
    FitParameter& _gamma;
    FitParameter& _k;
    
    FitParameter& _tau;
    FitParameter& _dGamma;
    FitParameter& _dm;
    FitParameter& _eff_tag;
    FitParameter& _w;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2));
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
    
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        const double C = (1.-(double)_r*(double)_r)/(1.+(double)_r*(double)_r);
        const double D = (double)_r * cos(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
        const double S = (double)_r * sin(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
        
        const double val =  exp(-fabs(t)/(double)_tau) *
        (
         (2.-fabs(q))*cosh((double)_dGamma/2.*t)
         +f*q*(1.-2.*w)* C *cos((double)_dm*t)
         -(2.-fabs(q))*2.0*_k* D *sinh((double)_dGamma/2.*t)
         -f*2.0*q*(1.-2.*w)*_k* S *sin((double)_dm*t)
         )*e_eff;
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double D = (double)_r * cos(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
        
        const double Gamma = 1./((double) _tau);
        
        double val = un_normalised(evt);
        
        double norm = 2.* ((4.*Gamma/(4.*Gamma*Gamma-_dGamma*_dGamma)) - 2.0*_k* D * (2.*_dGamma/(4.*Gamma*Gamma-_dGamma*_dGamma)));
        
        return val/norm;
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }
    
    TimePdf_mod(MINT::FitParameter& r,MINT::FitParameter& delta, MINT::FitParameter& gamma,MINT::FitParameter& k, MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm, MINT::FitParameter& eff_tag,MINT::FitParameter& w ):
    _r(r),_delta(delta),_gamma(gamma),_k(k),_tau(tau),_dGamma(dGamma),_dm(dm),_eff_tag(eff_tag), _w(w) {;}
    
};

enum basisType { noBasis=0  ,  expBasis= 3
    , sinBasis=13,  cosBasis=23
    , sinhBasis=63, coshBasis=53 };

// Full DsK like time PDF with additional coherence factor 
class FullTimePdf : public MINT::PdfBase<IDalitzEvent>
{
protected:
    FitParameter& _C;
    FitParameter& _D;
    FitParameter& _D_bar;
    FitParameter& _S;
    FitParameter& _S_bar;
    FitParameter& _k;
    
    FitParameter& _tau;
    FitParameter& _dGamma;
    FitParameter& _dm;
    FitParameter& _eff_tag;
    FitParameter& _w;
    
    // cast of MINT parameters to B2DX parameters
    RooRealVar* _r_t;
    RooRealVar* _r_dt;
    
    RooRealVar* _r_scale_mean_dt;
    RooRealVar* _r_scale_sigma_dt;
    RooAbsPdf* _pdf_sigma_t;
    
    RooCubicSplineFun* _spline;
    RooGaussEfficiencyModel* _efficiency;
    
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        _r_t->setVal(t);
        _r_dt->setVal(dt);
        
        if(t < _min_TAU || t > _max_TAU )return 0.;
        
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        const double D = 1./2. * ((1.+f) * _D + (1.-f) * _D_bar);
        const double S = 1./2. * ((1.+f) * _S + (1.-f) * _S_bar);
        
        const double cosh_term = _efficiency->evaluate(coshBasis,_tau,_dm,_dGamma);
        const double cos_term = _efficiency->evaluate(cosBasis,_tau,_dm,_dGamma);
        const double sinh_term = _efficiency->evaluate(sinhBasis,_tau,_dm,_dGamma);
        const double sin_term = _efficiency->evaluate(sinBasis,_tau,_dm,_dGamma);
        
        const double val = 
        (
         (2.-fabs(q))*cosh_term
         +f*q*(1.-2.*w)* _C *cos_term
         -(2.-fabs(q))*2.0*_k* D *sinh_term
         -f*2.0*q*(1.-2.*w)*_k* S *sin_term
         )*e_eff * _pdf_sigma_t->getVal();
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double D = 1./2. * ((1.+f) * _D + (1.-f) * _D_bar);
                
        double val = un_normalised(evt);
        
        double norm = 2. * (_efficiency->analyticalIntegral(coshBasis,_tau,_dm,_dGamma) - 2. * _k * D * _efficiency->analyticalIntegral(sinhBasis,_tau,_dm,_dGamma) );
        
        return val/norm;
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }
    
    FullTimePdf(MINT::FitParameter& C,MINT::FitParameter& D, MINT::FitParameter& D_bar,MINT::FitParameter& S, MINT::FitParameter& S_bar,MINT::FitParameter& k, MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm, MINT::FitParameter& eff_tag,MINT::FitParameter& w ):
    _C(C),_D(D),_D_bar(D_bar),_S(S),_S_bar(S_bar),_k(k),_tau(tau),_dGamma(dGamma),_dm(dm),_eff_tag(eff_tag), _w(w), _min_TAU("min_TAU", 0.4), _max_TAU("max_TAU", 10.)
    {
        // Init B2DX parameters
        _r_t = new RooRealVar("t", "time", _min_TAU, _max_TAU);
        _r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.25);
        
        _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", 0);
        _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", 1.2);
        
        //SPLINE KNOTS
        NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        vector<double> myBinning = knot_positions.getVector();
        
        NamedParameter<double> knot_values("knot_values", 0.38,0.63,0.86,1.05,1.14,1.24,1.22);
        vector<double> values = knot_values.getVector() ;
        
        //SPLINE COEFFICIENTS
        RooArgList tacc_list;
        for(int i= 0; i< values.size(); i++){
            tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 10.0)));
        }
        tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), 1.0)));
        
        RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
        
        tacc_list.add(*coeff_last);	
        
        //CUBIC SPLINE FUNCTION 
        _spline = new RooCubicSplineFun("splinePdf", "splinePdf", *_r_t, myBinning, tacc_list);        
        _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt, *_r_scale_mean_dt, *_r_scale_sigma_dt );
        
        // Plot acceptance
        TH1F *h_spline = new TH1F("", "", 100, _min_TAU, _max_TAU);
        for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
            _r_t->setVal(h_spline->GetXaxis()->GetBinCenter(i));
            h_spline->SetBinContent(i,_spline->getVal());
        }
        
        TCanvas* c = new TCanvas();
        h_spline->SetLineColor(kRed);
        h_spline->Draw("histc");
        c->Print("spline.eps");
        c->Print("spline.pdf");
        
        // Init pdf for sigma_t
        _pdf_sigma_t = new RooGenericPdf("pdf_sigma_t","pow(7. / @1, 7) / 720. * pow(@0, 6) * exp(-7. * @0 / @1)",RooArgList(*_r_dt, RooRealConstant::value(0.04)));
        
    }
    
};

// Full Time PDF in term of r, gamma, delta instead of CP coefficients
class FullTimePdf_mod : public MINT::PdfBase<IDalitzEvent>
{
protected:
    // actual fit parameters
    FitParameter& _r;
    FitParameter& _delta;
    FitParameter& _gamma;
    FitParameter& _k;
    
    FitParameter& _tau;
    FitParameter& _dGamma;
    FitParameter& _dm;
    FitParameter& _eff_tag;
    FitParameter& _w;
    
    // cast of MINT parameters to B2DX parameters
    RooRealVar* _r_t;
    RooRealVar* _r_dt;
        
    RooRealVar* _r_scale_mean_dt;
    RooRealVar* _r_scale_sigma_dt;
    RooAbsPdf* _pdf_sigma_t;
    
    RooCubicSplineFun* _spline;
    RooGaussEfficiencyModel* _efficiency;
    
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2));
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        _r_t->setVal(t);
        _r_dt->setVal(dt);
        
        if(t < _min_TAU || t > _max_TAU )return 0.;
        
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        const double C = (1.-(double)_r*(double)_r)/(1.+(double)_r*(double)_r);
        const double D = (double)_r * cos(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
        const double S = (double)_r * sin(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
        
        const double cosh_term = _efficiency->evaluate(coshBasis,_tau,_dm,_dGamma);
        const double cos_term = _efficiency->evaluate(cosBasis,_tau,_dm,_dGamma);
        const double sinh_term = _efficiency->evaluate(sinhBasis,_tau,_dm,_dGamma);
        const double sin_term = _efficiency->evaluate(sinBasis,_tau,_dm,_dGamma);
        
        const double val = 
        (
         (2.-fabs(q))*cosh_term
         +f*q*(1.-2.*w)* C *cos_term
         -(2.-fabs(q))*2.0*_k* D *sinh_term
         -f*2.0*q*(1.-2.*w)*_k* S *sin_term
         )*e_eff * _pdf_sigma_t->getVal();
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double D = (double)_r * cos(((double)_delta-(double)_gamma*f)/360.*2*pi)/(1.+(double)_r*(double)_r);
                
        double val = un_normalised(evt);
        
        double norm = 2. * (_efficiency->analyticalIntegral(coshBasis,_tau,_dm,_dGamma) - 2. * _k * D * _efficiency->analyticalIntegral(sinhBasis,_tau,_dm,_dGamma) );
        
        return val/norm;
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }
    
    FullTimePdf_mod(MINT::FitParameter& r,MINT::FitParameter& delta, MINT::FitParameter& gamma,MINT::FitParameter& k, MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm, MINT::FitParameter& eff_tag,MINT::FitParameter& w ):
    _r(r),_delta(delta),_gamma(gamma),_k(k),_tau(tau),_dGamma(dGamma),_dm(dm),_eff_tag(eff_tag), _w(w), _min_TAU("min_TAU", 0.4), _max_TAU("max_TAU", 10.)
    {
        
        // Init B2DX parameters
        _r_t = new RooRealVar("t", "time", _min_TAU, _max_TAU);
        _r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.25);
        
        _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", 0);
        _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", 1.2);
        
        //SPLINE KNOTS
        NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        vector<double> myBinning = knot_positions.getVector();
        
        NamedParameter<double> knot_values("knot_values", 0.38,0.63,0.86,1.05,1.14,1.24,1.22);
        vector<double> values = knot_values.getVector() ;
        
        //SPLINE COEFFICIENTS
        RooArgList tacc_list;
        for(int i= 0; i< values.size(); i++){
            tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 10.0)));
        }
        tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), 1.0)));
        
        RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
        
        tacc_list.add(*coeff_last);	
        
        //CUBIC SPLINE FUNCTION 
        _spline = new RooCubicSplineFun("splinePdf", "splinePdf", *_r_t, myBinning, tacc_list);        
        _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt, *_r_scale_mean_dt, *_r_scale_sigma_dt );

        // Plot acceptance
        TH1F *h_spline = new TH1F("", "", 100, _min_TAU, _max_TAU);
        for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
            _r_t->setVal(h_spline->GetXaxis()->GetBinCenter(i));
            h_spline->SetBinContent(i,_spline->getVal());
        }
        
        TCanvas* c = new TCanvas();
        h_spline->SetLineColor(kRed);
        h_spline->Draw("histc");
        c->Print("spline.eps");
        c->Print("spline.pdf");
    
        // Init pdf for sigma_t
        _pdf_sigma_t = new RooGenericPdf("pdf_sigma_t","pow(7. / @1, 7) / 720. * pow(@0, 6) * exp(-7. * @0 / @1)",RooArgList(*_r_dt, RooRealConstant::value(0.04)));
       
    }
    
};

void fullTimeFit(){

    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    NamedParameter<int>  Nevents("Nevents", 100);
    NamedParameter<double>  pdf_max("pdf_max", 100);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<int> nIterations("nIterations", 1);

    
    
    // Calculate pulls
    gDirectory->cd();
    TFile* paraFile = new TFile(((string)OutputDir+"pull.root").c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;
    TTree* tree = new TTree("pull","pull");
    double C_pull,D_pull,Dbar_pull,S_pull,Sbar_pull;
    TBranch* br_C_pull = tree->Branch( "C_pull", &C_pull, "C_pull/D" );
    TBranch* br_D_pull = tree->Branch( "D_pull", &D_pull, "D_pull/D" );
    TBranch* br_Dbar_pull = tree->Branch( "Dbar_pull", &Dbar_pull, "Dbar_pull/D" );
    TBranch* br_S_pull = tree->Branch( "S_pull", &S_pull, "S_pull/D" );
    TBranch* br_Sbar_pull = tree->Branch( "Sbar_pull", &Sbar_pull, "Sbar_pull/D" );
    
    for(int i = 0; i < nIterations; i++){
    
        FitParameter  C("C",0,1,0.1);
        FitParameter  D("D",0,-1.3,0.1);
        FitParameter  D_bar("D_bar",0,-0.8,0.1);
        FitParameter  S("S",0,-1.25,0.1);
        FitParameter  S_bar("S_bar",0,0.08,0.1);
        FitParameter  k("k");
        
        FitParameter  tau("tau");
        FitParameter  dGamma("dGamma");
        FitParameter  dm("dm");
        FitParameter  eff_tag("eff_tag");
        FitParameter  mistag("mistag");
        
        FullTimePdf t_pdf(C,D,D_bar,S,S_bar,k,tau,dGamma,dm,eff_tag,mistag );

        DalitzEventList eventList;
        DalitzEventList eventList_f, eventList_fbar;
        
        DalitzEventPattern _pat(pat);
        DalitzEvent evt(_pat);
        evt.generateThisToPhaseSpace();
        
        RooRealVar* _r_t = new RooRealVar("t", "time", min_TAU, max_TAU);
        RooRealVar* _r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.1);
        RooRealVar* _r_mistag = new RooRealVar("mistag", "mistag",0., 0.5);
        
        RooCategory* _r_f = new RooCategory("qf", "qf");
        _r_f->defineType("h+", +1);
        _r_f->defineType("h-", -1);
        
        RooCategory* _r_q = new RooCategory("qt", "qt");
        _r_q->defineType("B+", +1);
        _r_q->defineType("B-", -1) ;   
        _r_q->defineType("untagged", 0);    
        
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(*_r_t,*_r_dt,*_r_q,*_r_mistag,*_r_f));

        //simple hit and miss
        for(int i = 0; i < Nevents; i++){
            while(true){
                const double t = ranLux.Exp(tau);
                const double dt = ranLux.Uniform(0., 0.25);
                const double q_rand = ranLux.Uniform();
                int q = 0;
                if (q_rand < 1./3.) q = -1;
                if (q_rand > 2./3.) q = 1;
                int f;
                if(i<=Nevents/2)f = 1; 
                else f = -1;
                
                evt.setValueInVector(0, t);
                evt.setValueInVector(1, dt);
                evt.setValueInVector(2, q);
                evt.setValueInVector(3, mistag);
                evt.setValueInVector(4, f);
                
                const double pdfVal = t_pdf.un_normalised(evt);
                
                double maxVal = exp(-fabs(t)/(tau))/(tau)*pdf_max;
                const double height = ranLux.Uniform(0,maxVal);
                
                //Safety check on the maxmimum generated height
                if( pdfVal > maxVal ){
                    std::cout << "ERROR: PDF above determined maximum." << std::endl;
                    std::cout << pdfVal << " > " << maxVal << std::endl;
                    //exit(1);
                    pdf_max = pdf_max * 2.;
                }
                
                //Hit-and-miss
                if( height < pdfVal ){
                    eventList.Add(evt);
                    _r_t->setVal(t);
                    _r_dt->setVal(dt);
                    _r_q->setIndex(q);
                    _r_mistag->setVal(mistag);
                    _r_f->setIndex(f);
                    data->add(RooArgSet(*_r_t,*_r_dt,*_r_q,*_r_mistag,*_r_f));
                    if(f==1)eventList_f.Add(evt);
                    else eventList_fbar.Add(evt);
                    break;
                }
            }
        }

        // Fit with MINT Pdf
        Neg2LL fcn_t(t_pdf, eventList_f);
        Neg2LL fcn_t_bar(t_pdf, eventList_fbar);
        Neg2LLSum neg2LLSum_t(&fcn_t,&fcn_t_bar);
        
        Minimiser mini_t(&neg2LLSum_t);
        mini_t.doFit();
        mini_t.printResultVsInput();
        
        TCanvas* c = new TCanvas();
        if(i==0){
            int nBinst = 100;
            TH1D* h_t = new TH1D("",";t (ps);Events (norm.) ",nBinst,0,max_TAU);
            TH1D* h_dt = new TH1D("",";dt (ps);Events (norm.) ",nBinst,0,0.25);
            
            for (int i=0; i<eventList.size(); i++) {
                h_t->Fill(eventList[i].getValueFromVector(0));
                h_dt->Fill(eventList[i].getValueFromVector(1));
            }     
            
            TH1D* h_t_fit = new TH1D("",";t",nBinst,0,max_TAU);
            TH1D* h_dt_fit = new TH1D("",";dt (ps);Events (norm.) ",nBinst,0,0.25);
            
            for(int i = 0; i < 2000000; i++){
                
                const double t = ranLux.Exp(tau);
                const double dt = ranLux.Uniform(0., 0.25);
                const double q_rand = ranLux.Uniform();
                int q = 0;
                if (q_rand < 1./3.) q = -1;
                if (q_rand > 2./3.) q = 1;
                int f;
                if(i<=Nevents/2)f = 1; 
                else f = -1;
                
                evt.setValueInVector(0, t);
                evt.setValueInVector(1, dt);
                evt.setValueInVector(2, q);
                evt.setValueInVector(3, mistag);
                evt.setValueInVector(4, f);
                
                const double pdfVal = t_pdf.un_normalised(evt);
                double weight = pdfVal/exp(-fabs(t)/(tau));
                
                h_t_fit->Fill(t,weight);
                h_dt_fit->Fill(dt,weight);
                
            }
            
            h_t->SetLineColor(kBlack);
            h_t->DrawNormalized("e1",1);
            h_t_fit->SetLineColor(kBlue);
            h_t_fit->SetLineWidth(3);
            h_t_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_t.eps").c_str());
            gPad->SetLogy(1);
            c->Print(((string)OutputDir+"h_t_log.eps").c_str());
            gPad->SetLogy(0);
            
            h_dt->SetLineColor(kBlack);
            h_dt->DrawNormalized("e1",1);
            h_dt_fit->SetLineColor(kBlue);
            h_dt_fit->SetLineWidth(3);
            h_dt_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_dt.eps").c_str());
        }
        
        // Fit with B2DX Pdf
        RooRealVar* _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", 0);
        RooRealVar* _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", 1.2);
        
        NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        vector<double> myBinning = knot_positions.getVector();
        
        NamedParameter<double> knot_values("knot_values", 0.38,0.63,0.86,1.05,1.14,1.24,1.22);
        vector<double> values = knot_values.getVector() ;
        
        RooArgList tacc_list;
        for(int i= 0; i< values.size(); i++){
            tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i]*0.1)));
        }
        tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), .1)));
        
        RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
        
        tacc_list.add(*coeff_last);	
        
        RooCubicSplineFun* _spline = new RooCubicSplineFun("splinePdf", "splinePdf", *_r_t, myBinning, tacc_list);        
        RooGaussEfficiencyModel* _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt, *_r_scale_mean_dt, *_r_scale_sigma_dt );
        RooGenericPdf* _pdf_sigma_t = new RooGenericPdf("pdf_sigma_t","pow(7. / @1, 7) / 720. * pow(@0, 6) * exp(-7. * @0 / @1)",RooArgList(*_r_dt, RooRealConstant::value(0.04)));
        
        
        RooRealVar _r_tau("tau", "decay time", tau);    
        RooRealVar _r_dgamma("dgamma", "dgamma", dGamma);
        RooRealVar _r_dm("dm", "dm", dm);

        RooRealVar _r_C("C", "C",C,-4,4 );
        RooFormulaVar _r_Cbar("Cbar","-1. * @0",RooArgList(_r_C));
        RooRealVar _r_D("D", "D",D,-4,4  );
        RooRealVar _r_Dbar("Dbar", "Dbar",D_bar,-4,4 );
        RooRealVar _r_S("S", "S",S,-4,4 );
        RooRealVar _r_Sbar("Sbar", "Sbar",S_bar,-4,4 );
        
        DecRateCoeff_Bd cosh_coeff("cosh_coeff",
                        "cosh_coeff",
                        DecRateCoeff_Bd::kCosh ,
                        *_r_f,
                        RooRealConstant::value(1.),
                        RooRealConstant::value(1.),
                        *_r_q,
                        *_r_mistag,
                        RooRealConstant::value(0.),
                        RooRealConstant::value(1.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(eff_tag),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.)); 
                        
        DecRateCoeff_Bd cos_coeff("cos_coeff",
                                   "cos_coeff",
                                   DecRateCoeff_Bd::kCos ,
                                   *_r_f,
                                   _r_C,
                                   _r_Cbar,
                                   *_r_q,
                                  *_r_mistag,
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(1.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(eff_tag),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.));  
                                            
        DecRateCoeff_Bd sinh_coeff("sinh_coeff",
                                   "sinh_coeff",
                                   DecRateCoeff_Bd::kSinh ,
                                   *_r_f,
                                   _r_D,
                                   _r_Dbar,
                                   *_r_q,
                                   *_r_mistag,
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(1.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(eff_tag),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.));  
                                   
        DecRateCoeff_Bd sin_coeff("sin_coeff",
                                   "sin_coeff",
                                   DecRateCoeff_Bd::kSin,
                                   *_r_f,
                                   _r_S,
                                   _r_Sbar,
                                   *_r_q,
                                   *_r_mistag,
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(1.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(eff_tag),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.));       
                    
        
        /*                           
        RooFormulaVar sigma_dt("sigma_dt","@0*@1",RooArgList(*_r_dt,*_r_scale_sigma_dt));
        RooGaussModel resmodel("resmodel", "resolution model", *_r_t, RooRealConstant::value(0), sigma_dt);
         
        RooBDecay gendecay("gendecay", "decay time PDF for fitting",
                           *_r_t,_r_tau, _r_dgamma, 
                           cosh_coeff,sinh_coeff,cos_coeff,sin_coeff,
                           _r_dm, resmodel, RooBDecay::SingleSided);                                                                                                            
        
        RooEffProd genpdf_t("genpdf_t", "PDF for generation", gendecay, *_spline);
        
        RooProdPdf genpdf("genpdf", "Full PDF for generation",
                          RooArgSet(*_pdf_sigma_t),
                          RooFit::Conditional(RooArgSet(genpdf_t), RooArgSet(*_r_t,*_r_q,*_r_mistag,*_r_f)));                                                                                                                                                                               
                                          
        RooDataSet* data_2 = genpdf.generate(RooArgSet(*_r_t,*_r_dt,*_r_q,*_r_mistag,*_r_f), 300);
        */                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                            
        RooBDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting",
                           *_r_t,_r_tau, _r_dgamma, 
                           cosh_coeff,sinh_coeff,cos_coeff,sin_coeff,
                           _r_dm, *_efficiency, RooBDecay::SingleSided);
        
        RooProdPdf fitpdf("fitpdf", "Full PDF for fitting",
                          RooArgSet(*_pdf_sigma_t),
                          RooFit::Conditional(RooArgSet(fitpdf_t), RooArgSet(*_r_t,*_r_q,*_r_mistag,*_r_f)));
        
        
        fitpdf.fitTo(*data);
        if(i==0){
            RooPlot* tframefit = _r_t->frame();
            data->plotOn(tframefit);
            /// WTF is this doing ???
            fitpdf.plotOn(tframefit,ProjWData(RooArgSet(*_r_dt,*_r_q,*_r_mistag,*_r_f),*data));
            tframefit->Draw();
            c->Print("B2DX_fit.eps");        
            cout << " C = " << _r_C.getVal() << " ; Pull = " << (_r_C.getVal()-C)/_r_C.getError() << endl;
            cout << " D = " << _r_D.getVal() << " ; Pull = " << (_r_D.getVal()-D)/_r_D.getError() << endl;
            cout << " Dbar = " << _r_Dbar.getVal() << " ; Pull = " << (_r_Dbar.getVal()-D_bar)/_r_Dbar.getError() << endl;
            cout << " S = " << _r_S.getVal() << " ; Pull = " << (_r_S.getVal()-S)/_r_S.getError() << endl;
            cout << " Sbar = " << _r_Sbar.getVal() << " ; Pull = " << (_r_Sbar.getVal()-S_bar)/_r_Sbar.getError() << endl;
        }
        C_pull = (_r_C.getVal()-C)/_r_C.getError();
        D_pull = (_r_D.getVal()-D)/_r_D.getError();
        Dbar_pull =(_r_Dbar.getVal()-D_bar)/_r_Dbar.getError();
        S_pull = (_r_S.getVal()-S)/_r_S.getError();
        Sbar_pull = (_r_Sbar.getVal()-S_bar)/_r_Sbar.getError();
        
        tree->Fill();
    
    }
    
    paraFile->cd();
    tree->SetDirectory(paraFile);
    tree->Write();
    paraFile->Close();
    delete paraFile;
    
    return;
}

void timeFit_mod(){
    
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    NamedParameter<int>  Nevents("Nevents", 100);
    NamedParameter<double>  pdf_max("pdf_max", 100);
    
    FitParameter  r("r");
    FitParameter  delta("delta");
    FitParameter  gamma("gamma");
    FitParameter  k("k");
    
    FitParameter  tau("tau");
    FitParameter  dGamma("dGamma");
    FitParameter  dm("dm");
    FitParameter  eff_tag("eff_tag");
    FitParameter  mistag("mistag");
    
    TimePdf_mod t_pdf(r,delta,gamma,k,tau,dGamma,dm,eff_tag,mistag );
    
    DalitzEventList eventList;
    DalitzEventList eventList_f, eventList_fbar;
    
    DalitzEventPattern _pat(pat);
    DalitzEvent evt(_pat);
    evt.generateThisToPhaseSpace();
    
    //simple hit and miss
    for(int i = 0; i < Nevents; i++){
        while(true){
            const double t = ranLux.Exp(tau);
            const double q_rand = ranLux.Uniform();
            int q = 0;
            if (q_rand < 1./3.) q = -1;
            if (q_rand > 2./3.) q = 1;
            int f;
            if(i<=Nevents/2)f = 1; 
            else f = -1;
                        
            evt.setValueInVector(0, t);
            evt.setValueInVector(1, 0);
            evt.setValueInVector(2, q);
            evt.setValueInVector(3, mistag);
            evt.setValueInVector(4, f);
            
            const double pdfVal = t_pdf.un_normalised(evt);
            
            double maxVal = exp(-fabs(t)/(tau))/(tau)*pdf_max;
            const double height = ranLux.Uniform(0,maxVal);
            
            //Safety check on the maxmimum generated height
            if( pdfVal > maxVal ){
                std::cout << "ERROR: PDF above determined maximum." << std::endl;
                std::cout << pdfVal << " > " << maxVal << std::endl;
                //exit(1);
                pdf_max = pdf_max * 2.;
            }
            
            //Hit-and-miss
            if( height < pdfVal ){
                eventList.Add(evt);
                if(f==1)eventList_f.Add(evt);
                else eventList_fbar.Add(evt);
                break;
            }
        }
    }
    
    // Fit
    Neg2LL fcn_t(t_pdf, eventList_f);
    Neg2LL fcn_t_bar(t_pdf, eventList_fbar);
    Neg2LLSum neg2LLSum_t(&fcn_t,&fcn_t_bar);
    
    Minimiser mini_t(&neg2LLSum_t);
    mini_t.doFit();
    mini_t.printResultVsInput();
    return;
}

void fullTimeFit_mod(){
    
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    NamedParameter<int>  Nevents("Nevents", 100);
    NamedParameter<double>  pdf_max("pdf_max", 100);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    
    FitParameter  r("r");
    FitParameter  delta("delta");
    FitParameter  gamma("gamma");
    FitParameter  k("k");
    
    FitParameter  tau("tau");
    FitParameter  dGamma("dGamma");
    FitParameter  dm("dm");
    FitParameter  eff_tag("eff_tag");
    FitParameter  mistag("mistag");
    
    FullTimePdf_mod t_pdf(r,delta,gamma,k,tau,dGamma,dm,eff_tag,mistag );
        
    DalitzEventList eventList;
    DalitzEventList eventList_f, eventList_fbar;
    
    DalitzEventPattern _pat(pat);
    DalitzEvent evt(_pat);
    evt.generateThisToPhaseSpace();
    
    RooRealVar* _r_t = new RooRealVar("t", "time", min_TAU, max_TAU);
    RooRealVar* _r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.25);
    RooDataSet* data = new RooDataSet("data","data",RooArgSet(*_r_t,*_r_dt));
    
    //simple hit and miss
    for(int i = 0; i < Nevents; i++){
        while(true){
            const double t = ranLux.Exp(tau);
            const double dt = ranLux.Uniform(0., 0.25);
            const double q_rand = ranLux.Uniform();
            int q = 0;
            if (q_rand < 1./3.) q = -1;
            if (q_rand > 2./3.) q = 1;
            int f;
            if(i<=Nevents/2)f = 1; 
            else f = -1;
            
            evt.setValueInVector(0, t);
            evt.setValueInVector(1, dt);
            evt.setValueInVector(2, q);
            evt.setValueInVector(3, mistag);
            evt.setValueInVector(4, f);
            
            const double pdfVal = t_pdf.un_normalised(evt);
            
            double maxVal = exp(-fabs(t)/(tau))/(tau)*pdf_max;
            const double height = ranLux.Uniform(0,maxVal);
            
            //Safety check on the maxmimum generated height
            if( pdfVal > maxVal ){
                std::cout << "ERROR: PDF above determined maximum." << std::endl;
                std::cout << pdfVal << " > " << maxVal << std::endl;
                //exit(1);
                pdf_max = pdf_max * 2.;
            }
            
            //Hit-and-miss
            if( height < pdfVal ){
                eventList.Add(evt);
                _r_t->setVal(t);
                _r_dt->setVal(dt);
                data->add(RooArgSet(*_r_t,*_r_dt));
                if(f==1)eventList_f.Add(evt);
                else eventList_fbar.Add(evt);
                break;
            }
        }
    }
    
    // Fit
    Neg2LL fcn_t(t_pdf, eventList_f);
    Neg2LL fcn_t_bar(t_pdf, eventList_fbar);
    Neg2LLSum neg2LLSum_t(&fcn_t,&fcn_t_bar);
    
    Minimiser mini_t(&neg2LLSum_t);
    mini_t.doFit();
    mini_t.printResultVsInput();
    
    int nBinst = 100;
    TH1D* h_t = new TH1D("",";t (ps);Events (norm.) ",nBinst,0,max_TAU);
    TH1D* h_dt = new TH1D("",";dt (ps);Events (norm.) ",nBinst,0,0.25);

    for (int i=0; i<eventList.size(); i++) {
        h_t->Fill(eventList[i].getValueFromVector(0));
        h_dt->Fill(eventList[i].getValueFromVector(1));
    }     

    TH1D* h_t_fit = new TH1D("",";t",nBinst,0,max_TAU);
    TH1D* h_dt_fit = new TH1D("",";dt (ps);Events (norm.) ",nBinst,0,0.25);

    for(int i = 0; i < 20000000; i++){
    
        const double t = ranLux.Exp(tau);
        const double dt = ranLux.Uniform(0., 0.25);
        const double q_rand = ranLux.Uniform();
        int q = 0;
        if (q_rand < 1./3.) q = -1;
        if (q_rand > 2./3.) q = 1;
        int f;
        if(i<=Nevents/2)f = 1; 
        else f = -1;
        
        evt.setValueInVector(0, t);
        evt.setValueInVector(1, dt);
        evt.setValueInVector(2, q);
        evt.setValueInVector(3, mistag);
        evt.setValueInVector(4, f);
        
        const double pdfVal = t_pdf.un_normalised(evt);
        double weight = pdfVal/exp(-fabs(t)/(tau));

        h_t_fit->Fill(t,weight);
        h_dt_fit->Fill(dt,weight);

    }
    
    TCanvas* c = new TCanvas();
    h_t->SetLineColor(kBlack);
    h_t->DrawNormalized("e1",1);
    h_t_fit->SetLineColor(kBlue);
    h_t_fit->SetLineWidth(3);
    h_t_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_t.eps").c_str());
    gPad->SetLogy(1);
    c->Print(((string)OutputDir+"h_t_log.eps").c_str());
    gPad->SetLogy(0);

    h_dt->SetLineColor(kBlack);
    h_dt->DrawNormalized("e1",1);
    h_dt_fit->SetLineColor(kBlue);
    h_dt_fit->SetLineWidth(3);
    h_dt_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_dt.eps").c_str());
    
    
    ////
    RooRealVar* _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", 0);
    RooRealVar* _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", 1.2);
    
    NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
    vector<double> myBinning = knot_positions.getVector();
    
    NamedParameter<double> knot_values("knot_values", 0.38,0.63,0.86,1.05,1.14,1.24,1.22);
    vector<double> values = knot_values.getVector() ;
    
    RooArgList tacc_list;
    for(int i= 0; i< values.size(); i++){
        tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i])));
    }
    tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), 1.0)));
    
    RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
    
    tacc_list.add(*coeff_last);	
    
    RooCubicSplineFun* _spline = new RooCubicSplineFun("splinePdf", "splinePdf", *_r_t, myBinning, tacc_list);        
    RooGaussEfficiencyModel* _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt, *_r_scale_mean_dt, *_r_scale_sigma_dt );
    RooGenericPdf* _pdf_sigma_t = new RooGenericPdf("pdf_sigma_t","pow(7. / @1, 7) / 720. * pow(@0, 6) * exp(-7. * @0 / @1)",RooArgList(*_r_dt, RooRealConstant::value(0.04)));
    
    
    RooRealVar _r_tau("tau", "decay time", 1., 0., 2.);    
    RooRealVar _r_dgamma("dgamma", "dgamma", 0.1,0.,0.5);
    RooRealVar _r_dm("dm", "dm", 20);
    RooRealVar _r_f0("f0", "f0",1 );
    RooRealVar _r_f1("f1", "f1",0 );
    RooRealVar _r_f2("f2", "f2",0 );
    RooRealVar _r_f3("f3", "f3",0 );
    
     RooBDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting",
     *_r_t,_r_tau, _r_dgamma, 
     _r_f0, _r_f1, _r_f2, _r_f3,
     _r_dm, *_efficiency, RooBDecay::SingleSided);
     
     RooProdPdf fitpdf("fitpdf", "Full PDF for fitting",
     RooArgSet(*_pdf_sigma_t),
     RooFit::Conditional(RooArgSet(fitpdf_t), RooArgSet(*_r_t)));
     
    fitpdf.fitTo(*data);
    
    RooPlot* tframefit = _r_t->frame();
    data->plotOn(tframefit);
    fitpdf.plotOn(tframefit);
    tframefit->DrawClone();
    c->Print("B2DX_fit.eps");        
    return;
}

void spline_test(){

    TCanvas* c = new TCanvas();
    //c->Divide(1, 2);
    
    // define observables
    RooRealVar t("t", "time", 0., 14);
    RooRealVar dt("dt", "per-candidate time resolution estimate",0.001, 0.25);
    // define PDF parameters
    RooRealVar tau("tau", "decay time", 1., 0., 14);
    RooRealVar mean("mean", "mean", 0);
    RooRealVar sigma_t("sigma_t", "sigma_t", 0.044);
    
    RooRealVar dgamma("dgamma", "dgamma", 0.1);
    RooRealVar dm("dm", "dm", 20);
    RooRealVar f0("f0", "f0",1 );
    RooRealVar f1("f1", "f1",1 );
    RooRealVar f2("f2", "f2",0 );
    RooRealVar f3("f3", "f3",0 );
    
    // Generate random hist? Like 5 bins : from 1 to 0.5 to imitate
    // efficiency
    TH1F *hist = new TH1F("hist", "efficiency histogram", 10, 0., 14);
    for (int i = 0; i < 10; i++) {
        hist->SetBinContent(i + 1, 0.99 - 0.1 * i);
    }
        
    // define spline (from histogram)
    RooCubicSplineFun spline("spline", "spline from hist", t, hist);
        
    TH1F *h_spline = new TH1F("", "", 100, 0., 14);
    
    for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
            t.setVal(h_spline->GetXaxis()->GetBinCenter(i));
            h_spline->SetBinContent(i,spline.getVal());
    }

    hist->Draw();
    h_spline->SetLineColor(kRed);
    h_spline->Draw("histcsame");
    c->Print("spline.eps");
    c->Print("spline.pdf");

    // build PDF for per-candidate decay time resolution:
    // This gives the probability to observe a predicted time resolution of
    // dt; the chosen functional form pdf has roughly the shape one would see
    // in a forward detector like LHCb:
    //
    // P(dt) = ((n + 1) / sigma) ^ (n + 1) / n! *
    //         t^n * exp(-(n + 1) * dt / sigma_t)
    //
    // [This functional form is normalised, and has expectation value sigma_t.
    // The parameter n can be varied to influence the steepness of the
    // distribution.]
    RooGenericPdf sigmapdf("sigmapdf","pow(7. / @1, 7) / 720. * pow(@0, 6) * exp(-7. * @0 / @1)",RooArgList(dt, sigma_t));
    
    /*
     * build PDF used for the generation stage
     */
    // build resolution model
    RooGaussModel resmodel("resmodel", "resolution model", t, mean, dt);
    // build pdf
    //RooDecay gendecay("decay", "decay time PDF", t, tau, resmodel,RooDecay::SingleSided);
                      
    RooBDecay gendecay("gendecay", "gendecay",
                       t,tau, dgamma, 
                       f0, f1, f2, f3,
                       dm, resmodel, RooBDecay::SingleSided);                 
                      
                      
    // P(t|dt) - the decay time pdf for a given decay time resolution
    RooEffProd genpdf_t("genpdf_t", "PDF for generation", gendecay, spline);
    // build the full P(t|dt) * P(dt)
    RooProdPdf genpdf("genpdf", "Full PDF for generation",
                      RooArgSet(sigmapdf),
                      RooFit::Conditional(RooArgSet(genpdf_t), RooArgSet(t)));
    // Generate with RooEffProd
    RooDataSet* data = genpdf.generate(RooArgSet(t, dt), 200);
    // Fit using the generator PDF - very slow!
    //genpdf.fitTo(*data, Timer(), Verbose(kFALSE));
    
    
    RooPlot* tframegen = t.frame(Title("Decay Time p.d.f. (generation)"));
    data->plotOn(tframegen);
    genpdf.plotOn(tframegen);
    tframegen->DrawClone();
    
    c->Print("splineFit.pdf");
    c->Print("splineFit.eps");
    
    
    /*
     * build PDF used for the fit stage
     */
    // build spline-based resolution model with effects of acceptance
    RooGaussEfficiencyModel efficiency("efficiency",
                                       "spline convolved with resolution", t, spline, mean, sigma_t);
    // build pdf
    // P(t|dt) - the decay time pdf for a given decay time resolution, this
    // version can do analytical integrals
    //RooDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting", t, tau, efficiency, RooDecay::SingleSided);
                      
    RooBDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting",
                          t,tau, dgamma, 
                          f0, f1, f2, f3,
                          dm, efficiency, RooBDecay::SingleSided);
                      
    // build the full P(t|dt) * P(dt)
    RooProdPdf fitpdf("fitpdf", "Full PDF for fitting",
                      RooArgSet(sigmapdf),
                      RooFit::Conditional(RooArgSet(fitpdf_t), RooArgSet(t)));
    // fit - should be a lot faster!
    fitpdf.fitTo(*data, Timer(), Verbose(kFALSE));
    
    RooPlot* tframefit = t.frame(Title("Decay Time p.d.f."));
    data->plotOn(tframefit);
    fitpdf.plotOn(tframefit);
    tframefit->DrawClone();
    
    c->Print("splineFitFast.pdf");
    c->Print("splineFitFast.eps");
    
    /*
    t.setRange("range",0,14);
    //cout << "int = " << fitpdf_t.analyticalIntegral(1,"range") << endl;
    
    cout << "int = " << efficiency.analyticalIntegral(1,"range") << endl;

    sigma_t.setVal(1);
    cout << "int = " << efficiency.analyticalIntegral(1,"range") << endl;

    sigma_t.setVal(0.0001);
    cout << "int = " << efficiency.analyticalIntegral(1,"range") << endl;
    
    vector<double> knots;
    knots.push_back(1);
    knots.push_back(2);
    
    RooCubicSplineKnot k(knots);
    */
}


int main(int argc, char** argv){

  time_t startTime = time(0);

  gROOT->ProcessLine(".x ../lhcbStyle.C");

  NamedParameter<int>  doToyFit("doToyFit", 1);
  NamedParameter<int>  useCPcoefficients("useCPcoefficients", 0);
  
  NamedParameter<int>  doSplineTest("doSplineTest", 0);

  if(doToyFit){
      if(useCPcoefficients)fullTimeFit();
      else fullTimeFit_mod();
  }
  if(doSplineTest)spline_test();
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
