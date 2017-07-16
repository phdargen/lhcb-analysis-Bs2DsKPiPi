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

#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"

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
using namespace RooStats;
using namespace MINT;


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


void timeFit(){

    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    NamedParameter<int>  Nevents("Nevents", 10);
    NamedParameter<double>  pdf_max("pdf_max", 10);
    NamedParameter<double>  t_max("t_max", 10);
    
    NamedParameter<double>  r("r",0.64);
    NamedParameter<double>  delta("delta",1.63);
    NamedParameter<double>  gamma2("gamma2",1.2);
    
    FitParameter  C("C",0,1,0.1);
    FitParameter  D("D",0,-1.3,0.1);
    FitParameter  D_bar("D_bar",0,-0.8,0.1);
    FitParameter  S("S",0,-1.25,0.1);
    FitParameter  S_bar("S_bar",0,0.08,0.1);
    FitParameter  k("k");
    
    FitParameter  tau("tau");
    FitParameter  dGamma("dGamma");
    FitParameter  dm("dm");
    FitParameter  gamma("gamma");
    FitParameter  eff_tag("eff_tag");
    FitParameter  w("w");
    
    TimePdf t_pdf(C,D,D_bar,S,S_bar,k,tau,dGamma,dm,eff_tag,w );

    DalitzEventList eventList;
    DalitzEventList eventList_f, eventList_fbar;
    
    DalitzEventPattern _pat(pat);
    DalitzEvent evt(_pat);
    evt.generateThisToPhaseSpace();

    //simple hit and miss
    for(int i = 0; i < Nevents; i++){
        while(true){
            const double t = ranLux.Exp(tau*2.);
            const double q_rand = ranLux.Uniform();
            int q = 0;
            if (q_rand < 1./3.) q = -1;
            if (q_rand > 2./3.) q = 1;
            const int f = 1; //ranLux.Uniform() > 0.5 ? +1 : -1;
            
            evt.setValueInVector(0, t);
            //evt.setValueInVector(1, _dt);
            evt.setValueInVector(1, 0);
            evt.setValueInVector(2, q);
            evt.setValueInVector(3, w);
            evt.setValueInVector(4, f);
            
            const double pdfVal = t_pdf.un_normalised(evt);
            
            double maxVal = exp(-fabs(t)/(tau*2.))/(tau*2.)*pdf_max;
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
                if(f==1) eventList_f.Add(evt);
                else eventList_fbar.Add(evt);
                
                break;
            }
        }
        while(true){
            const double t = ranLux.Exp(tau*2.);
            const double q_rand = ranLux.Uniform();
            int q = 0;
            if (q_rand < 1./3.) q = -1;
            if (q_rand > 2./3.) q = 1;
            const int f = -1; //ranLux.Uniform() > 0.5 ? +1 : -1;
            
            evt.setValueInVector(0, t);
            //evt.setValueInVector(1, _dt);
            evt.setValueInVector(1, 0);
            evt.setValueInVector(2, q);
            evt.setValueInVector(3, w);
            evt.setValueInVector(4, f);
            
            const double pdfVal = t_pdf.un_normalised(evt);
            
            double maxVal = exp(-fabs(t)/(tau*2.))/(tau*2.)*pdf_max;
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
                if(f==1) eventList_f.Add(evt);
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

void timeFit_mod(){
    
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    NamedParameter<int>  Nevents("Nevents", 10);
    NamedParameter<double>  pdf_max("pdf_max", 10);
    NamedParameter<double>  t_max("t_max", 10);
    
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
            const int f = 1; 
            
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
                eventList_f.Add(evt);                
                break;
            }
        }
        while(true){
            const double t = ranLux.Exp(tau);
            const double q_rand = ranLux.Uniform();
            int q = 0;
            if (q_rand < 1./3.) q = -1;
            if (q_rand > 2./3.) q = 1;
            const int f = -1; 
            
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
                eventList_fbar.Add(evt);                
                break;
            }
        }
    }
    
    
    return;

    
    // Fit
    Neg2LL fcn_t(t_pdf, eventList_f);
    Neg2LL fcn_t_bar(t_pdf, eventList_fbar);
    Neg2LLSum neg2LLSum_t(&fcn_t,&fcn_t_bar);
    
    Minimiser mini_t(&neg2LLSum_t);
    mini_t.doFit();
    mini_t.printResultVsInput();
    return;
}

void TD_test(){

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
    RooRealVar f1("f1", "f1",0 );
    RooRealVar f2("f2", "f2",0 );
    RooRealVar f3("f3", "f3",0 );
    
    // Generate random hist? Like 5 bins : from 1 to 0.5 to imitate
    // efficiency
    TH1F *hist = new TH1F("hist", "efficiency histogram", 5, 0., 14);
    //for (int i = 0; i < 10; i++) {
        // linearly decreasing for this test
      //  hist->SetBinContent(i + 1, 0.99 - 0.1 * i);
    //}
    
    hist->SetBinContent(1, 1);
    hist->SetBinContent(2, 1);
    hist->SetBinContent(3, 1);
    hist->SetBinContent(4, 1);
    hist->SetBinContent(5, 1.);

    
    // define spline (from histogram)
    RooCubicSplineFun spline("spline", "spline from hist", t, hist);
    
    hist->Draw();
    c->Print("eff.pdf");
    
    TH1F *h_spline = new TH1F("", "", 100, 0., 14);
    
    for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
            t.setVal(h_spline->GetXaxis()->GetBinCenter(i));
            h_spline->SetBinContent(i+1,spline.getVal());
    }

    hist->Draw();
    h_spline->SetLineColor(kRed);
    h_spline->Draw("histcsame");
    c->Print("spline.pdf");

    //return;
    
    
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
    
    
    //c->cd(1);
    RooPlot* tframegen = t.frame(Title("Decay Time p.d.f. (generation)"));
    data->plotOn(tframegen);
    genpdf.plotOn(tframegen);
    tframegen->DrawClone();
    
    c->Print("spline2.pdf");
    
    
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
                      
    
    /*
    51   _f0("f0", "Cosh Coefficient", this, f0),
    52   _f1("f1", "Sinh Coefficient", this, f1),
    53   _f2("f2", "Cos Coefficient", this, f2),
    54   _f3("f3", "Sin Coefficient", this, f3),
     */
                  
    RooBDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting",
                          t,tau, dgamma, 
                          f0, f1, f2, f3,
                          dm, efficiency, RooBDecay::SingleSided);
                      
    // build the full P(t|dt) * P(dt)
    RooProdPdf fitpdf("fitpdf", "Full PDF for fitting",
                      RooArgSet(sigmapdf),
                      RooFit::Conditional(RooArgSet(fitpdf_t), RooArgSet(t)));
    // fit - should be a lot faster!
    //fitpdf.fitTo(*data, Timer(), Verbose(kFALSE));
    
    RooPlot* tframefit = t.frame(Title("Decay Time p.d.f."));
    data->plotOn(tframefit);
    fitpdf.plotOn(tframefit);
    tframefit->DrawClone();
    
    c->Print("spline3.pdf");
    
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
}


int main(int argc, char** argv){

  time_t startTime = time(0);

  gROOT->ProcessLine(".x ../lhcbStyle.C");

  //TD_test();
  timeFit_mod();
  
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
