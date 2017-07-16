// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk), Jack Benton (Jack.B.Benton@bristol.ac.uk)
// status:  Fri 28 Jun 2013 11:21:01 GMT
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
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include <omp.h>


using namespace std;
using namespace MINT;
using namespace RooFit ;

class AmpsPdfFlexiFast
: public DalitzPdfBaseFlexiFastInteg
{
protected:
    TRandom* _localRnd;
    SignalGenerator* _sgGen;
    FromFileGenerator* _fileGen;
    IEventGenerator<IDalitzEvent>* _chosenGen;
    NamedParameter<std::string> _integratorSource;
    std::string _integratorFileName;
public:
    double un_normalised_noPs(IDalitzEvent& evt){
        double ampSq =  _amps->RealVal(evt);
        return ampSq;// * getEvent()->phaseSpace();
    }
    
    std::complex<double> ComplexVal_un_normalised_noPs(IDalitzEvent& evt){
        return  _amps->ComplexVal(evt);
    }
    
    AmpsPdfFlexiFast(const DalitzEventPattern& pat
                     , IFastAmplitudeIntegrable* amps
                     , MinuitParameterSet* mps
                     , double precision=1.e-4
                     , std::string method="efficient"
                     , std::string fname =  "SignalIntegrationEvents.root", bool genMoreEvents = false
                     )
    : DalitzPdfBaseFlexiFastInteg(pat, 0, amps, precision, mps)
    , _localRnd(0)
    , _sgGen(0)
    , _fileGen(0)
    , _chosenGen(0)
    , _integratorSource("IntegratorSource", (std::string) "new", (char*) 0)
    , _integratorFileName(fname)
    {
        cout << " AmpsPdfFlexiFast with integ method " << method << endl;
        bool nonFlat = "efficient" == method;
        bool generateNew = ((string)_integratorSource == (string)"new");
        if(nonFlat){
            cout << "AmpsPdfFlexiFast uses nonFlat integration." << endl;
            if(generateNew){
                _sgGen =  new SignalGenerator(pat);
                _sgGen->setWeighted();
                //_sgGen->setSaveEvents(_integratorFileName);
                _chosenGen = _sgGen;
            }else{
                // here, SignalGenerator is used by FromFileGenerator, to fill
                // up missing events in case more are needed than found in the
                // file.  Since we don't know with which random seed the
                // events in the file were generated, we supply a random
                // number generator with randomised seed.
                _localRnd = new TRandom3(time(0));
                _sgGen =  new SignalGenerator(pat, _localRnd);
                _sgGen->setWeighted();
                _sgGen->dontSaveEvents();// saving events is done by FromFileGenerator
                if(genMoreEvents) _fileGen   = new FromFileGenerator(_integratorFileName, _sgGen);
                else{
                    _fileGen = new FromFileGenerator(_integratorFileName, 0, "UPDATE");
                    cout << "not going to generate any more events" << endl;
                }
                _chosenGen = _fileGen;
            }
            this->setEventGenerator(_chosenGen);
        }else{
            cout << "AmpsPdfFlexiFast uses flat integration." << endl;
        }
    }
    
    IFastAmplitudeIntegrable* getAmpSum(){ return _amps;}
    
    ~AmpsPdfFlexiFast(){
        if(0 != _fileGen)  delete _fileGen;
        if(0 != _sgGen)    delete _sgGen;
        if(0 != _localRnd) delete _localRnd;
    }
};


class AmpRatio : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    int _f;
public:
    AmpRatio(FitParameter& re, FitParameter& im,int f)
    : _re(re), _im(im), _f(f) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result(_re,static_cast<double>(_f) * _im); 
        return result;
    }
};

class CPV_amp : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    int _sign;
public:
    CPV_amp(FitParameter& re, FitParameter& im, int sign)
    : _re(re), _im(im), _sign(sign) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result((double) ( 1.+  _re * (double) _sign),(double) (_im * (double) _sign) ); 
        return result;
    }
};

class CPV_amp_polar : virtual public IReturnComplex{
    FitParameter& _r;
    FitParameter& _delta;
    int _sign;
public:
    CPV_amp_polar(FitParameter& r, FitParameter& delta, int sign)
    : _r(r), _delta(delta), _sign(sign) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result= polar((double) sqrt( 1.+  _r * (double) _sign),(double) (_delta/360.*2*pi * (double) _sign) ); 
        return result;
    }
};


void AddScaledAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, FitAmpSum& fasCC, std::string name, counted_ptr<IReturnComplex>& r_plus, counted_ptr<IReturnComplex>& r_minus){

    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    FitAmpSum fasCC_2(*List);
    fasCC_2.CPConjugateSameFitParameters();
    fasCC_2.CConjugateFinalStateSameFitParameters();
    fas_2.multiply(r_plus); 
    fasCC_2.multiply(r_minus); 
    fas.addAsList(fas_2,1.);
    fasCC.addAsList(fasCC_2,1.);
}


std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    //fas.print();
    //fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
        const std::complex<double> amp = fas.getVal(eventList[i]) ;
        const std::complex<double> amp_bar = fas_bar.getVal(eventList[i])*phase_diff ;
        valK += amp_bar*conj(amp);
        val1 += norm(amp);
        val2 += norm(amp_bar);
    }
    
    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);
    
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;
    
    return result;
}

class AmpsPdfFlexiFastCPV : public MINT::PdfBase<IDalitzEvent>
, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps1;
    AmpsPdfFlexiFast* _amps2;
    AmpsPdfFlexiFast* _ampsSum;

    FitParameter& _r;
    FitParameter& _delta;
    FitParameter& _gamma;

    FitParameter& _tau;
    FitParameter& _dGamma;
    FitParameter& _dm;
    FitParameter& _eff_tag;
    FitParameter& _w;
    
    double _intA;
    double _intAbar;
    complex<double> _intAAbar;

public:
    void parametersChanged(){
        _ampsSum->parametersChanged();
        _intA = _ampsSum->integralForMatchingPatterns(true,1);
        _intAbar = _ampsSum->integralForMatchingPatterns(true,-1);        
        _intAAbar = _ampsSum->ComplexSumForMatchingPatterns(false);
    }
    void beginFit(){
        _ampsSum->beginFit();
         printIntegralVals();
    }
    void endFit(){
        printIntegralVals();
        _ampsSum->endFit();
    }
    
    void printIntegralVals(){
        cout << "intSum = " << _ampsSum->getIntegralValue() << endl;
        cout << "intA = " << _intA << endl;
        cout << "intAbar = " << _intAbar << endl;
        cout << "intAAbar = " << _intAAbar.real() << endl;
    }
    
    inline double un_normalised_noPs(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));

        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);

        const std::complex<double> amp = _amps1->ComplexVal_un_normalised_noPs(evt) ;
        const std::complex<double> amp_bar = _amps2->ComplexVal_un_normalised_noPs(evt) * phase_diff;
        
        const double val =  exp(-fabs(t)/(double)_tau) *
        (
         (2.-fabs(q))*(norm(amp) + norm(amp_bar))*cosh((double)_dGamma/2.*t)
         +f*q*(1.-2.*w)*(norm(amp) - norm(amp_bar)) *cos((double)_dm*t)
         -(2.-fabs(q))*2.0*real(amp_bar*conj(amp))*sinh((double)_dGamma/2.*t)
         -f*2.0*q*(1.-2.*w)*imag(amp_bar*conj(amp))*sin((double)_dm*t)
         )*e_eff;
        
        return val;
    }
    
    inline double un_normalised_noPs_time_integrated(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        const double e_eff = fabs(q)*_eff_tag + (1.-fabs(q))*(1.-_eff_tag);
        
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);
        
        const std::complex<double> amp = _amps1->ComplexVal_un_normalised_noPs(evt) ;
        const std::complex<double> amp_bar = _amps2->ComplexVal_un_normalised_noPs(evt) * phase_diff;
        
        double Gamma = 1./((double) _tau);
        
        const double val =
        (
         (2.-fabs(q))*(norm(amp) + norm(amp_bar))*(4.*Gamma/(4.*Gamma*Gamma-_dGamma*_dGamma))
         +f*q*(1.-2.*w)*(norm(amp) - norm(amp_bar)) *(Gamma/(Gamma*Gamma+_dm*_dm))
         -(2.-fabs(q))*2.0*real(amp_bar*conj(amp))*(2.*_dGamma/(4.*Gamma*Gamma-_dGamma*_dGamma))
         -f*2.0*q*(1.-2.*w)*imag(amp_bar*conj(amp))*(_dm/(Gamma*Gamma+_dm*_dm))
         )*e_eff;
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){

        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        const double val = un_normalised_noPs(evt);

        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const double Gamma = 1./((double) _tau);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);
        const double int_interference = (phase_diff*_intAAbar).real();
        
        if(_intA == -1 ){
            cout << "AmpsPdfFlexiFastCPV:: _norm = -1, should not have happened." << endl;
            throw "can't deal with that";
        }
        
        return val/(2.* ((_intA + r* r * _intAbar) * 4.*Gamma - int_interference * 2. * _dGamma )/ (4.*Gamma*Gamma-_dGamma*_dGamma));
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
    
    virtual DalitzHistoSet histoSet(){return _ampsSum->histoSet();}
    
    
    void doFinalStatsAndSaveForAmp12(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults.txt", const std::string& fnameROOT="fitFractions.root"){
        _amps1->redoIntegrator();
        _amps2->redoIntegrator();
        _amps1->doFinalStatsAndSave(min,fname,fnameROOT);
        _amps2->doFinalStatsAndSave(min,fname,fnameROOT);        
    }
    
    AmpsPdfFlexiFastCPV(AmpsPdfFlexiFast* amps1, AmpsPdfFlexiFast* amps2, AmpsPdfFlexiFast* ampsSum, 
                        MINT::FitParameter& r, MINT::FitParameter& delta,MINT::FitParameter& gamma,
                        MINT::FitParameter& tau, MINT::FitParameter& dGamma, MINT::FitParameter& dm, MINT::FitParameter& eff_tag, MINT::FitParameter& w ):
        _amps1(amps1),_amps2(amps2),_ampsSum(ampsSum),_r(r),_delta(delta),_gamma(gamma),_tau(tau),_dGamma(dGamma),_dm(dm),_eff_tag(eff_tag),_w(w),
        _intA(-1),_intAbar(-1),_intAAbar(-1)
         {;}
};

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


int ampFit(int step=0){    
  TRandom3 ranLux;
  NamedParameter<int> RandomSeed("RandomSeed", 0);
  int seed = RandomSeed + step;
  ranLux.SetSeed((int)seed);
  gRandom = &ranLux;
  
  FitAmplitude::AutogenerateFitFile();
  
  NamedParameter<string> InputFileName("InputFileName", (std::string) "");
  std::string inputFile = InputFileName;
  bool generateNew = (std::string) InputFileName == "";
  std::cout << "InputFileName: " << InputFileName << std::endl;

  NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
  NamedParameter<string> OutputRootFile("OutputRootFile", (std::string) "OutputRootFile.root", (char*) 0);
    
  NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
  NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
  NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");

  NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
  DalitzEventPattern pat(EventPattern.getVector());
  cout << " got event pattern: " << pat << endl;
  DalitzEventPattern patBar(EventPattern.getVector());
  patBar[0].antiThis();

  NamedParameter<int>  Nevents("Nevents", 10);
  NamedParameter<double>  pdf_max("pdf_max", 10);
  NamedParameter<double>  t_max("t_max", 10);
  NamedParameter<int>  do2DScan("do2DScan", 0);
    
  FitParameter  tau("tau");
  FitParameter  dGamma("dGamma");
  FitParameter  dm("dm");

  FitParameter  r("r");
  FitParameter  delta("delta");
  FitParameter  gamma("gamma");

  FitParameter  eff_tag("eff_tag");
  FitParameter  mistag("mistag");
    
  DalitzEventList eventListPhsp,eventList;
  DalitzEventList eventList_f, eventList_fbar;

  eventListPhsp.generatePhaseSpaceEvents(100000,pat);

  FitAmpSum fas_tmp((DalitzEventPattern)pat);
  fas_tmp.getVal(eventListPhsp[0]);
  fas_tmp.normalizeAmps(eventListPhsp);
  
  counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K*(1410)+");
  FitAmpSum fas(*List_1);
  FitAmpSum fasCC(*List_1);
    //fas.getVal(eventListPhsp[0]);
    //fasCC.getVal(eventListPhsp[0]);

  fasCC.CPConjugateSameFitParameters();
  fasCC.CConjugateFinalStateSameFitParameters();
    //fas.getVal(eventListPhsp[0]);
    //fasCC.getVal(eventListPhsp[0]);
  FitParameter r_K1_re("r_K1_Re",2,0,0.01);
  FitParameter r_K1_im("r_K1_Im",2,0,0.01); 
  counted_ptr<IReturnComplex> r_K1_plus = new CPV_amp_polar(r_K1_re,r_K1_im,1);
  counted_ptr<IReturnComplex> r_K1_minus = new CPV_amp_polar(r_K1_re,r_K1_im,-1);
  fas.multiply(r_K1_plus); 
  fasCC.multiply(r_K1_minus);
   cout <<  fas.getVal(eventListPhsp[0]) << endl;
   cout << fasCC.getVal(eventListPhsp[0])<< endl;
    cout <<  fas.getVal(eventListPhsp[1]) << endl;
    cout << fasCC.getVal(eventListPhsp[1])<< endl;  
  
 /* 
  AddScaledAmpsToList(fas_tmp, fas, fasCC, "K(1)(1400)+", r_K1_plus, r_K1_minus );

  FitParameter r_2_re("r_2_Re",2,0,0.01);
  FitParameter r_2_im("r_2_Im",2,0,0.01); 
  counted_ptr<IReturnComplex> r_2_plus = new CPV_amp(r_2_re,r_2_im,1);
  counted_ptr<IReturnComplex> r_2_minus = new CPV_amp(r_2_re,r_2_im,-1);
  AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,K+),rho(770)0(->pi+,pi-)", r_2_plus, r_2_minus );
  AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResS0(->Ds-,K+),rho(770)0(->pi+,pi-)", r_2_plus, r_2_minus );


  FitParameter r_3_re("r_3_Re",2,0,0.01);
  FitParameter r_3_im("r_3_Im",2,0,0.01); 
  counted_ptr<IReturnComplex> r_3_plus = new CPV_amp(r_3_re,r_3_im,1);
  counted_ptr<IReturnComplex> r_3_minus = new CPV_amp(r_3_re,r_3_im,-1);
  AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
  //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
*/
  coherenceFactor(fas,fasCC,(double)r, (double)delta,eventListPhsp);
  return 0;
  counted_ptr<FitAmpList> sumList = fas.GetCloneSameFitParameters();
  FitAmpSum fas_sum(*sumList);
  
  fas_sum.addAsList(fasCC,1.);
    
  AmpsPdfFlexiFast ampsSig(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
  AmpsPdfFlexiFast ampsSigCC(pat, &fasCC, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
  AmpsPdfFlexiFast ampsSum(pat, &fas_sum, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
    
  AmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSigCC,&ampsSum, r, delta, gamma, tau, dGamma, dm, eff_tag, mistag );
  
  /*
    SignalGenerator sg(pat,&fas);
    sg.FillEventList(eventList,1000);
    Neg2LL fcnt(ampsSig, eventList);
    
    Neg2LLSum neg2LLSumt(&fcnt);
    //neg2LLSum.addConstraints();
    
    Minimiser minit(&neg2LLSumt);
    minit.doFit();
    minit.printResultVsInput();
    
    ampsSig.doFinalStatsAndSave(&minit,("FitAmpResults_rand_"+anythingToString((int)RandomSeed)+ ".txt").c_str(),("fitFractions_"+anythingToString((int)RandomSeed)+ ".root").c_str());
    DalitzHistoSet datH = eventList.histoSet();
    DalitzHistoSet fitH = ampsSig.histoSet();
    datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_"+anythingToString((int)RandomSeed)+"_").c_str(),"eps");
    return 0;
    */
    
  if(generateNew){
  time_t startTime = time(0); 

  
  TFile* fileTD= TFile::Open("dummy.root","RECREATE");
  fileTD->cd();
  
  TTree *tree = new TTree("TD","TD");
  double t,dt,w;
  int q,f;
  tree->Branch("t",&t,"t/D");
  tree->Branch("dt",&dt,"dt/D");
  tree->Branch("q",&q,"q/I");
  tree->Branch("w",&w,"w/D");
  tree->Branch("f",&f,"f/I");
  
  FitAmpIncoherentSum fasGen((DalitzEventPattern)pat); 
  fasGen.getVal(eventListPhsp[0]);
  fasGen.normalizeAmps(eventListPhsp);
  SignalGenerator sg(pat,&fasGen);
  fasGen.print();
 
  //simple hit and miss
  for(int i = 0; i < Nevents; i++){
      while(true){
          t = ranLux.Exp(tau);
          dt = 0;
          const double q_rand = ranLux.Uniform();
          q = 0;
          if (q_rand < 1./3.) q = -1;
          if (q_rand > 2./3.) q = 1;
          w = mistag;
          if(i<=Nevents/2)f = 1; 
          else f = -1;

          counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
          DalitzEvent evt(evtPtr.get());
          double maxVal = evt.getGeneratorPdfRelativeToPhaseSpace()*exp(-fabs(t)/(tau))/(tau)*pdf_max;

          evt.setValueInVector(0, t);
          evt.setValueInVector(1, dt);
          evt.setValueInVector(2, q);
          evt.setValueInVector(3, w);
          evt.setValueInVector(4, f);
          
          const double pdfVal = pdf.un_normalised_noPs(evt);
          
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
              tree->Fill();
              break;
          }
      }
  }
    
    cout << " Generated " << Nevents << " Events. Took " << (time(0) - startTime)/60 << " mins "<< endl;
    
    eventList.saveAsNtuple(OutputRootFile);
        
    TFile* file= TFile::Open(((std::string) OutputRootFile).c_str(),"UPDATE");
    file->cd();
    tree->CloneTree()->Write();
    file->Close();
    fileTD->Close();

  }
  else{

        cout << "reading events from file " << inputFile << endl;
      
        TFile *_InputFile =  TFile::Open(inputFile.c_str());
        TTree* in_tree, *in_treeTD;
        in_tree=dynamic_cast<TTree*>(_InputFile->Get("DalitzEventList"));
        eventList.fromNtuple(in_tree,1);
        in_treeTD=dynamic_cast<TTree*>(_InputFile->Get("TD"));
        double t,dt,w;
        int q,f;
        in_treeTD->SetBranchAddress("t",&t);
        in_treeTD->SetBranchAddress("dt",&dt);
        in_treeTD->SetBranchAddress("q",&q);
        in_treeTD->SetBranchAddress("w",&w);
        in_treeTD->SetBranchAddress("f",&f);
      
        if(in_treeTD->GetEntries() != in_tree->GetEntries()){cout << "ERROR: Different number of DalitzEvents and time information " << endl; throw "crash"; }
      
        for (unsigned int i = 0; i < in_treeTD->GetEntries(); i++) {
            in_treeTD->GetEntry(i);
            
            eventList[i].setValueInVector(0, t);
            eventList[i].setValueInVector(1, dt);
            eventList[i].setValueInVector(2, q);
            eventList[i].setValueInVector(3, w);
            eventList[i].setValueInVector(4, f);
            
            if(f==1)eventList_f.Add(eventList[i]);
            else eventList_fbar.Add(eventList[i]);
        }
      
        cout << " I've got " << eventList.size() << " events." << endl;
        _InputFile->Close(); 
  }
        
    // Fit
    
    /*
    FitParameter  k("k",2,1,0.1);
    TimePdf_mod t_pdf(r,delta,gamma,k,tau,dGamma,dm,eff_tag,mistag );
    Neg2LL fcn_t(t_pdf, eventList_f);
    Neg2LL fcn_t_bar(t_pdf, eventList_fbar);
    Neg2LLSum neg2LLSum_t(&fcn_t,&fcn_t_bar);

    Minimiser mini_t(&neg2LLSum_t);
    mini_t.doFit();
    mini_t.printResultVsInput();
    */
    
    Neg2LL fcn(pdf, eventList_f);
    Neg2LL fcn_bar(pdf, eventList_fbar);

    Neg2LLSum neg2LLSum(&fcn,&fcn_bar);
    //neg2LLSum.addConstraints();
    
    Minimiser mini(&neg2LLSum);
    mini.doFit();
    mini.printResultVsInput();
    
    gDirectory->cd();
    TFile* paraFile = new TFile(((string)OutputDir+"pull_"+anythingToString((int)seed)+".root").c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;
    TTree* tree = new TTree("Coherence","Coherence");
    double r_val,delta_val,gamma_val,k_val,n2ll;
    TBranch* br_r = tree->Branch( "r", &r_val, "r_val/D" );
    TBranch* br_delta = tree->Branch( "delta", &delta_val, "delta_val/D" );
    TBranch* br_gamma = tree->Branch( "gamma", &gamma_val, "gamma_val/D" );
    TBranch* br_k = tree->Branch( "k", &k_val, "k_val/D" );

    TBranch* br_n2ll = tree->Branch( "n2ll", &n2ll, "n2ll/D" );
    TBranch* br_seed = tree->Branch( "seed", &seed, "seed/I" );
    
    MinuitParameterSet::getDefaultSet()->fillNtp(paraFile, ntp);
    ntp->AutoSave();
    
    n2ll = neg2LLSum.getVal();
    
    for(unsigned int i=0; i < MinuitParameterSet::getDefaultSet()->size(); i++){
        if(0 == MinuitParameterSet::getDefaultSet()->getParPtr(i)) continue;
        if(A_is_in_B("gamma",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))gamma_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();
    }
    
    vector<double> k_fit = coherenceFactor(fas,fasCC,(double)r, (double)delta,eventListPhsp);
    r_val = k_fit[0];
    k_val = k_fit[1];
    delta_val = k_fit[2];
    
    tree->Fill();
    
    paraFile->cd();
    tree->SetDirectory(paraFile);
    tree->Write();
    paraFile->Close();
    delete paraFile;
    
    pdf.doFinalStatsAndSaveForAmp12(&mini,((string)OutputDir+"FitAmpResults_rand_"+anythingToString((int)seed)+ ".txt").c_str(),((string)OutputDir+"fitFractions_"+anythingToString((int)seed)+ ".root").c_str());
    //DalitzHistoSet datH = eventList.histoSet();
    //DalitzHistoSet fitH = ampsSig.histoSet();
    //datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_"+anythingToString((int)RandomSeed)+"_").c_str(),"eps");
    return 0;
    cout << "Now plotting:" << endl;

    int nBinst = 100;
    int nBins = 60;
    vector<int> s123;
    s123.push_back(1);
    s123.push_back(2);
    s123.push_back(3);
    
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);
    
    vector<int> s134;
    s134.push_back(1);
    s134.push_back(3);
    s134.push_back(4);
    
    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);
    
    TH1D* h_t = new TH1D("",";t",nBinst,0,t_max);
    TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
    TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,25);
    TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
        
    for (int i=0; i<eventList.size(); i++) {
        h_t->Fill(eventList[i].getValueFromVector(0));
        s_Kpipi->Fill(eventList[i].sij(s234)/(GeV*GeV));
        s_Kpi->Fill(eventList[i].s(2,4)/(GeV*GeV));
        s_pipi->Fill(eventList[i].s(3,4)/(GeV*GeV));
        s_Dspipi->Fill(eventList[i].sij(s134)/(GeV*GeV));
        s_DsK->Fill(eventList[i].s(1,2)/(GeV*GeV));
    }    
        
    TH1D* h_t_fit = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_mp = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_0p = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_pp = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_mm = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_0m = new TH1D("",";t",nBinst,0,t_max);
    TH1D* h_t_fit_pm = new TH1D("",";t",nBinst,0,t_max);

    TH1D* s_Kpipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
    TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,25);
    TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    
    TH1D* s_Kpipi_fitBs = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
    TH1D* s_Kpi_fitBs = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_pipi_fitBs = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi_fitBs = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,25);
    TH1D* s_DsK_fitBs = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    
    TH1D* s_Kpipi_fitBsbar = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
    TH1D* s_Kpi_fitBsbar = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_pipi_fitBsbar = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi_fitBsbar = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,25);
    TH1D* s_DsK_fitBsbar = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    
    TH1D* s_Kpipi_fit_notag = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
    TH1D* s_Kpi_fit_notag = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_pipi_fit_notag = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi_fit_notag = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,25);
    TH1D* s_DsK_fit_notag = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    
    SignalGenerator sg(pat,&fas);
    sg.setWeighted();

    for(int i = 0; i < 2000000; i++){

        const double t = ranLux.Exp(tau);
        double dt = 0;
        double w = mistag;
        
        counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
        DalitzEvent evt(evtPtr.get());
        
        evt.setValueInVector(0, t);
        evt.setValueInVector(1, dt);
        evt.setValueInVector(3, w);
        
        evt.setValueInVector(2, -1);
        evt.setValueInVector(4, 1);
        double weight_mp = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));

        evt.setValueInVector(2, 0);
        evt.setValueInVector(4, 1);
        double weight_0p = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));
        
        evt.setValueInVector(2, 1);
        evt.setValueInVector(4, 1);
        double weight_pp = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));
        
        evt.setValueInVector(2, -1);
        evt.setValueInVector(4, -1);
        double weight_mm = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));
        
        evt.setValueInVector(2, 0);
        evt.setValueInVector(4, -1);
        double weight_0m = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));
        
        evt.setValueInVector(2, 1);
        evt.setValueInVector(4, -1);
        double weight_pm = pdf.un_normalised_noPs(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace()/exp(-fabs(t)/(tau));

        double weight = weight_mp + weight_0p + weight_pp + weight_mm + weight_0m + weight_pm ;
        
        h_t_fit->Fill(t,weight);
        h_t_fit_mp->Fill(t,weight_mp);
        h_t_fit_0p->Fill(t,weight_0p);
        h_t_fit_pp->Fill(t,weight_pp);
        h_t_fit_mm->Fill(t,weight_mm);
        h_t_fit_0m->Fill(t,weight_0m);
        h_t_fit_pm->Fill(t,weight_pm);
        
        double weight_B = fas.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
        double weight_Bbar = fasCC.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
        double weight_notag = weight_Bbar+ (double)r * (double)r * weight_Bbar;

        const complex<double> phase_diff = polar((double)r, (double)delta/360.*2*pi);
        double weight_coherence = std::abs(conj(fas.getVal(evt)) * fasCC.getVal(evt)*phase_diff) *evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
        
        s_Kpipi_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
        s_Kpi_fit->Fill(evt.s(2,4)/(GeV*GeV),weight);
        s_pipi_fit->Fill(evt.s(3,4)/(GeV*GeV),weight);
        s_Dspipi_fit->Fill(evt.sij(s134)/(GeV*GeV),weight);
        s_DsK_fit->Fill(evt.s(1,2)/(GeV*GeV),weight);
        
        s_Kpipi_fitBs->Fill(evt.sij(s234)/(GeV*GeV),weight_B);
        s_Kpi_fitBs->Fill(evt.s(2,4)/(GeV*GeV),weight_B);
        s_pipi_fitBs->Fill(evt.s(3,4)/(GeV*GeV),weight_B);
        s_Dspipi_fitBs->Fill(evt.sij(s134)/(GeV*GeV),weight_B);
        s_DsK_fitBs->Fill(evt.s(1,2)/(GeV*GeV),weight_B);

        s_Kpipi_fitBsbar->Fill(evt.sij(s234)/(GeV*GeV),weight_Bbar);
        s_Kpi_fitBsbar->Fill(evt.s(2,4)/(GeV*GeV),weight_Bbar);
        s_pipi_fitBsbar->Fill(evt.s(3,4)/(GeV*GeV),weight_Bbar);
        s_Dspipi_fitBsbar->Fill(evt.sij(s134)/(GeV*GeV),weight_Bbar);
        s_DsK_fitBsbar->Fill(evt.s(1,2)/(GeV*GeV),weight_Bbar);

        s_Kpipi_fit_notag->Fill(evt.sij(s234)/(GeV*GeV),weight_notag);
        s_Kpi_fit_notag->Fill(evt.s(2,4)/(GeV*GeV),weight_notag);
        s_pipi_fit_notag->Fill(evt.s(3,4)/(GeV*GeV),weight_notag);
        s_Dspipi_fit_notag->Fill(evt.sij(s134)/(GeV*GeV),weight_notag);
        s_DsK_fit_notag->Fill(evt.s(1,2)/(GeV*GeV),weight_notag);
        
    }
    
    TCanvas* c = new TCanvas();
    
    h_t->SetLineColor(kBlack);
    h_t->DrawNormalized("e1",1);
    h_t_fit->SetLineColor(kBlack);
    h_t_fit->SetLineWidth(3);
    h_t_fit->DrawNormalized("histcsame",1);
    
    h_t_fit_pm->SetLineColor(kBlue);
    h_t_fit_pm->SetLineWidth(3);
    h_t_fit_pm->DrawNormalized("histcsame",0.5*0.5*eff_tag);
    h_t_fit_mm->SetLineColor(kBlue);
    h_t_fit_mm->SetLineWidth(3);
    h_t_fit_mm->SetLineStyle(kDashed);
    h_t_fit_mm->DrawNormalized("histcsame",0.5*0.5*eff_tag);

    h_t_fit_pp->SetLineColor(kRed);
    h_t_fit_pp->SetLineWidth(3);
    h_t_fit_pp->DrawNormalized("histcsame",0.5*0.5*eff_tag);
    h_t_fit_mp->SetLineColor(kRed);
    h_t_fit_mp->SetLineWidth(3);
    h_t_fit_mp->SetLineStyle(kDashed);
    h_t_fit_mp->DrawNormalized("histcsame",0.5*0.5*eff_tag); 
    
    h_t_fit_0p->SetLineColor(kGreen);
    h_t_fit_0p->SetLineWidth(3);
    h_t_fit_0p->DrawNormalized("histcsame",0.5*0.5*(1.-eff_tag));
    h_t_fit_0m->SetLineColor(kGreen);
    h_t_fit_0m->SetLineWidth(3);
    h_t_fit_0m->SetLineStyle(kDashed);
    h_t_fit_0m->DrawNormalized("histcsame",0.5*0.5*(1.-eff_tag));
    
    c->Print("h_t.pdf");
    gPad->SetLogy(1);
    c->Print("h_t_log.pdf");
    gPad->SetLogy(0);

    
    s_Kpipi->SetLineColor(kBlack);
    s_Kpipi->DrawNormalized("e1",1);
    s_Kpipi_fit->SetLineColor(kBlack);
    s_Kpipi_fit->SetLineWidth(3);
    s_Kpipi_fit->DrawNormalized("histcsame",1);
    s_Kpipi_fitBs->SetLineColor(kRed);
    s_Kpipi_fitBs->SetLineWidth(3);
    s_Kpipi_fitBs->SetLineStyle(kDashed);
    s_Kpipi_fitBs->DrawNormalized("histcsame",1./(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Kpipi_fitBsbar->SetLineWidth(3);
    s_Kpipi_fitBsbar->SetLineColor(kBlue);
    s_Kpipi_fitBsbar->SetLineStyle(kDashed);
    s_Kpipi_fitBsbar->DrawNormalized("histcsame",k_fit[0]*k_fit[0]/(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Kpipi_fit_notag->SetLineWidth(3);
    s_Kpipi_fit_notag->SetLineColor(kMagenta);
    s_Kpipi_fit_notag->SetLineStyle(kDashed);
    s_Kpipi_fit_notag->DrawNormalized("histcsame",(1.-eff_tag));
    c->Print("s_Kpipi.pdf");
    
    s_Kpi->SetLineColor(kBlack);
    s_Kpi->DrawNormalized("e1",1);
    s_Kpi_fit->SetLineColor(kBlack);
    s_Kpi_fit->SetLineWidth(3);
    s_Kpi_fit->DrawNormalized("histcsame",1);
    s_Kpi_fitBs->SetLineColor(kRed);
    s_Kpi_fitBs->SetLineWidth(3);
    s_Kpi_fitBs->SetLineStyle(kDashed);
    s_Kpi_fitBs->DrawNormalized("histcsame",1./(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Kpi_fitBsbar->SetLineWidth(3);
    s_Kpi_fitBsbar->SetLineColor(kBlue);
    s_Kpi_fitBsbar->SetLineStyle(kDashed);
    s_Kpi_fitBsbar->DrawNormalized("histcsame",k_fit[0]*k_fit[0]/(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Kpi_fit_notag->SetLineWidth(3);
    s_Kpi_fit_notag->SetLineColor(kMagenta);
    s_Kpi_fit_notag->SetLineStyle(kDashed);
    s_Kpi_fit_notag->DrawNormalized("histcsame",(1.-eff_tag));
    c->Print("s_Kpi.pdf");
    
    s_pipi->SetLineColor(kBlack);
    s_pipi->DrawNormalized("e1",1);
    s_pipi_fit->SetLineColor(kBlack);
    s_pipi_fit->SetLineWidth(3);
    s_pipi_fit->DrawNormalized("histcsame",1);
    s_pipi_fitBs->SetLineColor(kRed);
    s_pipi_fitBs->SetLineWidth(3);
    s_pipi_fitBs->SetLineStyle(kDashed);
    s_pipi_fitBs->DrawNormalized("histcsame",1./(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_pipi_fitBsbar->SetLineWidth(3);
    s_pipi_fitBsbar->SetLineColor(kBlue);
    s_pipi_fitBsbar->SetLineStyle(kDashed);
    s_pipi_fitBsbar->DrawNormalized("histcsame",k_fit[0]*k_fit[0]/(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_pipi_fit_notag->SetLineWidth(3);
    s_pipi_fit_notag->SetLineColor(kMagenta);
    s_pipi_fit_notag->SetLineStyle(kDashed);
    s_pipi_fit_notag->DrawNormalized("histcsame",(1.-eff_tag));
    c->Print("s_pipi.pdf");
    
    s_Dspipi->SetLineColor(kBlack);
    s_Dspipi->DrawNormalized("e1",1);
    s_Dspipi_fit->SetLineColor(kBlack);
    s_Dspipi_fit->SetLineWidth(3);
    s_Dspipi_fit->DrawNormalized("histcsame",1);
    s_Dspipi_fitBs->SetLineColor(kRed);
    s_Dspipi_fitBs->SetLineWidth(3);
    s_Dspipi_fitBs->SetLineStyle(kDashed);
    s_Dspipi_fitBs->DrawNormalized("histcsame",1./(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Dspipi_fitBsbar->SetLineWidth(3);
    s_Dspipi_fitBsbar->SetLineColor(kBlue);
    s_Dspipi_fitBsbar->SetLineStyle(kDashed);
    s_Dspipi_fitBsbar->DrawNormalized("histcsame",k_fit[0]*k_fit[0]/(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_Dspipi_fit_notag->SetLineWidth(3);
    s_Dspipi_fit_notag->SetLineColor(kMagenta);
    s_Dspipi_fit_notag->SetLineStyle(kDashed);
    s_Dspipi_fit_notag->DrawNormalized("histcsame",(1.-eff_tag));
    c->Print("s_Dspipi.pdf");
    
    s_DsK->SetLineColor(kBlack);
    s_DsK->DrawNormalized("e1",1);
    s_DsK_fit->SetLineColor(kBlack);
    s_DsK_fit->SetLineWidth(3);
    s_DsK_fit->DrawNormalized("histcsame",1);
    s_DsK_fitBs->SetLineColor(kRed);
    s_DsK_fitBs->SetLineWidth(3);
    s_DsK_fitBs->SetLineStyle(kDashed);
    s_DsK_fitBs->DrawNormalized("histcsame",1./(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_DsK_fitBsbar->SetLineWidth(3);
    s_DsK_fitBsbar->SetLineColor(kBlue);
    s_DsK_fitBsbar->SetLineStyle(kDashed);
    s_DsK_fitBsbar->DrawNormalized("histcsame",k_fit[0]*k_fit[0]/(1.+k_fit[0]*k_fit[0])*eff_tag);
    s_DsK_fit_notag->SetLineWidth(3);
    s_DsK_fit_notag->SetLineColor(kMagenta);
    s_DsK_fit_notag->SetLineStyle(kDashed);
    s_DsK_fit_notag->DrawNormalized("histcsame",(1.-eff_tag));
    c->Print("s_DsK.pdf");
    
        
    if(do2DScan == 1){
        cout << "Now doing 2D scan:" << endl;
        int scanBins=20;
        double scanMin=0, scanMax=360;
        double nSigmaZoom = 2;
        double scanMinGammaZoom=min(gamma.meanInit(), gamma.mean()) - nSigmaZoom*gamma.err();
        double scanMaxGammaZoom=max(gamma.meanInit(), gamma.mean()) + nSigmaZoom*gamma.err();
        double scanMinDeltaZoom=min(delta.meanInit(), delta.mean()) - nSigmaZoom*delta.err();
        double scanMaxDeltaZoom=max(delta.meanInit(), delta.mean()) + nSigmaZoom*delta.err();
        double gammaZoomRange = scanMaxGammaZoom - scanMinGammaZoom;
        double deltaZoomRange = scanMaxDeltaZoom - scanMinDeltaZoom;
        
        TFile* scanFile = new TFile("scan.root", "RECREATE");
        TH2D* scanHisto = new TH2D("scan", "; #gamma [deg]; #delta [deg]", scanBins,  scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanHistoP = new TH2D("scanP", "; #gamma [deg]; #delta [deg]", scanBins,  scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanHistoM = new TH2D("scanM" , "; #gamma [deg]; #delta [deg]", scanBins, scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanZoomHisto = new TH2D("scanZoom", "; #gamma [deg]; #delta [deg]", scanBins, scanMinGammaZoom, scanMaxGammaZoom, scanBins, scanMinDeltaZoom, scanMaxDeltaZoom);
        
        double scanMinLL=-9999;
        double scanMinLLP=-9999;
        double scanMinLLM=-9999;
        double scanMinLLZ=-9999;
        
        for(int i=0; i < scanBins; i++){
            double gamma_value = ((double)i+0.5)*360/((double)scanBins);
            gamma.setCurrentFitVal(gamma_value);
            for(int j=0; j < scanBins; j++){
                double delta_value = ((double)j+0.5)*360/((double)scanBins);
                delta.setCurrentFitVal(delta_value);
                
                double v = neg2LLSum.getNewVal();  
                if( (i==0 && j==0) || v < scanMinLL) scanMinLL=v;
                scanHisto->Fill(gamma_value, delta_value, v);
                
                double vP = fcn.getNewVal();  
                if( (i==0 && j==0) || vP < scanMinLLP) scanMinLLP=vP;
                scanHistoP->Fill(gamma_value, delta_value, vP);
                
                double vM = fcn_bar.getNewVal();  
                if( (i==0 && j==0) || vM < scanMinLLM) scanMinLLM=vM;
                scanHistoM->Fill(gamma_value, delta_value, vM);
            }
        }
        for(int i=0; i < scanBins; i++){
            double gamma_value = scanMinGammaZoom + ((double)i+0.5) * gammaZoomRange/((double)scanBins);
            gamma.setCurrentFitVal(gamma_value);
            for(int j=0; j < scanBins; j++){
                double delta_value = scanMinDeltaZoom + ((double)j+0.5) * deltaZoomRange/((double)scanBins);
                delta.setCurrentFitVal(delta_value);
                double v = neg2LLSum.getNewVal();
                
                if( (i==0 && j==0) || v < scanMinLLZ) scanMinLLZ=v;
                
                scanZoomHisto->Fill(gamma_value, delta_value, v);
            }
        }
        
        for(int i=0; i < scanBins; i++){
            double gamma_value = ((double)i+0.5)*360/((double)scanBins);
            for(int j=0; j < scanBins; j++){
                double delta_value = ((double)j+0.5)*360/((double)scanBins);
                scanHisto->Fill(gamma_value, delta_value, -scanMinLL);
                scanHistoP->Fill(gamma_value, delta_value, -scanMinLLP);
                scanHistoM->Fill(gamma_value, delta_value, -scanMinLLM);
            }
        }
        for(int i=0; i < scanBins; i++){
            double gamma_value = scanMinGammaZoom + ((double)i+0.5) * gammaZoomRange/((double)scanBins);
            for(int j=0; j < scanBins; j++){
                double delta_value = scanMinDeltaZoom + ((double)j+0.5) * deltaZoomRange/((double)scanBins);
                scanZoomHisto->Fill(gamma_value, delta_value, -scanMinLLZ);
            }
        }
        scanFile->cd();
        scanHisto->Write();
        scanHistoP->Write();
        scanHistoM->Write();
        scanZoomHisto->Write();
        scanFile->Close();
        
        cout<< "done 2-D scan" << endl;
    }

        
  return 0;
}

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

void makeIntegratorFile(){

    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventListPhsp,eventList;
    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.getVal(eventListPhsp[0]);
    fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);

    sg.FillEventList(eventList, IntegratorEvents);
    eventList.saveAsNtuple(IntegratorEventFile);
    
    return;
}

void calculateCoherence(int step=0){

    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    TFile* paraFile = new TFile(((string)OutputDir+"coherence_"+anythingToString((int)seed)+".root").c_str(), "RECREATE");
    paraFile->cd();
    TTree* tree = new TTree("coherence","coherence");
    double r_val,k,delta_val;
    double x,y;
    TBranch* br_r = tree->Branch( "r", &r_val, "r_val/D" );
    TBranch* br_k = tree->Branch( "k", &k, "k/D" );
    TBranch* br_delta = tree->Branch( "delta", &delta_val, "delta_val/D" );
    TBranch* br_x = tree->Branch( "x", &x, "x/D" );
    TBranch* br_y = tree->Branch( "y", &y, "y/D" );
    
    FitParameter  r("r");
    FitParameter  delta("delta");

    DalitzEventList eventListPhsp;
    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    // Bs -> f Amplitude
    FitAmpSum fas_tmp((DalitzEventPattern)pat);
    fas_tmp.getVal(eventListPhsp[0]);
    fas_tmp.normalizeAmps(eventListPhsp);
        
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1270)+");
    FitAmpSum fas(*List_1);
    FitAmpSum fasCC(*List_1);
    fasCC.CPConjugateSameFitParameters();
    fasCC.CConjugateFinalStateSameFitParameters();
    x = ranLux.Uniform(-1,1);
    y = ranLux.Uniform(-pi,pi);
    //y = sqrt(1. - x*x);
    //if(ranLux.Uniform() > 0.5) y= -y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl << endl;
    FitParameter r_K1("r_K1",2,x,0.1);
    FitParameter delta_K1("delta_K1",2,y,0.1); 
    counted_ptr<IReturnComplex> r_K1_plus = new CPV_amp_polar(r_K1,delta_K1,1);
    counted_ptr<IReturnComplex> r_K1_minus = new CPV_amp_polar(r_K1,delta_K1,-1);
    fas.multiply(r_K1_plus); 
    fasCC.multiply(r_K1_minus); 
    
    counted_ptr<FitAmpList> List_2 = fas_tmp.GetCloneOfSubsetSameFitParameters("NonResV0(->Ds-,K+),rho(770)0(->pi+,pi-)");
    FitAmpSum fas_2(*List_2);
    FitAmpSum fasCC_2(*List_2);
    fasCC_2.CPConjugateSameFitParameters();
    fasCC_2.CConjugateFinalStateSameFitParameters();
    x = ranLux.Uniform(-1,1);
    y = ranLux.Uniform(-pi,pi);
    //y = sqrt(1. - x*x);
    //if(ranLux.Uniform() > 0.5) y= -y;
    cout << "x = " << x << endl;
    cout << "y = " << y << endl << endl;
    FitParameter r_2("r_2",2,x,0.1);
    FitParameter delta_2("delta_2",2,y,0.1); 
    counted_ptr<IReturnComplex> r_2_plus = new CPV_amp_polar(r_2,delta_2,1);
    counted_ptr<IReturnComplex> r_2_minus = new CPV_amp_polar(r_2,delta_2,-1);
    fas_2.multiply(r_2_plus); 
    fasCC_2.multiply(r_2_minus); 
    fas.addAsList(fas_2,1.);
    fasCC.addAsList(fasCC_2,1.);
    
    vector<double> kappa = coherenceFactor(fas,fasCC,(double)r, (double)delta,eventListPhsp);
    
    r_val = kappa[0];
    k = kappa[1];
    delta_val = kappa[2];
    
    tree->Fill();
    paraFile->cd();
    tree->SetDirectory(paraFile);
    tree->Write();
    paraFile->Close();
    delete paraFile;
    
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


  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gStyle->SetOptStat(0);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleYSize(0.05);
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelOffset(0.01,"X");
  gStyle->SetTitleOffset(0.8,"X");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //TD_test();

  //makeIntegratorFile();
  timeFit_mod();
  //ampFit(atoi(argv[1]));
  
  //calculateCoherence(atoi(argv[1]));
  //for(int i=0;i<200;i++){
     // calculateCoherence(i);
  //}
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
