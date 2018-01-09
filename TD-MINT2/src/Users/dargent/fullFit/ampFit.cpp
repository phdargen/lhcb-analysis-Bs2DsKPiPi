// Full td amplitude fit
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
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooUniform.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "Mint/TimePdfMaster.h"
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
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include "TProfile.h"

using namespace std;
using namespace RooFit ;
using namespace MINT;

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
                     , std::string fname =  "SignalIntegrationEvents.root", bool generateNew = false, bool genMoreEvents = false
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
        //bool generateNew = ((string)_integratorSource == (string)"new");
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
        //std::complex<double> result= polar((double) sqrt( 1.+  _r * (double) _sign),(double) (_delta/360.*2.*pi * (double) _sign) ); 
	    std::complex<double> result= polar((double) ( 1.+  _r * (double) _sign),(double) (_delta/360.*2.*pi * (double) _sign) ); 
        return result;
    }
};

void AddScaledAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, FitAmpSum& fasCC, std::string name, counted_ptr<IReturnComplex>& r_plus, counted_ptr<IReturnComplex>& r_minus){
    
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    FitAmpSum fasCC_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas_2.multiply(r_plus); 
    fasCC_2.multiply(r_minus); 
    fas.addAsList(fas_2,1.);
    fasCC.addAsList(fasCC_2,1.);
}

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
        const std::complex<double> amp = fas.getVal(eventList[i]) ;
        const std::complex<double> amp_bar = fas_bar.getVal(eventList[i])*phase_diff ;
        valK += amp_bar*conj(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
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

std::vector<double> coherenceFactorCC(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;    
    //fas.CPConjugateSameFitParameters();
    //fas_bar.CPConjugateSameFitParameters();

    //fas.print();
    //fas_bar.print();
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
    
        DalitzEvent evt(eventList[i]);
        evt.P_conjugateYourself();
    
        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
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

// Full time-dependent PDF
class FullAmpsPdfFlexiFastCPV : public MINT::PdfBase<IDalitzEvent>
, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps1;
    AmpsPdfFlexiFast* _amps2;
    AmpsPdfFlexiFast* _ampsSum;
    AmpsPdfFlexiFast* _ampsSum_CP;

    double _intA;
    double _intAbar;
    double _intA_CP;
    double _intAbar_CP;
    complex<double> _intAAbar;   
    complex<double> _intAAbar_CP;   
 
    // Fit parameters
    FitParameter& _r;
    FitParameter& _delta;
    FitParameter& _gamma;
    
    // Time pdf master
    TimePdfMaster _timePdfMaster;

    // limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;

public:
    void parametersChanged(){
        _ampsSum->parametersChanged();
        _ampsSum_CP->parametersChanged();

        _intA = (_ampsSum->ComplexIntegralForTags(1,1)).real();
        _intAbar = (_ampsSum->ComplexIntegralForTags(-1,-1)).real();        
        _intAAbar = _ampsSum->ComplexIntegralForTags(1,-1);

        _intA_CP = (_ampsSum_CP->ComplexIntegralForTags(1,1)).real();
        _intAbar_CP = (_ampsSum_CP->ComplexIntegralForTags(-1,-1)).real();   
        _intAAbar_CP = _ampsSum_CP->ComplexIntegralForTags(1,-1);
    }
    void beginFit(){
        _ampsSum->beginFit();
        _ampsSum_CP->beginFit();
        _timePdfMaster.listFitParDependencies();
        printIntegralVals();
    }
    void endFit(){
        printIntegralVals();
        _ampsSum->endFit();
        _ampsSum_CP->endFit();
    }
    
    void printIntegralVals(){
        cout << "intSum = " << _ampsSum->getIntegralValue() << endl;
        cout << "intA = " << _intA << endl;
        cout << "intAbar = " << _intAbar << endl;
        cout << "intAAbar = " << _intAAbar << endl;

        cout << "intSum_CP = " << _ampsSum_CP->getIntegralValue() << endl;
        cout << "intA_CP = " << _intA_CP << endl;
        cout << "intAbar_CP = " << _intAbar_CP << endl;
        cout << "intAAbar_CP = " << _intAAbar_CP << endl;
    }

    /*    
    inline double getValForGeneration(IDalitzEvent& evt){
        const double t = (double) evt.getValueFromVector(0);
        const double dt = (double) evt.getValueFromVector(1);
        const double q = static_cast<double>((int)evt.getValueFromVector(2)); 
        const double w = (double) evt.getValueFromVector(3);
        const double f = static_cast<double>((int)evt.getValueFromVector(4));
        
        if(t < _min_TAU || t > _max_TAU )return 0.;        
        _r_t->setVal(t);
        _r_dt->setVal(dt);
        _r_q->setIndex(q);
        _r_mistag->setVal(w);
        _r_f->setIndex(f);
                
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma*f)/360.*2*pi);
        
        const std::complex<double> amp = _amps1->ComplexVal_un_normalised_noPs(evt) ;
        const std::complex<double> amp_bar = _amps2->ComplexVal_un_normalised_noPs(evt) * phase_diff;
        
        const double val =  exp(-fabs(t)/(double)_tau) *
         (
           (norm(amp) + norm(amp_bar)) *cosh((double)_dGamma/2.*t) * _cosh_coeff->evaluate()
          +(norm(amp) - norm(amp_bar)) *cos((double)_dm*t) * _cos_coeff->evaluate()
          +real(amp_bar*conj(amp)) *sinh((double)_dGamma/2.*t) * _sinh_coeff->evaluate()
          +imag(amp_bar*conj(amp)) *sin((double)_dm*t) * _sin_coeff->evaluate()
         ) *_pdf_sigma_t->getVal()*_spline->getVal();
        
        return val;
    }
 
    inline void getSmearedTime(double& t, double& dt, TRandom3& r){
        while(true){
                t = r.Exp(_tau);
                dt = r.Uniform(0.,0.25);
                t = t + r.Gaus(0,_r_scale_sigma_dt->getVal()*dt);
            if(_min_TAU< t && t < _max_TAU) return; 
        }
    }
*/
    inline double un_normalised_noPs(IDalitzEvent& evt){

        const double t = (double) evt.getValueFromVector(0);
        if(t < _min_TAU || t > _max_TAU )return 0.;
        _timePdfMaster.setAllObservablesAndFitParameters(evt);
                
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_bar = polar((double)r,((double) _delta +(double)_gamma)/360.*2*pi);

        const std::complex<double> amp = _amps1->ComplexVal_un_normalised_noPs(evt) ;
        const std::complex<double> amp_bar = _amps2->ComplexVal_un_normalised_noPs(evt) ;
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster.setCP_coeff(
        		(norm(amp) - norm(amp_bar))/(norm(amp) + norm(amp_bar)),
        		-(norm(amp) - norm(amp_bar))/(norm(amp) + norm(amp_bar)),
        		(-2.* real(amp_bar*conj(amp) * phase_diff) )/ (norm(amp) + norm(amp_bar)),
        		(-2.* real(amp_bar*conj(amp) * phase_diff_bar) )/ (norm(amp) + norm(amp_bar)),
        		(2.* imag(amp_bar*conj(amp) * phase_diff) )/ (norm(amp) + norm(amp_bar)),
        		(-2. * imag(amp_bar*conj(amp) * phase_diff_bar) )/ (norm(amp) + norm(amp_bar)) );
        
        const double val = 
	    (norm(amp) + norm(amp_bar)) *
             ( _timePdfMaster.get_cosh_term_Val(evt)
             +  _timePdfMaster.get_cos_term_Val(evt)
             +  _timePdfMaster.get_sinh_term_Val(evt)
             +  _timePdfMaster.get_sin_term_Val(evt)
             ) * _timePdfMaster.get_marginalPdfs_Val(evt);;
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double val = un_normalised_noPs(evt);
        double r = (double)_r; // * sqrt(_intA/_intAbar);
    
        /*
        if(_intA == -1 ){
            cout << "AmpsPdfFlexiFastCPV:: _norm = -1, should not have happened." << endl;
            throw "can't deal with that";
        }
	*/

        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_bar = polar((double)r,((double) _delta +(double)_gamma)/360.*2*pi);

        const complex<double> int_interference =  phase_diff * _intAAbar ;
        const complex<double> int_interference_bar = phase_diff_bar * _intAAbar_CP ;

        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster.setCP_coeff(
        		(_intA - r* r * _intAbar)/(_intA + r* r * _intAbar),
        		-(_intA_CP - r* r * _intAbar_CP)/(_intA + r* r * _intAbar) ,  /// sign ?
        		(- int_interference.real() )/ (_intA + r* r * _intAbar),
        		(- int_interference_bar.real() )/ (_intA + r* r * _intAbar),
        		(int_interference.imag() )/ (_intA + r* r * _intAbar),
        		(- int_interference_bar.imag() )/ (_intA + r* r * _intAbar) ) ;
        
        double norm = (_intA + r* r * _intAbar)
        *(     _timePdfMaster.get_cosh_term_Integral(evt)
            +  _timePdfMaster.get_cos_term_Integral(evt)
            +  _timePdfMaster.get_sinh_term_Integral(evt)
            +  _timePdfMaster.get_sin_term_Integral(evt) );

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

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_SS(evt);
    }
    
    virtual DalitzHistoSet histoSet(){return _ampsSum->histoSet();}
    
    void doFinalStatsAndSaveForAmp12(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults", const std::string& fnameROOT="fitFractions"){
        _amps1->redoIntegrator();
        _amps2->redoIntegrator();
        _amps1->doFinalStatsAndSave(min,((string)fname+".txt").c_str(),((string)fnameROOT+".root").c_str());
        _amps2->doFinalStatsAndSave(min,((string)fname+"_CC.txt").c_str(),((string)fnameROOT+"_CC.root").c_str());        
    }
    
    FullAmpsPdfFlexiFastCPV(
		AmpsPdfFlexiFast* amps1, AmpsPdfFlexiFast* amps2, AmpsPdfFlexiFast* ampsSum, AmpsPdfFlexiFast* ampsSum_CP, 
                MINT::FitParameter& r,MINT::FitParameter& delta, MINT::FitParameter& gamma
                ):
    _amps1(amps1),_amps2(amps2),_ampsSum(ampsSum),_ampsSum_CP(ampsSum_CP),
    _intA(-1),_intAbar(-1),_intAAbar(-1),_intAAbar_CP(-1),
    _r(r),_delta(delta),_gamma(gamma),
    _min_TAU("min_TAU", 0.4), _max_TAU("max_TAU", 10.)
    {
    }
};


void ampFit(int step=0){

    /// Options
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    
    TString prefix = "";

    // Generate list of amplitudes
    FitAmplitude::AutogenerateFitFile();
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    NamedParameter<double> w_max("w_max", 0.5);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  nBins("nBins", 40);
    NamedParameter<int>  nBinst("nBinst", 20);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 10);

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;
    NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
    NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");

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

    /// Define PDF 
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);

    /// Define amplitude model
    DalitzEventList eventListPhsp,eventList;
    DalitzEventList eventList_f, eventList_f_bar;
    
    eventListPhsp.generatePhaseSpaceEvents(200000,pat);
    
    FitAmpSum fas_tmp((DalitzEventPattern)pat);
    fas_tmp.getVal(eventListPhsp[0]);
    fas_tmp.normalizeAmps(eventListPhsp);
    
    //FitAmpSum fas_tmpCC(fas_tmp);
    //fas_tmpCC.CPConjugateSameFitParameters();
    //fas_tmpCC.CConjugateSameFitParameters();
    //fas_tmpCC.CConjugateInitialStateSameFitParameters();
    
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1270)+");
    FitAmpSum fas(*List_1);
    FitAmpSum fasCC(*List_1);
    //fasCC.CPConjugateSameFitParameters();
    //fasCC.CConjugateFinalStateSameFitParameters();
    FitParameter r_K1_re("r_K1_Re",2,0,0.01);
    FitParameter r_K1_im("r_K1_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_K1_plus = new CPV_amp_polar(r_K1_re,r_K1_im,1);
    counted_ptr<IReturnComplex> r_K1_minus = new CPV_amp_polar(r_K1_re,r_K1_im,-1);
    fas.multiply(r_K1_plus); 
    fasCC.multiply(r_K1_minus);
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "K(1)(1400)+", r_K1_plus, r_K1_minus );
    
    FitParameter r_2_re("r_2_Re",2,0,0.01);
    FitParameter r_2_im("r_2_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_2_plus = new CPV_amp_polar(r_2_re,r_2_im,1);
    counted_ptr<IReturnComplex> r_2_minus = new CPV_amp_polar(r_2_re,r_2_im,-1);
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "K*(1410)+", r_2_plus, r_2_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "K(1460)+(->K*(892)0(->K+,pi-),pi+),Ds-", r_2_plus, r_2_minus );

    FitParameter r_3_re("r_3_Re",2,0,0.01);
    FitParameter r_3_im("r_3_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_3_plus = new CPV_amp_polar(r_3_re,r_3_im,1);
    counted_ptr<IReturnComplex> r_3_minus = new CPV_amp_polar(r_3_re,r_3_im,-1);
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "BgSpinZeroBs0->NonResS0(->Ds-,K+),NonResS0(->pi+,pi-)", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "K(1460)+(->K*(892)0(->K+,pi-),pi+),Ds-", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,K+),sigma10(->pi+,pi-)", r_3_plus, r_3_minus );
     
    FitParameter r_4_re("r_4_Re",2,0,0.01);
    FitParameter r_4_im("r_4_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_4_plus = new CPV_amp_polar(r_4_re,r_4_im,1);
    counted_ptr<IReturnComplex> r_4_minus = new CPV_amp_polar(r_4_re,r_4_im,-1);
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResA0(->sigma10(->pi+,pi-),Ds-)", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResA0(->rho(770)0(->pi+,pi-),Ds-),K+", r_4_plus, r_4_minus );
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "BgSpinZeroBs0->NonResS0(->Ds-,K+),NonResS0(->pi+,pi-)", r_4_plus, r_4_minus );
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,K+),sigma10(->pi+,pi-)", r_4_plus, r_4_minus );
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_4_plus, r_4_minus );
    AddScaledAmpsToList(fas_tmp, fas, fasCC, "NonResV0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_4_plus, r_4_minus );
   

    /// Define B -> f amplitude        
    fas.setTag(1);
    /// Define Bbar -> f amplitude
    fasCC.setTag(-1);
    
    vector<double> k_gen = coherenceFactor(fas,fasCC,(double)r, (double)delta,eventListPhsp);
    
    counted_ptr<FitAmpList> sumList = fas.GetCloneSameFitParameters();
    FitAmpSum fas_sum(*sumList);
    fas_sum.addAsList(fasCC,1.);
    fas_sum.getVal(eventListPhsp[0]);
    
    AmpsPdfFlexiFast ampsSig(pat, &fas, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSigCC(pat, &fasCC, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSum(pat, &fas_sum, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSum_CP(pat, &fas_sum, 0, integPrecision,integMethod, (std::string) integratorEventFile.ReplaceAll(".root","_CP.root"));

    /// Make full time-dependent PDF
    FullAmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSigCC,&ampsSum,&ampsSum_CP, r,delta,gamma);

    /// Load data
    double t,dt;
    int f;
    int q_OS;
    Short_t q_SS;
    double eta_OS;
    Float_t eta_SS;
    double sw;
    int year,Ds_finalState;

    double K[4]; 
    double pip[4]; 
    double pim[4]; 
    double Ds_Kp[4],Ds_Km[4],Ds_pim[4];
    double mB;
    
    TChain* tree=new TChain("DecayTree");
    tree->Add(((string)InputDir+"Data/"+(string)channel+".root").c_str());
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("N_Bs_sw",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("*DEC",1);
    tree->SetBranchStatus("*PROB",1);
    tree->SetBranchStatus("*OS",1);
    tree->SetBranchStatus("*TAU*",1);
    tree->SetBranchStatus("*ID*",1);
    tree->SetBranchStatus("weight",1);
    tree->SetBranchStatus("Bs_DTF_MM",1);
    tree->SetBranchStatus("BsDTF_*P*",1);

    tree->SetBranchAddress("Bs_DTF_TAU",&t);
    tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
    tree->SetBranchAddress("Ds_ID",&f);
    tree->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
    tree->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&eta_OS);
    tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
    tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&eta_SS);
    tree->SetBranchAddress("N_Bs_sw",&sw);
    tree->SetBranchAddress("year",&year);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);    
    tree->SetBranchAddress("Bs_DTF_MM",&mB);
    
    tree->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
    tree->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
    tree->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]); 
    tree->SetBranchAddress("BsDTF_Kplus_PE",&K[3]); 
	
    tree->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
    tree->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
    tree->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]); 
    tree->SetBranchAddress("BsDTF_piplus_PE",&pip[3]); 
	
    tree->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
    tree->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
    tree->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]); 
    tree->SetBranchAddress("BsDTF_piminus_PE",&pim[3]); 
	
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PX",&Ds_Kp[0]);
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PY",&Ds_Kp[1]);
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PZ",&Ds_Kp[2]); 
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PE",&Ds_Kp[3]); 
    
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PX",&Ds_Km[0]);
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PY",&Ds_Km[1]);
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PZ",&Ds_Km[2]); 
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PE",&Ds_Km[3]); 

    tree->SetBranchAddress("BsDTF_Ds_piminus_PX",&Ds_pim[0]);
    tree->SetBranchAddress("BsDTF_Ds_piminus_PY",&Ds_pim[1]);
    tree->SetBranchAddress("BsDTF_Ds_piminus_PZ",&Ds_pim[2]); 
    tree->SetBranchAddress("BsDTF_Ds_piminus_PE",&Ds_pim[3]); 

    TRandom3 rndm;
    int badEvents = 0;

    for(int i=0; i< tree->GetEntries(); i++)
    {	
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
        tree->GetEntry(i);

        if(t < min_TAU || t > max_TAU )continue;
        if( dt < 0 || dt > 0.1 )continue;

        double sign = 1.;
        if(f > 0) sign = -1.;

        TLorentzVector K_p(sign*K[0],sign*K[1],sign*K[2],K[3]);
        TLorentzVector pip_p(sign*pip[0],sign*pip[1],sign*pip[2],pip[3]);
        TLorentzVector pim_p(sign*pim[0],sign*pim[1],sign*pim[2],pim[3]);
        TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
        TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
        TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
        TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
        TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        // array of vectors
        vector<TLorentzVector> vectorOfvectors; 

        if((string)channel=="norm"){
            TLorentzVector pip1_p, pip2_p;
            if(rndm.Rndm()<0.5) { 
                pip1_p = K_p;
                pip2_p = pip_p;
            } 
            else { 
                pip1_p = pip_p;
                pip2_p = K_p;
            } 
            K_p = pip1_p;
            pip_p = pip2_p;
        }

        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
        vectorOfvectors.push_back(pip_p*MeV);
        vectorOfvectors.push_back(pim_p*MeV);

        DalitzEvent evt(pat, vectorOfvectors);

        if(!(evt.phaseSpace() > 0.)){
            //cout << "evt " << i << " 0 phsp " << endl << evt << endl;
            badEvents++;
            continue;
        }
        if(TMath::IsNaN(norm(fas.getVal(evt)))){
            //cout << "evt " << i << " isNaN " << endl << evt << endl;
            badEvents++;
            continue;
        }

        evt.setWeight(sw);
        evt.setValueInVector(0, t);
        evt.setValueInVector(1, dt);
        if(f<0)evt.setValueInVector(2, 1);
        else if(f > 0)evt.setValueInVector(2, -1);
        else {
            cout << "ERROR:: Undefined final state";  
            throw "ERROR";
        }
        evt.setValueInVector(3, q_OS);
        evt.setValueInVector(4, eta_OS);
        evt.setValueInVector(5, q_SS);
        evt.setValueInVector(6, eta_SS);
        eventList.Add(evt);
        if(evt.getValueFromVector(2) == 1)eventList_f.Add(evt);
        else eventList_f_bar.Add(evt);
    }

   cout << endl << "bad events " << badEvents << " ( " << badEvents/(double) tree->GetEntries() * 100. << " %)" << endl << endl;

    /// Fit with MINT Pdf
    Neg2LL neg2LL(pdf, eventList);    
    //cout << "tau = " << endl << tau.mean() << endl <<  tau.blindedMean() << endl;
    //neg2LL.getVal();    
    Minimiser mini(&neg2LL);
    mini.doFit();
    mini.printResultVsInput();

     pdf.doFinalStatsAndSaveForAmp12(&mini,((string)OutputDir+"FitAmpResults_rand_"+anythingToString((int)seed)).c_str(),((string)OutputDir+"fitFractions_"+anythingToString((int)seed)).c_str());
    //cout << "tau = " << endl << tau.mean() << endl <<  tau.blindedMean() << endl;
    
    vector<double> k_fit = coherenceFactor(fas,fasCC,(double)r, (double)delta,eventListPhsp);

    /// Data histograms
    double tau = 1.509;
    double dm = 17.757;

    TH1D* h_t = new TH1D("h_t",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed = new TH1D("h_t_mixed",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed = new TH1D("h_t_unmixed",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged = new TH1D("h_t_untagegged",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_mp = new TH1D("h_t_mp",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0p = new TH1D("h_t_0p",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pp = new TH1D("h_t_pp",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_mm = new TH1D("h_t_mm",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0m = new TH1D("h_t_0m",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pm = new TH1D("h_t_pm",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_dt = new TH1D("h_dt",";#sigma_{t} (ps);Events (norm.) ",nBinst,0,0.15);
    TH1D* h_eta_OS = new TH1D("h_eta_OS",";#eta_{OS};Events (norm.) ",nBinst,0,0.5);
    TH1D* h_eta_SS = new TH1D("h_eta_SS",";#eta_{SS};Events (norm.) ",nBinst,0,0.5);

    TH1D* h_N_mixed = new TH1D("h_N_mixed",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed = (TH1D*) h_N_mixed->Clone("h_N_unmixed");

    TH1D* h_N_mixed_p = (TH1D*) h_N_mixed->Clone("h_N_mixed_p");
    TH1D* h_N_unmixed_p = (TH1D*) h_N_mixed->Clone("h_N_unmixed_p");
    TH1D* h_N_mixed_m = (TH1D*) h_N_mixed->Clone("h_N_mixed_m");
    TH1D* h_N_unmixed_m = (TH1D*) h_N_mixed->Clone("h_N_unmixed_m");

    TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,1,4);
    TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,2);
    TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30)  ;
    TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);

    TH1D* s_Kpipi_mixed_p = (TH1D*) s_Kpipi->Clone("s_Kpipi_mixed_p");
    TH1D* s_Kpipi_unmixed_p = (TH1D*) s_Kpipi->Clone("s_Kpipi_unmixed_p");
    TH1D* s_Kpipi_mixed_m = (TH1D*) s_Kpipi->Clone("s_Kpipi_mixed_m");
    TH1D* s_Kpipi_unmixed_m = (TH1D*) s_Kpipi->Clone("s_Kpipi_unmixed_m");

    /// Fit histograms
    TH1D* h_t_fit = new TH1D("h_t_fit",";t",nBinst,min_TAU,max_TAU);
    
    TH1D* h_t_mixed_fit = new TH1D("h_t_mixed_fit",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_fit = new TH1D("h_t_unmixed_fit",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged_fit = new TH1D("h_t_untagegged_fit",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_fit_mp = new TH1D("h_t_fit_mp",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0p = new TH1D("h_t_fit_0p",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pp = new TH1D("h_t_fit_pp",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_mm = new TH1D("h_t_fit_mm",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0m = new TH1D("h_t_fit_0m",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pm = new TH1D("h_t_fit_pm",";t (ps);Events (norm.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_dt_fit = new TH1D("h_dt_fit",";#sigma_{t} (ps);Events (norm.) ",nBinst,0,0.15);
    TH1D* h_eta_OS_fit = new TH1D("h_eta_OS_fit",";#eta_{OS};Events (norm.) ",nBinst,0,0.5);
    TH1D* h_eta_SS_fit = new TH1D("h_eta_SS_fit",";#eta_{SS};Events (norm.) ",nBinst,0,0.5);

    TH1D* h_N_mixed_fit = (TH1D*) h_N_mixed->Clone("h_N_mixed_fit");
    TH1D* h_N_unmixed_fit = (TH1D*) h_N_mixed->Clone("h_N_unmixed_fit");

    TH1D* h_N_mixed_p_fit = (TH1D*) h_N_mixed->Clone("h_N_mixed_p_fit");
    TH1D* h_N_unmixed_p_fit = (TH1D*) h_N_mixed->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit = (TH1D*) h_N_mixed->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit = (TH1D*) h_N_mixed->Clone("h_N_unmixed_m_fit");

    TH1D* s_Kpipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,1,4);
    TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,2);
    TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,2);
    TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsKpi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);

    TH1D* s_Kpipi_mixed_p_fit = (TH1D*) s_Kpipi->Clone("s_Kpipi_mixed_p_fit");
    TH1D* s_Kpipi_mixed_m_fit = (TH1D*) s_Kpipi->Clone("s_Kpipi_mixed_m_fit");
    TH1D* s_Kpipi_unmixed_p_fit = (TH1D*) s_Kpipi->Clone("s_Kpipi_unmixed_p_fit");
    TH1D* s_Kpipi_unmixed_m_fit = (TH1D*) s_Kpipi->Clone("s_Kpipi_unmixed_m_fit");

    TH1D* s_Kpipi_A = (TH1D*) s_Kpipi->Clone("s_Kpipi_A");
    TH1D* s_Kpipi_Abar = (TH1D*) s_Kpipi->Clone("s_Kpipi_Abar");
    TH1D* s_Kpi_A = (TH1D*) s_Kpi->Clone("s_Kpi_A");
    TH1D* s_Kpi_Abar = (TH1D*) s_Kpi->Clone("s_Kpi_Abar");
    TH1D* s_pipi_A = (TH1D*) s_pipi->Clone("s_pipi_A");
    TH1D* s_pipi_Abar = (TH1D*) s_pipi->Clone("s_pipi_Abar");
    TH1D* s_Dspipi_A = (TH1D*) s_Dspipi->Clone("s_Dspipi_A");
    TH1D* s_Dspipi_Abar = (TH1D*) s_Dspipi->Clone("s_Dspipi_Abar");
    TH1D* s_Dspi_A = (TH1D*) s_Dspi->Clone("s_Dspi_A");
    TH1D* s_Dspi_Abar = (TH1D*) s_Dspi->Clone("s_Dspi_Abar");
    TH1D* s_Dspim_A = (TH1D*) s_Dspim->Clone("s_Dspim_A");
    TH1D* s_Dspim_Abar = (TH1D*) s_Dspim->Clone("s_Dspim_Abar");
    TH1D* s_DsKpi_A = (TH1D*) s_DsKpi->Clone("s_DsKpi_A");
    TH1D* s_DsKpi_Abar = (TH1D*) s_DsKpi->Clone("s_DsKpi_Abar");
    TH1D* s_DsK_A = (TH1D*) s_DsK->Clone("s_DsK_A");
    TH1D* s_DsK_Abar = (TH1D*) s_DsK->Clone("s_DsK_Abar");

    TH1D* s_Kpipi_r = (TH1D*) s_Kpipi->Clone("s_Kpipi_r");

    /// Loop over data
    for (int i=0; i<eventList.size(); i++) {
        h_t->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        h_dt->Fill(eventList[i].getValueFromVector(1),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(3) != 0)h_eta_OS->Fill(eventList[i].getValueFromVector(4),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(5) != 0)h_eta_SS->Fill(eventList[i].getValueFromVector(6),eventList[i].getWeight());

        int f_evt = eventList[i].getValueFromVector(2);
        int q1 = eventList[i].getValueFromVector(3);
        int q2 = eventList[i].getValueFromVector(5);   
        int q_eff = 0;
        double w_eff = 0.5;

        if(q1 != 0 && q2 != 0){
            std::pair<double, double> calibrated_mistag_os = pdf.getCalibratedMistag_OS(eventList[i]);
            std::pair<double, double> calibrated_mistag_ss = pdf.getCalibratedMistag_SS(eventList[i]);
            
            double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
            double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
            
            if( p/(p+p_bar) > 0.5 ){ 
		q_eff = 1;
		w_eff = 1-p/(p+p_bar);
	    }
            else if( p/(p+p_bar) < 0.5 ){
		 q_eff = -1;
		 w_eff = p/(p+p_bar);
            }
            
        }
        else if( q1 != 0){
            q_eff = q1;
        }
        else if( q2 != 0){
            q_eff = q2;
        } 
        
        if((string)channel=="signal"){

            if(q_eff==-1 && f_evt == 1){ 
			h_t_mp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			if(w_eff<w_max)h_N_mixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
            }
	    else if(q_eff==0 && f_evt == 1)h_t_0p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == 1){
			h_t_pp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			if(w_eff<w_max)h_N_unmixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
	    }
	    else if(q_eff==-1 && f_evt == -1){
			h_t_mm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        	    	if(w_eff<w_max)h_N_unmixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
	    }
	    else if(q_eff==0 && f_evt == -1)h_t_0m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == -1){
			h_t_pm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			if(w_eff<w_max)h_N_mixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
	    }
        }
      
        else {
            if(q_eff == 0)h_t_untagegged->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff*f_evt > 0  ){
		h_t_mixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
		if(w_eff<w_max)h_N_mixed->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
	    }
            else {
		h_t_unmixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
	    	if(w_eff<w_max)h_N_unmixed->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
	    }
        }

            s_Kpipi->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
            s_Kpi->Fill(eventList[i].s(2,4)/(GeV*GeV),eventList[i].getWeight());
            s_pipi->Fill(eventList[i].s(3,4)/(GeV*GeV),eventList[i].getWeight());
            s_Dspipi->Fill(eventList[i].sij(s134)/(GeV*GeV),eventList[i].getWeight());
            s_DsK->Fill(eventList[i].s(1,2)/(GeV*GeV),eventList[i].getWeight());
            s_DsKpi->Fill(eventList[i].sij(s124)/(GeV*GeV),eventList[i].getWeight());
            s_Dspi->Fill(eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());
            s_Dspim->Fill(eventList[i].s(1,4)/(GeV*GeV),eventList[i].getWeight());

	    if(w_eff<w_max){
		//if(abs(fmod(eventList[i].getValueFromVector(0),2.*pi/dm)-0.18) < 2.*pi/dm/4.*0.5) s_Kpipi_mixed->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
		//else if(abs(fmod(eventList[i].getValueFromVector(0),2.*pi/dm)-0.18) > 2.*pi/dm/4. *1.5)s_Kpipi_unmixed->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
		if(cos(dm*eventList[i].getValueFromVector(0))>0){
			if(q_eff*f_evt < 0 )s_Kpipi_mixed_p->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
			else s_Kpipi_unmixed_p->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
	    	}
		else {
			if(q_eff*f_evt < 0 )s_Kpipi_mixed_m->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
			else s_Kpipi_unmixed_m->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
		}
	    }
    }     

    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");

    RooRealVar* r_t = new RooRealVar("t", "t",min_TAU,max_TAU);
    RooRealVar* r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.1);
    RooRealVar* r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5); 
    RooRealVar* r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5); 

    RooExponential gen_t("gen_t","gen_t", *r_t, RooRealConstant::value(1./tau));	
    RooGaussian gen_dt("gen_dt","gen_dt", *r_dt, RooRealConstant::value(0.03),RooRealConstant::value(0.015));	
    RooGaussian gen_eta_OS("gen_eta_OS","gen_eta_OS", *r_eta_OS, RooRealConstant::value(0.4),RooRealConstant::value(0.2));	
    RooGaussian gen_eta_SS("gen_eta_SS","gen_eta_SS", *r_eta_SS, RooRealConstant::value(0.5),RooRealConstant::value(0.1));

    /// Loop over MC
    for(int i = 0; i < eventListMC.size(); i++){
        
        double t_MC = 0.;
    	while(1) {
 		//double rand = ranLux.Uniform();
 		double tval = ranLux.Exp(tau); //- tau*log(rand);
  	    	if (tval< max_TAU && tval> min_TAU) {
         		t_MC = tval ;
			r_t->setVal(tval);
         		break ;
		}
       }
	
	//gen_t.generateEvent(1);
	//double t_MC = r_t->getVal();

	gen_dt.generateEvent(1);
	double dt_MC = r_dt->getVal(); //ranLux.Uniform(0., 0.1);
        
        //if(t_MC < min_TAU/2. || t_MC > max_TAU*1.2 )continue;
        //if( dt_MC < 0 || dt_MC > 0.1 )continue;
        
        double q_rand = ranLux.Uniform();
       
        int q_OS_MC = 0;
        if (q_rand < 1./3.) q_OS_MC = -1;
        if (q_rand > 2./3.) q_OS_MC = 1;
        
        q_rand = ranLux.Uniform();
        int q_SS_MC = 0;
        if (q_rand < 1./3.) q_SS_MC = -1;
        if (q_rand > 2./3.) q_SS_MC = 1;
        
	//const double eta_OS_MC = ranLux.Uniform(0., 0.5);
        //const double eta_SS_MC = ranLux.Uniform(0., 0.5);
	gen_eta_OS.generateEvent(1);
	double eta_OS_MC = r_eta_OS->getVal();        

	gen_eta_SS.generateEvent(1);
	double eta_SS_MC = r_eta_SS->getVal(); 

        q_rand = ranLux.Uniform();
        int f_MC = 0;
        if (q_rand > .5) f_MC = -1;
        else f_MC = 1;
     
	DalitzEvent evt(eventListMC.getEvent(i));
 	//if(!(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2))continue;

        evt.setValueInVector(0, t_MC);
        evt.setValueInVector(1, dt_MC);
        evt.setValueInVector(2, f_MC);
        evt.setValueInVector(3, q_OS_MC);
        evt.setValueInVector(4, eta_OS_MC);
        evt.setValueInVector(5, q_SS_MC);
        evt.setValueInVector(6, eta_SS_MC);
        
        //const double pdfVal = pdf.getVal(evt);
        const double pdfVal = pdf.un_normalised_noPs(evt);

        //double weight = pdfVal/exp(-fabs(t_MC)/(tau))*tau*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
        double weight = pdfVal*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	weight /=  exp(-t_MC/tau) / ( tau * ( exp(min_TAU/tau) - exp(max_TAU/tau) ) ) * gen_dt.getVal() * gen_eta_OS.getVal()* gen_eta_SS.getVal();

        h_t_fit->Fill(t_MC,weight);
        h_dt_fit->Fill(dt_MC,weight);
        if(evt.getValueFromVector(3) != 0)h_eta_OS_fit->Fill(evt.getValueFromVector(4),weight);
        if(evt.getValueFromVector(5) != 0)h_eta_SS_fit->Fill(evt.getValueFromVector(6),weight);
        
        int f_evt = evt.getValueFromVector(2);
        int q1 = evt.getValueFromVector(3);
        int q2 = evt.getValueFromVector(5);   
        int q_eff = 0;
	double w_eff = 0.5;
        
        if(q1 != 0 && q2 != 0){
            std::pair<double, double> calibrated_mistag_os = pdf.getCalibratedMistag_OS(evt);
            std::pair<double, double> calibrated_mistag_ss = pdf.getCalibratedMistag_SS(evt);
            
            double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
            double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
            
            if( p/(p+p_bar) > 0.5 ){ 
		q_eff = 1;
		w_eff = 1-p/(p+p_bar);
	    }
            else if( p/(p+p_bar) < 0.5 ){
		 q_eff = -1;
		 w_eff = p/(p+p_bar);
            }
            
        }
        else if( q1 != 0){
            q_eff = q1;
        }
        else if( q2 != 0){
            q_eff = q2;
        } 
        
        if((string)channel=="signal"){
            
            if(q_eff==-1 && f_evt == 1){
			h_t_fit_mp->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_mixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
	    }
            else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(evt.getValueFromVector(0),weight);
            else if(q_eff==1 && f_evt == 1){
			h_t_fit_pp->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_unmixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
            }
	    else if(q_eff==-1 && f_evt == -1){
			h_t_fit_mm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_unmixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
            }
	    else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(evt.getValueFromVector(0),weight);
            else if(q_eff==1 && f_evt == -1){
			h_t_fit_pm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_mixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
            }
        }
        
        else {
            if(q_eff == 0)h_t_untagegged_fit->Fill(evt.getValueFromVector(0),weight);
            else if(q_eff*f_evt > 0  ){
		h_t_mixed_fit->Fill(evt.getValueFromVector(0),weight);
		if(w_eff<w_max)h_N_mixed_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
	    }
            else{ 
		h_t_unmixed_fit->Fill(evt.getValueFromVector(0),weight);
		if(w_eff<w_max)h_N_unmixed_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
	    }
        }

        s_Kpipi_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
        s_Kpi_fit->Fill(evt.s(2,4)/(GeV*GeV),weight);
        s_pipi_fit->Fill(evt.s(3,4)/(GeV*GeV),weight);
        s_Dspipi_fit->Fill(evt.sij(s134)/(GeV*GeV),weight);
        s_DsK_fit->Fill(evt.s(1,2)/(GeV*GeV),weight);
        s_DsKpi_fit->Fill(evt.sij(s124)/(GeV*GeV),weight);
        s_Dspi_fit->Fill(evt.s(1,3)/(GeV*GeV),weight);
        s_Dspim_fit->Fill(evt.s(1,4)/(GeV*GeV),weight);

	double weight_A = fas.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	//weight_A /=  exp(-t_MC/tau) / ( tau * ( exp(min_TAU/tau) - exp(max_TAU/tau) ) ) * gen_dt.getVal() * gen_eta_OS.getVal()* gen_eta_SS.getVal();

	double weight_Abar = fasCC.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	//weight_Abar /=  exp(-t_MC/tau) / ( tau * ( exp(min_TAU/tau) - exp(max_TAU/tau) ) ) * gen_dt.getVal() * gen_eta_OS.getVal()* gen_eta_SS.getVal();

	s_Kpipi_A->Fill(evt.sij(s234)/(GeV*GeV),weight_A);
        s_Kpipi_Abar->Fill(evt.sij(s234)/(GeV*GeV),weight_Abar);
	s_Kpi_A->Fill(evt.s(2,4)/(GeV*GeV),weight_A);
        s_Kpi_Abar->Fill(evt.s(2,4)/(GeV*GeV),weight_Abar);
	s_pipi_A->Fill(evt.s(3,4)/(GeV*GeV),weight_A);
        s_pipi_Abar->Fill(evt.s(3,4)/(GeV*GeV),weight_Abar);
	s_Dspipi_A->Fill(evt.sij(s134)/(GeV*GeV),weight_A);
        s_Dspipi_Abar->Fill(evt.sij(s134)/(GeV*GeV),weight_Abar);
	s_Dspi_A->Fill(evt.s(1,3)/(GeV*GeV),weight_A);
        s_Dspi_Abar->Fill(evt.s(1,3)/(GeV*GeV),weight_Abar);
	s_Dspim_A->Fill(evt.s(1,4)/(GeV*GeV),weight_A);
        s_Dspim_Abar->Fill(evt.s(1,4)/(GeV*GeV),weight_Abar);
	s_DsK_A->Fill(evt.s(1,4)/(GeV*GeV),weight_A);
        s_DsK_Abar->Fill(evt.s(1,4)/(GeV*GeV),weight_Abar);
	s_DsKpi_A->Fill(evt.sij(s124)/(GeV*GeV),weight_A);
        s_DsKpi_Abar->Fill(evt.sij(s124)/(GeV*GeV),weight_Abar);

	if(w_eff<w_max){
		//if(abs(fmod(evt.getValueFromVector(0),2.*pi/dm)-0.18) < 2.*pi/dm/4.*0.5) s_Kpipi_mixed_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
		//else if(abs(fmod(evt.getValueFromVector(0),2.*pi/dm)-0.18) > 2.*pi/dm/4. *1.5)s_Kpipi_unmixed_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
		if(cos(dm*evt.getValueFromVector(0))>0){
			if(q_eff*f_evt < 0 )s_Kpipi_mixed_p_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
			else s_Kpipi_unmixed_p_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
		}
		else {
			if(q_eff*f_evt < 0 )s_Kpipi_mixed_m_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
			else s_Kpipi_unmixed_m_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
		}
	}

    }

    /// Plot
    TCanvas* c = new TCanvas();

    h_t->SetMinimum(0.1);    
    h_t->SetLineColor(kBlack);
    h_t->DrawNormalized("e",1);
        
    h_t_fit->SetLineColor(kBlue);
    h_t_fit->SetLineWidth(3);
    h_t_fit->SetMarkerColor(kBlue); 
    h_t_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_t.eps").c_str());
    gPad->SetLogy(1);
    c->Print(((string)OutputDir+"h_t_log.eps").c_str());
    gPad->SetLogy(0);
    
    h_dt->SetMinimum(0);        
    h_dt->SetLineColor(kBlack);
    h_dt->DrawNormalized("e1",1);
    h_dt_fit->SetLineColor(kBlue);
    h_dt_fit->SetLineWidth(3);
    h_dt_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_dt.eps").c_str());
    
    h_eta_OS->SetMinimum(0);        
    h_eta_OS->SetLineColor(kBlack);
    h_eta_OS->DrawNormalized("e1",1);
    h_eta_OS_fit->SetLineColor(kBlue);
    h_eta_OS_fit->SetLineWidth(3);
    h_eta_OS_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_eta_OS.eps").c_str());
    
    h_eta_SS->SetMinimum(0);        
    h_eta_SS->SetLineColor(kBlack);
    h_eta_SS->DrawNormalized("e1",1);
    h_eta_SS_fit->SetLineColor(kBlue);
    h_eta_SS_fit->SetLineWidth(3);
    h_eta_SS_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_eta_SS.eps").c_str());


    if((string)channel=="norm"){

        h_t_mixed->SetMarkerColor(kRed); 
        h_t_mixed->SetLineColor(kRed);
        h_t_mixed->DrawNormalized("e",1);
        
        h_t_unmixed->SetMarkerColor(kBlue); 
        h_t_unmixed->SetLineColor(kBlue);
        h_t_unmixed->DrawNormalized("esame",1);
        
        h_t_mixed_fit->SetMarkerColor(kRed); 
        h_t_mixed_fit->SetLineColor(kRed);
        h_t_mixed_fit->DrawNormalized("histcsame",1);
        
        h_t_unmixed_fit->SetMarkerColor(kBlue); 
        h_t_unmixed_fit->SetLineColor(kBlue);
        h_t_unmixed_fit->DrawNormalized("histcsame",1);
        
        h_t_untagegged->SetMarkerColor(kGreen); 
        h_t_untagegged->SetLineColor(kGreen);
        //h_t_untagegged->DrawNormalized("esame",1);
        
        h_t_untagegged_fit->SetMarkerColor(kGreen); 
        h_t_untagegged_fit->SetLineColor(kGreen);
        //h_t_untagegged_fit->DrawNormalized("histcsame",1);
        
        c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());

	TH1D* h_asym = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);	
        h_asym->SetMinimum(-0.25);
	h_asym->SetMaximum(0.25);
	TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
	h_asym_fit->SetLineColor(kRed);
	h_asym->Draw("e");
	h_asym_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());

    }

    else{

        h_t_mp->Scale(1./h_t_mp->Integral());
        double maxY= h_t_mp->GetMaximum()*1.3;        
        h_t_mp->SetMinimum(0.);  
        h_t_mp->SetMaximum(maxY);        
        h_t_mp->SetMarkerColor(kBlue); 
        h_t_mp->SetLineColor(kBlue);
        h_t_mp->DrawNormalized("e",1);
    
        h_t_pp->SetMarkerColor(kRed); 
        h_t_pp->SetLineColor(kRed);
        h_t_pp->DrawNormalized("esame",1);
        
        h_t_fit_mp->SetMarkerColor(kBlue); 
        h_t_fit_mp->SetLineColor(kBlue);
        h_t_fit_mp->DrawNormalized("histcsame",1);
        
        h_t_fit_pp->SetMarkerColor(kRed); 
        h_t_fit_pp->SetLineColor(kRed);
        h_t_fit_pp->DrawNormalized("histcsame",1);
        
        c->Print(((string)OutputDir+"h_t_mixed_p.eps").c_str());

        h_t_pm->Scale(1./h_t_pm->Integral());
        h_t_pm->SetMinimum(0.);     
        h_t_pm->SetMaximum(maxY);        
        h_t_pm->SetMarkerColor(kBlue); 
        h_t_pm->SetLineColor(kBlue);
        h_t_pm->DrawNormalized("e",1);
        
        h_t_mm->SetMarkerColor(kRed); 
        h_t_mm->SetLineColor(kRed);
        h_t_mm->DrawNormalized("esame",1);
        
        h_t_fit_pm->SetMarkerColor(kBlue); 
        h_t_fit_pm->SetLineColor(kBlue);
        h_t_fit_pm->DrawNormalized("histcsame",1);
        
        h_t_fit_mm->SetMarkerColor(kRed); 
        h_t_fit_mm->SetLineColor(kRed);
        h_t_fit_mm->DrawNormalized("histcsame",1);
        
        c->Print(((string)OutputDir+"h_t_mixed_m.eps").c_str());

	TH1D* h_asym_p = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_mixed_p);	
        //h_asym_p->SetMinimum(-20);
	//h_asym_p->SetMaximum(20);
	TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);	
	h_asym_p_fit->SetLineColor(kRed);
	h_asym_p->Draw("e");
	h_asym_p_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_p.eps").c_str());

	TH1D* h_asym_m = (TH1D*) h_N_unmixed_m->GetAsymmetry(h_N_mixed_m);	
        //h_asym_m->SetMinimum(-20);
	//h_asym_m->SetMaximum(20);
	TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);	
	h_asym_m_fit->SetLineColor(kRed);
	h_asym_m->Draw("e");
	h_asym_m_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_m.eps").c_str());

	h_asym_p_fit->Draw("histc");
	h_asym_m_fit->SetLineColor(kBlue);
 	h_asym_m_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());	
    }
    
            s_Kpipi->SetMinimum(0);
            s_Kpipi->SetLineColor(kBlack);
            s_Kpipi->DrawNormalized("e1",1);
            s_Kpipi_fit->SetLineColor(kBlue);
            s_Kpipi_fit->SetLineWidth(3);
            s_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi.eps").c_str());
            
            s_Kpi->SetMinimum(0);
            s_Kpi->SetLineColor(kBlack);
            s_Kpi->DrawNormalized("e1",1);
            s_Kpi_fit->SetLineColor(kBlue);
            s_Kpi_fit->SetLineWidth(3);
            s_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi.eps").c_str());
            
	    s_pipi->SetMinimum(0);            
            s_pipi->SetLineColor(kBlack);
            s_pipi->DrawNormalized("e1",1);
            s_pipi_fit->SetLineColor(kBlue);
            s_pipi_fit->SetLineWidth(3);
            s_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi.eps").c_str());
            
	    s_Dspipi->SetMinimum(0);            
            s_Dspipi->SetLineColor(kBlack);
            s_Dspipi->DrawNormalized("e1",1);
            s_Dspipi_fit->SetLineColor(kBlue);
            s_Dspipi_fit->SetLineWidth(3);
            s_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi.eps").c_str());
           
	    s_DsK->SetMinimum(0);
            s_DsK->SetLineColor(kBlack);
            s_DsK->DrawNormalized("e1",1);
            s_DsK_fit->SetLineColor(kBlue);
            s_DsK_fit->SetLineWidth(3);
            s_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK.eps").c_str());

	    s_DsKpi->SetMinimum(0);            
            s_DsKpi->SetLineColor(kBlack);
            s_DsKpi->DrawNormalized("e1",1);
            s_DsKpi_fit->SetLineColor(kBlue);
            s_DsKpi_fit->SetLineWidth(3);
            s_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi.eps").c_str());

	    s_Dspi->SetMinimum(0);
            s_Dspi->SetLineColor(kBlack);
            s_Dspi->DrawNormalized("e1",1);
            s_Dspi_fit->SetLineColor(kBlue);
            s_Dspi_fit->SetLineWidth(3);
            s_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi.eps").c_str());

	    s_Dspim->SetMinimum(0);
            s_Dspim->SetLineColor(kBlack);
            s_Dspim->DrawNormalized("e1",1);
            s_Dspim_fit->SetLineColor(kBlue);
            s_Dspim_fit->SetLineWidth(3);
            s_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim.eps").c_str());

/*
	    s_Kpipi_unmixed->Rebin(4);
	    s_Kpipi_mixed->Rebin(4);
            s_Kpipi_unmixed_fit->Rebin(4);
            s_Kpipi_mixed_fit->Rebin(4);
            s_Kpipi_unmixed->SetMinimum(0);
            s_Kpipi_unmixed->SetLineColor(kRed);
            s_Kpipi_unmixed->SetMarkerColor(kRed);
            s_Kpipi_unmixed->DrawNormalized("e1",1);
            s_Kpipi_unmixed_fit->SetLineColor(kRed);
            s_Kpipi_unmixed_fit->SetLineWidth(3);
            s_Kpipi_unmixed_fit->DrawNormalized("histcsame",1);
            s_Kpipi_mixed->SetMarkerColor(kBlue);
	    s_Kpipi_mixed->SetLineColor(kBlue);
            s_Kpipi_mixed->DrawNormalized("esame",1);
            s_Kpipi_mixed_fit->SetLineColor(kBlue);
            s_Kpipi_mixed_fit->SetLineWidth(3);
            s_Kpipi_mixed_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi_asym2.eps").c_str());
*/

	    //s_Kpipi_unmixed_p->Add(s_Kpipi_unmixed_m,-1.);	  
	    //s_Kpipi_mixed_p->Add(s_Kpipi_mixed_m,-1.);	  
	    //s_Kpipi_unmixed_p_fit->Add(s_Kpipi_unmixed_m_fit,-1.);	  
	    //s_Kpipi_mixed_p_fit->Add(s_Kpipi_mixed_m_fit,-1.);	  
/*
  	    for (int i=1; i<= s_Kpipi_unmixed_m->GetNbinsX(); i++) {
      		  s_Kpipi_mixed_m->SetBinContent(i,-s_Kpipi_mixed_m->GetBinContent(i));
      		  s_Kpipi_unmixed_m->SetBinContent(i,-s_Kpipi_unmixed_m->GetBinContent(i));
      		  s_Kpipi_mixed_m_fit->SetBinContent(i,-s_Kpipi_mixed_m_fit->GetBinContent(i));
      		  s_Kpipi_unmixed_m_fit->SetBinContent(i,-s_Kpipi_unmixed_m_fit->GetBinContent(i));
  	    }
*/

    	    s_Kpipi_unmixed_p->Rebin(4);
	    s_Kpipi_mixed_p->Rebin(4);
            s_Kpipi_unmixed_p_fit->Rebin(4);
            s_Kpipi_mixed_p_fit->Rebin(4);
        
    	    s_Kpipi_unmixed_m->Rebin(4);
	    s_Kpipi_mixed_m->Rebin(4);
            s_Kpipi_unmixed_m_fit->Rebin(4);
            s_Kpipi_mixed_m_fit->Rebin(4);
        
	    TH1D* s_Kpipi_asym_p = (TH1D*) s_Kpipi_unmixed_p->GetAsymmetry(s_Kpipi_mixed_p);	
	    //s_Kpipi_asym_p->Rebin(4);
	    TH1D* s_Kpipi_asym_m = (TH1D*) s_Kpipi_unmixed_m->GetAsymmetry(s_Kpipi_mixed_m);	
	    //s_Kpipi_asym_m->Rebin(4);
  	    for (int i=1; i<= s_Kpipi_asym_m->GetNbinsX(); i++) {
      		  s_Kpipi_asym_m->SetBinContent(i,-s_Kpipi_asym_m->GetBinContent(i));
	     }

	    s_Kpipi_asym_p->SetBit(TH1::kIsAverage);
      	    s_Kpipi_asym_m->SetBit(TH1::kIsAverage);
	    s_Kpipi_asym_p->Add(s_Kpipi_asym_m);

	    TH1D* s_Kpipi_asym_p_fit = (TH1D*) s_Kpipi_unmixed_p_fit->GetAsymmetry(s_Kpipi_mixed_p_fit);	
	    //s_Kpipi_asym_p_fit->Rebin(4);
	    TH1D* s_Kpipi_asym_m_fit = (TH1D*) s_Kpipi_unmixed_m_fit->GetAsymmetry(s_Kpipi_mixed_m_fit);	
	    //s_Kpipi_asym_m_fit->Rebin(4);
  	    for (int i=1; i<= s_Kpipi_asym_m_fit->GetNbinsX(); i++) {
      		  s_Kpipi_asym_m_fit->SetBinContent(i,-s_Kpipi_asym_m_fit->GetBinContent(i));
	     }

	    s_Kpipi_asym_p_fit->SetBit(TH1::kIsAverage);
      	    s_Kpipi_asym_m_fit->SetBit(TH1::kIsAverage);
	    s_Kpipi_asym_p_fit->Add(s_Kpipi_asym_m_fit);
	    s_Kpipi_asym_p_fit->SetLineColor(kRed);

            s_Kpipi_asym_p->SetMinimum(-.5);
	    s_Kpipi_asym_p->SetMaximum(.5);
	    s_Kpipi_asym_p->Draw("e");
	    s_Kpipi_asym_p_fit->Draw("histcsame");
	    
	    //s_Kpipi_asym_m_fit->SetLineColor(kBlue);
	    //s_Kpipi_asym_m->SetLineColor(kBlue);
	    //s_Kpipi_asym_m_fit->Draw("histcsame");
	    //s_Kpipi_asym_m->Draw("esame");


	    //s_Kpipi_asym_p_fit->Draw("histcsame");
            c->Print(((string)OutputDir+"s_Kpipi_asym.eps").c_str());
            c->Print(((string)OutputDir+"s_Kpipi_asym.C").c_str());
	    cout << "s_Kpipi_asym " << s_Kpipi_asym_p->Integral() << endl; 
	    cout << "s_Kpipi_asym_fit " << s_Kpipi_asym_p_fit->Integral() << endl; 

            s_Kpipi_A->SetMinimum(0);
            s_Kpipi_A->SetLineColor(kBlue);
            s_Kpipi_A->SetLineWidth(3);
	    s_Kpipi_A->DrawNormalized("histc",1);
            s_Kpipi_Abar->SetLineColor(kRed);
            s_Kpipi_Abar->SetLineWidth(3);
            s_Kpipi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi_AAbar.eps").c_str());

            s_Kpi_A->SetMinimum(0);
            s_Kpi_A->SetLineColor(kBlue);
            s_Kpi_A->SetLineWidth(3);
	    s_Kpi_A->DrawNormalized("histc",1);
            s_Kpi_Abar->SetLineColor(kRed);
            s_Kpi_Abar->SetLineWidth(3);
            s_Kpi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi_AAbar.eps").c_str());

            s_pipi_A->SetMinimum(0);
            s_pipi_A->SetLineColor(kBlue);
            s_pipi_A->SetLineWidth(3);
	    s_pipi_A->DrawNormalized("histc",1);
            s_pipi_Abar->SetLineColor(kRed);
            s_pipi_Abar->SetLineWidth(3);
            s_pipi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi_AAbar.eps").c_str());

            s_Dspipi_A->SetMinimum(0);
            s_Dspipi_A->SetLineColor(kBlue);
            s_Dspipi_A->SetLineWidth(3);
	    s_Dspipi_A->DrawNormalized("histc",1);
            s_Dspipi_Abar->SetLineColor(kRed);
            s_Dspipi_Abar->SetLineWidth(3);
            s_Dspipi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi_AAbar.eps").c_str());

            s_Dspi_A->SetMinimum(0);
            s_Dspi_A->SetLineColor(kBlue);
            s_Dspi_A->SetLineWidth(3);
	    s_Dspi_A->DrawNormalized("histc",1);
            s_Dspi_Abar->SetLineColor(kRed);
            s_Dspi_Abar->SetLineWidth(3);
            s_Dspi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi_AAbar.eps").c_str());

            s_Dspim_A->SetMinimum(0);
            s_Dspim_A->SetLineColor(kBlue);
            s_Dspim_A->SetLineWidth(3);
	    s_Dspim_A->DrawNormalized("histc",1);
            s_Dspim_Abar->SetLineColor(kRed);
            s_Dspim_Abar->SetLineWidth(3);
            s_Dspim_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim_AAbar.eps").c_str());

            s_DsK_A->SetMinimum(0);
            s_DsK_A->SetLineColor(kBlue);
            s_DsK_A->SetLineWidth(3);
	    s_DsK_A->DrawNormalized("histc",1);
            s_DsK_Abar->SetLineColor(kRed);
            s_DsK_Abar->SetLineWidth(3);
            s_DsK_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK_AAbar.eps").c_str());

            s_DsKpi_A->SetMinimum(0);
            s_DsKpi_A->SetLineColor(kBlue);
            s_DsKpi_A->SetLineWidth(3);
	    s_DsKpi_A->DrawNormalized("histc",1);
            s_DsKpi_Abar->SetLineColor(kRed);
            s_DsKpi_Abar->SetLineWidth(3);
            s_DsKpi_Abar->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi_AAbar.eps").c_str());


	    //getChi2(eventList,eventListMC);




    if(do2DScan == 1){
        cout << "Now doing 2D scan:" << endl;
        
        Neg2LL fcn(pdf, eventList_f);    
        Neg2LL fcn_bar(pdf, eventList_f_bar);    
        
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
                
                double v = neg2LL.getNewVal();  
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
                double v = neg2LL.getNewVal();
                
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

    
    return;
}

void produceMarginalPdfs(){

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);
    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);

    TString prefix = "";
    //TString prefix = "BsTaggingTool_";

    /// Load files
    
    // Data
    int q_OS;
    Short_t q_SS;
    double w_OS;
    Float_t w_SS;
    double sw;
    int year,Ds_finalState;
    double dt;
    double Bs_pt,Bs_eta,nTracks;
    
    TChain* tree=new TChain("DecayTree");
    tree->Add( ((string)InputDir + "Data/signal.root").c_str());
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("N_Bs_sw",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("*DEC",1);
    tree->SetBranchStatus("*PROB",1);
    tree->SetBranchStatus("*OS",1);
    tree->SetBranchStatus("*TAU*",1);
    tree->SetBranchStatus("*ETA",1);
    tree->SetBranchStatus("*PT",1);
    tree->SetBranchStatus("NTracks",1);

    tree->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
    tree->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&w_OS);
    tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
    tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&w_SS);
    tree->SetBranchAddress("N_Bs_sw",&sw);
    tree->SetBranchAddress("year",&year);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
    tree->SetBranchAddress("Bs_PT",&Bs_pt);
    tree->SetBranchAddress("Bs_ETA",&Bs_eta);
    tree->SetBranchAddress("NTracks",&nTracks);

    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm.root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*DEC",1);
    tree_norm->SetBranchStatus("*PROB",1);
    tree_norm->SetBranchStatus("*OS",1);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("*ETA",1);
    tree_norm->SetBranchStatus("*PT",1);
    tree_norm->SetBranchStatus("NTracks",1);
    
    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&w_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&w_SS);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree_norm->SetBranchAddress("Bs_DTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("Bs_PT",&Bs_pt);
    tree_norm->SetBranchAddress("Bs_ETA",&Bs_eta);
    tree_norm->SetBranchAddress("NTracks",&nTracks);

    // MC
    int q_OS_MC;
    Short_t q_SS_MC;
    double w_OS_MC;
    Float_t w_SS_MC;
    double w;
    int cat,yearMC,Ds_finalStateMC;
    double dt_MC;
    double Bs_pt_MC,Bs_eta_MC,nTracks_MC;

    TChain* treeMC =new TChain("DecayTree");
    treeMC->Add( ((string)InputDir + "MC/signal.root").c_str());
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("Bs_BKGCAT",1);
    treeMC->SetBranchStatus("weight",1);
    treeMC->SetBranchStatus("*DEC",1);
    treeMC->SetBranchStatus("*PROB",1);
    treeMC->SetBranchStatus("*OS",1);
    treeMC->SetBranchStatus("*TAU*",1);
    treeMC->SetBranchStatus("*ETA",1);
    treeMC->SetBranchStatus("*PT",1);
    treeMC->SetBranchStatus("NTracks",1);

    treeMC->SetBranchAddress("Bs_BKGCAT",&cat);
    treeMC->SetBranchAddress("year",&yearMC);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalStateMC);           
    treeMC->SetBranchAddress("weight",&w);
    treeMC->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS_MC);
    treeMC->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&w_OS_MC);
    treeMC->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS_MC);
    treeMC->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&w_SS_MC);
    treeMC->SetBranchAddress("Bs_DTF_TAUERR",&dt_MC);
    treeMC->SetBranchAddress("Bs_PT",&Bs_pt_MC);
    treeMC->SetBranchAddress("Bs_ETA",&Bs_eta_MC);
    treeMC->SetBranchAddress("NTracks",&nTracks_MC);
    
    TChain* treeMC_norm =new TChain("DecayTree");
    treeMC_norm->Add( ((string)InputDir + "MC/norm.root").c_str());
    treeMC_norm->SetBranchStatus("*",0);
    treeMC_norm->SetBranchStatus("Ds_finalState",1);
    treeMC_norm->SetBranchStatus("Bs_BKGCAT",1);
    treeMC_norm->SetBranchStatus("weight",1);
    treeMC_norm->SetBranchStatus("*DEC",1);
    treeMC_norm->SetBranchStatus("*PROB",1);
    treeMC_norm->SetBranchStatus("*OS",1);
    treeMC_norm->SetBranchStatus("*TAU*",1);
    treeMC_norm->SetBranchStatus("*ETA",1);
    treeMC_norm->SetBranchStatus("*PT",1);
    treeMC_norm->SetBranchStatus("NTracks",1);

    treeMC_norm->SetBranchAddress("Bs_BKGCAT",&cat);
    treeMC_norm->SetBranchAddress("year",&yearMC);
    treeMC_norm->SetBranchAddress("Ds_finalState",&Ds_finalStateMC);           
    treeMC_norm->SetBranchAddress("weight",&w);
    treeMC_norm->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS_MC);
    treeMC_norm->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&w_OS_MC);
    treeMC_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS_MC);
    treeMC_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&w_SS_MC);
    treeMC_norm->SetBranchAddress("Bs_DTF_TAUERR",&dt_MC);
    treeMC_norm->SetBranchAddress("Bs_PT",&Bs_pt_MC);
    treeMC_norm->SetBranchAddress("Bs_ETA",&Bs_eta_MC);
    treeMC_norm->SetBranchAddress("NTracks",&nTracks_MC);
    
    ///Make histograms
    int bins = 60;
    TH1D* h_w_OS = new TH1D("h_w_OS","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_MC = new TH1D("h_w_OS_MC","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_SS = new TH1D("h_w_SS","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_MC = new TH1D("h_w_SS_MC","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS = new TH1D("h_q_OS","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_MC = new TH1D("h_q_OS_MC","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_SS = new TH1D("h_q_SS","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_MC = new TH1D("h_q_SS_MC","; q_{SS}",3,-1.5,1.5);
    
    TH1D* h_dt = new TH1D("h_dt",";#sigma_{t} (ps);Events (norm.) ",bins,0,0.25);
    TH1D* h_dt_MC = new TH1D("h_dt_MC",";#sigma_{t} (ps);Events (norm.) ",bins,0,0.25);

    TH1D* h_w_OS_norm = new TH1D("h_w_OS_norm","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_MC_norm = new TH1D("h_w_OS_MC_norm","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_SS_norm = new TH1D("h_w_SS_norm","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_MC_norm = new TH1D("h_w_SS_MC_norm","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm = new TH1D("h_q_OS_norm","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_MC_norm = new TH1D("h_q_OS_MC_norm","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm = new TH1D("h_q_SS_norm","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_MC_norm = new TH1D("h_q_SS_MC_norm","; q_{SS}",3,-1.5,1.5);
    
    TH1D* h_dt_norm = new TH1D("h_dt_norm",";#sigma_{t} (ps);Events (norm.) ",bins,0,0.25);
    TH1D* h_dt_MC_norm = new TH1D("h_dt_MC_norm",";#sigma_{t} (ps);Events (norm.) ",bins,0,0.25);

    double eff_OS = 0; 
    double eff_SS = 0;
    double eff_OS_MC = 0; 
    double eff_SS_MC = 0;
    
    double eff_OS_norm = 0; 
    double eff_SS_norm = 0;
    double eff_OS_MC_norm = 0; 
    double eff_SS_MC_norm = 0;
    
    double sumw = 0;

    ///loop over data events
    for(int i=0; i< tree->GetEntries(); i++)
    {	
        //if (0ul == (i % 1000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
        tree->GetEntry(i);        
        //if(year > 13) continue;
        
        h_dt->Fill(dt,sw);
        h_q_OS->Fill((double)q_OS,sw);
        h_q_SS->Fill((double)q_SS,sw);

        if(q_OS != 0) {
            h_w_OS->Fill(w_OS,sw);
            eff_OS += sw;
        }
        if(q_SS != 0){
            h_w_SS->Fill(w_SS,sw);
            eff_SS += sw;
        }
        sumw += sw;
    }
    
    eff_OS /= sumw;
    eff_SS /= sumw;
    
    sumw = 0;
    for(int i=0; i< tree_norm->GetEntries(); i++)
    {	
        tree_norm->GetEntry(i);
        //if(year > 13) continue;
        
        h_dt_norm->Fill(dt,sw);
        h_q_OS_norm->Fill((double)q_OS,sw);
        h_q_SS_norm->Fill((double)q_SS,sw);
        
        if(q_OS != 0){
            h_w_OS_norm->Fill(w_OS,sw);
            eff_OS_norm += sw;
        }
        if(q_SS != 0){
            h_w_SS_norm->Fill(w_SS,sw);
            eff_SS_norm += sw;
            }
        sumw += sw;
    }
    
    eff_OS_norm /= sumw;    
    eff_SS_norm /= sumw;

    ///loop over MC events
    for(int i=0; i< treeMC->GetEntries(); i++)
    {	
        treeMC->GetEntry(i);
        
        h_dt_MC->Fill(dt_MC,sw);
        h_q_OS_MC->Fill((double)q_OS_MC,sw);
        h_q_SS_MC->Fill((double)q_SS_MC,sw);
        
        if(q_OS_MC != 0)h_w_OS_MC->Fill(w_OS_MC,w);
        if(q_SS_MC != 0)h_w_SS_MC->Fill(w_SS_MC,w);
    }
    
    for(int i=0; i< treeMC_norm->GetEntries(); i++)
    {	
        treeMC_norm->GetEntry(i);
        
        h_dt_MC_norm->Fill(dt_MC,sw);
        h_q_OS_MC_norm->Fill((double)q_OS_MC,sw);
        h_q_SS_MC_norm->Fill((double)q_SS_MC,sw);
        
        if(q_OS_MC != 0)h_w_OS_MC_norm->Fill(w_OS_MC,w);
        if(q_SS_MC != 0)h_w_SS_MC_norm->Fill(w_SS_MC,w);
    }
    
    ///Plot it
    TCanvas* c= new TCanvas();

    h_w_OS->SetMinimum(0);    
    h_w_OS->SetLineColor(kBlack);
    h_w_OS->DrawNormalized("e",1);
    h_w_OS_norm->SetMarkerColor(kRed);
    h_w_OS_norm->SetLineColor(kRed);
    h_w_OS_norm->DrawNormalized("esame",1);
    
    double KolmoTest = h_w_OS->KolmogorovTest(h_w_OS_norm);
    
    TLegend *leg = new TLegend(0.2,0.6,0.4,0.9,"");
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextFont(22);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    leg->SetTextAlign(12);
    
    leg->AddEntry(h_w_OS,"B_{s} #rightarrow D_{s}K#pi#pi  Data","LEP");
    
    stringstream ss1 ;
    TString leg_average = "<#omega_{OS}> = ";
    ss1 << std::fixed << std::setprecision(3) << h_w_OS->GetMean() ;
    leg_average += ss1.str();    
    leg->AddEntry((TObject*)0, leg_average, "");

    stringstream ss2 ;
    TString leg_eff = "#epsilon_{OS} = ";
    ss2 << std::fixed << std::setprecision(3) << eff_OS ;
    leg_eff += ss2.str();    
    leg->AddEntry((TObject*)0, leg_eff, "");

    leg->AddEntry(h_w_OS_norm,"B_{s} #rightarrow D_{s}#pi#pi#pi  Data","LEP");
    
    stringstream ss ;
    TString leg_kol = "Kolm.-Test : ";
    ss << std::fixed << std::setprecision(2) << KolmoTest ;
    leg_kol += ss.str();    
    //TLegendEntry* le = leg->AddEntry((TObject*)0, leg_kol, "");
    //le->SetTextColor(kRed);    
    
    //leg->Draw(); 
    c->Print(prefix+"w_OS.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/w_OS.pdf");

    h_w_SS->SetMinimum(0);    
    h_w_SS->SetLineColor(kBlack);
    h_w_SS->DrawNormalized("e",1);
    h_w_SS_norm->SetMarkerColor(kRed);
    h_w_SS_norm->SetLineColor(kRed);
    h_w_SS_norm->DrawNormalized("esame",1);
    c->Print(prefix+"w_SS.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/w_SS.pdf");

    h_w_OS_MC->SetMinimum(0);    
    h_w_OS_MC->SetLineColor(kBlack);
    h_w_OS_MC->DrawNormalized("e",1);
    h_w_OS_MC_norm->SetMarkerColor(kRed);
    h_w_OS_MC_norm->SetLineColor(kRed);
    h_w_OS_MC_norm->DrawNormalized("esame",1);
    c->Print(prefix+"w_OS_MC.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/w_OS_MC.pdf");

    h_w_SS_MC->SetMinimum(0);    
    h_w_SS_MC->SetLineColor(kBlack);
    h_w_SS_MC->DrawNormalized("e",1);
    h_w_SS_MC_norm->SetMarkerColor(kRed);
    h_w_SS_MC_norm->SetLineColor(kRed);
    h_w_SS_MC_norm->DrawNormalized("esame",1);
    c->Print(prefix+"w_SS_MC.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/w_SS_MC.pdf");

    h_q_OS->SetMinimum(0);    
    h_q_OS->SetLineColor(kBlack);
    h_q_OS->DrawNormalized("e",1);
    h_q_OS_norm->SetMarkerColor(kRed);
    h_q_OS_norm->SetLineColor(kRed);
    h_q_OS_norm->DrawNormalized("esame",1);
    c->Print(prefix+"q_OS.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/qOS.pdf");

    h_q_SS->SetMinimum(0);    
    h_q_SS->SetLineColor(kBlack);
    h_q_SS->DrawNormalized("e",1);
    h_q_SS_norm->SetMarkerColor(kRed);
    h_q_SS_norm->SetLineColor(kRed);
    h_q_SS_norm->DrawNormalized("esame",1);
    c->Print(prefix+"q_SS.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/q_SS.pdf");
    
    h_q_OS_MC->SetMinimum(0);    
    h_q_OS_MC->SetLineColor(kBlack);
    h_q_OS_MC->DrawNormalized("e",1);
    h_q_OS_MC_norm->SetMarkerColor(kRed);
    h_q_OS_MC_norm->SetLineColor(kRed);
    h_q_OS_MC_norm->DrawNormalized("esame",1);
    c->Print(prefix+"q_OS_MC.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/q_OS_MC.pdf");

    h_q_SS_MC->SetMinimum(0);    
    h_q_SS_MC->SetLineColor(kBlack);
    h_q_SS_MC->DrawNormalized("e",1);
    h_q_SS_MC_norm->SetMarkerColor(kRed);
    h_q_SS_MC_norm->SetLineColor(kRed);
    h_q_SS_MC_norm->DrawNormalized("esame",1);
    c->Print(prefix+"q_SS_MC.eps");
    if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Tagging/q_SS_MC.pdf");

    h_dt->SetMinimum(0);    
    h_dt->SetLineColor(kBlack);
    h_dt->DrawNormalized("e",1);
    h_dt_norm->SetMarkerColor(kRed);
    h_dt_norm->SetLineColor(kRed);
    h_dt_norm->DrawNormalized("esame",1);
    c->Print("dt.eps");
    
    h_dt_MC->SetMinimum(0);    
    h_dt_MC->SetLineColor(kBlack);
    h_dt_MC->DrawNormalized("e",1);
    h_dt_MC_norm->SetMarkerColor(kRed);
    h_dt_MC_norm->SetLineColor(kRed);
    h_dt_MC_norm->DrawNormalized("esame",1);
    c->Print("dt_MC.eps");
    
    TFile* out = new TFile("Mistag_pdfs.root","RECREATE");
    h_dt_norm->Write();
    h_q_OS_norm->Write();
    h_w_OS_norm->Write();
    h_q_SS_norm->Write();
    h_w_SS_norm->Write();
    out->Write();

}

void produceIntegratorFile_CP(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;

    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
    DiskResidentEventList eventList(((string) integratorEventFile.ReplaceAll(".root","_CP.root")).c_str(),"RECREATE");

    for(int i = 0; i < eventListMC.size(); i++){
	DalitzEvent evt(eventListMC.getEvent(i));
	evt.P_conjugateYourself();
	eventList.Add(evt);
    }

    eventList.save();
    return;
}


int main(int argc, char** argv){

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->ProcessLine(".x ../lhcbStyle.C");

  //produceMarginalPdfs();

  //produceIntegratorFile_CP();

  ampFit();
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
