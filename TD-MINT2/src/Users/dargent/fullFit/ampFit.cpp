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
#include "RooAddPdf.h"
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
#include "Mint/HyperHistogram.h"
#include "Mint/LASSO.h"
#include "Mint/LASSO_flexi.h"

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
    AmpRatio(FitParameter& re, FitParameter& im,int f = 1)
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
    
    CPV_amp(FitParameter& re, FitParameter& im, int sign, const IMinuitParameter* scale_re, const IMinuitParameter* scale_im)
    : _re(re), _im(im), _sign(sign) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result((double) ( 1.+  _re * (double) _sign),(double) (_im * (double) _sign) ); 
        return result;
    }
};

class CPV_amp_norm : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    const IMinuitParameter* _scale_re;
    const IMinuitParameter* _scale_im;
    int _sign;
public:
    CPV_amp_norm(FitParameter& re, FitParameter& im, int sign, const IMinuitParameter* scale_re, const IMinuitParameter* scale_im)
    : _re(re), _im(im), _sign(sign), _scale_re(scale_re), _scale_im(scale_im) {}
    
    complex<double> ComplexVal(){
        double norm = sqrt(pow(_scale_re->mean(),2)+pow(_scale_im->mean(),2));
        std::complex<double> result((double) ( 1.+  _re/norm * (double) _sign),(double) (_im/norm * (double) _sign) );
        return result;
    }
};

class CPV_amp_norm_scaled : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    const IMinuitParameter* _scale_re;
    const IMinuitParameter* _scale_im;
    int _sign;
public:
    CPV_amp_norm_scaled(FitParameter& re, FitParameter& im, int sign, const IMinuitParameter* scale_re, const IMinuitParameter* scale_im)
    : _re(re), _im(im), _sign(sign), _scale_re(scale_re), _scale_im(scale_im) {}
    
    complex<double> ComplexVal(){
        double norm = sqrt(pow(_scale_re->mean(),2)+pow(_scale_im->mean(),2));
        std::complex<double> scale(_scale_re->mean(),_scale_im->mean() );
        std::complex<double> result((double) ( 1.+  _re/norm * (double) _sign),(double) (_im/norm * (double) _sign) );
        return scale*result;
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
// 	    std::complex<double> result= polar((double) ( 1.+  _r * (double) _sign),(double) (_delta/360.*2.*pi * (double) _sign) ); 
	    std::complex<double> result= polar((double) ( _r ),(double) (_delta/360.*2.*pi ) ); 
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

void AddScaledAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name, counted_ptr<IReturnComplex>& scale){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.multiply(scale);
    fas.addAsList(fas_2,1.);
}

void AddAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.addAsList(fas_2,1.);
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

    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
    
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "phase = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "D = " << result[4] << endl;
    cout << "S = " << result[5] << endl;

    return result;
}

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
        DalitzEvent evt = eventList.getEvent(i);

        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
    }
    
    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);

    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
    
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "D = " << result[4] << endl;
    cout << "S = " << result[5] << endl;
    
    return result;
}

std::vector<double> coherenceFactor_CP(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
        DalitzEvent evt = eventList.getEvent(i);
        evt.CP_conjugateYourself();
        
        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
    }
    
    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);
    
    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );

    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "Dbar = " << result[4]<< endl;
    cout << "Sbar = " << result[5] << endl;
    
    return result;
}

std::vector<double> coherenceFactor_CP(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
        cout << "Calculating coherence factor ..." << endl << endl;
        fas.print();
        fas_bar.print();
    
        std::complex<double> valK(0,0);
        double val1 = 0;
        double val2 = 0;
    
        const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
        for(unsigned int i=0; i<eventList.size(); i++){
                DalitzEvent evt(eventList[i]);
                evt.CP_conjugateYourself();
                const std::complex<double> amp = fas.getVal(evt) ;
                const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
                valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
                val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
                val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
            }
	
	std::vector<double> result;
	result.push_back(sqrt(val2/val1));
	result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
	result.push_back(std::arg(valK)/(2.*pi)*360.);
	
	result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
	result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
	result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
	
	cout << "r = " << result[0] << endl;
	cout << "k = " << result[1] << endl;
	cout << "d = " << result[2] << " [deg]" << endl << endl;
	
	cout << "C = " << result[3] << endl;
	cout << "Dbar = " << result[4]<< endl;
	cout << "Sbar = " << result[5] << endl;
    
        return result;
}

// Full time-dependent PDF
class FullAmpsPdfFlexiFastCPV : public MINT::PdfBase<IDalitzEvent>
, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps;
    AmpsPdfFlexiFast* _amps_bar;

    AmpsPdfFlexiFast* _amps_CP;
    AmpsPdfFlexiFast* _amps_bar_CP;

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

    const MINT::FitParameter& _tau;
    const MINT::FitParameter& _dGamma;
    const MINT::FitParameter& _dm;
    
    const MINT::FitParameter& _offset_sigma_dt;    
    const MINT::FitParameter& _scale_mean_dt;
    const MINT::FitParameter& _scale_sigma_dt;
    const MINT::FitParameter& _scale_sigma_2_dt;

    const MINT::FitParameter& _c0;
    const MINT::FitParameter& _c1;
    const MINT::FitParameter& _c2;
    const MINT::FitParameter& _c3;
    const MINT::FitParameter& _c4;
    const MINT::FitParameter& _c5;
    const MINT::FitParameter& _c6;
    const MINT::FitParameter& _c7;
    const MINT::FitParameter& _c8;
    const MINT::FitParameter& _c9;
    
    const MINT::FitParameter& _p0_os;
    const MINT::FitParameter& _p1_os;
    const MINT::FitParameter& _delta_p0_os;
    const MINT::FitParameter& _delta_p1_os;
    const MINT::FitParameter& _avg_eta_os;
    const MINT::FitParameter& _tageff_os;
    const MINT::FitParameter& _tageff_asym_os;
    
    const MINT::FitParameter& _p0_ss;
    const MINT::FitParameter& _p1_ss;
    const MINT::FitParameter& _delta_p0_ss;
    const MINT::FitParameter& _delta_p1_ss;
    const MINT::FitParameter& _avg_eta_ss;
    const MINT::FitParameter& _tageff_ss;
    const MINT::FitParameter& _tageff_asym_ss;
    
    const MINT::FitParameter& _production_asym;
    const MINT::FitParameter& _detection_asym;
    
    string _marginalPdfsPrefix;
    // Time pdf master
    TimePdfMaster* _timePdfMaster;


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
        _timePdfMaster->listFitParDependencies();
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

    inline double getValForGeneration(IDalitzEvent& evt){

        const double t = (double) evt.getValueFromVector(0);
        const double f = (double) evt.getValueFromVector(2);
        //if(t < _min_TAU || t > _max_TAU )return 0.;
        _timePdfMaster->setAllObservablesAndFitParameters(evt);
        
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_CP = polar((double)r,((double) _delta +(double)_gamma)/360.*2*pi);
        
        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
        double norm_amp = 1.;
        
        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 1.;
        
        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt) ;
            amp_bar = _amps_bar->ComplexVal_un_normalised_noPs(evt) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt) ;
            amp_bar_CP = _amps_bar_CP->ComplexVal_un_normalised_noPs(evt) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
                                   norm_amp,
                                   norm_amp_CP,
                                   ( norm(amp) - norm(amp_bar) ) ,
                                   -( norm(amp_CP) - norm(amp_bar_CP) ),
                                   -2.* real(amp_bar*conj(amp) * phase_diff) ,
                                   2.* real(amp_bar_CP*conj(amp_CP) * phase_diff_CP) ,
                                   2. * imag(amp_bar*conj(amp) * phase_diff) ,
                                   2. * imag(amp_bar_CP*conj(amp_CP) * phase_diff_CP) );
        
        const double tau = _timePdfMaster->get_tau_Val();
        const double dGamma = _timePdfMaster->get_dGamma_Val();
        const double dm = _timePdfMaster->get_dm_Val();

        const double val =  exp(-fabs(t)/tau) *
        ( _timePdfMaster->get_cosh_coeff_Val(evt) *cosh(dGamma/2.*t)
         +  _timePdfMaster->get_cos_coeff_Val(evt) *cos(dm*t)
         +  _timePdfMaster->get_sinh_coeff_Val(evt) *sinh(dGamma/2.*t)
         +  _timePdfMaster->get_sin_coeff_Val(evt) *sin(dm*t)
         )
        * _timePdfMaster->get_spline_Val(evt)
        * _timePdfMaster->get_marginalPdfs_Val(evt);

        /*
        cout << _timePdfMaster.get_marginalPdfs_Val(evt) << endl;
        cout << _timePdfMaster.get_spline_Val(evt) << endl << endl;

        cout << norm_amp << endl;
        cout << _timePdfMaster.get_cosh_coeff_Val(evt) << endl;
        cout << _timePdfMaster.get_cos_coeff_Val(evt) << endl;
        cout << 2.* real(amp_bar*conj(amp) * phase_diff_CP) << endl;
        cout <<  2. * imag(amp_bar*conj(amp) * phase_diff_CP) << endl;

        throw "";
        */
        return val;
    }

    inline double un_normalised_noPs(IDalitzEvent& evt){

        //const double t = (double) evt.getValueFromVector(0);
        const double f = (double) evt.getValueFromVector(2);
        //if(t < _min_TAU || t > _max_TAU )return 0.;
        _timePdfMaster->setAllObservablesAndFitParameters(evt);
                
        double r = (double)_r; // * sqrt(_intA/_intAbar);
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_CP = polar((double)r,((double) _delta +(double)_gamma)/360.*2*pi);

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 1.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 1.;

        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt) ;
            amp_bar = _amps_bar->ComplexVal_un_normalised_noPs(evt) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt) ;
            amp_bar_CP = _amps_bar_CP->ComplexVal_un_normalised_noPs(evt) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }

        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp) * phase_diff) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP) * phase_diff_CP) ,
        		 2. * imag(amp_bar*conj(amp) * phase_diff) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP) * phase_diff_CP) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
             ( _timePdfMaster->get_cosh_term_Val(evt)
             +  _timePdfMaster->get_cos_term_Val(evt)
             +  _timePdfMaster->get_sinh_term_Val(evt)
             +  _timePdfMaster->get_sin_term_Val(evt)
             ) * _timePdfMaster->get_marginalPdfs_Val(evt);

        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double val = un_normalised_noPs(evt);
        double r = (double)_r; // * sqrt(_intA/_intAbar);
    
        const complex<double> phase_diff = polar((double)r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_CP = polar((double)r,((double) _delta +(double)_gamma)/360.*2*pi);

        const complex<double> int_interference =  phase_diff * _intAAbar ;
        const complex<double> int_interference_CP = phase_diff_CP * _intAAbar_CP ;

        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
			(_intA + r* r * _intAbar),
			(_intA_CP + r* r * _intAbar_CP),
        		(_intA - r* r * _intAbar),
        		-(_intA_CP - r* r * _intAbar_CP) , 
        		(- int_interference.real() ), /// *2 included in integration !!!
        		(- int_interference_CP.real() ), 
        		(int_interference.imag() ),  /// - sign included in DecRateCoeff !!!
        		(- int_interference_CP.imag() )   
		 ); 
        
        double norm =  // (_intA + r* r * _intAbar)  *
        (     _timePdfMaster->get_cosh_term_Integral(evt)
            +  _timePdfMaster->get_cos_term_Integral(evt)
            +  _timePdfMaster->get_sinh_term_Integral(evt)
            +  _timePdfMaster->get_sin_term_Integral(evt) );

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

    vector<double> coherenceFactor(bool print = true){
	
	if(_intA == -1){
		this->parametersChanged();
	}

	double val1 = _intA;
	double val2 =  _r* _r * _intAbar;

	double val1_CP = _intA_CP;
	double val2_CP =  _r* _r * _intAbar_CP;

        const complex<double> phase_diff = polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
        const complex<double> phase_diff_CP = polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);

        const complex<double> int_interference =  phase_diff * _intAAbar ;
        const complex<double> int_interference_CP = phase_diff_CP * _intAAbar_CP ;

	complex<double> valK = int_interference/2.;
	complex<double> valK_CP = int_interference_CP/2.;
	
	std::vector<double> result;

	result.push_back(sqrt(val2/val1));
	result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
	result.push_back(std::arg(valK)/(2.*pi)*360.);

	result.push_back(sqrt(val2_CP/val1_CP));
	result.push_back(std::abs(valK_CP)/sqrt(val1_CP)/sqrt(val2_CP));
	result.push_back(std::arg(valK_CP)/(2.*pi)*360.);

	result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
	result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
	result.push_back(2. * valK.imag()/val1 /(1.+pow(result[0],2)) );

	result.push_back(-(1.-pow(result[3],2))/(1.+pow(result[3],2)));
	result.push_back(-2. * valK_CP.real()/val1_CP /(1.+pow(result[3],2)) );
	result.push_back(-2. * valK_CP.imag()/val1_CP /(1.+pow(result[3],2)) );
	
	if(print){
		cout << "r = " << result[0] << endl;
		cout << "k = " << result[1] << endl;
		cout << "phase = " << result[2] << " [deg]" << endl << endl;
		cout << "r_CP = " << result[3] << endl;
		cout << "k_CP = " << result[4] << endl;
		cout << "phase_CP = " << result[5] << " [deg]" << endl << endl;
		
		cout << "C = " << result[6] << endl;
		cout << "D = " << result[7] << endl;
		cout << "S = " << result[8] << endl;
		cout << "Cbar = " << result[9] << endl;
		cout << "Dbar = " << result[10] << endl;
		cout << "Sbar = " << result[11] << endl;
	}
	return result;
    }

    double getCPcoeffChi2(vector<double>& coeff,double prec){
	vector<double> val = coherenceFactor(false);
	double chi2  = 0.;
	for(int i= 0; i< coeff.size();i++){
		chi2+= pow((coeff[i]-val[i+6])/prec,2);
	}
	return chi2;
    }

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_SS(evt);
    }
    
    double getCalibratedResolution(double& dt){
        return _timePdfMaster->getCalibratedResolution(dt);
    }
    
    virtual DalitzHistoSet histoSet(){return _ampsSum->histoSet();}
    
    void doFinalStatsAndSaveForAmp12(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults", const std::string& fnameROOT="fitFractions"){
        _amps->redoIntegrator();
        _amps_bar->redoIntegrator();
        _amps->doFinalStatsAndSave(min,((string)fname+".tex").c_str(),((string)fnameROOT+".root").c_str());
        _amps_bar->doFinalStatsAndSave(min,((string)fname+"_bar.tex").c_str(),((string)fnameROOT+"_Bar.root").c_str());        
    }
    
    FullAmpsPdfFlexiFastCPV(
		AmpsPdfFlexiFast* amps, AmpsPdfFlexiFast* amps_bar, 
		AmpsPdfFlexiFast* amps_CP, AmpsPdfFlexiFast* amps_bar_CP,
		AmpsPdfFlexiFast* ampsSum, AmpsPdfFlexiFast* ampsSum_CP, 
                MINT::FitParameter& r,MINT::FitParameter& delta, MINT::FitParameter& gamma,
		const MINT::FitParameter& tau, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
                ,const MINT::FitParameter& offset_sigma_dt, const MINT::FitParameter& scale_mean_dt 
		,const MINT::FitParameter& scale_sigma_dt, const MINT::FitParameter& scale_sigma_2_dt
                ,const MINT::FitParameter& c0, const MINT::FitParameter& c1, const MINT::FitParameter& c2
                ,const MINT::FitParameter& c3, const MINT::FitParameter& c4, const MINT::FitParameter& c5
                ,const MINT::FitParameter& c6, const MINT::FitParameter& c7, const MINT::FitParameter& c8
                ,const MINT::FitParameter& c9,
                const MINT::FitParameter& p0_os, const MINT::FitParameter& p1_os, const MINT::FitParameter& delta_p0_os, const MINT::FitParameter& delta_p1_os, 
                const MINT::FitParameter& avg_eta_os, const MINT::FitParameter& tageff_os, const MINT::FitParameter& tageff_asym_os, 
                const MINT::FitParameter& p0_ss, const MINT::FitParameter& p1_ss, const MINT::FitParameter& delta_p0_ss, const MINT::FitParameter& delta_p1_ss, 
                const MINT::FitParameter& avg_eta_ss, const MINT::FitParameter& tageff_ss, const MINT::FitParameter& tageff_asym_ss, 
                const MINT::FitParameter& production_asym, const MINT::FitParameter& detection_asym, string marginalPdfsPrefix = ""
                ):
    _amps(amps),_amps_bar(amps_bar),_amps_CP(amps_CP),_amps_bar_CP(amps_bar_CP),_ampsSum(ampsSum),_ampsSum_CP(ampsSum_CP),
    _intA(-1),_intAbar(-1),_intAAbar(-1),_intA_CP(-1),_intAbar_CP(-1),_intAAbar_CP(-1),
    _r(r),_delta(delta),_gamma(gamma),
    _tau(tau),
    _dGamma(dGamma),
    _dm(dm),
    _offset_sigma_dt(offset_sigma_dt),    
    _scale_mean_dt(scale_mean_dt),
    _scale_sigma_dt(scale_sigma_dt),
    _scale_sigma_2_dt(scale_sigma_2_dt),
    _c0(c0),
    _c1(c1),
    _c2(c2),
    _c3(c3),
    _c4(c4),
    _c5(c5),
    _c6(c6),
    _c7(c7),
    _c8(c8),
    _c9(c9),
    _p0_os(p0_os),
    _p1_os(p1_os),
    _delta_p0_os(delta_p0_os),
    _delta_p1_os(delta_p1_os),
    _avg_eta_os(avg_eta_os),
    _tageff_os(tageff_os),
    _tageff_asym_os(tageff_asym_os),
    _p0_ss(p0_ss),
    _p1_ss(p1_ss),
    _delta_p0_ss(delta_p0_ss),
    _delta_p1_ss(delta_p1_ss),
    _avg_eta_ss(avg_eta_ss),
    _tageff_ss(tageff_ss),
    _tageff_asym_ss(tageff_asym_ss),
    _production_asym(production_asym),
    _detection_asym(detection_asym),
    _marginalPdfsPrefix(marginalPdfsPrefix)
    {
		        _timePdfMaster = new TimePdfMaster(_tau, _dGamma, _dm
                                          ,_offset_sigma_dt, _scale_mean_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                          ,_c0, _c1, _c2
                                          ,_c3, _c4, _c5
                                          ,_c6, _c7, _c8
                                          ,_c9,
                                          _p0_os, _p1_os, _delta_p0_os, _delta_p1_os, 
                                          _avg_eta_os, _tageff_os, _tageff_asym_os, 
                                          _p0_ss, p1_ss, _delta_p0_ss, _delta_p1_ss, 
                                          _avg_eta_ss, _tageff_ss, _tageff_asym_ss, 
                                          _production_asym, _detection_asym,_marginalPdfsPrefix);
    }
};

class CPcoeffLL : public Minimisable{
  FullAmpsPdfFlexiFastCPV* _pdf;
  vector<double> coeff;
  double _prec;
  double _C,_Cbar,_S,_Sbar,_D,_Dbar;
  double _r,_k,_delta,_gamma;
public:
  CPcoeffLL(FullAmpsPdfFlexiFastCPV* pdf,double C,double Cbar,double S,double Sbar,double D,double Dbar, double prec) 
	: _pdf(pdf),_C(C),_Cbar(Cbar),_D(D),_Dbar(Dbar),_S(S),_Sbar(Sbar),_prec(prec){
	    _pdf->parametersChanged();//makes sure we are initialised
	    coeff.push_back(_C);
	    coeff.push_back(_Cbar);
	    coeff.push_back(_D);
	    coeff.push_back(_Dbar);
	    coeff.push_back(_S);
	    coeff.push_back(_Sbar);
  }
  double getVal(){
    _pdf->parametersChanged();
    return _pdf->getCPcoeffChi2(coeff,_prec);
  }
};

class FracLL : public Minimisable{
  FlexiFastAmplitudeIntegrator* _integ;
public:
  FracLL(FlexiFastAmplitudeIntegrator* integ) : _integ(integ){
    _integ->getVal();//makes sure we are initialised
  }
  double getVal(){
    return _integ->getFractionChi2();
  }
};

double getChi2(DalitzEventList& data, DiskResidentEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 5;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 25);
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4));
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4));
    //HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),4 * GeV * GeV);

    HyperCuboid limits(min, max );

    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI,
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (0)
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),                      
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")                 
                         );

//    dataHist.save("histData.root");
//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
	//if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //mcHist.save("histMC.root");
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    //mcHist.divide(dataHist);
    //mcHist.save("histMC_Data.root");

    return (double)chi2/(nBins-1.);
}

double getChi2_6D(DalitzEventList& data, DiskResidentEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 6;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
          
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);

    NamedParameter<int> minEventsPerBin("minEventsPerBin", 50);       
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4),(double)min_TAU);
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4),(double)max_TAU);
    //HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),4 * GeV * GeV);

    HyperCuboid limits(min, max );

    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.at(5)= evt.getValueFromVector(0);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART,  
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (4)                        
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),                      
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")                 
                         );

//    dataHist.save("histData.root");
//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.at(5)= evt.getValueFromVector(0);
	point.addWeight(evt.getWeight());
	//if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //mcHist.save("histMC.root");
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    //mcHist.divide(dataHist);
    //mcHist.save("histMC_Data.root");

    return (double)chi2/(nBins-1.);
}

double cosThetaAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}


void ampFit(int step=0){

    /// Options
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
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
    DalitzEventPattern pat_CP = pat.makeCPConjugate();

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);
    NamedParameter<string> InputGenMCFile("InputGenMCFile", (std::string) "/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_DsKpipi_CPV.root", (char*) 0);

    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<int>  generateNew("generateToys", 0);

    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);
    NamedParameter<double> w_max("w_max", 0.5);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  nBins("nBins", 50);
    NamedParameter<int>  nBinst("nBinst", 50);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 8);
    NamedParameter<int>  randomizeStartVals("randomizeStartVals", 0);
    NamedParameter<int>  doSimFit("doSimFit", 0);
    NamedParameter<int>  constrainFF("constrainFF", 0);
    NamedParameter<int>  useLASSO("useLASSO", 1);
    NamedParameter<double>  lambda("lambda", 1.);

    NamedParameter<int>  initCPcoeff("initCPcoeff", 0);
    NamedParameter<int>  fitGenMC("fitGenMC", 0);
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;
    TString integratorEventFile_CP = (string) IntegratorEventFile;
    integratorEventFile_CP.ReplaceAll(".root","_CP.root");

    NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
    NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");
    
    NamedParameter<int>  Nevents("Nevents", 100);
    NamedParameter<double>  pdf_max("pdf_max", 100);

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

    /// Define amplitude model
    DalitzEventList eventListPhsp,eventListPhsp_CP;
    eventListPhsp.generatePhaseSpaceEvents(200000,pat);
    eventListPhsp_CP.generatePhaseSpaceEvents(200000,pat_CP);

    FitAmpSum fas_tmp((DalitzEventPattern)pat);
    //fas_tmp.normalizeAmps(eventListPhsp);
    //if(randomizeStartVals)fas_tmp.randomizeStartVals(seed);
    //if(randomizeStartVals)fas_tmp.randomizePhaseStartVals(seed);
    
    /// Normalize amps
    {
    	DalitzEventList eventListNorm;
  	TFile *file =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
  	TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
  	eventListNorm.fromNtuple(tree,0.5);
  	fas_tmp.normalizeAmps(eventListNorm);
    }

    MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
    
    /// Choose reference amp
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1400)+");
    FitAmpSum fas(*List_1);
    FitAmpSum fas_bar(*List_1);
    FitAmpSum fas_tot(*List_1);
    FitParameter r_1_re("r_1_Re",2,0,0.01);
    FitParameter r_1_im("r_1_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_1_plus = new CPV_amp(r_1_re,r_1_im,1);
    counted_ptr<IReturnComplex> r_1_minus = new CPV_amp(r_1_re,r_1_im,-1);
    fas.multiply(r_1_plus); 
    fas_bar.multiply(r_1_minus);
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "K(1)(1270)+", r_K1_plus, r_K1_minus );

    /// Fix relative decay modes from charm input
    FitParameter a_K1_1270_re("a_K1_1270_Re",2,1,0.01);
    FitParameter a_K1_1270_im("a_K1_1270_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_K1_1270 = new AmpRatio(a_K1_1270_re,a_K1_1270_im);

    FitParameter a_Ks_1410_re("a_Ks_1410_Re",2,1,0.01);
    FitParameter a_Ks_1410_im("a_Ks_1410_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_Ks_1410 = new AmpRatio(a_Ks_1410_re,a_Ks_1410_im);

    FitParameter a_K_1460_re("a_K_1460_Re",2,1,0.01);
    FitParameter a_K_1460_im("a_K_1460_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_K_1460 = new AmpRatio(a_K_1460_re,a_K_1460_im);

    /// Add amps to A
    AddScaledAmpsToList(fas_tmp, fas, "K(1)(1270)+",a_K1_1270);
    AddScaledAmpsToList(fas_tmp, fas, "K*(1410)+",a_Ks_1410);
    //AddAmpsToList(fas_tmp, fas, "K*(1410)+");
    AddAmpsToList(fas_tmp, fas, "NonResS0(->Ds-,K+),sigma10(->pi+,pi-)");
    //AddAmpsToList(fas_tmp, fas, "NonResV0(->Ds-,K+),rho(770)0(->pi+,pi-)");

    /// Add amps to Abar
    AddScaledAmpsToList(fas_tmp, fas_bar, "K(1460)",a_K_1460);
    AddAmpsToList(fas_tmp, fas_bar, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)");

    /// Add amp to Atot
    AddScaledAmpsToList(fas_tmp, fas_tot, "K(1)(1270)+",a_K1_1270);
    AddScaledAmpsToList(fas_tmp, fas_tot, "K(1460)",a_K_1460);
    AddScaledAmpsToList(fas_tmp, fas_tot, "K*(1410)+",a_Ks_1410);
    AddAmpsToList(fas_tmp, fas_tot, "NonRes");

    /*
    const IMinuitParameter* mp_K1410_Re = mps->getParPtr("Bs0->K*(1410)+(->K*(892)0(->K+,pi-),pi+),Ds-_Re");
    const IMinuitParameter* mp_K1410_Im = mps->getParPtr("Bs0->K*(1410)+(->K*(892)0(->K+,pi-),pi+),Ds-_Im");
    FitParameter r_2_re("r_2_Re",2,0,0.01);
    FitParameter r_2_im("r_2_Im",2,0,0.01);
    counted_ptr<IReturnComplex> r_2_plus = new CPV_amp_norm(r_2_re,r_2_im,1,mp_K1410_Re,mp_K1410_Im);
    counted_ptr<IReturnComplex> r_2_minus = new CPV_amp_norm(r_2_re,r_2_im,-1,mp_K1410_Re,mp_K1410_Im);
    counted_ptr<IReturnComplex> r_2_plus_scaled = new CPV_amp_norm_scaled(r_2_re,r_2_im,1,mp_K1410_Re,mp_K1410_Im);
    counted_ptr<IReturnComplex> r_2_minus_scaled = new CPV_amp_norm_scaled(r_2_re,r_2_im,-1,mp_K1410_Re,mp_K1410_Im);
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "Bs0->K*(1410)+(->K*(892)0(->K+,pi-),pi+),Ds-", r_2_plus, r_2_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "Bs0->K*(1410)+(->rho(770)0(->pi+,pi-),K+),Ds-", r_2_plus_scaled, r_2_minus_scaled );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "K*(1410)+", r_2_plus, r_2_minus );

    const IMinuitParameter* mp_3_Re = mps->getParPtr("NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)_Re");
    const IMinuitParameter* mp_3_Im = mps->getParPtr("NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)_Im");
    FitParameter r_3_re("r_3_Re",2,0,0.01);
    FitParameter r_3_im("r_3_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_3_plus = new CPV_amp(r_3_re,r_3_im,1);
    counted_ptr<IReturnComplex> r_3_minus = new CPV_amp(r_3_re,r_3_im,-1);
    counted_ptr<IReturnComplex> r_3_plus_scaled = new CPV_amp_norm_scaled(r_3_re,r_3_im,1,mp_3_Re,mp_3_Im);
    counted_ptr<IReturnComplex> r_3_minus_scaled = new CPV_amp_norm_scaled(r_3_re,r_3_im,-1,mp_3_Re,mp_3_Im);
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "BgSpinZeroBs0->NonResS0(->Ds-,K+),NonResS0(->pi+,pi-)", r_3_plus_scaled, r_3_minus_scaled );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResS0(->Ds-,K+),sigma10(->pi+,pi-)", r_3_plus_scaled, r_3_minus_scaled );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResV0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "K(1460)+(->K*(892)0(->K+,pi-),pi+),Ds-", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "K*(1410)+", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResA0(->sigma10(->pi+,pi-),Ds-)", r_3_plus, r_3_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResA0(->rho(770)0(->pi+,pi-),Ds-),K+", r_3_plus, r_3_minus );   

    FitParameter r_4_re("r_4_Re",2,0,0.01);
    FitParameter r_4_im("r_4_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> r_4_plus = new CPV_amp(r_4_re,r_4_im,1);
    counted_ptr<IReturnComplex> r_4_minus = new CPV_amp(r_4_re,r_4_im,-1);
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "K(1)(1400)+", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResA0(->sigma10(->pi+,pi-),Ds-)", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResA0(->rho(770)0(->pi+,pi-),Ds-),K+", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "BgSpinZeroBs0->NonResS0(->Ds-,K+),NonResS0(->pi+,pi-)", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResV0(->Ds-,K+),sigma10(->pi+,pi-)", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResS0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_4_plus, r_4_minus );
    //AddScaledAmpsToList(fas_tmp, fas, fas_bar, "NonResV0(->Ds-,pi+),K*(892)0(->K+,pi-)", r_4_plus, r_4_minus );
    */
   
    /// Define B -> f amplitude        
    fas.setTag(1);
    /// Define Bbar -> f amplitude
    fas_bar.setTag(-1);

    /// CP conjugate amplitudes
    FitAmpSum fas_CP(fas);
    fas_CP.CPConjugateSameFitParameters();

    FitAmpSum fas_bar_CP(fas_bar);
    fas_bar_CP.CPConjugateSameFitParameters();

    /// Multiply r e^(i gamma)
    FitParameter x_m("xm",2,0,0.01);
    FitParameter y_m("ym",2,0,0.01); 
    FitParameter x_p("xp",2,0,0.01);
    FitParameter y_p("yp",2,0,0.01); 

    counted_ptr<IReturnComplex> xy_m = new AmpRatio(x_m,y_m,1);
    counted_ptr<IReturnComplex> xy_p = new AmpRatio(x_p,y_p,1);

    //counted_ptr<IReturnComplex> xy_m = new CPV_amp_polar(x_m,y_m,1);
    //counted_ptr<IReturnComplex> xy_p = new CPV_amp_polar(x_m,y_p,1);

    fas_bar.multiply(xy_m); 
    fas_bar_CP.multiply(xy_p); 

    /// Add amplitudes: A + r e^(i gamma) Abar
    counted_ptr<FitAmpList> sumList = fas.GetCloneSameFitParameters();
    FitAmpSum fas_sum(*sumList);
    fas_sum.addAsList(fas_bar,1.);
    fas_sum.getVal(eventListPhsp[0]);
    
    AmpsPdfFlexiFast ampsSig(pat, &fas, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSig_bar(pat, &fas_bar, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSum(pat, &fas_sum, 0, integPrecision,integMethod, (std::string) integratorEventFile);

    counted_ptr<FitAmpList> sumList_CP = fas_CP.GetCloneSameFitParameters();
    FitAmpSum fas_sum_CP(*sumList_CP);
    fas_sum_CP.addAsList(fas_bar_CP,1.);
    fas_sum_CP.getVal(eventListPhsp_CP[0]);
    
    AmpsPdfFlexiFast ampsSig_CP(pat_CP, &fas_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);
    AmpsPdfFlexiFast ampsSig_bar_CP(pat_CP, &fas_bar_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);
    AmpsPdfFlexiFast ampsSum_CP(pat_CP, &fas_sum_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);


    /// Fit parameters
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);

    FitParameter  tau("tau",2,1.509,0.1);
    FitParameter  dGamma("dGamma",2,0.09,0.1);
    FitParameter  dm("dm",2,17.757,0.1);
    
    FitParameter  scale_mean_dt("scale_mean_dt",1,1,0.1);
    FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
    FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.2,0.1);
    FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
    FitParameter  p0_os("p0_os",1,0.,0.);
    FitParameter  p1_os("p1_os",1,0.,0.);
    FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
    FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
    FitParameter  avg_eta_os("avg_eta_os",1,0.,0.);
    FitParameter  tageff_os("tageff_os",1,0.,0.);
    FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
    FitParameter  p0_ss("p0_ss",1,0.,0.);
    FitParameter  p1_ss("p1_ss",1,0.,0.);
    FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
    FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
    FitParameter  avg_eta_ss("avg_eta_ss",1,0.,0.);
    FitParameter  tageff_ss("tageff_ss",1,0.,0.);
    FitParameter  tageff_asym_ss("tageff_asym_ss",1,0.,0.);
    FitParameter  production_asym("production_asym",1,0.,0.);
    FitParameter  detection_asym("detection_asym",1,0.,0.);
    
    FitParameter  c0("c0",1,1,0.1);
    FitParameter  c1("c1",1,1,0.1);
    FitParameter  c2("c2",1,1,0.1);
    FitParameter  c3("c3",1,1,0.1);
    FitParameter  c4("c4",1,1,0.1);
    FitParameter  c5("c5",1,1,0.1);
    FitParameter  c6("c6",1,1,0.1);
    FitParameter  c7("c7",1,1,0.1);
    FitParameter  c8("c8",1,1,0.1);
    FitParameter  c9("c9",1,1,0.1);
    
    /// Fit parameters per Run    
    FitParameter  scale_mean_dt_Run1("scale_mean_dt_Run1",1,1,0.1);
    FitParameter  offset_sigma_dt_Run1("offset_sigma_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run1("scale_sigma_dt_Run1",1,1.2,0.1);
    FitParameter  scale_sigma_2_dt_Run1("scale_sigma_2_dt_Run1",1,0.,0.1);
    FitParameter  p0_os_Run1("p0_os_Run1",1,0.,0.);
    FitParameter  p1_os_Run1("p1_os_Run1",1,0.,0.);
    FitParameter  delta_p0_os_Run1("delta_p0_os_Run1",1,0.,0.);
    FitParameter  delta_p1_os_Run1("delta_p1_os_Run1",1,0.,0.);
    FitParameter  avg_eta_os_Run1("avg_eta_os_Run1",1,0.,0.);
    FitParameter  tageff_os_Run1("tageff_os_Run1",1,0.,0.);
    FitParameter  tageff_asym_os_Run1("tageff_asym_os_Run1",1,0.,0.);
    FitParameter  p0_ss_Run1("p0_ss_Run1",1,0.,0.);
    FitParameter  p1_ss_Run1("p1_ss_Run1",1,0.,0.);
    FitParameter  delta_p0_ss_Run1("delta_p0_ss_Run1",1,0.,0.);
    FitParameter  delta_p1_ss_Run1("delta_p1_ss_Run1",1,0.,0.);
    FitParameter  avg_eta_ss_Run1("avg_eta_ss_Run1",1,0.,0.);
    FitParameter  tageff_ss_Run1("tageff_ss_Run1",1,0.,0.);
    FitParameter  tageff_asym_ss_Run1("tageff_asym_ss_Run1",1,0.,0.);
    FitParameter  production_asym_Run1("production_asym_Run1",1,0.,0.);
    FitParameter  detection_asym_Run1("detection_asym_Run1",1,0.,0.);
    
    FitParameter  scale_mean_dt_Run2("scale_mean_dt_Run2",1,1,0.1);
    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.2,0.1);
    FitParameter  scale_sigma_2_dt_Run2("scale_sigma_2_dt_Run2",1,1.2,0.1);
    FitParameter  p0_os_Run2("p0_os_Run2",1,0.,0.);
    FitParameter  p1_os_Run2("p1_os_Run2",1,0.,0.);
    FitParameter  delta_p0_os_Run2("delta_p0_os_Run2",1,0.,0.);
    FitParameter  delta_p1_os_Run2("delta_p1_os_Run2",1,0.,0.);
    FitParameter  avg_eta_os_Run2("avg_eta_os_Run2",1,0.,0.);
    FitParameter  tageff_os_Run2("tageff_os_Run2",1,0.,0.);
    FitParameter  tageff_asym_os_Run2("tageff_asym_os_Run2",1,0.,0.);
    FitParameter  p0_ss_Run2("p0_ss_Run2",1,0.,0.);
    FitParameter  p1_ss_Run2("p1_ss_Run2",1,0.,0.);
    FitParameter  delta_p0_ss_Run2("delta_p0_ss_Run2",1,0.,0.);
    FitParameter  delta_p1_ss_Run2("delta_p1_ss_Run2",1,0.,0.);
    FitParameter  avg_eta_ss_Run2("avg_eta_ss_Run2",1,0.,0.);
    FitParameter  tageff_ss_Run2("tageff_ss_Run2",1,0.,0.);
    FitParameter  tageff_asym_ss_Run2("tageff_asym_ss_Run2",1,0.,0.);
    FitParameter  production_asym_Run2("production_asym_Run2",1,0.,0.);
    FitParameter  detection_asym_Run2("detection_asym_Run2",1,0.,0.);

    /// Fit parameters per run and trigger cat
    FitParameter  c0_Run1_t0("c0_Run1_t0",1,1,0.1);
    FitParameter  c1_Run1_t0("c1_Run1_t0",1,1,0.1);
    FitParameter  c2_Run1_t0("c2_Run1_t0",1,1,0.1);
    FitParameter  c3_Run1_t0("c3_Run1_t0",1,1,0.1);
    FitParameter  c4_Run1_t0("c4_Run1_t0",1,1,0.1);
    FitParameter  c5_Run1_t0("c5_Run1_t0",1,1,0.1);
    FitParameter  c6_Run1_t0("c6_Run1_t0",1,1,0.1);
    FitParameter  c7_Run1_t0("c7_Run1_t0",1,1,0.1);
    FitParameter  c8_Run1_t0("c8_Run1_t0",1,1,0.1);
    FitParameter  c9_Run1_t0("c9_Run1_t0",1,1,0.1);
    
    FitParameter  c0_Run1_t1("c0_Run1_t1",1,1,0.1);
    FitParameter  c1_Run1_t1("c1_Run1_t1",1,1,0.1);
    FitParameter  c2_Run1_t1("c2_Run1_t1",1,1,0.1);
    FitParameter  c3_Run1_t1("c3_Run1_t1",1,1,0.1);
    FitParameter  c4_Run1_t1("c4_Run1_t1",1,1,0.1);
    FitParameter  c5_Run1_t1("c5_Run1_t1",1,1,0.1);
    FitParameter  c6_Run1_t1("c6_Run1_t1",1,1,0.1);
    FitParameter  c7_Run1_t1("c7_Run1_t1",1,1,0.1);
    FitParameter  c8_Run1_t1("c8_Run1_t1",1,1,0.1);
    FitParameter  c9_Run1_t1("c9_Run1_t1",1,1,0.1);
    
    FitParameter  c0_Run2_t0("c0_Run2_t0",1,1,0.1);
    FitParameter  c1_Run2_t0("c1_Run2_t0",1,1,0.1);
    FitParameter  c2_Run2_t0("c2_Run2_t0",1,1,0.1);
    FitParameter  c3_Run2_t0("c3_Run2_t0",1,1,0.1);
    FitParameter  c4_Run2_t0("c4_Run2_t0",1,1,0.1);
    FitParameter  c5_Run2_t0("c5_Run2_t0",1,1,0.1);
    FitParameter  c6_Run2_t0("c6_Run2_t0",1,1,0.1);
    FitParameter  c7_Run2_t0("c7_Run2_t0",1,1,0.1);
    FitParameter  c8_Run2_t0("c8_Run2_t0",1,1,0.1);
    FitParameter  c9_Run2_t0("c9_Run2_t0",1,1,0.1);
    
    FitParameter  c0_Run2_t1("c0_Run2_t1",1,1,0.1);
    FitParameter  c1_Run2_t1("c1_Run2_t1",1,1,0.1);
    FitParameter  c2_Run2_t1("c2_Run2_t1",1,1,0.1);
    FitParameter  c3_Run2_t1("c3_Run2_t1",1,1,0.1);
    FitParameter  c4_Run2_t1("c4_Run2_t1",1,1,0.1);
    FitParameter  c5_Run2_t1("c5_Run2_t1",1,1,0.1);
    FitParameter  c6_Run2_t1("c6_Run2_t1",1,1,0.1);
    FitParameter  c7_Run2_t1("c7_Run2_t1",1,1,0.1);
    FitParameter  c8_Run2_t1("c8_Run2_t1",1,1,0.1);
    FitParameter  c9_Run2_t1("c9_Run2_t1",1,1,0.1);

    /// Calculate initial coherence factor
    vector<double> k_gen = coherenceFactor(fas,fas_bar,(double)r, (double)delta,eventListPhsp);
    //vector<double> k_gen_CP = coherenceFactor_CP(fas_CP,fas_bar_CP,(double)r, (double)delta,eventListPhsp_CP);
    vector<double> k_gen_CP = coherenceFactor_CP(fas_CP,fas_bar_CP,(double)r, (double)delta,eventListPhsp);
    
    /// Make full time-dependent PDF
    string marginalPdfsPrefix = "comb";
    if(fitGenMC)marginalPdfsPrefix = "Uniform";
    FullAmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP, r,delta,gamma,tau, dGamma, dm
                      ,offset_sigma_dt, scale_mean_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, marginalPdfsPrefix );

    // Simultaneous pdfs
    FullAmpsPdfFlexiFastCPV pdf_Run1_t0(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP, r,delta,gamma,
		      tau, dGamma, dm,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,c0_Run1_t0, c1_Run1_t0, c2_Run1_t0 ,c3_Run1_t0, c4_Run1_t0, c5_Run1_t0
                      ,c6_Run1_t0, c7_Run1_t0, c8_Run1_t0, c9_Run1_t0,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1_t0" );

    FullAmpsPdfFlexiFastCPV pdf_Run1_t1(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP, r,delta,gamma,
		      tau, dGamma, dm,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,c0_Run1_t1, c1_Run1_t1, c2_Run1_t1 ,c3_Run1_t1, c4_Run1_t1, c5_Run1_t1
                      ,c6_Run1_t1, c7_Run1_t1, c8_Run1_t1, c9_Run1_t1,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1_t1" );

    FullAmpsPdfFlexiFastCPV pdf_Run2_t0(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP, r,delta,gamma,
		      tau, dGamma, dm,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
                      ,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
                      ,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
                      p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                      avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                      p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                      avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                      production_asym_Run2, detection_asym_Run2, "Run2_t0" );

    FullAmpsPdfFlexiFastCPV pdf_Run2_t1(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP, r,delta,gamma,
		      tau, dGamma, dm,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
                      ,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
                      ,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
                      p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                      avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                      p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                      avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                      production_asym_Run2, detection_asym_Run2, "Run2_t1" );

    /*
    cout << fas.getVal(eventListPhsp[0]) << endl;
    cout << fas_bar.getVal(eventListPhsp[0]) << endl;
    DalitzEvent evt_test(eventListPhsp[0]);
    evt_test.CP_conjugateYourself();
    cout << fas_CP.getVal(evt_test) << endl;
    cout << fas_bar_CP.getVal(evt_test) << endl;
    throw "";
    */
    
    if(initCPcoeff){
	CPcoeffLL coeffLL(&pdf,0.741313,-0.741313,0.221461,0.53761,-0.333118,0.00910515,0.01); 
	Minimiser miniCP(&coeffLL);
	miniCP.doFit();
	miniCP.printResultVsInput();
	pdf.coherenceFactor();
	throw "";
    }

    /// Constrain total amp
    fas_tot.print();
    FromFileGenerator* fileGen = new FromFileGenerator(IntegratorEventFile, 0, "UPDATE");
    FlexiFastAmplitudeIntegrator integ_tot(pat, &fas_tot, fileGen, gRandom, integPrecision);
    FracLL fracLL(&integ_tot);
    //Minimiser mini2(&fracLL);
    //mini2.doFit();
    //mini2.printResultVsInput();
    //integ_tot.doFinalStats(&mini2);
    //throw "";

    /// Make marginal pdfs for MC generation and plotting
    // values only used for importance sampling:
    double eff_tag_OS = 0.3852;
    double eff_tag_SS = 0.6903;
    if(fitGenMC){
	eff_tag_OS = 1;
	eff_tag_SS = 1;
    }
    TFile* f_pdfs = new TFile("Mistag_pdfs.root","OPEN");
    
    RooRealVar* r_t = new RooRealVar("t", "t",min_TAU, max_TAU);
    RooRealVar* r_dt = new RooRealVar("dt", "dt",(double)min_TAUERR,(double) max_TAUERR);
    RooRealVar* r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5);
    RooRealVar* r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5);
    
    // time
    TH1D* h_t_norm = new TH1D( *((TH1D*) f_pdfs->Get("h_t_norm_comb")));
    RooDataHist* r_h_t = new RooDataHist("r_h_t","r_h_t",*r_t,h_t_norm);
    RooRealVar mean1_t("mean1_t","mean1_t", 1,0.,3.0);
    RooRealVar mean2_t("mean2_t","mean2_t", 0.,0.,10.0);
    RooRealVar sigma1_t("sigma1_t","sigma1_t", 1,0.,10.0);
    RooRealVar sigma2_t("sigma2_t","sigma2_t", 2,0.,20.0);
    RooRealVar f_t("f_t", "f_t", 0.5, 0., 1.);
    
    RooGaussian* gen1_t = new RooGaussian("gen1_t","gen1_t", *r_t, mean1_t,sigma1_t);
    RooGaussian* gen2_t = new RooGaussian("gen2_t","gen2_t", *r_t, mean2_t,sigma2_t);
    RooAddPdf* gen_t=new RooAddPdf("gen_dt", "gen_dt", RooArgList(*gen1_t, *gen2_t), RooArgList(f_t));
    gen_t->fitTo(*r_h_t,Save(kTRUE),SumW2Error(kTRUE));
    
    // time error
    TH1D* h_dt_norm = new TH1D( *((TH1D*) f_pdfs->Get("h_dt_norm_comb")));
    RooDataHist* r_h_dt = new RooDataHist("r_h_dt","r_h_dt",*r_dt,h_dt_norm);
    RooRealVar mean1_dt("mean1_dt","mean1_dt", 0.02,0.,1.0);
    RooRealVar mean2_dt("mean2_dt","mean2_dt", 0.03,0.,1.0);
    RooRealVar sigma1_dt("sigma1_dt","sigma1_dt", 0.015,0.,1.0);
    RooRealVar sigma2_dt("sigma2_dt","sigma2_dt", 0.015,0.,1.0);
    RooRealVar f_dt("f_dt", "f_dt", 0.5, 0., 1.);
    
    RooGaussian* gen1_dt = new RooGaussian("gen1_dt","gen1_dt", *r_dt, mean1_dt,sigma1_dt);
    RooGaussian* gen2_dt = new RooGaussian("gen2_dt","gen2_dt", *r_dt, mean2_dt,sigma2_dt);
    RooAddPdf* gen_dt=new RooAddPdf("gen_dt", "gen_dt", RooArgList(*gen1_dt, *gen2_dt), RooArgList(f_dt));
    gen_dt->fitTo(*r_h_dt,Save(kTRUE),SumW2Error(kTRUE));
    
    // mistag OS
    TH1D* h_w_OS_norm = new TH1D( *((TH1D*) f_pdfs->Get("h_w_OS_norm_comb")));
    RooDataHist* r_h_eta_OS = new RooDataHist("r_eta_OS","r_eta_OS",*r_eta_OS,h_w_OS_norm);
    RooRealVar mean1_eta_OS("mean1_eta_OS","mean1_eta_OS", 0.02,0.,1.0);
    RooRealVar mean2_eta_OS("mean2_eta_OS","mean2_eta_OS", 0.04,0.,1.0);
    RooRealVar sigma1_eta_OS("sigma1_eta_OS","sigma1_eta_OS", 0.015,0.,1.0);
    RooRealVar sigma2_eta_OS("sigma2_eta_OS","sigma2_eta_OS", 0.015,0.,1.0);
    RooRealVar f_eta_OS("f_eta_OS", "f_eta_OS", 0.5, 0., 1.);
    
    RooGaussian* gen1_eta_OS = new RooGaussian("gen1_eta_OS","gen1_eta_OS", *r_eta_OS, mean1_eta_OS,sigma1_eta_OS);
    RooGaussian* gen2_eta_OS = new RooGaussian("gen2_eta_OS","gen2_eta_OS", *r_eta_OS, mean2_eta_OS,sigma2_eta_OS);
    RooAddPdf* gen_eta_OS=new RooAddPdf("gen_eta_OS", "gen_eta_OS", RooArgList(*gen1_eta_OS, *gen2_eta_OS), RooArgList(f_eta_OS));
    gen_eta_OS->fitTo(*r_h_eta_OS,Save(kTRUE),SumW2Error(kTRUE));
    
    // mistag SS
    TH1D* h_w_SS_norm = new TH1D( *((TH1D*) f_pdfs->Get("h_w_SS_norm_comb")));
    RooDataHist* r_h_eta_SS = new RooDataHist("r_eta_SS","r_eta_SS",*r_eta_SS,h_w_SS_norm);
    RooRealVar mean1_eta_SS("mean1_eta_SS","mean1_eta_SS", 0.02,0.,1.0);
    RooRealVar mean2_eta_SS("mean2_eta_SS","mean2_eta_SS", 0.06,0.,1.0);
    RooRealVar sigma1_eta_SS("sigma1_eta_SS","sigma1_eta_SS", 0.015,0.,1.0);
    RooRealVar sigma2_eta_SS("sigma2_eta_SS","sigma2_eta_SS", 0.015,0.,1.0);
    RooRealVar f_eta_SS("f_eta_SS", "f_eta_SS", 0.5, 0., 1.);
    
    RooGaussian* gen1_eta_SS = new RooGaussian("gen1_eta_SS","gen1_eta_SS", *r_eta_SS, mean1_eta_SS,sigma1_eta_SS);
    RooGaussian* gen2_eta_SS = new RooGaussian("gen2_eta_SS","gen2_eta_SS", *r_eta_SS, mean2_eta_SS,sigma2_eta_SS);
    RooAddPdf* gen_eta_SS=new RooAddPdf("gen_eta_SS", "gen_eta_SS", RooArgList(*gen1_eta_SS, *gen2_eta_SS), RooArgList(f_eta_SS));
    gen_eta_SS->fitTo(*r_h_eta_SS,Save(kTRUE),SumW2Error(kTRUE));
    
    // plot
    TCanvas* c= new TCanvas("");
    RooPlot* frame_t= r_t->frame();
    r_h_t->plotOn(frame_t,Name("data"),MarkerSize(1));
    gen_t->plotOn(frame_t,Name("pdf"),LineColor(kBlue),LineWidth(3));
    frame_t->Draw();
    c->Print("samplingPdfs/t.eps");
    
    RooPlot* frame_dt= r_dt->frame();
    r_h_dt->plotOn(frame_dt,Name("data"),MarkerSize(1));
    gen_dt->plotOn(frame_dt,Name("pdf"),LineColor(kBlue),LineWidth(3));
    frame_dt->Draw();
    c->Print("samplingPdfs/dt.eps");
    
    RooPlot* frame_eta_OS= r_eta_OS->frame();
    r_h_eta_OS->plotOn(frame_eta_OS,MarkerSize(1));
    gen_eta_OS->plotOn(frame_eta_OS,LineColor(kBlue),LineWidth(3));
    frame_eta_OS->Draw();
    c->Print("samplingPdfs/eta_OS.eps");
    
    RooPlot* frame_eta_SS= r_eta_SS->frame();
    r_h_eta_SS->plotOn(frame_eta_SS,MarkerSize(1));
    gen_eta_SS->plotOn(frame_eta_SS,LineColor(kBlue),LineWidth(3));
    frame_eta_SS->Draw();
    c->Print("samplingPdfs/eta_SS.eps");
    
    f_pdfs->Close();
    
    /// Read data or generate toys
    DalitzEventList eventList, eventList_f, eventList_f_bar;
    DalitzEventList eventList_Run1_t0,eventList_Run1_t1,eventList_Run2_t0,eventList_Run2_t1;

    // Generate toys
    if(generateNew){
        time_t startTime = time(0);
        
        // simple amplitude model for importance sampling
        FitAmpIncoherentSum fasGen((DalitzEventPattern)pat);
        fasGen.getVal(eventListPhsp[0]);
        fasGen.normalizeAmps(eventListPhsp);
        SignalGenerator sg(pat,&fasGen);
        fasGen.print();
        
        //simple hit and miss
        for(int i = 0; i < Nevents; i++){
            while(true){
                
                double t_MC = 0.;
                double gaus_t = ranLux.Uniform();
                if(gaus_t < f_t.getVal()) {
                    while(true){
                        t_MC = ranLux.Gaus(mean1_t.getVal(),sigma1_t.getVal());
                        if(t_MC < max_TAU && t_MC > min_TAU)break;
                    }
                }
                else {
                    while(true){
                        t_MC = ranLux.Gaus(mean2_t.getVal(),sigma2_t.getVal());
                        if(t_MC < max_TAU && t_MC > min_TAU)break;
                    }
                }
                r_t->setVal(t_MC) ;
                
                double dt_MC = 0.;
                double gaus_f = ranLux.Uniform();
                if(gaus_f < f_dt.getVal()) {
                    while(true){
                        dt_MC = ranLux.Gaus(mean1_dt.getVal(),sigma1_dt.getVal());
                        if(dt_MC < max_TAUERR && dt_MC > min_TAUERR)break;
                    }
                }
                else {
                    while(true){
                        dt_MC = ranLux.Gaus(mean2_dt.getVal(),sigma2_dt.getVal());
                        if(dt_MC < max_TAUERR && dt_MC > min_TAUERR)break;
                    }
                }
                r_dt->setVal(dt_MC) ;
                
                
                double q_rand = ranLux.Uniform();
                int q_OS_MC = 0;
                if (q_rand < eff_tag_OS/2.  ) q_OS_MC = -1;
                if (q_rand > (1.-eff_tag_OS/2.) ) q_OS_MC = 1;
                
                double eta_OS_MC = 0;
                if(q_OS_MC == 0)eta_OS_MC = 0.5;
                else {
                    double gaus_eta_OS = ranLux.Uniform();
                    if(gaus_eta_OS < f_eta_OS.getVal()){
                        while(true){
                            eta_OS_MC = ranLux.Gaus(mean1_eta_OS.getVal(),sigma1_eta_OS.getVal());
                            if(eta_OS_MC > 0. && eta_OS_MC < 0.5) break;
                        }
                    }
                    else{
                        while(true){
                            eta_OS_MC = ranLux.Gaus(mean2_eta_OS.getVal(),sigma2_eta_OS.getVal());
                            if(eta_OS_MC > 0. && eta_OS_MC < 0.5) break;
                        }
                    }
                }
                r_eta_OS->setVal(eta_OS_MC) ;
                
                q_rand = ranLux.Uniform();
                int q_SS_MC = 0;
                if (q_rand < eff_tag_SS/2.  ) q_SS_MC = -1;
                if (q_rand > (1.-eff_tag_SS/2.) ) q_SS_MC = 1;
                
                double eta_SS_MC = 0;
                if(q_SS_MC == 0)eta_SS_MC = 0.5;
                else {
                    double gaus_eta_SS = ranLux.Uniform();
                    if(gaus_eta_SS < f_eta_SS.getVal()){
                        while(true){
                            eta_SS_MC = ranLux.Gaus(mean1_eta_SS.getVal(),sigma1_eta_SS.getVal());
                            if(eta_SS_MC > 0. && eta_SS_MC < 0.5) break;
                        }
                    }
                    else{
                        while(true){
                            eta_SS_MC = ranLux.Gaus(mean2_eta_SS.getVal(),sigma2_eta_SS.getVal());
                            if(eta_SS_MC > 0. && eta_SS_MC < 0.5) break;
                        }
                    }
                }
                r_eta_SS->setVal(eta_SS_MC) ;
                
                q_rand = ranLux.Uniform();
                int f_MC = 0;
                if (q_rand > .5) f_MC = -1;
                else f_MC = 1;
                
                counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
                DalitzEvent evt(evtPtr.get());
                if(!(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2))continue;
                
                if(f_MC<0){
                    evt.CP_conjugateYourself();
                }
                
                evt.setValueInVector(0, t_MC);
                evt.setValueInVector(1, dt_MC);
                evt.setValueInVector(2, f_MC);
                evt.setValueInVector(3, q_OS_MC);
                evt.setValueInVector(4, eta_OS_MC);
                evt.setValueInVector(5, q_SS_MC);
                evt.setValueInVector(6, eta_SS_MC);
                
                //const double pdfVal = pdf.getValForGeneration(evt);
                const double pdfVal = pdf.un_normalised_noPs(evt);
            
                double maxVal = evt.getGeneratorPdfRelativeToPhaseSpace(); //*exp(-t_MC/tau) / ( tau * ( exp(-min_TAU/tau) - exp(-max_TAU/tau) ) )
                maxVal *= gen_t->getVal(RooArgSet(*r_t))
                * gen_dt->getVal(RooArgSet(*r_dt)) 
                *  (abs(q_OS_MC)/2. * eff_tag_OS + ( 1. - abs(q_OS_MC)) * (1.-eff_tag_OS) )
                *  (abs(q_SS_MC)/2. * eff_tag_SS + ( 1. - abs(q_SS_MC)) * (1.-eff_tag_SS) )    
                *  pdf_max;
                if(q_OS_MC != 0) maxVal *= gen_eta_OS->getVal(RooArgSet(*r_eta_OS));
                if(q_SS_MC != 0) maxVal *= gen_eta_SS->getVal(RooArgSet(*r_eta_SS));
                
                
                const double height = ranLux.Uniform(0,maxVal);
                //Safety check on the maxmimum generated height
                if( pdfVal > maxVal ){
                    std::cout << "ERROR: PDF above determined maximum." << std::endl;
                    std::cout << pdfVal << " > " << maxVal << std::endl;
                    //exit(1);
                    pdf_max = pdf_max * 10.;
                }
                
                //Hit-and-miss
                if( height < pdfVal ){
                    eventList.Add(evt);
                    if(f_MC==1)eventList_f.Add(evt);
                    else eventList_f_bar.Add(evt);
                    break;
                }
            }
            if (0ul == (i % 500ul)) cout << "Generated " << i << " events" << endl;
        }
        cout << " Generated " << Nevents << " Events. Took " << (time(0) - startTime)/60. << " mins "<< endl;
    }
    else {
        
        /// Load data
        double t,dt;
        int f;
        double Bs_ID,Ds_ID;
        int q_OS;
        Short_t q_SS;
        double eta_OS;
        Float_t eta_SS;
        double sw;
        int year,run,Ds_finalState,trigger;

        double K[4];
        double pip[4];
        double pim[4];
        double Ds_Kp[4],Ds_Km[4],Ds_pim[4];
        double mB;
        
        TChain* tree;

	if(fitGenMC){
		tree=new TChain("MCDecayTreeTuple/MCDecayTree");
		tree->Add(((string)InputGenMCFile).c_str());
		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("*TAU*",1);
		tree->SetBranchStatus("*ID*",1);
		tree->SetBranchStatus("*P*",1);
		/*
		tree->SetBranchAddress("B_s0_TRUETAU",&t);
		tree->SetBranchAddress("D_sminus_TRUEID",&Ds_ID);
		tree->SetBranchAddress("B_s0_TRUEID",&Bs_ID);
		tree->SetBranchAddress("Kplus_TRUEP_X",&K[0]);
		tree->SetBranchAddress("Kplus_TRUEP_Y",&K[1]);
		tree->SetBranchAddress("Kplus_TRUEP_Z",&K[2]);
		tree->SetBranchAddress("Kplus_TRUEP_E",&K[3]);
		tree->SetBranchAddress("piplus_TRUEP_X",&pip[0]);
		tree->SetBranchAddress("piplus_TRUEP_Y",&pip[1]);
		tree->SetBranchAddress("piplus_TRUEP_Z",&pip[2]);
		tree->SetBranchAddress("piplus_TRUEP_E",&pip[3]);    
		tree->SetBranchAddress("piminus_TRUEP_X",&pim[0]);
		tree->SetBranchAddress("piminus_TRUEP_Y",&pim[1]);
		tree->SetBranchAddress("piminus_TRUEP_Z",&pim[2]);
		tree->SetBranchAddress("piminus_TRUEP_E",&pim[3]);    
		tree->SetBranchAddress("Kplus0_TRUEP_X",&Ds_Kp[0]);
		tree->SetBranchAddress("Kplus0_TRUEP_Y",&Ds_Kp[1]);
		tree->SetBranchAddress("Kplus0_TRUEP_Z",&Ds_Kp[2]);
		tree->SetBranchAddress("Kplus0_TRUEP_E",&Ds_Kp[3]);
		tree->SetBranchAddress("Kminus_TRUEP_X",&Ds_Km[0]);
		tree->SetBranchAddress("Kminus_TRUEP_Y",&Ds_Km[1]);
		tree->SetBranchAddress("Kminus_TRUEP_Z",&Ds_Km[2]);
		tree->SetBranchAddress("Kminus_TRUEP_E",&Ds_Km[3]);
		tree->SetBranchAddress("piminus0_TRUEP_X",&Ds_pim[0]);
		tree->SetBranchAddress("piminus0_TRUEP_Y",&Ds_pim[1]);
		tree->SetBranchAddress("piminus0_TRUEP_Z",&Ds_pim[2]);
		tree->SetBranchAddress("piminus0_TRUEP_E",&Ds_pim[3]);
		*/
		tree->SetBranchAddress("B_s0_TRUETAU",&t);
		tree->SetBranchAddress("D_splus_TRUEID",&Ds_ID);
		tree->SetBranchAddress("B_s0_TRUEID",&Bs_ID);
		tree->SetBranchAddress("Kminus_TRUEP_X",&K[0]);
		tree->SetBranchAddress("Kminus_TRUEP_Y",&K[1]);
		tree->SetBranchAddress("Kminus_TRUEP_Z",&K[2]);
		tree->SetBranchAddress("Kminus_TRUEP_E",&K[3]);
		tree->SetBranchAddress("piplus_TRUEP_X",&pim[0]);
		tree->SetBranchAddress("piplus_TRUEP_Y",&pim[1]);
		tree->SetBranchAddress("piplus_TRUEP_Z",&pim[2]);
		tree->SetBranchAddress("piplus_TRUEP_E",&pim[3]);    
		tree->SetBranchAddress("piminus_TRUEP_X",&pip[0]);
		tree->SetBranchAddress("piminus_TRUEP_Y",&pip[1]);
		tree->SetBranchAddress("piminus_TRUEP_Z",&pip[2]);
		tree->SetBranchAddress("piminus_TRUEP_E",&pip[3]);    
		tree->SetBranchAddress("Kplus_TRUEP_X",&Ds_Kp[0]);
		tree->SetBranchAddress("Kplus_TRUEP_Y",&Ds_Kp[1]);
		tree->SetBranchAddress("Kplus_TRUEP_Z",&Ds_Kp[2]);
		tree->SetBranchAddress("Kplus_TRUEP_E",&Ds_Kp[3]);
		tree->SetBranchAddress("Kminus0_TRUEP_X",&Ds_Km[0]);
		tree->SetBranchAddress("Kminus0_TRUEP_Y",&Ds_Km[1]);
		tree->SetBranchAddress("Kminus0_TRUEP_Z",&Ds_Km[2]);
		tree->SetBranchAddress("Kminus0_TRUEP_E",&Ds_Km[3]);
		tree->SetBranchAddress("piplus0_TRUEP_X",&Ds_pim[0]);
		tree->SetBranchAddress("piplus0_TRUEP_Y",&Ds_pim[1]);
		tree->SetBranchAddress("piplus0_TRUEP_Z",&Ds_pim[2]);
		tree->SetBranchAddress("piplus0_TRUEP_E",&Ds_pim[3]);
    	}
    	else {
		tree=new TChain("DecayTree");
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
		tree->SetBranchStatus("TriggerCat",1);
		tree->SetBranchStatus("run",1);
	
		tree->SetBranchAddress("Bs_DTF_TAU",&t);
		tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
		tree->SetBranchAddress("Ds_ID",&f);
		tree->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
		tree->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&eta_OS);
		tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
		tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&eta_SS);
		tree->SetBranchAddress("N_Bs_sw",&sw);
		tree->SetBranchAddress("year",&year);
		tree->SetBranchAddress("run",&run);
		tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
		tree->SetBranchAddress("TriggerCat",&trigger);
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
	}

	int N_sample = tree->GetEntries();
        vector<int> b_indices;
        while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(ranLux.Uniform(0,N_sample)));
        sort(b_indices.begin(), b_indices.end());
        if(doBootstrap)N_sample = b_indices.size();

	TRandom3 rndm;
	int badEvents = 0;
	for(int i=0; i< N_sample; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << N_sample << endl;
		if(doBootstrap) tree->GetEntry(b_indices[i]);
		else tree->GetEntry(i);
	
		if(fitGenMC){
			if(Ds_ID<0)f=-1;
	        	else if(Ds_ID > 0)f= 1;
		}
		double sign = 1.;
		//if(f > 0) sign = -1.;
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
	
		DalitzEvent evt;
	
		if(f < 0)evt = DalitzEvent(pat, vectorOfvectors);
		else evt = DalitzEvent(pat_CP, vectorOfvectors);
	
		if(!(evt.phaseSpace() > 0.)){
			/*
			cout << "evt " << i << " 0 phsp " << endl; // << evt << endl;
			cout << "t =" << t << endl;
			cout << "dt =" <<dt << endl;
			cout << "year =" <<year <<  endl;
			cout << "sw =" <<sw <<  endl;
			cout << "Dspi =" <<sqrt(evt.s(1,3))/GeV << endl;
			cout << "Kpi =" <<sqrt(evt.s(2,4))/GeV << endl;
			cout << "pipi =" <<sqrt(evt.s(3,4))/GeV << endl;
			cout << "DsKpi =" <<sqrt(evt.sij(s124))/GeV << endl;
			cout << "Kpipi =" <<sqrt(evt.sij(s234))/GeV << endl << endl;
			*/
			badEvents++;
			continue;
		}
		//if(TMath::IsNaN(norm(fas.getVal(evt)))){
			//cout << "evt " << i << " isNaN " << endl << evt << endl;
			//badEvents++;
			//continue;
		//}
	
		if(fitGenMC){		
			t = t*1000.+ranLux.Gaus(0.,0.04);
			if(t < min_TAU || t > max_TAU )continue;
			if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) continue;

			evt.setValueInVector(0, t);
			evt.setValueInVector(1, 0.04);
			evt.setValueInVector(2, -f);
			int q = 0;
			if(Bs_ID>0)q=1;
			else q = -1;
			evt.setValueInVector(3, q);
			evt.setValueInVector(4, 0.);
			evt.setValueInVector(5, q);
			evt.setValueInVector(6, 0.);
			evt.setValueInVector(7, 1);
			evt.setValueInVector(8, 0);
			eventList.Add(evt);
			continue;
		}
	
		if(t < min_TAU || t > max_TAU )continue;
		if( dt < min_TAUERR || dt > max_TAUERR )continue;
	
		evt.setWeight(sw);
		evt.setValueInVector(0, t);
		evt.setValueInVector(1, dt);
		//evt.setValueInVector(2, 1); /// ???
		if(f<0)evt.setValueInVector(2, 1);   /// ???
		else if(f > 0)evt.setValueInVector(2, -1);  /// ???
		else {
			cout << "ERROR:: Undefined final state";
			throw "ERROR";
		}
		evt.setValueInVector(3, sign*q_OS);
		evt.setValueInVector(4, eta_OS);
		evt.setValueInVector(5, sign*q_SS);
		evt.setValueInVector(6, eta_SS);
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);
	
		eventList.Add(evt);
		if(evt.getValueFromVector(2) == 1)eventList_f.Add(evt);
		else eventList_f_bar.Add(evt);
	
		if(run == 1 && trigger == 0) eventList_Run1_t0.Add(evt);
		else if(run == 1 && trigger == 1) eventList_Run1_t1.Add(evt);
		else if(run == 2 && trigger == 0) eventList_Run2_t0.Add(evt);
		else if(run == 2 && trigger == 1) eventList_Run2_t1.Add(evt);
	}
	cout << endl << "bad events " << badEvents << " ( " << badEvents/(double) N_sample * 100. << " %)" << endl << endl;
    }        

    /// Fit with MINT Pdf
    Neg2LL neg2LL(pdf, eventList);   

    Neg2LL neg2LL_Run1_t0(pdf_Run1_t0, eventList_Run1_t0);    
    Neg2LL neg2LL_Run1_t1(pdf_Run1_t1, eventList_Run1_t1);    
    Neg2LL neg2LL_Run2_t0(pdf_Run2_t0, eventList_Run2_t0);    
    Neg2LL neg2LL_Run2_t1(pdf_Run2_t1, eventList_Run2_t1);    
 
    Neg2LLSum neg2LL_sim;
    if(eventList_Run1_t0.size()>0)neg2LL_sim.add(&neg2LL_Run1_t0);
    if(eventList_Run1_t1.size()>0)neg2LL_sim.add(&neg2LL_Run1_t1);
    if(eventList_Run2_t0.size()>0)neg2LL_sim.add(&neg2LL_Run2_t0);
    if(eventList_Run2_t1.size()>0)neg2LL_sim.add(&neg2LL_Run2_t1);

    double stepSize = 1;
    lambda = lambda + (step-1) * stepSize;
    LASSO_flexi lasso(&ampsSig,lambda);
    LASSO_flexi lasso_bar(&ampsSig_bar,lambda);
    Neg2LLSum neg2LL_lasso(&neg2LL,&lasso,&lasso_bar);
    Neg2LLSum neg2LL_sim_lasso(&neg2LL_sim,&lasso,&lasso_bar);    

    if(constrainFF){
	neg2LL_sim.add(&fracLL);
    }
    
    Minimiser mini;
    if(useLASSO){
    	if(doSimFit)mini.attachFunction(&neg2LL_sim_lasso);
    	else mini.attachFunction(&neg2LL_lasso);    
    }
    else {
    	if(doSimFit)mini.attachFunction(&neg2LL_sim);
    	else mini.attachFunction(&neg2LL);    
    }
    mini.doFit();
    mini.printResultVsInput();

    /// Calculate pulls
    gDirectory->cd();
    TFile* paraFile = new TFile(((string)OutputDir+"pull_"+anythingToString((int)seed)+".root").c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;

    TTree* pull_tree = new TTree("Coherence","Coherence");
    double r_val,delta_val,gamma_val,k_val,n2ll;
    double C_val,D_val,S_val;
    double Cbar_val,Dbar_val,Sbar_val;
    double xp_val,xm_val,yp_val,ym_val;
    double chi2_val,chi2_6D_val;
    TBranch* br_r = pull_tree->Branch( "r", &r_val, "r_val/D" );
    TBranch* br_delta = pull_tree->Branch( "delta", &delta_val, "delta_val/D" );
    TBranch* br_gamma = pull_tree->Branch( "gamma", &gamma_val, "gamma_val/D" );
    TBranch* br_k = pull_tree->Branch( "k", &k_val, "k_val/D" );

    TBranch* br_C = pull_tree->Branch( "C", &C_val, "C_val/D" );
    TBranch* br_Cbar = pull_tree->Branch( "Cbar", &Cbar_val, "Cbar_val/D" );
    TBranch* br_D = pull_tree->Branch( "D", &D_val, "D_val/D" );
    TBranch* br_Dbar = pull_tree->Branch( "Dbar", &Dbar_val, "Dbar_val/D" );
    TBranch* br_S = pull_tree->Branch( "S", &S_val, "S_val/D" );
    TBranch* br_Sbar = pull_tree->Branch( "Sbar", &Sbar_val, "Sbar_val/D" );

    TBranch* br_xp = pull_tree->Branch( "xp", &xp_val, "xp_val/D" );
    TBranch* br_xm = pull_tree->Branch( "xm", &xm_val, "xm_val/D" );
    TBranch* br_yp = pull_tree->Branch( "yp", &yp_val, "yp_val/D" );
    TBranch* br_ym = pull_tree->Branch( "ym", &ym_val, "ym_val/D" );

    TBranch* br_chi2 = pull_tree->Branch( "chi2", &chi2_val, "chi2_val/D" );
    TBranch* br_chi2_6D = pull_tree->Branch( "chi2_6D", &chi2_6D_val, "chi2_6D_val/D" );
    TBranch* br_n2ll = pull_tree->Branch( "n2ll", &n2ll, "n2ll/D" );
    TBranch* br_seed = pull_tree->Branch( "seed", &seed, "seed/I" );
    
    MinuitParameterSet::getDefaultSet()->fillNtp(paraFile, ntp);
    ntp->AutoSave();
    
    if(doSimFit)n2ll = neg2LL_sim.getVal();
    else n2ll = neg2LL.getVal();

    for(unsigned int i=0; i < MinuitParameterSet::getDefaultSet()->size(); i++){
        if(0 == MinuitParameterSet::getDefaultSet()->getParPtr(i)) continue;
        if(A_is_in_B("gamma",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))gamma_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();

        if(A_is_in_B("xp",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))xp_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();
        if(A_is_in_B("xm",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))xm_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();
        if(A_is_in_B("yp",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))yp_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();
        if(A_is_in_B("ym",MinuitParameterSet::getDefaultSet()->getParPtr(i)->name()))ym_val = MinuitParameterSet::getDefaultSet()->getParPtr(i)->mean();
    }
    
    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
    vector<double> k_fit = coherenceFactor(fas,fas_bar,(double)r, (double)delta,eventListMC);
    r_val = k_fit[0];
    k_val = k_fit[1];
    delta_val = k_fit[2];
    C_val = k_fit[3];
    D_val = k_fit[4];
    S_val = k_fit[5];
    vector<double> k_fit_CP = coherenceFactor_CP(fas_CP,fas_bar_CP,(double)r, (double)delta,eventListMC);
    Cbar_val = -k_fit_CP[3];
    Dbar_val = k_fit_CP[4];
    Sbar_val = k_fit_CP[5];
    pdf.coherenceFactor();

    string outTableName = (string)OutputDir+"FitAmpResults_"+anythingToString((int)seed);
    if(updateAnaNote)outTableName = "../../../../../TD-AnaNote/latex/tables/fullFit/"+(string)OutputDir+"fitFractions";

    pdf.doFinalStatsAndSaveForAmp12(&mini,outTableName,((string)OutputDir+"fitFractions_"+anythingToString((int)seed)).c_str());
    if(constrainFF)integ_tot.doFinalStats(&mini);

    /// Create result tables
    if(updateAnaNote){

	ofstream resultsfile;
	resultsfile.open(("../../../../../TD-AnaNote/latex/tables/fullFit/"+(string)OutputDir+"result.tex").c_str(),std::ofstream::trunc);
	resultsfile << "\\begin{table}[h]" << "\n";
	resultsfile << "\\centering" << "\n";
// 	resultsfile << "\\small" << "\n";
	resultsfile << "\\caption{Result of the time-dependent amplitude fit to "; 
	resultsfile << "$B_s \\to D_s K \\pi \\pi$";
	resultsfile << " data.}\n";
	resultsfile << "\\begin{tabular}{c c}" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "Fit parameter & Value \\\\" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << std::fixed << std::setprecision(3) 
	<< "$x_{-}$ & " 
// 	<<   mps->getParPtr("xm")->mean() 
	<< " xx.xx " 
	<< " $\\pm$ " << mps->getParPtr("xm")->err() << "\\\\" << "\n";
	resultsfile << std::fixed << std::setprecision(3) 
	<< "$y_{-}$ & " 
// 	<<   mps->getParPtr("ym")->mean() 
	<< " xx.xx " 
	<< " $\\pm$ " << mps->getParPtr("ym")->err() << "\\\\" << "\n";
	resultsfile << std::fixed << std::setprecision(3) 
	<< "$x_{+}$ & " 
// 	<<   mps->getParPtr("xp")->mean() 
	<< " xx.xx " 
	<< " $\\pm$ " << mps->getParPtr("xp")->err() << "\\\\" << "\n";
	resultsfile << std::fixed << std::setprecision(3) 
	<< "$y_{+}$ & " 
// 	<<   mps->getParPtr("yp")->mean() 
	<< " xx.xx " 
	<< " $\\pm$ " << mps->getParPtr("yp")->err() << "\\\\" << "\n";
	
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\end{tabular}" << "\n";
	resultsfile << "\\label{table:fullFit_" << (string) channel << "}" << "\n";
	resultsfile << "\\end{table}";	
	resultsfile.close();

    } 
  
    /// Data histograms
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
    TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.3,1.6);
    TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,1.6);
    TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30)  ;
    TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);

    TH1D* m_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
    TH1D* m_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
    TH1D* m_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
    TH1D* m_Dspipi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_DsK = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
    TH1D* m_DsKpi = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_Dspi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* m_Dspim = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);

    TH1D* h_cosTheta_Kpi= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
    TH1D* h_cosTheta_Dspi= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    TH1D* h_phi_Kpi_Dspi= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

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
    TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.3,1.6);
    TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,1.6);
    TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsKpi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);

    TH1D* m_Kpipi_fit = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
    TH1D* m_Kpi_fit = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
    TH1D* m_pipi_fit = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
    TH1D* m_Dspipi_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_DsK_fit = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
    TH1D* m_DsKpi_fit = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_Dspi_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* m_Dspim_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* h_cosTheta_Kpi_fit= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1); 
    TH1D* h_cosTheta_Dspi_fit= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    TH1D* h_phi_Kpi_Dspi_fit= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

    TH1D* m_Kpipi_fit_A = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
    TH1D* m_Kpi_fit_A = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
    TH1D* m_pipi_fit_A = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
    TH1D* m_Dspipi_fit_A = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_DsK_fit_A = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
    TH1D* m_DsKpi_fit_A = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_Dspi_fit_A = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* m_Dspim_fit_A = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* h_cosTheta_Kpi_fit_A= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1); 
    TH1D* h_cosTheta_Dspi_fit_A= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    TH1D* h_phi_Kpi_Dspi_fit_A= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

    TH1D* m_Kpipi_fit_Abar = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
    TH1D* m_Kpi_fit_Abar = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
    TH1D* m_pipi_fit_Abar = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
    TH1D* m_Dspipi_fit_Abar = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_DsK_fit_Abar = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
    TH1D* m_DsKpi_fit_Abar = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* m_Dspi_fit_Abar = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* m_Dspim_fit_Abar = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* h_cosTheta_Kpi_fit_Abar= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1); 
    TH1D* h_cosTheta_Dspi_fit_Abar= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    TH1D* h_phi_Kpi_Dspi_fit_Abar= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

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

    double N = 0;
    double N_Run1_t0 = 0;
    double N_Run1_t1 = 0;
    double N_Run2_t0 = 0;
    double N_Run2_t1 = 0;
    double N_OS = 0;
    double N_SS = 0;
    double N_OS_SS = 0;
    double w_OS = 0;
    double w_SS = 0;
    double w_OS_SS = 0;
    double D_OS = 0;
    double D_SS = 0;
    double D_OS_SS = 0;
    double D_comb = 0;
    double N_OS_all = 0;
    double N_SS_all = 0;
    double w_OS_all = 0;
    double w_SS_all = 0;
    double D_OS_all = 0;
    double D_SS_all = 0;
    /// Loop over data
    for (int i=0; i<eventList.size(); i++) {

        N += eventList[i].getWeight();
        h_t->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        h_dt->Fill(eventList[i].getValueFromVector(1),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(3) != 0)h_eta_OS->Fill(eventList[i].getValueFromVector(4),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(5) != 0)h_eta_SS->Fill(eventList[i].getValueFromVector(6),eventList[i].getWeight());

        int f_evt = eventList[i].getValueFromVector(2);
        int q1 = eventList[i].getValueFromVector(3);
        int q2 = eventList[i].getValueFromVector(5);   
        int q_eff = 0;
        double w_eff = 0.5;
        int run_evt = eventList[i].getValueFromVector(7);   
        int trigger_evt = eventList[i].getValueFromVector(8);   
        
        if(run_evt==1 && trigger_evt == 0)N_Run1_t0 += eventList[i].getWeight();
        else if(run_evt==1 && trigger_evt == 1)N_Run1_t1 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 0)N_Run2_t0 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 1)N_Run2_t1 += eventList[i].getWeight();

        std::pair<double, double> calibrated_mistag_os;
        std::pair<double, double> calibrated_mistag_ss;

        if(doSimFit){
            if(run_evt==1){
                calibrated_mistag_os = pdf_Run1_t0.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = pdf_Run1_t0.getCalibratedMistag_SS(eventList[i]);
            }
            else{
                calibrated_mistag_os = pdf_Run2_t0.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = pdf_Run2_t0.getCalibratedMistag_SS(eventList[i]);                
            }
        }
        else{
            calibrated_mistag_os = pdf.getCalibratedMistag_OS(eventList[i]);
            calibrated_mistag_ss = pdf.getCalibratedMistag_SS(eventList[i]);        
        }

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

        if(q1 != 0 && q2 != 0){
            if(q_eff != 0){
                N_OS_SS += eventList[i].getWeight();
                w_OS_SS += w_eff * eventList[i].getWeight();
                D_OS_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
            }
        }
        else if( q1 != 0){
                //q_eff = q1;
                N_OS += eventList[i].getWeight();
                w_OS += w_eff * eventList[i].getWeight(); 
                D_OS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
        }
        else if( q2 != 0){
                //q_eff = q2;
                N_SS += eventList[i].getWeight();
	    		w_SS += w_eff * eventList[i].getWeight(); 
                D_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
        } 

        D_comb += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
        
        if(q1 != 0) N_OS_all += eventList[i].getWeight();
        if(q2 != 0)N_SS_all += eventList[i].getWeight();
            
        if(q1>0){
			w_OS_all +=  calibrated_mistag_os.first * eventList[i].getWeight();
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.first,2)* eventList[i].getWeight();
        } 
        else if(q1<0){
			w_OS_all +=  calibrated_mistag_os.second * eventList[i].getWeight();	
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.second,2)* eventList[i].getWeight();
        }

        if(q2>0){
			w_SS_all +=  calibrated_mistag_ss.first * eventList[i].getWeight();
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.first,2)* eventList[i].getWeight();
        } 
        else if(q2<0){
			w_SS_all +=  calibrated_mistag_ss.second * eventList[i].getWeight();	
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.second,2)* eventList[i].getWeight();
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

            m_Kpipi->Fill(sqrt(eventList[i].sij(s234)/(GeV*GeV)),eventList[i].getWeight());
            m_Kpi->Fill(sqrt(eventList[i].s(2,4)/(GeV*GeV)),eventList[i].getWeight());
            m_pipi->Fill(sqrt(eventList[i].s(3,4)/(GeV*GeV)),eventList[i].getWeight());
            m_Dspipi->Fill(sqrt(eventList[i].sij(s134)/(GeV*GeV)),eventList[i].getWeight());
            m_DsK->Fill(sqrt(eventList[i].s(1,2)/(GeV*GeV)),eventList[i].getWeight());
            m_DsKpi->Fill(sqrt(eventList[i].sij(s124)/(GeV*GeV)),eventList[i].getWeight());
            m_Dspi->Fill(sqrt(eventList[i].s(1,3)/(GeV*GeV)),eventList[i].getWeight());
            m_Dspim->Fill(sqrt(eventList[i].s(1,4)/(GeV*GeV)),eventList[i].getWeight());
            h_cosTheta_Kpi->Fill(cosThetaAngle(eventList[i],2,4,1,3),eventList[i].getWeight());
	    h_cosTheta_Dspi->Fill(cosThetaAngle(eventList[i],1,3,2,4),eventList[i].getWeight());
	    h_phi_Kpi_Dspi->Fill(acoplanarityAngle(eventList[i],2,4,1,3),eventList[i].getWeight());

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

    cout << "Tagging perfromance " << endl << endl;        
    cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;

    cout << "OS  | " << (N_OS+N_OS_SS)/N << " | " <<  (w_OS_all)/(N_OS+N_OS_SS) << " | " << D_OS_all/N << endl;
    cout << "SS  | " << (N_SS+N_OS_SS)/N << " | " <<  (w_SS_all)/(N_SS+N_OS_SS) << " | " << D_SS_all/N << endl << endl;

    cout << "OS only  | " << N_OS/N << " | " <<  w_OS/N_OS << " | " << N_OS/N * D_OS/N_OS << endl;
    cout << "SS only  | " << N_SS/N << " | " <<  w_SS/N_SS << " | " << N_SS/N * D_SS/N_SS << endl;
    cout << "OS+SS    | " << N_OS_SS/N << " | " <<  w_OS_SS/N_OS_SS << " | " << N_OS_SS/N * D_OS_SS/N_OS_SS << endl;
    cout << "Combined | " << (N_OS+N_SS+N_OS_SS)/N << " | "<<  (w_OS+w_SS+w_OS_SS)/(N_OS+N_SS+N_OS_SS) << " | " << (N_OS+N_SS+N_OS_SS)/N * D_comb/(N_OS+N_SS+N_OS_SS) << endl << endl ;

    /// Loop over MC
    DiskResidentEventList eventListMC_rw(pat,("dummy_"+anythingToString(step)+".root").c_str(),"RECREATE");
    for(int i = 0; i < eventListMC.size(); i++){
        
	double t_MC,dt_MC,eta_OS_MC,eta_SS_MC;
	int q_OS_MC,q_SS_MC,f_MC,trigger_MC,run_MC;

	if(fitGenMC){
		t_MC = ranLux.Exp(tau);
		if(t_MC > max_TAU && t_MC < min_TAU)continue;
	
		dt_MC = 0.04;
		
		double q_rand = ranLux.Uniform();
		q_OS_MC = 0;
		if (q_rand < eff_tag_OS/2.  ) q_OS_MC = -1;
		if (q_rand > (1.-eff_tag_OS/2.) ) q_OS_MC = 1;
		
		q_rand = ranLux.Uniform();
		q_SS_MC = q_OS_MC;
		
		eta_OS_MC = 0;
		eta_SS_MC = 0;
	
		q_rand = ranLux.Uniform();
		f_MC = 0;
		if (q_rand > .5) f_MC = -1;
		else f_MC = 1;
	}
	else {
	
		t_MC = 0.;
		double gaus_t = ranLux.Uniform();
		if(gaus_t < f_t.getVal()) {
		while(true){
			t_MC = ranLux.Gaus(mean1_t.getVal(),sigma1_t.getVal());
			if(t_MC < max_TAU && t_MC > min_TAU)break;
		}
		}
		else {
		while(true){
			t_MC = ranLux.Gaus(mean2_t.getVal(),sigma2_t.getVal());
			if(t_MC < max_TAU && t_MC > min_TAU)break;
		}
		}
		r_t->setVal(t_MC) ;
		
		dt_MC = 0.;
		double gaus_f = ranLux.Uniform();
		if(gaus_f < f_dt.getVal()) {
		while(true){
			dt_MC = ranLux.Gaus(mean1_dt.getVal(),sigma1_dt.getVal());
			if(dt_MC < max_TAUERR && dt_MC > min_TAUERR)break;
		}
		}
		else {
		while(true){
			dt_MC = ranLux.Gaus(mean2_dt.getVal(),sigma2_dt.getVal());
			if(dt_MC < max_TAUERR && dt_MC > min_TAUERR)break;
		}
		}
		r_dt->setVal(dt_MC) ;
		
		
		double q_rand = ranLux.Uniform();
		q_OS_MC = 0;
		if (q_rand < eff_tag_OS/2.  ) q_OS_MC = -1;
		if (q_rand > (1.-eff_tag_OS/2.) ) q_OS_MC = 1;
		
		eta_OS_MC = 0;
		if(q_OS_MC == 0)eta_OS_MC = 0.5;
		else {
		double gaus_eta_OS = ranLux.Uniform();
		if(gaus_eta_OS < f_eta_OS.getVal()){
			while(true){
			eta_OS_MC = ranLux.Gaus(mean1_eta_OS.getVal(),sigma1_eta_OS.getVal());
			if(eta_OS_MC > 0. && eta_OS_MC < 0.5) break;
			}
		}
		else{
			while(true){
			eta_OS_MC = ranLux.Gaus(mean2_eta_OS.getVal(),sigma2_eta_OS.getVal());
			if(eta_OS_MC > 0. && eta_OS_MC < 0.5) break;
			}
		}
		}
		r_eta_OS->setVal(eta_OS_MC) ;
		
		q_rand = ranLux.Uniform();
		q_SS_MC = 0;
		if (q_rand < eff_tag_SS/2.  ) q_SS_MC = -1;
		if (q_rand > (1.-eff_tag_SS/2.) ) q_SS_MC = 1;
		
		eta_SS_MC = 0;
		if(q_SS_MC == 0)eta_SS_MC = 0.5;
		else {
		double gaus_eta_SS = ranLux.Uniform();
		if(gaus_eta_SS < f_eta_SS.getVal()){
			while(true){
			eta_SS_MC = ranLux.Gaus(mean1_eta_SS.getVal(),sigma1_eta_SS.getVal());
			if(eta_SS_MC > 0. && eta_SS_MC < 0.5) break;
			}
		}
		else{
			while(true){
			eta_SS_MC = ranLux.Gaus(mean2_eta_SS.getVal(),sigma2_eta_SS.getVal());
			if(eta_SS_MC > 0. && eta_SS_MC < 0.5) break;
			}
		}
		}
		r_eta_SS->setVal(eta_SS_MC) ;
		
		q_rand = ranLux.Uniform();
		f_MC = 0;
		if (q_rand > .5) f_MC = -1;
		else f_MC = 1;
	
		q_rand = ranLux.Uniform();
		run_MC = 0;
		if (q_rand > .5) run_MC = 1;
		else run_MC = 2;
	
		q_rand = ranLux.Uniform();
		trigger_MC = 0;
		if (q_rand > .5) trigger_MC = 0;
		else trigger_MC = 1;
	}

        DalitzEvent evt(eventListMC.getEvent(i));
        if(f_MC<0){
            evt.CP_conjugateYourself();
            //evt.P_conjugateYourself();
        }
        //if(!(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2))continue;
        evt.setValueInVector(0, t_MC);
        evt.setValueInVector(1, dt_MC);
        evt.setValueInVector(2, f_MC);
        evt.setValueInVector(3, q_OS_MC);
        evt.setValueInVector(4, eta_OS_MC);
        evt.setValueInVector(5, q_SS_MC);
        evt.setValueInVector(6, eta_SS_MC);
        if(!fitGenMC)evt.setValueInVector(7, run_MC);
        if(!fitGenMC)evt.setValueInVector(8, trigger_MC);
        
        double pdfVal = 0;
        if(doSimFit) {
            if(run_MC==1 && trigger_MC == 0) pdfVal = pdf_Run1_t0.getVal(evt) * N_Run1_t0/N;
            else if(run_MC==1 && trigger_MC == 1) pdfVal = pdf_Run1_t1.getVal(evt) * N_Run1_t1/N;
            else if(run_MC==2 && trigger_MC == 0) pdfVal = pdf_Run2_t0.getVal(evt) * N_Run2_t0/N;
            else if(run_MC==2 && trigger_MC == 1) pdfVal = pdf_Run2_t1.getVal(evt) * N_Run2_t1/N;
        }
        else pdfVal = pdf.getVal(evt);
        //const double pdfVal = pdf.getValForGeneration(evt);
        //const double pdfVal = pdf.un_normalised_noPs(evt);

        //double weight = pdfVal/exp(-fabs(t_MC)/(tau))*tau*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
        double weight = pdfVal*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	if(fitGenMC){
		weight /=  exp(-t_MC/tau) / ( tau * ( exp(-min_TAU/tau) - exp(-max_TAU/tau) ) ) 
		*  (abs(q_OS_MC)/2. * eff_tag_OS + ( 1. - abs(q_OS_MC)) * (1.-eff_tag_OS) ) ;
	}
	else {
		weight /=  //exp(-t_MC/tau) / ( tau * ( exp(-min_TAU/tau) - exp(-max_TAU/tau) ) ) 
		gen_t->getVal(RooArgSet(*r_t))
		* gen_dt->getVal(RooArgSet(*r_dt)) 
		*  (abs(q_OS_MC)/2. * eff_tag_OS + ( 1. - abs(q_OS_MC)) * (1.-eff_tag_OS) )
		*  (abs(q_SS_MC)/2. * eff_tag_SS + ( 1. - abs(q_SS_MC)) * (1.-eff_tag_SS) )    ;
		if(q_OS_MC != 0) weight /= gen_eta_OS->getVal(RooArgSet(*r_eta_OS));
		if(q_SS_MC != 0) weight /= gen_eta_SS->getVal(RooArgSet(*r_eta_SS));
	}

        h_t_fit->Fill(t_MC,weight);
        h_dt_fit->Fill(dt_MC,weight);
        if(evt.getValueFromVector(3) != 0)h_eta_OS_fit->Fill(evt.getValueFromVector(4),weight);
        if(evt.getValueFromVector(5) != 0)h_eta_SS_fit->Fill(evt.getValueFromVector(6),weight);
        
        int f_evt = evt.getValueFromVector(2);
        int q1 = evt.getValueFromVector(3);
        int q2 = evt.getValueFromVector(5);   
        int q_eff = 0;
        double w_eff = 0.5;
        
	std::pair<double, double> calibrated_mistag_os;
	std::pair<double, double> calibrated_mistag_ss;
	if(doSimFit){
		if(run_MC==1){
			calibrated_mistag_os = pdf_Run1_t0.getCalibratedMistag_OS(evt);
			calibrated_mistag_ss = pdf_Run1_t0.getCalibratedMistag_SS(evt);
		}
		else{
			calibrated_mistag_os = pdf_Run2_t0.getCalibratedMistag_OS(evt);
			calibrated_mistag_ss = pdf_Run2_t0.getCalibratedMistag_SS(evt);                
		}
	}
	else{
		calibrated_mistag_os = pdf.getCalibratedMistag_OS(evt);
		calibrated_mistag_ss = pdf.getCalibratedMistag_SS(evt);        
	}
	
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

	m_Kpipi_fit->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight);
	m_Kpi_fit->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight);
	m_pipi_fit->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight);
	m_Dspipi_fit->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight);
	m_DsK_fit->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight);
	m_DsKpi_fit->Fill(sqrt(evt.sij(s124)/(GeV*GeV)),weight);
	m_Dspi_fit->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight);
	m_Dspim_fit->Fill(sqrt(evt.s(1,4)/(GeV*GeV)),weight);
	h_cosTheta_Kpi_fit->Fill(cosThetaAngle(evt,2,4,1,3),weight);
	h_cosTheta_Dspi_fit->Fill(cosThetaAngle(evt,1,3,2,4),weight);
	h_phi_Kpi_Dspi_fit->Fill(acoplanarityAngle(evt,2,4,1,3),weight);

	double weight_A = 0;
	double weight_Abar = 0;
	if(f_MC>0){
		weight_A = fas.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
		weight_Abar = fas_bar.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	else {
		weight_A = fas_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
		weight_Abar = fas_bar_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
	}

	m_Kpipi_fit_A->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight_A);
	m_Kpi_fit_A->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight_A);
	m_pipi_fit_A->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight_A);
	m_Dspipi_fit_A->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight_A);
	m_DsK_fit_A->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight_A);
	m_DsKpi_fit_A->Fill(sqrt(evt.sij(s124)/(GeV*GeV)),weight_A);
	m_Dspi_fit_A->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_A);
	m_Dspim_fit_A->Fill(sqrt(evt.s(1,4)/(GeV*GeV)),weight_A);
	h_cosTheta_Kpi_fit_A->Fill(cosThetaAngle(evt,2,4,1,3),weight_A);
	h_cosTheta_Dspi_fit_A->Fill(cosThetaAngle(evt,1,3,2,4),weight_A);
	h_phi_Kpi_Dspi_fit_A->Fill(acoplanarityAngle(evt,2,4,1,3),weight_A);

	m_Kpipi_fit_Abar->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight_Abar);
	m_Kpi_fit_Abar->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight_Abar);
	m_pipi_fit_Abar->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight_Abar);
	m_Dspipi_fit_Abar->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight_Abar);
	m_DsK_fit_Abar->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight_Abar);
	m_DsKpi_fit_Abar->Fill(sqrt(evt.sij(s124)/(GeV*GeV)),weight_Abar);
	m_Dspi_fit_Abar->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_Abar);
	m_Dspim_fit_Abar->Fill(sqrt(evt.s(1,4)/(GeV*GeV)),weight_Abar);
	h_cosTheta_Kpi_fit_Abar->Fill(cosThetaAngle(evt,2,4,1,3),weight_Abar);
	h_cosTheta_Dspi_fit_Abar->Fill(cosThetaAngle(evt,1,3,2,4),weight_Abar);
	h_phi_Kpi_Dspi_fit_Abar->Fill(acoplanarityAngle(evt,2,4,1,3),weight_Abar);

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

	evt.setWeight(weight);
	eventListMC_rw.Add(evt);
    }

    /// Plot
    h_t->SetMinimum(0.1);    
    h_t->SetLineColor(kBlack);
    h_t->DrawNormalized("e",1);
        
    h_t_fit->SetLineColor(kBlue);
    h_t_fit->SetLineWidth(3);
    h_t_fit->SetMarkerColor(kBlue); 
    h_t_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_t.eps").c_str());
    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t.pdf").c_str());
    gPad->SetLogy(1);
    c->Print(((string)OutputDir+"h_t_log.eps").c_str());
    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
    gPad->SetLogy(0);
    
    h_dt->SetMinimum(0);        
    h_dt->SetLineColor(kBlack);
    h_dt->DrawNormalized("e1",1);
    h_dt_fit->SetLineColor(kBlue);
    h_dt_fit->SetLineWidth(3);
    h_dt_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_dt.eps").c_str());
    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_dt.pdf").c_str());

    h_eta_OS->SetMinimum(0);        
    h_eta_OS->SetLineColor(kBlack);
    h_eta_OS->DrawNormalized("e1",1);
    h_eta_OS_fit->SetLineColor(kBlue);
    h_eta_OS_fit->SetLineWidth(3);
    h_eta_OS_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_eta_OS.eps").c_str());
    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_eta_OS.pdf").c_str());
    
    h_eta_SS->SetMinimum(0);        
    h_eta_SS->SetLineColor(kBlack);
    h_eta_SS->DrawNormalized("e1",1);
    h_eta_SS_fit->SetLineColor(kBlue);
    h_eta_SS_fit->SetLineWidth(3);
    h_eta_SS_fit->DrawNormalized("histcsame",1);
    c->Print(((string)OutputDir+"h_eta_SS.eps").c_str());
    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_eta_SS.pdf").c_str());

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
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed.pdf").c_str());

	TH1D* h_asym = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);	
        h_asym->SetMinimum(-0.25);
	h_asym->SetMaximum(0.25);
	TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
	h_asym_fit->SetLineColor(kRed);
	h_asym->Draw("e");
	h_asym_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
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
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_p.pdf").c_str());

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
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());

	cout << h_N_unmixed_p->GetEntries() << endl;
	cout << h_N_mixed_p->GetEntries() << endl;

	TH1D* h_asym_p = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_mixed_p);	
        //h_asym_p->SetMinimum(-20);
	//h_asym_p->SetMaximum(20);
	TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);	
	h_asym_p_fit->SetLineColor(kRed);
	h_asym_p->Draw("e");
	h_asym_p_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_p.eps").c_str());
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym_p.pdf").c_str());

	TH1D* h_asym_m = (TH1D*) h_N_unmixed_m->GetAsymmetry(h_N_mixed_m);	
        //h_asym_m->SetMinimum(-20);
	//h_asym_m->SetMaximum(20);
	TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);	
	h_asym_m_fit->SetLineColor(kRed);
	h_asym_m->Draw("e");
	h_asym_m_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_m.eps").c_str());
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym_m.pdf").c_str());

	h_asym_p->SetMaximum(max(h_asym_p->GetMaximum(),h_asym_m->GetMaximum())*1.25);
	h_asym_p->SetMarkerColor(kRed);
	h_asym_p->SetLineColor(kRed);
	h_asym_p->Draw("e");
	h_asym_m->SetLineColor(kBlue);
	h_asym_m->SetMarkerColor(kBlue);
	h_asym_m->Draw("esame");
	h_asym_p_fit->SetLineWidth(3);
	h_asym_m_fit->SetLineWidth(3);
	h_asym_p_fit->Draw("histcsame");
	h_asym_m_fit->SetLineColor(kBlue);
 	h_asym_m_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());	
        if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
    }
    
            s_Kpipi->SetMinimum(0);
            s_Kpipi->SetLineColor(kBlack);
            s_Kpipi->DrawNormalized("e1",1);
            s_Kpipi_fit->SetLineColor(kBlue);
            s_Kpipi_fit->SetLineWidth(3);
            s_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_Kpipi.pdf").c_str());

            s_Kpi->SetMinimum(0);
            s_Kpi->SetLineColor(kBlack);
            s_Kpi->DrawNormalized("e1",1);
            s_Kpi_fit->SetLineColor(kBlue);
            s_Kpi_fit->SetLineWidth(3);
            s_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_Kpi.pdf").c_str());

	    s_pipi->SetMinimum(0);            
            s_pipi->SetLineColor(kBlack);
            s_pipi->DrawNormalized("e1",1);
            s_pipi_fit->SetLineColor(kBlue);
            s_pipi_fit->SetLineWidth(3);
            s_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_pipi.pdf").c_str());

	    s_Dspipi->SetMinimum(0);            
            s_Dspipi->SetLineColor(kBlack);
            s_Dspipi->DrawNormalized("e1",1);
            s_Dspipi_fit->SetLineColor(kBlue);
            s_Dspipi_fit->SetLineWidth(3);
            s_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_Dspipi.pdf").c_str());

	    s_DsK->SetMinimum(0);
            s_DsK->SetLineColor(kBlack);
            s_DsK->DrawNormalized("e1",1);
            s_DsK_fit->SetLineColor(kBlue);
            s_DsK_fit->SetLineWidth(3);
            s_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_DsK.pdf").c_str());
	    
	    s_DsKpi->SetMinimum(0);            
            s_DsKpi->SetLineColor(kBlack);
            s_DsKpi->DrawNormalized("e1",1);
            s_DsKpi_fit->SetLineColor(kBlue);
            s_DsKpi_fit->SetLineWidth(3);
            s_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_DsKpi.pdf").c_str());
	    
	    s_Dspi->SetMinimum(0);
            s_Dspi->SetLineColor(kBlack);
            s_Dspi->DrawNormalized("e1",1);
            s_Dspi_fit->SetLineColor(kBlue);
            s_Dspi_fit->SetLineWidth(3);
            s_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_Dspi.pdf").c_str());
	    
	    s_Dspim->SetMinimum(0);
            s_Dspim->SetLineColor(kBlack);
            s_Dspim->DrawNormalized("e1",1);
            s_Dspim_fit->SetLineColor(kBlue);
            s_Dspim_fit->SetLineWidth(3);
            s_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"s_Dspim.pdf").c_str());
	    
            m_Kpipi->SetMinimum(0);
            m_Kpipi->SetLineColor(kBlack);
            m_Kpipi->DrawNormalized("e1",1);
            m_Kpipi_fit->SetLineColor(kBlue);
            m_Kpipi_fit->SetLineWidth(3);
            m_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Kpipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpipi.pdf").c_str());

            m_Kpi->SetMinimum(0);
            m_Kpi->SetLineColor(kBlack);
            m_Kpi->DrawNormalized("e1",1);
            m_Kpi_fit->SetLineColor(kBlue);
            m_Kpi_fit->SetLineWidth(3);
            m_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Kpi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpi.pdf").c_str());

	    m_pipi->SetMinimum(0);            
            m_pipi->SetLineColor(kBlack);
            m_pipi->DrawNormalized("e1",1);
            m_pipi_fit->SetLineColor(kBlue);
            m_pipi_fit->SetLineWidth(3);
            m_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_pipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_pipi.pdf").c_str());

	    m_Dspipi->SetMinimum(0);            
            m_Dspipi->SetLineColor(kBlack);
            m_Dspipi->DrawNormalized("e1",1);
            m_Dspipi_fit->SetLineColor(kBlue);
            m_Dspipi_fit->SetLineWidth(3);
            m_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspipi.pdf").c_str());

	    m_DsK->SetMinimum(0);
            m_DsK->SetLineColor(kBlack);
            m_DsK->DrawNormalized("e1",1);
            m_DsK_fit->SetLineColor(kBlue);
            m_DsK_fit->SetLineWidth(3);
            m_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_DsK.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsK.pdf").c_str());
	    
	    m_DsKpi->SetMinimum(0);            
            m_DsKpi->SetLineColor(kBlack);
            m_DsKpi->DrawNormalized("e1",1);
            m_DsKpi_fit->SetLineColor(kBlue);
            m_DsKpi_fit->SetLineWidth(3);
            m_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_DsKpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsKpi.pdf").c_str());
	    
	    m_Dspi->SetMinimum(0);
            m_Dspi->SetLineColor(kBlack);
            m_Dspi->DrawNormalized("e1",1);
            m_Dspi_fit->SetLineColor(kBlue);
            m_Dspi_fit->SetLineWidth(3);
            m_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspi.pdf").c_str());
	    
	    m_Dspim->SetMinimum(0);
            m_Dspim->SetLineColor(kBlack);
            m_Dspim->DrawNormalized("e1",1);
            m_Dspim_fit->SetLineColor(kBlue);
            m_Dspim_fit->SetLineWidth(3);
            m_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspim.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspim.pdf").c_str());

	    h_cosTheta_Kpi->SetMinimum(0);
            h_cosTheta_Kpi->SetLineColor(kBlack);
            h_cosTheta_Kpi->DrawNormalized("e1",1);
            h_cosTheta_Kpi_fit->SetLineColor(kBlue);
            h_cosTheta_Kpi_fit->SetLineWidth(3);
            h_cosTheta_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_cosTheta_Kpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Kpi.pdf").c_str());

	    h_cosTheta_Dspi->SetMinimum(0);
            h_cosTheta_Dspi->SetLineColor(kBlack);
            h_cosTheta_Dspi->DrawNormalized("e1",1);
            h_cosTheta_Dspi_fit->SetLineColor(kBlue);
            h_cosTheta_Dspi_fit->SetLineWidth(3);
            h_cosTheta_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_cosTheta_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Dspi.pdf").c_str());

	    h_phi_Kpi_Dspi->SetMinimum(0);
            h_phi_Kpi_Dspi->SetLineColor(kBlack);
            h_phi_Kpi_Dspi->DrawNormalized("e1",1);
            h_phi_Kpi_Dspi_fit->SetLineColor(kBlue);
            h_phi_Kpi_Dspi_fit->SetLineWidth(3);
            h_phi_Kpi_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_phi_Kpi_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_phi_Kpi_Dspi.pdf").c_str());
	    

	    double norm_A = m_Kpipi_fit_A->Integral()/(m_Kpipi_fit_A->Integral()+m_Kpipi_fit_Abar->Integral());
	    double norm_Abar = m_Kpipi_fit_Abar->Integral()/(m_Kpipi_fit_Abar->Integral()+m_Kpipi_fit_Abar->Integral());
	    cout << "ration = " << norm_Abar/norm_A << endl;

            m_Kpipi->SetMinimum(0.01);
            m_Kpipi->SetLineColor(kBlack);
            m_Kpipi->DrawNormalized("e1",1);
            m_Kpipi_fit->SetLineColor(kBlue);
            m_Kpipi_fit->SetLineWidth(3);
            m_Kpipi_fit->DrawNormalized("histcsame",1);
            m_Kpipi_fit_A->SetLineColor(kRed+1);
            m_Kpipi_fit_A->SetLineWidth(2);
            m_Kpipi_fit_A->SetFillColor(kRed+1);
            m_Kpipi_fit_A->SetFillStyle(3353);
            m_Kpipi_fit_A->DrawNormalized("histcsame",norm_A);
            m_Kpipi_fit_Abar->SetLineColor(kGreen+3);
            m_Kpipi_fit_Abar->SetLineWidth(2);
            m_Kpipi_fit_Abar->SetFillColor(kGreen+3);
            m_Kpipi_fit_Abar->SetFillStyle(3353);
            m_Kpipi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
            c->Print(((string)OutputDir+"m_Kpipi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpipi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Kpipi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpipi_mod_log.pdf").c_str());
            gPad->SetLogy(0);



	    chi2_val = getChi2(eventList,eventListMC_rw);
	    //chi2_6D_val = getChi2_6D(eventList,eventListMC_rw);

	    if(useLASSO){
		paraFile->cd();

		Double_t x[1],y[1];
		vector<double> thresholds;
		thresholds.push_back(0.001);
		thresholds.push_back(0.002);
		thresholds.push_back(0.003);
		thresholds.push_back(0.004);
		thresholds.push_back(0.005);
		thresholds.push_back(0.006);
		thresholds.push_back(0.007);
		thresholds.push_back(0.008);
		thresholds.push_back(0.009);
		thresholds.push_back(0.01);
		thresholds.push_back(0.02);
		thresholds.push_back(0.05);
		
		x[0]=lambda;
		for(int i = 0; i < thresholds.size() ; i++){
			if(doSimFit)y[0]=neg2LL_sim.getVal() + 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			else y[0]=neg2LL.getVal() + 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			y[0]+= 2. * lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* aic = new TGraph(1,x,y);
			aic->SetName( ("AIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			aic->SetTitle("");
			aic->GetXaxis()->SetTitle("#lambda");
			aic->GetXaxis()->SetTitleOffset(0.65);
			aic->GetYaxis()->SetTitle("AIC");
			aic->Draw("A*");
			aic->Write();
		}
		
		for(int i = 0; i < thresholds.size() ; i++){
			if(doSimFit)y[0]=neg2LL_sim.getVal() + lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(eventList.size());
			else y[0]=neg2LL.getVal() + lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(eventList.size());
			y[0]+= lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(eventList.size());
			TGraph* bic = new TGraph(1,x,y);
			bic->SetName( ("BIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			bic->SetTitle("");
			bic->GetXaxis()->SetTitle("#lambda");
			bic->GetXaxis()->SetTitleOffset(0.65);
			bic->GetYaxis()->SetTitle("BIC");
			bic->Draw("A*");
			bic->Write();
		}
		
		for(int i = 0; i < thresholds.size() ; i++){
			y[0]=lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			y[0]+= lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		for(int i = 0; i < thresholds.size() ; i++){
			y[0]=lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_A_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		for(int i = 0; i < thresholds.size() ; i++){
			y[0]=lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_Abar_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		if(doSimFit)y[0]=neg2LL_sim.getVal() ;
		else y[0]=neg2LL.getVal() ;
		TGraph* nll = new TGraph(1,x,y);
		nll->SetName("Neg2LL");
		nll->SetTitle("");
		nll->GetXaxis()->SetTitle("#lambda");
		nll->GetXaxis()->SetTitleOffset(0.65);
		nll->GetYaxis()->SetTitle("-2 logL");
		nll->Draw("A*");
		nll->Write();
		
		y[0]= lasso.sumOfFitFractions() ;
		TGraph* ff = new TGraph(1,x,y);
		ff->SetName("SumOfFitFractions A");
		ff->SetTitle("");
		ff->GetXaxis()->SetTitle("#lambda");
		ff->GetXaxis()->SetTitleOffset(0.65);
		ff->GetYaxis()->SetTitle("Total Fit Fraction");
		ff->Draw("A*");
		ff->Write();
		
		y[0]= lasso.absSumOfInterferenceFractions() ;
		TGraph* iff = new TGraph(1,x,y);
		iff->SetName("AbsSumOfInterferenceFractions A");
		iff->SetTitle("");
		iff->GetXaxis()->SetTitle("#lambda");
		iff->GetXaxis()->SetTitleOffset(0.65);
		iff->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
		iff->Draw("A*");
		iff->Write();

		y[0]= lasso_bar.sumOfFitFractions() ;
		TGraph* ff_bar = new TGraph(1,x,y);
		ff_bar->SetName("SumOfFitFractions Abar");
		ff_bar->SetTitle("");
		ff_bar->GetXaxis()->SetTitle("#lambda");
		ff_bar->GetXaxis()->SetTitleOffset(0.65);
		ff_bar->GetYaxis()->SetTitle("Total Fit Fraction");
		ff_bar->Draw("A*");
		ff_bar->Write();
		
		y[0]= lasso_bar.absSumOfInterferenceFractions() ;
		TGraph* iff_bar = new TGraph(1,x,y);
		iff_bar->SetName("AbsSumOfInterferenceFractions Abar");
		iff_bar->SetTitle("");
		iff_bar->GetXaxis()->SetTitle("#lambda");
		iff_bar->GetXaxis()->SetTitleOffset(0.65);
		iff_bar->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
		iff_bar->Draw("A*");
		iff_bar->Write();

		y[0]= chi2_val ;
		TGraph* chi2 = new TGraph(1,x,y);
		chi2->SetName("Chi2");
		chi2->SetTitle("");
		chi2->GetXaxis()->SetTitle("#lambda");
		chi2->GetXaxis()->SetTitleOffset(0.65);
		chi2->GetYaxis()->SetTitle("#Chi^{2}/Bin");
		chi2->Draw("A*");
		chi2->Write();
	
		// fill tree
		pull_tree->Fill();
		paraFile->cd();
		pull_tree->SetDirectory(paraFile);
		pull_tree->Write();
		paraFile->Close();
		delete paraFile;
	}

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
    TString prefix = "";
    //TString prefix = "BsTaggingTool_";
    NamedParameter<double> min_year("min_year", 11);
    NamedParameter<double> max_year("max_year", 16);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);    

    /// Load files
    // Data
    int q_OS,f,Ds_ID;
    Short_t q_SS;
    double w_OS;
    Float_t w_SS;
    double sw;
    int run,year,Ds_finalState,trigger;
    double t,dt;
    
    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm.root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*DEC",1);
    tree_norm->SetBranchStatus("*PROB",1);
    tree_norm->SetBranchStatus("*OS",1);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("run",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("Ds_ID",1);

    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&w_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&w_SS);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree_norm->SetBranchAddress("Bs_DTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_DTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("Ds_ID",&Ds_ID);

    ///Make histograms
    int bins = 60;
    TH1D* h_w_OS_norm = new TH1D("h_w_OS_norm_comb","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1 = new TH1D("h_w_OS_norm_Run1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2 = new TH1D("h_w_OS_norm_Run2","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t0 = new TH1D("h_w_OS_norm_Run1_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t0 = new TH1D("h_w_OS_norm_Run2_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t1 = new TH1D("h_w_OS_norm_Run1_t1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t1 = new TH1D("h_w_OS_norm_Run2_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm = new TH1D("h_w_SS_norm_comb","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1 = new TH1D("h_w_SS_norm_Run1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2 = new TH1D("h_w_SS_norm_Run2","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t0 = new TH1D("h_w_SS_norm_Run1_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t0 = new TH1D("h_w_SS_norm_Run2_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t1 = new TH1D("h_w_SS_norm_Run1_t1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t1 = new TH1D("h_w_SS_norm_Run2_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm = new TH1D("h_q_OS_norm_comb","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1 = new TH1D("h_q_OS_norm_Run1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2 = new TH1D("h_q_OS_norm_Run2","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t0 = new TH1D("h_q_OS_norm_Run1_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t0 = new TH1D("h_q_OS_norm_Run2_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t1 = new TH1D("h_q_OS_norm_Run1_t1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t1 = new TH1D("h_q_OS_norm_Run2_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm = new TH1D("h_q_SS_norm_comb","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1 = new TH1D("h_q_SS_norm_Run1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2 = new TH1D("h_q_SS_norm_Run2","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t0 = new TH1D("h_q_SS_norm_Run1_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t0 = new TH1D("h_q_SS_norm_Run2_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t1 = new TH1D("h_q_SS_norm_Run1_t1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t1 = new TH1D("h_q_SS_norm_Run2_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm = new TH1D("h_q_f_norm_comb","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1 = new TH1D("h_q_f_norm_Run1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2 = new TH1D("h_q_f_norm_Run2","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t0 = new TH1D("h_q_f_norm_Run1_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t0 = new TH1D("h_q_f_norm_Run2_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t1 = new TH1D("h_q_f_norm_Run1_t1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t1 = new TH1D("h_q_f_norm_Run2_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm = new TH1D("h_t_norm_comb",";t (ps);Events (norm.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run1 = new TH1D("h_t_norm_Run1",";t (ps);Events (norm.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run2 = new TH1D("h_t_norm_Run2",";t (ps);Events (norm.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm = new TH1D("h_dt_norm_comb",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1 = new TH1D("h_dt_norm_Run1",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2 = new TH1D("h_dt_norm_Run2",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t0 = new TH1D("h_dt_norm_Run1_t0",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t0 = new TH1D("h_dt_norm_Run2_t0",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t1 = new TH1D("h_dt_norm_Run1_t1",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t1 = new TH1D("h_dt_norm_Run2_t1",";#sigma_{t} (ps);Events (norm.) ",bins,min_TAUERR,max_TAUERR);
    
    ///loop over data events
    for(int i=0; i< tree_norm->GetEntries(); i++)
    {    
        tree_norm->GetEntry(i);
        if(year < min_year || year > max_year) continue;
	if(Ds_ID>0)f = -1;
	else f = 1;

        h_t_norm->Fill(t,sw);
        h_dt_norm->Fill(dt,sw);
        h_q_OS_norm->Fill((double)q_OS,sw);
        h_q_SS_norm->Fill((double)q_SS,sw);
        h_q_f_norm->Fill((double)f,sw);
        if(q_OS != 0)h_w_OS_norm->Fill(w_OS,sw);
        if(q_SS != 0)h_w_SS_norm->Fill(w_SS,sw);
            
        if(run==1){
            h_t_norm_Run1->Fill(t,sw);
            h_dt_norm_Run1->Fill(dt,sw);
            h_q_OS_norm_Run1->Fill((double)q_OS,sw);
            h_q_SS_norm_Run1->Fill((double)q_SS,sw);
	    h_q_f_norm_Run1->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run1->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run1->Fill(w_SS,sw);
	    if(trigger == 0){
		h_dt_norm_Run1_t0->Fill(dt,sw);
            	h_q_OS_norm_Run1_t0->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t0->Fill((double)q_SS,sw);
        	h_q_f_norm_Run1_t0->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t0->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t0->Fill(w_SS,sw);
	    }
	    else if(trigger == 1){
		h_dt_norm_Run1_t1->Fill(dt,sw);
            	h_q_OS_norm_Run1_t1->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t1->Fill((double)q_SS,sw);
       		h_q_f_norm_Run1_t1->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t1->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t1->Fill(w_SS,sw);
	    }
        }
        else if(run==2){
            h_t_norm_Run2->Fill(t,sw);
            h_dt_norm_Run2->Fill(dt,sw);
            h_q_OS_norm_Run2->Fill((double)q_OS,sw);
            h_q_SS_norm_Run2->Fill((double)q_SS,sw);
       	    h_q_f_norm_Run2->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run2->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run2->Fill(w_SS,sw);
	    if(trigger == 0){
		h_dt_norm_Run2_t0->Fill(dt,sw);
            	h_q_OS_norm_Run2_t0->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run2_t0->Fill((double)q_SS,sw);
        	h_q_f_norm_Run2_t0->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run2_t0->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run2_t0->Fill(w_SS,sw);
	    }
	    else if(trigger == 1){
		h_dt_norm_Run2_t1->Fill(dt,sw);
            	h_q_OS_norm_Run2_t1->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run2_t1->Fill((double)q_SS,sw);
        	h_q_f_norm_Run2_t1->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run2_t1->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run2_t1->Fill(w_SS,sw);
	    }
        }
       
    }
    
    TFile* out = new TFile("Mistag_pdfs.root","RECREATE");
    h_t_norm->Write();
    h_dt_norm->Write();
    h_q_OS_norm->Write();
    h_w_OS_norm->Write();
    h_q_SS_norm->Write();
    h_w_SS_norm->Write();
    h_q_f_norm->Write();

    h_t_norm_Run1->Write();
    h_dt_norm_Run1->Write();
    h_q_OS_norm_Run1->Write();
    h_w_OS_norm_Run1->Write();
    h_q_SS_norm_Run1->Write();
    h_w_SS_norm_Run1->Write();
    h_q_f_norm_Run1->Write();

    h_dt_norm_Run1_t0->Write();
    h_q_OS_norm_Run1_t0->Write();
    h_w_OS_norm_Run1_t0->Write();
    h_q_SS_norm_Run1_t0->Write();
    h_w_SS_norm_Run1_t0->Write();
    h_q_f_norm_Run1_t0->Write();

    h_dt_norm_Run1_t1->Write();
    h_q_OS_norm_Run1_t1->Write();
    h_w_OS_norm_Run1_t1->Write();
    h_q_SS_norm_Run1_t1->Write();
    h_w_SS_norm_Run1_t1->Write();
    h_q_f_norm_Run1_t1->Write();

    h_t_norm_Run2->Write();
    h_dt_norm_Run2->Write();
    h_q_OS_norm_Run2->Write();
    h_w_OS_norm_Run2->Write();
    h_q_SS_norm_Run2->Write();
    h_w_SS_norm_Run2->Write();
    h_q_f_norm_Run2->Write();

    h_dt_norm_Run2_t0->Write();
    h_q_OS_norm_Run2_t0->Write();
    h_w_OS_norm_Run2_t0->Write();
    h_q_SS_norm_Run2_t0->Write();
    h_w_SS_norm_Run2_t0->Write();
    h_q_f_norm_Run2_t0->Write();

    h_dt_norm_Run2_t1->Write();
    h_q_OS_norm_Run2_t1->Write();
    h_w_OS_norm_Run2_t1->Write();
    h_q_SS_norm_Run2_t1->Write();
    h_w_SS_norm_Run2_t1->Write();
    h_q_f_norm_Run2_t1->Write();

    out->Write();
}

void produceIntegratorFile_CP(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;

    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
    DiskResidentEventList eventList(((string) integratorEventFile.ReplaceAll(".root","_CP.root")).c_str(),"RECREATE");

    for(int i = 0; i < eventListMC.size(); i++){
	DalitzEvent evt(eventListMC.getEvent(i));
	evt.CP_conjugateYourself();
	//evt.P_conjugateYourself();
	eventList.Add(evt);
    }

    eventList.save();
    return;
}

void makeIntegratorFileForToys(int step = 0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    DalitzEventList eventListPhsp,eventList,eventList_cut,eventList_cut_CP;    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.print();
    fas.getVal(eventListPhsp[0]);
    //fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);
    
    sg.FillEventList(eventList, IntegratorEvents);
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);
    
    for(int i = 0; i < eventList.size(); i++){
        if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 1.95 && sqrt(eventList[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventList[i].s(3,4)/(GeV*GeV)) < 1.2){
            eventList_cut.Add(eventList[i]);
            DalitzEvent evt(eventList[i]);
            evt.CP_conjugateYourself();
            eventList_cut_CP.Add(evt);
        }
    }
    
    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    TString outputName = (string)IntegratorEventFile;
    if(step>0) outputName.ReplaceAll(".root",("_" + anythingToString(step) + ".root").c_str());
    eventList_cut.saveAsNtuple((string)outputName);
    eventList_cut_CP.saveAsNtuple((string)outputName.ReplaceAll(".root","_CP.root"));
    return;
}

int main(int argc, char** argv){

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->ProcessLine(".x ../lhcbStyle.C");

  //produceMarginalPdfs();
  //produceIntegratorFile_CP();
  //makeIntegratorFileForToys(atoi(argv[1]));
    
  //for(int i = 0; i < 200; i++)ampFit(atoi(argv[1])+i);
  ampFit(atoi(argv[1]));

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
