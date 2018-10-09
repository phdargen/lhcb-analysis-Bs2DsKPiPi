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
#include "Neg2LLMultiConstraint.h"
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
#include "RooTrace.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooRandom.h"
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
#include "Mint/FullTimePdf.h"
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

        std::complex<double> result= polar((double) ( _re ),(double) (_im/360.*2.*pi ) ); 
//        std::complex<double> result(_re,static_cast<double>(_f) * _im); 
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
//         std::complex<double> result= polar((double) sqrt( 1.+  _r * (double) _sign),(double) (_delta/360.*2.*pi * (double) _sign) ); 
//  	    std::complex<double> result= polar((double) ( 1.+  _r * (double) _sign),(double) (_delta/360.*2.*pi * (double) _sign) ); 
//	    std::complex<double> result= polar((double) ( _r ),(double) (_delta/360.*2.*pi ) ); 

//  	complex<double> result= polar((double) ( 1.+  _r * (double) _sign)/sqrt(2.*(1.+ _r * _r)),(double) (_delta/360.*2.*pi * (double) _sign) ); 

// 	cout << (double)_r << "," << (double)_delta << "," << _sign << endl;

 	complex<double> result(1.,0.);

	result += (double) _sign * polar((double) _r,(double) (_delta/360.*2.*pi * (double) _sign) ); 

// 	cout << polar((double) _r,(double) (_delta/360.*2.*pi * (double) _sign) ) << endl;
// 	cout << result << endl;
// 	cout << result/sqrt(2.*(1.+ _r * _r)) << endl << endl;


	return result/sqrt(2.*(1.+ _r * _r));

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
    fas_2.multiply(scale);
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

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList, DalitzEventList& eventListData){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    int nBins = 20;
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

    TH1D* m_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
    TH1D* r_Kpi = (TH1D*)m_Kpi->Clone();
    TH1D* k_Kpi = (TH1D*)m_Kpi->Clone();
    TH1D* rk_Kpi = (TH1D*)m_Kpi->Clone();

    TH1D* m_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
    TH1D* r_Kpipi = (TH1D*)m_Kpipi->Clone();
    TH1D* k_Kpipi = (TH1D*)m_Kpipi->Clone();
    TH1D* rk_Kpipi = (TH1D*)m_Kpipi->Clone();

    TH1D* m_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
    TH1D* r_pipi = (TH1D*)m_pipi->Clone();
    TH1D* k_pipi = (TH1D*)m_pipi->Clone();
    TH1D* rk_pipi = (TH1D*)m_pipi->Clone();

    TH1D* m_Dspipi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
    TH1D* r_Dspipi = (TH1D*)m_Dspipi->Clone();
    TH1D* k_Dspipi = (TH1D*)m_Dspipi->Clone();
    TH1D* rk_Dspipi = (TH1D*)m_Dspipi->Clone();

    TH1D* m_Dspi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
    TH1D* r_Dspi = (TH1D*)m_Dspi->Clone();
    TH1D* k_Dspi = (TH1D*)m_Dspi->Clone();
    TH1D* rk_Dspi = (TH1D*)m_Dspi->Clone();

    vector< complex<double> > valK_Kpipi(nBins+1,0);
    vector<double> val1_Kpipi(nBins+1,0);
    vector<double> val2_Kpipi(nBins+1,0);

    vector< complex<double> > valK_Kpi(nBins+1,0);
    vector<double> val1_Kpi(nBins+1,0);
    vector<double> val2_Kpi(nBins+1,0);

    vector< complex<double> > valK_pipi(nBins+1,0);
    vector<double> val1_pipi(nBins+1,0);
    vector<double> val2_pipi(nBins+1,0);

    vector< complex<double> > valK_Dspipi(nBins+1,0);
    vector<double> val1_Dspipi(nBins+1,0);
    vector<double> val2_Dspipi(nBins+1,0);

    vector< complex<double> > valK_Dspi(nBins+1,0);
    vector<double> val1_Dspi(nBins+1,0);
    vector<double> val2_Dspi(nBins+1,0);

    for(unsigned int i=0; i<eventList.size(); i++){
        DalitzEvent evt = eventList.getEvent(i);

        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();

	int bin_Kpi = r_Kpi->FindBin(sqrt(evt.s(2,4)/(GeV*GeV))); 
	if(!(bin_Kpi < 0 || bin_Kpi > nBins)) {
	        valK_Kpi[bin_Kpi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Kpi[bin_Kpi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val2_Kpi[bin_Kpi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Kpipi = r_Kpipi->FindBin(sqrt(evt.sij(s234)/(GeV*GeV))); 
	if(!(bin_Kpipi < 0 || bin_Kpipi > nBins)) {
	        valK_Kpipi[bin_Kpipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Kpipi[bin_Kpipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Kpipi[bin_Kpipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_pipi = r_pipi->FindBin(sqrt(evt.s(3,4)/(GeV*GeV))); 
	if(!(bin_pipi < 0 || bin_pipi > nBins)) {
	        valK_pipi[bin_pipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_pipi[bin_pipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_pipi[bin_pipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Dspipi = r_Dspipi->FindBin(sqrt(evt.sij(s134)/(GeV*GeV))); 
	if(!(bin_Dspipi < 0 || bin_Dspipi > nBins)) {
	        valK_Dspipi[bin_Dspipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Dspipi[bin_Dspipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Dspipi[bin_Dspipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Dspi = r_Dspi->FindBin(sqrt(evt.s(1,3)/(GeV*GeV))); 
	if(!(bin_Dspi < 0 || bin_Dspi > nBins)) {
	        valK_Dspi[bin_Dspi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Dspi[bin_Dspi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Dspi[bin_Dspi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}

    }

    for(unsigned int i=0; i<eventListData.size(); i++){
	m_Kpipi->Fill(sqrt(eventListData[i].sij(s234)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Kpi->Fill(sqrt(eventListData[i].s(2,4)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_pipi->Fill(sqrt(eventListData[i].s(3,4)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Dspipi->Fill(sqrt(eventListData[i].sij(s134)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Dspi->Fill(sqrt(eventListData[i].s(1,3)/(GeV*GeV)),eventListData[i].getWeight()) ;
    }

   for(int i = 1; i <= nBins; i++){
	if(val1_Kpipi[i]>0){
		r_Kpipi->SetBinContent(i,sqrt(val2_Kpipi[i]/val1_Kpipi[i]));
		k_Kpipi->SetBinContent(i,abs(valK_Kpipi[i])/sqrt(val1_Kpipi[i])/sqrt(val2_Kpipi[i]));
		rk_Kpipi->SetBinContent(i,r_Kpipi->GetBinContent(i)*k_Kpipi->GetBinContent(i));
	}
	if(val1_Kpi[i]>0){
		r_Kpi->SetBinContent(i,sqrt(val2_Kpi[i]/val1_Kpi[i]));
		k_Kpi->SetBinContent(i,abs(valK_Kpi[i])/sqrt(val1_Kpi[i])/sqrt(val2_Kpi[i]));
		rk_Kpi->SetBinContent(i,r_Kpi->GetBinContent(i)*k_Kpi->GetBinContent(i));
	}
	if(val1_pipi[i]>0){
		r_pipi->SetBinContent(i,sqrt(val2_pipi[i]/val1_pipi[i]));
		k_pipi->SetBinContent(i,abs(valK_pipi[i])/sqrt(val1_pipi[i])/sqrt(val2_pipi[i]));
		rk_pipi->SetBinContent(i,r_pipi->GetBinContent(i)*k_pipi->GetBinContent(i));
	}
	if(val1_Dspipi[i]>0){
		r_Dspipi->SetBinContent(i,sqrt(val2_Dspipi[i]/val1_Dspipi[i]));
		k_Dspipi->SetBinContent(i,abs(valK_Dspipi[i])/sqrt(val1_Dspipi[i])/sqrt(val2_Dspipi[i]));
		rk_Dspipi->SetBinContent(i,r_Dspipi->GetBinContent(i)*k_Dspipi->GetBinContent(i));
	}	
	if(val1_Dspi[i]>0){
		r_Dspi->SetBinContent(i,sqrt(val2_Dspi[i]/val1_Dspi[i]));
		k_Dspi->SetBinContent(i,abs(valK_Dspi[i])/sqrt(val1_Dspi[i])/sqrt(val2_Dspi[i]));
		rk_Dspi->SetBinContent(i,r_Dspi->GetBinContent(i)*k_Dspi->GetBinContent(i));
	}
    }	
    
    TCanvas* c = new TCanvas();

    m_Kpipi->SetFillColor(kGray+1);
    m_Kpipi->SetLineColor(kGray+1);
    m_Kpipi->Scale(1./m_Kpipi->GetMaximum());
    m_Kpipi->SetMaximum(1.2);
    m_Kpipi->SetMinimum(0);
    m_Kpipi->Draw("hist");

    r_Kpipi->SetLineColor(kRed);
    r_Kpipi->Scale(1./r_Kpipi->GetMaximum());
    r_Kpipi->Draw("histcsame");

    k_Kpipi->SetLineColor(kBlue);
    k_Kpipi->Scale(1./k_Kpipi->GetMaximum());
    k_Kpipi->Draw("histcsame");

    rk_Kpipi->SetLineColor(kGreen+3);
    rk_Kpipi->Scale(1./rk_Kpipi->GetMaximum());
    rk_Kpipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Kpipi.eps").c_str());

    m_Kpi->SetFillColor(kGray+1);
    m_Kpi->SetLineColor(kGray+1);
    m_Kpi->Scale(1./m_Kpi->GetMaximum());
    m_Kpi->SetMaximum(1.2);
    m_Kpi->SetMinimum(0);    
    m_Kpi->Draw("hist");

    r_Kpi->SetLineColor(kRed);
    r_Kpi->Scale(1./r_Kpi->GetMaximum());
    r_Kpi->Draw("histcsame");

    k_Kpi->SetLineColor(kBlue);
    k_Kpi->Scale(1./k_Kpi->GetMaximum());
    k_Kpi->Draw("histcsame");

    rk_Kpi->SetLineColor(kGreen+3);
    rk_Kpi->Scale(1./rk_Kpi->GetMaximum());
    rk_Kpi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Kpi.eps").c_str());

    m_pipi->SetFillColor(kGray+1);
    m_pipi->SetLineColor(kGray+1);
    m_pipi->Scale(1./m_pipi->GetMaximum());
    m_pipi->SetMaximum(1.2);
    m_pipi->SetMinimum(0);
    m_pipi->Draw("hist");

    r_pipi->SetLineColor(kRed);
    r_pipi->Scale(1./r_pipi->GetMaximum());
    r_pipi->Draw("histcsame");

    k_pipi->SetLineColor(kBlue);
    k_pipi->Scale(1./k_pipi->GetMaximum());
    k_pipi->Draw("histcsame");

    rk_pipi->SetLineColor(kGreen+3);
    rk_pipi->Scale(1./rk_pipi->GetMaximum());
    rk_pipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_pipi.eps").c_str());

    m_Dspipi->SetFillColor(kGray+1);
    m_Dspipi->SetLineColor(kGray+1);
    m_Dspipi->Scale(1./m_Dspipi->GetMaximum());
    m_Dspipi->SetMaximum(1.2);
    m_Dspipi->SetMinimum(0);
    m_Dspipi->Draw("hist");

    r_Dspipi->SetLineColor(kRed);
    r_Dspipi->Scale(1./r_Dspipi->GetMaximum());
    r_Dspipi->Draw("histcsame");

    k_Dspipi->SetLineColor(kBlue);
    k_Dspipi->Scale(1./k_Dspipi->GetMaximum());
    k_Dspipi->Draw("histcsame");

    rk_Dspipi->SetLineColor(kGreen+3);
    rk_Dspipi->Scale(1./rk_Dspipi->GetMaximum());
    rk_Dspipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Dspipi.eps").c_str());

    m_Dspi->SetFillColor(kGray+1);
    m_Dspi->SetLineColor(kGray+1);
    m_Dspi->Scale(1./m_Dspi->GetMaximum());
    m_Dspi->SetMaximum(1.2);
    m_Dspi->SetMinimum(0);
    m_Dspi->Draw("hist");

    r_Dspi->SetLineColor(kRed);
    r_Dspi->Scale(1./r_Dspi->GetMaximum());
    r_Dspi->Draw("histcsame");

    k_Dspi->SetLineColor(kBlue);
    k_Dspi->Scale(1./k_Dspi->GetMaximum());
    k_Dspi->Draw("histcsame");

    rk_Dspi->SetLineColor(kGreen+3);
    rk_Dspi->Scale(1./rk_Dspi->GetMaximum());
    rk_Dspi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Dspi.eps").c_str());

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

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList, int N = -1){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    //fas.print();
    //fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    N = (N > 0) ? N : eventList.size();
    for(unsigned int i=0; i< N; i++){
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
    result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
/*
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "Dbar = " << result[4]<< endl;
    cout << "Sbar = " << result[5] << endl;
*/  
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
class FullAmpsPdfFlexiFastCPV : 
	public MINT::PdfBase<IDalitzEvent>, virtual public IDalitzPdf{
    
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
    const MINT::FitParameter& _r;
    const MINT::FitParameter& _delta;
    const MINT::FitParameter& _gamma;

    const MINT::FitParameter& _xm;
    const MINT::FitParameter& _ym;

    const MINT::FitParameter& _xp;
    const MINT::FitParameter& _yp;

    const MINT::FitParameter& _Gamma;
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
    NamedParameter<int> _useCartCoord;
    NamedParameter<int> _directCPV;

    // Limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    NamedParameter<double> _min_TAUERR;
    NamedParameter<double> _max_TAUERR;

    // Toy generation
    FitAmpIncoherentSum* _fasGen;
    SignalGenerator* _sg;
    FitAmpIncoherentSum* _fasGen_CP;
    SignalGenerator* _sg_CP;
    vector<int> _s234;


public:
    void parametersChanged(){
        _ampsSum->parametersChanged();
        _intA = (_ampsSum->ComplexIntegralForTags(1,1)).real();
        _intAbar = (_ampsSum->ComplexIntegralForTags(-1,-1)).real();        
        _intAAbar = _ampsSum->ComplexIntegralForTags(1,-1);

	if(_directCPV){
		_ampsSum_CP->parametersChanged();
		_intA_CP = (_ampsSum_CP->ComplexIntegralForTags(1,1)).real();
		_intAbar_CP = (_ampsSum_CP->ComplexIntegralForTags(-1,-1)).real();   
		_intAAbar_CP = _ampsSum_CP->ComplexIntegralForTags(1,-1);
	}
	else {
		_intA_CP = _intA;
		_intAbar_CP = _intAbar;   
		_intAAbar_CP = _intAAbar;
	}
    }
    void beginFit(){
        _ampsSum->beginFit();
        _ampsSum->redoIntegrator();

	if(_directCPV){      
		_ampsSum_CP->beginFit();
		_ampsSum_CP->redoIntegrator();
	}
	parametersChanged(); 
        printIntegralVals();
	_timePdfMaster->listFitParDependencies();
    }
    void endFit(){
        printIntegralVals();
        _ampsSum->endFit();
        if(_directCPV)_ampsSum_CP->endFit();
    }
    
    void printIntegralVals(){
        cout << "intSum = " << _ampsSum->getIntegralValue() << endl;
        cout << "intA = " << _intA << endl;
        cout << "intAbar = " << _intAbar << endl;
        cout << "mag intAAbar = " << std::abs(_intAAbar) << endl;
	cout << "phase intAAbar = " << std::arg(_intAAbar)/(2.*pi)*360. << endl;

	if(_directCPV){
		cout << "intSum_CP = " << _ampsSum_CP->getIntegralValue() << endl;
		cout << "intA_CP = " << _intA_CP << endl;
		cout << "intAbar_CP = " << _intAbar_CP << endl;
		cout << "mag intAAbar_CP = " << std::abs(_intAAbar_CP) << endl;
		cout << "phase intAAbar_CP = " << std::arg(_intAAbar_CP)/(2.*pi)*360. << endl;
	}
    }

    inline double un_normalised_noPs(IDalitzEvent& evt){

        const double f = (double) evt.getValueFromVector(2);
        _timePdfMaster->setAllObservablesAndFitParameters(evt);
                
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = (_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = (_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 0.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 0.;

        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA) ;
            amp_bar = phase_diff * _amps_bar->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA_CP) ;
            amp_bar_CP = phase_diff_CP * _amps_bar_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar_CP) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }

        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp)) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP)) ,
        		 2. * imag(amp_bar*conj(amp)) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP)) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
             ( _timePdfMaster->get_cosh_term_Val(evt)
             +  _timePdfMaster->get_cos_term_Val(evt)
             +  _timePdfMaster->get_sinh_term_Val(evt)
             +  _timePdfMaster->get_sin_term_Val(evt)
             )  ; //* _timePdfMaster->get_marginalPdfs_product(evt); //* _timePdfMaster->get_marginalPdfs_Val(evt);

        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double val = un_normalised_noPs(evt);
    
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        const complex<double> int_interference =  _intAAbar  * phase_diff/sqrt(_intA)/sqrt(_intAbar) ;
        const complex<double> int_interference_CP = _intAAbar_CP  * phase_diff_CP/sqrt(_intA_CP)/sqrt(_intAbar_CP) ;

        // C,Cbar,D,Dbar,S,Sbar
/*
  	double r =(double)_r;
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
*/

         _timePdfMaster->setCP_coeff(
			(1. +  norm(phase_diff)),
			(1. +  norm(phase_diff_CP)),
        		(1. -  norm(phase_diff)),
        		-(1. -  norm(phase_diff_CP)) , 
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
    
    inline double getVal_timeIntegrated(IDalitzEvent& evt){

        const double f = (double) evt.getValueFromVector(2);
        _timePdfMaster->setAllObservablesToMean(evt);
                
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 0.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 0.;

        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA) ;
            amp_bar = phase_diff * _amps_bar->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA_CP) ;
            amp_bar_CP = phase_diff_CP * _amps_bar_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar_CP) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }
       
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp)) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP)) ,
        		 2. * imag(amp_bar*conj(amp)) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP)) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
        (     _timePdfMaster->get_cosh_term_Integral(evt)
            +  _timePdfMaster->get_cos_term_Integral(evt)
            +  _timePdfMaster->get_sinh_term_Integral(evt)
            +  _timePdfMaster->get_sin_term_Integral(evt) );

        return val;
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

    vector<double> calculateCP_coeff(bool print = true){

	vector<double> result;
	if(_intA == -1)this->parametersChanged();
    
	double delta_0 = std::arg(_intAAbar)/(2.*pi)*360.;
	double delta_0_CP = std::arg(_intAAbar_CP)/(2.*pi)*360.;

	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        const complex<double> int_interference =  phase_diff * _intAAbar/sqrt(_intA)/sqrt(_intAbar) ;
        const complex<double> int_interference_CP = phase_diff_CP * _intAAbar_CP/sqrt(_intA_CP)/sqrt(_intAbar_CP) ;

        result.push_back((1. -  norm(phase_diff))/(1. +  norm(phase_diff)));
        result.push_back(-(1. -  norm(phase_diff_CP))/(1. +  norm(phase_diff_CP)));
 
        result.push_back((- int_interference.real() )/(1. +  norm(phase_diff)));
        result.push_back((- int_interference_CP.real() )/(1. +  norm(phase_diff_CP))); 

        result.push_back((int_interference.imag())/(1. +  norm(phase_diff)));
        result.push_back((- int_interference_CP.imag() )/(1. +  norm(phase_diff_CP)));   

        result.push_back(abs(_intAAbar/sqrt(_intA)/sqrt(_intAbar))/2.);   

	if(print){
		cout << "CP coeffs:: (C,Cbar,D,Dbar,S,Sbar,k) " << endl << result << endl;
		cout <<  abs(phase_diff) << endl;
		cout <<  arg(phase_diff)/(2*pi)*360. << endl;
		cout <<  abs(phase_diff_CP) << endl;
		cout <<  arg(phase_diff_CP)/(2*pi)*360. << endl<< endl;
	}
	return result;
    }

    double getCPcoeffChi2(vector<double>& coeff,double prec){
	vector<double> val = calculateCP_coeff(false);
	double chi2  = 0.;
	for(int i= 0; i< coeff.size();i++){
		chi2+= pow((coeff[i]-val[i])/prec,2);
	}
	return chi2;
    }

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_SS(evt);
    }

    std::pair<double, double> getCalibratedMistag(double eta,double avg_eta,double p0,double p1,double delta_p0,double delta_p1 ){
        return _timePdfMaster->getCalibratedMistag(eta, avg_eta, p0, p1, delta_p0, delta_p1);
    }
    
    double getCalibratedResolution(double& dt){
        return _timePdfMaster->getCalibratedResolution(dt);
    }

    RooDataSet* sampleEvents(int N = 10000){
	return _timePdfMaster->sampleEvents(N);
    }

    inline double getSampledPdfVal(IDalitzEvent& evt){
	return _timePdfMaster->getSamplingPdfVal(evt);
    }

    void generateBkgToys(int N, DalitzEventList& eventListData){

	DalitzEventList eventList,eventListDataSideband;

	for(int i=0; i< eventListData.size(); i++)
	{	
		DalitzEvent evt = eventListData[i];
		if(abs(evt.getValueFromVector(9) -5370) > 60)eventListDataSideband.Add(evt);
	}
	int N_sample = eventListDataSideband.size();

	vector<int> b_indices;
	while( b_indices.size() < N )b_indices.push_back(TMath::Nint(gRandom->Uniform(0,N_sample)));
	sort(b_indices.begin(), b_indices.end());
	
	TRandom3 rndm;
	for(int i=0; i< N; i++)
	{	
		DalitzEvent evt = eventListDataSideband[b_indices[i]];
		eventList.Add(evt);
	}
	saveEventListToFile(eventList);
    }

    DalitzEvent generateWeightedEvent(){

	while(true){
		double t_MC = gRandom->Exp(1./_Gamma);
                if(t_MC > _max_TAU || t_MC < _min_TAU)continue;

		int f_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;		

		counted_ptr<IDalitzEvent> evtPtr;

		if(f_MC > 0)evtPtr = _sg->newEvent();
		else evtPtr = _sg_CP->newEvent();

		DalitzEvent evt(evtPtr.get());
                if(!((sqrt(evt.sij(_s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2)))continue;

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		double dt_MC = marginal_vals[0] ;
		double eta_OS_MC = marginal_vals[1] ;
		double eta_SS_MC = marginal_vals[2] ;
	
		// true flavor
		int q_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;

 	        int q_SS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
         	int q_OS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
 
		q_OS_MC = (gRandom->Uniform() < 0.5) ? - q_OS_MC : q_OS_MC;
		q_SS_MC = (gRandom->Uniform() < 0.5) ? - q_SS_MC : q_SS_MC;
		
		eta_OS_MC = (q_OS_MC == 0) ? 0.5 : eta_OS_MC;
		eta_SS_MC = (q_SS_MC == 0) ? 0.5 : eta_SS_MC;

// 		if(f_MC<0)evt.CP_conjugateYourself();
		evt.setValueInVector(0, t_MC);
		evt.setValueInVector(1, dt_MC);
		evt.setValueInVector(2, f_MC);
		evt.setValueInVector(3, q_OS_MC);
		evt.setValueInVector(4, eta_OS_MC);
		evt.setValueInVector(5, q_SS_MC);
		evt.setValueInVector(6, eta_SS_MC);

		evt.setGeneratorPdfRelativeToPhaseSpace(_Gamma * evt.getGeneratorPdfRelativeToPhaseSpace() * (exp(-t_MC*_Gamma) / ( ( exp(-_min_TAU*_Gamma) - exp(-_max_TAU*_Gamma) ))));
		return evt;
	}
    }

    DalitzEventList generateToys(int N = 10000, int run = -1 , int trigger = -1){

	time_t startTime = time(0);

	cout << "Generating " << N << " events" << endl;
	DalitzEventList eventList;
	if(_fasGen==0){
		        /// Simple amplitude model for importance sampling
   		        NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
        		DalitzEventPattern pat(EventPattern.getVector());
		        DalitzEventPattern pat_CP = pat.makeCPConjugate();     		

			_fasGen = new FitAmpIncoherentSum((DalitzEventPattern)pat);
        		_fasGen->print();
			{
				DalitzEventList eventListPhsp;
				eventListPhsp.generatePhaseSpaceEvents(200000,pat);
				_fasGen->normalizeAmps(eventListPhsp);
			}
	
			_fasGen_CP = new FitAmpIncoherentSum(*_fasGen);
			_fasGen_CP->CPConjugateSameFitParameters();

			{
				DalitzEventList eventListPhsp_CP;
				eventListPhsp_CP.generatePhaseSpaceEvents(2,pat_CP);
				_fasGen_CP->getVal(eventListPhsp_CP[0]);
			}

		        _sg = new SignalGenerator(pat,_fasGen);
		        _sg_CP = new SignalGenerator(pat_CP,_fasGen_CP);
	}

	/// Estimate max val
	vector<double> vals;	
	for(int i = 0; i < 100000; i++){
		DalitzEvent evt = generateWeightedEvent();
		double val = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); ///_timePdfMaster->get_marginalPdfs_product(evt);
		vals.push_back(val);	
	}

	cout << "Now calculating maximum val " << vals.size() << endl;
	double amax,pmax;
	generalisedPareto_estimateMaximum(vals,0.999,amax,pmax);
	
	double pdf_max = 1.;
	if(!TMath::IsNaN(pmax) && pmax > 0 && pmax < 100 * amax)pdf_max = pmax;
	else if(!TMath::IsNaN(amax))pdf_max = amax;
	if(amax > pmax)pdf_max = amax;
	// for safety
 	pdf_max *= 1.5;	

	cout << "pdf_max " << pdf_max << endl;

	int N_gen = 0;
	int N_tot = 0;
	while(true){
			DalitzEvent evt = generateWeightedEvent();
			double pdfVal = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); // /_timePdfMaster->get_marginalPdfs_product(evt) ;
			
			const double height = gRandom->Uniform(0,pdf_max);
			
			///Safety check on the maxmimum generated height
			if( pdfVal > pdf_max ){
				std::cout << "ERROR: PDF above determined maximum." << std::endl;
				std::cout << pdfVal << " > " << pdf_max << std::endl;
				pdf_max = pdf_max * 2.;
			}
			
			///Hit-and-miss
			if( height < pdfVal ) { 
				evt.setValueInVector(7, run);
				evt.setValueInVector(8, trigger);
				eventList.Add(evt);
				N_gen++;
				if (0ul == (N_gen % 500ul)) cout << "Generated event " << N_gen << "/" << N << endl;
			}		
			N_tot ++;
			if(N_gen == N)break;
	}

	cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
 	cout << "Generated " << N_gen << " events ! Efficiecy = " << (double)N_gen/(double)N_tot << endl;

	saveEventListToFile(eventList);
 	return eventList;
   }

    void saveEventListToFile(DalitzEventList& eventList, string name = "toys.root"){

	    TFile* out = new TFile(name.c_str(),"RECREATE");
	    TTree* tree = new TTree("DecayTree","DecayTree");
    		
	    double t,dt;
	    int Ds_ID;
	    int q_OS,q,q_SS;
	    double eta_OS;
	    double eta_SS;
	    double sw,w;
	    int run,trigger;
	    int year;
		
	    double K[4];
	    double pip[4];
	    double pim[4];
	    double Ds[4];
	    double mB,m_Kpipi,m_Kpi,m_pipi;

    	    TBranch* br_mB = tree->Branch( "Bs_DTF_MM", &mB, "Bs_DTF_MM/D" );
    	    TBranch* br_sw = tree->Branch( "N_Bs_sw", &sw, "N_Bs_sw/D" );
    	    TBranch* br_w = tree->Branch( "weight", &w, "weight/D" );

    	    TBranch* br_t = tree->Branch( "Bs_BsDTF_TAU", &t, "Bs_BsDTF_TAU/D" );
    	    TBranch* br_dt = tree->Branch( "Bs_BsDTF_TAUERR", &dt, "Bs_BsDTF_TAUERR/D" );

	    TBranch* br_Ds_ID = tree->Branch("Ds_ID",&Ds_ID,"Ds_ID/I");

	    TBranch* br_q_OS =tree->Branch("OS_Combination_DEC",&q_OS,"OS_Combination_DEC/I");
	    TBranch* br_eta_OS =  tree->Branch("OS_Combination_PROB",&eta_OS,"OS_Combination_PROB/D");
	    TBranch* br_q_SS = tree->Branch("SS_Kaon_DEC",&q_SS,"SS_Kaon_DEC/I");
	    TBranch* br_eta_SS = tree->Branch("SS_Kaon_PROB",&eta_SS,"SS_Kaon_PROB/D");

	    TBranch* br_run = tree->Branch("run",&run,"run/I");
    	    TBranch* br_year = tree->Branch( "year", &year, "year/I" );
	    TBranch* br_trigger = tree->Branch("TriggerCat",&trigger,"TriggerCat/I");
	    
	    TBranch* br_K0 = tree->Branch("BsDTF_Kplus_PX",&K[0],"BsDTF_Kplus_PX/D");
	    TBranch* br_K1 =tree->Branch("BsDTF_Kplus_PY",&K[1],"BsDTF_Kplus_PY/D");
	    TBranch* br_K2 =tree->Branch("BsDTF_Kplus_PZ",&K[2],"BsDTF_Kplus_PZ/D");
	    TBranch* br_K3 =tree->Branch("BsDTF_Kplus_PE",&K[3],"BsDTF_Kplus_PE/D");
	
            TBranch* br_pip0 =tree->Branch("BsDTF_piplus_PX",&pip[0],"BsDTF_piplus_PX/D");
	    TBranch* br_pip1 =tree->Branch("BsDTF_piplus_PY",&pip[1],"BsDTF_piplus_PY/D");
	    TBranch* br_pip2 =tree->Branch("BsDTF_piplus_PZ",&pip[2],"BsDTF_piplus_PZ/D");
	    TBranch* br_pip3 =tree->Branch("BsDTF_piplus_PE",&pip[3],"BsDTF_piplus_PE/D");
	
	    TBranch* br_pim0 =tree->Branch("BsDTF_piminus_PX",&pim[0],"BsDTF_piminus_PX/D");
	    TBranch* br_pim1 =tree->Branch("BsDTF_piminus_PY",&pim[1],"BsDTF_piminus_PY/D");
	    TBranch* br_pim2 =tree->Branch("BsDTF_piminus_PZ",&pim[2],"BsDTF_piminus_PZ/D");
	    TBranch* br_pim3 =tree->Branch("BsDTF_piminus_PE",&pim[3],"BsDTF_piminus_PE/D");
	
	    TBranch* br_Ds0 =tree->Branch("BsDTF_Ds_PX",&Ds[0],"BsDTF_Ds_PX/D");
	    TBranch* br_Ds1 =tree->Branch("BsDTF_Ds_PY",&Ds[1],"BsDTF_Ds_PY/D");
	    TBranch* br_Ds2 =tree->Branch("BsDTF_Ds_PZ",&Ds[2],"BsDTF_Ds_PZ/D");
	    TBranch* br_Ds3 =tree->Branch("BsDTF_Ds_PE",&Ds[3],"BsDTF_Ds_PE/D");

    	    TBranch* br_m_Kpipi = tree->Branch( "m_Kpipi", &m_Kpipi, "m_Kpipi/D" );
    	    TBranch* br_m_Kpi = tree->Branch( "m_Kpi", &m_Kpi, "m_Kpi/D" );
    	    TBranch* br_m_pipi = tree->Branch( "m_pipi", &m_pipi, "m_pipi/D" );

	    for(int i= 0; i < eventList.size(); i++){

		t = eventList[i].getValueFromVector(0);
		dt = eventList[i].getValueFromVector(1);

		Ds_ID = - eventList[i].getValueFromVector(2);
		
		q_OS = eventList[i].getValueFromVector(3);
		eta_OS = eventList[i].getValueFromVector(4);

		q_SS = eventList[i].getValueFromVector(5);
		eta_SS = eventList[i].getValueFromVector(6);

		run = eventList[i].getValueFromVector(7);
		if(run == 2) year = 16;
		else year = 12;
		trigger = eventList[i].getValueFromVector(8);

		mB = eventList[i].p(0).M(); //eventList[i].getValueFromVector(9);
		sw = eventList[i].getWeight();
		w = eventList[i].getWeight();

		Ds[0] = eventList[i].p(1).Px()/MeV;
		Ds[1] = eventList[i].p(1).Py()/MeV;
		Ds[2] = eventList[i].p(1).Pz()/MeV;
		Ds[3] = eventList[i].p(1).E()/MeV;

		K[0] = eventList[i].p(2).Px()/MeV;
		K[1] = eventList[i].p(2).Py()/MeV;
		K[2] = eventList[i].p(2).Pz()/MeV;
		K[3] = eventList[i].p(2).E()/MeV;

		pip[0] = eventList[i].p(3).Px()/MeV;
		pip[1] = eventList[i].p(3).Py()/MeV;
		pip[2] = eventList[i].p(3).Pz()/MeV;
		pip[3] = eventList[i].p(3).E()/MeV;

		pim[0] = eventList[i].p(4).Px()/MeV;
		pim[1] = eventList[i].p(4).Py()/MeV;
		pim[2] = eventList[i].p(4).Pz()/MeV;
		pim[3] = eventList[i].p(4).E()/MeV;

		m_Kpipi = sqrt(eventList[i].sij(_s234))/GeV;        
		m_Kpi = sqrt(eventList[i].s(2,4))/GeV;        
		m_pipi = sqrt(eventList[i].s(3,4))/GeV;        

		tree->Fill();
	     }

	     tree->Write();
	     out->Write();
	     out->Close();
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
                const MINT::FitParameter& r,const MINT::FitParameter& delta,const MINT::FitParameter& gamma,
                const MINT::FitParameter& xm,const MINT::FitParameter& ym,const MINT::FitParameter& xp,const MINT::FitParameter& yp,
		const MINT::FitParameter& Gamma, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
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
    _xm(xm),_ym(ym),_xp(xp),_yp(yp),
    _Gamma(Gamma),
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
    _marginalPdfsPrefix(marginalPdfsPrefix),
    _useCartCoord("FullAmpsPdfFlexiFastCPV::useCartCoord",1),
    _directCPV("FullAmpsPdfFlexiFastCPV::directCPV",0),
    _min_TAU("min_TAU", 0.4),
    _max_TAU("max_TAU", 10.),
    _min_TAUERR("min_TAUERR", 0.),
    _max_TAUERR("max_TAUERR", 0.1),
    _fasGen(0),
    _sg(0),
    _fasGen_CP(0),
    _sg_CP(0)
    {
		        _timePdfMaster = new TimePdfMaster(_Gamma, _dGamma, _dm
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

			_timePdfMaster->setAllFitParameters();

        		_s234.push_back(2);
		        _s234.push_back(3);
        		_s234.push_back(4);
    }
};


// Bkg PDF
class FullAmpsPdfFlexiFastBkg : 
	public MINT::PdfBase<IDalitzEvent>, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps;
    AmpsPdfFlexiFast* _amps_CP;

public:
    void parametersChanged(){
        _amps->parametersChanged();
        _amps_CP->parametersChanged();
    }
    void beginFit(){
        _amps->beginFit();
        _amps_CP->beginFit();
    }
    void endFit(){
        _amps->endFit();
        _amps_CP->endFit();
    }
    
    inline double un_normalised_noPs(IDalitzEvent& evt){
        const double f = (double) evt.getValueFromVector(2);
        if(f>0) return _amps->un_normalised_noPs(evt);
	else return _amps_CP->un_normalised_noPs(evt);
    }
    
    virtual double getVal(IDalitzEvent& evt){
        const double f = (double) evt.getValueFromVector(2);
        if(f>0) return _amps->getVal(evt);
	else return _amps_CP->getVal(evt);
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

    virtual DalitzHistoSet histoSet(){return _amps->histoSet();}

    FullAmpsPdfFlexiFastBkg(AmpsPdfFlexiFast* amps, AmpsPdfFlexiFast* amps_CP ):
    _amps(amps),_amps_CP(amps_CP)
    { ;  }
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


void ampFit(int step=0, string mode = "fit"){

    /// Generate list of amplitudes
    FitAmplitude::AutogenerateFitFile();

    /// Options
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    RooRandom::randomGenerator()->SetSeed(seed);
    TString prefix = "";
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    DalitzEventPattern pat_CP = pat.makeCPConjugate();

    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/signal_tagged.root", (char*) 0);
    NamedParameter<string> InputGenMCFile("InputGenMCFile", (std::string) "/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_DsKpipi_CPV.root", (char*) 0);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);

    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);
    NamedParameter<double> w_max("w_max", 0.5);

    NamedParameter<int>  doPlots("doPlots", 1);
    NamedParameter<int>  nBins("nBins", 50);
    NamedParameter<int>  nBinst("nBinst", 50);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 8);

    NamedParameter<int>  randomizeStartVals("randomizeStartVals", 0);
    NamedParameter<int>  doSimFit("doSimFit", 0);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  doScanParam("doScanParam", 0);
    NamedParameter<string> scanParam("scanParam", (std::string) "", (char*) 0);
    NamedParameter<double>  scanParam_min("scanParam_min", 0.);
    NamedParameter<double>  scanParam_max("scanParam_max", 1.);
    NamedParameter<double>  scanParam_steps("scanParam_steps", 100.);

    NamedParameter<int>  useLASSO("useLASSO", 0);
    NamedParameter<double>  lambda("lambda", 1.);

    NamedParameter<int>  initCPcoeff("initCPcoeff", 0);
    NamedParameter<int>  fitGenMC("fitGenMC", 0);
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);
    NamedParameter<int> doFractions("doFractions",0);
    NamedParameter<int> doFractionsErr("doFractionsErr",0);

    NamedParameter<int>  doToyStudy("doToyStudy", 0);
    NamedParameter<double> N_scale_toys("N_scale_toys", 1);

    NamedParameter<int>  useGaussConstrainsTagging("useGaussConstrainsTagging", 0);

    NamedParameter<int>  doAccSystematics("doAccSystematics", 0);
    NamedParameter<int>  useCholDec("useCholDec", 0);
    NamedParameter<int>  varPerParChol("varPerParChol", 100);
    int chol_index = (step-1)/varPerParChol ;

    NamedParameter<int>  doRadiusBWSystematic("doRadiusBWSystematic", 0);
    NamedParameter<string> doSystematic("doSystematic", (std::string) "", (char*) 0);
    NamedParameter<string> addAmpName("addAmpName", (std::string) "", (char*) 0);

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;
    TString integratorEventFile_CP = (string) IntegratorEventFile;
    integratorEventFile_CP.ReplaceAll(".root","_CP.root");
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

    /// Define amplitude model
    DalitzEventList eventListPhsp,eventListPhsp_CP;
    eventListPhsp.generatePhaseSpaceEvents(2,pat);
    eventListPhsp_CP.generatePhaseSpaceEvents(2,pat_CP);

    FitAmpSum fas_tmp((DalitzEventPattern)pat);
    //if(randomizeStartVals)fas_tmp.randomizeStartVals(seed);
    //if(randomizeStartVals)fas_tmp.randomizePhaseStartVals(seed);
    FitAmpIncoherentSum fasBkg(pat,"","Bkg_");

    /// Normalize amps
    {
    	DalitzEventList eventListNorm;
  	TFile *file =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
  	TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
  	eventListNorm.fromNtuple(tree,0.5);
  	fas_tmp.normalizeAmps(eventListNorm);
        fasBkg.normalizeAmps(eventListNorm);
	file->Close();
    }
    
    ///Choose reference amp
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1270)+");
    FitAmpSum fas(*List_1);
    FitAmpSum fas_bar(*List_1);
    
    /// Define relative decay modes    
    // A
    FitParameter a_K1_1400_Amp("a_K1_1400_Amp",1,1,0.01);
    FitParameter a_K1_1400_Phase("a_K1_1400_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_K1_1400 = new AmpRatio(a_K1_1400_Amp,a_K1_1400_Phase);

    FitParameter a_Ks_1410_Amp("a_Ks_1410_Amp",1,1,0.01);
    FitParameter a_Ks_1410_Phase("a_Ks_1410_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_Ks_1410 = new AmpRatio(a_Ks_1410_Amp,a_Ks_1410_Phase);

    FitParameter a_K_1460_Amp("a_K_1460_Amp",1,1,0.01);
    FitParameter a_K_1460_Phase("a_K_1460_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_K_1460 = new AmpRatio(a_K_1460_Amp,a_K_1460_Phase);

    FitParameter a_NS_Ks_Amp("a_NS_Ks_Amp",1,1,0.01);
    FitParameter a_NS_Ks_Phase("a_NS_Ks_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_Ks = new AmpRatio(a_NS_Ks_Amp,a_NS_Ks_Phase);

    FitParameter a_NS_sigma_Amp("a_NS_sigma_Amp",1,1,0.01);
    FitParameter a_NS_sigma_Phase("a_NS_sigma_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_sigma = new AmpRatio(a_NS_sigma_Amp,a_NS_sigma_Phase);

    FitParameter a_NS_rho_Amp("a_NS_rho_Amp",1,1,0.01);
    FitParameter a_NS_rho_Phase("a_NS_rho_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_rho = new AmpRatio(a_NS_rho_Amp,a_NS_rho_Phase);

    FitParameter a_sys_Amp("a_sys_Amp",1,1,0.01);
    FitParameter a_sys_Phase("a_sys_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_sys = new AmpRatio(a_sys_Amp,a_sys_Phase);

    // Abar
    FitParameter abar_K1_1400_Amp("abar_K1_1400_Amp",1,1,0.01);
    FitParameter abar_K1_1400_Phase("abar_K1_1400_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_K1_1400 = new AmpRatio(abar_K1_1400_Amp,abar_K1_1400_Phase);

    FitParameter abar_Ks_1410_Amp("abar_Ks_1410_Amp",1,1,0.01);
    FitParameter abar_Ks_1410_Phase("abar_Ks_1410_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_Ks_1410 = new AmpRatio(abar_Ks_1410_Amp,abar_Ks_1410_Phase);

    FitParameter abar_K_1460_Amp("abar_K_1460_Amp",1,1,0.01);
    FitParameter abar_K_1460_Phase("abar_K_1460_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_K_1460 = new AmpRatio(abar_K_1460_Amp,abar_K_1460_Phase);

    FitParameter abar_NS_Ks_Amp("abar_NS_Ks_Amp",1,1,0.01);
    FitParameter abar_NS_Ks_Phase("abar_NS_Ks_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_Ks = new AmpRatio(abar_NS_Ks_Amp,abar_NS_Ks_Phase);

    FitParameter abar_NS_sigma_Amp("abar_NS_sigma_Amp",1,1,0.01);
    FitParameter abar_NS_sigma_Phase("abar_NS_sigma_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_sigma = new AmpRatio(abar_NS_sigma_Amp,abar_NS_sigma_Phase);

    FitParameter abar_NS_rho_Amp("abar_NS_rho_Amp",1,1,0.01);
    FitParameter abar_NS_rho_Phase("abar_NS_rho_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_rho = new AmpRatio(abar_NS_rho_Amp,abar_NS_rho_Phase);

    FitParameter abar_sys_Amp("abar_sys_Amp",1,1,0.01);
    FitParameter abar_sys_Phase("abar_sys_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_sys = new AmpRatio(abar_sys_Amp,abar_sys_Phase);

    /// Randomize start vals
    MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
    if(randomizeStartVals){
	for(int i=0; i < mps->size(); i++){
		if(((FitParameter*)mps->getParPtr(i))->iFixInit()!=0)continue;
		if(A_is_in_B("abar_",((FitParameter*)mps->getParPtr(i))->name()) || A_is_in_B("a_",((FitParameter*)mps->getParPtr(i))->name()) ){
			double val = A_is_in_B("_Amp",((FitParameter*)mps->getParPtr(i))->name() ) ? gRandom->Uniform(0.5,1.5) * mps->getParPtr(i)->mean() : gRandom->Uniform(-180,180);
			mps->getParPtr(i)->setCurrentFitVal(val);
			((FitParameter*)mps->getParPtr(i))->setInit(val);	
			cout << "Setting " << ((FitParameter*)mps->getParPtr(i))->name() << "  to  " << val << endl; 
		}
	}
    } 

    /// Add amps to A and Abar
    AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K(1)(1400)+",a_K1_1400,abar_K1_1400);
    
    if(a_Ks_1410_Amp.iFixInit() != 1 && abar_Ks_1410_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K*(1410)+",a_Ks_1410,abar_Ks_1410);
    else if(a_Ks_1410_Amp.iFixInit() != 1 && abar_Ks_1410_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"K*(1410)+",a_Ks_1410);
    else if(a_Ks_1410_Amp.iFixInit() == 1 && abar_Ks_1410_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"K*(1410)+",abar_Ks_1410);
    
    if(a_K_1460_Amp.iFixInit() != 1 && abar_K_1460_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K(1460)+",a_K_1460,abar_K_1460);
    else if(a_K_1460_Amp.iFixInit() != 1 && abar_K_1460_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"K(1460)+",a_K_1460);
    else if(a_K_1460_Amp.iFixInit() == 1 && abar_K_1460_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp,fas_bar,"K(1460)+",abar_K_1460);
    
    if(a_NS_Ks_Amp.iFixInit() != 1 && abar_NS_Ks_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"(->Ds-,pi+),K*(892)0(->K+,pi-)",a_NS_Ks,abar_NS_Ks);
    else if(a_NS_Ks_Amp.iFixInit() != 1 && abar_NS_Ks_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"(->Ds-,pi+),K*(892)0(->K+,pi-)",a_NS_Ks);
    else if(a_NS_Ks_Amp.iFixInit() == 1 && abar_NS_Ks_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"(->Ds-,pi+),K*(892)0(->K+,pi-)",abar_NS_Ks);
    
    if(a_NS_rho_Amp.iFixInit() != 1 && abar_NS_rho_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"(->Ds-,K+),rho(770)0(->pi+,pi-)",a_NS_rho,abar_NS_rho);
    else if(a_NS_rho_Amp.iFixInit() != 1 && abar_NS_rho_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"(->Ds-,K+),rho(770)0(->pi+,pi-)",a_NS_rho);
    else if(a_NS_rho_Amp.iFixInit() == 1 && abar_NS_rho_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"(->Ds-,K+),rho(770)0(->pi+,pi-)",abar_NS_rho);

    if(a_sys_Amp.iFixInit() != 1 && abar_sys_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,(string) addAmpName,a_sys,abar_sys);
    else if(a_sys_Amp.iFixInit() != 1 && abar_sys_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,(string) addAmpName,a_sys);
    else if(a_sys_Amp.iFixInit() == 1 && abar_sys_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,(string) addAmpName,abar_sys);
   
    fas.print();
    fas_bar.print();

    /// Define B -> f amplitude        
    fas.setTag(1);
    /// Define Bbar -> f amplitude
    fas_bar.setTag(-1);

    /// CP conjugate amplitudes
    FitAmpSum fas_CP(fas);
    fas_CP.CPConjugateSameFitParameters();

    FitAmpSum fas_bar_CP(fas_bar);
    fas_bar_CP.CPConjugateSameFitParameters();

    FitAmpIncoherentSum fasBkg_CP(fasBkg);
    fasBkg_CP.CPConjugateSameFitParameters();

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

    /// Bkg
    FitParameter sigfraction("SigFraction",1,0.999999,0.01);
    AmpsPdfFlexiFast ampsBkg(pat, &fasBkg, 0, integPrecision*10,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsBkg_CP(pat_CP, &fasBkg_CP, 0, integPrecision*10,integMethod, (std::string) integratorEventFile_CP);
    FullAmpsPdfFlexiFastBkg pdf_bkg(&ampsBkg,&ampsBkg_CP);

    /// Fit parameters
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);

    FitParameter xm("xm",1,0,0.01);
    FitParameter ym("ym",1,0,0.01); 
    FitParameter xp("xp",1,0,0.01);
    FitParameter yp("yp",1,0,0.01); 

    if(doScanParam){
	for(int i=0; i < mps->size(); i++){
		if(A_is_in_B(scanParam,((FitParameter*)mps->getParPtr(i))->name())){
			double val = scanParam_min + step * (scanParam_max - scanParam_min)/scanParam_steps; 
			((FitParameter*)mps->getParPtr(i))->setInit(val);
			((FitParameter*)mps->getParPtr(i))->fix();		
			cout << "Setting " << ((FitParameter*)mps->getParPtr(i))->name() << "  to  " << val << endl; 
		}
	}
    } 

    /// Set blinded fit parameters to true values for toy generation
    if(mode == "gen"){
	double real_val = (double)r + r.blinding();
	r.setCurrentFitVal(real_val);
	r.setInit(real_val);

	real_val = (double)delta + delta.blinding();
	delta.setCurrentFitVal(real_val);
	delta.setInit(real_val);

	real_val = (double)gamma + gamma.blinding();
	gamma.setCurrentFitVal(real_val);
	gamma.setInit(real_val);

	real_val = (double)xm + xm.blinding();
	xm.setCurrentFitVal(real_val);
	xm.setInit(real_val);

	real_val = (double)ym + ym.blinding();
	ym.setCurrentFitVal(real_val);
	ym.setInit(real_val);

	real_val = (double)xp + xp.blinding();
	xp.setCurrentFitVal(real_val);
	xp.setInit(real_val);

	real_val = (double)yp + yp.blinding();
	yp.setCurrentFitVal(real_val);
	yp.setInit(real_val);
    }

    FitParameter  Gamma("Gamma",2,0.6629,0.0018);
    FitParameter  dGamma("dGamma",2,-0.088,0.006);
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

    /// Make full time-dependent PDF
    string marginalPdfsPrefix = "comb";
    if(fitGenMC)marginalPdfsPrefix = "Uniform";
    FullAmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma
		      , xm, ym, xp, yp
		      ,Gamma, dGamma, dm
                      ,offset_sigma_dt, scale_mean_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, marginalPdfsPrefix );

    // Simultaneous pdfs
    FullAmpsPdfFlexiFastCPV pdf_Run1_t0(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,c0_Run1_t0, c1_Run1_t0, c2_Run1_t0 ,c3_Run1_t0, c4_Run1_t0, c5_Run1_t0
                      ,c6_Run1_t0, c7_Run1_t0, c8_Run1_t0, c9_Run1_t0,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1_t0" );

    FullAmpsPdfFlexiFastCPV pdf_Run1_t1(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,c0_Run1_t1, c1_Run1_t1, c2_Run1_t1 ,c3_Run1_t1, c4_Run1_t1, c5_Run1_t1
                      ,c6_Run1_t1, c7_Run1_t1, c8_Run1_t1, c9_Run1_t1,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1_t1" );

    FullAmpsPdfFlexiFastCPV pdf_Run2_t0(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
                      ,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
                      ,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
                      p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                      avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                      p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                      avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                      production_asym_Run2, detection_asym_Run2, "Run2_t0" );

    FullAmpsPdfFlexiFastCPV pdf_Run2_t1(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
                      ,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
                      ,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
                      p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                      avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                      p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                      avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                      production_asym_Run2, detection_asym_Run2, "Run2_t1" );
    
	if(initCPcoeff){
		CPcoeffLL coeffLL(&pdf,0.741313,-0.741313,0.221461,0.53761,-0.333118,0.00910515,0.01); 
		Minimiser miniCP(&coeffLL);
		miniCP.doFit();
		miniCP.printResultVsInput();
		pdf.calculateCP_coeff();
		throw "";
	}

 	/// Read data
    	DalitzEventList eventList, eventList_f, eventList_f_bar;
    	DalitzEventList eventList_Run1_t0,eventList_Run1_t1,eventList_Run2_t0,eventList_Run2_t1;

	double t,dt;
        int f;
        double Bs_ID,Ds_ID;
        int q_OS;
// 	Short_t q_SS;
	int q_SS;
        double eta_OS;
// 	Float_t eta_SS;
	double eta_SS;
        double sw;
        int year,run,Ds_finalState,trigger;

        double K[4];
        double pip[4];
        double pim[4];
        double Ds_Kp[4],Ds_Km[4],Ds_pim[4],Ds[4];
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
		if(mode == "fit" && doToyStudy == 1)tree->Add(((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());
		else tree->Add(((string)InputFileName).c_str());
		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("N_Bs_sw",1);
		tree->SetBranchStatus("year",1);
		tree->SetBranchStatus("*DEC",1);
		tree->SetBranchStatus("*PROB",1);
		tree->SetBranchStatus("*OS*",1);
		tree->SetBranchStatus("*TAU*",1);
		tree->SetBranchStatus("*ID*",1);
		tree->SetBranchStatus("weight",1);
		tree->SetBranchStatus("Bs_DTF_MM",1);
		tree->SetBranchStatus("BsDTF_*P*",1);
		tree->SetBranchStatus("TriggerCat",1);
		tree->SetBranchStatus("run",1);
	
		tree->SetBranchAddress("Bs_BsDTF_TAU",&t);
		tree->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
		tree->SetBranchAddress("Ds_ID",&f);
// 		tree->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
// 		tree->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&eta_OS);
// 		tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
// 		tree->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&eta_SS);
		tree->SetBranchAddress("OS_Combination_DEC",&q_OS);
		tree->SetBranchAddress("OS_Combination_PROB",&eta_OS);
		tree->SetBranchAddress("SS_Kaon_DEC",&q_SS);
		tree->SetBranchAddress("SS_Kaon_PROB",&eta_SS);
		tree->SetBranchAddress("N_Bs_sw",&sw);
		tree->SetBranchAddress("year",&year);
		tree->SetBranchAddress("run",&run);
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
		
		if(doToyStudy && mode != "gen"){
			tree->SetBranchAddress("BsDTF_Ds_PX",&Ds[0]);
			tree->SetBranchAddress("BsDTF_Ds_PY",&Ds[1]);
			tree->SetBranchAddress("BsDTF_Ds_PZ",&Ds[2]);
			tree->SetBranchAddress("BsDTF_Ds_PE",&Ds[3]);    
		}
		else {
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
	}

	int N_sample = tree->GetEntries();
	if(N_bootstrap < 0)N_bootstrap = N_sample;
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
		TLorentzVector D_p;
		if(doToyStudy && mode == "fit"){
			D_p = TLorentzVector(sign*Ds[0],sign*Ds[1],sign*Ds[2],Ds[3]);
		}
		else {
			TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
			TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
			TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
			D_p = D_Kp_p + D_Km_p + D_pim_p;
		}
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
// 		if( sw > 1.4) continue;

		evt.setWeight(sw);
		evt.setValueInVector(0, t);
	        evt.setValueInVector(1, (year == 17) ? dt*1.1 : dt);   
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
    if(eventList.size()==0){
	cout << "ERROR: Have no data events !" << endl;
	throw "ERROR";
    }    


    /// Total model
    SumPdf<IDalitzEvent> tot_pdf_Run1_t0(sigfraction,pdf_Run1_t0,pdf_bkg);
    SumPdf<IDalitzEvent> tot_pdf_Run1_t1(sigfraction,pdf_Run1_t1,pdf_bkg);
    SumPdf<IDalitzEvent> tot_pdf_Run2_t0(sigfraction,pdf_Run2_t0,pdf_bkg);
    SumPdf<IDalitzEvent> tot_pdf_Run2_t1(sigfraction,pdf_Run2_t1,pdf_bkg);

    /// Likelihood
    Neg2LL neg2LL(pdf, eventList);   

    Neg2LL neg2LL_Run1_t0(tot_pdf_Run1_t0, eventList_Run1_t0);    
    Neg2LL neg2LL_Run1_t1(tot_pdf_Run1_t1, eventList_Run1_t1);    
    Neg2LL neg2LL_Run2_t0(tot_pdf_Run2_t0, eventList_Run2_t0);    
    Neg2LL neg2LL_Run2_t1(tot_pdf_Run2_t1, eventList_Run2_t1);    
 
    Neg2LLSum neg2LL_sim;
    if(eventList_Run1_t0.size()>0)neg2LL_sim.add(&neg2LL_Run1_t0);
    if(eventList_Run1_t1.size()>0)neg2LL_sim.add(&neg2LL_Run1_t1);
    if(eventList_Run2_t0.size()>0)neg2LL_sim.add(&neg2LL_Run2_t0);
    if(eventList_Run2_t1.size()>0)neg2LL_sim.add(&neg2LL_Run2_t1);

    /// Add tagging constraints
    Neg2LLMultiConstraint constrains_tagging_Run1(MinuitParameterSet::getDefaultSet(),"_Tagging_Run1");
    Neg2LLMultiConstraint constrains_tagging_Run2(MinuitParameterSet::getDefaultSet(),"_Tagging_Run2");
    if(useGaussConstrainsTagging){
 	    neg2LL_sim.add(&constrains_tagging_Run1);
 	    neg2LL_sim.add(&constrains_tagging_Run2);
    }

    /// Systematic studies
    Neg2LLMultiConstraint constrains_Acc(MinuitParameterSet::getDefaultSet(),"_Acc");
    if(doAccSystematics && mode == "fit"){
	if(useCholDec){
		if(chol_index > constrains_Acc.getNumberParams()-1){
			cout << "ERROR:: Invalid cholesky index ! "<< endl;
			throw "ERROR";
		}
		constrains_Acc.smearInputValuesChol(chol_index,(step-1) - chol_index * varPerParChol, 1);
	}
	else
	{ 
		constrains_Acc.smearInputValues();
	}
    }

    Neg2LLMultiConstraint constrains_sys(MinuitParameterSet::getDefaultSet(),("_" + (string)doSystematic).c_str());
    if((string)doSystematic != "" && mode == "fit")constrains_sys.smearInputValues();

    if(doRadiusBWSystematic && mode == "fit"){
				double val = gRandom->Uniform(0.,0.003);
				mps->getParPtr("BW_radius")->setCurrentFitVal(val);
				((FitParameter*)mps->getParPtr("BW_radius"))->setInit(val);
				cout << "Set parameter " << "BW_radius" << " to " << val << endl;  
    }

    /// LASSO
    double stepSize = 1;
    lambda = lambda + (step-1) * stepSize;
    LASSO_flexi lasso(&ampsSig,lambda);
    LASSO_flexi lasso_bar(&ampsSig_bar,lambda);
    
    Neg2LLSum neg2LL_lasso(&neg2LL,&lasso_bar);
    Neg2LLSum neg2LL_sim_lasso(&neg2LL_sim,&lasso_bar);    
    if(useLASSO == 1) { 
	neg2LL_lasso.add(&lasso);    
	neg2LL_sim_lasso.add(&lasso);
    }

    /// Fit
    Minimiser mini;
    if(useLASSO>0){
    	if(doSimFit)mini.attachFunction(&neg2LL_sim_lasso);
    	else mini.attachFunction(&neg2LL_lasso);    
    }
    else {
    	if(doSimFit)mini.attachFunction(&neg2LL_sim);
    	else mini.attachFunction(&neg2LL);    
    }
    if(mode == "fit"){
	mini.doFit();
    	mini.printResultVsInput();
// 	mini.CallMigrad();
//     	mini.printResultVsInput();
// 	mini.CallImprove();
//     	mini.printResultVsInput();
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
    TH1D* h_q_OS = new TH1D("h_q_OS",";q_{OS};Events (norm.) ",3,-1.5,1.5);
    TH1D* h_q_SS = new TH1D("h_q_SS",";q_{SS};Events (norm.) ",3,-1.5,1.5);
    TH1D* h_f = new TH1D("h_f",";q_{f};Events (norm.) ",2,-2,2);

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

        h_q_OS->Fill(q1,eventList[i].getWeight());
        h_q_SS->Fill(q2,eventList[i].getWeight());
        h_f->Fill(f_evt,eventList[i].getWeight());
        
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
        double D_tot = 1.;//(1.-2.*abs(w_eff)) * exp(-pow(t_pdf.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);

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
			h_t_mp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(w_eff<w_max)h_N_mixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
            }
	    else if(q_eff==0 && f_evt == 1)h_t_0p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == 1){
			h_t_pp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(w_eff<w_max)h_N_unmixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
	    }
	    else if(q_eff==-1 && f_evt == -1){
			h_t_mm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
        	    	if(w_eff<w_max)h_N_unmixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
	    }
	    else if(q_eff==0 && f_evt == -1)h_t_0m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == -1){
			h_t_pm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(w_eff<w_max)h_N_mixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
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

    /// Generate toys 
    DalitzEventList toys;
    if(mode == "gen"){
	if(doSimFit) {
		toys.Add(pdf_Run1_t0.generateToys(N_scale_toys * N_Run1_t0,1,0));
		toys.Add(pdf_Run1_t1.generateToys(N_scale_toys *N_Run1_t1,1,1));
		toys.Add(pdf_Run2_t0.generateToys(N_scale_toys *N_Run2_t0,2,0));
		toys.Add(pdf_Run2_t1.generateToys(N_scale_toys *N_Run2_t1,2,1));
	}
	else toys.Add(pdf.generateToys(N_scale_toys *N));

	pdf.saveEventListToFile(toys,((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());

	return;
    }

    /// Calculate pulls
    gDirectory->cd();
    string paraFileName = (string)OutputDir+"pull_"+ (string)doSystematic+ "_" + anythingToString((int)step)+".root";
    if(doAccSystematics) paraFileName = (string)OutputDir+"pullAcc_"+anythingToString((int)step)+".root";
    if(doAccSystematics && useCholDec) paraFileName = (string)OutputDir+"pullAccChol_"+anythingToString((int)step)+".root";
    
    TFile* paraFile = new TFile( paraFileName.c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;

    TTree* pull_tree = new TTree("Coherence","Coherence");
    double r_val,delta_val,gamma_val,k_val,n2ll;
    double r_err,delta_err,gamma_err,k_err;
    double sum_val,sum_bar_val,sum_err,sum_bar_err;
    double C_val,D_val,S_val;
    double Cbar_val,Dbar_val,Sbar_val;
    double xp_val,xm_val,yp_val,ym_val;
    double chi2_val,chi2_6D_val;
    double mass_K1_val,width_K1_val,mass_Ks_val,width_Ks_val;
    double mass_K1_err,width_K1_err,mass_Ks_err,width_Ks_err;

    TBranch* br_mass_K1 = pull_tree->Branch( "mass_K_1__1400_p_mean", &mass_K1_val, "mass_K1_val/D" );
    TBranch* br_mass_K1_err = pull_tree->Branch( "mass_K_1__1400_p_err", &mass_K1_err, "mass_K1_err/D" );
    TBranch* br_width_K1 = pull_tree->Branch( "width_K_1__1400_p_mean", &width_K1_val, "width_K1_val/D" );
    TBranch* br_width_K1_err = pull_tree->Branch( "width_K_1__1400_p_err", &width_K1_err, "width_K1_err/D" );

    TBranch* br_mass_Ks = pull_tree->Branch( "mass_Ks_1410_p_mean", &mass_Ks_val, "mass_Ks_val/D" );
    TBranch* br_mass_Ks_err = pull_tree->Branch( "mass_Ks_1410_p_err", &mass_Ks_err, "mass_Ks_err/D" );
    TBranch* br_width_Ks = pull_tree->Branch( "width_Ks_1410_p_mean", &width_Ks_val, "width_Ks_val/D" );
    TBranch* br_width_Ks_err = pull_tree->Branch( "width_Ks_1410_p_err", &width_Ks_err, "width_Ks_err/D" );

    TBranch* br_r = pull_tree->Branch( "r_mean", &r_val, "r_val/D" );
    TBranch* br_delta = pull_tree->Branch( "delta_mean", &delta_val, "delta_val/D" );
    TBranch* br_gamma = pull_tree->Branch( "gamma_mean", &gamma_val, "gamma_val/D" );
    TBranch* br_k = pull_tree->Branch( "k_mean", &k_val, "k_val/D" );

    TBranch* br_r_err = pull_tree->Branch( "r_err", &r_err, "r_err/D" );
    TBranch* br_delta_err = pull_tree->Branch( "delta_err", &delta_err, "delta_err/D" );
    TBranch* br_gamma_err = pull_tree->Branch( "gamma_err", &gamma_err, "gamma_err/D" );
    TBranch* br_k_err = pull_tree->Branch( "k_err", &k_err, "k_err/D" );

    TBranch* br_sum = pull_tree->Branch( "Sum_mean", &sum_val, "sum_val/D" );
    TBranch* br_sum_err = pull_tree->Branch( "Sum_err", &sum_err, "sum_err/D" );
    TBranch* br_sum_bar = pull_tree->Branch( "bar_Sum_mean", &sum_bar_val, "sum_bar_val/D" );
    TBranch* br_sum_bar_err = pull_tree->Branch( "bar_Sum_err", &sum_bar_err, "sum_bar_err/D" );

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

    r_val = r.blindedMean();
    delta_val = delta.blindedMean();
    gamma_val = gamma.blindedMean();

    r_err = r.err();
    delta_err = delta.err();
    gamma_err = gamma.err();

    xp_val = xp.blindedMean();
    yp_val = yp.blindedMean();
    xm_val = xm.blindedMean();
    ym_val = ym.blindedMean();

    for(int i=0; i < mps->size(); i++){
		if(A_is_in_B("mass_K_1__1400_p_mean",((FitParameter*)mps->getParPtr(i))->name()) ){
			 mass_K1_val = mps->getParPtr(i)->mean();
			 mass_K1_err = mps->getParPtr(i)->err();
		}
		if(A_is_in_B("width_K_1__1400_p_mean",((FitParameter*)mps->getParPtr(i))->name()) ){
			 width_K1_val = mps->getParPtr(i)->mean();
			 width_K1_err = mps->getParPtr(i)->err();
		}
		if(A_is_in_B("mass_Ks_1410_p_mean",((FitParameter*)mps->getParPtr(i))->name()) ){
			 mass_Ks_val = mps->getParPtr(i)->mean();
			 mass_Ks_err = mps->getParPtr(i)->err();
		}
		if(A_is_in_B("width_Ks_1410_p_mean",((FitParameter*)mps->getParPtr(i))->name()) ){
			 width_Ks_val = mps->getParPtr(i)->mean();
			 width_Ks_err = mps->getParPtr(i)->err();
		}
     }

    string outTableName = (string)OutputDir+"FitAmpResults_"+anythingToString((int)step);
    if(updateAnaNote)outTableName = "../../../../../TD-AnaNote/latex/tables/fullFit/"+(string)OutputDir+"fitFractions";

    string fractionFileName = (string)OutputDir+"fitFractions_"+ (string)doSystematic+ "_" + anythingToString((int)step);
    if(doAccSystematics) fractionFileName = (string)OutputDir+"fitFractionsAcc_"+anythingToString((int)step);
    if(doAccSystematics && useCholDec) fractionFileName = (string)OutputDir+"fitFractionsChol_"+anythingToString((int)step);

    vector<double> CP_coeff;
    if(mode == "fit"){
// 	pdf.doFinalStatsAndSaveForAmp12(&mini,outTableName,fractionFileName);
	CP_coeff = pdf_Run1_t0.calculateCP_coeff();
	C_val = CP_coeff[0];
	Cbar_val = CP_coeff[1];
	D_val = CP_coeff[2];
	Dbar_val = CP_coeff[3];
	S_val = CP_coeff[4];
	Sbar_val = CP_coeff[5];
	k_val = CP_coeff[6];

	if(doFractions){

		    vector<FitFractionList> ffl;
		    vector<FitFractionList> ffl_bar;

		    int N_frac_toys = (doFractionsErr ? 100 : 1);

		    RooArgList xvec, mu;
		    for(int i=0; i < mps->size(); i++){
					if(0 == mps->getParPtr(i)) continue;
					if(((FitParameter*)mps->getParPtr(i))->iFixInit() != 0)continue;
					double mean = mps->getParPtr(i)->mean();	
					double error = mps->getParPtr(i)->err();	
		
					RooRealVar* x = new RooRealVar(("x_"+anythingToString(i)).c_str(), ("x_"+anythingToString(i)).c_str(),mean-10.*error,mean+10.*error);
					xvec.add(*x);
					mu.add(RooRealConstant::value(mean));
		    }
	
		    TMatrixTSym<double> cov = mini.covMatrix();     
		    cov.Print();
		    xvec.Print();
		    mu.Print();
		    RooMultiVarGaussian gauss_cov("gauss_cov","gauss_cov",xvec, mu, cov);
		    RooDataSet* data_cov = gauss_cov.generate(xvec, N_frac_toys);

		    for(int i=0; i < N_frac_toys; i++){
	
			if(i>0){
					RooArgSet* xvec_cov= (RooArgSet*)data_cov->get(i);
	
   					for(int j=0; j < mps->size(); j++){
						if(0 == mps->getParPtr(j)) continue;
						if(((FitParameter*)mps->getParPtr(j))->iFixInit() != 0)continue;
						mps->getParPtr(j)->setCurrentFitVal(((RooRealVar*)xvec_cov->find(("x_"+anythingToString(j)).c_str()))->getVal());
					}
			}

			AmpsPdfFlexiFast ampsSig_noEff(pat, &fas, 0, integPrecision,integMethod, "SignalIntegrationEvents_toys_phspCut.root");
			AmpsPdfFlexiFast ampsSig_bar_noEff(pat, &fas_bar, 0, integPrecision,integMethod, "SignalIntegrationEvents_toys_phspCut.root");
	
			ampsSig_noEff.redoIntegrator();
			ampsSig_bar_noEff.redoIntegrator();
			if(i==0)ampsSig_noEff.doFinalStatsAndSave(&mini,(outTableName+".tex").c_str(),(fractionFileName+".root").c_str());
			if(i==0)ampsSig_bar_noEff.doFinalStatsAndSave(&mini,(outTableName+"_bar.tex").c_str(),(fractionFileName+"_Bar.root").c_str());

			ampsSig_noEff.doFractions();
			ampsSig_bar_noEff.doFractions();
			
			ffl.push_back(ampsSig_noEff.getFractions());
			ffl_bar.push_back(ampsSig_bar_noEff.getFractions());
		   }

		    vector<double> ffl_err;
		    for(int i=0; i < ffl[0].size(); i++){
			double ff_var = 0;
	
			for(int j=0; j < N_frac_toys; j++){
				for(int k=j+1; k < N_frac_toys; k++){
					ff_var += pow(ffl[j][i].frac() - ffl[k][i].frac(),2);
				}
			}
			ffl_err.push_back(sqrt(ff_var/((double)N_frac_toys*((double)N_frac_toys-.1))));
		    }

		    vector<double> ffl_sum;
		    for(int j=0; j < N_frac_toys; j++){
			    double sum = 0;
			    for(int i=0; i < ffl[0].size(); i++){
				sum += ffl[j][i].frac();	
			    }
			    ffl_sum.push_back(sum);
		    }

		   double sum_var = 0;
		   for(int j=0; j < N_frac_toys; j++){
				for(int k=j+1; k < N_frac_toys; k++){
					sum_var += pow(ffl_sum[j] - ffl_sum[k],2);
				}
		   }
		   sum_val = ffl_sum[0];
		   sum_err = sqrt(sum_var/((double)N_frac_toys*((double)N_frac_toys-.1)));

		    vector<double> ffl_bar_err;
		    for(int i=0; i < ffl_bar[0].size(); i++){
			double ff_var = 0;
	
			for(int j=0; j < N_frac_toys; j++){
				for(int k=j+1; k < N_frac_toys; k++){
					ff_var += pow(ffl_bar[j][i].frac() - ffl_bar[k][i].frac(),2);
				}
			}
			ffl_bar_err.push_back(sqrt(ff_var/((double)N_frac_toys*((double)N_frac_toys-.1))));
		    }

		    vector<double> ffl_bar_sum;
		    for(int j=0; j < N_frac_toys; j++){
			    double sum = 0;
			    for(int i=0; i < ffl_bar[0].size(); i++){
				sum += ffl_bar[j][i].frac();	
			    }
			    ffl_bar_sum.push_back(sum);
		    }

		   double sum_bar_var = 0;
		   for(int j=0; j < N_frac_toys; j++){
				for(int k=j+1; k < N_frac_toys; k++){
					sum_bar_var += pow(ffl_bar_sum[j] - ffl_bar_sum[k],2);
				}
		   }
		   sum_bar_val = ffl_bar_sum[0];
		   sum_bar_err = sqrt(sum_bar_var/((double)N_frac_toys*((double)N_frac_toys-.1)));

		    vector<double*> mean_vals;
		    vector<double*> err_vals;
		    for(int i=0; i < ffl[0].size(); i++){
				TString name = ffl[0][i].name() ;
				name.ReplaceAll("A(","");
				name.ReplaceAll("fit","");
				name.ReplaceAll("[","_");
				name.ReplaceAll("]","_");
				name.ReplaceAll("(","_");
				name.ReplaceAll(")","_");
				name.ReplaceAll("->","_");
				name.ReplaceAll("+","p");
				name.ReplaceAll("-","m");
				name.ReplaceAll("*","s");
				name.ReplaceAll(",","");
				name.ReplaceAll("GS","");
				name.ReplaceAll("Lass","");
				name.ReplaceAll("RhoOmega","");
				name.ReplaceAll("Bugg","");
				name.ReplaceAll("Flatte","");

			        double * mean = new double[1];
			        mean_vals.push_back(mean);
			        double * err = new double[1];
			        err_vals.push_back(err);

				*mean = ffl[0][i].frac();
				*err = ffl_err[i];
		
				TBranch* Bra_mean = pull_tree->Branch( ((string)name+"_mean").c_str(), mean, ((string)(name+"_mean/D")).c_str());
				TBranch* Bra_err = pull_tree->Branch( ((string)name+"_err").c_str(), err, ((string)(name+"_err/D")).c_str());
		   }

		    for(int i=0; i < ffl_bar[0].size(); i++){
				TString name = ffl_bar[0][i].name() ;
				name.ReplaceAll("A(","");
				name.ReplaceAll("fit","");
				name.ReplaceAll("[","_");
				name.ReplaceAll("]","_");
				name.ReplaceAll("(","_");
				name.ReplaceAll(")","_");
				name.ReplaceAll("->","_");
				name.ReplaceAll("+","p");
				name.ReplaceAll("-","m");
				name.ReplaceAll("*","s");
				name.ReplaceAll(",","");
				name.ReplaceAll("GS","");
				name.ReplaceAll("Lass","");
				name.ReplaceAll("RhoOmega","");
				name.ReplaceAll("Bugg","");
				name.ReplaceAll("Flatte","");

			        double * mean = new double[1];
			        mean_vals.push_back(mean);
			        double * err = new double[1];
			        err_vals.push_back(err);

				*mean = ffl_bar[0][i].frac();
				*err = ffl_bar_err[i];
		
				name = "bar_" + name;
				TBranch* Bra_mean = pull_tree->Branch( ((string)name+"_mean").c_str(), mean, ((string)(name+"_mean/D")).c_str());
				TBranch* Bra_err = pull_tree->Branch( ((string)name+"_err").c_str(), err, ((string)(name+"_err/D")).c_str());
		   }

		    vector<double> kl;
		    for(int i=0; i < N_frac_toys; i++){
	
			if(i>0){
					RooArgSet* xvec_cov= (RooArgSet*)data_cov->get(i);
	
   					for(int j=0; j < mps->size(); j++){
						if(0 == mps->getParPtr(j)) continue;
						if(((FitParameter*)mps->getParPtr(j))->iFixInit() != 0)continue;
						mps->getParPtr(j)->setCurrentFitVal(((RooRealVar*)xvec_cov->find(("x_"+anythingToString(j)).c_str()))->getVal());
					}
			}
		   
			DiskResidentEventList eventListMC("SignalIntegrationEvents_toys_phspCut.root","OPEN");
			vector<double> cf = coherenceFactor(fas, fas_bar, (double)r, (double)delta, eventListMC,200000);
			kl.push_back(cf[1]);
		   }

		   double k_var = 0;
		   for(int j=0; j < N_frac_toys; j++){
				for(int k=j+1; k < N_frac_toys; k++){
					k_var += pow(kl[j] - kl[k],2);
				}
		   }
		   k_val = kl[0];
		   k_err = sqrt(k_var/((double)N_frac_toys*((double)N_frac_toys-.1)));
	}

    }

    /// Loop over MC
    if(doPlots){

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
	TH1D* h_q_OS_fit = new TH1D("h_q_OS_fit",";q_{OS};Events (norm.) ",3,-1.5,1.5);
	TH1D* h_q_SS_fit = new TH1D("h_q_SS_fit",";q_{SS};Events (norm.) ",3,-1.5,1.5);
	TH1D* h_f_fit = new TH1D("h_f_fit",";q_{f};Events (norm.) ",2,-2,2);
	
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
	
	TH1D* h_N_mixed_p_fit_unfolded = new TH1D("h_N_mixed_p_fit_unfolded",";t/#tau;A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_N_unmixed_p_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_p_fit");
	TH1D* h_N_mixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_mixed_m_fit");
	TH1D* h_N_unmixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_m_fit");

	double N_OS_MC = 0;
	double N_SS_MC = 0;
	double N_OS_SS_MC = 0;
	double N_MC = 0;
	
	double w_OS_MC = 0;
	double w_SS_MC = 0;
	double w_OS_SS_MC = 0;
		
	double D_OS_MC = 0;
	double D_SS_MC = 0;
	double D_OS_SS_MC = 0;
	double D_comb_MC = 0;
	
	double N_OS_all_MC = 0;
	double N_SS_all_MC = 0;
	double w_OS_all_MC = 0;
	double w_SS_all_MC = 0;
	double D_OS_all_MC = 0;
	double D_SS_all_MC = 0;
	
	double N_Run1_t0_MC = 0;
	double N_Run1_t1_MC = 0;
	double N_Run2_t0_MC = 0;
	double N_Run2_t1_MC = 0;

	DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
	DiskResidentEventList eventListMC_rw(pat,("dummy_"+anythingToString(step)+".root").c_str(),"RECREATE");
	
	///Dalitz plots 
	for(int i = 0; i < eventListMC.size(); i++){
	
			DalitzEvent evt(eventListMC.getEvent(i));
			int f = gRandom->Uniform()>0.5 ? 1 : -1;
			if(f==-1)evt.CP_conjugateYourself();
			evt.setValueInVector(2, f);

			double pdfVal = 0;
			if(doSimFit) {
				pdfVal += pdf_Run1_t0.getVal_timeIntegrated(evt) * N_Run1_t0/N;
				pdfVal += pdf_Run1_t1.getVal_timeIntegrated(evt) * N_Run1_t1/N;
				pdfVal += pdf_Run2_t0.getVal_timeIntegrated(evt) * N_Run2_t0/N;
				pdfVal += pdf_Run2_t1.getVal_timeIntegrated(evt) * N_Run2_t1/N;
			}
			else pdfVal = pdf.getVal_timeIntegrated(evt);
			
			double weight = pdfVal*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();			
			
			double weight_A,weight_Abar;
			
			if(f==1){ 
				weight_A = fas.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
				weight_Abar = fas_bar.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			}
			else {		
				weight_A = fas_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
				weight_Abar = fas_bar_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
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
			
			evt.setWeight(weight);
			eventListMC_rw.Add(evt);
	}

	/// Time plots	
	FitParameter  C("C",1,CP_coeff[0],0.1);
	FitParameter  D("D",1,CP_coeff[2],0.1);
	FitParameter  D_bar("D_bar",1,CP_coeff[3],0.1);
	FitParameter  S("S",1,CP_coeff[4],0.1);
	FitParameter  S_bar("S_bar",1,CP_coeff[5],0.1);
	FitParameter  k("k",1,1.,0.1);
	
	FullTimePdf t_pdf(C, D, D_bar, S, S_bar, k,
			Gamma, dGamma, dm
			,offset_sigma_dt, scale_mean_dt, scale_sigma_dt, scale_sigma_2_dt
			,c0, c1, c2 ,c3, c4, c5
			,c6, c7, c8, c9,
			p0_os, p1_os, delta_p0_os, delta_p1_os, 
			avg_eta_os, tageff_os, tageff_asym_os, 
			p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
			avg_eta_ss, tageff_ss, tageff_asym_ss, 
			production_asym, detection_asym, marginalPdfsPrefix );
	
	FullTimePdf t_pdf_Run1_t0(C, D, D_bar, S, S_bar, k,
			Gamma, dGamma, dm
			,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
			,c0_Run1_t0, c1_Run1_t0, c2_Run1_t0 ,c3_Run1_t0, c4_Run1_t0, c5_Run1_t0
			,c6_Run1_t0, c7_Run1_t0, c8_Run1_t0, c9_Run1_t0,
			p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
			avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
			p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
			avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
			production_asym_Run1, detection_asym_Run1, "Run1_t0" );
	
	FullTimePdf t_pdf_Run1_t1(C, D, D_bar, S, S_bar, k,
				Gamma, dGamma, dm
				,offset_sigma_dt_Run1, scale_mean_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
				,c0_Run1_t1, c1_Run1_t1, c2_Run1_t1 ,c3_Run1_t1, c4_Run1_t1, c5_Run1_t1
				,c6_Run1_t1, c7_Run1_t1, c8_Run1_t1, c9_Run1_t1,
				p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
				avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
				p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
				avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
				production_asym_Run1, detection_asym_Run1, "Run1_t1" );
	
	FullTimePdf t_pdf_Run2_t0(C, D, D_bar, S, S_bar, k,
				Gamma, dGamma, dm
				,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
				,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
				,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
				p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
				avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
				p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
				avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
				production_asym_Run2, detection_asym_Run2, "Run2_t0" );
	
	FullTimePdf t_pdf_Run2_t1(C, D, D_bar, S, S_bar, k,
				Gamma, dGamma, dm
				,offset_sigma_dt_Run2, scale_mean_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
				,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
				,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
				p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
				avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
				p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
				avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
				production_asym_Run2, detection_asym_Run2, "Run2_t1" );


	for(int n = 0; n < 5; n++){   /// Multiple iterations needed to release memory 
		int N_sample = 250000;
		DalitzEventList sampleEvents;
   		if(doSimFit) {
			sampleEvents.Add(t_pdf_Run1_t0.generateToys(N_sample * N_Run1_t0/N,1,0));
			sampleEvents.Add(t_pdf_Run1_t1.generateToys(N_sample *N_Run1_t1/N,1,1));
			sampleEvents.Add(t_pdf_Run2_t0.generateToys(N_sample *N_Run2_t0/N,2,0));
			sampleEvents.Add(t_pdf_Run2_t1.generateToys(N_sample *N_Run2_t1/N,2,1));
		}
		else sampleEvents.Add(t_pdf_Run1_t0.generateToys(N_sample));	

		for(int i = 0; i < sampleEvents.size(); i++){

			DalitzEvent evt = sampleEvents[i];
			double t_MC = evt.getValueFromVector(0) ;
			double dt_MC = evt.getValueFromVector(1) ;
			int f_MC = evt.getValueFromVector(2) ;
			int q_OS_MC = evt.getValueFromVector(3) ;
			double eta_OS_MC = evt.getValueFromVector(4) ;
			int q_SS_MC = evt.getValueFromVector(5);
			double eta_SS_MC = evt.getValueFromVector(6);
			int run_MC = evt.getValueFromVector(7) ;
			int trigger_MC = evt.getValueFromVector(8) ;

			double weight = 1; //t_pdf_Run1_t0.getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace();
			N_MC += weight;
		
			h_t_fit->Fill(t_MC,weight);
			h_dt_fit->Fill(dt_MC,weight);
			if(q_OS_MC != 0)h_eta_OS_fit->Fill(eta_OS_MC,weight);
			if(q_SS_MC != 0)h_eta_SS_fit->Fill(eta_SS_MC,weight);
			
			int f_evt = f_MC;
			int q1 = q_OS_MC;
			int q2 = q_SS_MC;   
			int q_eff = 0;
			double w_eff = 0.5;
		
			h_q_OS_fit->Fill(q1,weight);
			h_q_SS_fit->Fill(q2,weight);
			h_f_fit->Fill(f_evt,weight);
			
			std::pair<double, double> calibrated_mistag_os;
			std::pair<double, double> calibrated_mistag_ss;
			if(doSimFit){
				if(run_MC==1){
					calibrated_mistag_os = t_pdf_Run1_t0.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run1_t0.getCalibratedMistag_SS(eta_SS_MC);
				}
				else{
					calibrated_mistag_os = t_pdf_Run2_t0.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run2_t0.getCalibratedMistag_SS(eta_SS_MC);                
				}
			}
			else{
				calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eta_OS_MC);
				calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eta_SS_MC);        
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
					N_OS_SS_MC += weight;
					w_OS_SS_MC += w_eff * weight;
					D_OS_SS_MC += pow(1.-2.*w_eff,2)* weight;
				}
			}
			else if( q1 != 0){
				//q_eff = q1;  flip tag ???
				N_OS_MC +=weight;
				w_OS_MC += w_eff * weight; 
				D_OS_MC += pow(1.-2.*w_eff,2)* weight;
			}
			else if( q2 != 0){
				//q_eff = q2;
				N_SS_MC += weight;
				w_SS_MC += w_eff * weight; 
				D_SS_MC += pow(1.-2.*w_eff,2)* weight; 
			} 
		
			D_comb_MC += pow(1.-2.*w_eff,2)* weight;
			
			if(q1 != 0) N_OS_all_MC += weight;
			if(q2 != 0)N_SS_all_MC += weight;
			
			if(q1>0){
					w_OS_all_MC +=  calibrated_mistag_os.first * weight;
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.first,2)* weight;
			} 
			else if(q1<0){
					w_OS_all_MC +=  calibrated_mistag_os.second * weight;	
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.second,2)* weight;
			}
		
			if(q2>0){
					w_SS_all_MC +=  calibrated_mistag_ss.first * weight;
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.first,2)* weight;
			} 
			else if(q2<0){
					w_SS_all_MC +=  calibrated_mistag_ss.second * weight;	
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.second,2)* weight;
			}
		
			if((string)channel=="signal"){
			
				if(q_eff==-1 && f_evt == 1){
					h_t_fit_mp->Fill(t_MC,weight);
					if(w_eff<w_max){
						h_N_mixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
						h_N_mixed_p_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight);
					}
				}
				else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight);
				else if(q_eff==1 && f_evt == 1){
					h_t_fit_pp->Fill(t_MC,weight);
					if(w_eff<w_max){
						h_N_unmixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
						h_N_unmixed_p_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight);
					}
				}
				else if(q_eff==-1 && f_evt == -1){
					h_t_fit_mm->Fill(t_MC,weight);
					if(w_eff<w_max){
						h_N_unmixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
						h_N_unmixed_m_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight);
					}	
				}
				else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(t_MC,weight);
				else if(q_eff==1 && f_evt == -1){
					h_t_fit_pm->Fill(t_MC,weight);
					if(w_eff<w_max){
						h_N_mixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
						h_N_mixed_m_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight);
					}
				}
			}
			else {   
				if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
				else if(q_eff*f_evt > 0  ){
					h_t_mixed_fit->Fill(t_MC,weight);
					if(w_eff<w_max)h_N_mixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
				}
				else{ 
					h_t_unmixed_fit->Fill(t_MC,weight);
					if(w_eff<w_max)h_N_unmixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
				}
			}
		}
	}
	
	cout << "N_MC = " << N_MC << endl << endl;
	cout << "N_OS = " << N_OS_all_MC << endl;
	cout << "N_SS = " << N_SS_all_MC << endl << endl;
	
	cout << "eff_OS = " <<(N_OS_all_MC)/N_MC << endl;
	cout << "eff_SS = " <<(N_SS_all_MC)/N_MC << endl;
		
	cout << "Tagging perfromance " << endl << endl;        
	cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
	
	cout << "OS  | " << (N_OS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_OS_all_MC)/(N_OS_MC+N_OS_SS_MC) << " | " << D_OS_all_MC/N_MC << endl;
	cout << "SS  | " << (N_SS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_SS_all_MC)/(N_SS_MC+N_OS_SS_MC) << " | " << D_SS_all_MC/N_MC << endl << endl;
	
	cout << "OS only  | " << N_OS_MC/N_MC << " | " <<  w_OS_MC/N_OS_MC << " | " << N_OS_MC/N_MC * D_OS_MC/N_OS_MC << endl;
	cout << "SS only  | " << N_SS_MC/N_MC << " | " <<  w_SS_MC/N_SS_MC << " | " << N_SS_MC/N_MC * D_SS_MC/N_SS_MC << endl;
	cout << "OS+SS    | " << N_OS_SS_MC/N_MC << " | " <<  w_OS_SS_MC/N_OS_SS_MC << " | " << N_OS_SS_MC/N_MC * D_OS_SS_MC/N_OS_SS_MC << endl;
	cout << "Combined | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC << " | "<<  (w_OS_MC+w_SS_MC+w_OS_SS_MC)/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << " | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC * D_comb_MC/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << endl << endl ;

	/// Plot
	TCanvas* c= new TCanvas("");
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

	h_q_OS->SetMinimum(0);        
	h_q_OS->SetLineColor(kBlack);
	h_q_OS->DrawNormalized("e1",1);
	h_q_OS_fit->SetLineColor(kBlue);
	h_q_OS_fit->SetLineWidth(3);
	h_q_OS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_q_OS.pdf").c_str());
	
	h_q_SS->SetMinimum(0);        
	h_q_SS->SetLineColor(kBlack);
	h_q_SS->DrawNormalized("e1",1);
	h_q_SS_fit->SetLineColor(kBlue);
	h_q_SS_fit->SetLineWidth(3);
	h_q_SS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_q_SS.pdf").c_str());
	
	h_f->SetMinimum(0);        
	h_f->SetLineColor(kBlack);
	h_f->DrawNormalized("e1",1);
	h_f_fit->SetLineColor(kBlue);
	h_f_fit->SetLineWidth(3);
	h_f_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_f.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_f.pdf").c_str());
	
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
		
		double max_asym = max(max(h_asym_p->GetMaximum(),h_asym_m->GetMaximum()),fabs(min(h_asym_p->GetMinimum(),h_asym_m->GetMinimum()))) *1.25;
		h_asym_p->SetMaximum(max_asym);
		h_asym_p->SetMinimum(-max_asym);
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

	//r_val = (double)r ;	
	double norm_A = 1./(1.+r_val*r_val); 
// 		m_Kpipi_fit_A->Integral()/(m_Kpipi_fit_A->Integral()+m_Kpipi_fit_Abar->Integral());
	double norm_Abar = r_val*r_val/(1.+r_val*r_val); 		 
// 		m_Kpipi_fit_Abar->Integral()/(m_Kpipi_fit_A->Integral()+m_Kpipi_fit_Abar->Integral());
	cout << "ratio = " << sqrt(norm_Abar/norm_A) << endl;

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

	m_Kpi->SetMinimum(0.01);
	m_Kpi->SetLineColor(kBlack);
	m_Kpi->DrawNormalized("e1",1);
	m_Kpi_fit->SetLineColor(kBlue);
	m_Kpi_fit->SetLineWidth(3);
	m_Kpi_fit->DrawNormalized("histcsame",1);
	m_Kpi_fit_A->SetLineColor(kRed+1);
	m_Kpi_fit_A->SetLineWidth(2);
	m_Kpi_fit_A->SetFillColor(kRed+1);
	m_Kpi_fit_A->SetFillStyle(3353);
	m_Kpi_fit_A->DrawNormalized("histcsame",norm_A);
	m_Kpi_fit_Abar->SetLineColor(kGreen+3);
	m_Kpi_fit_Abar->SetLineWidth(2);
	m_Kpi_fit_Abar->SetFillColor(kGreen+3);
	m_Kpi_fit_Abar->SetFillStyle(3353);
	m_Kpi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_Kpi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_Kpi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Kpi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_pipi->SetMinimum(0.01);
	m_pipi->SetLineColor(kBlack);
	m_pipi->DrawNormalized("e1",1);
	m_pipi_fit->SetLineColor(kBlue);
	m_pipi_fit->SetLineWidth(3);
	m_pipi_fit->DrawNormalized("histcsame",1);
	m_pipi_fit_A->SetLineColor(kRed+1);
	m_pipi_fit_A->SetLineWidth(2);
	m_pipi_fit_A->SetFillColor(kRed+1);
	m_pipi_fit_A->SetFillStyle(3353);
	m_pipi_fit_A->DrawNormalized("histcsame",norm_A);
	m_pipi_fit_Abar->SetLineColor(kGreen+3);
	m_pipi_fit_Abar->SetLineWidth(2);
	m_pipi_fit_Abar->SetFillColor(kGreen+3);
	m_pipi_fit_Abar->SetFillStyle(3353);
	m_pipi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_pipi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_pipi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_pipi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_pipi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_Dspipi->SetMinimum(0.01);
	m_Dspipi->SetLineColor(kBlack);
	m_Dspipi->DrawNormalized("e1",1);
	m_Dspipi_fit->SetLineColor(kBlue);
	m_Dspipi_fit->SetLineWidth(3);
	m_Dspipi_fit->DrawNormalized("histcsame",1);
	m_Dspipi_fit_A->SetLineColor(kRed+1);
	m_Dspipi_fit_A->SetLineWidth(2);
	m_Dspipi_fit_A->SetFillColor(kRed+1);
	m_Dspipi_fit_A->SetFillStyle(3353);
	m_Dspipi_fit_A->DrawNormalized("histcsame",norm_A);
	m_Dspipi_fit_Abar->SetLineColor(kGreen+3);
	m_Dspipi_fit_Abar->SetLineWidth(2);
	m_Dspipi_fit_Abar->SetFillColor(kGreen+3);
	m_Dspipi_fit_Abar->SetFillStyle(3353);
	m_Dspipi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_Dspipi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspipi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_Dspipi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspipi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_Dspi->SetMinimum(0.01);
	m_Dspi->SetLineColor(kBlack);
	m_Dspi->DrawNormalized("e1",1);
	m_Dspi_fit->SetLineColor(kBlue);
	m_Dspi_fit->SetLineWidth(3);
	m_Dspi_fit->DrawNormalized("histcsame",1);
	m_Dspi_fit_A->SetLineColor(kRed+1);
	m_Dspi_fit_A->SetLineWidth(2);
	m_Dspi_fit_A->SetFillColor(kRed+1);
	m_Dspi_fit_A->SetFillStyle(3353);
	m_Dspi_fit_A->DrawNormalized("histcsame",norm_A);
	m_Dspi_fit_Abar->SetLineColor(kGreen+3);
	m_Dspi_fit_Abar->SetLineWidth(2);
	m_Dspi_fit_Abar->SetFillColor(kGreen+3);
	m_Dspi_fit_Abar->SetFillStyle(3353);
	m_Dspi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_Dspi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_Dspi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_DsK->SetMinimum(0.01);
	m_DsK->SetLineColor(kBlack);
	m_DsK->DrawNormalized("e1",1);
	m_DsK_fit->SetLineColor(kBlue);
	m_DsK_fit->SetLineWidth(3);
	m_DsK_fit->DrawNormalized("histcsame",1);
	m_DsK_fit_A->SetLineColor(kRed+1);
	m_DsK_fit_A->SetLineWidth(2);
	m_DsK_fit_A->SetFillColor(kRed+1);
	m_DsK_fit_A->SetFillStyle(3353);
	m_DsK_fit_A->DrawNormalized("histcsame",norm_A);
	m_DsK_fit_Abar->SetLineColor(kGreen+3);
	m_DsK_fit_Abar->SetLineWidth(2);
	m_DsK_fit_Abar->SetFillColor(kGreen+3);
	m_DsK_fit_Abar->SetFillStyle(3353);
	m_DsK_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_DsK_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsK_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_DsK_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsK_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_DsKpi->SetMinimum(0.01);
	m_DsKpi->SetLineColor(kBlack);
	m_DsKpi->DrawNormalized("e1",1);
	m_DsKpi_fit->SetLineColor(kBlue);
	m_DsKpi_fit->SetLineWidth(3);
	m_DsKpi_fit->DrawNormalized("histcsame",1);
	m_DsKpi_fit_A->SetLineColor(kRed+1);
	m_DsKpi_fit_A->SetLineWidth(2);
	m_DsKpi_fit_A->SetFillColor(kRed+1);
	m_DsKpi_fit_A->SetFillStyle(3353);
	m_DsKpi_fit_A->DrawNormalized("histcsame",norm_A);
	m_DsKpi_fit_Abar->SetLineColor(kGreen+3);
	m_DsKpi_fit_Abar->SetLineWidth(2);
	m_DsKpi_fit_Abar->SetFillColor(kGreen+3);
	m_DsKpi_fit_Abar->SetFillStyle(3353);
	m_DsKpi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_DsKpi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsKpi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_DsKpi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_DsKpi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	m_Dspim->SetMinimum(0.01);
	m_Dspim->SetLineColor(kBlack);
	m_Dspim->DrawNormalized("e1",1);
	m_Dspim_fit->SetLineColor(kBlue);
	m_Dspim_fit->SetLineWidth(3);
	m_Dspim_fit->DrawNormalized("histcsame",1);
	m_Dspim_fit_A->SetLineColor(kRed+1);
	m_Dspim_fit_A->SetLineWidth(2);
	m_Dspim_fit_A->SetFillColor(kRed+1);
	m_Dspim_fit_A->SetFillStyle(3353);
	m_Dspim_fit_A->DrawNormalized("histcsame",norm_A);
	m_Dspim_fit_Abar->SetLineColor(kGreen+3);
	m_Dspim_fit_Abar->SetLineWidth(2);
	m_Dspim_fit_Abar->SetFillColor(kGreen+3);
	m_Dspim_fit_Abar->SetFillStyle(3353);
	m_Dspim_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"m_Dspim_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspim_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"m_Dspim_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"m_Dspim_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	h_cosTheta_Dspi->SetMinimum(0.01);
	h_cosTheta_Dspi->SetLineColor(kBlack);
	h_cosTheta_Dspi->DrawNormalized("e1",1);
	h_cosTheta_Dspi_fit->SetLineColor(kBlue);
	h_cosTheta_Dspi_fit->SetLineWidth(3);
	h_cosTheta_Dspi_fit->DrawNormalized("histcsame",1);
	h_cosTheta_Dspi_fit_A->SetLineColor(kRed+1);
	h_cosTheta_Dspi_fit_A->SetLineWidth(2);
	h_cosTheta_Dspi_fit_A->SetFillColor(kRed+1);
	h_cosTheta_Dspi_fit_A->SetFillStyle(3353);
	h_cosTheta_Dspi_fit_A->DrawNormalized("histcsame",norm_A);
	h_cosTheta_Dspi_fit_Abar->SetLineColor(kGreen+3);
	h_cosTheta_Dspi_fit_Abar->SetLineWidth(2);
	h_cosTheta_Dspi_fit_Abar->SetFillColor(kGreen+3);
	h_cosTheta_Dspi_fit_Abar->SetFillStyle(3353);
	h_cosTheta_Dspi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"h_cosTheta_Dspi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Dspi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_cosTheta_Dspi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Dspi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	h_cosTheta_Kpi->SetMinimum(0.01);
	h_cosTheta_Kpi->SetLineColor(kBlack);
	h_cosTheta_Kpi->DrawNormalized("e1",1);
	h_cosTheta_Kpi_fit->SetLineColor(kBlue);
	h_cosTheta_Kpi_fit->SetLineWidth(3);
	h_cosTheta_Kpi_fit->DrawNormalized("histcsame",1);
	h_cosTheta_Kpi_fit_A->SetLineColor(kRed+1);
	h_cosTheta_Kpi_fit_A->SetLineWidth(2);
	h_cosTheta_Kpi_fit_A->SetFillColor(kRed+1);
	h_cosTheta_Kpi_fit_A->SetFillStyle(3353);
	h_cosTheta_Kpi_fit_A->DrawNormalized("histcsame",norm_A);
	h_cosTheta_Kpi_fit_Abar->SetLineColor(kGreen+3);
	h_cosTheta_Kpi_fit_Abar->SetLineWidth(2);
	h_cosTheta_Kpi_fit_Abar->SetFillColor(kGreen+3);
	h_cosTheta_Kpi_fit_Abar->SetFillStyle(3353);
	h_cosTheta_Kpi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"h_cosTheta_Kpi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Kpi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_cosTheta_Kpi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_cosTheta_Kpi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	h_phi_Kpi_Dspi->SetMinimum(0.01);
	h_phi_Kpi_Dspi->SetLineColor(kBlack);
	h_phi_Kpi_Dspi->DrawNormalized("e1",1);
	h_phi_Kpi_Dspi_fit->SetLineColor(kBlue);
	h_phi_Kpi_Dspi_fit->SetLineWidth(3);
	h_phi_Kpi_Dspi_fit->DrawNormalized("histcsame",1);
	h_phi_Kpi_Dspi_fit_A->SetLineColor(kRed+1);
	h_phi_Kpi_Dspi_fit_A->SetLineWidth(2);
	h_phi_Kpi_Dspi_fit_A->SetFillColor(kRed+1);
	h_phi_Kpi_Dspi_fit_A->SetFillStyle(3353);
	h_phi_Kpi_Dspi_fit_A->DrawNormalized("histcsame",norm_A);
	h_phi_Kpi_Dspi_fit_Abar->SetLineColor(kGreen+3);
	h_phi_Kpi_Dspi_fit_Abar->SetLineWidth(2);
	h_phi_Kpi_Dspi_fit_Abar->SetFillColor(kGreen+3);
	h_phi_Kpi_Dspi_fit_Abar->SetFillStyle(3353);
	h_phi_Kpi_Dspi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
	c->Print(((string)OutputDir+"h_phi_Kpi_Dspi_mod.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_phi_Kpi_Dspi_mod.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_phi_Kpi_Dspi_mod_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_phi_Kpi_Dspi_mod_log.pdf").c_str());
	gPad->SetLogy(0);

	chi2_val = getChi2(eventList,eventListMC_rw);
	//chi2_6D_val = getChi2_6D(eventList,eventListMC_rw);
   }

   if(useLASSO>0){
		paraFile->cd();

		double N_sig = 0;
       		for (int i=0; i<eventList.size(); i++) N_sig += eventList[i].getWeight();

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
			if(doSimFit)y[0]=neg2LL_sim.getVal(); 
			else y[0]=neg2LL.getVal();
			if(useLASSO==1)y[0]+= 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
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
			if(doSimFit)y[0]=neg2LL_sim.getVal(); 
			else y[0]=neg2LL.getVal();
			if(useLASSO==1)y[0]+= lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
			y[0]+= lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
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
	}

	// fill tree
	paraFile->cd();
	pull_tree->Fill();
	pull_tree->SetDirectory(paraFile);
	pull_tree->Write();
	paraFile->Close();
	delete paraFile;

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

void animate(int step=0){
	TRandom3 ranLux;
	NamedParameter<int> RandomSeed("RandomSeed", 0);
	ranLux.SetSeed((int)RandomSeed);
	gRandom = &ranLux;
	
	FitAmplitude::AutogenerateFitFile();
	NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
	NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
	DalitzEventPattern pat(EventPattern.getVector());
	DalitzEventPattern pat_CP = pat.makeCPConjugate();

	NamedParameter<string> IntegratorEventFile("IntegratorEventFile" , (std::string) "SignalIntegrationEvents.root" , (char*) 0);
        TString integratorEventFile = (string) IntegratorEventFile;
        TString integratorEventFile_CP = (string) IntegratorEventFile;
        integratorEventFile_CP.ReplaceAll(".root","_CP.root");
        NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
        NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");
	DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
	
        NamedParameter<string> addAmpName("addAmpName", (std::string) "", (char*) 0);

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
		
	DalitzEventList eventListPhsp,eventListPhsp_CP;
	eventListPhsp.generatePhaseSpaceEvents(200,pat);
	eventListPhsp_CP.generatePhaseSpaceEvents(200,pat_CP);

	/// Define amplitude model
	FitAmpSum fas_tmp((DalitzEventPattern)pat);
	
	/// Normalize amps
	{
		DalitzEventList eventListNorm;
		TFile *file =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
		TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
		eventListNorm.fromNtuple(tree,0.2);
		fas_tmp.normalizeAmps(eventListNorm);
	}
	
	MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
	
    ///Choose reference amp
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1270)+");
    FitAmpSum fas(*List_1);
    FitAmpSum fas_bar(*List_1);
    
    /// Define relative decay modes    
    // A
    FitParameter a_K1_1400_Amp("a_K1_1400_Amp",1,1,0.01);
    FitParameter a_K1_1400_Phase("a_K1_1400_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_K1_1400 = new AmpRatio(a_K1_1400_Amp,a_K1_1400_Phase);

    FitParameter a_Ks_1410_Amp("a_Ks_1410_Amp",1,1,0.01);
    FitParameter a_Ks_1410_Phase("a_Ks_1410_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_Ks_1410 = new AmpRatio(a_Ks_1410_Amp,a_Ks_1410_Phase);

    FitParameter a_K_1460_Amp("a_K_1460_Amp",1,1,0.01);
    FitParameter a_K_1460_Phase("a_K_1460_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_K_1460 = new AmpRatio(a_K_1460_Amp,a_K_1460_Phase);

    FitParameter a_NS_Ks_Amp("a_NS_Ks_Amp",1,1,0.01);
    FitParameter a_NS_Ks_Phase("a_NS_Ks_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_Ks = new AmpRatio(a_NS_Ks_Amp,a_NS_Ks_Phase);

    FitParameter a_NS_sigma_Amp("a_NS_sigma_Amp",1,1,0.01);
    FitParameter a_NS_sigma_Phase("a_NS_sigma_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_sigma = new AmpRatio(a_NS_sigma_Amp,a_NS_sigma_Phase);

    FitParameter a_NS_rho_Amp("a_NS_rho_Amp",1,1,0.01);
    FitParameter a_NS_rho_Phase("a_NS_rho_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_NS_rho = new AmpRatio(a_NS_rho_Amp,a_NS_rho_Phase);

    FitParameter a_sys_Amp("a_sys_Amp",1,1,0.01);
    FitParameter a_sys_Phase("a_sys_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> a_sys = new AmpRatio(a_sys_Amp,a_sys_Phase);

    // Abar
    FitParameter abar_K1_1400_Amp("abar_K1_1400_Amp",1,1,0.01);
    FitParameter abar_K1_1400_Phase("abar_K1_1400_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_K1_1400 = new AmpRatio(abar_K1_1400_Amp,abar_K1_1400_Phase);

    FitParameter abar_Ks_1410_Amp("abar_Ks_1410_Amp",1,1,0.01);
    FitParameter abar_Ks_1410_Phase("abar_Ks_1410_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_Ks_1410 = new AmpRatio(abar_Ks_1410_Amp,abar_Ks_1410_Phase);

    FitParameter abar_K_1460_Amp("abar_K_1460_Amp",1,1,0.01);
    FitParameter abar_K_1460_Phase("abar_K_1460_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_K_1460 = new AmpRatio(abar_K_1460_Amp,abar_K_1460_Phase);

    FitParameter abar_NS_Ks_Amp("abar_NS_Ks_Amp",1,1,0.01);
    FitParameter abar_NS_Ks_Phase("abar_NS_Ks_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_Ks = new AmpRatio(abar_NS_Ks_Amp,abar_NS_Ks_Phase);

    FitParameter abar_NS_sigma_Amp("abar_NS_sigma_Amp",1,1,0.01);
    FitParameter abar_NS_sigma_Phase("abar_NS_sigma_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_sigma = new AmpRatio(abar_NS_sigma_Amp,abar_NS_sigma_Phase);

    FitParameter abar_NS_rho_Amp("abar_NS_rho_Amp",1,1,0.01);
    FitParameter abar_NS_rho_Phase("abar_NS_rho_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_NS_rho = new AmpRatio(abar_NS_rho_Amp,abar_NS_rho_Phase);

    FitParameter abar_sys_Amp("abar_sys_Amp",1,1,0.01);
    FitParameter abar_sys_Phase("abar_sys_Phase",1,0,0.01); 
    counted_ptr<IReturnComplex> abar_sys = new AmpRatio(abar_sys_Amp,abar_sys_Phase);


    /// Add amps to A and Abar
    AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K(1)(1400)+",a_K1_1400,abar_K1_1400);
    
    if(a_Ks_1410_Amp.iFixInit() != 1 && abar_Ks_1410_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K*(1410)+",a_Ks_1410,abar_Ks_1410);
    else if(a_Ks_1410_Amp.iFixInit() != 1 && abar_Ks_1410_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"K*(1410)+",a_Ks_1410);
    else if(a_Ks_1410_Amp.iFixInit() == 1 && abar_Ks_1410_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"K*(1410)+",abar_Ks_1410);
    
    if(a_K_1460_Amp.iFixInit() != 1 && abar_K_1460_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"K(1460)+",a_K_1460,abar_K_1460);
    else if(a_K_1460_Amp.iFixInit() != 1 && abar_K_1460_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"K(1460)+",a_K_1460);
    else if(a_K_1460_Amp.iFixInit() == 1 && abar_K_1460_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp,fas_bar,"K(1460)+",abar_K_1460);
    
    if(a_NS_Ks_Amp.iFixInit() != 1 && abar_NS_Ks_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"(->Ds-,pi+),K*(892)0(->K+,pi-)",a_NS_Ks,abar_NS_Ks);
    else if(a_NS_Ks_Amp.iFixInit() != 1 && abar_NS_Ks_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"(->Ds-,pi+),K*(892)0(->K+,pi-)",a_NS_Ks);
    else if(a_NS_Ks_Amp.iFixInit() == 1 && abar_NS_Ks_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"(->Ds-,pi+),K*(892)0(->K+,pi-)",abar_NS_Ks);
    
    if(a_NS_rho_Amp.iFixInit() != 1 && abar_NS_rho_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,"(->Ds-,K+),rho(770)0(->pi+,pi-)",a_NS_rho,abar_NS_rho);
    else if(a_NS_rho_Amp.iFixInit() != 1 && abar_NS_rho_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,"(->Ds-,K+),rho(770)0(->pi+,pi-)",a_NS_rho);
    else if(a_NS_rho_Amp.iFixInit() == 1 && abar_NS_rho_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,"(->Ds-,K+),rho(770)0(->pi+,pi-)",abar_NS_rho);

    if(a_sys_Amp.iFixInit() != 1 && abar_sys_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas, fas_bar,(string) addAmpName,a_sys,abar_sys);
    else if(a_sys_Amp.iFixInit() != 1 && abar_sys_Amp.iFixInit() == 1)AddScaledAmpsToList(fas_tmp, fas,(string) addAmpName,a_sys);
    else if(a_sys_Amp.iFixInit() == 1 && abar_sys_Amp.iFixInit() != 1)AddScaledAmpsToList(fas_tmp, fas_bar,(string) addAmpName,abar_sys);
   
    fas.print();
    fas_bar.print();

    /// Define B -> f amplitude        
    fas.setTag(1);
    /// Define Bbar -> f amplitude
    fas_bar.setTag(-1);

    /// CP conjugate amplitudes
    FitAmpSum fas_CP(fas);
    fas_CP.CPConjugateSameFitParameters();

    FitAmpSum fas_bar_CP(fas_bar);
    fas_bar_CP.CPConjugateSameFitParameters();

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
	
    FitParameter xm("xm",1,0,0.01);
    FitParameter ym("ym",1,0,0.01); 
    FitParameter xp("xp",1,0,0.01);
    FitParameter yp("yp",1,0,0.01); 


    FitParameter  Gamma("Gamma",2,0.6629,0.0018);
    FitParameter  dGamma("dGamma",2,-0.088,0.006);
    FitParameter  dm("dm",2,17.757,0.1);

	FitParameter  scale_mean_dt("scale_mean_dt",1,1,0.1);
	FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
	FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.,0.1);
	FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
	FitParameter  p0_os("p0_os",1,0.,0.);
	FitParameter  p1_os("p1_os",1,1.,0.);
	FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
	FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
	FitParameter  avg_eta_os("avg_eta_os",1,0.,0.);
	FitParameter  tageff_os("tageff_os",1,1.,0.);
	FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
	FitParameter  p0_ss("p0_ss",1,0.,0.);
	FitParameter  p1_ss("p1_ss",1,1.,0.);
	FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
	FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
	FitParameter  avg_eta_ss("avg_eta_ss",1,0.,0.);
	FitParameter  tageff_ss("tageff_ss",1,1.,0.);
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

        TH1D* h_t = new TH1D("h_t",";t",50,0,2*pi/dm);
        TH1D* m_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,1,2);
	TH1D* m_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.6,1.2);
	TH1D* m_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.2,1.2);
	TH1D* m_Dspipi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,2,5.5);
	TH1D* m_Dspi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});",50,1.5,5);

        TH1D* m_Kpipi_CP = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,1,2);
	TH1D* m_Kpi_CP = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.6,1.2);
	TH1D* m_pipi_CP = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.2,1.2);
	TH1D* m_Dspipi_CP = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,2,5.5);
	TH1D* m_Dspi_CP = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});",50,1.5,5);

        TH1D* m_Kpipi_A = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,1,2);
	TH1D* m_Kpi_A = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.6,1.2);
	TH1D* m_pipi_A = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.2,1.2);
	TH1D* m_Dspipi_A = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,2,5.5);
	TH1D* m_Dspi_A = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});",50,1.5,5);

        TH1D* m_Kpipi_Abar = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,1,2);
	TH1D* m_Kpi_Abar = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.6,1.2);
	TH1D* m_pipi_Abar = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,0.2,1.2);
	TH1D* m_Dspipi_Abar = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2}); ",50,2,5.5);
	TH1D* m_Dspi_Abar = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});",50,1.5,5);

	TH2D* dalitz = new TH2D("", ";m(K^{+} #pi^{-}) (GeV);m(#pi^{+} #pi^{-}) (GeV); ", 50, 0.6 ,1.2,50,0.25,1.2);
	dalitz->SetMarkerSize(0.2);


	/// Make full time-dependent PDF
   	FullAmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma
		      , xm, ym, xp, yp
		      ,Gamma, dGamma, dm
                      ,offset_sigma_dt, scale_mean_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, "Uniform" );

	pdf.beginFit();
//         ampsSig.beginFit();
//         ampsSig.redoIntegrator();
//         ampsSig_bar.beginFit();
//         ampsSig_bar.redoIntegrator();
//         vector<double> k_fit = coherenceFactor(fas,fas_bar,(double)r, (double)delta,eventListMC,eventListPhsp);
	
	cout << "Start loop " << endl;
	TCanvas* c = new TCanvas();
	TCanvas* c_1 = new TCanvas();
	c_1->Divide(3,2);
		
// 	for(int n = 1; n <= h_t->GetNbinsX(); n++){
// 		m_Kpipi->Clear();
// 		m_Kpi->Clear();
// 		m_pipi->Clear();
// 		m_Dspipi->Clear();
// 		m_Dspi->Clear();
// 		dalitz->Clear();		

		int n = step;
		double t =  h_t->GetBinLowEdge(n);
		cout << "t = " << t << endl; 

		double sumw = 0;
		double sumw_bar = 0;
		double sumw_CP = 0;

		for(int i = 0; i < eventListMC.size(); i++){
				DalitzEvent evt(eventListMC.getEvent(i));
				evt.setValueInVector(1, 0.0001);
				evt.setValueInVector(0,t);
				evt.setValueInVector(4, 0);
				evt.setValueInVector(6, 0);
				
				evt.setValueInVector(2, 1);
				evt.setValueInVector(3, 1);
				evt.setValueInVector(5, 1);

// 				r.setCurrentFitVal(0.45);
	
				double pdfVal  = pdf.getVal(evt);	
				double weight = pdfVal*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();	
	
				m_Kpipi->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight);
				m_Kpi->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight);
				m_pipi->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight);
				m_Dspipi->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight);
				m_Dspi->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight);

				dalitz->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),sqrt(evt.s(3,4)/(GeV*GeV)),weight);

				double weight_A = ampsSig.un_normalised_noPs(evt) *evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
				double weight_Abar = ampsSig_bar.un_normalised_noPs(evt) *evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();

				m_Kpipi_A->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight_A);
				m_Kpi_A->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight_A);
				m_pipi_A->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight_A);
				m_Dspipi_A->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight_A);
				m_Dspi_A->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_A);

				m_Kpipi_Abar->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight_Abar);
				m_Kpi_Abar->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight_Abar);
				m_pipi_Abar->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight_Abar);
				m_Dspipi_Abar->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight_Abar);
				m_Dspi_Abar->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_Abar);

// 				r.setCurrentFitVal(0);
				evt.setValueInVector(3, -1);
				evt.setValueInVector(5, -1);

				double pdfVal_bar  = pdf.getVal(evt);	
				double weight_bar = pdfVal_bar*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();	

				evt.setValueInVector(2, -1);
				evt.CP_conjugateYourself();

				double pdfVal_CP  = pdf.getVal(evt);	
				double weight_CP = pdfVal_CP*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();	

				m_Kpipi_CP->Fill(sqrt(evt.sij(s234)/(GeV*GeV)),weight_CP);
				m_Kpi_CP->Fill(sqrt(evt.s(2,4)/(GeV*GeV)),weight_CP);
				m_pipi_CP->Fill(sqrt(evt.s(3,4)/(GeV*GeV)),weight_CP);
				m_Dspipi_CP->Fill(sqrt(evt.sij(s134)/(GeV*GeV)),weight_CP);
				m_Dspi_CP->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_CP);

				sumw += weight;
				sumw_bar += weight_bar;
				sumw_CP += weight_CP;
		}	

		h_t->SetBinContent(n,sumw);	

		c->cd();

		double N_A = abs(sumw- sumw_bar)/(sumw + sumw_bar)/( (1-r*r)/(1+r*r) );
		double N_Abar = 1-N_A;//(sumw_bar)/(sumw + sumw_bar);

		cout << "N_A = " << N_A<< endl;
		cout << "N_Abar = " << N_Abar << endl;


		m_Kpipi->SetLineWidth(3);
		m_Kpipi->Scale(1./m_Kpipi->Integral());
		m_Kpipi->SetMaximum(0.1);
		m_Kpipi->Draw("histc");

		m_Kpipi_A->SetLineColor(kBlue+1);
		m_Kpipi_A->SetLineWidth(2);
		m_Kpipi_A->SetFillColor(kBlue+1);
		m_Kpipi_A->SetFillStyle(3353);
		m_Kpipi_A->Scale(1./m_Kpipi_A->Integral()*N_A);
		m_Kpipi_A->Draw("histcsame");

		m_Kpipi_Abar->SetLineWidth(2);
		m_Kpipi_Abar->SetLineColor(kRed+1);
		m_Kpipi_Abar->SetFillColor(kRed+1);
		m_Kpipi_Abar->SetFillStyle(3353);
		m_Kpipi_Abar->Scale(1./m_Kpipi_Abar->Integral()*N_Abar);
		m_Kpipi_Abar->Draw("histcsame");

		c->Print(((string)OutputDir+"m_Kpipi_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"m_Kpipi_"+anythingToString(n)+".png").c_str());

		m_Kpi_A->SetLineColor(kBlue+1);
		m_Kpi_A->SetLineWidth(2);
		m_Kpi_A->SetFillColor(kBlue+1);
		m_Kpi_A->SetFillStyle(3353);
		m_Kpi_A->Scale(1./m_Kpi_A->Integral()*N_A);
		m_Kpi_A->Draw("histcsame");

		m_Kpi_Abar->SetLineWidth(2);
		m_Kpi_Abar->SetLineColor(kRed+1);
		m_Kpi_Abar->SetFillColor(kRed+1);
		m_Kpi_Abar->SetFillStyle(3353);
		m_Kpi_Abar->Scale(1./m_Kpi_Abar->Integral()*N_Abar);
		m_Kpi_Abar->Draw("histcsame");

		m_Kpi->SetLineWidth(3);
		m_Kpi->Scale(1./m_Kpi->Integral());
		m_Kpi->SetMaximum(0.2);	
		m_Kpi->Draw("histc");
		c->Print(((string)OutputDir+"m_Kpi_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"m_Kpi_"+anythingToString(n)+".png").c_str());

		m_pipi_A->SetLineColor(kBlue+1);
		m_pipi_A->SetLineWidth(2);
		m_pipi_A->SetFillColor(kBlue+1);
		m_pipi_A->SetFillStyle(3353);
		m_pipi_A->Scale(1./m_pipi_A->Integral()*N_A);
		m_pipi_A->Draw("histcsame");

		m_pipi_Abar->SetLineWidth(2);
		m_pipi_Abar->SetLineColor(kRed+1);
		m_pipi_Abar->SetFillColor(kRed+1);
		m_pipi_Abar->SetFillStyle(3353);
		m_pipi_Abar->Scale(1./m_pipi_Abar->Integral()*N_Abar);
		m_pipi_Abar->Draw("histcsame");
	
		m_pipi->SetLineWidth(3);
		m_pipi->Scale(1./m_pipi->Integral());
		m_pipi->SetMaximum(0.1);
		m_pipi->Draw("histc");
		c->Print(((string)OutputDir+"m_pipi_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"m_pipi_"+anythingToString(n)+".png").c_str());
	
		m_Dspipi->SetLineWidth(3);
		m_Dspipi->Scale(1./m_Dspipi->Integral());
		m_Dspipi->SetMaximum(0.1);
		m_Dspipi->Draw("histc");
		c->Print(((string)OutputDir+"m_Dspipi_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"m_Dspipi_"+anythingToString(n)+".png").c_str());

		m_Dspipi_A->SetLineColor(kBlue+1);
		m_Dspipi_A->SetLineWidth(2);
		m_Dspipi_A->SetFillColor(kBlue+1);
		m_Dspipi_A->SetFillStyle(3353);
		m_Dspipi_A->Scale(1./m_Dspipi_A->Integral()*N_A);
		m_Dspipi_A->Draw("histcsame");

		m_Dspipi_Abar->SetLineWidth(2);
		m_Dspipi_Abar->SetLineColor(kRed+1);
		m_Dspipi_Abar->SetFillColor(kRed+1);
		m_Dspipi_Abar->SetFillStyle(3353);
		m_Dspipi_Abar->Scale(1./m_Dspipi_Abar->Integral()*N_Abar);
		m_Dspipi_Abar->Draw("histcsame");

		m_Dspi->SetLineWidth(3);
		m_Dspi->Scale(1./m_Dspi->Integral());
		m_Dspi->SetMaximum(0.1);	
		m_Dspi->Draw("histc");
		c->Print(((string)OutputDir+"m_Dspi_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"m_Dspi_"+anythingToString(n)+".png").c_str());

		m_Dspi_A->SetLineColor(kBlue+1);
		m_Dspi_A->SetLineWidth(2);
		m_Dspi_A->SetFillColor(kBlue+1);
		m_Dspi_A->SetFillStyle(3353);
		m_Dspi_A->Scale(1./m_Dspi_A->Integral()*N_A);
		m_Dspi_A->Draw("histcsame");

		m_Dspi_Abar->SetLineWidth(2);
		m_Dspi_Abar->SetLineColor(kRed+1);
		m_Dspi_Abar->SetFillColor(kRed+1);
		m_Dspi_Abar->SetFillStyle(3353);
		m_Dspi_Abar->Scale(1./m_Dspi_Abar->Integral()*N_Abar);
		m_Dspi_Abar->Draw("histcsame");


		dalitz->Draw();
		c->Print(((string)OutputDir+"dalitz_"+anythingToString(n)+".eps").c_str());
		c->Print(((string)OutputDir+"dalitz_"+anythingToString(n)+".png").c_str());
	
		TLegend leg(0.,0.,1,1,"");        
		leg.SetLineStyle(0);
		leg.SetLineColor(0);
		leg.SetFillColor(0);
		leg.SetTextFont(132);
// 		leg.SetTextColor(kRed);
		leg.SetTextSize(0.1);
		leg.SetTextAlign(12);

		stringstream ss ;
		TString label= "t = ";
		ss << std::fixed << std::setprecision(2) << t/(2*pi/dm);
		label += ss.str();
		label += "(2#pi/#Deltam_{s})";
	
		ss.str("");
		double N = (sumw - sumw_bar)/(sumw + sumw_bar);
		TString label_N= "N = ";
		ss << std::fixed << std::setprecision(2) << N;
		label_N += ss.str();
	
		leg.AddEntry((TObject*)0,label,"");
		TLegendEntry* le = leg.AddEntry(m_Kpipi_A,"B_{s}#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}","f");
		le->SetTextColor(kBlue);
		le = leg.AddEntry(m_Kpipi_Abar,"B_{s}#rightarrow#bar{B_{s}}#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}","f");
		le->SetTextColor(kRed);

		c_1->cd();
		c_1->cd(1);
		leg.Draw();
		c_1->cd(2);
		m_Dspipi->Draw("histc");
		m_Dspipi_A->Draw("histcsame");
		m_Dspipi_Abar->Draw("histcsame");
		c_1->cd(3);
		m_Dspi->Draw("histc");
		m_Dspi_A->Draw("histcsame");
		m_Dspi_Abar->Draw("histcsame");
		c_1->cd(4);
		m_Kpipi->Draw("histc");
		m_Kpipi_A->Draw("histcsame");
		m_Kpipi_Abar->Draw("histcsame");
		c_1->cd(5);
		m_Kpi->Draw("histc");
		m_Kpi_A->Draw("histcsame");
		m_Kpi_Abar->Draw("histcsame");
		c_1->cd(6);
		m_pipi->Draw("histc");
		m_pipi_A->Draw("histcsame");
		m_pipi_Abar->Draw("histcsame");
		
		c_1->Print(((string)OutputDir+"dalitz2_"+anythingToString(n)+".eps").c_str());
		c_1->Print(((string)OutputDir+"dalitz2_"+anythingToString(n)+".png").c_str());
	//}

	c->cd();
	h_t->DrawNormalized("histc",1);
	c->Print(((string)OutputDir+"h_t.eps").c_str());

		m_Kpipi_CP->SetLineWidth(3);
		m_Kpipi_CP->SetFillColor(kGreen+3);
		m_Kpipi_CP->SetLineColor(kGreen+3);
		m_Kpipi_CP->SetFillStyle(3353);
		m_Kpipi_CP->Scale(1./m_Kpipi_CP->Integral());

		m_Kpi_CP->SetLineWidth(3);
		m_Kpi_CP->SetFillColor(kGreen+3);
		m_Kpi_CP->SetLineColor(kGreen+3);
		m_Kpi_CP->SetFillStyle(3353);
		m_Kpi_CP->SetLineWidth(3);
		m_Kpi_CP->Scale(1./m_Kpi_CP->Integral());

		m_pipi_CP->SetLineWidth(3);
		m_pipi_CP->SetFillColor(kGreen+3);
		m_pipi_CP->SetLineColor(kGreen+3);
		m_pipi_CP->SetFillStyle(3353);
		m_pipi_CP->Scale(1./m_pipi_CP->Integral());

		m_Dspipi_CP->SetLineWidth(3);
		m_Dspipi_CP->SetFillColor(kGreen+3);
		m_Dspipi_CP->SetLineColor(kGreen+3);
		m_Dspipi_CP->SetFillStyle(3353);
		m_Dspipi_CP->Scale(1./m_Dspipi_CP->Integral());

		m_Dspi_CP->SetLineWidth(3);
		m_Dspi_CP->SetFillColor(kGreen+3);
		m_Dspi_CP->SetLineColor(kGreen+3);
		m_Dspi_CP->SetFillStyle(3353);
		m_Dspi_CP->Scale(1./m_Dspi_CP->Integral());

		TLegend leg_CP(0.,0.,1,1,"");        
		leg_CP.SetLineStyle(0);
		leg_CP.SetLineColor(0);
		leg_CP.SetFillColor(0);
		leg_CP.SetTextFont(132);
// 		leg_CP.SetTextColor(kRed);
		leg_CP.SetTextSize(0.1);
		leg_CP.SetTextAlign(12);
		
		leg_CP.AddEntry((TObject*)0,label,"");
		le = leg_CP.AddEntry(m_Kpipi,"B_{s}(#bar{B_{s}})#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}","l");
		le->SetTextColor(kBlack);
		le = leg_CP.AddEntry(m_Kpipi_CP,"#bar{B_{s}}(B_{s})#rightarrowD_{s}^{+}K^{-}#pi^{-}#pi^{+}","l");
		le->SetTextColor(kGreen+3);

		c_1->cd();
		c_1->cd(1);
		leg_CP.Draw();
		c_1->cd(2);
		m_Dspipi->Draw("histc");
		m_Dspipi_CP->Draw("histcsame");
		m_Dspipi->Draw("histcsame");
		c_1->cd(3);
		m_Dspi->Draw("histc");
		m_Dspi_CP->Draw("histcsame");
		m_Dspi->Draw("histcsame");
		c_1->cd(4);
		m_Kpipi->Draw("histc");
		m_Kpipi_CP->Draw("histcsame");
		m_Kpipi->Draw("histcsame");
		c_1->cd(5);
		m_Kpi->Draw("histc");
		m_Kpi_CP->Draw("histcsame");
		m_Kpi->Draw("histcsame");
		c_1->cd(6);
		m_pipi->Draw("histc");
		m_pipi_CP->Draw("histcsame");
		m_pipi->Draw("histcsame");
		
		c_1->Print(((string)OutputDir+"dalitz3_"+anythingToString(n)+".eps").c_str());
		c_1->Print(((string)OutputDir+"dalitz3_"+anythingToString(n)+".png").c_str());
	//}

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
    Int_t q_OS,f,Ds_ID,q_SS;
    Double_t w_OS,w_SS;
    double sw;
    int run,year,Ds_finalState,trigger;
    double t,dt;
    
    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm_tagged.root").c_str());
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

    tree_norm->SetBranchAddress("OS_Combination_DEC",&q_OS);
    tree_norm->SetBranchAddress("OS_Combination_PROB",&w_OS);
    tree_norm->SetBranchAddress("SS_Kaon_DEC",&q_SS);
    tree_norm->SetBranchAddress("SS_Kaon_PROB",&w_SS);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
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
	evt.P_conjugateYourself();
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

//   produceMarginalPdfs();
//   produceIntegratorFile_CP();
  //makeIntegratorFileForToys(atoi(argv[1]));
    
  //for(int i = 0; i < 200; i++)ampFit(atoi(argv[1])+i);
    ampFit(atoi(argv[1]),(string)argv[2]);
//   animate(atoi(argv[1]));

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
