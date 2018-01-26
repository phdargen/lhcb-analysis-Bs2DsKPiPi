// Time fits
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
#include "Mint/TimePdfMaster.h"
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
//using namespace RooStats;
using namespace MINT;


// Full DsK like time PDF with additional coherence factor 
class FullTimePdf : public MINT::PdfBase<IDalitzEvent>
{
protected:
    // fit parameters
    FitParameter _C;
    FitParameter _D;
    FitParameter _D_bar;
    FitParameter _S;
    FitParameter _S_bar;
    FitParameter _k;
   
    // Time pdf master
    TimePdfMaster _timePdfMaster;
    
    // limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
        _timePdfMaster.listFitParDependencies();
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        
        const double t = (double) evt.getValueFromVector(0);
        if(t < _min_TAU || t > _max_TAU )return 0.;
        _timePdfMaster.setAllObservablesAndFitParameters(evt);
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster.setCP_coeff(1.,
				   1.,
                                   _C,
                                   -_C,
                                   _k * _D,
                                   _k * _D_bar,
                                   _k * _S,
                                   -_k * _S_bar
                                   );
        
        double val =
        (
         _timePdfMaster.get_cosh_term_Val(evt)
         +  _timePdfMaster.get_cos_term_Val(evt)
         +  _timePdfMaster.get_sinh_term_Val(evt)
         +  _timePdfMaster.get_sin_term_Val(evt)
         ) * _timePdfMaster.get_marginalPdfs_Val(evt);
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        double val = un_normalised(evt);
        
        double norm =
        _timePdfMaster.get_cosh_term_Integral(evt)
        +  _timePdfMaster.get_cos_term_Integral(evt)
        +  _timePdfMaster.get_sinh_term_Integral(evt)
        +  _timePdfMaster.get_sin_term_Integral(evt);
        
        return val/norm;
    }
    
    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_SS(evt);
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
    
    FullTimePdf():
    _C("C",1,0.,0.1),
    _D("D",1,0.,0.1),
    _D_bar("D_bar",1,0.,0.1),
    _S("S",1,0.,0.1),
    _S_bar("S_bar",1,0.,0.1),
    _k("k",1,1.,0.1),
    _min_TAU("min_TAU", 0.4), _max_TAU("max_TAU", 10.)
    {
    }    
};
//

// Full DsK like time PDF in terms of r,gamma,delta with additional coherence factor 
class FullTimePdf_mod : public MINT::PdfBase<IDalitzEvent>
{
protected:
    // Fit parameters
    FitParameter& _r;
    FitParameter& _delta;
    FitParameter& _gamma;
    FitParameter& _k;
    
    // Time pdf master
    TimePdfMaster _timePdfMaster;

    // limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
        
public:
    void parametersChanged(){
    }
    void beginFit(){
        _timePdfMaster.listFitParDependencies();
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        
        const double t = (double) evt.getValueFromVector(0);
        if(t < _min_TAU || t > _max_TAU )return 0.;
        _timePdfMaster.setAllObservablesAndFitParameters(evt);
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster.setCP_coeff(
	    1.,
	    1.,
            (1.-_r*_r)/(1.+_r*_r),
           -(1.-_r*_r)/(1.+_r*_r),
            (-2.*_r * _k * cos( (_delta-_gamma)/360.*2*pi ))/(1.+_r*_r),
            (-2.*_r * _k * cos((_delta+_gamma)/360.*2*pi ))/(1.+_r*_r),
            (2.*_r * _k * sin((_delta-_gamma)/360.*2*pi))/(1.+_r*_r),
            (-2.*_r * _k * sin((_delta+_gamma)/360.*2*pi))/(1.+_r*_r)
        );
        
        double val = 
            (
                _timePdfMaster.get_cosh_term_Val(evt)
             +  _timePdfMaster.get_cos_term_Val(evt)
             +  _timePdfMaster.get_sinh_term_Val(evt)
             +  _timePdfMaster.get_sin_term_Val(evt)
             ) * _timePdfMaster.get_marginalPdfs_Val(evt);
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        double val = un_normalised(evt);
        
        double norm =
               _timePdfMaster.get_cosh_term_Integral(evt)
            +  _timePdfMaster.get_cos_term_Integral(evt)
            +  _timePdfMaster.get_sinh_term_Integral(evt)
            +  _timePdfMaster.get_sin_term_Integral(evt);
        
        return val/norm;
    }
    
    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster.getCalibratedMistag_SS(evt);
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
    
    FullTimePdf_mod( FitParameter& r, FitParameter& delta, FitParameter& gamma, FitParameter& k):
            _r(r),
            _delta(delta),
            _gamma(gamma),
            _k(k),
            _min_TAU("min_TAU", 0.4), _max_TAU("max_TAU", 10.)
    {
    }    
};
//

void fullTimeFit(){

    /// Options
    TString prefix = "";
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
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
    NamedParameter<int>  nBinst("nBinst", 40);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 10);
    
    NamedParameter<double> min_year("min_year", 11);
    NamedParameter<double> max_year("max_year", 16);

    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);
    FitParameter  k("k",1,1,1.);

    //FullTimePdf_mod t_pdf(r,delta,gamma,k);
    FullTimePdf t_pdf;

    /// Load data
    double t,dt;
    int f;
    int q_OS;
    Short_t q_SS;
    double eta_OS;
    Float_t eta_SS;
    double sw;
    int year,Ds_finalState;
    
    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add(((string)InputDir+"Data/"+(string)channel+".root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*DEC",1);
    tree_norm->SetBranchStatus("*PROB",1);
    tree_norm->SetBranchStatus("*OS",1);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("*ID*",1);
    tree_norm->SetBranchStatus("weight",1);

    tree_norm->SetBranchAddress("Bs_DTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_DTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("Ds_ID",&f);
    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGDECISION_OS",&q_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"TAGOMEGA_OS",&eta_OS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_DEC",&q_SS);
    tree_norm->SetBranchAddress("Bs_"+prefix+"SS_nnetKaon_PROB",&eta_SS);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);

    DalitzEventList eventList,eventList_f,eventList_f_bar;
    DalitzEventPattern _pat(pat);
    DalitzEvent evt_proto(_pat);
    evt_proto.generateThisToPhaseSpace();

    for(int i=0; i< tree_norm->GetEntries(); i++)
    {	
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree_norm->GetEntries() << endl;
        tree_norm->GetEntry(i);
        
        if(t < min_TAU || t > max_TAU )continue;
        if( dt < 0 || dt > 0.1 )continue;

        if(year < min_year || year > max_year) continue;
        
        DalitzEvent evt(evt_proto);
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

    /// Fit with MINT Pdf
    Neg2LL neg2LL(t_pdf, eventList);    
    
    //cout << "tau = " << endl << tau.mean() << endl <<  tau.blindedMean() << endl;
    neg2LL.getVal();    
    Minimiser mini_t(&neg2LL);
    mini_t.doFit();
    mini_t.printResultVsInput();
    //cout << "tau = " << endl << tau.mean() << endl <<  tau.blindedMean() << endl;
    
    /// Plot
    TCanvas* c = new TCanvas();
    
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

    TH1D* h_N_mixed_p_unfolded = new TH1D("h_N_mixed_p_unfolded",";t (ps);A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_m_fit");

    double N_OS = 0;
    double N_SS = 0;
    double N_OS_SS = 0;
    double N = 0;

    double w_OS = 0;
    double w_SS = 0;
    double w_OS_SS = 0;
 	
    double D_OS = 0;
    double D_SS = 0;
    double D_OS_SS = 0;
    double D_comb = 0;

    double w_OS_all = 0;
    double w_SS_all = 0;
    double D_OS_all = 0;
    double D_SS_all = 0;

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
 
        std::pair<double, double> calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eventList[i]);
        std::pair<double, double> calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eventList[i]);

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
		q_eff = q1;
	    	N_OS += eventList[i].getWeight();
	    	w_OS += w_eff * eventList[i].getWeight(); 
		D_OS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
        }
        else if( q2 != 0){
			q_eff = q2;
 			N_SS += eventList[i].getWeight();
	    		w_SS += w_eff * eventList[i].getWeight(); 
			D_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
        } 


	D_comb += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
     
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
			if(w_eff<w_max){
				h_N_mixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
				h_N_mixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			}

            }
	    else if(q_eff==0 && f_evt == 1)h_t_0p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == 1){
			h_t_pp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			if(w_eff<w_max){
				h_N_unmixed_p->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
				h_N_unmixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			}
	    }
	    else if(q_eff==-1 && f_evt == -1){
			h_t_mm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        	    	if(w_eff<w_max){
				h_N_unmixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
				h_N_unmixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			}
	    }
	    else if(q_eff==0 && f_evt == -1)h_t_0m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            else if(q_eff==1 && f_evt == -1){
			h_t_pm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			if(w_eff<w_max){
				h_N_mixed_m->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight());
				h_N_mixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
			}
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
    }     


    cout << "Tagging perfromance " << endl << endl;        
    cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;

    cout << "OS  | " << (N_OS+N_OS_SS)/N << " | " <<  (w_OS_all)/(N_OS+N_OS_SS) << " | " << D_OS_all/N << endl;
    cout << "SS  | " << (N_SS+N_OS_SS)/N << " | " <<  (w_SS_all)/(N_SS+N_OS_SS) << " | " << D_SS_all/N << endl << endl;

    cout << "OS only  | " << N_OS/N << " | " <<  w_OS/N_OS << " | " << N_OS/N * D_OS/N_OS << endl;
    cout << "SS only  | " << N_SS/N << " | " <<  w_SS/N_SS << " | " << N_SS/N * D_SS/N_SS << endl;
    cout << "OS+SS    | " << N_OS_SS/N << " | " <<  w_OS_SS/N_OS_SS << " | " << N_OS_SS/N * D_OS_SS/N_OS_SS << endl;
    cout << "Combined | " << (N_OS+N_SS+N_OS_SS)/N << " | "<<  (w_OS+w_SS+w_OS_SS)/(N_OS+N_SS+N_OS_SS) << " | " << (N_OS+N_SS+N_OS_SS)/N * D_comb/(N_OS+N_SS+N_OS_SS) << endl ;

    cout << endl << endl;        

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

    TH1D* h_N_mixed_p_fit_unfolded = new TH1D("h_N_mixed_p_fit_unfolded",";t/#tau;A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_m_fit");

    RooRealVar* r_t = new RooRealVar("t", "t",min_TAU,max_TAU);
    RooRealVar* r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.1);
    RooRealVar* r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5); 
    RooRealVar* r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5); 

    RooExponential gen_t("gen_t","gen_t", *r_t, RooRealConstant::value(1./tau));	
    RooGaussian gen_dt("gen_dt","gen_dt", *r_dt, RooRealConstant::value(0.03),RooRealConstant::value(0.015));	
    RooGaussian gen_eta_OS("gen_eta_OS","gen_eta_OS", *r_eta_OS, RooRealConstant::value(0.4),RooRealConstant::value(0.2));	
    RooGaussian gen_eta_SS("gen_eta_SS","gen_eta_SS", *r_eta_SS, RooRealConstant::value(0.5),RooRealConstant::value(0.1));

    for(int i = 0; i < 50000000; i++){
        
	/*
        const double t_MC = ranLux.Exp(tau);
        const double dt_MC = ranLux.Uniform(0., 0.1);
        
        if(t_MC < min_TAU/2. || t_MC > max_TAU*1.2 )continue;
        if( dt_MC < 0 || dt_MC > 0.1 )continue;
        
        double q_rand = ranLux.Uniform();
       
        int q_OS_MC = 0;
        if (q_rand < 1./3.) q_OS_MC = -1;
        if (q_rand > 2./3.) q_OS_MC = 1;
        
        q_rand = ranLux.Uniform();
        int q_SS_MC = 0;
        if (q_rand < 1./3.) q_SS_MC = -1;
        if (q_rand > 2./3.) q_SS_MC = 1;
        
        const double eta_OS_MC = ranLux.Uniform(0., 0.5);
        const double eta_SS_MC = ranLux.Uniform(0., 0.5);
        
        q_rand = ranLux.Uniform();
        int f_MC = 0;
        if (q_rand > .5) f_MC = -1;
        else f_MC = 1;
        */

        double t_MC = 0.;
    	while(1) {
 		double tval = ranLux.Exp(tau); 
  	    	if (tval< max_TAU && tval> min_TAU) {
         		t_MC = tval ;
			r_t->setVal(tval);
         		break ;
		}
       }
	
	gen_dt.generateEvent(1);
	double dt_MC = r_dt->getVal(); //ranLux.Uniform(0., 0.1);
        
        double q_rand = ranLux.Uniform();
       
        int q_OS_MC = 0;
        if (q_rand < 1./3.) q_OS_MC = -1;
        if (q_rand > 2./3.) q_OS_MC = 1;
        
        q_rand = ranLux.Uniform();
        int q_SS_MC = 0;
        if (q_rand < 1./3.) q_SS_MC = -1;
        if (q_rand > 2./3.) q_SS_MC = 1;
        
	gen_eta_OS.generateEvent(1);
	double eta_OS_MC = r_eta_OS->getVal();        

	gen_eta_SS.generateEvent(1);
	double eta_SS_MC = r_eta_SS->getVal(); 

        q_rand = ranLux.Uniform();
        int f_MC = 0;
        if (q_rand > .5) f_MC = -1;
        else f_MC = 1;

        DalitzEvent evt(evt_proto);

        evt.setWeight(1.);
        evt.setValueInVector(0, t_MC);
        evt.setValueInVector(1, dt_MC);
        evt.setValueInVector(2, f_MC);
        evt.setValueInVector(3, q_OS_MC);
        evt.setValueInVector(4, eta_OS_MC);
        evt.setValueInVector(5, q_SS_MC);
        evt.setValueInVector(6, eta_SS_MC);
        
        const double pdfVal = t_pdf.getVal(evt);
        double weight = pdfVal;
	weight /=  exp(-t_MC/tau) / ( tau * ( exp(min_TAU/tau) - exp(max_TAU/tau) ) );  
	weight /= gen_dt.getVal();
	weight /= gen_eta_OS.getVal();
	weight /= gen_eta_SS.getVal();
        
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
            std::pair<double, double> calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(evt);
            std::pair<double, double> calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(evt);
            
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
			if(w_eff<w_max){
				h_N_mixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
				h_N_mixed_p_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),tau),weight);
			}
		
	    }
            else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(evt.getValueFromVector(0),weight);
            else if(q_eff==1 && f_evt == 1){
			h_t_fit_pp->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
				h_N_unmixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
				h_N_unmixed_p_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),tau),weight);
			}
            }
	    else if(q_eff==-1 && f_evt == -1){
			h_t_fit_mm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
				h_N_unmixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
				h_N_unmixed_m_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),tau),weight);
			}	

            }
	    else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(evt.getValueFromVector(0),weight);
            else if(q_eff==1 && f_evt == -1){
			h_t_fit_pm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
				h_N_mixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
				h_N_mixed_m_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),tau),weight);
			}
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
        
    }

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

	h_asym_p_fit->SetLineWidth(3);
	h_asym_m_fit->SetLineWidth(3);
	h_asym_p_fit->Draw("histc");
	h_asym_m_fit->SetLineColor(kBlue);
 	h_asym_m_fit->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());	

	TH1D* h_asym_p_fit_unfolded = (TH1D*) h_N_unmixed_p_fit_unfolded->GetAsymmetry(h_N_mixed_p_fit_unfolded);	
	h_asym_p_fit_unfolded->SetLineWidth(3);
	h_asym_p_fit_unfolded->SetLineColor(kRed);
 	h_asym_p_fit_unfolded->Draw("histc");

	TH1D* h_asym_m_fit_unfolded = (TH1D*) h_N_unmixed_m_fit_unfolded->GetAsymmetry(h_N_mixed_m_fit_unfolded);
	h_asym_m_fit_unfolded->SetLineWidth(3);
	h_asym_m_fit_unfolded->SetLineColor(kBlue);
 	h_asym_m_fit_unfolded->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_unfolded.eps").c_str());	
        c->Print(((string)OutputDir+"h_asym_unfolded.pdf").c_str());	

	TH1D* h_asym_p_unfolded = (TH1D*) h_N_unmixed_p_unfolded->GetAsymmetry(h_N_mixed_p_unfolded);	
	h_asym_p_unfolded->SetLineWidth(3);
	h_asym_p_unfolded->SetLineColor(kRed);
 	h_asym_p_unfolded->Draw("e");

	TH1D* h_asym_m_unfolded = (TH1D*) h_N_unmixed_m_unfolded->GetAsymmetry(h_N_mixed_m_unfolded);
	h_asym_m_unfolded->SetLineWidth(3);
	h_asym_m_unfolded->SetLineColor(kBlue);
 	h_asym_m_unfolded->Draw("esame");

 	h_asym_p_fit_unfolded->Draw("histcsame");
 	h_asym_m_fit_unfolded->Draw("histcsame");
        c->Print(((string)OutputDir+"h_asym_unfolded_fit.eps").c_str());	

    }
    
    if(do2DScan == 1){
        cout << "Now doing 2D scan:" << endl;
        
        Neg2LL fcn(t_pdf, eventList_f);    
        Neg2LL fcn_bar(t_pdf, eventList_f_bar);    
        
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
    NamedParameter<double> min_year("min_year", 11);
    NamedParameter<double> max_year("max_year", 16);

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
        if(year < min_year || year > max_year) continue;
        
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
        if(year < min_year || year > max_year) continue;

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


int main(int argc, char** argv){

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->ProcessLine(".x ../lhcbStyle.C");

  produceMarginalPdfs();
  fullTimeFit();
  
  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
