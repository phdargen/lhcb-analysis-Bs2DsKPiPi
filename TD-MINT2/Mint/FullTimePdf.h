#ifndef FULLTIMEPDF_HH
#define FULLTIMEPDF_HH
// author: Philippe d'Argent
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/DalitzEventList.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/IDalitzEvent.h"
#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventPattern.h"
#include "Mint/IEventGenerator.h"
#include "Mint/TimePdfMaster.h"
#include "Mint/IReturnRealForEvent.h"
#include "Mint/IReturnComplexForEvent.h"

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
#include "RooRandom.h"
#include "RooGaussian.h"
#include "RooMultiVarGaussian.h"
#include "RooTruthModel.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include <TStyle.h>
#include <TROOT.h>
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
#include <TF2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;
using namespace RooFit ;
using namespace MINT;


// Full DsK like time PDF with additional coherence factor 
class FullTimePdf : public MINT::PdfBase<IDalitzEvent>
{
protected:
    /// fit parameters
    const FitParameter& _C;
    const FitParameter& _D;
    const FitParameter& _D_bar;
    const FitParameter& _S;
    const FitParameter& _S_bar;
    const FitParameter& _k;
   
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
    DalitzEventPattern _pat;
    DalitzEventPattern _pat_CP;

    // Limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    NamedParameter<double> _min_TAUERR;
    NamedParameter<double> _max_TAUERR;

    TF2 *_f_bkg;
    
public:
    void parametersChanged(){
    }
    void beginFit(){
        _timePdfMaster->listFitParDependencies();
    }
    void endFit(){
    }
    
    inline double un_normalised(IDalitzEvent& evt){
        
        //const double t = (double) evt.getValueFromVector(0);
        //if(t < _min_TAU || t > _max_TAU )return 0.;

        _timePdfMaster->setAllObservablesAndFitParameters(evt);
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(1.,
				   1.,
                                   _C,
                                   -_C,
                                   _k * _D,
                                   _k * _D_bar,
                                   _k * _S,
                                   _k * _S_bar
                                   );
        
        double val =
        (
         _timePdfMaster->get_cosh_term_Val(evt)
         +  _timePdfMaster->get_cos_term_Val(evt)
         +  _timePdfMaster->get_sinh_term_Val(evt)
         +  _timePdfMaster->get_sin_term_Val(evt)
         );// * _timePdfMaster->get_marginalPdfs_product(evt);
        
        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        double val = un_normalised(evt);
        
        double norm =
        _timePdfMaster->get_cosh_term_Integral(evt)
        +  _timePdfMaster->get_cos_term_Integral(evt)
        +  _timePdfMaster->get_sinh_term_Integral(evt)
        +  _timePdfMaster->get_sin_term_Integral(evt);
        
        return val/norm;
    }
    
    double getNorm(IDalitzEvent& evt){

        _timePdfMaster->setAllObservablesAndFitParameters(evt);
        
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(1.,
				   1.,
                                   _C,
                                   -_C,
                                   _k * _D,
                                   _k * _D_bar,
                                   _k * _S,
                                   _k * _S_bar
                                   );

        double norm =
        _timePdfMaster->get_cosh_term_Integral(evt)
        +  _timePdfMaster->get_cos_term_Integral(evt)
        +  _timePdfMaster->get_sinh_term_Integral(evt)
        +  _timePdfMaster->get_sin_term_Integral(evt);
        
        return norm;
    }

 inline double getSampledPdfVal(IDalitzEvent& evt){
	return _timePdfMaster->getSamplingPdfVal(evt);
    }

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_OS(evt);
    }

    std::pair<double, double> getCalibratedMistag_OS(double& eta_OS){
        return _timePdfMaster->getCalibratedMistag_OS(eta_OS);
    }

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt,double& avg_eta_os,double& p0_os,double& p1_os,double& delta_p0_os,double& delta_p1_os ){
        return _timePdfMaster->getCalibratedMistag_OS(evt, avg_eta_os, p0_os, p1_os, delta_p0_os, delta_p1_os);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_SS(evt);
    }

    std::pair<double, double> getCalibratedMistag_SS(double& eta_SS){
        return _timePdfMaster->getCalibratedMistag_SS(eta_SS);
    }

    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt,double& avg_eta_ss,double& p0_ss,double& p1_ss,double& delta_p0_ss,double& delta_p1_ss ){
        return _timePdfMaster->getCalibratedMistag_SS(evt, avg_eta_ss, p0_ss, p1_ss, delta_p0_ss, delta_p1_ss);
    }

    std::pair<double, double> getCalibratedMistag(double eta,double avg_eta,double p0,double p1,double delta_p0,double delta_p1 ){
        return _timePdfMaster->getCalibratedMistag(eta, avg_eta, p0, p1, delta_p0, delta_p1);
    }

    double getCalibratedResolution(double dt){
        return _timePdfMaster->getCalibratedResolution(dt);
    }

    double getCalibratedResolution(double& dt,double& scale_sigma_dt,double& scale_sigma_2_dt){
        return _timePdfMaster->getCalibratedResolution(dt,scale_sigma_dt,scale_sigma_2_dt);
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

    RooDataSet* sampleEvents(int N = 10000){
	return _timePdfMaster->sampleEvents(N);
    }

    DalitzEventList generateBkgToysBootstrap(int N, DalitzEventList& eventListData, int run = -1 , int trigger = -1){

	DalitzEventList eventList;
	int N_sample = eventListData.size();

	vector<int> b_indices;
	while( b_indices.size() < N )b_indices.push_back(TMath::Nint(gRandom->Uniform(0,N_sample-1)));
	sort(b_indices.begin(), b_indices.end());

	for(int i=0; i< N; i++)
	{	
		DalitzEvent evt = eventListData[b_indices[i]];

		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);
		//if(t < _min_TAU || t > _max_TAU )continue;

		eventList.Add(evt);
	}
	//saveEventListToFile(eventList);
	return eventList;
    }


    DalitzEventList generateBkgToys(int N, int run = -1, int trigger = -1,string input = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root"){

	DalitzEventList eventList;
	NamedParameter<int>  applyPhspCuts("generateBkgToys::applyPhspCuts", 0);

	if(_f_bkg==0){
		NamedParameter<int>  correlate("readBkgData::correlate", 1);
	
		double t_bkg,dt_bkg,m_bkg;
		
		TChain* tree;
		tree=new TChain("DecayTree");
		tree->Add(((string)input).c_str());
		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("*DTF*",1);
		
		tree->SetBranchAddress("Bs_DTF_MM",&m_bkg);
		tree->SetBranchAddress("Bs_DTF_TAU",&t_bkg);
		tree->SetBranchAddress("Bs_DTF_TAUERR",&dt_bkg);
		
		TH2D* h_bkg = new TH2D("h_bkg","h_bkg",100,5200,5700,100,_min_TAU,_max_TAU);

		for(int i=0; i< tree->GetEntries(); i++)
		{	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
			tree->GetEntry(i);

			h_bkg->Fill(m_bkg,t_bkg);	
		}

		TCanvas* c = new TCanvas();
		h_bkg->SetMarkerSize(0.1);
		h_bkg->Draw();
		TProfile* h_bkg_prof = h_bkg->ProfileX("p",1,-1,"");
		h_bkg_prof->SetLineColor(kRed);
		h_bkg_prof->SetMarkerColor(kRed);
		h_bkg_prof->Sumw2();
		h_bkg_prof->Rebin(4);
		cout << "bkg correlation = " << h_bkg->GetCorrelationFactor()*100 << endl;
		h_bkg_prof->Draw("histsame");
		c->Print("h_bkg.eps");

		 if(correlate)_f_bkg = new TF2("f_bkg","[0]*(1-exp(-[3]*(y-0.4)))*exp(-x*[1])*exp(-y*[2]*(1+(5200-x)*[4]))",5200,5700,_min_TAU,_max_TAU);
		 else _f_bkg = new TF2("f_bkg","[0]*(1-exp(-[3]*(y-0.4)))*exp(-x*[1])*exp(-y*[2])",5200,5700,_min_TAU,_max_TAU);
		_f_bkg->SetNpx(100);
		_f_bkg->SetNpy(100);

		_f_bkg->SetParLimits(0,1,10000000);
		_f_bkg->SetParLimits(1,0,1);
		_f_bkg->SetParLimits(2,0,10);
		_f_bkg->SetParLimits(3,-200,200);
		if(correlate)_f_bkg->SetParLimits(4,-0.1,0.1);

		_f_bkg->SetParameter(0,1000);
		_f_bkg->SetParameter(1,0.);
		_f_bkg->SetParameter(2,1.);
		_f_bkg->SetParameter(3,0.);
		if(correlate)_f_bkg->FixParameter(4,0.);

		h_bkg->Fit("f_bkg");
		if(correlate){ 
			_f_bkg->ReleaseParameter(4);
			h_bkg->Fit("f_bkg");
		}
		
		TH2D* h_fit = (TH2D*)_f_bkg->CreateHistogram();	
		h_fit->SetMarkerSize(0.1);
		h_fit->Draw();
// 		TProfile* h_fit_prof = h_fit->ProfileX("pfit",1,-1,"");
// 		h_fit_prof->SetLineColor(kRed);
// 		h_fit_prof->SetMarkerColor(kRed);
// 		h_fit_prof->Sumw2();
// 		h_fit_prof->Rebin(4);
// 		h_fit_prof->Draw("histsame");
		cout << "bkg pdf correlation = " << h_fit->GetCorrelationFactor()*100 << endl;
		c->Print("f_bkg.eps");
	}

	vector<int> s234;
   	s234.push_back(2);
        s234.push_back(3);
        s234.push_back(4);

	while(eventList.size()<N)
	{	

		double t,m;
		_f_bkg->GetRandom2(m,t);

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		double dt = marginal_vals[0] ;
		double eta_OS = marginal_vals[1] ;
		double eta_SS = marginal_vals[2] ;
		int f = (gRandom->Uniform() > 0.5) ? 1 : -1;

		DalitzEvent evt;
		if(f < 0)evt = DalitzEvent(_pat,gRandom);
		else evt = DalitzEvent(_pat_CP,gRandom);	
		
		if(!(evt.phaseSpace() > 0.))continue;
		if(applyPhspCuts)if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) continue;	
	
		evt.setWeight(1.);
		evt.setValueInVector(0, t);
		evt.setValueInVector(1, dt);   
		if(f<0)evt.setValueInVector(2, 1);
		else if(f > 0)evt.setValueInVector(2, -1);
		evt.setValueInVector(3, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(4, eta_OS);
		evt.setValueInVector(5, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(6, eta_SS);		
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);
		evt.setValueInVector(9, m);

		eventList.Add(evt);
	}
	//saveEventListToFile(eventList);
	return eventList;
    }

    DalitzEventList readBkgData(string input = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root"){

        NamedParameter<int>  correlate("readBkgData::correlate", 1);

	double t,dt,mB;
	int f;
	int q_OS;
	double Bs_ID,Ds_ID;
	Int_t q_SS;
	double eta_OS;
	Double_t eta_SS;
	double sw;
	int year,run,Ds_finalState,trigger;
	double K[4];
	double pip[4];
	double pim[4];
	double Ds_Kp[4],Ds_Km[4],Ds_pim[4],Ds[4];
	
	TChain* tree;

	tree=new TChain("DecayTree");
	tree->Add(((string)input).c_str());
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("*ID*",1);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);
	tree->SetBranchStatus("BsDTF_*P*",1);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	
	tree->SetBranchAddress("Bs_DTF_MM",&mB);
	tree->SetBranchAddress("Bs_DTF_TAU",&t);
	tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
	tree->SetBranchAddress("Ds_ID",&f);
	tree->SetBranchAddress("year",&year);
	tree->SetBranchAddress("run",&run);
	tree->SetBranchAddress("TriggerCat",&trigger);
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
	
	DalitzEventList eventList;

	for(int i=0; i< tree->GetEntries()-1; i++)
	{	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
		tree->GetEntry(i);

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		dt = marginal_vals[0] ;
		double eta_OS = marginal_vals[1] ;
		double eta_SS = marginal_vals[2] ;
		f = (gRandom->Uniform() > 0.5) ? 1 : -1;		
		
		double sign = 1.;
		TLorentzVector K_p(sign*K[0],sign*K[1],sign*K[2],K[3]);
		TLorentzVector pip_p(sign*pip[0],sign*pip[1],sign*pip[2],pip[3]);
		TLorentzVector pim_p(sign*pim[0],sign*pim[1],sign*pim[2],pim[3]);
		TLorentzVector D_p;
		TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
		TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
		TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
		D_p = D_Kp_p + D_Km_p + D_pim_p;
		
		TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
		// array of vectors
		vector<TLorentzVector> vectorOfvectors;
	
		vectorOfvectors.push_back(B_p*MeV);
		vectorOfvectors.push_back(D_p*MeV);
		vectorOfvectors.push_back(K_p*MeV);
		vectorOfvectors.push_back(pip_p*MeV);
		vectorOfvectors.push_back(pim_p*MeV);
		DalitzEvent evt;
	
		if(f < 0)evt = DalitzEvent(_pat, vectorOfvectors);
		else evt = DalitzEvent(_pat_CP, vectorOfvectors);
		
// 		t = gRandom->Exp(1./_Gamma);
// 		if(t < _min_TAU || t > _max_TAU )continue;
/*		if( dt < min_TAUERR || dt > max_TAUERR )continue;
		if(year < min_year || year > max_year) continue;
		if(eta_SS < w_min || eta_SS > w_max )continue;*/
	
		if(!(evt.phaseSpace() > 0.))continue;
			
		evt.setWeight(1.);
		evt.setValueInVector(0, t);
		evt.setValueInVector(1, dt);   
		if(f<0)evt.setValueInVector(2, 1);
		else if(f > 0)evt.setValueInVector(2, -1);
		else {
			cout << "ERROR:: Undefined final state " << f << endl;  
			throw "ERROR";
		}
		evt.setValueInVector(3, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(4, eta_OS);
		evt.setValueInVector(5, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(6, eta_SS);
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);

		if(!correlate)tree->GetEntry(i+1);
		if(mB <= 5200 || mB >= 5700 )continue;

		evt.setValueInVector(9, mB);
		eventList.Add(evt);
        }

	return eventList;
    }

    DalitzEventList generateToysRooFit(int N = 10000, int run = -1 , int trigger = -1){

	cout << "Generating " << N << " events" << endl;
	DalitzEventList eventList;

	RooDataSet* sample = sampleEvents(N);
	
	for(int i = 0; i < sample->numEntries(); i++){
	
		RooArgSet* sample_set= (RooArgSet*)sample->get(i);
	
		double t_MC = ((RooRealVar*)sample_set->find("t"))->getVal() ;
		double dt_MC = ((RooRealVar*)sample_set->find("dt"))->getVal() ;
		double eta_OS_MC = ((RooRealVar*)sample_set->find("eta_OS"))->getVal() ;
		double eta_SS_MC = ((RooRealVar*)sample_set->find("eta_SS"))->getVal() ;

		int f_MC = ((RooCategory*)sample_set->find("qf"))->getIndex() ;
		int q_OS_MC = ((RooCategory*)sample_set->find("q_OS"))->getIndex() ;
		int q_SS_MC = ((RooCategory*)sample_set->find("q_SS"))->getIndex() ; 	
	
		DalitzEvent evt(_pat,gRandom);
		evt.setValueInVector(0, t_MC);
		evt.setValueInVector(1, dt_MC);
		evt.setValueInVector(2, f_MC);
		evt.setValueInVector(3, q_OS_MC);
		evt.setValueInVector(4, eta_OS_MC);
		evt.setValueInVector(5, q_SS_MC);
		evt.setValueInVector(6, eta_SS_MC);
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);

		eventList.Add(evt);
	}
	return eventList;
   }

    DalitzEventList generateToys(int N = 10000, int run = -1 , int trigger = -1){

	time_t startTime = time(0);

	cout << "Generating " << N << " events" << endl;
	DalitzEventList eventList;

	/// Estimate max val
	vector<double> vals;	
	for(int i = 0; i < 100000; i++){
		DalitzEvent evt = generateWeightedEvent();
		double val = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); //_timePdfMaster->get_marginalPdfs_product(evt);
		vals.push_back(val);	
	}

	cout << "Now calculating maximum val " << vals.size() << endl;
	double amax,pmax;
	generalisedPareto_estimateMaximum(vals,0.999,amax,pmax);
	
	double pdf_max = 1.;
	if(!TMath::IsNaN(pmax) && pmax > 0 && pmax < 100 * amax)pdf_max = pmax;
	else if(!TMath::IsNaN(amax))pdf_max = amax;
	// for safety
 	pdf_max *= 1.5;	

	cout << "pdf_max " << pdf_max << endl;

	int N_gen = 0;
	int N_tot = 0;
	while(true){
			DalitzEvent evt = generateWeightedEvent();
			double pdfVal = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); // /_timePdfMaster->get_marginalPdfs_product(evt);
			
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
				//if (0ul == (N_gen % 500ul)) cout << "Generated event " << N_gen << "/" << N << endl;
			}		
			N_tot ++;
			if(N_gen == N)break;
	}

	cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
 	cout << "Generated " << N_gen << " events ! Efficiecy = " << (double)N_gen/(double)N_tot << endl;

// 	saveEventListToFile(eventList);
 	return eventList;
   }

    void saveEventListToFile(DalitzEventList& eventList, string name = "toys.root"){

	    NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

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
    	    TBranch* br_sw;
	    if(!addBkgToToys) br_sw = tree->Branch( "N_Bs_sw", &sw, "N_Bs_sw/D" );
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
        	
	    vector<int> s234;
   	    s234.push_back(2);
            s234.push_back(3);
            s234.push_back(4);

	    for(int i= 0; i < eventList.size(); i++){

		if(!(eventList[i].phaseSpace() > 0.))continue;

		t = eventList[i].getValueFromVector(0);
		dt = eventList[i].getValueFromVector(1);

		Ds_ID = - eventList[i].getValueFromVector(2);
		
		q_OS = eventList[i].getValueFromVector(3);
		eta_OS = eventList[i].getValueFromVector(4);

		q_SS = eventList[i].getValueFromVector(5);
		eta_SS = eventList[i].getValueFromVector(6);

		run = eventList[i].getValueFromVector(7);
		if(run == 2) year = 16;
		else if(run == 3) year = 17;
		else year = 12;
		trigger = eventList[i].getValueFromVector(8);

		//mB = eventList[i].p(0).M(); 
		mB = eventList[i].getValueFromVector(9);
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

		m_Kpipi = sqrt(eventList[i].sij(s234))/GeV;        
		m_Kpi = sqrt(eventList[i].s(2,4))/GeV;        
		m_pipi = sqrt(eventList[i].s(3,4))/GeV;        

		tree->Fill();
	     }

	     tree->Write();
	     out->Write();
	     out->Close();
    }

    DalitzEvent generateWeightedEvent(){

	while(true){
		double t_MC = gRandom->Exp(1./_Gamma);
// 		double t_MC = gRandom->Exp(1./_Gamma)+_min_TAU;
                if(t_MC > _max_TAU || t_MC < _min_TAU)continue;

		DalitzEvent evt(_pat,gRandom);

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		double dt_MC = marginal_vals[0] ;
		double eta_OS_MC = marginal_vals[1] ;
		double eta_SS_MC = marginal_vals[2] ;

		int f_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;		
	
		// true flavor
		int q_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;

 	        int q_SS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
         	int q_OS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
 
		q_OS_MC = (gRandom->Uniform() < 0.5) ? - q_OS_MC : q_OS_MC;
		q_SS_MC = (gRandom->Uniform() < 0.5) ? - q_SS_MC : q_SS_MC;
		
		eta_OS_MC = (q_OS_MC == 0) ? 0.5 : eta_OS_MC;
		eta_SS_MC = (q_SS_MC == 0) ? 0.5 : eta_SS_MC;

		if(f_MC<0)evt.CP_conjugateYourself();
		evt.setValueInVector(0, t_MC);
		evt.setValueInVector(1, dt_MC);
		evt.setValueInVector(2, f_MC);
		evt.setValueInVector(3, q_OS_MC);
		evt.setValueInVector(4, eta_OS_MC);
		evt.setValueInVector(5, q_SS_MC);
		evt.setValueInVector(6, eta_SS_MC);
		evt.setValueInVector(9, gRandom->Gaus(5379,20));

 		evt.setGeneratorPdfRelativeToPhaseSpace(_Gamma * exp(-t_MC*_Gamma) / ( ( exp(-_min_TAU*_Gamma) - exp(-_max_TAU * _Gamma) )));
// 		evt.setGeneratorPdfRelativeToPhaseSpace(_Gamma * ( exp(-t_MC*_Gamma) + _min_TAU )  / ( ( 1. - exp(-_max_TAU * _Gamma) ) + _Gamma * _min_TAU * ( _max_TAU - 0.)));

		return evt;
	}
    }

    TH1D* plotSpline() {
	 return _timePdfMaster->plotSpline();
    }
    
    FullTimePdf(const MINT::FitParameter& C, const MINT::FitParameter& D, const MINT::FitParameter& D_bar,
                const MINT::FitParameter& S, const MINT::FitParameter& S_bar, const MINT::FitParameter& k,
                const MINT::FitParameter& Gamma, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
                ,const MINT::FitParameter& offset_sigma_dt, const MINT::FitParameter& scale_mean_dt, const MINT::FitParameter& scale_sigma_dt, const MINT::FitParameter& scale_sigma_2_dt
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
    _C(C),
    _D(D),
    _D_bar(D_bar),
    _S(S),
    _S_bar(S_bar),
    _k(k),
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
    _min_TAU("min_TAU", 0.4),
    _max_TAU("max_TAU", 10.),
    _min_TAUERR("min_TAUERR", 0.),
    _max_TAUERR("max_TAUERR", 0.1)
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

        _timePdfMaster->setCP_coeff(1., 1.,
                                   _C,-_C,
                                   _k * _D, _k * _D_bar,
                                   _k * _S, _k * _S_bar  );  

        NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
        _pat = DalitzEventPattern(EventPattern.getVector());
        _pat_CP = _pat.makeCPConjugate();

	_f_bkg = 0;
    }
     
};
//

#endif
//
