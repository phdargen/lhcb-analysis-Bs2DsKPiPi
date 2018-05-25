// Fits Bs mass distribution and calculates sweights
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
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
#include <TFitResult.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include <RooDataSet.h>
#include <RooMCStudy.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooAddModel.h"
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooFitResult.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooMappedCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"
#include "RooBinning.h"
#include "RooBifurGauss.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooNDKeysPdf.h"
#include "RooKeysPdf.h"
#include "RooJohnsonSU.h"
#include "RooSimPdfBuilder.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
#include "Mint/RooHILLdini.h"
#include "Mint/RooHORNSdini.h"

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

vector<TString> str_year;
vector<TString> str_run;
vector<TString> str_Ds;
vector<TString> str_trigger;
static const double massKaon = 493.68;
static const double massPion = 139.57;

double prepareMisIdBkgShape(TString channel = "Dstar3pi"){

	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);

	TString fileName = "/auto/data/dargent/BsDsKpipi/Preselected/MC/";
	if(channel == "Dstar3pi")fileName += "norm_Ds2KKpi_12_Dstar_bkg.root";
	else if(channel == "Ds3pi")fileName += "norm_Ds2KKpi_bkg.root";
	else {
		cout << "ERROR::No channel specified" << endl;
		throw "ERROR";
	}
	TString calib_fileName = fileName;
	calib_fileName.ReplaceAll("bkg.root","bkg_PIDK_10.root");
	TString out_fileName = fileName;
	out_fileName.ReplaceAll("/Preselected/","/Final/");

	/// Load file
	TFile* calib_file = new TFile(calib_fileName);
	TTree* calib_tree = (TTree*) calib_file->Get("CalibTool_PIDCalibTree");
	float pi_plus1_PIDCalibEff,pi_plus2_PIDCalibEff;
	calib_tree->SetBranchAddress("pi_plus1_PIDCalibEff",&pi_plus1_PIDCalibEff);
	calib_tree->SetBranchAddress("pi_plus2_PIDCalibEff",&pi_plus2_PIDCalibEff);

	TFile* file = new TFile(fileName);
	TTree* tree = (TTree*) file->Get("DecayTree");
	
	TFile* out_file = new TFile(out_fileName,"RECREATE");
	TTree* new_tree = tree->CopyTree("3000 <= pi_plus1_P && 100000 > pi_plus1_P && 1.5 <= pi_plus1_ETA && 5 > pi_plus1_ETA && 3000 <= pi_plus2_P && 100000 > pi_plus2_P && 1.5 <= pi_plus2_ETA && 5 > pi_plus2_ETA");
	double fake_Bs_MM,EventWeight,fake_m_Kpipi,fake_m_Kpi,fake_m_pipi;
	TBranch* b_fake_Bs_MM = new_tree->Branch("fake_Bs_MM", &fake_Bs_MM, "fake_Bs_MM/D");
	TBranch* b_fake_m_Kpipi = new_tree->Branch("fake_m_Kpipi", &fake_m_Kpipi, "fake_m_Kpipi/D");
	TBranch* b_fake_m_Kpi = new_tree->Branch("fake_m_Kpi", &fake_m_Kpi, "fake_m_Kpi/D");
	TBranch* b_fake_m_pipi = new_tree->Branch("fake_m_pipi", &fake_m_pipi, "fake_m_pipi/D");
	TBranch* b_weight = new_tree->Branch("EventWeight", &EventWeight, "EventWeight/D");

	double Ds_PE,Ds_PX,Ds_PY,Ds_PZ,Bs_MM;
	double pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ;
	double pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ;
	double pi_minus_PX,pi_minus_PY,pi_minus_PZ;

        new_tree->SetBranchAddress("Bs_MM", &Bs_MM);
        new_tree->SetBranchAddress("Ds_PE", &Ds_PE);
        new_tree->SetBranchAddress("Ds_PX", &Ds_PX);
        new_tree->SetBranchAddress("Ds_PY", &Ds_PY);
        new_tree->SetBranchAddress("Ds_PZ", &Ds_PZ);
        new_tree->SetBranchAddress("pi_plus1_PX", &pi_plus1_PX);
        new_tree->SetBranchAddress("pi_plus1_PY", &pi_plus1_PY);
        new_tree->SetBranchAddress("pi_plus1_PZ", &pi_plus1_PZ);
        new_tree->SetBranchAddress("pi_plus2_PX", &pi_plus2_PX);
        new_tree->SetBranchAddress("pi_plus2_PY", &pi_plus2_PY);
        new_tree->SetBranchAddress("pi_plus2_PZ", &pi_plus2_PZ);
        new_tree->SetBranchAddress("pi_minus_PX", &pi_minus_PX);
        new_tree->SetBranchAddress("pi_minus_PY", &pi_minus_PY);
        new_tree->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ);

	int Ds_finalState;
	new_tree->SetBranchAddress("Ds_finalState", &Ds_finalState);

	if(calib_tree->GetEntries() != new_tree->GetEntries()){
		cout << "Error:: Event numbers don't match !" << endl;
		throw "ERROR";
	}

	double eff = 0;
	double eff_0 = 0;
	double eff_1 = 0;
	double eff_2 = 0;
	int n_0 = 0;
	int n_1 = 0;
	int n_2 = 0;
	TLorentzVector Ds;
	TLorentzVector pi_plus2;
	TLorentzVector pi_plus1;
	TLorentzVector pi_minus;
	for(int i= 0; i<new_tree->GetEntries();i++){
		calib_tree->GetEntry(i);
		new_tree->GetEntry(i);
		
		Ds.SetPxPyPzE(Ds_PX,Ds_PY,Ds_PZ,Ds_PE);
       		pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
		if(pi_plus1_PIDCalibEff > pi_plus2_PIDCalibEff){
			    pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massKaon);
      			    pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
			    fake_m_Kpi =  (pi_minus + pi_plus1).M() ;
			    fake_m_pipi =  (pi_minus + pi_plus2).M() ;
		}
		else {
			    pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
      			    pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon);
			    fake_m_pipi =  (pi_minus + pi_plus1).M() ;
			    fake_m_Kpi =  (pi_minus + pi_plus2).M() ;
		}

		fake_Bs_MM = (Ds + pi_minus + pi_plus1 + pi_plus2).M() ;
		fake_m_Kpipi =  (pi_minus + pi_plus1 + pi_plus2).M() ;
		
		EventWeight = max(pi_plus1_PIDCalibEff,pi_plus2_PIDCalibEff);

		if(Ds_finalState == 0) n_0 ++;
		else if(Ds_finalState == 1) n_1 ++;
		else if(Ds_finalState == 2) n_2 ++;

		if(fake_m_Kpipi < 1950. && fake_m_Kpi < 1200. && fake_m_pipi < 1200.) { 
			if(fake_Bs_MM > min_MM && fake_Bs_MM < max_MM){
				eff += EventWeight;
				if(Ds_finalState == 0) eff_0 += EventWeight;
				if(Ds_finalState == 1) eff_1 += EventWeight;
				if(Ds_finalState == 2) eff_2 += EventWeight;
			}
		}
		else fake_Bs_MM = -999;

		b_weight->Fill();
		b_fake_Bs_MM->Fill();
		b_fake_m_Kpipi->Fill();
		b_fake_m_Kpi->Fill();
		b_fake_m_pipi->Fill();
	}

	double fake_prob = eff/calib_tree->GetEntries();
	cout << "Fake prob. for " << channel << " = " << fake_prob * 100. << " % " << endl;
	cout << "Now for different Ds final states :" << endl;
	cout << eff_0/n_0 * 100. << " % " << endl;
	cout << eff_1/n_1 * 100. << " % " << endl;
	cout << eff_2/n_2 * 100. << " % " << endl;

	new_tree->Write();
	out_file->Close();

	return fake_prob;
}

vector<double> fitPartRecoBkgShape(){

        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);
	///define shape of Bs->Ds*pipipi BG as 3 bifurcated gaussians
	RooRealVar Bs_MM("Bs_DTF_MM", "m(D_{s}*K#pi#pi)", 5000., 5350.,"MeV/c^{2}");
	//mean of gaussians
	RooRealVar mean1("mean1","mu", 5059.,5040.,5070.);
	RooRealVar mean2("mean2","mu", 5182.,5140.,5205.);
	RooRealVar mean3("mean3","mu", 5285.,5270.,5300.);
	//width of gaussians
	RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", 25.9,15.,40.);
	RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", 99.4,50.,115.);
	RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", 13.1,5.,100.);
	RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", 49.5,25.,70.);
	RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", 107.,10.,125.);
	RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", 21.1,5.,33.);
	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);
	//fractions of gauss functions
	RooRealVar f_1("f_{1}", "fraction1", 0.405, 0., 1.);
	RooRealVar f_2("f_{2}", "fraction2", 0.1, 0., 1.);
	//add all gaussians
	RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2),kTRUE);
	
	/// Load file
	TFile* file = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12_Dstar_bkg.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs*MM",1);
	tree->SetBranchStatus("m_*",1);

	//TTree *new_tree = tree->CopyTree("m_Kpipi < 1950 && m_Kpi < 1200 && m_pipi < 1200");
	//Fill needed variable in RooDataSet
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));
	
	/// Fit
	RooFitResult *result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 
	//plot mass distribution and fit results
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");
        frame_m->GetXaxis()->SetTitle("m(D_{s}^{-}#pi^{+}#pi^{+}#pi^{-}) [MeV/c^{2}]");
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(50));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(BifGauss1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(BifGauss2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(BifGauss3),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstartpipipi.eps");
	if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/BkgShape/Bs2Dsstartpipipi.pdf");
	Bs_MM.setRange("signal_range",min_MM,max_MM);
	double eff_signal_range = pdf->createIntegral(Bs_MM,NormSet(Bs_MM),Range("signal_range"))->getVal();
	cout << "eff_signal_range" << eff_signal_range << endl;

	/// Return fit parameters
	vector<double> params;
	params.push_back(mean1.getVal());
	params.push_back(mean2.getVal());
	params.push_back(mean3.getVal());
	params.push_back(sigmaL1.getVal());
	params.push_back(sigmaR1.getVal());
	params.push_back(sigmaL2.getVal());
	params.push_back(sigmaR2.getVal());
	params.push_back(sigmaL3.getVal());
	params.push_back(sigmaR3.getVal());
	params.push_back(f_1.getVal());
	params.push_back(f_2.getVal());
	params.push_back(eff_signal_range);

	file->Close();
	
	return params;
}

vector<double> fitMisIdBkgShape_Ds3pi(){

        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
	/// Define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
	RooRealVar Bs_Mass("fake_Bs_MM", "m(D_{s}^{-} #pi_{K}^{+}#pi^{+}#pi^{-})", 5300., 6000.,"MeV/c^{2}");
	RooRealVar EventWeight("EventWeight", "EventWeight", 0.);
	RooRealVar Ds_finalState("Ds_finalState","Ds_finalState", 0.);
	RooRealVar Bs_MM("Bs_MM", "m(D_{s}^{-} K^{+}#pi^{-}#pi^{-})", 5320., 5420.,"MeV/c^{2}");

	//mean of crrystal balls
	RooRealVar mean1("mean1","mu", 5444.,5400.,5490.);
	RooRealVar mean2("mean2","mu", 5517.,5450.,5650.);
	// asymmetry parameter of crystsal balls
	RooRealVar a1("a1","a1",-1.5, -4.,3.);
	RooRealVar a2("a2","a2",-0.6, -1.5,1.5);
	RooRealVar n1("n1","n1",0.3, 0.,2.);
	RooRealVar n2("n2","n2",10., 0.,200.);
	//sigma of crystal balls
	RooRealVar sigma1("sigma_{1}", "sigma1", 24.,5.,300.);
	RooRealVar sigma2("sigma_{2}", "sigma2", 89.,5.,300.);
	//crystal Balls
	RooCBShape CB1("CB1", "CB1", Bs_Mass, mean1, sigma1, a1, n1);
	RooCBShape CB2("CB2", "CB2", Bs_Mass, mean2, sigma2, a2, n2);
	//fraction of crystal balls
	RooRealVar f_1("f_{1}", "fraction1", 0.5, 0., 1.);
	//add all gaussians
	RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1, CB2), RooArgList(f_1));
	
	///Load file
	TFile* file = new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm_Ds2KKpi_bkg.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("EventWeight",1);
	tree->SetBranchStatus("*Bs_MM",1);
	tree->SetBranchStatus("Ds_finalState",1);
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_Mass,Ds_finalState,Bs_MM,EventWeight),Import(*tree),WeightVar(EventWeight));
	
	/// Fit
	RooFitResult *result;
	result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3),SumW2Error(kTRUE));
	cout << "result is --------------- "<<endl;
	result->Print(); 
	
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
        frame_m->GetXaxis()->SetTitle("m(D_{s}^{-}#pi^{+}_{K}#pi^{+}#pi^{-}) [MeV/c^{2}]");

	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40),DataError(RooAbsData::SumW2));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi.eps");
	if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/BkgShape/Bs2Dspipipi_as_DsKpipi.pdf");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	RooDataSet* data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 0"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi_0.eps");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 1"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi_1.eps");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 2"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi_2.eps");

	/// Return fit params
	vector<double> params;
	params.push_back(mean1.getVal());
	params.push_back(mean2.getVal());
	params.push_back(a1.getVal());
	params.push_back(a2.getVal());
	params.push_back(n1.getVal());
	params.push_back(n2.getVal());
	params.push_back(sigma1.getVal());
	params.push_back(sigma2.getVal());
	params.push_back(f_1.getVal());

	return params;
}

vector<double> fitMisIdBkgShape_Dsstar3pi(){

        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
	/// Define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
	RooRealVar Bs_Mass("fake_Bs_MM", "m(D_{s}^{-} #pi_{K}^{+}#pi^{+}#pi^{-})", 4900., 6200.,"MeV/c^{2}");
	RooRealVar EventWeight("EventWeight","EventWeight", 0.);
	RooRealVar Ds_finalState("Ds_finalState","Ds_finalState", 0.);

	//mean of crrystal balls
	RooRealVar mean1("mean1","mu", 5350.,5200.,5490.);
	RooRealVar mean2("mean2","mu", 5517.,5400.,5650.);
	// asymmetry parameter of crystsal balls
	RooRealVar a1("a1","a1",-1.5, -3.,2.5);
	RooRealVar a2("a2","a2",-0.5, -2.5,2.5);
	RooRealVar n1("n1","n1",5.0, 0.,10.);
	RooRealVar n2("n2","n2",10., 0.,50.);
	//sigma of crystal balls
	RooRealVar sigma1("sigma_{1}", "sigma1", 70.,15.,180.);
	RooRealVar sigma2("sigma_{2}", "sigma2", 89.,15.,200.);
	//crystal Balls
	RooCBShape CB1("CB1", "CB1", Bs_Mass, mean1, sigma1, a1, n1);
	RooCBShape CB2("CB2", "CB2", Bs_Mass, mean2, sigma2, a2, n2);
	//fraction of crystal balls
	RooRealVar f_1("f_{1}", "fraction1", 0.5, 0.2, 0.8);
	//add all gaussians
	RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1, CB2), RooArgList(f_1));
	
	/// Load file
	TFile* file = new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm_Ds2KKpi_12_Dstar_bkg.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("EventWeight",1);
	tree->SetBranchStatus("fake_Bs_MM",1);
	tree->SetBranchStatus("Ds_finalState",1);
	//Fill needed variable in RooDataSet
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_Mass,Ds_finalState,EventWeight),Import(*tree),WeightVar(EventWeight));
	
	/// Fit
	RooFitResult *result;
	result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3),SumW2Error(kTRUE));
	cout << "result is --------------- "<<endl;
	result->Print(); 
	
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
        frame_m->GetXaxis()->SetTitle("m(D_{s}^{-}#pi^{+}_{K}#pi^{+}#pi^{-}) [MeV/c^{2}]");
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.eps");
	if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.pdf");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	RooDataSet* data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 0"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi_0.eps");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 1"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi_1.eps");

	frame_m= Bs_Mass.frame();
	frame_m->SetTitle("");
	data_slice = (RooDataSet*)data->reduce(Cut("Ds_finalState  == 2"));
	data_slice->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi_2.eps");

	/// Return fit params
	vector<double> params; 
	params.push_back(mean1.getVal());
	params.push_back(mean2.getVal());
	params.push_back(a1.getVal());
	params.push_back(a2.getVal());
	params.push_back(n1.getVal());
	params.push_back(n2.getVal());
	params.push_back(sigma1.getVal());
	params.push_back(sigma2.getVal());
	params.push_back(f_1.getVal());
	
return params;
}

vector<double> fitSignalShape(TString channel = "signal"){

        /// Options
	NamedParameter<double> cut_BDT("cut_BDT_MC",0.);
        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
	NamedParameter<int> fitPreselected("fitPreselected", 0);
	NamedParameter<int> sWeight("sWeightMC", 0);
	double min_MM = 5320. ;
	double max_MM = 5420. ;

	/// Load file
	TString inFileName = "/auto/data/dargent/BsDsKpipi/BDT/MC/"+channel+".root";
	TChain* tree = new TChain("DecayTree");	
	if(fitPreselected) {
		tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/"+channel+"_Ds2*_11.root");
		tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/"+channel+"_Ds2*_12.root");
	}
	else tree->Add(inFileName);

	if(!sWeight){
		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("Bs_BKGCAT",1);
		tree->SetBranchStatus("Bs_TRUEID",1);
		tree->SetBranchStatus("Bs_DTF_MM",1);
		if(!fitPreselected)tree->SetBranchStatus("BDTG_response",1);
		tree->SetBranchStatus("Ds_finalState",1);
	}
	tree->SetBranchStatus("weight",0);

	TFile* output;
	if(!sWeight || fitPreselected) output = new TFile("dummy.root","RECREATE");
	else output = new TFile(inFileName.ReplaceAll("/BDT/","/Final/"),"RECREATE");

	TTree* out_tree;
	if(fitPreselected)out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) ).c_str() );
	else out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) + " && BDTG_response > " + anythingToString((double)cut_BDT) ).c_str() );

	int bkgCAT,BKGCAT,Bs_TRUEID;
	double sw;
        TBranch* b_bkgCAT = out_tree->Branch("bkgCAT",&bkgCAT,"bkgCAT/I");
	out_tree->SetBranchAddress("Bs_BKGCAT",&BKGCAT);
	out_tree->SetBranchAddress("Bs_TRUEID",&Bs_TRUEID);
    	TBranch* b_w = out_tree->Branch("weight", &sw, "weight/D");

	for(int i= 0; i< out_tree->GetEntries();i++){
		out_tree->GetEntry(i);
		if(BKGCAT == 20 )bkgCAT= 0;
		//else if(BKGCAT == 60 || ( BKGCAT == 50 && abs(Bs_TRUEID) == 531) )bkgCAT= 1;
		else bkgCAT= 1;
		b_bkgCAT->Fill();
	}

	TString channelString;
        if(channel == "norm") channelString = "m(D_{s}#pi#pi#pi)" ;
        if(channel == "signal") channelString = "m(D_{s}K#pi#pi)" ;

	RooRealVar DTF_Bs_M("Bs_DTF_MM", channelString, min_MM, max_MM,"MeV/c^{2}");
	RooCategory Bs_BKGCAT("bkgCAT","bkgCAT");
	Bs_BKGCAT.defineType("signal",0);
	Bs_BKGCAT.defineType("ghost",1);

	RooArgList list =  RooArgList(DTF_Bs_M,Bs_BKGCAT);
        RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree));
	
	/// Signal pdf
	RooRealVar mean("mean", "mean", 5366.89,5350.,5390.); 
	RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
	RooRealVar gamma("gamma", "gamma", -0.5,-5,5.); 
	RooRealVar delta("delta", "delta", 0.5,-5,5.); 
	RooJohnsonSU* signal= new RooJohnsonSU("signal","signal",DTF_Bs_M, mean,sigma,gamma,delta);

	RooRealVar mean_ghost("mean_ghost", "mean_ghost", 5366.89,5350.,5390.); 
	RooRealVar sigma_ghost("sigma_ghost", "sigma_ghost", 20.,0.,80.); 
	RooRealVar gamma_ghost("gamma_ghost", "gamma_ghost", -0.5,-5,5.); 
	RooRealVar delta_ghost("delta_ghost", "delta_ghost", 0.5,-5,5.); 
	RooJohnsonSU* signal_ghost= new RooJohnsonSU("signal","signal_ghost",DTF_Bs_M, mean,sigma_ghost,gamma,delta);

	/// Bkg pdf
	RooRealVar c0_ghost("c0_ghost", "c0_ghost", .0,-10,10); 
	RooRealVar c1_ghost("c1_ghost", "c1_ghost", .0,-10,10); 
	RooRealVar c2_ghost("c2_ghost", "c2_ghost", .0,-10,10); 
	RooChebychev* bkg_ghost= new RooChebychev("bkg_ghost","bkg_ghost",DTF_Bs_M, RooArgList(c0_ghost));

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.9, 0., data->numEntries());
	RooRealVar n_sig_ghost("n_sig_ghost", "n_sig_ghost", data->numEntries()*0.5, 0., data->numEntries());
	RooRealVar n_bkg_ghost("n_bkg_ghost", "n_bkg_ghost", data->numEntries()/10., 0., data->numEntries());

	RooAddPdf* pdf_signal = new RooAddPdf("pdf_signal", "pdf_signal", RooArgList(*signal), RooArgList(n_sig));
	RooAddPdf* pdf_ghost = new RooAddPdf("pdf_ghost", "pdf_ghost", RooArgList(*signal_ghost, *bkg_ghost), RooArgList(n_sig_ghost, n_bkg_ghost));

	RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",Bs_BKGCAT);
	simPdf->addPdf(*pdf_signal,"signal");
	simPdf->addPdf(*pdf_ghost,"ghost");

	/// Fit
	RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),NumCPU(3),Extended(kTRUE));
	result->Print();

	if(sWeight){
		/// Calculate ghost weights
		RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"bkgCAT==bkgCAT::ghost");
		SPlot sPlot("sPlot","sPlot",*data_slice,pdf_ghost,RooArgList(n_sig_ghost,n_bkg_ghost)); 
	
		int n_ij = 0;  /// labels entry number of data slice
		for(int n = 0; n < out_tree->GetEntries(); n++){
			b_bkgCAT->GetEntry(n);
			if(bkgCAT == 1){
				sw = sPlot.GetSWeight(n_ij,"n_sig_ghost_sw");
				n_ij++;
			}
			else sw = 1.;
			b_w->Fill();
		}
	}
	/// Plotting
	TCanvas* c = new TCanvas();

	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
 	data->plotOn(frame,Name("data"),Binning(50),Cut("bkgCAT==bkgCAT::signal"));
	simPdf->plotOn(frame,Name("signal"),ProjWData(Bs_BKGCAT,*data),Slice(Bs_BKGCAT,"signal"));
	frame->Draw();
	c->Print("eps/SignalShape/"+channel+"MC.eps");

        RooPlot* frame2= DTF_Bs_M.frame();
 	data->plotOn(frame2,Name("data"),Binning(50),Cut("bkgCAT==bkgCAT::ghost"));
 	simPdf->plotOn(frame2,Name("signal"),Slice(Bs_BKGCAT,"ghost"),ProjWData(Bs_BKGCAT,*data));
	simPdf->plotOn(frame2,Name("bkg"),ProjWData(Bs_BKGCAT,*data),Slice(Bs_BKGCAT,"ghost"),LineColor(kRed),LineStyle(kDashed),Components("bkg_ghost"));
	frame2->Draw();
	c->Print("eps/SignalShape/"+channel+"MC_ghost.eps");

	TCanvas* canvas = new TCanvas();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);
        
        double max = 5.0 ;
        double min = -5.0 ;
        double rangeX = max-min;
        double zero = max/rangeX;
        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(max);
        graph->SetMinimum(min);
        graph->SetPoint(1,min_MM,0);
        graph->SetPoint(2,max_MM,0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(max);
        graph2->SetMinimum(min);
        graph2->SetPoint(1,min_MM,-3);
        graph2->SetPoint(2,max_MM,-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(max);
        graph3->SetMinimum(min);
        graph3->SetPoint(1,min_MM,3);
        graph3->SetPoint(2,max_MM,3);
        graph3->SetLineColor(kRed);
       
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame->GetYaxis()->SetRangeUser(0.01,frame->GetMaximum()*1.);
        frame->Draw();
        
        canvas->cd();
        TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2->SetBorderMode(0);
        pad2->SetBorderSize(-1);
        pad2->SetFillStyle(0);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        RooPlot* frame_p = DTF_Bs_M.frame();
	frame_p->SetTitle("");
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
        frame_p->GetXaxis()->SetTitle( channelString + "[MeV/c^{2}]");
        
        RooHist* hpull  = frame->pullHist("data","signal");
	hpull->SetTitle("");
        frame_p->addPlotable(hpull,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->SaveAs("eps/SignalShape/"+channel+"MC_pull.eps");
	if(updateAnaNotePlots && !fitPreselected) canvas->Print("../../../../../TD-AnaNote/latex/figs/MassFit/"+channel+"MC_pull.pdf");

	TCanvas* canvas_ghost = new TCanvas();
        canvas_ghost->SetTopMargin(0.05);
        canvas_ghost->SetBottomMargin(0.05);

        TPad* pad1_ghost = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1_ghost->SetBorderMode(0);
        pad1_ghost->SetBorderSize(-1);
        pad1_ghost->SetBottomMargin(0.);
        pad1_ghost->Draw();
        pad1_ghost->cd();
        frame2->GetYaxis()->SetRangeUser(0.01,frame2->GetMaximum()*1.);
        frame2->Draw();
        
        canvas_ghost->cd();
        TPad* pad2_ghost = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2_ghost->SetBorderMode(0);
        pad2_ghost->SetBorderSize(-1);
        pad2_ghost->SetFillStyle(0);
        pad2_ghost->SetTopMargin(0.);
        pad2_ghost->SetBottomMargin(0.35);
        pad2_ghost->Draw();
        pad2_ghost->cd();
        
        RooPlot* frame_p_ghost = DTF_Bs_M.frame();
	frame_p_ghost->SetTitle("");
        frame_p_ghost->GetYaxis()->SetNdivisions(5);
        frame_p_ghost->GetYaxis()->SetLabelSize(0.12);
        frame_p_ghost->GetXaxis()->SetLabelSize(0.12);
        frame_p_ghost->GetXaxis()->SetTitleOffset(0.75);
        frame_p_ghost->GetXaxis()->SetTitleSize(0.2);
        frame_p_ghost->GetXaxis()->SetTitle( channelString + "[MeV/c^{2}]");
        
        RooHist* hpull_ghost  = frame2->pullHist("data","signal");
	hpull_ghost->SetTitle("");
        frame_p_ghost->addPlotable(hpull_ghost,"BX");
        frame_p_ghost->GetYaxis()->SetRangeUser(min,max);
        
        frame_p_ghost->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2_ghost->Update();
        canvas_ghost->Update();
        canvas_ghost->SaveAs("eps/SignalShape/"+channel+"MC_ghost_pull.eps");
	if(updateAnaNotePlots && !fitPreselected) canvas_ghost->Print("../../../../../TD-AnaNote/latex/figs/MassFit/"+channel+"MC_ghost_pull.pdf");

	/// Return fit params
	vector<double> params;
	params.push_back(mean.getVal());
	params.push_back(sigma.getVal());
	params.push_back(gamma.getVal());
	params.push_back(delta.getVal());

	cout << endl << "Fraction of signal classified as ghosts = " << n_sig_ghost.getVal()/(n_sig.getVal()+n_sig_ghost.getVal()) << endl;

	out_tree->Write();
	output->Close();
	return params;
}

vector<double> fitPartRecoBkgShapeHILLHORN(){

        RooRealVar Bs_MM("Bs_DTF_MM", "m(D_{s}*K#pi#pi)", 5000., 5350.,"MeV/c^{2}");

	///define shape of Bs->Ds*pipipi BG as RooHILL and RooHORN
	RooRealVar a("a","a", 5059.,3040.,5100.);
	RooRealVar b("b","b", 5182.,4140.,5300.);
	RooRealVar a_HORNS("a_HORNS","a_HORNS", 5059.,3040.,5100.);
	RooRealVar b_HORNS("b_HORNS","b_HORNS", 5182.,4140.,5300.);
	RooRealVar csi("csi","csi", 1.,0.,5.);
	RooRealVar csi_HORNS("csi_HORNS","csi_HORNS", 1.,0.,5.);
	RooRealVar shift("shift","shift", 1.,0.,500.);
	RooRealVar sigma_HILL("sigma_HILL", "sigma_HILL", 25.9,0.,100.);
	RooRealVar sigma_HORNS("sigma_HORNS", "sigma_HORNS", 25.9,0.,100.);
	//RooRealVar ratio_sigma("ratio_sigma", "ratio_sigma", 1.,0.,50.);
	RooRealVar ratio_sigma("ratio_sigma", "ratio_sigma", 4.81);
	//RooRealVar fraction_sigma("fraction_sigma", "fraction_sigma", 0.5,0.,1.);
	RooRealVar fraction_sigma("fraction_sigma", "fraction_sigma", 1);

	RooRealVar f_1("f_{1}", "fraction1", 0.405, 0., 1.);

	RooHILLdini RooHILLBkgShape("RooHILLBkgShape", "RooHILLBkgShape", Bs_MM, a, b, csi, shift, sigma_HILL, ratio_sigma, fraction_sigma );
	RooHORNSdini RooHORNSBkgShape("RooHORNSBkgShape", "RooHORNSBkgShape", Bs_MM, a_HORNS, b_HORNS, csi_HORNS, shift, sigma_HORNS, ratio_sigma, fraction_sigma );
	RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(RooHILLBkgShape, RooHORNSBkgShape), RooArgList(f_1),kTRUE);

	
	/// Load file
	TFile* file = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12_Dstar_bkg.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs*MM",1);
	tree->SetBranchStatus("m_*",1);

	//TTree *new_tree = tree->CopyTree("m_Kpipi < 1950 && m_Kpi < 1200 && m_pipi < 1200");
	//Fill needed variable in RooDataSet
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));
	
	/// Fit
	RooFitResult *result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 
	//plot mass distribution and fit results
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(50));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(RooHILLBkgShape),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(RooHORNSBkgShape),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstartpipipi_HORNHILL.eps");

	/// Return fit parameters
	vector<double> params;
	params.push_back(a.getVal());
	params.push_back(b.getVal());
	params.push_back(csi.getVal());
	params.push_back(shift.getVal());
	params.push_back(sigma_HILL.getVal());
	params.push_back(ratio_sigma.getVal());
	params.push_back(fraction_sigma.getVal());
	params.push_back(a_HORNS.getVal());
	params.push_back(b_HORNS.getVal());
	params.push_back(csi_HORNS.getVal());
	params.push_back(sigma_HORNS.getVal());

	file->Close();
	
	return params;
}

vector< vector<double> > fitNorm(){

	/// Options
	NamedParameter<int> fixExpBkgFromSidebands("fixExpBkgFromSidebands", 0);
        NamedParameter<int> ignorePartRecoBkg("ignorePartRecoBkg", 0);
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> sWeight("sWeightNorm", 0);
	NamedParameter<int> nBins("nBins", 80);
	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);
	NamedParameter<string> cut_BDT("cut_BDT",(string)"");
	NamedParameter<string> inFileName("inFileNameNorm",(string)"/auto/data/dargent/BsDsKpipi/BDT/Data/norm.root");
	NamedParameter<string> outFileName("outFileNameNorm",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
        NamedParameter<int> altPartBkg("altPartBkg", 0);
        NamedParameter<int> fixSignalShapeFromMC("fixSignalShapeFromMC", 0);
	NamedParameter<int> useB0("useB0", 0);
	NamedParameter<int> useTriggerCat("useTriggerCat", 1);
	NamedParameter<int> fitPreselected("fitPreselected", 0);

	/// Define categories
	RooCategory year("year","year") ;
	year.defineType("y11",11);
	year.defineType("y12",12);
	year.defineType("y15",15);
	year.defineType("y16",16);
	year.defineType("y17",17);

  	RooCategory run("run","run") ;
  	run.defineType("Run1",1) ;
  	run.defineType("Run2",2) ;

	RooCategory Ds_finalState_mod("Ds_finalState_mod","Ds_finalState_mod") ;
	Ds_finalState_mod.defineType("phipi",0);
	Ds_finalState_mod.defineType("KsK",1);
	Ds_finalState_mod.defineType("KKpi_NR",2);
	Ds_finalState_mod.defineType("pipipi",3);
	//Ds_finalState.defineType("Kpipi",4);

	RooCategory TriggerCat("TriggerCat","TriggerCat") ;
	TriggerCat.defineType("t0",0);
	TriggerCat.defineType("t1",1);

	/// Load file
	TFile *file = 0;
	TTree* tree;	
   	if(fitPreselected){
		TChain* chain = new TChain("DecayTree");
		chain->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_11.root");
		chain->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_12.root");
		chain->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_15.root");
		chain->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_16.root");
		tree = (TTree*) chain;
	}
	else{
		 file = new TFile(((string)inFileName).c_str());
		 tree = (TTree*) file->Get("DecayTree");
	}
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	if(!fitPreselected)tree->SetBranchStatus("BDTG",1);
	tree->SetBranchStatus("Ds_finalState*",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);
	tree->SetBranchStatus("*PIDK",1);

        RooRealVar DTF_Bs_M("Bs_DTF_MM", "m(D_{s}^{-}#pi^{+}#pi^{+}#pi^{-})", min_MM, max_MM,"MeV/c^{2}");
        RooRealVar BDTG("BDTG", "BDTG", 0.);              
        RooRealVar pi_plus1_PIDK("pi_plus1_PIDK", "pi_plus1_PIDK", 0.);              
        RooRealVar pi_plus2_PIDK("pi_plus2_PIDK", "pi_plus2_PIDK", 0.);              

	RooArgList list =  RooArgList(DTF_Bs_M,Ds_finalState_mod,year,run,TriggerCat,pi_plus1_PIDK,pi_plus2_PIDK);
	if(!fitPreselected)list.add(BDTG);	
	RooDataSet*  data;
	if(!fitPreselected)data = new RooDataSet("data","data",tree,list, ((string)cut_BDT).c_str() );
	else data = new RooDataSet("data","data",tree,list);

	/// Signal Pdf
	vector<double> sig_params = fitSignalShape("norm");
	RooRealVar mean_MC("mean_MC", "#mu MC", sig_params[0]); 
	RooRealVar sigma_MC("sigma_MC", "#sigma MC", sig_params[1]);
	RooRealVar scale_mean("scale_mean", "scale #mu",1.,0.5,2.); 
	RooRealVar scale_sigma("scale_sigma", "scale #sigma", 1.2, 0.5,2.);
	RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
	RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
	RooRealVar alpha("alpha", "#alpha", sig_params[2],-5.,5.); 
	RooRealVar beta("beta", "#beta", sig_params[3]);
	if(fixSignalShapeFromMC){
		alpha.setConstant();
		beta.setConstant();
	} 
	RooJohnsonSU signal("signal","signal",DTF_Bs_M, mean,sigma,alpha,beta);

	/// add B⁰ -> D_s 3pi
	RooRealVar mean_MC_B0("mean_MC_B0", " #mu MC_B0", (sig_params[0] - 87.42)); 
	RooFormulaVar mean_B0("mean_B0","@0 * @1", RooArgSet(scale_mean,mean_MC_B0)); 
	RooJohnsonSU signal_B0("signal_B0","signal_B0",DTF_Bs_M, mean_B0,sigma,alpha,beta);

	/// Combinatorial bkg pdf
// 	RooRealVar exp_par("exp_par","#lambda",-1.6508e-03,-10.,0.);	
// 	RooExponential bkg_exp("bkg_exp","exponential bkg",DTF_Bs_M,exp_par);
// 	bkg_exp.fitTo(*data,Save(kTRUE),Range(5600.,5800.));

	RooRealVar exp_par1("exp_par1","exp_par1",-1.6508e-03,-10.,10.);	
	RooRealVar exp_par("exp_par","exp_par",0.);	
	RooExponential bkg_exp1("bkg_exp1","bkg_exp1",DTF_Bs_M,exp_par);
	RooRealVar exp_par2("exp_par2","exp_par2",0,-10.,10.);	
	RooRealVar exp_par3("exp_par3","exp_par3",0,-10.,10.);	

	RooChebychev bkg_exp2("bkg_exp2","bkg_exp2",DTF_Bs_M,RooArgList(exp_par2,exp_par3));
	RooRealVar f_bkg_exp("f_bkg_exp", "f_bkg_exp", 0.);
	RooAddPdf bkg_exp("bkg_exp", "bkg_exp", RooArgList(bkg_exp1, bkg_exp2), RooArgList(f_bkg_exp));

	/// Part. reco bkg
	vector<double> bkg_partReco_params(20,0);
	if(!altPartBkg) bkg_partReco_params = fitPartRecoBkgShape();
	RooRealVar mean1("mean1","mu", bkg_partReco_params[0]);
	RooRealVar mean2("mean2","mu", bkg_partReco_params[1]);
	RooRealVar mean3("mean3","mu", bkg_partReco_params[2]);
	RooRealVar sigmaL1("sigmaL1", "sigmaL1",  bkg_partReco_params[3]);
	RooRealVar sigmaR1("sigmaR1", "sigmaR1",  bkg_partReco_params[4]);
	RooRealVar sigmaL2("sigmaL2", "sigmaL2",  bkg_partReco_params[5]);//, bkg_partReco_params[5]*0.,bkg_partReco_params[5]*5);
	RooRealVar sigmaR2("sigmaR2", "sigmaR2",  bkg_partReco_params[6]);//, bkg_partReco_params[6]*0.,bkg_partReco_params[6]*2.);
	RooRealVar sigmaL3("sigmaL3", "sigmaL3",  bkg_partReco_params[7]);//,bkg_partReco_params[7]*0.5,bkg_partReco_params[7]*2.);
	RooRealVar sigmaR3("sigmaR3", "sigmaR3",  bkg_partReco_params[8]);// ,bkg_partReco_params[8]*0.5,bkg_partReco_params[8]*2.);
	RooRealVar f_1("f_1", "f_1", bkg_partReco_params[9]);//,0,1);
	RooRealVar f_2("f_2", "f_2", bkg_partReco_params[10]);//,0,1);

	RooRealVar scale_mean_partReco("scale_mean_partReco", "scale_mean_partReco", 1.);
	RooRealVar scale_sigma_partReco("scale_sigma_partReco", "scale_sigma_partReco", 1.,0.5,2.);
	RooFormulaVar mean1_scaled("mean1_scaled","@0 * @1", RooArgSet(scale_mean_partReco,mean1)); 
	RooFormulaVar mean2_scaled("mean2_scaled","@0 * @1", RooArgSet(scale_mean_partReco,mean2)); 
	RooFormulaVar mean3_scaled("mean3_scaled","@0 * @1", RooArgSet(scale_mean_partReco,mean3)); 
	RooFormulaVar sigmaL1_scaled("sigmaL1_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaL1)); 
	RooFormulaVar sigmaR1_scaled("sigmaR1_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaR1)); 
	RooFormulaVar sigmaL2_scaled("sigmaL2_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaL2)); 
	RooFormulaVar sigmaR2_scaled("sigmaR2_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaR2)); 
	RooFormulaVar sigmaL3_scaled("sigmaL3_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaL3)); 
	RooFormulaVar sigmaR3_scaled("sigmaR3_scaled","@0 * @1", RooArgSet(scale_sigma_partReco,sigmaR3)); 

	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", DTF_Bs_M, mean1_scaled, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", DTF_Bs_M, mean2_scaled, sigmaL2,sigmaR2_scaled);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", DTF_Bs_M, mean3_scaled, sigmaL3,sigmaR3_scaled);
	RooAddPdf* bkg_partReco_alt1= new RooAddPdf("bkg_partReco_alt1", "bkg_partReco_alt1", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2),kTRUE);

	///alternative modelling using RooHILLdini & RooHORNSdini 
	if(altPartBkg) bkg_partReco_params = fitPartRecoBkgShapeHILLHORN();

	RooRealVar a("a","a", bkg_partReco_params[0]);
	RooRealVar b("b","b", bkg_partReco_params[1]);
	RooRealVar a_HORNS("a_HORNS","a_HORNS", bkg_partReco_params[7]);
	RooRealVar b_HORNS("b_HORNS","b_HORNS", bkg_partReco_params[8]);
	RooRealVar csi("csi","csi", bkg_partReco_params[2]);
	RooRealVar csi_HORNS("csi_HORNS","csi_HORNS", bkg_partReco_params[9]);
	RooRealVar shift("shift","shift", bkg_partReco_params[3]);
	RooRealVar sigma_HILL("sigma_HILL", "sigma_HILL", bkg_partReco_params[4]);
	RooRealVar sigma_HORNS("sigma_HORNS", "sigma_HORNS", bkg_partReco_params[10]);
	RooRealVar ratio_sigma("ratio_sigma", "ratio_sigma",bkg_partReco_params[5]);
	RooRealVar fraction_sigma("fraction_sigma", "fraction_sigma", bkg_partReco_params[6]);

	RooHILLdini RooHILLBkgShape("RooHILLBkgShape", "RooHILLBkgShape", DTF_Bs_M, a, b, csi, shift, sigma_HILL, ratio_sigma, fraction_sigma );
	RooHORNSdini RooHORNSBkgShape("RooHORNSBkgShape", "RooHORNSBkgShape", DTF_Bs_M, a_HORNS, b_HORNS, csi_HORNS, shift, sigma_HORNS, ratio_sigma, fraction_sigma );
	RooAddPdf* bkg_partReco_alt2= new RooAddPdf("bkg_partReco_alt2", "bkg_partReco_alt2", RooArgList(RooHILLBkgShape, RooHORNSBkgShape), RooArgList(f_1),kTRUE);
	
	/// used Part.reco bkg
	RooAddPdf* bkg_partReco;
	if(!altPartBkg)bkg_partReco = bkg_partReco_alt1;
	else bkg_partReco = bkg_partReco_alt2;
	bkg_partReco->SetName("bkg_partReco");
	bkg_partReco->SetTitle("bkg_partReco");

        if(ignorePartRecoBkg){
                RooArgSet* fitParamsPartRecoBkg = bkg_partReco->getParameters(data);
                RooFIter iterat = fitParamsPartRecoBkg->fwdIterator();  
                RooAbsArg * next = 0;
                while(next=iterat.next()) ((RooRealVar*)next)->setConstant();
        }

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_sig_B0("n_sig_B0", "n_sig_B0", data->numEntries()/10., 10., data->numEntries());
	RooRealVar n_exp_bkg("n_exp_bkg", "n_exp_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()/5., 50., data->numEntries() );
        if(ignorePartRecoBkg){
                n_partReco_bkg.setVal(0.);
                n_partReco_bkg.setConstant();
        }

	RooAddPdf* pdf;
	if(useB0 == 1) pdf = new RooAddPdf("pdf", "pdf", RooArgList(signal, bkg_exp, *bkg_partReco, signal_B0), RooArgList(n_sig, n_exp_bkg, n_partReco_bkg, n_sig_B0));
	if(useB0 == 0) pdf = new RooAddPdf("pdf", "pdf", RooArgList(signal, bkg_exp, *bkg_partReco), RooArgList(n_sig, n_exp_bkg, n_partReco_bkg));

	/// Generate simultaneous pdf out of prototype pdf 
  	RooSimPdfBuilder* mgr = new RooSimPdfBuilder(*pdf) ;
  	RooArgSet* config = mgr->createProtoBuildConfig() ;
  	config->setStringValue("physModels","pdf") ;
  	if(useTriggerCat)config->setStringValue("splitCats" ,"run Ds_finalState_mod TriggerCat") ;
	else config->setStringValue("splitCats" ,"year Ds_finalState_mod") ;

	if(useTriggerCat){
		if (useB0 == 0) config->setStringValue("pdf", "run            : scale_sigma, scale_mean "
					"Ds_finalState_mod :  exp_par "  
					"run,Ds_finalState_mod,TriggerCat : n_sig, n_exp_bkg, n_partReco_bkg") ;  
	
		if (useB0 == 1) config->setStringValue("pdf", "run            : scale_mean "
					"run,TriggerCat :	scale_sigma "
					"Ds_finalState_mod :  exp_par "  
					"run,Ds_finalState_mod,TriggerCat : n_sig, n_exp_bkg, n_partReco_bkg, n_sig_B0") ; 
	}
	else {
		if (useB0 == 0) config->setStringValue("pdf", "year            : scale_sigma, scale_mean "
					"Ds_finalState :  exp_par "  
					"year,Ds_finalState_mod : n_sig, n_exp_bkg, n_partReco_bkg") ;  
	
		if (useB0 == 1) config->setStringValue("pdf", "year            : scale_sigma, scale_mean "
					"Ds_finalState_mod :  exp_par "  
					"year,Ds_finalState : n_sig, n_exp_bkg, n_partReco_bkg, n_sig_B0") ; 
	}
	
  	RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
	simPdf->Print("v") ;
	RooArgSet* fitParams = simPdf->getParameters(data);

	/// Fix Exp from Sidebands
	if(useTriggerCat){		
		for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			/// Get data slice
			RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat == TriggerCat::" + str_trigger[k]);
			/// Fit
			//bkg_exp.fitTo(*data_slice,Save(kTRUE),Range(5500.,5700.));
			/// Fix parameters
			//((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setVal(exp_par.getVal());
			//if(fixExpBkgFromSidebands)((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setConstant();
	
			/// Set start values for yields
			((RooRealVar*) fitParams->find("n_sig_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/2.);
			((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/2.);
			if(!ignorePartRecoBkg)((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/2.);
			//if(!ignorePartRecoBkg)((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setRange(data_slice->numEntries()*0.02,data_slice->numEntries()*0.3);
			if(useB0)((RooRealVar*) fitParams->find("n_sig_B0_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()*0.01);
			if(useB0)((RooRealVar*) fitParams->find("n_sig_B0_{"+str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}"))->setRange(data_slice->numEntries()*0.005,data_slice->numEntries()*0.1);
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
			/// Get data slice
			RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]);
			/// Fit
			bkg_exp.fitTo(*data_slice,Save(kTRUE),Range(5500.,5700.));
			/// Fix parameters
			//((RooRealVar*) fitParams->find("exp_par_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(exp_par.getVal());
			//if(fixExpBkgFromSidebands)((RooRealVar*) fitParams->find("exp_par_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
			((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setVal(exp_par.getVal());
			if(fixExpBkgFromSidebands)((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setConstant();
	
			/// Set start values for yields
			((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/2.);
			((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/2.);
			if(!ignorePartRecoBkg)((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/10.);
			if(useB0)((RooRealVar*) fitParams->find("n_sig_B0_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/20.);
		}
	}

	/// Perform fit
	fitParams->Print("v") ;
	RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
	cout << "result is --------------- "<<endl;
	result->Print();

	/// Plot combined data and fit
	TCanvas* c = new TCanvas();
	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
	data->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
	simPdf->plotOn(frame,Name("pdf"),ProjWData(year,*data),LineColor(kBlue+1),LineWidth(3));

	TLegend leg(0.6,0.4,0.9,0.9,"");
    	leg.SetLineStyle(0);
    	leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.05);
	leg.SetTextAlign(12);
	leg.AddEntry(frame->findObject("data"),"LHCb data","ep");
	leg.AddEntry(frame->findObject("pdf"),"Fit","l");

	/// Plot components 
	/// argh fuck RooFit, does it really need to be so complicated ? 
	TString last_name_signal,last_name_signal_B0,last_name_exp_bkg,last_name_partReco_bkg;
	
	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}");
			RooArgList pdf_slice_comp( pdf_slice->pdfList());
	
			TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_signal_B0("signal_B0_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_exp_bkg("exp_bkg_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
	
			if(i==0 && j == 0 && k == 0){
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			}
			else{
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_B0));
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
			}
	
			if(i== str_run.size()-1 && j == str_Ds.size() -1 && k == str_trigger.size() -1) {
	
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
	
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),FillColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0),DrawOption("F"),FillStyle(3353),LineColor(kGreen+1));
				if(useB0)pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),LineColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
				
				leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s}#rightarrowD_{s}^{-}#pi^{+}#pi^{+}#pi^{-}}","f");
				if(useB0)leg.AddEntry(frame->findObject(name_signal_B0),"#font[132]{B^{0}#rightarrowD_{s}^{-}#pi^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
				if(!ignorePartRecoBkg)leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
			}
			else {
				last_name_signal = name_signal;
				if(useB0)last_name_signal_B0 = name_signal_B0;
				last_name_exp_bkg = name_exp_bkg;
				last_name_partReco_bkg = name_partReco_bkg;
			}
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
			RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_year[i] + ";" + str_Ds[j] + "}");
			RooArgList pdf_slice_comp( pdf_slice->pdfList());
	
			TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_signal_B0("signal_B0_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_exp_bkg("exp_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
	
			if(i==0 && j == 0){
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			}
			else{
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_B0));
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
			}
	
			if(i== str_year.size()-1 && j == str_Ds.size() -1) {
	
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
	
				if(useB0)pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),FillColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0),DrawOption("F"),FillStyle(3353),LineColor(kGreen+1));
				if(useB0)pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),LineColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
				
				leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s}#rightarrowD_{s}^{-}#pi^{+}#pi^{+}#pi^{-}}","f");
				if(useB0)leg.AddEntry(frame->findObject(name_signal_B0),"#font[132]{B^{0}#rightarrowD_{s}^{-}#pi^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
				if(!ignorePartRecoBkg)leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
			}
			else {
				last_name_signal = name_signal;
				if(useB0)last_name_signal_B0 = name_signal_B0;
				last_name_exp_bkg = name_exp_bkg;
				last_name_partReco_bkg = name_partReco_bkg;
			}
		}
	}

	if(useTriggerCat)simPdf->plotOn(frame,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));
	else simPdf->plotOn(frame,Name("pdf"),ProjWData(year,*data),LineColor(kBlue+1),LineWidth(3));
	frame->Draw();
	leg.Draw();
	c->Print("eps/norm.eps");
        //if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm.pdf");

	double chi2 = 0.;
	double covmatr = result->covQual();
	double edm = result->edm();

	TCanvas* canvas = new TCanvas();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);
        
        double max = 5.0 ;
        double min = -5.0 ;
        double rangeX = max-min;
        double zero = max/rangeX;
        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(max);
        graph->SetMinimum(min);
        graph->SetPoint(1,min_MM,0);
        graph->SetPoint(2,max_MM,0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(max);
        graph2->SetMinimum(min);
        graph2->SetPoint(1,min_MM,-3);
        graph2->SetPoint(2,max_MM,-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(max);
        graph3->SetMinimum(min);
        graph3->SetPoint(1,min_MM,3);
        graph3->SetPoint(2,max_MM,3);
        graph3->SetLineColor(kRed);
       
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame->GetYaxis()->SetRangeUser(0.01,frame->GetMaximum()*1.);
        frame->Draw();
        leg.Draw();
        
        canvas->cd();
        TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2->SetBorderMode(0);
        pad2->SetBorderSize(-1);
        pad2->SetFillStyle(0);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        RooPlot* frame_p = DTF_Bs_M.frame();
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
        frame_p->GetXaxis()->SetTitle("m(D_{s}^{-}#pi^{+}#pi^{+}#pi^{-}) [MeV/c^{2}]");
        
        RooHist* hpull  = frame->pullHist("data","pdf");
        frame_p->addPlotable(hpull,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->SaveAs("eps/norm_pull.eps");
        if(updateAnaNotePlots && !fitPreselected)canvas->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_pull.pdf");
        if(updateAnaNotePlots && fitPreselected)canvas->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_preselected_pull.pdf");

	/// Output file
	TFile *output;
	TTree* out_tree = 0;
	double sw;
	Int_t t_year, t_run, t_Ds_finalState, t_TriggerCat;
	TBranch *b_sw, *b_w, *b_year, *b_run, *b_Ds_finalState, *b_TriggerCat;

	if(sWeight){
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);
		tree->SetBranchStatus("weight",0);

		if(fitPreselected)out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) ).c_str() );
		else out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT).c_str() );

    		b_sw = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");
    		b_w = out_tree->Branch("weight", &sw, "weight/D");

        	out_tree->SetBranchAddress("year", &t_year, &b_year);
        	out_tree->SetBranchAddress("run", &t_run, &b_run);
        	out_tree->SetBranchAddress("Ds_finalState_mod", &t_Ds_finalState, &b_Ds_finalState);
        	out_tree->SetBranchAddress("TriggerCat", &t_TriggerCat, &b_TriggerCat);
		if(out_tree->GetEntries() != data->numEntries()) {
			cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
			cout << out_tree->GetEntries() << endl;
			cout << data->numEntries() << endl;
			throw "ERROR";
		} 
	}
	double weights[(int)data->numEntries()];

	/// Calculate total signal yield
	double yield = 0.;
	vector< double> signal_yields;
	vector< double> signal_yieldsB0;
	vector< double> signal_yieldsB0_err;
	vector< double> partReco_yields;
	vector< double> signal_params;
	vector< double> signal_yields_err;
	vector< double> partReco_yields_err;
	vector< double> expBkg_yields;
	vector< double> expBkg_yields_err;


	/// Loop over pdf slices
	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++){
	
			TCanvas* c1 = new TCanvas();
			TLatex* lhcbtext = new TLatex();
			lhcbtext->SetTextFont(132);
			lhcbtext->SetTextColor(1);
			lhcbtext->SetTextSize(0.07);
			lhcbtext->SetTextAlign(13);
			lhcbtext->SetNDC(1);
	
			for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
				/// Get pdf slice
				RooAbsPdf* pdf_slice = simPdf->getPdf("{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}");
	
				/// Plot data and pdf slices
				frame=DTF_Bs_M.frame();
				data->plotOn(frame,Name("data_slice2"),Cut("run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat==TriggerCat::" + str_trigger[k]),MarkerSize(1),Binning(nBins));
				pdf_slice->plotOn(frame,LineColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				if(useB0)pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				if(useB0)pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+3),Components("signal_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
				//data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				frame->Draw();
				chi2 += frame->chiSquare("pdf_slice2","data_slice2");
				cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;	

				/// Print label
				TString label = "LHCb " + str_run[i];
				label.ReplaceAll("1","I");
				label.ReplaceAll("2","II");
				lhcbtext->SetTextFont(22);
				lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("Run","Run-"));
				if(str_Ds[j]=="phipi")label = "D_{s}^{-}#rightarrow #phi^{0}(1020)#pi^{-}";
				else if(str_Ds[j]=="KsK")label = "D_{s}^{-}#rightarrow K^{*0}(892)K^{-}";
				//else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-}#rightarrow (K^{+}K^{-}#pi^{-})_{NR}";
				else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-}#rightarrow (K^{-}h^{+}#pi^{-})";
				else if(str_Ds[j]=="pipipi")label = "D_{s}^{-}#rightarrow #pi^{+}#pi^{-}#pi^{-}";
				lhcbtext->SetTextFont(132);
				lhcbtext->DrawLatex(0.6,0.78,label);
				if(str_trigger[k]=="t0")lhcbtext->DrawLatex(0.6,0.68,"L0-TOS");
				else lhcbtext->DrawLatex(0.6,0.68,"L0-TIS");
				c1->Print("eps/norm_" + str_run[i] + "_" + str_Ds[j]+ "_" + str_trigger[k]  + ".eps");
				if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k]  + ".pdf");
				hpull = frame->pullHist("data_slice2","pdf_slice2") ;
				frame= DTF_Bs_M.frame();
				frame->addPlotable(hpull,"P") ;
				frame->Draw();
				c1->Print("eps/norm_pull_" + str_run[i] + "_" + str_Ds[j]+ "_" + str_trigger[k]  + ".eps");
				if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_pull_" + str_run[i] + "_" + str_Ds[j]+ "_" + str_trigger[k]  + ".pdf");
	
				/// Get signal yield
				yield += ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
				signal_yields.push_back(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal());
				signal_yields_err.push_back(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError());
				if(useB0)signal_yieldsB0.push_back(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal());
				if(useB0)signal_yieldsB0_err.push_back(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError());
				partReco_yields.push_back(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal());
				partReco_yields_err.push_back(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError());
				expBkg_yields.push_back(((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal());
				expBkg_yields_err.push_back(((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError());
	
				/// Calculate sWeights
				if(sWeight){
					RooArgList yield_list(
					*((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")),
					*((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))
					);
					if(useB0)yield_list.add(*((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")));
					if(!ignorePartRecoBkg)yield_list.add(*((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")));
					RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat==TriggerCat::" + str_trigger[k]);
					pdf_slice->Print();
					SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
		
					/// Plot the sWeight distributions as a function of mass
					TH2 * swHist = (TH2*)data_slice->createHistogram("Bs_DTF_MM,n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}" + "_sw");
					swHist->GetYaxis()->SetTitle("Signal sWeights");
					swHist->Draw();
					c1->Print("eps/norm_sweight_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".eps");
	
					/// Save sWeights
					/// Messy and dangerous hack but works for now
					int n_ij = 0;  /// labels entry number of data slice
					for(int n = 0; n < out_tree->GetEntries(); n++){
						b_run->GetEntry(n);
						b_Ds_finalState->GetEntry(n);
						b_TriggerCat->GetEntry(n);
						if(t_run == i+1 && t_Ds_finalState == j && t_TriggerCat == k){
							weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}" + "_sw");
							n_ij++;
						}
					}	
				}
			}
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++){
	
			TCanvas* c1 = new TCanvas();
			TLatex* lhcbtext = new TLatex();
			lhcbtext->SetTextFont(132);
			lhcbtext->SetTextColor(1);
			lhcbtext->SetTextSize(0.07);
			lhcbtext->SetTextAlign(13);
			lhcbtext->SetNDC(1);
	
			signal_params.push_back(((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->getVal());
			signal_params.push_back(((RooRealVar*) fitParams->find("scale_mean_"+str_year[i]))->getVal());

			// Doesn't really work
			/*
			frame= DTF_Bs_M.frame();
			data->plotOn(frame,Name("data_slice"),Cut("year==year::" + str_year[i]),MarkerSize(1),Binning(nBins));
			simPdf->plotOn(frame,Name("pdf_slice"),Slice(year,str_year[i]),ProjWData(year,*data),LineColor(kBlue),LineWidth(3));
			frame->Draw();
			c->Print("eps/norm_" + str_year[i] + ".eps");
			hpull = frame->pullHist("data_slice","pdf_slice") ;
			frame= DTF_Bs_M.frame();
			frame->addPlotable(hpull,"P") ;
			frame->Draw();
			c->Print("eps/norm_pull_" + str_year[i] + ".eps");
			*/
			for(int j=0; j<str_Ds.size(); j++){
				/// Get pdf slice
				RooAbsPdf* pdf_slice = simPdf->getPdf("{" + str_year[i] + ";" + str_Ds[j] + "}");
	
				/// Plot data and pdf slices
				frame=DTF_Bs_M.frame();
				data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				pdf_slice->plotOn(frame,LineColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				if(useB0)pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_B0_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				if(useB0)pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+3),Components("signal_B0_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
				//data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				frame->Draw();
				chi2 += frame->chiSquare("pdf_slice2","data_slice2");
				cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;	

				/// Print label
				TString label = "LHCb data " + str_year[i];
				lhcbtext->SetTextFont(22);
				lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("y",""));
				if(str_Ds[j]=="phipi")label = "D_{s}^{-}#rightarrow#phi^{0}(1020) #pi^{-}";
				else if(str_Ds[j]=="KsK")label = "D_{s}^{-}#rightarrowK^{*0}(892) K^{-}";
				else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-}#rightarrow(K^{+} K^{-} #pi^{-})_{NR}";
				else if(str_Ds[j]=="pipipi")label = "D_{s}^{-}#rightarrow#pi^{+} #pi^{-} #pi^{-}";
				lhcbtext->SetTextFont(132);
				lhcbtext->DrawLatex(0.6,0.78,label);
				c1->Print("eps/norm_" + str_year[i] + "_" + str_Ds[j] + ".eps");
				if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_" + str_year[i] + "_" + str_Ds[j] + ".pdf");
				hpull = frame->pullHist("data_slice2","pdf_slice2") ;
				frame= DTF_Bs_M.frame();
				frame->addPlotable(hpull,"P") ;
				frame->Draw();
				c1->Print("eps/norm_pull_" + str_year[i] + "_" + str_Ds[j] + ".eps");
				if(updateAnaNotePlots)c1->Print("../../../../../TD-AnaNote/latex/figs/MassFit/norm_pull_" + str_year[i] + "_" + str_Ds[j] + ".pdf");
	
				/// Get signal yield
				yield += ((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal();
				signal_yields.push_back(((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal());
				signal_yields_err.push_back(((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getError());
				if(useB0)signal_yieldsB0.push_back(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_year[i] + ";" + str_Ds[j] +  "}"))->getVal());
				if(useB0)signal_yieldsB0_err.push_back(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_year[i] + ";" + str_Ds[j] +  "}"))->getError());
				partReco_yields.push_back(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal());
				partReco_yields_err.push_back(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getError());
				expBkg_yields.push_back(((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal());
				expBkg_yields_err.push_back(((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getError());
	
				/// Calculate sWeights
				if(sWeight){
					RooArgList yield_list(
					*((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}")),
					*((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))
					);
					if(useB0)yield_list.add(*((RooRealVar*) fitParams->find("n_sig_B0_{"+str_year[i] + ";" + str_Ds[j] + "}")));
					if(!ignorePartRecoBkg)yield_list.add(*((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}")));
					RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]);
					pdf_slice->Print();
					SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
		
					/// Plot the sWeight distributions as a function of mass
					TH2 * swHist = (TH2*)data_slice->createHistogram("Bs_DTF_MM,n_sig_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
					swHist->GetYaxis()->SetTitle("Signal sWeights");
					swHist->Draw();
					c1->Print("eps/norm_sweight_" + str_year[i] + "_" + str_Ds[j] + ".eps");
	
					/// Save sWeights
					/// Messy and dangerous hack but works for now
					int n_ij = 0;  /// labels entry number of data slice
					for(int n = 0; n < out_tree->GetEntries(); n++){
						b_year->GetEntry(n);
						b_Ds_finalState->GetEntry(n);
						if(A_is_in_B(anythingToString(t_year), (string) str_year[i]) && t_Ds_finalState == j){
							weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
							n_ij++;
						}
					}	
				}
			}
		}
	}

	if(sWeight){
		for(int n = 0; n < out_tree->GetEntries(); n++){
			sw = weights[n];
			b_sw->Fill();
			b_w->Fill();
		}
	 	out_tree->Write();
   		output->Close();
		cout << endl;
		cout << "Created file " << (string)outFileName  << endl << endl;
	}
	
	chi2 = chi2*nBins/(str_year.size()*str_Ds.size()*nBins-fitParams->selectByAttrib("Constant",kFALSE)->getSize());
	cout << endl <<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
	cout << "Total signal yield = " << yield << endl;

	/// Return fit params
	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++) for(int k=0; k<str_trigger.size(); k++)
			signal_params.push_back(((RooRealVar*) fitParams->find("scale_sigma_{" + str_run[i] + ";" + str_trigger[k] + "}"))->getVal());
		for(int i=0; i<str_run.size(); i++) signal_params.push_back(((RooRealVar*) fitParams->find("scale_mean_"+str_run[i]))->getVal());
	}
	signal_params.push_back(alpha.getVal());
	signal_params.push_back(beta.getVal());

	vector<double> partReco_params;
	partReco_params.push_back(mean1.getVal());
	partReco_params.push_back(mean2.getVal());
	partReco_params.push_back(mean3.getVal());
	partReco_params.push_back(sigmaL1.getVal());
	partReco_params.push_back(sigmaR1.getVal());
	partReco_params.push_back(sigmaL2.getVal());
	partReco_params.push_back(sigmaR2.getVal());
	partReco_params.push_back(sigmaL3.getVal());
	partReco_params.push_back(sigmaR3.getVal());
	partReco_params.push_back(f_1.getVal());
	partReco_params.push_back(f_2.getVal());
	partReco_params.push_back(bkg_partReco_params[11]);

	vector< vector <double> > return_vec;
	return_vec.push_back(signal_yields);
	return_vec.push_back(signal_params);
	return_vec.push_back(partReco_yields);
	return_vec.push_back(partReco_params);
	return_vec.push_back(signal_yields_err);
	return_vec.push_back(partReco_yields_err);
	return_vec.push_back(expBkg_yields);
	return_vec.push_back(expBkg_yields_err);
	return_vec.push_back(signal_yieldsB0);
	return_vec.push_back(signal_yieldsB0_err);

	if(updateAnaNotePlots){

		double yield_sig = 0;
		double yield_B0 = 0;
		double yield_partReco = 0;
		double yield_comb = 0;
		double yield_sig_err = 0;
		double yield_B0_err = 0;
		double yield_partReco_err = 0;
		double yield_comb_err = 0;
		for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			yield_sig += ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_sig_err += pow(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);

			yield_B0 += ((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_B0_err += pow(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);	

			yield_partReco += ((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_partReco_err += pow(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);

			yield_comb += ((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_comb_err += pow(((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);
		}

		ofstream datafile;
		datafile.open("../../../../../TD-AnaNote/latex/tables/MassFit/yields_norm.tex",std::ofstream::trunc);
		//datafile << "\\begin{table}[h]" << "\n";
		//datafile << "\\centering" << "\n";
		datafile << " \\begin{tabular}{l r }" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "Component & Yield\\" << " \\\\" << "\n";
		datafile << "\\hline" << "\n";
		datafile << std::fixed << std::setprecision(0) << "$B_s \\to D_s \\pi \\pi \\pi$" << " & "<< yield_sig << " $\\pm$ " << sqrt(yield_sig_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "$B^{0} \\to D_s \\pi \\pi \\pi$" << " & "<< yield_B0 << " $\\pm$ " << sqrt(yield_B0_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "Partially reconstructed bkg." << " & "<< yield_partReco << " $\\pm$ " << sqrt(yield_partReco_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "Combinatorial bkg." << " & "<< yield_comb << " $\\pm$ " << sqrt(yield_comb_err) << " \\\\" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "\\end{tabular}" << "\n";
		//datafile << "\\caption{Total signal and background yields for the $B_s \\to D_s \\pi \\pi \\pi$ sample.}" << "\n";
		datafile << "\\label{table:normYields}" << "\n";
		//datafile << "\\end{table}" << "\n";
		datafile.close();

		vector<double> yield_sig_Ds(4,0);
		vector<double> yield_sig_Ds_err(4,0);

		for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			yield_sig_Ds[j] += ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_sig_Ds_err[j] += pow(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);
		}

		ofstream datafileDs;
		datafile.open("../../../../../TD-AnaNote/latex/tables/MassFit/yields_normDs.tex",std::ofstream::trunc);
		//datafile << "\\begin{table}[h]" << "\n";
		//datafile << "\\centering" << "\n";
		datafile << " \\begin{tabular}{l r }" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "$D_s$ final state  & Signal yield\\" << " \\\\" << "\n";
		datafile << "\\hline" << "\n";
		for(int j=0; j<str_Ds.size(); j++){ 
			string caption;
			if(str_Ds[j]=="phipi")caption = "$D_{s}^{-} \\to \\phi^{0}(1020)\\pi^{-}$";
			else if(str_Ds[j]=="KsK")caption = "$D_{s}^{-}\\to K^{*0}(892)K^{-}$";
			else if(str_Ds[j]=="KKpi_NR")caption = "$D_{s}^{-}\\to (K^{-}h^{+}\\pi^{-})$";
			else if(str_Ds[j]=="pipipi")caption = "$D_{s}^{-}\\to \\pi^{+}\\pi^{-}\\pi^{-}$";

			datafile << std::fixed << std::setprecision(0) << caption << " & "<< yield_sig_Ds[j] << " $\\pm$ " << sqrt(yield_sig_Ds_err[j]) << " \\\\" << "\n";
		}
		datafile << "\\hline\\hline" << "\n";
		datafile << "\\end{tabular}" << "\n";
		//datafile << "\\caption{Signal yield for the different $D_s$ final states contributing to the $B_s \\to D_s \\pi \\pi \\pi$ sample.}" << "\n";
		datafile << "\\label{table:normYieldsDs}" << "\n";
		//datafile << "\\end{table}" << "\n";
		datafile.close();
	}


	if(file!=0)file->Close();
	return return_vec;
}

void fitSignal(){

	///Options
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> sWeight("sWeightSignal", 0);
	NamedParameter<int> nBins("nBins", 80);
	NamedParameter<double> min_MM("min_MM",5100.);
	//min_MM.setVal((double)min_MM-25.);
	NamedParameter<double> max_MM("max_MM",5700.);
	NamedParameter<string> cut_BDT("cut_BDT",(string)"");
	NamedParameter<int> fixMisIDyields("fixMisIDyields",1);
	NamedParameter<int> useNormScaleFactors("useNormScaleFactors",1);
	NamedParameter<int> useNormSignalShape("useNormSignalShape",1);
	NamedParameter<int> useTriggerCat("useTriggerCat", 0);
	NamedParameter<int> constrainCombBkgFromSS("constrainCombBkgFromSS",0);
	NamedParameter<int> optimizeBDT("optimizeBDT",0);
	NamedParameter<int> newTable("newTable",0);
	NamedParameter<int> useB0("useB0", 0);
        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
	NamedParameter<string> inFileName("inFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");
	NamedParameter<string> outFileName("outFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");

	///define categories
	RooCategory year("year","year") ;
	year.defineType("y11",11);
	year.defineType("y12",12);
	year.defineType("y15",15);
	year.defineType("y16",16);
	year.defineType("y17",17);
  	RooCategory run("run","run") ;
  	run.defineType("Run1",1) ;
  	run.defineType("Run2",2) ;
	RooCategory Ds_finalState_mod("Ds_finalState_mod","Ds_finalState_mod") ;
	Ds_finalState_mod.defineType("phipi",0);
	Ds_finalState_mod.defineType("KsK",1);
	Ds_finalState_mod.defineType("KKpi_NR",2);
	Ds_finalState_mod.defineType("pipipi",3);
	//Ds_finalState.defineType("Kpipi",4);
	RooCategory TriggerCat("TriggerCat","TriggerCat") ;
	TriggerCat.defineType("t0",0);
	TriggerCat.defineType("t1",1);

	///Load file
	TFile *file= new TFile(((string)inFileName).c_str());
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	tree->SetBranchStatus("BDTG",1);
	tree->SetBranchStatus("Ds_finalState*",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*_PIDK",1);
	tree->SetBranchStatus("pi_plus_isMuon",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);

        RooRealVar DTF_Bs_M("Bs_DTF_MM", "m(D_{s}^{-}K^{+}#pi^{+}#pi^{-})", min_MM, max_MM,"MeV/c^{2}");
        RooRealVar BDTG_response("BDTG", "BDTG", 0.);        
        RooRealVar K_plus_PIDK("K_plus_PIDK", "K_plus_PIDK", 0.);        
        RooRealVar pi_plus_PIDK("pi_plus_PIDK", "pi_plus_PIDK", 0.);        
        RooRealVar pi_minus_PIDK("pi_minus_PIDK", "pi_minus_PIDK", 0.);        

	RooArgList list =  RooArgList(DTF_Bs_M,BDTG_response,Ds_finalState_mod,year,TriggerCat,run,pi_minus_PIDK);
	RooDataSet*  data = new RooDataSet("data","data",tree,list,((string)cut_BDT + " && pi_minus_PIDK < 0").c_str() );	

	/// Fit normalization mode first
	vector< vector <double> > norm_paramSet = fitNorm();

	/// Signal Pdf
	vector<double> sig_params = fitSignalShape("signal");
	RooRealVar mean_MC("mean_MC", "mean_MC", sig_params[0]); 
	RooRealVar sigma_MC("sigma_MC", "sigma_MC", sig_params[1]);
	RooRealVar scale_mean("scale_mean", "scale_mean",1.,0.,2. ); 
	RooRealVar scale_sigma("scale_sigma", "scale_sigma", 1.,0.,2.);
	RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
	RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
	RooRealVar alpha("alpha", "alpha", sig_params[2],-1,1); 
	RooRealVar beta("beta", "beta", sig_params[3]); 

	RooJohnsonSU signal("signal","signal",DTF_Bs_M, mean,sigma,alpha,beta);

	/// B0 pdf
	RooFormulaVar mean_B0("mean_B0","@0 - @1", RooArgSet(mean,RooConst(87.33))); 
	RooRealVar scale_sigma_B0("scale_sigma_B0", "scale_sigma_B0", 1.,0.,2.);
	RooFormulaVar sigma_B0("sigma_B0","@0 * @1", RooArgSet(scale_sigma_B0,sigma)); 
	RooRealVar alpha_B0("alpha_B0", "alpha_B0", sig_params[2],-1,1); 
	RooRealVar beta_B0("beta_B0", "beta_B0", sig_params[3]); 

	RooJohnsonSU signal_B0("signal_B0","signal_B0",DTF_Bs_M, mean_B0,sigma_B0,alpha,beta);

	/// Combinatorial bkg pdf
	RooRealVar exp_par1("exp_par1","exp_par1",-1.6508e-03,-10.,10.);	
	RooRealVar exp_par("exp_par","exp_par",0.);	
	RooExponential bkg_exp1("bkg_exp1","bkg_exp1",DTF_Bs_M,exp_par);
	RooRealVar exp_par2("exp_par2","exp_par2",0,-10.,10.);	
	RooRealVar exp_par3("exp_par3","exp_par3",0,-1,1);	

	RooChebychev bkg_exp2("bkg_exp2","bkg_exp2",DTF_Bs_M,RooArgList(exp_par2,exp_par3));
	RooRealVar f_bkg_exp("f_bkg_exp", "f_bkg_exp", 0.);
	RooAddPdf bkg_exp("bkg_exp", "bkg_exp", RooArgList(bkg_exp1, bkg_exp2), RooArgList(f_bkg_exp));
	//if(!constrainCombBkgFromSS)bkg_exp.fitTo(*data,Save(kTRUE),Range(5500.,5700.));

	/// Part. reco bkg
	vector<double> bkg_partReco_params = norm_paramSet[3];
	RooRealVar mean1("mean1","mean1", bkg_partReco_params[0]);
	RooRealVar mean2("mean2","mean2", bkg_partReco_params[1]);
	RooRealVar mean3("mean3","mean3", bkg_partReco_params[2]);
	RooRealVar mean1Shifted("mean1Shifted","mean1Shifted", mean1.getVal() - 87.33 );
	RooRealVar mean2Shifted("mean2Shifted","mean2Shifted", mean2.getVal() - 87.33 );
	RooRealVar mean3Shifted("mean3Shifted","mean3Shifted", mean3.getVal() - 87.33 );
	RooRealVar sigmaL1("sigmaL1", "sigmaL1",  bkg_partReco_params[3]);
	RooRealVar sigmaR1("sigmaR1", "sigmaR1",  bkg_partReco_params[4]);
	RooRealVar sigmaL2("sigmaR1", "sigmaL2",  bkg_partReco_params[5]);
	RooRealVar sigmaR2("sigmaR2", "sigmaR2",  bkg_partReco_params[6]);
	RooRealVar sigmaL3("sigmaL3", "sigmaL3",  bkg_partReco_params[7]);
	RooRealVar sigmaR3("sigmaR3", "sigmaR3",  bkg_partReco_params[8]);
	RooRealVar f_1("f_1", "f_1", bkg_partReco_params[9]);
	RooRealVar f_2("f_2", "f_2", bkg_partReco_params[10]);
	double eff_bkg_partReco_inSigRange = bkg_partReco_params[11];

	RooBifurGauss BifGauss1_Bs("BifGauss1_Bs","BifGauss1_Bs", DTF_Bs_M, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2_Bs("BifGauss2_Bs","BifGauss2_Bs", DTF_Bs_M, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3_Bs("BifGauss3_Bs","BifGauss3_Bs", DTF_Bs_M, mean3, sigmaL3,sigmaR3);
	RooAddPdf bkg_partReco_Bs("bkg_partReco_Bs", "bkg_partReco_Bs", RooArgList(BifGauss1_Bs, BifGauss2_Bs, BifGauss3_Bs), RooArgList(f_1,f_2));

	RooBifurGauss BifGauss1_B0("BifGauss1_B0","BifGauss1_B0", DTF_Bs_M, mean1Shifted, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2_B0("BifGauss2_B0","BifGauss2_B0", DTF_Bs_M, mean2Shifted, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3_B0("BifGauss3_B0","BifGauss3_B0", DTF_Bs_M, mean3Shifted, sigmaL3,sigmaR3);
	RooAddPdf bkg_partReco_B0("bkg_partReco_B0", "bkg_partReco_B0", RooArgList(BifGauss1_B0, BifGauss2_B0, BifGauss3_B0), RooArgList(f_1,f_2));

	RooRealVar partReco_f("partReco_f", "partReco_f", 0.9,0,1);
	RooAddPdf bkg_partReco("bkg_partReco", "bkg_partReco", RooArgList(bkg_partReco_Bs, bkg_partReco_B0), RooArgList(partReco_f));

	/// MisID bkg
	double fake_prob_Ds = prepareMisIdBkgShape("Ds3pi");
	double fake_prob_Dstar = prepareMisIdBkgShape("Dstar3pi");
	vector<double> bkg_misID_Ds3pi_params = fitMisIdBkgShape_Ds3pi();
	
	RooRealVar misID_Ds3pi_mean1("misID_Ds3pi_mean1","misID_Ds3pi_mean1", bkg_misID_Ds3pi_params[0]);//,bkg_misID_Ds3pi_params[0]-50,bkg_misID_Ds3pi_params[0]+50);
	RooRealVar misID_Ds3pi_mean2("misID_Ds3pi_mean2","misID_Ds3pi_mean2", bkg_misID_Ds3pi_params[1]);
	RooRealVar misID_Ds3pi_a1("misID_Ds3pi_a1","misID_Ds3pi_a1",bkg_misID_Ds3pi_params[2]);
	RooRealVar misID_Ds3pi_a2("misID_Ds3pi_a2","misID_Ds3pi_a2",bkg_misID_Ds3pi_params[3]);
	RooRealVar misID_Ds3pi_n1("misID_Ds3pi_n1","misID_Ds3pi_n1",bkg_misID_Ds3pi_params[4]);
	RooRealVar misID_Ds3pi_n2("misID_Ds3pi_n2","misID_Ds3pi_n2",bkg_misID_Ds3pi_params[5]);
	RooRealVar misID_Ds3pi_sigma1("misID_Ds3pi_sigma1", "misID_Ds3pi_sigma1", bkg_misID_Ds3pi_params[6]);//,  bkg_misID_Ds3pi_params[6] * 0.5,  bkg_misID_Ds3pi_params[6]*1.5);
	RooRealVar misID_Ds3pi_sigma2("misID_Ds3pi_sigma2", "misID_Ds3pi_sigma2", bkg_misID_Ds3pi_params[7]);
	RooRealVar misID_Ds3pi_f("misID_Ds3pi_f", "misID_Ds3pi_f", bkg_misID_Ds3pi_params[8]);

	RooCBShape misID_Ds3pi_CB1("misID_Ds3pi_CB1", "misID_Ds3pi_CB1", DTF_Bs_M, misID_Ds3pi_mean1, misID_Ds3pi_sigma1, misID_Ds3pi_a1, misID_Ds3pi_n1);
	RooCBShape misID_Ds3pi_CB2("misID_Ds3pi_CB2", "misID_Ds3pi_CB2", DTF_Bs_M, misID_Ds3pi_mean2, misID_Ds3pi_sigma2, misID_Ds3pi_a2, misID_Ds3pi_n2);
	RooAddPdf bkg_misID_Ds3pi("bkg_misID_Ds3pi", "bkg_misID_Ds3pi", RooArgList(misID_Ds3pi_CB1, misID_Ds3pi_CB2), RooArgList(misID_Ds3pi_f));

	vector<double> bkg_misID_Dsstar3pi_params = fitMisIdBkgShape_Dsstar3pi();
	RooRealVar misID_Dsstar3pi_mean1("misID_Dsstar3pi_mean1","misID_Dsstar3pi_mean1", bkg_misID_Dsstar3pi_params[0]);
	RooRealVar misID_Dsstar3pi_mean2("misID_Dsstar3pi_mean2","misID_Dsstar3pi_mean2", bkg_misID_Dsstar3pi_params[1]);
	RooRealVar misID_Dsstar3pi_a1("misID_Dsstar3pi_a1","misID_Dstars3pi_a1",bkg_misID_Dsstar3pi_params[2]);
	RooRealVar misID_Dsstar3pi_a2("misID_Dsstar3pi_a2","misID_Dsstar3pi_a2",bkg_misID_Dsstar3pi_params[3]);
	RooRealVar misID_Dsstar3pi_n1("misID_Dsstar3pi_n1","misID_Dsstar3pi_n1",bkg_misID_Dsstar3pi_params[4]);
	RooRealVar misID_Dsstar3pi_n2("misID_Dsstar3pi_n2","misID_Dsstar3pi_n2",bkg_misID_Dsstar3pi_params[5]);
	RooRealVar misID_Dsstar3pi_sigma1("misID_Dsstar3pi_sigma1", "misID_Dsstar3pi_sigma1", bkg_misID_Dsstar3pi_params[6]);
	RooRealVar misID_Dsstar3pi_sigma2("misID_Dsstar3pi_sigma2", "misID_Dsstar3pi_sigma2", bkg_misID_Dsstar3pi_params[7]);
	RooRealVar misID_Dsstar3pi_f("misID_Dsstar3pi_f", "misID_Dsstar3pi_f", bkg_misID_Dsstar3pi_params[8]);

	RooCBShape misID_Dsstar3pi_CB1("misID_Dsstar3pi_CB1", "misID_Dsstar3pi_CB1", DTF_Bs_M, misID_Dsstar3pi_mean1, misID_Dsstar3pi_sigma1, misID_Dsstar3pi_a1, misID_Dsstar3pi_n1);
	RooCBShape misID_Dsstar3pi_CB2("misID_Dsstar3pi_CB2", "misID_Dsstar3pi_CB2", DTF_Bs_M, misID_Dsstar3pi_mean2, misID_Dsstar3pi_sigma2, misID_Dsstar3pi_a2, misID_Dsstar3pi_n2);
	RooAddPdf bkg_misID_Dsstar3pi("bkg_misID_Dsstar3pi", "bkg_misID_Dsstar3pi", RooArgList(misID_Dsstar3pi_CB1, misID_Dsstar3pi_CB2), RooArgList(misID_Dsstar3pi_f));

	RooRealVar misID_f("misID_f", "misID_f", 0.5,0.,1.);
	RooAddPdf bkg_misID("bkg_misID", "bkg_misID", RooArgList(bkg_misID_Ds3pi, bkg_misID_Dsstar3pi), RooArgList(misID_f));

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/4., 0., data->numEntries());
	RooRealVar n_sig_B0("n_sig_B0", "n_sig_B0", data->numEntries()/4., 0., data->numEntries());
	RooRealVar n_exp_bkg("n_exp_bkg", "n_exp_bkg", data->numEntries()/2., 0., data->numEntries());
	//RooRealVar n_partReco_Bs_bkg("n_partReco_Bs_bkg", "n_partReco_Bs_bkg", data->numEntries()/2., 0., data->numEntries());
	//RooRealVar n_partReco_B0_bkg("n_partReco_B0_bkg", "n_partReco_B0_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()/6., 100., data->numEntries());
	RooRealVar n_misID_bkg("n_misID_bkg", "n_misID_bkg", data->numEntries()*0.01, 0., data->numEntries()*0.25);

	RooAddPdf pdf("pdf", "pdf", RooArgList(signal, signal_B0, bkg_exp, bkg_partReco, bkg_misID), RooArgList(n_sig, n_sig_B0, n_exp_bkg, n_partReco_bkg, n_misID_bkg));

	/// Generate simultaneous pdf out of prototype pdf 
  	RooSimPdfBuilder* mgr = new RooSimPdfBuilder(pdf) ;
  	RooArgSet* config = mgr->createProtoBuildConfig() ;
  	config->setStringValue("physModels","pdf") ;
  	if(useTriggerCat)config->setStringValue("splitCats" ,"run Ds_finalState_mod TriggerCat") ;
	else config->setStringValue("splitCats" ,"year Ds_finalState_mod") ;

	if(useTriggerCat){
		config->setStringValue("pdf", "run            : scale_mean "
					"run,TriggerCat :	scale_sigma "
					"Ds_finalState_mod :  exp_par2, exp_par3 "  
					"run,Ds_finalState_mod,TriggerCat : n_sig, n_sig_B0, n_exp_bkg, n_partReco_bkg, n_misID_bkg, misID_f") ; 
	}
	else {
		config->setStringValue("pdf", "year            : scale_sigma, scale_mean "
					"Ds_finalState_mod :  exp_par "  
					"year,Ds_finalState_mod : n_sig, n_sig_B0, n_exp_bkg, n_partReco_bkg, n_misID_bkg, misID_f") ; 
	}

  	RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
	simPdf->Print("v") ;
	RooArgSet* fitParams = simPdf->getParameters(data);
  	fitParams->Print("v") ;

	/// Fix scale factors from normalization mode
	vector<double> scaleFactors = norm_paramSet[1];
	int scaleFactor_counter = 0;

	if(useTriggerCat){	
		for(int i=0; i<str_run.size(); i++) for(int k=0; k<str_trigger.size(); k++){
			((RooRealVar*) fitParams->find("scale_sigma_{" + str_run[i] + ";" + str_trigger[k] + "}"))->setVal( scaleFactors[scaleFactor_counter] );
			if(useNormScaleFactors)((RooRealVar*) fitParams->find("scale_sigma_{" + str_run[i] + ";" + str_trigger[k] + "}"))->setConstant();
			scaleFactor_counter ++;
		}
		for(int i=0; i<str_run.size(); i++) {
			((RooRealVar*) fitParams->find("scale_mean_"+str_run[i]))->setVal( scaleFactors[scaleFactor_counter] );
			if(useNormScaleFactors)((RooRealVar*) fitParams->find("scale_mean_"+str_run[i]))->setConstant();
			scaleFactor_counter ++;
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++){
			((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->setVal( scaleFactors[scaleFactor_counter] );
			if(useNormScaleFactors)((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->setConstant();
			scaleFactor_counter ++;
			((RooRealVar*) fitParams->find("scale_mean_"+str_year[i]))->setVal( scaleFactors[scaleFactor_counter] );
			if(useNormScaleFactors)((RooRealVar*) fitParams->find("scale_mean_"+str_year[i]))->setConstant();
			scaleFactor_counter ++;
		}
	}
	if(useNormSignalShape){
		alpha.setVal(scaleFactors[scaleFactors.size()-2]);
		beta.setVal(scaleFactors[scaleFactors.size()-1]);
		alpha.setConstant();
		beta.setConstant();
	}

	/// Fix misID yields 
	vector<double> Ds_yields = norm_paramSet[0];
	vector<double> Dstar_yields = norm_paramSet[2];
	
	int counter = 0;

	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_Ds.size(); j++) for(int k=0; k<str_trigger.size(); k++){
			double val = fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]/eff_bkg_partReco_inSigRange   ;
			((RooRealVar*) fitParams->find("n_misID_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal(val);
			if(fixMisIDyields)((RooRealVar*) fitParams->find("n_misID_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setConstant();
			else ((RooRealVar*) fitParams->find("n_misID_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setRange(val * 0.5 , val * 3.);
	
			((RooRealVar*) fitParams->find("misID_f_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal( fake_prob_Ds * Ds_yields[counter] / (fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]/eff_bkg_partReco_inSigRange   ));
			((RooRealVar*) fitParams->find("misID_f_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setConstant();
	
			/// Set start values for yields
			RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat==TriggerCat::" + str_trigger[k]);
			((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/4.);
			((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/4.);
			((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/4.);
			((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setRange(data_slice->numEntries()*0.03,data_slice->numEntries()*0.3);
			((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->setVal(data_slice->numEntries()/4.);
	
			counter++;
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
			double val = fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]/eff_bkg_partReco_inSigRange   ;
			((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(val);
			if(fixMisIDyields)((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
			else ((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setRange(val * 0. , val * 2.);
	
			((RooRealVar*) fitParams->find("misID_f_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal( fake_prob_Ds * Ds_yields[counter] / (fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]/eff_bkg_partReco_inSigRange   ));
			((RooRealVar*) fitParams->find("misID_f_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
	
			/// Set start values for yields
			RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]);
			((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/2.);
			((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/2.);
			((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/10.);
			((RooRealVar*) fitParams->find("n_sig_B0_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(data_slice->numEntries()/20.);
	
			counter++;
		}
	}


	/// Constrain Exp from SS line
	RooArgSet constr_set;
	if(constrainCombBkgFromSS){		
		TFile *file_SS= new TFile("/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root");
		TTree* tree_SS = (TTree*) file_SS->Get("DecayTree");	
		tree_SS->SetBranchStatus("*",0);
		tree_SS->SetBranchStatus("Bs_DTF_MM",1);
		tree_SS->SetBranchStatus("Ds_finalState*",1);
		tree_SS->SetBranchStatus("year",1);
		tree_SS->SetBranchStatus("TriggerCat",1);
		tree_SS->SetBranchStatus("run",1);      
		tree_SS->SetBranchStatus("BDTG",1);      	

		RooArgList list_SS =  RooArgList(DTF_Bs_M,Ds_finalState_mod,year,run,TriggerCat);
		RooDataSet*  data_SS = new RooDataSet("data_SS","data_SS",tree_SS,list,((string)cut_BDT).c_str() );

		for(int j=0; j<str_Ds.size(); j++){
			/// Get data slice
			RooDataSet* data_slice = new RooDataSet("data_SS_slice","data_SS_slice",data_SS,list,"Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]);
			/// Fit
			bkg_exp.fitTo(*data_slice,Save(kTRUE));
			/// Constrain parameters
			((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setVal(exp_par.getVal());
			//((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setConstant();
			RooGaussian* constr = new RooGaussian("constr_"+str_Ds[j],"constr_"+str_Ds[j],*((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j])),RooRealConstant::value(exp_par.getVal()),RooRealConstant::value(exp_par.getError()));
			constr_set.add(*constr);
			/// Plot
			TCanvas* c = new TCanvas();
			RooPlot* frame= DTF_Bs_M.frame();
			frame->SetTitle("");
			data_slice->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
			bkg_exp.plotOn(frame,Name("bkg_exp"),LineColor(kBlue+1),LineWidth(3));
			frame->Draw();
			c->Print("eps/comb_"+str_Ds[j]+".eps");
		}
	}
	RooProdPdf* constr_exp; 
	if(constrainCombBkgFromSS) constr_exp = new RooProdPdf("constr_exp","constr_exp",constr_set);

	/// Perform fit
	RooFitResult* result;
	if(constrainCombBkgFromSS)result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU),ExternalConstraints(*constr_exp));
	else result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
	cout << "result is --------------- "<<endl;
	result->Print("v");

	/// Plot combined data and fit
	TCanvas* c = new TCanvas();
	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
	data->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
	simPdf->plotOn(frame,Name("pdf"),ProjWData(year,*data),LineColor(kBlue+1),LineWidth(3));

	TLegend leg(0.6,0.4,0.9,0.9,"");
    	leg.SetLineStyle(0);
    	leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.05);
	leg.SetTextAlign(12);
	leg.AddEntry(frame->findObject("data"),"LHCb data","ep");
	leg.AddEntry(frame->findObject("pdf"),"Fit","l");

	/// Plot components 
	/// argh fuck RooFit, is this seriously so complicated ? 
	TString last_name_signal,last_name_signal_B0,last_name_exp_bkg,last_name_partReco_bkg,last_name_misID_bkg;

	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_Ds[j] + ";" + str_trigger[k] + "}");
			RooArgList pdf_slice_comp( pdf_slice->pdfList());
	
			TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_signal_B0("signal_B0_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_exp_bkg("exp_bkg_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
			TString name_misID_bkg("misID_bkg_"+ anythingToString(i)+ "_" + anythingToString(j)+ "_" + anythingToString(k));
	
			if(i==0 && j == 0  && k == 0){
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
					
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			}
			else{
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
				
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
		
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_misID_bkg));		
			}
	
			if(i== str_run.size()-1 && j == str_Ds.size() -1 && k == str_trigger.size() -1) {
	
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),FillColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg),DrawOption("F"),FillStyle(3344),LineColor(kMagenta+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),LineColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg));
	
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
	
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),FillColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0),DrawOption("F"),FillStyle(3335),LineColor(kGreen+3));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),LineColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
	
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),FillColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg),DrawOption("F"),FillStyle(3344),LineColor(kMagenta+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),LineColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg));
	
				leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s}#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_signal_B0),"#font[132]{B^{0}#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
				leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
				leg.AddEntry(frame->findObject(name_misID_bkg),"MisID bkg.","f");
			}
			else {
				last_name_signal = name_signal;
				last_name_signal_B0 = name_signal_B0;
				last_name_exp_bkg = name_exp_bkg;
				last_name_partReco_bkg = name_partReco_bkg;
				last_name_misID_bkg = name_misID_bkg;
			}
		}
	}

	else {
		for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
			RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_year[i] + ";" + str_Ds[j] + "}");
			RooArgList pdf_slice_comp( pdf_slice->pdfList());
	
			TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_signal_B0("signal_B0_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_exp_bkg("exp_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
			TString name_misID_bkg("misID_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
	
			if(i==0 && j == 0){
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
					
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
				
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			}
			else{
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
				
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
		
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_misID_bkg));		
			}
	
			if(i== str_year.size()-1 && j == str_Ds.size() -1) {
	
				pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[3])),ProjWData(year,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
	
				pdf_slice->plotOn(frame,Name(name_misID_bkg),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),FillColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg),DrawOption("F"),FillStyle(3344),LineColor(kMagenta+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),LineColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg));
	
				pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
	
				pdf_slice->plotOn(frame,Name(name_signal_B0),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),FillColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0),DrawOption("F"),FillStyle(3335),LineColor(kGreen+3));
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),LineColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_B0));
	
				pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
	
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),FillColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg),DrawOption("F"),FillStyle(3344),LineColor(kMagenta+3));			
				pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[4])),ProjWData(year,*data),LineColor(kMagenta+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_misID_bkg));
	
				leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s}#rightarrowD_{s}^{-} K^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_signal_B0),"#font[132]{B^{0}#rightarrowD_{s}^{-}K^{+}#pi^{+}#pi^{-}}","f");
				leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
				leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
				leg.AddEntry(frame->findObject(name_misID_bkg),"MisID bkg.","f");
			}
			else {
				last_name_signal = name_signal;
				last_name_signal_B0 = name_signal_B0;
				last_name_exp_bkg = name_exp_bkg;
				last_name_partReco_bkg = name_partReco_bkg;
				last_name_misID_bkg = name_misID_bkg;
			}
		}
	}

	if(useTriggerCat)simPdf->plotOn(frame,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));  //,VisualizeError(*result,1,kTRUE)
	else simPdf->plotOn(frame,Name("pdf"),ProjWData(year,*data),LineColor(kBlue+1),LineWidth(3));
	frame->Draw();
	leg.Draw();
	c->Print("eps/signal.eps");
	//if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal.pdf");

	double chi2 = 0.;
	double covmatr = result->covQual();
	double edm = result->edm();

	TCanvas* canvas = new TCanvas();
	canvas->cd();
	canvas->UseCurrentStyle();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);
        
        double max = 5.0 ;
        double min = -5.0 ;
        double rangeX = max-min;
        double zero = max/rangeX;
        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(max);
        graph->SetMinimum(min);
        graph->SetPoint(1,min_MM,0);
        graph->SetPoint(2,max_MM,0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(max);
        graph2->SetMinimum(min);
        graph2->SetPoint(1,min_MM,-3);
        graph2->SetPoint(2,max_MM,-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(max);
        graph3->SetMinimum(min);
        graph3->SetPoint(1,min_MM,3);
        graph3->SetPoint(2,max_MM,3);
        graph3->SetLineColor(kRed);
       
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame->GetYaxis()->SetRangeUser(0.01,frame->GetMaximum()*1.);
        frame->Draw();
        leg.Draw();
        
        canvas->cd();
        TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2->SetBorderMode(0);
        pad2->SetBorderSize(-1);
        pad2->SetFillStyle(0);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        RooPlot* frame_p = DTF_Bs_M.frame();
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
        frame_p->GetXaxis()->SetTitle("m(D_{s}^{-}K^{+}#pi^{+}#pi^{-}) [MeV/c^{2}]");
        
        RooHist* hpull  = frame->pullHist("data","pdf");
        frame_p->addPlotable(hpull,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->Print("eps/signal_pull.eps");
        if(updateAnaNotePlots)canvas->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal_pull.pdf");
	c->cd();

	/// Output file
	TFile *output;
	TTree* out_tree;
	double sw,sw_B0;
	Int_t t_year, t_run, t_Ds_finalState, t_TriggerCat;
	TBranch *b_sw, *b_w, *b_year, *b_run, *b_Ds_finalState, *b_TriggerCat;

	if(sWeight){
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);
		tree->SetBranchStatus("weight",0);

		out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT ).c_str() );
    		b_sw = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");
    		b_w = out_tree->Branch("weight", &sw_B0, "weight/D");

        	out_tree->SetBranchAddress("year", &t_year, &b_year);
        	out_tree->SetBranchAddress("run", &t_run, &b_run);
        	out_tree->SetBranchAddress("Ds_finalState_mod", &t_Ds_finalState, &b_Ds_finalState);
        	out_tree->SetBranchAddress("TriggerCat", &t_TriggerCat, &b_TriggerCat);
		if(out_tree->GetEntries() != data->numEntries()) {
			cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
			cout << out_tree->GetEntries() << endl;
			cout << data->numEntries() << endl;
			throw "ERROR";
		} 
	}
	double weights[(int)data->numEntries()];
	double weights_B0[(int)data->numEntries()];

	/// Calculate total signal yield
	double signal_yield = 0.;
	double comb_bkg_yield = 0.;
	DTF_Bs_M.setRange("signal_range",mean.getVal()-45.,mean.getVal()+45.);

	/// Loop over pdf slices
	if(useTriggerCat){
		for(int i=0; i<str_run.size(); i++){

			TLatex* lhcbtext = new TLatex();
			lhcbtext->SetTextFont(22);
			lhcbtext->SetTextColor(1);
			lhcbtext->SetTextSize(0.07);
			lhcbtext->SetTextAlign(13);
			lhcbtext->SetNDC(1);
		
			for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
				/// Get pdf slice
				RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}");
				RooArgList pdf_slice_comp( pdf_slice->pdfList());
				pdf_slice->Print();
				/// Plot data and pdf slices
				frame=DTF_Bs_M.frame();
				data->plotOn(frame,Name("data_slice2"),Cut("run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat==TriggerCat::" + str_trigger[k]),MarkerSize(1),Binning(nBins));
				pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kMagenta+3),Components("bkg_misID_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3344),FillColor(kMagenta+3),Components("bkg_misID_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+3),Components("signal_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				//data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				frame->Draw();
				chi2 += frame->chiSquare("pdf_slice2","data_slice2");
				cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;	
	
				/// Print label
				TString label = "LHCb " + str_run[i];
				label.ReplaceAll("1","I");
				label.ReplaceAll("2","II");
				lhcbtext->SetTextFont(22);
				lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("Run","Run-"));
				if(str_Ds[j]=="phipi")label = "D_{s}^{-}#rightarrow #phi^{0}(1020)#pi^{-}";
				else if(str_Ds[j]=="KsK")label = "D_{s}^{-}#rightarrow K^{*0}(892)K^{-}";
				else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-}#rightarrow K^{-}h^{+}#pi^{-}";
	 			//else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-}#rightarrow (K^{+} K^{-}#pi^{-})_{NR}";
				else if(str_Ds[j]=="pipipi")label = "D_{s}^{-}#rightarrow #pi^{+}#pi^{-}#pi^{-}";
				lhcbtext->SetTextFont(132);
				lhcbtext->DrawLatex(0.6,0.78,label);
				if(str_trigger[k]=="t0")lhcbtext->DrawLatex(0.6,0.68,"L0-TOS");
				else lhcbtext->DrawLatex(0.6,0.68,"L0-TIS");
				c->Print("eps/signal_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".eps");
				if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".pdf");
				hpull = frame->pullHist("data_slice2","pdf_slice2") ;
				frame= DTF_Bs_M.frame();
				frame->addPlotable(hpull,"P") ;
				frame->Draw();
				c->Print("eps/signal_pull_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".eps");
				if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal_pull_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".pdf");
	
				/// Get signal yield
				signal_yield += ((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
				RooAbsPdf* pdf_slice_comb_bkg = (RooAbsPdf*) pdf_slice_comp.find("bkg_exp_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}");
				comb_bkg_yield += pdf_slice_comb_bkg->createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("signal_range"))->getVal()*((RooRealVar*) fitParams->find("n_exp_bkg_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
	
				/// Calculate sWeights
				if(sWeight){
					RooArgList yield_list(
					*((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")),
					*((RooRealVar*) fitParams->find("n_sig_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")),
					*((RooRealVar*) fitParams->find("n_exp_bkg_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")),
					*((RooRealVar*) fitParams->find("n_partReco_bkg_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}")),
					*((RooRealVar*) fitParams->find("n_misID_bkg_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))
					);
					RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j] + " && TriggerCat==TriggerCat::" + str_trigger[k]);
					SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
		
					/// Plot the sWeight distributions as a function of mass
					TH2 * swHist = (TH2*)data_slice->createHistogram("Bs_DTF_MM,n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}" + "_sw");
					swHist->GetYaxis()->SetTitle("Signal sWeights");
					swHist->Draw();
					c->Print("eps/signal_sweight_" + str_run[i] + "_" + str_Ds[j] + "_" + str_trigger[k] + ".eps");
	
					/// Save sWeights
					/// Messy and dangerous hack but works for now
					int n_ij = 0;  /// labels entry number of data slice
					for(int n = 0; n < out_tree->GetEntries(); n++){
						b_run->GetEntry(n);
						b_Ds_finalState->GetEntry(n);
						b_TriggerCat->GetEntry(n);
						if(t_run == i+1 && t_Ds_finalState == j && t_TriggerCat == k){
							weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}" + "_sw");
							weights_B0[n] = sPlot.GetSWeight(n_ij,"n_sig_B0_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}" + "_sw");
							n_ij++;
						}
					}	
				}
			}
		}
	}
	else {
		for(int i=0; i<str_year.size(); i++){
	
			TLatex* lhcbtext = new TLatex();
			lhcbtext->SetTextFont(22);
			lhcbtext->SetTextColor(1);
			lhcbtext->SetTextSize(0.07);
			lhcbtext->SetTextAlign(13);
			lhcbtext->SetNDC(1);
	
			frame= DTF_Bs_M.frame();
			data->plotOn(frame,Name("data_slice"),Cut("year==year::" + str_year[i]),MarkerSize(1),Binning(nBins));
			simPdf->plotOn(frame,Name("pdf_slice"),Slice(year,str_year[i]),ProjWData(year,*data),LineColor(kBlue),LineWidth(3));
			frame->Draw();
			c->Print("eps/signal_" + str_year[i] + ".eps");
	
			hpull = frame->pullHist("data_slice","pdf_slice") ;
			frame= DTF_Bs_M.frame();
			frame->addPlotable(hpull,"P") ;
			frame->Draw();
			c->Print("eps/signal_pull_" + str_year[i] + ".eps");
	
			for(int j=0; j<str_Ds.size(); j++){
				/// Get pdf slice
				RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_year[i] + ";" + str_Ds[j] + "}");
				RooArgList pdf_slice_comp( pdf_slice->pdfList());
				pdf_slice->Print();
				/// Plot data and pdf slices
				frame=DTF_Bs_M.frame();
				data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kMagenta+3),Components("bkg_misID_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3344),FillColor(kMagenta+3),Components("bkg_misID_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_B0_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+3),Components("signal_B0_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
	
				pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
				//data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
				frame->Draw();
				chi2 += frame->chiSquare("pdf_slice2","data_slice2");
				cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;	
	
				/// Print label
				TString label = "LHCb data " + str_year[i];
				lhcbtext->SetTextFont(22);
				lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("y",""));
				if(str_Ds[j]=="phipi")label = "D_{s}^{-} #rightarrow #phi^{0}(1020) #pi^{-}";
				else if(str_Ds[j]=="KsK")label = "D_{s}^{-} #rightarrow K^{*0}(892) K^{-}";
				else if(str_Ds[j]=="KKpi_NR")label = "D_{s}^{-} #rightarrow (K^{+} K^{-} #pi^{-})_{NR}";
				else if(str_Ds[j]=="pipipi")label = "D_{s}^{-} #rightarrow #pi^{+} #pi^{-} #pi^{-}";
				lhcbtext->SetTextFont(132);
				lhcbtext->DrawLatex(0.6,0.78,label);
				c->Print("eps/signal_" + str_year[i] + "_" + str_Ds[j] + ".eps");
				if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal_" + str_year[i] + "_" + str_Ds[j] + ".pdf");
				hpull = frame->pullHist("data_slice2","pdf_slice2") ;
				frame= DTF_Bs_M.frame();
				frame->addPlotable(hpull,"P") ;
				frame->Draw();
				c->Print("eps/signal_pull_" + str_year[i] + "_" + str_Ds[j] + ".eps");
				if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/MassFit/signal_pull_" + str_year[i] + "_" + str_Ds[j] + ".pdf");
	
				/// Get signal yield
				signal_yield += ((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal();
				RooAbsPdf* pdf_slice_comb_bkg = (RooAbsPdf*) pdf_slice_comp.find("bkg_exp_{" + str_year[i] + ";" + str_Ds[j] + "}");
				comb_bkg_yield += pdf_slice_comb_bkg->createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("signal_range"))->getVal()*((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal();
	
				/// Calculate sWeights
				if(sWeight){
					RooArgList yield_list(
					*((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}")),
					*((RooRealVar*) fitParams->find("n_sig_B0_{"+str_year[i] + ";" + str_Ds[j] + "}")),
					*((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}")),
					*((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}")),
					*((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))
					);
					RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState_mod == Ds_finalState_mod::" + str_Ds[j]);
					SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
		
					/// Plot the sWeight distributions as a function of mass
					TH2 * swHist = (TH2*)data_slice->createHistogram("Bs_DTF_MM,n_sig_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
					swHist->GetYaxis()->SetTitle("Signal sWeights");
					swHist->Draw();
					c->Print("eps/signal_sweight_" + str_year[i] + "_" + str_Ds[j] + ".eps");
	
					/// Save sWeights
					/// Messy and dangerous hack but works for now
					int n_ij = 0;  /// labels entry number of data slice
					for(int n = 0; n < out_tree->GetEntries(); n++){
						b_year->GetEntry(n);
						b_Ds_finalState->GetEntry(n);
						if(A_is_in_B(anythingToString(t_year), (string) str_year[i]) && t_Ds_finalState == j){
							weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
							weights_B0[n] = sPlot.GetSWeight(n_ij,"n_sig_B0_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
							n_ij++;
						}
					}
				}
			}
		}
	}

	if(sWeight){
		for(int n = 0; n < out_tree->GetEntries(); n++){
			sw = weights[n];
			sw_B0 = weights_B0[n];
			b_sw->Fill();
			b_w->Fill();
		}

	 	out_tree->Write();
   		output->Close();
		cout << endl;
		cout << "Created file " << "/auto/data/dargent/BsDsKpipi/Final/Data/signal.root"  << endl << endl;
	}

	chi2 = chi2*nBins/(str_year.size()*str_Ds.size()*nBins-fitParams->selectByAttrib("Constant",kFALSE)->getSize());
	cout << endl; cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
	cout << "Total signal yield = " << signal_yield << endl;

	if(optimizeBDT){

		TString cut_string = (string) cut_BDT;
		double BDT = atoi(cut_string.ReplaceAll("BDTG > 0.",""))/10.;		

		TCanvas* c_1 = new TCanvas();
		TCanvas* c_2 = new TCanvas();
		c_1->SetGrid(1);
		c_2->SetGrid(1);

		TH1D* h_BDT_s = new TH1D("h_BDT_s",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_b = new TH1D("h_BDT_b",";BDTG; Yield [norm.]",40,-1,1);

		TH1D* h_BDT_Run1_t0_s = new TH1D("h_BDT_Run1_t0_s",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run1_t1_s = new TH1D("h_BDT_Run1_t1_s",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run2_t0_s = new TH1D("h_BDT_Run2_t0_s",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run2_t1_s = new TH1D("h_BDT_Run2_t1_s",";BDTG; Yield [norm.]",40,-1,1);

		TH1D* h_BDT_Run1_t0_b = new TH1D("h_BDT_Run1_t0_b",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run1_t1_b = new TH1D("h_BDT_Run1_t1_b",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run2_t0_b = new TH1D("h_BDT_Run2_t0_b",";BDTG; Yield [norm.]",40,-1,1);
		TH1D* h_BDT_Run2_t1_b = new TH1D("h_BDT_Run2_t1_b",";BDTG; Yield [norm.]",40,-1,1);

		if(useTriggerCat){

			double signal_yield_opt = 0;
			double bkg_yield_opt = 0;
			double signal_yield_opt_eff = 0;
			double bkg_yield_opt_eff = 0;

			float t_BDTG, t_weight;
			int t_id;

			for(int i=0; i<str_run.size(); i++)for(int k=0; k<str_trigger.size(); k++){

				cout << "***********" << endl << "Optimizing BDT cut for " << str_run[i] << " ; " << str_trigger[k] << endl << endl ;

				/// Get yields for categories
				double signal_yield_cat = 0;
				double bkg_yield_cat = 0;

				for(int j=0; j<str_Ds.size(); j++){
	
					RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}");
					RooArgList pdf_slice_comp( pdf_slice->pdfList());
			
					signal_yield_cat += ((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
					RooAbsPdf* pdf_slice_comb_bkg = (RooAbsPdf*) pdf_slice_comp.find("bkg_exp_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}");
					bkg_yield_cat += pdf_slice_comb_bkg->createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("signal_range"))->getVal()*((RooRealVar*) fitParams->find("n_exp_bkg_{" + str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
				}
	
				TString effFile = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/Selection/TMVA_Bs2DsKpipi_BDTG_Data_"+str_run[i]+"_"+str_trigger[k]+"_even.root";
				effFile.ReplaceAll("Run","run");				

				TFile* f_even= TFile::Open(effFile);
				TTree* tree_even=(TTree*)f_even->Get("TestTree");

				tree_even->SetBranchAddress("classID",&t_id);
				tree_even->SetBranchAddress("BDTG",&t_BDTG);
				tree_even->SetBranchAddress("weight",&t_weight);
				for(int n = 0; n < tree_even->GetEntries(); n++){
					tree_even->GetEntry(n);
					if(i==0 && k ==0){
						if(t_id==0)h_BDT_Run1_t0_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run1_t0_b->Fill(t_BDTG,t_weight);
					}
					if(i==0 && k ==1){
						if(t_id==0)h_BDT_Run1_t1_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run1_t1_b->Fill(t_BDTG,t_weight);
					}
					if(i==1 && k ==0){
						if(t_id==0)h_BDT_Run2_t0_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run2_t0_b->Fill(t_BDTG,t_weight);
					}
					if(i==1 && k ==1){
						if(t_id==0)h_BDT_Run2_t1_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run2_t1_b->Fill(t_BDTG,t_weight);
					}

					if(t_id==0)h_BDT_s->Fill(t_BDTG,t_weight);	
					else h_BDT_b->Fill(t_BDTG,t_weight);
				}

				TH1D* effS_even=(TH1D*)f_even->Get("Method_BDT/BDTG/MVA_BDTG_effS");
				TH1D* effB_even=(TH1D*)f_even->Get("Method_BDT/BDTG/MVA_BDTG_effB");

				effS_even->Rebin(100);
				effB_even->Rebin(100);
				effS_even->Scale(1./100.);
				effB_even->Scale(1./100.);

				effFile.ReplaceAll("even","odd");				
				TFile* f_odd= TFile::Open(effFile);
				TTree* tree_odd=(TTree*)f_odd->Get("TestTree");

				tree_odd->SetBranchAddress("classID",&t_id);
				tree_odd->SetBranchAddress("BDTG",&t_BDTG);
				tree_odd->SetBranchAddress("weight",&t_weight);
				for(int n = 0; n < tree_odd->GetEntries(); n++){
					tree_odd->GetEntry(n);
					if(i==0 && k ==0){
						if(t_id==0)h_BDT_Run1_t0_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run1_t0_b->Fill(t_BDTG,t_weight);
					}
					if(i==0 && k ==1){
						if(t_id==0)h_BDT_Run1_t1_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run1_t1_b->Fill(t_BDTG,t_weight);
					}
					if(i==1 && k ==0){
						if(t_id==0)h_BDT_Run2_t0_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run2_t0_b->Fill(t_BDTG,t_weight);
					}
					if(i==1 && k ==1){
						if(t_id==0)h_BDT_Run2_t1_s->Fill(t_BDTG,t_weight);	
						else h_BDT_Run2_t1_b->Fill(t_BDTG,t_weight);
					}
					if(t_id==0)h_BDT_s->Fill(t_BDTG,t_weight);	
					else h_BDT_b->Fill(t_BDTG,t_weight);
				}

				TH1D* effS_odd=(TH1D*)f_odd->Get("Method_BDT/BDTG/MVA_BDTG_effS");
				TH1D* effB_odd=(TH1D*)f_odd->Get("Method_BDT/BDTG/MVA_BDTG_effB");

				effS_odd->Rebin(100);
				effB_odd->Rebin(100);
				effS_odd->Scale(1./100.);
				effB_odd->Scale(1./100.);

				cout << "Yields with BDT > " << BDT << ":" << endl;
				cout << "N_s = " << signal_yield_cat << endl;
				cout << "N_b = " << bkg_yield_cat << endl;
				
				/// Correct yields for BDT efficency
				double eff_s = ( effS_even->GetBinContent(effS_even->FindBin(BDT)) + effS_odd->GetBinContent(effS_odd->FindBin(BDT)))/2.;
				double eff_b = ( effB_even->GetBinContent(effB_even->FindBin(BDT)) + effB_odd->GetBinContent(effB_odd->FindBin(BDT)))/2.;
				double n_s =  signal_yield_cat/eff_s ;
				double n_b =  bkg_yield_cat/eff_b ;
				
				cout << endl << "BDT eff.corrected yields :" << endl;
				cout << "N_s = " << n_s << endl;
				cout << "N_b = " << n_b << endl;
		
				/// Find optimal BDT cut
				TH1D* h_sig = (TH1D*)effS_even->Clone("h_sig");
				TH1D* h_eff_s = (TH1D*)effS_even->Clone("h_eff_s");
				TH1D* h_eff_b = (TH1D*)effS_even->Clone("h_eff_b");

				for(int n = 1; n <= h_sig->GetNbinsX(); n++){
					double e_s = (effS_even->GetBinContent(n)+effS_odd->GetBinContent(n))/2.;
					double e_b = (effB_even->GetBinContent(n)+effB_odd->GetBinContent(n))/2.;
					h_sig->SetBinContent(n, e_s * n_s / sqrt( e_s * n_s + e_b * n_b )  );
					h_eff_s->SetBinContent(n, e_s  );
					h_eff_b->SetBinContent(n, e_b );
				}	
		
				c_1->cd();
				h_sig->SetTitle(";BDTG; Signal Significance [norm.]");
				h_sig->Scale(1./h_sig->GetMaximum());
				h_sig->SetMaximum(1.05);
				h_sig->SetMinimum(0.8);
				h_sig->SetLineWidth(4);
				if(k==0)h_sig->SetLineColor(kBlue);
				if(k==1){h_sig->SetLineColor(kRed);h_sig->SetLineStyle(7);}
				//if(k==1)h_sig->SetMarkerStyle(4);
				if(k==0)h_sig->Draw("histc");
				else h_sig->Draw("histcsame");
				if(i==0 && k==1){
					c_1->Print("eps/BDT_scan_Run1.eps");
					c_1->Print("../../../../../TD-AnaNote/latex/figs/TMVA/BDT_scan_Run1.pdf");
				}
				if(i==1 && k==1){
					c_1->Print("eps/BDT_scan_Run2.eps");
					c_1->Print("../../../../../TD-AnaNote/latex/figs/TMVA/BDT_scan_Run2.pdf");
				}

				c_2->cd();
				h_eff_s->SetTitle(";BDTG; Efficiency");
				h_eff_s->SetLineWidth(4);
				if(k==0)h_eff_s->SetLineColor(kBlue);
				if(k==1){h_eff_s->SetLineColor(kRed);h_sig->SetLineStyle(7);}
				//if(k==1)h_sig->SetMarkerStyle(4);
				if(k==0)h_eff_s->Draw("histc");
				else h_eff_s->Draw("histcsame");
				if(i==0 && k==1)c_2->Print("eps/BDT_eff_Run1.eps");
				if(i==1 && k==1)c_2->Print("eps/BDT_eff_Run2.eps");
			
				/// Update values with optimal cut
				int max_bin = h_sig->GetMaximumBin();
				double cut_BDT_opt = h_sig->GetBinCenter(max_bin);	

				eff_s = (effS_even->GetBinContent(max_bin)+effS_odd->GetBinContent(max_bin))/2.;
				eff_b = (effB_even->GetBinContent(max_bin)+effB_even->GetBinContent(max_bin))/2.;
				n_s *=  eff_s ;
				n_b *=  eff_b ;
				
				cout << endl << "Optimal BDT cut : BDT > " << cut_BDT_opt << endl; 
				cout << "effS = " << eff_s << endl;
				cout << "effB = " << eff_b << endl;
				cout <<  "N_s = " << n_s << endl;
				cout <<  "N_b = " << n_b << endl;
				cout << "Sig = " << n_s/sqrt(n_s+n_b) << endl;	
		
				signal_yield_opt += n_s;
				bkg_yield_opt += n_b;

				/// Update values with eff = 0.9 cut
				n_s /=  eff_s ;
				n_b /=  eff_b ;
				max_bin = h_eff_s->FindLastBinAbove(0.9);
				cut_BDT_opt = h_eff_s->GetBinCenter(max_bin);	

				eff_s = (effS_even->GetBinContent(max_bin)+effS_odd->GetBinContent(max_bin))/2.;
				eff_b = (effB_even->GetBinContent(max_bin)+effB_even->GetBinContent(max_bin))/2.;
				n_s *=  eff_s ;
				n_b *=  eff_b ;
				
				cout << endl << "BDT cut with eff = 0.9 : BDT > " << cut_BDT_opt << endl; 
				cout << "effS = " << eff_s << endl;
				cout << "effB = " << eff_b << endl;
				cout <<  "N_s = " << n_s << endl;
				cout <<  "N_b = " << n_b << endl;
				cout << "Sig = " << n_s/sqrt(n_s+n_b) << endl;	
		
				signal_yield_opt_eff += n_s;
				bkg_yield_opt_eff += n_b;
			}

			cout << "***********" << endl << "Summary :" << endl;
			cout << "Optimal BDT cut " << endl;
			cout <<  "N_s = " << signal_yield_opt << endl;
			cout <<  "N_b = " << bkg_yield_opt << endl;
			cout << "Sig = " << signal_yield_opt/sqrt(signal_yield_opt+bkg_yield_opt) << endl;	

			cout << endl << "BDT cut with eff = 0.9" << endl;
			cout <<  "N_s = " << signal_yield_opt_eff << endl;
			cout <<  "N_b = " << bkg_yield_opt_eff << endl;
			cout << "Sig = " << signal_yield_opt_eff/sqrt(signal_yield_opt_eff+bkg_yield_opt_eff) << endl;	
				
			c->cd();
			gPad->SetLogy();

			h_BDT_s->SetLineColor(kBlue+1);
			h_BDT_s->SetMarkerColor(kBlue+1);
			h_BDT_s->SetFillColor(kBlue-4);

			h_BDT_b->SetLineColor(kRed+1);
			h_BDT_b->SetMarkerColor(kRed+1);
			h_BDT_b->SetFillColor(kRed+1);
			h_BDT_b->SetFillStyle(3335);

			h_BDT_Run1_t0_s->SetLineColor(kRed+1);
			h_BDT_Run1_t0_s->SetMarkerColor(kRed+1);
			h_BDT_Run1_t0_s->SetFillColor(kRed+1);
			h_BDT_Run1_t0_s->SetFillStyle(3335);

			h_BDT_Run1_t1_s->SetLineColor(kRed);
			h_BDT_Run1_t1_s->SetMarkerColor(kRed);

			h_BDT_Run1_t0_b->SetLineColor(kGray+1);
			h_BDT_Run1_t0_b->SetMarkerColor(kGray+1);
			h_BDT_Run1_t0_b->SetFillColor(kGray+1);

			h_BDT_Run2_t0_s->SetLineColor(kRed+1);
			h_BDT_Run2_t0_s->SetMarkerColor(kRed+1);
			h_BDT_Run2_t0_s->SetFillColor(kRed+1);
			h_BDT_Run2_t0_s->SetFillStyle(3335);

			h_BDT_Run2_t1_s->SetLineColor(kRed);
			h_BDT_Run2_t1_s->SetMarkerColor(kRed);

			h_BDT_Run2_t0_b->SetLineColor(kGray+1);
			h_BDT_Run2_t0_b->SetMarkerColor(kGray+1);
			h_BDT_Run2_t0_b->SetFillColor(kGray+1);

			h_BDT_b->SetMinimum(10);
			h_BDT_b->Draw("e1");
			h_BDT_s->Draw("e1same");
			c->Print("eps/BDTG.eps");
			c->Print("../../../../../TD-AnaNote/latex/figs/TMVA/BDTG.pdf");

			h_BDT_Run1_t0_b->SetMinimum(1);
			h_BDT_Run1_t0_b->DrawNormalized("hist",1);
			h_BDT_Run1_t1_b->DrawNormalized("e1same",1);
			h_BDT_Run1_t0_s->DrawNormalized("histsame",1);
			h_BDT_Run1_t1_s->DrawNormalized("e1same",1);
			c->Print("eps/BDTG_Run1.eps");
			c->Print("../../../../../TD-AnaNote/latex/figs/TMVA/BDTG_Run1.pdf");

			h_BDT_Run2_t0_b->SetMinimum(1);
			h_BDT_Run2_t0_b->DrawNormalized("hist",1);
			h_BDT_Run2_t1_b->DrawNormalized("e1same",1);
			h_BDT_Run2_t0_s->DrawNormalized("histsame",1);
			h_BDT_Run2_t1_s->DrawNormalized("e1same",1);
			c->Print("eps/BDTG_Run2.eps");
			c->Print("../../../../../TD-AnaNote/latex/figs/TMVA/BDTG_Run2.pdf");
		}
// 		else {
// 			TFile* f= TFile::Open("/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/Selection/TMVA_Bs2DsKpipi.root");
// 			TH1D* effS=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effS");
// 			TH1D* effB=(TH1D*)f->Get("Method_BDT/BDT/MVA_BDT_effB");
// 			
// 			cout << endl << "Optimizing BDT cut" << endl << endl ;
// 			
// 			cout << "Yields with BDT > " << (double)cut_BDT << ":" << endl;
// 			cout << "N_s = " << signal_yield << endl;
// 			cout << "N_b = " << comb_bkg_yield << endl;
// 			
// 			/// Correct yields for BDT efficency
// 			double eff_s = effS->GetBinContent(effS->FindBin((double)cut_BDT));
// 			double eff_b = effB->GetBinContent(effB->FindBin((double)cut_BDT));
// 			double n_s =  signal_yield/eff_s ;
// 			double n_b =  comb_bkg_yield/eff_b ;
// 			
// 			cout << endl << "BDT eff.corrected yields :" << endl;
// 			cout << "N_s = " << n_s << endl;
// 			cout << "N_b = " << n_b << endl;
// 	
// 			/// Find optimal BDT cut
// 			TH1D* h_sig = (TH1D*)effS->Clone("h_sig");
// 			for(int i = 1; i <= h_sig->GetNbinsX(); i++){
// 				double e_s = effS->GetBinContent(i);
// 				double e_b = effB->GetBinContent(i);
// 				h_sig->SetBinContent(i, e_s * n_s / sqrt( e_s * n_s + e_b * n_b )  );
// 			}	
// 	
// 			h_sig->Draw("histc");
// 			c->Print("eps/BDT_scan.eps");
// 	
// 			int max_bin = h_sig->GetMaximumBin();
// 			double cut_BDT_opt = h_sig->GetBinCenter(max_bin);		
// 			
// 			/// Update values with optimal cut
// 			eff_s = effS->GetBinContent(max_bin);
// 			eff_b = effB->GetBinContent(max_bin);
// 			n_s *=  eff_s ;
// 			n_b *=  eff_b ;
// 			
// 			cout << endl << "Optimal BDT cut : BDT > " << cut_BDT_opt << endl; 
// 			cout << "effS = " << eff_s << endl;
// 			cout << "effB = " << eff_b << endl;
// 			cout <<  "N_s = " << n_s << endl;
// 			cout <<  "N_b = " << n_b << endl;
// 			//cout << "S/B = " << n_s/n_b << endl;
// 			//cout << "S/(S+B) = " << n_s/(n_s+n_b) << endl;
// 			cout << "Sig = " << n_s/sqrt(n_s+n_b) << endl;			
// 		}

	}

	///create a new table for Ana Note
	if(updateAnaNotePlots && useTriggerCat){

		double yield_sig = 0;
		double yield_B0 = 0;
		double yield_partReco = 0;
		double yield_comb = 0;
		double yield_misID = 0;
		double yield_sig_err = 0;
		double yield_B0_err = 0;
		double yield_partReco_err = 0;
		double yield_comb_err = 0;
		double yield_misID_err = 0;

		for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			yield_sig += ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_sig_err += pow(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);

			yield_B0 += ((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_B0_err += pow(((RooRealVar*) fitParams->find("n_sig_B0_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);	

			yield_partReco += ((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_partReco_err += pow(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);

			yield_comb += ((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_comb_err += pow(((RooRealVar*) fitParams->find("n_exp_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);

			yield_misID += ((RooRealVar*) fitParams->find("n_misID_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_misID_err += pow(((RooRealVar*) fitParams->find("n_misID_bkg_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);
		}

		ofstream datafile;
		datafile.open("../../../../../TD-AnaNote/latex/tables/MassFit/yields_signal.tex",std::ofstream::trunc);
		//datafile << "\\begin{table}[h]" << "\n";
		//datafile << "\\centering" << "\n";
		datafile << " \\begin{tabular}{l r }" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "Component & Yield\\" << " \\\\" << "\n";
		datafile << "\\hline" << "\n";
		datafile << std::fixed << std::setprecision(0) << "$B_s \\to D_s K \\pi \\pi$" << " & "<< yield_sig << " $\\pm$ " << sqrt(yield_sig_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "$B^{0} \\to D_s K \\pi \\pi$" << " & "<< yield_B0 << " $\\pm$ " << sqrt(yield_B0_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "Partially reconstructed bkg." << " & "<< yield_partReco << " $\\pm$ " << sqrt(yield_partReco_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "Misidentified bkg." << " & "<< yield_misID << " $\\pm$ " << sqrt(yield_misID_err) << " \\\\" << "\n";
		datafile << std::fixed<< std::setprecision(0) << "Combinatorial bkg." << " & "<< yield_comb << " $\\pm$ " << sqrt(yield_comb_err) << " \\\\" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "\\end{tabular}" << "\n";
		//datafile << "\\caption{Total signal and background yields for the $B_s \\to D_s \\pi \\pi \\pi$ sample.}" << "\n";
		datafile << "\\label{table:signalYields}" << "\n";
		//datafile << "\\end{table}" << "\n";
		datafile.close();

		vector<double> yield_sig_Ds(4,0);
		vector<double> yield_sig_Ds_err(4,0);

		for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_Ds.size(); j++)for(int k=0; k<str_trigger.size(); k++){
			yield_sig_Ds[j] += ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getVal();
			yield_sig_Ds_err[j] += pow(((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_Ds[j]+ ";" + str_trigger[k] + "}"))->getError(),2);
		}

		ofstream datafileDs;
		datafile.open("../../../../../TD-AnaNote/latex/tables/MassFit/yields_signalDs.tex",std::ofstream::trunc);
		//datafile << "\\begin{table}[h]" << "\n";
		//datafile << "\\centering" << "\n";
		datafile << " \\begin{tabular}{l r }" << "\n";
		datafile << "\\hline\\hline" << "\n";
		datafile << "$D_s$ final state  & Signal yield\\" << " \\\\" << "\n";
		datafile << "\\hline" << "\n";
		for(int j=0; j<str_Ds.size(); j++){ 
			string caption;
			if(str_Ds[j]=="phipi")caption = "$D_{s}^{-} \\to \\phi^{0}(1020)\\pi^{-}$";
			else if(str_Ds[j]=="KsK")caption = "$D_{s}^{-}\\to K^{*0}(892)K^{-}$";
			else if(str_Ds[j]=="KKpi_NR")caption = "$D_{s}^{-}\\to (K^{-}h^{+}\\pi^{-})$";
			else if(str_Ds[j]=="pipipi")caption = "$D_{s}^{-}\\to \\pi^{+}\\pi^{-}\\pi^{-}$";

			datafile << std::fixed << std::setprecision(0) << caption << " & "<< yield_sig_Ds[j] << " $\\pm$ " << sqrt(yield_sig_Ds_err[j]) << " \\\\" << "\n";
		}
		datafile << "\\hline\\hline" << "\n";
		datafile << "\\end{tabular}" << "\n";
		//datafile << "\\caption{Signal yield for the different $D_s$ final states contributing to the $B_s \\to D_s \\pi \\pi \\pi$ sample.}" << "\n";
		datafile << "\\label{table:signalYieldsDs}" << "\n";
		//datafile << "\\end{table}" << "\n";
		datafile.close();
	}
}

void AnalyticPartBkg(){


	///define shape of Bs->Ds*pipipi BG using analytical shape RooHILLdini
	
	RooRealVar Bs_MM("Bs_DTF_MM", "m(D_{s}*#pi#pi#pi)", 5000., 5350.,"MeV/c^{2}");

	RooRealVar a("a","a", 5.3108e+03, 5040., 5970.);
	RooRealVar b("b","b", 5.1718e+03, 4140., 5905.);
	RooRealVar a_HORNS("a_HORNS","a_HORNS", 5.0591e+03, 5000.,5200.);
	RooRealVar b_HORNS("b_HORNS","b_HORNS", 5.1741e+03, 5100.,5300.);
	RooRealVar csi("csi","csi", 1.,0.,10.);
	RooRealVar csi_HORNS("csi_HORNS","csi_HORNS", 1.,0.,10.);
	RooRealVar shift("shift","shift", 0., 0., 200.);

	RooRealVar sigma("sigma", "sigma", 20.,0.,100.);

	///fixed from othere LHCb analysis 
	RooRealVar ratio_sigma("ratio_sigma", "ratio_sigma", 1);
	RooRealVar fraction_sigma("fraction_sigma", "fraction_sigma", 1);

	RooRealVar f_1("f_1", "f_1", 0.935);//, 0., 1.);


	RooHILLdini RooHILLBkgShape("RooHILLBkgShape", "RooHILLBkgShape", Bs_MM, a, b, csi, shift, sigma, ratio_sigma, fraction_sigma );
	RooHORNSdini RooHORNSBkgShape("RooHORNSBkgShape", "RooHORNSBkgShape", Bs_MM, a_HORNS, b_HORNS, csi_HORNS, shift, sigma, ratio_sigma, fraction_sigma );
	RooAddPdf* bkg_partReco_alt = new RooAddPdf("bkg_partReco_alt", "bkg_partReco_alt", RooArgList(RooHILLBkgShape, RooHORNSBkgShape), RooArgList(f_1),kTRUE);

	/// Load file
	TFile* file = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12_Dstar_bkg.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs*MM",1);
	tree->SetBranchStatus("m_*",1);

	//Fill needed variable in RooDataSet
	RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));
	
	/// Fit
	RooFitResult *result = bkg_partReco_alt->fitTo(*data,Save(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 
	//plot mass distribution and fit results
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(50));
	bkg_partReco_alt->plotOn(frame_m,Name("RooHILLBkgShape"),Components(RooHILLBkgShape),LineStyle(kDashed),LineColor(kGreen),LineWidth(1));
	bkg_partReco_alt->plotOn(frame_m,Name("RooHORNSBkgShape"),Components(RooHORNSBkgShape),LineStyle(kDashed),LineColor(kRed),LineWidth(1));
	bkg_partReco_alt->plotOn(frame_m,Name("bkg_partReco_alt"),LineColor(kBlue),LineWidth(2));
	frame_m->Draw();
	c1->Print("eps/BkgShape/RooHILLBkgShape_Bs2Dsstartpipipi.eps");

	file->Close();
	

}


int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");

    NamedParameter<string> Channel("channel", (std::string) "Norm" , (char*) 0);
    string channel = Channel;

    str_year.push_back(TString("y11"));
    str_year.push_back(TString("y12"));
    str_year.push_back(TString("y15"));
    str_year.push_back(TString("y16"));
    str_year.push_back(TString("y17"));
    
    str_run.push_back(TString("Run1"));
    str_run.push_back(TString("Run2"));

    str_Ds.push_back(TString("phipi"));
    str_Ds.push_back(TString("KsK"));
    str_Ds.push_back(TString("KKpi_NR"));
    str_Ds.push_back(TString("pipipi"));
    //str_Ds.push_back(TString("Kpipi"));

    str_trigger.push_back(TString("t0"));
    str_trigger.push_back(TString("t1"));
    
    if(channel == "Norm") fitNorm();
    else if(channel == "Signal") fitSignal();
    else if(channel == "PartBkg") AnalyticPartBkg();
    else{
	cout << "*********************************************************" << endl;
    	cout << "please specify 'Signal' or 'Norm' in options file" << endl;
	cout << "*********************************************************" << endl;
    }
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
