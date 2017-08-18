// Fits Bs mass distribution and calculates sweights
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
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
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
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

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

vector<TString> str_year;
vector<TString> str_Ds;
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
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(50));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(BifGauss1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(BifGauss2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(BifGauss3),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstartpipipi.eps");

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

	file->Close();
	
	return params;
}

vector<double> fitMisIdBkgShape_Ds3pi(){

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
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40),DataError(RooAbsData::SumW2));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi.eps");

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
	data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(40));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
	pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	frame_m->Draw();
	c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.eps");

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
	NamedParameter<double> cut_BDT("cut_BDT",0.);

	/// Load file
	TString inFileName = "/auto/data/dargent/BsDsKpipi/BDT/MC/"+channel+".root";
	TFile* file = new TFile(inFileName);	
	TTree* tree = (TTree*) file->Get("DecayTree");	

	TFile* output = new TFile(inFileName.ReplaceAll("/BDT/","/Final/"),"RECREATE");
	TTree* out_tree = tree->CopyTree(("Bs_DTF_MM >= 5320. && Bs_DTF_MM <= 5420 && Bs_BKGCAT <= 20 && BDTG_response > " + anythingToString((double)cut_BDT) ).c_str() );

        RooRealVar DTF_Bs_M("Bs_DTF_MM", "m(D_{s} K #pi #pi)", 5320., 5420.,"MeV");
        RooRealVar BDTG_response("BDTG_response", "BDTG_response", 0.);
        RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
        RooRealVar Bs_BKGCAT("Bs_BKGCAT", "Bs_BKGCAT", 0.);
	RooRealVar EventWeight("weight","weight", 0.);

	RooArgList list =  RooArgList(DTF_Bs_M,BDTG_response,Ds_finalState,Bs_BKGCAT,EventWeight);
        RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree),WeightVar(EventWeight));
	
	/// Signal pdf
	RooRealVar mean("mean", "#mu", 5366.89,5350.,5390.); 
	RooRealVar sigma("sigma", "#sigma", 20.,0.,80.); 
	RooRealVar gamma("gamma", "#gamma", -0.5,-5,5.); 
	RooRealVar delta("delta", "#delta", 0.5,-5,5.); 
	RooJohnsonSU* signal= new RooJohnsonSU("signal","signal",DTF_Bs_M, mean,sigma,gamma,delta);

	/// Fit
	RooFitResult* result = signal->fitTo(*data,Save(kTRUE),NumCPU(3),SumW2Error(kTRUE));
	result->Print();

	/// Plotting
	TCanvas* c = new TCanvas();
	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
	data->plotOn(frame,Name("data"),Binning(50));
	signal->plotOn(frame,Name("signal"));
	frame->Draw();
	c->Print("eps/SignalShape/"+channel+"MC.eps");

	double chi2 = 0;
	chi2 = frame->chiSquare("signal","data",4);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

	RooHist* hpull = frame->pullHist("data","signal") ;
	frame= DTF_Bs_M.frame();
	frame->addPlotable(hpull,"P") ;
	frame->Draw();
        c->Print("eps/SignalShape/"+channel+"MC_pull.eps");

	/// Return fit params
	vector<double> params;
	params.push_back(mean.getVal());
	params.push_back(sigma.getVal());
	params.push_back(gamma.getVal());
	params.push_back(delta.getVal());

	out_tree->Write();
	file->Close();
	output->Close();
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
	NamedParameter<double> cut_BDT("cut_BDT",0.);
	NamedParameter<string> inFileName("inFileNameNorm",(string)"/auto/data/dargent/BsDsKpipi/BDT/Data/norm.root");
	NamedParameter<string> outFileName("outFileNameNorm",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");

	/// Define categories
	RooCategory year("year","year") ;
	year.defineType("y11",11);
	year.defineType("y12",12);
	year.defineType("y15",15);
	year.defineType("y16",16);
	RooCategory Ds_finalState("Ds_finalState","Ds_finalState") ;
	Ds_finalState.defineType("phipi",0);
	Ds_finalState.defineType("KsK",1);
	Ds_finalState.defineType("KKpi_NR",2);
	Ds_finalState.defineType("pipipi",3);

	/// Load file
	TFile *file= new TFile(((string)inFileName).c_str());
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	tree->SetBranchStatus("BDTG_response",1);
	tree->SetBranchStatus("Ds_finalState",1);
	tree->SetBranchStatus("year",1);

        RooRealVar DTF_Bs_M("Bs_DTF_MM", "m(D_{s}^{-} #pi^{+}#pi^{+}#pi^{-})", min_MM, max_MM,"MeV/c^{2}");
        RooRealVar BDTG_response("BDTG_response", "BDTG_response", 0.);              
	RooArgList list =  RooArgList(DTF_Bs_M,BDTG_response,Ds_finalState,year);
	
	RooDataSet*  data = new RooDataSet("data","data",tree,list,("BDTG_response > " + anythingToString((double)cut_BDT)).c_str() );	

	/// Signal Pdf
	vector<double> sig_params = fitSignalShape("norm");
	RooRealVar mean_MC("mean_MC", "#mu MC", sig_params[0]); 
	RooRealVar sigma_MC("sigma_MC", "#sigma MC", sig_params[1]);
	RooRealVar scale_mean("scale_mean", "scale #mu",1.,0.9,1.); 
	RooRealVar scale_sigma("scale_sigma", "scale #sigma", 1.2, 0.9,1.5);
	RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
	RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
	RooRealVar alpha("alpha", "#alpha", sig_params[2]); 
	RooRealVar beta("beta", "#beta", sig_params[3]); 

	RooJohnsonSU signal("signal","signal",DTF_Bs_M, mean,sigma,alpha,beta);

	/// Combinatorial bkg pdf
	RooRealVar exp_par("exp_par","#lambda",-1.6508e-03,-10.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential bkg",DTF_Bs_M,exp_par);
	bkg_exp.fitTo(*data,Save(kTRUE),Range(5600.,5800.));

	/// Part. reco bkg
	vector<double> bkg_partReco_params = fitPartRecoBkgShape();
	RooRealVar mean1("mean1","mu", bkg_partReco_params[0]);
	RooRealVar mean2("mean2","mu", bkg_partReco_params[1]);
	RooRealVar mean3("mean3","mu", bkg_partReco_params[2]);
	RooRealVar sigmaL1("sigmaL1", "sigmaL1",  bkg_partReco_params[3]);
	RooRealVar sigmaR1("sigmaR1", "sigmaR1",  bkg_partReco_params[4]);
	RooRealVar sigmaL2("sigmaL2", "sigmaL2",  bkg_partReco_params[5],bkg_partReco_params[5]*0.,bkg_partReco_params[5]*5);
	RooRealVar sigmaR2("sigmaR2", "sigmaR2",  bkg_partReco_params[6],bkg_partReco_params[6]*0.,bkg_partReco_params[6]*1.4);
	RooRealVar sigmaL3("sigmaL3", "sigmaL3",  bkg_partReco_params[7]);//,bkg_partReco_params[7]*0.5,bkg_partReco_params[7]*1.5);
	RooRealVar sigmaR3("sigmaR3", "sigmaR3",  bkg_partReco_params[8],bkg_partReco_params[8]*0.5,bkg_partReco_params[8]*1.5);
	RooRealVar f_1("f_1", "f_1", bkg_partReco_params[9],0,1);
	RooRealVar f_2("f_2", "f_2", bkg_partReco_params[10],0,1);

	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", DTF_Bs_M, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", DTF_Bs_M, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", DTF_Bs_M, mean3, sigmaL3,sigmaR3);
	RooAddPdf bkg_partReco("bkg_partReco", "bkg_partReco", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2),kTRUE);

	if(ignorePartRecoBkg){
		RooArgSet* fitParamsPartRecoBkg = bkg_partReco.getParameters(data);
		RooFIter iterat = fitParamsPartRecoBkg->fwdIterator();	
		RooAbsArg * next = 0;
		while(next=iterat.next()) ((RooRealVar*)next)->setConstant();
	}

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/4., 0., data->numEntries());
	RooRealVar n_exp_bkg("n_exp_bkg", "n_exp_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()/2., 0., data->numEntries() );
	if(ignorePartRecoBkg){
		n_partReco_bkg.setVal(0.);
		n_partReco_bkg.setConstant();
	}

	RooAddPdf pdf("pdf", "pdf", RooArgList(signal, bkg_exp, bkg_partReco), RooArgList(n_sig, n_exp_bkg, n_partReco_bkg));

	/// Generate simultaneous pdf out of prototype pdf 
  	RooSimPdfBuilder* mgr = new RooSimPdfBuilder(pdf) ;
  	RooArgSet* config = mgr->createProtoBuildConfig() ;
  	config->setStringValue("physModels","pdf") ;
  	config->setStringValue("splitCats" ,"year Ds_finalState") ;
  	config->setStringValue("pdf", "year            : scale_sigma "
				      "Ds_finalState :  exp_par "  
                                      "year,Ds_finalState : n_sig, n_exp_bkg, n_partReco_bkg") ;  
	
  	RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
	simPdf->Print("v") ;
	RooArgSet* fitParams = simPdf->getParameters(data);

	/// Fix Exp from Sidebands
	for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
		/// Get data slice
		RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]);
		/// Fit
		bkg_exp.fitTo(*data_slice,Save(kTRUE),Range(5500.,5700.));
		/// Fix parameters
		//((RooRealVar*) fitParams->find("exp_par_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(exp_par.getVal());
		//if(fixExpBkgFromSidebands)((RooRealVar*) fitParams->find("exp_par_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
		((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setVal(exp_par.getVal());
		if(fixExpBkgFromSidebands)((RooRealVar*) fitParams->find("exp_par_"+ str_Ds[j]))->setConstant();
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
	TString last_name_signal,last_name_exp_bkg,last_name_partReco_bkg;
	for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
		RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_year[i] + ";" + str_Ds[j] + "}");
		RooArgList pdf_slice_comp( pdf_slice->pdfList());

		TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
		TString name_exp_bkg("exp_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
		TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));

		if(i==0 && j == 0){
			pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			
			pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
			
			pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
		}
		else{
			pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
			
			pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
	
			pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
		}

		if(i== str_year.size()-1 && j == str_Ds.size() -1) {

			pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));			
			pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[2])),ProjWData(year,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));

			pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
			pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(year,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));

			pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(year,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
			
			leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s} #rightarrow D_{s}^{-} #pi^{+}#pi^{+}#pi^{-}}","f");
			leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
			if(!ignorePartRecoBkg)leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
		}
		else {
			last_name_signal = name_signal;
			last_name_exp_bkg = name_exp_bkg;
			last_name_partReco_bkg = name_partReco_bkg;
		}
	}

	simPdf->plotOn(frame,Name("pdf"),ProjWData(year,*data),LineColor(kBlue+1),LineWidth(3));
	frame->Draw();
	leg.Draw();
	c->Print("eps/norm.eps");

	double chi2 = 0.;
	double covmatr = result->covQual();
	double edm = result->edm();
	RooHist* hpull = frame->pullHist("data","pdf") ;
	frame= DTF_Bs_M.frame();
	frame->addPlotable(hpull,"P") ;
        frame->Draw();
        c->Print("eps/norm_pull.eps");

	/// Output file
	TFile *output;
	TTree* out_tree = 0;
	double sw;
	Int_t t_year, t_Ds_finalState;
	TBranch *b_sw, *b_year,*b_Ds_finalState;

	if(sWeight){
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);

		out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) + " && BDTG_response > " + anythingToString((double)cut_BDT) ).c_str() );
    		b_sw = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");
        	out_tree->SetBranchAddress("year", &t_year, &b_year);
        	out_tree->SetBranchAddress("Ds_finalState", &t_Ds_finalState, &b_Ds_finalState);
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
	vector< double> partReco_yields;
	vector< double> scaleFactors;

	/// Loop over pdf slices
	for(int i=0; i<str_year.size(); i++){

		TLatex* lhcbtext = new TLatex();
		lhcbtext->SetTextFont(132);
		lhcbtext->SetTextColor(1);
		lhcbtext->SetTextSize(0.07);
		lhcbtext->SetTextAlign(13);
		lhcbtext->SetNDC(1);

		scaleFactors.push_back(((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->getVal());

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
			data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
			pdf_slice->plotOn(frame,LineColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
			pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components("bkg_partReco"),Normalization(1.,RooAbsReal::RelativeExpected));
			pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
			pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
			pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_year[i] + ";" + str_Ds[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
			pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
			//data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
			frame->Draw();
			chi2 += frame->chiSquare("pdf_slice2","data_slice2");

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
			c->Print("eps/norm_" + str_year[i] + "_" + str_Ds[j] + ".eps");
			hpull = frame->pullHist("data_slice2","pdf_slice2") ;
			frame= DTF_Bs_M.frame();
			frame->addPlotable(hpull,"P") ;
        		frame->Draw();
        		c->Print("eps/norm_pull_" + str_year[i] + "_" + str_Ds[j] + ".eps");

			/// Get signal yield
			yield += ((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal();
			signal_yields.push_back(((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal());
			partReco_yields.push_back(((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->getVal());

			/// Calculate sWeights
			if(sWeight){
				RooArgList yield_list(
				*((RooRealVar*) fitParams->find("n_sig_{"+str_year[i] + ";" + str_Ds[j] + "}")),
				*((RooRealVar*) fitParams->find("n_exp_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))
				);
				if(!ignorePartRecoBkg)yield_list.add(*((RooRealVar*) fitParams->find("n_partReco_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}")));
				RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]);
				pdf_slice->Print();
				SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
	
				/// Plot the sWeight distributions as a function of mass
				TH2 * swHist = (TH2*)data_slice->createHistogram("Bs_DTF_MM,n_sig_{" +str_year[i] + ";" + str_Ds[j] + "}" + "_sw");
				swHist->GetYaxis()->SetTitle("Signal sWeights");
				swHist->Draw();
				c->Print("eps/norm_sweight_" + str_year[i] + "_" + str_Ds[j] + ".eps");

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

	if(sWeight){
		for(int n = 0; n < out_tree->GetEntries(); n++){
			sw = weights[n];
			b_sw->Fill();
		}
	 	out_tree->Write();
   		output->Close();
		cout << endl;
		cout << "Created file " << "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root"  << endl << endl;
	}
	
	chi2 = chi2*nBins/(str_year.size()*str_Ds.size()*nBins-fitParams->selectByAttrib("Constant",kFALSE)->getSize());
	cout << endl <<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
	cout << "Total signal yield = " << yield << endl;

	/// Return fit params
	scaleFactors.push_back(scale_mean.getVal());

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

	vector< vector <double> > return_vec;
	return_vec.push_back(signal_yields);
	return_vec.push_back(scaleFactors);
	return_vec.push_back(partReco_yields);
	return_vec.push_back(partReco_params);

	return return_vec;
}

void fitSignal(){

	///Options
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> sWeight("sWeightSignal", 0);
	NamedParameter<int> nBins("nBins", 80);
	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);
	NamedParameter<double> cut_BDT("cut_BDT",0.);
	NamedParameter<int> fixMisIDyields("fixMisIDyields",1);
	NamedParameter<int> useNormScaleFactors("useNormScaleFactors",1);
	NamedParameter<int> optimizeBDT("optimizeBDT",0);

	NamedParameter<string> inFileName("inFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");
	NamedParameter<string> outFileName("outFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");

	///define categories
	RooCategory year("year","year") ;
	year.defineType("y11",11);
	year.defineType("y12",12);
	year.defineType("y15",15);
	year.defineType("y16",16);
	RooCategory Ds_finalState("Ds_finalState","Ds_finalState") ;
	Ds_finalState.defineType("phipi",0);
	Ds_finalState.defineType("KsK",1);
	Ds_finalState.defineType("KKpi_NR",2);
	Ds_finalState.defineType("pipipi",3);

	///Load file
	TFile *file= new TFile(((string)inFileName).c_str());
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	tree->SetBranchStatus("BDTG_response",1);
	tree->SetBranchStatus("Ds_finalState",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*_PIDK",1);

        RooRealVar DTF_Bs_M("Bs_DTF_MM", "m(D_{s}^{-} K^{+}#pi^{+}#pi^{-})", min_MM, max_MM,"MeV/c^{2}");
        RooRealVar BDTG_response("BDTG_response", "BDTG_response", 0.);        
        //RooRealVar K_plus_PIDK("K_plus_PIDK", "K_plus_PIDK", 0.);             
	RooArgList list =  RooArgList(DTF_Bs_M,BDTG_response,Ds_finalState,year);
	
	RooDataSet*  data = new RooDataSet("data","data",tree,list,("BDTG_response > " + anythingToString((double)cut_BDT)).c_str() );	

	/// Fit normalization mode first
	vector< vector <double> > norm_paramSet = fitNorm();

	/// Signal Pdf
	vector<double> sig_params = fitSignalShape("signal");
	RooRealVar mean_MC("mean_MC", "#mu MC", sig_params[0]); 
	RooRealVar sigma_MC("sigma_MC", "#sigma MC", sig_params[1]);
	RooRealVar scale_mean("scale_mean", "scale #mu",1.,0.,2. ); 
	RooRealVar scale_sigma("scale_sigma", "scale #sigma", 1.,0.,2.);
	RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
	RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
	RooRealVar alpha("alpha", "#alpha", sig_params[2]); 
	RooRealVar beta("beta", "#beta", sig_params[3]); 

	RooJohnsonSU signal("signal","signal",DTF_Bs_M, mean,sigma,alpha,beta);

	/// B0 pdf
	RooFormulaVar mean_B0("mean_B0","@0 - @1", RooArgSet(mean,RooConst(87.33))); 
	RooRealVar scale_sigma_B0("scale_sigma_B0", "scale #sigma_B0", 1.,0.,2.);
	RooFormulaVar sigma_B0("sigma_B0","@0 * @1", RooArgSet(scale_sigma_B0,sigma_MC)); 

	RooJohnsonSU signal_B0("signal_B0","signal_B0",DTF_Bs_M, mean_B0,sigma_B0,alpha,beta);

	/// Combinatorial bkg pdf
	RooRealVar exp_par("exp_par","#lambda",-1.6508e-03,-10.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential bkg",DTF_Bs_M,exp_par);
	bkg_exp.fitTo(*data,Save(kTRUE),Range(5600.,5800.));

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

	RooBifurGauss BifGauss1_Bs("BifGauss1_Bs","BifGauss1_Bs", DTF_Bs_M, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2_Bs("BifGauss2_Bs","BifGauss2_Bs", DTF_Bs_M, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3_Bs("BifGauss3_Bs","BifGauss3_Bs", DTF_Bs_M, mean3, sigmaL3,sigmaR3);
	RooAddPdf bkg_partReco_Bs("bkg_partReco_Bs", "bkg_partReco_Bs", RooArgList(BifGauss1_Bs, BifGauss2_Bs, BifGauss3_Bs), RooArgList(f_1,f_2));

	RooBifurGauss BifGauss1_B0("BifGauss1_B0","BifGauss1_B0", DTF_Bs_M, mean1Shifted, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2_B0("BifGauss2_B0","BifGauss2_B0", DTF_Bs_M, mean2Shifted, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3_B0("BifGauss3_B0","BifGauss3_B0", DTF_Bs_M, mean3Shifted, sigmaL3,sigmaR3);
	RooAddPdf bkg_partReco_B0("bkg_partReco_B0", "bkg_partReco_B0", RooArgList(BifGauss1_B0, BifGauss2_B0, BifGauss3_B0), RooArgList(f_1,f_2));

	RooRealVar partReco_f("partReco_f", "partReco_f", .5,0,1.);
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
	RooRealVar n_partReco_Bs_bkg("n_partReco_Bs_bkg", "n_partReco_Bs_bkg", data->numEntries()/2., 0., data->numEntries());
	//RooRealVar n_partReco_B0_bkg("n_partReco_B0_bkg", "n_partReco_B0_bkg", data->numEntries()/2., 0., data->numEntries());
	//RooRealVar n_partReco_Bs_bkg("n_partReco_Bs_bkg", "n_partReco_Bs_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_misID_bkg("n_misID_bkg", "n_misID_bkg", data->numEntries()*0.01, 0., data->numEntries()*0.25);

	RooAddPdf pdf("pdf", "pdf", RooArgList(signal, signal_B0, bkg_exp, bkg_partReco, bkg_misID), RooArgList(n_sig, n_sig_B0, n_exp_bkg, n_partReco_bkg, n_misID_bkg));

	/// Generate simultaneous pdf out of prototype pdf 
  	RooSimPdfBuilder* mgr = new RooSimPdfBuilder(pdf) ;
  	RooArgSet* config = mgr->createProtoBuildConfig() ;
  	config->setStringValue("physModels","pdf") ;
  	config->setStringValue("splitCats" ,"year Ds_finalState") ;
  	config->setStringValue("pdf", "year            : scale_sigma "
				      "Ds_finalState : exp_par "  
                                      "year,Ds_finalState : n_sig, n_sig_B0, n_exp_bkg, n_partReco_bkg, n_misID_bkg, misID_f") ;  
	
  	RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
	simPdf->Print("v") ;
	RooArgSet* fitParams = simPdf->getParameters(data);
  	fitParams->Print("v") ;

	/// Fix scale factors from normalization mode
	vector<double> scaleFactors = norm_paramSet[1];
	for(int i=0; i<str_year.size(); i++){
		((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->setVal( scaleFactors[i] );
		if(useNormScaleFactors)((RooRealVar*) fitParams->find("scale_sigma_"+str_year[i]))->setConstant();
	}
	scale_mean.setVal(scaleFactors[scaleFactors.size()-1]);
	if(useNormScaleFactors)scale_mean.setConstant();

	/// Fix misID yields 
	vector<double> Ds_yields = norm_paramSet[0];
	vector<double> Dstar_yields = norm_paramSet[2];
	
	int counter = 0;
	for(int i=0; i<str_year.size(); i++) for(int j=0; j<str_Ds.size(); j++){
		double val = fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]   ;
		((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal(val);
		if(fixMisIDyields)((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
		else ((RooRealVar*) fitParams->find("n_misID_bkg_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setRange(val * 0. , val * 2.);

		((RooRealVar*) fitParams->find("misID_f_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setVal( fake_prob_Ds * Ds_yields[counter] / (fake_prob_Ds * Ds_yields[counter] + fake_prob_Dstar * Dstar_yields[counter]   ));
		((RooRealVar*) fitParams->find("misID_f_{"+str_year[i] + ";" + str_Ds[j] + "}"))->setConstant();
		counter++;
	}

	/// Perform fit
	RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
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

			leg.AddEntry(frame->findObject(name_signal),"#font[132]{B_{s} #rightarrow D_{s}^{-} K^{+}#pi^{+}#pi^{-}}","f");
			leg.AddEntry(frame->findObject(name_signal_B0),"#font[132]{B^{0} #rightarrow D_{s}^{-} K^{+}#pi^{+}#pi^{-}}","f");
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

	frame->Draw();
	leg.Draw();
	c->Print("eps/signal.eps");

	double chi2 = 0.;
	double covmatr = result->covQual();
	double edm = result->edm();
	RooHist* hpull = frame->pullHist("data","pdf") ;
	frame= DTF_Bs_M.frame();
	frame->addPlotable(hpull,"P") ;
        frame->Draw();
        c->Print("eps/signal_pull.eps");

	/// Output file
	TFile *output;
	TTree* out_tree;
	double sw;
	Int_t t_year, t_Ds_finalState;
	TBranch *b_sw, *b_year,*b_Ds_finalState;

	if(sWeight){
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);

		out_tree = tree->CopyTree(("Bs_DTF_MM >= " + anythingToString((double)min_MM) + " && Bs_DTF_MM <= " + anythingToString((double)max_MM) + " && BDTG_response > " + anythingToString((double)cut_BDT) ).c_str() );
    		b_sw = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");
        	out_tree->SetBranchAddress("year", &t_year, &b_year);
        	out_tree->SetBranchAddress("Ds_finalState", &t_Ds_finalState, &b_Ds_finalState);
		if(out_tree->GetEntries() != data->numEntries()) {
			cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
			cout << out_tree->GetEntries() << endl;
			cout << data->numEntries() << endl;
			throw "ERROR";
		} 
	}
	double weights[(int)data->numEntries()];

	/// Calculate total signal yield
	double signal_yield = 0.;
	double comb_bkg_yield = 0.;
	DTF_Bs_M.setRange("signal_range",mean.getVal()-45.,mean.getVal()+45.);

	/// Loop over pdf slices
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
			//pdf_slice->Print();
			/// Plot data and pdf slices
			frame=DTF_Bs_M.frame();
			data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
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
			hpull = frame->pullHist("data_slice2","pdf_slice2") ;
			frame= DTF_Bs_M.frame();
			frame->addPlotable(hpull,"P") ;
        		frame->Draw();
        		c->Print("eps/signal_pull_" + str_year[i] + "_" + str_Ds[j] + ".eps");

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
				RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]);
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
						n_ij++;
					}
				}
			}
		}
	}

	if(sWeight){
		for(int n = 0; n < out_tree->GetEntries(); n++){
			sw = weights[n];
			b_sw->Fill();
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
		TFile* f= TFile::Open("/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/Selection/TMVA_Bs2DsKpipi.root");
		TH1D* effS=(TH1D*)f->Get("Method_BDT/BDTG/MVA_BDTG_effS");
		TH1D* effB=(TH1D*)f->Get("Method_BDT/BDTG/MVA_BDTG_effB");
		
		cout << endl << "Optimizing BDT cut" << endl << endl ;
		
		cout << "Yields with BDT > " << (double)cut_BDT << ":" << endl;
		cout << "N_s = " << signal_yield << endl;
		cout << "N_b = " << comb_bkg_yield << endl;
		
		/// Correct yields for BDT efficency
		double eff_s = effS->GetBinContent(effS->FindBin((double)cut_BDT));
		double eff_b = effB->GetBinContent(effB->FindBin((double)cut_BDT));
		double n_s =  signal_yield/eff_s ;
		double n_b =  comb_bkg_yield/eff_b ;
		
		cout << endl << "BDT eff.corrected yields :" << endl;
		cout << "N_s = " << n_s << endl;
		cout << "N_b = " << n_b << endl;

		/// Find optimal BDT cut
		TH1D* h_sig = (TH1D*)effS->Clone("h_sig");
		for(int i = 1; i <= h_sig->GetNbinsX(); i++){
			double e_s = effS->GetBinContent(i);
			double e_b = effB->GetBinContent(i);
			h_sig->SetBinContent(i, e_s * n_s / sqrt( e_s * n_s + e_b * n_b )  );
		}	

		h_sig->Draw("histc");
		c->Print("eps/BDT_scan.eps");

		int max_bin = h_sig->GetMaximumBin();
		double cut_BDT_opt = h_sig->GetBinCenter(max_bin);		
		
		/// Update values with optimal cut
		eff_s = effS->GetBinContent(max_bin);
		eff_b = effB->GetBinContent(max_bin);
		n_s *=  eff_s ;
		n_b *=  eff_b ;
		
		cout << endl << "Optimal BDT cut : BDT > " << cut_BDT_opt << endl; 
		cout << "effS = " << eff_s << endl;
		cout << "effB = " << eff_b << endl;
		cout <<  "N_s = " << n_s << endl;
		cout <<  "N_b = " << n_b << endl;
		//cout << "S/B = " << n_s/n_b << endl;
		//cout << "S/(S+B) = " << n_s/(n_s+n_b) << endl;
		cout << "Sig = " << n_s/sqrt(n_s+n_b) << endl;			
	}
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
    
    str_Ds.push_back(TString("phipi"));
    str_Ds.push_back(TString("KsK"));
    str_Ds.push_back(TString("KKpi_NR"));
    str_Ds.push_back(TString("pipipi"));
    
    if(channel == "Norm") fitNorm();
    else if(channel == "Signal") fitSignal();
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
