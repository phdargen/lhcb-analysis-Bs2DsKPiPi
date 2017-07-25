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
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"


using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

//define global arrays to be filled with params from auxiliary fits
double fitValues_BGShapeNorm[11];
double fitValues_BGShapethreePi[9];
double fitValues_BGShapethreePiDstar[9];
double fitvalues_BDTNorm[15];

double *fitBGShapeNorm(){

//define shape of Bs->Ds*pipipi BG as 3 bifurcated gaussians
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 5000., 5350.,"MeV/c^{2}");

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
RooRealVar f_1("f_{1}", "fraction1", 0.405, 0., 0.6);
RooRealVar f_2("f_{2}", "fraction2", 0.1, 0., 0.6);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2));

//Load file
TFile* file;
file= new TFile("/auto/data/dargent/old_Bs2DsKpipi/MC/Norm/Bkg/Dsstpipipi.root");
TTree* tree = (TTree*) file->Get("DecayTree");
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_MM",1);

//Fill needed variable in RooDataSet
RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

//Fit
RooFitResult *result;
result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3));
cout << "result is --------------- "<<endl;
result->Print(); 

//fill array to return the fit values
fitValues_BGShapeNorm[0] = mean1.getVal();
fitValues_BGShapeNorm[1] = mean2.getVal();
fitValues_BGShapeNorm[2] = mean3.getVal();
fitValues_BGShapeNorm[3] = sigmaL1.getVal();
fitValues_BGShapeNorm[4] = sigmaR1.getVal();
fitValues_BGShapeNorm[5] = sigmaL2.getVal();
fitValues_BGShapeNorm[6] = sigmaR2.getVal();
fitValues_BGShapeNorm[7] = sigmaL3.getVal();
fitValues_BGShapeNorm[8] = sigmaR3.getVal();
fitValues_BGShapeNorm[9] = f_1.getVal();
fitValues_BGShapeNorm[10] = f_2.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_MM.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(BifGauss1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(BifGauss2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(BifGauss3),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dsstartpipipi.eps");

return fitValues_BGShapeNorm;
}

double *fitBGShapethreePi(){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_Mass("Bs_Mass", "m(D_{s}K#pi#pi) candidate", 5200., 6000.,"MeV/c^{2}");
RooRealVar EventWeight("EventWeight", "EventWeight", 0., 1.);

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


//Load file
TFile* file;
file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/Ds3piBKgShape.root");
TTree* tree = (TTree*) file->Get("DecayTree");
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_Mass",1);
tree->SetBranchStatus("EventWeight",1);

//Fill needed variable in RooDataSet
RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_Mass,EventWeight),Import(*tree),WeightVar(EventWeight));

//Fit
RooFitResult *result;
result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3),SumW2Error(kTRUE));
cout << "result is --------------- "<<endl;
result->Print(); 

//fill array to return the fit values
fitValues_BGShapethreePi[0] = mean1.getVal();
fitValues_BGShapethreePi[1] = mean2.getVal();
fitValues_BGShapethreePi[2] = a1.getVal();
fitValues_BGShapethreePi[3] = a2.getVal();
fitValues_BGShapethreePi[4] = n1.getVal();
fitValues_BGShapethreePi[5] = n2.getVal();
fitValues_BGShapethreePi[6] = sigma1.getVal();
fitValues_BGShapethreePi[7] = sigma2.getVal();
fitValues_BGShapethreePi[8] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_Mass.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi.eps");

return fitValues_BGShapethreePi;
}

double *fitBGShapethreePiDstar(){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_Mass("Bs_Mass", "m(D_{s}*K#pi#pi)", 4900., 6200.,"MeV/c^{2}");
RooRealVar EventWeight("EventWeight","EventWeight", 0., 1.);

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


//Load file
TFile* file;
file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/Dstar3piBKgShape.root");
TTree* tree = (TTree*) file->Get("DecayTree");
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("EventWeight",1);
tree->SetBranchStatus("Bs_Mass",1);

//Fill needed variable in RooDataSet
RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_Mass,EventWeight),Import(*tree),WeightVar(EventWeight));

//Fit
RooFitResult *result;
result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3),SumW2Error(kTRUE));
cout << "result is --------------- "<<endl;
result->Print(); 

//fill array to return the fit values
fitValues_BGShapethreePiDstar[0] = mean1.getVal();
fitValues_BGShapethreePiDstar[1] = mean2.getVal();
fitValues_BGShapethreePiDstar[2] = a1.getVal();
fitValues_BGShapethreePiDstar[3] = a2.getVal();
fitValues_BGShapethreePiDstar[4] = n1.getVal();
fitValues_BGShapethreePiDstar[5] = n2.getVal();
fitValues_BGShapethreePiDstar[6] = sigma1.getVal();
fitValues_BGShapethreePiDstar[7] = sigma2.getVal();
fitValues_BGShapethreePiDstar[8] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_Mass.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.eps");

return fitValues_BGShapethreePiDstar;
}

double getExpBkgShape(string fileName){

TFile* file;
TTree* tree;
file= new TFile(fileName.c_str());
tree = (TTree*) file->Get("DecayTree");	
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("DTF_Bs_M",1);

RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} #pi #pi #pi)", 5500., 6000.,"MeV/c^{2}");
RooArgList list =  RooArgList(DTF_Bs_M);
RooDataSet* data = 0;
data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));

RooRealVar exp_par("exp_par","#lambda_{1}",-1.6508e-03,-10.,0.);	
RooExponential bkg_exp("bkg_exp1","exponential background1",DTF_Bs_M,exp_par);

RooFitResult *result;
result = bkg_exp.fitTo(*data,Save(kTRUE),Extended(kFALSE));

cout << "result is --------------- "<<endl;
result->Print();

TCanvas* c1= new TCanvas("");
RooPlot* frame_m= DTF_Bs_M.frame();
data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(60));
bkg_exp.plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
frame_m->Draw();
c1->Print("eps/ExpBkg_from_Sideband.eps");

return exp_par.getVal();

}

double *fitBDTNorm(){

        NamedParameter<int> DsKKpi("DsKKpi", 1);
	NamedParameter<int> fitSimultaneous("fitSimultaneous", 1);
	NamedParameter<int> sWeighting("sWeighting", 0);
        NamedParameter<int> makePlot("makePlot", 0);
        NamedParameter<int> BDTScan("BDTScan", 0);
        NamedParameter<int> fitCombined("fitCombined", 0);	
        NamedParameter<double> Year("Year", 2011);


   	gStyle->SetOptStat(0);
    	gStyle->SetTitleXSize(0.03);
    	gStyle->SetTitleYSize(0.03);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();

	///fit options
	bool BDTscan = false;
	bool sWeight = false;
	bool Ds2KKpi = true;
        bool sevenTeV = true;
	bool fitSimultan = true;
        bool combineYears = false;	


	if(BDTScan == 1) BDTscan = true;
        if(sWeighting == 1) sWeight = true;	
        if(DsKKpi == 0) Ds2KKpi = false;
        if(Year == 2012) sevenTeV = false;
        if(fitSimultaneous == 0) fitSimultan = false;
        if(fitCombined == 1) combineYears = true;	


	string filename11="/auto/data/kecke/B2DPiPiPi/Data2011/data11_Bs2Dspipipi_fullSelectionBDTG_DsK3fb_Selection.root";
	string filename12="/auto/data/kecke/B2DPiPiPi/Data2012/data12_Bs2Dspipipi_fullSelectionBDTG.root";

	///Load file

	///seperate fits for 7 & 8 TeV data (fitSimultan = false)
	TFile* file;
	TTree* tree;
	if(!fitSimultan){
		if(sevenTeV && Ds2KKpi) file= new TFile("/auto/data/kecke/B2DPiPiPi/forMaster/Data2011/data11_Bs2Dspipipi_fullSelectionBDTG.root");
		if(!sevenTeV && Ds2KKpi) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
		if((!sevenTeV) && (!Ds2KKpi)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data12_Bs2Dspipipi_Ds2pipipi_fullSelectionBDTG.root");
		if(sevenTeV && (!Ds2KKpi)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data11_Bs2Dspipipi_Ds2pipipi_fullSelectionBDTG.root");
		tree = (TTree*) file->Get("DecayTree");	
   		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("DTF_Bs_M",1);
	}

	///one simultaneous fit for 7 & 8 TeV data (fitSimultan = true) 
	TFile* file11;
	TFile* file12;
	TTree* tree11;
	TTree* tree12;
	if(fitSimultan && Ds2KKpi && (!combineYears)){
		file11 = new TFile(filename11.c_str());
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile(filename12.c_str());
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}
	if(fitSimultan && (!Ds2KKpi) && (!combineYears)){
		file11 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data11_Bs2Dspipipi_Ds2pipipi_fullSelectionBDTG.root");
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data12_Bs2Dspipipi_Ds2pipipi_fullSelectionBDTG.root");
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}

	/// combine years and fit Ds->KKpi & Ds->pipipi simultaneous
	TFile* fileDs2KKpi;
	TFile* fileDs2pipipi;
	TTree* treeDs2KKpi;
	TTree* treeDs2pipipi;
	if(fitSimultan && combineYears){
		fileDs2KKpi = new TFile("/auto/data/kecke/B2DPiPiPi/Bs2Dspipipi_Ds2KKpi_fullSelection_combined.root");
		treeDs2KKpi = (TTree*) fileDs2KKpi->Get("DecayTree");
		treeDs2KKpi->SetBranchStatus("*",0);
		treeDs2KKpi->SetBranchStatus("DTF_Bs_M",1);
		fileDs2pipipi = new TFile("/auto/data/kecke/B2DPiPiPi/Bs2Dspipipi_Ds2pipipi_fullSelection_combined.root");
		treeDs2pipipi = (TTree*) fileDs2pipipi->Get("DecayTree");
		treeDs2pipipi->SetBranchStatus("*",0);
		treeDs2pipipi->SetBranchStatus("DTF_Bs_M",1);
	}
	
	

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} #pi #pi #pi)", 4975., 5800.,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);

        RooDataSet* data = 0;
	RooDataSet* data11 = 0;
	RooDataSet* data12 = 0;
	RooDataSet* dataDs2KKpi = 0;
	RooDataSet* dataDs2pipipi = 0;
	if(!fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree));
		dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*tree));
		dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*tree));
	}
	if(fitSimultan && (!combineYears)){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree11));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree11));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree12));
		dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*tree11));
		dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*tree11));
	}
	if(fitSimultan && combineYears){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
		dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
		dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*treeDs2pipipi));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*treeDs2pipipi));
	}

	///define category to distinguish between 2011 and 2012 data
	RooCategory sample_year("sample_year","sample_year") ;
	sample_year.defineType("y11");
	sample_year.defineType("y12");

	///define category to distinguish between Ds->KKpi and Ds->pipipi
	RooCategory sample_Ds("sample_Ds","sample_Ds") ;
	sample_Ds.defineType("Ds_kaonkaonpion");
	sample_Ds.defineType("Ds_pionpionpion");

	RooDataSet* combData;
	if(fitSimultan && (!combineYears)) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_year),Import("y11",*data11),Import("y12",*data12));
	if(fitSimultan && combineYears) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_Ds),Import("Ds_kaonkaonpion",*dataDs2KKpi),Import("Ds_pionpionpion",*dataDs2pipipi));


	///import Bkg shapes from MC fits--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//2) Ds*pipipi in Dspipipi
	// 11 params: mean1, mean2, mean3, sigmaL1, sigmaR1, sigmaL2, sigmaR2, sigmaL3, sigmaR3, f_1, f_2
	//double fillarr2[11];
	double *DstarpipipiNorm;
	DstarpipipiNorm = fitBGShapeNorm(/*fillarr2*/);

	///Define fit model----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///Signal model - two double gaussians for B0 and Bs signal 

	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5366.7,5345.,5380.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 32.56);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 13.54,10.,20.);
	RooRealVar a1("a1","a1", -1.8098e+00);
	RooRealVar n1("n1","n1", 5.5870e+00);
	RooRealVar a2("a2","a2", 9.7197e-01);
	RooRealVar n2("n2","n2", 1.7471e+00);
	RooRealVar sigmaCB1("sigmaCB1", "B_{s} CB #sigma_{1}", 1.3120e+01,0.,50.56);
	RooRealVar sigmaCB2("sigmaCB2", "B_{s} CB #sigma_{2}", 10.,0.,20.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanBs1, sigmaCB1, a1, n1);
	RooCBShape CB2("CB2", "CB2", DTF_Bs_M, meanBs1, sigmaCB2, a2, n2);
	RooRealVar f_CB("f_CB" , "f__{B_{s}}", 7.0763e-01);
	RooAddPdf DoubleCBs("DoubleCBs", "DoubleCBs", RooArgList(CB1,CB2),RooArgList(f_CB));
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.18);
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_CB));
	RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));

	//signal pdf
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential for first comb. BG. component
	RooRealVar exp_par("exp_par","#lambda_{1}",-1.6508e-03,-10.,0.);	
	RooRealVar exp_par_11("exp_par_11","#lambda_{1} 11", getExpBkgShape(filename11.c_str()));
	RooRealVar exp_par_12("exp_par_12","#lambda_{1} 12", getExpBkgShape(filename12.c_str()));
	RooExponential bkg_exp("bkg_exp1","exponential background1",DTF_Bs_M,exp_par);
	RooExponential bkg_exp_11("bkg_exp_11","exponential background 11",DTF_Bs_M,exp_par_11);
	RooExponential bkg_exp_12("bkg_exp 12","exponential background 12",DTF_Bs_M,exp_par_12);

	///assign imported shapes

	//1) Ds*pipipi in Dspipipi-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	double mean1Variation_low = DstarpipipiNorm[0] - (0.05 * DstarpipipiNorm[0]);
	double mean1Variation_high = DstarpipipiNorm[0] + (0.05 * DstarpipipiNorm[0]);
	double mean2Variation_low = DstarpipipiNorm[1] - (0.05 * DstarpipipiNorm[1]);
	double mean2Variation_high = DstarpipipiNorm[1] + (0.05 * DstarpipipiNorm[1]);
	double mean3Variation_low = DstarpipipiNorm[2] - (0.01 * DstarpipipiNorm[2]);
	double mean3Variation_high = DstarpipipiNorm[2] + (0.01 * DstarpipipiNorm[2]);

	double sigmaL1Variation_low = DstarpipipiNorm[3] -  DstarpipipiNorm[3];
	double sigmaL1Variation_high = DstarpipipiNorm[3] + DstarpipipiNorm[3];
	double sigmaR1Variation_low = DstarpipipiNorm[4] - DstarpipipiNorm[4];
	double sigmaR1Variation_high = DstarpipipiNorm[4] + DstarpipipiNorm[4];
	double sigmaL2Variation_low = DstarpipipiNorm[5] - DstarpipipiNorm[5];
	double sigmaL2Variation_high = DstarpipipiNorm[5] + DstarpipipiNorm[5];
	double sigmaR2Variation_low = DstarpipipiNorm[6] -  DstarpipipiNorm[6];
	double sigmaR2Variation_high = DstarpipipiNorm[6] +  DstarpipipiNorm[6];
	double sigmaL3Variation_low = DstarpipipiNorm[7] - DstarpipipiNorm[7];
	double sigmaL3Variation_high = DstarpipipiNorm[7] +  DstarpipipiNorm[7];
	double sigmaR3Variation_low = DstarpipipiNorm[8] - DstarpipipiNorm[8];
	double sigmaR3Variation_high = DstarpipipiNorm[8] + DstarpipipiNorm[8];

	RooRealVar* sigmaL1;
	RooRealVar* sigmaR1;
	RooRealVar* sigmaL2;
	RooRealVar* sigmaR2;
	RooRealVar* sigmaL3;
	RooRealVar* sigmaR3;

	RooRealVar* f_1;
	RooRealVar* f_2;

	///mean of gaussians
	RooRealVar* mean1;
	RooRealVar* mean2;
	RooRealVar* mean3;

	///width and fractions of gaussians 2011

	if(sevenTeV){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0]);//, mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1]);//, mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2]);//, mean3Variation_low, mean3Variation_high);
	}

	///fractions of gauss functions for 2011

	if(sevenTeV){
		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	///fractions of gauss functions for 2012

	if(!sevenTeV){
		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	///width and fractions of gaussians 2012

	if(!sevenTeV){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0]);//, mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1]);//, mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2]);//, mean3Variation_low, mean3Variation_high);
	}


	if(fitSimultan && Ds2KKpi && (!combineYears)){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0]);//, mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1]);//, mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2]);//, mean3Variation_low, mean3Variation_high);
	}

	if(fitSimultan && (!Ds2KKpi) && (!combineYears)){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 3*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0]);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1]);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2]);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	if(combineYears){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 3*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", DTF_Bs_M, *mean1, *sigmaL1,*sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", DTF_Bs_M, *mean2, *sigmaL2,*sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", DTF_Bs_M, *mean3, *sigmaL3,*sigmaR3);

	//add functions
	RooAddPdf Dstarpipipi_as_Dspipipi("Dstarpipipi_as_Dspipipi", "Dstarpipipi_as_Dspipipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(*f_1,*f_2));

	//sum all background shapes
	//yields
	RooRealVar N_Dstarpipipi("N_Dstarpipipi","N_Dstarpipipi", data->numEntries()/2., 0., data->numEntries());
	RooRealVar N_comb("N_comb","N_comb", data->numEntries()/2., 0., data->numEntries());

        //sum background pdfs for BDT estimation 
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar f_bg("f_bg","f_bg", 0.5,0.,1.);
	RooAddPdf bkg("bkg", "bkg", RooArgList(Dstarpipipi_as_Dspipipi, bkg_exp), RooArgList(f_bg));

	///add yields for possible simultaneous fit
	RooRealVar N_Bs_11("N_Bs_11", "2011 #B_{s}", data11->numEntries()/2., 0., data11->numEntries());
	RooRealVar N_Bs_12("N_Bs_12", "2012 #B_{s}", data12->numEntries()/2., 0., data12->numEntries());
	RooRealVar N_Dstarpipipi_11("N_Dstarpipipi_11","N_Dstarpipipi_11", data11->numEntries()/2., 0., data11->numEntries());
	RooRealVar N_Dstarpipipi_12("N_Dstarpipipi_12","N_Dstarpipipi_12", data12->numEntries()/2., 0., data12->numEntries());
	RooRealVar N_comb_11("N_comb_11","N_comb_11", data11->numEntries()/2., 0., data11->numEntries());
	RooRealVar N_comb_12("N_comb_12","N_comb_12", data12->numEntries()/2., 0., data12->numEntries());


	RooRealVar N_Bs_Ds2KKpi("N_Bs_Ds2KKpi", "Ds2KKpi #B_{s}", dataDs2KKpi->numEntries()/2., 0., dataDs2KKpi->numEntries());
	RooRealVar N_Bs_Ds2pipipi("N_Bs_Ds2pipipi", "20Ds2pipipi #B_{s}", dataDs2pipipi->numEntries()/2., 0., dataDs2pipipi->numEntries());
	RooRealVar N_Dstarpipipi_Ds2KKpi("N_Dstarpipipi_Ds2KKpi","N_Dstarpipipi_Ds2KKpi", dataDs2KKpi->numEntries()/2., 0., dataDs2KKpi->numEntries());
	RooRealVar N_Dstarpipipi_Ds2pipipi("N_Dstarpipipi_Ds2pipipi","N_Dstarpipipi_Ds2pipipi", dataDs2pipipi->numEntries()/2., 0., dataDs2pipipi->numEntries());
	RooRealVar N_comb_Ds2KKpi("N_comb_Ds2KKpi","N_comb_Ds2KKpi", dataDs2KKpi->numEntries()/2., 0., dataDs2KKpi->numEntries());
	RooRealVar N_comb_Ds2pipipi("N_comb_Ds2pipipi","N_comb_Ds2pipipi", dataDs2pipipi->numEntries()/2., 0., dataDs2pipipi->numEntries());

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///total pdf
	///----------------------

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(DoubleCBs, Dstarpipipi_as_Dspipipi ,bkg_exp), RooArgList(N_Bs, N_Dstarpipipi, N_comb));
	RooAbsPdf* pdfScan=new RooAddPdf("pdfScan", "pdfScan", RooArgList(DoubleGaussBs, bkg), RooArgList(N_Bs, n_bkg));

	///construct a simultaneous pdf using category sample as index
	RooSimultaneous* simPdf;
	RooAbsPdf* pdf11;
	RooAbsPdf* pdf12;
	RooAbsPdf* pdfDs2KKpi;
	RooAbsPdf* pdfDs2pipipi;
	if(fitSimultan && (!combineYears)){
		simPdf = new RooSimultaneous("simPdf","simultaneous pdf",sample_year);
		pdf11=new RooAddPdf("pdf11", "pdf11", RooArgList(DoubleCBs, Dstarpipipi_as_Dspipipi ,bkg_exp_11), RooArgList(N_Bs_11, N_Dstarpipipi_11, N_comb_11));
		pdf12=new RooAddPdf("pdf12", "pdf12", RooArgList(DoubleCBs, Dstarpipipi_as_Dspipipi ,bkg_exp_12), RooArgList(N_Bs_12, N_Dstarpipipi_12, N_comb_12));
		simPdf->addPdf(*pdf11,"y11");
		simPdf->addPdf(*pdf12,"y12");
	}
	if(fitSimultan && combineYears){
		simPdf = new RooSimultaneous("simPdf","simultaneous pdf",sample_Ds);
		pdfDs2KKpi=new RooAddPdf("pdfDs2KKpi", "pdfDs2KKpi", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_11), RooArgList(N_Bs_Ds2KKpi, N_Dstarpipipi_Ds2KKpi, N_comb_Ds2KKpi));
		pdfDs2pipipi=new RooAddPdf("pdfDs2pipipi", "pdfDs2pipipi", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_12), RooArgList(N_Bs_Ds2pipipi, N_Dstarpipipi_Ds2pipipi, N_comb_Ds2pipipi));
		simPdf->addPdf(*pdfDs2KKpi,"Ds_kaonkaonpion");
		simPdf->addPdf(*pdfDs2pipipi,"Ds_pionpionpion");
	}

	///Fit
	RooFitResult *result;
	if(!fitSimultan){
		if(!BDTscan) result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
		if(BDTscan) result = pdfScan->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
	}

	if(fitSimultan) result = simPdf->fitTo(*combData,Save(kTRUE),Extended(kTRUE),NumCPU(3));

	cout << "result is --------------- "<<endl;
	result->Print();

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= DTF_Bs_M.frame();
	RooPlot* frame_m_11= DTF_Bs_M.frame();
	RooPlot* frame_m_12= DTF_Bs_M.frame();
	frame_m->SetTitle("");
	frame_m_11->SetTitle("");
	frame_m_12->SetTitle("");

	string BsMassDistribution = "eps/Final/3pi_BmassFit";
	string BsMassResidual = "eps/Final/3pi_residual";
	string BsMassPull = "eps/Final/3pi_pull";
	string BsMassSweight = "eps/Final/3pi_Bs_sWeight";
	string simfit11 = "_sim11.eps";
	string simfit12 = "_sim12.eps";
	string StringDs2pipipi = "_Ds2pipipi";
	string StringDs2KKpi = "_Ds2KKpi";
	string eps = ".eps";

	string eleven = "_11.eps";
	string twelve = "_12.eps";

	string BsMassDistribution11 = "eps/Final/3pi_BmassFit";
	string BsMassDistribution12 = "eps/Final/3pi_BmassFit";
	string BsMassResidual11 = "eps/Final/3pi_residual";
	string BsMassResidual12 = "eps/Final/3pi_residual";
	string BsMassPull11 = "eps/Final/3pi_pull";
	string BsMassPull12 = "eps/Final/3pi_pull";
	string BsMassSweight11 = "eps/Final/3pi_Bs_sWeight";
	string BsMassSweight12 = "eps/Final/3pi_Bs_sWeight";

	// append Ds2pipipi string to all names in case this channel is fitted
	if(!Ds2KKpi){
		BsMassDistribution.append(StringDs2pipipi);
		BsMassResidual.append(StringDs2pipipi);
		BsMassPull.append(StringDs2pipipi);
		BsMassSweight.append(StringDs2pipipi);
	}

	// 2011 results
	if(sevenTeV && (!fitSimultan)){
		BsMassDistribution.append(eleven);
		BsMassResidual.append(eleven);
		BsMassPull.append(eleven);
		BsMassSweight.append(eleven);
	}

	//2012 results
	if(!sevenTeV && (!fitSimultan)){
		if(!Ds2KKpi) BsMassDistribution.append(StringDs2pipipi);
		BsMassDistribution.append(twelve);
		BsMassResidual.append(twelve);
		BsMassPull.append(twelve);
		BsMassSweight.append(twelve);
	}

	if(fitSimultan && Ds2KKpi && (!combineYears)){
	BsMassDistribution11.append(simfit11);
	BsMassResidual11.append(simfit11);
	BsMassPull11.append(simfit11);
	BsMassSweight11.append(simfit11);

	BsMassDistribution12.append(simfit12);
	BsMassResidual12.append(simfit12);
	BsMassPull12.append(simfit12);
	BsMassSweight12.append(simfit12);
	}

	if(fitSimultan && (!Ds2KKpi) && (!combineYears)){
	BsMassDistribution11.append(StringDs2pipipi);
	BsMassResidual11.append(StringDs2pipipi);
	BsMassPull11.append(StringDs2pipipi);
	BsMassSweight11.append(StringDs2pipipi);
	BsMassDistribution11.append(simfit11);
	BsMassResidual11.append(simfit11);
	BsMassPull11.append(simfit11);
	BsMassSweight11.append(simfit11);

	BsMassDistribution12.append(StringDs2pipipi);
	BsMassResidual12.append(StringDs2pipipi);
	BsMassPull12.append(StringDs2pipipi);
	BsMassSweight12.append(StringDs2pipipi);
	BsMassDistribution12.append(simfit12);
	BsMassResidual12.append(simfit12);
	BsMassPull12.append(simfit12);
	BsMassSweight12.append(simfit12);
	}

	if(fitSimultan && combineYears){
	BsMassDistribution12.append(StringDs2pipipi);
	BsMassResidual12.append(StringDs2pipipi);
	BsMassPull12.append(StringDs2pipipi);
	BsMassSweight12.append(StringDs2pipipi);
	BsMassDistribution12.append(eps);
	BsMassResidual12.append(eps);
	BsMassPull12.append(eps);
	BsMassSweight12.append(eps);

	BsMassDistribution11.append(StringDs2KKpi);
	BsMassResidual11.append(StringDs2KKpi);
	BsMassPull11.append(StringDs2KKpi);
	BsMassSweight11.append(StringDs2KKpi);
	BsMassDistribution11.append(eps);
	BsMassResidual11.append(eps);
	BsMassPull11.append(eps);
	BsMassSweight11.append(eps);

	}


	if((!BDTscan) && (!fitSimultan)){	

		data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(60));
		pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
		pdf->plotOn(frame_m,Components(DoubleGaussBs),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(Dstarpipipi_as_Dspipipi),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(bkg_exp),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution.c_str());
	}

	if((!BDTscan) && fitSimultan && (!combineYears)){
		combData->plotOn(frame_m_11,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(120));
		simPdf->plotOn(frame_m_11,Name("pdf11"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleCBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution11.c_str());

		combData->plotOn(frame_m_12,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(120));
		simPdf->plotOn(frame_m_12,Name("pdf12"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleCBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution12.c_str());
	}

	if((!BDTscan) && fitSimultan && combineYears){
		combData->plotOn(frame_m_11,Name("dataDs2KKpi"),Cut("sample_Ds==sample_Ds::Ds_kaonkaonpion"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_11,Name("pdfDs2KKpi"),Slice(sample_Ds,"Ds_kaonkaonpion"),ProjWData(sample_Ds,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(DoubleGaussBs),ProjWData(sample_Ds,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_Ds,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(bkg_exp_11),ProjWData(sample_Ds,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution11.c_str());

		combData->plotOn(frame_m_12,Name("dataDs2pipipi"),Cut("sample_Ds==sample_Ds::Ds_pionpionpion"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_12,Name("pdfDs2pipipi"),Slice(sample_Ds,"Ds_pionpionpion"),ProjWData(sample_Ds,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(DoubleGaussBs),ProjWData(sample_Ds,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_Ds,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(bkg_exp_12),ProjWData(sample_Ds,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution12.c_str());
	}


	///fit results
	double chi2 = 0;
	if(!fitSimultan) chi2 = frame_m->chiSquare("pdf","data",15);
	if(fitSimultan && (!combineYears)) chi2 = (frame_m_11->chiSquare("pdf11","data11",15) + frame_m_12->chiSquare("pdf12","data12",15)) / 2;
	if(fitSimultan && combineYears) chi2 = (frame_m_11->chiSquare("pdfDs2KKpi","dataDs2KKpi",15) + frame_m_12->chiSquare("pdfDs2pipipi","dataDs2pipipi",15)) / 2;
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid;
	RooHist* hresid11;
	RooHist* hresid12;
	RooHist* hpull;
	RooHist* hpull11;
	RooHist* hpull12;

	if(!fitSimultan){
		hresid = frame_m->residHist("data","pdf") ;
  		hpull = frame_m->pullHist("data","pdf") ;
	}
	if(fitSimultan && (!combineYears)){
		hresid11 = frame_m_11->residHist("data11","pdf11") ;
  		hpull11 = frame_m_11->pullHist("data11","pdf11") ;
		hresid12 = frame_m_12->residHist("data12","pdf12") ;
  		hpull12 = frame_m_12->pullHist("data12","pdf12") ;
	}
	if(fitSimultan && combineYears){
		hresid11 = frame_m_11->residHist("dataDs2KKpi","pdfDs2KKpi") ;
  		hpull11 = frame_m_11->pullHist("dataDs2KKpi","pdfDs2KKpi") ;
		hresid12 = frame_m_12->residHist("dataDs2pipipi","pdfDs2pipipi") ;
  		hpull12 = frame_m_12->pullHist("dataDs2pipipi","pdfDs2pipipi") ;
	}

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = DTF_Bs_M.frame(Title("Residual Distribution")) ;
	RooPlot* frame2_11 = DTF_Bs_M.frame(Title("Residual Distribution 2011"));
	RooPlot* frame2_12 = DTF_Bs_M.frame(Title("Residual Distribution 2012"));

	if((!BDTscan) && (!fitSimultan)){	
		frame2->SetTitle("");
		frame2->addPlotable(hresid,"P") ;
		frame2->Draw();
		if(makePlot == 1) c1->Print(BsMassResidual.c_str());
	}
	if((!BDTscan) && fitSimultan){
		frame2_11->SetTitle("");
		frame2_11->addPlotable(hresid11,"P") ;
		frame2_11->Draw();
		if(makePlot == 1) c1->Print(BsMassResidual11.c_str());

		frame2_12->SetTitle("");
		frame2_12->addPlotable(hresid12,"P") ;
		frame2_12->Draw();
		if(makePlot == 1) c1->Print(BsMassResidual12.c_str());
	}

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = DTF_Bs_M.frame(Title("Pull Distribution")) ;
	RooPlot* frame3_11 = DTF_Bs_M.frame(Title("Pull Distribution 2011")) ;
	RooPlot* frame3_12 = DTF_Bs_M.frame(Title("Pull Distribution 2012")) ;

	if((!BDTscan) && (!fitSimultan)){	
		frame3->SetTitle("");	
		frame3->SetLabelFont(62,"Y");
		frame3->addPlotable(hpull,"P") ;
		frame3->Draw();
		if(makePlot == 1) c1->Print(BsMassPull.c_str());
	}
	if((!BDTscan) && fitSimultan){
		frame3_11->SetTitle("");	
		frame3_11->SetLabelFont(62,"Y");
		frame3_11->addPlotable(hpull11,"P") ;
		frame3_11->SetFillColor(kBlue);
		frame3_11->Draw("hist");
		if(makePlot == 1) c1->Print(BsMassPull11.c_str());

		frame3_12->SetTitle("");	
		frame3_12->SetLabelFont(62,"Y");
		frame3_12->addPlotable(hpull12,"P") ;
		frame3_12->SetFillColor(kBlue);
		frame3_12->Draw("hist");
		if(makePlot == 1) c1->Print(BsMassPull12.c_str());
	}


	///save fitavalues in array to be use in fitBDT()
	fitvalues_BDTNorm[0] = mean1->getVal();
	fitvalues_BDTNorm[1] = mean2->getVal();
	fitvalues_BDTNorm[2] = mean3->getVal();
	fitvalues_BDTNorm[3] = sigmaL1->getVal();
	fitvalues_BDTNorm[4] = sigmaR1->getVal();
	fitvalues_BDTNorm[5] = sigmaL2->getVal();
	fitvalues_BDTNorm[6] = sigmaR2->getVal();
	fitvalues_BDTNorm[7] = sigmaL3->getVal();
	fitvalues_BDTNorm[8] = sigmaR3->getVal();
	fitvalues_BDTNorm[9] = f_1->getVal();
	fitvalues_BDTNorm[10] = f_2->getVal();

	if(!fitSimultan){ 
		fitvalues_BDTNorm[11] = N_Dstarpipipi.getVal();
		fitvalues_BDTNorm[12] = N_Bs.getVal();
		fitvalues_BDTNorm[13] = 0;
		fitvalues_BDTNorm[14] = 0;
	}
	if(fitSimultan && (!combineYears)){
		fitvalues_BDTNorm[11] = N_Dstarpipipi_11.getVal();
		fitvalues_BDTNorm[12] = N_Bs_11.getVal();
		fitvalues_BDTNorm[13] = N_Dstarpipipi_12.getVal();
		fitvalues_BDTNorm[14] = N_Bs_12.getVal();
	}
	if(fitSimultan && combineYears){
		fitvalues_BDTNorm[11] = N_Dstarpipipi_Ds2KKpi.getVal();
		fitvalues_BDTNorm[12] = N_Bs_Ds2KKpi.getVal();
		fitvalues_BDTNorm[13] = N_Dstarpipipi_Ds2pipipi.getVal();
		fitvalues_BDTNorm[14] = N_Bs_Ds2pipipi.getVal();
	}

if(sWeight){

		sigmaL1->setConstant();
 		sigmaR1->setConstant();
 		sigmaL2->setConstant();
 		sigmaR2->setConstant();
 		sigmaL3->setConstant();
 		sigmaR3->setConstant();
		sigmaCB1.setConstant();
		sigmaCB2.setConstant();
 		f_1->setConstant();
 		f_2->setConstant();
 		f_CB.setConstant();
		mean1->setConstant();
		mean2->setConstant();
		mean3->setConstant();
		meanBs1.setConstant();
		sigmaBs1.setConstant();
		sigmaBs2.setConstant();
		f_GaussBs.setConstant();
		exp_par_11.setConstant();
		exp_par_12.setConstant();

	
		SPlot* sData12 = new SPlot("sData12","An SPlot for 2012 data",*data12,pdf12,RooArgList(N_Bs_12, N_Dstarpipipi_12, N_comb_12));
		SPlot* sData11 = new SPlot("sData11","An SPlot for 2011 data",*data11,pdf11,RooArgList(N_Bs_11, N_Dstarpipipi_11, N_comb_11)); 

		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs11 = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist11 = (TH2*)data11->createHistogram("DTF_Bs_M,N_Bs_11_sw");
		SwDsHist11->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist11->SetTitle("");
		SwDsHist11->Draw();
		if(makePlot == 1) SwDs11->Print(BsMassSweight11.c_str());

		TCanvas* SwDs12 = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist12 = (TH2*)data12->createHistogram("DTF_Bs_M,N_Bs_12_sw");
		SwDsHist12->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist12->SetTitle("");
		SwDsHist12->Draw();
		if(makePlot == 1) SwDs12->Print(BsMassSweight12.c_str());


    		///Create output file
   		 TFile* output_11;
   		 TFile* output_12;
		 if(Ds2KKpi){
			output_11 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_final_sweight.root","RECREATE");
			output_12 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root","RECREATE");
		}
		 if(!Ds2KKpi){
			output_11 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_Ds2pipipi_11_final_sweight.root","RECREATE");
			output_12 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_Ds2pipipi_12_final_sweight.root","RECREATE");
		}

		 tree11->SetBranchStatus("*",1);
   		 TTree* new_tree11 = tree11->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
    		 double w_11;
    		 TBranch* Bra_sw11 = new_tree11->Branch("N_Bs_sw", &w_11, "N_Bs_sw/D");

		 tree12->SetBranchStatus("*",1);
   		 TTree* new_tree12 = tree12->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
    		 double w_12;
    		 TBranch* Bra_sw12 = new_tree12->Branch("N_Bs_sw", &w_12, "N_Bs_sw/D");



  		  ///loop over events
    		  int numEvents_11 = new_tree11->GetEntries();
    		  int numEvents_12 = new_tree12->GetEntries();

    		  for(int i=0; i< numEvents_11; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_11 << endl;
			tree11->GetEntry(i);
			w_11=sData11->GetSWeight(i,"N_Bs_11_sw");
			Bra_sw11->Fill();
  		  }
	cout << "loop finished!!!" << endl;
   		 new_tree11->Write();
   		 output_11->Close();

    		  for(int i=0; i< numEvents_12; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_12 << endl;
			tree12->GetEntry(i);
			w_12=sData12->GetSWeight(i,"N_Bs_12_sw");
			Bra_sw12->Fill();
  		  }
	cout << "loop finished!!!" << endl;
   		 new_tree12->Write();
   		 output_12->Close();
	}


return fitvalues_BDTNorm;

}

void fitBDT(){

        NamedParameter<int> DsKKpi("DsKKpi", 1);
	NamedParameter<int> fitSimultaneous("fitSimultaneous", 1);
	NamedParameter<int> sWeighting("sWeighting", 0);
        NamedParameter<int> makePlot("makePlot", 0);
        NamedParameter<int> fitCombined("fitCombined", 0);	
        NamedParameter<double> Year("Year", 2011);

	bool sWeight=false;
	bool sevenTeV=true;
	bool fitSimultan=true;
	bool Ds2KKpi = true;
        bool combineYears = false;

        if(sWeighting == 1) sWeight = true;	
        if(DsKKpi == 0) Ds2KKpi = false;
        if(Year == 2012) sevenTeV = false;
        if(fitSimultaneous == 0) fitSimultan = false;
        if(fitCombined == 1) combineYears = true;

   	gStyle->SetOptStat(0);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();

	string filename11="/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG.root";
	string filename12="/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG.root";

	
	///seperate fits for 7 & 8 TeV data (fitSimultan = false)
	TFile* file;
	TTree* tree;
	if(!fitSimultan && Ds2KKpi){
		if(sevenTeV) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root");	
		if(!sevenTeV) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree = (TTree*) file->Get("DecayTree");	
   		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("DTF_Bs_M",1);
	}

	///one simultaneous fit for 7 & 8 TeV data (fitSimultan = true) 
	TFile* file11;
	TFile* file12;
	TTree* tree11;
	TTree* tree12;
	if(fitSimultan && Ds2KKpi && (!combineYears)){
		file11 = new TFile(filename11.c_str());
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile(filename12.c_str());
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}
	if(fitSimultan && (!Ds2KKpi) && (!combineYears)){
		file11 = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2pipipi_fullSelectionBDTG.root");
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_fullSelectionBDTG.root");
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}

	/// combine years and fit Ds->KKpi & Ds->pipipi simultaneous
	TFile* fileDs2KKpi;
	TFile* fileDs2pipipi;
	TTree* treeDs2KKpi;
	TTree* treeDs2pipipi;
	if(fitSimultan && combineYears){
		fileDs2KKpi = new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_Ds2KKpi_fullSelection_combined.root");
		treeDs2KKpi = (TTree*) fileDs2KKpi->Get("DecayTree");
		treeDs2KKpi->SetBranchStatus("*",0);
		treeDs2KKpi->SetBranchStatus("DTF_Bs_M",1);
		fileDs2pipipi = new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_Ds2pipipi_fullSelection_combined.root");
		treeDs2pipipi = (TTree*) fileDs2pipipi->Get("DecayTree");
		treeDs2pipipi->SetBranchStatus("*",0);
		treeDs2pipipi->SetBranchStatus("DTF_Bs_M",1);
	}

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} K #pi #pi)", 4975., 5800.,"MeV");
	RooArgList list =  RooArgList(DTF_Bs_M);

        RooDataSet* data = 0;
	RooDataSet* data11 = 0;
	RooDataSet* data12 = 0;
        RooDataSet* dataDs2KKpi = 0;
        RooDataSet* dataDs2pipipi = 0;

	if(!fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree));
                dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*tree));
                dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*tree));

	}
	if(fitSimultan && (!combineYears)){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree11));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree11));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree12));
                dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*tree11));
                dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*tree11));

	}
        if(fitSimultan && combineYears){
                data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
                dataDs2KKpi = new RooDataSet("dataDs2KKpi","dataDs2KKpi",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
                dataDs2pipipi = new RooDataSet("dataDs2pipipi","dataDs2pipipi",RooArgSet(DTF_Bs_M),Import(*treeDs2pipipi));
                data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*treeDs2KKpi));
                data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*treeDs2pipipi));
        }

	///define category to distinguish between 2011 and 2012 data
	RooCategory sample_year("sample_year","sample_year") ;
	sample_year.defineType("y11");
	sample_year.defineType("y12");

	RooDataSet* combData;
	if(fitSimultan && (!combineYears)) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_year),Import("y11",*data11),Import("y12",*data12));
	if(fitSimultan && combineYears) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_year),Import("y11",*dataDs2KKpi),Import("y12",*dataDs2pipipi));

	///import Bkg shapes from Normalization fit and MC---------------------------------------------------------------------------------------------------------------------------------------------------------------



	//4) Ds*Kpipi in DsKpipi
	// 13 params: mean1-3, sigma1-6, f_1 & f_2, N_Dstar, N_Bs
	//double fillarr4[13];
	double *DsstarKpipifromNorm;
	DsstarKpipifromNorm = fitBDTNorm();

	//5) Dspipipi in DsKpipi
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
	//double fillarr5[9];
	double *DspipipiSig;
	DspipipiSig = fitBGShapethreePi(/*fillarr5*/);

	//6) Ds*pipipi in DsKpipi
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
	//double fillarr6[9];
	double *DstarpipipiSig;
	DstarpipipiSig = fitBGShapethreePiDstar(/*fillarr6*/);

	///Define fit model----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///Signal model - two double gaussians for B0 and Bs signal 

	//B0 signal shape
	RooRealVar meanB01("meanB01", "B_{0} #mu", 5281.71,5250.,5320.); 
	RooRealVar sigmaB1("sigmaB1", "B_{0} #sigma_{1}",  32.56,5.,200.);	
	RooRealVar sigmaB2("sigmaB2", "B_{0} sigma_{2}", 13.,10.,50.);
	RooGaussian GaussB01("GaussB01", "GaussB01", DTF_Bs_M, meanB01, sigmaB1);
	RooGaussian GaussB02("GaussB02", "GaussB02", DTF_Bs_M, meanB01, sigmaB2);
	RooRealVar f_GaussB("f_GaussB" , "f__{B_{0}}", 0.066);//, 0.2, 0.75);
        RooAddPdf DoubleGaussB0("DoubleGaussB0", "DoubleGaussB0", RooArgList(GaussB01,GaussB02),RooArgList(f_GaussB));

	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 1.4186e+01);//,10.,90.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 15.,10.,55.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaB1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0.2, 0.75);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussB));

	//double Crystal Ball Shape
	RooRealVar a1_Sig("a1_Sig","a1_Sig", -2.1977e+00);//,-10.,10.);
	RooRealVar n1_Sig("n1_Sig","n1_Sig", 3.1187e+00);//,0.,100.);
	RooRealVar a2_Sig("a2_Sig","a2_Sig", 1.9229e+00);
	RooRealVar n2_Sig("n2_Sig","n2_Sig", 9.0758e-01);
	RooRealVar f_CB_Sig("f_CB_Sig","f_CB_Sig",5.5008e-01);
	RooCBShape CB1_Sig("CB1_Sig", "CB1_Sig", DTF_Bs_M, meanBs1, sigmaBs1, a1_Sig, n1_Sig);
	RooCBShape CB2_Sig("CB2_Sig", "CB2_Sig", DTF_Bs_M, meanBs1, sigmaBs2, a2_Sig, n2_Sig);
	RooAddPdf DoubleCBBs("DoubleCBBs", "DoubleCBBs", RooArgList(CB1_Sig,CB2_Sig),RooArgList(f_CB_Sig));

	//signal pdf
	RooRealVar N_B0("N_B0", "#B_{0}", 500, 0, 5000);
	RooRealVar N_Bs("N_Bs", "#B_{s}", 500, 0, 5000);
	RooRealVar N_B0_11("N_B0_11", "#B_{0} 11", 500, 0, 5000);
	RooRealVar N_Bs_11("N_Bs_11", "#B_{s} 11", 500, 0, 5000);
	RooRealVar N_B0_12("N_B0_12", "#B_{0} 12", 500, 0, 5000);
	RooRealVar N_Bs_12("N_Bs_12", "#B_{s} 12", 500, 0, 5000);
//	RooAddPdf signal("signal", "signal", RooArgList(DoubleGaussB0, DoubleGaussBs),RooArgList(N_B0, N_Bs));	

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential
        RooRealVar exp_par("exp_par","#lambda_{1}",-1.6508e-03,-10.,0.);
        RooRealVar exp_par_11("exp_par_11","#lambda_{1} 11",-1.6508e-03,-10.,0.);
        RooRealVar exp_par_12("exp_par_12","#lambda_{1} 12",-1.6508e-03,-10.,0.);
        RooExponential bkg_exp("bkg_exp1","exponential background1",DTF_Bs_M,exp_par);
        RooExponential bkg_exp_11("bkg_exp_11","exponential background 11",DTF_Bs_M,exp_par_11);
        RooExponential bkg_exp_12("bkg_exp 12","exponential background 12",DTF_Bs_M,exp_par_12);


	///assign imported shapes

	//1) Ds*Kpipi in DsKpipi-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//fix shape from dspipipi normalization fit

	RooRealVar* mean1;
	RooRealVar* mean2;
	RooRealVar* mean3;

	RooRealVar* sigmaL1;
	RooRealVar* sigmaR1;
	RooRealVar* sigmaL2;
	RooRealVar* sigmaR2;
	RooRealVar* sigmaL3;
	RooRealVar* sigmaR3;	

	RooRealVar* f_1;
	RooRealVar* f_2;

	///mean of gaussians from normalization fit

	mean1 = new RooRealVar("mean1","mu", DsstarKpipifromNorm[0]);
	mean2 = new RooRealVar("mean2","mu", DsstarKpipifromNorm[1]);
	mean3 = new RooRealVar("mean3","mu", DsstarKpipifromNorm[2]);


	///width of gaussians from normalization fit
	
	sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1", DsstarKpipifromNorm[3]);
	sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1", DsstarKpipifromNorm[4]);
	sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2", DsstarKpipifromNorm[5]);
	sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2", DsstarKpipifromNorm[6]);
	sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3", DsstarKpipifromNorm[7]);
	sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3", DsstarKpipifromNorm[8]);
	

	///fractions of gauss functions for 2011 data from normalization fit

	f_1 = new RooRealVar("f_{1}", "fraction1", DsstarKpipifromNorm[9]);
	f_2 = new RooRealVar("f_{2}", "fraction2", DsstarKpipifromNorm[10]);

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", DTF_Bs_M, *mean1, *sigmaL1,*sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", DTF_Bs_M, *mean2, *sigmaL2,*sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", DTF_Bs_M, *mean3, *sigmaL3,*sigmaR3);

	//add functions
	RooAddPdf DstarKpipi_as_DsKpipi("DstarKpipi_as_DsKpipi", "DstarKpipi_as_DsKpipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(*f_1,*f_2));


	///same shape shifted down for m(B^0)
	///mean of gaussians shifted by m_Bs - m_B0
	RooRealVar mean1Shifted("mean1Shifted","mu", mean1->getVal() - 87.33 );
	RooRealVar mean2Shifted("mean2Shifted","mu", mean2->getVal() - 87.33 );
	RooRealVar mean3Shifted("mean3Shifted","mu", mean3->getVal() - 87.33 );

	RooBifurGauss BifGauss1Shifted("BifGauss1Shifted","BifGauss1Shifted", DTF_Bs_M, mean1Shifted, *sigmaL1,*sigmaR1);
	RooBifurGauss BifGauss2Shifted("BifGauss2Shifted","BifGauss2Shifted", DTF_Bs_M, mean2Shifted, *sigmaL2,*sigmaR2);
	RooBifurGauss BifGauss3Shifted("BifGauss3Shifted","BifGauss3Shifted", DTF_Bs_M, mean3Shifted, *sigmaL3,*sigmaR3);

	RooAddPdf DstarKpipi_as_DsKpipi_Shifted("DstarKpipi_as_DsKpipi_Shifted", "DstarKpipi_as_DsKpipi_Shifted", RooArgList(BifGauss1Shifted, BifGauss2Shifted, BifGauss3Shifted), RooArgList(*f_1,*f_2));


	//2) Dspipipi in DsKpipi--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	double meanDspipipi1Variation_low = DspipipiSig[0] - (0.00001 * DspipipiSig[0]);
	double meanDspipipi1Variation_high = DspipipiSig[0] + (0.01 * DspipipiSig[0]);
	double meanDspipipi2Variation_low = DspipipiSig[1] - (0.0001 * DspipipiSig[1]);
	double meanDspipipi2Variation_high = DspipipiSig[1] + (0.01 * DspipipiSig[1]);


	//mean of crrystal balls
	RooRealVar meanDspipipi1("meanDspipipi1","mu", DspipipiSig[0]);//, meanDspipipi1Variation_low, meanDspipipi1Variation_high);
	RooRealVar meanDspipipi2("meanDspipipi2","mu", DspipipiSig[1]);//, meanDspipipi2Variation_low, meanDspipipi2Variation_high);

	// asymmetry parameter of crystsal balls
	RooRealVar a1("a1","a1", DspipipiSig[2]);
	RooRealVar a2("a2","a2", DspipipiSig[3]);
	RooRealVar n1("n1","n1", DspipipiSig[4]);
	RooRealVar n2("n2","n2", DspipipiSig[5]);

	//sigma of crystal balls
	RooRealVar sigma1("sigma_{1}", "sigma1", DspipipiSig[6]);
	RooRealVar sigma2("sigma_{2}", "sigma2", DspipipiSig[7]);

	//crystal Balls
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanDspipipi1, sigma1, a1, n1);
	RooCBShape CB2("CB2", "CB2", DTF_Bs_M, meanDspipipi2, sigma2, a2, n2);

	//fraction of crystal balls
	RooRealVar f_3("f_{3}", "fraction3", DspipipiSig[8]);

	//add all functions
	RooAddPdf Dspipipi_as_DsKpipi("Dspipipi_as_DsKpipi", "Dspipipi_as_DsKpipi", RooArgList(CB1, CB2), RooArgList(f_3));


	//3) Ds*pipipi in DsKpipi-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
/*
	double meanDstarThreePi1Variation_low = DstarpipipiSig[0] - (0.01 * DstarpipipiSig[0]);
	double meanDstarThreePi1Variation_high = DstarpipipiSig[0] + (0.01 * DstarpipipiSig[0]);
	double meanDstarThreePi2Variation_low = DstarpipipiSig[1] - (0.01 * DstarpipipiSig[1]);
	double meanDstarThreePi2Variation_high = DstarpipipiSig[1] + (0.01 * DstarpipipiSig[1]);
*/

	//mean of crrystal balls
	RooRealVar meanDstarThreePi1("meanDstarThreePi1","mu", DstarpipipiSig[0]/*, meanDstarThreePi1Variation_low, meanDstarThreePi1Variation_high*/);
	RooRealVar meanDstarThreePi2("meanDstarThreePi2","mu", DstarpipipiSig[1]/*, meanDstarThreePi2Variation_low, meanDstarThreePi2Variation_high*/);

	// asymmetry parameter of crystsal balls
	RooRealVar aDstarThreePi1("aDstarThreePi1","a1", DstarpipipiSig[2]);
	RooRealVar aDstarThreePi2("aDstarThreePi2","a2", DstarpipipiSig[3]);
	RooRealVar nDstarThreePi1("nDstarThreePi1","n1", DstarpipipiSig[4]);
	RooRealVar nDstarThreePi2("nDstarThreePi2","n2", DstarpipipiSig[5]);


	//sigma of crystal balls
	RooRealVar sigmaDstarThreePi1("DstarThreePi sigma_{1}", "sigma1", DstarpipipiSig[6]);
	RooRealVar sigmaDstarThreePi2("DstarThreePi sigma_{2}", "sigma2", DstarpipipiSig[7]);

	//crystal Balls
	RooCBShape CBDstarThreePi1("CBDstarThreePi1", "CB1", DTF_Bs_M, meanDstarThreePi1, sigmaDstarThreePi1, aDstarThreePi1, nDstarThreePi1);
	RooCBShape CBDstarThreePi2("CBDstarThreePi2", "CB2", DTF_Bs_M, meanDstarThreePi2, sigmaDstarThreePi2, aDstarThreePi2, nDstarThreePi2);

	//fraction of crystal balls
	RooRealVar f_4("f_{4}", "fraction4", DstarpipipiSig[8]);

	//add all gaussians
	RooAddPdf Dstarpipipi_as_DsKpipi("Dstarpipipi_as_DsKpipi", "Dstarpipipi_as_DsKpipi", RooArgList(CBDstarThreePi1, CBDstarThreePi2), RooArgList(f_4));


	//sum all background shapes
	//yields
	RooRealVar N_DstarKpipi("N_DstarKpipi","N_DstarKpipi", 1275, 0, data->numEntries());
	RooRealVar N_DstarKpipiShifted("N_DstarKpipiShifted","N_DstarKpipiShifted", 1275, 0, data->numEntries());
	RooRealVar N_DstarKpipi_11("N_DstarKpipi_11","N_DstarKpipi 11", 1275, 0, data11->numEntries());

	RooRealVar *ratioDstar;
	if(Ds2KKpi) ratioDstar = new RooRealVar("ratioDstar","ratioDstar",1.,0.,10.);

	//fixed from Bs->Ds(->KKpi)Kpipi massfit
	if(!Ds2KKpi) ratioDstar = new RooRealVar("ratioDstar", "ratioDstar", 2.6771e+00);

	RooRealVar N_DstarKpipi_12("N_DstarKpipi_12","N_DstarKpipi 12", 1275, 0, data12->numEntries());
	RooGenericPdf N_DstarKpipiShifted_11("N_DstarKpipiShifted_11","N_DstarKpipiShifted 11", "@0*@1",RooArgList(*ratioDstar,N_DstarKpipi_11));
	RooGenericPdf N_DstarKpipiShifted_12("N_DstarKpipiShifted_12","N_DstarKpipiShifted 12", "@0*@1",RooArgList(*ratioDstar,N_DstarKpipi_12));

	RooRealVar* N_Dspipipi;
	RooRealVar* N_Dstarpipipi;
	RooRealVar* N_Dspipipi_11;
	RooRealVar* N_Dstarpipipi_11;
	RooRealVar* N_Dspipipi_12;
	RooRealVar* N_Dstarpipipi_12;


	N_Dspipipi = new RooRealVar("N_Dspipipi","N_Dspipipi", 0.0335 * DsstarKpipifromNorm[12]);//, 0., data->numEntries());
	N_Dstarpipipi = new RooRealVar("N_Dstarpipipi","N_Dstarpipipi", 0.036 * DsstarKpipifromNorm[11]);//, 0., data->numEntries());

	if(fitSimultan){
		N_Dstarpipipi_11 = new RooRealVar("N_Dstarpipipi_11","N_Dstarpipipi 11", 0.036 * DsstarKpipifromNorm[11]);
		N_Dspipipi_11 = new RooRealVar("N_Dspipipi_11","N_Dspipipi 11", 0.0335 * DsstarKpipifromNorm[12]);
		N_Dstarpipipi_12 = new RooRealVar("N_Dstarpipipi_12","N_Dstarpipipi 12", 0.036 * DsstarKpipifromNorm[13]);
		N_Dspipipi_12 = new RooRealVar("N_Dspipipi_12","N_Dspipipi 12", 0.0335 * DsstarKpipifromNorm[14]);
	}

	RooAddPdf bkg_model("bkg_model", "bkg_model", RooArgList(DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi), RooArgList(N_DstarKpipi, *N_Dspipipi, *N_Dstarpipipi));

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());
	RooRealVar N_comb("N_comb","N_comb",100, 0, data->numEntries());
	RooRealVar N_comb_11("N_comb_11","N_comb 11",100, 0, data11->numEntries());
	RooRealVar N_comb_12("N_comb_12","N_comb 12",100, 0, data12->numEntries());

cout << "put together pdfs" << endl;

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0, N_Bs, N_DstarKpipi, *N_Dspipipi, *N_Dstarpipipi, N_comb, N_DstarKpipiShifted));

cout << "put together pdfs 1.5" << endl;

        RooSimultaneous* simPdf;
        RooAbsPdf* pdf11;
        RooAbsPdf* pdf12;

        if(fitSimultan){
        simPdf = new RooSimultaneous("simPdf","simultaneous pdf",sample_year);	

	pdf11 = new RooAddPdf("pdf11", "pdf11", RooArgList(DoubleGaussB0, DoubleGaussBs/*DoubleCBBs*/, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_11, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11, *N_Dstarpipipi_11, N_comb_11, N_DstarKpipiShifted_11));

	pdf12 = new RooAddPdf("pdf12", "pdf12", RooArgList(DoubleGaussB0, DoubleGaussBs/*DoubleCBBs*/, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_12, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12));

	simPdf->addPdf(*pdf11,"y11");
	simPdf->addPdf(*pdf12,"y12");
	}
	///Fit
	RooFitResult *result;
	if(!fitSimultan)result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
        if(fitSimultan) result = simPdf->fitTo(*combData,Save(kTRUE),Extended(kTRUE),NumCPU(3));

	cout << "result is --------------- "<<endl;
	result->Print();

	//cout << "sigma2 =" << sigma1.getVal()*scale.getVal() << " pm " << sigma1.getError()*scale.getVal() << endl;

	///calculate # (signal)background events in signal region
	//cout << endl;
	//cout << endl;
/*
	Bplus_MM.setRange("signal",mean1.getVal()-60.,mean1.getVal()+60.);
	
	RooAbsReal* S_fr= signal.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t S = S_fr->getVal()*n_sig.getVal();
	cout<<"S= " << S << endl;
	RooAbsReal* B_fr= bkg.createIntegral(Bplus_MM,NormSet(Bplus_MM),Range("signal"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	cout<<"B= " << B << endl;
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
   	cout<<"S/B= " << S/B<< endl;

	cout << endl;
	cout << endl;
*/
	///Plot 
	//define strings to be used when saving plots
	
	string BsMassDistribution = "eps/Final/BmassFit";
	string BsMassResidual = "eps/Final/residual";
	string BsMassPull = "eps/Final/pull";
	string BsMassSweight = "eps/Final/Bs_sWeight";
        string simfit11 = "_sim11.eps";
        string simfit12 = "_sim12.eps";
        string StringDs2pipipi = "_Ds2pipipi";
        string StringDs2KKpi = "_Ds2KKpi";
        string eps = ".eps";


	string eleven = "_11.eps";
	string twelve = "_12.eps";

        string BsMassDistribution11 = "eps/Final/BmassFit";
        string BsMassDistribution12 = "eps/Final/BmassFit";
        string BsMassResidual11 = "eps/Final/residual";
        string BsMassResidual12 = "eps/Final/residual";
        string BsMassPull11 = "eps/Final/pull";
        string BsMassPull12 = "eps/Final/pull";
        string BsMassSweight11 = "eps/Final/Bs_sWeight";
        string BsMassSweight12 = "eps/Final/Bs_sWeight";

	// add Ds2pipipi tag to names where appropriate
	if(!Ds2KKpi){
		BsMassDistribution.append(StringDs2pipipi);
		BsMassResidual.append(StringDs2pipipi);
		BsMassPull.append(StringDs2pipipi);
		BsMassSweight.append(StringDs2pipipi);
	}

	// 2011 results
	if(sevenTeV && (!fitSimultan)){
		BsMassDistribution.append(eleven);
		BsMassResidual.append(eleven);
		BsMassPull.append(eleven);
		BsMassSweight.append(eleven);
	}

	//2012 results
	if((!sevenTeV) && (!fitSimultan)){
		BsMassDistribution.append(twelve);
		BsMassResidual.append(twelve);
		BsMassPull.append(twelve);
		BsMassSweight.append(twelve);
	}

	// results for simfit
        if(fitSimultan && Ds2KKpi  && (!combineYears)){
        BsMassDistribution11.append(simfit11);
        BsMassResidual11.append(simfit11);
        BsMassPull11.append(simfit11);
        BsMassSweight11.append(simfit11);

        BsMassDistribution12.append(simfit12);
        BsMassResidual12.append(simfit12);
        BsMassPull12.append(simfit12);
        BsMassSweight12.append(simfit12);
        }

        if(fitSimultan && (!Ds2KKpi)  && (!combineYears)){
        BsMassDistribution11.append(StringDs2pipipi);
        BsMassResidual11.append(StringDs2pipipi);
        BsMassPull11.append(StringDs2pipipi);
        BsMassSweight11.append(StringDs2pipipi);
        BsMassDistribution11.append(simfit11);
        BsMassResidual11.append(simfit11);
        BsMassPull11.append(simfit11);
        BsMassSweight11.append(simfit11);

        BsMassDistribution12.append(StringDs2pipipi);
        BsMassResidual12.append(StringDs2pipipi);
        BsMassPull12.append(StringDs2pipipi);
        BsMassSweight12.append(StringDs2pipipi);
        BsMassDistribution12.append(simfit12);
        BsMassResidual12.append(simfit12);
        BsMassPull12.append(simfit12);
        BsMassSweight12.append(simfit12);
        }

        if(fitSimultan && combineYears){
        BsMassDistribution12.append(StringDs2pipipi);
        BsMassResidual12.append(StringDs2pipipi);
        BsMassPull12.append(StringDs2pipipi);
        BsMassSweight12.append(StringDs2pipipi);
        BsMassDistribution12.append(eps);
        BsMassResidual12.append(eps);
        BsMassPull12.append(eps);
        BsMassSweight12.append(eps);

        BsMassDistribution11.append(StringDs2KKpi);
        BsMassResidual11.append(StringDs2KKpi);
        BsMassPull11.append(StringDs2KKpi);
        BsMassSweight11.append(StringDs2KKpi);
        BsMassDistribution11.append(eps);
        BsMassResidual11.append(eps);
        BsMassPull11.append(eps);
        BsMassSweight11.append(eps);

        }


	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= DTF_Bs_M.frame();
        RooPlot* frame_m_11= DTF_Bs_M.frame();
        RooPlot* frame_m_12= DTF_Bs_M.frame();
	frame_m->SetTitle("");
        frame_m_11->SetTitle("");
        frame_m_12->SetTitle("");

	TLatex* lhcbtext = new TLatex();
	lhcbtext->SetTextFont(132);
	lhcbtext->SetTextColor(1);
	lhcbtext->SetTextSize(0.05);
	lhcbtext->SetTextAlign(12);


	if(!fitSimultan){
		data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
		pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
		pdf->plotOn(frame_m,Components(DoubleGaussB0),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(DoubleGaussBs),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(RooArgSet(DstarKpipi_as_DsKpipi)),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(Dspipipi_as_DsKpipi),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(Dstarpipipi_as_DsKpipi),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(bkg_exp),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
		frame_m->Draw();
		if(makePlot == 1) c1->Print(BsMassDistribution.c_str());
	}

        if(fitSimultan){
                combData->plotOn(frame_m_11,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(90));
                simPdf->plotOn(frame_m_11,Name("pdf11"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussBs/*DoubleCBBs*/),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_11->Draw();
                if(makePlot == 1)  c1->Print(BsMassDistribution11.c_str());

                combData->plotOn(frame_m_12,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(90));
                simPdf->plotOn(frame_m_12,Name("pdf12"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussBs/*DoubleCBBs*/),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_12->Draw();
		//lhcbtext->DrawLatex(5480.,570.,"LHCb data");
		//lhcbtext->DrawLatex(5480.,500.,"L_{int} = 2 fb^{-1}");
		//lhcbtext->DrawLatex(5480.,430.,"(my analysis)");
                if(makePlot == 1) c1->Print(BsMassDistribution12.c_str());
        }

	///fit results
	double chi2 = 0;
	if(!fitSimultan)chi2 = frame_m->chiSquare("pdf","data",9);
        if(fitSimultan) chi2 = (frame_m_11->chiSquare("pdf11","data11",9) + frame_m_12->chiSquare("pdf12","data12",9)) / 2;
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
        RooHist* hresid;
        RooHist* hresid11;
        RooHist* hresid12;
        RooHist* hpull;
        RooHist* hpull11;
        RooHist* hpull12;

        if(!fitSimultan){
                hresid = frame_m->residHist("data","pdf") ;
                hpull = frame_m->pullHist("data","pdf") ;
        }
        if(fitSimultan){
                hresid11 = frame_m_11->residHist("data11","pdf11") ;
                hpull11 = frame_m_11->pullHist("data11","pdf11") ;
                hresid12 = frame_m_12->residHist("data12","pdf12") ;
                hpull12 = frame_m_12->pullHist("data12","pdf12") ;
        }

	// Create a new frame to draw the residual distribution and add the distribution to the frame
        RooPlot* frame2 = DTF_Bs_M.frame(Title("Residual Distribution")) ;
        RooPlot* frame2_11 = DTF_Bs_M.frame(Title("Residual Distribution 2011"));
        RooPlot* frame2_12 = DTF_Bs_M.frame(Title("Residual Distribution 2012"));

        if(!fitSimultan){
                frame2->SetTitle("");
                frame2->addPlotable(hresid,"P") ;
                frame2->Draw();
                if(makePlot == 1) c1->Print(BsMassResidual.c_str());
        }
        if(fitSimultan){
                frame2_11->SetTitle("");
                frame2_11->addPlotable(hresid11,"P") ;
                frame2_11->Draw();
                if(makePlot == 1) c1->Print(BsMassResidual11.c_str());

                frame2_12->SetTitle("");
                frame2_12->addPlotable(hresid12,"P") ;
                frame2_12->Draw();
                if(makePlot == 1) c1->Print(BsMassResidual12.c_str());
        }

	// Create a new frame to draw the pull distribution and add the distribution to the frame
        RooPlot* frame3 = DTF_Bs_M.frame(Title("Pull Distribution")) ;
        RooPlot* frame3_11 = DTF_Bs_M.frame(Title("Pull Distribution 2011")) ;
        RooPlot* frame3_12 = DTF_Bs_M.frame(Title("Pull Distribution 2012")) ;

        if(!fitSimultan){
                frame3->SetTitle("");
                frame3->SetLabelFont(62,"Y");
                frame3->addPlotable(hpull,"P") ;
                frame3->Draw();
                if(makePlot == 1) c1->Print(BsMassPull.c_str());
        }
        if(fitSimultan){
                frame3_11->SetTitle("");
                frame3_11->SetLabelFont(62,"Y");
                frame3_11->addPlotable(hpull11,"P") ;
                frame3_11->Draw();
                if(makePlot == 1) c1->Print(BsMassPull11.c_str());

                frame3_12->SetTitle("");
                frame3_12->SetLabelFont(62,"Y");
                frame3_12->addPlotable(hpull12,"P") ;
                frame3_12->Draw();
                if(makePlot == 1) c1->Print(BsMassPull12.c_str());
        }


	if(sWeight){

		meanB01.setConstant();
		meanBs1.setConstant();
		sigmaB2.setConstant();
		sigmaBs2.setConstant();
		exp_par.setConstant();
		meanDspipipi1.setConstant();
		meanDspipipi2.setConstant();


		///cloning DstarKpipiShifted into new pdf
		RooRealVar N_DstarKpipiShifted_11_var("N_DstarKpipiShifted_11_var", "N_DstarKpipiShifted_11_var", N_DstarKpipiShifted_11.getVal());
		RooRealVar N_DstarKpipiShifted_12_var("N_DstarKpipiShifted_12_var", "N_DstarKpipiShifted_12_var", N_DstarKpipiShifted_12.getVal());

        	RooAbsPdf* pdf11_clone;
        	RooAbsPdf* pdf12_clone;

		pdf11_clone = new RooAddPdf("pdf11_clone", "pdf11_clone", RooArgList(DoubleGaussB0, DoubleCBBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_11, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11, *N_Dstarpipipi_11, N_comb_11, N_DstarKpipiShifted_11_var));

		pdf12_clone = new RooAddPdf("pdf12_clone", "pdf12_clone", RooArgList(DoubleGaussB0, DoubleCBBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_12, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12_var));



		SPlot* sData11 = new SPlot("sData11","An SPlot for 2011",*data11,pdf11_clone,RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11, *N_Dstarpipipi_11, N_comb_11, N_DstarKpipiShifted_11_var)); 
		SPlot* sData12 = new SPlot("sData12","An SPlot for 2012",*data12,pdf12_clone,RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12_var));
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs11 = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist11 = (TH2*)data11->createHistogram("DTF_Bs_M,N_Bs_11_sw");
		SwDsHist11->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist11->SetTitle("");
		SwDsHist11->Draw();
		if(makePlot == 1) SwDs11->Print(BsMassSweight11.c_str());

		TCanvas* SwDs12 = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist12 = (TH2*)data12->createHistogram("DTF_Bs_M,N_Bs_12_sw");
		SwDsHist12->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist12->SetTitle("");
		SwDsHist12->Draw();
		if(makePlot == 1) SwDs12->Print(BsMassSweight12.c_str());


    		///Create output file
   		 TFile* output_11;
   		 TFile* output_12;
		 if(Ds2KKpi){
			output_11 = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs2DsKpipi_11_final_sweight.root","RECREATE");
			output_12 = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs2DsKpipi_12_final_sweight.root","RECREATE");
		}
		 if(!Ds2KKpi){
			output_11 = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs2DsKpipi_Ds2pipipi_11_final_sweight.root","RECREATE");
			output_12 = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs2DsKpipi_Ds2pipipi_12_final_sweight.root","RECREATE");
		}

		 tree11->SetBranchStatus("*",1);
   		 TTree* new_tree11 = tree11->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
    		 double w_11;
    		 TBranch* Bra_sw11 = new_tree11->Branch("N_Bs_sw", &w_11, "N_Bs_sw/D");

		 tree12->SetBranchStatus("*",1);
   		 TTree* new_tree12 = tree12->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
    		 double w_12;
    		 TBranch* Bra_sw12 = new_tree12->Branch("N_Bs_sw", &w_12, "N_Bs_sw/D");

  		  ///loop over events
    		  int numEvents_11 = new_tree11->GetEntries();
    		  int numEvents_12 = new_tree12->GetEntries();

    		  for(int i=0; i< numEvents_11; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_11 << endl;
			tree11->GetEntry(i);
			w_11=sData11->GetSWeight(i,"N_Bs_11_sw");
			Bra_sw11->Fill();
  		  }
	cout << "loop finished!!!" << endl;
   		 new_tree11->Write();
   		 output_11->Close();

    		  for(int i=0; i< numEvents_12; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_12 << endl;
			tree12->GetEntry(i);
			w_12=sData12->GetSWeight(i,"N_Bs_12_sw");
			Bra_sw12->Fill();
  		  }
	cout << "loop finished!!!" << endl;
   		 new_tree12->Write();
   		 output_12->Close();
	}

}

int main(int argc, char** argv){

    time_t startTime = time(0);
    //gROOT->ProcessLine(".x ../lhcbStyle.C");

    NamedParameter<string> Channel("Channel", (std::string) "");
    string Normalization = "Normalization";
    string Signal = "Signal";
    string Channel_convert = Channel;
    double *DsstarKpipifromNorm;


    if(Channel_convert.compare(Normalization) == 0) DsstarKpipifromNorm = fitBDTNorm();
    if(Channel_convert.compare(Signal) == 0) fitBDT();

    if((Channel_convert.compare(Normalization) != 0) && (Channel_convert.compare(Signal) != 0)){

	cout << "*********************************************************" << endl;
	cout << "*********************************************************" << endl;
    	cout << "please specify 'Signal' or 'Normalization' in options file" << endl;
	cout << "*********************************************************" << endl;
	cout << "*********************************************************" << endl;

    }
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
