#include <cmath>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include "TPaveText.h"
#include <sstream>
#include <RooDataSet.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
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
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

using namespace std;
using namespace RooFit ;
using namespace RooStats;

void plotEventVars(string Branch,string TitleX, int bins, double min, double max, bool useWeights=false,  bool reweight = false ){
    if(useWeights) reweight = false; 

   ///Load files
   TChain* tree=new TChain("DecayTree");
   tree->Add("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs2Dspipipi_11_final_sweight.root");
  // tree->Add("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_afterPreSel_sweight.root");	
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus(Branch.c_str(),1);
   tree->SetBranchStatus("N_Bs_sw",1);
   float var;
   double sw;
   tree->SetBranchAddress(Branch.c_str(),&var);
   tree->SetBranchAddress("N_Bs_sw",&sw);

   TString fileNameMC;
   if(useWeights)fileNameMC="/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_DsK3fb_Selection.root";
   else fileNameMC= "/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_forBDT_Reco14.root";
  // else fileNameMC= "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root";



   TFile* fileMC= new TFile(fileNameMC);
   TTree* treeMC = (TTree*) fileMC->Get("DecayTree");	
   treeMC->SetBranchStatus("*",0);
   treeMC->SetBranchStatus(Branch.c_str(),1);
   if(useWeights)treeMC->SetBranchStatus("weight",1);
   float varMC;
   double w;
   treeMC->SetBranchAddress(Branch.c_str(),&varMC);
   if(useWeights)treeMC->SetBranchAddress("weight",&w);
   else w=1.;

   ///Make histograms
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TString title= ";"+TitleX+";Yield [norm.]";
   TH1D* h= new TH1D(Branch.c_str(),title,bins,min,max);
   TH1D* h_MC= new TH1D((Branch+"_MC").c_str(),title,bins,min,max);
   TH1D* h_MC_rw= new TH1D((Branch+"_MC_rw").c_str(),title,bins,min,max);

    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);
	h->Fill(var,sw);
    }

    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
	treeMC->GetEntry(i);
	h_MC->Fill(varMC);
	h_MC_rw->Fill(varMC,w);
    }

    ///Plot it
    TCanvas* c= new TCanvas();

    h->Scale(1./h->Integral());
    h_MC->Scale(1./h_MC->Integral());
    h_MC_rw->Scale(1./h_MC_rw->Integral());
    double maxY= h->GetMaximum();
    if(h_MC->GetMaximum()>maxY)maxY=h_MC->GetMaximum();
    h->SetMinimum(0.);
    h->SetMaximum(maxY*1.2);
    h->SetLineColor(kBlack);
    h->Draw("");
    h_MC->SetLineColor(kRed);
    h_MC->Draw("same");
    h_MC_rw->SetLineColor(kBlue);
    if(useWeights)h_MC_rw->Draw("histsame");

    double KolmoTest = h->KolmogorovTest(h_MC);
    TPaveText *KolmOut= new TPaveText(0.6,0.55,0.85,0.7,"NDC");
    KolmOut->AddText(Form("Kolmogorov Test : %2f ", KolmoTest));
    KolmOut->SetLineColor(kWhite);
    KolmOut->SetFillColor(kWhite);
    KolmOut->SetShadowColor(0);
    KolmOut->SetTextSize(0.03);

    double KolmoTest_rw = h->KolmogorovTest(h_MC_rw);
    TPaveText *KolmOut_rw= new TPaveText(0.55,0.575,0.85,0.625,"NDC");
    KolmOut_rw->AddText(Form("Kolmogorov Test : %2f ", KolmoTest_rw));
    KolmOut_rw->SetLineColor(kWhite);
    KolmOut_rw->SetTextColor(kBlue);
    KolmOut_rw->SetFillColor(kWhite);
    KolmOut_rw->SetShadowColor(0);
    KolmOut_rw->SetTextSize(0.03);

    TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    leg->SetHeader(" ");
    leg->AddEntry(h,"Data","LEP");
    leg->AddEntry(h_MC,"MC","LEP");
    leg->SetLineColor(kWhite);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.05);
    KolmOut->Draw();
    leg->Draw(); 
    if(useWeights)KolmOut_rw->Draw();

    if(useWeights)c->Print(("DataVsReweightedMC/"+Branch+".eps").c_str());
    else c->Print(("DataVsMC/"+Branch+".eps").c_str());

    ///calculate weights
    if(reweight){
	TFile* output=new TFile((Branch+"_weights_norm.root").c_str(),"RECREATE");
	TH1D *h_weight = (TH1D*)h->Clone();
	h_weight->SetName((Branch+"_weight").c_str());
	h_weight->Divide(h,h_MC);
	h_weight->SetMinimum(0.);
  	h_weight->SetMaximum(h_weight->GetMaximum()*1.2);
	h_weight->GetYaxis()->SetTitle("Weight");
	h_weight->Draw();
	c->Print(("DataVsMC/"+Branch+"_weights.eps").c_str());
	h_weight->Write();
	output->Write();
	output->Close();
    }
}

void plot(string Branch,string TitleX, int bins, double min, double max, bool useWeights=false, bool reweight = false){
   if(useWeights) reweight = false; 

   ///Load files
   TChain* tree=new TChain("DecayTree");
   //tree->Add("/auto/data/dargent/Bs2DsKpipi/final/data2012_Ds2KKpi_sweight.root");
   tree->Add("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_afterPreSel_sweight.root");
   //tree->Add("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root");

   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus(Branch.c_str(),1);
   tree->SetBranchStatus("N_Bs_sw",1);
   int var;
   double sw;
   tree->SetBranchAddress(Branch.c_str(),&var);
   tree->SetBranchAddress("N_Bs_sw",&sw);

   TString fileNameMC;
   if(useWeights)fileNameMC="/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_DsK3fb_Selection.root";
   //else fileNameMC= "/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_forBDT_Reco14.root";
   else fileNameMC= "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root";

   TChain* treeMC =new TChain("DecayTree");
   treeMC->Add(fileNameMC);
   treeMC->SetBranchStatus("*",0);
   treeMC->SetBranchStatus(Branch.c_str(),1);
  // treeMC->SetBranchStatus("Bplus_BKGCAT",1);
   if(useWeights)treeMC->SetBranchStatus("weight",1);
   int varMC;
   double w;
  // int cat;
   treeMC->SetBranchAddress(Branch.c_str(),&varMC);
   //treeMC->SetBranchAddress("Bplus_BKGCAT",&cat);
   if(useWeights)treeMC->SetBranchAddress("weight",&w);
   else w=1.;

   ///Make histograms
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TString title= ";"+TitleX+";Yield [norm.]";
   TH1D* h= new TH1D(Branch.c_str(),title,bins,min,max);
   TH1D* h_MC= new TH1D((Branch+"_MC").c_str(),title,bins,min,max);
   TH1D* h_MC_rw= new TH1D((Branch+"_MC_rw").c_str(),title,bins,min,max);

    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);
	h->Fill(var,sw);
    }

    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
	treeMC->GetEntry(i);
        //if(cat>10)continue;
	h_MC->Fill(varMC);
	h_MC_rw->Fill(varMC,w);
    }

    ///Plot it
    TCanvas* c= new TCanvas();

    h->Scale(1./h->Integral());
    h_MC->Scale(1./h_MC->Integral());
    h_MC_rw->Scale(1./h_MC_rw->Integral());
    double maxY= h->GetMaximum();
    if(h_MC->GetMaximum()>maxY)maxY=h_MC->GetMaximum();
    h->SetMinimum(0.);
    h->SetMaximum(maxY*1.2);
    h->SetLineColor(kBlack);
    h->Draw("");
    h_MC->SetLineColor(kRed);
    h_MC->Draw("same");
    h_MC_rw->SetLineColor(kBlue);
    if(useWeights)h_MC_rw->Draw("histsame");

    double KolmoTest = h->KolmogorovTest(h_MC);
    TPaveText *KolmOut= new TPaveText(0.55,0.6,0.85,0.7,"NDC");
    KolmOut->AddText(Form("Kolmogorov Test : %2f ", KolmoTest));
    KolmOut->SetLineColor(kWhite);
    KolmOut->SetFillColor(kWhite);
    KolmOut->SetShadowColor(0);
    KolmOut->SetTextSize(0.03);
    KolmOut->SetTextColor(kRed);

    double KolmoTest_rw = h->KolmogorovTest(h_MC_rw);
    TPaveText *KolmOut_rw= new TPaveText(0.55,0.575,0.85,0.625,"NDC");
    KolmOut_rw->AddText(Form("Kolmogorov Test : %2f ", KolmoTest_rw));
    KolmOut_rw->SetLineColor(kWhite);
    KolmOut_rw->SetTextColor(kBlue);
    KolmOut_rw->SetFillColor(kWhite);
    KolmOut_rw->SetShadowColor(0);
    KolmOut_rw->SetTextSize(0.03);

    TLegend *leg = new TLegend(0.55,0.7,0.85,0.85);
    leg->SetHeader(" ");
    leg->AddEntry(h,"Data","LEP");
    leg->AddEntry(h_MC,"MC","LEP");
    if(useWeights)leg->AddEntry(h_MC_rw,"MC (reweighted)","LEP");
    leg->SetLineColor(kWhite);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.05);
    KolmOut->Draw();
    if(useWeights)KolmOut_rw->Draw();
    leg->Draw(); 

    if(useWeights)c->Print(("DataVsReweightedMC/"+Branch+".eps").c_str());
    else c->Print(("DataVsMC/"+Branch+".eps").c_str());

    ///calculate weights
    if(reweight){
	TFile* output=new TFile((Branch+"_weights_B2D3Pi12.root").c_str(),"RECREATE");
	TH1D *h_weight = (TH1D*)h->Clone();
	h_weight->SetName((Branch+"_weight").c_str());
	h_weight->Divide(h,h_MC);
	h_weight->SetMinimum(0.);
  	h_weight->SetMaximum(h_weight->GetMaximum()*1.2);
	h_weight->GetYaxis()->SetTitle("Weight");
	h_weight->Draw();
	c->Print(("DataVsMC/"+Branch+"_weights.eps").c_str());
	h_weight->Write();
	output->Write();
	output->Close();
    }
}

void reweight(){
	///Get weight histos
	//TFile* vtxFile=new TFile("");	
	//TH1D* h_vtx=(TH1D*)vtxFile->Get("Bplus_ENDVERTEX_CHI2_weight");
	
	//TFile* etaFile=new TFile("Bs_ETA_weights.root");	
	//TH1D* h_eta=(TH1D*)etaFile->Get("Bs_ETA_weight");

	TFile* ghostFile=new TFile("max_ghostProb_weights_B2D3Pi12.root");	
	TH1D* h_ghost=(TH1D*)ghostFile->Get("max_ghostProb_weight");
	
	TFile* nTFile=new TFile("nTracks_weights_B2D3Pi12.root");	
	TH1D* h_nT=(TH1D*)nTFile->Get("nTracks_weight");

	//TFile* isMuFile=new TFile("");	
	//TH2D* h_isMu=(TH2D*)isMuFile->Get("hist");

	//TFile* trackFile=new TFile("");	
	//TH2D* h_track=(TH2D*)trackFile->Get("Ratio");

        /*
	///Draw tables
	TCanvas* c=new TCanvas();
	h_track->UseCurrentStyle();
	h_track->SetTitle("");
	h_track->GetXaxis()->SetTitle("p [GeV]");
	h_track->GetYaxis()->SetTitleOffset(1.);
	h_track->Draw("colz");
	c->Print("DataVsMC/trackEff.eps");
	h_isMu->UseCurrentStyle();
	h_isMu->SetTitle("");
	h_isMu->GetXaxis()->SetTitle("p_{T} [MeV]");
	h_isMu->GetYaxis()->SetTitle("p [MeV]");
	h_isMu->Draw("colz");
	c->Print("DataVsMC/isMuEff.eps");
        */
	///Load MC file
	TFile* fileMC= new TFile("/auto/data/kecke/B2DPiPiPi/forMaster/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root");
	//TFile* fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_forBDT_Reco14.root");
   	TTree* treeMC = (TTree*) fileMC->Get("DecayTree");
	float ghostProb;
	int nT;
	//double K_p,K_eta;
	//double pip_p,pip_eta;
	//double pim_p,pim_eta;
	//double mup_p,mup_eta,mup_pt;
	//double mum_p,mum_eta,mum_pt;

   	//treeMC->SetBranchAddress("Bs_ENDVERTEX_CHI2",&chi);
   	//treeMC->SetBranchAddress("Bs_P",&p);
   	treeMC->SetBranchAddress("max_ghostProb",&ghostProb);
   	treeMC->SetBranchAddress("nTracks",&nT);
   	//treeMC->SetBranchAddress("Kplus_P",&K_p);
   	//treeMC->SetBranchAddress("Kplus_ETA",&K_eta);
   	//treeMC->SetBranchAddress("piplus_P",&pip_p);
   	//treeMC->SetBranchAddress("piplus_ETA",&pip_eta);
   	//treeMC->SetBranchAddress("piminus_P",&pim_p);
   	//treeMC->SetBranchAddress("piminus_ETA",&pim_eta);
	///Create new tree
	double w;
    	//TFile* output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_DsK3fb_Selection.root","RECREATE");
	TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/forMaster/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root","RECREATE");
    	TTree* new_tree = treeMC->CloneTree();//CopyTree();    
    	TBranch* Bra_w = new_tree->Branch("weight",&w,"weight/D");

        treeMC->SetBranchStatus("*",0);
   	treeMC->SetBranchStatus("*ghostProb*",1);
   	treeMC->SetBranchStatus("nTracks",1);

   	///loop over MC events
   	int numEventsMC = treeMC->GetEntries();
   	for(int i=0; i< numEventsMC; i++)
  	{	
		if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
		treeMC->GetEntry(i);
		w=1.;
		double tmp=w;
		/*
                tmp=h_vtx->GetBinContent(h_vtx->FindBin(chi));
		if(tmp==0)tmp=1.;
		w*=tmp;
		*/
                tmp=h_ghost->GetBinContent(h_ghost->FindBin(ghostProb));
		if(tmp==0)tmp=1.;
		w*=tmp;
		/*
                tmp=h_eta->GetBinContent(h_eta->FindBin(eta));
		if(tmp==0)tmp=1.;
		//w*=tmp;
		*/
		
		tmp=h_nT->GetBinContent(h_nT->FindBin(nT));
		if(tmp==0)tmp=1.;
		w*=tmp;
		/*
		tmp=h_track->GetBinContent(h_track->FindBin(K_p/1000.,K_eta));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_track->GetBinContent(h_track->FindBin(pip_p/1000.,pip_eta));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_track->GetBinContent(h_track->FindBin(pim_p/1000.,pim_eta));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_track->GetBinContent(h_track->FindBin(mup_p/1000.,mup_eta));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_track->GetBinContent(h_track->FindBin(mum_p/1000.,mum_eta));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_isMu->GetBinContent(h_isMu->FindBin(mup_pt,mup_p));
		if(tmp==0)tmp=1.;
		w*=tmp;
		tmp=h_isMu->GetBinContent(h_isMu->FindBin(mum_pt,mum_p));
		if(tmp==0)tmp=1.;
		w*=tmp;
		*/
		Bra_w->Fill();
    	}

	new_tree->Write();

	//etaFile->Close();
	ghostFile->Close();
	nTFile->Close();
	//isMuFile->Close();
	//trackFile->Close();
	fileMC->Close();
	output->Close();
}

void dataVsMC(bool useWeights=false){
	///B
	plot("Bs_P","p(B) [MeV]",40,0,500000,useWeights);
	plot("Bs_PT","p_{T}(B) [MeV]",40,0,40000,useWeights);
	plot("Bs_ETA","#eta(B)",40,1,6,useWeights);
	plot("Bs_IPCHI2_OWNPV","#chi^{2}_{IP}(B)",40,0,16,useWeights);
	plot("Bs_FDCHI2_OWNPV","#chi^{2}_{FD}(B)",40,0,100000,useWeights);
	plot("Bs_ENDVERTEX_CHI2","#chi^{2}_{vtx}(B)",40,0,35,useWeights);
	plot("Bs_TAU","#tau(B) [ns]",40,0,0.012,useWeights);
	plot("Bs_DIRA_OWNPV","cos(DIRA)",40,0.9999,1,useWeights);
	
	///K
	plot("K_plus_P","p(K^{+}) [MeV]",40,0,180000,useWeights);
	plot("K_plus_PT","p_{T}(K^{+}) [MeV]",40,0,10000,useWeights);
	plot("K_plus_ETA","#eta(K^{+})",40,1,6,useWeights);
	plot("K_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",40,0,13000,useWeights);
	//plot("Kplus_PIDK","DLL_{K#pi}(K^{+}) ",100,-100,100);
	plot("angK","#theta_{D_{s} K^{+}}",40,0,3.141,useWeights);
	plot("K_plus_ptasy_1.00","pt cone asymmetry (K^{+})",40,-1,1,useWeights);
	plot("K_plus_TRACK_GhostProb","ghost prob (K^{+})",40,0,0.4,useWeights);
	///pi+
	plot("pi_plus_P","p(#pi^{+}) [MeV]",40,0,180000,useWeights);
	plot("pi_plus_PT","p_{T}(#pi^{+}) [MeV]",40,0,10000,useWeights);
	plot("pi_plus_ETA","#eta(#pi^{+})",40,1,6,useWeights);
	plot("pi_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",40,0,13000,useWeights);
	//plot("piplus_PIDK","DLL_{K#pi}(#pi^{+}) ",100,-100,100);
	plot("angPip","#theta_{D_{s} #pi^{+}}",40,0,3.141,useWeights);
	plot("pi_plus_ptasy_1.00","pt cone asymmetry (#pi^{+})",40,0,1,useWeights);
	plot("pi_plus_TRACK_GhostProb","ghost prob (#pi^{+})",40,0,0.4,useWeights);
	///pi-
	plot("pi_minus_P","p(#pi^{-}) [MeV]",40,0,180000,useWeights);
	plot("pi_minus_PT","p_{T}(#pi^{-}) [MeV]",40,0,10000,useWeights);
	plot("pi_minus_ETA","#eta(#pi^{-})",40,1,6,useWeights);
	plot("pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",40,0,13000,useWeights);
	//plot("piminus_PIDK","DLL_{K#pi}(#pi^{-}) ",100,-100,100);
	plot("angPim","#theta_{D_{s} #pi^{-}}",40,0,3.141,useWeights);
	plot("pi_minus_ptasy_1.00","pt cone asymmetry (#pi^{-})",40,-1,1,useWeights);
	plot("pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",40,0,0.4,useWeights);

	///Xs
        plot("Xs_max_DOCA","X_{s} max DOCA [mm]",40, 0, 0.4 ,useWeights);
        plot("XsDaughters_min_IPCHI2","X_{s} min(#chi^{2}_{IP})",40, 0, 10 ,useWeights);

	///Ds
	plot("Ds_P","p(D_{s}) [MeV]",40,0,400000,useWeights);
	plot("Ds_PT","p_{T}(D_{s}) [MeV]",40,0,40000,useWeights);
	plot("Ds_ETA","#eta(D_{s})",40,1,6,useWeights);
	plot("Ds_FDCHI2_ORIVX","#chi^{2}_{FD}(D_{s})",40,0,40000,useWeights);
	plot("Ds_DIRA_OWNPV","cos(DIRA) (D_{s})",40,0.9999,1,useWeights);
        plot("DsDaughters_min_IPCHI2","D_{s} min(#chi^{2}_{IP})",40, 0, 10 ,useWeights);

	///K+ from Ds
	plot("K_plus_fromDs_P","p(K^{+} from D_{s}) [MeV]",40,0,180000,useWeights);
	plot("K_plus_fromDs_PT","p_{T}(K^{+} from D_{s}) [MeV]",40,0,10000,useWeights);
	plot("K_plus_fromDs_ETA","#eta(K^{+} from D_{s})",40,1,6,useWeights);
	plot("K_plus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+} from D_{s})",40,0,13000,useWeights);
	plot("K_plus_fromDs_ptasy_1.00","pt cone asymmetry (K^{+} from D_{s})",40,-1,1,useWeights);
	plot("K_plus_fromDs_TRACK_GhostProb","ghost prob (K^{+} from D_{s})",40,0,0.4,useWeights);
	///pi- from Ds
	plot("pi_minus_fromDs_P","p(#pi^{-} from D_{s}) [MeV]",40,0,180000,useWeights);
	plot("pi_minus_fromDs_PT","p_{T}(#pi^{-} from D_{s}) [MeV]",40,0,10000,useWeights);
	plot("pi_minus_fromDs_ETA","#eta(#pi^{-} from D_{s})",40,1,6,useWeights);
	plot("pi_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-} from D_{s})",40,0,13000,useWeights);
	plot("pi_minus_fromDs_ptasy_1.00","pt cone asymmetry (#pi^{-} from D_{s})",40,-1,1,useWeights);
	plot("pi_minus_fromDs_TRACK_GhostProb","ghost prob (#pi^{-} from D_{s})",40,0,0.4,useWeights);
	///K- from Ds
	plot("K_minus_fromDs_P","p(K^{-} from D_{s}) [MeV]",40,0,180000,useWeights);
	plot("K_minus_fromDs_PT","p_{T}(K^{-} from D_{s}) [MeV]",40,0,10000,useWeights);
	plot("K_minus_fromDs_ETA","#eta(K^{-} from D_{s})",40,1,6,useWeights);
	plot("K_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-} from D_{s})",40,0,13000,useWeights);
	plot("K_minus_fromDs_ptasy_1.00","pt cone asymmetry (K^{-} from D_{s})",40,-1,1,useWeights);
	plot("K_minus_fromDs_TRACK_GhostProb","ghost prob (K^{-} from D_{s})",40,0,0.4,useWeights);
	///event
	//plotEventVars("nPV","N_{PV}",10,0,10,useWeights);
      //  plotEventVars("max_ghostProb","max(Track_ghostProb)",25,0.,0.375,useWeights);
      //  plotEventVars("nTracks","# of tracks",40,0,450,useWeights);
}

void dataVsReweightedMC(){
   //reweight();
   dataVsMC(true);
}
 
void plotPID(string Branch,string TitleX, int bins, double min, double max, bool useWeights=true) {
    ///Set plot options
   gStyle->SetOptStat(0);
   gStyle->SetTitleXSize(0.05);
   gStyle->SetTitleYSize(0.05);
   gStyle->SetTitleXOffset(0.8);
   //gStyle->SetTitleYOffset(1.3);
   gStyle->SetTitleFont(42,"X");
   gStyle->SetTitleFont(42,"Y");
   gStyle->SetLabelFont(42,"X");
   gStyle->SetLabelFont(42,"Y");
   gStyle->SetLabelOffset(0.005,"X");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   
   ///Load files
   TFile* file= new TFile("/auto/data/dargent/Bu2JpsiKpipi/data/data_bdt_PIDK.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus(Branch.c_str(),1);
   tree->SetBranchStatus("n_sig_sw",1);
   double pid_data,sw;
   tree->SetBranchAddress(Branch.c_str(),&pid_data);
   tree->SetBranchAddress("n_sig_sw",&sw);

   TString fileNameMCcorr="/auto/data/dargent/Bu2JpsiKpipi/MC/MC_reweighted_forBDT_PIDcorr_cat10.root";
   TFile* fileMCcorr= new TFile(fileNameMCcorr);
   TTree* treeMCcorr = (TTree*) fileMCcorr->Get("DecayTree");	
   treeMCcorr->SetBranchStatus("*",0);
   treeMCcorr->SetBranchStatus(Branch.c_str(),1);
   treeMCcorr->SetBranchStatus((Branch+"_corr").c_str(),1);
   treeMCcorr->SetBranchStatus("weight",1);
   double pid_MC_corr, w;
   treeMCcorr->SetBranchAddress((Branch+"_corr").c_str(),&pid_MC_corr);
   treeMCcorr->SetBranchAddress("weight",&w);

   TString fileNameMC="/auto/data/dargent/Bu2JpsiKpipi/MC/MC_reweighted_forBDT.root";
   TFile* fileMC= new TFile(fileNameMC);
   TTree* treeMC = (TTree*) fileMC->Get("DecayTree");	
   treeMC->SetBranchStatus("*",0);
   treeMC->SetBranchStatus(Branch.c_str(),1);
   treeMC->SetBranchStatus("weight",1);
   double pid_MC;
   treeMC->SetBranchAddress(Branch.c_str(),&pid_MC);

   ///Make histograms
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TString title= ";"+TitleX+";Yield [norm.]";
   TH1D* h= new TH1D(Branch.c_str(),title,bins,min,max);
   TH1D* h_MC= new TH1D((Branch+"_MC").c_str(),title,bins,min,max);
   TH1D* h_MC_rw= new TH1D((Branch+"_MC_rw").c_str(),title,bins,min,max);

    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);
	h->Fill(pid_data,sw);
    }

    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
	treeMC->GetEntry(i);
	h_MC->Fill(pid_MC);
    }
    int numEventsMCcorr = treeMCcorr->GetEntries();
    
    for(int i=0; i< numEventsMCcorr; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMCcorr << endl;
	treeMCcorr->GetEntry(i);
	h_MC_rw->Fill(pid_MC_corr,w);
    }

    ///Plot it
    TCanvas* c= new TCanvas();

    h->Scale(1./h->Integral());
    h_MC->Scale(1./h_MC->Integral());
    h_MC_rw->Scale(1./h_MC_rw->Integral());
    double maxY= h->GetMaximum();
    if(h_MC->GetMaximum()>maxY)maxY=h_MC->GetMaximum();
    h->SetMinimum(0.);
    h->SetMaximum(maxY*1.2);
    h->SetLineColor(kBlack);
    h->Draw("");
    h_MC->SetLineColor(kRed);
    h_MC->Draw("histsame");
    h_MC_rw->SetLineColor(kBlue);
    h_MC_rw->Draw("histsame");
    c->Print(("DataVsReweightedMC/"+Branch+".eps").c_str());
}

void applyBDTcut(){
   ///Load file
   TFile* file= new TFile("/auto/data/dargent/Bu2JpsiKpipi/MC/MC_bdt_PIDcorr_cat10.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus("Bplus_MM",1);
   tree->SetBranchStatus("*PX",1);
   tree->SetBranchStatus("*PY",1);
   tree->SetBranchStatus("*PZ",1);
   tree->SetBranchStatus("*PE",1);
   tree->SetBranchStatus("m*pi",1);
   tree->SetBranchStatus("DTF_m*pi",1);
   tree->SetBranchStatus("BDT*",1);
   tree->SetBranchStatus("sample",1);
   tree->SetBranchStatus("*TRUE*",1);
   tree->SetBranchStatus("weight",1);

   TFile* output=new TFile("/auto/data/dargent/Bu2JpsiKpipi/MC/MC_final_PIDcorr_cat10.root","RECREATE");
   TTree* new_tree = tree->CopyTree("BDT_response2>0.1");
   new_tree->Write();
   
   file->Close();
   output->Close();	
}

void compareBDTresponse(){
   ///Load files
   TFile* file= new TFile("/auto/data/dargent/Bu2JpsiKpipi/data/data_bdt_PIDK.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus("BDT_response2",1);
   tree->SetBranchStatus("n_sig_sw",1);
   tree->SetBranchStatus("DTF_mKpipi",1);
   tree->SetBranchStatus("Bplus_MM",1);
   tree->SetBranchStatus("sample",1);

   double var,sw;
   double mKpipi,Bplus_MM;
   int sample;

   tree->SetBranchAddress("BDT_response2",&var);
   tree->SetBranchAddress("n_sig_sw",&sw);
   tree->SetBranchAddress("DTF_mKpipi",&mKpipi);
   tree->SetBranchAddress("Bplus_MM",&Bplus_MM);
   tree->SetBranchAddress("sample",&sample);

   TString fileNameMC;
   fileNameMC="/auto/data/dargent/Bu2JpsiKpipi/MC/MC_bdt_PIDcorr_cat10.root";
   TFile* fileMC= new TFile(fileNameMC);
   TTree* treeMC = (TTree*) fileMC->Get("DecayTree");	
   treeMC->SetBranchStatus("*",0);
   treeMC->SetBranchStatus("BDT_response2",1);
   treeMC->SetBranchStatus("weight",1);
   treeMC->SetBranchStatus("DTF_mKpipi",1);
   treeMC->SetBranchStatus("Bplus_MM",1);
   double varMC, w, mKpipiMC, mB_mc;
   treeMC->SetBranchAddress("BDT_response2",&varMC);
   treeMC->SetBranchAddress("weight",&w);
   treeMC->SetBranchAddress("DTF_mKpipi",&mKpipiMC);
   treeMC->SetBranchAddress("Bplus_MM",&mB_mc);

   ///Make histograms
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TString title= "; BDT ; Yield [norm.]";
   TH1D* h= new TH1D("BDT",title,100,-1,1);
   TH1D* h_MC= new TH1D("BDT_MC",title,100,-1,1);
   TH1D* h_MC_rw= new TH1D("BDT_MC_rw",title,100,-1,1);
   
   TProfile* prof_S = new TProfile("prof_S",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; <BDT>",30,0.6,4.8,-0.5,1.);
   TProfile* prof_B = new TProfile("prof_B",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; <BDT>",30,0.6,4.8,-0.5,1.);
   TProfile* prof_S_mB = new TProfile("prof_S_mB",";m(J/#psi K^{+} #pi^{+} #pi^{-}) [MeV]; <BDT>",20,5200,5600,-0.5,1.);
   TProfile* prof_B_mB = new TProfile("prof_S_mB",";m(J/#psi K^{+} #pi^{+} #pi^{-}) [MeV]; <BDT>",20,5200,5600,-0.5,1.);

    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);
	h->Fill(var,sw);
	if(sample==1 && Bplus_MM>5350.) {
		prof_B->Fill(mKpipi/1000000.,var);
		prof_B_mB->Fill(Bplus_MM,var);
	}
	//prof_S->Fill(mKpipi/1000000.,var,sw);
    }

    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
	treeMC->GetEntry(i);
	h_MC->Fill(varMC);
	h_MC_rw->Fill(varMC,w);
	prof_S->Fill(mKpipiMC/1000000.,varMC,w);
	prof_S_mB->Fill(mB_mc,varMC,w);
    }

    ///Plot it
    TCanvas* c= new TCanvas();

    h->Scale(1./h->Integral());
    h_MC->Scale(1./h_MC->Integral());
    h_MC_rw->Scale(1./h_MC_rw->Integral());
    double maxY= h->GetMaximum();
    if(h_MC->GetMaximum()>maxY)maxY=h_MC->GetMaximum();
    h->SetMinimum(0.);
    h->SetMaximum(maxY*1.2);
    h->SetLineColor(kBlack);
    h->Draw("");
    h_MC->SetLineColor(kRed);
    h_MC->Draw("histsame");
    h_MC_rw->SetLineColor(kBlue);
    h_MC_rw->Draw("histsame");
    c->Print("DataVsReweightedMC/BDT.eps");

    prof_S->SetLineColor(kBlue);
    prof_B->SetLineColor(kRed);
    prof_S->SetMinimum(-0.4);
    prof_S->SetMaximum(0.6);
    prof_S->Draw();
    prof_B->Draw("SAME");
    c->Print("BDT_corr_prof.eps");

    prof_S_mB->SetLineColor(kBlue);
    prof_B_mB->SetLineColor(kRed);
    prof_S_mB->SetMinimum(-0.4);
    prof_S_mB->SetMaximum(0.6);
    prof_S_mB->Draw();
    prof_B_mB->Draw("SAME");
    c->Print("BDT_corr_prof_mB.eps");

}



void applyPIDweights(){
	///Get weight histos
	TFile* DLLKFile=new TFile("/home/he/dargent/cmtuser/Urania_v4r0/PIDCalib/PIDPerfScripts/scripts/python/MultiTrack/bkg_Ds3pi_K_PID_weights.root");	
   	TTree* treeDLLK = (TTree*) DLLKFile->Get("CalibTool_PIDCalibTree");
	float w_DLLK;
  	treeDLLK->SetBranchAddress("Event_PIDCalibEff",&w_DLLK);
 
	///Load MC file
	TFile* fileMC= new TFile("/auto/data/dargent/Bs2DsKpipi/preselection/bkg_Ds3pi.root");
   	TTree* treeMC = (TTree*) fileMC->Get("DecayTree");
	double mB;
  	treeMC->SetBranchAddress("Bs_MM",&mB);
 
  	///Create new tree
	double w;
    	TFile* output = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection/bkg_Ds3pi_PID.root","RECREATE");
    	TTree* new_tree = treeMC->CopyTree("3000 <= K_plus_P && 100000 > K_plus_P && 1.5 <= K_plus_ETA && 5 > K_plus_ETA && 0 <= nTracks && 500 > nTracks && 3000 <= pi_plus_P && 100000 > pi_plus_P && 1.5 <= pi_plus_ETA && 5 > pi_plus_ETA && 3000 <= pi_minus_P && 100000 > pi_minus_P && 1.5 <= pi_minus_ETA && 5 > pi_minus_ETA && 3000 <= K_plus_fromDs_P && 100000 > K_plus_fromDs_P && 1.5 <= K_plus_fromDs_ETA && 5 > K_plus_fromDs_ETA && 3000 <= K_minus_fromDs_P && 100000 > K_minus_fromDs_P && 1.5 <= K_minus_fromDs_ETA && 5 > K_minus_fromDs_ETA && 3000 <= pi_minus_fromDs_P && 100000 > pi_minus_fromDs_P && 1.5 <= pi_minus_fromDs_ETA && 5 > pi_minus_fromDs_ETA" );

    	TBranch* Bra_w = new_tree->Branch("PID_eff",&w,"PID_eff/D");

   	///loop over MC events
   	int numEventsMC = new_tree->GetEntries();
	if(numEventsMC != treeDLLK->GetEntries()) { 
		cout << "Event numbers incompatible" << endl;
		cout << numEventsMC << endl;
		cout << treeDLLK->GetEntries() << endl;
		return;
	}


	TH1D* h_mB = new TH1D("","",50,4800,6500);

   	for(int i=0; i< numEventsMC; i++)
  	{	
		if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
		new_tree->GetEntry(i);
		treeDLLK->GetEntry(i);
		w=w_DLLK;
		Bra_w->Fill();
		h_mB->Fill(mB,w);
    	}

	TCanvas* c = new TCanvas();
	h_mB->Draw();
	c->Print("mB.eps");

	new_tree->Write();

	fileMC->Close();
	DLLKFile->Close();
	output->Close();
}


int main(){
   time_t startTime = time(0);
   
   ///Set plot options
   gStyle->SetOptStat(0);
   gStyle->SetTitleXSize(0.05);
   gStyle->SetTitleYSize(0.05);
   gStyle->SetTitleXOffset(0.8);
   //gStyle->SetTitleYOffset(1.3);
   gStyle->SetTitleFont(42,"X");
   gStyle->SetTitleFont(42,"Y");
   gStyle->SetLabelFont(42,"X");
   gStyle->SetLabelFont(42,"Y");
   gStyle->SetLabelOffset(0.005,"X");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();

  // plotEventVars("K_minus_fromDs_PIDK","#Delta LLK of K^{+}",25,-4,4,true,false);
  // plotEventVars("K_minus_fromDs_PIDp","#Delta LLp of K^{+}",25,-4,4,true,false);

  // plot("nTracks","# of tracks",40,0,450,false, true);
  // plot("nPV","# of PV",10,0,10,false,true);
  // plot("max_ghostProb","max(Track_ghostProb)",25,0.,0.375,false,true);

  // plot("max_ghostProb","max(Track_ghostProb)",25,0.,0.375,true);
   // dataVsMC(); 

    reweight();
   // dataVsReweightedMC(); 

    //compareBDTresponse();
    //applyBDTcut();

    //applyPIDweights();
    
    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;

return 0;
}

