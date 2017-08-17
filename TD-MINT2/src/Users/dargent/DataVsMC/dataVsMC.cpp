// MC studies
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
#include <TLegend.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include "TProfile.h"
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"

using namespace std;
using namespace MINT;

void plot(string Branch,string TitleX, int bins, double min, double max, int Year = 11, TString finalState = "KKpi", bool useWeights=false, TString Decay = "norm"){
    
    ///Load files
    TChain* tree=new TChain("DecayTree");
    tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_sweight.root");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus(Branch.c_str(),1);
    tree->SetBranchStatus("N_Bs_sw",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("Ds_finalState",1);
    double var;
    double sw;
    int year,Ds_finalState;
    tree->SetBranchAddress(Branch.c_str(),&var);
    tree->SetBranchAddress("N_Bs_sw",&sw);
    tree->SetBranchAddress("year",&year);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
           
    TString fileNameMC;
    if(useWeights)fileNameMC="";

    TChain* treeMC =new TChain("DecayTree");
    treeMC->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11.root");
    treeMC->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12.root");
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus(Branch.c_str(),1);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("Bs_BKGCAT",1);
    if(useWeights)treeMC->SetBranchStatus("weight",1);
   
    double varMC;
    double w;
    int cat,yearMC,Ds_finalStateMC;
    treeMC->SetBranchAddress(Branch.c_str(),&varMC);
    treeMC->SetBranchAddress("Bs_BKGCAT",&cat);
    treeMC->SetBranchAddress("year",&yearMC);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalStateMC);           
    if(useWeights)treeMC->SetBranchAddress("weight",&w);
    else w=1.;
    
    ///Make histograms
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
	if(year != Year) continue;
	if(finalState == "KKpi" && Ds_finalState == 3) continue;
        h->Fill(var,sw);
    }
    
    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
        treeMC->GetEntry(i);
	if(yearMC != Year) continue;
	if(finalState == "KKpi" && Ds_finalStateMC == 3) continue;
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
    h_MC->SetMarkerColor(kRed);
    h_MC->SetLineColor(kRed);
    h_MC->Draw("esame");
    h_MC_rw->SetLineColor(kBlue);
    h_MC_rw->SetMarkerColor(kBlue);
    if(useWeights)h_MC_rw->Draw("esame");
    
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
    
    /*
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
    */
}

void dataVsMC(int Year = 11, TString finalState = "KKpi", bool useWeights=false, TString Decay = "norm"){
    
    ///B
    plot("Bs_P","p(B) [MeV]",40,0,900000,Year, finalState,useWeights);
    plot("Bs_PT","p_{T}(B) [MeV]",40,0,40000,Year, finalState,useWeights);
    plot("Bs_ETA","#eta(B)",40,1,6,Year, finalState,useWeights);
    plot("Bs_FDCHI2_OWNPV","#chi^{2}_{FD}(B)",40,0,100000,Year, finalState,useWeights);
    plot("Bs_ENDVERTEX_CHI2","#chi^{2}_{vtx}(B)",40,0,35,Year, finalState,useWeights);
    plot("Bs_TAU","#tau(B) [ns]",40,0,0.012,Year, finalState,useWeights);

    /// BDT
    plot("DTF_CHI2NDOF","DTF CHI2",40,0.,7,Year, finalState,useWeights);    
    plot("Bs_IPCHI2_OWNPV","#chi^{2}_{IP}(B)",40,0,16,Year, finalState,useWeights);
    plot("Bs_DIRA_OWNPV","#chi^{2}_{IP}(B)",40,0.99997,1,Year, finalState,useWeights);

    plot("XsDaughters_min_IPCHI2","X_{s} min(#chi^{2}_{IP})",40, 0, 10 ,Year, finalState,useWeights);
    plot("a_1_1260_plus_ptasy_1.00","Xs_ptasy_1.00",40, -1, 1 ,Year, finalState,useWeights);
    plot("Xs_max_DOCA","X_{s} max DOCA [mm]",40, 0, 0.4 ,Year, finalState,useWeights);

    plot("DsDaughters_min_IPCHI2","D_{s} min(#chi^{2}_{IP})",40, 0, 10 ,Year, finalState,useWeights);
    plot("Ds_ptasy_1.00","Ds_ptasy_1.00",40, -1, 1 ,Year, finalState,useWeights);
    plot("Ds_FDCHI2_ORIVX","#chi^{2}_{FD}(D_{s})",40,0,40000,Year, finalState,useWeights);
    plot("Ds_RFD","Ds RFD",40,0,10,Year, finalState,useWeights);

    plot("maxCos","maxCos",40,-1,1,Year, finalState,useWeights);    
    plot("max_ghostProb","max(Track_ghostProb)",40,0.,0.375,Year, finalState,useWeights);

    plot("Ds_m12","Ds_m12",40,900,1900,Year, finalState,useWeights);
    plot("Ds_m13","Ds_m13",40,600,1600,Year, finalState,useWeights);


    /*
    ///K
    plot("K_plus_P","p(K^{+}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("K_plus_PT","p_{T}(K^{+}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("K_plus_ETA","#eta(K^{+})",40,1,6,Year, finalState,useWeights);
    plot("K_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",40,0,13000,Year, finalState,useWeights);
    //plot("Kplus_PIDK","DLL_{K#pi}(K^{+}) ",100,-100,100);
    plot("angK","#theta_{D_{s} K^{+}}",40,0,3.141,Year, finalState,useWeights);
    plot("K_plus_ptasy_1.00","pt cone asymmetry (K^{+})",40,-1,1,Year, finalState,useWeights);
    plot("K_plus_TRACK_GhostProb","ghost prob (K^{+})",40,0,0.4,Year, finalState,useWeights);
    ///pi+
    plot("pi_plus_P","p(#pi^{+}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("pi_plus_PT","p_{T}(#pi^{+}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("pi_plus_ETA","#eta(#pi^{+})",40,1,6,Year, finalState,useWeights);
    plot("pi_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",40,0,13000,Year, finalState,useWeights);
    //plot("piplus_PIDK","DLL_{K#pi}(#pi^{+}) ",100,-100,100);
    plot("angPip","#theta_{D_{s} #pi^{+}}",40,0,3.141,Year, finalState,useWeights);
    plot("pi_plus_ptasy_1.00","pt cone asymmetry (#pi^{+})",40,0,1,Year, finalState,useWeights);
    plot("pi_plus_TRACK_GhostProb","ghost prob (#pi^{+})",40,0,0.4,Year, finalState,useWeights);
    ///pi-
    plot("pi_minus_P","p(#pi^{-}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("pi_minus_PT","p_{T}(#pi^{-}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("pi_minus_ETA","#eta(#pi^{-})",40,1,6,Year, finalState,useWeights);
    plot("pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",40,0,13000,Year, finalState,useWeights);
    //plot("piminus_PIDK","DLL_{K#pi}(#pi^{-}) ",100,-100,100);
    plot("angPim","#theta_{D_{s} #pi^{-}}",40,0,3.141,Year, finalState,useWeights);
    plot("pi_minus_ptasy_1.00","pt cone asymmetry (#pi^{-})",40,-1,1,Year, finalState,useWeights);
    plot("pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",40,0,0.4,Year, finalState,useWeights);
   
    ///Ds
    plot("Ds_P","p(D_{s}) [MeV]",40,0,400000,Year, finalState,useWeights);
    plot("Ds_PT","p_{T}(D_{s}) [MeV]",40,0,40000,Year, finalState,useWeights);
    plot("Ds_ETA","#eta(D_{s})",40,1,6,Year, finalState,useWeights);
    plot("Ds_DIRA_OWNPV","cos(DIRA) (D_{s})",40,0.9999,1,Year, finalState,useWeights);
     */
    /*
    ///K+ from Ds
    plot("K_plus_fromDs_P","p(K^{+} from D_{s}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("K_plus_fromDs_PT","p_{T}(K^{+} from D_{s}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("K_plus_fromDs_ETA","#eta(K^{+} from D_{s})",40,1,6,Year, finalState,useWeights);
    plot("K_plus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+} from D_{s})",40,0,13000,Year, finalState,useWeights);
    plot("K_plus_fromDs_ptasy_1.00","pt cone asymmetry (K^{+} from D_{s})",40,-1,1,Year, finalState,useWeights);
    plot("K_plus_fromDs_TRACK_GhostProb","ghost prob (K^{+} from D_{s})",40,0,0.4,Year, finalState,useWeights);
    ///pi- from Ds
    plot("pi_minus_fromDs_P","p(#pi^{-} from D_{s}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("pi_minus_fromDs_PT","p_{T}(#pi^{-} from D_{s}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("pi_minus_fromDs_ETA","#eta(#pi^{-} from D_{s})",40,1,6,Year, finalState,useWeights);
    plot("pi_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-} from D_{s})",40,0,13000,Year, finalState,useWeights);
    plot("pi_minus_fromDs_ptasy_1.00","pt cone asymmetry (#pi^{-} from D_{s})",40,-1,1,Year, finalState,useWeights);
    plot("pi_minus_fromDs_TRACK_GhostProb","ghost prob (#pi^{-} from D_{s})",40,0,0.4,Year, finalState,useWeights);
    ///K- from Ds
    plot("K_minus_fromDs_P","p(K^{-} from D_{s}) [MeV]",40,0,180000,Year, finalState,useWeights);
    plot("K_minus_fromDs_PT","p_{T}(K^{-} from D_{s}) [MeV]",40,0,10000,Year, finalState,useWeights);
    plot("K_minus_fromDs_ETA","#eta(K^{-} from D_{s})",40,1,6,Year, finalState,useWeights);
    plot("K_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-} from D_{s})",40,0,13000,Year, finalState,useWeights);
    plot("K_minus_fromDs_ptasy_1.00","pt cone asymmetry (K^{-} from D_{s})",40,-1,1,Year, finalState,useWeights);
    plot("K_minus_fromDs_TRACK_GhostProb","ghost prob (K^{-} from D_{s})",40,0,0.4,Year, finalState,useWeights);
    */
    
    //plotEventVars("nPV","N_{PV}",10,0,10,Year, finalState,useWeights);
    plot("NTracks","# of tracks",40,0,450,Year, finalState,useWeights);
}

void plotPID(string Branch,string TitleX, int bins, double min, double max, bool useWeights=true) {
    
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
   
void applyCorrectionHisto(vector<TString> vars, vector<double> min, vector<double> max, int Year, TString FinalState = "KKpi", TString Decay = "norm"){

    const int dim = vars.size();
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;
    HyperHistogram hist_weights("correctionHistos/weights" + label + "_" + FinalState + "_" + anythingToString(Year) + ".root");

    TFile* f = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/"+ Decay + "_Ds2" + FinalState + "_" + anythingToString(Year) + ".root","UPDATE");
    TTree* treeMC =dynamic_cast<TTree*>(f->Get("DecayTree"));

    vector<double> var_MC(dim,0.);
    double weight;
    int Ds_finalState_MC, year_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    treeMC->SetBranchAddress("weight",&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    treeMC->SetBranchAddress("year",&year_MC);

    vector<double> old_weights;

    for (int i = 0; i < treeMC->GetEntries(); i++){
	treeMC->GetEntry(i);
    	old_weights.push_back(weight);
    }

    treeMC->SetBranchStatus("weight",0);
    TTree* summary_tree = treeMC->CloneTree();
    TBranch* b_w = summary_tree->Branch("weight",&weight,"weight/D"); 

    HyperPointSet points_MC( dim );

    for (int i = 0; i < treeMC->GetEntries(); i++){
    
        treeMC->GetEntry(i);

        HyperPoint point( dim );
	for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 

	double w = 0.;
        int bin = hist_weights.getBinning().getBinNum(point);
            if(hist_weights.checkBinNumber(bin)!= bin){
		 w = 0; //? should't happen
		 cout << "ERROR:: Event outside limits" << endl;
            }else w = hist_weights.getBinContent(bin);
	
	weight = w * old_weights[i];
	if(w < 0) {
		w = 0.;
		cout << "ERROR:: Negative weight" << endl;
	}

	b_w->Fill();
    }

   summary_tree->Write();

   f->Close();

   return;
}

void resetWeights(int Year, TString FinalState = "KKpi"){

    TFile* f;
    if( FinalState == "KKpi" ){    
	if(Year == 11)f = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11.root","UPDATE");
	else if (Year == 12)f = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12.root","UPDATE");
    }
    else {
	cout << "Error::No file found";
	throw "EROOR";
    }

    TTree* treeMC =dynamic_cast<TTree*>(f->Get("DecayTree"));
    treeMC->SetBranchStatus("weight",0);

    TTree* summary_tree = treeMC->CloneTree();
    double weight = 1.;
    TBranch* b_w = summary_tree->Branch("weight",&weight,"weight/D"); 

    for (int i = 0; i < treeMC->GetEntries(); i++) b_w->Fill();
    
    summary_tree->Write();

    f->Close();

    return;
}


void produceCorrectionHisto(vector<TString> vars, vector<double> min, vector<double> max, int Year, TString FinalState = "KKpi"){

    /// Options
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 60); 
    NamedParameter<int> maxBinsPerDim("maxBinsPerDim", 200); 

    /// Check vector sizes
    if(vars.size() * min.size() * max.size() != pow(vars.size(),3) ){
	cout << "ERROR:: Different number of vars and limits";
	throw "ERROR";
    }

    /// Get dimension and minimum bin width
    const int dim = vars.size();
    vector<double> vec_minBinWidths(dim,0.);
    for(int i = 0; i < dim; i++)vec_minBinWidths[i]= (max[i]-min[i])/(double)maxBinsPerDim;
    HyperPoint minBinWidths(vec_minBinWidths);    

    HyperPoint Min(min);
    HyperPoint Max(max);
    HyperCuboid limits(Min, Max );

    /// Get data           
    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_sweight.root");
    std::string inputFile = InputFileName;
    TFile *File =  TFile::Open(inputFile.c_str());
    TTree* tree =dynamic_cast<TTree*>(File->Get("DecayTree"));
    cout << "reading events from file " << inputFile.c_str() << endl;
    cout << " I've got " << tree->GetEntries() << " events." << endl;
    
    tree->SetBranchStatus("*",0);
    for(int i = 0; i < dim; i++)tree->SetBranchStatus(vars[i],1);
    tree->SetBranchStatus("N_Bs_sw",1);
    tree->SetBranchStatus("Ds_finalState",1);
    tree->SetBranchStatus("year",1);

    vector<double> var_data(dim,0.);
    double sw;
    int Ds_finalState, year;
    for(int i = 0; i < dim; i++)tree->SetBranchAddress(vars[i],&var_data[i]);
    tree->SetBranchAddress("N_Bs_sw",&sw);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree->SetBranchAddress("year",&year);

    HyperPointSet points( dim );
    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);
	
	if(Year != year) continue;
    	if(FinalState == "KKpi" && Ds_finalState == 3) continue;

        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_data[j];
        point.addWeight(sw);
        points.push_back(point);
    }
    
    /// Get MC
    TChain* treeMC =new TChain("DecayTree");
    treeMC->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11.root");
    treeMC->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12.root");
    
    treeMC->SetBranchStatus("*",0);
    for(int j = 0; j < dim; j++)treeMC->SetBranchStatus(vars[j],1);
    treeMC->SetBranchStatus("weight",1);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("year",1);

    vector<double> var_MC(dim,0.);
    double weight;
    int Ds_finalState_MC, year_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    treeMC->SetBranchAddress("weight",&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    treeMC->SetBranchAddress("year",&year_MC);

    HyperPointSet points_MC( dim );
    for (int i = 0; i < treeMC->GetEntries(); i++){
    
        treeMC->GetEntry(i);

	if(Year != year_MC) continue;
    	if(FinalState == "KKpi" && Ds_finalState_MC == 3) continue;

        HyperPoint point( dim );
	for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 
        point.addWeight(weight);
        points_MC.push_back(point);
    }

    /// Define binning based on MC
    HyperHistogram histMC(limits, points_MC, 
                         
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI, 
                         
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),    
                                             
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (minBinWidths),
                                                 
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

    /// Draw binning
    histMC.setNames(HyperName(vars));
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;

    if(dim < 3)histMC.draw("correctionHistos/binning" + label + "_" + FinalState + "_" + anythingToString(Year) );

    /// Draw density
    HyperHistogram hist( histMC.getBinning() );
    hist.fill(points); 
    hist.normalise(1);
    histMC.normalise(1);
       
    if(dim < 3)histMC.drawDensity("correctionHistos/density_mc" + label + "_" + FinalState + "_" + anythingToString(Year) );
    if(dim < 3)hist.drawDensity("correctionHistos/density_data" + label + "_" + FinalState + "_" + anythingToString(Year) );

    /// Produce MC correction histo 
    hist.divide(histMC);
    if(dim < 3)hist.drawDensity("correctionHistos/weights" + label + "_" + FinalState + "_" + anythingToString(Year) );
    hist.save("correctionHistos/weights" + label + "_" + FinalState + "_" + anythingToString(Year) + ".root" );
   
    File->Close();

  return;
}
 


int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
    /// Options
    NamedParameter<int> nIterations("nIterations", 1); 

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);
    
    vector<int> years;
    years.push_back(11);
    //year.push_back(12);
    //year.push_back(15);
    //year.push_back(16);

    vector<TString> Ds_finalStates;
    Ds_finalStates.push_back("KKpi");
    //Ds_finalStates.push_back("pipipi");

    /// Define reweighting vars
    vector<TString>  vars_1, vars_2, vars_3, vars_4;
    vector<double>   min_1, min_2, min_3, min_4;
    vector<double>   max_1, max_2, max_3, max_4;
     
    vars_1.push_back("Bs_PT");
    min_1.push_back(0.);
    max_1.push_back(60000.);
    vars_1.push_back("Bs_ETA");
    min_1.push_back(1.5);
    max_1.push_back(5.5);
    //vars_1.push_back("NTracks");
    //min_1.push_back(0.);
    //max_1.push_back(600.);
    
    vars_2.push_back("NTracks");
    min_2.push_back(0.);
    max_2.push_back(600.);
    vars_2.push_back("max_ghostProb");
    min_2.push_back(0.);
    max_2.push_back(0.4);

    vars_3.push_back("NTracks");
    min_3.push_back(0.);
    max_3.push_back(600.);
    vars_3.push_back("Bs_PT");
    min_3.push_back(0.);
    max_3.push_back(60000.);

    vars_4.push_back("DTF_CHI2NDOF");
    min_4.push_back(0.);
    max_4.push_back(10.);

    vector< vector<TString> > vars_set;
    vector< vector<double> > min_set;
    vector< vector<double> > max_set;

    vars_set.push_back(vars_1);
    vars_set.push_back(vars_2);
    vars_set.push_back(vars_3);
    vars_set.push_back(vars_4);

    min_set.push_back(min_1);
    min_set.push_back(min_2);
    min_set.push_back(min_3);
    min_set.push_back(min_4);

    max_set.push_back(max_1);
    max_set.push_back(max_2);
    max_set.push_back(max_3);
    max_set.push_back(max_4);

    /// Reset weights to 1
    for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++)resetWeights(years[i],Ds_finalStates[j]);

    /// Produce MC correction histos and apply weights
    /// Weights are applied on top of each other with the previous weighting applied
    for(int n= 0; n < nIterations; n++)
	for(int i= 0; i < years.size(); i++) 
		for(int j= 0; j < Ds_finalStates.size(); j++)		
			for(int k =0; k < vars_set.size(); k++)	{	
				produceCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j]);
				applyCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j],"norm");    
				//applyCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j],"signal");    
			}
	
    /// Draw comparison plots 
    for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++){ 
		dataVsMC(years[i],Ds_finalStates[j], true,"norm");
		//dataVsMC(years[i],Ds_finalStates[j], true,"signal");
    }

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
