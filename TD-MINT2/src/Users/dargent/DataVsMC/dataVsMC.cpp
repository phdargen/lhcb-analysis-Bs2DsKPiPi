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
#include <TLegendEntry.h>
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

NamedParameter<string> ReweightFromA("ReweightFromA", (std::string) "/auto/data/dargent/BsDsKpipi/Preselected/MC/norm.root");
NamedParameter<string> ReweightToB("ReweightToB", (std::string) "/auto/data/dargent/BsDsKpipi/Preselected/Data/norm.root");
NamedParameter<string> ApplyWeightToC("ApplyWeightToC", (std::string) "");

NamedParameter<string> weightVarA("weightVarA", (std::string) "weight");
NamedParameter<string> newWeightVarA("newWeightVarA", (std::string) "weight");
NamedParameter<string> weightVarB("weightVarB", (std::string) "N_Bs_sw");
NamedParameter<string> weightVarC("weightVarC", (std::string) "noweight");
NamedParameter<string> newWeightVarC("newWeightVarC", (std::string) "noweight");

NamedParameter<string> legTitle("legTitle", (std::string) "");
NamedParameter<string> nameA("nameA", (std::string) "MC");
NamedParameter<string> nameB("nameB", (std::string) "Data");

NamedParameter<string> cutA("cutA", (std::string) "");
NamedParameter<string> cutB("cutB", (std::string) "");

NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);
NamedParameter<int> updateAnaNoteHistos("updateAnaNoteHistos", 0);

void plot(TTree* tree, TTree* treeMC, TString Branch,TString TitleX, int bins, double min, double max, TString weightA, TString weightB, TString newWeightB, TString label, bool log = false, bool legendLeft = false){
        
    cout << "Plotting " << Branch << endl;

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus(Branch,1);
    tree->SetBranchStatus(weightA,1);
    double var;
    float varF[100];
    Short_t varS;
    int varI;
    double sw;
    if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")tree->SetBranchAddress(Branch,&varF);
    else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" )tree->SetBranchAddress(Branch,&varI);
    else if( A_is_in_B("_DEC",(string)Branch))tree->SetBranchAddress(Branch,&varS);
    else tree->SetBranchAddress(Branch,&var);
    tree->SetBranchAddress(weightA,&sw);
    
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus(Branch,1);
    if(weightB != "noweight")treeMC->SetBranchStatus(weightB,1);
    if(newWeightB != "noweight")treeMC->SetBranchStatus(newWeightB,1);
    double varMC;
    float varMCF[100];
    Short_t varMCS;
    int varMCI;
    double w = 1;
    double new_w = 1;
    if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")treeMC->SetBranchAddress(Branch,&varMCF);
    else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" )treeMC->SetBranchAddress(Branch,&varMCI);
    else if( A_is_in_B("_DEC",(string)Branch))treeMC->SetBranchAddress(Branch,&varMCS);
    else treeMC->SetBranchAddress(Branch,&varMC);         
    if(weightB != "noweight")treeMC->SetBranchAddress(weightB,&w);
    if(newWeightB != "noweight")treeMC->SetBranchAddress(newWeightB,&new_w);

    ///Make histograms
    TString title= ";"+TitleX+";Yield [norm.]";
    TH1D h(Branch,title,bins,min,max);
    TH1D h_MC(Branch+"_MC",title,bins,min,max);
    TH1D h_MC_rw(Branch+"_MC_rw",title,bins,min,max);
    
    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
        if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")var = (double)varF[0];
        else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" )var = (double)varI;
        else if( A_is_in_B("_DEC",(string)Branch))var = (double)varS;
        h.Fill(var,sw);
    }
    
    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
        treeMC->GetEntry(i);
        if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB"|| Branch == "Bs_DTF_MERR")varMC = (double)varMCF[0];
        else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" )varMC = (double)varMCI;
        else if( A_is_in_B("_DEC",(string)Branch))varMC = (double)varMCS;
        h_MC.Fill(varMC,w);
        h_MC_rw.Fill(varMC,new_w);
    }
    
    ///Plot it
    TCanvas c;
    
    h.Scale(1./h.Integral());
    h_MC.Scale(1./h_MC.Integral());
    h_MC_rw.Scale(1./h_MC_rw.Integral());
    double maxY= h.GetMaximum();
    if(h_MC.GetMaximum()>maxY)maxY=h_MC.GetMaximum();
    h.SetMinimum(0.);
    if(log){
        h.SetMinimum(0.0001);
        gPad->SetLogy(1);
    }
    else gPad->SetLogy(0);
    h.SetMaximum(maxY*1.4);
    h.SetLineColor(kBlack);
    h.Draw("");
    h_MC.SetMarkerColor(kRed);
    h_MC.SetLineColor(kRed);
    h_MC.Draw("esame");
    h_MC_rw.SetLineColor(kBlue);
    h_MC_rw.SetMarkerColor(kBlue);
    if(newWeightB != weightB && newWeightB != "noweight")h_MC_rw.Draw("esame");
    
    double KolmoTest = h.KolmogorovTest(&h_MC);
    double KolmoTest_rw = h.KolmogorovTest(&h_MC_rw);

    TLegend* leg;
    if(legendLeft)leg = new TLegend(0.15,0.6,0.45,0.9,"");
    else leg = new TLegend(0.55,0.6,0.85,0.9,"");
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextFont(22);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    leg->SetTextAlign(12);

    if((string)legTitle != "")leg->AddEntry((TObject*)0,((string)legTitle).c_str(), "");
    leg->AddEntry(&h,((string)nameB).c_str(),"LEP");
    
    TLegendEntry* le = leg->AddEntry(&h_MC,((string)nameA).c_str(),"LEP");
    le->SetTextColor(kRed);    

    stringstream ss ;
    TString leg_kol = "KS-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest ;
    leg_kol += ss.str();    
    le = leg->AddEntry((TObject*)0, leg_kol, "");
    le->SetTextColor(kRed);    

    if(newWeightB != weightB && newWeightB != "noweight"){
        le = leg->AddEntry(&h_MC_rw, ((string)nameA + " (reweighted)").c_str(),"LEP");
        le->SetTextColor(kBlue);  
    }
    ss.str("");
    leg_kol = "KS-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest_rw ;
    leg_kol += ss.str();    
    if(newWeightB != weightB && newWeightB != "noweight"){
        TLegendEntry* le = leg->AddEntry((TObject*)0, leg_kol, "");
        le->SetTextColor(kBlue);    
    }
    leg->Draw(); 
    
    cout << endl;
    c.Print((string)OutputDir + label + "_"+Branch+".eps");
    if(updateAnaNotePlots)c.Print("../../../../../TD-AnaNote/latex/figs/dataVsMC/" + (string)OutputDir + label + "_"+Branch+".pdf" );

}

void compare(TString fileA, TString fileB, TString weightA, TString weightB, TString newWeightB, TString CutA = "", TString CutB = "", int Year = -1, TString finalState = "all", int Trigger = -1, TString label = ""){
    
    // Cuts
    TString Cut;
    if(Year>10)Cut += " year == " + anythingToString(Year);
    else if(Year == -1)Cut += " year > " + anythingToString(Year);
    else Cut += " run == " + anythingToString(Year);
    if(finalState == "KKpi")Cut += " && Ds_finalState < 3 ";
    else if(finalState == "pipipi")Cut += " && Ds_finalState == 3 ";
    else if(finalState == "Kpipi")Cut += " && Ds_finalState == 4 ";
    if(Trigger != -1) Cut += " && TriggerCat == " + anythingToString(Trigger);
    
    if(CutA != "")CutA += " && ";
    CutA += Cut;
    if(CutB != "")CutB += " && ";
    CutB += Cut;
    
    cout << endl << "Comparing file " << endl << fileA << " ( " << CutA << " ) " << endl; 
    cout << " to " << endl << fileB << " ( " << CutB << " ) " << endl << endl;
    
    ///Load files
    TChain* treeA = new TChain("DecayTree");
    treeA->Add(fileA);
    
    TChain* treeB= new TChain("DecayTree");
    treeB->Add(fileB);
    
    TFile* output = new TFile("dummy.root","RECREATE");
    TTree* new_treeA = treeA->CopyTree(CutA);
    TTree* new_treeB = treeB->CopyTree(CutB);

    /// Options
    NamedParameter<int> nBins("nBins", 40); 
    
    label += "Ds2";
    if(finalState != "")label +=  finalState;
    else label += "all";
    if(Year>-1) label += "_" + anythingToString(Year);
    if(Trigger>-1) label += "_t" + anythingToString(Trigger);
    
    TString Decay, selection;
    if(A_is_in_B("signal",(string)fileA) && A_is_in_B("signal",(string)fileB)) Decay = "signal";
    if(A_is_in_B("norm",(string)fileA) && A_is_in_B("norm",(string)fileB)) Decay = "norm";
    if(A_is_in_B("Final",(string)fileA) && A_is_in_B("Final",(string)fileB)) selection = "Final";

    /// Bs
    plot(new_treeA,new_treeB,"Bs_PT","p_{T}(B) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_P","p(B) [MeV]",nBins,0,900000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_ETA","#eta(B)",nBins,1,6,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_FDCHI2_OWNPV","#chi^{2}_{FD}(B)",nBins,0,100000,weightA, weightB, newWeightB, label,true);
    plot(new_treeA,new_treeB,"Bs_ENDVERTEX_CHI2","#chi^{2}_{vtx}(B)",nBins,0,35,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_DTF_TAU","t(B) [ns]",nBins,0.,10.,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_DTF_TAUERR","#sigma_{t}(B) [ns]",nBins,0,0.15,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_ptasy_1.00","B_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"Bs_DTF_MERR","#sigma_{m}",nBins,4.,25.,weightA, weightB, newWeightB, label);

    /// BDT
    plot(new_treeA,new_treeB,"DTF_CHI2NDOF","DTF #chi^{2}",nBins,0.,7,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"Bs_IPCHI2_OWNPV","#chi^{2}_{IP}(B)",nBins,0,20,weightA, weightB, newWeightB, label,true);
    plot(new_treeA,new_treeB,"Bs_DIRA_OWNPV","DIRA(B)",nBins,0.99997,1,weightA, weightB, newWeightB, label,true,true);

    plot(new_treeA,new_treeB,"XsDaughters_min_IPCHI2","X_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,weightA, weightB, newWeightB, label,true);
    if(Decay == "norm")plot(new_treeA,new_treeB,"a_1_1260_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 1. ,weightA, weightB, newWeightB, label,false,true);
    else plot(new_treeA,new_treeB,"K_1_1270_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"Xs_max_DOCA","X_{s} max DOCA [mm]",nBins, 0, 0.4 ,weightA, weightB, newWeightB, label);

    plot(new_treeA,new_treeB,"DsDaughters_min_IPCHI2","D_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,weightA, weightB, newWeightB, label,true);
    plot(new_treeA,new_treeB,"Ds_ptasy_1.00","Ds_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"Ds_FDCHI2_ORIVX","#chi^{2}_{FD}(D_{s})",nBins,0,40000,weightA, weightB, newWeightB, label,true);
    plot(new_treeA,new_treeB,"Ds_RFD","Ds RFD",nBins,0,10,weightA, weightB, newWeightB, label);

    plot(new_treeA,new_treeB,"maxCos","maxCos",nBins,-1,1,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"max_ghostProb","max(Track_ghostProb)",nBins,0,0.4,weightA, weightB, newWeightB, label);
    
    /// Tagging
    plot(new_treeA,new_treeB,"Bs_TAGDECISION_OS","q_{OS}",8,-1.5,6.5,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"Bs_TAGOMEGA_OS","#eta_{OS}",nBins,0.,0.499999,weightA, weightB, newWeightB, label, false, true);    
    plot(new_treeA,new_treeB,"Bs_SS_nnetKaon_DEC","q_{SS}",8,-1.5,6.5,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"Bs_SS_nnetKaon_PROB","#eta_{SS}",nBins,0.,0.499999,weightA, weightB, newWeightB, label, false, true);  
    
    /// Ds
    plot(new_treeA,new_treeB,"Ds_m12","Ds_m12",nBins,900,1900,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_m13","Ds_m13",nBins,600,1600,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_PT","p_{T}(D_{s}) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_ETA","#eta(D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_finalState","D_{s} final state",8,-0.5,7.5,weightA, weightB, newWeightB, label);

    /// Trigger
    plot(new_treeA,new_treeB,"TriggerCat","Trigger category",8,-0.5,7.5,weightA, weightB, newWeightB, label);
    
    /// Dalitz
    plot(new_treeA,new_treeB,"m_Kpipi","m(K^{+}#pi^{+}#pi^{-})[MeV]",nBins,0,1950,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_Kpi","m(K^{+}#pi^{-})[MeV]",nBins,0,1200,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_pipi","m(#pi^{+}#pi^{-})[MeV]",nBins,0,1200,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_Dspipi","m(D_{s}^{-}#pi^{+}#pi^{-})[MeV]",nBins,1900,5550,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_Dspi","m(D_{s}^{-}#pi^{+})[MeV]",nBins,0,5500,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_DsK","m(D_{s}^{-}K^{+})[MeV]",nBins,0,5500,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_DsKpi","m(D_{s}^{-}K^{+}#pi^{-})[MeV]",nBins,1900,5500,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_DsKpip","m(D_{s}^{-}K^{+}#pi^{+})[MeV]",nBins,1900,5500,weightA, weightB, newWeightB, label,false,true);

    if(Decay == "signal"){
        plot(new_treeA,new_treeB,"K_plus_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_ETA","#eta(K^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"K_plus_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"K_plus_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_plus_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus_ETA","#eta(#pi^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"pi_plus_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_plus_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,weightA, weightB, newWeightB, label);
    }   
    else if(Decay == "norm") {
        plot(new_treeA,new_treeB,"pi_plus1_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus1_ETA","#eta(K^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus1_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"pi_plus1_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
       // plot(new_treeA,new_treeB,"pi_plus1_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_plus2_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus2_ETA","#eta(#pi^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_plus2_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"pi_plus2_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_plus2_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        plot(new_treeA,new_treeB,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,weightA, weightB, newWeightB, label);
    }
    if(finalState == "KKpi") {
        plot(new_treeA,new_treeB,"K_plus_fromDs_PT","p_{T}(K^{+} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_ETA","#eta(K^{+} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"K_plus_fromDs_TRACK_GhostProb","ghost prob (K^{+} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_minus_fromDs_PT","p_{T}(#pi^{-} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_ETA","#eta(#pi^{-} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"pi_minus_fromDs_TRACK_GhostProb","ghost prob (#pi^{-} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_PIDK","DLL_{K#pi}(#pi^{-} from D_{s}) ",nBins,-100,100,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"K_minus_fromDs_PT","p_{T}(K^{-} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_ETA","#eta(K^{-} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"K_minus_fromDs_TRACK_GhostProb","ghost prob (K^{-} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_PIDK","DLL_{K#pi}(K^{-} from D_{s}) ",nBins,-20,100,weightA, weightB, newWeightB, label);
    }    
    plot(new_treeA,new_treeB,"NTracks","# of tracks",nBins,0,550,weightA, weightB, newWeightB, label);
    if(selection == "Final") plot(new_treeA,new_treeB,"BDTG_response","BDTG",nBins,0,1.,weightA, weightB, newWeightB, label,false,true);
}

/*
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
*/
   
void applyCorrectionHisto(vector<TString> vars, int Year, TString FinalState, int Trigger, TString ApplyTo, TString& weightVar, TString NewWeightVar){

    NamedParameter<double> maxWeight("maxWeight",100.);
    const int dim = vars.size();

    double sumw = 0;
    double sumw2 = 0;
    double sumw_old = 0;
    double sumw2_old = 0;
    
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;
    label += "_Ds2" + FinalState;
    if(Year != -1)label += "_" + anythingToString(Year);
    if(Trigger != -1)label+= "_t" + anythingToString(Trigger);
    
    HyperHistogram hist_weights((string)OutputDir + "weights/weights" + label + ".root");

    TFile* f = new TFile(ApplyTo,"UPDATE");
    TTree* treeMC =dynamic_cast<TTree*>(f->Get("DecayTree"));

    cout << "Apply correction histo " << (string)OutputDir + "weights/weights" + label + ".root" << endl;
    cout << "to file " << ApplyTo << endl;
    
    vector<double> var_MC(dim,0.);
    double weight = 1;
    int Ds_finalState_MC, year_MC, trigger_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    if(weightVar != "noweight")treeMC->SetBranchAddress(weightVar,&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    if(Year>10)treeMC->SetBranchAddress("year",&year_MC);
    else treeMC->SetBranchAddress("run",&year_MC);
    treeMC->SetBranchAddress("TriggerCat",&trigger_MC);

    vector<double> old_weights;
    for (int i = 0; i < treeMC->GetEntries(); i++){
        treeMC->GetEntry(i);
    	old_weights.push_back(weight);
        sumw_old += weight;
        sumw2_old += weight*weight;
    }

    TBranch* br = (TBranch*)treeMC->GetListOfBranches()->FindObject(NewWeightVar);
    if(br != 0)treeMC->SetBranchStatus(NewWeightVar,0);
    if(NewWeightVar != weightVar)weightVar = NewWeightVar;
    TTree* summary_tree = treeMC->CloneTree();
    TBranch* b_w = summary_tree->Branch(NewWeightVar,&weight,NewWeightVar+"/D"); 

    HyperPointSet points_MC( dim );

    for (int i = 0; i < treeMC->GetEntries(); i++){
        treeMC->GetEntry(i);

        if( (Year != year_MC && Year != -1) 
           || (FinalState == "KKpi" && Ds_finalState_MC > 2) || (FinalState == "pipipi" && Ds_finalState_MC > 3) 
           || (Trigger != trigger_MC && Trigger != -1) ) {
            weight = old_weights[i];
            b_w->Fill();
            sumw += weight;
            sumw2 += weight*weight;
            continue;
        }
        else if(FinalState == "Kpipi" && Ds_finalState_MC > 4) throw "undefined final state";
        
        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 

        double w = 1.;
        int bin = hist_weights.getBinning().getBinNum(point);
        if(hist_weights.checkBinNumber(bin)!= bin){
            w = 1; //? should't happen
            cout << "ERROR:: Event outside limits" << endl;
        }else w = std::fmin(maxWeight,hist_weights.getBinContent(bin));
	
        if(w < 0) {
            w = 0.;
            cout << "ERROR:: Negative weight" << endl;
        }
        
        weight = w * old_weights[i];
        b_w->Fill();
        sumw += weight;
        sumw2 += weight*weight;
    }

   cout << "Effective weight before reweighting = " << sumw_old/sumw2_old << endl; 
   cout << "Effective weight after reweighting = " << sumw/sumw2 << endl; 

   summary_tree->Write();
   f->Close();

   return;
}

void produceCorrectionHisto(vector<TString> vars, vector<double> min, vector<double> max, int Year, TString FinalState, int Trigger, TString weightA){

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
    TChain* tree = new TChain("DecayTree");
    tree->Add(((string)ReweightToB).c_str());
    
    tree->SetBranchStatus("*",0);
    for(int i = 0; i < dim; i++)tree->SetBranchStatus(vars[i],1);
    if((string)weightVarB != "noweight")tree->SetBranchStatus(((string)weightVarB).c_str(),1);
    tree->SetBranchStatus("Ds_finalState",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("TriggerCat",1);

    vector<double> var_data(dim,0.);
    double sw = 1;
    int Ds_finalState, year, trigger;
    for(int i = 0; i < dim; i++)tree->SetBranchAddress(vars[i],&var_data[i]);
    if((string)weightVarB != "noweight")tree->SetBranchAddress(((string)weightVarB).c_str(),&sw);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
    if(Year>10)tree->SetBranchAddress("year",&year);
    else tree->SetBranchAddress("run",&year);
    tree->SetBranchAddress("TriggerCat",&trigger);

    HyperPointSet points( dim );
    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);
	
        if(Year != year && Year != -1) continue;
    	if(FinalState == "KKpi" && Ds_finalState > 2) continue;
        else if(FinalState == "pipipi" && Ds_finalState > 3) continue;
        else if(FinalState == "Kpipi" && Ds_finalState > 4) throw "undefined final state";
        if(Trigger != trigger && Trigger != -1) continue;

        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_data[j];
        point.addWeight(sw);
        points.push_back(point);
    }
    
    /// Get MC
    TChain* treeMC =new TChain("DecayTree");
    treeMC->Add(((string)ReweightFromA).c_str());
    
    treeMC->SetBranchStatus("*",0);
    for(int j = 0; j < dim; j++)treeMC->SetBranchStatus(vars[j],1);
    if(weightA != "noweight")treeMC->SetBranchStatus(weightA,1);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("year",1);
    treeMC->SetBranchStatus("run",1);
    treeMC->SetBranchStatus("TriggerCat",1);

    vector<double> var_MC(dim,0.);
    double weight = 1;
    int Ds_finalState_MC, year_MC, trigger_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    if(weightA != "noweight")treeMC->SetBranchAddress(weightA,&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    if(Year>10)treeMC->SetBranchAddress("year",&year_MC);
    else treeMC->SetBranchAddress("run",&year_MC);
    treeMC->SetBranchAddress("TriggerCat",&trigger_MC);

    HyperPointSet points_MC( dim );
    for (int i = 0; i < treeMC->GetEntries(); i++){
    
        treeMC->GetEntry(i);

        if(Year != year_MC && Year != -1 ) continue;
        if(FinalState == "KKpi" && Ds_finalState_MC > 2) continue;
        else if(FinalState == "pipipi" && Ds_finalState_MC > 3) continue;
        else if(FinalState == "Kpipi" && Ds_finalState_MC > 4) throw "undefined final state";
        if(Trigger != trigger_MC && Trigger != -1) continue;

        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 
        point.addWeight(weight);
        points_MC.push_back(point);
    }

    /// Define binning based on sample with smaller statistic
    HyperHistogram* histMC;
    HyperHistogram* hist;
    
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;
    label += "_Ds2" + FinalState;
    if(Year != -1)label += "_" + anythingToString(Year);
    if(Trigger != -1)label+= "_t" + anythingToString(Trigger);
    
    if(points.getSumW()>points_MC.getSumW()){
         histMC= new HyperHistogram(limits, points_MC, 
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
            histMC->setNames(HyperName(vars));
            hist= new HyperHistogram( histMC->getBinning() );
            hist->fill(points); 
            /// Draw binning
            if(dim < 3)histMC->draw((string)OutputDir + "weights/binning" + label);
    }
    else {
        hist= new HyperHistogram(limits, points, 
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
        hist->setNames(HyperName(vars));
        histMC= new HyperHistogram( hist->getBinning() );
        histMC->fill(points_MC); 
        /// Draw binning
        if(dim < 3)hist->draw((string)OutputDir + "weights/binning" + label );
    }
        
    /// Draw density
    hist->normalise(1);
    histMC->normalise(1);
    if(dim < 3)histMC->drawDensity((string)OutputDir + "weights/density_mc" + label );
    if(dim < 3)hist->drawDensity((string)OutputDir + "weights/density_data" + label );

    /// Produce MC correction histo 
    hist->divide(*histMC);
    if(dim < 3){
        hist->draw((string)OutputDir + "weights/weights" + label );
        if(updateAnaNoteHistos)hist->draw("../../../../../TD-AnaNote/latex/figs/dataVsMC/" +(string)OutputDir + "weights/weights" + label );
    }
    hist->save((string)OutputDir + "weights/weights" + label + ".root" );
    
    return;
}

void createSubset(TString file, TString newfile, TString Cut){
    TChain* tree = new TChain("DecayTree");
    tree->Add(file);
    TFile* output = new TFile(newfile,"RECREATE");
    TTree* new_tree = tree->CopyTree(Cut);
    new_tree->Write();
    output->Close();
}


int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_t0.root","TriggerCat == 0");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_t1.root","TriggerCat == 1");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_r1.root","run == 1");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_r2.root","run == 2");

    /// Options
    NamedParameter<int> nIterations("nIterations", 1); 
    NamedParameter<int> reweight("reweight", 1); 
    NamedParameter<int> reweightInBinsOfRun("reweightInBinsOfRun", 1); 
    NamedParameter<int> reweightInBinsOfFinalState("reweightInBinsOfFinalState", 1); 
    NamedParameter<int> reweightInBinsOfTrigger("reweightInBinsOfTrigger", 1); 

    NamedParameter<int> reweightVarSet1("reweightVarSet1", 1); 
    NamedParameter<int> reweightVarSet2("reweightVarSet2", 1); 
    NamedParameter<int> reweightVarSet3("reweightVarSet3", 1); 
    NamedParameter<int> reweightVarSet4("reweightVarSet4", 1); 
    NamedParameter<int> reweightVarSet5("reweightVarSet5", 1); 

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);
    
    vector<int> years;
    if(reweightInBinsOfRun==1){
        years.push_back(1); // Means run 1
        //years.push_back(2); // Means run 2
    }
    else if(reweightInBinsOfRun==0){
        years.push_back(11);
        years.push_back(12);
        //years.push_back(15);
        //years.push_back(16);
    }
    else years.push_back(-1); // Means all 
    
    vector<TString> Ds_finalStates;
    if(reweightInBinsOfFinalState){
        Ds_finalStates.push_back("KKpi");
        //Ds_finalStates.push_back("pipipi");
        //Ds_finalStates.push_back("Kpipi");
    }
    else Ds_finalStates.push_back("all");
    
    vector<int> trigger;
    if(reweightInBinsOfTrigger){
        trigger.push_back(0);
        trigger.push_back(1);
    }
    else trigger.push_back(-1);

    /// Define reweighting vars
    vector<TString>  vars_1, vars_2, vars_3, vars_4, vars_5;
    vector<double>   min_1, min_2, min_3, min_4, min_5;
    vector<double>   max_1, max_2, max_3, max_4, max_5;
     
    vars_1.push_back("Bs_PT");
    min_1.push_back(0.);
    max_1.push_back(200000.);
    vars_1.push_back("Bs_ETA");
    min_1.push_back(1.5);
    max_1.push_back(5.5);
    //vars_1.push_back("NTracks");
    //min_1.push_back(0.);
    //max_1.push_back(1000.);
    
    vars_2.push_back("NTracks");
    min_2.push_back(0.);
    max_2.push_back(1000.);
    vars_2.push_back("max_ghostProb");
    min_2.push_back(0.);
    max_2.push_back(0.4);

    vars_3.push_back("Bs_DTF_TAUERR");
    min_3.push_back(0.);
    max_3.push_back(0.2);
    //vars_3.push_back("Bs_PT");
    //min_3.push_back(0.);
    //max_3.push_back(200000.);

    vars_4.push_back("DTF_CHI2NDOF");
    min_4.push_back(0.);
    max_4.push_back(10.);
    
    vars_5.push_back("Bs_IPCHI2_OWNPV");
    min_5.push_back(0.);
    max_5.push_back(20.);
    
    vector< vector<TString> > vars_set;
    vector< vector<double> > min_set;
    vector< vector<double> > max_set;

    if(reweightVarSet1)vars_set.push_back(vars_1);
    if(reweightVarSet2)vars_set.push_back(vars_2);
    if(reweightVarSet3)vars_set.push_back(vars_3);
    if(reweightVarSet4)vars_set.push_back(vars_4);
    if(reweightVarSet5)vars_set.push_back(vars_5);
    
    if(reweightVarSet1)min_set.push_back(min_1);
    if(reweightVarSet2)min_set.push_back(min_2);
    if(reweightVarSet3)min_set.push_back(min_3);
    if(reweightVarSet4)min_set.push_back(min_4);
    if(reweightVarSet5)min_set.push_back(min_5);

    if(reweightVarSet1)max_set.push_back(max_1);
    if(reweightVarSet2)max_set.push_back(max_2);
    if(reweightVarSet3)max_set.push_back(max_3);
    if(reweightVarSet4)max_set.push_back(max_4);
    if(reweightVarSet5)max_set.push_back(max_5);

    /// Produce MC correction histos and apply weights
    /// Weights are applied on top of each other with the previous weighting applied
    TString weightA = (string) weightVarA;
    TString weightC = (string) weightVarC;
    
    if(reweight)for(int n= 0; n < nIterations; n++)
			for(int i= 0; i < years.size(); i++) 
				for(int j= 0; j < Ds_finalStates.size(); j++)	
                    for(int k= 0; k < trigger.size(); k++)    
                        for(int l =0; l < vars_set.size(); l++)	{	
                            produceCorrectionHisto(vars_set[l],min_set[l],max_set[l],years[i],Ds_finalStates[j],trigger[k], weightA);
                            applyCorrectionHisto(vars_set[l],years[i],Ds_finalStates[j],trigger[k],(string)ReweightFromA, weightA, (string)newWeightVarA);    
                            if((string)ApplyWeightToC != "")applyCorrectionHisto(vars_set[l],years[i],Ds_finalStates[j],trigger[k],(string)ApplyWeightToC, weightC, (string)newWeightVarC);
					}
			
    /// Draw comparison plots 
    for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++)for(int k= 0; k < trigger.size(); k++){ 
        compare((string) ReweightToB, (string) ReweightFromA, (string) weightVarB, (string) weightVarA, (string) newWeightVarA, (string) cutB, (string) cutA, years[i], Ds_finalStates[j],trigger[k]);
    }

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
