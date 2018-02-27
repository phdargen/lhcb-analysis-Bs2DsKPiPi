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

TString plot(TString Branch,TString TitleX, int bins, double min, double max, int Year = 11, TString finalState = "KKpi", bool useWeights=false, TString Decay = "norm", TString selection = "Preselected"){
    
    cout << "Plotting " << Branch << endl;

    ///Load files
    TChain tree("DecayTree");
    TString decay = Decay;
    if(Decay == "signal") selection = "Final" ;
    if(selection == "Preselected") decay.ReplaceAll(Decay,Decay+"_sweight");
    tree.Add("/auto/data/dargent/BsDsKpipi/" + selection + "/Data/" + decay + ".root");
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus(Branch,1);
    tree.SetBranchStatus("N_Bs_sw",1);
    tree.SetBranchStatus("year",1);
    tree.SetBranchStatus("Ds_finalState",1);
    double var;
    float varF;
    double sw;
    int year,Ds_finalState;
    if(Branch == "BDTG_response")tree.SetBranchAddress(Branch,&varF);
    else tree.SetBranchAddress(Branch,&var);
    tree.SetBranchAddress("N_Bs_sw",&sw);
    tree.SetBranchAddress("year",&year);
    tree.SetBranchAddress("Ds_finalState",&Ds_finalState);
    
    TChain treeMC("DecayTree");
    TString fileNameMC = "/auto/data/dargent/BsDsKpipi/" + selection + "/MC/" + Decay + "_Ds2" + finalState + "_" + anythingToString(Year) + ".root";
    if(selection == "Final") fileNameMC = "/auto/data/dargent/BsDsKpipi/" + selection + "/MC/" + Decay + ".root";
    treeMC.Add(fileNameMC);
    treeMC.SetBranchStatus("*",0);
    treeMC.SetBranchStatus(Branch,1);
    treeMC.SetBranchStatus("Ds_finalState",1);
    treeMC.SetBranchStatus("Bs_BKGCAT",1);
    if(useWeights)treeMC.SetBranchStatus("weight",1);
   
    double varMC;
    float varMCF;
    double w;
    int cat,yearMC,Ds_finalStateMC;
    if(Branch == "BDTG_response"){treeMC.SetBranchAddress(Branch,&varMCF); cout << "true" << endl;}
    else treeMC.SetBranchAddress(Branch,&varMC);
    treeMC.SetBranchAddress("Bs_BKGCAT",&cat);
    treeMC.SetBranchAddress("year",&yearMC);
    treeMC.SetBranchAddress("Ds_finalState",&Ds_finalStateMC);           
    if(useWeights)treeMC.SetBranchAddress("weight",&w);
    else w=1.;
    
    ///Make histograms
    TString title= ";"+TitleX+";Yield [norm.]";
    TH1D h(Branch,title,bins,min,max);
    TH1D h_MC(Branch+"_MC",title,bins,min,max);
    TH1D h_MC_rw(Branch+"_MC_rw",title,bins,min,max);
    
    ///loop over data events
    int numEvents = tree.GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree.GetEntry(i);
	if(year != Year) continue;
	if(finalState == "KKpi" && Ds_finalState == 3) continue;
	if(Branch == "BDTG_response")var = (double)varF;
        h.Fill(var,sw);
    }
    
    ///loop over MC events
    int numEventsMC = treeMC.GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
        treeMC.GetEntry(i);
	if(yearMC != Year) continue;
	if(finalState == "KKpi" && Ds_finalStateMC == 3) continue;
	if(Branch == "BDTG_response")varMC = (double)varMCF;
        h_MC.Fill(varMC);
        h_MC_rw.Fill(varMC,w);
    }
    
    ///Plot it
    TCanvas c;
    
    h.Scale(1./h.Integral());
    h_MC.Scale(1./h_MC.Integral());
    h_MC_rw.Scale(1./h_MC_rw.Integral());
    double maxY= h.GetMaximum();
    if(h_MC.GetMaximum()>maxY)maxY=h_MC.GetMaximum();
    h.SetMinimum(0.);
    h.SetMaximum(maxY*1.4);
    h.SetLineColor(kBlack);
    h.Draw("");
    h_MC.SetMarkerColor(kRed);
    h_MC.SetLineColor(kRed);
    h_MC.Draw("esame");
    h_MC_rw.SetLineColor(kBlue);
    h_MC_rw.SetMarkerColor(kBlue);
    if(useWeights)h_MC_rw.Draw("esame");
    
    double KolmoTest = h.KolmogorovTest(&h_MC);
    double KolmoTest_rw = h.KolmogorovTest(&h_MC_rw);

    TLegend leg(0.6,0.6,0.9,0.9,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(22);
    leg.SetTextColor(1);
    leg.SetTextSize(0.04);
    leg.SetTextAlign(12);

    leg.AddEntry(&h,"Data","LEP");
    leg.AddEntry(&h_MC,"MC","LEP");

    stringstream ss ;
    TString leg_kol = "Kolm.-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest ;
    leg_kol += ss.str();    
    TLegendEntry* le = leg.AddEntry((TObject*)0, leg_kol, "");
    le->SetTextColor(kRed);    

    if(useWeights)leg.AddEntry(&h_MC_rw,"MC (reweighted)","LEP");
    ss.str("");
    leg_kol = "Kolm.-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest_rw ;
    leg_kol += ss.str();    
    if(useWeights){
	TLegendEntry* le = leg.AddEntry((TObject*)0, leg_kol, "");
	le->SetTextColor(kBlue);    
    }
    leg.Draw(); 
    
    if(useWeights)c.Print("DataVsReweightedMC/"+ Decay + "/" + finalState + "/" + selection + "/" + Branch+ "_" + anythingToString(Year) + ".eps");
    else c.Print("DataVsMC/"+Branch+".eps");
}

void dataVsMC(int Year = 11, TString finalState = "KKpi", bool useWeights=false, TString Decay = "norm", TString selection = "Preselected"){
    
    /// Options
    NamedParameter<int> nBins("nBins", 40); 

    /// Bs
    plot("Bs_P","p(B) [MeV]",nBins,0,900000,Year, finalState,useWeights, Decay, selection);
    plot("Bs_PT","p_{T}(B) [MeV]",nBins,0,40000,Year, finalState,useWeights, Decay, selection);
    plot("Bs_ETA","#eta(B)",nBins,1,6,Year, finalState,useWeights, Decay, selection);
    plot("Bs_FDCHI2_OWNPV","#chi^{2}_{FD}(B)",nBins,0,100000,Year, finalState,useWeights, Decay, selection);
    plot("Bs_ENDVERTEX_CHI2","#chi^{2}_{vtx}(B)",nBins,0,35,Year, finalState,useWeights, Decay, selection);
    plot("Bs_TAU","#tau(B) [ns]",nBins,0,16.,Year, finalState,useWeights, Decay, selection);

    /// BDT
    plot("DTF_CHI2NDOF","DTF CHI2",nBins,0.,7,Year, finalState,useWeights, Decay, selection);    
    plot("Bs_IPCHI2_OWNPV","#chi^{2}_{IP}(B)",nBins,0,16,Year, finalState,useWeights, Decay, selection);
    plot("Bs_DIRA_OWNPV","#chi^{2}_{IP}(B)",nBins,0.99997,1,Year, finalState,useWeights, Decay, selection);

    plot("XsDaughters_min_IPCHI2","X_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,Year, finalState,useWeights, Decay, selection);
    if(Decay == "norm")plot("a_1_1260_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 2.5 ,Year, finalState,useWeights, Decay, selection);
    else plot("K_1_1270_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 2.5 ,Year, finalState,useWeights, Decay, selection);
    plot("Xs_max_DOCA","X_{s} max DOCA [mm]",nBins, 0, 0.4 ,Year, finalState,useWeights, Decay, selection);

    plot("DsDaughters_min_IPCHI2","D_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,Year, finalState,useWeights, Decay, selection);
    plot("Ds_ptasy_1.00","Ds_ptasy_1.00",nBins, -1, 2.5 ,Year, finalState,useWeights, Decay, selection);
    plot("Ds_FDCHI2_ORIVX","#chi^{2}_{FD}(D_{s})",nBins,0,40000,Year, finalState,useWeights, Decay, selection);
    plot("Ds_RFD","Ds RFD",nBins,0,10,Year, finalState,useWeights, Decay, selection);

    plot("maxCos","maxCos",nBins,-1,1,Year, finalState,useWeights, Decay, selection);    
    plot("max_ghostProb","max(Track_ghostProb)",nBins,0.,0.375,Year, finalState,useWeights, Decay, selection);

    /// Ds
    plot("Ds_m12","Ds_m12",nBins,900,1900,Year, finalState,useWeights, Decay, selection);
    plot("Ds_m13","Ds_m13",nBins,600,1600,Year, finalState,useWeights, Decay, selection);
    plot("Ds_PT","p_{T}(D_{s}) [MeV]",nBins,0,40000,Year, finalState,useWeights, Decay, selection);
    plot("Ds_ETA","#eta(D_{s})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
    
    if(Decay == "signal"){
	plot("K_plus_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_ETA","#eta(K^{+})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);

	plot("pi_plus_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus_ETA","#eta(#pi^{+})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);

	plot("pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);
    }   
    
    else {
	plot("pi_plus1_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus1_ETA","#eta(K^{+})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus1_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus1_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus1_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);

	plot("pi_plus2_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus2_ETA","#eta(#pi^{+})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus2_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus2_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);
	plot("pi_plus2_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);

	plot("pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);
    }

    if(finalState == "KKpi") {
	plot("K_plus_fromDs_PT","p_{T}(K^{+} from D_{s}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_fromDs_ETA","#eta(K^{+} from D_{s})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+} from D_{s})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_fromDs_TRACK_GhostProb","ghost prob (K^{+} from D_{s})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);
	plot("K_plus_fromDs_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,Year, finalState,useWeights, Decay, selection);

	plot("pi_minus_fromDs_PT","p_{T}(#pi^{-} from D_{s}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_fromDs_ETA","#eta(#pi^{-} from D_{s})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-} from D_{s})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_fromDs_TRACK_GhostProb","ghost prob (#pi^{-} from D_{s})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);
	plot("pi_minus_fromDs_PIDK","DLL_{K#pi}(#pi^{-} from D_{s}) ",nBins,-100,100,Year, finalState,useWeights, Decay, selection);

	plot("K_minus_fromDs_PT","p_{T}(K^{-} from D_{s}) [MeV]",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_minus_fromDs_ETA","#eta(K^{-} from D_{s})",nBins,1,6,Year, finalState,useWeights, Decay, selection);
	plot("K_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-} from D_{s})",nBins,0,10000,Year, finalState,useWeights, Decay, selection);
	plot("K_minus_fromDs_TRACK_GhostProb","ghost prob (K^{-} from D_{s})",nBins,0,0.4,Year, finalState,useWeights, Decay, selection);
	plot("K_minus_fromDs_PIDK","DLL_{K#pi}(K^{-} from D_{s}) ",nBins,-20,100,Year, finalState,useWeights, Decay, selection);
    }    

    plot("NTracks","# of tracks",nBins,0,450,Year, finalState,useWeights, Decay, selection);
    if(selection == "Final") plot("BDTG_response","BDTG",nBins,0,1.2,Year, finalState,useWeights, Decay, selection);
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

void resetWeights(int Year, TString FinalState = "KKpi", TString Decay = "norm"){

    TFile* f = new TFile("/auto/data/dargent/BsDsKpipi/Preselected/MC/" + Decay + "_Ds2" + FinalState + "_" + anythingToString(Year) + ".root","UPDATE");
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
    NamedParameter<int> reweight("reweight", 1); 

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);
    
    vector<int> years;
    years.push_back(11);
    years.push_back(12);
    //years.push_back(15);
    //years.push_back(16);

    vector<TString> Ds_finalStates;
    Ds_finalStates.push_back("KKpi");
    //Ds_finalStates.push_back("pipipi");

    /// Define reweighting vars
    vector<TString>  vars_1, vars_2, vars_3, vars_4;
    vector<double>   min_1, min_2, min_3, min_4;
    vector<double>   max_1, max_2, max_3, max_4;
     
    vars_1.push_back("Bs_PT");
    min_1.push_back(0.);
    max_1.push_back(200000.);
    vars_1.push_back("Bs_ETA");
    min_1.push_back(1.5);
    max_1.push_back(5.5);
    //vars_1.push_back("NTracks");
    //min_1.push_back(0.);
    //max_1.push_back(600.);
    
    vars_2.push_back("NTracks");
    min_2.push_back(0.);
    max_2.push_back(1000.);
    vars_2.push_back("max_ghostProb");
    min_2.push_back(0.);
    max_2.push_back(0.4);

    vars_3.push_back("NTracks");
    min_3.push_back(0.);
    max_3.push_back(1000.);
    vars_3.push_back("Bs_PT");
    min_3.push_back(0.);
    max_3.push_back(200000.);

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
    if(reweight)for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++){ 
			resetWeights(years[i],Ds_finalStates[j],"norm");
			resetWeights(years[i],Ds_finalStates[j],"signal");
    }
    /// Produce MC correction histos and apply weights
    /// Weights are applied on top of each other with the previous weighting applied
    if(reweight)for(int n= 0; n < nIterations; n++)
			for(int i= 0; i < years.size(); i++) 
				for(int j= 0; j < Ds_finalStates.size(); j++)		
					for(int k =0; k < vars_set.size(); k++)	{	
						produceCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j]);
						applyCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j],"norm");    
						applyCorrectionHisto(vars_set[k],min_set[k],max_set[k],years[i],Ds_finalStates[j],"signal");    
					}
			
    /// Draw comparison plots 
    for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++){ 
		dataVsMC(years[i],Ds_finalStates[j], true,"norm","Preselected");
		//dataVsMC(years[i],Ds_finalStates[j], true,"norm","Final");
		//dataVsMC(years[i],Ds_finalStates[j], true,"signal","Final");
    }

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
