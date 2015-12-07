//main file for Bs->DsKpipi Analysis
//philippe d'argent & matthieu kecke

#include <boost/lexical_cast.hpp>
#include <cmath>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TNtuple.h>
#include "TRandom3.h"
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
#include "RooBifurGauss.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooNDKeysPdf.h"
#include "RooKeysPdf.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

using namespace std;
using namespace RooFit ;
using namespace RooStats;

int preselect(bool MC=false , bool preselection =true) {

    //specifiy decay channel of Ds
    bool Ds2KKpi = true;
    bool Ds2Kpipi = false;
    bool Ds2pipipi = false;

    //specification for normalisation channel
    bool BsDspipipi = false;

    // extended offline selection ?
    bool offline_selection = true;

    ///Load file
    TChain* tree = 0;
	
    //Ds2KKpi case
    if(Ds2KKpi){
	tree=new TChain("Bs2DsKpipi_Ds2KKpi_Tuple/DecayTree");
	tree->Add("/auto/data/kecke/B2DKPiPi/12D/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12U/*.root");
	//tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D/*.root");
	//tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/11U/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12D-PID/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/11U-MC/*.root");
    }

    //Ds2pipipi case
    if(Ds2Kpipi){
    	tree=new TChain("Bs2DsKpipi_Ds2Kpipi_Tuple/DecayTree");
	tree->Add("/auto/data/kecke/B2DKPiPi/12D/*.root");
	tree->Add("/auto/data/kecke/B2DKPiPi/12U/*.root");
    }

    //Ds2pipipi case
    if(Ds2pipipi){
    	tree=new TChain("Bs2DsKpipi_Ds2pipipi_Tuple/DecayTree");
	tree->Add("/auto/data/kecke/B2DKPiPi/12D/*.root");
	tree->Add("/auto/data/kecke/B2DKPiPi/12U/*.root");
    }

    //Bs2Dspipipi normalisation case
    if(BsDspipipi){
    tree=new TChain("Bs2Dspipipi_Ds2KKpi_Tuple/DecayTree");
  //  tree->Add("/auto/data/kecke/B2DKPiPi/12D-3pi/*.root");
   // tree->Add("/auto/data/kecke/B2DKPiPi/12U-3pi/*.root");
      tree->Add("/auto/data/kecke/B2DKPiPi/12D-3pi-PID/*.root");
      //tree->Add("/auto/data/dargent/Bs2DsKpipi/3pi/11D/*.root");
    }

    int N = tree->GetEntries();
    cout << "Old file contains " << N << " events" <<  endl;
    
    //Disable all branches
    tree->SetBranchStatus("*",0);
        
    //activate needed branches
    tree->SetBranchStatus("Bs*Decision*",1) ;

    tree->SetBranchStatus("*PID*",1) ;
    tree->SetBranchStatus("*ProbNN*",1) ;
    //tree->SetBranchStatus("*GhostProb",1) ;

	tree->SetBranchStatus("nCandidate",1) ;
	tree->SetBranchStatus("nTracks",1) ;
	tree->SetBranchStatus("nPV",1) ;
	tree->SetBranchStatus("eventNumber",1) ;
	tree->SetBranchStatus("runNumber",1) ;
	tree->SetBranchStatus("EventInSequence",1) ;
	tree->SetBranchStatus("totCandidates",1) ;
	tree->SetBranchStatus("Polarity",1) ;

	//write cone variables
	tree->SetBranchStatus("*_cp*",1);
	tree->SetBranchStatus("*deltaEta*",1);
	tree->SetBranchStatus("*deltaPhi*",1);
	tree->SetBranchStatus("*cmult*",1);
	tree->SetBranchStatus("*pxasy*",1);
	tree->SetBranchStatus("*pyasy*",1);
	tree->SetBranchStatus("*pzasy*",1);
	tree->SetBranchStatus("*ptasy*",1);
	tree->SetBranchStatus("*pasy*",1);


    tree->SetBranchStatus("*M",1) ;
	tree->SetBranchStatus("*MM",1) ;
    tree->SetBranchStatus( "*P", 1 );
	tree->SetBranchStatus( "*PX", 1 );
	tree->SetBranchStatus( "*PY", 1);
	tree->SetBranchStatus( "*PZ", 1);
	tree->SetBranchStatus( "*PE", 1);
	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*TAU", 1 );
    tree->SetBranchStatus( "*ETA", 1 );
    tree->SetBranchStatus( "*FD*", 1 );
	tree->SetBranchStatus( "*IP*", 1 );

	tree->SetBranchStatus( "*IPCHI2_OWNPV", 1 );
	tree->SetBranchStatus( "*FDCHI2_OWNPV",1 );
	tree->SetBranchStatus( "*OWNPV*", 1 );
	tree->SetBranchStatus( "*DIRA_OWNPV",1);
	tree->SetBranchStatus( "*ENDVERTEX_CHI2",1 );
        tree->SetBranchStatus( "*ENDVERTEX*",1 );
	tree->SetBranchStatus( "*DOCA*",1 );
	tree->SetBranchStatus("*TRACK_CHI2NDOF",1) ;
        tree->SetBranchStatus("Bs_L0Global_TIS",1) ;
    
    tree->SetBranchStatus( "*chi2",1 );
    tree->SetBranchStatus( "*nDOF",1 );
    tree->SetBranchStatus( "*status",1 );
    tree->SetBranchStatus( "*ctau",1 );

	if(!MC)tree->SetBranchStatus("*TRUE*",0) ;
    else tree->SetBranchStatus("*TRUE*",1) ;
    
	///Create output file
	TFile* output = 0;

	// Ds2KKpi case
	if(Ds2KKpi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_preselected_S21_PIDLine_MagDown.root","RECREATE");

	//Ds2Kpipi case
	if(Ds2Kpipi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2Kpipi_preselected.root","RECREATE");

	// Ds2pipipi case
	if(Ds2pipipi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_preselected.root","RECREATE");

	//Bs->Dspipipi normalisation case
	if(BsDspipipi)output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Bs2Dspipipi_preselected_S21PID_MagDown.root","RECREATE");
		
	///cuts
	string cuts;
	string kin_cut_pt;
        string kin_cut_p;
	string track_cuts;
	string IP_cuts;
	string Ds_cuts;
	string pid_cuts;
	string kin_cut_sumpt;

	string Bs_selection;
	string K_selection;
	
	//Bs mass window
	string mass_cuts = "Bs_MM > 4000 && Bs_MM <6000 ";

	// Trigger Selection
	string trigger_cuts_L0 = "Bs_L0HadronDecision_TOS==1 || Bs_L0Global_TIS==1";
	string trigger_cuts_HLT1 = "Bs_Hlt1TrackAllL0Decision_TOS==1";
	string trigger_cuts_HLT2 = "Bs_Hlt2Topo2BodyBBDTDecision_TOS==1 || Bs_Hlt2Topo3BodyBBDTDecision_TOS==1 || Bs_Hlt2Topo4BodyBBDTDecision_TOS==1";

	//Preselection for Bs->Ds(2KKpi)Kpipi channel 
	if(Ds2KKpi){
        	kin_cut_pt = "Ds_PT>100 && pi_minus_fromDs_PT>100 && K_minus_fromDs_PT>100 && K_plus_fromDs_PT>100 && K_plus_PT>100 && pi_plus_PT>100 && pi_minus_PT>100" ;
		kin_cut_p = "Ds_P>1000 && pi_minus_fromDs_P>1000 && K_minus_fromDs_P>1000 && K_plus_fromDs_P>1000 && K_plus_P>1000 && pi_plus_P>1000 && pi_minus_P>1000" ;
		track_cuts = "K_plus_fromDs_TRACK_CHI2NDOF < 4 && K_minus_fromDs_TRACK_CHI2NDOF < 4 && pi_minus_fromDs_TRACK_CHI2NDOF < 4 && K_plus_TRACK_CHI2NDOF < 4 && pi_plus_TRACK_CHI2NDOF < 4 && pi_minus_TRACK_CHI2NDOF < 4 && (pi_minus_PT + K_plus_PT + pi_plus_PT)>1250 && K_1_1270_plus_DOCA1<0.4 && K_1_1270_plus_DOCA2<0.4 && K_1_1270_plus_DOCA3<0.4 && abs(K_1_1270_plus_OWNPV_Z - K_1_1270_plus_ENDVERTEX_Z)>2.0 && abs(K_1_1270_plus_OWNPV_X - K_1_1270_plus_ENDVERTEX_X)>0.1 && abs(K_1_1270_plus_OWNPV_Y - K_1_1270_plus_ENDVERTEX_Y)>0.1 && (K_1_1270_plus_ENDVERTEX_CHI2/K_1_1270_plus_ENDVERTEX_NDOF)<8 && K_1_1270_plus_FDCHI2_OWNPV>16 && K_1_1270_plus_FDCHI2_ORIVX>16 && K_1_1270_plus_DIRA_OWNPV > 0.98";

		IP_cuts = "K_plus_fromDs_IPCHI2_OWNPV > 4 && K_minus_fromDs_IPCHI2_OWNPV > 4 && pi_minus_fromDs_IPCHI2_OWNPV > 4 && K_plus_IPCHI2_OWNPV > 4 && pi_plus_IPCHI2_OWNPV > 4 && pi_minus_IPCHI2_OWNPV > 4";

        	Ds_cuts = "(pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT)>1800 && abs(Ds_MM - 1968.30)<50 && (Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)<10 && Ds_FDCHI2_OWNPV>36 && Ds_FDCHI2_ORIVX>36 && Ds_DOCA1<0.5 && Ds_DOCA2<0.5 && Ds_DOCA3<0.5";
		pid_cuts = "K_minus_fromDs_PIDK>-10 && K_plus_fromDs_PIDK>-10 && K_plus_PIDK>-10 && pi_minus_fromDs_PIDK<20 && pi_plus_PIDK<20 && pi_minus_PIDK<20";
	}

	//Preselection for Bs->Ds(2KKpi)pipipi channel 
	if(BsDspipipi){
        	kin_cut_pt = "Ds_PT>100 && pi_minus_fromDs_PT>100 && K_minus_fromDs_PT>100 && K_plus_fromDs_PT>100 && pi_plus1_PT>100 && pi_plus2_PT>100 && pi_minus_PT>100" ;
		kin_cut_p = "Ds_P>1000 && pi_minus_fromDs_P>1000 && K_minus_fromDs_P>1000 && K_plus_fromDs_P>1000 && pi_plus1_P>1000 && pi_plus2_P>1000 && pi_minus_P>1000" ;
		track_cuts = "K_plus_fromDs_TRACK_CHI2NDOF < 4 && K_minus_fromDs_TRACK_CHI2NDOF < 4 && pi_minus_fromDs_TRACK_CHI2NDOF < 4 && pi_plus1_TRACK_CHI2NDOF < 4 && pi_plus2_TRACK_CHI2NDOF < 4 && pi_minus_TRACK_CHI2NDOF < 4 && (pi_minus_PT + pi_plus1_PT + pi_plus2_PT)>1250 && (pi_minus_P + pi_plus1_P + pi_plus2_P)>2000 && a_1_1260_plus_DOCA1<0.4 && a_1_1260_plus_DOCA2<0.4 && a_1_1260_plus_DOCA3<0.4 && abs(a_1_1260_plus_OWNPV_Z - a_1_1260_plus_ENDVERTEX_Z)>2.0 && abs(a_1_1260_plus_OWNPV_X - a_1_1260_plus_ENDVERTEX_X)>0.1 && abs(a_1_1260_plus_OWNPV_Y - a_1_1260_plus_ENDVERTEX_Y)>0.1 && (a_1_1260_plus_ENDVERTEX_CHI2/a_1_1260_plus_ENDVERTEX_NDOF)<8 && a_1_1260_plus_FDCHI2_OWNPV>16 && a_1_1260_plus_FDCHI2_ORIVX>16 && a_1_1260_plus_DIRA_OWNPV > 0.98";

		IP_cuts = "K_plus_fromDs_IPCHI2_OWNPV > 4 && K_minus_fromDs_IPCHI2_OWNPV > 4 && pi_minus_fromDs_IPCHI2_OWNPV > 4 && pi_plus1_IPCHI2_OWNPV > 4 && pi_plus2_IPCHI2_OWNPV > 4 && pi_minus_IPCHI2_OWNPV > 4";

        	Ds_cuts = "(pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT)>1800 && abs(Ds_MM - 1968.30)<50 && (Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)<10 && Ds_FDCHI2_OWNPV>36 && Ds_FDCHI2_ORIVX>36 && Ds_DOCA1<0.5 && Ds_DOCA2<0.5 && Ds_DOCA3<0.5";
		pid_cuts = "K_minus_fromDs_PIDK>-5 && K_plus_fromDs_PIDK>-5 && pi_plus1_PIDK<10 && pi_minus_fromDs_PIDK<10 && pi_plus2_PIDK<10 && pi_minus_PIDK<10";
	}

	//Preselection for Ds2Kpipi channel
	if(Ds2Kpipi){
        	kin_cut_pt = "Ds_PT>100 && pi_minus_fromDs_PT>100 && K_minus_fromDs_PT>100 && pi_plus_fromDs_PT>100 && K_plus_PT>100 && pi_plus_PT>100 && pi_minus_PT>100" ;
		kin_cut_p = "Ds_P>1000 && pi_minus_fromDs_P>1000 && K_minus_fromDs_P>1000 && pi_plus_fromDs_P>1000 && K_plus_P>1000 && pi_plus_P>1000 && pi_minus_P>1000" ;
		track_cuts = "pi_plus_fromDs_TRACK_CHI2NDOF < 4 && K_minus_fromDs_TRACK_CHI2NDOF < 4 && pi_minus_fromDs_TRACK_CHI2NDOF < 4 && K_plus_TRACK_CHI2NDOF < 4 && pi_plus_TRACK_CHI2NDOF < 4 && pi_minus_TRACK_CHI2NDOF < 4";
		IP_cuts = "pi_plus_fromDs_IPCHI2_OWNPV > 4 && K_minus_fromDs_IPCHI2_OWNPV > 4 && pi_minus_fromDs_IPCHI2_OWNPV > 4 && K_plus_IPCHI2_OWNPV > 4 && pi_plus_IPCHI2_OWNPV > 4 && pi_minus_IPCHI2_OWNPV > 4";

        	Ds_cuts = "(pi_minus_fromDs_PT + K_minus_fromDs_PT + pi_plus_fromDs_PT)>1800 && abs(Ds_MM - 1968.30)<50 && (Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)<8 && Ds_FDCHI2_OWNPV>16 && Ds_FDCHI2_ORIVX>16 && Ds_DIRA_OWNPV>0.98 && abs(Ds_OWNPV_Z - Ds_ENDVERTEX_Z)>2.0 && abs(Ds_OWNPV_X - Ds_ENDVERTEX_X)>0.1 && abs(Ds_OWNPV_Y - Ds_ENDVERTEX_Y)>0.1 && Ds_DOCA<0.5 && Ds_DOCA2<0.5 && Ds_DOCA3<0.5";
		pid_cuts = "K_minus_fromDs_PIDK>-5 && pi_plus_fromDs_PIDK<10 && pi_minus_fromDs_PIDK<10 && K_plus_PIDK>-10 && pi_plus_PIDK<20 && pi_minus_PIDK<20";
	}

	//Preselection for Ds2pipipi channel
	if(Ds2pipipi){
        	kin_cut_pt = "Ds_PT>100 && pi_minus_fromDs_PT>100 && pi_minus2_fromDs_PT>100 && pi_plus_fromDs_PT>100 && K_plus_PT>100 && pi_plus_PT>100 && pi_minus_PT>100" ;
		kin_cut_p = "Ds_P>1000 && pi_minus_fromDs_P>2000 && pi_minus2_fromDs_P>2000 && pi_plus_fromDs_P>2000 && K_plus_P>1000 && pi_plus_P>1000 && pi_minus_P>1000" ;
		track_cuts = "pi_plus_fromDs_TRACK_CHI2NDOF < 4 && pi_minus2_fromDs_TRACK_CHI2NDOF < 4 && pi_minus_fromDs_TRACK_CHI2NDOF < 4 && K_plus_TRACK_CHI2NDOF < 4 && pi_plus_TRACK_CHI2NDOF < 4 && pi_minus_TRACK_CHI2NDOF < 4";
		IP_cuts = "pi_plus_fromDs_IPCHI2_OWNPV > 4 && pi_minus2_fromDs_IPCHI2_OWNPV > 4 && pi_minus_fromDs_IPCHI2_OWNPV > 4 && K_plus_IPCHI2_OWNPV > 4 && pi_plus_IPCHI2_OWNPV > 4 && pi_minus_IPCHI2_OWNPV > 4";

        	Ds_cuts = "(pi_minus_fromDs_PT + pi_minus2_fromDs_PT + pi_plus_fromDs_PT)>1800 && abs(Ds_MM - 1968.30)<50 && (Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)<8 && Ds_FDCHI2_OWNPV>16 && Ds_FDCHI2_ORIVX>16 && Ds_DIRA_OWNPV>0.98 && abs(Ds_OWNPV_Z - Ds_ENDVERTEX_Z)>2.0 && abs(Ds_OWNPV_X - Ds_ENDVERTEX_X)>0.1 && abs(Ds_OWNPV_Y - Ds_ENDVERTEX_Y)>0.1 && Ds_DOCA<0.5 && Ds_DOCA2<0.5 && Ds_DOCA3<0.5";
		pid_cuts = "pi_minus2_fromDs_PIDK<10 && pi_plus_fromDs_PIDK<10 && pi_minus_fromDs_PIDK<10 && K_plus_PIDK>-10 && pi_plus_PIDK<20 && pi_minus_PIDK<20";
		kin_cut_sumpt = "(pi_minus_fromDs_PT>300 && pi_minus2_fromDs_PT>300) || (pi_minus_fromDs_PT>300 && pi_plus_fromDs_PT>300) || (pi_minus2_fromDs_PT>300 && pi_plus_fromDs_PT>300)";
	}


	//further offline selection
	if(offline_selection){
		Bs_selection = "Bs_DIRA_OWNPV>0.99994 && Bs_IPCHI2_OWNPV<20 && Bs_FDCHI2_OWNPV>100 && (Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)<40";
		K_selection = "K_plus_PIDK>8 && pi_plus_PIDK<10 && pi_minus_PIDK<10";
	}/*K_plus_PIDK>8 && */
	 

    string MCTRUE_cut = "abs(Bs_TRUEID)==531 && abs(Ds_TRUEID)==431 && abs(K_plus_TRUEID)==321 && abs(pi_plus_TRUEID)==211 && abs(pi_minus_TRUEID)==211" ;
    
    
	cuts.append(mass_cuts);
    if (preselection) {
        cuts.append("&&");
        cuts.append(kin_cut_pt);
        cuts.append("&&");
        cuts.append(kin_cut_p);
        cuts.append("&&");
        cuts.append(trigger_cuts_L0);
        cuts.append("&&");
        cuts.append(trigger_cuts_HLT1);
        cuts.append("&&");
        cuts.append(trigger_cuts_HLT2);
        cuts.append("&&");
        cuts.append(Ds_cuts);
        cuts.append("&&");
        cuts.append(track_cuts);
        cuts.append("&&");
        cuts.append(IP_cuts);
    	if(!MC){
        	cuts.append("&&");
        	cuts.append(pid_cuts);
	}
    }
    if (Ds2pipipi) {
        cuts.append("&&");
        cuts.append(kin_cut_sumpt);
    }
    if (offline_selection) {
        cuts.append("&&");
        cuts.append(Bs_selection);
	if(!MC){
        	cuts.append("&&");
        	cuts.append(K_selection);
	}
    }
    if (MC) {
        cuts.append("&&");
        cuts.append(MCTRUE_cut);
    }
    
    TTree* tree_sel = tree->CopyTree(cuts.c_str());    
    cout << "New file contains " << tree_sel->GetEntries() << " events" <<  endl;
    if(MC)cout << "Signal Eff = " << (double) ((double) tree_sel->GetEntries() / (double) tree->GetEntries()) << endl;
    
    output->WriteTObject(tree_sel);
    output->Close();
    delete output;
    return 0;
}
/*
void chooseBestPV(string fit="psiDTF", string input="/auto/data/kecke/B2DKPiPi/data2011_preselected.root"){
 
    double unit=1000000.;

    ///Read file
    TChain* tree=new TChain("DecayTree");
    tree->Add(input.c_str());
    
    ///momenta
    int nPV;
    float chi2[100], ndof[100];

    float B_M[100] ,K1_M[100] , B_TAU[100];
    float K_PX[100], K_PY[100], K_PZ[100], K_PE[100];
    float pip_PX[100], pip_PY[100], pip_PZ[100], pip_PE[100];
    float pim_PX[100], pim_PY[100], pim_PZ[100], pim_PE[100];
    float mup_PX[100], mup_PY[100], mup_PZ[100], mup_PE[100];
    float mum_PX[100], mum_PY[100], mum_PZ[100], mum_PE[100];

    tree->SetBranchAddress(("Bplus_"+fit+"_nPV").c_str(),&nPV);
    tree->SetBranchAddress(("Bplus_"+fit+"_chi2").c_str(),&chi2);
    tree->SetBranchAddress(("Bplus_"+fit+"_nDOF").c_str(),&ndof);

    tree->SetBranchAddress(("Bplus_"+fit+"_M").c_str(),&B_M);
    tree->SetBranchAddress(("Bplus_"+fit+"_ctau").c_str(),&B_TAU);

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PX").c_str(),&K_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PY").c_str(),&K_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PZ").c_str(),&K_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_Kplus_PE").c_str(),&K_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PX").c_str(),&pip_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PY").c_str(),&pip_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PZ").c_str(),&pip_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_PE").c_str(),&pip_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PX").c_str(),&pim_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PY").c_str(),&pim_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PZ").c_str(),&pim_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_K_1_1270_plus_piplus_0_PE").c_str(),&pim_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PX").c_str(),&mup_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PY").c_str(),&mup_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PZ").c_str(),&mup_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_0_PE").c_str(),&mup_PE); 

    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PX").c_str(),&mum_PX);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PY").c_str(),&mum_PY);
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PZ").c_str(),&mum_PZ); 
    tree->SetBranchAddress(("Bplus_"+fit+"_psi_2S_muminus_PE").c_str(),&mum_PE); 

    ///create new tree
    TFile* f = new TFile(("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_"+fit+".root").c_str(),"RECREATE");
    TTree* new_tree = tree->CloneTree();//CopyTree();    
    double DTF_Bplus_M, DTF_psi_M, DTF_chi2, DTF_Bplus_TAU; 
    double DTF_Bplus_PX, DTF_Bplus_PY, DTF_Bplus_PZ, DTF_Bplus_PE;      
    double DTF_Kplus_PX, DTF_Kplus_PY, DTF_Kplus_PZ, DTF_Kplus_PE;      
    double DTF_piplus_PX, DTF_piplus_PY, DTF_piplus_PZ, DTF_piplus_PE;      
    double DTF_piminus_PX, DTF_piminus_PY, DTF_piminus_PZ, DTF_piminus_PE;      
    double DTF_muplus_PX, DTF_muplus_PY, DTF_muplus_PZ, DTF_muplus_PE;      
    double DTF_muminus_PX, DTF_muminus_PY, DTF_muminus_PZ, DTF_muminus_PE; 
    double DTF_Jpsi_1S_PX, DTF_Jpsi_1S_PY, DTF_Jpsi_1S_PZ, DTF_Jpsi_1S_PE;      

    TBranch* Bra_DTF_Bplus_M = new_tree->Branch((fit+"_Bplus_M").c_str(), &DTF_Bplus_M, (fit+"_Bplus_M/D").c_str());
    TBranch* Bra_DTF_Bplus_TAU = new_tree->Branch((fit+"_TAU").c_str(), &DTF_Bplus_TAU, (fit+"_TAU/D").c_str());
    TBranch* Bra_DTF_psi_M = new_tree->Branch((fit+"_psi_M").c_str(), &DTF_psi_M,  (fit+"_psi_M/D").c_str());
    TBranch* Bra_DTF_chi2 = new_tree->Branch((fit+"_chi2").c_str(), &DTF_chi2,  (fit+"_chi2/D").c_str());

    TBranch* Bra_DTF_Bplus_PX = new_tree->Branch((fit+"_Bplus_PX").c_str(), &DTF_Bplus_PX, (fit+"_Bplus_PX/D").c_str());
    TBranch* Bra_DTF_Bplus_PY = new_tree->Branch((fit+"_Bplus_PY").c_str(), &DTF_Bplus_PY, (fit+"_Bplus_PY/D").c_str());
    TBranch* Bra_DTF_Bplus_PZ = new_tree->Branch((fit+"_Bplus_PZ").c_str(), &DTF_Bplus_PZ, (fit+"_Bplus_PZ/D").c_str());
    TBranch* Bra_DTF_Bplus_PE = new_tree->Branch((fit+"_Bplus_PE").c_str(), &DTF_Bplus_PE, (fit+"_Bplus_PE/D").c_str());

    TBranch* Bra_DTF_Kplus_PX = new_tree->Branch((fit+"_Kplus_PX").c_str(), &DTF_Kplus_PX, (fit+"_Kplus_PX/D").c_str());
    TBranch* Bra_DTF_Kplus_PY = new_tree->Branch((fit+"_Kplus_PY").c_str(), &DTF_Kplus_PY, (fit+"_Kplus_PY/D").c_str());
    TBranch* Bra_DTF_Kplus_PZ = new_tree->Branch((fit+"_Kplus_PZ").c_str(), &DTF_Kplus_PZ, (fit+"_Kplus_PZ/D").c_str());
    TBranch* Bra_DTF_Kplus_PE = new_tree->Branch((fit+"_Kplus_PE").c_str(), &DTF_Kplus_PE, (fit+"_Kplus_PE/D").c_str());

    TBranch* Bra_DTF_piplus_PX = new_tree->Branch((fit+"_piplus_PX").c_str(), &DTF_piplus_PX, (fit+"_piplus_PX/D").c_str());
    TBranch* Bra_DTF_piplus_PY = new_tree->Branch((fit+"_piplus_PY").c_str(), &DTF_piplus_PY, (fit+"_piplus_PY/D").c_str());
    TBranch* Bra_DTF_piplus_PZ = new_tree->Branch((fit+"_piplus_PZ").c_str(), &DTF_piplus_PZ, (fit+"_piplus_PZ/D").c_str());
    TBranch* Bra_DTF_piplus_PE = new_tree->Branch((fit+"_piplus_PE").c_str(), &DTF_piplus_PE, (fit+"_piplus_PE/D").c_str());

    TBranch* Bra_DTF_piminus_PX = new_tree->Branch((fit+"_piminus_PX").c_str(), &DTF_piminus_PX, (fit+"_piminus_PX/D").c_str());
    TBranch* Bra_DTF_piminus_PY = new_tree->Branch((fit+"_piminus_PY").c_str(), &DTF_piminus_PY, (fit+"_piminus_PY/D").c_str());
    TBranch* Bra_DTF_piminus_PZ = new_tree->Branch((fit+"_piminus_PZ").c_str(), &DTF_piminus_PZ, (fit+"_piminus_PZ/D").c_str());
    TBranch* Bra_DTF_piminus_PE = new_tree->Branch((fit+"_piminus_PE").c_str(), &DTF_piminus_PE, (fit+"_piminus_PE/D").c_str());

    TBranch* Bra_DTF_muplus_PX = new_tree->Branch((fit+"_muplus_PX").c_str(), &DTF_muplus_PX, (fit+"_muplus_PX/D").c_str());
    TBranch* Bra_DTF_muplus_PY = new_tree->Branch((fit+"_muplus_PY").c_str(), &DTF_muplus_PY, (fit+"_muplus_PY/D").c_str());
    TBranch* Bra_DTF_muplus_PZ = new_tree->Branch((fit+"_muplus_PZ").c_str(), &DTF_muplus_PZ, (fit+"_muplus_PZ/D").c_str());
    TBranch* Bra_DTF_muplus_PE = new_tree->Branch((fit+"_muplus_PE").c_str(), &DTF_muplus_PE, (fit+"_muplus_PE/D").c_str());

    TBranch* Bra_DTF_muminus_PX = new_tree->Branch((fit+"_muminus_PX").c_str(), &DTF_muminus_PX, (fit+"_muminus_PX/D").c_str());
    TBranch* Bra_DTF_muminus_PY = new_tree->Branch((fit+"_muminus_PY").c_str(), &DTF_muminus_PY, (fit+"_muminus_PY/D").c_str());
    TBranch* Bra_DTF_muminus_PZ = new_tree->Branch((fit+"_muminus_PZ").c_str(), &DTF_muminus_PZ, (fit+"_muminus_PZ/D").c_str());
    TBranch* Bra_DTF_muminus_PE = new_tree->Branch((fit+"_muminus_PE").c_str(), &DTF_muminus_PE, (fit+"_muminus_PE/D").c_str());

    TBranch* Bra_DTF_Jpsi_1S_PX = new_tree->Branch((fit+"_Jpsi_1S_PX").c_str(), &DTF_Jpsi_1S_PX, (fit+"_Jpsi_1S_PX/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PY = new_tree->Branch((fit+"_Jpsi_1S_PY").c_str(), &DTF_Jpsi_1S_PY, (fit+"_Jpsi_1S_PY/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PZ = new_tree->Branch((fit+"_Jpsi_1S_PZ").c_str(), &DTF_Jpsi_1S_PZ, (fit+"_Jpsi_1S_PZ/D").c_str());
    TBranch* Bra_DTF_Jpsi_1S_PE = new_tree->Branch((fit+"_Jpsi_1S_PE").c_str(), &DTF_Jpsi_1S_PE, (fit+"_Jpsi_1S_PE/D").c_str());

   ///loop over events
    int numEvents = tree->GetEntries();
    
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);

	///find best PV
	int best=0;
	double tmp=chi2[0]/ndof[0];
	if (0ul == (i % 10000ul))cout << "nPV= " << nPV << endl;
	if (0ul == (i % 10000ul))cout << "first= " << tmp << endl;

	for(int j=1; j<nPV ;j++){
                if (0ul == (i % 10000ul))cout << j << "_chi2= " << chi2[j]/ndof[j] << endl;
		if(chi2[j]/ndof[j]< tmp){
			 tmp = chi2[j]/ndof[j];
			 best=j;
		}
	}

	if (0ul == (i % 10000ul)){
		cout << "best= " << tmp << endl;
		cout << endl ;
		cout << endl ;
	}

        TLorentzVector K_P(K_PX[best],K_PY[best],K_PZ[best],K_PE[best]);
        TLorentzVector pip_P(pip_PX[best],pip_PY[best],pip_PZ[best],pip_PE[best]);
        TLorentzVector pim_P(pim_PX[best],pim_PY[best],pim_PZ[best],pim_PE[best]);
        TLorentzVector mup_P(mup_PX[best],mup_PY[best],mup_PZ[best],mup_PE[best]);
        TLorentzVector mum_P(mum_PX[best],mum_PY[best],mum_PZ[best],mum_PE[best]);

	TLorentzVector B_P = K_P + pip_P + pim_P + mup_P + mum_P;
	TLorentzVector psi_P = mup_P + mum_P;
	

	///Fill branches
	DTF_Bplus_M = B_M[best];
	if(DTF_Bplus_M<5200. || DTF_Bplus_M>5600.)continue;
	DTF_Bplus_TAU = B_TAU[best];
	DTF_chi2 = tmp;
	DTF_psi_M = psi_P.M();

	DTF_Bplus_PX = B_P.X();
	DTF_Bplus_PY = B_P.Y();
	DTF_Bplus_PZ = B_P.Z();
	DTF_Bplus_PE = B_P.E();

	DTF_Jpsi_1S_PX = psi_P.X();
	DTF_Jpsi_1S_PY = psi_P.Y();
	DTF_Jpsi_1S_PZ = psi_P.Z();
	DTF_Jpsi_1S_PE = psi_P.E();

	DTF_Kplus_PX = K_PX[best];
	DTF_Kplus_PY = K_PY[best];
	DTF_Kplus_PZ = K_PZ[best];
	DTF_Kplus_PE = K_PE[best];

	DTF_piplus_PX = pip_PX[best];
	DTF_piplus_PY = pip_PY[best];
	DTF_piplus_PZ = pip_PZ[best];
	DTF_piplus_PE = pip_PE[best];

	DTF_piminus_PX = pim_PX[best];
	DTF_piminus_PY = pim_PY[best];
	DTF_piminus_PZ = pim_PZ[best];
	DTF_piminus_PE = pim_PE[best];

	DTF_muplus_PX = mup_PX[best];
	DTF_muplus_PY = mup_PY[best];
	DTF_muplus_PZ = mup_PZ[best];
	DTF_muplus_PE = mup_PE[best];

	DTF_muminus_PX = mum_PX[best];
	DTF_muminus_PY = mum_PY[best];
	DTF_muminus_PZ = mum_PZ[best];
	DTF_muminus_PE = mum_PE[best];

        Bra_DTF_Bplus_M->Fill();
        Bra_DTF_Bplus_TAU->Fill();
	Bra_DTF_psi_M->Fill();
        Bra_DTF_chi2->Fill();

        Bra_DTF_Bplus_PX->Fill();
        Bra_DTF_Bplus_PY->Fill();
        Bra_DTF_Bplus_PZ->Fill();
        Bra_DTF_Bplus_PE->Fill();

        Bra_DTF_Jpsi_1S_PX->Fill();
        Bra_DTF_Jpsi_1S_PY->Fill();
        Bra_DTF_Jpsi_1S_PZ->Fill();
        Bra_DTF_Jpsi_1S_PE->Fill();

        Bra_DTF_Kplus_PX->Fill();
        Bra_DTF_Kplus_PY->Fill();
        Bra_DTF_Kplus_PZ->Fill();
        Bra_DTF_Kplus_PE->Fill();

        Bra_DTF_piplus_PX->Fill();
        Bra_DTF_piplus_PY->Fill();
        Bra_DTF_piplus_PZ->Fill();
        Bra_DTF_piplus_PE->Fill();

        Bra_DTF_piminus_PX->Fill();
        Bra_DTF_piminus_PY->Fill();
        Bra_DTF_piminus_PZ->Fill();
        Bra_DTF_piminus_PE->Fill();

        Bra_DTF_muplus_PX->Fill();
        Bra_DTF_muplus_PY->Fill();
        Bra_DTF_muplus_PZ->Fill();
        Bra_DTF_muplus_PE->Fill();

        Bra_DTF_muminus_PX->Fill();
        Bra_DTF_muminus_PY->Fill();
        Bra_DTF_muminus_PZ->Fill();
        Bra_DTF_muminus_PE->Fill();

     }

    new_tree->Write();
    f->Close();
}
*/
/*
void fitPreselected(){

	bool binned=false;
	bool sWeight=true;

   	 gStyle->SetOptStat(0);
    	//gStyle->SetTitleXOffset(1.1);
    	//gStyle->SetTitleYOffset(1.3);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	//gStyle->SetLabelSize(0.033,"X");
    	//gStyle->SetLabelSize(0.033,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();
	///Load file
	TFile* file;
	file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_psiDTF.root");	
	//file= new TFile("data_test.root");
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bplus_MM",1);
	tree->SetBranchStatus("Bplus_TAU",1);
	tree->SetBranchStatus("psiFit_Bplus_M",1);
   	tree->SetBranchStatus("*PT",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///B+
        RooRealVar Bplus_MM("Bplus_MM", "m(#psi K #pi #pi)", 5200., 5600.,"MeV");
	RooRealVar Bplus_PT("Bplus_PT","Bplus_PT",0.);
	RooRealVar Bplus_TAU("Bplus_TAU","Bplus_TAU",0.);
  	RooRealVar J_psi_1S_PT("J_psi_1S_PT", "J_psi_1S_PT", 0.);
  	RooRealVar muminus_PT("muminus_PT", "muminus_PT", 0.);
	RooRealVar muplus_PT("muplus_PT", "muplus_PT", 0.);
  	RooRealVar piplus_PT("piplus_PT","piplus_PT",0.);
	RooRealVar piminus_PT("piminus_PT","piminus_PT",0.);
  	RooRealVar Kplus_PT("Kplus_PT", "Kplus_PT", 0.);

	RooArgList list =  RooArgList(Bplus_MM,Bplus_PT,Bplus_TAU,J_psi_1S_PT,Kplus_PT,piplus_PT,piminus_PT,muplus_PT,muminus_PT);
	RooDataSet* data = new RooDataSet("data", "data", tree, list, "");
	 //"Bplus_PT>2000 && Kplus_PT>200 && piplus_PT>200 && piminus_PT>200 && Bplus_TAU >0.00035");
	
	RooDataSet* data_small = (RooDataSet*) data->reduce(SelectVars(RooArgSet(Bplus_MM)));
	RooDataHist* data_binned = data_small->binnedClone();

	///Define fit model
	///----------------

	///Signal model
	///-----------------------

	RooRealVar mean1("mu", "mean1", 5281.71,5150.,5350.); 
	RooRealVar sigma1("sigma_{1}", "sigma1", 15.45,10.,20.);	
	RooRealVar sigma2("sigma_{2}", "sigma2", 29.,20.,40.);
	RooRealVar scale("s", "scale", 1.924,1.,3.);
	RooFormulaVar sigma3("sigma3", "s*sigma_{1}", RooArgSet(scale,sigma1));

	//cout<<sigma3.getVal()<<endl;

	///Gaussian
	RooGaussian Gauss1("Gauss1", "Gauss1", Bplus_MM, mean1, sigma1);
	RooGaussian Gauss2("Gauss2", "Gauss2", Bplus_MM, mean1, sigma2);
	//RooGaussian Gauss3("Gauss3", "Gauss3", DeltaM, mean2, sigma3);
	//RooRealVar f_sig("f_{sig}", "signal fraction", 0.56, 0., 1.);

	/// Crystal Ball
	RooRealVar alpha("alpha","alpha",1.41,0.5,10.);
	RooRealVar n("n","n",34.4);
	RooCBShape crystal1("crystal1","CB PDF1",Bplus_MM,mean1,sigma1, alpha, n);
        RooCBShape crystal2("crystal2","CB PDF2",Bplus_MM,mean1,sigma3, alpha, n);	

	///signal pdf
	RooRealVar f_sig("f_{sig}", "signal fraction", 0.812,0.,0.95);
	//RooAddPdf signal("signal", "signal", RooArgList(Gauss1,Gauss2),RooArgList(f_sig));
	RooAddPdf signal("signal", "signal", RooArgList(crystal1,crystal2),RooArgList(f_sig));

	///Background model
	///-------------------------

	///Exponential
	RooRealVar exp_par("lambda","exp",-0.003535,-1.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bplus_MM,exp_par);
	
	///Polynomial
	RooRealVar c_0("c_{0}","c_0",0.);
	//c_0.setConstant(0.);
	RooRealVar c_1("c_{1}","c_1",-.356,-1.,0.);
	RooRealVar c_2("c_{2}","c_2",-0.01,-1.,1.);

	RooChebychev bkg_Chebychev("bkg_Chebychev","bkg_Chebychev", Bplus_MM, RooArgList(c_1));
	RooPolynomial bkg_Poly("bkg_Poly","bkg_Poly",Bplus_MM,RooArgList(c_0));

	///Gaussian
	RooRealVar mean_bkg("mu_{bkg}", "mean_bkg", 5063.,5000.,5100.); 
	RooRealVar sigma_bkg("sigma_{bkg}", "sigma_bkg", 48.5, 40., 80.);
	RooRealVar sigma2_bkg("sigma_{2}", "sigma2", 19.2, 0.1, 100.5);
	RooGaussian bkg_Gauss("bkg_Gauss", "bkg_Gauss", Bplus_MM, mean_bkg, sigma_bkg);
	
	///bkg pdf
	RooRealVar f_bkg("f_{bkg}", "background fraction", 0.11, 0., 1.);
	RooRealVar f2_bkg("f2_{bkg}", "background fraction2", 0.11, 0., 1.);

	//RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_exp, bkg_Chebychev), RooArgList(f_bkg));
	RooAddPdf bkg("bkg", "bkg", RooArgList(bkg_Gauss, bkg_Chebychev), RooArgList(f_bkg));

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(signal, bkg_Chebychev), RooArgList(n_sig,n_bkg));

	///Fit
	RooFitResult *result;
	if(binned) result = pdf->fitTo(*data_binned,Save(kTRUE),Extended());
	else result = pdf->fitTo(*data,Save(kTRUE),Extended(),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print(); 

	cout << "sigma2 =" << sigma1.getVal()*scale.getVal() << " pm " << sigma1.getError()*scale.getVal() << endl;

	///calculate # (signal)background events in signal region
	cout << endl;
	cout << endl;

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

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bplus_MM.frame();
	frame_m->SetTitle("");

	data->plotOn(frame_m,Name("data"),MarkerSize(0.1),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack));
	pdf->plotOn(frame_m,Components(signal),LineColor(kBlue),LineStyle(kDashed));
	pdf->plotOn(frame_m,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed));
	pdf->paramOn(frame_m,Layout(0.6));
	frame_m->Draw();
	c1->Print("Preselection/BmassFit.eps");

	RooPlot* frame_m2= Bplus_MM.frame();
	frame_m2->SetTitle("");
	//data->plotOn(frame_m2,MarkerSize(0.2),Binning(100),DrawOption("e"),DataError(RooAbsData::SumW2));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m2,LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m2,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m2,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	frame_m2->Draw();
	c1->Print("Preselection/BmassFit2.eps");

	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Preselection/BmassFit_log.eps");


	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bplus_MM.frame(Title("Residual Distribution")) ;
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("Preselection/residual.eps");

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bplus_MM.frame(Title("Pull Distribution")) ;
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("Preselection/pull.eps");

	if(sWeight){
		mean1.setConstant();
		sigma1.setConstant();
		sigma2.setConstant();
		f_sig.setConstant();
		f_bkg.setConstant();
		exp_par.setConstant();
		c_0.setConstant();
		alpha.setConstant();
		n.setConstant();
		//sigma3.setConstant();
		c_1.setConstant();
		scale.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bd sWeight","Bd sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("Bplus_MM,n_sig_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print("Preselection/Bd_sWeight.eps");

    		///Create output file
   		 TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CloneTree();//CopyTree();    
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("n_sig_sw", &w, "n_sig_sw/D");

  		  ///loop over events
    		  int numEvents = tree->GetEntries();
    
    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			w=sData->GetSWeight(i,"n_sig_sw");
			Bra_sw->Fill();
  		  }
   		 new_tree->Write();
   		 output->Close();
	}
}
*/
/*
void addVarsForBDT(){
	
    ///Load file
    TFile* file;
    file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight.root");	
    TTree* tree = (TTree*) file->Get("DecayTree");	

    ///Reconstructed momenta
    double K_rec[4]; 
    double pip_rec[4]; 
    double pim_rec[4]; 
    double psi_rec[4];

    tree->SetBranchAddress("Kplus_PX",&K_rec[0]);
    tree->SetBranchAddress("Kplus_PY",&K_rec[1]);
    tree->SetBranchAddress("Kplus_PZ",&K_rec[2]); 
    tree->SetBranchAddress("Kplus_PE",&K_rec[3]); 

    tree->SetBranchAddress("piplus_PX",&pip_rec[0]);
    tree->SetBranchAddress("piplus_PY",&pip_rec[1]);
    tree->SetBranchAddress("piplus_PZ",&pip_rec[2]); 
    tree->SetBranchAddress("piplus_PE",&pip_rec[3]); 

    tree->SetBranchAddress("piminus_PX",&pim_rec[0]);
    tree->SetBranchAddress("piminus_PY",&pim_rec[1]);
    tree->SetBranchAddress("piminus_PZ",&pim_rec[2]); 
    tree->SetBranchAddress("piminus_PE",&pim_rec[3]); 
    
    tree->SetBranchAddress("J_psi_1S_PX",&psi_rec[0]);
    tree->SetBranchAddress("J_psi_1S_PY",&psi_rec[1]);
    tree->SetBranchAddress("J_psi_1S_PZ",&psi_rec[2]); 
    tree->SetBranchAddress("J_psi_1S_PE",&psi_rec[3]);

    ///Refitted momenta
    double K_dtf[4]; 
    double pip_dtf[4]; 
    double pim_dtf[4]; 
    double psi_dtf[4]; 

    tree->SetBranchAddress("psiDTF_Kplus_PX",&K_dtf[0]);
    tree->SetBranchAddress("psiDTF_Kplus_PY",&K_dtf[1]);
    tree->SetBranchAddress("psiDTF_Kplus_PZ",&K_dtf[2]); 
    tree->SetBranchAddress("psiDTF_Kplus_PE",&K_dtf[3]); 

    tree->SetBranchAddress("psiDTF_piplus_PX",&pip_dtf[0]);
    tree->SetBranchAddress("psiDTF_piplus_PY",&pip_dtf[1]);
    tree->SetBranchAddress("psiDTF_piplus_PZ",&pip_dtf[2]); 
    tree->SetBranchAddress("psiDTF_piplus_PE",&pip_dtf[3]); 

    tree->SetBranchAddress("psiDTF_piminus_PX",&pim_dtf[0]);
    tree->SetBranchAddress("psiDTF_piminus_PY",&pim_dtf[1]);
    tree->SetBranchAddress("psiDTF_piminus_PZ",&pim_dtf[2]); 
    tree->SetBranchAddress("psiDTF_piminus_PE",&pim_dtf[3]); 

    tree->SetBranchAddress("psiDTF_Jpsi_1S_PX",&psi_dtf[0]);
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PY",&psi_dtf[1]);
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PZ",&psi_dtf[2]); 
    tree->SetBranchAddress("psiDTF_Jpsi_1S_PE",&psi_dtf[3]); 
   		
    ///Create output file
    TFile* output = new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_preselected_sweight_bdt.root","RECREATE");

    TTree* new_tree = tree->CloneTree();//CopyTree();    
    int sample;
    double angK,angPip, angPim;
    double mKpi, mpipi, mKpipi, mJpsipi, mJpsipipi;
    double DTF_mKpi, DTF_mpipi, DTF_mKpipi, DTF_mJpsipi, DTF_mJpsipipi;

    TBranch* Bra_s = new_tree->Branch("sample", &sample, "sample/I");
    TBranch* Bra_angK = new_tree->Branch("angK", &angK, "angK/D");
    TBranch* Bra_angPip = new_tree->Branch("angPip", &angPip, "angPip/D");
    TBranch* Bra_angPim = new_tree->Branch("angPim", &angPim, "angPim/D");
    TBranch* Bra_mKpi = new_tree->Branch("mKpi", &mKpi, "mKpi/D");
    TBranch* Bra_mpipi = new_tree->Branch("mpipi", &mpipi, "mpipi/D");
    TBranch* Bra_mKpipi = new_tree->Branch("mKpipi", &mKpipi, "mKpipi/D");
    TBranch* Bra_mJpsipi = new_tree->Branch("mJpsipi", &mJpsipi, "mJpsipi/D");
    TBranch* Bra_mJpsipipi = new_tree->Branch("mJpsipipi", &mJpsipipi, "mJpsipipi/D");
    TBranch* Bra_DTF_mKpi = new_tree->Branch("DTF_mKpi", &DTF_mKpi, "DTF_mKpi/D");
    TBranch* Bra_DTF_mpipi = new_tree->Branch("DTF_mpipi", &DTF_mpipi, "DTF_mpipi/D");
    TBranch* Bra_DTF_mKpipi = new_tree->Branch("DTF_mKpipi", &DTF_mKpipi, "DTF_mKpipi/D");
    TBranch* Bra_DTF_mJpsipi = new_tree->Branch("DTF_mJpsipi", &DTF_mJpsipi, "DTF_mJpsipi/D");
    TBranch* Bra_DTF_mJpsipipi = new_tree->Branch("DTF_mJpsipipi", &DTF_mJpsipipi, "DTF_mJpsipipi/D");

    TRandom3 r;
    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			if(r.Rndm()<0.5)sample=1;
			else sample=2;
			Bra_s->Fill();
			
			TVector3 psi(psi_rec[0],psi_rec[1],0.);
			TVector3 K(K_rec[0],K_rec[1],0.);
			TVector3 pip(pip_rec[0],pip_rec[1],0.);
			TVector3 pim(pim_rec[0],pim_rec[1],0.);
			angK= psi.Angle(K);
			angPip= psi.Angle(pip);
			angPim= psi.Angle(pim);
			Bra_angK->Fill();
			Bra_angPip->Fill();
			Bra_angPim->Fill();

    		        TLorentzVector K_recP(K_rec[0],K_rec[1],K_rec[2],K_rec[3]);
        		TLorentzVector pip_recP(pip_rec[0],pip_rec[1],pip_rec[2],pip_rec[3]);
			TLorentzVector pim_recP(pim_rec[0],pim_rec[1],pim_rec[2],pim_rec[3]);
        		TLorentzVector psi_recP(psi_rec[0],psi_rec[1],psi_rec[2],psi_rec[3]);

			mKpi= (K_recP+pim_recP).M2();
			mpipi= (pip_recP+pim_recP).M2();
			mKpipi= (K_recP+pim_recP+pip_recP).M2();
			mJpsipi= (pip_recP+psi_recP).M2();
			mJpsipipi= (psi_recP+pim_recP+pip_recP).M2();

			Bra_mKpi->Fill();
			Bra_mpipi->Fill();
			Bra_mKpipi->Fill();
			Bra_mJpsipi->Fill();
			Bra_mJpsipipi->Fill();

        		TLorentzVector K_dtfP(K_dtf[0],K_dtf[1],K_dtf[2],K_dtf[3]);
       			TLorentzVector pip_dtfP(pip_dtf[0],pip_dtf[1],pip_dtf[2],pip_dtf[3]);
			TLorentzVector pim_dtfP(pim_dtf[0],pim_dtf[1],pim_dtf[2],pim_dtf[3]);
        		TLorentzVector psi_dtfP(psi_dtf[0],psi_dtf[1],psi_dtf[2],psi_dtf[3]);

			DTF_mKpi= (K_dtfP+pim_dtfP).M2();
			DTF_mpipi= (pip_dtfP+pim_dtfP).M2();
			DTF_mKpipi= (K_dtfP+pim_dtfP+pip_dtfP).M2();
			DTF_mJpsipi= (pip_dtfP+psi_dtfP).M2();
			DTF_mJpsipipi= (psi_dtfP+pim_dtfP+pip_dtfP).M2();

			Bra_DTF_mKpi->Fill();
			Bra_DTF_mpipi->Fill();
			Bra_DTF_mKpipi->Fill();
			Bra_DTF_mJpsipi->Fill();
			Bra_DTF_mJpsipipi->Fill();
			
			
     }
     new_tree->Write();
     output->Close();	
}
*/

void addCut(){

TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_with_BDT_variables_S21_PID.root");
TTree* tree = (TTree*) file->Get("DecayTree");


   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_with_BDT_variables_S21_PID_P_ETA.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree("pi_plus1_ETA > 1.5 && pi_plus2_ETA > 1.5 && pi_minus_ETA > 1.5 && pi_plus1_ETA < 5.0 && pi_plus2_ETA < 5.0 && pi_minus_ETA < 5.0 && pi_plus1_P < 150000 && pi_plus2_P < 150000 && pi_minus_P < 150000");
   new_tree_BDTG->Write();
   output_BDTG->Close();

   //close file at the end
   file->Close();

}


void applyBDTcut(string cutoff){
   ///Load file
   //TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/Bs2DsKpipi_Ds2KKpi_BDTtrained_S21PID.root");
   //TTree* tree = (TTree*) file->Get("DecayTree");	
TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_BDTtrained_S21PID.root");
TTree* tree = (TTree*) file->Get("DecayTree");

//string cutstring = "Ds_MM > 1950 && Ds_MM < 1990 && BDTG_response > ";
string cutstring = "BDTG_response > ";
cutstring.append(cutoff);
const char* cstringcutstring = cutstring.c_str();

   //BDT application
   //TFile* output_BDT=new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDT.root","RECREATE");
   //TTree* new_tree_BDT = tree->CopyTree("BDT_response>0.0521");
   //new_tree_BDT->Write();
   //output_BDT->Close();


    //BDTG application for 2011 Signal Channel, 0.8300 for maximum S/sqrt(S+B) with S = 286 and B = 7668
/*
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2011/Bs2DsKpipi_Ds2pipipi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

    //BDTG application for 2012 Signal Channel (9500BG/450S), 0.8040 for maximum S/sqrt(S+B) , with S = 726 , B = 18751 
/*
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/Bs2DsKpipi_Ds2KKpi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

    //BDTG application for 2012 normalization Channel (9500BG/450S) , 0.1458 for maximum S/sqrt(S+B) with S = 13867 , B = 23588
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();

/*
    //BDTG application for 2011 normalization Channel (32.500BG/6500S) ------> BDTG > 0.2084 with S=5455 and B=8827
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

   //close file at the end
   file->Close();
}

double *fitBGShape(double fitValues[11]){

//define shape of Bs->Ds*Kpipi BG as 3 bifurcated gaussians
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 5000., 5350.,"MeV/c^{2}");

//mean of gaussians
RooRealVar mean1("mean1","mu", 5059.,5040.,5070.);
RooRealVar mean2("mean2","mu", 5182.,5165.,5205.);
RooRealVar mean3("mean3","mu", 5285.,5270.,5300.);

//width of gaussians
RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", 25.9,15.,40.);
RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", 99.4,80.,115.);
RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", 13.1,5.,25.);
RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", 49.5,25.,70.);
RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", 107.,90.,125.);
RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", 21.1,5.,33.);

//bifurcated gaussians
RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);

//fractions of gauss functions
RooRealVar f_1("f_{1}", "fraction1", 0.405, 0.2, 0.6);
RooRealVar f_2("f_{2}", "fraction2", 0.1, 0., 0.3);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2));

//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/DsstKpipi.root");
TTree* tree = (TTree*) file->Get("DecayTree");
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_MM",1);

//Fill needed variable in RooDataSet
RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

//try RooKernel
RooKeysPdf kest("kest","kest",Bs_MM,*data,RooKeysPdf::NoMirror);

//Fit
RooFitResult *result;
result = pdf->fitTo(*data,Save(kTRUE),NumCPU(3));
cout << "result is --------------- "<<endl;
result->Print(); 

//fill array to return the fit values
fitValues[0] = mean1.getVal();
fitValues[1] = mean2.getVal();
fitValues[2] = mean3.getVal();
fitValues[3] = sigmaL1.getVal();
fitValues[4] = sigmaR1.getVal();
fitValues[5] = sigmaL2.getVal();
fitValues[6] = sigmaR2.getVal();
fitValues[7] = sigmaL3.getVal();
fitValues[8] = sigmaR3.getVal();
fitValues[9] = f_1.getVal();
fitValues[10] = f_2.getVal();

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
c1->Print("eps/BkgShape/Bs2DsstartKpipi.eps");

return fitValues;
}

double *fitBGShapeNorm(double fitValues[11]){

//define shape of Bs->Ds*Kpipi BG as 3 bifurcated gaussians
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
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/Dsstpipipi.root");
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
fitValues[0] = mean1.getVal();
fitValues[1] = mean2.getVal();
fitValues[2] = mean3.getVal();
fitValues[3] = sigmaL1.getVal();
fitValues[4] = sigmaR1.getVal();
fitValues[5] = sigmaL2.getVal();
fitValues[6] = sigmaR2.getVal();
fitValues[7] = sigmaL3.getVal();
fitValues[8] = sigmaR3.getVal();
fitValues[9] = f_1.getVal();
fitValues[10] = f_2.getVal();

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

return fitValues;
}


double *fitBGShapeNormDstKpipi(double fitValues[9]){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 4800., 5300.,"MeV/c^{2}");

//mean of crrystal balls
RooRealVar mean1("mean1","mu", 5050.,4900.,5100.);
RooRealVar mean2("mean2","mu", 5200.,5100.,5300.);

// asymmetry parameter of crystsal balls
RooRealVar a1("a1","a1",5., 0.,10.);
RooRealVar a2("a2","a2",8.47);
RooRealVar n1("n1","n1",0.3, 0.,5.);
RooRealVar n2("n2","n2",2.6);

//sigma of crystal balls
RooRealVar sigma1("sigma_{1}", "sigma1", 24.,5.,300.);
RooRealVar sigma2("sigma_{2}", "sigma2", 89.,5.,300.);

//crystal Balls
RooCBShape CB1("CB1", "CB1", Bs_MM, mean1, sigma1, a1, n1);
RooCBShape CB2("CB2", "CB2", Bs_MM, mean2, sigma2, a2, n2);

//fraction of crystal balls
RooRealVar f_1("f_{1}", "fraction1", 0.5, 0., 1.);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1, CB2), RooArgList(f_1));


//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsstKpipi.root");
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
fitValues[0] = mean1.getVal();
fitValues[1] = mean2.getVal();
fitValues[2] = a1.getVal();
fitValues[3] = a2.getVal();
fitValues[4] = n1.getVal();
fitValues[5] = n2.getVal();
fitValues[6] = sigma1.getVal();
fitValues[7] = sigma2.getVal();
fitValues[8] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_MM.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2DsstarKpipi_as_Dspipipi.eps");

return fitValues;
}


double *fitBGShapeNormKpipi(double fitValues[6]){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 4800., 5800.,"MeV/c^{2}");

//mean of crrystal balls
RooRealVar mean1("mean1","mu", 5300.,5200.,5400.);

// asymmetry parameter of crystsal balls
RooRealVar a1("a1","a1",5., 0.,10.);
RooRealVar n1("n1","n1",0.3, 0.,50.);

//sigma of crystal balls
RooRealVar sigma1("sigma_{1}", "sigma1", 24.,5.,120.);

//exp background
RooRealVar alpha("alpha", "alpha", 0.,-10.,10.);
RooExponential expBkg("expBkg", "expBkg", Bs_MM, alpha);

//crystal Balls
RooCBShape CB1("CB1", "CB1", Bs_MM, mean1, sigma1, a1, n1);

//fraction of crystal balls
RooRealVar f_1("f_{1}", "fraction1", 0.5, 0., 1.);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1,expBkg), RooArgList(f_1));


//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsKpipi.root");
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
fitValues[0] = mean1.getVal();
fitValues[1] = a1.getVal();
fitValues[2] = n1.getVal();
fitValues[3] = alpha.getVal();
fitValues[4] = sigma1.getVal();
fitValues[5] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_MM.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(expBkg),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2DsKpipi_as_Dspipipi.eps");

return fitValues;
}

double *fitBGShapethreePi(double fitValues[9]){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 5200., 6000.,"MeV/c^{2}");

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
RooCBShape CB1("CB1", "CB1", Bs_MM, mean1, sigma1, a1, n1);
RooCBShape CB2("CB2", "CB2", Bs_MM, mean2, sigma2, a2, n2);

//fraction of crystal balls
RooRealVar f_1("f_{1}", "fraction1", 0.5, 0., 1.);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1, CB2), RooArgList(f_1));


//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dspipipi.root");
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
fitValues[0] = mean1.getVal();
fitValues[1] = mean2.getVal();
fitValues[2] = a1.getVal();
fitValues[3] = a2.getVal();
fitValues[4] = n1.getVal();
fitValues[5] = n2.getVal();
fitValues[6] = sigma1.getVal();
fitValues[7] = sigma2.getVal();
fitValues[8] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_MM.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dspipipi_as_DsKpipi.eps");

return fitValues;
}


double *fitBGShapethreePiDstar(double fitValues[9]){

//define shape of Bs->Ds(*)pipipi BG as 2 crystal balls
//observable
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 4980., 6000.,"MeV/c^{2}");

//mean of crrystal balls
RooRealVar mean1("mean1","mu", 5350.,5200.,5490.);
RooRealVar mean2("mean2","mu", 5517.,5300.,5650.);

// asymmetry parameter of crystsal balls
RooRealVar a1("a1","a1",-1.5, -4.,3.);
RooRealVar a2("a2","a2",-0.6, -3.5,1.5);
RooRealVar n1("n1","n1",0.3, 0.,5.);
RooRealVar n2("n2","n2",10., 0.,200.);


//sigma of crystal balls
RooRealVar sigma1("sigma_{1}", "sigma1", 24.,5.,300.);
RooRealVar sigma2("sigma_{2}", "sigma2", 89.,5.,300.);

//crystal Balls
RooCBShape CB1("CB1", "CB1", Bs_MM, mean1, sigma1, a1, n1);
RooCBShape CB2("CB2", "CB2", Bs_MM, mean2, sigma2, a2, n2);

//fraction of crystal balls
RooRealVar f_1("f_{1}", "fraction1", 0.5, 0., 1.);

//add all gaussians
RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(CB1, CB2), RooArgList(f_1));


//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dsstpipipi.root");
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
fitValues[0] = mean1.getVal();
fitValues[1] = mean2.getVal();
fitValues[2] = a1.getVal();
fitValues[3] = a2.getVal();
fitValues[4] = n1.getVal();
fitValues[5] = n2.getVal();
fitValues[6] = sigma1.getVal();
fitValues[7] = sigma2.getVal();
fitValues[8] = f_1.getVal();

//plot mass distribution and fit results
TCanvas* c1= new TCanvas("");
RooPlot* frame_m= Bs_MM.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.eps");

return fitValues;
}


/*
void chooseMultCand(){
   ///Load file
   TFile* file= new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final_multCand_PIDK.root");
   TTree* tree = (TTree*) file->Get("DecayTree");	
   ///Create output
   TFile* output=new TFile("/auto/data/dargent/Bu2psiKpipi/data/data_final_chosenCand_PIDK_new.root","RECREATE");
   TTree* new_tree = tree->CopyTree("isSelectedMultipleCand==1 && Bplus_MM<5500");
   new_tree->Write();
   file->Close();
   output->Close();	
}*/

void fitBDT(){

	bool sWeight=false;

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
	
	///Load file
	TFile* file;
	file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/Bs2DsKpipi_Ds2KKpi_fullSelectionBDTG.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_MM",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar Bs_MM("Bs_MM", "m(D_{s} K #pi #pi)", 4800., 5800.,"MeV");
	RooArgList list =  RooArgList(Bs_MM);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

	///import Bkg shapes from MC fits--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//1) Ds*Kpipi in DsKpipi
	// 11 params: mean1, mean2, mean3, sigmaL1, sigmaR1, sigmaL2, sigmaR2, sigmaL3, sigmaR3, f_1, f_2
/*	double fillarr1[11];
	double *DstarKpipiSig;
	DstarKpipiSig = fitBGShape(fillarr1);
*/
/*
	//2) Ds*pipipi in Dspipipi
	// 11 params: mean1, mean2, mean3, sigmaL1, sigmaR1, sigmaL2, sigmaR2, sigmaL3, sigmaR3, f_1, f_2
	double fillarr2[11];
	double *DstarpipipiNorm;
	DstarpipipiNorm = fitBGShapeNorm(fillarr2);

	//3) Ds*Kpipi in Dspipipi
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
	double fillarr3[9];
	double *DstarKpipiNorm;
	DstarKpipiNorm = fitBGShapeNormDstKpipi(fillarr3);

	//4) DsKpipi in Dspipipi
	// 6 params: mean1 a1, n1, alpha, sigma1, f_1
	double fillarr4[6];
	double *DsKpipiNorm;
	DsKpipiNorm = fitBGShapeNormKpipi(fillarr4);
*/
	//5) Dspipipi in DsKpipi
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
	double fillarr5[9];
	double *DspipipiSig;
	DspipipiSig = fitBGShapethreePi(fillarr5);

	//6) Ds*pipipi in DsKpipi
	// 9 params: mean1, mean2, a1, a2, n1, n2, sigma1, sigma2, f_1
	double fillarr6[9];
	double *DstarpipipiSig;
	DstarpipipiSig = fitBGShapethreePiDstar(fillarr6);

	///Define fit model----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///Signal model - two double gaussians for B0 and Bs signal 

	//B0 signal shape
	RooRealVar meanB01("meanB01", "B_{0} #mu", 5281.71,5250.,5320.); 
	RooRealVar sigmaB1("sigmaB1", "B_{0} #sigma_{1}",  32.56);	
	RooRealVar sigmaB2("sigmaB2", "B_{0} sigma_{2}", 13.,10.,20.);
	RooGaussian GaussB01("GaussB01", "GaussB01", Bs_MM, meanB01, sigmaB1);
	RooGaussian GaussB02("GaussB02", "GaussB02", Bs_MM, meanB01, sigmaB2);
	RooRealVar f_GaussB("f_GaussB" , "f__{B_{0}}", 0.066);//, 0.2, 0.75);
        RooAddPdf DoubleGaussB0("DoubleGaussB0", "DoubleGaussB0", RooArgList(GaussB01,GaussB02),RooArgList(f_GaussB));

	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	//RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 15.45,10.,40.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 15.,10.,55.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_MM, meanBs1, sigmaB1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_MM, meanBs1, sigmaBs2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0.2, 0.75);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussB));

	//signal pdf
	RooRealVar N_B0("N_B0", "#B_{0}", 500, 0, 5000);
	RooRealVar N_Bs("N_Bs", "#B_{s}", 500, 0, 5000);
//	RooAddPdf signal("signal", "signal", RooArgList(DoubleGaussB0, DoubleGaussBs),RooArgList(N_B0, N_Bs));	

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential
	RooRealVar exp_par("exp_par","#lambda",0.,-10.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bs_MM,exp_par);


	///assign imported shapes

	//1) Ds*Kpipi in DsKpipi-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
	double mean1Variation_low = DstarKpipiSig[0] - (0.01 * DstarKpipiSig[0]);
	double mean1Variation_high = DstarKpipiSig[0] + (0.01 * DstarKpipiSig[0]);
	double mean2Variation_low = DstarKpipiSig[1] - (0.01 * DstarKpipiSig[1]);
	double mean2Variation_high = DstarKpipiSig[1] + (0.01 * DstarKpipiSig[1]);
	double mean3Variation_low = DstarKpipiSig[2] - (0.11 * DstarKpipiSig[2]);
	double mean3Variation_high = DstarKpipiSig[2] + (0.01 * DstarKpipiSig[2]);


	//mean of gaussians
	RooRealVar mean1("mean1","mu", DstarKpipiSig[0], mean1Variation_low, mean1Variation_high);
	RooRealVar mean2("mean2","mu", DstarKpipiSig[1], mean2Variation_low, mean2Variation_high);
	RooRealVar mean3("mean3","mu", DstarKpipiSig[2], mean3Variation_low, mean3Variation_high);

	//width of gaussians
	RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", DstarKpipiSig[3]);
	RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", DstarKpipiSig[4]);
	RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", DstarKpipiSig[5]);
	RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", DstarKpipiSig[6]);
	RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", DstarKpipiSig[7]);
	RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", DstarKpipiSig[8]);

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);

	//fractions of gauss functions
	RooRealVar f_1("f_{1}", "fraction1", DstarKpipiSig[9] );
	RooRealVar f_2("f_{2}", "fraction2", DstarKpipiSig[10] );

	//add functions
	RooAddPdf DstarKpipi_as_DsKpipi("DstarKpipi_as_DsKpipi", "DstarKpipi_as_DsKpipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2));

*/
	//fix shape from dspipipi normalization fit

	//mean of gaussians for 2011 data
	/*
	RooRealVar mean1("mean1","mu", 4.90193e+03 );
	RooRealVar mean2("mean2","mu", 5.17709e+03 );
	RooRealVar mean3("mean3","mu", 5.27608e+03 );
	*/

	//mean of gaussians for 2012 data
	RooRealVar mean1("mean1","mu", 4.9082e+03 );
	RooRealVar mean2("mean2","mu", 5.1797e+03 );
	RooRealVar mean3("mean3","mu", 5.2767e+03 );

	//width of gaussians for 2011 data
	/*
	RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", 4.45211e+01 );
	RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", 9.93543e+01 );
	RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", 1.07812e+02 );
	RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", 5.03657e+01 );
	RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", 1.04149e+02 );
	RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", 3.47862e+01 );
	*/
	//width of gaussians for 2012 data
	
	RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", 4.4246e+01 );
	RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", 1.4536e+02 );
	RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", 9.3027e+01 );
	RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", 4.9408e+01 );
	RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", 5.1916e+01 );
	RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", 2.9898e+01 );

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);

	//fractions of gauss functions for 2011 data
	/*
	RooRealVar f_1("f_{1}", "fraction1", 1.54511e-01 );
	RooRealVar f_2("f_{2}", "fraction2", 5.72217e-01 );
	*/

	//fractions of gauss functions for 2012 data

	RooRealVar f_1("f_{1}", "fraction1", 2.0851e-01 );
	RooRealVar f_2("f_{2}", "fraction2", 6.2242e-01 );

	//add functions
	RooAddPdf DstarKpipi_as_DsKpipi("DstarKpipi_as_DsKpipi", "DstarKpipi_as_DsKpipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2));


	///same shape shifted down for m(B^0)
	///mean of gaussians shifted by m_Bs - m_B0
	RooRealVar mean1Shifted("mean1Shifted","mu", mean1.getVal() - 87.33 );
	RooRealVar mean2Shifted("mean2Shifted","mu", mean2.getVal() - 87.33 );
	RooRealVar mean3Shifted("mean3Shifted","mu", mean3.getVal() - 87.33 );

	RooBifurGauss BifGauss1Shifted("BifGauss1Shifted","BifGauss1Shifted", Bs_MM, mean1Shifted, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2Shifted("BifGauss2Shifted","BifGauss2Shifted", Bs_MM, mean2Shifted, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3Shifted("BifGauss3Shifted","BifGauss3Shifted", Bs_MM, mean3Shifted, sigmaL3,sigmaR3);

	RooAddPdf DstarKpipi_as_DsKpipi_Shifted("DstarKpipi_as_DsKpipi_Shifted", "DstarKpipi_as_DsKpipi_Shifted", RooArgList(BifGauss1Shifted, BifGauss2Shifted, BifGauss3Shifted), RooArgList(f_1,f_2));


	//2) Dspipipi in DsKpipi--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	double meanDspipipi1Variation_low = DspipipiSig[0] - (0.00001 * DspipipiSig[0]);
	double meanDspipipi1Variation_high = DspipipiSig[0] + (0.01 * DspipipiSig[0]);
	double meanDspipipi2Variation_low = DspipipiSig[1] - (0.0001 * DspipipiSig[1]);
	double meanDspipipi2Variation_high = DspipipiSig[1] + (0.01 * DspipipiSig[1]);


	//mean of crrystal balls
	RooRealVar meanDspipipi1("meanDspipipi1","mu",5480. /*DspipipiSig[0]*/, 5450. /*meanDspipipi1Variation_low*/, 5550./*meanDspipipi1Variation_high*/);
	RooRealVar meanDspipipi2("meanDspipipi2","mu",5480. /*DspipipiSig[1]*/, 5450. /*meanDspipipi2Variation_low*/, 5550./*meanDspipipi2Variation_high*/);

	// asymmetry parameter of crystsal balls
	RooRealVar a1("a1","a1", DspipipiSig[2]);
	RooRealVar a2("a2","a2", DspipipiSig[3]);
	RooRealVar n1("n1","n1", DspipipiSig[4]);
	RooRealVar n2("n2","n2", DspipipiSig[5]);

	//sigma of crystal balls
	RooRealVar sigma1("sigma_{1}", "sigma1", DspipipiSig[6]);
	RooRealVar sigma2("sigma_{2}", "sigma2", DspipipiSig[7]);

	//crystal Balls
	RooCBShape CB1("CB1", "CB1", Bs_MM, meanDspipipi1, sigma1, a1, n1);
	RooCBShape CB2("CB2", "CB2", Bs_MM, meanDspipipi2, sigma2, a2, n2);

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
	RooCBShape CBDstarThreePi1("CBDstarThreePi1", "CB1", Bs_MM, meanDstarThreePi1, sigmaDstarThreePi1, aDstarThreePi1, nDstarThreePi1);
	RooCBShape CBDstarThreePi2("CBDstarThreePi2", "CB2", Bs_MM, meanDstarThreePi2, sigmaDstarThreePi2, aDstarThreePi2, nDstarThreePi2);

	//fraction of crystal balls
	RooRealVar f_4("f_{4}", "fraction4", DstarpipipiSig[8]);

	//add all gaussians
	RooAddPdf Dstarpipipi_as_DsKpipi("Dstarpipipi_as_DsKpipi", "Dstarpipipi_as_DsKpipi", RooArgList(CBDstarThreePi1, CBDstarThreePi2), RooArgList(f_4));


	//sum all background shapes
	//yields
	RooRealVar N_DstarKpipi("N_DstarKpipi","N_DstarKpipi", 1275, 0, data->numEntries());
	RooRealVar N_DstarKpipiShifted("N_DstarKpipiShifted","N_DstarKpipiShifted", 1275, 0, data->numEntries());
	/// 2011 bg yield
	//RooRealVar N_Dspipipi("N_Dspipipi","N_Dspipipi", 228, 0., data->numEntries());
	//RooRealVar N_Dstarpipipi("N_Dstarpipipi","N_Dstarpipipi", 444);//, 0., data->numEntries());
	/// 2012 bg yield
	RooRealVar N_Dspipipi("N_Dspipipi","N_Dspipipi", 673, 0., data->numEntries());
	RooRealVar N_Dstarpipipi("N_Dstarpipipi","N_Dstarpipipi", 1297);//, 0., data->numEntries());

	RooAddPdf bkg_model("bkg_model", "bkg_model", RooArgList(DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi), RooArgList(N_DstarKpipi, N_Dspipipi, N_Dstarpipipi));

	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///total pdf
	///----------------------
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2. , 0., data->numEntries());
	RooRealVar N_comb("N_comb","N_comb",100, 0, data->numEntries());

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0, N_Bs, N_DstarKpipi, N_Dspipipi, N_Dstarpipipi, N_comb, N_DstarKpipiShifted));

	///Fit
	RooFitResult *result;
	result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
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
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");


	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(60));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
	//pdf->plotOn(frame_m,Components(signal),LineColor(kBlue),LineWidth(1));
	pdf->plotOn(frame_m,Components(DoubleGaussB0),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(DoubleGaussBs),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	//pdf->plotOn(frame_m,Components(bkg_model),LineColor(kMagenta),LineWidth(1));
	pdf->plotOn(frame_m,Components(DstarKpipi_as_DsKpipi),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(Dspipipi_as_DsKpipi),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(Dstarpipipi_as_DsKpipi),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(bkg_exp),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(DstarKpipi_as_DsKpipi_Shifted),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
	//pdf->paramOn(frame_m,Layout(0.6));
	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
	frame_m->Draw();
	c1->Print("eps/Final/BmassFit_12.eps");
/*
	RooPlot* frame_m2= Bplus_MM.frame();
	frame_m2->SetTitle("");
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m2,LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m2,Components(signal),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m2,Components(bkg_Chebychev),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	data->plotOn(frame_m2,MarkerSize(0.5),Binning(50));
	frame_m2->Draw();
	c1->Print("Final/BmassFit2.eps");

	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Final/BmassFit_log.eps");
*/

	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bs_MM.frame(Title("Residual Distribution")) ;
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("eps/Final/residual_12.eps");

	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bs_MM.frame(Title("Pull Distribution")) ;
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("eps/Final/pull_12.eps");

	if(sWeight){

		meanB01.setConstant();
		meanBs1.setConstant();
		sigmaB2.setConstant();
		sigmaBs2.setConstant();
		exp_par.setConstant();
		meanDspipipi1.setConstant();
		meanDspipipi2.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(N_B0, N_Bs, N_DstarKpipi, N_Dspipipi, N_Dstarpipipi, N_comb, N_DstarKpipiShifted)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("Bs_MM,N_Bs_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print("eps/Final/Bs_12_sWeight.eps");

    		///Create output file
   		 TFile* output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root","RECREATE");
		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("Bs_MM > 4800 && Bs_MM < 5800");
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("N_Bs_sw", &w, "N_Bs_sw/D");

  		  ///loop over events
    		  int numEvents = new_tree->GetEntries();

    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);
			w=sData->GetSWeight(i,"N_Bs_sw");
			Bra_sw->Fill();
  		  }
	cout << "loop finished!!!" << endl;
   		 new_tree->Write();
   		 output->Close();
	}
}


double fitBDTNorm(){

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

	bool BDTscan = false;

	///Load file
	TFile* file;
	file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG.root");
	//file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Bs2Dspipipi_with_BDT_variables_S21_PID.root");
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_MM",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar Bs_MM("Bs_MM", "m(D_{s} #pi #pi #pi)", 4800., 5800.,"MeV/c^{2}");
  	//RooRealVar Bs_MM("Bs_MM", "m(D_{s} #pi #pi #pi)", 5300., 5450.,"MeV/c^{2}");
	RooArgList list =  RooArgList(Bs_MM);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

	///import Bkg shapes from MC fits--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//2) Ds*pipipi in Dspipipi
	// 11 params: mean1, mean2, mean3, sigmaL1, sigmaR1, sigmaL2, sigmaR2, sigmaL3, sigmaR3, f_1, f_2
	double fillarr2[11];
	double *DstarpipipiNorm;
	DstarpipipiNorm = fitBGShapeNorm(fillarr2);

	///Define fit model----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///Signal model - two double gaussians for B0 and Bs signal 


/*
	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5366.7,5345.,5380.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 32.56);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 13.54,10.,20.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_MM, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_MM, meanBs1, sigmaBs2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.066);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
*/

	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5366.7,5345.,5380.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 32.56);//, 5.,60.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 13.54,10.,20.);
	RooRealVar a1("a1","a1", 1.83);
	RooRealVar n1("n1","n1", 1.48);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_MM, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_MM, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", Bs_MM, meanBs1, sigmaBs2, a1, n1);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.18);//, 0.,0.5);
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_GaussBs));
	RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));

	//signal pdf
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential
	RooRealVar exp_par("exp_par","#lambda",-1.6508e-03,-10.,0.);	
	RooExponential bkg_exp("bkg_exp","exponential background",Bs_MM,exp_par);


	///assign imported shapes

	//1) Ds*pipipi in Dspipipi-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	double mean1Variation_low = DstarpipipiNorm[0] - (0.05 * DstarpipipiNorm[0]);
	double mean1Variation_high = DstarpipipiNorm[0] + (0.05 * DstarpipipiNorm[0]);
	double mean2Variation_low = DstarpipipiNorm[1] - (0.05 * DstarpipipiNorm[1]);
	double mean2Variation_high = DstarpipipiNorm[1] + (0.05 * DstarpipipiNorm[1]);
	double mean3Variation_low = DstarpipipiNorm[2] - (0.01 * DstarpipipiNorm[2]);
	double mean3Variation_high = DstarpipipiNorm[2] + (0.01 * DstarpipipiNorm[2]);

	double sigmaL1Variation_low = DstarpipipiNorm[3] - (/*0.75 **/ DstarpipipiNorm[3]);
	double sigmaL1Variation_high = DstarpipipiNorm[3] + (/*0.75 **/ DstarpipipiNorm[3]);
	double sigmaR1Variation_low = DstarpipipiNorm[4] - (/*0.75 **/ DstarpipipiNorm[4]);
	double sigmaR1Variation_high = DstarpipipiNorm[4] + (/*0.80 **/ DstarpipipiNorm[4]);
	double sigmaL2Variation_low = DstarpipipiNorm[5] - (/*0.40 **/ DstarpipipiNorm[5]);
	double sigmaL2Variation_high = DstarpipipiNorm[5] + (/*0.40 **/ DstarpipipiNorm[5]);
	double sigmaR2Variation_low = DstarpipipiNorm[6] - (/*0.40 **/ DstarpipipiNorm[6]);
	double sigmaR2Variation_high = DstarpipipiNorm[6] + (/*0.40 **/ DstarpipipiNorm[6]);
	double sigmaL3Variation_low = DstarpipipiNorm[7] - (/*0.40 **/ DstarpipipiNorm[7]);
	double sigmaL3Variation_high = DstarpipipiNorm[7] + (/*0.40 **/ DstarpipipiNorm[7]);
	double sigmaR3Variation_low = DstarpipipiNorm[8] - (/*0.1 **/ DstarpipipiNorm[8]);
	double sigmaR3Variation_high = DstarpipipiNorm[8] + (/*0.1 **/ DstarpipipiNorm[8]);


	//mean of gaussians
	RooRealVar mean1("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
	RooRealVar mean2("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
	RooRealVar mean3("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);

	//width of gaussians
	RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", /*4.4521e+01*/ DstarpipipiNorm[3], sigmaL1Variation_low, sigmaL1Variation_high);
	RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", /*1.4708e+02*/ DstarpipipiNorm[4], sigmaR1Variation_low, sigmaR1Variation_high);
	RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", /*9.3833e+01*/ DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
	RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", /*5.9429e+01*/ DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
	RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", /*1.2325e+01*/ DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
	RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", /*2.6268e+01*/ DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);

	//fractions of gauss functions
	RooRealVar f_1("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
	RooRealVar f_2("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);

	//add functions
	RooAddPdf Dstarpipipi_as_Dspipipi("Dstarpipipi_as_Dspipipi", "Dstarpipipi_as_Dspipipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2));

	//sum all background shapes
	//yields
	RooRealVar N_Dstarpipipi("N_Dstarpipipi","N_Dstarpipipi", 9500., 0., data->numEntries());

	RooRealVar N_comb("N_comb","N_comb", data->numEntries()/2., 0., data->numEntries());

        //sum background pdfs for BDT estimation 
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar f_bg("f_bg","f_bg", 0.5,0.,1.);
	RooAddPdf bkg("bkg", "bkg", RooArgList(Dstarpipipi_as_Dspipipi, bkg_exp), RooArgList(f_bg));
	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///total pdf
	///----------------------

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp), RooArgList(N_Bs, N_Dstarpipipi, N_comb));
	RooAbsPdf* pdfScan=new RooAddPdf("pdfScan", "pdfScan", RooArgList(DoubleGaussBs, bkg), RooArgList(N_Bs, n_bkg));

	///Fit
	RooFitResult *result;
if(!BDTscan)	result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
if(BDTscan)	result = pdfScan->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print();

	///calculate # (signal)background events in signal region

	Bs_MM.setRange("SigRange",meanBs1.getVal()-50.,meanBs1.getVal()+50.);
	
	RooAbsReal* S_fr= DoubleGaussBs.createIntegral(Bs_MM,NormSet(Bs_MM),Range("SigRange"));
	Double_t S = S_fr->getVal()*N_Bs.getVal();
	RooAbsReal* B_fr= bkg.createIntegral(Bs_MM,NormSet(Bs_MM),Range("SigRange"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
	cout<<"S/B= " << S/B<< endl;
	cout<<"S= " << S<< endl;
	cout<<"B= " << B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");



if(!BDTscan){	

	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(60));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m,Components(DoubleGaussBs),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(Dstarpipipi_as_Dspipipi),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(bkg_exp),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	//pdf->paramOn(frame_m,Layout(0.6));
	frame_m->Draw();
	c1->Print("eps/Final/3pi_BmassFit_12.eps");
}
/*
	RooPlot* frame_m3= Bplus_MM.frame("");
	data->plotOn(frame_m3);
	pdf->plotOn(frame_m3);
	gPad->SetLogy();
	frame_m3->Draw();
	c1->Print("Final/BmassFit_log.eps");
*/

	gPad->SetLogy(0);
	///fit results
	double chi2 = frame_m->chiSquare("pdf","data",9);
	double covmatr = result->covQual();
	double edm = result->edm();
	cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl;

  	// Construct a histogram with the pulls of the data w.r.t the curve
	RooHist* hresid = frame_m->residHist("data","pdf") ;
  	RooHist* hpull = frame_m->pullHist("data","pdf") ;

	// Create a new frame to draw the residual distribution and add the distribution to the frame
	RooPlot* frame2 = Bs_MM.frame(Title("Residual Distribution")) ;
if(!BDTscan){	
	frame2->SetTitle("");
	frame2->addPlotable(hresid,"P") ;
	frame2->Draw();
	c1->Print("eps/Final/3pi_residual_12.eps");
}
	// Create a new frame to draw the pull distribution and add the distribution to the frame
	RooPlot* frame3 = Bs_MM.frame(Title("Pull Distribution")) ;
if(!BDTscan){	
	frame3->SetTitle("");	
	frame3->SetLabelFont(62,"Y");
	frame3->addPlotable(hpull,"P") ;
	frame3->Draw();
	c1->Print("eps/Final/3pi_pull_12.eps");
}

return S/sqrt(S+B);

}



void makePlots(){
/*
	gStyle->SetOptStat(0);
    	gStyle->SetTitleXOffset(0.9);
    	//gStyle->SetTitleYOffset(1.3);
    	gStyle->SetTitleXSize(0.05);
    	gStyle->SetTitleYSize(0.05);
    	gStyle->SetTitleFont(42,"X");
    	gStyle->SetTitleFont(42,"Y");
    	//gStyle->SetLabelSize(0.033,"X");
    	//gStyle->SetLabelSize(0.033,"Y");
    	gStyle->SetLabelFont(42,"X");
    	gStyle->SetLabelFont(42,"Y");
   	gStyle->SetLabelOffset(0.01,"X");
    	gStyle->SetPadTickX(1);
    	gStyle->SetPadTickY(1);
    
    	TH1::SetDefaultSumw2();
    	TH2::SetDefaultSumw2();*/
    
    ///Load file
    TFile* file;
   file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/Bs2DsKpipi_Ds2pipipi_fullSelectionBDTG.root");
 //   file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT.root");    
    TTree* tree = (TTree*) file->Get("DecayTree");	

    double Bs_MM, Ds_MM;
    tree->SetBranchAddress("Bs_MM",&Bs_MM);
    tree->SetBranchAddress("Ds_MM",&Ds_MM);
    TH1D* massBs = new TH1D("mass of B_{s} candidates",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 5000., 5600.);
  //  TH1D* massBs = new TH1D("mass of B_{s} candidates",";m(D_{s}^{+} #pi^{-} #pi^{+} #pi^{-}) [MeV];Entries",40, 5000., 5600.);
    TH1D* massDs = new TH1D("mass of D_{s} candidates",";m(K^{-} K^{+} #pi^{+}) [MeV];Entries",40, 1900., 2050.);
/*
    double DTF_mKpipi, mKpipi, DTF_mKpi, mKpi, mPsi;
    double mB, sw;
    tree->SetBranchAddress("DTF_mKpipi",&DTF_mKpipi);
    tree->SetBranchAddress("mKpipi",&mKpipi);
    tree->SetBranchAddress("DTF_mKpi",&DTF_mKpi);
    tree->SetBranchAddress("mKpi",&mKpi);
    tree->SetBranchAddress("Bplus_MM",&mB);
    tree->SetBranchAddress("n_sig_sw",&sw);
    tree->SetBranchAddress("J_psi_1S_MM",&mPsi);*/
/*
    TH1D* hKpi = new TH1D("hKpi","",100,0.,4.5);
    TH1D* hKpi_DTF = new TH1D("hKpi_DTF",";m^{2}(K #pi) [GeV^{2}]; Yield [norm.]",100,0.,4.5);
    
    TH1D* hKpipi = new TH1D("hKpipi","",100,0.5,6);
    TH1D* hKpipi_DTF = new TH1D("hKpipi_DTF",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    TH1D* hKpipi_s = new TH1D("hKpipi_s",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; Candidates ",100,0.5,6);
    TH1D* hKpipi_DTF_s = new TH1D("hKpipi_DTF_s",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}]; Candidates",100,0.5,6);
    TH1D* hKpipi_b = new TH1D("hKpipi_b",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    TH1D* hKpipi_DTF_b = new TH1D("hKpipi_DTF_b",";m^{2}(K #pi #pi) [GeV^{2}]; Yield [norm.]",100,0.5,6);
    
    TH2D* mB_mKpipi = new TH2D("mB_mKpipi",";m^{2}(K^{+} #pi^{+} #pi^{-}) [GeV^{2}];m(B^{+}) [MeV]",50,0,6.5,50,5200,5600);
    TH1D* hpsi = new TH1D("mPsi","",100,3640,3740);
*/
    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		tree->GetEntry(i);

		massBs->Fill(Bs_MM);
		massDs->Fill(Ds_MM);
/*
		hJpsipipi->Fill(mJpsipipi/1000000,sw);
		hJpsipipi_DTF->Fill(DTF_mJpsipipi/1000000,sw);
		hKpipi->Fill(mKpipi/1000000,sw);
		hKpipi_DTF->Fill(DTF_mKpipi/1000000,sw);
		hKpi->Fill(mKpi/1000000,sw);
		hKpi_DTF->Fill(DTF_mKpi/1000000,sw);
		mB_mKpipi->Fill(mKpipi/1000000,mB);
		hpsi->Fill(mPsi,sw);
		if(abs(mB-5283.3)<60.){
			hKpipi_s->Fill(mKpipi/1000000);
			hKpipi_DTF_s->Fill(DTF_mKpipi/1000000);
		}
		if(mB>5400.){
			hKpipi_b->Fill(mKpipi/1000000);
			hKpipi_DTF_b->Fill(DTF_mKpipi/1000000);
		}*/
    }

    TCanvas* c= new TCanvas();
 //   massBs->Draw("e1"); c->Print("eps/mass_Bs_fullSel_12.eps");
  //  massDs->Draw("e1"); c->Print("eps/mass_Ds_fullSel_12.eps");

    massBs->Draw("e1"); c->Print("eps/mass_Bs_fullSel_12_Ds2pipipi.eps");
    massDs->Draw("e1"); c->Print("eps/mass_Ds_fullSel_12_Ds2pipipi.eps");
/*
    hpsi->Draw("e");
    c->Print("plots/mpsi.eps");

    hKpi_DTF->SetLineColor(kRed);
    hKpi_DTF->Draw("");   
    hKpi->Draw("esame");
    c->Print("plots/mKpi.eps");

    hKpipi_DTF->SetLineColor(kRed);
    hKpipi_DTF->Draw("");   
    hKpipi->Draw("esame");
    c->Print("plots/mKpipi.eps");

    hKpipi_b->SetLineColor(kRed);
    hKpipi_s->Draw("");   
    hKpipi_b->DrawNormalized("esame",51495);
    c->Print("plots/mKpipi_sb.eps");

    hKpipi_DTF_b->SetLineColor(kRed);
    hKpipi_DTF_s->Draw("");   
    hKpipi_DTF_b->DrawNormalized("esame",51495);
    c->Print("plots/mKpipi_DTF_sb.eps");

    mB_mKpipi->Draw();
    c->Print("plots/mB_mKpipi.eps");
    
    hJpsipipi_DTF->SetLineColor(kBlue);
    hJpsipipi_DTF->DrawNormalized("hist",1);   
    hJpsipipi->SetLineColor(kRed);
    hJpsipipi->DrawNormalized("histsame",1);
    c->Print("plots/mJpsipipi.eps");
    gPad->SetLogy(1);
    c->Print("plots/mJpsipipi_log.eps");*/
}

void quickFit(){

///Load file
	TFile* file;
	file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Bs2Dspipipi_with_BDT_variables_S21_PID.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_MM",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar Bs_MM("Bs_MM", "m(D_{s} #pi #pi #pi)", 5300., 5450.,"MeV/c^{2}");
	RooArgList list =  RooArgList(Bs_MM);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 39.27,5.,45.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 11.,5.,45.);
	RooRealVar a1("a1","a1", 1.49,0.,5.);
	RooRealVar n1("n1","n1", 1., 0., 100.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_MM, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_MM, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", Bs_MM, meanBs1, sigmaBs2, a1, n1);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0., 1.);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf signal("signal", "signal", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_GaussBs));

	//yields
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());
	RooRealVar n_bkg("n_bkg","n_bkg", data->numEntries()/2., 0., data->numEntries());

	//Exponential
	RooRealVar f_comb("f_comb" , "f_comb", 0.5, 0., 1.);
	RooRealVar exp_par("exp_par","#lambda",0.,-1.,1.);	
	RooExponential bkg("bkg_exp","exponential background",Bs_MM,exp_par);

	//add pdfs
	//RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, GaussBs2, bkg_exp), RooArgList(f_GaussBs, f_comb));
	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, bkg), RooArgList(N_Bs, n_bkg));


	///Fit
	RooFitResult *result;
	result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print();


	///calculate # (signal)background events in signal region

	Bs_MM.setRange("SigRange",meanBs1.getVal()-50.,meanBs1.getVal()+50.);
	
	RooAbsReal* S_fr= GaussBs1.createIntegral(Bs_MM,NormSet(Bs_MM),Range("SigRange"));
	Double_t S = S_fr->getVal()*N_Bs.getVal();
	RooAbsReal* B_fr= bkg.createIntegral(Bs_MM,NormSet(Bs_MM),Range("SigRange"));
	Double_t B = B_fr->getVal()*n_bkg.getVal();
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
	cout<<"S/B= " << S/B<< endl;
	cout<<"S= " << S<< endl;
	cout<<"B= " << B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= Bs_MM.frame();
	frame_m->SetTitle("");

	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
	pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
	pdf->plotOn(frame_m,Components(GaussBs1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	//pdf->plotOn(frame_m,Components(GaussBs2),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	//pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
	pdf->plotOn(frame_m,Components(bkg),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
	//pdf->paramOn(frame_m,Layout(0.75));
	data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
	//gPad->SetLogy();
	frame_m->Draw();
	c1->Print("eps/3pi_BmassShape_afterPreSel_12.eps");
}


void quickSignalEstimate(){

///Load file
	TFile* file;
	file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2Kpipi_with_BDT_variables_S21_PID.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("Bs_MM",1);

	///Fill all needed variables in RooDataSet

	///Bs
        RooRealVar Bs_MM("Bs_MM", "m(D_{s} K #pi #pi)", 5320., 5420.,"MeV/c^{2}");
	RooArgList list =  RooArgList(Bs_MM);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

	int expectedYield_12 = 842;
	int expectedYield_11 = 341;

	int expectedYield_11_Ds2pipipi = expectedYield_11 * 0.2 ;
	int expectedYield_11_Ds2Kpipi = expectedYield_11 * 0.12 ;

	//int expectedYield_11

	int all = data->numEntries();
	cout << "all Events in Signal Range:  " << all << endl;

	cout<<"expected Signal Yield:  " <<  expectedYield_11_Ds2Kpipi << endl;
	cout<<"Background Yield in Signal Region:  " <<  all - expectedYield_11_Ds2Kpipi << endl;

}


void iterateBDT(double startvalue, double stopvalue, double steps){

int numScans= (stopvalue - startvalue) / steps;
int Scan = 0;

double FOMarray[numScans];

for(double i = startvalue; i < (stopvalue + steps); i = i + steps){

	std::string cutAsString = boost::lexical_cast<std::string>(i);

	applyBDTcut(cutAsString);
	FOMarray[Scan] = fitBDTNorm();
	Scan++;
}

int readout = 0;

for(double i = startvalue; i < (stopvalue + steps); i = i + steps){

	cout << "<Figure of merit at BDTG-cut= " << i << "  is S/sqrt(S+B)= " <<  FOMarray[readout] << endl;
	readout++;
}
	
}

void MCStudies(){

///Load files to study peaking bg

//load mc file
TFile* file;
file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_with_BDT_variables_S21_PID_P_ETA.root");
TTree* tree = (TTree*) file->Get("DecayTree");	

//load weight files
TFile* fileUp1_w;
fileUp1_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/PIDEfficiencies_peakBG_Up1.root");
TTree* treeUp1_w = (TTree*) fileUp1_w->Get("CalibTool_PIDCalibTree");	

TFile* fileUp2_w;
fileUp2_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/PIDEfficiencies_peakBG_Up2.root");
TTree* treeUp2_w = (TTree*) fileUp2_w->Get("CalibTool_PIDCalibTree");

TFile* fileDown2_w;
fileDown2_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/PIDEfficiencies_peakBG_Down2.root");
TTree* treeDown2_w = (TTree*) fileDown2_w->Get("CalibTool_PIDCalibTree");

TFile* fileDown1_w;
fileDown1_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/PIDEfficiencies_peakBG_Down1.root");
TTree* treeDown1_w = (TTree*) fileDown1_w->Get("CalibTool_PIDCalibTree");


///define variables
//4-vectors
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;
TLorentzVector pi_plus2;
TLorentzVector pi_plus1;
TLorentzVector pi_minus;

//momentas
Double_t K_plus_fromDs_PX;
Double_t K_plus_fromDs_PY;
Double_t K_plus_fromDs_PZ;
Double_t pi_plus1_PX;
Double_t pi_plus1_PY;
Double_t pi_plus1_PZ;
Double_t K_minus_fromDs_PX;
Double_t K_minus_fromDs_PY;
Double_t K_minus_fromDs_PZ;
Double_t pi_minus_fromDs_PX;
Double_t pi_minus_fromDs_PY;
Double_t pi_minus_fromDs_PZ;
Double_t pi_minus_PX;
Double_t pi_minus_PY;
Double_t pi_minus_PZ;
Double_t pi_plus2_PX;
Double_t pi_plus2_PY;
Double_t pi_plus2_PZ;

//weights
Float_t Event_PIDCalibEffWeight_Up1;
Float_t Event_PIDCalibEffWeight_Up2;
Float_t Event_PIDCalibEffWeight_Down1;
Float_t Event_PIDCalibEffWeight_Down2;

//masses
double massKaon = 493.68;
double massPion = 139.57;


//link variables to tree
treeUp1_w -> SetBranchAddress( "Event_PIDCalibEffWeight" , &Event_PIDCalibEffWeight_Up1 );
treeUp2_w -> SetBranchAddress( "Event_PIDCalibEffWeight" , &Event_PIDCalibEffWeight_Up2 );
treeDown1_w -> SetBranchAddress( "Event_PIDCalibEffWeight" , &Event_PIDCalibEffWeight_Down1 );
treeDown2_w -> SetBranchAddress( "Event_PIDCalibEffWeight" , &Event_PIDCalibEffWeight_Down2 );


tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );

tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );

tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

tree -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
tree -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
tree -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ )
;
tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );

tree -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
tree -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
tree -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );


TH1D* massBs_Up1 = new TH1D("mass of B_{s} candidates up1",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Up1_fitRange = new TH1D("mass of B_{s} candidates up1 in fitRange",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Up2 = new TH1D("mass of B_{s} candidates up2",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Up2_fitRange = new TH1D("mass of B_{s} candidates up2 in fitRange",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Down1 = new TH1D("mass of B_{s} candidates Down1",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Down1_fitRange = new TH1D("mass of B_{s} candidates Down1 in fitRange",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Down2 = new TH1D("mass of B_{s} candidates Down2",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);
TH1D* massBs_Down2_fitRange = new TH1D("mass of B_{s} candidates Down2 in fitRange",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",40, 4500., 6500.);


double Bs_MM = 0;

///loop over events

///Up1
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	treeUp1_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massKaon); //flip mass hypothesis here
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Up1->Fill(Bs_MM,Event_PIDCalibEffWeight_Up1);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Up1_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Up1);
	}

///Up2
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	treeUp2_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion); 
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon); //flip mass hypothesis here

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Up2->Fill(Bs_MM,Event_PIDCalibEffWeight_Up2);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Up2_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Up2);
	}

///Down1
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	treeDown1_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massKaon); //flip mass hypothesis here
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Down1->Fill(Bs_MM,Event_PIDCalibEffWeight_Down1);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Down1_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Down1);
	}

///Down2
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	treeDown2_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion); 
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon); //flip mass hypothesis here

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Down2->Fill(Bs_MM,Event_PIDCalibEffWeight_Down2);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Down2_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Down2);
	}


TCanvas* c= new TCanvas();

massBs_Up1->Sumw2(); massBs_Up1->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up1.eps");
massBs_Up1_fitRange->Sumw2(); massBs_Up1_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up1_inRange.eps");

massBs_Up2->Sumw2(); massBs_Up2->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up2.eps");
massBs_Up2_fitRange->Sumw2(); massBs_Up2_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up2_inRange.eps");

massBs_Down1->Sumw2(); massBs_Down1->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down1.eps");
massBs_Down1_fitRange->Sumw2(); massBs_Down1_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down1_inRange.eps");

massBs_Down2->Sumw2(); massBs_Down2->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down2.eps");
massBs_Down2_fitRange->Sumw2(); massBs_Down2_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down2_inRange.eps");

}



int main(){
    time_t startTime = time(0);

 //  quickSignalEstimate();
  //  iterateBDT(-0.9,0.9,0.05);
  //  preselect();
  //  applyBDTcut("0.1949");
    //makePlots();
  //  //chooseBestPV("BFit","/auto/data/dargent/Bu2psiKpipi/data/data_preselected.root");
    //chooseBestPV("psiFit","/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_BFit.root");
    //chooseBestPV("psiDTF","/auto/data/dargent/Bu2psiKpipi/data/data_preselected_bestPV_psiFit.root");
    //fitPreselected();  
    //fitPreselected_psiConstrained();      
    //addVarsForBDT();
    //chooseMultCand();
/*
RooRealVar Bs_MM("Bs_MM", "m(D_{s}*K#pi#pi)", 5000., 5350.,"MeV/c^{2}");

//Load file
TFile* file;
file= new TFile("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/DsstKpipi.root");
TTree* tree = (TTree*) file->Get("DecayTree");
tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_MM",1);

//Fill needed variable in RooDataSet
RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*tree));

//try RooKernel
RooKeysPdf kest("kest","kest",Bs_MM,*data,RooKeysPdf::NoMirror);

TCanvas* c1= new TCanvas("");

// Plot kernel estimation pdfs with and without mirroring over data
RooPlot* frame = Bs_MM.frame(Title("Adaptive kernel estimation pdf w/o mirroring"),Bins(20)) ;
data->plotOn(frame);
kest.plotOn(frame,LineStyle(kDashed),LineColor(kRed)) ;
frame->Draw();
c1->Print("eps/BkgShape/RooKeyKernelstimator.eps");
*/
  // addCut();
    fitBDT();
 // fitBDTNorm();
 //  fitBGShapeNorm();
  // fitBGShapeNormKpipi();
  // fitBGShapeNormDstKpipi();
   //fitBGShapethreePi();
   //fitBGShapethreePiDstar();
 //  quickFit();
 //   MCStudies();

    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;

return 0;
}

