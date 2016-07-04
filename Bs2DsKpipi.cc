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

int preselect(bool MC=true , bool preselection =false) {

    //specifiy decay channel of Ds
    bool Ds2KKpi = false;
    bool Ds2Kpipi = false;
    bool Ds2pipipi = false;

    //specification for normalisation channel
    bool BsDspipipi = true;

    // extended offline selection ?
    bool offline_selection = false;

    ///Load file
    TChain* tree = 0;
	
    //Ds2KKpi case
    if(Ds2KKpi){
	tree=new TChain("Bs2DsKpipi_Ds2KKpi_Tuple/DecayTree");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12D/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12U/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12D-MC/*.root");
	//tree->Add("/auto/data/kecke/B2DKPiPi/12U-MC/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/11U/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/11D/*.root");
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
      tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-U/*.root");
      tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-D/*.root");
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
	if(Ds2KKpi) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_TriggeredandTruthMatched.root","RECREATE");

	//Ds2Kpipi case
	if(Ds2Kpipi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2Kpipi_preselected.root","RECREATE");

	// Ds2pipipi case
	if(Ds2pipipi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_preselected.root","RECREATE");

	//Bs->Dspipipi normalisation case
	if(BsDspipipi)output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_TriggeredandTruthMatched.root","RECREATE");
		
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
	string trigger_cuts_L0 = "(Bs_L0HadronDecision_TOS==1 || Bs_L0Global_TIS==1)";
	string trigger_cuts_HLT1 = "(Bs_Hlt1TrackAllL0Decision_TOS==1)";
	string trigger_cuts_HLT2 = "(Bs_Hlt2Topo2BodyBBDTDecision_TOS==1 || Bs_Hlt2Topo3BodyBBDTDecision_TOS==1 || Bs_Hlt2Topo4BodyBBDTDecision_TOS==1)";

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
	 

    string MCTRUE_cut = "abs(Bs_TRUEID)==531 && abs(Ds_TRUEID)==431 && abs(K_plus_TRUEID)==321 && abs(pi_plus_TRUEID)==211 && abs(pi_minus_TRUEID)==211 && abs(K_plus_fromDs_TRUEID)==321 && abs(K_minus_fromDs_TRUEID)==321 && abs(pi_minus_fromDs_TRUEID)==211" ;

    string MCTRUE_cut_Norm = "abs(Bs_TRUEID)==531 && abs(Ds_TRUEID)==431 && abs(pi_plus1_TRUEID)==211 && abs(pi_plus2_TRUEID)==211 && abs(pi_minus_TRUEID)==211 && abs(K_plus_fromDs_TRUEID)==321 && abs(K_minus_fromDs_TRUEID)==321 && abs(pi_minus_fromDs_TRUEID)==211" ;
    
/*    
	cuts.append(mass_cuts);
	cuts.append("&&");
	cuts.append(pid_cuts);

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
    }*/
    if (MC) {
        cuts.append(MCTRUE_cut_Norm);
        cuts.append("&&");
        cuts.append(trigger_cuts_L0);
        cuts.append("&&");
        cuts.append(trigger_cuts_HLT1);
        cuts.append("&&");
        cuts.append(trigger_cuts_HLT2);
    }
    //string cuts_noPID;
    //cuts_noPID.append(MCTRUE_cut);
    
    TTree* tree_sel = tree->CopyTree(cuts.c_str());
    //TTree* tree_noPID = tree->CopyTree(cuts_noPID.c_str());
    cout << "New file contains " << tree_sel->GetEntries() << " events" <<  endl;
    //if(MC)cout << "Signal Eff = " << (double) ((double) tree_sel->GetEntries() / (double) tree_noPID->GetEntries()) << endl;
    
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

TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
TTree* tree = (TTree*) file->Get("DecayTree");


   //TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/mc2011_Bs2Dspipipi_Ds2KKpi_forBDT_P_ETA_pi_plus2.root","RECREATE");
/*
   TTree* new_tree_BDTG = tree->CopyTree("(Bs_L0HadronDecision_TOS==1 || Bs_L0Global_TIS==1) && (Bs_Hlt1TrackAllL0Decision_TOS==1) && K_plus_fromDs_P>1600 && K_minus_fromDs_P>1600 && pi_minus_fromDs_P>1600");

   TTree* new_tree_BDTG = tree->CopyTree("pi_plus2_ETA > 1.5 && pi_plus2_ETA < 5.0  && pi_plus2_P < 100000 && pi_plus2_P > 3000");
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ_harderLbVeto.root","RECREATE");
   //TTree* new_tree_BDTG = tree->CopyTree("Ds_MM > 1950 && Ds_MM < 1990");
   //TTree* new_tree_BDTG = tree->CopyTree("K_plus_PIDK > 8");
   //TTree* new_tree_BDTG = tree->CopyTree("(Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) > 0");
   TTree* new_tree_BDTG = tree->CopyTree("(K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > -1");
   new_tree_BDTG->Write();
   output_BDTG->Close();


   //close file at the end
   file->Close();

}


void applyBDTcut(string cutoff){
   ///Load file

  //TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data12_Ds2KKpi_BDT_tmp.root");
 // TTree* tree = (TTree*) file->Get("DecayTree");	
TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data12_Ds2KKpi_BDT_tmp.root");
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


    //BDTG application for 2011 Signal Channel, 0.8300 for maximum S/sqrt(S+B) with S = 617 and B = 37476
/*
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

    //BDTG application for 2012 Signal Channel (9500BG/450S), 0.8040 for maximum S/sqrt(S+B) , with S = 726 , B = 18751 
/*
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

    //BDTG application for 2012 normalization Channel (9500BG/450S) , 0.1458 for maximum S/sqrt(S+B) with S = 13867 , B = 23588
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();

/*
    //BDTG application for 2011 normalization Channel (32.500BG/6500S) ------> BDTG > 0.5571 with S=14511 and B=98053
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root","RECREATE");
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
RooPlot* frame_m= Bs_Mass.frame();
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
RooPlot* frame_m= Bs_Mass.frame();
frame_m->SetTitle("");

data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
pdf->plotOn(frame_m,Components(CB2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
frame_m->Draw();
c1->Print("eps/BkgShape/Bs2Dsstarpipipi_as_DsKpipi.eps");

return fitValues;
}

double *fitBDTNorm(double fitvalues[15], bool sevenTeV, bool fitSimultan, bool doToys=false, bool altModel = false){

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

	///Load file

	///seperate fits for 7 & 8 TeV data (fitSimultan = false)
	TFile* file;
	TTree* tree;
	if(!fitSimultan){
		if(sevenTeV) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
		if(!sevenTeV) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree = (TTree*) file->Get("DecayTree");	
   		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("DTF_Bs_M",1);
	}

	///one simultaneous fit for 7 & 8 TeV data (fitSimultan = true) 
	TFile* file11;
	TFile* file12;
	TTree* tree11;
	TTree* tree12;
	if(fitSimultan){
		file11 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}
	
	

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} #pi #pi #pi)", 4800., 5800.,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);

        RooDataSet* data = 0;
	RooDataSet* data11 = 0;
	RooDataSet* data12 = 0;
	if(!fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree));
	}
	if(fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree11));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree11));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree12));
	}

	///define category to distinguish between 2011 and 2012 data
	RooCategory sample_year("sample_year","sample_year") ;
	sample_year.defineType("y11");
	sample_year.defineType("y12");

	RooDataSet* combData;
	if(fitSimultan) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_year),Import("y11",*data11),Import("y12",*data12));


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
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.066);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
*/

	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5366.7,5345.,5380.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 32.56);//, 5.,60.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 13.54,10.,20.);
	RooRealVar a1("a1","a1", -1.83);
	RooRealVar n1("n1","n1", 1.48);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanBs1, sigmaBs2, a1, n1);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.18);//, 0.,0.5);
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_GaussBs));
	RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));

	//signal pdf
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential for first comb. BG. component
	RooRealVar exp_par("exp_par","#lambda_{1}",-1.6508e-03,-10.,0.);	
	RooRealVar exp_par_11("exp_par_11","#lambda_{1} 11",-1.6508e-03,-10.,0.);
	RooRealVar exp_par_12("exp_par_12","#lambda_{1} 12",-1.6508e-03,-10.,0.);
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

	RooRealVar* sigmaL1;
	RooRealVar* sigmaR1;
	RooRealVar* sigmaL2;
	RooRealVar* sigmaR2;
	RooRealVar* sigmaL3;
	RooRealVar* sigmaR3;

	RooRealVar* f_1;
	RooRealVar* f_2;

	///mean of gaussians
	RooRealVar mean1("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
	RooRealVar mean2("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
	RooRealVar mean3("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);

	///width and fractions of gaussians 2011

	if(sevenTeV){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);// sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);
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
	sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
	sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
	sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
	sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);//, sigmaL3Variation_low, sigmaL3Variation_high);
	sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);
	}

	if(fitSimultan){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);//, sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	//bifurcated gaussians
	RooBifurGauss BifGauss1("BifGauss1","BifGauss1", DTF_Bs_M, mean1, *sigmaL1,*sigmaR1);
	RooBifurGauss BifGauss2("BifGauss2","BifGauss2", DTF_Bs_M, mean2, *sigmaL2,*sigmaR2);
	RooBifurGauss BifGauss3("BifGauss3","BifGauss3", DTF_Bs_M, mean3, *sigmaL3,*sigmaR3);


	//add functions
	RooAddPdf Dstarpipipi_as_Dspipipi("Dstarpipipi_as_Dspipipi", "Dstarpipipi_as_Dspipipi", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(*f_1,*f_2));

	//sum all background shapes
	//yields
	RooRealVar N_Dstarpipipi("N_Dstarpipipi","N_Dstarpipipi", 9500., 0., data->numEntries());

	RooRealVar N_comb("N_comb","N_comb", data->numEntries()/2., 0., data->numEntries());

        //sum background pdfs for BDT estimation 
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()/2., 0., data->numEntries());
	RooRealVar f_bg("f_bg","f_bg", 0.5,0.,1.);
	RooAddPdf bkg("bkg", "bkg", RooArgList(Dstarpipipi_as_Dspipipi, bkg_exp), RooArgList(f_bg));

	///add yields for possible simultaneous fit
	RooRealVar N_Bs_11("N_Bs_11", "2011 #B_{s}", data11->numEntries()/2., 0., data11->numEntries());
	RooRealVar N_Bs_12("N_Bs_12", "2012 #B_{s}", data12->numEntries()/2., 0., data12->numEntries());
	RooRealVar N_Dstarpipipi_11("N_Dstarpipipi_11","N_Dstarpipipi_11", 9500., 0., data11->numEntries());
	RooRealVar N_Dstarpipipi_12("N_Dstarpipipi_12","N_Dstarpipipi_12", 9500., 0., data12->numEntries());
	RooRealVar N_comb_11("N_comb_11","N_comb_11", data11->numEntries()/2., 0., data11->numEntries());
	RooRealVar N_comb_12("N_comb_12","N_comb_12", data12->numEntries()/2., 0., data12->numEntries());


	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///total pdf
	///----------------------

	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp), RooArgList(N_Bs, N_Dstarpipipi, N_comb));
	RooAbsPdf* pdfScan=new RooAddPdf("pdfScan", "pdfScan", RooArgList(DoubleGaussBs, bkg), RooArgList(N_Bs, n_bkg));

	///construct a simultaneous pdf using category sample as index
	RooSimultaneous* simPdf;
	RooAbsPdf* pdf11;
	RooAbsPdf* pdf12;
	if(fitSimultan){
		simPdf = new RooSimultaneous("simPdf","simultaneous pdf",sample_year);
		pdf11=new RooAddPdf("pdf11", "pdf11", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_11), RooArgList(N_Bs_11, N_Dstarpipipi_11, N_comb_11));
		pdf12=new RooAddPdf("pdf12", "pdf12", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_12), RooArgList(N_Bs_12, N_Dstarpipipi_12, N_comb_12));
		simPdf->addPdf(*pdf11,"y11");
		simPdf->addPdf(*pdf12,"y12");
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

	string eleven = "_11.eps";
	string twelve = "_12.eps";

	// 2011 results
	if(sevenTeV){
		BsMassDistribution.append(eleven);
		BsMassResidual.append(eleven);
		BsMassPull.append(eleven);
		BsMassSweight.append(eleven);
	}

	//2012 results
	if(!sevenTeV){
		BsMassDistribution.append(twelve);
		BsMassResidual.append(twelve);
		BsMassPull.append(twelve);
		BsMassSweight.append(twelve);
	}


	if((!BDTscan) && (!fitSimultan)){	

		data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(60));
		pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlack),LineWidth(2));
		pdf->plotOn(frame_m,Components(DoubleGaussBs),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(Dstarpipipi_as_Dspipipi),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		pdf->plotOn(frame_m,Components(bkg_exp),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m->Draw();
		c1->Print(BsMassDistribution.c_str());
	}

	if((!BDTscan) && fitSimultan){
		combData->plotOn(frame_m_11,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_11,Name("pdf11"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11->Draw();
		c1->Print("eps/Final/3pi_BmassFit_sim11.eps");

		combData->plotOn(frame_m_12,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_12,Name("pdf12"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12->Draw();
		c1->Print("eps/Final/3pi_BmassFit_sim12.eps");
	}


	///fit results
	double chi2 = 0;
	if(!fitSimultan) chi2 = frame_m->chiSquare("pdf","data",15);
	if(fitSimultan) chi2 = (frame_m_11->chiSquare("pdf11","data11",15) + frame_m_12->chiSquare("pdf12","data12",15)) / 2;
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

	if((!BDTscan) && (!fitSimultan)){	
		frame2->SetTitle("");
		frame2->addPlotable(hresid,"P") ;
		frame2->Draw();
		c1->Print(BsMassResidual.c_str());
	}
	if((!BDTscan) && fitSimultan){
		frame2_11->SetTitle("");
		frame2_11->addPlotable(hresid11,"P") ;
		frame2_11->Draw();
		c1->Print("eps/Final/3pi_residual_sim11.eps");

		frame2_12->SetTitle("");
		frame2_12->addPlotable(hresid12,"P") ;
		frame2_12->Draw();
		c1->Print("eps/Final/3pi_residual_sim12.eps");
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
		c1->Print(BsMassPull.c_str());
	}
	if((!BDTscan) && fitSimultan){
		frame3_11->SetTitle("");	
		frame3_11->SetLabelFont(62,"Y");
		frame3_11->addPlotable(hpull11,"P") ;
		frame3_11->Draw();
		c1->Print("eps/Final/3pi_pull_sim11.eps");

		frame3_12->SetTitle("");	
		frame3_12->SetLabelFont(62,"Y");
		frame3_12->addPlotable(hpull12,"P") ;
		frame3_12->Draw();
		c1->Print("eps/Final/3pi_pull_sim12.eps");
	}


	///save fitavalues in array to be use in fitBDT()
	fitvalues[0] = mean1.getVal();
	fitvalues[1] = mean2.getVal();
	fitvalues[2] = mean3.getVal();
	fitvalues[3] = sigmaL1->getVal();
	fitvalues[4] = sigmaR1->getVal();
	fitvalues[5] = sigmaL2->getVal();
	fitvalues[6] = sigmaR2->getVal();
	fitvalues[7] = sigmaL3->getVal();
	fitvalues[8] = sigmaR3->getVal();
	fitvalues[9] = f_1->getVal();
	fitvalues[10] = f_2->getVal();

	if(!fitSimultan){ 
		fitvalues[11] = N_Dstarpipipi.getVal();
		fitvalues[12] = N_Bs.getVal();
		fitvalues[13] = 0;
		fitvalues[14] = 0;
	}
	if(fitSimultan){
		fitvalues[11] = N_Dstarpipipi_11.getVal();
		fitvalues[12] = N_Bs_11.getVal();
		fitvalues[13] = N_Dstarpipipi_12.getVal();
		fitvalues[14] = N_Bs_12.getVal();
	}

if(doToys && (!fitSimultan)){

	Int_t SampleEvents;
	//if(sevenTeV) SampleEvents = 41663;
	//if(!sevenTeV) SampleEvents = 93604;
	SampleEvents = 250000;


	TFile* toy_output;
	if(sevenTeV) toy_output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Ds3pi_11_final_toy.root","RECREATE");
	if(!sevenTeV) toy_output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Ds3pi_12_final_toy.root","RECREATE");

	RooMCStudy* Toystudy = new RooMCStudy(*pdf,RooArgSet(DTF_Bs_M),Binned(kFALSE),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
	Toystudy->generateAndFit(1000,SampleEvents,"kTRUE");
	RooDataSet toy_data = Toystudy->fitParDataSet();

	// Make plots of the distributions of mean, the error on mean and the pull of mean
	
	RooPlot* frame_Bs_value = Toystudy->plotParam(N_Bs,Bins(40));
	RooPlot* frame_Bs_error = Toystudy->plotError(N_Bs,Bins(40));
	RooPlot* frame_Bs_pull = Toystudy->plotPull(N_Bs,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_Bs_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_Distribution_12.eps");
 	frame_Bs_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_Error_12.eps");
 	frame_Bs_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Bs_PullDistribution_12.eps");


	
	RooPlot* frame_N_Dstarpipipi_value = Toystudy->plotParam(N_Dstarpipipi,Bins(40));
	RooPlot* frame_N_Dstarpipipi_error = Toystudy->plotError(N_Dstarpipipi,Bins(40));
	RooPlot* frame_N_Dstarpipipi_pull = Toystudy->plotPull(N_Dstarpipipi,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_N_Dstarpipipi_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_Distribution_12.eps");
	frame_N_Dstarpipipi_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_Error_12.eps");
	frame_N_Dstarpipipi_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_Dstarpipipi_PullDistribution_12.eps");



	RooPlot* frame_N_comb_value = Toystudy->plotParam(N_comb,Bins(40));
	RooPlot* frame_N_comb_error = Toystudy->plotError(N_comb,Bins(40));
	RooPlot* frame_N_comb_pull = Toystudy->plotPull(N_comb,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_N_comb_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_Distribution_12.eps");
	frame_N_comb_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_Error_12.eps");
	frame_N_comb_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/N_comb_PullDistribution_12.eps");



	RooPlot* frame_exp_par_value = Toystudy->plotParam(exp_par,Bins(40));
	RooPlot* frame_exp_par_error = Toystudy->plotError(exp_par,Bins(40));
	RooPlot* frame_exp_par_pull = Toystudy->plotPull(exp_par,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_exp_par_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_Distribution_12.eps");
	frame_exp_par_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_Error_12.eps");
	frame_exp_par_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/exp_par_PullDistribution_12.eps");



	RooPlot* frame_f_1_value = Toystudy->plotParam(*f_1,Bins(40));
	RooPlot* frame_f_1_error = Toystudy->plotError(*f_1,Bins(40));
	RooPlot* frame_f_1_pull = Toystudy->plotPull(*f_1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_f_1_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/f_1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/f_1_Distribution_12.eps");
	frame_f_1_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/f_1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/f_1_Error_12.eps");
	frame_f_1_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/f_1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/f_1_PullDistribution_12.eps");



	RooPlot* frame_mean1_value = Toystudy->plotParam(mean1,Bins(40));
	RooPlot* frame_mean1_error = Toystudy->plotError(mean1,Bins(40));
	RooPlot* frame_mean1_pull = Toystudy->plotPull(mean1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_mean1_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Distribution_12.eps");
	frame_mean1_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Error_12.eps");
	frame_mean1_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_PullDistribution_12.eps");



	RooPlot* frame_mean2_value = Toystudy->plotParam(mean2,Bins(40));
	RooPlot* frame_mean2_error = Toystudy->plotError(mean2,Bins(40));
	RooPlot* frame_mean2_pull = Toystudy->plotPull(mean2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_mean2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Distribution_12.eps");
	frame_mean2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Error_12.eps");
	frame_mean2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_PullDistribution_12.eps");



	RooPlot* frame_mean3_value = Toystudy->plotParam(mean3,Bins(40));
	RooPlot* frame_mean3_error = Toystudy->plotError(mean3,Bins(40));
	RooPlot* frame_mean3_pull = Toystudy->plotPull(mean3,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_mean3_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean3_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean3_Distribution_12.eps");
	frame_mean3_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean3_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean3_Error_12.eps");
	frame_mean3_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean3_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean3_PullDistribution_12.eps");


	RooPlot* frame_meanBs1_value = Toystudy->plotParam(meanBs1,Bins(40));
	RooPlot* frame_meanBs1_error = Toystudy->plotError(meanBs1,Bins(40));
	RooPlot* frame_meanBs1_pull = Toystudy->plotPull(meanBs1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_meanBs1_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_Distribution_12.eps");
	frame_meanBs1_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_Error_12.eps");
	frame_meanBs1_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/meanBs1_PullDistribution_12.eps");



	RooPlot* frame_sigmaBs2_value = Toystudy->plotParam(sigmaBs2,Bins(40));
	RooPlot* frame_sigmaBs2_error = Toystudy->plotError(sigmaBs2,Bins(40));
	RooPlot* frame_sigmaBs2_pull = Toystudy->plotPull(sigmaBs2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaBs2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_Distribution_12.eps");
	frame_sigmaBs2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_Error_12.eps");
	frame_sigmaBs2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaBs2_PullDistribution_12.eps");



	RooPlot* frame_sigmaL1_value = Toystudy->plotParam(*sigmaL1,Bins(40));
	RooPlot* frame_sigmaL1_error = Toystudy->plotError(*sigmaL1,Bins(40));
	RooPlot* frame_sigmaL1_pull = Toystudy->plotPull(*sigmaL1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaL1_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_Distribution_12.eps");
	frame_sigmaL1_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_Error_12.eps");
	frame_sigmaL1_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL1_PullDistribution_12.eps");



	RooPlot* frame_sigmaL2_value = Toystudy->plotParam(*sigmaL2,Bins(40));
	RooPlot* frame_sigmaL2_error = Toystudy->plotError(*sigmaL2,Bins(40));
	RooPlot* frame_sigmaL2_pull = Toystudy->plotPull(*sigmaL2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaL2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_Distribution_12.eps");
	frame_sigmaL2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_Error_12.eps");
	frame_sigmaL2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL2_PullDistribution_12.eps");


	RooPlot* frame_sigmaR2_value = Toystudy->plotParam(*sigmaR2,Bins(40));
	RooPlot* frame_sigmaR2_error = Toystudy->plotError(*sigmaR2,Bins(40));
	RooPlot* frame_sigmaR2_pull = Toystudy->plotPull(*sigmaR2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaR2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_Distribution_12.eps");
	frame_sigmaR2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_Error_12.eps");
	frame_sigmaR2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR2_PullDistribution_12.eps");



	RooPlot* frame_sigmaR3_value = Toystudy->plotParam(*sigmaR3,Bins(40));
	RooPlot* frame_sigmaR3_error = Toystudy->plotError(*sigmaR3,Bins(40));
	RooPlot* frame_sigmaR3_pull = Toystudy->plotPull(*sigmaR3,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaR3_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_Distribution_12.eps");
	frame_sigmaR3_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_Error_12.eps");
	frame_sigmaR3_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaR3_PullDistribution_12.eps");



	RooPlot* frame_sigmaL3_value = Toystudy->plotParam(*sigmaL3,Bins(40));
	RooPlot* frame_sigmaL3_error = Toystudy->plotError(*sigmaL3,Bins(40));
	RooPlot* frame_sigmaL3_pull = Toystudy->plotPull(*sigmaL3,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaL3_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_Distribution_12.eps");
	frame_sigmaL3_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_Error_12.eps");
	frame_sigmaL3_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/sigmaL3_PullDistribution_12.eps");

	toy_data.Write();
	toy_output->Close();
}


if(sWeight && (!fitSimultan)){

		sigmaL1->setConstant();
 		sigmaR1->setConstant();
 		sigmaL2->setConstant();
 		sigmaR2->setConstant();
 		sigmaL3->setConstant();
 		sigmaR3->setConstant();
 		f_1->setConstant();
 		f_2->setConstant();
		mean1.setConstant();
		mean2.setConstant();
		mean3.setConstant();
		meanBs1.setConstant();
		sigmaBs1.setConstant();
		sigmaBs2.setConstant();
		f_GaussBs.setConstant();
		exp_par.setConstant();

	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(N_Bs, N_Dstarpipipi, N_comb)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("DTF_Bs_M,N_Bs_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print(BsMassSweight.c_str());

    		///Create output file
   		 TFile* output;
		 if(sevenTeV) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_final_sweight.root","RECREATE");
		 if(!sevenTeV) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("DTF_Bs_M > 4800 && DTF_Bs_M < 5800");
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

	if(altModel && fitSimultan){

		// set polynomial background
		RooRealVar a1_11("a1_11","a1_11",5.,-100.,100.);
		RooRealVar a2_11("a2_11","a2_11",0.,-10.,10.);
		RooRealVar a3_11("a3_11","a3_11",0.,-1.,1.);
		RooPolynomial polyBG_11("polyBG_11","polyBG_11",DTF_Bs_M,RooArgList(a1_11,a2_11),0);

		RooRealVar a1_12("a1_12","a1_12",5.,-100.,100.);
		RooRealVar a2_12("a2_12","a2_12",0.,-10.,10.);
		RooRealVar a3_12("a3_12","a3_12",0.,-1.,1.);
		RooPolynomial polyBG_12("polyBG_12","polyBG_12",DTF_Bs_M,RooArgList(a1_12,a2_12),0);

		mean1.setVal(DstarpipipiNorm[0]);
		mean2.setVal(DstarpipipiNorm[1]);
		mean3.setVal(DstarpipipiNorm[2]);
/*
		sigmaL1->setVal(DstarpipipiNorm[3]);
		sigmaR1->setVal(DstarpipipiNorm[4]);
		sigmaL2->setVal(DstarpipipiNorm[5]);
		sigmaR2->setVal(DstarpipipiNorm[6]);
		sigmaL3->setVal(DstarpipipiNorm[7]);
		sigmaR3->setVal(DstarpipipiNorm[8]);
*/
		f_1->setVal(DstarpipipiNorm[9]);
		f_2->setVal(DstarpipipiNorm[10]);
/*
		sigmaL1->setConstant();
 		sigmaR1->setConstant();
 		sigmaL2->setConstant();
 		sigmaR2->setConstant();
 		sigmaL3->setConstant();
 		sigmaR3->setConstant();
*/
		RooKeysPdf altBGModel_11("altBGModel_11","altBGModel_11",DTF_Bs_M,*combData);
		RooKeysPdf altBGModel_12("altBGModel_12","altBGModel_12",DTF_Bs_M,*combData);

 		f_1->setConstant();
 		f_2->setConstant();
		mean1.setConstant();
		mean2.setConstant();
		mean3.setConstant();

		RooSimultaneous* simPdf_alt;
		RooAbsPdf* pdf11_alt;
		RooAbsPdf* pdf12_alt;
		simPdf_alt = new RooSimultaneous("simPdf_alt","alternative simultaneous pdf",sample_year);
		pdf11_alt=new RooAddPdf("pdf11_alt", "pdf11_alt", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_11), RooArgList(N_Bs_11, N_Dstarpipipi_11, N_comb_11));
		pdf12_alt=new RooAddPdf("pdf12_alt", "pdf12_alt", RooArgList(DoubleGaussBs, Dstarpipipi_as_Dspipipi ,bkg_exp_12), RooArgList(N_Bs_12, N_Dstarpipipi_12, N_comb_12));
		simPdf_alt->addPdf(*pdf11_alt,"y11");
		simPdf_alt->addPdf(*pdf12_alt,"y12");

		//fit alt. models
		RooFitResult *result_alt;
		result_alt = simPdf_alt->fitTo(*combData,Save(kTRUE),Extended(kTRUE),NumCPU(3));

		cout << "result is --------------- "<<endl;
		result_alt->Print();

		//draw alt. models
        	RooPlot* frame_m_11_alt= DTF_Bs_M.frame();
        	RooPlot* frame_m_12_alt= DTF_Bs_M.frame();

      	  	frame_m_11_alt->SetTitle("");
       		frame_m_12_alt->SetTitle("");

		combData->plotOn(frame_m_11_alt,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(60));
		simPdf_alt->plotOn(frame_m_11_alt,Name("pdf11_alt"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11_alt->Draw();
		c1->Print("eps/altModels/3pi_BmassFit_sim11_fixedBG.eps");

		combData->plotOn(frame_m_12_alt,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(60));
		simPdf_alt->plotOn(frame_m_12_alt,Name("pdf12_alt"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12_alt->Draw();
		c1->Print("eps/altModels/3pi_BmassFit_sim12_fixedBG.eps");

	}

return fitvalues;

}


void fitBDT(){

	bool sWeight=false;
	bool sevenTeV=false; // 2011 = true, 2012 = false
	bool doToys=true;
	bool fitSimultan=false;
	bool altModel=false;

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
	
	///seperate fits for 7 & 8 TeV data (fitSimultan = false)
	TFile* file;
	TTree* tree;
	if(!fitSimultan){
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
	if(fitSimultan){
		file11 = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree11 = (TTree*) file11->Get("DecayTree");
		tree11->SetBranchStatus("*",0);
		tree11->SetBranchStatus("DTF_Bs_M",1);
		file12 = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG_tightDCUT_DZ.root");
		tree12 = (TTree*) file12->Get("DecayTree");
		tree12->SetBranchStatus("*",0);
		tree12->SetBranchStatus("DTF_Bs_M",1);
	}

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} K #pi #pi)", 4800., 5800.,"MeV");
	RooArgList list =  RooArgList(DTF_Bs_M);

        RooDataSet* data = 0;
	RooDataSet* data11 = 0;
	RooDataSet* data12 = 0;
	if(!fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree));
	}
	if(fitSimultan){
		data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree11));
		data11 = new RooDataSet("data11","data11",RooArgSet(DTF_Bs_M),Import(*tree11));
		data12 = new RooDataSet("data12","data12",RooArgSet(DTF_Bs_M),Import(*tree12));
	}

	///define category to distinguish between 2011 and 2012 data
	RooCategory sample_year("sample_year","sample_year") ;
	sample_year.defineType("y11");
	sample_year.defineType("y12");

	RooDataSet* combData;
	if(fitSimultan) combData = new RooDataSet("combData","combined data",RooArgSet(DTF_Bs_M),Index(sample_year),Import("y11",*data11),Import("y12",*data12));

	///import Bkg shapes from Normalization fit and MC---------------------------------------------------------------------------------------------------------------------------------------------------------------



	//4) Ds*Kpipi in DsKpipi
	// 13 params: mean1-3, sigma1-6, f_1 & f_2, N_Dstar, N_Bs
	double fillarr4[13];
	double *DsstarKpipifromNorm;
	DsstarKpipifromNorm = fitBDTNorm(fillarr4, sevenTeV, fitSimultan);

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
	RooRealVar sigmaB1("sigmaB1", "B_{0} #sigma_{1}",  32.56,5.,200.);	
	RooRealVar sigmaB2("sigmaB2", "B_{0} sigma_{2}", 13.,10.,50.);
	RooGaussian GaussB01("GaussB01", "GaussB01", DTF_Bs_M, meanB01, sigmaB1);
	RooGaussian GaussB02("GaussB02", "GaussB02", DTF_Bs_M, meanB01, sigmaB2);
	RooRealVar f_GaussB("f_GaussB" , "f__{B_{0}}", 0.066);//, 0.2, 0.75);
        RooAddPdf DoubleGaussB0("DoubleGaussB0", "DoubleGaussB0", RooArgList(GaussB01,GaussB02),RooArgList(f_GaussB));

	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 15.45,10.,90.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 15.,10.,55.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaB1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0.2, 0.75);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussB));

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
	RooRealVar ratioDstar("ratioDstar","ratioDstar",1.,0.,10.);
	RooRealVar N_DstarKpipi_12("N_DstarKpipi_12","N_DstarKpipi 12", 1275, 0, data12->numEntries());
	RooGenericPdf N_DstarKpipiShifted_11("N_DstarKpipiShifted_11","N_DstarKpipiShifted 11", "@0*@1",RooArgList(ratioDstar,N_DstarKpipi_11));
	RooGenericPdf N_DstarKpipiShifted_12("N_DstarKpipiShifted_12","N_DstarKpipiShifted 12", "@0*@1",RooArgList(ratioDstar,N_DstarKpipi_12));

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
	pdf11 = new RooAddPdf("pdf11", "pdf11", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_11, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11, *N_Dstarpipipi_11, N_comb_11, N_DstarKpipiShifted_11));

	pdf12 = new RooAddPdf("pdf12", "pdf12", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_12, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12));

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

	string eleven = "_11.eps";
	string twelve = "_12.eps";

	// 2011 results
	if(sevenTeV){
		BsMassDistribution.append(eleven);
		BsMassResidual.append(eleven);
		BsMassPull.append(eleven);
		BsMassSweight.append(eleven);
	}

	//2012 results
	if(!sevenTeV){
		BsMassDistribution.append(twelve);
		BsMassResidual.append(twelve);
		BsMassPull.append(twelve);
		BsMassSweight.append(twelve);
	}
	
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= DTF_Bs_M.frame();
        RooPlot* frame_m_11= DTF_Bs_M.frame();
        RooPlot* frame_m_12= DTF_Bs_M.frame();
	frame_m->SetTitle("");
        frame_m_11->SetTitle("");
        frame_m_12->SetTitle("");

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
		c1->Print(BsMassDistribution.c_str());
	}

        if(fitSimultan){
                combData->plotOn(frame_m_11,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(60));
                simPdf->plotOn(frame_m_11,Name("pdf11"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_11->Draw();
                c1->Print("eps/Final/BmassFit_sim11.eps");

                combData->plotOn(frame_m_12,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(60));
                simPdf->plotOn(frame_m_12,Name("pdf12"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_12->Draw();
                c1->Print("eps/Final/BmassFit_sim12.eps");
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
                c1->Print(BsMassResidual.c_str());
        }
        if(fitSimultan){
                frame2_11->SetTitle("");
                frame2_11->addPlotable(hresid11,"P") ;
                frame2_11->Draw();
                c1->Print("eps/Final/residual_sim11.eps");

                frame2_12->SetTitle("");
                frame2_12->addPlotable(hresid12,"P") ;
                frame2_12->Draw();
                c1->Print("eps/Final/residual_sim12.eps");
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
                c1->Print(BsMassPull.c_str());
        }
        if(fitSimultan){
                frame3_11->SetTitle("");
                frame3_11->SetLabelFont(62,"Y");
                frame3_11->addPlotable(hpull11,"P") ;
                frame3_11->Draw();
                c1->Print("eps/Final/pull_sim11.eps");

                frame3_12->SetTitle("");
                frame3_12->SetLabelFont(62,"Y");
                frame3_12->addPlotable(hpull12,"P") ;
                frame3_12->Draw();
                c1->Print("eps/Final/pull_sim12.eps");
        }

if(doToys && (!fitSimultan)){

	Int_t SampleEvents;
	//if(sevenTeV) SampleEvents = 10377;
	//if(!sevenTeV) SampleEvents = 22299;
	SampleEvents = 250000;

	//const RooDataSet& toy_data;
	//RooDataSet toy_data; 
	TFile* toy_output;
	if(sevenTeV) toy_output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs_11_final_toy.root","RECREATE");
	if(!sevenTeV) toy_output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_toy.root","RECREATE");

	RooMCStudy* Toystudy = new RooMCStudy(*pdf,RooArgSet(DTF_Bs_M),Binned(kFALSE),Silence(),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
	Toystudy->generateAndFit(1000,SampleEvents,"kTRUE");
	RooDataSet toy_data = Toystudy->fitParDataSet();

	// Make plots of the distributions of mean, the error on mean and the pull of mean
	
	RooPlot* frame_Bs_value = Toystudy->plotParam(N_Bs,Bins(40));
	RooPlot* frame_Bs_error = Toystudy->plotError(N_Bs,Bins(40));
	RooPlot* frame_Bs_pull = Toystudy->plotPull(N_Bs,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_Bs_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_Distribution_12.eps");
 	frame_Bs_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_Error_12.eps");
 	frame_Bs_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_Bs_PullDistribution_12.eps");



	RooPlot* frame_B0_value = Toystudy->plotParam(N_B0,Bins(40));
	RooPlot* frame_B0_error = Toystudy->plotError(N_B0,Bins(40));
	RooPlot* frame_B0_pull = Toystudy->plotPull(N_B0,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_B0_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_B0_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_B0_Distribution_12.eps");
 	frame_B0_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_B0_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_B0_Error_12.eps");
 	frame_B0_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_B0_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_B0_PullDistribution_12.eps");


	
	RooPlot* frame_N_DstarKpipi_value = Toystudy->plotParam(N_DstarKpipi,Bins(40));
	RooPlot* frame_N_DstarKpipi_error = Toystudy->plotError(N_DstarKpipi,Bins(40));
	RooPlot* frame_N_DstarKpipi_pull = Toystudy->plotPull(N_DstarKpipi,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_N_DstarKpipi_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_Distribution_12.eps");
	frame_N_DstarKpipi_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_Error_12.eps");
	frame_N_DstarKpipi_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipi_PullDistribution_12.eps");



	RooPlot* frame_N_comb_value = Toystudy->plotParam(N_comb,Bins(40));
	RooPlot* frame_N_comb_error = Toystudy->plotError(N_comb,Bins(40));
	RooPlot* frame_N_comb_pull = Toystudy->plotPull(N_comb,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_N_comb_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_comb_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_comb_Distribution_12.eps");
	frame_N_comb_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_comb_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_comb_Error_12.eps");
	frame_N_comb_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_comb_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_comb_PullDistribution_12.eps");



	RooPlot* frame_exp_par_value = Toystudy->plotParam(exp_par,Bins(40));
	RooPlot* frame_exp_par_error = Toystudy->plotError(exp_par,Bins(40));
	RooPlot* frame_exp_par_pull = Toystudy->plotPull(exp_par,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_exp_par_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/exp_par_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/exp_par_Distribution_12.eps");
	frame_exp_par_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/exp_par_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/exp_par_Error_12.eps");
	frame_exp_par_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/exp_par_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/exp_par_PullDistribution_12.eps");



	RooPlot* frame_N_DstarKpipiShifted_value = Toystudy->plotParam(N_DstarKpipiShifted,Bins(40));
	RooPlot* frame_N_DstarKpipiShifted_error = Toystudy->plotError(N_DstarKpipiShifted,Bins(40));
	RooPlot* frame_N_DstarKpipiShifted_pull = Toystudy->plotPull(N_DstarKpipiShifted,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_N_DstarKpipiShifted_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_Distribution_12.eps");
	frame_N_DstarKpipiShifted_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_Error_12.eps");
	frame_N_DstarKpipiShifted_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/N_DstarKpipiShifted_PullDistribution_12.eps");



	RooPlot* frame_meanB01_value = Toystudy->plotParam(meanB01,Bins(40));
	RooPlot* frame_meanB01_error = Toystudy->plotError(meanB01,Bins(40));
	RooPlot* frame_meanB01_pull = Toystudy->plotPull(meanB01,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_meanB01_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanB01_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanB01_Distribution_12.eps");
	frame_meanB01_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanB01_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanB01_Error_12.eps");
	frame_meanB01_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanB01_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanB01_PullDistribution_12.eps");



	RooPlot* frame_meanBs1_value = Toystudy->plotParam(meanBs1,Bins(40));
	RooPlot* frame_meanBs1_error = Toystudy->plotError(meanBs1,Bins(40));
	RooPlot* frame_meanBs1_pull = Toystudy->plotPull(meanBs1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_meanBs1_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_Distribution_12.eps");
	frame_meanBs1_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_Error_12.eps");
	frame_meanBs1_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/meanBs1_PullDistribution_12.eps");



	RooPlot* frame_sigmaBs2_value = Toystudy->plotParam(sigmaBs2,Bins(40));
	RooPlot* frame_sigmaBs2_error = Toystudy->plotError(sigmaBs2,Bins(40));
	RooPlot* frame_sigmaBs2_pull = Toystudy->plotPull(sigmaBs2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaBs2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_Distribution_12.eps");
	frame_sigmaBs2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_Error_12.eps");
	frame_sigmaBs2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaBs2_PullDistribution_12.eps");



	RooPlot* frame_sigmaB1_value = Toystudy->plotParam(sigmaB1,Bins(40));
	RooPlot* frame_sigmaB1_error = Toystudy->plotError(sigmaB1,Bins(40));
	RooPlot* frame_sigmaB1_pull = Toystudy->plotPull(sigmaB1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaB1_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_Distribution_12.eps");
	frame_sigmaB1_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_Error_12.eps");
	frame_sigmaB1_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB1_PullDistribution_12.eps");



	RooPlot* frame_sigmaB2_value = Toystudy->plotParam(sigmaB2,Bins(40));
	RooPlot* frame_sigmaB2_error = Toystudy->plotError(sigmaB2,Bins(40));
	RooPlot* frame_sigmaB2_pull = Toystudy->plotPull(sigmaB2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_sigmaB2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_Distribution_12.eps");
	frame_sigmaB2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_Error_12.eps");
	frame_sigmaB2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Signal/sigmaB2_PullDistribution_12.eps");

	toy_data.Write();
	toy_output->Close();
}

	if(sWeight && (!fitSimultan)){

		meanB01.setConstant();
		meanBs1.setConstant();
		sigmaB2.setConstant();
		sigmaBs2.setConstant();
		exp_par.setConstant();
		meanDspipipi1.setConstant();
		meanDspipipi2.setConstant();
	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(N_B0, N_Bs, N_DstarKpipi, *N_Dspipipi, *N_Dstarpipipi, N_comb, N_DstarKpipiShifted)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("DTF_Bs_M,N_Bs_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print(BsMassSweight.c_str());

    		///Create output file
   		 TFile* output;
		 if(sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs_11_final_sweight.root","RECREATE");
		 if(!sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("DTF_Bs_M > 4800 && DTF_Bs_M < 5800");
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


	if(altModel && fitSimultan){


		// set polynomial background
		RooRealVar a1_11("a1_11","a1_11",5.,-100.,100.);
		RooRealVar a2_11("a2_11","a2_11",0.,-10.,10.);
		RooRealVar a3_11("a3_11","a3_11",0.,-1.,1.);
		RooPolynomial polyBG_11("polyBG_11","polyBG_11",DTF_Bs_M,RooArgList(a1_11,a2_11,a3_11),0);

		RooRealVar a1_12("a1_12","a1_12",5.,-100.,100.);
		RooRealVar a2_12("a2_12","a2_12",0.,-10.,10.);
		RooRealVar a3_12("a3_12","a3_12",0.,-1.,1.);
		RooPolynomial polyBG_12("polyBG_12","polyBG_12",DTF_Bs_M,RooArgList(a1_12,a2_12,a3_12),0);

		// vary peaking background within PIDCalib uncertainties 

		double DstarVariation = -0.00315;
		double DsVariation = 0.00276;
		//N_Dstarpipipi_11->setVal((0.036 + DstarVariation) * DsstarKpipifromNorm[11]);
		//N_Dspipipi_11->setVal((0.0335 + DsVariation) * DsstarKpipifromNorm[12]);
		RooRealVar* N_Dspipipi_11_alt = new RooRealVar("N_Dspipipi_11_alt","N_Dspipipi_11_alt",DsstarKpipifromNorm[12],0,10000);
		RooRealVar* N_Dstarpipipi_11_alt = new RooRealVar("N_Dstarpipipi_11_alt","N_Dstarpipipi_11_alt",DsstarKpipifromNorm[11],0,100000);
		//N_Dstarpipipi_12->setVal((0.036 + DstarVariation) * DsstarKpipifromNorm[13]);
		//N_Dspipipi_12->setVal((0.0335 + DsVariation) * DsstarKpipifromNorm[14]);
		RooRealVar* N_Dspipipi_12_alt = new RooRealVar("N_Dspipipi_12_alt","N_Dspipipi_12_alt",DsstarKpipifromNorm[14],0,10000);
		RooRealVar* N_Dstarpipipi_12_alt = new RooRealVar("N_Dstarpipipi_12_alt","N_Dstarpipipi_12_alt",DsstarKpipifromNorm[13],0,100000);


        	RooSimultaneous* simPdf_alt;
        	RooAbsPdf* pdf11_alt;
        	RooAbsPdf* pdf12_alt;

        	simPdf_alt = new RooSimultaneous("simPdf_alt","alternative simultaneous pdf",sample_year);	
		pdf11_alt = new RooAddPdf("pdf11_alt", "pdf11_alt", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_11, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11_alt, *N_Dstarpipipi_11_alt, N_comb_11, N_DstarKpipiShifted_11));

		pdf12_alt = new RooAddPdf("pdf12_alt", "pdf12_alt", RooArgList(DoubleGaussB0, DoubleGaussBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_12, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12_alt, *N_Dstarpipipi_12_alt, N_comb_12, N_DstarKpipiShifted_12));

		simPdf_alt->addPdf(*pdf11_alt,"y11");
		simPdf_alt->addPdf(*pdf12_alt,"y12");

		//fit alt. models
		RooFitResult *result_alt;
		result_alt = simPdf_alt->fitTo(*combData,Save(kTRUE),Extended(kTRUE),NumCPU(3));

		cout << "result is --------------- "<<endl;
		result_alt->Print();

		//draw alt. models
        	RooPlot* frame_m_11_alt= DTF_Bs_M.frame();
        	RooPlot* frame_m_12_alt= DTF_Bs_M.frame();

      	  	frame_m_11_alt->SetTitle("");
       		frame_m_12_alt->SetTitle("");

                combData->plotOn(frame_m_11_alt,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(60));
                simPdf_alt->plotOn(frame_m_11_alt,Name("pdf11_alt"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_11_alt,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_11_alt->Draw();
                c1->Print("eps/altModels/BmassFit_sim11_allBGfloat.eps");

                combData->plotOn(frame_m_12_alt,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(60));
                simPdf_alt->plotOn(frame_m_12_alt,Name("pdf12_alt"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(DoubleGaussB0),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(RooArgSet(DstarKpipi_as_DsKpipi_Shifted)),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(Dspipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(Dstarpipipi_as_DsKpipi),ProjWData(sample_year,*combData),LineColor(kOrange),LineStyle(kDashed),LineWidth(1));
                simPdf_alt->plotOn(frame_m_12_alt,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
                frame_m_12_alt->Draw();
                c1->Print("eps/altModels/BmassFit_sim12_allBGfloat.eps");

	}


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

int quickFit(bool sevenTeV=true ){

	bool sWeight=false;

	///Load file
	TFile* file;
	if(sevenTeV) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_forBDT_tmp.root");
	if(!sevenTeV)file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_forBDT_tmp.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("DTF_Bs_M",1);

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} #pi #pi #pi)", 5200., 5500.,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));

	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 39.27,5.,45.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 11.,5.,45.);
	RooRealVar a1("a1","a1", 1.49,0.,5.);
	RooRealVar n1("n1","n1", 1., 0., 100.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanBs1, sigmaBs2, a1, n1);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0., 1.);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf signal("signal", "signal", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_GaussBs));

	//yields
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());
	RooRealVar N_comb("N_comb","N_comb", data->numEntries()/2., 0., data->numEntries());

	//Exponential
	RooRealVar f_comb("f_comb" , "f_comb", 0.5, 0., 1.);
	RooRealVar exp_par("exp_par","#lambda",0.,-1.,1.);	
	RooExponential bkg("bkg_exp","exponential background",DTF_Bs_M,exp_par);

	//add pdfs
	//RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, GaussBs2, bkg_exp), RooArgList(f_GaussBs, f_comb));
	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, bkg), RooArgList(N_Bs, N_comb));


	///Fit
	RooFitResult *result;
	result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
	cout << "result is --------------- "<<endl;
	result->Print();


	///calculate # (signal)background events in signal region

	DTF_Bs_M.setRange("SigRange",meanBs1.getVal()-50.,meanBs1.getVal()+50.);
	
	RooAbsReal* S_fr= GaussBs1.createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("SigRange"));
	Double_t S = S_fr->getVal()*N_Bs.getVal();
	RooAbsReal* B_fr= bkg.createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("SigRange"));
	Double_t B = B_fr->getVal()*N_comb.getVal();
	
	cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
	cout<<"S/B= " << S/B<< endl;
	cout<<"S= " << S<< endl;
	cout<<"B= " << B<< endl;

	cout << endl;
	cout << endl;

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= DTF_Bs_M.frame();
	frame_m->SetTitle("");

	string BsMassDistribution = "eps/3pi_BmassShape_afterPreSel";
	string BsMassSweight = "eps/afterPresel/3pi_Bs_sWeight";

	string eleven = "_11.eps";
	string twelve = "_12.eps";

	// 2011 results
	if(sevenTeV){
		BsMassDistribution.append(eleven);
		BsMassSweight.append(eleven);
	}

	//2012 results
	if(!sevenTeV){
		BsMassDistribution.append(twelve);
		BsMassSweight.append(twelve);
	}



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
	c1->Print(BsMassDistribution.c_str());

if(sWeight){

		meanBs1.setConstant();
		sigmaBs1.setConstant();
		exp_par.setConstant();

	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(N_Bs, N_comb));
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("DTF_Bs_M,N_Bs_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print(BsMassSweight.c_str());

    		///Create output file
   		 TFile* output;
		 if(sevenTeV) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_afterPreSel_sweight.root","RECREATE");
		 if(!sevenTeV) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_afterPreSel_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("DTF_Bs_M > 5200 && DTF_Bs_M < 5500");
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

return S;

}


void quickSignalEstimate(){

	bool sevenTeV=false;


///Load file
	TFile* file;
	if(sevenTeV) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_forBDT_tmp.root");
	if(!sevenTeV) file=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT_tmp.root");
	//file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_forBDT.root");	
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("DTF_Bs_M",1);

	///Fill all needed variables in RooDataSet

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} K #pi #pi)", 5320., 5420.,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));

	//int expectedYield_12 =1494;
	//int expectedYield_12 =1212;
	//int expectedYield_11 = 506;
	//int expectedYield_11 = 623;
	//int expectedYield = 0.05238 * 0.95122 * quickFit(sevenTeV);
	//int expectedYield = 0.05238 * 0.595 * quickFit(sevenTeV);
	int expectedYield = 0.05238 * 0.863 * quickFit(sevenTeV);

	int expectedYield_Ds2pipipi = expectedYield * 0.2 ;
	int expectedYield_Ds2Kpipi = expectedYield * 0.12 ;

	//int expectedYield_11

	int all = data->numEntries();
	cout << "all Events in Signal Range:  " << all << endl;

	cout<<"expected Signal Yield:  " <<  expectedYield << endl;
	cout<<"Background Yield in Signal Region:  " <<  all - expectedYield << endl;

}

/*
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
	
}*/

void MCStudies(){

///Load files to study peaking bg

//load mc file
TFile* file;
file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/mc2011_Bs2Dspipipi_Ds2KKpi_forBDT_P_ETA_pi_plus1.root");
TTree* tree = (TTree*) file->Get("DecayTree");

TFile* file2;
file2= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/mc2011_Bs2Dspipipi_Ds2KKpi_forBDT_P_ETA_pi_plus2.root");
TTree* tree2 = (TTree*) file2->Get("DecayTree");

//load weight files
TFile* fileUp1_w;
fileUp1_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/PIDEfficiencies_peakBG_pi_plus1_BstoDs3pi_MagUp_11.root");
TTree* treeUp1_w = (TTree*) fileUp1_w->Get("CalibTool_PIDCalibTree");	

TFile* fileUp2_w;
fileUp2_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/PIDEfficiencies_peakBG_pi_plus2_BstoDs3pi_MagUp_11.root");
TTree* treeUp2_w = (TTree*) fileUp2_w->Get("CalibTool_PIDCalibTree");

TFile* fileDown2_w;
fileDown2_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/PIDEfficiencies_peakBG_pi_plus2_BstoDs3pi_MagDown_11.root");
TTree* treeDown2_w = (TTree*) fileDown2_w->Get("CalibTool_PIDCalibTree");

TFile* fileDown1_w;
fileDown1_w= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/PIDEfficiencies_peakBG_pi_plus1_BstoDs3pi_MagDown_11.root");
TTree* treeDown1_w = (TTree*) fileDown1_w->Get("CalibTool_PIDCalibTree");


//the output trees
TFile *Dstar3piBkg_Up1 = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/peakBG/Ds3piBKgShape.root","recreate");
TTree *output_Up1= tree->CloneTree(0);

double EventWeight = 0;
double Bs_Mass = 0;
output_Up1->Branch("EventWeight",&EventWeight,"EventWeight/D");
output_Up1->Branch("Bs_Mass",&Bs_Mass,"Bs_Mass/D");

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
tree -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );

tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );

tree -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
tree -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
tree -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );

//-----------------------------------------------------------------------------

tree2 -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
tree2 -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
tree2 -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );

tree2 -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
tree2 -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
tree2 -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );

tree2 -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
tree2 -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
tree2 -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

tree2 -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
tree2 -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
tree2 -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );

tree2 -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree2 -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree2 -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );

tree2 -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
tree2 -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
tree2 -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );


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
int numEvents2 = tree2->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	treeUp1_w->GetEntry(i);
	treeUp2_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
	
	//flip the mass hypothesis of pion with bigger miss-ID weight
        if(Event_PIDCalibEffWeight_Up1 > Event_PIDCalibEffWeight_Up2){ 
		pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massKaon); //flip mass hypothesis here
        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        	pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
	}
        if(Event_PIDCalibEffWeight_Up2 > Event_PIDCalibEffWeight_Up1){ 
		pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        	pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon); //flip mass hypothesis here
	}

	if((pi_plus1 + pi_minus + pi_plus2).M() > 3000.) continue;

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	if(Event_PIDCalibEffWeight_Up1 > Event_PIDCalibEffWeight_Up2) massBs_Up1->Fill(Bs_MM,Event_PIDCalibEffWeight_Up1);
        if(Event_PIDCalibEffWeight_Up1 < Event_PIDCalibEffWeight_Up2) massBs_Up1->Fill(Bs_MM,Event_PIDCalibEffWeight_Up2);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Up1_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Up1);
	EventWeight = TMath::Max(Event_PIDCalibEffWeight_Up1,Event_PIDCalibEffWeight_Up2);
	Bs_Mass = Bs_MM;
	output_Up1->Fill();
	}


///Up2
for(int i=0; i< numEvents2; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents2 << endl;
        tree2->GetEntry(i);
	treeUp2_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion); 
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon); //flip mass hypothesis here

	if((pi_plus1 + pi_minus + pi_plus2).M() > 3000.) continue;

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Up2->Fill(Bs_MM,Event_PIDCalibEffWeight_Up2);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Up2_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Up2);
	EventWeight = Event_PIDCalibEffWeight_Up2;
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

	if((pi_plus1 + pi_minus + pi_plus2).M() > 3000.) continue;

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Down1->Fill(Bs_MM,Event_PIDCalibEffWeight_Down1);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Down1_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Down1);
	}

///Down2
for(int i=0; i< numEvents2; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents2 << endl;
        tree2->GetEntry(i);
	treeDown2_w->GetEntry(i);

        //define the Lorentz vectors
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion); 
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon); //flip mass hypothesis here

	if((pi_plus1 + pi_minus + pi_plus2).M() > 3000.) continue;

	Bs_MM = (pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs + pi_plus1 + pi_minus + pi_plus2).M();

	massBs_Down2->Fill(Bs_MM,Event_PIDCalibEffWeight_Down2);
	if(Bs_MM > 4800. && Bs_MM < 5800.) massBs_Down2_fitRange->Fill(Bs_MM,Event_PIDCalibEffWeight_Down2);
	}


TCanvas* c= new TCanvas();

double Integral_Up1 = 0;
double Integral_Up1_inRange = 0;
double Integral_Up2 = 0;
double Integral_Up2_inRange = 0;
double Integral_Down1 = 0;
double Integral_Down1_inRange = 0;
double Integral_Down2 = 0;
double Integral_Down2_inRange = 0;

Integral_Up1 = massBs_Up1->Integral();
Integral_Up1_inRange = massBs_Up1_fitRange->Integral();
Integral_Up2 = massBs_Up2->Integral();
Integral_Up2_inRange = massBs_Up2_fitRange->Integral();
Integral_Down1 = massBs_Down1->Integral();
Integral_Down1_inRange = massBs_Down1_fitRange->Integral();
Integral_Down2 = massBs_Down2->Integral();
Integral_Down2_inRange = massBs_Down2_fitRange->Integral();

cout << "Integral of MagUp, Pion1: " << Integral_Up1 << endl;
cout << "Integral of MagUp, Pion1, only (4800 - 5800) MeV: " << Integral_Up1_inRange << endl;
cout << "The ratio is " << (Integral_Up1_inRange/Integral_Up1) << endl;

cout << "Integral of MagUp, Pion2: " << Integral_Up2 << endl;
cout << "Integral of MagUp, Pion2, only (4800 - 5800) MeV: " << Integral_Up2_inRange << endl;
cout << "The ratio is " << (Integral_Up2_inRange/Integral_Up2) << endl;

cout << "Integral of MagDown, Pion1: " << Integral_Down1 << endl;
cout << "Integral of MagDown, Pion1, only (4800 - 5800) MeV: " << Integral_Down1_inRange << endl;
cout << "The ratio is " << (Integral_Down1_inRange/Integral_Down1) << endl;

cout << "Integral of MagDown, Pion2: " << Integral_Down2 << endl;
cout << "Integral of MagDown, Pion2, only (4800 - 5800) MeV: " << Integral_Down2_inRange << endl;
cout << "The ratio is " << (Integral_Down2_inRange/Integral_Down2) << endl;

massBs_Up1->Sumw2(); massBs_Up1->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up1.eps");
massBs_Up1_fitRange->Sumw2(); massBs_Up1_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up1_inRange.eps");

massBs_Up2->Sumw2(); massBs_Up2->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up2.eps");
massBs_Up2_fitRange->Sumw2(); massBs_Up2_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Up2_inRange.eps");

massBs_Down1->Sumw2(); massBs_Down1->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down1.eps");
massBs_Down1_fitRange->Sumw2(); massBs_Down1_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down1_inRange.eps");

massBs_Down2->Sumw2(); massBs_Down2->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down2.eps");
massBs_Down2_fitRange->Sumw2(); massBs_Down2_fitRange->Draw("e1"); c->Print("eps/BkgShape/PIDCalib_forDs3pi_Down2_inRange.eps");

output_Up1->Write();
Dstar3piBkg_Up1->Close();
}

void splitMC(){

TFile* file;
file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
TTree* tree = (TTree*) file->Get("DecayTree");

int N = tree->GetEntries();
cout << "full MC file contains " << N << " events" <<  endl;

Double_t K_plus_fromDs_PX;
Double_t K_plus_fromDs_PY;
Double_t K_plus_fromDs_PZ;

Double_t K_minus_fromDs_PX;
Double_t K_minus_fromDs_PY;
Double_t K_minus_fromDs_PZ;

Double_t pi_minus_fromDs_PX;
Double_t pi_minus_fromDs_PY;
Double_t pi_minus_fromDs_PZ;

TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;

double massKaon = 493.68;
double massPion = 139.57;
double massPhi = 1019.46;
double massKstar = 895.81;

tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );

tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );

tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );


TFile* output_PhiResonant = 0;
TFile* output_KstarResonant = 0;
TFile* output_nonResonant = 0;

output_PhiResonant = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Ds2KKpi_Bs2Dspipipi_forBDT_PhiResonant.root","RECREATE");
TTree* phi_tree = tree->CloneTree(0);
int numEvents = tree->GetEntries();


for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

        //define the Lorentz vectors
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);

	if(TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) < 20) phi_tree->Fill();
	}

cout << "New file for Phi resonant MC contains " << phi_tree->GetEntries() << " events" <<  endl;
phi_tree->Write();
output_PhiResonant->Close();



output_KstarResonant = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Ds2KKpi_Bs2Dspipipi_forBDT_KstarResonant.root","RECREATE");
TTree* Kstar_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

        //define the Lorentz vectors
	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);

	if(TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) < 75) Kstar_tree->Fill();
	}

cout << "New file for K* resonant MC contains " << Kstar_tree->GetEntries() << " events" <<  endl;
Kstar_tree->Write();
output_KstarResonant->Close();



output_nonResonant = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Ds2KKpi_Bs2Dspipipi_forBDT_nonResonant.root","RECREATE");
TTree* nonRes_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

        //define the Lorentz vectors
	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);

	if((TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20) && (TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) > 75)) nonRes_tree->Fill();
	}

cout << "New file for non-resonant MC contains " << nonRes_tree->GetEntries() << " events" <<  endl;
nonRes_tree->Write();
output_nonResonant->Close();

}

void DataStructure(bool sevenTeV, bool SignalMode){

	TFile* file;
	if(sevenTeV && SignalMode) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs_11_final_sweight.root");
	if((!sevenTeV) && SignalMode) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root");
	if(sevenTeV && (!SignalMode)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_final_sweight.root");
	if((!sevenTeV) && (!SignalMode)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("DTF_Bs_M",1);
	tree->SetBranchStatus("N_Bs_sw",1);
	tree->SetBranchStatus( "*PX", 1 );
	tree->SetBranchStatus( "*PY", 1);
	tree->SetBranchStatus( "*PZ", 1);

	double massKaon = 493.68;
	double massPion = 139.57;
	double massPhi = 1019.46;
	double massKstar = 895.81;

	Double_t K_plus_fromDs_PX;
	Double_t K_plus_fromDs_PY;
	Double_t K_plus_fromDs_PZ;
	Double_t K_minus_fromDs_PX;
	Double_t K_minus_fromDs_PY;
	Double_t K_minus_fromDs_PZ;
	Double_t pi_minus_fromDs_PX;
	Double_t pi_minus_fromDs_PY;
	Double_t pi_minus_fromDs_PZ;
	Double_t N_Bs_sw;

	TLorentzVector K_plus_fromDs;
	TLorentzVector K_minus_fromDs;
	TLorentzVector pi_minus_fromDs;

	TH1D* resonant_structure = new TH1D("resonant structure", ";resonance category;Counts", 3, -0.5, 2.5);

        tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
        tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
        tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
        tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
        tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
        tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
        tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
        tree -> SetBranchAddress( "N_Bs_sw" , &N_Bs_sw );


int numEvents = tree->GetEntries();
int resonance = 0;

	for(int i=0; i< numEvents; i++)
        {
        	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        	tree->GetEntry(i);

		pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
		K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);


                resonance = 0;
                if(TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) < 20) resonance = 1;
                else if(TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) < 75) {
                        resonance = 2;
                }
                else{
                        resonance = 0;
		}

	resonant_structure->Fill(resonance,N_Bs_sw);

	}
	TCanvas* c = new TCanvas();

	resonant_structure->Draw("E1");
	c->Print("eps/resonant_structure_S11.eps");

	cout << "PhiPi: " << resonant_structure->GetBinContent(2) << " p/m " << resonant_structure->GetBinError(2) << endl;
	cout << "KstarK: " << resonant_structure->GetBinContent(3) << " p/m " << resonant_structure->GetBinError(3) << endl;
	cout << "nonRes: " << resonant_structure->GetBinContent(1) << " p/m " << resonant_structure->GetBinError(1) << endl;
}

int main(){
    time_t startTime = time(0);

 // quickSignalEstimate();
  //  iterateBDT(-0.9,0.9,0.05);
    //preselect();
 // applyBDTcut("0.7012");
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
   //addCut();
  // fitBDT();
  // quickFit(true);
  //  MCStudies();
   //splitMC();
/// bool sevenTeV , bool signalMode
  //DataStructure(true, true);

	bool sevenTeV=false;
	bool fitSimultan=false;
	bool doToys=true;
	double fillarr4[13];
	double *DsstarKpipifromNorm;
	DsstarKpipifromNorm = fitBDTNorm(fillarr4,sevenTeV,fitSimultan,doToys,false);

    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;

return 0;
}

