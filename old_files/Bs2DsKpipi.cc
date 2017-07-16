
//main file for Bs->DsKpipi Analysis
//philippe d'argent & matthieu kecke

#include <boost/lexical_cast.hpp>
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

using namespace std;
using namespace RooFit ;
using namespace RooStats;

int preselect(bool MC=false , bool preselection =false) {

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
    //tree=new TChain("Bs2Dspipipi_Ds2pipipi_Tuple/DecayTree");
tree=new TChain("DecayTree");
  //  tree->Add("/auto/data/kecke/B2DKPiPi/12D-3pi/*.root");
   // tree->Add("/auto/data/kecke/B2DKPiPi/12U-3pi/*.root");
      //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-U/*.root");
      //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-D/*.root");
      tree->Add("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_forBDT.root");
    }

    int N = tree->GetEntries();
    cout << "Old file contains " << N << " events" <<  endl;
    
    //Disable all branches
    tree->SetBranchStatus("*",1);
        
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
	if(BsDspipipi)output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_forBDT_tmp.root","RECREATE");
		
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
    
    TTree* tree_sel = tree->CopyTree(/*cuts.c_str()*/"(Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) > 1.5");
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

TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT.root");
TTree* tree = (TTree*) file->Get("DecayTree");


TFile* output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT_DsK3fb_Selection.root","RECREATE");
TTree* summary_tree = tree->CloneTree(0);

   double massKaon = 493.68;
   double massPion = 139.57;
   double massProton = 938.27;

   double massDs = 1968.30;
   double massDminus = 1869.61;

   double massLambda_c = 2286.46;
   double massLambda_b = 5619.51;

   Double_t pi_plus1_PX;
   Double_t pi_plus1_PY;
   Double_t pi_plus1_PZ;
   Double_t pi_minus_PX;
   Double_t pi_minus_PY;
   Double_t pi_minus_PZ;
   Double_t pi_plus2_PX;
   Double_t pi_plus2_PY;
   Double_t pi_plus2_PZ;
   Double_t Ds_PX;
   Double_t Ds_PY;
   Double_t Ds_PZ;

   Double_t Ds_ENDVERTEX_Z;
   Double_t Bs_ENDVERTEX_Z;
   Double_t Ds_FDCHI2_ORIVX;

   Double_t K_plus_fromDs_PX;
   Double_t K_plus_fromDs_PY;
   Double_t K_plus_fromDs_PZ;
   Double_t K_minus_fromDs_PX;
   Double_t K_minus_fromDs_PY;
   Double_t K_minus_fromDs_PZ;
   Double_t pi_minus_fromDs_PX;
   Double_t pi_minus_fromDs_PY;
   Double_t pi_minus_fromDs_PZ;

   Double_t K_minus_fromDs_PIDK;
   Double_t K_minus_fromDs_PIDp;

  //look for D_s1
  Double_t DTF_Bs_M;

TLorentzVector pi_plus1;
TLorentzVector pi_plus2;
TLorentzVector pi_minus;
TLorentzVector Ds;
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;
TLorentzVector Kminus_asProton_MissID;
TLorentzVector Kminus_asPiminus_MissID;

TLorentzVector Lambda_c;

/*
    tree->SetBranchAddress("pi_plus1_PX",&pi_plus1_PX);
    tree->SetBranchAddress("pi_plus1_PY",&pi_plus1_PY);
    tree->SetBranchAddress("pi_plus1_PZ",&pi_plus1_PZ);
    tree->SetBranchAddress("pi_plus2_PX",&pi_plus2_PX);
    tree->SetBranchAddress("pi_plus2_PY",&pi_plus2_PY);
    tree->SetBranchAddress("pi_plus2_PZ",&pi_plus2_PZ);
    tree->SetBranchAddress("pi_minus_PX",&pi_minus_PX);
    tree->SetBranchAddress("pi_minus_PY",&pi_minus_PY);
    tree->SetBranchAddress("pi_minus_PZ",&pi_minus_PZ);
*/
    tree->SetBranchAddress("Ds_ENDVERTEX_Z",&Ds_ENDVERTEX_Z);
    tree->SetBranchAddress("Bs_ENDVERTEX_Z",&Bs_ENDVERTEX_Z);
    tree->SetBranchAddress("Ds_FDCHI2_ORIVX",&Ds_FDCHI2_ORIVX);

    tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
    tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
    tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
    tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
    tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
    tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
    tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
    tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
    tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

    tree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
    tree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );



    tree->SetBranchAddress("Ds_PX",&Ds_PX);
    tree->SetBranchAddress("Ds_PY",&Ds_PY);
    tree->SetBranchAddress("Ds_PZ",&Ds_PZ);

    tree->SetBranchAddress("DTF_Bs_M",&DTF_Bs_M);

    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		tree->GetEntry(i);
/*
		pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        	pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
*/
		pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);

		Kminus_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);
		Kminus_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);

		Ds.SetXYZM(Ds_PX,Ds_PY,Ds_PZ,massDs);
		Lambda_c = K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MissID;

		///vetos from Bs->DsK 3 fb^-1 analysis
		//if(TMath::Abs((Ds + pi_minus + pi_plus1).M() - 2535.11) < 25) continue;
		//if(TMath::Abs((Ds + pi_minus + pi_plus2).M() - 2535.11) < 25) continue;
		//D veto
		if(TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 30. && (K_minus_fromDs_PIDK < 10) ) continue;
		//D0 veto
		if((K_plus_fromDs + K_minus_fromDs).M() > 1840.) continue;
		// Lc veto
        	if(TMath::Abs((K_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 30. && (K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) continue;
		// Ds vertex seperation
		if(Ds_FDCHI2_ORIVX < 2) continue;
		//if(TMath::Abs((Lambda_c + pi_minus + pi_plus1 + pi_plus2).M() - massLambda_b) < 15) continue;
		//if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 1.) continue;
		summary_tree->Fill();
    }
summary_tree->Write();
output->Close();

   //TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/Bs2Dspipipi_fullSelectionBDTG_tightDCUT_DZ_Ds1Veto.root","RECREATE");
   //TTree* new_tree_BDTG = tree->CopyTree("Ds_MM > 1950 && Ds_MM < 1990");
   //TTree* new_tree_BDTG = tree->CopyTree("K_plus_PIDK > 8");
   //TTree* new_tree_BDTG = tree->CopyTree("(Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) > 3");
   //TTree* new_tree_BDTG = tree->CopyTree("(K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > -1");
   //new_tree_BDTG->Write();
   //output_BDTG->Close();


   //close file at the end
   file->Close();

}

void addVars(){

TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");
TTree* tree = (TTree*) file->Get("DecayTree");

/*
Double_t Bs_TAU;
Double_t Bs_TAUERR;
*/
Double_t Bs_TRUETAU;

/*
tree -> SetBranchAddress("Bs_TAU" , &Bs_TAU );
tree -> SetBranchAddress("Bs_TAUERR" , &Bs_TAUERR );
*/

tree -> SetBranchAddress("Bs_TRUETAU" , &Bs_TRUETAU );
/*
Int_t Bs_TRUEID;
Int_t K_plus_TRUEID;
*/
//tree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
//tree -> SetBranchAddress( "K_plus_TRUEID" , &K_plus_TRUEID );


TFile* output=new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined_wTrueCT.root","RECREATE");
TTree* new_tree = tree->CloneTree();


//double Bs_ct = 0;
double Bs_cttrue = 0;
//double Bs_cterr = 0;

//int qt = 0;
//int qf = 0;

//TBranch* Bs_ct_Branch = new_tree->Branch("Bs_ct",&Bs_ct,"Bs_ct/D");
TBranch* Bs_cttrue_Branch = new_tree->Branch("Bs_cttrue",&Bs_cttrue,"Bs_cttrue/D");
//TBranch* Bs_cterr_Branch = new_tree->Branch("Bs_cterr",&Bs_cterr,"Bs_cterr/D");

//TBranch* qt_Branch = new_tree->Branch("qt",&qt,"qt/I");
//TBranch* qf_Branch = new_tree->Branch("qf",&qf,"qf/I");

int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

        //Bs_ct = Bs_TAU*1000;
        Bs_cttrue = Bs_TRUETAU*1000;
        //Bs_cterr = Bs_TAUERR*1000;
	
	//Bs_ct_Branch->Fill();
	Bs_cttrue_Branch->Fill();
	//Bs_cterr_Branch->Fill();
/*
	if(Bs_TRUEID > 0) qt = 1;
	if(Bs_TRUEID < 0) qt = -1;
	if(K_plus_TRUEID > 0) qf = 1;
	if(K_plus_TRUEID < 0) qf = -1;

	qt_Branch->Fill();
	qf_Branch->Fill();
*/
}


new_tree->Write();
output->Close();

file->Close();


}

void checkForDuplicates(){

TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/12/b2dhhh_285.root");
TTree* tree = (TTree*) file->Get("Bs2Dspipipi_Ds2KKpi_Tuple/DecayTree");

//TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_BDT.root");
//TTree* tree = (TTree*) file->Get("DecayTree");

cout << "check 1" << endl;

ULong64_t eventNumber;
UInt_t runNumber;
Double_t Bs_MM;
Double_t Ds_MM;
Double_t Bs_P;
Double_t Ds_P;


tree -> SetBranchAddress("eventNumber" , &eventNumber );
tree -> SetBranchAddress("runNumber" , &runNumber );
tree -> SetBranchAddress("Bs_MM" , &Bs_MM );
tree -> SetBranchAddress("Ds_MM" , &Ds_MM );
tree -> SetBranchAddress("Bs_P" , &Bs_P );
tree -> SetBranchAddress("Ds_P" , &Ds_P );

cout << "check 2" << endl;

int numEvents = tree->GetEntries();

long allNumbers[numEvents];
int  allRuns[numEvents];

cout << "check 3" << endl;

long double Bs_mass[numEvents];

cout << "check 3" << endl;

long double Bs_momentum[numEvents];

cout << "check 3" << endl;

long double Ds_mass[numEvents];

cout << "check 3" << endl;

long double Ds_momentum[numEvents];



for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	allNumbers[i] = eventNumber;
        allRuns[i] = runNumber;
        Bs_mass[i] = Bs_MM;
        Bs_momentum[i] = Bs_P;
        Ds_mass[i] = Ds_MM;
	Ds_momentum[i] = Ds_P;
        }

/*
for(int i=0; i< 500; i++)
        {
	
        cout << "event number in event " << i << "  is: " << allNumbers[i] << endl;
        cout << "run number number in event  " << i << "  is: " << allRuns[i] << endl;

	}
*/
	
  // Check for duplicate numbers in user inputted data
  // the kecke way

int m; // Need to declare m here so that it can be accessed by the 'inner' loop 
for(m = 0;m < numEvents; m++) { // Check each other number in the array
    for(int n = m; n < numEvents; n++) { // Check the rest of the numbers
      if(n != m) { // Makes sure don't check number against itself
        if((allNumbers[m]==allNumbers[n])&&(allRuns[m] == allRuns[n])&&(Bs_mass[m] == Bs_mass[n])&&(Ds_mass[m]==Ds_mass[n])&&(Bs_momentum[m]==Bs_momentum[n])&&(Ds_momentum[m]==Ds_momentum[n])){
	      	
		cout << "following event number is present twice.  " << m +1 << ". Duplicate numbers are not cool." << endl;
	}    
      }
                
    } // Comparison loop
     // Reset the boolean after each number entered has been checked
} // Main check loop  


    //the elegant way
//int* end = allNumbers + numEvents;
//std::sort(allNumbers, numEvents);
//bool containsDuplicates = (std::unique(allNumbers, numEvents) != numEvents);
 

file->Close();


}


void applyBDTcut(string cutoff){
   ///Load file

//  TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data12_Ds2KKpi_BDT.root");
//  TTree* tree = (TTree*) file->Get("DecayTree");	
TFile* file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_BDT.root");
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
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG_DsK3fb_Selection.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();

*/
    //BDTG application for 2012 Signal Channel (9500BG/450S), 0.8040 for maximum S/sqrt(S+B) , with S = 726 , B = 18751 
/*
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();
*/

    //BDTG application for 2012 normalization Channel (9500BG/450S) , 0.1458 for maximum S/sqrt(S+B) with S = 13867 , B = 23588
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data12_Bs2Dspipipi_fullSelectionBDTG.root","RECREATE");
   TTree* new_tree_BDTG = tree->CopyTree(cstringcutstring);
   new_tree_BDTG->Write();
   output_BDTG->Close();

/*
    //BDTG application for 2011 normalization Channel (32.500BG/6500S) ------> BDTG > 0.5571 with S=14511 and B=98053
   TFile* output_BDTG=new TFile("/auto/data/kecke/B2DPiPiPi/forMaster/Data2011/data11_Bs2Dspipipi_fullSelectionBDTG.root","RECREATE");
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


double *fitBDTNorm(double fitvalues[15], bool sevenTeV, bool fitSimultan, bool Ds2KKpi, bool combineYears, bool doToys=false, bool altModel = false){

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

	string filename11="/auto/data/kecke/B2DPiPiPi/Data2011/data11_Bs2Dspipipi_fullSelectionBDTG_DsK3fb_Selection.root";
	string filename12="/auto/data/kecke/B2DPiPiPi/Data2012/data12_Bs2Dspipipi_fullSelectionBDTG_DsK3fb_Selection.root";

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
	RooRealVar a1("a1","a1", -1.8098e+00);//,-10.,10.);
	RooRealVar n1("n1","n1", 5.5870e+00);//,0.,100.);
	RooRealVar a2("a2","a2", 9.7197e-01);
	RooRealVar n2("n2","n2", 1.7471e+00);
	RooRealVar sigmaCB1("sigmaCB1", "B_{s} CB #sigma_{1}", 1.3120e+01,0.,50.56); //32.56
	RooRealVar sigmaCB2("sigmaCB2", "B_{s} CB #sigma_{2}", 10.,0.,20.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanBs1, sigmaCB1, a1, n1);
	RooCBShape CB2("CB2", "CB2", DTF_Bs_M, meanBs1, sigmaCB2, a2, n2);
	RooRealVar f_CB("f_CB" , "f__{B_{s}}", 7.0763e-01);//,0.,1.);
	RooAddPdf DoubleCBs("DoubleCBs", "DoubleCBs", RooArgList(CB1,CB2),RooArgList(f_CB));
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.18);//, 0.,0.5);
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_CB));
	RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));

	//signal pdf
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());

	///Background model - exponential for combinatorial Bkg + fixed BG shape from peaking Bkg

	//Exponential for first comb. BG. component
	RooRealVar exp_par("exp_par","#lambda_{1}",-1.6508e-03,-10.,0.);	
	RooRealVar exp_par_11("exp_par_11","#lambda_{1} 11", getExpBkgShape(filename11.c_str())/*-1.6508e-03,-10.,0.*/);
	RooRealVar exp_par_12("exp_par_12","#lambda_{1} 12", getExpBkgShape(filename12.c_str())/*-1.6508e-03,-10.,0.*/);
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
	RooRealVar* mean1;
	RooRealVar* mean2;
	RooRealVar* mean3;

	///width and fractions of gaussians 2011

	if(sevenTeV){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);// sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);
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

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);
	}


	if(fitSimultan && Ds2KKpi && (!combineYears)){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 2*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);//, sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8], sigmaR3Variation_low, sigmaR3Variation_high);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0], mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1], mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2], mean3Variation_low, mean3Variation_high);
	}

	if(fitSimultan && (!Ds2KKpi) && (!combineYears)){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 3*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7], sigmaL3Variation_low, sigmaL3Variation_high);
		sigmaR3 = new RooRealVar("sigma_{3R}", "sigmaR3",  DstarpipipiNorm[8]);//, sigmaR3Variation_low, sigmaR3Variation_high);

		mean1 = new RooRealVar("mean1","mu", DstarpipipiNorm[0]);//, mean1Variation_low, mean1Variation_high);
		mean2 = new RooRealVar("mean2","mu", DstarpipipiNorm[1]);//, mean2Variation_low, mean2Variation_high);
		mean3 = new RooRealVar("mean3","mu", DstarpipipiNorm[2]);//, mean3Variation_low, mean3Variation_high);

		f_1 = new RooRealVar("f_{1}", "fraction1", DstarpipipiNorm[9], 0.,1.);
		f_2 = new RooRealVar("f_{2}", "fraction2", DstarpipipiNorm[10], 0.,1.);
	}

	if(combineYears){
		sigmaL1 = new RooRealVar("sigma_{1L}", "sigmaL1",  DstarpipipiNorm[3], sigmaL1Variation_low, 3*sigmaL1Variation_high);
		sigmaR1 = new RooRealVar("sigma_{1R}", "sigmaR1",  DstarpipipiNorm[4]);//, sigmaR1Variation_low, sigmaR1Variation_high);
		sigmaL2 = new RooRealVar("sigma_{2L}", "sigmaL2",  DstarpipipiNorm[5], sigmaL2Variation_low, sigmaL2Variation_high);
		sigmaR2 = new RooRealVar("sigma_{2R}", "sigmaR2",  DstarpipipiNorm[6], sigmaR2Variation_low, sigmaR2Variation_high);
		sigmaL3 = new RooRealVar("sigma_{3L}", "sigmaL3",  DstarpipipiNorm[7]);//, sigmaL3Variation_low, sigmaL3Variation_high);
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
		pdf11=new RooAddPdf("pdf11", "pdf11", RooArgList(/*DoubleGaussBs*/DoubleCBs, Dstarpipipi_as_Dspipipi ,bkg_exp_11), RooArgList(N_Bs_11, N_Dstarpipipi_11, N_comb_11));
		pdf12=new RooAddPdf("pdf12", "pdf12", RooArgList(/*DoubleGaussBs*/DoubleCBs, Dstarpipipi_as_Dspipipi ,bkg_exp_12), RooArgList(N_Bs_12, N_Dstarpipipi_12, N_comb_12));
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
	string simfit11 = "_sim11.pdf";
	string simfit12 = "_sim12.pdf";
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
		c1->Print(BsMassDistribution.c_str());
	}

	if((!BDTscan) && fitSimultan && (!combineYears)){
		combData->plotOn(frame_m_11,Name("data11"),Cut("sample_year==sample_year::y11"),MarkerSize(0.5),Binning(120));
		simPdf->plotOn(frame_m_11,Name("pdf11"),Slice(sample_year,"y11"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		//simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(DoubleCBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_year,"y11"),Components(bkg_exp_11),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11->Draw();
		c1->Print(BsMassDistribution11.c_str());

		combData->plotOn(frame_m_12,Name("data12"),Cut("sample_year==sample_year::y12"),MarkerSize(0.5),Binning(120));
		simPdf->plotOn(frame_m_12,Name("pdf12"),Slice(sample_year,"y12"),ProjWData(sample_year,*combData),LineColor(kBlack),LineWidth(2));
		//simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleGaussBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(DoubleCBs),ProjWData(sample_year,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_year,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_year,"y12"),Components(bkg_exp_12),ProjWData(sample_year,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12->Draw();
		c1->Print(BsMassDistribution12.c_str());
	}

	if((!BDTscan) && fitSimultan && combineYears){
		combData->plotOn(frame_m_11,Name("dataDs2KKpi"),Cut("sample_Ds==sample_Ds::Ds_kaonkaonpion"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_11,Name("pdfDs2KKpi"),Slice(sample_Ds,"Ds_kaonkaonpion"),ProjWData(sample_Ds,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(DoubleGaussBs),ProjWData(sample_Ds,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_Ds,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_11,Slice(sample_Ds,"Ds_kaonkaonpion"),Components(bkg_exp_11),ProjWData(sample_Ds,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_11->Draw();
		c1->Print(BsMassDistribution11.c_str());

		combData->plotOn(frame_m_12,Name("dataDs2pipipi"),Cut("sample_Ds==sample_Ds::Ds_pionpionpion"),MarkerSize(0.5),Binning(60));
		simPdf->plotOn(frame_m_12,Name("pdfDs2pipipi"),Slice(sample_Ds,"Ds_pionpionpion"),ProjWData(sample_Ds,*combData),LineColor(kBlack),LineWidth(2));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(DoubleGaussBs),ProjWData(sample_Ds,*combData),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(Dstarpipipi_as_Dspipipi),ProjWData(sample_Ds,*combData),LineColor(kMagenta),LineStyle(kDashed),LineWidth(1));
		simPdf->plotOn(frame_m_12,Slice(sample_Ds,"Ds_pionpionpion"),Components(bkg_exp_12),ProjWData(sample_Ds,*combData),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		frame_m_12->Draw();
		c1->Print(BsMassDistribution12.c_str());
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
		c1->Print(BsMassResidual.c_str());
	}
	if((!BDTscan) && fitSimultan){
		frame2_11->SetTitle("");
		frame2_11->addPlotable(hresid11,"P") ;
		frame2_11->Draw();
		c1->Print(BsMassResidual11.c_str());

		frame2_12->SetTitle("");
		frame2_12->addPlotable(hresid12,"P") ;
		frame2_12->Draw();
		c1->Print(BsMassResidual12.c_str());
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
		frame3_11->SetFillColor(kBlue);
		frame3_11->Draw("hist");
		c1->Print(BsMassPull11.c_str());

		frame3_12->SetTitle("");	
		frame3_12->SetLabelFont(62,"Y");
		frame3_12->addPlotable(hpull12,"P") ;
		frame3_12->SetFillColor(kBlue);
		frame3_12->Draw("hist");
		c1->Print(BsMassPull12.c_str());
	}


	///save fitavalues in array to be use in fitBDT()
	fitvalues[0] = mean1->getVal();
	fitvalues[1] = mean2->getVal();
	fitvalues[2] = mean3->getVal();
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
	if(fitSimultan && (!combineYears)){
		fitvalues[11] = N_Dstarpipipi_11.getVal();
		fitvalues[12] = N_Bs_11.getVal();
		fitvalues[13] = N_Dstarpipipi_12.getVal();
		fitvalues[14] = N_Bs_12.getVal();
	}
	if(fitSimultan && combineYears){
		fitvalues[11] = N_Dstarpipipi_Ds2KKpi.getVal();
		fitvalues[12] = N_Bs_Ds2KKpi.getVal();
		fitvalues[13] = N_Dstarpipipi_Ds2pipipi.getVal();
		fitvalues[14] = N_Bs_Ds2pipipi.getVal();
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



	RooPlot* frame_mean1_value = Toystudy->plotParam(*mean1,Bins(40));
	RooPlot* frame_mean1_error = Toystudy->plotError(*mean1,Bins(40));
	RooPlot* frame_mean1_pull = Toystudy->plotPull(*mean1,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_mean1_value->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Distribution_12.eps");
	frame_mean1_error->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_Error_12.eps");
	frame_mean1_pull->Draw("E1");
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean1_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean1_PullDistribution_12.eps");



	RooPlot* frame_mean2_value = Toystudy->plotParam(*mean2,Bins(40));
	RooPlot* frame_mean2_error = Toystudy->plotError(*mean2,Bins(40));
	RooPlot* frame_mean2_pull = Toystudy->plotPull(*mean2,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

	frame_mean2_value->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Distribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Distribution_12.eps");
	frame_mean2_error->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Error_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_Error_12.eps");
	frame_mean2_pull->Draw("E1"); 
	if(sevenTeV) c1->Print("eps/Toys/Normalization/mean2_PullDistribution_11.eps");
	if(!sevenTeV) c1->Print("eps/Toys/Normalization/mean2_PullDistribution_12.eps");



	RooPlot* frame_mean3_value = Toystudy->plotParam(*mean3,Bins(40));
	RooPlot* frame_mean3_error = Toystudy->plotError(*mean3,Bins(40));
	RooPlot* frame_mean3_pull = Toystudy->plotPull(*mean3,Bins(30),FrameRange(-5., 5.),FitGauss(kTRUE));

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


if(sWeight/* && (!fitSimultan)*/){

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

	
		SPlot* sData = new SPlot("sData","An SPlot",*data,pdf,RooArgList(N_Bs, N_Dstarpipipi, N_comb)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data->createHistogram("DTF_Bs_M,N_Bs_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print(BsMassSweight11.c_str());

    		///Create output file
   		 TFile* output;
		 if(sevenTeV && Ds2KKpi) output = new TFile("/auto/data/kecke/B2DPiPiPi/forMaster/Data2011/data_Bs2Dspipipi_11_final_sweight.root","RECREATE");
		// if(!sevenTeV && Ds2KKpi) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root","RECREATE");
		// if(sevenTeV && (!Ds2KKpi)) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_Ds2pipipi_11_final_sweight.root","RECREATE");
		// if(!sevenTeV && (!Ds2KKpi)) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_Ds2pipipi_12_final_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
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

		mean1->setVal(DstarpipipiNorm[0]);
		mean2->setVal(DstarpipipiNorm[1]);
		mean3->setVal(DstarpipipiNorm[2]);
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
		mean1->setConstant();
		mean2->setConstant();
		mean3->setConstant();

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


void createRooDataSet(string fileName,bool sevenTeV = true){

bool data = false;
bool SigMode = false;


string root = ".root";

string Bs2Dspipipi = "B2DPiPiPi/";
string Bs2DsKpipi = "B2DKPiPi/";

string eleven = "Data2011/";
string twelve = "Data2012/";
string eleven_mc = "MC2011/";
string twelve_mc = "MC2012/";

string autoDataLink = "/auto/data/kecke/";
string newLink = "/auto/data/kecke/";

// Bs2DsKpipi data
if(sevenTeV && SigMode && data){
autoDataLink.append(Bs2DsKpipi);
autoDataLink.append(eleven);
newLink.append(Bs2DsKpipi);
newLink.append(eleven);
}
if((!sevenTeV) && SigMode && data){
autoDataLink.append(Bs2DsKpipi);
autoDataLink.append(twelve);
newLink.append(Bs2DsKpipi);
newLink.append(twelve);
}

//Bs2Dspipipi data
if((!sevenTeV) && (!SigMode) && data){
autoDataLink.append(Bs2Dspipipi);
autoDataLink.append(twelve);
newLink.append(Bs2Dspipipi);
newLink.append(twelve);
}

if((sevenTeV) && (!SigMode) && data){
autoDataLink.append(Bs2Dspipipi);
autoDataLink.append(eleven);
newLink.append(Bs2Dspipipi);
newLink.append(eleven);
}

//Bs2Dspipipi mc
if((!sevenTeV) && (!SigMode) && (!data)){
autoDataLink.append(Bs2Dspipipi);
autoDataLink.append(twelve_mc);
newLink.append(Bs2Dspipipi);
newLink.append(twelve_mc);
}
if((sevenTeV) && (!SigMode) && (!data)){
autoDataLink.append(Bs2Dspipipi);
autoDataLink.append(eleven_mc);
newLink.append(Bs2Dspipipi);
newLink.append(eleven_mc);
}

// Bs2DsKpipi mc 
if(sevenTeV && SigMode && (!data)){
autoDataLink.append(Bs2DsKpipi);
autoDataLink.append(eleven_mc);
newLink.append(Bs2DsKpipi);
newLink.append(eleven_mc);
}
if((!sevenTeV) && SigMode && (!data)){
autoDataLink.append(Bs2DsKpipi);
autoDataLink.append(twelve_mc);
newLink.append(Bs2DsKpipi);
newLink.append(twelve_mc);
}


string completeLink = autoDataLink.append(fileName);
completeLink.append(root);

// load the file to be converted to RooDataSet
TFile* file = new TFile(completeLink.c_str());
TTree* tree = (TTree*) file->Get("DecayTree");

string saveAs = "_asRooDataSet.root";
string newFileName = fileName.append(saveAs);

string saveLink = newLink.append(newFileName);

// define RooRealVar for observables
RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar DTF_TAU("DTF_TAU", "DTF_TAU", 0., 15., "ps");
RooRealVar Bs_BsDTF_ctauErr("Bs_BsDTF_ctauErr", "Bs_BsDTF_ctauErr", 0.,15.,"ps");
RooRealVar N_Bs_sw("N_Bs_sw", "N_Bs_sw", -0.3, 1.2);
RooRealVar BDTG_response("BDTG_response", "BDTG_response", 0., 1.);

RooArgSet observables;

if(data) observables = RooArgSet(DTF_Bs_M, DTF_TAU, N_Bs_sw, BDTG_response);
if(!data) observables =RooArgSet(DTF_Bs_M, DTF_TAU,Bs_BsDTF_ctauErr);

RooDataSet *dataset = new RooDataSet("dataset","dataset",tree, observables);

TFile* output = new TFile(saveLink.c_str(),"RECREATE");
RooWorkspace* workspace = new RooWorkspace("workspace", "workspace");


workspace->import(*dataset);

workspace->Write();
output->Close();

}

void fitBDT(){

	bool sWeight=false;
	bool sevenTeV=false; // 2011 = true, 2012 = false
	bool doToys=false;
	bool fitSimultan=true;
	bool altModel=false;
	bool Ds2KKpi = false;
        bool combineYears = false;

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

	string filename11="/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_fullSelectionBDTG_DsK3fb_Selection.root";
	string filename12="/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_fullSelectionBDTG_DsK3fb_Selection.root";

	
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
	double fillarr4[13];
	double *DsstarKpipifromNorm;
	DsstarKpipifromNorm = fitBDTNorm(fillarr4, sevenTeV, fitSimultan, Ds2KKpi, combineYears);

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
        string simfit11 = "_sim11.pdf";
        string simfit12 = "_sim12.pdf";
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
		c1->Print(BsMassDistribution.c_str());
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
                c1->Print(BsMassDistribution11.c_str());

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
                c1->Print(BsMassDistribution12.c_str());
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
                c1->Print(BsMassResidual11.c_str());

                frame2_12->SetTitle("");
                frame2_12->addPlotable(hresid12,"P") ;
                frame2_12->Draw();
                c1->Print(BsMassResidual12.c_str());
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
                c1->Print(BsMassPull11.c_str());

                frame3_12->SetTitle("");
                frame3_12->SetLabelFont(62,"Y");
                frame3_12->addPlotable(hpull12,"P") ;
                frame3_12->Draw();
                c1->Print(BsMassPull12.c_str());
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

		pdf11_clone = new RooAddPdf("pdf11_clone", "pdf11_clone", RooArgList(DoubleGaussB0, /*DoubleGaussBs*/DoubleCBBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_11, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_11, N_Bs_11, N_DstarKpipi_11, *N_Dspipipi_11, *N_Dstarpipipi_11, N_comb_11, N_DstarKpipiShifted_11_var));

		pdf12_clone = new RooAddPdf("pdf12_clone", "pdf12_clone", RooArgList(DoubleGaussB0, /*DoubleGaussBs*/DoubleCBBs, DstarKpipi_as_DsKpipi, Dspipipi_as_DsKpipi, Dstarpipipi_as_DsKpipi ,bkg_exp_12, DstarKpipi_as_DsKpipi_Shifted), RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12_var));



		SPlot* sData = new SPlot("sData","An SPlot",*data12,pdf12_clone,RooArgList(N_B0_12, N_Bs_12, N_DstarKpipi_12, *N_Dspipipi_12, *N_Dstarpipipi_12, N_comb_12, N_DstarKpipiShifted_12_var)); 
		gStyle->SetOptStat(0);
	
		///Plot the sWeight distributions as a function of mass
		TCanvas* SwDs = new TCanvas("Bs sWeight","Bs sWeight distribution");
		TH2 * SwDsHist = (TH2*)data12->createHistogram("DTF_Bs_M,N_Bs_12_sw");
		SwDsHist->GetYaxis()->SetTitle("Signal sWeights");
		SwDsHist->SetTitle("");
		//SwDs->Write();
		SwDsHist->Draw();
		SwDs->Print(BsMassSweight12.c_str());

    		///Create output file
   		 TFile* output;
		 //if(sevenTeV && Ds2KKpi) 
		output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root","RECREATE");
		// if(!sevenTeV && Ds2KKpi) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root","RECREATE");
		/*if(sevenTeV && (!Ds2KKpi))*/ //output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_Ds2pipipi_12_final_sweight.root","RECREATE");
		 //if(!sevenTeV && (!Ds2KKpi)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_Ds2pipipi_12_final_sweight.root","RECREATE");

		 tree12->SetBranchStatus("*",1);
   		 TTree* new_tree = tree12->CopyTree("DTF_Bs_M > 4975 && DTF_Bs_M < 5800");
    		 double w;
    		 TBranch* Bra_sw = new_tree->Branch("N_Bs_sw", &w, "N_Bs_sw/D");

  		  ///loop over events
    		  int numEvents = new_tree->GetEntries();

    		  for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree12->GetEntry(i);
			w=sData->GetSWeight(i,"N_Bs_12_sw");
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
    //TFile* fileMC;
  // fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_Reco14.root");
   //file= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_XdReweighted_Reco14.root");
    //file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT_tmp.root");
    file= new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");
   //fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_withTriggerWeights_Reco14.root");
    TTree* tree = (TTree*) file->Get("DecayTree");
  //  TTree* treeMC = (TTree*) fileMC->Get("DecayTree");

   // TH1D* massBs = new TH1D("mass of B_{s}^{0} candidates",";m(D_{s}^{+} K^{-} #pi^{+} #pi^{-}) [MeV];Entries",75, 4750., 6000.);
   // TH1D* massBssyst_L0 = new TH1D("mass of B_{s} candidates L0 syst",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",75, 4750., 6000.);
   // TH1D* massBssyst_HLT1 = new TH1D("mass of B_{s} candidates HLT1 syst",";m(D_{s}^{+} K #pi^{+} #pi^{-}) [MeV];Entries",75, 4750., 6000.);
  //  TH1D* massBs = new TH1D("mass of B_{s} candidates",";m(D_{s}^{+} #pi^{-} #pi^{+} #pi^{-}) [MeV];Entries",40, 5000., 5600.);
   // TH1D* massDs = new TH1D("mass of D_{s} candidates",";m(K^{-} K^{+} #pi^{+}) [MeV];Entries",40, 1900., 2050.);

    //TH2D* nTracksVsGhostProb = new TH2D("nTracks vs max(ghostProb)",";nTracks;max(ghostProb)",10, 50., 300.,10,0.,0.1);
   //plot("nTracks","# of tracks",40,0,450,false,true);
  // plot("max_ghostProb","max(Track_ghostProb)",25,0.,0.375,false,true);

    ///pairs
    TH1D* massDs_Pi1 = new TH1D("mass of D_{s} + #pi 1 candidates",";m(D_{s}^{-} #pi^{+}) [MeV];Entries",75, 1500., 5500.);
    TH1D* massDs_Pi2 = new TH1D("mass of D_{s} + #pi 2 candidates",";m(D_{s}^{-} #pi^{+}) [MeV];Entries",75, 1500., 5500.);
    TH1D* massPiPi1 = new TH1D("mass of #pi + #pi first candidates",";m(#pi^{+} #pi^{-}) [MeV];Entries",75, 200., 2500.);
    TH1D* massPiPi2 = new TH1D("mass of #pi + #pi second candidates",";m(#pi^{+} #pi^{-}) [MeV];Entries",75, 200., 2500.);

   ///triplets
   TH1D* massDs_PiPi1 = new TH1D("mass of D_{s} + #pi +#pi 1 candidates",";m(D_{s}^{-} #pi^{+} #pi^{-}) [MeV];Entries",75, 1500., 6000.);
   TH1D* massDs_PiPi2 = new TH1D("mass of D_{s} + #pi +#pi 2 candidates",";m(D_{s}^{-} #pi^{-} #pi^{+}) [MeV];Entries",75, 1500., 6000.);
   TH1D* massDs_PiPi3 = new TH1D("mass of D_{s} + #pi +#pi 3 candidates",";m(D_{s}^{-} #pi^{+} #pi^{+}) [MeV];Entries",75, 1500., 6000.);
   TH1D* mass_PiPiPi = new TH1D("mass of #pi + #pi +#pi candidates",";m(#pi^{-} #pi^{+} #pi^{+}) [MeV];Entries",75, 500., 4000.);

   ///Bs->D_s1 + pi
   TH1D* mass_D_s1_signal = new TH1D("#Delta m of D_{s1} candidates",";m(D_{s}^{-}#pi^{+}#pi^{-}) - m(D_{s}^{-}) [MeV];Entries",25, 450., 700.);
   /// D_s1-> Ds + pi + pi 
   TH1D* mass_D_s1_peakBG1 = new TH1D("mass of D_{s1} candidates",";m(D_{s}^{-}#pi^{+}#pi^{-}) - m(D_{s}^{-}) [MeV];Entries",50, 2000., 3000.);
   TH1D* mass_D_s1_peakBG2 = new TH1D("mass of D_{s1} candidates ",";m(D_{s}^{-}#pi^{+}#pi^{-}) - m(D_{s}^{-}) [MeV];Entries",50, 2000., 3000.);

   /// Lambda_c -> K + P + pi miss ID? 
   TH1D* mass_Lambda_c_peakBG = new TH1D("mass of #Lambda_{c} candidates",";m(K P #pi) [MeV];Entries",50, 1800., 3000.);
   TH1D* mass_Lambda_b_to_Lambda_c_peakBG = new TH1D("mass of #Lambda_{b} candidates",";m(#Lambda_{c} #pi #pi #pi) [MeV];Entries",50, 5200., 6000.);

   TH1D* cterr = new TH1D("distribution of #sigma(t) (MC)",";#sigma(t) [fs];Entries",30, 0., 150.);

    double mB;
    int nTracks;
    float max_ghostProb;
    //double w;
    double w_L0;
    double w_HLT1;

   double massKaon = 493.68;
   double massPion = 139.57;
   double massDs = 1968.30;
   double massProton = 938.27;


   Double_t pi_plus1_PX;
   Double_t pi_plus1_PY;
   Double_t pi_plus1_PZ;
   Double_t pi_minus_PX;
   Double_t pi_minus_PY;
   Double_t pi_minus_PZ;
   Double_t pi_plus2_PX;
   Double_t pi_plus2_PY;
   Double_t pi_plus2_PZ;
   Double_t Ds_PX;
   Double_t Ds_PY;
   Double_t Ds_PZ;

   Double_t K_plus_fromDs_PX;
   Double_t K_plus_fromDs_PY;
   Double_t K_plus_fromDs_PZ;
   Double_t K_minus_fromDs_PX;
   Double_t K_minus_fromDs_PY;
   Double_t K_minus_fromDs_PZ;
   Double_t pi_minus_fromDs_PX;
   Double_t pi_minus_fromDs_PY;
   Double_t pi_minus_fromDs_PZ;

   Double_t Bs_cterr;

  //look for D_s1
  Double_t DTF_Bs_M;

TLorentzVector pi_plus1;
TLorentzVector pi_plus2;
TLorentzVector pi_minus;
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;
TLorentzVector Kminus_asProton_MissID;
TLorentzVector Ds;
TLorentzVector Lambda_c;

    tree->SetBranchAddress("Bs_cterr",&Bs_cterr);

    tree->SetBranchAddress("pi_plus1_PX",&pi_plus1_PX);
    tree->SetBranchAddress("pi_plus1_PY",&pi_plus1_PY);
    tree->SetBranchAddress("pi_plus1_PZ",&pi_plus1_PZ);
    tree->SetBranchAddress("pi_plus2_PX",&pi_plus2_PX);
    tree->SetBranchAddress("pi_plus2_PY",&pi_plus2_PY);
    tree->SetBranchAddress("pi_plus2_PZ",&pi_plus2_PZ);
    tree->SetBranchAddress("pi_minus_PX",&pi_minus_PX);
    tree->SetBranchAddress("pi_minus_PY",&pi_minus_PY);
    tree->SetBranchAddress("pi_minus_PZ",&pi_minus_PZ);

    tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
    tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
    tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
    tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
    tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
    tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
    tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
    tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
    tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

    tree->SetBranchAddress("Ds_PX",&Ds_PX);
    tree->SetBranchAddress("Ds_PY",&Ds_PY);
    tree->SetBranchAddress("Ds_PZ",&Ds_PZ);

    tree->SetBranchAddress("DTF_Bs_M",&DTF_Bs_M);

    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		tree->GetEntry(i);

		cterr->Fill(Bs_cterr*1000);

		pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        	pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
		Ds.SetXYZM(Ds_PX,Ds_PY,Ds_PZ,massDs);

        	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
		pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		Kminus_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);

		Lambda_c = K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MissID;


    		massDs_Pi1->Fill((Ds + pi_plus1).M());
   		massDs_Pi2->Fill((Ds + pi_plus2).M());
    		massPiPi1->Fill((pi_minus + pi_plus1).M());
    		massPiPi2->Fill((pi_minus + pi_plus2).M());
   		massDs_PiPi1->Fill((Ds + pi_minus + pi_plus1).M());
   		massDs_PiPi2->Fill((Ds + pi_minus + pi_plus2).M());
   		massDs_PiPi3->Fill((Ds + pi_plus2 + pi_plus1).M());
   		mass_PiPiPi->Fill((pi_minus + pi_plus1 + pi_plus2).M());

		//look for Bs->Ds1 + pi signal
		if(DTF_Bs_M > 5327. && DTF_Bs_M < 5407.){
		mass_D_s1_signal->Fill(((Ds + pi_minus + pi_plus1).M() - Ds.M()));
		}

		//look for D_s1 in the uper mass sideband of Bs
		//if(DTF_Bs_M > 5367.){
			mass_D_s1_peakBG1->Fill((Ds + pi_minus + pi_plus1).M());
			mass_D_s1_peakBG2->Fill((Ds + pi_minus + pi_plus2).M());
		//}

		//Lambda_c and Lambda_b->Lambda_c K pi pi 
		mass_Lambda_c_peakBG->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MissID).M());
		mass_Lambda_b_to_Lambda_c_peakBG->Fill((Lambda_c + pi_minus + pi_plus1 + pi_plus2).M());

    }
   TCanvas* c= new TCanvas();
   // pairs
   cterr->Draw("e1"); c->Print("eps/Bs_cterr_MC.pdf");
//   massDs_Pi1->Draw("e1"); c->Print("eps/CombBG_DsPi1.eps");
//   massDs_Pi2->Draw("e1"); c->Print("eps/CombBG_DsPi2.eps");
//   massPiPi1->Draw("e1"); c->Print("eps/CombBG_PiPi1.eps");
//   massPiPi2->Draw("e1"); c->Print("eps/CombBG_PiPi2.eps");

   // triplets
//   massDs_PiPi1->Draw("e1"); c->Print("eps/CombBG_DsPiPi1.eps");
//   massDs_PiPi2->Draw("e1"); c->Print("eps/CombBG_DsPiPi2.eps");
//   massDs_PiPi3->Draw("e1"); c->Print("eps/CombBG_DsPiPi3.eps");
//   mass_PiPiPi->Draw("e1");  c->Print("eps/CombBG_PiPiPi.eps");

   //D_s1 studies
//   mass_D_s1_peakBG1->Draw("E1"); c->Print("eps/CombBG_Ds1_1.eps");
//   mass_D_s1_peakBG2->Draw("E1"); c->Print("eps/CombBG_Ds1_2.eps");

//   mass_D_s1_signal->Draw("E1"); c->Print("eps/Final/DeltaMass_Ds1.eps");
//   mass_Lambda_c_peakBG->Draw("E1"); c->Print("eps/CombBG_Lambda_c.eps");
//   mass_Lambda_b_to_Lambda_c_peakBG->Draw("E1"); c->Print("eps/CombBG_Lambda_b.eps");


/*
    treeMC->SetBranchAddress("DTF_Bs_M",&mB);
    treeMC->SetBranchAddress("weight",&w);
    treeMC->SetBranchAddress("weight_L0",&w_L0);
    treeMC->SetBranchAddress("weight_HLT1",&w_HLT1);
    treeMC->SetBranchAddress("nTracks",&nTracks);
    treeMC->SetBranchAddress("max_ghostProb",&max_ghostProb);

    ///loop over events
    int numEvents = treeMC->GetEntries();
    for(int i=0; i< numEvents; i++){	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		treeMC->GetEntry(i);

		massBs->Fill(mB);
		massBssyst_L0->Fill(mB,w_L0);
		massBssyst_HLT1->Fill(mB,w_HLT1);
		//nTracksVsGhostProb->Fill(nTracks,max_ghostProb,w);

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
		}
    }
*/
 //   massBs->Draw("e1"); c->Print("eps/mass_Bs_fullSel_12.eps");
  //  massDs->Draw("e1"); c->Print("eps/mass_Ds_fullSel_12.eps");


   /// nTracksVsGhostProb->SetStats(0);
   /// nTracksVsGhostProb->Draw("COLZ"); 
   /// c->Print("DataVsMC/nTracksVsGhostProb.eps"); 
   /// c->Print("DataVsMC/nTracksVsGhostProb.pdf");


    //massBs->SetStats(0);
    //massBs->Draw("e1"); c->Print("eps/afterPresel/mass_Bs_forBDT_12.png");
    //c->Print("eps/afterPresel/mass_Bs_forBDT_12.root");
   // massBssyst_L0->Draw("e1"); c->Print("eps/Systematics/mass_Bs_systL0_12.eps");

	//cout << "integral without weights is: " << massBs->Integral() << endl;
	//cout << "integral with L0 weight is: " << massBssyst_L0->Integral() << endl;
	//cout << "integral with HLT1 weight is: " << massBssyst_HLT1->Integral() << endl;
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

void makeSweightedPlots(){

//TFile* file = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight_withCT.root");
//TFile* file = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_DsK3fb_Selection.root");
TFile* file = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_with_BDT_variables_Reco14_LTUB.root");
TTree* tree = (TTree*) file->Get("DecayTree");



RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar Ds_ct("Ds_ct", "Ds_ct", 0., 1.,"ps");
//RooRealVar Bs_ct("Bs_ct", "Bs_ct", 0.15, 10., "ps");
//RooRealVar Bs_cterr("Bs_cterr", "Bs_cterr", 0.0001, 1.,"ps");
//RooRealVar N_Bs_sw("N_Bs_sw", "N_Bs_sw", -1., 2.);

RooArgSet observables( Ds_ct/*, Bs_cterr, N_Bs_sw*/);

RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree)/*, WeightVar(N_Bs_sw.GetName())*/);

TCanvas* canvas_kecko = new TCanvas("canvas", "canvas", 1200, 800);
RooPlot* frame_kecke = Ds_ct.frame();
frame_kecke->SetTitle("");
dataset->plotOn(frame_kecke, Binning(50));
frame_kecke->Draw();
canvas_kecko->Print("eps/kecko_testPlot.eps");

file->Close();
}


int quickFit(bool sevenTeV=false, bool Ds2KKpi=true){

	bool sWeight = false;
	bool Systematics = false;
	bool Simulation = true;

	bool L0 = true;
	bool HLT1 = false;
	bool HLT2 = false;

	if(Systematics) sWeight = false;

	///Load file
	TFile* file;
	if(sevenTeV && Ds2KKpi && (!Simulation)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_forBDT_tmp_DsK3fb_Selection.root");
	if(sevenTeV && (!Ds2KKpi) && (!Simulation)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2pipipi_forBDT.root");
	if((!sevenTeV) && Ds2KKpi && (!Simulation))file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_forBDT_tmp_DsK3fb_Selection.root");
	if((!sevenTeV) && (!Ds2KKpi) && (!Simulation))file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_forBDT.root");
	//if(Simulation) file= new TFile("/auto/data/kecke/B2DPiPiPi/Bs2Dspipipi_MC_fullSel_reweighted_combined.root");
	if(Simulation) file= new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("DTF_Bs_M",1);

	///Load files for TISTOS data studies
	/// load L0 efficiencies
	//TISTOS files
	TFile* fileData_1_TISTOS;
	if(L0) fileData_1_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_pt_1525.root");
	else if(HLT1) fileData_1_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_1525.root");
	TTree* treeData_1_TISTOS = (TTree*) fileData_1_TISTOS->Get("DecayTree");
   	treeData_1_TISTOS->SetBranchStatus("*",0);
	treeData_1_TISTOS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_2_TISTOS;
	if(L0) fileData_2_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_pt_2535.root");
	else if(HLT1) fileData_2_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_2535.root");
	TTree* treeData_2_TISTOS = (TTree*) fileData_2_TISTOS->Get("DecayTree");
   	treeData_2_TISTOS->SetBranchStatus("*",0);
	treeData_2_TISTOS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_3_TISTOS;
	if(L0) fileData_3_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_pt_3545.root");
	else if(HLT1) fileData_3_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_3545.root");
	TTree* treeData_3_TISTOS = (TTree*) fileData_3_TISTOS->Get("DecayTree");
   	treeData_3_TISTOS->SetBranchStatus("*",0);
	treeData_3_TISTOS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_4_TISTOS;
	if(L0) fileData_4_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_pt_4555.root");
	else if(HLT1) fileData_4_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_4555.root");
	TTree* treeData_4_TISTOS = (TTree*) fileData_4_TISTOS->Get("DecayTree");
   	treeData_4_TISTOS->SetBranchStatus("*",0);
	treeData_4_TISTOS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_5_TISTOS;
	if(L0) fileData_5_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_alt_pt_5565.root");
	else if(HLT1) fileData_5_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_5565.root");
	TTree* treeData_5_TISTOS = (TTree*) fileData_5_TISTOS->Get("DecayTree");
   	treeData_5_TISTOS->SetBranchStatus("*",0);
	treeData_5_TISTOS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_6_TISTOS;
	if(L0) fileData_6_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TISTOS_pt_6575.root");
	else if(HLT1) fileData_6_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TISTOS_pt_6575.root");
	TTree* treeData_6_TISTOS = (TTree*) fileData_6_TISTOS->Get("DecayTree");
   	treeData_6_TISTOS->SetBranchStatus("*",0);
	treeData_6_TISTOS->SetBranchStatus("Bs_MM",1);


	//TIS only files
	TFile* fileData_1_TIS;
	if(L0) fileData_1_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_1525.root");
	else if(HLT1) fileData_1_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_1525.root");
	TTree* treeData_1_TIS = (TTree*) fileData_1_TIS->Get("DecayTree");
   	treeData_1_TIS->SetBranchStatus("*",0);
	treeData_1_TIS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_2_TIS;
	if(L0) fileData_2_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_2535.root");
	else if(HLT1) fileData_2_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_2535.root");
	TTree* treeData_2_TIS = (TTree*) fileData_2_TIS->Get("DecayTree");
   	treeData_2_TIS->SetBranchStatus("*",0);
	treeData_2_TIS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_3_TIS;
	if(L0) fileData_3_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_3545.root");
	else if(HLT1) fileData_3_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_3545.root");
	TTree* treeData_3_TIS = (TTree*) fileData_3_TIS->Get("DecayTree");
   	treeData_3_TIS->SetBranchStatus("*",0);
	treeData_3_TIS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_4_TIS;
	if(L0) fileData_4_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_4555.root");
	else if(HLT1) fileData_4_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_4555.root");
	TTree* treeData_4_TIS = (TTree*) fileData_4_TIS->Get("DecayTree");
   	treeData_4_TIS->SetBranchStatus("*",0);
	treeData_4_TIS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_5_TIS;
	if(L0) fileData_5_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_5565.root");
	else if(HLT1) fileData_5_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_5565.root");
	TTree* treeData_5_TIS = (TTree*) fileData_5_TIS->Get("DecayTree");
   	treeData_5_TIS->SetBranchStatus("*",0);
	treeData_5_TIS->SetBranchStatus("Bs_MM",1);

	TFile* fileData_6_TIS;
	if(L0) fileData_6_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_6575.root");
	else if(HLT1) fileData_6_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_HLT1_TIS_pt_6575.root");
	TTree* treeData_6_TIS = (TTree*) fileData_6_TIS->Get("DecayTree");
   	treeData_6_TIS->SetBranchStatus("*",0);
	treeData_6_TIS->SetBranchStatus("Bs_MM",1);

	

	///Fill all needed variables in RooDataSet
  	///---------------------------------------

	///Bs
        RooRealVar DTF_Bs_M("DTF_Bs_M", "m(D_{s} K #pi #pi)", 5300., 5450.,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",RooArgSet(DTF_Bs_M),Import(*tree));

	///for TISTOS
	RooRealVar Bs_MM("Bs_MM", "m(D_{s} #pi #pi #pi)", 5200., 5500.,"MeV/c^{2}");
	RooArgList list_TISTOS =  RooArgList(Bs_MM);
	RooDataSet* data_TISTOS_1;
	RooDataSet* data_TISTOS_2;
	RooDataSet* data_TISTOS_3;
	RooDataSet* data_TISTOS_4;
	RooDataSet* data_TISTOS_5;
	RooDataSet* data_TISTOS_6;
	RooDataSet* data_TIS_1;
	RooDataSet* data_TIS_2;
	RooDataSet* data_TIS_3;
	RooDataSet* data_TIS_4;
	RooDataSet* data_TIS_5;
	RooDataSet* data_TIS_6;

	if(Systematics){
		data_TISTOS_1 = new RooDataSet("data_TISTOS_1","data_TISTOS_1",RooArgSet(Bs_MM),Import(*treeData_1_TISTOS));
		data_TISTOS_2 = new RooDataSet("data_TISTOS_2","data_TISTOS_2",RooArgSet(Bs_MM),Import(*treeData_2_TISTOS));
		data_TISTOS_3 = new RooDataSet("data_TISTOS_3","data_TISTOS_3",RooArgSet(Bs_MM),Import(*treeData_3_TISTOS));
		data_TISTOS_4 = new RooDataSet("data_TISTOS_4","data_TISTOS_4",RooArgSet(Bs_MM),Import(*treeData_4_TISTOS));
		data_TISTOS_5 = new RooDataSet("data_TISTOS_5","data_TISTOS_5",RooArgSet(Bs_MM),Import(*treeData_5_TISTOS));
		data_TISTOS_6 = new RooDataSet("data_TISTOS_6","data_TISTOS_6",RooArgSet(Bs_MM),Import(*treeData_6_TISTOS));

		data_TIS_1 = new RooDataSet("data_TIS_1","data_TIS_1",RooArgSet(Bs_MM),Import(*treeData_1_TIS));
		data_TIS_2 = new RooDataSet("data_TIS_2","data_TIS_2",RooArgSet(Bs_MM),Import(*treeData_2_TIS));
		data_TIS_3 = new RooDataSet("data_TIS_3","data_TIS_3",RooArgSet(Bs_MM),Import(*treeData_3_TIS));
		data_TIS_4 = new RooDataSet("data_TIS_4","data_TIS_4",RooArgSet(Bs_MM),Import(*treeData_4_TIS));
		data_TIS_5 = new RooDataSet("data_TIS_5","data_TIS_5",RooArgSet(Bs_MM),Import(*treeData_5_TIS));
		data_TIS_6 = new RooDataSet("data_TIS_6","data_TIS_6",RooArgSet(Bs_MM),Import(*treeData_6_TIS));
	}

	double N_TISTOS_1 = 0;
	double N_TISTOS_2 = 0;
	double N_TISTOS_3 = 0;
	double N_TISTOS_4 = 0;
	double N_TISTOS_5 = 0;
	double N_TISTOS_6 = 0;

	double N_TIS_1 = 0;
	double N_TIS_2 = 0;
	double N_TIS_3 = 0;
	double N_TIS_4 = 0;
	double N_TIS_5 = 0;
	double N_TIS_6 = 0;

	///put together the roofit pdf
	//Bs signal shape
	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 5370.,5320.,5420.); 
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 39.27,5.,45.);	
	RooRealVar sigmaBs2("sigmaBs2", "B_{s} sigma_{2}", 11.,5.,45.);
	RooRealVar a1("a1","a1", -1.49,-5.,0.);
	RooRealVar n1("n1","n1", 1., 0., 100.);
	RooRealVar a2("a2","a2", 1.49,0.,5.);
	RooRealVar n2("n2","n2", 1., 0., 100.);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", DTF_Bs_M, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", DTF_Bs_M, meanBs1, sigmaBs2);
	RooGaussian GaussBs("GaussBs", "GaussBs", Bs_MM, meanBs1, sigmaBs1);
	RooCBShape CB1("CB1", "CB1", DTF_Bs_M, meanBs1, sigmaBs1, a1, n1);
	RooCBShape CB2("CB2", "CB2", DTF_Bs_M, meanBs1, sigmaBs2, a2, n2);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0., 1.);
        RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf signal("signal", "signal", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
        RooAddPdf GaussCBBs("GaussCBBs", "GaussCBBs", RooArgList(GaussBs1,CB1),RooArgList(f_GaussBs));
        RooAddPdf DoubleCBs("DoubleCBs", "DoubleCBs", RooArgList(CB1,CB2),RooArgList(f_GaussBs));

	//yields
	RooRealVar N_Bs("N_Bs", "#B_{s}", data->numEntries()/2., 0., data->numEntries());
	RooRealVar N_comb("N_comb","N_comb", data->numEntries()/2., 0., data->numEntries());

	//Exponential
	RooRealVar f_comb("f_comb" , "f_comb", 0.5, 0., 1.);
	RooRealVar exp_par("exp_par","#lambda",0.,-1.,0.);	
	RooExponential bkg("bkg_exp","exponential background",DTF_Bs_M,exp_par);
	RooExponential exp_bkg("exp_bkg","exponential background function",Bs_MM,exp_par);

	//add pdfs
	//RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, GaussBs2, bkg_exp), RooArgList(f_GaussBs, f_comb));
	RooAbsPdf* pdf=new RooAddPdf("pdf", "pdf", RooArgList(GaussBs1, bkg), RooArgList(N_Bs, N_comb));

	RooAbsPdf* pdf_CB=new RooAddPdf("pdf_CB", "pdf_CB", RooArgList(DoubleCBs, bkg), RooArgList(f_comb));

	RooAbsPdf* pdf_TISTOS = new RooAddPdf("pdf_TISTOS", "pdf_TISTOS", RooArgList(GaussBs, exp_bkg), RooArgList(N_Bs, N_comb));


	///Fit
	RooFitResult *result;
	if(!Systematics){
		if(!Simulation) result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(3));
		if(Simulation) result = DoubleCBs.fitTo(*data,Save(kTRUE),Extended(kFALSE),NumCPU(3));
		cout << "result is --------------- "<<endl;
		result->Print();
	}

	///Fit TIS & TISTOS datasets
	RooFitResult *result_TISTOS_1;
	RooFitResult *result_TISTOS_2;
	RooFitResult *result_TISTOS_3;
	RooFitResult *result_TISTOS_4;
	RooFitResult *result_TISTOS_5;
	RooFitResult *result_TISTOS_6;

	RooFitResult *result_TIS_1;
	RooFitResult *result_TIS_2;
	RooFitResult *result_TIS_3;
	RooFitResult *result_TIS_4;
	RooFitResult *result_TIS_5;
	RooFitResult *result_TIS_6;

	if(Systematics){
		result_TISTOS_1 = pdf_TISTOS->fitTo(*data_TISTOS_1,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_1 Bin is --------------- "<<endl;
		result_TISTOS_1->Print();
		N_TISTOS_1 = N_Bs.getVal();
	

		result_TISTOS_2 = pdf_TISTOS->fitTo(*data_TISTOS_2,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_2 Bin is --------------- "<<endl;
		result_TISTOS_2->Print();
		N_TISTOS_2 = N_Bs.getVal();

		result_TISTOS_3 = pdf_TISTOS->fitTo(*data_TISTOS_3,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_3 Bin is --------------- "<<endl;
		result_TISTOS_3->Print();
		N_TISTOS_3 = N_Bs.getVal();

		result_TISTOS_4 = pdf_TISTOS->fitTo(*data_TISTOS_4,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_4 Bin is --------------- "<<endl;
		result_TISTOS_4->Print();
		N_TISTOS_4 = N_Bs.getVal();

		result_TISTOS_5 = pdf_TISTOS->fitTo(*data_TISTOS_5,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_5 Bin is --------------- "<<endl;
		result_TISTOS_5->Print();
		N_TISTOS_5 = N_Bs.getVal();

		result_TISTOS_6 = pdf_TISTOS->fitTo(*data_TISTOS_6,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TISTOS_6 Bin is --------------- "<<endl;
		result_TISTOS_6->Print();
		N_TISTOS_6 = N_Bs.getVal();



		result_TIS_1 = pdf_TISTOS->fitTo(*data_TIS_1,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_1 Bin is --------------- "<<endl;
		result_TIS_1->Print();
		N_TIS_1 = N_Bs.getVal();


		result_TIS_2 = pdf_TISTOS->fitTo(*data_TIS_2,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_2 Bin is --------------- "<<endl;
		result_TIS_2->Print();
		N_TIS_2 = N_Bs.getVal();


		result_TIS_3 = pdf_TISTOS->fitTo(*data_TIS_3,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_3 Bin is --------------- "<<endl;
		result_TIS_3->Print();
		N_TIS_3 = N_Bs.getVal();


		result_TIS_4 = pdf_TISTOS->fitTo(*data_TIS_4,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_4 Bin is --------------- "<<endl;
		result_TIS_4->Print();
		N_TIS_4 = N_Bs.getVal();


		result_TIS_5 = pdf_TISTOS->fitTo(*data_TIS_5,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_5 Bin is --------------- "<<endl;
		result_TIS_5->Print();
		N_TIS_5 = N_Bs.getVal();


		result_TIS_6 = pdf_TISTOS->fitTo(*data_TIS_6,Save(kTRUE),Extended(kTRUE),NumCPU(2));
		cout << "result for TIS_6 Bin is --------------- "<<endl;
		result_TIS_6->Print();
		N_TIS_6 = N_Bs.getVal();
	}

	///calculate # (signal)background events in signal region

	DTF_Bs_M.setRange("SigRange",meanBs1.getVal()-50.,meanBs1.getVal()+50.);
	
	RooAbsReal* S_fr= GaussBs1.createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("SigRange"));
	Double_t S = S_fr->getVal()*N_Bs.getVal();
	RooAbsReal* B_fr= bkg.createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("SigRange"));
	Double_t B = B_fr->getVal()*N_comb.getVal();
	
	if(!Systematics){
		cout<< "S/sqrt(S+B)= " << S/sqrt(S+B) << endl;
		cout<<"S/B= " << S/B<< endl;
		cout<<"S= " << S<< endl;
		cout<<"B= " << B<< endl;

		cout << endl;
		cout << endl;
	}

	///Plot 
	///----------
	TCanvas* c1= new TCanvas("");
	RooPlot* frame_m= DTF_Bs_M.frame();
	frame_m->SetTitle("");

	string BsMassDistribution = "eps/3pi_BmassShape_afterPreSel";
	string BsMassSweight = "eps/afterPresel/3pi_Bs_sWeight";
	string Dsto3pi = "_Ds23pi";

	string eleven = "_11.eps";
	string twelve = "_12.eps";

	if(!Systematics){
		// 2011 results
		if(sevenTeV && Ds2KKpi){
			BsMassDistribution.append(eleven);
			BsMassSweight.append(eleven);
		}

		//2012 results
		if((!sevenTeV) && Ds2KKpi){
			BsMassDistribution.append(twelve);
			BsMassSweight.append(twelve);
		}

		if((!sevenTeV) && (!Ds2KKpi)){
			BsMassDistribution.append(Dsto3pi);
			BsMassSweight.append(Dsto3pi);
			BsMassDistribution.append(twelve);
			BsMassSweight.append(twelve);
		}
		if(sevenTeV && (!Ds2KKpi)){
			BsMassDistribution.append(Dsto3pi);
			BsMassSweight.append(Dsto3pi);
			BsMassDistribution.append(eleven);
			BsMassSweight.append(eleven);
		}


		data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
		DoubleCBs.plotOn(frame_m,Name(""),LineColor(kBlack),LineWidth(2));
		DoubleCBs.plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		DoubleCBs.plotOn(frame_m,Components(CB2),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		//pdf->plotOn(frame_m,Components(GaussBs1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		//pdf->plotOn(frame_m,LineColor(kBlack),LineWidth(1));
		//pdf->plotOn(frame_m,Components(CB1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
		//pdf->plotOn(frame_m,Components(bkg),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
		//pdf->paramOn(frame_m,Layout(0.75));
		//data->plotOn(frame_m,Name("data"),MarkerSize(0.5),Binning(50));
		//gPad->SetLogy();
		frame_m->Draw();
		c1->Print(BsMassDistribution.c_str());
	}


	//setup Histograms
	TH1D* e_data;
	if(L0) e_data  = new TH1D("L0 efficiency",";p_{t} [GeV/c];#epsilon(L0)",6, 1.5, 7.5);
	if(HLT1) e_data = new TH1D("HLT1 efficiency",";p_{t} [GeV/c];#epsilon(HLT1)",6, 1.5, 7.5);
	if(HLT2) e_data = new TH1D("HLT2 efficiency",";p_{t} [GeV/c];#epsilon(HLT2)",6, 1.5, 7.5);


	double e_data_1525 = ( N_TISTOS_1 / N_TIS_1);
	double e_data_1525_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_1)/N_TISTOS_1),2) + TMath::Power((TMath::Sqrt(N_TIS_1)/N_TIS_1),2) ) * e_data_1525;

	double e_data_2535 = ( N_TISTOS_2 / N_TIS_2 );
	double e_data_2535_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_2)/N_TISTOS_2),2) + TMath::Power((TMath::Sqrt(N_TIS_2)/N_TIS_2),2) ) * e_data_2535;

	double e_data_3545 = ( N_TISTOS_3 / N_TIS_3 );
	double e_data_3545_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_3)/N_TISTOS_3),2) + TMath::Power((TMath::Sqrt(N_TIS_3)/N_TIS_3),2) ) * e_data_3545;

	double e_data_4555 = ( N_TISTOS_4 / N_TIS_4 );
	double e_data_4555_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_4)/N_TISTOS_4),2) + TMath::Power((TMath::Sqrt(N_TIS_4)/N_TIS_4),2) ) * e_data_4555;

	double e_data_5565 = ( N_TISTOS_5 / N_TIS_5 );
	double e_data_5565_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_5)/N_TISTOS_5),2) + TMath::Power((TMath::Sqrt(N_TIS_5)/N_TIS_5),2) ) * e_data_5565;

	double e_data_6575 = ( N_TISTOS_6 / N_TIS_6 );
	double e_data_6575_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_TISTOS_6)/N_TISTOS_6),2) + TMath::Power((TMath::Sqrt(N_TIS_6)/N_TIS_6),2) ) * e_data_6575;

	//fill the histos
	if(Systematics){
		e_data->SetBinContent(1,e_data_1525); e_data->SetBinError(1,e_data_1525_error);
		e_data->SetBinContent(2,e_data_2535); e_data->SetBinError(2,e_data_2535_error);
		e_data->SetBinContent(3,e_data_3545); e_data->SetBinError(3,e_data_3545_error);
		e_data->SetBinContent(4,e_data_4555); e_data->SetBinError(4,e_data_4555_error);
		e_data->SetBinContent(5,e_data_5565); e_data->SetBinError(5,e_data_5565_error);
		e_data->SetBinContent(6,e_data_6575); e_data->SetBinError(6,e_data_6575_error);
	}

	//draw histos
	if(L0 && Systematics){
		e_data->SetStats(0);
		e_data->SetAxisRange(0.,1.2,"Y");
		e_data->Draw("E1"); c1->Print("eps/efficiency_L0_data.eps");
	}
	if(HLT1 && Systematics){
		e_data->SetStats(0);
		e_data->SetAxisRange(0.,1.2,"Y");
		e_data->Draw("E1"); c1->Print("eps/efficiency_HLT1_data.eps");
	}



	TFile* output_TISTOS;
	if(L0 && Systematics){
		output_TISTOS = new TFile("rootHistos/L0_data_Hist.root","RECREATE");
		e_data->Draw();
		e_data->Write();
		output_TISTOS->Write();
		output_TISTOS->Close();
	}
	if(HLT1 && Systematics){
		output_TISTOS = new TFile("rootHistos/HLT1_data_Hist.root","RECREATE");
		e_data->Draw();
		e_data->Write();
		output_TISTOS->Write();
		output_TISTOS->Close();
	}

	if(Systematics) output_TISTOS->Close();


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
		 if(sevenTeV && Ds2KKpi) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_afterPreSel_sweight.root","RECREATE");
		 if((!sevenTeV) && Ds2KKpi) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_afterPreSel_sweight.root","RECREATE");
		 if((!sevenTeV) && (!Ds2KKpi)) output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_Ds2pipipi_12_afterPreSel_sweight.root","RECREATE");

		 tree->SetBranchStatus("*",1);
   		 TTree* new_tree = tree->CopyTree("DTF_Bs_M > 5300 && DTF_Bs_M < 5450");
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
	bool Ds2KKpi=true;

///Load file
	TFile* file;
	if(sevenTeV && Ds2KKpi) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_forBDT_DsK3fb_Selection.root");
	if(!sevenTeV && Ds2KKpi) file=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT_DsK3fb_Selection.root");
	if(sevenTeV && (!Ds2KKpi)) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2pipipi_forBDT.root");
	if(!sevenTeV && (!Ds2KKpi)) file=new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_forBDT.root");
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
	int expectedYield = 0.05238 * 0.863 * quickFit(sevenTeV,Ds2KKpi);

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

void taggingCategories(){

TFile* file= new TFile("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");
TTree* tree = (TTree*) file->Get("DecayTree");

int N = tree->GetEntries();
cout << "full MC file contains " << N << " events" <<  endl;

Int_t Bs_TRUEID;
Int_t Bs_TAGDECISION;

Int_t K_plus_TRUEID;

tree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
tree -> SetBranchAddress( "Bs_TAGDECISION" , &Bs_TAGDECISION );

tree -> SetBranchAddress( "K_plus_TRUEID" , &K_plus_TRUEID );

int numEvents = tree->GetEntries();


TFile* output_TagBs_hplus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TagBs_hplus.root","RECREATE");
TTree* TagBs_hplus_tree = tree->CloneTree(0);


for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION > 0){
		if(K_plus_TRUEID > 0) TagBs_hplus_tree->Fill();
	}
}
TagBs_hplus_tree->Write();
output_TagBs_hplus->Close();


TFile* output_TagBsbar_hplus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TagBsbar_hplus.root","RECREATE");
TTree* TagBsbar_hplus_tree = tree->CloneTree(0);


for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION < 0){
		if(K_plus_TRUEID > 0) TagBsbar_hplus_tree->Fill();
	}
}
TagBsbar_hplus_tree->Write();
output_TagBsbar_hplus->Close();


TFile* output_noTag_hplus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_noTag_hplus.root","RECREATE");
TTree* noTag_hplus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION == 0){
		if(K_plus_TRUEID > 0) noTag_hplus_tree->Fill();	
	}

}
noTag_hplus_tree->Write();
output_noTag_hplus->Close();


TFile* output_TagBs_hminus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TagBs_hminus.root","RECREATE");
TTree* TagBs_hminus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION > 0){
		if(K_plus_TRUEID < 0) TagBs_hminus_tree->Fill();
	}
}
TagBs_hminus_tree->Write();
output_TagBs_hminus->Close();


TFile* output_TagBsbar_hminus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TagBsbar_hminus.root","RECREATE");
TTree* TagBsbar_hminus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION < 0){
		if(K_plus_TRUEID < 0) TagBsbar_hminus_tree->Fill();
	}
}
TagBsbar_hminus_tree->Write();
output_TagBsbar_hminus->Close();



TFile* output_noTag_hminus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_noTag_hminus.root","RECREATE");
TTree* noTag_hminus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TAGDECISION == 0){
		if(K_plus_TRUEID < 0) noTag_hminus_tree->Fill();	
	}
}
noTag_hminus_tree->Write();
output_noTag_hminus->Close();


TFile* output_TrueBs_hplus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TrueBs_hplus.root","RECREATE");
TTree* TrueBs_hplus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TRUEID == 531){
		if(K_plus_TRUEID > 0) TrueBs_hplus_tree->Fill();
	}
}
TrueBs_hplus_tree->Write();
output_TrueBs_hplus->Close();

TFile* output_TrueBsbar_hplus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TrueBsbar_hplus.root","RECREATE");
TTree* TrueBsbar_hplus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TRUEID == -531){
		if(K_plus_TRUEID > 0) TrueBsbar_hplus_tree->Fill();
	}
}
TrueBsbar_hplus_tree->Write();
output_TrueBsbar_hplus->Close();



TFile* output_TrueBsbar_hminus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TrueBsbar_hminus.root","RECREATE");
TTree* TrueBsbar_hminus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TRUEID == -531){
		if(K_plus_TRUEID < 0) TrueBsbar_hminus_tree->Fill();
	}
}
TrueBsbar_hminus_tree->Write();
output_TrueBsbar_hminus->Close();



TFile* output_TrueBs_hminus = new TFile("/auto/data/kecke/B2DKPiPi/TaggingCategories_MC/mc_Bs2DsKpipi_Run1_reweighted_TrueBs_hminus.root","RECREATE");
TTree* TrueBs_hminus_tree = tree->CloneTree(0);

for(int i=0; i< numEvents; i++)
{
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

	if(Bs_TRUEID == 531){
		if(K_plus_TRUEID < 0) TrueBs_hminus_tree->Fill();
	}
}


TrueBs_hminus_tree->Write();
output_TrueBs_hminus->Close();

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

void Systematics(bool sevenTeV, bool SignalMode){

	bool useWeights = false;
	bool reweight = true;
	if(useWeights) reweight = false;

	/// data from signal & normalization mode
	TFile* file;
	if(sevenTeV && SignalMode) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs_11_final_sweight.root");
	if((!sevenTeV) && SignalMode) file= new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root");
	if(sevenTeV && (!SignalMode)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data_Bs2Dspipipi_11_final_sweight.root");
	if((!sevenTeV) && (!SignalMode)) file= new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root");
	TTree* tree = (TTree*) file->Get("DecayTree");

	/// mc from signal & normalization mode
	TFile* fileMC;
	if(sevenTeV && SignalMode && (!useWeights)) fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_Reco14.root");
	if((!sevenTeV) && SignalMode && (!useWeights)) fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_Reco14.root");
	if(sevenTeV && (!SignalMode)&& (!useWeights)) fileMC= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
	if((!sevenTeV) && (!SignalMode) && (!useWeights)) fileMC= new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
	if((!SignalMode) && useWeights && sevenTeV) fileMC = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_XdReweighted_Reco14.root");
	if(SignalMode && useWeights && sevenTeV) fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_XsReweighted_Reco14.root");
	if((!SignalMode) && useWeights && (!sevenTeV)) fileMC = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_XdReweighted_Reco14.root");
	if(SignalMode && useWeights && (!sevenTeV)) fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_XsReweighted_Reco14.root");
	TTree* treeMC = (TTree*) fileMC->Get("DecayTree");

	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("N_Bs_sw",1);
	tree->SetBranchStatus( "*PX", 1 );
	tree->SetBranchStatus( "*PY", 1);
	tree->SetBranchStatus( "*PZ", 1);

	treeMC->SetBranchStatus("*",0);
	treeMC->SetBranchStatus( "*PX", 1 );
	treeMC->SetBranchStatus( "*PY", 1);
	treeMC->SetBranchStatus( "*PZ", 1);
	treeMC->SetBranchStatus( "*weight", 1);

	double massKaon = 493.68;
	double massPion = 139.57;
	double w;
	if(!useWeights) w = 1.;
	else if(useWeights) treeMC -> SetBranchAddress( "weight" , &w );

	//signal data
	Double_t K_plus_PX;
	Double_t K_plus_PY;
	Double_t K_plus_PZ;
	Double_t pi_plus_PX;
	Double_t pi_plus_PY;
	Double_t pi_plus_PZ;
	Double_t pi_minus_PX;
	Double_t pi_minus_PY;
	Double_t pi_minus_PZ;
	Double_t N_Bs_sw;

	//signal mc
	Double_t K_plus_PX_MC;
	Double_t K_plus_PY_MC;
	Double_t K_plus_PZ_MC;
	Double_t pi_plus_PX_MC;
	Double_t pi_plus_PY_MC;
	Double_t pi_plus_PZ_MC;
	Double_t pi_minus_PX_MC;
	Double_t pi_minus_PY_MC;
	Double_t pi_minus_PZ_MC;

	//normalization data
	Double_t pi_plus1_PX;
	Double_t pi_plus1_PY;
	Double_t pi_plus1_PZ;
	Double_t pi_plus2_PX;
	Double_t pi_plus2_PY;
	Double_t pi_plus2_PZ;
	Double_t pi_minusN_PX;
	Double_t pi_minusN_PY;
	Double_t pi_minusN_PZ;

	//normalization mc
	Double_t pi_plus1_PX_MC;
	Double_t pi_plus1_PY_MC;
	Double_t pi_plus1_PZ_MC;
	Double_t pi_plus2_PX_MC;
	Double_t pi_plus2_PY_MC;
	Double_t pi_plus2_PZ_MC;
	Double_t pi_minusN_PX_MC;
	Double_t pi_minusN_PY_MC;
	Double_t pi_minusN_PZ_MC;

	TLorentzVector K_plus;
	TLorentzVector pi_minus;
	TLorentzVector pi_plus;

	TLorentzVector K_plus_MC;
	TLorentzVector pi_minus_MC;
	TLorentzVector pi_plus_MC;

	TLorentzVector pi_plus1;
	TLorentzVector pi_minusN;
	TLorentzVector pi_plus2;

	TLorentzVector pi_plus1_MC;
	TLorentzVector pi_minusN_MC;
	TLorentzVector pi_plus2_MC;

	if(SignalMode){

        	tree -> SetBranchAddress( "K_plus_PX" , &K_plus_PX );
        	tree -> SetBranchAddress( "K_plus_PY" , &K_plus_PY );
        	tree -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );
        	tree -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
        	tree -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
        	tree -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );
        	tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
        	tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
        	tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
        	tree -> SetBranchAddress( "N_Bs_sw" , &N_Bs_sw );

        	treeMC -> SetBranchAddress( "K_plus_PX" , &K_plus_PX_MC );
        	treeMC -> SetBranchAddress( "K_plus_PY" , &K_plus_PY_MC );
        	treeMC -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ_MC );
        	treeMC -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX_MC );
        	treeMC -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY_MC );
        	treeMC -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ_MC );
	}

	if(!SignalMode){

        	tree -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
        	tree -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
        	tree -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );
        	tree -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
        	tree -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
        	tree -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );
        	tree -> SetBranchAddress( "pi_minus_PX" , &pi_minusN_PX );
        	tree -> SetBranchAddress( "pi_minus_PY" , &pi_minusN_PY );
        	tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minusN_PZ );
        	tree -> SetBranchAddress( "N_Bs_sw" , &N_Bs_sw );

        	treeMC -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX_MC );
        	treeMC -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY_MC );
        	treeMC -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ_MC );
        	treeMC -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX_MC );
        	treeMC -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY_MC );
        	treeMC -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PX" , &pi_minusN_PX_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PY" , &pi_minusN_PY_MC );
        	treeMC -> SetBranchAddress( "pi_minus_PZ" , &pi_minusN_PZ_MC );
	}

	string TitleX;
	if(SignalMode) TitleX = "X_{s} mass";
	else if (!SignalMode) TitleX = "X_{d} mass";
	string Branch;
	if(SignalMode) Branch = "Xs_Mass";
	if(!SignalMode) Branch = "Xd_Mass";
	int bins = 25;
	double min;
	double max;
	if(SignalMode){
		min = 700.;
		max = 3400.;
	}
	if(!SignalMode){
		min = 500.;
		max = 3100.;
	}

	///Make histograms
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	TString title= ";"+TitleX+";Yield [norm.]";
	TH1D* h= new TH1D(Branch.c_str(),title,bins,min,max);
	TH1D* h_MC= new TH1D((Branch+"_MC").c_str(),title,bins,min,max);
	TH1D* h_MC_rw= new TH1D((Branch+"_MC_rw").c_str(),title,bins,min,max);
	h->SetStats(0);
	h_MC->SetStats(0);
	h_MC_rw->SetStats(0);



	///loop over data events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
	{
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		tree->GetEntry(i);
		if(SignalMode){
			K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
			pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
			pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);
			h->Fill((K_plus + pi_minus + pi_plus).M(),N_Bs_sw);
		}
		if(!SignalMode){
			pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
			pi_minusN.SetXYZM(pi_minusN_PX,pi_minusN_PY,pi_minusN_PZ,massPion);
			pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
			h->Fill((pi_plus1 + pi_minusN + pi_plus2).M(),N_Bs_sw);
		}
	}

	///loop over MC events
	int numEventsMC = treeMC->GetEntries();
	for(int i=0; i< numEventsMC; i++)
	{
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
		treeMC->GetEntry(i);
		if(SignalMode){
			K_plus_MC.SetXYZM(K_plus_PX_MC,K_plus_PY_MC,K_plus_PZ_MC,massKaon);
			pi_minus_MC.SetXYZM(pi_minus_PX_MC,pi_minus_PY_MC,pi_minus_PZ_MC,massPion);
			pi_plus_MC.SetXYZM(pi_plus_PX_MC,pi_plus_PY_MC,pi_plus_PZ_MC,massPion);
			h_MC->Fill((K_plus_MC + pi_minus_MC + pi_plus_MC).M());
			h_MC_rw->Fill((K_plus_MC + pi_minus_MC + pi_plus_MC).M(),w);
		}
		if(!SignalMode){
			pi_plus1_MC.SetXYZM(pi_plus1_PX_MC,pi_plus1_PY_MC,pi_plus1_PZ_MC,massPion);
			pi_minusN_MC.SetXYZM(pi_minusN_PX_MC,pi_minusN_PY_MC,pi_minusN_PZ_MC,massPion);
			pi_plus2_MC.SetXYZM(pi_plus2_PX_MC,pi_plus2_PY_MC,pi_plus2_PZ_MC,massPion);
			h_MC->Fill((pi_plus2_MC + pi_minusN_MC + pi_plus1_MC).M());
			h_MC_rw->Fill((pi_plus2_MC + pi_minusN_MC + pi_plus1_MC).M(),w);
		}
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
		TFile* output=new TFile((Branch+"_weights.root").c_str(),"RECREATE");
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

void reweight(bool SignalMode, bool Trigger){

	///Get weight histos
	TFile* massFile;
	if(!SignalMode && (!Trigger)) massFile = new TFile("Xd_Mass_weights.root");
	if(SignalMode && (!Trigger)) massFile = new TFile("Xs_Mass_weights.root");	
	TH1D* h_mass;
	if(!SignalMode && (!Trigger)) h_mass = (TH1D*)massFile->Get("Xd_Mass_weight");
	if(SignalMode && (!Trigger)) h_mass = (TH1D*)massFile->Get("Xs_Mass_weight");

	TFile* L0_File;
	TFile* HLT1_File;
	TH1D* h_L0;
	TH1D* h_HLT1;
	if(Trigger){
		L0_File = new TFile("rootHistos/L0_weights_Hist.root");
		HLT1_File = new TFile("rootHistos/HLT1_weights_Hist.root");
		h_L0 = (TH1D*)L0_File->Get("L0_weights");
		h_HLT1 = (TH1D*)HLT1_File->Get("HLT1_weights");
	}

	///Load MC file
	TFile* fileMC;
	if(!SignalMode) fileMC= new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root");
	if(SignalMode) fileMC= new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_forBDT_Reco14.root");
   	TTree* treeMC = (TTree*) fileMC->Get("DecayTree");

	Double_t K_plus_PX;
	Double_t K_plus_PY;
	Double_t K_plus_PZ;
	Double_t K_plus_PT;

	Double_t pi_plus_PX;
	Double_t pi_plus_PY;
	Double_t pi_plus_PZ;
	Double_t pi_plus_PT;

	Double_t pi_minus_PX;
	Double_t pi_minus_PY;
	Double_t pi_minus_PZ;
	Double_t pi_minus_PT;

	Double_t pi_plus1_PX;
	Double_t pi_plus1_PY;
	Double_t pi_plus1_PZ;
	Double_t pi_plus1_PT;

	Double_t pi_plus2_PX;
	Double_t pi_plus2_PY;
	Double_t pi_plus2_PZ;
	Double_t pi_plus2_PT;

	Double_t K_plus_fromDs_PT;
	Double_t K_minus_fromDs_PT;
	Double_t pi_minus_fromDs_PT;


	treeMC -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
	treeMC -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
	treeMC -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );

	if(SignalMode){
        	treeMC -> SetBranchAddress( "K_plus_PX" , &K_plus_PX );
        	treeMC -> SetBranchAddress( "K_plus_PY" , &K_plus_PY );
        	treeMC -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );
        	treeMC -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
        	treeMC -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
        	treeMC -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
        	treeMC -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );
        	treeMC -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
        	treeMC -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
        	treeMC -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
        	treeMC -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
        	treeMC -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
	}
	if(!SignalMode){
        	treeMC -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
        	treeMC -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
        	treeMC -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );
        	treeMC -> SetBranchAddress( "pi_plus2_PT" , &pi_plus2_PT );
        	treeMC -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
        	treeMC -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
        	treeMC -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );
        	treeMC -> SetBranchAddress( "pi_plus1_PT" , &pi_plus1_PT );
        	treeMC -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
        	treeMC -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
        	treeMC -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
        	treeMC -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
	}

	TLorentzVector K_plus;
	TLorentzVector pi_minus;
	TLorentzVector pi_plus;
	TLorentzVector pi_plus1;
	TLorentzVector pi_plus2;

	double massKaon = 493.68;
	double massPion = 139.57;
	double massX = 0;
	double maxPt = 0;

	///Create new tree
	double w;
	double w_L0;
	double w_HLT1;
    	///TFile* output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_Reco14.root","RECREATE");
	TFile* output;
	if(!SignalMode && (!Trigger)) output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_XdReweighted_Reco14.root","RECREATE");
	if(SignalMode && (!Trigger)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_XsReweighted_Reco14.root","RECREATE");

	if(!SignalMode && Trigger) output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_withTriggerWeights_Reco14.root","RECREATE");
	if(SignalMode && Trigger) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_withTriggerWeights_Reco14.root","RECREATE");

    	TTree* new_tree = treeMC->CloneTree();
    	TBranch* Bra_w;
	TBranch* Bra_w_L0;
	TBranch* Bra_w_HLT1;
	if(!Trigger) Bra_w = new_tree->Branch("weight",&w,"weight/D");
	if(Trigger){
		Bra_w_L0 = new_tree->Branch("weight_L0",&w_L0,"weight_L0/D");
		Bra_w_HLT1 = new_tree->Branch("weight_HLT1",&w_HLT1,"weight_HLT1/D");
	}

	treeMC->SetBranchStatus("*",0);
	treeMC->SetBranchStatus( "*PX", 1 );
	treeMC->SetBranchStatus( "*PY", 1);
	treeMC->SetBranchStatus( "*PZ", 1);
	treeMC->SetBranchStatus( "*PT", 1);

   	///loop over MC events
   	int numEventsMC = treeMC->GetEntries();
   	for(int i=0; i< numEventsMC; i++)
  	{	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
		treeMC->GetEntry(i);
		w=1.;
		w_L0=1.;
		w_HLT1=1.;

		double tmp=w;
		double tmp_L0=w_L0;
		double tmp_HLT1=w_HLT1;


		if(!Trigger){
			if(SignalMode){
				K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
				pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
				pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);
				massX = (pi_plus + pi_minus + K_plus).M();
			}
			if(!SignalMode){
				pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
				pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
				pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
				massX = (pi_plus1 + pi_plus2 + pi_minus).M();
			}
		
                	tmp=h_mass->GetBinContent(h_mass->FindBin(massX));
			if(tmp==0)tmp=1.;
			w*=tmp;

			Bra_w->Fill();
		}
		if(Trigger){
			if(!SignalMode){
				//maxPt in GeV units
				maxPt = TMath::Max(pi_minus_PT,TMath::Max(pi_plus2_PT,TMath::Max(pi_plus1_PT,TMath::Max(pi_minus_fromDs_PT,TMath::Max(K_plus_fromDs_PT,K_minus_fromDs_PT))))) / 1000;
			}
			if(SignalMode){
				//maxPt in GeV units
				maxPt = TMath::Max(pi_minus_PT,TMath::Max(K_plus_PT,TMath::Max(pi_plus_PT,TMath::Max(pi_minus_fromDs_PT,TMath::Max(K_plus_fromDs_PT,K_minus_fromDs_PT))))) / 1000;
			}

                	tmp_L0=h_L0->GetBinContent(h_L0->FindBin(maxPt));
			if(tmp_L0==0)tmp_L0=1.;
			w_L0*=tmp_L0;

                	tmp_HLT1=h_HLT1->GetBinContent(h_HLT1->FindBin(maxPt));
			if(tmp_HLT1==0)tmp_HLT1=1.;
			w_HLT1*=tmp_HLT1;

			Bra_w_L0->Fill();
			Bra_w_HLT1->Fill();
		}
    	}
	new_tree->Write();
	if(!Trigger) massFile->Close();
	fileMC->Close();
	output->Close();
	if(Trigger){
		L0_File->Close();
		HLT1_File->Close();
	}
}

void TriggerSyst(bool reweight = false){

///baseline for estimation , max pt of tracks between 1.5 - 7.5 GeV , using 6 bins
/*
TFile* fileMC;
fileMC= new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
TTree* treeMC = (TTree*) fileMC->Get("DecayTree");

int N = treeMC->GetEntries();
cout << "full MC file contains " << N << " events" <<  endl;

Double_t K_plus_fromDs_PT;
Double_t K_minus_fromDs_PT;
Double_t pi_minus_fromDs_PT;
Double_t pi_plus1_PT;
Double_t pi_plus2_PT;
Double_t pi_minus_PT;

treeMC -> SetBranchAddress( "pi_plus1_PT" , &pi_plus1_PT );
treeMC -> SetBranchAddress( "pi_plus2_PT" , &pi_plus2_PT );
treeMC -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );

treeMC -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
treeMC -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
treeMC -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );


double maxPt = 0;


TH1D* max_track_pt = new TH1D("max track p_{t}",";p_{t} [MeV];Entries",6, 1500., 7500.);

	///loop over data events
int numEvents = treeMC->GetEntries();
for(int i=0; i< numEvents; i++)
{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	treeMC->GetEntry(i);

	maxPt = TMath::Max(pi_minus_PT,TMath::Max(pi_plus2_PT,TMath::Max(pi_plus1_PT,TMath::Max(pi_minus_fromDs_PT,TMath::Max(K_plus_fromDs_PT,K_minus_fromDs_PT)))));
	max_track_pt->Fill(maxPt);

	}

TCanvas* c= new TCanvas();
max_track_pt->Draw("E1"); c->Print("eps/max_track_pt_Norm_MC.eps");
*/

///load and compute all MC TISTOS efficiencies first

/// load L0 efficiencies
//TISTOS files
TFile* fileMC_L0_1_TISTOS;
fileMC_L0_1_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_alt_pt_1525.root");
TTree* treeMC_L0_1_TISTOS = (TTree*) fileMC_L0_1_TISTOS->Get("DecayTree");
TTree* treeMC_L0_1_TISTOS_inRange = treeMC_L0_1_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_2_TISTOS;
fileMC_L0_2_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_pt_2535.root");
TTree* treeMC_L0_2_TISTOS = (TTree*) fileMC_L0_2_TISTOS->Get("DecayTree");
TTree* treeMC_L0_2_TISTOS_inRange = treeMC_L0_2_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_3_TISTOS;
fileMC_L0_3_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_pt_3545.root");
TTree* treeMC_L0_3_TISTOS = (TTree*) fileMC_L0_3_TISTOS->Get("DecayTree");
TTree* treeMC_L0_3_TISTOS_inRange = treeMC_L0_3_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_4_TISTOS;
fileMC_L0_4_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_pt_4555.root");
TTree* treeMC_L0_4_TISTOS = (TTree*) fileMC_L0_4_TISTOS->Get("DecayTree");
TTree* treeMC_L0_4_TISTOS_inRange = treeMC_L0_4_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_5_TISTOS;
fileMC_L0_5_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_pt_5565.root");
TTree* treeMC_L0_5_TISTOS = (TTree*) fileMC_L0_5_TISTOS->Get("DecayTree");
TTree* treeMC_L0_5_TISTOS_inRange = treeMC_L0_5_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_6_TISTOS;
fileMC_L0_6_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TISTOS_pt_6575.root");
TTree* treeMC_L0_6_TISTOS = (TTree*) fileMC_L0_6_TISTOS->Get("DecayTree");
TTree* treeMC_L0_6_TISTOS_inRange = treeMC_L0_6_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

//TIS only files
TFile* fileMC_L0_1_TIS;
fileMC_L0_1_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_alt_pt_1525.root");
TTree* treeMC_L0_1_TIS = (TTree*) fileMC_L0_1_TIS->Get("DecayTree");
TTree* treeMC_L0_1_TIS_inRange = treeMC_L0_1_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_2_TIS;
fileMC_L0_2_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_pt_2535.root");
TTree* treeMC_L0_2_TIS = (TTree*) fileMC_L0_2_TIS->Get("DecayTree");
TTree* treeMC_L0_2_TIS_inRange = treeMC_L0_2_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_3_TIS;
fileMC_L0_3_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_pt_3545.root");
TTree* treeMC_L0_3_TIS = (TTree*) fileMC_L0_3_TIS->Get("DecayTree");
TTree* treeMC_L0_3_TIS_inRange = treeMC_L0_3_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_4_TIS;
fileMC_L0_4_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_pt_4555.root");
TTree* treeMC_L0_4_TIS = (TTree*) fileMC_L0_4_TIS->Get("DecayTree");
TTree* treeMC_L0_4_TIS_inRange = treeMC_L0_4_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_5_TIS;
fileMC_L0_5_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_pt_5565.root");
TTree* treeMC_L0_5_TIS = (TTree*) fileMC_L0_5_TIS->Get("DecayTree");
TTree* treeMC_L0_5_TIS_inRange = treeMC_L0_5_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_L0_6_TIS;
fileMC_L0_6_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_L0_TIS_pt_6575.root");
TTree* treeMC_L0_6_TIS = (TTree*) fileMC_L0_6_TIS->Get("DecayTree");
TTree* treeMC_L0_6_TIS_inRange = treeMC_L0_6_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

/// load HLT1 efficiencies
//TISTOS files
TFile* fileMC_HLT1_1_TISTOS;
fileMC_HLT1_1_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_1525.root");
TTree* treeMC_HLT1_1_TISTOS = (TTree*) fileMC_HLT1_1_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_1_TISTOS_inRange = treeMC_HLT1_1_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_2_TISTOS;
fileMC_HLT1_2_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_2535.root");
TTree* treeMC_HLT1_2_TISTOS = (TTree*) fileMC_HLT1_2_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_2_TISTOS_inRange = treeMC_HLT1_2_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_3_TISTOS;
fileMC_HLT1_3_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_3545.root");
TTree* treeMC_HLT1_3_TISTOS = (TTree*) fileMC_HLT1_3_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_3_TISTOS_inRange = treeMC_HLT1_3_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_4_TISTOS;
fileMC_HLT1_4_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_4555.root");
TTree* treeMC_HLT1_4_TISTOS = (TTree*) fileMC_HLT1_4_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_4_TISTOS_inRange = treeMC_HLT1_4_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_5_TISTOS;
fileMC_HLT1_5_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_5565.root");
TTree* treeMC_HLT1_5_TISTOS = (TTree*) fileMC_HLT1_5_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_5_TISTOS_inRange = treeMC_HLT1_5_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_6_TISTOS;
fileMC_HLT1_6_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TISTOS_pt_6575.root");
TTree* treeMC_HLT1_6_TISTOS = (TTree*) fileMC_HLT1_6_TISTOS->Get("DecayTree");
TTree* treeMC_HLT1_6_TISTOS_inRange = treeMC_HLT1_6_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");


//TIS only files
TFile* fileMC_HLT1_1_TIS;
fileMC_HLT1_1_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_1525.root");
TTree* treeMC_HLT1_1_TIS = (TTree*) fileMC_HLT1_1_TIS->Get("DecayTree");
TTree* treeMC_HLT1_1_TIS_inRange = treeMC_HLT1_1_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_2_TIS;
fileMC_HLT1_2_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_2535.root");
TTree* treeMC_HLT1_2_TIS = (TTree*) fileMC_HLT1_2_TIS->Get("DecayTree");
TTree* treeMC_HLT1_2_TIS_inRange = treeMC_HLT1_2_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_3_TIS;
fileMC_HLT1_3_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_3545.root");
TTree* treeMC_HLT1_3_TIS = (TTree*) fileMC_HLT1_3_TIS->Get("DecayTree");
TTree* treeMC_HLT1_3_TIS_inRange = treeMC_HLT1_3_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_4_TIS;
fileMC_HLT1_4_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_4555.root");
TTree* treeMC_HLT1_4_TIS = (TTree*) fileMC_HLT1_4_TIS->Get("DecayTree");
TTree* treeMC_HLT1_4_TIS_inRange = treeMC_HLT1_4_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_5_TIS;
fileMC_HLT1_5_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_5565.root");
TTree* treeMC_HLT1_5_TIS = (TTree*) fileMC_HLT1_5_TIS->Get("DecayTree");
TTree* treeMC_HLT1_5_TIS_inRange = treeMC_HLT1_5_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT1_6_TIS;
fileMC_HLT1_6_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT1_TIS_pt_6575.root");
TTree* treeMC_HLT1_6_TIS = (TTree*) fileMC_HLT1_6_TIS->Get("DecayTree");
TTree* treeMC_HLT1_6_TIS_inRange = treeMC_HLT1_6_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

/// load HLT2 efficiencies
//TISTOS files
TFile* fileMC_HLT2_1_TISTOS;
fileMC_HLT2_1_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_1525.root");
TTree* treeMC_HLT2_1_TISTOS = (TTree*) fileMC_HLT2_1_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_1_TISTOS_inRange = treeMC_HLT2_1_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_2_TISTOS;
fileMC_HLT2_2_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_2535.root");
TTree* treeMC_HLT2_2_TISTOS = (TTree*) fileMC_HLT2_2_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_2_TISTOS_inRange = treeMC_HLT2_2_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_3_TISTOS;
fileMC_HLT2_3_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_3545.root");
TTree* treeMC_HLT2_3_TISTOS = (TTree*) fileMC_HLT2_3_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_3_TISTOS_inRange = treeMC_HLT2_3_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_4_TISTOS;
fileMC_HLT2_4_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_4555.root");
TTree* treeMC_HLT2_4_TISTOS = (TTree*) fileMC_HLT2_4_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_4_TISTOS_inRange = treeMC_HLT2_4_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_5_TISTOS;
fileMC_HLT2_5_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_5565.root");
TTree* treeMC_HLT2_5_TISTOS = (TTree*) fileMC_HLT2_5_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_5_TISTOS_inRange = treeMC_HLT2_5_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_6_TISTOS;
fileMC_HLT2_6_TISTOS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TISTOS_pt_6575.root");
TTree* treeMC_HLT2_6_TISTOS = (TTree*) fileMC_HLT2_6_TISTOS->Get("DecayTree");
TTree* treeMC_HLT2_6_TISTOS_inRange = treeMC_HLT2_6_TISTOS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");


//TIS only files
TFile* fileMC_HLT2_1_TIS;
fileMC_HLT2_1_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_1525.root");
TTree* treeMC_HLT2_1_TIS = (TTree*) fileMC_HLT2_1_TIS->Get("DecayTree");
TTree* treeMC_HLT2_1_TIS_inRange = treeMC_HLT2_1_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_2_TIS;
fileMC_HLT2_2_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_2535.root");
TTree* treeMC_HLT2_2_TIS = (TTree*) fileMC_HLT2_2_TIS->Get("DecayTree");
TTree* treeMC_HLT2_2_TIS_inRange = treeMC_HLT2_2_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_3_TIS;
fileMC_HLT2_3_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_3545.root");
TTree* treeMC_HLT2_3_TIS = (TTree*) fileMC_HLT2_3_TIS->Get("DecayTree");
TTree* treeMC_HLT2_3_TIS_inRange = treeMC_HLT2_3_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_4_TIS;
fileMC_HLT2_4_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_4555.root");
TTree* treeMC_HLT2_4_TIS = (TTree*) fileMC_HLT2_4_TIS->Get("DecayTree");
TTree* treeMC_HLT2_4_TIS_inRange = treeMC_HLT2_4_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_5_TIS;
fileMC_HLT2_5_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_5565.root");
TTree* treeMC_HLT2_5_TIS = (TTree*) fileMC_HLT2_5_TIS->Get("DecayTree");
TTree* treeMC_HLT2_5_TIS_inRange = treeMC_HLT2_5_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");

TFile* fileMC_HLT2_6_TIS;
fileMC_HLT2_6_TIS = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_6575.root");
TTree* treeMC_HLT2_6_TIS = (TTree*) fileMC_HLT2_6_TIS->Get("DecayTree");
TTree* treeMC_HLT2_6_TIS_inRange = treeMC_HLT2_6_TIS->CopyTree("Bs_MM > 5300 && Bs_MM < 5420");


//setup Histograms
TH1D* e_L0_MC = new TH1D("L0 efficiency",";p_{t} [GeV/c];#epsilon(L0)",6, 1.5, 7.5);
TH1D* e_HLT1_MC = new TH1D("HLT1 efficiency",";p_{t} [GeV/c];#epsilon(HLT1)",6, 1.5, 7.5);
TH1D* e_HLT2_MC = new TH1D("HLT2 efficiency",";p_{t} [GeV/c];#epsilon(HLT2)",6, 1.5, 7.5);

//store number of events in doubles
double N_L0_TISTOS_1 = treeMC_L0_1_TISTOS_inRange->GetEntries();
double N_L0_TISTOS_2 = treeMC_L0_2_TISTOS_inRange->GetEntries();
double N_L0_TISTOS_3 = treeMC_L0_3_TISTOS_inRange->GetEntries();
double N_L0_TISTOS_4 = treeMC_L0_4_TISTOS_inRange->GetEntries();
double N_L0_TISTOS_5 = treeMC_L0_5_TISTOS_inRange->GetEntries();
double N_L0_TISTOS_6 = treeMC_L0_6_TISTOS_inRange->GetEntries();

double N_L0_TIS_1 = treeMC_L0_1_TIS_inRange->GetEntries();
double N_L0_TIS_2 = treeMC_L0_2_TIS_inRange->GetEntries();
double N_L0_TIS_3 = treeMC_L0_3_TIS_inRange->GetEntries();
double N_L0_TIS_4 = treeMC_L0_4_TIS_inRange->GetEntries();
double N_L0_TIS_5 = treeMC_L0_5_TIS_inRange->GetEntries();
double N_L0_TIS_6 = treeMC_L0_6_TIS_inRange->GetEntries();

double N_HLT1_TISTOS_1 = treeMC_HLT1_1_TISTOS_inRange->GetEntries();
double N_HLT1_TISTOS_2 = treeMC_HLT1_2_TISTOS_inRange->GetEntries();
double N_HLT1_TISTOS_3 = treeMC_HLT1_3_TISTOS_inRange->GetEntries();
double N_HLT1_TISTOS_4 = treeMC_HLT1_4_TISTOS_inRange->GetEntries();
double N_HLT1_TISTOS_5 = treeMC_HLT1_5_TISTOS_inRange->GetEntries();
double N_HLT1_TISTOS_6 = treeMC_HLT1_6_TISTOS_inRange->GetEntries();

double N_HLT1_TIS_1 = treeMC_HLT1_1_TIS_inRange->GetEntries();
double N_HLT1_TIS_2 = treeMC_HLT1_2_TIS_inRange->GetEntries();
double N_HLT1_TIS_3 = treeMC_HLT1_3_TIS_inRange->GetEntries();
double N_HLT1_TIS_4 = treeMC_HLT1_4_TIS_inRange->GetEntries();
double N_HLT1_TIS_5 = treeMC_HLT1_5_TIS_inRange->GetEntries();
double N_HLT1_TIS_6 = treeMC_HLT1_6_TIS_inRange->GetEntries();

double N_HLT2_TISTOS_1 = treeMC_HLT2_1_TISTOS_inRange->GetEntries();
double N_HLT2_TISTOS_2 = treeMC_HLT2_2_TISTOS_inRange->GetEntries();
double N_HLT2_TISTOS_3 = treeMC_HLT2_3_TISTOS_inRange->GetEntries();
double N_HLT2_TISTOS_4 = treeMC_HLT2_4_TISTOS_inRange->GetEntries();
double N_HLT2_TISTOS_5 = treeMC_HLT2_5_TISTOS_inRange->GetEntries();
double N_HLT2_TISTOS_6 = treeMC_HLT2_6_TISTOS_inRange->GetEntries();

double N_HLT2_TIS_1 = treeMC_HLT2_1_TIS_inRange->GetEntries();
double N_HLT2_TIS_2 = treeMC_HLT2_2_TIS_inRange->GetEntries();
double N_HLT2_TIS_3 = treeMC_HLT2_3_TIS_inRange->GetEntries();
double N_HLT2_TIS_4 = treeMC_HLT2_4_TIS_inRange->GetEntries();
double N_HLT2_TIS_5 = treeMC_HLT2_5_TIS_inRange->GetEntries();
double N_HLT2_TIS_6 = treeMC_HLT2_6_TIS_inRange->GetEntries();

//compute efficiencies for different p_t bins
double e_L0_MC_1525 = ( N_L0_TISTOS_1 / N_L0_TIS_1);
double e_L0_MC_1525_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_1)/N_L0_TISTOS_1),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_1)/N_L0_TIS_1),2) ) * e_L0_MC_1525;

double e_L0_MC_2535 = ( N_L0_TISTOS_2 / N_L0_TIS_2 );
double e_L0_MC_2535_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_2)/N_L0_TISTOS_2),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_2)/N_L0_TIS_2),2) ) * e_L0_MC_2535;

double e_L0_MC_3545 = ( N_L0_TISTOS_3 / N_L0_TIS_3 );
double e_L0_MC_3545_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_3)/N_L0_TISTOS_3),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_3)/N_L0_TIS_3),2) ) * e_L0_MC_3545;

double e_L0_MC_4555 = ( N_L0_TISTOS_4 / N_L0_TIS_4 );
double e_L0_MC_4555_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_4)/N_L0_TISTOS_4),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_4)/N_L0_TIS_4),2) ) * e_L0_MC_4555;

double e_L0_MC_5565 = ( N_L0_TISTOS_5 / N_L0_TIS_5 );
double e_L0_MC_5565_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_5)/N_L0_TISTOS_5),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_5)/N_L0_TIS_5),2) ) * e_L0_MC_5565;

double e_L0_MC_6575 = ( N_L0_TISTOS_6 / N_L0_TIS_6 );
double e_L0_MC_6575_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_L0_TISTOS_6)/N_L0_TISTOS_6),2) + TMath::Power((TMath::Sqrt(N_L0_TIS_6)/N_L0_TIS_6),2) ) * e_L0_MC_6575;



double e_HLT1_MC_1525 = ( N_HLT1_TISTOS_1 / N_HLT1_TIS_1);
double e_HLT1_MC_1525_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_1)/N_HLT1_TISTOS_1),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_1)/N_HLT1_TIS_1),2) ) * e_HLT1_MC_1525;

double e_HLT1_MC_2535 = ( N_HLT1_TISTOS_2 / N_HLT1_TIS_2 );
double e_HLT1_MC_2535_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_2)/N_HLT1_TISTOS_2),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_2)/N_HLT1_TIS_2),2) ) * e_HLT1_MC_2535;

double e_HLT1_MC_3545 = ( N_HLT1_TISTOS_3 / N_HLT1_TIS_3 );
double e_HLT1_MC_3545_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_3)/N_HLT1_TISTOS_3),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_3)/N_HLT1_TIS_3),2) ) * e_HLT1_MC_3545;

double e_HLT1_MC_4555 = ( N_HLT1_TISTOS_4 / N_HLT1_TIS_4 );
double e_HLT1_MC_4555_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_4)/N_HLT1_TISTOS_4),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_4)/N_HLT1_TIS_4),2) ) * e_HLT1_MC_4555;

double e_HLT1_MC_5565 = ( N_HLT1_TISTOS_5 / N_HLT1_TIS_5 );
double e_HLT1_MC_5565_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_5)/N_HLT1_TISTOS_5),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_5)/N_HLT1_TIS_5),2) ) * e_HLT1_MC_5565;

double e_HLT1_MC_6575 = ( N_HLT1_TISTOS_6 / N_HLT1_TIS_6 );
double e_HLT1_MC_6575_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT1_TISTOS_6)/N_HLT1_TISTOS_6),2) + TMath::Power((TMath::Sqrt(N_HLT1_TIS_6)/N_HLT1_TIS_6),2) ) * e_HLT1_MC_6575;



double e_HLT2_MC_1525 = ( N_HLT2_TISTOS_1 / N_HLT2_TIS_1);
double e_HLT2_MC_1525_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_1)/N_HLT2_TISTOS_1),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_1)/N_HLT2_TIS_1),2) ) * e_HLT2_MC_1525;

double e_HLT2_MC_2535 = ( N_HLT2_TISTOS_2 / N_HLT2_TIS_2 );
double e_HLT2_MC_2535_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_2)/N_HLT2_TISTOS_2),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_2)/N_HLT2_TIS_2),2) ) * e_HLT2_MC_2535;

double e_HLT2_MC_3545 = ( N_HLT2_TISTOS_3 / N_HLT2_TIS_3 );
double e_HLT2_MC_3545_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_3)/N_HLT2_TISTOS_3),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_3)/N_HLT2_TIS_3),2) ) * e_HLT2_MC_3545;

double e_HLT2_MC_4555 = ( N_HLT2_TISTOS_4 / N_HLT2_TIS_4 );
double e_HLT2_MC_4555_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_4)/N_HLT2_TISTOS_4),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_4)/N_HLT2_TIS_4),2) ) * e_HLT2_MC_4555;

double e_HLT2_MC_5565 = ( N_HLT2_TISTOS_5 / N_HLT2_TIS_5 );
double e_HLT2_MC_5565_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_5)/N_HLT2_TISTOS_5),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_5)/N_HLT2_TIS_5),2) ) * e_HLT2_MC_5565;

double e_HLT2_MC_6575 = ( N_HLT2_TISTOS_6 / N_HLT2_TIS_6 );
double e_HLT2_MC_6575_error = TMath::Sqrt( TMath::Power((TMath::Sqrt(N_HLT2_TISTOS_6)/N_HLT2_TISTOS_6),2) + TMath::Power((TMath::Sqrt(N_HLT2_TIS_6)/N_HLT2_TIS_6),2) ) * e_HLT2_MC_6575;



//fill the histos
e_L0_MC->SetBinContent(1,e_L0_MC_1525); e_L0_MC->SetBinError(1,e_L0_MC_1525_error);
e_L0_MC->SetBinContent(2,e_L0_MC_2535); e_L0_MC->SetBinError(2,e_L0_MC_2535_error);
e_L0_MC->SetBinContent(3,e_L0_MC_3545); e_L0_MC->SetBinError(3,e_L0_MC_3545_error);
e_L0_MC->SetBinContent(4,e_L0_MC_4555); e_L0_MC->SetBinError(4,e_L0_MC_4555_error);
e_L0_MC->SetBinContent(5,e_L0_MC_5565); e_L0_MC->SetBinError(5,e_L0_MC_5565_error);
e_L0_MC->SetBinContent(6,e_L0_MC_6575); e_L0_MC->SetBinError(6,e_L0_MC_6575_error);

e_HLT1_MC->SetBinContent(1,e_HLT1_MC_1525); e_HLT1_MC->SetBinError(1,e_HLT1_MC_1525_error);
e_HLT1_MC->SetBinContent(2,e_HLT1_MC_2535); e_HLT1_MC->SetBinError(2,e_HLT1_MC_2535_error);
e_HLT1_MC->SetBinContent(3,e_HLT1_MC_3545); e_HLT1_MC->SetBinError(3,e_HLT1_MC_3545_error);
e_HLT1_MC->SetBinContent(4,e_HLT1_MC_4555); e_HLT1_MC->SetBinError(4,e_HLT1_MC_4555_error);
e_HLT1_MC->SetBinContent(5,e_HLT1_MC_5565); e_HLT1_MC->SetBinError(5,e_HLT1_MC_5565_error);
e_HLT1_MC->SetBinContent(6,e_HLT1_MC_6575); e_HLT1_MC->SetBinError(6,e_HLT1_MC_6575_error);

e_HLT2_MC->SetBinContent(1,e_HLT2_MC_1525); e_HLT2_MC->SetBinError(1,e_HLT2_MC_1525_error);
e_HLT2_MC->SetBinContent(2,e_HLT2_MC_2535); e_HLT2_MC->SetBinError(2,e_HLT2_MC_2535_error);
e_HLT2_MC->SetBinContent(3,e_HLT2_MC_3545); e_HLT2_MC->SetBinError(3,e_HLT2_MC_3545_error);
e_HLT2_MC->SetBinContent(4,e_HLT2_MC_4555); e_HLT2_MC->SetBinError(4,e_HLT2_MC_4555_error);
e_HLT2_MC->SetBinContent(5,e_HLT2_MC_5565); e_HLT2_MC->SetBinError(5,e_HLT2_MC_5565_error);
e_HLT2_MC->SetBinContent(6,e_HLT2_MC_6575); e_HLT2_MC->SetBinError(6,e_HLT2_MC_6575_error);

//draw histos
TCanvas* c= new TCanvas();
e_L0_MC->SetStats(0);
e_L0_MC->SetAxisRange(0.,1.2,"Y");
e_L0_MC->Draw("E1"); c->Print("eps/efficiency_L0_MC.eps");
e_HLT1_MC->SetStats(0);
e_HLT1_MC->SetAxisRange(0.,1.2,"Y");
e_HLT1_MC->Draw("E1"); c->Print("eps/efficiency_HLT1_MC.eps");
e_HLT2_MC->SetStats(0);
e_HLT2_MC->SetAxisRange(0.,1.2,"Y");
e_HLT2_MC->Draw("E1"); c->Print("eps/efficiency_HLT2_MC.eps");


///create the MC vs Data Histos
//load L0 data efficiency
TFile* L0_data = new TFile("rootHistos/L0_data_Hist.root");
TH1D* e_L0_data = (TH1D*) L0_data->Get("L0 efficiency");

double maxY= e_L0_data->GetMaximum();
if(e_L0_MC->GetMaximum()>maxY)maxY=e_L0_MC->GetMaximum();
e_L0_data->SetMinimum(0.);
e_L0_data->SetMaximum(maxY*1.2);
e_L0_data->SetLineColor(kBlack);
e_L0_data->Draw("");
e_L0_MC->SetLineColor(kRed);
e_L0_MC->Draw("same");

TLegend *leg = new TLegend(0.55,0.7,0.85,0.85);
leg->SetHeader(" ");
leg->AddEntry(e_L0_data,"Data","LEP");
leg->AddEntry(e_L0_MC,"MC","LEP");
leg->SetLineColor(kWhite);
leg->SetFillColor(kWhite);
leg->SetTextSize(0.05);
leg->Draw();

//calculate weights
if(reweight){
        TFile* output_L0_weights=new TFile("rootHistos/L0_weights_Hist.root","RECREATE");
        TH1D *h_weight_L0 = (TH1D*)e_L0_MC->Clone();
        h_weight_L0->SetName("L0_weights");
        h_weight_L0->Divide(e_L0_data,e_L0_MC);
        h_weight_L0->SetMinimum(0.);
        h_weight_L0->SetMaximum(h_weight_L0->GetMaximum()*1.2);
        h_weight_L0->GetYaxis()->SetTitle("Weight");
        h_weight_L0->Draw();
        c->Print("eps/TriggerSyst/L0_efficiency_weights.eps");
        h_weight_L0->Write();
        output_L0_weights->Write();
        output_L0_weights->Close();
}


c->Print("eps/TriggerSyst/L0_efficiency_comparison.eps");
L0_data->Close();


//load HLT1 data efficiency
TFile* HLT1_data = new TFile("rootHistos/HLT1_data_Hist.root");
TH1D* e_HLT1_data = (TH1D*) HLT1_data->Get("HLT1 efficiency");

double maxY_HLT1= e_HLT1_data->GetMaximum();
if(e_HLT1_MC->GetMaximum()>maxY)maxY_HLT1=e_HLT1_MC->GetMaximum();
e_HLT1_data->SetMinimum(0.);
e_HLT1_data->SetMaximum(maxY_HLT1*1.2);
e_HLT1_data->SetLineColor(kBlack);
e_HLT1_data->Draw("");
e_HLT1_MC->SetLineColor(kRed);
e_HLT1_MC->Draw("same");

TLegend *leg_HLT1 = new TLegend(0.55,0.7,0.85,0.85);
leg_HLT1->SetHeader(" ");
leg_HLT1->AddEntry(e_HLT1_data,"Data","LEP");
leg_HLT1->AddEntry(e_HLT1_MC,"MC","LEP");
leg_HLT1->SetLineColor(kWhite);
leg_HLT1->SetFillColor(kWhite);
leg_HLT1->SetTextSize(0.05);
leg_HLT1->Draw();

//calculate weights
if(reweight){
        TFile* output_HLT1_weights=new TFile("rootHistos/HLT1_weights_Hist.root","RECREATE");
        TH1D *h_weight_HLT1 = (TH1D*)e_HLT1_MC->Clone();
        h_weight_HLT1->SetName("HLT1_weights");
        h_weight_HLT1->Divide(e_HLT1_data,e_HLT1_MC);
        h_weight_HLT1->SetMinimum(0.);
        h_weight_HLT1->SetMaximum(h_weight_HLT1->GetMaximum()*1.2);
        h_weight_HLT1->GetYaxis()->SetTitle("Weight");
        h_weight_HLT1->Draw();
        c->Print("eps/TriggerSyst/HLT1_efficiency_weights.eps");
        h_weight_HLT1->Write();
        output_HLT1_weights->Write();
        output_HLT1_weights->Close();
}

c->Print("eps/TriggerSyst/HLT1_efficiency_comparison.eps");
HLT1_data->Close();


//close all files at the end
fileMC_L0_1_TISTOS->Close();
fileMC_L0_2_TISTOS->Close();
fileMC_L0_3_TISTOS->Close();
fileMC_L0_4_TISTOS->Close();
fileMC_L0_5_TISTOS->Close();
fileMC_L0_6_TISTOS->Close();

fileMC_L0_1_TIS->Close();
fileMC_L0_2_TIS->Close();
fileMC_L0_3_TIS->Close();
fileMC_L0_4_TIS->Close();
fileMC_L0_5_TIS->Close();
fileMC_L0_6_TIS->Close();

fileMC_HLT1_1_TISTOS->Close();
fileMC_HLT1_2_TISTOS->Close();
fileMC_HLT1_3_TISTOS->Close();
fileMC_HLT1_4_TISTOS->Close();
fileMC_HLT1_5_TISTOS->Close();
fileMC_HLT1_6_TISTOS->Close();

fileMC_HLT1_1_TIS->Close();
fileMC_HLT1_2_TIS->Close();
fileMC_HLT1_3_TIS->Close();
fileMC_HLT1_4_TIS->Close();
fileMC_HLT1_5_TIS->Close();
fileMC_HLT1_6_TIS->Close();

fileMC_HLT2_1_TISTOS->Close();
fileMC_HLT2_2_TISTOS->Close();
fileMC_HLT2_3_TISTOS->Close();
fileMC_HLT2_4_TISTOS->Close();
fileMC_HLT2_5_TISTOS->Close();
fileMC_HLT2_6_TISTOS->Close();

fileMC_HLT2_1_TIS->Close();
fileMC_HLT2_2_TIS->Close();
fileMC_HLT2_3_TIS->Close();
fileMC_HLT2_4_TIS->Close();
fileMC_HLT2_5_TIS->Close();
fileMC_HLT2_6_TIS->Close();
}


void TrackingEff(){ 



///compute tracking efficiency ratio from MC

bool sevenTeV = true;

//load tracking eff histo from tracking group
TFile* tracking;
if(sevenTeV) tracking = new TFile("rootHistos/ratio2011S20MC17.root");
if(!sevenTeV) tracking = new TFile("rootHistos/ratio2012S20.root");
TH2D* effHisto = (TH2D*) tracking->Get("Ratio");


TFile* fileMCNorm;
if(sevenTeV) fileMCNorm= new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
if(!sevenTeV) fileMCNorm= new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");
TTree* treeMC_Norm = (TTree*) fileMCNorm->Get("DecayTree");


TFile* fileMCSig;
if(sevenTeV) fileMCSig= new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_Reco14.root");
if(!sevenTeV) fileMCSig= new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_Reco14.root");
TTree* treeMC_Sig = (TTree*) fileMCSig->Get("DecayTree");


Double_t K_plus_fromDs_P_Norm;
Double_t K_minus_fromDs_P_Norm;
Double_t pi_minus_fromDs_P_Norm;
Double_t pi_plus1_P_Norm;
Double_t pi_plus2_P_Norm;
Double_t pi_minus_P_Norm;
Double_t K_plus_fromDs_ETA_Norm;
Double_t K_minus_fromDs_ETA_Norm;
Double_t pi_minus_fromDs_ETA_Norm;
Double_t pi_plus1_ETA_Norm;
Double_t pi_plus2_ETA_Norm;
Double_t pi_minus_ETA_Norm;

Double_t Bs_DTF_M_Norm;
Double_t weight_Norm;

Double_t K_plus_fromDs_P_Sig;
Double_t K_minus_fromDs_P_Sig;
Double_t pi_minus_fromDs_P_Sig;
Double_t K_plus_P_Sig;
Double_t pi_plus_P_Sig;
Double_t pi_minus_P_Sig;
Double_t K_plus_fromDs_ETA_Sig;
Double_t K_minus_fromDs_ETA_Sig;
Double_t pi_minus_fromDs_ETA_Sig;
Double_t K_plus_ETA_Sig;
Double_t pi_plus_ETA_Sig;
Double_t pi_minus_ETA_Sig;

Double_t Bs_DTF_M_Sig;
Double_t weight_Sig;

treeMC_Norm -> SetBranchAddress( "pi_plus1_P" , &pi_plus1_P_Norm );
treeMC_Norm -> SetBranchAddress( "pi_plus2_P" , &pi_plus2_P_Norm );
treeMC_Norm -> SetBranchAddress( "pi_minus_P" , &pi_minus_P_Norm );
treeMC_Norm -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P_Norm );
treeMC_Norm -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P_Norm );
treeMC_Norm -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P_Norm );
treeMC_Norm -> SetBranchAddress( "pi_plus1_ETA" , &pi_plus1_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "pi_plus2_ETA" , &pi_plus2_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "pi_minus_ETA" , &pi_minus_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "K_plus_fromDs_ETA" , &K_plus_fromDs_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "K_minus_fromDs_ETA" , &K_minus_fromDs_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "pi_minus_fromDs_ETA" , &pi_minus_fromDs_ETA_Norm );
treeMC_Norm -> SetBranchAddress( "DTF_Bs_M" , &Bs_DTF_M_Norm );
treeMC_Norm -> SetBranchAddress( "weight" , &weight_Norm );

treeMC_Sig -> SetBranchAddress( "pi_plus_P" , &pi_plus_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_P" , &K_plus_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_P" , &pi_minus_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_plus_ETA" , &pi_plus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "pi_plus_ETA" , &pi_plus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_ETA" , &pi_minus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_fromDs_ETA" , &K_plus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "K_minus_fromDs_ETA" , &K_minus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_fromDs_ETA" , &pi_minus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "DTF_Bs_M" , &Bs_DTF_M_Sig );
treeMC_Sig -> SetBranchAddress( "weight" , &weight_Sig );


//define efficiencies
double eff_pi_plus1_Norm = 0;
double eff_pi_plus2_Norm = 0;
double eff_pi_minus_Norm = 0;
double eff_K_plus_fromDs_Norm = 0;
double eff_K_minus_fromDs_Norm = 0;
double eff_pi_minus_fromDs_Norm = 0;

double eff_pi_plus1_error_Norm = 0;
double eff_pi_plus2_error_Norm = 0;
double eff_pi_minus_error_Norm = 0;
double eff_K_plus_fromDs_error_Norm = 0;
double eff_K_minus_fromDs_error_Norm = 0;
double eff_pi_minus_fromDs_error_Norm = 0;

int normalization_Norm = 0;
double trackingWeight_Norm = 0;
double trackingWeight_error_Norm = 0;

double weighted_efficiency_Norm = 0;
double weighted_efficiency_error_Norm = 0;

double eff_pi_plus_Sig = 0;
double eff_pi_minus_Sig = 0;
double eff_K_plus_Sig = 0;
double eff_K_plus_fromDs_Sig = 0;
double eff_K_minus_fromDs_Sig = 0;
double eff_pi_minus_fromDs_Sig = 0;

double eff_pi_plus_error_Sig = 0;
double eff_pi_minus_error_Sig = 0;
double eff_K_plus_error_Sig = 0;
double eff_K_plus_fromDs_error_Sig = 0;
double eff_K_minus_fromDs_error_Sig = 0;
double eff_pi_minus_fromDs_error_Sig = 0;

int normalization_Sig = 0;
double trackingWeight_Sig = 0;
double trackingWeight_error_Sig = 0;

double weighted_efficiency_Sig = 0;
double weighted_efficiency_error_Sig = 0;

/*
	TFile* output_Norm;
	if(sevenTeV) output_Norm = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc11_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_wTrackingWeights_SigRegion_Reco14.root","RECREATE");
	if(!sevenTeV) output_Norm = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_wTrackingWeights_SigRegion_Reco14.root","RECREATE");

	TFile* output_Sig;
	if(sevenTeV) output_Sig = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BDT_reweighted_wTrackingWeights_SigRegion_Reco14.root","RECREATE");
	if(!sevenTeV) output_Sig = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BDT_reweighted_wTrackingWeights_SigRegion_Reco14.root","RECREATE");
*/

///loop over Norm MC
int numEvents_Norm = treeMC_Norm->GetEntries();
for(int i=0; i< numEvents_Norm; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Norm << endl;
	treeMC_Norm->GetEntry(i);

	//only signal region
	if(Bs_DTF_M_Norm < 5300 || Bs_DTF_M_Norm > 5420) continue;
	

	//get efficiency ratio for every track from Histo
	eff_pi_plus1_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus1_P_Norm/1000),pi_plus1_ETA_Norm));
	eff_pi_plus1_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus1_P_Norm/1000),pi_plus1_ETA_Norm));
	if(eff_pi_plus1_Norm == 0){
		eff_pi_plus1_Norm = 1;
		eff_pi_plus1_error_Norm = 0;
	}

	eff_pi_plus2_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus2_P_Norm/1000),pi_plus2_ETA_Norm));
	eff_pi_plus2_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus2_P_Norm/1000),pi_plus2_ETA_Norm));
	if(eff_pi_plus2_Norm == 0){
		eff_pi_plus2_Norm = 1;
		eff_pi_plus2_error_Norm = 0;
	}

	eff_pi_minus_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Norm/1000),pi_minus_ETA_Norm));
	eff_pi_minus_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Norm/1000),pi_minus_ETA_Norm));
	if(eff_pi_minus_Norm == 0){ 
		eff_pi_minus_Norm = 1;
		eff_pi_minus_error_Norm = 0;
	}

	eff_K_plus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Norm/1000),K_plus_fromDs_ETA_Norm));
	eff_K_plus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Norm/1000),K_plus_fromDs_ETA_Norm));
	if(eff_K_plus_fromDs_Norm == 0){ 
		eff_K_plus_fromDs_Norm =1;
		eff_K_plus_fromDs_error_Norm =0;
	}

	eff_K_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Norm/1000),K_minus_fromDs_ETA_Norm));
	eff_K_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Norm/1000),K_minus_fromDs_ETA_Norm));
	if(eff_K_minus_fromDs_Norm == 0){
		eff_K_minus_fromDs_Norm =1;
		eff_K_minus_fromDs_error_Norm =0;
	}

	eff_pi_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Norm/1000),pi_minus_fromDs_ETA_Norm));
	eff_pi_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Norm/1000),pi_minus_fromDs_ETA_Norm));
	if(eff_pi_minus_fromDs_Norm == 0){ 
		eff_pi_minus_fromDs_Norm =1;
		eff_pi_minus_fromDs_error_Norm =0;
	}

	//compute event tracking efficiency as product of single track efficiencies
	trackingWeight_Norm = eff_pi_plus1_Norm * eff_pi_plus2_Norm * eff_pi_minus_Norm * eff_K_plus_fromDs_Norm * eff_K_minus_fromDs_Norm * eff_pi_minus_fromDs_Norm;
	trackingWeight_error_Norm = TMath::Sqrt(TMath::Power(eff_pi_plus1_error_Norm,2) + TMath::Power(eff_pi_plus2_error_Norm,2) + TMath::Power(eff_pi_minus_error_Norm,2) + TMath::Power(eff_K_plus_fromDs_error_Norm,2) + TMath::Power(eff_K_minus_fromDs_error_Norm,2) + TMath::Power(eff_pi_minus_fromDs_error_Norm,2));

	//compute weighted efficiency ratio, using the MC weights
	weighted_efficiency_Norm = weighted_efficiency_Norm + (weight_Norm * trackingWeight_Norm);
	weighted_efficiency_error_Norm = weighted_efficiency_error_Norm + trackingWeight_error_Norm;
	//conuter for normalization
	normalization_Norm++;
	}



///loop over Sig MC
int numEvents_Sig = treeMC_Sig->GetEntries();
for(int i=0; i< numEvents_Sig; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Sig << endl;
	treeMC_Sig->GetEntry(i);

	//only signal region
	if(Bs_DTF_M_Sig < 5300 || Bs_DTF_M_Sig > 5420) continue;
	

	//get efficiency ratio for every track from Histo
	eff_pi_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_plus_P_Sig/1000),pi_plus_ETA_Sig));
	eff_pi_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_plus_P_Sig/1000),pi_plus_ETA_Sig));
	if(eff_pi_plus_Sig == 0){
		 eff_pi_plus_Sig = 1;
		 eff_pi_plus_error_Sig = 0;
	}

	eff_K_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_P_Sig/1000),K_plus_ETA_Sig));
	eff_K_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_P_Sig/1000),K_plus_ETA_Sig));
	if(eff_K_plus_Sig == 0){
		eff_K_plus_Sig = 1;
		eff_K_plus_error_Sig = 0;
	}

	eff_pi_minus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Sig/1000),pi_minus_ETA_Sig));
	eff_pi_minus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Sig/1000),pi_minus_ETA_Sig));
	if(eff_pi_minus_Sig == 0){
		eff_pi_minus_Sig = 1;
		eff_pi_minus_error_Sig = 0;		
	}

	eff_K_plus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Sig/1000),K_plus_fromDs_ETA_Sig));
	eff_K_plus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Sig/1000),K_plus_fromDs_ETA_Sig));
	if(eff_K_plus_fromDs_Sig == 0){
		eff_K_plus_fromDs_Sig =1;
		eff_K_plus_fromDs_error_Sig =0;
	}

	eff_K_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Sig/1000),K_minus_fromDs_ETA_Sig));
	eff_K_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Sig/1000),K_minus_fromDs_ETA_Sig));
	if(eff_K_minus_fromDs_Sig == 0){
		eff_K_minus_fromDs_Sig =1;
		eff_K_minus_fromDs_error_Sig =0;
	}


	eff_pi_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Sig/1000),pi_minus_fromDs_ETA_Sig));
	eff_pi_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Sig/1000),pi_minus_fromDs_ETA_Sig));
	if(eff_pi_minus_fromDs_Sig == 0){
		eff_pi_minus_fromDs_Sig =1;
		eff_pi_minus_fromDs_error_Sig =0;
	}

	//compute event tracking efficiency as product of single track efficiencies
	trackingWeight_Sig = eff_pi_plus_Sig * eff_K_plus_Sig * eff_pi_minus_Sig * eff_K_plus_fromDs_Sig * eff_K_minus_fromDs_Sig * eff_pi_minus_fromDs_Sig;
	trackingWeight_error_Sig = TMath::Sqrt(TMath::Power(eff_pi_plus_error_Sig,2) + TMath::Power(eff_K_plus_error_Sig,2) + TMath::Power(eff_pi_minus_error_Sig,2) + TMath::Power(eff_K_plus_fromDs_error_Sig,2) + TMath::Power(eff_K_minus_fromDs_error_Sig,2) + TMath::Power(eff_pi_minus_fromDs_error_Sig,2));

	//compute weighted efficiency ratio, using the MC weights
	weighted_efficiency_Sig = weighted_efficiency_Sig + (weight_Sig * trackingWeight_Sig);
	weighted_efficiency_error_Sig = weighted_efficiency_error_Sig + trackingWeight_error_Sig;

	//conuter for normalization
	normalization_Sig++;
	}

TCanvas* c= new TCanvas();
gPad->SetLogx(1);
effHisto->Draw("COLZ");
if(sevenTeV)c->Print("eps/TrackingEff_MCvData_11.eps");
if(!sevenTeV)c->Print("eps/TrackingEff_MCvData_12.eps");
gPad->SetLogx(0);

cout << "****************************************************************************************************" << endl;

//cout << "weighted tracking eff ratio for normalization channel:  " << weighted_efficiency_Norm << endl;
//cout << "weighted tracking eff ratio error for normalization channel:  " << weighted_efficiency_error_Norm << endl;
//cout << "normalization for normalization channel:  " << normalization_Norm << endl;
cout << "normalized weighted tracking eff ratio for normalization channel:  " << weighted_efficiency_Norm / normalization_Norm << endl;
cout << "normalized weighted tracking eff ratio error for normalization channel:  " << weighted_efficiency_error_Norm / normalization_Norm << endl;

cout << "****************************************************************************************************" << endl;

//cout << "weighted tracking eff ratio for signal channel:  " << weighted_efficiency_Sig << endl;
//cout << "normalization for signal channel:  " << normalization_Sig << endl;
cout << "normalized weighted tracking eff ratio for signal channel:  " << weighted_efficiency_Sig / normalization_Sig << endl;
cout << "normalized weighted tracking eff ratio error for signal channel:  " << weighted_efficiency_error_Sig / normalization_Sig << endl;

cout << "****************************************************************************************************" << endl;
cout << "****************************************************************************************************" << endl;

cout << "Systematic uncertainty due to tracking efficiency:  " << TMath::Abs(1. - ((weighted_efficiency_Norm / normalization_Norm)/ (weighted_efficiency_Sig / normalization_Sig))) << endl;
}

void getDecayTimeAcceptance(string Branch,string TitleX, int bins, double min, double max){

   ///Load files
   TChain* tree=new TChain("DecayTree");
   tree->Add("/auto/data/kecke/B2DPiPiPi/Data2012/data_Bs2Dspipipi_12_final_sweight.root");
   tree->SetBranchStatus("*",0);
   tree->SetBranchStatus(Branch.c_str(),1);
   tree->SetBranchStatus("N_Bs_sw",1);
   double var;
   double sw;
   tree->SetBranchAddress(Branch.c_str(),&var);
   tree->SetBranchAddress("N_Bs_sw",&sw);

   //TString fileNameMC;
   //if(useWeights)fileNameMC="/auto/data/kecke/B2DKPiPi/MC2011/mc_Ds2KKpi_BDT_reweighted.root";
   //else fileNameMC= "/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_forBDT.root";
   //fileNameMC= "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root";

   ///generate true lifetime distribution
   RooRealVar true_tau("true_tau","#tau",-0.64935,"ps^{-1}");
   RooRealVar time("time","time",0.,10.,"ps");
   RooExponential decayTime("decayTime","decayTime",time,true_tau);

   RooDataSet *generatorSet = decayTime.generate(RooArgSet(time),1000000);
   TH1D* h_MC = new TH1D("#tau toy genrator level", ";#tau [ps];Entries",bins,min,max);
   generatorSet->fillHistogram(h_MC,time);


   ///Make histograms
   TH1::SetDefaultSumw2();
   TH2::SetDefaultSumw2();
   TString title= ";"+TitleX+";Yield [norm.]";
   TH1D* h= new TH1D(Branch.c_str(),title,bins,min,max);

   // getRandom check
   TH1D* h_random = new TH1D("random check alta",title,bins,min,max);

   double ct = 0;

    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
	if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);
	ct = var * 1000;
	h->Fill(ct,sw);
    }

    ///Plot it

    TCanvas* c= new TCanvas();

    h->Scale(1./h->Integral());
    h_MC->Scale(1./h_MC->Integral());
    double maxY= h->GetMaximum();
    if(h_MC->GetMaximum()>maxY)maxY=h_MC->GetMaximum();
    h->SetMinimum(0.);
    h->SetMaximum(maxY*1.2);
    h->SetLineColor(kBlack);
    h->Draw("");
    h_MC->SetLineColor(kRed);
    h_MC->Draw("same");

    c->Print("eps/DataVStheory_tau_from_Bs2Ds3pi.eps");

    TFile* output=new TFile("/auto/data/kecke/B2DKPiPi/timeAcc/timeAcc_from_Bs2Ds3pi_var.root","RECREATE");
    TTree* summary_tree = new TTree("DecayTree","DecayTree");
    double timeAcc;
    summary_tree->Branch("timeAcc",&timeAcc,"timeAcc/D");
    TH1D *h_weight = (TH1D*)h->Clone();
    h_weight->SetName((Branch+"_acceptance").c_str());
    h_weight->Divide(h,h_MC);
    h_weight->SetMinimum(0.);
    h_weight->SetMaximum(h_weight->GetMaximum()*1.2);
    h_weight->GetYaxis()->SetTitle("Acceptance");
    h_weight->Draw();
    c->Print(("eps/"+Branch+"_acceptance.eps").c_str());
   // h_weight->Write();
   // output->Write();
   // output->Close();


	for( int i = 0; i < 1000000; i++){
		timeAcc = h_weight->GetRandom();
		h_random->Fill(timeAcc);
		summary_tree->Fill();
	}

   h_random->Draw("e1"); c->Print("eps/randomCheck.eps");
   summary_tree->Write();
   output->Write();
   output->Close();


}


void getBR(){

/// acceptance efficiency

double epsilon_acc_11_Sig = 0.1137;
double epsilon_acc_11_Sig_error = 0.0002;

double epsilon_acc_12_Sig = 0.1163;
double epsilon_acc_12_Sig_error = 0.0002;

double epsilon_acc_11_Norm = 0.1066;
double epsilon_acc_11_Norm_error = 0.0002;

double epsilon_acc_12_Norm = 0.1090; 
double epsilon_acc_12_Norm_error = 0.0002;


/// selection efficiency

double epsilon_sel_11_Sig = 0.0107; 
double epsilon_sel_11_Sig_error = 0.0001; 

double epsilon_sel_12_Sig = 0.0095; 
double epsilon_sel_12_Sig_error = 0.0001; 

double epsilon_sel_11_Norm = 0.0112; 
double epsilon_sel_11_Norm_error = 0.0001; 

double epsilon_sel_12_Norm = 0.0096; 
double epsilon_sel_12_Norm_error = 0.0001; 


/// pid efficiency

double epsilon_pid_11_Sig = 0.7325;
double epsilon_pid_11_Sig_error = 0.0088; 
 
double epsilon_pid_12_Sig = 0.7196; 
double epsilon_pid_12_Sig_error = 0.0090; 

double epsilon_pid_11_Norm = 0.8850; 
double epsilon_pid_11_Norm_error = 0.0059; 

double epsilon_pid_12_Norm = 0.8839; 
double epsilon_pid_12_Norm_error = 0.0059; 


/// compute overall efficiency

double epsilon_2011_Sig = epsilon_acc_11_Sig * epsilon_sel_11_Sig * epsilon_pid_11_Sig;
double epsilon_2012_Sig = epsilon_acc_12_Sig * epsilon_sel_12_Sig * epsilon_pid_12_Sig;

double epsilon_2011_Sig_error = TMath::Sqrt( TMath::Power((epsilon_acc_11_Sig_error/epsilon_acc_11_Sig),2) + TMath::Power((epsilon_sel_11_Sig_error/epsilon_sel_11_Sig),2) + TMath::Power((epsilon_pid_11_Sig_error/epsilon_pid_11_Sig),2)) * epsilon_2011_Sig;
double epsilon_2012_Sig_error = TMath::Sqrt( TMath::Power((epsilon_acc_12_Sig_error/epsilon_acc_12_Sig),2) + TMath::Power((epsilon_sel_12_Sig_error/epsilon_sel_12_Sig),2) + TMath::Power((epsilon_pid_12_Sig_error/epsilon_pid_12_Sig),2)) * epsilon_2012_Sig;

double epsilon_2011_Norm = epsilon_acc_11_Norm * epsilon_sel_11_Norm * epsilon_pid_11_Norm;
double epsilon_2012_Norm = epsilon_acc_12_Norm * epsilon_sel_12_Norm * epsilon_pid_12_Norm;

double epsilon_2011_Norm_error = TMath::Sqrt( TMath::Power((epsilon_acc_11_Norm_error/epsilon_acc_11_Norm),2) + TMath::Power((epsilon_sel_11_Norm_error/epsilon_sel_11_Norm),2) + TMath::Power((epsilon_pid_11_Norm_error/epsilon_pid_11_Norm),2)) * epsilon_2011_Norm;
double epsilon_2012_Norm_error = TMath::Sqrt( TMath::Power((epsilon_acc_12_Norm_error/epsilon_acc_12_Norm),2) + TMath::Power((epsilon_sel_12_Norm_error/epsilon_sel_12_Norm),2) + TMath::Power((epsilon_pid_12_Norm_error/epsilon_pid_12_Norm),2)) * epsilon_2012_Norm;


///get the yields 

double N_Bs2DsKpipi_11 = 351.;
double N_Bs2DsKpipi_11_error = 26.;

double N_Bs2DsKpipi_12 = 858.;
double N_Bs2DsKpipi_12_error = 40.;

double N_Bs2Dspipipi_11 = 7671.;
double N_Bs2Dspipipi_11_error = 96.;

double N_Bs2Dspipipi_12 = 17379.;
double N_Bs2Dspipipi_12_error = 148.;


///get the systematic errors 

double syst_pid = 0.004;
double syst_fit = 0.04;
double syst_sel = 0.009;
double syst_trigger = 0.001;
double syst_tracking = 0.015;
double syst_BDTG = 0.019;
double syst_MC = 0.013;

double syst_combined = TMath::Sqrt(TMath::Power(syst_pid,2) + TMath::Power(syst_fit,2) + TMath::Power(syst_sel,2) + TMath::Power(syst_trigger,2) + TMath::Power(syst_tracking,2) + TMath::Power(syst_BDTG,2) + TMath::Power(syst_MC,2));


///compute the BR

double BR_11 = (N_Bs2DsKpipi_11/N_Bs2Dspipipi_11) * (epsilon_2011_Norm/epsilon_2011_Sig);
double BR_11_error_stat = TMath::Sqrt(TMath::Power((N_Bs2DsKpipi_11_error/N_Bs2DsKpipi_11),2) + TMath::Power((N_Bs2Dspipipi_11_error/N_Bs2Dspipipi_11),2)) * BR_11;

double BR_12 = (N_Bs2DsKpipi_12/N_Bs2Dspipipi_12) * (epsilon_2012_Norm/epsilon_2012_Sig);
double BR_12_error_stat = TMath::Sqrt(TMath::Power((N_Bs2DsKpipi_12_error/N_Bs2DsKpipi_12),2) + TMath::Power((N_Bs2Dspipipi_12_error/N_Bs2Dspipipi_12),2)) * BR_12;

// combine 7 & 8 TeV branching ratios
double BR = (BR_11 + BR_12) / 2.;
double BR_error_stat = (BR_11_error_stat + BR_12_error_stat) / 2.;
double BR_error_syst = BR * syst_combined;

///print efficiencies

cout << "total efficiency for 2011 normalization channel: " << "e=  " << epsilon_2011_Norm*100 << "+/-" << epsilon_2011_Norm_error*100 << " %" << endl;
cout << "total efficiency for 2012 normalization channel: " << "e=  " << epsilon_2012_Norm*100 << "+/-" << epsilon_2012_Norm_error*100 << " %" << endl;

cout << "total efficiency for 2011 signal channel: " << "e=  " << epsilon_2011_Sig*100 << "+/-" << epsilon_2011_Sig_error*100 << " %" << endl;
cout << "total efficiency for 2012 signal channel: " << "e=  " << epsilon_2012_Sig*100 << "+/-" << epsilon_2012_Sig_error*100 << " %" << endl;


cout << "**********************************" << endl;
cout << "**********************************" << endl;
cout << "**********************************" << endl;

cout << "Branching Ratio:    " << BR << "  +/-  " << BR_error_stat << "  +/-  " << BR_error_syst << endl;

}


void getDecayTimeError(){

///Load files
TChain* tree=new TChain("DecayTree");
tree->Add("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");

tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_TAU",1);
tree->SetBranchStatus("Bs_cterr",1);
tree->SetBranchStatus("Bs_TRUETAU",1);


Double_t Bs_TAU;
Double_t Bs_TRUETAU;
Double_t Bs_cterr;

tree -> SetBranchAddress( "Bs_TAU" , &Bs_TAU );
tree -> SetBranchAddress( "Bs_TRUETAU" , &Bs_TRUETAU );
tree -> SetBranchAddress( "Bs_cterr" , &Bs_cterr );

TFile* output = 0;
output = new TFile("/auto/data/kecke/B2DKPiPi/TimeRes/Bs2DsKpipi_MCcombined_9to19.root","RECREATE");

TTree* summary_tree = tree->CloneTree(0);

double Bs_DeltaTau = 0;
summary_tree->Branch("Bs_DeltaTau",&Bs_DeltaTau,"Bs_DeltaTau/D");

double Sum_cterr = 0;
int Nevents = 0;

//loop over events
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	if(Bs_TAU < 0. || Bs_TRUETAU < 0.) continue;
	if( Bs_cterr > 0.009 && Bs_cterr < 0.019){
		Bs_DeltaTau = (1000*Bs_TAU) - (1000*Bs_TRUETAU);
		summary_tree->Fill();
		}

	Sum_cterr = Sum_cterr + Bs_cterr;
	Nevents++;
	}


cout << "average sigma_t:    " <<  (Sum_cterr/Nevents)*1000 << " fs" << endl;


summary_tree->Write();
output->Close();
}


void FitTimeRes(string ResoBin, string BinName){

string Directory = "/auto/data/kecke/B2DKPiPi/TimeRes/";

string FullName = Directory.append(ResoBin);

TFile* file= new TFile(FullName.c_str());
TTree* tree= (TTree*) file->Get("DecayTree");

RooRealVar Bs_DeltaTau("Bs_DeltaTau", "#Delta(t)", -0.2, 0.2,"[ps]");
RooArgList list =  RooArgList(Bs_DeltaTau);

RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_DeltaTau),Import(*tree));


RooRealVar meanBs1("meanBs1", "B_{s} #mu", 0., -0.2, 0.2);

RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 0.010,0.,0.2);
RooRealVar sigmaBs2("sigmaBs2", "B_{s} #sigma_{2}", 0.045,0.,0.2);
RooRealVar sigmaBs3("sigmaBs3", "B_{s} #sigma_{3}", 0.040,0.,0.2);

RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_DeltaTau, meanBs1, sigmaBs1);
RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_DeltaTau, meanBs1, sigmaBs2);
RooGaussian GaussBs3("GaussBs3", "GaussBs3", Bs_DeltaTau, meanBs1, sigmaBs3);
RooGaussian GaussBs("GaussBs", "GaussBs", Bs_DeltaTau, meanBs1, sigmaBs1);
RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0., 1.);
RooRealVar f_GaussBs2("f_GaussBs2" , "2f__{B_{s}}", 0.5, 0., 1.);
RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
RooAddPdf TripleGaussBs("TripleGaussBs", "TripleGaussBs", RooArgList(GaussBs1,GaussBs2,GaussBs3),RooArgList(f_GaussBs,f_GaussBs2));

RooFitResult *result;
result = DoubleGaussBs.fitTo(*data,Save(kTRUE),Extended(kFALSE),NumCPU(3));
cout << "result is --------------- "<<endl;
result->Print();


double range_dw = Bs_DeltaTau.getMin();
double range_up = Bs_DeltaTau.getMax();
int bin = 75;


RooPlot* frame_m= Bs_DeltaTau.frame();
frame_m->SetTitle("");

frame_m->GetXaxis()->SetLabelSize( 0.06 );
frame_m->GetYaxis()->SetLabelSize( 0.06 );
frame_m->GetXaxis()->SetLabelFont( 132 );
frame_m->GetYaxis()->SetLabelFont( 132 );
frame_m->GetXaxis()->SetLabelOffset( 0.006 );
frame_m->GetYaxis()->SetLabelOffset( 0.006 );
frame_m->GetXaxis()->SetLabelColor( kWhite);

frame_m->GetXaxis()->SetTitleSize( 0.06 );
frame_m->GetYaxis()->SetTitleSize( 0.06 );
frame_m->GetYaxis()->SetNdivisions(512);
frame_m->GetXaxis()->SetTitleOffset( 1.00 );
frame_m->GetYaxis()->SetTitleOffset( 1.00 );


data->plotOn(frame_m,Name("dataSetCut"),MarkerSize(0.5),Binning(bin));
DoubleGaussBs.plotOn(frame_m,Name("FullPdf"),LineColor(kBlack),LineWidth(2));
DoubleGaussBs.plotOn(frame_m,Components(GaussBs1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
DoubleGaussBs.plotOn(frame_m,Components(GaussBs2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));

TLatex* lhcbtext = new TLatex();
lhcbtext->SetTextFont(132);
lhcbtext->SetTextColor(1);
lhcbtext->SetTextSize(0.07);
lhcbtext->SetTextAlign(12);


TCanvas* canvas = new TCanvas("canvas", "canvas", 1200, 800);
canvas->cd();
canvas->SetLeftMargin(0.01);
canvas->SetRightMargin(0.01);
canvas->SetTopMargin(0.05);
canvas->SetBottomMargin(0.05);
TPad* pad1 = new TPad("upperPad", "upperPad", .050, .22, 1.0, 1.0);
pad1->SetBorderMode(0);
pad1->SetBorderSize(-1);
pad1->SetFillStyle(0);
pad1->SetTickx(0);
pad1->SetLeftMargin(0.115);
pad1->SetRightMargin(0.05);
pad1->Draw();
pad1->cd();
frame_m->GetYaxis()->SetRangeUser(0.1,frame_m->GetMaximum()*1.0);
frame_m->Draw();

lhcbtext->DrawTextNDC(0.70,0.85,"LHCb simulation");
lhcbtext->DrawTextNDC(0.70,0.75,"preliminary");
lhcbtext->DrawTextNDC(0.70,0.65,"Bin(19 - 24) fs");

canvas->cd();
TPad* pad2 = new TPad("lowerPad", "lowerPad", .050, .005, 1.0, .3275);
pad2->SetBorderMode(0);
pad2->SetBorderSize(-1);
pad2->SetFillStyle(0);
pad2->SetBottomMargin(0.35);
pad2->SetLeftMargin(0.115);
pad2->SetRightMargin(0.05);
pad2->SetTickx(0);
pad2->Draw();
pad2->SetLogy(0);
pad2->cd();


frame_m->Print("v");
RooPlot* frame_p = Bs_DeltaTau.frame(Title("pull_frame"));
frame_p->Print("v");
frame_p->SetTitle("");
frame_p->GetYaxis()->SetTitle("");
frame_p->GetYaxis()->SetTitleSize(0.09);
frame_p->GetYaxis()->SetTitleOffset(0.26);
frame_p->GetYaxis()->SetTitleFont(62);
frame_p->GetYaxis()->SetNdivisions(106);
frame_p->GetYaxis()->SetLabelSize(0.12);
frame_p->GetYaxis()->SetLabelOffset(0.006);
frame_p->GetXaxis()->SetTitleSize(0.15);
frame_p->GetXaxis()->SetTitleFont(132);
frame_p->GetXaxis()->SetTitleOffset(0.85);
frame_p->GetXaxis()->SetNdivisions(515);
frame_p->GetYaxis()->SetNdivisions(5);
frame_p->GetXaxis()->SetLabelSize(0.12);
frame_p->GetXaxis()->SetLabelFont( 132 );
frame_p->GetYaxis()->SetLabelFont( 132 );
frame_p->GetXaxis()->SetTitle("#font[132]{#Delta t(B_{s}) [ps]}");

TString* obsTS = new TString(Bs_DeltaTau.GetName());
TString* pullnameTS = new TString("FullPdf");
TString* pullname2TS = new TString("dataSetCut");
RooHist* pullHist  = frame_m->pullHist(pullname2TS->Data(),pullnameTS->Data());
frame_p->addPlotable(pullHist,"P");

double chi2 = frame_m->chiSquare();
double chi22 = frame_m->chiSquare(pullnameTS->Data(),pullname2TS->Data());



TAxis* axisX = pullHist->GetXaxis();
RooBinning* Bin = new RooBinning(range_dw,range_up,"P");
Bin->addUniform(bin, range_dw, range_up);
axisX->Set(Bin->numBins(), Bin->array());

TAxis* axisY = pullHist->GetYaxis();
double max = 5.0 ;
double min = -5.0 ;
axisY->SetLabelSize(0.12);
axisY->SetNdivisions(5);
axisX->SetLabelSize(0.12);

double rangeX = max-min;
double zero = max/rangeX;

TGraph* graph = new TGraph(2);
graph->SetMaximum(max);
graph->SetMinimum(min);
graph->SetPoint(1,range_dw,0);
graph->SetPoint(2,range_up,0);

TGraph* graph2 = new TGraph(2);
graph2->SetMaximum(max);
graph2->SetMinimum(min);
graph2->SetPoint(1,range_dw,-3);
graph2->SetPoint(2,range_up,-3);
graph2->SetLineColor(kRed);

TGraph* graph3 = new TGraph(2);
graph3->SetMaximum(max);
graph3->SetMinimum(min);
graph3->SetPoint(1,range_dw,3);
graph3->SetPoint(2,range_up,3);
graph3->SetLineColor(kRed);

pullHist->GetXaxis()->SetLabelFont( 132 );
pullHist->GetYaxis()->SetLabelFont( 132 );
pullHist->SetTitle("");

frame_p->GetYaxis()->SetRangeUser(-5.0,5.0);
frame_p->Draw();

graph->Draw("same");
//graph2->Draw("same");
//graph3->Draw("same");

pad2->Update();
canvas->Update();

string epsName = "eps/TimeRes/";
string epsFullName = epsName.append(BinName);

canvas->SaveAs(epsFullName.c_str());


double f1 = f_GaussBs.getVal();
double df1 = f_GaussBs.getError();
double f2 = 1. - f1;
double df2 = f_GaussBs.getError();

double sig1 = sigmaBs1.getVal();
double dsig1 = sigmaBs1.getError();
double sig2 = sigmaBs2.getVal();
double dsig2 = sigmaBs2.getError();

///compute resolution from combination of Gauss sigmas
double resolution = TMath::Sqrt( (f1*sig1*sig1) + (f2*sig2*sig2) );
double dresolution = TMath::Sqrt(TMath::Power(((sig1*f1*dsig1)/resolution),2)+TMath::Power(((sig2*f2*dsig2)/resolution),2)+TMath::Power(((sig1*sig1*df1)/(2*resolution)),2)+TMath::Power(((sig2*sig2*df2)/(2*resolution)),2));

cout << "Measured resolution from GaussComb:   " << resolution*1000 << " +/- " << dresolution*1000 <<" fs" << endl; 
cout << "***********************************************************************" << endl;



///compute resolution from dilution
double dms = 17.757; // [1/ps]

double dilution = f1*TMath::Exp(-1*(sig1*sig1*dms*dms)/2.) + f2*TMath::Exp(-1*(sig2*sig2*dms*dms)/2.);

double Term1 = f1*TMath::Exp(-1*(sig1*sig1*dms*dms)/2.);
double Term2 = f2*TMath::Exp(-1*(sig2*sig2*dms*dms)/2.);

double ddilution = TMath::Sqrt( TMath::Power(df1*(Term1/f1),2) + TMath::Power(df2*(Term2/f2),2) + TMath::Power((Term1*dms*dms*sig1)*dsig1,2) + TMath::Power((Term2*dms*dms*sig2)*dsig2,2) );

double resolution_eff = TMath::Sqrt((-2/(dms*dms))*log(dilution));
double dresolution_eff = ((2/(dms*dms))/(2*dilution*TMath::Sqrt((-2/(dms*dms))*log(dilution)))) * ddilution;


cout << "Measured resolution from dilution:   " << resolution_eff*1000 << " +/- " << dresolution_eff*1000 <<" fs" << endl;


file->Close();

}


void FitResoRelation(){

//define sigma bins
double sigma_t_BinCenter[9] = {0., 19. , 24., 29. , 34. , 39. , 44. , 49., 150.};
//double sigma_t_BinWidth[8] = {9.5 , 2.5 , 2.5 , 2.5, 2.5 , 2.5 , 2.5 , 50.5};

//define reso t bins
double reso_t_BinCenter[8] = {27.4569 , 30.639 , 34.6615 , 39.0903 , 44.7556 , 50.9845 ,  54.8914 , 68.9104};
double reso_t_BinWidth[8] = {8.81782 , 8.48237 , 5.27469 , 10.4688 , 11.7751 , 15.1059 , 1.6030 ,12.5501};


//define polynom

//TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 45.);
TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 150.);
fitFunc->SetParNames("c0","s");
fitFunc->SetParameters(10.,1.2);
fitFunc->SetParLimits(0,0.,30.);
fitFunc->SetParLimits(1,0.,10.);

fitFunc->FixParameter(0,0.);
//fitFunc->FixParameter(1,1.280);

//fill histo
TH1D* ResoRelation = new TH1D("ResoRelation" ,"ResoRelation", 8, sigma_t_BinCenter);

for(int i = 0; i < 8; i++)
	{
	ResoRelation->SetBinContent((i+1), reso_t_BinCenter[i]);
	ResoRelation->SetBinError((i+1), reso_t_BinWidth[i]);
	}

ResoRelation->SetXTitle("decay time error #sigma(t) [fs]");
ResoRelation->SetYTitle("decay time resolution #Delta(t) [fs]");

ResoRelation->SetAxisRange(0.,180.);
ResoRelation->SetAxisRange(0.,120.,"Y");

TCanvas* canvas = new TCanvas();


ResoRelation->Fit(fitFunc,"RL");

ResoRelation->Draw("e1");
fitFunc->Draw("same");
canvas->Print("eps/TimeRes/ResoRelation_MC.pdf");

}

int main(){
    time_t startTime = time(0);


//	FitResoRelation();
//	getDecayTimeError();
	//taggingCategories();
	//getBR();
//	FitTimeRes("Bs2DsKpipi_MCcombined_19to24.root", "SignalMC_19to24.pdf");
	//getExpBkgShape();


	//applyBDTcut("0.7012");
/*
	bool sevenTeV=true;
	bool fitSimultan=true;
	bool Ds2KKpi = false;
	bool combineYears = false;
	double fillarr4[13];
	double *DsstarKpipifromNorm;
	DsstarKpipifromNorm = fitBDTNorm(fillarr4,sevenTeV,fitSimultan,Ds2KKpi,combineYears);
*/

   //getDecayTimeAcceptance("Bs_TAU","#tau_{Bs}",10,0.,10.);
 // addVars();
//makeSweightedPlots();
 //quickFit(true,true);
   //createRooDataSet("mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14",false);
 // TriggerSyst(false);
 // reweight(false, true);
//TrackingEff();

 // quickSignalEstimate();
  //  iterateBDT(-0.9,0.9,0.05);
    //preselect();

   //nominal value for Ds->KKpi final state
   // applyBDTcut("0.7012");

    checkForDuplicates();

   ///nominal value for Ds->KKpi final state, DsK ,3 fb^-1 selection
  //applyBDTcut("0.6891");

   //troll value for Ds->KKpi final state
 // applyBDTcut("0.9082");

//   makePlots();
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
 //  getBR();
  //  MCStudies();
   //splitMC();
// bool sevenTeV , bool signalMode
  //DataStructure(true, true);
 // Systematics(false,true);


    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;

return 0;
}

