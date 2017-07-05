//offline selection of BDT variables for Bs->DsKpipi
//philippe d'argent & matthieu kecke

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
#include <THStack.h>
#include "TRandom3.h"
#include <sstream>


using namespace std;

int preselect() {

    //specifiy decay channel of Ds
    bool Ds2KKpi = true;
    bool Ds2Kpipi = false;
    bool Ds2pipipi = false;

    //MC or Data?
    bool MC = false;

    //Bkg sample?
    bool Ds3piBkg = false;
    bool Dsstar3piBkg = false;

    //which year?
    bool sevenTeV = false; //for 2011 = true , for 2012 = false 

    TChain* tree = 0;
    tree=new TChain("Bs2DsKpipi_Ds2KKpi_Tuple/DecayTree");
    //TString loc = "/auto/data/dargent/gangadir/workspace/dargent/LocalXML/55/";
    //TString file = "/output/b2dhhh.root";
    tree->Add("/auto/data/kecke/B2DKPiPi/12/*.root");
/*
cout << "yolo" << endl;
    for(int i = 0; i<1000; i++){
            TString dir_i = TString::Format("%d",i);
            tree->Add(loc + dir_i + file);
    }
cout << "yolo oho" << endl;
*/	
    //Ds2KKpi case
/// CAREFULL ! ONLY TEMPORARELY DISABLE
/*
    if(Ds2KKpi){
	tree=new TChain("Bs2DsKpipi_Ds2KKpi_Tuple/DecayTree");
	if((!MC) && (sevenTeV)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U/*.root");
	}
	else if((!MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-PID/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-PID/*.root");
	}
	else if((MC) && (Ds3piBkg)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Ds3pi-U/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Ds3pi-D/*.root");
	}
	else if((MC) && (Dsstar3piBkg)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dst3pi-U/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dst3pi-D/*.root");
	}
	else if((MC) && (!sevenTeV)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/12D/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/12U/*.root");
	}
	else if((MC) && (sevenTeV)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/11U/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/11D/*.root");
	}
    }
*/
    //Ds2Kpipi case
    if(Ds2Kpipi){
    	tree=new TChain("Bs2DsKpipi_Ds2Kpipi_Tuple/DecayTree");
	if((!MC) && (sevenTeV)){
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U/*.root");
	}
	if((!MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-PID/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-PID/*.root");
	}
    }

    //Ds2pipipi case
    if(Ds2pipipi){
    	tree=new TChain("Bs2DsKpipi_Ds2pipipi_Tuple/DecayTree");
	if((!MC) && (sevenTeV)){
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U/*.root");
	}
	else if((!MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-d23pi/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-d23pi/*.root");
	}
	else if((MC) && (!sevenTeV)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/12D/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/12U/*.root");
	}
	else if((MC) && (sevenTeV)){
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/11U/*.root");
		tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Signal/11D/*.root");
	}
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
    tree->SetBranchStatus( "*TAU*", 1 );
    tree->SetBranchStatus( "*ERR*", 1 );
    tree->SetBranchStatus( "*DTF*", 1 );    
    tree->SetBranchStatus( "*ETA", 1 );
    tree->SetBranchStatus( "*FD*", 1 );
    tree->SetBranchStatus( "*IP*", 1 );
    tree->SetBranchStatus( "*TAG*", 1);
    tree->SetBranchStatus("*_SS_*",1);
    tree->SetBranchStatus("*_OS_*",1);
    tree->SetBranchStatus("*VtxCharge*",1);
    tree->SetBranchStatus("*DEC*",1);
    tree->SetBranchStatus("*PROB*",1);

    tree->SetBranchStatus( "*IPCHI2_OWNPV", 1 );
    tree->SetBranchStatus( "*FDCHI2_OWNPV",1 );
    tree->SetBranchStatus( "*OWNPV*", 1 );
    tree->SetBranchStatus( "*ORIVX*", 1 );
    tree->SetBranchStatus( "*DIRA_OWNPV",1);
    tree->SetBranchStatus( "*ENDVERTEX_CHI2",1 );
    tree->SetBranchStatus( "*ENDVERTEX*",1 );
    tree->SetBranchStatus( "*DOCA*",1 );
    tree->SetBranchStatus("*TRACK*",1) ;
    tree->SetBranchStatus("Bs_L0Global_TIS",1) ;

    tree->SetBranchStatus( "*chi2",1 );
    tree->SetBranchStatus( "*nDOF",1 );
    tree->SetBranchStatus( "*status",1 );
    tree->SetBranchStatus( "*ctau",1 );
    tree->SetBranchStatus( "*ID*",1 );

    if(!MC)tree->SetBranchStatus("*TRUE*",0) ;
    else {
    	tree->SetBranchStatus("*TRUE*",1) ;
    	tree->SetBranchStatus("Bs_BKGCAT",1);
    }
	

//define variables
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;
TLorentzVector K_plus;
TLorentzVector pi_plus;
TLorentzVector pi_minus;

//Ds->Kpipi case
TLorentzVector pi_plus_fromDs;

//Ds->pipipi case
TLorentzVector pi_minus2_fromDs;

TLorentzVector Kminus_asPiminus_MissID;
TLorentzVector Kminus_asProton_MissID;

double massKaon = 493.68;
double massPion = 139.57;
double massProton = 938.27;
double massPhi = 1019.46;
double massKstar = 895.81;
double massDs = 1968.30;
double massDminus = 1869.61;
double massLambda_c = 2286.46;
//double massJPsi =  3096.92;

Double_t Bs_MM;
Double_t Ds_MM;
Int_t Bs_BKGCAT;
Double_t K_plus_fromDs_PX;
Double_t K_plus_fromDs_PY;
Double_t K_plus_fromDs_PZ;
Double_t K_plus_PX;
Double_t K_plus_PY;
Double_t K_plus_PZ;
Double_t K_minus_fromDs_PX;
Double_t K_minus_fromDs_PY;
Double_t K_minus_fromDs_PZ;
Double_t pi_minus_fromDs_PX;
Double_t pi_minus_fromDs_PY;
Double_t pi_minus_fromDs_PZ;
Double_t pi_minus_PX;
Double_t pi_minus_PY;
Double_t pi_minus_PZ;
Double_t pi_plus_PX;
Double_t pi_plus_PY;
Double_t pi_plus_PZ;

// Ds->Kpipi case------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Double_t pi_plus_fromDs_PX;
Double_t pi_plus_fromDs_PY;
Double_t pi_plus_fromDs_PZ;
Double_t pi_plus_fromDs_PT;
Double_t pi_plus_fromDs_P;
Double_t pi_plus_fromDs_PIDK;
Double_t pi_plus_fromDs_IPCHI2_OWNPV;
Double_t pi_plus_fromDs_TRACK_CHI2NDOF;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Ds->pipipi case------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Double_t pi_minus2_fromDs_PX;
Double_t pi_minus2_fromDs_PY;
Double_t pi_minus2_fromDs_PZ;
Double_t pi_minus2_fromDs_PT;
Double_t pi_minus2_fromDs_P;
Double_t pi_minus2_fromDs_PIDK;
Double_t pi_minus2_fromDs_IPCHI2_OWNPV;
Double_t pi_minus2_fromDs_TRACK_CHI2NDOF;
Double_t pi_plus_fromDs_TRACK_GhostProb;
Double_t pi_minus2_fromDs_TRACK_GhostProb;
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Double_t K_plus_P;
Double_t pi_plus_P;
Double_t pi_minus_P;
Double_t pi_minus_fromDs_P;
Double_t K_minus_fromDs_P;
Double_t K_plus_fromDs_P;


Double_t K_plus_fromDs_PIDK;
Double_t K_minus_fromDs_PIDK;
Double_t pi_minus_fromDs_PIDK;
Double_t K_plus_PIDK;
Double_t pi_minus_PIDK;
Double_t pi_plus_PIDK;
Double_t K_minus_fromDs_PIDp;

Double_t Ds_FDCHI2_ORIVX;
Double_t Ds_FDCHI2_OWNPV;

//BDT variables
Double_t Bs_IPCHI2_OWNPV;
Double_t K_plus_fromDs_IPCHI2_OWNPV;
Double_t K_minus_fromDs_IPCHI2_OWNPV;
Double_t pi_minus_fromDs_IPCHI2_OWNPV;
Double_t K_1_1270_plus_IPCHI2_OWNPV;
Double_t K_1_1270_plus_DOCA1;
Double_t K_1_1270_plus_DOCA;
Double_t K_1_1270_plus_DOCA2;
Double_t K_1_1270_plus_DOCA3;
Double_t K_1_1270_plus_ENDVERTEX_X;
Double_t K_1_1270_plus_ENDVERTEX_Y;
Double_t K_1_1270_plus_ENDVERTEX_Z;
Double_t K_1_1270_plus_OWNPV_X;
Double_t K_1_1270_plus_OWNPV_Y;
Double_t K_1_1270_plus_OWNPV_Z;
Double_t Ds_DOCA;
Double_t Ds_DOCA1;
Double_t Ds_DOCA2;
Double_t Ds_DOCA3;
Double_t K_plus_IPCHI2_OWNPV;
Double_t pi_plus_IPCHI2_OWNPV;
Double_t pi_minus_IPCHI2_OWNPV;
Double_t K_plus_fromDs_PT;
Double_t K_minus_fromDs_PT;
Double_t pi_minus_fromDs_PT;
Double_t K_plus_PT;
Double_t pi_plus_PT;
Double_t pi_minus_PT;
Double_t K_plus_fromDs_TRACK_CHI2NDOF;
Double_t K_minus_fromDs_TRACK_CHI2NDOF;
Double_t pi_minus_fromDs_TRACK_CHI2NDOF;
Double_t K_plus_TRACK_CHI2NDOF;
Double_t pi_plus_TRACK_CHI2NDOF;
Double_t pi_minus_TRACK_CHI2NDOF;

Double_t K_plus_fromDs_TRACK_GhostProb;
Double_t K_minus_fromDs_TRACK_GhostProb;
Double_t pi_minus_fromDs_TRACK_GhostProb;
Double_t K_plus_TRACK_GhostProb;
Double_t pi_plus_TRACK_GhostProb;
Double_t pi_minus_TRACK_GhostProb;

Double_t Ds_ENDVERTEX_CHI2;
Int_t Ds_ENDVERTEX_NDOF;
Double_t Ds_ENDVERTEX_Z;
Double_t K_1_1270_plus_ENDVERTEX_CHI2;
Int_t K_1_1270_plus_ENDVERTEX_NDOF;
Double_t K_1_1270_plus_FDCHI2_OWNPV;
Double_t Bs_OWNPV_X; 
Double_t Bs_ENDVERTEX_X;
Double_t Bs_OWNPV_Y; 
Double_t Bs_ENDVERTEX_Y;
Double_t Bs_OWNPV_Z; 
Double_t Bs_ENDVERTEX_Z;
Double_t K_1_1270_plus_DIRA_OWNPV;

Double_t Bs_DIRA_OWNPV;
Double_t Bs_FDCHI2_OWNPV;
Double_t Bs_ENDVERTEX_CHI2;
Int_t Bs_ENDVERTEX_NDOF;

Double_t Bs_TAU;
Double_t Bs_TAUERR;
Double_t Ds_TAU;
Double_t Ds_TAUERR;

Int_t Bs_TRUEID;
Int_t Ds_TRUEID;
Int_t K_plus_TRUEID;
Int_t pi_plus_TRUEID;
Int_t pi_minus_TRUEID;
Int_t K_plus_fromDs_TRUEID;
Int_t K_minus_fromDs_TRUEID;
Int_t pi_minus_fromDs_TRUEID;

Int_t pi_plus_fromDs_TRUEID;
Int_t pi_minus2_fromDs_TRUEID;

Bool_t Bs_L0Global_TIS;
Bool_t Bs_L0HadronDecision_TOS;
Bool_t Bs_Hlt1TrackAllL0Decision_TOS;
Bool_t Bs_Hlt2Topo2BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo3BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo4BodyBBDTDecision_TOS;

Double_t pi_plus_fromDs_PIDp;
Double_t pi_minus_fromDs_PIDp;
Double_t pi_minus2_fromDs_PIDp;

Int_t Bs_ID;
Int_t K_plus_ID;


//trigger decisions
tree -> SetBranchAddress( "Bs_L0Global_TIS" , &Bs_L0Global_TIS );
tree -> SetBranchAddress( "Bs_L0HadronDecision_TOS" , &Bs_L0HadronDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt1TrackAllL0Decision_TOS" , &Bs_Hlt1TrackAllL0Decision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo2BodyBBDTDecision_TOS" , &Bs_Hlt2Topo2BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo3BodyBBDTDecision_TOS" , &Bs_Hlt2Topo3BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo4BodyBBDTDecision_TOS" , &Bs_Hlt2Topo4BodyBBDTDecision_TOS );

tree -> SetBranchAddress( "Bs_BKGCAT" , &Bs_BKGCAT );
tree -> SetBranchAddress( "Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z );
tree -> SetBranchAddress( "Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF );
tree -> SetBranchAddress( "Ds_FDCHI2_ORIVX", &Ds_FDCHI2_ORIVX );
tree -> SetBranchAddress( "Ds_FDCHI2_OWNPV", &Ds_FDCHI2_OWNPV );
tree -> SetBranchAddress( "K_1_1270_plus_ENDVERTEX_CHI2", &K_1_1270_plus_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "K_1_1270_plus_ENDVERTEX_NDOF", &K_1_1270_plus_ENDVERTEX_NDOF );
tree -> SetBranchAddress( "K_1_1270_plus_FDCHI2_OWNPV", &K_1_1270_plus_FDCHI2_OWNPV );
tree -> SetBranchAddress("K_1_1270_plus_ENDVERTEX_X", &K_1_1270_plus_ENDVERTEX_X);
tree -> SetBranchAddress("K_1_1270_plus_ENDVERTEX_Y", &K_1_1270_plus_ENDVERTEX_Y);
tree -> SetBranchAddress("K_1_1270_plus_ENDVERTEX_Z", &K_1_1270_plus_ENDVERTEX_Z);
tree -> SetBranchAddress("K_1_1270_plus_OWNPV_X", &K_1_1270_plus_OWNPV_X);
tree -> SetBranchAddress("K_1_1270_plus_OWNPV_Y", &K_1_1270_plus_OWNPV_Y);
tree -> SetBranchAddress("K_1_1270_plus_OWNPV_Z", &K_1_1270_plus_OWNPV_Z);
tree -> SetBranchAddress( "Bs_OWNPV_X", &Bs_OWNPV_X ); 
tree -> SetBranchAddress( "Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X );
tree -> SetBranchAddress( "Bs_OWNPV_Y", &Bs_OWNPV_Y ); 
tree -> SetBranchAddress( "Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y );
tree -> SetBranchAddress( "Bs_OWNPV_Z", &Bs_OWNPV_Z ); 
tree -> SetBranchAddress( "Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z );
tree -> SetBranchAddress( "K_1_1270_plus_DIRA_OWNPV", &K_1_1270_plus_DIRA_OWNPV );

tree -> SetBranchAddress( "Bs_ID" , &Bs_ID );
tree -> SetBranchAddress( "K_plus_ID" , &K_plus_ID );

//Ds->KKpi case----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Ds2KKpi){
	tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
	tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
	tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
	tree -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P );
	tree -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
	tree -> SetBranchAddress( "K_plus_fromDs_TRACK_CHI2NDOF" , &K_plus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "K_plus_fromDs_IPCHI2_OWNPV" ,&K_plus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "K_plus_fromDs_PIDK" , &K_plus_fromDs_PIDK );
	tree -> SetBranchAddress( "K_plus_fromDs_TRUEID" , &K_plus_fromDs_TRUEID );

	tree -> SetBranchAddress("K_plus_fromDs_TRACK_GhostProb", &K_plus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress("K_plus_TRACK_GhostProb", &K_plus_TRACK_GhostProb );
	tree -> SetBranchAddress("pi_plus_TRACK_GhostProb", &pi_plus_TRACK_GhostProb );
	tree -> SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb );
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//overlap---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Ds2Kpipi || Ds2KKpi){
	tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
	tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
	tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
	tree -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
	tree -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P );
	tree -> SetBranchAddress( "K_minus_fromDs_IPCHI2_OWNPV" ,&K_minus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "K_minus_fromDs_TRACK_CHI2NDOF" , &K_minus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
	tree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );
	tree -> SetBranchAddress( "K_minus_fromDs_TRUEID" , &K_minus_fromDs_TRUEID );
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//overlap----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Ds2Kpipi || Ds2pipipi){
	tree -> SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
	tree -> SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
	tree -> SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
	tree -> SetBranchAddress( "pi_plus_fromDs_PT" , &pi_plus_fromDs_PT );
	tree -> SetBranchAddress( "pi_plus_fromDs_P" , &pi_plus_fromDs_P );
	tree -> SetBranchAddress( "pi_plus_fromDs_IPCHI2_OWNPV" ,&pi_plus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRACK_CHI2NDOF" , &pi_plus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "pi_plus_fromDs_PIDK" , &pi_plus_fromDs_PIDK );
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//Ds->pipipi case----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Ds2pipipi){
	tree -> SetBranchAddress( "pi_minus2_fromDs_PX" , &pi_minus2_fromDs_PX );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PY" , &pi_minus2_fromDs_PY );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PZ" , &pi_minus2_fromDs_PZ );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PT" , &pi_minus2_fromDs_PT );
	tree -> SetBranchAddress( "pi_minus2_fromDs_P" , &pi_minus2_fromDs_P );
	tree -> SetBranchAddress( "pi_minus2_fromDs_IPCHI2_OWNPV" ,&pi_minus2_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRACK_CHI2NDOF" , &pi_minus2_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PIDK" , &pi_minus2_fromDs_PIDK );
	tree -> SetBranchAddress( "pi_plus_fromDs_PIDp" , &pi_plus_fromDs_PIDp );
	tree -> SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PIDp" , &pi_minus2_fromDs_PIDp );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRUEID" , &pi_minus2_fromDs_TRUEID );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRUEID" , &pi_plus_fromDs_TRUEID );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRACK_GhostProb", &pi_plus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRACK_GhostProb", &pi_minus2_fromDs_TRACK_GhostProb );
}
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


tree -> SetBranchAddress( "Bs_DIRA_OWNPV", &Bs_DIRA_OWNPV );
tree -> SetBranchAddress( "Bs_FDCHI2_OWNPV", &Bs_FDCHI2_OWNPV );
tree -> SetBranchAddress( "Bs_IPCHI2_OWNPV", &Bs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF );

tree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
tree -> SetBranchAddress( "Ds_TRUEID" , &Ds_TRUEID );
tree -> SetBranchAddress( "K_plus_TRUEID" , &K_plus_TRUEID );
tree -> SetBranchAddress( "pi_plus_TRUEID" , &pi_plus_TRUEID );
tree -> SetBranchAddress( "pi_minus_TRUEID" , &pi_minus_TRUEID );
tree -> SetBranchAddress( "pi_minus_fromDs_TRUEID" , &pi_minus_fromDs_TRUEID );

tree -> SetBranchAddress("Bs_TAU" , &Bs_TAU );
tree -> SetBranchAddress("Bs_TAUERR" , &Bs_TAUERR );
tree -> SetBranchAddress("Ds_TAU" , &Ds_TAU );
tree -> SetBranchAddress("Ds_TAUERR" , &Ds_TAUERR );

//link variables to tree branches
tree -> SetBranchAddress( "Bs_MM" , &Bs_MM );
tree -> SetBranchAddress( "Ds_MM" , &Ds_MM );
tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
tree -> SetBranchAddress( "K_plus_PX" , &K_plus_PX );
tree -> SetBranchAddress( "K_plus_PY" , &K_plus_PY );
tree -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );
tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
tree -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
tree -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
tree -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );

tree -> SetBranchAddress( "pi_plus_P" , &pi_plus_P );
tree -> SetBranchAddress( "pi_minus_P" , &pi_minus_P );
tree -> SetBranchAddress( "K_plus_P" , &K_plus_P );
tree -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P );

//BDT variables
tree -> SetBranchAddress( "Bs_IPCHI2_OWNPV" ,&Bs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_fromDs_IPCHI2_OWNPV" ,&pi_minus_fromDs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "K_1_1270_plus_IPCHI2_OWNPV" ,&K_1_1270_plus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "K_plus_IPCHI2_OWNPV" ,&K_plus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus_IPCHI2_OWNPV" ,&pi_plus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_IPCHI2_OWNPV" ,&pi_minus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );
tree -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
tree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
tree -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
if((Ds2Kpipi && sevenTeV) || Ds2KKpi || (Ds2pipipi && sevenTeV))  tree -> SetBranchAddress( "K_1_1270_plus_DOCA1" ,&K_1_1270_plus_DOCA1 );
if(Ds2pipipi && (!sevenTeV))  tree -> SetBranchAddress( "K_1_1270_plus_DOCA" ,&K_1_1270_plus_DOCA );
if(Ds2Kpipi && (!sevenTeV))  tree -> SetBranchAddress( "K_1_1270_plus_DOCA" ,&K_1_1270_plus_DOCA );
tree -> SetBranchAddress( "K_1_1270_plus_DOCA2" ,&K_1_1270_plus_DOCA2 );
tree -> SetBranchAddress( "K_1_1270_plus_DOCA3" ,&K_1_1270_plus_DOCA3 );
if((Ds2Kpipi && sevenTeV) || Ds2KKpi || (Ds2pipipi && sevenTeV))   tree -> SetBranchAddress( "Ds_DOCA1" ,&Ds_DOCA1);
if(Ds2pipipi && (!sevenTeV))tree -> SetBranchAddress( "Ds_DOCA" ,&Ds_DOCA);
if(Ds2Kpipi && (!sevenTeV))tree -> SetBranchAddress( "Ds_DOCA" ,&Ds_DOCA);
tree -> SetBranchAddress( "Ds_DOCA2" ,&Ds_DOCA2);
tree -> SetBranchAddress( "Ds_DOCA3" ,&Ds_DOCA3);
tree -> SetBranchAddress( "pi_minus_fromDs_TRACK_CHI2NDOF" , &pi_minus_fromDs_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "K_plus_TRACK_CHI2NDOF" , &K_plus_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_minus_TRACK_CHI2NDOF" , &pi_minus_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_plus_TRACK_CHI2NDOF" , &pi_plus_TRACK_CHI2NDOF );

//PID variables
tree -> SetBranchAddress( "pi_minus_fromDs_PIDK" , &pi_minus_fromDs_PIDK );
tree -> SetBranchAddress( "K_plus_PIDK" , &K_plus_PIDK);
tree -> SetBranchAddress( "pi_minus_PIDK" , &pi_minus_PIDK);
tree -> SetBranchAddress( "pi_plus_PIDK" , &pi_plus_PIDK);



//define histograms--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TCanvas* c = new TCanvas();

//general histos
TH1D* mass_peak_Bs = new TH1D("B_{s} candidates", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
TH1D* mass_peak_Ds = new TH1D("D_{s} candidates", ";m(KK#pi) [MeV];Entries", 75, 1920., 2020.);
TH1D* ct_Bs = new TH1D("decay time of B_{s} candidates",";time [ps]; Entries", 75, 0., 10.);

//different intermediate states
//Ds->Phipi->KKpi
TH1D* mass_Bs_Phipi_Dsideband = new TH1D("D_{s}->#Phi#pi back", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_Phipi_Dsideband->SetFillColor(kBlue);
mass_Bs_Phipi_Dsideband->SetMarkerStyle(21);
mass_Bs_Phipi_Dsideband->SetMarkerColor(kBlue);
TH1D* mass_Bs_Phipi_Dsignal = new TH1D("D_{s}->'Phi#pi sig", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_Phipi_Dsignal->SetFillColor(kWhite);
mass_Bs_Phipi_Dsignal->SetMarkerStyle(21);
mass_Bs_Phipi_Dsignal->SetMarkerColor(kWhite);

//Ds->K*K->KKpi
TH1D* mass_Bs_KstarK_Dsideband = new TH1D("D_{s}->K^{*}K back ", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_KstarK_Dsideband->SetFillColor(kBlue);
mass_Bs_KstarK_Dsideband->SetMarkerStyle(21);
mass_Bs_KstarK_Dsideband->SetMarkerColor(kBlue);
TH1D* mass_Bs_KstarK_Dsignal = new TH1D("D_{s}->K^{*}K sig", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_KstarK_Dsignal->SetFillColor(kWhite);
mass_Bs_KstarK_Dsignal->SetMarkerStyle(21);
mass_Bs_KstarK_Dsignal->SetMarkerColor(kWhite);

//non resonant
TH1D* mass_Bs_nonRes_Dsideband = new TH1D("D_{s} nonRes back", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_nonRes_Dsideband->SetFillColor(kBlue);
mass_Bs_nonRes_Dsideband->SetMarkerStyle(21);
mass_Bs_nonRes_Dsideband->SetMarkerColor(kBlue);
TH1D* mass_Bs_nonRes_Dsignal = new TH1D("D_{s} nonRes sig", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_nonRes_Dsignal->SetFillColor(kWhite);
mass_Bs_nonRes_Dsignal->SetMarkerStyle(21);
mass_Bs_nonRes_Dsignal->SetMarkerColor(kWhite);

//all
TH1D* mass_Bs_all_Dsideband = new TH1D("B_{s} all candidates back", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_all_Dsideband->SetFillColor(kBlue);
mass_Bs_all_Dsideband->SetMarkerStyle(21);
mass_Bs_all_Dsideband->SetMarkerColor(kBlue);
TH1D* mass_Bs_all_Dsignal = new TH1D("B_{s} all candidates sig", ";m(D_{s}K#pi#pi) [MeV];Entries", 75, 5000., 5600.);
mass_Bs_all_Dsignal->SetFillColor(kWhite);
mass_Bs_all_Dsignal->SetMarkerStyle(21);
mass_Bs_all_Dsignal->SetMarkerColor(kWhite);


//the stacked plots
THStack *PhiPi = new THStack("PhiPi","D_{s}->#Phi#pi->KK#pi candidates");
THStack *KstK = new THStack("KstK","D_{s}->K*K->KK#pi candidates");
THStack *nonRes = new THStack("nonRes","non resonant candidates");
THStack *all = new THStack("all","all candidates");
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//Create output file
	TFile* output = 0;

	// Ds2KKpi case
	if(Ds2KKpi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_with_BDT_variables_S21_PID_temporary.root","RECREATE");
	else if(Ds2KKpi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_with_BDT_variables_S21_PID.root","RECREATE");
	else if(Ds2KKpi && MC && Ds3piBkg) output = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection/bkg_Ds3pi.root","RECREATE");
	else if(Ds2KKpi && MC && Dsstar3piBkg) output = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection/bkg_Dsstar3pi.root","RECREATE");
	else if(Ds2KKpi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_with_BDT_variables_Reco14.root","RECREATE");
	//else if(Ds2KKpi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_fullSelectionL0HLT1Triggered.root","RECREATE");
	else if(Ds2KKpi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_with_BDT_variables_Reco14.root","RECREATE");
	//else if(Ds2KKpi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_fullSelectionL0HLT1Triggered.root","RECREATE");

	//else if(Ds2KKpi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/HLT2efficiency/withHLT2.root","RECREATE");
	//else if(Ds2KKpi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/HLT2efficiency/withHLT2.root","RECREATE");


	//Ds2Kpipi case
	if(Ds2Kpipi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2Kpipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2Kpipi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2Kpipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2Kpipi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2Kpipi_with_BDT_variables.root","RECREATE");
	if(Ds2Kpipi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2Kpipi_with_BDT_variables.root","RECREATE");


	// Ds2pipipi case
	if(Ds2pipipi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2pipipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2pipipi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2pipipi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2pipipi_with_BDT_variables_Reco14.root","RECREATE");
	if(Ds2pipipi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2pipipi_with_BDT_variables_Reco14.root","RECREATE");


	TTree* summary_tree = tree->CloneTree(0);

	float DsDaughters_min_IPCHI2 = 0;
	float XsDaughters_min_IPCHI2 = 0;
	float DsDaughters_max_IPCHI2 = 0;
	float XsDaughters_max_IPCHI2 = 0;
	float DsDaughters_min_PT = 0;
	float XsDaughters_min_PT = 0;
	float Xs_max_DOCA = 0;
	float max_TrackChi2 = 0;
	float min_TrackChi2 = 0;
	float max_ghostProb = 0;
	int resonant = 0 ;
        double Bs_ct = 0;
        double Bs_cterr = 0;
        double Ds_ct = 0;
        double Ds_cterr = 0;
	int qt = 0;
	int qf = 0;

	summary_tree->Branch("DsDaughters_min_IPCHI2",&DsDaughters_min_IPCHI2,"DsDaughters_min_IPCHI2/F");
        summary_tree->Branch("XsDaughters_min_IPCHI2",&XsDaughters_min_IPCHI2,"XsDaughters_min_IPCHI2/F");
	summary_tree->Branch("DsDaughters_max_IPCHI2",&DsDaughters_max_IPCHI2,"DsDaughters_max_IPCHI2/F");
	summary_tree->Branch("XsDaughters_max_IPCHI2",&XsDaughters_max_IPCHI2,"XsDaughters_max_IPCHI2/F");
	summary_tree->Branch("DsDaughters_min_PT",&DsDaughters_min_PT,"DsDaughters_min_PT/F");
	summary_tree->Branch("XsDaughters_min_PT",&XsDaughters_min_PT,"XsDaughters_min_PT/F");
	summary_tree->Branch("Xs_max_DOCA",&Xs_max_DOCA,"Xs_max_DOCA/F");
	summary_tree->Branch("max_TrackChi2",&max_TrackChi2,"max_Track_Chi2/F");
	summary_tree->Branch("min_TrackChi2",&min_TrackChi2,"min_Track_Chi2/F");
	summary_tree->Branch("max_ghostProb",&max_ghostProb,"max_ghostProb/F");
	summary_tree->Branch("resonant",&resonant,"resonant/I");
        summary_tree->Branch("Bs_ct",&Bs_ct,"Bs_ct/D");
        summary_tree->Branch("Bs_cterr",&Bs_cterr,"Bs_cterr/D");
        summary_tree->Branch("Ds_ct",&Ds_ct,"Ds_ct/D");
        summary_tree->Branch("Ds_cterr",&Ds_cterr,"Ds_cterr/D");
	summary_tree->Branch("qt",&qt,"qt/I");
	summary_tree->Branch("qf",&qf,"qf/I");

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//define some intermediate floats to compute max/min track chi2
float interMin12 = 0;
float interMin34 = 0;
float interMin56 = 0;
float interMin1to4 = 0;

float interMax12 = 0;
float interMax34 = 0;
float interMax56 = 0;
float interMax1to4 = 0;


//loop over events
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
	tree->GetEntry(i);

	//define the Lorentz vectors
	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
	pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);

	if(Ds2KKpi){
		K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
		Kminus_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
	}
	if(Ds2Kpipi || Ds2pipipi){
		pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
	}
	if(Ds2KKpi || Ds2Kpipi){
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
		Kminus_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);
	}

	if(Ds2pipipi){
		pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);

	}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//trigger selection
	//L0 stage
	if((!Bs_L0Global_TIS) && (!Bs_L0HadronDecision_TOS)) continue;

	//HLT 1 stage
	if(!Bs_Hlt1TrackAllL0Decision_TOS) continue;

	//HLT2 stage
	if((!Bs_Hlt2Topo2BodyBBDTDecision_TOS) &&  (!Bs_Hlt2Topo3BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo4BodyBBDTDecision_TOS)) continue;
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	if(Ds2KKpi){
        	if(K_plus_fromDs_PT<100) continue;
		if(K_plus_fromDs_P<1600) continue;
		if(K_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_plus_fromDs_IPCHI2_OWNPV< 4) continue; 
		if((pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT) < 1800) continue;
		if(!MC)if(K_plus_fromDs_PIDK<-10) continue;
	}

	if(Ds2Kpipi){
		if((pi_minus_fromDs_PT + K_minus_fromDs_PT + pi_plus_fromDs_PT) < 1800) continue;
	}

	if(Ds2Kpipi || Ds2pipipi){
      		if(pi_plus_fromDs_PT<100) continue;
		if(pi_plus_fromDs_P<1600) continue;
		if(pi_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_plus_fromDs_IPCHI2_OWNPV< 4) continue;
		if(!MC)if(pi_plus_fromDs_PIDK>10) continue;
	}

	if(Ds2KKpi || Ds2Kpipi){
		if(K_minus_fromDs_PT<100) continue;
		if(K_minus_fromDs_P<1600) continue;
		if(K_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_minus_fromDs_IPCHI2_OWNPV< 4) continue;
		if(!MC)if(K_minus_fromDs_PIDK<-10) continue;
	}

	if(Ds2pipipi){
      		if(pi_minus2_fromDs_PT<100) continue;
		if(pi_minus2_fromDs_P<1600) continue;
		if(pi_minus2_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_minus2_fromDs_IPCHI2_OWNPV< 4) continue;
		if(!MC)if(pi_minus2_fromDs_PIDK>10) continue;
		if((pi_minus_fromDs_PT + pi_minus2_fromDs_PT + pi_plus_fromDs_PT) < 1800) continue;
	}

	//track pt cut
 	if(pi_minus_fromDs_PT<100) continue;
	if(K_plus_PT<100) continue;
	if(pi_plus_PT<100) continue;
	if(pi_minus_PT<100) continue;

	//track p cut
	if(pi_minus_fromDs_P<1600) continue;
	if(K_plus_P<1000) continue;
	if(pi_plus_P<1000) continue;
	if(pi_minus_P<1000) continue;

	//track chi2 cut
	if(pi_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
	if(K_plus_TRACK_CHI2NDOF> 4) continue;
	if(pi_plus_TRACK_CHI2NDOF> 4) continue;
	if(pi_minus_TRACK_CHI2NDOF> 4) continue;

	//track ip chi2 cut
	if(pi_minus_fromDs_IPCHI2_OWNPV< 4) continue;
	if(K_plus_IPCHI2_OWNPV< 4) continue;
	if(pi_plus_IPCHI2_OWNPV< 4) continue;
	if(pi_minus_IPCHI2_OWNPV< 4) continue;
	

	//ds daughters pt cut

	//ds daughters doca cut
	if((Ds2Kpipi && sevenTeV) || Ds2KKpi || (Ds2pipipi && sevenTeV)) if(Ds_DOCA1> 0.5) continue;
	if((Ds2pipipi && (!sevenTeV)) || (Ds2Kpipi && (!sevenTeV))) if(Ds_DOCA> 0.5) continue;
	if(Ds_DOCA2> 0.5) continue;
	if(Ds_DOCA3> 0.5) continue;

	//ds mass window of 40 MeV
	if( Ds_MM < 1950 || Ds_MM > 1990) continue;

	//ds vertex chi2 cut
	if((Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF) > 10) continue;

	//ds fd chi2 cut
	if(Ds_FDCHI2_OWNPV<36) continue;

	//loose pid requirements
if(!MC){
	if(K_plus_PIDK<8) continue;
	if(pi_minus_fromDs_PIDK>10) continue;
	if(pi_plus_PIDK>10) continue;
	if(pi_minus_PIDK>10) continue;
}
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//cuts on Xs candidates
	
	//Xs DOCA cut
	if((Ds2Kpipi && sevenTeV) || Ds2KKpi || (Ds2pipipi && sevenTeV)) if(K_1_1270_plus_DOCA1> 0.4) continue;
	if((Ds2pipipi && (!sevenTeV)) || (Ds2Kpipi && (!sevenTeV))) if(K_1_1270_plus_DOCA> 0.4) continue;
	if(K_1_1270_plus_DOCA2> 0.4) continue;
	if(K_1_1270_plus_DOCA3> 0.4) continue;

	//Xs pt sum cut
	if((pi_minus_PT + K_plus_PT + pi_plus_PT) < 1250) continue;

	//Xs vertex chi2 cut
	if((K_1_1270_plus_ENDVERTEX_CHI2/K_1_1270_plus_ENDVERTEX_NDOF) > 8) continue;

	//Xs fd chi2
	if(K_1_1270_plus_FDCHI2_OWNPV < 16) continue;

	//vertex displacement 
	 if(TMath::Sqrt(((K_1_1270_plus_ENDVERTEX_X - K_1_1270_plus_OWNPV_X)*(K_1_1270_plus_ENDVERTEX_X - K_1_1270_plus_OWNPV_X))+((K_1_1270_plus_ENDVERTEX_Y - K_1_1270_plus_OWNPV_Y)*(K_1_1270_plus_ENDVERTEX_Y - K_1_1270_plus_OWNPV_Y)))<0.1) continue;

	if(TMath::Abs(K_1_1270_plus_ENDVERTEX_Z - K_1_1270_plus_OWNPV_Z)<2.0) continue;

	//Xs DIRA
	if(K_1_1270_plus_DIRA_OWNPV < 0.98) continue;


	//selection cuts specific to Ds->Kpipi--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if(Ds2Kpipi){
		//D^0 veto
		if((pi_plus_fromDs + K_minus_fromDs).M() > 1750) continue;

		//Lambda_c veto
		if((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 0 && TMath::Abs((pi_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 30) continue;

		//PID requirements
		if(K_minus_fromDs_PIDK < 10) continue;
		if(pi_plus_fromDs_PIDK > 5 || pi_minus_fromDs_PIDK > 5) continue;
	}

	//selection cuts specific to Ds->pipipi------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	if(Ds2pipipi){
		//D^0 veto
		if((pi_plus_fromDs + pi_minus_fromDs).M() > 1700 || (pi_plus_fromDs + pi_minus2_fromDs).M() > 1700) continue;

		if(Ds_FDCHI2_ORIVX < 9) continue;

		//PID requirements
		if(!MC) if(pi_plus_fromDs_PIDK > 10 || pi_minus_fromDs_PIDK > 10 || pi_minus2_fromDs_PIDK > 10) continue;
		if(!MC) if(pi_plus_fromDs_PIDp > 10 || pi_minus_fromDs_PIDp > 10 || pi_minus2_fromDs_PIDp > 10) continue;
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//additional loose requirements on b-hadron

	//dira cut
	if(Bs_DIRA_OWNPV<0.99994) continue;

	//ip chi2 cut
	if(Bs_IPCHI2_OWNPV>20) continue;

	//fd chi2 cut
	if(Bs_FDCHI2_OWNPV<100) continue;

	//vertex fit quality 
	if((Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)> 8) continue;

	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(!MC){
	//implement PID requirements on Ds daughters
	if(Ds2KKpi){
		resonant = 0;
		if( TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) < 20) resonant = 1; 
		else if(TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) < 75) {
			if((K_plus_fromDs_PIDK>0 && K_minus_fromDs_PIDK>0) ) resonant = 1;
			else continue;
		}  		
		else{
			if((K_plus_fromDs_PIDK>5 && K_minus_fromDs_PIDK>5) ) resonant = 0;
			else continue;
		}
	}
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
}
	//rejection cuts for peaking background
	if(Ds2KKpi){
		if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) continue;
		if(Ds_FDCHI2_ORIVX < 2) continue;
	}


	//Bs->DsDs suppression
	if(TMath::Abs((K_plus + pi_plus + pi_minus).M() - massDs) < 20) continue;

	if(!MC){

	//Bs->DsKKpi suppression
		if(pi_minus_PIDK > 5) continue;
	}

	if(Ds2KKpi){
		//if(TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20){

		//B0->D^-Kpipi suppression
		if( TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 30. && K_minus_fromDs_PIDK < 10 ) continue;

		//D⁰ veto
		if((K_plus_fromDs + K_minus_fromDs).M() > 1840.) continue;

		//Lambda_b->Lambda_c Kpipi suppression
		if( TMath::Abs((K_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 30. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) continue;

		//}
	}
	//MC Truth matching

	if(MC && Ds2KKpi){
		//if(Bs_BKGCAT == 50) continue; // remove ghosts
		if(TMath::Abs(Bs_TRUEID) != 531) continue;
		if(TMath::Abs(Ds_TRUEID) != 431) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(K_plus_TRUEID) != 321) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(pi_plus_TRUEID) != 211) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(pi_minus_TRUEID) != 211) continue;
		if(TMath::Abs(K_plus_fromDs_TRUEID) != 321) continue;
		if(TMath::Abs(pi_minus_fromDs_TRUEID) != 211) continue;
		if(TMath::Abs(K_minus_fromDs_TRUEID) != 321) continue;
	}

	if(MC && Ds2pipipi){
		//if(Bs_BKGCAT == 50) continue; // remove ghosts
		if(TMath::Abs(Bs_TRUEID) != 531) continue;
		if(TMath::Abs(Ds_TRUEID) != 431) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(K_plus_TRUEID) != 321) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(pi_plus_TRUEID) != 211) continue;
		if(!Ds3piBkg)if(!Dsstar3piBkg)if(TMath::Abs(pi_minus_TRUEID) != 211) continue;
		if(TMath::Abs(pi_plus_fromDs_TRUEID) != 211) continue;
		if(TMath::Abs(pi_minus_fromDs_TRUEID) != 211) continue;
		if(TMath::Abs(pi_minus2_fromDs_TRUEID) != 211) continue;
	}

	//fill all sorts of histrograms
	ct_Bs->Fill(Bs_TAU*1000);
	mass_peak_Bs->Fill(Bs_MM);
	mass_peak_Ds->Fill(Ds_MM);
        Bs_ct = Bs_TAU*1000;
        Bs_cterr = Bs_TAUERR*1000;
        Ds_ct = Ds_TAU*1000;
        Ds_cterr = Ds_TAUERR*1000;

	if(Bs_ID > 0) qt = 1;
	if(Bs_ID < 0) qt = -1;
	if(Bs_ID == 0) qt = 0;
	if(K_plus_ID > 0) qf = 1;
	if(K_plus_ID < 0) qf = -1;


	if(Ds2KKpi){
	//Ds->PhiPi->KKpi case
		if(TMath::Abs((K_plus_fromDs + K_minus_fromDs).M() - massPhi) < 20 ){
			if((Ds_MM < (massDs - 25)) || (Ds_MM > (massDs + 45))) mass_Bs_Phipi_Dsideband->Fill(Bs_MM);
			if( TMath::Abs(Ds_MM - massDs) < 20 ) mass_Bs_Phipi_Dsignal->Fill(Bs_MM);
		}

	//Ds->K*K
		if(TMath::Abs((pi_minus_fromDs + K_plus_fromDs).M() - massKstar) < 75){
			if((Ds_MM < (massDs - 25)) || (Ds_MM > (massDs + 45))) mass_Bs_KstarK_Dsideband->Fill(Bs_MM);
			if( TMath::Abs(Ds_MM - massDs) < 20 ) mass_Bs_KstarK_Dsignal->Fill(Bs_MM);
		}

	//non res-case
		if( TMath::Abs((K_plus_fromDs + K_minus_fromDs).M() - massPhi) > 20 && TMath::Abs((pi_minus_fromDs + K_plus_fromDs).M() - massKstar) > 75 ){
			if((Ds_MM < (massDs - 25)) || (Ds_MM > (massDs + 45))) mass_Bs_nonRes_Dsideband->Fill(Bs_MM);
			if( TMath::Abs(Ds_MM - massDs) < 20 ) mass_Bs_nonRes_Dsignal->Fill(Bs_MM);
		}

	//all
	 	if((Ds_MM < (massDs - 25)) || (Ds_MM > (massDs + 45))) mass_Bs_all_Dsideband->Fill(Bs_MM);
	 	if( TMath::Abs(Ds_MM - massDs) < 20 ) mass_Bs_all_Dsignal->Fill(Bs_MM);
	}
	//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	//compute the variables used for BDT training ------> only for Ds->KKpi
	if(Ds2KKpi){
		//min IP chi² of Ds daughters
		if( (K_plus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
		if( (K_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
		if( (pi_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

		//min IP chi² of Xs daughters
		if( (K_plus_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) && (K_plus_IPCHI2_OWNPV < pi_plus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = K_plus_IPCHI2_OWNPV;
		if( (pi_minus_IPCHI2_OWNPV < K_plus_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV < pi_plus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_minus_IPCHI2_OWNPV;
		if( (pi_plus_IPCHI2_OWNPV < K_plus_IPCHI2_OWNPV) && (pi_plus_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus_IPCHI2_OWNPV;

		//max IP chi² of Xs daughters
		if( (K_plus_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) && (K_plus_IPCHI2_OWNPV > pi_plus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = K_plus_IPCHI2_OWNPV;
		if( (pi_minus_IPCHI2_OWNPV > K_plus_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV > pi_plus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_minus_IPCHI2_OWNPV;
		if( (pi_plus_IPCHI2_OWNPV > K_plus_IPCHI2_OWNPV) && (pi_plus_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus_IPCHI2_OWNPV;

		//max IP chi² of Ds daughters
		if( (K_plus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
		if( (K_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
		if( (pi_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;		

		//min p_t of Ds daughters
		if( (K_plus_fromDs_PT < K_minus_fromDs_PT) && (K_plus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_plus_fromDs_PT;
		if( (K_minus_fromDs_PT < K_plus_fromDs_PT) && (K_minus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_minus_fromDs_PT;
		if( (pi_minus_fromDs_PT < K_plus_fromDs_PT) && (pi_minus_fromDs_PT < K_minus_fromDs_PT) ) DsDaughters_min_PT = pi_minus_fromDs_PT;

		//min p_t of Xs daughters
		if( (K_plus_PT < pi_minus_PT) && (K_plus_PT < pi_plus_PT) ) XsDaughters_min_PT = K_plus_PT;
		if( (pi_minus_PT < K_plus_PT) && (pi_minus_PT < pi_plus_PT) ) XsDaughters_min_PT = pi_minus_PT;
		if( (pi_plus_PT < K_plus_PT) && (pi_plus_PT < pi_minus_PT) ) XsDaughters_min_PT = pi_plus_PT;


		//max DOCA of Xs
		if( (K_1_1270_plus_DOCA3 > K_1_1270_plus_DOCA2) && (K_1_1270_plus_DOCA3 > K_1_1270_plus_DOCA1) ) Xs_max_DOCA =  K_1_1270_plus_DOCA3;
		if( (K_1_1270_plus_DOCA2 > K_1_1270_plus_DOCA3) && (K_1_1270_plus_DOCA2 > K_1_1270_plus_DOCA1) ) Xs_max_DOCA =  K_1_1270_plus_DOCA2;
		if( (K_1_1270_plus_DOCA1 > K_1_1270_plus_DOCA2) && (K_1_1270_plus_DOCA1 > K_1_1270_plus_DOCA3) ) Xs_max_DOCA =  K_1_1270_plus_DOCA1;

		//min Track chi2
		interMin12 = TMath::Min(K_plus_fromDs_TRACK_CHI2NDOF,K_minus_fromDs_TRACK_CHI2NDOF);
		interMin34 = TMath::Min(pi_minus_fromDs_TRACK_CHI2NDOF,K_plus_TRACK_CHI2NDOF);
		interMin56 = TMath::Min(pi_minus_TRACK_CHI2NDOF,pi_plus_TRACK_CHI2NDOF);
		interMin1to4 = TMath::Min(interMin12,interMin34);
		min_TrackChi2 = TMath::Min(interMin1to4,interMin56);

		//max Track chi2
		interMax12 = TMath::Max(K_plus_fromDs_TRACK_CHI2NDOF,K_minus_fromDs_TRACK_CHI2NDOF);
		interMax34 = TMath::Max(pi_minus_fromDs_TRACK_CHI2NDOF,K_plus_TRACK_CHI2NDOF);
		interMax56 = TMath::Max(pi_minus_TRACK_CHI2NDOF,pi_plus_TRACK_CHI2NDOF);
		interMax1to4 = TMath::Max(interMax12,interMax34);
		max_TrackChi2 = TMath::Max(interMax1to4,interMax56);

		//max GhostProb
		max_ghostProb = TMath::Max(pi_plus_TRACK_GhostProb,TMath::Max(K_plus_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(K_plus_fromDs_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb)))));
	}
	//compute the variables used for BDT training ------> only for Ds->pipipi
	if(Ds2pipipi){
		//min IP chi² of Ds daughters
		if( (pi_plus_fromDs_IPCHI2_OWNPV < pi_minus2_fromDs_IPCHI2_OWNPV) && (pi_plus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_plus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus2_fromDs_IPCHI2_OWNPV < pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus2_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus2_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV < pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV < pi_minus2_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

		//min IP chi² of Xs daughters
		if( (K_plus_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) && (K_plus_IPCHI2_OWNPV < pi_plus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = K_plus_IPCHI2_OWNPV;
		if( (pi_minus_IPCHI2_OWNPV < K_plus_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV < pi_plus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_minus_IPCHI2_OWNPV;
		if( (pi_plus_IPCHI2_OWNPV < K_plus_IPCHI2_OWNPV) && (pi_plus_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus_IPCHI2_OWNPV;

		//max IP chi² of Xs daughters
		if( (K_plus_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) && (K_plus_IPCHI2_OWNPV > pi_plus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = K_plus_IPCHI2_OWNPV;
		if( (pi_minus_IPCHI2_OWNPV > K_plus_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV > pi_plus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_minus_IPCHI2_OWNPV;
		if( (pi_plus_IPCHI2_OWNPV > K_plus_IPCHI2_OWNPV) && (pi_plus_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus_IPCHI2_OWNPV;

		//max IP chi² of Ds daughters
                if( (pi_plus_fromDs_IPCHI2_OWNPV > pi_minus2_fromDs_IPCHI2_OWNPV) && (pi_plus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_plus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus2_fromDs_IPCHI2_OWNPV > pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus2_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus2_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV > pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV > pi_minus2_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

		//min p_t of Ds daughters
                if( (pi_plus_fromDs_PT < pi_minus2_fromDs_PT) && (pi_plus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = pi_plus_fromDs_PT;
                if( (pi_minus2_fromDs_PT < pi_plus_fromDs_PT) && (pi_minus2_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = pi_minus2_fromDs_PT;
                if( (pi_minus_fromDs_PT < pi_plus_fromDs_PT) && (pi_minus_fromDs_PT < pi_minus2_fromDs_PT) ) DsDaughters_min_PT = pi_minus_fromDs_PT;

		//min p_t of Xs daughters
		if( (K_plus_PT < pi_minus_PT) && (K_plus_PT < pi_plus_PT) ) XsDaughters_min_PT = K_plus_PT;
		if( (pi_minus_PT < K_plus_PT) && (pi_minus_PT < pi_plus_PT) ) XsDaughters_min_PT = pi_minus_PT;
		if( (pi_plus_PT < K_plus_PT) && (pi_plus_PT < pi_minus_PT) ) XsDaughters_min_PT = pi_plus_PT;

		//max DOCA of Xs
		if( (K_1_1270_plus_DOCA3 > K_1_1270_plus_DOCA2) && (K_1_1270_plus_DOCA3 > K_1_1270_plus_DOCA1) ) Xs_max_DOCA =  K_1_1270_plus_DOCA3;
		if( (K_1_1270_plus_DOCA2 > K_1_1270_plus_DOCA3) && (K_1_1270_plus_DOCA2 > K_1_1270_plus_DOCA1) ) Xs_max_DOCA =  K_1_1270_plus_DOCA2;
		if( (K_1_1270_plus_DOCA1 > K_1_1270_plus_DOCA2) && (K_1_1270_plus_DOCA1 > K_1_1270_plus_DOCA3) ) Xs_max_DOCA =  K_1_1270_plus_DOCA1;

		//min Track chi2
		interMin12 = TMath::Min(pi_plus_fromDs_TRACK_CHI2NDOF,pi_minus2_fromDs_TRACK_CHI2NDOF);
		interMin34 = TMath::Min(pi_minus_fromDs_TRACK_CHI2NDOF,K_plus_TRACK_CHI2NDOF);
		interMin56 = TMath::Min(pi_minus_TRACK_CHI2NDOF,pi_plus_TRACK_CHI2NDOF);
		interMin1to4 = TMath::Min(interMin12,interMin34);
		min_TrackChi2 = TMath::Min(interMin1to4,interMin56);

		//max Track chi2
		interMax12 = TMath::Max(pi_plus_fromDs_TRACK_CHI2NDOF,pi_minus2_fromDs_TRACK_CHI2NDOF);
		interMax34 = TMath::Max(pi_minus_fromDs_TRACK_CHI2NDOF,K_plus_TRACK_CHI2NDOF);
		interMax56 = TMath::Max(pi_minus_TRACK_CHI2NDOF,pi_plus_TRACK_CHI2NDOF);
		interMax1to4 = TMath::Max(interMax12,interMax34);
		max_TrackChi2 = TMath::Max(interMax1to4,interMax56);

		//max GhostProb
		max_ghostProb = TMath::Max(pi_plus_TRACK_GhostProb,TMath::Max(K_plus_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(pi_plus_fromDs_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,pi_minus2_fromDs_TRACK_GhostProb)))));
	}
		summary_tree->Fill();

}

//draw histograms for different years and final states

string BsMassDistribution = "eps/mass_Bs_preSel";
string DsMassDistribution = "eps/mass_Ds_preSel";
string BsLifetimeDistribution = "eps/decaytime";
string PhiPiDistribution = "eps/mass_contributions_PhiPi";
string KstKDistribution = "eps/mass_contributions_KstK";
string NonResDistribution = "eps/mass_contributions_nonRes";
string AllDistribution = "eps/mass_contributions_all";

//which final state? (default is Ds->KKpi)
string Kpipi = "_Ds2Kpipi";
string pipipi = "_Ds2pipipi";

//which year? data or mc?
string elevenData = "_11Data.eps";
string twelveData = "_12Data.eps";
string elevenMC = "_11MC.eps";
string twelveMC = "_12MC.eps";

if(Ds2KKpi){
	if((!MC) && (!sevenTeV)){
		AllDistribution.append(twelveData);
		NonResDistribution.append(twelveData);
		KstKDistribution.append(twelveData);
		PhiPiDistribution.append(twelveData);
		BsMassDistribution.append(twelveData);
		DsMassDistribution.append(twelveData);
		BsLifetimeDistribution.append(twelveData);
	}
	if((!MC) && (sevenTeV)){
		AllDistribution.append(elevenData);
		NonResDistribution.append(elevenData);
		KstKDistribution.append(elevenData);
		PhiPiDistribution.append(elevenData);
		BsMassDistribution.append(elevenData);
		DsMassDistribution.append(elevenData);
		BsLifetimeDistribution.append(elevenData);
	}
	if((MC) && (!sevenTeV)){
		AllDistribution.append(twelveMC);
		NonResDistribution.append(twelveMC);
		KstKDistribution.append(twelveMC);
		PhiPiDistribution.append(twelveMC);
		BsMassDistribution.append(twelveMC);
		DsMassDistribution.append(twelveMC);
		BsLifetimeDistribution.append(twelveMC);
	}
	if((MC) && (sevenTeV)){
		AllDistribution.append(elevenMC);
		NonResDistribution.append(elevenMC);
		KstKDistribution.append(elevenMC);
		PhiPiDistribution.append(elevenMC);
		BsMassDistribution.append(elevenMC);
		DsMassDistribution.append(elevenMC);
		BsLifetimeDistribution.append(elevenMC);
	}
}
if(Ds2Kpipi){
	BsMassDistribution.append(Kpipi);
	DsMassDistribution.append(Kpipi);
	BsLifetimeDistribution.append(Kpipi);

	if((!MC) && (!sevenTeV)){
		BsMassDistribution.append(twelveData);
		DsMassDistribution.append(twelveData);
		BsLifetimeDistribution.append(twelveData);
	}
	if((!MC) && (sevenTeV)){
		BsMassDistribution.append(elevenData);
		DsMassDistribution.append(elevenData);
		BsLifetimeDistribution.append(elevenData);
	}
	if((MC) && (!sevenTeV)){
		BsMassDistribution.append(twelveMC);
		DsMassDistribution.append(twelveMC);
		BsLifetimeDistribution.append(twelveMC);
	}
	if((MC) && (sevenTeV)){ 
		BsMassDistribution.append(elevenMC);
		DsMassDistribution.append(elevenMC);
		BsLifetimeDistribution.append(elevenMC);
	}
}

if(Ds2pipipi){
	BsMassDistribution.append(pipipi);
	DsMassDistribution.append(pipipi);
	BsLifetimeDistribution.append(pipipi);

	if((!MC) && (!sevenTeV)){
		BsMassDistribution.append(twelveData);
		DsMassDistribution.append(twelveData);
		BsLifetimeDistribution.append(twelveData);
	}
	if((!MC) && (sevenTeV)){
		BsMassDistribution.append(elevenData);
		DsMassDistribution.append(elevenData);
		BsLifetimeDistribution.append(elevenData);
	}
	if((MC) && (!sevenTeV)){
		BsMassDistribution.append(twelveMC);
		DsMassDistribution.append(twelveMC);
		BsLifetimeDistribution.append(twelveMC);
	}
	if((MC) && (sevenTeV)){ 
		BsMassDistribution.append(elevenMC);
		DsMassDistribution.append(elevenMC);
		BsLifetimeDistribution.append(elevenMC);
	}
}

//convert strings to chars in order to use in c->Print();
const char* cstringBMass = BsMassDistribution.c_str();
const char* cstringDMass = DsMassDistribution.c_str();
const char* cstringBLifetime = BsLifetimeDistribution.c_str();
const char* cstringPhiPiDistribution = PhiPiDistribution.c_str();
const char* cstringKstKDistribution = KstKDistribution.c_str();
const char* cstringNonResDistribution = NonResDistribution.c_str();
const char* cstringAllDistribution = AllDistribution.c_str();

mass_peak_Bs->Draw("E1"); c->Print(cstringBMass);
mass_peak_Ds->Draw("E1"); c->Print(cstringDMass);
ct_Bs->Draw("E1"); c->Print(cstringBLifetime);

if(Ds2KKpi){
	PhiPi->Add(mass_Bs_Phipi_Dsideband); PhiPi->Add(mass_Bs_Phipi_Dsignal);
	PhiPi->Draw();
	PhiPi->GetXaxis()->SetTitle("m(D_{s}K#pi#pi) [MeV]");
	PhiPi->GetYaxis()->SetTitle("Candidates  ");
	c->Print(cstringPhiPiDistribution);

	KstK->Add(mass_Bs_KstarK_Dsideband); KstK->Add(mass_Bs_KstarK_Dsignal);
	KstK->Draw();
	KstK->GetXaxis()->SetTitle("m(D_{s}K#pi#pi) [MeV]");
	KstK->GetYaxis()->SetTitle("Candidates ");
	c->Print(cstringKstKDistribution);

	nonRes->Add(mass_Bs_nonRes_Dsideband); nonRes->Add(mass_Bs_nonRes_Dsignal);
	nonRes->Draw();
	nonRes->GetXaxis()->SetTitle("m(D_{s}K#pi#pi) [MeV]");
	nonRes->GetYaxis()->SetTitle("Candidates  ");
	c->Print(cstringNonResDistribution);

	all->Add(mass_Bs_all_Dsideband); all->Add(mass_Bs_all_Dsignal);
	all->Draw();
	all->GetXaxis()->SetTitle("m(D_{s}K#pi#pi) [MeV]");
	all->GetYaxis()->SetTitle("Candidates   ");
	c->Print(cstringAllDistribution);
}

cout << "New file contains " << summary_tree->GetEntries() << " events" <<  endl;
summary_tree->Write();
output->Close();	

return 0;
}

void chooseBestPV(string fit, string input, string output){
 
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
    float Ds_Kp_PX[100], Ds_Kp_PY[100], Ds_Kp_PZ[100], Ds_Kp_PE[100];
    float Ds_pip_PX[100], Ds_pip_PY[100], Ds_pip_PZ[100], Ds_pip_PE[100];
    float Ds_Km_PX[100], Ds_Km_PY[100], Ds_Km_PZ[100], Ds_Km_PE[100];

    tree->SetBranchAddress(("Bs_"+fit+"_nPV").c_str(),&nPV);
    tree->SetBranchAddress(("Bs_"+fit+"_chi2").c_str(),&chi2);
    tree->SetBranchAddress(("Bs_"+fit+"_nDOF").c_str(),&ndof);

    tree->SetBranchAddress(("Bs_"+fit+"_M").c_str(),&B_M);
    tree->SetBranchAddress(("Bs_"+fit+"_ctau").c_str(),&B_TAU);

    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_Kplus_PX").c_str(),&K_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_Kplus_PY").c_str(),&K_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_Kplus_PZ").c_str(),&K_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_Kplus_PE").c_str(),&K_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_PX").c_str(),&pip_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_PY").c_str(),&pip_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_PZ").c_str(),&pip_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_PE").c_str(),&pip_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_0_PX").c_str(),&pim_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_0_PY").c_str(),&pim_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_0_PZ").c_str(),&pim_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_K_1_1270_plus_piplus_0_PE").c_str(),&pim_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_PX").c_str(),&Ds_Kp_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_PY").c_str(),&Ds_Kp_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_PZ").c_str(),&Ds_Kp_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_PE").c_str(),&Ds_Kp_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_piplus_PX").c_str(),&Ds_pip_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_piplus_PY").c_str(),&Ds_pip_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_piplus_PZ").c_str(),&Ds_pip_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_piplus_PE").c_str(),&Ds_pip_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_0_PX").c_str(),&Ds_Km_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_0_PY").c_str(),&Ds_Km_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_0_PZ").c_str(),&Ds_Km_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_D_splus_Kplus_0_PE").c_str(),&Ds_Km_PE); 

    ///create new tree
    TFile* f = new TFile(output.c_str(),"RECREATE");
    TTree* new_tree = tree->CloneTree();;    
    double DTF_Bs_M, DTF_Bs_MM, DTF_Ds_M, DTF_chi2, DTF_Bs_TAU; 
    double DTF_Bs_PX, DTF_Bs_PY, DTF_Bs_PZ, DTF_Bs_PE;      
    double DTF_Kplus_PX, DTF_Kplus_PY, DTF_Kplus_PZ, DTF_Kplus_PE;      
    double DTF_piplus_PX, DTF_piplus_PY, DTF_piplus_PZ, DTF_piplus_PE;      
    double DTF_piminus_PX, DTF_piminus_PY, DTF_piminus_PZ, DTF_piminus_PE;

    double DTF_Ds_Kplus_PX, DTF_Ds_Kplus_PY, DTF_Ds_Kplus_PZ, DTF_Ds_Kplus_PE;      
    double DTF_Ds_piplus_PX, DTF_Ds_piplus_PY, DTF_Ds_piplus_PZ, DTF_Ds_piplus_PE;      
    double DTF_Ds_Kminus_PX, DTF_Ds_Kminus_PY, DTF_Ds_Kminus_PZ, DTF_Ds_Kminus_PE; 
    double DTF_Ds_PX, DTF_Ds_PY, DTF_Ds_PZ, DTF_Ds_PE;      

    TBranch* Bra_DTF_Bs_M = new_tree->Branch((fit+"_Bs_M").c_str(), &DTF_Bs_M, (fit+"_Bs_M/D").c_str());
    TBranch* Bra_DTF_Bs_MM = new_tree->Branch((fit+"_Bs_MM").c_str(), &DTF_Bs_MM, (fit+"_Bs_MM/D").c_str());
    TBranch* Bra_DTF_Bs_TAU = new_tree->Branch((fit+"_TAU").c_str(), &DTF_Bs_TAU, (fit+"_TAU/D").c_str());
    TBranch* Bra_DTF_Ds_M = new_tree->Branch((fit+"_Ds_M").c_str(), &DTF_Ds_M,  (fit+"_Ds_M/D").c_str());
    TBranch* Bra_DTF_chi2 = new_tree->Branch((fit+"_chi2").c_str(), &DTF_chi2,  (fit+"_chi2/D").c_str());

    TBranch* Bra_DTF_Bs_PX = new_tree->Branch((fit+"_Bs_PX").c_str(), &DTF_Bs_PX, (fit+"_Bs_PX/D").c_str());
    TBranch* Bra_DTF_Bs_PY = new_tree->Branch((fit+"_Bs_PY").c_str(), &DTF_Bs_PY, (fit+"_Bs_PY/D").c_str());
    TBranch* Bra_DTF_Bs_PZ = new_tree->Branch((fit+"_Bs_PZ").c_str(), &DTF_Bs_PZ, (fit+"_Bs_PZ/D").c_str());
    TBranch* Bra_DTF_Bs_PE = new_tree->Branch((fit+"_Bs_PE").c_str(), &DTF_Bs_PE, (fit+"_Bs_PE/D").c_str());

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

    TBranch* Bra_DTF_Ds_Kplus_PX = new_tree->Branch((fit+"_Ds_Kplus_PX").c_str(), &DTF_Ds_Kplus_PX, (fit+"_Ds_Kplus_PX/D").c_str());
    TBranch* Bra_DTF_Ds_Kplus_PY = new_tree->Branch((fit+"_Ds_Kplus_PY").c_str(), &DTF_Ds_Kplus_PY, (fit+"_Ds_Kplus_PY/D").c_str());
    TBranch* Bra_DTF_Ds_Kplus_PZ = new_tree->Branch((fit+"_Ds_Kplus_PZ").c_str(), &DTF_Ds_Kplus_PZ, (fit+"_Ds_Kplus_PZ/D").c_str());
    TBranch* Bra_DTF_Ds_Kplus_PE = new_tree->Branch((fit+"_Ds_Kplus_PE").c_str(), &DTF_Ds_Kplus_PE, (fit+"_Ds_Kplus_PE/D").c_str());

    TBranch* Bra_DTF_Ds_piplus_PX = new_tree->Branch((fit+"_Ds_piplus_PX").c_str(), &DTF_Ds_piplus_PX, (fit+"_Ds_piplus_PX/D").c_str());
    TBranch* Bra_DTF_Ds_piplus_PY = new_tree->Branch((fit+"_Ds_piplus_PY").c_str(), &DTF_Ds_piplus_PY, (fit+"_Ds_piplus_PY/D").c_str());
    TBranch* Bra_DTF_Ds_piplus_PZ = new_tree->Branch((fit+"_Ds_piplus_PZ").c_str(), &DTF_Ds_piplus_PZ, (fit+"_Ds_piplus_PZ/D").c_str());
    TBranch* Bra_DTF_Ds_piplus_PE = new_tree->Branch((fit+"_Ds_piplus_PE").c_str(), &DTF_Ds_piplus_PE, (fit+"_Ds_piplus_PE/D").c_str());

    TBranch* Bra_DTF_Ds_Kminus_PX = new_tree->Branch((fit+"_Ds_Kminus_PX").c_str(), &DTF_Ds_Kminus_PX, (fit+"_Ds_Kminus_PX/D").c_str());
    TBranch* Bra_DTF_Ds_Kminus_PY = new_tree->Branch((fit+"_Ds_Kminus_PY").c_str(), &DTF_Ds_Kminus_PY, (fit+"_Ds_Kminus_PY/D").c_str());
    TBranch* Bra_DTF_Ds_Kminus_PZ = new_tree->Branch((fit+"_Ds_Kminus_PZ").c_str(), &DTF_Ds_Kminus_PZ, (fit+"_Ds_Kminus_PZ/D").c_str());
    TBranch* Bra_DTF_Ds_Kminus_PE = new_tree->Branch((fit+"_Ds_Kminus_PE").c_str(), &DTF_Ds_Kminus_PE, (fit+"_Ds_Kminus_PE/D").c_str());

    TBranch* Bra_DTF_Ds_PX = new_tree->Branch((fit+"_Ds_PX").c_str(), &DTF_Ds_PX, (fit+"_Ds_PX/D").c_str());
    TBranch* Bra_DTF_Ds_PY = new_tree->Branch((fit+"_Ds_PY").c_str(), &DTF_Ds_PY, (fit+"_Ds_PY/D").c_str());
    TBranch* Bra_DTF_Ds_PZ = new_tree->Branch((fit+"_Ds_PZ").c_str(), &DTF_Ds_PZ, (fit+"_Ds_PZ/D").c_str());
    TBranch* Bra_DTF_Ds_PE = new_tree->Branch((fit+"_Ds_PE").c_str(), &DTF_Ds_PE, (fit+"_Ds_PE/D").c_str());

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

        TLorentzVector Ds_Kp_P(Ds_Kp_PX[best],Ds_Kp_PY[best],Ds_Kp_PZ[best],Ds_Kp_PE[best]);
        TLorentzVector Ds_pip_P(Ds_pip_PX[best],Ds_pip_PY[best],Ds_pip_PZ[best],Ds_pip_PE[best]);
        TLorentzVector Ds_Km_P(Ds_Km_PX[best],Ds_Km_PY[best],Ds_Km_PZ[best],Ds_Km_PE[best]);

	TLorentzVector Ds_P = Ds_Kp_P + Ds_Km_P + Ds_pip_P;
	TLorentzVector B_P = K_P + pip_P + pim_P + Ds_P;


	///Fill branches
	DTF_Bs_M = B_M[best];
	DTF_Bs_MM = B_P.M();
	//if(DTF_Bs_M<4000. || DTF_Bs_M>6000.)continue;
	DTF_Bs_TAU = B_TAU[best];
	DTF_chi2 = tmp;
	DTF_Ds_M = Ds_P.M();

	DTF_Bs_PX = B_P.X();
	DTF_Bs_PY = B_P.Y();
	DTF_Bs_PZ = B_P.Z();
	DTF_Bs_PE = B_P.E();

	DTF_Ds_PX = Ds_P.X();
	DTF_Ds_PY = Ds_P.Y();
	DTF_Ds_PZ = Ds_P.Z();
	DTF_Ds_PE = Ds_P.E();

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

	DTF_Ds_Kplus_PX = Ds_Kp_PX[best];
	DTF_Ds_Kplus_PY = Ds_Kp_PY[best];
	DTF_Ds_Kplus_PZ = Ds_Kp_PZ[best];
	DTF_Ds_Kplus_PE = Ds_Kp_PE[best];

	DTF_Ds_piplus_PX = Ds_pip_PX[best];
	DTF_Ds_piplus_PY = Ds_pip_PY[best];
	DTF_Ds_piplus_PZ = Ds_pip_PZ[best];
	DTF_Ds_piplus_PE = Ds_pip_PE[best];

	DTF_Ds_Kminus_PX = Ds_Km_PX[best];
	DTF_Ds_Kminus_PY = Ds_Km_PY[best];
	DTF_Ds_Kminus_PZ = Ds_Km_PZ[best];
	DTF_Ds_Kminus_PE = Ds_Km_PE[best];

        Bra_DTF_Bs_M->Fill();
        Bra_DTF_Bs_MM->Fill();
        Bra_DTF_Bs_TAU->Fill();
	Bra_DTF_Ds_M->Fill();
        Bra_DTF_chi2->Fill();

        Bra_DTF_Bs_PX->Fill();
        Bra_DTF_Bs_PY->Fill();
        Bra_DTF_Bs_PZ->Fill();
        Bra_DTF_Bs_PE->Fill();

        Bra_DTF_Ds_PX->Fill();
        Bra_DTF_Ds_PY->Fill();
        Bra_DTF_Ds_PZ->Fill();
        Bra_DTF_Ds_PE->Fill();

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

        Bra_DTF_Ds_Kplus_PX->Fill();
        Bra_DTF_Ds_Kplus_PY->Fill();
        Bra_DTF_Ds_Kplus_PZ->Fill();
        Bra_DTF_Ds_Kplus_PE->Fill();

        Bra_DTF_Ds_piplus_PX->Fill();
        Bra_DTF_Ds_piplus_PY->Fill();
        Bra_DTF_Ds_piplus_PZ->Fill();
        Bra_DTF_Ds_piplus_PE->Fill();

        Bra_DTF_Ds_Kminus_PX->Fill();
        Bra_DTF_Ds_Kminus_PY->Fill();
        Bra_DTF_Ds_Kminus_PZ->Fill();
        Bra_DTF_Ds_Kminus_PE->Fill();

     }

    new_tree->Write();
    f->Close();
}

void addVarsForBDT(string input, string output){
	
    ///Load file
    TChain* tree=new TChain("DecayTree");
    tree->Add(input.c_str());
	
    ///Reconstructed momenta
    double K_rec[4]; 
    double pip_rec[4]; 
    double pim_rec[4]; 
    double Ds_rec[4];
    double gp[6];

    tree->SetBranchAddress("K_plus_PX",&K_rec[0]);
    tree->SetBranchAddress("K_plus_PY",&K_rec[1]);
    tree->SetBranchAddress("K_plus_PZ",&K_rec[2]); 
    tree->SetBranchAddress("K_plus_PE",&K_rec[3]); 

    tree->SetBranchAddress("pi_plus_PX",&pip_rec[0]);
    tree->SetBranchAddress("pi_plus_PY",&pip_rec[1]);
    tree->SetBranchAddress("pi_plus_PZ",&pip_rec[2]); 
    tree->SetBranchAddress("pi_plus_PE",&pip_rec[3]); 

    tree->SetBranchAddress("pi_minus_PX",&pim_rec[0]);
    tree->SetBranchAddress("pi_minus_PY",&pim_rec[1]);
    tree->SetBranchAddress("pi_minus_PZ",&pim_rec[2]); 
    tree->SetBranchAddress("pi_minus_PE",&pim_rec[3]); 
    
    tree->SetBranchAddress("Ds_PX",&Ds_rec[0]);
    tree->SetBranchAddress("Ds_PY",&Ds_rec[1]);
    tree->SetBranchAddress("Ds_PZ",&Ds_rec[2]); 
    tree->SetBranchAddress("Ds_PE",&Ds_rec[3]);

    tree->SetBranchAddress( "K_plus_fromDs_TRACK_GhostProb", &gp[0] );
    tree->SetBranchAddress( "K_minus_fromDs_TRACK_GhostProb", &gp[1] );
    tree->SetBranchAddress( "pi_minus_fromDs_TRACK_GhostProb", &gp[2] );
    tree->SetBranchAddress( "K_plus_TRACK_GhostProb", &gp[3] );
    tree->SetBranchAddress( "pi_plus_TRACK_GhostProb", &gp[4] );
    tree->SetBranchAddress( "pi_minus_TRACK_GhostProb", &gp[5] );

    ///Create output file
    TFile* outputFile = new TFile(output.c_str(),"RECREATE");

    TTree* new_tree = tree->CloneTree();//CopyTree();    
    double angK,angPip, angPim, maxCos, maxGP;

    TBranch* Bra_angK = new_tree->Branch("angK", &angK, "angK/D");
    TBranch* Bra_angPip = new_tree->Branch("angPip", &angPip, "angPip/D");
    TBranch* Bra_angPim = new_tree->Branch("angPim", &angPim, "angPim/D");
    TBranch* Bra_maxCos = new_tree->Branch("maxCos", &maxCos, "maxCos/D");
    TBranch* Bra_maxGP = new_tree->Branch("maxGP", &maxGP, "maxGP/D");

    ///loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++){	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
			tree->GetEntry(i);

			TVector3 Ds(Ds_rec[0],Ds_rec[1],0.);
			TVector3 K(K_rec[0],K_rec[1],0.);
			TVector3 pip(pip_rec[0],pip_rec[1],0.);
			TVector3 pim(pim_rec[0],pim_rec[1],0.);
			angK= Ds.Angle(K);
			angPip= Ds.Angle(pip);
			angPim= Ds.Angle(pim);
			maxCos = cos(max(angK,max(angPip,angPim)));
			maxGP = max(gp[0],max(gp[1],max(gp[2],max(gp[3],max(gp[4],gp[5])))));
			Bra_maxGP->Fill();
			Bra_angK->Fill();
			Bra_angPip->Fill();
			Bra_angPim->Fill();
			Bra_maxCos->Fill();		
     }
     new_tree->Write();
     outputFile->Close();	
}


int main() {
	time_t startTime = time(0);

	/// This applies all preselection cuts except:
	/// PIDK_K_plus > 8 INN 
	/// Ds_MM signal region INN
	/// LO + HLT1 Trigger ? INN
	/// Bkg Cat for MC ?
	preselect();

        /// Add refitted momenta and variables used for BDT training
        /// Output goes to TMVAClassificationDsKpipi.C
/*
        chooseBestPV("DTF","/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_with_BDT_variables_S21_PID_temporary.root", "/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_DTF.root");
        chooseBestPV("BsDTF","/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_BsDTF.root");
        addVarsForBDT("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_forBDT.root");
*/

        chooseBestPV("DTF","/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_with_BDT_variables_S21_PID.root", "/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_DTF.root");
        chooseBestPV("BsDTF","/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_BsDTF.root");
        addVarsForBDT("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_forBDT.root");


/*
        chooseBestPV("DTF","/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_with_BDT_variables_Reco14.root",   "/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_DTF.root");
        chooseBestPV("BsDTF","/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BsDTF.root");
        addVarsForBDT("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_forBDT_Reco14.root");
*/

/*
        chooseBestPV("DTF","/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_with_BDT_variables_Reco14.root", "/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_DTF.root");
        chooseBestPV("BsDTF","/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BsDTF.root");
        addVarsForBDT("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_forBDT_Reco14.root");
*/

	/*
   	 chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Bkg/DsstKpipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Bkg/DsstKpipi_DTF.root");
   	 chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dspipipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dspipipi_DTF.root");
   	 chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dsstpipipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Bkg/Dsstpipipi_DTF.root");
	*/

    	 cout << "==============================================" << endl;
   	 cout << " Done " 
   	 << " \n Time since start " << (time(0) - startTime)/60.0
   	 << " min." << endl;
   	 cout << "==============================================" << endl; 
	
	 return 0;
}
