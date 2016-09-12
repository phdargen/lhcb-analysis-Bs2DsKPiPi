//offline selection of BDT variables for Bs->Dspipipi
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

void preselect() {

    bool MC = false;
    bool Ds2KKpi = false;
    bool Ds2pipipi = true;

    TChain* tree = 0;
    if(Ds2pipipi) tree=new TChain("Bs2Dspipipi_Ds2pipipi_Tuple/DecayTree");
    if(Ds2KKpi) tree=new TChain("Bs2Dspipipi_Ds2KKpi_Tuple/DecayTree");
    tree->Add("/auto/data/kecke/B2DKPiPi/12U-3pi-PID/*.root");
    tree->Add("/auto/data/kecke/B2DKPiPi/12D-3pi-PID/*.root");
    //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Norm/11U/*.root");
   // tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Norm/11D/*.root");
   // tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Norm/12U/*.root");
    //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC_Reco14/Norm/12D/*.root");
    //tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D_3pi/*.root");
   // tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U_3pi/*.root");
   //tree->Add("/auto/data/kecke/B2DPiPiPi/MC2012/mc12_Bs2Dspipipi_Ds2KKpi_BDT_reweighted_Reco14.root");

    //TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_with_BDT_variables_S21_PID.root","RECREATE");
    //TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_with_BDT_variables_S21_PID_tmp.root","RECREATE");
    TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_with_BDT_variables_S21_PID.root","RECREATE");
    //TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_with_BDT_variables_S21_PID_Reco14.root","RECREATE");
     // TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_with_BDT_variables_S21_PID_Reco14.root","RECREATE");
    //TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/TriggerTisTos/mc_HLT2_TIS_pt_3545.root","RECREATE");
  // TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/MC2012/HLT2efficiency/noHLT2.root","RECREATE");
   // TFile* output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/TriggerTisTos/data_L0_TIS_pt_6575.root","RECREATE");


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
TLorentzVector pi_plus_fromDs;
TLorentzVector pi_minus2_fromDs;
TLorentzVector pi_plus1;
TLorentzVector pi_plus2;
TLorentzVector pi_minus;

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
double massLambda_b = 5619.51;


Double_t Bs_MM;
Double_t Ds_MM;
Int_t Bs_BKGCAT;
Double_t K_plus_fromDs_PX;
Double_t K_plus_fromDs_PY;
Double_t K_plus_fromDs_PZ;
Double_t pi_plus_fromDs_PX;
Double_t pi_plus_fromDs_PY;
Double_t pi_plus_fromDs_PZ;
Double_t pi_plus1_PX;
Double_t pi_plus1_PY;
Double_t pi_plus1_PZ;
Double_t K_minus_fromDs_PX;
Double_t K_minus_fromDs_PY;
Double_t K_minus_fromDs_PZ;
Double_t pi_minus_fromDs_PX;
Double_t pi_minus_fromDs_PY;
Double_t pi_minus_fromDs_PZ;
Double_t pi_minus2_fromDs_PX;
Double_t pi_minus2_fromDs_PY;
Double_t pi_minus2_fromDs_PZ;
Double_t pi_minus_PX;
Double_t pi_minus_PY;
Double_t pi_minus_PZ;
Double_t pi_plus2_PX;
Double_t pi_plus2_PY;
Double_t pi_plus2_PZ;
Double_t pi_minus_fromDs_P;
Double_t K_minus_fromDs_P;
Double_t K_plus_fromDs_P;
Double_t pi_minus2_fromDs_P;
Double_t pi_plus_fromDs_P;
Double_t pi_plus1_P;
Double_t pi_plus2_P;
Double_t pi_minus_P;

Double_t K_plus_fromDs_PIDK;
Double_t K_minus_fromDs_PIDK;
Double_t K_minus_fromDs_PIDp;
Double_t pi_minus_fromDs_PIDK;
Double_t pi_minus_PIDK;
Double_t pi_plus2_PIDK;
Double_t pi_plus1_PIDK;


Double_t Ds_FDCHI2_ORIVX;
Double_t Ds_FDCHI2_OWNPV;
Double_t a_1_1260_plus_FDCHI2_OWNPV;

Double_t Ds_ENDVERTEX_CHI2;
Double_t Ds_ENDVERTEX_Z;
Int_t Ds_ENDVERTEX_NDOF;
Double_t a_1_1260_plus_ENDVERTEX_CHI2;
Int_t a_1_1260_plus_ENDVERTEX_NDOF;
Double_t Bs_ENDVERTEX_Z;
Double_t Bs_OWNPV_Z;
Double_t Bs_ENDVERTEX_X;
Double_t Bs_OWNPV_X;
Double_t Bs_ENDVERTEX_Y;
Double_t Bs_OWNPV_Y;
Double_t Bs_FDCHI2_OWNPV;
Double_t Bs_ENDVERTEX_CHI2;
Int_t Bs_ENDVERTEX_NDOF;

Double_t pi_plus_fromDs_TRACK_GhostProb;
Double_t pi_minus2_fromDs_TRACK_GhostProb;
Double_t pi_plus_fromDs_PIDK;
Double_t pi_minus2_fromDs_PIDK;
Double_t pi_minus2_fromDs_PIDp;
Double_t pi_plus_fromDs_PIDp;
Double_t pi_minus_fromDs_PIDp;

//BDT variables
Double_t Bs_IPCHI2_OWNPV;
Double_t Bs_DIRA_OWNPV;
Double_t K_plus_fromDs_IPCHI2_OWNPV;
Double_t K_minus_fromDs_IPCHI2_OWNPV;
Double_t pi_minus_fromDs_IPCHI2_OWNPV;
Double_t pi_plus_fromDs_IPCHI2_OWNPV;
Double_t pi_minus2_fromDs_IPCHI2_OWNPV;
Double_t a_1_1260_plus_IPCHI2_OWNPV;
Double_t a_1_1260_plus_DOCA1;
Double_t a_1_1260_plus_DOCA2;
Double_t a_1_1260_plus_DOCA3;
Double_t a_1_1260_plus_ENDVERTEX_X;
Double_t a_1_1260_plus_ENDVERTEX_Y;
Double_t a_1_1260_plus_ENDVERTEX_Z;
Double_t a_1_1260_plus_OWNPV_X;
Double_t a_1_1260_plus_OWNPV_Y;
Double_t a_1_1260_plus_OWNPV_Z;
Double_t a_1_1260_plus_DIRA_OWNPV;
Double_t pi_plus1_IPCHI2_OWNPV;
Double_t pi_plus2_IPCHI2_OWNPV;
Double_t pi_minus_IPCHI2_OWNPV;
Double_t K_plus_fromDs_PT;
Double_t K_minus_fromDs_PT;
Double_t pi_minus_fromDs_PT;
Double_t pi_plus_fromDs_PT;
Double_t pi_minus2_fromDs_PT;
Double_t pi_plus1_PT;
Double_t pi_plus2_PT;
Double_t pi_minus_PT;
Double_t K_plus_fromDs_TRACK_CHI2NDOF;
Double_t K_minus_fromDs_TRACK_CHI2NDOF;
Double_t pi_minus_fromDs_TRACK_CHI2NDOF;
Double_t pi_plus_fromDs_TRACK_CHI2NDOF;
Double_t pi_minus2_fromDs_TRACK_CHI2NDOF;
Double_t pi_plus1_TRACK_CHI2NDOF;
Double_t pi_plus2_TRACK_CHI2NDOF;
Double_t pi_minus_TRACK_CHI2NDOF;

Double_t K_plus_fromDs_TRACK_GhostProb;
Double_t K_minus_fromDs_TRACK_GhostProb;
Double_t pi_minus_fromDs_TRACK_GhostProb;
Double_t pi_plus1_TRACK_GhostProb;
Double_t pi_plus2_TRACK_GhostProb;
Double_t pi_minus_TRACK_GhostProb;


//MC truth variables
Int_t Bs_TRUEID;
Int_t Ds_TRUEID;
Int_t pi_plus1_TRUEID;
Int_t pi_plus2_TRUEID;
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

//for TisTos
Bool_t Bs_L0HadronDecision_TIS;
Bool_t Bs_Hlt1TrackAllL0Decision_TIS;
Bool_t Bs_Hlt2Topo2BodyBBDTDecision_TIS;
Bool_t Bs_Hlt2Topo3BodyBBDTDecision_TIS;
Bool_t Bs_Hlt2Topo4BodyBBDTDecision_TIS;

Double_t Ds_DOCA1;
Double_t Ds_DOCA2;
Double_t Ds_DOCA3;

//trigger decisions
tree -> SetBranchAddress( "Bs_L0Global_TIS" , &Bs_L0Global_TIS );
tree -> SetBranchAddress( "Bs_L0HadronDecision_TOS" , &Bs_L0HadronDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt1TrackAllL0Decision_TOS" , &Bs_Hlt1TrackAllL0Decision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo2BodyBBDTDecision_TOS" , &Bs_Hlt2Topo2BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo3BodyBBDTDecision_TOS" , &Bs_Hlt2Topo3BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo4BodyBBDTDecision_TOS" , &Bs_Hlt2Topo4BodyBBDTDecision_TOS );

//for TisTos
tree -> SetBranchAddress( "Bs_L0HadronDecision_TIS" , &Bs_L0HadronDecision_TIS );
tree -> SetBranchAddress( "Bs_Hlt1TrackAllL0Decision_TIS" , &Bs_Hlt1TrackAllL0Decision_TIS );
tree -> SetBranchAddress( "Bs_Hlt2Topo2BodyBBDTDecision_TIS" , &Bs_Hlt2Topo2BodyBBDTDecision_TIS );
tree -> SetBranchAddress( "Bs_Hlt2Topo3BodyBBDTDecision_TIS" , &Bs_Hlt2Topo3BodyBBDTDecision_TIS );
tree -> SetBranchAddress( "Bs_Hlt2Topo4BodyBBDTDecision_TIS" , &Bs_Hlt2Topo4BodyBBDTDecision_TIS );

//set branch addresses
tree -> SetBranchAddress( "Bs_BKGCAT" , &Bs_BKGCAT );
tree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
tree -> SetBranchAddress( "Ds_TRUEID" , &Ds_TRUEID );
tree -> SetBranchAddress( "pi_plus1_TRUEID" , &pi_plus1_TRUEID );
tree -> SetBranchAddress( "pi_plus2_TRUEID" , &pi_plus2_TRUEID );
tree -> SetBranchAddress( "pi_minus_TRUEID" , &pi_minus_TRUEID );

tree -> SetBranchAddress( "Ds_DOCA1" , &Ds_DOCA1 );
tree -> SetBranchAddress( "Ds_DOCA2" , &Ds_DOCA2 );
tree -> SetBranchAddress( "Ds_DOCA3" , &Ds_DOCA3 );
tree -> SetBranchAddress( "Bs_DIRA_OWNPV" ,&Bs_DIRA_OWNPV );
tree -> SetBranchAddress( "Bs_FDCHI2_OWNPV" ,&Bs_FDCHI2_OWNPV );
tree -> SetBranchAddress( "Bs_ENDVERTEX_CHI2" ,&Bs_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "Bs_ENDVERTEX_NDOF" ,&Bs_ENDVERTEX_NDOF );
tree -> SetBranchAddress( "a_1_1260_plus_ENDVERTEX_CHI2" ,&a_1_1260_plus_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "a_1_1260_plus_ENDVERTEX_NDOF" ,&a_1_1260_plus_ENDVERTEX_NDOF );
tree -> SetBranchAddress( "a_1_1260_plus_FDCHI2_OWNPV" ,&a_1_1260_plus_FDCHI2_OWNPV );

tree -> SetBranchAddress( "Bs_MM" , &Bs_MM );
tree -> SetBranchAddress( "Ds_MM" , &Ds_MM );


if(Ds2KKpi)
{
	tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
	tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
	tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
	tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
	tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
	tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
	tree -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P );
	tree -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P );
	tree -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
	tree -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
	tree -> SetBranchAddress( "K_plus_fromDs_TRUEID" , &K_plus_fromDs_TRUEID );
	tree -> SetBranchAddress( "K_minus_fromDs_TRUEID" , &K_minus_fromDs_TRUEID );
	tree -> SetBranchAddress( "K_plus_fromDs_IPCHI2_OWNPV" ,&K_plus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "K_minus_fromDs_IPCHI2_OWNPV" ,&K_minus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "K_plus_fromDs_TRACK_CHI2NDOF" , &K_plus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "K_minus_fromDs_TRACK_CHI2NDOF" , &K_minus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "K_plus_fromDs_TRACK_GhostProb" , &K_plus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "K_minus_fromDs_TRACK_GhostProb" , &K_minus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "K_plus_fromDs_PIDK" , &K_plus_fromDs_PIDK );
	tree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
	tree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );
};
if(Ds2KKpi || Ds2pipipi)
{
	tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
	tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
	tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
	tree -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P );
	tree -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );
	tree -> SetBranchAddress( "pi_minus_fromDs_TRUEID" , &pi_minus_fromDs_TRUEID );
	tree -> SetBranchAddress( "pi_minus_fromDs_IPCHI2_OWNPV" ,&pi_minus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "pi_minus_fromDs_TRACK_CHI2NDOF" , &pi_minus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "pi_minus_fromDs_TRACK_GhostProb" , &pi_minus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "pi_minus_fromDs_PIDK" , &pi_minus_fromDs_PIDK );
	tree -> SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
}
if(Ds2pipipi)
{
	tree -> SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
	tree -> SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
	tree -> SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PX" , &pi_minus2_fromDs_PX );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PY" , &pi_minus2_fromDs_PY );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PZ" , &pi_minus2_fromDs_PZ );
	tree -> SetBranchAddress( "pi_plus_fromDs_P" , &pi_plus_fromDs_P );
	tree -> SetBranchAddress( "pi_minus2_fromDs_P" , &pi_minus2_fromDs_P );
	tree -> SetBranchAddress( "pi_plus_fromDs_PT" , &pi_plus_fromDs_PT );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PT" , &pi_minus2_fromDs_PT );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRUEID" , &pi_plus_fromDs_TRUEID );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRUEID" , &pi_minus2_fromDs_TRUEID );
	tree -> SetBranchAddress( "pi_plus_fromDs_IPCHI2_OWNPV" ,&pi_plus_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "pi_minus2_fromDs_IPCHI2_OWNPV" ,&pi_minus2_fromDs_IPCHI2_OWNPV );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRACK_CHI2NDOF" , &pi_plus_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRACK_CHI2NDOF" , &pi_minus2_fromDs_TRACK_CHI2NDOF );
	tree -> SetBranchAddress( "pi_plus_fromDs_TRACK_GhostProb" , &pi_plus_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "pi_minus2_fromDs_TRACK_GhostProb" , &pi_minus2_fromDs_TRACK_GhostProb );
	tree -> SetBranchAddress( "pi_plus_fromDs_PIDK" , &pi_plus_fromDs_PIDK );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PIDK" , &pi_minus2_fromDs_PIDK );
	tree -> SetBranchAddress( "pi_minus2_fromDs_PIDp" , &pi_minus2_fromDs_PIDp );
	tree -> SetBranchAddress( "pi_plus_fromDs_PIDp" , &pi_plus_fromDs_PIDp );
}

tree -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
tree -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
tree -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );
tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
tree -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
tree -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
tree -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );
tree -> SetBranchAddress( "pi_minus_P" , &pi_minus_P );
tree -> SetBranchAddress( "pi_plus1_P" , &pi_plus1_P );
tree -> SetBranchAddress( "pi_plus2_P" , &pi_plus2_P );
tree -> SetBranchAddress( "Bs_ENDVERTEX_X" , &Bs_ENDVERTEX_X );
tree -> SetBranchAddress( "Bs_OWNPV_X" , &Bs_OWNPV_X );
tree -> SetBranchAddress( "Bs_ENDVERTEX_Y" , &Bs_ENDVERTEX_Y );
tree -> SetBranchAddress( "Bs_OWNPV_Y" , &Bs_OWNPV_Y );
tree -> SetBranchAddress( "Bs_ENDVERTEX_Z" , &Bs_ENDVERTEX_Z );
tree -> SetBranchAddress( "Bs_OWNPV_Z" , &Bs_OWNPV_Z );
tree -> SetBranchAddress( "a_1_1260_plus_DIRA_OWNPV" ,&a_1_1260_plus_DIRA_OWNPV );
tree -> SetBranchAddress( "Ds_ENDVERTEX_Z" , &Ds_ENDVERTEX_Z );

//BDT variables
tree -> SetBranchAddress( "Bs_IPCHI2_OWNPV" ,&Bs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "a_1_1260_plus_IPCHI2_OWNPV" ,&a_1_1260_plus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus1_IPCHI2_OWNPV" ,&pi_plus1_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus2_IPCHI2_OWNPV" ,&pi_plus2_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_IPCHI2_OWNPV" ,&pi_minus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus1_PT" , &pi_plus1_PT );
tree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
tree -> SetBranchAddress( "pi_plus2_PT" , &pi_plus2_PT );
tree -> SetBranchAddress( "a_1_1260_plus_DOCA1" ,&a_1_1260_plus_DOCA1 );
tree -> SetBranchAddress( "a_1_1260_plus_DOCA2" ,&a_1_1260_plus_DOCA2 );
tree -> SetBranchAddress( "a_1_1260_plus_DOCA3" ,&a_1_1260_plus_DOCA3 );
tree -> SetBranchAddress( "a_1_1260_plus_ENDVERTEX_X" ,&a_1_1260_plus_ENDVERTEX_X);
tree -> SetBranchAddress( "a_1_1260_plus_ENDVERTEX_Y" ,&a_1_1260_plus_ENDVERTEX_Y);
tree -> SetBranchAddress( "a_1_1260_plus_ENDVERTEX_Z" ,&a_1_1260_plus_ENDVERTEX_Z);
tree -> SetBranchAddress( "a_1_1260_plus_OWNPV_X" ,&a_1_1260_plus_OWNPV_X);
tree -> SetBranchAddress( "a_1_1260_plus_OWNPV_Y" ,&a_1_1260_plus_OWNPV_Y);
tree -> SetBranchAddress( "a_1_1260_plus_OWNPV_Z" ,&a_1_1260_plus_OWNPV_Z);

tree -> SetBranchAddress( "pi_plus1_TRACK_CHI2NDOF" , &pi_plus1_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_minus_TRACK_CHI2NDOF" , &pi_minus_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_plus2_TRACK_CHI2NDOF" , &pi_plus2_TRACK_CHI2NDOF );

tree -> SetBranchAddress( "pi_plus1_TRACK_GhostProb" , &pi_plus1_TRACK_GhostProb );
tree -> SetBranchAddress( "pi_minus_TRACK_GhostProb" , &pi_minus_TRACK_GhostProb );
tree -> SetBranchAddress( "pi_plus2_TRACK_GhostProb" , &pi_plus2_TRACK_GhostProb );

//PID variables
tree -> SetBranchAddress( "pi_minus_PIDK" , &pi_minus_PIDK);
tree -> SetBranchAddress( "pi_plus1_PIDK" , &pi_plus1_PIDK);
tree -> SetBranchAddress( "pi_plus2_PIDK" , &pi_plus2_PIDK);

tree -> SetBranchAddress( "Ds_FDCHI2_ORIVX" , &Ds_FDCHI2_ORIVX );
tree -> SetBranchAddress( "Ds_FDCHI2_OWNPV" , &Ds_FDCHI2_OWNPV );
tree -> SetBranchAddress( "Ds_ENDVERTEX_CHI2" , &Ds_ENDVERTEX_CHI2 );
tree -> SetBranchAddress( "Ds_ENDVERTEX_NDOF" , &Ds_ENDVERTEX_NDOF );


TCanvas* c = new TCanvas();

//general histos
TH1D* mass_peak_Bs = new TH1D("B_{s} candidates", ";m(D_{s}#pi#pi#pi) [MeV];Entries", 75, 5000., 5600.);
TH1D* mass_peak_Ds = new TH1D("D_{s} candidates", ";m(KK#pi) [MeV];Entries", 75, 1880., 2060.);


//create outputfile
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//define some intermediate floats to compute max/min track chi2
float interMin12 = 0;
float interMin34 = 0;
float interMin56 = 0;
float interMin1to4 = 0;

float interMax12 = 0;
float interMax34 = 0;
float interMax56 = 0;
float interMax1to4 = 0;

int ptscore1 = 0;
int ptscore2 = 0;
int ptscore3 = 0;

double maxPt = 0;

//loop over events
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);

        //define the Lorentz vectors
	if(Ds2KKpi)
	{
        	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
	}
	if(Ds2KKpi || Ds2pipipi)
	{
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	}
	if(Ds2pipipi)
	{
	        pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
        	pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);
	}
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);

        Kminus_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
        Kminus_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);

	//trigger selection
	//L0 stage
	if((!Bs_L0Global_TIS) && (!Bs_L0HadronDecision_TOS)) continue;

	//HLT 1 stage
	if(!Bs_Hlt1TrackAllL0Decision_TOS) continue;

	//HLT2 stage
	if( (!Bs_Hlt2Topo2BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo3BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo4BodyBBDTDecision_TOS) ) continue;
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	///pt bins for TisTosing
	//maxPt = TMath::Max(pi_minus_PT,TMath::Max(pi_plus2_PT,TMath::Max(pi_plus1_PT,TMath::Max(pi_minus_fromDs_PT,TMath::Max(K_plus_fromDs_PT,K_minus_fromDs_PT)))));
	//if(maxPt < 3500 || maxPt > 4500) continue;
	//if(maxPt < 6500) continue;

	if(Ds2KKpi || Ds2pipipi)
	{
 		if(pi_minus_fromDs_PT<100) continue;
		if(pi_minus_fromDs_P<1600) continue;
		if(pi_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_minus_fromDs_IPCHI2_OWNPV< 4) continue;
	}
	if(Ds2KKpi)
	{
		if(K_minus_fromDs_PT<100) continue;
 		if(K_plus_fromDs_PT<100) continue;
		if((pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT) < 1800) continue;
		if(K_minus_fromDs_P<1600) continue;
		if(K_plus_fromDs_P<1600) continue;
		if(K_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_plus_fromDs_IPCHI2_OWNPV< 4) continue; 
		if(K_minus_fromDs_IPCHI2_OWNPV< 4) continue;
	}
	if(Ds2pipipi)
	{
		if(pi_minus2_fromDs_PT<100) continue;
 		if(pi_plus_fromDs_PT<100) continue;
		if((pi_minus_fromDs_PT + pi_minus2_fromDs_PT + pi_plus_fromDs_PT) < 1800) continue;
		if(pi_minus2_fromDs_P<1600) continue;
		if(pi_plus_fromDs_P<1600) continue;
		if(pi_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_minus2_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_plus_fromDs_IPCHI2_OWNPV< 4) continue; 
		if(pi_minus2_fromDs_IPCHI2_OWNPV< 4) continue;
	}

	if(pi_plus2_PT<100) continue;
	if(pi_plus1_PT<100) continue;
	if(pi_minus_PT<100) continue;

	if(pi_plus2_TRACK_CHI2NDOF> 4) continue;
	if(pi_plus1_TRACK_CHI2NDOF> 4) continue;
	if(pi_minus_TRACK_CHI2NDOF> 4) continue;

	if(pi_plus2_IPCHI2_OWNPV< 4) continue;
	if(pi_plus1_IPCHI2_OWNPV< 4) continue;
	if(pi_minus_IPCHI2_OWNPV< 4) continue;
	
	if(Ds_DOCA1> 0.5) continue;
	if(Ds_DOCA2> 0.5) continue;
	if(Ds_DOCA3> 0.5) continue;

	//ds mass window of 40 MeV
	if( Ds_MM < 1950 || Ds_MM > 1990) continue;

	//ds vertex chi2 cut
	if((Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF) > 10) continue;

	//ds fd chi2 cut
	if(Ds_FDCHI2_OWNPV<36) continue;

if(!MC){
	//loose pid requirements
	if(Ds2KKpi)
	{
		if(K_minus_fromDs_PIDK<-10) continue;
		if(K_plus_fromDs_PIDK<-10) continue;
	}
	if(Ds2KKpi || Ds2pipipi)
	{
		if(pi_minus_fromDs_PIDK>10) continue;
	}
	if(pi_plus2_PIDK>10) continue;
	if(pi_plus1_PIDK>10) continue;
	if(pi_minus_PIDK>10) continue;
}
	//cuts on Xd candidates-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//Xd daughter tracks
	if(pi_plus2_P<2000) continue;
	if(pi_plus1_P<2000) continue;
	if(pi_minus_P<2000) continue;

	//Xd min track pt

	if(pi_minus_PT < 300) ptscore1 = 1; if(pi_minus_PT > 300) ptscore1 = 0; 
	if(pi_plus1_PT < 300) ptscore2 = 1; if(pi_plus1_PT > 300) ptscore2 = 0;
	if(pi_plus2_PT < 300) ptscore3 = 1; if(pi_plus2_PT > 300) ptscore3 = 0;

	if((ptscore1 + ptscore2 + ptscore3) > 1) continue;

	//Xd DOCA cut
	if(a_1_1260_plus_DOCA1> 0.4) continue;
	if(a_1_1260_plus_DOCA2> 0.4) continue;
	if(a_1_1260_plus_DOCA3> 0.4) continue;

	//Xd pt sum cut
	if((pi_minus_PT + pi_plus1_PT + pi_plus2_PT) < 1250) continue;

	//Xd vertex chi2 cut
	if((a_1_1260_plus_ENDVERTEX_CHI2/a_1_1260_plus_ENDVERTEX_NDOF) > 8) continue;

	//Xs fd chi2
	if(a_1_1260_plus_FDCHI2_OWNPV < 16) continue;

	//vertex displacement 
	if(TMath::Abs(a_1_1260_plus_ENDVERTEX_Z - a_1_1260_plus_OWNPV_Z)<2.0) continue;
	
	if(TMath::Sqrt(((a_1_1260_plus_ENDVERTEX_X - a_1_1260_plus_OWNPV_X)*(a_1_1260_plus_ENDVERTEX_X - a_1_1260_plus_OWNPV_X))+((a_1_1260_plus_ENDVERTEX_Y - a_1_1260_plus_OWNPV_Y)*(a_1_1260_plus_ENDVERTEX_Y - a_1_1260_plus_OWNPV_Y)))<0.1) continue;
	
	//Xd DIRA
	if(a_1_1260_plus_DIRA_OWNPV < 0.98) continue;

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

	//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if((!MC) && Ds2KKpi){
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
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if(Ds2pipipi)
{
                //D^0 veto
                if((pi_plus_fromDs + pi_minus_fromDs).M() > 1700 || (pi_plus_fromDs + pi_minus2_fromDs).M() > 1700) continue;

		if(Ds_FDCHI2_ORIVX  < 9) continue;

                //PID requirements 
                if(!MC) if(pi_plus_fromDs_PIDK > 10 || pi_minus_fromDs_PIDK > 10 || pi_minus2_fromDs_PIDK > 10) continue;
                if(!MC) if(pi_plus_fromDs_PIDp > 10 || pi_minus_fromDs_PIDp > 10 || pi_minus2_fromDs_PIDp > 10) continue;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 //rejection cuts for peaking background

	if(Ds2KKpi) if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0)continue;

	//Bs->DsDs suppression
	if(TMath::Abs((pi_plus1 + pi_plus2 + pi_minus).M() - massDs) < 20)continue; 

if(!MC && Ds2KKpi){
	if(TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20){

        	//Bs->D^-Kpipi suppression
       		 if(TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 20 && (K_minus_fromDs_PIDK < 10) ) continue;

       		 ///Lambda_b->Lambda_c Kpipi suppression
        	 if(TMath::Abs((K_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 15 && (K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 0) continue;
	}

}

        //MC Truth matching
        if(MC){
		//if(Bs_BKGCAT == 50) continue; //remove ghosts
                if(TMath::Abs(Bs_TRUEID) != 531) continue;
                if(TMath::Abs(Ds_TRUEID) != 431) continue;
                if(TMath::Abs(pi_plus1_TRUEID) != 211) continue;
                if(TMath::Abs(pi_plus2_TRUEID) != 211) continue;
                if(TMath::Abs(pi_minus_TRUEID) != 211) continue;
		if(Ds2KKpi)
		{
			if(TMath::Abs(K_plus_fromDs_TRUEID) != 321) continue;
			if(TMath::Abs(K_minus_fromDs_TRUEID) != 321) continue;
		}
		if(Ds2pipipi)
		{
			if(TMath::Abs(pi_plus_fromDs_TRUEID) != 211) continue;
			if(TMath::Abs(pi_minus2_fromDs_TRUEID) != 211) continue;		
		}
		if(TMath::Abs(pi_minus_fromDs_TRUEID) != 211) continue;
        }

        //fill all sorts of histrograms
        mass_peak_Bs->Fill(Bs_MM);
        mass_peak_Ds->Fill(Ds_MM);



	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 //compute the variables used for BDT training
	if(Ds2KKpi){
                //min IP chi² of Ds daughters
                if( (K_plus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
                if( (K_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

                //min IP chi² of Xs daughters
                if( (pi_plus2_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Xs daughters
                if( (pi_plus2_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Ds daughters
                if( (K_plus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
                if( (K_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;      

                //min p_t of Ds daughters
                if( (K_plus_fromDs_PT < K_minus_fromDs_PT) && (K_plus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_plus_fromDs_PT;
                if( (K_minus_fromDs_PT < K_plus_fromDs_PT) && (K_minus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_minus_fromDs_PT;
                if( (pi_minus_fromDs_PT < K_plus_fromDs_PT) && (pi_minus_fromDs_PT < K_minus_fromDs_PT) ) DsDaughters_min_PT = pi_minus_fromDs_PT;

                //min p_t of Xs daughters
                if( (pi_plus2_PT < pi_minus_PT) && (pi_plus2_PT < pi_plus1_PT) ) XsDaughters_min_PT = pi_plus2_PT;
                if( (pi_minus_PT < pi_plus2_PT) && (pi_minus_PT < pi_plus1_PT) ) XsDaughters_min_PT = pi_minus_PT;
                if( (pi_plus1_PT < pi_plus2_PT) && (pi_plus1_PT < pi_minus_PT) ) XsDaughters_min_PT = pi_plus1_PT;


                //max DOCA of Xs
                if( (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA1) ) Xs_max_DOCA =  a_1_1260_plus_DOCA3;
                if( (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA3) && (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA1) ) Xs_max_DOCA =  a_1_1260_plus_DOCA2;
                if( (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA3) ) Xs_max_DOCA =  a_1_1260_plus_DOCA1;

                //min Track chi2
                interMin12 = TMath::Min(K_plus_fromDs_TRACK_CHI2NDOF,K_minus_fromDs_TRACK_CHI2NDOF);
                interMin34 = TMath::Min(pi_minus_fromDs_TRACK_CHI2NDOF,pi_plus2_TRACK_CHI2NDOF);
                interMin56 = TMath::Min(pi_minus_TRACK_CHI2NDOF,pi_plus1_TRACK_CHI2NDOF);
                interMin1to4 = TMath::Min(interMin12,interMin34);
                min_TrackChi2 = TMath::Min(interMin1to4,interMin56);

                //max Track chi2
                interMax12 = TMath::Max(K_plus_fromDs_TRACK_CHI2NDOF,K_minus_fromDs_TRACK_CHI2NDOF);
                interMax34 = TMath::Max(pi_minus_fromDs_TRACK_CHI2NDOF,pi_plus2_TRACK_CHI2NDOF);
                interMax56 = TMath::Max(pi_minus_TRACK_CHI2NDOF,pi_plus1_TRACK_CHI2NDOF);
                interMax1to4 = TMath::Max(interMax12,interMax34);
		max_TrackChi2 = TMath::Max(interMax1to4,interMax56);

		//max ghostProb
		max_ghostProb = TMath::Max(pi_plus1_TRACK_GhostProb,TMath::Max(pi_plus2_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(K_plus_fromDs_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb)))));
	}
if(Ds2pipipi){
                //min IP chi² of Ds daughters
                if( (pi_plus_fromDs_IPCHI2_OWNPV < pi_minus2_fromDs_IPCHI2_OWNPV) && (pi_plus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_plus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus2_fromDs_IPCHI2_OWNPV < pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus2_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus2_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV < pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV < pi_minus2_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

                //min IP chi² of Xs daughters
                if( (pi_plus2_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) ) XsDaughters_min_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Xs daughters
                if( (pi_plus2_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) ) XsDaughters_max_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Ds daughters
                if( (pi_plus_fromDs_IPCHI2_OWNPV > pi_minus2_fromDs_IPCHI2_OWNPV) && (pi_plus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_plus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus2_fromDs_IPCHI2_OWNPV > pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus2_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus2_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV > pi_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV > pi_minus2_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;      

                //min p_t of Ds daughters
                if( (pi_plus_fromDs_PT < pi_minus2_fromDs_PT) && (pi_plus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = pi_plus_fromDs_PT;
                if( (pi_minus2_fromDs_PT < pi_plus_fromDs_PT) && (pi_minus2_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = pi_minus2_fromDs_PT;
                if( (pi_minus_fromDs_PT < pi_plus_fromDs_PT) && (pi_minus_fromDs_PT < pi_minus2_fromDs_PT) ) DsDaughters_min_PT = pi_minus_fromDs_PT;

                //min p_t of Xs daughters
                if( (pi_plus2_PT < pi_minus_PT) && (pi_plus2_PT < pi_plus1_PT) ) XsDaughters_min_PT = pi_plus2_PT;
                if( (pi_minus_PT < pi_plus2_PT) && (pi_minus_PT < pi_plus1_PT) ) XsDaughters_min_PT = pi_minus_PT;
                if( (pi_plus1_PT < pi_plus2_PT) && (pi_plus1_PT < pi_minus_PT) ) XsDaughters_min_PT = pi_plus1_PT;


                //max DOCA of Xs
                if( (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA1) ) Xs_max_DOCA =  a_1_1260_plus_DOCA3;
                if( (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA3) && (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA1) ) Xs_max_DOCA =  a_1_1260_plus_DOCA2;
                if( (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA3) ) Xs_max_DOCA =  a_1_1260_plus_DOCA1;

                //min Track chi2
                interMin12 = TMath::Min(pi_plus_fromDs_TRACK_CHI2NDOF,pi_minus2_fromDs_TRACK_CHI2NDOF);
                interMin34 = TMath::Min(pi_minus_fromDs_TRACK_CHI2NDOF,pi_plus2_TRACK_CHI2NDOF);
                interMin56 = TMath::Min(pi_minus_TRACK_CHI2NDOF,pi_plus1_TRACK_CHI2NDOF);
                interMin1to4 = TMath::Min(interMin12,interMin34);
                min_TrackChi2 = TMath::Min(interMin1to4,interMin56);

                //max Track chi2
                interMax12 = TMath::Max(pi_plus_fromDs_TRACK_CHI2NDOF,pi_minus2_fromDs_TRACK_CHI2NDOF);
                interMax34 = TMath::Max(pi_minus_fromDs_TRACK_CHI2NDOF,pi_plus2_TRACK_CHI2NDOF);
                interMax56 = TMath::Max(pi_minus_TRACK_CHI2NDOF,pi_plus1_TRACK_CHI2NDOF);
                interMax1to4 = TMath::Max(interMax12,interMax34);
		max_TrackChi2 = TMath::Max(interMax1to4,interMax56);

		//max ghostProb
		max_ghostProb = TMath::Max(pi_plus1_TRACK_GhostProb,TMath::Max(pi_plus2_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(pi_plus_fromDs_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,pi_minus2_fromDs_TRACK_GhostProb)))));
	}


                summary_tree->Fill();
}


mass_peak_Bs->Draw("E1");
//c->Print("eps/mass_Bs_preSel_12MC_3pi.eps");
mass_peak_Ds->Draw("E1");
//c->Print("eps/mass_Ds_preSel_12MC_3pi.eps");

cout << "New file contains " << summary_tree->GetEntries() << " events" <<  endl;
summary_tree->Write();
output->Close();
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

    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_1_PX").c_str(),&K_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_1_PY").c_str(),&K_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_1_PZ").c_str(),&K_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_1_PE").c_str(),&K_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_PX").c_str(),&pip_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_PY").c_str(),&pip_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_PZ").c_str(),&pip_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_PE").c_str(),&pip_PE); 

    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_0_PX").c_str(),&pim_PX);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_0_PY").c_str(),&pim_PY);
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_0_PZ").c_str(),&pim_PZ); 
    tree->SetBranchAddress(("Bs_"+fit+"_a_1_1260_plus_piplus_0_PE").c_str(),&pim_PE); 

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
    double DTF_Bs_M, DTF_Ds_M, DTF_chi2, DTF_Bs_TAU; 
    double DTF_Bs_PX, DTF_Bs_PY, DTF_Bs_PZ, DTF_Bs_PE;      
    double DTF_Kplus_PX, DTF_Kplus_PY, DTF_Kplus_PZ, DTF_Kplus_PE;      
    double DTF_piplus_PX, DTF_piplus_PY, DTF_piplus_PZ, DTF_piplus_PE;      
    double DTF_piminus_PX, DTF_piminus_PY, DTF_piminus_PZ, DTF_piminus_PE;

    double DTF_Ds_Kplus_PX, DTF_Ds_Kplus_PY, DTF_Ds_Kplus_PZ, DTF_Ds_Kplus_PE;      
    double DTF_Ds_piplus_PX, DTF_Ds_piplus_PY, DTF_Ds_piplus_PZ, DTF_Ds_piplus_PE;      
    double DTF_Ds_Kminus_PX, DTF_Ds_Kminus_PY, DTF_Ds_Kminus_PZ, DTF_Ds_Kminus_PE; 
    double DTF_Ds_PX, DTF_Ds_PY, DTF_Ds_PZ, DTF_Ds_PE;      

    TBranch* Bra_DTF_Bs_M = new_tree->Branch((fit+"_Bs_M").c_str(), &DTF_Bs_M, (fit+"_Bs_M/D").c_str());
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

    tree->SetBranchAddress("pi_plus1_PX",&K_rec[0]);
    tree->SetBranchAddress("pi_plus1_PY",&K_rec[1]);
    tree->SetBranchAddress("pi_plus1_PZ",&K_rec[2]); 
    tree->SetBranchAddress("pi_plus1_PE",&K_rec[3]); 

    tree->SetBranchAddress("pi_plus2_PX",&pip_rec[0]);
    tree->SetBranchAddress("pi_plus2_PY",&pip_rec[1]);
    tree->SetBranchAddress("pi_plus2_PZ",&pip_rec[2]); 
    tree->SetBranchAddress("pi_plus2_PE",&pip_rec[3]); 

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
    tree->SetBranchAddress( "pi_plus1_TRACK_GhostProb", &gp[3] );
    tree->SetBranchAddress( "pi_plus2_TRACK_GhostProb", &gp[4] );
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
			Bra_maxCos->Fill();	
			Bra_angK->Fill();
			Bra_angPip->Fill();
			Bra_angPim->Fill();		
     }
     new_tree->Write();
     outputFile->Close();	
}

int main(){
     time_t startTime = time(0);

     /// This applies all preselection cuts except:
     /// Ds_MM signal region
     /// LO + HLT1 Trigger ?
     /// Bkg Cat for MC ?
     preselect();

     /// Add refitted momenta and variables used for BDT training
     /*
     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/Data2011/data2011_with_BDT_variables_S21_PID_temporary.root", "/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_DTF_tmp.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_DTF_tmp.root", "/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_BsDTF_tmp.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_BsDTF_tmp.root", "/auto/data/kecke/B2DPiPiPi/Data2011/data2011_Ds2KKpi_forBDT_tmp.root");
     */
/*
     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/Data2012/data2012_with_BDT_variables_S21_PID_tmp.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_DTF_tmp.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_DTF_tmp.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_BsDTF_tmp.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_BsDTF_tmp.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2KKpi_forBDT_tmp.root");
 */

     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_with_BDT_variables_S21_PID.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_DTF.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_DTF.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_BsDTF.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_BsDTF.root", "/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Ds2pipipi_forBDT.root");

/*	
     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_with_BDT_variables_S21_PID_Reco14.root", "/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_Ds2KKpi_DTF.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_Ds2KKpi_BsDTF.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DPiPiPi/MC2011/mc2011_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root");
*/
/*
     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_with_BDT_variables_S21_PID_Reco14.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_DTF.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_DTF.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_BsDTF.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_BsDTF.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2KKpi_forBDT_Reco14.root");
*/	
/*
     chooseBestPV("DTF","/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_with_BDT_variables_S21_PID_Reco14.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_DTF.root");
     chooseBestPV("BsDTF","/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_DTF.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_BsDTF.root");
     addVarsForBDT("/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_BsDTF.root", "/auto/data/kecke/B2DPiPiPi/MC2012/mc2012_Bs2Dspipipi_Ds2pipipi_forBDT_Reco14.root");
*/
     /*
     chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsstKpipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsstKpipi_DTF.root");
     chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsKpipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/DsKpipi_DTF.root");
     chooseBestPV("DTF","/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/Dsstpipipi.root", "/auto/data/dargent/Bs2DsKpipi/MC/Norm/Bkg/Dsstpipipi_DTF.root");
     */
    cout << "==============================================" << endl;
    cout << " Done " 
    << " \n Time since start " << (time(0) - startTime)/60.0
    << " min." << endl;
    cout << "==============================================" << endl;
    return 0;
}



















