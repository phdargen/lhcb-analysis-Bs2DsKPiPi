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
// using namespace RooFit
using namespace RooStats;

int main() {

    //specifiy decay channel of Ds
    bool Ds2KKpi = true;
    bool Ds2Kpipi = false;
    bool Ds2pipipi = false;

    //MC or Data?
    bool MC = true;

    //which year?
    bool sevenTeV = true; //for 2011 = true , for 2012 = false 

    TChain* tree = 0;
	
    //Ds2KKpi case
    if(Ds2KKpi){
	tree=new TChain("Bs2DsKpipi_Ds2KKpi_Tuple/DecayTree");
	if((!MC) && (sevenTeV)){
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U/*.root");
	}
	if((!MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-PID/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-PID/*.root");
	}
	if((MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-MC/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-MC/*.root");
	}
	if((MC) && (sevenTeV)){
	tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/11U/*.root");
	tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/11D/*.root");
	}
    }

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
	if((!MC) && (!sevenTeV)){
		tree->Add("/auto/data/kecke/B2DKPiPi/12D-d23pi/*.root");
		tree->Add("/auto/data/kecke/B2DKPiPi/12U-d23pi/*.root");
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
    tree->SetBranchStatus("*TRACK*",1) ;
    tree->SetBranchStatus("Bs_L0Global_TIS",1) ;

    tree->SetBranchStatus( "*chi2",1 );
    tree->SetBranchStatus( "*nDOF",1 );
    tree->SetBranchStatus( "*status",1 );
    tree->SetBranchStatus( "*ctau",1 );

    if(!MC)tree->SetBranchStatus("*TRUE*",0) ;
    else tree->SetBranchStatus("*TRUE*",1) ;



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

Double_t Ds_ENDVERTEX_CHI2;
Int_t Ds_ENDVERTEX_NDOF;
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

Int_t Bs_TRUEID;
Int_t Ds_TRUEID;
Int_t K_plus_TRUEID;
Int_t pi_plus_TRUEID;
Int_t pi_minus_TRUEID;

Bool_t Bs_L0Global_TIS;
Bool_t Bs_L0HadronDecision_TOS;
Bool_t Bs_Hlt1TrackAllL0Decision_TOS;
Bool_t Bs_Hlt2Topo2BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo3BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo4BodyBBDTDecision_TOS;

Double_t pi_plus_fromDs_PIDp;
Double_t pi_minus_fromDs_PIDp;
Double_t pi_minus2_fromDs_PIDp;


//trigger decisions
tree -> SetBranchAddress( "Bs_L0Global_TIS" , &Bs_L0Global_TIS );
tree -> SetBranchAddress( "Bs_L0HadronDecision_TOS" , &Bs_L0HadronDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt1TrackAllL0Decision_TOS" , &Bs_Hlt1TrackAllL0Decision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo2BodyBBDTDecision_TOS" , &Bs_Hlt2Topo2BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo3BodyBBDTDecision_TOS" , &Bs_Hlt2Topo3BodyBBDTDecision_TOS );
tree -> SetBranchAddress( "Bs_Hlt2Topo4BodyBBDTDecision_TOS" , &Bs_Hlt2Topo4BodyBBDTDecision_TOS );

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

tree -> SetBranchAddress("Bs_TAU" , &Bs_TAU );

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
	if(Ds2KKpi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2KKpi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2KKpi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2KKpi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2KKpi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2KKpi_with_BDT_variables.root","RECREATE");
	if(Ds2KKpi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2KKpi_with_BDT_variables.root","RECREATE");

	//Ds2Kpipi case
	if(Ds2Kpipi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2Kpipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2Kpipi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2Kpipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2Kpipi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2Kpipi_with_BDT_variables.root","RECREATE");
	if(Ds2Kpipi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2Kpipi_with_BDT_variables.root","RECREATE");


	// Ds2pipipi case
	if(Ds2pipipi && (!MC) && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data2011_Ds2pipipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2pipipi && (!MC) && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/Data2012/data2012_Ds2pipipi_with_BDT_variables_S21_PID.root","RECREATE");
	if(Ds2pipipi && MC && (!sevenTeV)) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2012/mc12_Ds2pipipi_with_BDT_variables.root","RECREATE");
	if(Ds2pipipi && MC && sevenTeV) output = new TFile("/auto/data/kecke/B2DKPiPi/MC2011/mc11_Ds2pipipi_with_BDT_variables.root","RECREATE");


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


	summary_tree->Branch("DsDaughters_min_IPCHI2",&DsDaughters_min_IPCHI2,"DsDaughters_min_IPCHI2/F");
        summary_tree->Branch("XsDaughters_min_IPCHI2",&XsDaughters_min_IPCHI2,"XsDaughters_min_IPCHI2/F");
	summary_tree->Branch("DsDaughters_max_IPCHI2",&DsDaughters_max_IPCHI2,"DsDaughters_max_IPCHI2/F");
	summary_tree->Branch("XsDaughters_max_IPCHI2",&XsDaughters_max_IPCHI2,"XsDaughters_max_IPCHI2/F");
	summary_tree->Branch("DsDaughters_min_PT",&DsDaughters_min_PT,"DsDaughters_min_PT/F");
	summary_tree->Branch("XsDaughters_min_PT",&XsDaughters_min_PT,"XsDaughters_min_PT/F");
	summary_tree->Branch("Xs_max_DOCA",&Xs_max_DOCA,"Xs_max_DOCA/F");
	summary_tree->Branch("max_TrackChi2",&max_TrackChi2,"max_Track_Chi2/F");
	summary_tree->Branch("min_TrackChi2",&min_TrackChi2,"min_Track_Chi2/F");

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
		if(K_plus_fromDs_P<1000) continue;
		if(K_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_plus_fromDs_IPCHI2_OWNPV< 4) continue; 
		if((pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT) < 1800) continue;
		if(K_plus_fromDs_PIDK<-10) continue;
	}

	if(Ds2Kpipi){
		if((pi_minus_fromDs_PT + K_minus_fromDs_PT + pi_plus_fromDs_PT) < 1800) continue;
	}

	if(Ds2Kpipi || Ds2pipipi){
      		if(pi_plus_fromDs_PT<100) continue;
		if(pi_plus_fromDs_P<1000) continue;
		if(pi_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_plus_fromDs_IPCHI2_OWNPV< 4) continue;
		if(pi_plus_fromDs_PIDK>10) continue;
	}

	if(Ds2KKpi || Ds2Kpipi){
		if(K_minus_fromDs_PT<100) continue;
		if(K_minus_fromDs_P<1000) continue;
		if(K_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(K_minus_fromDs_IPCHI2_OWNPV< 4) continue;
		if(K_minus_fromDs_PIDK<-10) continue;
	}

	if(Ds2pipipi){
      		if(pi_minus2_fromDs_PT<100) continue;
		if(pi_minus2_fromDs_P<1000) continue;
		if(pi_minus2_fromDs_TRACK_CHI2NDOF> 4) continue;
		if(pi_minus2_fromDs_IPCHI2_OWNPV< 4) continue;
		if(pi_minus2_fromDs_PIDK>10) continue;
		if((pi_minus_fromDs_PT + pi_minus2_fromDs_PT + pi_plus_fromDs_PT) < 1800) continue;
	}

	//track pt cut
 	if(pi_minus_fromDs_PT<100) continue;
	if(K_plus_PT<100) continue;
	if(pi_plus_PT<100) continue;
	if(pi_minus_PT<100) continue;

	//track p cut
	if(pi_minus_fromDs_P<1000) continue;
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

	//ds mass window of 100 MeV
	if( Ds_MM < 1945 || Ds_MM > 1995) continue;

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

		//PID requirements
		if(pi_plus_fromDs_PIDK > 10 || pi_minus_fromDs_PIDK > 10 || pi_minus2_fromDs_PIDK > 10) continue;
		if(pi_plus_fromDs_PIDp > 10 || pi_minus_fromDs_PIDp > 10 || pi_minus2_fromDs_PIDp > 10) continue;
	}
	//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
		//cut for non PhiPi, but K*K candidates
		if( TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20 && TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) < 75 && (K_plus_fromDs_PIDK<0 || K_minus_fromDs_PIDK<0) ) continue;

		//cut for non resonant candidates
		if( TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20 && TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) > 75 && (K_plus_fromDs_PIDK<5 || K_minus_fromDs_PIDK<5) ) continue;
	}
	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
}
	//rejection cuts for peaking background

	//charmless suppression
	if(Ds_FDCHI2_ORIVX < 9) continue;

	//Bs->DsDs suppression
	if(TMath::Abs((K_plus + pi_plus + pi_minus).M() - massDs) < 20) continue;

if(!MC){

	//Bs->DsKKpi suppression
	if(pi_minus_PIDK > 0) continue;

	if(Ds2KKpi){
		if(TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) > 20){

			//Bs->D^-Kpipi suppression
			if( TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 20 && K_minus_fromDs_PIDK < 10 ) continue;

			//Lambda_b->Lambda_c Kpipi suppression
			if( TMath::Abs((K_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 15 && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 0) ) continue;
		}
	}
}
	//MC Truth matching
	if(MC){
		if(TMath::Abs(Bs_TRUEID) != 531) continue;
		if(TMath::Abs(Ds_TRUEID) != 431) continue;
		if(TMath::Abs(K_plus_TRUEID) != 321) continue;
		if(TMath::Abs(pi_plus_TRUEID) != 211) continue;
		if(TMath::Abs(pi_minus_TRUEID) != 211) continue;
	}

	//fill all sorts of histrograms
	ct_Bs->Fill(Bs_TAU*1000);
	mass_peak_Bs->Fill(Bs_MM);
	mass_peak_Ds->Fill(Ds_MM);


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

summary_tree->Write();
output->Close();	

return 0;
}
