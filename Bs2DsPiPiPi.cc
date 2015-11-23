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


int main() {

bool MC = false; 

 TChain* tree = 0;

        tree=new TChain("Bs2Dspipipi_Ds2KKpi_Tuple/DecayTree");
//tree=new TChain("DecayTree");
        //tree->Add("/auto/data/kecke/B2DKPiPi/Data2012/Bs2DsKpipi_fullSelectionBDT.root");
        tree->Add("/auto/data/kecke/B2DKPiPi/12U-3pi-PID/*.root");
        tree->Add("/auto/data/kecke/B2DKPiPi/12D-3pi-PID/*.root");
        //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-U/*.root");
        //tree->Add("/auto/data/dargent/Bs2DsKpipi/MC/Norm/11-D/*.root");
	 //tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11D_3pi/*.root");
 	 //tree->Add("/auto/data/dargent/Bs2DsKpipi/Data_PID/11U_3pi/*.root");
	//tree->Add("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Bs2Dspipipi_with_BDT_variables_S21_PID.root");
        //tree->Add("/auto/data/kecke/B2DKPiPi/MC2012/mc2012_Ds2KKpi_preselected.root");
        //tree->Add("/auto/data/kecke/B2DKPiPi/MC2011/mc2011_Ds2KKpi_preselected.root");



 int N = tree->GetEntries();
    cout << "Old file contains " << N << " events" <<  endl;

//define variables
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;
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


Double_t Bs_MM;
Double_t Ds_MM;
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
Double_t pi_minus_fromDs_P;
Double_t K_minus_fromDs_P;
Double_t K_plus_fromDs_P;
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

//BDT variables
Double_t Bs_IPCHI2_OWNPV;
Double_t Bs_DIRA_OWNPV;
Double_t K_plus_fromDs_IPCHI2_OWNPV;
Double_t K_minus_fromDs_IPCHI2_OWNPV;
Double_t pi_minus_fromDs_IPCHI2_OWNPV;
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
Double_t pi_plus1_PT;
Double_t pi_plus2_PT;
Double_t pi_minus_PT;
Double_t K_plus_fromDs_TRACK_CHI2NDOF;
Double_t K_minus_fromDs_TRACK_CHI2NDOF;
Double_t pi_minus_fromDs_TRACK_CHI2NDOF;
Double_t pi_plus1_TRACK_CHI2NDOF;
Double_t pi_plus2_TRACK_CHI2NDOF;
Double_t pi_minus_TRACK_CHI2NDOF;


//MC truth variables
Int_t Bs_TRUEID;
Int_t Ds_TRUEID;
Int_t pi_plus1_TRUEID;
Int_t pi_plus2_TRUEID;
Int_t pi_minus_TRUEID;

Bool_t Bs_L0Global_TIS;
Bool_t Bs_L0HadronDecision_TOS;
Bool_t Bs_Hlt1TrackAllL0Decision_TOS;
Bool_t Bs_Hlt2Topo2BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo3BodyBBDTDecision_TOS;
Bool_t Bs_Hlt2Topo4BodyBBDTDecision_TOS;

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

//set branch addresses
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
tree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
tree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
tree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
tree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
tree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
tree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
tree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
tree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
tree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
tree -> SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
tree -> SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
tree -> SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );
tree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
tree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
tree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
tree -> SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
tree -> SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
tree -> SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );
tree -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P );
tree -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P );
tree -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P );
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


//BDT variables
tree -> SetBranchAddress( "Bs_IPCHI2_OWNPV" ,&Bs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "K_plus_fromDs_IPCHI2_OWNPV" ,&K_plus_fromDs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "K_minus_fromDs_IPCHI2_OWNPV" ,&K_minus_fromDs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_fromDs_IPCHI2_OWNPV" ,&pi_minus_fromDs_IPCHI2_OWNPV );
tree -> SetBranchAddress( "a_1_1260_plus_IPCHI2_OWNPV" ,&a_1_1260_plus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus1_IPCHI2_OWNPV" ,&pi_plus1_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_plus2_IPCHI2_OWNPV" ,&pi_plus2_IPCHI2_OWNPV );
tree -> SetBranchAddress( "pi_minus_IPCHI2_OWNPV" ,&pi_minus_IPCHI2_OWNPV );
tree -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
tree -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
tree -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );
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

tree -> SetBranchAddress( "K_plus_fromDs_TRACK_CHI2NDOF" , &K_plus_fromDs_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "K_minus_fromDs_TRACK_CHI2NDOF" , &K_minus_fromDs_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_minus_fromDs_TRACK_CHI2NDOF" , &pi_minus_fromDs_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_plus1_TRACK_CHI2NDOF" , &pi_plus1_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_minus_TRACK_CHI2NDOF" , &pi_minus_TRACK_CHI2NDOF );
tree -> SetBranchAddress( "pi_plus2_TRACK_CHI2NDOF" , &pi_plus2_TRACK_CHI2NDOF );


//PID variables
tree -> SetBranchAddress( "K_plus_fromDs_PIDK" , &K_plus_fromDs_PIDK );
tree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
tree -> SetBranchAddress( "pi_minus_fromDs_PIDK" , &pi_minus_fromDs_PIDK );
tree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );
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
TFile* output = 0;

        //Bs2Dspipipi normalization case
        output = new TFile("/auto/data/kecke/B2DPiPiPi/Data2012/data2012_Bs2Dspipipi_with_BDT_variables_S21_PID.root","RECREATE");


        TTree* summary_tree = tree->CloneTree(0);

        float DsDaughters_min_IPCHI2 = 0;
        float XdDaughters_min_IPCHI2 = 0;
        float DsDaughters_max_IPCHI2 = 0;
        float XdDaughters_max_IPCHI2 = 0;
        float DsDaughters_min_PT = 0;
        float XdDaughters_min_PT = 0;
        float Xd_max_DOCA = 0;
        float max_TrackChi2 = 0;
        float min_TrackChi2 = 0;


        summary_tree->Branch("DsDaughters_min_IPCHI2",&DsDaughters_min_IPCHI2,"DsDaughters_min_IPCHI2/F");
        summary_tree->Branch("XdDaughters_min_IPCHI2",&XdDaughters_min_IPCHI2,"XdDaughters_min_IPCHI2/F");
        summary_tree->Branch("DsDaughters_max_IPCHI2",&DsDaughters_max_IPCHI2,"DsDaughters_max_IPCHI2/F");
        summary_tree->Branch("XdDaughters_max_IPCHI2",&XdDaughters_max_IPCHI2,"XdDaughters_max_IPCHI2/F");
        summary_tree->Branch("DsDaughters_min_PT",&DsDaughters_min_PT,"DsDaughters_min_PT/F");
        summary_tree->Branch("XdDaughters_min_PT",&XdDaughters_min_PT,"XdDaughters_min_PT/F");
        summary_tree->Branch("Xd_max_DOCA",&Xd_max_DOCA,"Xd_max_DOCA/F");
        summary_tree->Branch("max_TrackChi2",&max_TrackChi2,"max_Track_Chi2/F");
        summary_tree->Branch("min_TrackChi2",&min_TrackChi2,"min_Track_Chi2/F");
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

//loop over events
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);


        //define the Lorentz vectors
        K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
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
	if((!Bs_Hlt2Topo2BodyBBDTDecision_TOS) &&  (!Bs_Hlt2Topo3BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo4BodyBBDTDecision_TOS)) continue;
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

	//track pt cut
 	if(pi_minus_fromDs_PT<100) continue;
	if(K_minus_fromDs_PT<100) continue;
 	if(K_plus_fromDs_PT<100) continue;
	if(pi_plus2_PT<100) continue;
	if(pi_plus1_PT<100) continue;
	if(pi_minus_PT<100) continue;

	//track p cut
	if(pi_minus_fromDs_P<1000) continue;
	if(K_minus_fromDs_P<1000) continue;
	if(K_plus_fromDs_P<1000) continue;

	//track chi2 cut
	if(K_plus_fromDs_TRACK_CHI2NDOF> 4) continue;
	if(K_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
	if(pi_minus_fromDs_TRACK_CHI2NDOF> 4) continue;
	if(pi_plus2_TRACK_CHI2NDOF> 4) continue;
	if(pi_plus1_TRACK_CHI2NDOF> 4) continue;
	if(pi_minus_TRACK_CHI2NDOF> 4) continue;

	//track ip chi2 cut
	if(K_plus_fromDs_IPCHI2_OWNPV< 4) continue; 
	if(K_minus_fromDs_IPCHI2_OWNPV< 4) continue;
	if(pi_minus_fromDs_IPCHI2_OWNPV< 4) continue;
	if(pi_plus2_IPCHI2_OWNPV< 4) continue;
	if(pi_plus1_IPCHI2_OWNPV< 4) continue;
	if(pi_minus_IPCHI2_OWNPV< 4) continue;
	

	//ds daughters pt cut
	if((pi_minus_fromDs_PT + K_minus_fromDs_PT + K_plus_fromDs_PT) < 1800) continue;

	//ds daughters doca cut
	if(Ds_DOCA1> 0.5) continue;
	if(Ds_DOCA2> 0.5) continue;
	if(Ds_DOCA3> 0.5) continue;

	//ds mass window of 100 MeV
	if( Ds_MM < 1945 || Ds_MM > 1995) continue;

	//ds vertex chi2 cut
	if((Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF) > 10) continue;

	//ds fd chi2 cut
	if(Ds_FDCHI2_OWNPV<36) continue;

	//loose pid requirements
	if(K_minus_fromDs_PIDK<-10) continue;
	if(K_plus_fromDs_PIDK<-10) continue;
	if(pi_plus2_PIDK>10) continue;
	if(pi_minus_fromDs_PIDK>10) continue;
	if(pi_plus1_PIDK>10) continue;
	if(pi_minus_PIDK>10) continue;

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

        //implement PID requirements on Ds daughters
        //cut for non PhiPi, but K*K candidates
        if( TMath::Abs((K_plus_fromDs + K_minus_fromDs).M() - massPhi) > 20 && TMath::Abs((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)< 75 && (K_plus_fromDs_PIDK<0 || K_minus_fromDs_PIDK<0) ) continue;

        //cut for non resonant candidates
        if( TMath::Abs((K_plus_fromDs + K_minus_fromDs).M() - massPhi) > 20 && TMath::Abs((pi_minus_fromDs + K_plus_fromDs).M() - massKstar) > 75 && (K_plus_fromDs_PIDK<5 || K_minus_fromDs_PIDK<5) ) continue;
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 //rejection cuts for peaking background

        //charmless suppression
        if(Ds_FDCHI2_ORIVX < 9) continue;

        //Bs->DsDs suppression
        if(TMath::Abs((pi_plus1 + pi_plus2 + pi_minus).M() - massDs) < 20) continue;

        //Bs->DsKKpi suppression
        if(pi_minus_PIDK > 0) continue;
        if(pi_plus1_PIDK > 0) continue;
        if(pi_plus2_PIDK > 0) continue;

        //Bs->D^-Kpipi suppression
        if( TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 20 && (K_minus_fromDs_PIDK < 10) ) continue;

        //Lambda_b->Lambda_c Kpipi suppression
        if( TMath::Abs((K_plus_fromDs + Kminus_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 15 && (K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 0) continue;


        //MC Truth matching
        if(MC){
                if(TMath::Abs(Bs_TRUEID) != 531) continue;
                if(TMath::Abs(Ds_TRUEID) != 431) continue;
                if(TMath::Abs(pi_plus1_TRUEID) != 211) continue;
                if(TMath::Abs(pi_plus2_TRUEID) != 211) continue;
                if(TMath::Abs(pi_minus_TRUEID) != 211) continue;
        }

        //fill all sorts of histrograms
        mass_peak_Bs->Fill(Bs_MM);
        mass_peak_Ds->Fill(Ds_MM);



	//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

 //compute the variables used for BDT training

                //min IP chi² of Ds daughters
                if( (K_plus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
                if( (K_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV < pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV < K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV < K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_min_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;

                //min IP chi² of Xd daughters
                if( (pi_plus2_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XdDaughters_min_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV < pi_plus1_IPCHI2_OWNPV) ) XdDaughters_min_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV < pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV < pi_minus_IPCHI2_OWNPV) ) XdDaughters_min_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Xd daughters
                if( (pi_plus2_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) && (pi_plus2_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XdDaughters_max_IPCHI2 = pi_plus2_IPCHI2_OWNPV;
                if( (pi_minus_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_minus_IPCHI2_OWNPV > pi_plus1_IPCHI2_OWNPV) ) XdDaughters_max_IPCHI2 = pi_minus_IPCHI2_OWNPV;
                if( (pi_plus1_IPCHI2_OWNPV > pi_plus2_IPCHI2_OWNPV) && (pi_plus1_IPCHI2_OWNPV > pi_minus_IPCHI2_OWNPV) ) XdDaughters_max_IPCHI2 = pi_plus1_IPCHI2_OWNPV;

                //max IP chi² of Ds daughters
                if( (K_plus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) && (K_plus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_plus_fromDs_IPCHI2_OWNPV;
                if( (K_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (K_minus_fromDs_IPCHI2_OWNPV > pi_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = K_minus_fromDs_IPCHI2_OWNPV;
                if( (pi_minus_fromDs_IPCHI2_OWNPV > K_plus_fromDs_IPCHI2_OWNPV) && (pi_minus_fromDs_IPCHI2_OWNPV > K_minus_fromDs_IPCHI2_OWNPV) ) DsDaughters_max_IPCHI2 = pi_minus_fromDs_IPCHI2_OWNPV;      

                //min p_t of Ds daughters
                if( (K_plus_fromDs_PT < K_minus_fromDs_PT) && (K_plus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_plus_fromDs_PT;
                if( (K_minus_fromDs_PT < K_plus_fromDs_PT) && (K_minus_fromDs_PT < pi_minus_fromDs_PT) ) DsDaughters_min_PT = K_minus_fromDs_PT;
                if( (pi_minus_fromDs_PT < K_plus_fromDs_PT) && (pi_minus_fromDs_PT < K_minus_fromDs_PT) ) DsDaughters_min_PT = pi_minus_fromDs_PT;

                //min p_t of Xd daughters
                if( (pi_plus2_PT < pi_minus_PT) && (pi_plus2_PT < pi_plus1_PT) ) XdDaughters_min_PT = pi_plus2_PT;
                if( (pi_minus_PT < pi_plus2_PT) && (pi_minus_PT < pi_plus1_PT) ) XdDaughters_min_PT = pi_minus_PT;
                if( (pi_plus1_PT < pi_plus2_PT) && (pi_plus1_PT < pi_minus_PT) ) XdDaughters_min_PT = pi_plus1_PT;


                //max DOCA of Xd
                if( (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA3 > a_1_1260_plus_DOCA1) ) Xd_max_DOCA =  a_1_1260_plus_DOCA3;
                if( (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA3) && (a_1_1260_plus_DOCA2 > a_1_1260_plus_DOCA1) ) Xd_max_DOCA =  a_1_1260_plus_DOCA2;
                if( (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA2) && (a_1_1260_plus_DOCA1 > a_1_1260_plus_DOCA3) ) Xd_max_DOCA =  a_1_1260_plus_DOCA1;

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

                summary_tree->Fill();
}


mass_peak_Bs->Draw("E1");
c->Print("eps/mass_Bs_preSel_12Data_3pi.eps");
mass_peak_Ds->Draw("E1");
c->Print("eps/mass_Ds_preSel_12Data_3pi.eps");

summary_tree->Write();
output->Close();

return 0;
}





















