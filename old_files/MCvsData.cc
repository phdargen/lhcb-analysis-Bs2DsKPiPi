//Comparing Data and MC for Bs->DsKpipi
//philippe d'argent & matthieu kecke

#include <boost/shared_ptr.hpp>
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
#include <TPaveText.h> 
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
using namespace TMath;
using namespace RooFit;

int main(int argc, char** argv)
{

//load data and mc tuples
/*
TFile *datafile = new TFile("/auto/data/kecke/B2DKPiPi/Data2011/data_Bs_11_final_sweight.root");
TTree *datatree = (TTree*)datafile->Get("DecayTree");
*/

TChain* datatree = 0;
datatree=new TChain("DecayTree");
datatree->Add("/auto/data/kecke/B2DKPiPi/Data2012/data_Bs_12_final_sweight.root");

TChain* mctree = 0;
mctree=new TChain("DecayTree");
mctree->Add("/auto/data/kecke/B2DKPiPi/mc_Ds2KKpi_with_BDT_variables.root");


//initiate variables of interest
TLorentzVector K_plus;
TLorentzVector pi_plus;
TLorentzVector pi_minus;
TLorentzVector K_plus_fromDs;
TLorentzVector K_minus_fromDs;
TLorentzVector pi_minus_fromDs;

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
Double_t Ds_PT;

Double_t K_plus_fromDs_TRACK_GhostProb;
Double_t K_minus_fromDs_TRACK_GhostProb;
Double_t pi_minus_fromDs_TRACK_GhostProb;
Double_t pi_minus_TRACK_GhostProb;
Double_t pi_plus_TRACK_GhostProb;
Double_t K_plus_TRACK_GhostProb;

Double_t K_plus_fromDs_PIDK;
Double_t K_minus_fromDs_PIDK;
Double_t K_minus_fromDs_PIDp;
Double_t pi_minus_PIDK;

Double_t Ds_FDCHI2_ORIVX;

//BDT variables
Double_t Bs_IPCHI2_OWNPV;
Double_t K_plus_fromDs_IPCHI2_OWNPV;
Double_t K_minus_fromDs_IPCHI2_OWNPV;
Double_t pi_minus_fromDs_IPCHI2_OWNPV;
Double_t K_1_1270_plus_IPCHI2_OWNPV;
Double_t K_1_1270_plus_DOCA1;
Double_t K_1_1270_plus_DOCA2;
Double_t K_1_1270_plus_DOCA3;
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

Double_t Bs_TAU;

Int_t Bs_TRUEID;
Int_t Ds_TRUEID;
Int_t K_plus_TRUEID;
Int_t pi_plus_TRUEID;
Int_t pi_minus_TRUEID;

//sWeights
Double_t N_Bs_sw;
datatree -> SetBranchAddress("N_Bs_sw" , &N_Bs_sw);



datatree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
datatree -> SetBranchAddress( "Ds_TRUEID" , &Ds_TRUEID );
datatree -> SetBranchAddress( "K_plus_TRUEID" , &K_plus_TRUEID );
datatree -> SetBranchAddress( "pi_plus_TRUEID" , &pi_plus_TRUEID );
datatree -> SetBranchAddress( "pi_minus_TRUEID" , &pi_minus_TRUEID );

datatree -> SetBranchAddress("Bs_TAU" , &Bs_TAU );

//link variables to datatree branches
datatree -> SetBranchAddress( "Bs_MM" , &Bs_MM );
datatree -> SetBranchAddress( "Ds_MM" , &Ds_MM );
datatree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
datatree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
datatree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
datatree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
datatree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
datatree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
datatree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
datatree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
datatree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
datatree -> SetBranchAddress( "K_plus_PX" , &K_plus_PX );
datatree -> SetBranchAddress( "K_plus_PY" , &K_plus_PY );
datatree -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );
datatree -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
datatree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
datatree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
datatree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
datatree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
datatree -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
datatree -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
datatree -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );
datatree -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
datatree -> SetBranchAddress( "Ds_PT" , &Ds_PT );

//BDT variables
datatree -> SetBranchAddress( "Bs_IPCHI2_OWNPV" ,&Bs_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "K_plus_fromDs_IPCHI2_OWNPV" ,&K_plus_fromDs_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "K_minus_fromDs_IPCHI2_OWNPV" ,&K_minus_fromDs_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "pi_minus_fromDs_IPCHI2_OWNPV" ,&pi_minus_fromDs_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "K_1_1270_plus_IPCHI2_OWNPV" ,&K_1_1270_plus_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "K_plus_IPCHI2_OWNPV" ,&K_plus_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "pi_plus_IPCHI2_OWNPV" ,&pi_plus_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "pi_minus_IPCHI2_OWNPV" ,&pi_minus_IPCHI2_OWNPV );
datatree -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
datatree -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
datatree -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );
datatree -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
datatree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
datatree -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
datatree -> SetBranchAddress( "K_1_1270_plus_DOCA1" ,&K_1_1270_plus_DOCA1 );
datatree -> SetBranchAddress( "K_1_1270_plus_DOCA2" ,&K_1_1270_plus_DOCA2 );
datatree -> SetBranchAddress( "K_1_1270_plus_DOCA3" ,&K_1_1270_plus_DOCA3 );
datatree -> SetBranchAddress( "K_plus_fromDs_TRACK_CHI2NDOF" , &K_plus_fromDs_TRACK_CHI2NDOF );
datatree -> SetBranchAddress( "K_minus_fromDs_TRACK_CHI2NDOF" , &K_minus_fromDs_TRACK_CHI2NDOF );
datatree -> SetBranchAddress( "pi_minus_fromDs_TRACK_CHI2NDOF" , &pi_minus_fromDs_TRACK_CHI2NDOF );
datatree -> SetBranchAddress( "K_plus_TRACK_CHI2NDOF" , &K_plus_TRACK_CHI2NDOF );
datatree -> SetBranchAddress( "pi_minus_TRACK_CHI2NDOF" , &pi_minus_TRACK_CHI2NDOF );
datatree -> SetBranchAddress( "pi_plus_TRACK_CHI2NDOF" , &pi_plus_TRACK_CHI2NDOF );

//PID variables
datatree -> SetBranchAddress( "K_plus_fromDs_PIDK" , &K_plus_fromDs_PIDK );
datatree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
datatree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );
datatree -> SetBranchAddress( "pi_minus_PIDK" , &pi_minus_PIDK);
datatree -> SetBranchAddress( "Ds_FDCHI2_ORIVX" , &Ds_FDCHI2_ORIVX );

datatree -> SetBranchAddress( "K_plus_fromDs_TRACK_GhostProb" , &K_plus_fromDs_TRACK_GhostProb);
datatree -> SetBranchAddress( "K_minus_fromDs_TRACK_GhostProb" , &K_minus_fromDs_TRACK_GhostProb);
datatree -> SetBranchAddress( "pi_minus_fromDs_TRACK_GhostProb" , &pi_minus_fromDs_TRACK_GhostProb);
datatree -> SetBranchAddress( "pi_minus_TRACK_GhostProb" , &pi_minus_TRACK_GhostProb);
datatree -> SetBranchAddress( "pi_plus_TRACK_GhostProb" , &pi_plus_TRACK_GhostProb);
datatree -> SetBranchAddress( "K_plus_TRACK_GhostProb" , &K_plus_TRACK_GhostProb);

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mctree -> SetBranchAddress( "Bs_TRUEID" , &Bs_TRUEID );
mctree -> SetBranchAddress( "Ds_TRUEID" , &Ds_TRUEID );
mctree -> SetBranchAddress( "K_plus_TRUEID" , &K_plus_TRUEID );
mctree -> SetBranchAddress( "pi_plus_TRUEID" , &pi_plus_TRUEID );
mctree -> SetBranchAddress( "pi_minus_TRUEID" , &pi_minus_TRUEID );

mctree -> SetBranchAddress("Bs_TAU" , &Bs_TAU );

//link variables to mctree branches
mctree -> SetBranchAddress( "Bs_MM" , &Bs_MM );
mctree -> SetBranchAddress( "Ds_MM" , &Ds_MM );
mctree -> SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
mctree -> SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
mctree -> SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
mctree -> SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
mctree -> SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
mctree -> SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
mctree -> SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
mctree -> SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
mctree -> SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
mctree -> SetBranchAddress( "K_plus_PX" , &K_plus_PX );
mctree -> SetBranchAddress( "K_plus_PY" , &K_plus_PY );
mctree -> SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );
mctree -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
mctree -> SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
mctree -> SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
mctree -> SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
mctree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
mctree -> SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
mctree -> SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
mctree -> SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );
mctree -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
mctree -> SetBranchAddress( "Ds_PT" , &Ds_PT );

mctree -> SetBranchAddress( "K_plus_fromDs_TRACK_GhostProb" , &K_plus_fromDs_TRACK_GhostProb);
mctree -> SetBranchAddress( "K_minus_fromDs_TRACK_GhostProb" , &K_minus_fromDs_TRACK_GhostProb);
mctree -> SetBranchAddress( "pi_minus_fromDs_TRACK_GhostProb" , &pi_minus_fromDs_TRACK_GhostProb);
mctree -> SetBranchAddress( "pi_minus_TRACK_GhostProb" , &pi_minus_TRACK_GhostProb);
mctree -> SetBranchAddress( "pi_plus_TRACK_GhostProb" , &pi_plus_TRACK_GhostProb);
mctree -> SetBranchAddress( "K_plus_TRACK_GhostProb" , &K_plus_TRACK_GhostProb);


//BDT variables
mctree -> SetBranchAddress( "Bs_IPCHI2_OWNPV" ,&Bs_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "K_plus_fromDs_IPCHI2_OWNPV" ,&K_plus_fromDs_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "K_minus_fromDs_IPCHI2_OWNPV" ,&K_minus_fromDs_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "pi_minus_fromDs_IPCHI2_OWNPV" ,&pi_minus_fromDs_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "K_1_1270_plus_IPCHI2_OWNPV" ,&K_1_1270_plus_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "K_plus_IPCHI2_OWNPV" ,&K_plus_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "pi_plus_IPCHI2_OWNPV" ,&pi_plus_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "pi_minus_IPCHI2_OWNPV" ,&pi_minus_IPCHI2_OWNPV );
mctree -> SetBranchAddress( "K_plus_fromDs_PT" , &K_plus_fromDs_PT );
mctree -> SetBranchAddress( "K_minus_fromDs_PT" , &K_minus_fromDs_PT );
mctree -> SetBranchAddress( "pi_minus_fromDs_PT" , &pi_minus_fromDs_PT );
mctree -> SetBranchAddress( "K_plus_PT" , &K_plus_PT );
mctree -> SetBranchAddress( "pi_minus_PT" , &pi_minus_PT );
mctree -> SetBranchAddress( "pi_plus_PT" , &pi_plus_PT );
mctree -> SetBranchAddress( "K_1_1270_plus_DOCA1" ,&K_1_1270_plus_DOCA1 );
mctree -> SetBranchAddress( "K_1_1270_plus_DOCA2" ,&K_1_1270_plus_DOCA2 );
mctree -> SetBranchAddress( "K_1_1270_plus_DOCA3" ,&K_1_1270_plus_DOCA3 );
mctree -> SetBranchAddress( "K_plus_fromDs_TRACK_CHI2NDOF" , &K_plus_fromDs_TRACK_CHI2NDOF );
mctree -> SetBranchAddress( "K_minus_fromDs_TRACK_CHI2NDOF" , &K_minus_fromDs_TRACK_CHI2NDOF );
mctree -> SetBranchAddress( "pi_minus_fromDs_TRACK_CHI2NDOF" , &pi_minus_fromDs_TRACK_CHI2NDOF );
mctree -> SetBranchAddress( "K_plus_TRACK_CHI2NDOF" , &K_plus_TRACK_CHI2NDOF );
mctree -> SetBranchAddress( "pi_minus_TRACK_CHI2NDOF" , &pi_minus_TRACK_CHI2NDOF );
mctree -> SetBranchAddress( "pi_plus_TRACK_CHI2NDOF" , &pi_plus_TRACK_CHI2NDOF );

//PID variables
mctree -> SetBranchAddress( "K_plus_fromDs_PIDK" , &K_plus_fromDs_PIDK );
mctree -> SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
mctree -> SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );
mctree -> SetBranchAddress( "pi_minus_PIDK" , &pi_minus_PIDK);

mctree -> SetBranchAddress( "Ds_FDCHI2_ORIVX" , &Ds_FDCHI2_ORIVX );


//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TCanvas* c = new TCanvas();
 boost::shared_ptr<TCanvas> c(new TCanvas);
gStyle->SetOptStat(1000000001);

//histos to compare distributions

//data
TH1D* m_Kpipi_data = new TH1D("K#pi#pi candidates", ";m(K#pi#pi) [MeV];Normalized Entries", 25, 500., 3000.);
TH1D* m_pipi_data = new TH1D("#pi#pi candidates", ";m(#pi#pi) [MeV];Normalized Entries", 25, 0., 3000.);
TH1D* m_Kpi_data = new TH1D("high mass K#pi candidates",";m(K#pi);Normalized Entries", 25, 0., 3000.);

TH1D* pt_Ds_data = new TH1D("D^{-}_{s} p_{t}",";D_{s} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* pt_Xs_data = new TH1D("X_{s} p_{t}",";X_{s} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);

TH1D* p_piPlus_data = new TH1D("#pi^{+} p_{t}",";#pi^{+} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* p_piMinus_data = new TH1D("#pi^{-} p_{t}",";#pi^{-} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* p_KPlus_data = new TH1D("K^{+} p_{t}",";#K^{+} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);

TH1D* max_ghostProb_data = new TH1D("max ghost prob. of tracks",";ghost prob.;Normalized Entries", 25, 0.,0.4);
TH1D* min_ghostProb_data = new TH1D("min ghost prob. of tracks",";ghost prob.;Normalized Entries", 25, 0.,1.);

//mc
TH1D* m_Kpipi_mc = new TH1D("K#pi#pi ", ";m(K#pi#pi) [MeV];Normalized Entries", 25, 500., 3000.);
TH1D* m_pipi_mc = new TH1D("#pi#pi", ";m(#pi#pi) [MeV];Normalized Entries", 25, 0., 3000.);
TH1D* m_Kpi_mc = new TH1D("high mass K#pi",";m(K#pi);Normalized Entries", 25, 0., 3000.);

TH1D* pt_Ds_mc = new TH1D("D_{s} - p_{t}",";D_{s} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* pt_Xs_mc = new TH1D("X_{s} - p_{t}",";X_{s} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);

TH1D* p_piPlus_mc = new TH1D("#pi^{+} - p_{t}",";#pi^{+} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* p_piMinus_mc = new TH1D("#pi^{-} - p_{t}",";#pi^{-} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);
TH1D* p_KPlus_mc = new TH1D("K - p_{t}",";#K^{+} p_{t} [MeV];Normalized Entries", 25, 0., 15000.);

TH1D* max_ghostProb_mc = new TH1D("max ghost prob. tracks",";ghost prob.;Normalized Entries", 25, 0.,0.4);
TH1D* min_ghostProb_mc = new TH1D("min ghost prob. tracks",";ghost prob.;Normalized Entries", 25, 0.,1.);

//define some auxilliary variables
Double_t m_Kpi_highM;

//data loop
int numEvents_data = datatree->GetEntries();
for(int i=0; i< numEvents_data; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_data << endl;
	datatree->GetEntry(i);


	//define the Lorentz vectors
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
	pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);

	m_Kpi_highM = TMath::Max((K_plus + pi_minus).M(),(K_plus + pi_plus).M());

	//fill histograms with data, weight by signal sWeight

	//mass of X_s
	m_Kpipi_data->Fill((K_plus + pi_minus + pi_plus).M(),N_Bs_sw);
	m_pipi_data->Fill((pi_minus + pi_plus).M(),N_Bs_sw);
	m_Kpi_data->Fill(m_Kpi_highM,N_Bs_sw);

	//pt of X_s
	pt_Ds_data->Fill(Ds_PT,N_Bs_sw);
	pt_Xs_data->Fill((K_plus + pi_minus + pi_plus).Pt(),N_Bs_sw);

	//pt of X_s daughters
	p_piPlus_data->Fill(pi_plus_PT,N_Bs_sw);
	p_piMinus_data->Fill(pi_minus_PT,N_Bs_sw);
	p_KPlus_data->Fill(K_plus_PT,N_Bs_sw);

	//ghost track probability
	max_ghostProb_data->Fill(TMath::Max(K_plus_TRACK_GhostProb,TMath::Max( pi_plus_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,TMath::Max( K_plus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))));
 
	min_ghostProb_data->Fill(TMath::Min(K_plus_TRACK_GhostProb,TMath::Min( pi_plus_TRACK_GhostProb,TMath::Min(pi_minus_TRACK_GhostProb,TMath::Min(pi_minus_fromDs_TRACK_GhostProb,TMath::Min( K_plus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))));

	}


//mc loop
int numEvents_mc = mctree->GetEntries();
for(int i=0; i< numEvents_mc; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_mc << endl;
	mctree->GetEntry(i);


	//define the Lorentz vectors
	K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
	K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
	K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
	pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);

	m_Kpi_highM = TMath::Max((K_plus + pi_minus).M(),(K_plus + pi_plus).M());

	//fill histograms with data

	//mass of X_s
	m_Kpipi_mc->Fill((K_plus + pi_minus + pi_plus).M());
	m_pipi_mc->Fill((pi_minus + pi_plus).M());
	m_Kpi_mc->Fill(m_Kpi_highM);

	//pt of X_s
	pt_Ds_mc->Fill(Ds_PT);
	pt_Xs_mc->Fill((K_plus + pi_minus + pi_plus).Pt());

	//pt of X_s daughters
	p_piPlus_mc->Fill(pi_plus_PT);
	p_piMinus_mc->Fill(pi_minus_PT);
	p_KPlus_mc->Fill(K_plus_PT);

	//ghost track probability
	max_ghostProb_mc->Fill(TMath::Max(K_plus_TRACK_GhostProb,TMath::Max( pi_plus_TRACK_GhostProb,TMath::Max(pi_minus_TRACK_GhostProb,TMath::Max(pi_minus_fromDs_TRACK_GhostProb,TMath::Max( K_plus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))));
 
	min_ghostProb_mc->Fill(TMath::Min(K_plus_TRACK_GhostProb,TMath::Min( pi_plus_TRACK_GhostProb,TMath::Min(pi_minus_TRACK_GhostProb,TMath::Min(pi_minus_fromDs_TRACK_GhostProb,TMath::Min( K_plus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))));

	}

//plot distrubutions and compare by Kolmogorov test 

//M(K pi pi)
double m_Kpipi_KolmoTest = m_Kpipi_data->KolmogorovTest(m_Kpipi_mc);

TPaveText *m_Kpipi_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
m_Kpipi_KolmOut->AddText(Form("Kolmogorov Test : %2f ", m_Kpipi_KolmoTest));
TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
leg->SetHeader(" ");
m_Kpipi_data->SetLineColor(kBlack);
m_Kpipi_data->Sumw2();
m_Kpipi_data->DrawNormalized("e1");
m_Kpipi_mc->SetLineColor(kBlue);
m_Kpipi_mc->Sumw2();
m_Kpipi_mc->DrawNormalized("e1SAME");
leg->AddEntry(m_Kpipi_data,"Data","LEP");
leg->AddEntry(m_Kpipi_mc,"Simulation","LEP");
m_Kpipi_KolmOut->Draw();
leg->Draw(); c->Print("eps/MC-vs-Data/m_Kpipi_comp.eps");

//M(pi pi)
double m_pipi_KolmoTest = m_pipi_data->KolmogorovTest(m_pipi_mc);

TPaveText *m_pipi_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
m_pipi_KolmOut->AddText(Form("Kolmogorov Test : %2f ", m_pipi_KolmoTest));
TLegend *leg2 = new TLegend(0.7,0.75,0.9,0.9);
leg2->SetHeader(" ");
m_pipi_data->SetLineColor(kBlack);
m_pipi_data->Sumw2();
m_pipi_data->DrawNormalized("e1");
m_pipi_mc->SetLineColor(kBlue);
m_pipi_mc->Sumw2();
m_pipi_mc->DrawNormalized("e1SAME");
leg2->AddEntry(m_pipi_data,"Data","LEP");
leg2->AddEntry(m_pipi_mc,"Simulation","LEP");
m_pipi_KolmOut->Draw();
leg2->Draw(); c->Print("eps/MC-vs-Data/m_pipi_comp.eps");

//M(K pi)
double m_Kpi_KolmoTest = m_Kpi_data->KolmogorovTest(m_Kpi_mc);

TPaveText *m_Kpi_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
m_Kpi_KolmOut->AddText(Form("Kolmogorov Test : %2f ", m_Kpi_KolmoTest));
TLegend *leg3 = new TLegend(0.7,0.75,0.9,0.9);
leg3->SetHeader(" ");
m_Kpi_data->SetLineColor(kBlack);
m_Kpi_data->Sumw2();
m_Kpi_data->DrawNormalized("e1");
m_Kpi_mc->SetLineColor(kBlue);
m_Kpi_mc->Sumw2();
m_Kpi_mc->DrawNormalized("e1SAME");
leg3->AddEntry(m_Kpi_data,"Data","LEP");
leg3->AddEntry(m_Kpi_mc,"Simulation","LEP");
m_Kpi_KolmOut->Draw();
leg3->Draw(); c->Print("eps/MC-vs-Data/m_Kpi_comp.eps");

//pt of Ds
double pt_Ds_KolmoTest = pt_Ds_data->KolmogorovTest(pt_Ds_mc);

TPaveText *pt_Ds_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
pt_Ds_KolmOut->AddText(Form("Kolmogorov Test : %2f ", pt_Ds_KolmoTest));
TLegend *leg4 = new TLegend(0.7,0.75,0.9,0.9);
leg4->SetHeader(" ");
pt_Ds_data->SetLineColor(kBlack);
pt_Ds_data->Sumw2();
pt_Ds_data->DrawNormalized("e1");
pt_Ds_mc->SetLineColor(kBlue);
pt_Ds_mc->Sumw2();
pt_Ds_mc->DrawNormalized("e1SAME");
leg4->AddEntry(pt_Ds_data,"Data","LEP");
leg4->AddEntry(pt_Ds_mc,"Simulation","LEP");
pt_Ds_KolmOut->Draw();
leg4->Draw(); c->Print("eps/MC-vs-Data/pt_Ds_comp.eps");

//pt of Xs
double pt_Xs_KolmoTest = pt_Xs_data->KolmogorovTest(pt_Xs_mc);

TPaveText *pt_Xs_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
pt_Xs_KolmOut->AddText(Form("Kolmogorov Test : %2f ", pt_Xs_KolmoTest));
TLegend *leg5 = new TLegend(0.7,0.75,0.9,0.9);
leg5->SetHeader(" ");
pt_Xs_data->SetLineColor(kBlack);
pt_Xs_data->Sumw2();
pt_Xs_data->DrawNormalized("e1");
pt_Xs_mc->SetLineColor(kBlue);
pt_Xs_mc->Sumw2();
pt_Xs_mc->DrawNormalized("e1SAME");
leg5->AddEntry(pt_Xs_data,"Data","LEP");
leg5->AddEntry(pt_Xs_mc,"Simulation","LEP");
pt_Xs_KolmOut->Draw();
leg5->Draw(); c->Print("eps/MC-vs-Data/pt_Xs_comp.eps");


//p of pi plus from Xs
double p_piPlus_KolmoTest = p_piPlus_data->KolmogorovTest(p_piPlus_mc);

TPaveText *p_piPlus_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
p_piPlus_KolmOut->AddText(Form("Kolmogorov Test : %2f ", p_piPlus_KolmoTest));
TLegend *leg6 = new TLegend(0.7,0.75,0.9,0.9);
leg6->SetHeader(" ");
p_piPlus_data->SetLineColor(kBlack);
p_piPlus_data->Sumw2();
p_piPlus_data->DrawNormalized("e1");
p_piPlus_mc->SetLineColor(kBlue);
p_piPlus_mc->Sumw2();
p_piPlus_mc->DrawNormalized("e1SAME");
leg6->AddEntry(p_piPlus_data,"Data","LEP");
leg6->AddEntry(p_piPlus_mc,"Simulation","LEP");
p_piPlus_KolmOut->Draw();
leg6->Draw(); c->Print("eps/MC-vs-Data/p_piPlus_comp.eps");

//p of pi minus from Xs
double p_piMinus_KolmoTest = p_piMinus_data->KolmogorovTest(p_piMinus_mc);

TPaveText *p_piMinus_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
p_piMinus_KolmOut->AddText(Form("Kolmogorov Test : %2f ", p_piMinus_KolmoTest));
TLegend *leg7 = new TLegend(0.7,0.75,0.9,0.9);
leg7->SetHeader(" ");
p_piMinus_data->SetLineColor(kBlack);
p_piMinus_data->Sumw2();
p_piMinus_data->DrawNormalized("e1");
p_piMinus_mc->SetLineColor(kBlue);
p_piMinus_mc->Sumw2();
p_piMinus_mc->DrawNormalized("e1SAME");
leg7->AddEntry(p_piMinus_data,"Data","LEP");
leg7->AddEntry(p_piMinus_mc,"Simulation","LEP");
p_piMinus_KolmOut->Draw();
leg7->Draw(); c->Print("eps/MC-vs-Data/p_piMinus_comp.eps");

//p of K plus from Xs
double p_KPlus_KolmoTest = p_KPlus_data->KolmogorovTest(p_KPlus_mc);

TPaveText *p_KPlus_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
p_KPlus_KolmOut->AddText(Form("Kolmogorov Test : %2f ", p_KPlus_KolmoTest));
TLegend *leg8 = new TLegend(0.7,0.75,0.9,0.9);
leg8->SetHeader(" ");
p_KPlus_data->SetLineColor(kBlack);
p_KPlus_data->Sumw2();
p_KPlus_data->DrawNormalized("e1");
p_KPlus_mc->SetLineColor(kBlue);
p_KPlus_mc->Sumw2();
p_KPlus_mc->DrawNormalized("e1SAME");
leg8->AddEntry(p_KPlus_data,"Data","LEP");
leg8->AddEntry(p_KPlus_mc,"Simulation","LEP");
p_KPlus_KolmOut->Draw();
leg8->Draw(); c->Print("eps/MC-vs-Data/p_KPlus_comp.eps");

//max track prob
double max_ghostProb_KolmoTest = max_ghostProb_data->KolmogorovTest(max_ghostProb_mc);

TPaveText *max_ghostProb_KolmOut= new TPaveText(0.5,0.75,0.7,0.9,"NDC");
max_ghostProb_KolmOut->AddText(Form("Kolmogorov Test : %2f ", max_ghostProb_KolmoTest));
TLegend *leg9 = new TLegend(0.7,0.75,0.9,0.9);
leg9->SetHeader(" ");
max_ghostProb_data->SetLineColor(kBlack);
max_ghostProb_data->Sumw2();
max_ghostProb_data->DrawNormalized("e1");
max_ghostProb_mc->SetLineColor(kBlue);
max_ghostProb_mc->Sumw2();
max_ghostProb_mc->DrawNormalized("e1SAME");
leg9->AddEntry(max_ghostProb_data,"Data","LEP");
leg9->AddEntry(max_ghostProb_mc,"Simulation","LEP");
max_ghostProb_KolmOut->Draw();
leg9->Draw(); c->Print("eps/MC-vs-Data/max_ghostProb_comp.eps");


//min track prob
TLegend *leg10 = new TLegend(0.7,0.75,0.9,0.9);
leg10->SetHeader(" ");
min_ghostProb_data->SetLineColor(kBlack);
min_ghostProb_data->Sumw2();
min_ghostProb_data->DrawNormalized("e1");
min_ghostProb_mc->SetLineColor(kBlue);
min_ghostProb_mc->Sumw2();
min_ghostProb_mc->DrawNormalized("e1SAME");
leg10->AddEntry(min_ghostProb_data,"Data","LEP");
leg10->AddEntry(min_ghostProb_mc,"Simulation","LEP");
leg10->Draw(); c->Print("eps/MC-vs-Data/min_ghostProb_comp.eps");


return 0;
}
