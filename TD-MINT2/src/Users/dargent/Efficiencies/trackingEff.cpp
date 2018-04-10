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
#include <TMath.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
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
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
#include "Mint/RooHILLdini.h"
#include "Mint/RooHORNSdini.h"



using namespace std;
using namespace MINT;


void TrackingEff(int Year){ 

///compute tracking efficiency ratio from MC

//load tracking eff histo from tracking group
TFile* tracking;
if(Year == 11) tracking = new TFile("rootHistos/ratio2011S20MC17.root");
if(Year == 12) tracking = new TFile("rootHistos/ratio2012S20.root");
if(Year == 15) tracking = new TFile("rootHistos/Ratio_2015_25ns_Sim09b_Final_P-ETA_Final_FineInP.root");
if(Year == 16) tracking = new TFile("rootHistos/Ratio_2016_25ns_Sim09b_Long_P-ETA_FineInP.root");
TH2D* effHisto = (TH2D*) tracking->Get("Ratio");


TFile* fileMCNorm;
fileMCNorm= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm.root");
TTree* treeMC_Norm = (TTree*) fileMCNorm->Get("DecayTree");


TFile* fileMCSig;
fileMCSig= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
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
Int_t year_Norm;

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
Int_t year_Sig;

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
treeMC_Norm -> SetBranchAddress( "Bs_DTF_MM" , &Bs_DTF_M_Norm );
treeMC_Norm -> SetBranchAddress( "weight" , &weight_Norm );
treeMC_Norm -> SetBranchAddress( "year" , &year_Norm );

treeMC_Sig -> SetBranchAddress( "pi_plus_P" , &pi_plus_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_P" , &K_plus_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_P" , &pi_minus_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_fromDs_P" , &K_plus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "K_minus_fromDs_P" , &K_minus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_fromDs_P" , &pi_minus_fromDs_P_Sig );
treeMC_Sig -> SetBranchAddress( "pi_plus_ETA" , &pi_plus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_ETA" , &K_plus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_ETA" , &pi_minus_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "K_plus_fromDs_ETA" , &K_plus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "K_minus_fromDs_ETA" , &K_minus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "pi_minus_fromDs_ETA" , &pi_minus_fromDs_ETA_Sig );
treeMC_Sig -> SetBranchAddress( "Bs_DTF_MM" , &Bs_DTF_M_Sig );
treeMC_Sig -> SetBranchAddress( "weight" , &weight_Sig );
treeMC_Sig -> SetBranchAddress( "year" , &year_Sig );

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
	//only correct year
        if(year_Norm != Year) continue;	

	//get efficiency ratio for every track from Histo
	/// x-Axis units switch from GeV to MeV between Run I and II ......... 
	if(Year == 11 || Year == 12){
		eff_pi_plus1_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus1_P_Norm/1000),pi_plus1_ETA_Norm));
		eff_pi_plus1_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus1_P_Norm/1000),pi_plus1_ETA_Norm));
	}
	else{
		eff_pi_plus1_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus1_P_Norm),pi_plus1_ETA_Norm));
		eff_pi_plus1_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus1_P_Norm),pi_plus1_ETA_Norm));
	}
	if(eff_pi_plus1_Norm == 0){
		eff_pi_plus1_Norm = 1;
		eff_pi_plus1_error_Norm = 0;
	}


	if(Year == 11 || Year == 12){
		eff_pi_plus2_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus2_P_Norm/1000),pi_plus2_ETA_Norm));
		eff_pi_plus2_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus2_P_Norm/1000),pi_plus2_ETA_Norm));
	}
	else{
		eff_pi_plus2_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_plus2_P_Norm),pi_plus2_ETA_Norm));
		eff_pi_plus2_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_plus2_P_Norm),pi_plus2_ETA_Norm));
	}
	if(eff_pi_plus2_Norm == 0){
		eff_pi_plus2_Norm = 1;
		eff_pi_plus2_error_Norm = 0;
	}


	if(Year == 11 || Year == 12){	
		eff_pi_minus_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Norm/1000),pi_minus_ETA_Norm));
		eff_pi_minus_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Norm/1000),pi_minus_ETA_Norm));
	}
	else{
		eff_pi_minus_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Norm),pi_minus_ETA_Norm));
		eff_pi_minus_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Norm),pi_minus_ETA_Norm));
	}
	if(eff_pi_minus_Norm == 0){ 
		eff_pi_minus_Norm = 1;
		eff_pi_minus_error_Norm = 0;
	}


	if(Year == 11 || Year == 12){
		eff_K_plus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Norm/1000),K_plus_fromDs_ETA_Norm));
		eff_K_plus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Norm/1000),K_plus_fromDs_ETA_Norm));
	}
	else{
		eff_K_plus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Norm),K_plus_fromDs_ETA_Norm));
		eff_K_plus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Norm),K_plus_fromDs_ETA_Norm));
	}
	if(eff_K_plus_fromDs_Norm == 0){ 
		eff_K_plus_fromDs_Norm =1;
		eff_K_plus_fromDs_error_Norm =0;
	}


	if(Year == 11 || Year == 12){
		eff_K_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Norm/1000),K_minus_fromDs_ETA_Norm));
		eff_K_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Norm/1000),K_minus_fromDs_ETA_Norm));
	}
	else{
		eff_K_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Norm),K_minus_fromDs_ETA_Norm));
		eff_K_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Norm),K_minus_fromDs_ETA_Norm));
	}
	if(eff_K_minus_fromDs_Norm == 0){
		eff_K_minus_fromDs_Norm =1;
		eff_K_minus_fromDs_error_Norm =0;
	}


	if(Year == 11 || Year == 12){
		eff_pi_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Norm/1000),pi_minus_fromDs_ETA_Norm));
		eff_pi_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Norm/1000),pi_minus_fromDs_ETA_Norm));
	}
	else{
		eff_pi_minus_fromDs_Norm = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Norm),pi_minus_fromDs_ETA_Norm));
		eff_pi_minus_fromDs_error_Norm = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Norm),pi_minus_fromDs_ETA_Norm));
	}
	if(eff_pi_minus_fromDs_Norm == 0){ 
		eff_pi_minus_fromDs_Norm =1;
		eff_pi_minus_fromDs_error_Norm =0;
	}

	//compute event tracking efficiency as product of single track efficiencies
	trackingWeight_Norm = eff_pi_plus1_Norm * eff_pi_plus2_Norm * eff_pi_minus_Norm * eff_K_plus_fromDs_Norm * eff_K_minus_fromDs_Norm * eff_pi_minus_fromDs_Norm;
	trackingWeight_error_Norm = TMath::Sqrt(TMath::Power(eff_pi_plus1_error_Norm,2) + TMath::Power(eff_pi_plus2_error_Norm,2) + TMath::Power(eff_pi_minus_error_Norm,2) + TMath::Power(eff_K_plus_fromDs_error_Norm,2) + TMath::Power(eff_K_minus_fromDs_error_Norm,2) + TMath::Power(eff_pi_minus_fromDs_error_Norm,2));

	//compute weighted efficiency ratio, using the MC weights
	weighted_efficiency_Norm = weighted_efficiency_Norm + (weight_Norm * trackingWeight_Norm);
	weighted_efficiency_error_Norm = weighted_efficiency_error_Norm + TMath::Power(trackingWeight_error_Norm,2);
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
	//only correct year
        if(year_Sig != Year) continue;	

	//get efficiency ratio for every track from Histo
	/// x-Axis units switch from GeV to MeV between Run I and II ......... 

	if(Year == 11 || Year == 12){
		eff_pi_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_plus_P_Sig/1000),pi_plus_ETA_Sig));
		eff_pi_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_plus_P_Sig/1000),pi_plus_ETA_Sig));
	}
	else{
		eff_pi_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_plus_P_Sig),pi_plus_ETA_Sig));
		eff_pi_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_plus_P_Sig),pi_plus_ETA_Sig));

	}
	if(eff_pi_plus_Sig == 0){
		 eff_pi_plus_Sig = 1;
		 eff_pi_plus_error_Sig = 0;
	}

	if(Year == 11 || Year == 12){
		eff_K_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_P_Sig/1000),K_plus_ETA_Sig));
		eff_K_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_P_Sig/1000),K_plus_ETA_Sig));
	}
	else{
		eff_K_plus_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_P_Sig),K_plus_ETA_Sig));
		eff_K_plus_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_P_Sig),K_plus_ETA_Sig));

	}
	if(eff_K_plus_Sig == 0){
		eff_K_plus_Sig = 1;
		eff_K_plus_error_Sig = 0;
	}

	if(Year == 11 || Year == 12){
		eff_pi_minus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Sig/1000),pi_minus_ETA_Sig));
		eff_pi_minus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Sig/1000),pi_minus_ETA_Sig));
	}
	else{
		eff_pi_minus_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_P_Sig),pi_minus_ETA_Sig));
		eff_pi_minus_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_P_Sig),pi_minus_ETA_Sig));

	}
	if(eff_pi_minus_Sig == 0){
		eff_pi_minus_Sig = 1;
		eff_pi_minus_error_Sig = 0;		
	}

	if(Year == 11 || Year == 12){
		eff_K_plus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Sig/1000),K_plus_fromDs_ETA_Sig));
		eff_K_plus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Sig/1000),K_plus_fromDs_ETA_Sig));
	}
	else{
		eff_K_plus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_plus_fromDs_P_Sig),K_plus_fromDs_ETA_Sig));
		eff_K_plus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_plus_fromDs_P_Sig),K_plus_fromDs_ETA_Sig));

	}
	if(eff_K_plus_fromDs_Sig == 0){
		eff_K_plus_fromDs_Sig =1;
		eff_K_plus_fromDs_error_Sig =0;
	}

	if(Year == 11 || Year == 12){
		eff_K_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Sig/1000),K_minus_fromDs_ETA_Sig));
		eff_K_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Sig/1000),K_minus_fromDs_ETA_Sig));
	}
	else{
		eff_K_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((K_minus_fromDs_P_Sig),K_minus_fromDs_ETA_Sig));
		eff_K_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((K_minus_fromDs_P_Sig),K_minus_fromDs_ETA_Sig));

	}
	if(eff_K_minus_fromDs_Sig == 0){
		eff_K_minus_fromDs_Sig =1;
		eff_K_minus_fromDs_error_Sig =0;
	}

	if(Year == 11 || Year == 12){
		eff_pi_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Sig/1000),pi_minus_fromDs_ETA_Sig));
		eff_pi_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Sig/1000),pi_minus_fromDs_ETA_Sig));
	}
	else{
		eff_pi_minus_fromDs_Sig = effHisto->GetBinContent(effHisto->FindBin((pi_minus_fromDs_P_Sig),pi_minus_fromDs_ETA_Sig));
		eff_pi_minus_fromDs_error_Sig = effHisto->GetBinError(effHisto->FindBin((pi_minus_fromDs_P_Sig),pi_minus_fromDs_ETA_Sig));

	}
	if(eff_pi_minus_fromDs_Sig == 0){
		eff_pi_minus_fromDs_Sig =1;
		eff_pi_minus_fromDs_error_Sig =0;
	}

	//compute event tracking efficiency as product of single track efficiencies
	trackingWeight_Sig = eff_pi_plus_Sig * eff_K_plus_Sig * eff_pi_minus_Sig * eff_K_plus_fromDs_Sig * eff_K_minus_fromDs_Sig * eff_pi_minus_fromDs_Sig;
	trackingWeight_error_Sig = TMath::Sqrt(TMath::Power(eff_pi_plus_error_Sig,2) + TMath::Power(eff_K_plus_error_Sig,2) + TMath::Power(eff_pi_minus_error_Sig,2) + TMath::Power(eff_K_plus_fromDs_error_Sig,2) + TMath::Power(eff_K_minus_fromDs_error_Sig,2) + TMath::Power(eff_pi_minus_fromDs_error_Sig,2));

	//compute weighted efficiency ratio, using the MC weights
	weighted_efficiency_Sig = weighted_efficiency_Sig + (weight_Sig * trackingWeight_Sig);
	weighted_efficiency_error_Sig = weighted_efficiency_error_Sig + TMath::Power(trackingWeight_error_Sig,2);

	//conuter for normalization
	normalization_Sig++;
	}

TCanvas* c= new TCanvas();
gPad->SetLogx(1);
effHisto->Draw("COLZ");
if(Year == 11)c->Print("eps/TrackingEff_MCvData_11.eps");
if(Year == 12)c->Print("eps/TrackingEff_MCvData_12.eps");
if(Year == 15)c->Print("eps/TrackingEff_MCvData_15.eps");
if(Year == 16)c->Print("eps/TrackingEff_MCvData_16.eps");
gPad->SetLogx(0);

cout << "****************************************************************************************************" << endl;

cout << "normalized weighted tracking eff ratio for normalization channel:  " << weighted_efficiency_Norm / normalization_Norm << endl;
cout << "normalized weighted tracking eff ratio error for normalization channel:  " << TMath::Sqrt(weighted_efficiency_error_Norm) / normalization_Norm << endl;

cout << "****************************************************************************************************" << endl;

cout << "normalized weighted tracking eff ratio for signal channel:  " << weighted_efficiency_Sig / normalization_Sig << endl;
cout << "normalized weighted tracking eff ratio error for signal channel:  " << TMath::Sqrt(weighted_efficiency_error_Sig) / normalization_Sig << endl;

cout << "****************************************************************************************************" << endl;
cout << "****************************************************************************************************" << endl;

//cout << "Systematic uncertainty due to tracking efficiency:  " << TMath::Abs(1. - ((weighted_efficiency_Norm / normalization_Norm)/ (weighted_efficiency_Sig / normalization_Sig))) << endl;
}


int main(int argc, char** argv){

NamedParameter<int> year("year", 11);

TrackingEff(year);

return 0;
}




