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

void make_BsAsymrootTable(){

TCanvas* c= new TCanvas();
gStyle->SetOptStat(0000000000);

double pt_bins[5] = {2000., 7000., 9500., 12000., 30000.};
double y_bins[4] = {2.1, 3., 3.3, 4.5};

//2011 prod asymmetry
TH2D* Bs_prod_asym_11 = new TH2D("Bs_prod_asym", "B_{s} prod. asymmetry 2011", 4, pt_bins, 3, y_bins);
Bs_prod_asym_11->GetXaxis()->SetTitle("p_{t} [Gev/c]");
Bs_prod_asym_11->GetYaxis()->SetTitle("#eta");

Bs_prod_asym_11->SetBinContent(1, 1,  0.0166); 
Bs_prod_asym_11->SetBinError(1, 1,  TMath::Sqrt(TMath::Power(0.0632,2) + TMath::Power(0.0125,2)));

Bs_prod_asym_11->SetBinContent(1, 2, 0.0311); 
Bs_prod_asym_11->SetBinError(1, 2, TMath::Sqrt(TMath::Power(0.0773,2) + TMath::Power(0.0151,2)));

Bs_prod_asym_11->SetBinContent(1, 3, -0.0833); 
Bs_prod_asym_11->SetBinError(1, 3, TMath::Sqrt(TMath::Power(0.0558,2) + TMath::Power(0.0132,2)));


Bs_prod_asym_11->SetBinContent(2, 1, 0.0364); 
Bs_prod_asym_11->SetBinError(2, 1, TMath::Sqrt(TMath::Power(0.0479,2) + TMath::Power(0.0068,2)));

Bs_prod_asym_11->SetBinContent(2, 2, 0.0206); 
Bs_prod_asym_11->SetBinError(2, 2, TMath::Sqrt(TMath::Power(0.0682,2) + TMath::Power(0.0127,2)));

Bs_prod_asym_11->SetBinContent(2, 3, 0.0058); 
Bs_prod_asym_11->SetBinError(2, 3, TMath::Sqrt(TMath::Power(0.0584,2) + TMath::Power(0.0089,2)));


Bs_prod_asym_11->SetBinContent(3, 1, -0.0039); 
Bs_prod_asym_11->SetBinError(3, 1, TMath::Sqrt(TMath::Power(0.0456,2) + TMath::Power(0.0121,2)));

Bs_prod_asym_11->SetBinContent(3, 2, 0.1095); 
Bs_prod_asym_11->SetBinError(3, 2, TMath::Sqrt(TMath::Power(0.0723,2) + TMath::Power(0.0179,2)));

Bs_prod_asym_11->SetBinContent(3, 3, 0.1539); 
Bs_prod_asym_11->SetBinError(3, 3, TMath::Sqrt(TMath::Power(0.0722,2) + TMath::Power(0.0212,2)));


Bs_prod_asym_11->SetBinContent(4, 1, -0.0271); 
Bs_prod_asym_11->SetBinError(4, 1, TMath::Sqrt(TMath::Power(0.0336,2) + TMath::Power(0.0061,2)));

Bs_prod_asym_11->SetBinContent(4, 2, -0.0542); 
Bs_prod_asym_11->SetBinError(4, 2, TMath::Sqrt(TMath::Power(0.0612,2) + TMath::Power(0.0106,2)));

Bs_prod_asym_11->SetBinContent(4, 3, -0.0586); 
Bs_prod_asym_11->SetBinError(4, 3, TMath::Sqrt(TMath::Power(0.0648,2) + TMath::Power(0.0150,2)));

Bs_prod_asym_11->SaveAs("AsymmetryHistos/Bs_prodAsym_11.root");
Bs_prod_asym_11->Draw("COLZ");
c->Print("AsymmetryHistos/Bs_prodAsym_11.pdf");


//2012 prod asymmetry
TH2D* Bs_prod_asym_12 = new TH2D("Bs_prod_asym", "B_{s} prod. asymmetry 2012", 4, pt_bins, 3, y_bins);
Bs_prod_asym_12->GetXaxis()->SetTitle("p_{t} [Gev/c]");
Bs_prod_asym_12->GetYaxis()->SetTitle("#eta");

Bs_prod_asym_12->SetBinContent(1, 1, 0.0412); 
Bs_prod_asym_12->SetBinError(1, 1,  TMath::Sqrt(TMath::Power(0.0416,2) + TMath::Power(0.0150,2)));

Bs_prod_asym_12->SetBinContent(1, 2, -0.0241); 
Bs_prod_asym_12->SetBinError(1, 2,  TMath::Sqrt(TMath::Power(0.0574,2) + TMath::Power(0.0079,2)));

Bs_prod_asym_12->SetBinContent(1, 3, 0.0166); 
Bs_prod_asym_12->SetBinError(1, 3,  TMath::Sqrt(TMath::Power(0.0391,2) + TMath::Power(0.0092,2)));


Bs_prod_asym_12->SetBinContent(2, 1, 0.0482); 
Bs_prod_asym_12->SetBinError(2, 1,  TMath::Sqrt(TMath::Power(0.0320,2) + TMath::Power(0.0067,2)));

Bs_prod_asym_12->SetBinContent(2, 2, 0.0983); 
Bs_prod_asym_12->SetBinError(2, 2,  TMath::Sqrt(TMath::Power(0.0470,2) + TMath::Power(0.0155,2)));

Bs_prod_asym_12->SetBinContent(2, 3, -0.0430); 
Bs_prod_asym_12->SetBinError(2, 3,  TMath::Sqrt(TMath::Power(0.0386,2) + TMath::Power(0.0079,2)));


Bs_prod_asym_12->SetBinContent(3, 1, 0.0067); 
Bs_prod_asym_12->SetBinError(3, 1,  TMath::Sqrt(TMath::Power(0.0303,2) + TMath::Power(0.0063,2)));

Bs_prod_asym_12->SetBinContent(3, 2, -0.1283); 
Bs_prod_asym_12->SetBinError(3, 2,  TMath::Sqrt(TMath::Power(0.0503,2) + TMath::Power(0.0171,2)));

Bs_prod_asym_12->SetBinContent(3, 3, -0.0500); 
Bs_prod_asym_12->SetBinError(3, 3,  TMath::Sqrt(TMath::Power(0.0460,2) + TMath::Power(0.0104,2)));


Bs_prod_asym_12->SetBinContent(4, 1, -0.0012); 
Bs_prod_asym_12->SetBinError(4, 1,  TMath::Sqrt(TMath::Power(0.0222,2) + TMath::Power(0.0050,2)));

Bs_prod_asym_12->SetBinContent(4, 2, 0.0421); 
Bs_prod_asym_12->SetBinError(4, 2,  TMath::Sqrt(TMath::Power(0.0416,2) + TMath::Power(0.0162,2)));

Bs_prod_asym_12->SetBinContent(4, 3, 0.0537); 
Bs_prod_asym_12->SetBinError(4, 3,  TMath::Sqrt(TMath::Power(0.0447,2) + TMath::Power(0.0124,2)));

Bs_prod_asym_12->SaveAs("AsymmetryHistos/Bs_prodAsym_12.root");
Bs_prod_asym_12->Draw("COLZ");
c->Print("AsymmetryHistos/Bs_prodAsym_12.pdf");

}

void ComputeBsAsym(int Year){ 

//load asymmetry histo
TFile* BsAsym;
if(Year == 11) BsAsym = new TFile("AsymmetryHistos/Bs_prodAsym_11.root");
else BsAsym = new TFile("AsymmetryHistos/Bs_prodAsym_12.root");
TH2D* H_BsAsym = (TH2D*) BsAsym->Get("Bs_prod_asym");


TFile* fileNorm;
fileNorm= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
TTree* tree_Norm = (TTree*) fileNorm->Get("DecayTree");


TFile* fileSig;
fileSig= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
TTree* tree_Sig = (TTree*) fileSig->Get("DecayTree");


tree_Norm->SetBranchStatus("*",0);
tree_Sig->SetBranchStatus("*",0);

tree_Norm->SetBranchStatus("Bs_ETA",1);
tree_Norm->SetBranchStatus("Bs_PT",1);
tree_Norm->SetBranchStatus("weight",1);
tree_Norm->SetBranchStatus("N_Bs_sw",1);
tree_Norm->SetBranchStatus("year",1);
tree_Sig->SetBranchStatus("Bs_ETA",1);
tree_Sig->SetBranchStatus("Bs_PT",1);
tree_Sig->SetBranchStatus("weight",1);
tree_Sig->SetBranchStatus("N_Bs_sw",1);
tree_Sig->SetBranchStatus("year",1);


Double_t Bs_PT_Norm;
Double_t Bs_ETA_Norm;
Double_t weight_Norm;
Int_t year_Norm;

Double_t Bs_PT_Sig;
Double_t Bs_ETA_Sig;
Double_t weight_Sig;
Int_t year_Sig;


tree_Norm -> SetBranchAddress( "Bs_PT" , &Bs_PT_Norm );
tree_Norm -> SetBranchAddress( "Bs_ETA" , &Bs_ETA_Norm );
tree_Norm -> SetBranchAddress( "N_Bs_sw" , &weight_Norm );
tree_Norm -> SetBranchAddress( "year" , &year_Norm );


tree_Sig -> SetBranchAddress( "Bs_PT" , &Bs_PT_Sig );
tree_Sig -> SetBranchAddress( "Bs_ETA" , &Bs_ETA_Sig );
tree_Sig -> SetBranchAddress( "N_Bs_sw" , &weight_Sig );
tree_Sig -> SetBranchAddress( "year" , &year_Sig );


//define asymmetries for norm and signal
double A_Bs_Norm = 0;
double A_Bs_error_Norm = 0;
double weighted_A_Bs_Norm = 0;
double weighted_A_Bs_error_Norm = 0;

double A_Bs_Sig = 0;
double A_Bs_error_Sig = 0;
double weighted_A_Bs_Sig = 0;
double weighted_A_Bs_error_Sig = 0;


///loop over Norm 
int numEvents_Norm = tree_Norm->GetEntries();
for(int i=0; i< numEvents_Norm; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Norm << endl;
	tree_Norm->GetEntry(i);
	if(year_Norm != Year) continue;
	

	//get asymmetry from Histo
	A_Bs_Norm = H_BsAsym->GetBinContent(H_BsAsym->FindBin(Bs_PT_Norm,Bs_ETA_Norm));
	A_Bs_error_Norm = H_BsAsym->GetBinError(H_BsAsym->FindBin(Bs_PT_Norm,Bs_ETA_Norm));
	if(A_Bs_Norm == 0) A_Bs_error_Norm = 0;	

	//compute weighted asymmetry using sWeights
	weighted_A_Bs_Norm = weighted_A_Bs_Norm + (weight_Norm * A_Bs_Norm);
	weighted_A_Bs_error_Norm = weighted_A_Bs_error_Norm + TMath::Power((weight_Norm * A_Bs_error_Norm),2);
	}



///loop over Sig 
int numEvents_Sig = tree_Sig->GetEntries();
for(int i=0; i< numEvents_Sig; i++)
	{
	if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Sig << endl;
	tree_Sig->GetEntry(i);
	if(year_Sig != Year) continue;
	

	//get asymmetry from Histo
	A_Bs_Sig = H_BsAsym->GetBinContent(H_BsAsym->FindBin(Bs_PT_Sig,Bs_ETA_Sig));
	A_Bs_error_Sig = H_BsAsym->GetBinError(H_BsAsym->FindBin(Bs_PT_Sig,Bs_ETA_Sig));
	if(A_Bs_Sig == 0) A_Bs_error_Sig = 0;	

	//compute weighted asymmetry using sWeights
	weighted_A_Bs_Sig = weighted_A_Bs_Sig + (weight_Sig * A_Bs_Sig);
	weighted_A_Bs_error_Sig = weighted_A_Bs_error_Sig + TMath::Power((weight_Sig * A_Bs_error_Sig),2);
	}

cout << "****************************************************************************************************" << endl;

cout << " Yield of weighted Bs production asymmetry for normalization channel in your chosen dataset:  " << weighted_A_Bs_Norm << endl;
cout << " Yield of weighted Bs production asymmetry error for normalization channel in your chosen dataset:  " << TMath::Sqrt(weighted_A_Bs_error_Norm) << endl;

cout << "****************************************************************************************************" << endl;

cout << "Yield of weighted Bs production asymmetry for signal channel in your chosen dataset:  " << weighted_A_Bs_Sig  << endl;
cout << "Yield of weighted Bs production asymmetry error for signal channel in your chosen dataset:  " << TMath::Sqrt(weighted_A_Bs_error_Sig) << endl;

cout << "****************************************************************************************************" << endl;
cout << "****************************************************************************************************" << endl;

}

int main(int argc, char** argv){

//set parameters
NamedParameter<int> makeNewHistos("makeNewHistos", 0);
NamedParameter<int> year("year", 11);



if(makeNewHistos == 1)  make_BsAsymrootTable();

ComputeBsAsym(year);


return 0;
}






