
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
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;


void BkgStudies_norm_Ds2KKpi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector K_plus_fromDs;
	TLorentzVector K_minus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector K_plus;
	TLorentzVector pi_plus;
	TLorentzVector pi_minus;

	//misIDs
	TLorentzVector Kminus_asPiminus_MisID;
	TLorentzVector Kminus_asProton_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
	Double_t Bs_MM;
	Double_t Ds_MM;
	Double_t K_plus_PX;
	Double_t K_plus_PY;
	Double_t K_plus_PZ;
	Double_t pi_minus_PX;
	Double_t pi_minus_PY;
	Double_t pi_minus_PZ;
	Double_t pi_plus_PX;
	Double_t pi_plus_PY;
	Double_t pi_plus_PZ;

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

	//link to tree
        tree->SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
        tree->SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
        tree->SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
        tree->SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
        tree->SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
        tree->SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

        tree->SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
        tree->SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );


	//misID histos
	TH1D* mass_Ds2KKpi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K^{+}K^{-}_{#pi}#pi^{-}) [MeV/c^{2}];Events", 50, 1700., 1950.);
	TH1D* mass_Ds2KKpi_as_D2Kpipi_afterVeto = new TH1D("   ", ";m(K^{+}K^{-}_{#pi}#pi^{-}) [MeV/c^{2}];Events", 50, 1700., 1950.);

	TH1D* mass_Ds2KKpi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(K^{+} K^{-}_{p} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);
	TH1D* mass_Ds2KKpi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(K^{+} K^{-}_{p} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);

	TH1D* mass_Ds2KKpi_as_D02KK_beforeVeto = new TH1D(" ", ";m(K^{+} K^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2KKpi_as_D02KK_afterVeto = new TH1D("   ", ";m(K^{+} K^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);


	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);


		//misID Lorentz vectors
		Kminus_asPiminus_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
 		Kminus_asProton_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);

		//fill histos without veto cuts
		mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asPiminus_MisID).M());
		mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MisID).M());
		mass_Ds2KKpi_as_D02KK_beforeVeto->Fill((K_plus_fromDs + K_minus_fromDs).M());



		//fill histos with veto cuts
 		if( TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MisID + pi_minus_fromDs).M() - massDminus) > 30. || K_minus_fromDs_PIDK > 10 ){
			mass_Ds2KKpi_as_D2Kpipi_afterVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asPiminus_MisID).M());
		}

		if( TMath::Abs((K_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M() - massLambda_c) > 30. || ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > 5) ){
			mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MisID).M());
		}

		if((K_plus_fromDs + K_minus_fromDs).M() < 1840.){
			mass_Ds2KKpi_as_D02KK_afterVeto->Fill((K_plus_fromDs + K_minus_fromDs).M());
		}

	}

	TCanvas* c = new TCanvas();

	//plots
	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D2Kpipi_beforeVeto.eps");
	mass_Ds2KKpi_as_D2Kpipi_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D2Kpipi_afterVeto.eps");

    	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.2,0.7,0.45,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2KKpi_as_D2Kpipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2KKpi_as_D2Kpipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D2Kpipi_compareVeto.eps");



	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_Lc2KPpi_beforeVeto.eps");
	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_Lc2KPpi_afterVeto.eps");

    	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Draw("same");

    	TLegend *leg_2 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2KKpi_as_Lc2KPpi_beforeVeto,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2KKpi_as_Lc2KPpi_afterVeto,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_Lc2KPpi_compareVeto.eps");



	mass_Ds2KKpi_as_D02KK_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D02KK_beforeVeto.eps");
	mass_Ds2KKpi_as_D02KK_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D02KK_afterVeto.eps");

    	mass_Ds2KKpi_as_D02KK_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_D02KK_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_D02KK_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_D02KK_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_D02KK_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_D02KK_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.2,0.2,0.45,0.35);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2KKpi_as_D02KK_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2KKpi_as_D02KK_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2KKpi/norm_Ds2KKpi_as_D02KK_compareVeto.eps");
}

void BkgStudies_signal_Ds2KKpi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector K_plus_fromDs;
	TLorentzVector K_minus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector K_plus;
	TLorentzVector pi_plus;
	TLorentzVector pi_minus;

	//misIDs
	TLorentzVector Kminus_asPiminus_MisID;
	TLorentzVector Kminus_asProton_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
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

	//link to tree
        tree->SetBranchAddress( "K_plus_fromDs_PX" , &K_plus_fromDs_PX );
        tree->SetBranchAddress( "K_plus_fromDs_PY" , &K_plus_fromDs_PY );
        tree->SetBranchAddress( "K_plus_fromDs_PZ" , &K_plus_fromDs_PZ );
        tree->SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
        tree->SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
        tree->SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );

        tree->SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
        tree->SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );


	//misID histos
	TH1D* mass_Ds2KKpi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K^{+}K^{-}_{#pi}#pi^{-}) [MeV/c^{2}];Events", 50, 1700., 1950.);
	TH1D* mass_Ds2KKpi_as_D2Kpipi_afterVeto = new TH1D("   ", ";m(K^{+}K^{-}_{#pi}#pi^{-}) [MeV/c^{2}];Events", 50, 1700., 1950.);

	TH1D* mass_Ds2KKpi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(K^{+} K^{-}_{p} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);
	TH1D* mass_Ds2KKpi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(K^{+} K^{-}_{p} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);

	TH1D* mass_Ds2KKpi_as_D02KK_beforeVeto = new TH1D(" ", ";m(K^{+} K^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2KKpi_as_D02KK_afterVeto = new TH1D("   ", ";m(K^{+} K^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);

	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);


		//misID Lorentz vectors
		Kminus_asPiminus_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
 		Kminus_asProton_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);

		//fill histos without veto cuts
		mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asPiminus_MisID).M());
		mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MisID).M());
		mass_Ds2KKpi_as_D02KK_beforeVeto->Fill((K_plus_fromDs + K_minus_fromDs).M());



		//fill histos with veto cuts
 		if( TMath::Abs((K_plus_fromDs + Kminus_asPiminus_MisID + pi_minus_fromDs).M() - massDminus) > 30. || K_minus_fromDs_PIDK > 10 ){
			mass_Ds2KKpi_as_D2Kpipi_afterVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asPiminus_MisID).M());
		}

		if( TMath::Abs((K_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M() - massLambda_c) > 30. || ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > 5) ){
		mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Fill((K_plus_fromDs + pi_minus_fromDs + Kminus_asProton_MisID).M());
		}

		if((K_plus_fromDs + K_minus_fromDs).M() < 1840.){
		mass_Ds2KKpi_as_D02KK_afterVeto->Fill((K_plus_fromDs + K_minus_fromDs).M());
		}

	}

	TCanvas* c = new TCanvas();

	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D2Kpipi_beforeVeto.eps");
	mass_Ds2KKpi_as_D2Kpipi_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D2Kpipi_afterVeto.eps");

    	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_D2Kpipi_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_D2Kpipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.2,0.7,0.45,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2KKpi_as_D2Kpipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2KKpi_as_D2Kpipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D2Kpipi_compareVeto.eps");

	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_Lc2KPpi_beforeVeto.eps");
	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_Lc2KPpi_afterVeto.eps");

    	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_Lc2KPpi_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_Lc2KPpi_afterVeto->Draw("same");

    	TLegend *leg_2 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2KKpi_as_Lc2KPpi_beforeVeto,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2KKpi_as_Lc2KPpi_afterVeto,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_Lc2KPpi_compareVeto.eps");

	mass_Ds2KKpi_as_D02KK_beforeVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D02KK_beforeVeto.eps");
	mass_Ds2KKpi_as_D02KK_afterVeto->Draw("E1"); c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D02KK_afterVeto.eps");

    	mass_Ds2KKpi_as_D02KK_beforeVeto->SetLineColor(kRed);
    	mass_Ds2KKpi_as_D02KK_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2KKpi_as_D02KK_beforeVeto->Draw("");
    	mass_Ds2KKpi_as_D02KK_afterVeto->SetLineColor(kBlack);
    	mass_Ds2KKpi_as_D02KK_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2KKpi_as_D02KK_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.2,0.2,0.45,0.35);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2KKpi_as_D02KK_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2KKpi_as_D02KK_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2KKpi/signal_Ds2KKpi_as_D02KK_compareVeto.eps");
}

void BkgStudies_norm_Ds2Kpipi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector pi_plus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector K_minus_fromDs;

	//misIDs
	TLorentzVector piplus_asKaon_MisID;
	TLorentzVector Kminus_asProton_MisID;
	TLorentzVector Kminus_asPion_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
	Double_t pi_plus_fromDs_PX;
	Double_t pi_plus_fromDs_PY;
	Double_t pi_plus_fromDs_PZ;
	Double_t pi_minus_fromDs_PX;
	Double_t pi_minus_fromDs_PY;
	Double_t pi_minus_fromDs_PZ;
	Double_t K_minus_fromDs_PX;
	Double_t K_minus_fromDs_PY;
	Double_t K_minus_fromDs_PZ;

	Double_t pi_minus_fromDs_PIDp;
	Double_t K_minus_fromDs_PIDK;
	Double_t K_minus_fromDs_PIDp;

	//link to tree
        tree->SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
        tree->SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
        tree->SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
        tree->SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
        tree->SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
        tree->SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );

        tree->SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
        tree->SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
        tree->SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );

	//misID histos
	TH1D* mass_Ds2Kpipi_as_D2pipipi_beforeVeto = new TH1D(" ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Ds2Kpipi_as_D2pipipi_afterVeto = new TH1D("  ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}) [MeV/c^{2}];Events", 50, 1600., 2200.);

	TH1D* mass_Ds2Kpipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}_{K}) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Ds2Kpipi_as_D2Kpipi_afterVeto = new TH1D("  ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}_{K}) [MeV/c^{2}];Events", 50, 1600., 2200.);

	TH1D* mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 3000.);
	TH1D* mass_Ds2Kpipi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(#pi^{+}_{K} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 3000.);

	TH1D* mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 2800.);
	TH1D* mass_Ds2Kpipi_as_Lc2piPpi_afterVeto = new TH1D("   ", ";m(#pi^{+} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 2800.);

	TH1D* mass_Ds2Kpipi_as_D02Kpi_beforeVeto = new TH1D(" ", ";m(K^{-} #pi^{+}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2Kpipi_as_D02Kpi_afterVeto = new TH1D("   ", ";m(K^{-} #pi^{+}) [MeV/c^{2}];Events", 50, 1600., 1900.);



	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);


		//misID Lorentz vectors
		piplus_asKaon_MisID.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massKaon);
		Kminus_asProton_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massProton);
		Kminus_asPion_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);


		//fill histos without veto cuts
		mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Fill((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		mass_Ds2Kpipi_as_D2Kpipi_beforeVeto->Fill((Kminus_asPion_MisID + piplus_asKaon_MisID +  pi_minus_fromDs).M());
		mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto->Fill((Kminus_asProton_MisID + piplus_asKaon_MisID + pi_minus_fromDs).M());
		mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Fill((Kminus_asProton_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Fill((K_minus_fromDs + pi_plus_fromDs).M());

		//fill histos with veto cuts
		if(TMath::Abs((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M() - massDminus ) > 30. || K_minus_fromDs_PIDK > 20.){
		mass_Ds2Kpipi_as_D2pipipi_afterVeto->Fill((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M() - massLambda_c) > 30. || (K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > 5.){
		mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Fill((pi_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M());
		}
		if((pi_plus_fromDs + K_minus_fromDs).M() < 1750.){
		mass_Ds2Kpipi_as_D02Kpi_afterVeto->Fill((K_minus_fromDs + pi_plus_fromDs).M());
		}

	}

	TCanvas* c = new TCanvas();

	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D2pipipi_beforeVeto.eps");
	mass_Ds2Kpipi_as_D2pipipi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D2pipipi_afterVeto.eps");

    	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2Kpipi_as_D2pipipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2Kpipi_as_D2pipipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D2pipipi_compareVeto.eps");


	mass_Ds2Kpipi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D2Kpipi_beforeVeto.eps");

	mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_Lc2KPpi_beforeVeto.eps");


	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_Lc2piPpi_beforeVeto.eps");
	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_Lc2piPpi_afterVeto.eps");

    	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Draw("same");

    	TLegend *leg_2 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2Kpipi_as_Lc2piPpi_afterVeto,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_Lc2piPpi_compareVeto.eps");



	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D02Kpi_beforeVeto.eps");
	mass_Ds2Kpipi_as_D02Kpi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D02Kpi_afterVeto.eps");

    	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2Kpipi_as_D02Kpi_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2Kpipi_as_D02Kpi_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2Kpipi/norm_Ds2Kpipi_as_D02Kpi_compareVeto.eps");

}

void BkgStudies_signal_Ds2Kpipi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector pi_plus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector K_minus_fromDs;

	//misIDs
	TLorentzVector piplus_asKaon_MisID;
	TLorentzVector Kminus_asProton_MisID;
	TLorentzVector Kminus_asPion_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
	Double_t pi_plus_fromDs_PX;
	Double_t pi_plus_fromDs_PY;
	Double_t pi_plus_fromDs_PZ;
	Double_t pi_minus_fromDs_PX;
	Double_t pi_minus_fromDs_PY;
	Double_t pi_minus_fromDs_PZ;
	Double_t K_minus_fromDs_PX;
	Double_t K_minus_fromDs_PY;
	Double_t K_minus_fromDs_PZ;

	Double_t pi_minus_fromDs_PIDp;
	Double_t K_minus_fromDs_PIDK;
	Double_t K_minus_fromDs_PIDp;

	//link to tree
        tree->SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
        tree->SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
        tree->SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
        tree->SetBranchAddress( "K_minus_fromDs_PX" , &K_minus_fromDs_PX );
        tree->SetBranchAddress( "K_minus_fromDs_PY" , &K_minus_fromDs_PY );
        tree->SetBranchAddress( "K_minus_fromDs_PZ" , &K_minus_fromDs_PZ );

        tree->SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
        tree->SetBranchAddress( "K_minus_fromDs_PIDK" , &K_minus_fromDs_PIDK );
        tree->SetBranchAddress( "K_minus_fromDs_PIDp" , &K_minus_fromDs_PIDp );

	//misID histos
	TH1D* mass_Ds2Kpipi_as_D2pipipi_beforeVeto = new TH1D(" ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Ds2Kpipi_as_D2pipipi_afterVeto = new TH1D("  ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}) [MeV/c^{2}];Events", 50, 1600., 2200.);

	TH1D* mass_Ds2Kpipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}_{K}) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Ds2Kpipi_as_D2Kpipi_afterVeto = new TH1D("  ", ";m(K^{-}_{#pi}#pi^{-}#pi^{+}_{K}) [MeV/c^{2}];Events", 50, 1600., 2200.);

	TH1D* mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 3000.);
	TH1D* mass_Ds2Kpipi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(#pi^{+}_{K} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 3000.);

	TH1D* mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 2800.);
	TH1D* mass_Ds2Kpipi_as_Lc2piPpi_afterVeto = new TH1D("   ", ";m(#pi^{+} K^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2000., 2800.);

	TH1D* mass_Ds2Kpipi_as_D02Kpi_beforeVeto = new TH1D(" ", ";m(K^{-} #pi^{+}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2Kpipi_as_D02Kpi_afterVeto = new TH1D("   ", ";m(K^{-} #pi^{+}) [MeV/c^{2}];Events", 50, 1600., 1900.);



	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
		K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);


		//misID Lorentz vectors
		piplus_asKaon_MisID.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massKaon);
		Kminus_asProton_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massProton);
		Kminus_asPion_MisID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);


		//fill histos without veto cuts
		mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Fill((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		mass_Ds2Kpipi_as_D2Kpipi_beforeVeto->Fill((Kminus_asPion_MisID + piplus_asKaon_MisID +  pi_minus_fromDs).M());
		mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto->Fill((Kminus_asProton_MisID + piplus_asKaon_MisID + pi_minus_fromDs).M());
		mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Fill((Kminus_asProton_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Fill((K_minus_fromDs + pi_plus_fromDs).M());

		//fill histos with veto cuts
		if(TMath::Abs((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M() - massDminus ) > 30. || K_minus_fromDs_PIDK > 20.){
		mass_Ds2Kpipi_as_D2pipipi_afterVeto->Fill((Kminus_asPion_MisID + pi_minus_fromDs + pi_plus_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M() - massLambda_c) > 30. || (K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) > 5.){
		mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Fill((pi_plus_fromDs + Kminus_asProton_MisID + pi_minus_fromDs).M());
		}
		if((pi_plus_fromDs + K_minus_fromDs).M() < 1750.){
		mass_Ds2Kpipi_as_D02Kpi_afterVeto->Fill((K_minus_fromDs + pi_plus_fromDs).M());
		}

	}

	TCanvas* c = new TCanvas();

	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D2pipipi_beforeVeto.eps");
	mass_Ds2Kpipi_as_D2pipipi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D2pipipi_afterVeto.eps");

    	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_D2pipipi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_D2pipipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2Kpipi_as_D2pipipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2Kpipi_as_D2pipipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D2pipipi_compareVeto.eps");



	mass_Ds2Kpipi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D2Kpipi_beforeVeto.eps");

	mass_Ds2Kpipi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_Lc2KPpi_beforeVeto.eps");



	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_Lc2piPpi_beforeVeto.eps");
	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_Lc2piPpi_afterVeto.eps");

    	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_Lc2piPpi_afterVeto->Draw("same");

    	TLegend *leg_2 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2Kpipi_as_Lc2piPpi_beforeVeto,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2Kpipi_as_Lc2piPpi_afterVeto,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_Lc2piPpi_compareVeto.eps");



	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D02Kpi_beforeVeto.eps");
	mass_Ds2Kpipi_as_D02Kpi_afterVeto->Draw("E1"); c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D02Kpi_afterVeto.eps");

    	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2Kpipi_as_D02Kpi_beforeVeto->Draw("");
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2Kpipi_as_D02Kpi_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2Kpipi_as_D02Kpi_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2Kpipi_as_D02Kpi_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2Kpipi/signal_Ds2Kpipi_as_D02Kpi_compareVeto.eps");

}

void BkgStudies_norm_Ds2pipipi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector pi_plus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector pi_minus2_fromDs;

	//misIDs
	TLorentzVector piplus_asKaon_MisID;
	TLorentzVector piminus_asProton_MisID;
	TLorentzVector piminus2_asProton_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
	Double_t pi_plus_fromDs_PX;
	Double_t pi_plus_fromDs_PY;
	Double_t pi_plus_fromDs_PZ;
	Double_t pi_minus_fromDs_PX;
	Double_t pi_minus_fromDs_PY;
	Double_t pi_minus_fromDs_PZ;
	Double_t pi_minus2_fromDs_PX;
	Double_t pi_minus2_fromDs_PY;
	Double_t pi_minus2_fromDs_PZ;

	Double_t pi_minus_fromDs_PIDp;
	Double_t pi_minus2_fromDs_PIDp;

	//link to tree
        tree->SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
        tree->SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
        tree->SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus2_fromDs_PX" , &pi_minus2_fromDs_PX );
        tree->SetBranchAddress( "pi_minus2_fromDs_PY" , &pi_minus2_fromDs_PY );
        tree->SetBranchAddress( "pi_minus2_fromDs_PZ" , &pi_minus2_fromDs_PZ );

        tree->SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
        tree->SetBranchAddress( "pi_minus2_fromDs_PIDp" , &pi_minus2_fromDs_PIDp );


	//misID histos
	TH1D* mass_Ds2pipipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 1900., 3000.);

	TH1D* mass_Ds2pipipi_as_D02pipi_beforeVeto = new TH1D(" ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_afterVeto = new TH1D("   ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_beforeVeto_2 = new TH1D("     ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_afterVeto_2 = new TH1D("        ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);

	TH1D* mass_Ds2pipipi_as_Lc2piPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_afterVeto = new TH1D("   ", ";m(#pi^{+} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2 = new TH1D("     ", ";m(#pi^{+}  #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2 = new TH1D("       ", ";m(#pi^{+} #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2150., 2600.);

	TH1D* mass_Ds2pipipi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2400., 2800.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(#pi^{+}_{K} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2 = new TH1D("    ", ";m(#pi^{+}_{K} #pi^{-}#pi^{-}_{#bar{p}} ) [MeV/c^{2}];Events", 50, 2400., 2800.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_afterVeto_2 = new TH1D("    ", ";m(#pi^{+}_{K} #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2100., 2600.);

	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
		pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);


		//misID Lorentz vectors
		piplus_asKaon_MisID.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massKaon);
		piminus_asProton_MisID.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massProton);
		piminus2_asProton_MisID.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massProton);


		//fill histos without veto cuts
		mass_Ds2pipipi_as_D2Kpipi_beforeVeto->Fill((piplus_asKaon_MisID + pi_minus_fromDs + pi_minus2_fromDs).M());

		mass_Ds2pipipi_as_D02pipi_beforeVeto->Fill((pi_plus_fromDs + pi_minus_fromDs).M());
		mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Fill((pi_plus_fromDs + pi_minus2_fromDs).M());

		mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Fill((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M());
		mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Fill((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M());

		mass_Ds2pipipi_as_Lc2KPpi_beforeVeto->Fill((piplus_asKaon_MisID + piminus_asProton_MisID + pi_minus2_fromDs).M());
		mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2->Fill((piplus_asKaon_MisID + pi_minus_fromDs + piminus2_asProton_MisID).M());


		//fill histos with veto cuts
		if((pi_plus_fromDs + pi_minus_fromDs).M() < 1700.){
		mass_Ds2pipipi_as_D02pipi_afterVeto->Fill((pi_plus_fromDs + pi_minus_fromDs).M());
		}
		if((pi_plus_fromDs + pi_minus2_fromDs).M() < 1700.){
		mass_Ds2pipipi_as_D02pipi_afterVeto_2->Fill((pi_plus_fromDs + pi_minus2_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M() - massLambda_c) > 30. || pi_minus_fromDs_PIDp < 0.){
		mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Fill((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M() - massLambda_c) > 30. || pi_minus2_fromDs_PIDp < 0.){
		mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Fill((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M());
		}

	}

	TCanvas* c = new TCanvas();

	mass_Ds2pipipi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D2Kpipi_beforeVeto.eps");

	mass_Ds2pipipi_as_D02pipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_beforeVeto.eps");
	mass_Ds2pipipi_as_D02pipi_afterVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_afterVeto.eps");

    	mass_Ds2pipipi_as_D02pipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2pipipi_as_D02pipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_D02pipi_beforeVeto->Draw("");
    	mass_Ds2pipipi_as_D02pipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.2,0.7,0.45,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2pipipi_as_D02pipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2pipipi_as_D02pipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_compareVeto.eps");


	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_beforeVeto_2.eps");
	mass_Ds2pipipi_as_D02pipi_afterVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_afterVeto_2.eps");

    	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->SetLineColor(kRed);
    	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Draw("");
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->Draw("same");

    	TLegend *leg_2 = new TLegend(0.2,0.7,0.45,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2pipipi_as_D02pipi_beforeVeto_2,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2pipipi_as_D02pipi_afterVeto_2,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_D02pipi_compareVeto_2.eps");



	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_beforeVeto.eps");
	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_afterVeto.eps");

    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Draw("");
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_compareVeto.eps");



	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_beforeVeto_2.eps");
	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_afterVeto_2.eps");

    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->SetLineColor(kRed);
    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Draw("");
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Draw("same");

    	TLegend *leg_4 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_4->SetHeader(" ");
    	leg_4->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2,"Data without veto","LEP");
    	leg_4->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2,"Data with veto","LEP");
    	leg_4->SetLineColor(kWhite);
    	leg_4->SetFillColor(kWhite);
    	leg_4->SetTextSize(0.05);
    	leg_4->Draw();
	c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2piPpi_compareVeto_2.eps");



	mass_Ds2pipipi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2KPpi_beforeVeto.eps");
	mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/norm_Ds2pipipi_as_Lc2KPpi_beforeVeto_2.eps");
}

void BkgStudies_signal_Ds2pipipi(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector pi_plus_fromDs;
	TLorentzVector pi_minus_fromDs;
	TLorentzVector pi_minus2_fromDs;

	//misIDs
	TLorentzVector piplus_asKaon_MisID;
	TLorentzVector piminus_asProton_MisID;
	TLorentzVector piminus2_asProton_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;

	//define observables
	Double_t pi_plus_fromDs_PX;
	Double_t pi_plus_fromDs_PY;
	Double_t pi_plus_fromDs_PZ;
	Double_t pi_minus_fromDs_PX;
	Double_t pi_minus_fromDs_PY;
	Double_t pi_minus_fromDs_PZ;
	Double_t pi_minus2_fromDs_PX;
	Double_t pi_minus2_fromDs_PY;
	Double_t pi_minus2_fromDs_PZ;

	Double_t pi_minus_fromDs_PIDK;
	Double_t pi_minus_fromDs_PIDp;
	Double_t pi_minus2_fromDs_PIDp;

	//link to tree
        tree->SetBranchAddress( "pi_plus_fromDs_PX" , &pi_plus_fromDs_PX );
        tree->SetBranchAddress( "pi_plus_fromDs_PY" , &pi_plus_fromDs_PY );
        tree->SetBranchAddress( "pi_plus_fromDs_PZ" , &pi_plus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus_fromDs_PX" , &pi_minus_fromDs_PX );
        tree->SetBranchAddress( "pi_minus_fromDs_PY" , &pi_minus_fromDs_PY );
        tree->SetBranchAddress( "pi_minus_fromDs_PZ" , &pi_minus_fromDs_PZ );
        tree->SetBranchAddress( "pi_minus2_fromDs_PX" , &pi_minus2_fromDs_PX );
        tree->SetBranchAddress( "pi_minus2_fromDs_PY" , &pi_minus2_fromDs_PY );
        tree->SetBranchAddress( "pi_minus2_fromDs_PZ" , &pi_minus2_fromDs_PZ );

        tree->SetBranchAddress( "pi_minus_fromDs_PIDK" , &pi_minus_fromDs_PIDK );
        tree->SetBranchAddress( "pi_minus_fromDs_PIDp" , &pi_minus_fromDs_PIDp );
        tree->SetBranchAddress( "pi_minus2_fromDs_PIDp" , &pi_minus2_fromDs_PIDp );


	//misID histos
	TH1D* mass_Ds2pipipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 1900., 3000.);

	TH1D* mass_Ds2pipipi_as_D02pipi_beforeVeto = new TH1D(" ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_afterVeto = new TH1D("   ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_beforeVeto_2 = new TH1D("     ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);
	TH1D* mass_Ds2pipipi_as_D02pipi_afterVeto_2 = new TH1D("        ", ";m(#pi^{+} #pi^{-}) [MeV/c^{2}];Events", 50, 1600., 1900.);

	TH1D* mass_Ds2pipipi_as_Lc2piPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_afterVeto = new TH1D("   ", ";m(#pi^{+} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2 = new TH1D("     ", ";m(#pi^{+}  #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2150., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2 = new TH1D("       ", ";m(#pi^{+} #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2150., 2600.);

	TH1D* mass_Ds2pipipi_as_Lc2KPpi_beforeVeto = new TH1D(" ", ";m(#pi^{+}_{K} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2400., 2800.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_afterVeto = new TH1D("   ", ";m(#pi^{+}_{K} #pi^{-}_{#bar{p}} #pi^{-}) [MeV/c^{2}];Events", 50, 2100., 2600.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2 = new TH1D("    ", ";m(#pi^{+}_{K} #pi^{-}#pi^{-}_{#bar{p}} ) [MeV/c^{2}];Events", 50, 2400., 2800.);
	TH1D* mass_Ds2pipipi_as_Lc2KPpi_afterVeto_2 = new TH1D("    ", ";m(#pi^{+}_{K} #pi^{-} #pi^{-}_{#bar{p}}) [MeV/c^{2}];Events", 50, 2100., 2600.);


	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
        	pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);
		pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
		pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);


		//misID Lorentz vectors
		piplus_asKaon_MisID.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massKaon);
		piminus_asProton_MisID.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massProton);
		piminus2_asProton_MisID.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massProton);

		//fill histos without veto cuts
		mass_Ds2pipipi_as_D2Kpipi_beforeVeto->Fill((piplus_asKaon_MisID + pi_minus_fromDs + pi_minus2_fromDs).M());
		mass_Ds2pipipi_as_D02pipi_beforeVeto->Fill((pi_plus_fromDs + pi_minus_fromDs).M());
		mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Fill((pi_plus_fromDs + pi_minus2_fromDs).M());

		mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Fill((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M());
		mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Fill((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M());

		mass_Ds2pipipi_as_Lc2KPpi_beforeVeto->Fill((piplus_asKaon_MisID + piminus_asProton_MisID + pi_minus2_fromDs).M());
		mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2->Fill((piplus_asKaon_MisID + pi_minus_fromDs + piminus2_asProton_MisID).M());


		//fill histos with veto cuts
		if((pi_plus_fromDs + pi_minus_fromDs).M() < 1700.){
		mass_Ds2pipipi_as_D02pipi_afterVeto->Fill((pi_plus_fromDs + pi_minus_fromDs).M());
		}
		if((pi_plus_fromDs + pi_minus2_fromDs).M() < 1700.){
		mass_Ds2pipipi_as_D02pipi_afterVeto_2->Fill((pi_plus_fromDs + pi_minus2_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M() - massLambda_c) > 30. || pi_minus_fromDs_PIDp < 0.){
		mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Fill((pi_plus_fromDs + piminus_asProton_MisID + pi_minus2_fromDs).M());
		}
		if(TMath::Abs((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M() - massLambda_c) > 30. || pi_minus2_fromDs_PIDp < 0.){
		mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Fill((pi_plus_fromDs + pi_minus_fromDs + piminus2_asProton_MisID).M());
		}
	}

	TCanvas* c = new TCanvas();

	mass_Ds2pipipi_as_D2Kpipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D2Kpipi_beforeVeto.eps");

	mass_Ds2pipipi_as_D02pipi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_beforeVeto.eps");
	mass_Ds2pipipi_as_D02pipi_afterVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_afterVeto.eps");

    	mass_Ds2pipipi_as_D02pipi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2pipipi_as_D02pipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_D02pipi_beforeVeto->Draw("");
    	mass_Ds2pipipi_as_D02pipi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.2,0.3,0.45,0.45);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Ds2pipipi_as_D02pipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Ds2pipipi_as_D02pipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_compareVeto.eps");



	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_beforeVeto_2.eps");
	mass_Ds2pipipi_as_D02pipi_afterVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_afterVeto_2.eps");

    	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->SetLineColor(kRed);
    	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_D02pipi_beforeVeto_2->Draw("");
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_D02pipi_afterVeto_2->Draw("same");

    	TLegend *leg_2 = new TLegend(0.2,0.3,0.45,0.45);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Ds2pipipi_as_D02pipi_beforeVeto_2,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Ds2pipipi_as_D02pipi_afterVeto_2,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_D02pipi_compareVeto_2.eps");



	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_beforeVeto.eps");
	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_afterVeto.eps");

    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->SetLineColor(kRed);
    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto->Draw("");
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_compareVeto.eps");



	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_beforeVeto_2.eps");
	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_afterVeto_2.eps");

    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->SetLineColor(kRed);
    	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2->Draw("");
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->SetLineColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2->Draw("same");

    	TLegend *leg_4 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_4->SetHeader(" ");
    	leg_4->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_beforeVeto_2,"Data without veto","LEP");
    	leg_4->AddEntry(mass_Ds2pipipi_as_Lc2piPpi_afterVeto_2,"Data with veto","LEP");
    	leg_4->SetLineColor(kWhite);
    	leg_4->SetFillColor(kWhite);
    	leg_4->SetTextSize(0.05);
    	leg_4->Draw();
	c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2piPpi_compareVeto_2.eps");


	mass_Ds2pipipi_as_Lc2KPpi_beforeVeto->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2KPpi_beforeVeto.eps");
	mass_Ds2pipipi_as_Lc2KPpi_beforeVeto_2->Draw("E1"); c->Print("eps/Ds2pipipi/signal_Ds2pipipi_as_Lc2KPpi_beforeVeto_2.eps");

}

void BkgStudies_Xd(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_16_noVetoes.root");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_16_noVetoes.root");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector Ds;
	TLorentzVector pi_plus;
	TLorentzVector pi_plus2;
	TLorentzVector pi_minus;

	//misIDs
	TLorentzVector piplus_asKaon_MisID;
	TLorentzVector piplus2_asKaon_MisID;

	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;


	//define observables
	Double_t Ds_PX;
	Double_t Ds_PY;
	Double_t Ds_PZ;

	Double_t pi_plus2_PX;
	Double_t pi_plus2_PY;
	Double_t pi_plus2_PZ;
	Double_t pi_minus_PX;
	Double_t pi_minus_PY;
	Double_t pi_minus_PZ;
	Double_t pi_plus1_PX;
	Double_t pi_plus1_PY;
	Double_t pi_plus1_PZ;

	Double_t pi_plus1_PIDK;
	Double_t pi_plus2_PIDK;


	//link to tree
        tree->SetBranchAddress( "Ds_PX" , &Ds_PX );
        tree->SetBranchAddress( "Ds_PY" , &Ds_PY );
        tree->SetBranchAddress( "Ds_PZ" , &Ds_PZ );

        tree->SetBranchAddress( "pi_plus1_PX" , &pi_plus1_PX );
        tree->SetBranchAddress( "pi_plus1_PY" , &pi_plus1_PY );
        tree->SetBranchAddress( "pi_plus1_PZ" , &pi_plus1_PZ );
        tree->SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
        tree->SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
        tree->SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
        tree->SetBranchAddress( "pi_plus2_PX" , &pi_plus2_PX );
        tree->SetBranchAddress( "pi_plus2_PY" , &pi_plus2_PY );
        tree->SetBranchAddress( "pi_plus2_PZ" , &pi_plus2_PZ );

        tree->SetBranchAddress( "pi_plus1_PIDK" , &pi_plus1_PIDK );
        tree->SetBranchAddress( "pi_plus2_PIDK" , &pi_plus2_PIDK );


	//misID histos
	TH1D* mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto = new TH1D(" ", ";m(D_{s}^{-}#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto = new TH1D("    ", ";m(D_{s}^{-}#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2 = new TH1D("     ", ";m(D_{s}^{-}#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2 = new TH1D("       ", ";m(D_{s}^{-}#pi^{+}_{K}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);

	TH1D* mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto = new TH1D(" ", ";m(B^{0}_{s}#rightarrow D^{-}_{s} (D_{s}^{+}#rightarrow#pi^{+}_{K}#pi^{-}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto = new TH1D("    ", ";m(B^{0}_{s}#rightarrow D^{-}_{s} (D_{s}^{+}#rightarrow#pi^{+}_{K}#pi^{-}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2 = new TH1D(" ", ";m(B^{0}_{s}#rightarrow D^{-}_{s} (D_{s}^{+}#rightarrow#pi^{+}_{K}#pi^{-}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2200.);
	TH1D* mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2 = new TH1D("    ", ";m(B^{0}_{s}#rightarrow D^{-}_{s} (D_{s}^{+}#rightarrow#pi^{+}_{K}#pi^{-}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2200.);

	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
		Ds.SetXYZM(Ds_PX, Ds_PY, Ds_PZ, massDs);

        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
		pi_plus.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
		pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);


		//misID Lorentz vectors
		piplus_asKaon_MisID.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massKaon);
		piplus2_asKaon_MisID.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massKaon);

		//fill histos without veto cuts
		mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto->Fill((Ds + piplus_asKaon_MisID + pi_minus + pi_plus2).M());
		mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2->Fill((Ds + piplus2_asKaon_MisID + pi_minus + pi_plus).M());

		mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto->Fill((piplus_asKaon_MisID + pi_minus + pi_plus2).M());
		mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2->Fill((piplus2_asKaon_MisID + pi_minus + pi_plus).M());


		//fill histo after cut
		if(pi_plus1_PIDK < 5.){
		mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto->Fill((Ds + piplus_asKaon_MisID + pi_minus + pi_plus2).M());
		}
		if(pi_plus2_PIDK < 5.){
		mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2->Fill((Ds + pi_plus + pi_minus + piplus2_asKaon_MisID).M());
		}
		if(TMath::Abs((piplus_asKaon_MisID + pi_minus + pi_plus2).M() - massDs) > 30. ||  pi_plus1_PIDK < 0.){
		mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto->Fill((piplus_asKaon_MisID + pi_minus + pi_plus2).M());
		}
		if(TMath::Abs((pi_plus + pi_minus + piplus2_asKaon_MisID).M() - massDs) > 30. ||  pi_plus2_PIDK < 0.){
		mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2->Fill((piplus2_asKaon_MisID + pi_minus + pi_plus).M());
		}
	}

	TCanvas* c = new TCanvas();

	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto.eps");
	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto->Draw("E1");c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto.eps");

    	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto->SetLineColor(kRed);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto->Draw("");
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto->SetLineColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_compareVeto.eps");



	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2.eps");
	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2->Draw("E1");c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2.eps");

    	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2->SetLineColor(kRed);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2->Draw("");
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2->SetLineColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2->Draw("same");

    	TLegend *leg_2 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Bs2Dspipipi_as_Bs2DsKpipi_beforeVeto_2,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Bs2Dspipipi_as_Bs2DsKpipi_afterVeto_2,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsKpipi_compareVeto_2.eps");



	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_beforeVeto.eps");
	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_afterVeto.eps");

    	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto->SetLineColor(kRed);
    	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto->SetMarkerColor(kRed);
   	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto->Draw("");
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto->SetLineColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto->SetMarkerColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto->Draw("same");

    	TLegend *leg_3 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_3->SetHeader(" ");
    	leg_3->AddEntry(mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto,"Data without veto","LEP");
    	leg_3->AddEntry(mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto,"Data with veto","LEP");
    	leg_3->SetLineColor(kWhite);
    	leg_3->SetFillColor(kWhite);
    	leg_3->SetTextSize(0.05);
    	leg_3->Draw();
	c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_compareVeto.eps");



	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2.eps");
	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2->Draw("E1"); c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2.eps");

    	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2->SetLineColor(kRed);
    	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2->SetMarkerColor(kRed);
   	mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2->Draw("");
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2->SetLineColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2->SetMarkerColor(kBlack);
    	mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2->Draw("same");

    	TLegend *leg_4 = new TLegend(0.6,0.7,0.85,0.85);
    	leg_4->SetHeader(" ");
    	leg_4->AddEntry(mass_Bs2Dspipipi_as_Bs2DsDs_beforeVeto_2,"Data without veto","LEP");
    	leg_4->AddEntry(mass_Bs2Dspipipi_as_Bs2DsDs_afterVeto_2,"Data with veto","LEP");
    	leg_4->SetLineColor(kWhite);
    	leg_4->SetFillColor(kWhite);
    	leg_4->SetTextSize(0.05);
    	leg_4->Draw();
	c->Print("eps/Xd/norm_Bs2Dspipipi_as_Bs2DsDs_compareVeto_2.eps");
}


void BkgStudies_Xs(){

    	TChain* tree = 0;
   	tree =new TChain("DecayTree");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_noVetoes.root");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_16_noVetoes.root");

    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_11_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_12_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_15_noVetoes.root");
    	tree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_16_noVetoes.root");

   	tree->SetBranchStatus("*",0);
    	tree->SetBranchStatus("*M",1) ;
    	tree->SetBranchStatus("*MM",1) ;
    	tree->SetBranchStatus( "*P", 1 );
   	tree->SetBranchStatus( "*PX", 1 );
    	tree->SetBranchStatus( "*PY", 1);
   	tree->SetBranchStatus( "*PZ", 1);
    	tree->SetBranchStatus( "*PE", 1);
    	tree->SetBranchStatus( "*PT", 1 );
	tree->SetBranchStatus( "*PIDK*", 1);
	tree->SetBranchStatus( "*PIDp*", 1);

	//define Lorentz vectors
	TLorentzVector Ds;
	TLorentzVector pi_plus;
	TLorentzVector K_plus;
	TLorentzVector pi_minus;

	//misIDs
	TLorentzVector Kplus_asPion_MisID;
	TLorentzVector piminus_asKaon_MisID;
	//masses
	double massKaon = 493.68;
	double massPion = 139.57;
	double massProton = 938.27;
	double massPhi = 1019.46;
	double massKstar = 895.81;
	double massDs = 1968.30;
	double massDminus = 1869.61;
	double massLambda_c = 2286.46;
	double massD0 = 1864.83;


	//define observables
	Double_t Ds_PX;
	Double_t Ds_PY;
	Double_t Ds_PZ;

	Double_t K_plus_PX;
	Double_t K_plus_PY;
	Double_t K_plus_PZ;
	Double_t pi_minus_PX;
	Double_t pi_minus_PY;
	Double_t pi_minus_PZ;
	Double_t pi_plus_PX;
	Double_t pi_plus_PY;
	Double_t pi_plus_PZ;

	Double_t K_plus_PIDK;
	Double_t pi_minus_PIDK;

	//link to tree
        tree->SetBranchAddress( "Ds_PX" , &Ds_PX );
        tree->SetBranchAddress( "Ds_PY" , &Ds_PY );
        tree->SetBranchAddress( "Ds_PZ" , &Ds_PZ );

        tree->SetBranchAddress( "pi_plus_PX" , &pi_plus_PX );
        tree->SetBranchAddress( "pi_plus_PY" , &pi_plus_PY );
        tree->SetBranchAddress( "pi_plus_PZ" , &pi_plus_PZ );
        tree->SetBranchAddress( "pi_minus_PX" , &pi_minus_PX );
        tree->SetBranchAddress( "pi_minus_PY" , &pi_minus_PY );
        tree->SetBranchAddress( "pi_minus_PZ" , &pi_minus_PZ );
        tree->SetBranchAddress( "K_plus_PX" , &K_plus_PX );
        tree->SetBranchAddress( "K_plus_PY" , &K_plus_PY );
        tree->SetBranchAddress( "K_plus_PZ" , &K_plus_PZ );

        tree->SetBranchAddress( "K_plus_PIDK" , &K_plus_PIDK );
        tree->SetBranchAddress( "pi_minus_PIDK" , &pi_minus_PIDK );


	//misID histos
	TH1D* mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto = new TH1D(" ", ";m(D_{s}^{-}K^{+}_{#pi}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);
	TH1D* mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto = new TH1D("   ", ";m(D_{s}^{-}K^{+}_{#pi}#pi^{-}#pi^{-}) [MeV/c^{2}];Events", 50, 5000., 6000.);

	TH1D* mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto = new TH1D(" ", ";m(B^{0}_{s}#rightarrow D^{-}_{s}(D_{s}^{+}#rightarrow K^{+}#pi^{-}_{K}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2100.);
	TH1D* mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto = new TH1D("   ", ";m(B^{0}_{s}#rightarrow D^{-}_{s}(D_{s}^{+}#rightarrow K^{+}#pi^{-}_{K}#pi^{+})) [MeV/c^{2}];Events", 50, 1600., 2100.);

	//loop over events
	int numEvents = tree->GetEntries();
	for(int i=0; i< numEvents; i++)
        	{
        	if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        	tree->GetEntry(i);

        	//fill the Lorentz vectors
		Ds.SetXYZM(Ds_PX, Ds_PY, Ds_PZ, massDs);

        	pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
		pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);
		K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);


		//misID Lorentz vectors
		Kplus_asPion_MisID.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massPion);
		piminus_asKaon_MisID.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massKaon);

		//fill histos without veto cuts
		mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto->Fill((Ds + Kplus_asPion_MisID + pi_minus + pi_plus).M());
		mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto->Fill((K_plus + piminus_asKaon_MisID + pi_plus).M());


		//fill histos with vetoes applied
		if(K_plus_PIDK > 10.){
		mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto->Fill((Ds + Kplus_asPion_MisID + pi_minus + pi_plus).M());
		}
		if(TMath::Abs((K_plus + piminus_asKaon_MisID + pi_plus).M() - massDs) > 30. || pi_minus_PIDK < 0.){
		mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto->Fill((K_plus + piminus_asKaon_MisID + pi_plus).M());
		}

	}

	TCanvas* c = new TCanvas();
	mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto->Draw("E1"); c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto.eps");
	mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto->Draw("E1"); c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto.eps");

   	mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto->SetLineColor(kRed);
   	mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto->SetMarkerColor(kRed);
   	mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto->Draw("");
    	mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto->SetLineColor(kBlack);
    	mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto->SetMarkerColor(kBlack);
    	mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto->Draw("same");

    	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
    	leg->SetHeader(" ");
    	leg->AddEntry(mass_Bs2DsKpipi_as_Bs2Dspipipi_beforeVeto,"Data without veto","LEP");
    	leg->AddEntry(mass_Bs2DsKpipi_as_Bs2Dspipipi_afterVeto,"Data with veto","LEP");
    	leg->SetLineColor(kWhite);
    	leg->SetFillColor(kWhite);
    	leg->SetTextSize(0.05);
    	leg->Draw();
	c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2Dspipipi_compareVeto.eps");




	mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto->Draw("E1"); c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2DsDs_beforeVeto.eps");
	mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto->Draw("E1"); c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2DsDs_afterVeto.eps");

   	mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto->SetLineColor(kRed);
   	mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto->SetMarkerColor(kRed);
   	mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto->Draw("");
    	mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto->SetLineColor(kBlack);
    	mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto->SetMarkerColor(kBlack);
    	mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto->Draw("same");

    	TLegend *leg_2 = new TLegend(0.2,0.3,0.45,0.45);
    	leg_2->SetHeader(" ");
    	leg_2->AddEntry(mass_Bs2DsKpipi_as_Bs2DsDs_beforeVeto,"Data without veto","LEP");
    	leg_2->AddEntry(mass_Bs2DsKpipi_as_Bs2DsDs_afterVeto,"Data with veto","LEP");
    	leg_2->SetLineColor(kWhite);
    	leg_2->SetFillColor(kWhite);
    	leg_2->SetTextSize(0.05);
    	leg_2->Draw();
	c->Print("eps/Xs/signal_Bs2DsKpipi_as_Bs2DsDs_compareVeto.eps");
}


int main(int argc, char** argv){

	time_t startTime = time(0);
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gROOT->ProcessLine(".x ../lhcbStyle.C");

	NamedParameter<string> Channel("channel", (std::string) "Norm" , (char*) 0);
	string channel = Channel;
	NamedParameter<int>  do_Ds2pipipi("do_Ds2pipipi", 1);
	NamedParameter<int>  do_Ds2Kpipi("do_Ds2Kpipi", 1);
	NamedParameter<int>  do_Ds2KKpi("do_Ds2KKpi", 1);



	if(channel == "Norm") {
		BkgStudies_Xd();
		if(do_Ds2pipipi == 1) BkgStudies_norm_Ds2pipipi();
		if(do_Ds2Kpipi == 1) BkgStudies_norm_Ds2Kpipi();
		if(do_Ds2KKpi == 1) BkgStudies_norm_Ds2KKpi();
	}
  	else if(channel == "Signal"){
		BkgStudies_Xs();
		if(do_Ds2pipipi == 1) BkgStudies_signal_Ds2pipipi();
		if(do_Ds2Kpipi == 1) BkgStudies_signal_Ds2Kpipi();
		if(do_Ds2KKpi == 1) BkgStudies_signal_Ds2KKpi();		
	}

    	else{
        	std::cout << "*********************************************************" << std::endl;
        	std::cout << "please specify 'Signal' or 'Norm' in options file" << std::endl;
        	std::cout << "*********************************************************" << std::endl;
    	}


	std::cout << "==============================================" << std::endl;
	std::cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << std::endl;
	std::cout << "==============================================" << std::endl;

	return 0;

}