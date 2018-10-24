/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <TChain.h>
#include "TStopwatch.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;

void TMVAClassificationApplication(TString decay = "Signal", TString dataType = "Data", TString myMethod = "BDTG", TString trainedOn = "Data" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------
   TChain* theTree = new TChain("DecayTree");


   TString outFileName = "/auto/data/dargent/BsDsKpipi/BDT/";

   if(decay == "Signal" && dataType == "Data"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_17.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_16.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_17.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_16.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2Kpipi_17.root");	

	outFileName += "Data/signal.root";
   }

   else if(decay == "Signal" && dataType == "Data_SS"){ 	  
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12_SS.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_SS.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12_SS.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_16_SS.root");
	
	outFileName += "Data/signal_SS.root";
   }

   else if(decay == "Signal" && dataType == "MC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_16.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_16.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_16.root");

	outFileName += "MC/signal.root";
   }
   else if(decay == "Signal" && dataType == "MC_PIDGen"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12_PIDGen.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_16_PIDGen.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_12_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_16_PIDGen.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_12_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_16_PIDGen.root");

	outFileName += "MC/signal_PIDGen.root";
   }
   else if(decay == "Signal" && dataType == "MC_PIDMC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12_PIDMC.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_16_PIDMC.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_12_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2pipipi_16_PIDMC.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_12_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2Kpipi_16_PIDMC.root");

	outFileName += "MC/signal_PIDMC.root";
   }


   else if(decay == "Norm" && dataType == "Data"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_16.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_17.root");
   	
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_16.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_17.root");

	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_16.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2Kpipi_17.root");

	outFileName += "Data/norm.root";
   }

   else if(decay == "Norm" && dataType == "MC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_16.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_16.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_16.root");

	outFileName += "MC/norm.root";
   }
   else if(decay == "Norm" && dataType == "MC_PIDGen"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12_PIDGen.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_16_PIDGen.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_12_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_16_PIDGen.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_11_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_12_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_15_PIDGen.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_16_PIDGen.root");

	outFileName += "MC/norm_PIDGen.root";
   }
   else if(decay == "Norm" && dataType == "MC_PIDMC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12_PIDMC.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_16_PIDMC.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_12_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_16_PIDMC.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_11_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_12_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_15_PIDMC.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2Kpipi_16_PIDMC.root");

	outFileName += "MC/norm_PIDMC.root";
   }

   else {
	cout << "Unknown options, I'll crash now." << endl;
	throw "ERROR";
   }

   // Disable not needed branches
   //theTree->SetBranchStatus("*PARTICLES*",0); 
   theTree->SetBranchStatus("*MC12Tune*",0);
   //if(decay == "Norm")theTree->SetBranchStatus("*a_1_1260*",0);
   //else theTree->SetBranchStatus("*K_1_1270*",0);
   //theTree->SetBranchStatus("*SS_Proton*",0);
   //theTree->SetBranchStatus("*SS_Pion*",0);
   /*
   theTree->SetBranchStatus("*",0); 
   theTree->SetBranchStatus("*CHI2*",1); 
   theTree->SetBranchStatus("*DOCA*",1);
   theTree->SetBranchStatus("*DIRA*",1);
   theTree->SetBranchStatus("*PT*",1);
   theTree->SetBranchStatus("*RFD*",1);
   theTree->SetBranchStatus("*max*",1);
   theTree->SetBranchStatus("*ptasy*",1);
   theTree->SetBranchStatus("*MM*",1);
   theTree->SetBranchStatus("*Trigger*",1);
   theTree->SetBranchStatus("*State*",1);
   theTree->SetBranchStatus("year",1);
   theTree->SetBranchStatus("run",1);
   theTree->SetBranchStatus("*PIDK",1);
   theTree->SetBranchStatus("pi_plus_isMuon",1);
   theTree->SetBranchStatus("weight",1);
   theTree->SetBranchStatus("eventNumber",1);
   */
   // Ouput tree
   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   Float_t r_DTF_CHI2NDOF;
   Float_t r_log_Bs_IPCHI2_OWNPV;
   Float_t r_log_Bs_DIRA;
   Float_t r_log_XsDaughters_min_IPCHI2;
   Float_t r_Xs_ptasy;
   Float_t r_Xs_max_DOCA;
   Float_t r_log_DsDaughters_min_IPCHI2;
   Float_t r_Ds_ptasy;
   Float_t r_log_Ds_FDCHI2_ORIVX;
   Float_t r_log_Ds_RFD;
   Float_t r_maxCos;
   Float_t r_max_ghostProb;
   Float_t r_Bs_ptasy;
   Float_t r_log_Bs_SmallestDeltaChi2OneTrack;

   reader->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)",&r_log_Bs_IPCHI2_OWNPV );
   reader->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)",&r_log_Bs_DIRA );
   reader->AddVariable( "PV_CHI2NDOF", &r_DTF_CHI2NDOF );
   reader->AddVariable( "log_Bs_SmallestDeltaChi2OneTrack:= log(Bs_SmallestDeltaChi2OneTrack)",&r_log_Bs_SmallestDeltaChi2OneTrack );
   reader->AddVariable( "Bs_ptasy_1.00",&r_Bs_ptasy );
   reader->AddVariable("max_ghostProb",&r_max_ghostProb);

   reader->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)",&r_log_XsDaughters_min_IPCHI2 );
   //reader->AddVariable( "Xs_ptasy_1.00",&r_Xs_ptasy );
   reader->AddVariable( "Xs_max_DOCA",&r_Xs_max_DOCA);

   reader->AddVariable( "maxCos", &r_maxCos );

   reader->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)",&r_log_DsDaughters_min_IPCHI2);
   //reader->AddVariable( "Ds_ptasy_1.00",&r_Ds_ptasy);
   reader->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)",&r_log_Ds_FDCHI2_ORIVX);
   reader->AddVariable( "log_Ds_RFD:=log(Ds_RFD)",&r_log_Ds_RFD);


   // --- Book the MVA methods
   TString prefix = "weights/TMVAClassification_"+trainedOn+ "_";

   std::vector<TString> weightFiles;
   weightFiles.push_back("run1_t0_odd");
   weightFiles.push_back("run1_t0_even");
   weightFiles.push_back("run1_t1_odd");
   weightFiles.push_back("run1_t1_even");
   weightFiles.push_back("run2_t0_odd");
   weightFiles.push_back("run2_t0_even");
   weightFiles.push_back("run2_t1_odd");
   weightFiles.push_back("run2_t1_even");
 
   for(int i= 0 ; i < weightFiles.size(); i++) 
	reader->BookMVA( myMethod + weightFiles[i], prefix + weightFiles[i] + "_" + myMethod + ".weights.xml" ); 

   Double_t DTF_CHI2NDOF;
   Double_t Bs_IPCHI2_OWNPV;
   Double_t Bs_DIRA_OWNPV;
   Double_t XsDaughters_min_IPCHI2;
   Double_t Xs_ptasy;
   Double_t Xs_max_DOCA;
   Double_t DsDaughters_min_IPCHI2;
   Double_t Ds_ptasy;
   Double_t Ds_FDCHI2_ORIVX;
   Double_t Ds_RFD;
   Double_t maxCos;
   Double_t max_ghostProb;
   Double_t Bs_ptasy;
   Double_t Bs_SmallestDeltaChi2OneTrack;

   theTree->SetBranchAddress( "PV_CHI2NDOF", &DTF_CHI2NDOF );
   theTree->SetBranchAddress( "Bs_IPCHI2_OWNPV", &Bs_IPCHI2_OWNPV );
   theTree->SetBranchAddress( "Bs_DIRA_OWNPV", &Bs_DIRA_OWNPV );
   theTree->SetBranchAddress( "XsDaughters_min_IPCHI2", &XsDaughters_min_IPCHI2 );
   theTree->SetBranchAddress( "Xs_ptasy_1.00", &Xs_ptasy );
   theTree->SetBranchAddress( "Xs_max_DOCA", &Xs_max_DOCA );
   theTree->SetBranchAddress( "DsDaughters_min_IPCHI2", &DsDaughters_min_IPCHI2 );
   theTree->SetBranchAddress( "Ds_ptasy_1.00", &Ds_ptasy );
   theTree->SetBranchAddress( "Ds_FDCHI2_ORIVX", &Ds_FDCHI2_ORIVX );
   theTree->SetBranchAddress( "Ds_RFD", &Ds_RFD );
   theTree->SetBranchAddress( "maxCos", &maxCos );
   theTree->SetBranchAddress( "max_ghostProb", &max_ghostProb );
   theTree->SetBranchAddress( "Bs_ptasy_1.00", &Bs_ptasy );
   theTree->SetBranchAddress( "Bs_SmallestDeltaChi2OneTrack", &Bs_SmallestDeltaChi2OneTrack );
   
   Int_t year, run, Ds_finalState, TriggerCat; 
   ULong64_t eventNumber;

   theTree->SetBranchAddress( "year", &year );
   theTree->SetBranchAddress( "run", &run );
   theTree->SetBranchAddress( "Ds_finalState", &Ds_finalState );
   theTree->SetBranchAddress( "TriggerCat", &TriggerCat );
   theTree->SetBranchAddress( "eventNumber", &eventNumber );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   Float_t BDTG_response;
   double BDTG;
   tree->Branch("BDTG_response",&BDTG_response, "BDTG_response/F");
   tree->Branch("BDTG",&BDTG, "BDTG/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);

	if(Bs_SmallestDeltaChi2OneTrack<=5)continue;
	if(maxCos <= -0.95) continue;
	if(Bs_IPCHI2_OWNPV >= 16) continue;
	if(DTF_CHI2NDOF >= 15)continue;

	r_DTF_CHI2NDOF= float(DTF_CHI2NDOF);
        r_log_Bs_IPCHI2_OWNPV = float(log(Bs_IPCHI2_OWNPV));
        r_log_Bs_DIRA = float(log(1.-Bs_DIRA_OWNPV));
        r_log_XsDaughters_min_IPCHI2 = float(log(XsDaughters_min_IPCHI2));
        //r_Xs_ptasy = float(Xs_ptasy);
        r_Xs_max_DOCA = float(Xs_max_DOCA);
        r_log_DsDaughters_min_IPCHI2 = float(log(DsDaughters_min_IPCHI2));
        //r_Ds_ptasy = float(Ds_ptasy);
        r_log_Ds_FDCHI2_ORIVX = float(log(Ds_FDCHI2_ORIVX));
        r_log_Ds_RFD = float(log(Ds_RFD));
   	r_maxCos = float(maxCos);
        r_max_ghostProb = float(max_ghostProb);
        r_Bs_ptasy = float(Bs_ptasy);
	r_log_Bs_SmallestDeltaChi2OneTrack = float(log(Bs_SmallestDeltaChi2OneTrack));

	TString methodName = myMethod + "run";
	methodName += run;
	methodName += "_t"; 
	methodName += TriggerCat;
	// apply BDT trained on even sample to odd sample and viceversa
	if(eventNumber % 2 == 0) methodName += "_odd" ;
	else methodName += "_even";

        BDTG_response=reader->EvaluateMVA(methodName);
	BDTG = double(BDTG_response);
        tree->Fill();    
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   tree->Write();
   delete reader;
   hFile->Close();
    
   cout << "Wrote to file: " << outFileName << endl;
   cout << "==> TMVAClassificationApplication is done!" << endl << endl;
} 

void applyToAll(TString myMethod = "BDTG", TString trainedOn = "Data" ){

//  	TMVAClassificationApplication("Signal", "Data", myMethod, trainedOn );
// 	TMVAClassificationApplication("Norm", "Data", myMethod, trainedOn );

	TMVAClassificationApplication("Signal", "MC", myMethod, trainedOn );
	TMVAClassificationApplication("Signal", "MC_PIDGen", myMethod, trainedOn );
	TMVAClassificationApplication("Signal", "MC_PIDMC", myMethod, trainedOn );

	TMVAClassificationApplication("Norm", "MC", myMethod, trainedOn );
	TMVAClassificationApplication("Norm", "MC_PIDGen", myMethod, trainedOn );
	TMVAClassificationApplication("Norm", "MC_PIDMC", myMethod, trainedOn );
}