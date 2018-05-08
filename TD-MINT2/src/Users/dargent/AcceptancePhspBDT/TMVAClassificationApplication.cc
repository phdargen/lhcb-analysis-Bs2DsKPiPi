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

void TMVAClassificationApplication(TString myMethod = "BDTG", TString applyTo = "Gen") 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------
   TChain* theTree;
   TString outFileName;

   if(applyTo == "Gen"){
   	theTree = new TChain("MCDecayTree");
   	outFileName = "GenMC_BDT.root";
   	theTree->Add("GenMC.root");
   }
   else if(applyTo == "Sel"){
   	theTree = new TChain("DecayTree");
   	outFileName = "SelMC_BDT.root";
   	theTree->Add("SelMC.root");
   }
   else if(applyTo == "MINT"){
   	theTree = new TChain("DalitzEventList");
   	outFileName = "MintMC_BDT.root";
   	theTree->Add("MintMC.root");
   }

   // Disable not needed branches
   // Ouput tree
   theTree->SetBranchStatus("*",0);
   theTree->SetBranchStatus("s_*",1);
   theTree->SetBranchStatus("cos_*",1);
   theTree->SetBranchStatus("phi_*",1);
   if(applyTo == "Gen")theTree->SetBranchStatus("*TRUEP*",1);
   if(applyTo == "MINT")theTree->SetBranchStatus("*",1);

   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   Float_t r_s_Kpipi;
   Float_t r_s_Kpi;
   Float_t r_s_pipi;
   Float_t r_s_Dspipi;
   Float_t r_s_Dspi;
   Float_t r_s_DsK;
   Float_t r_cos_theta_Kpi;
   Float_t r_cos_theta_Dspi;
   Float_t r_phi_Kpi_Dspi;
   Float_t r_cos_theta_pipi;
   Float_t r_cos_theta_DsK;
   Float_t r_theta_DsK;
   Float_t r_phi_pipi_DsK;
 
   //reader->AddVariable( "s_Kpipi", &r_s_Kpipi );
   reader->AddVariable( "s_Kpi", &r_s_Kpi );
   //reader->AddVariable( "s_pipi", &r_s_pipi );
   reader->AddVariable( "s_Dspi", &r_s_Dspi );
   //reader->AddVariable( "s_Dspipi", &r_s_Dspipi );
   //reader->AddVariable( "s_DsK", &r_s_DsK );   
   reader->AddVariable( "cos_theta_Kpi", &r_cos_theta_Kpi );
   reader->AddVariable( "cos_theta_Dspi", &r_cos_theta_Dspi );
   reader->AddVariable( "phi_Kpi_Dspi", &r_phi_Kpi_Dspi );
   //reader->AddVariable( "cos_theta_pipi", &r_cos_theta_pipi );
   //reader->AddVariable( "theta_DsK:=acos(cos_theta_DsK)", &r_theta_DsK );
   //reader->AddVariable( "phi_pipi_DsK", &r_phi_pipi_DsK );

   // --- Book the MVA methods
   TString prefix = "weights/TMVAClassification_run1_t0";
   reader->BookMVA( myMethod, prefix + "_" + myMethod + ".weights.xml" ); 

   Double_t s_Kpipi;
   Double_t s_Kpi;
   Double_t s_pipi;
   Double_t s_Dspipi;
   Double_t s_Dspi;
   Double_t s_DsK;
   Double_t cos_theta_Kpi;
   Double_t cos_theta_Dspi;
   Double_t phi_Kpi_Dspi;
   Double_t cos_theta_pipi;
   Double_t cos_theta_DsK;
   Double_t phi_pipi_DsK;

   theTree->SetBranchAddress( "s_Kpipi", &s_Kpipi );
   theTree->SetBranchAddress( "s_Kpi", &s_Kpi );
   theTree->SetBranchAddress( "s_pipi", &s_pipi );
   theTree->SetBranchAddress( "s_Dspipi", &s_Dspipi );
   theTree->SetBranchAddress( "s_Dspi", &s_Dspi );
   theTree->SetBranchAddress( "s_DsK", &s_DsK );
   theTree->SetBranchAddress( "cos_theta_Kpi", &cos_theta_Kpi );
   theTree->SetBranchAddress( "cos_theta_Dspi", &cos_theta_Dspi );
   theTree->SetBranchAddress( "phi_Kpi_Dspi", &phi_Kpi_Dspi );
   //theTree->SetBranchAddress( "cos_theta_pipi", &cos_theta_pipi );
   //theTree->SetBranchAddress( "cos_theta_DsK", &cos_theta_DsK );
   //theTree->SetBranchAddress( "phi_pipi_DsK", &phi_pipi_DsK );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   double BDTG;
   tree->Branch("BDTG",&BDTG, "BDTG/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);

	r_s_Kpipi= float(s_Kpipi);
	r_s_Kpi= float(s_Kpi);
	r_s_pipi= float(s_pipi);
	r_s_Dspipi= float(s_Dspipi);
	r_s_Dspi= float(s_Dspi);
	r_s_DsK= float(s_DsK);

	r_cos_theta_Kpi= float(cos_theta_Kpi);
	r_cos_theta_Dspi= float(cos_theta_Dspi);
	r_phi_Kpi_Dspi= float(phi_Kpi_Dspi);

	//r_cos_theta_pipi= float(cos_theta_pipi);
	//r_cos_theta_DsK= float(cos_theta_DsK);
	//r_theta_DsK= float(acos(cos_theta_DsK));
	//r_phi_pipi_DsK= float(phi_pipi_DsK);

	TString methodName = myMethod;
	BDTG = double(reader->EvaluateMVA(methodName));
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

void applyToAll(TString myMethod = "BDTG"){
	TMVAClassificationApplication(myMethod, "Gen" );
	TMVAClassificationApplication(myMethod, "Sel" );
	TMVAClassificationApplication(myMethod, "MINT" );
}