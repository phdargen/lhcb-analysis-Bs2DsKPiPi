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

void TMVAClassificationApplication(TString decay = "Signal", TString dataType = "Data", TString myMethodList = "BDTG" ) 
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
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_up.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_down.root");

   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_16.root");
	
	outFileName += "Data/signal.root";
   }

   else if(decay == "Signal" && dataType == "MC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12.root");
   	//theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_11.root");
   	//theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2pipipi_12.root");
	outFileName += "MC/signal.root";
   }

   else if(decay == "Norm" && dataType == "Data"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_12.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2KKpi_16.root");
   	
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_12.root");
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_15.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2pipipi_16.root");

	outFileName += "Data/norm.root";
   }

   else if(decay == "Norm" && dataType == "MC"){ 	  
	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_11.root");
   	theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2KKpi_12.root");
   	//theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_11.root");
   	//theTree->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/norm_Ds2pipipi_12.root");

	outFileName += "MC/norm.root";
   }

   else {
	cout << "Unknown options, I'll crash now." << endl;
	throw "ERROR";
   }

   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   Float_t r_DTF_CHI2NDOF;
   Float_t r_log_Bs_IPCHI2_OWNPV;
   Float_t r_log_Bs_DIRA;
   Float_t r_log_XsDaughters_min_IPCHI2;
   Float_t r_K_1_1270_plus_ptasy;
   Float_t r_Xs_max_DOCA;
   Float_t r_log_DsDaughters_min_IPCHI2;
   Float_t r_Ds_ptasy;
   Float_t r_log_Ds_FDCHI2_ORIVX;
   Float_t r_log_Ds_RFD;
   Float_t r_maxCos;
   Float_t r_max_ghostProb;

   reader->AddVariable( "DTF_CHI2NDOF", &r_DTF_CHI2NDOF );
   reader->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)",&r_log_Bs_IPCHI2_OWNPV );
   reader->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)",&r_log_Bs_DIRA );
  
   reader->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)",&r_log_XsDaughters_min_IPCHI2 );
   reader->AddVariable( "K_1_1270_plus_ptasy_1.00",&r_K_1_1270_plus_ptasy );
   reader->AddVariable( "Xs_max_DOCA",&r_Xs_max_DOCA);

   reader->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)",&r_log_DsDaughters_min_IPCHI2);
   reader->AddVariable( "Ds_ptasy_1.00",&r_Ds_ptasy);
   reader->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)",&r_log_Ds_FDCHI2_ORIVX);
   reader->AddVariable( "log_Ds_RFD:=log(Ds_RFD)",&r_log_Ds_RFD);
   reader->AddVariable( "maxCos", &r_maxCos );

   reader->AddVariable("max_ghostProb",&r_max_ghostProb);

   // --- Book the MVA methods
   TString dir    = "weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   Double_t DTF_CHI2NDOF;
   Double_t Bs_IPCHI2_OWNPV;
   Double_t Bs_DIRA_OWNPV;
   Double_t XsDaughters_min_IPCHI2;
   Double_t K_1_1270_plus_ptasy;
   Double_t Xs_max_DOCA;
   Double_t DsDaughters_min_IPCHI2;
   Double_t Ds_ptasy;
   Double_t Ds_FDCHI2_ORIVX;
   Double_t Ds_RFD;
   Double_t maxCos;
   Double_t max_ghostProb;

   theTree->SetBranchAddress( "DTF_CHI2NDOF", &DTF_CHI2NDOF );
   theTree->SetBranchAddress( "Bs_IPCHI2_OWNPV", &Bs_IPCHI2_OWNPV );
   theTree->SetBranchAddress( "Bs_DIRA_OWNPV", &Bs_DIRA_OWNPV );
   theTree->SetBranchAddress( "XsDaughters_min_IPCHI2", &XsDaughters_min_IPCHI2 );
   if(decay == "Signal")theTree->SetBranchAddress( "K_1_1270_plus_ptasy_1.00", &K_1_1270_plus_ptasy );
   else theTree->SetBranchAddress( "a_1_1260_plus_ptasy_1.00", &K_1_1270_plus_ptasy );
   theTree->SetBranchAddress( "Xs_max_DOCA", &Xs_max_DOCA );
   theTree->SetBranchAddress( "DsDaughters_min_IPCHI2", &DsDaughters_min_IPCHI2 );
   theTree->SetBranchAddress( "Ds_ptasy_1.00", &Ds_ptasy );
   theTree->SetBranchAddress( "Ds_FDCHI2_ORIVX", &Ds_FDCHI2_ORIVX );
   theTree->SetBranchAddress( "Ds_RFD", &Ds_RFD );
   theTree->SetBranchAddress( "maxCos", &maxCos );
   theTree->SetBranchAddress( "max_ghostProb", &max_ghostProb );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   Float_t BDTG_response;
   tree->Branch("BDTG_response",&BDTG_response, "BDTG_response/F");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);

	r_DTF_CHI2NDOF= float(DTF_CHI2NDOF);
        r_log_Bs_IPCHI2_OWNPV = float(log(Bs_IPCHI2_OWNPV));
        r_log_Bs_DIRA = float(log(1.-Bs_DIRA_OWNPV));
        r_log_XsDaughters_min_IPCHI2 = float(log(XsDaughters_min_IPCHI2));
        r_K_1_1270_plus_ptasy = float(K_1_1270_plus_ptasy);
        r_Xs_max_DOCA = float(Xs_max_DOCA);
        r_log_DsDaughters_min_IPCHI2 = float(log(DsDaughters_min_IPCHI2));
        r_Ds_ptasy = float(Ds_ptasy);
        r_log_Ds_FDCHI2_ORIVX = float(log(Ds_FDCHI2_ORIVX));
        r_log_Ds_RFD = float(log(Ds_RFD));
   	r_maxCos = float(maxCos);
        r_max_ghostProb = float(max_ghostProb);

        BDTG_response=reader->EvaluateMVA("BDTG method");
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

int main(int argc, char** argv){
	TMVAClassificationApplication(TString((string)argv[1]),TString((string)argv[2]));
	return 0;
}
