#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "BDTG", TString run = "Run1" )
{
   TChain* background = new TChain("DecayTree");
   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_11.root");
   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12.root");

   TChain* signal = new TChain("DecayTree");
   signal->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11.root");
   signal->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12.root");

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA_Bs2DsKpipi.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


   signal->SetBranchStatus("*",0);  // disable all branches
   signal->SetBranchStatus("*CHI2*",1); 
   signal->SetBranchStatus("*DOCA*",1);
   signal->SetBranchStatus("*DIRA*",1);
   signal->SetBranchStatus("*PT*",1);
   signal->SetBranchStatus("*RFD*",1);
   signal->SetBranchStatus("*max*",1);
   signal->SetBranchStatus("*ptasy*",1);
   signal->SetBranchStatus("*MM*",1);
   signal->SetBranchStatus("*TAU*",1);

   background->SetBranchStatus("*",0);  // disable all branches
   background->SetBranchStatus("*CHI2*",1); 
   background->SetBranchStatus("*DOCA*",1);
   background->SetBranchStatus("*DIRA*",1);
   background->SetBranchStatus("*PT*",1);
   background->SetBranchStatus("*RFD*",1);
   background->SetBranchStatus("*max*",1);
   background->SetBranchStatus("*ptasy*",1);
   background->SetBranchStatus("*MM*",1);

   // Define the input variables that shall be used for the MVA training
   factory->AddVariable( "DTF_CHI2NDOF", "DTF Fit chi2", "", 'F' );
   factory->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)","B_{s} ln(IP #chi^{2})", "", 'F' );
   factory->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)","ln(1 - B_{s} DIRA)","", 'F' );
  
   factory->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)","X_{s} daughters min ln(IP#chi^{2})", "", 'F' );
   factory->AddVariable( "K_1_1270_plus_ptasy_1.00","K1^{+}  cone p_{t} asy","", 'F' );
   factory->AddVariable( "Xs_max_DOCA","X_{s} max DOCA", "mm", 'F' );

   factory->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)","D_{s} daughters min ln(IP#chi^{2})", "", 'F' );
   factory->AddVariable( "Ds_ptasy_1.00","Ds  cone p_{t} asy","", 'F' );
   factory->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)","D_{s} FD significance", "", 'F' );
   factory->AddVariable( "log_Ds_RFD:=log(Ds_RFD)","D_{s} RFD", "", 'F' );
   factory->AddVariable( "maxCos", "cos(max[#theta_{Ds h}])", "", 'F' );

   factory->AddVariable("max_ghostProb","max(ghostProb)","",'F');

   // Additional variables for testing
   //factory->AddVariable( "log_Bs_RFD:=log(Bs_RFD)","B_{s} RFD", "", 'F' );
   //factory->AddVariable( "Bs_PT","Bs p_t","MeV", 'D' );
   //factory->AddVariable( "Bs_ENDVERTEX_CHI2", "Bs Vertex fit", "", 'D' );
   //factory->AddVariable( "PV_CHI2NDOF", "PV Fit chi2", "", 'D' );
   //factory->AddVariable( "log_XsDaughters_min_PT := log(XsDaughters_min_PT)","X_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   //factory->AddVariable( "log_DsDaughters_min_PT := log(DsDaughters_min_PT)","D_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   //factory->AddVariable( "log_K_1_1270_plus_IPCHI2_OWNPV := log(K_1_1270_plus_IPCHI2_OWNPV)","X_{s} ln(IP #chi^{2})", "", 'D' );   
   //factory->AddVariable( "log_XsDaughters_max_IPCHI2 := log(XsDaughters_max_IPCHI2)","X_{s} daughters max ln(IP#chi^{2})", "", 'F' );
   //factory->AddVariable( "log_DsDaughters_max_IPCHI2 := log(DsDaughters_max_IPCHI2)","D_{s} daughters max ln(IP#chi^{2})", "", 'F' );
   //factory->AddVariable( "Ds_max_DOCA","D_{s} max DOCA", "mm", 'F' );
   //factory->AddVariable( "log_Bs_FDCHI2_OWNPV := log(Bs_FDCHI2_OWNPV)","B_{s} ln(FD #chi^{2})", "", 'D' );
   //factory->AddVariable( "log_Ds_DIRA := log(1-Ds_DIRA_OWNPV)","ln(1 - D_{s} DIRA)","", 'D' );
   //factory->AddVariable( "K_plus_fromDs_ptasy_1.00","K^{+} (from D_{s}) cone p_{t} asy","", 'D' );
   //factory->AddVariable( "K_minus_fromDs_ptasy_1.00","K^{-} (from D_{s}) cone p_{t} asy","", 'D' );
   //factory->AddVariable( "pi_minus_fromDs_ptasy_1.00","#pi^{+} (from D_{s}) cone p_{t} asy","", 'D' );
   //factory->AddVariable( "K_plus_ptasy_1.00","K^{+} cone p_{t} asy","", 'D' );
   //factory->AddVariable( "pi_plus_ptasy_1.00","#pi^{+} cone p_{t}","", 'D' );
   //factory->AddVariable( "pi_minus_ptasy_1.00","#pi^{-} cone p_{t} asy","", 'D' );
   
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   factory->SetSignalWeightExpression("weight");

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "Bs_MM > 5300 && Bs_MM < 5420"  ; 
   TCut mycutb = "Bs_MM > 5650";
   
   // Tell the factory how to use the training and testing events
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
