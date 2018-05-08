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
/*
#include "TMVAGui.C"
#include "/work/dargent/TMVA-v4.2.0/test/variables.C"
#include "/work/dargent/TMVA-v4.2.0/test/efficiencies.C"
#include "/work/dargent/TMVA-v4.2.0/test/mvas.C"
#include "/work/dargent/TMVA-v4.2.0/test/correlations.C"*/

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "BDTG", TString run = "run1", TString trigger = "t0")
{
   TChain* background = new TChain("MCDecayTree");
   background->Add("GenMC.root");
  
   TChain* signal = new TChain("DecayTree");
   signal->Add("SelMC.root");

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

   //TString outDir = "plots";
   TString outDir = "../TD-AnaNote/latex/figs/TMVA/";
   outDir +=  myMethodList + "_" + run + "_" + trigger;
 
   TString outfileName = "TMVA_Bs2DsKpipi_" +  myMethodList + "_" + run + "_" + trigger + ".root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_" + run + "_" + trigger, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


   signal->SetBranchStatus("*",0);  // disable all branches   signal->SetBranchStatus("s_*",1); 
   signal->SetBranchStatus("s_*",1); 
   signal->SetBranchStatus("cos_*",1); 
   signal->SetBranchStatus("phi_*",1); 
   //signal->SetBranchStatus("weight",1); 

   background->SetBranchStatus("*",0);  // disable all branches
   background->SetBranchStatus("s_*",1);
   background->SetBranchStatus("cos_*",1); 
   background->SetBranchStatus("phi_*",1); 

   // Define the input variables that shall be used for the MVA training
   //factory->AddVariable( "s_Kpipi", "m(K#pi#pi)", "MeV", 'F' );
   factory->AddVariable( "s_Kpi", "m(K#pi)", "MeV", 'F' );
   //factory->AddVariable( "s_pipi", "m(#pi#pi)", "MeV", 'F' );
   factory->AddVariable( "s_Dspi", "m(D_{s}#pi)", "MeV", 'F' );
   //factory->AddVariable( "s_Dspipi", "m(D_{s}#pi#pi)", "MeV", 'F' );
   
   //factory->AddVariable( "s_DsK", "m(D_{s}K)", "MeV", 'F' );
   //factory->AddVariable( "s_DsKpi", "m(D_{s}K#pi)", "MeV", 'F' );
   //factory->AddVariable( "s_Kpip", "m(Kpi)", "MeV", 'F' );

   factory->AddVariable( "cos_theta_Kpi", "cos(theta_{K#pi})", "", 'F' );
   factory->AddVariable( "cos_theta_Dspi", "cos(theta_{D_{s}#pi})", "", 'F' );
   factory->AddVariable( "phi_Kpi_Dspi", "phi_{K#pi,D_{s}#pi}", "", 'F' );

   //factory->AddVariable( "cos_theta_pipi", "cos_theta_pipi", "", 'F' );
   //factory->AddVariable( "theta_DsK:=acos(cos_theta_DsK)", "theta_DsK", "", 'F' );
   //factory->AddVariable( "phi_pipi_DsK", "phi_pipi_DsK", "", 'F' );

   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "s_Kpipi > 1000 && s_Kpipi < 1950 && s_Kpi < 1200 && s_pipi < 1200 && s_Kpi > 200 && s_Dspi > 1900 && s_Dspi < 5000 && s_Dspipi > 2000 && s_Dspipi < 5100 && cos_theta_Dspi > 0"; 
   //mycuts = "run == " + run.ReplaceAll("run","")  ;
   //mycuts += "TriggerCat == " + trigger.ReplaceAll("t","")  ;
   TCut mycutb = "s_Kpipi > 1000 && s_Kpipi < 1950 && s_Kpi < 1200 && s_pipi < 1200 && s_Kpi > 200 && s_Dspi > 1900 && s_Dspi < 5000 && s_Dspipi > 2000 && s_Dspipi < 5100 && cos_theta_Dspi > 0 ";
   
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=12000:nTrain_Background=150000:nTest_Background=15000:SplitMode=Random:NormMode=NumEvents:!V" );

   //factory->SetSignalWeightExpression("weight");
   // ---- Book MVA methods

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
	//       factory->BookMethod( TMVA::Types::kBDT, "BDTG",
	//                            "!H:!V:NTrees=400:MinNodeSize=0%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=2:NegWeightTreatment=Pray:nEventsMin=500" );
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=500:MinNodeSize=0.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=2:NegWeightTreatment=Pray" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=500:MinNodeSize=3.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=80" );

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

   //variables(outfileName,"InputVariables_Id", "TMVA Input Variables",kFALSE, kTRUE, outDir);
   //correlations( outfileName,  kFALSE, kFALSE, kTRUE ,outDir);
   //efficiencies( outfileName,  2, kTRUE ,outDir);
   //mvas( outfileName, CompareType,  kTRUE , outDir, true);

   // Launch the GUI for the root macros
//    if (!gROOT->IsBatch()) TMVAGui( outfileName );
 

}


void trainAll( TString myMethodList = "BDTG", TString trainOn = "MC") {

	gROOT->SetBatch(true);
	TMVAClassification( myMethodList, "run1",  "t0" );


}
