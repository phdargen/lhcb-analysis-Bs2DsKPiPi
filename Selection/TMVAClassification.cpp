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
#include "/work/dargent/TMVA-v4.2.0/test/variables.C"
#include "/work/dargent/TMVA-v4.2.0/test/efficiencies.C"
#include "/work/dargent/TMVA-v4.2.0/test/mvas.C"
#include "/work/dargent/TMVA-v4.2.0/test/correlations.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "BDTG", TString trainOn = "MC", TString run = "run1", TString trigger = "t0", TString sample = "even" )
{
   TChain* background = new TChain("DecayTree");
   if(run == "run1" || run == "all"){
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_11.root");
  	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_12.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_11.root");
  	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_12.root");
// 	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_12_SS.root");
   }
   if(run == "run2" || run == "all") {
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_15.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_16.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_17.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2*_18.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_15.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_16.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_17.root");
	   background->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_Ds2*_18.root");
   }

   TChain* signal = new TChain("DecayTree");
   //signal->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_11.root");
   //signal->Add("/auto/data/dargent/BsDsKpipi/Preselected/MC/signal_Ds2KKpi_12.root");
   if(trainOn == "MC")signal->Add("/auto/data/dargent/BsDsKpipi/BDT/MC/signal.root");
   else signal->Add("/auto/data/dargent/BsDsKpipi/Preselected/Data/norm_sweight.root");

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
   outDir +=  myMethodList + "_" + trainOn + "_" + run + "_" + trigger + "_" + sample;
 
   TString outfileName = "TMVA_Bs2DsKpipi_" +  myMethodList + "_" + trainOn + "_" + run + "_" + trigger + "_" + sample + ".root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_" + trainOn + "_" + run + "_" + trigger + "_" + sample, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


   signal->SetBranchStatus("*",0);  // disable all branches
   signal->SetBranchStatus("*CHI2*",1); 
   signal->SetBranchStatus("*Chi2*",1); 
   signal->SetBranchStatus("*DOCA*",1);
   signal->SetBranchStatus("*DIRA*",1);
   signal->SetBranchStatus("*PT*",1);
   signal->SetBranchStatus("*RFD*",1);
   signal->SetBranchStatus("*max*",1);
   signal->SetBranchStatus("*ptasy*",1);
   signal->SetBranchStatus("*MM*",1);
   signal->SetBranchStatus("*TAU*",1);
   signal->SetBranchStatus("*Trigger*",1);
   signal->SetBranchStatus("*PT",1);
   signal->SetBranchStatus("*m1*",1);
   signal->SetBranchStatus("*PIDK",1);
   signal->SetBranchStatus("*BKG*",1);
   signal->SetBranchStatus("*State*",1);
   signal->SetBranchStatus("m*",1);
   signal->SetBranchStatus("weight",1);
   signal->SetBranchStatus("eventNumber",1);
   signal->SetBranchStatus("run",1);
   signal->SetBranchStatus("Ds_finalState",1);
   signal->SetBranchStatus("Ds_finalState",1);
   signal->SetBranchStatus("Bs_PV_nPV",1);
   signal->SetBranchStatus("Bs_MINIPCHI2NEXTBEST",1);
   signal->SetBranchStatus("Bs_IPCHI2_OWNPV",1);
   signal->SetBranchStatus("Bs_BsDTF_TAU",1);
   //signal->SetBranchStatus("K_1_1270_plus_MM",1);

   background->SetBranchStatus("*",0);  // disable all branches
   background->SetBranchStatus("*CHI2*",1); 
   background->SetBranchStatus("*Chi2*",1); 
   background->SetBranchStatus("*DOCA*",1);
   background->SetBranchStatus("*DIRA*",1);
   background->SetBranchStatus("*PT*",1);
   background->SetBranchStatus("*RFD*",1);
   background->SetBranchStatus("*max*",1);
   background->SetBranchStatus("*ptasy*",1);
   background->SetBranchStatus("*MM*",1);
   background->SetBranchStatus("*Trigger*",1);
   background->SetBranchStatus("*PT",1);
   background->SetBranchStatus("*m1*",1);
   background->SetBranchStatus("*PIDK",1);
   background->SetBranchStatus("*State*",1);
   background->SetBranchStatus("*TAU*",1);
   background->SetBranchStatus("m*",1);
   background->SetBranchStatus("weight",1);
   background->SetBranchStatus("eventNumber",1);
   background->SetBranchStatus("run",1);
   background->SetBranchStatus("Ds_finalState",1);
   background->SetBranchStatus("Bs_PV_nPV",1);
   background->SetBranchStatus("Bs_MINIPCHI2NEXTBEST",1);
   background->SetBranchStatus("Bs_IPCHI2_OWNPV",1);
   background->SetBranchStatus("Bs_BsDTF_TAU",1);
   //background->SetBranchStatus("K_1_1270_plus_MM",1);

   // Define the input variables that shall be used for the MVA training

   //factory->AddVariable( "DTF_CHI2NDOF", "chi^{2}_{DTF}/#nu", "", 'F' );
   factory->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)","B_{s} ln(IP #chi^{2})", "", 'F' );
   factory->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)","B_{s} ln(1 - DIRA)","", 'F' );
   factory->AddVariable( "PV_CHI2NDOF", "#chi^{2}_{DTF}/ndf", "", 'F' );
   factory->AddVariable( "log_Bs_SmallestDeltaChi2OneTrack:= log(Bs_SmallestDeltaChi2OneTrack)","#Delta#chi^{2}_{add-track}","", 'F' );
   factory->AddVariable( "Bs_ptasy_1.00","B_{s} A_{p_{t}}^{cone}","", 'F' );
   factory->AddVariable("max_ghostProb","max[ghostProb]","",'F');
  
   factory->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)","X_{s} daughters min[ln(IP#chi^{2})]", "", 'F' );
   //factory->AddVariable( "Xs_ptasy_1.00","X_{s} cone p_{t} asy","", 'F' );
   factory->AddVariable( "Xs_max_DOCA","X_{s} max[DOCA]", "mm", 'F' );

   factory->AddVariable( "maxCos", "cos(max[#theta_{Ds h}])", "", 'F' );

   factory->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)","D_{s} daughters min[ln(IP#chi^{2})]", "", 'F' );
   //factory->AddVariable( "Ds_ptasy_1.00","D_{s} cone p_{t} asy","", 'F' );
   factory->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)","D_{s} ln(#chi^{2}_{FD})", "", 'F' );
   factory->AddVariable( "log_Ds_RFD:=log(Ds_RFD)","D_{s} log(RFD)", "", 'F' );


   // Additional variables for testing
  
   //factory->AddVariable("max_ProbNNghost","max_ProbNNghost","",'F');
   //factory->AddVariable( "log_track_min_PT := log(track_min_PT)","track_min_PT","", 'F' );
   //factory->AddVariable("Bs_DTF_TAU","Bs_DTF_TAU","",'F');
   //factory->AddVariable("m_Kpipi","m_Kpipi","",'F');
   //factory->AddVariable("m_Dspipi","m_Dspipi","",'F');
   //factory->AddVariable( "Bs_MINIPCHI2NEXTBEST := log(Bs_MINIPCHI2NEXTBEST)","Bs_MINIPCHI2NEXTBEST","", 'F' );
   //factory->AddVariable( "Bs_SmallestDeltaChi2MassOneTrack:= log(Bs_SmallestDeltaChi2MassOneTrack)","Bs_SmallestDeltaChi2MassOneTrack","", 'F' );
   //factory->AddVariable( "Ds_SmallestDeltaChi2OneTrack:= log(Ds_SmallestDeltaChi2OneTrack)","Bs_SmallestDeltaChi2OneTrack","", 'F' );
   //factory->AddVariable( "Ds_SmallestDeltaChi2MassOneTrack:= log(Ds_SmallestDeltaChi2MassOneTrack)","Bs_SmallestDeltaChi2MassOneTrack","", 'F' );
   //factory->AddVariable( "Bs_SmallestDeltaChi2TwoTracks:= log(Bs_SmallestDeltaChi2TwoTracks)","Bs_SmallestDeltaChi2TwoTrack","", 'F' );
   //factory->AddVariable( "Bs_SmallestDeltaChi2MassTwoTracks:= log(Bs_SmallestDeltaChi2MassTwoTracks)","Bs_SmallestDeltaChi2MassTwoTrack","", 'F' );
   //factory->AddVariable("Bs_DTF_MM","Bs_DTF_MM","",'F');
   //factory->AddVariable( "log_track_min_IPCHI2 := log(track_min_IPCHI2)","track_min_IPCHI2","", 'F' );
   //factory->AddVariable( "log_Bs_PT := log(Bs_PT)","Bs p_t","MeV", 'D' );
   //factory->AddVariable( "log_Ds_PT := log(Ds_PT)","Ds p_t","MeV", 'D' );
   //factory->AddVariable( "Bs_ptasy_1.00","Bs  cone p_{t} asy","", 'F' );
   //factory->AddVariable( "Ds_max_DOCA","D_{s} max DOCA", "mm", 'F' );
   //factory->AddVariable( "Ds_finalState","f","", 'I' );
   //factory->AddVariable( "Ds_m12","Ds_m12","", 'F' );
   //factory->AddVariable( "Ds_m13","Ds_m13","", 'F' );
   //factory->AddVariable( "K_plus_fromDs_PIDK","K_plus_fromDs_PIDK","", 'F' );
   //factory->AddVariable( "K_minus_fromDs_PIDK","K_minus_fromDs_PIDK","", 'F' );
   //factory->AddVariable( "pi_minus_fromDs_PIDK","pi_minus_fromDs_PIDK","", 'F' );
   //factory->AddVariable( "K_plus_PIDK","K_plus_PIDK","", 'F' );
   //factory->AddVariable( "log_Bs_RFD:=log(Bs_RFD)","B_{s} RFD", "", 'F' );
   //factory->AddVariable( "Ds_ENDVERTEX_CHI2", "Ds Vertex fit", "", 'D' );
   //factory->AddVariable( "PV_CHI2NDOF", "PV Fit chi2", "", 'D' );
   //factory->AddVariable( "log_XsDaughters_min_PT := log(XsDaughters_min_PT)","X_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   //factory->AddVariable( "log_DsDaughters_min_PT := log(DsDaughters_min_PT)","D_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   //factory->AddVariable( "log_Ds_IPCHI2_OWNPV := log(Ds_IPCHI2_OWNPV)","D_{s} ln(IP #chi^{2})", "", 'D' );   
   //factory->AddVariable( "log_XsDaughters_max_IPCHI2 := log(XsDaughters_max_IPCHI2)","X_{s} daughters max ln(IP#chi^{2})", "", 'F' );
   //factory->AddVariable( "log_DsDaughters_max_IPCHI2 := log(DsDaughters_max_IPCHI2)","D_{s} daughters max ln(IP#chi^{2})", "", 'F' );
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

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts; 
   mycuts = "maxCos > -0.95 && Bs_IPCHI2_OWNPV < 16 && Bs_SmallestDeltaChi2OneTrack > 5 && Ds_finalState < 5 && PV_CHI2NDOF < 15 && Xs_max_DOCA < 0.5 && m_pipi < 1200 && m_Kpi < 1200 && (Bs_PV_nPV ==1 || (Bs_MINIPCHI2NEXTBEST-Bs_IPCHI2_OWNPV) > 20) && Bs_BsDTF_TAU >= 0.4";
   if(run != "all")mycuts += "run == " + run.ReplaceAll("run","");
   if(trigger != "all")mycuts += "TriggerCat == " + trigger.ReplaceAll("t","");
   if(trainOn == "MC") mycuts += "Bs_MM > 5300 && Bs_MM < 5420 && Bs_BKGCAT == 20";

   TCut mycutb = "Bs_MM > 5500";
   mycutb += "maxCos > -0.95 && Bs_IPCHI2_OWNPV < 16 && Bs_SmallestDeltaChi2OneTrack > 5 && Ds_finalState < 5 && PV_CHI2NDOF < 15 && Xs_max_DOCA < 0.5 && m_pipi < 1200 && m_Kpi < 1200 && (Bs_PV_nPV ==1 || (Bs_MINIPCHI2NEXTBEST-Bs_IPCHI2_OWNPV) > 20) && Bs_BsDTF_TAU >= 0.4" ;
   if(run != "all")mycutb += "run == " + run;
   if(trigger != "all")mycutb += "TriggerCat == " + trigger ;
   
   TFile* dummy;

   if(sample == "all"){
	   factory->AddSignalTree    ( signal,     signalWeight     );
   	   factory->AddBackgroundTree( background, backgroundWeight );
           factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   }
   else {
	dummy = new TFile("dummy.root","RECREATE");
	TCut cut_even = " eventNumber % 2 == 0";
	TCut cut_odd =  " eventNumber % 2 == 1";

        TTree* signal_training;
        TTree* signal_testing;
        TTree* background_training;
        TTree* background_testing;

	if(sample == "even"){
		signal_training = signal->CopyTree(mycuts+cut_even);
		signal_testing = signal->CopyTree(mycuts+cut_odd);
		background_training = background->CopyTree(mycutb+cut_even);
		background_testing = background->CopyTree(mycutb+cut_odd);
	}
	else {
		signal_training = signal->CopyTree(mycuts+cut_odd);
		signal_testing = signal->CopyTree(mycuts+cut_even);
		background_training = background->CopyTree(mycutb+cut_odd);
		background_testing = background->CopyTree(mycutb+cut_even);
        }
	factory->AddSignalTree(signal_training,signalWeight,TMVA::Types::kTraining);
	factory->AddSignalTree(signal_testing,signalWeight,TMVA::Types::kTesting);
	factory->AddBackgroundTree(background_training,backgroundWeight,TMVA::Types::kTraining);
	factory->AddBackgroundTree(background_testing,backgroundWeight,TMVA::Types::kTesting);
   }

   factory->SetSignalWeightExpression("weight");

   // ---- Book MVA methods

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=500:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=3:NegWeightTreatment=Pray" );

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

   variables(outfileName,"InputVariables_Id", "TMVA Input Variables",kFALSE, kTRUE, outDir);
   correlations( outfileName,  kFALSE, kFALSE, kTRUE ,outDir);
   efficiencies( outfileName,  2, kTRUE ,outDir);
   mvas( outfileName, CompareType,  kTRUE , outDir, true);

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
 

}


void trainAll( TString myMethodList = "BDTG", TString trainOn = "Data") {

	gROOT->SetBatch(true);
	TMVAClassification( myMethodList, trainOn , "all",  "all", "all" );

	TMVAClassification( myMethodList, trainOn , "run1",  "t0", "even" );
 	TMVAClassification( myMethodList, trainOn , "run1",  "t1", "even" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t0", "even" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t1", "even" );

	TMVAClassification( myMethodList, trainOn , "run1",  "t0", "odd" );
 	TMVAClassification( myMethodList, trainOn , "run1",  "t1", "odd" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t0", "odd" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t1", "odd" );

	TMVAClassification( myMethodList, trainOn , "run1",  "t0", "all" );
 	TMVAClassification( myMethodList, trainOn , "run1",  "t1", "all" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t0", "all" );
 	TMVAClassification( myMethodList, trainOn , "run2",  "t1", "all" );
}
