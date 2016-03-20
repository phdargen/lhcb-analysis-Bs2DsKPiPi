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

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <TChain.h>
#include "TStopwatch.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void TMVAClassificationApplication( TString myMethodList = "BDTG" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   TChain* theTree = new TChain("DecayTree");
   theTree->Add("/auto/data/dargent/Bs2DsKpipi/preselection_norm/data2011_Ds2KKpi_forBDT.root");
   //theTree->Add("/auto/data/dargent/Bs2DsKpipi/preselection_norm/data2012_Ds2KKpi_forBDT.root");
   //theTree->Add("/auto/data/dargent/Bs2DsKpipi/preselection_norm/mc11_Ds2KKpi_forBDT.root");
   //theTree->Add("/auto/data/dargent/Bs2DsKpipi/preselection_norm/mc12_Ds2KKpi_forBDT.root");

   TFile *hFile = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection_norm/data11_Ds2KKpi_BDT.root","RECREATE");
   //TFile *hFile = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection_norm/data12_Ds2KKpi_BDT.root","RECREATE");
   //TFile *hFile = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection_norm/mc11_Ds2KKpi_BDT.root","RECREATE");
   //TFile *hFile = new TFile("/auto/data/dargent/Bs2DsKpipi/preselection_norm/mc12_Ds2KKpi_BDT.root","RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Cut optimisation
   Use["Cuts"]            = 1;
   Use["CutsD"]           = 1;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 1;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 1;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 1;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 1; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 1; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 1;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

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
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  /* Float_t var1, var2;
   Float_t var3, var4;
   reader->AddVariable( "myvar1 := var1+var2", &var1 );
   reader->AddVariable( "myvar2 := var1-var2", &var2 );
   reader->AddVariable( "var3",                &var3 );
   reader->AddVariable( "var4",                &var4 );*/

   Float_t var[100];

   reader->AddVariable( "Bs_ENDVERTEX_CHI2", &var[0]);
   reader->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)", &var[1]);
   reader->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)", &var[2]);

   reader->AddVariable( "Xs_max_DOCA", &var[3]);
   reader->AddVariable( "log_Bs_FDCHI2_OWNPV := log(Bs_FDCHI2_OWNPV)", &var[4]);
   reader->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)", &var[5]);
   reader->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)", &var[6]);

   reader->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)", &var[7] );
   reader->AddVariable( "log_Ds_DIRA := log(1-Ds_DIRA_OWNPV)", &var[8] );

   reader->AddVariable( "K_plus_fromDs_ptasy_1.00", &var[9] );
   reader->AddVariable( "K_minus_fromDs_ptasy_1.00", &var[10] );
   reader->AddVariable( "pi_minus_fromDs_ptasy_1.00", &var[11] );
   reader->AddVariable( "K_plus_ptasy_1.00", &var[12] );
   reader->AddVariable( "pi_plus_ptasy_1.00", &var[13] );
   reader->AddVariable( "pi_minus_ptasy_1.00", &var[14] );

   reader->AddVariable( "max_ghostProb := max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb)))))", &var[15] );
   reader->AddVariable( "cos := cos( max(max(angK,angPip),angPim))   ", &var[16] );

   //reader->AddVariable( "log(Bs_PT)", &var1);
   //reader->AddVariable( "log_XsDaughters_min_PT := log(XsDaughters_min_PT)", &var3);
   //reader->AddVariable( "log_DsDaughters_min_PT := log(DsDaughters_min_PT)", &var4);
   //reader->AddVariable( "log_XsDaughters_max_IPCHI2 := log(XsDaughters_max_IPCHI2)", &var7);
   //reader->AddVariable( "log_DsDaughters_max_IPCHI2 := log(DsDaughters_max_IPCHI2)", &var8);
   //reader->AddVariable( "log_K_1_1270_plus_IPCHI2_OWNPV := log(K_1_1270_plus_IPCHI2_OWNPV)", &var13 );
   //reader->AddVariable( "min_ghostProb := min(pi_minus_TRACK_GhostProb,min(pi_plus_TRACK_GhostProb,min(K_plus_TRACK_GhostProb,min(K_plus_fromDs_TRACK_GhostProb,min(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb)))))", &var33 );

/*
   reader->AddVariable( "Bs_ETA", &var16 );
   reader->AddVariable( "Ds_ETA", &var17 );
   reader->AddVariable( "K_plus_fromDs_ETA", &var18 );
   reader->AddVariable( "K_minus_fromDs_ETA", &var19 );
   reader->AddVariable( "pi_minus_fromDs_ETA", &var20 );
   reader->AddVariable( "K_plus_ETA", &var21 );
   reader->AddVariable( "pi_plus_ETA", &var22 );
   reader->AddVariable( "pi_minus_ETA", &var23 );
*/
   //reader->AddVariable( "max_TrackChi2", &var24 );
   //reader->AddVariable( "min_TrackChi2", &var25 );

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
   
   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   Double_t myVar1, myVar2;
   Float_t  myVar3, myVar4, myVar5, myVar6, myVar7, myVar8, myVar9;
   Double_t myVar10, myVar11, myVar12, myVar13 , myVar14 ,myVar15, myVar16, myVar17, myVar18, myVar19, myVar20;
   Double_t myVar21, myVar22, myVar23 ,myVar25, myVar26, myVar27, myVar28, myVar29, myVar30;
   Float_t myVar24;
   Double_t myVar31, myVar32, myVar33, myVar34, myVar35, myVar36, myVar37, myVar38, myVar39, myVar40;

   theTree->SetBranchAddress( "Bs_PT", &myVar1 );
   theTree->SetBranchAddress( "Bs_ENDVERTEX_CHI2",&myVar2 );
   theTree->SetBranchAddress( "XdDaughters_min_PT", &myVar3 );
   theTree->SetBranchAddress( "DsDaughters_min_PT", &myVar4 );
   theTree->SetBranchAddress( "XdDaughters_min_IPCHI2", &myVar5 );
   theTree->SetBranchAddress( "DsDaughters_min_IPCHI2", &myVar6 );
   theTree->SetBranchAddress( "XdDaughters_max_IPCHI2", &myVar7 );
   theTree->SetBranchAddress( "DsDaughters_max_IPCHI2", &myVar8 );
   theTree->SetBranchAddress( "Xd_max_DOCA", &myVar9 );
   theTree->SetBranchAddress( "Bs_FDCHI2_OWNPV", &myVar10 );
   theTree->SetBranchAddress( "Ds_FDCHI2_ORIVX", &myVar11 );
   theTree->SetBranchAddress( "Bs_IPCHI2_OWNPV", &myVar12 );
   theTree->SetBranchAddress( "a_1_1260_plus_IPCHI2_OWNPV", &myVar13 );
   theTree->SetBranchAddress( "Bs_DIRA_OWNPV", &myVar14 );
   theTree->SetBranchAddress( "Ds_DIRA_OWNPV", &myVar15 );
   theTree->SetBranchAddress( "Bs_ETA", &myVar16 );
   theTree->SetBranchAddress( "Ds_ETA", &myVar17 );
   theTree->SetBranchAddress( "K_plus_fromDs_ETA", &myVar18 );
   theTree->SetBranchAddress( "K_minus_fromDs_ETA", &myVar19 );
   theTree->SetBranchAddress( "pi_minus_fromDs_ETA", &myVar20 );
   theTree->SetBranchAddress( "pi_plus1_ETA", &myVar21 );
   theTree->SetBranchAddress( "pi_plus2_ETA", &myVar22 );
   theTree->SetBranchAddress( "pi_minus_ETA", &myVar23 );
   theTree->SetBranchAddress( "max_TrackChi2", &myVar24 );
   theTree->SetBranchAddress( "K_plus_fromDs_ptasy_1.00", &myVar26 );
   theTree->SetBranchAddress( "K_minus_fromDs_ptasy_1.00", &myVar27 );
   theTree->SetBranchAddress( "pi_minus_fromDs_ptasy_1.00", &myVar28 );
   theTree->SetBranchAddress( "pi_plus1_ptasy_1.00", &myVar29 );
   theTree->SetBranchAddress( "pi_plus2_ptasy_1.00", &myVar30 );
   theTree->SetBranchAddress( "pi_minus_ptasy_1.00", &myVar31 );
   theTree->SetBranchAddress( "K_plus_fromDs_TRACK_GhostProb", &myVar32 );
   theTree->SetBranchAddress( "K_minus_fromDs_TRACK_GhostProb", &myVar33 );
   theTree->SetBranchAddress( "pi_minus_fromDs_TRACK_GhostProb", &myVar34 );
   theTree->SetBranchAddress( "pi_plus1_TRACK_GhostProb", &myVar35 );
   theTree->SetBranchAddress( "pi_plus2_TRACK_GhostProb", &myVar36 );
   theTree->SetBranchAddress( "pi_minus_TRACK_GhostProb", &myVar37 );
   theTree->SetBranchAddress( "angK", &myVar38 );
   theTree->SetBranchAddress( "angPip", &myVar39 );
   theTree->SetBranchAddress( "angPim", &myVar40 );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   Float_t BDTG_response;
   tree->Branch("BDTG_response",&BDTG_response, "BDTG_response/F");

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;

   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);

	var[0] = myVar2;
        var[1] = log(myVar5);
        var[2] = log(myVar6);
        var[3] = myVar9;
        var[4] = log(myVar10);
        var[5] = log(myVar11);
        var[6] = log(myVar12);
        var[7] = log(1.-myVar14);
        var[8] = log(1.-myVar15);
        var[9] = myVar26;
        var[10] = myVar27;
        var[11] = myVar28;
        var[12] = myVar29;
        var[13] = myVar30;
        var[14] = myVar31;
        var[15] = max(myVar32,max(myVar33,max(myVar34,max(myVar35,max(myVar36,myVar37)))));
	var[16] = cos(max(max(myVar38,myVar39),myVar40));

        BDTG_response=reader->EvaluateMVA("BDTG method");
        tree->Fill();


      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
     // if (Use["BDT"          ])   histBdt    ->Fill( BDT_response );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( BDTG_response );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // --- Write histograms

   //TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();

   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   
   tree->Write();
   delete reader;
 //  input->Close();
   hFile->Close(); 
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
} 
