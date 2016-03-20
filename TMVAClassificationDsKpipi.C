// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

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

void TMVAClassification( TString myMethodList = "BDTG" )
{
   TChain* background = new TChain("DecayTree");
   background->Add("/auto/data/dargent/Bs2DsKpipi/preselection/data2011_Ds2KKpi_forBDT.root");
   background->Add("/auto/data/dargent/Bs2DsKpipi/preselection/data2012_Ds2KKpi_forBDT.root");

   TChain* signal = new TChain("DecayTree");
   signal->Add("/auto/data/dargent/Bs2DsKpipi/preselection/mc11_Ds2KKpi_forBDT.root");
   signal->Add("/auto/data/dargent/Bs2DsKpipi/preselection/mc12_Ds2KKpi_forBDT.root");

   //---------------------------------------------------------------
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
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 1;
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
   TString outfileName( "TMVA_Bs2DsKpipi_2012Ana.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   // Define the input variables that shall be used for the MVA training
   //factory->AddVariable( "Bs_PT","Bs p_t","MeV", 'D' );
   factory->AddVariable( "Bs_ENDVERTEX_CHI2", "Bs Vertex fit", "", 'D' );
   //factory->AddVariable( "log_XsDaughters_min_PT := log(XsDaughters_min_PT)","X_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   //factory->AddVariable( "log_DsDaughters_min_PT := log(DsDaughters_min_PT)","D_{s} daughters ln(p_{t})", "ln(MeV)", 'F' );
   factory->AddVariable( "log_XsDaughters_min_IPCHI2 := log(XsDaughters_min_IPCHI2)","X_{s} daughters min ln(IP#chi^{2})", "", 'F' );
   factory->AddVariable( "log_DsDaughters_min_IPCHI2 := log(DsDaughters_min_IPCHI2)","D_{s} daughters min ln(IP#chi^{2})", "", 'F' );
   //factory->AddVariable( "log_XsDaughters_max_IPCHI2 := log(XsDaughters_max_IPCHI2)","X_{s} daughters max ln(IP#chi^{2})", "", 'F' );
   //factory->AddVariable( "log_DsDaughters_max_IPCHI2 := log(DsDaughters_max_IPCHI2)","D_{s} daughters max ln(IP#chi^{2})", "", 'F' );
   factory->AddVariable( "Xs_max_DOCA","X_{s} max DOCA", "mm", 'F' );
   factory->AddVariable( "log_Bs_FDCHI2_OWNPV := log(Bs_FDCHI2_OWNPV)","B_{s} ln(FD #chi^{2})", "", 'D' );
   factory->AddVariable( "log_Ds_FDCHI2_ORIVX := log(Ds_FDCHI2_ORIVX)","D_{s} FD significance", "", 'D' );

   factory->AddVariable( "log_Bs_IPCHI2_OWNPV := log(Bs_IPCHI2_OWNPV)","B_{s} ln(IP #chi^{2})", "", 'D' );
   //factory->AddVariable( "log_K_1_1270_plus_IPCHI2_OWNPV := log(K_1_1270_plus_IPCHI2_OWNPV)","X_{s} ln(IP #chi^{2})", "", 'D' );



   //investigate possible variables

   //DIRA variables
   factory->AddVariable( "log_Bs_DIRA := log(1-Bs_DIRA_OWNPV)","ln(1 - B_{s} DIRA)","", 'D' );
   factory->AddVariable( "log_Ds_DIRA := log(1-Ds_DIRA_OWNPV)","ln(1 - D_{s} DIRA)","", 'D' );

/*
*/
/*
   //pseudo rapidity variables
   factory->AddVariable( "Bs_ETA","B_{s} #eta","", 'D' );
   factory->AddVariable( "Ds_ETA","D{s} #eta","", 'D' );
   factory->AddVariable( "K_plus_fromDs_ETA","K^{+} (from D_{s}) #eta","", 'D' );
   factory->AddVariable( "K_minus_fromDs_ETA","K^{-} (from D_{s}) #eta","", 'D' );
   factory->AddVariable( "pi_minus_fromDs_ETA","#pi^{+} (from D_{s}) #eta","", 'D' );
   factory->AddVariable( "K_plus_ETA","K^{+} #eta","", 'D' );
   factory->AddVariable( "pi_plus_ETA","#pi^{+} #eta","", 'D' );
   factory->AddVariable( "pi_minus_ETA","#pi^{-} #eta","", 'D' );
*/

   //track chi2
   //factory->AddVariable( "max_TrackChi2","max(track #chi^{2})","", 'F' );
  // factory->AddVariable( "min_TrackChi2","min(track #chi^{2})","", 'F' );

   //cone pt asy
   factory->AddVariable( "K_plus_fromDs_ptasy_1.00","K^{+} (from D_{s}) cone p_{t} asy","", 'D' );
   factory->AddVariable( "K_minus_fromDs_ptasy_1.00","K^{-} (from D_{s}) cone p_{t} asy","", 'D' );
   factory->AddVariable( "pi_minus_fromDs_ptasy_1.00","#pi^{+} (from D_{s}) cone p_{t} asy","", 'D' );
   factory->AddVariable( "K_plus_ptasy_1.00","K^{+} cone p_{t} asy","", 'D' );
   factory->AddVariable( "pi_plus_ptasy_1.00","#pi^{+} cone p_{t}","", 'D' );
   factory->AddVariable( "pi_minus_ptasy_1.00","#pi^{-} cone p_{t} asy","", 'D' );

   factory->AddVariable("max_ghostProb := max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb)))))","","",'D');
   
   factory->AddVariable( "cos := cos( max(max(angK,angPip),angPim))   ", "cos(max[#theta_{Ds h}])", "", 'F' );

   //factory->AddVariable( "K_plus_PIDK","K^{+} PIDK","", 'D' );

   //factory->AddVariable( "Ds_FDsig := (Ds_ENDVERTEX_Z - Ds_ORIVX_Z) / sqrt( pow(Ds_ENDVERTEX_ZERR,2) + pow(Ds_ORIVX_ZERR,2) )","D_{s} FD significance", "", 'D' );
  //factory->AddVariable( "Bs_RFD := sqrt((Bs_ENDVERTEX_X - Bs_OWNPV_X)*(Bs_ENDVERTEX_X - Bs_OWNPV_X) + (Bs_ENDVERTEX_Y - Bs_OWNPV_Y) * (Bs_ENDVERTEX_Y - Bs_OWNPV_Y))","B_{s} RFD", "", 'D' );

  //factory->AddVariable( "Ds_RFD := sqrt((Ds_ENDVERTEX_X - Ds_OWNPV_X)*(Ds_ENDVERTEX_X - Ds_OWNPV_X) + (Ds_ENDVERTEX_Y - Ds_OWNPV_Y) * (Ds_ENDVERTEX_Y - Ds_OWNPV_Y))","D_{s} RFD", "", 'D' );


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;
   
   // You can add an arbitrary number of signal or background trees
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
  
   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = "Bs_MM > 5300 && Bs_MM < 5420 && max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) > 0 && max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) < 1"; 
  //&& min(pi_minus_TRACK_GhostProb,min(pi_plus_TRACK_GhostProb,min(K_plus_TRACK_GhostProb,min(K_plus_fromDs_TRACK_GhostProb,min(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) > 0 &&  min(pi_minus_TRACK_GhostProb,min(pi_plus_TRACK_GhostProb,min(K_plus_TRACK_GhostProb,min(K_plus_fromDs_TRACK_GhostProb,min(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) < 1";
   TCut mycutb = "Bs_MM > 5600 && max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) > 0 && max(pi_minus_TRACK_GhostProb,max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,max(K_plus_fromDs_TRACK_GhostProb,max(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) < 1 && abs(Ds_MM - 1969) < 40   ";
   //&& min(pi_minus_TRACK_GhostProb,min(pi_plus_TRACK_GhostProb,min(K_plus_TRACK_GhostProb,min(K_plus_fromDs_TRACK_GhostProb,min(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) > 0 &&  min(pi_minus_TRACK_GhostProb,min(pi_plus_TRACK_GhostProb,min(K_plus_TRACK_GhostProb,min(K_plus_fromDs_TRACK_GhostProb,min(pi_minus_fromDs_TRACK_GhostProb,K_minus_fromDs_TRACK_GhostProb))))) < 1";

   // Tell the factory how to use the training and testing events
   factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods

   // Cut optimisation
   if (Use["Cuts"])
      factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

   if (Use["CutsD"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

   if (Use["CutsPCA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

   if (Use["CutsGA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

   if (Use["CutsSA"])
      factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
      factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

   // Decorrelated likelihood
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

   if (Use["PDERSD"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

   if (Use["PDERSPCA"])
      factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

   if (Use["PDEFoamBoost"])
      factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"])
      factory->BookMethod( TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"])
      factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"])
      factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher discriminant (same as LD)
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"])
      factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"])
      factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

   if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

   if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

   if (Use["FDA_MT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

   if (Use["FDA_GAMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

   if (Use["FDA_MCMT"])
      factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   if (Use["MLPBFGS"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

   if (Use["MLPBNN"])
      factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"])
      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

   // Tmlp(Root)ANN
   if (Use["TMlpANN"])
      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

   // Support Vector Machine
   if (Use["SVM"])
      factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

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

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"])
      factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's ! 
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

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
