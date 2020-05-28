// Time fits
// author: Philippe d'Argent
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
#include "Mint/Neg2LLMultiConstraint.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/DalitzPdfBase.h"
#include "Mint/DalitzPdfBaseFastInteg.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/AmpRatios.h"
#include "Mint/IEventGenerator.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/SignalGenerator.h"
#include "Mint/FromFileGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/cexp.h"
#include "Mint/DalitzPdfNormChecker.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/DalitzPdfSaveInteg.h"
#include "Mint/Chi2Binning.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpList.h"
#include "Mint/DalitzPdfBaseMCInteg.h"
#include "Mint/TimePdfMaster.h"
#include "Mint/FullTimePdf.h"

#include "RooRealConstant.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDecay.h"
#include "RooBDecay.h"
#include "RooPlot.h"
#include "RooEffProd.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooUniform.h"
#include "RooExponential.h"
#include "RooRandom.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooMultiVarGaussian.h"
#include "RooTruthModel.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
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
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDecompChol.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace MINT;

void fullTimeFit(int step=0, string mode = "fit"){

    /// Options
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    RooRandom::randomGenerator()->SetSeed(seed);

    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    DalitzEventPattern pat_CP = pat.makeCPConjugate();

    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/signal_tagged.root", (char*) 0);
    NamedParameter<string> InputGenMCFile("InputGenMCFile", (std::string) "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/EvtGen/GenLevMC/Gen_DsK.root", (char*) 0);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);

    NamedParameter<int>  useTrueTau("useTrueTau", 0);
    NamedParameter<int>  useTrueTagging("useTrueTagging", 0);

    NamedParameter<double> w_min("w_min", 0.);
    NamedParameter<double> w_max("w_max", 1.);

    NamedParameter<double> min_year("min_year", 10);
    NamedParameter<double> max_year("max_year", 20);
    NamedParameter<double> min_trigger("min_trigger", -10);
    NamedParameter<double> max_trigger("max_trigger", 10);
    NamedParameter<double> min_Ds_finalState("min_Ds_finalState", -10);
    NamedParameter<double> max_Ds_finalState("max_Ds_finalState", 10);

    NamedParameter<int>  nBinst("nBinst", 40);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 10);
    NamedParameter<int>  doPlots("doPlots", 1);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  doScanParam("doScanParam", 0);
    NamedParameter<string> scanParam("scanParam", (std::string) "", (char*) 0);
    NamedParameter<double>  scanParam_min("scanParam_min", 0.);
    NamedParameter<double>  scanParam_max("scanParam_max", 1.);
    NamedParameter<double>  scanParam_steps("scanParam_steps", 100.);

    NamedParameter<int>  doSimFit("doSimFit", 0);
    NamedParameter<int>  doSimFitInBins("doSimFitInBins", 0);

    NamedParameter<int>  fitGenMC("fitGenMC", 0);
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);

    NamedParameter<int>  doToyStudy("doToyStudy", 0);
    NamedParameter<double> N_scale_toys("N_scale_toys", 1);
    NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

    NamedParameter<int>  useGaussConstrainsTagging("useGaussConstrainsTagging", 0);
    NamedParameter<int>  useGaussConstrainsBias("useGaussConstrainsBias", 0);

    NamedParameter<int>  doAccSystematics("doAccSystematics", 0);
    NamedParameter<int>  useCholDec("useCholDec", 0);
    NamedParameter<int>  varPerParChol("varPerParChol", 100);
    int chol_index = (step-1)/varPerParChol ;

    NamedParameter<string> doSystematic("doSystematic", (std::string) "", (char*) 0);
    NamedParameter<int>  randomizeStartVals("randomizeStartVals", 0);

    NamedParameter<string> weightName("weightName", (std::string) "N_Bs_sw", (char*) 0);

    NamedParameter<int>  dilutionWeight("dilutionWeight", 1);
    NamedParameter<int>  N_plot_it("N_plot_it",10);
    NamedParameter<int> scale_asym("scale_asym", 1);


    /// Common fit parameters
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);
    FitParameter  k("k",1,1,1.);
        
    FitParameter  C("C",1,0.,0.1);
    FitParameter  D("D",1,0.,0.1);
    FitParameter  D_bar("D_bar",1,0.,0.1);
    FitParameter  S("S",1,0.,0.1);
    FitParameter  S_bar("S_bar",1,0.,0.1);
    
    FitParameter  C_1("C_1",1,0.,0.1);
    FitParameter  D_1("D_1",1,0.,0.1);
    FitParameter  D_bar_1("D_bar_1",1,0.,0.1);
    FitParameter  S_1("S_1",1,0.,0.1);
    FitParameter  S_bar_1("S_bar_1",1,0.,0.1);
    
    FitParameter  C_2("C_2",1,0.,0.1);
    FitParameter  D_2("D_2",1,0.,0.1);
    FitParameter  D_bar_2("D_bar_2",1,0.,0.1);
    FitParameter  S_2("S_2",1,0.,0.1);
    FitParameter  S_bar_2("S_bar_2",1,0.,0.1);
    
    FitParameter  C_3("C_3",1,0.,0.1);
    FitParameter  D_3("D_3",1,0.,0.1);
    FitParameter  D_bar_3("D_bar_3",1,0.,0.1);
    FitParameter  S_3("S_3",1,0.,0.1);
    FitParameter  S_bar_3("S_bar_3",1,0.,0.1);
    
    FitParameter  Gamma("Gamma",2,0.6629,0.0018);
    FitParameter  dGamma("dGamma",2,-0.088,0.006);
    FitParameter  dm("dm",2,17.757,0.1);
    
    FitParameter  offset_mean_dt("offset_mean_dt",1,0,0.1);
    FitParameter  scale_mean_dt("scale_mean_dt",1,0,0.1);
    FitParameter  scale_mean_2_dt("scale_mean_2_dt",1,0,0.1);
    FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
    FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.,0.1);
    FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
    FitParameter  offset_sigma2_dt("offset_sigma2_dt",1,0.,0.1);
    FitParameter  scale_sigma2_dt("scale_sigma2_dt",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt("scale_sigma2_2_dt",1,0.,0.1);
    FitParameter  offset_sigma3_dt("offset_sigma3_dt",1,0.,0.1);
    FitParameter  scale_sigma3_dt("scale_sigma3_dt",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt("scale_sigma3_2_dt",1,0.,0.1);
    FitParameter  offset_f_dt("offset_f_dt",1,1,0.1);
    FitParameter  scale_f_dt("scale_f_dt",1,0.,0.1);
    FitParameter  scale_f_2_dt("scale_f_2_dt",1,0.,0.1);
    FitParameter  offset_f2_dt("offset_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_dt("scale_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_2_dt("scale_f2_2_dt",1,0.,0.1);

    FitParameter  p0_os("p0_os",1,0.,0.);
    FitParameter  p1_os("p1_os",1,1.,0.);
    FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
    FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
    FitParameter  avg_eta_os("avg_eta_os",1,0.4,0.);
    FitParameter  tageff_os("tageff_os",1,0.4,0.);
    FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
    FitParameter  p0_ss("p0_ss",1,0.,0.);
    FitParameter  p1_ss("p1_ss",1,1.,0.);
    FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
    FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
    FitParameter  avg_eta_ss("avg_eta_ss",1,0.4,0.);
    FitParameter  tageff_ss("tageff_ss",1,0.7,0.);
    FitParameter  tageff_asym_ss("tageff_asym_ss",1,0.,0.);
    FitParameter  production_asym("production_asym",1,0.,0.);
    FitParameter  detection_asym("detection_asym",1,0.,0.);
    
    FitParameter  c0("c0",1,1,0.1);
    FitParameter  c1("c1",1,1,0.1);
    FitParameter  c2("c2",1,1,0.1);
    FitParameter  c3("c3",1,1,0.1);
    FitParameter  c4("c4",1,1,0.1);
    FitParameter  c5("c5",1,1,0.1);
    FitParameter  c6("c6",1,1,0.1);
    FitParameter  c7("c7",1,1,0.1);
    FitParameter  c8("c8",1,1,0.1);
    FitParameter  c9("c9",1,1,0.1);
    
    /// Fit parameters per Run    
    FitParameter  offset_mean_dt_Run1("offset_mean_dt_Run1",1,0,0.1);
    FitParameter  scale_mean_dt_Run1("scale_mean_dt_Run1",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run1("scale_mean_2_dt_Run1",1,0,0.1);
    FitParameter  offset_sigma_dt_Run1("offset_sigma_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run1("scale_sigma_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run1("scale_sigma_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run1("offset_sigma2_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run1("scale_sigma2_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run1("scale_sigma2_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run1("offset_sigma3_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run1("scale_sigma3_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run1("scale_sigma3_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_f_dt_Run1("offset_f_dt_Run1",1,1,0.1);
    FitParameter  scale_f_dt_Run1("scale_f_dt_Run1",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run1("scale_f_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_f2_dt_Run1("offset_f2_dt_Run1",1,0.,0.1);
    FitParameter  scale_f2_dt_Run1("scale_f2_dt_Run1",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run1("scale_f2_2_dt_Run1",1,0.,0.1);

    FitParameter  p0_os_Run1("p0_os_Run1",1,0.,0.);
    FitParameter  p1_os_Run1("p1_os_Run1",1,0.,0.);
    FitParameter  delta_p0_os_Run1("delta_p0_os_Run1",1,0.,0.);
    FitParameter  delta_p1_os_Run1("delta_p1_os_Run1",1,0.,0.);
    FitParameter  avg_eta_os_Run1("avg_eta_os_Run1",1,0.,0.);
    FitParameter  tageff_os_Run1("tageff_os_Run1",1,0.,0.);
    FitParameter  tageff_asym_os_Run1("tageff_asym_os_Run1",1,0.,0.);
    FitParameter  p0_ss_Run1("p0_ss_Run1",1,0.,0.);
    FitParameter  p1_ss_Run1("p1_ss_Run1",1,0.,0.);
    FitParameter  delta_p0_ss_Run1("delta_p0_ss_Run1",1,0.,0.);
    FitParameter  delta_p1_ss_Run1("delta_p1_ss_Run1",1,0.,0.);
    FitParameter  avg_eta_ss_Run1("avg_eta_ss_Run1",1,0.,0.);
    FitParameter  tageff_ss_Run1("tageff_ss_Run1",1,0.,0.);
    FitParameter  tageff_asym_ss_Run1("tageff_asym_ss_Run1",1,0.,0.);
    FitParameter  production_asym_Run1("production_asym_Run1",1,0.,0.);
    FitParameter  detection_asym_Run1("detection_asym_Run1",1,0.,0.);
    
    FitParameter  offset_mean_dt_Run2("offset_mean_dt_Run2",1,0,0.1);
    FitParameter  scale_mean_dt_Run2("scale_mean_dt_Run2",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run2("scale_mean_2_dt_Run2",1,0,0.1);
    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run2("scale_sigma_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run2("offset_sigma2_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run2("scale_sigma2_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run2("scale_sigma2_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run2("offset_sigma3_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run2("scale_sigma3_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run2("scale_sigma3_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_f_dt_Run2("offset_f_dt_Run2",1,1,0.1);
    FitParameter  scale_f_dt_Run2("scale_f_dt_Run2",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run2("scale_f_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_f2_dt_Run2("offset_f2_dt_Run2",1,0.,0.1);
    FitParameter  scale_f2_dt_Run2("scale_f2_dt_Run2",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run2("scale_f2_2_dt_Run2",1,0.,0.1);

    FitParameter  p0_os_Run2("p0_os_Run2",1,0.,0.);
    FitParameter  p1_os_Run2("p1_os_Run2",1,0.,0.);
    FitParameter  delta_p0_os_Run2("delta_p0_os_Run2",1,0.,0.);
    FitParameter  delta_p1_os_Run2("delta_p1_os_Run2",1,0.,0.);
    FitParameter  avg_eta_os_Run2("avg_eta_os_Run2",1,0.,0.);
    FitParameter  tageff_os_Run2("tageff_os_Run2",1,0.,0.);
    FitParameter  tageff_asym_os_Run2("tageff_asym_os_Run2",1,0.,0.);
    FitParameter  p0_ss_Run2("p0_ss_Run2",1,0.,0.);
    FitParameter  p1_ss_Run2("p1_ss_Run2",1,0.,0.);
    FitParameter  delta_p0_ss_Run2("delta_p0_ss_Run2",1,0.,0.);
    FitParameter  delta_p1_ss_Run2("delta_p1_ss_Run2",1,0.,0.);
    FitParameter  avg_eta_ss_Run2("avg_eta_ss_Run2",1,0.,0.);
    FitParameter  tageff_ss_Run2("tageff_ss_Run2",1,0.,0.);
    FitParameter  tageff_asym_ss_Run2("tageff_asym_ss_Run2",1,0.,0.);
    FitParameter  production_asym_Run2("production_asym_Run2",1,0.,0.);
    FitParameter  detection_asym_Run2("detection_asym_Run2",1,0.,0.);

    FitParameter  offset_mean_dt_Run2_17("offset_mean_dt_Run2_17",1,0,0.1);
    FitParameter  scale_mean_dt_Run2_17("scale_mean_dt_Run2_17",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run2_17("scale_mean_2_dt_Run2_17",1,0,0.1);
    FitParameter  offset_sigma_dt_Run2_17("offset_sigma_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2_17("scale_sigma_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run2_17("scale_sigma_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run2_17("offset_sigma2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run2_17("scale_sigma2_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run2_17("scale_sigma2_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run2_17("offset_sigma3_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run2_17("scale_sigma3_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run2_17("scale_sigma3_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_f_dt_Run2_17("offset_f_dt_Run2_17",1,1,0.1);
    FitParameter  scale_f_dt_Run2_17("scale_f_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run2_17("scale_f_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_f2_dt_Run2_17("offset_f2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f2_dt_Run2_17("scale_f2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run2_17("scale_f2_2_dt_Run2_17",1,0.,0.1);

    FitParameter  offset_mean_dt_Run2_18("offset_mean_dt_Run2_18",1,0,0.1);
    FitParameter  scale_mean_dt_Run2_18("scale_mean_dt_Run2_18",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run2_18("scale_mean_2_dt_Run2_18",1,0,0.1);
    FitParameter  offset_sigma_dt_Run2_18("offset_sigma_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2_18("scale_sigma_dt_Run2_18",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run2_18("scale_sigma_2_dt_Run2_18",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run2_18("offset_sigma2_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run2_18("scale_sigma2_dt_Run2_18",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run2_18("scale_sigma2_2_dt_Run2_18",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run2_18("offset_sigma3_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run2_18("scale_sigma3_dt_Run2_18",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run2_18("scale_sigma3_2_dt_Run2_18",1,0.,0.1);
    FitParameter  offset_f_dt_Run2_18("offset_f_dt_Run2_18",1,1,0.1);
    FitParameter  scale_f_dt_Run2_18("scale_f_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run2_18("scale_f_2_dt_Run2_18",1,0.,0.1);
    FitParameter  offset_f2_dt_Run2_18("offset_f2_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_f2_dt_Run2_18("scale_f2_dt_Run2_18",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run2_18("scale_f2_2_dt_Run2_18",1,0.,0.1);

    /// Fit parameters per run and trigger cat
    FitParameter  c0_Run1_t0("c0_Run1_t0",1,1,0.1);
    FitParameter  c1_Run1_t0("c1_Run1_t0",1,1,0.1);
    FitParameter  c2_Run1_t0("c2_Run1_t0",1,1,0.1);
    FitParameter  c3_Run1_t0("c3_Run1_t0",1,1,0.1);
    FitParameter  c4_Run1_t0("c4_Run1_t0",1,1,0.1);
    FitParameter  c5_Run1_t0("c5_Run1_t0",1,1,0.1);
    FitParameter  c6_Run1_t0("c6_Run1_t0",1,1,0.1);
    FitParameter  c7_Run1_t0("c7_Run1_t0",1,1,0.1);
    FitParameter  c8_Run1_t0("c8_Run1_t0",1,1,0.1);
    FitParameter  c9_Run1_t0("c9_Run1_t0",1,1,0.1);
    
    FitParameter  c0_Run1_t1("c0_Run1_t1",1,1,0.1);
    FitParameter  c1_Run1_t1("c1_Run1_t1",1,1,0.1);
    FitParameter  c2_Run1_t1("c2_Run1_t1",1,1,0.1);
    FitParameter  c3_Run1_t1("c3_Run1_t1",1,1,0.1);
    FitParameter  c4_Run1_t1("c4_Run1_t1",1,1,0.1);
    FitParameter  c5_Run1_t1("c5_Run1_t1",1,1,0.1);
    FitParameter  c6_Run1_t1("c6_Run1_t1",1,1,0.1);
    FitParameter  c7_Run1_t1("c7_Run1_t1",1,1,0.1);
    FitParameter  c8_Run1_t1("c8_Run1_t1",1,1,0.1);
    FitParameter  c9_Run1_t1("c9_Run1_t1",1,1,0.1);
    
    FitParameter  c0_Run2_t0("c0_Run2_t0",1,1,0.1);
    FitParameter  c1_Run2_t0("c1_Run2_t0",1,1,0.1);
    FitParameter  c2_Run2_t0("c2_Run2_t0",1,1,0.1);
    FitParameter  c3_Run2_t0("c3_Run2_t0",1,1,0.1);
    FitParameter  c4_Run2_t0("c4_Run2_t0",1,1,0.1);
    FitParameter  c5_Run2_t0("c5_Run2_t0",1,1,0.1);
    FitParameter  c6_Run2_t0("c6_Run2_t0",1,1,0.1);
    FitParameter  c7_Run2_t0("c7_Run2_t0",1,1,0.1);
    FitParameter  c8_Run2_t0("c8_Run2_t0",1,1,0.1);
    FitParameter  c9_Run2_t0("c9_Run2_t0",1,1,0.1);
    
    FitParameter  c0_Run2_t1("c0_Run2_t1",1,1,0.1);
    FitParameter  c1_Run2_t1("c1_Run2_t1",1,1,0.1);
    FitParameter  c2_Run2_t1("c2_Run2_t1",1,1,0.1);
    FitParameter  c3_Run2_t1("c3_Run2_t1",1,1,0.1);
    FitParameter  c4_Run2_t1("c4_Run2_t1",1,1,0.1);
    FitParameter  c5_Run2_t1("c5_Run2_t1",1,1,0.1);
    FitParameter  c6_Run2_t1("c6_Run2_t1",1,1,0.1);
    FitParameter  c7_Run2_t1("c7_Run2_t1",1,1,0.1);
    FitParameter  c8_Run2_t1("c8_Run2_t1",1,1,0.1);
    FitParameter  c9_Run2_t1("c9_Run2_t1",1,1,0.1);

    FitParameter  c0_Run2_17_t0("c0_Run2_17_t0",1,1,0.1);
    FitParameter  c1_Run2_17_t0("c1_Run2_17_t0",1,1,0.1);
    FitParameter  c2_Run2_17_t0("c2_Run2_17_t0",1,1,0.1);
    FitParameter  c3_Run2_17_t0("c3_Run2_17_t0",1,1,0.1);
    FitParameter  c4_Run2_17_t0("c4_Run2_17_t0",1,1,0.1);
    FitParameter  c5_Run2_17_t0("c5_Run2_17_t0",1,1,0.1);
    FitParameter  c6_Run2_17_t0("c6_Run2_17_t0",1,1,0.1);
    FitParameter  c7_Run2_17_t0("c7_Run2_17_t0",1,1,0.1);
    FitParameter  c8_Run2_17_t0("c8_Run2_17_t0",1,1,0.1);
    FitParameter  c9_Run2_17_t0("c9_Run2_17_t0",1,1,0.1);
    
    FitParameter  c0_Run2_17_t1("c0_Run2_17_t1",1,1,0.1);
    FitParameter  c1_Run2_17_t1("c1_Run2_17_t1",1,1,0.1);
    FitParameter  c2_Run2_17_t1("c2_Run2_17_t1",1,1,0.1);
    FitParameter  c3_Run2_17_t1("c3_Run2_17_t1",1,1,0.1);
    FitParameter  c4_Run2_17_t1("c4_Run2_17_t1",1,1,0.1);
    FitParameter  c5_Run2_17_t1("c5_Run2_17_t1",1,1,0.1);
    FitParameter  c6_Run2_17_t1("c6_Run2_17_t1",1,1,0.1);
    FitParameter  c7_Run2_17_t1("c7_Run2_17_t1",1,1,0.1);
    FitParameter  c8_Run2_17_t1("c8_Run2_17_t1",1,1,0.1);
    FitParameter  c9_Run2_17_t1("c9_Run2_17_t1",1,1,0.1);

    FitParameter  c0_Run2_18_t0("c0_Run2_18_t0",1,1,0.1);
    FitParameter  c1_Run2_18_t0("c1_Run2_18_t0",1,1,0.1);
    FitParameter  c2_Run2_18_t0("c2_Run2_18_t0",1,1,0.1);
    FitParameter  c3_Run2_18_t0("c3_Run2_18_t0",1,1,0.1);
    FitParameter  c4_Run2_18_t0("c4_Run2_18_t0",1,1,0.1);
    FitParameter  c5_Run2_18_t0("c5_Run2_18_t0",1,1,0.1);
    FitParameter  c6_Run2_18_t0("c6_Run2_18_t0",1,1,0.1);
    FitParameter  c7_Run2_18_t0("c7_Run2_18_t0",1,1,0.1);
    FitParameter  c8_Run2_18_t0("c8_Run2_18_t0",1,1,0.1);
    FitParameter  c9_Run2_18_t0("c9_Run2_18_t0",1,1,0.1);
    
    FitParameter  c0_Run2_18_t1("c0_Run2_18_t1",1,1,0.1);
    FitParameter  c1_Run2_18_t1("c1_Run2_18_t1",1,1,0.1);
    FitParameter  c2_Run2_18_t1("c2_Run2_18_t1",1,1,0.1);
    FitParameter  c3_Run2_18_t1("c3_Run2_18_t1",1,1,0.1);
    FitParameter  c4_Run2_18_t1("c4_Run2_18_t1",1,1,0.1);
    FitParameter  c5_Run2_18_t1("c5_Run2_18_t1",1,1,0.1);
    FitParameter  c6_Run2_18_t1("c6_Run2_18_t1",1,1,0.1);
    FitParameter  c7_Run2_18_t1("c7_Run2_18_t1",1,1,0.1);
    FitParameter  c8_Run2_18_t1("c8_Run2_18_t1",1,1,0.1);
    FitParameter  c9_Run2_18_t1("c9_Run2_18_t1",1,1,0.1);

    //FullTimePdf_mod t_pdf(r,delta,gamma,k);
    string marginalPdfsPrefix = "comb";
    if(fitGenMC)marginalPdfsPrefix = "Uniform";
    FullTimePdf t_pdf(C, D, D_bar, S, S_bar, k,
                      Gamma, dGamma, dm
		      ,offset_mean_dt,scale_mean_dt,scale_mean_2_dt
                      ,offset_sigma_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,offset_sigma2_dt, scale_sigma2_dt, scale_sigma2_2_dt
                      ,offset_sigma3_dt, scale_sigma3_dt, scale_sigma3_2_dt
                      ,offset_f_dt, scale_f_dt, scale_f_2_dt
                      ,offset_f2_dt, scale_f2_dt, scale_f2_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, marginalPdfsPrefix );

    /// Simultaneous pdfs
    FullTimePdf t_pdf_Run1_t0(C, D, D_bar, S, S_bar, k,
                      Gamma, dGamma, dm
		      ,offset_mean_dt_Run1,scale_mean_dt_Run1,scale_mean_2_dt_Run1
                      ,offset_sigma_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,offset_sigma2_dt_Run1, scale_sigma2_dt_Run1, scale_sigma2_2_dt_Run1
                      ,offset_sigma3_dt_Run1, scale_sigma3_dt_Run1, scale_sigma3_2_dt_Run1
                      ,offset_f_dt_Run1, scale_f_dt_Run1, scale_f_2_dt_Run1
                      ,offset_f2_dt_Run1, scale_f2_dt_Run1, scale_f2_2_dt_Run1
                      ,c0_Run1_t0, c1_Run1_t0, c2_Run1_t0 ,c3_Run1_t0, c4_Run1_t0, c5_Run1_t0
                      ,c6_Run1_t0, c7_Run1_t0, c8_Run1_t0, c9_Run1_t0,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1_t0" );
    
    FullTimePdf t_pdf_Run1_t1(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run1,scale_mean_dt_Run1,scale_mean_2_dt_Run1
			      ,offset_sigma_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
			      ,offset_sigma2_dt_Run1, scale_sigma2_dt_Run1, scale_sigma2_2_dt_Run1
			      ,offset_sigma3_dt_Run1, scale_sigma3_dt_Run1, scale_sigma3_2_dt_Run1
			      ,offset_f_dt_Run1, scale_f_dt_Run1, scale_f_2_dt_Run1
			      ,offset_f2_dt_Run1, scale_f2_dt_Run1, scale_f2_2_dt_Run1
                              ,c0_Run1_t1, c1_Run1_t1, c2_Run1_t1 ,c3_Run1_t1, c4_Run1_t1, c5_Run1_t1
                              ,c6_Run1_t1, c7_Run1_t1, c8_Run1_t1, c9_Run1_t1,
                              p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                              avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                              p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                              avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                              production_asym_Run1, detection_asym_Run1, "Run1_t1" );
    
    FullTimePdf t_pdf_Run2_t0(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2,scale_mean_dt_Run2,scale_mean_2_dt_Run2
			      ,offset_sigma_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
			      ,offset_sigma2_dt_Run2, scale_sigma2_dt_Run2, scale_sigma2_2_dt_Run2
			      ,offset_sigma3_dt_Run2, scale_sigma3_dt_Run2, scale_sigma3_2_dt_Run2
			      ,offset_f_dt_Run2, scale_f_dt_Run2, scale_f_2_dt_Run2
			      ,offset_f2_dt_Run2, scale_f2_dt_Run2, scale_f2_2_dt_Run2
                              ,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
                              ,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_t0" );
    
    FullTimePdf t_pdf_Run2_t1(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2,scale_mean_dt_Run2,scale_mean_2_dt_Run2
			      ,offset_sigma_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
			      ,offset_sigma2_dt_Run2, scale_sigma2_dt_Run2, scale_sigma2_2_dt_Run2
			      ,offset_sigma3_dt_Run2, scale_sigma3_dt_Run2, scale_sigma3_2_dt_Run2
			      ,offset_f_dt_Run2, scale_f_dt_Run2, scale_f_2_dt_Run2
			      ,offset_f2_dt_Run2, scale_f2_dt_Run2, scale_f2_2_dt_Run2
                              ,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
                              ,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_t1" );


    FullTimePdf t_pdf_Run2_17_t0(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2_17,scale_mean_dt_Run2_17,scale_mean_2_dt_Run2_17
			      ,offset_sigma_dt_Run2_17, scale_sigma_dt_Run2_17, scale_sigma_2_dt_Run2_17
			      ,offset_sigma2_dt_Run2_17, scale_sigma2_dt_Run2_17, scale_sigma2_2_dt_Run2_17
			      ,offset_sigma3_dt_Run2_17, scale_sigma3_dt_Run2_17, scale_sigma3_2_dt_Run2_17
			      ,offset_f_dt_Run2_17, scale_f_dt_Run2_17, scale_f_2_dt_Run2_17
			      ,offset_f2_dt_Run2_17, scale_f2_dt_Run2_17, scale_f2_2_dt_Run2_17
                              ,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
                              ,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_17_t0" );


    FullTimePdf t_pdf_Run2_17_t1(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2_17,scale_mean_dt_Run2_17,scale_mean_2_dt_Run2_17
			      ,offset_sigma_dt_Run2_17, scale_sigma_dt_Run2_17, scale_sigma_2_dt_Run2_17
			      ,offset_sigma2_dt_Run2_17, scale_sigma2_dt_Run2_17, scale_sigma2_2_dt_Run2_17
			      ,offset_sigma3_dt_Run2_17, scale_sigma3_dt_Run2_17, scale_sigma3_2_dt_Run2_17
			      ,offset_f_dt_Run2_17, scale_f_dt_Run2_17, scale_f_2_dt_Run2_17
			      ,offset_f2_dt_Run2_17, scale_f2_dt_Run2_17, scale_f2_2_dt_Run2_17
                              ,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
                              ,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_17_t1" );
    

    FullTimePdf t_pdf_Run2_18_t0(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2_18,scale_mean_dt_Run2_18,scale_mean_2_dt_Run2_18
			      ,offset_sigma_dt_Run2_18, scale_sigma_dt_Run2_18, scale_sigma_2_dt_Run2_18
			      ,offset_sigma2_dt_Run2_18, scale_sigma2_dt_Run2_18, scale_sigma2_2_dt_Run2_18
			      ,offset_sigma3_dt_Run2_18, scale_sigma3_dt_Run2_18, scale_sigma3_2_dt_Run2_18
			      ,offset_f_dt_Run2_18, scale_f_dt_Run2_18, scale_f_2_dt_Run2_18
			      ,offset_f2_dt_Run2_18, scale_f2_dt_Run2_18, scale_f2_2_dt_Run2_18
                              ,c0_Run2_t0, c1_Run2_t0, c2_Run2_t0 ,c3_Run2_t0, c4_Run2_t0, c5_Run2_t0
                              ,c6_Run2_t0, c7_Run2_t0, c8_Run2_t0, c9_Run2_t0,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_18_t0" );


    FullTimePdf t_pdf_Run2_18_t1(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2_18,scale_mean_dt_Run2_18,scale_mean_2_dt_Run2_18
			      ,offset_sigma_dt_Run2_18, scale_sigma_dt_Run2_18, scale_sigma_2_dt_Run2_18
			      ,offset_sigma2_dt_Run2_18, scale_sigma2_dt_Run2_18, scale_sigma2_2_dt_Run2_18
			      ,offset_sigma3_dt_Run2_18, scale_sigma3_dt_Run2_18, scale_sigma3_2_dt_Run2_18
			      ,offset_f_dt_Run2_18, scale_f_dt_Run2_18, scale_f_2_dt_Run2_18
			      ,offset_f2_dt_Run2_18, scale_f2_dt_Run2_18, scale_f2_2_dt_Run2_18
                              ,c0_Run2_t1, c1_Run2_t1, c2_Run2_t1 ,c3_Run2_t1, c4_Run2_t1, c5_Run2_t1
                              ,c6_Run2_t1, c7_Run2_t1, c8_Run2_t1, c9_Run2_t1,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2_18_t1" );


    /// Randomize start vals
    MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
    if(randomizeStartVals && (string)channel=="signal"){
	double val = gRandom->Uniform(0,1);
	mps->getParPtr("C")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("C"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("D")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("D"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("D_bar")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("D_bar"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("S")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("S"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("S_bar")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("S_bar"))->setInit(val);	
    }

    if(doScanParam){
	for(int i=0; i < mps->size(); i++){
		if((string)scanParam==((FitParameter*)mps->getParPtr(i))->name()){
			double val = scanParam_min + step * (scanParam_max - scanParam_min)/scanParam_steps; 
			((FitParameter*)mps->getParPtr(i))->setInit(val);
			((FitParameter*)mps->getParPtr(i))->fix();		
			cout << "Setting " << ((FitParameter*)mps->getParPtr(i))->name() << "  to  " << val << endl; 
		}
	}
    } 


    /// Load data
    double t,dt,mB;
    int f;
    int q_OS;
    double Bs_ID,Ds_ID;
    int Bs_TRUEID,bkgCAT;
    Int_t q_SS;
    double eta_OS;
    Double_t eta_SS;
    double sw;
    int year,run,Ds_finalState,trigger;
    double K[4];
    double pip[4];
    double pim[4];
    double Ds_Kp[4],Ds_Km[4],Ds_pim[4],Ds[4];
    
    TChain* tree_norm;

    if(fitGenMC){
	tree_norm=new TChain("MCDecayTreeTuple/MCDecayTree");
	tree_norm->Add(((string)InputGenMCFile).c_str());
	tree_norm->SetBranchStatus("*",0);
	tree_norm->SetBranchStatus("*TAU*",1);
	tree_norm->SetBranchStatus("*ID*",1);
	tree_norm->SetBranchStatus("*P*",1);
	
	tree_norm->SetBranchAddress("B_s0_TRUETAU",&t);
// 	tree_norm->SetBranchAddress("D_sminus_TRUEID",&Ds_ID);
	tree_norm->SetBranchAddress("D_splus_TRUEID",&Ds_ID);
	tree_norm->SetBranchAddress("B_s0_TRUEID",&Bs_ID);

	tree_norm->SetBranchAddress("Kminus_TRUEP_X",&K[0]);
	tree_norm->SetBranchAddress("Kminus_TRUEP_Y",&K[1]);
	tree_norm->SetBranchAddress("Kminus_TRUEP_Z",&K[2]);
	tree_norm->SetBranchAddress("Kminus_TRUEP_E",&K[3]);
	tree_norm->SetBranchAddress("piplus_TRUEP_X",&pim[0]);
	tree_norm->SetBranchAddress("piplus_TRUEP_Y",&pim[1]);
	tree_norm->SetBranchAddress("piplus_TRUEP_Z",&pim[2]);
	tree_norm->SetBranchAddress("piplus_TRUEP_E",&pim[3]);    
	tree_norm->SetBranchAddress("piminus_TRUEP_X",&pip[0]);
	tree_norm->SetBranchAddress("piminus_TRUEP_Y",&pip[1]);
	tree_norm->SetBranchAddress("piminus_TRUEP_Z",&pip[2]);
	tree_norm->SetBranchAddress("piminus_TRUEP_E",&pip[3]);    
	tree_norm->SetBranchAddress("Kplus_TRUEP_X",&Ds_Kp[0]);
	tree_norm->SetBranchAddress("Kplus_TRUEP_Y",&Ds_Kp[1]);
	tree_norm->SetBranchAddress("Kplus_TRUEP_Z",&Ds_Kp[2]);
	tree_norm->SetBranchAddress("Kplus_TRUEP_E",&Ds_Kp[3]);
	tree_norm->SetBranchAddress("Kminus0_TRUEP_X",&Ds_Km[0]);
	tree_norm->SetBranchAddress("Kminus0_TRUEP_Y",&Ds_Km[1]);
	tree_norm->SetBranchAddress("Kminus0_TRUEP_Z",&Ds_Km[2]);
	tree_norm->SetBranchAddress("Kminus0_TRUEP_E",&Ds_Km[3]);
	tree_norm->SetBranchAddress("piplus0_TRUEP_X",&Ds_pim[0]);
	tree_norm->SetBranchAddress("piplus0_TRUEP_Y",&Ds_pim[1]);
	tree_norm->SetBranchAddress("piplus0_TRUEP_Z",&Ds_pim[2]);
	tree_norm->SetBranchAddress("piplus0_TRUEP_E",&Ds_pim[3]);
    }
    else {
	tree_norm=new TChain("DecayTree");
	if(mode == "fit" && doToyStudy == 1){
		if(addBkgToToys)tree_norm->Add(((string)OutputDir+"sw_toys_"+anythingToString((int)step)+".root").c_str());
		else tree_norm->Add(((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());
	}
	else tree_norm->Add(((string)InputFileName).c_str());
	tree_norm->SetBranchStatus("*",0);
	tree_norm->SetBranchStatus("N_Bs_sw*",1);
	tree_norm->SetBranchStatus("year",1);
	tree_norm->SetBranchStatus("*DEC",1);
	tree_norm->SetBranchStatus("*PROB",1);
	tree_norm->SetBranchStatus("*OS*",1);
	tree_norm->SetBranchStatus("*TAU*",1);
	tree_norm->SetBranchStatus("*ID*",1);
	tree_norm->SetBranchStatus("weight",1);
	tree_norm->SetBranchStatus("TriggerCat",1);
	tree_norm->SetBranchStatus("run",1);
	tree_norm->SetBranchStatus("*finalState*",1);
	tree_norm->SetBranchStatus("BsDTF_*P*",1);
	tree_norm->SetBranchStatus("Bs_DTF_MM",1);
	if(useTrueTagging)tree_norm->SetBranchStatus("Bs_TRUEID",1);
	if(useTrueTagging)tree_norm->SetBranchStatus("bkgCAT",1);
	
	tree_norm->SetBranchAddress("Bs_DTF_MM",&mB);
	if(useTrueTau)tree_norm->SetBranchAddress("Bs_TRUETAU",&t);
	else tree_norm->SetBranchAddress("Bs_BsDTF_TAU",&t);
	tree_norm->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
	if(useTrueTagging)tree_norm->SetBranchAddress("Bs_TRUEID",&Bs_TRUEID);
	if(useTrueTagging)tree_norm->SetBranchAddress("bkgCAT",&bkgCAT);
	tree_norm->SetBranchAddress("Ds_ID",&f);
	tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
	tree_norm->SetBranchAddress("OS_Combination_DEC",&q_OS);
	tree_norm->SetBranchAddress("OS_Combination_PROB",&eta_OS);
	tree_norm->SetBranchAddress("SS_Kaon_DEC",&q_SS);
	tree_norm->SetBranchAddress("SS_Kaon_PROB",&eta_SS);
	tree_norm->SetBranchAddress(((string)weightName).c_str(),&sw);
	tree_norm->SetBranchAddress("year",&year);
	tree_norm->SetBranchAddress("run",&run);
	tree_norm->SetBranchAddress("TriggerCat",&trigger);
	tree_norm->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
	tree_norm->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
	tree_norm->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]);
	tree_norm->SetBranchAddress("BsDTF_Kplus_PE",&K[3]);
	tree_norm->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
	tree_norm->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
	tree_norm->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]);
	tree_norm->SetBranchAddress("BsDTF_piplus_PE",&pip[3]);    
	tree_norm->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
	tree_norm->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
	tree_norm->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]);
	tree_norm->SetBranchAddress("BsDTF_piminus_PE",&pim[3]);    
	if(doToyStudy && mode == "fit"){
		tree_norm->SetBranchAddress("BsDTF_Ds_PX",&Ds[0]);
		tree_norm->SetBranchAddress("BsDTF_Ds_PY",&Ds[1]);
		tree_norm->SetBranchAddress("BsDTF_Ds_PZ",&Ds[2]);
		tree_norm->SetBranchAddress("BsDTF_Ds_PE",&Ds[3]);    
	}
	else {
		tree_norm->SetBranchAddress("BsDTF_Ds_Kplus_PX",&Ds_Kp[0]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kplus_PY",&Ds_Kp[1]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kplus_PZ",&Ds_Kp[2]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kplus_PE",&Ds_Kp[3]);    
		tree_norm->SetBranchAddress("BsDTF_Ds_Kminus_PX",&Ds_Km[0]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kminus_PY",&Ds_Km[1]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kminus_PZ",&Ds_Km[2]);
		tree_norm->SetBranchAddress("BsDTF_Ds_Kminus_PE",&Ds_Km[3]);
		tree_norm->SetBranchAddress("BsDTF_Ds_piminus_PX",&Ds_pim[0]);
		tree_norm->SetBranchAddress("BsDTF_Ds_piminus_PY",&Ds_pim[1]);
		tree_norm->SetBranchAddress("BsDTF_Ds_piminus_PZ",&Ds_pim[2]);
		tree_norm->SetBranchAddress("BsDTF_Ds_piminus_PE",&Ds_pim[3]);
	}
    }
    DalitzEventList eventList,eventList_f,eventList_f_bar;
    DalitzEventList eventList_Run1_t0,eventList_Run1_t1,eventList_Run2_t0,eventList_Run2_t1,eventList_Run2_17_t0,eventList_Run2_17_t1,eventList_Run2_18_t0,eventList_Run2_18_t1;

    TH2D* h_Kpi_pipi_bin1 = new TH2D("h_Kpi_pipi_bin1",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});#left[m(#pi^{+} #pi^{+})#right] (GeV/c^{2}); Events (a.u.)",40,0.6,1.2,40,0.2,1.2);
    TH2D* h_Kpi_pipi_bin2 = new TH2D("h_Kpi_pipi_bin2",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});#left[m(#pi^{+} #pi^{+})#right] (GeV/c^{2}); Events (a.u.)",40,0.6,1.2,40,0.2,1.2);
    TH2D* h_Kpi_pipi_bin3 = new TH2D("h_Kpi_pipi_bin3",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});#left[m(#pi^{+} #pi^{+})#right] (GeV/c^{2}); Events (a.u.)",40,0.6,1.2,40,0.2,1.2);
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    RooRealVar* r_t = new RooRealVar("t", "t", min_TAU, max_TAU);
    RooRealVar* r_dt = new RooRealVar("dt", "dt",min_TAUERR, max_TAUERR);
    RooRealVar* r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5);
    RooRealVar* r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5);

    RooRealVar* r_mistag = new RooRealVar("mistag", "mistag",0., 0.5);

    RooCategory* r_f = new RooCategory("qf", "qf");
    r_f->defineType("h+", +1);
    r_f->defineType("h-", -1);

    RooCategory* r_q = new RooCategory("qt", "qt");
    r_q->defineType("B+", +1);
    r_q->defineType("B-", -1) ;   
    r_q->defineType("untagged", 0);    

    RooCategory* r_q_OS = new RooCategory("q_OS", "q_OS");
    r_q_OS->defineType("B+", +1);
    r_q_OS->defineType("B-", -1) ;   
    r_q_OS->defineType("untagged", 0);    

    RooCategory* r_q_SS = new RooCategory("q_SS", "q_SS");
    r_q_SS->defineType("B+", +1);
    r_q_SS->defineType("B-", -1) ;   
    r_q_SS->defineType("untagged", 0);    

    RooDataSet* data = new RooDataSet("data","data",RooArgSet(*r_t,*r_dt,*r_q,*r_mistag,*r_f));

    int N_sample = tree_norm->GetEntries();
    if(N_bootstrap == -1)N_bootstrap = N_sample;

    vector<int> b_indices;
    while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(ranLux.Uniform(0,N_sample-1)));
    sort(b_indices.begin(), b_indices.end());
    if(doBootstrap)N_sample = b_indices.size();

    TRandom3 rndm;
    TRandom3 randRes(seed);

    for(int i=0; i< N_sample; i++)
    {	
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << N_sample << endl;
	if(doBootstrap) tree_norm->GetEntry(b_indices[i]);
	else tree_norm->GetEntry(i);
        
	if(useTrueTau)t= t*1000.;//+ gRandom->Gaus(0.,offset_sigma_dt_Run2_17+scale_sigma_dt_Run2_17*dt);

	//t = -1.*(-1.85772e-04 + t*9.54051e-05)+t;
	//dt = 3.76166e-02 + 2.16574e-03 * t;

	if(fitGenMC){
			if(Ds_ID<0)f=-1;
	        	else if(Ds_ID > 0)f= 1;
	}
	double sign = 1.;
	//if(f > 0) sign = -1.;
	TLorentzVector K_p(sign*K[0],sign*K[1],sign*K[2],K[3]);
	TLorentzVector pip_p(sign*pip[0],sign*pip[1],sign*pip[2],pip[3]);
	TLorentzVector pim_p(sign*pim[0],sign*pim[1],sign*pim[2],pim[3]);
	TLorentzVector D_p;
	if(doToyStudy && mode == "fit"){
		D_p = TLorentzVector(sign*Ds[0],sign*Ds[1],sign*Ds[2],Ds[3]);
	}
	else {
		TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
		TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
		TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
		D_p = D_Kp_p + D_Km_p + D_pim_p;
	}
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
	// array of vectors
	vector<TLorentzVector> vectorOfvectors;

	if((string)channel=="norm"){
		TLorentzVector pip1_p, pip2_p;
		if(rndm.Rndm()<0.5) {
			pip1_p = K_p;
			pip2_p = pip_p;
		}
		else {
			pip1_p = pip_p;
			pip2_p = K_p;
		}
		K_p = pip1_p;
		pip_p = pip2_p;
	}

	vectorOfvectors.push_back(B_p*MeV);
	vectorOfvectors.push_back(D_p*MeV);
	vectorOfvectors.push_back(K_p*MeV);
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);
	DalitzEvent evt;

	if(f < 0)evt = DalitzEvent(pat, vectorOfvectors);
	else evt = DalitzEvent(pat_CP, vectorOfvectors);

	if(fitGenMC){
// 		if(i>4000)break;
		if(i<5000)continue;

		//DalitzEvent evt;
        	//if(Ds_ID<0)f=1;
        	//else if(Ds_ID > 0)f= -1;
		//if(f < 0)evt = DalitzEvent(pat);
		//else evt = DalitzEvent(pat_CP);
		
		t = t*1000.+randRes.Gaus(0.,offset_sigma_dt);
		if(i<5020)cout << "t = " << t << endl;

                if(t < min_TAU || t > max_TAU )continue;
		if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) continue;

   		evt.setValueInVector(0, t);
        	evt.setValueInVector(1, 0.);
        	evt.setValueInVector(2, -f);
        	int q = 0;
		if(Bs_ID>0)q=1;
		else q = -1;
	        evt.setValueInVector(3, q);
        	evt.setValueInVector(4, 0.);
        	evt.setValueInVector(5, q);
        	evt.setValueInVector(6, 0.);
        	evt.setValueInVector(7, 1);
        	evt.setValueInVector(8, 0);
                eventList.Add(evt);

                r_t->setVal(t);
                r_dt->setVal(0.04);
                r_q->setIndex(q);
                r_mistag->setVal(0.);
                r_f->setIndex(evt.getValueFromVector(2));
                data->add(RooArgSet(*r_t,*r_dt,*r_q,*r_mistag,*r_f));
		
		if(eventList.size()==5000)break;
		continue;
	}

	//if(Ds_finalState>0)continue;

	if(useTrueTagging){
		if(bkgCAT==1)continue;
		if(Bs_TRUEID == 531){ 
			q_OS=1;
			q_SS=1;
		}
		else if(Bs_TRUEID == -531) {
		 	q_OS = -1;
		 	q_SS = -1;
		}
		else continue;
		eta_OS = 0.;
		eta_SS = 0.;

 		//if( abs(mB-5366.770) > 15)continue;
		//t = t/mB*5366.770;
	}
        if(t < min_TAU || t > max_TAU )continue;
        if( dt < min_TAUERR || dt > max_TAUERR )continue;
        if(year < min_year || year > max_year) continue;
        if(trigger < min_trigger || trigger > max_trigger) continue;
        //if(Ds_finalState < min_Ds_finalState || Ds_finalState > max_Ds_finalState) continue;
        if(eta_SS < w_min || eta_SS > w_max )continue;
 	if((string)channel=="signal")if(!(evt.phaseSpace() > 0.))continue;

    	if(!(doToyStudy && mode == "fit"))run = (run == 2 && year == 17) ? 3 : run;
    	if(!(doToyStudy && mode == "fit"))run = (run == 2 && year == 18) ? 4 : run;
        evt.setWeight(sw);
        evt.setValueInVector(0, t);
        evt.setValueInVector(1, dt);   
        if(f<0)evt.setValueInVector(2, 1);
        else if(f > 0)evt.setValueInVector(2, -1);
        else {
            cout << "ERROR:: Undefined final state " << f << endl;  
            throw "ERROR";
        }
        evt.setValueInVector(3, q_OS);
        evt.setValueInVector(4, eta_OS);
        evt.setValueInVector(5, q_SS);
        evt.setValueInVector(6, eta_SS);
        evt.setValueInVector(7, run);
        evt.setValueInVector(8, trigger);
        evt.setValueInVector(9, mB);
        eventList.Add(evt);

        if(evt.getValueFromVector(2) == 1)eventList_f.Add(evt);
        else eventList_f_bar.Add(evt);
        
        if(run == 1 && trigger == 0) eventList_Run1_t0.Add(evt);
        else if(run == 1 && trigger == 1) eventList_Run1_t1.Add(evt);
        else if(run == 2 && trigger == 0) eventList_Run2_t0.Add(evt);
        else if(run == 2 && trigger == 1) eventList_Run2_t1.Add(evt);
        else if(run == 3 && trigger == 0) eventList_Run2_17_t0.Add(evt);
        else if(run == 3 && trigger == 1) eventList_Run2_17_t1.Add(evt);
        else if(run == 4 && trigger == 0) eventList_Run2_18_t0.Add(evt);
        else if(run == 4 && trigger == 1) eventList_Run2_18_t1.Add(evt);
    }
    
    /// Fit with MINT Pdf
    Neg2LL neg2LL(t_pdf, eventList);    

    Neg2LL neg2LL_Run1_t0(t_pdf_Run1_t0, eventList_Run1_t0);    
    Neg2LL neg2LL_Run1_t1(t_pdf_Run1_t1, eventList_Run1_t1);    
    Neg2LL neg2LL_Run2_t0(t_pdf_Run2_t0, eventList_Run2_t0);    
    Neg2LL neg2LL_Run2_t1(t_pdf_Run2_t1, eventList_Run2_t1);    
    Neg2LL neg2LL_Run2_17_t0(t_pdf_Run2_17_t0, eventList_Run2_17_t0);    
    Neg2LL neg2LL_Run2_17_t1(t_pdf_Run2_17_t1, eventList_Run2_17_t1);    
    Neg2LL neg2LL_Run2_18_t0(t_pdf_Run2_18_t0, eventList_Run2_18_t0);    
    Neg2LL neg2LL_Run2_18_t1(t_pdf_Run2_18_t1, eventList_Run2_18_t1);    

    Neg2LLSum neg2LL_sim;
    if(eventList_Run1_t0.size()>0)neg2LL_sim.add(&neg2LL_Run1_t0);
    if(eventList_Run1_t1.size()>0)neg2LL_sim.add(&neg2LL_Run1_t1);
    if(eventList_Run2_t0.size()>0)neg2LL_sim.add(&neg2LL_Run2_t0);
    if(eventList_Run2_t1.size()>0)neg2LL_sim.add(&neg2LL_Run2_t1);
    if(eventList_Run2_17_t0.size()>0)neg2LL_sim.add(&neg2LL_Run2_17_t0);
    if(eventList_Run2_17_t1.size()>0)neg2LL_sim.add(&neg2LL_Run2_17_t1);
    if(eventList_Run2_18_t0.size()>0)neg2LL_sim.add(&neg2LL_Run2_18_t0);
    if(eventList_Run2_18_t1.size()>0)neg2LL_sim.add(&neg2LL_Run2_18_t1);

     Neg2LLMultiConstraint constrains_tagging_Run1(MinuitParameterSet::getDefaultSet(),"_Tagging_Run1");
     Neg2LLMultiConstraint constrains_tagging_Run2(MinuitParameterSet::getDefaultSet(),"_Tagging_Run2");
     if(useGaussConstrainsTagging){
 	    neg2LL_sim.add(&constrains_tagging_Run1);
 	    neg2LL_sim.add(&constrains_tagging_Run2);
     }

    Neg2LLMultiConstraint constrains_Acc(MinuitParameterSet::getDefaultSet(),"_Acc");
    if(doAccSystematics && mode == "fit"){
	if(useCholDec){
		if(chol_index > constrains_Acc.getNumberParams()-1){
			cout << "ERROR:: Invalid cholesky index ! "<< endl;
			throw "ERROR";
		}
		constrains_Acc.smearInputValuesChol(chol_index,(step-1) - chol_index * varPerParChol, 1);
	}
	else
	{ 
		constrains_Acc.smearInputValues();
	}
    }

    Neg2LLMultiConstraint constrains_sys(MinuitParameterSet::getDefaultSet(),("_" + (string)doSystematic).c_str());
    if((string)doSystematic != "" && mode == "fit")constrains_sys.smearInputValues();


    Neg2LLMultiConstraint constrains_bias(MinuitParameterSet::getDefaultSet(),"_bias");
    if(useGaussConstrainsBias){
 	    neg2LL_sim.add(&constrains_bias);
     }

    Minimiser mini;
    if(doSimFit)mini.attachFunction(&neg2LL_sim);
    else mini.attachFunction(&neg2LL);
    if(mode == "fit"){
	mini.doFit();
    	mini.printResultVsInput();
//  	return;
    }

    /// Plot
    TCanvas* c = new TCanvas();
    h_Kpi_pipi_bin1->SetMarkerSize(0.1);
    h_Kpi_pipi_bin1->SetMarkerColor(kRed);
    h_Kpi_pipi_bin1->Draw();
    h_Kpi_pipi_bin2->SetMarkerSize(0.1);
    h_Kpi_pipi_bin2->SetMarkerColor(kBlue);
    h_Kpi_pipi_bin2->Draw("same");
    h_Kpi_pipi_bin3->SetMarkerSize(0.1);
    h_Kpi_pipi_bin3->SetMarkerColor(kBlack);
    h_Kpi_pipi_bin3->Draw("same");
    c->Print(((string)OutputDir+"h_Kpi_pipi.eps").c_str());
        
    TH1D* h_t = new TH1D("h_t",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed = new TH1D("h_t_mixed",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed = new TH1D("h_t_unmixed",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged = new TH1D("h_t_untagegged",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_OS = new TH1D("h_t_OS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_OS = new TH1D("h_t_mixed_OS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_OS = new TH1D("h_t_unmixed_OS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_SS = new TH1D("h_t_SS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_SS = new TH1D("h_t_mixed_SS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_SS = new TH1D("h_t_unmixed_SS",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_mp = new TH1D("h_t_mp",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0p = new TH1D("h_t_0p",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pp = new TH1D("h_t_pp",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_mm = new TH1D("h_t_mm",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0m = new TH1D("h_t_0m",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pm = new TH1D("h_t_pm",";t [ps];Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_p = new TH1D("h_t_p",";t [ps];Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_m = new TH1D("h_t_m",";t [ps];Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_N = new TH1D("h_t_N",";t [ps];Yield  (a.u.) ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_t_Nbar = new TH1D("h_t_Nbar",";t [ps];Yield  (a.u.) ",nBinsAsym,0.,2.*pi/dm);

    double mint = 0.;//min_TAU;
    double binsl[5] = {0.3,1.,2.,5.,max_TAU};
//     h_t_p->SetBins(4,binsl);
//     h_t_m->SetBins(4,binsl);

    TH1D* h_dt = new TH1D("h_dt",";#sigma_{t} [ps];Yield  (a.u.) ",nBinst,0,0.15);
    TH1D* h_eta_OS = new TH1D("h_eta_OS",";#eta_{OS};Yield  (a.u.) ",nBinst,0,0.5);
    TH1D* h_eta_SS = new TH1D("h_eta_SS",";#eta_{SS};Yield  (a.u.) ",nBinst,0,0.5);
    TH1D* h_q_OS = new TH1D("h_q_OS",";q_{OS};Yield  (a.u.) ",3,-1.5,1.5);
    TH1D* h_q_SS = new TH1D("h_q_SS",";q_{SS};Yield  (a.u.) ",3,-1.5,1.5);
    TH1D* h_f = new TH1D("h_f",";q_{f};Yield  (a.u.) ",2,-2,2);

    TH1D* h_N_mixed = new TH1D("h_N_mixed",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed = (TH1D*) h_N_mixed->Clone("h_N_unmixed");
    TH1D* h_N_mixed_OS = new TH1D("h_N_mixed_OS",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_OS = (TH1D*) h_N_mixed->Clone("h_N_unmixed_OS");
    TH1D* h_N_mixed_SS = new TH1D("h_N_mixed_SS",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_SS = (TH1D*) h_N_mixed->Clone("h_N_unmixed_SS");

    TH1D* h_N_mixed_p = (TH1D*) h_N_mixed->Clone("h_N_mixed_p");
    TH1D* h_N_unmixed_p = (TH1D*) h_N_mixed->Clone("h_N_unmixed_p");
    TH1D* h_N_mixed_m = (TH1D*) h_N_mixed->Clone("h_N_mixed_m");
    TH1D* h_N_unmixed_m = (TH1D*) h_N_mixed->Clone("h_N_unmixed_m");

    TH1D* h_N_mixed_p_unfolded = new TH1D("h_N_mixed_p_unfolded",";t [ps];A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_m_fit");

    double N_OS = 0;
    double N_SS = 0;
    double N_OS_SS = 0;
    double N = 0;

    double w_OS = 0;
    double w_SS = 0;
    double w_OS_SS = 0;
 	
    double D_OS = 0;
    double D_SS = 0;
    double D_OS_SS = 0;
    double D_comb = 0;

    double N_OS_all = 0;
    double N_SS_all = 0;
    double w_OS_all = 0;
    double w_SS_all = 0;
    double D_OS_all = 0;
    double D_SS_all = 0;
    
    double N_Run1_t0 = 0;
    double N_Run1_t1 = 0;
    double N_Run2_t0 = 0;
    double N_Run2_t1 = 0;
    double N_Run2_17_t0 = 0;
    double N_Run2_17_t1 = 0;
    double N_Run2_18_t0 = 0;
    double N_Run2_18_t1 = 0;

    double sigma_t_eff = 0;

    for (unsigned int i=0; i<eventList.size(); i++) {
 
        N += eventList[i].getWeight();
	if(!doSimFit)sigma_t_eff += t_pdf.getCalibratedResolution(eventList[i].getValueFromVector(1)) * eventList[i].getWeight();

        h_t->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        h_dt->Fill(eventList[i].getValueFromVector(1),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(3) != 0)h_eta_OS->Fill(eventList[i].getValueFromVector(4),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(5) != 0)h_eta_SS->Fill(eventList[i].getValueFromVector(6),eventList[i].getWeight());

        int f_evt = eventList[i].getValueFromVector(2);
        int q1 = eventList[i].getValueFromVector(3);
        int q2 = eventList[i].getValueFromVector(5);   
        int q_eff = 0;
        double w_eff = 0.5;
        int run_evt = eventList[i].getValueFromVector(7);   
        int trigger_evt = eventList[i].getValueFromVector(8);   
        
        h_q_OS->Fill(q1,eventList[i].getWeight());
        h_q_SS->Fill(q2,eventList[i].getWeight());
        h_f->Fill(f_evt,eventList[i].getWeight());

	if( q1 != 0)h_t_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
	if( q2 != 0)h_t_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());

        if(run_evt==1 && trigger_evt == 0)N_Run1_t0 += eventList[i].getWeight();
        else if(run_evt==1 && trigger_evt == 1)N_Run1_t1 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 0)N_Run2_t0 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 1)N_Run2_t1 += eventList[i].getWeight();
        else if(run_evt==3 && trigger_evt == 0)N_Run2_17_t0 += eventList[i].getWeight();
        else if(run_evt==3 && trigger_evt == 1)N_Run2_17_t1 += eventList[i].getWeight();
        else if(run_evt==4 && trigger_evt == 0)N_Run2_18_t0 += eventList[i].getWeight();
        else if(run_evt==4 && trigger_evt == 1)N_Run2_18_t1 += eventList[i].getWeight();

        std::pair<double, double> calibrated_mistag_os;
        std::pair<double, double> calibrated_mistag_ss;
        if(doSimFit){
            if(run_evt==1){
                calibrated_mistag_os = t_pdf_Run1_t0.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = t_pdf_Run1_t0.getCalibratedMistag_SS(eventList[i]);
            }
            else{
                calibrated_mistag_os = t_pdf_Run2_t0.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = t_pdf_Run2_t0.getCalibratedMistag_SS(eventList[i]);                
            }
        }
        else{
            calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eventList[i]);
            calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eventList[i]);        
        }

        double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
        double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
            
        if( p/(p+p_bar) > 0.5 ){ 
            q_eff = 1;
            w_eff = 1-p/(p+p_bar);
        }
        else if( p/(p+p_bar) < 0.5 ){
            q_eff = -1;
            w_eff = p/(p+p_bar);
        }

        if(q1 != 0 && q2 != 0){
            if(q_eff != 0){
                N_OS_SS += eventList[i].getWeight();
                w_OS_SS += w_eff * eventList[i].getWeight();
                D_OS_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
            }
        }
        else if( q1 != 0){
                //q_eff = q1;  flip tag ???
                N_OS += eventList[i].getWeight();
                w_OS += w_eff * eventList[i].getWeight(); 
                D_OS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
        }
        else if( q2 != 0){
                //q_eff = q2;
                N_SS += eventList[i].getWeight();
	    	w_SS += w_eff * eventList[i].getWeight(); 
                D_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
        } 

        D_comb += pow(1.-2.*w_eff,2)* eventList[i].getWeight();

	double D_res = 1.;
        if(doSimFit){
            if(run_evt==1){
		D_res = exp(-pow(t_pdf_Run1_t0.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
            else if(run_evt==2){
		D_res = exp(-pow(t_pdf_Run2_t0.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
            else if(run_evt==3){
		D_res = exp(-pow(t_pdf_Run2_17_t0.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
            else if(run_evt==4){
		D_res = exp(-pow(t_pdf_Run2_18_t0.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
	}

	double D_tag = 0.;
	if(q_eff != 0) D_tag = (1.-2.*(w_eff));
	double D_tot = D_tag * D_res;
	if(!dilutionWeight)D_tot = 1.;
	//cout << D_tot << "; t= " << eventList[i].getValueFromVector(1) << ";" << D_tag << ";" << D_res << endl;
// 	double D_tot = D_res;
        
        if(q1 != 0) N_OS_all += eventList[i].getWeight();
        if(q2 != 0){
                //if(q2 > 0 && calibrated_mistag_ss.first < 0.5) 
		N_SS_all += eventList[i].getWeight();
                //else if(q2 < 0 && calibrated_mistag_ss.second < 0.5) N_SS_all += eventList[i].getWeight();
        }    
            
        if(q1>0){
			w_OS_all +=  calibrated_mistag_os.first * eventList[i].getWeight();
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.first,2)* eventList[i].getWeight();
        } 
        else if(q1<0){
			w_OS_all +=  calibrated_mistag_os.second * eventList[i].getWeight();	
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.second,2)* eventList[i].getWeight();
        }

        if(q2>0){
			w_SS_all +=  calibrated_mistag_ss.first * eventList[i].getWeight();
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.first,2)* eventList[i].getWeight();
        } 
        else if(q2<0){
			w_SS_all +=  calibrated_mistag_ss.second * eventList[i].getWeight();	
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.second,2)* eventList[i].getWeight();
        }
   
        if((string)channel=="signal"){

	    if(f_evt == 1){ 
		h_t_p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
	    }
	    else h_t_m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
	
	    if(q_eff == 1)h_t_N->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);
	    else if(q_eff == -1)h_t_Nbar->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);

            if(q_eff==-1 && f_evt == 1){ 
			h_t_mp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(w_eff<w_max){
				h_N_mixed_p->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);
				h_N_mixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			}

            }
            else if(q_eff==0 && f_evt == 1)h_t_0p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
            else if(q_eff==1 && f_evt == 1){
                        h_t_pp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                            if(w_eff<w_max){
                                    h_N_unmixed_p->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);
                                    h_N_unmixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                            }
            }
            else if(q_eff==-1 && f_evt == -1){
                    h_t_mm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
        	    	if(w_eff<w_max){
                            h_N_unmixed_m->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);
                            h_N_unmixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    }
           }
           else if(q_eff==0 && f_evt == -1)h_t_0m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
           else if(q_eff==1 && f_evt == -1){
                    h_t_pm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    if(w_eff<w_max){
                        h_N_mixed_m->Fill(fmod(eventList[i].getValueFromVector(0)-mint,2.*pi/dm),eventList[i].getWeight()*D_tot);
                        h_N_mixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    }
                }
        }
//         else { 	
            if(q_eff == 0)h_t_untagegged->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
                else if(q_eff*f_evt > 0  ){
                        h_t_mixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        //if(w_eff<w_max)
			h_N_mixed->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);

                        if(q1 != 0)h_t_mixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(q1 != 0)h_N_mixed_OS->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);

                        if(q2 != 0)h_t_mixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(q2 != 0)h_N_mixed_SS->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
                    }
                else {
                    h_t_unmixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    //if(w_eff<w_max)
		    h_N_unmixed->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);

                    if(q1 != 0)h_t_unmixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
		    if(q1 != 0)h_N_unmixed_OS->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);

                    if(q2 != 0)h_t_unmixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
	            if(q2 != 0)h_N_unmixed_SS->Fill(fmod(eventList[i].getValueFromVector(0),2.*pi/dm),eventList[i].getWeight()*D_tot);
                }
//         }
  
    }     

    cout << "tree size = " << eventList.size() << endl;
    cout << "sumw = " << N << endl << endl;
    cout << "N_Run1_t0 =" << N_Run1_t0/N <<  endl;
    cout << "N_Run1_t1 =" << N_Run1_t1/N <<  endl;
    cout << "N_Run2_t0 =" << N_Run2_t0/N <<  endl;
    cout << "N_Run2_t1 =" << N_Run2_t1/N <<  endl;
    cout << "N_Run2_17_t0 =" << N_Run2_17_t0/N <<  endl;
    cout << "N_Run2_17_t1 =" << N_Run2_17_t1/N <<  endl;
    cout << "N_Run2_18_t0 =" << N_Run2_18_t0/N <<  endl;
    cout << "N_Run2_18_t1 =" << N_Run2_18_t1/N <<  endl;

    cout << "sigma_t_eff = " << sigma_t_eff/N << endl << endl;

    cout << "N_OS = " << N_OS_all << endl;
    cout << "N_SS = " << N_SS_all << endl << endl;

    cout << "eff_OS = " <<(N_OS_all)/N << endl;
    cout << "eff_SS = " <<(N_SS_all)/N << endl;
        
    cout << "Tagging perfromance " << endl << endl;        
    cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;

    cout << "OS  | " << (N_OS+N_OS_SS)/N << " | " <<  (w_OS_all)/(N_OS+N_OS_SS) << " | " << D_OS_all/N << endl;
    cout << "SS  | " << (N_SS+N_OS_SS)/N << " | " <<  (w_SS_all)/(N_SS+N_OS_SS) << " | " << D_SS_all/N << endl << endl;

    cout << "OS only  | " << N_OS/N << " | " <<  w_OS/N_OS << " | " << N_OS/N * D_OS/N_OS << endl;
    cout << "SS only  | " << N_SS/N << " | " <<  w_SS/N_SS << " | " << N_SS/N * D_SS/N_SS << endl;
    cout << "OS+SS    | " << N_OS_SS/N << " | " <<  w_OS_SS/N_OS_SS << " | " << N_OS_SS/N * D_OS_SS/N_OS_SS << endl;
    cout << "Combined | " << (N_OS+N_SS+N_OS_SS)/N << " | "<<  (w_OS+w_SS+w_OS_SS)/(N_OS+N_SS+N_OS_SS) << " | " << (N_OS+N_SS+N_OS_SS)/N * D_comb/(N_OS+N_SS+N_OS_SS) << endl << endl ;

    /// Generate toys 
    DalitzEventList toys;
    if(mode == "gen"){
	if(doSimFit) {
		toys.Add(t_pdf_Run1_t0.generateToys(N_scale_toys * N_Run1_t0,1,0));
		toys.Add(t_pdf_Run1_t1.generateToys(N_scale_toys *N_Run1_t1,1,1));
		toys.Add(t_pdf_Run2_t0.generateToys(N_scale_toys *N_Run2_t0,2,0));
		toys.Add(t_pdf_Run2_t1.generateToys(N_scale_toys *N_Run2_t1,2,1));
		toys.Add(t_pdf_Run2_17_t0.generateToys(N_scale_toys *N_Run2_17_t0,3,0));
		toys.Add(t_pdf_Run2_17_t1.generateToys(N_scale_toys *N_Run2_17_t1,3,1));
		toys.Add(t_pdf_Run2_18_t0.generateToys(N_scale_toys *N_Run2_18_t0,4,0));
		toys.Add(t_pdf_Run2_18_t1.generateToys(N_scale_toys *N_Run2_18_t1,4,1));
	}
	else  toys.Add(t_pdf.generateToys(N_scale_toys *N));
	
	if(addBkgToToys == 1){
			t_pdf.saveEventListToFile(toys,((string)OutputDir+"signal_toys_"+anythingToString((int)step)+".root").c_str());

 			string bkg_template_input = "signal_SS.root";

			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run1_t0/N,1,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run1_t1/N,1,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_t0/N,2,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_t1/N,2,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_17_t0/N,3,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_17_t1/N,3,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_18_t0/N,4,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToys(N_scale_toys * (eventList.size()-N) * N_Run2_18_t1/N,4,1,bkg_template_input));
	}
	if(addBkgToToys == 2){
			t_pdf.saveEventListToFile(toys,((string)OutputDir+"signal_toys_"+anythingToString((int)step)+".root").c_str());

			string bkg_template_input = "bkg_template_test2.root";

			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run1_t0/N,1,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run1_t1/N,1,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_t0/N,2,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_t1/N,2,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_17_t0/N,3,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_17_t1/N,3,1,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_18_t0/N,4,0,bkg_template_input));
			toys.Add(t_pdf.generateBkgToysBootstrap(N_scale_toys * (eventList.size()-N) * N_Run2_18_t1/N,4,1,bkg_template_input));
	}	

	
	t_pdf.saveEventListToFile(toys,((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());
	return;
    }

    /// Save results 
    TMatrixTSym<double> cov_full = mini.covMatrixFull();     
    //cov_full.Print();
	
    ofstream fitResults;
    fitResults.open(((string)OutputDir+"FitResults.txt").c_str(),std::ofstream::trunc);

    vector<string> prefix;
    if(doSimFit){
		prefix.push_back("_Run1");
    		prefix.push_back("_Run2");
    }
    else prefix.push_back("");

    for(int p = 0 ; p < prefix.size(); p++){
	
	fitResults << "\"" << "ConstrainMulti_Tagging" << prefix[p] << "\"" << "     " << "\"" ;
	for(int i = 0 ; i < mps->size(); i++){	
		if(mps->getParPtr(i)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(i)->name()))
			fitResults << mps->getParPtr(i)->name() << " ";
	}	
	fitResults << "\"";
	fitResults<< "\n";
	fitResults<< "\n";
	fitResults << "\"" << "ConstrainMulti_Tagging" << prefix[p] << "_corr" << "\"" << "     " << "\"" ;
	for(int i= 0; i< cov_full.GetNcols(); i++){
			if(!(mps->getParPtr(i)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(i)->name())))continue;
			for(int j= i; j< cov_full.GetNcols(); j++){
				if(!(mps->getParPtr(j)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(j)->name())))continue;
				if(i == j) fitResults << 1. << " " ;
				else if(cov_full[i][i] == 0. || cov_full[j][j] == 0.) fitResults << 0. << " " ;
				//else if(cov_full[i][j]/sqrt(cov_full[i][i])/sqrt(cov_full[j][j]) < 0.00001) fitResults << 0. << " " ; 
				else fitResults << (cov_full[i][j]/sqrt(cov_full[i][i])/sqrt(cov_full[j][j])) << " ";
			}
	}
	fitResults << "\"";
	fitResults<< "\n";
	fitResults<< "\n";
    }

     for(int i = 0 ; i < mps->size(); i++)	    
		fitResults << "\"" << mps->getParPtr(i)->name() << "\"" << "    " << mps->getParPtr(i)->iFixInit() << "    " << mps->getParPtr(i)->mean() << "    " << mps->getParPtr(i)->err() <<  "\n";
     fitResults.close();

    
    /// Save pulls
    gDirectory->cd();
    string paraFileName = (string)OutputDir+"pull_"+ (string)doSystematic+ "_" + anythingToString((int)step)+".root";
    if(doAccSystematics) paraFileName = (string)OutputDir+"pullAcc_"+anythingToString((int)step)+".root";
    if(doAccSystematics && useCholDec) paraFileName = (string)OutputDir+"pullAccChol_"+anythingToString((int)step)+".root";
    cout << paraFileName << endl;

    TFile* paraFile = new TFile( paraFileName.c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;
    mps->fillNtp(paraFile, ntp);
    ntp->AutoSave();

    TTree* pull_tree = new TTree("Coherence","Coherence");
    double n2ll;
    double C_val,D_val,S_val;
    double Cbar_val,Dbar_val,Sbar_val;
    double C_err,D_err,S_err;
    double Cbar_err,Dbar_err,Sbar_err;
   
    TBranch* br_C = pull_tree->Branch( "C", &C_val, "C_val/D" );
    TBranch* br_Cbar = pull_tree->Branch( "Cbar", &Cbar_val, "Cbar_val/D" );
    TBranch* br_D = pull_tree->Branch( "D", &D_val, "D_val/D" );
    TBranch* br_Dbar = pull_tree->Branch( "Dbar", &Dbar_val, "Dbar_val/D" );
    TBranch* br_S = pull_tree->Branch( "S", &S_val, "S_val/D" );
    TBranch* br_Sbar = pull_tree->Branch( "Sbar", &Sbar_val, "Sbar_val/D" );

    TBranch* br_C_err = pull_tree->Branch( "C_err", &C_err, "C_err/D" );
    TBranch* br_Cbar_err = pull_tree->Branch( "Cbar_err", &Cbar_err, "Cbar_err/D" );
    TBranch* br_D_err = pull_tree->Branch( "D_err", &D_err, "D_err/D" );
    TBranch* br_Dbar_err = pull_tree->Branch( "Dbar_err", &Dbar_err, "Dbar_err/D" );
    TBranch* br_S_err = pull_tree->Branch( "S_err", &S_err, "S_err/D" );
    TBranch* br_Sbar_err = pull_tree->Branch( "Sbar_err", &Sbar_err, "Sbar_err/D" );

    TBranch* br_n2ll = pull_tree->Branch( "n2ll_mean", &n2ll, "n2ll/D" );
    TBranch* br_seed = pull_tree->Branch( "seed", &seed, "seed/I" );
        
    if(doSimFit)n2ll = neg2LL_sim.getVal();
    else n2ll = neg2LL.getVal();

    C_val = C.blindedMean();
    C_err = C.err();
    Cbar_val = -C.blindedMean();
    Cbar_err = C.err();
    D_val = D.blindedMean();
    D_err = D.err();
    Dbar_val = D_bar.blindedMean();
    Dbar_err = D_bar.err();
    S_val = S.blindedMean();
    S_err = S.err();
    Sbar_val = S_bar.blindedMean();
    Sbar_err = S_bar.err();

    paraFile->cd();
    pull_tree->Fill();
    pull_tree->SetDirectory(paraFile); 
    pull_tree->Write();
    paraFile->Close();
    delete paraFile;

    /// Create tagging perfromance tables
    if(updateAnaNote){

	ofstream resultsfile;
	resultsfile.open(("../../../../../TD-AnaNote/latex/tables/timeFit/"+(string)OutputDir+"result.tex").c_str(),std::ofstream::trunc);
	resultsfile << "\\begin{table}[h]" << "\n";
	resultsfile << "\\centering" << "\n";
// 	resultsfile << "\\small" << "\n";
	resultsfile << "\\caption{Result of the phase-space integrated fit to "; 
	if((string)channel == "norm")resultsfile << "$B_s \\to D_s \\pi \\pi \\pi$";
	else if((string)channel == "signal")resultsfile << "$B_s \\to D_s K \\pi \\pi$";
	resultsfile << " data.}\n";
	resultsfile << "\\begin{tabular}{c c c}" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "& Fit parameter & Value \\\\" << "\n";
	resultsfile << "\\hline" << "\n";

	if((string)channel == "norm"){
		resultsfile << std::fixed << std::setprecision(4) 
		<< "Run-I & $p_{0}^{\\text{OS}}$ & " <<  mps->getParPtr("p0_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p0_os_Run1")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("p1_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p1_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p0_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p1_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_asym_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_os_Run1")->err() << "\\\\" << "\n"
		
		<< "& $p_{0}^{\\text{SS}}$ & " <<  mps->getParPtr("p0_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p0_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("p1_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p1_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p0_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p1_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_asym_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_ss_Run1")->err() << "\\\\" << "\n"
	
		<< "&$A_{p}$ & " <<  mps->getParPtr("production_asym_Run1")->mean() << " $\\pm$ " << mps->getParPtr("production_asym_Run1")->err() << "\\\\" << "\n";
	
		resultsfile << "\\\\" << "\n" ;
		resultsfile << std::fixed << std::setprecision(4) 
		<< "Run-II & $p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("p0_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p0_os_Run2")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("p1_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p1_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p0_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p1_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_asym_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_os_Run2")->err() << "\\\\" << "\n"
	
		<< "& $p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("p0_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p0_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("p1_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p1_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p0_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p1_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_asym_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_ss_Run2")->err() << "\\\\" << "\n"
	
		<< "&$A_{p}$ & " <<  mps->getParPtr("production_asym_Run2")->mean() << " $\\pm$ " << mps->getParPtr("production_asym_Run2")->err() << "\\\\" << "\n";
	
		resultsfile << "\\\\" << "\n" ;
		resultsfile << std::fixed << std::setprecision(4) 
		<< "&$\\Delta m_{s}$ & " 
		//<<   mps->getParPtr("dm")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("dm")->err() << "\\\\" << "\n";
	}

	else {
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $C$ & " 
	// 	<<   mps->getParPtr("C")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("C")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "&$D$ & " 
	// 	<<   mps->getParPtr("D")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("D")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "&$\\bar D$ & " 
	// 	<<   mps->getParPtr("D_bar")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("D_bar")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $S$ & " 
	// 	<<   mps->getParPtr("S")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("S")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $\\bar S$ & " 
	// 	<<   mps->getParPtr("S_bar")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("S_bar")->err() << "\\\\" << "\n";
	}

	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\end{tabular}" << "\n";
	resultsfile << "\\label{table:timeFit_" << (string) channel << "}" << "\n";
	resultsfile << "\\end{table}";	
	resultsfile.close();

	for(int p = 0 ; p < prefix.size(); p++){

		vector<string> cov_params;
//     		if(!mps->getParPtr("Gamma")->iFixInit())cov_params.push_back("Gamma");
//     		if(!mps->getParPtr("dGamma")->iFixInit())cov_params.push_back("dGamma");
//     		if(!mps->getParPtr("dm")->iFixInit())cov_params.push_back("dm");

//     		if(!mps->getParPtr("offset_sigma_dt"+prefix[p])->iFixInit())cov_params.push_back("offset_sigma_dt"+prefix[p]);
//     		if(!mps->getParPtr("scale_sigma_dt"+prefix[p])->iFixInit())cov_params.push_back("scale_sigma_dt"+prefix[p]);

		if(!mps->getParPtr("p0_os"+prefix[p])->iFixInit())cov_params.push_back("p0_os"+prefix[p]);
    		if(!mps->getParPtr("p1_os"+prefix[p])->iFixInit())cov_params.push_back("p1_os"+prefix[p]);
    		if(!mps->getParPtr("delta_p0_os"+prefix[p])->iFixInit())cov_params.push_back("delta_p0_os"+prefix[p]);
    		if(!mps->getParPtr("delta_p1_os"+prefix[p])->iFixInit())cov_params.push_back("delta_p1_os"+prefix[p]);
//     		if(!mps->getParPtr("avg_eta_os"+prefix[p])->iFixInit())cov_params.push_back("avg_eta_os"+prefix[p]);
    		if(!mps->getParPtr("tageff_os"+prefix[p])->iFixInit())cov_params.push_back("tageff_os"+prefix[p]);
    		//if(!mps->getParPtr("tageff_asym_os"+prefix[p])->iFixInit())cov_params.push_back("tageff_asym_os"+prefix[p]);

		if(!mps->getParPtr("p0_ss"+prefix[p])->iFixInit())cov_params.push_back("p0_ss"+prefix[p]);
    		if(!mps->getParPtr("p1_ss"+prefix[p])->iFixInit())cov_params.push_back("p1_ss"+prefix[p]);
    		if(!mps->getParPtr("delta_p0_ss"+prefix[p])->iFixInit())cov_params.push_back("delta_p0_ss"+prefix[p]);
    		if(!mps->getParPtr("delta_p1_ss"+prefix[p])->iFixInit())cov_params.push_back("delta_p1_ss"+prefix[p]);
//     		if(!mps->getParPtr("avg_eta_ss"+prefix[p])->iFixInit())cov_params.push_back("avg_eta_ss"+prefix[p]);
    		if(!mps->getParPtr("tageff_ss"+prefix[p])->iFixInit())cov_params.push_back("tageff_ss"+prefix[p]);
    		//if(!mps->getParPtr("tageff_asym_ss"+prefix[p])->iFixInit())cov_params.push_back("tageff_asym_ss"+prefix[p]);

//     		if(!mps->getParPtr("production_asym"+prefix[p])->iFixInit())cov_params.push_back("production_asym"+prefix[p]);
//     		if(!mps->getParPtr("detection_asym"+prefix[p])->iFixInit())cov_params.push_back("detection_asym"+prefix[p]);

		vector<int> cov_params_id;
		RooArgList xvec, mu;
		

		for(int i = 0; i < cov_params.size(); i++){
			cov_params_id.push_back(mps->findParPtr(cov_params[i]));
			double mean = mps->getParPtr(cov_params[i])->mean();	
			double error = mps->getParPtr(cov_params[i])->err();	
		
			RooRealVar* x = new RooRealVar(("x_"+cov_params[i]).c_str(), ("x_"+cov_params[i]).c_str(),mean-10.*error,mean+10.*error);
			xvec.add(*x);
			mu.add(RooRealConstant::value(mean));
		}
		
		if(cov_params.size()>0){

			TMatrixTSym<double> cov(cov_params.size());	
			for(int i = 0; i < cov_params_id.size(); i++)for(int j = 0; j < cov_params_id.size(); j++) cov[i][j] = cov_full[cov_params_id[i]][cov_params_id[j]]; 
				
			cov.Print();
			xvec.Print();
			mu.Print();
			RooMultiVarGaussian gauss_cov("gauss_cov","gauss_cov",xvec, mu, cov);
	
			const int N_toys_cov = 1000; 
			RooDataSet* data_cov = gauss_cov.generate(xvec, N_toys_cov);
	
			double N_tot = 0;
			vector<double> v_N_OS(N_toys_cov,0.);
			vector<double> v_N_SS(N_toys_cov,0.);
			vector<double> v_N_OS_SS(N_toys_cov,0.);
			
			vector<double> v_w_OS(N_toys_cov,0.);
			vector<double> v_w_SS(N_toys_cov,0.);
			vector<double> v_w_OS_SS(N_toys_cov,0.);
				
			vector<double> v_D_OS(N_toys_cov,0.);
			vector<double> v_D_SS(N_toys_cov,0.);
			vector<double> v_D_OS_SS(N_toys_cov,0.);
			vector<double> v_D_comb(N_toys_cov,0.);
			
			for (unsigned int i=0; i<eventList.size(); i++) {
			
				int f_evt = eventList[i].getValueFromVector(2);
				int q1 = eventList[i].getValueFromVector(3);
				int q2 = eventList[i].getValueFromVector(5);   
				int q_eff = 0;
				double w_eff = 0.5;
				int run_evt = eventList[i].getValueFromVector(7);   
				int trigger_evt = eventList[i].getValueFromVector(8);   
	
				if(A_is_in_B("Run1",prefix[p]) && run_evt > 1) continue;
				else if(A_is_in_B("Run2",prefix[p]) && run_evt == 1) continue;
	
				N_tot += eventList[i].getWeight();
	
				for(int j = 0 ; j < N_toys_cov; j++){
					RooArgSet* xvec_cov= (RooArgSet*)data_cov->get(j);
	
					double x_avg_eta_ss = xvec_cov->find(("x_avg_eta_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_avg_eta_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("avg_eta_ss"+prefix[p])->mean(); 
					double x_p0_ss = xvec_cov->find(("x_p0_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p0_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p0_ss"+prefix[p])->mean();
					double x_p1_ss = xvec_cov->find(("x_p1_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p1_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p1_ss"+prefix[p])->mean(); 
					double x_delta_p0_ss = xvec_cov->find(("x_delta_p0_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p0_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p0_ss"+prefix[p])->mean();
					double x_delta_p1_ss = xvec_cov->find(("x_delta_p1_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p1_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p1_ss"+prefix[p])->mean(); 
	
					double x_avg_eta_os = xvec_cov->find(("x_avg_eta_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_avg_eta_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("avg_eta_os"+prefix[p])->mean(); 
					double x_p0_os = xvec_cov->find(("x_p0_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p0_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p0_os"+prefix[p])->mean();
					double x_p1_os = xvec_cov->find(("x_p1_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p1_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p1_os"+prefix[p])->mean(); 
					double x_delta_p0_os = xvec_cov->find(("x_delta_p0_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p0_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p0_os"+prefix[p])->mean();
					double x_delta_p1_os = xvec_cov->find(("x_delta_p1_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p1_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p1_os"+prefix[p])->mean(); 
	
					std::pair<double, double> calibrated_mistag_os;
					std::pair<double, double> calibrated_mistag_ss;
		
					calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eventList[i],x_avg_eta_ss,x_p0_ss,x_p1_ss,x_delta_p0_ss,x_delta_p1_ss);
					calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eventList[i],x_avg_eta_os,x_p0_os,x_p1_os,x_delta_p0_os,x_delta_p1_os);    
		
					double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
					double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
					
					if( p/(p+p_bar) > 0.5 ){ 
						q_eff = 1;
						w_eff = 1-p/(p+p_bar);
					}
					else if( p/(p+p_bar) < 0.5 ){
						q_eff = -1;
						w_eff = p/(p+p_bar);
					}
				
					if(q1 != 0 && q2 != 0){
					if(q_eff != 0){
						v_N_OS_SS[j] += eventList[i].getWeight();
						v_w_OS_SS[j] += w_eff * eventList[i].getWeight();
						v_D_OS_SS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
					}
					}
					else if( q1 != 0){
						//q_eff = q1;  flip tag ???
						v_N_OS[j] += eventList[i].getWeight();
						v_w_OS[j] += w_eff * eventList[i].getWeight(); 
						v_D_OS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
					}
					else if( q2 != 0){
						//q_eff = q2;
						v_N_SS[j] += eventList[i].getWeight();
						v_w_SS[j] += w_eff * eventList[i].getWeight(); 
						v_D_SS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
					} 
				}
			}
	
			double eff_OS = 0;
			double eff_OS_err = 0;
			double w_OS = 0;
			double w_OS_err = 0;
			double e_OS = 0;
			double e_OS_err = 0;
	
			double eff_SS = 0;
			double eff_SS_err = 0;
			double w_SS = 0;
			double w_SS_err = 0;
			double e_SS = 0;
			double e_SS_err = 0;
	
			double eff_OS_SS = 0;
			double eff_OS_SS_err = 0;
			double w_OS_SS = 0;
			double w_OS_SS_err = 0;
			double e_OS_SS = 0;
			double e_OS_SS_err = 0;
		
			for(int j = 0 ; j < N_toys_cov; j++){
				eff_OS += (v_N_OS[j])/N_tot/N_toys_cov;
				w_OS += (v_w_OS[j])/(v_N_OS[j])/N_toys_cov;
				e_OS += (v_D_OS[j])/N_tot/N_toys_cov;
	
				eff_SS += (v_N_SS[j])/N_tot/N_toys_cov;
				w_SS += (v_w_SS[j])/(v_N_SS[j])/N_toys_cov;
				e_SS += (v_D_SS[j])/N_tot/N_toys_cov;
	
				eff_OS_SS += (v_N_OS_SS[j])/N_tot/N_toys_cov;
				w_OS_SS += (v_w_OS_SS[j])/(v_N_OS_SS[j])/N_toys_cov;
				e_OS_SS += (v_D_OS_SS[j])/N_tot/N_toys_cov;	
			}
			for(int j = 0 ; j < N_toys_cov; j++){
				eff_OS_err += pow(((v_N_OS[j])/N_tot-eff_OS),2)/(N_toys_cov-1.);
				w_OS_err += pow(((v_w_OS[j])/(v_N_OS[j]) - w_OS),2)/(N_toys_cov-1.);
				e_OS_err += pow((v_D_OS[j]/N_tot-e_OS),2)/(N_toys_cov-1.);
	
				eff_SS_err += pow(((v_N_SS[j])/N_tot-eff_SS),2)/(N_toys_cov-1.);
				w_SS_err += pow(((v_w_SS[j])/(v_N_SS[j]) - w_SS),2)/(N_toys_cov-1.);
				e_SS_err += pow((v_D_SS[j]/N_tot-e_SS),2)/(N_toys_cov-1.);
	
				eff_OS_SS_err += pow(((v_N_OS_SS[j])/N_tot-eff_OS_SS),2)/(N_toys_cov-1.);
				w_OS_SS_err += pow(((v_w_OS_SS[j])/(v_N_OS_SS[j]) - w_OS_SS),2)/(N_toys_cov-1.);
				e_OS_SS_err += pow((v_D_OS_SS[j]/N_tot-e_OS_SS),2)/(N_toys_cov-1.);
			}
			
			double rel_eff_err_OS = mps->getParPtr("tageff_os"+prefix[p])->iFixInit() ? 0. : mps->getParPtr("tageff_os"+prefix[p])->err() / mps->getParPtr("tageff_os"+prefix[p])->mean(); 
	
			eff_OS_err += pow( rel_eff_err_OS * eff_OS,2);
			w_OS_err += pow( rel_eff_err_OS * w_OS,2);
			e_OS_err += pow( rel_eff_err_OS * e_OS,2);
	
			double rel_eff_err_SS = mps->getParPtr("tageff_ss"+prefix[p])->iFixInit() ? 0. : mps->getParPtr("tageff_ss"+prefix[p])->err() / mps->getParPtr("tageff_ss"+prefix[p])->mean(); 
	
			eff_SS_err += pow( rel_eff_err_SS * eff_SS,2);
			w_SS_err += pow( rel_eff_err_SS * w_SS,2);
			e_SS_err += pow( rel_eff_err_SS * e_SS,2);
	
			eff_OS_SS_err += pow( rel_eff_err_OS * eff_OS_SS,2) + pow( rel_eff_err_SS * eff_OS_SS,2);
			w_OS_SS_err += pow( rel_eff_err_OS * w_OS_SS,2) + pow( rel_eff_err_SS * w_OS_SS,2);
			e_OS_SS_err += pow( rel_eff_err_OS * e_OS_SS,2) + pow( rel_eff_err_SS * e_OS_SS,2);
	
			double eff_tot = eff_OS + eff_SS + eff_OS_SS;
			double eff_tot_err = eff_OS_err + eff_SS_err + eff_OS_SS_err;
			double w_tot = (eff_OS * w_OS + eff_SS * w_SS + eff_OS_SS * w_OS_SS)/eff_tot;
			double w_tot_err = (eff_OS * w_OS_err + eff_SS * w_SS_err + eff_OS_SS * w_OS_SS_err)/eff_tot;
			double e_tot = e_OS + e_SS + e_OS_SS;
			double e_tot_err = e_OS_err + e_SS_err + e_OS_SS_err;
	
			ofstream datafile;
			datafile.open(("../../../../../TD-AnaNote/latex/tables/Tagging/"+(string)OutputDir + "tagPower"+ prefix[p] + ".tex").c_str(),std::ofstream::trunc);
			//datafile << "\\begin{table}[h]" << "\n";
			//datafile << "\\centering" << "\n";
			//datafile << "\\caption{The flavour tagging performances for only OS tagged, only SS tagged and both OS and SS tagged events";
			//if(A_is_in_B("Run1", prefix[p])) datafile << " for Run-I data"; 
			//else if(A_is_in_B("Run2", prefix[p])) datafile << " for Run-II data"; 
			//datafile << ".}\n";;
			datafile << "\\begin{tabular}{c c c c}" << "\n";
			datafile << "\\hline" << "\n";
			datafile << "\\hline" << "\n";
			if((string)channel == "norm")datafile << "$ B_s \\to D_s \\pi \\pi \\pi$";
			else if((string)channel == "signal")datafile << "$ B_s \\to D_s K \\pi \\pi$";
			datafile << " & $\\epsilon_{tag} [\\%]$ & $\\langle \\omega \\rangle [\\%] $ & $\\epsilon_{eff} [\\%]$ \\\\" << "\n";
			datafile << "\\hline" << "\n";
			datafile << std::fixed << std::setprecision(2) << "Only OS & " 
			<< eff_OS * 100. << " $\\pm$ " << sqrt(eff_OS_err) * 100. << " & " 
			<< w_OS * 100. << " $\\pm$ " << sqrt(w_OS_err) * 100. << " & "
			<< e_OS * 100.<< " $\\pm$ " << sqrt(e_OS_err) * 100. << "\\\\" << "\n";
	
			datafile << std::fixed << std::setprecision(2) << "Only SS & " 
			<< eff_SS * 100.<< " $\\pm$ " << sqrt(eff_SS_err) * 100. << " & " 
			<< w_SS * 100.<< " $\\pm$ " << sqrt(w_SS_err) * 100. << " & "
			<< e_SS * 100.<< " $\\pm$ " << sqrt(e_SS_err) * 100. << "\\\\" << "\n";
	
			datafile << std::fixed << std::setprecision(2) << "Both OS-SS & " 
			<< eff_OS_SS * 100.<< " $\\pm$ " << sqrt(eff_OS_SS_err) * 100. << " & " 
			<< w_OS_SS * 100.<< " $\\pm$ " << sqrt(w_OS_SS_err) * 100. << " & "
			<< e_OS_SS * 100.<< " $\\pm$ " << sqrt(e_OS_SS_err) * 100. << "\\\\" << "\n";
	
			datafile << "\\hline" << "\n" << std::fixed << std::setprecision(2) << "Combined & " 
			<< eff_tot * 100.<< " $\\pm$ " << sqrt(eff_tot_err) * 100. << " & " 
			<< w_tot * 100.<< " $\\pm$ " << sqrt(w_tot_err) * 100. << " & "
			<< e_tot * 100.<< " $\\pm$ " << sqrt(e_tot_err) * 100. << "\\\\" << "\n";
	
			datafile << "\\hline" << "\n";
			datafile << "\\hline" << "\n";
			datafile << "\\end{tabular}" << "\n";
			//datafile << "\\label{table:tagging" << prefix[p] << "}" << "\n";
			//datafile << "\\end{table}";	
		}
	}
    }

//     nBinst *= 2;
    nBinsAsym *= scale_asym;

    TH1D* h_t_fit = new TH1D("h_t_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_p_fit = new TH1D("h_t_p",";t [ps];Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_m_fit = new TH1D("h_t_m",";t [ps];Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_N_fit = new TH1D("h_t_N_fit",";t [ps];Events (a.u.) ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_t_Nbar_fit = new TH1D("h_t_Nbar_fit",";t [ps];Events (a.u.) ",nBinsAsym,0.,2.*pi/dm);
//     h_t_p_fit->SetBins(4,binsl);
//     h_t_m_fit->SetBins(4,binsl);

    TH1D* h_t_mixed_fit = new TH1D("h_t_mixed_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_fit = new TH1D("h_t_unmixed_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged_fit = new TH1D("h_t_untagegged_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_OS_fit = new TH1D("h_t_OS_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_OS_fit = new TH1D("h_t_mixed_OS_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_OS_fit = new TH1D("h_t_unmixed_OS_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_SS_fit = new TH1D("h_t_SS_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_SS_fit = new TH1D("h_t_mixed_SS_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_SS_fit = new TH1D("h_t_unmixed_SS_fit",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_fit_mp = new TH1D("h_t_fit_mp",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0p = new TH1D("h_t_fit_0p",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pp = new TH1D("h_t_fit_pp",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_mm = new TH1D("h_t_fit_mm",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0m = new TH1D("h_t_fit_0m",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pm = new TH1D("h_t_fit_pm",";t [ps];Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_dt_fit = new TH1D("h_dt_fit",";#sigma_{t} [ps];Events (a.u.) ",nBinst,0,0.15);
    TH1D* h_eta_OS_fit = new TH1D("h_eta_OS_fit",";#eta_{OS};Events (a.u.) ",nBinst,0,0.5);
    TH1D* h_eta_SS_fit = new TH1D("h_eta_SS_fit",";#eta_{SS};Events (a.u.) ",nBinst,0,0.5);
    TH1D* h_q_OS_fit = new TH1D("h_q_OS_fit",";q_{OS};Events (a.u.) ",3,-1.5,1.5);
    TH1D* h_q_SS_fit = new TH1D("h_q_SS_fit",";q_{SS};Events (a.u.) ",3,-1.5,1.5);
    TH1D* h_f_fit = new TH1D("h_f_fit",";q_{f};Events (a.u.) ",2,-2,2);

    TH1D* h_N_mixed_fit = new TH1D("h_N_mixed_fit",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_fit");

    TH1D* h_N_mixed_OS_fit = new TH1D("h_N_mixed_OS_fit",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_OS_fit = (TH1D*) h_N_mixed_OS_fit->Clone("h_N_unmixed_OS_fit");
    TH1D* h_N_mixed_SS_fit = new TH1D("h_N_mixed_SS_fit",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_SS_fit = (TH1D*) h_N_mixed_SS_fit->Clone("h_N_unmixed_SS_fit");

    TH1D* h_N_mixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_p_fit");
    TH1D* h_N_unmixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_m_fit");

    TH1D* h_N_mixed_p_fit_unfolded = new TH1D("h_N_mixed_p_fit_unfolded",";t/#tau;A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_m_fit");

    double N_OS_MC = 0;
    double N_SS_MC = 0;
    double N_OS_SS_MC = 0;
    double N_MC = 0;

    double w_OS_MC = 0;
    double w_SS_MC = 0;
    double w_OS_SS_MC = 0;
 	
    double D_OS_MC = 0;
    double D_SS_MC = 0;
    double D_OS_SS_MC = 0;
    double D_comb_MC = 0;

    double N_OS_all_MC = 0;
    double N_SS_all_MC = 0;
    double w_OS_all_MC = 0;
    double w_SS_all_MC = 0;
    double D_OS_all_MC = 0;
    double D_SS_all_MC = 0;
    
    double N_Run1_t0_MC = 0;
    double N_Run1_t1_MC = 0;
    double N_Run2_t0_MC = 0;
    double N_Run2_t1_MC = 0;

   if(fitGenMC){

	/// Fit with B2DX Pdf	
	//RooDataSet* data = new RooDataSet("data","data",RooArgSet(*r_t,*r_dt,*r_q,*r_mistag,*r_f));
	RooDataSet* protoData = new RooDataSet("protoData","protoData",RooArgSet(*r_dt,*r_q_OS,*r_q_SS,*r_f,*r_eta_OS,*r_eta_SS));
	//RooDataSet* protoData = new RooDataSet("protoData","protoData",RooArgSet(*r_dt,*r_eta_OS,*r_eta_SS));

	NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
	vector<double> myBinning = knot_positions.getVector();
	NamedParameter<double> knot_values("knot_values", 0.38,0.63,0.86,1.05,1.14,1.24,1.22);
	vector<double> values = knot_values.getVector() ;
	
	RooArgList tacc_list;
	for(int i= 0; i< values.size(); i++){
		tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i])));
	}
	tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString((int)values.size())).c_str(), ("coeff_"+anythingToString((int)values.size())).c_str(), 1.)));
	RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString((int)values.size()+1)).c_str(),("coeff_"+anythingToString((int)values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString((int)values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(r_t->getMax()) ));
	tacc_list.add(*coeff_last);	
	
	RooCubicSplineFun* spline = new RooCubicSplineFun("splinePdf", "splinePdf", *r_t, myBinning, tacc_list);        
	RooGaussEfficiencyModel* sample_efficiency = new RooGaussEfficiencyModel("sample_efficiency", "sample_efficiency", *r_t, *spline, RooRealConstant::value(0.),*r_dt, RooRealConstant::value(0.),RooRealConstant::value(1.));
			
	RooRealVar r_tau("tau", "decay time", 1./Gamma);    
	RooRealVar r_dgamma("dgamma", "dgamma", dGamma);
	RooRealVar r_dm("dm", "dm", dm);

        RooRealVar r_C("C", "C",C,0,1);
        RooFormulaVar r_Cbar("Cbar","-1. * @0",RooArgList(r_C));
        RooRealVar r_D("D", "D",D,-1,1);
        RooRealVar r_Dbar("Dbar", "Dbar",D_bar,-1,1);
        RooRealVar r_S("S", "S",S,-1,1);
        RooRealVar r_Sbar("Sbar", "Sbar",S_bar,-1,1);

	RooGaussEfficiencyModel* efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *r_t, *spline, RooRealConstant::value(0.), RooRealConstant::value(1.), RooRealConstant::value(0.), RooRealConstant::value(1.) );

	RooTruthModel* deltaFunc = new RooTruthModel("resmodelDelta", "resmodelDelta", *r_t);

        DecRateCoeff_Bd cosh_coeff_gen("cosh_coeff_gen",
                        "cosh_coeff_gen",
                        DecRateCoeff_Bd::kCosh ,
                        *r_f,
                        RooRealConstant::value(1.),
                        RooRealConstant::value(1.),
                        *r_q,
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(1.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.),
                        RooRealConstant::value(0.)); 
                        
        DecRateCoeff_Bd cos_coeff_gen("cos_coeff_gen",
                                   "cos_coeff_gen",
                                   DecRateCoeff_Bd::kCos ,
                                   *r_f,
                                   r_C,
                                   r_Cbar,
                                   *r_q,
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(1.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.),
                                  RooRealConstant::value(0.));  
                                            
        DecRateCoeff_Bd sinh_coeff_gen("sinh_coeff_gen",
                                   "sinh_coeff_gen",
                                   DecRateCoeff_Bd::kSinh ,
                                   *r_f,
                                   r_D,
                                   r_Dbar,
                                   *r_q,
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(1.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.));  
                                   
        DecRateCoeff_Bd sin_coeff_gen("sin_coeff_gen",
                                   "sin_coeff_gen",
                                   DecRateCoeff_Bd::kSin,
                                   *r_f,
                                   r_S,
                                   r_Sbar,
                                   *r_q,
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(1.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.),
                                   RooRealConstant::value(0.));       

        RooBDecay fitpdf_t("fitpdf_t", "decay time PDF for fitting",
                           *r_t,r_tau, r_dgamma, 
                           //RooRealConstant::value(1.),RooRealConstant::value(0.),RooRealConstant::value(0.),RooRealConstant::value(0.),
                           cosh_coeff_gen,sinh_coeff_gen,cos_coeff_gen,sin_coeff_gen,
                           r_dm, *efficiency, RooBDecay::SingleSided); 
        
        fitpdf_t.fitTo(*data,NumCPU(4));
	RooDataSet* toy = fitpdf_t.generate(RooArgSet(*r_t,*r_f,*r_q),1000);

	RooPlot* tframefit = r_t->frame();
        toy->plotOn(tframefit,Binning(nBinst));
        fitpdf_t.plotOn(tframefit);
        tframefit->Draw();
        c->Print("B2DX_fit.eps");      
        cout << " C = " << r_C.getVal() << " ; Pull = " << (r_C.getVal()-C)/r_C.getError() << endl;
        cout << " D = " << r_D.getVal() << " ; Pull = " << (r_D.getVal()-D)/r_D.getError() << endl;
        cout << " Dbar = " << r_Dbar.getVal() << " ; Pull = " << (r_Dbar.getVal()-D_bar)/r_Dbar.getError() << endl;
        cout << " S = " << r_S.getVal() << " ; Pull = " << (r_S.getVal()-S)/r_S.getError() << endl;
        cout << " Sbar = " << r_Sbar.getVal() << " ; Pull = " << (r_Sbar.getVal()-S_bar)/r_Sbar.getError() << endl;  

	DalitzEventPattern _pat(pat);
    	DalitzEvent evt_proto(_pat);
    	evt_proto.generateThisToPhaseSpace();

	for(int i = 0; i < 100000; i++){
		
		double t_MC = ranLux.Exp(1./Gamma);
		if(t_MC > max_TAU && t_MC < min_TAU)continue;
	
		double dt_MC = 0.04;
		
		double q_rand = ranLux.Uniform();
		int q_OS_MC = 0;
		if (q_rand < 1./2.  ) q_OS_MC = -1;
		if (q_rand > (1.-1./2.) ) q_OS_MC = 1;
		
		q_rand = ranLux.Uniform();
		int q_SS_MC = q_OS_MC;
		
		double eta_OS_MC = 0;
		double eta_SS_MC = 0;
	
		q_rand = ranLux.Uniform();
		int f_MC = 0;
		if (q_rand > .5) f_MC = -1;
		else f_MC = 1;
			
		DalitzEvent evt(evt_proto);
	
		evt.setWeight(1.);
		evt.setValueInVector(0, t_MC);
		evt.setValueInVector(1, dt_MC);
		evt.setValueInVector(2, f_MC);
		evt.setValueInVector(3, q_OS_MC);
		evt.setValueInVector(4, eta_OS_MC);
		evt.setValueInVector(5, q_SS_MC);
		evt.setValueInVector(6, eta_SS_MC);
		
		double pdfVal = t_pdf.getVal(evt);
	
		double weight = pdfVal;
		weight /=  exp(-t_MC*Gamma) / ( 1./Gamma * ( exp(-min_TAU*Gamma) - exp(-max_TAU*Gamma) ) ) ;
		//*  (abs(q_OS_MC)/2. * eff_tag_OS + ( 1. - abs(q_OS_MC)) * (1.-eff_tag_OS) ) ;
	
		h_t_fit->Fill(t_MC,weight);
		h_dt_fit->Fill(dt_MC,weight);
		if(evt.getValueFromVector(3) != 0)h_eta_OS_fit->Fill(evt.getValueFromVector(4),weight);
		if(evt.getValueFromVector(5) != 0)h_eta_SS_fit->Fill(evt.getValueFromVector(6),weight);
		
		int f_evt = evt.getValueFromVector(2);
		int q1 = evt.getValueFromVector(3);
		int q2 = evt.getValueFromVector(5);   
		int q_eff = 0;
		double w_eff = 0.5;
		
		std::pair<double, double> calibrated_mistag_os;
		std::pair<double, double> calibrated_mistag_ss;
		calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(evt);
		calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(evt);        
		
		double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
		double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
		
		if( p/(p+p_bar) > 0.5 ){ 
			q_eff = 1;
			w_eff = 1-p/(p+p_bar);
		}
		else if( p/(p+p_bar) < 0.5 ){
			q_eff = -1;
			w_eff = p/(p+p_bar);
		}
		
		if((string)channel=="signal"){
		
		if(q_eff==-1 && f_evt == 1){
			h_t_fit_mp->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
			h_N_mixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
			h_N_mixed_p_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),1./Gamma),weight);
			}
		}
		else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(evt.getValueFromVector(0),weight);
		else if(q_eff==1 && f_evt == 1){
			h_t_fit_pp->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
				h_N_unmixed_p_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
				h_N_unmixed_p_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),1./Gamma),weight);
			}
		}
		else if(q_eff==-1 && f_evt == -1){
			h_t_fit_mm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
			h_N_unmixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
			h_N_unmixed_m_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),1./Gamma),weight);
			}	
		}
		else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(evt.getValueFromVector(0),weight);
		else if(q_eff==1 && f_evt == -1){
			h_t_fit_pm->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max){
			h_N_mixed_m_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
			h_N_mixed_m_fit_unfolded->Fill(fmod(evt.getValueFromVector(0),1./Gamma),weight);
			}
		}
		}
		else {
		if(q_eff == 0)h_t_untagegged_fit->Fill(evt.getValueFromVector(0),weight);
		else if(q_eff*f_evt > 0  ){
			h_t_mixed_fit->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_mixed_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
		}
		else{ 
			h_t_unmixed_fit->Fill(evt.getValueFromVector(0),weight);
			if(w_eff<w_max)h_N_unmixed_fit->Fill(fmod(evt.getValueFromVector(0),2.*pi/dm),weight);
		}
	   	}
	}
    }
   else if (doPlots){ 
	for(int n = 0; n < N_plot_it; n++){   /// Multiple iterations needed to release memory 
		int N_sample = 500000;
		DalitzEventList sampleEvents;
// 		mps->getParPtr("Gamma")->setCurrentFitVal(1./15.);
		if(doSimFit) {
			if(N_Run1_t0>0)sampleEvents.Add(t_pdf_Run1_t0.generateToys(N_sample * N_Run1_t0/N,1,0));
			if(N_Run1_t1>0)sampleEvents.Add(t_pdf_Run1_t1.generateToys(N_sample *N_Run1_t1/N,1,1));
			if(N_Run2_t0>0)sampleEvents.Add(t_pdf_Run2_t0.generateToys(N_sample *N_Run2_t0/N,2,0));
			if(N_Run2_t1>0)sampleEvents.Add(t_pdf_Run2_t1.generateToys(N_sample *N_Run2_t1/N,2,1));
			if(N_Run2_17_t0>0)sampleEvents.Add(t_pdf_Run2_17_t0.generateToys(N_sample *N_Run2_17_t0/N,3,0));
			if(N_Run2_17_t1>0)sampleEvents.Add(t_pdf_Run2_17_t1.generateToys(N_sample *N_Run2_17_t1/N,3,1));
			if(N_Run2_18_t0>0)sampleEvents.Add(t_pdf_Run2_18_t0.generateToys(N_sample *N_Run2_18_t0/N,4,0));
			if(N_Run2_18_t1>0)sampleEvents.Add(t_pdf_Run2_18_t1.generateToys(N_sample *N_Run2_18_t1/N,4,1));
		}
		else sampleEvents.Add(t_pdf_Run1_t0.generateToys(N_sample));	

		for(int i = 0; i < sampleEvents.size(); i++){

			DalitzEvent evt = sampleEvents[i];
			double t_MC = evt.getValueFromVector(0) ;
			double dt_MC = evt.getValueFromVector(1) ;
			int f_MC = evt.getValueFromVector(2) ;
			int q_OS_MC = evt.getValueFromVector(3) ;
			double eta_OS_MC = evt.getValueFromVector(4) ;
			int q_SS_MC = evt.getValueFromVector(5);
			double eta_SS_MC = evt.getValueFromVector(6);
			int run_MC = evt.getValueFromVector(7) ;
			int trigger_MC = evt.getValueFromVector(8) ;

			double weight = 1; //t_pdf_Run1_t0.getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace();
			N_MC += weight;
		
			h_t_fit->Fill(t_MC,weight);
			h_dt_fit->Fill(dt_MC,weight);
			if(q_OS_MC != 0)h_eta_OS_fit->Fill(eta_OS_MC,weight);
			if(q_SS_MC != 0)h_eta_SS_fit->Fill(eta_SS_MC,weight);
			
			int f_evt = f_MC;
			int q1 = q_OS_MC;
			int q2 = q_SS_MC;   
			int q_eff = 0;
			double w_eff = 0.5;
		
			h_q_OS_fit->Fill(q1,weight);
			h_q_SS_fit->Fill(q2,weight);
			h_f_fit->Fill(f_evt,weight);
			
			if( q1 != 0)h_t_OS_fit->Fill(t_MC,weight);
			if( q2 != 0)h_t_SS_fit->Fill(t_MC,weight);

			std::pair<double, double> calibrated_mistag_os;
			std::pair<double, double> calibrated_mistag_ss;
			if(doSimFit){
				if(run_MC==1){
					calibrated_mistag_os = t_pdf_Run1_t0.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run1_t0.getCalibratedMistag_SS(eta_SS_MC);
				}
				else{
					calibrated_mistag_os = t_pdf_Run2_t0.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run2_t0.getCalibratedMistag_SS(eta_SS_MC);  
				}
			}
			else{
				calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eta_OS_MC);
				calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eta_SS_MC);        
			}
			
			double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
			double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
			
			if( p/(p+p_bar) > 0.5 ){ 
				q_eff = 1;
				w_eff = 1-p/(p+p_bar);
			}
			else if( p/(p+p_bar) < 0.5 ){
				q_eff = -1;
				w_eff = p/(p+p_bar);
			}
			
			if(q1 != 0 && q2 != 0){
				if(q_eff != 0){
					N_OS_SS_MC += weight;
					w_OS_SS_MC += w_eff * weight;
					D_OS_SS_MC += pow(1.-2.*w_eff,2)* weight;
				}
			}
			else if( q1 != 0){
				//q_eff = q1;  flip tag ???
				N_OS_MC +=weight;
				w_OS_MC += w_eff * weight; 
				D_OS_MC += pow(1.-2.*w_eff,2)* weight;
			}
			else if( q2 != 0){
				//q_eff = q2;
				N_SS_MC += weight;
				w_SS_MC += w_eff * weight; 
				D_SS_MC += pow(1.-2.*w_eff,2)* weight; 
			} 
		
			D_comb_MC += pow(1.-2.*w_eff,2)* weight;
			double D_res = 1.;
			if(doSimFit){
				if(run_MC==1){
					D_res = exp(-pow(t_pdf_Run1_t0.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
				else if(run_MC==2){
					D_res = exp(-pow(t_pdf_Run2_t0.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
				else if(run_MC==3){
					D_res = exp(-pow(t_pdf_Run2_17_t0.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
				else if(run_MC==4){
					D_res = exp(-pow(t_pdf_Run2_18_t0.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
			}
			double D_tag = 0.;
			if(q_eff != 0) D_tag = (1.-2.*(w_eff));
			double D_tot = D_tag * D_res;
// 			double D_tot =  D_res;
			if(!dilutionWeight)D_tot = 1.;
// 			double D_tot = (1.-2.*abs(w_eff)) * exp(-pow(t_pdf.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);

			if(q1 != 0) N_OS_all_MC += weight;
			if(q2 != 0)N_SS_all_MC += weight;
			
			if(q1>0){
					w_OS_all_MC +=  calibrated_mistag_os.first * weight;
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.first,2)* weight;
			} 
			else if(q1<0){
					w_OS_all_MC +=  calibrated_mistag_os.second * weight;	
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.second,2)* weight;
			}
		
			if(q2>0){
					w_SS_all_MC +=  calibrated_mistag_ss.first * weight;
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.first,2)* weight;
			} 
			else if(q2<0){
					w_SS_all_MC +=  calibrated_mistag_ss.second * weight;	
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.second,2)* weight;
			}
		
			if((string)channel=="signal"){

				if(f_evt == 1){ 
					h_t_p_fit->Fill(t_MC,weight);
				}
				else h_t_m_fit->Fill(t_MC,weight);
			
	   			if(q_eff == 1)h_t_N_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);
	    			else if(q_eff == -1)h_t_Nbar_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);

				if(q_eff==-1 && f_evt == 1){
					h_t_fit_mp->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_mixed_p_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);
						h_N_mixed_p_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight*D_tot);
					}
				}
				else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight*D_tot);
				else if(q_eff==1 && f_evt == 1){
					h_t_fit_pp->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_unmixed_p_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);
						h_N_unmixed_p_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight*D_tot);
					}
				}
				else if(q_eff==-1 && f_evt == -1){
					h_t_fit_mm->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_unmixed_m_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);
						h_N_unmixed_m_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight*D_tot);
					}	
				}
				else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(t_MC,weight*D_tot);
				else if(q_eff==1 && f_evt == -1){
					h_t_fit_pm->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_mixed_m_fit->Fill(fmod(t_MC-mint,2.*pi/dm),weight*D_tot);
						h_N_mixed_m_fit_unfolded->Fill(fmod(t_MC,1./Gamma),weight*D_tot);
					}
				}
			}
//  			else {   
				if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
				else if(q_eff*f_evt > 0  ){
					h_t_mixed_fit->Fill(t_MC,weight*D_tot);
					//if(w_eff<w_max)
					h_N_mixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);

					if( q1 != 0)h_t_mixed_OS_fit->Fill(t_MC,weight*D_tot);
					if( q1 != 0)h_N_mixed_OS_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);

					if( q2 != 0)h_t_mixed_SS_fit->Fill(t_MC,weight*D_tot);
					if( q2 != 0)h_N_mixed_SS_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);
				}
				else{ 
					h_t_unmixed_fit->Fill(t_MC,weight*D_tot);
					//if(w_eff<w_max)
					h_N_unmixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);

					if( q1 != 0)h_t_unmixed_OS_fit->Fill(t_MC,weight*D_tot);
					if( q1 != 0)h_N_unmixed_OS_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);

					if( q2 != 0)h_t_unmixed_SS_fit->Fill(t_MC,weight*D_tot);
					if( q2 != 0)h_N_unmixed_SS_fit->Fill(fmod(t_MC,2.*pi/dm),weight*D_tot);
				}
 			}
// 		}
	}
	
	cout << "N_MC = " << N_MC << endl << endl;
	cout << "N_OS = " << N_OS_all_MC << endl;
	cout << "N_SS = " << N_SS_all_MC << endl << endl;
	
	cout << "eff_OS = " <<(N_OS_all_MC)/N_MC << endl;
	cout << "eff_SS = " <<(N_SS_all_MC)/N_MC << endl;
		
	cout << "Tagging perfromance " << endl << endl;        
	cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
	
	cout << "OS  | " << (N_OS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_OS_all_MC)/(N_OS_MC+N_OS_SS_MC) << " | " << D_OS_all_MC/N_MC << endl;
	cout << "SS  | " << (N_SS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_SS_all_MC)/(N_SS_MC+N_OS_SS_MC) << " | " << D_SS_all_MC/N_MC << endl << endl;
	
	cout << "OS only  | " << N_OS_MC/N_MC << " | " <<  w_OS_MC/N_OS_MC << " | " << N_OS_MC/N_MC * D_OS_MC/N_OS_MC << endl;
	cout << "SS only  | " << N_SS_MC/N_MC << " | " <<  w_SS_MC/N_SS_MC << " | " << N_SS_MC/N_MC * D_SS_MC/N_SS_MC << endl;
	cout << "OS+SS    | " << N_OS_SS_MC/N_MC << " | " <<  w_OS_SS_MC/N_OS_SS_MC << " | " << N_OS_SS_MC/N_MC * D_OS_SS_MC/N_OS_SS_MC << endl;
	cout << "Combined | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC << " | "<<  (w_OS_MC+w_SS_MC+w_OS_SS_MC)/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << " | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC * D_comb_MC/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << endl << endl ;
	
	/// Plots
	TGraph* graph = new TGraph(2);
	graph->SetPoint(1,min_TAU,0);
	graph->SetPoint(2,max_TAU,0);
	graph->SetLineStyle(kDashed);

        TH1D* h_t_acc = new TH1D("h_t_acc",";t",nBinst,min_TAU,max_TAU);    
        for (int i = 1; i<=h_t_acc->GetNbinsX(); i++) {
		double val_acc = N_Run1_t0*t_pdf_Run1_t0.get_spline_Val(h_t_acc->GetXaxis()->GetBinCenter(i));
		val_acc += N_Run1_t1*t_pdf_Run1_t1.get_spline_Val(h_t_acc->GetXaxis()->GetBinCenter(i));
		val_acc += (N_Run2_t0+N_Run2_17_t0)*t_pdf_Run2_t0.get_spline_Val(h_t_acc->GetXaxis()->GetBinCenter(i));
		val_acc += (N_Run2_t1+N_Run2_17_t1)*t_pdf_Run2_t1.get_spline_Val(h_t_acc->GetXaxis()->GetBinCenter(i));
		h_t_acc->SetBinContent(i,val_acc/N);
        }

	TLegend leg_lhcb(0.75,0.8,0.9,0.9);
    	leg_lhcb.SetLineStyle(0);
    	leg_lhcb.SetLineColor(0);
	leg_lhcb.SetFillColor(0);
	leg_lhcb.SetTextFont(132);
	leg_lhcb.SetTextColor(1);
	leg_lhcb.SetTextSize(0.06);
	leg_lhcb.SetTextAlign(12);
	leg_lhcb.AddEntry((TObject*)0,"#font[22]{LHCb}","");
        stringstream ss ;

	TString title_t = "Yield / (";
	ss << std::fixed << std::setprecision(1) << (h_t->GetXaxis()->GetXmax() -h_t->GetXaxis()->GetXmin())/h_t->GetNbinsX();  
	title_t += ss.str();
	title_t += " ps)";
	h_t->GetYaxis()->SetTitle(title_t);
	h_t->Scale(N/h_t->Integral());

	h_t->SetMinimum(0.01);    
	h_t->SetLineColor(kBlack);
	h_t->Draw("e");
	h_t_fit->SetLineColor(kBlue);
	h_t_fit->SetLineWidth(3);
	h_t_fit->SetMarkerColor(kBlue); 
	h_t_fit->DrawNormalized("histcsame",N);

	h_t_acc->SetLineColor(kRed);
	h_t_acc->SetLineWidth(3);
	h_t_acc->SetLineStyle(kDashed);
	h_t_acc->DrawNormalized("histcsame",N);

	leg_lhcb.Draw();
	c->Print(((string)OutputDir+"h_t.eps").c_str());
	c->Print(((string)OutputDir+"h_t.pdf").c_str());
	c->Print(((string)OutputDir+"h_t.png").c_str());
	c->Print(((string)OutputDir+"h_t.C").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t.pdf").c_str());
	
	h_t_acc->DrawNormalized("histcsame",N/2);
	c->Print(((string)OutputDir+"h_acc.eps").c_str());


	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_t_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
	gPad->SetLogy(0);
	
	h_t_OS->SetMinimum(0.1);    
	h_t_OS->SetLineColor(kBlack);
	h_t_OS->DrawNormalized("e",1);
	h_t_OS_fit->SetLineColor(kBlue);
	h_t_OS_fit->SetLineWidth(3);
	h_t_OS_fit->SetMarkerColor(kBlue); 
	h_t_OS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_OS.eps").c_str());

	h_t_SS->SetMinimum(0.1);    
	h_t_SS->SetLineColor(kBlack);
	h_t_SS->DrawNormalized("e",1);
	h_t_SS_fit->SetLineColor(kBlue);
	h_t_SS_fit->SetLineWidth(3);
	h_t_SS_fit->SetMarkerColor(kBlue); 
	h_t_SS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_SS.eps").c_str());

	h_t_untagegged->SetMinimum(0.1);    
	h_t_untagegged->SetLineColor(kBlack);
	h_t_untagegged->DrawNormalized("e",1);
	h_t_untagegged_fit->SetLineColor(kBlue);
	h_t_untagegged_fit->SetLineWidth(3);
	h_t_untagegged_fit->SetMarkerColor(kBlue); 
	h_t_untagegged_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_untagged.eps").c_str());

	h_dt->SetMinimum(0);        
	h_dt->SetLineColor(kBlack);
	h_dt->DrawNormalized("e1",1);
	h_dt_fit->SetLineColor(kBlue);
	h_dt_fit->SetLineWidth(3);
	h_dt_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_dt.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_dt.pdf").c_str());
	
	h_eta_OS->SetMinimum(0);        
	h_eta_OS->SetLineColor(kBlack);
	h_eta_OS->DrawNormalized("e1",1);
	h_eta_OS_fit->SetLineColor(kBlue);
	h_eta_OS_fit->SetLineWidth(3);
	h_eta_OS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_eta_OS.pdf").c_str());
	
	h_eta_SS->SetMinimum(0);        
	h_eta_SS->SetLineColor(kBlack);
	h_eta_SS->DrawNormalized("e1",1);
	h_eta_SS_fit->SetLineColor(kBlue);
	h_eta_SS_fit->SetLineWidth(3);
	h_eta_SS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_eta_SS.pdf").c_str());
	
	h_q_OS->SetMinimum(0);        
	h_q_OS->SetLineColor(kBlack);
	h_q_OS->DrawNormalized("e1",1);
	h_q_OS_fit->SetLineColor(kBlue);
	h_q_OS_fit->SetLineWidth(3);
	h_q_OS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_q_OS.pdf").c_str());
	
	h_q_SS->SetMinimum(0);        
	h_q_SS->SetLineColor(kBlack);
	h_q_SS->DrawNormalized("e1",1);
	h_q_SS_fit->SetLineColor(kBlue);
	h_q_SS_fit->SetLineWidth(3);
	h_q_SS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_q_SS.pdf").c_str());
	
	h_f->SetMinimum(0);        
	h_f->SetLineColor(kBlack);
	h_f->DrawNormalized("e1",1);
	h_f_fit->SetLineColor(kBlue);
	h_f_fit->SetLineWidth(3);
	h_f_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_f.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_f.pdf").c_str());
	
	if((string)channel=="norm"){
	
		ss.str("");
		title_t = "Yield / (";
		ss << std::fixed << std::setprecision(1) << (h_t_mixed->GetXaxis()->GetXmax() -h_t_mixed->GetXaxis()->GetXmin())/h_t_mixed->GetNbinsX();  
		title_t += ss.str();
		title_t += " ps)";
		h_t_mixed->GetYaxis()->SetTitle(title_t);
		h_t_mixed->Scale((N_OS+N_SS+N_OS_SS)/2./h_t_mixed->Integral());

		h_t_mixed->SetMarkerColor(kRed); 
		h_t_mixed->SetLineColor(kRed);

		//h_t_mixed->SetMaximum(h_t_mixed->GetMaximum()*1.1);
		//h_t_mixed->SetMinimum(h_t_mixed->GetMinimum()*0.5);

		h_t_mixed->Draw("e");
		
		h_t_unmixed->SetMarkerColor(kBlue); 
		h_t_unmixed->SetLineColor(kBlue);
		h_t_unmixed->DrawNormalized("esame",(N_OS+N_SS+N_OS_SS)/2.);
		
		h_t_mixed_fit->SetMarkerColor(kRed); 
		h_t_mixed_fit->SetLineColor(kRed);
		h_t_mixed_fit->DrawNormalized("histcsame",(N_OS+N_SS+N_OS_SS)/2.);
		
		h_t_unmixed_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_fit->SetLineColor(kBlue);
		h_t_unmixed_fit->DrawNormalized("histcsame",(N_OS+N_SS+N_OS_SS)/2.);
		
		h_t_untagegged->SetMarkerColor(kGreen); 
		h_t_untagegged->SetLineColor(kGreen);
		//h_t_untagegged->DrawNormalized("esame",1);
		
		h_t_untagegged_fit->SetMarkerColor(kGreen); 
		h_t_untagegged_fit->SetLineColor(kGreen);
		//h_t_untagegged_fit->DrawNormalized("histcsame",1);

		TLegend leg(0.6,0.65,0.9,0.9,"");
		leg.SetLineStyle(0);
		leg.SetLineColor(0);
		leg.SetFillColor(0);
		leg.SetTextFont(132);
		leg.SetTextColor(1);
		leg.SetTextSize(0.06);
		leg.SetTextAlign(12);

		leg.AddEntry((TObject*)0,"#font[22]{LHCb}","h");
		TLegendEntry* le = leg.AddEntry(h_t_mixed,"B^{0}_{s}(t)#rightarrow D_{s}^{-}#pi^{+}#pi^{-}#pi^{+}","l");
		le->SetTextColor(kRed);
		le->SetTextSize(0.05);
		le = leg.AddEntry(h_t_unmixed,"#bar{B^{0}_{s}}(t)#rightarrow D_{s}^{-}#pi^{+}#pi^{-}#pi^{+}","l");
	        le->SetTextColor(kBlue);    
		le->SetTextSize(0.05);
		leg.Draw();
	
		c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
		c->Print(((string)OutputDir+"h_t_mixed.C").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed.pdf").c_str());
	
		TLegend leg2(0.475,0.75,0.625,0.9,"");
		leg2.SetLineStyle(0);
		leg2.SetLineColor(0);
		leg2.SetFillColor(0);
		leg2.SetTextFont(132);
		leg2.SetTextColor(1);
		leg2.SetTextSize(0.06);
		leg2.SetTextAlign(12);

		TLegendEntry* le_lhcb = leg2.AddEntry((TObject*)0,"#font[22]{LHCb}","h");
		le_lhcb->SetTextAlign(22);
		leg2.AddEntry((TObject*)0,"","");
		leg2.AddEntry((TObject*)0,"","");

		
		TH1D* h_asym = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);	
		double max_asym = 0.5; //max(h_asym->GetMaximum(),fabs(h_asym->GetMinimum())) *1.25;
		h_asym->SetMaximum(max_asym);
		h_asym->SetMinimum(-max_asym);

		TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
		h_asym_fit->SetLineColor(kBlue);
		h_asym_fit->SetLineWidth(3);

		h_asym->SetName("h_asym");
		h_asym_fit->SetName("h_asym_fit");

		h_asym->Draw("e");
		h_asym_fit->Draw("histcsame");
 		leg2.Draw();
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());
		c->Print(((string)OutputDir+"h_asym.pdf").c_str());
		c->Print(((string)OutputDir+"h_asym.png").c_str());
		c->Print(((string)OutputDir+"h_asym.C").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym.pdf").c_str());

		h_t_mixed_OS->SetMarkerColor(kRed); 
		h_t_mixed_OS->SetLineColor(kRed);
		h_t_mixed_OS->DrawNormalized("e",1);
		
		h_t_unmixed_OS->SetMarkerColor(kBlue); 
		h_t_unmixed_OS->SetLineColor(kBlue);
		h_t_unmixed_OS->DrawNormalized("esame",1);
		
		h_t_mixed_OS_fit->SetMarkerColor(kRed); 
		h_t_mixed_OS_fit->SetLineColor(kRed);
		h_t_mixed_OS_fit->DrawNormalized("histcsame",1);
		
		h_t_unmixed_OS_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_OS_fit->SetLineColor(kBlue);
		h_t_unmixed_OS_fit->DrawNormalized("histcsame",1);
		c->Print(((string)OutputDir+"h_t_mixed_OS.eps").c_str());

		h_t_mixed_SS->SetMarkerColor(kRed); 
		h_t_mixed_SS->SetLineColor(kRed);
		h_t_mixed_SS->DrawNormalized("e",1);
		
		h_t_unmixed_SS->SetMarkerColor(kBlue); 
		h_t_unmixed_SS->SetLineColor(kBlue);
		h_t_unmixed_SS->DrawNormalized("esame",1);
		
		h_t_mixed_SS_fit->SetMarkerColor(kRed); 
		h_t_mixed_SS_fit->SetLineColor(kRed);
		h_t_mixed_SS_fit->DrawNormalized("histcsame",1);
		
		h_t_unmixed_SS_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_SS_fit->SetLineColor(kBlue);
		h_t_unmixed_SS_fit->DrawNormalized("histcsame",1);
		c->Print(((string)OutputDir+"h_t_mixed_SS.eps").c_str());
	}
	
	else{
	
		h_t_mp->Scale(1./h_t_mp->Integral());
		double maxY= h_t_mp->GetMaximum()*1.3;        
		h_t_mp->SetMinimum(0.);  
		h_t_mp->SetMaximum(maxY);        
		h_t_mp->SetMarkerColor(kBlue); 
		h_t_mp->SetLineColor(kBlue);
		h_t_mp->DrawNormalized("e",1);
	
		h_t_pp->SetMarkerColor(kRed); 
		h_t_pp->SetLineColor(kRed);
		h_t_pp->DrawNormalized("esame",1);
		
		h_t_fit_mp->SetMarkerColor(kBlue); 
		h_t_fit_mp->SetLineColor(kBlue);
		h_t_fit_mp->DrawNormalized("histcsame",1);
		
		h_t_fit_pp->SetMarkerColor(kRed); 
		h_t_fit_pp->SetLineColor(kRed);
		h_t_fit_pp->DrawNormalized("histcsame",1);
		
		c->Print(((string)OutputDir+"h_t_mixed_p.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed_p.pdf").c_str());
	
		h_t_pm->Scale(1./h_t_pm->Integral());
		h_t_pm->SetMinimum(0.);     
		h_t_pm->SetMaximum(maxY);        
		h_t_pm->SetMarkerColor(kBlue); 
		h_t_pm->SetLineColor(kBlue);
		h_t_pm->DrawNormalized("e",1);
		
		h_t_mm->SetMarkerColor(kRed); 
		h_t_mm->SetLineColor(kRed);
		h_t_mm->DrawNormalized("esame",1);
		
		h_t_fit_pm->SetMarkerColor(kBlue); 
		h_t_fit_pm->SetLineColor(kBlue);
		h_t_fit_pm->DrawNormalized("histcsame",1);
		
		h_t_fit_mm->SetMarkerColor(kRed); 
		h_t_fit_mm->SetLineColor(kRed);
		h_t_fit_mm->DrawNormalized("histcsame",1);
		
		c->Print(((string)OutputDir+"h_t_mixed_m.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());
	
		cout << h_N_unmixed_p->GetEntries() << endl;
		cout << h_N_mixed_p->GetEntries() << endl;
	
		TH1D* h_asym_p = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_mixed_p);	
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);	
		h_asym_p_fit->SetLineColor(kRed);
		h_asym_p->Draw("e");
		h_asym_p_fit->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_p.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_p.pdf").c_str());

	
	
		TH1D* h_asym_m = (TH1D*) h_N_unmixed_m->GetAsymmetry(h_N_mixed_m);	
		//h_asym_m->SetMinimum(-20);
		//h_asym_m->SetMaximum(20);
		TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);	
		h_asym_m_fit->SetLineColor(kRed);
		h_asym_m->Draw("e");
		h_asym_m_fit->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_m.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_m.pdf").c_str());
	

		h_asym_m->SetName("h_asym_m");
		h_asym_p->SetName("h_asym_p");
		h_asym_p_fit->SetName("h_asym_p_fit");
		h_asym_m_fit->SetName("h_asym_m_fit");

		double max_asym = 0.5; // max(max(h_asym_p->GetMaximum(),h_asym_m->GetMaximum()),fabs(min(h_asym_p->GetMinimum(),h_asym_m->GetMinimum()))) *1.25;
		h_asym_p->SetMaximum(max_asym);
		h_asym_p->SetMinimum(-max_asym);	
// 		h_asym_p->SetMarkerColor(kRed);
// 		h_asym_p->SetLineColor(kRed);

		h_asym_m->SetMarkerSize(1.15);
		h_asym_p->SetMarkerSize(1.15);
		h_asym_m->SetMarkerStyle(21);

		h_asym_p->SetLineColor(kGray+3);
		h_asym_p->SetMarkerColor(kGray+3);
		h_asym_p->Draw("e");
		h_asym_m->SetLineColor(kGreen+3);
		h_asym_m->SetMarkerColor(kGreen+3);
		h_asym_m->Draw("esame");
		h_asym_p_fit->SetLineColor(kGray+3);
		h_asym_p_fit->SetLineWidth(5);
		h_asym_m_fit->SetLineWidth(5);
		h_asym_p_fit->Draw("histcsame");
		h_asym_m_fit->SetLineColor(kGreen+3);
		h_asym_m_fit->SetLineStyle(kDashed);
		h_asym_m_fit->Draw("histcsame");

		TLegend leg2(0.4,0.7,0.6,0.9,"");
		leg2.SetLineStyle(0);
		leg2.SetLineColor(0);
		leg2.SetFillColor(0);
		leg2.SetTextFont(132);
		leg2.SetTextColor(1);
		leg2.SetTextSize(0.06);
		leg2.SetTextAlign(12);

		TLegendEntry* le_lhcb = leg2.AddEntry((TObject*)0,"#font[22]{LHCb}","hc");
		le_lhcb->SetTextAlign(22);
		TLegendEntry* le_p = leg2.AddEntry(h_asym_p_fit,"f=D_{s}^{-}K^{+}#pi^{-}#pi^{+}","l");
		le_p->SetTextColor(kGray+3); 
		le_p->SetTextSize(0.055);                   
		TLegendEntry* le_m  = leg2.AddEntry(h_asym_m_fit,"#bar{f}=D_{s}^{+}K^{-}#pi^{+}#pi^{-}","l");
		le_m->SetTextColor(kGreen+3); 
		le_m->SetTextSize(0.055);                   
		leg2.Draw();                   
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());	
		c->Print(((string)OutputDir+"h_asym.C").c_str());	
		c->Print(((string)OutputDir+"h_asym.pdf").c_str());	
		c->Print(((string)OutputDir+"h_asym.png").c_str());	

		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym.pdf").c_str());
	
		TH1D* h_asym_p_fit_unfolded = (TH1D*) h_N_unmixed_p_fit_unfolded->GetAsymmetry(h_N_mixed_p_fit_unfolded);	
		h_asym_p_fit_unfolded->SetLineWidth(3);
		h_asym_p_fit_unfolded->SetLineColor(kRed);
		h_asym_p_fit_unfolded->Draw("histc");
	
		TH1D* h_asym_m_fit_unfolded = (TH1D*) h_N_unmixed_m_fit_unfolded->GetAsymmetry(h_N_mixed_m_fit_unfolded);
		h_asym_m_fit_unfolded->SetLineWidth(3);
		h_asym_m_fit_unfolded->SetLineColor(kBlue);
		h_asym_m_fit_unfolded->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_unfolded.eps").c_str());	
		c->Print(((string)OutputDir+"h_asym_unfolded.pdf").c_str());	
	
		TH1D* h_asym_p_unfolded = (TH1D*) h_N_unmixed_p_unfolded->GetAsymmetry(h_N_mixed_p_unfolded);	
		h_asym_p_unfolded->SetLineWidth(3);
		h_asym_p_unfolded->SetLineColor(kRed);
		h_asym_p_unfolded->Draw("e");
	
		TH1D* h_asym_m_unfolded = (TH1D*) h_N_unmixed_m_unfolded->GetAsymmetry(h_N_mixed_m_unfolded);
		h_asym_m_unfolded->SetLineWidth(3);
		h_asym_m_unfolded->SetLineColor(kBlue);
		h_asym_m_unfolded->Draw("esame");
	
		h_asym_p_fit_unfolded->Draw("histcsame");
		h_asym_m_fit_unfolded->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_unfolded_fit.eps").c_str());	


// 		TH1D* h_asym_CP_f = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);
// 		h_asym_CP_f->SetTitle(";t modulo (2#pi/#Deltam_{s}) [ps]; A_{CP}^{#LTf#GT} ");	
// 		h_asym_CP_f->GetYaxis()->SetTitleOffset(0.9);	
// 		//h_asym_p->SetMinimum(-20);
// 		//h_asym_p->SetMaximum(20);
// 		TH1D* h_asym_CP_f_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
// 		h_asym_CP_f_fit->SetLineColor(kBlue+1);
// 		h_asym_CP_f_fit->SetLineWidth(5);
// 		h_asym_CP_f->Draw("e");
// 		h_asym_CP_f_fit->Draw("histcsame");
// 		graph->Draw("same");
// 		c->Print(((string)OutputDir+"h_asym_CP_f.eps").c_str());
// 		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_f.pdf").c_str());

		TH1D* h_asym_CP_f = (TH1D*) h_t_N->GetAsymmetry(h_t_Nbar);
		h_asym_CP_f->SetTitle(";t modulo (2#pi/#Deltam_{s}) [ps]; A_{CP}^{#LTf#GT} ");	
		h_asym_CP_f->GetYaxis()->SetTitleOffset(0.9);	
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_CP_f_fit = (TH1D*) h_t_N_fit->GetAsymmetry(h_t_Nbar_fit);	
		h_asym_CP_f_fit->SetLineColor(kBlue+1);
		h_asym_CP_f_fit->SetLineWidth(5);
		h_asym_CP_f->Draw("e");
		h_asym_CP_f_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym_CP_f.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_f.pdf").c_str());


		TH1D* h_asym_CP_i = (TH1D*) h_t_p->GetAsymmetry(h_t_m);	
		h_asym_CP_i->SetTitle(";t [ps]; A_{CP}^{#LTi#GT} ");	
		h_asym_CP_i->GetYaxis()->SetTitleOffset(0.9);		
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_CP_i_fit = (TH1D*) h_t_p_fit->GetAsymmetry(h_t_m_fit);	
		h_asym_CP_i_fit->SetLineColor(kBlue+1);
		h_asym_CP_i_fit->SetLineWidth(5);
		h_asym_CP_i->Draw("e");
		h_asym_CP_i_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym_CP_i.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_i.pdf").c_str());


/*		TH1D* h_asymCP = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_unmixed_m);	
		TH1D* h_asymCP_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);	
		h_asymCP_fit->SetLineColor(kBlack);
		h_asymCP_fit->SetLineWidth(5);
// 		h_asymCP_fit->SetMinimum(-1.);
// 		h_asymCP_fit->SetMaximum(1.);
		h_asymCP->Draw("e1");
		h_asymCP_fit->Draw("histc");

		TH1D* h_asymCP2 = (TH1D*) h_N_mixed_m->GetAsymmetry(h_N_mixed_p);	
		TH1D* h_asymCP2_fit = (TH1D*) h_N_mixed_m_fit->GetAsymmetry(h_N_mixed_p_fit);	
		h_asymCP2_fit->SetLineColor(kRed);
		h_asymCP2_fit->SetLineWidth(5);
// 		h_asymCP_fit->SetMinimum(-1.);
// 		h_asymCP_fit->SetMaximum(1.);
		h_asymCP2->Draw("e1");
		h_asymCP2_fit->Draw("histcsame");

		c->Print(((string)OutputDir+"h_asymCP.eps").c_str());	

		h_asym_p->Add(h_asym_m,-1.);	
		h_asym_p_fit->Add(h_asym_m_fit,-1.);	
		h_asym_p_fit->SetLineColor(kRed);
		h_asym_p_fit->SetLineWidth(5);
// 		h_asym_asym_fit->SetMinimum(-1.);
// 		h_asym_asym_fit->SetMaximum(1.);
		h_asym_p_fit->SetYTitle("A_{CP}");
		h_asym_p->Draw("e1");
		h_asym_p_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_delta_asym.eps").c_str());	


// 		h_t_pp->Add(h_t_mp);
// 		h_t_mm->Add(h_t_pm);
// 		h_t_pp->Rebin(10);
// 		h_t_mm->Rebin(10);
		TH1D* h_asym2 = (TH1D*) h_t_m->GetAsymmetry(h_t_p);	

		h_t_fit_pp->Add(h_t_fit_mp);
		h_t_fit_mm->Add(h_t_fit_pm);
		h_t_fit_pp->Rebin(10);
		h_t_fit_mm->Rebin(10);
		TH1D* h_asym2_fit = (TH1D*) h_t_fit_mm->GetAsymmetry(h_t_fit_pp);	
		h_asym2_fit->SetLineColor(kBlack);
		h_asym2_fit->SetLineWidth(5);
// 		h_asym2_fit->SetMinimum(-1.);
// 		h_asym2_fit->SetMaximum(1.);
		h_asym2->Draw("e1");
		h_asym2_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym2.eps").c_str());


		h_N_unmixed_p->Add(h_N_mixed_m,1.);	
		h_N_unmixed_m->Add(h_N_mixed_p,1.);	
		TH1D* h_asym3 = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_unmixed_m);	

		h_N_unmixed_p_fit->Add(h_N_mixed_m_fit,1.);	
		h_N_unmixed_m_fit->Add(h_N_mixed_p_fit,1.);	
		TH1D* h_asym3_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);	
		h_asym3->Draw("e1");
		h_asym3_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym3.eps").c_str());*/
	
	}

    }
	
    if(do2DScan == 1){
        cout << "Now doing 2D scan:" << endl;
        
        Neg2LL fcn(t_pdf, eventList_f);    
        Neg2LL fcn_bar(t_pdf, eventList_f_bar);    
        
        int scanBins=20;
        double scanMin=0, scanMax=360;
        double nSigmaZoom = 2;
        double scanMinGammaZoom=min(gamma.meanInit(), gamma.mean()) - nSigmaZoom*gamma.err();
        double scanMaxGammaZoom=max(gamma.meanInit(), gamma.mean()) + nSigmaZoom*gamma.err();
        double scanMinDeltaZoom=min(delta.meanInit(), delta.mean()) - nSigmaZoom*delta.err();
        double scanMaxDeltaZoom=max(delta.meanInit(), delta.mean()) + nSigmaZoom*delta.err();
        double gammaZoomRange = scanMaxGammaZoom - scanMinGammaZoom;
        double deltaZoomRange = scanMaxDeltaZoom - scanMinDeltaZoom;
        
        TFile* scanFile = new TFile("scan.root", "RECREATE");
        TH2D* scanHisto = new TH2D("scan", "; #gamma [deg]; #delta [deg]", scanBins,  scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanHistoP = new TH2D("scanP", "; #gamma [deg]; #delta [deg]", scanBins,  scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanHistoM = new TH2D("scanM" , "; #gamma [deg]; #delta [deg]", scanBins, scanMin, scanMax, scanBins, scanMin, scanMax);
        TH2D* scanZoomHisto = new TH2D("scanZoom", "; #gamma [deg]; #delta [deg]", scanBins, scanMinGammaZoom, scanMaxGammaZoom, scanBins, scanMinDeltaZoom, scanMaxDeltaZoom);
        
        double scanMinLL=-9999;
        double scanMinLLP=-9999;
        double scanMinLLM=-9999;
        double scanMinLLZ=-9999;
        
        for(int i=0; i < scanBins; i++){
            double gamma_value = ((double)i+0.5)*360/((double)scanBins);
            gamma.setCurrentFitVal(gamma_value);
            for(int j=0; j < scanBins; j++){
                double delta_value = ((double)j+0.5)*360/((double)scanBins);
                delta.setCurrentFitVal(delta_value);
                
                double v = neg2LL.getNewVal();  
                if( (i==0 && j==0) || v < scanMinLL) scanMinLL=v;
                scanHisto->Fill(gamma_value, delta_value, v);
                
                double vP = fcn.getNewVal();  
                if( (i==0 && j==0) || vP < scanMinLLP) scanMinLLP=vP;
                scanHistoP->Fill(gamma_value, delta_value, vP);
                
                double vM = fcn_bar.getNewVal();  
                if( (i==0 && j==0) || vM < scanMinLLM) scanMinLLM=vM;
                scanHistoM->Fill(gamma_value, delta_value, vM);
            }
        }
        for(int i=0; i < scanBins; i++){
            double gamma_value = scanMinGammaZoom + ((double)i+0.5) * gammaZoomRange/((double)scanBins);
            gamma.setCurrentFitVal(gamma_value);
            for(int j=0; j < scanBins; j++){
                double delta_value = scanMinDeltaZoom + ((double)j+0.5) * deltaZoomRange/((double)scanBins);
                delta.setCurrentFitVal(delta_value);
                double v = neg2LL.getNewVal();
                
                if( (i==0 && j==0) || v < scanMinLLZ) scanMinLLZ=v;
                
                scanZoomHisto->Fill(gamma_value, delta_value, v);
            }
        }
        
        for(int i=0; i < scanBins; i++){
            double gamma_value = ((double)i+0.5)*360/((double)scanBins);
            for(int j=0; j < scanBins; j++){
                double delta_value = ((double)j+0.5)*360/((double)scanBins);
                scanHisto->Fill(gamma_value, delta_value, -scanMinLL);
                scanHistoP->Fill(gamma_value, delta_value, -scanMinLLP);
                scanHistoM->Fill(gamma_value, delta_value, -scanMinLLM);
            }
        }
        for(int i=0; i < scanBins; i++){
            double gamma_value = scanMinGammaZoom + ((double)i+0.5) * gammaZoomRange/((double)scanBins);
            for(int j=0; j < scanBins; j++){
                double delta_value = scanMinDeltaZoom + ((double)j+0.5) * deltaZoomRange/((double)scanBins);
                scanZoomHisto->Fill(gamma_value, delta_value, -scanMinLLZ);
            }
        }
        scanFile->cd();
        scanHisto->Write();
        scanHistoP->Write();
        scanHistoM->Write();
        scanZoomHisto->Write();
        scanFile->Close();
        
        cout<< "done 2-D scan" << endl;
    }
    
//     if(plotAcceptance){
//     	Neg2LLMultiConstraint gauss_acc_run1_t0(MinuitParameterSet::getDefaultSet(),"_Run1");
// 	TH1D* t_acc = t_pdf_Run1_t0.plotSpline();
// 	t_acc->Draw("histc");
// 
// 	for(int i = 0; i < 100; i++){
// 		gauss_acc_run1_t0.smearInputValues();
// 		TH1D* t_acc_i = t_pdf_Run1_t0.plotSpline();
// 		t_acc->SetLineColor(kBlue);
// 		t_acc_i->Draw("histcsame");
// 	}
// 	t_acc->SetLineColor(kRed);
// 	t_acc->Draw("histcsame");
// 	c->Print("smearedAcc.eps");
//     }

    return;
}

void produceMarginalPdfs(){
    
    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);
    TString prefix = "";
    //TString prefix = "BsTaggingTool_";
    NamedParameter<double> min_year("min_year", 11);
    NamedParameter<double> max_year("max_year", 18);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);    

    /// Load files
    // Data
    Int_t q_OS,f,Ds_ID,q_SS;
    Double_t w_OS,w_SS;
    double sw;
    int run,year,Ds_finalState,trigger;
    double t,dt;
    
    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm_18.root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*DEC",1);
    tree_norm->SetBranchStatus("*PROB",1);
    tree_norm->SetBranchStatus("*OS",1);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("run",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("Ds_ID",1);

    tree_norm->SetBranchAddress("OS_Combination_DEC",&q_OS);
    tree_norm->SetBranchAddress("OS_Combination_PROB",&w_OS);
    tree_norm->SetBranchAddress("SS_Kaon_DEC",&q_SS);
    tree_norm->SetBranchAddress("SS_Kaon_PROB",&w_SS);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("Ds_ID",&Ds_ID);

    ///Make histograms
    int bins = 60;
    TH1D* h_w_OS_norm = new TH1D("h_w_OS_norm_comb","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1 = new TH1D("h_w_OS_norm_Run1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2 = new TH1D("h_w_OS_norm_Run2","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t0 = new TH1D("h_w_OS_norm_Run1_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t0 = new TH1D("h_w_OS_norm_Run2_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t1 = new TH1D("h_w_OS_norm_Run1_t1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t1 = new TH1D("h_w_OS_norm_Run2_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm = new TH1D("h_w_SS_norm_comb","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1 = new TH1D("h_w_SS_norm_Run1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2 = new TH1D("h_w_SS_norm_Run2","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t0 = new TH1D("h_w_SS_norm_Run1_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t0 = new TH1D("h_w_SS_norm_Run2_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t1 = new TH1D("h_w_SS_norm_Run1_t1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t1 = new TH1D("h_w_SS_norm_Run2_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm = new TH1D("h_q_OS_norm_comb","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1 = new TH1D("h_q_OS_norm_Run1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2 = new TH1D("h_q_OS_norm_Run2","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t0 = new TH1D("h_q_OS_norm_Run1_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t0 = new TH1D("h_q_OS_norm_Run2_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t1 = new TH1D("h_q_OS_norm_Run1_t1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t1 = new TH1D("h_q_OS_norm_Run2_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm = new TH1D("h_q_SS_norm_comb","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1 = new TH1D("h_q_SS_norm_Run1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2 = new TH1D("h_q_SS_norm_Run2","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t0 = new TH1D("h_q_SS_norm_Run1_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t0 = new TH1D("h_q_SS_norm_Run2_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t1 = new TH1D("h_q_SS_norm_Run1_t1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t1 = new TH1D("h_q_SS_norm_Run2_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm = new TH1D("h_q_f_norm_comb","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1 = new TH1D("h_q_f_norm_Run1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2 = new TH1D("h_q_f_norm_Run2","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t0 = new TH1D("h_q_f_norm_Run1_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t0 = new TH1D("h_q_f_norm_Run2_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t1 = new TH1D("h_q_f_norm_Run1_t1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t1 = new TH1D("h_q_f_norm_Run2_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm = new TH1D("h_t_norm_comb",";t [ps];Events (a.u.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run1 = new TH1D("h_t_norm_Run1",";t [ps];Events (a.u.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run2 = new TH1D("h_t_norm_Run2",";t [ps];Events (a.u.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm = new TH1D("h_dt_norm_comb",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1 = new TH1D("h_dt_norm_Run1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2 = new TH1D("h_dt_norm_Run2",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t0 = new TH1D("h_dt_norm_Run1_t0",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t0 = new TH1D("h_dt_norm_Run2_t0",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t1 = new TH1D("h_dt_norm_Run1_t1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t1 = new TH1D("h_dt_norm_Run2_t1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);

    TH1D* h_w_OS_norm_Run2_17 = new TH1D("h_w_OS_norm_Run2_17","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_17_t0 = new TH1D("h_w_OS_norm_Run2_17_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_17_t1 = new TH1D("h_w_OS_norm_Run2_17_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm_Run2_17 = new TH1D("h_w_SS_norm_Run2_17","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_17_t0 = new TH1D("h_w_SS_norm_Run2_17_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_17_t1 = new TH1D("h_w_SS_norm_Run2_17_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm_Run2_17 = new TH1D("h_q_OS_norm_Run2_17","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_17_t0 = new TH1D("h_q_OS_norm_Run2_17_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_17_t1 = new TH1D("h_q_OS_norm_Run2_17_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm_Run2_17 = new TH1D("h_q_SS_norm_Run2_17","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_17_t0 = new TH1D("h_q_SS_norm_Run2_17_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_17_t1 = new TH1D("h_q_SS_norm_Run2_17_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm_Run2_17 = new TH1D("h_q_f_norm_Run2_17","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_17_t0 = new TH1D("h_q_f_norm_Run2_17_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_17_t1 = new TH1D("h_q_f_norm_Run2_17_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm_Run2_17 = new TH1D("h_t_norm_Run2_17",";t [ps];Events (a.u.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm_Run2_17 = new TH1D("h_dt_norm_Run2_17",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17_t0 = new TH1D("h_dt_norm_Run2_17_t0",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17_t1 = new TH1D("h_dt_norm_Run2_17_t1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    

  TH1D* h_w_OS_norm_Run2_18 = new TH1D("h_w_OS_norm_Run2_18","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_18_t0 = new TH1D("h_w_OS_norm_Run2_18_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_18_t1 = new TH1D("h_w_OS_norm_Run2_18_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm_Run2_18 = new TH1D("h_w_SS_norm_Run2_18","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_18_t0 = new TH1D("h_w_SS_norm_Run2_18_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_18_t1 = new TH1D("h_w_SS_norm_Run2_18_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm_Run2_18 = new TH1D("h_q_OS_norm_Run2_18","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_18_t0 = new TH1D("h_q_OS_norm_Run2_18_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_18_t1 = new TH1D("h_q_OS_norm_Run2_18_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm_Run2_18 = new TH1D("h_q_SS_norm_Run2_18","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_18_t0 = new TH1D("h_q_SS_norm_Run2_18_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_18_t1 = new TH1D("h_q_SS_norm_Run2_18_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm_Run2_18 = new TH1D("h_q_f_norm_Run2_18","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_18_t0 = new TH1D("h_q_f_norm_Run2_18_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_18_t1 = new TH1D("h_q_f_norm_Run2_18_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm_Run2_18 = new TH1D("h_t_norm_Run2_18",";t [ps];Events (a.u.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm_Run2_18 = new TH1D("h_dt_norm_Run2_18",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_18_t0 = new TH1D("h_dt_norm_Run2_18_t0",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_18_t1 = new TH1D("h_dt_norm_Run2_18_t1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);


    ///loop over data events
    for(int i=0; i< tree_norm->GetEntries(); i++)
    {    
        tree_norm->GetEntry(i);
        if(year < min_year || year > max_year) continue;
	if(Ds_ID>0)f = -1;
	else f = 1;

        h_t_norm->Fill(t,sw);
        h_dt_norm->Fill(dt,sw);
        h_q_OS_norm->Fill((double)q_OS,sw);
        h_q_SS_norm->Fill((double)q_SS,sw);
        h_q_f_norm->Fill((double)f,sw);
        if(q_OS != 0)h_w_OS_norm->Fill(w_OS,sw);
        if(q_SS != 0)h_w_SS_norm->Fill(w_SS,sw);
            
        if(run==1){
            h_t_norm_Run1->Fill(t,sw);
            h_dt_norm_Run1->Fill(dt,sw);
            h_q_OS_norm_Run1->Fill((double)q_OS,sw);
            h_q_SS_norm_Run1->Fill((double)q_SS,sw);
	    h_q_f_norm_Run1->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run1->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run1->Fill(w_SS,sw);
	    if(trigger == 0){
		h_dt_norm_Run1_t0->Fill(dt,sw);
            	h_q_OS_norm_Run1_t0->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t0->Fill((double)q_SS,sw);
        	h_q_f_norm_Run1_t0->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t0->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t0->Fill(w_SS,sw);
	    }
	    else if(trigger == 1){
		h_dt_norm_Run1_t1->Fill(dt,sw);
            	h_q_OS_norm_Run1_t1->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t1->Fill((double)q_SS,sw);
       		h_q_f_norm_Run1_t1->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t1->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t1->Fill(w_SS,sw);
	    }
        }
        else if(run==2){
            h_t_norm_Run2->Fill(t,sw);
            h_dt_norm_Run2->Fill(dt,sw);
            h_q_OS_norm_Run2->Fill((double)q_OS,sw);
            h_q_SS_norm_Run2->Fill((double)q_SS,sw);
       	    h_q_f_norm_Run2->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run2->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run2->Fill(w_SS,sw);
	    if(trigger == 0){
		if(year < 17){
			h_dt_norm_Run2_t0->Fill(dt,sw);
			h_q_OS_norm_Run2_t0->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_t0->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_t0->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_t0->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_t0->Fill(w_SS,sw);
		}
		else if(year == 17) {
			h_dt_norm_Run2_17_t0->Fill(dt,sw);
			h_q_OS_norm_Run2_17_t0->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_17_t0->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_17_t0->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_17_t0->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_17_t0->Fill(w_SS,sw);
		}
		else if(year == 18) {
			h_dt_norm_Run2_18_t0->Fill(dt,sw);
			h_q_OS_norm_Run2_18_t0->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_18_t0->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_18_t0->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_18_t0->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_18_t0->Fill(w_SS,sw);
		}
	    }
	    else if(trigger == 1){
		if(year < 17){
			h_dt_norm_Run2_t1->Fill(dt,sw);
			h_q_OS_norm_Run2_t1->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_t1->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_t1->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_t1->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_t1->Fill(w_SS,sw);
		}
		else if(year == 17) {
			h_dt_norm_Run2_17_t1->Fill(dt,sw);
			h_q_OS_norm_Run2_17_t1->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_17_t1->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_17_t1->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_17_t1->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_17_t1->Fill(w_SS,sw);
		}
		else if(year == 18) {
			h_dt_norm_Run2_18_t1->Fill(dt,sw);
			h_q_OS_norm_Run2_18_t1->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_18_t1->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_18_t1->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_18_t1->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_18_t1->Fill(w_SS,sw);
		}
	    }
        }
    }
    
    TFile* out = new TFile("Mistag_pdfs.root","RECREATE");
    h_t_norm->Write();
    h_dt_norm->Write();
    h_q_OS_norm->Write();
    h_w_OS_norm->Write();
    h_q_SS_norm->Write();
    h_w_SS_norm->Write();
    h_q_f_norm->Write();

    h_t_norm_Run1->Write();
    h_dt_norm_Run1->Write();
    h_q_OS_norm_Run1->Write();
    h_w_OS_norm_Run1->Write();
    h_q_SS_norm_Run1->Write();
    h_w_SS_norm_Run1->Write();
    h_q_f_norm_Run1->Write();

    h_dt_norm_Run1_t0->Write();
    h_q_OS_norm_Run1_t0->Write();
    h_w_OS_norm_Run1_t0->Write();
    h_q_SS_norm_Run1_t0->Write();
    h_w_SS_norm_Run1_t0->Write();
    h_q_f_norm_Run1_t0->Write();

    h_dt_norm_Run1_t1->Write();
    h_q_OS_norm_Run1_t1->Write();
    h_w_OS_norm_Run1_t1->Write();
    h_q_SS_norm_Run1_t1->Write();
    h_w_SS_norm_Run1_t1->Write();
    h_q_f_norm_Run1_t1->Write();

    h_t_norm_Run2->Write();
    h_dt_norm_Run2->Write();
    h_q_OS_norm_Run2->Write();
    h_w_OS_norm_Run2->Write();
    h_q_SS_norm_Run2->Write();
    h_w_SS_norm_Run2->Write();
    h_q_f_norm_Run2->Write();

    h_dt_norm_Run2_t0->Write();
    h_q_OS_norm_Run2_t0->Write();
    h_w_OS_norm_Run2_t0->Write();
    h_q_SS_norm_Run2_t0->Write();
    h_w_SS_norm_Run2_t0->Write();
    h_q_f_norm_Run2_t0->Write();

    h_dt_norm_Run2_t1->Write();
    h_q_OS_norm_Run2_t1->Write();
    h_w_OS_norm_Run2_t1->Write();
    h_q_SS_norm_Run2_t1->Write();
    h_w_SS_norm_Run2_t1->Write();
    h_q_f_norm_Run2_t1->Write();

    h_dt_norm_Run2_17_t0->Write();
    h_q_OS_norm_Run2_17_t0->Write();
    h_w_OS_norm_Run2_17_t0->Write();
    h_q_SS_norm_Run2_17_t0->Write();
    h_w_SS_norm_Run2_17_t0->Write();
    h_q_f_norm_Run2_17_t0->Write();

    h_dt_norm_Run2_17_t1->Write();
    h_q_OS_norm_Run2_17_t1->Write();
    h_w_OS_norm_Run2_17_t1->Write();
    h_q_SS_norm_Run2_17_t1->Write();
    h_w_SS_norm_Run2_17_t1->Write();
    h_q_f_norm_Run2_17_t1->Write();

    h_dt_norm_Run2_18_t0->Write();
    h_q_OS_norm_Run2_18_t0->Write();
    h_w_OS_norm_Run2_18_t0->Write();
    h_q_SS_norm_Run2_18_t0->Write();
    h_w_SS_norm_Run2_18_t0->Write();
    h_q_f_norm_Run2_18_t0->Write();

    h_dt_norm_Run2_18_t1->Write();
    h_q_OS_norm_Run2_18_t1->Write();
    h_w_OS_norm_Run2_18_t1->Write();
    h_q_SS_norm_Run2_18_t1->Write();
    h_w_SS_norm_Run2_18_t1->Write();
    h_q_f_norm_Run2_18_t1->Write();

    out->Write();
}


void test_multiGaussConstraints(){

    //time reso part

    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.2,0.1);

    Neg2LLMultiConstraint gauss_constrains(MinuitParameterSet::getDefaultSet(),"_Run2");

    RooDataSet* data_cov = gauss_constrains.generateToys(100);

    TCanvas* c = new TCanvas();
    TF1 *nominalFunc = new TF1("nominalFunc", "[0]+[1]*x ", 0., 0.15);
    nominalFunc->SetParameters(0, offset_sigma_dt_Run2);
    nominalFunc->SetParameters(1, scale_sigma_dt_Run2);
    nominalFunc->SetLineWidth(2.5);
    nominalFunc->Draw();

    for(int i = 0 ; i < 100; i++){
					RooArgSet* xvec_cov= (RooArgSet*)data_cov->get(i);

					double p0 = ((RooRealVar*)xvec_cov->find("offset_sigma_dt_Run2"))->getVal(); 
					double p1 = ((RooRealVar*)xvec_cov->find("scale_sigma_dt_Run2"))->getVal(); 

					cout << "scaling function : " << p0 << " +/- " << p1 << " * dt "  << endl ;

					//plot it 
					TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 0.15);
					fitFunc->SetParameters(0, p0);
					fitFunc->SetParameters(1, p1);
					fitFunc->SetLineColor(i);
					fitFunc->Draw("LSAME");

     }

    nominalFunc->Draw("LSAME");
    c->Print("ScalingFunctions.eps");
    c->Close();


    TCanvas* c_new = new TCanvas();

    //time acceptance part

    FitParameter  c0_Run1_t0("c0_Run1_t0",1,0.,2.);
    FitParameter  c1_Run1_t0("c1_Run1_t0",1,0.,2.);
    FitParameter  c2_Run1_t0("c2_Run1_t0",1,0.,2.);
    FitParameter  c3_Run1_t0("c3_Run1_t0",1,0.,2.);

    Neg2LLMultiConstraint tagging_constrains_Run1_t0(MinuitParameterSet::getDefaultSet(),"_Tagging_Run1_t0");

    RooDataSet* tagging_cov_Run1_t0 = tagging_constrains_Run1_t0.generateToys(100);

    Double_t xAxis[4]  = {0.8, 1.6, 2.5 , 6.5};
    Double_t yNominal[4] = {5.7696e-01, 7.5715e-01, 8.8174e-01, 1.0844e+00};

   TGraphErrors *NominalSpline = new TGraphErrors(4,xAxis,yNominal);
   NominalSpline->SetTitle("Spline Coefficients c_{i}");
   NominalSpline->GetXaxis()->SetTitle("t [ps]");
   NominalSpline->GetYaxis()->SetTitle("c_{i}");
   NominalSpline->Draw();

    for(int i = 0 ; i < 100; i++){
					RooArgSet* xvec_cov_Run1_t0= (RooArgSet*)tagging_cov_Run1_t0->get(i);

                                        double c0 = ((RooRealVar*)xvec_cov_Run1_t0->find("c0_Run1_t0"))->getVal();
					double c1 = ((RooRealVar*)xvec_cov_Run1_t0->find("c1_Run1_t0"))->getVal();
                                        double c2 = ((RooRealVar*)xvec_cov_Run1_t0->find("c2_Run1_t0"))->getVal();
					double c3 = ((RooRealVar*)xvec_cov_Run1_t0->find("c3_Run1_t0"))->getVal();

                                        cout << "spline : c0 = " << c0 << " , c1 = " << c1 << " , c2 = " << c2 << " , c3= " << c3  << endl ;

   					Double_t yAxis[4]  = {c0, c1, c2, c3};
   					TGraphErrors *gr = new TGraphErrors(4,xAxis,yAxis);
   					gr->SetMarkerColor(i);
					gr->SetLineColor(i);
   					gr->Draw("SAME");

   }

   TGraphErrors *NominalSpline2 = new TGraphErrors(4,xAxis,yNominal);
   NominalSpline2->SetLineWidth(2.5);
   NominalSpline2->GetXaxis()->SetTitle("t [ps]");
   NominalSpline2->GetYaxis()->SetTitle("c_{i}");
   NominalSpline2->Draw("SAME");

   c_new->Print("SplineCoeffs.eps");



   TMatrixDSym* cov = tagging_constrains_Run1_t0.getCovMatrix();

   tagging_constrains_Run1_t0.smearInputValuesChol(0,0);
   tagging_constrains_Run1_t0.smearInputValuesChol(0,0);
   tagging_constrains_Run1_t0.smearInputValuesChol(1,1);
   tagging_constrains_Run1_t0.smearInputValuesChol(2,2);
   tagging_constrains_Run1_t0.smearInputValuesChol(3,3);

   tagging_constrains_Run1_t0.smearInputValuesChol(0,4);
   tagging_constrains_Run1_t0.smearInputValuesChol(1,5);


     throw "";

   TDecompChol tdc(*cov);
   tdc.Decompose();
   TMatrixD U = tdc.GetU();
   TMatrixD UT(TMatrixD::kTransposed,U);

   cov->Print();
   U.Print();
   UT.Print();

   int N = 10000;
   vector< vector<double> > results;

   for(int n = 0 ; n < N; n++){
	vector<double> v;
	vector<double> r;

	for(int i = 0 ; i < UT.GetNcols(); i++){
		double val = tagging_constrains_Run1_t0.smearInputValuesChol(i,n);
		v.push_back(val);
	}
	results.push_back(v);
   }

   vector<double> mean(results[0].size(),0.);
   vector<double> sigma(results[0].size(),0.);
   TMatrixD cov_col(results[0].size(),results[0].size());

   for(int n = 0 ; n < N; n++){
	for(int i = 0 ; i < results[0].size(); i++){
		mean[i] += results[n][i];
		sigma[i] += pow(results[n][i],2);	
	}
	for(int i = 0 ; i < results[0].size(); i++)for(int j = 0 ; j < results[0].size(); j++){
		cov_col(i,j) += results[n][i]*results[n][j]/N;
	}

   }

   for(int i = 0 ; i < results[0].size(); i++){
		   cout << "mean = " << mean[i]/N << " pm " << sqrt(sigma[i]/N) << endl;
   }

   cov_col.Print();


}


void animate(int step=0){
	TRandom3 ranLux;
	NamedParameter<int> RandomSeed("RandomSeed", 0);
	ranLux.SetSeed((int)RandomSeed);
	gRandom = &ranLux;
	
	NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
	
	FitParameter  r("r",1,0.,0.1);
	FitParameter  delta("delta",1,100.,1.);
	FitParameter  gamma("gamma",1,70,1.);
	FitParameter  k("k",1,1,1.);
	
	FitParameter  C("C",1,0.,0.1);
	FitParameter  D("D",1,0.,0.1);
	FitParameter  D_bar("D_bar",1,0.,0.1);
	FitParameter  S("S",1,0.,0.1);
	FitParameter  S_bar("S_bar",1,0.,0.1);
	
	C.setCurrentFitVal((1.-r*r)/(1.+r*r));

	D.setCurrentFitVal(-2.*r*k*cos((delta-gamma)/360.*2.*pi)/(1.+r*r));
	D_bar.setCurrentFitVal(-2.*r*k*cos((delta+gamma)/360.*2.*pi)/(1.+r*r));

	S.setCurrentFitVal(2.*r*k*sin((delta-gamma)/360.*2.*pi)/(1.+r*r));
	S_bar.setCurrentFitVal(-2.*r*k*sin((delta+gamma)/360.*2.*pi)/(1.+r*r));

	cout << "C = " << (double) C << endl;
	cout << "D = " << (double) D << endl;
	cout << "Dbar = " << (double) D_bar << endl;
	cout << "S = " << (double) S << endl;
	cout << "Sbar = " << (double) S_bar << endl << endl;

   	FitParameter  Gamma("Gamma",2,0.6629,0.0018);
	FitParameter  dGamma("dGamma",2,0.09,0.1);
	FitParameter  dm("dm",2,17.757,0.1);
	
	FitParameter  offset_mean_dt("offset_mean_dt",1,0,0.1);
	FitParameter  scale_mean_dt("scale_mean_dt",1,0,0.1);
	FitParameter  scale_mean_2_dt("scale_mean_2_dt",1,0,0.1);
	FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
	FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.,0.1);
	FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
	FitParameter  offset_sigma2_dt("offset_sigma2_dt",1,0.,0.1);
	FitParameter  scale_sigma2_dt("scale_sigma2_dt",1,1.,0.1);
	FitParameter  scale_sigma2_2_dt("scale_sigma2_2_dt",1,0.,0.1);
	FitParameter  offset_sigma3_dt("offset_sigma3_dt",1,0.,0.1);
	FitParameter  scale_sigma3_dt("scale_sigma3_dt",1,1.,0.1);
	FitParameter  scale_sigma3_2_dt("scale_sigma3_2_dt",1,0.,0.1);
	FitParameter  offset_f_dt("offset_f_dt",1,1,0.1);
	FitParameter  scale_f_dt("scale_f_dt",1,0.,0.1);
	FitParameter  scale_f_2_dt("scale_f_2_dt",1,0.,0.1);
	FitParameter  offset_f2_dt("offset_f2_dt",1,0.,0.1);
	FitParameter  scale_f2_dt("scale_f2_dt",1,0.,0.1);
	FitParameter  scale_f2_2_dt("scale_f2_2_dt",1,0.,0.1);

	FitParameter  p0_os("p0_os",1,0.,0.);
	FitParameter  p1_os("p1_os",1,1.,0.);
	FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
	FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
	FitParameter  avg_eta_os("avg_eta_os",1,0.,0.);
	FitParameter  tageff_os("tageff_os",1,1.,0.);
	FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
	FitParameter  p0_ss("p0_ss",1,0.,0.);
	FitParameter  p1_ss("p1_ss",1,1.,0.);
	FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
	FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
	FitParameter  avg_eta_ss("avg_eta_ss",1,0.,0.);
	FitParameter  tageff_ss("tageff_ss",1,1.,0.);
	FitParameter  tageff_asym_ss("tageff_asym_ss",1,0.,0.);
	FitParameter  production_asym("production_asym",1,0.,0.);
	FitParameter  detection_asym("detection_asym",1,0.1,0.);
	
	FitParameter  c0("c0",1,1,0.1);
	FitParameter  c1("c1",1,1,0.1);
	FitParameter  c2("c2",1,1,0.1);
	FitParameter  c3("c3",1,1,0.1);
	FitParameter  c4("c4",1,1,0.1);
	FitParameter  c5("c5",1,1,0.1);
	FitParameter  c6("c6",1,1,0.1);
	FitParameter  c7("c7",1,1,0.1);
	FitParameter  c8("c8",1,1,0.1);
	FitParameter  c9("c9",1,1,0.1);
	
	FullTimePdf t_pdf(C, D, D_bar, S, S_bar, k,
			Gamma, dGamma, dm
		      ,offset_mean_dt,scale_mean_dt,scale_mean_2_dt
                      ,offset_sigma_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,offset_sigma2_dt, scale_sigma2_dt, scale_sigma2_2_dt
                      ,offset_sigma3_dt, scale_sigma3_dt, scale_sigma3_2_dt
                      ,offset_f_dt, scale_f_dt, scale_f_2_dt
                      ,offset_f2_dt, scale_f2_dt, scale_f2_2_dt
			,c0, c1, c2 ,c3, c4, c5
			,c6, c7, c8, c9,
			p0_os, p1_os, delta_p0_os, delta_p1_os, 
			avg_eta_os, tageff_os, tageff_asym_os, 
			p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
			avg_eta_ss, tageff_ss, tageff_asym_ss, 
			production_asym, detection_asym, "comb" );

    	NamedParameter<int>  nBinst("nBinst", 50);
        NamedParameter<int>  nBinsAsym("nBinsAsym", 10);
	NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    	NamedParameter<double> min_TAU("min_TAU", 0.4);
   	NamedParameter<double> max_TAU("max_TAU", 10.);
   	NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);

 	min_TAU.setVal(-0.2);
	TGraph* graph = new TGraph(2);
	graph->SetPoint(1,-100,0);
	graph->SetPoint(2,100,0);
	graph->SetLineStyle(kDashed);
        
	TH1D* h_t_fit = new TH1D("h_t",";t [ps];Events ",nBinst,min_TAU,max_TAU);    

	TH1D* h_t_mixed_fit = new TH1D("h_t_mixed_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_unmixed_fit = new TH1D("h_t_unmixed_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_untagegged_fit = new TH1D("h_t_untagegged_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);

// 	double xBins[7] = {-1,1,3,5,7,9,11};
	TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",11,-0.5,10.5);
	TH1D* h_N_m_fit = new TH1D("h_N_m_fit",";t [ps];Events  ",11,-0.5,10.5);
// 	TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);
// 	TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);


	TH1D* h_t_fit_p = new TH1D("h_t_fit_p",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_m = new TH1D("h_t_fit_m",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);

	TH1D* h_t_fit_mp = new TH1D("h_t_fit_mp",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_0p = new TH1D("h_t_fit_0p",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_pp = new TH1D("h_t_fit_pp",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_mm = new TH1D("h_t_fit_mm",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_0m = new TH1D("h_t_fit_0m",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
	TH1D* h_t_fit_pm = new TH1D("h_t_fit_pm",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);

	TH1D* h_N_mixed_fit = new TH1D("h_N_mixed",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
	TH1D* h_N_unmixed_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed");	
	TH1D* h_N_mixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_p");
	TH1D* h_N_unmixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_p");
	TH1D* h_N_mixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_m");
	TH1D* h_N_unmixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_m");

	t_pdf.plotSpline();
	double N_OS_MC = 0;
	double N_SS_MC = 0;
	double N_OS_SS_MC = 0;
	double N_MC = 0;
	
	double w_OS_MC = 0;
	double w_SS_MC = 0;
	double w_OS_SS_MC = 0;
		
	double D_OS_MC = 0;
	double D_SS_MC = 0;
	double D_OS_SS_MC = 0;
	double D_comb_MC = 0;
	
	double N_OS_all_MC = 0;
	double N_SS_all_MC = 0;
	double w_OS_all_MC = 0;
	double w_SS_all_MC = 0;
	double D_OS_all_MC = 0;
	double D_SS_all_MC = 0;
	
	double N_Run1_t0_MC = 0;
	double N_Run1_t1_MC = 0;
	double N_Run2_t0_MC = 0;
	double N_Run2_t1_MC = 0;


	for(int n = 0 ; n < 5; n++){
		DalitzEventList sampleEvents = t_pdf.generateToys(1000000,1,0);
		//t_pdf.saveEventListToFile(sampleEvents);
	
		for(int i = 0; i < sampleEvents.size(); i++){
	
				DalitzEvent evt = sampleEvents[i];
				double t_MC = evt.getValueFromVector(0) ;
// 				double t_MC = abs(evt.getValueFromVector(0)) ;
				double dt_MC = evt.getValueFromVector(1) ;
				int f_MC = evt.getValueFromVector(2) ;
				int q_OS_MC = evt.getValueFromVector(3) ;
				double eta_OS_MC = evt.getValueFromVector(4) ;
				int q_SS_MC = evt.getValueFromVector(5);
				double eta_SS_MC = evt.getValueFromVector(6);
				int run_MC = evt.getValueFromVector(7) ;
				int trigger_MC = evt.getValueFromVector(8) ;
	
				double weight = 1; //t_pdf_Run1_t0.getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace();
				N_MC += weight;
			
				h_t_fit->Fill(t_MC,weight);
				
				int f_evt = f_MC;
				int q1 = q_OS_MC;
				int q2 = q_SS_MC;   
				int q_eff = 0;
				double w_eff = 0.5;
						
				std::pair<double, double> calibrated_mistag_os;
				std::pair<double, double> calibrated_mistag_ss;
				calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eta_OS_MC);
				calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eta_SS_MC);        
				
				double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
				double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
				
				if( p/(p+p_bar) > 0.5 ){ 
					q_eff = 1;
					w_eff = 1-p/(p+p_bar);
				}
				else if( p/(p+p_bar) < 0.5 ){
					q_eff = -1;
					w_eff = p/(p+p_bar);
				}
				
				if(q1 != 0 && q2 != 0){
					if(q_eff != 0){
						N_OS_SS_MC += weight;
						w_OS_SS_MC += w_eff * weight;
						D_OS_SS_MC += pow(1.-2.*w_eff,2)* weight;
					}
				}
				else if( q1 != 0){
					//q_eff = q1;  flip tag ???
					N_OS_MC +=weight;
					w_OS_MC += w_eff * weight; 
					D_OS_MC += pow(1.-2.*w_eff,2)* weight;
				}
				else if( q2 != 0){
					//q_eff = q2;
					N_SS_MC += weight;
					w_SS_MC += w_eff * weight; 
					D_SS_MC += pow(1.-2.*w_eff,2)* weight; 
				} 
			
				D_comb_MC += pow(1.-2.*w_eff,2)* weight;
				
				if(q1 != 0) N_OS_all_MC += weight;
				if(q2 != 0)N_SS_all_MC += weight;
				
				if(q1>0){
						w_OS_all_MC +=  calibrated_mistag_os.first * weight;
						D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.first,2)* weight;
				} 
				else if(q1<0){
						w_OS_all_MC +=  calibrated_mistag_os.second * weight;	
						D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.second,2)* weight;
				}
			
				if(q2>0){
						w_SS_all_MC +=  calibrated_mistag_ss.first * weight;
						D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.first,2)* weight;
				} 
				else if(q2<0){
						w_SS_all_MC +=  calibrated_mistag_ss.second * weight;	
						D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.second,2)* weight;
				}
			
				if((string)channel=="signal"){
	
					if(f_evt == 1){
						 h_t_fit_p->Fill(t_MC,weight);
					  	 h_N_p_fit->Fill(t_MC,weight);
					}
					else{
						 h_t_fit_m->Fill(t_MC,weight);
   						 h_N_m_fit->Fill(t_MC,weight);
					}
					if(q_eff==-1 && f_evt == 1){
						h_t_fit_mp->Fill(t_MC,weight);
						h_N_mixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
					else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight);
					else if(q_eff==1 && f_evt == 1){
						h_t_fit_pp->Fill(t_MC,weight);
						h_N_unmixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
					else if(q_eff==-1 && f_evt == -1){
						h_t_fit_mm->Fill(t_MC,weight);
						h_N_unmixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
					else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(t_MC,weight);
					else if(q_eff==1 && f_evt == -1){
						h_t_fit_pm->Fill(t_MC,weight);
						h_N_mixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
				}
				else {   
					if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
					else if(q_eff*f_evt > 0  ){
						h_t_mixed_fit->Fill(t_MC,weight);
						h_N_mixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
					else{ 
						h_t_unmixed_fit->Fill(t_MC,weight);
						h_N_unmixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
					}
				}
		}

	}
	
	cout << "N_MC = " << N_MC << endl << endl;
	cout << "N_OS = " << N_OS_all_MC << endl;
	cout << "N_SS = " << N_SS_all_MC << endl << endl;
	
	cout << "eff_OS = " <<(N_OS_all_MC)/N_MC << endl;
	cout << "eff_SS = " <<(N_SS_all_MC)/N_MC << endl;
		
	cout << "Tagging perfromance " << endl << endl;        
	cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
	
	cout << "OS  | " << (N_OS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_OS_all_MC)/(N_OS_MC+N_OS_SS_MC) << " | " << D_OS_all_MC/N_MC << endl;
	cout << "SS  | " << (N_SS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_SS_all_MC)/(N_SS_MC+N_OS_SS_MC) << " | " << D_SS_all_MC/N_MC << endl << endl;
	
	cout << "OS only  | " << N_OS_MC/N_MC << " | " <<  w_OS_MC/N_OS_MC << " | " << N_OS_MC/N_MC * D_OS_MC/N_OS_MC << endl;
	cout << "SS only  | " << N_SS_MC/N_MC << " | " <<  w_SS_MC/N_SS_MC << " | " << N_SS_MC/N_MC * D_SS_MC/N_SS_MC << endl;
	cout << "OS+SS    | " << N_OS_SS_MC/N_MC << " | " <<  w_OS_SS_MC/N_OS_SS_MC << " | " << N_OS_SS_MC/N_MC * D_OS_SS_MC/N_OS_SS_MC << endl;
	cout << "Combined | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC << " | "<<  (w_OS_MC+w_SS_MC+w_OS_SS_MC)/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << " | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC * D_comb_MC/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << endl << endl ;

	/// Plot
	TCanvas* c= new TCanvas("");

	TLegend leg(0.6,0.6,0.925,0.9,"");
	leg.SetLineStyle(0);
	leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.06);
	leg.SetTextAlign(12);
// 	leg.AddEntry((TObject*)0,"#font[22]{LHCb unofficial}","");
	leg.SetNColumns(2);

	stringstream ss ;
	TString label_r = "r = ";
	ss << std::fixed << std::setprecision(1) << (double)r;
	label_r += ss.str();

	ss.str("");
	TString label_k = "#kappa = ";
	ss << std::fixed << std::setprecision(1) << (double)k;
	label_k += ss.str();

	ss.str("");
	TString label_d = "#delta = ";
	ss << std::fixed << std::setprecision(0) << (double)delta;
	label_d += ss.str();

	ss.str("");
	TString label_g = "#gamma = ";
	ss << std::fixed << std::setprecision(0) << (double)gamma;
	label_g += ss.str();


	h_t_fit->SetLineColor(kBlue);
	h_t_fit->SetLineWidth(5);
	h_t_fit->SetMarkerColor(kBlue); 
	h_t_fit->Draw("histc");
// 	leg.Draw();

	c->Print(((string)OutputDir+"h_t.eps").c_str());
	c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_t_log.eps").c_str());
	c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
	gPad->SetLogy(0);

	if((string)channel=="norm"){
		h_t_mixed_fit->SetMarkerColor(kRed); 
		h_t_mixed_fit->SetLineColor(kRed);
		h_t_mixed_fit->DrawNormalized("histc",1);
		
		h_t_unmixed_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_fit->SetLineColor(kBlue);
		h_t_unmixed_fit->DrawNormalized("histcsame",1);
				
		h_t_untagegged_fit->SetMarkerColor(kGreen); 
		h_t_untagegged_fit->SetLineColor(kGreen);
		//h_t_untagegged_fit->DrawNormalized("histcsame",1);
		
		c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
		c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed.pdf").c_str());
		
		TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
		h_asym_fit->SetMinimum(-1);
		h_asym_fit->SetMaximum(1);
		h_asym_fit->SetLineColor(kRed);
		h_asym_fit->Draw("histc");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());
		c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
	}
	
	else{		

// 		h_t_fit_p->SetAxisRange(0.4, max_TAU*2,"X");
// 		h_t_fit_p->Smooth(10,"R");
// 		h_t_fit_p->SetAxisRange(min_TAU, max_TAU,"X");

		h_t_fit_m->Scale(1./10.);
		h_t_fit_p->Scale(1./10.);
		h_t_fit_p->SetMinimum(0);
		h_t_fit_p->SetMaximum(4000);
		h_t_fit_p->SetLineWidth(4);
		h_t_fit_p->SetLineColor(kBlack);
		h_t_fit_p->Draw("e1");

		h_t_fit_mp->Scale(1./10.);
// 		h_t_fit_mp->SetBinContent(h_t_fit_pp->FindBin(0.),0.);
// 		h_t_fit_mp->Smooth(10);
		h_t_fit_mp->SetLineWidth(4);
		h_t_fit_mp->SetMarkerColor(kBlue); 
		h_t_fit_mp->SetLineColor(kBlue);
		h_t_fit_mp->Draw("histcsame");


		h_t_fit_pp->Scale(1./10.);
// 		h_t_fit_pp->SetBinContent(h_t_fit_pp->FindBin(0.),h_t_fit_p->GetBinContent(h_t_fit_pp->FindBin(0.)));
// 		h_t_fit_pp->Smooth(10);
		h_t_fit_pp->SetLineWidth(4);		
		h_t_fit_pp->SetMarkerColor(kRed); 
		h_t_fit_pp->SetLineColor(kRed);
		h_t_fit_pp->Draw("histcsame");

// 		leg.AddEntry(h_t_fit_p,"Total","l");
		TLegendEntry* le = leg.AddEntry(h_t_fit_pp,"B_{s}#rightarrow f","l");
		le->SetTextColor(kRed);    				
		le = leg.AddEntry(h_t_fit_mp,"#bar{B_{s}}#rightarrow f","l");
		le->SetTextColor(kBlue);    	
		leg.AddEntry((TObject*)0,label_r,"");
		leg.AddEntry((TObject*)0,label_k,"");
		
		leg.AddEntry((TObject*)0,label_d,"");
		leg.AddEntry((TObject*)0,label_g,"");
		leg.Draw();

		c->Print(((string)OutputDir+"h_t_mixed_p.C").c_str());
		c->Print(((string)OutputDir+"h_t_mixed_p.root").c_str());
		c->Print(((string)OutputDir+"h_t_mixed_p.eps").c_str());


// 		h_t_fit_mp->SetMaximum(1600);
		h_t_fit_mp->Draw("histc");
		h_t_fit_pp->Draw("histcsame");

		h_t_fit_pm->Scale(1./10.);
// 		h_t_fit_pm->Smooth(10);
		h_t_fit_pm->SetLineWidth(4);
		h_t_fit_pm->SetLineStyle(kDashed); 			
		h_t_fit_pm->SetMarkerColor(kBlue+1); 
		h_t_fit_pm->SetLineColor(kBlue+1);
		h_t_fit_pm->Draw("histcsame");

 		h_t_fit_mm->Scale(1./10.);
// 		h_t_fit_mm->Smooth(10);
		h_t_fit_mm->SetLineWidth(4);
		h_t_fit_mm->SetLineStyle(kDashed); 					
		h_t_fit_mm->SetMarkerColor(kRed+1); 
		h_t_fit_mm->SetLineColor(kRed+1);
		h_t_fit_mm->Draw("histcsame");

		TLegend leg3(0.6,0.6,0.925,0.9,"");
		leg3.SetLineStyle(0);
		leg3.SetLineColor(0);
		leg3.SetFillColor(0);
		leg3.SetTextFont(22);
		leg3.SetTextColor(1);
		leg3.SetTextSize(0.06);
		leg3.SetTextAlign(12);
		leg3.SetNColumns(2);

		le = leg3.AddEntry(h_t_fit_pp,"B_{s}#rightarrowf","l");
		le->SetTextColor(kRed);    				
		le = leg3.AddEntry(h_t_fit_mp,"#bar{B_{s}}#rightarrowf","l");
		le->SetTextColor(kBlue);    

		le = leg3.AddEntry(h_t_fit_mm,"#bar{B_{s}}#rightarrow#bar{f}","l");
		le->SetTextColor(kRed+1);    				
		le = leg3.AddEntry(h_t_fit_pm,"B_{s}#rightarrow#bar{f}","l");
		le->SetTextColor(kBlue+1);    

		leg3.AddEntry((TObject*)0,label_r,"");
		leg3.AddEntry((TObject*)0,label_k,"");
		
		leg3.AddEntry((TObject*)0,label_d,"");
		leg3.AddEntry((TObject*)0,label_g,"");
		
		leg3.Draw();

		c->Print(((string)OutputDir+"h_t_mixed.C").c_str());
		c->Print(((string)OutputDir+"h_t_mixed.root").c_str());
		c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
		c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());



		TH1D* h_asymCP_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);	
		h_asymCP_fit->SetLineColor(kBlack);
		h_asymCP_fit->SetLineWidth(5);
// 		h_asymCP_fit->SetMinimum(-1.);
// 		h_asymCP_fit->SetMaximum(1.);
		h_asymCP_fit->SetYTitle("A_{CP}");
		h_asymCP_fit->Draw("histc");

		TH1D* h_asymCP2_fit = (TH1D*) h_N_mixed_m_fit->GetAsymmetry(h_N_mixed_p_fit);	
		h_asymCP2_fit->SetLineColor(kRed);
		h_asymCP2_fit->SetLineWidth(5);
// 		h_asymCP_fit->SetMinimum(-1.);
// 		h_asymCP_fit->SetMaximum(1.);
		h_asymCP2_fit->Draw("histcsame");
		graph->Draw("same");

		c->Print(((string)OutputDir+"h_asymCP.eps").c_str());	
		c->Print(((string)OutputDir+"h_asymCP.C").c_str());	
		c->Print(((string)OutputDir+"h_asymCP.root").c_str());	


		TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);	
		h_asym_p_fit->SetLineColor(kBlack);
		h_asym_p_fit->SetLineWidth(5);
		h_asym_p_fit->SetMinimum(-1.);
		h_asym_p_fit->SetMaximum(1.);
		h_asym_p_fit->Draw("histc");
		
		TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);	
		h_asym_m_fit->SetLineWidth(5);
		h_asym_m_fit->SetLineColor(kGreen+3);
		h_asym_m_fit->SetLineStyle(kDashed);
		h_asym_m_fit->Draw("histcsame");

		TLegend leg2(0.45,0.75,0.6,0.92,"");
		leg2.SetLineStyle(0);
		leg2.SetLineColor(0);
		leg2.SetFillColor(0);
		leg2.SetTextFont(22);
		leg2.SetTextColor(1);
		leg2.SetTextSize(0.06);
		leg2.SetTextAlign(12);

		le = leg2.AddEntry(h_asym_p_fit,"f=D_{s}^{-}K^{+}#pi#pi","l");
		le->SetTextColor(kBlack);    				
		le = leg2.AddEntry(h_asym_m_fit,"#bar{f}=D_{s}^{+}K^{-}#pi#pi","l");
		le->SetTextColor(kGreen+3);    				
		leg2.Draw();

		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());	
		c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
		c->Print(((string)OutputDir+"h_asym.C").c_str());	
		c->Print(((string)OutputDir+"h_asym.root").c_str());	

		h_asym_p_fit->Add(h_asym_m_fit,-1.);	
		h_asym_p_fit->Smooth(2);
		h_asym_p_fit->SetLineColor(kMagenta+3);
		h_asym_p_fit->SetLineWidth(5);
		h_asym_p_fit->SetMinimum(-0.4);
		h_asym_p_fit->SetMaximum(0.4);
		h_asym_p_fit->SetYTitle("A_{CP}");
// 		h_asym_p_fit->SetYTitle("#DeltaA_{flav}");
		h_asym_p_fit->Draw("histc");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_delta_asym.eps").c_str());

		h_t_fit_p->Draw("e1");
		h_t_fit_m->Draw("histcsame");
		c->Print(((string)OutputDir+"h_t_f.eps").c_str());

// 		h_t_fit_p->Rebin(20);
// 		h_t_fit_m->Rebin(20);

		TH1D* h_asym2_fit = (TH1D*) h_N_m_fit->GetAsymmetry(h_N_p_fit);	
		h_asym2_fit->SetLineColor(kBlack);
		h_asym2_fit->SetLineWidth(5);
/*		h_asym2_fit->Rebin(40);*/
// 		h_asym2_fit->SetMinimum(-1.);
// 		h_asym2_fit->SetMaximum(1.);

// 		h_asym2_fit->GetXaxis()->SetRangeUser(0.,10.);
// 		h_asym2_fit->Draw("e1");

	 	TH1D* h_N_fake = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);
		h_N_fake->SetMaximum(h_asym2_fit->GetMaximum());
		h_N_fake->SetMinimum(h_asym2_fit->GetMinimum());
		h_N_fake->Draw("e");
		h_asym2_fit->SetLineColor(kRed);
		h_asym2_fit->Draw("histcsame");
// 		h_asym2_fit->Rebin(2);
// 		h_asym2_fit->SetLineColor(kBlue);
// 		h_asym2_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym2.eps").c_str());
		c->Print(((string)OutputDir+"h_asym2.C").c_str());

		h_N_unmixed_p_fit->Add(h_N_mixed_m_fit,1.);	
		h_N_unmixed_m_fit->Add(h_N_mixed_p_fit,1.);	
		TH1D* h_asym3_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);	
		h_asym3_fit->Draw("histc");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym3.eps").c_str());


	}


	
}


void animate2(int step=0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);
    FitParameter  k("k",1,1,1.);
    
    FitParameter  C("C",1,0.,0.1);
    FitParameter  D("D",1,0.,0.1);
    FitParameter  D_bar("D_bar",1,0.,0.1);
    FitParameter  S("S",1,0.,0.1);
    FitParameter  S_bar("S_bar",1,0.,0.1);
    
    C.setCurrentFitVal((1.-r*r)/(1.+r*r));
    
    D.setCurrentFitVal(-2.*r*k*cos((delta-gamma)/360.*2.*pi)/(1.+r*r));
    D_bar.setCurrentFitVal(-2.*r*k*cos((delta+gamma)/360.*2.*pi)/(1.+r*r));
    
    S.setCurrentFitVal(2.*r*k*sin((delta-gamma)/360.*2.*pi)/(1.+r*r));
    S_bar.setCurrentFitVal(-2.*r*k*sin((delta+gamma)/360.*2.*pi)/(1.+r*r));
    
    cout << "C = " << (double) C << endl;
    cout << "D = " << (double) D << endl;
    cout << "Dbar = " << (double) D_bar << endl;
    cout << "S = " << (double) S << endl;
    cout << "Sbar = " << (double) S_bar << endl << endl;
    
    FitParameter  Gamma("Gamma",2,0.6629,0.0018);
    FitParameter  dGamma("dGamma",2,0.09,0.1);
    FitParameter  dm("dm",2,17.757,0.1);

    FitParameter  offset_mean_dt("offset_mean_dt",1,0,0.1);
    FitParameter  scale_mean_dt("scale_mean_dt",1,0,0.1);
    FitParameter  scale_mean_2_dt("scale_mean_2_dt",1,0,0.1);
    FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
    FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.,0.1);
    FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
    FitParameter  offset_sigma2_dt("offset_sigma2_dt",1,0.,0.1);
    FitParameter  scale_sigma2_dt("scale_sigma2_dt",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt("scale_sigma2_2_dt",1,0.,0.1);
    FitParameter  offset_sigma3_dt("offset_sigma3_dt",1,0.,0.1);
    FitParameter  scale_sigma3_dt("scale_sigma3_dt",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt("scale_sigma3_2_dt",1,0.,0.1);
    FitParameter  offset_f_dt("offset_f_dt",1,1,0.1);
    FitParameter  scale_f_dt("scale_f_dt",1,0.,0.1);
    FitParameter  scale_f_2_dt("scale_f_2_dt",1,0.,0.1);
    FitParameter  offset_f2_dt("offset_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_dt("scale_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_2_dt("scale_f2_2_dt",1,0.,0.1);

    FitParameter  p0_os("p0_os",1,0.,0.);
    FitParameter  p1_os("p1_os",1,1.,0.);
    FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
    FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
    FitParameter  avg_eta_os("avg_eta_os",1,0.,0.);
    FitParameter  tageff_os("tageff_os",1,1.,0.);
    FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
    FitParameter  p0_ss("p0_ss",1,0.,0.);
    FitParameter  p1_ss("p1_ss",1,1.,0.);
    FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
    FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
    FitParameter  avg_eta_ss("avg_eta_ss",1,0.,0.);
    FitParameter  tageff_ss("tageff_ss",1,1.,0.);
    FitParameter  tageff_asym_ss("tageff_asym_ss",1,0.,0.);
    FitParameter  production_asym("production_asym",1,0.,0.);
    FitParameter  detection_asym("detection_asym",1,0.1,0.);
    
    FitParameter  c0("c0",1,1,0.1);
    FitParameter  c1("c1",1,1,0.1);
    FitParameter  c2("c2",1,1,0.1);
    FitParameter  c3("c3",1,1,0.1);
    FitParameter  c4("c4",1,1,0.1);
    FitParameter  c5("c5",1,1,0.1);
    FitParameter  c6("c6",1,1,0.1);
    FitParameter  c7("c7",1,1,0.1);
    FitParameter  c8("c8",1,1,0.1);
    FitParameter  c9("c9",1,1,0.1);
    
    FullTimePdf t_pdf(C, D, D_bar, S, S_bar, k,
                      Gamma, dGamma, dm
		      ,offset_mean_dt,scale_mean_dt,scale_mean_2_dt
                      ,offset_sigma_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,offset_sigma2_dt, scale_sigma2_dt, scale_sigma2_2_dt
                      ,offset_sigma3_dt, scale_sigma3_dt, scale_sigma3_2_dt
                      ,offset_f_dt, scale_f_dt, scale_f_2_dt
                      ,offset_f2_dt, scale_f2_dt, scale_f2_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, "comb" );
    
    NamedParameter<int>  nBinst("nBinst", 50);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 10);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    
    min_TAU.setVal(-0.2);
    TGraph* graph = new TGraph(2);
    graph->SetPoint(1,-100,0);
    graph->SetPoint(2,100,0);
    graph->SetLineStyle(kDashed);
    
    TH1D* h_t_fit = new TH1D("h_t",";t [ps];Events ",nBinst,min_TAU,max_TAU);    
    
    TH1D* h_t_mixed_fit = new TH1D("h_t_mixed_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_fit = new TH1D("h_t_unmixed_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged_fit = new TH1D("h_t_untagegged_fit",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    
    //     double xBins[7] = {-1,1,3,5,7,9,11};
    TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",11,-0.5,10.5);
    TH1D* h_N_m_fit = new TH1D("h_N_m_fit",";t [ps];Events  ",11,-0.5,10.5);
    //     TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);
    //     TH1D* h_N_p_fit = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);
    
    
    TH1D* h_t_fit_p = new TH1D("h_t_fit_p",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_m = new TH1D("h_t_fit_m",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    
    TH1D* h_t_fit_mp = new TH1D("h_t_fit_mp",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0p = new TH1D("h_t_fit_0p",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pp = new TH1D("h_t_fit_pp",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_mm = new TH1D("h_t_fit_mm",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0m = new TH1D("h_t_fit_0m",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pm = new TH1D("h_t_fit_pm",";t [ps];Events  ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    
    TH1D* h_N_mixed_fit = new TH1D("h_N_mixed",";t modulo (2#pi/#Deltam_{s}) [ps]; A_{mix} ",nBinsAsym,0.,2.*pi/dm);
    TH1D* h_N_unmixed_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed");    
    TH1D* h_N_mixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_p");
    TH1D* h_N_unmixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_p");
    TH1D* h_N_mixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_m");
    TH1D* h_N_unmixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_m");
    
    t_pdf.plotSpline();
    double N_OS_MC = 0;
    double N_SS_MC = 0;
    double N_OS_SS_MC = 0;
    double N_MC = 0;
    
    double w_OS_MC = 0;
    double w_SS_MC = 0;
    double w_OS_SS_MC = 0;
    
    double D_OS_MC = 0;
    double D_SS_MC = 0;
    double D_OS_SS_MC = 0;
    double D_comb_MC = 0;
    
    double N_OS_all_MC = 0;
    double N_SS_all_MC = 0;
    double w_OS_all_MC = 0;
    double w_SS_all_MC = 0;
    double D_OS_all_MC = 0;
    double D_SS_all_MC = 0;
    
    double N_Run1_t0_MC = 0;
    double N_Run1_t1_MC = 0;
    double N_Run2_t0_MC = 0;
    double N_Run2_t1_MC = 0;
    
    
    TF1 *N_f = new TF1("N_f"," [5] * [6] * (cosh([3]/2*x)+ [0] * cos([4]*x) + [1]  * sinh([3]/2*x) - [2] * sin([4]*x) )",min_TAU,max_TAU);
    N_f->SetParameter(0,(double) C);
    N_f->SetParameter(1,(double) D);
    N_f->SetParameter(2,(double) S);
    N_f->SetParameter(3,(double) dGamma);
    N_f->SetParameter(4,(double) dm);
    N_f->SetParameter(5,1- (double) production_asym);
    N_f->SetParameter(6,1- (double) detection_asym);
    
    TF1 *Nbar_f = new TF1("Nbar_f","[5] * [6] * ( cosh([3]/2*x) - [0] * cos([4]*x) + [1]  * sinh([3]/2*x) + [2]* sin([4]*x) )",min_TAU,max_TAU);
    Nbar_f->SetParameter(0,(double) C);
    Nbar_f->SetParameter(1,(double) D);
    Nbar_f->SetParameter(2,(double) S);
    Nbar_f->SetParameter(3,(double) dGamma);
    Nbar_f->SetParameter(4,(double) dm);
    Nbar_f->SetParameter(5,1+ (double) production_asym);
    Nbar_f->SetParameter(6,1- (double) detection_asym);
    
    TF1 *N_fbar = new TF1("N_fbar","[5] * [6] * ( cosh([3]/2*x)- [0] * cos([4]*x) + [1]  * sinh([3]/2*x) - [2] *sin([4]*x) )",min_TAU,max_TAU);
    N_fbar->SetParameter(0,(double) C);
    N_fbar->SetParameter(1,(double) D_bar);
    N_fbar->SetParameter(2,(double) S_bar);
    N_fbar->SetParameter(3,(double) dGamma);
    N_fbar->SetParameter(4,(double) dm);
    N_fbar->SetParameter(5,1- (double) production_asym);
    N_fbar->SetParameter(6,1+ (double) detection_asym);
    
    TF1 *Nbar_fbar = new TF1("Nbar_fbar","[5] * [6] * ( cosh([3]/2*x)+ [0] * cos([4]*x) + [1]  * sinh([3]/2*x) + [2] * sin([4]*x) )",min_TAU,max_TAU);
    Nbar_fbar->SetParameter(0,(double) C);
    Nbar_fbar->SetParameter(1,(double) D_bar);
    Nbar_fbar->SetParameter(2,(double) S_bar);
    Nbar_fbar->SetParameter(3,(double) dGamma);
    Nbar_fbar->SetParameter(4,(double) dm);
    Nbar_fbar->SetParameter(5,1+ (double) production_asym);
    Nbar_fbar->SetParameter(6,1+ (double) detection_asym);
    
    TF1 *A_mix = new TF1("A_mix", "fmod( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) , 2*3.141/17.757)",0,2*3.141/17.757);
    TF1 *A_q= new TF1("A_q", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);

    TF1 *A_mix_true = new TF1("A_mix", "fmod( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) , 2*3.141/17.757)",0,2*3.141/17.757);
    TF1 *A_q_true= new TF1("A_q", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);

    
    N_f->SetParameter(5,1- 0.05);
    N_fbar->SetParameter(5,1- 0.05);
    Nbar_f->SetParameter(5,1+ 0.05);
    Nbar_fbar->SetParameter(5,1+ 0.05);    
    TF1 *A_mix_2 = new TF1("A_mix_2","( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) )",0,2*3.141/17.757);
    TF1 *A_q_2= new TF1("A_q_2", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);
    
    N_f->SetParameter(5,1- 0.1);
    N_fbar->SetParameter(5,1- 0.1);
    Nbar_f->SetParameter(5,1+ 0.1);
    Nbar_fbar->SetParameter(5,1+ 0.1);    
    TF1 *A_mix_3 = new TF1("A_mix_3", "( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) )",0,2*3.141/17.757);
    TF1 *A_q_3= new TF1("A_q_3", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);
    
    N_f->SetParameter(5,1);
    N_fbar->SetParameter(5,1.);
    Nbar_f->SetParameter(5,1.);
    Nbar_fbar->SetParameter(5,1.); 
    N_f->SetParameter(6,1-0.05);
    N_fbar->SetParameter(6,1+0.05);
    Nbar_f->SetParameter(6,1-0.05);
    Nbar_fbar->SetParameter(6,1+0.05); 
    TF1 *A_mix_4 = new TF1("A_mix_4", "fmod( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) , 2*3.141/17.757)",0,2*3.141/17.757);
    TF1 *A_q_4= new TF1("A_q_4", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);
    
    N_f->SetParameter(6,1-0.1);
    N_fbar->SetParameter(6,1+0.1);
    Nbar_f->SetParameter(6,1-0.1);
    Nbar_fbar->SetParameter(6,1+0.1);  
    TF1 *A_mix_5 = new TF1("A_mix_5", "fmod( ( (N_f+N_fbar) - (Nbar_fbar + Nbar_f) )/((N_f+N_fbar) + (Nbar_fbar + Nbar_f) ) , 2*3.141/17.757)",0,2*3.141/17.757);
    TF1 *A_q_5= new TF1("A_q_5", "( (N_f+Nbar_f) - (N_fbar+Nbar_fbar) )/((N_f+Nbar_f) + (N_fbar+Nbar_fbar))",0,max_TAU);
    
    
    gStyle->SetTitleOffset(0.95,"Y");

    A_mix->SetTitle("; t modulo (2#pi/#Deltam_{s}) [ps]; A_{raw}^{#LTf#GT}");
    A_mix->SetMinimum(-0.5);
    A_mix->SetMaximum(0.5);
    A_mix->SetNpx(500);
    A_mix_2->SetNpx(500);
    A_mix_3->SetNpx(500);
    A_mix_4->SetNpx(500);
    A_mix_5->SetNpx(500);

    A_q->SetTitle("; t [ps] ; A_{raw}^{#LTi#GT}");
    A_q->SetMinimum(-0.1);
    A_q->SetMaximum(0.1);
    A_q->SetNpx(1000);
    A_q_2->SetNpx(1000);
    A_q_3->SetNpx(1000);
    A_q_4->SetNpx(1000);
    A_q_5->SetNpx(1000);
    
    
    A_mix_true->SetTitle("; t modulo (2#pi/#Deltam_{s}) [ps]; A_{CP}^{#LTf#GT}");
    A_mix_true->SetMinimum(-0.5);
    A_mix_true->SetMaximum(0.5);
    A_mix_true->SetNpx(500);
    A_q_true->SetTitle("; t [ps] ; A_{CP}^{#LTi#GT}");
    A_q_true->SetMinimum(-0.1);
    A_q_true->SetMaximum(0.1);
    A_q_true->SetNpx(1000);
    
    for(int n = 0 ; n < 5; n++){
        DalitzEventList sampleEvents = t_pdf.generateToys(1000000,1,0);
        //t_pdf.saveEventListToFile(sampleEvents);
        
        for(int i = 0; i < sampleEvents.size(); i++){
            
            DalitzEvent evt = sampleEvents[i];
            double t_MC = evt.getValueFromVector(0) ;
            //                 double t_MC = abs(evt.getValueFromVector(0)) ;
            double dt_MC = evt.getValueFromVector(1) ;
            int f_MC = evt.getValueFromVector(2) ;
            int q_OS_MC = evt.getValueFromVector(3) ;
            double eta_OS_MC = evt.getValueFromVector(4) ;
            int q_SS_MC = evt.getValueFromVector(5);
            double eta_SS_MC = evt.getValueFromVector(6);
            int run_MC = evt.getValueFromVector(7) ;
            int trigger_MC = evt.getValueFromVector(8) ;
            
            double weight = 1; //t_pdf_Run1_t0.getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace();
            N_MC += weight;
            
            h_t_fit->Fill(t_MC,weight);
            
            int f_evt = f_MC;
            int q1 = q_OS_MC;
            int q2 = q_SS_MC;   
            int q_eff = 0;
            double w_eff = 0.5;
            
            std::pair<double, double> calibrated_mistag_os;
            std::pair<double, double> calibrated_mistag_ss;
            calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eta_OS_MC);
            calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eta_SS_MC);        
            
            double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
            double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
            
            if( p/(p+p_bar) > 0.5 ){ 
                q_eff = 1;
                w_eff = 1-p/(p+p_bar);
            }
            else if( p/(p+p_bar) < 0.5 ){
                q_eff = -1;
                w_eff = p/(p+p_bar);
            }
            
            if(q1 != 0 && q2 != 0){
                if(q_eff != 0){
                    N_OS_SS_MC += weight;
                    w_OS_SS_MC += w_eff * weight;
                    D_OS_SS_MC += pow(1.-2.*w_eff,2)* weight;
                }
            }
            else if( q1 != 0){
                //q_eff = q1;  flip tag ???
                N_OS_MC +=weight;
                w_OS_MC += w_eff * weight; 
                D_OS_MC += pow(1.-2.*w_eff,2)* weight;
            }
            else if( q2 != 0){
                //q_eff = q2;
                N_SS_MC += weight;
                w_SS_MC += w_eff * weight; 
                D_SS_MC += pow(1.-2.*w_eff,2)* weight; 
            } 
            
            D_comb_MC += pow(1.-2.*w_eff,2)* weight;
            
            if(q1 != 0) N_OS_all_MC += weight;
            if(q2 != 0)N_SS_all_MC += weight;
            
            if(q1>0){
                w_OS_all_MC +=  calibrated_mistag_os.first * weight;
                D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.first,2)* weight;
            } 
            else if(q1<0){
                w_OS_all_MC +=  calibrated_mistag_os.second * weight;    
                D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.second,2)* weight;
            }
            
            if(q2>0){
                w_SS_all_MC +=  calibrated_mistag_ss.first * weight;
                D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.first,2)* weight;
            } 
            else if(q2<0){
                w_SS_all_MC +=  calibrated_mistag_ss.second * weight;    
                D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.second,2)* weight;
            }
            
            if((string)channel=="signal"){
                
                if(f_evt == 1){
                    h_t_fit_p->Fill(t_MC,weight);
                    h_N_p_fit->Fill(t_MC,weight);
                }
                else{
                    h_t_fit_m->Fill(t_MC,weight);
                    h_N_m_fit->Fill(t_MC,weight);
                }
                if(q_eff==-1 && f_evt == 1){
                    h_t_fit_mp->Fill(t_MC,weight);
                    h_N_mixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
                else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight);
                else if(q_eff==1 && f_evt == 1){
                    h_t_fit_pp->Fill(t_MC,weight);
                    h_N_unmixed_p_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
                else if(q_eff==-1 && f_evt == -1){
                    h_t_fit_mm->Fill(t_MC,weight);
                    h_N_unmixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
                else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(t_MC,weight);
                else if(q_eff==1 && f_evt == -1){
                    h_t_fit_pm->Fill(t_MC,weight);
                    h_N_mixed_m_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
            }
            else {   
                if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
                else if(q_eff*f_evt > 0  ){
                    h_t_mixed_fit->Fill(t_MC,weight);
                    h_N_mixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
                else{ 
                    h_t_unmixed_fit->Fill(t_MC,weight);
                    h_N_unmixed_fit->Fill(fmod(t_MC,2.*pi/dm),weight);
                }
            }
        }
        
    }
    
    cout << "N_MC = " << N_MC << endl << endl;
    cout << "N_OS = " << N_OS_all_MC << endl;
    cout << "N_SS = " << N_SS_all_MC << endl << endl;
    
    cout << "eff_OS = " <<(N_OS_all_MC)/N_MC << endl;
    cout << "eff_SS = " <<(N_SS_all_MC)/N_MC << endl;
    
    cout << "Tagging perfromance " << endl << endl;        
    cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
    
    cout << "OS  | " << (N_OS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_OS_all_MC)/(N_OS_MC+N_OS_SS_MC) << " | " << D_OS_all_MC/N_MC << endl;
    cout << "SS  | " << (N_SS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_SS_all_MC)/(N_SS_MC+N_OS_SS_MC) << " | " << D_SS_all_MC/N_MC << endl << endl;
    
    cout << "OS only  | " << N_OS_MC/N_MC << " | " <<  w_OS_MC/N_OS_MC << " | " << N_OS_MC/N_MC * D_OS_MC/N_OS_MC << endl;
    cout << "SS only  | " << N_SS_MC/N_MC << " | " <<  w_SS_MC/N_SS_MC << " | " << N_SS_MC/N_MC * D_SS_MC/N_SS_MC << endl;
    cout << "OS+SS    | " << N_OS_SS_MC/N_MC << " | " <<  w_OS_SS_MC/N_OS_SS_MC << " | " << N_OS_SS_MC/N_MC * D_OS_SS_MC/N_OS_SS_MC << endl;
    cout << "Combined | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC << " | "<<  (w_OS_MC+w_SS_MC+w_OS_SS_MC)/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << " | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC * D_comb_MC/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << endl << endl ;
    
    /// Plot
    TCanvas* c= new TCanvas("");
    
    TLegend leg(0.5,0.6,0.925,0.9,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(22);
    leg.SetTextColor(1);
    leg.SetTextSize(0.06);
    leg.SetTextAlign(12);
    //     leg.AddEntry((TObject*)0,"#font[22]{LHCb unofficial}","");
    leg.SetNColumns(2);
    
    stringstream ss ;
    TString label_r = "r = ";
    ss << std::fixed << std::setprecision(1) << (double)r;
    label_r += ss.str();
    
    ss.str("");
    TString label_k = "#kappa = ";
    ss << std::fixed << std::setprecision(1) << (double)k;
    label_k += ss.str();
    
    ss.str("");
    TString label_d = "#delta = ";
    ss << std::fixed << std::setprecision(0) << (double)delta;
    label_d += ss.str();
    label_d += "#circ";
    
    ss.str("");
    TString label_g = "#gamma = ";
    ss << std::fixed << std::setprecision(0) << (double)gamma;
    label_g += ss.str();
    label_g += "#circ";
    
    
    
    TLegend legA(0.6,0.6,0.925,0.9,"");
    legA.SetLineStyle(0);
    legA.SetLineColor(0);
    legA.SetFillColor(0);
    legA.SetTextFont(22);
    legA.SetTextColor(1);
    legA.SetTextSize(0.06);
    legA.SetTextAlign(12);
    
    TLegend legAd(0.6,0.6,0.925,0.9,"");
    legAd.SetLineStyle(0);
    legAd.SetLineColor(0);
    legAd.SetFillColor(0);
    legAd.SetTextFont(22);
    legAd.SetTextColor(1);
    legAd.SetTextSize(0.06);
    legAd.SetTextAlign(12);

    A_mix->SetFillColor(kBlue);    
    A_mix->SetLineColor(kBlue);
    A_mix->SetLineWidth(5);
    A_mix->Draw("");
    
    A_mix_2->SetLineColor(kRed);
    A_mix_2->SetLineWidth(5);
    A_mix_2->SetLineStyle(kDashed);
    A_mix_2->Draw("same");

    A_mix_3->SetLineColor(kGreen+3);
    A_mix_3->SetLineWidth(5);
    A_mix_3->SetLineStyle(kDashed);
    A_mix_3->Draw("same");
    
    TLegendEntry* leA = legA.AddEntry(A_mix,"A_{P} = 0 %","");
    leA->SetTextColor(kBlue);
    leA = legA.AddEntry(A_mix_2,"A_{P} = 5 %","");
    leA->SetTextColor(kRed);
    leA = legA.AddEntry(A_mix_3,"A_{P} = 10 %","");
    leA->SetTextColor(kGreen+3);

    legA.Draw();
    graph->Draw("same");
    c->Print(((string)OutputDir+"a_mix_Ap.pdf").c_str());
    c->Print(((string)OutputDir+"a_mix_Ap.eps").c_str());

    A_mix->Draw("");
    A_mix_4->SetLineColor(kRed);
    A_mix_4->SetLineWidth(5);
    A_mix_4->SetLineStyle(kDashed);
    A_mix_4->Draw("same");
    
    A_mix_5->SetLineColor(kGreen+3);
    A_mix_5->SetLineWidth(5);
    A_mix_5->SetLineStyle(kDashed);
    A_mix_5->Draw("same");
    
    TLegendEntry* leAd = legAd.AddEntry(A_mix,"A_{D} = 0 %","");
    leAd->SetTextColor(kBlue);
    leAd = legAd.AddEntry(A_mix_2,"A_{D} = 5 %","");
    leAd->SetTextColor(kRed);
    leAd = legAd.AddEntry(A_mix_3,"A_{D} = 10 %","");
    leAd->SetTextColor(kGreen+3);
    
    legAd.Draw();
    graph->Draw("same");
    c->Print(((string)OutputDir+"a_mix_Ad.pdf").c_str());
    c->Print(((string)OutputDir+"a_mix_Ad.eps").c_str());
    

    A_q->SetLineColor(kBlue);
    A_q->SetLineWidth(5);
    A_q->Draw("C");
    
    A_q_2->SetLineColor(kRed);
    A_q_2->SetLineWidth(5);
    //A_q_2->SetLineStyle(kDashed);
    A_q_2->Draw("same");

    A_q_3->SetLineColor(kGreen+3);
    A_q_3->SetLineWidth(3);
    A_q_3->SetLineStyle(kDashed);
    A_q_3->Draw("same");

    graph->Draw("same");
    c->Print(((string)OutputDir+"a_q_Ap.pdf").c_str());
    c->Print(((string)OutputDir+"a_q_Ap.eps").c_str());
    
    A_q->Draw("C");
    A_q_4->SetLineColor(kRed);
    A_q_4->SetLineWidth(5);
    A_q_4->SetLineStyle(kDashed);
    A_q_4->Draw("same");
    
    A_q_5->SetLineColor(kGreen+3);
    A_q_5->SetLineWidth(3);
    A_q_5->SetLineStyle(kDashed);
    A_q_5->Draw("same");
    
    graph->Draw("same");
    c->Print(((string)OutputDir+"a_q_Ad.pdf").c_str());
    c->Print(((string)OutputDir+"a_q_Ad.eps").c_str());
    
    

    A_mix_true->SetLineColor(kBlue+1);
    A_mix_true->SetLineWidth(5);
    A_mix_true->Draw("");
    graph->Draw("same");
    c->Print(((string)OutputDir+"a_f.pdf").c_str());
    c->Print(((string)OutputDir+"a_f.eps").c_str());
    
    A_q_true->SetLineColor(kBlue+1);
    A_q_true->SetLineWidth(5);
    A_q_true->Draw("C");
    graph->Draw("same");
    c->Print(((string)OutputDir+"a_i.pdf").c_str());
    c->Print(((string)OutputDir+"a_i.eps").c_str());
    
    //return;
    
    h_t_fit->SetLineColor(kBlue);
    h_t_fit->SetLineWidth(5);
    h_t_fit->SetMarkerColor(kBlue); 
    h_t_fit->Draw("histc");
    //     leg.Draw();
    
    c->Print(((string)OutputDir+"h_t.eps").c_str());
    c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t.pdf").c_str());
    gPad->SetLogy(1);
    c->Print(((string)OutputDir+"h_t_log.eps").c_str());
    c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
    gPad->SetLogy(0);
    
    if((string)channel=="norm"){
        h_t_mixed_fit->SetMarkerColor(kRed); 
        h_t_mixed_fit->SetLineColor(kRed);
        h_t_mixed_fit->DrawNormalized("histc",1);
        
        h_t_unmixed_fit->SetMarkerColor(kBlue); 
        h_t_unmixed_fit->SetLineColor(kBlue);
        h_t_unmixed_fit->DrawNormalized("histcsame",1);
        
        h_t_untagegged_fit->SetMarkerColor(kGreen); 
        h_t_untagegged_fit->SetLineColor(kGreen);
        //h_t_untagegged_fit->DrawNormalized("histcsame",1);
        
        c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
        c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed.pdf").c_str());
        
        TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);    
        h_asym_fit->SetMinimum(-1);
        h_asym_fit->SetMaximum(1);
        h_asym_fit->SetLineColor(kRed);
        h_asym_fit->Draw("histc");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());
        c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
    }
    
    else{        
        
        //         h_t_fit_p->SetAxisRange(0.4, max_TAU*2,"X");
        //         h_t_fit_p->Smooth(10,"R");
        //         h_t_fit_p->SetAxisRange(min_TAU, max_TAU,"X");
        
        h_t_fit_m->Scale(1./10.);
        h_t_fit_p->Scale(1./10.);
        h_t_fit_p->SetMinimum(0);
        h_t_fit_p->SetMaximum(4000);
        h_t_fit_p->SetLineWidth(4);
        h_t_fit_p->SetLineColor(kBlack);
        h_t_fit_p->Draw("e1");
        
        h_t_fit_mp->Scale(1./10.);
        //         h_t_fit_mp->SetBinContent(h_t_fit_pp->FindBin(0.),0.);
        //         h_t_fit_mp->Smooth(10);
        h_t_fit_mp->SetLineWidth(4);
        h_t_fit_mp->SetMarkerColor(kBlue); 
        h_t_fit_mp->SetLineColor(kBlue);
        h_t_fit_mp->Draw("histcsame");
        
        
        h_t_fit_pp->Scale(1./10.);
        //         h_t_fit_pp->SetBinContent(h_t_fit_pp->FindBin(0.),h_t_fit_p->GetBinContent(h_t_fit_pp->FindBin(0.)));
        //         h_t_fit_pp->Smooth(10);
        h_t_fit_pp->SetLineWidth(4);        
        h_t_fit_pp->SetMarkerColor(kRed); 
        h_t_fit_pp->SetLineColor(kRed);
        h_t_fit_pp->Draw("histcsame");
        
        //         leg.AddEntry(h_t_fit_p,"Total","l");
        TLegendEntry* le = leg.AddEntry(h_t_fit_pp,"B^{0}_{s}(t)#rightarrow f","l");
        le->SetTextColor(kRed);                    
        le = leg.AddEntry(h_t_fit_mp,"#bar{B^{0}_{s}}(t)#rightarrow f","l");
        le->SetTextColor(kBlue);        
        leg.AddEntry((TObject*)0,label_r,"");
        leg.AddEntry((TObject*)0,label_k,"");
        
        leg.AddEntry((TObject*)0,label_d,"");
        leg.AddEntry((TObject*)0,label_g,"");
        leg.Draw();
        
        c->Print(((string)OutputDir+"h_t_mixed_p.C").c_str());
        c->Print(((string)OutputDir+"h_t_mixed_p.root").c_str());
        c->Print(((string)OutputDir+"h_t_mixed_p.eps").c_str());
        
        
        //         h_t_fit_mp->SetMaximum(1600);
        h_t_fit_mp->Draw("histc");
        h_t_fit_pp->Draw("histcsame");
        
        h_t_fit_pm->Scale(1./10.);
        //         h_t_fit_pm->Smooth(10);
        h_t_fit_pm->SetLineWidth(4);
        h_t_fit_pm->SetLineStyle(kDashed);             
        h_t_fit_pm->SetMarkerColor(kBlue+1); 
        h_t_fit_pm->SetLineColor(kBlue+1);
        h_t_fit_pm->Draw("histcsame");
        
        h_t_fit_mm->Scale(1./10.);
        //         h_t_fit_mm->Smooth(10);
        h_t_fit_mm->SetLineWidth(4);
        h_t_fit_mm->SetLineStyle(kDashed);                     
        h_t_fit_mm->SetMarkerColor(kRed+1); 
        h_t_fit_mm->SetLineColor(kRed+1);
        h_t_fit_mm->Draw("histcsame");
        
        TLegend leg3(0.6,0.6,0.925,0.9,"");
        leg3.SetLineStyle(0);
        leg3.SetLineColor(0);
        leg3.SetFillColor(0);
        leg3.SetTextFont(22);
        leg3.SetTextColor(1);
        leg3.SetTextSize(0.06);
        leg3.SetTextAlign(12);
        leg3.SetNColumns(2);
        
        le = leg3.AddEntry(h_t_fit_pp,"B^{0}_{s}(t)#rightarrowf","l");
        le->SetTextColor(kRed);                    
        le = leg3.AddEntry(h_t_fit_mp,"#bar{B^{0}_{s}}(t)#rightarrowf","l");
        le->SetTextColor(kBlue);    
        
        le = leg3.AddEntry(h_t_fit_mm,"#bar{B^{0}_{s}}(t)#rightarrow#bar{f}","l");
        le->SetTextColor(kRed+1);                    
        le = leg3.AddEntry(h_t_fit_pm,"B^{0}_{s}(t)#rightarrow#bar{f}","l");
        le->SetTextColor(kBlue+1);    
        
        leg3.AddEntry((TObject*)0,label_r,"");
        leg3.AddEntry((TObject*)0,label_k,"");
        
        leg3.AddEntry((TObject*)0,label_d,"");
        leg3.AddEntry((TObject*)0,label_g,"");
        
        leg3.Draw();
        
        c->Print(((string)OutputDir+"h_t_mixed.C").c_str());
        c->Print(((string)OutputDir+"h_t_mixed.root").c_str());
        c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
        c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());
        
        
        
        TH1D* h_asymCP_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);    
        h_asymCP_fit->SetLineColor(kBlack);
        h_asymCP_fit->SetLineWidth(5);
        //         h_asymCP_fit->SetMinimum(-1.);
        //         h_asymCP_fit->SetMaximum(1.);
        h_asymCP_fit->SetYTitle("A_{CP}");
        h_asymCP_fit->Draw("histc");
        
        TH1D* h_asymCP2_fit = (TH1D*) h_N_mixed_m_fit->GetAsymmetry(h_N_mixed_p_fit);    
        h_asymCP2_fit->SetLineColor(kRed);
        h_asymCP2_fit->SetLineWidth(5);
        //         h_asymCP_fit->SetMinimum(-1.);
        //         h_asymCP_fit->SetMaximum(1.);
        h_asymCP2_fit->Draw("histcsame");
        graph->Draw("same");
        
        c->Print(((string)OutputDir+"h_asymCP.eps").c_str());    
        c->Print(((string)OutputDir+"h_asymCP.C").c_str());    
        c->Print(((string)OutputDir+"h_asymCP.root").c_str());    
        
        
        TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);    
        h_asym_p_fit->SetLineColor(kBlack);
        h_asym_p_fit->SetLineWidth(5);
        h_asym_p_fit->SetMinimum(-1.);
        h_asym_p_fit->SetMaximum(1.);
        h_asym_p_fit->Draw("histc");
        
        TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);    
        h_asym_m_fit->SetLineWidth(5);
        h_asym_m_fit->SetLineColor(kGreen+3);
        h_asym_m_fit->SetLineStyle(kDashed);
        h_asym_m_fit->Draw("histcsame");
        
        TLegend leg2(0.45,0.75,0.6,0.92,"");
        leg2.SetLineStyle(0);
        leg2.SetLineColor(0);
        leg2.SetFillColor(0);
        leg2.SetTextFont(22);
        leg2.SetTextColor(1);
        leg2.SetTextSize(0.06);
        leg2.SetTextAlign(12);

//         le = leg2.AddEntry(h_asym_p_fit,"f=D_{s}^{-}#pi^{+}#pi^{+}#pi^{-}","l");        
        le = leg2.AddEntry(h_asym_p_fit,"f=D_{s}^{-}K^{+}#pi^{+}#pi^{-}","l");
        le->SetTextColor(kBlack);                    
//         le = leg2.AddEntry(h_asym_m_fit,"#bar{f}=D_{s}^{+}#pi^{-}#pi^{-}#pi^{+}","l");
        le = leg2.AddEntry(h_asym_m_fit,"#bar{f}=D_{s}^{+}K^{-}#pi^{-}#pi^{+}","l");
        le->SetTextColor(kGreen+3);                    
        leg2.Draw();
        
        graph->Draw("same");
        c->Print(((string)OutputDir+"h_asym.eps").c_str());    
        c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());
        c->Print(((string)OutputDir+"h_asym.C").c_str());    
        c->Print(((string)OutputDir+"h_asym.root").c_str());    
        
        h_asym_p_fit->Add(h_asym_m_fit,-1.);    
        h_asym_p_fit->Smooth(2);
        h_asym_p_fit->SetLineColor(kMagenta+3);
        h_asym_p_fit->SetLineWidth(5);
        h_asym_p_fit->SetMinimum(-0.4);
        h_asym_p_fit->SetMaximum(0.4);
        h_asym_p_fit->SetYTitle("A_{CP}");
        //         h_asym_p_fit->SetYTitle("#DeltaA_{flav}");
        h_asym_p_fit->Draw("histc");
        graph->Draw("same");
        c->Print(((string)OutputDir+"h_delta_asym.eps").c_str());
        
        h_t_fit_p->Draw("e1");
        h_t_fit_m->Draw("histcsame");
        c->Print(((string)OutputDir+"h_t_f.eps").c_str());
        
        //         h_t_fit_p->Rebin(20);
        //         h_t_fit_m->Rebin(20);
        
        TH1D* h_asym2_fit = (TH1D*) h_N_m_fit->GetAsymmetry(h_N_p_fit);    
        h_asym2_fit->SetLineColor(kBlack);
        h_asym2_fit->SetLineWidth(5);
        /*        h_asym2_fit->Rebin(40);*/
        //         h_asym2_fit->SetMinimum(-1.);
        //         h_asym2_fit->SetMaximum(1.);
        
        //         h_asym2_fit->GetXaxis()->SetRangeUser(0.,10.);
        //         h_asym2_fit->Draw("e1");
        
        TH1D* h_N_fake = new TH1D("h_N_p_fit",";t [ps];Events  ",100,0.,10.);
        h_N_fake->SetMaximum(h_asym2_fit->GetMaximum());
        h_N_fake->SetMinimum(h_asym2_fit->GetMinimum());
        h_N_fake->Draw("e");
        h_asym2_fit->SetLineColor(kRed);
        h_asym2_fit->Draw("histcsame");
        //         h_asym2_fit->Rebin(2);
        //         h_asym2_fit->SetLineColor(kBlue);
        //         h_asym2_fit->Draw("histcsame");
        graph->Draw("same");
        c->Print(((string)OutputDir+"h_asym2.eps").c_str());
        c->Print(((string)OutputDir+"h_asym2.C").c_str());
        
        h_N_unmixed_p_fit->Add(h_N_mixed_m_fit,1.);    
        h_N_unmixed_m_fit->Add(h_N_mixed_p_fit,1.);    
        TH1D* h_asym3_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_unmixed_m_fit);    
        h_asym3_fit->Draw("histc");
        graph->Draw("same");
        c->Print(((string)OutputDir+"h_asym3.eps").c_str());
        
        
    }
    
    
    
}


void calib(){

	const int nBins = 4;
 	double x[nBins]; 
        double xerr[nBins]; 
        double xerrL[nBins]; 
        double xerrH[nBins]; 
        double y[nBins]; 
        double yerr[nBins]; 


	x[0] = 0.225;
	xerr[0] = 0.;
	xerrL[0]= 0.225- 0. ;
	xerrH[0]= 0.3 - 0.225 ;

	x[1] = 0.353;
	xerr[1] = 0.;
	xerrL[1]= 0.353- 0.3 ;
	xerrH[1]= 0.4 - 0.353 ;

	x[2] = 0.444;
	xerr[2] = 0.;
	xerrL[2]= 0.444- 0.4 ;
	xerrH[2]= 0.48 - 0.444 ;

	x[3] = 0.489;
	xerr[3] = 0.;
	xerrL[3]= 0.489- 0.48 ;
	xerrH[3]= 0.5 - 0.489 ;

        y[0] = 0.298069;
	yerr[0] = 0.013;

        y[1] = 0.386296;
	yerr[1] = 0.0125;

        y[2] = 0.46;
	yerr[2] = 0.01;

        y[3] = 0.49;
	yerr[3] = 0.013;

        TGraphErrors *ResoRelation_g = new TGraphErrors(nBins, x,y,xerr,yerr);
	TGraphAsymmErrors *ResoRelation_ga = new TGraphAsymmErrors(nBins, x,y,xerrL,xerrH,yerr,yerr);

	ResoRelation_ga->SetTitle(";#eta_{SS};#omega");
	ResoRelation_ga->SetMinimum(0);
 	ResoRelation_ga->SetMaximum(0.55);

	//define polynom for fit
	TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*(x-0.40398) ", 0., 0.5);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetParNames("p0","p1");
	fitFunc->SetParameters(0.428425,0.787313);
	fitFunc->SetParLimits(0,0.,1);
	fitFunc->SetParLimits(1,0.,1.5);

	TFitResultPtr result = ResoRelation_g->Fit(fitFunc,"RS");
	TCanvas* c = new TCanvas();
	ResoRelation_ga->Draw("AP");
	fitFunc->Draw("same");
	c->Print("norm_calib/tagCalib.eps");
}


void compareCP(){

    RooArgList x;
    RooArgList m;
    
    RooRealVar* r = new  RooRealVar("r", "r", 0.4, 0.,1.);
    RooRealVar* delta = new  RooRealVar("delta", "delta", 0, -90.,90.);
    RooRealVar* gamma = new  RooRealVar("gamma", "gamma", 0, 0.,180.);
    RooRealVar *k = new  RooRealVar("k", "k", 0, 0.,1.);

    x.add(*r); 
    x.add(*delta);
    x.add(*gamma);
    x.add(*k);

    m.add(RooRealConstant::value(0.4));
    m.add(RooRealConstant::value(-10));
    m.add(RooRealConstant::value(70));
    m.add(RooRealConstant::value(0.5));
    
    
    vector<double> sigma;
    //sigma.push_back(0.03);
    //sigma.push_back(15);
    //sigma.push_back(16.4);
    //sigma.push_back(0.05);
    
    sigma.push_back(0.033);
    sigma.push_back(12);
    sigma.push_back(12);
    sigma.push_back(0.054);
    
    
    double rho_r_delta = -0.03775 ;
    double rho_r_gamma = -0.1892;
    double rho_r_kappa = 0.597;
    
    double rho_delta_gamma = 0.1602;    
    double rho_delta_kappa = -0.1277;
    
    double rho_gamma_kappa = -0.04696;
    
    TMatrixDSym cov(4);
    cov(0,0) = 1.;
    cov(0,1) = cov(1,0) = rho_r_delta;
    cov(0,2) = cov(2,0) = rho_r_gamma;
    cov(0,3) = cov(3,0) = rho_r_kappa;
    
    cov(1,1) = 1.;
    cov(1,2) = cov(2,1) = rho_delta_gamma;
    cov(1,3) = cov(3,1) = rho_delta_kappa;
    
    cov(2,2) = 1.;
    cov(2,3) = cov(3,2) = rho_gamma_kappa;
    
    cov(3,3) = 1.;
    
    //
    //1     -0.03775     -0.1892       0.597 
    // -0.03775           1      0.1602     -0.1277 
    //-0.1892      0.1602           1    -0.04696 
    //0.597     -0.1277    -0.04696           1 
    
    //1.000 -0.109  0.052
    //-0.109  1.000  0.092
    //0.052  0.092  1.000
    
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov(i,j) = cov(i,j) * sigma[i] * sigma[j];
        }
    }
    
    cov.Print();

    int N_toys = 10000;
    RooMultiVarGaussian* gauss_cov = new RooMultiVarGaussian("gauss_cov","gauss_cov",x, m, cov);
    RooDataSet* data = gauss_cov->generate(x,N_toys);
    
    vector<double> C,D,Dbar,S,Sbar;
    
    for(Int_t i=0;i<N_toys;i++) 
    {
        RooArgList* l = (RooArgList*) data->get(i);
    
        double r_val = ((RooRealVar *) l->at(0))->getVal();
        double delta_val = ((RooRealVar *) l->at(1))->getVal()/360.*2.*pi;
        double gamma_val = ((RooRealVar *) l->at(2))->getVal()/360.*2.*pi;
        double k_val = ((RooRealVar *) l->at(3))->getVal();
        
        C.push_back( (1.-pow(r_val,2))/(1.+pow(r_val,2)) );

        D.push_back( (-2.*r_val*k_val*cos(delta_val-gamma_val)) /(1.+pow(r_val,2)) );
        Dbar.push_back( (-2.*r_val*k_val*cos(delta_val+gamma_val)) /(1.+pow(r_val,2)) );
       
        S.push_back( (2.*r_val*k_val*sin(delta_val-gamma_val)) /(1.+pow(r_val,2)) );
        Sbar.push_back( (-2.*r_val*k_val*sin(delta_val+gamma_val)) /(1.+pow(r_val,2)) );
    }
    
    double C_var = 0;
    double D_var = 0;
    double Dbar_var = 0;
    double S_var = 0;
    double Sbar_var = 0;
    
    for(int j=0; j < N_toys; j++){
        for(int k=j+1; k < N_toys; k++){
            C_var += pow(C[j] - C[k],2);
            D_var += pow(D[j] - D[k],2);
            Dbar_var += pow(Dbar[j] - Dbar[k],2);
            S_var += pow(S[j] - S[k],2);
            Sbar_var += pow(Sbar[j] - Sbar[k],2);
        }
    }
    double C_err = sqrt(C_var/((double)N_toys*((double)N_toys-.1)));
    double D_err = sqrt(D_var/((double)N_toys*((double)N_toys-.1)));
    double Dbar_err = sqrt(Dbar_var/((double)N_toys*((double)N_toys-.1)));
    double S_err = sqrt(S_var/((double)N_toys*((double)N_toys-.1)));
    double Sbar_err = sqrt(Sbar_var/((double)N_toys*((double)N_toys-.1)));

    
    double r_val = ((RooRealVar *) m.at(0))->getVal();
    double delta_val = ((RooRealVar *) m.at(1))->getVal()/360.*2.*pi;
    double gamma_val = ((RooRealVar *) m.at(2))->getVal()/360.*2.*pi;
    double k_val = ((RooRealVar *) m.at(3))->getVal();
    
    double C_val =( (1.-pow(r_val,2))/(1.+pow(r_val,2)) );
    double D_val =( (-2.*r_val*k_val*cos(delta_val-gamma_val)) /(1.+pow(r_val,2)) );
    double Dbar_val =( (-2.*r_val*k_val*cos(delta_val+gamma_val)) /(1.+pow(r_val,2)) );
    double S_val =( (2.*r_val*k_val*sin(delta_val-gamma_val)) /(1.+pow(r_val,2)) );
    double Sbar_val =( (-2.*r_val*k_val*sin(delta_val+gamma_val)) /(1.+pow(r_val,2)) );

    cout << "C = " << C_val << " +- " << C_err << endl;
    cout << "D = " << D_val << " +- " << D_err << endl;
    cout << "Dbar = " << Dbar_val << " +- " << Dbar_err << endl;
    cout << "S = " << S_val << " +- " << S_err << endl;
    cout << "Sbar = " << Sbar_val << " +- " <<  Sbar_err << endl;
}    

void calculateAverageReso(){

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);

    int run,year,Ds_finalState,trigger;
    double t,dt;
    double sw;
    double K_P;

    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm_tagged.root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("run",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("K_minus_fromDs_P",1);
//     tree_norm->SetBranchStatus("K_plus_P",1);
    tree_norm->SetBranchStatus("*finalState*",1);

    tree_norm->SetBranchAddress("Bs_BsDTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("K_minus_fromDs_P",&K_P);
//     tree_norm->SetBranchAddress("K_plus_P",&K_P);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);

    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);    
    FitParameter  offset_sigma_dt_Run1("offset_sigma_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run1("scale_sigma_dt_Run1",1,1.2,0.1);
    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.2,0.1);
    FitParameter  offset_sigma_dt_Run2_17("offset_sigma_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2_17("scale_sigma_dt_Run2_17",1,1.2,0.1);

    int bins = 60;
    TH1D* h_dt_norm = new TH1D("h_dt_norm_comb",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1 = new TH1D("h_dt_norm_Run1",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2 = new TH1D("h_dt_norm_Run2",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17 = new TH1D("h_dt_norm_Run2_17",";#sigma_{t} [ps];Events (a.u.) ",bins,min_TAUERR,max_TAUERR);

    TFile* f_det_asym_Run1 = new TFile("../Asymmetries/AsymmetryHistos/det_asym_Run1.root");
    TH1D* det_asym_Run1 = (TH1D*) f_det_asym_Run1->Get("det_asym_Run1");
    TH1D* det_asym_Run1_toy = (TH1D*) f_det_asym_Run1->Get("det_asym_Run1");
    TFile* f_det_asym_Run2 = new TFile("../Asymmetries/AsymmetryHistos/det_asym_Run2.root");
    TH1D* det_asym_Run2 = (TH1D*) f_det_asym_Run2->Get("det_asym_Run2");
    TH1D* det_asym_Run2_toy = (TH1D*) f_det_asym_Run2->Get("det_asym_Run2");

    double sigma_dt = 0.;
    double sigma_dt_Run1 = 0.;
    double sigma_dt_Run2 = 0.;
    double sigma_dt_Run2_17 = 0.;
    double N_Run1 = 0.;
    double N_Run2 = 0.;
    double N_Run2_17 = 0.;

    double det_asy_Run1= 0.;
    double det_asy_Run2= 0.;
    double det_asy_N_Run1= 0.;
    double det_asy_N_Run2= 0.;

    double n_4 = 0.;

    for(int i=0; i< tree_norm->GetEntries(); i++)
    {    
        tree_norm->GetEntry(i);
	if(t < min_TAU || t > max_TAU )continue;
	if( dt < min_TAUERR || dt > max_TAUERR )continue;

	if(run==1){
		sigma_dt += sw*(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt);
		sigma_dt_Run1 += sw*(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt);
		N_Run1 += sw;
		h_dt_norm_Run1->Fill(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt,sw);
	}
	else if(year < 17){
		sigma_dt += sw*(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt);
		sigma_dt_Run2 += sw*(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt);
		N_Run2 += sw;
		h_dt_norm_Run2->Fill(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt,sw);
	}
	else{
		sigma_dt += sw*(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt);
		sigma_dt_Run2_17 += sw*(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt);
		N_Run2_17 += sw;
		h_dt_norm_Run2_17->Fill(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt,sw);
	}

	if(Ds_finalState == 4)n_4 += sw;

	double det_asy= 0.;
	if(run==2){
			double kaon_p = min(69.99,max(3.01,K_P/1000.));
			det_asy = -0.0091;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy_Run2 += det_asy * sw;
			det_asy_N_Run2 += sw;
	}
	else{
			double kaon_p = min(99.99,max(5.01,K_P/1000.));
			det_asy = -0.01;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			det_asy_Run1 += det_asy * sw;
			det_asy_N_Run1 += sw;
	}
	

    }

    cout << "n4 = " << n_4/(N_Run1+N_Run2+N_Run2_17) <<  endl;

    cout << "det_asy_Run1 = " << det_asy_Run1/det_asy_N_Run1 <<  endl;
    cout << "det_asy_Run2 = " << det_asy_Run2/det_asy_N_Run2 <<  endl;

    cout << "Reso Run1 = " << sigma_dt_Run1/N_Run1 *1000.<< endl;
    cout << "Reso Run2 = " << sigma_dt_Run2/N_Run2 *1000.<< endl;
    cout << "Reso Run2_17 = " << sigma_dt_Run2_17/N_Run2_17 *1000.<< endl;
    cout << "Reso tot = " << sigma_dt/(N_Run1+N_Run2+N_Run2_17) *1000.<< endl;

return;

    TCanvas* c = new TCanvas();
    h_dt_norm_Run1->DrawNormalized("e1",1);
    h_dt_norm_Run2->SetMarkerColor(kRed);
    h_dt_norm_Run2->SetLineColor(kRed);
    h_dt_norm_Run2_17->SetMarkerColor(kBlue);
    h_dt_norm_Run2_17->SetLineColor(kBlue);
    h_dt_norm_Run2->DrawNormalized("histcsame",1);
    h_dt_norm_Run2_17->DrawNormalized("histcsame",1);
    c->Print("dt_calib.eps");

// 	return;

    int N_toys = 100;
    double rho = -0.9 ;    
    TMatrixDSym cov(2);
    cov(0,0) = 1.;
    cov(0,1) = cov(1,0) = rho;    
    cov(1,1) = 1.;    

    RooArgList x_Run1;
    RooArgList m_Run1;
    RooRealVar* off_Run1 = new  RooRealVar("off_Run1", "off_Run1", offset_sigma_dt_Run1);
    RooRealVar* scale_Run1 = new  RooRealVar("scale_Run1", "scale_Run1", scale_sigma_dt_Run1);
    x_Run1.add(*off_Run1); 
    x_Run1.add(*scale_Run1);
    m_Run1.add(RooRealConstant::value(offset_sigma_dt_Run1));
    m_Run1.add(RooRealConstant::value(scale_sigma_dt_Run1));
        
    vector<double> sigma_Run1;    
    sigma_Run1.push_back(1.5/1000.);
    sigma_Run1.push_back(0.042);    
    
    TMatrixDSym cov_Run1(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run1(i,j) = cov(i,j) * sigma_Run1[i] * sigma_Run1[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run1 = new RooMultiVarGaussian("gauss_cov_Run1","gauss_cov_Run1",x_Run1, m_Run1, cov_Run1);
    RooDataSet* data_Run1 = gauss_cov_Run1->generate(x_Run1,N_toys);


    RooArgList x_Run2;
    RooArgList m_Run2;
    RooRealVar* off_Run2 = new  RooRealVar("off_Run2", "off_Run2", offset_sigma_dt_Run2);
    RooRealVar* scale_Run2 = new  RooRealVar("scale_Run2", "scale_Run2", scale_sigma_dt_Run2);
    x_Run2.add(*off_Run2); 
    x_Run2.add(*scale_Run2);
    m_Run2.add(RooRealConstant::value(offset_sigma_dt_Run2));
    m_Run2.add(RooRealConstant::value(scale_sigma_dt_Run2));
        
    vector<double> sigma_Run2;    
    sigma_Run2.push_back(1.5/1000.);
    sigma_Run2.push_back(0.042);    
    
    TMatrixDSym cov_Run2(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run2(i,j) = cov(i,j) * sigma_Run2[i] * sigma_Run2[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run2 = new RooMultiVarGaussian("gauss_cov_Run2","gauss_cov_Run2",x_Run2, m_Run2, cov_Run2);
    RooDataSet* data_Run2 = gauss_cov_Run2->generate(x_Run2,N_toys);


    RooArgList x_Run2_17;
    RooArgList m_Run2_17;
    RooRealVar* off_Run2_17 = new  RooRealVar("off_Run2_17", "off_Run2_17", offset_sigma_dt_Run2_17);
    RooRealVar* scale_Run2_17 = new  RooRealVar("scale_Run2_17", "scale_Run2_17", scale_sigma_dt_Run2_17);
    x_Run2_17.add(*off_Run2_17); 
    x_Run2_17.add(*scale_Run2_17);
    m_Run2_17.add(RooRealConstant::value(offset_sigma_dt_Run2_17));
    m_Run2_17.add(RooRealConstant::value(scale_sigma_dt_Run2_17));
        
    vector<double> sigma_Run2_17;    
    sigma_Run2_17.push_back(1.5/1000.);
    sigma_Run2_17.push_back(0.042);    
    
    TMatrixDSym cov_Run2_17(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run2_17(i,j) = cov(i,j) * sigma_Run2_17[i] * sigma_Run2_17[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run2_17 = new RooMultiVarGaussian("gauss_cov_Run2_17","gauss_cov_Run2_17",x_Run2_17, m_Run2_17, cov_Run2_17);
    RooDataSet* data_Run2_17 = gauss_cov_Run2_17->generate(x_Run2_17,N_toys);

    
    vector<double> sigma_dt_vec,sigma_dt_Run1_vec,sigma_dt_Run2_vec,sigma_dt_Run2_17_vec;
    vector<double> det_asy_Run1_vec,det_asy_Run2_vec;
    sigma_dt_vec.push_back(sigma_dt/(N_Run1+N_Run2+N_Run2_17));       
    sigma_dt_Run1_vec.push_back(sigma_dt_Run1/N_Run1);       
    sigma_dt_Run2_vec.push_back(sigma_dt_Run2/N_Run2);       
    sigma_dt_Run2_17_vec.push_back(sigma_dt_Run2_17/N_Run2_17);   
    det_asy_Run1_vec.push_back(det_asy_Run1/det_asy_N_Run1);
    det_asy_Run2_vec.push_back(det_asy_Run2/det_asy_N_Run2);  

    for(Int_t j=0;j<N_toys;j++) 
    {
        RooArgList* l_Run1 = (RooArgList*) data_Run1->get(j);    
        double off_Run1_val = ((RooRealVar *) l_Run1->at(0))->getVal();
        double scale_Run1_val = ((RooRealVar *) l_Run1->at(1))->getVal();

        RooArgList* l_Run2 = (RooArgList*) data_Run2->get(j);    
        double off_Run2_val = ((RooRealVar *) l_Run2->at(0))->getVal();
        double scale_Run2_val = ((RooRealVar *) l_Run2->at(1))->getVal();

        RooArgList* l_Run2_17 = (RooArgList*) data_Run2_17->get(j);    
        double off_Run2_17_val = ((RooRealVar *) l_Run2_17->at(0))->getVal();
        double scale_Run2_17_val = ((RooRealVar *) l_Run2_17->at(1))->getVal();

	sigma_dt = 0.;
	sigma_dt_Run1 = 0.;
	sigma_dt_Run2 = 0.;
	sigma_dt_Run2_17 = 0.;

	det_asy_Run1= 0.;
    	det_asy_Run2= 0.;
    	det_asy_N_Run1= 0.;
    	det_asy_N_Run2= 0.;

	for(int n = 1; n <= det_asym_Run1->GetNbinsX(); n++){
		det_asym_Run1_toy->SetBinContent(n, det_asym_Run1->GetBinContent(n) + gRandom->Gaus(0.,det_asym_Run1->GetBinError(n)) );
	}

	for(int n = 1; n <= det_asym_Run2->GetNbinsX(); n++){
		det_asym_Run2_toy->SetBinContent(n, det_asym_Run2->GetBinContent(n) + gRandom->Gaus(0.,det_asym_Run2->GetBinError(n)) );
	}


	
	for(int i=0; i< tree_norm->GetEntries(); i++)
	{    
		tree_norm->GetEntry(i);
		if(t < min_TAU || t > max_TAU )continue;
		if( dt < min_TAUERR || dt > max_TAUERR )continue;
	
		if(run==1){
			sigma_dt += sw*(off_Run1_val + scale_Run1_val * dt);
			sigma_dt_Run1 += sw*(off_Run1_val + scale_Run1_val * dt);
		}
		else if(year < 17){
			sigma_dt += sw*(off_Run2_val + scale_Run2_val * dt);
			sigma_dt_Run2 += sw*(off_Run2_val + scale_Run2_val * dt);
		}
		else{
			sigma_dt += sw*(off_Run2_17_val + scale_Run2_17_val * dt);
			sigma_dt_Run2_17 += sw*(off_Run2_17_val + scale_Run2_17_val * dt);
		}

		double det_asy= 0.;
		if(run==2){
			double kaon_p = min(69.99,max(3.01,K_P/1000.));
			det_asy = -0.0091;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			if(Ds_finalState == 4)det_asy_Run2 += det_asy * sw;
			if(Ds_finalState == 4)det_asy_N_Run2 += sw;
		}
		else{
			double kaon_p = min(99.99,max(5.01,K_P/1000.));
			det_asy = -0.01;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			if(Ds_finalState == 4)det_asy_Run1 += det_asy * sw;
			if(Ds_finalState == 4)det_asy_N_Run1 += sw;
		}

    	}
	sigma_dt_vec.push_back(sigma_dt/(N_Run1+N_Run2+N_Run2_17));       
	sigma_dt_Run1_vec.push_back(sigma_dt_Run1/N_Run1);       
	sigma_dt_Run2_vec.push_back(sigma_dt_Run2/N_Run2);       
	sigma_dt_Run2_17_vec.push_back(sigma_dt_Run2_17/N_Run2_17);   
	det_asy_Run1_vec.push_back(det_asy_Run1/det_asy_N_Run1);
	det_asy_Run2_vec.push_back(det_asy_Run2/det_asy_N_Run2);     
    }
    
    double sigma_dt_var = 0;
    double sigma_dt_var_Run1 = 0;
    double sigma_dt_var_Run2 = 0;
    double sigma_dt_var_Run2_17 = 0;

    double sigma_det_asy_Run1 = 0;
    double sigma_det_asy_Run2 = 0;

    for(int j=0; j < N_toys; j++){
        for(int k=j+1; k < N_toys; k++){
            sigma_dt_var += pow(sigma_dt_vec[j] - sigma_dt_vec[k],2);
            sigma_dt_var_Run1 += pow(sigma_dt_Run1_vec[j] - sigma_dt_Run1_vec[k],2);
            sigma_dt_var_Run2 += pow(sigma_dt_Run2_vec[j] - sigma_dt_Run2_vec[k],2);
            sigma_dt_var_Run2_17 += pow(sigma_dt_Run2_17_vec[j] - sigma_dt_Run2_17_vec[k],2);

            sigma_det_asy_Run1 += pow(det_asy_Run1_vec[j] - det_asy_Run1_vec[k],2);
            sigma_det_asy_Run2 += pow(det_asy_Run2_vec[j] - det_asy_Run2_vec[k],2);
        }
    }
    cout << "Reso Run1 = " << sigma_dt_Run1/N_Run1 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run1/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso Run2 = " << sigma_dt_Run2/N_Run2 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run2/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso Run2_17 = " << sigma_dt_Run2_17/N_Run2_17 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run2_17/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso tot = " << sigma_dt/(N_Run1+N_Run2+N_Run2_17) *1000.<< " +- " <<  sqrt(sigma_dt_var/((double)N_toys*((double)N_toys-.1))) *1000 << endl;

    cout << "DA Run1 = " <<  det_asy_Run1_vec[0] *100.<< " +- " <<  sqrt(sigma_det_asy_Run1/((double)N_toys*((double)N_toys-.1))) *100 << endl;
    cout << "DA Run2 = " <<  det_asy_Run2_vec[0] *100.<< " +- " <<  sqrt(sigma_det_asy_Run2/((double)N_toys*((double)N_toys-.1))) *100 << endl;
}

void calculateSweightsForToys(int step = 0){

	/// Load file
        NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
	TString inFileName = ((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str();
	TString outFileName = ((string)OutputDir+"sw_toys_"+anythingToString((int)step)+".root").c_str();

	TFile* f = new TFile(inFileName,"Open");
	TTree* tree = (TTree*) f->Get("DecayTree");
	tree->SetBranchStatus("N_Bs_sw",0);

	TFile* out = new TFile(outFileName,"RECREATE");
	TTree* out_tree = tree->CopyTree("");

 	double sw;
     	TBranch* b_w = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");

	TString channelString = "m(D_{s}K#pi#pi)" ;
	RooRealVar DTF_Bs_M("Bs_DTF_MM", channelString, 5200, 5700,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree));
	
	/// Signal pdf
	RooRealVar mean("mean", "mean", 5367,5350.,5390.); 
	RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
	RooGaussian* signal = new RooGaussian("signal","signal",DTF_Bs_M, mean,sigma);
	mean.setConstant();
	sigma.setConstant();

	/// Bkg pdf
	RooRealVar c0("c0", "c0", .0,-10,10); 
	RooRealVar c1("c1", "c1", .0,-10,10); 
	RooRealVar c2("c2", "c2", .0,-10,10); 
	RooChebychev* bkg= new RooChebychev("bkg","bkg",DTF_Bs_M, RooArgList(c0,c1));

	RooRealVar exp_par("exp_par","exp_par",-1.6508e-03,-10.,10.);
	RooExponential bkg_exp("bkg_exp","bkg_exp",DTF_Bs_M,exp_par);

        NamedParameter<double>  f_bkg("readBkgData::f_bkg", 0.);
	RooRealVar mean_bkg("mean_bkg", "mean_bkg", 5241.,5350.,5390.); 
	RooRealVar sigma_bkg("sigma_bkg", "sigma_bkg", 37.,0.,80.); 
	RooGaussian* partReco = new RooGaussian("partReco","partReco",DTF_Bs_M, mean_bkg,sigma_bkg);
// 	if(f_bkg==0.){
		mean_bkg.setConstant();
		sigma_bkg.setConstant();
// 	}
	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.75, 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.25, 0., data->numEntries());
	RooRealVar n_bkg2("n_bkg2", "n_bkg2", data->numEntries()*0.25, 0., data->numEntries());

	RooAddPdf* pdf;
	if(f_bkg==0.) pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, *bkg), RooArgList(n_sig, n_bkg));
	else pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, bkg_exp, *partReco), RooArgList(n_sig, n_bkg,n_bkg2));

	/// Fit
	RooFitResult* result = pdf->fitTo(*data,Save(kTRUE),NumCPU(6),Extended(kTRUE));
	result->Print();

	/// Calculate sweights
	SPlot *sPlot;
	if(f_bkg==0.)sPlot= new SPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
	else sPlot= new SPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg,n_bkg2));	

	for(int n = 0; n < tree->GetEntries(); n++){
			sw = sPlot->GetSWeight(n,"n_sig_sw");
			b_w->Fill();
	}
	
	/// Plotting
	TCanvas* c = new TCanvas();

	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
 	data->plotOn(frame,Name("data"),Binning(200));
	pdf->plotOn(frame,Name("pdf"),LineColor(kBlue));
	pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("signal"));
	pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("bkg_exp"));
	pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("bkg"));
	frame->Draw();
	c->Print("toyFit.eps");

	out_tree->Write();
	out->Close();
	f->Close();
}

int main(int argc, char** argv){

  NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->ProcessLine(".x ../lhcbStyle.C");
  gStyle->SetTitleOffset(1.,"Y");

//  calib();
//     test_multiGaussConstraints();
   produceMarginalPdfs();
  //for(int i = 0; i < 100; i++) fullTimeFit(atoi(argv[1])+i);
   fullTimeFit(atoi(argv[1]),(string)argv[2]);
   if((string)argv[2] == "gen" && addBkgToToys)calculateSweightsForToys(atoi(argv[1]));

  //for(int i = 1; i <= 100; i++)calculateSweightsForToys(i);

         // animate2(atoi(argv[1]));
// calculateAverageReso();

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
