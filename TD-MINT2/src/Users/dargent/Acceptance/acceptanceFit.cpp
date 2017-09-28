// Fits the time acceptance
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TPaveText.h>
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
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooNDKeysPdf.h"
#include "RooKeysPdf.h"
#include "RooBDecay.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/RooSplineProduct.h"
#include "Mint/Utils.h"
#include <fstream>
#include "Mint/HyperHistogram.h"
#include "Mint/HyperBinningPainter1D.h"
using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

/// HFLAV summer 17 values
double tau = 1.509;
double dgamma = 0.09; 
double deltaMs = 17.757;

/// MC values
double tau_MC = 1.512; // ?
double dgamma_MC = .109; // ?
double deltaMs_MC = 17.8; // or 20.0 ??? 

TH1D* createBinning(){

	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");		
	NamedParameter<double> TAU_min("TAU_min", 0.4);		
	NamedParameter<double> TAU_max("TAU_max", 10);		
        NamedParameter<int> minEventsPerBin("minEventsPerBin", 1000); 
	int dim = 1;

	TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
	TTree* tree= (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus(((string)Bs_TAU_Var).c_str(),1);
	double weight,t;
	tree->SetBranchAddress("weight",&weight);
	tree->SetBranchAddress(((string)Bs_TAU_Var).c_str(),&t);
	
	HyperPoint Min((double)TAU_min);
    	HyperPoint Max((double)TAU_max);
    	HyperCuboid limits(Min, Max );
	HyperPointSet points( dim );

	for (int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		
		HyperPoint point( dim );
		point.at(0)= t;
		point.addWeight(weight);
		points.push_back(point);
	}
	
   	/// Define binning based on MC
    	HyperHistogram hist(limits, points,                     
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI, 
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.001),
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (0)
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")
                         );
	/// Draw binning
// 	hist.setNames(HyperName(vars));
	hist.draw("Plot/binning");
	hist.drawDensity("Plot/density");

	TCanvas* c = new TCanvas();
        HyperBinningPainter1D painter(&hist);
 	TH1D* h = painter.getHistogram("binning");
	//h->Draw();
	//c->Print("test.eps");
	return h;
}


vector< vector<double> > fitSplineAcc(string CutString){

	// Options
	NamedParameter<string> BinningName("BinningName",(string)"default");
        NamedParameter<int> makePlots("makePlots", 0);
	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");
	NamedParameter<double> min_TAU("min_TAU", 0.4);
	NamedParameter<double> max_TAU("max_TAU", 10.);
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> useAdaptiveBinningKnots("useAdaptiveBinningKnots", 0);

	// Read Dataset
	TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
	TTree* tree = (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("*sw",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*finalState",1);
	
	//Define RooRealVar for observables
	RooRealVar Bs_TAU(((string)Bs_TAU_Var).c_str(), ((string)Bs_TAU_Var).c_str(), min_TAU, max_TAU, "ps");
	RooRealVar Bs_TAUERR(((string)Bs_TAU_Var+"ERR").c_str(), ((string)Bs_TAU_Var+"ERR").c_str(), 0.00001, 1.,"ps");
	RooRealVar N_Bs_sw("N_Bs_sw", "N_Bs_sw", 0.);
	RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
	RooRealVar year("year", "year", 0.);
	
	RooArgSet observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, N_Bs_sw);
	RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree), WeightVar(N_Bs_sw.GetName()), Cut(CutString.c_str()));
	
	///SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
	
	//SPLINE KNOTS
 	NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        NamedParameter<double> knot_values("knot_values", 3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.0, 1.7757e+00);
        vector<double> myBinning;    
        vector<double> values;

    	if(useAdaptiveBinningKnots){
		TH1D* binning = createBinning();	
		cout << endl << "knot positios: " << endl;
		for(int i = 1; i <= binning->GetNbinsX(); i++){
			cout << binning->GetBinCenter(i) << " , ";
			myBinning.push_back(binning->GetBinCenter(i));
			values.push_back(1.);
		}
		cout << endl << endl;
    	}
    	else { 
		myBinning = knot_positions.getVector();
    		values = knot_values.getVector() ;
    	}

	//SPLINE COEFFICIENTS
	RooArgList tacc_list;
        for(int i= 0; i< values.size(); i++){
		tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 10.0)));
	}

	tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), 1.0)));

	RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));

	tacc_list.add(*coeff_last);	

	//CUBIC SPLINE FUNCTION 
 	RooCubicSplineFun* spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_TAU, myBinning, tacc_list);
		
	RooRealVar trm_mean( "trm_mean" , "Gaussian resolution model mean", 0.0, "ps" );
	//RooRealVar trm_sigma( "trm_sigma" , "Gaussian resolution model sigma" , 0.040 , "ps");
	RooRealVar trm_scale( "trm_scale", "Gaussian resolution model scale factor", 1.20);
	RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *spl, trm_mean, Bs_TAUERR, trm_mean, trm_scale );
	
	RooBDecay* timePDF = new RooBDecay("Bdecay", "Bdecay", Bs_TAU, RooRealConstant::value(tau),
			RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
			RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs),
			trm, RooBDecay::SingleSided);
	
	
	///Fit and Print
	RooFitResult *myfitresult = timePDF->fitTo(*dataset, Save(1), Optimize(2), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE),NumCPU(numCPU));
	myfitresult->Print("v");
	
	//put coefficients into vector
	vector<double> myCoeffs,myCoeffsErr;
	for(int i= 0; i< values.size()+2; i++){
		myCoeffs.push_back(((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal());
        myCoeffsErr.push_back(((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError());
	}

	/// Plot
	int bin = 80;

	TCanvas* canvas = new TCanvas();
	canvas->cd();
	canvas->SetTopMargin(0.05);
	canvas->SetBottomMargin(0.05);

	TLegend leg(0.65,0.65,0.9,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.06);
	leg.SetTextAlign(12);

	RooPlot* frame_m = Bs_TAU.frame();	
	frame_m->GetXaxis()->SetLabelColor( kWhite);
	frame_m->GetYaxis()->SetTitleOffset(0.95);

	dataset->plotOn(frame_m, Binning(bin), Name("data"));
	timePDF->plotOn(frame_m, LineColor(kBlue+1),  Name("pdf"));
	spl->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline"));
	
	leg.AddEntry(frame_m->findObject("data"),"LHCb data","ep");
	leg.AddEntry(frame_m->findObject("pdf"),"Fit","l");
	leg.AddEntry(frame_m->findObject("spline"),"Acceptance","l");
	
	TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
	pad1->SetBorderMode(0);
	pad1->SetBorderSize(-1);
	pad1->SetBottomMargin(0.);
	pad1->Draw();
	pad1->cd();
	frame_m->GetYaxis()->SetRangeUser(0.011,frame_m->GetMaximum()*1.1);
	frame_m->Draw();
	leg.Draw();

	canvas->cd();
	TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
	pad2->SetBorderMode(0);
	pad2->SetBorderSize(-1);
	pad2->SetFillStyle(0);
	pad2->SetTopMargin(0.);
	pad2->SetBottomMargin(0.35);
	pad2->Draw();
	pad2->cd();
	
	RooPlot* frame_p = Bs_TAU.frame();
	frame_p->GetYaxis()->SetNdivisions(5);
	frame_p->GetYaxis()->SetLabelSize(0.12);
	frame_p->GetXaxis()->SetLabelSize(0.12);
	frame_p->GetXaxis()->SetTitleOffset(0.75);
	frame_p->GetXaxis()->SetTitleSize(0.2);
	frame_p->GetXaxis()->SetTitle("#font[132]{t(B_{s}) [ps]}");
	
	RooHist* pullHist  = frame_m->pullHist("data","pdf");
	frame_p->addPlotable(pullHist,"BX");
	
	double max = 5.0 ;
	double min = -5.0 ;

	double rangeX = max-min;
	double zero = max/rangeX;
	frame_p->GetYaxis()->SetRangeUser(min,max);

	TGraph* graph = new TGraph(2);
	graph->SetMaximum(max);
	graph->SetMinimum(min);
	graph->SetPoint(1,min_TAU,0);
	graph->SetPoint(2,max_TAU,0);
	
	TGraph* graph2 = new TGraph(2);
	graph2->SetMaximum(max);
	graph2->SetMinimum(min);
	graph2->SetPoint(1,min_TAU,-3);
	graph2->SetPoint(2,max_TAU,-3);
	graph2->SetLineColor(kRed);
	
	TGraph* graph3 = new TGraph(2);
	graph3->SetMaximum(max);
	graph3->SetMinimum(min);
	graph3->SetPoint(1,min_TAU,3);
	graph3->SetPoint(2,max_TAU,3);
	graph3->SetLineColor(kRed);
	
	frame_p->Draw();
	graph->Draw("same");
	graph2->Draw("same");
	graph3->Draw("same");
		
	pad2->Update();
	canvas->Update();
	canvas->SaveAs(("Plot/timeAccFit_"+(string)BinningName+ ".eps").c_str());
	
	pad1->SetLogy(1);
	pad1->Update();
	canvas->Update();
	canvas->SaveAs(("Plot/timeAccFit_"+(string)BinningName+ "_log.eps").c_str());
	
	// ???
	double chi2 = frame_m->chiSquare("pdf","data",values.size());
	cout << "chi2 = " << chi2 << endl;
	cout << "used datasets:     "<< CutString.c_str() << endl;
	
	file->Close();
    
   	vector< vector<double> > myCoeffsAndErr;
    	myCoeffsAndErr.push_back(myCoeffs);
    	myCoeffsAndErr.push_back(myCoeffsErr);    
    
	return myCoeffsAndErr;
}


void fitSplineAccRatio(string CutString){
    
    // Options
    NamedParameter<string> BinningName("BinningName",(string)"default");
    NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<int> numCPU("numCPU", 6);
    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 1);
    NamedParameter<int> useAdaptiveBinningKnots("useAdaptiveBinningKnots", 0);

    // Read Datasets
    TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
    TTree* tree = (TTree*) file->Get("DecayTree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("*TAU*",1);
    tree->SetBranchStatus("weight",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("*finalState",1);
    
    TFile* file_norm= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
    TTree* tree_norm = (TTree*) file_norm->Get("DecayTree");
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("weight",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*finalState",1);

    TFile* file_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
    TTree* tree_mc = (TTree*) file_mc->Get("DecayTree");
    tree_mc->SetBranchStatus("*",0);
    tree_mc->SetBranchStatus("*TAU*",1);
    tree_mc->SetBranchStatus("weight",1);
    tree_mc->SetBranchStatus("year",1);
    tree_mc->SetBranchStatus("*finalState",1);
    
    TFile* file_norm_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm.root");
    TTree* tree_norm_mc = (TTree*) file_norm_mc->Get("DecayTree");
    tree_norm_mc->SetBranchStatus("*",0);
    tree_norm_mc->SetBranchStatus("*TAU*",1);
    tree_norm_mc->SetBranchStatus("weight",1);
    tree_norm_mc->SetBranchStatus("year",1);
    tree_norm_mc->SetBranchStatus("*finalState",1);
    
    // Define observables
    RooRealVar Bs_TAU(((string)Bs_TAU_Var).c_str(), ((string)Bs_TAU_Var).c_str(), min_TAU, max_TAU, "ps");
    RooRealVar Bs_TAUERR(((string)Bs_TAU_Var+"ERR").c_str(), ((string)Bs_TAU_Var+"ERR").c_str(), 0.00001, 1.,"ps");
    RooRealVar weight("weight" , "weight", 0.);
    RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
    RooRealVar year("year", "year", 0.);
    RooArgSet observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, weight);
    
    // Define category to distinguish between singal and norm data
    RooCategory decay("decay","decay") ;
    decay.defineType("norm");
    decay.defineType("signal_mc");
    decay.defineType("norm_mc");
    
    RooDataSet* data_signal = new RooDataSet("data_signal","data_signal", observables,Import(*tree), WeightVar(weight.GetName()), Cut(CutString.c_str()));
    RooDataSet* data_norm = new RooDataSet("data_norm","data_norm", observables,Import(*tree_norm), WeightVar(weight.GetName()), Cut(CutString.c_str()));
    RooDataSet* data_signal_mc = new RooDataSet("data_signal_mc","data_signal_mc", observables,Import(*tree_mc), WeightVar(weight.GetName()), Cut(""));
    RooDataSet* data_norm_mc = new RooDataSet("data_norm_mc","data_norm_mc", observables,Import(*tree_norm_mc), WeightVar(weight.GetName()), Cut(""));
    
    RooDataSet* dataset = new RooDataSet("dataset","dataset",observables,Index(decay),Import("signal_mc",*data_signal_mc),Import("norm_mc",*data_norm_mc),Import("norm",*data_norm), WeightVar(weight.GetName()));
    
    /// SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
    
    // SPLINE KNOTS
    NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
    NamedParameter<double> knot_values("knot_values", 3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.0, 1.7757e+00);

    vector<double> myBinning;    
    vector<double> values;

    if(useAdaptiveBinningKnots){
	TH1D* binning = createBinning();	
	cout << endl << "knot positios: " << endl;
	for(int i = 1; i <= binning->GetNbinsX(); i++){
		cout << binning->GetBinCenter(i) << " , " ;
		myBinning.push_back(binning->GetBinCenter(i));
		values.push_back(1.);
	}
	
	cout << endl << endl;
    }
    else { 
	myBinning = knot_positions.getVector();
    	values = knot_values.getVector() ;
    }

    // Spline for DsKpipi acceptance    
    RooArgList tacc_list;
    for(int i= 0; i< values.size(); i++){
        tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 10.0)));
    }    
    tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(values.size())).c_str(), ("coeff_"+anythingToString(values.size())).c_str(), 1.0)));
    
    RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));
    
    tacc_list.add(*coeff_last);
    
    cout << myBinning.size() << endl;
    cout << values.size() << endl;
    tacc_list.Print();

    RooCubicSplineFun* spline_signal = new RooCubicSplineFun("spline_signal", "spline_signal", Bs_TAU, myBinning, tacc_list);
    
    // Spline for DsKpipi MC acceptance
    RooArgList mc_tacc_list;
    for(int i= 0; i< values.size(); i++){
        mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(i)).c_str(), ("mc_coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 10.0)));
    }    
    mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(values.size())).c_str(), ("mc_coeff_"+anythingToString(values.size())).c_str(), 1.0)));
    
    RooFormulaVar* mc_coeff_last = new RooFormulaVar(("mc_coeff_"+anythingToString(values.size()+1)).c_str(),("mc_coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *mc_tacc_list.find(("mc_coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));
    
    mc_tacc_list.add(*mc_coeff_last);
    
    RooCubicSplineFun* spline_signal_mc = new RooCubicSplineFun("spline_signal_mc", "spline_signal_mc", Bs_TAU, myBinning, mc_tacc_list);
    
    // Spline for Ds3pi/DsKpipi ratio
    vector<double> ratio_values(values.size(),1.) ;
    
    RooArgList ratio_tacc_list;
    for(int i= 0; i< ratio_values.size(); i++){
        ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(i)).c_str(), ("ratio_"+anythingToString(i)).c_str(), ratio_values[i], 0.0, 10.0)));
    }    
    ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(ratio_values.size())).c_str(), ("ratio_"+anythingToString(ratio_values.size())).c_str(), 1.0)));
    
    RooFormulaVar* ratio_last = new RooFormulaVar(("ratio_"+anythingToString(ratio_values.size()+1)).c_str(),("ratio_"+anythingToString(ratio_values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *ratio_tacc_list.find(("ratio_"+anythingToString(ratio_values.size()-1)).c_str()), RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));
    
    ratio_tacc_list.add(*ratio_last);
    
    RooCubicSplineFun* spline_ratio = new RooCubicSplineFun("spline_ratio", "spline_ratio", Bs_TAU, myBinning, ratio_tacc_list);
    
    // Spline for Ds3pi MC = ratio * DsKpipi MC acceptance
    RooSplineProduct* spline_norm_mc = new RooSplineProduct("spline_norm_mc","spline_norm_mc", Bs_TAU, *spline_signal_mc, *spline_ratio);

    // Spline for Ds3pi data = ratio * DsKpipi data acceptance
    RooSplineProduct* spline_norm = new RooSplineProduct("spline_norm","spline_norm", Bs_TAU, *spline_signal, *spline_ratio);
    
    /// Build simultaneous pdf
    RooRealVar trm_mean( "trm_mean" , "trm_mean", 0.0, "ps" );
    RooRealVar trm_scale( "trm_scale", "trm_scale", 1.20);
    RooRealVar trm_mean_mc( "trm_mean_mc" , "trm_mean_mc", 0.0, "ps" );
    RooRealVar trm_scale_mc( "trm_scale_mc", "trm_scale_mc", 1.20);
    
    RooGaussEfficiencyModel trm_signal("trm_signal", "trm_signal", Bs_TAU, *spline_signal, trm_mean, Bs_TAUERR, trm_mean, trm_scale );
    RooGaussEfficiencyModel trm_norm("trm_norm", "trm_norm", Bs_TAU, *spline_norm, trm_mean, Bs_TAUERR, trm_mean, trm_scale );
    RooGaussEfficiencyModel trm_signal_mc("trm_signal_mc", "trm_signal_mc", Bs_TAU, *spline_signal_mc, trm_mean, Bs_TAUERR, trm_mean_mc, trm_scale_mc );
    RooGaussEfficiencyModel trm_norm_mc("trm_norm_mc", "trm_norm_mc", Bs_TAU, *spline_norm_mc, trm_mean, Bs_TAUERR, trm_mean_mc, trm_scale_mc );

    RooBDecay* pdf_signal = new RooBDecay("pdf_signal", "pdf_signal", Bs_TAU, RooRealConstant::value(tau),
                                          RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                          RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs),
                                          trm_signal, RooBDecay::SingleSided);
    
    RooBDecay* pdf_norm = new RooBDecay("pdf_norm", "pdf_norm", Bs_TAU, RooRealConstant::value(tau),
                                        RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                        RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs),
                                        trm_norm, RooBDecay::SingleSided);
                                        
    RooBDecay* pdf_signal_mc = new RooBDecay("pdf_signal_mc", "pdf_signal_mc", Bs_TAU, RooRealConstant::value(tau_MC),
                                          RooRealConstant::value(dgamma_MC), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                          RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs_MC),
                                          trm_signal_mc, RooBDecay::SingleSided);
    
    RooBDecay* pdf_norm_mc = new RooBDecay("pdf_norm_mc", "pdf_norm_mc", Bs_TAU, RooRealConstant::value(tau_MC),
                                        RooRealConstant::value(dgamma_MC), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                        RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs_MC),
                                        trm_norm_mc, RooBDecay::SingleSided);
    
    RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simultaneous pdf",decay);
    simPdf->addPdf(*pdf_signal_mc,"signal_mc");
    simPdf->addPdf(*pdf_norm_mc,"norm_mc");
    simPdf->addPdf(*pdf_norm,"norm");

    /// Fit
    
    // Fit only signal MC to set reasonable start parameters
    pdf_signal_mc->fitTo(*data_signal_mc, Save(1), SumW2Error(kTRUE), NumCPU(numCPU),Optimize(2), Strategy(2),Extended(kFALSE));
    
    // Perform simulataneous fit
    RooFitResult *myfitresult = simPdf->fitTo(*dataset, Save(1), SumW2Error(kTRUE), NumCPU(numCPU),Optimize(2), Strategy(2),Extended(kFALSE));
    myfitresult->Print("v");
    
    /// Plot
    int bin = 80;
    
    vector<TString> decays;
    decays.push_back("signal_mc");
    decays.push_back("norm_mc");
    decays.push_back("signal");
    decays.push_back("norm");
        
    for(int i = 0 ; i < decays.size(); i++){
        TCanvas* canvas = new TCanvas();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);
        
        TLegend leg(0.65,0.5,0.9,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.06);
        leg.SetTextAlign(12);
        
        double max = 5.0 ;
        double min = -5.0 ;
        double rangeX = max-min;
        double zero = max/rangeX;
        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(max);
        graph->SetMinimum(min);
        graph->SetPoint(1,min_TAU,0);
        graph->SetPoint(2,max_TAU,0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(max);
        graph2->SetMinimum(min);
        graph2->SetPoint(1,min_TAU,-3);
        graph2->SetPoint(2,max_TAU,-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(max);
        graph3->SetMinimum(min);
        graph3->SetPoint(1,min_TAU,3);
        graph3->SetPoint(2,max_TAU,3);
        graph3->SetLineColor(kRed);
        
        RooPlot* frame_m = Bs_TAU.frame();	
        frame_m->GetXaxis()->SetLabelColor( kWhite);
        frame_m->GetYaxis()->SetTitleOffset(0.95);
        
        if(decays[i]=="signal"){
            data_signal->plotOn(frame_m, Binning(bin), Name("data_"+decays[i]));
            pdf_signal->plotOn(frame_m, LineColor(kBlue+1),  Name("pdf_"+decays[i]));
            spline_signal->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
        }
        else {
            dataset->plotOn(frame_m, Binning(bin), Name("data_"+decays[i]),Cut("decay==decay::"+decays[i]));
            simPdf->plotOn(frame_m, LineColor(kBlue+1),  Name("pdf_"+decays[i]),Slice(decay,decays[i]),ProjWData(decay,*dataset));
        }
        
        if(decays[i]=="signal_mc")spline_signal_mc->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
        
        else if(decays[i]=="norm_mc"){ 
            spline_norm_mc->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
            spline_ratio->plotOn(frame_m, LineColor(kMagenta+3), LineWidth(3), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_ratio_"+decays[i]));
            spline_signal_mc->plotOn(frame_m, LineColor(kGreen+3), LineStyle(kDashed), LineWidth(3), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_signal_"+decays[i]));
        }
        
        else if(decays[i]=="norm"){ 
            spline_norm->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
            spline_ratio->plotOn(frame_m, LineColor(kMagenta+3), LineWidth(3), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_ratio_"+decays[i]));
            spline_signal->plotOn(frame_m, LineColor(kGreen+3), LineWidth(3), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_signal_"+decays[i]));
        }
        
        if(decays[i]=="signal")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrow D_{s}K#pi#pi Data","ep");
        else if(decays[i]=="signal_mc")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrow D_{s}K#pi#pi MC","ep");
        else if(decays[i]=="norm")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrow D_{s}#pi#pi#pi Data","ep");
        else if(decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrow D_{s}#pi#pi#pi MC","ep");

        if(decays[i]=="signal")leg.AddEntry(frame_m->findObject("pdf_"+decays[i]),"Prediction","l");
        else leg.AddEntry(frame_m->findObject("pdf_"+decays[i]),"Fit","l");
        
        leg.AddEntry(frame_m->findObject("spline_"+decays[i]),"Acceptance","l");
        if(decays[i]=="norm")leg.AddEntry(frame_m->findObject("spline_signal_"+decays[i]),"Acc(t)_{D_{s}K#pi#pi}^{Data}","l");
        else if(decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("spline_signal_"+decays[i]),"Acc(t)_{D_{s}K#pi#pi}^{MC}","l");
        if(decays[i]=="norm" || decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("spline_ratio_"+decays[i]),"Ratio","l");
        
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame_m->GetYaxis()->SetRangeUser(0.011,frame_m->GetMaximum()*1.1);
        frame_m->Draw();
        leg.Draw();
        
        canvas->cd();
        TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2->SetBorderMode(0);
        pad2->SetBorderSize(-1);
        pad2->SetFillStyle(0);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        RooPlot* frame_p = Bs_TAU.frame();
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
        frame_p->GetXaxis()->SetTitle("#font[132]{t(B_{s}) [ps]}");
        
        RooHist* pullHist  = frame_m->pullHist("data_"+decays[i],"pdf_"+decays[i]);
        frame_p->addPlotable(pullHist,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->SaveAs("Plot/timeAccRatioFit_"+decays[i]+"_"+(string)BinningName+ ".eps");
        if(updateAnaNotePlots)canvas->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAccRatioFit_"+decays[i]+"_"+(string)BinningName+ ".pdf");
        
        pad1->SetLogy(1);
        pad1->Update();
        canvas->Update();
        canvas->SaveAs("Plot/timeAccRatioFit_"+decays[i]+"_"+(string)BinningName+ "_log.eps");
        //if(updateAnaNotePlots)canvas->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAccRatioFit_"+decays[i]+"_"+(string)BinningName+ "_log.pdf");
        pad1->SetLogy(0);
        
        // ???
        double chi2 = frame_m->chiSquare("pdf_"+decays[i],"data_"+decays[i],values.size());
        cout << "chi2 = " << chi2 << endl;
    }
    
    
    //put coefficients into table    
    ofstream datafile;
    if(updateAnaNotePlots) datafile.open ("../../../../../TD-AnaNote/latex/tables/SplineCoeffs_Bs2DsKPiPi_Data.tex");
    else datafile.open ("SplineCoeffs_Bs2DsKPiPi_Data.tex");
    datafile << "\\begin{table}[h]" << "\n";
    datafile << "\\centering" << "\n";
    datafile << "\\begin{tabular}{l l}" << "\n";
    datafile << "Parameter & Fit to $\\Bs\\to\\Ds\\pion\\pion\\pion$ data \\\\" << "\n";
    datafile << "\\hline" << "\n";

    for(int i= 0; i< values.size()+2; i++){        
        datafile << std::fixed << std::setprecision(3) << ("$v_{"+anythingToString(i)).c_str()<<"}$ & " << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal() << " $\\pm$ "  << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError() << "\\\\" << "\n"; 
    }

    datafile << "\\end{tabular}" << "\n";
    datafile << "\\caption{Summary of the obtained parameters from the acceptance fit to " <<CutString.c_str()<<  ".} " << "\n";
    datafile << "\\label{table: Splines}" << "\n";
    datafile << "\\end{table}" << "\n";
}


void compareAcceptance(){

    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 1);
    
    NamedParameter<string> BinningName("BinningName",(string)"default");
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
    vector<double> myBinning = knot_positions.getVector();
    const int nBins = myBinning.size()+2;
    
    // Compare different years
    vector<TString> cuts_year;
    cuts_year.push_back("(year == 11)");
    cuts_year.push_back("(year == 12)");
    cuts_year.push_back("(year == 15)");
    cuts_year.push_back("(year == 16)");  
    
    vector<TString> legend_year;
    legend_year.push_back("Year 11");
    legend_year.push_back("Year 12");
    legend_year.push_back("Year 15");
    legend_year.push_back("Year 16");
    
    vector<TString> legend_year_2;
    legend_year_2.push_back("Others");
    legend_year_2.push_back("Others");
    legend_year_2.push_back("Others");
    legend_year_2.push_back("Others");
    
    // Compare different Ds final states
    vector<TString> cuts_Ds;
    cuts_Ds.push_back("(Ds_finalState == 0)");
    cuts_Ds.push_back("(Ds_finalState == 1)");
    cuts_Ds.push_back("(Ds_finalState == 2)");
    cuts_Ds.push_back("(Ds_finalState == 3)");    
    
    vector<TString> legend_Ds;
    legend_Ds.push_back("D_{s}^{-} #rightarrow #phi^{0}(1020) #pi^{-}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow K^{*0}(892) K^{-}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow (K^{+} K^{-} #pi^{-})_{NR}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow #pi^{+} #pi^{-} #pi^{-}"); 
    
    vector<TString> legend_Ds_2 = legend_year_2;
    
    // Compare different runs
    vector<TString> cuts_run;
    cuts_run.push_back("(year == 11 || year == 12) ");
    
    vector<TString> legend_run;
    legend_run.push_back("Run 1");
    
    vector<TString> legend_run_2;
    legend_run_2.push_back("Run 2");
    
    // Compare different years per run 
    vector<TString> cuts_run_year;
    cuts_run_year.push_back("(year == 11) && year != 15 && year != 16 ");
    cuts_run_year.push_back("(year == 15) && year != 11 && year != 12 ");
    
    vector<TString> legend_run_year;
    legend_run_year.push_back("Year 11");
    legend_run_year.push_back("Year 15");
    
    vector<TString> legend_run_year_2;
    legend_run_year_2.push_back("Year 12");
    legend_run_year_2.push_back("Year 16");
    
    // Combine cuts into vector to iterate over
    vector< vector<TString> > cut_set;    
    cut_set.push_back(cuts_year);    
    cut_set.push_back(cuts_Ds);
    cut_set.push_back(cuts_run);    
    cut_set.push_back(cuts_run_year);    
    
    vector< vector<TString> > legend_title_set;
    legend_title_set.push_back(legend_year);
    legend_title_set.push_back(legend_Ds);
    legend_title_set.push_back(legend_run);
    legend_title_set.push_back(legend_run_year);
    
    vector< vector<TString> > legend_title_set_2;
    legend_title_set_2.push_back(legend_year_2);
    legend_title_set_2.push_back(legend_Ds_2);
    legend_title_set_2.push_back(legend_run_2);
    legend_title_set_2.push_back(legend_run_year_2);
    
    vector<TString> plot_titles;
    plot_titles.push_back("year");
    plot_titles.push_back("DsFinalState");
    plot_titles.push_back("run");
    plot_titles.push_back("run_year");
    
    for (int a = 0; a < cut_set.size(); a++) {
        
        TCanvas* c = new TCanvas();
        TCanvas *combinedCanvas = new TCanvas();
        
        vector<TString> cuts = cut_set[a];
        
        if(cuts.size()==2)combinedCanvas->Divide(2,1);
        else if(cuts.size()==3)combinedCanvas->Divide(3,1);
        else if(cuts.size()==4)combinedCanvas->Divide(2,2);
        
        for (int n = 0; n < cuts.size(); n++) {
            
            vector< vector<double> > myCoeffsAndErr = fitSplineAcc((string)cuts[n]);
            vector< vector<double> > myCoeffsAndErr_reversed = fitSplineAcc((string)(cuts[n].ReplaceAll("(","!(")));
            
            double x[nBins]; 
            double xerr[nBins]; 
            double y[nBins]; 
            double yerr[nBins]; 
            double y_reversed[nBins]; 
            double yerr_reversed[nBins]; 
            
            double chi2 = 0;
            
            for (int i= 0; i < nBins; i++) {
                
                if(i==0) x[0] = min_TAU;
                else if(i==nBins-1) x[nBins-1] = max_TAU;
                else x[i] = myBinning[i-1];
                
                xerr[i] = 0;        
                y[i] = myCoeffsAndErr[0][i];
                yerr[i] = myCoeffsAndErr[1][i];
                
                y_reversed[i] = myCoeffsAndErr_reversed[0][i];
                yerr_reversed[i] = myCoeffsAndErr_reversed[1][i];
                
                if(i < nBins - 2)chi2 += pow(y[i]-y_reversed[i],2)/(pow(yerr[i],2)+pow(yerr_reversed[i],2));
            }
            
            chi2 = chi2/((double)nBins-2.);
            
            TGraphErrors *KnotsVsCoeffs = new TGraphErrors(nBins, x,y,xerr,yerr);
            KnotsVsCoeffs->SetMaximum(2.5);
            KnotsVsCoeffs->SetTitle("; t(B_{s}) [ps]; v_{i}");
            
            TGraphErrors *KnotsVsCoeffs_reversed = new TGraphErrors(nBins, x,y_reversed,xerr,yerr_reversed);
            KnotsVsCoeffs_reversed->SetMarkerColor(kRed+1);
            KnotsVsCoeffs_reversed->SetLineColor(kRed+1);
            KnotsVsCoeffs_reversed->SetFillColor(kRed+1);
            KnotsVsCoeffs_reversed->SetFillStyle(3001);
            
            TLegend* leg = new TLegend(0.5,0.65,0.9,0.9,"");
            leg->SetLineStyle(0);
            leg->SetLineColor(0);
            leg->SetFillColor(0);
            leg->SetTextFont(22);
            leg->SetTextColor(1);
            leg->SetTextSize(0.06);
            leg->SetTextAlign(12);
            
            if(cuts.size()>1)combinedCanvas->cd(n+1);
            else combinedCanvas->cd();
            KnotsVsCoeffs->Draw("AP");
            KnotsVsCoeffs_reversed->Draw("PSAME");
            
            leg->AddEntry(KnotsVsCoeffs,legend_title_set[a][n],"ep");
            leg->AddEntry(KnotsVsCoeffs_reversed,legend_title_set_2[a][n],"ep");
            
            stringstream ss ;
            TString leg_chi2 = "#chi^{2} = ";
            ss << std::fixed << std::setprecision(2) << chi2 ;
            leg_chi2 += ss.str();    
            leg->AddEntry((TObject*)0, leg_chi2, "");
            
            leg->Draw();
            
            if(n == 0){
                KnotsVsCoeffs->SetMarkerColor(kBlack);
                KnotsVsCoeffs->SetLineColor(kBlack);
            }
            else if(n == 1){
                KnotsVsCoeffs->SetMarkerColor(kBlue);        
                KnotsVsCoeffs->SetLineColor(kBlue);        
            }
            else if(n == 2){
                KnotsVsCoeffs->SetMarkerColor(kGreen+3);
                KnotsVsCoeffs->SetLineColor(kGreen+3);
            }
            else if(n == 3){
                KnotsVsCoeffs->SetMarkerColor(kMagenta+3);
                KnotsVsCoeffs->SetLineColor(kMagenta+3);
            }
            else if(n == 4){
                KnotsVsCoeffs->SetMarkerColor(kGray);
                KnotsVsCoeffs->SetLineColor(kGray);
            }
            c->cd();
            if(n==0)KnotsVsCoeffs->Draw("AP");
            else KnotsVsCoeffs->Draw("PSAME");
            
        }
        
        c->Print("Plot/timeAcc_comparison_by_"+plot_titles[a]+"_"+(string)BinningName+ ".eps");
        //if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAcc_comparison_by_"+plot_titles[a]+"_"+(string)BinningName+ ".pdf");
        
        combinedCanvas->Print("Plot/timeAcc_combined_by_"+plot_titles[a]+"_"+(string)BinningName+ ".eps");
        if(updateAnaNotePlots)combinedCanvas->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAcc_combined_by_"+plot_titles[a]+"_"+(string)BinningName+ ".pdf");
    }

}


int main(int argc, char** argv){
    
    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    
    NamedParameter<int> CompareAcceptance("CompareAcceptance", 1);
    NamedParameter<int> FitSplineAccRatio("FitSplineAccRatio", 1);
    
    if(CompareAcceptance)compareAcceptance();

    if(FitSplineAccRatio)fitSplineAccRatio(" year == 11 || year == 12");

        
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
