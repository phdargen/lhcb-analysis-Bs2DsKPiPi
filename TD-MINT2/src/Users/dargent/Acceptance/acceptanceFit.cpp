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
#include <TMatrixDSymfwd.h>
#include <TDecompChol.h>
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
#include "RooHistPdf.h"
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
#include "RooNumIntConfig.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooNDKeysPdf.h"
#include "RooKeysPdf.h"
#include "RooBDecay.h"
#include "RooProdPdf.h"
#include "RooMultiVarGaussian.h"
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

/// HFLAV 19 values
double tau = 1.510;
double sigma_tau = 0.004;
double Gamma = 0.6624;
double sigma_Gamma = 0.0018;
double dgamma = 0.090; 
double sigma_dgamma = 0.005;
double deltaMs = 17.757;

double tau_B0 = 1.519;
double sigma_tau_B0 = 0.004;
double dgamma_B0 = 0.0; 
double deltaMd = 0.0;


/// HFLAV summer 17 values
// double tau = 1.509;
// double sigma_tau = 0.004;
// double Gamma = 0.6629;
// double sigma_Gamma = 0.0018;
// double dgamma = -0.088; 
// double sigma_dgamma = 0.006;
// double deltaMs = 17.757;

// double tau_B0 = 1.520;
// double sigma_tau_B0 = 0.004;
// double dgamma_B0 = 0.0; 
// double deltaMd = 0.0;

/// MC values (old decFile)
// double tau_MC = 1.510; 
// double dgamma_MC = .09166; 
// double deltaMs_MC = 17.8;  

/// MC values (new decFile)
double tau_MC = 1.512; 
double dgamma_MC = .1097; 
double deltaMs_MC = 17.8; 

TH1D* createBinning(){

	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");		
	NamedParameter<double> TAU_min("TAU_min", 0.4);		
	NamedParameter<double> TAU_max("TAU_max", 10);		
        NamedParameter<int> minEventsPerBin("minEventsPerBin", 1000); 
	int dim = 1;

	TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal_scaled.root");
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
                         HyperBinningAlgorithms::LIKELIHOOD, 
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

vector< TGraph* > fitSplineAcc(string CutString, string marginalPdfsPrefix = "", double offset_dt = 0. , double scale_dt = 1.2, bool plot = false, double offset_mean_dt = 0. ){

	// Options
    	NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);
	NamedParameter<string> BinningName("BinningName",(string)"default");
   	NamedParameter<int> makePlots("makePlots", 0);
	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");
	NamedParameter<double> min_TAU("min_TAU", 0.4);
	NamedParameter<double> max_TAU("max_TAU", 10.);
    	NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    	NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);
	NamedParameter<int> nBins("nBins", 100);
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> useAdaptiveBinningKnots("useAdaptiveBinningKnots", 0);
	NamedParameter<int> fixFirstKnot("fixFirstKnot", 0);

	// Read Dataset
    	TChain* tree=new TChain("DecayTree");
    	tree->Add( ((string)InputDir + "Data/norm_18.root").c_str());
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("*sw",1);
	tree->SetBranchStatus("weight*",1);
	//tree->SetBranchStatus("BDT*",1);
        //tree->SetBranchStatus("*BKG*",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*finalState",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);
//     	tree->SetBranchStatus("Bs_L0Global_TIS",1);
//    	tree->SetBranchStatus("Bs_*_TOS",1);
//    	tree->SetBranchStatus("*PT",1);
//     	tree->SetBranchStatus("m_*",1);

	//Define RooRealVar for observables
	RooRealVar Bs_TAU(((string)Bs_TAU_Var).c_str(), ((string)Bs_TAU_Var).c_str(), min_TAU, max_TAU, "ps");
	RooRealVar Bs_TAUERR(((string)Bs_TAU_Var+"ERR").c_str(), ((string)Bs_TAU_Var+"ERR").c_str(), min_TAUERR,max_TAUERR,"ps");
	RooRealVar N_Bs_sw("weight", "weight", 0.);
	RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
	RooRealVar year("year", "year", 0.);
	RooRealVar TriggerCat("TriggerCat", "TriggerCat", 0.);
    
	RooRealVar m_Kpipi("m_Kpipi", "m_Kpipi", 0.);
	RooRealVar DsDaughters_min_PT("DsDaughters_min_PT", "DsDaughters_min_PT", 0.);
	RooRealVar XsDaughters_min_PT("XsDaughters_min_PT", "XsDaughters_min_PT", 0.);
	
	RooRealVar run("run", "run", 0.);
	RooRealVar Bs_L0Global_TIS("Bs_L0Global_TIS", "Bs_L0Global_TIS", 0.);
	RooRealVar Bs_L0HadronDecision_TOS("Bs_L0HadronDecision_TOS", "Bs_L0HadronDecision_TOS", 0.);
	
	RooRealVar Bs_Hlt1TrackMVADecision_TOS("Bs_Hlt1TrackMVADecision_TOS", "Bs_Hlt1TrackMVADecision_TOS", 0.);
	RooRealVar Bs_Hlt1TwoTrackMVADecision_TOS("Bs_Hlt1TwoTrackMVADecision_TOS", "Bs_Hlt1TwoTrackMVADecision_TOS", 0.);
	RooRealVar Bs_Hlt1TrackAllL0Decision_TOS("Bs_Hlt1TrackAllL0Decision_TOS", "Bs_Hlt1TrackAllL0Decision_TOS", 0.);
	
	RooRealVar Bs_Hlt2Topo2BodyDecision_TOS("Bs_Hlt2Topo2BodyDecision_TOS", "Bs_Hlt2Topo2BodyDecision_TOS", 0.);
	RooRealVar Bs_Hlt2Topo3BodyDecision_TOS("Bs_Hlt2Topo3BodyDecision_TOS", "Bs_Hlt2Topo3BodyDecision_TOS", 0.);
	RooRealVar Bs_Hlt2Topo4BodyDecision_TOS("Bs_Hlt2Topo4BodyDecision_TOS", "Bs_Hlt2Topo4BodyDecision_TOS", 0.);
	RooRealVar Bs_Hlt2PhiIncPhiDecision_TOS("Bs_Hlt2PhiIncPhiDecision_TOS", "Bs_Hlt2PhiIncPhiDecision_TOS", 0.);
	
	RooRealVar Bs_Hlt2Topo2BodyBBDTDecision_TOS("Bs_Hlt2Topo2BodyBBDTDecision_TOS", "Bs_Hlt2Topo2BodyBBDTDecision_TOS", 0.);
	RooRealVar Bs_Hlt2Topo3BodyBBDTDecision_TOS("Bs_Hlt2Topo3BodyBBDTDecision_TOS", "Bs_Hlt2Topo3BodyBBDTDecision_TOS", 0.);
	RooRealVar Bs_Hlt2Topo4BodyBBDTDecision_TOS("Bs_Hlt2Topo4BodyBBDTDecision_TOS", "Bs_Hlt2Topo4BodyBBDTDecision_TOS", 0.);
	RooRealVar Bs_Hlt2IncPhiDecision_TOS("Bs_Hlt2IncPhiDecision_TOS", "Bs_Hlt2IncPhiDecision_TOS", 0.);
	
	RooArgList observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, N_Bs_sw, run, TriggerCat);
	RooArgList observables2(Bs_L0Global_TIS, Bs_L0HadronDecision_TOS, Bs_Hlt1TrackMVADecision_TOS, Bs_Hlt1TwoTrackMVADecision_TOS, Bs_Hlt1TrackAllL0Decision_TOS);
	RooArgList observables3(Bs_Hlt2Topo2BodyDecision_TOS, Bs_Hlt2Topo3BodyDecision_TOS, Bs_Hlt2Topo4BodyDecision_TOS, Bs_Hlt2Topo2BodyBBDTDecision_TOS, Bs_Hlt2Topo3BodyBBDTDecision_TOS, Bs_Hlt2Topo4BodyBBDTDecision_TOS,Bs_Hlt2IncPhiDecision_TOS,Bs_Hlt2PhiIncPhiDecision_TOS);
	RooArgList observables4(m_Kpipi,DsDaughters_min_PT,XsDaughters_min_PT);
	
// 	observables.add(observables2);
// 	observables.add(observables3);
// 	observables.add(observables4);

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
        for(int i= 0; i<= values.size(); i++){
		if(fixFirstKnot){
			if(i==0)tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), 1.0)));
			else tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i-1], 0.0, 5.0)));
		}
		else{
			if(i==values.size())tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), 1.0)));
			else tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 5.0)));
		}
	}
	
	RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(*tacc_list.find(("coeff_"+anythingToString(values.size())).c_str()), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));

	tacc_list.add(*coeff_last);	

	//CUBIC SPLINE FUNCTION 
 	RooCubicSplineFun* spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_TAU, myBinning, tacc_list);
		
	RooRealVar trm_mean( "trm_mean" , "trm_mean", offset_mean_dt );
	RooRealVar trm_offset( "trm_offset", "trm_offset", offset_dt);
	RooRealVar trm_scale( "trm_scale", "trm_scale", scale_dt);
	//RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *spl, trm_mean, Bs_TAUERR, trm_mean, trm_scale );
        RooFormulaVar dt_scaled( "dt_scaled","dt_scaled", "@0+@1*@2",RooArgList(trm_offset,trm_scale,Bs_TAUERR));
        RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *spl, trm_mean, dt_scaled, RooRealConstant::value(1.), RooRealConstant::value(1.) );
	
	RooRealVar GammaVal("Gamma","Gamma",Gamma);
	RooRealVar DeltaGamma("DeltaGamma","DeltaGamma",dgamma);
        RooFormulaVar Tau( "Tau","Tau", "1./@0",RooArgList(GammaVal));

	TMatrixDSym cov(2);
	cov[0][0] = pow(0.0020,2) ;
	cov[0][1] = -0.239 * 0.0020 * 0.006 ;
 	cov[1][0] = -0.239 * 0.0020 * 0.006 ;
	cov[1][1] = pow(0.006,2) ;

	RooMultiVarGaussian* gauss_cov = new RooMultiVarGaussian("gauss_cov","gauss_cov",RooArgList(GammaVal,DeltaGamma), RooArgList(RooRealConstant::value(Gamma),RooRealConstant::value(dgamma)), cov);


	RooBDecay* timePdf = new RooBDecay("Bdecay", "Bdecay", Bs_TAU, Tau, DeltaGamma,
			RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
			RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs),
			trm, RooBDecay::SingleSided);
	
        /// Marginal pdfs
	TFile* f_pdfs = new TFile("Mistag_pdfs.root","OPEN");
	TH1D* h_dt;
	if((string)marginalPdfsPrefix != "")h_dt = new TH1D( *((TH1D*) f_pdfs->Get(("h_dt_norm_"+(string)marginalPdfsPrefix).c_str())));
	else h_dt = new TH1D( *((TH1D*) f_pdfs->Get("h_dt_norm")));
	// RooHistPdf doesn't like negative or 0 bins, set them to a small positive number
	h_dt->Smooth();
	for(int i= 1 ; i<=h_dt->GetNbinsX(); i++){
	if(h_dt->GetBinContent(i) <= 0.)h_dt->SetBinContent(i,0.000000001*h_dt->GetMaximum());
	}
	RooDataHist* r_h_dt = new RooDataHist("r_h_dt","r_h_dt",Bs_TAUERR,h_dt);
	RooHistPdf* pdf_sigma_t = new RooHistPdf("pdf_sigma_t","pdf_sigma_t",Bs_TAUERR,*r_h_dt);
	f_pdfs->Close();

	RooProdPdf* totPdf= new RooProdPdf("totPdf","totPdf",RooArgSet(*pdf_sigma_t),Conditional(RooArgSet(*timePdf),RooArgSet(Bs_TAU)));
        
	///Fit and Print
// 	RooFitResult *myfitresult = totPdf->fitTo(*dataset, Save(1), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE),NumCPU(numCPU),ExternalConstraints(*gauss_cov));
	RooFitResult *myfitresult = totPdf->fitTo(*dataset, Save(1), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE),NumCPU(numCPU));
	myfitresult->Print("v");

	//put coefficients into vector
	vector<double> myCoeffs,myCoeffsErr;
	for(int i= 0; i< values.size()+2; i++){
		myCoeffs.push_back(((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal());
       		myCoeffsErr.push_back(((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError());
	}
        
	if(plot){
		/// Plot
		TCanvas* canvas = new TCanvas();
		canvas->cd();
		canvas->SetTopMargin(0.05);
		canvas->SetBottomMargin(0.05);
	
		TLegend leg(0.5,0.65,0.9,0.9,"");
		leg.SetLineStyle(0);
		leg.SetLineColor(0);
		leg.SetFillColor(0);
		leg.SetTextFont(22);
		leg.SetTextColor(1);
		leg.SetTextSize(0.07);
		leg.SetTextAlign(12);
	
		RooPlot* frame_m = Bs_TAU.frame();	
		frame_m->GetXaxis()->SetLabelColor( kWhite);
		frame_m->GetYaxis()->SetTitleOffset(0.95);
	
		dataset->plotOn(frame_m, Binning(nBins), Name("data"));
		totPdf->plotOn(frame_m, LineColor(kBlue+1),  Name("pdf"));
		//totPdf->plotOn(frame_m, LineColor(kBlue+1),  Name("pdf"),ProjWData(Bs_TAUERR,*dataset));
		spl->plotOn(frame_m, LineColor(kRed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline"));//,VisualizeError(*myfitresult));
		
		leg.AddEntry(frame_m->findObject("data"),"LHCb data unofficial","ep");
		leg.AddEntry(frame_m->findObject("pdf"),"Fit","l");
		leg.AddEntry(frame_m->findObject("spline"),"Acceptance","l");
		
		TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
		pad1->SetBorderMode(0);
		pad1->SetBorderSize(-1);
		pad1->SetBottomMargin(0.);
		pad1->Draw();
		pad1->cd();
		frame_m->GetYaxis()->SetRangeUser(0.00011,frame_m->GetMaximum()*1.1);
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
		canvas->SaveAs(("Plot/timeAccFit_"+(string)BinningName+ "_" + marginalPdfsPrefix + ".eps").c_str());
		
		pad1->SetLogy(1);
		pad1->Update();
		canvas->Update();
		canvas->SaveAs(("Plot/timeAccFit_"+(string)BinningName+ "_" + marginalPdfsPrefix + "_log.eps").c_str());
		
		// ???
		double chi2 = frame_m->chiSquare("pdf","data",values.size());
		cout << "chi2 = " << chi2 << endl;
		cout << "used datasets:     "<< CutString.c_str() << endl;
		
		ofstream resultsFile;
		resultsFile.open(("results_norm_" + marginalPdfsPrefix + "_" + (string)BinningName+ ".txt").c_str(),std::ofstream::trunc);
		resultsFile << "knot_positions " ;
		for(int i= 0; i< myBinning.size(); i++){
			resultsFile << myBinning[i] << " " ;
		}
		resultsFile << endl;
		for(int i= 0; i< myBinning.size(); i++){
			resultsFile << "c" + anythingToString(i) + "_" + marginalPdfsPrefix << "  " << 2 << "  " << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal() 
			<< "  " <<  ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError() << endl;
		}
	}

   	vector< vector<double> > myCoeffsAndErr;
    	myCoeffsAndErr.push_back(myCoeffs);
    	myCoeffsAndErr.push_back(myCoeffsErr);    

        // Plot acceptance
        TH1F *h_spline = new TH1F("", "", 100, min_TAU, max_TAU);
        for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
            Bs_TAU.setVal(h_spline->GetXaxis()->GetBinCenter(i));
            h_spline->SetBinContent(i,spl->getVal());
        }
        
        TCanvas* c = new TCanvas();
	c->cd();
        h_spline->SetLineColor(kRed);
        h_spline->Draw("histc");
        c->Print("spline.eps");
        c->Print("spline.pdf");


	TMatrixDSym covM = myfitresult->covarianceMatrix();
	covM.Print();

	    int nbins = myBinning.size()+2;
            double x[nbins]; 
            double xerr[nbins]; 
            double y[nbins]; 
            double yerr[nbins]; 
            
            for (int i= 0; i < nbins; i++) {
                
                if(i==0) x[0] = min_TAU;
                else if(i==nbins-1) x[nbins-1] = max_TAU;
                else x[i] = myBinning[i-1];
                
                xerr[i] = 0;        
                y[i] = myCoeffsAndErr[0][i];
                yerr[i] = myCoeffsAndErr[1][i];
            }
            
            TGraphErrors *KnotsVsCoeffs = new TGraphErrors(nbins, x,y,xerr,yerr);
            //KnotsVsCoeffs->SetMinimum(0.3 );
            //KnotsVsCoeffs->SetMaximum(1.5 );
            KnotsVsCoeffs->SetTitle("; t(B_{s}) [ps]; v_{i}");

            KnotsVsCoeffs->Draw("AP");


	    Bs_TAU.setVal(x[nbins-2]);
	    double y_fit_norm = spl->getVal();

	    nbins *= 100;

            double x_fit[nbins]; 
            double xerr_fit[nbins]; 
            double y_fit[nbins]; 
            double yerr_fit_h[nbins]; 
            double yerr_fit_l[nbins]; 
           
            for (int i= 0; i < nbins; i++) {
                
                x_fit[i] = min_TAU + (max_TAU - min_TAU)/(double)(nbins-1.) * (double)i;
                xerr_fit[i] = 0;        

	        Bs_TAU.setVal(x_fit[i]);
		y_fit[i] = spl->getVal() - y_fit_norm +1.;
            }
            
        int nIter = 100000;
	vector< TH1D* > randomHists;
	for (int i= 0; i < nbins; i++) randomHists.push_back(new TH1D("","",nIter,-1,1));

	for(int n=0 ; n < nIter; n++){

	    RooArgList randomPars = myfitresult->randomizePars();
	    for (int i= 0; i < myBinning.size(); i++)
		((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->setVal(((RooRealVar*)randomPars.at(i))->getVal());

	    for (int i= 0; i < nbins; i++) {
			Bs_TAU.setVal(x_fit[i]);
			randomHists[i]->Fill(spl->getVal() - y_fit_norm +1. - y_fit[i]);
	    }
	}
	
	for (int i= 0; i < nbins; i++) {
		double* integral = randomHists[i]->GetIntegral();
		for(int n= 0; n < randomHists[i]->GetNbinsX();n++){
			if(integral[n] > (0.5-0.341)) { 
				yerr_fit_l[i] = abs(randomHists[i]->GetBinCenter(n));
				break;
			}
		}
		for(int n= 0; n < randomHists[i]->GetNbinsX();n++){
			if(integral[n] > 1.-(0.5-0.341)) { 
				yerr_fit_h[i] = abs(randomHists[i]->GetBinCenter(n));
				break;
			}
		}
	}


	TGraphAsymmErrors *KnotsVsCoeffs_fit = new TGraphAsymmErrors(nbins, x_fit,y_fit,xerr_fit,xerr_fit,yerr_fit_l,yerr_fit_h);
	KnotsVsCoeffs_fit->SetMarkerColor(kRed);
	KnotsVsCoeffs_fit->SetLineColor(kRed);
	KnotsVsCoeffs_fit->SetFillColor(kRed);
	KnotsVsCoeffs_fit->SetFillStyle(3001);
	KnotsVsCoeffs_fit->Draw("3SAME");

	TGraph *KnotsVsCoeffs_fit2 = new TGraph(nbins, x_fit,y_fit);
	KnotsVsCoeffs_fit2->SetMarkerColor(kRed);
	KnotsVsCoeffs_fit2->SetLineColor(kRed);
	KnotsVsCoeffs_fit2->SetLineWidth(5);
//             KnotsVsCoeffs_fit->SetFillColor(kRed+1);
//             KnotsVsCoeffs_fit->SetFillStyle(3001);
	KnotsVsCoeffs_fit2->Draw("LSAME");


	c->Print("Plot/timeAcc_test.eps");
	c->Print("Plot/timeAcc_test.pdf");
	c->Print("Plot/timeAcc_test.png");

	
	vector< TGraph* > ret_vec;
	ret_vec.push_back(KnotsVsCoeffs);
	ret_vec.push_back(KnotsVsCoeffs_fit2);
	ret_vec.push_back(KnotsVsCoeffs_fit);

	return ret_vec;
}

RooFitResult * fitSplineAccRatio(string CutString, string CutStringMC, string marginalPdfsPrefix, double offset_dt, double scale_dt , double offset_dt_MC, double scale_dt_MC, double Tau, double DeltaGamma, bool plot = true, double offset_mean_dt = 0., double offset_mean_dt_MC = 0. ){
    
    // Options
    NamedParameter<string> BinningName("BinningName",(string)"default");
    NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);
    NamedParameter<int> nBins("nBins", 100);
    NamedParameter<int> numCPU("numCPU", 6);
    NamedParameter<int> fitB0("fitB0", 0);
    NamedParameter<int> fixRatio("fixRatio", 0);
    NamedParameter<int> fixFirstKnot("fixFirstKnot", 0);
    NamedParameter<int> useAdaptiveBinningKnots("useAdaptiveBinningKnots", 0);
    NamedParameter<int> updateAnaNote("updateAnaNote", 1);

    // Read Datasets
    TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/signal_18.root");
    TTree* tree = (TTree*) file->Get("DecayTree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("*TAU*",1);
    tree->SetBranchStatus("weight",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("*finalState",1);
    tree->SetBranchStatus("TriggerCat",1);
    tree->SetBranchStatus("run",1);    

    TFile* file_norm= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/norm_18.root");
    TTree* tree_norm = (TTree*) file_norm->Get("DecayTree");
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("weight",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*finalState",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("run",1);    

   // TFile* file_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal_18_newBDT.root");
    TFile* file_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal_scaled.root");
    TTree* tree_mc = (TTree*) file_mc->Get("DecayTree");
    tree_mc->SetBranchStatus("*",0);
    tree_mc->SetBranchStatus("*TAU*",1);
    tree_mc->SetBranchStatus("weight",1);
    tree_mc->SetBranchStatus("year",1);
    tree_mc->SetBranchStatus("*finalState",1);
    tree_mc->SetBranchStatus("TriggerCat",1);
    tree_mc->SetBranchStatus("run",1);    
    
    //TFile* file_norm_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm_18_newBDT.root");
    TFile* file_norm_mc= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm_scaled.root");
    TTree* tree_norm_mc = (TTree*) file_norm_mc->Get("DecayTree");
    tree_norm_mc->SetBranchStatus("*",0);
    tree_norm_mc->SetBranchStatus("*TAU*",1);
    tree_norm_mc->SetBranchStatus("weight",1);
    tree_norm_mc->SetBranchStatus("year",1);
    tree_norm_mc->SetBranchStatus("*finalState",1);
    tree_norm_mc->SetBranchStatus("TriggerCat",1);
    tree_norm_mc->SetBranchStatus("run",1);    
    
    // Define observables
    RooRealVar Bs_TAU(((string)Bs_TAU_Var).c_str(), "t", min_TAU, max_TAU, "ps");
    RooRealVar Bs_TAUERR(((string)Bs_TAU_Var+"ERR").c_str(), ((string)Bs_TAU_Var+"ERR").c_str(), min_TAUERR,max_TAUERR,"ps");

    RooRealVar B0_TAU("Bs_DTF_TAU", "t", min_TAU, max_TAU, "ps");
    RooRealVar B0_TAUERR("Bs_DTF_TAUERR", "Bs_DTF_TAUERR", min_TAUERR,max_TAUERR,"ps");

    RooRealVar weight("weight" , "weight", 0.);
    RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
    RooRealVar year("year", "year", 0.);
    RooRealVar run("run", "run", 0.);
    RooRealVar TriggerCat("TriggerCat", "TriggerCat", 0.);
    RooArgSet observables(Bs_TAU, Bs_TAUERR,B0_TAU, B0_TAUERR, Ds_finalState, year, run, weight, TriggerCat);
    
    // Define category to distinguish between singal and norm data
    RooCategory decay("decay","decay") ;
    decay.defineType("norm");
    decay.defineType("signal_B0");
    decay.defineType("signal_mc");
    decay.defineType("norm_mc");
    
    RooDataSet* data = new RooDataSet("data","data", observables,Import(*tree), WeightVar(weight.GetName()), Cut(CutString.c_str()));
    RooDataSet* data_norm = new RooDataSet("data_norm","data_norm", observables,Import(*tree_norm), WeightVar(weight.GetName()), Cut(CutString.c_str()));
    RooDataSet* data_signal_mc = new RooDataSet("data_signal_mc","data_signal_mc", observables,Import(*tree_mc), WeightVar(weight.GetName()), Cut(CutStringMC.c_str()));
    RooDataSet* data_norm_mc = new RooDataSet("data_norm_mc","data_norm_mc", observables,Import(*tree_norm_mc), WeightVar(weight.GetName()), Cut(CutStringMC.c_str()));
    
    RooDataSet* dataset = new RooDataSet("dataset","dataset",observables,Index(decay),Import("signal_B0",*data),Import("signal_mc",*data_signal_mc),Import("norm_mc",*data_norm_mc),Import("norm",*data_norm), WeightVar(weight.GetName()));
    
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
    for(int i= 0; i<= values.size(); i++){
	if(fixFirstKnot){
		if(i==0)tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), 1.0)));
		else tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i-1], 0.0, 5.0)));
    	}
	else {
		if(i==values.size())tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), 1.0)));
		else tacc_list.add(*(new RooRealVar(("coeff_"+anythingToString(i)).c_str(), ("coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 5.0)));
	}
    }
    RooFormulaVar* coeff_last = new RooFormulaVar(("coeff_"+anythingToString(values.size()+1)).c_str(),("coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(*tacc_list.find(("coeff_"+anythingToString(values.size())).c_str()), *tacc_list.find(("coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));

    tacc_list.add(*coeff_last);
    
    RooCubicSplineFun* spline_signal = new RooCubicSplineFun("spline_signal", "spline_signal", Bs_TAU, myBinning, tacc_list);
    RooCubicSplineFun* spline_signal_B0 = new RooCubicSplineFun("spline_signal_B0", "spline_signal_B0", B0_TAU, myBinning, tacc_list);
    
    // Spline for DsKpipi MC acceptance
    RooArgList mc_tacc_list;
    for(int i= 0; i<= values.size(); i++){
	if(fixFirstKnot){
	        if(i==0)mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(i)).c_str(), ("mc_coeff_"+anythingToString(i)).c_str(), 1.0)));
        	else mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(i)).c_str(), ("mc_coeff_"+anythingToString(i)).c_str(), values[i-1], 0.0, 5.0)));
	}
	else {
		if(i==values.size())mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(i)).c_str(), ("mc_coeff_"+anythingToString(i)).c_str(), 1.0)));
        	else mc_tacc_list.add(*(new RooRealVar(("mc_coeff_"+anythingToString(i)).c_str(), ("mc_coeff_"+anythingToString(i)).c_str(), values[i], 0.0, 5.0)));
	}
    }    
    RooFormulaVar* mc_coeff_last = new RooFormulaVar(("mc_coeff_"+anythingToString(values.size()+1)).c_str(),("mc_coeff_"+anythingToString(values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(*mc_tacc_list.find(("mc_coeff_"+anythingToString(values.size())).c_str()), *mc_tacc_list.find(("mc_coeff_"+anythingToString(values.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));
    
    mc_tacc_list.add(*mc_coeff_last);
    
    RooCubicSplineFun* spline_signal_mc = new RooCubicSplineFun("spline_signal_mc", "spline_signal_mc", Bs_TAU, myBinning, mc_tacc_list);
    
    // Spline for Ds3pi/DsKpipi ratio
    vector<double> ratio_values(values.size(),1.) ;
    
    RooArgList ratio_tacc_list;
    for(int i= 0; i<= ratio_values.size(); i++){
	if(fixFirstKnot){
	        if(i==0)ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(i)).c_str(), ("ratio_"+anythingToString(i)).c_str(), 1.0)));
        	else ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(i)).c_str(), ("ratio_"+anythingToString(i)).c_str(), ratio_values[i-1],0.,2.)));
	}
	else {
        	if(i==values.size())ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(i)).c_str(), ("ratio_"+anythingToString(i)).c_str(), 1.0)));
        	else ratio_tacc_list.add(*(new RooRealVar(("ratio_"+anythingToString(i)).c_str(), ("ratio_"+anythingToString(i)).c_str(), ratio_values[i],0.,2.)));
	}
    }        
    RooFormulaVar* ratio_last = new RooFormulaVar(("ratio_"+anythingToString(ratio_values.size()+1)).c_str(),("ratio_"+anythingToString(ratio_values.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(*ratio_tacc_list.find(("ratio_"+anythingToString(ratio_values.size())).c_str()), *ratio_tacc_list.find(("ratio_"+anythingToString(ratio_values.size()-1)).c_str()), RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(Bs_TAU.getMax()) ));
  
    ratio_tacc_list.add(*ratio_last);

    if(fixRatio){
		RooFIter iterat = ratio_tacc_list.fwdIterator();  
                RooAbsArg * next = 0;
                while(next=iterat.next()) ((RooRealVar*)next)->setConstant();
    }
    RooCubicSplineFun* spline_ratio = new RooCubicSplineFun("spline_ratio", "spline_ratio", Bs_TAU, myBinning, ratio_tacc_list);
    
    // Spline for Ds3pi MC = ratio * DsKpipi MC acceptance
    RooSplineProduct* spline_norm_mc = new RooSplineProduct("spline_norm_mc","spline_norm_mc", Bs_TAU, *spline_signal_mc, *spline_ratio);

    // Spline for Ds3pi data = ratio * DsKpipi data acceptance
    RooSplineProduct* spline_norm = new RooSplineProduct("spline_norm","spline_norm", Bs_TAU, *spline_signal, *spline_ratio);
    
    /// Build simultaneous pdf
    RooRealVar trm_mean( "trm_mean" , "trm_mean", offset_mean_dt );
    RooRealVar trm_offset( "trm_offset", "trm_offset", offset_dt);
    RooRealVar trm_scale( "trm_scale", "trm_scale", scale_dt);
    RooFormulaVar dt_scaled( "dt_scaled","dt_scaled", "@0+@1*@2",RooArgList(trm_offset,trm_scale,Bs_TAUERR));
    RooFormulaVar dt_scaled_B0( "dt_scaled_B0","dt_scaled_B0", "@0+@1*@2",RooArgList(trm_offset,trm_scale,B0_TAUERR));

    RooRealVar trm_mean_mc( "trm_mean_mc" , "trm_mean_mc", offset_mean_dt_MC );
    RooRealVar trm_offset_mc( "trm_offset_mc", "trm_offset_mc", offset_dt_MC);
    RooRealVar trm_scale_mc( "trm_scale_mc", "trm_scale_mc", scale_dt_MC);
    RooFormulaVar dt_scaled_mc( "dt_scaled_mc","dt_scaled_mc", "@0+@1*@2",RooArgList(trm_offset_mc,trm_scale_mc,Bs_TAUERR));
    
    RooGaussEfficiencyModel trm_signal_B0("trm_signal_B0", "trm_signal_B0", B0_TAU, *spline_signal_B0, trm_mean, dt_scaled_B0,  RooRealConstant::value(1.), RooRealConstant::value(1.) );
    RooGaussEfficiencyModel trm_norm("trm_norm", "trm_norm", Bs_TAU, *spline_norm,trm_mean, dt_scaled,  RooRealConstant::value(1.), RooRealConstant::value(1.) );
    RooGaussEfficiencyModel trm_signal_mc("trm_signal_mc", "trm_signal_mc", Bs_TAU, *spline_signal_mc, trm_mean_mc, dt_scaled_mc,  RooRealConstant::value(1.), RooRealConstant::value(1.) );
    RooGaussEfficiencyModel trm_norm_mc("trm_norm_mc", "trm_norm_mc", Bs_TAU, *spline_norm_mc, trm_mean_mc, dt_scaled_mc,  RooRealConstant::value(1.), RooRealConstant::value(1.) );

    // time pdfs
    RooBDecay* time_pdf_signal_B0 = new RooBDecay("time_pdf_signal_B0", "time_pdf_signal_B0", B0_TAU, RooRealConstant::value(tau_B0),
                                        RooRealConstant::value(dgamma_B0), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                        RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMd),
                                        trm_signal_B0, RooBDecay::SingleSided);
                                        
    RooBDecay* time_pdf_norm = new RooBDecay("time_pdf_norm", "time_pdf_norm", Bs_TAU, RooRealConstant::value(Tau),
                                        RooRealConstant::value(DeltaGamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                        RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs),
                                        trm_norm, RooBDecay::SingleSided);
                                        
    RooBDecay* time_pdf_signal_mc = new RooBDecay("time_pdf_signal_mc", "time_pdf_signal_mc", Bs_TAU, RooRealConstant::value(tau_MC),
                                          RooRealConstant::value(dgamma_MC), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                          RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs_MC),
                                          trm_signal_mc, RooBDecay::SingleSided);
    
    RooBDecay* time_pdf_norm_mc = new RooBDecay("time_pdf_norm_mc", "time_pdf_norm_mc", Bs_TAU, RooRealConstant::value(tau_MC),
                                        RooRealConstant::value(dgamma_MC), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                                        RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(deltaMs_MC),
                                        trm_norm_mc, RooBDecay::SingleSided);

    /// Marginal pdfs
    TFile* f_pdfs = new TFile("Mistag_pdfs.root","OPEN");
    TH1D* h_dt;
    if((string)marginalPdfsPrefix != "")h_dt = new TH1D( *((TH1D*) f_pdfs->Get(("h_dt_norm_"+(string)marginalPdfsPrefix).c_str())));
    else h_dt = new TH1D( *((TH1D*) f_pdfs->Get("h_dt_norm")));
    // RooHistPdf doesn't like negative or 0 bins, set them to a small positive number
    h_dt->Smooth();
    for(int i= 1 ; i<=h_dt->GetNbinsX(); i++){
	if(h_dt->GetBinContent(i) <= 0.)h_dt->SetBinContent(i,0.000000001*h_dt->GetMaximum());
    }
    RooDataHist* r_h_dt = new RooDataHist("r_h_dt","r_h_dt",Bs_TAUERR,h_dt);
    RooHistPdf* pdf_sigma_t = new RooHistPdf("pdf_sigma_t","pdf_sigma_t",Bs_TAUERR,*r_h_dt);

    RooDataHist* r_h_dt_B0 = new RooDataHist("r_h_dt_B0","r_h_dt_B0",B0_TAUERR,h_dt);
    RooHistPdf* pdf_sigma_t_B0 = new RooHistPdf("pdf_sigma_t_B0","pdf_sigma_t_B0",B0_TAUERR,*r_h_dt_B0);
    f_pdfs->Close();

    // total pdfs
    RooProdPdf* pdf_signal_B0= new RooProdPdf("pdf_signal_B0","pdf_signal_B0",RooArgSet(*pdf_sigma_t_B0),Conditional(RooArgSet(*time_pdf_signal_B0),RooArgSet(B0_TAU)));
    RooProdPdf* pdf_norm= new RooProdPdf("pdf_norm","pdf_norm",RooArgSet(*pdf_sigma_t),Conditional(RooArgSet(*time_pdf_norm),RooArgSet(Bs_TAU)));
    RooProdPdf* pdf_signal_mc= new RooProdPdf("pdf_signal_mc","pdf_signal_mc",RooArgSet(*pdf_sigma_t),Conditional(RooArgSet(*time_pdf_signal_mc),RooArgSet(Bs_TAU)));
    RooProdPdf* pdf_norm_mc= new RooProdPdf("pdf_norm_mc","pdf_norm_mc",RooArgSet(*pdf_sigma_t),Conditional(RooArgSet(*time_pdf_norm_mc),RooArgSet(Bs_TAU)));
    
    RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simultaneous pdf",decay);
    simPdf->addPdf(*pdf_signal_mc,"signal_mc");
    simPdf->addPdf(*pdf_norm_mc,"norm_mc");
    simPdf->addPdf(*pdf_norm,"norm");
    if(fitB0)simPdf->addPdf(*pdf_signal_B0,"signal_B0");

    /// Fit
    RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-7);
    RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-7);
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("extrapolation","WynnEpsilon");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("maxSteps","1000");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("minSteps","0");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","21Points");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 1000);

    /// Fit only signal MC to set reasonable start parameters
    pdf_signal_mc->fitTo(*data_signal_mc, Save(1), SumW2Error(kTRUE), NumCPU(numCPU),Extended(kFALSE));
//     pdf_signal_B0->fitTo(*data, Save(1), SumW2Error(kTRUE), NumCPU(numCPU),Extended(kFALSE));
    
    /// Perform simulataneous fit
    RooFitResult *myfitresult = simPdf->fitTo(*dataset, Save(1), SumW2Error(kTRUE), NumCPU(numCPU),Extended(kFALSE));
    myfitresult->Print("v");
    
    if(!plot)return myfitresult;

    /// Plot    
    vector<TString> decays;
    decays.push_back("signal_mc");
    decays.push_back("norm_mc");
    decays.push_back("norm");
    decays.push_back("signal_B0");
            
    for(int i = 0 ; i < decays.size(); i++){
        TCanvas* canvas = new TCanvas();
//         canvas->SetTopMargin(0.05);
//         canvas->SetBottomMargin(0.05);
	double legY = 0.6;
	if(A_is_in_B("norm",(string)decays[i])) legY = 0.4;
        TLegend leg(0.5,legY,0.9,0.9,"");
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
	if(decays[i]=="signal_B0") frame_m = B0_TAU.frame();	
//         frame_m->GetXaxis()->SetLabelColor( kWhite);
/*        frame_m->GetYaxis()->SetTitleOffset(0.95);*/
        frame_m->GetYaxis()->SetTitleOffset(1.5);

	if(decays[i]=="signal_B0" && !fitB0){
        	data->plotOn(frame_m, Binning(nBins/3), Name("data_"+decays[i]));
		pdf_signal_B0->plotOn(frame_m, LineColor(kBlue+1), LineWidth(4),  Name("pdf_"+decays[i]));
	}
	else{
		if(decays[i]=="signal_B0") dataset->plotOn(frame_m, Binning(25), Name("data_"+decays[i]),Cut("decay==decay::"+decays[i]));
		else dataset->plotOn(frame_m, Binning(nBins), Name("data_"+decays[i]),Cut("decay==decay::"+decays[i]));
        	simPdf->plotOn(frame_m, LineColor(kBlue+1), LineWidth(4),  Name("pdf_"+decays[i]),Slice(decay,decays[i]),ProjWData(decay,*dataset));
	}
        if(decays[i]=="signal_mc")spline_signal_mc->plotOn(frame_m, LineWidth(4), LineColor(kGreen+2), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));

        else if(decays[i]=="signal_B0")spline_signal_B0->plotOn(frame_m, LineWidth(4), LineColor(kRed+1), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
     
        else if(decays[i]=="norm_mc"){ 
            spline_norm_mc->plotOn(frame_m, LineColor(kBlack), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
            spline_ratio->plotOn(frame_m, LineColor(kMagenta+3), LineWidth(4), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_ratio_"+decays[i]));
            spline_signal_mc->plotOn(frame_m, LineColor(kGreen+2), LineStyle(kDashed), LineWidth(4), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_signal_"+decays[i]));
        }
        
        else if(decays[i]=="norm"){ 
            spline_norm->plotOn(frame_m, LineWidth(4), LineColor(kOrange+7), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_"+decays[i]));
            spline_ratio->plotOn(frame_m, LineColor(kMagenta+3), LineWidth(4), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_ratio_"+decays[i]));
            spline_signal->plotOn(frame_m, LineColor(kRed+1), LineWidth(4), LineStyle(kDashed), Normalization(frame_m->GetMaximum()*0.25, RooAbsReal::NumEvent),Name("spline_signal_"+decays[i]));
        }
        
        if(decays[i]=="signal_mc")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrowD_{s}K#pi#pi MC","ep");
        else if(decays[i]=="signal_B0")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{d}#rightarrowD_{s}K#pi#pi Data","ep");
	else if(decays[i]=="norm")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrowD_{s}#pi#pi#pi Data","ep");
        else if(decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("data_"+decays[i]),"B_{s}#rightarrowD_{s}#pi#pi#pi MC","ep");

        leg.AddEntry(frame_m->findObject("pdf_"+decays[i]),"Fit","l");
        
        if(decays[i]=="signal_mc")leg.AddEntry(frame_m->findObject("spline_"+decays[i]),"#varepsilon_{D_{s}K#pi#pi}^{MC}(t)","l");
        else if(decays[i]=="signal_B0")leg.AddEntry(frame_m->findObject("spline_"+decays[i]),"#varepsilon_{D_{s}K#pi#pi}^{Data}(t)","l");
	else if(decays[i]=="norm")leg.AddEntry(frame_m->findObject("spline_"+decays[i]),"#varepsilon_{D_{s}#pi#pi#pi}^{Data}(t)","l");
        else if(decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("spline_"+decays[i]),"#varepsilon_{D_{s}#pi#pi#pi}^{MC}(t)","l");

        if(decays[i]=="norm")leg.AddEntry(frame_m->findObject("spline_signal_"+decays[i]),"#varepsilon_{D_{s}K#pi#pi}^{Data}(t)","l");
        else if(decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("spline_signal_"+decays[i]),"#varepsilon_{D_{s}K#pi#pi}^{MC}(t)","l");
        if(decays[i]=="norm" || decays[i]=="norm_mc")leg.AddEntry(frame_m->findObject("spline_ratio_"+decays[i]),"R(t)","l");


        frame_m->SetMinimum(0.0);
        frame_m->GetYaxis()->SetTitleSize(0.05);

        frame_m->Draw();
        leg.Draw();
        canvas->SaveAs("Plot/timeAccRatioFit2_"+decays[i]+"_"+ marginalPdfsPrefix + "_" + (string)BinningName+ ".eps");

        double chi2 = frame_m->chiSquare("pdf_"+decays[i],"data_"+decays[i],values.size());
        cout << "chi2 = " << chi2 << endl;
        
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame_m->GetYaxis()->SetRangeUser(0.00011,frame_m->GetMaximum()*1.1);
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
        frame_p->GetXaxis()->SetTitle("#font[132]{t[ps]}");
        
        RooHist* pullHist  = frame_m->pullHist("data_"+decays[i],"pdf_"+decays[i]);
        frame_p->addPlotable(pullHist,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->SaveAs("Plot/timeAccRatioFit_"+decays[i]+"_"+ marginalPdfsPrefix + "_" + (string)BinningName+ ".eps");
        if(updateAnaNote)canvas->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/"+(string)BinningName+"/timeAccRatioFit_"+decays[i]+"_"+ marginalPdfsPrefix + ".pdf");
        
        pad1->SetLogy(1);
        pad1->Update();
        canvas->Update();
        canvas->SaveAs("Plot/timeAccRatioFit_"+decays[i]+"_"+marginalPdfsPrefix + "_" +(string)BinningName+ "_log.eps");
        //if(updateAnaNote)canvas->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAccRatioFit_"+decays[i]+"_"+(string)BinningName+ "_log.pdf");
        pad1->SetLogy(0);
    }
    
    //put coefficients into table    
    ofstream datafile;
    if(updateAnaNote) datafile.open(("../../../../../TD-AnaNote/latex/tables/Acceptance/"+(string)BinningName+"/splineCoeffs_"+ marginalPdfsPrefix + ".tex").c_str(),std::ofstream::trunc);
    else datafile.open(("splineCoeffs_"+ marginalPdfsPrefix + "_" + (string)BinningName+ ".tex").c_str(),std::ofstream::trunc);
//     datafile << "\\begin{table}[hp!]" << "\n";
//     datafile << "\\centering" << "\n";
//     datafile << "\\small" << "\n";
//     datafile << "\\caption{Time acceptance parameters for ";
//     if(CutString == "") datafile << " all events. }" << "\n";
//     else datafile << "events in category [";
//     if(A_is_in_B("run == 1", CutString)) datafile << "\\textsf{Run-I}"; 
//     else if(A_is_in_B("run == 2", CutString)) datafile << "\\textsf{Run-II}"; 
//     if(A_is_in_B("&&", CutString)) datafile << ","; 
//     if(A_is_in_B("TriggerCat == 0", CutString)) datafile << "\\textsf{L0-TOS}"; 
//     else if(A_is_in_B("TriggerCat == 1", CutString)) datafile << "\\textsf{L0-TIS}"; 
//     if(CutString != "") datafile << "].}" << "\n";
    datafile << "\\begin{tabular}{c c c c c}" << "\n";
    datafile << "\\hline" << "\n";
    datafile << "\\hline" << "\n";
    datafile << "Knot position & Coefficient & $\\Bs\\to\\Ds\\kaon\\pion\\pion$ data & $\\Bs\\to\\Ds\\kaon\\pion\\pion$ MC & Ratio \\\\" << "\n";
    datafile << "\\hline" << "\n";
    for(int i= 0; i< values.size()+2; i++){        
	double knot_pos;
	if(i==0)knot_pos = Bs_TAU.getMin();
	else if(i==values.size()+1)knot_pos= Bs_TAU.getMax();
	else knot_pos = myBinning[i-1];
        datafile << std::setprecision(1) << knot_pos << " & ";
	datafile << std::fixed << std::setprecision(3) << ("$v_{"+anythingToString(i)).c_str()<<"}$ & ";
	if(i < values.size()){
		datafile << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal() << " $\\pm$ "  << 	((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError()  <<  " & " <<
		((RooRealVar*)mc_tacc_list.find(("mc_coeff_"+anythingToString(i)).c_str()))->getVal() << " $\\pm$ "  << ((RooRealVar*)mc_tacc_list.find(("mc_coeff_"+anythingToString(i)).c_str()))->getError()  <<  " & " <<
		((RooRealVar*)ratio_tacc_list.find(("ratio_"+anythingToString(i)).c_str()))->getVal() << " $\\pm$ "  << 			((RooRealVar*)ratio_tacc_list.find(("ratio_"+anythingToString(i)).c_str()))->getError() 	
		<< "\\\\" << "\n";
	} 
	else if (i==values.size()){
		datafile << " 1.0 (fixed) & 1.0 (fixed) & 1.0 (fixed)" << "\\\\" << "\n";
	}
	else {
		datafile << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal() << " (interpolated) & "  << 
		((RooRealVar*)mc_tacc_list.find(("mc_coeff_"+anythingToString(i)).c_str()))->getVal() << " (interpolated) & "  << 
		((RooRealVar*)ratio_tacc_list.find(("ratio_"+anythingToString(i)).c_str()))->getVal() << " (interpolated) "  <<  "\\\\" << "\n";	
		}
    }
    datafile << "\\hline" << "\n";
    datafile << "\\hline" << "\n";
    datafile << "\\end{tabular}" << "\n";
//     datafile << "\\label{table:splines}" << "\n";
//     datafile << "\\end{table}";

    ofstream resultsFile;
    resultsFile.open(("results_" + marginalPdfsPrefix + "_" + (string)BinningName+ ".txt").c_str(),std::ofstream::trunc);
    resultsFile << "knot_positions " ;
    for(int i= 0; i< myBinning.size(); i++){
	resultsFile << myBinning[i] << " " ;
    }
    resultsFile << endl;
    for(int i= 0; i< myBinning.size(); i++){
	resultsFile << "c" + anythingToString(i) + "_" + marginalPdfsPrefix << "  " << 2 << "  " << ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getVal() 
	<< "  " <<  ((RooRealVar*)tacc_list.find(("coeff_"+anythingToString(i)).c_str()))->getError() << endl;
    }

    return myfitresult;
}

void compareAcceptance(){

    NamedParameter<int> updateAnaNote("updateAnaNote", 1);
    
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
    cuts_year.push_back("(year == 17)");  
    
    vector<TString> legend_year;
    legend_year.push_back("Year 11");
    legend_year.push_back("Year 12");
    legend_year.push_back("Year 15");
    legend_year.push_back("Year 16");
    legend_year.push_back("Year 17");


    vector<TString> cuts_year_t0;
    cuts_year_t0.push_back("(year == 11)&& TriggerCat == 0");
    cuts_year_t0.push_back("(year == 12)&& TriggerCat == 0");
    cuts_year_t0.push_back("(year == 15)&& TriggerCat == 0");
    cuts_year_t0.push_back("(year == 16) && TriggerCat == 0 ");  
    cuts_year_t0.push_back("(year == 17) && TriggerCat == 0 ");  

    vector<TString> cuts_year_t1;
    cuts_year_t1.push_back("(year == 11)&& TriggerCat == 1");
    cuts_year_t1.push_back("(year == 12)&& TriggerCat == 1");
    cuts_year_t1.push_back("(year == 15)&& TriggerCat == 1");
    cuts_year_t1.push_back("(year == 16) && TriggerCat == 1 ");  
    cuts_year_t1.push_back("(year == 17) && TriggerCat == 1 ");  
    

    // Compare different Ds final states
    vector<TString> cuts_Ds;
    cuts_Ds.push_back("(Ds_finalState == 0)");
    cuts_Ds.push_back("(Ds_finalState == 1)");
    cuts_Ds.push_back("(Ds_finalState == 2)");
    cuts_Ds.push_back("(Ds_finalState == 3)");    
    cuts_Ds.push_back("(Ds_finalState == 4)");    

    vector<TString> legend_Ds;
    legend_Ds.push_back("D_{s}^{-} #rightarrow #phi^{0}(1020) #pi^{-}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow K^{*0}(892) K^{-}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow (K^{+} K^{-} #pi^{-})_{NR}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow #pi^{+} #pi^{-} #pi^{-}"); 
    legend_Ds.push_back("D_{s}^{-} #rightarrow K^{-} #pi^{+} #pi^{-}"); 


    vector<TString> cuts_Ds_mod;
    cuts_Ds_mod.push_back("(Ds_finalState < 3)");
    cuts_Ds_mod.push_back("(Ds_finalState == 3)");    
    cuts_Ds_mod.push_back("(Ds_finalState == 4)");    

    vector<TString> legend_Ds_mod;
    legend_Ds_mod.push_back("D_{s}^{-} #rightarrow (K^{+} K^{-} #pi^{-})"); 
    legend_Ds_mod.push_back("D_{s}^{-} #rightarrow #pi^{+} #pi^{-} #pi^{-}"); 
    legend_Ds_mod.push_back("D_{s}^{-} #rightarrow K^{-} #pi^{+} #pi^{-}"); 
    

    // Compare different runs
    vector<TString> cuts_run;
    cuts_run.push_back("( run ==1 ) ");
    cuts_run.push_back("( run ==2 ) ");
    
    vector<TString> legend_run;
    legend_run.push_back("Run-I");
    legend_run.push_back("Run-II");
    

    // Compare different triggers
    vector<TString> cuts_trigger;
    cuts_trigger.push_back("(TriggerCat == 0)");
    cuts_trigger.push_back("(TriggerCat == 1)");
    
    vector<TString> legend_trigger;
    legend_trigger.push_back("LO-TOS");
    legend_trigger.push_back("LO-TIS");


    /// Combine cuts into vector to iterate over
    vector< vector<TString> > cut_set;    
//     cut_set.push_back(cuts_year);
//     cut_set.push_back(cuts_year_t0);
    cut_set.push_back(cuts_year_t1);    
//      cut_set.push_back(cuts_Ds);
//     cut_set.push_back(cuts_Ds_mod);
//     cut_set.push_back(cuts_run);    
//     cut_set.push_back(cuts_trigger);    
    
    vector< vector<TString> > legend_title_set;
//     legend_title_set.push_back(legend_year);
//      legend_title_set.push_back(legend_year);
     legend_title_set.push_back(legend_year);
//     legend_title_set.push_back(legend_Ds);
//     legend_title_set.push_back(legend_Ds_mod);
//     legend_title_set.push_back(legend_run);
//     legend_title_set.push_back(legend_trigger);
        
    vector<TString> plot_titles;
//      plot_titles.push_back("year");
//     plot_titles.push_back("year_t0");
     plot_titles.push_back("year_t1");
//     plot_titles.push_back("DsFinalState");
//     plot_titles.push_back("DsFinalState_mod");
//     plot_titles.push_back("run");
//     plot_titles.push_back("trigger");    

    for (int a = 0; a < cut_set.size(); a++) {
        
        TCanvas* c = new TCanvas();    
	c->cd();    
        vector<TString> cuts = cut_set[a];
        
    	TLegend* leg = new TLegend(0.6,0.7,0.9,0.9,"");
        leg->SetLineStyle(0);
	leg->SetLineColor(0);
	leg->SetFillColor(0);
	leg->SetTextFont(22);
	leg->SetTextColor(1);
	leg->SetTextSize(0.045);
	leg->SetTextAlign(12);

	vector< vector< TGraph* > > myGraphs_vec;
        for (int n = 0; n < cuts.size(); n++) {
            
            vector< TGraph* > myGraphs = fitSplineAcc((string)cuts[n],"comb",0.,1.);
	    myGraphs_vec.push_back(myGraphs);
	    c->cd();                          
            
	    int color = kRed;
	    if(n == 1) color = kBlue;
	    if(n == 2) color = kGreen;
	    if(n == 3) color = kMagenta;
	    if(n == 4) color = kGray+3;

	    myGraphs[0]->SetMinimum(0.);
	    myGraphs[0]->SetMaximum(2.);
            myGraphs[0]->SetMarkerColor(color+2);
            myGraphs[0]->SetLineColor(color+2);
            if(n==0)myGraphs[0]->Draw("AP");
	    else myGraphs[0]->Draw("PSAME");

            myGraphs[2]->SetMarkerColor(color);
            myGraphs[2]->SetLineColor(color);
            myGraphs[2]->SetFillColor(color);
            myGraphs[2]->SetFillStyle(1001);
	    myGraphs[2]->Draw("3SAME");

            myGraphs[1]->SetMarkerColor(color);
            myGraphs[1]->SetLineColor(color);
            myGraphs[1]->SetLineWidth(5);
	    myGraphs[1]->Draw("LSAME");

            leg->AddEntry(myGraphs[1],legend_title_set[a][n],"ep");                
        }
	for (int n = 0; n < cuts.size(); n++)myGraphs_vec[n][0]->Draw("PSAME");       

        leg->Draw();
        c->Print("Plot/timeAcc_comparison_by_"+plot_titles[a]+"_"+(string)BinningName+ ".eps");
        if(updateAnaNote)c->Print("../../../../../TD-AnaNote/latex/figs/Acceptance/timeAcc_comparison_by_"+plot_titles[a]+"_"+(string)BinningName+ ".pdf");
    }

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

void checkPV(){

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/BDT/", (char*) 0);
    NamedParameter<int> nBins("nBins", 10);

    int run,year,Ds_finalState,trigger,Bs_DTF_nPV;
    double t;
    double w;
    int cat,yearMC,Ds_finalStateMC;
    float Bs_DTF_chi2[100];
    double Bs_TRUEORIGINVERTEX_Z,Bs_OWNPV_Z,Bs_OWNPV_ZERR;
    double Bs_TRUEORIGINVERTEX_Y,Bs_OWNPV_Y,Bs_OWNPV_YERR;
    double Bs_TRUEORIGINVERTEX_X,Bs_OWNPV_X,Bs_OWNPV_XERR;

    TChain* treeMC =new TChain("DecayTree");
    treeMC->Add( ((string)InputDir + "MC/signal_PIDMC.root").c_str());
    //treeMC->Add( ((string)InputDir + "MC/norm.root").c_str());
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("Bs_BKGCAT",1);
    treeMC->SetBranchStatus("weight",1);
    treeMC->SetBranchStatus("*TAU*",1);
    treeMC->SetBranchStatus("*Bs_TRUEORIGINVERTEX_Z",1);
    treeMC->SetBranchStatus("*PV*",1);
    
    treeMC->SetBranchAddress("Bs_DTF_TAU",&t);
    treeMC->SetBranchAddress("Bs_BKGCAT",&cat);
    treeMC->SetBranchAddress("year",&yearMC);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalStateMC);           
    treeMC->SetBranchAddress("weight",&w);
    treeMC->SetBranchAddress("Bs_TRUEORIGINVERTEX_X",&Bs_TRUEORIGINVERTEX_X);
    treeMC->SetBranchAddress("Bs_OWNPV_X",&Bs_OWNPV_X);
    treeMC->SetBranchAddress("Bs_OWNPV_XERR",&Bs_OWNPV_XERR);
    treeMC->SetBranchAddress("Bs_TRUEORIGINVERTEX_Y",&Bs_TRUEORIGINVERTEX_Y);
    treeMC->SetBranchAddress("Bs_OWNPV_Y",&Bs_OWNPV_Y);
    treeMC->SetBranchAddress("Bs_OWNPV_YERR",&Bs_OWNPV_YERR);
    treeMC->SetBranchAddress("Bs_TRUEORIGINVERTEX_Z",&Bs_TRUEORIGINVERTEX_Z);
    treeMC->SetBranchAddress("Bs_OWNPV_Z",&Bs_OWNPV_Z);
    treeMC->SetBranchAddress("Bs_OWNPV_ZERR",&Bs_OWNPV_ZERR);
    treeMC->SetBranchAddress("Bs_DTF_chi2",&Bs_DTF_chi2);
    treeMC->SetBranchAddress("Bs_DTF_nPV",&Bs_DTF_nPV);


    TH1D* h_t = new TH1D("h_t",";t [ps];Events (norm.) ",50,0,15);
    TH1D* h_t_cut = new TH1D("h_t_cut",";t [ps];Events (norm.) ",50,0,15);

    TH1D* h_t_nPV1 = new TH1D("h_t_nPV1",";t [ps];Events (norm.) ",nBins,0,15);
    TH1D* h_t_nPV2 = new TH1D("h_t_nPV2",";t [ps];Events (norm.) ",nBins,0,15);

    TH1D* h_t_wrongPV = new TH1D("h_t_wrongPV",";t [ps];Events (norm.) ",nBins,0,15);
    TH1D* h_t_rightPV = new TH1D("h_t_rightPV",";t [ps];Events (norm.) ",nBins,0,15);

    TH1D* h_t_wrongPV_cut = new TH1D("h_t_wrongPV",";t [ps];Events (norm.) ",nBins,0,15);
    TH1D* h_t_rightPV_cut = new TH1D("h_t_rightPV",";t [ps];Events (norm.) ",nBins,0,15);
    
    TH1D* h_chi2_right = new TH1D("h_chi2_right",";#Delta #chi^{2}_{DTF};Events (norm.) ",50,0,500);
    TH1D* h_chi2_wrong = new TH1D("h_chi2_wrong",";#Delta #chi^{2}_{DTF};Events (norm.) ",50,0,500);

    ///loop over data events
    for(int i=0; i< treeMC->GetEntries(); i++)
    {    
        //if (0ul == (i % 1000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
        treeMC->GetEntry(i);        
        
        if(Bs_DTF_nPV == 1) h_t_nPV1->Fill(t,w);
        else{ 
            h_t->Fill(t,w/(exp(-t/tau_MC)*cosh(dgamma_MC/2.*t)));
            if(Bs_DTF_chi2[1]-Bs_DTF_chi2[0]>1000)h_t_cut->Fill(t,w/(exp(-t/tau_MC)*cosh(dgamma_MC/2.*t)));

            h_t_nPV2->Fill(t,w);
            if( abs(Bs_TRUEORIGINVERTEX_Z-Bs_OWNPV_Z)/Bs_OWNPV_ZERR > 5 && abs(Bs_TRUEORIGINVERTEX_X-Bs_OWNPV_X)/Bs_OWNPV_XERR > 5 && abs(Bs_TRUEORIGINVERTEX_Y-Bs_OWNPV_Y)/Bs_OWNPV_YERR > 5) {
                h_t_wrongPV->Fill(t,w);
                h_chi2_wrong->Fill(Bs_DTF_chi2[1]-Bs_DTF_chi2[0],w);
                if(Bs_DTF_chi2[1]-Bs_DTF_chi2[0]>15)h_t_wrongPV_cut->Fill(t,w);
            }
            else {
                h_t_rightPV->Fill(t,w);
                h_chi2_right->Fill(Bs_DTF_chi2[1]-Bs_DTF_chi2[0],w);
                if(Bs_DTF_chi2[1]-Bs_DTF_chi2[0]>15)h_t_rightPV_cut->Fill(t,w);
                }
        }
    }
    
    double N_right =  h_t_rightPV->Integral();
    double N_right_cut =  h_t_rightPV_cut->Integral();
    double N_wrong =  h_t_wrongPV->Integral();
    double N_wrong_cut =  h_t_wrongPV_cut->Integral();

    cout << "N_right = " << N_right << endl;
    cout << "N_wrong = " << N_wrong << endl;

    cout << "eff signal = " << N_right_cut/N_right << endl;
    cout << "bkg rejection = " << 1.-N_wrong_cut/N_wrong << endl;

    TCanvas*c = new TCanvas();
    
    h_chi2_right->Draw();
    h_chi2_wrong->SetLineColor(kBlue);
    h_chi2_wrong->Draw("e1same");
    c->Print("h_chi2.eps");
    
    h_t_wrongPV->Divide(h_t_wrongPV,h_t_rightPV);
    h_t_wrongPV->SetMinimum(0);
    h_t_wrongPV->Draw("e1");
    c->Print("h_t_PV.eps");
    
    h_t_wrongPV_cut->SetLineColor(kBlue);
    h_t_wrongPV_cut->Divide(h_t_wrongPV_cut,h_t_rightPV_cut);
    h_t_wrongPV_cut->SetMinimum(0);
    h_t_wrongPV_cut->Draw("e1");
    c->Print("h_t_PV_cut.eps");
    
    h_t_nPV2->Divide(h_t_nPV2,h_t_nPV1);
    h_t_nPV2->SetMinimum(0);
    h_t_nPV2->Draw("e1");
    c->Print("h_t_nPV2.eps");
    
    h_t->Draw("e1");
    c->Print("h_t.eps");
    h_t_cut->Draw("e1");
    c->Print("h_t_cut.eps");
    
    h_t_cut->Divide(h_t_cut,h_t);
    h_t_cut->SetMinimum(0);
    h_t_cut->Draw("e1");
    c->Print("h_t_ratio.eps");
    
    throw "";
    
}

int main(int argc, char** argv){
    
    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    
    NamedParameter<string> BinningName("BinningName",(string)"default");
    NamedParameter<int> CompareAcceptance("CompareAcceptance", 1);
    NamedParameter<int> FitSplineAccRatio("FitSplineAccRatio", 1);
    NamedParameter<int> FitSplineNorm("FitSplineNorm", 1);
    NamedParameter<int> doSystematics("doSystematics", 0);

    NamedParameter<double> offset_mean_dt("offset_mean_dt", 0.0);
    NamedParameter<double> offset_sigma_dt("offset_sigma_dt", 0.0);
    NamedParameter<double> scale_sigma_dt("scale_sigma_dt", 1.2);
    NamedParameter<double> offset_sigma_dt_MC("offset_sigma_dt_MC", 0.0);
    NamedParameter<double> scale_sigma_dt_MC("scale_sigma_dt_MC", 1.2);

    NamedParameter<double> offset_mean_dt_Run1("offset_mean_dt_Run1", 0.0);
    NamedParameter<double> offset_sigma_dt_Run1("offset_sigma_dt_Run1", 0.0);
    NamedParameter<double> scale_sigma_dt_Run1("scale_sigma_dt_Run1", 1.2);
    NamedParameter<double> offset_sigma_dt_Run1_MC("offset_sigma_dt_Run1_MC", 0.0);
    NamedParameter<double> scale_sigma_dt_Run1_MC("scale_sigma_dt_Run1_MC", 1.2);

    NamedParameter<double> offset_mean_dt_Run2("offset_mean_dt_Run2", 0.0);
    NamedParameter<double> offset_sigma_dt_Run2("offset_sigma_dt_Run2", 0.0);
    NamedParameter<double> scale_sigma_dt_Run2("scale_sigma_dt_Run2", 1.);
    NamedParameter<double> offset_sigma_dt_Run2_MC("offset_sigma_dt_Run2_MC", 0.0);
    NamedParameter<double> scale_sigma_dt_Run2_MC("scale_sigma_dt_Run2_MC", 1.);

    NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);

    //checkPV(); 
    produceMarginalPdfs();

    if(CompareAcceptance)compareAcceptance();

    if(FitSplineAccRatio){

	vector<string> dataCuts,mcCuts,marginalPdfs;
	dataCuts.push_back(" run == 1 && TriggerCat == 0 ");
	mcCuts.push_back(" run == 1 && TriggerCat == 0 ");
	marginalPdfs.push_back("Run1_t0");

	dataCuts.push_back(" run == 1 && TriggerCat == 1 ");
	mcCuts.push_back(" run == 1 && TriggerCat == 1 ");
	marginalPdfs.push_back("Run1_t1");

	dataCuts.push_back(" run == 2 && TriggerCat == 0 ");
	mcCuts.push_back("run == 2 && TriggerCat == 0");
	marginalPdfs.push_back("Run2_t0");

	dataCuts.push_back(" run == 2 && TriggerCat == 1 ");
	mcCuts.push_back("run == 2 && TriggerCat == 1");
	marginalPdfs.push_back("Run2_t1");

// 	dataCuts.push_back(" (year == 15 || year == 16) && TriggerCat == 0 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 0");
// 	marginalPdfs.push_back("Run2_t0");
// 	dataCuts.push_back(" (year == 15 || year == 16) && TriggerCat == 1 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 1");
// 	marginalPdfs.push_back("Run2_t1");
// 
// 	dataCuts.push_back(" (year == 17) && TriggerCat == 0 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 0");
// 	marginalPdfs.push_back("Run2_17_t0");
// 	dataCuts.push_back(" (year == 17) && TriggerCat == 1 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 1");
// 	marginalPdfs.push_back("Run2_17_t1");
// 
// 	dataCuts.push_back(" (year == 18) && TriggerCat == 0 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 0");
// 	marginalPdfs.push_back("Run2_18_t0");
// 	dataCuts.push_back(" (year == 18) && TriggerCat == 1 ");
// 	mcCuts.push_back("run == 2 && TriggerCat == 1");
// 	marginalPdfs.push_back("Run2_18_t1");


	//dataCuts.push_back("");
	//mcCuts.push_back("");
	//marginalPdfs.push_back("comb");

	vector<TMatrixDSym*> cors_fit;
	vector< vector<double>  > cors_Gamma;
	vector< vector<double>  > cors_deltaGamma;

	for(int opt = 0 ; opt < dataCuts.size(); opt++){
		/// Baseline fit	
		RooFitResult* result = fitSplineAccRatio(dataCuts[opt],mcCuts[opt], marginalPdfs[opt], (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, (double) offset_sigma_dt_Run2_MC, (double) scale_sigma_dt_Run2_MC,1./Gamma,dgamma);
	
		TMatrixDSym cor(knot_positions.getVector().size());
		vector<double> vals,errs;
	
		for(int i = 0; i < knot_positions.getVector().size(); i++){
			vals.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getVal());
			errs.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getError());
			for(int j = 0; j < knot_positions.getVector().size(); j++){
				cor[i][j] = result->correlation(("coeff_"+anythingToString(i)).c_str(),("coeff_"+anythingToString(j)).c_str());
			}
		}
		cors_fit.push_back(new TMatrixDSym(cor));
	
		/// Gamma +/- 1 sigma fits
		vector<double> vals_sigma_Gamma_p,vals_sigma_Gamma_m;
	
		if(doSystematics)result = fitSplineAccRatio(dataCuts[opt],mcCuts[opt], marginalPdfs[opt], (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, (double) offset_sigma_dt_Run2_MC, (double) scale_sigma_dt_Run2_MC,1./(Gamma+sigma_Gamma),dgamma,false);
	
		for(int i = 0; i < knot_positions.getVector().size(); i++)
			vals_sigma_Gamma_p.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getVal());
		
	
	
		if(doSystematics)result = fitSplineAccRatio(dataCuts[opt],mcCuts[opt], marginalPdfs[opt], (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, (double) offset_sigma_dt_Run2_MC, (double) scale_sigma_dt_Run2_MC,1./(Gamma-sigma_Gamma),dgamma,false);
	
		for(int i = 0; i < knot_positions.getVector().size(); i++)
			vals_sigma_Gamma_m.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getVal());
		
	
		/// deltaGamma +/- 1 sigma fits
		vector<double> vals_sigma_deltaGamma_p,vals_sigma_deltaGamma_m;
	
		if(doSystematics)result = fitSplineAccRatio(dataCuts[opt],mcCuts[opt], marginalPdfs[opt], (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, (double) offset_sigma_dt_Run2_MC, (double) scale_sigma_dt_Run2_MC,1./Gamma,dgamma+sigma_dgamma,false);
	
		for(int i = 0; i < knot_positions.getVector().size(); i++)
			vals_sigma_deltaGamma_p.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getVal());
		
	
	
		if(doSystematics)result = fitSplineAccRatio(dataCuts[opt],mcCuts[opt], marginalPdfs[opt], (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, (double) offset_sigma_dt_Run2_MC, (double) scale_sigma_dt_Run2_MC,1./Gamma,dgamma-sigma_dgamma,false);
	
		for(int i = 0; i < knot_positions.getVector().size(); i++)
			vals_sigma_deltaGamma_m.push_back( ((RooRealVar*)result->floatParsFinal().find(("coeff_"+anythingToString(i)).c_str()))->getVal());
	

		/// Calculate correlation
		vector<double> cor_Gamma;
		for(int i = 0; i < vals.size(); i++)cor_Gamma.push_back((vals_sigma_Gamma_p[i]-vals_sigma_Gamma_m[i])/(2.*errs[i]));
		cors_Gamma.push_back(cor_Gamma);

		vector<double> cor_deltaGamma;
		for(int i = 0; i < vals.size(); i++)cor_deltaGamma.push_back((vals_sigma_deltaGamma_p[i]-vals_sigma_deltaGamma_m[i])/(2.*errs[i]));
		cors_deltaGamma.push_back(cor_deltaGamma);

		cor.Print();
		cout << cor_Gamma << endl;
		cout << cor_deltaGamma << endl << endl;	
	}
		
	TMatrixDSym cor(knot_positions.getVector().size()*dataCuts.size()+2);
	cor[0][0] = 1.;
	cor[1][1] = 1.;
// 	cor[0][1] = 0.11; /// sign ???
// 	cor[1][0] = 0.11;
	cor[0][1] = -0.08; /// sign ???
	cor[1][0] = -0.08;

	for(int opt = 0 ; opt < dataCuts.size(); opt++){ 
			for(int i = 2+opt*knot_positions.getVector().size(); i < 2+(opt+1)*knot_positions.getVector().size(); i++){
				cor[0][i] = cors_Gamma[opt][i-(2+opt*knot_positions.getVector().size())];
				cor[1][i] = cors_deltaGamma[opt][i-(2+opt*knot_positions.getVector().size())];
				cor[i][0] = cors_Gamma[opt][i-(2+opt*knot_positions.getVector().size())];
				cor[i][1] = cors_deltaGamma[opt][i-(2+opt*knot_positions.getVector().size())];
			}
	}

	for(int opt = 0 ; opt < dataCuts.size(); opt++)
			for(int i = 2+opt*knot_positions.getVector().size(); i < 2+(opt+1)*knot_positions.getVector().size(); i++)
				for(int j = 2+opt*knot_positions.getVector().size(); j < 2+(opt+1)*knot_positions.getVector().size(); j++)
					cor[i][j] = (*cors_fit[opt])[i-(2+opt*knot_positions.getVector().size())][j-(2+opt*knot_positions.getVector().size())];

	cor.Print();	

	/// Save correlations
	if(doSystematics){ 
		ofstream correlationFile;
		correlationFile.open(("Correlations_"+ (string)BinningName+".txt").c_str(),std::ofstream::trunc);
		correlationFile << "\"" << "ConstrainMulti_Acc" << "\"" << "     " << "\"" << "Gamma dGamma ";
		for(int opt = 0 ; opt < dataCuts.size(); opt++)	
			for(int i = 0; i < knot_positions.getVector().size(); i++)
				correlationFile << "c" << i << "_" + marginalPdfs[opt] << " ";
	
		correlationFile << "\"";
		correlationFile<< "\n";
		correlationFile<< "\n";
		correlationFile << "\"" << "ConstrainMulti_Acc_corr" << "\"" << "     " << "\"" ;
		for(int i= 0; i< cor.GetNcols(); i++){
			for(int j= i; j< cor.GetNcols(); j++){
				correlationFile << cor[i][j] << " ";
			}
		}
		correlationFile << "\"";
	}

    }

    if(FitSplineNorm){
	//fitSplineAcc("" , "comb", 0.01 , 1.,true);
 //	fitSplineAcc(" run == 2" , "Run2", (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, true, offset_mean_dt_Run2);
 	fitSplineAcc(" run == 1 && TriggerCat == 0 " , "Run1_t0", (double) offset_sigma_dt_Run1, (double) scale_sigma_dt_Run1, true, offset_mean_dt_Run1);
 	fitSplineAcc(" run == 1 && TriggerCat == 1 " , "Run1_t1", (double) offset_sigma_dt_Run1, (double) scale_sigma_dt_Run1, true, offset_mean_dt_Run1);
 	fitSplineAcc(" run == 2 && TriggerCat == 0 " , "Run2_t0", (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, true,offset_mean_dt_Run2);
 	fitSplineAcc(" run == 2 && TriggerCat == 1 " , "Run2_t1", (double) offset_sigma_dt_Run2, (double) scale_sigma_dt_Run2, true, offset_mean_dt_Run2);
    }

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
