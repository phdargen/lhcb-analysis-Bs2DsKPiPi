// Resolution studies
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCut.h>
#include <TFile.h>
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
#include <TGraphErrors.h>
#include "TGraphAsymmErrors.h"
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
#include "RooConstVar.h"
#include "RooRealConstant.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooJohnsonSU.h"
#include "RooKeysPdf.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"
#include "Mint/HyperBinningPainter1D.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace MINT;

TFile* file_res = 0;
TTree* tree_res = 0;

vector<double> FitTimeRes(double min, double max, string binName = "", TString Bs_TAU_Var = "Bs_DTF_TAU", TString dataType = "MC"){
	
	/// Options
	NamedParameter<int> updateAnaNote("updateAnaNote", 1);
	NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
        NamedParameter<double> TAUERR_max("TAUERR_max", 0.15);
	NamedParameter<int> useTransformedSigma("useTransformedSigma", 0);

	NamedParameter<int> nGauss("nGauss", 2);
        NamedParameter<int> fixMean("fixMean", 0);
        NamedParameter<int> fixGaussFraction("fixGaussFraction", 0);
        NamedParameter<double> gaussFraction("gaussFraction", 0.75);

	/// Load file
        RooRealVar Bs_TAU(Bs_TAU_Var, Bs_TAU_Var, -20.,20.);
        RooRealVar Bs_TAUERR(Bs_TAU_Var+"ERR", Bs_TAU_Var+"ERR", min, max,"ps");
	RooRealVar Bs_TRUETAU("Bs_TRUETAU", "Bs_TRUETAU", 0.,20.);
	RooRealVar weight("weight" , "weight", 0.);
        RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
        RooRealVar year("year", "year", 0.);	

	RooArgList list =  RooArgList(Bs_TAU,Bs_TAUERR,weight,Ds_finalState,year);
	if(dataType == "MC")list.add(Bs_TRUETAU);
	RooDataSet* data = new RooDataSet("data","data",list,Import(*tree_res),WeightVar("weight"));

        /// Add residuals to dataset
	RooFormulaVar* Bs_DeltaTau_func;
	if(dataType=="MC")Bs_DeltaTau_func = new RooFormulaVar("Bs_DeltaTau_func","#Deltat","@0 - @1 * 1000.",RooArgList(Bs_TAU,Bs_TRUETAU));
	else Bs_DeltaTau_func = new RooFormulaVar("Bs_DeltaTau_func","t","@0",RooArgList(Bs_TAU));	
	RooRealVar* Bs_DeltaTau = (RooRealVar*) data->addColumn(*Bs_DeltaTau_func);

	/// Pdf
	RooRealVar* mean1 = new RooRealVar("mean1", "mean1", 0., -0.1*data->mean(Bs_TAUERR), 0.1*data->mean(Bs_TAUERR));
	if(fixMean)mean1->setConstant();
	RooRealVar* f = new RooRealVar("f" , "f", (double)gaussFraction,0.,1.);
	if(fixGaussFraction)f->setConstant();
	RooRealVar* f2 = new RooRealVar("f2" , "f2", 0.5, 0., 1.);

	RooRealVar *sigma1, *sigma2;
	RooRealVar *scale, *sigma_average, *sigma_RMS;
	if(useTransformedSigma){
		sigma_average = new RooRealVar("sigma_average", "sigma_average", (max-min)/2.,min,max);
		sigma_RMS = new RooRealVar("sigma_RMS", "sigma_RMS", (max-min),0,2.*max);
		sigma1 = (RooRealVar*) new RooFormulaVar("sigma1","@0 - @1 * sqrt(@2/(1.-@2))", RooArgList(*sigma_average,*sigma_RMS,*f));
		sigma2 = (RooRealVar*) new RooFormulaVar("sigma2","@0 + @1 * sqrt((1.-@2)/@2)", RooArgList(*sigma_average,*sigma_RMS,*f));
	}
	else {
		sigma1 = new RooRealVar("sigma1", "sigma1", 0.020,0.,0.1);
		scale  = new RooRealVar("scale", "scale", 2.,1.,10.);
		sigma2 = (RooRealVar*) new RooFormulaVar("sigma2", "@0*@1", RooArgList(*scale,*sigma1));
		//RooRealVar sigma2("sigma2", "sigma2", 0.045,0.,0.2);
	}
	RooRealVar* sigma3 = new RooRealVar("sigma3", "sigma3", 0.040,0.,0.2);

	RooGaussian Gauss1("Gauss1", "Gauss1", *Bs_DeltaTau, *mean1, *sigma1);
	RooGaussian Gauss2("Gauss2", "Gauss2", *Bs_DeltaTau, *mean1, *sigma2);
	RooGaussian Gauss3("Gauss3", "Gauss3", *Bs_DeltaTau, *mean1, *sigma3);
	
	RooAddPdf* pdf;
	if(nGauss<3)pdf = new RooAddPdf("pdf", "pdf", RooArgList(Gauss1,Gauss2),RooArgList(*f));
	else pdf = new RooAddPdf("pdf", "pdf", RooArgList(Gauss1,Gauss2,Gauss3),RooArgList(*f,*f2));
	if(nGauss==1){
		f->setVal(1);
		f->setConstant();
		sigma2->setConstant();
		scale->setConstant();
	}	

        /// Fit
	double fitRange_min, fitRange_max;
	if(dataType=="MC"){
		fitRange_min = -0.2;
		fitRange_max = 0.2;
	}
	else {
		fitRange_min = -0.25;
		fitRange_max = 0.;
		//fitRange_min = -4.*data->mean(Bs_TAUERR);
		//fitRange_max = 0.5*data->mean(Bs_TAUERR);
	}
	Bs_DeltaTau->setRange(fitRange_min,fitRange_max);

	RooFitResult *result = pdf->fitTo(*data,Save(kTRUE),SumW2Error(kTRUE),Extended(kFALSE),NumCPU(3),Range(fitRange_min,fitRange_max));
	cout << "result is --------------- "<<endl;
	result->Print();
	double covmatr = result->covQual();
	double edm = result->edm();
	double status = result->status();
	cout << endl <<"Status = " << status << " Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;

	/// Calculate effective resolution
	double dms = 17.757;
	RooFormulaVar dilution("dilution","@0 * exp(-@1*@1*@2*@2/2.) + (1. - @0) * exp(-@3*@3*@2*@2/2.)",RooArgList(*f,*sigma1,RooRealConstant::value(dms),*sigma2));
	double dilution_val = dilution.getVal();
	double dilution_error = dilution.getPropagatedError(*result);
	double f1 = f->getVal();
	double df1 = f->getError();
	double sig1 = sigma1->getVal();
	double dsig1 = sigma1->getError();
	double sig2 = sigma2->getVal();
	double dsig2 = sigma2->getPropagatedError(*result);

	RooFormulaVar resolution_eff("resolution_eff","sqrt(-2./@0/@0*log(@1))",RooArgList(RooRealConstant::value(dms),dilution)); 
	//double resolution_eff = sqrt(-2./pow(dms,2)*log(dilution_val));
	//double resolution_eff_error = -1./(dms*dilution_val*sqrt(log(1./pow(dilution_val,2))))*dilution_error;
	//double resolution_eff_error = ((2/(dms*dms))/(2*dilution_val*TMath::Sqrt((-2/(dms*dms))*log(dilution_val)))) * dilution_error;
	cout << "Measured resolution from dilution:   " << resolution_eff.getVal()*1000 << " +/- " << resolution_eff.getPropagatedError(*result)*1000 <<" fs" << endl;

	/*
	if(dilution_error > dilution_val/20.){
		cout << endl << "ERROR:: ERROR suspicously high, fit probably failed, will repeat fit fixed fraction of gaussians" << endl << endl;
		f_GaussBs.setConstant();
		result = DoubleGaussBs.fitTo(*data,Save(kTRUE),SumW2Error(kTRUE),Extended(kFALSE),NumCPU(3),Range(-0.2,0.2));
		cout << "result is --------------- "<<endl;
		result->Print();
		covmatr = result->covQual();
		edm = result->edm();
		status = result->status();
		cout << endl <<"Status = " << status << " Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
		dilution_val = dilution.getVal();
	        dilution_error = dilution.getPropagatedError(*result);
	}
	*/

	/// Plot
	int bin = 50;	
	RooPlot* frame_m= Bs_DeltaTau->frame();
	frame_m->SetTitle("");
	
	frame_m->GetXaxis()->SetLabelSize( 0.06 );
	frame_m->GetYaxis()->SetLabelSize( 0.06 );
	frame_m->GetXaxis()->SetLabelFont( 132 );
	frame_m->GetYaxis()->SetLabelFont( 132 );
	frame_m->GetXaxis()->SetLabelOffset( 0.006 );
	frame_m->GetYaxis()->SetLabelOffset( 0.006 );
	frame_m->GetXaxis()->SetLabelColor( kWhite);
	
	frame_m->GetXaxis()->SetTitleSize( 0.06 );
	frame_m->GetYaxis()->SetTitleSize( 0.1 );
	//frame_m->GetYaxis()->SetNdivisions(512);
	frame_m->GetXaxis()->SetTitleOffset( 1.00 );
	frame_m->GetYaxis()->SetTitleOffset( 0.5 );
	frame_m->GetYaxis()->SetTitle("Yield [norm.]");

	data->plotOn(frame_m,Name("dataSetCut"),Binning(bin));
	pdf->plotOn(frame_m,Name("FullPdf"),LineColor(kBlue),LineWidth(2));
	//DoubleGaussBs.plotOn(frame_m,Components(GaussBs1),LineColor(kRed+1),LineStyle(kDashed),LineWidth(1));
	//DoubleGaussBs.plotOn(frame_m,Components(GaussBs2),LineColor(kMagenta+3),LineStyle(kDashed),LineWidth(1));
	
	TCanvas* canvas = new TCanvas();
	canvas->cd();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);

        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame_m->GetYaxis()->SetRangeUser(0.01,frame_m->GetMaximum()*1.2);
        frame_m->Draw();

	stringstream ss ;
    	TString leg_min = " fs";
    	ss << std::fixed << std::setprecision(1) << min*1000. ;
    	leg_min = ss.str() + leg_min; 
	ss.str("");
    	TString leg_max = " fs";
    	ss << std::fixed << std::setprecision(1) << max*1000. ;
    	leg_max = ss.str() + leg_max; 

	ss.str("");
    	TString leg_sigma = "#sigma_{eff} = ";
    	ss << std::fixed << std::setprecision(1) << resolution_eff.getVal()*1000 ;
    	leg_sigma += ss.str(); 
	ss.str("");
    	ss << std::fixed << std::setprecision(1) << resolution_eff.getPropagatedError(*result)*1000 ;
	leg_sigma += " #pm " + ss.str() + " fs";
        
	TLegend leg(0.15,0.5,0.4,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.065);
        leg.SetTextAlign(12);
	TString label = dataType;
	if(dataType=="MC")label = "Simulation";
	leg.AddEntry((TObject*)0,"#font[22]{LHCb " + label+ "}","");
	leg.AddEntry("dataSetCut",leg_min + " < #sigma_{t} < " + leg_max,"ep");
	//leg.AddEntry("FullPdf","Fit","l");
	leg.AddEntry("FullPdf",leg_sigma,"l");
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
        
	RooPlot* frame_p = Bs_DeltaTau->frame();
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
	if(dataType=="MC")frame_p->GetXaxis()->SetTitle("#font[132]{#Deltat [ps]}");
	else frame_p->GetXaxis()->SetTitle("#font[132]{t [ps]}");        

        RooHist* pullHist  = frame_m->pullHist("dataSetCut","FullPdf");
        frame_p->addPlotable(pullHist,"BX");
        
        double maxPull = 5.0 ;
        double minPull = -5.0 ;        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(maxPull);
        graph->SetMinimum(minPull);
        graph->SetPoint(0,Bs_DeltaTau->getMin(),0);
        graph->SetPoint(1,Bs_DeltaTau->getMax(),0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(maxPull);
        graph2->SetMinimum(minPull);
        graph2->SetPoint(0,Bs_DeltaTau->getMin(),-3);
        graph2->SetPoint(1,Bs_DeltaTau->getMax(),-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(maxPull);
        graph3->SetMinimum(minPull);
        graph3->SetPoint(0,Bs_DeltaTau->getMin(),3);
        graph3->SetPoint(1,Bs_DeltaTau->getMax(),3);
        graph3->SetLineColor(kRed);
	
	pullHist->GetXaxis()->SetLabelFont( 132 );
	pullHist->GetYaxis()->SetLabelFont( 132 );
	pullHist->SetTitle("");
	
	frame_p->GetYaxis()->SetRangeUser(minPull,maxPull);
	frame_p->Draw();
	
	graph->Draw("sameL");
	graph2->Draw("sameL");
	graph3->Draw("sameL");
	
	pad2->Update();
	canvas->Update();
	
	canvas->Print("Plots/Signal"+dataType+"_bin_"+binName+".eps");
	if(updateAnaNote)canvas->Print("../../../../../TD-AnaNote/latex/figs/Resolution/Signal"+dataType+"_bin_"+binName+".pdf");

	///create a new table for Ana Note
	ofstream datafile;
	if(updateAnaNote)datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ios_base::app);
	else datafile.open("ResoTable_"+dataType+".txt", std::ios_base::app);
	datafile << std::setprecision(3) << leg_min.ReplaceAll("fs","") + " - " + leg_max.ReplaceAll("fs","") << " & "<< sig1 * 1000 << " $\\pm$ " << dsig1 * 1000 << " & " << sig2 * 1000 << " $\\pm$ " << dsig2 * 1000 << " & " << f1 << " $\\pm$ " << df1 << " & " << dilution_val << " $\\pm$ " <<  dilution_error << " & " << resolution_eff.getVal() * 1000 << " $\\pm$ " << resolution_eff.getPropagatedError(*result)* 1000 << " \\\\" << "\n";
	datafile.close();

	vector<double> resoValues;
	resoValues.push_back(resolution_eff.getVal());
	resoValues.push_back(resolution_eff.getPropagatedError(*result));
	resoValues.push_back(data->mean(Bs_TAUERR));
	return resoValues;
}

TH1D* createBinning(TString Bs_TAU_Var = "Bs_DTF_TAU"){

	NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
	NamedParameter<double> TAUERR_max("TAUERR_max", 0.15);		
        NamedParameter<int> minEventsPerBin("minEventsPerBin", 1000); 
	int dim = 1;

	double weight,dt;
	tree_res->SetBranchAddress("weight",&weight);
	tree_res->SetBranchAddress(Bs_TAU_Var+"ERR",&dt);
	
	HyperPoint Min((double)TAUERR_min);
    	HyperPoint Max((double)TAUERR_max);
    	HyperCuboid limits(Min, Max );
	HyperPointSet points( dim );

	for (int i = 0; i < tree_res->GetEntries(); i++){
		tree_res->GetEntry(i);
		
		HyperPoint point( dim );
		point.at(0)= dt;
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
	hist.draw("Plots/binning");
	hist.drawDensity("Plots/density");

	TCanvas* c = new TCanvas();
        HyperBinningPainter1D painter(&hist);
 	TH1D* h = painter.getHistogram("binning");
	//h->Draw();
	//c->Print("test.eps");
	return h;
}

void FitResoRelation(TString Bs_TAU_Var = "Bs_DTF_TAU",TString dataType = "MC"){

        NamedParameter<int> updateAnaNote("updateAnaNote", 1);

	TH1D* binning = createBinning(Bs_TAU_Var);	
	TH1D* ResoRelation = (TH1D*) binning->Clone("ResoRelation");
	
	ofstream datafile;
	if(updateAnaNote)datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ofstream::trunc);
	else datafile.open("ResoTable_"+dataType+".txt", std::ofstream::trunc);		
	datafile << "\\begin{table}[h]" << "\n";
	datafile << "\\centering" << "\n";
	datafile << " \\begin{tabular}{l || l l l | l l}" << "\n";
	datafile << "$\\sigma_{t}$ Bin [fs] & $\\sigma_{1}$ [fs] & $\\sigma_{2}$ [fs] & $f_{1}$ & D & $\\sigma_{eff}$ [fs]" << " \\\\" << "\n";
	datafile << "\\hline" << "\n";
	datafile.close();

	const int nBins = binning->GetNbinsX();
 	double x[nBins]; 
        double xerr[nBins]; 
        double xerrL[nBins]; 
        double xerrH[nBins]; 
        double y[nBins]; 
        double yerr[nBins]; 

	for(int i = 1; i <= binning->GetNbinsX(); i++){
		vector<double> reso_bin = FitTimeRes(binning->GetBinLowEdge(i),binning->GetBinLowEdge(i+1),anythingToString((int)i),Bs_TAU_Var,dataType);
		
		ResoRelation->SetBinContent(i, reso_bin[0]);
		ResoRelation->SetBinError(i, reso_bin[1]);

		x[i-1] = reso_bin[2];
		xerr[i-1] = 0.;
		xerrL[i-1]= reso_bin[2]- binning->GetBinLowEdge(i) ;
		xerrH[i-1]= binning->GetBinLowEdge(i+1) - reso_bin[2] ;
		y[i-1] = reso_bin[0];
		yerr[i-1] = reso_bin[1];
	}

        TGraphErrors *ResoRelation_g = new TGraphErrors(nBins, x,y,xerr,yerr);
	TGraphAsymmErrors *ResoRelation_ga = new TGraphAsymmErrors(nBins, x,y,xerrL,xerrH,yerr,yerr);

	if(updateAnaNote)datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ios_base::app);
	else datafile.open("ResoTable_"+dataType+".txt", std::ios_base::app);		
	datafile << "\\hline" << "\n";
	datafile << "\\end{tabular}" << "\n";
	datafile << "\\caption{Summary of the obtained parameters from the resolution fits described above.}" << "\n";
	datafile << "\\label{table:ResoParams}" << "\n";
	datafile << "\\end{table}" << "\n";
	datafile.close();

	ResoRelation_ga->SetTitle(";#sigma_{t} [ps];#sigma_{eff} [ps]");
	ResoRelation_ga->SetMinimum(0);
	ResoRelation_ga->SetMaximum(0.12);

	//define polynom for fit
	TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 0.15);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetParNames("c0","s");
	fitFunc->SetParameters(0.,1.2);
	fitFunc->SetParLimits(0,-0.5,0.5);
	fitFunc->SetParLimits(1,0.,3.);
	if(dataType=="MC")fitFunc->FixParameter(0,0.);
	//fitFunc->FixParameter(1,1.280);
	
	TF1 *fitFunc2 = new TF1("fitFunc2", "[0]+[1]*x+[2]*x*x", 0., 0.15);
	fitFunc2->SetLineColor(kGreen+3);
	fitFunc2->SetLineStyle(kDotted);
	fitFunc2->SetParNames("c0","s","s2");
	fitFunc2->SetParameters(0.,1.2,0.);
	fitFunc2->SetParLimits(0,-0.5,0.5);
	fitFunc2->SetParLimits(1,-5.,5.);
	fitFunc2->SetParLimits(2,-10.,10.);
	//fitFunc2->FixParameter(0,0.);
	//fitFunc->FixParameter(1,1.280);

	// draw polynom from DsK analysis for comparison
	TF1 *fitFunc_DsK_data = new TF1("fitFunc_DsK_data", "[0]+[1]*x ", 0., 0.15);
	fitFunc_DsK_data->SetParNames("c0_data","s_data");
	fitFunc_DsK_data->SetLineColor(kMagenta+3);
	fitFunc_DsK_data->SetLineStyle(kDotted);
	fitFunc_DsK_data->SetParameters(10.,1.2);
	fitFunc_DsK_data->FixParameter(0,0.010262);
	fitFunc_DsK_data->FixParameter(1,1.280);
	
	TF1 *fitFunc_DsK_mc = new TF1("fitFunc_DsK_mc", "[0]+[1]*x ", 0., 0.15);
	fitFunc_DsK_mc->SetParNames("c0_mc","s_mc");
	fitFunc_DsK_mc->SetLineColor(kRed);
	fitFunc_DsK_mc->SetLineStyle(kDotted);
	fitFunc_DsK_mc->SetParameters(10.,1.2);
	fitFunc_DsK_mc->FixParameter(0,0.);
	fitFunc_DsK_mc->FixParameter(1,1.201);

	TCanvas* c = new TCanvas();

        TLegend leg(0.15,0.65,0.45,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.05);
        leg.SetTextAlign(12);

	ResoRelation_g->Fit(fitFunc,"R");
	ResoRelation_g->Fit(fitFunc2,"R");

	ResoRelation_ga->Draw("AP");
	fitFunc->Draw("same");
	//fitFunc_DsK_data->Draw("same");
	if(dataType=="MC")fitFunc_DsK_mc->Draw("same");
	fitFunc2->Draw("same");
	//ResoRelation_ga->Draw("Psame");

        //leg.AddEntry((TObject*)0,"LHCb Simulation","");
        if(dataType=="MC")leg.AddEntry(ResoRelation_ga,"B_{s} #rightarrow D_{s}K#pi#pi MC","ep");
	else leg.AddEntry(ResoRelation_ga,"Prompt-D_{s} Data","ep");
        leg.AddEntry(fitFunc,"Linear Fit","l");
	leg.AddEntry(fitFunc2,"Quadratic Fit","l");
        if(dataType=="MC")leg.AddEntry(fitFunc_DsK_mc,"B_{s} #rightarrow D_{s}K MC","l");
// 	else leg.AddEntry(fitFunc2,"Quadratic Fit","l");
	leg.Draw();
	
	c->Print("Plots/ScaleFactor_"+dataType+".eps");
        if(updateAnaNote)c->Print("../../../../../TD-AnaNote/latex/figs/Resolution/ProperTimeReso_"+dataType+".pdf");
}

void fitSignalShape(TCut cut = ""){

        /// Options
        NamedParameter<int> updateAnaNote("updateAnaNote", 0);
	NamedParameter<int> sWeight("sWeight", 0);
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<double> min_MM("min_MM",1925.);
	NamedParameter<double> max_MM("max_MM",2015.);
	NamedParameter<string> InFileName("inFileNameForDsMassFit",(string)"/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_LTU.root");
	TString inFileName = (TString)((string)InFileName);
 	//NamedParameter<string> outFileName("outFileNameNorm",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
	
	/// Load file
	TFile *file = new TFile(inFileName);
	TTree* tree = (TTree*) file->Get("DecayTree");	
	tree->SetBranchStatus("weight",0);

	TFile* output;
	if(!sWeight) output = new TFile("dummy.root","RECREATE");
	else output = new TFile(inFileName.ReplaceAll("/Preselected/","/Final/"),"RECREATE");

	cut += ("Ds_MM >= " + anythingToString((double)min_MM) + " && Ds_MM <= " + anythingToString((double)max_MM)).c_str();
	TTree* out_tree = tree->CopyTree(cut);

	double sw;
    	TBranch* b_w = out_tree->Branch("weight", &sw, "weight/D");

	RooRealVar DTF_Bs_M("Ds_MM", "m(D_{s})", min_MM, max_MM,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree));
	
	/// Signal pdf
	RooRealVar mean("mean", "mean", 1968.,1960.,1980.); 
	RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
	RooRealVar gamma("gamma", "gamma", -0.5,-5,5.); 
	RooRealVar delta("delta", "delta", 0.5,-5,5.); 
	RooJohnsonSU* signal= new RooJohnsonSU("signal","signal",DTF_Bs_M, mean,sigma,gamma,delta);

	/// Bkg pdf
	RooRealVar c0("c0", "c0", .0); 
	RooRealVar c1("c1", "c1", .0,-10,10); 
	RooRealVar c2("c2", "c2", .0,-10,10); 
	RooChebychev* bkg= new RooChebychev("bkg","bkg",DTF_Bs_M, RooArgList(c0));

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.8, 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.2, 0., data->numEntries());
	RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, *bkg), RooArgList(n_sig, n_bkg));

	/// Fit
	RooFitResult* result = pdf->fitTo(*data,Save(kTRUE),NumCPU(numCPU),Extended(kTRUE));
	result->Print();

	if(sWeight){
		/// Calculate sWeights
		SPlot sPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
		int N = out_tree->GetEntries(); 
		for(int n = 0; n < N; n++){
			sw = sPlot.GetSWeight(n,"n_sig_sw");
			b_w->Fill();
		}
	}
	/// Plotting
	TCanvas* c = new TCanvas();

	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
 	data->plotOn(frame,Name("data"),Binning(100));
	pdf->plotOn(frame,Name("pdf"));
	pdf->plotOn(frame,Name("signal"),LineColor(kBlue),LineStyle(kDashed),Components("signal"));
	pdf->plotOn(frame,Name("bkg"),LineColor(kRed),LineStyle(kDashed),Components("bkg"));
	frame->Draw();
	c->Print("Ds_M.eps");

	TCanvas* canvas = new TCanvas();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);
        
        double max = 5.0 ;
        double min = -5.0 ;
        double rangeX = max-min;
        double zero = max/rangeX;
        
        TGraph* graph = new TGraph(2);
        graph->SetMaximum(max);
        graph->SetMinimum(min);
        graph->SetPoint(1,min_MM,0);
        graph->SetPoint(2,max_MM,0);
        
        TGraph* graph2 = new TGraph(2);
        graph2->SetMaximum(max);
        graph2->SetMinimum(min);
        graph2->SetPoint(1,min_MM,-3);
        graph2->SetPoint(2,max_MM,-3);
        graph2->SetLineColor(kRed);
        
        TGraph* graph3 = new TGraph(2);
        graph3->SetMaximum(max);
        graph3->SetMinimum(min);
        graph3->SetPoint(1,min_MM,3);
        graph3->SetPoint(2,max_MM,3);
        graph3->SetLineColor(kRed);
       
        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame->GetYaxis()->SetRangeUser(0.01,frame->GetMaximum()*1.);
        frame->Draw();
        
        canvas->cd();
        TPad* pad2 = new TPad("lowerPad", "lowerPad", .0, .005, 1.0, .3);
        pad2->SetBorderMode(0);
        pad2->SetBorderSize(-1);
        pad2->SetFillStyle(0);
        pad2->SetTopMargin(0.);
        pad2->SetBottomMargin(0.35);
        pad2->Draw();
        pad2->cd();
        
        RooPlot* frame_p = DTF_Bs_M.frame();
	frame_p->SetTitle("");
        frame_p->GetYaxis()->SetNdivisions(5);
        frame_p->GetYaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetLabelSize(0.12);
        frame_p->GetXaxis()->SetTitleOffset(0.75);
        frame_p->GetXaxis()->SetTitleSize(0.2);
//         frame_p->GetXaxis()->SetTitle( channelString + "[MeV/c^{2}]");
        
        RooHist* hpull  = frame->pullHist("data","pdf");
	hpull->SetTitle("");
        frame_p->addPlotable(hpull,"BX");
        frame_p->GetYaxis()->SetRangeUser(min,max);
        
        frame_p->Draw();
        graph->Draw("same");
        graph2->Draw("same");
        graph3->Draw("same");
        
        pad2->Update();
        canvas->Update();
        canvas->SaveAs("Ds_M_pull.eps");
 	if(updateAnaNote) c->Print("../../../../../TD-AnaNote/latex/figs/Resolution/Ds_M_pull.pdf");

	out_tree->Write();
	output->Close();
	return;
}

int main(int argc, char** argv){

    time_t startTime = time(0);
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
 
    /// Options
    NamedParameter<string> DataType("dataType",(string)"MC");
    TString dataType = TString((string) DataType);
    NamedParameter<string> inFileName("inFileName", (string)"/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
    NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
    NamedParameter<double> TAUERR_max("TAUERR_max", 0.12);	
    NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");		

    NamedParameter<int> fitDsMass("fitDsMass", 0);	
    NamedParameter<int> fitIntegratedResolution("fitIntegratedResolution", 0);	
    NamedParameter<int> fitResoRelation("fitResoRelation", 0);	

    if(fitDsMass)fitSignalShape("!isRejectedMultipleCandidate");

    if(fitIntegratedResolution || fitResoRelation){
	  /// Load file
  	  TFile* file = new TFile(((string)inFileName).c_str());
  	  TTree* tree = (TTree*) file->Get("DecayTree");	
	  tree->SetBranchStatus("*",0);
       	  tree->SetBranchStatus("*TAU*",1);
	  tree->SetBranchStatus("weight",1);
	  tree->SetBranchStatus("year",1);
	  tree->SetBranchStatus("*finalState",1);
  	  file_res = new TFile("dummy_res.root","RECREATE");
  	  tree_res = tree->CopyTree("");
	  file->Close();
    }

    if(fitIntegratedResolution)FitTimeRes(TAUERR_min, TAUERR_max, "all", TString((string) Bs_TAU_Var), dataType);
    if(fitResoRelation)FitResoRelation(TString((string) Bs_TAU_Var),dataType);

    //if(!file_res)file_res->Close();
  
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
