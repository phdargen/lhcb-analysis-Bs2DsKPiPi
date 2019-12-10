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

vector<double> FitTimeRes(double min, double max, TString varName = "Bs_DTF_TAUERR", string binName = "", TString Bs_TAU_Var = "Bs_DTF_TAU", TString dataType = "MC", int year = 16, int doSystematics= 0, int ngauss = 2){
	
	/// Options
        NamedParameter<string> weightVar("weightVar", (string)"weight");
	NamedParameter<int> updateAnaNote("updateAnaNote", 1);
	NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
        NamedParameter<double> TAUERR_max("TAUERR_max", 0.15);
	NamedParameter<int> useTransformedSigma("useTransformedSigma", 0);

        NamedParameter<int> fixMean("fixMean", 0);
        NamedParameter<int> fixGaussFraction("fixGaussFraction", 0);
        NamedParameter<double> gaussFraction("gaussFraction", 0.75);

	if(dataType == "Data" && year == 0) dataType.Append("_all");
	if(dataType == "Data" && year == 16) dataType.Append("_16");
	if(dataType == "Data" && year == 17) dataType.Append("_17");

	/// Load file
        RooRealVar Bs_TAU(Bs_TAU_Var, Bs_TAU_Var, -20.,20.);
        RooRealVar Bs_TAUERR(Bs_TAU_Var+"ERR", Bs_TAU_Var+"ERR", (double)TAUERR_min, (double)TAUERR_max,"ps");
	RooRealVar Bs_TRUETAU("Bs_TRUETAU", "Bs_TRUETAU", 0.,20.);
	RooRealVar weight(((string)weightVar).c_str() , ((string)weightVar).c_str(), 0.);
	RooArgList list =  RooArgList(Bs_TAU,Bs_TAUERR,weight);

        RooRealVar *var;
	if(varName != Bs_TAU_Var && varName != Bs_TAU_Var+"ERR"){ 
		var = new RooRealVar(varName, varName, 0.);
		list.add(*var);		
	}
	if(dataType == "MC")list.add(Bs_TRUETAU);
	TString cut = varName + " >= " + anythingToString(min) + " && " + varName + "<=" + anythingToString(max);

	RooDataSet* data = new RooDataSet("data","data",list,Import(*tree_res),WeightVar(((string)weightVar).c_str()),Cut(cut));

        /// Add residuals to dataset
	RooFormulaVar* Bs_DeltaTau_func;
	if(dataType=="MC")Bs_DeltaTau_func = new RooFormulaVar("Bs_DeltaTau_func","#Deltat","@0 - @1 * 1000.",RooArgList(Bs_TAU,Bs_TRUETAU));
	else Bs_DeltaTau_func = new RooFormulaVar("Bs_DeltaTau_func","t","@0",RooArgList(Bs_TAU));	
	RooRealVar* Bs_DeltaTau = (RooRealVar*) data->addColumn(*Bs_DeltaTau_func);
	Bs_DeltaTau->setUnit("ps");
	
	/// Pdf
	RooRealVar* mean1 = new RooRealVar("mean1", "mean1", 0., -0.05, 0.05);
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
	if(ngauss<3)pdf = new RooAddPdf("pdf", "pdf", RooArgList(Gauss1,Gauss2),RooArgList(*f));
	else pdf = new RooAddPdf("pdf", "pdf", RooArgList(Gauss1,Gauss2,Gauss3),RooArgList(*f,*f2));
	if(ngauss==1){
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
		fitRange_min = -0.2;
		fitRange_max = 0.2;
		//fitRange_min = -4.*data->mean(Bs_TAUERR);
		//fitRange_max = 0.5*data->mean(Bs_TAUERR);
	}
	if(dataType=="MC")Bs_DeltaTau->setRange(fitRange_min,fitRange_max);
	else Bs_DeltaTau->setRange(-0.2,0.2);

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
	TCanvas* canvas = new TCanvas();
	int bin = 50;	

	stringstream ss ;
    	TString leg_min = " fs";
    	ss << std::fixed << std::setprecision(1) << min*1000. ;
    	leg_min = ss.str() + leg_min; 
	ss.str("");
    	TString leg_max = " fs";
    	ss << std::fixed << std::setprecision(1) << max*1000. ;
    	leg_max = ss.str() + leg_max; 

	ss.str("");
/*    	TString leg_sigma = "#sigma_{eff} = ";*/
    	TString leg_sigma = "#LT#sigma_{t}#GT = ";  
  	ss << std::fixed << std::setprecision(1) << resolution_eff.getVal()*1000 ;
    	leg_sigma += ss.str(); 
	ss.str("");
    	ss << std::fixed << std::setprecision(1) << resolution_eff.getPropagatedError(*result)*1000 ;
	leg_sigma += " #pm " + ss.str() + " fs";
        
	ss.str("");
    	TString leg_mean = "#LT#mu_{t}#GT = ";  
  	ss << std::fixed << std::setprecision(2) << mean1->getVal()*1000. ;
    	leg_mean += ss.str(); 
	ss.str("");
    	ss << std::fixed << std::setprecision(2) << mean1->getError()*1000. ;
	leg_mean += " #pm " + ss.str() + " fs";

	cout << leg_sigma << endl;
	cout << leg_mean << endl;

	RooPlot* frame_m= Bs_DeltaTau->frame();
	frame_m->SetTitle("");
	frame_m->GetXaxis()->SetTitle("#font[132]{t (ps)}");
        frame_m->SetMinimum(0.0);
	data->plotOn(frame_m,Name("dataSetCut"),Binning(bin));
	pdf->plotOn(frame_m,Name("FullPdf"),LineColor(kBlue),LineWidth(3));
	//DoubleGaussBs.plotOn(frame_m,Components(GaussBs1),LineColor(kRed+1),LineStyle(kDashed),LineWidth(1));
	//DoubleGaussBs.plotOn(frame_m,Components(GaussBs2),LineColor(kMagenta+3),LineStyle(kDashed),LineWidth(1));
//         frame_m->GetYaxis()->SetRangeUser(0.01,frame_m->GetMaximum()*1.2);
        frame_m->SetMinimum(0.0);
        frame_m->Draw();	
        
	TLegend leg(0.15,0.5,0.4,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.05);
        leg.SetTextAlign(12);
	TString label = dataType;
	if(dataType=="MC")label = "Simulation";
// 	leg.AddEntry((TObject*)0,"#font[22]{LHCb " + label+ "}","");
// 	leg.AddEntry("dataSetCut",leg_min + " < #sigma_{t} < " + leg_max,"ep");
	leg.AddEntry("dataSetCut",leg_min + " < #delta_{t} < " + leg_max,"ep");
	//leg.AddEntry("FullPdf","Fit","l");
	leg.AddEntry("FullPdf",leg_sigma,"l");
 	leg.AddEntry((TObject*)0,leg_mean,"");


	leg.Draw();
	canvas->Print("Plots/Signal2"+dataType+"_bin_"+binName+".eps");
	

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
// 	frame_m->GetYaxis()->SetTitle("Yield (norm.)");


	canvas->cd();
        canvas->SetTopMargin(0.05);
        canvas->SetBottomMargin(0.05);

        TPad* pad1 = new TPad("upperPad", "upperPad", .0, .3, 1.0, 1.0);
        pad1->SetBorderMode(0);
        pad1->SetBorderSize(-1);
        pad1->SetBottomMargin(0.);
        pad1->Draw();
        pad1->cd();
        frame_m->Draw();


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
	else frame_p->GetXaxis()->SetTitle("#font[132]{t (ps)}");        

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
	if(updateAnaNote) datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ios_base::app);

	else datafile.open("ResoTable_"+dataType+".txt", std::ios_base::app);

	if(!doSystematics) datafile << std::setprecision(3) << leg_min.ReplaceAll("fs","") + " - " + leg_max.ReplaceAll("fs","") << " & "<< sig1 * 1000 << " $\\pm$ " << dsig1 * 1000 << " & " << sig2 * 1000 << " $\\pm$ " << dsig2 * 1000 << " & " << f1 << " $\\pm$ " << df1 << " & " << dilution_val << " $\\pm$ " <<  dilution_error << " & " << resolution_eff.getVal() * 1000 << " $\\pm$ " << resolution_eff.getPropagatedError(*result)* 1000 << " \\\\" << "\n";
	if(doSystematics) datafile << std::setprecision(3) << leg_min.ReplaceAll("fs","") + " - " + leg_max.ReplaceAll("fs","") << " & "<< sig1 * 1000 << " $\\pm$ " << dsig1 * 1000 << " \\\\" << "\n";

	datafile.close();

	vector<double> resoValues;
	if(!doSystematics){
		resoValues.push_back(resolution_eff.getVal());
		resoValues.push_back(resolution_eff.getPropagatedError(*result));
		if(varName == Bs_TAU_Var)resoValues.push_back(data->mean(Bs_TAU));
		else if(varName == Bs_TAU_Var+"ERR")resoValues.push_back(data->mean(Bs_TAUERR));
		else resoValues.push_back(data->mean(*var));
		resoValues.push_back(mean1->getVal());
		resoValues.push_back(mean1->getError());

		resoValues.push_back(sigma1->getVal());
		resoValues.push_back(sigma1->getError());

		resoValues.push_back(sigma2->getVal());
		resoValues.push_back(sigma2->getError());

		resoValues.push_back(f->getVal());
		resoValues.push_back(f->getError());
	}

	if(doSystematics){
		resoValues.push_back(sig1);
		resoValues.push_back(dsig1);
		if(varName == Bs_TAU_Var)resoValues.push_back(data->mean(Bs_TAU));
		else if(varName == Bs_TAU_Var+"ERR")resoValues.push_back(data->mean(Bs_TAUERR));
		else resoValues.push_back(data->mean(*var));
	}

	return resoValues;
}

TH1D* createBinning(TString var = "Bs_DTF_TAUERR"){

	NamedParameter<double> Binning_min("Binning_min", 0.);		
	NamedParameter<double> Binning_max("Binning_max", 0.1);		
        NamedParameter<int> minEventsPerBin("minEventsPerBin", 1000); 
	int dim = 1;

	double weight,dt;
	tree_res->SetBranchAddress("weight",&weight);
	tree_res->SetBranchAddress(var,&dt);
	
	HyperPoint Min((double)Binning_min);
    	HyperPoint Max((double)Binning_max);
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

void FitResoRelation(TString varName = "Bs_DTF_TAUERR", TString Bs_TAU_Var = "Bs_DTF_TAU", TString dataType = "MC", int year = 16, int doSystematics = 0, int ngauss = 2){

        NamedParameter<int> updateAnaNote("updateAnaNote", 1);

	NamedParameter<double> Binning_min("Binning_min", 0.);		
	NamedParameter<double> Binning_max("Binning_max", 0.1);	
	TH1D* binning = createBinning(varName);	
	TH1D* ResoRelation = (TH1D*) binning->Clone("ResoRelation");

	if(dataType == "Data" && year == 0) dataType.Append("_all");
	if(dataType == "Data" && year == 16) dataType.Append("_16");
	if(dataType == "Data" && year == 17) dataType.Append("_17");
	
	ofstream datafile;
	if(updateAnaNote)datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ofstream::trunc);
	else datafile.open("ResoTable_"+dataType+".txt", std::ofstream::trunc);		
	datafile << " \\begin{tabular}{c | c c c | c c}" << "\n";
	datafile << "\\hline" << "\n";
	datafile << "\\hline" << "\n";
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
        double y_mean[nBins]; 
        double yerr_mean[nBins]; 
        double y_sigma1[nBins]; 
        double yerr_sigma1[nBins]; 
        double y_sigma2[nBins]; 
        double yerr_sigma2[nBins]; 
        double y_f[nBins]; 
        double yerr_f[nBins]; 
	double scale = 1000.;

	for(int i = 1; i <= binning->GetNbinsX(); i++){
		vector<double> reso_bin = FitTimeRes(binning->GetBinLowEdge(i),binning->GetBinLowEdge(i+1),varName, anythingToString((int)i),Bs_TAU_Var,dataType, year, doSystematics, ngauss);
		
		ResoRelation->SetBinContent(i, reso_bin[0]);
		ResoRelation->SetBinError(i, reso_bin[1]);

		x[i-1] = reso_bin[2]*scale;
		xerr[i-1] = 0.;
		xerrL[i-1]= (reso_bin[2]- binning->GetBinLowEdge(i))*scale ;
		xerrH[i-1]= (binning->GetBinLowEdge(i+1) - reso_bin[2]) *scale;

		y[i-1] = reso_bin[0]*scale;
		yerr[i-1] = reso_bin[1]*scale;

		y_mean[i-1] = reso_bin[3]*scale;
		yerr_mean[i-1] = reso_bin[4]*scale;

		y_sigma1[i-1] = reso_bin[5]*scale;
		yerr_sigma1[i-1] = reso_bin[6]*scale;

		y_sigma2[i-1] = reso_bin[7]*scale;
		yerr_sigma2[i-1] = reso_bin[8]*scale;

		y_f[i-1] = reso_bin[9];
		yerr_f[i-1] = reso_bin[10];
	}

        TGraphErrors *ResoRelation_g = new TGraphErrors(nBins, x,y,xerr,yerr);
	TGraphAsymmErrors *ResoRelation_ga = new TGraphAsymmErrors(nBins, x,y,xerrL,xerrH,yerr,yerr);
	
	TGraphErrors *ResoRelation_mean = new TGraphErrors(nBins, x,y_mean,xerr,yerr_mean);
	TGraphErrors *ResoRelation_sigma1 = new TGraphErrors(nBins, x,y_sigma1,xerr,yerr_sigma1);
	TGraphErrors *ResoRelation_sigma2 = new TGraphErrors(nBins, x,y_sigma2,xerr,yerr_sigma2);
	TGraphErrors *ResoRelation_f = new TGraphErrors(nBins, x,y_f,xerr,yerr_f);

	TGraphAsymmErrors *ResoRelation_mean_a = new TGraphAsymmErrors(nBins, x,y_mean,xerrL,xerrH,yerr_mean,yerr_mean);
	TGraphAsymmErrors *ResoRelation_sigma1_a = new TGraphAsymmErrors(nBins, x,y_sigma1,xerrL,xerrH,yerr_sigma1,yerr_sigma1);
	TGraphAsymmErrors *ResoRelation_sigma2_a = new TGraphAsymmErrors(nBins, x,y_sigma2,xerrL,xerrH,yerr_sigma2,yerr_sigma2);
	TGraphAsymmErrors *ResoRelation_f_a = new TGraphAsymmErrors(nBins, x,y_f,xerrL,xerrH,yerr_f,yerr_f);

	ResoRelation_g->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_ga->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	
        ResoRelation_mean->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_sigma1->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_sigma2->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_f->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);

        ResoRelation_mean_a->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_sigma1_a->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_sigma2_a->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);
	ResoRelation_f_a->GetXaxis()->SetLimits((double)Binning_min*scale,(double)Binning_max*scale);

	if(updateAnaNote)datafile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ResoTable_"+dataType+".txt",std::ios_base::app);
	else datafile.open("ResoTable_"+dataType+".txt", std::ios_base::app);		
	datafile << "\\hline" << "\n";
	datafile << "\\hline" << "\n";
	datafile << "\\end{tabular}" << "\n";
	datafile.close();

// 	ResoRelation_ga->SetTitle(";#sigma_{t} [ps];#sigma_{eff} [ps]");
	ResoRelation_ga->SetTitle(";#delta_{t} (fs);#sigma_{t} (fs)");
	ResoRelation_ga->SetMinimum(0/1000.*scale);
	ResoRelation_ga->SetMaximum(120/1000.*scale);

	//define polynom for fit
	TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", (double)Binning_min*scale, (double)Binning_max*scale);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetParNames("c0","s");
	fitFunc->SetParameters(0.,1.2);
	fitFunc->SetParLimits(0,-50,50);
	fitFunc->SetParLimits(1,-3,3.);
	//if(dataType=="MC")fitFunc->FixParameter(0,0.);
	//fitFunc->FixParameter(1,1.280);
	
	TF1 *fitFunc2 = new TF1("fitFunc2", "[0]+[1]*x+[2]*x*x", (double)Binning_min*scale, (double)Binning_max*scale);
	fitFunc2->SetLineColor(kGreen+3);
	fitFunc2->SetLineStyle(kDotted);
	fitFunc2->SetParNames("c0","s","s2");
	fitFunc2->SetParameters(0.,1.2,0.);
	fitFunc2->SetParLimits(0,-50,50);
	fitFunc2->SetParLimits(1,-5.,5.);
	fitFunc2->SetParLimits(2,-20.,20.);
	//fitFunc2->FixParameter(0,0.);
	//fitFunc->FixParameter(1,1.280);

	// draw polynom from DsK analysis for comparison
	TF1 *fitFunc_DsK_data = new TF1("fitFunc_DsK_data", "[0]+[1]*x ", (double)Binning_min*scale, (double)Binning_max*scale);
	fitFunc_DsK_data->SetParNames("c0_data","s_data");
	fitFunc_DsK_data->SetLineColor(kMagenta+3);
	fitFunc_DsK_data->SetLineStyle(kDotted);
	fitFunc_DsK_data->SetParameters(10.,1.2);
	fitFunc_DsK_data->FixParameter(0,0.010262);
	fitFunc_DsK_data->FixParameter(1,1.280);
	
	TF1 *fitFunc_DsK_mc = new TF1("fitFunc_DsK_mc", "[0]+[1]*x ", (double)Binning_min*scale, (double)Binning_max*scale);
	fitFunc_DsK_mc->SetParNames("c0_mc","s_mc");
	fitFunc_DsK_mc->SetLineColor(kRed);
	fitFunc_DsK_mc->SetLineStyle(kDotted);
	fitFunc_DsK_mc->SetParameters(10.,1.2);
	fitFunc_DsK_mc->FixParameter(0,0.);
	fitFunc_DsK_mc->FixParameter(1,1.201);

	TF1 *fitFunc_jpsiPhi_data = new TF1("fitFunc_jpsiPhi_data", "[0]+[1]*x ", (double)Binning_min*scale, (double)Binning_max*scale);
	fitFunc_jpsiPhi_data->SetParNames("c0_data","s_data");
	fitFunc_jpsiPhi_data->SetLineColor(kRed);
	fitFunc_jpsiPhi_data->SetLineStyle(kDotted);
	fitFunc_jpsiPhi_data->SetParameters(10.,1.2);
	fitFunc_jpsiPhi_data->FixParameter(0,0.01206);
	fitFunc_jpsiPhi_data->FixParameter(1,0.8793);

	TCanvas* c = new TCanvas();
        TLegend leg(0.15,0.65,0.45,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.05);
        leg.SetTextAlign(12);

	TFitResultPtr result = ResoRelation_g->Fit(fitFunc,"RS");
	ResoRelation_g->Fit(fitFunc2,"R");

	//put fitresult in TFitResultPtr to get covariance
	TMatrixTSym<double> CovMatrix = result->GetCovarianceMatrix();
	ofstream ResoFitCov;
	ResoFitCov.open("Resolution_CovarianceMatrix_"+dataType+".txt",std::ofstream::trunc);
	ResoFitCov << "\"" << CovMatrix(0,0)/(fitFunc->GetParError(0) * fitFunc->GetParError(0)) << " " << CovMatrix(0,1)/(fitFunc->GetParError(1) * fitFunc->GetParError(0)) << " " << CovMatrix(1,0)/(fitFunc->GetParError(1) * fitFunc->GetParError(0)) << " " << CovMatrix(1,1)/(fitFunc->GetParError(1) * fitFunc->GetParError(1)) << " \"";
	ResoFitCov.close();

	ResoRelation_ga->Draw("AP");
	fitFunc->Draw("same");
	//fitFunc_DsK_data->Draw("same");
	if(dataType=="MC")fitFunc_DsK_mc->Draw("same");
// 	if(dataType!="MC" && year == 16) fitFunc_jpsiPhi_data->Draw("same");
// // 	fitFunc2->Draw("same");
	//ResoRelation_ga->Draw("Psame");

        //leg.AddEntry((TObject*)0,"LHCb Simulation","");
        if(dataType=="MC")leg.AddEntry(ResoRelation_ga,"B_{s} #rightarrow D_{s}K#pi#pi MC","ep");
	else leg.AddEntry(ResoRelation_ga,"Prompt-D_{s} Data","ep");
        leg.AddEntry(fitFunc,"Linear Fit","l");
	leg.AddEntry(fitFunc2,"Quadratic Fit","l");
        if(dataType=="MC")leg.AddEntry(fitFunc_DsK_mc,"B_{s} #rightarrow D_{s}K MC","l");
// 	if(dataType!="MC" && year == 16) leg.AddEntry(fitFunc_jpsiPhi_data,"B_{s} #rightarrow J/#psi #phi Data","l");
// 	else leg.AddEntry(fitFunc2,"Quadratic Fit","l");
// 	leg.Draw();
	
	c->Print("Plots/ScaleFactor_"+dataType+".eps");
        if(updateAnaNote) c->Print("../../../../../TD-AnaNote/latex/figs/Resolution/ScaleFactor_"+dataType+".pdf");


	fitFunc2->SetParameters(0.,-1.,0.);
	fitFunc2->SetParLimits(0,-50,50);
	fitFunc2->SetParLimits(1,-100.,100.);
	fitFunc2->SetParLimits(2,-200.,200.);
	fitFunc2->FixParameter(2,0);

	ResoRelation_mean->Fit(fitFunc,"RS");
	ResoRelation_mean_a->SetTitle(";#delta_{t} (fs);#mu_{t} (fs)");
	ResoRelation_mean_a->SetMinimum(-10/1000.*scale);
	ResoRelation_mean_a->SetMaximum(10/1000.*scale);
	ResoRelation_mean_a->Draw("AP");
	fitFunc->Draw("same");
	fitFunc2->ReleaseParameter(2);
	ResoRelation_mean->Fit(fitFunc2,"RS");
        fitFunc2->SetLineColor(kRed);
	//fitFunc2->Draw("same");
	c->Print("Plots/Bias_"+dataType+".eps");


	fitFunc2->SetParameters(0.,1.,0.);
	fitFunc2->FixParameter(2,0);

	ResoRelation_sigma1->Fit(fitFunc2,"RS");
	ResoRelation_sigma1_a->SetTitle(";#delta_{t} (fs);#sigma_{1} (fs)");
	ResoRelation_sigma1_a->SetMinimum(0/1000.*scale);
	ResoRelation_sigma1_a->SetMaximum(150/1000.*scale);
	ResoRelation_sigma1_a->Draw("AP");
	fitFunc2->Draw("same");
	c->Print("Plots/sigma1_"+dataType+".eps");

	ResoRelation_sigma2->Fit(fitFunc2,"RS");
	ResoRelation_sigma2_a->SetTitle(";#delta_{t} (fs);#sigma_{2} (fs)");
	ResoRelation_sigma2_a->SetMinimum(0/1000.*scale);
	ResoRelation_sigma2_a->SetMaximum(150/1000.*scale);
	ResoRelation_sigma2_a->Draw("AP");
	fitFunc2->Draw("same");
	c->Print("Plots/sigma2_"+dataType+".eps");

	fitFunc2->SetParameters(0.9,0.,0.);
	fitFunc2->ReleaseParameter(2);

	ResoRelation_f->Fit(fitFunc2,"RS");
	ResoRelation_f_a->SetTitle(";#delta_{t} (fs); f");
	ResoRelation_f_a->SetMinimum(0/1000.*scale);
	ResoRelation_f_a->SetMaximum(1/1000.*scale);
	ResoRelation_f_a->Draw("AP");
	fitFunc2->Draw("same");
	c->Print("Plots/f_"+dataType+".eps");


	TString dataTypeLable = dataType.ReplaceAll("_",",");
	if(updateAnaNote){
		ofstream eqfile;
		eqfile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ScaleFactor_"+dataType+".txt",std::ofstream::trunc);
		eqfile << "\\begin{equation}" << "\n";
		eqfile <<  "\\sigma_{eff}^{"+dataTypeLable+"}(\\sigma_t) = \\left( " ;
		eqfile << std::setprecision(1) << std::fixed <<  fitFunc->GetParameter(0) * 1000. << " \\pm " << fitFunc->GetParError(0) * 1000. << " \\right) \\text{fs} + \\left( ";
		eqfile << std::setprecision(3) << std::fixed  <<  fitFunc->GetParameter(1) << " \\pm " << fitFunc->GetParError(1) ;
		eqfile << " \\right) \\sigma_t" << "\n";
		eqfile << "\\label{eq:scaleFactor"+dataType+"}" << "\n";
		eqfile << "\\end{equation}" << "\n";
		eqfile.close();
	}


	if(doSystematics == 1 && ngauss == 2){
		ofstream eqfile;
		eqfile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ScaleFactor_"+dataType+"_coreGauss.txt",std::ofstream::trunc);
		eqfile << "\\begin{equation}" << "\n";
		eqfile <<  "\\sigma_{eff, "+dataTypeLable+"}^{core-Gauss}(\\sigma_t) = \\left( " ;
		eqfile << std::setprecision(1) << std::fixed <<  fitFunc->GetParameter(0) * 1000. << " \\pm " << fitFunc->GetParError(0) * 1000. << " \\right) \\text{fs} + \\left( ";
		eqfile << std::setprecision(3) << std::fixed  <<  fitFunc->GetParameter(1) << " \\pm " << fitFunc->GetParError(1) ;
		eqfile << " \\right) \\sigma_t" << "\n";
		eqfile << "\\label{eq:scaleFactor"+dataType+",core}" << "\n";
		eqfile << "\\end{equation}" << "\n";
		eqfile.close();
	}
	if(ngauss == 1){
		ofstream eqfile;
		eqfile.open("../../../../../TD-AnaNote/latex/tables/Resolution/ScaleFactor_"+dataType+"_singleGauss.txt",std::ofstream::trunc);
		eqfile << "\\begin{equation}" << "\n";
		eqfile <<  "\\sigma_{eff "+dataTypeLable+"}^{single-Gauss}(\\sigma_t) = \\left( " ;
		eqfile << std::setprecision(1) << std::fixed <<  fitFunc->GetParameter(0) * 1000. << " \\pm " << fitFunc->GetParError(0) * 1000. << " \\right) \\text{fs} + \\left( ";
		eqfile << std::setprecision(3) << std::fixed  <<  fitFunc->GetParameter(1) << " \\pm " << fitFunc->GetParError(1) ;
		eqfile << " \\right) \\sigma_t" << "\n";
		eqfile << "\\label{eq:scaleFactor"+dataType+",single}" << "\n";
		eqfile << "\\end{equation}" << "\n";
		eqfile.close();
	}
}

void fitSignalShape(TCut cut = "", int year = 16){

        /// Options
        NamedParameter<int> updateAnaNote("updateAnaNote", 0);
	NamedParameter<int> sWeight("sWeight", 0);
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<double> min_MM("min_MM",1925.);
	NamedParameter<double> max_MM("max_MM",2015.);
	NamedParameter<string> InFileName("inFileNameForDsMassFit",(string)"/auto/data/dargent/BsDsKpipi/Preselected/Data/signal_Ds2KKpi_16_LTU.root");
	TString inFileName = (TString)((string)InFileName);
 	NamedParameter<string> outFileNameForDsMassFit("outFileNameForDsMassFit",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/signal_LTU.root");
	
	/// Load file
	//TFile *file = new TFile(inFileName);
	//TTree* tree = (TTree*) file->Get("DecayTree");	
	TChain* tree = new TChain("DecayTree");
	tree->Add(inFileName);
	tree->SetBranchStatus("weight",0);

	TFile* output;
	if(!sWeight) output = new TFile("dummy.root","RECREATE");
	//else output = new TFile((inFileName.ReplaceAll("/Preselected/","/Final/")).ReplaceAll(),"RECREATE");
	else output = new TFile(((string)outFileNameForDsMassFit).c_str(),"RECREATE");

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
	RooRealVar c0("c0", "c0", .0,-1,1); 
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
        canvas->SaveAs(("Ds_M_pull_"+anythingToString(int(year))+".eps").c_str());
 	if(updateAnaNote)canvas->Print(("../../../../../TD-AnaNote/latex/figs/Resolution/Ds_M_pull_"+anythingToString(int(year))+".pdf").c_str());

	out_tree->Write();
	output->Close();
	return;
}

void reformatTuple(string input, string output, string tauName = "BeautyTime", string tauErrName = "BeautyTimeErr", double scale = 1., bool useFloat = false ){

	float t_float,dt_float;
	double t,dt,mB,Bs_TRUETAU;
	int Bs_ID,Ds_ID;
	Int_t q;
	double eta;
	double sw = 1;
	int year,run,trigger;
	double pt,Bs_eta,px,py,pz;
	double Ds_pt,Ds_eta,Ds_px,Ds_py,Ds_pz,mD;
	double pi_pt,pi_eta,pi_px,pi_py,pi_pz;
	double Ds_FDCHI2_ORIVX,Ds_PV_TAU;

	int lab1_TRUEID,lab3_TRUEID,lab4_TRUEID,lab5_TRUEID;
        double lab1_PIDK;
	int nPV;
	
	TChain* in_tree = new TChain("DecayTree");
	in_tree->Add(((string)input).c_str());
	
	TString cut = "";//"abs(lab2_TRUEID)==431 &&abs(lab1_TRUEID)==211 &&abs(lab2_MC_MOTHER_ID)!=531 &&abs(lab2_MC_MOTHER_ID)!=511 && abs(lab2_MC_MOTHER_ID)!=521 && abs(lab2_MC_MOTHER_ID)!=541 && abs(lab2_MC_MOTHER_ID)<1000 && abs(lab2_MC_GD_MOTHER_ID)!=531 &&abs(lab2_MC_GD_MOTHER_ID)!=511 &&abs(lab2_MC_GD_MOTHER_ID)!=521 &&abs(lab2_MC_GD_MOTHER_ID)!=541 &&abs(lab2_MC_GD_MOTHER_ID)<1000 &&abs(lab3_TRUEID)==321 &&abs(lab4_TRUEID)==321 &&abs(lab5_TRUEID)==211 &&abs(lab3_MC_MOTHER_ID)==431&&abs(lab4_MC_MOTHER_ID)==431 &&abs(lab5_MC_MOTHER_ID)==431 &&abs(lab3_MC_GD_MOTHER_ID)!=531&&abs(lab4_MC_GD_MOTHER_ID)!=531&&abs(lab5_MC_GD_MOTHER_ID)!=531&&abs(lab3_MC_GD_MOTHER_ID)!=511 &&abs(lab4_MC_GD_MOTHER_ID)!=511&&abs(lab5_MC_GD_MOTHER_ID)!=511&&abs(lab3_MC_GD_MOTHER_ID)!=521 &&abs(lab4_MC_GD_MOTHER_ID)!=521&&abs(lab5_MC_GD_MOTHER_ID)!=521&&abs(lab3_MC_GD_MOTHER_ID)!=541 &&abs(lab4_MC_GD_MOTHER_ID)!=541&&abs(lab5_MC_GD_MOTHER_ID)!=541&&abs(lab3_MC_GD_MOTHER_ID)<1000 &&abs(lab4_MC_GD_MOTHER_ID)<1000&&abs(lab5_MC_GD_MOTHER_ID)<1000 &&abs(lab3_MC_GD_GD_MOTHER_ID)!=53&&abs(lab4_MC_GD_GD_MOTHER_ID)!=531 &&abs(lab5_MC_GD_GD_MOTHER_ID)!=531&&abs(lab3_MC_GD_GD_MOTHER_ID)!=511 && abs(lab4_MC_GD_GD_MOTHER_ID)!=511 && abs(lab5_MC_GD_GD_MOTHER_ID)!=511&&abs(lab3_MC_GD_GD_MOTHER_ID)!=521 &&abs(lab4_MC_GD_GD_MOTHER_ID)!=521&&abs(lab5_MC_GD_GD_MOTHER_ID)!=521&&abs(lab3_MC_GD_GD_MOTHER_ID)!=541 &&abs(lab4_MC_GD_GD_MOTHER_ID)!=541&&abs(lab5_MC_GD_GD_MOTHER_ID)!=541&&abs(lab3_MC_GD_GD_MOTHER_ID)<1000 &&abs(lab4_MC_GD_GD_MOTHER_ID)<1000&&abs(lab5_MC_GD_GD_MOTHER_ID)<1000";
	TTree* tree = in_tree->CopyTree(cut);

	tree->SetBranchAddress("lab0_MM",&mB);
	tree->SetBranchAddress("lab0_PX",&px);
	tree->SetBranchAddress("lab0_PY",&py);
	tree->SetBranchAddress("lab0_PZ",&pz);

	tree->SetBranchAddress("lab2_MM",&mD);
	tree->SetBranchAddress("lab2_PX",&Ds_px);
	tree->SetBranchAddress("lab2_PY",&Ds_py);
	tree->SetBranchAddress("lab2_PZ",&Ds_pz);

	tree->SetBranchAddress("lab1_PX",&pi_px);
	tree->SetBranchAddress("lab1_PY",&pi_py);
	tree->SetBranchAddress("lab1_PZ",&pi_pz);


	tree->SetBranchAddress("lab0_TRUETAU",&Bs_TRUETAU);
	if(useFloat){
		tree->SetBranchAddress(((string)tauName).c_str(),&t_float);
		tree->SetBranchAddress(((string)tauErrName).c_str(),&dt_float);
	}
	else{
		tree->SetBranchAddress(((string)tauName).c_str(),&t);
		tree->SetBranchAddress(((string)tauErrName).c_str(),&dt);
	}
	tree->SetBranchAddress("lab2_TRUEID",&Ds_ID);
	tree->SetBranchAddress("lab0_TRUEID",&Bs_ID);
	tree->SetBranchAddress("lab1_TRUEID",&lab1_TRUEID);
	tree->SetBranchAddress("lab3_TRUEID",&lab3_TRUEID);
	tree->SetBranchAddress("lab4_TRUEID",&lab4_TRUEID);
	tree->SetBranchAddress("lab5_TRUEID",&lab5_TRUEID);
	tree->SetBranchAddress("lab1_PIDK",&lab1_PIDK);
	tree->SetBranchAddress("nPV",&nPV);
	tree->SetBranchAddress("lab2_TAU",&Ds_PV_TAU);
	tree->SetBranchAddress("lab2_FDCHI2_ORIVX",&Ds_FDCHI2_ORIVX);

	TFile* out = new TFile(((string)output).c_str(),"RECREATE");
	TTree* out_tree = new TTree("DecayTree","DecayTree");

	TBranch* br_mB = out_tree->Branch( "Bs_DTF_MM", &mB, "Bs_DTF_MM/D" );
	TBranch* br_sw= out_tree->Branch( "N_Bs_sw", &sw, "N_Bs_sw/D" );
	TBranch* br_w = out_tree->Branch( "weight", &sw, "weight/D" );

	TBranch* br_pt = out_tree->Branch( "Bs_PT", &pt, "Bs_PT/D" );
	TBranch* br_eta = out_tree->Branch( "Bs_ETA", &Bs_eta, "Bs_ETA/D" );

	TBranch* br_Dspt = out_tree->Branch( "Ds_PT", &Ds_pt, "Ds_PT/D" );
	TBranch* br_Dseta = out_tree->Branch( "Ds_ETA", &Ds_eta, "Ds_ETA/D" );

	TBranch* br_Ds_FDCHI2_ORIVX = out_tree->Branch( "Ds_FDCHI2_ORIVX", &Ds_FDCHI2_ORIVX, "Ds_FDCHI2_ORIVX/D" );
	TBranch* br_Ds_TAU = out_tree->Branch( "Ds_PV_TAU", &Ds_PV_TAU, "Ds_PV_TAU/D" );

	TBranch* br_pipt = out_tree->Branch( "pi_PT", &pi_pt, "pi_PT/D" );
	TBranch* br_pieta = out_tree->Branch( "pi_ETA", &pi_eta, "pi_ETA/D" );

	TBranch* br_t = out_tree->Branch( "Bs_BsDTF_TAU", &t, "Bs_BsDTF_TAU/D" );
	TBranch* br_dt = out_tree->Branch( "Bs_BsDTF_TAUERR", &dt, "Bs_BsDTF_TAUERR/D" );
	TBranch* br_trueTAU = out_tree->Branch( "Bs_TRUETAU", &Bs_TRUETAU, "Bs_TRUETAU/D" );

	TBranch* br_Ds_ID = out_tree->Branch("Ds_ID",&Ds_ID,"Ds_ID/I");

	TBranch* br_q_OS =out_tree->Branch("OS_Combination_DEC",&q,"OS_Combination_DEC/I");
	TBranch* br_eta_OS =  out_tree->Branch("OS_Combination_PROB",&eta,"OS_Combination_PROB/D");
	TBranch* br_q_SS = out_tree->Branch("SS_Kaon_DEC",&q,"SS_Kaon_DEC/I");
	TBranch* br_eta_SS = out_tree->Branch("SS_Kaon_PROB",&eta,"SS_Kaon_PROB/D");

	TBranch* br_run = out_tree->Branch("run",&run,"run/I");
	TBranch* br_year = out_tree->Branch( "year", &year, "year/I" );
	TBranch* br_trigger = out_tree->Branch("TriggerCat",&trigger,"TriggerCat/I");

	for(int i=0; i< tree->GetEntries(); i++)
	{	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
		tree->GetEntry(i);

// 		if(abs(lab1_TRUEID)!=211)continue;
// 		if(abs(Ds_ID)!=431)continue;
// 		if(abs(lab3_TRUEID)!=321)continue;
// 		if(abs(lab4_TRUEID)!=321)continue;
// 		if(abs(lab5_TRUEID)!=211)continue;
// 
		//if(lab1_PIDK>0)continue;
		//if(nPV>1)continue;

		if(useFloat){
			t = (double)t_float*scale;
			dt = (double)dt_float*scale;
		}
		else{
			t = t*scale;
			dt = dt*scale;		
		}	

		t = t/mB*5366.770;
	
		//if(abs(t-Bs_TRUETAU*1000.)>0.02)continue;

		run = 2;
		year = 17;
		trigger = 0;

		if(Bs_ID>0)q=1;
		else q = -1;
		eta = 0.;

		TLorentzVector p;
		p.SetXYZM(px,py,pz,mB);

		pt = p.Pt();
		Bs_eta = p.PseudoRapidity();

		TLorentzVector p_Ds;
		p_Ds.SetXYZM(Ds_px,Ds_py,Ds_pz,mD);

		Ds_pt = p_Ds.Pt();
		Ds_eta = p_Ds.PseudoRapidity();

		TLorentzVector p_pi;
		p_pi.SetXYZM(pi_px,pi_py,pi_pz, 139.57);

		pi_pt = p_pi.Pt();
		pi_eta = p_pi.PseudoRapidity();

		out_tree->Fill();
	}

	out_tree->Write();
	out->Write();
	out->Close();
}

int main(int argc, char** argv){

//        reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/0micron.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/norm_0micron_scaled.root","BeautyTime","BeautyTimeErr",1.,false);
//        reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/4micron.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/norm_4micron_scaled.root","BeautyTime","BeautyTimeErr",1.,false);

  //   reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/PromptDsMCLRMisaligned_0micron_m*_All.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/prompt_0micron.root","lab0_LifetimeFit_ctau0","lab0_LifetimeFit_ctauErr0",3.33564095,true);
  //   reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/PromptDsMCLRMisaligned_4micron_m*_All.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/prompt_4micron.root","lab0_LifetimeFit_ctau0","lab0_LifetimeFit_ctauErr0",3.33564095,true);
// 
  //   reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/PromptDsMCLRMisaligned_0micron_m*_All.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/prompt_0micron_TAU.root","lab0_TAU","lab0_TAUERR",1000);
//     reformatTuple("/auto/data/dargent/BsDsKpipi/decayTimeBias/PromptDsMCLRMisaligned_4micron_m*_All.root","/auto/data/dargent/BsDsKpipi/decayTimeBias/prompt_4micron_TAU.root","lab0_TAU","lab0_TAUERR",1000);

//       return 0;


    time_t startTime = time(0);
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
 
    /// Options
    NamedParameter<int> Year("Year", 16);
    int year = Year;
    NamedParameter<string> DataType("dataType",(string)"MC");
    TString dataType = TString((string) DataType);
    NamedParameter<string> inFileName("inFileName", (string)"/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
    NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
    NamedParameter<double> TAUERR_max("TAUERR_max", 0.12);	
    NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");	
    NamedParameter<int> DoSystematics("DoSystematics", 0);
    int doSystematics = DoSystematics;
    NamedParameter<int> nGauss("nGauss", 2);
    int ngauss = nGauss;

    NamedParameter<int> fitDsMass("fitDsMass", 0);	
    NamedParameter<int> fitIntegratedResolution("fitIntegratedResolution", 0);	
    NamedParameter<int> fitResoRelation("fitResoRelation", 0);	
    NamedParameter<string> cut("cut",(string)"");

    if(fitDsMass)fitSignalShape("!isRejectedMultipleCandidate", year);
    //if(fitDsMass)fitSignalShape("", year);


    if(fitIntegratedResolution || fitResoRelation){
	  /// Load file
  	  TFile* file = new TFile(((string)inFileName).c_str());
  	  TTree* tree = (TTree*) file->Get("DecayTree");	
	  tree->SetBranchStatus("*",0);
       	  tree->SetBranchStatus("*TAU*",1);
	  tree->SetBranchStatus("weight*",1);
	  tree->SetBranchStatus("year",1);
	  tree->SetBranchStatus("run",1);
	  tree->SetBranchStatus("Ds_FDCHI2_ORIVX",1);
	  tree->SetBranchStatus("Bs_PT",1);
	  tree->SetBranchStatus("Bs_DTF_MM",1);
	  tree->SetBranchStatus("*finalState*",1);
	  tree->SetBranchStatus("*bkg*",1);


  	  file_res = new TFile("dummy_res.root","RECREATE");
	 // if(dataType == "Data"){
           //     if(year == 0)tree_res = tree->CopyTree("");
           //     else tree_res = tree->CopyTree(("year == "+anythingToString((int)year)).c_str());
         // }
//  	  else 
	  tree_res = tree->CopyTree(((string)cut).c_str());
	  file->Close();
    }
    if(fitIntegratedResolution)FitTimeRes(TAUERR_min, TAUERR_max, TString((string) Bs_TAU_Var)+"ERR", "all", TString((string) Bs_TAU_Var), dataType, year, doSystematics, ngauss);
    if(fitResoRelation)FitResoRelation(TString((string) Bs_TAU_Var)+"ERR",TString((string) Bs_TAU_Var),dataType, year, doSystematics, ngauss);
   //if(fitResoRelation)FitResoRelation(TString("Bs_DTF_MM"),TString((string) Bs_TAU_Var),dataType, year, doSystematics, ngauss);

    if(!file_res)file_res->Close();
  
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
