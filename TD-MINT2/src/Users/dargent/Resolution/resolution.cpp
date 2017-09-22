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
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooNDKeysPdf.h"
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

vector<double> FitTimeRes(double min, double max, string binName = ""){
	
	NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 1);
	NamedParameter<int> updateTable("updateTable", 1);
	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");

	TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
	TTree* tree= (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*finalState",1);
	
        RooRealVar Bs_TAU(((string)Bs_TAU_Var).c_str(), ((string)Bs_TAU_Var).c_str(), 0.,20.);
        RooRealVar Bs_TAUERR(((string)Bs_TAU_Var+"ERR").c_str(), ((string)Bs_TAU_Var+"ERR").c_str(), min, max,"ps");
	RooRealVar Bs_TRUETAU("Bs_TRUETAU", "Bs_TRUETAU", 0.,20.);
	RooRealVar weight("weight" , "weight", 0.);
        RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0.);
        RooRealVar year("year", "year", 0.);	

	RooArgList list =  RooArgList(Bs_TAU,Bs_TAUERR,Bs_TRUETAU,weight,Ds_finalState,year);
	RooDataSet* data = new RooDataSet("data","data",list,Import(*tree),WeightVar("weight"));

	RooFormulaVar Bs_DeltaTau_func("Bs_DeltaTau_func","#Deltat","@0 - @1 * 1000.",RooArgList(Bs_TAU,Bs_TRUETAU));
	RooRealVar* Bs_DeltaTau = (RooRealVar*) data->addColumn(Bs_DeltaTau_func);
	Bs_DeltaTau->setRange(-0.2,0.2);

	RooRealVar meanBs1("meanBs1", "B_{s} #mu", 0., -0.2, 0.2);
	RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 0.020,0.,0.1);
	RooRealVar scale("scale", "scale", 2.,1.,10.);
	RooFormulaVar sigmaBs2("sigmaBs2", "@0*@1", RooArgList(scale,sigmaBs1));
	//RooRealVar sigmaBs2("sigmaBs2", "B_{s} #sigma_{2}", 0.045,0.,0.2);
	RooRealVar sigmaBs3("sigmaBs3", "B_{s} #sigma_{3}", 0.040,0.,0.2);
	RooGaussian GaussBs1("GaussBs1", "GaussBs1", *Bs_DeltaTau, meanBs1, sigmaBs1);
	RooGaussian GaussBs2("GaussBs2", "GaussBs2", *Bs_DeltaTau, meanBs1, sigmaBs2);
	RooGaussian GaussBs3("GaussBs3", "GaussBs3", *Bs_DeltaTau, meanBs1, sigmaBs3);
	RooGaussian GaussBs("GaussBs", "GaussBs", *Bs_DeltaTau, meanBs1, sigmaBs1);
	RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.8, 0., 1.);
	RooRealVar f_GaussBs2("f_GaussBs2" , "2f__{B_{s}}", 0.5, 0., 1.);
	RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
	RooAddPdf TripleGaussBs("TripleGaussBs", "TripleGaussBs", RooArgList(GaussBs1,GaussBs2,GaussBs3),RooArgList(f_GaussBs,f_GaussBs2));
	
	RooFitResult *result = DoubleGaussBs.fitTo(*data,Save(kTRUE),SumW2Error(kTRUE),Extended(kFALSE),NumCPU(3),Range(-0.2,0.2));
	cout << "result is --------------- "<<endl;
	result->Print();
	
	int bin = 75;
	
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
	frame_m->GetYaxis()->SetTitleSize( 0.06 );
	frame_m->GetYaxis()->SetNdivisions(512);
	frame_m->GetXaxis()->SetTitleOffset( 1.00 );
	frame_m->GetYaxis()->SetTitleOffset( 1.00 );
	
	data->plotOn(frame_m,Name("dataSetCut"),Binning(bin));
	DoubleGaussBs.plotOn(frame_m,Name("FullPdf"),LineColor(kBlue),LineWidth(2));
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
        frame_m->GetYaxis()->SetRangeUser(0.01,frame_m->GetMaximum()*1.);
        frame_m->Draw();

	stringstream ss ;
    	TString leg_min = " fs";
    	ss << std::fixed << std::setprecision(1) << min*1000. ;
    	leg_min = ss.str() + leg_min; 
	ss.str("");
    	TString leg_max = " fs";
    	ss << std::fixed << std::setprecision(1) << max*1000. ;
    	leg_max = ss.str() + leg_max; 

        TLegend leg(0.6,0.5,0.9,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.06);
        leg.SetTextAlign(12);
	leg.AddEntry((TObject*)0,"#font[22]{LHCb Simulation}","");
	leg.AddEntry("dataSetCut",leg_min + " < #sigma_{t} < " + leg_max,"ep");
	leg.AddEntry("FullPdf","Fit","l");
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
	frame_p->GetXaxis()->SetTitle("#font[132]{#Deltat(B_{s}) [ps]}");
        
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
	
	canvas->Print(("Plots/SignalMC_bin_"+binName+".eps").c_str());
	if(updateAnaNotePlots)canvas->Print(("../../../../../TD-AnaNote/latex/figs/Resolution/+SignalMC_bin_"+binName+".pdf").c_str());

	double f1 = f_GaussBs.getVal();
	double df1 = f_GaussBs.getError();
	double sig1 = sigmaBs1.getVal();
	double dsig1 = sigmaBs1.getError();
	double sig2 = sigmaBs2.getVal();
	double dsig2 = sigmaBs2.getPropagatedError(*result);
	double dms = 17.757;

	RooFormulaVar dilution("dilution","@0 * exp(-@1*@1*@2*@2/2.) + (1. - @0) * exp(-@3*@3*@2*@2/2.)",RooArgList(f_GaussBs,sigmaBs1,RooRealConstant::value(dms),sigmaBs2));

	double dilution_val = dilution.getVal();
	double dilution_error = dilution.getPropagatedError(*result);

	RooFormulaVar resolution_eff("resolution_eff","sqrt(-2./@0/@0*log(@1))",RooArgList(RooRealConstant::value(dms),dilution)); 
	//double resolution_eff = sqrt(-2./pow(dms,2)*log(dilution_val));
	//double resolution_eff_error = -1./(dms*dilution_val*sqrt(log(1./pow(dilution_val,2))))*dilution_error;
	//double resolution_eff_error = ((2/(dms*dms))/(2*dilution_val*TMath::Sqrt((-2/(dms*dms))*log(dilution_val)))) * dilution_error;

	cout << "Measured resolution from dilution:   " << resolution_eff.getVal()*1000 << " +/- " << resolution_eff.getPropagatedError(*result)*1000 <<" fs" << endl;

	vector<double> resoValues;
	resoValues.push_back(resolution_eff.getVal());
	resoValues.push_back(resolution_eff.getPropagatedError(*result));
	resoValues.push_back(data->mean(Bs_TAUERR));
	
	file->Close();
	///create a new table for Ana Note
	if(updateTable){
		ofstream datafile;
		if(updateAnaNotePlots)datafile.open("../../../../../TD-AnaNote/latex/tables/ResoTable.txt", std::ios_base::app);
		else datafile.open("ResoTable.txt", std::ios_base::app);
		datafile << std::setprecision(3) << leg_min.ReplaceAll("fs","") + " - " + leg_max.ReplaceAll("fs","") << " & "<< sig1 * 1000 << " $\\pm$ " << dsig1 * 1000 << " & " << sig2 * 1000 << " $\\pm$ " << dsig2 * 1000 << " & " << f1 << " $\\pm$ " << df1 << " & " << dilution_val << " $\\pm$ " <<  dilution_error << " & " << resolution_eff.getVal() * 1000 << " $\\pm$ " << resolution_eff.getPropagatedError(*result)* 1000 << " \\\\" << "\n";
		datafile.close();
	}

	return resoValues;
}

TH1D* createBinning(){

	NamedParameter<string> Bs_TAU_Var("Bs_TAU_Var",(string)"Bs_TAU");		
	NamedParameter<double> TAUERR_min("TAUERR_min", 0.);		
	NamedParameter<double> TAUERR_max("TAUERR_max", 0.15);		
        NamedParameter<int> minEventsPerBin("minEventsPerBin", 1000); 
	int dim = 1;

	TFile* file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
	TTree* tree= (TTree*) file->Get("DecayTree");
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus(((string)Bs_TAU_Var+"ERR").c_str(),1);
	double weight,dt;
	tree->SetBranchAddress("weight",&weight);
	tree->SetBranchAddress(((string)Bs_TAU_Var+"ERR").c_str(),&dt);
	
	HyperPoint Min((double)TAUERR_min);
    	HyperPoint Max((double)TAUERR_max);
    	HyperCuboid limits(Min, Max );
	HyperPointSet points( dim );

	for (int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		
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

void FitResoRelation(){

        NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 1);
	NamedParameter<int> updateTable("updateTable", 1);

	TH1D* binning = createBinning();	
	TH1D* ResoRelation = (TH1D*) binning->Clone("ResoRelation");
	
	ofstream datafile;
	if(updateTable){
		if(updateAnaNotePlots)datafile.open("../../../../../TD-AnaNote/latex/tables/ResoTable.txt", std::ios_base::app);
		else datafile.open("ResoTable.txt", std::ios_base::app);		datafile << "\\begin{table}[h]" << "\n";
		datafile << "\\centering" << "\n";
		datafile << " \\begin{tabular}{l || l l l | l l}" << "\n";
		datafile << "$\\sigma_{t}$ Bin [fs] & $\\sigma_{1}$ [fs] & $\\sigma_{2}$ [fs] & $f_{1}$ & D & $\\sigma_{eff}$ [fs]" << " \\\\" << "\n";
		datafile << "\\hline" << "\n";
		datafile.close();
	}

	const int nBins = binning->GetNbinsX();
 	double x[nBins]; 
        double xerr[nBins]; 
        double xerrL[nBins]; 
        double xerrH[nBins]; 
        double y[nBins]; 
        double yerr[nBins]; 

	for(int i = 1; i <= binning->GetNbinsX(); i++){
		vector<double> reso_bin = FitTimeRes(binning->GetBinLowEdge(i),binning->GetBinLowEdge(i+1),anythingToString((int)i));
		
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

	if(updateTable){
		if(updateAnaNotePlots)datafile.open("../../../../../TD-AnaNote/latex/tables/ResoTable.txt", std::ios_base::app);
		else datafile.open("ResoTable.txt", std::ios_base::app);		datafile << "\\hline" << "\n";
		datafile << "\\end{tabular}" << "\n";
		datafile << "\\caption{Summary of the obtained parameters from the resolution fits described above.}" << "\n";
		datafile << "\\label{table:ResoParams}" << "\n";
		datafile << "\\end{table}" << "\n";
		datafile.close();
	}

	ResoRelation_ga->SetTitle(";#sigma_{t} [ps];#sigma_{eff} [ps]");
	ResoRelation_ga->SetMinimum(0);
	ResoRelation_ga->SetMaximum(0.12);

	//define polynom for fit
	TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 0.15);
	fitFunc->SetLineColor(kBlue);
	fitFunc->SetParNames("c0","s");
	fitFunc->SetParameters(0.,1.2);
	fitFunc->SetParLimits(0,0.,0.1);
	fitFunc->SetParLimits(1,0.,3.);
	fitFunc->FixParameter(0,0.);
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

        TLegend leg(0.6,0.65,0.9,0.9,"");
        leg.SetLineStyle(0);
        leg.SetLineColor(0);
        leg.SetFillColor(0);
        leg.SetTextFont(132);
        leg.SetTextColor(1);
        leg.SetTextSize(0.05);
        leg.SetTextAlign(12);

	ResoRelation_g->Fit(fitFunc,"R");
	ResoRelation_ga->Draw("AP");
	fitFunc->Draw("same");
	//fitFunc_DsK_data->Draw("same");
	fitFunc_DsK_mc->Draw("same");
	ResoRelation_ga->Draw("Psame");

        //leg.AddEntry((TObject*)0,"LHCb Simulation","");
        leg.AddEntry(ResoRelation_ga,"B_{s} #rightarrow D_{s}K#pi#pi MC","ep");
        leg.AddEntry(fitFunc,"Fit","l");
        leg.AddEntry(fitFunc_DsK_mc,"B_{s} #rightarrow D_{s}K MC","l");
	leg.Draw();
	
	c->Print("Plots/ProperTimeReso_MC.eps");
        if(updateAnaNotePlots)c->Print("../../../../../TD-AnaNote/latex/figs/Resolution/ProperTimeReso_MC.pdf");
}

int main(int argc, char** argv){

    time_t startTime = time(0);
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();

    FitResoRelation();

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
