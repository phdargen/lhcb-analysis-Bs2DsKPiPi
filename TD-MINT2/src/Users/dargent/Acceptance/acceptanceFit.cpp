// Fits the time acceptance
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
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

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;


void fitSplineAcc(){

// Read Dataset
TFile* file;
TTree* tree;

file= new TFile("/auto/data/kecke/B2DPiPiPi/Bs2Dspipipi_Ds2KKpi_fullSelection_combined.root");
tree = (TTree*) file->Get("DecayTree");



//Define RooRealVar for observables
RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar Bs_ct("Bs_ct", "Bs_ct", 0.4, 10., "ps");
RooRealVar Bs_cterr("Bs_cterr", "Bs_cterr", 0.0001, 1.,"ps");
RooRealVar N_Bs_sw("N_Bs_sw", "N_Bs_sw", -0.3, 1.3);

RooArgSet observables(Bs_ct, Bs_cterr, N_Bs_sw);

RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree), WeightVar(N_Bs_sw.GetName()));



///SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
//SPLINE KNOTS

double myknots [6] = { 0.5,1., 1.5, 2., 3., 9.5};
cout << "Bs_ct_min:   " << Bs_ct.getMin() << endl;
cout << "Bs_ct_max:   " << Bs_ct.getMax() << endl;


std::vector<double> myBinning;

myBinning.push_back(myknots[0]);
myBinning.push_back(myknots[1]);
myBinning.push_back(myknots[2]);
myBinning.push_back(myknots[3]);
myBinning.push_back(myknots[4]);
myBinning.push_back(myknots[5]);

double values [6] =  {3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.7757e+00};
int numKnots = 6;


//SPLINE COEFFICIENTS
RooRealVar coeff_0("coeff_0", "coeff_0", values[0], 0.0, 10.0);
RooRealVar coeff_1("coeff_1", "coeff_1", values[1], 0.0, 10.0);
RooRealVar coeff_2("coeff_2", "coeff_2", values[2], 0.0, 10.0);
RooRealVar coeff_3("coeff_3", "coeff_3", values[3], 0.0, 10.0);
RooRealVar coeff_4("coeff_4", "coeff_4", values[4], 0.0, 10.0);
RooRealVar coeff_5("coeff_5", "coeff_5", values[5], 0.0, 10.0);

RooRealVar coeff_6("coeff_6", "coeff_6", 1.0);


//#############Implement GetCoefffomBinning directly###################
//Do linear extrapolation for last coeff
/*
double deltaY = values[6] - values[5];
double deltaX= myknots[5] - myknots[4];
double gradient = deltaY / deltaX;
double lastcoeff = 1 + gradient * (Bs_ct.getMax() - myknots[5]);
*/

RooRealVar knot_5("knot_5","knot_5",myknots[5]);
RooRealVar knot_4("knot_4","knot_4",myknots[4]);
RooRealVar Bs_ct_max("Bs_ct_max","Bs_ct_max", Bs_ct.getMax());

RooGenericPdf coeff_7("coeff_7","coeff_7", "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(coeff_6, coeff_5, knot_5, knot_4, Bs_ct_max) );


//RooArgList* tacc_list;
RooArgList tacc_list(coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_6, coeff_7);

double tauH = 1.661000;
double tauL = 1.405000;
double gammaH = 1.0/tauH;
double gammaL = 1.0/tauL;
double gamma = (gammaH + gammaL) / 2.0;
double tau = 1.512;
double dgamma = -0.082; //gammaH - gammaL
double DeltaMs = 17.768;


//CUBIC SPLINE FUNCTION 
RooCubicSplineFun* spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_ct, myBinning, tacc_list);


RooRealVar trm_mean( "trm_mean" , "Gaussian resolution model mean", 0.0, "ps" );
//RooRealVar trm_sigma( "trm_sigma" , "Gaussian resolution model sigma" , 0.040 , "ps");
RooRealVar trm_scale( "trm_scale", "Gaussian resolution model scale factor", 1.20);
RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_ct, *spl, trm_mean, Bs_cterr, trm_mean, trm_scale );


RooBDecay* timePDF = new RooBDecay("Bdecay", "Bdecay", Bs_ct, RooRealConstant::value(tau),
                     RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                     RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(DeltaMs),
                     trm, RooBDecay::SingleSided);


///Fit and Print
//Fit
RooFitResult *myfitresult;
myfitresult = timePDF->fitTo(*dataset, Save(1), Optimize(2), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE));
myfitresult->Print("v");


//Print
double range_dw = Bs_ct.getMin();
double range_up = Bs_ct.getMax();

TLegend* legend = new TLegend( 0.62, 0.70, 0.88, 0.88 );

legend->SetTextSize(0.05);
legend->SetTextFont(12);
legend->SetFillColor(4000);
legend->SetShadowColor(0);
legend->SetBorderSize(0);
legend->SetTextFont(132);

TLine* l1 = new TLine();
l1->SetLineColor(kBlue+3);
l1->SetLineWidth(4);
l1->SetLineStyle(kSolid);
legend->AddEntry(l1, "decay", "L");

TLine* l2 = new TLine();
l2->SetLineColor(kRed);
l2->SetLineWidth(4);
l2->SetLineStyle(kSolid);
legend->AddEntry(l2, "acceptance", "L");

RooPlot* frame_m = Bs_ct.frame();
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

int bin = 150;
dataset->plotOn(frame_m, Binning(bin), Name("dataSetCut"));
timePDF->plotOn(frame_m, LineColor(kBlue+3),  Name("FullPdf"));
spl->plotOn(frame_m, LineColor(kRed), Normalization(50, RooAbsReal::Relative));

TLatex* lhcbtext = new TLatex();
lhcbtext->SetTextFont(132);
lhcbtext->SetTextColor(1);
lhcbtext->SetTextSize(0.07);
lhcbtext->SetTextAlign(12);

TCanvas* canvas = new TCanvas("canvas", "canvas", 1200, 800);
canvas->cd();
canvas->SetLeftMargin(0.01);
canvas->SetRightMargin(0.01);
canvas->SetTopMargin(0.05);
canvas->SetBottomMargin(0.05);
TPad* pad1 = new TPad("upperPad", "upperPad", .050, .22, 1.0, 1.0);
pad1->SetBorderMode(0);
pad1->SetBorderSize(-1);
pad1->SetFillStyle(0);
pad1->SetTickx(0);
pad1->SetLeftMargin(0.115);
pad1->SetRightMargin(0.05);
pad1->Draw();
pad1->cd();
frame_m->GetYaxis()->SetRangeUser(0.1,frame_m->GetMaximum()*1.0);
frame_m->Draw();

lhcbtext->DrawTextNDC(0.70,0.85,"LHCb Run I data");
lhcbtext->DrawTextNDC(0.70,0.75,"preliminary");

canvas->cd();
TPad* pad2 = new TPad("lowerPad", "lowerPad", .050, .005, 1.0, .3275);
pad2->SetBorderMode(0);
pad2->SetBorderSize(-1);
pad2->SetFillStyle(0);
pad2->SetBottomMargin(0.35);
pad2->SetLeftMargin(0.115);
pad2->SetRightMargin(0.05);
pad2->SetTickx(0);
pad2->Draw();
pad2->SetLogy(0);
pad2->cd();

frame_m->Print("v");
RooPlot* frame_p = Bs_ct.frame(Title("pull_frame"));
frame_p->Print("v");
frame_p->SetTitle("");
frame_p->GetYaxis()->SetTitle("");
frame_p->GetYaxis()->SetTitleSize(0.09);
frame_p->GetYaxis()->SetTitleOffset(0.26);
frame_p->GetYaxis()->SetTitleFont(62);
frame_p->GetYaxis()->SetNdivisions(106);
frame_p->GetYaxis()->SetLabelSize(0.12);
frame_p->GetYaxis()->SetLabelOffset(0.006);
frame_p->GetXaxis()->SetTitleSize(0.15);
frame_p->GetXaxis()->SetTitleFont(132);
frame_p->GetXaxis()->SetTitleOffset(0.85);
frame_p->GetXaxis()->SetNdivisions(515);
frame_p->GetYaxis()->SetNdivisions(5);
frame_p->GetXaxis()->SetLabelSize(0.12);
frame_p->GetXaxis()->SetLabelFont( 132 );
frame_p->GetYaxis()->SetLabelFont( 132 );
frame_p->GetXaxis()->SetTitle("#font[132]{t(B_{s}) [ps]}");

TString* obsTS = new TString(Bs_ct.GetName());
TString* pullnameTS = new TString("FullPdf");
TString* pullname2TS = new TString("dataSetCut");
RooHist* pullHist  = frame_m->pullHist(pullname2TS->Data(),pullnameTS->Data());
frame_p->addPlotable(pullHist,"P");

double chi2 = frame_m->chiSquare();
double chi22 = frame_m->chiSquare(pullnameTS->Data(),pullname2TS->Data());

TAxis* axisX = pullHist->GetXaxis();
RooBinning* Bin = new RooBinning(range_dw,range_up,"P");
Bin->addUniform(bin, range_dw, range_up);
axisX->Set(Bin->numBins(), Bin->array());

TAxis* axisY = pullHist->GetYaxis();
double max = 5.0 ;
double min = -5.0 ;
axisY->SetLabelSize(0.12);
axisY->SetNdivisions(5);
axisX->SetLabelSize(0.12);

double rangeX = max-min;
double zero = max/rangeX;

TGraph* graph = new TGraph(2);
graph->SetMaximum(max);
graph->SetMinimum(min);
graph->SetPoint(1,range_dw,0);
graph->SetPoint(2,range_up,0);

TGraph* graph2 = new TGraph(2);
graph2->SetMaximum(max);
graph2->SetMinimum(min);
graph2->SetPoint(1,range_dw,-3);
graph2->SetPoint(2,range_up,-3);
graph2->SetLineColor(kRed);

TGraph* graph3 = new TGraph(2);
graph3->SetMaximum(max);
graph3->SetMinimum(min);
graph3->SetPoint(1,range_dw,3);
graph3->SetPoint(2,range_up,3);
graph3->SetLineColor(kRed);

pullHist->GetXaxis()->SetLabelFont( 132 );
pullHist->GetYaxis()->SetLabelFont( 132 );
pullHist->SetTitle("");

frame_p->GetYaxis()->SetRangeUser(-5.0,5.0);
frame_p->Draw();

graph->Draw("same");

pad2->Update();
canvas->Update();

string epsName = "Plot/";
string PlotName = "timeAccFit.eps";
string epsFullName = epsName.append(PlotName);

canvas->SaveAs(epsFullName.c_str());

file->Close();

}


int main(int argc, char** argv){
    
    time_t startTime = time(0);
    //gROOT->ProcessLine(".x ../lhcbStyle.C");


    fitSplineAcc();
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
