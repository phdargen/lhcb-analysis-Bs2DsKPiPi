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
#include "Mint/RooSplineProduct.h"

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;


vector<double> fitSplineAcc(string CutString){

// Read Dataset
TFile* file;
TTree* tree;

file= new TFile("/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
tree = (TTree*) file->Get("DecayTree");

//Define RooRealVar for observables
RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar Bs_TAU("Bs_TAU", "Bs_TAU", 0.4, 10., "ps");
RooRealVar Bs_TAUERR("Bs_TAUERR", "Bs_TAUERR", 0.0001, 1.,"ps");
RooRealVar N_Bs_sw("N_Bs_sw", "N_Bs_sw", -0.3, 1.3);
RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0, 3);
RooRealVar year("year", "year", 11, 16);

RooArgSet observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, N_Bs_sw);

RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree), WeightVar(N_Bs_sw.GetName()), Cut(CutString.c_str()));


///SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
//SPLINE KNOTS

double myknots [7] = { 0.5, 1., 1.5, 2., 3., 6., 9.5};
cout << "Bs_TAU_min:   " << Bs_TAU.getMin() << endl;
cout << "Bs_TAU_max:   " << Bs_TAU.getMax() << endl;


std::vector<double> myBinning;

myBinning.push_back(myknots[0]);
myBinning.push_back(myknots[1]);
myBinning.push_back(myknots[2]);
myBinning.push_back(myknots[3]);
myBinning.push_back(myknots[4]);
myBinning.push_back(myknots[5]);
myBinning.push_back(myknots[6]);

double values [6] =  {3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.7757e+00};
int numKnots = 6;

RooRealVar coeff_55("coeff_55", "coeff_55", 1., 0.0, 10.0);

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
double lastcoeff = 1 + gradient * (Bs_TAU.getMax() - myknots[5]);
*/

RooRealVar knot_4("knot_4","knot_4",myknots[4]);
RooRealVar knot_5("knot_5","knot_5",myknots[5]);
RooRealVar knot_6("knot_6","knot_6",myknots[6]);
RooRealVar Bs_TAU_max("Bs_TAU_max","Bs_TAU_max", Bs_TAU.getMax());

RooGenericPdf coeff_7("coeff_7","coeff_7", "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(coeff_6, coeff_55, knot_6, knot_5, Bs_TAU_max) );


//RooArgList* tacc_list;
RooArgList tacc_list(coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_55, coeff_6, coeff_7);

double tauH = 1.661000;
double tauL = 1.405000;
double gammaH = 1.0/tauH;
double gammaL = 1.0/tauL;
double gamma = (gammaH + gammaL) / 2.0;
double tau = 1.512;
double dgamma = -0.082; //gammaH - gammaL
double DeltaMs = 17.768;


//CUBIC SPLINE FUNCTION 
RooCubicSplineFun* spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_TAU, myBinning, tacc_list);


RooRealVar trm_mean( "trm_mean" , "Gaussian resolution model mean", 0.0, "ps" );
//RooRealVar trm_sigma( "trm_sigma" , "Gaussian resolution model sigma" , 0.040 , "ps");
RooRealVar trm_scale( "trm_scale", "Gaussian resolution model scale factor", 1.20);
RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *spl, trm_mean, Bs_TAUERR, trm_mean, trm_scale );


RooBDecay* timePDF = new RooBDecay("Bdecay", "Bdecay", Bs_TAU, RooRealConstant::value(tau),
                     RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                     RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(DeltaMs),
                     trm, RooBDecay::SingleSided);


///Fit and Print
//Fit
RooFitResult *myfitresult;
myfitresult = timePDF->fitTo(*dataset, Save(1), Optimize(2), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE));
myfitresult->Print("v");

//put coefficients into vector
std::vector<double> myCoeffs;
//values
myCoeffs.push_back(coeff_0.getValV());
myCoeffs.push_back(coeff_1.getValV());
myCoeffs.push_back(coeff_2.getValV());
myCoeffs.push_back(coeff_3.getValV());
myCoeffs.push_back(coeff_4.getValV());
myCoeffs.push_back(coeff_5.getValV());
myCoeffs.push_back(coeff_55.getValV());
myCoeffs.push_back(coeff_6.getVal());
myCoeffs.push_back(coeff_7.getVal());
//errors
myCoeffs.push_back(coeff_0.getError());
myCoeffs.push_back(coeff_1.getError());
myCoeffs.push_back(coeff_2.getError());
myCoeffs.push_back(coeff_3.getError());
myCoeffs.push_back(coeff_4.getError());
myCoeffs.push_back(coeff_5.getError());
myCoeffs.push_back(coeff_55.getError());

//Print
double range_dw = Bs_TAU.getMin();
double range_up = Bs_TAU.getMax();

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
legend->AddEntry(l1, "decay PDF", "L");

TLine* l2 = new TLine();
l2->SetLineColor(kRed);
l2->SetLineWidth(4);
l2->SetLineStyle(kSolid);
legend->AddEntry(l2, "spline acceptance", "L");

RooPlot* frame_m = Bs_TAU.frame();
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
spl->plotOn(frame_m, LineColor(kRed), Normalization(350, RooAbsReal::Relative));

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

legend->Draw("same");
lhcbtext->DrawTextNDC(0.70,0.55,"LHCb Run 1&2 data");
lhcbtext->DrawTextNDC(0.70,0.45,"preliminary");

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
RooPlot* frame_p = Bs_TAU.frame(Title("pull_frame"));
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

TString* obsTS = new TString(Bs_TAU.GetName());
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

cout << "used datasets:     "<< CutString.c_str() << endl;

file->Close();

return myCoeffs;

}

vector<double> fitSplineAccMC(string CutString){

TFile* file;
TTree* tree;

file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/norm.root");
tree = (TTree*) file->Get("DecayTree");


//Define RooRealVar for observables
RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar Bs_TAU("Bs_TAU", "Bs_TAU", 0.4, 10., "ps");
RooRealVar Bs_TAUERR("Bs_TAUERR", "Bs_TAUERR", 0.0001, 1.,"ps");
RooRealVar weight("weight" , "weight", 0., 12.);
RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0, 3);
RooRealVar year("year", "year", 11, 16);

RooArgSet observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, weight);

RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree), WeightVar(weight.GetName()), Cut(CutString.c_str()));


///SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
//SPLINE KNOTS

double myknots [7] = { 0.5, 1., 1.5, 2., 3., 6., 9.5};
cout << "Bs_TAU_min:   " << Bs_TAU.getMin() << endl;
cout << "Bs_TAU_max:   " << Bs_TAU.getMax() << endl;


std::vector<double> myBinning;

myBinning.push_back(myknots[0]);
myBinning.push_back(myknots[1]);
myBinning.push_back(myknots[2]);
myBinning.push_back(myknots[3]);
myBinning.push_back(myknots[4]);
myBinning.push_back(myknots[5]);
myBinning.push_back(myknots[6]);

//fix the first spline from Bs->Dspipipi data fit 

vector<double> fixed_coeffs;
fixed_coeffs = fitSplineAcc(CutString.c_str());

RooRealVar fixedcoeff_0("fixedcoeff_0", "fixedcoeff_0", fixed_coeffs[0]);
RooRealVar fixedcoeff_1("fixedcoeff_1", "fixedcoeff_1", fixed_coeffs[1]);
RooRealVar fixedcoeff_2("fixedcoeff_2", "fixedcoeff_2", fixed_coeffs[2]);
RooRealVar fixedcoeff_3("fixedcoeff_3", "fixedcoeff_3", fixed_coeffs[3]);
RooRealVar fixedcoeff_4("fixedcoeff_4", "fixedcoeff_4", fixed_coeffs[4]);
RooRealVar fixedcoeff_5("fixedcoeff_5", "fixedcoeff_5", fixed_coeffs[5]);
RooRealVar fixedcoeff_6("fixedcoeff_6", "fixedcoeff_6", fixed_coeffs[6]);
RooRealVar fixedcoeff_7("fixedcoeff_7", "fixedcoeff_7", fixed_coeffs[7]);
RooRealVar fixedcoeff_8("fixedcoeff_8", "fixedcoeff_8", fixed_coeffs[8]);

RooArgList tacc_list_fixed(fixedcoeff_0, fixedcoeff_1, fixedcoeff_2, fixedcoeff_3, fixedcoeff_4, fixedcoeff_5, fixedcoeff_6, fixedcoeff_7, fixedcoeff_8);


double tau = 1.512;
double dgamma = -0.082; //gammaH - gammaL
double DeltaMs = 17.768;

//Setup fixed SPLINE FUNCTION 
RooCubicSplineFun* fixed_spl = new RooCubicSplineFun("fixed splinePdf", "fixed splinePdf", Bs_TAU, myBinning, tacc_list_fixed);



//define floating correction spline

double values [6] =  {3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.7757e+00};


//SPLINE COEFFICIENTS
RooRealVar coeff_0("coeff_0", "coeff_0", values[0], 0.0, 10.0);
RooRealVar coeff_1("coeff_1", "coeff_1", values[1], 0.0, 10.0);
RooRealVar coeff_2("coeff_2", "coeff_2", values[2], 0.0, 10.0);
RooRealVar coeff_3("coeff_3", "coeff_3", values[3], 0.0, 10.0);
RooRealVar coeff_4("coeff_4", "coeff_4", values[4], 0.0, 10.0);
RooRealVar coeff_5("coeff_5", "coeff_5", values[5], 0.0, 10.0);
RooRealVar coeff_55("coeff_55", "coeff_55", 1., 0.0, 10.0);

RooRealVar coeff_6("coeff_6", "coeff_6", 1.0);


RooRealVar knot_5("knot_5","knot_5",myknots[5]);
RooRealVar knot_6("knot_6","knot_6",myknots[6]);
RooRealVar Bs_TAU_max("Bs_TAU_max","Bs_TAU_max", Bs_TAU.getMax());

RooGenericPdf coeff_7("coeff_7","coeff_7", "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(coeff_6, coeff_55, knot_6, knot_5, Bs_TAU_max) );


RooArgList tacc_list(coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_55, coeff_6, coeff_7);


//setup floated spline 
RooCubicSplineFun* R_spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_TAU, myBinning, tacc_list);



///setup product of fixed and floated spline
RooSplineProduct* splineProd = new RooSplineProduct("spline Product","spline Product", Bs_TAU, *fixed_spl, *R_spl);


//combine spline product and time resolution
RooRealVar trm_mean( "trm_mean" , "Gaussian resolution model mean", 0.0, "ps" );
RooRealVar trm_scale( "trm_scale", "Gaussian resolution model scale factor", 1.20);
RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *splineProd, trm_mean, Bs_TAUERR, trm_mean, trm_scale );


RooBDecay* timePDF = new RooBDecay("Bdecay", "Bdecay", Bs_TAU, RooRealConstant::value(tau),
                     RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                     RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(DeltaMs),
                     trm, RooBDecay::SingleSided);

///Fit and Print
//Fit
RooFitResult *myfitresult;
myfitresult = timePDF->fitTo(*dataset, Save(1), Optimize(2), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE));
myfitresult->Print("v");


std::vector<double> myCoeffs;
//values
myCoeffs.push_back(coeff_0.getValV());
myCoeffs.push_back(coeff_1.getValV());
myCoeffs.push_back(coeff_2.getValV());
myCoeffs.push_back(coeff_3.getValV());
myCoeffs.push_back(coeff_4.getValV());
myCoeffs.push_back(coeff_5.getValV());
myCoeffs.push_back(coeff_55.getValV());
myCoeffs.push_back(coeff_6.getVal());
myCoeffs.push_back(coeff_7.getVal());


double range_dw = Bs_TAU.getMin();
double range_up = Bs_TAU.getMax();

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
legend->AddEntry(l1, "decay PDF", "L");

TLine* l3 = new TLine();
l3->SetLineColor(kGreen);
l3->SetLineWidth(4);
l3->SetLineStyle(kSolid);
legend->AddEntry(l3, "fixed spline from datafit", "L");

TLine* l4 = new TLine();
l4->SetLineColor(kMagenta);
l4->SetLineWidth(4);
l4->SetLineStyle(kSolid);
legend->AddEntry(l4, "floated spline", "L");

TLine* l2 = new TLine();
l2->SetLineColor(kRed);
l2->SetLineWidth(4);
l2->SetLineStyle(kSolid);
legend->AddEntry(l2, "spline #otimes spline", "L");


RooPlot* frame_m = Bs_TAU.frame();
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
splineProd->plotOn(frame_m, LineColor(kRed), Normalization(50, RooAbsReal::Relative));
fixed_spl->plotOn(frame_m, LineColor(kGreen), Normalization(50, RooAbsReal::Relative));
R_spl->plotOn(frame_m, LineColor(kMagenta), Normalization(50, RooAbsReal::Relative));

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

legend->Draw("same");
lhcbtext->DrawTextNDC(0.70,0.55,"LHCb Simulation");
lhcbtext->DrawTextNDC(0.70,0.45,"preliminary");

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
RooPlot* frame_p = Bs_TAU.frame(Title("pull_frame"));
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

TString* obsTS = new TString(Bs_TAU.GetName());
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
string PlotName = "timeAccFit_MC_Norm.eps";;

string epsFullName = epsName.append(PlotName);

canvas->SaveAs(epsFullName.c_str());

cout << "used datasets:     "<< CutString.c_str() << endl;

return myCoeffs;

}

void SplineAccCorrection(string CutString){

TFile* file;
TTree* tree;

file= new TFile("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
tree = (TTree*) file->Get("DecayTree");


//Define RooRealVar for observables
RooRealVar DTF_Bs_M("DTF_Bs_M", "DTF_Bs_M", 4975., 5800., "MeV");
RooRealVar Bs_TAU("Bs_TAU", "Bs_TAU", 0.4, 10., "ps");
RooRealVar Bs_TAUERR("Bs_TAUERR", "Bs_TAUERR", 0.0001, 1.,"ps");
RooRealVar weight("weight" , "weight", 0., 12.);
RooRealVar Ds_finalState("Ds_finalState", "Ds_finalState", 0, 3);
RooRealVar year("year", "year", 11, 16);

RooArgSet observables(Bs_TAU, Bs_TAUERR, Ds_finalState, year, weight);

RooDataSet* dataset = new RooDataSet("dataset","dataset", observables, Import(*tree), WeightVar(weight.GetName()), Cut(CutString.c_str()));


///SETUP FITTER AND FIT TO DECAYTIME DISTRIBUTION
//SPLINE KNOTS

double myknots [7] = { 0.5, 1., 1.5, 2., 3., 6., 9.5};
cout << "Bs_TAU_min:   " << Bs_TAU.getMin() << endl;
cout << "Bs_TAU_max:   " << Bs_TAU.getMax() << endl;


std::vector<double> SplineBinning;

SplineBinning.push_back(myknots[0]);
SplineBinning.push_back(myknots[1]);
SplineBinning.push_back(myknots[2]);
SplineBinning.push_back(myknots[3]);
SplineBinning.push_back(myknots[4]);
SplineBinning.push_back(myknots[5]);
SplineBinning.push_back(myknots[6]);

//fix the first spline from Bs->Dspipipi data fit 

vector<double> fixed_coeffs;
fixed_coeffs = fitSplineAccMC(CutString.c_str());

RooRealVar fixedcoeff_0("fixedcoeff_0", "fixedcoeff_0", fixed_coeffs[0]);
RooRealVar fixedcoeff_1("fixedcoeff_1", "fixedcoeff_1", fixed_coeffs[1]);
RooRealVar fixedcoeff_2("fixedcoeff_2", "fixedcoeff_2", fixed_coeffs[2]);
RooRealVar fixedcoeff_3("fixedcoeff_3", "fixedcoeff_3", fixed_coeffs[3]);
RooRealVar fixedcoeff_4("fixedcoeff_4", "fixedcoeff_4", fixed_coeffs[4]);
RooRealVar fixedcoeff_5("fixedcoeff_5", "fixedcoeff_5", fixed_coeffs[5]);
RooRealVar fixedcoeff_6("fixedcoeff_6", "fixedcoeff_6", fixed_coeffs[6]);
RooRealVar fixedcoeff_7("fixedcoeff_7", "fixedcoeff_7", fixed_coeffs[7]);
RooRealVar fixedcoeff_8("fixedcoeff_8", "fixedcoeff_8", fixed_coeffs[8]);

RooArgList tacc_list_fixed(fixedcoeff_0, fixedcoeff_1, fixedcoeff_2, fixedcoeff_3, fixedcoeff_4, fixedcoeff_5, fixedcoeff_6, fixedcoeff_7, fixedcoeff_8);

double tau = 1.512;
double dgamma = -0.082; //gammaH - gammaL
double DeltaMs = 17.768;

//Setup fixed SPLINE FUNCTION 
RooCubicSplineFun* fixed_spl = new RooCubicSplineFun("fixed splinePdf", "fixed splinePdf", Bs_TAU, SplineBinning, tacc_list_fixed);


//define floating correction spline

double values [6] =  {3.1692e-01, 5.9223e-01, 1.1015e+00, 1.3984e+00, 1.7174e+00, 1.7757e+00};


//SPLINE COEFFICIENTS
RooRealVar coeff_0("coeff_0", "coeff_0", values[0], 0.0, 10.0);
RooRealVar coeff_1("coeff_1", "coeff_1", values[1], 0.0, 10.0);
RooRealVar coeff_2("coeff_2", "coeff_2", values[2], 0.0, 10.0);
RooRealVar coeff_3("coeff_3", "coeff_3", values[3], 0.0, 10.0);
RooRealVar coeff_4("coeff_4", "coeff_4", values[4], 0.0, 10.0);
RooRealVar coeff_5("coeff_5", "coeff_5", values[5], 0.0, 10.0);
RooRealVar coeff_55("coeff_55", "coeff_55", 1., 0.0, 10.0);

RooRealVar coeff_6("coeff_6", "coeff_6", 1.0);


RooRealVar knot_5("knot_5","knot_5",myknots[5]);
RooRealVar knot_6("knot_6","knot_6",myknots[6]);
RooRealVar Bs_TAU_max("Bs_TAU_max","Bs_TAU_max", Bs_TAU.getMax());

RooGenericPdf coeff_7("coeff_7","coeff_7", "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(coeff_6, coeff_55, knot_6, knot_5, Bs_TAU_max) );


RooArgList tacc_list(coeff_0, coeff_1, coeff_2, coeff_3, coeff_4, coeff_5, coeff_55, coeff_6, coeff_7);

//setup floated spline 
RooCubicSplineFun* aDataSig_spl = new RooCubicSplineFun("splinePdf", "splinePdf", Bs_TAU, SplineBinning, tacc_list);



///setup product of fixed and floated spline
RooSplineProduct* splineProd = new RooSplineProduct("spline Product","spline Product", Bs_TAU, *fixed_spl, *aDataSig_spl);


//combine spline product and time resolution
RooRealVar trm_mean( "trm_mean" , "Gaussian resolution model mean", 0.0, "ps" );
RooRealVar trm_scale( "trm_scale", "Gaussian resolution model scale factor", 1.20);
RooGaussEfficiencyModel trm("resmodel", "resmodel", Bs_TAU, *splineProd, trm_mean, Bs_TAUERR, trm_mean, trm_scale );


RooBDecay* timePDF = new RooBDecay("Bdecay", "Bdecay", Bs_TAU, RooRealConstant::value(tau),
                     RooRealConstant::value(dgamma), RooRealConstant::value(1.0),  RooRealConstant::value(0.0),
                     RooRealConstant::value(0.0),  RooRealConstant::value(0.0),  RooRealConstant::value(DeltaMs),
                     trm, RooBDecay::SingleSided);


///Fit and Print
//Fit
RooFitResult *myfitresult;
myfitresult = timePDF->fitTo(*dataset, Save(1), Optimize(2), Strategy(2), Verbose(kFALSE), SumW2Error(kTRUE), Extended(kFALSE), Offset(kTRUE));
myfitresult->Print("v");


double range_dw = Bs_TAU.getMin();
double range_up = Bs_TAU.getMax();

TLegend* legend = new TLegend( 0.45, 0.70, 0.88, 0.88 );

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
legend->AddEntry(l1, "decay PDF", "L");

TLine* l3 = new TLine();
l3->SetLineColor(kGreen);
l3->SetLineWidth(4);
l3->SetLineStyle(kSolid);
legend->AddEntry(l3, "spline from B_{s} #rightarrow D_{s}#pi#pi#pi simulation", "L");

TLine* l4 = new TLine();
l4->SetLineColor(kMagenta);
l4->SetLineWidth(4);
l4->SetLineStyle(kSolid);
legend->AddEntry(l4, "floated spline", "L");

TLine* l2 = new TLine();
l2->SetLineColor(kRed);
l2->SetLineWidth(4);
l2->SetLineStyle(kSolid);
legend->AddEntry(l2, "spline #otimes spline", "L");


RooPlot* frame_m = Bs_TAU.frame();
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
splineProd->plotOn(frame_m, LineColor(kRed), Normalization(80, RooAbsReal::Relative));
fixed_spl->plotOn(frame_m, LineColor(kGreen), Normalization(80, RooAbsReal::Relative));
aDataSig_spl->plotOn(frame_m, LineColor(kMagenta), Normalization(80, RooAbsReal::Relative));

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

legend->Draw("same");
lhcbtext->DrawTextNDC(0.60,0.55,"LHCb Simulation");
lhcbtext->DrawTextNDC(0.60,0.45,"preliminary");

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
RooPlot* frame_p = Bs_TAU.frame(Title("pull_frame"));
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

TString* obsTS = new TString(Bs_TAU.GetName());
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
string PlotName = "timeAccFit_MC_Sig.eps";;

string epsFullName = epsName.append(PlotName);

canvas->SaveAs(epsFullName.c_str());

cout << "used datasets:     "<< CutString.c_str() << endl;

}

int main(int argc, char** argv){
    
    time_t startTime = time(0);
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetOptStat(1000000001);

    NamedParameter<int> makePlots("makePlots", 0);
    NamedParameter<string> Ds_finalState("Ds_finalState", (std::string) "all");
    NamedParameter<string> year("year", (std::string) "all");
    NamedParameter<int> CorrectFromMC("CorrectFromMC", 0);
    NamedParameter<int> applyCorrection("applyCorrection", 0);
    string Run1 = "Run1";
    string Run2 = "Run2";
    string all = "all";
    string Ds2KKpi = "Ds2KKpi";
    string Ds2pipipi = "Ds2pipipi";
    string Ds_finalState_convert = Ds_finalState;
    string year_convert = year;

    if(applyCorrection == 1) CorrectFromMC = 0;

if(makePlots == 1)
{
    int numKnots = 9;
    double positions[10] = {0., 0.25, 0.75, 1.25, 1.75, 2.25, 3.75, 8.25, 10.75, 11.};
    double *KnotBinning;
    KnotBinning = positions;

    TH1D *KnotsVsCoeffs_1 = new TH1D("D_{s} -> KK#pi, 2011", "D_{s} -> KK#pi, 2011", numKnots, KnotBinning);
    KnotsVsCoeffs_1->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_1->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_2 = new TH1D("D_{s} -> KK#pi, 2012", "D_{s} -> KK#pi, 2012", numKnots, KnotBinning);
    KnotsVsCoeffs_2->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_2->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_3 = new TH1D("D_{s} -> KK#pi, 2015", "D_{s} -> KK#pi, 2015", numKnots, KnotBinning);
    KnotsVsCoeffs_3->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_3->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_4 = new TH1D("D_{s} -> KK#pi, 2016", "D_{s} -> KK#pi, 2016", numKnots, KnotBinning);
    KnotsVsCoeffs_4->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_4->GetYaxis()->SetTitle("v_{i} values");

    TH1D *KnotsVsCoeffs_5 = new TH1D("D_{s} -> #phi #pi, Run1&2", "D_{s} -> #phi #pi, Run1&2", numKnots, KnotBinning);
    KnotsVsCoeffs_5->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_5->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_6 = new TH1D("D_{s} -> K^{*}K, Run1&2", "D_{s} -> K^{*}K, Run1&2", numKnots, KnotBinning);
    KnotsVsCoeffs_6->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_6->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_7 = new TH1D("D_{s} -> non-resonant, Run1&2", "D_{s} -> non-resonant, Run1&2", numKnots, KnotBinning);
    KnotsVsCoeffs_7->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_7->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_8 = new TH1D("D_{s} -> #pi#pi#pi, Run1&2", "D_{s} -> #pi#pi#pi, Run1&2", numKnots, KnotBinning);
    KnotsVsCoeffs_8->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_8->GetYaxis()->SetTitle("v_{i} values");

    TH1D *KnotsVsCoeffs_9 = new TH1D("D_{s} -> KK#pi, Run1", "D_{s} -> KK#pi, Run1", numKnots, KnotBinning);
    KnotsVsCoeffs_9->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_9->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_10 = new TH1D("D_{s} -> #pi#pi#pi, Run1", "D_{s} -> #pi#pi#pi, Run1", numKnots, KnotBinning);
    KnotsVsCoeffs_10->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_10->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_11 = new TH1D("D_{s} -> KK#pi, Run2", "D_{s} -> KK#pi, Run2", numKnots, KnotBinning);
    KnotsVsCoeffs_11->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_11->GetYaxis()->SetTitle("v_{i} values");
    TH1D *KnotsVsCoeffs_12 = new TH1D("D_{s} -> #pi#pi#pi, Run2", "D_{s} -> #pi#pi#pi, Run2", numKnots, KnotBinning);
    KnotsVsCoeffs_12->GetXaxis()->SetTitle("Knot position [ps]");
    KnotsVsCoeffs_12->GetYaxis()->SetTitle("v_{i} values");

    vector<double> Ds2KKpi_11;
    vector<double> Ds2KKpi_12;
    vector<double> Ds2KKpi_15;
    vector<double> Ds2KKpi_16;

    vector<double> Ds2phiPi;
    vector<double> Ds2KstarK;
    vector<double> Ds2nonRes;
    vector<double> Ds2PiPiPi;

    vector<double> Ds2KKpi_Run1;
    vector<double> Ds2pipipi_Run1;
    vector<double> Ds2KKpi_Run2;
    vector<double> Ds2pipipi_Run2;

    Ds2KKpi_11 = fitSplineAcc("year == 11 && Ds_finalState !=3");
    Ds2KKpi_12 = fitSplineAcc("year == 12 && Ds_finalState !=3");
    Ds2KKpi_15 = fitSplineAcc("year == 15 && Ds_finalState !=3");
    Ds2KKpi_16 = fitSplineAcc("year == 16 && Ds_finalState !=3");

    Ds2phiPi = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && Ds_finalState == 0");
    Ds2KstarK = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && Ds_finalState == 1");
    Ds2nonRes = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && Ds_finalState == 2");
    Ds2PiPiPi = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && Ds_finalState == 3");

    Ds2KKpi_Run1 = fitSplineAcc("(year == 12 || year == 11) && Ds_finalState !=3");
    Ds2pipipi_Run1 = fitSplineAcc("(year == 12 || year == 11) && Ds_finalState ==3");
    Ds2KKpi_Run2 = fitSplineAcc("(year == 15 || year == 16) && Ds_finalState !=3");
    Ds2pipipi_Run2 = fitSplineAcc("(year == 15 || year == 16) && Ds_finalState ==3");

    for(int i = 0; i < numKnots; i++)
    {

	KnotsVsCoeffs_1->SetBinContent((i+1), Ds2KKpi_11[i]);
	KnotsVsCoeffs_2->SetBinContent((i+1), Ds2KKpi_12[i]);
	KnotsVsCoeffs_3->SetBinContent((i+1), Ds2KKpi_15[i]);
	KnotsVsCoeffs_4->SetBinContent((i+1), Ds2KKpi_16[i]);

	KnotsVsCoeffs_5->SetBinContent((i+1), Ds2phiPi[i]);
	KnotsVsCoeffs_6->SetBinContent((i+1), Ds2KstarK[i]);
	KnotsVsCoeffs_7->SetBinContent((i+1), Ds2nonRes[i]);
	KnotsVsCoeffs_8->SetBinContent((i+1), Ds2PiPiPi[i]);

	KnotsVsCoeffs_9->SetBinContent((i+1), Ds2KKpi_Run1[i]);
	KnotsVsCoeffs_10->SetBinContent((i+1), Ds2pipipi_Run1[i]);
	KnotsVsCoeffs_11->SetBinContent((i+1), Ds2KKpi_Run2[i]);
	KnotsVsCoeffs_12->SetBinContent((i+1), Ds2pipipi_Run2[i]);


	if((i+numKnots) != 15 && (i+numKnots) != 16){
		KnotsVsCoeffs_1->SetBinError((i+1), Ds2KKpi_11[i+numKnots]);
		KnotsVsCoeffs_2->SetBinError((i+1), Ds2KKpi_12[i+numKnots]);
		KnotsVsCoeffs_3->SetBinError((i+1), Ds2KKpi_15[i+numKnots]);
		KnotsVsCoeffs_4->SetBinError((i+1), Ds2KKpi_16[i+numKnots]);

		KnotsVsCoeffs_5->SetBinError((i+1), Ds2phiPi[i+numKnots]);
		KnotsVsCoeffs_6->SetBinError((i+1), Ds2KstarK[i+numKnots]);
		KnotsVsCoeffs_7->SetBinError((i+1), Ds2nonRes[i+numKnots]);
		KnotsVsCoeffs_8->SetBinError((i+1), Ds2PiPiPi[i+numKnots]);

		KnotsVsCoeffs_9->SetBinError((i+1), Ds2KKpi_Run1[i+numKnots]);
		KnotsVsCoeffs_10->SetBinError((i+1), Ds2pipipi_Run1[i+numKnots]);
		KnotsVsCoeffs_11->SetBinError((i+1), Ds2KKpi_Run2[i+numKnots]);
		KnotsVsCoeffs_12->SetBinError((i+1), Ds2pipipi_Run2[i+numKnots]);
	}

	else if((i+numKnots) == 15 && (i+numKnots) == 16){
		KnotsVsCoeffs_1->SetBinError((i+1), 0.);
		KnotsVsCoeffs_2->SetBinError((i+1), 0.);
		KnotsVsCoeffs_3->SetBinError((i+1), 0.);
		KnotsVsCoeffs_4->SetBinError((i+1), 0.);

		KnotsVsCoeffs_5->SetBinError((i+1), 0.);
		KnotsVsCoeffs_6->SetBinError((i+1), 0.);
		KnotsVsCoeffs_7->SetBinError((i+1), 0.);
		KnotsVsCoeffs_8->SetBinError((i+1), 0.);

		KnotsVsCoeffs_9->SetBinError((i+1), 0.);
		KnotsVsCoeffs_10->SetBinError((i+1), 0.);
		KnotsVsCoeffs_11->SetBinError((i+1), 0.);
		KnotsVsCoeffs_12->SetBinError((i+1), 0.);
	}
    }


    TCanvas *combinedCanvas_years = new TCanvas("Summary of acceptance analysis", "combinedCanvas_years");
    combinedCanvas_years->Divide(2,2);

    combinedCanvas_years->cd(1);
    KnotsVsCoeffs_1->Draw("E1");
    combinedCanvas_years->cd(2);
    KnotsVsCoeffs_2->Draw("E1");
    combinedCanvas_years->cd(3);
    KnotsVsCoeffs_3->Draw("E1");
    combinedCanvas_years->cd(4);
    KnotsVsCoeffs_4->Draw("E1");

    combinedCanvas_years->SaveAs("Plot/Ds2KKpi_Canvas_altBinning.eps");
    combinedCanvas_years->SaveAs("Plot/Ds2KKpi_Canvas_altBinning.root");



    TCanvas *combinedCanvas_state = new TCanvas("Summary of acceptance analysis", "combinedCanvas_state");
    combinedCanvas_state->Divide(2,2);


    combinedCanvas_state->cd(1);
    KnotsVsCoeffs_5->Draw("E1");
    combinedCanvas_state->cd(2);
    KnotsVsCoeffs_6->Draw("E1");
    combinedCanvas_state->cd(3);
    KnotsVsCoeffs_7->Draw("E1");
    combinedCanvas_state->cd(4);
    KnotsVsCoeffs_8->Draw("E1");

    combinedCanvas_state->SaveAs("Plot/Run12_Canvas_altBinning.eps");
    combinedCanvas_state->SaveAs("Plot/Run12_Canvas_altBinning.root");



    TCanvas *combinedCanvas_state_year = new TCanvas("Summary of acceptance analysis", "combinedCanvas_state_year");
    combinedCanvas_state_year->Divide(2,2);

    combinedCanvas_state_year->cd(1);
    KnotsVsCoeffs_9->Draw("E1");
    combinedCanvas_state_year->cd(2);
    KnotsVsCoeffs_10->Draw("E1");
    combinedCanvas_state_year->cd(3);
    KnotsVsCoeffs_11->Draw("E1");
    combinedCanvas_state_year->cd(4);
    KnotsVsCoeffs_12->Draw("E1");

    combinedCanvas_state_year->SaveAs("Plot/stateRun_Canvas_altBinning.eps");
    combinedCanvas_state_year->SaveAs("Plot/stateRun_Canvas_altBinning.root");

}


    if(makePlots != 1)
    {

    	vector<double> dummy;

    	if(year_convert.compare(all) == 0)
    	{
		if(Ds_finalState_convert.compare(all) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("year == 16 || year == 15 || year == 12 || year == 11");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("year == 16 || year == 15 || year == 12 || year == 11");
		if(applyCorrection == 1) SplineAccCorrection("year == 16 || year == 15 || year == 12 || year == 11");
		}
		else if(Ds_finalState_convert.compare(Ds2KKpi) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState !=3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState !=3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState !=3)");
		}
		else if(Ds_finalState_convert.compare(Ds2pipipi) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState ==3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState ==3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 16 || year == 15 || year == 12 || year == 11) && (Ds_finalState ==3)");
		}
    	}

    	if(Ds_finalState_convert.compare(all) == 0)
    	{
		if(year_convert.compare(all) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("year == 16 || year == 15 || year == 12 || year == 11");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("year == 16 || year == 15 || year == 12 || year == 11");
		if(applyCorrection == 1) SplineAccCorrection("year == 16 || year == 15 || year == 12 || year == 11");
		}
		else if(year_convert.compare(Run1) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("year == 12 || year == 11");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("year == 12 || year == 11");
		if(applyCorrection == 1) SplineAccCorrection("year == 12 || year == 11");
		}
		else if(year_convert.compare(Run2) == 0)
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("year == 15 || year == 16");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("year == 15 || year == 16");
		if(applyCorrection == 1) SplineAccCorrection("year == 15 || year == 16");
		}
    	}

   	if((Ds_finalState_convert.compare(all) != 0) && (year_convert.compare(all) != 0))
   	{
		if((year_convert.compare(Run1) == 0) && (Ds_finalState_convert.compare(Ds2KKpi) == 0))
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 12 || year == 11) && (Ds_finalState !=3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 12 || year == 11) && (Ds_finalState !=3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 12 || year == 11) && (Ds_finalState !=3)");
		}
		else if((year_convert.compare(Run1) == 0) && (Ds_finalState_convert.compare(Ds2pipipi) == 0))
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 12 || year == 11) && (Ds_finalState ==3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 12 || year == 11) && (Ds_finalState ==3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 12 || year == 11) && (Ds_finalState ==3)");
		}
		else if((year_convert.compare(Run2) == 0) && (Ds_finalState_convert.compare(Ds2KKpi) == 0))
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 15 || year == 16) && (Ds_finalState !=3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 15 || year == 16) && (Ds_finalState !=3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 15 || year == 16) && (Ds_finalState !=3)");
		}
		else if((year_convert.compare(Run2) == 0) && (Ds_finalState_convert.compare(Ds2pipipi) == 0))
		{
		if(CorrectFromMC != 1 && applyCorrection != 1) dummy = fitSplineAcc("(year == 15 || year == 16) && (Ds_finalState ==3)");
		if(CorrectFromMC == 1) dummy = fitSplineAccMC("(year == 15 || year == 16) && (Ds_finalState ==3)");
		if(applyCorrection == 1) SplineAccCorrection("(year == 15 || year == 16) && (Ds_finalState ==3)");
		}
   	}
    }
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
