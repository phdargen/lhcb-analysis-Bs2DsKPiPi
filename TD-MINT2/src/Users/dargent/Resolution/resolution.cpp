// Resolution studies
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
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/NamedParameter.h"

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

void getDecayTimeErrorDistribution(string bin, double lowLimit, double highLimit){

///Load files
TChain* tree=new TChain("DecayTree");
tree->Add("/auto/data/kecke/B2DKPiPi/Bs2DsKpipi_MC_fullSel_reweighted_combined.root");

tree->SetBranchStatus("*",0);
tree->SetBranchStatus("Bs_TAU",1);
tree->SetBranchStatus("Bs_cterr",1);
tree->SetBranchStatus("Bs_TRUETAU",1);


Double_t Bs_TAU;
Double_t Bs_TRUETAU;
Double_t Bs_cterr;

tree -> SetBranchAddress( "Bs_TAU" , &Bs_TAU );
tree -> SetBranchAddress( "Bs_TRUETAU" , &Bs_TRUETAU );
tree -> SetBranchAddress( "Bs_cterr" , &Bs_cterr );

string address = "/auto/data/kecke/B2DKPiPi/TimeRes/Bs2DsKpipi_MCcombined";
string line = "_";
string root = ".root";
address.append(line);
address.append(bin);
address.append(root);

TFile* output = 0;
output = new TFile(address.c_str(),"RECREATE");

TTree* summary_tree = tree->CloneTree(0);

double Bs_DeltaTau = 0;
summary_tree->Branch("Bs_DeltaTau",&Bs_DeltaTau,"Bs_DeltaTau/D");

double Sum_cterr = 0;
int Nevents = 0;

//loop over events
int numEvents = tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
	if(Bs_TAU < 0. || Bs_TRUETAU < 0.) continue;
	if( Bs_cterr > lowLimit && Bs_cterr < highLimit){
		Bs_DeltaTau = (1000*Bs_TAU) - (1000*Bs_TRUETAU);
		summary_tree->Fill();
		}

	Sum_cterr = Sum_cterr + Bs_cterr;
	Nevents++;
	}


cout << "average sigma_t:    " <<  (Sum_cterr/Nevents)*1000 << " fs" << endl;


summary_tree->Write();
output->Close();
}

double *FitTimeRes(string ResoBin, string BinName, double resoValues[2]){


string Directory = "/auto/data/kecke/B2DKPiPi/TimeRes/";

string FullName = Directory.append(ResoBin);

TFile* file= new TFile(FullName.c_str());
TTree* tree= (TTree*) file->Get("DecayTree");

RooRealVar Bs_DeltaTau("Bs_DeltaTau", "#Delta(t)", -0.2, 0.2,"[ps]");
RooArgList list =  RooArgList(Bs_DeltaTau);

RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_DeltaTau),Import(*tree));


RooRealVar meanBs1("meanBs1", "B_{s} #mu", 0., -0.2, 0.2);

RooRealVar sigmaBs1("sigmaBs1", "B_{s} #sigma_{1}", 0.010,0.,0.2);
RooRealVar sigmaBs2("sigmaBs2", "B_{s} #sigma_{2}", 0.045,0.,0.2);
RooRealVar sigmaBs3("sigmaBs3", "B_{s} #sigma_{3}", 0.040,0.,0.2);

RooGaussian GaussBs1("GaussBs1", "GaussBs1", Bs_DeltaTau, meanBs1, sigmaBs1);
RooGaussian GaussBs2("GaussBs2", "GaussBs2", Bs_DeltaTau, meanBs1, sigmaBs2);
RooGaussian GaussBs3("GaussBs3", "GaussBs3", Bs_DeltaTau, meanBs1, sigmaBs3);
RooGaussian GaussBs("GaussBs", "GaussBs", Bs_DeltaTau, meanBs1, sigmaBs1);
RooRealVar f_GaussBs("f_GaussBs" , "f__{B_{s}}", 0.5, 0., 1.);
RooRealVar f_GaussBs2("f_GaussBs2" , "2f__{B_{s}}", 0.5, 0., 1.);
RooAddPdf DoubleGaussBs("DoubleGaussBs", "DoubleGaussBs", RooArgList(GaussBs1,GaussBs2),RooArgList(f_GaussBs));
RooAddPdf TripleGaussBs("TripleGaussBs", "TripleGaussBs", RooArgList(GaussBs1,GaussBs2,GaussBs3),RooArgList(f_GaussBs,f_GaussBs2));

RooFitResult *result;
result = DoubleGaussBs.fitTo(*data,Save(kTRUE),Extended(kFALSE),NumCPU(3));
cout << "result is --------------- "<<endl;
result->Print();


double range_dw = Bs_DeltaTau.getMin();
double range_up = Bs_DeltaTau.getMax();
int bin = 75;


RooPlot* frame_m= Bs_DeltaTau.frame();
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


data->plotOn(frame_m,Name("dataSetCut"),MarkerSize(0.5),Binning(bin));
DoubleGaussBs.plotOn(frame_m,Name("FullPdf"),LineColor(kBlack),LineWidth(2));
DoubleGaussBs.plotOn(frame_m,Components(GaussBs1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
DoubleGaussBs.plotOn(frame_m,Components(GaussBs2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));

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

lhcbtext->DrawTextNDC(0.70,0.85,"LHCb simulation");
lhcbtext->DrawTextNDC(0.70,0.75,"preliminary");
//lhcbtext->DrawTextNDC(0.70,0.65,"Bin(19 - 24) fs");

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
RooPlot* frame_p = Bs_DeltaTau.frame(Title("pull_frame"));
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
frame_p->GetXaxis()->SetTitle("#font[132]{#Delta t(B_{s}) [ps]}");

TString* obsTS = new TString(Bs_DeltaTau.GetName());
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
//graph2->Draw("same");
//graph3->Draw("same");

pad2->Update();
canvas->Update();

string epsName = "Plots/";
string epsFullName = epsName.append(BinName);

canvas->SaveAs(epsFullName.c_str());


double f1 = f_GaussBs.getVal();
double df1 = f_GaussBs.getError();
double f2 = 1. - f1;
double df2 = f_GaussBs.getError();

double sig1 = sigmaBs1.getVal();
double dsig1 = sigmaBs1.getError();
double sig2 = sigmaBs2.getVal();
double dsig2 = sigmaBs2.getError();

///compute resolution from combination of Gauss sigmas
double resolution = TMath::Sqrt( (f1*sig1*sig1) + (f2*sig2*sig2) );
double dresolution = TMath::Sqrt(TMath::Power(((sig1*f1*dsig1)/resolution),2)+TMath::Power(((sig2*f2*dsig2)/resolution),2)+TMath::Power(((sig1*sig1*df1)/(2*resolution)),2)+TMath::Power(((sig2*sig2*df2)/(2*resolution)),2));

cout << "Measured resolution from GaussComb:   " << resolution*1000 << " +/- " << dresolution*1000 <<" fs" << endl; 
cout << "***********************************************************************" << endl;



///compute resolution from dilution
double dms = 17.757; // [1/ps]

double dilution = f1*TMath::Exp(-1*(sig1*sig1*dms*dms)/2.) + f2*TMath::Exp(-1*(sig2*sig2*dms*dms)/2.);

double Term1 = f1*TMath::Exp(-1*(sig1*sig1*dms*dms)/2.);
double Term2 = f2*TMath::Exp(-1*(sig2*sig2*dms*dms)/2.);

double ddilution = TMath::Sqrt( TMath::Power(df1*(Term1/f1),2) + TMath::Power(df2*(Term2/f2),2) + TMath::Power((Term1*dms*dms*sig1)*dsig1,2) + TMath::Power((Term2*dms*dms*sig2)*dsig2,2) );

double resolution_eff = TMath::Sqrt((-2/(dms*dms))*log(dilution));
double dresolution_eff = ((2/(dms*dms))/(2*dilution*TMath::Sqrt((-2/(dms*dms))*log(dilution)))) * ddilution;


//cout << "Measured resolution from dilution:   " << resolution_eff*1000 << " +/- " << dresolution_eff*1000 <<" fs" << endl;

resoValues[0] = (resolution_eff*1000);
resoValues[1] = (dresolution_eff*1000);

file->Close();

return resoValues;

}

void FitResoRelation(){

double bins[9] = {0., 19. , 24., 29. , 34. , 39. , 44. , 49., 150.};

//define sigma bins
double *sigma_t_BinCenter;
sigma_t_BinCenter = bins;


/// get the resolutions and errors 
double reso_bin_0_fill[2];
double reso_bin_1_fill[2];
double reso_bin_2_fill[2];
double reso_bin_3_fill[2];
double reso_bin_4_fill[2];
double reso_bin_5_fill[2];
double reso_bin_6_fill[2];
double reso_bin_7_fill[2];

double *reso_bin_0;
double *reso_bin_1;
double *reso_bin_2;
double *reso_bin_3;
double *reso_bin_4;
double *reso_bin_5;
double *reso_bin_6;
double *reso_bin_7;

reso_bin_0 = FitTimeRes("Bs2DsKpipi_MCcombined_0to19.root", "SignalMC_0to19.eps", reso_bin_0_fill);
reso_bin_1 = FitTimeRes("Bs2DsKpipi_MCcombined_19to24.root", "SignalMC_19to24.eps", reso_bin_1_fill);
reso_bin_2 = FitTimeRes("Bs2DsKpipi_MCcombined_24to29.root", "SignalMC_24to29.eps", reso_bin_2_fill);
reso_bin_3 = FitTimeRes("Bs2DsKpipi_MCcombined_29to34.root", "SignalMC_29to34.eps", reso_bin_3_fill);
reso_bin_4 = FitTimeRes("Bs2DsKpipi_MCcombined_34to39.root", "SignalMC_34to39.eps", reso_bin_4_fill);
reso_bin_5 = FitTimeRes("Bs2DsKpipi_MCcombined_39to44.root", "SignalMC_39to44.eps", reso_bin_5_fill);
reso_bin_6 = FitTimeRes("Bs2DsKpipi_MCcombined_44to49.root", "SignalMC_44to49.eps", reso_bin_6_fill);
reso_bin_7 = FitTimeRes("Bs2DsKpipi_MCcombined_49to150.root", "SignalMC_49to150.eps", reso_bin_7_fill);

//pass resolution in bin center and its error
double reso_t_BinCenter[8] = {reso_bin_0[0], reso_bin_1[0], reso_bin_2[0], reso_bin_3[0], reso_bin_4[0], reso_bin_5[0], reso_bin_6[0], reso_bin_7[0]};
double reso_t_BinWidth[8] = {reso_bin_0[1], reso_bin_1[1], reso_bin_2[1], reso_bin_3[1], reso_bin_4[1], reso_bin_5[1], reso_bin_6[1], reso_bin_7[1]};

//define polynom

//TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 45.);
TF1 *fitFunc = new TF1("fitFunc", "[0]+[1]*x ", 0., 150.);
fitFunc->SetParNames("c0","s");
fitFunc->SetParameters(10.,1.2);
fitFunc->SetParLimits(0,0.,30.);
fitFunc->SetParLimits(1,0.,10.);

fitFunc->FixParameter(0,0.);
//fitFunc->FixParameter(1,1.280);

//fill histo
TH1D* ResoRelation = new TH1D("ResoRelation" ,"ResoRelation", 8, sigma_t_BinCenter);

for(int i = 0; i < 8; i++)
	{
	ResoRelation->SetBinContent((i+1), reso_t_BinCenter[i]);
	ResoRelation->SetBinError((i+1), reso_t_BinWidth[i]);
	}

ResoRelation->SetXTitle("decay time error #sigma(t) [fs]");
ResoRelation->SetYTitle("decay time resolution #Delta(t) [fs]");

ResoRelation->SetAxisRange(0.,180.);
ResoRelation->SetAxisRange(0.,120.,"Y");

TCanvas* canvas = new TCanvas();


ResoRelation->Fit(fitFunc,"RL");

ResoRelation->Draw("e1");
fitFunc->Draw("same");
canvas->Print("Plots/ProperTimeReso_MC.eps");

}

int main(int argc, char** argv){

    time_t startTime = time(0);
    gROOT->ProcessLine(".x ../lhcbStyle.C");

    NamedParameter<int> FitReso("FitReso", 1);
    NamedParameter<int> CreateTimeErrorDis("CreateTimeErrorDis", 0);
    NamedParameter<double> lowLimit("lowLimit", 0.);
    NamedParameter<double> highLimit("highLimit", 0.019);
    NamedParameter<string> bin("bin", (std::string) "0to19");

    if(CreateTimeErrorDis == 1) getDecayTimeErrorDistribution(bin, lowLimit, highLimit);

    if(FitReso == 1) FitResoRelation();

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
