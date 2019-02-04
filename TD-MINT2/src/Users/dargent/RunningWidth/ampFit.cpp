// author: Philippe d'Argent (p.dargent@cern.ch)
// status:  Fri 28 Jun 2013 11:21:01 GMT
#include "Mint/NamedParameter.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/SignalGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/FitAmpList.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TF1.h"
#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include "RunningWidthCalculator.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include "TPaveText.h"

using namespace std;
using namespace MINT;

class BW_plotter{
    std::vector<TF1*> _gamma;
    std::vector<double> _c;
    double _m0;
    double _gamma0;
public:
    BW_plotter(std::vector<double> couplings , std::vector<TF1*> f, double m0, double gamma0): 
    _c(couplings), _gamma(f), _m0(m0/GeV), _gamma0(gamma0/GeV){
    }
    double Eval(double x){ 
        double val = 0.;
        double val0 = 0.;
        for (int j=0; j<_gamma.size(); j++) {
            val += _c[j]*_gamma[j]->Eval(x);
            val0 += _c[j]*_gamma[j]->Eval(_m0*_m0);
        }
        double gamma = val/val0*_gamma0;
        return _m0*_gamma0/((_m0*_m0-x)*(_m0*_m0-x)+(_m0*gamma)*(_m0*gamma));
    }
    double operator()(double* x, double* p){
        return Eval(x[0]);
    }

};

int makeHistosForD4pi(){  
  TCanvas* c = new TCanvas();  
    
  NamedParameter<int> EventPattern("Event Pattern", 20213, 211, 211, -211);
  DalitzEventPattern pat(EventPattern.getVector());
  cout << " got event pattern: " << pat << endl;

  NamedParameter<int> nBins("nBins",2000);
  NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
  NamedParameter<double> max_s_inGeV2("max_s_inGeV2",3.1);  
      
  RunningWidthCalculator rwcalc(pat);
  
  // Integrate over a1 Dalitz plot   
  FitAmpSum fas(pat);
  fas.print();   
  //TH1D* h_dalitz = rwcalc.makeHisto_dalitz(nBins, max_s_inGeV2, nIntegrationEvents, &fas);
    
  // Integrate over flat phase space   
  //TH1D* h_phaseSpace = rwcalc.makeHisto_phaseSpace(nBins, max_s_inGeV2);  
   //return 0;    


  // Calculate a1->rho pi 3-body width 
  DecayTree dk_a1_rhopi = DecayTree(20213);
  dk_a1_rhopi.addDgtr(211, 113)->addDgtr(211,-211);  
  TH1D* h_a1_rhopi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_rhopi,"RunningWidth_a1_rhopi.root");
    
  dk_a1_rhopi.getVal().setL(2);  
  TH1D* h_a1_rhopi_D = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_rhopi,"RunningWidth_a1_rhopi_D.root");
//return 0;
  DecayTree dk_a1_sigmapi = DecayTree(20213);
  dk_a1_sigmapi.addDgtr(211, 999001)->addDgtr(211,-211);   
  TH1D* h_a1_sigmapi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_sigmapi,"RunningWidth_a1_sigmapi.root");
  
  // Calculate a1->K*(892) K 3-body width 
  DecayTree dk_a1_KstarK = DecayTree(20213);
  dk_a1_KstarK.addDgtr(321, 323)->addDgtr(321,-211); 
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK));  
  TH1D* h_a1_KstarK = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_KstarK, "RunningWidth_a1_KstarK.root");
    
  // Calculate a1->K*(892) K 2-body width 
  DecayTree dk_a1_KstarK_2body = DecayTree(20213);
  dk_a1_KstarK_2body.addDgtr(321, 323); 
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK_2body));  
  TH1D* h_a1_KstarK_2body = rwcalc.makeHisto_2body(nBins, max_s_inGeV2, dk_a1_KstarK_2body, "RunningWidth_a1_KstarK_2body.root");  

  dk_a1_rhopi.getVal().setL(2);    
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_rhopi));  
  TF1* f_a1_rhopi_D = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_rhopi);  
  f_a1_rhopi_D->SetName("f_a1_rhopi_D");
    
  dk_a1_rhopi.getVal().setL(0);  
  TF1* f_a1_rhopi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_rhopi);  
  f_a1_rhopi->SetName("f_a1_rhopi");
 
  TF1* f_a1_sigmapi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_sigmapi);  
  f_a1_sigmapi->SetName("f_a1_sigmapi");
    
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK));  
  TF1* f_a1_KstarK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_KstarK);
  f_a1_KstarK->SetName("f_a1_KstarK");
    
  vector<TF1*> gamma_tot;
  gamma_tot.push_back(f_a1_rhopi);
  gamma_tot.push_back(f_a1_rhopi_D);
  gamma_tot.push_back(f_a1_sigmapi);
  gamma_tot.push_back(f_a1_KstarK);
    
  vector<double> BF;
  // 1. iteration
  /*  
  BF.push_back(0.77);
  BF.push_back(0.05);
  BF.push_back(0.15);
  BF.push_back(0.03);  
  double m0 = 1255;
  double gamma0 = 367.;
  */
    
  // 2. iteration  
  /*  
  BF.push_back(33.6/(33.6+2.4+13.8)*0.967);
  BF.push_back(2.4/(33.6+2.4+13.8)*0.967);
  BF.push_back(13.8/(33.6+2.4+13.8)*0.967);
  BF.push_back(0.033);  
  double m0 = 1222;
  double gamma0 = 504.;
  */
  
  /*
  // 3. iteration  
  BF.push_back(41.46/(41.46+6.54)*0.967);
  BF.push_back(0./(41.46+6.54)*0.967);
  BF.push_back(6.54/(41.46+6.54)*0.967);
  BF.push_back(0.033);  
  double m0 = 1220;
  double gamma0 = 412.;   
  */
    
  /*
  // 4. iteration  
  BF.push_back(38.25/(38.25+6.63)*0.967);
  BF.push_back(0./(38.25+6.63)*0.967);
  BF.push_back(6.63/(38.25+6.63)*0.967);
  BF.push_back(0.033);  
  double m0 = 1210;
  double gamma0 = 408.;
  */
  /*
  // 5. iteration  
  BF.push_back(38.98/(38.98+6.62)*0.967);
  BF.push_back(0./(38.98+6.62)*0.967);
  BF.push_back(6.62/(38.98+6.62)*0.967);
  BF.push_back(0.033);  
  double m0 = 1210.;
  double gamma0 = 415.;
  cout << "BF = " <<  BF  << endl;
  */
  /*
  // 6. iteration  
  BF.push_back(38.84/(38.84+6.57)*0.967);
  BF.push_back(0./(38.84+6.57)*0.967);
  BF.push_back(6.57/(38.84+6.57)*0.967);
  BF.push_back(0.033);  
  double m0 = 1211.;
  double gamma0 = 416.;
  */
  /*
  // 7. iteration  
  BF.push_back(39.45/(39.45+6.38)*0.967);
  BF.push_back(0./(39.45+6.38)*0.967);
  BF.push_back(6.38/(39.45+6.38)*0.967);
  BF.push_back(0.033);  
  double m0 = 1212.;
  double gamma0 = 417.;
  */
  /*
  // 8. iteration  
  BF.push_back(39.94/(39.94+6.32)*0.967);
  BF.push_back(0./(39.94+6.32)*0.967);
  BF.push_back(6.32/(39.94+6.32)*0.967);
  BF.push_back(0.033);  
  double m0 = 1210.;
  double gamma0 = 414.;
  */
  // 9. it
  /*
  BF.push_back(39.26/(39.26+6.38)*0.967);
  BF.push_back(0./(39.26+6.38)*0.967);
  BF.push_back(6.32/(39.26+6.38)*0.967);
  BF.push_back(0.033);  
  double m0 = 1210.;
  double gamma0 = 415.;
  */
  // 10. it
  /*
  BF.push_back(38.9/(38.9+7.4)*0.967);
  BF.push_back(0./(38.9+7.4)*0.967);
  BF.push_back(7.4/(38.9+7.4)*0.967);
  BF.push_back(0.033);  
  double m0 = 1226.;
  double gamma0 = 442.;
  */
  // 11. it
 /* 
  BF.push_back(39.7/(39.7+6.9)*0.967);
  BF.push_back(0./(39.7+6.9)*0.967);
  BF.push_back(6.9/(39.7+6.9)*0.967);
  BF.push_back(0.033);  
  double m0 = 1227.;
  double gamma0 = 440.;
*/
  /*
  // 11. it
  BF.push_back(36.06/(36.06+11.73)*0.967);
  BF.push_back(0./(39.7+6.9)*0.967);
  BF.push_back(11.73/(36.06+11.73)*0.967);
  BF.push_back(0.033);  
  double m0 = 1230.;
  double gamma0 = 452.;
  */
  // 12. it
  BF.push_back(38.1/(38.1+10.24)*0.967);
  BF.push_back(0./(39.7+6.9)*0.967);
  BF.push_back(10.24/(38.1+10.24)*0.967);
  BF.push_back(0.033);  
  double m0 = 1225.;
  double gamma0 = 430.;
  

  cout << "BF = " <<  BF  << endl;
  
  vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
  //vector<double> couplings;
  /*couplings.push_back(0.848592);
  couplings.push_back(0.0712339);
  couplings.push_back(0.201538);
  couplings.push_back(0.122192);
  */
  /*  
  couplings.push_back(1.09458);
  couplings.push_back(0.0730657);
  couplings.push_back(0.522786);
  couplings.push_back(0.174616);  
   */
  cout << "couplings = " <<  couplings  << endl;
    
  TH1D* gamma = (TH1D*) h_a1_rhopi->Clone();
  TH1D* gamma_rhopi = (TH1D*) h_a1_rhopi->Clone();
  TH1D* gamma_rhopi_D = (TH1D*) h_a1_rhopi->Clone(); 
  TH1D* gamma_sigmapi = (TH1D*) h_a1_rhopi->Clone();
  TH1D* gamma_KstarK = (TH1D*) h_a1_rhopi->Clone();
  TH1D* gamma_3pi = (TH1D*) h_a1_rhopi->Clone();
    
  for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        gamma_rhopi->SetBinContent(i,couplings[0]*f_a1_rhopi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_rhopi_D->SetBinContent(i,couplings[1]*f_a1_rhopi_D->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_sigmapi->SetBinContent(i,couplings[2]*f_a1_sigmapi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_KstarK->SetBinContent(i,couplings[3]*f_a1_KstarK->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_3pi->SetBinContent(i,couplings[0]*f_a1_rhopi->Eval(gamma->GetXaxis()->GetBinCenter(i))
                                 + couplings[1]*f_a1_rhopi_D->Eval(gamma->GetXaxis()->GetBinCenter(i))
                                 + couplings[2]*f_a1_sigmapi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
   }  
    
   TFile* f = new TFile("RunningWidth_a(1)(1260)+_12it2.root","RECREATE");
   gamma->Write();
   //gamma_KstarK->Write();   
    
   gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
   gamma->SetLineColor(kBlack);
   gamma->Draw("histc");
   gamma_3pi->SetLineStyle(kDashed);
   gamma_3pi->SetLineColor(kBlue);
   gamma_3pi->Draw("histcsame");
   gamma_rhopi->SetLineStyle(kDashed);
   gamma_rhopi->SetLineColor(kBlue);
   //gamma_rhopi->Draw("same");
   gamma_rhopi_D->SetLineStyle(kDashed);
   gamma_rhopi_D->SetLineColor(kGreen);
   //gamma_rhopi_D->Draw("same");
   gamma_sigmapi->SetLineStyle(kDashed);
   gamma_sigmapi->SetLineColor(kMagenta);
   //gamma_sigmapi->Draw("same");
   gamma_KstarK->SetLineStyle(kDotted);
   gamma_KstarK->SetLineColor(kRed);
   gamma_KstarK->Draw("histcsame");
   //c->Print("gamma_a1_1260.eps");    

        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(a)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    //c->Print("plots/gamma_a1_1260.eps");

    gamma->Draw("");
    //c->Print("plots/runningWidth_a1_1260.C");


   
   rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_rhopi));    
   TH1D* h_m = rwcalc.makeRunningMassHisto_3body(20, max_s_inGeV2, gamma_tot, couplings);  
    
   double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
      
   for (int i=1; i<= h_m->GetNbinsX(); i++) {
     h_m->SetBinContent(i,((h_m->GetBinContent(i) - delta_m0) *m0/GeV + m0/GeV*m0/GeV) );
   }    
     
    h_m->Write();
    
   h_m->Draw("histc");
   c->Print("m_a1.eps");  
       
    
  f->Write(); 
  f->Close(); 
  
  // Calculate pi->rho pi 3-body width 
  DecayTree dk_pi_rhopi = DecayTree(100211);
  dk_pi_rhopi.addDgtr(211, 113)->addDgtr(211,-211); 
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_pi_rhopi));  
  TH1D* h_pi_rhopi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_pi_rhopi,"RunningWidth_pi_rhopi.root");
    
  // Calculate pi->sigma pi 3-body width 
  DecayTree dk_pi_sigmapi = DecayTree(100211);
  dk_pi_sigmapi.addDgtr(211, 999001)->addDgtr(211,-211);   
  TH1D* h_pi_sigmapi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_pi_sigmapi,"RunningWidth_pi_sigmapi.root");  
    
  // Add pi histos
  TH1D* gamma_pi = (TH1D*) h_pi_rhopi->Clone();
  double b_pi = gamma_pi->FindBin(1.69);
  double norm_pi_1 = h_pi_rhopi->GetBinContent(b_pi);
  double norm_pi_2 = h_pi_sigmapi->GetBinContent(b_pi);

  for (int i=1; i<= gamma_pi->GetNbinsX(); i++) {
        h_pi_rhopi->SetBinContent(i,h_pi_rhopi->GetBinContent(i)/norm_pi_1*0.4*0.5);
        h_pi_sigmapi->SetBinContent(i,h_pi_sigmapi->GetBinContent(i)/norm_pi_2*0.4*0.5);
  }
    
  for (int i=1; i<= gamma_pi->GetNbinsX(); i++) {
        double val = h_pi_rhopi->GetBinContent(i) +  h_pi_sigmapi->GetBinContent(i);
        gamma_pi->SetBinContent(i,val);
  }
    
  gamma_pi->SetMinimum(0.);
  gamma_pi->SetLineColor(kRed);
  gamma_pi->SetTitle("; s [GeV^{2}] ; #Gamma(s) [GeV]");  
  gamma_pi->Draw("");
  h_pi_rhopi->SetLineStyle(kDashed);
  h_pi_sigmapi->SetLineStyle(kDashed);  
  h_pi_rhopi->Draw("same");
  h_pi_sigmapi->Draw("same");
  c->Print("gamma_pi.eps");  

  TFile* f_pi = new TFile("RunningWidth_pi(1300)+.root","RECREATE");
  gamma_pi->Write();  
  f_pi->Write(); 
  f_pi->Close();
    
  // Calculate pi->rho pi 3-body width 
  DecayTree dk_a2_rhopi = DecayTree(215);
  dk_a2_rhopi.addDgtr(211, 113)->addDgtr(211,-211); 
  rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a2_rhopi));  
  TH1D* h_a2_rhopi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a2_rhopi,"RunningWidth_a(2)(1320)+.root");  
    
  return 0;
}

int make_pi_1300_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 100211, 211, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",3.1);  
    
    RunningWidthCalculator rwcalc(pat);
    
    /*
    // Integrate over pi Dalitz plot   
    FitAmpSum fas(pat);
    fas.print();   
    DalitzEventList eventTest;
    eventTest.generatePhaseSpaceEvents(200000,pat);
    fas.normalizeAmps(eventTest);
    TH1D* h_dalitz = rwcalc.makeHisto_dalitz(nBins, max_s_inGeV2, nIntegrationEvents, &fas);
    return 0;
    */

    // Calculate pi->rho pi 3-body width 
    DecayTree dk_1 = DecayTree(100211);
    dk_1.addDgtr(211, 113)->addDgtr(211,-211);  
    TF1* f_1 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_1);
    
    // Calculate pi->sigma pi 3-body width 
    DecayTree dk_2 = DecayTree(100211);
    dk_2.addDgtr(211, 999001)->addDgtr(211,-211); 
    TF1* f_2 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_2);
    
    // Calculate pi->NV pi 3-body width 
    DecayTree dk_3 = DecayTree(100211);
    dk_3.addDgtr(211, 9993)->addDgtr(211,-211);  
    TF1* f_3 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_3);


    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_3);
    gamma_tot.push_back(f_2);
    vector<double> BF;

    //pdg
    /*
    //BF.push_back(1./(1.+2.2));
    //BF.push_back(2.2/(1.+2.2));
    BF.push_back(0.4);
    BF.push_back(0.6);
    double m0 = 1200.;
    double gamma0 = 470.;
    */

    //1.it
    /*
    BF.push_back(5.4/(3.6+5.4));
    BF.push_back(3.6/(3.6+5.4));
    double m0 = 1182.;
    double gamma0 = 313.;
    */

    //2.it
  /*  BF.push_back(7.3/(4.3+7.3));
    BF.push_back(4.3/(4.3+7.3));
    double m0 = 1173.;
    double gamma0 = 290.;
*/
/*
    //3.it
    BF.push_back(5.5/(5.5+3.7));
    BF.push_back(3.7/(5.5+3.7));
    double m0 = 1180.;
    double gamma0 = 308.;
*/
    //4.it
    BF.push_back(0.);
    BF.push_back(1.);
    double m0 = 1128.;
    double gamma0 = 314.;


    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
     couplings.push_back(0.166291);
     couplings.push_back(0.678194);
    */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_1 = new TH1D("RunningWidth pi(1300)","; s (GeV^{2}) ; #Gamma(s) (GeV)", nBins, 0., max_s_inGeV2);
    TH1D* h_2 = (TH1D*) h_1->Clone();    
    TH1D* gamma = (TH1D*) h_1->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        h_1->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_2->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }  

    //cout << "Gamma(sigma pi)/Gamma(rho pi) = " << h_2->GetBinContent(gamma->FindBin(pow(m0/GeV,2)))/h_1->GetBinContent(gamma->FindBin(pow(m0/GeV,2))) << endl;
    
    gamma->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma->SetLineColor(kBlack);
    gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
    gamma->Draw("");
    h_1->SetLineColor(kBlue);
    h_1->SetLineStyle(kDashed);
    //h_1->Draw("same");
    h_2->SetLineColor(kBlack);
    h_2->SetLineStyle(kDashed);
    //h_2->Draw("same");
    
    /*
     BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
     TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
     TH1D* h_bw = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw->SetLineColor(kBlack);
     h_bw->SetLineStyle(kDashed);
     h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
     h_bw->Draw("same");
     
     */
    
    /*
     TF1* const_bw = new TF1("const_bw","1",0,max_s_inGeV2,0);
     vector<TF1*> gamma_const;
     gamma_const.push_back(const_bw);
     
     BW_plotter bw_plotter_const(vector<double>(1),gamma_const,m0,gamma0);
     TF1 *bw_const = new TF1("bw_const",bw_plotter_const,0,max_s_inGeV2,0);
     TH1D* h_bw_const = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw_const->SetBinContent(i,bw_const->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw_const->SetLineColor(kBlack);
     h_bw_const->SetLineStyle(kDashed);
     h_bw_const->Scale(gamma->GetMaximum()*0.8/bw_const->Eval(m0*m0/(GeV*GeV)));
     h_bw_const->Draw("same");
     */
    
        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(a)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_pi_1300.eps");

    gamma->Draw("");
    c->Print("plots/runningWidth_pi_1300.C");
    
    TFile* out = new TFile("RunningWidth_pi(1300)+_5it_mod.root","RECREATE");
    gamma->Write();  
    out->Write(); 
    out->Close();
    
    /*
     
     //rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
     TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
     
     double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
     
     for (int i=1; i<= h_m->GetNbinsX(); i++) {
     h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
     }    
     
     h_m->Draw("");
     c->Print("m_K1s_1680.pdf");
     
     */
    
    return 0;
}

int make_a1_1640_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 20213, 211, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4.1);  
    
    RunningWidthCalculator rwcalc(pat);
    
    
    // Calculate a1->rho pi 3-body width 
    DecayTree dk_a1_rhopi = DecayTree(20213);
    dk_a1_rhopi.addDgtr(211, 113)->addDgtr(211,-211);  
    TH1D* h_a1_rhopi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_rhopi);
    
    dk_a1_rhopi.getVal().setL(2);  
    TH1D* h_a1_rhopi_D = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_rhopi,"RunningWidth_a1_rhopi_D.root");
    
    DecayTree dk_a1_sigmapi = DecayTree(20213);
    dk_a1_sigmapi.addDgtr(211, 999001)->addDgtr(211,-211);   
    TH1D* h_a1_sigmapi = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_sigmapi,"RunningWidth_a1_sigmapi.root");
    
    // Calculate a1->K*(892) K 3-body width 
    DecayTree dk_a1_KstarK = DecayTree(20213);
    dk_a1_KstarK.addDgtr(321, 323)->addDgtr(321,-211); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK));  
    TH1D* h_a1_KstarK = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_a1_KstarK, "RunningWidth_a1_KstarK.root");
    
    // Calculate a1->K*(892) K 2-body width 
    DecayTree dk_a1_KstarK_2body = DecayTree(20213);
    dk_a1_KstarK_2body.addDgtr(321, 323); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK_2body));  
    TH1D* h_a1_KstarK_2body = rwcalc.makeHisto_2body(nBins, max_s_inGeV2, dk_a1_KstarK_2body, "RunningWidth_a1_KstarK_2body.root");  
    
    dk_a1_rhopi.getVal().setL(2);    
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_rhopi));  
    TF1* f_a1_rhopi_D = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_rhopi);  
    f_a1_rhopi_D->SetName("f_a1_rhopi_D");
    
    dk_a1_rhopi.getVal().setL(0);  
    TF1* f_a1_rhopi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_rhopi);  
    f_a1_rhopi->SetName("f_a1_rhopi");
    
    TF1* f_a1_sigmapi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_sigmapi);  
    f_a1_sigmapi->SetName("f_a1_sigmapi");
    
    DecayTree dk_a1_f2pi = DecayTree(20213);
    dk_a1_f2pi.addDgtr(211, 225)->addDgtr(211,-211);  
    TF1* f_a1_f2pi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_f2pi);  
    f_a1_f2pi->SetName("f_a1_f2pi");
    
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_KstarK));  
    TF1* f_a1_KstarK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_a1_KstarK);
    f_a1_KstarK->SetName("f_a1_KstarK");
    
    vector<TF1*> gamma_tot;
    //gamma_tot.push_back(f_a1_rhopi);
    gamma_tot.push_back(f_a1_rhopi_D);
    gamma_tot.push_back(f_a1_sigmapi);
    gamma_tot.push_back(f_a1_f2pi);
    //gamma_tot.push_back(f_a1_KstarK);
    
    vector<double> BF;

     // a1(1640)
     //BF.push_back(2.9/(2.9+2.9));
     //BF.push_back(2.9/(2.9+2.9));
     //BF.push_back(0.033);  
    
     /*
     BF.push_back(4.17/(4.17+1.86));
     BF.push_back(1.86/(4.17+1.86));
     double m0 = 1647;
     double gamma0 = 254;  
     */

     /*
     //2.it
     BF.push_back(3.9/(3.9+1.89));
     BF.push_back(1.89/(3.9+1.89));
     double m0 = 1647;
     double gamma0 = 254;  
     */

     /*
     //3.it
     BF.push_back(3.52/(3.52+1.86));
     BF.push_back(1.86/(3.52+1.86));
     double m0 = 1647;
     double gamma0 = 254;  
 
     //4.it
     BF.push_back(2.97/(2.97+1.91));
     BF.push_back(1.91/(2.97+1.91));
     double m0 = 1647;
     double gamma0 = 254;  
	*/
     /*
     //5.it
     BF.push_back(3.68/(3.68+1.72));
     BF.push_back(1.72/(3.68+1.72));
     double m0 = 1647;
     double gamma0 = 254;  
     */
	/*
    //5.it mod
     BF.push_back(5.29/(5.29+2.88));
     BF.push_back(2.88/(5.29+2.88));
     double m0 = 1729.;
     double gamma0 =161.;  
*/
     /*
    //5.it mod 2
     BF.push_back(0.68/(3.68+1.72+0.68+1.));
     BF.push_back(3.68/(3.68+1.72+0.68+1.));
     BF.push_back(1.72/(3.68+1.72+0.68+1.));
     BF.push_back(1./(3.68+1.72+0.68+1.));
     double m0 = 1647;
     double gamma0 = 254;  
     cout << "BF = " <<  BF  << endl;
    */
   
     /*
     //6.it
     BF.push_back(3.70/(3.7+1.81));
     BF.push_back(1.81/(3.7+1.81));
     double m0 = 1647;
     double gamma0 = 254;  
     */
     /*
     //7a.it
     BF.push_back(3.80/(3.8+1.9));
     BF.push_back(0.9/(3.8+.9));
     double m0 = 1647;
     double gamma0 = 254;  
     */
     
     //7b.it
     /*
     BF.push_back(4.2/(4.2+1.6));
     BF.push_back(1.6/(4.2+1.6));
     double m0 = 1707;
     double gamma0 = 138;  
     */
     
     //8.it
     /*
     BF.push_back(3.9/(3.9+.1));
     BF.push_back(0.1/(3.9+0.1));
     double m0 = 1656;
     double gamma0 = 277;  
     */
     //9.it
 /*    BF.push_back(3.9/(3.9+1.+1.));
     BF.push_back(1./(3.9+1.+1.));
     BF.push_back(1./(3.9+1.+1.));
     double m0 = 1656;
     double gamma0 = 254;  
*/
/*
     //10.it
     BF.push_back(3.4/(3.4+1.1));
     BF.push_back(1.1/(3.4+1.1));
     BF.push_back(0./(3.9+1.1));
     double m0 = 1647;
*/
     //11.it
     BF.push_back(4.2/(4.2+2.4));
     BF.push_back(2.4/(4.2+2.4));
     BF.push_back(0./(3.9+1.1));
     double m0 = 1691.;
     double gamma0 = 171.;            

     cout << "BF = " <<  BF  << endl;

     vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //vector<double> couplings;
    /*couplings.push_back(0.848592);
     couplings.push_back(0.0712339);
     couplings.push_back(0.201538);
     couplings.push_back(0.122192);
     */
    /*
    couplings.push_back(1.09458);
    couplings.push_back(0.0730657);
    couplings.push_back(0.522786);
    couplings.push_back(0.174616);  
    */
    cout << "couplings = " <<  couplings  << endl;
    
    TH1D* gamma = (TH1D*) h_a1_rhopi->Clone();
    TH1D* gamma_rhopi = (TH1D*) h_a1_rhopi->Clone();
    TH1D* gamma_rhopi_D = (TH1D*) h_a1_rhopi->Clone(); 
    TH1D* gamma_sigmapi = (TH1D*) h_a1_rhopi->Clone();
    TH1D* gamma_f2pi = (TH1D*) h_a1_rhopi->Clone();
    TH1D* gamma_KstarK = (TH1D*) h_a1_rhopi->Clone();
    TH1D* gamma_3pi = (TH1D*) h_a1_rhopi->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        //gamma_rhopi->SetBinContent(i,couplings[0]*f_a1_rhopi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_rhopi_D->SetBinContent(i,couplings[0]*f_a1_rhopi_D->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_f2pi->SetBinContent(i,couplings[1]*f_a1_f2pi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        //gamma_sigmapi->SetBinContent(i,couplings[2]*f_a1_sigmapi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        //gamma_KstarK->SetBinContent(i,couplings[3]*f_a1_KstarK->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        gamma_3pi->SetBinContent(i,//couplings[0]*f_a1_rhopi->Eval(gamma->GetXaxis()->GetBinCenter(i))
                                  couplings[0]*f_a1_rhopi_D->Eval(gamma->GetXaxis()->GetBinCenter(i))
                                + couplings[1]*f_a1_f2pi->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }  
    
    TFile* f = new TFile("RunningWidth_a(1)(1640)+_11it.root","RECREATE");
    gamma->Write();
    //gamma_KstarK->Write();   
    
    gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
    gamma->SetLineColor(kBlack);
    gamma->Draw("histc");
    gamma_3pi->SetLineStyle(kDashed);
    gamma_3pi->SetLineColor(kBlue);
    //gamma_3pi->Draw("histcsame");
    //gamma_rhopi->SetLineStyle(kDashed);
    //gamma_rhopi->SetLineColor(kBlue);
    //gamma_rhopi->Draw("same");
    gamma_rhopi_D->SetLineStyle(kDashed);
    gamma_rhopi_D->SetLineColor(kGreen);
    //gamma_rhopi_D->Draw("same");
    //gamma_sigmapi->SetLineStyle(kDashed);
    //gamma_sigmapi->SetLineColor(kMagenta);
    //gamma_sigmapi->Draw("same");
    gamma_f2pi->SetLineStyle(kDashed);
    gamma_f2pi->SetLineColor(kMagenta);
    //gamma_f2pi->Draw("same");
    //gamma_KstarK->SetLineStyle(kDotted);
    //gamma_KstarK->SetLineColor(kRed);
    //gamma_KstarK->Draw("histcsame");

/*
     BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
     TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
     TH1D* h_bw = (TH1D*) gamma->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw->SetLineColor(kBlack);
     h_bw->SetLineStyle(kDashed);
     h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
     h_bw->Draw("same");
     
*/


        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(b)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_a1_1640.eps");

    gamma->Draw("");
    c->Print("plots/runningWidth_a1_1640.C");
    
    /*
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_a1_rhopi));    
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(50, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,sqrt((h_m->GetBinContent(i) - delta_m0) *m0/GeV + m0/GeV*m0/GeV) );
    }    
    
    h_m->Write();
    
    h_m->Draw("histc");
    c->Print("m_a1_1640.pdf");  
    
    */
    
    f->Write(); 
    f->Close(); 
    
    return 0;
}

int make_pi2_1670_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 10215, 211, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4.1);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate pi2->rho pi 3-body width 
    DecayTree dk_1 = DecayTree(10215);
    dk_1.addDgtr(211, 113)->addDgtr(211,-211);  
    TF1* f_1 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_1);
    
    // Calculate pi2->f2 pi 3-body width 
    DecayTree dk_2 = DecayTree(10215);
    dk_2.addDgtr(211, 225)->addDgtr(211,-211); 
    TF1* f_2 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_2);
    
    // Calculate pi2->sigma pi 3-body width 
    DecayTree dk_3 = DecayTree(10215);
    dk_3.addDgtr(211, 999001)->addDgtr(211,-211); 
    TF1* f_3 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_3);
    
    // Calculate pi2->rho omega 3-body width 
    DecayTree dk_4 = DecayTree(10215);
    dk_4.addDgtr(223, 213)->addDgtr(211,111); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_4));  
    TF1* f_4 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_4);
    
    // Calculate pi2->K*(892) K 3-body width     
    DecayTree dk_5 = DecayTree(10215);
    dk_5.addDgtr(321, 323)->addDgtr(321,-211); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_5));  
    TF1* f_5 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_5);
    
    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_1);
    gamma_tot.push_back(f_2);
    gamma_tot.push_back(f_3);
    gamma_tot.push_back(f_4);
    gamma_tot.push_back(f_5);
    
    vector<double> BF;

    
    //pdg
    BF.push_back(0.310/(0.31 +0.563 + 0.109 )* (1.- (0.027+0.042)));
    BF.push_back(0.563/(0.31 +0.563 + 0.109 )* (1.- (0.027+0.042)));
    BF.push_back(0.109/(0.31 +0.563 + 0.109 )* (1.- (0.027+0.042))); 
    //BF.push_back((1.- (0.027+0.042)));  
    BF.push_back(0.027);  
    BF.push_back(0.042);  
    double m0 = 1672.2;
    double gamma0 = 260.;
    
    /*
    //2.it
    //BF.push_back(0.0/(0. +0.564 + 0.311 )* (1.- (0.027+0.042)));
    BF.push_back(0.516/(0. +0.516 + 0.318 )* (1.- (0.027+0.042)));
    BF.push_back(0.318/(0. +0.564 + 0.318 )* (1.- (0.027+0.042))); 
    //BF.push_back((1.- (0.027+0.042)));  
    BF.push_back(0.027);  
    BF.push_back(0.042);  
    double m0 = 1672.2;
    double gamma0 = 260.;
    */
    /*
    //4.it
    //BF.push_back(0.0/(0. +0.564 + 0.311 )* (1.- (0.027+0.042)));
    BF.push_back(0.519/(0. +0.519 + 0.311 )* (1.- (0.027+0.042)));
    BF.push_back(0.311/(0. +0.519 + 0.311 )* (1.- (0.027+0.042))); 
    //BF.push_back((1.- (0.027+0.042)));  
    BF.push_back(0.027);  
    BF.push_back(0.042);  
    double m0 = 1672.2;
    double gamma0 = 260.;
    */

    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
     couplings.push_back(0.166291);
     couplings.push_back(0.678194);
     couplings.push_back(0.0406787);
    */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_1 = new TH1D("RunningWidth pi_2(1670)","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_2 = (TH1D*) h_1->Clone();
    TH1D* h_3 = (TH1D*) h_1->Clone();
    TH1D* h_4 = (TH1D*) h_1->Clone();
    TH1D* h_5 = (TH1D*) h_1->Clone();
    TH1D* h_3pi = (TH1D*) h_1->Clone();

    TH1D* gamma = (TH1D*) h_1->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        h_1->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_2->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_3->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_4->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_5->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_3pi->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i))+couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i))+couplings[2]*gamma_tot[2]->Eval(gamma->GetXaxis()->GetBinCenter(i)));

    }  
    
    
    gamma->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma->SetLineColor(kBlack);
    gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma->Draw("histc");
    h_1->SetLineColor(kBlue);
    h_1->SetLineStyle(kDashed);
    //h_1->Draw("same");
    h_2->SetLineColor(kBlack);
    h_2->SetLineStyle(kDashed);
    //h_2->Draw("same");
    h_3->SetLineColor(kMagenta);
    h_3->SetLineStyle(kDashed);
    //h_3->Draw("same");  
    h_3pi->SetLineColor(kBlue);
    h_3pi->SetLineStyle(kDashed);
    h_3pi->Draw("same");
    h_4->SetLineColor(kRed);
    h_4->SetLineStyle(kDotted);
    h_4->Draw("same");  
    h_5->SetLineColor(kGreen);
    h_5->SetLineStyle(10);
    h_5->Draw("same");  
    
    
    /*
     BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
     TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
     TH1D* h_bw = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw->SetLineColor(kBlack);
     h_bw->SetLineStyle(kDashed);
     h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
     h_bw->Draw("same");
     
     */
    
    /*
     TF1* const_bw = new TF1("const_bw","1",0,max_s_inGeV2,0);
     vector<TF1*> gamma_const;
     gamma_const.push_back(const_bw);
     
     BW_plotter bw_plotter_const(vector<double>(1),gamma_const,m0,gamma0);
     TF1 *bw_const = new TF1("bw_const",bw_plotter_const,0,max_s_inGeV2,0);
     TH1D* h_bw_const = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw_const->SetBinContent(i,bw_const->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw_const->SetLineColor(kBlack);
     h_bw_const->SetLineStyle(kDashed);
     h_bw_const->Scale(gamma->GetMaximum()*0.8/bw_const->Eval(m0*m0/(GeV*GeV)));
     h_bw_const->Draw("same");
     */
    
        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(a)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	//text->Draw();

    c->Print("plots/gamma_pi2_1670.eps");

    gamma->Draw("histc");
    c->Print("plots/runningWidth_pi2_1670.C");

    //c->Print("gamma_pi2_1670.eps");
    
    TFile* out = new TFile("RunningWidth_pi(2)(1670)+.root","RECREATE");
    gamma->Write();  
    out->Write(); 
    out->Close();
    
    /*
     
     //rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
     TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
     
     double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
     
     for (int i=1; i<= h_m->GetNbinsX(); i++) {
     h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
     }    
     
     h_m->Draw("");
     c->Print("m_K1s_1680.pdf");
     
     */
    
    return 0;
}

int makeHistosForPsiKpipi(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 10323, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Integrate over a1 Dalitz plot   
    FitAmpSum fas(pat);
    fas.print();   
    //TH1D* h_dalitz = rwcalc.makeHisto_dalitz(nBins, max_s_inGeV2, nIntegrationEvents, &fas);
    
    // Integrate over flat phase space   
    //TH1D* h_phaseSpace = rwcalc.makeHisto_phaseSpace(nBins, max_s_inGeV2);  
    
    // Calculate K1->rho K 3-body width 
    DecayTree dk_K1_rhoK = DecayTree(10323);
    dk_K1_rhoK.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_K1_rhoK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_rhoK);
    
    // Calculate K1->K*(892) pi 3-body width 
    DecayTree dk_K1_Kstarpi = DecayTree(10323);
    dk_K1_Kstarpi.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_K1_Kstarpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);

    dk_K1_Kstarpi.getVal().setL(2);  
    TF1* f_K1_Kstarpi_D = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);

    // Calculate K1->K0*(1430) pi 3-body width 
    DecayTree dk_K1_K0starpi = DecayTree(10323);
    dk_K1_K0starpi.addDgtr(211, 10311)->addDgtr(321,-211); 
    TF1* f_K1_K0starpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_K0starpi);
   
    // Calculate K1->f0(1370) pi 3-body width 
    DecayTree dk_K1_f0K = DecayTree(10323);
    dk_K1_f0K.addDgtr(321, 30221)->addDgtr(211,-211); 
    TF1* f_K1_f0K = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_f0K);
    
    // Calculate K1->omega pi 2-body width 
    DecayTree dk_K1_omegapi = DecayTree(10323);
    dk_K1_omegapi.addDgtr(223, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_omegapi));  
    TF1* f_K1_omegapi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_K1_omegapi);  

    // Add K1 contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_K1_rhoK);
    gamma_tot.push_back(f_K1_Kstarpi);
    gamma_tot.push_back(f_K1_omegapi);
    gamma_tot.push_back(f_K1_K0starpi);
    gamma_tot.push_back(f_K1_f0K);
    //gamma_tot.push_back(f_K1_Kstarpi_D);

    vector<double> BF;
    BF.push_back(0.42);
    BF.push_back(0.16);
    BF.push_back(0.11);  
    BF.push_back(0.28);
    BF.push_back(0.03);
    //BF.push_back(0.3);
    double m0 = 1270;
    double gamma0 = 90.;
    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
    couplings.push_back(0.624354);
    couplings.push_back(0.0713955);
    couplings.push_back(0.334468);
    couplings.push_back(6.78019);
    */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_K1_rhoK = new TH1D("RunningWidth K_1(1270)","; s (GeV^{2}) ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_K1_Kstarpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kstarpi_D = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_K0starpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_omegapi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_f0K = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kpipi = (TH1D*) h_K1_rhoK->Clone();

    TH1D* gamma_K1 = (TH1D*) h_K1_rhoK->Clone();
    
    for (int i=1; i<= gamma_K1->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i));
        }
        gamma_K1->SetBinContent(i,val);
        h_K1_rhoK->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_Kstarpi->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_omegapi->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_K0starpi->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_f0K->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
	h_K1_Kpipi->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
    }  

    
    gamma_K1->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma_K1->SetLineColor(kBlack);
    gamma_K1->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_K1->Draw("histc");
    h_K1_Kpipi->SetLineColor(kBlue);
    h_K1_Kpipi->SetLineStyle(kDashed);
    h_K1_Kpipi->Draw("histcsame");

    h_K1_rhoK->SetLineColor(kBlue);
    h_K1_rhoK->SetLineStyle(kDashed);
    //h_K1_rhoK->Draw("histcsame");
    h_K1_Kstarpi->SetLineColor(kGreen);
    h_K1_Kstarpi->SetLineStyle(kDashed);
    //h_K1_Kstarpi->Draw("same");
    //h_K1_Kstarpi_D->SetLineColor(kBlue);
    //h_K1_Kstarpi_D->SetLineStyle(kDashed);
    //h_K1_Kstarpi_D->Draw("same");
    h_K1_omegapi->SetLineColor(kRed);
    h_K1_omegapi->SetLineStyle(kDotted);
    h_K1_omegapi->Draw("histcsame");
    h_K1_K0starpi->SetLineColor(kBlack);
    h_K1_K0starpi->SetLineStyle(kDashed);
    //h_K1_K0starpi->Draw("histcsame");
    h_K1_f0K->SetLineColor(kYellow);
    h_K1_f0K->SetLineStyle(kDashed);
    //h_K1_f0K->Draw("histcsame");

        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(a)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_K1_1270.eps");

    gamma_K1->Draw("histc");
    c->Print("plots/runningWidth_K1_1270.C");
    
    TFile* f_K1 = new TFile("RunningWidth_K(1)(1270)+.root","RECREATE");
    gamma_K1->Write();  
    f_K1->Write(); 
    f_K1->Close();
    
    /*
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_K1_1270.pdf");  
    */
    return 0;
}

int make_K1_1270_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 10323, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Integrate over a1 Dalitz plot   
    FitAmpSum fas(pat);
    fas.print();   
    //TH1D* h_dalitz = rwcalc.makeHisto_dalitz(nBins, max_s_inGeV2, nIntegrationEvents, &fas);
    
    // Integrate over flat phase space   
    //TH1D* h_phaseSpace = rwcalc.makeHisto_phaseSpace(nBins, max_s_inGeV2);  

    // Calculate K1->rho K 3-body width 
    DecayTree dk_K1_rhoK = DecayTree(10323);
    dk_K1_rhoK.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_K1_rhoK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_rhoK);
    
    DecayTree dk_K1_rhoPrimeK = DecayTree(10323);
    dk_K1_rhoPrimeK.addDgtr(321, 100113)->addDgtr(211,-211);  
    TF1* f_K1_rhoPrimeK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_rhoPrimeK);
    TH1D* h_1 = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_K1_rhoPrimeK,"RunningWidth_K1_rhoPrimeK.root");
    
    // Calculate K1->K*(892) pi 3-body width 
    DecayTree dk_K1_Kstarpi = DecayTree(10323);
    dk_K1_Kstarpi.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_K1_Kstarpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);

    dk_K1_Kstarpi.getVal().setL(2);  
    TF1* f_K1_Kstarpi_D = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);

    // Calculate K1->K0*(1430) pi 3-body width 
    DecayTree dk_K1_K0starpi = DecayTree(10323);
    dk_K1_K0starpi.addDgtr(211, 10311)->addDgtr(321,-211); 
    TF1* f_K1_K0starpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_K0starpi);
    TH1D* h_2 = rwcalc.makeHisto_3body(nBins, max_s_inGeV2, dk_K1_K0starpi,"RunningWidth_K1_K0starpi.root");

    // Calculate K1->f0(1370) pi 3-body width 
    DecayTree dk_K1_f0K = DecayTree(10323);
    dk_K1_f0K.addDgtr(321, 30221)->addDgtr(211,-211); 
    TF1* f_K1_f0K = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_f0K);
    
    // Calculate K1->omega pi 2-body width 
    DecayTree dk_K1_omegapi = DecayTree(10323);
    dk_K1_omegapi.addDgtr(223, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_omegapi));  
    TF1* f_K1_omegapi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_K1_omegapi);  

    // Add K1 contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_K1_rhoK);
    gamma_tot.push_back(f_K1_rhoPrimeK);
    //gamma_tot.push_back(f_K1_omegapi);
    gamma_tot.push_back(f_K1_Kstarpi);
    gamma_tot.push_back(f_K1_Kstarpi_D);
    gamma_tot.push_back(f_K1_K0starpi);
    //gamma_tot.push_back(f_K1_f0K);

    vector<double> BF;
    
    // K3pi LHCb
    //BF.push_back(0.963);
    //BF.push_back(0.4909);
    //BF.push_back(0.2708);
    //BF.push_back(0.0347);
    //BF.push_back(0.229);
    
    // K3pi LHCb alt
    //BF.push_back(0.5066);
    //BF.push_back(0.2525);
    //BF.push_back(0.0273);
    //BF.push_back(0.0597);
    //double m0 = 1285.03;
    //double gamma0 = 90.79;
    
    // KKpipi LHCb
/*    BF.push_back(0.4791);
    BF.push_back(0.0024);
    BF.push_back(0.3627);
    BF.push_back(0.0205);
    BF.push_back(0.0803);
    double m0 = 1289.81;
    double gamma0 = 116.11;*/
    
    // sys
    BF.push_back(1.0);
    double m0 = 1289.81;
    double gamma0 = 116.11;


    cout << "Have BF = " << BF << endl;
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
    couplings.push_back(0.624354);
    couplings.push_back(0.0713955);
    couplings.push_back(0.334468);
    couplings.push_back(6.78019);
    */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_K1_rhoK = new TH1D("RunningWidth","; s (GeV^{2}) ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_K1_rhoPrimeK = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kstarpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kstarpi_D = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_K0starpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_omegapi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_f0K = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kpipi = (TH1D*) h_K1_rhoK->Clone();

    TH1D* gamma_K1 = (TH1D*) h_K1_rhoK->Clone();
    
    for (int i=1; i<= gamma_K1->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i));
        }
        gamma_K1->SetBinContent(i,val);
/*        h_K1_rhoK->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_rhoPrimeK->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_Kstarpi->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_Kstarpi_D->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        //h_K1_omegapi->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_K0starpi->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));*/
        //h_K1_f0K->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
	
        /*h_K1_Kpipi->SetBinContent(i,
		couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));*/
    }  

    
    gamma_K1->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma_K1->SetLineColor(kBlack);
    gamma_K1->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_K1->Draw("histc");
    
    h_K1_Kpipi->SetLineColor(kBlue);
    h_K1_Kpipi->SetLineStyle(kDashed);
    //h_K1_Kpipi->Draw("histcsame");

    h_K1_rhoK->SetLineColor(kBlue);
    h_K1_rhoK->SetLineStyle(kDashed);
    //h_K1_rhoK->Draw("histcsame");
    h_K1_Kstarpi->SetLineColor(kGreen);
    h_K1_Kstarpi->SetLineStyle(kDashed);
    //h_K1_Kstarpi->Draw("same");
    //h_K1_Kstarpi_D->SetLineColor(kBlue);
    //h_K1_Kstarpi_D->SetLineStyle(kDashed);
    //h_K1_Kstarpi_D->Draw("same");
    h_K1_omegapi->SetLineColor(kRed);
    h_K1_omegapi->SetLineStyle(kDotted);
    //h_K1_omegapi->Draw("histcsame");
    h_K1_K0starpi->SetLineColor(kBlack);
    h_K1_K0starpi->SetLineStyle(kDashed);
    //h_K1_K0starpi->Draw("histcsame");
    h_K1_f0K->SetLineColor(kYellow);
    h_K1_f0K->SetLineStyle(kDashed);
    //h_K1_f0K->Draw("histcsame");

//         TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
//    	text->AddText("(a)");
//    	text->SetLineColor(kWhite);
//    	text->SetFillColor(kWhite);
//   	text->SetShadowColor(0);
//    	text->SetTextSize(0.075);
//    	text->SetTextFont(40);
//    	text->SetTextColor(kBlack);	
// 	text->Draw();

    c->Print("plots/gamma_K1_1270_alt.eps");

    gamma_K1->Draw("histc");
    c->Print("plots/runningWidth_K1_1270_alt.C");
    
    TFile* f_K1 = new TFile("RunningWidth_K(1)(1270)+_alt.root","RECREATE");
    gamma_K1->Write();  
    f_K1->Write(); 
    f_K1->Close();
    
    return 0;
}

int make_K1_1400_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 20213, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4.);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Integrate over a1 Dalitz plot   
    //FitAmpSum fas(pat);
    //fas.print();   
    //TH1D* h_dalitz = rwcalc.makeHisto_dalitz(nBins, max_s_inGeV2, nIntegrationEvents, &fas);
    
    // Integrate over flat phase space   
    //TH1D* h_phaseSpace = rwcalc.makeHisto_phaseSpace(nBins, max_s_inGeV2);  
    
    // Calculate K1->rho K 3-body width 
    DecayTree dk_K1_rhoK = DecayTree(20213);
    dk_K1_rhoK.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_K1_rhoK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_rhoK);
    
    // Calculate K1->K*(892) pi 3-body width 
    DecayTree dk_K1_Kstarpi = DecayTree(20213);
    dk_K1_Kstarpi.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_K1_Kstarpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);
    
    dk_K1_Kstarpi.getVal().setL(2);  
    TF1* f_K1_Kstarpi_D = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);
    
    // Calculate K1->K0*(1430) pi 3-body width 
    DecayTree dk_K1_K0starpi = DecayTree(20213);
    dk_K1_K0starpi.addDgtr(211, 10311)->addDgtr(321,-211); 
    TF1* f_K1_K0starpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_K0starpi);
    
    // Calculate K1->f0(1370) pi 3-body width 
    DecayTree dk_K1_f0K = DecayTree(20213);
    dk_K1_f0K.addDgtr(321, 30221)->addDgtr(211,-211); 
    TF1* f_K1_f0K = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_f0K);
    
    // Calculate K1->omega pi 2-body width 
    DecayTree dk_K1_omegapi = DecayTree(20213);
    dk_K1_omegapi.addDgtr(223, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_omegapi));  
    TF1* f_K1_omegapi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_K1_omegapi);  
    
    // Calculate K1->phi pi 2-body width 
    DecayTree dk_K1_phipi = DecayTree(20213);
    dk_K1_phipi.addDgtr(321, 333); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_phipi));  
    TF1* f_K1_phipi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_K1_phipi);  
    
    // Add K1 contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_K1_rhoK);
    gamma_tot.push_back(f_K1_Kstarpi);
    //gamma_tot.push_back(f_K1_Kstarpi_D);
    gamma_tot.push_back(f_K1_omegapi);
    //gamma_tot.push_back(f_K1_K0starpi);
    gamma_tot.push_back(f_K1_f0K);
    //gamma_tot.push_back(f_K1_phipi);

    
    vector<double> BF;
    // CLEO
//     BF.push_back(0.03);
//     BF.push_back(0.94);
//     //BF.push_back(0.3);
//     BF.push_back(0.01);  
//     BF.push_back(0.02);
//     //BF.push_back(0.06);
//     double m0 = 1403;
//     double gamma0 = 174.;

    // sys
    BF.push_back(1.);
    double m0 = 1403;
    double gamma0 = 174.;

    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //vector<double> couplings;
     //couplings.push_back(0.418994);
     //couplings.push_back(0.385653);
     //couplings.push_back(0.334468);
     //couplings.push_back(6.78019);
     
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_K1_rhoK = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_K1_Kstarpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kstarpi_D = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_K0starpi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_omegapi = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_f0K = (TH1D*) h_K1_rhoK->Clone();
    TH1D* h_K1_Kpipi = (TH1D*) h_K1_rhoK->Clone();
    
    TH1D* gamma_K1 = (TH1D*) h_K1_rhoK->Clone();
    
    for (int i=1; i<= gamma_K1->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i));
        }
        gamma_K1->SetBinContent(i,val);
        h_K1_rhoK->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_Kstarpi->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        //h_K1_Kstarpi_D->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_omegapi->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        //h_K1_K0starpi->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
        h_K1_f0K->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
	h_K1_Kpipi->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[1]*gamma_tot[1]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		//+couplings[4]*gamma_tot[4]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i))
		+couplings[3]*gamma_tot[3]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
    }  
    
    
    gamma_K1->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma_K1->SetLineColor(kBlack);
    gamma_K1->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_K1->Draw("histc");

    h_K1_Kpipi->SetLineColor(kBlue);
    h_K1_Kpipi->SetLineStyle(kDashed);
    h_K1_Kpipi->Draw("same");

    h_K1_rhoK->SetLineColor(kBlue);
    h_K1_rhoK->SetLineStyle(kDashed);
    //h_K1_rhoK->Draw("same");
    h_K1_Kstarpi->SetLineColor(kBlue);
    h_K1_Kstarpi->SetLineStyle(kDashed);
    //h_K1_Kstarpi->Draw("same");
    h_K1_omegapi->SetLineColor(kRed);
    h_K1_omegapi->SetLineStyle(kDotted);
    h_K1_omegapi->Draw("same");
    //h_K1_K0starpi->SetLineColor(kBlue);
    //h_K1_K0starpi->SetLineStyle(kDashed);
    //h_K1_K0starpi->Draw("same");
    h_K1_f0K->SetLineColor(kBlue);
    h_K1_f0K->SetLineStyle(kDashed);
    //h_K1_f0K->Draw("same");
    
    /*
    BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
    TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
    TH1D* h_bw = (TH1D*) h_K1_rhoK->Clone();
    for (int i=1; i<= gamma_K1->GetNbinsX(); i++) {
        h_bw->SetBinContent(i,bw->Eval(gamma_K1->GetXaxis()->GetBinCenter(i)));
    }
    h_bw->SetLineColor(kBlack);
    h_bw->SetLineStyle(kDashed);
    h_bw->Scale(gamma_K1->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
    h_bw->Draw("same");
    */    

//         TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
//    	text->AddText("(b)");
//    	text->SetLineColor(kWhite);
//    	text->SetFillColor(kWhite);
//   	text->SetShadowColor(0);
//    	text->SetTextSize(0.075);
//    	text->SetTextFont(40);
//    	text->SetTextColor(kBlack);	
// 	text->Draw();

    c->Print("plots/gamma_K1_1400_alt.eps");

    gamma_K1->Draw("histc");
    c->Print("plots/runningWidth_K1_1400_alt.C");    
    TFile* f_K1 = new TFile("RunningWidth_K(1)(1400)+_alt.root","RECREATE");
    gamma_K1->Write();  
    f_K1->Write(); 
    f_K1->Close();
    
    return 0;

    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV + m0/GeV*m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_K1_1400.pdf");  
    
    return 0;
}

int make_K_1460_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 20213, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4.);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate K1->K*(892) pi 3-body width 
    DecayTree dk_K1_Kstarpi = DecayTree(20213);
    dk_K1_Kstarpi.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_K1_Kstarpi = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_Kstarpi);
        
    // Calculate K1->f0(1370) pi 3-body width 
    DecayTree dk_K1_sigmaK = DecayTree(20213);
    dk_K1_sigmaK.addDgtr(321, 999001)->addDgtr(211,-211); 
    TF1* f_K1_sigmaK = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_K1_sigmaK);
   
    // Add K1 contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_K1_Kstarpi);
    gamma_tot.push_back(f_K1_sigmaK);
    
    vector<double> BF;
//     BF.push_back(0.5139);
//     BF.push_back(0.3123);

    BF.push_back(0.999);
    BF.push_back(0.001);
   
    double m0 = 1482.4;
    double gamma0 = 335.6;
    
    cout << "Have BF = " << BF << endl;
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //vector<double> couplings;
     //couplings.push_back(0.418994);
     //couplings.push_back(0.385653);
     //couplings.push_back(0.334468);
     //couplings.push_back(6.78019);
     
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* gamma_K1 = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    
    for (int i=1; i<= gamma_K1->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_K1->GetXaxis()->GetBinCenter(i));
        }
        gamma_K1->SetBinContent(i,val);
    }  
    
    
    gamma_K1->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma_K1->SetLineColor(kBlack);
    gamma_K1->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_K1->Draw("histc");
    c->Print("plots/gamma_K_1460_alt.eps");
    c->Print("plots/runningWidth_K_1460_alt.C");    
    
    TFile* f_K1 = new TFile("RunningWidth_K(1460)+_alt.root","RECREATE");
    gamma_K1->Write();  
    f_K1->Write(); 
    f_K1->Close();
    
    return 0;
}

int make_Ks_1410_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 100323, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate Ks->rho K 3-body width 
    DecayTree dk_1 = DecayTree(100323);
    dk_1.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_1 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_1);
    
    // Calculate Ks->K*(892) pi 3-body width 
    DecayTree dk_2 = DecayTree(100323);
    dk_2.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_2 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_2);
    
    // Calculate Ks->K pi 2-body width 
    DecayTree dk_3 = DecayTree(100323);
    dk_3.addDgtr(111, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_3));  
    TF1* f_3 = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_3);  
    
    // Add contributions    
    vector<TF1*> gamma_tot;
//     gamma_tot.push_back(f_1);
    gamma_tot.push_back(f_2);
//     gamma_tot.push_back(f_3);
    
    vector<double> BF;
    // CLEO
//     BF.push_back(0.034);
//     BF.push_back(0.80);
//     BF.push_back(0.066);  
//     double m0 = 1414;
//     double gamma0 = 232.;

    // 2. it
//     BF.push_back(0.0585/(0.0585+0.146)*(1.-0.066));
//     BF.push_back(0.146/(0.0585+0.146)*(1.-0.066));
//     BF.push_back(0.066);  
//     double m0 = 1436.61;
//     double gamma0 = 386.654;

    // sys
    BF.push_back(1.);
    double m0 = 1436.61;
    double gamma0 = 386.654;

    
    cout << "Have BF = " << BF << endl;    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_1 = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_2 = (TH1D*) h_1->Clone();
    TH1D* h_3 = (TH1D*) h_1->Clone();
    TH1D* h_4 = (TH1D*) h_1->Clone();
    
    TH1D* gamma = (TH1D*) h_1->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
//         h_1->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
//         h_2->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
//         h_3->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
//         h_4->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i))+couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }  
    
    
    gamma->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma->SetLineColor(kBlack);
    gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma->Draw("histc");
    h_4->SetLineColor(kBlue);
    h_4->SetLineStyle(kDashed);
    h_4->Draw("histcsame");
    h_2->SetLineColor(kBlack);
    h_2->SetLineStyle(kDashed);
    //h_2->Draw("same");
    h_3->SetLineColor(kRed);
    h_3->SetLineStyle(kDotted);
    h_3->Draw("histcsame");    
    
    
    /*
    BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
    TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
    TH1D* h_bw = (TH1D*) h_1->Clone();
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }
    h_bw->SetLineColor(kBlack);
    h_bw->SetLineStyle(kDashed);
    h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
    h_bw->Draw("same");
    
    */
    
    /*
     TF1* const_bw = new TF1("const_bw","1",0,max_s_inGeV2,0);
     vector<TF1*> gamma_const;
     gamma_const.push_back(const_bw);
     
     BW_plotter bw_plotter_const(vector<double>(1),gamma_const,m0,gamma0);
     TF1 *bw_const = new TF1("bw_const",bw_plotter_const,0,max_s_inGeV2,0);
     TH1D* h_bw_const = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw_const->SetBinContent(i,bw_const->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw_const->SetLineColor(kBlack);
     h_bw_const->SetLineStyle(kDashed);
     h_bw_const->Scale(gamma->GetMaximum()*0.8/bw_const->Eval(m0*m0/(GeV*GeV)));
     h_bw_const->Draw("same");
     */
    
//         TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
//    	text->AddText("(a)");
//    	text->SetLineColor(kWhite);
//    	text->SetFillColor(kWhite);
//   	text->SetShadowColor(0);
//    	text->SetTextSize(0.075);
//    	text->SetTextFont(40);
//    	text->SetTextColor(kBlack);	
// 	text->Draw();

    c->Print("plots/gamma_Ks_1410_alt.eps");

    gamma->Draw("histc");
    c->Print("plots/runningWidth_Ks_1410_alt.C"); 
    
    TFile* out = new TFile("RunningWidth_K*(1410)+_alt.root","RECREATE");
    gamma->Write();  
    out->Write(); 
    out->Close();
    
    /*
    
    //rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_K1s_1680.pdf");
    
    */
    
    return 0;
}

int make_K2s_1430_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 325, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",5);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate K2->rho K 3-body width 
    DecayTree dk_1 = DecayTree(325);
    dk_1.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_1 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_1);
    
    // Calculate K2->K*(892) pi 3-body width 
    DecayTree dk_2 = DecayTree(325);
    dk_2.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_2 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_2);
    
    // Calculate K2->K pi 2-body width 
    DecayTree dk_3 = DecayTree(325);
    dk_3.addDgtr(111, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_3));  
    TF1* f_3 = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_3);  
    
    DecayTree dk_4 = DecayTree(325);
    dk_4.addDgtr(111, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_4));  
    TF1* f_4= rwcalc.getRunningWidthFunction_phaseSpace(max_s_inGeV2);  

    
    // Calculate K2->K omega 2-body width 
    DecayTree dk_5 = DecayTree(325);
    dk_5.addDgtr(313, 211, -211); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_5));  
    TF1* f_5 = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_5); 
    
    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_1);
    gamma_tot.push_back(f_2);
    gamma_tot.push_back(f_3);
    gamma_tot.push_back(f_4);
    gamma_tot.push_back(f_5);
    
    vector<double> BF;
    BF.push_back(0.087/(0.087 +0.247 + 0.499 + 0.134 + 0.029));
    BF.push_back(0.247/(0.087 +0.247 + 0.499 + 0.134 + 0.029));
    BF.push_back(0.499/(0.087 +0.247 + 0.499 + 0.134 + 0.029));  
    BF.push_back(0.134/(0.087 +0.247 + 0.499 + 0.134 + 0.029));  
    BF.push_back(0.029/(0.087 +0.247 + 0.499 + 0.134 + 0.029));  
    double m0 = 1425.6;
    double gamma0 = 98.5;
    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_1 = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_2 = (TH1D*) h_1->Clone();
    TH1D* h_3 = (TH1D*) h_1->Clone();
    TH1D* h_4 = (TH1D*) h_1->Clone();
    TH1D* h_5 = (TH1D*) h_1->Clone();
    
    TH1D* gamma = (TH1D*) h_1->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        h_1->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_2->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_3->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_4->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_5->SetBinContent(i,couplings[4]*gamma_tot[4]->Eval(gamma->GetXaxis()->GetBinCenter(i)));

    }  
    
    
    gamma->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma->SetLineColor(kRed);
    gamma->SetTitle("; s [GeV^{2}] ; #Gamma(s) [GeV]");  
    gamma->Draw("");
    h_1->SetLineColor(kBlue);
    h_1->SetLineStyle(kDashed);
    h_1->Draw("same");
    h_2->SetLineColor(kBlack);
    h_2->SetLineStyle(kDashed);
    h_2->Draw("same");
    h_3->SetLineColor(kMagenta);
    h_3->SetLineStyle(kDashed);
    h_3->Draw("same");  
    h_4->SetLineColor(kGreen);
    h_4->SetLineStyle(kDashed);
    h_4->Draw("same");  
    h_5->SetLineColor(kYellow);
    h_5->SetLineStyle(kDashed);
    h_5->Draw("same");  
    
    
    /*
     BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
     TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
     TH1D* h_bw = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw->SetLineColor(kBlack);
     h_bw->SetLineStyle(kDashed);
     h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
     h_bw->Draw("same");
     
     */
    
    /*
     TF1* const_bw = new TF1("const_bw","1",0,max_s_inGeV2,0);
     vector<TF1*> gamma_const;
     gamma_const.push_back(const_bw);
     
     BW_plotter bw_plotter_const(vector<double>(1),gamma_const,m0,gamma0);
     TF1 *bw_const = new TF1("bw_const",bw_plotter_const,0,max_s_inGeV2,0);
     TH1D* h_bw_const = (TH1D*) h_1->Clone();
     for (int i=1; i<= gamma->GetNbinsX(); i++) {
     h_bw_const->SetBinContent(i,bw_const->Eval(gamma->GetXaxis()->GetBinCenter(i)));
     }
     h_bw_const->SetLineColor(kBlack);
     h_bw_const->SetLineStyle(kDashed);
     h_bw_const->Scale(gamma->GetMaximum()*0.8/bw_const->Eval(m0*m0/(GeV*GeV)));
     h_bw_const->Draw("same");
     */
    
    c->Print("gamma_K2s_1430.pdf");
    
    TFile* out = new TFile("RunningWidth_K(2)*(1430)+.root","RECREATE");
    gamma->Write();  
    out->Write(); 
    out->Close();
    
    /*
     
     //rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
     TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
     
     double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
     
     for (int i=1; i<= h_m->GetNbinsX(); i++) {
     h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
     }    
     
     h_m->Draw("");
     c->Print("m_K1s_1680.pdf");
     
     */
    
    return 0;
}

int make_Ks_1680_Histo(){  
    TCanvas* c = new TCanvas();  
    
    NamedParameter<int> EventPattern("Event Pattern", 30323, 321, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<int> nIntegrationEvents("nIntegrationEvents",500000);  
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4);  
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate Ks->rho K 3-body width 
    DecayTree dk_1 = DecayTree(30323);
    dk_1.addDgtr(321, 113)->addDgtr(211,-211);  
    TF1* f_1 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_1);
    
    // Calculate Ks->K*(892) pi 3-body width 
    DecayTree dk_2 = DecayTree(30323);
    dk_2.addDgtr(211, 323)->addDgtr(321,-211); 
    TF1* f_2 = rwcalc.getRunningWidthFunction_3body(max_s_inGeV2, dk_2);
    
    // Calculate Ks->K pi 2-body width 
    DecayTree dk_3 = DecayTree(30323);
    dk_3.addDgtr(111, 321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_3));  
    TF1* f_3 = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_3);  
    
    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_1);
    gamma_tot.push_back(f_2);
    gamma_tot.push_back(f_3);
    
    vector<double> BF;
    BF.push_back(0.314);
    BF.push_back(0.299);
    BF.push_back(0.387);  
    double m0 = 1717;
    double gamma0 = 322.;
    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    /*vector<double> couplings;
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     couplings.push_back(0.2);
     */
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_1 = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_2 = (TH1D*) h_1->Clone();
    TH1D* h_3 = (TH1D*) h_1->Clone();
    TH1D* h_4 = (TH1D*) h_1->Clone();
    
    TH1D* gamma = (TH1D*) h_1->Clone();
    
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma->GetXaxis()->GetBinCenter(i));
        }
        gamma->SetBinContent(i,val);
        h_1->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_2->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
        h_3->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
 	h_4->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma->GetXaxis()->GetBinCenter(i))+couplings[1]*gamma_tot[1]->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }  
    
    
    gamma->SetMinimum(0.);
    //gamma_K1->SetMaximum(1.7);
    gamma->SetLineColor(kBlack);
    gamma->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma->Draw("histc");
    h_4->SetLineColor(kBlue);
    h_4->SetLineStyle(kDashed);
    h_4->Draw("histcsame");
    h_2->SetLineColor(kBlack);
    h_2->SetLineStyle(kDashed);
    //h_2->Draw("same");
    h_3->SetLineColor(kRed);
    h_3->SetLineStyle(kDotted);
    h_3->Draw("histcsame");    
    
    /*
    BW_plotter bw_plotter(couplings,gamma_tot,m0,gamma0);
    TF1 *bw = new TF1("bw",bw_plotter,0,max_s_inGeV2,0);
    TH1D* h_bw = (TH1D*) h_1->Clone();
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        h_bw->SetBinContent(i,bw->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }
    h_bw->SetLineColor(kBlack);
    h_bw->SetLineStyle(kDashed);
    h_bw->Scale(gamma->GetMaximum()*0.8/bw->Eval(m0*m0/(GeV*GeV)));
    h_bw->Draw("same");
    */

    /*
    TF1* const_bw = new TF1("const_bw","1",0,max_s_inGeV2,0);
    vector<TF1*> gamma_const;
    gamma_const.push_back(const_bw);
    
    BW_plotter bw_plotter_const(vector<double>(1),gamma_const,m0,gamma0);
    TF1 *bw_const = new TF1("bw_const",bw_plotter_const,0,max_s_inGeV2,0);
    TH1D* h_bw_const = (TH1D*) h_1->Clone();
    for (int i=1; i<= gamma->GetNbinsX(); i++) {
        h_bw_const->SetBinContent(i,bw_const->Eval(gamma->GetXaxis()->GetBinCenter(i)));
    }
    h_bw_const->SetLineColor(kBlack);
    h_bw_const->SetLineStyle(kDashed);
    h_bw_const->Scale(gamma->GetMaximum()*0.8/bw_const->Eval(m0*m0/(GeV*GeV)));
    h_bw_const->Draw("same");
    */
    
        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(b)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_Ks_1680.eps");

    gamma->Draw("histc");
    c->Print("plots/runningWidth_Ks_1680.C"); 

    
    TFile* out = new TFile("RunningWidth_K*(1680)+.root","RECREATE");
    gamma->Write();  
    out->Write(); 
    out->Close();
    
    return 0;

    //rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_K1_rhoK));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_K1s_1680.pdf");
    
    
    return 0;
}

int make_f2_histo(){      
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",3.);  
    NamedParameter<int> EventPattern("Event Pattern", 225, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;

    RunningWidthCalculator rwcalc(pat);
    
    // Calculate f2-> 2pi width 
    DecayTree dk_f2_2pi = DecayTree(225);
    dk_f2_2pi.addDgtr(211, -211); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f2_2pi));
    TF1* f_f2_2pi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f2_2pi); 
    
    // Calculate f2-> 4pi fermi width  
    TF1* f_f2_4pi = rwcalc.getRunningWidthFunction_Fermi(max_s_inGeV2,3.39, 3.238);
    
    // Calculate f2-> 2K width 
    DecayTree dk_f2_2K = DecayTree(225);
    dk_f2_2K.addDgtr(321, -321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f2_2K));
    TF1* f_f2_2K = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f2_2K); 
    
    // Calculate f2-> 2eta width 
    DecayTree dk_f2_2eta = DecayTree(225);
    dk_f2_2eta.addDgtr(221, 221); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f2_2eta));
    TF1* f_f2_2eta = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f2_2eta); 
    
    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_f2_2pi);
    gamma_tot.push_back(f_f2_4pi);
    gamma_tot.push_back(f_f2_2K);
    gamma_tot.push_back(f_f2_2eta);

    vector<double> BF;
    BF.push_back(0.842);
    BF.push_back(0.108);
    BF.push_back(0.046);
    BF.push_back(0.004);

    double m0 = 1275.5;
    double gamma0 = 186.7;
    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //couplings.push_back(0.457455);
    //couplings.push_back(0.674893);
    //couplings.push_back(0.100314);
    //couplings.push_back(0.00872282);
    
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_f2_2pi = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_f2_4pi = (TH1D*) h_f2_2pi->Clone();
    TH1D* h_f2_2K = (TH1D*) h_f2_2pi->Clone();
    TH1D* h_f2_2eta = (TH1D*) h_f2_2pi->Clone();
    TH1D* h_f2_2K2eta = (TH1D*) h_f2_2pi->Clone();

    TH1D* gamma_f2 = (TH1D*) h_f2_2pi->Clone();
    
    for (int i=1; i<= gamma_f2->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i));
        }
        gamma_f2->SetBinContent(i,val);
        h_f2_2pi->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i)));
        h_f2_4pi->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i)));
        h_f2_2K->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i)));
        h_f2_2eta->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i)));
        h_f2_2K2eta->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i))+couplings[3]*gamma_tot[3]->Eval(gamma_f2->GetXaxis()->GetBinCenter(i)));

    }  
    
    TCanvas* c = new TCanvas();  
    gamma_f2->SetMinimum(0.);
    gamma_f2->SetLineColor(kBlack);
    gamma_f2->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_f2->Draw("");
    h_f2_2pi->SetLineColor(kBlue);
    h_f2_2pi->SetLineStyle(kDashed);
    h_f2_2pi->Draw("same");
    h_f2_4pi->SetLineColor(kRed);
    h_f2_4pi->SetLineStyle(kDotted);
    h_f2_4pi->Draw("same");
    h_f2_2K->SetLineColor(kGreen);
    h_f2_2K->SetLineStyle(kDashed);
    //h_f2_2K->Draw("same");
    h_f2_2eta->SetLineColor(kBlue);
    h_f2_2eta->SetLineStyle(kDashed);
    //h_f2_2eta->Draw("same");
    h_f2_2K2eta->SetLineColor(kGreen);
    h_f2_2K2eta->SetLineStyle(10);
    h_f2_2K2eta->Draw("same");

        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(b)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_f2.eps");

    gamma_f2->Draw("");
    c->Print("plots/runningWidth_f2_1270.C");


    //c->Print("gamma_f2.eps");
    
    TFile* f_f2 = new TFile("RunningWidth_f(2)(1270)+.root","RECREATE");
    gamma_f2->Write();  
    f_f2->Write(); 
    f_f2->Close();

    return 0;

    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f2_2pi));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_f2_1270.eps");
    
    return 0;
}

void make_f0_1370_histo_Bugg_test(){

    NamedParameter<int> nBins("nBins",20);
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",4.);  
    NamedParameter<int> EventPattern("Event Pattern", 30221, 211, -211, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());

    RunningWidthCalculator rwcalc(pat);
    
    vector<double> BF;
    BF.push_back(1);
    
    double m0 = 1315;
    double gamma0 = 54;
    
    cout << "Have BF = " << BF << endl;
    
    TF1* f_f0_4pi = rwcalc.getRunningWidthFunction_Fermi(max_s_inGeV2,3.39, 3.238);
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_f0_4pi);
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //vector<double> couplings;  
    //couplings.push_back(1);
    
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings); 
    TH1D* h_m_Bugg = (TH1D*) h_m->Clone();
    TH1D* h_f0_4pi = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(h_m->GetXaxis()->GetBinCenter(i));
        }
        h_m->SetBinContent(i,-(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
        h_f0_4pi->SetBinContent(i,val);
        
        
        vector<double> a;
        a.push_back(3.5320);
        a.push_back(-0.0427);
        vector<double> s;
        s.push_back(2.9876);
        s.push_back(-0.4619);
        vector<double> w;
        w.push_back(0.8804);
        w.push_back(-0.0036);
        
        double bugg = 0.;
        
        for (int j=0; j<2; j++) {
            bugg += a[j]/( pow(h_m->GetXaxis()->GetBinCenter(i)-s[j],2) + w[j]*w[j]);
            bugg += -a[j]/( pow(m0*m0/GeV/GeV,2) + w[j]*w[j]);
        }
        
        h_m_Bugg->SetBinContent(i,bugg );

    }  
    
    cout << m0*m0/GeV/GeV<< endl;
    cout << couplings[0]*gamma_tot[0]->Eval(m0*m0/GeV/GeV) << endl;
    
    TCanvas* c = new TCanvas();  
    h_f0_4pi->Draw("");
    c->Print("gamma_f0.pdf");
    
    h_m->Draw("");
    h_m_Bugg->SetLineColor(kBlack);
    h_m_Bugg->Draw("same");
    c->Print("m_f0.pdf");

}

int make_f0_histo(){      
    NamedParameter<int> nBins("nBins",2000);
    NamedParameter<double> max_s_inGeV2("max_s_inGeV2",3);  
    NamedParameter<int> EventPattern("Event Pattern", 30221, 211, -211);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    RunningWidthCalculator rwcalc(pat);
    
    // Calculate f0-> 2pi width 
    DecayTree dk_f0_2pi = DecayTree(30221);
    dk_f0_2pi.addDgtr(211, -211); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f0_2pi));
    TF1* f_f0_2pi = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f0_2pi); 
    
    // Calculate f0-> 4pi fermi width  
    TF1* f_f0_4pi = rwcalc.getRunningWidthFunction_Fermi(max_s_inGeV2,3.39, 3.238);
    
    // Calculate f0-> 2K width 
    DecayTree dk_f0_2K = DecayTree(30221);
    dk_f0_2K.addDgtr(321, -321); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f0_2K));
    TF1* f_f0_2K = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f0_2K); 
    
    // Calculate f0-> 2eta width 
    DecayTree dk_f0_2eta = DecayTree(30221);
    dk_f0_2eta.addDgtr(221, 221); 
    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f0_2eta));
    TF1* f_f0_2eta = rwcalc.getRunningWidthFunction_2body(max_s_inGeV2, dk_f0_2eta); 
    
    // Add contributions    
    vector<TF1*> gamma_tot;
    gamma_tot.push_back(f_f0_2pi);
    gamma_tot.push_back(f_f0_4pi);
    gamma_tot.push_back(f_f0_2K);
    //gamma_tot.push_back(f_f0_2eta);
    
    vector<double> BF;
    //BF.push_back(1);
    BF.push_back(0.26);
    BF.push_back(0.62);
    BF.push_back(0.12);
    //BF.push_back(0.004);
    
    double m0 = 1475;
    double gamma0 = 113;
    
    cout << "Have BF = " << BF << endl;
    
    vector<double> couplings = rwcalc.getPartialWidthCouplingsFromBF(BF , gamma_tot, m0, gamma0);  
    //couplings.push_back(0.457455);
    //couplings.push_back(0.674893);
    //couplings.push_back(0.100314);
    //couplings.push_back(0.00872282);
    
    cout << "Couplings = " << couplings  << endl;
    
    TH1D* h_f0_2pi = new TH1D("RunningWidth","; s [GeV^{2}] ; #Gamma(s) [GeV]", nBins, 0., max_s_inGeV2);
    TH1D* h_f0_4pi = (TH1D*) h_f0_2pi->Clone();
    TH1D* h_f0_2K = (TH1D*) h_f0_2pi->Clone();
    TH1D* h_f0_2eta = (TH1D*) h_f0_2pi->Clone();
    
    TH1D* gamma_f0 = (TH1D*) h_f0_2pi->Clone();
    
    for (int i=1; i<= gamma_f0->GetNbinsX(); i++) {
        double val = 0.;
        for (int j=0; j<gamma_tot.size(); j++) {
            val += couplings[j]*gamma_tot[j]->Eval(gamma_f0->GetXaxis()->GetBinCenter(i));
        }
        gamma_f0->SetBinContent(i,val);
        h_f0_2pi->SetBinContent(i,couplings[0]*gamma_tot[0]->Eval(gamma_f0->GetXaxis()->GetBinCenter(i)));
        h_f0_4pi->SetBinContent(i,couplings[1]*gamma_tot[1]->Eval(gamma_f0->GetXaxis()->GetBinCenter(i)));
        h_f0_2K->SetBinContent(i,couplings[2]*gamma_tot[2]->Eval(gamma_f0->GetXaxis()->GetBinCenter(i)));
        //h_f0_2eta->SetBinContent(i,couplings[3]*gamma_tot[3]->Eval(gamma_f0->GetXaxis()->GetBinCenter(i)));
    }  
    
    TCanvas* c = new TCanvas();  
    gamma_f0->SetMinimum(0.);
    gamma_f0->SetLineColor(kBlack);
    gamma_f0->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)");  
    gamma_f0->Draw("");
    h_f0_2pi->SetLineColor(kBlue);
    h_f0_2pi->SetLineStyle(kDashed);
    h_f0_2pi->Draw("same");
    h_f0_4pi->SetLineColor(kRed);
    h_f0_4pi->SetLineStyle(kDotted);
    h_f0_4pi->Draw("same");
    h_f0_2K->SetLineColor(kGreen);
    h_f0_2K->SetLineStyle(10);
    h_f0_2K->Draw("same");
    //h_f0_2eta->SetLineColor(kBlue);
    //h_f0_2eta->SetLineStyle(kDashed);
    //h_f0_2eta->Draw("same");
    //c->Print("gamma_f0_1370.eps");
    
        TPaveText *text= new TPaveText(0.17,0.85,0.3,0.9,"NDC");
   	text->AddText("(a)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();

    c->Print("plots/gamma_f0_1370.eps");

    gamma_f0->Draw("");
    c->Print("plots/runningWidth_f0_1370.C");

    TFile* f_f0 = new TFile("RunningWidth_f(0)(1370)+.root","RECREATE");
    gamma_f0->Write();  
    f_f0->Write(); 
    f_f0->Close();
    
    return 0;

    rwcalc.setDalitzEventPattern(DalitzEventPattern(dk_f0_2pi));  
    TH1D* h_m = rwcalc.makeRunningMassHisto_3body(nBins, max_s_inGeV2, gamma_tot, couplings);  
    
    double delta_m0 = h_m->Interpolate(m0*m0/GeV/GeV);
    
    for (int i=1; i<= h_m->GetNbinsX(); i++) {
        h_m->SetBinContent(i,(h_m->GetBinContent(i) - delta_m0) *m0/GeV );
    }    
    
    h_m->Draw("");
    c->Print("m_f0_1370.pdf");
    
    return 0;
}


void plotIterations(){
    TFile *f1 =  TFile::Open("RunningWidth_a(1)(1260)+_12it.root");
    TH1D* h1;
    h1=dynamic_cast<TH1D*>(f1->Get("RunningWidth"));

    TFile *f2 =  TFile::Open("RunningWidth_a(1)(1260)+_10it.root");
    TH1D* h2;
    h2=dynamic_cast<TH1D*>(f2->Get("RunningWidth"));

    TFile *f3 =  TFile::Open("RunningWidth_a(1)(1260)+_8it.root");
    TH1D* h3;
    h3=dynamic_cast<TH1D*>(f3->Get("RunningWidth"));

    TFile *f4 =  TFile::Open("RunningWidth_a1+_PhaseSpace.root");
    TH1D* h4;
    h4=dynamic_cast<TH1D*>(f4->Get("RunningWidth"));

    TFile *f5 =  TFile::Open("RunningWidth_a1_rhopi.root");
    TH1D* h5;
    h5=dynamic_cast<TH1D*>(f5->Get("RunningWidth"));

    TFile *f6 =  TFile::Open("RunningWidth_a1_rhopi_D.root");
    TH1D* h6;
    h6=dynamic_cast<TH1D*>(f6->Get("RunningWidth"));

    TFile *f7 =  TFile::Open("RunningMass_a(1)(1260)+_3it.root");
    //TFile *f7 =  TFile::Open("RunningMass_a1+_3body.root");
    //TFile *f7 =  TFile::Open("Histos_1/RunningMass_a(1)(1260)+.root");
    TH1D* h7;
    h7=dynamic_cast<TH1D*>(f7->Get("RunningMass"));

    double m0 = 1225.;
    double gamma0 = 430.;
    double norm = h4->Interpolate(m0*m0/GeV/GeV);
    double norm5 = h5->Interpolate(m0*m0/GeV/GeV);
    double norm6 = h6->Interpolate(m0*m0/GeV/GeV);

    for (int i=1; i<= h4->GetNbinsX(); i++) {
    	 h4->SetBinContent(i,h4->GetBinContent(i)/norm * gamma0/GeV );
    	 h5->SetBinContent(i,h5->GetBinContent(i)/norm5 * gamma0/GeV );
    	 h6->SetBinContent(i,h6->GetBinContent(i)/norm6 * gamma0/GeV );
    }    
     
    TCanvas *c = new TCanvas();

    h1->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(3);
    h1->Draw("histc");

    h2->SetLineColor(kBlue);
    h2->SetLineStyle(kDashed);
    h2->Draw("histcsame");

    h3->SetLineColor(kRed);
    h3->SetLineStyle(kDashed);
    h3->Draw("histcsame");

    h4->SetLineColor(kMagenta);
    h4->SetLineStyle(kDashed);
    h4->Draw("histcsame");

    //h1->SetLineColor(kMagenta);
    //h1->SetLineStyle(kDashed);
    //h1->Draw("histcsame");

    c->Print("rw_iterations.eps");

    h5->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
    h5->SetLineColor(kBlack);
    h5->SetLineWidth(3);
    h5->Draw("histc");
    c->Print("a1rho.eps");

    h6->SetTitle("; s (GeV^{2}) ; #Gamma(s) (GeV)"); 
    h6->SetLineColor(kBlack);
    h6->SetLineWidth(3);
    h6->Draw("histc");
    c->Print("a1rhoD.eps");


    m0 = 1225.;
    gamma0 = 430.;
    	double d_m02 = h7->Interpolate(m0*m0/GeV/GeV);

    for (int i=1; i<= h7->GetNbinsX(); i++) {
    	double d_m2  = h7->GetBinContent(i);

	h7->SetBinContent(i,  (m0*m0/GeV/GeV + m0/GeV*(d_m2 - d_m02) * gamma0/GeV));
	//h7->SetBinContent(i,  ((d_m2 - d_m02)/0.43));

    }    

    
    //TFile *f =  new TFile("RunningMass_a(1)(1260)+_4it.root","RECREATE");    
    //h7->Write();
    //f->Write();    

    h7->Draw("histc");
    
        TPaveText *text= new TPaveText(0.2,0.85,0.33,0.9,"NDC");
   	text->AddText("(b)");
   	text->SetLineColor(kWhite);
   	text->SetFillColor(kWhite);
  	text->SetShadowColor(0);
   	text->SetTextSize(0.075);
   	text->SetTextFont(40);
   	text->SetTextColor(kBlack);	
	text->Draw();
    //c->Print("plots/gamma_a1_1260.eps");

    //gamma->Draw("");

    c->Print("a1_rw.eps");
    c->Print("plots/runningWidth_a1_1260_new.C");

}






void plot(string nameA, string nameB, double m0, double g0, string out){
    TFile *f1 =  TFile::Open(nameA.c_str());
    TH1D* h1=dynamic_cast<TH1D*>(f1->Get("RunningWidth"));
    TFile *f2 =  TFile::Open(nameB.c_str());
    TH1D* h2 = dynamic_cast<TH1D*>(f2->Get("RunningWidth"));
     
    double norm1 = h1->Interpolate(m0*m0/GeV/GeV);
    double norm2 = h2->Interpolate(m0*m0/GeV/GeV);

    TH1D* h1_new= new TH1D("","; #sqrt{s} (GeV) ; #Gamma(s) (GeV)",100,0.7,2);
    TH1D* h2_new= new TH1D("","; #sqrt{s} (GeV) ; #Gamma(s) (GeV)",100,0.7,2);

    for (int i=1; i<= h1_new->GetNbinsX(); i++) {
    	 h1_new->SetBinContent(i,h1->Interpolate(pow(h1_new->GetBinCenter(i),2))/norm1 * g0/GeV );
    	 h2_new->SetBinContent(i,h2->Interpolate(pow(h2_new->GetBinCenter(i),2))/norm2 * g0/GeV );
    }


    TCanvas *c = new TCanvas();
    h1_new->SetLineColor(kBlue+1);
    h1_new->SetLineWidth(5);
    h1_new->Draw("histc");

    h2_new->SetLineColor(kRed+1);
    h2_new->SetLineStyle(kDashed);
    h2_new->SetLineWidth(3);
    h2_new->Draw("histcsame");

    c->Print((out+".eps").c_str());
}



int main(){
    time_t startTime = time(0);
    
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    gROOT->ProcessLine(".x lhcbStyle.C");
    

    plot("Histos_KKpipi_LHCb/RunningWidth_K(1)(1270)+.root","Histos_alt/RunningWidth_K(1)(1270)+_alt.root",1289.81,116.11,"rw_K1_1270");
    plot("Histos_CLEO/RunningWidth_K(1)(1400)+.root","Histos_alt/RunningWidth_K(1)(1400)+_alt.root",1397.87,204.377,"rw_K1_1400");
    plot("Histos_CLEO/RunningWidth_K*(1410+.root","Histos_alt/RunningWidth_K*(1410+_alt.root",1432.07,344.149,"rw_Ks");
    plot("Histos_K3pi_LHCb/RunningWidth_K(1460)+.root","Histos_alt/RunningWidth_K(1460)+_alt.root",1482.4,335.6 ,"rw_K");

    //plotIterations();
 
    //   lhcbStyle();
  
    /*
    gStyle->SetTitleXSize(0.05);
    gStyle->SetTitleYSize(0.05);
    gStyle->SetTitleFont(42,"X");
    gStyle->SetTitleFont(42,"Y");
    gStyle->SetLabelFont(42,"X");
    gStyle->SetLabelFont(42,"Y");
    gStyle->SetLabelOffset(0.01,"X");
    gStyle->SetTitleOffset(0.8,"X");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    */
    
    //makeHistosForD4pi();
    //make_a1_1640_Histo();
    //make_pi_1300_Histo();
    //make_pi2_1670_Histo();
    
    //makeHistosForPsiKpipi();
//     make_K1_1270_Histo();
//     make_K1_1400_Histo();
//     make_Ks_1410_Histo();
//     make_K_1460_Histo();
    //make_K2s_1430_Histo();
    //make_Ks_1680_Histo();
    //make_f2_histo();
    //make_f0_1370_histo_Bugg_test();
    //make_f0_histo();

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
//
