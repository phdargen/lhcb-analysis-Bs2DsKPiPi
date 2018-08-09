#define pull_cxx
#include "pull.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void pull::Loop(string parName)
{
//   In a ROOT session, you can do:
//      Root > .L pull.C
//      Root > pull t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//
   if (fChain == 0) return;

   TH1D* Plot_C_shift = new TH1D("Plot_C_shift","C_{init}-C_{fit}", 10, -0.5,0.5);
   TH1D* Plot_D_shift = new TH1D("Plot_D_shift","D_{init}-D_{fit}", 10, -0.5,0.5);
   TH1D* Plot_S_shift = new TH1D("Plot_S_shift","S_{init}-S_{fit}", 10, -0.5,0.5);
   TH1D* Plot_D_bar_shift = new TH1D("Plot_D_bar_shift","D_bar_{init}-D_bar_{fit}", 10, -0.5,0.5);
   TH1D* Plot_S_bar_shift = new TH1D("Plot_S_bar_shift","S_bar_{init}-S_bar_{fit}", 10, -0.5,0.5);

   TH1D* pullPlot_C = new TH1D("pullPlot_C","Pull C", 10, -3.,3.);
   TH1D* pullPlot_D = new TH1D("pullPlot_D","Pull D", 10, -3.,3.);
   TH1D* pullPlot_S = new TH1D("pullPlot_S","Pull S", 10, -3.,3.);
   TH1D* pullPlot_D_bar = new TH1D("pullPlot_D_bar","Pull D_bar", 10, -3.,3.);
   TH1D* pullPlot_S_bar = new TH1D("pullPlot_S_bar","Pull S_bar", 10, -3.,3.);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

	Plot_C_shift->Fill((C_init-C_mean));
	Plot_D_shift->Fill((D_init-D_mean));
	Plot_S_shift->Fill((S_init-S_mean));
	Plot_D_bar_shift->Fill((D_bar_init-D_bar_mean));
	Plot_S_bar_shift->Fill((S_bar_init-S_bar_mean));

        pullPlot_C->Fill(C_pull);
        pullPlot_D->Fill(D_pull);
        pullPlot_S->Fill(S_pull);
        pullPlot_D_bar->Fill(D_bar_pull);
        pullPlot_S_bar->Fill(S_bar_pull);

   }

   Plot_C_shift->SaveAs(("pull_result_new/C_shift_"+ parName +".root").c_str());
   Plot_D_shift->SaveAs(("pull_result_new/D_shift_"+ parName +".root").c_str());
   Plot_S_shift->SaveAs(("pull_result_new/S_shift_"+ parName +".root").c_str());
   Plot_D_bar_shift->SaveAs(("pull_result_new/D_bar_shift_"+ parName +".root").c_str());
   Plot_S_bar_shift->SaveAs(("pull_result_new/S_bar_shift_"+ parName +".root").c_str());

   pullPlot_C->SaveAs(("pull_result_new/pull_C_"+ parName +".root").c_str());
   pullPlot_D->SaveAs(("pull_result_new/pull_D_"+ parName +".root").c_str());
   pullPlot_S->SaveAs(("pull_result_new/pull_S_"+ parName +".root").c_str());
   pullPlot_D_bar->SaveAs(("pull_result_new/pull_D_bar_"+ parName +".root").c_str());
   pullPlot_S_bar->SaveAs(("pull_result_new/pull_S_bar_"+ parName +".root").c_str());

  //perform single-Gauss fit to pull distributions
   TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
   gaussian->SetParLimits(0, 1., 1000.);
   gaussian->SetParLimits(1,-1., 1.);
   gaussian->SetParLimits(2, 0., 2.);

   ofstream SummaryFile;
   SummaryFile.open(("pull_result_new/PullFile_" + parName + ".tex").c_str(),std::ofstream::trunc);
   SummaryFile << "\\begin{table}[hp!]" << "\n";
   SummaryFile << "\\centering" << "\n";
   SummaryFile << "\\caption{Pull parameters for CP coefficients from the toy studies for the time-dependent fit.}" << "\n";
   SummaryFile << "\\begin{tabular}{l | c | c}" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\hline" << "\n";


   double shift_C_mean = 0;
   double shift_D_mean = 0;
   double shift_S_mean = 0;
   double shift_D_bar_mean = 0;
   double shift_S_bar_mean = 0;

   double shift_C_sigma = 0;
   double shift_D_sigma = 0;
   double shift_S_sigma = 0;
   double shift_D_bar_sigma = 0;
   double shift_S_bar_sigma = 0;


   TCanvas* c = new TCanvas();
   gStyle->SetOptStat(11);
   gStyle->SetOptFit(111);

   Plot_C_shift->Fit(gaussian);
   cout << "C shift mean  :   " << gaussian->GetParameter(1) << " +/- " << gaussian->GetParError(1) <<"  C shift sigma:   " << gaussian->GetParameter(2) << " +/- " <<  gaussian->GetParError(2) << endl;
   shift_C_mean =  gaussian->GetParameter(1);
   shift_C_sigma =  gaussian->GetParameter(2);
   Plot_C_shift->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/C_shift_"+ parName + ".eps").c_str());

   Plot_D_shift->Fit(gaussian);
   cout << "D shift mean  :   " << gaussian->GetParameter(1) << " +/- " << gaussian->GetParError(1) <<"  D shift sigma:   " << gaussian->GetParameter(2) << " +/- " <<  gaussian->GetParError(2) << endl; 
   shift_D_mean =  gaussian->GetParameter(1);
   shift_D_sigma =  gaussian->GetParameter(2);
   Plot_D_shift->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/D_shift_"+ parName + ".eps").c_str());

   Plot_S_shift->Fit(gaussian);
   cout << "S shift mean  :   " << gaussian->GetParameter(1) << " +/- " << gaussian->GetParError(1) <<"  S shift sigma:   " << gaussian->GetParameter(2) << " +/- " <<  gaussian->GetParError(2) << endl; 
   shift_S_mean =  gaussian->GetParameter(1);
   shift_S_sigma =  gaussian->GetParameter(2);
   Plot_S_shift->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/S_shift_"+ parName + ".eps").c_str());

   Plot_D_bar_shift->Fit(gaussian);
   cout << "D_bar shift mean  :   " << gaussian->GetParameter(1) << " +/- " << gaussian->GetParError(1) <<"  D_bar shift sigma:   " << gaussian->GetParameter(2) << " +/- " <<  gaussian->GetParError(2) << endl; 
   shift_D_bar_mean =  gaussian->GetParameter(1);
   shift_D_bar_sigma =  gaussian->GetParameter(2);
   Plot_D_bar_shift->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/D_bar_shift_"+ parName + ".eps").c_str());

   Plot_S_bar_shift->Fit(gaussian);
   cout << "S_bar shift mean  :   " << gaussian->GetParameter(1) << " +/- " << gaussian->GetParError(1) <<"  S_bar shift sigma:   " << gaussian->GetParameter(2) << " +/- " <<  gaussian->GetParError(2) << endl;
   shift_S_bar_mean =  gaussian->GetParameter(1);
   shift_S_bar_sigma =  gaussian->GetParameter(2);
   Plot_S_bar_shift->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/S_bar_shift_"+ parName + ".eps").c_str());




   pullPlot_C->Fit(gaussian);
   SummaryFile << "C & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/pull_C_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D->Fit(gaussian);
   SummaryFile << "D & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/pull_D_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S->Fit(gaussian);
   SummaryFile << "S & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/pull_S_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar->Fit(gaussian);
   SummaryFile << "$\\bar{D}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/pull_D_bar_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar->Fit(gaussian);
   SummaryFile << "$\\bar{S}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result_new/pull_S_bar_"+ parName + "_Gaussfit.eps").c_str());

   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\end{tabular}" << "\n";
   SummaryFile << "\\label{table:Pulls_tDFit}" << "\n";
   SummaryFile << "\\end{table}";



  cout << "syst. error on C:  " << TMath::Sqrt( (shift_C_mean*shift_C_mean) + (shift_C_sigma*shift_C_sigma) ) << endl;
  cout << "syst. error on D:  " << TMath::Sqrt( (shift_D_mean*shift_D_mean) + (shift_D_sigma*shift_D_sigma) ) << endl;
  cout << "syst. error on S:  " << TMath::Sqrt( (shift_S_mean*shift_S_mean) + (shift_S_sigma*shift_S_sigma) ) << endl;
  cout << "syst. error on D_bar:  " << TMath::Sqrt( (shift_D_bar_mean*shift_D_bar_mean) + (shift_D_bar_sigma*shift_D_bar_sigma) ) << endl;
  cout << "syst. error on S_bar:  " << TMath::Sqrt( (shift_S_bar_mean*shift_S_bar_mean) + (shift_S_bar_sigma*shift_S_bar_sigma) ) << endl;

}

void pull::LoopSyst(string parName)
{
if (fChain == 0) return;

   TH1D* pullPlot_C0 = new TH1D("pullPlot_C","Pull C, varying c_{0}", 10, -0.2,0.2);
   TH1D* pullPlot_D0 = new TH1D("pullPlot_D","Pull D, varying c_{0}", 10, -1.,1.);
   TH1D* pullPlot_S0 = new TH1D("pullPlot_S","Pull S, varying c_{0}", 10, -0.2,0.2);
   TH1D* pullPlot_D_bar0 = new TH1D("pullPlot_D_bar","Pull D_bar, varying c_{0}", 10, -1.,1.);
   TH1D* pullPlot_S_bar0 = new TH1D("pullPlot_S_bar","Pull S_bar, varying c_{0}", 10, -0.2,0.2);

   TH1D* pullPlot_C1 = new TH1D(" pullPlot_C","Pull C, varying c_{1}", 10, -0.2,0.2);
   TH1D* pullPlot_D1 = new TH1D(" pullPlot_D","Pull D, varying c_{1}", 10, -1.,1.);
   TH1D* pullPlot_S1 = new TH1D(" pullPlot_S","Pull S, varying c_{1}", 10, -0.2,0.2);
   TH1D* pullPlot_D_bar1 = new TH1D(" pullPlot_D_bar","Pull D_bar, varying c_{1}", 11, -1.,1.);
   TH1D* pullPlot_S_bar1 = new TH1D(" pullPlot_S_bar","Pull S_bar, varying c_{1}", 11, -0.2,0.2);

   TH1D* pullPlot_C2 = new TH1D("pullPlot_C ","Pull C, varying c_{2}", 10, -0.2,0.2);
   TH1D* pullPlot_D2 = new TH1D("pullPlot_D ","Pull D, varying c_{2}", 10, -1.,1.);
   TH1D* pullPlot_S2 = new TH1D("pullPlot_S ","Pull S, varying c_{2}", 10, -0.2,0.2);
   TH1D* pullPlot_D_bar2 = new TH1D("pullPlot_D_bar ","Pull D_bar, varying c_{2}", 10, -1.,1.);
   TH1D* pullPlot_S_bar2 = new TH1D("pullPlot_S_bar ","Pull S_bar, varying c_{2}", 10, -0.2,0.2);

   TH1D* pullPlot_C3 = new TH1D(" pullPlot_C ","Pull C, varying c_{3}", 10, -0.2,0.2);
   TH1D* pullPlot_D3 = new TH1D(" pullPlot_D ","Pull D, varying c_{3}", 10, -1.,1.);
   TH1D* pullPlot_S3 = new TH1D(" pullPlot_S ","Pull S, varying c_{3}", 10, -0.2,0.2);
   TH1D* pullPlot_D_bar3 = new TH1D(" pullPlot_D_bar ","Pull D_bar, varying c_{3}", 10, -1.,1.);
   TH1D* pullPlot_S_bar3 = new TH1D(" pullPlot_S_bar ","Pull S_bar, varying c_{3}", 10, -0.2,0.2);

/*
   TH1D* Plot_C_shift = new TH1D("Plot_C_shift","C_{init}-C_{fit}", 10, -0.5,0.5);
   TH1D* Plot_D_shift = new TH1D("Plot_D_shift","D_{init}-D_{fit}", 10, -0.5,0.5);
   TH1D* Plot_S_shift = new TH1D("Plot_S_shift","S_{init}-S_{fit}", 10, -0.5,0.5);
   TH1D* Plot_D_bar_shift = new TH1D("Plot_D_bar_shift","D_bar_{init}-D_bar_{fit}", 10, -0.5,0.5);
   TH1D* Plot_S_bar_shift = new TH1D("Plot_S_bar_shift","S_bar_{init}-S_bar_{fit}", 10, -0.5,0.5);
*/
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

        pullPlot_C0->Fill(delta_pull_C0);
        pullPlot_D0->Fill(delta_pull_D0);
        pullPlot_S0->Fill(delta_pull_S0);
        pullPlot_D_bar0->Fill(delta_pull_D_bar0);
        pullPlot_S_bar0->Fill(delta_pull_S_bar0);

        pullPlot_C1->Fill(delta_pull_C1);
        pullPlot_D1->Fill(delta_pull_D1);
        pullPlot_S1->Fill(delta_pull_S1);
        pullPlot_D_bar1->Fill(delta_pull_D_bar1);
        pullPlot_S_bar1->Fill(delta_pull_S_bar1);

        pullPlot_C2->Fill(delta_pull_C2);
        pullPlot_D2->Fill(delta_pull_D2);
        pullPlot_S2->Fill(delta_pull_S2);
        pullPlot_D_bar2->Fill(delta_pull_D_bar2);
        pullPlot_S_bar2->Fill(delta_pull_S_bar2);

        pullPlot_C3->Fill(delta_pull_C3);
        pullPlot_D3->Fill(delta_pull_D3);
        pullPlot_S3->Fill(delta_pull_S3);
        pullPlot_D_bar3->Fill(delta_pull_D_bar3);
        pullPlot_S_bar3->Fill(delta_pull_S_bar3);

   }

   pullPlot_C0->SaveAs(("pull_result/pull_C_0_"+ parName +".root").c_str());
   pullPlot_D0->SaveAs(("pull_result/pull_D_0_"+ parName +".root").c_str());
   pullPlot_S0->SaveAs(("pull_result/pull_S_0_"+ parName +".root").c_str());
   pullPlot_D_bar0->SaveAs(("pull_result/pull_D_bar_0_"+ parName +".root").c_str());
   pullPlot_S_bar0->SaveAs(("pull_result/pull_S_bar_0_"+ parName +".root").c_str());

   pullPlot_C1->SaveAs(("pull_result/pull_C_1_"+ parName +".root").c_str());
   pullPlot_D1->SaveAs(("pull_result/pull_D_1_"+ parName +".root").c_str());
   pullPlot_S1->SaveAs(("pull_result/pull_S_1_"+ parName +".root").c_str());
   pullPlot_D_bar1->SaveAs(("pull_result/pull_D_bar_1_"+ parName +".root").c_str());
   pullPlot_S_bar1->SaveAs(("pull_result/pull_S_bar_1_"+ parName +".root").c_str());

   pullPlot_C2->SaveAs(("pull_result/pull_C_2_"+ parName +".root").c_str());
   pullPlot_D2->SaveAs(("pull_result/pull_D_2_"+ parName +".root").c_str());
   pullPlot_S2->SaveAs(("pull_result/pull_S_2_"+ parName +".root").c_str());
   pullPlot_D_bar2->SaveAs(("pull_result/pull_D_bar_2_"+ parName +".root").c_str());
   pullPlot_S_bar2->SaveAs(("pull_result/pull_S_bar_2_"+ parName +".root").c_str());

   pullPlot_C3->SaveAs(("pull_result/pull_C_3_"+ parName +".root").c_str());
   pullPlot_D3->SaveAs(("pull_result/pull_D_3_"+ parName +".root").c_str());
   pullPlot_S3->SaveAs(("pull_result/pull_S_3_"+ parName +".root").c_str());
   pullPlot_D_bar3->SaveAs(("pull_result/pull_D_bar_3_"+ parName +".root").c_str());
   pullPlot_S_bar3->SaveAs(("pull_result/pull_S_bar_3_"+ parName +".root").c_str());

  //perform single-Gauss fit to pull distributions
   TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
   gaussian->SetParLimits(0, 1., 1000.);
   gaussian->SetParLimits(1,-1., 1.);
   gaussian->SetParLimits(2, 0., 2.);

   ofstream SummaryFile;
   SummaryFile.open(("pull_result/PullFile_" + parName + ".tex").c_str(),std::ofstream::trunc);
   SummaryFile << "\\begin{table}[hp!]" << "\n";
   SummaryFile << "\\centering" << "\n";
   SummaryFile << "\\caption{Pull parameters for CP coefficients from the toy studies for the time-dependent fit.}" << "\n";
   SummaryFile << "\\begin{tabular}{l | c | c}" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\hline" << "\n";


   gStyle->SetOptStat(11);
   gStyle->SetOptFit(111);
   TCanvas* c = new TCanvas();

   double C0_Gaussmean = 0;
   double C1_Gaussmean = 0;
   double C2_Gaussmean = 0;
   double C3_Gaussmean = 0;

   double D0_Gaussmean = 0;
   double D1_Gaussmean = 0;
   double D2_Gaussmean = 0;
   double D3_Gaussmean = 0;

   double S0_Gaussmean = 0;
   double S1_Gaussmean = 0;
   double S2_Gaussmean = 0;
   double S3_Gaussmean = 0;

   double D_bar0_Gaussmean = 0;
   double D_bar1_Gaussmean = 0;
   double D_bar2_Gaussmean = 0;
   double D_bar3_Gaussmean = 0;

   double S_bar0_Gaussmean = 0;
   double S_bar1_Gaussmean = 0;
   double S_bar2_Gaussmean = 0;
   double S_bar3_Gaussmean = 0;


   double C0_Gausserror = 0;
   double C1_Gausserror = 0;
   double C2_Gausserror = 0;
   double C3_Gausserror = 0;

   double D0_Gausserror = 0;
   double D1_Gausserror = 0;
   double D2_Gausserror = 0;
   double D3_Gausserror = 0;

   double S0_Gausserror = 0;
   double S1_Gausserror = 0;
   double S2_Gausserror = 0;
   double S3_Gausserror = 0;

   double D_bar0_Gausserror = 0;
   double D_bar1_Gausserror = 0;
   double D_bar2_Gausserror = 0;
   double D_bar3_Gausserror = 0;

   double S_bar0_Gausserror = 0;
   double S_bar1_Gausserror = 0;
   double S_bar2_Gausserror = 0;
   double S_bar3_Gausserror = 0;

   pullPlot_C0->Fit(gaussian);
   C0_Gaussmean = gaussian->GetParameter(1);
   C0_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "C, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C1->Fit(gaussian);
   C1_Gaussmean = gaussian->GetParameter(1);
   C1_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "C, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C2->Fit(gaussian);
   C2_Gaussmean = gaussian->GetParameter(1);
   C2_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "C, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C3->Fit(gaussian);
   C3_Gaussmean = gaussian->GetParameter(1);
   C3_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "C, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_D0->Fit(gaussian);
   D0_Gaussmean = gaussian->GetParameter(1);
   D0_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "D, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D1->Fit(gaussian);
   D1_Gaussmean = gaussian->GetParameter(1);
   D1_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "D, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D2->Fit(gaussian);
   D2_Gaussmean = gaussian->GetParameter(1);
   D2_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "D, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D3->Fit(gaussian);
   D3_Gaussmean = gaussian->GetParameter(1);
   D3_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "D, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_S0->Fit(gaussian);
   S0_Gaussmean = gaussian->GetParameter(1);
   S0_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "S, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S1->Fit(gaussian);
   S1_Gaussmean = gaussian->GetParameter(1);
   S1_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "S, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S2->Fit(gaussian);
   S2_Gaussmean = gaussian->GetParameter(1);
   S2_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "S, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S3->Fit(gaussian);
   S3_Gaussmean = gaussian->GetParameter(1);
   S3_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "S, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_D_bar0->Fit(gaussian);
   D_bar0_Gaussmean = gaussian->GetParameter(1);
   D_bar0_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{D}$, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar1->Fit(gaussian);
   D_bar1_Gaussmean = gaussian->GetParameter(1);
   D_bar1_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{D}$, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar2->Fit(gaussian);
   D_bar2_Gaussmean = gaussian->GetParameter(1);
   D_bar2_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{D}$, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar3->Fit(gaussian);
   D_bar3_Gaussmean = gaussian->GetParameter(1);
   D_bar3_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{D}$, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_S_bar0->Fit(gaussian);
   S_bar0_Gaussmean = gaussian->GetParameter(1);
   S_bar0_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{S}$, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar1->Fit(gaussian);
   S_bar1_Gaussmean = gaussian->GetParameter(1);
   S_bar1_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{S}$, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar2->Fit(gaussian);
   S_bar2_Gaussmean = gaussian->GetParameter(1);
   S_bar2_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{S}$, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar3->Fit(gaussian);
   S_bar3_Gaussmean = gaussian->GetParameter(1);
   S_bar3_Gausserror = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{S}$, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_3_"+ parName + "_Gaussfit.eps").c_str());



   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\end{tabular}" << "\n";
   SummaryFile << "\\label{table:Pulls_tDFit}" << "\n";
   SummaryFile << "\\end{table}";

   cout << "combined shift in C:   "  <<  TMath::Sqrt((C0_Gaussmean * C0_Gaussmean) + (C0_Gausserror * C0_Gausserror) + (C1_Gausserror * C1_Gausserror) + (C1_Gaussmean * C1_Gaussmean) + (C2_Gaussmean * C2_Gaussmean) + (C2_Gausserror * C2_Gausserror) + (C3_Gaussmean * C3_Gaussmean) + (C3_Gausserror * C3_Gausserror)) << endl;

   cout << "combined shift in D:   "  <<  TMath::Sqrt((D0_Gaussmean * D0_Gaussmean) + (D0_Gausserror * D0_Gausserror) + (D1_Gausserror * D1_Gausserror) + (D1_Gaussmean * D1_Gaussmean) + (D2_Gaussmean * D2_Gaussmean) + (D2_Gausserror * D2_Gausserror) + (D3_Gaussmean * D3_Gaussmean) + (D3_Gausserror * D3_Gausserror)) << endl;

   cout << "combined shift in S:   "  <<  TMath::Sqrt((S0_Gaussmean * S0_Gaussmean) + (S0_Gausserror * S0_Gausserror) + (S1_Gausserror * S1_Gausserror) + (S1_Gaussmean * S1_Gaussmean) + (S2_Gaussmean * S2_Gaussmean) + (S2_Gausserror * S2_Gausserror) + (S3_Gaussmean * S3_Gaussmean) + (S3_Gausserror * S3_Gausserror)) << endl;

   cout << "combined shift in D:   "  <<  TMath::Sqrt((D_bar0_Gaussmean * D_bar0_Gaussmean) + (D_bar0_Gausserror * D_bar0_Gausserror) + (D_bar1_Gausserror * D_bar1_Gausserror) + (D_bar1_Gaussmean * D_bar1_Gaussmean) + (D_bar2_Gaussmean * D_bar2_Gaussmean) + (D_bar2_Gausserror * D_bar2_Gausserror) + (D_bar3_Gaussmean * D_bar3_Gaussmean) + (D_bar3_Gausserror * D_bar3_Gausserror)) << endl;

   cout << "combined shift in D:   "  <<  TMath::Sqrt((S_bar0_Gaussmean * S_bar0_Gaussmean) + (S_bar0_Gausserror * S_bar0_Gausserror) + (S_bar1_Gausserror * S_bar1_Gausserror) + (S_bar1_Gaussmean * S_bar1_Gaussmean) + (S_bar2_Gaussmean * S_bar2_Gaussmean) + (S_bar2_Gausserror * S_bar2_Gausserror) + (S_bar3_Gaussmean * S_bar3_Gaussmean) + (S_bar3_Gausserror * S_bar3_Gausserror)) << endl;

/*
   cout << "combined shift in D:   "  <<  TMath::Sqrt((D0_Gaussmean * D0_Gaussmean) + (D1_Gaussmean * D1_Gaussmean) + (D2_Gaussmean * D2_Gaussmean) + (D3_Gaussmean * D3_Gaussmean)) << endl;
   cout << "combined shift in S:   "  <<  TMath::Sqrt((S0_Gaussmean * S0_Gaussmean) + (S1_Gaussmean * S1_Gaussmean) + (S2_Gaussmean * S2_Gaussmean) + (S3_Gaussmean * S3_Gaussmean)) << endl;
   cout << "combined shift in D_bar:   "  <<  TMath::Sqrt((D_bar0_Gaussmean * D_bar0_Gaussmean) + (D_bar1_Gaussmean * D_bar1_Gaussmean) + (D_bar2_Gaussmean * D_bar2_Gaussmean) + (D_bar3_Gaussmean * D_bar3_Gaussmean)) << endl;
   cout << "combined shift in S_bar:   "  <<  TMath::Sqrt((S_bar0_Gaussmean * S_bar0_Gaussmean) + (S_bar1_Gaussmean * S_bar1_Gaussmean) + (S_bar2_Gaussmean * S_bar2_Gaussmean) + (S_bar3_Gaussmean * S_bar3_Gaussmean)) << endl;
*/


}

void pull::LoopSyst_noChol(string parName)
{
   if (fChain == 0) return;

   TH1D* pullPlot_C = new TH1D("pullPlot_C","Pull C", 10, -0.2,0.2);
   TH1D* pullPlot_D = new TH1D("pullPlot_D","Pull D", 10, -1.,1.);
   TH1D* pullPlot_S = new TH1D("pullPlot_S","Pull S", 10, -0.2,0.2);
   TH1D* pullPlot_D_bar = new TH1D("pullPlot_D_bar","Pull D_bar", 10, -1.,1.);
   TH1D* pullPlot_S_bar = new TH1D("pullPlot_S_bar","Pull S_bar", 10, -0.2,0.2);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

        pullPlot_C->Fill(delta_pull_C);
        pullPlot_D->Fill(delta_pull_D);
        pullPlot_S->Fill(delta_pull_S);
        pullPlot_D_bar->Fill(delta_pull_D_bar);
        pullPlot_S_bar->Fill(delta_pull_S_bar);

   }

   pullPlot_C->SaveAs(("pull_result/pull_C_"+ parName +".root").c_str());
   pullPlot_D->SaveAs(("pull_result/pull_D_"+ parName +".root").c_str());
   pullPlot_S->SaveAs(("pull_result/pull_S_"+ parName +".root").c_str());
   pullPlot_D_bar->SaveAs(("pull_result/pull_D_bar_"+ parName +".root").c_str());
   pullPlot_S_bar->SaveAs(("pull_result/pull_S_bar_"+ parName +".root").c_str());

  //perform single-Gauss fit to pull distributions
   TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
   gaussian->SetParLimits(0, 1., 1000.);
   gaussian->SetParLimits(1,-1., 1.);
   gaussian->SetParLimits(2, 0., 2.);

   ofstream SummaryFile;
   SummaryFile.open(("pull_result/PullFile_" + parName + ".tex").c_str(),std::ofstream::trunc);
   SummaryFile << "\\begin{table}[hp!]" << "\n";
   SummaryFile << "\\centering" << "\n";
   SummaryFile << "\\caption{Pull parameters for CP coefficients from the toy studies for the time-dependent fit.}" << "\n";
   SummaryFile << "\\begin{tabular}{l | c | c}" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\hline" << "\n";


   gStyle->SetOptStat(11);
   gStyle->SetOptFit(111);
   TCanvas* c = new TCanvas();

   double C_Gaussmean = 0;
   double D_Gaussmean = 0;
   double S_Gaussmean = 0;
   double D_bar_Gaussmean = 0;
   double S_bar_Gaussmean = 0;

   double C_Gausswidth = 0;
   double D_Gausswidth = 0;
   double S_Gausswidth = 0;
   double D_bar_Gausswidth = 0;
   double S_bar_Gausswidth = 0;

   pullPlot_C->Fit(gaussian);
   C_Gaussmean = gaussian->GetParameter(1);
   C_Gausswidth = gaussian->GetParameter(2);
   SummaryFile << "C & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D->Fit(gaussian);
   D_Gaussmean = gaussian->GetParameter(1);
   D_Gausswidth = gaussian->GetParameter(2);
   SummaryFile << "D & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S->Fit(gaussian);
   S_Gaussmean = gaussian->GetParameter(1);
   S_Gausswidth = gaussian->GetParameter(2);
   SummaryFile << "S & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar->Fit(gaussian);
   D_bar_Gaussmean = gaussian->GetParameter(1);
   D_bar_Gausswidth = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{D}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar->Fit(gaussian);
   S_bar_Gaussmean = gaussian->GetParameter(1);
   S_bar_Gausswidth = gaussian->GetParameter(2);
   SummaryFile << "$\\bar{S}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_"+ parName + "_Gaussfit.eps").c_str());

   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\end{tabular}" << "\n";
   SummaryFile << "\\label{table:Pulls_tDFit}" << "\n";
   SummaryFile << "\\end{table}";


   cout << "combined syst. error in C:   "  << TMath::Sqrt( (C_Gaussmean*C_Gaussmean) + (C_Gausswidth*C_Gausswidth) ) << endl;
   cout << "combined syst. error in D:   "  << TMath::Sqrt( (D_Gaussmean*D_Gaussmean) + (D_Gausswidth*D_Gausswidth) ) << endl;
   cout << "combined syst. error in S:   "  << TMath::Sqrt( (S_Gaussmean*S_Gaussmean) + (S_Gausswidth*S_Gausswidth) )  << endl;
   cout << "combined syst. error in D_bar:   "  << TMath::Sqrt( (D_bar_Gaussmean*D_bar_Gaussmean) + (D_bar_Gausswidth*D_bar_Gausswidth) ) << endl;
   cout << "combined syst. error in S_bar:   "  << TMath::Sqrt( (S_bar_Gaussmean*S_bar_Gaussmean) + (S_bar_Gausswidth*S_bar_Gausswidth) )  << endl;
}

void pull::getShift(string parName)
{
   if (fChain == 0) return;


   Long64_t nentries = fChain->GetEntriesFast();

   double shift_C = 0;
   double shift_D = 0;
   double shift_S = 0;
   double shift_D_bar = 0;
   double shift_S_bar = 0;

   int nen = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;


     shift_C += shift_C + (C_init-C_mean);
     shift_D += shift_D + (D_init-D_mean);
     shift_S += shift_S + (S_init-S_mean);
     shift_D_bar += shift_D_bar + (D_bar_init-D_bar_mean);
     shift_S_bar += shift_S_bar + (S_bar_init-S_bar_mean);

     nen++;
   }

   cout << "mean shift for C_" << parName.c_str() << ":   " << shift_C/nen << endl;
   cout << "mean shift for D_" << parName.c_str() << ":   " << shift_D/nen << endl;
   cout << "mean shift for S_" << parName.c_str() << ":   " << shift_S/nen << endl;
   cout << "mean shift for D_bar_" << parName.c_str() << ":   " << shift_D_bar/nen << endl;
   cout << "mean shift for S_bar_" << parName.c_str() << ":   " << shift_S_bar/nen << endl;


   ofstream SummaryFile;
   SummaryFile.open(("pull_result_new/ShiftFile_" + parName + ".txt").c_str(),std::ofstream::trunc);
   SummaryFile << "C_shift = " << shift_C/nen  << "\n";
   SummaryFile << "D_shift = " << shift_D/nen  << "\n";
   SummaryFile << "S_shift = " << shift_S/nen  << "\n";
   SummaryFile << "D_bar_shift = " << shift_D_bar/nen  << "\n";
   SummaryFile << "S_bar_shift = " << shift_S_bar/nen;
}


void pull::makeDeltaPulls_Col()
{

 TChain* chain_noSyst =  new TChain("MinuitParameterSetNtp");

 for(int i=1; i < 2001; i++){
	stringstream number;
	number << i;
 	chain_noSyst->Add(("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_" + number.str() + ".root").c_str());
 }

 TChain* chain_c0 =  new TChain("MinuitParameterSetNtp");
 TChain* chain_c1 =  new TChain("MinuitParameterSetNtp");
 TChain* chain_c2 =  new TChain("MinuitParameterSetNtp");
 TChain* chain_c3 =  new TChain("MinuitParameterSetNtp");

 for(int i=1; i < 501; i++){
	stringstream number;
	number << i;
 	chain_c0->Add(("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_par0_" + number.str() + ".root").c_str());
 	chain_c1->Add(("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_par1_" + number.str() + ".root").c_str());
	chain_c2->Add(("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_par2_" + number.str() + ".root").c_str());
 	chain_c3->Add(("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_par3_" + number.str() + ".root").c_str());
}

 Double_t C_pull_c0;
 Double_t C_pull_c1;
 Double_t C_pull_c2;
 Double_t C_pull_c3;

 Double_t D_pull_c0;
 Double_t D_pull_c1;
 Double_t D_pull_c2;
 Double_t D_pull_c3;

 Double_t S_pull_c0;
 Double_t S_pull_c1;
 Double_t S_pull_c2;
 Double_t S_pull_c3;

 Double_t D_bar_pull_c0;
 Double_t D_bar_pull_c1;
 Double_t D_bar_pull_c2;
 Double_t D_bar_pull_c3;

 Double_t S_bar_pull_c0;
 Double_t S_bar_pull_c1;
 Double_t S_bar_pull_c2;
 Double_t S_bar_pull_c3;

 Double_t C_pull_noSyst;
 Double_t D_pull_noSyst;
 Double_t S_pull_noSyst;
 Double_t D_bar_pull_noSyst;
 Double_t S_bar_pull_noSyst;


 Double_t C_mean_noSyst;
 Double_t D_mean_noSyst;
 Double_t S_mean_noSyst;
 Double_t D_bar_mean_noSyst;
 Double_t S_bar_mean_noSyst;

 Double_t C_err_noSyst;
 Double_t D_err_noSyst;
 Double_t S_err_noSyst;
 Double_t D_bar_err_noSyst;
 Double_t S_bar_err_noSyst;

 Double_t C_mean_c0;
 Double_t C_mean_c1;
 Double_t C_mean_c2;
 Double_t C_mean_c3;

 Double_t D_mean_c0;
 Double_t D_mean_c1;
 Double_t D_mean_c2;
 Double_t D_mean_c3;

 Double_t S_mean_c0;
 Double_t S_mean_c1;
 Double_t S_mean_c2;
 Double_t S_mean_c3;

 Double_t D_bar_mean_c0;
 Double_t D_bar_mean_c1;
 Double_t D_bar_mean_c2;
 Double_t D_bar_mean_c3;

 Double_t S_bar_mean_c0;
 Double_t S_bar_mean_c1;
 Double_t S_bar_mean_c2;
 Double_t S_bar_mean_c3;


 chain_c0 -> SetBranchAddress("C_pull" , &C_pull_c0 );
 chain_c1 -> SetBranchAddress("C_pull" , &C_pull_c1 );
 chain_c2 -> SetBranchAddress("C_pull" , &C_pull_c2 );
 chain_c3 -> SetBranchAddress("C_pull" , &C_pull_c3 );

 chain_c0 -> SetBranchAddress("D_pull" , &D_pull_c0 );
 chain_c1 -> SetBranchAddress("D_pull" , &D_pull_c1 );
 chain_c2 -> SetBranchAddress("D_pull" , &D_pull_c2 );
 chain_c3 -> SetBranchAddress("D_pull" , &D_pull_c3 );

 chain_c0 -> SetBranchAddress("S_pull" , &S_pull_c0 );
 chain_c1 -> SetBranchAddress("S_pull" , &S_pull_c1 );
 chain_c2 -> SetBranchAddress("S_pull" , &S_pull_c2 );
 chain_c3 -> SetBranchAddress("S_pull" , &S_pull_c3 );

 chain_c0 -> SetBranchAddress("D_bar_pull" , &D_bar_pull_c0 );
 chain_c1 -> SetBranchAddress("D_bar_pull" , &D_bar_pull_c1 );
 chain_c2 -> SetBranchAddress("D_bar_pull" , &D_bar_pull_c2 );
 chain_c3 -> SetBranchAddress("D_bar_pull" , &D_bar_pull_c3 );

 chain_c0 -> SetBranchAddress("S_bar_pull" , &S_bar_pull_c0 );
 chain_c1 -> SetBranchAddress("S_bar_pull" , &S_bar_pull_c1 );
 chain_c2 -> SetBranchAddress("S_bar_pull" , &S_bar_pull_c2 );
 chain_c3 -> SetBranchAddress("S_bar_pull" , &S_bar_pull_c3 );


 chain_noSyst -> SetBranchAddress("C_pull" , &C_pull_noSyst );
 chain_noSyst -> SetBranchAddress("D_pull" , &D_pull_noSyst );
 chain_noSyst -> SetBranchAddress("S_pull" , &S_pull_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_pull" , &D_bar_pull_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_pull" , &S_bar_pull_noSyst );



 chain_noSyst -> SetBranchAddress("C_mean" , &C_mean_noSyst );
 chain_noSyst -> SetBranchAddress("D_mean" , &D_mean_noSyst );
 chain_noSyst -> SetBranchAddress("S_mean" , &S_mean_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_mean" , &D_bar_mean_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_mean" , &S_bar_mean_noSyst );

 chain_noSyst -> SetBranchAddress("C_err" , &C_err_noSyst );
 chain_noSyst -> SetBranchAddress("D_err" , &D_err_noSyst );
 chain_noSyst -> SetBranchAddress("S_err" , &S_err_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_err" , &D_bar_err_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_err" , &S_bar_err_noSyst );

 chain_c0 -> SetBranchAddress("C_mean" , &C_mean_c0 );
 chain_c1 -> SetBranchAddress("C_mean" , &C_mean_c1 );
 chain_c2 -> SetBranchAddress("C_mean" , &C_mean_c2 );
 chain_c3 -> SetBranchAddress("C_mean" , &C_mean_c3 );

 chain_c0 -> SetBranchAddress("D_mean" , &D_mean_c0 );
 chain_c1 -> SetBranchAddress("D_mean" , &D_mean_c1 );
 chain_c2 -> SetBranchAddress("D_mean" , &D_mean_c2 );
 chain_c3 -> SetBranchAddress("D_mean" , &D_mean_c3 );

 chain_c0 -> SetBranchAddress("S_mean" , &S_mean_c0 );
 chain_c1 -> SetBranchAddress("S_mean" , &S_mean_c1 );
 chain_c2 -> SetBranchAddress("S_mean" , &S_mean_c2 );
 chain_c3 -> SetBranchAddress("S_mean" , &S_mean_c3 );

 chain_c0 -> SetBranchAddress("D_bar_mean" , &D_bar_mean_c0 );
 chain_c1 -> SetBranchAddress("D_bar_mean" , &D_bar_mean_c1 );
 chain_c2 -> SetBranchAddress("D_bar_mean" , &D_bar_mean_c2 );
 chain_c3 -> SetBranchAddress("D_bar_mean" , &D_bar_mean_c3 );

 chain_c0 -> SetBranchAddress("S_bar_mean" , &S_bar_mean_c0 );
 chain_c1 -> SetBranchAddress("S_bar_mean" , &S_bar_mean_c1 );
 chain_c2 -> SetBranchAddress("S_bar_mean" , &S_bar_mean_c2 );
 chain_c3 -> SetBranchAddress("S_bar_mean" , &S_bar_mean_c3 );

 int entries = chain_c0->GetEntries();

 Double_t delta_pull_C0;
 Double_t delta_pull_D0;
 Double_t delta_pull_S0;
 Double_t delta_pull_D_bar0;
 Double_t delta_pull_S_bar0;

 Double_t delta_pull_C1;
 Double_t delta_pull_D1;
 Double_t delta_pull_S1;
 Double_t delta_pull_D_bar1;
 Double_t delta_pull_S_bar1;

 Double_t delta_pull_C2;
 Double_t delta_pull_D2;
 Double_t delta_pull_S2;
 Double_t delta_pull_D_bar2;
 Double_t delta_pull_S_bar2;

 Double_t delta_pull_C3;
 Double_t delta_pull_D3;
 Double_t delta_pull_S3;
 Double_t delta_pull_D_bar3;
 Double_t delta_pull_S_bar3;

 TFile* output = new TFile("tdfit_deltaPulls_Chol.root","RECREATE");
 TTree* new_tree = new TTree("MinuitParameterSetNtp","MinuitParameterSetNtp");

 new_tree->Branch("delta_pull_C0",&delta_pull_C0,"delta_pull_C0/D");
 new_tree->Branch("delta_pull_D0",&delta_pull_D0,"delta_pull_D0/D");
 new_tree->Branch("delta_pull_S0",&delta_pull_S0,"delta_pull_S0/D");
 new_tree->Branch("delta_pull_D_bar0",&delta_pull_D_bar0,"delta_pull_D_bar0/D");
 new_tree->Branch("delta_pull_S_bar0",&delta_pull_S_bar0,"delta_pull_S_bar0/D");

 new_tree->Branch("delta_pull_C1",&delta_pull_C1,"delta_pull_C1/D");
 new_tree->Branch("delta_pull_D1",&delta_pull_D1,"delta_pull_D1/D");
 new_tree->Branch("delta_pull_S1",&delta_pull_S1,"delta_pull_S1/D");
 new_tree->Branch("delta_pull_D_bar1",&delta_pull_D_bar1,"delta_pull_D_bar1/D");
 new_tree->Branch("delta_pull_S_bar1",&delta_pull_S_bar1,"delta_pull_S_bar1/D");

 new_tree->Branch("delta_pull_C2",&delta_pull_C2,"delta_pull_C2/D");
 new_tree->Branch("delta_pull_D2",&delta_pull_D2,"delta_pull_D2/D");
 new_tree->Branch("delta_pull_S2",&delta_pull_S2,"delta_pull_S2/D");
 new_tree->Branch("delta_pull_D_bar2",&delta_pull_D_bar2,"delta_pull_D_bar2/D");
 new_tree->Branch("delta_pull_S_bar2",&delta_pull_S_bar2,"delta_pull_S_bar2/D");

 new_tree->Branch("delta_pull_C3",&delta_pull_C3,"delta_pull_C3/D");
 new_tree->Branch("delta_pull_D3",&delta_pull_D3,"delta_pull_D3/D");
 new_tree->Branch("delta_pull_S3",&delta_pull_S3,"delta_pull_S3/D");
 new_tree->Branch("delta_pull_D_bar3",&delta_pull_D_bar3,"delta_pull_D_bar3/D");
 new_tree->Branch("delta_pull_S_bar3",&delta_pull_S_bar3,"delta_pull_S_bar3/D");

 for(int i=0; i < entries; i++){

 	chain_c0->GetEntry(i);
 	chain_noSyst->GetEntry(i);

 	delta_pull_C0 = (C_mean_noSyst - C_mean_c0) / C_err_noSyst;
 	delta_pull_D0 = (D_mean_noSyst - D_mean_c0) / D_err_noSyst;
 	delta_pull_S0 = (S_mean_noSyst - S_mean_c0) / S_err_noSyst;
 	delta_pull_D_bar0 = (D_bar_mean_noSyst - D_bar_mean_c0) / D_bar_err_noSyst;
 	delta_pull_S_bar0 = (S_bar_mean_noSyst - S_bar_mean_c0) / S_bar_err_noSyst;


        chain_c1->GetEntry(i);
        chain_noSyst->GetEntry(i+entries);

 	delta_pull_C1 = (C_mean_noSyst - C_mean_c1) / C_err_noSyst;
 	delta_pull_D1 = (D_mean_noSyst - D_mean_c1) / D_err_noSyst;
 	delta_pull_S1 = (S_mean_noSyst - S_mean_c1) / S_err_noSyst;
 	delta_pull_D_bar1 = (D_bar_mean_noSyst - D_bar_mean_c1) / D_bar_err_noSyst;
 	delta_pull_S_bar1 = (S_bar_mean_noSyst - S_bar_mean_c1) / S_bar_err_noSyst;


        chain_c2->GetEntry(i);
        chain_noSyst->GetEntry(i+(2*entries));


 	delta_pull_C2 = (C_mean_noSyst - C_mean_c2) / C_err_noSyst;
 	delta_pull_D2 = (D_mean_noSyst - D_mean_c2) / D_err_noSyst;
 	delta_pull_S2 = (S_mean_noSyst - S_mean_c2) / S_err_noSyst;
 	delta_pull_D_bar2 = (D_bar_mean_noSyst - D_bar_mean_c2) / D_bar_err_noSyst;
 	delta_pull_S_bar2 = (S_bar_mean_noSyst - S_bar_mean_c2) / S_bar_err_noSyst;


        chain_c3->GetEntry(i);
        chain_noSyst->GetEntry(i+(3*entries));


 	delta_pull_C3 = (C_mean_noSyst - C_mean_c3) / C_err_noSyst;
 	delta_pull_D3 = (D_mean_noSyst - D_mean_c3) / D_err_noSyst;
 	delta_pull_S3 = (S_mean_noSyst - S_mean_c3) / S_err_noSyst;
 	delta_pull_D_bar3 = (D_bar_mean_noSyst - D_bar_mean_c3) / D_bar_err_noSyst;
 	delta_pull_S_bar3 = (S_bar_mean_noSyst - S_bar_mean_c3) / S_bar_err_noSyst;

 	new_tree->Fill();

 }

  new_tree->Write();
  output->Close();

}


void pull::makeDeltaPulls_noCol()
{

 TChain* chain_noSyst =  new TChain("MinuitParameterSetNtp");
 chain_noSyst->Add("/work/kecke/Promotion/Bs2DsKpipi/Bs2DsKpipi_repository/TD-MINT2/src/Users/dargent/timeFit/signal_toy_fitterValid/pull_*.root"); //1



 TChain* chain_syst =  new TChain("MinuitParameterSetNtp");
 chain_syst->Add("/auto/data/kecke/BsDsKpipi/toys/signal_noSyst/pull_*.root");


 Double_t C_pull_syst;
 Double_t D_pull_syst;
 Double_t S_pull_syst;
 Double_t D_bar_pull_syst;
 Double_t S_bar_pull_syst;

 Double_t C_pull_noSyst;
 Double_t D_pull_noSyst;
 Double_t S_pull_noSyst;
 Double_t D_bar_pull_noSyst;
 Double_t S_bar_pull_noSyst;


 Double_t C_mean_syst;
 Double_t D_mean_syst;
 Double_t S_mean_syst;
 Double_t D_bar_mean_syst;
 Double_t S_bar_mean_syst;

 Double_t C_mean_noSyst;
 Double_t D_mean_noSyst;
 Double_t S_mean_noSyst;
 Double_t D_bar_mean_noSyst;
 Double_t S_bar_mean_noSyst;

 Double_t C_err_noSyst;
 Double_t D_err_noSyst;
 Double_t S_err_noSyst;
 Double_t D_bar_err_noSyst;
 Double_t S_bar_err_noSyst;


 chain_syst -> SetBranchAddress("C_pull" , &C_pull_syst );
 chain_syst -> SetBranchAddress("D_pull" , &D_pull_syst );
 chain_syst -> SetBranchAddress("S_pull" , &S_pull_syst );
 chain_syst -> SetBranchAddress("D_bar_pull" , &D_bar_pull_syst );
 chain_syst -> SetBranchAddress("S_bar_pull" , &S_bar_pull_syst );


 chain_noSyst -> SetBranchAddress("C_pull" , &C_pull_noSyst );
 chain_noSyst -> SetBranchAddress("D_pull" , &D_pull_noSyst );
 chain_noSyst -> SetBranchAddress("S_pull" , &S_pull_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_pull" , &D_bar_pull_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_pull" , &S_bar_pull_noSyst );


 chain_syst -> SetBranchAddress("C_mean" , &C_mean_syst );
 chain_syst -> SetBranchAddress("D_mean" , &D_mean_syst );
 chain_syst -> SetBranchAddress("S_mean" , &S_mean_syst );
 chain_syst -> SetBranchAddress("D_bar_mean" , &D_bar_mean_syst );
 chain_syst -> SetBranchAddress("S_bar_mean" , &S_bar_mean_syst );

 chain_noSyst -> SetBranchAddress("C_mean" , &C_mean_noSyst );
 chain_noSyst -> SetBranchAddress("D_mean" , &D_mean_noSyst );
 chain_noSyst -> SetBranchAddress("S_mean" , &S_mean_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_mean" , &D_bar_mean_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_mean" , &S_bar_mean_noSyst );

 chain_noSyst -> SetBranchAddress("C_err" , &C_err_noSyst );
 chain_noSyst -> SetBranchAddress("D_err" , &D_err_noSyst );
 chain_noSyst -> SetBranchAddress("S_err" , &S_err_noSyst );
 chain_noSyst -> SetBranchAddress("D_bar_err" , &D_bar_err_noSyst );
 chain_noSyst -> SetBranchAddress("S_bar_err" , &S_bar_err_noSyst );

 TFile* output = new TFile("tdfit_deltaPulls_noChol.root","RECREATE");
 TTree* new_tree = new TTree("MinuitParameterSetNtp","MinuitParameterSetNtp");


 Double_t delta_pull_C;
 Double_t delta_pull_D;
 Double_t delta_pull_S;
 Double_t delta_pull_D_bar;
 Double_t delta_pull_S_bar;


 new_tree->Branch("delta_pull_C",&delta_pull_C,"delta_pull_C/D");
 new_tree->Branch("delta_pull_D",&delta_pull_D,"delta_pull_D/D");
 new_tree->Branch("delta_pull_S",&delta_pull_S,"delta_pull_S/D");
 new_tree->Branch("delta_pull_D_bar",&delta_pull_D_bar,"delta_pull_D_bar/D");
 new_tree->Branch("delta_pull_S_bar",&delta_pull_S_bar,"delta_pull_S_bar/D");



 int entries = chain_syst->GetEntries();
 
 for(int i=0; i < entries; i++){

        chain_syst->GetEntry(i);
        chain_noSyst->GetEntry(i);

 delta_pull_C = (C_mean_noSyst - C_mean_syst) / C_err_noSyst;
 delta_pull_D = (D_mean_noSyst - D_mean_syst) / D_err_noSyst;
 delta_pull_S = (S_mean_noSyst - S_mean_syst) / S_err_noSyst;
 delta_pull_D_bar = (D_bar_mean_noSyst - D_bar_mean_syst) / D_bar_err_noSyst;
 delta_pull_S_bar = (S_bar_mean_noSyst - S_bar_mean_syst) / S_bar_err_noSyst;

  new_tree->Fill();
 }

  new_tree->Write();
  output->Close();
}