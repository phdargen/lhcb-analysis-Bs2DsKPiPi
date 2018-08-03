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

   TH1D* Plot_C_shift = new TH1D("Plot_C_shift","C_{init}-C_{fit}", 10, -0.4,0.4);
   TH1D* Plot_D_shift = new TH1D("Plot_D_shift","D_{init}-D_{fit}", 10, -0.4,0.4);
   TH1D* Plot_S_shift = new TH1D("Plot_S_shift","S_{init}-S_{fit}", 10, -0.4,0.4);
   TH1D* Plot_D_bar_shift = new TH1D("Plot_D_bar_shift","D_bar_{init}-D_bar_{fit}", 10, -0.4,0.4);
   TH1D* Plot_S_bar_shift = new TH1D("Plot_S_bar_shift","S_bar_{init}-S_bar_{fit}", 10, -0.4,0.4);

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


   gStyle->SetOptStat(11);
   gStyle->SetOptFit(111);
   TCanvas* c = new TCanvas();

   Plot_C_shift->Draw("E1");
   c->Print(("pull_result_new/C_shift_"+ parName + ".eps").c_str());

   Plot_D_shift->Draw("E1");
   c->Print(("pull_result_new/D_shift_"+ parName + ".eps").c_str());

   Plot_S_shift->Draw("E1");
   c->Print(("pull_result_new/S_shift_"+ parName + ".eps").c_str());

   Plot_D_bar_shift->Draw("E1");
   c->Print(("pull_result_new/D_bar_shift_"+ parName + ".eps").c_str());

   Plot_S_bar_shift->Draw("E1");
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


}

void pull::LoopSyst(string parName)
{
if (fChain == 0) return;

   TH1D* pullPlot_C0 = new TH1D("pullPlot_C","Pull C, varying c_{0}", 10, -3.,3.);
   TH1D* pullPlot_D0 = new TH1D("pullPlot_D","Pull D, varying c_{0}", 10, -3.,3.);
   TH1D* pullPlot_S0 = new TH1D("pullPlot_S","Pull S, varying c_{0}", 10, -3.,3.);
   TH1D* pullPlot_D_bar0 = new TH1D("pullPlot_D_bar","Pull D_bar, varying c_{0}", 10, -3.,3.);
   TH1D* pullPlot_S_bar0 = new TH1D("pullPlot_S_bar","Pull S_bar, varying c_{0}", 10, -3.,3.);

   TH1D* pullPlot_C1 = new TH1D(" pullPlot_C","Pull C, varying c_{1}", 10, -3.,3.);
   TH1D* pullPlot_D1 = new TH1D(" pullPlot_D","Pull D, varying c_{1}", 10, -3.,3.);
   TH1D* pullPlot_S1 = new TH1D(" pullPlot_S","Pull S, varying c_{1}", 10, -3.,3.);
   TH1D* pullPlot_D_bar1 = new TH1D(" pullPlot_D_bar","Pull D_bar, varying c_{1}", 11, -3.,3.);
   TH1D* pullPlot_S_bar1 = new TH1D(" pullPlot_S_bar","Pull S_bar, varying c_{1}", 11, -3.,3.);

   TH1D* pullPlot_C2 = new TH1D("pullPlot_C ","Pull C, varying c_{2}", 10, -3.,3.);
   TH1D* pullPlot_D2 = new TH1D("pullPlot_D ","Pull D, varying c_{2}", 10, -3.,3.);
   TH1D* pullPlot_S2 = new TH1D("pullPlot_S ","Pull S, varying c_{2}", 10, -3.,3.);
   TH1D* pullPlot_D_bar2 = new TH1D("pullPlot_D_bar ","Pull D_bar, varying c_{2}", 10, -3.,3.);
   TH1D* pullPlot_S_bar2 = new TH1D("pullPlot_S_bar ","Pull S_bar, varying c_{2}", 10, -3.,3.);

   TH1D* pullPlot_C3 = new TH1D(" pullPlot_C ","Pull C, varying c_{3}", 10, -3.,3.);
   TH1D* pullPlot_D3 = new TH1D(" pullPlot_D ","Pull D, varying c_{3}", 10, -3.,3.);
   TH1D* pullPlot_S3 = new TH1D(" pullPlot_S ","Pull S, varying c_{3}", 10, -3.,3.);
   TH1D* pullPlot_D_bar3 = new TH1D(" pullPlot_D_bar ","Pull D_bar, varying c_{3}", 10, -3.,3.);
   TH1D* pullPlot_S_bar3 = new TH1D(" pullPlot_S_bar ","Pull S_bar, varying c_{3}", 10, -3.,3.);

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

   pullPlot_C0->Fit(gaussian);
   SummaryFile << "C, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C1->Fit(gaussian);
   SummaryFile << "C, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C2->Fit(gaussian);
   SummaryFile << "C, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_C3->Fit(gaussian);
   SummaryFile << "C, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_D0->Fit(gaussian);
   SummaryFile << "D, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D1->Fit(gaussian);
   SummaryFile << "D, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D2->Fit(gaussian);
   SummaryFile << "D, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D3->Fit(gaussian);
   SummaryFile << "D, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_S0->Fit(gaussian);
   SummaryFile << "S, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S1->Fit(gaussian);
   SummaryFile << "S, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S2->Fit(gaussian);
   SummaryFile << "S, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S3->Fit(gaussian);
   SummaryFile << "S, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_D_bar0->Fit(gaussian);
   SummaryFile << "$\\bar{D}$, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar1->Fit(gaussian);
   SummaryFile << "$\\bar{D}$, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar2->Fit(gaussian);
   SummaryFile << "$\\bar{D}$, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar3->Fit(gaussian);
   SummaryFile << "$\\bar{D}$, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_3_"+ parName + "_Gaussfit.eps").c_str());


   SummaryFile << "\\hline" << "\n";


   pullPlot_S_bar0->Fit(gaussian);
   SummaryFile << "$\\bar{S}$, varying $\\c_{0}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar0->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_0_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar1->Fit(gaussian);
   SummaryFile << "$\\bar{S}$, varying $\\c_{1}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar1->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_1_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar2->Fit(gaussian);
   SummaryFile << "$\\bar{S}$, varying $\\c_{2}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar2->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_2_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar3->Fit(gaussian);
   SummaryFile << "$\\bar{S}$, varying $\\c_{3}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar3->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_3_"+ parName + "_Gaussfit.eps").c_str());



   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\end{tabular}" << "\n";
   SummaryFile << "\\label{table:Pulls_tDFit}" << "\n";
   SummaryFile << "\\end{table}";



}

void pull::LoopSyst_noChol(string parName)
{
   if (fChain == 0) return;

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

   pullPlot_C->Fit(gaussian);
   SummaryFile << "C & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_C->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D->Fit(gaussian);
   SummaryFile << "D & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S->Fit(gaussian);
   SummaryFile << "S & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar->Fit(gaussian);
   SummaryFile << "$\\bar{D}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_D_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar->Fit(gaussian);
   SummaryFile << "$\\bar{S}$ & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
   pullPlot_S_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_"+ parName + "_Gaussfit.eps").c_str());

   SummaryFile << "\\hline" << "\n";
   SummaryFile << "\\end{tabular}" << "\n";
   SummaryFile << "\\label{table:Pulls_tDFit}" << "\n";
   SummaryFile << "\\end{table}";


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
