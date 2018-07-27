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

   Plot_C_shift->SaveAs(("pull_result/C_shift_"+ parName +".root").c_str());
   Plot_D_shift->SaveAs(("pull_result/D_shift_"+ parName +".root").c_str());
   Plot_S_shift->SaveAs(("pull_result/S_shift_"+ parName +".root").c_str());
   Plot_D_bar_shift->SaveAs(("pull_result/D_bar_shift_"+ parName +".root").c_str());
   Plot_S_bar_shift->SaveAs(("pull_result/S_bar_shift_"+ parName +".root").c_str());

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

   gStyle->SetOptStat(11);
   gStyle->SetOptFit(111);
   TCanvas* c = new TCanvas();

   Plot_C_shift->Draw("E1");
   c->Print(("pull_result/C_shift_"+ parName + ".eps").c_str());

   Plot_D_shift->Draw("E1");
   c->Print(("pull_result/D_shift_"+ parName + ".eps").c_str());

   Plot_S_shift->Draw("E1");
   c->Print(("pull_result/S_shift_"+ parName + ".eps").c_str());

   Plot_D_bar_shift->Draw("E1");
   c->Print(("pull_result/D_bar_shift_"+ parName + ".eps").c_str());

   Plot_S_bar_shift->Draw("E1");
   c->Print(("pull_result/S_bar_shift_"+ parName + ".eps").c_str());


   pullPlot_C->Fit(gaussian);
   pullPlot_C->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_C_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D->Fit(gaussian);
   pullPlot_D->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S->Fit(gaussian);
   pullPlot_S->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_D_bar->Fit(gaussian);
   pullPlot_D_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_D_bar_"+ parName + "_Gaussfit.eps").c_str());

   pullPlot_S_bar->Fit(gaussian);
   pullPlot_S_bar->Draw("E1");
   gaussian->Draw("SAME");
   c->Print(("pull_result/pull_S_bar_"+ parName + "_Gaussfit.eps").c_str());

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


     shift_C = shift_C + TMath::Abs(C_init-C_mean);
     shift_D = shift_D + TMath::Abs(D_init-D_mean);
     shift_S = shift_S + TMath::Abs(S_init-S_mean);
     shift_D_bar = shift_D_bar + TMath::Abs(D_bar_init-D_bar_mean);
     shift_S_bar = shift_S_bar + TMath::Abs(S_bar_init-S_bar_mean);

     nen++;
   }

   cout << "mean shift for C_" << parName.c_str() << ":   " << shift_C/nen << endl;
   cout << "mean shift for D_" << parName.c_str() << ":   " << shift_D/nen << endl;
   cout << "mean shift for S_" << parName.c_str() << ":   " << shift_S/nen << endl;
   cout << "mean shift for D_bar_" << parName.c_str() << ":   " << shift_D_bar/nen << endl;
   cout << "mean shift for S_bar_" << parName.c_str() << ":   " << shift_S_bar/nen << endl;


   ofstream SummaryFile;
   SummaryFile.open(("pull_result/ShiftFile_" + parName + ".txt").c_str(),std::ofstream::trunc);
   SummaryFile << "C_shift = " << shift_C/nen  << "\n";
   SummaryFile << "D_shift = " << shift_D/nen  << "\n";
   SummaryFile << "S_shift = " << shift_S/nen  << "\n";
   SummaryFile << "D_bar_shift = " << shift_D_bar/nen  << "\n";
   SummaryFile << "S_bar_shift = " << shift_S_bar/nen;
}