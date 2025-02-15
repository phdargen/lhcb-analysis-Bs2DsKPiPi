void spectrum_trigger_pi_trigger_pi_pt()
{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Thu Apr  5 15:09:37 2018) by ROOT version6.08/06
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",13,60,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(-3456.79,-0.10369,21234.57,0.5443724);
   c1_n2->SetFillColor(0);
   c1_n2->SetBorderMode(0);
   c1_n2->SetBorderSize(2);
   c1_n2->SetTickx(1);
   c1_n2->SetTicky(1);
   c1_n2->SetLeftMargin(0.14);
   c1_n2->SetRightMargin(0.05);
   c1_n2->SetTopMargin(0.05);
   c1_n2->SetBottomMargin(0.16);
   c1_n2->SetFrameLineWidth(2);
   c1_n2->SetFrameBorderMode(0);
   c1_n2->SetFrameLineWidth(2);
   c1_n2->SetFrameBorderMode(0);
   Double_t xAxis25[14] = {0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 7500, 10000, 13000, 16000, 20000}; 
   
   TH1F *ks_to_kpipi_pi_spectra_from_pt__25 = new TH1F("ks_to_kpipi_pi_spectra_from_pt__25","ks_to_kpipi_pi_spectra_from_pt",13, xAxis25);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(3,4.643485e-06);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(4,0.3656924);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(5,0.3071082);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(6,0.1847682);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(7,0.09172343);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(8,0.03788222);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(9,0.01073773);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(10,0.001823828);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(11,0.0002542861);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinContent(12,5.085722e-06);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(3,1.755072e-06);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(4,0.000492528);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(5,0.0004513555);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(6,0.0003500956);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(7,0.0001744207);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(8,0.0001120923);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(9,3.774368e-05);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(10,1.555536e-05);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(11,5.302232e-06);
   ks_to_kpipi_pi_spectra_from_pt__25->SetBinError(12,7.498488e-07);
   ks_to_kpipi_pi_spectra_from_pt__25->SetMinimum(0);
   ks_to_kpipi_pi_spectra_from_pt__25->SetMaximum(0.5119693);
   ks_to_kpipi_pi_spectra_from_pt__25->SetEntries(1780581);
   ks_to_kpipi_pi_spectra_from_pt__25->SetStats(0);
   ks_to_kpipi_pi_spectra_from_pt__25->SetFillColor(4);
   ks_to_kpipi_pi_spectra_from_pt__25->SetFillStyle(3004);
   ks_to_kpipi_pi_spectra_from_pt__25->SetLineColor(4);
   ks_to_kpipi_pi_spectra_from_pt__25->SetLineWidth(2);
   ks_to_kpipi_pi_spectra_from_pt__25->SetMarkerStyle(22);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetTitle("#it{p_{T}} [MeV/#it{c}]");
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_from_pt__25->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetTitle("Candidates (arb. units)");
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_from_pt__25->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_pt__25->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_pt__25->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_pi_spectra_from_pt__25->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_pt__25->Draw("HIST");
   Double_t xAxis26[14] = {0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 7500, 10000, 13000, 16000, 20000}; 
   
   TH1F *ks_to_kpipi_pi_spectra_to_pt__26 = new TH1F("ks_to_kpipi_pi_spectra_to_pt__26","ks_to_kpipi_pi_spectra_to_pt",13, xAxis26);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(3,1.383222e-06);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(4,0.3144674);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(5,0.3144973);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(6,0.2047359);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(7,0.1076074);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(8,0.04529556);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(9,0.01211654);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(10,0.001230849);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(11,4.768475e-05);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinContent(12,2.426705e-08);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(3,3.173328e-07);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(4,0.0001513063);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(5,0.0001513134);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(6,0.000122086);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(7,6.258571e-05);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(8,4.060522e-05);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(9,1.32823e-05);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(10,4.233372e-06);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(11,7.60647e-07);
   ks_to_kpipi_pi_spectra_to_pt__26->SetBinError(12,1.715939e-08);
   ks_to_kpipi_pi_spectra_to_pt__26->SetEntries(1.657298e+07);
   ks_to_kpipi_pi_spectra_to_pt__26->SetStats(0);
   ks_to_kpipi_pi_spectra_to_pt__26->SetFillColor(2);
   ks_to_kpipi_pi_spectra_to_pt__26->SetFillStyle(3005);
   ks_to_kpipi_pi_spectra_to_pt__26->SetLineColor(2);
   ks_to_kpipi_pi_spectra_to_pt__26->SetLineWidth(2);
   ks_to_kpipi_pi_spectra_to_pt__26->SetMarkerStyle(23);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_to_pt__26->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_to_pt__26->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_pt__26->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_pt__26->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_pi_spectra_to_pt__26->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_pt__26->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.555,0.73,0.89,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("ks_to_kpipi_pi_spectra_from_pt","[#it{D^{+}#rightarrowK_{s}^{0}#pi^{+}}] #pi^{+}","f");
   entry->SetFillColor(4);
   entry->SetFillStyle(3004);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("ks_to_kpipi_pi_spectra_to_pt","[#it{D^{#pm}#rightarrowK^{#mp}#pi^{#pm}#pi^{#pm}}] #pi^{#pm}","f");
   entry->SetFillColor(2);
   entry->SetFillStyle(3005);
   entry->SetLineColor(2);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   leg->Draw();
   c1_n2->Modified();
   c1_n2->cd();
   c1_n2->SetSelected(c1_n2);
}
