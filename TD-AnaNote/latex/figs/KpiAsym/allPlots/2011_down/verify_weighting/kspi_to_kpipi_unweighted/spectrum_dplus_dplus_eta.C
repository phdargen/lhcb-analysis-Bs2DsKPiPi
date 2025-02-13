void spectrum_dplus_dplus_eta()
{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Thu Apr  5 16:15:30 2018) by ROOT version6.08/06
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",13,60,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(1.346914,-0.06650582,5.297531,0.3491556);
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
   Double_t xAxis31[9] = {1.9, 2.4, 2.65, 3, 3.3, 3.5, 3.8, 4.5, 5.1}; 
   
   TH1F *ks_to_kpipi_dplus_spectra_from_eta__31 = new TH1F("ks_to_kpipi_dplus_spectra_from_eta__31","ks_to_kpipi_dplus_spectra_from_eta",8, xAxis31);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(1,0.007661807);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(2,0.07674833);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(3,0.1668932);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(4,0.2173068);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(5,0.2330927);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(6,0.2103008);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(7,0.08503354);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinContent(8,0.00296282);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(1,8.900057e-05);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(2,0.0003983609);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(3,0.0004964751);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(4,0.0006119109);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(5,0.0007761784);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(6,0.0006019661);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(7,0.0002505869);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetBinError(8,5.052301e-05);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetMinimum(0);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetMaximum(0.3283725);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetEntries(614478);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetStats(0);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetFillColor(4);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetFillStyle(3004);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetLineColor(4);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetLineWidth(2);
   ks_to_kpipi_dplus_spectra_from_eta__31->SetMarkerStyle(22);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetTitle("#it{#eta}");
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetTitle("Candidates (arb. units)");
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_dplus_spectra_from_eta__31->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_from_eta__31->Draw("HIST");
   Double_t xAxis32[9] = {1.9, 2.4, 2.65, 3, 3.3, 3.5, 3.8, 4.5, 5.1}; 
   
   TH1F *ks_to_kpipi_dplus_spectra_to_eta__32 = new TH1F("ks_to_kpipi_dplus_spectra_to_eta__32","ks_to_kpipi_dplus_spectra_to_eta",8, xAxis32);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(1,0.02183996);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(2,0.1215061);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(3,0.204777);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(4,0.2345518);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(5,0.2141876);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(6,0.1586675);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(7,0.04392317);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinContent(8,0.0005469174);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(1,4.810518e-05);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(2,0.0001604649);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(3,0.0001760587);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(4,0.0002035212);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(5,0.0002381952);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(6,0.0001673918);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(7,5.765656e-05);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetBinError(8,6.949221e-06);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetEntries(5754128);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetStats(0);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetFillColor(2);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetFillStyle(3005);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetLineColor(2);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetLineWidth(2);
   ks_to_kpipi_dplus_spectra_to_eta__32->SetMarkerStyle(23);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_dplus_spectra_to_eta__32->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_dplus_spectra_to_eta__32->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.555,0.73,0.89,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("ks_to_kpipi_dplus_spectra_from_eta","[#it{D^{+}#rightarrowK_{s}^{0}#pi^{+}}] D^{+}","f");
   entry->SetFillColor(4);
   entry->SetFillStyle(3004);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("ks_to_kpipi_dplus_spectra_to_eta","[#it{D^{#pm}#rightarrowK^{#mp}#pi^{#pm}#pi^{#pm}}] D^{+}","f");
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
