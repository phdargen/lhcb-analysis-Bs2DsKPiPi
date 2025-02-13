void spectrum_piplus_piplus_eta()
{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Thu Apr  5 15:30:49 2018) by ROOT version6.08/06
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",13,60,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(1.346914,-0.07058323,5.297531,0.370562);
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
   Double_t xAxis7[9] = {1.9, 2.4, 2.65, 3, 3.3, 3.5, 3.8, 4.5, 5.1}; 
   
   TH1F *kpipi_to_signal_pi_spectra_from_eta__7 = new TH1F("kpipi_to_signal_pi_spectra_from_eta__7","kpipi_to_signal_pi_spectra_from_eta",8, xAxis7);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(0,0.0001333472);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(1,0.03739998);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(2,0.1677002);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(3,0.2489319);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(4,0.1891015);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(5,0.1843667);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(6,0.1254098);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(7,0.04539694);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(8,0.001693006);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinContent(9,5.568256e-06);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(0,2.228686e-06);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(1,0.0001382651);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(2,0.0004100057);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(3,0.0003733279);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(4,0.0003027134);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(5,0.0003845346);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(6,0.000259765);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(7,0.0001192507);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(8,1.062752e-05);
   kpipi_to_signal_pi_spectra_from_eta__7->SetBinError(9,2.76823e-07);
   kpipi_to_signal_pi_spectra_from_eta__7->SetMinimum(0);
   kpipi_to_signal_pi_spectra_from_eta__7->SetMaximum(0.3485047);
   kpipi_to_signal_pi_spectra_from_eta__7->SetEntries(3.196529e+07);
   kpipi_to_signal_pi_spectra_from_eta__7->SetStats(0);
   kpipi_to_signal_pi_spectra_from_eta__7->SetFillColor(4);
   kpipi_to_signal_pi_spectra_from_eta__7->SetFillStyle(3004);
   kpipi_to_signal_pi_spectra_from_eta__7->SetLineColor(4);
   kpipi_to_signal_pi_spectra_from_eta__7->SetLineWidth(2);
   kpipi_to_signal_pi_spectra_from_eta__7->SetMarkerStyle(22);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetTitle("#it{#eta}");
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetNdivisions(505);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetLabelOffset(0.01);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetTitleOffset(0.95);
   kpipi_to_signal_pi_spectra_from_eta__7->GetXaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetTitle("Candidates (arb. units)");
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetLabelOffset(0.01);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetTitleOffset(0.95);
   kpipi_to_signal_pi_spectra_from_eta__7->GetYaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->GetZaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->GetZaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_from_eta__7->GetZaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_from_eta__7->GetZaxis()->SetTitleOffset(1.2);
   kpipi_to_signal_pi_spectra_from_eta__7->GetZaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_from_eta__7->Draw("HIST");
   Double_t xAxis8[9] = {1.9, 2.4, 2.65, 3, 3.3, 3.5, 3.8, 4.5, 5.1}; 
   
   TH1F *kpipi_to_signal_pi_spectra_to_eta__8 = new TH1F("kpipi_to_signal_pi_spectra_to_eta__8","kpipi_to_signal_pi_spectra_to_eta",8, xAxis8);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(0,-0.0003807946);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(1,0.06531555);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(2,0.1648966);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(3,0.2324655);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(4,0.1465496);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(5,0.2003434);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(6,0.12209);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(7,0.06203899);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinContent(8,0.006300311);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(0,0.0003807946);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(1,0.0216781);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(2,0.05313266);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(3,0.05242182);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(4,0.05299425);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(5,0.07209606);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(6,0.03910497);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(7,0.01914817);
   kpipi_to_signal_pi_spectra_to_eta__8->SetBinError(8,0.005285425);
   kpipi_to_signal_pi_spectra_to_eta__8->SetEntries(535);
   kpipi_to_signal_pi_spectra_to_eta__8->SetStats(0);
   kpipi_to_signal_pi_spectra_to_eta__8->SetFillColor(2);
   kpipi_to_signal_pi_spectra_to_eta__8->SetFillStyle(3005);
   kpipi_to_signal_pi_spectra_to_eta__8->SetLineColor(2);
   kpipi_to_signal_pi_spectra_to_eta__8->SetLineWidth(2);
   kpipi_to_signal_pi_spectra_to_eta__8->SetMarkerStyle(23);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetNdivisions(505);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetLabelOffset(0.01);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetTitleOffset(0.95);
   kpipi_to_signal_pi_spectra_to_eta__8->GetXaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetLabelOffset(0.01);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetTitleOffset(0.95);
   kpipi_to_signal_pi_spectra_to_eta__8->GetYaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->GetZaxis()->SetLabelFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->GetZaxis()->SetLabelSize(0.06);
   kpipi_to_signal_pi_spectra_to_eta__8->GetZaxis()->SetTitleSize(0.072);
   kpipi_to_signal_pi_spectra_to_eta__8->GetZaxis()->SetTitleOffset(1.2);
   kpipi_to_signal_pi_spectra_to_eta__8->GetZaxis()->SetTitleFont(132);
   kpipi_to_signal_pi_spectra_to_eta__8->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.555,0.73,0.89,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("kpipi_to_signal_pi_spectra_from_eta","[#it{D^{#pm}#rightarrowK^{#mp}#pi^{#pm}#pi^{#pm}}] #pi^{#pm}","f");
   entry->SetFillColor(4);
   entry->SetFillStyle(3004);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("kpipi_to_signal_pi_spectra_to_eta","[#it{Signal}] #pi^{#pm}","f");
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
