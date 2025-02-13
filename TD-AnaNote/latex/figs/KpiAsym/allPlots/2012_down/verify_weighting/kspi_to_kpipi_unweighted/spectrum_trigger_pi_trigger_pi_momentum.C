void spectrum_trigger_pi_trigger_pi_momentum()
{
//=========Macro generated from canvas: c1_n2/c1_n2
//=========  (Mon Mar 19 23:46:06 2018) by ROOT version6.08/06
   TCanvas *c1_n2 = new TCanvas("c1_n2", "c1_n2",13,60,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1_n2->SetHighLightColor(2);
   c1_n2->Range(-19.25926,-0.03875851,132.5926,0.2034822);
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
   Double_t xAxis49[18] = {2, 6, 14, 18, 22, 26, 30, 34, 38, 40, 44, 48, 52, 56, 60, 75, 100, 125}; 
   
   TH1F *ks_to_kpipi_pi_spectra_from_momentum__49 = new TH1F("ks_to_kpipi_pi_spectra_from_momentum__49","ks_to_kpipi_pi_spectra_from_momentum",17, xAxis49);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(2,0.04962444);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(3,0.1203394);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(4,0.1219774);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(5,0.1145404);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(6,0.1026633);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(7,0.0895745);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(8,0.07712238);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(9,0.06814416);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(10,0.06025588);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(11,0.05131158);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(12,0.04362121);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(13,0.03682631);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(14,0.03124866);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(15,0.0204999);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(16,0.008979817);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(17,0.003270632);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinContent(18,0.05813208);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(2,0.000131036);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(3,0.0002885772);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(4,0.0002905345);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(5,0.0002815381);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(6,0.0002665419);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(7,0.0002489718);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(8,0.0002310192);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(9,0.0003071053);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(10,0.0002042009);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(11,0.0001884368);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(12,0.0001737428);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(13,0.0001596383);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(14,0.0001470529);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(15,6.15061e-05);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(16,3.153201e-05);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(17,1.902978e-05);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetBinError(18,0.0004011399);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetMinimum(0);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetMaximum(0.1913701);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetEntries(1662949);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetStats(0);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetFillColor(4);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetFillStyle(3004);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetLineColor(4);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetLineWidth(2);
   ks_to_kpipi_pi_spectra_from_momentum__49->SetMarkerStyle(22);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetTitle("#it{momentum} [GeV/#it{c}]");
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetTitle("Candidates (arb. units)");
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_pi_spectra_from_momentum__49->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_from_momentum__49->Draw("HIST");
   Double_t xAxis50[18] = {2, 6, 14, 18, 22, 26, 30, 34, 38, 40, 44, 48, 52, 56, 60, 75, 100, 125}; 
   
   TH1F *ks_to_kpipi_pi_spectra_to_momentum__50 = new TH1F("ks_to_kpipi_pi_spectra_to_momentum__50","ks_to_kpipi_pi_spectra_to_momentum",17, xAxis50);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(1,1.12808e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(2,0.0505552);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(3,0.1332267);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(4,0.136693);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(5,0.1253567);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(6,0.1085242);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(7,0.09113029);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(8,0.07560547);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(9,0.06522574);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(10,0.05609616);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(11,0.04575145);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(12,0.03709817);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(13,0.0299405);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(14,0.02412267);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(15,0.01444123);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(16,0.004926842);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(17,0.00129449);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinContent(18,0.01325508);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(1,8.863035e-07);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(2,4.195468e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(3,9.631815e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(4,9.75631e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(5,9.342997e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(6,8.693121e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(7,7.966065e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(8,7.255863e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(9,9.530968e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(10,6.249986e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(11,5.644363e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(12,5.082633e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(13,4.566064e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(14,4.098503e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(15,1.637566e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(16,7.408952e-06);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(17,3.797711e-06);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetBinError(18,6.076223e-05);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetEntries(1.570531e+07);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetStats(0);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetFillColor(2);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetFillStyle(3005);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetLineColor(2);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetLineWidth(2);
   ks_to_kpipi_pi_spectra_to_momentum__50->SetMarkerStyle(23);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetNdivisions(505);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetXaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetLabelOffset(0.01);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetTitleOffset(0.95);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetYaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetZaxis()->SetLabelFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetZaxis()->SetLabelSize(0.06);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetZaxis()->SetTitleSize(0.072);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetZaxis()->SetTitleOffset(1.2);
   ks_to_kpipi_pi_spectra_to_momentum__50->GetZaxis()->SetTitleFont(132);
   ks_to_kpipi_pi_spectra_to_momentum__50->Draw("HIST SAME");
   
   TLegend *leg = new TLegend(0.555,0.73,0.89,0.92,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetTextSize(0.05);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("ks_to_kpipi_pi_spectra_from_momentum","[#it{D^{+}#rightarrowK_{s}^{0}#pi^{+}}] #pi^{+}","f");
   entry->SetFillColor(4);
   entry->SetFillStyle(3004);
   entry->SetLineColor(4);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("ks_to_kpipi_pi_spectra_to_momentum","[#it{D^{#pm}#rightarrowK^{#mp}#pi^{#pm}#pi^{#pm}}] #pi^{#pm}","f");
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
