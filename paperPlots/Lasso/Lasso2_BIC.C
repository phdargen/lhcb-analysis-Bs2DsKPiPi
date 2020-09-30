{
//=========Macro generated from canvas: c1/c1
//=========  (Thu Jul  2 21:35:22 2020) by ROOT version5.34/10
   TCanvas *c1 = new TCanvas("c1", "c1",10,32,700,500);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-1.576268,-31.24362,2.213091,101.529);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.14);
   c1->SetRightMargin(0.05);
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.16);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineWidth(2);
   c1->SetFrameBorderMode(0);
   
   TGraph *graph = new TGraph(40);
   graph->SetName("Graph0");
   graph->SetTitle("");
   graph->SetFillColor(1);
   graph->SetLineWidth(3);

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#0000ff");
   graph->SetMarkerColor(ci);
   graph->SetMarkerStyle(20);
   graph->SetMarkerSize(1.5);
   graph->SetPoint(0,0.1,8.998137404);
   graph->SetPoint(1,0.5,8.998137404);
   graph->SetPoint(2,1,8.998137404);
   graph->SetPoint(3,2,9.723722427);
   graph->SetPoint(4,3,10.69423394);
   graph->SetPoint(5,4,11.85440457);
   graph->SetPoint(6,5,13.17693655);
   graph->SetPoint(7,6,15.35288096);
   graph->SetPoint(8,7,16.92866606);
   graph->SetPoint(9,8,18.59733789);
   graph->SetPoint(10,9,2.494031344);
   graph->SetPoint(11,10,4.290410924);
   graph->SetPoint(12,11,6.139988172);
   graph->SetPoint(13,12,11.97350651);
   graph->SetPoint(14,13,15.12299547);
   graph->SetPoint(15,14,0);
   graph->SetPoint(16,15,2.562998228);
   graph->SetPoint(17,16,4.999080142);
   graph->SetPoint(18,17,7.340920342);
   graph->SetPoint(19,18,9.616979253);
   graph->SetPoint(20,19,11.84168199);
   graph->SetPoint(21,20,14.05435751);
   graph->SetPoint(22,22,18.56624692);
   graph->SetPoint(23,24,20.9591334);
   graph->SetPoint(24,25,22.52937366);
   graph->SetPoint(25,28,28.98402994);
   graph->SetPoint(26,30,13.58491434);
   graph->SetPoint(27,33,24.26794277);
   graph->SetPoint(28,39,27.31918075);
   graph->SetPoint(29,45,44.48514974);
   graph->SetPoint(30,50,53.9646776);
   graph->SetPoint(31,55,66.78365442);
   graph->SetPoint(32,60,71.39323007);
   graph->SetPoint(33,65,79.70280394);
   graph->SetPoint(34,70,50.90600433);
   graph->SetPoint(35,75,56.73498857);
   graph->SetPoint(36,81,65.51869643);
   graph->SetPoint(37,85,70.70530027);
   graph->SetPoint(38,93,82.18910966);
   graph->SetPoint(39,96,86.26397608);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0.09,165);
   Graph_Graph1->SetMinimum(-10);
   Graph_Graph1->SetMaximum(94.89037);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineWidth(2);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetTitle("#it{#lambda}");
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(132);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.005);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.065);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.9);
   Graph_Graph1->GetXaxis()->SetTitleFont(132);
   Graph_Graph1->GetYaxis()->SetTitle("#it{#DeltaBIC}");
   Graph_Graph1->GetYaxis()->SetLabelFont(132);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.005);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.065);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.9);
   Graph_Graph1->GetYaxis()->SetTitleFont(132);
   Graph_Graph1->GetZaxis()->SetLabelFont(132);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph1->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph1);
   
   graph->Draw("apl");
   
   TPaveText *pt = new TPaveText(0.2,0.8,0.5,0.9,"brNDC");
   pt->SetFillColor(0);
   pt->SetLineColor(0);
   pt->SetLineWidth(2);
   pt->SetBorderSize(0);
   pt->SetTextFont(132);
   pt->SetTextSize(0.08);
    pt->SetTextAlign(12);

   TText *text = pt->AddText("LHCb");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
    
    c1->Print("Lasso2_BIC.pdf");
    c1->Print("Lasso2_BIC.png");
    c1->Print("Lasso2_BIC.eps");

    c1->Print("Fig11b.pdf");
    c1->Print("Fig11b.png");
    c1->Print("Fig11b.eps");
    c1->Print("Fig11b.C");
}
