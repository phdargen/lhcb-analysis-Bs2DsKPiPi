// Phaspace acceptance studies
// author: Philippe d'Argent, Matthieu Kecke
#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventList.h"
#include "TTree.h"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Mint/DalitzEventPattern.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "Mint/PlotSet.h"
#include "TChain.h"
#include "Mint/NamedParameter.h"
#include <TROOT.h>
#include <sstream>
#include "TProfile.h"
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"

using namespace std;
using namespace MINT;

int efficiencyPlots(){

    DalitzEventPattern pdg(531, -431, 321, 211, -211);

    TChain* tree=new TChain("DecayTree");
    tree->Add("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
   
    TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
    tree_gen->Add("../GenMC.root");

    HyperHistogram hist_weights("plots/weights.root");
    int dim = 2;

    double K[5]; 
    double pip[5]; 
    double pim[5]; 
    double Ds_Kp[5],Ds_Km[5],Ds_pim[5];
    double weight;

    tree->SetBranchAddress("weight",&weight);
    
    tree->SetBranchAddress("K_plus_TRUEP_X",&K[0]);
    tree->SetBranchAddress("K_plus_TRUEP_Y",&K[1]);
    tree->SetBranchAddress("K_plus_TRUEP_Z",&K[2]); 
    tree->SetBranchAddress("K_plus_TRUEP_E",&K[3]); 
    tree->SetBranchAddress("K_plus_TRUEPT",&K[4]); 
	
    tree->SetBranchAddress("pi_plus_TRUEP_X",&pip[0]);
    tree->SetBranchAddress("pi_plus_TRUEP_Y",&pip[1]);
    tree->SetBranchAddress("pi_plus_TRUEP_Z",&pip[2]); 
    tree->SetBranchAddress("pi_plus_TRUEP_E",&pip[3]); 
    tree->SetBranchAddress("pi_plus_TRUEPT",&pip[4]); 

    tree->SetBranchAddress("pi_minus_TRUEP_X",&pim[0]);
    tree->SetBranchAddress("pi_minus_TRUEP_Y",&pim[1]);
    tree->SetBranchAddress("pi_minus_TRUEP_Z",&pim[2]); 
    tree->SetBranchAddress("pi_minus_TRUEP_E",&pim[3]); 
    tree->SetBranchAddress("pi_minus_TRUEPT",&pim[4]); 
	
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_X",&Ds_Kp[0]);
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_Y",&Ds_Kp[1]);
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_Z",&Ds_Kp[2]); 
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_E",&Ds_Kp[3]); 
    tree->SetBranchAddress("K_plus_fromDs_TRUEPT",&Ds_Kp[4]); 
    
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_X",&Ds_Km[0]);
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_Y",&Ds_Km[1]);
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_Z",&Ds_Km[2]); 
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_E",&Ds_Km[3]); 
    tree->SetBranchAddress("K_minus_fromDs_TRUEPT",&Ds_Km[4]); 

    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_X",&Ds_pim[0]);
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_Y",&Ds_pim[1]);
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_Z",&Ds_pim[2]); 
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_E",&Ds_pim[3]); 
    tree->SetBranchAddress("pi_minus_fromDs_TRUEPT",&Ds_pim[4]); 

            vector<int> s123;
            s123.push_back(1);
            s123.push_back(2);
            s123.push_back(3);
            
            vector<int> s234;
            s234.push_back(2);
            s234.push_back(3);
            s234.push_back(4);
            
            vector<int> s134;
            s134.push_back(1);
            s134.push_back(3);
            s134.push_back(4);
            
            vector<int> s124;
            s124.push_back(1);
            s124.push_back(2);
            s124.push_back(4);

	int nBins = 25;
        TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.8,4);
        TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ", nBins,0.,2);
        TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,2);
            TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
            TH1D* s_Dspi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            TH1D* s_Dspim = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            
    TH1D* min_pt = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);
        TH1D* min_Ds_pt = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);


    int numEvents = tree->GetEntries();
    //loop over tree and fill eventList
    for(int i=0; i< numEvents; i++)
    {
	tree->GetEntry(i);
        
        // Lorentz vectors: P=(Px,Py,Pz,E)
        TLorentzVector K_p(K[0],K[1],K[2],K[3]);
        TLorentzVector pip_p(pip[0],pip[1],pip[2],pip[3]);
	TLorentzVector pim_p(pim[0],pim[1],pim[2],pim[3]);
        TLorentzVector D_Kp_p(Ds_Kp[0],Ds_Kp[1],Ds_Kp[2],Ds_Kp[3]);
        TLorentzVector D_Km_p(Ds_Km[0],Ds_Km[1],Ds_Km[2],Ds_Km[3]);
        TLorentzVector D_pim_p(Ds_pim[0],Ds_pim[1],Ds_pim[2],Ds_pim[3]);
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);

	DalitzEvent evt(pdg, vectorOfvectors);
        evt.setWeight(weight);
	//eventList.Add(evt); // this fills the event list	

	s_Kpipi->Fill((evt.sij(s234)/(GeV*GeV)),evt.getWeight());
        s_Kpi->Fill((evt.s(2,4)/(GeV*GeV)),evt.getWeight());
        s_pipi->Fill((evt.s(3,4)/(GeV*GeV)),evt.getWeight());
                s_Dspipi->Fill(evt.sij(s134)/(GeV*GeV),evt.getWeight());
                s_DsK->Fill(evt.s(1,2)/(GeV*GeV),evt.getWeight());
                s_DsKpi->Fill(evt.sij(s124)/(GeV*GeV),evt.getWeight());
                s_Dspi->Fill(evt.s(1,3)/(GeV*GeV),evt.getWeight());
                s_Dspim->Fill(evt.s(1,4)/(GeV*GeV),evt.getWeight());

	min_pt->Fill( min(K[4],min(pip[4],min(pim[4],min(Ds_Kp[4],min(Ds_Km[4],Ds_pim[4]))))), weight  );
	min_Ds_pt->Fill( min(Ds_Kp[4],min(Ds_Km[4],Ds_pim[4])), weight  );

    }

    double K_gen[5]; 
    double pip_gen[5]; 
    double pim_gen[5]; 
    double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];
    
    tree_gen->SetBranchAddress("Kplus_TRUEP_X",&K_gen[0]);
    tree_gen->SetBranchAddress("Kplus_TRUEP_Y",&K_gen[1]);
    tree_gen->SetBranchAddress("Kplus_TRUEP_Z",&K_gen[2]); 
    tree_gen->SetBranchAddress("Kplus_TRUEP_E",&K_gen[3]); 
    tree_gen->SetBranchAddress("Kplus_TRUEPT",&K_gen[4]); 
	
    tree_gen->SetBranchAddress("piplus_TRUEP_X",&pip_gen[0]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Y",&pip_gen[1]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Z",&pip_gen[2]); 
    tree_gen->SetBranchAddress("piplus_TRUEP_E",&pip_gen[3]); 
    tree_gen->SetBranchAddress("piplus_TRUEPT",&pip_gen[4]); 

    tree_gen->SetBranchAddress("piminus_TRUEP_X",&pim_gen[0]);
    tree_gen->SetBranchAddress("piminus_TRUEP_Y",&pim_gen[1]);
    tree_gen->SetBranchAddress("piminus_TRUEP_Z",&pim_gen[2]); 
    tree_gen->SetBranchAddress("piminus_TRUEP_E",&pim_gen[3]); 
    tree_gen->SetBranchAddress("piminus_TRUEPT",&pim_gen[4]); 
	
    tree_gen->SetBranchAddress("Kplus0_TRUEP_X",&Ds_Kp_gen[0]);
    tree_gen->SetBranchAddress("Kplus0_TRUEP_Y",&Ds_Kp_gen[1]);
    tree_gen->SetBranchAddress("Kplus0_TRUEP_Z",&Ds_Kp_gen[2]); 
    tree_gen->SetBranchAddress("Kplus0_TRUEP_E",&Ds_Kp_gen[3]); 
    tree_gen->SetBranchAddress("Kplus0_TRUEPT",&Ds_Kp_gen[4]); 
    
    tree_gen->SetBranchAddress("Kminus_TRUEP_X",&Ds_Km_gen[0]);
    tree_gen->SetBranchAddress("Kminus_TRUEP_Y",&Ds_Km_gen[1]);
    tree_gen->SetBranchAddress("Kminus_TRUEP_Z",&Ds_Km_gen[2]); 
    tree_gen->SetBranchAddress("Kminus_TRUEP_E",&Ds_Km_gen[3]); 
    tree_gen->SetBranchAddress("Kminus_TRUEPT",&Ds_Km_gen[4]); 

    tree_gen->SetBranchAddress("piminus0_TRUEP_X",&Ds_pim_gen[0]);
    tree_gen->SetBranchAddress("piminus0_TRUEP_Y",&Ds_pim_gen[1]);
    tree_gen->SetBranchAddress("piminus0_TRUEP_Z",&Ds_pim_gen[2]); 
    tree_gen->SetBranchAddress("piminus0_TRUEP_E",&Ds_pim_gen[3]); 
    tree_gen->SetBranchAddress("piminus0_TRUEPT",&Ds_pim_gen[4]); 

    TH1D* s_Kpipi_gen = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",nBins,0.8,4);
    TH1D* s_Kpi_gen = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV);Events (norm.) ", nBins,0.,2);
    TH1D* s_pipi_gen = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",nBins,0.,2);
    TH1D* s_Dspipi_gen = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK_gen = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsKpi_gen = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi_gen = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim_gen = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* min_pt_gen = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);
    TH1D* min_Ds_pt_gen = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);

    TH1D* s_Kpipi_gen_rw = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",nBins,0.8,4);
    TH1D* s_Kpi_gen_rw = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV);Events (norm.) ", nBins,0.,2);
    TH1D* s_pipi_gen_rw = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",nBins,0.,2);
    TH1D* s_Dspipi_gen_rw = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsK_gen_rw = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
    TH1D* s_DsKpi_gen_rw = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
    TH1D* s_Dspi_gen_rw = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* s_Dspim_gen_rw = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
    TH1D* min_pt_gen_rw = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);
    TH1D* min_Ds_pt_gen_rw = new TH1D("",";min pt (GeV);Events (norm.) ",nBins,0.,4000);


    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
	tree_gen->GetEntry(i);
        
	double min_pt = min(K_gen[4],min(pip_gen[4],min(pim_gen[4],min(Ds_Kp_gen[4],min(Ds_Km_gen[4],Ds_pim_gen[4])))));
	double max_pt = max(K_gen[4],max(pip_gen[4],max(pim_gen[4],max(Ds_Kp_gen[4],max(Ds_Km_gen[4],Ds_pim_gen[4])))));

        // Lorentz vectors: P=(Px,Py,Pz,E)
        TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
        TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
	TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
        TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
        TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
        TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);
	DalitzEvent evt(pdg, vectorOfvectors);

	double w = 1.;
        //evt.setWeight(w);
	//eventList.Add(evt); // this fills the event list	
	s_Kpipi_gen->Fill((evt.sij(s234)/(GeV*GeV)),evt.getWeight());
        s_Kpi_gen->Fill((evt.s(2,4)/(GeV*GeV)),evt.getWeight());
        s_pipi_gen->Fill((evt.s(3,4)/(GeV*GeV)),evt.getWeight());
        s_Dspipi_gen->Fill(evt.sij(s134)/(GeV*GeV),evt.getWeight());
        s_DsK_gen->Fill(evt.s(1,2)/(GeV*GeV),evt.getWeight());
        s_DsKpi_gen->Fill(evt.sij(s124)/(GeV*GeV),evt.getWeight());
        s_Dspi_gen->Fill(evt.s(1,3)/(GeV*GeV),evt.getWeight());
        s_Dspim_gen->Fill(evt.s(1,4)/(GeV*GeV),evt.getWeight());
	min_pt_gen->Fill( min_pt,w  );
	min_Ds_pt_gen->Fill( min(Ds_Kp_gen[4],min(Ds_Km_gen[4],Ds_pim_gen[4])),w );

        HyperPoint point( dim );
        point.at(0)= 	log( min_pt  );
        if(K_gen[4] == min_pt) point.at(1)= 1/2. * log( (K_gen[3]+K_gen[2])/(K_gen[3]-K_gen[2]));
        else if(pip_gen[4] == min_pt) point.at(1)= 1/2. * log( (pip_gen[3]+pip_gen[2])/(pip_gen[3]-pip_gen[2]));
        else if(pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (pim_gen[3]+pim_gen[2])/(pim_gen[3]-pim_gen[2]));
        else if(Ds_Kp_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Kp_gen[3]+Ds_Kp_gen[2])/(Ds_Kp_gen[3]-Ds_Kp_gen[2]));
        else if(Ds_Km_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Km_gen[3]+Ds_Km_gen[2])/(Ds_Km_gen[3]-Ds_Km_gen[2]));
        else if(Ds_pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_pim_gen[3]+Ds_pim_gen[2])/(Ds_pim_gen[3]-Ds_pim_gen[2]));

	w = 0.;
        int bin = hist_weights.getBinning().getBinNum(point);
        if(hist_weights.checkBinNumber(bin)!= bin){
		 w = 0; //? should't happen
		 cout << "ERROR:: Event outside limits" << endl;	
		cout << point.at(0) << endl;
		cout << point.at(1) << endl << endl;
        }else w = hist_weights.getBinContent(bin);
 
	if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) w = 0.;
	if( min_pt < 100  ) w = 0.;
	if( max_pt < 1700  ) w = 0.;

	if( min(D_Kp_p.P(),min(D_Km_p.P(),D_pim_p.P())) < 1000 ) w = 0.;
	if( D_Kp_p.Pt()+D_Km_p.Pt()+D_pim_p.Pt() < 1800 ) w = 0.;

	if( min(K_p.P(),min(pip_p.P(),pim_p.P())) < 2000 ) w = 0.;
	if( K_p.Pt()+pip_p.Pt()+pim_p.Pt() < 1250 ) w = 0.;

	int pt_counter = 0;
	if(K_p.Pt() > 300)pt_counter++;
	if(pip_p.Pt() > 300)pt_counter++;
	if(pim_p.Pt() > 300)pt_counter++;
	if(pt_counter < 2) w = 0.;

        evt.setWeight(w);
	s_Kpipi_gen_rw->Fill((evt.sij(s234)/(GeV*GeV)),evt.getWeight());
        s_Kpi_gen_rw->Fill((evt.s(2,4)/(GeV*GeV)),evt.getWeight());
        s_pipi_gen_rw->Fill((evt.s(3,4)/(GeV*GeV)),evt.getWeight());
        s_Dspipi_gen_rw->Fill(evt.sij(s134)/(GeV*GeV),evt.getWeight());
        s_DsK_gen_rw->Fill(evt.s(1,2)/(GeV*GeV),evt.getWeight());
        s_DsKpi_gen_rw->Fill(evt.sij(s124)/(GeV*GeV),evt.getWeight());
        s_Dspi_gen_rw->Fill(evt.s(1,3)/(GeV*GeV),evt.getWeight());
        s_Dspim_gen_rw->Fill(evt.s(1,4)/(GeV*GeV),evt.getWeight());
	min_pt_gen_rw->Fill( min_pt,w  );
	min_Ds_pt_gen_rw->Fill( min(Ds_Kp_gen[4],min(Ds_Km_gen[4],Ds_pim_gen[4])),w );
    }

    TCanvas* c = new TCanvas();

    s_Kpipi->SetMinimum(0);
    s_Kpipi->DrawNormalized("e",1);
    s_Kpipi_gen->SetLineColor(kRed);
    s_Kpipi_gen->DrawNormalized("histsame",1);
    s_Kpipi_gen_rw->SetLineColor(kBlue);
    s_Kpipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Kpipi.eps");
    s_Kpipi->Scale(1./s_Kpipi->Integral());
    s_Kpipi_gen->Scale(1./s_Kpipi_gen->Integral());
    s_Kpipi_gen_rw->Scale(1./s_Kpipi_gen_rw->Integral());
    s_Kpipi->Divide(s_Kpipi,s_Kpipi_gen);
    s_Kpipi->Draw("e");
    s_Kpipi_gen_rw->Divide(s_Kpipi_gen_rw,s_Kpipi_gen);
    s_Kpipi_gen_rw->Draw("histsame");
    c->Print("eff_Kpipi.eps");

    s_Kpi->SetMinimum(0);
    s_Kpi->DrawNormalized("e",1);
    s_Kpi_gen->SetLineColor(kRed);
    s_Kpi_gen->DrawNormalized("histsame",1);
    s_Kpi_gen_rw->SetLineColor(kBlue);
    s_Kpi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Kpi.eps");
    s_Kpi->Scale(1./s_Kpi->Integral());
    s_Kpi_gen->Scale(1./s_Kpi_gen->Integral());
    s_Kpi_gen_rw->Scale(1./s_Kpi_gen_rw->Integral());
    s_Kpi->Divide(s_Kpi,s_Kpi_gen);
    s_Kpi->Draw("e");
    s_Kpi_gen_rw->Divide(s_Kpi_gen_rw,s_Kpi_gen);
    s_Kpi_gen_rw->Draw("histsame");
    c->Print("eff_Kpi.eps");

 s_pipi->SetMinimum(0);
    s_pipi->DrawNormalized("e",1);
    s_pipi_gen->SetLineColor(kRed);
    s_pipi_gen->DrawNormalized("histsame",1);
    s_pipi_gen_rw->SetLineColor(kBlue);
    s_pipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_pipi.eps");
    s_pipi->Scale(1./s_pipi->Integral());
    s_pipi_gen->Scale(1./s_pipi_gen->Integral());
    s_pipi_gen_rw->Scale(1./s_pipi_gen_rw->Integral());
    s_pipi->Divide(s_pipi,s_pipi_gen);
    s_pipi->Draw("e");
    s_pipi_gen_rw->Divide(s_pipi_gen_rw,s_pipi_gen);
    s_pipi_gen_rw->Draw("histsame");
    c->Print("eff_pipi.eps");

 s_Dspipi->SetMinimum(0);
    s_Dspipi->DrawNormalized("e",1);
    s_Dspipi_gen->SetLineColor(kRed);
    s_Dspipi_gen->DrawNormalized("histsame",1);
    s_Dspipi_gen_rw->SetLineColor(kBlue);
    s_Dspipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Dspipi.eps");
    s_Dspipi->Scale(1./s_Dspipi->Integral());
    s_Dspipi_gen->Scale(1./s_Dspipi_gen->Integral());
    s_Dspipi_gen_rw->Scale(1./s_Dspipi_gen_rw->Integral());
    s_Dspipi->Divide(s_Dspipi,s_Dspipi_gen);
    s_Dspipi->Draw("e");
    s_Dspipi_gen_rw->Divide(s_Dspipi_gen_rw,s_Dspipi_gen);
    s_Dspipi_gen_rw->Draw("histsame");
    c->Print("eff_Dspipi.eps");

 s_Dspi->SetMinimum(0);
    s_Dspi->DrawNormalized("e",1);
    s_Dspi_gen->SetLineColor(kRed);
    s_Dspi_gen->DrawNormalized("histsame",1);
    s_Dspi_gen_rw->SetLineColor(kBlue);
    s_Dspi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Dspi.eps");
    s_Dspi->Scale(1./s_Dspi->Integral());
    s_Dspi_gen->Scale(1./s_Dspi_gen->Integral());
    s_Dspi_gen_rw->Scale(1./s_Dspi_gen_rw->Integral());
    s_Dspi->Divide(s_Dspi,s_Dspi_gen);
    s_Dspi->Draw("e");
    s_Dspi_gen_rw->Divide(s_Dspi_gen_rw,s_Dspi_gen);
    s_Dspi_gen_rw->Draw("histsame");
    c->Print("eff_Dspi.eps");

 s_Dspim->SetMinimum(0);
    s_Dspim->DrawNormalized("e",1);
    s_Dspim_gen->SetLineColor(kRed);
    s_Dspim_gen->DrawNormalized("histsame",1);
    s_Dspim_gen_rw->SetLineColor(kBlue);
    s_Dspim_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Dspim.eps");
    s_Dspim->Scale(1./s_Dspim->Integral());
    s_Dspim_gen->Scale(1./s_Dspim_gen->Integral());
    s_Dspim_gen_rw->Scale(1./s_Dspim_gen_rw->Integral());
    s_Dspim->Divide(s_Dspim,s_Dspim_gen);
    s_Dspim->Draw("e");
    s_Dspim_gen_rw->Divide(s_Dspim_gen_rw,s_Dspim_gen);
    s_Dspim_gen_rw->Draw("histsame");
    c->Print("eff_Dspim.eps");

 s_DsK->SetMinimum(0);
    s_DsK->DrawNormalized("e",1);
    s_DsK_gen->SetLineColor(kRed);
    s_DsK_gen->DrawNormalized("histsame",1);
    s_DsK_gen_rw->SetLineColor(kBlue);
    s_DsK_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_DsK.eps");
    s_DsK->Scale(1./s_DsK->Integral());
    s_DsK_gen->Scale(1./s_DsK_gen->Integral());
    s_DsK_gen_rw->Scale(1./s_DsK_gen_rw->Integral());
    s_DsK->Divide(s_DsK,s_DsK_gen);
    s_DsK->Draw("e");
    s_DsK_gen_rw->Divide(s_DsK_gen_rw,s_DsK_gen);
    s_DsK_gen_rw->Draw("histsame");
    c->Print("eff_DsK.eps");

 s_DsKpi->SetMinimum(0);
    s_DsKpi->DrawNormalized("e",1);
    s_DsKpi_gen->SetLineColor(kRed);
    s_DsKpi_gen->DrawNormalized("histsame",1);
    s_DsKpi_gen_rw->SetLineColor(kBlue);
    s_DsKpi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_DsKpi.eps");
    s_DsKpi->Scale(1./s_DsKpi->Integral());
    s_DsKpi_gen->Scale(1./s_DsKpi_gen->Integral());
    s_DsKpi_gen_rw->Scale(1./s_DsKpi_gen_rw->Integral());
    s_DsKpi->Divide(s_DsKpi,s_DsKpi_gen);
    s_DsKpi->Draw("e");
    s_DsKpi_gen_rw->Divide(s_DsKpi_gen_rw,s_DsKpi_gen);
    s_DsKpi_gen_rw->Draw("histsame");
    c->Print("eff_DsKpi.eps");

min_pt->SetMinimum(0);
   min_pt->DrawNormalized("e",1);
   min_pt_gen->SetLineColor(kRed);
   min_pt_gen->DrawNormalized("histsame",1);
   min_pt_gen_rw->SetLineColor(kBlue);
   min_pt_gen_rw->DrawNormalized("histsame",1);
    c->Print("min_pt.eps");
   min_pt->Scale(1./min_pt->Integral());
   min_pt_gen->Scale(1./min_pt_gen->Integral());
   min_pt_gen_rw->Scale(1./min_pt_gen_rw->Integral());
   min_pt->Divide(min_pt,min_pt_gen);
   min_pt->Draw("e");
   min_pt_gen_rw->Divide(min_pt_gen_rw,min_pt_gen);
   min_pt_gen_rw->Draw("histsame");
    c->Print("eff_min_pt.eps");

min_Ds_pt->SetMinimum(0);
   min_Ds_pt->DrawNormalized("e",1);
   min_Ds_pt_gen->SetLineColor(kRed);
   min_Ds_pt_gen->DrawNormalized("histsame",1);
   min_Ds_pt_gen_rw->SetLineColor(kBlue);
   min_Ds_pt_gen_rw->DrawNormalized("histsame",1);
    c->Print("min_Ds_pt.eps");
   min_Ds_pt->Scale(1./min_Ds_pt->Integral());
   min_Ds_pt_gen->Scale(1./min_Ds_pt_gen->Integral());
   min_Ds_pt_gen_rw->Scale(1./min_Ds_pt_gen_rw->Integral());
   min_Ds_pt->Divide(min_Ds_pt,min_Ds_pt_gen);
   min_Ds_pt->Draw("e");
   min_Ds_pt_gen_rw->Divide(min_Ds_pt_gen_rw,min_Ds_pt_gen);
   min_Ds_pt_gen_rw->Draw("histsame");
    c->Print("eff_min_Ds_pt.eps");


}


void produceCorrectionHisto(){

    /// Options
    DalitzEventPattern pdg(531, -431, 321, 211, -211);
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 100); 
    NamedParameter<int> maxBinsPerDim("maxBinsPerDim", 200); 

    /// Get dimension and minimum bin width
    const int dim = 2;
   
    HyperPoint Min(4.,1.79);
    HyperPoint Max(11.,5.21);
    HyperCuboid limits(Min, Max );

    /// Get data           
    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/MC/signal.root");
    std::string inputFile = InputFileName;
    TFile *File =  TFile::Open(inputFile.c_str());
    TTree* tree =dynamic_cast<TTree*>(File->Get("DecayTree"));
    cout << "reading events from file " << inputFile.c_str() << endl;
    cout << " I've got " << tree->GetEntries() << " events." << endl;
   
    double K[5]; 
    double pip[5]; 
    double pim[5]; 
    double Ds_Kp[5],Ds_Km[5],Ds_pim[5];
    double weight;

    tree->SetBranchAddress("weight",&weight);
    
    tree->SetBranchAddress("K_plus_TRUEP_X",&K[0]);
    tree->SetBranchAddress("K_plus_TRUEP_Y",&K[1]);
    tree->SetBranchAddress("K_plus_TRUEP_Z",&K[2]); 
    tree->SetBranchAddress("K_plus_TRUEP_E",&K[3]); 
    tree->SetBranchAddress("K_plus_TRUEPT",&K[4]); 
	
    tree->SetBranchAddress("pi_plus_TRUEP_X",&pip[0]);
    tree->SetBranchAddress("pi_plus_TRUEP_Y",&pip[1]);
    tree->SetBranchAddress("pi_plus_TRUEP_Z",&pip[2]); 
    tree->SetBranchAddress("pi_plus_TRUEP_E",&pip[3]); 
    tree->SetBranchAddress("pi_plus_TRUEPT",&pip[4]); 

    tree->SetBranchAddress("pi_minus_TRUEP_X",&pim[0]);
    tree->SetBranchAddress("pi_minus_TRUEP_Y",&pim[1]);
    tree->SetBranchAddress("pi_minus_TRUEP_Z",&pim[2]); 
    tree->SetBranchAddress("pi_minus_TRUEP_E",&pim[3]); 
    tree->SetBranchAddress("pi_minus_TRUEPT",&pim[4]); 
	
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_X",&Ds_Kp[0]);
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_Y",&Ds_Kp[1]);
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_Z",&Ds_Kp[2]); 
    tree->SetBranchAddress("K_plus_fromDs_TRUEP_E",&Ds_Kp[3]); 
    tree->SetBranchAddress("K_plus_fromDs_TRUEPT",&Ds_Kp[4]); 
    
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_X",&Ds_Km[0]);
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_Y",&Ds_Km[1]);
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_Z",&Ds_Km[2]); 
    tree->SetBranchAddress("K_minus_fromDs_TRUEP_E",&Ds_Km[3]); 
    tree->SetBranchAddress("K_minus_fromDs_TRUEPT",&Ds_Km[4]); 

    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_X",&Ds_pim[0]);
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_Y",&Ds_pim[1]);
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_Z",&Ds_pim[2]); 
    tree->SetBranchAddress("pi_minus_fromDs_TRUEP_E",&Ds_pim[3]); 
    tree->SetBranchAddress("pi_minus_fromDs_TRUEPT",&Ds_pim[4]); 


    HyperPointSet points( dim );
    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);

        HyperPoint point( dim );
	double min_pt = min(K[4],min(pip[4],min(pim[4],min(Ds_Kp[4],min(Ds_Km[4],Ds_pim[4])))));
        point.at(0)= 	log( min_pt  );
        if(K[4] == min_pt) point.at(1)= 1/2. * log( (K[3]+K[2])/(K[3]-K[2]));
        else if(pip[4] == min_pt) point.at(1)= 1/2. * log( (pip[3]+pip[2])/(pip[3]-pip[2]));
        else if(pim[4] == min_pt) point.at(1)= 1/2. * log( (pim[3]+pim[2])/(pim[3]-pim[2]));
        else if(Ds_Kp[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Kp[3]+Ds_Kp[2])/(Ds_Kp[3]-Ds_Kp[2]));
        else if(Ds_Km[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Km[3]+Ds_Km[2])/(Ds_Km[3]-Ds_Km[2]));
        else if(Ds_pim[4] == min_pt) point.at(1)= 1/2. * log( (Ds_pim[3]+Ds_pim[2])/(Ds_pim[3]-Ds_pim[2]));

        point.addWeight(weight);
        points.push_back(point);
    }
    
    /// Get MC
    TChain* treeMC =new TChain("MCDecayTreeTuple/MCDecayTree");
    treeMC->Add("../GenMC.root");

    double K_gen[5]; 
    double pip_gen[5]; 
    double pim_gen[5]; 
    double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];
    
    treeMC->SetBranchAddress("Kplus_TRUEP_X",&K_gen[0]);
    treeMC->SetBranchAddress("Kplus_TRUEP_Y",&K_gen[1]);
    treeMC->SetBranchAddress("Kplus_TRUEP_Z",&K_gen[2]); 
    treeMC->SetBranchAddress("Kplus_TRUEP_E",&K_gen[3]); 
    treeMC->SetBranchAddress("Kplus_TRUEPT",&K_gen[4]); 
	
    treeMC->SetBranchAddress("piplus_TRUEP_X",&pip_gen[0]);
    treeMC->SetBranchAddress("piplus_TRUEP_Y",&pip_gen[1]);
    treeMC->SetBranchAddress("piplus_TRUEP_Z",&pip_gen[2]); 
    treeMC->SetBranchAddress("piplus_TRUEP_E",&pip_gen[3]); 
    treeMC->SetBranchAddress("piplus_TRUEPT",&pip_gen[4]); 

    treeMC->SetBranchAddress("piminus_TRUEP_X",&pim_gen[0]);
    treeMC->SetBranchAddress("piminus_TRUEP_Y",&pim_gen[1]);
    treeMC->SetBranchAddress("piminus_TRUEP_Z",&pim_gen[2]); 
    treeMC->SetBranchAddress("piminus_TRUEP_E",&pim_gen[3]); 
    treeMC->SetBranchAddress("piminus_TRUEPT",&pim_gen[4]); 
	
    treeMC->SetBranchAddress("Kplus0_TRUEP_X",&Ds_Kp_gen[0]);
    treeMC->SetBranchAddress("Kplus0_TRUEP_Y",&Ds_Kp_gen[1]);
    treeMC->SetBranchAddress("Kplus0_TRUEP_Z",&Ds_Kp_gen[2]); 
    treeMC->SetBranchAddress("Kplus0_TRUEP_E",&Ds_Kp_gen[3]); 
    treeMC->SetBranchAddress("Kplus0_TRUEPT",&Ds_Kp_gen[4]); 
    
    treeMC->SetBranchAddress("Kminus_TRUEP_X",&Ds_Km_gen[0]);
    treeMC->SetBranchAddress("Kminus_TRUEP_Y",&Ds_Km_gen[1]);
    treeMC->SetBranchAddress("Kminus_TRUEP_Z",&Ds_Km_gen[2]); 
    treeMC->SetBranchAddress("Kminus_TRUEP_E",&Ds_Km_gen[3]); 
    treeMC->SetBranchAddress("Kminus_TRUEPT",&Ds_Km_gen[4]); 

    treeMC->SetBranchAddress("piminus0_TRUEP_X",&Ds_pim_gen[0]);
    treeMC->SetBranchAddress("piminus0_TRUEP_Y",&Ds_pim_gen[1]);
    treeMC->SetBranchAddress("piminus0_TRUEP_Z",&Ds_pim_gen[2]); 
    treeMC->SetBranchAddress("piminus0_TRUEP_E",&Ds_pim_gen[3]); 
    treeMC->SetBranchAddress("piminus0_TRUEPT",&Ds_pim_gen[4]); 

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    HyperPointSet points_MC( dim );
    for (int i = 0; i < treeMC->GetEntries(); i++){
    
        treeMC->GetEntry(i);
        
	double min_pt = min(K_gen[4],min(pip_gen[4],min(pim_gen[4],min(Ds_Kp_gen[4],min(Ds_Km_gen[4],Ds_pim_gen[4])))));
	double max_pt = max(K_gen[4],max(pip_gen[4],max(pim_gen[4],max(Ds_Kp_gen[4],max(Ds_Km_gen[4],Ds_pim_gen[4])))));

        // Lorentz vectors: P=(Px,Py,Pz,E)
        TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
        TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
	TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
        TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
        TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
        TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);
	DalitzEvent evt(pdg, vectorOfvectors);

	double w = 1.;
	if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) w = 0.;
	if( min_pt < 100  ) w = 0.;
	if( max_pt < 1700  ) w = 0.;

	if( min(D_Kp_p.P(),min(D_Km_p.P(),D_pim_p.P())) < 1000 ) w = 0.;
	if( D_Kp_p.Pt()+D_Km_p.Pt()+D_pim_p.Pt() < 1800 ) w = 0.;

	if( min(K_p.P(),min(pip_p.P(),pim_p.P())) < 2000 ) w = 0.;
	if( K_p.Pt()+pip_p.Pt()+pim_p.Pt() < 1250 ) w = 0.;

	int pt_counter = 0;
	if(K_p.Pt() > 300)pt_counter++;
	if(pip_p.Pt() > 300)pt_counter++;
	if(pim_p.Pt() > 300)pt_counter++;
	if(pt_counter < 2) w = 0.;

	if( w == 0.) continue; 

        HyperPoint point( dim );
        point.at(0)= 	log( min_pt  );
        if(K_gen[4] == min_pt) point.at(1)= 1/2. * log( (K_gen[3]+K_gen[2])/(K_gen[3]-K_gen[2]));
        else if(pip_gen[4] == min_pt) point.at(1)= 1/2. * log( (pip_gen[3]+pip_gen[2])/(pip_gen[3]-pip_gen[2]));
        else if(pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (pim_gen[3]+pim_gen[2])/(pim_gen[3]-pim_gen[2]));
        else if(Ds_Kp_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Kp_gen[3]+Ds_Kp_gen[2])/(Ds_Kp_gen[3]-Ds_Kp_gen[2]));
        else if(Ds_Km_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Km_gen[3]+Ds_Km_gen[2])/(Ds_Km_gen[3]-Ds_Km_gen[2]));
        else if(Ds_pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_pim_gen[3]+Ds_pim_gen[2])/(Ds_pim_gen[3]-Ds_pim_gen[2]));

 	//point.addWeight(weight);
        points_MC.push_back(point);
    }

    cout << points_MC.size() << endl;
    cout << points_MC.size()/(double)treeMC->GetEntries() << endl;

    /// Define binning based on MC
    HyperHistogram hist(limits, points, 
                         
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI, 
                         
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),    
                                             
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (HyperPoint(0.01,0.5)),
                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),
                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),
                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (0)
                         
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),
                         
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")
                         
                         );

    /// Draw binning
    if(dim < 3)hist.draw("plots/binning");

    /// Draw density
    HyperHistogram histMC( hist.getBinning() );
    histMC.fill(points_MC); 


    cout << histMC.integral() << endl;
    cout << histMC.integral()/(double)treeMC->GetEntries() << endl;


    hist.normalise(1);
    histMC.normalise(1);
       
    if(dim < 3)histMC.drawDensity("plots/density_gen" );
    if(dim < 3)hist.drawDensity("plots/density_mc" );

    /// Produce MC correction histo 
    hist.divide(histMC);
    if(dim < 3)hist.drawDensity("plots/weights" );
    hist.save("plots/weights.root" );
   
    File->Close();

  return;
}





int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

    produceCorrectionHisto();
//     efficiencyPlots();
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
