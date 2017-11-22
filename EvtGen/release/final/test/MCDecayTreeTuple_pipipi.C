#define MCDecayTreeTuple_pipipi_cxx
#include "MCDecayTreeTuple_pipipi.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <algorithm>

void MCDecayTreeTuple_pipipi::Loop()
{

   if (fChain == 0) return;

    TFile* output = new TFile("Gen_pipipi.root","RECREATE");
    TTree* summary_tree = fChain->CloneTree(0);
    
    Long64_t nentries = fChain->GetEntries();
    cout << "Have " << nentries << " events" <<  endl;
    
    double m_DsK,m_Dspi,m_Dspipi,m_DsKpi,m_DsKpip,m_Kpipi,m_Kpi,m_pipi;
    summary_tree->Branch("m_DsK", &m_DsK, "m_DsK/D");
    summary_tree->Branch("m_Dspi", &m_Dspi, "m_Dspi/D");
    summary_tree->Branch("m_Dspipi", &m_Dspipi, "m_Dspipi/D");
    summary_tree->Branch("m_DsKpi", &m_DsKpi, "m_DsKpi/D");
    summary_tree->Branch("m_DsKpip", &m_DsKpip, "m_DsKpip/D");
    summary_tree->Branch("m_Kpipi", &m_Kpipi, "m_Kpipi/D");
    summary_tree->Branch("m_Kpi", &m_Kpi, "m_Kpi/D");
    summary_tree->Branch("m_pipi", &m_pipi, "m_pipi/D");
    
    double Ds_m12,Ds_m13;
    summary_tree->Branch("Ds_m12", &Ds_m12, "Ds_m12/D");
    summary_tree->Branch("Ds_m13", &Ds_m13, "Ds_m13/D");
    
    double Ds_m,Ds_SumPT,Xs_SumPT,Bs_P,minPT,minP;
    summary_tree->Branch("Ds_m", &Ds_m, "Ds_m/D");
    summary_tree->Branch("Ds_SumPT", &Ds_SumPT, "Ds_SumPT/D");
    summary_tree->Branch("Xs_SumPT", &Xs_SumPT, "Xs_SumPT/D");
    summary_tree->Branch("Bs_P", &Bs_P, "Bs_P/D");
    summary_tree->Branch("minPT", &minPT, "minPT/D");
    summary_tree->Branch("minP", &minP, "minP/D");

    TLorentzVector K_plus_fromDs;
    TLorentzVector K_minus_fromDs;
    TLorentzVector pi_minus_fromDs;
    TLorentzVector Ds;
    TLorentzVector K_plus;
    TLorentzVector pi_plus;
    TLorentzVector pi_minus;
    TLorentzVector pi_plus_fromDs;
    TLorentzVector pi_minus2_fromDs;
    TLorentzVector pi_plus1;
    TLorentzVector pi_plus2;
    
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);   
      
       K_plus.SetXYZM(Kplus_TRUEP_X,Kplus_TRUEP_Y,Kplus_TRUEP_Z,massKaon);
       pi_minus.SetXYZM(piminus_TRUEP_X,piminus_TRUEP_Y,piminus_TRUEP_Z,massPion);
       pi_plus.SetXYZM(piplus_TRUEP_X,piplus_TRUEP_Y,piplus_TRUEP_Z,massPion);   
      
       pi_minus_fromDs.SetXYZM(piminus0_TRUEP_X,piminus0_TRUEP_Y,piminus0_TRUEP_Z,massPion);        
       K_plus_fromDs.SetXYZM(piplus0_TRUEP_X,piplus0_TRUEP_Y,piplus0_TRUEP_Z,massPion);
       K_minus_fromDs.SetXYZM(piminus1_TRUEP_X,piminus1_TRUEP_Y,piminus1_TRUEP_Z,massPion);
       Ds = pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs;
      
       m_DsK = (Ds+K_plus).M();
       m_Dspi = (Ds+pi_plus).M();
       m_Dspipi = (Ds+pi_plus+pi_minus).M();
       m_DsKpi = (Ds+K_plus+pi_minus).M();
       m_DsKpip = (Ds+K_plus+pi_plus).M();
       m_Kpipi = (K_plus+pi_plus+pi_minus).M();
       m_Kpi = (K_plus+pi_minus).M();
       m_pipi = (pi_plus+pi_minus).M();  
       
       Ds_m12 = (K_plus_fromDs + K_minus_fromDs).M(); 
       Ds_m13 = (K_plus_fromDs + pi_minus_fromDs).M();   
       Ds_m = (Ds).M();   
 
       Ds_SumPT = piminus1_TRUEPT + piminus0_TRUEPT +  piplus0_TRUEPT ;
       Xs_SumPT = piminus_TRUEPT + piplus_TRUEPT +  Kplus_TRUEPT ;

	Bs_P = (Ds + K_plus + pi_plus + pi_minus).P();
	minPT = min(K_plus.Pt(),min(pi_minus.Pt(),min(pi_plus.Pt(),min(pi_minus_fromDs.Pt(),min(K_plus_fromDs.Pt(),K_minus_fromDs.Pt())))));
	minP = min(K_plus.P(),min(pi_minus.P(),min(pi_plus.P(),min(pi_minus_fromDs.P(),min(K_plus_fromDs.P(),K_minus_fromDs.P())))));

       summary_tree->Fill();
   }

    summary_tree->Draw("m_Kpi");
    summary_tree->Draw("m_pipi");
    summary_tree->Draw("m_Kpipi");
    
    summary_tree->Write();
    output->Close();
}
