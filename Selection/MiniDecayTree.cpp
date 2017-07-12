#define MiniDecayTree_cxx
#include "MiniDecayTree.h"
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include "TLorentzVector.h"

using namespace std; 

inline Bool_t MiniDecayTree::PID_Cuts(){

    if(_Ds_finalState == Ds_finalState::phipi){
    
        if(K_plus_fromDs_PIDK < -10) return false;
        else if(K_minus_fromDs_PIDK < -10) return false;
        else if(pi_minus_fromDs_PIDK > 20) return false;

    }

    else if(_Ds_finalState == Ds_finalState::KsK){
        
        if(K_plus_fromDs_PIDK < -10) return false;
        else if(K_minus_fromDs_PIDK < -5) return false;
        else if(pi_minus_fromDs_PIDK > 20) return false;        
    }

    else if(_Ds_finalState == Ds_finalState::KKpi_NR){
        
        if(K_plus_fromDs_PIDK < 0) return false;
        else if(K_minus_fromDs_PIDK < 0) return false;
        else if(pi_minus_fromDs_PIDK > 10) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::pipipi){
        
        if(pi_plus_fromDs_PIDK > 10) return false;
        else if(pi_minus_fromDs_PIDK > 10) return false;
        else if(pi_minus2_fromDs_PIDK > 10) return false;
        
        if(pi_plus_fromDs_PIDp > 10) return false;
        else if(pi_minus_fromDs_PIDp > 10) return false;
        else if(pi_minus2_fromDs_PIDp > 10) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
        
        
    }
    

    if(_decay == Decay::signal){
        
        if(K_plus_PIDK < 5) return false;
        else if(pi_plus_PIDK > 10) return false;
        else if(pi_minus_PIDK > 10) return false;
        
    }

    else if(_decay == Decay::norm){

        if(pi_plus1_PIDK > 10) return false;
        else if(pi_plus2_PIDK > 10) return false;
        else if(pi_minus_PIDK > 10) return false;
        
    }

    
    return true;
}


inline Bool_t MiniDecayTree::Veto_Cuts(){

    if(_Ds_finalState == Ds_finalState::phipi || _Ds_finalState == Ds_finalState::KsK || _Ds_finalState == Ds_finalState::KKpi_NR){
        
        //D- veto
        if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus) < 30. && K_minus_fromDs_PIDK < 10 ) return false;
        
        //D0 veto
        if((K_plus_fromDs + K_minus_fromDs).M() > 1840.) return false;
        
        //Lambda_c veto
        if( TMath::Abs((K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c) < 30. && ((K_minus_fromDs_PIDK - K_minus_fromDs_PIDp) < 5) ) return false;
        
        // Charmless veto
        if(_Ds_finalState == Ds_finalState::phipi){
            if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < -1) return false;
            if(Ds_FDCHI2_ORIVX < 0) return false;
        }
        else if(_Ds_finalState == Ds_finalState::KsK){
            if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
            if(Ds_FDCHI2_ORIVX < 0) return false;
        }
        else if(_Ds_finalState == Ds_finalState::KKpi_NR){
            if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
            if(Ds_FDCHI2_ORIVX < 2) return false;
        }
    }

    else if(_Ds_finalState == Ds_finalState::pipipi){
               
        // D0 veto
        if((pi_plus_fromDs + pi_minus_fromDs).M() > 1700 || (pi_plus_fromDs + pi_minus2_fromDs).M() > 1700) return false;
        
        // Charmless veto
        if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < 0) return false;
        if(Ds_FDCHI2_ORIVX < 9) return false;
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
    }
    
    if(_decay == Decay::signal){
            
        // Ds veto
        if(TMath::Abs((K_plus + pi_plus + pi_minus).M() - massDs) < 20) return false;
        
    }
    
    else if(_decay == Decay::norm){
               
        // Ds veto
        if(TMath::Abs((pi_plus1 + pi_plus2 + pi_minus).M() - massDs) < 20) return false;

    }

    return true;
}

  
inline Bool_t MiniDecayTree::Preselection_Cuts(){
    
        if(fabs(Ds_MM-massDs) > 25 ) return false;
    
        return true;
}

inline Ds_finalState::Type MiniDecayTree::get_Ds_finalState(){

    if(_Ds_finalState == Ds_finalState::pipipi || _Ds_finalState == Ds_finalState::Kpipi) return _Ds_finalState;

    else {
    
        K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
    
        if( TMath::Abs(((K_plus_fromDs + K_minus_fromDs).M() - massPhi)) < 20) return Ds_finalState::phipi;
        
        else if (TMath::Abs(((pi_minus_fromDs + K_plus_fromDs).M() - massKstar)) < 75) return Ds_finalState::KsK;
        
        else return Ds_finalState::KKpi_NR;
    
    }
}

inline void MiniDecayTree::set_LorentzVectors(){

    if(_Ds_finalState == Ds_finalState::phipi || _Ds_finalState == Ds_finalState::KsK || _Ds_finalState == Ds_finalState::KKpi_NR){
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);        
        K_plus_fromDs.SetXYZM(K_plus_fromDs_PX,K_plus_fromDs_PY,K_plus_fromDs_PZ,massKaon);
        K_minus_fromDs.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massKaon);
        Ds = pi_minus_fromDs + K_plus_fromDs + K_minus_fromDs;
        
        Kminus_fromDs_asProton_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ, massProton);
        Kminus_fromDs_asPiminus_MissID.SetXYZM(K_minus_fromDs_PX,K_minus_fromDs_PY,K_minus_fromDs_PZ,massPion);
    }
    
    else if(_Ds_finalState == Ds_finalState::pipipi){
        pi_plus_fromDs.SetXYZM(pi_plus_fromDs_PX,pi_plus_fromDs_PY,pi_plus_fromDs_PZ,massPion);
        pi_minus_fromDs.SetXYZM(pi_minus_fromDs_PX,pi_minus_fromDs_PY,pi_minus_fromDs_PZ,massPion);        
        pi_minus2_fromDs.SetXYZM(pi_minus2_fromDs_PX,pi_minus2_fromDs_PY,pi_minus2_fromDs_PZ,massPion);
        Ds = pi_minus_fromDs + pi_plus_fromDs + pi_minus2_fromDs;
    }
    
    else if(_Ds_finalState == Ds_finalState::Kpipi){
    }
    
    if(_decay == Decay::signal){
        K_plus.SetXYZM(K_plus_PX,K_plus_PY,K_plus_PZ,massKaon);
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus.SetXYZM(pi_plus_PX,pi_plus_PY,pi_plus_PZ,massPion);        
    }
    
    else if(_decay == Decay::norm){
        pi_minus.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ,massPion);
        pi_plus1.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ,massPion);
        pi_plus2.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ,massPion);
    }
    
    return;
}


TTree* MiniDecayTree::GetInputTree(){
        
    cout << "Read from file " << _inFileName << endl;    
        
    TChain* chain= new TChain("DecayTree");
    chain->Add(_inFileName);
    //chain->Add("Signal/11D/*.root");
    if(chain==0){
        cout << "ERROR: No file found" << endl;
        throw "ERROR";
    }
    
    return (TTree*)chain;
}

void MiniDecayTree::Loop()
{

   time_t startTime = time(0);

   Init();
   if (fChain == 0) return;

   TFile* output = new TFile(_outFileName,"RECREATE");
   TTree* summary_tree = fChain->CloneTree(0);
    
   Long64_t nentries = fChain->GetEntries();
   cout << "Have " << nentries << " events" <<  endl;
   
    // Add new branches to output tree
    int Ds_finalState = 0 ;
    summary_tree->Branch("Ds_finalState",&Ds_finalState,"Ds_finalState/I");
    
    float DsDaughters_min_IPCHI2 = 0;
    float XsDaughters_min_IPCHI2 = 0;
    float DsDaughters_max_IPCHI2 = 0;
    float XsDaughters_max_IPCHI2 = 0;
    float DsDaughters_min_PT = 0;
    float XsDaughters_min_PT = 0;
    float Xs_max_DOCA = 0;
    float Ds_max_DOCA = 0;
    float max_TrackChi2 = 0;
    float min_TrackChi2 = 0;
    float Xs_max_ghostProb = 0;
    float Ds_max_ghostProb = 0;
    float max_ghostProb = 0;
    double Bs_ct = 0;
    double Bs_cterr = 0;
    double Ds_ct = 0;
    double Ds_cterr = 0;
    int qt = 0;
    int qf = 0;
    summary_tree->Branch("DsDaughters_min_IPCHI2",&DsDaughters_min_IPCHI2,"DsDaughters_min_IPCHI2/F");
    summary_tree->Branch("XsDaughters_min_IPCHI2",&XsDaughters_min_IPCHI2,"XsDaughters_min_IPCHI2/F");
    summary_tree->Branch("DsDaughters_max_IPCHI2",&DsDaughters_max_IPCHI2,"DsDaughters_max_IPCHI2/F");
    summary_tree->Branch("XsDaughters_max_IPCHI2",&XsDaughters_max_IPCHI2,"XsDaughters_max_IPCHI2/F");
    summary_tree->Branch("DsDaughters_min_PT",&DsDaughters_min_PT,"DsDaughters_min_PT/F");
    summary_tree->Branch("XsDaughters_min_PT",&XsDaughters_min_PT,"XsDaughters_min_PT/F");
    summary_tree->Branch("Xs_max_DOCA",&Xs_max_DOCA,"Xs_max_DOCA/F");
    summary_tree->Branch("Ds_max_DOCA",&Ds_max_DOCA,"Ds_max_DOCA/F");
    //summary_tree->Branch("max_TrackChi2",&max_TrackChi2,"max_Track_Chi2/F");
    //summary_tree->Branch("min_TrackChi2",&min_TrackChi2,"min_Track_Chi2/F");
    summary_tree->Branch("Xs_max_ghostProb",&Xs_max_ghostProb,"Xs_max_ghostProb/F");
    summary_tree->Branch("Ds_max_ghostProb",&Ds_max_ghostProb,"Ds_max_ghostProb/F");
    summary_tree->Branch("max_ghostProb",&max_ghostProb,"max_ghostProb/F");
    summary_tree->Branch("Bs_ct",&Bs_ct,"Bs_ct/D");
    summary_tree->Branch("Bs_cterr",&Bs_cterr,"Bs_cterr/D");
    summary_tree->Branch("Ds_ct",&Ds_ct,"Ds_ct/D");
    summary_tree->Branch("Ds_cterr",&Ds_cterr,"Ds_cterr/D");
    summary_tree->Branch("qt",&qt,"qt/I");
    summary_tree->Branch("qf",&qf,"qf/I");
    
    double angK,angPip, angPim, maxCos, maxGP;
    summary_tree->Branch("angK", &angK, "angK/D");
    summary_tree->Branch("angPip", &angPip, "angPip/D");
    summary_tree->Branch("angPim", &angPim, "angPim/D");
    summary_tree->Branch("maxCos", &maxCos, "maxCos/D");

    double Bs_RFD,Ds_RFD,Ds_FDsig,Ds_z;
    summary_tree->Branch("Bs_RFD", &Bs_RFD, "Bs_RFD/D");
    summary_tree->Branch("Ds_RFD", &Ds_RFD, "Ds_RFD/D");
    summary_tree->Branch("Ds_FDsig", &Ds_FDsig, "Ds_FDsig/D");
    summary_tree->Branch("Ds_z", &Ds_z, "Ds_z/D");
    
    double Ds_m12,Ds_m13;
    summary_tree->Branch("Ds_m12", &Ds_m12, "Ds_m12/D");
    summary_tree->Branch("Ds_m13", &Ds_m13, "Ds_m13/D");
    
    double bkg_D_as_Ds_m, bkg_Lambdac_as_Ds_m;
    summary_tree->Branch("bkg_D_as_Ds_m", &bkg_D_as_Ds_m, "bkg_D_as_Ds_m/D");
    summary_tree->Branch("bkg_Lambdac_as_Ds_m", &bkg_Lambdac_as_Ds_m, "bkg_Lambdac_as_Ds_m/D");
    
    double bkg_KKpi_as_Xs_m, bkg_3pi_as_Xs_m;
    summary_tree->Branch("bkg_KKpi_as_Xs_m", &bkg_KKpi_as_Xs_m, "bkg_KKpi_as_Xs_m/D");
    summary_tree->Branch("bkg_3pi_as_Xs_m", &bkg_3pi_as_Xs_m, "bkg_3pi_as_Xs_m/D");

    double m_DsK,m_Dspi,m_Dspipi,m_DsKpi,m_DsKpip;
    summary_tree->Branch("m_DsK", &m_DsK, "m_DsK/D");
    summary_tree->Branch("m_Dspi", &m_Dspi, "m_Dspi/D");
    summary_tree->Branch("m_Dspipi", &m_Dspipi, "m_Dspipi/D");
    summary_tree->Branch("m_DsKpi", &m_DsKpi, "m_DsKpi/D");
    summary_tree->Branch("m_DsKpip", &m_DsKpip, "m_DsKpip/D");
    
    double beta_K_plus;
    double beta_pi_plus;
    double beta_pi_minus;
    double beta_Ds;
    summary_tree->Branch("beta_K_plus", &beta_K_plus, "beta_K_plus/D");
    summary_tree->Branch("beta_pi_plus", &beta_pi_plus, "beta_pi_plus/D");
    summary_tree->Branch("beta_pi_minus", &beta_pi_minus, "beta_pi_minus/D");
    summary_tree->Branch("beta_Ds", &beta_Ds, "beta_Ds/D");

    double beta_K_plus_fromDs;
    double beta_pi_minus_fromDs;
    double beta_K_minus_fromDs;
    summary_tree->Branch("beta_K_plus_fromDs", &beta_K_plus_fromDs, "beta_K_plus_fromDs/D");
    summary_tree->Branch("beta_pi_minus_fromDs", &beta_pi_minus_fromDs, "beta_pi_minus_fromDs/D");
    summary_tree->Branch("beta_K_minus_fromDs", &beta_K_minus_fromDs, "beta_K_minus_fromDs/D");

    for (Long64_t i=0; i<nentries;i++) {
        
        if(0ul == (i % 10000ul)) cout << "Read event " << i << "/" << nentries <<
        "  ( " << i/(double)nentries * 100. << " % )" << endl;

        fChain->GetEntry(i);           
        //Long64_t j = LoadTree(i);
        //if (j < 0) break;
        
        // Apply preselection cust
        if(!Preselection_Cuts()) continue;
    
        _Ds_finalState = get_Ds_finalState();
        if(_data)if(!PID_Cuts()) continue;
    
        set_LorentzVectors();
        if(!Veto_Cuts()) continue;
        
        // Add new variables
        Ds_finalState = get_Ds_finalState();

        Bs_RFD = sqrt(pow(Bs_ENDVERTEX_X-Bs_OWNPV_X,2)+pow(Bs_ENDVERTEX_Y-Bs_OWNPV_Y,2));
        Ds_RFD = sqrt(pow(Ds_ENDVERTEX_X-Ds_OWNPV_X,2)+pow(Ds_ENDVERTEX_Y-Ds_OWNPV_Y,2));
        Ds_FDsig = (Ds_ENDVERTEX_Z-Ds_ORIVX_Z)/sqrt(pow(Ds_ENDVERTEX_ZERR,2)+pow(Ds_ORIVX_ZERR,2));
        Ds_z = Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z;
        
        if(_decay == Decay::signal){
            TVector3 v_Ds(Ds_PX,Ds_PY,0.);
            TVector3 v_K(K_plus_PX,K_plus_PY,0.);
            TVector3 v_pip(pi_plus_PX,pi_plus_PY,0.);
            TVector3 v_pim(pi_minus_PX,pi_minus_PY,0.);
            angK= v_Ds.Angle(v_K);
            angPip= v_Ds.Angle(v_pip);
            angPim= v_Ds.Angle(v_pim);
            maxCos = cos(max(angK,max(angPip,angPim)));
            
            XsDaughters_min_IPCHI2 = min(K_plus_IPCHI2_OWNPV,min(pi_minus_IPCHI2_OWNPV,pi_plus_IPCHI2_OWNPV));            
            XsDaughters_max_IPCHI2 = max(K_plus_IPCHI2_OWNPV,max(pi_minus_IPCHI2_OWNPV,pi_plus_IPCHI2_OWNPV));            
            XsDaughters_min_PT = min(K_plus_PT,min(pi_minus_PT,pi_plus_PT));          
            Xs_max_DOCA = max(K_1_1270_plus_DOCA1,max(K_1_1270_plus_DOCA2,K_1_1270_plus_DOCA3));
            Xs_max_ghostProb = max(pi_plus_TRACK_GhostProb,max(K_plus_TRACK_GhostProb,pi_minus_TRACK_GhostProb)); 
            
            m_DsK = (Ds+K_plus).M();
            m_Dspi = (Ds+pi_plus).M();
            m_Dspipi = (Ds+pi_plus+pi_minus).M();
            m_DsKpi = (Ds+K_plus+pi_minus).M();
            m_DsKpip = (Ds+K_plus+pi_plus).M();
            
            TLorentzVector pi_minus_asK_MissID; 
            pi_minus_asK_MissID.SetXYZM(pi_minus_PX,pi_minus_PY,pi_minus_PZ, massKaon);
            bkg_KKpi_as_Xs_m = (K_plus+pi_plus+pi_minus_asK_MissID).M();
            
            beta_K_plus = (-K_plus_P + pi_plus_P + pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_pi_plus = (K_plus_P - pi_plus_P + pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_pi_minus = (K_plus_P + pi_plus_P - pi_minus_P + Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);
            beta_Ds = (K_plus_P + pi_plus_P + pi_minus_P - Ds_P)/(K_plus_P + pi_plus_P + pi_minus_P + Ds_P);

        }
        else {
            TVector3 v_Ds(Ds_PX,Ds_PY,0.);
            TVector3 v_K(pi_plus1_PX,pi_plus1_PY,0.);
            TVector3 v_pip(pi_plus2_PX,pi_plus2_PY,0.);
            TVector3 v_pim(pi_minus_PX,pi_minus_PY,0.);
            angK= v_Ds.Angle(v_K);
            angPip= v_Ds.Angle(v_pip);
            angPim= v_Ds.Angle(v_pim);
            maxCos = cos(max(angK,max(angPip,angPim)));
            
            XsDaughters_min_IPCHI2 = min(pi_plus1_IPCHI2_OWNPV,min(pi_minus_IPCHI2_OWNPV,pi_plus2_IPCHI2_OWNPV));            
            XsDaughters_max_IPCHI2 = max(pi_plus1_IPCHI2_OWNPV,max(pi_minus_IPCHI2_OWNPV,pi_plus2_IPCHI2_OWNPV));            
            XsDaughters_min_PT = min(pi_plus1_PT,min(pi_minus_PT,pi_plus2_PT));          
            Xs_max_DOCA = max(a_1_1260_plus_DOCA1,max(a_1_1260_plus_DOCA2,a_1_1260_plus_DOCA3));
            Xs_max_ghostProb = max(pi_plus1_TRACK_GhostProb,max(pi_plus2_TRACK_GhostProb,pi_minus_TRACK_GhostProb));   
            
            m_DsK = (Ds+pi_plus1).M();
            m_Dspi = (Ds+pi_plus2).M();
            m_Dspipi = (Ds+pi_plus1+pi_minus).M();
            m_DsKpi = (Ds+pi_plus2+pi_minus).M();
            m_DsKpip = (Ds+pi_plus1+pi_plus2).M(); 
            
            TLorentzVector pi_plus1_asK_MissID; 
            TLorentzVector pi_plus2_asK_MissID; 

            pi_plus1_asK_MissID.SetXYZM(pi_plus1_PX,pi_plus1_PY,pi_plus1_PZ, massKaon);
            pi_plus2_asK_MissID.SetXYZM(pi_plus2_PX,pi_plus2_PY,pi_plus2_PZ, massKaon);

            if(pi_plus1_PIDK > pi_plus2_PIDK )bkg_3pi_as_Xs_m = (pi_minus+pi_plus2+pi_plus1_asK_MissID).M();
            else bkg_3pi_as_Xs_m = (pi_minus+pi_plus1+pi_plus2_asK_MissID).M();
            
            beta_K_plus = (-pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_pi_plus = (pi_plus1_P - pi_plus2_P + pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_pi_minus = (pi_plus1_P + pi_plus2_P - pi_minus_P + Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
            beta_Ds = (pi_plus1_P + pi_plus2_P + pi_minus_P - Ds_P)/(pi_plus1_P + pi_plus2_P + pi_minus_P + Ds_P);
        }
        
        if(_Ds_finalState == Ds_finalState::pipipi){
            DsDaughters_min_IPCHI2 = min(pi_plus_fromDs_IPCHI2_OWNPV,min(pi_minus_fromDs_IPCHI2_OWNPV,pi_minus2_fromDs_IPCHI2_OWNPV));            
            DsDaughters_max_IPCHI2 = max(pi_plus_fromDs_IPCHI2_OWNPV,max(pi_minus_fromDs_IPCHI2_OWNPV,pi_minus2_fromDs_IPCHI2_OWNPV));            
            DsDaughters_min_PT = min(pi_plus_fromDs_PT,min(pi_minus_fromDs_PT,pi_minus2_fromDs_PT));            
            Ds_max_DOCA = max(Ds_DOCA1,max(Ds_DOCA2,Ds_DOCA3));
            Ds_max_ghostProb = max(pi_plus_fromDs_TRACK_GhostProb,max(pi_minus2_fromDs_TRACK_GhostProb,pi_minus_fromDs_TRACK_GhostProb));  
            
            Ds_m12 = (pi_plus_fromDs + pi_minus_fromDs).M(); 
            Ds_m13 = (pi_plus_fromDs + pi_minus2_fromDs).M(); 
            
            beta_K_minus_fromDs = (pi_plus_fromDs_P + pi_minus_fromDs_P - pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
            beta_K_plus_fromDs = (-pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
            beta_pi_minus_fromDs = (pi_plus_fromDs_P - pi_minus_fromDs_P + pi_minus2_fromDs_P) / (pi_plus_fromDs_P + pi_minus_fromDs_P + pi_minus2_fromDs_P);
        }
        
        else if(_Ds_finalState == Ds_finalState::Kpipi){
        
        }
        
        else {
            DsDaughters_min_IPCHI2 = min(K_plus_fromDs_IPCHI2_OWNPV,min(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_max_IPCHI2 = max(K_plus_fromDs_IPCHI2_OWNPV,max(pi_minus_fromDs_IPCHI2_OWNPV,K_minus_fromDs_IPCHI2_OWNPV));            
            DsDaughters_min_PT = min(K_plus_fromDs_PT,min(pi_minus_fromDs_PT,K_minus_fromDs_PT));            
            Ds_max_DOCA = max(Ds_DOCA1,max(Ds_DOCA2,Ds_DOCA3));
            Ds_max_ghostProb = max(K_plus_fromDs_TRACK_GhostProb,max(K_minus_fromDs_TRACK_GhostProb,pi_minus_fromDs_TRACK_GhostProb));  
            
            Ds_m12 = (K_plus_fromDs + K_minus_fromDs).M(); 
            Ds_m13 = (K_plus_fromDs + pi_minus_fromDs).M(); 
            
            bkg_D_as_Ds_m= (K_plus_fromDs + Kminus_fromDs_asPiminus_MissID + pi_minus_fromDs).M() - massDminus;
            bkg_Lambdac_as_Ds_m=  (K_plus_fromDs + Kminus_fromDs_asProton_MissID + pi_minus_fromDs).M() - massLambda_c;  
            
            beta_K_minus_fromDs = (K_plus_fromDs_P + pi_minus_fromDs_P - K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_K_plus_fromDs = (-K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);
            beta_pi_minus_fromDs = (K_plus_fromDs_P - pi_minus_fromDs_P + K_minus_fromDs_P) / (K_plus_fromDs_P + pi_minus_fromDs_P + K_minus_fromDs_P);

        }
        
        max_ghostProb = max(Xs_max_ghostProb,Ds_max_ghostProb);

        Bs_ct = Bs_TAU*1000;
        Bs_cterr = Bs_TAUERR*1000;
        Ds_ct = Ds_TAU*1000;
        Ds_cterr = Ds_TAUERR*1000;

        // ???
        if(Bs_ID > 0) qt = 1;
        if(Bs_ID < 0) qt = -1;
        if(Bs_ID == 0) qt = 0;
        if(Ds_ID > 0) qf = -1;
        if(Ds_ID < 0) qf = +1;

        
        // Fill tree
        summary_tree->Fill();
        if(0ul == (i % 100000ul))summary_tree->AutoSave();
    }
   
    cout << "Selected " << summary_tree->GetEntries() << " events" <<  endl;
    cout << "Efficiency = " << summary_tree->GetEntries()/(double)nentries * 100. << " %" <<  endl;
    
    summary_tree->Write();
    output->Close();
    
    cout << "Created new file " << _outFileName << endl;
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;

}
