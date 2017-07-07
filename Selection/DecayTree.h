//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  4 15:02:48 2017 by ROOT version 5.34/36
// from TTree DecayTree/DecayTree
// found on file: b2hhh_4.root
//////////////////////////////////////////////////////////

#ifndef DecayTree_h
#define DecayTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;

struct Ds_finalState{
    enum  Type { phipi, Kspi, KKpi_NR, pipipi, Kpipi };
};

struct Decay{
    enum  Type { signal, norm };
};

struct Year{
    enum  Type { y11 = 11, y12 = 12, y15 = 15, y16 = 16 };
};

struct DataType{
    enum  Type { mc, data };
};

static const double massKaon = 493.68;
static const double massPion = 139.57;
static const double massProton = 938.27;
static const double massPhi = 1019.46;
static const double massKstar = 895.81;
static const double massDs = 1968.30;
static const double massDminus = 1869.61;
static const double massLambda_c = 2286.46;
static const double c_light = 0.299792458;

class DecayTree {
public :
    TString _outFileName;
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    Decay::Type _decay;
    Year::Type _year;
    Ds_finalState::Type _Ds_finalState;
    DataType::Type _data;

    DecayTree(Decay::Type decay, Year::Type year, Ds_finalState::Type finalState, DataType::Type dataType );
    virtual ~DecayTree();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual TTree* GetInputTree();
    
    virtual  Bool_t TriggerCuts();
    virtual  Bool_t LooseCuts();

    // Declaration of leaf types
    Double_t        Bs_DOCA1;
    Double_t        Bs_DOCA2;
    Double_t        Bs_DOCA3;
    Double_t        Bs_ETA;
    Double_t        Bs_ENDVERTEX_X;
    Double_t        Bs_ENDVERTEX_Y;
    Double_t        Bs_ENDVERTEX_Z;
    Double_t        Bs_ENDVERTEX_XERR;
    Double_t        Bs_ENDVERTEX_YERR;
    Double_t        Bs_ENDVERTEX_ZERR;
    Double_t        Bs_ENDVERTEX_CHI2;
    Int_t           Bs_ENDVERTEX_NDOF;
    Float_t         Bs_ENDVERTEX_COV_[3][3];
    Double_t        Bs_OWNPV_X;
    Double_t        Bs_OWNPV_Y;
    Double_t        Bs_OWNPV_Z;
    Double_t        Bs_OWNPV_XERR;
    Double_t        Bs_OWNPV_YERR;
    Double_t        Bs_OWNPV_ZERR;
    Double_t        Bs_OWNPV_CHI2;
    Int_t           Bs_OWNPV_NDOF;
    Float_t         Bs_OWNPV_COV_[3][3];
    Double_t        Bs_IP_OWNPV;
    Double_t        Bs_IPCHI2_OWNPV;
    Double_t        Bs_FD_OWNPV;
    Double_t        Bs_FDCHI2_OWNPV;
    Double_t        Bs_DIRA_OWNPV;
    Double_t        Bs_P;
    Double_t        Bs_PT;
    Double_t        Bs_PE;
    Double_t        Bs_PX;
    Double_t        Bs_PY;
    Double_t        Bs_PZ;
    Double_t        Bs_MM;
    Double_t        Bs_MMERR;
    Double_t        Bs_M;
    Int_t           Bs_ID;
    Double_t        Bs_TAU;
    Double_t        Bs_TAUERR;
    Double_t        Bs_TAUCHI2;
    Bool_t          Bs_L0Global_Dec;
    Bool_t          Bs_L0Global_TIS;
    Bool_t          Bs_L0Global_TOS;
    Bool_t          Bs_Hlt1Global_Dec;
    Bool_t          Bs_Hlt1Global_TIS;
    Bool_t          Bs_Hlt1Global_TOS;
    Bool_t          Bs_Hlt1Phys_Dec;
    Bool_t          Bs_Hlt1Phys_TIS;
    Bool_t          Bs_Hlt1Phys_TOS;
    Bool_t          Bs_Hlt2Global_Dec;
    Bool_t          Bs_Hlt2Global_TIS;
    Bool_t          Bs_Hlt2Global_TOS;
    Bool_t          Bs_Hlt2Phys_Dec;
    Bool_t          Bs_Hlt2Phys_TIS;
    Bool_t          Bs_Hlt2Phys_TOS;
    Bool_t          Bs_L0HadronDecision_Dec;
    Bool_t          Bs_L0HadronDecision_TIS;
    Bool_t          Bs_L0HadronDecision_TOS;
    Bool_t          Bs_L0MuonDecision_Dec;
    Bool_t          Bs_L0MuonDecision_TIS;
    Bool_t          Bs_L0MuonDecision_TOS;
    Bool_t          Bs_L0GlobalDecision_Dec;
    Bool_t          Bs_L0GlobalDecision_TIS;
    Bool_t          Bs_L0GlobalDecision_TOS;
    Bool_t          Bs_Hlt1TrackAllL0Decision_Dec;
    Bool_t          Bs_Hlt1TrackAllL0Decision_TIS;
    Bool_t          Bs_Hlt1TrackAllL0Decision_TOS;
    Bool_t          Bs_Hlt1TrackMVADecision_Dec;
    Bool_t          Bs_Hlt1TrackMVADecision_TIS;
    Bool_t          Bs_Hlt1TrackMVADecision_TOS;
    Bool_t          Bs_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          Bs_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          Bs_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          Bs_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          Bs_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          Bs_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          Bs_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          Bs_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          Bs_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          Bs_Hlt2IncPhiDecision_Dec;
    Bool_t          Bs_Hlt2IncPhiDecision_TIS;
    Bool_t          Bs_Hlt2IncPhiDecision_TOS;
    Bool_t          Bs_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          Bs_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          Bs_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          Bs_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          Bs_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          Bs_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          Bs_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          Bs_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          Bs_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          Bs_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          Bs_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          Bs_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          Bs_Hlt2Topo2BodyDecision_Dec;
    Bool_t          Bs_Hlt2Topo2BodyDecision_TIS;
    Bool_t          Bs_Hlt2Topo2BodyDecision_TOS;
    Bool_t          Bs_Hlt2Topo3BodyDecision_Dec;
    Bool_t          Bs_Hlt2Topo3BodyDecision_TIS;
    Bool_t          Bs_Hlt2Topo3BodyDecision_TOS;
    Bool_t          Bs_Hlt2Topo4BodyDecision_Dec;
    Bool_t          Bs_Hlt2Topo4BodyDecision_TIS;
    Bool_t          Bs_Hlt2Topo4BodyDecision_TOS;
    Int_t           Bs_TAGDECISION;
    Double_t        Bs_TAGOMEGA;
    Int_t           Bs_TAGDECISION_OS;
    Double_t        Bs_TAGOMEGA_OS;
    Int_t           Bs_TAGGER;
    Short_t         Bs_OS_Muon_DEC;
    Float_t         Bs_OS_Muon_PROB;
    Short_t         Bs_OS_Electron_DEC;
    Float_t         Bs_OS_Electron_PROB;
    Short_t         Bs_OS_Kaon_DEC;
    Float_t         Bs_OS_Kaon_PROB;
    Short_t         Bs_SS_Kaon_DEC;
    Float_t         Bs_SS_Kaon_PROB;
    Short_t         Bs_SS_Pion_DEC;
    Float_t         Bs_SS_Pion_PROB;
    Short_t         Bs_SS_PionBDT_DEC;
    Float_t         Bs_SS_PionBDT_PROB;
    Short_t         Bs_VtxCharge_DEC;
    Float_t         Bs_VtxCharge_PROB;
    Short_t         Bs_OS_nnetKaon_DEC;
    Float_t         Bs_OS_nnetKaon_PROB;
    Short_t         Bs_SS_nnetKaon_DEC;
    Float_t         Bs_SS_nnetKaon_PROB;
    Short_t         Bs_SS_Proton_DEC;
    Float_t         Bs_SS_Proton_PROB;
    Short_t         Bs_OS_Charm_DEC;
    Float_t         Bs_OS_Charm_PROB;
    Double_t        Bs_cpx_0_50;
    Double_t        Bs_cpy_0_50;
    Double_t        Bs_cpz_0_50;
    Double_t        Bs_cpt_0_50;
    Double_t        Bs_cp_0_50;
    Int_t           Bs_cmult_0_50;
    Double_t        Bs_deltaEta_0_50;
    Double_t        Bs_deltaPhi_0_50;
    Double_t        Bs_pxasy_0_50;
    Double_t        Bs_pyasy_0_50;
    Double_t        Bs_pzasy_0_50;
    Double_t        Bs_pasy_0_50;
    Double_t        Bs_ptasy_0_50;
    Double_t        Bs_cpx_0_60;
    Double_t        Bs_cpy_0_60;
    Double_t        Bs_cpz_0_60;
    Double_t        Bs_cpt_0_60;
    Double_t        Bs_cp_0_60;
    Int_t           Bs_cmult_0_60;
    Double_t        Bs_deltaEta_0_60;
    Double_t        Bs_deltaPhi_0_60;
    Double_t        Bs_pxasy_0_60;
    Double_t        Bs_pyasy_0_60;
    Double_t        Bs_pzasy_0_60;
    Double_t        Bs_pasy_0_60;
    Double_t        Bs_ptasy_0_60;
    Double_t        Bs_cpx_0_70;
    Double_t        Bs_cpy_0_70;
    Double_t        Bs_cpz_0_70;
    Double_t        Bs_cpt_0_70;
    Double_t        Bs_cp_0_70;
    Int_t           Bs_cmult_0_70;
    Double_t        Bs_deltaEta_0_70;
    Double_t        Bs_deltaPhi_0_70;
    Double_t        Bs_pxasy_0_70;
    Double_t        Bs_pyasy_0_70;
    Double_t        Bs_pzasy_0_70;
    Double_t        Bs_pasy_0_70;
    Double_t        Bs_ptasy_0_70;
    Double_t        Bs_cpx_0_80;
    Double_t        Bs_cpy_0_80;
    Double_t        Bs_cpz_0_80;
    Double_t        Bs_cpt_0_80;
    Double_t        Bs_cp_0_80;
    Int_t           Bs_cmult_0_80;
    Double_t        Bs_deltaEta_0_80;
    Double_t        Bs_deltaPhi_0_80;
    Double_t        Bs_pxasy_0_80;
    Double_t        Bs_pyasy_0_80;
    Double_t        Bs_pzasy_0_80;
    Double_t        Bs_pasy_0_80;
    Double_t        Bs_ptasy_0_80;
    Double_t        Bs_cpx_0_90;
    Double_t        Bs_cpy_0_90;
    Double_t        Bs_cpz_0_90;
    Double_t        Bs_cpt_0_90;
    Double_t        Bs_cp_0_90;
    Int_t           Bs_cmult_0_90;
    Double_t        Bs_deltaEta_0_90;
    Double_t        Bs_deltaPhi_0_90;
    Double_t        Bs_pxasy_0_90;
    Double_t        Bs_pyasy_0_90;
    Double_t        Bs_pzasy_0_90;
    Double_t        Bs_pasy_0_90;
    Double_t        Bs_ptasy_0_90;
    Double_t        Bs_cpx_1_00;
    Double_t        Bs_cpy_1_00;
    Double_t        Bs_cpz_1_00;
    Double_t        Bs_cpt_1_00;
    Double_t        Bs_cp_1_00;
    Int_t           Bs_cmult_1_00;
    Double_t        Bs_deltaEta_1_00;
    Double_t        Bs_deltaPhi_1_00;
    Double_t        Bs_pxasy_1_00;
    Double_t        Bs_pyasy_1_00;
    Double_t        Bs_pzasy_1_00;
    Double_t        Bs_pasy_1_00;
    Double_t        Bs_ptasy_1_00;
    Int_t           Bs_B0DTF_nPV;
    Float_t         Bs_B0DTF_D_splus_Kplus_0_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_0_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_0_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_0_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_0_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_Kplus_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_M[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_MERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_P[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_PERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_ctau[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_ctauErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_decayLength[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_decayLengthErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_Kplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_Kplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_Kplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_Kplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_Kplus_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_M[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_MERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_P[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_PERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_ctau[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_ctauErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_decayLength[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_decayLengthErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_0_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_0_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_0_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_0_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_0_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_K_1_1270_plus_piplus_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_M[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_MERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_P[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_PERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_PV_X[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_PV_Y[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_PV_Z[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_PV_key[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_chi2[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_ctau[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_ctauErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_decayLength[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_decayLengthErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_nDOF[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_nIter[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_status[100];   //[Bs_B0DTF_nPV]
    Int_t           Bs_BsDTF_nPV;
    Float_t         Bs_BsDTF_D_splus_Kplus_0_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_0_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_0_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_0_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_0_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_Kplus_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_M[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_MERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_P[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_PERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_ctau[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_ctauErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_decayLength[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_decayLengthErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_Kplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_Kplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_Kplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_Kplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_Kplus_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_M[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_MERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_P[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_PERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_ctau[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_ctauErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_decayLength[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_decayLengthErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_0_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_0_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_0_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_0_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_0_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_K_1_1270_plus_piplus_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_M[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_MERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_P[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_PERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_PV_X[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_PV_Y[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_PV_Z[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_PV_key[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_chi2[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_ctau[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_ctauErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_decayLength[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_decayLengthErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_nDOF[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_nIter[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_status[100];   //[Bs_BsDTF_nPV]
    Int_t           Bs_DTF_nPV;
    Float_t         Bs_DTF_D_splus_Kplus_0_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_0_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_0_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_0_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_0_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_Kplus_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_M[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_MERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_P[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_PERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_ctau[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_ctauErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_decayLength[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_decayLengthErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_Kplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_Kplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_Kplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_Kplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_Kplus_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_M[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_MERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_P[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_PERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_ctau[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_ctauErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_decayLength[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_decayLengthErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_0_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_0_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_0_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_0_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_0_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_K_1_1270_plus_piplus_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_M[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_MERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_P[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_PERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_PV_X[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_PV_Y[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_PV_Z[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_PV_key[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_chi2[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_ctau[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_ctauErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_decayLength[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_decayLengthErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_nDOF[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_nIter[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_status[100];   //[Bs_DTF_nPV]
    Int_t           Bs_PV_nPV;
    Float_t         Bs_PV_Dplus_Kplus_0_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_0_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_0_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_0_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_0_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_Kplus_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_M[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_MERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_P[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_PERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_ctau[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_ctauErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_decayLength[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_decayLengthErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_Kplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_Kplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_Kplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_Kplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_Kplus_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_M[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_MERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_P[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_PERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_ctau[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_ctauErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_decayLength[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_decayLengthErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_0_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_0_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_0_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_0_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_0_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_K_1_1270_plus_piplus_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_M[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_MERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_P[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_PERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_PV_X[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_PV_Y[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_PV_Z[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_PV_key[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_chi2[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_ctau[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_ctauErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_decayLength[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_decayLengthErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_nDOF[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_nIter[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_status[100];   //[Bs_PV_nPV]
    Int_t           Bs_BsTaggingTool_TAGDECISION;
    Double_t        Bs_BsTaggingTool_TAGOMEGA;
    Int_t           Bs_BsTaggingTool_TAGDECISION_OS;
    Double_t        Bs_BsTaggingTool_TAGOMEGA_OS;
    Int_t           Bs_BsTaggingTool_TAGGER;
    Short_t         Bs_BsTaggingTool_OS_Muon_DEC;
    Float_t         Bs_BsTaggingTool_OS_Muon_PROB;
    Short_t         Bs_BsTaggingTool_OS_Electron_DEC;
    Float_t         Bs_BsTaggingTool_OS_Electron_PROB;
    Short_t         Bs_BsTaggingTool_OS_Kaon_DEC;
    Float_t         Bs_BsTaggingTool_OS_Kaon_PROB;
    Short_t         Bs_BsTaggingTool_SS_Kaon_DEC;
    Float_t         Bs_BsTaggingTool_SS_Kaon_PROB;
    Short_t         Bs_BsTaggingTool_SS_Pion_DEC;
    Float_t         Bs_BsTaggingTool_SS_Pion_PROB;
    Short_t         Bs_BsTaggingTool_SS_PionBDT_DEC;
    Float_t         Bs_BsTaggingTool_SS_PionBDT_PROB;
    Short_t         Bs_BsTaggingTool_VtxCharge_DEC;
    Float_t         Bs_BsTaggingTool_VtxCharge_PROB;
    Short_t         Bs_BsTaggingTool_OS_nnetKaon_DEC;
    Float_t         Bs_BsTaggingTool_OS_nnetKaon_PROB;
    Short_t         Bs_BsTaggingTool_SS_nnetKaon_DEC;
    Float_t         Bs_BsTaggingTool_SS_nnetKaon_PROB;
    Short_t         Bs_BsTaggingTool_SS_Proton_DEC;
    Float_t         Bs_BsTaggingTool_SS_Proton_PROB;
    Short_t         Bs_BsTaggingTool_OS_Charm_DEC;
    Float_t         Bs_BsTaggingTool_OS_Charm_PROB;
    Double_t        Ds_DOCA1;
    Double_t        Ds_DOCA2;
    Double_t        Ds_DOCA3;
    Double_t        Ds_ETA;
    Double_t        Ds_CosTheta;
    Double_t        Ds_ENDVERTEX_X;
    Double_t        Ds_ENDVERTEX_Y;
    Double_t        Ds_ENDVERTEX_Z;
    Double_t        Ds_ENDVERTEX_XERR;
    Double_t        Ds_ENDVERTEX_YERR;
    Double_t        Ds_ENDVERTEX_ZERR;
    Double_t        Ds_ENDVERTEX_CHI2;
    Int_t           Ds_ENDVERTEX_NDOF;
    Float_t         Ds_ENDVERTEX_COV_[3][3];
    Double_t        Ds_OWNPV_X;
    Double_t        Ds_OWNPV_Y;
    Double_t        Ds_OWNPV_Z;
    Double_t        Ds_OWNPV_XERR;
    Double_t        Ds_OWNPV_YERR;
    Double_t        Ds_OWNPV_ZERR;
    Double_t        Ds_OWNPV_CHI2;
    Int_t           Ds_OWNPV_NDOF;
    Float_t         Ds_OWNPV_COV_[3][3];
    Double_t        Ds_IP_OWNPV;
    Double_t        Ds_IPCHI2_OWNPV;
    Double_t        Ds_FD_OWNPV;
    Double_t        Ds_FDCHI2_OWNPV;
    Double_t        Ds_DIRA_OWNPV;
    Double_t        Ds_ORIVX_X;
    Double_t        Ds_ORIVX_Y;
    Double_t        Ds_ORIVX_Z;
    Double_t        Ds_ORIVX_XERR;
    Double_t        Ds_ORIVX_YERR;
    Double_t        Ds_ORIVX_ZERR;
    Double_t        Ds_ORIVX_CHI2;
    Int_t           Ds_ORIVX_NDOF;
    Float_t         Ds_ORIVX_COV_[3][3];
    Double_t        Ds_FD_ORIVX;
    Double_t        Ds_FDCHI2_ORIVX;
    Double_t        Ds_DIRA_ORIVX;
    Double_t        Ds_P;
    Double_t        Ds_PT;
    Double_t        Ds_PE;
    Double_t        Ds_PX;
    Double_t        Ds_PY;
    Double_t        Ds_PZ;
    Double_t        Ds_MM;
    Double_t        Ds_MMERR;
    Double_t        Ds_M;
    Int_t           Ds_ID;
    Double_t        Ds_TAU;
    Double_t        Ds_TAUERR;
    Double_t        Ds_TAUCHI2;
    Bool_t          Ds_L0Global_Dec;
    Bool_t          Ds_L0Global_TIS;
    Bool_t          Ds_L0Global_TOS;
    Bool_t          Ds_Hlt1Global_Dec;
    Bool_t          Ds_Hlt1Global_TIS;
    Bool_t          Ds_Hlt1Global_TOS;
    Bool_t          Ds_Hlt1Phys_Dec;
    Bool_t          Ds_Hlt1Phys_TIS;
    Bool_t          Ds_Hlt1Phys_TOS;
    Bool_t          Ds_Hlt2Global_Dec;
    Bool_t          Ds_Hlt2Global_TIS;
    Bool_t          Ds_Hlt2Global_TOS;
    Bool_t          Ds_Hlt2Phys_Dec;
    Bool_t          Ds_Hlt2Phys_TIS;
    Bool_t          Ds_Hlt2Phys_TOS;
    Bool_t          Ds_L0HadronDecision_Dec;
    Bool_t          Ds_L0HadronDecision_TIS;
    Bool_t          Ds_L0HadronDecision_TOS;
    Bool_t          Ds_L0MuonDecision_Dec;
    Bool_t          Ds_L0MuonDecision_TIS;
    Bool_t          Ds_L0MuonDecision_TOS;
    Bool_t          Ds_L0GlobalDecision_Dec;
    Bool_t          Ds_L0GlobalDecision_TIS;
    Bool_t          Ds_L0GlobalDecision_TOS;
    Bool_t          Ds_Hlt1TrackAllL0Decision_Dec;
    Bool_t          Ds_Hlt1TrackAllL0Decision_TIS;
    Bool_t          Ds_Hlt1TrackAllL0Decision_TOS;
    Bool_t          Ds_Hlt1TrackMVADecision_Dec;
    Bool_t          Ds_Hlt1TrackMVADecision_TIS;
    Bool_t          Ds_Hlt1TrackMVADecision_TOS;
    Bool_t          Ds_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          Ds_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          Ds_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          Ds_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          Ds_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          Ds_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          Ds_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          Ds_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          Ds_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          Ds_Hlt2IncPhiDecision_Dec;
    Bool_t          Ds_Hlt2IncPhiDecision_TIS;
    Bool_t          Ds_Hlt2IncPhiDecision_TOS;
    Bool_t          Ds_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          Ds_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          Ds_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          Ds_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          Ds_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          Ds_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          Ds_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          Ds_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          Ds_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          Ds_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          Ds_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          Ds_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          Ds_Hlt2Topo2BodyDecision_Dec;
    Bool_t          Ds_Hlt2Topo2BodyDecision_TIS;
    Bool_t          Ds_Hlt2Topo2BodyDecision_TOS;
    Bool_t          Ds_Hlt2Topo3BodyDecision_Dec;
    Bool_t          Ds_Hlt2Topo3BodyDecision_TIS;
    Bool_t          Ds_Hlt2Topo3BodyDecision_TOS;
    Bool_t          Ds_Hlt2Topo4BodyDecision_Dec;
    Bool_t          Ds_Hlt2Topo4BodyDecision_TIS;
    Bool_t          Ds_Hlt2Topo4BodyDecision_TOS;
    Double_t        Ds_cpx_0_50;
    Double_t        Ds_cpy_0_50;
    Double_t        Ds_cpz_0_50;
    Double_t        Ds_cpt_0_50;
    Double_t        Ds_cp_0_50;
    Int_t           Ds_cmult_0_50;
    Double_t        Ds_deltaEta_0_50;
    Double_t        Ds_deltaPhi_0_50;
    Double_t        Ds_pxasy_0_50;
    Double_t        Ds_pyasy_0_50;
    Double_t        Ds_pzasy_0_50;
    Double_t        Ds_pasy_0_50;
    Double_t        Ds_ptasy_0_50;
    Double_t        Ds_cpx_0_60;
    Double_t        Ds_cpy_0_60;
    Double_t        Ds_cpz_0_60;
    Double_t        Ds_cpt_0_60;
    Double_t        Ds_cp_0_60;
    Int_t           Ds_cmult_0_60;
    Double_t        Ds_deltaEta_0_60;
    Double_t        Ds_deltaPhi_0_60;
    Double_t        Ds_pxasy_0_60;
    Double_t        Ds_pyasy_0_60;
    Double_t        Ds_pzasy_0_60;
    Double_t        Ds_pasy_0_60;
    Double_t        Ds_ptasy_0_60;
    Double_t        Ds_cpx_0_70;
    Double_t        Ds_cpy_0_70;
    Double_t        Ds_cpz_0_70;
    Double_t        Ds_cpt_0_70;
    Double_t        Ds_cp_0_70;
    Int_t           Ds_cmult_0_70;
    Double_t        Ds_deltaEta_0_70;
    Double_t        Ds_deltaPhi_0_70;
    Double_t        Ds_pxasy_0_70;
    Double_t        Ds_pyasy_0_70;
    Double_t        Ds_pzasy_0_70;
    Double_t        Ds_pasy_0_70;
    Double_t        Ds_ptasy_0_70;
    Double_t        Ds_cpx_0_80;
    Double_t        Ds_cpy_0_80;
    Double_t        Ds_cpz_0_80;
    Double_t        Ds_cpt_0_80;
    Double_t        Ds_cp_0_80;
    Int_t           Ds_cmult_0_80;
    Double_t        Ds_deltaEta_0_80;
    Double_t        Ds_deltaPhi_0_80;
    Double_t        Ds_pxasy_0_80;
    Double_t        Ds_pyasy_0_80;
    Double_t        Ds_pzasy_0_80;
    Double_t        Ds_pasy_0_80;
    Double_t        Ds_ptasy_0_80;
    Double_t        Ds_cpx_0_90;
    Double_t        Ds_cpy_0_90;
    Double_t        Ds_cpz_0_90;
    Double_t        Ds_cpt_0_90;
    Double_t        Ds_cp_0_90;
    Int_t           Ds_cmult_0_90;
    Double_t        Ds_deltaEta_0_90;
    Double_t        Ds_deltaPhi_0_90;
    Double_t        Ds_pxasy_0_90;
    Double_t        Ds_pyasy_0_90;
    Double_t        Ds_pzasy_0_90;
    Double_t        Ds_pasy_0_90;
    Double_t        Ds_ptasy_0_90;
    Double_t        Ds_cpx_1_00;
    Double_t        Ds_cpy_1_00;
    Double_t        Ds_cpz_1_00;
    Double_t        Ds_cpt_1_00;
    Double_t        Ds_cp_1_00;
    Int_t           Ds_cmult_1_00;
    Double_t        Ds_deltaEta_1_00;
    Double_t        Ds_deltaPhi_1_00;
    Double_t        Ds_pxasy_1_00;
    Double_t        Ds_pyasy_1_00;
    Double_t        Ds_pzasy_1_00;
    Double_t        Ds_pasy_1_00;
    Double_t        Ds_ptasy_1_00;
    Double_t        K_plus_fromDs_DOCA1;
    Double_t        K_plus_fromDs_DOCA2;
    Double_t        K_plus_fromDs_DOCA3;
    Double_t        K_plus_fromDs_ETA;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNe;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNmu;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNpi;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNk;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNp;
    Double_t        K_plus_fromDs_MC12TuneV2_ProbNNghost;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNe;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNmu;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNpi;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNk;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNp;
    Double_t        K_plus_fromDs_MC12TuneV3_ProbNNghost;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNe;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNmu;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNpi;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNk;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNp;
    Double_t        K_plus_fromDs_MC12TuneV4_ProbNNghost;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNe;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNmu;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNpi;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNk;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNp;
    Double_t        K_plus_fromDs_MC15TuneV1_ProbNNghost;
    Double_t        K_plus_fromDs_CosTheta;
    Double_t        K_plus_fromDs_OWNPV_X;
    Double_t        K_plus_fromDs_OWNPV_Y;
    Double_t        K_plus_fromDs_OWNPV_Z;
    Double_t        K_plus_fromDs_OWNPV_XERR;
    Double_t        K_plus_fromDs_OWNPV_YERR;
    Double_t        K_plus_fromDs_OWNPV_ZERR;
    Double_t        K_plus_fromDs_OWNPV_CHI2;
    Int_t           K_plus_fromDs_OWNPV_NDOF;
    Float_t         K_plus_fromDs_OWNPV_COV_[3][3];
    Double_t        K_plus_fromDs_IP_OWNPV;
    Double_t        K_plus_fromDs_IPCHI2_OWNPV;
    Double_t        K_plus_fromDs_ORIVX_X;
    Double_t        K_plus_fromDs_ORIVX_Y;
    Double_t        K_plus_fromDs_ORIVX_Z;
    Double_t        K_plus_fromDs_ORIVX_XERR;
    Double_t        K_plus_fromDs_ORIVX_YERR;
    Double_t        K_plus_fromDs_ORIVX_ZERR;
    Double_t        K_plus_fromDs_ORIVX_CHI2;
    Int_t           K_plus_fromDs_ORIVX_NDOF;
    Float_t         K_plus_fromDs_ORIVX_COV_[3][3];
    Double_t        K_plus_fromDs_P;
    Double_t        K_plus_fromDs_PT;
    Double_t        K_plus_fromDs_PE;
    Double_t        K_plus_fromDs_PX;
    Double_t        K_plus_fromDs_PY;
    Double_t        K_plus_fromDs_PZ;
    Double_t        K_plus_fromDs_M;
    Int_t           K_plus_fromDs_ID;
    Double_t        K_plus_fromDs_PIDe;
    Double_t        K_plus_fromDs_PIDmu;
    Double_t        K_plus_fromDs_PIDK;
    Double_t        K_plus_fromDs_PIDp;
    Double_t        K_plus_fromDs_ProbNNe;
    Double_t        K_plus_fromDs_ProbNNk;
    Double_t        K_plus_fromDs_ProbNNp;
    Double_t        K_plus_fromDs_ProbNNpi;
    Double_t        K_plus_fromDs_ProbNNmu;
    Double_t        K_plus_fromDs_ProbNNghost;
    Bool_t          K_plus_fromDs_hasMuon;
    Bool_t          K_plus_fromDs_isMuon;
    Bool_t          K_plus_fromDs_hasRich;
    Bool_t          K_plus_fromDs_UsedRichAerogel;
    Bool_t          K_plus_fromDs_UsedRich1Gas;
    Bool_t          K_plus_fromDs_UsedRich2Gas;
    Bool_t          K_plus_fromDs_RichAboveElThres;
    Bool_t          K_plus_fromDs_RichAboveMuThres;
    Bool_t          K_plus_fromDs_RichAbovePiThres;
    Bool_t          K_plus_fromDs_RichAboveKaThres;
    Bool_t          K_plus_fromDs_RichAbovePrThres;
    Bool_t          K_plus_fromDs_hasCalo;
    Bool_t          K_plus_fromDs_L0Global_Dec;
    Bool_t          K_plus_fromDs_L0Global_TIS;
    Bool_t          K_plus_fromDs_L0Global_TOS;
    Bool_t          K_plus_fromDs_Hlt1Global_Dec;
    Bool_t          K_plus_fromDs_Hlt1Global_TIS;
    Bool_t          K_plus_fromDs_Hlt1Global_TOS;
    Bool_t          K_plus_fromDs_Hlt1Phys_Dec;
    Bool_t          K_plus_fromDs_Hlt1Phys_TIS;
    Bool_t          K_plus_fromDs_Hlt1Phys_TOS;
    Bool_t          K_plus_fromDs_Hlt2Global_Dec;
    Bool_t          K_plus_fromDs_Hlt2Global_TIS;
    Bool_t          K_plus_fromDs_Hlt2Global_TOS;
    Bool_t          K_plus_fromDs_Hlt2Phys_Dec;
    Bool_t          K_plus_fromDs_Hlt2Phys_TIS;
    Bool_t          K_plus_fromDs_Hlt2Phys_TOS;
    Bool_t          K_plus_fromDs_L0HadronDecision_Dec;
    Bool_t          K_plus_fromDs_L0HadronDecision_TIS;
    Bool_t          K_plus_fromDs_L0HadronDecision_TOS;
    Bool_t          K_plus_fromDs_L0MuonDecision_Dec;
    Bool_t          K_plus_fromDs_L0MuonDecision_TIS;
    Bool_t          K_plus_fromDs_L0MuonDecision_TOS;
    Bool_t          K_plus_fromDs_L0GlobalDecision_Dec;
    Bool_t          K_plus_fromDs_L0GlobalDecision_TIS;
    Bool_t          K_plus_fromDs_L0GlobalDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt1TrackAllL0Decision_Dec;
    Bool_t          K_plus_fromDs_Hlt1TrackAllL0Decision_TIS;
    Bool_t          K_plus_fromDs_Hlt1TrackAllL0Decision_TOS;
    Bool_t          K_plus_fromDs_Hlt1TrackMVADecision_Dec;
    Bool_t          K_plus_fromDs_Hlt1TrackMVADecision_TIS;
    Bool_t          K_plus_fromDs_Hlt1TrackMVADecision_TOS;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          K_plus_fromDs_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2IncPhiDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2IncPhiDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2IncPhiDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo2BodyDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo3BodyDecision_TOS;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyDecision_Dec;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyDecision_TIS;
    Bool_t          K_plus_fromDs_Hlt2Topo4BodyDecision_TOS;
    Int_t           K_plus_fromDs_TRACK_Type;
    Int_t           K_plus_fromDs_TRACK_Key;
    Double_t        K_plus_fromDs_TRACK_CHI2NDOF;
    Double_t        K_plus_fromDs_TRACK_PCHI2;
    Double_t        K_plus_fromDs_TRACK_MatchCHI2;
    Double_t        K_plus_fromDs_TRACK_GhostProb;
    Double_t        K_plus_fromDs_TRACK_CloneDist;
    Double_t        K_plus_fromDs_TRACK_Likelihood;
    Double_t        K_plus_fromDs_cpx_0_50;
    Double_t        K_plus_fromDs_cpy_0_50;
    Double_t        K_plus_fromDs_cpz_0_50;
    Double_t        K_plus_fromDs_cpt_0_50;
    Double_t        K_plus_fromDs_cp_0_50;
    Int_t           K_plus_fromDs_cmult_0_50;
    Double_t        K_plus_fromDs_deltaEta_0_50;
    Double_t        K_plus_fromDs_deltaPhi_0_50;
    Double_t        K_plus_fromDs_pxasy_0_50;
    Double_t        K_plus_fromDs_pyasy_0_50;
    Double_t        K_plus_fromDs_pzasy_0_50;
    Double_t        K_plus_fromDs_pasy_0_50;
    Double_t        K_plus_fromDs_ptasy_0_50;
    Double_t        K_plus_fromDs_cpx_0_60;
    Double_t        K_plus_fromDs_cpy_0_60;
    Double_t        K_plus_fromDs_cpz_0_60;
    Double_t        K_plus_fromDs_cpt_0_60;
    Double_t        K_plus_fromDs_cp_0_60;
    Int_t           K_plus_fromDs_cmult_0_60;
    Double_t        K_plus_fromDs_deltaEta_0_60;
    Double_t        K_plus_fromDs_deltaPhi_0_60;
    Double_t        K_plus_fromDs_pxasy_0_60;
    Double_t        K_plus_fromDs_pyasy_0_60;
    Double_t        K_plus_fromDs_pzasy_0_60;
    Double_t        K_plus_fromDs_pasy_0_60;
    Double_t        K_plus_fromDs_ptasy_0_60;
    Double_t        K_plus_fromDs_cpx_0_70;
    Double_t        K_plus_fromDs_cpy_0_70;
    Double_t        K_plus_fromDs_cpz_0_70;
    Double_t        K_plus_fromDs_cpt_0_70;
    Double_t        K_plus_fromDs_cp_0_70;
    Int_t           K_plus_fromDs_cmult_0_70;
    Double_t        K_plus_fromDs_deltaEta_0_70;
    Double_t        K_plus_fromDs_deltaPhi_0_70;
    Double_t        K_plus_fromDs_pxasy_0_70;
    Double_t        K_plus_fromDs_pyasy_0_70;
    Double_t        K_plus_fromDs_pzasy_0_70;
    Double_t        K_plus_fromDs_pasy_0_70;
    Double_t        K_plus_fromDs_ptasy_0_70;
    Double_t        K_plus_fromDs_cpx_0_80;
    Double_t        K_plus_fromDs_cpy_0_80;
    Double_t        K_plus_fromDs_cpz_0_80;
    Double_t        K_plus_fromDs_cpt_0_80;
    Double_t        K_plus_fromDs_cp_0_80;
    Int_t           K_plus_fromDs_cmult_0_80;
    Double_t        K_plus_fromDs_deltaEta_0_80;
    Double_t        K_plus_fromDs_deltaPhi_0_80;
    Double_t        K_plus_fromDs_pxasy_0_80;
    Double_t        K_plus_fromDs_pyasy_0_80;
    Double_t        K_plus_fromDs_pzasy_0_80;
    Double_t        K_plus_fromDs_pasy_0_80;
    Double_t        K_plus_fromDs_ptasy_0_80;
    Double_t        K_plus_fromDs_cpx_0_90;
    Double_t        K_plus_fromDs_cpy_0_90;
    Double_t        K_plus_fromDs_cpz_0_90;
    Double_t        K_plus_fromDs_cpt_0_90;
    Double_t        K_plus_fromDs_cp_0_90;
    Int_t           K_plus_fromDs_cmult_0_90;
    Double_t        K_plus_fromDs_deltaEta_0_90;
    Double_t        K_plus_fromDs_deltaPhi_0_90;
    Double_t        K_plus_fromDs_pxasy_0_90;
    Double_t        K_plus_fromDs_pyasy_0_90;
    Double_t        K_plus_fromDs_pzasy_0_90;
    Double_t        K_plus_fromDs_pasy_0_90;
    Double_t        K_plus_fromDs_ptasy_0_90;
    Double_t        K_plus_fromDs_cpx_1_00;
    Double_t        K_plus_fromDs_cpy_1_00;
    Double_t        K_plus_fromDs_cpz_1_00;
    Double_t        K_plus_fromDs_cpt_1_00;
    Double_t        K_plus_fromDs_cp_1_00;
    Int_t           K_plus_fromDs_cmult_1_00;
    Double_t        K_plus_fromDs_deltaEta_1_00;
    Double_t        K_plus_fromDs_deltaPhi_1_00;
    Double_t        K_plus_fromDs_pxasy_1_00;
    Double_t        K_plus_fromDs_pyasy_1_00;
    Double_t        K_plus_fromDs_pzasy_1_00;
    Double_t        K_plus_fromDs_pasy_1_00;
    Double_t        K_plus_fromDs_ptasy_1_00;
    Double_t        K_minus_fromDs_DOCA1;
    Double_t        K_minus_fromDs_DOCA2;
    Double_t        K_minus_fromDs_DOCA3;
    Double_t        K_minus_fromDs_ETA;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNe;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNmu;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNpi;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNk;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNp;
    Double_t        K_minus_fromDs_MC12TuneV2_ProbNNghost;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNe;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNmu;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNpi;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNk;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNp;
    Double_t        K_minus_fromDs_MC12TuneV3_ProbNNghost;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNe;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNmu;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNpi;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNk;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNp;
    Double_t        K_minus_fromDs_MC12TuneV4_ProbNNghost;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNe;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNmu;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNpi;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNk;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNp;
    Double_t        K_minus_fromDs_MC15TuneV1_ProbNNghost;
    Double_t        K_minus_fromDs_CosTheta;
    Double_t        K_minus_fromDs_OWNPV_X;
    Double_t        K_minus_fromDs_OWNPV_Y;
    Double_t        K_minus_fromDs_OWNPV_Z;
    Double_t        K_minus_fromDs_OWNPV_XERR;
    Double_t        K_minus_fromDs_OWNPV_YERR;
    Double_t        K_minus_fromDs_OWNPV_ZERR;
    Double_t        K_minus_fromDs_OWNPV_CHI2;
    Int_t           K_minus_fromDs_OWNPV_NDOF;
    Float_t         K_minus_fromDs_OWNPV_COV_[3][3];
    Double_t        K_minus_fromDs_IP_OWNPV;
    Double_t        K_minus_fromDs_IPCHI2_OWNPV;
    Double_t        K_minus_fromDs_ORIVX_X;
    Double_t        K_minus_fromDs_ORIVX_Y;
    Double_t        K_minus_fromDs_ORIVX_Z;
    Double_t        K_minus_fromDs_ORIVX_XERR;
    Double_t        K_minus_fromDs_ORIVX_YERR;
    Double_t        K_minus_fromDs_ORIVX_ZERR;
    Double_t        K_minus_fromDs_ORIVX_CHI2;
    Int_t           K_minus_fromDs_ORIVX_NDOF;
    Float_t         K_minus_fromDs_ORIVX_COV_[3][3];
    Double_t        K_minus_fromDs_P;
    Double_t        K_minus_fromDs_PT;
    Double_t        K_minus_fromDs_PE;
    Double_t        K_minus_fromDs_PX;
    Double_t        K_minus_fromDs_PY;
    Double_t        K_minus_fromDs_PZ;
    Double_t        K_minus_fromDs_M;
    Int_t           K_minus_fromDs_ID;
    Double_t        K_minus_fromDs_PIDe;
    Double_t        K_minus_fromDs_PIDmu;
    Double_t        K_minus_fromDs_PIDK;
    Double_t        K_minus_fromDs_PIDp;
    Double_t        K_minus_fromDs_ProbNNe;
    Double_t        K_minus_fromDs_ProbNNk;
    Double_t        K_minus_fromDs_ProbNNp;
    Double_t        K_minus_fromDs_ProbNNpi;
    Double_t        K_minus_fromDs_ProbNNmu;
    Double_t        K_minus_fromDs_ProbNNghost;
    Bool_t          K_minus_fromDs_hasMuon;
    Bool_t          K_minus_fromDs_isMuon;
    Bool_t          K_minus_fromDs_hasRich;
    Bool_t          K_minus_fromDs_UsedRichAerogel;
    Bool_t          K_minus_fromDs_UsedRich1Gas;
    Bool_t          K_minus_fromDs_UsedRich2Gas;
    Bool_t          K_minus_fromDs_RichAboveElThres;
    Bool_t          K_minus_fromDs_RichAboveMuThres;
    Bool_t          K_minus_fromDs_RichAbovePiThres;
    Bool_t          K_minus_fromDs_RichAboveKaThres;
    Bool_t          K_minus_fromDs_RichAbovePrThres;
    Bool_t          K_minus_fromDs_hasCalo;
    Bool_t          K_minus_fromDs_L0Global_Dec;
    Bool_t          K_minus_fromDs_L0Global_TIS;
    Bool_t          K_minus_fromDs_L0Global_TOS;
    Bool_t          K_minus_fromDs_Hlt1Global_Dec;
    Bool_t          K_minus_fromDs_Hlt1Global_TIS;
    Bool_t          K_minus_fromDs_Hlt1Global_TOS;
    Bool_t          K_minus_fromDs_Hlt1Phys_Dec;
    Bool_t          K_minus_fromDs_Hlt1Phys_TIS;
    Bool_t          K_minus_fromDs_Hlt1Phys_TOS;
    Bool_t          K_minus_fromDs_Hlt2Global_Dec;
    Bool_t          K_minus_fromDs_Hlt2Global_TIS;
    Bool_t          K_minus_fromDs_Hlt2Global_TOS;
    Bool_t          K_minus_fromDs_Hlt2Phys_Dec;
    Bool_t          K_minus_fromDs_Hlt2Phys_TIS;
    Bool_t          K_minus_fromDs_Hlt2Phys_TOS;
    Bool_t          K_minus_fromDs_L0HadronDecision_Dec;
    Bool_t          K_minus_fromDs_L0HadronDecision_TIS;
    Bool_t          K_minus_fromDs_L0HadronDecision_TOS;
    Bool_t          K_minus_fromDs_L0MuonDecision_Dec;
    Bool_t          K_minus_fromDs_L0MuonDecision_TIS;
    Bool_t          K_minus_fromDs_L0MuonDecision_TOS;
    Bool_t          K_minus_fromDs_L0GlobalDecision_Dec;
    Bool_t          K_minus_fromDs_L0GlobalDecision_TIS;
    Bool_t          K_minus_fromDs_L0GlobalDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt1TrackAllL0Decision_Dec;
    Bool_t          K_minus_fromDs_Hlt1TrackAllL0Decision_TIS;
    Bool_t          K_minus_fromDs_Hlt1TrackAllL0Decision_TOS;
    Bool_t          K_minus_fromDs_Hlt1TrackMVADecision_Dec;
    Bool_t          K_minus_fromDs_Hlt1TrackMVADecision_TIS;
    Bool_t          K_minus_fromDs_Hlt1TrackMVADecision_TOS;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          K_minus_fromDs_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2IncPhiDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2IncPhiDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2IncPhiDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo2BodyDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo3BodyDecision_TOS;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyDecision_Dec;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyDecision_TIS;
    Bool_t          K_minus_fromDs_Hlt2Topo4BodyDecision_TOS;
    Int_t           K_minus_fromDs_TRACK_Type;
    Int_t           K_minus_fromDs_TRACK_Key;
    Double_t        K_minus_fromDs_TRACK_CHI2NDOF;
    Double_t        K_minus_fromDs_TRACK_PCHI2;
    Double_t        K_minus_fromDs_TRACK_MatchCHI2;
    Double_t        K_minus_fromDs_TRACK_GhostProb;
    Double_t        K_minus_fromDs_TRACK_CloneDist;
    Double_t        K_minus_fromDs_TRACK_Likelihood;
    Double_t        K_minus_fromDs_cpx_0_50;
    Double_t        K_minus_fromDs_cpy_0_50;
    Double_t        K_minus_fromDs_cpz_0_50;
    Double_t        K_minus_fromDs_cpt_0_50;
    Double_t        K_minus_fromDs_cp_0_50;
    Int_t           K_minus_fromDs_cmult_0_50;
    Double_t        K_minus_fromDs_deltaEta_0_50;
    Double_t        K_minus_fromDs_deltaPhi_0_50;
    Double_t        K_minus_fromDs_pxasy_0_50;
    Double_t        K_minus_fromDs_pyasy_0_50;
    Double_t        K_minus_fromDs_pzasy_0_50;
    Double_t        K_minus_fromDs_pasy_0_50;
    Double_t        K_minus_fromDs_ptasy_0_50;
    Double_t        K_minus_fromDs_cpx_0_60;
    Double_t        K_minus_fromDs_cpy_0_60;
    Double_t        K_minus_fromDs_cpz_0_60;
    Double_t        K_minus_fromDs_cpt_0_60;
    Double_t        K_minus_fromDs_cp_0_60;
    Int_t           K_minus_fromDs_cmult_0_60;
    Double_t        K_minus_fromDs_deltaEta_0_60;
    Double_t        K_minus_fromDs_deltaPhi_0_60;
    Double_t        K_minus_fromDs_pxasy_0_60;
    Double_t        K_minus_fromDs_pyasy_0_60;
    Double_t        K_minus_fromDs_pzasy_0_60;
    Double_t        K_minus_fromDs_pasy_0_60;
    Double_t        K_minus_fromDs_ptasy_0_60;
    Double_t        K_minus_fromDs_cpx_0_70;
    Double_t        K_minus_fromDs_cpy_0_70;
    Double_t        K_minus_fromDs_cpz_0_70;
    Double_t        K_minus_fromDs_cpt_0_70;
    Double_t        K_minus_fromDs_cp_0_70;
    Int_t           K_minus_fromDs_cmult_0_70;
    Double_t        K_minus_fromDs_deltaEta_0_70;
    Double_t        K_minus_fromDs_deltaPhi_0_70;
    Double_t        K_minus_fromDs_pxasy_0_70;
    Double_t        K_minus_fromDs_pyasy_0_70;
    Double_t        K_minus_fromDs_pzasy_0_70;
    Double_t        K_minus_fromDs_pasy_0_70;
    Double_t        K_minus_fromDs_ptasy_0_70;
    Double_t        K_minus_fromDs_cpx_0_80;
    Double_t        K_minus_fromDs_cpy_0_80;
    Double_t        K_minus_fromDs_cpz_0_80;
    Double_t        K_minus_fromDs_cpt_0_80;
    Double_t        K_minus_fromDs_cp_0_80;
    Int_t           K_minus_fromDs_cmult_0_80;
    Double_t        K_minus_fromDs_deltaEta_0_80;
    Double_t        K_minus_fromDs_deltaPhi_0_80;
    Double_t        K_minus_fromDs_pxasy_0_80;
    Double_t        K_minus_fromDs_pyasy_0_80;
    Double_t        K_minus_fromDs_pzasy_0_80;
    Double_t        K_minus_fromDs_pasy_0_80;
    Double_t        K_minus_fromDs_ptasy_0_80;
    Double_t        K_minus_fromDs_cpx_0_90;
    Double_t        K_minus_fromDs_cpy_0_90;
    Double_t        K_minus_fromDs_cpz_0_90;
    Double_t        K_minus_fromDs_cpt_0_90;
    Double_t        K_minus_fromDs_cp_0_90;
    Int_t           K_minus_fromDs_cmult_0_90;
    Double_t        K_minus_fromDs_deltaEta_0_90;
    Double_t        K_minus_fromDs_deltaPhi_0_90;
    Double_t        K_minus_fromDs_pxasy_0_90;
    Double_t        K_minus_fromDs_pyasy_0_90;
    Double_t        K_minus_fromDs_pzasy_0_90;
    Double_t        K_minus_fromDs_pasy_0_90;
    Double_t        K_minus_fromDs_ptasy_0_90;
    Double_t        K_minus_fromDs_cpx_1_00;
    Double_t        K_minus_fromDs_cpy_1_00;
    Double_t        K_minus_fromDs_cpz_1_00;
    Double_t        K_minus_fromDs_cpt_1_00;
    Double_t        K_minus_fromDs_cp_1_00;
    Int_t           K_minus_fromDs_cmult_1_00;
    Double_t        K_minus_fromDs_deltaEta_1_00;
    Double_t        K_minus_fromDs_deltaPhi_1_00;
    Double_t        K_minus_fromDs_pxasy_1_00;
    Double_t        K_minus_fromDs_pyasy_1_00;
    Double_t        K_minus_fromDs_pzasy_1_00;
    Double_t        K_minus_fromDs_pasy_1_00;
    Double_t        K_minus_fromDs_ptasy_1_00;
    Double_t        pi_minus_fromDs_DOCA1;
    Double_t        pi_minus_fromDs_DOCA2;
    Double_t        pi_minus_fromDs_DOCA3;
    Double_t        pi_minus_fromDs_ETA;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNe;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNmu;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNpi;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNk;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNp;
    Double_t        pi_minus_fromDs_MC12TuneV2_ProbNNghost;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNe;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNmu;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNpi;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNk;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNp;
    Double_t        pi_minus_fromDs_MC12TuneV3_ProbNNghost;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNe;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNmu;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNpi;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNk;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNp;
    Double_t        pi_minus_fromDs_MC12TuneV4_ProbNNghost;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNe;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNmu;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNpi;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNk;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNp;
    Double_t        pi_minus_fromDs_MC15TuneV1_ProbNNghost;
    Double_t        pi_minus_fromDs_CosTheta;
    Double_t        pi_minus_fromDs_OWNPV_X;
    Double_t        pi_minus_fromDs_OWNPV_Y;
    Double_t        pi_minus_fromDs_OWNPV_Z;
    Double_t        pi_minus_fromDs_OWNPV_XERR;
    Double_t        pi_minus_fromDs_OWNPV_YERR;
    Double_t        pi_minus_fromDs_OWNPV_ZERR;
    Double_t        pi_minus_fromDs_OWNPV_CHI2;
    Int_t           pi_minus_fromDs_OWNPV_NDOF;
    Float_t         pi_minus_fromDs_OWNPV_COV_[3][3];
    Double_t        pi_minus_fromDs_IP_OWNPV;
    Double_t        pi_minus_fromDs_IPCHI2_OWNPV;
    Double_t        pi_minus_fromDs_ORIVX_X;
    Double_t        pi_minus_fromDs_ORIVX_Y;
    Double_t        pi_minus_fromDs_ORIVX_Z;
    Double_t        pi_minus_fromDs_ORIVX_XERR;
    Double_t        pi_minus_fromDs_ORIVX_YERR;
    Double_t        pi_minus_fromDs_ORIVX_ZERR;
    Double_t        pi_minus_fromDs_ORIVX_CHI2;
    Int_t           pi_minus_fromDs_ORIVX_NDOF;
    Float_t         pi_minus_fromDs_ORIVX_COV_[3][3];
    Double_t        pi_minus_fromDs_P;
    Double_t        pi_minus_fromDs_PT;
    Double_t        pi_minus_fromDs_PE;
    Double_t        pi_minus_fromDs_PX;
    Double_t        pi_minus_fromDs_PY;
    Double_t        pi_minus_fromDs_PZ;
    Double_t        pi_minus_fromDs_M;
    Int_t           pi_minus_fromDs_ID;
    Double_t        pi_minus_fromDs_PIDe;
    Double_t        pi_minus_fromDs_PIDmu;
    Double_t        pi_minus_fromDs_PIDK;
    Double_t        pi_minus_fromDs_PIDp;
    Double_t        pi_minus_fromDs_ProbNNe;
    Double_t        pi_minus_fromDs_ProbNNk;
    Double_t        pi_minus_fromDs_ProbNNp;
    Double_t        pi_minus_fromDs_ProbNNpi;
    Double_t        pi_minus_fromDs_ProbNNmu;
    Double_t        pi_minus_fromDs_ProbNNghost;
    Bool_t          pi_minus_fromDs_hasMuon;
    Bool_t          pi_minus_fromDs_isMuon;
    Bool_t          pi_minus_fromDs_hasRich;
    Bool_t          pi_minus_fromDs_UsedRichAerogel;
    Bool_t          pi_minus_fromDs_UsedRich1Gas;
    Bool_t          pi_minus_fromDs_UsedRich2Gas;
    Bool_t          pi_minus_fromDs_RichAboveElThres;
    Bool_t          pi_minus_fromDs_RichAboveMuThres;
    Bool_t          pi_minus_fromDs_RichAbovePiThres;
    Bool_t          pi_minus_fromDs_RichAboveKaThres;
    Bool_t          pi_minus_fromDs_RichAbovePrThres;
    Bool_t          pi_minus_fromDs_hasCalo;
    Bool_t          pi_minus_fromDs_L0Global_Dec;
    Bool_t          pi_minus_fromDs_L0Global_TIS;
    Bool_t          pi_minus_fromDs_L0Global_TOS;
    Bool_t          pi_minus_fromDs_Hlt1Global_Dec;
    Bool_t          pi_minus_fromDs_Hlt1Global_TIS;
    Bool_t          pi_minus_fromDs_Hlt1Global_TOS;
    Bool_t          pi_minus_fromDs_Hlt1Phys_Dec;
    Bool_t          pi_minus_fromDs_Hlt1Phys_TIS;
    Bool_t          pi_minus_fromDs_Hlt1Phys_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Global_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Global_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Global_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Phys_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Phys_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Phys_TOS;
    Bool_t          pi_minus_fromDs_L0HadronDecision_Dec;
    Bool_t          pi_minus_fromDs_L0HadronDecision_TIS;
    Bool_t          pi_minus_fromDs_L0HadronDecision_TOS;
    Bool_t          pi_minus_fromDs_L0MuonDecision_Dec;
    Bool_t          pi_minus_fromDs_L0MuonDecision_TIS;
    Bool_t          pi_minus_fromDs_L0MuonDecision_TOS;
    Bool_t          pi_minus_fromDs_L0GlobalDecision_Dec;
    Bool_t          pi_minus_fromDs_L0GlobalDecision_TIS;
    Bool_t          pi_minus_fromDs_L0GlobalDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt1TrackAllL0Decision_Dec;
    Bool_t          pi_minus_fromDs_Hlt1TrackAllL0Decision_TIS;
    Bool_t          pi_minus_fromDs_Hlt1TrackAllL0Decision_TOS;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVADecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVADecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVADecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2IncPhiDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2IncPhiDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2IncPhiDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo2BodyDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo3BodyDecision_TOS;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyDecision_Dec;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyDecision_TIS;
    Bool_t          pi_minus_fromDs_Hlt2Topo4BodyDecision_TOS;
    Int_t           pi_minus_fromDs_TRACK_Type;
    Int_t           pi_minus_fromDs_TRACK_Key;
    Double_t        pi_minus_fromDs_TRACK_CHI2NDOF;
    Double_t        pi_minus_fromDs_TRACK_PCHI2;
    Double_t        pi_minus_fromDs_TRACK_MatchCHI2;
    Double_t        pi_minus_fromDs_TRACK_GhostProb;
    Double_t        pi_minus_fromDs_TRACK_CloneDist;
    Double_t        pi_minus_fromDs_TRACK_Likelihood;
    Double_t        pi_minus_fromDs_cpx_0_50;
    Double_t        pi_minus_fromDs_cpy_0_50;
    Double_t        pi_minus_fromDs_cpz_0_50;
    Double_t        pi_minus_fromDs_cpt_0_50;
    Double_t        pi_minus_fromDs_cp_0_50;
    Int_t           pi_minus_fromDs_cmult_0_50;
    Double_t        pi_minus_fromDs_deltaEta_0_50;
    Double_t        pi_minus_fromDs_deltaPhi_0_50;
    Double_t        pi_minus_fromDs_pxasy_0_50;
    Double_t        pi_minus_fromDs_pyasy_0_50;
    Double_t        pi_minus_fromDs_pzasy_0_50;
    Double_t        pi_minus_fromDs_pasy_0_50;
    Double_t        pi_minus_fromDs_ptasy_0_50;
    Double_t        pi_minus_fromDs_cpx_0_60;
    Double_t        pi_minus_fromDs_cpy_0_60;
    Double_t        pi_minus_fromDs_cpz_0_60;
    Double_t        pi_minus_fromDs_cpt_0_60;
    Double_t        pi_minus_fromDs_cp_0_60;
    Int_t           pi_minus_fromDs_cmult_0_60;
    Double_t        pi_minus_fromDs_deltaEta_0_60;
    Double_t        pi_minus_fromDs_deltaPhi_0_60;
    Double_t        pi_minus_fromDs_pxasy_0_60;
    Double_t        pi_minus_fromDs_pyasy_0_60;
    Double_t        pi_minus_fromDs_pzasy_0_60;
    Double_t        pi_minus_fromDs_pasy_0_60;
    Double_t        pi_minus_fromDs_ptasy_0_60;
    Double_t        pi_minus_fromDs_cpx_0_70;
    Double_t        pi_minus_fromDs_cpy_0_70;
    Double_t        pi_minus_fromDs_cpz_0_70;
    Double_t        pi_minus_fromDs_cpt_0_70;
    Double_t        pi_minus_fromDs_cp_0_70;
    Int_t           pi_minus_fromDs_cmult_0_70;
    Double_t        pi_minus_fromDs_deltaEta_0_70;
    Double_t        pi_minus_fromDs_deltaPhi_0_70;
    Double_t        pi_minus_fromDs_pxasy_0_70;
    Double_t        pi_minus_fromDs_pyasy_0_70;
    Double_t        pi_minus_fromDs_pzasy_0_70;
    Double_t        pi_minus_fromDs_pasy_0_70;
    Double_t        pi_minus_fromDs_ptasy_0_70;
    Double_t        pi_minus_fromDs_cpx_0_80;
    Double_t        pi_minus_fromDs_cpy_0_80;
    Double_t        pi_minus_fromDs_cpz_0_80;
    Double_t        pi_minus_fromDs_cpt_0_80;
    Double_t        pi_minus_fromDs_cp_0_80;
    Int_t           pi_minus_fromDs_cmult_0_80;
    Double_t        pi_minus_fromDs_deltaEta_0_80;
    Double_t        pi_minus_fromDs_deltaPhi_0_80;
    Double_t        pi_minus_fromDs_pxasy_0_80;
    Double_t        pi_minus_fromDs_pyasy_0_80;
    Double_t        pi_minus_fromDs_pzasy_0_80;
    Double_t        pi_minus_fromDs_pasy_0_80;
    Double_t        pi_minus_fromDs_ptasy_0_80;
    Double_t        pi_minus_fromDs_cpx_0_90;
    Double_t        pi_minus_fromDs_cpy_0_90;
    Double_t        pi_minus_fromDs_cpz_0_90;
    Double_t        pi_minus_fromDs_cpt_0_90;
    Double_t        pi_minus_fromDs_cp_0_90;
    Int_t           pi_minus_fromDs_cmult_0_90;
    Double_t        pi_minus_fromDs_deltaEta_0_90;
    Double_t        pi_minus_fromDs_deltaPhi_0_90;
    Double_t        pi_minus_fromDs_pxasy_0_90;
    Double_t        pi_minus_fromDs_pyasy_0_90;
    Double_t        pi_minus_fromDs_pzasy_0_90;
    Double_t        pi_minus_fromDs_pasy_0_90;
    Double_t        pi_minus_fromDs_ptasy_0_90;
    Double_t        pi_minus_fromDs_cpx_1_00;
    Double_t        pi_minus_fromDs_cpy_1_00;
    Double_t        pi_minus_fromDs_cpz_1_00;
    Double_t        pi_minus_fromDs_cpt_1_00;
    Double_t        pi_minus_fromDs_cp_1_00;
    Int_t           pi_minus_fromDs_cmult_1_00;
    Double_t        pi_minus_fromDs_deltaEta_1_00;
    Double_t        pi_minus_fromDs_deltaPhi_1_00;
    Double_t        pi_minus_fromDs_pxasy_1_00;
    Double_t        pi_minus_fromDs_pyasy_1_00;
    Double_t        pi_minus_fromDs_pzasy_1_00;
    Double_t        pi_minus_fromDs_pasy_1_00;
    Double_t        pi_minus_fromDs_ptasy_1_00;
    Double_t        K_1_1270_plus_DOCA1;
    Double_t        K_1_1270_plus_DOCA2;
    Double_t        K_1_1270_plus_DOCA3;
    Double_t        K_1_1270_plus_ETA;
    Double_t        K_1_1270_plus_CosTheta;
    Double_t        K_1_1270_plus_ENDVERTEX_X;
    Double_t        K_1_1270_plus_ENDVERTEX_Y;
    Double_t        K_1_1270_plus_ENDVERTEX_Z;
    Double_t        K_1_1270_plus_ENDVERTEX_XERR;
    Double_t        K_1_1270_plus_ENDVERTEX_YERR;
    Double_t        K_1_1270_plus_ENDVERTEX_ZERR;
    Double_t        K_1_1270_plus_ENDVERTEX_CHI2;
    Int_t           K_1_1270_plus_ENDVERTEX_NDOF;
    Float_t         K_1_1270_plus_ENDVERTEX_COV_[3][3];
    Double_t        K_1_1270_plus_OWNPV_X;
    Double_t        K_1_1270_plus_OWNPV_Y;
    Double_t        K_1_1270_plus_OWNPV_Z;
    Double_t        K_1_1270_plus_OWNPV_XERR;
    Double_t        K_1_1270_plus_OWNPV_YERR;
    Double_t        K_1_1270_plus_OWNPV_ZERR;
    Double_t        K_1_1270_plus_OWNPV_CHI2;
    Int_t           K_1_1270_plus_OWNPV_NDOF;
    Float_t         K_1_1270_plus_OWNPV_COV_[3][3];
    Double_t        K_1_1270_plus_IP_OWNPV;
    Double_t        K_1_1270_plus_IPCHI2_OWNPV;
    Double_t        K_1_1270_plus_FD_OWNPV;
    Double_t        K_1_1270_plus_FDCHI2_OWNPV;
    Double_t        K_1_1270_plus_DIRA_OWNPV;
    Double_t        K_1_1270_plus_ORIVX_X;
    Double_t        K_1_1270_plus_ORIVX_Y;
    Double_t        K_1_1270_plus_ORIVX_Z;
    Double_t        K_1_1270_plus_ORIVX_XERR;
    Double_t        K_1_1270_plus_ORIVX_YERR;
    Double_t        K_1_1270_plus_ORIVX_ZERR;
    Double_t        K_1_1270_plus_ORIVX_CHI2;
    Int_t           K_1_1270_plus_ORIVX_NDOF;
    Float_t         K_1_1270_plus_ORIVX_COV_[3][3];
    Double_t        K_1_1270_plus_FD_ORIVX;
    Double_t        K_1_1270_plus_FDCHI2_ORIVX;
    Double_t        K_1_1270_plus_DIRA_ORIVX;
    Double_t        K_1_1270_plus_P;
    Double_t        K_1_1270_plus_PT;
    Double_t        K_1_1270_plus_PE;
    Double_t        K_1_1270_plus_PX;
    Double_t        K_1_1270_plus_PY;
    Double_t        K_1_1270_plus_PZ;
    Double_t        K_1_1270_plus_MM;
    Double_t        K_1_1270_plus_MMERR;
    Double_t        K_1_1270_plus_M;
    Int_t           K_1_1270_plus_ID;
    Double_t        K_1_1270_plus_TAU;
    Double_t        K_1_1270_plus_TAUERR;
    Double_t        K_1_1270_plus_TAUCHI2;
    Bool_t          K_1_1270_plus_L0Global_Dec;
    Bool_t          K_1_1270_plus_L0Global_TIS;
    Bool_t          K_1_1270_plus_L0Global_TOS;
    Bool_t          K_1_1270_plus_Hlt1Global_Dec;
    Bool_t          K_1_1270_plus_Hlt1Global_TIS;
    Bool_t          K_1_1270_plus_Hlt1Global_TOS;
    Bool_t          K_1_1270_plus_Hlt1Phys_Dec;
    Bool_t          K_1_1270_plus_Hlt1Phys_TIS;
    Bool_t          K_1_1270_plus_Hlt1Phys_TOS;
    Bool_t          K_1_1270_plus_Hlt2Global_Dec;
    Bool_t          K_1_1270_plus_Hlt2Global_TIS;
    Bool_t          K_1_1270_plus_Hlt2Global_TOS;
    Bool_t          K_1_1270_plus_Hlt2Phys_Dec;
    Bool_t          K_1_1270_plus_Hlt2Phys_TIS;
    Bool_t          K_1_1270_plus_Hlt2Phys_TOS;
    Bool_t          K_1_1270_plus_L0HadronDecision_Dec;
    Bool_t          K_1_1270_plus_L0HadronDecision_TIS;
    Bool_t          K_1_1270_plus_L0HadronDecision_TOS;
    Bool_t          K_1_1270_plus_L0MuonDecision_Dec;
    Bool_t          K_1_1270_plus_L0MuonDecision_TIS;
    Bool_t          K_1_1270_plus_L0MuonDecision_TOS;
    Bool_t          K_1_1270_plus_L0GlobalDecision_Dec;
    Bool_t          K_1_1270_plus_L0GlobalDecision_TIS;
    Bool_t          K_1_1270_plus_L0GlobalDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt1TrackAllL0Decision_Dec;
    Bool_t          K_1_1270_plus_Hlt1TrackAllL0Decision_TIS;
    Bool_t          K_1_1270_plus_Hlt1TrackAllL0Decision_TOS;
    Bool_t          K_1_1270_plus_Hlt1TrackMVADecision_Dec;
    Bool_t          K_1_1270_plus_Hlt1TrackMVADecision_TIS;
    Bool_t          K_1_1270_plus_Hlt1TrackMVADecision_TOS;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          K_1_1270_plus_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2IncPhiDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2IncPhiDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2IncPhiDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo2BodyDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo3BodyDecision_TOS;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyDecision_Dec;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyDecision_TIS;
    Bool_t          K_1_1270_plus_Hlt2Topo4BodyDecision_TOS;
    Double_t        K_1_1270_plus_cpx_0_50;
    Double_t        K_1_1270_plus_cpy_0_50;
    Double_t        K_1_1270_plus_cpz_0_50;
    Double_t        K_1_1270_plus_cpt_0_50;
    Double_t        K_1_1270_plus_cp_0_50;
    Int_t           K_1_1270_plus_cmult_0_50;
    Double_t        K_1_1270_plus_deltaEta_0_50;
    Double_t        K_1_1270_plus_deltaPhi_0_50;
    Double_t        K_1_1270_plus_pxasy_0_50;
    Double_t        K_1_1270_plus_pyasy_0_50;
    Double_t        K_1_1270_plus_pzasy_0_50;
    Double_t        K_1_1270_plus_pasy_0_50;
    Double_t        K_1_1270_plus_ptasy_0_50;
    Double_t        K_1_1270_plus_cpx_0_60;
    Double_t        K_1_1270_plus_cpy_0_60;
    Double_t        K_1_1270_plus_cpz_0_60;
    Double_t        K_1_1270_plus_cpt_0_60;
    Double_t        K_1_1270_plus_cp_0_60;
    Int_t           K_1_1270_plus_cmult_0_60;
    Double_t        K_1_1270_plus_deltaEta_0_60;
    Double_t        K_1_1270_plus_deltaPhi_0_60;
    Double_t        K_1_1270_plus_pxasy_0_60;
    Double_t        K_1_1270_plus_pyasy_0_60;
    Double_t        K_1_1270_plus_pzasy_0_60;
    Double_t        K_1_1270_plus_pasy_0_60;
    Double_t        K_1_1270_plus_ptasy_0_60;
    Double_t        K_1_1270_plus_cpx_0_70;
    Double_t        K_1_1270_plus_cpy_0_70;
    Double_t        K_1_1270_plus_cpz_0_70;
    Double_t        K_1_1270_plus_cpt_0_70;
    Double_t        K_1_1270_plus_cp_0_70;
    Int_t           K_1_1270_plus_cmult_0_70;
    Double_t        K_1_1270_plus_deltaEta_0_70;
    Double_t        K_1_1270_plus_deltaPhi_0_70;
    Double_t        K_1_1270_plus_pxasy_0_70;
    Double_t        K_1_1270_plus_pyasy_0_70;
    Double_t        K_1_1270_plus_pzasy_0_70;
    Double_t        K_1_1270_plus_pasy_0_70;
    Double_t        K_1_1270_plus_ptasy_0_70;
    Double_t        K_1_1270_plus_cpx_0_80;
    Double_t        K_1_1270_plus_cpy_0_80;
    Double_t        K_1_1270_plus_cpz_0_80;
    Double_t        K_1_1270_plus_cpt_0_80;
    Double_t        K_1_1270_plus_cp_0_80;
    Int_t           K_1_1270_plus_cmult_0_80;
    Double_t        K_1_1270_plus_deltaEta_0_80;
    Double_t        K_1_1270_plus_deltaPhi_0_80;
    Double_t        K_1_1270_plus_pxasy_0_80;
    Double_t        K_1_1270_plus_pyasy_0_80;
    Double_t        K_1_1270_plus_pzasy_0_80;
    Double_t        K_1_1270_plus_pasy_0_80;
    Double_t        K_1_1270_plus_ptasy_0_80;
    Double_t        K_1_1270_plus_cpx_0_90;
    Double_t        K_1_1270_plus_cpy_0_90;
    Double_t        K_1_1270_plus_cpz_0_90;
    Double_t        K_1_1270_plus_cpt_0_90;
    Double_t        K_1_1270_plus_cp_0_90;
    Int_t           K_1_1270_plus_cmult_0_90;
    Double_t        K_1_1270_plus_deltaEta_0_90;
    Double_t        K_1_1270_plus_deltaPhi_0_90;
    Double_t        K_1_1270_plus_pxasy_0_90;
    Double_t        K_1_1270_plus_pyasy_0_90;
    Double_t        K_1_1270_plus_pzasy_0_90;
    Double_t        K_1_1270_plus_pasy_0_90;
    Double_t        K_1_1270_plus_ptasy_0_90;
    Double_t        K_1_1270_plus_cpx_1_00;
    Double_t        K_1_1270_plus_cpy_1_00;
    Double_t        K_1_1270_plus_cpz_1_00;
    Double_t        K_1_1270_plus_cpt_1_00;
    Double_t        K_1_1270_plus_cp_1_00;
    Int_t           K_1_1270_plus_cmult_1_00;
    Double_t        K_1_1270_plus_deltaEta_1_00;
    Double_t        K_1_1270_plus_deltaPhi_1_00;
    Double_t        K_1_1270_plus_pxasy_1_00;
    Double_t        K_1_1270_plus_pyasy_1_00;
    Double_t        K_1_1270_plus_pzasy_1_00;
    Double_t        K_1_1270_plus_pasy_1_00;
    Double_t        K_1_1270_plus_ptasy_1_00;
    Double_t        K_plus_DOCA1;
    Double_t        K_plus_DOCA2;
    Double_t        K_plus_DOCA3;
    Double_t        K_plus_ETA;
    Double_t        K_plus_MC12TuneV2_ProbNNe;
    Double_t        K_plus_MC12TuneV2_ProbNNmu;
    Double_t        K_plus_MC12TuneV2_ProbNNpi;
    Double_t        K_plus_MC12TuneV2_ProbNNk;
    Double_t        K_plus_MC12TuneV2_ProbNNp;
    Double_t        K_plus_MC12TuneV2_ProbNNghost;
    Double_t        K_plus_MC12TuneV3_ProbNNe;
    Double_t        K_plus_MC12TuneV3_ProbNNmu;
    Double_t        K_plus_MC12TuneV3_ProbNNpi;
    Double_t        K_plus_MC12TuneV3_ProbNNk;
    Double_t        K_plus_MC12TuneV3_ProbNNp;
    Double_t        K_plus_MC12TuneV3_ProbNNghost;
    Double_t        K_plus_MC12TuneV4_ProbNNe;
    Double_t        K_plus_MC12TuneV4_ProbNNmu;
    Double_t        K_plus_MC12TuneV4_ProbNNpi;
    Double_t        K_plus_MC12TuneV4_ProbNNk;
    Double_t        K_plus_MC12TuneV4_ProbNNp;
    Double_t        K_plus_MC12TuneV4_ProbNNghost;
    Double_t        K_plus_MC15TuneV1_ProbNNe;
    Double_t        K_plus_MC15TuneV1_ProbNNmu;
    Double_t        K_plus_MC15TuneV1_ProbNNpi;
    Double_t        K_plus_MC15TuneV1_ProbNNk;
    Double_t        K_plus_MC15TuneV1_ProbNNp;
    Double_t        K_plus_MC15TuneV1_ProbNNghost;
    Double_t        K_plus_CosTheta;
    Double_t        K_plus_OWNPV_X;
    Double_t        K_plus_OWNPV_Y;
    Double_t        K_plus_OWNPV_Z;
    Double_t        K_plus_OWNPV_XERR;
    Double_t        K_plus_OWNPV_YERR;
    Double_t        K_plus_OWNPV_ZERR;
    Double_t        K_plus_OWNPV_CHI2;
    Int_t           K_plus_OWNPV_NDOF;
    Float_t         K_plus_OWNPV_COV_[3][3];
    Double_t        K_plus_IP_OWNPV;
    Double_t        K_plus_IPCHI2_OWNPV;
    Double_t        K_plus_ORIVX_X;
    Double_t        K_plus_ORIVX_Y;
    Double_t        K_plus_ORIVX_Z;
    Double_t        K_plus_ORIVX_XERR;
    Double_t        K_plus_ORIVX_YERR;
    Double_t        K_plus_ORIVX_ZERR;
    Double_t        K_plus_ORIVX_CHI2;
    Int_t           K_plus_ORIVX_NDOF;
    Float_t         K_plus_ORIVX_COV_[3][3];
    Double_t        K_plus_P;
    Double_t        K_plus_PT;
    Double_t        K_plus_PE;
    Double_t        K_plus_PX;
    Double_t        K_plus_PY;
    Double_t        K_plus_PZ;
    Double_t        K_plus_M;
    Int_t           K_plus_ID;
    Double_t        K_plus_PIDe;
    Double_t        K_plus_PIDmu;
    Double_t        K_plus_PIDK;
    Double_t        K_plus_PIDp;
    Double_t        K_plus_ProbNNe;
    Double_t        K_plus_ProbNNk;
    Double_t        K_plus_ProbNNp;
    Double_t        K_plus_ProbNNpi;
    Double_t        K_plus_ProbNNmu;
    Double_t        K_plus_ProbNNghost;
    Bool_t          K_plus_hasMuon;
    Bool_t          K_plus_isMuon;
    Bool_t          K_plus_hasRich;
    Bool_t          K_plus_UsedRichAerogel;
    Bool_t          K_plus_UsedRich1Gas;
    Bool_t          K_plus_UsedRich2Gas;
    Bool_t          K_plus_RichAboveElThres;
    Bool_t          K_plus_RichAboveMuThres;
    Bool_t          K_plus_RichAbovePiThres;
    Bool_t          K_plus_RichAboveKaThres;
    Bool_t          K_plus_RichAbovePrThres;
    Bool_t          K_plus_hasCalo;
    Bool_t          K_plus_L0Global_Dec;
    Bool_t          K_plus_L0Global_TIS;
    Bool_t          K_plus_L0Global_TOS;
    Bool_t          K_plus_Hlt1Global_Dec;
    Bool_t          K_plus_Hlt1Global_TIS;
    Bool_t          K_plus_Hlt1Global_TOS;
    Bool_t          K_plus_Hlt1Phys_Dec;
    Bool_t          K_plus_Hlt1Phys_TIS;
    Bool_t          K_plus_Hlt1Phys_TOS;
    Bool_t          K_plus_Hlt2Global_Dec;
    Bool_t          K_plus_Hlt2Global_TIS;
    Bool_t          K_plus_Hlt2Global_TOS;
    Bool_t          K_plus_Hlt2Phys_Dec;
    Bool_t          K_plus_Hlt2Phys_TIS;
    Bool_t          K_plus_Hlt2Phys_TOS;
    Bool_t          K_plus_L0HadronDecision_Dec;
    Bool_t          K_plus_L0HadronDecision_TIS;
    Bool_t          K_plus_L0HadronDecision_TOS;
    Bool_t          K_plus_L0MuonDecision_Dec;
    Bool_t          K_plus_L0MuonDecision_TIS;
    Bool_t          K_plus_L0MuonDecision_TOS;
    Bool_t          K_plus_L0GlobalDecision_Dec;
    Bool_t          K_plus_L0GlobalDecision_TIS;
    Bool_t          K_plus_L0GlobalDecision_TOS;
    Bool_t          K_plus_Hlt1TrackAllL0Decision_Dec;
    Bool_t          K_plus_Hlt1TrackAllL0Decision_TIS;
    Bool_t          K_plus_Hlt1TrackAllL0Decision_TOS;
    Bool_t          K_plus_Hlt1TrackMVADecision_Dec;
    Bool_t          K_plus_Hlt1TrackMVADecision_TIS;
    Bool_t          K_plus_Hlt1TrackMVADecision_TOS;
    Bool_t          K_plus_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          K_plus_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          K_plus_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          K_plus_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          K_plus_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          K_plus_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          K_plus_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          K_plus_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          K_plus_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          K_plus_Hlt2IncPhiDecision_Dec;
    Bool_t          K_plus_Hlt2IncPhiDecision_TIS;
    Bool_t          K_plus_Hlt2IncPhiDecision_TOS;
    Bool_t          K_plus_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          K_plus_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          K_plus_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          K_plus_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          K_plus_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          K_plus_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          K_plus_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          K_plus_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          K_plus_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          K_plus_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          K_plus_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          K_plus_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          K_plus_Hlt2Topo2BodyDecision_Dec;
    Bool_t          K_plus_Hlt2Topo2BodyDecision_TIS;
    Bool_t          K_plus_Hlt2Topo2BodyDecision_TOS;
    Bool_t          K_plus_Hlt2Topo3BodyDecision_Dec;
    Bool_t          K_plus_Hlt2Topo3BodyDecision_TIS;
    Bool_t          K_plus_Hlt2Topo3BodyDecision_TOS;
    Bool_t          K_plus_Hlt2Topo4BodyDecision_Dec;
    Bool_t          K_plus_Hlt2Topo4BodyDecision_TIS;
    Bool_t          K_plus_Hlt2Topo4BodyDecision_TOS;
    Int_t           K_plus_TRACK_Type;
    Int_t           K_plus_TRACK_Key;
    Double_t        K_plus_TRACK_CHI2NDOF;
    Double_t        K_plus_TRACK_PCHI2;
    Double_t        K_plus_TRACK_MatchCHI2;
    Double_t        K_plus_TRACK_GhostProb;
    Double_t        K_plus_TRACK_CloneDist;
    Double_t        K_plus_TRACK_Likelihood;
    Double_t        K_plus_cpx_0_50;
    Double_t        K_plus_cpy_0_50;
    Double_t        K_plus_cpz_0_50;
    Double_t        K_plus_cpt_0_50;
    Double_t        K_plus_cp_0_50;
    Int_t           K_plus_cmult_0_50;
    Double_t        K_plus_deltaEta_0_50;
    Double_t        K_plus_deltaPhi_0_50;
    Double_t        K_plus_pxasy_0_50;
    Double_t        K_plus_pyasy_0_50;
    Double_t        K_plus_pzasy_0_50;
    Double_t        K_plus_pasy_0_50;
    Double_t        K_plus_ptasy_0_50;
    Double_t        K_plus_cpx_0_60;
    Double_t        K_plus_cpy_0_60;
    Double_t        K_plus_cpz_0_60;
    Double_t        K_plus_cpt_0_60;
    Double_t        K_plus_cp_0_60;
    Int_t           K_plus_cmult_0_60;
    Double_t        K_plus_deltaEta_0_60;
    Double_t        K_plus_deltaPhi_0_60;
    Double_t        K_plus_pxasy_0_60;
    Double_t        K_plus_pyasy_0_60;
    Double_t        K_plus_pzasy_0_60;
    Double_t        K_plus_pasy_0_60;
    Double_t        K_plus_ptasy_0_60;
    Double_t        K_plus_cpx_0_70;
    Double_t        K_plus_cpy_0_70;
    Double_t        K_plus_cpz_0_70;
    Double_t        K_plus_cpt_0_70;
    Double_t        K_plus_cp_0_70;
    Int_t           K_plus_cmult_0_70;
    Double_t        K_plus_deltaEta_0_70;
    Double_t        K_plus_deltaPhi_0_70;
    Double_t        K_plus_pxasy_0_70;
    Double_t        K_plus_pyasy_0_70;
    Double_t        K_plus_pzasy_0_70;
    Double_t        K_plus_pasy_0_70;
    Double_t        K_plus_ptasy_0_70;
    Double_t        K_plus_cpx_0_80;
    Double_t        K_plus_cpy_0_80;
    Double_t        K_plus_cpz_0_80;
    Double_t        K_plus_cpt_0_80;
    Double_t        K_plus_cp_0_80;
    Int_t           K_plus_cmult_0_80;
    Double_t        K_plus_deltaEta_0_80;
    Double_t        K_plus_deltaPhi_0_80;
    Double_t        K_plus_pxasy_0_80;
    Double_t        K_plus_pyasy_0_80;
    Double_t        K_plus_pzasy_0_80;
    Double_t        K_plus_pasy_0_80;
    Double_t        K_plus_ptasy_0_80;
    Double_t        K_plus_cpx_0_90;
    Double_t        K_plus_cpy_0_90;
    Double_t        K_plus_cpz_0_90;
    Double_t        K_plus_cpt_0_90;
    Double_t        K_plus_cp_0_90;
    Int_t           K_plus_cmult_0_90;
    Double_t        K_plus_deltaEta_0_90;
    Double_t        K_plus_deltaPhi_0_90;
    Double_t        K_plus_pxasy_0_90;
    Double_t        K_plus_pyasy_0_90;
    Double_t        K_plus_pzasy_0_90;
    Double_t        K_plus_pasy_0_90;
    Double_t        K_plus_ptasy_0_90;
    Double_t        K_plus_cpx_1_00;
    Double_t        K_plus_cpy_1_00;
    Double_t        K_plus_cpz_1_00;
    Double_t        K_plus_cpt_1_00;
    Double_t        K_plus_cp_1_00;
    Int_t           K_plus_cmult_1_00;
    Double_t        K_plus_deltaEta_1_00;
    Double_t        K_plus_deltaPhi_1_00;
    Double_t        K_plus_pxasy_1_00;
    Double_t        K_plus_pyasy_1_00;
    Double_t        K_plus_pzasy_1_00;
    Double_t        K_plus_pasy_1_00;
    Double_t        K_plus_ptasy_1_00;
    Double_t        pi_plus_DOCA1;
    Double_t        pi_plus_DOCA2;
    Double_t        pi_plus_DOCA3;
    Double_t        pi_plus_ETA;
    Double_t        pi_plus_MC12TuneV2_ProbNNe;
    Double_t        pi_plus_MC12TuneV2_ProbNNmu;
    Double_t        pi_plus_MC12TuneV2_ProbNNpi;
    Double_t        pi_plus_MC12TuneV2_ProbNNk;
    Double_t        pi_plus_MC12TuneV2_ProbNNp;
    Double_t        pi_plus_MC12TuneV2_ProbNNghost;
    Double_t        pi_plus_MC12TuneV3_ProbNNe;
    Double_t        pi_plus_MC12TuneV3_ProbNNmu;
    Double_t        pi_plus_MC12TuneV3_ProbNNpi;
    Double_t        pi_plus_MC12TuneV3_ProbNNk;
    Double_t        pi_plus_MC12TuneV3_ProbNNp;
    Double_t        pi_plus_MC12TuneV3_ProbNNghost;
    Double_t        pi_plus_MC12TuneV4_ProbNNe;
    Double_t        pi_plus_MC12TuneV4_ProbNNmu;
    Double_t        pi_plus_MC12TuneV4_ProbNNpi;
    Double_t        pi_plus_MC12TuneV4_ProbNNk;
    Double_t        pi_plus_MC12TuneV4_ProbNNp;
    Double_t        pi_plus_MC12TuneV4_ProbNNghost;
    Double_t        pi_plus_MC15TuneV1_ProbNNe;
    Double_t        pi_plus_MC15TuneV1_ProbNNmu;
    Double_t        pi_plus_MC15TuneV1_ProbNNpi;
    Double_t        pi_plus_MC15TuneV1_ProbNNk;
    Double_t        pi_plus_MC15TuneV1_ProbNNp;
    Double_t        pi_plus_MC15TuneV1_ProbNNghost;
    Double_t        pi_plus_CosTheta;
    Double_t        pi_plus_OWNPV_X;
    Double_t        pi_plus_OWNPV_Y;
    Double_t        pi_plus_OWNPV_Z;
    Double_t        pi_plus_OWNPV_XERR;
    Double_t        pi_plus_OWNPV_YERR;
    Double_t        pi_plus_OWNPV_ZERR;
    Double_t        pi_plus_OWNPV_CHI2;
    Int_t           pi_plus_OWNPV_NDOF;
    Float_t         pi_plus_OWNPV_COV_[3][3];
    Double_t        pi_plus_IP_OWNPV;
    Double_t        pi_plus_IPCHI2_OWNPV;
    Double_t        pi_plus_ORIVX_X;
    Double_t        pi_plus_ORIVX_Y;
    Double_t        pi_plus_ORIVX_Z;
    Double_t        pi_plus_ORIVX_XERR;
    Double_t        pi_plus_ORIVX_YERR;
    Double_t        pi_plus_ORIVX_ZERR;
    Double_t        pi_plus_ORIVX_CHI2;
    Int_t           pi_plus_ORIVX_NDOF;
    Float_t         pi_plus_ORIVX_COV_[3][3];
    Double_t        pi_plus_P;
    Double_t        pi_plus_PT;
    Double_t        pi_plus_PE;
    Double_t        pi_plus_PX;
    Double_t        pi_plus_PY;
    Double_t        pi_plus_PZ;
    Double_t        pi_plus_M;
    Int_t           pi_plus_ID;
    Double_t        pi_plus_PIDe;
    Double_t        pi_plus_PIDmu;
    Double_t        pi_plus_PIDK;
    Double_t        pi_plus_PIDp;
    Double_t        pi_plus_ProbNNe;
    Double_t        pi_plus_ProbNNk;
    Double_t        pi_plus_ProbNNp;
    Double_t        pi_plus_ProbNNpi;
    Double_t        pi_plus_ProbNNmu;
    Double_t        pi_plus_ProbNNghost;
    Bool_t          pi_plus_hasMuon;
    Bool_t          pi_plus_isMuon;
    Bool_t          pi_plus_hasRich;
    Bool_t          pi_plus_UsedRichAerogel;
    Bool_t          pi_plus_UsedRich1Gas;
    Bool_t          pi_plus_UsedRich2Gas;
    Bool_t          pi_plus_RichAboveElThres;
    Bool_t          pi_plus_RichAboveMuThres;
    Bool_t          pi_plus_RichAbovePiThres;
    Bool_t          pi_plus_RichAboveKaThres;
    Bool_t          pi_plus_RichAbovePrThres;
    Bool_t          pi_plus_hasCalo;
    Bool_t          pi_plus_L0Global_Dec;
    Bool_t          pi_plus_L0Global_TIS;
    Bool_t          pi_plus_L0Global_TOS;
    Bool_t          pi_plus_Hlt1Global_Dec;
    Bool_t          pi_plus_Hlt1Global_TIS;
    Bool_t          pi_plus_Hlt1Global_TOS;
    Bool_t          pi_plus_Hlt1Phys_Dec;
    Bool_t          pi_plus_Hlt1Phys_TIS;
    Bool_t          pi_plus_Hlt1Phys_TOS;
    Bool_t          pi_plus_Hlt2Global_Dec;
    Bool_t          pi_plus_Hlt2Global_TIS;
    Bool_t          pi_plus_Hlt2Global_TOS;
    Bool_t          pi_plus_Hlt2Phys_Dec;
    Bool_t          pi_plus_Hlt2Phys_TIS;
    Bool_t          pi_plus_Hlt2Phys_TOS;
    Bool_t          pi_plus_L0HadronDecision_Dec;
    Bool_t          pi_plus_L0HadronDecision_TIS;
    Bool_t          pi_plus_L0HadronDecision_TOS;
    Bool_t          pi_plus_L0MuonDecision_Dec;
    Bool_t          pi_plus_L0MuonDecision_TIS;
    Bool_t          pi_plus_L0MuonDecision_TOS;
    Bool_t          pi_plus_L0GlobalDecision_Dec;
    Bool_t          pi_plus_L0GlobalDecision_TIS;
    Bool_t          pi_plus_L0GlobalDecision_TOS;
    Bool_t          pi_plus_Hlt1TrackAllL0Decision_Dec;
    Bool_t          pi_plus_Hlt1TrackAllL0Decision_TIS;
    Bool_t          pi_plus_Hlt1TrackAllL0Decision_TOS;
    Bool_t          pi_plus_Hlt1TrackMVADecision_Dec;
    Bool_t          pi_plus_Hlt1TrackMVADecision_TIS;
    Bool_t          pi_plus_Hlt1TrackMVADecision_TOS;
    Bool_t          pi_plus_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          pi_plus_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          pi_plus_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          pi_plus_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          pi_plus_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          pi_plus_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          pi_plus_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          pi_plus_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          pi_plus_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          pi_plus_Hlt2IncPhiDecision_Dec;
    Bool_t          pi_plus_Hlt2IncPhiDecision_TIS;
    Bool_t          pi_plus_Hlt2IncPhiDecision_TOS;
    Bool_t          pi_plus_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          pi_plus_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          pi_plus_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo2BodyDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo2BodyDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo2BodyDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo3BodyDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo3BodyDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo3BodyDecision_TOS;
    Bool_t          pi_plus_Hlt2Topo4BodyDecision_Dec;
    Bool_t          pi_plus_Hlt2Topo4BodyDecision_TIS;
    Bool_t          pi_plus_Hlt2Topo4BodyDecision_TOS;
    Int_t           pi_plus_TRACK_Type;
    Int_t           pi_plus_TRACK_Key;
    Double_t        pi_plus_TRACK_CHI2NDOF;
    Double_t        pi_plus_TRACK_PCHI2;
    Double_t        pi_plus_TRACK_MatchCHI2;
    Double_t        pi_plus_TRACK_GhostProb;
    Double_t        pi_plus_TRACK_CloneDist;
    Double_t        pi_plus_TRACK_Likelihood;
    Double_t        pi_plus_cpx_0_50;
    Double_t        pi_plus_cpy_0_50;
    Double_t        pi_plus_cpz_0_50;
    Double_t        pi_plus_cpt_0_50;
    Double_t        pi_plus_cp_0_50;
    Int_t           pi_plus_cmult_0_50;
    Double_t        pi_plus_deltaEta_0_50;
    Double_t        pi_plus_deltaPhi_0_50;
    Double_t        pi_plus_pxasy_0_50;
    Double_t        pi_plus_pyasy_0_50;
    Double_t        pi_plus_pzasy_0_50;
    Double_t        pi_plus_pasy_0_50;
    Double_t        pi_plus_ptasy_0_50;
    Double_t        pi_plus_cpx_0_60;
    Double_t        pi_plus_cpy_0_60;
    Double_t        pi_plus_cpz_0_60;
    Double_t        pi_plus_cpt_0_60;
    Double_t        pi_plus_cp_0_60;
    Int_t           pi_plus_cmult_0_60;
    Double_t        pi_plus_deltaEta_0_60;
    Double_t        pi_plus_deltaPhi_0_60;
    Double_t        pi_plus_pxasy_0_60;
    Double_t        pi_plus_pyasy_0_60;
    Double_t        pi_plus_pzasy_0_60;
    Double_t        pi_plus_pasy_0_60;
    Double_t        pi_plus_ptasy_0_60;
    Double_t        pi_plus_cpx_0_70;
    Double_t        pi_plus_cpy_0_70;
    Double_t        pi_plus_cpz_0_70;
    Double_t        pi_plus_cpt_0_70;
    Double_t        pi_plus_cp_0_70;
    Int_t           pi_plus_cmult_0_70;
    Double_t        pi_plus_deltaEta_0_70;
    Double_t        pi_plus_deltaPhi_0_70;
    Double_t        pi_plus_pxasy_0_70;
    Double_t        pi_plus_pyasy_0_70;
    Double_t        pi_plus_pzasy_0_70;
    Double_t        pi_plus_pasy_0_70;
    Double_t        pi_plus_ptasy_0_70;
    Double_t        pi_plus_cpx_0_80;
    Double_t        pi_plus_cpy_0_80;
    Double_t        pi_plus_cpz_0_80;
    Double_t        pi_plus_cpt_0_80;
    Double_t        pi_plus_cp_0_80;
    Int_t           pi_plus_cmult_0_80;
    Double_t        pi_plus_deltaEta_0_80;
    Double_t        pi_plus_deltaPhi_0_80;
    Double_t        pi_plus_pxasy_0_80;
    Double_t        pi_plus_pyasy_0_80;
    Double_t        pi_plus_pzasy_0_80;
    Double_t        pi_plus_pasy_0_80;
    Double_t        pi_plus_ptasy_0_80;
    Double_t        pi_plus_cpx_0_90;
    Double_t        pi_plus_cpy_0_90;
    Double_t        pi_plus_cpz_0_90;
    Double_t        pi_plus_cpt_0_90;
    Double_t        pi_plus_cp_0_90;
    Int_t           pi_plus_cmult_0_90;
    Double_t        pi_plus_deltaEta_0_90;
    Double_t        pi_plus_deltaPhi_0_90;
    Double_t        pi_plus_pxasy_0_90;
    Double_t        pi_plus_pyasy_0_90;
    Double_t        pi_plus_pzasy_0_90;
    Double_t        pi_plus_pasy_0_90;
    Double_t        pi_plus_ptasy_0_90;
    Double_t        pi_plus_cpx_1_00;
    Double_t        pi_plus_cpy_1_00;
    Double_t        pi_plus_cpz_1_00;
    Double_t        pi_plus_cpt_1_00;
    Double_t        pi_plus_cp_1_00;
    Int_t           pi_plus_cmult_1_00;
    Double_t        pi_plus_deltaEta_1_00;
    Double_t        pi_plus_deltaPhi_1_00;
    Double_t        pi_plus_pxasy_1_00;
    Double_t        pi_plus_pyasy_1_00;
    Double_t        pi_plus_pzasy_1_00;
    Double_t        pi_plus_pasy_1_00;
    Double_t        pi_plus_ptasy_1_00;
    Double_t        pi_minus_DOCA1;
    Double_t        pi_minus_DOCA2;
    Double_t        pi_minus_DOCA3;
    Double_t        pi_minus_ETA;
    Double_t        pi_minus_MC12TuneV2_ProbNNe;
    Double_t        pi_minus_MC12TuneV2_ProbNNmu;
    Double_t        pi_minus_MC12TuneV2_ProbNNpi;
    Double_t        pi_minus_MC12TuneV2_ProbNNk;
    Double_t        pi_minus_MC12TuneV2_ProbNNp;
    Double_t        pi_minus_MC12TuneV2_ProbNNghost;
    Double_t        pi_minus_MC12TuneV3_ProbNNe;
    Double_t        pi_minus_MC12TuneV3_ProbNNmu;
    Double_t        pi_minus_MC12TuneV3_ProbNNpi;
    Double_t        pi_minus_MC12TuneV3_ProbNNk;
    Double_t        pi_minus_MC12TuneV3_ProbNNp;
    Double_t        pi_minus_MC12TuneV3_ProbNNghost;
    Double_t        pi_minus_MC12TuneV4_ProbNNe;
    Double_t        pi_minus_MC12TuneV4_ProbNNmu;
    Double_t        pi_minus_MC12TuneV4_ProbNNpi;
    Double_t        pi_minus_MC12TuneV4_ProbNNk;
    Double_t        pi_minus_MC12TuneV4_ProbNNp;
    Double_t        pi_minus_MC12TuneV4_ProbNNghost;
    Double_t        pi_minus_MC15TuneV1_ProbNNe;
    Double_t        pi_minus_MC15TuneV1_ProbNNmu;
    Double_t        pi_minus_MC15TuneV1_ProbNNpi;
    Double_t        pi_minus_MC15TuneV1_ProbNNk;
    Double_t        pi_minus_MC15TuneV1_ProbNNp;
    Double_t        pi_minus_MC15TuneV1_ProbNNghost;
    Double_t        pi_minus_CosTheta;
    Double_t        pi_minus_OWNPV_X;
    Double_t        pi_minus_OWNPV_Y;
    Double_t        pi_minus_OWNPV_Z;
    Double_t        pi_minus_OWNPV_XERR;
    Double_t        pi_minus_OWNPV_YERR;
    Double_t        pi_minus_OWNPV_ZERR;
    Double_t        pi_minus_OWNPV_CHI2;
    Int_t           pi_minus_OWNPV_NDOF;
    Float_t         pi_minus_OWNPV_COV_[3][3];
    Double_t        pi_minus_IP_OWNPV;
    Double_t        pi_minus_IPCHI2_OWNPV;
    Double_t        pi_minus_ORIVX_X;
    Double_t        pi_minus_ORIVX_Y;
    Double_t        pi_minus_ORIVX_Z;
    Double_t        pi_minus_ORIVX_XERR;
    Double_t        pi_minus_ORIVX_YERR;
    Double_t        pi_minus_ORIVX_ZERR;
    Double_t        pi_minus_ORIVX_CHI2;
    Int_t           pi_minus_ORIVX_NDOF;
    Float_t         pi_minus_ORIVX_COV_[3][3];
    Double_t        pi_minus_P;
    Double_t        pi_minus_PT;
    Double_t        pi_minus_PE;
    Double_t        pi_minus_PX;
    Double_t        pi_minus_PY;
    Double_t        pi_minus_PZ;
    Double_t        pi_minus_M;
    Int_t           pi_minus_ID;
    Double_t        pi_minus_PIDe;
    Double_t        pi_minus_PIDmu;
    Double_t        pi_minus_PIDK;
    Double_t        pi_minus_PIDp;
    Double_t        pi_minus_ProbNNe;
    Double_t        pi_minus_ProbNNk;
    Double_t        pi_minus_ProbNNp;
    Double_t        pi_minus_ProbNNpi;
    Double_t        pi_minus_ProbNNmu;
    Double_t        pi_minus_ProbNNghost;
    Bool_t          pi_minus_hasMuon;
    Bool_t          pi_minus_isMuon;
    Bool_t          pi_minus_hasRich;
    Bool_t          pi_minus_UsedRichAerogel;
    Bool_t          pi_minus_UsedRich1Gas;
    Bool_t          pi_minus_UsedRich2Gas;
    Bool_t          pi_minus_RichAboveElThres;
    Bool_t          pi_minus_RichAboveMuThres;
    Bool_t          pi_minus_RichAbovePiThres;
    Bool_t          pi_minus_RichAboveKaThres;
    Bool_t          pi_minus_RichAbovePrThres;
    Bool_t          pi_minus_hasCalo;
    Bool_t          pi_minus_L0Global_Dec;
    Bool_t          pi_minus_L0Global_TIS;
    Bool_t          pi_minus_L0Global_TOS;
    Bool_t          pi_minus_Hlt1Global_Dec;
    Bool_t          pi_minus_Hlt1Global_TIS;
    Bool_t          pi_minus_Hlt1Global_TOS;
    Bool_t          pi_minus_Hlt1Phys_Dec;
    Bool_t          pi_minus_Hlt1Phys_TIS;
    Bool_t          pi_minus_Hlt1Phys_TOS;
    Bool_t          pi_minus_Hlt2Global_Dec;
    Bool_t          pi_minus_Hlt2Global_TIS;
    Bool_t          pi_minus_Hlt2Global_TOS;
    Bool_t          pi_minus_Hlt2Phys_Dec;
    Bool_t          pi_minus_Hlt2Phys_TIS;
    Bool_t          pi_minus_Hlt2Phys_TOS;
    Bool_t          pi_minus_L0HadronDecision_Dec;
    Bool_t          pi_minus_L0HadronDecision_TIS;
    Bool_t          pi_minus_L0HadronDecision_TOS;
    Bool_t          pi_minus_L0MuonDecision_Dec;
    Bool_t          pi_minus_L0MuonDecision_TIS;
    Bool_t          pi_minus_L0MuonDecision_TOS;
    Bool_t          pi_minus_L0GlobalDecision_Dec;
    Bool_t          pi_minus_L0GlobalDecision_TIS;
    Bool_t          pi_minus_L0GlobalDecision_TOS;
    Bool_t          pi_minus_Hlt1TrackAllL0Decision_Dec;
    Bool_t          pi_minus_Hlt1TrackAllL0Decision_TIS;
    Bool_t          pi_minus_Hlt1TrackAllL0Decision_TOS;
    Bool_t          pi_minus_Hlt1TrackMVADecision_Dec;
    Bool_t          pi_minus_Hlt1TrackMVADecision_TIS;
    Bool_t          pi_minus_Hlt1TrackMVADecision_TOS;
    Bool_t          pi_minus_Hlt1TwoTrackMVADecision_Dec;
    Bool_t          pi_minus_Hlt1TwoTrackMVADecision_TIS;
    Bool_t          pi_minus_Hlt1TwoTrackMVADecision_TOS;
    Bool_t          pi_minus_Hlt1TrackMVALooseDecision_Dec;
    Bool_t          pi_minus_Hlt1TrackMVALooseDecision_TIS;
    Bool_t          pi_minus_Hlt1TrackMVALooseDecision_TOS;
    Bool_t          pi_minus_Hlt1TwoTrackMVALooseDecision_Dec;
    Bool_t          pi_minus_Hlt1TwoTrackMVALooseDecision_TIS;
    Bool_t          pi_minus_Hlt1TwoTrackMVALooseDecision_TOS;
    Bool_t          pi_minus_Hlt2IncPhiDecision_Dec;
    Bool_t          pi_minus_Hlt2IncPhiDecision_TIS;
    Bool_t          pi_minus_Hlt2IncPhiDecision_TOS;
    Bool_t          pi_minus_Hlt2PhiIncPhiDecision_Dec;
    Bool_t          pi_minus_Hlt2PhiIncPhiDecision_TIS;
    Bool_t          pi_minus_Hlt2PhiIncPhiDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo4BodyBBDTDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo2BodyDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo2BodyDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo2BodyDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo3BodyDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo3BodyDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo3BodyDecision_TOS;
    Bool_t          pi_minus_Hlt2Topo4BodyDecision_Dec;
    Bool_t          pi_minus_Hlt2Topo4BodyDecision_TIS;
    Bool_t          pi_minus_Hlt2Topo4BodyDecision_TOS;
    Int_t           pi_minus_TRACK_Type;
    Int_t           pi_minus_TRACK_Key;
    Double_t        pi_minus_TRACK_CHI2NDOF;
    Double_t        pi_minus_TRACK_PCHI2;
    Double_t        pi_minus_TRACK_MatchCHI2;
    Double_t        pi_minus_TRACK_GhostProb;
    Double_t        pi_minus_TRACK_CloneDist;
    Double_t        pi_minus_TRACK_Likelihood;
    Double_t        pi_minus_cpx_0_50;
    Double_t        pi_minus_cpy_0_50;
    Double_t        pi_minus_cpz_0_50;
    Double_t        pi_minus_cpt_0_50;
    Double_t        pi_minus_cp_0_50;
    Int_t           pi_minus_cmult_0_50;
    Double_t        pi_minus_deltaEta_0_50;
    Double_t        pi_minus_deltaPhi_0_50;
    Double_t        pi_minus_pxasy_0_50;
    Double_t        pi_minus_pyasy_0_50;
    Double_t        pi_minus_pzasy_0_50;
    Double_t        pi_minus_pasy_0_50;
    Double_t        pi_minus_ptasy_0_50;
    Double_t        pi_minus_cpx_0_60;
    Double_t        pi_minus_cpy_0_60;
    Double_t        pi_minus_cpz_0_60;
    Double_t        pi_minus_cpt_0_60;
    Double_t        pi_minus_cp_0_60;
    Int_t           pi_minus_cmult_0_60;
    Double_t        pi_minus_deltaEta_0_60;
    Double_t        pi_minus_deltaPhi_0_60;
    Double_t        pi_minus_pxasy_0_60;
    Double_t        pi_minus_pyasy_0_60;
    Double_t        pi_minus_pzasy_0_60;
    Double_t        pi_minus_pasy_0_60;
    Double_t        pi_minus_ptasy_0_60;
    Double_t        pi_minus_cpx_0_70;
    Double_t        pi_minus_cpy_0_70;
    Double_t        pi_minus_cpz_0_70;
    Double_t        pi_minus_cpt_0_70;
    Double_t        pi_minus_cp_0_70;
    Int_t           pi_minus_cmult_0_70;
    Double_t        pi_minus_deltaEta_0_70;
    Double_t        pi_minus_deltaPhi_0_70;
    Double_t        pi_minus_pxasy_0_70;
    Double_t        pi_minus_pyasy_0_70;
    Double_t        pi_minus_pzasy_0_70;
    Double_t        pi_minus_pasy_0_70;
    Double_t        pi_minus_ptasy_0_70;
    Double_t        pi_minus_cpx_0_80;
    Double_t        pi_minus_cpy_0_80;
    Double_t        pi_minus_cpz_0_80;
    Double_t        pi_minus_cpt_0_80;
    Double_t        pi_minus_cp_0_80;
    Int_t           pi_minus_cmult_0_80;
    Double_t        pi_minus_deltaEta_0_80;
    Double_t        pi_minus_deltaPhi_0_80;
    Double_t        pi_minus_pxasy_0_80;
    Double_t        pi_minus_pyasy_0_80;
    Double_t        pi_minus_pzasy_0_80;
    Double_t        pi_minus_pasy_0_80;
    Double_t        pi_minus_ptasy_0_80;
    Double_t        pi_minus_cpx_0_90;
    Double_t        pi_minus_cpy_0_90;
    Double_t        pi_minus_cpz_0_90;
    Double_t        pi_minus_cpt_0_90;
    Double_t        pi_minus_cp_0_90;
    Int_t           pi_minus_cmult_0_90;
    Double_t        pi_minus_deltaEta_0_90;
    Double_t        pi_minus_deltaPhi_0_90;
    Double_t        pi_minus_pxasy_0_90;
    Double_t        pi_minus_pyasy_0_90;
    Double_t        pi_minus_pzasy_0_90;
    Double_t        pi_minus_pasy_0_90;
    Double_t        pi_minus_ptasy_0_90;
    Double_t        pi_minus_cpx_1_00;
    Double_t        pi_minus_cpy_1_00;
    Double_t        pi_minus_cpz_1_00;
    Double_t        pi_minus_cpt_1_00;
    Double_t        pi_minus_cp_1_00;
    Int_t           pi_minus_cmult_1_00;
    Double_t        pi_minus_deltaEta_1_00;
    Double_t        pi_minus_deltaPhi_1_00;
    Double_t        pi_minus_pxasy_1_00;
    Double_t        pi_minus_pyasy_1_00;
    Double_t        pi_minus_pzasy_1_00;
    Double_t        pi_minus_pasy_1_00;
    Double_t        pi_minus_ptasy_1_00;
    UInt_t          nCandidate;
    ULong64_t       totCandidates;
    ULong64_t       EventInSequence;
    UInt_t          runNumber;
    ULong64_t       eventNumber;
    UInt_t          BCID;
    Int_t           BCType;
    UInt_t          OdinTCK;
    UInt_t          L0DUTCK;
    UInt_t          HLT1TCK;
    UInt_t          HLT2TCK;
    ULong64_t       GpsTime;
    Short_t         Polarity;
    Int_t           nPV;
    Float_t         PVX[100];   //[nPV]
    Float_t         PVY[100];   //[nPV]
    Float_t         PVZ[100];   //[nPV]
    Float_t         PVXERR[100];   //[nPV]
    Float_t         PVYERR[100];   //[nPV]
    Float_t         PVZERR[100];   //[nPV]
    Float_t         PVCHI2[100];   //[nPV]
    Float_t         PVNDOF[100];   //[nPV]
    Float_t         PVNTRACKS[100];   //[nPV]
    Int_t           nPVs;
    Int_t           nTracks;
    Int_t           nLongTracks;
    Int_t           nDownstreamTracks;
    Int_t           nUpstreamTracks;
    Int_t           nVeloTracks;
    Int_t           nTTracks;
    Int_t           nBackTracks;
    Int_t           nRich1Hits;
    Int_t           nRich2Hits;
    Int_t           nVeloClusters;
    Int_t           nITClusters;
    Int_t           nTTClusters;
    Int_t           nOTClusters;
    Int_t           nSPDHits;
    Int_t           nMuonCoordsS0;
    Int_t           nMuonCoordsS1;
    Int_t           nMuonCoordsS2;
    Int_t           nMuonCoordsS3;
    Int_t           nMuonCoordsS4;
    Int_t           nMuonTracks;
    
    // List of branches
    TBranch        *b_Bs_DOCA1;   //!
    TBranch        *b_Bs_DOCA2;   //!
    TBranch        *b_Bs_DOCA3;   //!
    TBranch        *b_Bs_ETA;   //!
    TBranch        *b_Bs_ENDVERTEX_X;   //!
    TBranch        *b_Bs_ENDVERTEX_Y;   //!
    TBranch        *b_Bs_ENDVERTEX_Z;   //!
    TBranch        *b_Bs_ENDVERTEX_XERR;   //!
    TBranch        *b_Bs_ENDVERTEX_YERR;   //!
    TBranch        *b_Bs_ENDVERTEX_ZERR;   //!
    TBranch        *b_Bs_ENDVERTEX_CHI2;   //!
    TBranch        *b_Bs_ENDVERTEX_NDOF;   //!
    TBranch        *b_Bs_ENDVERTEX_COV_;   //!
    TBranch        *b_Bs_OWNPV_X;   //!
    TBranch        *b_Bs_OWNPV_Y;   //!
    TBranch        *b_Bs_OWNPV_Z;   //!
    TBranch        *b_Bs_OWNPV_XERR;   //!
    TBranch        *b_Bs_OWNPV_YERR;   //!
    TBranch        *b_Bs_OWNPV_ZERR;   //!
    TBranch        *b_Bs_OWNPV_CHI2;   //!
    TBranch        *b_Bs_OWNPV_NDOF;   //!
    TBranch        *b_Bs_OWNPV_COV_;   //!
    TBranch        *b_Bs_IP_OWNPV;   //!
    TBranch        *b_Bs_IPCHI2_OWNPV;   //!
    TBranch        *b_Bs_FD_OWNPV;   //!
    TBranch        *b_Bs_FDCHI2_OWNPV;   //!
    TBranch        *b_Bs_DIRA_OWNPV;   //!
    TBranch        *b_Bs_P;   //!
    TBranch        *b_Bs_PT;   //!
    TBranch        *b_Bs_PE;   //!
    TBranch        *b_Bs_PX;   //!
    TBranch        *b_Bs_PY;   //!
    TBranch        *b_Bs_PZ;   //!
    TBranch        *b_Bs_MM;   //!
    TBranch        *b_Bs_MMERR;   //!
    TBranch        *b_Bs_M;   //!
    TBranch        *b_Bs_ID;   //!
    TBranch        *b_Bs_TAU;   //!
    TBranch        *b_Bs_TAUERR;   //!
    TBranch        *b_Bs_TAUCHI2;   //!
    TBranch        *b_Bs_L0Global_Dec;   //!
    TBranch        *b_Bs_L0Global_TIS;   //!
    TBranch        *b_Bs_L0Global_TOS;   //!
    TBranch        *b_Bs_Hlt1Global_Dec;   //!
    TBranch        *b_Bs_Hlt1Global_TIS;   //!
    TBranch        *b_Bs_Hlt1Global_TOS;   //!
    TBranch        *b_Bs_Hlt1Phys_Dec;   //!
    TBranch        *b_Bs_Hlt1Phys_TIS;   //!
    TBranch        *b_Bs_Hlt1Phys_TOS;   //!
    TBranch        *b_Bs_Hlt2Global_Dec;   //!
    TBranch        *b_Bs_Hlt2Global_TIS;   //!
    TBranch        *b_Bs_Hlt2Global_TOS;   //!
    TBranch        *b_Bs_Hlt2Phys_Dec;   //!
    TBranch        *b_Bs_Hlt2Phys_TIS;   //!
    TBranch        *b_Bs_Hlt2Phys_TOS;   //!
    TBranch        *b_Bs_L0HadronDecision_Dec;   //!
    TBranch        *b_Bs_L0HadronDecision_TIS;   //!
    TBranch        *b_Bs_L0HadronDecision_TOS;   //!
    TBranch        *b_Bs_L0MuonDecision_Dec;   //!
    TBranch        *b_Bs_L0MuonDecision_TIS;   //!
    TBranch        *b_Bs_L0MuonDecision_TOS;   //!
    TBranch        *b_Bs_L0GlobalDecision_Dec;   //!
    TBranch        *b_Bs_L0GlobalDecision_TIS;   //!
    TBranch        *b_Bs_L0GlobalDecision_TOS;   //!
    TBranch        *b_Bs_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_Bs_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_Bs_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_Bs_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_Bs_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_Bs_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_Bs_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_Bs_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_Bs_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_Bs_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_Bs_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_Bs_TAGDECISION;   //!
    TBranch        *b_Bs_TAGOMEGA;   //!
    TBranch        *b_Bs_TAGDECISION_OS;   //!
    TBranch        *b_Bs_TAGOMEGA_OS;   //!
    TBranch        *b_Bs_TAGGER;   //!
    TBranch        *b_Bs_OS_Muon_DEC;   //!
    TBranch        *b_Bs_OS_Muon_PROB;   //!
    TBranch        *b_Bs_OS_Electron_DEC;   //!
    TBranch        *b_Bs_OS_Electron_PROB;   //!
    TBranch        *b_Bs_OS_Kaon_DEC;   //!
    TBranch        *b_Bs_OS_Kaon_PROB;   //!
    TBranch        *b_Bs_SS_Kaon_DEC;   //!
    TBranch        *b_Bs_SS_Kaon_PROB;   //!
    TBranch        *b_Bs_SS_Pion_DEC;   //!
    TBranch        *b_Bs_SS_Pion_PROB;   //!
    TBranch        *b_Bs_SS_PionBDT_DEC;   //!
    TBranch        *b_Bs_SS_PionBDT_PROB;   //!
    TBranch        *b_Bs_VtxCharge_DEC;   //!
    TBranch        *b_Bs_VtxCharge_PROB;   //!
    TBranch        *b_Bs_OS_nnetKaon_DEC;   //!
    TBranch        *b_Bs_OS_nnetKaon_PROB;   //!
    TBranch        *b_Bs_SS_nnetKaon_DEC;   //!
    TBranch        *b_Bs_SS_nnetKaon_PROB;   //!
    TBranch        *b_Bs_SS_Proton_DEC;   //!
    TBranch        *b_Bs_SS_Proton_PROB;   //!
    TBranch        *b_Bs_OS_Charm_DEC;   //!
    TBranch        *b_Bs_OS_Charm_PROB;   //!
    TBranch        *b_Bs_cpx_0_50;   //!
    TBranch        *b_Bs_cpy_0_50;   //!
    TBranch        *b_Bs_cpz_0_50;   //!
    TBranch        *b_Bs_cpt_0_50;   //!
    TBranch        *b_Bs_cp_0_50;   //!
    TBranch        *b_Bs_cmult_0_50;   //!
    TBranch        *b_Bs_deltaEta_0_50;   //!
    TBranch        *b_Bs_deltaPhi_0_50;   //!
    TBranch        *b_Bs_pxasy_0_50;   //!
    TBranch        *b_Bs_pyasy_0_50;   //!
    TBranch        *b_Bs_pzasy_0_50;   //!
    TBranch        *b_Bs_pasy_0_50;   //!
    TBranch        *b_Bs_ptasy_0_50;   //!
    TBranch        *b_Bs_cpx_0_60;   //!
    TBranch        *b_Bs_cpy_0_60;   //!
    TBranch        *b_Bs_cpz_0_60;   //!
    TBranch        *b_Bs_cpt_0_60;   //!
    TBranch        *b_Bs_cp_0_60;   //!
    TBranch        *b_Bs_cmult_0_60;   //!
    TBranch        *b_Bs_deltaEta_0_60;   //!
    TBranch        *b_Bs_deltaPhi_0_60;   //!
    TBranch        *b_Bs_pxasy_0_60;   //!
    TBranch        *b_Bs_pyasy_0_60;   //!
    TBranch        *b_Bs_pzasy_0_60;   //!
    TBranch        *b_Bs_pasy_0_60;   //!
    TBranch        *b_Bs_ptasy_0_60;   //!
    TBranch        *b_Bs_cpx_0_70;   //!
    TBranch        *b_Bs_cpy_0_70;   //!
    TBranch        *b_Bs_cpz_0_70;   //!
    TBranch        *b_Bs_cpt_0_70;   //!
    TBranch        *b_Bs_cp_0_70;   //!
    TBranch        *b_Bs_cmult_0_70;   //!
    TBranch        *b_Bs_deltaEta_0_70;   //!
    TBranch        *b_Bs_deltaPhi_0_70;   //!
    TBranch        *b_Bs_pxasy_0_70;   //!
    TBranch        *b_Bs_pyasy_0_70;   //!
    TBranch        *b_Bs_pzasy_0_70;   //!
    TBranch        *b_Bs_pasy_0_70;   //!
    TBranch        *b_Bs_ptasy_0_70;   //!
    TBranch        *b_Bs_cpx_0_80;   //!
    TBranch        *b_Bs_cpy_0_80;   //!
    TBranch        *b_Bs_cpz_0_80;   //!
    TBranch        *b_Bs_cpt_0_80;   //!
    TBranch        *b_Bs_cp_0_80;   //!
    TBranch        *b_Bs_cmult_0_80;   //!
    TBranch        *b_Bs_deltaEta_0_80;   //!
    TBranch        *b_Bs_deltaPhi_0_80;   //!
    TBranch        *b_Bs_pxasy_0_80;   //!
    TBranch        *b_Bs_pyasy_0_80;   //!
    TBranch        *b_Bs_pzasy_0_80;   //!
    TBranch        *b_Bs_pasy_0_80;   //!
    TBranch        *b_Bs_ptasy_0_80;   //!
    TBranch        *b_Bs_cpx_0_90;   //!
    TBranch        *b_Bs_cpy_0_90;   //!
    TBranch        *b_Bs_cpz_0_90;   //!
    TBranch        *b_Bs_cpt_0_90;   //!
    TBranch        *b_Bs_cp_0_90;   //!
    TBranch        *b_Bs_cmult_0_90;   //!
    TBranch        *b_Bs_deltaEta_0_90;   //!
    TBranch        *b_Bs_deltaPhi_0_90;   //!
    TBranch        *b_Bs_pxasy_0_90;   //!
    TBranch        *b_Bs_pyasy_0_90;   //!
    TBranch        *b_Bs_pzasy_0_90;   //!
    TBranch        *b_Bs_pasy_0_90;   //!
    TBranch        *b_Bs_ptasy_0_90;   //!
    TBranch        *b_Bs_cpx_1_00;   //!
    TBranch        *b_Bs_cpy_1_00;   //!
    TBranch        *b_Bs_cpz_1_00;   //!
    TBranch        *b_Bs_cpt_1_00;   //!
    TBranch        *b_Bs_cp_1_00;   //!
    TBranch        *b_Bs_cmult_1_00;   //!
    TBranch        *b_Bs_deltaEta_1_00;   //!
    TBranch        *b_Bs_deltaPhi_1_00;   //!
    TBranch        *b_Bs_pxasy_1_00;   //!
    TBranch        *b_Bs_pyasy_1_00;   //!
    TBranch        *b_Bs_pzasy_1_00;   //!
    TBranch        *b_Bs_pasy_1_00;   //!
    TBranch        *b_Bs_ptasy_1_00;   //!
    TBranch        *b_Bs_B0DTF_nPV;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_0_ID;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_0_PE;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_0_PX;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_0_PY;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_0_PZ;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_ID;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_PE;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_PX;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_PY;   //!
    TBranch        *b_Bs_B0DTF_D_splus_Kplus_PZ;   //!
    TBranch        *b_Bs_B0DTF_D_splus_M;   //!
    TBranch        *b_Bs_B0DTF_D_splus_MERR;   //!
    TBranch        *b_Bs_B0DTF_D_splus_P;   //!
    TBranch        *b_Bs_B0DTF_D_splus_PERR;   //!
    TBranch        *b_Bs_B0DTF_D_splus_ctau;   //!
    TBranch        *b_Bs_B0DTF_D_splus_ctauErr;   //!
    TBranch        *b_Bs_B0DTF_D_splus_decayLength;   //!
    TBranch        *b_Bs_B0DTF_D_splus_decayLengthErr;   //!
    TBranch        *b_Bs_B0DTF_D_splus_piplus_ID;   //!
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PE;   //!
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PX;   //!
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PY;   //!
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PZ;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_Kplus_ID;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_Kplus_PE;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_Kplus_PX;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_Kplus_PY;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_Kplus_PZ;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_M;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_MERR;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_P;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_PERR;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_ctau;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_ctauErr;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_decayLength;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_decayLengthErr;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_0_ID;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_0_PE;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_0_PX;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_0_PY;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_0_PZ;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_ID;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_PE;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_PX;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_PY;   //!
    TBranch        *b_Bs_B0DTF_K_1_1270_plus_piplus_PZ;   //!
    TBranch        *b_Bs_B0DTF_M;   //!
    TBranch        *b_Bs_B0DTF_MERR;   //!
    TBranch        *b_Bs_B0DTF_P;   //!
    TBranch        *b_Bs_B0DTF_PERR;   //!
    TBranch        *b_Bs_B0DTF_PV_X;   //!
    TBranch        *b_Bs_B0DTF_PV_Y;   //!
    TBranch        *b_Bs_B0DTF_PV_Z;   //!
    TBranch        *b_Bs_B0DTF_PV_key;   //!
    TBranch        *b_Bs_B0DTF_chi2;   //!
    TBranch        *b_Bs_B0DTF_ctau;   //!
    TBranch        *b_Bs_B0DTF_ctauErr;   //!
    TBranch        *b_Bs_B0DTF_decayLength;   //!
    TBranch        *b_Bs_B0DTF_decayLengthErr;   //!
    TBranch        *b_Bs_B0DTF_nDOF;   //!
    TBranch        *b_Bs_B0DTF_nIter;   //!
    TBranch        *b_Bs_B0DTF_status;   //!
    TBranch        *b_Bs_BsDTF_nPV;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_0_ID;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_0_PE;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_0_PX;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_0_PY;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_0_PZ;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_ID;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_PE;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_PX;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_PY;   //!
    TBranch        *b_Bs_BsDTF_D_splus_Kplus_PZ;   //!
    TBranch        *b_Bs_BsDTF_D_splus_M;   //!
    TBranch        *b_Bs_BsDTF_D_splus_MERR;   //!
    TBranch        *b_Bs_BsDTF_D_splus_P;   //!
    TBranch        *b_Bs_BsDTF_D_splus_PERR;   //!
    TBranch        *b_Bs_BsDTF_D_splus_ctau;   //!
    TBranch        *b_Bs_BsDTF_D_splus_ctauErr;   //!
    TBranch        *b_Bs_BsDTF_D_splus_decayLength;   //!
    TBranch        *b_Bs_BsDTF_D_splus_decayLengthErr;   //!
    TBranch        *b_Bs_BsDTF_D_splus_piplus_ID;   //!
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PE;   //!
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PX;   //!
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PY;   //!
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PZ;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_Kplus_ID;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_Kplus_PE;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_Kplus_PX;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_Kplus_PY;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_Kplus_PZ;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_M;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_MERR;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_P;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_PERR;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_ctau;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_ctauErr;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_decayLength;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_decayLengthErr;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_0_ID;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_0_PE;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_0_PX;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_0_PY;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_0_PZ;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_ID;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_PE;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_PX;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_PY;   //!
    TBranch        *b_Bs_BsDTF_K_1_1270_plus_piplus_PZ;   //!
    TBranch        *b_Bs_BsDTF_M;   //!
    TBranch        *b_Bs_BsDTF_MERR;   //!
    TBranch        *b_Bs_BsDTF_P;   //!
    TBranch        *b_Bs_BsDTF_PERR;   //!
    TBranch        *b_Bs_BsDTF_PV_X;   //!
    TBranch        *b_Bs_BsDTF_PV_Y;   //!
    TBranch        *b_Bs_BsDTF_PV_Z;   //!
    TBranch        *b_Bs_BsDTF_PV_key;   //!
    TBranch        *b_Bs_BsDTF_chi2;   //!
    TBranch        *b_Bs_BsDTF_ctau;   //!
    TBranch        *b_Bs_BsDTF_ctauErr;   //!
    TBranch        *b_Bs_BsDTF_decayLength;   //!
    TBranch        *b_Bs_BsDTF_decayLengthErr;   //!
    TBranch        *b_Bs_BsDTF_nDOF;   //!
    TBranch        *b_Bs_BsDTF_nIter;   //!
    TBranch        *b_Bs_BsDTF_status;   //!
    TBranch        *b_Bs_DTF_nPV;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_0_ID;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_0_PE;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_0_PX;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_0_PY;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_0_PZ;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_ID;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_PE;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_PX;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_PY;   //!
    TBranch        *b_Bs_DTF_D_splus_Kplus_PZ;   //!
    TBranch        *b_Bs_DTF_D_splus_M;   //!
    TBranch        *b_Bs_DTF_D_splus_MERR;   //!
    TBranch        *b_Bs_DTF_D_splus_P;   //!
    TBranch        *b_Bs_DTF_D_splus_PERR;   //!
    TBranch        *b_Bs_DTF_D_splus_ctau;   //!
    TBranch        *b_Bs_DTF_D_splus_ctauErr;   //!
    TBranch        *b_Bs_DTF_D_splus_decayLength;   //!
    TBranch        *b_Bs_DTF_D_splus_decayLengthErr;   //!
    TBranch        *b_Bs_DTF_D_splus_piplus_ID;   //!
    TBranch        *b_Bs_DTF_D_splus_piplus_PE;   //!
    TBranch        *b_Bs_DTF_D_splus_piplus_PX;   //!
    TBranch        *b_Bs_DTF_D_splus_piplus_PY;   //!
    TBranch        *b_Bs_DTF_D_splus_piplus_PZ;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_Kplus_ID;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_Kplus_PE;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_Kplus_PX;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_Kplus_PY;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_Kplus_PZ;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_M;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_MERR;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_P;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_PERR;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_ctau;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_ctauErr;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_decayLength;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_decayLengthErr;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_0_ID;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_0_PE;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_0_PX;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_0_PY;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_0_PZ;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_ID;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_PE;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_PX;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_PY;   //!
    TBranch        *b_Bs_DTF_K_1_1270_plus_piplus_PZ;   //!
    TBranch        *b_Bs_DTF_M;   //!
    TBranch        *b_Bs_DTF_MERR;   //!
    TBranch        *b_Bs_DTF_P;   //!
    TBranch        *b_Bs_DTF_PERR;   //!
    TBranch        *b_Bs_DTF_PV_X;   //!
    TBranch        *b_Bs_DTF_PV_Y;   //!
    TBranch        *b_Bs_DTF_PV_Z;   //!
    TBranch        *b_Bs_DTF_PV_key;   //!
    TBranch        *b_Bs_DTF_chi2;   //!
    TBranch        *b_Bs_DTF_ctau;   //!
    TBranch        *b_Bs_DTF_ctauErr;   //!
    TBranch        *b_Bs_DTF_decayLength;   //!
    TBranch        *b_Bs_DTF_decayLengthErr;   //!
    TBranch        *b_Bs_DTF_nDOF;   //!
    TBranch        *b_Bs_DTF_nIter;   //!
    TBranch        *b_Bs_DTF_status;   //!
    TBranch        *b_Bs_PV_nPV;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_0_ID;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_0_PE;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_0_PX;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_0_PY;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_0_PZ;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_ID;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_PE;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_PX;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_PY;   //!
    TBranch        *b_Bs_PV_Dplus_Kplus_PZ;   //!
    TBranch        *b_Bs_PV_Dplus_M;   //!
    TBranch        *b_Bs_PV_Dplus_MERR;   //!
    TBranch        *b_Bs_PV_Dplus_P;   //!
    TBranch        *b_Bs_PV_Dplus_PERR;   //!
    TBranch        *b_Bs_PV_Dplus_ctau;   //!
    TBranch        *b_Bs_PV_Dplus_ctauErr;   //!
    TBranch        *b_Bs_PV_Dplus_decayLength;   //!
    TBranch        *b_Bs_PV_Dplus_decayLengthErr;   //!
    TBranch        *b_Bs_PV_Dplus_piplus_ID;   //!
    TBranch        *b_Bs_PV_Dplus_piplus_PE;   //!
    TBranch        *b_Bs_PV_Dplus_piplus_PX;   //!
    TBranch        *b_Bs_PV_Dplus_piplus_PY;   //!
    TBranch        *b_Bs_PV_Dplus_piplus_PZ;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_Kplus_ID;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_Kplus_PE;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_Kplus_PX;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_Kplus_PY;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_Kplus_PZ;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_M;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_MERR;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_P;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_PERR;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_ctau;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_ctauErr;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_decayLength;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_decayLengthErr;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_0_ID;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_0_PE;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_0_PX;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_0_PY;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_0_PZ;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_ID;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_PE;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_PX;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_PY;   //!
    TBranch        *b_Bs_PV_K_1_1270_plus_piplus_PZ;   //!
    TBranch        *b_Bs_PV_M;   //!
    TBranch        *b_Bs_PV_MERR;   //!
    TBranch        *b_Bs_PV_P;   //!
    TBranch        *b_Bs_PV_PERR;   //!
    TBranch        *b_Bs_PV_PV_X;   //!
    TBranch        *b_Bs_PV_PV_Y;   //!
    TBranch        *b_Bs_PV_PV_Z;   //!
    TBranch        *b_Bs_PV_PV_key;   //!
    TBranch        *b_Bs_PV_chi2;   //!
    TBranch        *b_Bs_PV_ctau;   //!
    TBranch        *b_Bs_PV_ctauErr;   //!
    TBranch        *b_Bs_PV_decayLength;   //!
    TBranch        *b_Bs_PV_decayLengthErr;   //!
    TBranch        *b_Bs_PV_nDOF;   //!
    TBranch        *b_Bs_PV_nIter;   //!
    TBranch        *b_Bs_PV_status;   //!
    TBranch        *b_Bs_BsTaggingTool_TAGDECISION;   //!
    TBranch        *b_Bs_BsTaggingTool_TAGOMEGA;   //!
    TBranch        *b_Bs_BsTaggingTool_TAGDECISION_OS;   //!
    TBranch        *b_Bs_BsTaggingTool_TAGOMEGA_OS;   //!
    TBranch        *b_Bs_BsTaggingTool_TAGGER;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Muon_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Muon_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Electron_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Electron_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Kaon_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Kaon_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Kaon_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Kaon_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Pion_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Pion_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_PionBDT_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_PionBDT_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_VtxCharge_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_VtxCharge_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_nnetKaon_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_nnetKaon_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_nnetKaon_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_nnetKaon_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Proton_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_SS_Proton_PROB;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Charm_DEC;   //!
    TBranch        *b_Bs_BsTaggingTool_OS_Charm_PROB;   //!
    TBranch        *b_Ds_DOCA1;   //!
    TBranch        *b_Ds_DOCA2;   //!
    TBranch        *b_Ds_DOCA3;   //!
    TBranch        *b_Ds_ETA;   //!
    TBranch        *b_Ds_CosTheta;   //!
    TBranch        *b_Ds_ENDVERTEX_X;   //!
    TBranch        *b_Ds_ENDVERTEX_Y;   //!
    TBranch        *b_Ds_ENDVERTEX_Z;   //!
    TBranch        *b_Ds_ENDVERTEX_XERR;   //!
    TBranch        *b_Ds_ENDVERTEX_YERR;   //!
    TBranch        *b_Ds_ENDVERTEX_ZERR;   //!
    TBranch        *b_Ds_ENDVERTEX_CHI2;   //!
    TBranch        *b_Ds_ENDVERTEX_NDOF;   //!
    TBranch        *b_Ds_ENDVERTEX_COV_;   //!
    TBranch        *b_Ds_OWNPV_X;   //!
    TBranch        *b_Ds_OWNPV_Y;   //!
    TBranch        *b_Ds_OWNPV_Z;   //!
    TBranch        *b_Ds_OWNPV_XERR;   //!
    TBranch        *b_Ds_OWNPV_YERR;   //!
    TBranch        *b_Ds_OWNPV_ZERR;   //!
    TBranch        *b_Ds_OWNPV_CHI2;   //!
    TBranch        *b_Ds_OWNPV_NDOF;   //!
    TBranch        *b_Ds_OWNPV_COV_;   //!
    TBranch        *b_Ds_IP_OWNPV;   //!
    TBranch        *b_Ds_IPCHI2_OWNPV;   //!
    TBranch        *b_Ds_FD_OWNPV;   //!
    TBranch        *b_Ds_FDCHI2_OWNPV;   //!
    TBranch        *b_Ds_DIRA_OWNPV;   //!
    TBranch        *b_Ds_ORIVX_X;   //!
    TBranch        *b_Ds_ORIVX_Y;   //!
    TBranch        *b_Ds_ORIVX_Z;   //!
    TBranch        *b_Ds_ORIVX_XERR;   //!
    TBranch        *b_Ds_ORIVX_YERR;   //!
    TBranch        *b_Ds_ORIVX_ZERR;   //!
    TBranch        *b_Ds_ORIVX_CHI2;   //!
    TBranch        *b_Ds_ORIVX_NDOF;   //!
    TBranch        *b_Ds_ORIVX_COV_;   //!
    TBranch        *b_Ds_FD_ORIVX;   //!
    TBranch        *b_Ds_FDCHI2_ORIVX;   //!
    TBranch        *b_Ds_DIRA_ORIVX;   //!
    TBranch        *b_Ds_P;   //!
    TBranch        *b_Ds_PT;   //!
    TBranch        *b_Ds_PE;   //!
    TBranch        *b_Ds_PX;   //!
    TBranch        *b_Ds_PY;   //!
    TBranch        *b_Ds_PZ;   //!
    TBranch        *b_Ds_MM;   //!
    TBranch        *b_Ds_MMERR;   //!
    TBranch        *b_Ds_M;   //!
    TBranch        *b_Ds_ID;   //!
    TBranch        *b_Ds_TAU;   //!
    TBranch        *b_Ds_TAUERR;   //!
    TBranch        *b_Ds_TAUCHI2;   //!
    TBranch        *b_Ds_L0Global_Dec;   //!
    TBranch        *b_Ds_L0Global_TIS;   //!
    TBranch        *b_Ds_L0Global_TOS;   //!
    TBranch        *b_Ds_Hlt1Global_Dec;   //!
    TBranch        *b_Ds_Hlt1Global_TIS;   //!
    TBranch        *b_Ds_Hlt1Global_TOS;   //!
    TBranch        *b_Ds_Hlt1Phys_Dec;   //!
    TBranch        *b_Ds_Hlt1Phys_TIS;   //!
    TBranch        *b_Ds_Hlt1Phys_TOS;   //!
    TBranch        *b_Ds_Hlt2Global_Dec;   //!
    TBranch        *b_Ds_Hlt2Global_TIS;   //!
    TBranch        *b_Ds_Hlt2Global_TOS;   //!
    TBranch        *b_Ds_Hlt2Phys_Dec;   //!
    TBranch        *b_Ds_Hlt2Phys_TIS;   //!
    TBranch        *b_Ds_Hlt2Phys_TOS;   //!
    TBranch        *b_Ds_L0HadronDecision_Dec;   //!
    TBranch        *b_Ds_L0HadronDecision_TIS;   //!
    TBranch        *b_Ds_L0HadronDecision_TOS;   //!
    TBranch        *b_Ds_L0MuonDecision_Dec;   //!
    TBranch        *b_Ds_L0MuonDecision_TIS;   //!
    TBranch        *b_Ds_L0MuonDecision_TOS;   //!
    TBranch        *b_Ds_L0GlobalDecision_Dec;   //!
    TBranch        *b_Ds_L0GlobalDecision_TIS;   //!
    TBranch        *b_Ds_L0GlobalDecision_TOS;   //!
    TBranch        *b_Ds_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_Ds_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_Ds_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_Ds_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_Ds_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_Ds_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_Ds_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_Ds_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_Ds_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_Ds_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_Ds_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_Ds_cpx_0_50;   //!
    TBranch        *b_Ds_cpy_0_50;   //!
    TBranch        *b_Ds_cpz_0_50;   //!
    TBranch        *b_Ds_cpt_0_50;   //!
    TBranch        *b_Ds_cp_0_50;   //!
    TBranch        *b_Ds_cmult_0_50;   //!
    TBranch        *b_Ds_deltaEta_0_50;   //!
    TBranch        *b_Ds_deltaPhi_0_50;   //!
    TBranch        *b_Ds_pxasy_0_50;   //!
    TBranch        *b_Ds_pyasy_0_50;   //!
    TBranch        *b_Ds_pzasy_0_50;   //!
    TBranch        *b_Ds_pasy_0_50;   //!
    TBranch        *b_Ds_ptasy_0_50;   //!
    TBranch        *b_Ds_cpx_0_60;   //!
    TBranch        *b_Ds_cpy_0_60;   //!
    TBranch        *b_Ds_cpz_0_60;   //!
    TBranch        *b_Ds_cpt_0_60;   //!
    TBranch        *b_Ds_cp_0_60;   //!
    TBranch        *b_Ds_cmult_0_60;   //!
    TBranch        *b_Ds_deltaEta_0_60;   //!
    TBranch        *b_Ds_deltaPhi_0_60;   //!
    TBranch        *b_Ds_pxasy_0_60;   //!
    TBranch        *b_Ds_pyasy_0_60;   //!
    TBranch        *b_Ds_pzasy_0_60;   //!
    TBranch        *b_Ds_pasy_0_60;   //!
    TBranch        *b_Ds_ptasy_0_60;   //!
    TBranch        *b_Ds_cpx_0_70;   //!
    TBranch        *b_Ds_cpy_0_70;   //!
    TBranch        *b_Ds_cpz_0_70;   //!
    TBranch        *b_Ds_cpt_0_70;   //!
    TBranch        *b_Ds_cp_0_70;   //!
    TBranch        *b_Ds_cmult_0_70;   //!
    TBranch        *b_Ds_deltaEta_0_70;   //!
    TBranch        *b_Ds_deltaPhi_0_70;   //!
    TBranch        *b_Ds_pxasy_0_70;   //!
    TBranch        *b_Ds_pyasy_0_70;   //!
    TBranch        *b_Ds_pzasy_0_70;   //!
    TBranch        *b_Ds_pasy_0_70;   //!
    TBranch        *b_Ds_ptasy_0_70;   //!
    TBranch        *b_Ds_cpx_0_80;   //!
    TBranch        *b_Ds_cpy_0_80;   //!
    TBranch        *b_Ds_cpz_0_80;   //!
    TBranch        *b_Ds_cpt_0_80;   //!
    TBranch        *b_Ds_cp_0_80;   //!
    TBranch        *b_Ds_cmult_0_80;   //!
    TBranch        *b_Ds_deltaEta_0_80;   //!
    TBranch        *b_Ds_deltaPhi_0_80;   //!
    TBranch        *b_Ds_pxasy_0_80;   //!
    TBranch        *b_Ds_pyasy_0_80;   //!
    TBranch        *b_Ds_pzasy_0_80;   //!
    TBranch        *b_Ds_pasy_0_80;   //!
    TBranch        *b_Ds_ptasy_0_80;   //!
    TBranch        *b_Ds_cpx_0_90;   //!
    TBranch        *b_Ds_cpy_0_90;   //!
    TBranch        *b_Ds_cpz_0_90;   //!
    TBranch        *b_Ds_cpt_0_90;   //!
    TBranch        *b_Ds_cp_0_90;   //!
    TBranch        *b_Ds_cmult_0_90;   //!
    TBranch        *b_Ds_deltaEta_0_90;   //!
    TBranch        *b_Ds_deltaPhi_0_90;   //!
    TBranch        *b_Ds_pxasy_0_90;   //!
    TBranch        *b_Ds_pyasy_0_90;   //!
    TBranch        *b_Ds_pzasy_0_90;   //!
    TBranch        *b_Ds_pasy_0_90;   //!
    TBranch        *b_Ds_ptasy_0_90;   //!
    TBranch        *b_Ds_cpx_1_00;   //!
    TBranch        *b_Ds_cpy_1_00;   //!
    TBranch        *b_Ds_cpz_1_00;   //!
    TBranch        *b_Ds_cpt_1_00;   //!
    TBranch        *b_Ds_cp_1_00;   //!
    TBranch        *b_Ds_cmult_1_00;   //!
    TBranch        *b_Ds_deltaEta_1_00;   //!
    TBranch        *b_Ds_deltaPhi_1_00;   //!
    TBranch        *b_Ds_pxasy_1_00;   //!
    TBranch        *b_Ds_pyasy_1_00;   //!
    TBranch        *b_Ds_pzasy_1_00;   //!
    TBranch        *b_Ds_pasy_1_00;   //!
    TBranch        *b_Ds_ptasy_1_00;   //!
    TBranch        *b_K_plus_fromDs_DOCA1;   //!
    TBranch        *b_K_plus_fromDs_DOCA2;   //!
    TBranch        *b_K_plus_fromDs_DOCA3;   //!
    TBranch        *b_K_plus_fromDs_ETA;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_K_plus_fromDs_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_K_plus_fromDs_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_K_plus_fromDs_CosTheta;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_X;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_Y;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_Z;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_XERR;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_YERR;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_ZERR;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_CHI2;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_NDOF;   //!
    TBranch        *b_K_plus_fromDs_OWNPV_COV_;   //!
    TBranch        *b_K_plus_fromDs_IP_OWNPV;   //!
    TBranch        *b_K_plus_fromDs_IPCHI2_OWNPV;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_X;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_Y;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_Z;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_XERR;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_YERR;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_ZERR;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_CHI2;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_NDOF;   //!
    TBranch        *b_K_plus_fromDs_ORIVX_COV_;   //!
    TBranch        *b_K_plus_fromDs_P;   //!
    TBranch        *b_K_plus_fromDs_PT;   //!
    TBranch        *b_K_plus_fromDs_PE;   //!
    TBranch        *b_K_plus_fromDs_PX;   //!
    TBranch        *b_K_plus_fromDs_PY;   //!
    TBranch        *b_K_plus_fromDs_PZ;   //!
    TBranch        *b_K_plus_fromDs_M;   //!
    TBranch        *b_K_plus_fromDs_ID;   //!
    TBranch        *b_K_plus_fromDs_PIDe;   //!
    TBranch        *b_K_plus_fromDs_PIDmu;   //!
    TBranch        *b_K_plus_fromDs_PIDK;   //!
    TBranch        *b_K_plus_fromDs_PIDp;   //!
    TBranch        *b_K_plus_fromDs_ProbNNe;   //!
    TBranch        *b_K_plus_fromDs_ProbNNk;   //!
    TBranch        *b_K_plus_fromDs_ProbNNp;   //!
    TBranch        *b_K_plus_fromDs_ProbNNpi;   //!
    TBranch        *b_K_plus_fromDs_ProbNNmu;   //!
    TBranch        *b_K_plus_fromDs_ProbNNghost;   //!
    TBranch        *b_K_plus_fromDs_hasMuon;   //!
    TBranch        *b_K_plus_fromDs_isMuon;   //!
    TBranch        *b_K_plus_fromDs_hasRich;   //!
    TBranch        *b_K_plus_fromDs_UsedRichAerogel;   //!
    TBranch        *b_K_plus_fromDs_UsedRich1Gas;   //!
    TBranch        *b_K_plus_fromDs_UsedRich2Gas;   //!
    TBranch        *b_K_plus_fromDs_RichAboveElThres;   //!
    TBranch        *b_K_plus_fromDs_RichAboveMuThres;   //!
    TBranch        *b_K_plus_fromDs_RichAbovePiThres;   //!
    TBranch        *b_K_plus_fromDs_RichAboveKaThres;   //!
    TBranch        *b_K_plus_fromDs_RichAbovePrThres;   //!
    TBranch        *b_K_plus_fromDs_hasCalo;   //!
    TBranch        *b_K_plus_fromDs_L0Global_Dec;   //!
    TBranch        *b_K_plus_fromDs_L0Global_TIS;   //!
    TBranch        *b_K_plus_fromDs_L0Global_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Global_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Global_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Global_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Phys_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Phys_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1Phys_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Global_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Global_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Global_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Phys_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Phys_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Phys_TOS;   //!
    TBranch        *b_K_plus_fromDs_L0HadronDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_L0HadronDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_L0HadronDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_L0MuonDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_L0MuonDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_L0MuonDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_L0GlobalDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_L0GlobalDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_L0GlobalDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_K_plus_fromDs_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_K_plus_fromDs_TRACK_Type;   //!
    TBranch        *b_K_plus_fromDs_TRACK_Key;   //!
    TBranch        *b_K_plus_fromDs_TRACK_CHI2NDOF;   //!
    TBranch        *b_K_plus_fromDs_TRACK_PCHI2;   //!
    TBranch        *b_K_plus_fromDs_TRACK_MatchCHI2;   //!
    TBranch        *b_K_plus_fromDs_TRACK_GhostProb;   //!
    TBranch        *b_K_plus_fromDs_TRACK_CloneDist;   //!
    TBranch        *b_K_plus_fromDs_TRACK_Likelihood;   //!
    TBranch        *b_K_plus_fromDs_cpx_0_50;   //!
    TBranch        *b_K_plus_fromDs_cpy_0_50;   //!
    TBranch        *b_K_plus_fromDs_cpz_0_50;   //!
    TBranch        *b_K_plus_fromDs_cpt_0_50;   //!
    TBranch        *b_K_plus_fromDs_cp_0_50;   //!
    TBranch        *b_K_plus_fromDs_cmult_0_50;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_0_50;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_0_50;   //!
    TBranch        *b_K_plus_fromDs_pxasy_0_50;   //!
    TBranch        *b_K_plus_fromDs_pyasy_0_50;   //!
    TBranch        *b_K_plus_fromDs_pzasy_0_50;   //!
    TBranch        *b_K_plus_fromDs_pasy_0_50;   //!
    TBranch        *b_K_plus_fromDs_ptasy_0_50;   //!
    TBranch        *b_K_plus_fromDs_cpx_0_60;   //!
    TBranch        *b_K_plus_fromDs_cpy_0_60;   //!
    TBranch        *b_K_plus_fromDs_cpz_0_60;   //!
    TBranch        *b_K_plus_fromDs_cpt_0_60;   //!
    TBranch        *b_K_plus_fromDs_cp_0_60;   //!
    TBranch        *b_K_plus_fromDs_cmult_0_60;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_0_60;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_0_60;   //!
    TBranch        *b_K_plus_fromDs_pxasy_0_60;   //!
    TBranch        *b_K_plus_fromDs_pyasy_0_60;   //!
    TBranch        *b_K_plus_fromDs_pzasy_0_60;   //!
    TBranch        *b_K_plus_fromDs_pasy_0_60;   //!
    TBranch        *b_K_plus_fromDs_ptasy_0_60;   //!
    TBranch        *b_K_plus_fromDs_cpx_0_70;   //!
    TBranch        *b_K_plus_fromDs_cpy_0_70;   //!
    TBranch        *b_K_plus_fromDs_cpz_0_70;   //!
    TBranch        *b_K_plus_fromDs_cpt_0_70;   //!
    TBranch        *b_K_plus_fromDs_cp_0_70;   //!
    TBranch        *b_K_plus_fromDs_cmult_0_70;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_0_70;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_0_70;   //!
    TBranch        *b_K_plus_fromDs_pxasy_0_70;   //!
    TBranch        *b_K_plus_fromDs_pyasy_0_70;   //!
    TBranch        *b_K_plus_fromDs_pzasy_0_70;   //!
    TBranch        *b_K_plus_fromDs_pasy_0_70;   //!
    TBranch        *b_K_plus_fromDs_ptasy_0_70;   //!
    TBranch        *b_K_plus_fromDs_cpx_0_80;   //!
    TBranch        *b_K_plus_fromDs_cpy_0_80;   //!
    TBranch        *b_K_plus_fromDs_cpz_0_80;   //!
    TBranch        *b_K_plus_fromDs_cpt_0_80;   //!
    TBranch        *b_K_plus_fromDs_cp_0_80;   //!
    TBranch        *b_K_plus_fromDs_cmult_0_80;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_0_80;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_0_80;   //!
    TBranch        *b_K_plus_fromDs_pxasy_0_80;   //!
    TBranch        *b_K_plus_fromDs_pyasy_0_80;   //!
    TBranch        *b_K_plus_fromDs_pzasy_0_80;   //!
    TBranch        *b_K_plus_fromDs_pasy_0_80;   //!
    TBranch        *b_K_plus_fromDs_ptasy_0_80;   //!
    TBranch        *b_K_plus_fromDs_cpx_0_90;   //!
    TBranch        *b_K_plus_fromDs_cpy_0_90;   //!
    TBranch        *b_K_plus_fromDs_cpz_0_90;   //!
    TBranch        *b_K_plus_fromDs_cpt_0_90;   //!
    TBranch        *b_K_plus_fromDs_cp_0_90;   //!
    TBranch        *b_K_plus_fromDs_cmult_0_90;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_0_90;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_0_90;   //!
    TBranch        *b_K_plus_fromDs_pxasy_0_90;   //!
    TBranch        *b_K_plus_fromDs_pyasy_0_90;   //!
    TBranch        *b_K_plus_fromDs_pzasy_0_90;   //!
    TBranch        *b_K_plus_fromDs_pasy_0_90;   //!
    TBranch        *b_K_plus_fromDs_ptasy_0_90;   //!
    TBranch        *b_K_plus_fromDs_cpx_1_00;   //!
    TBranch        *b_K_plus_fromDs_cpy_1_00;   //!
    TBranch        *b_K_plus_fromDs_cpz_1_00;   //!
    TBranch        *b_K_plus_fromDs_cpt_1_00;   //!
    TBranch        *b_K_plus_fromDs_cp_1_00;   //!
    TBranch        *b_K_plus_fromDs_cmult_1_00;   //!
    TBranch        *b_K_plus_fromDs_deltaEta_1_00;   //!
    TBranch        *b_K_plus_fromDs_deltaPhi_1_00;   //!
    TBranch        *b_K_plus_fromDs_pxasy_1_00;   //!
    TBranch        *b_K_plus_fromDs_pyasy_1_00;   //!
    TBranch        *b_K_plus_fromDs_pzasy_1_00;   //!
    TBranch        *b_K_plus_fromDs_pasy_1_00;   //!
    TBranch        *b_K_plus_fromDs_ptasy_1_00;   //!
    TBranch        *b_K_minus_fromDs_DOCA1;   //!
    TBranch        *b_K_minus_fromDs_DOCA2;   //!
    TBranch        *b_K_minus_fromDs_DOCA3;   //!
    TBranch        *b_K_minus_fromDs_ETA;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_K_minus_fromDs_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_K_minus_fromDs_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_K_minus_fromDs_CosTheta;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_X;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_Y;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_Z;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_XERR;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_YERR;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_ZERR;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_CHI2;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_NDOF;   //!
    TBranch        *b_K_minus_fromDs_OWNPV_COV_;   //!
    TBranch        *b_K_minus_fromDs_IP_OWNPV;   //!
    TBranch        *b_K_minus_fromDs_IPCHI2_OWNPV;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_X;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_Y;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_Z;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_XERR;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_YERR;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_ZERR;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_CHI2;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_NDOF;   //!
    TBranch        *b_K_minus_fromDs_ORIVX_COV_;   //!
    TBranch        *b_K_minus_fromDs_P;   //!
    TBranch        *b_K_minus_fromDs_PT;   //!
    TBranch        *b_K_minus_fromDs_PE;   //!
    TBranch        *b_K_minus_fromDs_PX;   //!
    TBranch        *b_K_minus_fromDs_PY;   //!
    TBranch        *b_K_minus_fromDs_PZ;   //!
    TBranch        *b_K_minus_fromDs_M;   //!
    TBranch        *b_K_minus_fromDs_ID;   //!
    TBranch        *b_K_minus_fromDs_PIDe;   //!
    TBranch        *b_K_minus_fromDs_PIDmu;   //!
    TBranch        *b_K_minus_fromDs_PIDK;   //!
    TBranch        *b_K_minus_fromDs_PIDp;   //!
    TBranch        *b_K_minus_fromDs_ProbNNe;   //!
    TBranch        *b_K_minus_fromDs_ProbNNk;   //!
    TBranch        *b_K_minus_fromDs_ProbNNp;   //!
    TBranch        *b_K_minus_fromDs_ProbNNpi;   //!
    TBranch        *b_K_minus_fromDs_ProbNNmu;   //!
    TBranch        *b_K_minus_fromDs_ProbNNghost;   //!
    TBranch        *b_K_minus_fromDs_hasMuon;   //!
    TBranch        *b_K_minus_fromDs_isMuon;   //!
    TBranch        *b_K_minus_fromDs_hasRich;   //!
    TBranch        *b_K_minus_fromDs_UsedRichAerogel;   //!
    TBranch        *b_K_minus_fromDs_UsedRich1Gas;   //!
    TBranch        *b_K_minus_fromDs_UsedRich2Gas;   //!
    TBranch        *b_K_minus_fromDs_RichAboveElThres;   //!
    TBranch        *b_K_minus_fromDs_RichAboveMuThres;   //!
    TBranch        *b_K_minus_fromDs_RichAbovePiThres;   //!
    TBranch        *b_K_minus_fromDs_RichAboveKaThres;   //!
    TBranch        *b_K_minus_fromDs_RichAbovePrThres;   //!
    TBranch        *b_K_minus_fromDs_hasCalo;   //!
    TBranch        *b_K_minus_fromDs_L0Global_Dec;   //!
    TBranch        *b_K_minus_fromDs_L0Global_TIS;   //!
    TBranch        *b_K_minus_fromDs_L0Global_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Global_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Global_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Global_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Phys_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Phys_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1Phys_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Global_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Global_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Global_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Phys_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Phys_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Phys_TOS;   //!
    TBranch        *b_K_minus_fromDs_L0HadronDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_L0HadronDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_L0HadronDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_L0MuonDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_L0MuonDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_L0MuonDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_L0GlobalDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_L0GlobalDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_L0GlobalDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_K_minus_fromDs_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_K_minus_fromDs_TRACK_Type;   //!
    TBranch        *b_K_minus_fromDs_TRACK_Key;   //!
    TBranch        *b_K_minus_fromDs_TRACK_CHI2NDOF;   //!
    TBranch        *b_K_minus_fromDs_TRACK_PCHI2;   //!
    TBranch        *b_K_minus_fromDs_TRACK_MatchCHI2;   //!
    TBranch        *b_K_minus_fromDs_TRACK_GhostProb;   //!
    TBranch        *b_K_minus_fromDs_TRACK_CloneDist;   //!
    TBranch        *b_K_minus_fromDs_TRACK_Likelihood;   //!
    TBranch        *b_K_minus_fromDs_cpx_0_50;   //!
    TBranch        *b_K_minus_fromDs_cpy_0_50;   //!
    TBranch        *b_K_minus_fromDs_cpz_0_50;   //!
    TBranch        *b_K_minus_fromDs_cpt_0_50;   //!
    TBranch        *b_K_minus_fromDs_cp_0_50;   //!
    TBranch        *b_K_minus_fromDs_cmult_0_50;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_0_50;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_0_50;   //!
    TBranch        *b_K_minus_fromDs_pxasy_0_50;   //!
    TBranch        *b_K_minus_fromDs_pyasy_0_50;   //!
    TBranch        *b_K_minus_fromDs_pzasy_0_50;   //!
    TBranch        *b_K_minus_fromDs_pasy_0_50;   //!
    TBranch        *b_K_minus_fromDs_ptasy_0_50;   //!
    TBranch        *b_K_minus_fromDs_cpx_0_60;   //!
    TBranch        *b_K_minus_fromDs_cpy_0_60;   //!
    TBranch        *b_K_minus_fromDs_cpz_0_60;   //!
    TBranch        *b_K_minus_fromDs_cpt_0_60;   //!
    TBranch        *b_K_minus_fromDs_cp_0_60;   //!
    TBranch        *b_K_minus_fromDs_cmult_0_60;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_0_60;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_0_60;   //!
    TBranch        *b_K_minus_fromDs_pxasy_0_60;   //!
    TBranch        *b_K_minus_fromDs_pyasy_0_60;   //!
    TBranch        *b_K_minus_fromDs_pzasy_0_60;   //!
    TBranch        *b_K_minus_fromDs_pasy_0_60;   //!
    TBranch        *b_K_minus_fromDs_ptasy_0_60;   //!
    TBranch        *b_K_minus_fromDs_cpx_0_70;   //!
    TBranch        *b_K_minus_fromDs_cpy_0_70;   //!
    TBranch        *b_K_minus_fromDs_cpz_0_70;   //!
    TBranch        *b_K_minus_fromDs_cpt_0_70;   //!
    TBranch        *b_K_minus_fromDs_cp_0_70;   //!
    TBranch        *b_K_minus_fromDs_cmult_0_70;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_0_70;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_0_70;   //!
    TBranch        *b_K_minus_fromDs_pxasy_0_70;   //!
    TBranch        *b_K_minus_fromDs_pyasy_0_70;   //!
    TBranch        *b_K_minus_fromDs_pzasy_0_70;   //!
    TBranch        *b_K_minus_fromDs_pasy_0_70;   //!
    TBranch        *b_K_minus_fromDs_ptasy_0_70;   //!
    TBranch        *b_K_minus_fromDs_cpx_0_80;   //!
    TBranch        *b_K_minus_fromDs_cpy_0_80;   //!
    TBranch        *b_K_minus_fromDs_cpz_0_80;   //!
    TBranch        *b_K_minus_fromDs_cpt_0_80;   //!
    TBranch        *b_K_minus_fromDs_cp_0_80;   //!
    TBranch        *b_K_minus_fromDs_cmult_0_80;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_0_80;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_0_80;   //!
    TBranch        *b_K_minus_fromDs_pxasy_0_80;   //!
    TBranch        *b_K_minus_fromDs_pyasy_0_80;   //!
    TBranch        *b_K_minus_fromDs_pzasy_0_80;   //!
    TBranch        *b_K_minus_fromDs_pasy_0_80;   //!
    TBranch        *b_K_minus_fromDs_ptasy_0_80;   //!
    TBranch        *b_K_minus_fromDs_cpx_0_90;   //!
    TBranch        *b_K_minus_fromDs_cpy_0_90;   //!
    TBranch        *b_K_minus_fromDs_cpz_0_90;   //!
    TBranch        *b_K_minus_fromDs_cpt_0_90;   //!
    TBranch        *b_K_minus_fromDs_cp_0_90;   //!
    TBranch        *b_K_minus_fromDs_cmult_0_90;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_0_90;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_0_90;   //!
    TBranch        *b_K_minus_fromDs_pxasy_0_90;   //!
    TBranch        *b_K_minus_fromDs_pyasy_0_90;   //!
    TBranch        *b_K_minus_fromDs_pzasy_0_90;   //!
    TBranch        *b_K_minus_fromDs_pasy_0_90;   //!
    TBranch        *b_K_minus_fromDs_ptasy_0_90;   //!
    TBranch        *b_K_minus_fromDs_cpx_1_00;   //!
    TBranch        *b_K_minus_fromDs_cpy_1_00;   //!
    TBranch        *b_K_minus_fromDs_cpz_1_00;   //!
    TBranch        *b_K_minus_fromDs_cpt_1_00;   //!
    TBranch        *b_K_minus_fromDs_cp_1_00;   //!
    TBranch        *b_K_minus_fromDs_cmult_1_00;   //!
    TBranch        *b_K_minus_fromDs_deltaEta_1_00;   //!
    TBranch        *b_K_minus_fromDs_deltaPhi_1_00;   //!
    TBranch        *b_K_minus_fromDs_pxasy_1_00;   //!
    TBranch        *b_K_minus_fromDs_pyasy_1_00;   //!
    TBranch        *b_K_minus_fromDs_pzasy_1_00;   //!
    TBranch        *b_K_minus_fromDs_pasy_1_00;   //!
    TBranch        *b_K_minus_fromDs_ptasy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_DOCA1;   //!
    TBranch        *b_pi_minus_fromDs_DOCA2;   //!
    TBranch        *b_pi_minus_fromDs_DOCA3;   //!
    TBranch        *b_pi_minus_fromDs_ETA;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_pi_minus_fromDs_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_pi_minus_fromDs_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_pi_minus_fromDs_CosTheta;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_X;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_Y;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_Z;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_XERR;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_YERR;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_ZERR;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_CHI2;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_NDOF;   //!
    TBranch        *b_pi_minus_fromDs_OWNPV_COV_;   //!
    TBranch        *b_pi_minus_fromDs_IP_OWNPV;   //!
    TBranch        *b_pi_minus_fromDs_IPCHI2_OWNPV;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_X;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_Y;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_Z;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_XERR;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_YERR;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_ZERR;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_CHI2;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_NDOF;   //!
    TBranch        *b_pi_minus_fromDs_ORIVX_COV_;   //!
    TBranch        *b_pi_minus_fromDs_P;   //!
    TBranch        *b_pi_minus_fromDs_PT;   //!
    TBranch        *b_pi_minus_fromDs_PE;   //!
    TBranch        *b_pi_minus_fromDs_PX;   //!
    TBranch        *b_pi_minus_fromDs_PY;   //!
    TBranch        *b_pi_minus_fromDs_PZ;   //!
    TBranch        *b_pi_minus_fromDs_M;   //!
    TBranch        *b_pi_minus_fromDs_ID;   //!
    TBranch        *b_pi_minus_fromDs_PIDe;   //!
    TBranch        *b_pi_minus_fromDs_PIDmu;   //!
    TBranch        *b_pi_minus_fromDs_PIDK;   //!
    TBranch        *b_pi_minus_fromDs_PIDp;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNe;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNk;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNp;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNpi;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNmu;   //!
    TBranch        *b_pi_minus_fromDs_ProbNNghost;   //!
    TBranch        *b_pi_minus_fromDs_hasMuon;   //!
    TBranch        *b_pi_minus_fromDs_isMuon;   //!
    TBranch        *b_pi_minus_fromDs_hasRich;   //!
    TBranch        *b_pi_minus_fromDs_UsedRichAerogel;   //!
    TBranch        *b_pi_minus_fromDs_UsedRich1Gas;   //!
    TBranch        *b_pi_minus_fromDs_UsedRich2Gas;   //!
    TBranch        *b_pi_minus_fromDs_RichAboveElThres;   //!
    TBranch        *b_pi_minus_fromDs_RichAboveMuThres;   //!
    TBranch        *b_pi_minus_fromDs_RichAbovePiThres;   //!
    TBranch        *b_pi_minus_fromDs_RichAboveKaThres;   //!
    TBranch        *b_pi_minus_fromDs_RichAbovePrThres;   //!
    TBranch        *b_pi_minus_fromDs_hasCalo;   //!
    TBranch        *b_pi_minus_fromDs_L0Global_Dec;   //!
    TBranch        *b_pi_minus_fromDs_L0Global_TIS;   //!
    TBranch        *b_pi_minus_fromDs_L0Global_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Global_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Global_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Global_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Phys_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Phys_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1Phys_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Global_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Global_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Global_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Phys_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Phys_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Phys_TOS;   //!
    TBranch        *b_pi_minus_fromDs_L0HadronDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_L0HadronDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_L0HadronDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_L0MuonDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_L0MuonDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_L0MuonDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_L0GlobalDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_L0GlobalDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_L0GlobalDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_fromDs_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_Type;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_Key;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_CHI2NDOF;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_PCHI2;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_MatchCHI2;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_GhostProb;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_CloneDist;   //!
    TBranch        *b_pi_minus_fromDs_TRACK_Likelihood;   //!
    TBranch        *b_pi_minus_fromDs_cpx_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cpy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cpz_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cpt_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cp_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cmult_0_50;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_0_50;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_0_50;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_pasy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_0_50;   //!
    TBranch        *b_pi_minus_fromDs_cpx_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cpy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cpz_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cpt_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cp_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cmult_0_60;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_0_60;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_0_60;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_pasy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_0_60;   //!
    TBranch        *b_pi_minus_fromDs_cpx_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cpy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cpz_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cpt_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cp_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cmult_0_70;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_0_70;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_0_70;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_pasy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_0_70;   //!
    TBranch        *b_pi_minus_fromDs_cpx_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cpy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cpz_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cpt_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cp_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cmult_0_80;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_0_80;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_0_80;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_pasy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_0_80;   //!
    TBranch        *b_pi_minus_fromDs_cpx_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cpy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cpz_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cpt_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cp_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cmult_0_90;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_0_90;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_0_90;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_pasy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_0_90;   //!
    TBranch        *b_pi_minus_fromDs_cpx_1_00;   //!
    TBranch        *b_pi_minus_fromDs_cpy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_cpz_1_00;   //!
    TBranch        *b_pi_minus_fromDs_cpt_1_00;   //!
    TBranch        *b_pi_minus_fromDs_cp_1_00;   //!
    TBranch        *b_pi_minus_fromDs_cmult_1_00;   //!
    TBranch        *b_pi_minus_fromDs_deltaEta_1_00;   //!
    TBranch        *b_pi_minus_fromDs_deltaPhi_1_00;   //!
    TBranch        *b_pi_minus_fromDs_pxasy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_pyasy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_pzasy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_pasy_1_00;   //!
    TBranch        *b_pi_minus_fromDs_ptasy_1_00;   //!
    TBranch        *b_K_1_1270_plus_DOCA1;   //!
    TBranch        *b_K_1_1270_plus_DOCA2;   //!
    TBranch        *b_K_1_1270_plus_DOCA3;   //!
    TBranch        *b_K_1_1270_plus_ETA;   //!
    TBranch        *b_K_1_1270_plus_CosTheta;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_X;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_Y;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_Z;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_XERR;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_YERR;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_ZERR;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_CHI2;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_NDOF;   //!
    TBranch        *b_K_1_1270_plus_ENDVERTEX_COV_;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_X;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_Y;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_Z;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_XERR;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_YERR;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_ZERR;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_CHI2;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_NDOF;   //!
    TBranch        *b_K_1_1270_plus_OWNPV_COV_;   //!
    TBranch        *b_K_1_1270_plus_IP_OWNPV;   //!
    TBranch        *b_K_1_1270_plus_IPCHI2_OWNPV;   //!
    TBranch        *b_K_1_1270_plus_FD_OWNPV;   //!
    TBranch        *b_K_1_1270_plus_FDCHI2_OWNPV;   //!
    TBranch        *b_K_1_1270_plus_DIRA_OWNPV;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_X;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_Y;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_Z;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_XERR;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_YERR;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_ZERR;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_CHI2;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_NDOF;   //!
    TBranch        *b_K_1_1270_plus_ORIVX_COV_;   //!
    TBranch        *b_K_1_1270_plus_FD_ORIVX;   //!
    TBranch        *b_K_1_1270_plus_FDCHI2_ORIVX;   //!
    TBranch        *b_K_1_1270_plus_DIRA_ORIVX;   //!
    TBranch        *b_K_1_1270_plus_P;   //!
    TBranch        *b_K_1_1270_plus_PT;   //!
    TBranch        *b_K_1_1270_plus_PE;   //!
    TBranch        *b_K_1_1270_plus_PX;   //!
    TBranch        *b_K_1_1270_plus_PY;   //!
    TBranch        *b_K_1_1270_plus_PZ;   //!
    TBranch        *b_K_1_1270_plus_MM;   //!
    TBranch        *b_K_1_1270_plus_MMERR;   //!
    TBranch        *b_K_1_1270_plus_M;   //!
    TBranch        *b_K_1_1270_plus_ID;   //!
    TBranch        *b_K_1_1270_plus_TAU;   //!
    TBranch        *b_K_1_1270_plus_TAUERR;   //!
    TBranch        *b_K_1_1270_plus_TAUCHI2;   //!
    TBranch        *b_K_1_1270_plus_L0Global_Dec;   //!
    TBranch        *b_K_1_1270_plus_L0Global_TIS;   //!
    TBranch        *b_K_1_1270_plus_L0Global_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Global_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Global_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Global_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Phys_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Phys_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1Phys_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Global_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Global_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Global_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Phys_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Phys_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Phys_TOS;   //!
    TBranch        *b_K_1_1270_plus_L0HadronDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_L0HadronDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_L0HadronDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_L0MuonDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_L0MuonDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_L0MuonDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_L0GlobalDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_L0GlobalDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_L0GlobalDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_K_1_1270_plus_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_K_1_1270_plus_cpx_0_50;   //!
    TBranch        *b_K_1_1270_plus_cpy_0_50;   //!
    TBranch        *b_K_1_1270_plus_cpz_0_50;   //!
    TBranch        *b_K_1_1270_plus_cpt_0_50;   //!
    TBranch        *b_K_1_1270_plus_cp_0_50;   //!
    TBranch        *b_K_1_1270_plus_cmult_0_50;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_0_50;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_0_50;   //!
    TBranch        *b_K_1_1270_plus_pxasy_0_50;   //!
    TBranch        *b_K_1_1270_plus_pyasy_0_50;   //!
    TBranch        *b_K_1_1270_plus_pzasy_0_50;   //!
    TBranch        *b_K_1_1270_plus_pasy_0_50;   //!
    TBranch        *b_K_1_1270_plus_ptasy_0_50;   //!
    TBranch        *b_K_1_1270_plus_cpx_0_60;   //!
    TBranch        *b_K_1_1270_plus_cpy_0_60;   //!
    TBranch        *b_K_1_1270_plus_cpz_0_60;   //!
    TBranch        *b_K_1_1270_plus_cpt_0_60;   //!
    TBranch        *b_K_1_1270_plus_cp_0_60;   //!
    TBranch        *b_K_1_1270_plus_cmult_0_60;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_0_60;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_0_60;   //!
    TBranch        *b_K_1_1270_plus_pxasy_0_60;   //!
    TBranch        *b_K_1_1270_plus_pyasy_0_60;   //!
    TBranch        *b_K_1_1270_plus_pzasy_0_60;   //!
    TBranch        *b_K_1_1270_plus_pasy_0_60;   //!
    TBranch        *b_K_1_1270_plus_ptasy_0_60;   //!
    TBranch        *b_K_1_1270_plus_cpx_0_70;   //!
    TBranch        *b_K_1_1270_plus_cpy_0_70;   //!
    TBranch        *b_K_1_1270_plus_cpz_0_70;   //!
    TBranch        *b_K_1_1270_plus_cpt_0_70;   //!
    TBranch        *b_K_1_1270_plus_cp_0_70;   //!
    TBranch        *b_K_1_1270_plus_cmult_0_70;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_0_70;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_0_70;   //!
    TBranch        *b_K_1_1270_plus_pxasy_0_70;   //!
    TBranch        *b_K_1_1270_plus_pyasy_0_70;   //!
    TBranch        *b_K_1_1270_plus_pzasy_0_70;   //!
    TBranch        *b_K_1_1270_plus_pasy_0_70;   //!
    TBranch        *b_K_1_1270_plus_ptasy_0_70;   //!
    TBranch        *b_K_1_1270_plus_cpx_0_80;   //!
    TBranch        *b_K_1_1270_plus_cpy_0_80;   //!
    TBranch        *b_K_1_1270_plus_cpz_0_80;   //!
    TBranch        *b_K_1_1270_plus_cpt_0_80;   //!
    TBranch        *b_K_1_1270_plus_cp_0_80;   //!
    TBranch        *b_K_1_1270_plus_cmult_0_80;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_0_80;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_0_80;   //!
    TBranch        *b_K_1_1270_plus_pxasy_0_80;   //!
    TBranch        *b_K_1_1270_plus_pyasy_0_80;   //!
    TBranch        *b_K_1_1270_plus_pzasy_0_80;   //!
    TBranch        *b_K_1_1270_plus_pasy_0_80;   //!
    TBranch        *b_K_1_1270_plus_ptasy_0_80;   //!
    TBranch        *b_K_1_1270_plus_cpx_0_90;   //!
    TBranch        *b_K_1_1270_plus_cpy_0_90;   //!
    TBranch        *b_K_1_1270_plus_cpz_0_90;   //!
    TBranch        *b_K_1_1270_plus_cpt_0_90;   //!
    TBranch        *b_K_1_1270_plus_cp_0_90;   //!
    TBranch        *b_K_1_1270_plus_cmult_0_90;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_0_90;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_0_90;   //!
    TBranch        *b_K_1_1270_plus_pxasy_0_90;   //!
    TBranch        *b_K_1_1270_plus_pyasy_0_90;   //!
    TBranch        *b_K_1_1270_plus_pzasy_0_90;   //!
    TBranch        *b_K_1_1270_plus_pasy_0_90;   //!
    TBranch        *b_K_1_1270_plus_ptasy_0_90;   //!
    TBranch        *b_K_1_1270_plus_cpx_1_00;   //!
    TBranch        *b_K_1_1270_plus_cpy_1_00;   //!
    TBranch        *b_K_1_1270_plus_cpz_1_00;   //!
    TBranch        *b_K_1_1270_plus_cpt_1_00;   //!
    TBranch        *b_K_1_1270_plus_cp_1_00;   //!
    TBranch        *b_K_1_1270_plus_cmult_1_00;   //!
    TBranch        *b_K_1_1270_plus_deltaEta_1_00;   //!
    TBranch        *b_K_1_1270_plus_deltaPhi_1_00;   //!
    TBranch        *b_K_1_1270_plus_pxasy_1_00;   //!
    TBranch        *b_K_1_1270_plus_pyasy_1_00;   //!
    TBranch        *b_K_1_1270_plus_pzasy_1_00;   //!
    TBranch        *b_K_1_1270_plus_pasy_1_00;   //!
    TBranch        *b_K_1_1270_plus_ptasy_1_00;   //!
    TBranch        *b_K_plus_DOCA1;   //!
    TBranch        *b_K_plus_DOCA2;   //!
    TBranch        *b_K_plus_DOCA3;   //!
    TBranch        *b_K_plus_ETA;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_K_plus_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_K_plus_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_K_plus_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_K_plus_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_K_plus_CosTheta;   //!
    TBranch        *b_K_plus_OWNPV_X;   //!
    TBranch        *b_K_plus_OWNPV_Y;   //!
    TBranch        *b_K_plus_OWNPV_Z;   //!
    TBranch        *b_K_plus_OWNPV_XERR;   //!
    TBranch        *b_K_plus_OWNPV_YERR;   //!
    TBranch        *b_K_plus_OWNPV_ZERR;   //!
    TBranch        *b_K_plus_OWNPV_CHI2;   //!
    TBranch        *b_K_plus_OWNPV_NDOF;   //!
    TBranch        *b_K_plus_OWNPV_COV_;   //!
    TBranch        *b_K_plus_IP_OWNPV;   //!
    TBranch        *b_K_plus_IPCHI2_OWNPV;   //!
    TBranch        *b_K_plus_ORIVX_X;   //!
    TBranch        *b_K_plus_ORIVX_Y;   //!
    TBranch        *b_K_plus_ORIVX_Z;   //!
    TBranch        *b_K_plus_ORIVX_XERR;   //!
    TBranch        *b_K_plus_ORIVX_YERR;   //!
    TBranch        *b_K_plus_ORIVX_ZERR;   //!
    TBranch        *b_K_plus_ORIVX_CHI2;   //!
    TBranch        *b_K_plus_ORIVX_NDOF;   //!
    TBranch        *b_K_plus_ORIVX_COV_;   //!
    TBranch        *b_K_plus_P;   //!
    TBranch        *b_K_plus_PT;   //!
    TBranch        *b_K_plus_PE;   //!
    TBranch        *b_K_plus_PX;   //!
    TBranch        *b_K_plus_PY;   //!
    TBranch        *b_K_plus_PZ;   //!
    TBranch        *b_K_plus_M;   //!
    TBranch        *b_K_plus_ID;   //!
    TBranch        *b_K_plus_PIDe;   //!
    TBranch        *b_K_plus_PIDmu;   //!
    TBranch        *b_K_plus_PIDK;   //!
    TBranch        *b_K_plus_PIDp;   //!
    TBranch        *b_K_plus_ProbNNe;   //!
    TBranch        *b_K_plus_ProbNNk;   //!
    TBranch        *b_K_plus_ProbNNp;   //!
    TBranch        *b_K_plus_ProbNNpi;   //!
    TBranch        *b_K_plus_ProbNNmu;   //!
    TBranch        *b_K_plus_ProbNNghost;   //!
    TBranch        *b_K_plus_hasMuon;   //!
    TBranch        *b_K_plus_isMuon;   //!
    TBranch        *b_K_plus_hasRich;   //!
    TBranch        *b_K_plus_UsedRichAerogel;   //!
    TBranch        *b_K_plus_UsedRich1Gas;   //!
    TBranch        *b_K_plus_UsedRich2Gas;   //!
    TBranch        *b_K_plus_RichAboveElThres;   //!
    TBranch        *b_K_plus_RichAboveMuThres;   //!
    TBranch        *b_K_plus_RichAbovePiThres;   //!
    TBranch        *b_K_plus_RichAboveKaThres;   //!
    TBranch        *b_K_plus_RichAbovePrThres;   //!
    TBranch        *b_K_plus_hasCalo;   //!
    TBranch        *b_K_plus_L0Global_Dec;   //!
    TBranch        *b_K_plus_L0Global_TIS;   //!
    TBranch        *b_K_plus_L0Global_TOS;   //!
    TBranch        *b_K_plus_Hlt1Global_Dec;   //!
    TBranch        *b_K_plus_Hlt1Global_TIS;   //!
    TBranch        *b_K_plus_Hlt1Global_TOS;   //!
    TBranch        *b_K_plus_Hlt1Phys_Dec;   //!
    TBranch        *b_K_plus_Hlt1Phys_TIS;   //!
    TBranch        *b_K_plus_Hlt1Phys_TOS;   //!
    TBranch        *b_K_plus_Hlt2Global_Dec;   //!
    TBranch        *b_K_plus_Hlt2Global_TIS;   //!
    TBranch        *b_K_plus_Hlt2Global_TOS;   //!
    TBranch        *b_K_plus_Hlt2Phys_Dec;   //!
    TBranch        *b_K_plus_Hlt2Phys_TIS;   //!
    TBranch        *b_K_plus_Hlt2Phys_TOS;   //!
    TBranch        *b_K_plus_L0HadronDecision_Dec;   //!
    TBranch        *b_K_plus_L0HadronDecision_TIS;   //!
    TBranch        *b_K_plus_L0HadronDecision_TOS;   //!
    TBranch        *b_K_plus_L0MuonDecision_Dec;   //!
    TBranch        *b_K_plus_L0MuonDecision_TIS;   //!
    TBranch        *b_K_plus_L0MuonDecision_TOS;   //!
    TBranch        *b_K_plus_L0GlobalDecision_Dec;   //!
    TBranch        *b_K_plus_L0GlobalDecision_TIS;   //!
    TBranch        *b_K_plus_L0GlobalDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_K_plus_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_K_plus_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_K_plus_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_K_plus_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_K_plus_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_K_plus_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_K_plus_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_K_plus_TRACK_Type;   //!
    TBranch        *b_K_plus_TRACK_Key;   //!
    TBranch        *b_K_plus_TRACK_CHI2NDOF;   //!
    TBranch        *b_K_plus_TRACK_PCHI2;   //!
    TBranch        *b_K_plus_TRACK_MatchCHI2;   //!
    TBranch        *b_K_plus_TRACK_GhostProb;   //!
    TBranch        *b_K_plus_TRACK_CloneDist;   //!
    TBranch        *b_K_plus_TRACK_Likelihood;   //!
    TBranch        *b_K_plus_cpx_0_50;   //!
    TBranch        *b_K_plus_cpy_0_50;   //!
    TBranch        *b_K_plus_cpz_0_50;   //!
    TBranch        *b_K_plus_cpt_0_50;   //!
    TBranch        *b_K_plus_cp_0_50;   //!
    TBranch        *b_K_plus_cmult_0_50;   //!
    TBranch        *b_K_plus_deltaEta_0_50;   //!
    TBranch        *b_K_plus_deltaPhi_0_50;   //!
    TBranch        *b_K_plus_pxasy_0_50;   //!
    TBranch        *b_K_plus_pyasy_0_50;   //!
    TBranch        *b_K_plus_pzasy_0_50;   //!
    TBranch        *b_K_plus_pasy_0_50;   //!
    TBranch        *b_K_plus_ptasy_0_50;   //!
    TBranch        *b_K_plus_cpx_0_60;   //!
    TBranch        *b_K_plus_cpy_0_60;   //!
    TBranch        *b_K_plus_cpz_0_60;   //!
    TBranch        *b_K_plus_cpt_0_60;   //!
    TBranch        *b_K_plus_cp_0_60;   //!
    TBranch        *b_K_plus_cmult_0_60;   //!
    TBranch        *b_K_plus_deltaEta_0_60;   //!
    TBranch        *b_K_plus_deltaPhi_0_60;   //!
    TBranch        *b_K_plus_pxasy_0_60;   //!
    TBranch        *b_K_plus_pyasy_0_60;   //!
    TBranch        *b_K_plus_pzasy_0_60;   //!
    TBranch        *b_K_plus_pasy_0_60;   //!
    TBranch        *b_K_plus_ptasy_0_60;   //!
    TBranch        *b_K_plus_cpx_0_70;   //!
    TBranch        *b_K_plus_cpy_0_70;   //!
    TBranch        *b_K_plus_cpz_0_70;   //!
    TBranch        *b_K_plus_cpt_0_70;   //!
    TBranch        *b_K_plus_cp_0_70;   //!
    TBranch        *b_K_plus_cmult_0_70;   //!
    TBranch        *b_K_plus_deltaEta_0_70;   //!
    TBranch        *b_K_plus_deltaPhi_0_70;   //!
    TBranch        *b_K_plus_pxasy_0_70;   //!
    TBranch        *b_K_plus_pyasy_0_70;   //!
    TBranch        *b_K_plus_pzasy_0_70;   //!
    TBranch        *b_K_plus_pasy_0_70;   //!
    TBranch        *b_K_plus_ptasy_0_70;   //!
    TBranch        *b_K_plus_cpx_0_80;   //!
    TBranch        *b_K_plus_cpy_0_80;   //!
    TBranch        *b_K_plus_cpz_0_80;   //!
    TBranch        *b_K_plus_cpt_0_80;   //!
    TBranch        *b_K_plus_cp_0_80;   //!
    TBranch        *b_K_plus_cmult_0_80;   //!
    TBranch        *b_K_plus_deltaEta_0_80;   //!
    TBranch        *b_K_plus_deltaPhi_0_80;   //!
    TBranch        *b_K_plus_pxasy_0_80;   //!
    TBranch        *b_K_plus_pyasy_0_80;   //!
    TBranch        *b_K_plus_pzasy_0_80;   //!
    TBranch        *b_K_plus_pasy_0_80;   //!
    TBranch        *b_K_plus_ptasy_0_80;   //!
    TBranch        *b_K_plus_cpx_0_90;   //!
    TBranch        *b_K_plus_cpy_0_90;   //!
    TBranch        *b_K_plus_cpz_0_90;   //!
    TBranch        *b_K_plus_cpt_0_90;   //!
    TBranch        *b_K_plus_cp_0_90;   //!
    TBranch        *b_K_plus_cmult_0_90;   //!
    TBranch        *b_K_plus_deltaEta_0_90;   //!
    TBranch        *b_K_plus_deltaPhi_0_90;   //!
    TBranch        *b_K_plus_pxasy_0_90;   //!
    TBranch        *b_K_plus_pyasy_0_90;   //!
    TBranch        *b_K_plus_pzasy_0_90;   //!
    TBranch        *b_K_plus_pasy_0_90;   //!
    TBranch        *b_K_plus_ptasy_0_90;   //!
    TBranch        *b_K_plus_cpx_1_00;   //!
    TBranch        *b_K_plus_cpy_1_00;   //!
    TBranch        *b_K_plus_cpz_1_00;   //!
    TBranch        *b_K_plus_cpt_1_00;   //!
    TBranch        *b_K_plus_cp_1_00;   //!
    TBranch        *b_K_plus_cmult_1_00;   //!
    TBranch        *b_K_plus_deltaEta_1_00;   //!
    TBranch        *b_K_plus_deltaPhi_1_00;   //!
    TBranch        *b_K_plus_pxasy_1_00;   //!
    TBranch        *b_K_plus_pyasy_1_00;   //!
    TBranch        *b_K_plus_pzasy_1_00;   //!
    TBranch        *b_K_plus_pasy_1_00;   //!
    TBranch        *b_K_plus_ptasy_1_00;   //!
    TBranch        *b_pi_plus_DOCA1;   //!
    TBranch        *b_pi_plus_DOCA2;   //!
    TBranch        *b_pi_plus_DOCA3;   //!
    TBranch        *b_pi_plus_ETA;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_pi_plus_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_pi_plus_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_pi_plus_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_pi_plus_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_pi_plus_CosTheta;   //!
    TBranch        *b_pi_plus_OWNPV_X;   //!
    TBranch        *b_pi_plus_OWNPV_Y;   //!
    TBranch        *b_pi_plus_OWNPV_Z;   //!
    TBranch        *b_pi_plus_OWNPV_XERR;   //!
    TBranch        *b_pi_plus_OWNPV_YERR;   //!
    TBranch        *b_pi_plus_OWNPV_ZERR;   //!
    TBranch        *b_pi_plus_OWNPV_CHI2;   //!
    TBranch        *b_pi_plus_OWNPV_NDOF;   //!
    TBranch        *b_pi_plus_OWNPV_COV_;   //!
    TBranch        *b_pi_plus_IP_OWNPV;   //!
    TBranch        *b_pi_plus_IPCHI2_OWNPV;   //!
    TBranch        *b_pi_plus_ORIVX_X;   //!
    TBranch        *b_pi_plus_ORIVX_Y;   //!
    TBranch        *b_pi_plus_ORIVX_Z;   //!
    TBranch        *b_pi_plus_ORIVX_XERR;   //!
    TBranch        *b_pi_plus_ORIVX_YERR;   //!
    TBranch        *b_pi_plus_ORIVX_ZERR;   //!
    TBranch        *b_pi_plus_ORIVX_CHI2;   //!
    TBranch        *b_pi_plus_ORIVX_NDOF;   //!
    TBranch        *b_pi_plus_ORIVX_COV_;   //!
    TBranch        *b_pi_plus_P;   //!
    TBranch        *b_pi_plus_PT;   //!
    TBranch        *b_pi_plus_PE;   //!
    TBranch        *b_pi_plus_PX;   //!
    TBranch        *b_pi_plus_PY;   //!
    TBranch        *b_pi_plus_PZ;   //!
    TBranch        *b_pi_plus_M;   //!
    TBranch        *b_pi_plus_ID;   //!
    TBranch        *b_pi_plus_PIDe;   //!
    TBranch        *b_pi_plus_PIDmu;   //!
    TBranch        *b_pi_plus_PIDK;   //!
    TBranch        *b_pi_plus_PIDp;   //!
    TBranch        *b_pi_plus_ProbNNe;   //!
    TBranch        *b_pi_plus_ProbNNk;   //!
    TBranch        *b_pi_plus_ProbNNp;   //!
    TBranch        *b_pi_plus_ProbNNpi;   //!
    TBranch        *b_pi_plus_ProbNNmu;   //!
    TBranch        *b_pi_plus_ProbNNghost;   //!
    TBranch        *b_pi_plus_hasMuon;   //!
    TBranch        *b_pi_plus_isMuon;   //!
    TBranch        *b_pi_plus_hasRich;   //!
    TBranch        *b_pi_plus_UsedRichAerogel;   //!
    TBranch        *b_pi_plus_UsedRich1Gas;   //!
    TBranch        *b_pi_plus_UsedRich2Gas;   //!
    TBranch        *b_pi_plus_RichAboveElThres;   //!
    TBranch        *b_pi_plus_RichAboveMuThres;   //!
    TBranch        *b_pi_plus_RichAbovePiThres;   //!
    TBranch        *b_pi_plus_RichAboveKaThres;   //!
    TBranch        *b_pi_plus_RichAbovePrThres;   //!
    TBranch        *b_pi_plus_hasCalo;   //!
    TBranch        *b_pi_plus_L0Global_Dec;   //!
    TBranch        *b_pi_plus_L0Global_TIS;   //!
    TBranch        *b_pi_plus_L0Global_TOS;   //!
    TBranch        *b_pi_plus_Hlt1Global_Dec;   //!
    TBranch        *b_pi_plus_Hlt1Global_TIS;   //!
    TBranch        *b_pi_plus_Hlt1Global_TOS;   //!
    TBranch        *b_pi_plus_Hlt1Phys_Dec;   //!
    TBranch        *b_pi_plus_Hlt1Phys_TIS;   //!
    TBranch        *b_pi_plus_Hlt1Phys_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Global_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Global_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Global_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Phys_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Phys_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Phys_TOS;   //!
    TBranch        *b_pi_plus_L0HadronDecision_Dec;   //!
    TBranch        *b_pi_plus_L0HadronDecision_TIS;   //!
    TBranch        *b_pi_plus_L0HadronDecision_TOS;   //!
    TBranch        *b_pi_plus_L0MuonDecision_Dec;   //!
    TBranch        *b_pi_plus_L0MuonDecision_TIS;   //!
    TBranch        *b_pi_plus_L0MuonDecision_TOS;   //!
    TBranch        *b_pi_plus_L0GlobalDecision_Dec;   //!
    TBranch        *b_pi_plus_L0GlobalDecision_TIS;   //!
    TBranch        *b_pi_plus_L0GlobalDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_pi_plus_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_pi_plus_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_pi_plus_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_pi_plus_TRACK_Type;   //!
    TBranch        *b_pi_plus_TRACK_Key;   //!
    TBranch        *b_pi_plus_TRACK_CHI2NDOF;   //!
    TBranch        *b_pi_plus_TRACK_PCHI2;   //!
    TBranch        *b_pi_plus_TRACK_MatchCHI2;   //!
    TBranch        *b_pi_plus_TRACK_GhostProb;   //!
    TBranch        *b_pi_plus_TRACK_CloneDist;   //!
    TBranch        *b_pi_plus_TRACK_Likelihood;   //!
    TBranch        *b_pi_plus_cpx_0_50;   //!
    TBranch        *b_pi_plus_cpy_0_50;   //!
    TBranch        *b_pi_plus_cpz_0_50;   //!
    TBranch        *b_pi_plus_cpt_0_50;   //!
    TBranch        *b_pi_plus_cp_0_50;   //!
    TBranch        *b_pi_plus_cmult_0_50;   //!
    TBranch        *b_pi_plus_deltaEta_0_50;   //!
    TBranch        *b_pi_plus_deltaPhi_0_50;   //!
    TBranch        *b_pi_plus_pxasy_0_50;   //!
    TBranch        *b_pi_plus_pyasy_0_50;   //!
    TBranch        *b_pi_plus_pzasy_0_50;   //!
    TBranch        *b_pi_plus_pasy_0_50;   //!
    TBranch        *b_pi_plus_ptasy_0_50;   //!
    TBranch        *b_pi_plus_cpx_0_60;   //!
    TBranch        *b_pi_plus_cpy_0_60;   //!
    TBranch        *b_pi_plus_cpz_0_60;   //!
    TBranch        *b_pi_plus_cpt_0_60;   //!
    TBranch        *b_pi_plus_cp_0_60;   //!
    TBranch        *b_pi_plus_cmult_0_60;   //!
    TBranch        *b_pi_plus_deltaEta_0_60;   //!
    TBranch        *b_pi_plus_deltaPhi_0_60;   //!
    TBranch        *b_pi_plus_pxasy_0_60;   //!
    TBranch        *b_pi_plus_pyasy_0_60;   //!
    TBranch        *b_pi_plus_pzasy_0_60;   //!
    TBranch        *b_pi_plus_pasy_0_60;   //!
    TBranch        *b_pi_plus_ptasy_0_60;   //!
    TBranch        *b_pi_plus_cpx_0_70;   //!
    TBranch        *b_pi_plus_cpy_0_70;   //!
    TBranch        *b_pi_plus_cpz_0_70;   //!
    TBranch        *b_pi_plus_cpt_0_70;   //!
    TBranch        *b_pi_plus_cp_0_70;   //!
    TBranch        *b_pi_plus_cmult_0_70;   //!
    TBranch        *b_pi_plus_deltaEta_0_70;   //!
    TBranch        *b_pi_plus_deltaPhi_0_70;   //!
    TBranch        *b_pi_plus_pxasy_0_70;   //!
    TBranch        *b_pi_plus_pyasy_0_70;   //!
    TBranch        *b_pi_plus_pzasy_0_70;   //!
    TBranch        *b_pi_plus_pasy_0_70;   //!
    TBranch        *b_pi_plus_ptasy_0_70;   //!
    TBranch        *b_pi_plus_cpx_0_80;   //!
    TBranch        *b_pi_plus_cpy_0_80;   //!
    TBranch        *b_pi_plus_cpz_0_80;   //!
    TBranch        *b_pi_plus_cpt_0_80;   //!
    TBranch        *b_pi_plus_cp_0_80;   //!
    TBranch        *b_pi_plus_cmult_0_80;   //!
    TBranch        *b_pi_plus_deltaEta_0_80;   //!
    TBranch        *b_pi_plus_deltaPhi_0_80;   //!
    TBranch        *b_pi_plus_pxasy_0_80;   //!
    TBranch        *b_pi_plus_pyasy_0_80;   //!
    TBranch        *b_pi_plus_pzasy_0_80;   //!
    TBranch        *b_pi_plus_pasy_0_80;   //!
    TBranch        *b_pi_plus_ptasy_0_80;   //!
    TBranch        *b_pi_plus_cpx_0_90;   //!
    TBranch        *b_pi_plus_cpy_0_90;   //!
    TBranch        *b_pi_plus_cpz_0_90;   //!
    TBranch        *b_pi_plus_cpt_0_90;   //!
    TBranch        *b_pi_plus_cp_0_90;   //!
    TBranch        *b_pi_plus_cmult_0_90;   //!
    TBranch        *b_pi_plus_deltaEta_0_90;   //!
    TBranch        *b_pi_plus_deltaPhi_0_90;   //!
    TBranch        *b_pi_plus_pxasy_0_90;   //!
    TBranch        *b_pi_plus_pyasy_0_90;   //!
    TBranch        *b_pi_plus_pzasy_0_90;   //!
    TBranch        *b_pi_plus_pasy_0_90;   //!
    TBranch        *b_pi_plus_ptasy_0_90;   //!
    TBranch        *b_pi_plus_cpx_1_00;   //!
    TBranch        *b_pi_plus_cpy_1_00;   //!
    TBranch        *b_pi_plus_cpz_1_00;   //!
    TBranch        *b_pi_plus_cpt_1_00;   //!
    TBranch        *b_pi_plus_cp_1_00;   //!
    TBranch        *b_pi_plus_cmult_1_00;   //!
    TBranch        *b_pi_plus_deltaEta_1_00;   //!
    TBranch        *b_pi_plus_deltaPhi_1_00;   //!
    TBranch        *b_pi_plus_pxasy_1_00;   //!
    TBranch        *b_pi_plus_pyasy_1_00;   //!
    TBranch        *b_pi_plus_pzasy_1_00;   //!
    TBranch        *b_pi_plus_pasy_1_00;   //!
    TBranch        *b_pi_plus_ptasy_1_00;   //!
    TBranch        *b_pi_minus_DOCA1;   //!
    TBranch        *b_pi_minus_DOCA2;   //!
    TBranch        *b_pi_minus_DOCA3;   //!
    TBranch        *b_pi_minus_ETA;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNe;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_pi_minus_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNe;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_pi_minus_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNe;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_pi_minus_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNe;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_pi_minus_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_pi_minus_CosTheta;   //!
    TBranch        *b_pi_minus_OWNPV_X;   //!
    TBranch        *b_pi_minus_OWNPV_Y;   //!
    TBranch        *b_pi_minus_OWNPV_Z;   //!
    TBranch        *b_pi_minus_OWNPV_XERR;   //!
    TBranch        *b_pi_minus_OWNPV_YERR;   //!
    TBranch        *b_pi_minus_OWNPV_ZERR;   //!
    TBranch        *b_pi_minus_OWNPV_CHI2;   //!
    TBranch        *b_pi_minus_OWNPV_NDOF;   //!
    TBranch        *b_pi_minus_OWNPV_COV_;   //!
    TBranch        *b_pi_minus_IP_OWNPV;   //!
    TBranch        *b_pi_minus_IPCHI2_OWNPV;   //!
    TBranch        *b_pi_minus_ORIVX_X;   //!
    TBranch        *b_pi_minus_ORIVX_Y;   //!
    TBranch        *b_pi_minus_ORIVX_Z;   //!
    TBranch        *b_pi_minus_ORIVX_XERR;   //!
    TBranch        *b_pi_minus_ORIVX_YERR;   //!
    TBranch        *b_pi_minus_ORIVX_ZERR;   //!
    TBranch        *b_pi_minus_ORIVX_CHI2;   //!
    TBranch        *b_pi_minus_ORIVX_NDOF;   //!
    TBranch        *b_pi_minus_ORIVX_COV_;   //!
    TBranch        *b_pi_minus_P;   //!
    TBranch        *b_pi_minus_PT;   //!
    TBranch        *b_pi_minus_PE;   //!
    TBranch        *b_pi_minus_PX;   //!
    TBranch        *b_pi_minus_PY;   //!
    TBranch        *b_pi_minus_PZ;   //!
    TBranch        *b_pi_minus_M;   //!
    TBranch        *b_pi_minus_ID;   //!
    TBranch        *b_pi_minus_PIDe;   //!
    TBranch        *b_pi_minus_PIDmu;   //!
    TBranch        *b_pi_minus_PIDK;   //!
    TBranch        *b_pi_minus_PIDp;   //!
    TBranch        *b_pi_minus_ProbNNe;   //!
    TBranch        *b_pi_minus_ProbNNk;   //!
    TBranch        *b_pi_minus_ProbNNp;   //!
    TBranch        *b_pi_minus_ProbNNpi;   //!
    TBranch        *b_pi_minus_ProbNNmu;   //!
    TBranch        *b_pi_minus_ProbNNghost;   //!
    TBranch        *b_pi_minus_hasMuon;   //!
    TBranch        *b_pi_minus_isMuon;   //!
    TBranch        *b_pi_minus_hasRich;   //!
    TBranch        *b_pi_minus_UsedRichAerogel;   //!
    TBranch        *b_pi_minus_UsedRich1Gas;   //!
    TBranch        *b_pi_minus_UsedRich2Gas;   //!
    TBranch        *b_pi_minus_RichAboveElThres;   //!
    TBranch        *b_pi_minus_RichAboveMuThres;   //!
    TBranch        *b_pi_minus_RichAbovePiThres;   //!
    TBranch        *b_pi_minus_RichAboveKaThres;   //!
    TBranch        *b_pi_minus_RichAbovePrThres;   //!
    TBranch        *b_pi_minus_hasCalo;   //!
    TBranch        *b_pi_minus_L0Global_Dec;   //!
    TBranch        *b_pi_minus_L0Global_TIS;   //!
    TBranch        *b_pi_minus_L0Global_TOS;   //!
    TBranch        *b_pi_minus_Hlt1Global_Dec;   //!
    TBranch        *b_pi_minus_Hlt1Global_TIS;   //!
    TBranch        *b_pi_minus_Hlt1Global_TOS;   //!
    TBranch        *b_pi_minus_Hlt1Phys_Dec;   //!
    TBranch        *b_pi_minus_Hlt1Phys_TIS;   //!
    TBranch        *b_pi_minus_Hlt1Phys_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Global_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Global_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Global_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Phys_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Phys_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Phys_TOS;   //!
    TBranch        *b_pi_minus_L0HadronDecision_Dec;   //!
    TBranch        *b_pi_minus_L0HadronDecision_TIS;   //!
    TBranch        *b_pi_minus_L0HadronDecision_TOS;   //!
    TBranch        *b_pi_minus_L0MuonDecision_Dec;   //!
    TBranch        *b_pi_minus_L0MuonDecision_TIS;   //!
    TBranch        *b_pi_minus_L0MuonDecision_TOS;   //!
    TBranch        *b_pi_minus_L0GlobalDecision_Dec;   //!
    TBranch        *b_pi_minus_L0GlobalDecision_TIS;   //!
    TBranch        *b_pi_minus_L0GlobalDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_pi_minus_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_pi_minus_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVADecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVADecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVADecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVADecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVADecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVADecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt1TrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVALooseDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVALooseDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt1TwoTrackMVALooseDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2IncPhiDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2IncPhiDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2IncPhiDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2PhiIncPhiDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2PhiIncPhiDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2PhiIncPhiDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo2BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo3BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyDecision_Dec;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyDecision_TIS;   //!
    TBranch        *b_pi_minus_Hlt2Topo4BodyDecision_TOS;   //!
    TBranch        *b_pi_minus_TRACK_Type;   //!
    TBranch        *b_pi_minus_TRACK_Key;   //!
    TBranch        *b_pi_minus_TRACK_CHI2NDOF;   //!
    TBranch        *b_pi_minus_TRACK_PCHI2;   //!
    TBranch        *b_pi_minus_TRACK_MatchCHI2;   //!
    TBranch        *b_pi_minus_TRACK_GhostProb;   //!
    TBranch        *b_pi_minus_TRACK_CloneDist;   //!
    TBranch        *b_pi_minus_TRACK_Likelihood;   //!
    TBranch        *b_pi_minus_cpx_0_50;   //!
    TBranch        *b_pi_minus_cpy_0_50;   //!
    TBranch        *b_pi_minus_cpz_0_50;   //!
    TBranch        *b_pi_minus_cpt_0_50;   //!
    TBranch        *b_pi_minus_cp_0_50;   //!
    TBranch        *b_pi_minus_cmult_0_50;   //!
    TBranch        *b_pi_minus_deltaEta_0_50;   //!
    TBranch        *b_pi_minus_deltaPhi_0_50;   //!
    TBranch        *b_pi_minus_pxasy_0_50;   //!
    TBranch        *b_pi_minus_pyasy_0_50;   //!
    TBranch        *b_pi_minus_pzasy_0_50;   //!
    TBranch        *b_pi_minus_pasy_0_50;   //!
    TBranch        *b_pi_minus_ptasy_0_50;   //!
    TBranch        *b_pi_minus_cpx_0_60;   //!
    TBranch        *b_pi_minus_cpy_0_60;   //!
    TBranch        *b_pi_minus_cpz_0_60;   //!
    TBranch        *b_pi_minus_cpt_0_60;   //!
    TBranch        *b_pi_minus_cp_0_60;   //!
    TBranch        *b_pi_minus_cmult_0_60;   //!
    TBranch        *b_pi_minus_deltaEta_0_60;   //!
    TBranch        *b_pi_minus_deltaPhi_0_60;   //!
    TBranch        *b_pi_minus_pxasy_0_60;   //!
    TBranch        *b_pi_minus_pyasy_0_60;   //!
    TBranch        *b_pi_minus_pzasy_0_60;   //!
    TBranch        *b_pi_minus_pasy_0_60;   //!
    TBranch        *b_pi_minus_ptasy_0_60;   //!
    TBranch        *b_pi_minus_cpx_0_70;   //!
    TBranch        *b_pi_minus_cpy_0_70;   //!
    TBranch        *b_pi_minus_cpz_0_70;   //!
    TBranch        *b_pi_minus_cpt_0_70;   //!
    TBranch        *b_pi_minus_cp_0_70;   //!
    TBranch        *b_pi_minus_cmult_0_70;   //!
    TBranch        *b_pi_minus_deltaEta_0_70;   //!
    TBranch        *b_pi_minus_deltaPhi_0_70;   //!
    TBranch        *b_pi_minus_pxasy_0_70;   //!
    TBranch        *b_pi_minus_pyasy_0_70;   //!
    TBranch        *b_pi_minus_pzasy_0_70;   //!
    TBranch        *b_pi_minus_pasy_0_70;   //!
    TBranch        *b_pi_minus_ptasy_0_70;   //!
    TBranch        *b_pi_minus_cpx_0_80;   //!
    TBranch        *b_pi_minus_cpy_0_80;   //!
    TBranch        *b_pi_minus_cpz_0_80;   //!
    TBranch        *b_pi_minus_cpt_0_80;   //!
    TBranch        *b_pi_minus_cp_0_80;   //!
    TBranch        *b_pi_minus_cmult_0_80;   //!
    TBranch        *b_pi_minus_deltaEta_0_80;   //!
    TBranch        *b_pi_minus_deltaPhi_0_80;   //!
    TBranch        *b_pi_minus_pxasy_0_80;   //!
    TBranch        *b_pi_minus_pyasy_0_80;   //!
    TBranch        *b_pi_minus_pzasy_0_80;   //!
    TBranch        *b_pi_minus_pasy_0_80;   //!
    TBranch        *b_pi_minus_ptasy_0_80;   //!
    TBranch        *b_pi_minus_cpx_0_90;   //!
    TBranch        *b_pi_minus_cpy_0_90;   //!
    TBranch        *b_pi_minus_cpz_0_90;   //!
    TBranch        *b_pi_minus_cpt_0_90;   //!
    TBranch        *b_pi_minus_cp_0_90;   //!
    TBranch        *b_pi_minus_cmult_0_90;   //!
    TBranch        *b_pi_minus_deltaEta_0_90;   //!
    TBranch        *b_pi_minus_deltaPhi_0_90;   //!
    TBranch        *b_pi_minus_pxasy_0_90;   //!
    TBranch        *b_pi_minus_pyasy_0_90;   //!
    TBranch        *b_pi_minus_pzasy_0_90;   //!
    TBranch        *b_pi_minus_pasy_0_90;   //!
    TBranch        *b_pi_minus_ptasy_0_90;   //!
    TBranch        *b_pi_minus_cpx_1_00;   //!
    TBranch        *b_pi_minus_cpy_1_00;   //!
    TBranch        *b_pi_minus_cpz_1_00;   //!
    TBranch        *b_pi_minus_cpt_1_00;   //!
    TBranch        *b_pi_minus_cp_1_00;   //!
    TBranch        *b_pi_minus_cmult_1_00;   //!
    TBranch        *b_pi_minus_deltaEta_1_00;   //!
    TBranch        *b_pi_minus_deltaPhi_1_00;   //!
    TBranch        *b_pi_minus_pxasy_1_00;   //!
    TBranch        *b_pi_minus_pyasy_1_00;   //!
    TBranch        *b_pi_minus_pzasy_1_00;   //!
    TBranch        *b_pi_minus_pasy_1_00;   //!
    TBranch        *b_pi_minus_ptasy_1_00;   //!
    TBranch        *b_nCandidate;   //!
    TBranch        *b_totCandidates;   //!
    TBranch        *b_EventInSequence;   //!
    TBranch        *b_runNumber;   //!
    TBranch        *b_eventNumber;   //!
    TBranch        *b_BCID;   //!
    TBranch        *b_BCType;   //!
    TBranch        *b_OdinTCK;   //!
    TBranch        *b_L0DUTCK;   //!
    TBranch        *b_HLT1TCK;   //!
    TBranch        *b_HLT2TCK;   //!
    TBranch        *b_GpsTime;   //!
    TBranch        *b_Polarity;   //!
    TBranch        *b_nPV;   //!
    TBranch        *b_PVX;   //!
    TBranch        *b_PVY;   //!
    TBranch        *b_PVZ;   //!
    TBranch        *b_PVXERR;   //!
    TBranch        *b_PVYERR;   //!
    TBranch        *b_PVZERR;   //!
    TBranch        *b_PVCHI2;   //!
    TBranch        *b_PVNDOF;   //!
    TBranch        *b_PVNTRACKS;   //!
    TBranch        *b_nPVs;   //!
    TBranch        *b_nTracks;   //!
    TBranch        *b_nLongTracks;   //!
    TBranch        *b_nDownstreamTracks;   //!
    TBranch        *b_nUpstreamTracks;   //!
    TBranch        *b_nVeloTracks;   //!
    TBranch        *b_nTTracks;   //!
    TBranch        *b_nBackTracks;   //!
    TBranch        *b_nRich1Hits;   //!
    TBranch        *b_nRich2Hits;   //!
    TBranch        *b_nVeloClusters;   //!
    TBranch        *b_nITClusters;   //!
    TBranch        *b_nTTClusters;   //!
    TBranch        *b_nOTClusters;   //!
    TBranch        *b_nSPDHits;   //!
    TBranch        *b_nMuonCoordsS0;   //!
    TBranch        *b_nMuonCoordsS1;   //!
    TBranch        *b_nMuonCoordsS2;   //!
    TBranch        *b_nMuonCoordsS3;   //!
    TBranch        *b_nMuonCoordsS4;   //!
    TBranch        *b_nMuonTracks;   //!
   
};

#endif

#ifdef DecayTree_cxx
DecayTree::~DecayTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DecayTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DecayTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DecayTree::Init()
{

    TTree* tree = this->GetInputTree();
    cout << "Found files, now init" << endl;

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
    
    fChain->SetBranchAddress("Bs_DOCA1", &Bs_DOCA1, &b_Bs_DOCA1);
    fChain->SetBranchAddress("Bs_DOCA2", &Bs_DOCA2, &b_Bs_DOCA2);
    fChain->SetBranchAddress("Bs_DOCA3", &Bs_DOCA3, &b_Bs_DOCA3);
    fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
    fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
    fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
    fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
    fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("Bs_ENDVERTEX_COV_", Bs_ENDVERTEX_COV_, &b_Bs_ENDVERTEX_COV_);
    fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
    fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
    fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
    fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
    fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
    fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
    fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
    fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
    fChain->SetBranchAddress("Bs_OWNPV_COV_", Bs_OWNPV_COV_, &b_Bs_OWNPV_COV_);
    fChain->SetBranchAddress("Bs_IP_OWNPV", &Bs_IP_OWNPV, &b_Bs_IP_OWNPV);
    fChain->SetBranchAddress("Bs_IPCHI2_OWNPV", &Bs_IPCHI2_OWNPV, &b_Bs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("Bs_FD_OWNPV", &Bs_FD_OWNPV, &b_Bs_FD_OWNPV);
    fChain->SetBranchAddress("Bs_FDCHI2_OWNPV", &Bs_FDCHI2_OWNPV, &b_Bs_FDCHI2_OWNPV);
    fChain->SetBranchAddress("Bs_DIRA_OWNPV", &Bs_DIRA_OWNPV, &b_Bs_DIRA_OWNPV);
    fChain->SetBranchAddress("Bs_P", &Bs_P, &b_Bs_P);
    fChain->SetBranchAddress("Bs_PT", &Bs_PT, &b_Bs_PT);
    fChain->SetBranchAddress("Bs_PE", &Bs_PE, &b_Bs_PE);
    fChain->SetBranchAddress("Bs_PX", &Bs_PX, &b_Bs_PX);
    fChain->SetBranchAddress("Bs_PY", &Bs_PY, &b_Bs_PY);
    fChain->SetBranchAddress("Bs_PZ", &Bs_PZ, &b_Bs_PZ);
    fChain->SetBranchAddress("Bs_MM", &Bs_MM, &b_Bs_MM);
    fChain->SetBranchAddress("Bs_MMERR", &Bs_MMERR, &b_Bs_MMERR);
    fChain->SetBranchAddress("Bs_M", &Bs_M, &b_Bs_M);
    fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
    fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
    fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
    fChain->SetBranchAddress("Bs_TAUCHI2", &Bs_TAUCHI2, &b_Bs_TAUCHI2);
    fChain->SetBranchAddress("Bs_L0Global_Dec", &Bs_L0Global_Dec, &b_Bs_L0Global_Dec);
    fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
    fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);
    fChain->SetBranchAddress("Bs_Hlt1Global_Dec", &Bs_Hlt1Global_Dec, &b_Bs_Hlt1Global_Dec);
    fChain->SetBranchAddress("Bs_Hlt1Global_TIS", &Bs_Hlt1Global_TIS, &b_Bs_Hlt1Global_TIS);
    fChain->SetBranchAddress("Bs_Hlt1Global_TOS", &Bs_Hlt1Global_TOS, &b_Bs_Hlt1Global_TOS);
    fChain->SetBranchAddress("Bs_Hlt1Phys_Dec", &Bs_Hlt1Phys_Dec, &b_Bs_Hlt1Phys_Dec);
    fChain->SetBranchAddress("Bs_Hlt1Phys_TIS", &Bs_Hlt1Phys_TIS, &b_Bs_Hlt1Phys_TIS);
    fChain->SetBranchAddress("Bs_Hlt1Phys_TOS", &Bs_Hlt1Phys_TOS, &b_Bs_Hlt1Phys_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Global_Dec", &Bs_Hlt2Global_Dec, &b_Bs_Hlt2Global_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Global_TIS", &Bs_Hlt2Global_TIS, &b_Bs_Hlt2Global_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Global_TOS", &Bs_Hlt2Global_TOS, &b_Bs_Hlt2Global_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Phys_Dec", &Bs_Hlt2Phys_Dec, &b_Bs_Hlt2Phys_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Phys_TIS", &Bs_Hlt2Phys_TIS, &b_Bs_Hlt2Phys_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Phys_TOS", &Bs_Hlt2Phys_TOS, &b_Bs_Hlt2Phys_TOS);
    fChain->SetBranchAddress("Bs_L0HadronDecision_Dec", &Bs_L0HadronDecision_Dec, &b_Bs_L0HadronDecision_Dec);
    fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
    fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
    fChain->SetBranchAddress("Bs_L0MuonDecision_Dec", &Bs_L0MuonDecision_Dec, &b_Bs_L0MuonDecision_Dec);
    fChain->SetBranchAddress("Bs_L0MuonDecision_TIS", &Bs_L0MuonDecision_TIS, &b_Bs_L0MuonDecision_TIS);
    fChain->SetBranchAddress("Bs_L0MuonDecision_TOS", &Bs_L0MuonDecision_TOS, &b_Bs_L0MuonDecision_TOS);
    fChain->SetBranchAddress("Bs_L0GlobalDecision_Dec", &Bs_L0GlobalDecision_Dec, &b_Bs_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("Bs_L0GlobalDecision_TIS", &Bs_L0GlobalDecision_TIS, &b_Bs_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("Bs_L0GlobalDecision_TOS", &Bs_L0GlobalDecision_TOS, &b_Bs_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_Dec", &Bs_Hlt1TrackAllL0Decision_Dec, &b_Bs_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TIS", &Bs_Hlt1TrackAllL0Decision_TIS, &b_Bs_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TOS", &Bs_Hlt1TrackAllL0Decision_TOS, &b_Bs_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_Dec", &Bs_Hlt1TrackMVADecision_Dec, &b_Bs_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TIS", &Bs_Hlt1TrackMVADecision_TIS, &b_Bs_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TOS", &Bs_Hlt1TrackMVADecision_TOS, &b_Bs_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_Dec", &Bs_Hlt1TwoTrackMVADecision_Dec, &b_Bs_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TIS", &Bs_Hlt1TwoTrackMVADecision_TIS, &b_Bs_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TOS", &Bs_Hlt1TwoTrackMVADecision_TOS, &b_Bs_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_Dec", &Bs_Hlt1TrackMVALooseDecision_Dec, &b_Bs_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TIS", &Bs_Hlt1TrackMVALooseDecision_TIS, &b_Bs_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TOS", &Bs_Hlt1TrackMVALooseDecision_TOS, &b_Bs_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_Dec", &Bs_Hlt1TwoTrackMVALooseDecision_Dec, &b_Bs_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TIS", &Bs_Hlt1TwoTrackMVALooseDecision_TIS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TOS", &Bs_Hlt1TwoTrackMVALooseDecision_TOS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_Dec", &Bs_Hlt2IncPhiDecision_Dec, &b_Bs_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TIS", &Bs_Hlt2IncPhiDecision_TIS, &b_Bs_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TOS", &Bs_Hlt2IncPhiDecision_TOS, &b_Bs_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_Dec", &Bs_Hlt2PhiIncPhiDecision_Dec, &b_Bs_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TIS", &Bs_Hlt2PhiIncPhiDecision_TIS, &b_Bs_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TOS", &Bs_Hlt2PhiIncPhiDecision_TOS, &b_Bs_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_Dec", &Bs_Hlt2Topo2BodyBBDTDecision_Dec, &b_Bs_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TIS", &Bs_Hlt2Topo2BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TOS", &Bs_Hlt2Topo2BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_Dec", &Bs_Hlt2Topo3BodyBBDTDecision_Dec, &b_Bs_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TIS", &Bs_Hlt2Topo3BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TOS", &Bs_Hlt2Topo3BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_Dec", &Bs_Hlt2Topo4BodyBBDTDecision_Dec, &b_Bs_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TIS", &Bs_Hlt2Topo4BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TOS", &Bs_Hlt2Topo4BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_Dec", &Bs_Hlt2Topo2BodyDecision_Dec, &b_Bs_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TIS", &Bs_Hlt2Topo2BodyDecision_TIS, &b_Bs_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TOS", &Bs_Hlt2Topo2BodyDecision_TOS, &b_Bs_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_Dec", &Bs_Hlt2Topo3BodyDecision_Dec, &b_Bs_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TIS", &Bs_Hlt2Topo3BodyDecision_TIS, &b_Bs_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TOS", &Bs_Hlt2Topo3BodyDecision_TOS, &b_Bs_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_Dec", &Bs_Hlt2Topo4BodyDecision_Dec, &b_Bs_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TIS", &Bs_Hlt2Topo4BodyDecision_TIS, &b_Bs_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TOS", &Bs_Hlt2Topo4BodyDecision_TOS, &b_Bs_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_TAGDECISION", &Bs_TAGDECISION, &b_Bs_TAGDECISION);
    fChain->SetBranchAddress("Bs_TAGOMEGA", &Bs_TAGOMEGA, &b_Bs_TAGOMEGA);
    fChain->SetBranchAddress("Bs_TAGDECISION_OS", &Bs_TAGDECISION_OS, &b_Bs_TAGDECISION_OS);
    fChain->SetBranchAddress("Bs_TAGOMEGA_OS", &Bs_TAGOMEGA_OS, &b_Bs_TAGOMEGA_OS);
    fChain->SetBranchAddress("Bs_TAGGER", &Bs_TAGGER, &b_Bs_TAGGER);
    fChain->SetBranchAddress("Bs_OS_Muon_DEC", &Bs_OS_Muon_DEC, &b_Bs_OS_Muon_DEC);
    fChain->SetBranchAddress("Bs_OS_Muon_PROB", &Bs_OS_Muon_PROB, &b_Bs_OS_Muon_PROB);
    fChain->SetBranchAddress("Bs_OS_Electron_DEC", &Bs_OS_Electron_DEC, &b_Bs_OS_Electron_DEC);
    fChain->SetBranchAddress("Bs_OS_Electron_PROB", &Bs_OS_Electron_PROB, &b_Bs_OS_Electron_PROB);
    fChain->SetBranchAddress("Bs_OS_Kaon_DEC", &Bs_OS_Kaon_DEC, &b_Bs_OS_Kaon_DEC);
    fChain->SetBranchAddress("Bs_OS_Kaon_PROB", &Bs_OS_Kaon_PROB, &b_Bs_OS_Kaon_PROB);
    fChain->SetBranchAddress("Bs_SS_Kaon_DEC", &Bs_SS_Kaon_DEC, &b_Bs_SS_Kaon_DEC);
    fChain->SetBranchAddress("Bs_SS_Kaon_PROB", &Bs_SS_Kaon_PROB, &b_Bs_SS_Kaon_PROB);
    fChain->SetBranchAddress("Bs_SS_Pion_DEC", &Bs_SS_Pion_DEC, &b_Bs_SS_Pion_DEC);
    fChain->SetBranchAddress("Bs_SS_Pion_PROB", &Bs_SS_Pion_PROB, &b_Bs_SS_Pion_PROB);
    fChain->SetBranchAddress("Bs_SS_PionBDT_DEC", &Bs_SS_PionBDT_DEC, &b_Bs_SS_PionBDT_DEC);
    fChain->SetBranchAddress("Bs_SS_PionBDT_PROB", &Bs_SS_PionBDT_PROB, &b_Bs_SS_PionBDT_PROB);
    fChain->SetBranchAddress("Bs_VtxCharge_DEC", &Bs_VtxCharge_DEC, &b_Bs_VtxCharge_DEC);
    fChain->SetBranchAddress("Bs_VtxCharge_PROB", &Bs_VtxCharge_PROB, &b_Bs_VtxCharge_PROB);
    fChain->SetBranchAddress("Bs_OS_nnetKaon_DEC", &Bs_OS_nnetKaon_DEC, &b_Bs_OS_nnetKaon_DEC);
    fChain->SetBranchAddress("Bs_OS_nnetKaon_PROB", &Bs_OS_nnetKaon_PROB, &b_Bs_OS_nnetKaon_PROB);
    fChain->SetBranchAddress("Bs_SS_nnetKaon_DEC", &Bs_SS_nnetKaon_DEC, &b_Bs_SS_nnetKaon_DEC);
    fChain->SetBranchAddress("Bs_SS_nnetKaon_PROB", &Bs_SS_nnetKaon_PROB, &b_Bs_SS_nnetKaon_PROB);
    fChain->SetBranchAddress("Bs_SS_Proton_DEC", &Bs_SS_Proton_DEC, &b_Bs_SS_Proton_DEC);
    fChain->SetBranchAddress("Bs_SS_Proton_PROB", &Bs_SS_Proton_PROB, &b_Bs_SS_Proton_PROB);
    fChain->SetBranchAddress("Bs_OS_Charm_DEC", &Bs_OS_Charm_DEC, &b_Bs_OS_Charm_DEC);
    fChain->SetBranchAddress("Bs_OS_Charm_PROB", &Bs_OS_Charm_PROB, &b_Bs_OS_Charm_PROB);
    fChain->SetBranchAddress("Bs_cpx_0.50", &Bs_cpx_0_50, &b_Bs_cpx_0_50);
    fChain->SetBranchAddress("Bs_cpy_0.50", &Bs_cpy_0_50, &b_Bs_cpy_0_50);
    fChain->SetBranchAddress("Bs_cpz_0.50", &Bs_cpz_0_50, &b_Bs_cpz_0_50);
    fChain->SetBranchAddress("Bs_cpt_0.50", &Bs_cpt_0_50, &b_Bs_cpt_0_50);
    fChain->SetBranchAddress("Bs_cp_0.50", &Bs_cp_0_50, &b_Bs_cp_0_50);
    fChain->SetBranchAddress("Bs_cmult_0.50", &Bs_cmult_0_50, &b_Bs_cmult_0_50);
    fChain->SetBranchAddress("Bs_deltaEta_0.50", &Bs_deltaEta_0_50, &b_Bs_deltaEta_0_50);
    fChain->SetBranchAddress("Bs_deltaPhi_0.50", &Bs_deltaPhi_0_50, &b_Bs_deltaPhi_0_50);
    fChain->SetBranchAddress("Bs_pxasy_0.50", &Bs_pxasy_0_50, &b_Bs_pxasy_0_50);
    fChain->SetBranchAddress("Bs_pyasy_0.50", &Bs_pyasy_0_50, &b_Bs_pyasy_0_50);
    fChain->SetBranchAddress("Bs_pzasy_0.50", &Bs_pzasy_0_50, &b_Bs_pzasy_0_50);
    fChain->SetBranchAddress("Bs_pasy_0.50", &Bs_pasy_0_50, &b_Bs_pasy_0_50);
    fChain->SetBranchAddress("Bs_ptasy_0.50", &Bs_ptasy_0_50, &b_Bs_ptasy_0_50);
    fChain->SetBranchAddress("Bs_cpx_0.60", &Bs_cpx_0_60, &b_Bs_cpx_0_60);
    fChain->SetBranchAddress("Bs_cpy_0.60", &Bs_cpy_0_60, &b_Bs_cpy_0_60);
    fChain->SetBranchAddress("Bs_cpz_0.60", &Bs_cpz_0_60, &b_Bs_cpz_0_60);
    fChain->SetBranchAddress("Bs_cpt_0.60", &Bs_cpt_0_60, &b_Bs_cpt_0_60);
    fChain->SetBranchAddress("Bs_cp_0.60", &Bs_cp_0_60, &b_Bs_cp_0_60);
    fChain->SetBranchAddress("Bs_cmult_0.60", &Bs_cmult_0_60, &b_Bs_cmult_0_60);
    fChain->SetBranchAddress("Bs_deltaEta_0.60", &Bs_deltaEta_0_60, &b_Bs_deltaEta_0_60);
    fChain->SetBranchAddress("Bs_deltaPhi_0.60", &Bs_deltaPhi_0_60, &b_Bs_deltaPhi_0_60);
    fChain->SetBranchAddress("Bs_pxasy_0.60", &Bs_pxasy_0_60, &b_Bs_pxasy_0_60);
    fChain->SetBranchAddress("Bs_pyasy_0.60", &Bs_pyasy_0_60, &b_Bs_pyasy_0_60);
    fChain->SetBranchAddress("Bs_pzasy_0.60", &Bs_pzasy_0_60, &b_Bs_pzasy_0_60);
    fChain->SetBranchAddress("Bs_pasy_0.60", &Bs_pasy_0_60, &b_Bs_pasy_0_60);
    fChain->SetBranchAddress("Bs_ptasy_0.60", &Bs_ptasy_0_60, &b_Bs_ptasy_0_60);
    fChain->SetBranchAddress("Bs_cpx_0.70", &Bs_cpx_0_70, &b_Bs_cpx_0_70);
    fChain->SetBranchAddress("Bs_cpy_0.70", &Bs_cpy_0_70, &b_Bs_cpy_0_70);
    fChain->SetBranchAddress("Bs_cpz_0.70", &Bs_cpz_0_70, &b_Bs_cpz_0_70);
    fChain->SetBranchAddress("Bs_cpt_0.70", &Bs_cpt_0_70, &b_Bs_cpt_0_70);
    fChain->SetBranchAddress("Bs_cp_0.70", &Bs_cp_0_70, &b_Bs_cp_0_70);
    fChain->SetBranchAddress("Bs_cmult_0.70", &Bs_cmult_0_70, &b_Bs_cmult_0_70);
    fChain->SetBranchAddress("Bs_deltaEta_0.70", &Bs_deltaEta_0_70, &b_Bs_deltaEta_0_70);
    fChain->SetBranchAddress("Bs_deltaPhi_0.70", &Bs_deltaPhi_0_70, &b_Bs_deltaPhi_0_70);
    fChain->SetBranchAddress("Bs_pxasy_0.70", &Bs_pxasy_0_70, &b_Bs_pxasy_0_70);
    fChain->SetBranchAddress("Bs_pyasy_0.70", &Bs_pyasy_0_70, &b_Bs_pyasy_0_70);
    fChain->SetBranchAddress("Bs_pzasy_0.70", &Bs_pzasy_0_70, &b_Bs_pzasy_0_70);
    fChain->SetBranchAddress("Bs_pasy_0.70", &Bs_pasy_0_70, &b_Bs_pasy_0_70);
    fChain->SetBranchAddress("Bs_ptasy_0.70", &Bs_ptasy_0_70, &b_Bs_ptasy_0_70);
    fChain->SetBranchAddress("Bs_cpx_0.80", &Bs_cpx_0_80, &b_Bs_cpx_0_80);
    fChain->SetBranchAddress("Bs_cpy_0.80", &Bs_cpy_0_80, &b_Bs_cpy_0_80);
    fChain->SetBranchAddress("Bs_cpz_0.80", &Bs_cpz_0_80, &b_Bs_cpz_0_80);
    fChain->SetBranchAddress("Bs_cpt_0.80", &Bs_cpt_0_80, &b_Bs_cpt_0_80);
    fChain->SetBranchAddress("Bs_cp_0.80", &Bs_cp_0_80, &b_Bs_cp_0_80);
    fChain->SetBranchAddress("Bs_cmult_0.80", &Bs_cmult_0_80, &b_Bs_cmult_0_80);
    fChain->SetBranchAddress("Bs_deltaEta_0.80", &Bs_deltaEta_0_80, &b_Bs_deltaEta_0_80);
    fChain->SetBranchAddress("Bs_deltaPhi_0.80", &Bs_deltaPhi_0_80, &b_Bs_deltaPhi_0_80);
    fChain->SetBranchAddress("Bs_pxasy_0.80", &Bs_pxasy_0_80, &b_Bs_pxasy_0_80);
    fChain->SetBranchAddress("Bs_pyasy_0.80", &Bs_pyasy_0_80, &b_Bs_pyasy_0_80);
    fChain->SetBranchAddress("Bs_pzasy_0.80", &Bs_pzasy_0_80, &b_Bs_pzasy_0_80);
    fChain->SetBranchAddress("Bs_pasy_0.80", &Bs_pasy_0_80, &b_Bs_pasy_0_80);
    fChain->SetBranchAddress("Bs_ptasy_0.80", &Bs_ptasy_0_80, &b_Bs_ptasy_0_80);
    fChain->SetBranchAddress("Bs_cpx_0.90", &Bs_cpx_0_90, &b_Bs_cpx_0_90);
    fChain->SetBranchAddress("Bs_cpy_0.90", &Bs_cpy_0_90, &b_Bs_cpy_0_90);
    fChain->SetBranchAddress("Bs_cpz_0.90", &Bs_cpz_0_90, &b_Bs_cpz_0_90);
    fChain->SetBranchAddress("Bs_cpt_0.90", &Bs_cpt_0_90, &b_Bs_cpt_0_90);
    fChain->SetBranchAddress("Bs_cp_0.90", &Bs_cp_0_90, &b_Bs_cp_0_90);
    fChain->SetBranchAddress("Bs_cmult_0.90", &Bs_cmult_0_90, &b_Bs_cmult_0_90);
    fChain->SetBranchAddress("Bs_deltaEta_0.90", &Bs_deltaEta_0_90, &b_Bs_deltaEta_0_90);
    fChain->SetBranchAddress("Bs_deltaPhi_0.90", &Bs_deltaPhi_0_90, &b_Bs_deltaPhi_0_90);
    fChain->SetBranchAddress("Bs_pxasy_0.90", &Bs_pxasy_0_90, &b_Bs_pxasy_0_90);
    fChain->SetBranchAddress("Bs_pyasy_0.90", &Bs_pyasy_0_90, &b_Bs_pyasy_0_90);
    fChain->SetBranchAddress("Bs_pzasy_0.90", &Bs_pzasy_0_90, &b_Bs_pzasy_0_90);
    fChain->SetBranchAddress("Bs_pasy_0.90", &Bs_pasy_0_90, &b_Bs_pasy_0_90);
    fChain->SetBranchAddress("Bs_ptasy_0.90", &Bs_ptasy_0_90, &b_Bs_ptasy_0_90);
    fChain->SetBranchAddress("Bs_cpx_1.00", &Bs_cpx_1_00, &b_Bs_cpx_1_00);
    fChain->SetBranchAddress("Bs_cpy_1.00", &Bs_cpy_1_00, &b_Bs_cpy_1_00);
    fChain->SetBranchAddress("Bs_cpz_1.00", &Bs_cpz_1_00, &b_Bs_cpz_1_00);
    fChain->SetBranchAddress("Bs_cpt_1.00", &Bs_cpt_1_00, &b_Bs_cpt_1_00);
    fChain->SetBranchAddress("Bs_cp_1.00", &Bs_cp_1_00, &b_Bs_cp_1_00);
    fChain->SetBranchAddress("Bs_cmult_1.00", &Bs_cmult_1_00, &b_Bs_cmult_1_00);
    fChain->SetBranchAddress("Bs_deltaEta_1.00", &Bs_deltaEta_1_00, &b_Bs_deltaEta_1_00);
    fChain->SetBranchAddress("Bs_deltaPhi_1.00", &Bs_deltaPhi_1_00, &b_Bs_deltaPhi_1_00);
    fChain->SetBranchAddress("Bs_pxasy_1.00", &Bs_pxasy_1_00, &b_Bs_pxasy_1_00);
    fChain->SetBranchAddress("Bs_pyasy_1.00", &Bs_pyasy_1_00, &b_Bs_pyasy_1_00);
    fChain->SetBranchAddress("Bs_pzasy_1.00", &Bs_pzasy_1_00, &b_Bs_pzasy_1_00);
    fChain->SetBranchAddress("Bs_pasy_1.00", &Bs_pasy_1_00, &b_Bs_pasy_1_00);
    fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
    fChain->SetBranchAddress("Bs_B0DTF_nPV", &Bs_B0DTF_nPV, &b_Bs_B0DTF_nPV);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_0_ID", Bs_B0DTF_D_splus_Kplus_0_ID, &b_Bs_B0DTF_D_splus_Kplus_0_ID);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_0_PE", Bs_B0DTF_D_splus_Kplus_0_PE, &b_Bs_B0DTF_D_splus_Kplus_0_PE);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_0_PX", Bs_B0DTF_D_splus_Kplus_0_PX, &b_Bs_B0DTF_D_splus_Kplus_0_PX);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_0_PY", Bs_B0DTF_D_splus_Kplus_0_PY, &b_Bs_B0DTF_D_splus_Kplus_0_PY);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_0_PZ", Bs_B0DTF_D_splus_Kplus_0_PZ, &b_Bs_B0DTF_D_splus_Kplus_0_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_ID", Bs_B0DTF_D_splus_Kplus_ID, &b_Bs_B0DTF_D_splus_Kplus_ID);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_PE", Bs_B0DTF_D_splus_Kplus_PE, &b_Bs_B0DTF_D_splus_Kplus_PE);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_PX", Bs_B0DTF_D_splus_Kplus_PX, &b_Bs_B0DTF_D_splus_Kplus_PX);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_PY", Bs_B0DTF_D_splus_Kplus_PY, &b_Bs_B0DTF_D_splus_Kplus_PY);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_Kplus_PZ", Bs_B0DTF_D_splus_Kplus_PZ, &b_Bs_B0DTF_D_splus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_M", Bs_B0DTF_D_splus_M, &b_Bs_B0DTF_D_splus_M);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_MERR", Bs_B0DTF_D_splus_MERR, &b_Bs_B0DTF_D_splus_MERR);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_P", Bs_B0DTF_D_splus_P, &b_Bs_B0DTF_D_splus_P);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_PERR", Bs_B0DTF_D_splus_PERR, &b_Bs_B0DTF_D_splus_PERR);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctau", Bs_B0DTF_D_splus_ctau, &b_Bs_B0DTF_D_splus_ctau);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctauErr", Bs_B0DTF_D_splus_ctauErr, &b_Bs_B0DTF_D_splus_ctauErr);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLength", Bs_B0DTF_D_splus_decayLength, &b_Bs_B0DTF_D_splus_decayLength);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLengthErr", Bs_B0DTF_D_splus_decayLengthErr, &b_Bs_B0DTF_D_splus_decayLengthErr);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_ID", Bs_B0DTF_D_splus_piplus_ID, &b_Bs_B0DTF_D_splus_piplus_ID);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PE", Bs_B0DTF_D_splus_piplus_PE, &b_Bs_B0DTF_D_splus_piplus_PE);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PX", Bs_B0DTF_D_splus_piplus_PX, &b_Bs_B0DTF_D_splus_piplus_PX);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PY", Bs_B0DTF_D_splus_piplus_PY, &b_Bs_B0DTF_D_splus_piplus_PY);
    fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PZ", Bs_B0DTF_D_splus_piplus_PZ, &b_Bs_B0DTF_D_splus_piplus_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_Kplus_ID", Bs_B0DTF_K_1_1270_plus_Kplus_ID, &b_Bs_B0DTF_K_1_1270_plus_Kplus_ID);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_Kplus_PE", Bs_B0DTF_K_1_1270_plus_Kplus_PE, &b_Bs_B0DTF_K_1_1270_plus_Kplus_PE);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_Kplus_PX", Bs_B0DTF_K_1_1270_plus_Kplus_PX, &b_Bs_B0DTF_K_1_1270_plus_Kplus_PX);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_Kplus_PY", Bs_B0DTF_K_1_1270_plus_Kplus_PY, &b_Bs_B0DTF_K_1_1270_plus_Kplus_PY);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_Kplus_PZ", Bs_B0DTF_K_1_1270_plus_Kplus_PZ, &b_Bs_B0DTF_K_1_1270_plus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_M", Bs_B0DTF_K_1_1270_plus_M, &b_Bs_B0DTF_K_1_1270_plus_M);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_MERR", Bs_B0DTF_K_1_1270_plus_MERR, &b_Bs_B0DTF_K_1_1270_plus_MERR);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_P", Bs_B0DTF_K_1_1270_plus_P, &b_Bs_B0DTF_K_1_1270_plus_P);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_PERR", Bs_B0DTF_K_1_1270_plus_PERR, &b_Bs_B0DTF_K_1_1270_plus_PERR);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_ctau", Bs_B0DTF_K_1_1270_plus_ctau, &b_Bs_B0DTF_K_1_1270_plus_ctau);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_ctauErr", Bs_B0DTF_K_1_1270_plus_ctauErr, &b_Bs_B0DTF_K_1_1270_plus_ctauErr);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_decayLength", Bs_B0DTF_K_1_1270_plus_decayLength, &b_Bs_B0DTF_K_1_1270_plus_decayLength);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_decayLengthErr", Bs_B0DTF_K_1_1270_plus_decayLengthErr, &b_Bs_B0DTF_K_1_1270_plus_decayLengthErr);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_0_ID", Bs_B0DTF_K_1_1270_plus_piplus_0_ID, &b_Bs_B0DTF_K_1_1270_plus_piplus_0_ID);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_0_PE", Bs_B0DTF_K_1_1270_plus_piplus_0_PE, &b_Bs_B0DTF_K_1_1270_plus_piplus_0_PE);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_0_PX", Bs_B0DTF_K_1_1270_plus_piplus_0_PX, &b_Bs_B0DTF_K_1_1270_plus_piplus_0_PX);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_0_PY", Bs_B0DTF_K_1_1270_plus_piplus_0_PY, &b_Bs_B0DTF_K_1_1270_plus_piplus_0_PY);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_0_PZ", Bs_B0DTF_K_1_1270_plus_piplus_0_PZ, &b_Bs_B0DTF_K_1_1270_plus_piplus_0_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_ID", Bs_B0DTF_K_1_1270_plus_piplus_ID, &b_Bs_B0DTF_K_1_1270_plus_piplus_ID);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_PE", Bs_B0DTF_K_1_1270_plus_piplus_PE, &b_Bs_B0DTF_K_1_1270_plus_piplus_PE);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_PX", Bs_B0DTF_K_1_1270_plus_piplus_PX, &b_Bs_B0DTF_K_1_1270_plus_piplus_PX);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_PY", Bs_B0DTF_K_1_1270_plus_piplus_PY, &b_Bs_B0DTF_K_1_1270_plus_piplus_PY);
    fChain->SetBranchAddress("Bs_B0DTF_K_1_1270_plus_piplus_PZ", Bs_B0DTF_K_1_1270_plus_piplus_PZ, &b_Bs_B0DTF_K_1_1270_plus_piplus_PZ);
    fChain->SetBranchAddress("Bs_B0DTF_M", Bs_B0DTF_M, &b_Bs_B0DTF_M);
    fChain->SetBranchAddress("Bs_B0DTF_MERR", Bs_B0DTF_MERR, &b_Bs_B0DTF_MERR);
    fChain->SetBranchAddress("Bs_B0DTF_P", Bs_B0DTF_P, &b_Bs_B0DTF_P);
    fChain->SetBranchAddress("Bs_B0DTF_PERR", Bs_B0DTF_PERR, &b_Bs_B0DTF_PERR);
    fChain->SetBranchAddress("Bs_B0DTF_PV_X", Bs_B0DTF_PV_X, &b_Bs_B0DTF_PV_X);
    fChain->SetBranchAddress("Bs_B0DTF_PV_Y", Bs_B0DTF_PV_Y, &b_Bs_B0DTF_PV_Y);
    fChain->SetBranchAddress("Bs_B0DTF_PV_Z", Bs_B0DTF_PV_Z, &b_Bs_B0DTF_PV_Z);
    fChain->SetBranchAddress("Bs_B0DTF_PV_key", Bs_B0DTF_PV_key, &b_Bs_B0DTF_PV_key);
    fChain->SetBranchAddress("Bs_B0DTF_chi2", Bs_B0DTF_chi2, &b_Bs_B0DTF_chi2);
    fChain->SetBranchAddress("Bs_B0DTF_ctau", Bs_B0DTF_ctau, &b_Bs_B0DTF_ctau);
    fChain->SetBranchAddress("Bs_B0DTF_ctauErr", Bs_B0DTF_ctauErr, &b_Bs_B0DTF_ctauErr);
    fChain->SetBranchAddress("Bs_B0DTF_decayLength", Bs_B0DTF_decayLength, &b_Bs_B0DTF_decayLength);
    fChain->SetBranchAddress("Bs_B0DTF_decayLengthErr", Bs_B0DTF_decayLengthErr, &b_Bs_B0DTF_decayLengthErr);
    fChain->SetBranchAddress("Bs_B0DTF_nDOF", Bs_B0DTF_nDOF, &b_Bs_B0DTF_nDOF);
    fChain->SetBranchAddress("Bs_B0DTF_nIter", Bs_B0DTF_nIter, &b_Bs_B0DTF_nIter);
    fChain->SetBranchAddress("Bs_B0DTF_status", Bs_B0DTF_status, &b_Bs_B0DTF_status);
    fChain->SetBranchAddress("Bs_BsDTF_nPV", &Bs_BsDTF_nPV, &b_Bs_BsDTF_nPV);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_0_ID", Bs_BsDTF_D_splus_Kplus_0_ID, &b_Bs_BsDTF_D_splus_Kplus_0_ID);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_0_PE", Bs_BsDTF_D_splus_Kplus_0_PE, &b_Bs_BsDTF_D_splus_Kplus_0_PE);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_0_PX", Bs_BsDTF_D_splus_Kplus_0_PX, &b_Bs_BsDTF_D_splus_Kplus_0_PX);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_0_PY", Bs_BsDTF_D_splus_Kplus_0_PY, &b_Bs_BsDTF_D_splus_Kplus_0_PY);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_0_PZ", Bs_BsDTF_D_splus_Kplus_0_PZ, &b_Bs_BsDTF_D_splus_Kplus_0_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_ID", Bs_BsDTF_D_splus_Kplus_ID, &b_Bs_BsDTF_D_splus_Kplus_ID);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_PE", Bs_BsDTF_D_splus_Kplus_PE, &b_Bs_BsDTF_D_splus_Kplus_PE);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_PX", Bs_BsDTF_D_splus_Kplus_PX, &b_Bs_BsDTF_D_splus_Kplus_PX);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_PY", Bs_BsDTF_D_splus_Kplus_PY, &b_Bs_BsDTF_D_splus_Kplus_PY);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_Kplus_PZ", Bs_BsDTF_D_splus_Kplus_PZ, &b_Bs_BsDTF_D_splus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_M", Bs_BsDTF_D_splus_M, &b_Bs_BsDTF_D_splus_M);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_MERR", Bs_BsDTF_D_splus_MERR, &b_Bs_BsDTF_D_splus_MERR);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_P", Bs_BsDTF_D_splus_P, &b_Bs_BsDTF_D_splus_P);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_PERR", Bs_BsDTF_D_splus_PERR, &b_Bs_BsDTF_D_splus_PERR);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctau", Bs_BsDTF_D_splus_ctau, &b_Bs_BsDTF_D_splus_ctau);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctauErr", Bs_BsDTF_D_splus_ctauErr, &b_Bs_BsDTF_D_splus_ctauErr);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLength", Bs_BsDTF_D_splus_decayLength, &b_Bs_BsDTF_D_splus_decayLength);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLengthErr", Bs_BsDTF_D_splus_decayLengthErr, &b_Bs_BsDTF_D_splus_decayLengthErr);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_ID", Bs_BsDTF_D_splus_piplus_ID, &b_Bs_BsDTF_D_splus_piplus_ID);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PE", Bs_BsDTF_D_splus_piplus_PE, &b_Bs_BsDTF_D_splus_piplus_PE);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PX", Bs_BsDTF_D_splus_piplus_PX, &b_Bs_BsDTF_D_splus_piplus_PX);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PY", Bs_BsDTF_D_splus_piplus_PY, &b_Bs_BsDTF_D_splus_piplus_PY);
    fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PZ", Bs_BsDTF_D_splus_piplus_PZ, &b_Bs_BsDTF_D_splus_piplus_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_Kplus_ID", Bs_BsDTF_K_1_1270_plus_Kplus_ID, &b_Bs_BsDTF_K_1_1270_plus_Kplus_ID);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_Kplus_PE", Bs_BsDTF_K_1_1270_plus_Kplus_PE, &b_Bs_BsDTF_K_1_1270_plus_Kplus_PE);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_Kplus_PX", Bs_BsDTF_K_1_1270_plus_Kplus_PX, &b_Bs_BsDTF_K_1_1270_plus_Kplus_PX);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_Kplus_PY", Bs_BsDTF_K_1_1270_plus_Kplus_PY, &b_Bs_BsDTF_K_1_1270_plus_Kplus_PY);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_Kplus_PZ", Bs_BsDTF_K_1_1270_plus_Kplus_PZ, &b_Bs_BsDTF_K_1_1270_plus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_M", Bs_BsDTF_K_1_1270_plus_M, &b_Bs_BsDTF_K_1_1270_plus_M);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_MERR", Bs_BsDTF_K_1_1270_plus_MERR, &b_Bs_BsDTF_K_1_1270_plus_MERR);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_P", Bs_BsDTF_K_1_1270_plus_P, &b_Bs_BsDTF_K_1_1270_plus_P);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_PERR", Bs_BsDTF_K_1_1270_plus_PERR, &b_Bs_BsDTF_K_1_1270_plus_PERR);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_ctau", Bs_BsDTF_K_1_1270_plus_ctau, &b_Bs_BsDTF_K_1_1270_plus_ctau);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_ctauErr", Bs_BsDTF_K_1_1270_plus_ctauErr, &b_Bs_BsDTF_K_1_1270_plus_ctauErr);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_decayLength", Bs_BsDTF_K_1_1270_plus_decayLength, &b_Bs_BsDTF_K_1_1270_plus_decayLength);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_decayLengthErr", Bs_BsDTF_K_1_1270_plus_decayLengthErr, &b_Bs_BsDTF_K_1_1270_plus_decayLengthErr);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_0_ID", Bs_BsDTF_K_1_1270_plus_piplus_0_ID, &b_Bs_BsDTF_K_1_1270_plus_piplus_0_ID);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_0_PE", Bs_BsDTF_K_1_1270_plus_piplus_0_PE, &b_Bs_BsDTF_K_1_1270_plus_piplus_0_PE);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_0_PX", Bs_BsDTF_K_1_1270_plus_piplus_0_PX, &b_Bs_BsDTF_K_1_1270_plus_piplus_0_PX);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_0_PY", Bs_BsDTF_K_1_1270_plus_piplus_0_PY, &b_Bs_BsDTF_K_1_1270_plus_piplus_0_PY);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_0_PZ", Bs_BsDTF_K_1_1270_plus_piplus_0_PZ, &b_Bs_BsDTF_K_1_1270_plus_piplus_0_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_ID", Bs_BsDTF_K_1_1270_plus_piplus_ID, &b_Bs_BsDTF_K_1_1270_plus_piplus_ID);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_PE", Bs_BsDTF_K_1_1270_plus_piplus_PE, &b_Bs_BsDTF_K_1_1270_plus_piplus_PE);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_PX", Bs_BsDTF_K_1_1270_plus_piplus_PX, &b_Bs_BsDTF_K_1_1270_plus_piplus_PX);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_PY", Bs_BsDTF_K_1_1270_plus_piplus_PY, &b_Bs_BsDTF_K_1_1270_plus_piplus_PY);
    fChain->SetBranchAddress("Bs_BsDTF_K_1_1270_plus_piplus_PZ", Bs_BsDTF_K_1_1270_plus_piplus_PZ, &b_Bs_BsDTF_K_1_1270_plus_piplus_PZ);
    fChain->SetBranchAddress("Bs_BsDTF_M", Bs_BsDTF_M, &b_Bs_BsDTF_M);
    fChain->SetBranchAddress("Bs_BsDTF_MERR", Bs_BsDTF_MERR, &b_Bs_BsDTF_MERR);
    fChain->SetBranchAddress("Bs_BsDTF_P", Bs_BsDTF_P, &b_Bs_BsDTF_P);
    fChain->SetBranchAddress("Bs_BsDTF_PERR", Bs_BsDTF_PERR, &b_Bs_BsDTF_PERR);
    fChain->SetBranchAddress("Bs_BsDTF_PV_X", Bs_BsDTF_PV_X, &b_Bs_BsDTF_PV_X);
    fChain->SetBranchAddress("Bs_BsDTF_PV_Y", Bs_BsDTF_PV_Y, &b_Bs_BsDTF_PV_Y);
    fChain->SetBranchAddress("Bs_BsDTF_PV_Z", Bs_BsDTF_PV_Z, &b_Bs_BsDTF_PV_Z);
    fChain->SetBranchAddress("Bs_BsDTF_PV_key", Bs_BsDTF_PV_key, &b_Bs_BsDTF_PV_key);
    fChain->SetBranchAddress("Bs_BsDTF_chi2", Bs_BsDTF_chi2, &b_Bs_BsDTF_chi2);
    fChain->SetBranchAddress("Bs_BsDTF_ctau", Bs_BsDTF_ctau, &b_Bs_BsDTF_ctau);
    fChain->SetBranchAddress("Bs_BsDTF_ctauErr", Bs_BsDTF_ctauErr, &b_Bs_BsDTF_ctauErr);
    fChain->SetBranchAddress("Bs_BsDTF_decayLength", Bs_BsDTF_decayLength, &b_Bs_BsDTF_decayLength);
    fChain->SetBranchAddress("Bs_BsDTF_decayLengthErr", Bs_BsDTF_decayLengthErr, &b_Bs_BsDTF_decayLengthErr);
    fChain->SetBranchAddress("Bs_BsDTF_nDOF", Bs_BsDTF_nDOF, &b_Bs_BsDTF_nDOF);
    fChain->SetBranchAddress("Bs_BsDTF_nIter", Bs_BsDTF_nIter, &b_Bs_BsDTF_nIter);
    fChain->SetBranchAddress("Bs_BsDTF_status", Bs_BsDTF_status, &b_Bs_BsDTF_status);
    fChain->SetBranchAddress("Bs_DTF_nPV", &Bs_DTF_nPV, &b_Bs_DTF_nPV);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_0_ID", Bs_DTF_D_splus_Kplus_0_ID, &b_Bs_DTF_D_splus_Kplus_0_ID);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_0_PE", Bs_DTF_D_splus_Kplus_0_PE, &b_Bs_DTF_D_splus_Kplus_0_PE);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_0_PX", Bs_DTF_D_splus_Kplus_0_PX, &b_Bs_DTF_D_splus_Kplus_0_PX);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_0_PY", Bs_DTF_D_splus_Kplus_0_PY, &b_Bs_DTF_D_splus_Kplus_0_PY);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_0_PZ", Bs_DTF_D_splus_Kplus_0_PZ, &b_Bs_DTF_D_splus_Kplus_0_PZ);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_ID", Bs_DTF_D_splus_Kplus_ID, &b_Bs_DTF_D_splus_Kplus_ID);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_PE", Bs_DTF_D_splus_Kplus_PE, &b_Bs_DTF_D_splus_Kplus_PE);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_PX", Bs_DTF_D_splus_Kplus_PX, &b_Bs_DTF_D_splus_Kplus_PX);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_PY", Bs_DTF_D_splus_Kplus_PY, &b_Bs_DTF_D_splus_Kplus_PY);
    fChain->SetBranchAddress("Bs_DTF_D_splus_Kplus_PZ", Bs_DTF_D_splus_Kplus_PZ, &b_Bs_DTF_D_splus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_DTF_D_splus_M", Bs_DTF_D_splus_M, &b_Bs_DTF_D_splus_M);
    fChain->SetBranchAddress("Bs_DTF_D_splus_MERR", Bs_DTF_D_splus_MERR, &b_Bs_DTF_D_splus_MERR);
    fChain->SetBranchAddress("Bs_DTF_D_splus_P", Bs_DTF_D_splus_P, &b_Bs_DTF_D_splus_P);
    fChain->SetBranchAddress("Bs_DTF_D_splus_PERR", Bs_DTF_D_splus_PERR, &b_Bs_DTF_D_splus_PERR);
    fChain->SetBranchAddress("Bs_DTF_D_splus_ctau", Bs_DTF_D_splus_ctau, &b_Bs_DTF_D_splus_ctau);
    fChain->SetBranchAddress("Bs_DTF_D_splus_ctauErr", Bs_DTF_D_splus_ctauErr, &b_Bs_DTF_D_splus_ctauErr);
    fChain->SetBranchAddress("Bs_DTF_D_splus_decayLength", Bs_DTF_D_splus_decayLength, &b_Bs_DTF_D_splus_decayLength);
    fChain->SetBranchAddress("Bs_DTF_D_splus_decayLengthErr", Bs_DTF_D_splus_decayLengthErr, &b_Bs_DTF_D_splus_decayLengthErr);
    fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_ID", Bs_DTF_D_splus_piplus_ID, &b_Bs_DTF_D_splus_piplus_ID);
    fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PE", Bs_DTF_D_splus_piplus_PE, &b_Bs_DTF_D_splus_piplus_PE);
    fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PX", Bs_DTF_D_splus_piplus_PX, &b_Bs_DTF_D_splus_piplus_PX);
    fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PY", Bs_DTF_D_splus_piplus_PY, &b_Bs_DTF_D_splus_piplus_PY);
    fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PZ", Bs_DTF_D_splus_piplus_PZ, &b_Bs_DTF_D_splus_piplus_PZ);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_Kplus_ID", Bs_DTF_K_1_1270_plus_Kplus_ID, &b_Bs_DTF_K_1_1270_plus_Kplus_ID);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_Kplus_PE", Bs_DTF_K_1_1270_plus_Kplus_PE, &b_Bs_DTF_K_1_1270_plus_Kplus_PE);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_Kplus_PX", Bs_DTF_K_1_1270_plus_Kplus_PX, &b_Bs_DTF_K_1_1270_plus_Kplus_PX);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_Kplus_PY", Bs_DTF_K_1_1270_plus_Kplus_PY, &b_Bs_DTF_K_1_1270_plus_Kplus_PY);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_Kplus_PZ", Bs_DTF_K_1_1270_plus_Kplus_PZ, &b_Bs_DTF_K_1_1270_plus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_M", Bs_DTF_K_1_1270_plus_M, &b_Bs_DTF_K_1_1270_plus_M);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_MERR", Bs_DTF_K_1_1270_plus_MERR, &b_Bs_DTF_K_1_1270_plus_MERR);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_P", Bs_DTF_K_1_1270_plus_P, &b_Bs_DTF_K_1_1270_plus_P);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_PERR", Bs_DTF_K_1_1270_plus_PERR, &b_Bs_DTF_K_1_1270_plus_PERR);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_ctau", Bs_DTF_K_1_1270_plus_ctau, &b_Bs_DTF_K_1_1270_plus_ctau);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_ctauErr", Bs_DTF_K_1_1270_plus_ctauErr, &b_Bs_DTF_K_1_1270_plus_ctauErr);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_decayLength", Bs_DTF_K_1_1270_plus_decayLength, &b_Bs_DTF_K_1_1270_plus_decayLength);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_decayLengthErr", Bs_DTF_K_1_1270_plus_decayLengthErr, &b_Bs_DTF_K_1_1270_plus_decayLengthErr);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_0_ID", Bs_DTF_K_1_1270_plus_piplus_0_ID, &b_Bs_DTF_K_1_1270_plus_piplus_0_ID);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_0_PE", Bs_DTF_K_1_1270_plus_piplus_0_PE, &b_Bs_DTF_K_1_1270_plus_piplus_0_PE);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_0_PX", Bs_DTF_K_1_1270_plus_piplus_0_PX, &b_Bs_DTF_K_1_1270_plus_piplus_0_PX);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_0_PY", Bs_DTF_K_1_1270_plus_piplus_0_PY, &b_Bs_DTF_K_1_1270_plus_piplus_0_PY);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_0_PZ", Bs_DTF_K_1_1270_plus_piplus_0_PZ, &b_Bs_DTF_K_1_1270_plus_piplus_0_PZ);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_ID", Bs_DTF_K_1_1270_plus_piplus_ID, &b_Bs_DTF_K_1_1270_plus_piplus_ID);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_PE", Bs_DTF_K_1_1270_plus_piplus_PE, &b_Bs_DTF_K_1_1270_plus_piplus_PE);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_PX", Bs_DTF_K_1_1270_plus_piplus_PX, &b_Bs_DTF_K_1_1270_plus_piplus_PX);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_PY", Bs_DTF_K_1_1270_plus_piplus_PY, &b_Bs_DTF_K_1_1270_plus_piplus_PY);
    fChain->SetBranchAddress("Bs_DTF_K_1_1270_plus_piplus_PZ", Bs_DTF_K_1_1270_plus_piplus_PZ, &b_Bs_DTF_K_1_1270_plus_piplus_PZ);
    fChain->SetBranchAddress("Bs_DTF_M", Bs_DTF_M, &b_Bs_DTF_M);
    fChain->SetBranchAddress("Bs_DTF_MERR", Bs_DTF_MERR, &b_Bs_DTF_MERR);
    fChain->SetBranchAddress("Bs_DTF_P", Bs_DTF_P, &b_Bs_DTF_P);
    fChain->SetBranchAddress("Bs_DTF_PERR", Bs_DTF_PERR, &b_Bs_DTF_PERR);
    fChain->SetBranchAddress("Bs_DTF_PV_X", Bs_DTF_PV_X, &b_Bs_DTF_PV_X);
    fChain->SetBranchAddress("Bs_DTF_PV_Y", Bs_DTF_PV_Y, &b_Bs_DTF_PV_Y);
    fChain->SetBranchAddress("Bs_DTF_PV_Z", Bs_DTF_PV_Z, &b_Bs_DTF_PV_Z);
    fChain->SetBranchAddress("Bs_DTF_PV_key", Bs_DTF_PV_key, &b_Bs_DTF_PV_key);
    fChain->SetBranchAddress("Bs_DTF_chi2", Bs_DTF_chi2, &b_Bs_DTF_chi2);
    fChain->SetBranchAddress("Bs_DTF_ctau", Bs_DTF_ctau, &b_Bs_DTF_ctau);
    fChain->SetBranchAddress("Bs_DTF_ctauErr", Bs_DTF_ctauErr, &b_Bs_DTF_ctauErr);
    fChain->SetBranchAddress("Bs_DTF_decayLength", Bs_DTF_decayLength, &b_Bs_DTF_decayLength);
    fChain->SetBranchAddress("Bs_DTF_decayLengthErr", Bs_DTF_decayLengthErr, &b_Bs_DTF_decayLengthErr);
    fChain->SetBranchAddress("Bs_DTF_nDOF", Bs_DTF_nDOF, &b_Bs_DTF_nDOF);
    fChain->SetBranchAddress("Bs_DTF_nIter", Bs_DTF_nIter, &b_Bs_DTF_nIter);
    fChain->SetBranchAddress("Bs_DTF_status", Bs_DTF_status, &b_Bs_DTF_status);
    fChain->SetBranchAddress("Bs_PV_nPV", &Bs_PV_nPV, &b_Bs_PV_nPV);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_0_ID", Bs_PV_Dplus_Kplus_0_ID, &b_Bs_PV_Dplus_Kplus_0_ID);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_0_PE", Bs_PV_Dplus_Kplus_0_PE, &b_Bs_PV_Dplus_Kplus_0_PE);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_0_PX", Bs_PV_Dplus_Kplus_0_PX, &b_Bs_PV_Dplus_Kplus_0_PX);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_0_PY", Bs_PV_Dplus_Kplus_0_PY, &b_Bs_PV_Dplus_Kplus_0_PY);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_0_PZ", Bs_PV_Dplus_Kplus_0_PZ, &b_Bs_PV_Dplus_Kplus_0_PZ);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_ID", Bs_PV_Dplus_Kplus_ID, &b_Bs_PV_Dplus_Kplus_ID);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_PE", Bs_PV_Dplus_Kplus_PE, &b_Bs_PV_Dplus_Kplus_PE);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_PX", Bs_PV_Dplus_Kplus_PX, &b_Bs_PV_Dplus_Kplus_PX);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_PY", Bs_PV_Dplus_Kplus_PY, &b_Bs_PV_Dplus_Kplus_PY);
    fChain->SetBranchAddress("Bs_PV_Dplus_Kplus_PZ", Bs_PV_Dplus_Kplus_PZ, &b_Bs_PV_Dplus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_PV_Dplus_M", Bs_PV_Dplus_M, &b_Bs_PV_Dplus_M);
    fChain->SetBranchAddress("Bs_PV_Dplus_MERR", Bs_PV_Dplus_MERR, &b_Bs_PV_Dplus_MERR);
    fChain->SetBranchAddress("Bs_PV_Dplus_P", Bs_PV_Dplus_P, &b_Bs_PV_Dplus_P);
    fChain->SetBranchAddress("Bs_PV_Dplus_PERR", Bs_PV_Dplus_PERR, &b_Bs_PV_Dplus_PERR);
    fChain->SetBranchAddress("Bs_PV_Dplus_ctau", Bs_PV_Dplus_ctau, &b_Bs_PV_Dplus_ctau);
    fChain->SetBranchAddress("Bs_PV_Dplus_ctauErr", Bs_PV_Dplus_ctauErr, &b_Bs_PV_Dplus_ctauErr);
    fChain->SetBranchAddress("Bs_PV_Dplus_decayLength", Bs_PV_Dplus_decayLength, &b_Bs_PV_Dplus_decayLength);
    fChain->SetBranchAddress("Bs_PV_Dplus_decayLengthErr", Bs_PV_Dplus_decayLengthErr, &b_Bs_PV_Dplus_decayLengthErr);
    fChain->SetBranchAddress("Bs_PV_Dplus_piplus_ID", Bs_PV_Dplus_piplus_ID, &b_Bs_PV_Dplus_piplus_ID);
    fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PE", Bs_PV_Dplus_piplus_PE, &b_Bs_PV_Dplus_piplus_PE);
    fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PX", Bs_PV_Dplus_piplus_PX, &b_Bs_PV_Dplus_piplus_PX);
    fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PY", Bs_PV_Dplus_piplus_PY, &b_Bs_PV_Dplus_piplus_PY);
    fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PZ", Bs_PV_Dplus_piplus_PZ, &b_Bs_PV_Dplus_piplus_PZ);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_Kplus_ID", Bs_PV_K_1_1270_plus_Kplus_ID, &b_Bs_PV_K_1_1270_plus_Kplus_ID);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_Kplus_PE", Bs_PV_K_1_1270_plus_Kplus_PE, &b_Bs_PV_K_1_1270_plus_Kplus_PE);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_Kplus_PX", Bs_PV_K_1_1270_plus_Kplus_PX, &b_Bs_PV_K_1_1270_plus_Kplus_PX);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_Kplus_PY", Bs_PV_K_1_1270_plus_Kplus_PY, &b_Bs_PV_K_1_1270_plus_Kplus_PY);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_Kplus_PZ", Bs_PV_K_1_1270_plus_Kplus_PZ, &b_Bs_PV_K_1_1270_plus_Kplus_PZ);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_M", Bs_PV_K_1_1270_plus_M, &b_Bs_PV_K_1_1270_plus_M);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_MERR", Bs_PV_K_1_1270_plus_MERR, &b_Bs_PV_K_1_1270_plus_MERR);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_P", Bs_PV_K_1_1270_plus_P, &b_Bs_PV_K_1_1270_plus_P);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_PERR", Bs_PV_K_1_1270_plus_PERR, &b_Bs_PV_K_1_1270_plus_PERR);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_ctau", Bs_PV_K_1_1270_plus_ctau, &b_Bs_PV_K_1_1270_plus_ctau);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_ctauErr", Bs_PV_K_1_1270_plus_ctauErr, &b_Bs_PV_K_1_1270_plus_ctauErr);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_decayLength", Bs_PV_K_1_1270_plus_decayLength, &b_Bs_PV_K_1_1270_plus_decayLength);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_decayLengthErr", Bs_PV_K_1_1270_plus_decayLengthErr, &b_Bs_PV_K_1_1270_plus_decayLengthErr);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_0_ID", Bs_PV_K_1_1270_plus_piplus_0_ID, &b_Bs_PV_K_1_1270_plus_piplus_0_ID);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_0_PE", Bs_PV_K_1_1270_plus_piplus_0_PE, &b_Bs_PV_K_1_1270_plus_piplus_0_PE);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_0_PX", Bs_PV_K_1_1270_plus_piplus_0_PX, &b_Bs_PV_K_1_1270_plus_piplus_0_PX);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_0_PY", Bs_PV_K_1_1270_plus_piplus_0_PY, &b_Bs_PV_K_1_1270_plus_piplus_0_PY);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_0_PZ", Bs_PV_K_1_1270_plus_piplus_0_PZ, &b_Bs_PV_K_1_1270_plus_piplus_0_PZ);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_ID", Bs_PV_K_1_1270_plus_piplus_ID, &b_Bs_PV_K_1_1270_plus_piplus_ID);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_PE", Bs_PV_K_1_1270_plus_piplus_PE, &b_Bs_PV_K_1_1270_plus_piplus_PE);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_PX", Bs_PV_K_1_1270_plus_piplus_PX, &b_Bs_PV_K_1_1270_plus_piplus_PX);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_PY", Bs_PV_K_1_1270_plus_piplus_PY, &b_Bs_PV_K_1_1270_plus_piplus_PY);
    fChain->SetBranchAddress("Bs_PV_K_1_1270_plus_piplus_PZ", Bs_PV_K_1_1270_plus_piplus_PZ, &b_Bs_PV_K_1_1270_plus_piplus_PZ);
    fChain->SetBranchAddress("Bs_PV_M", Bs_PV_M, &b_Bs_PV_M);
    fChain->SetBranchAddress("Bs_PV_MERR", Bs_PV_MERR, &b_Bs_PV_MERR);
    fChain->SetBranchAddress("Bs_PV_P", Bs_PV_P, &b_Bs_PV_P);
    fChain->SetBranchAddress("Bs_PV_PERR", Bs_PV_PERR, &b_Bs_PV_PERR);
    fChain->SetBranchAddress("Bs_PV_PV_X", Bs_PV_PV_X, &b_Bs_PV_PV_X);
    fChain->SetBranchAddress("Bs_PV_PV_Y", Bs_PV_PV_Y, &b_Bs_PV_PV_Y);
    fChain->SetBranchAddress("Bs_PV_PV_Z", Bs_PV_PV_Z, &b_Bs_PV_PV_Z);
    fChain->SetBranchAddress("Bs_PV_PV_key", Bs_PV_PV_key, &b_Bs_PV_PV_key);
    fChain->SetBranchAddress("Bs_PV_chi2", Bs_PV_chi2, &b_Bs_PV_chi2);
    fChain->SetBranchAddress("Bs_PV_ctau", Bs_PV_ctau, &b_Bs_PV_ctau);
    fChain->SetBranchAddress("Bs_PV_ctauErr", Bs_PV_ctauErr, &b_Bs_PV_ctauErr);
    fChain->SetBranchAddress("Bs_PV_decayLength", Bs_PV_decayLength, &b_Bs_PV_decayLength);
    fChain->SetBranchAddress("Bs_PV_decayLengthErr", Bs_PV_decayLengthErr, &b_Bs_PV_decayLengthErr);
    fChain->SetBranchAddress("Bs_PV_nDOF", Bs_PV_nDOF, &b_Bs_PV_nDOF);
    fChain->SetBranchAddress("Bs_PV_nIter", Bs_PV_nIter, &b_Bs_PV_nIter);
    fChain->SetBranchAddress("Bs_PV_status", Bs_PV_status, &b_Bs_PV_status);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION", &Bs_BsTaggingTool_TAGDECISION, &b_Bs_BsTaggingTool_TAGDECISION);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA", &Bs_BsTaggingTool_TAGOMEGA, &b_Bs_BsTaggingTool_TAGOMEGA);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION_OS", &Bs_BsTaggingTool_TAGDECISION_OS, &b_Bs_BsTaggingTool_TAGDECISION_OS);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA_OS", &Bs_BsTaggingTool_TAGOMEGA_OS, &b_Bs_BsTaggingTool_TAGOMEGA_OS);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGGER", &Bs_BsTaggingTool_TAGGER, &b_Bs_BsTaggingTool_TAGGER);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Muon_DEC", &Bs_BsTaggingTool_OS_Muon_DEC, &b_Bs_BsTaggingTool_OS_Muon_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Muon_PROB", &Bs_BsTaggingTool_OS_Muon_PROB, &b_Bs_BsTaggingTool_OS_Muon_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Electron_DEC", &Bs_BsTaggingTool_OS_Electron_DEC, &b_Bs_BsTaggingTool_OS_Electron_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Electron_PROB", &Bs_BsTaggingTool_OS_Electron_PROB, &b_Bs_BsTaggingTool_OS_Electron_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Kaon_DEC", &Bs_BsTaggingTool_OS_Kaon_DEC, &b_Bs_BsTaggingTool_OS_Kaon_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Kaon_PROB", &Bs_BsTaggingTool_OS_Kaon_PROB, &b_Bs_BsTaggingTool_OS_Kaon_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Kaon_DEC", &Bs_BsTaggingTool_SS_Kaon_DEC, &b_Bs_BsTaggingTool_SS_Kaon_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Kaon_PROB", &Bs_BsTaggingTool_SS_Kaon_PROB, &b_Bs_BsTaggingTool_SS_Kaon_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Pion_DEC", &Bs_BsTaggingTool_SS_Pion_DEC, &b_Bs_BsTaggingTool_SS_Pion_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Pion_PROB", &Bs_BsTaggingTool_SS_Pion_PROB, &b_Bs_BsTaggingTool_SS_Pion_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_PionBDT_DEC", &Bs_BsTaggingTool_SS_PionBDT_DEC, &b_Bs_BsTaggingTool_SS_PionBDT_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_PionBDT_PROB", &Bs_BsTaggingTool_SS_PionBDT_PROB, &b_Bs_BsTaggingTool_SS_PionBDT_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_VtxCharge_DEC", &Bs_BsTaggingTool_VtxCharge_DEC, &b_Bs_BsTaggingTool_VtxCharge_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_VtxCharge_PROB", &Bs_BsTaggingTool_VtxCharge_PROB, &b_Bs_BsTaggingTool_VtxCharge_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_nnetKaon_DEC", &Bs_BsTaggingTool_OS_nnetKaon_DEC, &b_Bs_BsTaggingTool_OS_nnetKaon_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_nnetKaon_PROB", &Bs_BsTaggingTool_OS_nnetKaon_PROB, &b_Bs_BsTaggingTool_OS_nnetKaon_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_nnetKaon_DEC", &Bs_BsTaggingTool_SS_nnetKaon_DEC, &b_Bs_BsTaggingTool_SS_nnetKaon_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_nnetKaon_PROB", &Bs_BsTaggingTool_SS_nnetKaon_PROB, &b_Bs_BsTaggingTool_SS_nnetKaon_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Proton_DEC", &Bs_BsTaggingTool_SS_Proton_DEC, &b_Bs_BsTaggingTool_SS_Proton_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_SS_Proton_PROB", &Bs_BsTaggingTool_SS_Proton_PROB, &b_Bs_BsTaggingTool_SS_Proton_PROB);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Charm_DEC", &Bs_BsTaggingTool_OS_Charm_DEC, &b_Bs_BsTaggingTool_OS_Charm_DEC);
    fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Charm_PROB", &Bs_BsTaggingTool_OS_Charm_PROB, &b_Bs_BsTaggingTool_OS_Charm_PROB);
    fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
    fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
    fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
    fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
    fChain->SetBranchAddress("Ds_CosTheta", &Ds_CosTheta, &b_Ds_CosTheta);
    fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
    fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
    fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
    fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("Ds_ENDVERTEX_COV_", Ds_ENDVERTEX_COV_, &b_Ds_ENDVERTEX_COV_);
    fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
    fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
    fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
    fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
    fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
    fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
    fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
    fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
    fChain->SetBranchAddress("Ds_OWNPV_COV_", Ds_OWNPV_COV_, &b_Ds_OWNPV_COV_);
    fChain->SetBranchAddress("Ds_IP_OWNPV", &Ds_IP_OWNPV, &b_Ds_IP_OWNPV);
    fChain->SetBranchAddress("Ds_IPCHI2_OWNPV", &Ds_IPCHI2_OWNPV, &b_Ds_IPCHI2_OWNPV);
    fChain->SetBranchAddress("Ds_FD_OWNPV", &Ds_FD_OWNPV, &b_Ds_FD_OWNPV);
    fChain->SetBranchAddress("Ds_FDCHI2_OWNPV", &Ds_FDCHI2_OWNPV, &b_Ds_FDCHI2_OWNPV);
    fChain->SetBranchAddress("Ds_DIRA_OWNPV", &Ds_DIRA_OWNPV, &b_Ds_DIRA_OWNPV);
    fChain->SetBranchAddress("Ds_ORIVX_X", &Ds_ORIVX_X, &b_Ds_ORIVX_X);
    fChain->SetBranchAddress("Ds_ORIVX_Y", &Ds_ORIVX_Y, &b_Ds_ORIVX_Y);
    fChain->SetBranchAddress("Ds_ORIVX_Z", &Ds_ORIVX_Z, &b_Ds_ORIVX_Z);
    fChain->SetBranchAddress("Ds_ORIVX_XERR", &Ds_ORIVX_XERR, &b_Ds_ORIVX_XERR);
    fChain->SetBranchAddress("Ds_ORIVX_YERR", &Ds_ORIVX_YERR, &b_Ds_ORIVX_YERR);
    fChain->SetBranchAddress("Ds_ORIVX_ZERR", &Ds_ORIVX_ZERR, &b_Ds_ORIVX_ZERR);
    fChain->SetBranchAddress("Ds_ORIVX_CHI2", &Ds_ORIVX_CHI2, &b_Ds_ORIVX_CHI2);
    fChain->SetBranchAddress("Ds_ORIVX_NDOF", &Ds_ORIVX_NDOF, &b_Ds_ORIVX_NDOF);
    fChain->SetBranchAddress("Ds_ORIVX_COV_", Ds_ORIVX_COV_, &b_Ds_ORIVX_COV_);
    fChain->SetBranchAddress("Ds_FD_ORIVX", &Ds_FD_ORIVX, &b_Ds_FD_ORIVX);
    fChain->SetBranchAddress("Ds_FDCHI2_ORIVX", &Ds_FDCHI2_ORIVX, &b_Ds_FDCHI2_ORIVX);
    fChain->SetBranchAddress("Ds_DIRA_ORIVX", &Ds_DIRA_ORIVX, &b_Ds_DIRA_ORIVX);
    fChain->SetBranchAddress("Ds_P", &Ds_P, &b_Ds_P);
    fChain->SetBranchAddress("Ds_PT", &Ds_PT, &b_Ds_PT);
    fChain->SetBranchAddress("Ds_PE", &Ds_PE, &b_Ds_PE);
    fChain->SetBranchAddress("Ds_PX", &Ds_PX, &b_Ds_PX);
    fChain->SetBranchAddress("Ds_PY", &Ds_PY, &b_Ds_PY);
    fChain->SetBranchAddress("Ds_PZ", &Ds_PZ, &b_Ds_PZ);
    fChain->SetBranchAddress("Ds_MM", &Ds_MM, &b_Ds_MM);
    fChain->SetBranchAddress("Ds_MMERR", &Ds_MMERR, &b_Ds_MMERR);
    fChain->SetBranchAddress("Ds_M", &Ds_M, &b_Ds_M);
    fChain->SetBranchAddress("Ds_ID", &Ds_ID, &b_Ds_ID);
    fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
    fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);
    fChain->SetBranchAddress("Ds_TAUCHI2", &Ds_TAUCHI2, &b_Ds_TAUCHI2);
    fChain->SetBranchAddress("Ds_L0Global_Dec", &Ds_L0Global_Dec, &b_Ds_L0Global_Dec);
    fChain->SetBranchAddress("Ds_L0Global_TIS", &Ds_L0Global_TIS, &b_Ds_L0Global_TIS);
    fChain->SetBranchAddress("Ds_L0Global_TOS", &Ds_L0Global_TOS, &b_Ds_L0Global_TOS);
    fChain->SetBranchAddress("Ds_Hlt1Global_Dec", &Ds_Hlt1Global_Dec, &b_Ds_Hlt1Global_Dec);
    fChain->SetBranchAddress("Ds_Hlt1Global_TIS", &Ds_Hlt1Global_TIS, &b_Ds_Hlt1Global_TIS);
    fChain->SetBranchAddress("Ds_Hlt1Global_TOS", &Ds_Hlt1Global_TOS, &b_Ds_Hlt1Global_TOS);
    fChain->SetBranchAddress("Ds_Hlt1Phys_Dec", &Ds_Hlt1Phys_Dec, &b_Ds_Hlt1Phys_Dec);
    fChain->SetBranchAddress("Ds_Hlt1Phys_TIS", &Ds_Hlt1Phys_TIS, &b_Ds_Hlt1Phys_TIS);
    fChain->SetBranchAddress("Ds_Hlt1Phys_TOS", &Ds_Hlt1Phys_TOS, &b_Ds_Hlt1Phys_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Global_Dec", &Ds_Hlt2Global_Dec, &b_Ds_Hlt2Global_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Global_TIS", &Ds_Hlt2Global_TIS, &b_Ds_Hlt2Global_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Global_TOS", &Ds_Hlt2Global_TOS, &b_Ds_Hlt2Global_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Phys_Dec", &Ds_Hlt2Phys_Dec, &b_Ds_Hlt2Phys_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Phys_TIS", &Ds_Hlt2Phys_TIS, &b_Ds_Hlt2Phys_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Phys_TOS", &Ds_Hlt2Phys_TOS, &b_Ds_Hlt2Phys_TOS);
    fChain->SetBranchAddress("Ds_L0HadronDecision_Dec", &Ds_L0HadronDecision_Dec, &b_Ds_L0HadronDecision_Dec);
    fChain->SetBranchAddress("Ds_L0HadronDecision_TIS", &Ds_L0HadronDecision_TIS, &b_Ds_L0HadronDecision_TIS);
    fChain->SetBranchAddress("Ds_L0HadronDecision_TOS", &Ds_L0HadronDecision_TOS, &b_Ds_L0HadronDecision_TOS);
    fChain->SetBranchAddress("Ds_L0MuonDecision_Dec", &Ds_L0MuonDecision_Dec, &b_Ds_L0MuonDecision_Dec);
    fChain->SetBranchAddress("Ds_L0MuonDecision_TIS", &Ds_L0MuonDecision_TIS, &b_Ds_L0MuonDecision_TIS);
    fChain->SetBranchAddress("Ds_L0MuonDecision_TOS", &Ds_L0MuonDecision_TOS, &b_Ds_L0MuonDecision_TOS);
    fChain->SetBranchAddress("Ds_L0GlobalDecision_Dec", &Ds_L0GlobalDecision_Dec, &b_Ds_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("Ds_L0GlobalDecision_TIS", &Ds_L0GlobalDecision_TIS, &b_Ds_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("Ds_L0GlobalDecision_TOS", &Ds_L0GlobalDecision_TOS, &b_Ds_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt1TrackAllL0Decision_Dec", &Ds_Hlt1TrackAllL0Decision_Dec, &b_Ds_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("Ds_Hlt1TrackAllL0Decision_TIS", &Ds_Hlt1TrackAllL0Decision_TIS, &b_Ds_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("Ds_Hlt1TrackAllL0Decision_TOS", &Ds_Hlt1TrackAllL0Decision_TOS, &b_Ds_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVADecision_Dec", &Ds_Hlt1TrackMVADecision_Dec, &b_Ds_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVADecision_TIS", &Ds_Hlt1TrackMVADecision_TIS, &b_Ds_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVADecision_TOS", &Ds_Hlt1TrackMVADecision_TOS, &b_Ds_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVADecision_Dec", &Ds_Hlt1TwoTrackMVADecision_Dec, &b_Ds_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVADecision_TIS", &Ds_Hlt1TwoTrackMVADecision_TIS, &b_Ds_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVADecision_TOS", &Ds_Hlt1TwoTrackMVADecision_TOS, &b_Ds_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVALooseDecision_Dec", &Ds_Hlt1TrackMVALooseDecision_Dec, &b_Ds_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVALooseDecision_TIS", &Ds_Hlt1TrackMVALooseDecision_TIS, &b_Ds_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt1TrackMVALooseDecision_TOS", &Ds_Hlt1TrackMVALooseDecision_TOS, &b_Ds_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVALooseDecision_Dec", &Ds_Hlt1TwoTrackMVALooseDecision_Dec, &b_Ds_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVALooseDecision_TIS", &Ds_Hlt1TwoTrackMVALooseDecision_TIS, &b_Ds_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt1TwoTrackMVALooseDecision_TOS", &Ds_Hlt1TwoTrackMVALooseDecision_TOS, &b_Ds_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2IncPhiDecision_Dec", &Ds_Hlt2IncPhiDecision_Dec, &b_Ds_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2IncPhiDecision_TIS", &Ds_Hlt2IncPhiDecision_TIS, &b_Ds_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2IncPhiDecision_TOS", &Ds_Hlt2IncPhiDecision_TOS, &b_Ds_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2PhiIncPhiDecision_Dec", &Ds_Hlt2PhiIncPhiDecision_Dec, &b_Ds_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2PhiIncPhiDecision_TIS", &Ds_Hlt2PhiIncPhiDecision_TIS, &b_Ds_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2PhiIncPhiDecision_TOS", &Ds_Hlt2PhiIncPhiDecision_TOS, &b_Ds_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyBBDTDecision_Dec", &Ds_Hlt2Topo2BodyBBDTDecision_Dec, &b_Ds_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyBBDTDecision_TIS", &Ds_Hlt2Topo2BodyBBDTDecision_TIS, &b_Ds_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyBBDTDecision_TOS", &Ds_Hlt2Topo2BodyBBDTDecision_TOS, &b_Ds_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyBBDTDecision_Dec", &Ds_Hlt2Topo3BodyBBDTDecision_Dec, &b_Ds_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyBBDTDecision_TIS", &Ds_Hlt2Topo3BodyBBDTDecision_TIS, &b_Ds_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyBBDTDecision_TOS", &Ds_Hlt2Topo3BodyBBDTDecision_TOS, &b_Ds_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyBBDTDecision_Dec", &Ds_Hlt2Topo4BodyBBDTDecision_Dec, &b_Ds_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyBBDTDecision_TIS", &Ds_Hlt2Topo4BodyBBDTDecision_TIS, &b_Ds_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyBBDTDecision_TOS", &Ds_Hlt2Topo4BodyBBDTDecision_TOS, &b_Ds_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyDecision_Dec", &Ds_Hlt2Topo2BodyDecision_Dec, &b_Ds_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyDecision_TIS", &Ds_Hlt2Topo2BodyDecision_TIS, &b_Ds_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo2BodyDecision_TOS", &Ds_Hlt2Topo2BodyDecision_TOS, &b_Ds_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyDecision_Dec", &Ds_Hlt2Topo3BodyDecision_Dec, &b_Ds_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyDecision_TIS", &Ds_Hlt2Topo3BodyDecision_TIS, &b_Ds_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo3BodyDecision_TOS", &Ds_Hlt2Topo3BodyDecision_TOS, &b_Ds_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyDecision_Dec", &Ds_Hlt2Topo4BodyDecision_Dec, &b_Ds_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyDecision_TIS", &Ds_Hlt2Topo4BodyDecision_TIS, &b_Ds_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("Ds_Hlt2Topo4BodyDecision_TOS", &Ds_Hlt2Topo4BodyDecision_TOS, &b_Ds_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("Ds_cpx_0.50", &Ds_cpx_0_50, &b_Ds_cpx_0_50);
    fChain->SetBranchAddress("Ds_cpy_0.50", &Ds_cpy_0_50, &b_Ds_cpy_0_50);
    fChain->SetBranchAddress("Ds_cpz_0.50", &Ds_cpz_0_50, &b_Ds_cpz_0_50);
    fChain->SetBranchAddress("Ds_cpt_0.50", &Ds_cpt_0_50, &b_Ds_cpt_0_50);
    fChain->SetBranchAddress("Ds_cp_0.50", &Ds_cp_0_50, &b_Ds_cp_0_50);
    fChain->SetBranchAddress("Ds_cmult_0.50", &Ds_cmult_0_50, &b_Ds_cmult_0_50);
    fChain->SetBranchAddress("Ds_deltaEta_0.50", &Ds_deltaEta_0_50, &b_Ds_deltaEta_0_50);
    fChain->SetBranchAddress("Ds_deltaPhi_0.50", &Ds_deltaPhi_0_50, &b_Ds_deltaPhi_0_50);
    fChain->SetBranchAddress("Ds_pxasy_0.50", &Ds_pxasy_0_50, &b_Ds_pxasy_0_50);
    fChain->SetBranchAddress("Ds_pyasy_0.50", &Ds_pyasy_0_50, &b_Ds_pyasy_0_50);
    fChain->SetBranchAddress("Ds_pzasy_0.50", &Ds_pzasy_0_50, &b_Ds_pzasy_0_50);
    fChain->SetBranchAddress("Ds_pasy_0.50", &Ds_pasy_0_50, &b_Ds_pasy_0_50);
    fChain->SetBranchAddress("Ds_ptasy_0.50", &Ds_ptasy_0_50, &b_Ds_ptasy_0_50);
    fChain->SetBranchAddress("Ds_cpx_0.60", &Ds_cpx_0_60, &b_Ds_cpx_0_60);
    fChain->SetBranchAddress("Ds_cpy_0.60", &Ds_cpy_0_60, &b_Ds_cpy_0_60);
    fChain->SetBranchAddress("Ds_cpz_0.60", &Ds_cpz_0_60, &b_Ds_cpz_0_60);
    fChain->SetBranchAddress("Ds_cpt_0.60", &Ds_cpt_0_60, &b_Ds_cpt_0_60);
    fChain->SetBranchAddress("Ds_cp_0.60", &Ds_cp_0_60, &b_Ds_cp_0_60);
    fChain->SetBranchAddress("Ds_cmult_0.60", &Ds_cmult_0_60, &b_Ds_cmult_0_60);
    fChain->SetBranchAddress("Ds_deltaEta_0.60", &Ds_deltaEta_0_60, &b_Ds_deltaEta_0_60);
    fChain->SetBranchAddress("Ds_deltaPhi_0.60", &Ds_deltaPhi_0_60, &b_Ds_deltaPhi_0_60);
    fChain->SetBranchAddress("Ds_pxasy_0.60", &Ds_pxasy_0_60, &b_Ds_pxasy_0_60);
    fChain->SetBranchAddress("Ds_pyasy_0.60", &Ds_pyasy_0_60, &b_Ds_pyasy_0_60);
    fChain->SetBranchAddress("Ds_pzasy_0.60", &Ds_pzasy_0_60, &b_Ds_pzasy_0_60);
    fChain->SetBranchAddress("Ds_pasy_0.60", &Ds_pasy_0_60, &b_Ds_pasy_0_60);
    fChain->SetBranchAddress("Ds_ptasy_0.60", &Ds_ptasy_0_60, &b_Ds_ptasy_0_60);
    fChain->SetBranchAddress("Ds_cpx_0.70", &Ds_cpx_0_70, &b_Ds_cpx_0_70);
    fChain->SetBranchAddress("Ds_cpy_0.70", &Ds_cpy_0_70, &b_Ds_cpy_0_70);
    fChain->SetBranchAddress("Ds_cpz_0.70", &Ds_cpz_0_70, &b_Ds_cpz_0_70);
    fChain->SetBranchAddress("Ds_cpt_0.70", &Ds_cpt_0_70, &b_Ds_cpt_0_70);
    fChain->SetBranchAddress("Ds_cp_0.70", &Ds_cp_0_70, &b_Ds_cp_0_70);
    fChain->SetBranchAddress("Ds_cmult_0.70", &Ds_cmult_0_70, &b_Ds_cmult_0_70);
    fChain->SetBranchAddress("Ds_deltaEta_0.70", &Ds_deltaEta_0_70, &b_Ds_deltaEta_0_70);
    fChain->SetBranchAddress("Ds_deltaPhi_0.70", &Ds_deltaPhi_0_70, &b_Ds_deltaPhi_0_70);
    fChain->SetBranchAddress("Ds_pxasy_0.70", &Ds_pxasy_0_70, &b_Ds_pxasy_0_70);
    fChain->SetBranchAddress("Ds_pyasy_0.70", &Ds_pyasy_0_70, &b_Ds_pyasy_0_70);
    fChain->SetBranchAddress("Ds_pzasy_0.70", &Ds_pzasy_0_70, &b_Ds_pzasy_0_70);
    fChain->SetBranchAddress("Ds_pasy_0.70", &Ds_pasy_0_70, &b_Ds_pasy_0_70);
    fChain->SetBranchAddress("Ds_ptasy_0.70", &Ds_ptasy_0_70, &b_Ds_ptasy_0_70);
    fChain->SetBranchAddress("Ds_cpx_0.80", &Ds_cpx_0_80, &b_Ds_cpx_0_80);
    fChain->SetBranchAddress("Ds_cpy_0.80", &Ds_cpy_0_80, &b_Ds_cpy_0_80);
    fChain->SetBranchAddress("Ds_cpz_0.80", &Ds_cpz_0_80, &b_Ds_cpz_0_80);
    fChain->SetBranchAddress("Ds_cpt_0.80", &Ds_cpt_0_80, &b_Ds_cpt_0_80);
    fChain->SetBranchAddress("Ds_cp_0.80", &Ds_cp_0_80, &b_Ds_cp_0_80);
    fChain->SetBranchAddress("Ds_cmult_0.80", &Ds_cmult_0_80, &b_Ds_cmult_0_80);
    fChain->SetBranchAddress("Ds_deltaEta_0.80", &Ds_deltaEta_0_80, &b_Ds_deltaEta_0_80);
    fChain->SetBranchAddress("Ds_deltaPhi_0.80", &Ds_deltaPhi_0_80, &b_Ds_deltaPhi_0_80);
    fChain->SetBranchAddress("Ds_pxasy_0.80", &Ds_pxasy_0_80, &b_Ds_pxasy_0_80);
    fChain->SetBranchAddress("Ds_pyasy_0.80", &Ds_pyasy_0_80, &b_Ds_pyasy_0_80);
    fChain->SetBranchAddress("Ds_pzasy_0.80", &Ds_pzasy_0_80, &b_Ds_pzasy_0_80);
    fChain->SetBranchAddress("Ds_pasy_0.80", &Ds_pasy_0_80, &b_Ds_pasy_0_80);
    fChain->SetBranchAddress("Ds_ptasy_0.80", &Ds_ptasy_0_80, &b_Ds_ptasy_0_80);
    fChain->SetBranchAddress("Ds_cpx_0.90", &Ds_cpx_0_90, &b_Ds_cpx_0_90);
    fChain->SetBranchAddress("Ds_cpy_0.90", &Ds_cpy_0_90, &b_Ds_cpy_0_90);
    fChain->SetBranchAddress("Ds_cpz_0.90", &Ds_cpz_0_90, &b_Ds_cpz_0_90);
    fChain->SetBranchAddress("Ds_cpt_0.90", &Ds_cpt_0_90, &b_Ds_cpt_0_90);
    fChain->SetBranchAddress("Ds_cp_0.90", &Ds_cp_0_90, &b_Ds_cp_0_90);
    fChain->SetBranchAddress("Ds_cmult_0.90", &Ds_cmult_0_90, &b_Ds_cmult_0_90);
    fChain->SetBranchAddress("Ds_deltaEta_0.90", &Ds_deltaEta_0_90, &b_Ds_deltaEta_0_90);
    fChain->SetBranchAddress("Ds_deltaPhi_0.90", &Ds_deltaPhi_0_90, &b_Ds_deltaPhi_0_90);
    fChain->SetBranchAddress("Ds_pxasy_0.90", &Ds_pxasy_0_90, &b_Ds_pxasy_0_90);
    fChain->SetBranchAddress("Ds_pyasy_0.90", &Ds_pyasy_0_90, &b_Ds_pyasy_0_90);
    fChain->SetBranchAddress("Ds_pzasy_0.90", &Ds_pzasy_0_90, &b_Ds_pzasy_0_90);
    fChain->SetBranchAddress("Ds_pasy_0.90", &Ds_pasy_0_90, &b_Ds_pasy_0_90);
    fChain->SetBranchAddress("Ds_ptasy_0.90", &Ds_ptasy_0_90, &b_Ds_ptasy_0_90);
    fChain->SetBranchAddress("Ds_cpx_1.00", &Ds_cpx_1_00, &b_Ds_cpx_1_00);
    fChain->SetBranchAddress("Ds_cpy_1.00", &Ds_cpy_1_00, &b_Ds_cpy_1_00);
    fChain->SetBranchAddress("Ds_cpz_1.00", &Ds_cpz_1_00, &b_Ds_cpz_1_00);
    fChain->SetBranchAddress("Ds_cpt_1.00", &Ds_cpt_1_00, &b_Ds_cpt_1_00);
    fChain->SetBranchAddress("Ds_cp_1.00", &Ds_cp_1_00, &b_Ds_cp_1_00);
    fChain->SetBranchAddress("Ds_cmult_1.00", &Ds_cmult_1_00, &b_Ds_cmult_1_00);
    fChain->SetBranchAddress("Ds_deltaEta_1.00", &Ds_deltaEta_1_00, &b_Ds_deltaEta_1_00);
    fChain->SetBranchAddress("Ds_deltaPhi_1.00", &Ds_deltaPhi_1_00, &b_Ds_deltaPhi_1_00);
    fChain->SetBranchAddress("Ds_pxasy_1.00", &Ds_pxasy_1_00, &b_Ds_pxasy_1_00);
    fChain->SetBranchAddress("Ds_pyasy_1.00", &Ds_pyasy_1_00, &b_Ds_pyasy_1_00);
    fChain->SetBranchAddress("Ds_pzasy_1.00", &Ds_pzasy_1_00, &b_Ds_pzasy_1_00);
    fChain->SetBranchAddress("Ds_pasy_1.00", &Ds_pasy_1_00, &b_Ds_pasy_1_00);
    fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_DOCA1", &K_plus_fromDs_DOCA1, &b_K_plus_fromDs_DOCA1);
    fChain->SetBranchAddress("K_plus_fromDs_DOCA2", &K_plus_fromDs_DOCA2, &b_K_plus_fromDs_DOCA2);
    fChain->SetBranchAddress("K_plus_fromDs_DOCA3", &K_plus_fromDs_DOCA3, &b_K_plus_fromDs_DOCA3);
    fChain->SetBranchAddress("K_plus_fromDs_ETA", &K_plus_fromDs_ETA, &b_K_plus_fromDs_ETA);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNe", &K_plus_fromDs_MC12TuneV2_ProbNNe, &b_K_plus_fromDs_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNmu", &K_plus_fromDs_MC12TuneV2_ProbNNmu, &b_K_plus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNpi", &K_plus_fromDs_MC12TuneV2_ProbNNpi, &b_K_plus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNk", &K_plus_fromDs_MC12TuneV2_ProbNNk, &b_K_plus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNp", &K_plus_fromDs_MC12TuneV2_ProbNNp, &b_K_plus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNghost", &K_plus_fromDs_MC12TuneV2_ProbNNghost, &b_K_plus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNe", &K_plus_fromDs_MC12TuneV3_ProbNNe, &b_K_plus_fromDs_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNmu", &K_plus_fromDs_MC12TuneV3_ProbNNmu, &b_K_plus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNpi", &K_plus_fromDs_MC12TuneV3_ProbNNpi, &b_K_plus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNk", &K_plus_fromDs_MC12TuneV3_ProbNNk, &b_K_plus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNp", &K_plus_fromDs_MC12TuneV3_ProbNNp, &b_K_plus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNghost", &K_plus_fromDs_MC12TuneV3_ProbNNghost, &b_K_plus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNe", &K_plus_fromDs_MC12TuneV4_ProbNNe, &b_K_plus_fromDs_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNmu", &K_plus_fromDs_MC12TuneV4_ProbNNmu, &b_K_plus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNpi", &K_plus_fromDs_MC12TuneV4_ProbNNpi, &b_K_plus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNk", &K_plus_fromDs_MC12TuneV4_ProbNNk, &b_K_plus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNp", &K_plus_fromDs_MC12TuneV4_ProbNNp, &b_K_plus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNghost", &K_plus_fromDs_MC12TuneV4_ProbNNghost, &b_K_plus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNe", &K_plus_fromDs_MC15TuneV1_ProbNNe, &b_K_plus_fromDs_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNmu", &K_plus_fromDs_MC15TuneV1_ProbNNmu, &b_K_plus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNpi", &K_plus_fromDs_MC15TuneV1_ProbNNpi, &b_K_plus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNk", &K_plus_fromDs_MC15TuneV1_ProbNNk, &b_K_plus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNp", &K_plus_fromDs_MC15TuneV1_ProbNNp, &b_K_plus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNghost", &K_plus_fromDs_MC15TuneV1_ProbNNghost, &b_K_plus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_CosTheta", &K_plus_fromDs_CosTheta, &b_K_plus_fromDs_CosTheta);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_X", &K_plus_fromDs_OWNPV_X, &b_K_plus_fromDs_OWNPV_X);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_Y", &K_plus_fromDs_OWNPV_Y, &b_K_plus_fromDs_OWNPV_Y);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_Z", &K_plus_fromDs_OWNPV_Z, &b_K_plus_fromDs_OWNPV_Z);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_XERR", &K_plus_fromDs_OWNPV_XERR, &b_K_plus_fromDs_OWNPV_XERR);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_YERR", &K_plus_fromDs_OWNPV_YERR, &b_K_plus_fromDs_OWNPV_YERR);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_ZERR", &K_plus_fromDs_OWNPV_ZERR, &b_K_plus_fromDs_OWNPV_ZERR);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_CHI2", &K_plus_fromDs_OWNPV_CHI2, &b_K_plus_fromDs_OWNPV_CHI2);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_NDOF", &K_plus_fromDs_OWNPV_NDOF, &b_K_plus_fromDs_OWNPV_NDOF);
    fChain->SetBranchAddress("K_plus_fromDs_OWNPV_COV_", K_plus_fromDs_OWNPV_COV_, &b_K_plus_fromDs_OWNPV_COV_);
    fChain->SetBranchAddress("K_plus_fromDs_IP_OWNPV", &K_plus_fromDs_IP_OWNPV, &b_K_plus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("K_plus_fromDs_IPCHI2_OWNPV", &K_plus_fromDs_IPCHI2_OWNPV, &b_K_plus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_X", &K_plus_fromDs_ORIVX_X, &b_K_plus_fromDs_ORIVX_X);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_Y", &K_plus_fromDs_ORIVX_Y, &b_K_plus_fromDs_ORIVX_Y);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_Z", &K_plus_fromDs_ORIVX_Z, &b_K_plus_fromDs_ORIVX_Z);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_XERR", &K_plus_fromDs_ORIVX_XERR, &b_K_plus_fromDs_ORIVX_XERR);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_YERR", &K_plus_fromDs_ORIVX_YERR, &b_K_plus_fromDs_ORIVX_YERR);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_ZERR", &K_plus_fromDs_ORIVX_ZERR, &b_K_plus_fromDs_ORIVX_ZERR);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_CHI2", &K_plus_fromDs_ORIVX_CHI2, &b_K_plus_fromDs_ORIVX_CHI2);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_NDOF", &K_plus_fromDs_ORIVX_NDOF, &b_K_plus_fromDs_ORIVX_NDOF);
    fChain->SetBranchAddress("K_plus_fromDs_ORIVX_COV_", K_plus_fromDs_ORIVX_COV_, &b_K_plus_fromDs_ORIVX_COV_);
    fChain->SetBranchAddress("K_plus_fromDs_P", &K_plus_fromDs_P, &b_K_plus_fromDs_P);
    fChain->SetBranchAddress("K_plus_fromDs_PT", &K_plus_fromDs_PT, &b_K_plus_fromDs_PT);
    fChain->SetBranchAddress("K_plus_fromDs_PE", &K_plus_fromDs_PE, &b_K_plus_fromDs_PE);
    fChain->SetBranchAddress("K_plus_fromDs_PX", &K_plus_fromDs_PX, &b_K_plus_fromDs_PX);
    fChain->SetBranchAddress("K_plus_fromDs_PY", &K_plus_fromDs_PY, &b_K_plus_fromDs_PY);
    fChain->SetBranchAddress("K_plus_fromDs_PZ", &K_plus_fromDs_PZ, &b_K_plus_fromDs_PZ);
    fChain->SetBranchAddress("K_plus_fromDs_M", &K_plus_fromDs_M, &b_K_plus_fromDs_M);
    fChain->SetBranchAddress("K_plus_fromDs_ID", &K_plus_fromDs_ID, &b_K_plus_fromDs_ID);
    fChain->SetBranchAddress("K_plus_fromDs_PIDe", &K_plus_fromDs_PIDe, &b_K_plus_fromDs_PIDe);
    fChain->SetBranchAddress("K_plus_fromDs_PIDmu", &K_plus_fromDs_PIDmu, &b_K_plus_fromDs_PIDmu);
    fChain->SetBranchAddress("K_plus_fromDs_PIDK", &K_plus_fromDs_PIDK, &b_K_plus_fromDs_PIDK);
    fChain->SetBranchAddress("K_plus_fromDs_PIDp", &K_plus_fromDs_PIDp, &b_K_plus_fromDs_PIDp);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNe", &K_plus_fromDs_ProbNNe, &b_K_plus_fromDs_ProbNNe);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNk", &K_plus_fromDs_ProbNNk, &b_K_plus_fromDs_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNp", &K_plus_fromDs_ProbNNp, &b_K_plus_fromDs_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNpi", &K_plus_fromDs_ProbNNpi, &b_K_plus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNmu", &K_plus_fromDs_ProbNNmu, &b_K_plus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNghost", &K_plus_fromDs_ProbNNghost, &b_K_plus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_hasMuon", &K_plus_fromDs_hasMuon, &b_K_plus_fromDs_hasMuon);
    fChain->SetBranchAddress("K_plus_fromDs_isMuon", &K_plus_fromDs_isMuon, &b_K_plus_fromDs_isMuon);
    fChain->SetBranchAddress("K_plus_fromDs_hasRich", &K_plus_fromDs_hasRich, &b_K_plus_fromDs_hasRich);
    fChain->SetBranchAddress("K_plus_fromDs_UsedRichAerogel", &K_plus_fromDs_UsedRichAerogel, &b_K_plus_fromDs_UsedRichAerogel);
    fChain->SetBranchAddress("K_plus_fromDs_UsedRich1Gas", &K_plus_fromDs_UsedRich1Gas, &b_K_plus_fromDs_UsedRich1Gas);
    fChain->SetBranchAddress("K_plus_fromDs_UsedRich2Gas", &K_plus_fromDs_UsedRich2Gas, &b_K_plus_fromDs_UsedRich2Gas);
    fChain->SetBranchAddress("K_plus_fromDs_RichAboveElThres", &K_plus_fromDs_RichAboveElThres, &b_K_plus_fromDs_RichAboveElThres);
    fChain->SetBranchAddress("K_plus_fromDs_RichAboveMuThres", &K_plus_fromDs_RichAboveMuThres, &b_K_plus_fromDs_RichAboveMuThres);
    fChain->SetBranchAddress("K_plus_fromDs_RichAbovePiThres", &K_plus_fromDs_RichAbovePiThres, &b_K_plus_fromDs_RichAbovePiThres);
    fChain->SetBranchAddress("K_plus_fromDs_RichAboveKaThres", &K_plus_fromDs_RichAboveKaThres, &b_K_plus_fromDs_RichAboveKaThres);
    fChain->SetBranchAddress("K_plus_fromDs_RichAbovePrThres", &K_plus_fromDs_RichAbovePrThres, &b_K_plus_fromDs_RichAbovePrThres);
    fChain->SetBranchAddress("K_plus_fromDs_hasCalo", &K_plus_fromDs_hasCalo, &b_K_plus_fromDs_hasCalo);
    fChain->SetBranchAddress("K_plus_fromDs_L0Global_Dec", &K_plus_fromDs_L0Global_Dec, &b_K_plus_fromDs_L0Global_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_L0Global_TIS", &K_plus_fromDs_L0Global_TIS, &b_K_plus_fromDs_L0Global_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_L0Global_TOS", &K_plus_fromDs_L0Global_TOS, &b_K_plus_fromDs_L0Global_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Global_Dec", &K_plus_fromDs_Hlt1Global_Dec, &b_K_plus_fromDs_Hlt1Global_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Global_TIS", &K_plus_fromDs_Hlt1Global_TIS, &b_K_plus_fromDs_Hlt1Global_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Global_TOS", &K_plus_fromDs_Hlt1Global_TOS, &b_K_plus_fromDs_Hlt1Global_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Phys_Dec", &K_plus_fromDs_Hlt1Phys_Dec, &b_K_plus_fromDs_Hlt1Phys_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Phys_TIS", &K_plus_fromDs_Hlt1Phys_TIS, &b_K_plus_fromDs_Hlt1Phys_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1Phys_TOS", &K_plus_fromDs_Hlt1Phys_TOS, &b_K_plus_fromDs_Hlt1Phys_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Global_Dec", &K_plus_fromDs_Hlt2Global_Dec, &b_K_plus_fromDs_Hlt2Global_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Global_TIS", &K_plus_fromDs_Hlt2Global_TIS, &b_K_plus_fromDs_Hlt2Global_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Global_TOS", &K_plus_fromDs_Hlt2Global_TOS, &b_K_plus_fromDs_Hlt2Global_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Phys_Dec", &K_plus_fromDs_Hlt2Phys_Dec, &b_K_plus_fromDs_Hlt2Phys_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Phys_TIS", &K_plus_fromDs_Hlt2Phys_TIS, &b_K_plus_fromDs_Hlt2Phys_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Phys_TOS", &K_plus_fromDs_Hlt2Phys_TOS, &b_K_plus_fromDs_Hlt2Phys_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_L0HadronDecision_Dec", &K_plus_fromDs_L0HadronDecision_Dec, &b_K_plus_fromDs_L0HadronDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_L0HadronDecision_TIS", &K_plus_fromDs_L0HadronDecision_TIS, &b_K_plus_fromDs_L0HadronDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_L0HadronDecision_TOS", &K_plus_fromDs_L0HadronDecision_TOS, &b_K_plus_fromDs_L0HadronDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_L0MuonDecision_Dec", &K_plus_fromDs_L0MuonDecision_Dec, &b_K_plus_fromDs_L0MuonDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_L0MuonDecision_TIS", &K_plus_fromDs_L0MuonDecision_TIS, &b_K_plus_fromDs_L0MuonDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_L0MuonDecision_TOS", &K_plus_fromDs_L0MuonDecision_TOS, &b_K_plus_fromDs_L0MuonDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_L0GlobalDecision_Dec", &K_plus_fromDs_L0GlobalDecision_Dec, &b_K_plus_fromDs_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_L0GlobalDecision_TIS", &K_plus_fromDs_L0GlobalDecision_TIS, &b_K_plus_fromDs_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_L0GlobalDecision_TOS", &K_plus_fromDs_L0GlobalDecision_TOS, &b_K_plus_fromDs_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackAllL0Decision_Dec", &K_plus_fromDs_Hlt1TrackAllL0Decision_Dec, &b_K_plus_fromDs_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackAllL0Decision_TIS", &K_plus_fromDs_Hlt1TrackAllL0Decision_TIS, &b_K_plus_fromDs_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackAllL0Decision_TOS", &K_plus_fromDs_Hlt1TrackAllL0Decision_TOS, &b_K_plus_fromDs_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVADecision_Dec", &K_plus_fromDs_Hlt1TrackMVADecision_Dec, &b_K_plus_fromDs_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVADecision_TIS", &K_plus_fromDs_Hlt1TrackMVADecision_TIS, &b_K_plus_fromDs_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVADecision_TOS", &K_plus_fromDs_Hlt1TrackMVADecision_TOS, &b_K_plus_fromDs_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVADecision_Dec", &K_plus_fromDs_Hlt1TwoTrackMVADecision_Dec, &b_K_plus_fromDs_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVADecision_TIS", &K_plus_fromDs_Hlt1TwoTrackMVADecision_TIS, &b_K_plus_fromDs_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVADecision_TOS", &K_plus_fromDs_Hlt1TwoTrackMVADecision_TOS, &b_K_plus_fromDs_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVALooseDecision_Dec", &K_plus_fromDs_Hlt1TrackMVALooseDecision_Dec, &b_K_plus_fromDs_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVALooseDecision_TIS", &K_plus_fromDs_Hlt1TrackMVALooseDecision_TIS, &b_K_plus_fromDs_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TrackMVALooseDecision_TOS", &K_plus_fromDs_Hlt1TrackMVALooseDecision_TOS, &b_K_plus_fromDs_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec", &K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec, &b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS", &K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS, &b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS", &K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS, &b_K_plus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2IncPhiDecision_Dec", &K_plus_fromDs_Hlt2IncPhiDecision_Dec, &b_K_plus_fromDs_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2IncPhiDecision_TIS", &K_plus_fromDs_Hlt2IncPhiDecision_TIS, &b_K_plus_fromDs_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2IncPhiDecision_TOS", &K_plus_fromDs_Hlt2IncPhiDecision_TOS, &b_K_plus_fromDs_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2PhiIncPhiDecision_Dec", &K_plus_fromDs_Hlt2PhiIncPhiDecision_Dec, &b_K_plus_fromDs_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2PhiIncPhiDecision_TIS", &K_plus_fromDs_Hlt2PhiIncPhiDecision_TIS, &b_K_plus_fromDs_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2PhiIncPhiDecision_TOS", &K_plus_fromDs_Hlt2PhiIncPhiDecision_TOS, &b_K_plus_fromDs_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec", &K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec, &b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS", &K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS, &b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS", &K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS, &b_K_plus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec", &K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec, &b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS", &K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS, &b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS", &K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS, &b_K_plus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec", &K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec, &b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS", &K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS, &b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS", &K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS, &b_K_plus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyDecision_Dec", &K_plus_fromDs_Hlt2Topo2BodyDecision_Dec, &b_K_plus_fromDs_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyDecision_TIS", &K_plus_fromDs_Hlt2Topo2BodyDecision_TIS, &b_K_plus_fromDs_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo2BodyDecision_TOS", &K_plus_fromDs_Hlt2Topo2BodyDecision_TOS, &b_K_plus_fromDs_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyDecision_Dec", &K_plus_fromDs_Hlt2Topo3BodyDecision_Dec, &b_K_plus_fromDs_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyDecision_TIS", &K_plus_fromDs_Hlt2Topo3BodyDecision_TIS, &b_K_plus_fromDs_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo3BodyDecision_TOS", &K_plus_fromDs_Hlt2Topo3BodyDecision_TOS, &b_K_plus_fromDs_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyDecision_Dec", &K_plus_fromDs_Hlt2Topo4BodyDecision_Dec, &b_K_plus_fromDs_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyDecision_TIS", &K_plus_fromDs_Hlt2Topo4BodyDecision_TIS, &b_K_plus_fromDs_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_fromDs_Hlt2Topo4BodyDecision_TOS", &K_plus_fromDs_Hlt2Topo4BodyDecision_TOS, &b_K_plus_fromDs_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_Type", &K_plus_fromDs_TRACK_Type, &b_K_plus_fromDs_TRACK_Type);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_Key", &K_plus_fromDs_TRACK_Key, &b_K_plus_fromDs_TRACK_Key);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_CHI2NDOF", &K_plus_fromDs_TRACK_CHI2NDOF, &b_K_plus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_PCHI2", &K_plus_fromDs_TRACK_PCHI2, &b_K_plus_fromDs_TRACK_PCHI2);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_MatchCHI2", &K_plus_fromDs_TRACK_MatchCHI2, &b_K_plus_fromDs_TRACK_MatchCHI2);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_GhostProb", &K_plus_fromDs_TRACK_GhostProb, &b_K_plus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_CloneDist", &K_plus_fromDs_TRACK_CloneDist, &b_K_plus_fromDs_TRACK_CloneDist);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_Likelihood", &K_plus_fromDs_TRACK_Likelihood, &b_K_plus_fromDs_TRACK_Likelihood);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_0.50", &K_plus_fromDs_cpx_0_50, &b_K_plus_fromDs_cpx_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_0.50", &K_plus_fromDs_cpy_0_50, &b_K_plus_fromDs_cpy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_0.50", &K_plus_fromDs_cpz_0_50, &b_K_plus_fromDs_cpz_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_0.50", &K_plus_fromDs_cpt_0_50, &b_K_plus_fromDs_cpt_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cp_0.50", &K_plus_fromDs_cp_0_50, &b_K_plus_fromDs_cp_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_0.50", &K_plus_fromDs_cmult_0_50, &b_K_plus_fromDs_cmult_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_0.50", &K_plus_fromDs_deltaEta_0_50, &b_K_plus_fromDs_deltaEta_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_0.50", &K_plus_fromDs_deltaPhi_0_50, &b_K_plus_fromDs_deltaPhi_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_0.50", &K_plus_fromDs_pxasy_0_50, &b_K_plus_fromDs_pxasy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_0.50", &K_plus_fromDs_pyasy_0_50, &b_K_plus_fromDs_pyasy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_0.50", &K_plus_fromDs_pzasy_0_50, &b_K_plus_fromDs_pzasy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_0.50", &K_plus_fromDs_pasy_0_50, &b_K_plus_fromDs_pasy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_0.50", &K_plus_fromDs_ptasy_0_50, &b_K_plus_fromDs_ptasy_0_50);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_0.60", &K_plus_fromDs_cpx_0_60, &b_K_plus_fromDs_cpx_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_0.60", &K_plus_fromDs_cpy_0_60, &b_K_plus_fromDs_cpy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_0.60", &K_plus_fromDs_cpz_0_60, &b_K_plus_fromDs_cpz_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_0.60", &K_plus_fromDs_cpt_0_60, &b_K_plus_fromDs_cpt_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cp_0.60", &K_plus_fromDs_cp_0_60, &b_K_plus_fromDs_cp_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_0.60", &K_plus_fromDs_cmult_0_60, &b_K_plus_fromDs_cmult_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_0.60", &K_plus_fromDs_deltaEta_0_60, &b_K_plus_fromDs_deltaEta_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_0.60", &K_plus_fromDs_deltaPhi_0_60, &b_K_plus_fromDs_deltaPhi_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_0.60", &K_plus_fromDs_pxasy_0_60, &b_K_plus_fromDs_pxasy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_0.60", &K_plus_fromDs_pyasy_0_60, &b_K_plus_fromDs_pyasy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_0.60", &K_plus_fromDs_pzasy_0_60, &b_K_plus_fromDs_pzasy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_0.60", &K_plus_fromDs_pasy_0_60, &b_K_plus_fromDs_pasy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_0.60", &K_plus_fromDs_ptasy_0_60, &b_K_plus_fromDs_ptasy_0_60);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_0.70", &K_plus_fromDs_cpx_0_70, &b_K_plus_fromDs_cpx_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_0.70", &K_plus_fromDs_cpy_0_70, &b_K_plus_fromDs_cpy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_0.70", &K_plus_fromDs_cpz_0_70, &b_K_plus_fromDs_cpz_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_0.70", &K_plus_fromDs_cpt_0_70, &b_K_plus_fromDs_cpt_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cp_0.70", &K_plus_fromDs_cp_0_70, &b_K_plus_fromDs_cp_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_0.70", &K_plus_fromDs_cmult_0_70, &b_K_plus_fromDs_cmult_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_0.70", &K_plus_fromDs_deltaEta_0_70, &b_K_plus_fromDs_deltaEta_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_0.70", &K_plus_fromDs_deltaPhi_0_70, &b_K_plus_fromDs_deltaPhi_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_0.70", &K_plus_fromDs_pxasy_0_70, &b_K_plus_fromDs_pxasy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_0.70", &K_plus_fromDs_pyasy_0_70, &b_K_plus_fromDs_pyasy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_0.70", &K_plus_fromDs_pzasy_0_70, &b_K_plus_fromDs_pzasy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_0.70", &K_plus_fromDs_pasy_0_70, &b_K_plus_fromDs_pasy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_0.70", &K_plus_fromDs_ptasy_0_70, &b_K_plus_fromDs_ptasy_0_70);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_0.80", &K_plus_fromDs_cpx_0_80, &b_K_plus_fromDs_cpx_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_0.80", &K_plus_fromDs_cpy_0_80, &b_K_plus_fromDs_cpy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_0.80", &K_plus_fromDs_cpz_0_80, &b_K_plus_fromDs_cpz_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_0.80", &K_plus_fromDs_cpt_0_80, &b_K_plus_fromDs_cpt_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cp_0.80", &K_plus_fromDs_cp_0_80, &b_K_plus_fromDs_cp_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_0.80", &K_plus_fromDs_cmult_0_80, &b_K_plus_fromDs_cmult_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_0.80", &K_plus_fromDs_deltaEta_0_80, &b_K_plus_fromDs_deltaEta_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_0.80", &K_plus_fromDs_deltaPhi_0_80, &b_K_plus_fromDs_deltaPhi_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_0.80", &K_plus_fromDs_pxasy_0_80, &b_K_plus_fromDs_pxasy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_0.80", &K_plus_fromDs_pyasy_0_80, &b_K_plus_fromDs_pyasy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_0.80", &K_plus_fromDs_pzasy_0_80, &b_K_plus_fromDs_pzasy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_0.80", &K_plus_fromDs_pasy_0_80, &b_K_plus_fromDs_pasy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_0.80", &K_plus_fromDs_ptasy_0_80, &b_K_plus_fromDs_ptasy_0_80);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_0.90", &K_plus_fromDs_cpx_0_90, &b_K_plus_fromDs_cpx_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_0.90", &K_plus_fromDs_cpy_0_90, &b_K_plus_fromDs_cpy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_0.90", &K_plus_fromDs_cpz_0_90, &b_K_plus_fromDs_cpz_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_0.90", &K_plus_fromDs_cpt_0_90, &b_K_plus_fromDs_cpt_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cp_0.90", &K_plus_fromDs_cp_0_90, &b_K_plus_fromDs_cp_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_0.90", &K_plus_fromDs_cmult_0_90, &b_K_plus_fromDs_cmult_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_0.90", &K_plus_fromDs_deltaEta_0_90, &b_K_plus_fromDs_deltaEta_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_0.90", &K_plus_fromDs_deltaPhi_0_90, &b_K_plus_fromDs_deltaPhi_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_0.90", &K_plus_fromDs_pxasy_0_90, &b_K_plus_fromDs_pxasy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_0.90", &K_plus_fromDs_pyasy_0_90, &b_K_plus_fromDs_pyasy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_0.90", &K_plus_fromDs_pzasy_0_90, &b_K_plus_fromDs_pzasy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_0.90", &K_plus_fromDs_pasy_0_90, &b_K_plus_fromDs_pasy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_0.90", &K_plus_fromDs_ptasy_0_90, &b_K_plus_fromDs_ptasy_0_90);
    fChain->SetBranchAddress("K_plus_fromDs_cpx_1.00", &K_plus_fromDs_cpx_1_00, &b_K_plus_fromDs_cpx_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_cpy_1.00", &K_plus_fromDs_cpy_1_00, &b_K_plus_fromDs_cpy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_cpz_1.00", &K_plus_fromDs_cpz_1_00, &b_K_plus_fromDs_cpz_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_cpt_1.00", &K_plus_fromDs_cpt_1_00, &b_K_plus_fromDs_cpt_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_cp_1.00", &K_plus_fromDs_cp_1_00, &b_K_plus_fromDs_cp_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_cmult_1.00", &K_plus_fromDs_cmult_1_00, &b_K_plus_fromDs_cmult_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_deltaEta_1.00", &K_plus_fromDs_deltaEta_1_00, &b_K_plus_fromDs_deltaEta_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_deltaPhi_1.00", &K_plus_fromDs_deltaPhi_1_00, &b_K_plus_fromDs_deltaPhi_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_pxasy_1.00", &K_plus_fromDs_pxasy_1_00, &b_K_plus_fromDs_pxasy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_pyasy_1.00", &K_plus_fromDs_pyasy_1_00, &b_K_plus_fromDs_pyasy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_pzasy_1.00", &K_plus_fromDs_pzasy_1_00, &b_K_plus_fromDs_pzasy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_pasy_1.00", &K_plus_fromDs_pasy_1_00, &b_K_plus_fromDs_pasy_1_00);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_1.00", &K_plus_fromDs_ptasy_1_00, &b_K_plus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_DOCA1", &K_minus_fromDs_DOCA1, &b_K_minus_fromDs_DOCA1);
    fChain->SetBranchAddress("K_minus_fromDs_DOCA2", &K_minus_fromDs_DOCA2, &b_K_minus_fromDs_DOCA2);
    fChain->SetBranchAddress("K_minus_fromDs_DOCA3", &K_minus_fromDs_DOCA3, &b_K_minus_fromDs_DOCA3);
    fChain->SetBranchAddress("K_minus_fromDs_ETA", &K_minus_fromDs_ETA, &b_K_minus_fromDs_ETA);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNe", &K_minus_fromDs_MC12TuneV2_ProbNNe, &b_K_minus_fromDs_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNmu", &K_minus_fromDs_MC12TuneV2_ProbNNmu, &b_K_minus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNpi", &K_minus_fromDs_MC12TuneV2_ProbNNpi, &b_K_minus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNk", &K_minus_fromDs_MC12TuneV2_ProbNNk, &b_K_minus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNp", &K_minus_fromDs_MC12TuneV2_ProbNNp, &b_K_minus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNghost", &K_minus_fromDs_MC12TuneV2_ProbNNghost, &b_K_minus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNe", &K_minus_fromDs_MC12TuneV3_ProbNNe, &b_K_minus_fromDs_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNmu", &K_minus_fromDs_MC12TuneV3_ProbNNmu, &b_K_minus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNpi", &K_minus_fromDs_MC12TuneV3_ProbNNpi, &b_K_minus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNk", &K_minus_fromDs_MC12TuneV3_ProbNNk, &b_K_minus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNp", &K_minus_fromDs_MC12TuneV3_ProbNNp, &b_K_minus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNghost", &K_minus_fromDs_MC12TuneV3_ProbNNghost, &b_K_minus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNe", &K_minus_fromDs_MC12TuneV4_ProbNNe, &b_K_minus_fromDs_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNmu", &K_minus_fromDs_MC12TuneV4_ProbNNmu, &b_K_minus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNpi", &K_minus_fromDs_MC12TuneV4_ProbNNpi, &b_K_minus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNk", &K_minus_fromDs_MC12TuneV4_ProbNNk, &b_K_minus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNp", &K_minus_fromDs_MC12TuneV4_ProbNNp, &b_K_minus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNghost", &K_minus_fromDs_MC12TuneV4_ProbNNghost, &b_K_minus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNe", &K_minus_fromDs_MC15TuneV1_ProbNNe, &b_K_minus_fromDs_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNmu", &K_minus_fromDs_MC15TuneV1_ProbNNmu, &b_K_minus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNpi", &K_minus_fromDs_MC15TuneV1_ProbNNpi, &b_K_minus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNk", &K_minus_fromDs_MC15TuneV1_ProbNNk, &b_K_minus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNp", &K_minus_fromDs_MC15TuneV1_ProbNNp, &b_K_minus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNghost", &K_minus_fromDs_MC15TuneV1_ProbNNghost, &b_K_minus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_CosTheta", &K_minus_fromDs_CosTheta, &b_K_minus_fromDs_CosTheta);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_X", &K_minus_fromDs_OWNPV_X, &b_K_minus_fromDs_OWNPV_X);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_Y", &K_minus_fromDs_OWNPV_Y, &b_K_minus_fromDs_OWNPV_Y);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_Z", &K_minus_fromDs_OWNPV_Z, &b_K_minus_fromDs_OWNPV_Z);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_XERR", &K_minus_fromDs_OWNPV_XERR, &b_K_minus_fromDs_OWNPV_XERR);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_YERR", &K_minus_fromDs_OWNPV_YERR, &b_K_minus_fromDs_OWNPV_YERR);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_ZERR", &K_minus_fromDs_OWNPV_ZERR, &b_K_minus_fromDs_OWNPV_ZERR);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_CHI2", &K_minus_fromDs_OWNPV_CHI2, &b_K_minus_fromDs_OWNPV_CHI2);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_NDOF", &K_minus_fromDs_OWNPV_NDOF, &b_K_minus_fromDs_OWNPV_NDOF);
    fChain->SetBranchAddress("K_minus_fromDs_OWNPV_COV_", K_minus_fromDs_OWNPV_COV_, &b_K_minus_fromDs_OWNPV_COV_);
    fChain->SetBranchAddress("K_minus_fromDs_IP_OWNPV", &K_minus_fromDs_IP_OWNPV, &b_K_minus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("K_minus_fromDs_IPCHI2_OWNPV", &K_minus_fromDs_IPCHI2_OWNPV, &b_K_minus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_X", &K_minus_fromDs_ORIVX_X, &b_K_minus_fromDs_ORIVX_X);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_Y", &K_minus_fromDs_ORIVX_Y, &b_K_minus_fromDs_ORIVX_Y);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_Z", &K_minus_fromDs_ORIVX_Z, &b_K_minus_fromDs_ORIVX_Z);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_XERR", &K_minus_fromDs_ORIVX_XERR, &b_K_minus_fromDs_ORIVX_XERR);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_YERR", &K_minus_fromDs_ORIVX_YERR, &b_K_minus_fromDs_ORIVX_YERR);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_ZERR", &K_minus_fromDs_ORIVX_ZERR, &b_K_minus_fromDs_ORIVX_ZERR);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_CHI2", &K_minus_fromDs_ORIVX_CHI2, &b_K_minus_fromDs_ORIVX_CHI2);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_NDOF", &K_minus_fromDs_ORIVX_NDOF, &b_K_minus_fromDs_ORIVX_NDOF);
    fChain->SetBranchAddress("K_minus_fromDs_ORIVX_COV_", K_minus_fromDs_ORIVX_COV_, &b_K_minus_fromDs_ORIVX_COV_);
    fChain->SetBranchAddress("K_minus_fromDs_P", &K_minus_fromDs_P, &b_K_minus_fromDs_P);
    fChain->SetBranchAddress("K_minus_fromDs_PT", &K_minus_fromDs_PT, &b_K_minus_fromDs_PT);
    fChain->SetBranchAddress("K_minus_fromDs_PE", &K_minus_fromDs_PE, &b_K_minus_fromDs_PE);
    fChain->SetBranchAddress("K_minus_fromDs_PX", &K_minus_fromDs_PX, &b_K_minus_fromDs_PX);
    fChain->SetBranchAddress("K_minus_fromDs_PY", &K_minus_fromDs_PY, &b_K_minus_fromDs_PY);
    fChain->SetBranchAddress("K_minus_fromDs_PZ", &K_minus_fromDs_PZ, &b_K_minus_fromDs_PZ);
    fChain->SetBranchAddress("K_minus_fromDs_M", &K_minus_fromDs_M, &b_K_minus_fromDs_M);
    fChain->SetBranchAddress("K_minus_fromDs_ID", &K_minus_fromDs_ID, &b_K_minus_fromDs_ID);
    fChain->SetBranchAddress("K_minus_fromDs_PIDe", &K_minus_fromDs_PIDe, &b_K_minus_fromDs_PIDe);
    fChain->SetBranchAddress("K_minus_fromDs_PIDmu", &K_minus_fromDs_PIDmu, &b_K_minus_fromDs_PIDmu);
    fChain->SetBranchAddress("K_minus_fromDs_PIDK", &K_minus_fromDs_PIDK, &b_K_minus_fromDs_PIDK);
    fChain->SetBranchAddress("K_minus_fromDs_PIDp", &K_minus_fromDs_PIDp, &b_K_minus_fromDs_PIDp);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNe", &K_minus_fromDs_ProbNNe, &b_K_minus_fromDs_ProbNNe);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNk", &K_minus_fromDs_ProbNNk, &b_K_minus_fromDs_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNp", &K_minus_fromDs_ProbNNp, &b_K_minus_fromDs_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNpi", &K_minus_fromDs_ProbNNpi, &b_K_minus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNmu", &K_minus_fromDs_ProbNNmu, &b_K_minus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNghost", &K_minus_fromDs_ProbNNghost, &b_K_minus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_hasMuon", &K_minus_fromDs_hasMuon, &b_K_minus_fromDs_hasMuon);
    fChain->SetBranchAddress("K_minus_fromDs_isMuon", &K_minus_fromDs_isMuon, &b_K_minus_fromDs_isMuon);
    fChain->SetBranchAddress("K_minus_fromDs_hasRich", &K_minus_fromDs_hasRich, &b_K_minus_fromDs_hasRich);
    fChain->SetBranchAddress("K_minus_fromDs_UsedRichAerogel", &K_minus_fromDs_UsedRichAerogel, &b_K_minus_fromDs_UsedRichAerogel);
    fChain->SetBranchAddress("K_minus_fromDs_UsedRich1Gas", &K_minus_fromDs_UsedRich1Gas, &b_K_minus_fromDs_UsedRich1Gas);
    fChain->SetBranchAddress("K_minus_fromDs_UsedRich2Gas", &K_minus_fromDs_UsedRich2Gas, &b_K_minus_fromDs_UsedRich2Gas);
    fChain->SetBranchAddress("K_minus_fromDs_RichAboveElThres", &K_minus_fromDs_RichAboveElThres, &b_K_minus_fromDs_RichAboveElThres);
    fChain->SetBranchAddress("K_minus_fromDs_RichAboveMuThres", &K_minus_fromDs_RichAboveMuThres, &b_K_minus_fromDs_RichAboveMuThres);
    fChain->SetBranchAddress("K_minus_fromDs_RichAbovePiThres", &K_minus_fromDs_RichAbovePiThres, &b_K_minus_fromDs_RichAbovePiThres);
    fChain->SetBranchAddress("K_minus_fromDs_RichAboveKaThres", &K_minus_fromDs_RichAboveKaThres, &b_K_minus_fromDs_RichAboveKaThres);
    fChain->SetBranchAddress("K_minus_fromDs_RichAbovePrThres", &K_minus_fromDs_RichAbovePrThres, &b_K_minus_fromDs_RichAbovePrThres);
    fChain->SetBranchAddress("K_minus_fromDs_hasCalo", &K_minus_fromDs_hasCalo, &b_K_minus_fromDs_hasCalo);
    fChain->SetBranchAddress("K_minus_fromDs_L0Global_Dec", &K_minus_fromDs_L0Global_Dec, &b_K_minus_fromDs_L0Global_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_L0Global_TIS", &K_minus_fromDs_L0Global_TIS, &b_K_minus_fromDs_L0Global_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_L0Global_TOS", &K_minus_fromDs_L0Global_TOS, &b_K_minus_fromDs_L0Global_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Global_Dec", &K_minus_fromDs_Hlt1Global_Dec, &b_K_minus_fromDs_Hlt1Global_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Global_TIS", &K_minus_fromDs_Hlt1Global_TIS, &b_K_minus_fromDs_Hlt1Global_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Global_TOS", &K_minus_fromDs_Hlt1Global_TOS, &b_K_minus_fromDs_Hlt1Global_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Phys_Dec", &K_minus_fromDs_Hlt1Phys_Dec, &b_K_minus_fromDs_Hlt1Phys_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Phys_TIS", &K_minus_fromDs_Hlt1Phys_TIS, &b_K_minus_fromDs_Hlt1Phys_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1Phys_TOS", &K_minus_fromDs_Hlt1Phys_TOS, &b_K_minus_fromDs_Hlt1Phys_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Global_Dec", &K_minus_fromDs_Hlt2Global_Dec, &b_K_minus_fromDs_Hlt2Global_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Global_TIS", &K_minus_fromDs_Hlt2Global_TIS, &b_K_minus_fromDs_Hlt2Global_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Global_TOS", &K_minus_fromDs_Hlt2Global_TOS, &b_K_minus_fromDs_Hlt2Global_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Phys_Dec", &K_minus_fromDs_Hlt2Phys_Dec, &b_K_minus_fromDs_Hlt2Phys_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Phys_TIS", &K_minus_fromDs_Hlt2Phys_TIS, &b_K_minus_fromDs_Hlt2Phys_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Phys_TOS", &K_minus_fromDs_Hlt2Phys_TOS, &b_K_minus_fromDs_Hlt2Phys_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_L0HadronDecision_Dec", &K_minus_fromDs_L0HadronDecision_Dec, &b_K_minus_fromDs_L0HadronDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_L0HadronDecision_TIS", &K_minus_fromDs_L0HadronDecision_TIS, &b_K_minus_fromDs_L0HadronDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_L0HadronDecision_TOS", &K_minus_fromDs_L0HadronDecision_TOS, &b_K_minus_fromDs_L0HadronDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_L0MuonDecision_Dec", &K_minus_fromDs_L0MuonDecision_Dec, &b_K_minus_fromDs_L0MuonDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_L0MuonDecision_TIS", &K_minus_fromDs_L0MuonDecision_TIS, &b_K_minus_fromDs_L0MuonDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_L0MuonDecision_TOS", &K_minus_fromDs_L0MuonDecision_TOS, &b_K_minus_fromDs_L0MuonDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_L0GlobalDecision_Dec", &K_minus_fromDs_L0GlobalDecision_Dec, &b_K_minus_fromDs_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_L0GlobalDecision_TIS", &K_minus_fromDs_L0GlobalDecision_TIS, &b_K_minus_fromDs_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_L0GlobalDecision_TOS", &K_minus_fromDs_L0GlobalDecision_TOS, &b_K_minus_fromDs_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackAllL0Decision_Dec", &K_minus_fromDs_Hlt1TrackAllL0Decision_Dec, &b_K_minus_fromDs_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackAllL0Decision_TIS", &K_minus_fromDs_Hlt1TrackAllL0Decision_TIS, &b_K_minus_fromDs_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackAllL0Decision_TOS", &K_minus_fromDs_Hlt1TrackAllL0Decision_TOS, &b_K_minus_fromDs_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVADecision_Dec", &K_minus_fromDs_Hlt1TrackMVADecision_Dec, &b_K_minus_fromDs_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVADecision_TIS", &K_minus_fromDs_Hlt1TrackMVADecision_TIS, &b_K_minus_fromDs_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVADecision_TOS", &K_minus_fromDs_Hlt1TrackMVADecision_TOS, &b_K_minus_fromDs_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVADecision_Dec", &K_minus_fromDs_Hlt1TwoTrackMVADecision_Dec, &b_K_minus_fromDs_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVADecision_TIS", &K_minus_fromDs_Hlt1TwoTrackMVADecision_TIS, &b_K_minus_fromDs_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVADecision_TOS", &K_minus_fromDs_Hlt1TwoTrackMVADecision_TOS, &b_K_minus_fromDs_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVALooseDecision_Dec", &K_minus_fromDs_Hlt1TrackMVALooseDecision_Dec, &b_K_minus_fromDs_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVALooseDecision_TIS", &K_minus_fromDs_Hlt1TrackMVALooseDecision_TIS, &b_K_minus_fromDs_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TrackMVALooseDecision_TOS", &K_minus_fromDs_Hlt1TrackMVALooseDecision_TOS, &b_K_minus_fromDs_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec", &K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec, &b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS", &K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS, &b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS", &K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS, &b_K_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2IncPhiDecision_Dec", &K_minus_fromDs_Hlt2IncPhiDecision_Dec, &b_K_minus_fromDs_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2IncPhiDecision_TIS", &K_minus_fromDs_Hlt2IncPhiDecision_TIS, &b_K_minus_fromDs_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2IncPhiDecision_TOS", &K_minus_fromDs_Hlt2IncPhiDecision_TOS, &b_K_minus_fromDs_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2PhiIncPhiDecision_Dec", &K_minus_fromDs_Hlt2PhiIncPhiDecision_Dec, &b_K_minus_fromDs_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2PhiIncPhiDecision_TIS", &K_minus_fromDs_Hlt2PhiIncPhiDecision_TIS, &b_K_minus_fromDs_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2PhiIncPhiDecision_TOS", &K_minus_fromDs_Hlt2PhiIncPhiDecision_TOS, &b_K_minus_fromDs_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec", &K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec, &b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS", &K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS, &b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS", &K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS, &b_K_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec", &K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec, &b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS", &K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS, &b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS", &K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS, &b_K_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec", &K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec, &b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS", &K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS, &b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS", &K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS, &b_K_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyDecision_Dec", &K_minus_fromDs_Hlt2Topo2BodyDecision_Dec, &b_K_minus_fromDs_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyDecision_TIS", &K_minus_fromDs_Hlt2Topo2BodyDecision_TIS, &b_K_minus_fromDs_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo2BodyDecision_TOS", &K_minus_fromDs_Hlt2Topo2BodyDecision_TOS, &b_K_minus_fromDs_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyDecision_Dec", &K_minus_fromDs_Hlt2Topo3BodyDecision_Dec, &b_K_minus_fromDs_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyDecision_TIS", &K_minus_fromDs_Hlt2Topo3BodyDecision_TIS, &b_K_minus_fromDs_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo3BodyDecision_TOS", &K_minus_fromDs_Hlt2Topo3BodyDecision_TOS, &b_K_minus_fromDs_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyDecision_Dec", &K_minus_fromDs_Hlt2Topo4BodyDecision_Dec, &b_K_minus_fromDs_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyDecision_TIS", &K_minus_fromDs_Hlt2Topo4BodyDecision_TIS, &b_K_minus_fromDs_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("K_minus_fromDs_Hlt2Topo4BodyDecision_TOS", &K_minus_fromDs_Hlt2Topo4BodyDecision_TOS, &b_K_minus_fromDs_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_Type", &K_minus_fromDs_TRACK_Type, &b_K_minus_fromDs_TRACK_Type);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_Key", &K_minus_fromDs_TRACK_Key, &b_K_minus_fromDs_TRACK_Key);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_CHI2NDOF", &K_minus_fromDs_TRACK_CHI2NDOF, &b_K_minus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_PCHI2", &K_minus_fromDs_TRACK_PCHI2, &b_K_minus_fromDs_TRACK_PCHI2);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_MatchCHI2", &K_minus_fromDs_TRACK_MatchCHI2, &b_K_minus_fromDs_TRACK_MatchCHI2);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb, &b_K_minus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_CloneDist", &K_minus_fromDs_TRACK_CloneDist, &b_K_minus_fromDs_TRACK_CloneDist);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_Likelihood", &K_minus_fromDs_TRACK_Likelihood, &b_K_minus_fromDs_TRACK_Likelihood);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_0.50", &K_minus_fromDs_cpx_0_50, &b_K_minus_fromDs_cpx_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_0.50", &K_minus_fromDs_cpy_0_50, &b_K_minus_fromDs_cpy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_0.50", &K_minus_fromDs_cpz_0_50, &b_K_minus_fromDs_cpz_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_0.50", &K_minus_fromDs_cpt_0_50, &b_K_minus_fromDs_cpt_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cp_0.50", &K_minus_fromDs_cp_0_50, &b_K_minus_fromDs_cp_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_0.50", &K_minus_fromDs_cmult_0_50, &b_K_minus_fromDs_cmult_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_0.50", &K_minus_fromDs_deltaEta_0_50, &b_K_minus_fromDs_deltaEta_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_0.50", &K_minus_fromDs_deltaPhi_0_50, &b_K_minus_fromDs_deltaPhi_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_0.50", &K_minus_fromDs_pxasy_0_50, &b_K_minus_fromDs_pxasy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_0.50", &K_minus_fromDs_pyasy_0_50, &b_K_minus_fromDs_pyasy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_0.50", &K_minus_fromDs_pzasy_0_50, &b_K_minus_fromDs_pzasy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_0.50", &K_minus_fromDs_pasy_0_50, &b_K_minus_fromDs_pasy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_0.50", &K_minus_fromDs_ptasy_0_50, &b_K_minus_fromDs_ptasy_0_50);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_0.60", &K_minus_fromDs_cpx_0_60, &b_K_minus_fromDs_cpx_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_0.60", &K_minus_fromDs_cpy_0_60, &b_K_minus_fromDs_cpy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_0.60", &K_minus_fromDs_cpz_0_60, &b_K_minus_fromDs_cpz_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_0.60", &K_minus_fromDs_cpt_0_60, &b_K_minus_fromDs_cpt_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cp_0.60", &K_minus_fromDs_cp_0_60, &b_K_minus_fromDs_cp_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_0.60", &K_minus_fromDs_cmult_0_60, &b_K_minus_fromDs_cmult_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_0.60", &K_minus_fromDs_deltaEta_0_60, &b_K_minus_fromDs_deltaEta_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_0.60", &K_minus_fromDs_deltaPhi_0_60, &b_K_minus_fromDs_deltaPhi_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_0.60", &K_minus_fromDs_pxasy_0_60, &b_K_minus_fromDs_pxasy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_0.60", &K_minus_fromDs_pyasy_0_60, &b_K_minus_fromDs_pyasy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_0.60", &K_minus_fromDs_pzasy_0_60, &b_K_minus_fromDs_pzasy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_0.60", &K_minus_fromDs_pasy_0_60, &b_K_minus_fromDs_pasy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_0.60", &K_minus_fromDs_ptasy_0_60, &b_K_minus_fromDs_ptasy_0_60);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_0.70", &K_minus_fromDs_cpx_0_70, &b_K_minus_fromDs_cpx_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_0.70", &K_minus_fromDs_cpy_0_70, &b_K_minus_fromDs_cpy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_0.70", &K_minus_fromDs_cpz_0_70, &b_K_minus_fromDs_cpz_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_0.70", &K_minus_fromDs_cpt_0_70, &b_K_minus_fromDs_cpt_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cp_0.70", &K_minus_fromDs_cp_0_70, &b_K_minus_fromDs_cp_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_0.70", &K_minus_fromDs_cmult_0_70, &b_K_minus_fromDs_cmult_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_0.70", &K_minus_fromDs_deltaEta_0_70, &b_K_minus_fromDs_deltaEta_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_0.70", &K_minus_fromDs_deltaPhi_0_70, &b_K_minus_fromDs_deltaPhi_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_0.70", &K_minus_fromDs_pxasy_0_70, &b_K_minus_fromDs_pxasy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_0.70", &K_minus_fromDs_pyasy_0_70, &b_K_minus_fromDs_pyasy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_0.70", &K_minus_fromDs_pzasy_0_70, &b_K_minus_fromDs_pzasy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_0.70", &K_minus_fromDs_pasy_0_70, &b_K_minus_fromDs_pasy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_0.70", &K_minus_fromDs_ptasy_0_70, &b_K_minus_fromDs_ptasy_0_70);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_0.80", &K_minus_fromDs_cpx_0_80, &b_K_minus_fromDs_cpx_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_0.80", &K_minus_fromDs_cpy_0_80, &b_K_minus_fromDs_cpy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_0.80", &K_minus_fromDs_cpz_0_80, &b_K_minus_fromDs_cpz_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_0.80", &K_minus_fromDs_cpt_0_80, &b_K_minus_fromDs_cpt_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cp_0.80", &K_minus_fromDs_cp_0_80, &b_K_minus_fromDs_cp_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_0.80", &K_minus_fromDs_cmult_0_80, &b_K_minus_fromDs_cmult_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_0.80", &K_minus_fromDs_deltaEta_0_80, &b_K_minus_fromDs_deltaEta_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_0.80", &K_minus_fromDs_deltaPhi_0_80, &b_K_minus_fromDs_deltaPhi_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_0.80", &K_minus_fromDs_pxasy_0_80, &b_K_minus_fromDs_pxasy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_0.80", &K_minus_fromDs_pyasy_0_80, &b_K_minus_fromDs_pyasy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_0.80", &K_minus_fromDs_pzasy_0_80, &b_K_minus_fromDs_pzasy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_0.80", &K_minus_fromDs_pasy_0_80, &b_K_minus_fromDs_pasy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_0.80", &K_minus_fromDs_ptasy_0_80, &b_K_minus_fromDs_ptasy_0_80);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_0.90", &K_minus_fromDs_cpx_0_90, &b_K_minus_fromDs_cpx_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_0.90", &K_minus_fromDs_cpy_0_90, &b_K_minus_fromDs_cpy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_0.90", &K_minus_fromDs_cpz_0_90, &b_K_minus_fromDs_cpz_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_0.90", &K_minus_fromDs_cpt_0_90, &b_K_minus_fromDs_cpt_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cp_0.90", &K_minus_fromDs_cp_0_90, &b_K_minus_fromDs_cp_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_0.90", &K_minus_fromDs_cmult_0_90, &b_K_minus_fromDs_cmult_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_0.90", &K_minus_fromDs_deltaEta_0_90, &b_K_minus_fromDs_deltaEta_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_0.90", &K_minus_fromDs_deltaPhi_0_90, &b_K_minus_fromDs_deltaPhi_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_0.90", &K_minus_fromDs_pxasy_0_90, &b_K_minus_fromDs_pxasy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_0.90", &K_minus_fromDs_pyasy_0_90, &b_K_minus_fromDs_pyasy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_0.90", &K_minus_fromDs_pzasy_0_90, &b_K_minus_fromDs_pzasy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_0.90", &K_minus_fromDs_pasy_0_90, &b_K_minus_fromDs_pasy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_0.90", &K_minus_fromDs_ptasy_0_90, &b_K_minus_fromDs_ptasy_0_90);
    fChain->SetBranchAddress("K_minus_fromDs_cpx_1.00", &K_minus_fromDs_cpx_1_00, &b_K_minus_fromDs_cpx_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_cpy_1.00", &K_minus_fromDs_cpy_1_00, &b_K_minus_fromDs_cpy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_cpz_1.00", &K_minus_fromDs_cpz_1_00, &b_K_minus_fromDs_cpz_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_cpt_1.00", &K_minus_fromDs_cpt_1_00, &b_K_minus_fromDs_cpt_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_cp_1.00", &K_minus_fromDs_cp_1_00, &b_K_minus_fromDs_cp_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_cmult_1.00", &K_minus_fromDs_cmult_1_00, &b_K_minus_fromDs_cmult_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_deltaEta_1.00", &K_minus_fromDs_deltaEta_1_00, &b_K_minus_fromDs_deltaEta_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_deltaPhi_1.00", &K_minus_fromDs_deltaPhi_1_00, &b_K_minus_fromDs_deltaPhi_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_pxasy_1.00", &K_minus_fromDs_pxasy_1_00, &b_K_minus_fromDs_pxasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_pyasy_1.00", &K_minus_fromDs_pyasy_1_00, &b_K_minus_fromDs_pyasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_pzasy_1.00", &K_minus_fromDs_pzasy_1_00, &b_K_minus_fromDs_pzasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_pasy_1.00", &K_minus_fromDs_pasy_1_00, &b_K_minus_fromDs_pasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_1.00", &K_minus_fromDs_ptasy_1_00, &b_K_minus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_DOCA1", &pi_minus_fromDs_DOCA1, &b_pi_minus_fromDs_DOCA1);
    fChain->SetBranchAddress("pi_minus_fromDs_DOCA2", &pi_minus_fromDs_DOCA2, &b_pi_minus_fromDs_DOCA2);
    fChain->SetBranchAddress("pi_minus_fromDs_DOCA3", &pi_minus_fromDs_DOCA3, &b_pi_minus_fromDs_DOCA3);
    fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNe", &pi_minus_fromDs_MC12TuneV2_ProbNNe, &b_pi_minus_fromDs_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNe", &pi_minus_fromDs_MC12TuneV3_ProbNNe, &b_pi_minus_fromDs_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNe", &pi_minus_fromDs_MC12TuneV4_ProbNNe, &b_pi_minus_fromDs_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNmu", &pi_minus_fromDs_MC12TuneV4_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNpi", &pi_minus_fromDs_MC12TuneV4_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNk", &pi_minus_fromDs_MC12TuneV4_ProbNNk, &b_pi_minus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNp", &pi_minus_fromDs_MC12TuneV4_ProbNNp, &b_pi_minus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNghost", &pi_minus_fromDs_MC12TuneV4_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNe", &pi_minus_fromDs_MC15TuneV1_ProbNNe, &b_pi_minus_fromDs_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNmu", &pi_minus_fromDs_MC15TuneV1_ProbNNmu, &b_pi_minus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNpi", &pi_minus_fromDs_MC15TuneV1_ProbNNpi, &b_pi_minus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNk", &pi_minus_fromDs_MC15TuneV1_ProbNNk, &b_pi_minus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNp", &pi_minus_fromDs_MC15TuneV1_ProbNNp, &b_pi_minus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNghost", &pi_minus_fromDs_MC15TuneV1_ProbNNghost, &b_pi_minus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_CosTheta", &pi_minus_fromDs_CosTheta, &b_pi_minus_fromDs_CosTheta);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_X", &pi_minus_fromDs_OWNPV_X, &b_pi_minus_fromDs_OWNPV_X);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_Y", &pi_minus_fromDs_OWNPV_Y, &b_pi_minus_fromDs_OWNPV_Y);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_Z", &pi_minus_fromDs_OWNPV_Z, &b_pi_minus_fromDs_OWNPV_Z);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_XERR", &pi_minus_fromDs_OWNPV_XERR, &b_pi_minus_fromDs_OWNPV_XERR);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_YERR", &pi_minus_fromDs_OWNPV_YERR, &b_pi_minus_fromDs_OWNPV_YERR);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_ZERR", &pi_minus_fromDs_OWNPV_ZERR, &b_pi_minus_fromDs_OWNPV_ZERR);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_CHI2", &pi_minus_fromDs_OWNPV_CHI2, &b_pi_minus_fromDs_OWNPV_CHI2);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_NDOF", &pi_minus_fromDs_OWNPV_NDOF, &b_pi_minus_fromDs_OWNPV_NDOF);
    fChain->SetBranchAddress("pi_minus_fromDs_OWNPV_COV_", pi_minus_fromDs_OWNPV_COV_, &b_pi_minus_fromDs_OWNPV_COV_);
    fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_X", &pi_minus_fromDs_ORIVX_X, &b_pi_minus_fromDs_ORIVX_X);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_Y", &pi_minus_fromDs_ORIVX_Y, &b_pi_minus_fromDs_ORIVX_Y);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_Z", &pi_minus_fromDs_ORIVX_Z, &b_pi_minus_fromDs_ORIVX_Z);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_XERR", &pi_minus_fromDs_ORIVX_XERR, &b_pi_minus_fromDs_ORIVX_XERR);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_YERR", &pi_minus_fromDs_ORIVX_YERR, &b_pi_minus_fromDs_ORIVX_YERR);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_ZERR", &pi_minus_fromDs_ORIVX_ZERR, &b_pi_minus_fromDs_ORIVX_ZERR);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_CHI2", &pi_minus_fromDs_ORIVX_CHI2, &b_pi_minus_fromDs_ORIVX_CHI2);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_NDOF", &pi_minus_fromDs_ORIVX_NDOF, &b_pi_minus_fromDs_ORIVX_NDOF);
    fChain->SetBranchAddress("pi_minus_fromDs_ORIVX_COV_", pi_minus_fromDs_ORIVX_COV_, &b_pi_minus_fromDs_ORIVX_COV_);
    fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
    fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
    fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
    fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
    fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
    fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
    fChain->SetBranchAddress("pi_minus_fromDs_M", &pi_minus_fromDs_M, &b_pi_minus_fromDs_M);
    fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDe", &pi_minus_fromDs_PIDe, &b_pi_minus_fromDs_PIDe);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNe", &pi_minus_fromDs_ProbNNe, &b_pi_minus_fromDs_ProbNNe);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_hasMuon", &pi_minus_fromDs_hasMuon, &b_pi_minus_fromDs_hasMuon);
    fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
    fChain->SetBranchAddress("pi_minus_fromDs_hasRich", &pi_minus_fromDs_hasRich, &b_pi_minus_fromDs_hasRich);
    fChain->SetBranchAddress("pi_minus_fromDs_UsedRichAerogel", &pi_minus_fromDs_UsedRichAerogel, &b_pi_minus_fromDs_UsedRichAerogel);
    fChain->SetBranchAddress("pi_minus_fromDs_UsedRich1Gas", &pi_minus_fromDs_UsedRich1Gas, &b_pi_minus_fromDs_UsedRich1Gas);
    fChain->SetBranchAddress("pi_minus_fromDs_UsedRich2Gas", &pi_minus_fromDs_UsedRich2Gas, &b_pi_minus_fromDs_UsedRich2Gas);
    fChain->SetBranchAddress("pi_minus_fromDs_RichAboveElThres", &pi_minus_fromDs_RichAboveElThres, &b_pi_minus_fromDs_RichAboveElThres);
    fChain->SetBranchAddress("pi_minus_fromDs_RichAboveMuThres", &pi_minus_fromDs_RichAboveMuThres, &b_pi_minus_fromDs_RichAboveMuThres);
    fChain->SetBranchAddress("pi_minus_fromDs_RichAbovePiThres", &pi_minus_fromDs_RichAbovePiThres, &b_pi_minus_fromDs_RichAbovePiThres);
    fChain->SetBranchAddress("pi_minus_fromDs_RichAboveKaThres", &pi_minus_fromDs_RichAboveKaThres, &b_pi_minus_fromDs_RichAboveKaThres);
    fChain->SetBranchAddress("pi_minus_fromDs_RichAbovePrThres", &pi_minus_fromDs_RichAbovePrThres, &b_pi_minus_fromDs_RichAbovePrThres);
    fChain->SetBranchAddress("pi_minus_fromDs_hasCalo", &pi_minus_fromDs_hasCalo, &b_pi_minus_fromDs_hasCalo);
    fChain->SetBranchAddress("pi_minus_fromDs_L0Global_Dec", &pi_minus_fromDs_L0Global_Dec, &b_pi_minus_fromDs_L0Global_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_L0Global_TIS", &pi_minus_fromDs_L0Global_TIS, &b_pi_minus_fromDs_L0Global_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0Global_TOS", &pi_minus_fromDs_L0Global_TOS, &b_pi_minus_fromDs_L0Global_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Global_Dec", &pi_minus_fromDs_Hlt1Global_Dec, &b_pi_minus_fromDs_Hlt1Global_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Global_TIS", &pi_minus_fromDs_Hlt1Global_TIS, &b_pi_minus_fromDs_Hlt1Global_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Global_TOS", &pi_minus_fromDs_Hlt1Global_TOS, &b_pi_minus_fromDs_Hlt1Global_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Phys_Dec", &pi_minus_fromDs_Hlt1Phys_Dec, &b_pi_minus_fromDs_Hlt1Phys_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Phys_TIS", &pi_minus_fromDs_Hlt1Phys_TIS, &b_pi_minus_fromDs_Hlt1Phys_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1Phys_TOS", &pi_minus_fromDs_Hlt1Phys_TOS, &b_pi_minus_fromDs_Hlt1Phys_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Global_Dec", &pi_minus_fromDs_Hlt2Global_Dec, &b_pi_minus_fromDs_Hlt2Global_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Global_TIS", &pi_minus_fromDs_Hlt2Global_TIS, &b_pi_minus_fromDs_Hlt2Global_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Global_TOS", &pi_minus_fromDs_Hlt2Global_TOS, &b_pi_minus_fromDs_Hlt2Global_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Phys_Dec", &pi_minus_fromDs_Hlt2Phys_Dec, &b_pi_minus_fromDs_Hlt2Phys_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Phys_TIS", &pi_minus_fromDs_Hlt2Phys_TIS, &b_pi_minus_fromDs_Hlt2Phys_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Phys_TOS", &pi_minus_fromDs_Hlt2Phys_TOS, &b_pi_minus_fromDs_Hlt2Phys_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0HadronDecision_Dec", &pi_minus_fromDs_L0HadronDecision_Dec, &b_pi_minus_fromDs_L0HadronDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_L0HadronDecision_TIS", &pi_minus_fromDs_L0HadronDecision_TIS, &b_pi_minus_fromDs_L0HadronDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0HadronDecision_TOS", &pi_minus_fromDs_L0HadronDecision_TOS, &b_pi_minus_fromDs_L0HadronDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0MuonDecision_Dec", &pi_minus_fromDs_L0MuonDecision_Dec, &b_pi_minus_fromDs_L0MuonDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_L0MuonDecision_TIS", &pi_minus_fromDs_L0MuonDecision_TIS, &b_pi_minus_fromDs_L0MuonDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0MuonDecision_TOS", &pi_minus_fromDs_L0MuonDecision_TOS, &b_pi_minus_fromDs_L0MuonDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0GlobalDecision_Dec", &pi_minus_fromDs_L0GlobalDecision_Dec, &b_pi_minus_fromDs_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_L0GlobalDecision_TIS", &pi_minus_fromDs_L0GlobalDecision_TIS, &b_pi_minus_fromDs_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_L0GlobalDecision_TOS", &pi_minus_fromDs_L0GlobalDecision_TOS, &b_pi_minus_fromDs_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackAllL0Decision_Dec", &pi_minus_fromDs_Hlt1TrackAllL0Decision_Dec, &b_pi_minus_fromDs_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackAllL0Decision_TIS", &pi_minus_fromDs_Hlt1TrackAllL0Decision_TIS, &b_pi_minus_fromDs_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackAllL0Decision_TOS", &pi_minus_fromDs_Hlt1TrackAllL0Decision_TOS, &b_pi_minus_fromDs_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVADecision_Dec", &pi_minus_fromDs_Hlt1TrackMVADecision_Dec, &b_pi_minus_fromDs_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVADecision_TIS", &pi_minus_fromDs_Hlt1TrackMVADecision_TIS, &b_pi_minus_fromDs_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVADecision_TOS", &pi_minus_fromDs_Hlt1TrackMVADecision_TOS, &b_pi_minus_fromDs_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVADecision_Dec", &pi_minus_fromDs_Hlt1TwoTrackMVADecision_Dec, &b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVADecision_TIS", &pi_minus_fromDs_Hlt1TwoTrackMVADecision_TIS, &b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVADecision_TOS", &pi_minus_fromDs_Hlt1TwoTrackMVADecision_TOS, &b_pi_minus_fromDs_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVALooseDecision_Dec", &pi_minus_fromDs_Hlt1TrackMVALooseDecision_Dec, &b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVALooseDecision_TIS", &pi_minus_fromDs_Hlt1TrackMVALooseDecision_TIS, &b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TrackMVALooseDecision_TOS", &pi_minus_fromDs_Hlt1TrackMVALooseDecision_TOS, &b_pi_minus_fromDs_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec", &pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec, &b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS", &pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS, &b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS", &pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS, &b_pi_minus_fromDs_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2IncPhiDecision_Dec", &pi_minus_fromDs_Hlt2IncPhiDecision_Dec, &b_pi_minus_fromDs_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2IncPhiDecision_TIS", &pi_minus_fromDs_Hlt2IncPhiDecision_TIS, &b_pi_minus_fromDs_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2IncPhiDecision_TOS", &pi_minus_fromDs_Hlt2IncPhiDecision_TOS, &b_pi_minus_fromDs_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2PhiIncPhiDecision_Dec", &pi_minus_fromDs_Hlt2PhiIncPhiDecision_Dec, &b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2PhiIncPhiDecision_TIS", &pi_minus_fromDs_Hlt2PhiIncPhiDecision_TIS, &b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2PhiIncPhiDecision_TOS", &pi_minus_fromDs_Hlt2PhiIncPhiDecision_TOS, &b_pi_minus_fromDs_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec", &pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS", &pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS", &pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec", &pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS", &pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS", &pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec", &pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS", &pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS", &pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyDecision_Dec", &pi_minus_fromDs_Hlt2Topo2BodyDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyDecision_TIS", &pi_minus_fromDs_Hlt2Topo2BodyDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo2BodyDecision_TOS", &pi_minus_fromDs_Hlt2Topo2BodyDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyDecision_Dec", &pi_minus_fromDs_Hlt2Topo3BodyDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyDecision_TIS", &pi_minus_fromDs_Hlt2Topo3BodyDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo3BodyDecision_TOS", &pi_minus_fromDs_Hlt2Topo3BodyDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyDecision_Dec", &pi_minus_fromDs_Hlt2Topo4BodyDecision_Dec, &b_pi_minus_fromDs_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyDecision_TIS", &pi_minus_fromDs_Hlt2Topo4BodyDecision_TIS, &b_pi_minus_fromDs_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_fromDs_Hlt2Topo4BodyDecision_TOS", &pi_minus_fromDs_Hlt2Topo4BodyDecision_TOS, &b_pi_minus_fromDs_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_Type", &pi_minus_fromDs_TRACK_Type, &b_pi_minus_fromDs_TRACK_Type);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_Key", &pi_minus_fromDs_TRACK_Key, &b_pi_minus_fromDs_TRACK_Key);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CHI2NDOF", &pi_minus_fromDs_TRACK_CHI2NDOF, &b_pi_minus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_PCHI2", &pi_minus_fromDs_TRACK_PCHI2, &b_pi_minus_fromDs_TRACK_PCHI2);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_MatchCHI2", &pi_minus_fromDs_TRACK_MatchCHI2, &b_pi_minus_fromDs_TRACK_MatchCHI2);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CloneDist", &pi_minus_fromDs_TRACK_CloneDist, &b_pi_minus_fromDs_TRACK_CloneDist);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_Likelihood", &pi_minus_fromDs_TRACK_Likelihood, &b_pi_minus_fromDs_TRACK_Likelihood);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_0.50", &pi_minus_fromDs_cpx_0_50, &b_pi_minus_fromDs_cpx_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_0.50", &pi_minus_fromDs_cpy_0_50, &b_pi_minus_fromDs_cpy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_0.50", &pi_minus_fromDs_cpz_0_50, &b_pi_minus_fromDs_cpz_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_0.50", &pi_minus_fromDs_cpt_0_50, &b_pi_minus_fromDs_cpt_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_0.50", &pi_minus_fromDs_cp_0_50, &b_pi_minus_fromDs_cp_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_0.50", &pi_minus_fromDs_cmult_0_50, &b_pi_minus_fromDs_cmult_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_0.50", &pi_minus_fromDs_deltaEta_0_50, &b_pi_minus_fromDs_deltaEta_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_0.50", &pi_minus_fromDs_deltaPhi_0_50, &b_pi_minus_fromDs_deltaPhi_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_0.50", &pi_minus_fromDs_pxasy_0_50, &b_pi_minus_fromDs_pxasy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_0.50", &pi_minus_fromDs_pyasy_0_50, &b_pi_minus_fromDs_pyasy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_0.50", &pi_minus_fromDs_pzasy_0_50, &b_pi_minus_fromDs_pzasy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_0.50", &pi_minus_fromDs_pasy_0_50, &b_pi_minus_fromDs_pasy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_0.50", &pi_minus_fromDs_ptasy_0_50, &b_pi_minus_fromDs_ptasy_0_50);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_0.60", &pi_minus_fromDs_cpx_0_60, &b_pi_minus_fromDs_cpx_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_0.60", &pi_minus_fromDs_cpy_0_60, &b_pi_minus_fromDs_cpy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_0.60", &pi_minus_fromDs_cpz_0_60, &b_pi_minus_fromDs_cpz_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_0.60", &pi_minus_fromDs_cpt_0_60, &b_pi_minus_fromDs_cpt_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_0.60", &pi_minus_fromDs_cp_0_60, &b_pi_minus_fromDs_cp_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_0.60", &pi_minus_fromDs_cmult_0_60, &b_pi_minus_fromDs_cmult_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_0.60", &pi_minus_fromDs_deltaEta_0_60, &b_pi_minus_fromDs_deltaEta_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_0.60", &pi_minus_fromDs_deltaPhi_0_60, &b_pi_minus_fromDs_deltaPhi_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_0.60", &pi_minus_fromDs_pxasy_0_60, &b_pi_minus_fromDs_pxasy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_0.60", &pi_minus_fromDs_pyasy_0_60, &b_pi_minus_fromDs_pyasy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_0.60", &pi_minus_fromDs_pzasy_0_60, &b_pi_minus_fromDs_pzasy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_0.60", &pi_minus_fromDs_pasy_0_60, &b_pi_minus_fromDs_pasy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_0.60", &pi_minus_fromDs_ptasy_0_60, &b_pi_minus_fromDs_ptasy_0_60);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_0.70", &pi_minus_fromDs_cpx_0_70, &b_pi_minus_fromDs_cpx_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_0.70", &pi_minus_fromDs_cpy_0_70, &b_pi_minus_fromDs_cpy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_0.70", &pi_minus_fromDs_cpz_0_70, &b_pi_minus_fromDs_cpz_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_0.70", &pi_minus_fromDs_cpt_0_70, &b_pi_minus_fromDs_cpt_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_0.70", &pi_minus_fromDs_cp_0_70, &b_pi_minus_fromDs_cp_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_0.70", &pi_minus_fromDs_cmult_0_70, &b_pi_minus_fromDs_cmult_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_0.70", &pi_minus_fromDs_deltaEta_0_70, &b_pi_minus_fromDs_deltaEta_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_0.70", &pi_minus_fromDs_deltaPhi_0_70, &b_pi_minus_fromDs_deltaPhi_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_0.70", &pi_minus_fromDs_pxasy_0_70, &b_pi_minus_fromDs_pxasy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_0.70", &pi_minus_fromDs_pyasy_0_70, &b_pi_minus_fromDs_pyasy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_0.70", &pi_minus_fromDs_pzasy_0_70, &b_pi_minus_fromDs_pzasy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_0.70", &pi_minus_fromDs_pasy_0_70, &b_pi_minus_fromDs_pasy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_0.70", &pi_minus_fromDs_ptasy_0_70, &b_pi_minus_fromDs_ptasy_0_70);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_0.80", &pi_minus_fromDs_cpx_0_80, &b_pi_minus_fromDs_cpx_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_0.80", &pi_minus_fromDs_cpy_0_80, &b_pi_minus_fromDs_cpy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_0.80", &pi_minus_fromDs_cpz_0_80, &b_pi_minus_fromDs_cpz_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_0.80", &pi_minus_fromDs_cpt_0_80, &b_pi_minus_fromDs_cpt_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_0.80", &pi_minus_fromDs_cp_0_80, &b_pi_minus_fromDs_cp_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_0.80", &pi_minus_fromDs_cmult_0_80, &b_pi_minus_fromDs_cmult_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_0.80", &pi_minus_fromDs_deltaEta_0_80, &b_pi_minus_fromDs_deltaEta_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_0.80", &pi_minus_fromDs_deltaPhi_0_80, &b_pi_minus_fromDs_deltaPhi_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_0.80", &pi_minus_fromDs_pxasy_0_80, &b_pi_minus_fromDs_pxasy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_0.80", &pi_minus_fromDs_pyasy_0_80, &b_pi_minus_fromDs_pyasy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_0.80", &pi_minus_fromDs_pzasy_0_80, &b_pi_minus_fromDs_pzasy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_0.80", &pi_minus_fromDs_pasy_0_80, &b_pi_minus_fromDs_pasy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_0.80", &pi_minus_fromDs_ptasy_0_80, &b_pi_minus_fromDs_ptasy_0_80);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_0.90", &pi_minus_fromDs_cpx_0_90, &b_pi_minus_fromDs_cpx_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_0.90", &pi_minus_fromDs_cpy_0_90, &b_pi_minus_fromDs_cpy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_0.90", &pi_minus_fromDs_cpz_0_90, &b_pi_minus_fromDs_cpz_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_0.90", &pi_minus_fromDs_cpt_0_90, &b_pi_minus_fromDs_cpt_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_0.90", &pi_minus_fromDs_cp_0_90, &b_pi_minus_fromDs_cp_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_0.90", &pi_minus_fromDs_cmult_0_90, &b_pi_minus_fromDs_cmult_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_0.90", &pi_minus_fromDs_deltaEta_0_90, &b_pi_minus_fromDs_deltaEta_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_0.90", &pi_minus_fromDs_deltaPhi_0_90, &b_pi_minus_fromDs_deltaPhi_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_0.90", &pi_minus_fromDs_pxasy_0_90, &b_pi_minus_fromDs_pxasy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_0.90", &pi_minus_fromDs_pyasy_0_90, &b_pi_minus_fromDs_pyasy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_0.90", &pi_minus_fromDs_pzasy_0_90, &b_pi_minus_fromDs_pzasy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_0.90", &pi_minus_fromDs_pasy_0_90, &b_pi_minus_fromDs_pasy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_0.90", &pi_minus_fromDs_ptasy_0_90, &b_pi_minus_fromDs_ptasy_0_90);
    fChain->SetBranchAddress("pi_minus_fromDs_cpx_1.00", &pi_minus_fromDs_cpx_1_00, &b_pi_minus_fromDs_cpx_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_cpy_1.00", &pi_minus_fromDs_cpy_1_00, &b_pi_minus_fromDs_cpy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_cpz_1.00", &pi_minus_fromDs_cpz_1_00, &b_pi_minus_fromDs_cpz_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_cpt_1.00", &pi_minus_fromDs_cpt_1_00, &b_pi_minus_fromDs_cpt_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_cp_1.00", &pi_minus_fromDs_cp_1_00, &b_pi_minus_fromDs_cp_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_cmult_1.00", &pi_minus_fromDs_cmult_1_00, &b_pi_minus_fromDs_cmult_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaEta_1.00", &pi_minus_fromDs_deltaEta_1_00, &b_pi_minus_fromDs_deltaEta_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_deltaPhi_1.00", &pi_minus_fromDs_deltaPhi_1_00, &b_pi_minus_fromDs_deltaPhi_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_pxasy_1.00", &pi_minus_fromDs_pxasy_1_00, &b_pi_minus_fromDs_pxasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_pyasy_1.00", &pi_minus_fromDs_pyasy_1_00, &b_pi_minus_fromDs_pyasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_pzasy_1.00", &pi_minus_fromDs_pzasy_1_00, &b_pi_minus_fromDs_pzasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_pasy_1.00", &pi_minus_fromDs_pasy_1_00, &b_pi_minus_fromDs_pasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_1.00", &pi_minus_fromDs_ptasy_1_00, &b_pi_minus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA1", &K_1_1270_plus_DOCA1, &b_K_1_1270_plus_DOCA1);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA2", &K_1_1270_plus_DOCA2, &b_K_1_1270_plus_DOCA2);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA3", &K_1_1270_plus_DOCA3, &b_K_1_1270_plus_DOCA3);
    fChain->SetBranchAddress("K_1_1270_plus_ETA", &K_1_1270_plus_ETA, &b_K_1_1270_plus_ETA);
    fChain->SetBranchAddress("K_1_1270_plus_CosTheta", &K_1_1270_plus_CosTheta, &b_K_1_1270_plus_CosTheta);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_X", &K_1_1270_plus_ENDVERTEX_X, &b_K_1_1270_plus_ENDVERTEX_X);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Y", &K_1_1270_plus_ENDVERTEX_Y, &b_K_1_1270_plus_ENDVERTEX_Y);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Z", &K_1_1270_plus_ENDVERTEX_Z, &b_K_1_1270_plus_ENDVERTEX_Z);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_XERR", &K_1_1270_plus_ENDVERTEX_XERR, &b_K_1_1270_plus_ENDVERTEX_XERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_YERR", &K_1_1270_plus_ENDVERTEX_YERR, &b_K_1_1270_plus_ENDVERTEX_YERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_ZERR", &K_1_1270_plus_ENDVERTEX_ZERR, &b_K_1_1270_plus_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_CHI2", &K_1_1270_plus_ENDVERTEX_CHI2, &b_K_1_1270_plus_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_NDOF", &K_1_1270_plus_ENDVERTEX_NDOF, &b_K_1_1270_plus_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_COV_", K_1_1270_plus_ENDVERTEX_COV_, &b_K_1_1270_plus_ENDVERTEX_COV_);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_X", &K_1_1270_plus_OWNPV_X, &b_K_1_1270_plus_OWNPV_X);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Y", &K_1_1270_plus_OWNPV_Y, &b_K_1_1270_plus_OWNPV_Y);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Z", &K_1_1270_plus_OWNPV_Z, &b_K_1_1270_plus_OWNPV_Z);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_XERR", &K_1_1270_plus_OWNPV_XERR, &b_K_1_1270_plus_OWNPV_XERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_YERR", &K_1_1270_plus_OWNPV_YERR, &b_K_1_1270_plus_OWNPV_YERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_ZERR", &K_1_1270_plus_OWNPV_ZERR, &b_K_1_1270_plus_OWNPV_ZERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_CHI2", &K_1_1270_plus_OWNPV_CHI2, &b_K_1_1270_plus_OWNPV_CHI2);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_NDOF", &K_1_1270_plus_OWNPV_NDOF, &b_K_1_1270_plus_OWNPV_NDOF);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_COV_", K_1_1270_plus_OWNPV_COV_, &b_K_1_1270_plus_OWNPV_COV_);
    fChain->SetBranchAddress("K_1_1270_plus_IP_OWNPV", &K_1_1270_plus_IP_OWNPV, &b_K_1_1270_plus_IP_OWNPV);
    fChain->SetBranchAddress("K_1_1270_plus_IPCHI2_OWNPV", &K_1_1270_plus_IPCHI2_OWNPV, &b_K_1_1270_plus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_1_1270_plus_FD_OWNPV", &K_1_1270_plus_FD_OWNPV, &b_K_1_1270_plus_FD_OWNPV);
    fChain->SetBranchAddress("K_1_1270_plus_FDCHI2_OWNPV", &K_1_1270_plus_FDCHI2_OWNPV, &b_K_1_1270_plus_FDCHI2_OWNPV);
    fChain->SetBranchAddress("K_1_1270_plus_DIRA_OWNPV", &K_1_1270_plus_DIRA_OWNPV, &b_K_1_1270_plus_DIRA_OWNPV);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_X", &K_1_1270_plus_ORIVX_X, &b_K_1_1270_plus_ORIVX_X);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_Y", &K_1_1270_plus_ORIVX_Y, &b_K_1_1270_plus_ORIVX_Y);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_Z", &K_1_1270_plus_ORIVX_Z, &b_K_1_1270_plus_ORIVX_Z);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_XERR", &K_1_1270_plus_ORIVX_XERR, &b_K_1_1270_plus_ORIVX_XERR);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_YERR", &K_1_1270_plus_ORIVX_YERR, &b_K_1_1270_plus_ORIVX_YERR);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_ZERR", &K_1_1270_plus_ORIVX_ZERR, &b_K_1_1270_plus_ORIVX_ZERR);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_CHI2", &K_1_1270_plus_ORIVX_CHI2, &b_K_1_1270_plus_ORIVX_CHI2);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_NDOF", &K_1_1270_plus_ORIVX_NDOF, &b_K_1_1270_plus_ORIVX_NDOF);
    fChain->SetBranchAddress("K_1_1270_plus_ORIVX_COV_", K_1_1270_plus_ORIVX_COV_, &b_K_1_1270_plus_ORIVX_COV_);
    fChain->SetBranchAddress("K_1_1270_plus_FD_ORIVX", &K_1_1270_plus_FD_ORIVX, &b_K_1_1270_plus_FD_ORIVX);
    fChain->SetBranchAddress("K_1_1270_plus_FDCHI2_ORIVX", &K_1_1270_plus_FDCHI2_ORIVX, &b_K_1_1270_plus_FDCHI2_ORIVX);
    fChain->SetBranchAddress("K_1_1270_plus_DIRA_ORIVX", &K_1_1270_plus_DIRA_ORIVX, &b_K_1_1270_plus_DIRA_ORIVX);
    fChain->SetBranchAddress("K_1_1270_plus_P", &K_1_1270_plus_P, &b_K_1_1270_plus_P);
    fChain->SetBranchAddress("K_1_1270_plus_PT", &K_1_1270_plus_PT, &b_K_1_1270_plus_PT);
    fChain->SetBranchAddress("K_1_1270_plus_PE", &K_1_1270_plus_PE, &b_K_1_1270_plus_PE);
    fChain->SetBranchAddress("K_1_1270_plus_PX", &K_1_1270_plus_PX, &b_K_1_1270_plus_PX);
    fChain->SetBranchAddress("K_1_1270_plus_PY", &K_1_1270_plus_PY, &b_K_1_1270_plus_PY);
    fChain->SetBranchAddress("K_1_1270_plus_PZ", &K_1_1270_plus_PZ, &b_K_1_1270_plus_PZ);
    fChain->SetBranchAddress("K_1_1270_plus_MM", &K_1_1270_plus_MM, &b_K_1_1270_plus_MM);
    fChain->SetBranchAddress("K_1_1270_plus_MMERR", &K_1_1270_plus_MMERR, &b_K_1_1270_plus_MMERR);
    fChain->SetBranchAddress("K_1_1270_plus_M", &K_1_1270_plus_M, &b_K_1_1270_plus_M);
    fChain->SetBranchAddress("K_1_1270_plus_ID", &K_1_1270_plus_ID, &b_K_1_1270_plus_ID);
    fChain->SetBranchAddress("K_1_1270_plus_TAU", &K_1_1270_plus_TAU, &b_K_1_1270_plus_TAU);
    fChain->SetBranchAddress("K_1_1270_plus_TAUERR", &K_1_1270_plus_TAUERR, &b_K_1_1270_plus_TAUERR);
    fChain->SetBranchAddress("K_1_1270_plus_TAUCHI2", &K_1_1270_plus_TAUCHI2, &b_K_1_1270_plus_TAUCHI2);
    fChain->SetBranchAddress("K_1_1270_plus_L0Global_Dec", &K_1_1270_plus_L0Global_Dec, &b_K_1_1270_plus_L0Global_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_L0Global_TIS", &K_1_1270_plus_L0Global_TIS, &b_K_1_1270_plus_L0Global_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_L0Global_TOS", &K_1_1270_plus_L0Global_TOS, &b_K_1_1270_plus_L0Global_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Global_Dec", &K_1_1270_plus_Hlt1Global_Dec, &b_K_1_1270_plus_Hlt1Global_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Global_TIS", &K_1_1270_plus_Hlt1Global_TIS, &b_K_1_1270_plus_Hlt1Global_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Global_TOS", &K_1_1270_plus_Hlt1Global_TOS, &b_K_1_1270_plus_Hlt1Global_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Phys_Dec", &K_1_1270_plus_Hlt1Phys_Dec, &b_K_1_1270_plus_Hlt1Phys_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Phys_TIS", &K_1_1270_plus_Hlt1Phys_TIS, &b_K_1_1270_plus_Hlt1Phys_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1Phys_TOS", &K_1_1270_plus_Hlt1Phys_TOS, &b_K_1_1270_plus_Hlt1Phys_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Global_Dec", &K_1_1270_plus_Hlt2Global_Dec, &b_K_1_1270_plus_Hlt2Global_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Global_TIS", &K_1_1270_plus_Hlt2Global_TIS, &b_K_1_1270_plus_Hlt2Global_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Global_TOS", &K_1_1270_plus_Hlt2Global_TOS, &b_K_1_1270_plus_Hlt2Global_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Phys_Dec", &K_1_1270_plus_Hlt2Phys_Dec, &b_K_1_1270_plus_Hlt2Phys_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Phys_TIS", &K_1_1270_plus_Hlt2Phys_TIS, &b_K_1_1270_plus_Hlt2Phys_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Phys_TOS", &K_1_1270_plus_Hlt2Phys_TOS, &b_K_1_1270_plus_Hlt2Phys_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_L0HadronDecision_Dec", &K_1_1270_plus_L0HadronDecision_Dec, &b_K_1_1270_plus_L0HadronDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_L0HadronDecision_TIS", &K_1_1270_plus_L0HadronDecision_TIS, &b_K_1_1270_plus_L0HadronDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_L0HadronDecision_TOS", &K_1_1270_plus_L0HadronDecision_TOS, &b_K_1_1270_plus_L0HadronDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_L0MuonDecision_Dec", &K_1_1270_plus_L0MuonDecision_Dec, &b_K_1_1270_plus_L0MuonDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_L0MuonDecision_TIS", &K_1_1270_plus_L0MuonDecision_TIS, &b_K_1_1270_plus_L0MuonDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_L0MuonDecision_TOS", &K_1_1270_plus_L0MuonDecision_TOS, &b_K_1_1270_plus_L0MuonDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_L0GlobalDecision_Dec", &K_1_1270_plus_L0GlobalDecision_Dec, &b_K_1_1270_plus_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_L0GlobalDecision_TIS", &K_1_1270_plus_L0GlobalDecision_TIS, &b_K_1_1270_plus_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_L0GlobalDecision_TOS", &K_1_1270_plus_L0GlobalDecision_TOS, &b_K_1_1270_plus_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackAllL0Decision_Dec", &K_1_1270_plus_Hlt1TrackAllL0Decision_Dec, &b_K_1_1270_plus_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackAllL0Decision_TIS", &K_1_1270_plus_Hlt1TrackAllL0Decision_TIS, &b_K_1_1270_plus_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackAllL0Decision_TOS", &K_1_1270_plus_Hlt1TrackAllL0Decision_TOS, &b_K_1_1270_plus_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVADecision_Dec", &K_1_1270_plus_Hlt1TrackMVADecision_Dec, &b_K_1_1270_plus_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVADecision_TIS", &K_1_1270_plus_Hlt1TrackMVADecision_TIS, &b_K_1_1270_plus_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVADecision_TOS", &K_1_1270_plus_Hlt1TrackMVADecision_TOS, &b_K_1_1270_plus_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVADecision_Dec", &K_1_1270_plus_Hlt1TwoTrackMVADecision_Dec, &b_K_1_1270_plus_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVADecision_TIS", &K_1_1270_plus_Hlt1TwoTrackMVADecision_TIS, &b_K_1_1270_plus_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVADecision_TOS", &K_1_1270_plus_Hlt1TwoTrackMVADecision_TOS, &b_K_1_1270_plus_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVALooseDecision_Dec", &K_1_1270_plus_Hlt1TrackMVALooseDecision_Dec, &b_K_1_1270_plus_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVALooseDecision_TIS", &K_1_1270_plus_Hlt1TrackMVALooseDecision_TIS, &b_K_1_1270_plus_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TrackMVALooseDecision_TOS", &K_1_1270_plus_Hlt1TrackMVALooseDecision_TOS, &b_K_1_1270_plus_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_Dec", &K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_Dec, &b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TIS", &K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TIS, &b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TOS", &K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TOS, &b_K_1_1270_plus_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2IncPhiDecision_Dec", &K_1_1270_plus_Hlt2IncPhiDecision_Dec, &b_K_1_1270_plus_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2IncPhiDecision_TIS", &K_1_1270_plus_Hlt2IncPhiDecision_TIS, &b_K_1_1270_plus_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2IncPhiDecision_TOS", &K_1_1270_plus_Hlt2IncPhiDecision_TOS, &b_K_1_1270_plus_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2PhiIncPhiDecision_Dec", &K_1_1270_plus_Hlt2PhiIncPhiDecision_Dec, &b_K_1_1270_plus_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2PhiIncPhiDecision_TIS", &K_1_1270_plus_Hlt2PhiIncPhiDecision_TIS, &b_K_1_1270_plus_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2PhiIncPhiDecision_TOS", &K_1_1270_plus_Hlt2PhiIncPhiDecision_TOS, &b_K_1_1270_plus_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_Dec", &K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_Dec, &b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TIS", &K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TIS, &b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TOS", &K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TOS, &b_K_1_1270_plus_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_Dec", &K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_Dec, &b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TIS", &K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TIS, &b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TOS", &K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TOS, &b_K_1_1270_plus_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_Dec", &K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_Dec, &b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TIS", &K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TIS, &b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TOS", &K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TOS, &b_K_1_1270_plus_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyDecision_Dec", &K_1_1270_plus_Hlt2Topo2BodyDecision_Dec, &b_K_1_1270_plus_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyDecision_TIS", &K_1_1270_plus_Hlt2Topo2BodyDecision_TIS, &b_K_1_1270_plus_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo2BodyDecision_TOS", &K_1_1270_plus_Hlt2Topo2BodyDecision_TOS, &b_K_1_1270_plus_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyDecision_Dec", &K_1_1270_plus_Hlt2Topo3BodyDecision_Dec, &b_K_1_1270_plus_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyDecision_TIS", &K_1_1270_plus_Hlt2Topo3BodyDecision_TIS, &b_K_1_1270_plus_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo3BodyDecision_TOS", &K_1_1270_plus_Hlt2Topo3BodyDecision_TOS, &b_K_1_1270_plus_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyDecision_Dec", &K_1_1270_plus_Hlt2Topo4BodyDecision_Dec, &b_K_1_1270_plus_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyDecision_TIS", &K_1_1270_plus_Hlt2Topo4BodyDecision_TIS, &b_K_1_1270_plus_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("K_1_1270_plus_Hlt2Topo4BodyDecision_TOS", &K_1_1270_plus_Hlt2Topo4BodyDecision_TOS, &b_K_1_1270_plus_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_0.50", &K_1_1270_plus_cpx_0_50, &b_K_1_1270_plus_cpx_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_0.50", &K_1_1270_plus_cpy_0_50, &b_K_1_1270_plus_cpy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_0.50", &K_1_1270_plus_cpz_0_50, &b_K_1_1270_plus_cpz_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_0.50", &K_1_1270_plus_cpt_0_50, &b_K_1_1270_plus_cpt_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cp_0.50", &K_1_1270_plus_cp_0_50, &b_K_1_1270_plus_cp_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_0.50", &K_1_1270_plus_cmult_0_50, &b_K_1_1270_plus_cmult_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_0.50", &K_1_1270_plus_deltaEta_0_50, &b_K_1_1270_plus_deltaEta_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_0.50", &K_1_1270_plus_deltaPhi_0_50, &b_K_1_1270_plus_deltaPhi_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_0.50", &K_1_1270_plus_pxasy_0_50, &b_K_1_1270_plus_pxasy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_0.50", &K_1_1270_plus_pyasy_0_50, &b_K_1_1270_plus_pyasy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_0.50", &K_1_1270_plus_pzasy_0_50, &b_K_1_1270_plus_pzasy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_0.50", &K_1_1270_plus_pasy_0_50, &b_K_1_1270_plus_pasy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_0.50", &K_1_1270_plus_ptasy_0_50, &b_K_1_1270_plus_ptasy_0_50);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_0.60", &K_1_1270_plus_cpx_0_60, &b_K_1_1270_plus_cpx_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_0.60", &K_1_1270_plus_cpy_0_60, &b_K_1_1270_plus_cpy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_0.60", &K_1_1270_plus_cpz_0_60, &b_K_1_1270_plus_cpz_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_0.60", &K_1_1270_plus_cpt_0_60, &b_K_1_1270_plus_cpt_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cp_0.60", &K_1_1270_plus_cp_0_60, &b_K_1_1270_plus_cp_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_0.60", &K_1_1270_plus_cmult_0_60, &b_K_1_1270_plus_cmult_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_0.60", &K_1_1270_plus_deltaEta_0_60, &b_K_1_1270_plus_deltaEta_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_0.60", &K_1_1270_plus_deltaPhi_0_60, &b_K_1_1270_plus_deltaPhi_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_0.60", &K_1_1270_plus_pxasy_0_60, &b_K_1_1270_plus_pxasy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_0.60", &K_1_1270_plus_pyasy_0_60, &b_K_1_1270_plus_pyasy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_0.60", &K_1_1270_plus_pzasy_0_60, &b_K_1_1270_plus_pzasy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_0.60", &K_1_1270_plus_pasy_0_60, &b_K_1_1270_plus_pasy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_0.60", &K_1_1270_plus_ptasy_0_60, &b_K_1_1270_plus_ptasy_0_60);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_0.70", &K_1_1270_plus_cpx_0_70, &b_K_1_1270_plus_cpx_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_0.70", &K_1_1270_plus_cpy_0_70, &b_K_1_1270_plus_cpy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_0.70", &K_1_1270_plus_cpz_0_70, &b_K_1_1270_plus_cpz_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_0.70", &K_1_1270_plus_cpt_0_70, &b_K_1_1270_plus_cpt_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cp_0.70", &K_1_1270_plus_cp_0_70, &b_K_1_1270_plus_cp_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_0.70", &K_1_1270_plus_cmult_0_70, &b_K_1_1270_plus_cmult_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_0.70", &K_1_1270_plus_deltaEta_0_70, &b_K_1_1270_plus_deltaEta_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_0.70", &K_1_1270_plus_deltaPhi_0_70, &b_K_1_1270_plus_deltaPhi_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_0.70", &K_1_1270_plus_pxasy_0_70, &b_K_1_1270_plus_pxasy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_0.70", &K_1_1270_plus_pyasy_0_70, &b_K_1_1270_plus_pyasy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_0.70", &K_1_1270_plus_pzasy_0_70, &b_K_1_1270_plus_pzasy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_0.70", &K_1_1270_plus_pasy_0_70, &b_K_1_1270_plus_pasy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_0.70", &K_1_1270_plus_ptasy_0_70, &b_K_1_1270_plus_ptasy_0_70);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_0.80", &K_1_1270_plus_cpx_0_80, &b_K_1_1270_plus_cpx_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_0.80", &K_1_1270_plus_cpy_0_80, &b_K_1_1270_plus_cpy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_0.80", &K_1_1270_plus_cpz_0_80, &b_K_1_1270_plus_cpz_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_0.80", &K_1_1270_plus_cpt_0_80, &b_K_1_1270_plus_cpt_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cp_0.80", &K_1_1270_plus_cp_0_80, &b_K_1_1270_plus_cp_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_0.80", &K_1_1270_plus_cmult_0_80, &b_K_1_1270_plus_cmult_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_0.80", &K_1_1270_plus_deltaEta_0_80, &b_K_1_1270_plus_deltaEta_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_0.80", &K_1_1270_plus_deltaPhi_0_80, &b_K_1_1270_plus_deltaPhi_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_0.80", &K_1_1270_plus_pxasy_0_80, &b_K_1_1270_plus_pxasy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_0.80", &K_1_1270_plus_pyasy_0_80, &b_K_1_1270_plus_pyasy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_0.80", &K_1_1270_plus_pzasy_0_80, &b_K_1_1270_plus_pzasy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_0.80", &K_1_1270_plus_pasy_0_80, &b_K_1_1270_plus_pasy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_0.80", &K_1_1270_plus_ptasy_0_80, &b_K_1_1270_plus_ptasy_0_80);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_0.90", &K_1_1270_plus_cpx_0_90, &b_K_1_1270_plus_cpx_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_0.90", &K_1_1270_plus_cpy_0_90, &b_K_1_1270_plus_cpy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_0.90", &K_1_1270_plus_cpz_0_90, &b_K_1_1270_plus_cpz_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_0.90", &K_1_1270_plus_cpt_0_90, &b_K_1_1270_plus_cpt_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cp_0.90", &K_1_1270_plus_cp_0_90, &b_K_1_1270_plus_cp_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_0.90", &K_1_1270_plus_cmult_0_90, &b_K_1_1270_plus_cmult_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_0.90", &K_1_1270_plus_deltaEta_0_90, &b_K_1_1270_plus_deltaEta_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_0.90", &K_1_1270_plus_deltaPhi_0_90, &b_K_1_1270_plus_deltaPhi_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_0.90", &K_1_1270_plus_pxasy_0_90, &b_K_1_1270_plus_pxasy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_0.90", &K_1_1270_plus_pyasy_0_90, &b_K_1_1270_plus_pyasy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_0.90", &K_1_1270_plus_pzasy_0_90, &b_K_1_1270_plus_pzasy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_0.90", &K_1_1270_plus_pasy_0_90, &b_K_1_1270_plus_pasy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_0.90", &K_1_1270_plus_ptasy_0_90, &b_K_1_1270_plus_ptasy_0_90);
    fChain->SetBranchAddress("K_1_1270_plus_cpx_1.00", &K_1_1270_plus_cpx_1_00, &b_K_1_1270_plus_cpx_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_cpy_1.00", &K_1_1270_plus_cpy_1_00, &b_K_1_1270_plus_cpy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_cpz_1.00", &K_1_1270_plus_cpz_1_00, &b_K_1_1270_plus_cpz_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_cpt_1.00", &K_1_1270_plus_cpt_1_00, &b_K_1_1270_plus_cpt_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_cp_1.00", &K_1_1270_plus_cp_1_00, &b_K_1_1270_plus_cp_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_cmult_1.00", &K_1_1270_plus_cmult_1_00, &b_K_1_1270_plus_cmult_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_deltaEta_1.00", &K_1_1270_plus_deltaEta_1_00, &b_K_1_1270_plus_deltaEta_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_deltaPhi_1.00", &K_1_1270_plus_deltaPhi_1_00, &b_K_1_1270_plus_deltaPhi_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_pxasy_1.00", &K_1_1270_plus_pxasy_1_00, &b_K_1_1270_plus_pxasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_pyasy_1.00", &K_1_1270_plus_pyasy_1_00, &b_K_1_1270_plus_pyasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_pzasy_1.00", &K_1_1270_plus_pzasy_1_00, &b_K_1_1270_plus_pzasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_pasy_1.00", &K_1_1270_plus_pasy_1_00, &b_K_1_1270_plus_pasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_1.00", &K_1_1270_plus_ptasy_1_00, &b_K_1_1270_plus_ptasy_1_00);
    fChain->SetBranchAddress("K_plus_DOCA1", &K_plus_DOCA1, &b_K_plus_DOCA1);
    fChain->SetBranchAddress("K_plus_DOCA2", &K_plus_DOCA2, &b_K_plus_DOCA2);
    fChain->SetBranchAddress("K_plus_DOCA3", &K_plus_DOCA3, &b_K_plus_DOCA3);
    fChain->SetBranchAddress("K_plus_ETA", &K_plus_ETA, &b_K_plus_ETA);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNe", &K_plus_MC12TuneV2_ProbNNe, &b_K_plus_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNmu", &K_plus_MC12TuneV2_ProbNNmu, &b_K_plus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNpi", &K_plus_MC12TuneV2_ProbNNpi, &b_K_plus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNk", &K_plus_MC12TuneV2_ProbNNk, &b_K_plus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNp", &K_plus_MC12TuneV2_ProbNNp, &b_K_plus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNghost", &K_plus_MC12TuneV2_ProbNNghost, &b_K_plus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNe", &K_plus_MC12TuneV3_ProbNNe, &b_K_plus_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNmu", &K_plus_MC12TuneV3_ProbNNmu, &b_K_plus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNpi", &K_plus_MC12TuneV3_ProbNNpi, &b_K_plus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNk", &K_plus_MC12TuneV3_ProbNNk, &b_K_plus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNp", &K_plus_MC12TuneV3_ProbNNp, &b_K_plus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNghost", &K_plus_MC12TuneV3_ProbNNghost, &b_K_plus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNe", &K_plus_MC12TuneV4_ProbNNe, &b_K_plus_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNmu", &K_plus_MC12TuneV4_ProbNNmu, &b_K_plus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNpi", &K_plus_MC12TuneV4_ProbNNpi, &b_K_plus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNk", &K_plus_MC12TuneV4_ProbNNk, &b_K_plus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNp", &K_plus_MC12TuneV4_ProbNNp, &b_K_plus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNghost", &K_plus_MC12TuneV4_ProbNNghost, &b_K_plus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNe", &K_plus_MC15TuneV1_ProbNNe, &b_K_plus_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNmu", &K_plus_MC15TuneV1_ProbNNmu, &b_K_plus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNpi", &K_plus_MC15TuneV1_ProbNNpi, &b_K_plus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNk", &K_plus_MC15TuneV1_ProbNNk, &b_K_plus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNp", &K_plus_MC15TuneV1_ProbNNp, &b_K_plus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNghost", &K_plus_MC15TuneV1_ProbNNghost, &b_K_plus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_plus_CosTheta", &K_plus_CosTheta, &b_K_plus_CosTheta);
    fChain->SetBranchAddress("K_plus_OWNPV_X", &K_plus_OWNPV_X, &b_K_plus_OWNPV_X);
    fChain->SetBranchAddress("K_plus_OWNPV_Y", &K_plus_OWNPV_Y, &b_K_plus_OWNPV_Y);
    fChain->SetBranchAddress("K_plus_OWNPV_Z", &K_plus_OWNPV_Z, &b_K_plus_OWNPV_Z);
    fChain->SetBranchAddress("K_plus_OWNPV_XERR", &K_plus_OWNPV_XERR, &b_K_plus_OWNPV_XERR);
    fChain->SetBranchAddress("K_plus_OWNPV_YERR", &K_plus_OWNPV_YERR, &b_K_plus_OWNPV_YERR);
    fChain->SetBranchAddress("K_plus_OWNPV_ZERR", &K_plus_OWNPV_ZERR, &b_K_plus_OWNPV_ZERR);
    fChain->SetBranchAddress("K_plus_OWNPV_CHI2", &K_plus_OWNPV_CHI2, &b_K_plus_OWNPV_CHI2);
    fChain->SetBranchAddress("K_plus_OWNPV_NDOF", &K_plus_OWNPV_NDOF, &b_K_plus_OWNPV_NDOF);
    fChain->SetBranchAddress("K_plus_OWNPV_COV_", K_plus_OWNPV_COV_, &b_K_plus_OWNPV_COV_);
    fChain->SetBranchAddress("K_plus_IP_OWNPV", &K_plus_IP_OWNPV, &b_K_plus_IP_OWNPV);
    fChain->SetBranchAddress("K_plus_IPCHI2_OWNPV", &K_plus_IPCHI2_OWNPV, &b_K_plus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_plus_ORIVX_X", &K_plus_ORIVX_X, &b_K_plus_ORIVX_X);
    fChain->SetBranchAddress("K_plus_ORIVX_Y", &K_plus_ORIVX_Y, &b_K_plus_ORIVX_Y);
    fChain->SetBranchAddress("K_plus_ORIVX_Z", &K_plus_ORIVX_Z, &b_K_plus_ORIVX_Z);
    fChain->SetBranchAddress("K_plus_ORIVX_XERR", &K_plus_ORIVX_XERR, &b_K_plus_ORIVX_XERR);
    fChain->SetBranchAddress("K_plus_ORIVX_YERR", &K_plus_ORIVX_YERR, &b_K_plus_ORIVX_YERR);
    fChain->SetBranchAddress("K_plus_ORIVX_ZERR", &K_plus_ORIVX_ZERR, &b_K_plus_ORIVX_ZERR);
    fChain->SetBranchAddress("K_plus_ORIVX_CHI2", &K_plus_ORIVX_CHI2, &b_K_plus_ORIVX_CHI2);
    fChain->SetBranchAddress("K_plus_ORIVX_NDOF", &K_plus_ORIVX_NDOF, &b_K_plus_ORIVX_NDOF);
    fChain->SetBranchAddress("K_plus_ORIVX_COV_", K_plus_ORIVX_COV_, &b_K_plus_ORIVX_COV_);
    fChain->SetBranchAddress("K_plus_P", &K_plus_P, &b_K_plus_P);
    fChain->SetBranchAddress("K_plus_PT", &K_plus_PT, &b_K_plus_PT);
    fChain->SetBranchAddress("K_plus_PE", &K_plus_PE, &b_K_plus_PE);
    fChain->SetBranchAddress("K_plus_PX", &K_plus_PX, &b_K_plus_PX);
    fChain->SetBranchAddress("K_plus_PY", &K_plus_PY, &b_K_plus_PY);
    fChain->SetBranchAddress("K_plus_PZ", &K_plus_PZ, &b_K_plus_PZ);
    fChain->SetBranchAddress("K_plus_M", &K_plus_M, &b_K_plus_M);
    fChain->SetBranchAddress("K_plus_ID", &K_plus_ID, &b_K_plus_ID);
    fChain->SetBranchAddress("K_plus_PIDe", &K_plus_PIDe, &b_K_plus_PIDe);
    fChain->SetBranchAddress("K_plus_PIDmu", &K_plus_PIDmu, &b_K_plus_PIDmu);
    fChain->SetBranchAddress("K_plus_PIDK", &K_plus_PIDK, &b_K_plus_PIDK);
    fChain->SetBranchAddress("K_plus_PIDp", &K_plus_PIDp, &b_K_plus_PIDp);
    fChain->SetBranchAddress("K_plus_ProbNNe", &K_plus_ProbNNe, &b_K_plus_ProbNNe);
    fChain->SetBranchAddress("K_plus_ProbNNk", &K_plus_ProbNNk, &b_K_plus_ProbNNk);
    fChain->SetBranchAddress("K_plus_ProbNNp", &K_plus_ProbNNp, &b_K_plus_ProbNNp);
    fChain->SetBranchAddress("K_plus_ProbNNpi", &K_plus_ProbNNpi, &b_K_plus_ProbNNpi);
    fChain->SetBranchAddress("K_plus_ProbNNmu", &K_plus_ProbNNmu, &b_K_plus_ProbNNmu);
    fChain->SetBranchAddress("K_plus_ProbNNghost", &K_plus_ProbNNghost, &b_K_plus_ProbNNghost);
    fChain->SetBranchAddress("K_plus_hasMuon", &K_plus_hasMuon, &b_K_plus_hasMuon);
    fChain->SetBranchAddress("K_plus_isMuon", &K_plus_isMuon, &b_K_plus_isMuon);
    fChain->SetBranchAddress("K_plus_hasRich", &K_plus_hasRich, &b_K_plus_hasRich);
    fChain->SetBranchAddress("K_plus_UsedRichAerogel", &K_plus_UsedRichAerogel, &b_K_plus_UsedRichAerogel);
    fChain->SetBranchAddress("K_plus_UsedRich1Gas", &K_plus_UsedRich1Gas, &b_K_plus_UsedRich1Gas);
    fChain->SetBranchAddress("K_plus_UsedRich2Gas", &K_plus_UsedRich2Gas, &b_K_plus_UsedRich2Gas);
    fChain->SetBranchAddress("K_plus_RichAboveElThres", &K_plus_RichAboveElThres, &b_K_plus_RichAboveElThres);
    fChain->SetBranchAddress("K_plus_RichAboveMuThres", &K_plus_RichAboveMuThres, &b_K_plus_RichAboveMuThres);
    fChain->SetBranchAddress("K_plus_RichAbovePiThres", &K_plus_RichAbovePiThres, &b_K_plus_RichAbovePiThres);
    fChain->SetBranchAddress("K_plus_RichAboveKaThres", &K_plus_RichAboveKaThres, &b_K_plus_RichAboveKaThres);
    fChain->SetBranchAddress("K_plus_RichAbovePrThres", &K_plus_RichAbovePrThres, &b_K_plus_RichAbovePrThres);
    fChain->SetBranchAddress("K_plus_hasCalo", &K_plus_hasCalo, &b_K_plus_hasCalo);
    fChain->SetBranchAddress("K_plus_L0Global_Dec", &K_plus_L0Global_Dec, &b_K_plus_L0Global_Dec);
    fChain->SetBranchAddress("K_plus_L0Global_TIS", &K_plus_L0Global_TIS, &b_K_plus_L0Global_TIS);
    fChain->SetBranchAddress("K_plus_L0Global_TOS", &K_plus_L0Global_TOS, &b_K_plus_L0Global_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1Global_Dec", &K_plus_Hlt1Global_Dec, &b_K_plus_Hlt1Global_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1Global_TIS", &K_plus_Hlt1Global_TIS, &b_K_plus_Hlt1Global_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1Global_TOS", &K_plus_Hlt1Global_TOS, &b_K_plus_Hlt1Global_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1Phys_Dec", &K_plus_Hlt1Phys_Dec, &b_K_plus_Hlt1Phys_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1Phys_TIS", &K_plus_Hlt1Phys_TIS, &b_K_plus_Hlt1Phys_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1Phys_TOS", &K_plus_Hlt1Phys_TOS, &b_K_plus_Hlt1Phys_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Global_Dec", &K_plus_Hlt2Global_Dec, &b_K_plus_Hlt2Global_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Global_TIS", &K_plus_Hlt2Global_TIS, &b_K_plus_Hlt2Global_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Global_TOS", &K_plus_Hlt2Global_TOS, &b_K_plus_Hlt2Global_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Phys_Dec", &K_plus_Hlt2Phys_Dec, &b_K_plus_Hlt2Phys_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Phys_TIS", &K_plus_Hlt2Phys_TIS, &b_K_plus_Hlt2Phys_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Phys_TOS", &K_plus_Hlt2Phys_TOS, &b_K_plus_Hlt2Phys_TOS);
    fChain->SetBranchAddress("K_plus_L0HadronDecision_Dec", &K_plus_L0HadronDecision_Dec, &b_K_plus_L0HadronDecision_Dec);
    fChain->SetBranchAddress("K_plus_L0HadronDecision_TIS", &K_plus_L0HadronDecision_TIS, &b_K_plus_L0HadronDecision_TIS);
    fChain->SetBranchAddress("K_plus_L0HadronDecision_TOS", &K_plus_L0HadronDecision_TOS, &b_K_plus_L0HadronDecision_TOS);
    fChain->SetBranchAddress("K_plus_L0MuonDecision_Dec", &K_plus_L0MuonDecision_Dec, &b_K_plus_L0MuonDecision_Dec);
    fChain->SetBranchAddress("K_plus_L0MuonDecision_TIS", &K_plus_L0MuonDecision_TIS, &b_K_plus_L0MuonDecision_TIS);
    fChain->SetBranchAddress("K_plus_L0MuonDecision_TOS", &K_plus_L0MuonDecision_TOS, &b_K_plus_L0MuonDecision_TOS);
    fChain->SetBranchAddress("K_plus_L0GlobalDecision_Dec", &K_plus_L0GlobalDecision_Dec, &b_K_plus_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("K_plus_L0GlobalDecision_TIS", &K_plus_L0GlobalDecision_TIS, &b_K_plus_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("K_plus_L0GlobalDecision_TOS", &K_plus_L0GlobalDecision_TOS, &b_K_plus_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackAllL0Decision_Dec", &K_plus_Hlt1TrackAllL0Decision_Dec, &b_K_plus_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1TrackAllL0Decision_TIS", &K_plus_Hlt1TrackAllL0Decision_TIS, &b_K_plus_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackAllL0Decision_TOS", &K_plus_Hlt1TrackAllL0Decision_TOS, &b_K_plus_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVADecision_Dec", &K_plus_Hlt1TrackMVADecision_Dec, &b_K_plus_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVADecision_TIS", &K_plus_Hlt1TrackMVADecision_TIS, &b_K_plus_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVADecision_TOS", &K_plus_Hlt1TrackMVADecision_TOS, &b_K_plus_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVADecision_Dec", &K_plus_Hlt1TwoTrackMVADecision_Dec, &b_K_plus_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVADecision_TIS", &K_plus_Hlt1TwoTrackMVADecision_TIS, &b_K_plus_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVADecision_TOS", &K_plus_Hlt1TwoTrackMVADecision_TOS, &b_K_plus_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVALooseDecision_Dec", &K_plus_Hlt1TrackMVALooseDecision_Dec, &b_K_plus_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVALooseDecision_TIS", &K_plus_Hlt1TrackMVALooseDecision_TIS, &b_K_plus_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1TrackMVALooseDecision_TOS", &K_plus_Hlt1TrackMVALooseDecision_TOS, &b_K_plus_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVALooseDecision_Dec", &K_plus_Hlt1TwoTrackMVALooseDecision_Dec, &b_K_plus_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVALooseDecision_TIS", &K_plus_Hlt1TwoTrackMVALooseDecision_TIS, &b_K_plus_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt1TwoTrackMVALooseDecision_TOS", &K_plus_Hlt1TwoTrackMVALooseDecision_TOS, &b_K_plus_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2IncPhiDecision_Dec", &K_plus_Hlt2IncPhiDecision_Dec, &b_K_plus_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2IncPhiDecision_TIS", &K_plus_Hlt2IncPhiDecision_TIS, &b_K_plus_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2IncPhiDecision_TOS", &K_plus_Hlt2IncPhiDecision_TOS, &b_K_plus_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2PhiIncPhiDecision_Dec", &K_plus_Hlt2PhiIncPhiDecision_Dec, &b_K_plus_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2PhiIncPhiDecision_TIS", &K_plus_Hlt2PhiIncPhiDecision_TIS, &b_K_plus_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2PhiIncPhiDecision_TOS", &K_plus_Hlt2PhiIncPhiDecision_TOS, &b_K_plus_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyBBDTDecision_Dec", &K_plus_Hlt2Topo2BodyBBDTDecision_Dec, &b_K_plus_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyBBDTDecision_TIS", &K_plus_Hlt2Topo2BodyBBDTDecision_TIS, &b_K_plus_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyBBDTDecision_TOS", &K_plus_Hlt2Topo2BodyBBDTDecision_TOS, &b_K_plus_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyBBDTDecision_Dec", &K_plus_Hlt2Topo3BodyBBDTDecision_Dec, &b_K_plus_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyBBDTDecision_TIS", &K_plus_Hlt2Topo3BodyBBDTDecision_TIS, &b_K_plus_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyBBDTDecision_TOS", &K_plus_Hlt2Topo3BodyBBDTDecision_TOS, &b_K_plus_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyBBDTDecision_Dec", &K_plus_Hlt2Topo4BodyBBDTDecision_Dec, &b_K_plus_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyBBDTDecision_TIS", &K_plus_Hlt2Topo4BodyBBDTDecision_TIS, &b_K_plus_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyBBDTDecision_TOS", &K_plus_Hlt2Topo4BodyBBDTDecision_TOS, &b_K_plus_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyDecision_Dec", &K_plus_Hlt2Topo2BodyDecision_Dec, &b_K_plus_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyDecision_TIS", &K_plus_Hlt2Topo2BodyDecision_TIS, &b_K_plus_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo2BodyDecision_TOS", &K_plus_Hlt2Topo2BodyDecision_TOS, &b_K_plus_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyDecision_Dec", &K_plus_Hlt2Topo3BodyDecision_Dec, &b_K_plus_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyDecision_TIS", &K_plus_Hlt2Topo3BodyDecision_TIS, &b_K_plus_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo3BodyDecision_TOS", &K_plus_Hlt2Topo3BodyDecision_TOS, &b_K_plus_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyDecision_Dec", &K_plus_Hlt2Topo4BodyDecision_Dec, &b_K_plus_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyDecision_TIS", &K_plus_Hlt2Topo4BodyDecision_TIS, &b_K_plus_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("K_plus_Hlt2Topo4BodyDecision_TOS", &K_plus_Hlt2Topo4BodyDecision_TOS, &b_K_plus_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("K_plus_TRACK_Type", &K_plus_TRACK_Type, &b_K_plus_TRACK_Type);
    fChain->SetBranchAddress("K_plus_TRACK_Key", &K_plus_TRACK_Key, &b_K_plus_TRACK_Key);
    fChain->SetBranchAddress("K_plus_TRACK_CHI2NDOF", &K_plus_TRACK_CHI2NDOF, &b_K_plus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_plus_TRACK_PCHI2", &K_plus_TRACK_PCHI2, &b_K_plus_TRACK_PCHI2);
    fChain->SetBranchAddress("K_plus_TRACK_MatchCHI2", &K_plus_TRACK_MatchCHI2, &b_K_plus_TRACK_MatchCHI2);
    fChain->SetBranchAddress("K_plus_TRACK_GhostProb", &K_plus_TRACK_GhostProb, &b_K_plus_TRACK_GhostProb);
    fChain->SetBranchAddress("K_plus_TRACK_CloneDist", &K_plus_TRACK_CloneDist, &b_K_plus_TRACK_CloneDist);
    fChain->SetBranchAddress("K_plus_TRACK_Likelihood", &K_plus_TRACK_Likelihood, &b_K_plus_TRACK_Likelihood);
    fChain->SetBranchAddress("K_plus_cpx_0.50", &K_plus_cpx_0_50, &b_K_plus_cpx_0_50);
    fChain->SetBranchAddress("K_plus_cpy_0.50", &K_plus_cpy_0_50, &b_K_plus_cpy_0_50);
    fChain->SetBranchAddress("K_plus_cpz_0.50", &K_plus_cpz_0_50, &b_K_plus_cpz_0_50);
    fChain->SetBranchAddress("K_plus_cpt_0.50", &K_plus_cpt_0_50, &b_K_plus_cpt_0_50);
    fChain->SetBranchAddress("K_plus_cp_0.50", &K_plus_cp_0_50, &b_K_plus_cp_0_50);
    fChain->SetBranchAddress("K_plus_cmult_0.50", &K_plus_cmult_0_50, &b_K_plus_cmult_0_50);
    fChain->SetBranchAddress("K_plus_deltaEta_0.50", &K_plus_deltaEta_0_50, &b_K_plus_deltaEta_0_50);
    fChain->SetBranchAddress("K_plus_deltaPhi_0.50", &K_plus_deltaPhi_0_50, &b_K_plus_deltaPhi_0_50);
    fChain->SetBranchAddress("K_plus_pxasy_0.50", &K_plus_pxasy_0_50, &b_K_plus_pxasy_0_50);
    fChain->SetBranchAddress("K_plus_pyasy_0.50", &K_plus_pyasy_0_50, &b_K_plus_pyasy_0_50);
    fChain->SetBranchAddress("K_plus_pzasy_0.50", &K_plus_pzasy_0_50, &b_K_plus_pzasy_0_50);
    fChain->SetBranchAddress("K_plus_pasy_0.50", &K_plus_pasy_0_50, &b_K_plus_pasy_0_50);
    fChain->SetBranchAddress("K_plus_ptasy_0.50", &K_plus_ptasy_0_50, &b_K_plus_ptasy_0_50);
    fChain->SetBranchAddress("K_plus_cpx_0.60", &K_plus_cpx_0_60, &b_K_plus_cpx_0_60);
    fChain->SetBranchAddress("K_plus_cpy_0.60", &K_plus_cpy_0_60, &b_K_plus_cpy_0_60);
    fChain->SetBranchAddress("K_plus_cpz_0.60", &K_plus_cpz_0_60, &b_K_plus_cpz_0_60);
    fChain->SetBranchAddress("K_plus_cpt_0.60", &K_plus_cpt_0_60, &b_K_plus_cpt_0_60);
    fChain->SetBranchAddress("K_plus_cp_0.60", &K_plus_cp_0_60, &b_K_plus_cp_0_60);
    fChain->SetBranchAddress("K_plus_cmult_0.60", &K_plus_cmult_0_60, &b_K_plus_cmult_0_60);
    fChain->SetBranchAddress("K_plus_deltaEta_0.60", &K_plus_deltaEta_0_60, &b_K_plus_deltaEta_0_60);
    fChain->SetBranchAddress("K_plus_deltaPhi_0.60", &K_plus_deltaPhi_0_60, &b_K_plus_deltaPhi_0_60);
    fChain->SetBranchAddress("K_plus_pxasy_0.60", &K_plus_pxasy_0_60, &b_K_plus_pxasy_0_60);
    fChain->SetBranchAddress("K_plus_pyasy_0.60", &K_plus_pyasy_0_60, &b_K_plus_pyasy_0_60);
    fChain->SetBranchAddress("K_plus_pzasy_0.60", &K_plus_pzasy_0_60, &b_K_plus_pzasy_0_60);
    fChain->SetBranchAddress("K_plus_pasy_0.60", &K_plus_pasy_0_60, &b_K_plus_pasy_0_60);
    fChain->SetBranchAddress("K_plus_ptasy_0.60", &K_plus_ptasy_0_60, &b_K_plus_ptasy_0_60);
    fChain->SetBranchAddress("K_plus_cpx_0.70", &K_plus_cpx_0_70, &b_K_plus_cpx_0_70);
    fChain->SetBranchAddress("K_plus_cpy_0.70", &K_plus_cpy_0_70, &b_K_plus_cpy_0_70);
    fChain->SetBranchAddress("K_plus_cpz_0.70", &K_plus_cpz_0_70, &b_K_plus_cpz_0_70);
    fChain->SetBranchAddress("K_plus_cpt_0.70", &K_plus_cpt_0_70, &b_K_plus_cpt_0_70);
    fChain->SetBranchAddress("K_plus_cp_0.70", &K_plus_cp_0_70, &b_K_plus_cp_0_70);
    fChain->SetBranchAddress("K_plus_cmult_0.70", &K_plus_cmult_0_70, &b_K_plus_cmult_0_70);
    fChain->SetBranchAddress("K_plus_deltaEta_0.70", &K_plus_deltaEta_0_70, &b_K_plus_deltaEta_0_70);
    fChain->SetBranchAddress("K_plus_deltaPhi_0.70", &K_plus_deltaPhi_0_70, &b_K_plus_deltaPhi_0_70);
    fChain->SetBranchAddress("K_plus_pxasy_0.70", &K_plus_pxasy_0_70, &b_K_plus_pxasy_0_70);
    fChain->SetBranchAddress("K_plus_pyasy_0.70", &K_plus_pyasy_0_70, &b_K_plus_pyasy_0_70);
    fChain->SetBranchAddress("K_plus_pzasy_0.70", &K_plus_pzasy_0_70, &b_K_plus_pzasy_0_70);
    fChain->SetBranchAddress("K_plus_pasy_0.70", &K_plus_pasy_0_70, &b_K_plus_pasy_0_70);
    fChain->SetBranchAddress("K_plus_ptasy_0.70", &K_plus_ptasy_0_70, &b_K_plus_ptasy_0_70);
    fChain->SetBranchAddress("K_plus_cpx_0.80", &K_plus_cpx_0_80, &b_K_plus_cpx_0_80);
    fChain->SetBranchAddress("K_plus_cpy_0.80", &K_plus_cpy_0_80, &b_K_plus_cpy_0_80);
    fChain->SetBranchAddress("K_plus_cpz_0.80", &K_plus_cpz_0_80, &b_K_plus_cpz_0_80);
    fChain->SetBranchAddress("K_plus_cpt_0.80", &K_plus_cpt_0_80, &b_K_plus_cpt_0_80);
    fChain->SetBranchAddress("K_plus_cp_0.80", &K_plus_cp_0_80, &b_K_plus_cp_0_80);
    fChain->SetBranchAddress("K_plus_cmult_0.80", &K_plus_cmult_0_80, &b_K_plus_cmult_0_80);
    fChain->SetBranchAddress("K_plus_deltaEta_0.80", &K_plus_deltaEta_0_80, &b_K_plus_deltaEta_0_80);
    fChain->SetBranchAddress("K_plus_deltaPhi_0.80", &K_plus_deltaPhi_0_80, &b_K_plus_deltaPhi_0_80);
    fChain->SetBranchAddress("K_plus_pxasy_0.80", &K_plus_pxasy_0_80, &b_K_plus_pxasy_0_80);
    fChain->SetBranchAddress("K_plus_pyasy_0.80", &K_plus_pyasy_0_80, &b_K_plus_pyasy_0_80);
    fChain->SetBranchAddress("K_plus_pzasy_0.80", &K_plus_pzasy_0_80, &b_K_plus_pzasy_0_80);
    fChain->SetBranchAddress("K_plus_pasy_0.80", &K_plus_pasy_0_80, &b_K_plus_pasy_0_80);
    fChain->SetBranchAddress("K_plus_ptasy_0.80", &K_plus_ptasy_0_80, &b_K_plus_ptasy_0_80);
    fChain->SetBranchAddress("K_plus_cpx_0.90", &K_plus_cpx_0_90, &b_K_plus_cpx_0_90);
    fChain->SetBranchAddress("K_plus_cpy_0.90", &K_plus_cpy_0_90, &b_K_plus_cpy_0_90);
    fChain->SetBranchAddress("K_plus_cpz_0.90", &K_plus_cpz_0_90, &b_K_plus_cpz_0_90);
    fChain->SetBranchAddress("K_plus_cpt_0.90", &K_plus_cpt_0_90, &b_K_plus_cpt_0_90);
    fChain->SetBranchAddress("K_plus_cp_0.90", &K_plus_cp_0_90, &b_K_plus_cp_0_90);
    fChain->SetBranchAddress("K_plus_cmult_0.90", &K_plus_cmult_0_90, &b_K_plus_cmult_0_90);
    fChain->SetBranchAddress("K_plus_deltaEta_0.90", &K_plus_deltaEta_0_90, &b_K_plus_deltaEta_0_90);
    fChain->SetBranchAddress("K_plus_deltaPhi_0.90", &K_plus_deltaPhi_0_90, &b_K_plus_deltaPhi_0_90);
    fChain->SetBranchAddress("K_plus_pxasy_0.90", &K_plus_pxasy_0_90, &b_K_plus_pxasy_0_90);
    fChain->SetBranchAddress("K_plus_pyasy_0.90", &K_plus_pyasy_0_90, &b_K_plus_pyasy_0_90);
    fChain->SetBranchAddress("K_plus_pzasy_0.90", &K_plus_pzasy_0_90, &b_K_plus_pzasy_0_90);
    fChain->SetBranchAddress("K_plus_pasy_0.90", &K_plus_pasy_0_90, &b_K_plus_pasy_0_90);
    fChain->SetBranchAddress("K_plus_ptasy_0.90", &K_plus_ptasy_0_90, &b_K_plus_ptasy_0_90);
    fChain->SetBranchAddress("K_plus_cpx_1.00", &K_plus_cpx_1_00, &b_K_plus_cpx_1_00);
    fChain->SetBranchAddress("K_plus_cpy_1.00", &K_plus_cpy_1_00, &b_K_plus_cpy_1_00);
    fChain->SetBranchAddress("K_plus_cpz_1.00", &K_plus_cpz_1_00, &b_K_plus_cpz_1_00);
    fChain->SetBranchAddress("K_plus_cpt_1.00", &K_plus_cpt_1_00, &b_K_plus_cpt_1_00);
    fChain->SetBranchAddress("K_plus_cp_1.00", &K_plus_cp_1_00, &b_K_plus_cp_1_00);
    fChain->SetBranchAddress("K_plus_cmult_1.00", &K_plus_cmult_1_00, &b_K_plus_cmult_1_00);
    fChain->SetBranchAddress("K_plus_deltaEta_1.00", &K_plus_deltaEta_1_00, &b_K_plus_deltaEta_1_00);
    fChain->SetBranchAddress("K_plus_deltaPhi_1.00", &K_plus_deltaPhi_1_00, &b_K_plus_deltaPhi_1_00);
    fChain->SetBranchAddress("K_plus_pxasy_1.00", &K_plus_pxasy_1_00, &b_K_plus_pxasy_1_00);
    fChain->SetBranchAddress("K_plus_pyasy_1.00", &K_plus_pyasy_1_00, &b_K_plus_pyasy_1_00);
    fChain->SetBranchAddress("K_plus_pzasy_1.00", &K_plus_pzasy_1_00, &b_K_plus_pzasy_1_00);
    fChain->SetBranchAddress("K_plus_pasy_1.00", &K_plus_pasy_1_00, &b_K_plus_pasy_1_00);
    fChain->SetBranchAddress("K_plus_ptasy_1.00", &K_plus_ptasy_1_00, &b_K_plus_ptasy_1_00);
    fChain->SetBranchAddress("pi_plus_DOCA1", &pi_plus_DOCA1, &b_pi_plus_DOCA1);
    fChain->SetBranchAddress("pi_plus_DOCA2", &pi_plus_DOCA2, &b_pi_plus_DOCA2);
    fChain->SetBranchAddress("pi_plus_DOCA3", &pi_plus_DOCA3, &b_pi_plus_DOCA3);
    fChain->SetBranchAddress("pi_plus_ETA", &pi_plus_ETA, &b_pi_plus_ETA);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNe", &pi_plus_MC12TuneV2_ProbNNe, &b_pi_plus_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNmu", &pi_plus_MC12TuneV2_ProbNNmu, &b_pi_plus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNpi", &pi_plus_MC12TuneV2_ProbNNpi, &b_pi_plus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNk", &pi_plus_MC12TuneV2_ProbNNk, &b_pi_plus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNp", &pi_plus_MC12TuneV2_ProbNNp, &b_pi_plus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNghost", &pi_plus_MC12TuneV2_ProbNNghost, &b_pi_plus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNe", &pi_plus_MC12TuneV3_ProbNNe, &b_pi_plus_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNmu", &pi_plus_MC12TuneV3_ProbNNmu, &b_pi_plus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNpi", &pi_plus_MC12TuneV3_ProbNNpi, &b_pi_plus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNk", &pi_plus_MC12TuneV3_ProbNNk, &b_pi_plus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNp", &pi_plus_MC12TuneV3_ProbNNp, &b_pi_plus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNghost", &pi_plus_MC12TuneV3_ProbNNghost, &b_pi_plus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNe", &pi_plus_MC12TuneV4_ProbNNe, &b_pi_plus_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNmu", &pi_plus_MC12TuneV4_ProbNNmu, &b_pi_plus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNpi", &pi_plus_MC12TuneV4_ProbNNpi, &b_pi_plus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNk", &pi_plus_MC12TuneV4_ProbNNk, &b_pi_plus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNp", &pi_plus_MC12TuneV4_ProbNNp, &b_pi_plus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNghost", &pi_plus_MC12TuneV4_ProbNNghost, &b_pi_plus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNe", &pi_plus_MC15TuneV1_ProbNNe, &b_pi_plus_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNmu", &pi_plus_MC15TuneV1_ProbNNmu, &b_pi_plus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNpi", &pi_plus_MC15TuneV1_ProbNNpi, &b_pi_plus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNk", &pi_plus_MC15TuneV1_ProbNNk, &b_pi_plus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNp", &pi_plus_MC15TuneV1_ProbNNp, &b_pi_plus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNghost", &pi_plus_MC15TuneV1_ProbNNghost, &b_pi_plus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_CosTheta", &pi_plus_CosTheta, &b_pi_plus_CosTheta);
    fChain->SetBranchAddress("pi_plus_OWNPV_X", &pi_plus_OWNPV_X, &b_pi_plus_OWNPV_X);
    fChain->SetBranchAddress("pi_plus_OWNPV_Y", &pi_plus_OWNPV_Y, &b_pi_plus_OWNPV_Y);
    fChain->SetBranchAddress("pi_plus_OWNPV_Z", &pi_plus_OWNPV_Z, &b_pi_plus_OWNPV_Z);
    fChain->SetBranchAddress("pi_plus_OWNPV_XERR", &pi_plus_OWNPV_XERR, &b_pi_plus_OWNPV_XERR);
    fChain->SetBranchAddress("pi_plus_OWNPV_YERR", &pi_plus_OWNPV_YERR, &b_pi_plus_OWNPV_YERR);
    fChain->SetBranchAddress("pi_plus_OWNPV_ZERR", &pi_plus_OWNPV_ZERR, &b_pi_plus_OWNPV_ZERR);
    fChain->SetBranchAddress("pi_plus_OWNPV_CHI2", &pi_plus_OWNPV_CHI2, &b_pi_plus_OWNPV_CHI2);
    fChain->SetBranchAddress("pi_plus_OWNPV_NDOF", &pi_plus_OWNPV_NDOF, &b_pi_plus_OWNPV_NDOF);
    fChain->SetBranchAddress("pi_plus_OWNPV_COV_", pi_plus_OWNPV_COV_, &b_pi_plus_OWNPV_COV_);
    fChain->SetBranchAddress("pi_plus_IP_OWNPV", &pi_plus_IP_OWNPV, &b_pi_plus_IP_OWNPV);
    fChain->SetBranchAddress("pi_plus_IPCHI2_OWNPV", &pi_plus_IPCHI2_OWNPV, &b_pi_plus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_plus_ORIVX_X", &pi_plus_ORIVX_X, &b_pi_plus_ORIVX_X);
    fChain->SetBranchAddress("pi_plus_ORIVX_Y", &pi_plus_ORIVX_Y, &b_pi_plus_ORIVX_Y);
    fChain->SetBranchAddress("pi_plus_ORIVX_Z", &pi_plus_ORIVX_Z, &b_pi_plus_ORIVX_Z);
    fChain->SetBranchAddress("pi_plus_ORIVX_XERR", &pi_plus_ORIVX_XERR, &b_pi_plus_ORIVX_XERR);
    fChain->SetBranchAddress("pi_plus_ORIVX_YERR", &pi_plus_ORIVX_YERR, &b_pi_plus_ORIVX_YERR);
    fChain->SetBranchAddress("pi_plus_ORIVX_ZERR", &pi_plus_ORIVX_ZERR, &b_pi_plus_ORIVX_ZERR);
    fChain->SetBranchAddress("pi_plus_ORIVX_CHI2", &pi_plus_ORIVX_CHI2, &b_pi_plus_ORIVX_CHI2);
    fChain->SetBranchAddress("pi_plus_ORIVX_NDOF", &pi_plus_ORIVX_NDOF, &b_pi_plus_ORIVX_NDOF);
    fChain->SetBranchAddress("pi_plus_ORIVX_COV_", pi_plus_ORIVX_COV_, &b_pi_plus_ORIVX_COV_);
    fChain->SetBranchAddress("pi_plus_P", &pi_plus_P, &b_pi_plus_P);
    fChain->SetBranchAddress("pi_plus_PT", &pi_plus_PT, &b_pi_plus_PT);
    fChain->SetBranchAddress("pi_plus_PE", &pi_plus_PE, &b_pi_plus_PE);
    fChain->SetBranchAddress("pi_plus_PX", &pi_plus_PX, &b_pi_plus_PX);
    fChain->SetBranchAddress("pi_plus_PY", &pi_plus_PY, &b_pi_plus_PY);
    fChain->SetBranchAddress("pi_plus_PZ", &pi_plus_PZ, &b_pi_plus_PZ);
    fChain->SetBranchAddress("pi_plus_M", &pi_plus_M, &b_pi_plus_M);
    fChain->SetBranchAddress("pi_plus_ID", &pi_plus_ID, &b_pi_plus_ID);
    fChain->SetBranchAddress("pi_plus_PIDe", &pi_plus_PIDe, &b_pi_plus_PIDe);
    fChain->SetBranchAddress("pi_plus_PIDmu", &pi_plus_PIDmu, &b_pi_plus_PIDmu);
    fChain->SetBranchAddress("pi_plus_PIDK", &pi_plus_PIDK, &b_pi_plus_PIDK);
    fChain->SetBranchAddress("pi_plus_PIDp", &pi_plus_PIDp, &b_pi_plus_PIDp);
    fChain->SetBranchAddress("pi_plus_ProbNNe", &pi_plus_ProbNNe, &b_pi_plus_ProbNNe);
    fChain->SetBranchAddress("pi_plus_ProbNNk", &pi_plus_ProbNNk, &b_pi_plus_ProbNNk);
    fChain->SetBranchAddress("pi_plus_ProbNNp", &pi_plus_ProbNNp, &b_pi_plus_ProbNNp);
    fChain->SetBranchAddress("pi_plus_ProbNNpi", &pi_plus_ProbNNpi, &b_pi_plus_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_ProbNNmu", &pi_plus_ProbNNmu, &b_pi_plus_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_ProbNNghost", &pi_plus_ProbNNghost, &b_pi_plus_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_hasMuon", &pi_plus_hasMuon, &b_pi_plus_hasMuon);
    fChain->SetBranchAddress("pi_plus_isMuon", &pi_plus_isMuon, &b_pi_plus_isMuon);
    fChain->SetBranchAddress("pi_plus_hasRich", &pi_plus_hasRich, &b_pi_plus_hasRich);
    fChain->SetBranchAddress("pi_plus_UsedRichAerogel", &pi_plus_UsedRichAerogel, &b_pi_plus_UsedRichAerogel);
    fChain->SetBranchAddress("pi_plus_UsedRich1Gas", &pi_plus_UsedRich1Gas, &b_pi_plus_UsedRich1Gas);
    fChain->SetBranchAddress("pi_plus_UsedRich2Gas", &pi_plus_UsedRich2Gas, &b_pi_plus_UsedRich2Gas);
    fChain->SetBranchAddress("pi_plus_RichAboveElThres", &pi_plus_RichAboveElThres, &b_pi_plus_RichAboveElThres);
    fChain->SetBranchAddress("pi_plus_RichAboveMuThres", &pi_plus_RichAboveMuThres, &b_pi_plus_RichAboveMuThres);
    fChain->SetBranchAddress("pi_plus_RichAbovePiThres", &pi_plus_RichAbovePiThres, &b_pi_plus_RichAbovePiThres);
    fChain->SetBranchAddress("pi_plus_RichAboveKaThres", &pi_plus_RichAboveKaThres, &b_pi_plus_RichAboveKaThres);
    fChain->SetBranchAddress("pi_plus_RichAbovePrThres", &pi_plus_RichAbovePrThres, &b_pi_plus_RichAbovePrThres);
    fChain->SetBranchAddress("pi_plus_hasCalo", &pi_plus_hasCalo, &b_pi_plus_hasCalo);
    fChain->SetBranchAddress("pi_plus_L0Global_Dec", &pi_plus_L0Global_Dec, &b_pi_plus_L0Global_Dec);
    fChain->SetBranchAddress("pi_plus_L0Global_TIS", &pi_plus_L0Global_TIS, &b_pi_plus_L0Global_TIS);
    fChain->SetBranchAddress("pi_plus_L0Global_TOS", &pi_plus_L0Global_TOS, &b_pi_plus_L0Global_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1Global_Dec", &pi_plus_Hlt1Global_Dec, &b_pi_plus_Hlt1Global_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1Global_TIS", &pi_plus_Hlt1Global_TIS, &b_pi_plus_Hlt1Global_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1Global_TOS", &pi_plus_Hlt1Global_TOS, &b_pi_plus_Hlt1Global_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1Phys_Dec", &pi_plus_Hlt1Phys_Dec, &b_pi_plus_Hlt1Phys_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1Phys_TIS", &pi_plus_Hlt1Phys_TIS, &b_pi_plus_Hlt1Phys_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1Phys_TOS", &pi_plus_Hlt1Phys_TOS, &b_pi_plus_Hlt1Phys_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Global_Dec", &pi_plus_Hlt2Global_Dec, &b_pi_plus_Hlt2Global_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Global_TIS", &pi_plus_Hlt2Global_TIS, &b_pi_plus_Hlt2Global_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Global_TOS", &pi_plus_Hlt2Global_TOS, &b_pi_plus_Hlt2Global_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Phys_Dec", &pi_plus_Hlt2Phys_Dec, &b_pi_plus_Hlt2Phys_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Phys_TIS", &pi_plus_Hlt2Phys_TIS, &b_pi_plus_Hlt2Phys_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Phys_TOS", &pi_plus_Hlt2Phys_TOS, &b_pi_plus_Hlt2Phys_TOS);
    fChain->SetBranchAddress("pi_plus_L0HadronDecision_Dec", &pi_plus_L0HadronDecision_Dec, &b_pi_plus_L0HadronDecision_Dec);
    fChain->SetBranchAddress("pi_plus_L0HadronDecision_TIS", &pi_plus_L0HadronDecision_TIS, &b_pi_plus_L0HadronDecision_TIS);
    fChain->SetBranchAddress("pi_plus_L0HadronDecision_TOS", &pi_plus_L0HadronDecision_TOS, &b_pi_plus_L0HadronDecision_TOS);
    fChain->SetBranchAddress("pi_plus_L0MuonDecision_Dec", &pi_plus_L0MuonDecision_Dec, &b_pi_plus_L0MuonDecision_Dec);
    fChain->SetBranchAddress("pi_plus_L0MuonDecision_TIS", &pi_plus_L0MuonDecision_TIS, &b_pi_plus_L0MuonDecision_TIS);
    fChain->SetBranchAddress("pi_plus_L0MuonDecision_TOS", &pi_plus_L0MuonDecision_TOS, &b_pi_plus_L0MuonDecision_TOS);
    fChain->SetBranchAddress("pi_plus_L0GlobalDecision_Dec", &pi_plus_L0GlobalDecision_Dec, &b_pi_plus_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("pi_plus_L0GlobalDecision_TIS", &pi_plus_L0GlobalDecision_TIS, &b_pi_plus_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("pi_plus_L0GlobalDecision_TOS", &pi_plus_L0GlobalDecision_TOS, &b_pi_plus_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackAllL0Decision_Dec", &pi_plus_Hlt1TrackAllL0Decision_Dec, &b_pi_plus_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackAllL0Decision_TIS", &pi_plus_Hlt1TrackAllL0Decision_TIS, &b_pi_plus_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackAllL0Decision_TOS", &pi_plus_Hlt1TrackAllL0Decision_TOS, &b_pi_plus_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVADecision_Dec", &pi_plus_Hlt1TrackMVADecision_Dec, &b_pi_plus_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVADecision_TIS", &pi_plus_Hlt1TrackMVADecision_TIS, &b_pi_plus_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVADecision_TOS", &pi_plus_Hlt1TrackMVADecision_TOS, &b_pi_plus_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVADecision_Dec", &pi_plus_Hlt1TwoTrackMVADecision_Dec, &b_pi_plus_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVADecision_TIS", &pi_plus_Hlt1TwoTrackMVADecision_TIS, &b_pi_plus_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVADecision_TOS", &pi_plus_Hlt1TwoTrackMVADecision_TOS, &b_pi_plus_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVALooseDecision_Dec", &pi_plus_Hlt1TrackMVALooseDecision_Dec, &b_pi_plus_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVALooseDecision_TIS", &pi_plus_Hlt1TrackMVALooseDecision_TIS, &b_pi_plus_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1TrackMVALooseDecision_TOS", &pi_plus_Hlt1TrackMVALooseDecision_TOS, &b_pi_plus_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVALooseDecision_Dec", &pi_plus_Hlt1TwoTrackMVALooseDecision_Dec, &b_pi_plus_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVALooseDecision_TIS", &pi_plus_Hlt1TwoTrackMVALooseDecision_TIS, &b_pi_plus_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt1TwoTrackMVALooseDecision_TOS", &pi_plus_Hlt1TwoTrackMVALooseDecision_TOS, &b_pi_plus_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2IncPhiDecision_Dec", &pi_plus_Hlt2IncPhiDecision_Dec, &b_pi_plus_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2IncPhiDecision_TIS", &pi_plus_Hlt2IncPhiDecision_TIS, &b_pi_plus_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2IncPhiDecision_TOS", &pi_plus_Hlt2IncPhiDecision_TOS, &b_pi_plus_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2PhiIncPhiDecision_Dec", &pi_plus_Hlt2PhiIncPhiDecision_Dec, &b_pi_plus_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2PhiIncPhiDecision_TIS", &pi_plus_Hlt2PhiIncPhiDecision_TIS, &b_pi_plus_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2PhiIncPhiDecision_TOS", &pi_plus_Hlt2PhiIncPhiDecision_TOS, &b_pi_plus_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyBBDTDecision_Dec", &pi_plus_Hlt2Topo2BodyBBDTDecision_Dec, &b_pi_plus_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyBBDTDecision_TIS", &pi_plus_Hlt2Topo2BodyBBDTDecision_TIS, &b_pi_plus_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyBBDTDecision_TOS", &pi_plus_Hlt2Topo2BodyBBDTDecision_TOS, &b_pi_plus_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyBBDTDecision_Dec", &pi_plus_Hlt2Topo3BodyBBDTDecision_Dec, &b_pi_plus_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyBBDTDecision_TIS", &pi_plus_Hlt2Topo3BodyBBDTDecision_TIS, &b_pi_plus_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyBBDTDecision_TOS", &pi_plus_Hlt2Topo3BodyBBDTDecision_TOS, &b_pi_plus_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyBBDTDecision_Dec", &pi_plus_Hlt2Topo4BodyBBDTDecision_Dec, &b_pi_plus_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyBBDTDecision_TIS", &pi_plus_Hlt2Topo4BodyBBDTDecision_TIS, &b_pi_plus_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyBBDTDecision_TOS", &pi_plus_Hlt2Topo4BodyBBDTDecision_TOS, &b_pi_plus_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyDecision_Dec", &pi_plus_Hlt2Topo2BodyDecision_Dec, &b_pi_plus_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyDecision_TIS", &pi_plus_Hlt2Topo2BodyDecision_TIS, &b_pi_plus_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo2BodyDecision_TOS", &pi_plus_Hlt2Topo2BodyDecision_TOS, &b_pi_plus_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyDecision_Dec", &pi_plus_Hlt2Topo3BodyDecision_Dec, &b_pi_plus_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyDecision_TIS", &pi_plus_Hlt2Topo3BodyDecision_TIS, &b_pi_plus_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo3BodyDecision_TOS", &pi_plus_Hlt2Topo3BodyDecision_TOS, &b_pi_plus_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyDecision_Dec", &pi_plus_Hlt2Topo4BodyDecision_Dec, &b_pi_plus_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyDecision_TIS", &pi_plus_Hlt2Topo4BodyDecision_TIS, &b_pi_plus_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("pi_plus_Hlt2Topo4BodyDecision_TOS", &pi_plus_Hlt2Topo4BodyDecision_TOS, &b_pi_plus_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("pi_plus_TRACK_Type", &pi_plus_TRACK_Type, &b_pi_plus_TRACK_Type);
    fChain->SetBranchAddress("pi_plus_TRACK_Key", &pi_plus_TRACK_Key, &b_pi_plus_TRACK_Key);
    fChain->SetBranchAddress("pi_plus_TRACK_CHI2NDOF", &pi_plus_TRACK_CHI2NDOF, &b_pi_plus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_plus_TRACK_PCHI2", &pi_plus_TRACK_PCHI2, &b_pi_plus_TRACK_PCHI2);
    fChain->SetBranchAddress("pi_plus_TRACK_MatchCHI2", &pi_plus_TRACK_MatchCHI2, &b_pi_plus_TRACK_MatchCHI2);
    fChain->SetBranchAddress("pi_plus_TRACK_GhostProb", &pi_plus_TRACK_GhostProb, &b_pi_plus_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_plus_TRACK_CloneDist", &pi_plus_TRACK_CloneDist, &b_pi_plus_TRACK_CloneDist);
    fChain->SetBranchAddress("pi_plus_TRACK_Likelihood", &pi_plus_TRACK_Likelihood, &b_pi_plus_TRACK_Likelihood);
    fChain->SetBranchAddress("pi_plus_cpx_0.50", &pi_plus_cpx_0_50, &b_pi_plus_cpx_0_50);
    fChain->SetBranchAddress("pi_plus_cpy_0.50", &pi_plus_cpy_0_50, &b_pi_plus_cpy_0_50);
    fChain->SetBranchAddress("pi_plus_cpz_0.50", &pi_plus_cpz_0_50, &b_pi_plus_cpz_0_50);
    fChain->SetBranchAddress("pi_plus_cpt_0.50", &pi_plus_cpt_0_50, &b_pi_plus_cpt_0_50);
    fChain->SetBranchAddress("pi_plus_cp_0.50", &pi_plus_cp_0_50, &b_pi_plus_cp_0_50);
    fChain->SetBranchAddress("pi_plus_cmult_0.50", &pi_plus_cmult_0_50, &b_pi_plus_cmult_0_50);
    fChain->SetBranchAddress("pi_plus_deltaEta_0.50", &pi_plus_deltaEta_0_50, &b_pi_plus_deltaEta_0_50);
    fChain->SetBranchAddress("pi_plus_deltaPhi_0.50", &pi_plus_deltaPhi_0_50, &b_pi_plus_deltaPhi_0_50);
    fChain->SetBranchAddress("pi_plus_pxasy_0.50", &pi_plus_pxasy_0_50, &b_pi_plus_pxasy_0_50);
    fChain->SetBranchAddress("pi_plus_pyasy_0.50", &pi_plus_pyasy_0_50, &b_pi_plus_pyasy_0_50);
    fChain->SetBranchAddress("pi_plus_pzasy_0.50", &pi_plus_pzasy_0_50, &b_pi_plus_pzasy_0_50);
    fChain->SetBranchAddress("pi_plus_pasy_0.50", &pi_plus_pasy_0_50, &b_pi_plus_pasy_0_50);
    fChain->SetBranchAddress("pi_plus_ptasy_0.50", &pi_plus_ptasy_0_50, &b_pi_plus_ptasy_0_50);
    fChain->SetBranchAddress("pi_plus_cpx_0.60", &pi_plus_cpx_0_60, &b_pi_plus_cpx_0_60);
    fChain->SetBranchAddress("pi_plus_cpy_0.60", &pi_plus_cpy_0_60, &b_pi_plus_cpy_0_60);
    fChain->SetBranchAddress("pi_plus_cpz_0.60", &pi_plus_cpz_0_60, &b_pi_plus_cpz_0_60);
    fChain->SetBranchAddress("pi_plus_cpt_0.60", &pi_plus_cpt_0_60, &b_pi_plus_cpt_0_60);
    fChain->SetBranchAddress("pi_plus_cp_0.60", &pi_plus_cp_0_60, &b_pi_plus_cp_0_60);
    fChain->SetBranchAddress("pi_plus_cmult_0.60", &pi_plus_cmult_0_60, &b_pi_plus_cmult_0_60);
    fChain->SetBranchAddress("pi_plus_deltaEta_0.60", &pi_plus_deltaEta_0_60, &b_pi_plus_deltaEta_0_60);
    fChain->SetBranchAddress("pi_plus_deltaPhi_0.60", &pi_plus_deltaPhi_0_60, &b_pi_plus_deltaPhi_0_60);
    fChain->SetBranchAddress("pi_plus_pxasy_0.60", &pi_plus_pxasy_0_60, &b_pi_plus_pxasy_0_60);
    fChain->SetBranchAddress("pi_plus_pyasy_0.60", &pi_plus_pyasy_0_60, &b_pi_plus_pyasy_0_60);
    fChain->SetBranchAddress("pi_plus_pzasy_0.60", &pi_plus_pzasy_0_60, &b_pi_plus_pzasy_0_60);
    fChain->SetBranchAddress("pi_plus_pasy_0.60", &pi_plus_pasy_0_60, &b_pi_plus_pasy_0_60);
    fChain->SetBranchAddress("pi_plus_ptasy_0.60", &pi_plus_ptasy_0_60, &b_pi_plus_ptasy_0_60);
    fChain->SetBranchAddress("pi_plus_cpx_0.70", &pi_plus_cpx_0_70, &b_pi_plus_cpx_0_70);
    fChain->SetBranchAddress("pi_plus_cpy_0.70", &pi_plus_cpy_0_70, &b_pi_plus_cpy_0_70);
    fChain->SetBranchAddress("pi_plus_cpz_0.70", &pi_plus_cpz_0_70, &b_pi_plus_cpz_0_70);
    fChain->SetBranchAddress("pi_plus_cpt_0.70", &pi_plus_cpt_0_70, &b_pi_plus_cpt_0_70);
    fChain->SetBranchAddress("pi_plus_cp_0.70", &pi_plus_cp_0_70, &b_pi_plus_cp_0_70);
    fChain->SetBranchAddress("pi_plus_cmult_0.70", &pi_plus_cmult_0_70, &b_pi_plus_cmult_0_70);
    fChain->SetBranchAddress("pi_plus_deltaEta_0.70", &pi_plus_deltaEta_0_70, &b_pi_plus_deltaEta_0_70);
    fChain->SetBranchAddress("pi_plus_deltaPhi_0.70", &pi_plus_deltaPhi_0_70, &b_pi_plus_deltaPhi_0_70);
    fChain->SetBranchAddress("pi_plus_pxasy_0.70", &pi_plus_pxasy_0_70, &b_pi_plus_pxasy_0_70);
    fChain->SetBranchAddress("pi_plus_pyasy_0.70", &pi_plus_pyasy_0_70, &b_pi_plus_pyasy_0_70);
    fChain->SetBranchAddress("pi_plus_pzasy_0.70", &pi_plus_pzasy_0_70, &b_pi_plus_pzasy_0_70);
    fChain->SetBranchAddress("pi_plus_pasy_0.70", &pi_plus_pasy_0_70, &b_pi_plus_pasy_0_70);
    fChain->SetBranchAddress("pi_plus_ptasy_0.70", &pi_plus_ptasy_0_70, &b_pi_plus_ptasy_0_70);
    fChain->SetBranchAddress("pi_plus_cpx_0.80", &pi_plus_cpx_0_80, &b_pi_plus_cpx_0_80);
    fChain->SetBranchAddress("pi_plus_cpy_0.80", &pi_plus_cpy_0_80, &b_pi_plus_cpy_0_80);
    fChain->SetBranchAddress("pi_plus_cpz_0.80", &pi_plus_cpz_0_80, &b_pi_plus_cpz_0_80);
    fChain->SetBranchAddress("pi_plus_cpt_0.80", &pi_plus_cpt_0_80, &b_pi_plus_cpt_0_80);
    fChain->SetBranchAddress("pi_plus_cp_0.80", &pi_plus_cp_0_80, &b_pi_plus_cp_0_80);
    fChain->SetBranchAddress("pi_plus_cmult_0.80", &pi_plus_cmult_0_80, &b_pi_plus_cmult_0_80);
    fChain->SetBranchAddress("pi_plus_deltaEta_0.80", &pi_plus_deltaEta_0_80, &b_pi_plus_deltaEta_0_80);
    fChain->SetBranchAddress("pi_plus_deltaPhi_0.80", &pi_plus_deltaPhi_0_80, &b_pi_plus_deltaPhi_0_80);
    fChain->SetBranchAddress("pi_plus_pxasy_0.80", &pi_plus_pxasy_0_80, &b_pi_plus_pxasy_0_80);
    fChain->SetBranchAddress("pi_plus_pyasy_0.80", &pi_plus_pyasy_0_80, &b_pi_plus_pyasy_0_80);
    fChain->SetBranchAddress("pi_plus_pzasy_0.80", &pi_plus_pzasy_0_80, &b_pi_plus_pzasy_0_80);
    fChain->SetBranchAddress("pi_plus_pasy_0.80", &pi_plus_pasy_0_80, &b_pi_plus_pasy_0_80);
    fChain->SetBranchAddress("pi_plus_ptasy_0.80", &pi_plus_ptasy_0_80, &b_pi_plus_ptasy_0_80);
    fChain->SetBranchAddress("pi_plus_cpx_0.90", &pi_plus_cpx_0_90, &b_pi_plus_cpx_0_90);
    fChain->SetBranchAddress("pi_plus_cpy_0.90", &pi_plus_cpy_0_90, &b_pi_plus_cpy_0_90);
    fChain->SetBranchAddress("pi_plus_cpz_0.90", &pi_plus_cpz_0_90, &b_pi_plus_cpz_0_90);
    fChain->SetBranchAddress("pi_plus_cpt_0.90", &pi_plus_cpt_0_90, &b_pi_plus_cpt_0_90);
    fChain->SetBranchAddress("pi_plus_cp_0.90", &pi_plus_cp_0_90, &b_pi_plus_cp_0_90);
    fChain->SetBranchAddress("pi_plus_cmult_0.90", &pi_plus_cmult_0_90, &b_pi_plus_cmult_0_90);
    fChain->SetBranchAddress("pi_plus_deltaEta_0.90", &pi_plus_deltaEta_0_90, &b_pi_plus_deltaEta_0_90);
    fChain->SetBranchAddress("pi_plus_deltaPhi_0.90", &pi_plus_deltaPhi_0_90, &b_pi_plus_deltaPhi_0_90);
    fChain->SetBranchAddress("pi_plus_pxasy_0.90", &pi_plus_pxasy_0_90, &b_pi_plus_pxasy_0_90);
    fChain->SetBranchAddress("pi_plus_pyasy_0.90", &pi_plus_pyasy_0_90, &b_pi_plus_pyasy_0_90);
    fChain->SetBranchAddress("pi_plus_pzasy_0.90", &pi_plus_pzasy_0_90, &b_pi_plus_pzasy_0_90);
    fChain->SetBranchAddress("pi_plus_pasy_0.90", &pi_plus_pasy_0_90, &b_pi_plus_pasy_0_90);
    fChain->SetBranchAddress("pi_plus_ptasy_0.90", &pi_plus_ptasy_0_90, &b_pi_plus_ptasy_0_90);
    fChain->SetBranchAddress("pi_plus_cpx_1.00", &pi_plus_cpx_1_00, &b_pi_plus_cpx_1_00);
    fChain->SetBranchAddress("pi_plus_cpy_1.00", &pi_plus_cpy_1_00, &b_pi_plus_cpy_1_00);
    fChain->SetBranchAddress("pi_plus_cpz_1.00", &pi_plus_cpz_1_00, &b_pi_plus_cpz_1_00);
    fChain->SetBranchAddress("pi_plus_cpt_1.00", &pi_plus_cpt_1_00, &b_pi_plus_cpt_1_00);
    fChain->SetBranchAddress("pi_plus_cp_1.00", &pi_plus_cp_1_00, &b_pi_plus_cp_1_00);
    fChain->SetBranchAddress("pi_plus_cmult_1.00", &pi_plus_cmult_1_00, &b_pi_plus_cmult_1_00);
    fChain->SetBranchAddress("pi_plus_deltaEta_1.00", &pi_plus_deltaEta_1_00, &b_pi_plus_deltaEta_1_00);
    fChain->SetBranchAddress("pi_plus_deltaPhi_1.00", &pi_plus_deltaPhi_1_00, &b_pi_plus_deltaPhi_1_00);
    fChain->SetBranchAddress("pi_plus_pxasy_1.00", &pi_plus_pxasy_1_00, &b_pi_plus_pxasy_1_00);
    fChain->SetBranchAddress("pi_plus_pyasy_1.00", &pi_plus_pyasy_1_00, &b_pi_plus_pyasy_1_00);
    fChain->SetBranchAddress("pi_plus_pzasy_1.00", &pi_plus_pzasy_1_00, &b_pi_plus_pzasy_1_00);
    fChain->SetBranchAddress("pi_plus_pasy_1.00", &pi_plus_pasy_1_00, &b_pi_plus_pasy_1_00);
    fChain->SetBranchAddress("pi_plus_ptasy_1.00", &pi_plus_ptasy_1_00, &b_pi_plus_ptasy_1_00);
    fChain->SetBranchAddress("pi_minus_DOCA1", &pi_minus_DOCA1, &b_pi_minus_DOCA1);
    fChain->SetBranchAddress("pi_minus_DOCA2", &pi_minus_DOCA2, &b_pi_minus_DOCA2);
    fChain->SetBranchAddress("pi_minus_DOCA3", &pi_minus_DOCA3, &b_pi_minus_DOCA3);
    fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNe", &pi_minus_MC12TuneV2_ProbNNe, &b_pi_minus_MC12TuneV2_ProbNNe);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNe", &pi_minus_MC12TuneV3_ProbNNe, &b_pi_minus_MC12TuneV3_ProbNNe);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNe", &pi_minus_MC12TuneV4_ProbNNe, &b_pi_minus_MC12TuneV4_ProbNNe);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNmu", &pi_minus_MC12TuneV4_ProbNNmu, &b_pi_minus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNpi", &pi_minus_MC12TuneV4_ProbNNpi, &b_pi_minus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNk", &pi_minus_MC12TuneV4_ProbNNk, &b_pi_minus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNp", &pi_minus_MC12TuneV4_ProbNNp, &b_pi_minus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNghost", &pi_minus_MC12TuneV4_ProbNNghost, &b_pi_minus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNe", &pi_minus_MC15TuneV1_ProbNNe, &b_pi_minus_MC15TuneV1_ProbNNe);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNmu", &pi_minus_MC15TuneV1_ProbNNmu, &b_pi_minus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNpi", &pi_minus_MC15TuneV1_ProbNNpi, &b_pi_minus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNk", &pi_minus_MC15TuneV1_ProbNNk, &b_pi_minus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNp", &pi_minus_MC15TuneV1_ProbNNp, &b_pi_minus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNghost", &pi_minus_MC15TuneV1_ProbNNghost, &b_pi_minus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_CosTheta", &pi_minus_CosTheta, &b_pi_minus_CosTheta);
    fChain->SetBranchAddress("pi_minus_OWNPV_X", &pi_minus_OWNPV_X, &b_pi_minus_OWNPV_X);
    fChain->SetBranchAddress("pi_minus_OWNPV_Y", &pi_minus_OWNPV_Y, &b_pi_minus_OWNPV_Y);
    fChain->SetBranchAddress("pi_minus_OWNPV_Z", &pi_minus_OWNPV_Z, &b_pi_minus_OWNPV_Z);
    fChain->SetBranchAddress("pi_minus_OWNPV_XERR", &pi_minus_OWNPV_XERR, &b_pi_minus_OWNPV_XERR);
    fChain->SetBranchAddress("pi_minus_OWNPV_YERR", &pi_minus_OWNPV_YERR, &b_pi_minus_OWNPV_YERR);
    fChain->SetBranchAddress("pi_minus_OWNPV_ZERR", &pi_minus_OWNPV_ZERR, &b_pi_minus_OWNPV_ZERR);
    fChain->SetBranchAddress("pi_minus_OWNPV_CHI2", &pi_minus_OWNPV_CHI2, &b_pi_minus_OWNPV_CHI2);
    fChain->SetBranchAddress("pi_minus_OWNPV_NDOF", &pi_minus_OWNPV_NDOF, &b_pi_minus_OWNPV_NDOF);
    fChain->SetBranchAddress("pi_minus_OWNPV_COV_", pi_minus_OWNPV_COV_, &b_pi_minus_OWNPV_COV_);
    fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
    fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_minus_ORIVX_X", &pi_minus_ORIVX_X, &b_pi_minus_ORIVX_X);
    fChain->SetBranchAddress("pi_minus_ORIVX_Y", &pi_minus_ORIVX_Y, &b_pi_minus_ORIVX_Y);
    fChain->SetBranchAddress("pi_minus_ORIVX_Z", &pi_minus_ORIVX_Z, &b_pi_minus_ORIVX_Z);
    fChain->SetBranchAddress("pi_minus_ORIVX_XERR", &pi_minus_ORIVX_XERR, &b_pi_minus_ORIVX_XERR);
    fChain->SetBranchAddress("pi_minus_ORIVX_YERR", &pi_minus_ORIVX_YERR, &b_pi_minus_ORIVX_YERR);
    fChain->SetBranchAddress("pi_minus_ORIVX_ZERR", &pi_minus_ORIVX_ZERR, &b_pi_minus_ORIVX_ZERR);
    fChain->SetBranchAddress("pi_minus_ORIVX_CHI2", &pi_minus_ORIVX_CHI2, &b_pi_minus_ORIVX_CHI2);
    fChain->SetBranchAddress("pi_minus_ORIVX_NDOF", &pi_minus_ORIVX_NDOF, &b_pi_minus_ORIVX_NDOF);
    fChain->SetBranchAddress("pi_minus_ORIVX_COV_", pi_minus_ORIVX_COV_, &b_pi_minus_ORIVX_COV_);
    fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
    fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
    fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
    fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
    fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
    fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
    fChain->SetBranchAddress("pi_minus_M", &pi_minus_M, &b_pi_minus_M);
    fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
    fChain->SetBranchAddress("pi_minus_PIDe", &pi_minus_PIDe, &b_pi_minus_PIDe);
    fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
    fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
    fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
    fChain->SetBranchAddress("pi_minus_ProbNNe", &pi_minus_ProbNNe, &b_pi_minus_ProbNNe);
    fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
    fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
    fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_hasMuon", &pi_minus_hasMuon, &b_pi_minus_hasMuon);
    fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
    fChain->SetBranchAddress("pi_minus_hasRich", &pi_minus_hasRich, &b_pi_minus_hasRich);
    fChain->SetBranchAddress("pi_minus_UsedRichAerogel", &pi_minus_UsedRichAerogel, &b_pi_minus_UsedRichAerogel);
    fChain->SetBranchAddress("pi_minus_UsedRich1Gas", &pi_minus_UsedRich1Gas, &b_pi_minus_UsedRich1Gas);
    fChain->SetBranchAddress("pi_minus_UsedRich2Gas", &pi_minus_UsedRich2Gas, &b_pi_minus_UsedRich2Gas);
    fChain->SetBranchAddress("pi_minus_RichAboveElThres", &pi_minus_RichAboveElThres, &b_pi_minus_RichAboveElThres);
    fChain->SetBranchAddress("pi_minus_RichAboveMuThres", &pi_minus_RichAboveMuThres, &b_pi_minus_RichAboveMuThres);
    fChain->SetBranchAddress("pi_minus_RichAbovePiThres", &pi_minus_RichAbovePiThres, &b_pi_minus_RichAbovePiThres);
    fChain->SetBranchAddress("pi_minus_RichAboveKaThres", &pi_minus_RichAboveKaThres, &b_pi_minus_RichAboveKaThres);
    fChain->SetBranchAddress("pi_minus_RichAbovePrThres", &pi_minus_RichAbovePrThres, &b_pi_minus_RichAbovePrThres);
    fChain->SetBranchAddress("pi_minus_hasCalo", &pi_minus_hasCalo, &b_pi_minus_hasCalo);
    fChain->SetBranchAddress("pi_minus_L0Global_Dec", &pi_minus_L0Global_Dec, &b_pi_minus_L0Global_Dec);
    fChain->SetBranchAddress("pi_minus_L0Global_TIS", &pi_minus_L0Global_TIS, &b_pi_minus_L0Global_TIS);
    fChain->SetBranchAddress("pi_minus_L0Global_TOS", &pi_minus_L0Global_TOS, &b_pi_minus_L0Global_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1Global_Dec", &pi_minus_Hlt1Global_Dec, &b_pi_minus_Hlt1Global_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1Global_TIS", &pi_minus_Hlt1Global_TIS, &b_pi_minus_Hlt1Global_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1Global_TOS", &pi_minus_Hlt1Global_TOS, &b_pi_minus_Hlt1Global_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1Phys_Dec", &pi_minus_Hlt1Phys_Dec, &b_pi_minus_Hlt1Phys_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1Phys_TIS", &pi_minus_Hlt1Phys_TIS, &b_pi_minus_Hlt1Phys_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1Phys_TOS", &pi_minus_Hlt1Phys_TOS, &b_pi_minus_Hlt1Phys_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Global_Dec", &pi_minus_Hlt2Global_Dec, &b_pi_minus_Hlt2Global_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Global_TIS", &pi_minus_Hlt2Global_TIS, &b_pi_minus_Hlt2Global_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Global_TOS", &pi_minus_Hlt2Global_TOS, &b_pi_minus_Hlt2Global_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Phys_Dec", &pi_minus_Hlt2Phys_Dec, &b_pi_minus_Hlt2Phys_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Phys_TIS", &pi_minus_Hlt2Phys_TIS, &b_pi_minus_Hlt2Phys_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Phys_TOS", &pi_minus_Hlt2Phys_TOS, &b_pi_minus_Hlt2Phys_TOS);
    fChain->SetBranchAddress("pi_minus_L0HadronDecision_Dec", &pi_minus_L0HadronDecision_Dec, &b_pi_minus_L0HadronDecision_Dec);
    fChain->SetBranchAddress("pi_minus_L0HadronDecision_TIS", &pi_minus_L0HadronDecision_TIS, &b_pi_minus_L0HadronDecision_TIS);
    fChain->SetBranchAddress("pi_minus_L0HadronDecision_TOS", &pi_minus_L0HadronDecision_TOS, &b_pi_minus_L0HadronDecision_TOS);
    fChain->SetBranchAddress("pi_minus_L0MuonDecision_Dec", &pi_minus_L0MuonDecision_Dec, &b_pi_minus_L0MuonDecision_Dec);
    fChain->SetBranchAddress("pi_minus_L0MuonDecision_TIS", &pi_minus_L0MuonDecision_TIS, &b_pi_minus_L0MuonDecision_TIS);
    fChain->SetBranchAddress("pi_minus_L0MuonDecision_TOS", &pi_minus_L0MuonDecision_TOS, &b_pi_minus_L0MuonDecision_TOS);
    fChain->SetBranchAddress("pi_minus_L0GlobalDecision_Dec", &pi_minus_L0GlobalDecision_Dec, &b_pi_minus_L0GlobalDecision_Dec);
    fChain->SetBranchAddress("pi_minus_L0GlobalDecision_TIS", &pi_minus_L0GlobalDecision_TIS, &b_pi_minus_L0GlobalDecision_TIS);
    fChain->SetBranchAddress("pi_minus_L0GlobalDecision_TOS", &pi_minus_L0GlobalDecision_TOS, &b_pi_minus_L0GlobalDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackAllL0Decision_Dec", &pi_minus_Hlt1TrackAllL0Decision_Dec, &b_pi_minus_Hlt1TrackAllL0Decision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackAllL0Decision_TIS", &pi_minus_Hlt1TrackAllL0Decision_TIS, &b_pi_minus_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackAllL0Decision_TOS", &pi_minus_Hlt1TrackAllL0Decision_TOS, &b_pi_minus_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVADecision_Dec", &pi_minus_Hlt1TrackMVADecision_Dec, &b_pi_minus_Hlt1TrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVADecision_TIS", &pi_minus_Hlt1TrackMVADecision_TIS, &b_pi_minus_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVADecision_TOS", &pi_minus_Hlt1TrackMVADecision_TOS, &b_pi_minus_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVADecision_Dec", &pi_minus_Hlt1TwoTrackMVADecision_Dec, &b_pi_minus_Hlt1TwoTrackMVADecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVADecision_TIS", &pi_minus_Hlt1TwoTrackMVADecision_TIS, &b_pi_minus_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVADecision_TOS", &pi_minus_Hlt1TwoTrackMVADecision_TOS, &b_pi_minus_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVALooseDecision_Dec", &pi_minus_Hlt1TrackMVALooseDecision_Dec, &b_pi_minus_Hlt1TrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVALooseDecision_TIS", &pi_minus_Hlt1TrackMVALooseDecision_TIS, &b_pi_minus_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1TrackMVALooseDecision_TOS", &pi_minus_Hlt1TrackMVALooseDecision_TOS, &b_pi_minus_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVALooseDecision_Dec", &pi_minus_Hlt1TwoTrackMVALooseDecision_Dec, &b_pi_minus_Hlt1TwoTrackMVALooseDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVALooseDecision_TIS", &pi_minus_Hlt1TwoTrackMVALooseDecision_TIS, &b_pi_minus_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt1TwoTrackMVALooseDecision_TOS", &pi_minus_Hlt1TwoTrackMVALooseDecision_TOS, &b_pi_minus_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2IncPhiDecision_Dec", &pi_minus_Hlt2IncPhiDecision_Dec, &b_pi_minus_Hlt2IncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2IncPhiDecision_TIS", &pi_minus_Hlt2IncPhiDecision_TIS, &b_pi_minus_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2IncPhiDecision_TOS", &pi_minus_Hlt2IncPhiDecision_TOS, &b_pi_minus_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2PhiIncPhiDecision_Dec", &pi_minus_Hlt2PhiIncPhiDecision_Dec, &b_pi_minus_Hlt2PhiIncPhiDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2PhiIncPhiDecision_TIS", &pi_minus_Hlt2PhiIncPhiDecision_TIS, &b_pi_minus_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2PhiIncPhiDecision_TOS", &pi_minus_Hlt2PhiIncPhiDecision_TOS, &b_pi_minus_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyBBDTDecision_Dec", &pi_minus_Hlt2Topo2BodyBBDTDecision_Dec, &b_pi_minus_Hlt2Topo2BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyBBDTDecision_TIS", &pi_minus_Hlt2Topo2BodyBBDTDecision_TIS, &b_pi_minus_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyBBDTDecision_TOS", &pi_minus_Hlt2Topo2BodyBBDTDecision_TOS, &b_pi_minus_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyBBDTDecision_Dec", &pi_minus_Hlt2Topo3BodyBBDTDecision_Dec, &b_pi_minus_Hlt2Topo3BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyBBDTDecision_TIS", &pi_minus_Hlt2Topo3BodyBBDTDecision_TIS, &b_pi_minus_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyBBDTDecision_TOS", &pi_minus_Hlt2Topo3BodyBBDTDecision_TOS, &b_pi_minus_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyBBDTDecision_Dec", &pi_minus_Hlt2Topo4BodyBBDTDecision_Dec, &b_pi_minus_Hlt2Topo4BodyBBDTDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyBBDTDecision_TIS", &pi_minus_Hlt2Topo4BodyBBDTDecision_TIS, &b_pi_minus_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyBBDTDecision_TOS", &pi_minus_Hlt2Topo4BodyBBDTDecision_TOS, &b_pi_minus_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyDecision_Dec", &pi_minus_Hlt2Topo2BodyDecision_Dec, &b_pi_minus_Hlt2Topo2BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyDecision_TIS", &pi_minus_Hlt2Topo2BodyDecision_TIS, &b_pi_minus_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo2BodyDecision_TOS", &pi_minus_Hlt2Topo2BodyDecision_TOS, &b_pi_minus_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyDecision_Dec", &pi_minus_Hlt2Topo3BodyDecision_Dec, &b_pi_minus_Hlt2Topo3BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyDecision_TIS", &pi_minus_Hlt2Topo3BodyDecision_TIS, &b_pi_minus_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo3BodyDecision_TOS", &pi_minus_Hlt2Topo3BodyDecision_TOS, &b_pi_minus_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyDecision_Dec", &pi_minus_Hlt2Topo4BodyDecision_Dec, &b_pi_minus_Hlt2Topo4BodyDecision_Dec);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyDecision_TIS", &pi_minus_Hlt2Topo4BodyDecision_TIS, &b_pi_minus_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("pi_minus_Hlt2Topo4BodyDecision_TOS", &pi_minus_Hlt2Topo4BodyDecision_TOS, &b_pi_minus_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("pi_minus_TRACK_Type", &pi_minus_TRACK_Type, &b_pi_minus_TRACK_Type);
    fChain->SetBranchAddress("pi_minus_TRACK_Key", &pi_minus_TRACK_Key, &b_pi_minus_TRACK_Key);
    fChain->SetBranchAddress("pi_minus_TRACK_CHI2NDOF", &pi_minus_TRACK_CHI2NDOF, &b_pi_minus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_minus_TRACK_PCHI2", &pi_minus_TRACK_PCHI2, &b_pi_minus_TRACK_PCHI2);
    fChain->SetBranchAddress("pi_minus_TRACK_MatchCHI2", &pi_minus_TRACK_MatchCHI2, &b_pi_minus_TRACK_MatchCHI2);
    fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_minus_TRACK_CloneDist", &pi_minus_TRACK_CloneDist, &b_pi_minus_TRACK_CloneDist);
    fChain->SetBranchAddress("pi_minus_TRACK_Likelihood", &pi_minus_TRACK_Likelihood, &b_pi_minus_TRACK_Likelihood);
    fChain->SetBranchAddress("pi_minus_cpx_0.50", &pi_minus_cpx_0_50, &b_pi_minus_cpx_0_50);
    fChain->SetBranchAddress("pi_minus_cpy_0.50", &pi_minus_cpy_0_50, &b_pi_minus_cpy_0_50);
    fChain->SetBranchAddress("pi_minus_cpz_0.50", &pi_minus_cpz_0_50, &b_pi_minus_cpz_0_50);
    fChain->SetBranchAddress("pi_minus_cpt_0.50", &pi_minus_cpt_0_50, &b_pi_minus_cpt_0_50);
    fChain->SetBranchAddress("pi_minus_cp_0.50", &pi_minus_cp_0_50, &b_pi_minus_cp_0_50);
    fChain->SetBranchAddress("pi_minus_cmult_0.50", &pi_minus_cmult_0_50, &b_pi_minus_cmult_0_50);
    fChain->SetBranchAddress("pi_minus_deltaEta_0.50", &pi_minus_deltaEta_0_50, &b_pi_minus_deltaEta_0_50);
    fChain->SetBranchAddress("pi_minus_deltaPhi_0.50", &pi_minus_deltaPhi_0_50, &b_pi_minus_deltaPhi_0_50);
    fChain->SetBranchAddress("pi_minus_pxasy_0.50", &pi_minus_pxasy_0_50, &b_pi_minus_pxasy_0_50);
    fChain->SetBranchAddress("pi_minus_pyasy_0.50", &pi_minus_pyasy_0_50, &b_pi_minus_pyasy_0_50);
    fChain->SetBranchAddress("pi_minus_pzasy_0.50", &pi_minus_pzasy_0_50, &b_pi_minus_pzasy_0_50);
    fChain->SetBranchAddress("pi_minus_pasy_0.50", &pi_minus_pasy_0_50, &b_pi_minus_pasy_0_50);
    fChain->SetBranchAddress("pi_minus_ptasy_0.50", &pi_minus_ptasy_0_50, &b_pi_minus_ptasy_0_50);
    fChain->SetBranchAddress("pi_minus_cpx_0.60", &pi_minus_cpx_0_60, &b_pi_minus_cpx_0_60);
    fChain->SetBranchAddress("pi_minus_cpy_0.60", &pi_minus_cpy_0_60, &b_pi_minus_cpy_0_60);
    fChain->SetBranchAddress("pi_minus_cpz_0.60", &pi_minus_cpz_0_60, &b_pi_minus_cpz_0_60);
    fChain->SetBranchAddress("pi_minus_cpt_0.60", &pi_minus_cpt_0_60, &b_pi_minus_cpt_0_60);
    fChain->SetBranchAddress("pi_minus_cp_0.60", &pi_minus_cp_0_60, &b_pi_minus_cp_0_60);
    fChain->SetBranchAddress("pi_minus_cmult_0.60", &pi_minus_cmult_0_60, &b_pi_minus_cmult_0_60);
    fChain->SetBranchAddress("pi_minus_deltaEta_0.60", &pi_minus_deltaEta_0_60, &b_pi_minus_deltaEta_0_60);
    fChain->SetBranchAddress("pi_minus_deltaPhi_0.60", &pi_minus_deltaPhi_0_60, &b_pi_minus_deltaPhi_0_60);
    fChain->SetBranchAddress("pi_minus_pxasy_0.60", &pi_minus_pxasy_0_60, &b_pi_minus_pxasy_0_60);
    fChain->SetBranchAddress("pi_minus_pyasy_0.60", &pi_minus_pyasy_0_60, &b_pi_minus_pyasy_0_60);
    fChain->SetBranchAddress("pi_minus_pzasy_0.60", &pi_minus_pzasy_0_60, &b_pi_minus_pzasy_0_60);
    fChain->SetBranchAddress("pi_minus_pasy_0.60", &pi_minus_pasy_0_60, &b_pi_minus_pasy_0_60);
    fChain->SetBranchAddress("pi_minus_ptasy_0.60", &pi_minus_ptasy_0_60, &b_pi_minus_ptasy_0_60);
    fChain->SetBranchAddress("pi_minus_cpx_0.70", &pi_minus_cpx_0_70, &b_pi_minus_cpx_0_70);
    fChain->SetBranchAddress("pi_minus_cpy_0.70", &pi_minus_cpy_0_70, &b_pi_minus_cpy_0_70);
    fChain->SetBranchAddress("pi_minus_cpz_0.70", &pi_minus_cpz_0_70, &b_pi_minus_cpz_0_70);
    fChain->SetBranchAddress("pi_minus_cpt_0.70", &pi_minus_cpt_0_70, &b_pi_minus_cpt_0_70);
    fChain->SetBranchAddress("pi_minus_cp_0.70", &pi_minus_cp_0_70, &b_pi_minus_cp_0_70);
    fChain->SetBranchAddress("pi_minus_cmult_0.70", &pi_minus_cmult_0_70, &b_pi_minus_cmult_0_70);
    fChain->SetBranchAddress("pi_minus_deltaEta_0.70", &pi_minus_deltaEta_0_70, &b_pi_minus_deltaEta_0_70);
    fChain->SetBranchAddress("pi_minus_deltaPhi_0.70", &pi_minus_deltaPhi_0_70, &b_pi_minus_deltaPhi_0_70);
    fChain->SetBranchAddress("pi_minus_pxasy_0.70", &pi_minus_pxasy_0_70, &b_pi_minus_pxasy_0_70);
    fChain->SetBranchAddress("pi_minus_pyasy_0.70", &pi_minus_pyasy_0_70, &b_pi_minus_pyasy_0_70);
    fChain->SetBranchAddress("pi_minus_pzasy_0.70", &pi_minus_pzasy_0_70, &b_pi_minus_pzasy_0_70);
    fChain->SetBranchAddress("pi_minus_pasy_0.70", &pi_minus_pasy_0_70, &b_pi_minus_pasy_0_70);
    fChain->SetBranchAddress("pi_minus_ptasy_0.70", &pi_minus_ptasy_0_70, &b_pi_minus_ptasy_0_70);
    fChain->SetBranchAddress("pi_minus_cpx_0.80", &pi_minus_cpx_0_80, &b_pi_minus_cpx_0_80);
    fChain->SetBranchAddress("pi_minus_cpy_0.80", &pi_minus_cpy_0_80, &b_pi_minus_cpy_0_80);
    fChain->SetBranchAddress("pi_minus_cpz_0.80", &pi_minus_cpz_0_80, &b_pi_minus_cpz_0_80);
    fChain->SetBranchAddress("pi_minus_cpt_0.80", &pi_minus_cpt_0_80, &b_pi_minus_cpt_0_80);
    fChain->SetBranchAddress("pi_minus_cp_0.80", &pi_minus_cp_0_80, &b_pi_minus_cp_0_80);
    fChain->SetBranchAddress("pi_minus_cmult_0.80", &pi_minus_cmult_0_80, &b_pi_minus_cmult_0_80);
    fChain->SetBranchAddress("pi_minus_deltaEta_0.80", &pi_minus_deltaEta_0_80, &b_pi_minus_deltaEta_0_80);
    fChain->SetBranchAddress("pi_minus_deltaPhi_0.80", &pi_minus_deltaPhi_0_80, &b_pi_minus_deltaPhi_0_80);
    fChain->SetBranchAddress("pi_minus_pxasy_0.80", &pi_minus_pxasy_0_80, &b_pi_minus_pxasy_0_80);
    fChain->SetBranchAddress("pi_minus_pyasy_0.80", &pi_minus_pyasy_0_80, &b_pi_minus_pyasy_0_80);
    fChain->SetBranchAddress("pi_minus_pzasy_0.80", &pi_minus_pzasy_0_80, &b_pi_minus_pzasy_0_80);
    fChain->SetBranchAddress("pi_minus_pasy_0.80", &pi_minus_pasy_0_80, &b_pi_minus_pasy_0_80);
    fChain->SetBranchAddress("pi_minus_ptasy_0.80", &pi_minus_ptasy_0_80, &b_pi_minus_ptasy_0_80);
    fChain->SetBranchAddress("pi_minus_cpx_0.90", &pi_minus_cpx_0_90, &b_pi_minus_cpx_0_90);
    fChain->SetBranchAddress("pi_minus_cpy_0.90", &pi_minus_cpy_0_90, &b_pi_minus_cpy_0_90);
    fChain->SetBranchAddress("pi_minus_cpz_0.90", &pi_minus_cpz_0_90, &b_pi_minus_cpz_0_90);
    fChain->SetBranchAddress("pi_minus_cpt_0.90", &pi_minus_cpt_0_90, &b_pi_minus_cpt_0_90);
    fChain->SetBranchAddress("pi_minus_cp_0.90", &pi_minus_cp_0_90, &b_pi_minus_cp_0_90);
    fChain->SetBranchAddress("pi_minus_cmult_0.90", &pi_minus_cmult_0_90, &b_pi_minus_cmult_0_90);
    fChain->SetBranchAddress("pi_minus_deltaEta_0.90", &pi_minus_deltaEta_0_90, &b_pi_minus_deltaEta_0_90);
    fChain->SetBranchAddress("pi_minus_deltaPhi_0.90", &pi_minus_deltaPhi_0_90, &b_pi_minus_deltaPhi_0_90);
    fChain->SetBranchAddress("pi_minus_pxasy_0.90", &pi_minus_pxasy_0_90, &b_pi_minus_pxasy_0_90);
    fChain->SetBranchAddress("pi_minus_pyasy_0.90", &pi_minus_pyasy_0_90, &b_pi_minus_pyasy_0_90);
    fChain->SetBranchAddress("pi_minus_pzasy_0.90", &pi_minus_pzasy_0_90, &b_pi_minus_pzasy_0_90);
    fChain->SetBranchAddress("pi_minus_pasy_0.90", &pi_minus_pasy_0_90, &b_pi_minus_pasy_0_90);
    fChain->SetBranchAddress("pi_minus_ptasy_0.90", &pi_minus_ptasy_0_90, &b_pi_minus_ptasy_0_90);
    fChain->SetBranchAddress("pi_minus_cpx_1.00", &pi_minus_cpx_1_00, &b_pi_minus_cpx_1_00);
    fChain->SetBranchAddress("pi_minus_cpy_1.00", &pi_minus_cpy_1_00, &b_pi_minus_cpy_1_00);
    fChain->SetBranchAddress("pi_minus_cpz_1.00", &pi_minus_cpz_1_00, &b_pi_minus_cpz_1_00);
    fChain->SetBranchAddress("pi_minus_cpt_1.00", &pi_minus_cpt_1_00, &b_pi_minus_cpt_1_00);
    fChain->SetBranchAddress("pi_minus_cp_1.00", &pi_minus_cp_1_00, &b_pi_minus_cp_1_00);
    fChain->SetBranchAddress("pi_minus_cmult_1.00", &pi_minus_cmult_1_00, &b_pi_minus_cmult_1_00);
    fChain->SetBranchAddress("pi_minus_deltaEta_1.00", &pi_minus_deltaEta_1_00, &b_pi_minus_deltaEta_1_00);
    fChain->SetBranchAddress("pi_minus_deltaPhi_1.00", &pi_minus_deltaPhi_1_00, &b_pi_minus_deltaPhi_1_00);
    fChain->SetBranchAddress("pi_minus_pxasy_1.00", &pi_minus_pxasy_1_00, &b_pi_minus_pxasy_1_00);
    fChain->SetBranchAddress("pi_minus_pyasy_1.00", &pi_minus_pyasy_1_00, &b_pi_minus_pyasy_1_00);
    fChain->SetBranchAddress("pi_minus_pzasy_1.00", &pi_minus_pzasy_1_00, &b_pi_minus_pzasy_1_00);
    fChain->SetBranchAddress("pi_minus_pasy_1.00", &pi_minus_pasy_1_00, &b_pi_minus_pasy_1_00);
    fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
    fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
    fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
    fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
    fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
    fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
    fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
    fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
    fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
    fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
    fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
    fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
    fChain->SetBranchAddress("PVX", PVX, &b_PVX);
    fChain->SetBranchAddress("PVY", PVY, &b_PVY);
    fChain->SetBranchAddress("PVZ", PVZ, &b_PVZ);
    fChain->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
    fChain->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
    fChain->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
    fChain->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
    fChain->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
    fChain->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
    fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
    fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
    fChain->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
    fChain->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
    fChain->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
    fChain->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
    fChain->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
    fChain->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
    fChain->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
    fChain->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
    fChain->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
    fChain->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
    fChain->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
    fChain->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
    fChain->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
    fChain->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
    fChain->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
    fChain->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
    fChain->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
    fChain->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
    fChain->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
    
   Notify();
}

Bool_t DecayTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DecayTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef DecayTree_cxx
