#ifndef MiniDecayTree_h
#define MiniDecayTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "DecayTree.h"
#include "TLorentzVector.h"


class MiniDecayTree: public DecayTree {
public :

    TString _inFileName;
    Bool_t _charmLess;
    TString _usePIDvar;
    MiniDecayTree(Decay::Type decay, Year::Type year, Ds_finalState::Type finalState, DataType::Type dataType, TString polarity, TString inFileLoc = "/auto/data/dargent/BsDsKpipi/", TString outFileLoc = "/auto/data/dargent/BsDsKpipi/", Bool_t charmLess = false, TString _usePIDvar = "Corr", Bool_t bkg = false, Bool_t ltu = false, Bool_t ss = false );
    virtual ~MiniDecayTree();
    
    virtual void  Init();
    virtual TTree* GetInputTree();
    virtual void Loop();

    void set_inFileName( TString s){ _inFileName = s; }
    void set_outFileName( TString s){ _outFileName = s; }
    const TString get_inFileName(){ return _inFileName; }
    const TString get_outFileName(){ return _outFileName; }

    inline Ds_finalState::Type get_Ds_finalState();
    inline void set_LorentzVectors();
    
    inline Bool_t Preselection_Cuts();    
    inline Bool_t PID_Cuts();
    inline Bool_t Veto_Cuts();
    inline Bool_t PhaseSpace_Cuts();
    inline Bool_t PIDCalib_Cuts();
    inline Bool_t MC_Cuts();
    inline Bool_t LTU_Cuts();

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
    // changed mass hypothesis
    TLorentzVector Kminus_fromDs_asPiminus_MissID;
    TLorentzVector Kminus_fromDs_asProton_MissID;
    
    Double_t        Bs_MINIPNEXTBEST;
    Double_t        Bs_MINIPCHI2NEXTBEST;

    Float_t         Bs_B0DTF_D_splus_piplus_0_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_0_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_0_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_0_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_0_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_1_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_1_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_1_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_1_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_1_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_D_splus_piplus_PZ[100];   //[Bs_B0DTF_nPV]

    Float_t         Bs_BsDTF_D_splus_piplus_0_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_0_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_0_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_0_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_0_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_1_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_1_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_1_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_1_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_1_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_D_splus_piplus_PZ[100];   //[Bs_BsDTF_nPV]

    Float_t         Bs_DTF_D_splus_piplus_0_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_0_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_0_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_0_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_0_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_1_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_1_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_1_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_1_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_1_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_D_splus_piplus_PZ[100];   //[Bs_DTF_nPV]

    Float_t         Bs_PV_Dplus_piplus_0_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_0_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_0_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_0_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_0_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_1_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_1_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_1_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_1_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_1_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_Dplus_piplus_PZ[100];   //[Bs_PV_nPV]

    Double_t        pi_plus_fromDs_ETA;
    Double_t        pi_plus_fromDs_MC12TuneV2_ProbNNmu;
    Double_t        pi_plus_fromDs_MC12TuneV2_ProbNNpi;
    Double_t        pi_plus_fromDs_MC12TuneV2_ProbNNk;
    Double_t        pi_plus_fromDs_MC12TuneV2_ProbNNp;
    Double_t        pi_plus_fromDs_MC12TuneV2_ProbNNghost;
    Double_t        pi_plus_fromDs_MC12TuneV3_ProbNNmu;
    Double_t        pi_plus_fromDs_MC12TuneV3_ProbNNpi;
    Double_t        pi_plus_fromDs_MC12TuneV3_ProbNNk;
    Double_t        pi_plus_fromDs_MC12TuneV3_ProbNNp;
    Double_t        pi_plus_fromDs_MC12TuneV3_ProbNNghost;
    Double_t        pi_plus_fromDs_MC12TuneV4_ProbNNmu;
    Double_t        pi_plus_fromDs_MC12TuneV4_ProbNNpi;
    Double_t        pi_plus_fromDs_MC12TuneV4_ProbNNk;
    Double_t        pi_plus_fromDs_MC12TuneV4_ProbNNp;
    Double_t        pi_plus_fromDs_MC12TuneV4_ProbNNghost;
    Double_t        pi_plus_fromDs_MC15TuneV1_ProbNNmu;
    Double_t        pi_plus_fromDs_MC15TuneV1_ProbNNpi;
    Double_t        pi_plus_fromDs_MC15TuneV1_ProbNNk;
    Double_t        pi_plus_fromDs_MC15TuneV1_ProbNNp;
    Double_t        pi_plus_fromDs_MC15TuneV1_ProbNNghost;
    Double_t        pi_plus_fromDs_IP_OWNPV;
    Double_t        pi_plus_fromDs_IPCHI2_OWNPV;
    Double_t        pi_plus_fromDs_P;
    Double_t        pi_plus_fromDs_PT;
    Double_t        pi_plus_fromDs_PE;
    Double_t        pi_plus_fromDs_PX;
    Double_t        pi_plus_fromDs_PY;
    Double_t        pi_plus_fromDs_PZ;
    Int_t           pi_plus_fromDs_ID;
    Double_t        pi_plus_fromDs_PIDmu;
    Double_t        pi_plus_fromDs_PIDK;
    Double_t        pi_plus_fromDs_PIDp;
    Double_t        pi_plus_fromDs_ProbNNk;
    Double_t        pi_plus_fromDs_ProbNNp;
    Double_t        pi_plus_fromDs_ProbNNpi;
    Double_t        pi_plus_fromDs_ProbNNmu;
    Double_t        pi_plus_fromDs_ProbNNghost;
    Bool_t          pi_plus_fromDs_isMuon;
    Double_t        pi_plus_fromDs_TRACK_CHI2NDOF;
    Double_t        pi_plus_fromDs_TRACK_GhostProb;
    Double_t        pi_plus_fromDs_ptasy_1_00;

    Double_t        pi_minus2_fromDs_ETA;
    Double_t        pi_minus2_fromDs_MC12TuneV2_ProbNNmu;
    Double_t        pi_minus2_fromDs_MC12TuneV2_ProbNNpi;
    Double_t        pi_minus2_fromDs_MC12TuneV2_ProbNNk;
    Double_t        pi_minus2_fromDs_MC12TuneV2_ProbNNp;
    Double_t        pi_minus2_fromDs_MC12TuneV2_ProbNNghost;
    Double_t        pi_minus2_fromDs_MC12TuneV3_ProbNNmu;
    Double_t        pi_minus2_fromDs_MC12TuneV3_ProbNNpi;
    Double_t        pi_minus2_fromDs_MC12TuneV3_ProbNNk;
    Double_t        pi_minus2_fromDs_MC12TuneV3_ProbNNp;
    Double_t        pi_minus2_fromDs_MC12TuneV3_ProbNNghost;
    Double_t        pi_minus2_fromDs_MC12TuneV4_ProbNNmu;
    Double_t        pi_minus2_fromDs_MC12TuneV4_ProbNNpi;
    Double_t        pi_minus2_fromDs_MC12TuneV4_ProbNNk;
    Double_t        pi_minus2_fromDs_MC12TuneV4_ProbNNp;
    Double_t        pi_minus2_fromDs_MC12TuneV4_ProbNNghost;
    Double_t        pi_minus2_fromDs_MC15TuneV1_ProbNNmu;
    Double_t        pi_minus2_fromDs_MC15TuneV1_ProbNNpi;
    Double_t        pi_minus2_fromDs_MC15TuneV1_ProbNNk;
    Double_t        pi_minus2_fromDs_MC15TuneV1_ProbNNp;
    Double_t        pi_minus2_fromDs_MC15TuneV1_ProbNNghost;
    Double_t        pi_minus2_fromDs_IP_OWNPV;
    Double_t        pi_minus2_fromDs_IPCHI2_OWNPV;
    Double_t        pi_minus2_fromDs_P;
    Double_t        pi_minus2_fromDs_PT;
    Double_t        pi_minus2_fromDs_PE;
    Double_t        pi_minus2_fromDs_PX;
    Double_t        pi_minus2_fromDs_PY;
    Double_t        pi_minus2_fromDs_PZ;
    Int_t           pi_minus2_fromDs_ID;
    Double_t        pi_minus2_fromDs_PIDmu;
    Double_t        pi_minus2_fromDs_PIDK;
    Double_t        pi_minus2_fromDs_PIDp;
    Double_t        pi_minus2_fromDs_ProbNNk;
    Double_t        pi_minus2_fromDs_ProbNNp;
    Double_t        pi_minus2_fromDs_ProbNNpi;
    Double_t        pi_minus2_fromDs_ProbNNmu;
    Double_t        pi_minus2_fromDs_ProbNNghost;
    Bool_t          pi_minus2_fromDs_isMuon;
    Double_t        pi_minus2_fromDs_TRACK_CHI2NDOF;
    Double_t        pi_minus2_fromDs_TRACK_GhostProb;
    Double_t        pi_minus2_fromDs_ptasy_1_00;

    Float_t         Bs_B0DTF_a_1_1260_plus_M[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_MERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_P[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_PERR[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_ctau[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_ctauErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_decayLength[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_decayLengthErr[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_0_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_0_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_0_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_0_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_0_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_1_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_1_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_1_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_1_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_1_PZ[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_ID[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_PE[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_PX[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_PY[100];   //[Bs_B0DTF_nPV]
    Float_t         Bs_B0DTF_a_1_1260_plus_piplus_PZ[100];   //[Bs_B0DTF_nPV]

    Float_t         Bs_BsDTF_a_1_1260_plus_M[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_MERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_P[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_PERR[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_ctau[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_ctauErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_decayLength[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_decayLengthErr[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_0_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_0_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_0_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_0_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_0_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_1_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_1_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_1_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_1_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_1_PZ[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_ID[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_PE[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_PX[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_PY[100];   //[Bs_BsDTF_nPV]
    Float_t         Bs_BsDTF_a_1_1260_plus_piplus_PZ[100];   //[Bs_BsDTF_nPV]

    Float_t         Bs_DTF_a_1_1260_plus_M[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_MERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_P[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_PERR[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_ctau[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_ctauErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_decayLength[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_decayLengthErr[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_0_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_0_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_0_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_0_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_0_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_1_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_1_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_1_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_1_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_1_PZ[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_ID[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_PE[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_PX[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_PY[100];   //[Bs_DTF_nPV]
    Float_t         Bs_DTF_a_1_1260_plus_piplus_PZ[100];   //[Bs_DTF_nPV]

    Float_t         Bs_PV_a_1_1260_plus_M[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_MERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_P[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_PERR[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_ctau[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_ctauErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_decayLength[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_decayLengthErr[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_0_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_0_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_0_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_0_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_0_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_1_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_1_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_1_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_1_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_1_PZ[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_ID[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_PE[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_PX[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_PY[100];   //[Bs_PV_nPV]
    Float_t         Bs_PV_a_1_1260_plus_piplus_PZ[100];   //[Bs_PV_nPV]

    Double_t        a_1_1260_plus_DOCA1;
    Double_t        a_1_1260_plus_DOCA2;
    Double_t        a_1_1260_plus_DOCA3;
    Double_t        a_1_1260_plus_ETA;
    Double_t        a_1_1260_plus_ENDVERTEX_X;
    Double_t        a_1_1260_plus_ENDVERTEX_Y;
    Double_t        a_1_1260_plus_ENDVERTEX_Z;
    Double_t        a_1_1260_plus_ENDVERTEX_XERR;
    Double_t        a_1_1260_plus_ENDVERTEX_YERR;
    Double_t        a_1_1260_plus_ENDVERTEX_ZERR;
    Double_t        a_1_1260_plus_ENDVERTEX_CHI2;
    Int_t           a_1_1260_plus_ENDVERTEX_NDOF;
    Double_t        a_1_1260_plus_OWNPV_X;
    Double_t        a_1_1260_plus_OWNPV_Y;
    Double_t        a_1_1260_plus_OWNPV_Z;
    Double_t        a_1_1260_plus_OWNPV_XERR;
    Double_t        a_1_1260_plus_OWNPV_YERR;
    Double_t        a_1_1260_plus_OWNPV_ZERR;
    Double_t        a_1_1260_plus_OWNPV_CHI2;
    Int_t           a_1_1260_plus_OWNPV_NDOF;
    Double_t        a_1_1260_plus_IP_OWNPV;
    Double_t        a_1_1260_plus_IPCHI2_OWNPV;
    Double_t        a_1_1260_plus_FD_OWNPV;
    Double_t        a_1_1260_plus_FDCHI2_OWNPV;
    Double_t        a_1_1260_plus_DIRA_OWNPV;
    Double_t        a_1_1260_plus_ORIVX_X;
    Double_t        a_1_1260_plus_ORIVX_Y;
    Double_t        a_1_1260_plus_ORIVX_Z;
    Double_t        a_1_1260_plus_ORIVX_XERR;
    Double_t        a_1_1260_plus_ORIVX_YERR;
    Double_t        a_1_1260_plus_ORIVX_ZERR;
    Double_t        a_1_1260_plus_ORIVX_CHI2;
    Int_t           a_1_1260_plus_ORIVX_NDOF;
    Double_t        a_1_1260_plus_FD_ORIVX;
    Double_t        a_1_1260_plus_FDCHI2_ORIVX;
    Double_t        a_1_1260_plus_DIRA_ORIVX;
    Double_t        a_1_1260_plus_P;
    Double_t        a_1_1260_plus_PT;
    Double_t        a_1_1260_plus_PE;
    Double_t        a_1_1260_plus_PX;
    Double_t        a_1_1260_plus_PY;
    Double_t        a_1_1260_plus_PZ;
    Double_t        a_1_1260_plus_MM;
    Double_t        a_1_1260_plus_MMERR;
    Int_t           a_1_1260_plus_ID;
    Double_t        a_1_1260_plus_TAU;
    Double_t        a_1_1260_plus_TAUERR;
    Double_t        a_1_1260_plus_TAUCHI2;
    Double_t        a_1_1260_plus_ptasy_1_00;
    
    Double_t        pi_plus1_ETA;
    Double_t        pi_plus1_MC12TuneV2_ProbNNmu;
    Double_t        pi_plus1_MC12TuneV2_ProbNNpi;
    Double_t        pi_plus1_MC12TuneV2_ProbNNk;
    Double_t        pi_plus1_MC12TuneV2_ProbNNp;
    Double_t        pi_plus1_MC12TuneV2_ProbNNghost;
    Double_t        pi_plus1_MC12TuneV3_ProbNNmu;
    Double_t        pi_plus1_MC12TuneV3_ProbNNpi;
    Double_t        pi_plus1_MC12TuneV3_ProbNNk;
    Double_t        pi_plus1_MC12TuneV3_ProbNNp;
    Double_t        pi_plus1_MC12TuneV3_ProbNNghost;
    Double_t        pi_plus1_IP_OWNPV;
    Double_t        pi_plus1_IPCHI2_OWNPV;
    Double_t        pi_plus1_P;
    Double_t        pi_plus1_PT;
    Double_t        pi_plus1_PE;
    Double_t        pi_plus1_PX;
    Double_t        pi_plus1_PY;
    Double_t        pi_plus1_PZ;
    Int_t           pi_plus1_ID;
    Double_t        pi_plus1_PIDmu;
    Double_t        pi_plus1_PIDK;
    Double_t        pi_plus1_PIDp;
    Double_t        pi_plus1_ProbNNk;
    Double_t        pi_plus1_ProbNNp;
    Double_t        pi_plus1_ProbNNpi;
    Double_t        pi_plus1_ProbNNmu;
    Double_t        pi_plus1_ProbNNghost;
    Bool_t          pi_plus1_isMuon;
    Double_t        pi_plus1_TRACK_CHI2NDOF;
    Double_t        pi_plus1_TRACK_GhostProb;
    Double_t        pi_plus1_ptasy_1_00;
    Double_t        pi_plus2_ETA;
    Double_t        pi_plus2_MC12TuneV2_ProbNNmu;
    Double_t        pi_plus2_MC12TuneV2_ProbNNpi;
    Double_t        pi_plus2_MC12TuneV2_ProbNNk;
    Double_t        pi_plus2_MC12TuneV2_ProbNNp;
    Double_t        pi_plus2_MC12TuneV2_ProbNNghost;
    Double_t        pi_plus2_MC12TuneV3_ProbNNmu;
    Double_t        pi_plus2_MC12TuneV3_ProbNNpi;
    Double_t        pi_plus2_MC12TuneV3_ProbNNk;
    Double_t        pi_plus2_MC12TuneV3_ProbNNp;
    Double_t        pi_plus2_MC12TuneV3_ProbNNghost;
    Double_t        pi_plus2_IP_OWNPV;
    Double_t        pi_plus2_IPCHI2_OWNPV;
    Double_t        pi_plus2_P;
    Double_t        pi_plus2_PT;
    Double_t        pi_plus2_PE;
    Double_t        pi_plus2_PX;
    Double_t        pi_plus2_PY;
    Double_t        pi_plus2_PZ;
    Int_t           pi_plus2_ID;
    Double_t        pi_plus2_PIDmu;
    Double_t        pi_plus2_PIDK;
    Double_t        pi_plus2_PIDp;
    Double_t        pi_plus2_ProbNNk;
    Double_t        pi_plus2_ProbNNp;
    Double_t        pi_plus2_ProbNNpi;
    Double_t        pi_plus2_ProbNNmu;
    Double_t        pi_plus2_ProbNNghost;
    Bool_t          pi_plus2_isMuon;
    Double_t        pi_plus2_TRACK_CHI2NDOF;
    Double_t        pi_plus2_TRACK_GhostProb;
    Double_t        pi_plus2_ptasy_1_00;

    Int_t	   Bs_TRUEID,Ds_TRUEID;
    Int_t	   Bs_BKGCAT;
    Int_t	   K_plus_TRUEID;
    Int_t	   pi_plus_TRUEID;
    Int_t	   pi_minus_TRUEID;
    Int_t	   pi_plus1_TRUEID;
    Int_t	   pi_plus2_TRUEID;
    Int_t	   K_plus_fromDs_TRUEID;
    Int_t	   K_minus_fromDs_TRUEID;
    Int_t	   pi_minus_fromDs_TRUEID;
    Int_t	   pi_plus_fromDs_TRUEID;
    Int_t	   pi_minus2_fromDs_TRUEID;

    Int_t	   Ds_MC_MOTHER_ID;
    Int_t	   K_plus_MC_MOTHER_ID;
    Int_t	   pi_plus_MC_MOTHER_ID;
    Int_t	   pi_minus_MC_MOTHER_ID;
    Int_t	   pi_plus1_MC_MOTHER_ID;
    Int_t	   pi_plus2_MC_MOTHER_ID;
    Int_t	   K_plus_fromDs_MC_MOTHER_ID;
    Int_t	   K_minus_fromDs_MC_MOTHER_ID;
    Int_t	   pi_minus_fromDs_MC_MOTHER_ID;
    Int_t	   pi_plus_fromDs_MC_MOTHER_ID;
    Int_t	   pi_minus2_fromDs_MC_MOTHER_ID;

    Double_t K_plus_PIDK_gen_MagDown,pi_plus_PIDK_gen_MagDown,pi_minus_PIDK_gen_MagDown,K_plus_fromDs_PIDK_gen_MagDown,K_minus_fromDs_PIDK_gen_MagDown,
pi_minus_fromDs_PIDK_gen_MagDown,pi_minus2_fromDs_PIDK_gen_MagDown,pi_plus_fromDs_PIDK_gen_MagDown;
    Double_t K_plus_PIDK_corr_MagDown,pi_plus_PIDK_corr_MagDown,pi_minus_PIDK_corr_MagDown,K_plus_fromDs_PIDK_corr_MagDown,K_minus_fromDs_PIDK_corr_MagDown,
pi_minus_fromDs_PIDK_corr_MagDown,pi_minus2_fromDs_PIDK_corr_MagDown,pi_plus_fromDs_PIDK_corr_MagDown;

    Double_t K_plus_PIDK_gen_MagUp,pi_plus_PIDK_gen_MagUp,pi_minus_PIDK_gen_MagUp,K_plus_fromDs_PIDK_gen_MagUp,K_minus_fromDs_PIDK_gen_MagUp,
pi_minus_fromDs_PIDK_gen_MagUp,pi_minus2_fromDs_PIDK_gen_MagUp,pi_plus_fromDs_PIDK_gen_MagUp;
    Double_t K_plus_PIDK_corr_MagUp,pi_plus_PIDK_corr_MagUp,pi_minus_PIDK_corr_MagUp,K_plus_fromDs_PIDK_corr_MagUp,K_minus_fromDs_PIDK_corr_MagUp,
pi_minus_fromDs_PIDK_corr_MagUp,pi_minus2_fromDs_PIDK_corr_MagUp,pi_plus_fromDs_PIDK_corr_MagUp;

   Double_t pi_plus1_PIDK_gen_MagDown, pi_plus1_PIDK_gen_MagUp, pi_plus1_PIDK_corr_MagDown, pi_plus1_PIDK_corr_MagUp;
   Double_t pi_plus2_PIDK_gen_MagDown, pi_plus2_PIDK_gen_MagUp, pi_plus2_PIDK_corr_MagDown, pi_plus2_PIDK_corr_MagUp;

    Double_t K_plus_PIDp_gen_MagDown,pi_plus_PIDp_gen_MagDown,pi_minus_PIDp_gen_MagDown,K_plus_fromDs_PIDp_gen_MagDown,K_minus_fromDs_PIDp_gen_MagDown,
pi_minus_fromDs_PIDp_gen_MagDown;
    Double_t K_plus_PIDp_corr_MagDown,pi_plus_PIDp_corr_MagDown,pi_minus_PIDp_corr_MagDown,K_plus_fromDs_PIDp_corr_MagDown,K_minus_fromDs_PIDp_corr_MagDown,
pi_minus_fromDs_PIDp_corr_MagDown;

    Double_t K_plus_PIDp_gen_MagUp,pi_plus_PIDp_gen_MagUp,pi_minus_PIDp_gen_MagUp,K_plus_fromDs_PIDp_gen_MagUp,K_minus_fromDs_PIDp_gen_MagUp,
pi_minus_fromDs_PIDp_gen_MagUp;
    Double_t K_plus_PIDp_corr_MagUp,pi_plus_PIDp_corr_MagUp,pi_minus_PIDp_corr_MagUp,K_plus_fromDs_PIDp_corr_MagUp,K_minus_fromDs_PIDp_corr_MagUp,
pi_minus_fromDs_PIDp_corr_MagUp;

   Double_t pi_plus1_PIDp_gen_MagDown, pi_plus1_PIDp_gen_MagUp, pi_plus1_PIDp_corr_MagDown, pi_plus1_PIDp_corr_MagUp;
   Double_t pi_plus2_PIDp_gen_MagDown, pi_plus2_PIDp_gen_MagUp, pi_plus2_PIDp_corr_MagDown, pi_plus2_PIDp_corr_MagUp;

/*
 Short_t         Bs_OS_Muon_DEC;
 Float_t         Bs_OS_Muon_PROB;
 Int_t           Bs_OS_Muon_PARTICLES_NUM;
 Float_t         Bs_OS_Muon_PARTICLES_ID[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_P[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PX[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PY[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PZ[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PT[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_THETA[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PIDmu[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PIDk[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_PIDp[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_IP_OWNPV[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_IP_BVertex[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Float_t         Bs_OS_Muon_PARTICLES_IPCHI2_BVertex[20];   //[Bs_OS_Muon_PARTICLES_NUM]
 Short_t         Bs_OS_Electron_DEC;
 Float_t         Bs_OS_Electron_PROB;
 Int_t           Bs_OS_Electron_PARTICLES_NUM;
 Float_t         Bs_OS_Electron_PARTICLES_ID[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_P[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PX[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PY[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PZ[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PT[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_THETA[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PIDmu[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PIDk[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_PIDp[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_IP_OWNPV[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_IP_BVertex[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Float_t         Bs_OS_Electron_PARTICLES_IPCHI2_BVertex[20];   //[Bs_OS_Electron_PARTICLES_NUM]
 Short_t         Bs_OS_Kaon_DEC;
 Float_t         Bs_OS_Kaon_PROB;
 Int_t           Bs_OS_Kaon_PARTICLES_NUM;
 Float_t         Bs_OS_Kaon_PARTICLES_ID[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_P[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PX[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PY[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PZ[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PT[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_THETA[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PIDmu[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PIDk[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_PIDp[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_IP_OWNPV[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_IP_BVertex[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Float_t         Bs_OS_Kaon_PARTICLES_IPCHI2_BVertex[20];   //[Bs_OS_Kaon_PARTICLES_NUM]
 Short_t         Bs_SS_Kaon_DEC;
 Float_t         Bs_SS_Kaon_PROB;
 Int_t           Bs_SS_Kaon_PARTICLES_NUM;
 Float_t         Bs_SS_Kaon_PARTICLES_ID[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_P[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PX[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PY[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PZ[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PT[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_THETA[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PIDmu[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PIDk[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_PIDp[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_IP_OWNPV[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_IP_BVertex[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Float_t         Bs_SS_Kaon_PARTICLES_IPCHI2_BVertex[20];   //[Bs_SS_Kaon_PARTICLES_NUM]
 Short_t         Bs_SS_Pion_DEC;
 Float_t         Bs_SS_Pion_PROB;
 Int_t           Bs_SS_Pion_PARTICLES_NUM;
 Float_t         Bs_SS_Pion_PARTICLES_ID[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_P[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PX[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PY[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PZ[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PT[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_THETA[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PIDmu[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PIDk[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_PIDp[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_IP_OWNPV[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_IP_BVertex[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Float_t         Bs_SS_Pion_PARTICLES_IPCHI2_BVertex[20];   //[Bs_SS_Pion_PARTICLES_NUM]
 Short_t         Bs_SS_PionBDT_DEC;
 Float_t         Bs_SS_PionBDT_PROB;
 Int_t           Bs_SS_PionBDT_PARTICLES_NUM;
 Float_t         Bs_SS_PionBDT_PARTICLES_ID[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_P[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PX[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PY[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PZ[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PT[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_THETA[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PIDmu[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PIDk[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_PIDp[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_IP_OWNPV[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_IP_BVertex[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Float_t         Bs_SS_PionBDT_PARTICLES_IPCHI2_BVertex[20];   //[Bs_SS_PionBDT_PARTICLES_NUM]
 Short_t         Bs_VtxCharge_DEC;
 Float_t         Bs_VtxCharge_PROB;
 Int_t           Bs_VtxCharge_PARTICLES_NUM;
 Float_t         Bs_VtxCharge_PARTICLES_ID[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_P[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PX[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PY[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PZ[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PT[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_THETA[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PIDmu[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PIDk[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_PIDp[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_IP_OWNPV[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_IP_BVertex[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Float_t         Bs_VtxCharge_PARTICLES_IPCHI2_BVertex[20];   //[Bs_VtxCharge_PARTICLES_NUM]
 Short_t         Bs_OS_nnetKaon_DEC;
 Float_t         Bs_OS_nnetKaon_PROB;
 Int_t           Bs_OS_nnetKaon_PARTICLES_NUM;
 Float_t         Bs_OS_nnetKaon_PARTICLES_ID[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_P[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PX[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PY[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PZ[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PT[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_THETA[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PIDmu[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PIDk[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_PIDp[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_IP_OWNPV[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_IP_BVertex[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_OS_nnetKaon_PARTICLES_IPCHI2_BVertex[20];   //[Bs_OS_nnetKaon_PARTICLES_NUM]
 Short_t         Bs_SS_nnetKaon_DEC;
 Float_t         Bs_SS_nnetKaon_PROB;
 Int_t           Bs_SS_nnetKaon_PARTICLES_NUM;
 Float_t         Bs_SS_nnetKaon_PARTICLES_ID[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_P[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PX[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PY[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PZ[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PT[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_THETA[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PIDmu[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PIDk[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_PIDp[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_IP_OWNPV[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_IP_BVertex[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Float_t         Bs_SS_nnetKaon_PARTICLES_IPCHI2_BVertex[20];   //[Bs_SS_nnetKaon_PARTICLES_NUM]
 Short_t         Bs_SS_Proton_DEC;
 Float_t         Bs_SS_Proton_PROB;
 Int_t           Bs_SS_Proton_PARTICLES_NUM;
 Float_t         Bs_SS_Proton_PARTICLES_ID[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_P[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PX[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PY[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PZ[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PT[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_THETA[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PIDmu[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PIDk[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_PIDp[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_IP_OWNPV[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_IP_BVertex[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Float_t         Bs_SS_Proton_PARTICLES_IPCHI2_BVertex[20];   //[Bs_SS_Proton_PARTICLES_NUM]
 Short_t         Bs_OS_Charm_DEC;
 Float_t         Bs_OS_Charm_PROB;
 Int_t           Bs_OS_Charm_PARTICLES_NUM;
 Float_t         Bs_OS_Charm_PARTICLES_ID[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_P[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PX[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PY[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PZ[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PT[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_THETA[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PIDmu[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PIDk[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_PIDp[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_IP_OWNPV[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_IPCHI2_OWNPV[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_IP_BVertex[20];   //[Bs_OS_Charm_PARTICLES_NUM]
 Float_t         Bs_OS_Charm_PARTICLES_IPCHI2_BVertex[20];   //[Bs_OS_Charm_PARTICLES_NUM]
    Int_t           Bs_BsTaggingTool_TAGDECISION_OS;
    Double_t        Bs_BsTaggingTool_TAGOMEGA_OS;
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
   */
   
   
    TBranch        *b_Bs_MINIPNEXTBEST;   //!
    TBranch        *b_Bs_MINIPCHI2NEXTBEST;   //!

    TBranch         *b_Bs_B0DTF_D_splus_piplus_0_ID;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_0_PE;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_0_PX;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_0_PY;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_0_PZ;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_1_ID;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_1_PE;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_1_PX;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_1_PY;  //[Bs_B0DTF_nPV]
    TBranch         *b_Bs_B0DTF_D_splus_piplus_1_PZ;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_D_splus_piplus_ID;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PE;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PX;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PY;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_D_splus_piplus_PZ;  //[Bs_B0DTF_nPV]
    
    TBranch        *b_Bs_BsDTF_D_splus_piplus_0_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_0_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_0_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_0_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_0_PZ;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_1_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_1_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_1_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_1_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_1_PZ;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_D_splus_piplus_PZ;  //[Bs_BsDTF_nPV]
    
    TBranch        *b_Bs_DTF_D_splus_piplus_0_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_0_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_0_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_0_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_0_PZ;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_1_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_1_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_1_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_1_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_1_PZ;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_D_splus_piplus_PZ;  //[Bs_DTF_nPV]
    
    TBranch        *b_Bs_PV_Dplus_piplus_0_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_0_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_0_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_0_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_0_PZ;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_1_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_1_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_1_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_1_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_1_PZ;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_Dplus_piplus_PZ;  //[Bs_PV_nPV]
    
    TBranch        *b_pi_plus_fromDs_ETA;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV2_ProbNNmu;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV2_ProbNNpi;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV2_ProbNNk;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV2_ProbNNp;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV2_ProbNNghost;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV3_ProbNNmu;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV3_ProbNNpi;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV3_ProbNNk;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV3_ProbNNp;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV3_ProbNNghost;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV4_ProbNNmu;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV4_ProbNNpi;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV4_ProbNNk;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV4_ProbNNp;   //!
    TBranch        *b_pi_plus_fromDs_MC12TuneV4_ProbNNghost;   //!
    TBranch        *b_pi_plus_fromDs_MC15TuneV1_ProbNNmu;   //!
    TBranch        *b_pi_plus_fromDs_MC15TuneV1_ProbNNpi;   //!
    TBranch        *b_pi_plus_fromDs_MC15TuneV1_ProbNNk;   //!
    TBranch        *b_pi_plus_fromDs_MC15TuneV1_ProbNNp;   //!
    TBranch        *b_pi_plus_fromDs_MC15TuneV1_ProbNNghost;   //!
    TBranch        *b_pi_plus_fromDs_IP_OWNPV;   //!
    TBranch        *b_pi_plus_fromDs_IPCHI2_OWNPV;   //!
    TBranch        *b_pi_plus_fromDs_P;   //!
    TBranch        *b_pi_plus_fromDs_PT;   //!
    TBranch        *b_pi_plus_fromDs_PE;   //!
    TBranch        *b_pi_plus_fromDs_PX;   //!
    TBranch        *b_pi_plus_fromDs_PY;   //!
    TBranch        *b_pi_plus_fromDs_PZ;   //!
    TBranch        *b_pi_plus_fromDs_ID;   //!
    TBranch        *b_pi_plus_fromDs_PIDmu;   //!
    TBranch        *b_pi_plus_fromDs_PIDK;   //!
    TBranch        *b_pi_plus_fromDs_PIDp;   //!
    TBranch        *b_pi_plus_fromDs_ProbNNk;   //!
    TBranch        *b_pi_plus_fromDs_ProbNNp;   //!
    TBranch        *b_pi_plus_fromDs_ProbNNpi;   //!
    TBranch        *b_pi_plus_fromDs_ProbNNmu;   //!
    TBranch        *b_pi_plus_fromDs_ProbNNghost;   //!
    TBranch        *b_pi_plus_fromDs_isMuon;   //!
    TBranch        *b_pi_plus_fromDs_TRACK_CHI2NDOF;   //!
    TBranch        *b_pi_plus_fromDs_TRACK_GhostProb;   //!
    TBranch        *b_pi_plus_fromDs_ptasy_1_00;   //!
    
    TBranch      *b_pi_minus2_fromDs_ETA;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV2_ProbNNmu;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV2_ProbNNpi;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV2_ProbNNk;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV2_ProbNNp;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV2_ProbNNghost;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV3_ProbNNmu;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV3_ProbNNpi;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV3_ProbNNk;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV3_ProbNNp;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV3_ProbNNghost;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV4_ProbNNmu;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV4_ProbNNpi;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV4_ProbNNk;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV4_ProbNNp;
    TBranch      *b_pi_minus2_fromDs_MC12TuneV4_ProbNNghost;
    TBranch      *b_pi_minus2_fromDs_MC15TuneV1_ProbNNmu;
    TBranch      *b_pi_minus2_fromDs_MC15TuneV1_ProbNNpi;
    TBranch      *b_pi_minus2_fromDs_MC15TuneV1_ProbNNk;
    TBranch      *b_pi_minus2_fromDs_MC15TuneV1_ProbNNp;
    TBranch      *b_pi_minus2_fromDs_MC15TuneV1_ProbNNghost;
    TBranch      *b_pi_minus2_fromDs_IP_OWNPV;
    TBranch      *b_pi_minus2_fromDs_IPCHI2_OWNPV;
    TBranch      *b_pi_minus2_fromDs_P;
    TBranch      *b_pi_minus2_fromDs_PT;
    TBranch      *b_pi_minus2_fromDs_PE;
    TBranch      *b_pi_minus2_fromDs_PX;
    TBranch      *b_pi_minus2_fromDs_PY;
    TBranch      *b_pi_minus2_fromDs_PZ;
    TBranch         *b_pi_minus2_fromDs_ID;
    TBranch      *b_pi_minus2_fromDs_PIDmu;
    TBranch      *b_pi_minus2_fromDs_PIDK;
    TBranch      *b_pi_minus2_fromDs_PIDp;
    TBranch      *b_pi_minus2_fromDs_ProbNNk;
    TBranch      *b_pi_minus2_fromDs_ProbNNp;
    TBranch      *b_pi_minus2_fromDs_ProbNNpi;
    TBranch      *b_pi_minus2_fromDs_ProbNNmu;
    TBranch      *b_pi_minus2_fromDs_ProbNNghost;
    TBranch        *b_pi_minus2_fromDs_isMuon;
    TBranch      *b_pi_minus2_fromDs_TRACK_CHI2NDOF;
    TBranch      *b_pi_minus2_fromDs_TRACK_GhostProb;
    TBranch      *b_pi_minus2_fromDs_ptasy_1_00;
    
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_M;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_MERR;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_P;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_PERR;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_ctau;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_ctauErr;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_decayLength;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_decayLengthErr;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_0_ID;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_0_PE;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_0_PX;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_0_PY;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_0_PZ;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_1_ID;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_1_PE;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_1_PX;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_1_PY;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_1_PZ;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_ID;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_PE;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_PX;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_PY;  //[Bs_B0DTF_nPV]
    TBranch        *b_Bs_B0DTF_a_1_1260_plus_piplus_PZ;  //[Bs_B0DTF_nPV]
    
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_M;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_MERR;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_P;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_PERR;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_ctau;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_ctauErr;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_decayLength;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_decayLengthErr;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_0_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_0_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_0_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_0_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_0_PZ;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_1_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_1_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_1_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_1_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_1_PZ;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_ID;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_PE;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_PX;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_PY;  //[Bs_BsDTF_nPV]
    TBranch        *b_Bs_BsDTF_a_1_1260_plus_piplus_PZ;  //[Bs_BsDTF_nPV]
    
    TBranch        *b_Bs_DTF_a_1_1260_plus_M;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_MERR;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_P;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_PERR;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_ctau;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_ctauErr;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_decayLength;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_decayLengthErr;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_0_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_0_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_0_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_0_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_0_PZ;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_1_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_1_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_1_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_1_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_1_PZ;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_ID;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_PE;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_PX;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_PY;  //[Bs_DTF_nPV]
    TBranch        *b_Bs_DTF_a_1_1260_plus_piplus_PZ;  //[Bs_DTF_nPV]
    
    TBranch        *b_Bs_PV_a_1_1260_plus_M;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_MERR;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_P;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_PERR;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_ctau;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_ctauErr;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_decayLength;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_decayLengthErr;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_0_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_0_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_0_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_0_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_0_PZ;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_1_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_1_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_1_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_1_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_1_PZ;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_ID;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_PE;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_PX;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_PY;  //[Bs_PV_nPV]
    TBranch        *b_Bs_PV_a_1_1260_plus_piplus_PZ;  //[Bs_PV_nPV]
    
    TBranch        *b_a_1_1260_plus_DOCA1;
    TBranch        *b_a_1_1260_plus_DOCA2;
    TBranch        *b_a_1_1260_plus_DOCA3;
    TBranch        *b_a_1_1260_plus_ETA;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_X;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_Y;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_Z;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_XERR;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_YERR;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_ZERR;
    TBranch        *b_a_1_1260_plus_ENDVERTEX_CHI2;
    TBranch           *b_a_1_1260_plus_ENDVERTEX_NDOF;
    TBranch        *b_a_1_1260_plus_OWNPV_X;
    TBranch        *b_a_1_1260_plus_OWNPV_Y;
    TBranch        *b_a_1_1260_plus_OWNPV_Z;
    TBranch        *b_a_1_1260_plus_OWNPV_XERR;
    TBranch        *b_a_1_1260_plus_OWNPV_YERR;
    TBranch        *b_a_1_1260_plus_OWNPV_ZERR;
    TBranch        *b_a_1_1260_plus_OWNPV_CHI2;
    TBranch           *b_a_1_1260_plus_OWNPV_NDOF;
    TBranch        *b_a_1_1260_plus_IP_OWNPV;
    TBranch        *b_a_1_1260_plus_IPCHI2_OWNPV;
    TBranch        *b_a_1_1260_plus_FD_OWNPV;
    TBranch        *b_a_1_1260_plus_FDCHI2_OWNPV;
    TBranch        *b_a_1_1260_plus_DIRA_OWNPV;
    TBranch        *b_a_1_1260_plus_ORIVX_X;
    TBranch        *b_a_1_1260_plus_ORIVX_Y;
    TBranch        *b_a_1_1260_plus_ORIVX_Z;
    TBranch        *b_a_1_1260_plus_ORIVX_XERR;
    TBranch        *b_a_1_1260_plus_ORIVX_YERR;
    TBranch        *b_a_1_1260_plus_ORIVX_ZERR;
    TBranch        *b_a_1_1260_plus_ORIVX_CHI2;
    TBranch           *b_a_1_1260_plus_ORIVX_NDOF;
    TBranch        *b_a_1_1260_plus_FD_ORIVX;
    TBranch        *b_a_1_1260_plus_FDCHI2_ORIVX;
    TBranch        *b_a_1_1260_plus_DIRA_ORIVX;
    TBranch        *b_a_1_1260_plus_P;
    TBranch        *b_a_1_1260_plus_PT;
    TBranch        *b_a_1_1260_plus_PE;
    TBranch        *b_a_1_1260_plus_PX;
    TBranch        *b_a_1_1260_plus_PY;
    TBranch        *b_a_1_1260_plus_PZ;
    TBranch        *b_a_1_1260_plus_MM;
    TBranch        *b_a_1_1260_plus_MMERR;
    TBranch           *b_a_1_1260_plus_ID;
    TBranch        *b_a_1_1260_plus_TAU;
    TBranch        *b_a_1_1260_plus_TAUERR;
    TBranch        *b_a_1_1260_plus_TAUCHI2;
    TBranch        *b_a_1_1260_plus_ptasy_1_00;
    
    TBranch      *b_pi_plus1_ETA;
    TBranch      *b_pi_plus1_MC12TuneV2_ProbNNmu;
    TBranch      *b_pi_plus1_MC12TuneV2_ProbNNpi;
    TBranch      *b_pi_plus1_MC12TuneV2_ProbNNk;
    TBranch      *b_pi_plus1_MC12TuneV2_ProbNNp;
    TBranch      *b_pi_plus1_MC12TuneV2_ProbNNghost;
    TBranch      *b_pi_plus1_MC12TuneV3_ProbNNmu;
    TBranch      *b_pi_plus1_MC12TuneV3_ProbNNpi;
    TBranch      *b_pi_plus1_MC12TuneV3_ProbNNk;
    TBranch      *b_pi_plus1_MC12TuneV3_ProbNNp;
    TBranch      *b_pi_plus1_MC12TuneV3_ProbNNghost;
    TBranch      *b_pi_plus1_IP_OWNPV;
    TBranch      *b_pi_plus1_IPCHI2_OWNPV;
    TBranch      *b_pi_plus1_P;
    TBranch      *b_pi_plus1_PT;
    TBranch      *b_pi_plus1_PE;
    TBranch      *b_pi_plus1_PX;
    TBranch      *b_pi_plus1_PY;
    TBranch      *b_pi_plus1_PZ;
    TBranch         *b_pi_plus1_ID;
    TBranch      *b_pi_plus1_PIDmu;
    TBranch      *b_pi_plus1_PIDK;
    TBranch      *b_pi_plus1_PIDp;
    TBranch      *b_pi_plus1_ProbNNk;
    TBranch      *b_pi_plus1_ProbNNp;
    TBranch      *b_pi_plus1_ProbNNpi;
    TBranch      *b_pi_plus1_ProbNNmu;
    TBranch      *b_pi_plus1_ProbNNghost;
    TBranch        *b_pi_plus1_isMuon;
    TBranch      *b_pi_plus1_TRACK_CHI2NDOF;
    TBranch      *b_pi_plus1_TRACK_GhostProb;
    TBranch      *b_pi_plus1_ptasy_1_00;
    TBranch      *b_pi_plus2_ETA;
    TBranch      *b_pi_plus2_MC12TuneV2_ProbNNmu;
    TBranch      *b_pi_plus2_MC12TuneV2_ProbNNpi;
    TBranch      *b_pi_plus2_MC12TuneV2_ProbNNk;
    TBranch      *b_pi_plus2_MC12TuneV2_ProbNNp;
    TBranch      *b_pi_plus2_MC12TuneV2_ProbNNghost;
    TBranch      *b_pi_plus2_MC12TuneV3_ProbNNmu;
    TBranch      *b_pi_plus2_MC12TuneV3_ProbNNpi;
    TBranch      *b_pi_plus2_MC12TuneV3_ProbNNk;
    TBranch      *b_pi_plus2_MC12TuneV3_ProbNNp;
    TBranch      *b_pi_plus2_MC12TuneV3_ProbNNghost;
    TBranch      *b_pi_plus2_IP_OWNPV;
    TBranch      *b_pi_plus2_IPCHI2_OWNPV;
    TBranch      *b_pi_plus2_P;
    TBranch      *b_pi_plus2_PT;
    TBranch      *b_pi_plus2_PE;
    TBranch      *b_pi_plus2_PX;
    TBranch      *b_pi_plus2_PY;
    TBranch      *b_pi_plus2_PZ;
    TBranch         *b_pi_plus2_ID;
    TBranch      *b_pi_plus2_PIDmu;
    TBranch      *b_pi_plus2_PIDK;
    TBranch      *b_pi_plus2_PIDp;
    TBranch      *b_pi_plus2_ProbNNk;
    TBranch      *b_pi_plus2_ProbNNp;
    TBranch      *b_pi_plus2_ProbNNpi;
    TBranch      *b_pi_plus2_ProbNNmu;
    TBranch      *b_pi_plus2_ProbNNghost;
    TBranch        *b_pi_plus2_isMuon;
    TBranch      *b_pi_plus2_TRACK_CHI2NDOF;
    TBranch      *b_pi_plus2_TRACK_GhostProb;
    TBranch      *b_pi_plus2_ptasy_1_00;
    
     TBranch        *b_Bs_OS_Muon_DEC;
     TBranch        *b_Bs_OS_Muon_PROB;
     TBranch        *b_Bs_OS_Electron_DEC;
     TBranch        *b_Bs_OS_Electron_PROB;
     TBranch        *b_Bs_OS_Kaon_DEC;
     TBranch        *b_Bs_OS_Kaon_PROB;
     TBranch        *b_Bs_OS_Kaon_PARTICLES_NUM;
     TBranch        *b_Bs_SS_Kaon_DEC;
     TBranch        *b_Bs_SS_Kaon_PROB;
     TBranch          *b_Bs_SS_Kaon_PARTICLES_NUM;
     TBranch        *b_Bs_SS_Pion_DEC;
     TBranch        *b_Bs_SS_Pion_PROB;
     TBranch          *b_Bs_SS_Pion_PARTICLES_NUM;
     TBranch        *b_Bs_SS_PionBDT_DEC;
     TBranch        *b_Bs_SS_PionBDT_PROB;
     TBranch        *b_Bs_VtxCharge_DEC;
     TBranch        *b_Bs_VtxCharge_PROB;
     TBranch        *b_Bs_OS_nnetKaon_DEC;
     TBranch        *b_Bs_OS_nnetKaon_PROB;
     TBranch        *b_Bs_SS_nnetKaon_DEC;
     TBranch        *b_Bs_SS_nnetKaon_PROB;
     TBranch        *b_Bs_SS_Proton_DEC;
     TBranch        *b_Bs_SS_Proton_PROB;
     TBranch        *b_Bs_OS_Charm_DEC;
     TBranch        *b_Bs_OS_Charm_PROB;
     
     /*
     TBranch          *b_Bs_BsTaggingTool_TAGDECISION_OS;
     TBranch       *b_Bs_BsTaggingTool_TAGOMEGA_OS;
     TBranch        *b_Bs_BsTaggingTool_OS_Muon_DEC;
     TBranch        *b_Bs_BsTaggingTool_OS_Muon_PROB;
     TBranch        *b_Bs_BsTaggingTool_OS_Electron_DEC;
     TBranch        *b_Bs_BsTaggingTool_OS_Electron_PROB;
     TBranch        *b_Bs_BsTaggingTool_OS_Kaon_DEC;
     TBranch        *b_Bs_BsTaggingTool_OS_Kaon_PROB;
     TBranch        *b_Bs_BsTaggingTool_SS_Kaon_DEC;
     TBranch        *b_Bs_BsTaggingTool_SS_Kaon_PROB;
     TBranch        *b_Bs_BsTaggingTool_SS_Pion_DEC;
     TBranch        *b_Bs_BsTaggingTool_SS_Pion_PROB;
     TBranch        *b_Bs_BsTaggingTool_SS_PionBDT_DEC;
     TBranch        *b_Bs_BsTaggingTool_SS_PionBDT_PROB;
     TBranch        *b_Bs_BsTaggingTool_VtxCharge_DEC;
     TBranch        *b_Bs_BsTaggingTool_VtxCharge_PROB;
     TBranch        *b_Bs_BsTaggingTool_OS_nnetKaon_DEC;
     TBranch        *b_Bs_BsTaggingTool_OS_nnetKaon_PROB;
     TBranch        *b_Bs_BsTaggingTool_SS_nnetKaon_DEC;
     TBranch        *b_Bs_BsTaggingTool_SS_nnetKaon_PROB;
     TBranch        *b_Bs_BsTaggingTool_SS_Proton_DEC;
     TBranch        *b_Bs_BsTaggingTool_SS_Proton_PROB;
     TBranch        *b_Bs_BsTaggingTool_OS_Charm_DEC;
     TBranch        *b_Bs_BsTaggingTool_OS_Charm_PROB;
     */
   
};

#endif

#ifdef MiniDecayTree_cxx
MiniDecayTree::MiniDecayTree(Decay::Type decay, Year::Type year, Ds_finalState::Type finalState, DataType::Type dataType, TString polarity, TString inFileLoc, TString outFileLoc, Bool_t charmLess, TString usePIDvar, Bool_t bkg, Bool_t ltu, Bool_t ss ) : DecayTree(decay,year,finalState,dataType, polarity, inFileLoc, outFileLoc, bkg, ltu, ss), _charmLess(charmLess), _usePIDvar(usePIDvar)
{
    _inFileName = _outFileName;
    if(!_data) _inFileName.ReplaceAll(".root","_PID.root");
    _outFileName = _outFileName.ReplaceAll(TString("Mini/"),TString("Preselected/"));   
    if(_charmLess) _outFileName.ReplaceAll(".root","_charmLess.root");
}

MiniDecayTree::~MiniDecayTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void MiniDecayTree::Init()
{
    
    TTree* tree = this->GetInputTree();
    cout << "Found files, now init" << endl;
    
    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);
        
    fChain->SetBranchAddress("Bs_MINIPNEXTBEST", &Bs_MINIPNEXTBEST, &b_Bs_MINIPNEXTBEST);
    fChain->SetBranchAddress("Bs_MINIPCHI2NEXTBEST", &Bs_MINIPCHI2NEXTBEST, &b_Bs_MINIPCHI2NEXTBEST);

    if(_Ds_finalState == Ds_finalState::phipi && _decay == Decay::signal){

    if(!_data){
	    fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
	    fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
    	    fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);

	    fChain->SetBranchAddress("K_plus_TRUEID", &K_plus_TRUEID);
    	    fChain->SetBranchAddress("pi_plus_TRUEID", &pi_plus_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
    	    fChain->SetBranchAddress("K_plus_fromDs_TRUEID", &K_plus_fromDs_TRUEID);
    	    fChain->SetBranchAddress("K_minus_fromDs_TRUEID", &K_minus_fromDs_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);

	    fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
	    fChain->SetBranchAddress("K_plus_MC_MOTHER_ID", &K_plus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus_MC_MOTHER_ID", &pi_plus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("K_plus_fromDs_MC_MOTHER_ID", &K_plus_fromDs_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("K_minus_fromDs_MC_MOTHER_ID", &K_minus_fromDs_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);

	    fChain->SetBranchAddress("K_plus_PIDK_gen_MagDown", &K_plus_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_plus_PIDK_gen_MagDown", &pi_plus_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagDown", &pi_minus_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_gen_MagDown", &K_plus_fromDs_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagDown", &K_minus_fromDs_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagDown", &pi_minus_fromDs_PIDK_gen_MagDown);

	    fChain->SetBranchAddress("K_plus_PIDK_gen_MagUp", &K_plus_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_plus_PIDK_gen_MagUp", &pi_plus_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagUp", &pi_minus_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_gen_MagUp", &K_plus_fromDs_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagUp", &K_minus_fromDs_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagUp", &pi_minus_fromDs_PIDK_gen_MagUp);

	    fChain->SetBranchAddress("K_plus_PIDK_corr_MagDown", &K_plus_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_plus_PIDK_corr_MagDown", &pi_plus_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagDown", &pi_minus_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_corr_MagDown", &K_plus_fromDs_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagDown", &K_minus_fromDs_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagDown", &pi_minus_fromDs_PIDK_corr_MagDown);

	    fChain->SetBranchAddress("K_plus_PIDK_corr_MagUp", &K_plus_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_plus_PIDK_corr_MagUp", &pi_plus_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagUp", &pi_minus_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_corr_MagUp", &K_plus_fromDs_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagUp", &K_minus_fromDs_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagUp", &pi_minus_fromDs_PIDK_corr_MagUp);

            /*
	    fChain->SetBranchAddress("K_plus_PIDp_gen_MagDown", &K_plus_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_plus_PIDp_gen_MagDown", &pi_plus_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDp_gen_MagDown", &pi_minus_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagDown", &K_plus_fromDs_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagDown", &K_minus_fromDs_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagDown", &pi_minus_fromDs_PIDp_gen_MagDown);

	    fChain->SetBranchAddress("K_plus_PIDp_gen_MagUp", &K_plus_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_plus_PIDp_gen_MagUp", &pi_plus_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDp_gen_MagUp", &pi_minus_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagUp", &K_plus_fromDs_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagUp", &K_minus_fromDs_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagUp", &pi_minus_fromDs_PIDp_gen_MagUp);

	    fChain->SetBranchAddress("K_plus_PIDp_corr_MagDown", &K_plus_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_plus_PIDp_corr_MagDown", &pi_plus_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDp_corr_MagDown", &pi_minus_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagDown", &K_plus_fromDs_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagDown", &K_minus_fromDs_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagDown", &pi_minus_fromDs_PIDp_corr_MagDown);

	    fChain->SetBranchAddress("K_plus_PIDp_corr_MagUp", &K_plus_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_plus_PIDp_corr_MagUp", &pi_plus_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDp_corr_MagUp", &pi_minus_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagUp", &K_plus_fromDs_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagUp", &K_minus_fromDs_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagUp", &pi_minus_fromDs_PIDp_corr_MagUp);
	    */
    }

    fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
    fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
    fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
    fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
    fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
    fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
    fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
    fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
    fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
    fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
    fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
    fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
    fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
    fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
    fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
    //fChain->SetBranchAddress("Bs_TAUCHI2", &Bs_TAUCHI2, &b_Bs_TAUCHI2);
    fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
    fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);
    fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
    fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
    //fChain->SetBranchAddress("Bs_L0GlobalDecision_TIS", &Bs_L0GlobalDecision_TIS, &b_Bs_L0GlobalDecision_TIS);
    //fChain->SetBranchAddress("Bs_L0GlobalDecision_TOS", &Bs_L0GlobalDecision_TOS, &b_Bs_L0GlobalDecision_TOS);
    /*fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TIS", &Bs_Hlt1TrackAllL0Decision_TIS, &b_Bs_Hlt1TrackAllL0Decision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TOS", &Bs_Hlt1TrackAllL0Decision_TOS, &b_Bs_Hlt1TrackAllL0Decision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TIS", &Bs_Hlt1TrackMVADecision_TIS, &b_Bs_Hlt1TrackMVADecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TOS", &Bs_Hlt1TrackMVADecision_TOS, &b_Bs_Hlt1TrackMVADecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TIS", &Bs_Hlt1TwoTrackMVADecision_TIS, &b_Bs_Hlt1TwoTrackMVADecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TOS", &Bs_Hlt1TwoTrackMVADecision_TOS, &b_Bs_Hlt1TwoTrackMVADecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TIS", &Bs_Hlt1TrackMVALooseDecision_TIS, &b_Bs_Hlt1TrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TOS", &Bs_Hlt1TrackMVALooseDecision_TOS, &b_Bs_Hlt1TrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TIS", &Bs_Hlt1TwoTrackMVALooseDecision_TIS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TOS", &Bs_Hlt1TwoTrackMVALooseDecision_TOS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TIS", &Bs_Hlt2IncPhiDecision_TIS, &b_Bs_Hlt2IncPhiDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TOS", &Bs_Hlt2IncPhiDecision_TOS, &b_Bs_Hlt2IncPhiDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TIS", &Bs_Hlt2PhiIncPhiDecision_TIS, &b_Bs_Hlt2PhiIncPhiDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TOS", &Bs_Hlt2PhiIncPhiDecision_TOS, &b_Bs_Hlt2PhiIncPhiDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TIS", &Bs_Hlt2Topo2BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TOS", &Bs_Hlt2Topo2BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TIS", &Bs_Hlt2Topo3BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TOS", &Bs_Hlt2Topo3BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TIS", &Bs_Hlt2Topo4BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TOS", &Bs_Hlt2Topo4BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TIS", &Bs_Hlt2Topo2BodyDecision_TIS, &b_Bs_Hlt2Topo2BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TOS", &Bs_Hlt2Topo2BodyDecision_TOS, &b_Bs_Hlt2Topo2BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TIS", &Bs_Hlt2Topo3BodyDecision_TIS, &b_Bs_Hlt2Topo3BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TOS", &Bs_Hlt2Topo3BodyDecision_TOS, &b_Bs_Hlt2Topo3BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TIS", &Bs_Hlt2Topo4BodyDecision_TIS, &b_Bs_Hlt2Topo4BodyDecision_TIS);
    fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TOS", &Bs_Hlt2Topo4BodyDecision_TOS, &b_Bs_Hlt2Topo4BodyDecision_TOS);
    fChain->SetBranchAddress("Bs_TAGDECISION", &Bs_TAGDECISION, &b_Bs_TAGDECISION);
    fChain->SetBranchAddress("Bs_TAGOMEGA", &Bs_TAGOMEGA, &b_Bs_TAGOMEGA);
    fChain->SetBranchAddress("Bs_TAGDECISION_OS", &Bs_TAGDECISION_OS, &b_Bs_TAGDECISION_OS);
    fChain->SetBranchAddress("Bs_TAGOMEGA_OS", &Bs_TAGOMEGA_OS, &b_Bs_TAGOMEGA_OS);
    fChain->SetBranchAddress("Bs_TAGGER", &Bs_TAGGER, &b_Bs_TAGGER);*/
    fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
/*    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION_OS", &Bs_BsTaggingTool_TAGDECISION_OS, &b_Bs_BsTaggingTool_TAGDECISION_OS);
    fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA_OS", &Bs_BsTaggingTool_TAGOMEGA_OS, &b_Bs_BsTaggingTool_TAGOMEGA_OS);*/
    
    fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
    fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
    fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
    fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
    fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
    fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
    fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
    fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
    fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
    fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
    fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
    fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
    fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
    fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
    fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
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
    fChain->SetBranchAddress("Ds_ID", &Ds_ID, &b_Ds_ID);
    //fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
    //fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);
    //fChain->SetBranchAddress("Ds_TAUCHI2", &Ds_TAUCHI2, &b_Ds_TAUCHI2);
    fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);

    fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
    fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
    fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
    fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
    fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
    fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
    
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
    
    fChain->SetBranchAddress("K_plus_fromDs_ETA", &K_plus_fromDs_ETA, &b_K_plus_fromDs_ETA);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNmu", &K_plus_fromDs_MC12TuneV2_ProbNNmu, &b_K_plus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNpi", &K_plus_fromDs_MC12TuneV2_ProbNNpi, &b_K_plus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNk", &K_plus_fromDs_MC12TuneV2_ProbNNk, &b_K_plus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNp", &K_plus_fromDs_MC12TuneV2_ProbNNp, &b_K_plus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNghost", &K_plus_fromDs_MC12TuneV2_ProbNNghost, &b_K_plus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNmu", &K_plus_fromDs_MC12TuneV3_ProbNNmu, &b_K_plus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNpi", &K_plus_fromDs_MC12TuneV3_ProbNNpi, &b_K_plus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNk", &K_plus_fromDs_MC12TuneV3_ProbNNk, &b_K_plus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNp", &K_plus_fromDs_MC12TuneV3_ProbNNp, &b_K_plus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNghost", &K_plus_fromDs_MC12TuneV3_ProbNNghost, &b_K_plus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNmu", &K_plus_fromDs_MC12TuneV4_ProbNNmu, &b_K_plus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNpi", &K_plus_fromDs_MC12TuneV4_ProbNNpi, &b_K_plus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNk", &K_plus_fromDs_MC12TuneV4_ProbNNk, &b_K_plus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNp", &K_plus_fromDs_MC12TuneV4_ProbNNp, &b_K_plus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV4_ProbNNghost", &K_plus_fromDs_MC12TuneV4_ProbNNghost, &b_K_plus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNmu", &K_plus_fromDs_MC15TuneV1_ProbNNmu, &b_K_plus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNpi", &K_plus_fromDs_MC15TuneV1_ProbNNpi, &b_K_plus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNk", &K_plus_fromDs_MC15TuneV1_ProbNNk, &b_K_plus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNp", &K_plus_fromDs_MC15TuneV1_ProbNNp, &b_K_plus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_MC15TuneV1_ProbNNghost", &K_plus_fromDs_MC15TuneV1_ProbNNghost, &b_K_plus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_IP_OWNPV", &K_plus_fromDs_IP_OWNPV, &b_K_plus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("K_plus_fromDs_IPCHI2_OWNPV", &K_plus_fromDs_IPCHI2_OWNPV, &b_K_plus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_plus_fromDs_P", &K_plus_fromDs_P, &b_K_plus_fromDs_P);
    fChain->SetBranchAddress("K_plus_fromDs_PT", &K_plus_fromDs_PT, &b_K_plus_fromDs_PT);
    fChain->SetBranchAddress("K_plus_fromDs_PE", &K_plus_fromDs_PE, &b_K_plus_fromDs_PE);
    fChain->SetBranchAddress("K_plus_fromDs_PX", &K_plus_fromDs_PX, &b_K_plus_fromDs_PX);
    fChain->SetBranchAddress("K_plus_fromDs_PY", &K_plus_fromDs_PY, &b_K_plus_fromDs_PY);
    fChain->SetBranchAddress("K_plus_fromDs_PZ", &K_plus_fromDs_PZ, &b_K_plus_fromDs_PZ);
    fChain->SetBranchAddress("K_plus_fromDs_ID", &K_plus_fromDs_ID, &b_K_plus_fromDs_ID);
    fChain->SetBranchAddress("K_plus_fromDs_PIDmu", &K_plus_fromDs_PIDmu, &b_K_plus_fromDs_PIDmu);
    fChain->SetBranchAddress("K_plus_fromDs_PIDK", &K_plus_fromDs_PIDK, &b_K_plus_fromDs_PIDK);
    fChain->SetBranchAddress("K_plus_fromDs_PIDp", &K_plus_fromDs_PIDp, &b_K_plus_fromDs_PIDp);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNk", &K_plus_fromDs_ProbNNk, &b_K_plus_fromDs_ProbNNk);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNp", &K_plus_fromDs_ProbNNp, &b_K_plus_fromDs_ProbNNp);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNpi", &K_plus_fromDs_ProbNNpi, &b_K_plus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNmu", &K_plus_fromDs_ProbNNmu, &b_K_plus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("K_plus_fromDs_ProbNNghost", &K_plus_fromDs_ProbNNghost, &b_K_plus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("K_plus_fromDs_isMuon", &K_plus_fromDs_isMuon, &b_K_plus_fromDs_isMuon);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_CHI2NDOF", &K_plus_fromDs_TRACK_CHI2NDOF, &b_K_plus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_plus_fromDs_TRACK_GhostProb", &K_plus_fromDs_TRACK_GhostProb, &b_K_plus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("K_plus_fromDs_ptasy_1.00", &K_plus_fromDs_ptasy_1_00, &b_K_plus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("K_minus_fromDs_ETA", &K_minus_fromDs_ETA, &b_K_minus_fromDs_ETA);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNmu", &K_minus_fromDs_MC12TuneV2_ProbNNmu, &b_K_minus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNpi", &K_minus_fromDs_MC12TuneV2_ProbNNpi, &b_K_minus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNk", &K_minus_fromDs_MC12TuneV2_ProbNNk, &b_K_minus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNp", &K_minus_fromDs_MC12TuneV2_ProbNNp, &b_K_minus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNghost", &K_minus_fromDs_MC12TuneV2_ProbNNghost, &b_K_minus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNmu", &K_minus_fromDs_MC12TuneV3_ProbNNmu, &b_K_minus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNpi", &K_minus_fromDs_MC12TuneV3_ProbNNpi, &b_K_minus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNk", &K_minus_fromDs_MC12TuneV3_ProbNNk, &b_K_minus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNp", &K_minus_fromDs_MC12TuneV3_ProbNNp, &b_K_minus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNghost", &K_minus_fromDs_MC12TuneV3_ProbNNghost, &b_K_minus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNmu", &K_minus_fromDs_MC12TuneV4_ProbNNmu, &b_K_minus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNpi", &K_minus_fromDs_MC12TuneV4_ProbNNpi, &b_K_minus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNk", &K_minus_fromDs_MC12TuneV4_ProbNNk, &b_K_minus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNp", &K_minus_fromDs_MC12TuneV4_ProbNNp, &b_K_minus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV4_ProbNNghost", &K_minus_fromDs_MC12TuneV4_ProbNNghost, &b_K_minus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNmu", &K_minus_fromDs_MC15TuneV1_ProbNNmu, &b_K_minus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNpi", &K_minus_fromDs_MC15TuneV1_ProbNNpi, &b_K_minus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNk", &K_minus_fromDs_MC15TuneV1_ProbNNk, &b_K_minus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNp", &K_minus_fromDs_MC15TuneV1_ProbNNp, &b_K_minus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_MC15TuneV1_ProbNNghost", &K_minus_fromDs_MC15TuneV1_ProbNNghost, &b_K_minus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_IP_OWNPV", &K_minus_fromDs_IP_OWNPV, &b_K_minus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("K_minus_fromDs_IPCHI2_OWNPV", &K_minus_fromDs_IPCHI2_OWNPV, &b_K_minus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_minus_fromDs_P", &K_minus_fromDs_P, &b_K_minus_fromDs_P);
    fChain->SetBranchAddress("K_minus_fromDs_PT", &K_minus_fromDs_PT, &b_K_minus_fromDs_PT);
    fChain->SetBranchAddress("K_minus_fromDs_PE", &K_minus_fromDs_PE, &b_K_minus_fromDs_PE);
    fChain->SetBranchAddress("K_minus_fromDs_PX", &K_minus_fromDs_PX, &b_K_minus_fromDs_PX);
    fChain->SetBranchAddress("K_minus_fromDs_PY", &K_minus_fromDs_PY, &b_K_minus_fromDs_PY);
    fChain->SetBranchAddress("K_minus_fromDs_PZ", &K_minus_fromDs_PZ, &b_K_minus_fromDs_PZ);
    fChain->SetBranchAddress("K_minus_fromDs_ID", &K_minus_fromDs_ID, &b_K_minus_fromDs_ID);
    fChain->SetBranchAddress("K_minus_fromDs_PIDmu", &K_minus_fromDs_PIDmu, &b_K_minus_fromDs_PIDmu);
    fChain->SetBranchAddress("K_minus_fromDs_PIDK", &K_minus_fromDs_PIDK, &b_K_minus_fromDs_PIDK);
    fChain->SetBranchAddress("K_minus_fromDs_PIDp", &K_minus_fromDs_PIDp, &b_K_minus_fromDs_PIDp);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNk", &K_minus_fromDs_ProbNNk, &b_K_minus_fromDs_ProbNNk);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNp", &K_minus_fromDs_ProbNNp, &b_K_minus_fromDs_ProbNNp);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNpi", &K_minus_fromDs_ProbNNpi, &b_K_minus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNmu", &K_minus_fromDs_ProbNNmu, &b_K_minus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("K_minus_fromDs_ProbNNghost", &K_minus_fromDs_ProbNNghost, &b_K_minus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("K_minus_fromDs_isMuon", &K_minus_fromDs_isMuon, &b_K_minus_fromDs_isMuon);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_CHI2NDOF", &K_minus_fromDs_TRACK_CHI2NDOF, &b_K_minus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb, &b_K_minus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("K_minus_fromDs_ptasy_1.00", &K_minus_fromDs_ptasy_1_00, &b_K_minus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNmu", &pi_minus_fromDs_MC12TuneV4_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNpi", &pi_minus_fromDs_MC12TuneV4_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNk", &pi_minus_fromDs_MC12TuneV4_ProbNNk, &b_pi_minus_fromDs_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNp", &pi_minus_fromDs_MC12TuneV4_ProbNNp, &b_pi_minus_fromDs_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNghost", &pi_minus_fromDs_MC12TuneV4_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNmu", &pi_minus_fromDs_MC15TuneV1_ProbNNmu, &b_pi_minus_fromDs_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNpi", &pi_minus_fromDs_MC15TuneV1_ProbNNpi, &b_pi_minus_fromDs_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNk", &pi_minus_fromDs_MC15TuneV1_ProbNNk, &b_pi_minus_fromDs_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNp", &pi_minus_fromDs_MC15TuneV1_ProbNNp, &b_pi_minus_fromDs_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNghost", &pi_minus_fromDs_MC15TuneV1_ProbNNghost, &b_pi_minus_fromDs_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
    fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
    fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
    fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
    fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
    fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
    fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
    fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
    fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CHI2NDOF", &pi_minus_fromDs_TRACK_CHI2NDOF, &b_pi_minus_fromDs_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_minus_fromDs_ptasy_1.00", &pi_minus_fromDs_ptasy_1_00, &b_pi_minus_fromDs_ptasy_1_00);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA1", &K_1_1270_plus_DOCA1, &b_K_1_1270_plus_DOCA1);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA2", &K_1_1270_plus_DOCA2, &b_K_1_1270_plus_DOCA2);
    fChain->SetBranchAddress("K_1_1270_plus_DOCA3", &K_1_1270_plus_DOCA3, &b_K_1_1270_plus_DOCA3);
    fChain->SetBranchAddress("K_1_1270_plus_ETA", &K_1_1270_plus_ETA, &b_K_1_1270_plus_ETA);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_X", &K_1_1270_plus_ENDVERTEX_X, &b_K_1_1270_plus_ENDVERTEX_X);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Y", &K_1_1270_plus_ENDVERTEX_Y, &b_K_1_1270_plus_ENDVERTEX_Y);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Z", &K_1_1270_plus_ENDVERTEX_Z, &b_K_1_1270_plus_ENDVERTEX_Z);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_XERR", &K_1_1270_plus_ENDVERTEX_XERR, &b_K_1_1270_plus_ENDVERTEX_XERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_YERR", &K_1_1270_plus_ENDVERTEX_YERR, &b_K_1_1270_plus_ENDVERTEX_YERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_ZERR", &K_1_1270_plus_ENDVERTEX_ZERR, &b_K_1_1270_plus_ENDVERTEX_ZERR);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_CHI2", &K_1_1270_plus_ENDVERTEX_CHI2, &b_K_1_1270_plus_ENDVERTEX_CHI2);
    fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_NDOF", &K_1_1270_plus_ENDVERTEX_NDOF, &b_K_1_1270_plus_ENDVERTEX_NDOF);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_X", &K_1_1270_plus_OWNPV_X, &b_K_1_1270_plus_OWNPV_X);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Y", &K_1_1270_plus_OWNPV_Y, &b_K_1_1270_plus_OWNPV_Y);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Z", &K_1_1270_plus_OWNPV_Z, &b_K_1_1270_plus_OWNPV_Z);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_XERR", &K_1_1270_plus_OWNPV_XERR, &b_K_1_1270_plus_OWNPV_XERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_YERR", &K_1_1270_plus_OWNPV_YERR, &b_K_1_1270_plus_OWNPV_YERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_ZERR", &K_1_1270_plus_OWNPV_ZERR, &b_K_1_1270_plus_OWNPV_ZERR);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_CHI2", &K_1_1270_plus_OWNPV_CHI2, &b_K_1_1270_plus_OWNPV_CHI2);
    fChain->SetBranchAddress("K_1_1270_plus_OWNPV_NDOF", &K_1_1270_plus_OWNPV_NDOF, &b_K_1_1270_plus_OWNPV_NDOF);
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
    fChain->SetBranchAddress("K_1_1270_plus_ID", &K_1_1270_plus_ID, &b_K_1_1270_plus_ID);
/*    fChain->SetBranchAddress("K_1_1270_plus_TAU", &K_1_1270_plus_TAU, &b_K_1_1270_plus_TAU);
    fChain->SetBranchAddress("K_1_1270_plus_TAUERR", &K_1_1270_plus_TAUERR, &b_K_1_1270_plus_TAUERR);*/
    //fChain->SetBranchAddress("K_1_1270_plus_TAUCHI2", &K_1_1270_plus_TAUCHI2, &b_K_1_1270_plus_TAUCHI2);
    fChain->SetBranchAddress("K_1_1270_plus_ptasy_1.00", &K_1_1270_plus_ptasy_1_00, &b_K_1_1270_plus_ptasy_1_00);
    fChain->SetBranchAddress("K_plus_ETA", &K_plus_ETA, &b_K_plus_ETA);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNmu", &K_plus_MC12TuneV2_ProbNNmu, &b_K_plus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNpi", &K_plus_MC12TuneV2_ProbNNpi, &b_K_plus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNk", &K_plus_MC12TuneV2_ProbNNk, &b_K_plus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNp", &K_plus_MC12TuneV2_ProbNNp, &b_K_plus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNghost", &K_plus_MC12TuneV2_ProbNNghost, &b_K_plus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNmu", &K_plus_MC12TuneV3_ProbNNmu, &b_K_plus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNpi", &K_plus_MC12TuneV3_ProbNNpi, &b_K_plus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNk", &K_plus_MC12TuneV3_ProbNNk, &b_K_plus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNp", &K_plus_MC12TuneV3_ProbNNp, &b_K_plus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNghost", &K_plus_MC12TuneV3_ProbNNghost, &b_K_plus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNmu", &K_plus_MC12TuneV4_ProbNNmu, &b_K_plus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNpi", &K_plus_MC12TuneV4_ProbNNpi, &b_K_plus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNk", &K_plus_MC12TuneV4_ProbNNk, &b_K_plus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNp", &K_plus_MC12TuneV4_ProbNNp, &b_K_plus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNghost", &K_plus_MC12TuneV4_ProbNNghost, &b_K_plus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNmu", &K_plus_MC15TuneV1_ProbNNmu, &b_K_plus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNpi", &K_plus_MC15TuneV1_ProbNNpi, &b_K_plus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNk", &K_plus_MC15TuneV1_ProbNNk, &b_K_plus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNp", &K_plus_MC15TuneV1_ProbNNp, &b_K_plus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNghost", &K_plus_MC15TuneV1_ProbNNghost, &b_K_plus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("K_plus_IP_OWNPV", &K_plus_IP_OWNPV, &b_K_plus_IP_OWNPV);
    fChain->SetBranchAddress("K_plus_IPCHI2_OWNPV", &K_plus_IPCHI2_OWNPV, &b_K_plus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("K_plus_P", &K_plus_P, &b_K_plus_P);
    fChain->SetBranchAddress("K_plus_PT", &K_plus_PT, &b_K_plus_PT);
    fChain->SetBranchAddress("K_plus_PE", &K_plus_PE, &b_K_plus_PE);
    fChain->SetBranchAddress("K_plus_PX", &K_plus_PX, &b_K_plus_PX);
    fChain->SetBranchAddress("K_plus_PY", &K_plus_PY, &b_K_plus_PY);
    fChain->SetBranchAddress("K_plus_PZ", &K_plus_PZ, &b_K_plus_PZ);
    fChain->SetBranchAddress("K_plus_ID", &K_plus_ID, &b_K_plus_ID);
    fChain->SetBranchAddress("K_plus_PIDmu", &K_plus_PIDmu, &b_K_plus_PIDmu);
    fChain->SetBranchAddress("K_plus_PIDK", &K_plus_PIDK, &b_K_plus_PIDK);
    fChain->SetBranchAddress("K_plus_PIDp", &K_plus_PIDp, &b_K_plus_PIDp);
    fChain->SetBranchAddress("K_plus_ProbNNk", &K_plus_ProbNNk, &b_K_plus_ProbNNk);
    fChain->SetBranchAddress("K_plus_ProbNNp", &K_plus_ProbNNp, &b_K_plus_ProbNNp);
    fChain->SetBranchAddress("K_plus_ProbNNpi", &K_plus_ProbNNpi, &b_K_plus_ProbNNpi);
    fChain->SetBranchAddress("K_plus_ProbNNmu", &K_plus_ProbNNmu, &b_K_plus_ProbNNmu);
    fChain->SetBranchAddress("K_plus_ProbNNghost", &K_plus_ProbNNghost, &b_K_plus_ProbNNghost);
    fChain->SetBranchAddress("K_plus_isMuon", &K_plus_isMuon, &b_K_plus_isMuon);
    fChain->SetBranchAddress("K_plus_TRACK_CHI2NDOF", &K_plus_TRACK_CHI2NDOF, &b_K_plus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("K_plus_TRACK_GhostProb", &K_plus_TRACK_GhostProb, &b_K_plus_TRACK_GhostProb);
    fChain->SetBranchAddress("K_plus_ptasy_1.00", &K_plus_ptasy_1_00, &b_K_plus_ptasy_1_00);
    fChain->SetBranchAddress("pi_plus_ETA", &pi_plus_ETA, &b_pi_plus_ETA);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNmu", &pi_plus_MC12TuneV2_ProbNNmu, &b_pi_plus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNpi", &pi_plus_MC12TuneV2_ProbNNpi, &b_pi_plus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNk", &pi_plus_MC12TuneV2_ProbNNk, &b_pi_plus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNp", &pi_plus_MC12TuneV2_ProbNNp, &b_pi_plus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNghost", &pi_plus_MC12TuneV2_ProbNNghost, &b_pi_plus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNmu", &pi_plus_MC12TuneV3_ProbNNmu, &b_pi_plus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNpi", &pi_plus_MC12TuneV3_ProbNNpi, &b_pi_plus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNk", &pi_plus_MC12TuneV3_ProbNNk, &b_pi_plus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNp", &pi_plus_MC12TuneV3_ProbNNp, &b_pi_plus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNghost", &pi_plus_MC12TuneV3_ProbNNghost, &b_pi_plus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNmu", &pi_plus_MC12TuneV4_ProbNNmu, &b_pi_plus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNpi", &pi_plus_MC12TuneV4_ProbNNpi, &b_pi_plus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNk", &pi_plus_MC12TuneV4_ProbNNk, &b_pi_plus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNp", &pi_plus_MC12TuneV4_ProbNNp, &b_pi_plus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNghost", &pi_plus_MC12TuneV4_ProbNNghost, &b_pi_plus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNmu", &pi_plus_MC15TuneV1_ProbNNmu, &b_pi_plus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNpi", &pi_plus_MC15TuneV1_ProbNNpi, &b_pi_plus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNk", &pi_plus_MC15TuneV1_ProbNNk, &b_pi_plus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNp", &pi_plus_MC15TuneV1_ProbNNp, &b_pi_plus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNghost", &pi_plus_MC15TuneV1_ProbNNghost, &b_pi_plus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_IP_OWNPV", &pi_plus_IP_OWNPV, &b_pi_plus_IP_OWNPV);
    fChain->SetBranchAddress("pi_plus_IPCHI2_OWNPV", &pi_plus_IPCHI2_OWNPV, &b_pi_plus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_plus_P", &pi_plus_P, &b_pi_plus_P);
    fChain->SetBranchAddress("pi_plus_PT", &pi_plus_PT, &b_pi_plus_PT);
    fChain->SetBranchAddress("pi_plus_PE", &pi_plus_PE, &b_pi_plus_PE);
    fChain->SetBranchAddress("pi_plus_PX", &pi_plus_PX, &b_pi_plus_PX);
    fChain->SetBranchAddress("pi_plus_PY", &pi_plus_PY, &b_pi_plus_PY);
    fChain->SetBranchAddress("pi_plus_PZ", &pi_plus_PZ, &b_pi_plus_PZ);
    fChain->SetBranchAddress("pi_plus_ID", &pi_plus_ID, &b_pi_plus_ID);
    fChain->SetBranchAddress("pi_plus_PIDmu", &pi_plus_PIDmu, &b_pi_plus_PIDmu);
    fChain->SetBranchAddress("pi_plus_PIDK", &pi_plus_PIDK, &b_pi_plus_PIDK);
    fChain->SetBranchAddress("pi_plus_PIDp", &pi_plus_PIDp, &b_pi_plus_PIDp);
    fChain->SetBranchAddress("pi_plus_ProbNNk", &pi_plus_ProbNNk, &b_pi_plus_ProbNNk);
    fChain->SetBranchAddress("pi_plus_ProbNNp", &pi_plus_ProbNNp, &b_pi_plus_ProbNNp);
    fChain->SetBranchAddress("pi_plus_ProbNNpi", &pi_plus_ProbNNpi, &b_pi_plus_ProbNNpi);
    fChain->SetBranchAddress("pi_plus_ProbNNmu", &pi_plus_ProbNNmu, &b_pi_plus_ProbNNmu);
    fChain->SetBranchAddress("pi_plus_ProbNNghost", &pi_plus_ProbNNghost, &b_pi_plus_ProbNNghost);
    fChain->SetBranchAddress("pi_plus_isMuon", &pi_plus_isMuon, &b_pi_plus_isMuon);
    fChain->SetBranchAddress("pi_plus_TRACK_CHI2NDOF", &pi_plus_TRACK_CHI2NDOF, &b_pi_plus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_plus_TRACK_GhostProb", &pi_plus_TRACK_GhostProb, &b_pi_plus_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_plus_ptasy_1.00", &pi_plus_ptasy_1_00, &b_pi_plus_ptasy_1_00);
    fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNmu", &pi_minus_MC12TuneV4_ProbNNmu, &b_pi_minus_MC12TuneV4_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNpi", &pi_minus_MC12TuneV4_ProbNNpi, &b_pi_minus_MC12TuneV4_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNk", &pi_minus_MC12TuneV4_ProbNNk, &b_pi_minus_MC12TuneV4_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNp", &pi_minus_MC12TuneV4_ProbNNp, &b_pi_minus_MC12TuneV4_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNghost", &pi_minus_MC12TuneV4_ProbNNghost, &b_pi_minus_MC12TuneV4_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNmu", &pi_minus_MC15TuneV1_ProbNNmu, &b_pi_minus_MC15TuneV1_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNpi", &pi_minus_MC15TuneV1_ProbNNpi, &b_pi_minus_MC15TuneV1_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNk", &pi_minus_MC15TuneV1_ProbNNk, &b_pi_minus_MC15TuneV1_ProbNNk);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNp", &pi_minus_MC15TuneV1_ProbNNp, &b_pi_minus_MC15TuneV1_ProbNNp);
    fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNghost", &pi_minus_MC15TuneV1_ProbNNghost, &b_pi_minus_MC15TuneV1_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
    fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
    fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
    fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
    fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
    fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
    fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
    fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
    fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
    fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
    fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
    fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
    fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
    fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
    fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
    fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
    fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
    fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
    fChain->SetBranchAddress("pi_minus_TRACK_CHI2NDOF", &pi_minus_TRACK_CHI2NDOF, &b_pi_minus_TRACK_CHI2NDOF);
    fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
    fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
    
    }

    if(_Ds_finalState == Ds_finalState::pipipi && _decay == Decay::signal){

	 if(!_data){
	    fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
	    fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
    	    fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);

	    fChain->SetBranchAddress("K_plus_TRUEID", &K_plus_TRUEID);
    	    fChain->SetBranchAddress("pi_plus_TRUEID", &pi_plus_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);
	    fChain->SetBranchAddress("pi_plus_fromDs_TRUEID", &pi_plus_fromDs_TRUEID);
   	    fChain->SetBranchAddress("pi_minus2_fromDs_TRUEID", &pi_minus2_fromDs_TRUEID);

	    fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
	    fChain->SetBranchAddress("K_plus_MC_MOTHER_ID", &K_plus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus_MC_MOTHER_ID", &pi_plus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);
	    fChain->SetBranchAddress("pi_plus_fromDs_MC_MOTHER_ID", &pi_plus_fromDs_MC_MOTHER_ID);
   	    fChain->SetBranchAddress("pi_minus2_fromDs_MC_MOTHER_ID", &pi_minus2_fromDs_MC_MOTHER_ID);
	}

        fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
        fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
        fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
        fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
        fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
        fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
        fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
        fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
        fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
        fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
        fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
        fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
        //fChain->SetBranchAddress("Bs_TAUCHI2", &Bs_TAUCHI2, &b_Bs_TAUCHI2);
        fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
        fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);/*
        fChain->SetBranchAddress("Bs_L0GlobalDecision_TIS", &Bs_L0GlobalDecision_TIS, &b_Bs_L0GlobalDecision_TIS);
        fChain->SetBranchAddress("Bs_L0GlobalDecision_TOS", &Bs_L0GlobalDecision_TOS, &b_Bs_L0GlobalDecision_TOS);*//*
        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TIS", &Bs_Hlt1TrackAllL0Decision_TIS, &b_Bs_Hlt1TrackAllL0Decision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TOS", &Bs_Hlt1TrackAllL0Decision_TOS, &b_Bs_Hlt1TrackAllL0Decision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TIS", &Bs_Hlt1TrackMVADecision_TIS, &b_Bs_Hlt1TrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TOS", &Bs_Hlt1TrackMVADecision_TOS, &b_Bs_Hlt1TrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TIS", &Bs_Hlt1TwoTrackMVADecision_TIS, &b_Bs_Hlt1TwoTrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TOS", &Bs_Hlt1TwoTrackMVADecision_TOS, &b_Bs_Hlt1TwoTrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TIS", &Bs_Hlt1TrackMVALooseDecision_TIS, &b_Bs_Hlt1TrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TOS", &Bs_Hlt1TrackMVALooseDecision_TOS, &b_Bs_Hlt1TrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TIS", &Bs_Hlt1TwoTrackMVALooseDecision_TIS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TOS", &Bs_Hlt1TwoTrackMVALooseDecision_TOS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TIS", &Bs_Hlt2IncPhiDecision_TIS, &b_Bs_Hlt2IncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TOS", &Bs_Hlt2IncPhiDecision_TOS, &b_Bs_Hlt2IncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TIS", &Bs_Hlt2PhiIncPhiDecision_TIS, &b_Bs_Hlt2PhiIncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TOS", &Bs_Hlt2PhiIncPhiDecision_TOS, &b_Bs_Hlt2PhiIncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TIS", &Bs_Hlt2Topo2BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TOS", &Bs_Hlt2Topo2BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TIS", &Bs_Hlt2Topo3BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TOS", &Bs_Hlt2Topo3BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TIS", &Bs_Hlt2Topo4BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TOS", &Bs_Hlt2Topo4BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TIS", &Bs_Hlt2Topo2BodyDecision_TIS, &b_Bs_Hlt2Topo2BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TOS", &Bs_Hlt2Topo2BodyDecision_TOS, &b_Bs_Hlt2Topo2BodyDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TIS", &Bs_Hlt2Topo3BodyDecision_TIS, &b_Bs_Hlt2Topo3BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TOS", &Bs_Hlt2Topo3BodyDecision_TOS, &b_Bs_Hlt2Topo3BodyDecision_TOS);
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
        fChain->SetBranchAddress("Bs_OS_Charm_PROB", &Bs_OS_Charm_PROB, &b_Bs_OS_Charm_PROB);*/
        fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
        fChain->SetBranchAddress("Bs_B0DTF_nPV", &Bs_B0DTF_nPV, &b_Bs_B0DTF_nPV);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_M", Bs_B0DTF_D_splus_M, &b_Bs_B0DTF_D_splus_M);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_MERR", Bs_B0DTF_D_splus_MERR, &b_Bs_B0DTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_P", Bs_B0DTF_D_splus_P, &b_Bs_B0DTF_D_splus_P);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_PERR", Bs_B0DTF_D_splus_PERR, &b_Bs_B0DTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctau", Bs_B0DTF_D_splus_ctau, &b_Bs_B0DTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctauErr", Bs_B0DTF_D_splus_ctauErr, &b_Bs_B0DTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLength", Bs_B0DTF_D_splus_decayLength, &b_Bs_B0DTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLengthErr", Bs_B0DTF_D_splus_decayLengthErr, &b_Bs_B0DTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_ID", Bs_B0DTF_D_splus_piplus_0_ID, &b_Bs_B0DTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PE", Bs_B0DTF_D_splus_piplus_0_PE, &b_Bs_B0DTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PX", Bs_B0DTF_D_splus_piplus_0_PX, &b_Bs_B0DTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PY", Bs_B0DTF_D_splus_piplus_0_PY, &b_Bs_B0DTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PZ", Bs_B0DTF_D_splus_piplus_0_PZ, &b_Bs_B0DTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_ID", Bs_B0DTF_D_splus_piplus_1_ID, &b_Bs_B0DTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PE", Bs_B0DTF_D_splus_piplus_1_PE, &b_Bs_B0DTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PX", Bs_B0DTF_D_splus_piplus_1_PX, &b_Bs_B0DTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PY", Bs_B0DTF_D_splus_piplus_1_PY, &b_Bs_B0DTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PZ", Bs_B0DTF_D_splus_piplus_1_PZ, &b_Bs_B0DTF_D_splus_piplus_1_PZ);
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
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_M", Bs_BsDTF_D_splus_M, &b_Bs_BsDTF_D_splus_M);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_MERR", Bs_BsDTF_D_splus_MERR, &b_Bs_BsDTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_P", Bs_BsDTF_D_splus_P, &b_Bs_BsDTF_D_splus_P);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_PERR", Bs_BsDTF_D_splus_PERR, &b_Bs_BsDTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctau", Bs_BsDTF_D_splus_ctau, &b_Bs_BsDTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctauErr", Bs_BsDTF_D_splus_ctauErr, &b_Bs_BsDTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLength", Bs_BsDTF_D_splus_decayLength, &b_Bs_BsDTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLengthErr", Bs_BsDTF_D_splus_decayLengthErr, &b_Bs_BsDTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_ID", Bs_BsDTF_D_splus_piplus_0_ID, &b_Bs_BsDTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PE", Bs_BsDTF_D_splus_piplus_0_PE, &b_Bs_BsDTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PX", Bs_BsDTF_D_splus_piplus_0_PX, &b_Bs_BsDTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PY", Bs_BsDTF_D_splus_piplus_0_PY, &b_Bs_BsDTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PZ", Bs_BsDTF_D_splus_piplus_0_PZ, &b_Bs_BsDTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_ID", Bs_BsDTF_D_splus_piplus_1_ID, &b_Bs_BsDTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PE", Bs_BsDTF_D_splus_piplus_1_PE, &b_Bs_BsDTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PX", Bs_BsDTF_D_splus_piplus_1_PX, &b_Bs_BsDTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PY", Bs_BsDTF_D_splus_piplus_1_PY, &b_Bs_BsDTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PZ", Bs_BsDTF_D_splus_piplus_1_PZ, &b_Bs_BsDTF_D_splus_piplus_1_PZ);
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
        fChain->SetBranchAddress("Bs_DTF_D_splus_M", Bs_DTF_D_splus_M, &b_Bs_DTF_D_splus_M);
        fChain->SetBranchAddress("Bs_DTF_D_splus_MERR", Bs_DTF_D_splus_MERR, &b_Bs_DTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_DTF_D_splus_P", Bs_DTF_D_splus_P, &b_Bs_DTF_D_splus_P);
        fChain->SetBranchAddress("Bs_DTF_D_splus_PERR", Bs_DTF_D_splus_PERR, &b_Bs_DTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_DTF_D_splus_ctau", Bs_DTF_D_splus_ctau, &b_Bs_DTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_DTF_D_splus_ctauErr", Bs_DTF_D_splus_ctauErr, &b_Bs_DTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_D_splus_decayLength", Bs_DTF_D_splus_decayLength, &b_Bs_DTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_DTF_D_splus_decayLengthErr", Bs_DTF_D_splus_decayLengthErr, &b_Bs_DTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_ID", Bs_DTF_D_splus_piplus_0_ID, &b_Bs_DTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PE", Bs_DTF_D_splus_piplus_0_PE, &b_Bs_DTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PX", Bs_DTF_D_splus_piplus_0_PX, &b_Bs_DTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PY", Bs_DTF_D_splus_piplus_0_PY, &b_Bs_DTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PZ", Bs_DTF_D_splus_piplus_0_PZ, &b_Bs_DTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_ID", Bs_DTF_D_splus_piplus_1_ID, &b_Bs_DTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PE", Bs_DTF_D_splus_piplus_1_PE, &b_Bs_DTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PX", Bs_DTF_D_splus_piplus_1_PX, &b_Bs_DTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PY", Bs_DTF_D_splus_piplus_1_PY, &b_Bs_DTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PZ", Bs_DTF_D_splus_piplus_1_PZ, &b_Bs_DTF_D_splus_piplus_1_PZ);
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
        fChain->SetBranchAddress("Bs_PV_Dplus_M", Bs_PV_Dplus_M, &b_Bs_PV_Dplus_M);
        fChain->SetBranchAddress("Bs_PV_Dplus_MERR", Bs_PV_Dplus_MERR, &b_Bs_PV_Dplus_MERR);
        fChain->SetBranchAddress("Bs_PV_Dplus_P", Bs_PV_Dplus_P, &b_Bs_PV_Dplus_P);
        fChain->SetBranchAddress("Bs_PV_Dplus_PERR", Bs_PV_Dplus_PERR, &b_Bs_PV_Dplus_PERR);
        fChain->SetBranchAddress("Bs_PV_Dplus_ctau", Bs_PV_Dplus_ctau, &b_Bs_PV_Dplus_ctau);
        fChain->SetBranchAddress("Bs_PV_Dplus_ctauErr", Bs_PV_Dplus_ctauErr, &b_Bs_PV_Dplus_ctauErr);
        fChain->SetBranchAddress("Bs_PV_Dplus_decayLength", Bs_PV_Dplus_decayLength, &b_Bs_PV_Dplus_decayLength);
        fChain->SetBranchAddress("Bs_PV_Dplus_decayLengthErr", Bs_PV_Dplus_decayLengthErr, &b_Bs_PV_Dplus_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_ID", Bs_PV_Dplus_piplus_0_ID, &b_Bs_PV_Dplus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PE", Bs_PV_Dplus_piplus_0_PE, &b_Bs_PV_Dplus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PX", Bs_PV_Dplus_piplus_0_PX, &b_Bs_PV_Dplus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PY", Bs_PV_Dplus_piplus_0_PY, &b_Bs_PV_Dplus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PZ", Bs_PV_Dplus_piplus_0_PZ, &b_Bs_PV_Dplus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_ID", Bs_PV_Dplus_piplus_1_ID, &b_Bs_PV_Dplus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PE", Bs_PV_Dplus_piplus_1_PE, &b_Bs_PV_Dplus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PX", Bs_PV_Dplus_piplus_1_PX, &b_Bs_PV_Dplus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PY", Bs_PV_Dplus_piplus_1_PY, &b_Bs_PV_Dplus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PZ", Bs_PV_Dplus_piplus_1_PZ, &b_Bs_PV_Dplus_piplus_1_PZ);
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
        fChain->SetBranchAddress("Bs_PV_status", Bs_PV_status, &b_Bs_PV_status);/*
        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION_OS", &Bs_BsTaggingTool_TAGDECISION_OS, &b_Bs_BsTaggingTool_TAGDECISION_OS);
        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA_OS", &Bs_BsTaggingTool_TAGOMEGA_OS, &b_Bs_BsTaggingTool_TAGOMEGA_OS);
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
        fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Charm_PROB", &Bs_BsTaggingTool_OS_Charm_PROB, &b_Bs_BsTaggingTool_OS_Charm_PROB);*/
        fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
        fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
        fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
        fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
        fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
        fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
        fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
        fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
        fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
        fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
        fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
        fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
        fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Ds_ID", &Ds_ID, &b_Ds_ID);/*
        fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
        fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);*/
        //fChain->SetBranchAddress("Ds_TAUCHI2", &Ds_TAUCHI2, &b_Ds_TAUCHI2);
        fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus_fromDs_ETA", &pi_plus_fromDs_ETA, &b_pi_plus_fromDs_ETA);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNmu", &pi_plus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNpi", &pi_plus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNk", &pi_plus_fromDs_MC12TuneV2_ProbNNk, &b_pi_plus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNp", &pi_plus_fromDs_MC12TuneV2_ProbNNp, &b_pi_plus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNghost", &pi_plus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNmu", &pi_plus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNpi", &pi_plus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNk", &pi_plus_fromDs_MC12TuneV3_ProbNNk, &b_pi_plus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNp", &pi_plus_fromDs_MC12TuneV3_ProbNNp, &b_pi_plus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNghost", &pi_plus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV4_ProbNNmu", &pi_plus_fromDs_MC12TuneV4_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV4_ProbNNpi", &pi_plus_fromDs_MC12TuneV4_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV4_ProbNNk", &pi_plus_fromDs_MC12TuneV4_ProbNNk, &b_pi_plus_fromDs_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV4_ProbNNp", &pi_plus_fromDs_MC12TuneV4_ProbNNp, &b_pi_plus_fromDs_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV4_ProbNNghost", &pi_plus_fromDs_MC12TuneV4_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_MC15TuneV1_ProbNNmu", &pi_plus_fromDs_MC15TuneV1_ProbNNmu, &b_pi_plus_fromDs_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC15TuneV1_ProbNNpi", &pi_plus_fromDs_MC15TuneV1_ProbNNpi, &b_pi_plus_fromDs_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC15TuneV1_ProbNNk", &pi_plus_fromDs_MC15TuneV1_ProbNNk, &b_pi_plus_fromDs_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC15TuneV1_ProbNNp", &pi_plus_fromDs_MC15TuneV1_ProbNNp, &b_pi_plus_fromDs_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC15TuneV1_ProbNNghost", &pi_plus_fromDs_MC15TuneV1_ProbNNghost, &b_pi_plus_fromDs_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_IP_OWNPV", &pi_plus_fromDs_IP_OWNPV, &b_pi_plus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus_fromDs_IPCHI2_OWNPV", &pi_plus_fromDs_IPCHI2_OWNPV, &b_pi_plus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus_fromDs_P", &pi_plus_fromDs_P, &b_pi_plus_fromDs_P);
        fChain->SetBranchAddress("pi_plus_fromDs_PT", &pi_plus_fromDs_PT, &b_pi_plus_fromDs_PT);
        fChain->SetBranchAddress("pi_plus_fromDs_PE", &pi_plus_fromDs_PE, &b_pi_plus_fromDs_PE);
        fChain->SetBranchAddress("pi_plus_fromDs_PX", &pi_plus_fromDs_PX, &b_pi_plus_fromDs_PX);
        fChain->SetBranchAddress("pi_plus_fromDs_PY", &pi_plus_fromDs_PY, &b_pi_plus_fromDs_PY);
        fChain->SetBranchAddress("pi_plus_fromDs_PZ", &pi_plus_fromDs_PZ, &b_pi_plus_fromDs_PZ);
        fChain->SetBranchAddress("pi_plus_fromDs_ID", &pi_plus_fromDs_ID, &b_pi_plus_fromDs_ID);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDmu", &pi_plus_fromDs_PIDmu, &b_pi_plus_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDK", &pi_plus_fromDs_PIDK, &b_pi_plus_fromDs_PIDK);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDp", &pi_plus_fromDs_PIDp, &b_pi_plus_fromDs_PIDp);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNk", &pi_plus_fromDs_ProbNNk, &b_pi_plus_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNp", &pi_plus_fromDs_ProbNNp, &b_pi_plus_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNpi", &pi_plus_fromDs_ProbNNpi, &b_pi_plus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNmu", &pi_plus_fromDs_ProbNNmu, &b_pi_plus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNghost", &pi_plus_fromDs_ProbNNghost, &b_pi_plus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_isMuon", &pi_plus_fromDs_isMuon, &b_pi_plus_fromDs_isMuon);
        fChain->SetBranchAddress("pi_plus_fromDs_TRACK_CHI2NDOF", &pi_plus_fromDs_TRACK_CHI2NDOF, &b_pi_plus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus_fromDs_TRACK_GhostProb", &pi_plus_fromDs_TRACK_GhostProb, &b_pi_plus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus_fromDs_ptasy_1.00", &pi_plus_fromDs_ptasy_1_00, &b_pi_plus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNmu", &pi_minus_fromDs_MC12TuneV4_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNpi", &pi_minus_fromDs_MC12TuneV4_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNk", &pi_minus_fromDs_MC12TuneV4_ProbNNk, &b_pi_minus_fromDs_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNp", &pi_minus_fromDs_MC12TuneV4_ProbNNp, &b_pi_minus_fromDs_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV4_ProbNNghost", &pi_minus_fromDs_MC12TuneV4_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNmu", &pi_minus_fromDs_MC15TuneV1_ProbNNmu, &b_pi_minus_fromDs_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNpi", &pi_minus_fromDs_MC15TuneV1_ProbNNpi, &b_pi_minus_fromDs_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNk", &pi_minus_fromDs_MC15TuneV1_ProbNNk, &b_pi_minus_fromDs_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNp", &pi_minus_fromDs_MC15TuneV1_ProbNNp, &b_pi_minus_fromDs_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC15TuneV1_ProbNNghost", &pi_minus_fromDs_MC15TuneV1_ProbNNghost, &b_pi_minus_fromDs_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
        fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
        fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
        fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
        fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
        fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
        fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CHI2NDOF", &pi_minus_fromDs_TRACK_CHI2NDOF, &b_pi_minus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_fromDs_ptasy_1.00", &pi_minus_fromDs_ptasy_1_00, &b_pi_minus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus2_fromDs_ETA", &pi_minus2_fromDs_ETA, &b_pi_minus2_fromDs_ETA);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNmu", &pi_minus2_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNpi", &pi_minus2_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNk", &pi_minus2_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNp", &pi_minus2_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNghost", &pi_minus2_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNmu", &pi_minus2_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNpi", &pi_minus2_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNk", &pi_minus2_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNp", &pi_minus2_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNghost", &pi_minus2_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV4_ProbNNmu", &pi_minus2_fromDs_MC12TuneV4_ProbNNmu, &b_pi_minus2_fromDs_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV4_ProbNNpi", &pi_minus2_fromDs_MC12TuneV4_ProbNNpi, &b_pi_minus2_fromDs_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV4_ProbNNk", &pi_minus2_fromDs_MC12TuneV4_ProbNNk, &b_pi_minus2_fromDs_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV4_ProbNNp", &pi_minus2_fromDs_MC12TuneV4_ProbNNp, &b_pi_minus2_fromDs_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV4_ProbNNghost", &pi_minus2_fromDs_MC12TuneV4_ProbNNghost, &b_pi_minus2_fromDs_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC15TuneV1_ProbNNmu", &pi_minus2_fromDs_MC15TuneV1_ProbNNmu, &b_pi_minus2_fromDs_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC15TuneV1_ProbNNpi", &pi_minus2_fromDs_MC15TuneV1_ProbNNpi, &b_pi_minus2_fromDs_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC15TuneV1_ProbNNk", &pi_minus2_fromDs_MC15TuneV1_ProbNNk, &b_pi_minus2_fromDs_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC15TuneV1_ProbNNp", &pi_minus2_fromDs_MC15TuneV1_ProbNNp, &b_pi_minus2_fromDs_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC15TuneV1_ProbNNghost", &pi_minus2_fromDs_MC15TuneV1_ProbNNghost, &b_pi_minus2_fromDs_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_IP_OWNPV", &pi_minus2_fromDs_IP_OWNPV, &b_pi_minus2_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus2_fromDs_IPCHI2_OWNPV", &pi_minus2_fromDs_IPCHI2_OWNPV, &b_pi_minus2_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus2_fromDs_P", &pi_minus2_fromDs_P, &b_pi_minus2_fromDs_P);
        fChain->SetBranchAddress("pi_minus2_fromDs_PT", &pi_minus2_fromDs_PT, &b_pi_minus2_fromDs_PT);
        fChain->SetBranchAddress("pi_minus2_fromDs_PE", &pi_minus2_fromDs_PE, &b_pi_minus2_fromDs_PE);
        fChain->SetBranchAddress("pi_minus2_fromDs_PX", &pi_minus2_fromDs_PX, &b_pi_minus2_fromDs_PX);
        fChain->SetBranchAddress("pi_minus2_fromDs_PY", &pi_minus2_fromDs_PY, &b_pi_minus2_fromDs_PY);
        fChain->SetBranchAddress("pi_minus2_fromDs_PZ", &pi_minus2_fromDs_PZ, &b_pi_minus2_fromDs_PZ);
        fChain->SetBranchAddress("pi_minus2_fromDs_ID", &pi_minus2_fromDs_ID, &b_pi_minus2_fromDs_ID);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDmu", &pi_minus2_fromDs_PIDmu, &b_pi_minus2_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDK", &pi_minus2_fromDs_PIDK, &b_pi_minus2_fromDs_PIDK);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDp", &pi_minus2_fromDs_PIDp, &b_pi_minus2_fromDs_PIDp);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNk", &pi_minus2_fromDs_ProbNNk, &b_pi_minus2_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNp", &pi_minus2_fromDs_ProbNNp, &b_pi_minus2_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNpi", &pi_minus2_fromDs_ProbNNpi, &b_pi_minus2_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNmu", &pi_minus2_fromDs_ProbNNmu, &b_pi_minus2_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNghost", &pi_minus2_fromDs_ProbNNghost, &b_pi_minus2_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_isMuon", &pi_minus2_fromDs_isMuon, &b_pi_minus2_fromDs_isMuon);
        fChain->SetBranchAddress("pi_minus2_fromDs_TRACK_CHI2NDOF", &pi_minus2_fromDs_TRACK_CHI2NDOF, &b_pi_minus2_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus2_fromDs_TRACK_GhostProb", &pi_minus2_fromDs_TRACK_GhostProb, &b_pi_minus2_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus2_fromDs_ptasy_1.00", &pi_minus2_fromDs_ptasy_1_00, &b_pi_minus2_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("K_1_1270_plus_DOCA1", &K_1_1270_plus_DOCA1, &b_K_1_1270_plus_DOCA1);
        fChain->SetBranchAddress("K_1_1270_plus_DOCA2", &K_1_1270_plus_DOCA2, &b_K_1_1270_plus_DOCA2);
        fChain->SetBranchAddress("K_1_1270_plus_DOCA3", &K_1_1270_plus_DOCA3, &b_K_1_1270_plus_DOCA3);
        fChain->SetBranchAddress("K_1_1270_plus_ETA", &K_1_1270_plus_ETA, &b_K_1_1270_plus_ETA);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_X", &K_1_1270_plus_ENDVERTEX_X, &b_K_1_1270_plus_ENDVERTEX_X);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Y", &K_1_1270_plus_ENDVERTEX_Y, &b_K_1_1270_plus_ENDVERTEX_Y);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_Z", &K_1_1270_plus_ENDVERTEX_Z, &b_K_1_1270_plus_ENDVERTEX_Z);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_XERR", &K_1_1270_plus_ENDVERTEX_XERR, &b_K_1_1270_plus_ENDVERTEX_XERR);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_YERR", &K_1_1270_plus_ENDVERTEX_YERR, &b_K_1_1270_plus_ENDVERTEX_YERR);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_ZERR", &K_1_1270_plus_ENDVERTEX_ZERR, &b_K_1_1270_plus_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_CHI2", &K_1_1270_plus_ENDVERTEX_CHI2, &b_K_1_1270_plus_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("K_1_1270_plus_ENDVERTEX_NDOF", &K_1_1270_plus_ENDVERTEX_NDOF, &b_K_1_1270_plus_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_X", &K_1_1270_plus_OWNPV_X, &b_K_1_1270_plus_OWNPV_X);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Y", &K_1_1270_plus_OWNPV_Y, &b_K_1_1270_plus_OWNPV_Y);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_Z", &K_1_1270_plus_OWNPV_Z, &b_K_1_1270_plus_OWNPV_Z);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_XERR", &K_1_1270_plus_OWNPV_XERR, &b_K_1_1270_plus_OWNPV_XERR);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_YERR", &K_1_1270_plus_OWNPV_YERR, &b_K_1_1270_plus_OWNPV_YERR);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_ZERR", &K_1_1270_plus_OWNPV_ZERR, &b_K_1_1270_plus_OWNPV_ZERR);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_CHI2", &K_1_1270_plus_OWNPV_CHI2, &b_K_1_1270_plus_OWNPV_CHI2);
        fChain->SetBranchAddress("K_1_1270_plus_OWNPV_NDOF", &K_1_1270_plus_OWNPV_NDOF, &b_K_1_1270_plus_OWNPV_NDOF);
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
        fChain->SetBranchAddress("K_1_1270_plus_ID", &K_1_1270_plus_ID, &b_K_1_1270_plus_ID);/*
        fChain->SetBranchAddress("K_1_1270_plus_TAU", &K_1_1270_plus_TAU, &b_K_1_1270_plus_TAU);
        fChain->SetBranchAddress("K_1_1270_plus_TAUERR", &K_1_1270_plus_TAUERR, &b_K_1_1270_plus_TAUERR);*/
        //fChain->SetBranchAddress("K_1_1270_plus_TAUCHI2", &K_1_1270_plus_TAUCHI2, &b_K_1_1270_plus_TAUCHI2);
        fChain->SetBranchAddress("K_1_1270_plus_ptasy_1.00", &K_1_1270_plus_ptasy_1_00, &b_K_1_1270_plus_ptasy_1_00);
        fChain->SetBranchAddress("K_plus_ETA", &K_plus_ETA, &b_K_plus_ETA);
        fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNmu", &K_plus_MC12TuneV2_ProbNNmu, &b_K_plus_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNpi", &K_plus_MC12TuneV2_ProbNNpi, &b_K_plus_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNk", &K_plus_MC12TuneV2_ProbNNk, &b_K_plus_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNp", &K_plus_MC12TuneV2_ProbNNp, &b_K_plus_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNghost", &K_plus_MC12TuneV2_ProbNNghost, &b_K_plus_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNmu", &K_plus_MC12TuneV3_ProbNNmu, &b_K_plus_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNpi", &K_plus_MC12TuneV3_ProbNNpi, &b_K_plus_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNk", &K_plus_MC12TuneV3_ProbNNk, &b_K_plus_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNp", &K_plus_MC12TuneV3_ProbNNp, &b_K_plus_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNghost", &K_plus_MC12TuneV3_ProbNNghost, &b_K_plus_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNmu", &K_plus_MC12TuneV4_ProbNNmu, &b_K_plus_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNpi", &K_plus_MC12TuneV4_ProbNNpi, &b_K_plus_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNk", &K_plus_MC12TuneV4_ProbNNk, &b_K_plus_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNp", &K_plus_MC12TuneV4_ProbNNp, &b_K_plus_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("K_plus_MC12TuneV4_ProbNNghost", &K_plus_MC12TuneV4_ProbNNghost, &b_K_plus_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNmu", &K_plus_MC15TuneV1_ProbNNmu, &b_K_plus_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNpi", &K_plus_MC15TuneV1_ProbNNpi, &b_K_plus_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNk", &K_plus_MC15TuneV1_ProbNNk, &b_K_plus_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNp", &K_plus_MC15TuneV1_ProbNNp, &b_K_plus_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("K_plus_MC15TuneV1_ProbNNghost", &K_plus_MC15TuneV1_ProbNNghost, &b_K_plus_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("K_plus_IP_OWNPV", &K_plus_IP_OWNPV, &b_K_plus_IP_OWNPV);
        fChain->SetBranchAddress("K_plus_IPCHI2_OWNPV", &K_plus_IPCHI2_OWNPV, &b_K_plus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("K_plus_P", &K_plus_P, &b_K_plus_P);
        fChain->SetBranchAddress("K_plus_PT", &K_plus_PT, &b_K_plus_PT);
        fChain->SetBranchAddress("K_plus_PE", &K_plus_PE, &b_K_plus_PE);
        fChain->SetBranchAddress("K_plus_PX", &K_plus_PX, &b_K_plus_PX);
        fChain->SetBranchAddress("K_plus_PY", &K_plus_PY, &b_K_plus_PY);
        fChain->SetBranchAddress("K_plus_PZ", &K_plus_PZ, &b_K_plus_PZ);
        fChain->SetBranchAddress("K_plus_ID", &K_plus_ID, &b_K_plus_ID);
        fChain->SetBranchAddress("K_plus_PIDmu", &K_plus_PIDmu, &b_K_plus_PIDmu);
        fChain->SetBranchAddress("K_plus_PIDK", &K_plus_PIDK, &b_K_plus_PIDK);
        fChain->SetBranchAddress("K_plus_PIDp", &K_plus_PIDp, &b_K_plus_PIDp);
        fChain->SetBranchAddress("K_plus_ProbNNk", &K_plus_ProbNNk, &b_K_plus_ProbNNk);
        fChain->SetBranchAddress("K_plus_ProbNNp", &K_plus_ProbNNp, &b_K_plus_ProbNNp);
        fChain->SetBranchAddress("K_plus_ProbNNpi", &K_plus_ProbNNpi, &b_K_plus_ProbNNpi);
        fChain->SetBranchAddress("K_plus_ProbNNmu", &K_plus_ProbNNmu, &b_K_plus_ProbNNmu);
        fChain->SetBranchAddress("K_plus_ProbNNghost", &K_plus_ProbNNghost, &b_K_plus_ProbNNghost);
        fChain->SetBranchAddress("K_plus_isMuon", &K_plus_isMuon, &b_K_plus_isMuon);
        fChain->SetBranchAddress("K_plus_TRACK_CHI2NDOF", &K_plus_TRACK_CHI2NDOF, &b_K_plus_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("K_plus_TRACK_GhostProb", &K_plus_TRACK_GhostProb, &b_K_plus_TRACK_GhostProb);
        fChain->SetBranchAddress("K_plus_ptasy_1.00", &K_plus_ptasy_1_00, &b_K_plus_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus_ETA", &pi_plus_ETA, &b_pi_plus_ETA);
        fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNmu", &pi_plus_MC12TuneV2_ProbNNmu, &b_pi_plus_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNpi", &pi_plus_MC12TuneV2_ProbNNpi, &b_pi_plus_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNk", &pi_plus_MC12TuneV2_ProbNNk, &b_pi_plus_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNp", &pi_plus_MC12TuneV2_ProbNNp, &b_pi_plus_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNghost", &pi_plus_MC12TuneV2_ProbNNghost, &b_pi_plus_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNmu", &pi_plus_MC12TuneV3_ProbNNmu, &b_pi_plus_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNpi", &pi_plus_MC12TuneV3_ProbNNpi, &b_pi_plus_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNk", &pi_plus_MC12TuneV3_ProbNNk, &b_pi_plus_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNp", &pi_plus_MC12TuneV3_ProbNNp, &b_pi_plus_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNghost", &pi_plus_MC12TuneV3_ProbNNghost, &b_pi_plus_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNmu", &pi_plus_MC12TuneV4_ProbNNmu, &b_pi_plus_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNpi", &pi_plus_MC12TuneV4_ProbNNpi, &b_pi_plus_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNk", &pi_plus_MC12TuneV4_ProbNNk, &b_pi_plus_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNp", &pi_plus_MC12TuneV4_ProbNNp, &b_pi_plus_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("pi_plus_MC12TuneV4_ProbNNghost", &pi_plus_MC12TuneV4_ProbNNghost, &b_pi_plus_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNmu", &pi_plus_MC15TuneV1_ProbNNmu, &b_pi_plus_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNpi", &pi_plus_MC15TuneV1_ProbNNpi, &b_pi_plus_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNk", &pi_plus_MC15TuneV1_ProbNNk, &b_pi_plus_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNp", &pi_plus_MC15TuneV1_ProbNNp, &b_pi_plus_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("pi_plus_MC15TuneV1_ProbNNghost", &pi_plus_MC15TuneV1_ProbNNghost, &b_pi_plus_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_IP_OWNPV", &pi_plus_IP_OWNPV, &b_pi_plus_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus_IPCHI2_OWNPV", &pi_plus_IPCHI2_OWNPV, &b_pi_plus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus_P", &pi_plus_P, &b_pi_plus_P);
        fChain->SetBranchAddress("pi_plus_PT", &pi_plus_PT, &b_pi_plus_PT);
        fChain->SetBranchAddress("pi_plus_PE", &pi_plus_PE, &b_pi_plus_PE);
        fChain->SetBranchAddress("pi_plus_PX", &pi_plus_PX, &b_pi_plus_PX);
        fChain->SetBranchAddress("pi_plus_PY", &pi_plus_PY, &b_pi_plus_PY);
        fChain->SetBranchAddress("pi_plus_PZ", &pi_plus_PZ, &b_pi_plus_PZ);
        fChain->SetBranchAddress("pi_plus_ID", &pi_plus_ID, &b_pi_plus_ID);
        fChain->SetBranchAddress("pi_plus_PIDmu", &pi_plus_PIDmu, &b_pi_plus_PIDmu);
        fChain->SetBranchAddress("pi_plus_PIDK", &pi_plus_PIDK, &b_pi_plus_PIDK);
        fChain->SetBranchAddress("pi_plus_PIDp", &pi_plus_PIDp, &b_pi_plus_PIDp);
        fChain->SetBranchAddress("pi_plus_ProbNNk", &pi_plus_ProbNNk, &b_pi_plus_ProbNNk);
        fChain->SetBranchAddress("pi_plus_ProbNNp", &pi_plus_ProbNNp, &b_pi_plus_ProbNNp);
        fChain->SetBranchAddress("pi_plus_ProbNNpi", &pi_plus_ProbNNpi, &b_pi_plus_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_ProbNNmu", &pi_plus_ProbNNmu, &b_pi_plus_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_ProbNNghost", &pi_plus_ProbNNghost, &b_pi_plus_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_isMuon", &pi_plus_isMuon, &b_pi_plus_isMuon);
        fChain->SetBranchAddress("pi_plus_TRACK_CHI2NDOF", &pi_plus_TRACK_CHI2NDOF, &b_pi_plus_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus_TRACK_GhostProb", &pi_plus_TRACK_GhostProb, &b_pi_plus_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus_ptasy_1.00", &pi_plus_ptasy_1_00, &b_pi_plus_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNmu", &pi_minus_MC12TuneV4_ProbNNmu, &b_pi_minus_MC12TuneV4_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNpi", &pi_minus_MC12TuneV4_ProbNNpi, &b_pi_minus_MC12TuneV4_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNk", &pi_minus_MC12TuneV4_ProbNNk, &b_pi_minus_MC12TuneV4_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNp", &pi_minus_MC12TuneV4_ProbNNp, &b_pi_minus_MC12TuneV4_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV4_ProbNNghost", &pi_minus_MC12TuneV4_ProbNNghost, &b_pi_minus_MC12TuneV4_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNmu", &pi_minus_MC15TuneV1_ProbNNmu, &b_pi_minus_MC15TuneV1_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNpi", &pi_minus_MC15TuneV1_ProbNNpi, &b_pi_minus_MC15TuneV1_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNk", &pi_minus_MC15TuneV1_ProbNNk, &b_pi_minus_MC15TuneV1_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNp", &pi_minus_MC15TuneV1_ProbNNp, &b_pi_minus_MC15TuneV1_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC15TuneV1_ProbNNghost", &pi_minus_MC15TuneV1_ProbNNghost, &b_pi_minus_MC15TuneV1_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
        fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
        fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
        fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
        fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
        fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
        fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
        fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
        fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
        fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
        fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
        fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
        fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
        fChain->SetBranchAddress("pi_minus_TRACK_CHI2NDOF", &pi_minus_TRACK_CHI2NDOF, &b_pi_minus_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
        fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
        fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
        fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
        fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
        fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
        fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
        fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
        fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
        fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);

    }

    if(_Ds_finalState == Ds_finalState::phipi && _decay == Decay::norm){

	 if(!_data){
	    fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
	    fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
    	    fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);
	  
            fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
    	    fChain->SetBranchAddress("pi_plus1_TRUEID", &pi_plus1_TRUEID);
    	    fChain->SetBranchAddress("pi_plus2_TRUEID", &pi_plus2_TRUEID);
    	    fChain->SetBranchAddress("K_plus_fromDs_TRUEID", &K_plus_fromDs_TRUEID);
    	    fChain->SetBranchAddress("K_minus_fromDs_TRUEID", &K_minus_fromDs_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);

            fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
            fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus1_MC_MOTHER_ID", &pi_plus1_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus2_MC_MOTHER_ID", &pi_plus2_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("K_plus_fromDs_MC_MOTHER_ID", &K_plus_fromDs_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("K_minus_fromDs_MC_MOTHER_ID", &K_minus_fromDs_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);

	    fChain->SetBranchAddress("pi_plus1_PIDK_gen_MagDown", &pi_plus1_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_gen_MagDown", &pi_plus2_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagDown", &pi_minus_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_gen_MagDown", &K_plus_fromDs_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagDown", &K_minus_fromDs_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagDown", &pi_minus_fromDs_PIDK_gen_MagDown);

	    fChain->SetBranchAddress("pi_plus1_PIDK_gen_MagUp", &pi_plus1_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_gen_MagUp", &pi_plus2_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagUp", &pi_minus_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_gen_MagUp", &K_plus_fromDs_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagUp", &K_minus_fromDs_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagUp", &pi_minus_fromDs_PIDK_gen_MagUp);

	    fChain->SetBranchAddress("pi_plus1_PIDK_corr_MagDown", &pi_plus1_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_corr_MagDown", &pi_plus2_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagDown", &pi_minus_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_corr_MagDown", &K_plus_fromDs_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagDown", &K_minus_fromDs_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagDown", &pi_minus_fromDs_PIDK_corr_MagDown);

	    fChain->SetBranchAddress("pi_plus1_PIDK_corr_MagUp", &pi_plus1_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_corr_MagUp", &pi_plus2_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagUp", &pi_minus_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDK_corr_MagUp", &K_plus_fromDs_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagUp", &K_minus_fromDs_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagUp", &pi_minus_fromDs_PIDK_corr_MagUp);

            /*
    	    fChain->SetBranchAddress("pi_plus1_PIDp_gen_MagDown", &pi_plus1_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDp_gen_MagDown", &pi_plus2_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDp_gen_MagDown", &pi_minus_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagDown", &K_plus_fromDs_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagDown", &K_minus_fromDs_PIDp_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagDown", &pi_minus_fromDs_PIDp_gen_MagDown);

	    fChain->SetBranchAddress("pi_plus1_PIDp_gen_MagUp", &pi_plus1_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDp_gen_MagUp", &pi_plus2_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDp_gen_MagUp", &pi_minus_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagUp", &K_plus_fromDs_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagUp", &K_minus_fromDs_PIDp_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagUp", &pi_minus_fromDs_PIDp_gen_MagUp);

	    fChain->SetBranchAddress("pi_plus1_PIDp_corr_MagDown", &pi_plus1_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDp_corr_MagDown", &pi_plus2_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDp_corr_MagDown", &pi_minus_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagDown", &K_plus_fromDs_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagDown", &K_minus_fromDs_PIDp_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagDown", &pi_minus_fromDs_PIDp_corr_MagDown);

	    fChain->SetBranchAddress("pi_plus1_PIDp_corr_MagUp", &pi_plus1_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDp_corr_MagUp", &pi_plus2_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDp_corr_MagUp", &pi_minus_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagUp", &K_plus_fromDs_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagUp", &K_minus_fromDs_PIDp_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagUp", &pi_minus_fromDs_PIDp_corr_MagUp);
	    */
	}

        fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
        fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
        fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
        fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
        fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
        fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
        fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
        fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
        fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
        fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
        fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
        fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
        //fChain->SetBranchAddress("Bs_TAUCHI2", &Bs_TAUCHI2, &b_Bs_TAUCHI2);
        fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
        fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
    /*    fChain->SetBranchAddress("Bs_L0GlobalDecision_TIS", &Bs_L0GlobalDecision_TIS, &b_Bs_L0GlobalDecision_TIS);
        fChain->SetBranchAddress("Bs_L0GlobalDecision_TOS", &Bs_L0GlobalDecision_TOS, &b_Bs_L0GlobalDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TIS", &Bs_Hlt1TrackAllL0Decision_TIS, &b_Bs_Hlt1TrackAllL0Decision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TOS", &Bs_Hlt1TrackAllL0Decision_TOS, &b_Bs_Hlt1TrackAllL0Decision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TIS", &Bs_Hlt1TrackMVADecision_TIS, &b_Bs_Hlt1TrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TOS", &Bs_Hlt1TrackMVADecision_TOS, &b_Bs_Hlt1TrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TIS", &Bs_Hlt1TwoTrackMVADecision_TIS, &b_Bs_Hlt1TwoTrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TOS", &Bs_Hlt1TwoTrackMVADecision_TOS, &b_Bs_Hlt1TwoTrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TIS", &Bs_Hlt1TrackMVALooseDecision_TIS, &b_Bs_Hlt1TrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TOS", &Bs_Hlt1TrackMVALooseDecision_TOS, &b_Bs_Hlt1TrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TIS", &Bs_Hlt1TwoTrackMVALooseDecision_TIS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TOS", &Bs_Hlt1TwoTrackMVALooseDecision_TOS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TIS", &Bs_Hlt2IncPhiDecision_TIS, &b_Bs_Hlt2IncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TOS", &Bs_Hlt2IncPhiDecision_TOS, &b_Bs_Hlt2IncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TIS", &Bs_Hlt2PhiIncPhiDecision_TIS, &b_Bs_Hlt2PhiIncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TOS", &Bs_Hlt2PhiIncPhiDecision_TOS, &b_Bs_Hlt2PhiIncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TIS", &Bs_Hlt2Topo2BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TOS", &Bs_Hlt2Topo2BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TIS", &Bs_Hlt2Topo3BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TOS", &Bs_Hlt2Topo3BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TIS", &Bs_Hlt2Topo4BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TOS", &Bs_Hlt2Topo4BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TIS", &Bs_Hlt2Topo2BodyDecision_TIS, &b_Bs_Hlt2Topo2BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TOS", &Bs_Hlt2Topo2BodyDecision_TOS, &b_Bs_Hlt2Topo2BodyDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TIS", &Bs_Hlt2Topo3BodyDecision_TIS, &b_Bs_Hlt2Topo3BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TOS", &Bs_Hlt2Topo3BodyDecision_TOS, &b_Bs_Hlt2Topo3BodyDecision_TOS);
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
        fChain->SetBranchAddress("Bs_OS_Charm_PROB", &Bs_OS_Charm_PROB, &b_Bs_OS_Charm_PROB);*/
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
        fChain->SetBranchAddress("Bs_B0DTF_M", Bs_B0DTF_M, &b_Bs_B0DTF_M);
        fChain->SetBranchAddress("Bs_B0DTF_MERR", Bs_B0DTF_MERR, &b_Bs_B0DTF_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_P", Bs_B0DTF_P, &b_Bs_B0DTF_P);
        fChain->SetBranchAddress("Bs_B0DTF_PERR", Bs_B0DTF_PERR, &b_Bs_B0DTF_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_PV_X", Bs_B0DTF_PV_X, &b_Bs_B0DTF_PV_X);
        fChain->SetBranchAddress("Bs_B0DTF_PV_Y", Bs_B0DTF_PV_Y, &b_Bs_B0DTF_PV_Y);
        fChain->SetBranchAddress("Bs_B0DTF_PV_Z", Bs_B0DTF_PV_Z, &b_Bs_B0DTF_PV_Z);
        fChain->SetBranchAddress("Bs_B0DTF_PV_key", Bs_B0DTF_PV_key, &b_Bs_B0DTF_PV_key);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_M", Bs_B0DTF_a_1_1260_plus_M, &b_Bs_B0DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_MERR", Bs_B0DTF_a_1_1260_plus_MERR, &b_Bs_B0DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_P", Bs_B0DTF_a_1_1260_plus_P, &b_Bs_B0DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_PERR", Bs_B0DTF_a_1_1260_plus_PERR, &b_Bs_B0DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctau", Bs_B0DTF_a_1_1260_plus_ctau, &b_Bs_B0DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctauErr", Bs_B0DTF_a_1_1260_plus_ctauErr, &b_Bs_B0DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLength", Bs_B0DTF_a_1_1260_plus_decayLength, &b_Bs_B0DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLengthErr", Bs_B0DTF_a_1_1260_plus_decayLengthErr, &b_Bs_B0DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_ID", Bs_B0DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PE", Bs_B0DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PX", Bs_B0DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PY", Bs_B0DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PZ", Bs_B0DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_ID", Bs_B0DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PE", Bs_B0DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PX", Bs_B0DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PY", Bs_B0DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PZ", Bs_B0DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_ID", Bs_B0DTF_a_1_1260_plus_piplus_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PE", Bs_B0DTF_a_1_1260_plus_piplus_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PX", Bs_B0DTF_a_1_1260_plus_piplus_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PY", Bs_B0DTF_a_1_1260_plus_piplus_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PZ", Bs_B0DTF_a_1_1260_plus_piplus_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_PZ);
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
        fChain->SetBranchAddress("Bs_BsDTF_M", Bs_BsDTF_M, &b_Bs_BsDTF_M);
        fChain->SetBranchAddress("Bs_BsDTF_MERR", Bs_BsDTF_MERR, &b_Bs_BsDTF_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_P", Bs_BsDTF_P, &b_Bs_BsDTF_P);
        fChain->SetBranchAddress("Bs_BsDTF_PERR", Bs_BsDTF_PERR, &b_Bs_BsDTF_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_PV_X", Bs_BsDTF_PV_X, &b_Bs_BsDTF_PV_X);
        fChain->SetBranchAddress("Bs_BsDTF_PV_Y", Bs_BsDTF_PV_Y, &b_Bs_BsDTF_PV_Y);
        fChain->SetBranchAddress("Bs_BsDTF_PV_Z", Bs_BsDTF_PV_Z, &b_Bs_BsDTF_PV_Z);
        fChain->SetBranchAddress("Bs_BsDTF_PV_key", Bs_BsDTF_PV_key, &b_Bs_BsDTF_PV_key);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_M", Bs_BsDTF_a_1_1260_plus_M, &b_Bs_BsDTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_MERR", Bs_BsDTF_a_1_1260_plus_MERR, &b_Bs_BsDTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_P", Bs_BsDTF_a_1_1260_plus_P, &b_Bs_BsDTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_PERR", Bs_BsDTF_a_1_1260_plus_PERR, &b_Bs_BsDTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctau", Bs_BsDTF_a_1_1260_plus_ctau, &b_Bs_BsDTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctauErr", Bs_BsDTF_a_1_1260_plus_ctauErr, &b_Bs_BsDTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLength", Bs_BsDTF_a_1_1260_plus_decayLength, &b_Bs_BsDTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLengthErr", Bs_BsDTF_a_1_1260_plus_decayLengthErr, &b_Bs_BsDTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_ID", Bs_BsDTF_a_1_1260_plus_piplus_0_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PE", Bs_BsDTF_a_1_1260_plus_piplus_0_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PX", Bs_BsDTF_a_1_1260_plus_piplus_0_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PY", Bs_BsDTF_a_1_1260_plus_piplus_0_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PZ", Bs_BsDTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_ID", Bs_BsDTF_a_1_1260_plus_piplus_1_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PE", Bs_BsDTF_a_1_1260_plus_piplus_1_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PX", Bs_BsDTF_a_1_1260_plus_piplus_1_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PY", Bs_BsDTF_a_1_1260_plus_piplus_1_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PZ", Bs_BsDTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_ID", Bs_BsDTF_a_1_1260_plus_piplus_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PE", Bs_BsDTF_a_1_1260_plus_piplus_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PX", Bs_BsDTF_a_1_1260_plus_piplus_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PY", Bs_BsDTF_a_1_1260_plus_piplus_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PZ", Bs_BsDTF_a_1_1260_plus_piplus_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_PZ);
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
        fChain->SetBranchAddress("Bs_DTF_M", Bs_DTF_M, &b_Bs_DTF_M);
        fChain->SetBranchAddress("Bs_DTF_MERR", Bs_DTF_MERR, &b_Bs_DTF_MERR);
        fChain->SetBranchAddress("Bs_DTF_P", Bs_DTF_P, &b_Bs_DTF_P);
        fChain->SetBranchAddress("Bs_DTF_PERR", Bs_DTF_PERR, &b_Bs_DTF_PERR);
        fChain->SetBranchAddress("Bs_DTF_PV_X", Bs_DTF_PV_X, &b_Bs_DTF_PV_X);
        fChain->SetBranchAddress("Bs_DTF_PV_Y", Bs_DTF_PV_Y, &b_Bs_DTF_PV_Y);
        fChain->SetBranchAddress("Bs_DTF_PV_Z", Bs_DTF_PV_Z, &b_Bs_DTF_PV_Z);
        fChain->SetBranchAddress("Bs_DTF_PV_key", Bs_DTF_PV_key, &b_Bs_DTF_PV_key);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_M", Bs_DTF_a_1_1260_plus_M, &b_Bs_DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_MERR", Bs_DTF_a_1_1260_plus_MERR, &b_Bs_DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_P", Bs_DTF_a_1_1260_plus_P, &b_Bs_DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_PERR", Bs_DTF_a_1_1260_plus_PERR, &b_Bs_DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctau", Bs_DTF_a_1_1260_plus_ctau, &b_Bs_DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctauErr", Bs_DTF_a_1_1260_plus_ctauErr, &b_Bs_DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLength", Bs_DTF_a_1_1260_plus_decayLength, &b_Bs_DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLengthErr", Bs_DTF_a_1_1260_plus_decayLengthErr, &b_Bs_DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_ID", Bs_DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PE", Bs_DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PX", Bs_DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PY", Bs_DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PZ", Bs_DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_ID", Bs_DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PE", Bs_DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PX", Bs_DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PY", Bs_DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PZ", Bs_DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_ID", Bs_DTF_a_1_1260_plus_piplus_ID, &b_Bs_DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PE", Bs_DTF_a_1_1260_plus_piplus_PE, &b_Bs_DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PX", Bs_DTF_a_1_1260_plus_piplus_PX, &b_Bs_DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PY", Bs_DTF_a_1_1260_plus_piplus_PY, &b_Bs_DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PZ", Bs_DTF_a_1_1260_plus_piplus_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_PZ);
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
        fChain->SetBranchAddress("Bs_PV_M", Bs_PV_M, &b_Bs_PV_M);
        fChain->SetBranchAddress("Bs_PV_MERR", Bs_PV_MERR, &b_Bs_PV_MERR);
        fChain->SetBranchAddress("Bs_PV_P", Bs_PV_P, &b_Bs_PV_P);
        fChain->SetBranchAddress("Bs_PV_PERR", Bs_PV_PERR, &b_Bs_PV_PERR);
        fChain->SetBranchAddress("Bs_PV_PV_X", Bs_PV_PV_X, &b_Bs_PV_PV_X);
        fChain->SetBranchAddress("Bs_PV_PV_Y", Bs_PV_PV_Y, &b_Bs_PV_PV_Y);
        fChain->SetBranchAddress("Bs_PV_PV_Z", Bs_PV_PV_Z, &b_Bs_PV_PV_Z);
        fChain->SetBranchAddress("Bs_PV_PV_key", Bs_PV_PV_key, &b_Bs_PV_PV_key);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_M", Bs_PV_a_1_1260_plus_M, &b_Bs_PV_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_MERR", Bs_PV_a_1_1260_plus_MERR, &b_Bs_PV_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_P", Bs_PV_a_1_1260_plus_P, &b_Bs_PV_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_PERR", Bs_PV_a_1_1260_plus_PERR, &b_Bs_PV_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctau", Bs_PV_a_1_1260_plus_ctau, &b_Bs_PV_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctauErr", Bs_PV_a_1_1260_plus_ctauErr, &b_Bs_PV_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLength", Bs_PV_a_1_1260_plus_decayLength, &b_Bs_PV_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLengthErr", Bs_PV_a_1_1260_plus_decayLengthErr, &b_Bs_PV_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_ID", Bs_PV_a_1_1260_plus_piplus_0_ID, &b_Bs_PV_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PE", Bs_PV_a_1_1260_plus_piplus_0_PE, &b_Bs_PV_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PX", Bs_PV_a_1_1260_plus_piplus_0_PX, &b_Bs_PV_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PY", Bs_PV_a_1_1260_plus_piplus_0_PY, &b_Bs_PV_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PZ", Bs_PV_a_1_1260_plus_piplus_0_PZ, &b_Bs_PV_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_ID", Bs_PV_a_1_1260_plus_piplus_1_ID, &b_Bs_PV_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PE", Bs_PV_a_1_1260_plus_piplus_1_PE, &b_Bs_PV_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PX", Bs_PV_a_1_1260_plus_piplus_1_PX, &b_Bs_PV_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PY", Bs_PV_a_1_1260_plus_piplus_1_PY, &b_Bs_PV_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PZ", Bs_PV_a_1_1260_plus_piplus_1_PZ, &b_Bs_PV_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_ID", Bs_PV_a_1_1260_plus_piplus_ID, &b_Bs_PV_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PE", Bs_PV_a_1_1260_plus_piplus_PE, &b_Bs_PV_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PX", Bs_PV_a_1_1260_plus_piplus_PX, &b_Bs_PV_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PY", Bs_PV_a_1_1260_plus_piplus_PY, &b_Bs_PV_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PZ", Bs_PV_a_1_1260_plus_piplus_PZ, &b_Bs_PV_a_1_1260_plus_piplus_PZ);
        fChain->SetBranchAddress("Bs_PV_chi2", Bs_PV_chi2, &b_Bs_PV_chi2);
        fChain->SetBranchAddress("Bs_PV_ctau", Bs_PV_ctau, &b_Bs_PV_ctau);
        fChain->SetBranchAddress("Bs_PV_ctauErr", Bs_PV_ctauErr, &b_Bs_PV_ctauErr);
        fChain->SetBranchAddress("Bs_PV_decayLength", Bs_PV_decayLength, &b_Bs_PV_decayLength);
        fChain->SetBranchAddress("Bs_PV_decayLengthErr", Bs_PV_decayLengthErr, &b_Bs_PV_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_nDOF", Bs_PV_nDOF, &b_Bs_PV_nDOF);
        fChain->SetBranchAddress("Bs_PV_nIter", Bs_PV_nIter, &b_Bs_PV_nIter);
        fChain->SetBranchAddress("Bs_PV_status", Bs_PV_status, &b_Bs_PV_status);
/*        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION_OS", &Bs_BsTaggingTool_TAGDECISION_OS, &b_Bs_BsTaggingTool_TAGDECISION_OS);
        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA_OS", &Bs_BsTaggingTool_TAGOMEGA_OS, &b_Bs_BsTaggingTool_TAGOMEGA_OS);
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
        fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Charm_PROB", &Bs_BsTaggingTool_OS_Charm_PROB, &b_Bs_BsTaggingTool_OS_Charm_PROB);*/
        fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
        fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
        fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
        fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
        fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
        fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
        fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
        fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
        fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
        fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
        fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
        fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
        fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Ds_ID", &Ds_ID, &b_Ds_ID);
/*        fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
        fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);*/
        //fChain->SetBranchAddress("Ds_TAUCHI2", &Ds_TAUCHI2, &b_Ds_TAUCHI2);
        fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
        fChain->SetBranchAddress("K_plus_fromDs_ETA", &K_plus_fromDs_ETA, &b_K_plus_fromDs_ETA);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNmu", &K_plus_fromDs_MC12TuneV2_ProbNNmu, &b_K_plus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNpi", &K_plus_fromDs_MC12TuneV2_ProbNNpi, &b_K_plus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNk", &K_plus_fromDs_MC12TuneV2_ProbNNk, &b_K_plus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNp", &K_plus_fromDs_MC12TuneV2_ProbNNp, &b_K_plus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV2_ProbNNghost", &K_plus_fromDs_MC12TuneV2_ProbNNghost, &b_K_plus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNmu", &K_plus_fromDs_MC12TuneV3_ProbNNmu, &b_K_plus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNpi", &K_plus_fromDs_MC12TuneV3_ProbNNpi, &b_K_plus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNk", &K_plus_fromDs_MC12TuneV3_ProbNNk, &b_K_plus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNp", &K_plus_fromDs_MC12TuneV3_ProbNNp, &b_K_plus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("K_plus_fromDs_MC12TuneV3_ProbNNghost", &K_plus_fromDs_MC12TuneV3_ProbNNghost, &b_K_plus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("K_plus_fromDs_IP_OWNPV", &K_plus_fromDs_IP_OWNPV, &b_K_plus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("K_plus_fromDs_IPCHI2_OWNPV", &K_plus_fromDs_IPCHI2_OWNPV, &b_K_plus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("K_plus_fromDs_P", &K_plus_fromDs_P, &b_K_plus_fromDs_P);
        fChain->SetBranchAddress("K_plus_fromDs_PT", &K_plus_fromDs_PT, &b_K_plus_fromDs_PT);
        fChain->SetBranchAddress("K_plus_fromDs_PE", &K_plus_fromDs_PE, &b_K_plus_fromDs_PE);
        fChain->SetBranchAddress("K_plus_fromDs_PX", &K_plus_fromDs_PX, &b_K_plus_fromDs_PX);
        fChain->SetBranchAddress("K_plus_fromDs_PY", &K_plus_fromDs_PY, &b_K_plus_fromDs_PY);
        fChain->SetBranchAddress("K_plus_fromDs_PZ", &K_plus_fromDs_PZ, &b_K_plus_fromDs_PZ);
        fChain->SetBranchAddress("K_plus_fromDs_ID", &K_plus_fromDs_ID, &b_K_plus_fromDs_ID);
        fChain->SetBranchAddress("K_plus_fromDs_PIDmu", &K_plus_fromDs_PIDmu, &b_K_plus_fromDs_PIDmu);
        fChain->SetBranchAddress("K_plus_fromDs_PIDK", &K_plus_fromDs_PIDK, &b_K_plus_fromDs_PIDK);
        fChain->SetBranchAddress("K_plus_fromDs_PIDp", &K_plus_fromDs_PIDp, &b_K_plus_fromDs_PIDp);
        fChain->SetBranchAddress("K_plus_fromDs_ProbNNk", &K_plus_fromDs_ProbNNk, &b_K_plus_fromDs_ProbNNk);
        fChain->SetBranchAddress("K_plus_fromDs_ProbNNp", &K_plus_fromDs_ProbNNp, &b_K_plus_fromDs_ProbNNp);
        fChain->SetBranchAddress("K_plus_fromDs_ProbNNpi", &K_plus_fromDs_ProbNNpi, &b_K_plus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("K_plus_fromDs_ProbNNmu", &K_plus_fromDs_ProbNNmu, &b_K_plus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("K_plus_fromDs_ProbNNghost", &K_plus_fromDs_ProbNNghost, &b_K_plus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("K_plus_fromDs_isMuon", &K_plus_fromDs_isMuon, &b_K_plus_fromDs_isMuon);
        fChain->SetBranchAddress("K_plus_fromDs_TRACK_CHI2NDOF", &K_plus_fromDs_TRACK_CHI2NDOF, &b_K_plus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("K_plus_fromDs_TRACK_GhostProb", &K_plus_fromDs_TRACK_GhostProb, &b_K_plus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("K_plus_fromDs_ptasy_1.00", &K_plus_fromDs_ptasy_1_00, &b_K_plus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("K_minus_fromDs_ETA", &K_minus_fromDs_ETA, &b_K_minus_fromDs_ETA);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNmu", &K_minus_fromDs_MC12TuneV2_ProbNNmu, &b_K_minus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNpi", &K_minus_fromDs_MC12TuneV2_ProbNNpi, &b_K_minus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNk", &K_minus_fromDs_MC12TuneV2_ProbNNk, &b_K_minus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNp", &K_minus_fromDs_MC12TuneV2_ProbNNp, &b_K_minus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNghost", &K_minus_fromDs_MC12TuneV2_ProbNNghost, &b_K_minus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNmu", &K_minus_fromDs_MC12TuneV3_ProbNNmu, &b_K_minus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNpi", &K_minus_fromDs_MC12TuneV3_ProbNNpi, &b_K_minus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNk", &K_minus_fromDs_MC12TuneV3_ProbNNk, &b_K_minus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNp", &K_minus_fromDs_MC12TuneV3_ProbNNp, &b_K_minus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNghost", &K_minus_fromDs_MC12TuneV3_ProbNNghost, &b_K_minus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("K_minus_fromDs_IP_OWNPV", &K_minus_fromDs_IP_OWNPV, &b_K_minus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("K_minus_fromDs_IPCHI2_OWNPV", &K_minus_fromDs_IPCHI2_OWNPV, &b_K_minus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("K_minus_fromDs_P", &K_minus_fromDs_P, &b_K_minus_fromDs_P);
        fChain->SetBranchAddress("K_minus_fromDs_PT", &K_minus_fromDs_PT, &b_K_minus_fromDs_PT);
        fChain->SetBranchAddress("K_minus_fromDs_PE", &K_minus_fromDs_PE, &b_K_minus_fromDs_PE);
        fChain->SetBranchAddress("K_minus_fromDs_PX", &K_minus_fromDs_PX, &b_K_minus_fromDs_PX);
        fChain->SetBranchAddress("K_minus_fromDs_PY", &K_minus_fromDs_PY, &b_K_minus_fromDs_PY);
        fChain->SetBranchAddress("K_minus_fromDs_PZ", &K_minus_fromDs_PZ, &b_K_minus_fromDs_PZ);
        fChain->SetBranchAddress("K_minus_fromDs_ID", &K_minus_fromDs_ID, &b_K_minus_fromDs_ID);
        fChain->SetBranchAddress("K_minus_fromDs_PIDmu", &K_minus_fromDs_PIDmu, &b_K_minus_fromDs_PIDmu);
        fChain->SetBranchAddress("K_minus_fromDs_PIDK", &K_minus_fromDs_PIDK, &b_K_minus_fromDs_PIDK);
        fChain->SetBranchAddress("K_minus_fromDs_PIDp", &K_minus_fromDs_PIDp, &b_K_minus_fromDs_PIDp);
        fChain->SetBranchAddress("K_minus_fromDs_ProbNNk", &K_minus_fromDs_ProbNNk, &b_K_minus_fromDs_ProbNNk);
        fChain->SetBranchAddress("K_minus_fromDs_ProbNNp", &K_minus_fromDs_ProbNNp, &b_K_minus_fromDs_ProbNNp);
        fChain->SetBranchAddress("K_minus_fromDs_ProbNNpi", &K_minus_fromDs_ProbNNpi, &b_K_minus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("K_minus_fromDs_ProbNNmu", &K_minus_fromDs_ProbNNmu, &b_K_minus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("K_minus_fromDs_ProbNNghost", &K_minus_fromDs_ProbNNghost, &b_K_minus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("K_minus_fromDs_isMuon", &K_minus_fromDs_isMuon, &b_K_minus_fromDs_isMuon);
        fChain->SetBranchAddress("K_minus_fromDs_TRACK_CHI2NDOF", &K_minus_fromDs_TRACK_CHI2NDOF, &b_K_minus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb, &b_K_minus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("K_minus_fromDs_ptasy_1.00", &K_minus_fromDs_ptasy_1_00, &b_K_minus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
        fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
        fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
        fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
        fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
        fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
        fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CHI2NDOF", &pi_minus_fromDs_TRACK_CHI2NDOF, &b_pi_minus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_fromDs_ptasy_1.00", &pi_minus_fromDs_ptasy_1_00, &b_pi_minus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA1", &a_1_1260_plus_DOCA1, &b_a_1_1260_plus_DOCA1);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA2", &a_1_1260_plus_DOCA2, &b_a_1_1260_plus_DOCA2);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA3", &a_1_1260_plus_DOCA3, &b_a_1_1260_plus_DOCA3);
        fChain->SetBranchAddress("a_1_1260_plus_ETA", &a_1_1260_plus_ETA, &b_a_1_1260_plus_ETA);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_X", &a_1_1260_plus_ENDVERTEX_X, &b_a_1_1260_plus_ENDVERTEX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Y", &a_1_1260_plus_ENDVERTEX_Y, &b_a_1_1260_plus_ENDVERTEX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Z", &a_1_1260_plus_ENDVERTEX_Z, &b_a_1_1260_plus_ENDVERTEX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_XERR", &a_1_1260_plus_ENDVERTEX_XERR, &b_a_1_1260_plus_ENDVERTEX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_YERR", &a_1_1260_plus_ENDVERTEX_YERR, &b_a_1_1260_plus_ENDVERTEX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_ZERR", &a_1_1260_plus_ENDVERTEX_ZERR, &b_a_1_1260_plus_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_CHI2", &a_1_1260_plus_ENDVERTEX_CHI2, &b_a_1_1260_plus_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_NDOF", &a_1_1260_plus_ENDVERTEX_NDOF, &b_a_1_1260_plus_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_X", &a_1_1260_plus_OWNPV_X, &b_a_1_1260_plus_OWNPV_X);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Y", &a_1_1260_plus_OWNPV_Y, &b_a_1_1260_plus_OWNPV_Y);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Z", &a_1_1260_plus_OWNPV_Z, &b_a_1_1260_plus_OWNPV_Z);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_XERR", &a_1_1260_plus_OWNPV_XERR, &b_a_1_1260_plus_OWNPV_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_YERR", &a_1_1260_plus_OWNPV_YERR, &b_a_1_1260_plus_OWNPV_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_ZERR", &a_1_1260_plus_OWNPV_ZERR, &b_a_1_1260_plus_OWNPV_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_CHI2", &a_1_1260_plus_OWNPV_CHI2, &b_a_1_1260_plus_OWNPV_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_NDOF", &a_1_1260_plus_OWNPV_NDOF, &b_a_1_1260_plus_OWNPV_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_IP_OWNPV", &a_1_1260_plus_IP_OWNPV, &b_a_1_1260_plus_IP_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_IPCHI2_OWNPV", &a_1_1260_plus_IPCHI2_OWNPV, &b_a_1_1260_plus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FD_OWNPV", &a_1_1260_plus_FD_OWNPV, &b_a_1_1260_plus_FD_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_OWNPV", &a_1_1260_plus_FDCHI2_OWNPV, &b_a_1_1260_plus_FDCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_OWNPV", &a_1_1260_plus_DIRA_OWNPV, &b_a_1_1260_plus_DIRA_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_X", &a_1_1260_plus_ORIVX_X, &b_a_1_1260_plus_ORIVX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Y", &a_1_1260_plus_ORIVX_Y, &b_a_1_1260_plus_ORIVX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Z", &a_1_1260_plus_ORIVX_Z, &b_a_1_1260_plus_ORIVX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_XERR", &a_1_1260_plus_ORIVX_XERR, &b_a_1_1260_plus_ORIVX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_YERR", &a_1_1260_plus_ORIVX_YERR, &b_a_1_1260_plus_ORIVX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_ZERR", &a_1_1260_plus_ORIVX_ZERR, &b_a_1_1260_plus_ORIVX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_CHI2", &a_1_1260_plus_ORIVX_CHI2, &b_a_1_1260_plus_ORIVX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_NDOF", &a_1_1260_plus_ORIVX_NDOF, &b_a_1_1260_plus_ORIVX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_FD_ORIVX", &a_1_1260_plus_FD_ORIVX, &b_a_1_1260_plus_FD_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_ORIVX", &a_1_1260_plus_FDCHI2_ORIVX, &b_a_1_1260_plus_FDCHI2_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_ORIVX", &a_1_1260_plus_DIRA_ORIVX, &b_a_1_1260_plus_DIRA_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_P", &a_1_1260_plus_P, &b_a_1_1260_plus_P);
        fChain->SetBranchAddress("a_1_1260_plus_PT", &a_1_1260_plus_PT, &b_a_1_1260_plus_PT);
        fChain->SetBranchAddress("a_1_1260_plus_PE", &a_1_1260_plus_PE, &b_a_1_1260_plus_PE);
        fChain->SetBranchAddress("a_1_1260_plus_PX", &a_1_1260_plus_PX, &b_a_1_1260_plus_PX);
        fChain->SetBranchAddress("a_1_1260_plus_PY", &a_1_1260_plus_PY, &b_a_1_1260_plus_PY);
        fChain->SetBranchAddress("a_1_1260_plus_PZ", &a_1_1260_plus_PZ, &b_a_1_1260_plus_PZ);
        fChain->SetBranchAddress("a_1_1260_plus_MM", &a_1_1260_plus_MM, &b_a_1_1260_plus_MM);
        fChain->SetBranchAddress("a_1_1260_plus_MMERR", &a_1_1260_plus_MMERR, &b_a_1_1260_plus_MMERR);
        fChain->SetBranchAddress("a_1_1260_plus_ID", &a_1_1260_plus_ID, &b_a_1_1260_plus_ID);
/*        fChain->SetBranchAddress("a_1_1260_plus_TAU", &a_1_1260_plus_TAU, &b_a_1_1260_plus_TAU);
        fChain->SetBranchAddress("a_1_1260_plus_TAUERR", &a_1_1260_plus_TAUERR, &b_a_1_1260_plus_TAUERR);*/
        //fChain->SetBranchAddress("a_1_1260_plus_TAUCHI2", &a_1_1260_plus_TAUCHI2, &b_a_1_1260_plus_TAUCHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ptasy_1.00", &a_1_1260_plus_ptasy_1_00, &b_a_1_1260_plus_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus1_ETA", &pi_plus1_ETA, &b_pi_plus1_ETA);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNmu", &pi_plus1_MC12TuneV2_ProbNNmu, &b_pi_plus1_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNpi", &pi_plus1_MC12TuneV2_ProbNNpi, &b_pi_plus1_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNk", &pi_plus1_MC12TuneV2_ProbNNk, &b_pi_plus1_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNp", &pi_plus1_MC12TuneV2_ProbNNp, &b_pi_plus1_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNghost", &pi_plus1_MC12TuneV2_ProbNNghost, &b_pi_plus1_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNmu", &pi_plus1_MC12TuneV3_ProbNNmu, &b_pi_plus1_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNpi", &pi_plus1_MC12TuneV3_ProbNNpi, &b_pi_plus1_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNk", &pi_plus1_MC12TuneV3_ProbNNk, &b_pi_plus1_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNp", &pi_plus1_MC12TuneV3_ProbNNp, &b_pi_plus1_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNghost", &pi_plus1_MC12TuneV3_ProbNNghost, &b_pi_plus1_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_IP_OWNPV", &pi_plus1_IP_OWNPV, &b_pi_plus1_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus1_IPCHI2_OWNPV", &pi_plus1_IPCHI2_OWNPV, &b_pi_plus1_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus1_P", &pi_plus1_P, &b_pi_plus1_P);
        fChain->SetBranchAddress("pi_plus1_PT", &pi_plus1_PT, &b_pi_plus1_PT);
        fChain->SetBranchAddress("pi_plus1_PE", &pi_plus1_PE, &b_pi_plus1_PE);
        fChain->SetBranchAddress("pi_plus1_PX", &pi_plus1_PX, &b_pi_plus1_PX);
        fChain->SetBranchAddress("pi_plus1_PY", &pi_plus1_PY, &b_pi_plus1_PY);
        fChain->SetBranchAddress("pi_plus1_PZ", &pi_plus1_PZ, &b_pi_plus1_PZ);
        fChain->SetBranchAddress("pi_plus1_ID", &pi_plus1_ID, &b_pi_plus1_ID);
        fChain->SetBranchAddress("pi_plus1_PIDmu", &pi_plus1_PIDmu, &b_pi_plus1_PIDmu);
        fChain->SetBranchAddress("pi_plus1_PIDK", &pi_plus1_PIDK, &b_pi_plus1_PIDK);
        fChain->SetBranchAddress("pi_plus1_PIDp", &pi_plus1_PIDp, &b_pi_plus1_PIDp);
        fChain->SetBranchAddress("pi_plus1_ProbNNk", &pi_plus1_ProbNNk, &b_pi_plus1_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_ProbNNp", &pi_plus1_ProbNNp, &b_pi_plus1_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_ProbNNpi", &pi_plus1_ProbNNpi, &b_pi_plus1_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_ProbNNmu", &pi_plus1_ProbNNmu, &b_pi_plus1_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_ProbNNghost", &pi_plus1_ProbNNghost, &b_pi_plus1_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_isMuon", &pi_plus1_isMuon, &b_pi_plus1_isMuon);
        fChain->SetBranchAddress("pi_plus1_TRACK_CHI2NDOF", &pi_plus1_TRACK_CHI2NDOF, &b_pi_plus1_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus1_TRACK_GhostProb", &pi_plus1_TRACK_GhostProb, &b_pi_plus1_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus1_ptasy_1.00", &pi_plus1_ptasy_1_00, &b_pi_plus1_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus2_ETA", &pi_plus2_ETA, &b_pi_plus2_ETA);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNmu", &pi_plus2_MC12TuneV2_ProbNNmu, &b_pi_plus2_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNpi", &pi_plus2_MC12TuneV2_ProbNNpi, &b_pi_plus2_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNk", &pi_plus2_MC12TuneV2_ProbNNk, &b_pi_plus2_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNp", &pi_plus2_MC12TuneV2_ProbNNp, &b_pi_plus2_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNghost", &pi_plus2_MC12TuneV2_ProbNNghost, &b_pi_plus2_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNmu", &pi_plus2_MC12TuneV3_ProbNNmu, &b_pi_plus2_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNpi", &pi_plus2_MC12TuneV3_ProbNNpi, &b_pi_plus2_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNk", &pi_plus2_MC12TuneV3_ProbNNk, &b_pi_plus2_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNp", &pi_plus2_MC12TuneV3_ProbNNp, &b_pi_plus2_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNghost", &pi_plus2_MC12TuneV3_ProbNNghost, &b_pi_plus2_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_IP_OWNPV", &pi_plus2_IP_OWNPV, &b_pi_plus2_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus2_IPCHI2_OWNPV", &pi_plus2_IPCHI2_OWNPV, &b_pi_plus2_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus2_P", &pi_plus2_P, &b_pi_plus2_P);
        fChain->SetBranchAddress("pi_plus2_PT", &pi_plus2_PT, &b_pi_plus2_PT);
        fChain->SetBranchAddress("pi_plus2_PE", &pi_plus2_PE, &b_pi_plus2_PE);
        fChain->SetBranchAddress("pi_plus2_PX", &pi_plus2_PX, &b_pi_plus2_PX);
        fChain->SetBranchAddress("pi_plus2_PY", &pi_plus2_PY, &b_pi_plus2_PY);
        fChain->SetBranchAddress("pi_plus2_PZ", &pi_plus2_PZ, &b_pi_plus2_PZ);
        fChain->SetBranchAddress("pi_plus2_ID", &pi_plus2_ID, &b_pi_plus2_ID);
        fChain->SetBranchAddress("pi_plus2_PIDmu", &pi_plus2_PIDmu, &b_pi_plus2_PIDmu);
        fChain->SetBranchAddress("pi_plus2_PIDK", &pi_plus2_PIDK, &b_pi_plus2_PIDK);
        fChain->SetBranchAddress("pi_plus2_PIDp", &pi_plus2_PIDp, &b_pi_plus2_PIDp);
        fChain->SetBranchAddress("pi_plus2_ProbNNk", &pi_plus2_ProbNNk, &b_pi_plus2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_ProbNNp", &pi_plus2_ProbNNp, &b_pi_plus2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_ProbNNpi", &pi_plus2_ProbNNpi, &b_pi_plus2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_ProbNNmu", &pi_plus2_ProbNNmu, &b_pi_plus2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_ProbNNghost", &pi_plus2_ProbNNghost, &b_pi_plus2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_isMuon", &pi_plus2_isMuon, &b_pi_plus2_isMuon);
        fChain->SetBranchAddress("pi_plus2_TRACK_CHI2NDOF", &pi_plus2_TRACK_CHI2NDOF, &b_pi_plus2_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus2_TRACK_GhostProb", &pi_plus2_TRACK_GhostProb, &b_pi_plus2_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus2_ptasy_1.00", &pi_plus2_ptasy_1_00, &b_pi_plus2_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
        fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
        fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
        fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
        fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
        fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
        fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
        fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
        fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
        fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
        fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
        fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
        fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
        fChain->SetBranchAddress("pi_minus_TRACK_CHI2NDOF", &pi_minus_TRACK_CHI2NDOF, &b_pi_minus_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
        fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
        fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
        fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
        fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
        fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
        fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
        fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
        fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
        fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);

    }

    if(_Ds_finalState == Ds_finalState::pipipi && _decay == Decay::norm){
	
 	if(!_data){
	    fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
	    fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
    	    fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);
	
	    fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
    	    fChain->SetBranchAddress("pi_plus1_TRUEID", &pi_plus1_TRUEID);
    	    fChain->SetBranchAddress("pi_plus2_TRUEID", &pi_plus2_TRUEID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);
	    fChain->SetBranchAddress("pi_plus_fromDs_TRUEID", &pi_plus_fromDs_TRUEID);
   	    fChain->SetBranchAddress("pi_minus2_fromDs_TRUEID", &pi_minus2_fromDs_TRUEID);

	    fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
	    fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus1_MC_MOTHER_ID", &pi_plus1_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus2_MC_MOTHER_ID", &pi_plus2_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);
	    fChain->SetBranchAddress("pi_plus_fromDs_MC_MOTHER_ID", &pi_plus_fromDs_MC_MOTHER_ID);
   	    fChain->SetBranchAddress("pi_minus2_fromDs_MC_MOTHER_ID", &pi_minus2_fromDs_MC_MOTHER_ID);
	}

        fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
        fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
        fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
        fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
        fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
        fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
        fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
        fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
        fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
        fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
        fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
        fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
        fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
        //fChain->SetBranchAddress("Bs_TAUCHI2", &Bs_TAUCHI2, &b_Bs_TAUCHI2);
        fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
        fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
        fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
        fChain->SetBranchAddress("Bs_L0GlobalDecision_TIS", &Bs_L0GlobalDecision_TIS, &b_Bs_L0GlobalDecision_TIS);
        fChain->SetBranchAddress("Bs_L0GlobalDecision_TOS", &Bs_L0GlobalDecision_TOS, &b_Bs_L0GlobalDecision_TOS);
/*        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TIS", &Bs_Hlt1TrackAllL0Decision_TIS, &b_Bs_Hlt1TrackAllL0Decision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackAllL0Decision_TOS", &Bs_Hlt1TrackAllL0Decision_TOS, &b_Bs_Hlt1TrackAllL0Decision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TIS", &Bs_Hlt1TrackMVADecision_TIS, &b_Bs_Hlt1TrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVADecision_TOS", &Bs_Hlt1TrackMVADecision_TOS, &b_Bs_Hlt1TrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TIS", &Bs_Hlt1TwoTrackMVADecision_TIS, &b_Bs_Hlt1TwoTrackMVADecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVADecision_TOS", &Bs_Hlt1TwoTrackMVADecision_TOS, &b_Bs_Hlt1TwoTrackMVADecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TIS", &Bs_Hlt1TrackMVALooseDecision_TIS, &b_Bs_Hlt1TrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TrackMVALooseDecision_TOS", &Bs_Hlt1TrackMVALooseDecision_TOS, &b_Bs_Hlt1TrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TIS", &Bs_Hlt1TwoTrackMVALooseDecision_TIS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt1TwoTrackMVALooseDecision_TOS", &Bs_Hlt1TwoTrackMVALooseDecision_TOS, &b_Bs_Hlt1TwoTrackMVALooseDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TIS", &Bs_Hlt2IncPhiDecision_TIS, &b_Bs_Hlt2IncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2IncPhiDecision_TOS", &Bs_Hlt2IncPhiDecision_TOS, &b_Bs_Hlt2IncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TIS", &Bs_Hlt2PhiIncPhiDecision_TIS, &b_Bs_Hlt2PhiIncPhiDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2PhiIncPhiDecision_TOS", &Bs_Hlt2PhiIncPhiDecision_TOS, &b_Bs_Hlt2PhiIncPhiDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TIS", &Bs_Hlt2Topo2BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyBBDTDecision_TOS", &Bs_Hlt2Topo2BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo2BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TIS", &Bs_Hlt2Topo3BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyBBDTDecision_TOS", &Bs_Hlt2Topo3BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo3BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TIS", &Bs_Hlt2Topo4BodyBBDTDecision_TIS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo4BodyBBDTDecision_TOS", &Bs_Hlt2Topo4BodyBBDTDecision_TOS, &b_Bs_Hlt2Topo4BodyBBDTDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TIS", &Bs_Hlt2Topo2BodyDecision_TIS, &b_Bs_Hlt2Topo2BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo2BodyDecision_TOS", &Bs_Hlt2Topo2BodyDecision_TOS, &b_Bs_Hlt2Topo2BodyDecision_TOS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TIS", &Bs_Hlt2Topo3BodyDecision_TIS, &b_Bs_Hlt2Topo3BodyDecision_TIS);
        fChain->SetBranchAddress("Bs_Hlt2Topo3BodyDecision_TOS", &Bs_Hlt2Topo3BodyDecision_TOS, &b_Bs_Hlt2Topo3BodyDecision_TOS);
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
        fChain->SetBranchAddress("Bs_OS_Charm_PROB", &Bs_OS_Charm_PROB, &b_Bs_OS_Charm_PROB);*/
        fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
        fChain->SetBranchAddress("Bs_B0DTF_nPV", &Bs_B0DTF_nPV, &b_Bs_B0DTF_nPV);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_M", Bs_B0DTF_D_splus_M, &b_Bs_B0DTF_D_splus_M);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_MERR", Bs_B0DTF_D_splus_MERR, &b_Bs_B0DTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_P", Bs_B0DTF_D_splus_P, &b_Bs_B0DTF_D_splus_P);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_PERR", Bs_B0DTF_D_splus_PERR, &b_Bs_B0DTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctau", Bs_B0DTF_D_splus_ctau, &b_Bs_B0DTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_ctauErr", Bs_B0DTF_D_splus_ctauErr, &b_Bs_B0DTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLength", Bs_B0DTF_D_splus_decayLength, &b_Bs_B0DTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_decayLengthErr", Bs_B0DTF_D_splus_decayLengthErr, &b_Bs_B0DTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_ID", Bs_B0DTF_D_splus_piplus_0_ID, &b_Bs_B0DTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PE", Bs_B0DTF_D_splus_piplus_0_PE, &b_Bs_B0DTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PX", Bs_B0DTF_D_splus_piplus_0_PX, &b_Bs_B0DTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PY", Bs_B0DTF_D_splus_piplus_0_PY, &b_Bs_B0DTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PZ", Bs_B0DTF_D_splus_piplus_0_PZ, &b_Bs_B0DTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_ID", Bs_B0DTF_D_splus_piplus_1_ID, &b_Bs_B0DTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PE", Bs_B0DTF_D_splus_piplus_1_PE, &b_Bs_B0DTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PX", Bs_B0DTF_D_splus_piplus_1_PX, &b_Bs_B0DTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PY", Bs_B0DTF_D_splus_piplus_1_PY, &b_Bs_B0DTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_1_PZ", Bs_B0DTF_D_splus_piplus_1_PZ, &b_Bs_B0DTF_D_splus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_ID", Bs_B0DTF_D_splus_piplus_ID, &b_Bs_B0DTF_D_splus_piplus_ID);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PE", Bs_B0DTF_D_splus_piplus_PE, &b_Bs_B0DTF_D_splus_piplus_PE);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PX", Bs_B0DTF_D_splus_piplus_PX, &b_Bs_B0DTF_D_splus_piplus_PX);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PY", Bs_B0DTF_D_splus_piplus_PY, &b_Bs_B0DTF_D_splus_piplus_PY);
        fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PZ", Bs_B0DTF_D_splus_piplus_PZ, &b_Bs_B0DTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_M", Bs_B0DTF_M, &b_Bs_B0DTF_M);
        fChain->SetBranchAddress("Bs_B0DTF_MERR", Bs_B0DTF_MERR, &b_Bs_B0DTF_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_P", Bs_B0DTF_P, &b_Bs_B0DTF_P);
        fChain->SetBranchAddress("Bs_B0DTF_PERR", Bs_B0DTF_PERR, &b_Bs_B0DTF_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_PV_X", Bs_B0DTF_PV_X, &b_Bs_B0DTF_PV_X);
        fChain->SetBranchAddress("Bs_B0DTF_PV_Y", Bs_B0DTF_PV_Y, &b_Bs_B0DTF_PV_Y);
        fChain->SetBranchAddress("Bs_B0DTF_PV_Z", Bs_B0DTF_PV_Z, &b_Bs_B0DTF_PV_Z);
        fChain->SetBranchAddress("Bs_B0DTF_PV_key", Bs_B0DTF_PV_key, &b_Bs_B0DTF_PV_key);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_M", Bs_B0DTF_a_1_1260_plus_M, &b_Bs_B0DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_MERR", Bs_B0DTF_a_1_1260_plus_MERR, &b_Bs_B0DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_P", Bs_B0DTF_a_1_1260_plus_P, &b_Bs_B0DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_PERR", Bs_B0DTF_a_1_1260_plus_PERR, &b_Bs_B0DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctau", Bs_B0DTF_a_1_1260_plus_ctau, &b_Bs_B0DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctauErr", Bs_B0DTF_a_1_1260_plus_ctauErr, &b_Bs_B0DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLength", Bs_B0DTF_a_1_1260_plus_decayLength, &b_Bs_B0DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLengthErr", Bs_B0DTF_a_1_1260_plus_decayLengthErr, &b_Bs_B0DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_ID", Bs_B0DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PE", Bs_B0DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PX", Bs_B0DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PY", Bs_B0DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PZ", Bs_B0DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_ID", Bs_B0DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PE", Bs_B0DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PX", Bs_B0DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PY", Bs_B0DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PZ", Bs_B0DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_ID", Bs_B0DTF_a_1_1260_plus_piplus_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PE", Bs_B0DTF_a_1_1260_plus_piplus_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PX", Bs_B0DTF_a_1_1260_plus_piplus_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PY", Bs_B0DTF_a_1_1260_plus_piplus_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PZ", Bs_B0DTF_a_1_1260_plus_piplus_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_chi2", Bs_B0DTF_chi2, &b_Bs_B0DTF_chi2);
        fChain->SetBranchAddress("Bs_B0DTF_ctau", Bs_B0DTF_ctau, &b_Bs_B0DTF_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_ctauErr", Bs_B0DTF_ctauErr, &b_Bs_B0DTF_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_decayLength", Bs_B0DTF_decayLength, &b_Bs_B0DTF_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_decayLengthErr", Bs_B0DTF_decayLengthErr, &b_Bs_B0DTF_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_nDOF", Bs_B0DTF_nDOF, &b_Bs_B0DTF_nDOF);
        fChain->SetBranchAddress("Bs_B0DTF_nIter", Bs_B0DTF_nIter, &b_Bs_B0DTF_nIter);
        fChain->SetBranchAddress("Bs_B0DTF_status", Bs_B0DTF_status, &b_Bs_B0DTF_status);
        fChain->SetBranchAddress("Bs_BsDTF_nPV", &Bs_BsDTF_nPV, &b_Bs_BsDTF_nPV);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_M", Bs_BsDTF_D_splus_M, &b_Bs_BsDTF_D_splus_M);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_MERR", Bs_BsDTF_D_splus_MERR, &b_Bs_BsDTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_P", Bs_BsDTF_D_splus_P, &b_Bs_BsDTF_D_splus_P);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_PERR", Bs_BsDTF_D_splus_PERR, &b_Bs_BsDTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctau", Bs_BsDTF_D_splus_ctau, &b_Bs_BsDTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_ctauErr", Bs_BsDTF_D_splus_ctauErr, &b_Bs_BsDTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLength", Bs_BsDTF_D_splus_decayLength, &b_Bs_BsDTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_decayLengthErr", Bs_BsDTF_D_splus_decayLengthErr, &b_Bs_BsDTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_ID", Bs_BsDTF_D_splus_piplus_0_ID, &b_Bs_BsDTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PE", Bs_BsDTF_D_splus_piplus_0_PE, &b_Bs_BsDTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PX", Bs_BsDTF_D_splus_piplus_0_PX, &b_Bs_BsDTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PY", Bs_BsDTF_D_splus_piplus_0_PY, &b_Bs_BsDTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PZ", Bs_BsDTF_D_splus_piplus_0_PZ, &b_Bs_BsDTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_ID", Bs_BsDTF_D_splus_piplus_1_ID, &b_Bs_BsDTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PE", Bs_BsDTF_D_splus_piplus_1_PE, &b_Bs_BsDTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PX", Bs_BsDTF_D_splus_piplus_1_PX, &b_Bs_BsDTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PY", Bs_BsDTF_D_splus_piplus_1_PY, &b_Bs_BsDTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_1_PZ", Bs_BsDTF_D_splus_piplus_1_PZ, &b_Bs_BsDTF_D_splus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_ID", Bs_BsDTF_D_splus_piplus_ID, &b_Bs_BsDTF_D_splus_piplus_ID);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PE", Bs_BsDTF_D_splus_piplus_PE, &b_Bs_BsDTF_D_splus_piplus_PE);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PX", Bs_BsDTF_D_splus_piplus_PX, &b_Bs_BsDTF_D_splus_piplus_PX);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PY", Bs_BsDTF_D_splus_piplus_PY, &b_Bs_BsDTF_D_splus_piplus_PY);
        fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PZ", Bs_BsDTF_D_splus_piplus_PZ, &b_Bs_BsDTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_M", Bs_BsDTF_M, &b_Bs_BsDTF_M);
        fChain->SetBranchAddress("Bs_BsDTF_MERR", Bs_BsDTF_MERR, &b_Bs_BsDTF_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_P", Bs_BsDTF_P, &b_Bs_BsDTF_P);
        fChain->SetBranchAddress("Bs_BsDTF_PERR", Bs_BsDTF_PERR, &b_Bs_BsDTF_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_PV_X", Bs_BsDTF_PV_X, &b_Bs_BsDTF_PV_X);
        fChain->SetBranchAddress("Bs_BsDTF_PV_Y", Bs_BsDTF_PV_Y, &b_Bs_BsDTF_PV_Y);
        fChain->SetBranchAddress("Bs_BsDTF_PV_Z", Bs_BsDTF_PV_Z, &b_Bs_BsDTF_PV_Z);
        fChain->SetBranchAddress("Bs_BsDTF_PV_key", Bs_BsDTF_PV_key, &b_Bs_BsDTF_PV_key);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_M", Bs_BsDTF_a_1_1260_plus_M, &b_Bs_BsDTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_MERR", Bs_BsDTF_a_1_1260_plus_MERR, &b_Bs_BsDTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_P", Bs_BsDTF_a_1_1260_plus_P, &b_Bs_BsDTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_PERR", Bs_BsDTF_a_1_1260_plus_PERR, &b_Bs_BsDTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctau", Bs_BsDTF_a_1_1260_plus_ctau, &b_Bs_BsDTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctauErr", Bs_BsDTF_a_1_1260_plus_ctauErr, &b_Bs_BsDTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLength", Bs_BsDTF_a_1_1260_plus_decayLength, &b_Bs_BsDTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLengthErr", Bs_BsDTF_a_1_1260_plus_decayLengthErr, &b_Bs_BsDTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_ID", Bs_BsDTF_a_1_1260_plus_piplus_0_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PE", Bs_BsDTF_a_1_1260_plus_piplus_0_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PX", Bs_BsDTF_a_1_1260_plus_piplus_0_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PY", Bs_BsDTF_a_1_1260_plus_piplus_0_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PZ", Bs_BsDTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_ID", Bs_BsDTF_a_1_1260_plus_piplus_1_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PE", Bs_BsDTF_a_1_1260_plus_piplus_1_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PX", Bs_BsDTF_a_1_1260_plus_piplus_1_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PY", Bs_BsDTF_a_1_1260_plus_piplus_1_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PZ", Bs_BsDTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_ID", Bs_BsDTF_a_1_1260_plus_piplus_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PE", Bs_BsDTF_a_1_1260_plus_piplus_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PX", Bs_BsDTF_a_1_1260_plus_piplus_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PY", Bs_BsDTF_a_1_1260_plus_piplus_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PZ", Bs_BsDTF_a_1_1260_plus_piplus_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_chi2", Bs_BsDTF_chi2, &b_Bs_BsDTF_chi2);
        fChain->SetBranchAddress("Bs_BsDTF_ctau", Bs_BsDTF_ctau, &b_Bs_BsDTF_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_ctauErr", Bs_BsDTF_ctauErr, &b_Bs_BsDTF_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_decayLength", Bs_BsDTF_decayLength, &b_Bs_BsDTF_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_decayLengthErr", Bs_BsDTF_decayLengthErr, &b_Bs_BsDTF_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_nDOF", Bs_BsDTF_nDOF, &b_Bs_BsDTF_nDOF);
        fChain->SetBranchAddress("Bs_BsDTF_nIter", Bs_BsDTF_nIter, &b_Bs_BsDTF_nIter);
        fChain->SetBranchAddress("Bs_BsDTF_status", Bs_BsDTF_status, &b_Bs_BsDTF_status);
        fChain->SetBranchAddress("Bs_DTF_nPV", &Bs_DTF_nPV, &b_Bs_DTF_nPV);
        fChain->SetBranchAddress("Bs_DTF_D_splus_M", Bs_DTF_D_splus_M, &b_Bs_DTF_D_splus_M);
        fChain->SetBranchAddress("Bs_DTF_D_splus_MERR", Bs_DTF_D_splus_MERR, &b_Bs_DTF_D_splus_MERR);
        fChain->SetBranchAddress("Bs_DTF_D_splus_P", Bs_DTF_D_splus_P, &b_Bs_DTF_D_splus_P);
        fChain->SetBranchAddress("Bs_DTF_D_splus_PERR", Bs_DTF_D_splus_PERR, &b_Bs_DTF_D_splus_PERR);
        fChain->SetBranchAddress("Bs_DTF_D_splus_ctau", Bs_DTF_D_splus_ctau, &b_Bs_DTF_D_splus_ctau);
        fChain->SetBranchAddress("Bs_DTF_D_splus_ctauErr", Bs_DTF_D_splus_ctauErr, &b_Bs_DTF_D_splus_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_D_splus_decayLength", Bs_DTF_D_splus_decayLength, &b_Bs_DTF_D_splus_decayLength);
        fChain->SetBranchAddress("Bs_DTF_D_splus_decayLengthErr", Bs_DTF_D_splus_decayLengthErr, &b_Bs_DTF_D_splus_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_ID", Bs_DTF_D_splus_piplus_0_ID, &b_Bs_DTF_D_splus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PE", Bs_DTF_D_splus_piplus_0_PE, &b_Bs_DTF_D_splus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PX", Bs_DTF_D_splus_piplus_0_PX, &b_Bs_DTF_D_splus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PY", Bs_DTF_D_splus_piplus_0_PY, &b_Bs_DTF_D_splus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PZ", Bs_DTF_D_splus_piplus_0_PZ, &b_Bs_DTF_D_splus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_ID", Bs_DTF_D_splus_piplus_1_ID, &b_Bs_DTF_D_splus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PE", Bs_DTF_D_splus_piplus_1_PE, &b_Bs_DTF_D_splus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PX", Bs_DTF_D_splus_piplus_1_PX, &b_Bs_DTF_D_splus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PY", Bs_DTF_D_splus_piplus_1_PY, &b_Bs_DTF_D_splus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_1_PZ", Bs_DTF_D_splus_piplus_1_PZ, &b_Bs_DTF_D_splus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_ID", Bs_DTF_D_splus_piplus_ID, &b_Bs_DTF_D_splus_piplus_ID);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PE", Bs_DTF_D_splus_piplus_PE, &b_Bs_DTF_D_splus_piplus_PE);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PX", Bs_DTF_D_splus_piplus_PX, &b_Bs_DTF_D_splus_piplus_PX);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PY", Bs_DTF_D_splus_piplus_PY, &b_Bs_DTF_D_splus_piplus_PY);
        fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PZ", Bs_DTF_D_splus_piplus_PZ, &b_Bs_DTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_DTF_M", Bs_DTF_M, &b_Bs_DTF_M);
        fChain->SetBranchAddress("Bs_DTF_MERR", Bs_DTF_MERR, &b_Bs_DTF_MERR);
        fChain->SetBranchAddress("Bs_DTF_P", Bs_DTF_P, &b_Bs_DTF_P);
        fChain->SetBranchAddress("Bs_DTF_PERR", Bs_DTF_PERR, &b_Bs_DTF_PERR);
        fChain->SetBranchAddress("Bs_DTF_PV_X", Bs_DTF_PV_X, &b_Bs_DTF_PV_X);
        fChain->SetBranchAddress("Bs_DTF_PV_Y", Bs_DTF_PV_Y, &b_Bs_DTF_PV_Y);
        fChain->SetBranchAddress("Bs_DTF_PV_Z", Bs_DTF_PV_Z, &b_Bs_DTF_PV_Z);
        fChain->SetBranchAddress("Bs_DTF_PV_key", Bs_DTF_PV_key, &b_Bs_DTF_PV_key);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_M", Bs_DTF_a_1_1260_plus_M, &b_Bs_DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_MERR", Bs_DTF_a_1_1260_plus_MERR, &b_Bs_DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_P", Bs_DTF_a_1_1260_plus_P, &b_Bs_DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_PERR", Bs_DTF_a_1_1260_plus_PERR, &b_Bs_DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctau", Bs_DTF_a_1_1260_plus_ctau, &b_Bs_DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctauErr", Bs_DTF_a_1_1260_plus_ctauErr, &b_Bs_DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLength", Bs_DTF_a_1_1260_plus_decayLength, &b_Bs_DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLengthErr", Bs_DTF_a_1_1260_plus_decayLengthErr, &b_Bs_DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_ID", Bs_DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PE", Bs_DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PX", Bs_DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PY", Bs_DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PZ", Bs_DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_ID", Bs_DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PE", Bs_DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PX", Bs_DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PY", Bs_DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PZ", Bs_DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_ID", Bs_DTF_a_1_1260_plus_piplus_ID, &b_Bs_DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PE", Bs_DTF_a_1_1260_plus_piplus_PE, &b_Bs_DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PX", Bs_DTF_a_1_1260_plus_piplus_PX, &b_Bs_DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PY", Bs_DTF_a_1_1260_plus_piplus_PY, &b_Bs_DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PZ", Bs_DTF_a_1_1260_plus_piplus_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_PZ);
        fChain->SetBranchAddress("Bs_DTF_chi2", Bs_DTF_chi2, &b_Bs_DTF_chi2);
        fChain->SetBranchAddress("Bs_DTF_ctau", Bs_DTF_ctau, &b_Bs_DTF_ctau);
        fChain->SetBranchAddress("Bs_DTF_ctauErr", Bs_DTF_ctauErr, &b_Bs_DTF_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_decayLength", Bs_DTF_decayLength, &b_Bs_DTF_decayLength);
        fChain->SetBranchAddress("Bs_DTF_decayLengthErr", Bs_DTF_decayLengthErr, &b_Bs_DTF_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_nDOF", Bs_DTF_nDOF, &b_Bs_DTF_nDOF);
        fChain->SetBranchAddress("Bs_DTF_nIter", Bs_DTF_nIter, &b_Bs_DTF_nIter);
        fChain->SetBranchAddress("Bs_DTF_status", Bs_DTF_status, &b_Bs_DTF_status);
        fChain->SetBranchAddress("Bs_PV_nPV", &Bs_PV_nPV, &b_Bs_PV_nPV);
        fChain->SetBranchAddress("Bs_PV_Dplus_M", Bs_PV_Dplus_M, &b_Bs_PV_Dplus_M);
        fChain->SetBranchAddress("Bs_PV_Dplus_MERR", Bs_PV_Dplus_MERR, &b_Bs_PV_Dplus_MERR);
        fChain->SetBranchAddress("Bs_PV_Dplus_P", Bs_PV_Dplus_P, &b_Bs_PV_Dplus_P);
        fChain->SetBranchAddress("Bs_PV_Dplus_PERR", Bs_PV_Dplus_PERR, &b_Bs_PV_Dplus_PERR);
        fChain->SetBranchAddress("Bs_PV_Dplus_ctau", Bs_PV_Dplus_ctau, &b_Bs_PV_Dplus_ctau);
        fChain->SetBranchAddress("Bs_PV_Dplus_ctauErr", Bs_PV_Dplus_ctauErr, &b_Bs_PV_Dplus_ctauErr);
        fChain->SetBranchAddress("Bs_PV_Dplus_decayLength", Bs_PV_Dplus_decayLength, &b_Bs_PV_Dplus_decayLength);
        fChain->SetBranchAddress("Bs_PV_Dplus_decayLengthErr", Bs_PV_Dplus_decayLengthErr, &b_Bs_PV_Dplus_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_ID", Bs_PV_Dplus_piplus_0_ID, &b_Bs_PV_Dplus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PE", Bs_PV_Dplus_piplus_0_PE, &b_Bs_PV_Dplus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PX", Bs_PV_Dplus_piplus_0_PX, &b_Bs_PV_Dplus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PY", Bs_PV_Dplus_piplus_0_PY, &b_Bs_PV_Dplus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PZ", Bs_PV_Dplus_piplus_0_PZ, &b_Bs_PV_Dplus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_ID", Bs_PV_Dplus_piplus_1_ID, &b_Bs_PV_Dplus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PE", Bs_PV_Dplus_piplus_1_PE, &b_Bs_PV_Dplus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PX", Bs_PV_Dplus_piplus_1_PX, &b_Bs_PV_Dplus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PY", Bs_PV_Dplus_piplus_1_PY, &b_Bs_PV_Dplus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_1_PZ", Bs_PV_Dplus_piplus_1_PZ, &b_Bs_PV_Dplus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_ID", Bs_PV_Dplus_piplus_ID, &b_Bs_PV_Dplus_piplus_ID);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PE", Bs_PV_Dplus_piplus_PE, &b_Bs_PV_Dplus_piplus_PE);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PX", Bs_PV_Dplus_piplus_PX, &b_Bs_PV_Dplus_piplus_PX);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PY", Bs_PV_Dplus_piplus_PY, &b_Bs_PV_Dplus_piplus_PY);
        fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PZ", Bs_PV_Dplus_piplus_PZ, &b_Bs_PV_Dplus_piplus_PZ);
        fChain->SetBranchAddress("Bs_PV_M", Bs_PV_M, &b_Bs_PV_M);
        fChain->SetBranchAddress("Bs_PV_MERR", Bs_PV_MERR, &b_Bs_PV_MERR);
        fChain->SetBranchAddress("Bs_PV_P", Bs_PV_P, &b_Bs_PV_P);
        fChain->SetBranchAddress("Bs_PV_PERR", Bs_PV_PERR, &b_Bs_PV_PERR);
        fChain->SetBranchAddress("Bs_PV_PV_X", Bs_PV_PV_X, &b_Bs_PV_PV_X);
        fChain->SetBranchAddress("Bs_PV_PV_Y", Bs_PV_PV_Y, &b_Bs_PV_PV_Y);
        fChain->SetBranchAddress("Bs_PV_PV_Z", Bs_PV_PV_Z, &b_Bs_PV_PV_Z);
        fChain->SetBranchAddress("Bs_PV_PV_key", Bs_PV_PV_key, &b_Bs_PV_PV_key);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_M", Bs_PV_a_1_1260_plus_M, &b_Bs_PV_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_MERR", Bs_PV_a_1_1260_plus_MERR, &b_Bs_PV_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_P", Bs_PV_a_1_1260_plus_P, &b_Bs_PV_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_PERR", Bs_PV_a_1_1260_plus_PERR, &b_Bs_PV_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctau", Bs_PV_a_1_1260_plus_ctau, &b_Bs_PV_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctauErr", Bs_PV_a_1_1260_plus_ctauErr, &b_Bs_PV_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLength", Bs_PV_a_1_1260_plus_decayLength, &b_Bs_PV_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLengthErr", Bs_PV_a_1_1260_plus_decayLengthErr, &b_Bs_PV_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_ID", Bs_PV_a_1_1260_plus_piplus_0_ID, &b_Bs_PV_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PE", Bs_PV_a_1_1260_plus_piplus_0_PE, &b_Bs_PV_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PX", Bs_PV_a_1_1260_plus_piplus_0_PX, &b_Bs_PV_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PY", Bs_PV_a_1_1260_plus_piplus_0_PY, &b_Bs_PV_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PZ", Bs_PV_a_1_1260_plus_piplus_0_PZ, &b_Bs_PV_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_ID", Bs_PV_a_1_1260_plus_piplus_1_ID, &b_Bs_PV_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PE", Bs_PV_a_1_1260_plus_piplus_1_PE, &b_Bs_PV_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PX", Bs_PV_a_1_1260_plus_piplus_1_PX, &b_Bs_PV_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PY", Bs_PV_a_1_1260_plus_piplus_1_PY, &b_Bs_PV_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PZ", Bs_PV_a_1_1260_plus_piplus_1_PZ, &b_Bs_PV_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_ID", Bs_PV_a_1_1260_plus_piplus_ID, &b_Bs_PV_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PE", Bs_PV_a_1_1260_plus_piplus_PE, &b_Bs_PV_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PX", Bs_PV_a_1_1260_plus_piplus_PX, &b_Bs_PV_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PY", Bs_PV_a_1_1260_plus_piplus_PY, &b_Bs_PV_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PZ", Bs_PV_a_1_1260_plus_piplus_PZ, &b_Bs_PV_a_1_1260_plus_piplus_PZ);
        fChain->SetBranchAddress("Bs_PV_chi2", Bs_PV_chi2, &b_Bs_PV_chi2);
        fChain->SetBranchAddress("Bs_PV_ctau", Bs_PV_ctau, &b_Bs_PV_ctau);
        fChain->SetBranchAddress("Bs_PV_ctauErr", Bs_PV_ctauErr, &b_Bs_PV_ctauErr);
        fChain->SetBranchAddress("Bs_PV_decayLength", Bs_PV_decayLength, &b_Bs_PV_decayLength);
        fChain->SetBranchAddress("Bs_PV_decayLengthErr", Bs_PV_decayLengthErr, &b_Bs_PV_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_nDOF", Bs_PV_nDOF, &b_Bs_PV_nDOF);
        fChain->SetBranchAddress("Bs_PV_nIter", Bs_PV_nIter, &b_Bs_PV_nIter);
        fChain->SetBranchAddress("Bs_PV_status", Bs_PV_status, &b_Bs_PV_status);
/*        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGDECISION_OS", &Bs_BsTaggingTool_TAGDECISION_OS, &b_Bs_BsTaggingTool_TAGDECISION_OS);
        fChain->SetBranchAddress("Bs_BsTaggingTool_TAGOMEGA_OS", &Bs_BsTaggingTool_TAGOMEGA_OS, &b_Bs_BsTaggingTool_TAGOMEGA_OS);
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
        fChain->SetBranchAddress("Bs_BsTaggingTool_OS_Charm_PROB", &Bs_BsTaggingTool_OS_Charm_PROB, &b_Bs_BsTaggingTool_OS_Charm_PROB);*/
        fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
        fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
        fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
        fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
        fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
        fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
        fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
        fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
        fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
        fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
        fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
        fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
        fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
        fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
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
        fChain->SetBranchAddress("Ds_ID", &Ds_ID, &b_Ds_ID);
/*        fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
        fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);*/
        //fChain->SetBranchAddress("Ds_TAUCHI2", &Ds_TAUCHI2, &b_Ds_TAUCHI2);
        fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus_fromDs_ETA", &pi_plus_fromDs_ETA, &b_pi_plus_fromDs_ETA);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNmu", &pi_plus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNpi", &pi_plus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNk", &pi_plus_fromDs_MC12TuneV2_ProbNNk, &b_pi_plus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNp", &pi_plus_fromDs_MC12TuneV2_ProbNNp, &b_pi_plus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNghost", &pi_plus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNmu", &pi_plus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNpi", &pi_plus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNk", &pi_plus_fromDs_MC12TuneV3_ProbNNk, &b_pi_plus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNp", &pi_plus_fromDs_MC12TuneV3_ProbNNp, &b_pi_plus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNghost", &pi_plus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_IP_OWNPV", &pi_plus_fromDs_IP_OWNPV, &b_pi_plus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus_fromDs_IPCHI2_OWNPV", &pi_plus_fromDs_IPCHI2_OWNPV, &b_pi_plus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus_fromDs_P", &pi_plus_fromDs_P, &b_pi_plus_fromDs_P);
        fChain->SetBranchAddress("pi_plus_fromDs_PT", &pi_plus_fromDs_PT, &b_pi_plus_fromDs_PT);
        fChain->SetBranchAddress("pi_plus_fromDs_PE", &pi_plus_fromDs_PE, &b_pi_plus_fromDs_PE);
        fChain->SetBranchAddress("pi_plus_fromDs_PX", &pi_plus_fromDs_PX, &b_pi_plus_fromDs_PX);
        fChain->SetBranchAddress("pi_plus_fromDs_PY", &pi_plus_fromDs_PY, &b_pi_plus_fromDs_PY);
        fChain->SetBranchAddress("pi_plus_fromDs_PZ", &pi_plus_fromDs_PZ, &b_pi_plus_fromDs_PZ);
        fChain->SetBranchAddress("pi_plus_fromDs_ID", &pi_plus_fromDs_ID, &b_pi_plus_fromDs_ID);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDmu", &pi_plus_fromDs_PIDmu, &b_pi_plus_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDK", &pi_plus_fromDs_PIDK, &b_pi_plus_fromDs_PIDK);
        fChain->SetBranchAddress("pi_plus_fromDs_PIDp", &pi_plus_fromDs_PIDp, &b_pi_plus_fromDs_PIDp);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNk", &pi_plus_fromDs_ProbNNk, &b_pi_plus_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNp", &pi_plus_fromDs_ProbNNp, &b_pi_plus_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNpi", &pi_plus_fromDs_ProbNNpi, &b_pi_plus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNmu", &pi_plus_fromDs_ProbNNmu, &b_pi_plus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_plus_fromDs_ProbNNghost", &pi_plus_fromDs_ProbNNghost, &b_pi_plus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_plus_fromDs_isMuon", &pi_plus_fromDs_isMuon, &b_pi_plus_fromDs_isMuon);
        fChain->SetBranchAddress("pi_plus_fromDs_TRACK_CHI2NDOF", &pi_plus_fromDs_TRACK_CHI2NDOF, &b_pi_plus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus_fromDs_TRACK_GhostProb", &pi_plus_fromDs_TRACK_GhostProb, &b_pi_plus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus_fromDs_ptasy_1.00", &pi_plus_fromDs_ptasy_1_00, &b_pi_plus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
        fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
        fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
        fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
        fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
        fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
        fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
        fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_CHI2NDOF", &pi_minus_fromDs_TRACK_CHI2NDOF, &b_pi_minus_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_fromDs_ptasy_1.00", &pi_minus_fromDs_ptasy_1_00, &b_pi_minus_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus2_fromDs_ETA", &pi_minus2_fromDs_ETA, &b_pi_minus2_fromDs_ETA);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNmu", &pi_minus2_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNpi", &pi_minus2_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNk", &pi_minus2_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNp", &pi_minus2_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV2_ProbNNghost", &pi_minus2_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus2_fromDs_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNmu", &pi_minus2_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNpi", &pi_minus2_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNk", &pi_minus2_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNp", &pi_minus2_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_MC12TuneV3_ProbNNghost", &pi_minus2_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus2_fromDs_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_IP_OWNPV", &pi_minus2_fromDs_IP_OWNPV, &b_pi_minus2_fromDs_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus2_fromDs_IPCHI2_OWNPV", &pi_minus2_fromDs_IPCHI2_OWNPV, &b_pi_minus2_fromDs_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus2_fromDs_P", &pi_minus2_fromDs_P, &b_pi_minus2_fromDs_P);
        fChain->SetBranchAddress("pi_minus2_fromDs_PT", &pi_minus2_fromDs_PT, &b_pi_minus2_fromDs_PT);
        fChain->SetBranchAddress("pi_minus2_fromDs_PE", &pi_minus2_fromDs_PE, &b_pi_minus2_fromDs_PE);
        fChain->SetBranchAddress("pi_minus2_fromDs_PX", &pi_minus2_fromDs_PX, &b_pi_minus2_fromDs_PX);
        fChain->SetBranchAddress("pi_minus2_fromDs_PY", &pi_minus2_fromDs_PY, &b_pi_minus2_fromDs_PY);
        fChain->SetBranchAddress("pi_minus2_fromDs_PZ", &pi_minus2_fromDs_PZ, &b_pi_minus2_fromDs_PZ);
        fChain->SetBranchAddress("pi_minus2_fromDs_ID", &pi_minus2_fromDs_ID, &b_pi_minus2_fromDs_ID);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDmu", &pi_minus2_fromDs_PIDmu, &b_pi_minus2_fromDs_PIDmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDK", &pi_minus2_fromDs_PIDK, &b_pi_minus2_fromDs_PIDK);
        fChain->SetBranchAddress("pi_minus2_fromDs_PIDp", &pi_minus2_fromDs_PIDp, &b_pi_minus2_fromDs_PIDp);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNk", &pi_minus2_fromDs_ProbNNk, &b_pi_minus2_fromDs_ProbNNk);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNp", &pi_minus2_fromDs_ProbNNp, &b_pi_minus2_fromDs_ProbNNp);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNpi", &pi_minus2_fromDs_ProbNNpi, &b_pi_minus2_fromDs_ProbNNpi);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNmu", &pi_minus2_fromDs_ProbNNmu, &b_pi_minus2_fromDs_ProbNNmu);
        fChain->SetBranchAddress("pi_minus2_fromDs_ProbNNghost", &pi_minus2_fromDs_ProbNNghost, &b_pi_minus2_fromDs_ProbNNghost);
        fChain->SetBranchAddress("pi_minus2_fromDs_isMuon", &pi_minus2_fromDs_isMuon, &b_pi_minus2_fromDs_isMuon);
        fChain->SetBranchAddress("pi_minus2_fromDs_TRACK_CHI2NDOF", &pi_minus2_fromDs_TRACK_CHI2NDOF, &b_pi_minus2_fromDs_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus2_fromDs_TRACK_GhostProb", &pi_minus2_fromDs_TRACK_GhostProb, &b_pi_minus2_fromDs_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus2_fromDs_ptasy_1.00", &pi_minus2_fromDs_ptasy_1_00, &b_pi_minus2_fromDs_ptasy_1_00);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA1", &a_1_1260_plus_DOCA1, &b_a_1_1260_plus_DOCA1);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA2", &a_1_1260_plus_DOCA2, &b_a_1_1260_plus_DOCA2);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA3", &a_1_1260_plus_DOCA3, &b_a_1_1260_plus_DOCA3);
        fChain->SetBranchAddress("a_1_1260_plus_ETA", &a_1_1260_plus_ETA, &b_a_1_1260_plus_ETA);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_X", &a_1_1260_plus_ENDVERTEX_X, &b_a_1_1260_plus_ENDVERTEX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Y", &a_1_1260_plus_ENDVERTEX_Y, &b_a_1_1260_plus_ENDVERTEX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Z", &a_1_1260_plus_ENDVERTEX_Z, &b_a_1_1260_plus_ENDVERTEX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_XERR", &a_1_1260_plus_ENDVERTEX_XERR, &b_a_1_1260_plus_ENDVERTEX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_YERR", &a_1_1260_plus_ENDVERTEX_YERR, &b_a_1_1260_plus_ENDVERTEX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_ZERR", &a_1_1260_plus_ENDVERTEX_ZERR, &b_a_1_1260_plus_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_CHI2", &a_1_1260_plus_ENDVERTEX_CHI2, &b_a_1_1260_plus_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_NDOF", &a_1_1260_plus_ENDVERTEX_NDOF, &b_a_1_1260_plus_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_X", &a_1_1260_plus_OWNPV_X, &b_a_1_1260_plus_OWNPV_X);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Y", &a_1_1260_plus_OWNPV_Y, &b_a_1_1260_plus_OWNPV_Y);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Z", &a_1_1260_plus_OWNPV_Z, &b_a_1_1260_plus_OWNPV_Z);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_XERR", &a_1_1260_plus_OWNPV_XERR, &b_a_1_1260_plus_OWNPV_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_YERR", &a_1_1260_plus_OWNPV_YERR, &b_a_1_1260_plus_OWNPV_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_ZERR", &a_1_1260_plus_OWNPV_ZERR, &b_a_1_1260_plus_OWNPV_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_CHI2", &a_1_1260_plus_OWNPV_CHI2, &b_a_1_1260_plus_OWNPV_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_NDOF", &a_1_1260_plus_OWNPV_NDOF, &b_a_1_1260_plus_OWNPV_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_IP_OWNPV", &a_1_1260_plus_IP_OWNPV, &b_a_1_1260_plus_IP_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_IPCHI2_OWNPV", &a_1_1260_plus_IPCHI2_OWNPV, &b_a_1_1260_plus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FD_OWNPV", &a_1_1260_plus_FD_OWNPV, &b_a_1_1260_plus_FD_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_OWNPV", &a_1_1260_plus_FDCHI2_OWNPV, &b_a_1_1260_plus_FDCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_OWNPV", &a_1_1260_plus_DIRA_OWNPV, &b_a_1_1260_plus_DIRA_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_X", &a_1_1260_plus_ORIVX_X, &b_a_1_1260_plus_ORIVX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Y", &a_1_1260_plus_ORIVX_Y, &b_a_1_1260_plus_ORIVX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Z", &a_1_1260_plus_ORIVX_Z, &b_a_1_1260_plus_ORIVX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_XERR", &a_1_1260_plus_ORIVX_XERR, &b_a_1_1260_plus_ORIVX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_YERR", &a_1_1260_plus_ORIVX_YERR, &b_a_1_1260_plus_ORIVX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_ZERR", &a_1_1260_plus_ORIVX_ZERR, &b_a_1_1260_plus_ORIVX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_CHI2", &a_1_1260_plus_ORIVX_CHI2, &b_a_1_1260_plus_ORIVX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_NDOF", &a_1_1260_plus_ORIVX_NDOF, &b_a_1_1260_plus_ORIVX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_FD_ORIVX", &a_1_1260_plus_FD_ORIVX, &b_a_1_1260_plus_FD_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_ORIVX", &a_1_1260_plus_FDCHI2_ORIVX, &b_a_1_1260_plus_FDCHI2_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_ORIVX", &a_1_1260_plus_DIRA_ORIVX, &b_a_1_1260_plus_DIRA_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_P", &a_1_1260_plus_P, &b_a_1_1260_plus_P);
        fChain->SetBranchAddress("a_1_1260_plus_PT", &a_1_1260_plus_PT, &b_a_1_1260_plus_PT);
        fChain->SetBranchAddress("a_1_1260_plus_PE", &a_1_1260_plus_PE, &b_a_1_1260_plus_PE);
        fChain->SetBranchAddress("a_1_1260_plus_PX", &a_1_1260_plus_PX, &b_a_1_1260_plus_PX);
        fChain->SetBranchAddress("a_1_1260_plus_PY", &a_1_1260_plus_PY, &b_a_1_1260_plus_PY);
        fChain->SetBranchAddress("a_1_1260_plus_PZ", &a_1_1260_plus_PZ, &b_a_1_1260_plus_PZ);
        fChain->SetBranchAddress("a_1_1260_plus_MM", &a_1_1260_plus_MM, &b_a_1_1260_plus_MM);
        fChain->SetBranchAddress("a_1_1260_plus_MMERR", &a_1_1260_plus_MMERR, &b_a_1_1260_plus_MMERR);
        fChain->SetBranchAddress("a_1_1260_plus_ID", &a_1_1260_plus_ID, &b_a_1_1260_plus_ID);
  /*      fChain->SetBranchAddress("a_1_1260_plus_TAU", &a_1_1260_plus_TAU, &b_a_1_1260_plus_TAU);
        fChain->SetBranchAddress("a_1_1260_plus_TAUERR", &a_1_1260_plus_TAUERR, &b_a_1_1260_plus_TAUERR);
  */      //fChain->SetBranchAddress("a_1_1260_plus_TAUCHI2", &a_1_1260_plus_TAUCHI2, &b_a_1_1260_plus_TAUCHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ptasy_1.00", &a_1_1260_plus_ptasy_1_00, &b_a_1_1260_plus_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus1_ETA", &pi_plus1_ETA, &b_pi_plus1_ETA);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNmu", &pi_plus1_MC12TuneV2_ProbNNmu, &b_pi_plus1_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNpi", &pi_plus1_MC12TuneV2_ProbNNpi, &b_pi_plus1_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNk", &pi_plus1_MC12TuneV2_ProbNNk, &b_pi_plus1_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNp", &pi_plus1_MC12TuneV2_ProbNNp, &b_pi_plus1_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNghost", &pi_plus1_MC12TuneV2_ProbNNghost, &b_pi_plus1_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNmu", &pi_plus1_MC12TuneV3_ProbNNmu, &b_pi_plus1_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNpi", &pi_plus1_MC12TuneV3_ProbNNpi, &b_pi_plus1_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNk", &pi_plus1_MC12TuneV3_ProbNNk, &b_pi_plus1_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNp", &pi_plus1_MC12TuneV3_ProbNNp, &b_pi_plus1_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNghost", &pi_plus1_MC12TuneV3_ProbNNghost, &b_pi_plus1_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_IP_OWNPV", &pi_plus1_IP_OWNPV, &b_pi_plus1_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus1_IPCHI2_OWNPV", &pi_plus1_IPCHI2_OWNPV, &b_pi_plus1_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus1_P", &pi_plus1_P, &b_pi_plus1_P);
        fChain->SetBranchAddress("pi_plus1_PT", &pi_plus1_PT, &b_pi_plus1_PT);
        fChain->SetBranchAddress("pi_plus1_PE", &pi_plus1_PE, &b_pi_plus1_PE);
        fChain->SetBranchAddress("pi_plus1_PX", &pi_plus1_PX, &b_pi_plus1_PX);
        fChain->SetBranchAddress("pi_plus1_PY", &pi_plus1_PY, &b_pi_plus1_PY);
        fChain->SetBranchAddress("pi_plus1_PZ", &pi_plus1_PZ, &b_pi_plus1_PZ);
        fChain->SetBranchAddress("pi_plus1_ID", &pi_plus1_ID, &b_pi_plus1_ID);
        fChain->SetBranchAddress("pi_plus1_PIDmu", &pi_plus1_PIDmu, &b_pi_plus1_PIDmu);
        fChain->SetBranchAddress("pi_plus1_PIDK", &pi_plus1_PIDK, &b_pi_plus1_PIDK);
        fChain->SetBranchAddress("pi_plus1_PIDp", &pi_plus1_PIDp, &b_pi_plus1_PIDp);
        fChain->SetBranchAddress("pi_plus1_ProbNNk", &pi_plus1_ProbNNk, &b_pi_plus1_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_ProbNNp", &pi_plus1_ProbNNp, &b_pi_plus1_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_ProbNNpi", &pi_plus1_ProbNNpi, &b_pi_plus1_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_ProbNNmu", &pi_plus1_ProbNNmu, &b_pi_plus1_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_ProbNNghost", &pi_plus1_ProbNNghost, &b_pi_plus1_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_isMuon", &pi_plus1_isMuon, &b_pi_plus1_isMuon);
        fChain->SetBranchAddress("pi_plus1_TRACK_CHI2NDOF", &pi_plus1_TRACK_CHI2NDOF, &b_pi_plus1_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus1_TRACK_GhostProb", &pi_plus1_TRACK_GhostProb, &b_pi_plus1_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus1_ptasy_1.00", &pi_plus1_ptasy_1_00, &b_pi_plus1_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus2_ETA", &pi_plus2_ETA, &b_pi_plus2_ETA);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNmu", &pi_plus2_MC12TuneV2_ProbNNmu, &b_pi_plus2_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNpi", &pi_plus2_MC12TuneV2_ProbNNpi, &b_pi_plus2_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNk", &pi_plus2_MC12TuneV2_ProbNNk, &b_pi_plus2_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNp", &pi_plus2_MC12TuneV2_ProbNNp, &b_pi_plus2_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNghost", &pi_plus2_MC12TuneV2_ProbNNghost, &b_pi_plus2_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNmu", &pi_plus2_MC12TuneV3_ProbNNmu, &b_pi_plus2_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNpi", &pi_plus2_MC12TuneV3_ProbNNpi, &b_pi_plus2_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNk", &pi_plus2_MC12TuneV3_ProbNNk, &b_pi_plus2_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNp", &pi_plus2_MC12TuneV3_ProbNNp, &b_pi_plus2_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNghost", &pi_plus2_MC12TuneV3_ProbNNghost, &b_pi_plus2_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_IP_OWNPV", &pi_plus2_IP_OWNPV, &b_pi_plus2_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus2_IPCHI2_OWNPV", &pi_plus2_IPCHI2_OWNPV, &b_pi_plus2_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus2_P", &pi_plus2_P, &b_pi_plus2_P);
        fChain->SetBranchAddress("pi_plus2_PT", &pi_plus2_PT, &b_pi_plus2_PT);
        fChain->SetBranchAddress("pi_plus2_PE", &pi_plus2_PE, &b_pi_plus2_PE);
        fChain->SetBranchAddress("pi_plus2_PX", &pi_plus2_PX, &b_pi_plus2_PX);
        fChain->SetBranchAddress("pi_plus2_PY", &pi_plus2_PY, &b_pi_plus2_PY);
        fChain->SetBranchAddress("pi_plus2_PZ", &pi_plus2_PZ, &b_pi_plus2_PZ);
        fChain->SetBranchAddress("pi_plus2_ID", &pi_plus2_ID, &b_pi_plus2_ID);
        fChain->SetBranchAddress("pi_plus2_PIDmu", &pi_plus2_PIDmu, &b_pi_plus2_PIDmu);
        fChain->SetBranchAddress("pi_plus2_PIDK", &pi_plus2_PIDK, &b_pi_plus2_PIDK);
        fChain->SetBranchAddress("pi_plus2_PIDp", &pi_plus2_PIDp, &b_pi_plus2_PIDp);
        fChain->SetBranchAddress("pi_plus2_ProbNNk", &pi_plus2_ProbNNk, &b_pi_plus2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_ProbNNp", &pi_plus2_ProbNNp, &b_pi_plus2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_ProbNNpi", &pi_plus2_ProbNNpi, &b_pi_plus2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_ProbNNmu", &pi_plus2_ProbNNmu, &b_pi_plus2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_ProbNNghost", &pi_plus2_ProbNNghost, &b_pi_plus2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_isMuon", &pi_plus2_isMuon, &b_pi_plus2_isMuon);
        fChain->SetBranchAddress("pi_plus2_TRACK_CHI2NDOF", &pi_plus2_TRACK_CHI2NDOF, &b_pi_plus2_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus2_TRACK_GhostProb", &pi_plus2_TRACK_GhostProb, &b_pi_plus2_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus2_ptasy_1.00", &pi_plus2_ptasy_1_00, &b_pi_plus2_ptasy_1_00);
        fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
        fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
        fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
        fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
        fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
        fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
        fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
        fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
        fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
        fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
        fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
        fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
        fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
        fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
        fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
        fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
        fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
        fChain->SetBranchAddress("pi_minus_TRACK_CHI2NDOF", &pi_minus_TRACK_CHI2NDOF, &b_pi_minus_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
        fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
        fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
        fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
        fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
        fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
        fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
        fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
        fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
        fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
    }

    if(_Ds_finalState == Ds_finalState::Kpipi && _decay == Decay::signal){

    	if(!_data){
		fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
		fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
		fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);
	
		fChain->SetBranchAddress("K_plus_TRUEID", &K_plus_TRUEID);
		fChain->SetBranchAddress("pi_plus_TRUEID", &pi_plus_TRUEID);
		fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
		fChain->SetBranchAddress("pi_plus_fromDs_TRUEID", &pi_plus_fromDs_TRUEID);
		fChain->SetBranchAddress("K_minus_fromDs_TRUEID", &K_minus_fromDs_TRUEID);
		fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);
	
		fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
		fChain->SetBranchAddress("K_plus_MC_MOTHER_ID", &K_plus_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_plus_MC_MOTHER_ID", &pi_plus_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_plus_fromDs_MC_MOTHER_ID", &pi_plus_fromDs_MC_MOTHER_ID);
		fChain->SetBranchAddress("K_minus_fromDs_MC_MOTHER_ID", &K_minus_fromDs_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);
	
		fChain->SetBranchAddress("K_plus_PIDK_gen_MagDown", &K_plus_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDK_gen_MagDown", &pi_plus_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDK_gen_MagDown", &pi_minus_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_gen_MagDown", &pi_plus_fromDs_PIDK_gen_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagDown", &K_minus_fromDs_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagDown", &pi_minus_fromDs_PIDK_gen_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDK_gen_MagUp", &K_plus_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDK_gen_MagUp", &pi_plus_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDK_gen_MagUp", &pi_minus_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_gen_MagUp", &pi_plus_fromDs_PIDK_gen_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagUp", &K_minus_fromDs_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagUp", &pi_minus_fromDs_PIDK_gen_MagUp);
	
		fChain->SetBranchAddress("K_plus_PIDK_corr_MagDown", &K_plus_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDK_corr_MagDown", &pi_plus_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDK_corr_MagDown", &pi_minus_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_corr_MagDown", &pi_plus_fromDs_PIDK_corr_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagDown", &K_minus_fromDs_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagDown", &pi_minus_fromDs_PIDK_corr_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDK_corr_MagUp", &K_plus_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDK_corr_MagUp", &pi_plus_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDK_corr_MagUp", &pi_minus_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_corr_MagUp", &pi_plus_fromDs_PIDK_corr_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagUp", &K_minus_fromDs_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagUp", &pi_minus_fromDs_PIDK_corr_MagUp);
	
		/*
		fChain->SetBranchAddress("K_plus_PIDp_gen_MagDown", &K_plus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDp_gen_MagDown", &pi_plus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDp_gen_MagDown", &pi_minus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagDown", &K_plus_fromDs_PIDp_gen_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagDown", &K_minus_fromDs_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagDown", &pi_minus_fromDs_PIDp_gen_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDp_gen_MagUp", &K_plus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDp_gen_MagUp", &pi_plus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDp_gen_MagUp", &pi_minus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagUp", &K_plus_fromDs_PIDp_gen_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagUp", &K_minus_fromDs_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagUp", &pi_minus_fromDs_PIDp_gen_MagUp);
	
		fChain->SetBranchAddress("K_plus_PIDp_corr_MagDown", &K_plus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDp_corr_MagDown", &pi_plus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDp_corr_MagDown", &pi_minus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagDown", &K_plus_fromDs_PIDp_corr_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagDown", &K_minus_fromDs_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagDown", &pi_minus_fromDs_PIDp_corr_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDp_corr_MagUp", &K_plus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDp_corr_MagUp", &pi_plus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDp_corr_MagUp", &pi_minus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagUp", &K_plus_fromDs_PIDp_corr_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagUp", &K_minus_fromDs_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagUp", &pi_minus_fromDs_PIDp_corr_MagUp);
		*/
	}

	fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
	fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
	fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
	fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
	fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
	fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
	fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
	fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
	fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
	fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
	fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
	fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
	fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
	fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
	fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
// 	fChain->SetBranchAddress("Bs_M", &Bs_M, &b_Bs_M);
	fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
// 	fChain->SetBranchAddress("Bs_L0Global_Dec", &Bs_L0Global_Dec, &b_Bs_L0Global_Dec);
	fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
	fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);/*
	fChain->SetBranchAddress("Bs_Hlt1Global_Dec", &Bs_Hlt1Global_Dec, &b_Bs_Hlt1Global_Dec);
	fChain->SetBranchAddress("Bs_Hlt1Global_TIS", &Bs_Hlt1Global_TIS, &b_Bs_Hlt1Global_TIS);
	fChain->SetBranchAddress("Bs_Hlt1Global_TOS", &Bs_Hlt1Global_TOS, &b_Bs_Hlt1Global_TOS);*/
/*	fChain->SetBranchAddress("Bs_Hlt1Phys_Dec", &Bs_Hlt1Phys_Dec, &b_Bs_Hlt1Phys_Dec);
	fChain->SetBranchAddress("Bs_Hlt1Phys_TIS", &Bs_Hlt1Phys_TIS, &b_Bs_Hlt1Phys_TIS);
	fChain->SetBranchAddress("Bs_Hlt1Phys_TOS", &Bs_Hlt1Phys_TOS, &b_Bs_Hlt1Phys_TOS);
	fChain->SetBranchAddress("Bs_Hlt2Global_Dec", &Bs_Hlt2Global_Dec, &b_Bs_Hlt2Global_Dec);
	fChain->SetBranchAddress("Bs_Hlt2Global_TIS", &Bs_Hlt2Global_TIS, &b_Bs_Hlt2Global_TIS);
	fChain->SetBranchAddress("Bs_Hlt2Global_TOS", &Bs_Hlt2Global_TOS, &b_Bs_Hlt2Global_TOS);
	fChain->SetBranchAddress("Bs_Hlt2Phys_Dec", &Bs_Hlt2Phys_Dec, &b_Bs_Hlt2Phys_Dec);
	fChain->SetBranchAddress("Bs_Hlt2Phys_TIS", &Bs_Hlt2Phys_TIS, &b_Bs_Hlt2Phys_TIS);
	fChain->SetBranchAddress("Bs_Hlt2Phys_TOS", &Bs_Hlt2Phys_TOS, &b_Bs_Hlt2Phys_TOS);
	fChain->SetBranchAddress("Bs_L0HadronDecision_Dec", &Bs_L0HadronDecision_Dec, &b_Bs_L0HadronDecision_Dec);*/
	fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
	fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
/*	fChain->SetBranchAddress("Bs_L0MuonDecision_Dec", &Bs_L0MuonDecision_Dec, &b_Bs_L0MuonDecision_Dec);
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
	fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TOS", &Bs_Hlt2Topo4BodyDecision_TOS, &b_Bs_Hlt2Topo4BodyDecision_TOS);*/
	fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
	fChain->SetBranchAddress("Bs_B0DTF_nPV", &Bs_B0DTF_nPV, &b_Bs_B0DTF_nPV);
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
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_ID", Bs_B0DTF_D_splus_piplus_0_ID, &b_Bs_B0DTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PE", Bs_B0DTF_D_splus_piplus_0_PE, &b_Bs_B0DTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PX", Bs_B0DTF_D_splus_piplus_0_PX, &b_Bs_B0DTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PY", Bs_B0DTF_D_splus_piplus_0_PY, &b_Bs_B0DTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PZ", Bs_B0DTF_D_splus_piplus_0_PZ, &b_Bs_B0DTF_D_splus_piplus_0_PZ);
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
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_ID", Bs_BsDTF_D_splus_piplus_0_ID, &b_Bs_BsDTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PE", Bs_BsDTF_D_splus_piplus_0_PE, &b_Bs_BsDTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PX", Bs_BsDTF_D_splus_piplus_0_PX, &b_Bs_BsDTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PY", Bs_BsDTF_D_splus_piplus_0_PY, &b_Bs_BsDTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PZ", Bs_BsDTF_D_splus_piplus_0_PZ, &b_Bs_BsDTF_D_splus_piplus_0_PZ);
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
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_ID", Bs_DTF_D_splus_piplus_0_ID, &b_Bs_DTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PE", Bs_DTF_D_splus_piplus_0_PE, &b_Bs_DTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PX", Bs_DTF_D_splus_piplus_0_PX, &b_Bs_DTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PY", Bs_DTF_D_splus_piplus_0_PY, &b_Bs_DTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PZ", Bs_DTF_D_splus_piplus_0_PZ, &b_Bs_DTF_D_splus_piplus_0_PZ);
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
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_ID", Bs_PV_Dplus_piplus_0_ID, &b_Bs_PV_Dplus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PE", Bs_PV_Dplus_piplus_0_PE, &b_Bs_PV_Dplus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PX", Bs_PV_Dplus_piplus_0_PX, &b_Bs_PV_Dplus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PY", Bs_PV_Dplus_piplus_0_PY, &b_Bs_PV_Dplus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PZ", Bs_PV_Dplus_piplus_0_PZ, &b_Bs_PV_Dplus_piplus_0_PZ);
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
	fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
	fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
	fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
	fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
/*	fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
	fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);*/
// 	fChain->SetBranchAddress("Ds_CosTheta", &Ds_CosTheta, &b_Ds_CosTheta);
	fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
	fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
	fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
	fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
	fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
// 	fChain->SetBranchAddress("Ds_ENDVERTEX_COV_", Ds_ENDVERTEX_COV_, &b_Ds_ENDVERTEX_COV_);
	fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
	fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
	fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
	fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
	fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
	fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
	fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
	fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
// 	fChain->SetBranchAddress("Ds_OWNPV_COV_", Ds_OWNPV_COV_, &b_Ds_OWNPV_COV_);
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
// 	fChain->SetBranchAddress("Ds_ORIVX_COV_", Ds_ORIVX_COV_, &b_Ds_ORIVX_COV_);
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
	fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
	fChain->SetBranchAddress("pi_plus_fromDs_ETA", &pi_plus_fromDs_ETA, &b_pi_plus_fromDs_ETA);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNe", &pi_plus_fromDs_MC12TuneV2_ProbNNe, &b_pi_plus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNmu", &pi_plus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNpi", &pi_plus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNk", &pi_plus_fromDs_MC12TuneV2_ProbNNk, &b_pi_plus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNp", &pi_plus_fromDs_MC12TuneV2_ProbNNp, &b_pi_plus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNghost", &pi_plus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNe", &pi_plus_fromDs_MC12TuneV3_ProbNNe, &b_pi_plus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNmu", &pi_plus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNpi", &pi_plus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNk", &pi_plus_fromDs_MC12TuneV3_ProbNNk, &b_pi_plus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNp", &pi_plus_fromDs_MC12TuneV3_ProbNNp, &b_pi_plus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNghost", &pi_plus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_fromDs_IP_OWNPV", &pi_plus_fromDs_IP_OWNPV, &b_pi_plus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("pi_plus_fromDs_IPCHI2_OWNPV", &pi_plus_fromDs_IPCHI2_OWNPV, &b_pi_plus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_plus_fromDs_P", &pi_plus_fromDs_P, &b_pi_plus_fromDs_P);
	fChain->SetBranchAddress("pi_plus_fromDs_PT", &pi_plus_fromDs_PT, &b_pi_plus_fromDs_PT);
	fChain->SetBranchAddress("pi_plus_fromDs_PE", &pi_plus_fromDs_PE, &b_pi_plus_fromDs_PE);
	fChain->SetBranchAddress("pi_plus_fromDs_PX", &pi_plus_fromDs_PX, &b_pi_plus_fromDs_PX);
	fChain->SetBranchAddress("pi_plus_fromDs_PY", &pi_plus_fromDs_PY, &b_pi_plus_fromDs_PY);
	fChain->SetBranchAddress("pi_plus_fromDs_PZ", &pi_plus_fromDs_PZ, &b_pi_plus_fromDs_PZ);
	fChain->SetBranchAddress("pi_plus_fromDs_ID", &pi_plus_fromDs_ID, &b_pi_plus_fromDs_ID);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDmu", &pi_plus_fromDs_PIDmu, &b_pi_plus_fromDs_PIDmu);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDK", &pi_plus_fromDs_PIDK, &b_pi_plus_fromDs_PIDK);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDp", &pi_plus_fromDs_PIDp, &b_pi_plus_fromDs_PIDp);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNk", &pi_plus_fromDs_ProbNNk, &b_pi_plus_fromDs_ProbNNk);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNp", &pi_plus_fromDs_ProbNNp, &b_pi_plus_fromDs_ProbNNp);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNpi", &pi_plus_fromDs_ProbNNpi, &b_pi_plus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNmu", &pi_plus_fromDs_ProbNNmu, &b_pi_plus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNghost", &pi_plus_fromDs_ProbNNghost, &b_pi_plus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_fromDs_isMuon", &pi_plus_fromDs_isMuon, &b_pi_plus_fromDs_isMuon);
	fChain->SetBranchAddress("pi_plus_fromDs_TRACK_GhostProb", &pi_plus_fromDs_TRACK_GhostProb, &b_pi_plus_fromDs_TRACK_GhostProb);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA1", &pi_minus_fromDs_DOCA1, &b_pi_minus_fromDs_DOCA1);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA2", &pi_minus_fromDs_DOCA2, &b_pi_minus_fromDs_DOCA2);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA3", &pi_minus_fromDs_DOCA3, &b_pi_minus_fromDs_DOCA3);
	fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNe", &pi_minus_fromDs_MC12TuneV2_ProbNNe, &b_pi_minus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNe", &pi_minus_fromDs_MC12TuneV3_ProbNNe, &b_pi_minus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
	fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
	fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
	fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
	fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
	fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
	fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
	fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
	fChain->SetBranchAddress("K_minus_fromDs_ETA", &K_minus_fromDs_ETA, &b_K_minus_fromDs_ETA);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNe", &K_minus_fromDs_MC12TuneV2_ProbNNe, &b_K_minus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNmu", &K_minus_fromDs_MC12TuneV2_ProbNNmu, &b_K_minus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNpi", &K_minus_fromDs_MC12TuneV2_ProbNNpi, &b_K_minus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNk", &K_minus_fromDs_MC12TuneV2_ProbNNk, &b_K_minus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNp", &K_minus_fromDs_MC12TuneV2_ProbNNp, &b_K_minus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNghost", &K_minus_fromDs_MC12TuneV2_ProbNNghost, &b_K_minus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNe", &K_minus_fromDs_MC12TuneV3_ProbNNe, &b_K_minus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNmu", &K_minus_fromDs_MC12TuneV3_ProbNNmu, &b_K_minus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNpi", &K_minus_fromDs_MC12TuneV3_ProbNNpi, &b_K_minus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNk", &K_minus_fromDs_MC12TuneV3_ProbNNk, &b_K_minus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNp", &K_minus_fromDs_MC12TuneV3_ProbNNp, &b_K_minus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNghost", &K_minus_fromDs_MC12TuneV3_ProbNNghost, &b_K_minus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("K_minus_fromDs_IP_OWNPV", &K_minus_fromDs_IP_OWNPV, &b_K_minus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("K_minus_fromDs_IPCHI2_OWNPV", &K_minus_fromDs_IPCHI2_OWNPV, &b_K_minus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("K_minus_fromDs_P", &K_minus_fromDs_P, &b_K_minus_fromDs_P);
	fChain->SetBranchAddress("K_minus_fromDs_PT", &K_minus_fromDs_PT, &b_K_minus_fromDs_PT);
	fChain->SetBranchAddress("K_minus_fromDs_PE", &K_minus_fromDs_PE, &b_K_minus_fromDs_PE);
	fChain->SetBranchAddress("K_minus_fromDs_PX", &K_minus_fromDs_PX, &b_K_minus_fromDs_PX);
	fChain->SetBranchAddress("K_minus_fromDs_PY", &K_minus_fromDs_PY, &b_K_minus_fromDs_PY);
	fChain->SetBranchAddress("K_minus_fromDs_PZ", &K_minus_fromDs_PZ, &b_K_minus_fromDs_PZ);
	fChain->SetBranchAddress("K_minus_fromDs_ID", &K_minus_fromDs_ID, &b_K_minus_fromDs_ID);
	fChain->SetBranchAddress("K_minus_fromDs_PIDmu", &K_minus_fromDs_PIDmu, &b_K_minus_fromDs_PIDmu);
	fChain->SetBranchAddress("K_minus_fromDs_PIDK", &K_minus_fromDs_PIDK, &b_K_minus_fromDs_PIDK);
	fChain->SetBranchAddress("K_minus_fromDs_PIDp", &K_minus_fromDs_PIDp, &b_K_minus_fromDs_PIDp);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNk", &K_minus_fromDs_ProbNNk, &b_K_minus_fromDs_ProbNNk);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNp", &K_minus_fromDs_ProbNNp, &b_K_minus_fromDs_ProbNNp);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNpi", &K_minus_fromDs_ProbNNpi, &b_K_minus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNmu", &K_minus_fromDs_ProbNNmu, &b_K_minus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNghost", &K_minus_fromDs_ProbNNghost, &b_K_minus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("K_minus_fromDs_isMuon", &K_minus_fromDs_isMuon, &b_K_minus_fromDs_isMuon);
// 	fChain->SetBranchAddress("K_minus_fromDs_hasRich", &K_minus_fromDs_hasRich, &b_K_minus_fromDs_hasRich);
// 	fChain->SetBranchAddress("K_minus_fromDs_hasCalo", &K_minus_fromDs_hasCalo, &b_K_minus_fromDs_hasCalo);
	fChain->SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb, &b_K_minus_fromDs_TRACK_GhostProb);
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
        fChain->SetBranchAddress("K_1_1270_plus_ptasy_1.00", &K_1_1270_plus_ptasy_1_00, &b_K_1_1270_plus_ptasy_1_00);
	fChain->SetBranchAddress("K_plus_ETA", &K_plus_ETA, &b_K_plus_ETA);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNe", &K_plus_MC12TuneV2_ProbNNe, &b_K_plus_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNmu", &K_plus_MC12TuneV2_ProbNNmu, &b_K_plus_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNpi", &K_plus_MC12TuneV2_ProbNNpi, &b_K_plus_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNk", &K_plus_MC12TuneV2_ProbNNk, &b_K_plus_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNp", &K_plus_MC12TuneV2_ProbNNp, &b_K_plus_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNghost", &K_plus_MC12TuneV2_ProbNNghost, &b_K_plus_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNe", &K_plus_MC12TuneV3_ProbNNe, &b_K_plus_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNmu", &K_plus_MC12TuneV3_ProbNNmu, &b_K_plus_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNpi", &K_plus_MC12TuneV3_ProbNNpi, &b_K_plus_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNk", &K_plus_MC12TuneV3_ProbNNk, &b_K_plus_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNp", &K_plus_MC12TuneV3_ProbNNp, &b_K_plus_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNghost", &K_plus_MC12TuneV3_ProbNNghost, &b_K_plus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("K_plus_IP_OWNPV", &K_plus_IP_OWNPV, &b_K_plus_IP_OWNPV);
	fChain->SetBranchAddress("K_plus_IPCHI2_OWNPV", &K_plus_IPCHI2_OWNPV, &b_K_plus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("K_plus_P", &K_plus_P, &b_K_plus_P);
	fChain->SetBranchAddress("K_plus_PT", &K_plus_PT, &b_K_plus_PT);
	fChain->SetBranchAddress("K_plus_PE", &K_plus_PE, &b_K_plus_PE);
	fChain->SetBranchAddress("K_plus_PX", &K_plus_PX, &b_K_plus_PX);
	fChain->SetBranchAddress("K_plus_PY", &K_plus_PY, &b_K_plus_PY);
	fChain->SetBranchAddress("K_plus_PZ", &K_plus_PZ, &b_K_plus_PZ);
	fChain->SetBranchAddress("K_plus_ID", &K_plus_ID, &b_K_plus_ID);
	fChain->SetBranchAddress("K_plus_PIDmu", &K_plus_PIDmu, &b_K_plus_PIDmu);
	fChain->SetBranchAddress("K_plus_PIDK", &K_plus_PIDK, &b_K_plus_PIDK);
	fChain->SetBranchAddress("K_plus_PIDp", &K_plus_PIDp, &b_K_plus_PIDp);
	fChain->SetBranchAddress("K_plus_ProbNNk", &K_plus_ProbNNk, &b_K_plus_ProbNNk);
	fChain->SetBranchAddress("K_plus_ProbNNp", &K_plus_ProbNNp, &b_K_plus_ProbNNp);
	fChain->SetBranchAddress("K_plus_ProbNNpi", &K_plus_ProbNNpi, &b_K_plus_ProbNNpi);
	fChain->SetBranchAddress("K_plus_ProbNNmu", &K_plus_ProbNNmu, &b_K_plus_ProbNNmu);
	fChain->SetBranchAddress("K_plus_ProbNNghost", &K_plus_ProbNNghost, &b_K_plus_ProbNNghost);
	fChain->SetBranchAddress("K_plus_isMuon", &K_plus_isMuon, &b_K_plus_isMuon);
	fChain->SetBranchAddress("K_plus_TRACK_GhostProb", &K_plus_TRACK_GhostProb, &b_K_plus_TRACK_GhostProb);
	fChain->SetBranchAddress("pi_plus_ETA", &pi_plus_ETA, &b_pi_plus_ETA);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNe", &pi_plus_MC12TuneV2_ProbNNe, &b_pi_plus_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNmu", &pi_plus_MC12TuneV2_ProbNNmu, &b_pi_plus_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNpi", &pi_plus_MC12TuneV2_ProbNNpi, &b_pi_plus_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNk", &pi_plus_MC12TuneV2_ProbNNk, &b_pi_plus_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNp", &pi_plus_MC12TuneV2_ProbNNp, &b_pi_plus_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV2_ProbNNghost", &pi_plus_MC12TuneV2_ProbNNghost, &b_pi_plus_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNe", &pi_plus_MC12TuneV3_ProbNNe, &b_pi_plus_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNmu", &pi_plus_MC12TuneV3_ProbNNmu, &b_pi_plus_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNpi", &pi_plus_MC12TuneV3_ProbNNpi, &b_pi_plus_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNk", &pi_plus_MC12TuneV3_ProbNNk, &b_pi_plus_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNp", &pi_plus_MC12TuneV3_ProbNNp, &b_pi_plus_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_MC12TuneV3_ProbNNghost", &pi_plus_MC12TuneV3_ProbNNghost, &b_pi_plus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_IP_OWNPV", &pi_plus_IP_OWNPV, &b_pi_plus_IP_OWNPV);
	fChain->SetBranchAddress("pi_plus_IPCHI2_OWNPV", &pi_plus_IPCHI2_OWNPV, &b_pi_plus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_plus_P", &pi_plus_P, &b_pi_plus_P);
	fChain->SetBranchAddress("pi_plus_PT", &pi_plus_PT, &b_pi_plus_PT);
	fChain->SetBranchAddress("pi_plus_PE", &pi_plus_PE, &b_pi_plus_PE);
	fChain->SetBranchAddress("pi_plus_PX", &pi_plus_PX, &b_pi_plus_PX);
	fChain->SetBranchAddress("pi_plus_PY", &pi_plus_PY, &b_pi_plus_PY);
	fChain->SetBranchAddress("pi_plus_PZ", &pi_plus_PZ, &b_pi_plus_PZ);
	fChain->SetBranchAddress("pi_plus_ID", &pi_plus_ID, &b_pi_plus_ID);
	fChain->SetBranchAddress("pi_plus_PIDmu", &pi_plus_PIDmu, &b_pi_plus_PIDmu);
	fChain->SetBranchAddress("pi_plus_PIDK", &pi_plus_PIDK, &b_pi_plus_PIDK);
	fChain->SetBranchAddress("pi_plus_PIDp", &pi_plus_PIDp, &b_pi_plus_PIDp);
	fChain->SetBranchAddress("pi_plus_ProbNNk", &pi_plus_ProbNNk, &b_pi_plus_ProbNNk);
	fChain->SetBranchAddress("pi_plus_ProbNNp", &pi_plus_ProbNNp, &b_pi_plus_ProbNNp);
	fChain->SetBranchAddress("pi_plus_ProbNNpi", &pi_plus_ProbNNpi, &b_pi_plus_ProbNNpi);
	fChain->SetBranchAddress("pi_plus_ProbNNmu", &pi_plus_ProbNNmu, &b_pi_plus_ProbNNmu);
	fChain->SetBranchAddress("pi_plus_ProbNNghost", &pi_plus_ProbNNghost, &b_pi_plus_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_isMuon", &pi_plus_isMuon, &b_pi_plus_isMuon);
	fChain->SetBranchAddress("pi_plus_TRACK_GhostProb", &pi_plus_TRACK_GhostProb, &b_pi_plus_TRACK_GhostProb);
	fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNe", &pi_minus_MC12TuneV2_ProbNNe, &b_pi_minus_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNe", &pi_minus_MC12TuneV3_ProbNNe, &b_pi_minus_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
	fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
	fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
	fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
	fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
	fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
	fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
	fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
	fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
	fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
	fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
	fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
	fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
	fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
	fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
	fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
	fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
	fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
	fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
	fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
	//    fChain->SetBranchAddress("nLong", &nLong, &b_nLong);
	fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
	fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    }

    if(_Ds_finalState == Ds_finalState::Kpipi && _decay == Decay::norm){

    	if(!_data){
		fChain->SetBranchAddress("Bs_TRUEID", &Bs_TRUEID);
		fChain->SetBranchAddress("Ds_TRUEID", &Ds_TRUEID);
		fChain->SetBranchAddress("Bs_BKGCAT", &Bs_BKGCAT);
	
                fChain->SetBranchAddress("pi_minus_TRUEID", &pi_minus_TRUEID);
    	        fChain->SetBranchAddress("pi_plus1_TRUEID", &pi_plus1_TRUEID);
    	        fChain->SetBranchAddress("pi_plus2_TRUEID", &pi_plus2_TRUEID);
		fChain->SetBranchAddress("pi_plus_fromDs_TRUEID", &pi_plus_fromDs_TRUEID);
		fChain->SetBranchAddress("K_minus_fromDs_TRUEID", &K_minus_fromDs_TRUEID);
		fChain->SetBranchAddress("pi_minus_fromDs_TRUEID", &pi_minus_fromDs_TRUEID);
	
		fChain->SetBranchAddress("Ds_MC_MOTHER_ID", &Ds_MC_MOTHER_ID);
            fChain->SetBranchAddress("pi_minus_MC_MOTHER_ID", &pi_minus_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus1_MC_MOTHER_ID", &pi_plus1_MC_MOTHER_ID);
    	    fChain->SetBranchAddress("pi_plus2_MC_MOTHER_ID", &pi_plus2_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_plus_fromDs_MC_MOTHER_ID", &pi_plus_fromDs_MC_MOTHER_ID);
		fChain->SetBranchAddress("K_minus_fromDs_MC_MOTHER_ID", &K_minus_fromDs_MC_MOTHER_ID);
		fChain->SetBranchAddress("pi_minus_fromDs_MC_MOTHER_ID", &pi_minus_fromDs_MC_MOTHER_ID);
	
	    fChain->SetBranchAddress("pi_plus1_PIDK_gen_MagDown", &pi_plus1_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_gen_MagDown", &pi_plus2_PIDK_gen_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagDown", &pi_minus_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_gen_MagDown", &pi_plus_fromDs_PIDK_gen_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagDown", &K_minus_fromDs_PIDK_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagDown", &pi_minus_fromDs_PIDK_gen_MagDown);
	
	    fChain->SetBranchAddress("pi_plus1_PIDK_gen_MagUp", &pi_plus1_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_gen_MagUp", &pi_plus2_PIDK_gen_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_gen_MagUp", &pi_minus_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_gen_MagUp", &pi_plus_fromDs_PIDK_gen_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_gen_MagUp", &K_minus_fromDs_PIDK_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_gen_MagUp", &pi_minus_fromDs_PIDK_gen_MagUp);
	
	    fChain->SetBranchAddress("pi_plus1_PIDK_corr_MagDown", &pi_plus1_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_corr_MagDown", &pi_plus2_PIDK_corr_MagDown);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagDown", &pi_minus_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_corr_MagDown", &pi_plus_fromDs_PIDK_corr_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagDown", &K_minus_fromDs_PIDK_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagDown", &pi_minus_fromDs_PIDK_corr_MagDown);
	
	    fChain->SetBranchAddress("pi_plus1_PIDK_corr_MagUp", &pi_plus1_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_plus2_PIDK_corr_MagUp", &pi_plus2_PIDK_corr_MagUp);
    	    fChain->SetBranchAddress("pi_minus_PIDK_corr_MagUp", &pi_minus_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_plus_fromDs_PIDK_corr_MagUp", &pi_plus_fromDs_PIDK_corr_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDK_corr_MagUp", &K_minus_fromDs_PIDK_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDK_corr_MagUp", &pi_minus_fromDs_PIDK_corr_MagUp);
	
		/*
		fChain->SetBranchAddress("K_plus_PIDp_gen_MagDown", &K_plus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDp_gen_MagDown", &pi_plus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDp_gen_MagDown", &pi_minus_PIDp_gen_MagDown);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagDown", &K_plus_fromDs_PIDp_gen_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagDown", &K_minus_fromDs_PIDp_gen_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagDown", &pi_minus_fromDs_PIDp_gen_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDp_gen_MagUp", &K_plus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDp_gen_MagUp", &pi_plus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDp_gen_MagUp", &pi_minus_PIDp_gen_MagUp);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_gen_MagUp", &K_plus_fromDs_PIDp_gen_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_gen_MagUp", &K_minus_fromDs_PIDp_gen_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_gen_MagUp", &pi_minus_fromDs_PIDp_gen_MagUp);
	
		fChain->SetBranchAddress("K_plus_PIDp_corr_MagDown", &K_plus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_plus_PIDp_corr_MagDown", &pi_plus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_PIDp_corr_MagDown", &pi_minus_PIDp_corr_MagDown);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagDown", &K_plus_fromDs_PIDp_corr_MagDown);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagDown", &K_minus_fromDs_PIDp_corr_MagDown);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagDown", &pi_minus_fromDs_PIDp_corr_MagDown);
	
		fChain->SetBranchAddress("K_plus_PIDp_corr_MagUp", &K_plus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_plus_PIDp_corr_MagUp", &pi_plus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_PIDp_corr_MagUp", &pi_minus_PIDp_corr_MagUp);
		fChain->SetBranchAddress("K_plus_fromDs_PIDp_corr_MagUp", &K_plus_fromDs_PIDp_corr_MagUp);
		fChain->SetBranchAddress("K_minus_fromDs_PIDp_corr_MagUp", &K_minus_fromDs_PIDp_corr_MagUp);
		fChain->SetBranchAddress("pi_minus_fromDs_PIDp_corr_MagUp", &pi_minus_fromDs_PIDp_corr_MagUp);
		*/
	}

	fChain->SetBranchAddress("Bs_ETA", &Bs_ETA, &b_Bs_ETA);
	fChain->SetBranchAddress("Bs_TAU", &Bs_TAU, &b_Bs_TAU);
	fChain->SetBranchAddress("Bs_TAUERR", &Bs_TAUERR, &b_Bs_TAUERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_X", &Bs_ENDVERTEX_X, &b_Bs_ENDVERTEX_X);
	fChain->SetBranchAddress("Bs_ENDVERTEX_Y", &Bs_ENDVERTEX_Y, &b_Bs_ENDVERTEX_Y);
	fChain->SetBranchAddress("Bs_ENDVERTEX_Z", &Bs_ENDVERTEX_Z, &b_Bs_ENDVERTEX_Z);
	fChain->SetBranchAddress("Bs_ENDVERTEX_XERR", &Bs_ENDVERTEX_XERR, &b_Bs_ENDVERTEX_XERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_YERR", &Bs_ENDVERTEX_YERR, &b_Bs_ENDVERTEX_YERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_ZERR", &Bs_ENDVERTEX_ZERR, &b_Bs_ENDVERTEX_ZERR);
	fChain->SetBranchAddress("Bs_ENDVERTEX_CHI2", &Bs_ENDVERTEX_CHI2, &b_Bs_ENDVERTEX_CHI2);
	fChain->SetBranchAddress("Bs_ENDVERTEX_NDOF", &Bs_ENDVERTEX_NDOF, &b_Bs_ENDVERTEX_NDOF);
	fChain->SetBranchAddress("Bs_OWNPV_X", &Bs_OWNPV_X, &b_Bs_OWNPV_X);
	fChain->SetBranchAddress("Bs_OWNPV_Y", &Bs_OWNPV_Y, &b_Bs_OWNPV_Y);
	fChain->SetBranchAddress("Bs_OWNPV_Z", &Bs_OWNPV_Z, &b_Bs_OWNPV_Z);
	fChain->SetBranchAddress("Bs_OWNPV_XERR", &Bs_OWNPV_XERR, &b_Bs_OWNPV_XERR);
	fChain->SetBranchAddress("Bs_OWNPV_YERR", &Bs_OWNPV_YERR, &b_Bs_OWNPV_YERR);
	fChain->SetBranchAddress("Bs_OWNPV_ZERR", &Bs_OWNPV_ZERR, &b_Bs_OWNPV_ZERR);
	fChain->SetBranchAddress("Bs_OWNPV_CHI2", &Bs_OWNPV_CHI2, &b_Bs_OWNPV_CHI2);
	fChain->SetBranchAddress("Bs_OWNPV_NDOF", &Bs_OWNPV_NDOF, &b_Bs_OWNPV_NDOF);
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
// 	fChain->SetBranchAddress("Bs_M", &Bs_M, &b_Bs_M);
	fChain->SetBranchAddress("Bs_ID", &Bs_ID, &b_Bs_ID);
// 	fChain->SetBranchAddress("Bs_L0Global_Dec", &Bs_L0Global_Dec, &b_Bs_L0Global_Dec);
	fChain->SetBranchAddress("Bs_L0Global_TIS", &Bs_L0Global_TIS, &b_Bs_L0Global_TIS);
	fChain->SetBranchAddress("Bs_L0Global_TOS", &Bs_L0Global_TOS, &b_Bs_L0Global_TOS);/*
	fChain->SetBranchAddress("Bs_Hlt1Global_Dec", &Bs_Hlt1Global_Dec, &b_Bs_Hlt1Global_Dec);
	fChain->SetBranchAddress("Bs_Hlt1Global_TIS", &Bs_Hlt1Global_TIS, &b_Bs_Hlt1Global_TIS);
	fChain->SetBranchAddress("Bs_Hlt1Global_TOS", &Bs_Hlt1Global_TOS, &b_Bs_Hlt1Global_TOS);*/
/*	fChain->SetBranchAddress("Bs_Hlt1Phys_Dec", &Bs_Hlt1Phys_Dec, &b_Bs_Hlt1Phys_Dec);
	fChain->SetBranchAddress("Bs_Hlt1Phys_TIS", &Bs_Hlt1Phys_TIS, &b_Bs_Hlt1Phys_TIS);
	fChain->SetBranchAddress("Bs_Hlt1Phys_TOS", &Bs_Hlt1Phys_TOS, &b_Bs_Hlt1Phys_TOS);
	fChain->SetBranchAddress("Bs_Hlt2Global_Dec", &Bs_Hlt2Global_Dec, &b_Bs_Hlt2Global_Dec);
	fChain->SetBranchAddress("Bs_Hlt2Global_TIS", &Bs_Hlt2Global_TIS, &b_Bs_Hlt2Global_TIS);
	fChain->SetBranchAddress("Bs_Hlt2Global_TOS", &Bs_Hlt2Global_TOS, &b_Bs_Hlt2Global_TOS);
	fChain->SetBranchAddress("Bs_Hlt2Phys_Dec", &Bs_Hlt2Phys_Dec, &b_Bs_Hlt2Phys_Dec);
	fChain->SetBranchAddress("Bs_Hlt2Phys_TIS", &Bs_Hlt2Phys_TIS, &b_Bs_Hlt2Phys_TIS);
	fChain->SetBranchAddress("Bs_Hlt2Phys_TOS", &Bs_Hlt2Phys_TOS, &b_Bs_Hlt2Phys_TOS);
	fChain->SetBranchAddress("Bs_L0HadronDecision_Dec", &Bs_L0HadronDecision_Dec, &b_Bs_L0HadronDecision_Dec);*/
	fChain->SetBranchAddress("Bs_L0HadronDecision_TIS", &Bs_L0HadronDecision_TIS, &b_Bs_L0HadronDecision_TIS);
	fChain->SetBranchAddress("Bs_L0HadronDecision_TOS", &Bs_L0HadronDecision_TOS, &b_Bs_L0HadronDecision_TOS);
/*	fChain->SetBranchAddress("Bs_L0MuonDecision_Dec", &Bs_L0MuonDecision_Dec, &b_Bs_L0MuonDecision_Dec);
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
	fChain->SetBranchAddress("Bs_Hlt2Topo4BodyDecision_TOS", &Bs_Hlt2Topo4BodyDecision_TOS, &b_Bs_Hlt2Topo4BodyDecision_TOS);*/
	fChain->SetBranchAddress("Bs_ptasy_1.00", &Bs_ptasy_1_00, &b_Bs_ptasy_1_00);
	fChain->SetBranchAddress("Bs_B0DTF_nPV", &Bs_B0DTF_nPV, &b_Bs_B0DTF_nPV);
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
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_ID", Bs_B0DTF_D_splus_piplus_0_ID, &b_Bs_B0DTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PE", Bs_B0DTF_D_splus_piplus_0_PE, &b_Bs_B0DTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PX", Bs_B0DTF_D_splus_piplus_0_PX, &b_Bs_B0DTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PY", Bs_B0DTF_D_splus_piplus_0_PY, &b_Bs_B0DTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_0_PZ", Bs_B0DTF_D_splus_piplus_0_PZ, &b_Bs_B0DTF_D_splus_piplus_0_PZ);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_ID", Bs_B0DTF_D_splus_piplus_ID, &b_Bs_B0DTF_D_splus_piplus_ID);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PE", Bs_B0DTF_D_splus_piplus_PE, &b_Bs_B0DTF_D_splus_piplus_PE);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PX", Bs_B0DTF_D_splus_piplus_PX, &b_Bs_B0DTF_D_splus_piplus_PX);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PY", Bs_B0DTF_D_splus_piplus_PY, &b_Bs_B0DTF_D_splus_piplus_PY);
	fChain->SetBranchAddress("Bs_B0DTF_D_splus_piplus_PZ", Bs_B0DTF_D_splus_piplus_PZ, &b_Bs_B0DTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_M", Bs_B0DTF_a_1_1260_plus_M, &b_Bs_B0DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_MERR", Bs_B0DTF_a_1_1260_plus_MERR, &b_Bs_B0DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_P", Bs_B0DTF_a_1_1260_plus_P, &b_Bs_B0DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_PERR", Bs_B0DTF_a_1_1260_plus_PERR, &b_Bs_B0DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctau", Bs_B0DTF_a_1_1260_plus_ctau, &b_Bs_B0DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_ctauErr", Bs_B0DTF_a_1_1260_plus_ctauErr, &b_Bs_B0DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLength", Bs_B0DTF_a_1_1260_plus_decayLength, &b_Bs_B0DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_decayLengthErr", Bs_B0DTF_a_1_1260_plus_decayLengthErr, &b_Bs_B0DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_ID", Bs_B0DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PE", Bs_B0DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PX", Bs_B0DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PY", Bs_B0DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_0_PZ", Bs_B0DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_ID", Bs_B0DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PE", Bs_B0DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PX", Bs_B0DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PY", Bs_B0DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_1_PZ", Bs_B0DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_ID", Bs_B0DTF_a_1_1260_plus_piplus_ID, &b_Bs_B0DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PE", Bs_B0DTF_a_1_1260_plus_piplus_PE, &b_Bs_B0DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PX", Bs_B0DTF_a_1_1260_plus_piplus_PX, &b_Bs_B0DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PY", Bs_B0DTF_a_1_1260_plus_piplus_PY, &b_Bs_B0DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_B0DTF_a_1_1260_plus_piplus_PZ", Bs_B0DTF_a_1_1260_plus_piplus_PZ, &b_Bs_B0DTF_a_1_1260_plus_piplus_PZ);
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
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_ID", Bs_BsDTF_D_splus_piplus_0_ID, &b_Bs_BsDTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PE", Bs_BsDTF_D_splus_piplus_0_PE, &b_Bs_BsDTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PX", Bs_BsDTF_D_splus_piplus_0_PX, &b_Bs_BsDTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PY", Bs_BsDTF_D_splus_piplus_0_PY, &b_Bs_BsDTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_0_PZ", Bs_BsDTF_D_splus_piplus_0_PZ, &b_Bs_BsDTF_D_splus_piplus_0_PZ);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_ID", Bs_BsDTF_D_splus_piplus_ID, &b_Bs_BsDTF_D_splus_piplus_ID);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PE", Bs_BsDTF_D_splus_piplus_PE, &b_Bs_BsDTF_D_splus_piplus_PE);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PX", Bs_BsDTF_D_splus_piplus_PX, &b_Bs_BsDTF_D_splus_piplus_PX);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PY", Bs_BsDTF_D_splus_piplus_PY, &b_Bs_BsDTF_D_splus_piplus_PY);
	fChain->SetBranchAddress("Bs_BsDTF_D_splus_piplus_PZ", Bs_BsDTF_D_splus_piplus_PZ, &b_Bs_BsDTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_M", Bs_BsDTF_a_1_1260_plus_M, &b_Bs_BsDTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_MERR", Bs_BsDTF_a_1_1260_plus_MERR, &b_Bs_BsDTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_P", Bs_BsDTF_a_1_1260_plus_P, &b_Bs_BsDTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_PERR", Bs_BsDTF_a_1_1260_plus_PERR, &b_Bs_BsDTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctau", Bs_BsDTF_a_1_1260_plus_ctau, &b_Bs_BsDTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_ctauErr", Bs_BsDTF_a_1_1260_plus_ctauErr, &b_Bs_BsDTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLength", Bs_BsDTF_a_1_1260_plus_decayLength, &b_Bs_BsDTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_decayLengthErr", Bs_BsDTF_a_1_1260_plus_decayLengthErr, &b_Bs_BsDTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_ID", Bs_BsDTF_a_1_1260_plus_piplus_0_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PE", Bs_BsDTF_a_1_1260_plus_piplus_0_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PX", Bs_BsDTF_a_1_1260_plus_piplus_0_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PY", Bs_BsDTF_a_1_1260_plus_piplus_0_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_0_PZ", Bs_BsDTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_ID", Bs_BsDTF_a_1_1260_plus_piplus_1_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PE", Bs_BsDTF_a_1_1260_plus_piplus_1_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PX", Bs_BsDTF_a_1_1260_plus_piplus_1_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PY", Bs_BsDTF_a_1_1260_plus_piplus_1_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_1_PZ", Bs_BsDTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_ID", Bs_BsDTF_a_1_1260_plus_piplus_ID, &b_Bs_BsDTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PE", Bs_BsDTF_a_1_1260_plus_piplus_PE, &b_Bs_BsDTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PX", Bs_BsDTF_a_1_1260_plus_piplus_PX, &b_Bs_BsDTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PY", Bs_BsDTF_a_1_1260_plus_piplus_PY, &b_Bs_BsDTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_BsDTF_a_1_1260_plus_piplus_PZ", Bs_BsDTF_a_1_1260_plus_piplus_PZ, &b_Bs_BsDTF_a_1_1260_plus_piplus_PZ);
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
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_ID", Bs_DTF_D_splus_piplus_0_ID, &b_Bs_DTF_D_splus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PE", Bs_DTF_D_splus_piplus_0_PE, &b_Bs_DTF_D_splus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PX", Bs_DTF_D_splus_piplus_0_PX, &b_Bs_DTF_D_splus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PY", Bs_DTF_D_splus_piplus_0_PY, &b_Bs_DTF_D_splus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_0_PZ", Bs_DTF_D_splus_piplus_0_PZ, &b_Bs_DTF_D_splus_piplus_0_PZ);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_ID", Bs_DTF_D_splus_piplus_ID, &b_Bs_DTF_D_splus_piplus_ID);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PE", Bs_DTF_D_splus_piplus_PE, &b_Bs_DTF_D_splus_piplus_PE);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PX", Bs_DTF_D_splus_piplus_PX, &b_Bs_DTF_D_splus_piplus_PX);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PY", Bs_DTF_D_splus_piplus_PY, &b_Bs_DTF_D_splus_piplus_PY);
	fChain->SetBranchAddress("Bs_DTF_D_splus_piplus_PZ", Bs_DTF_D_splus_piplus_PZ, &b_Bs_DTF_D_splus_piplus_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_M", Bs_DTF_a_1_1260_plus_M, &b_Bs_DTF_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_MERR", Bs_DTF_a_1_1260_plus_MERR, &b_Bs_DTF_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_P", Bs_DTF_a_1_1260_plus_P, &b_Bs_DTF_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_PERR", Bs_DTF_a_1_1260_plus_PERR, &b_Bs_DTF_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctau", Bs_DTF_a_1_1260_plus_ctau, &b_Bs_DTF_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_ctauErr", Bs_DTF_a_1_1260_plus_ctauErr, &b_Bs_DTF_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLength", Bs_DTF_a_1_1260_plus_decayLength, &b_Bs_DTF_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_decayLengthErr", Bs_DTF_a_1_1260_plus_decayLengthErr, &b_Bs_DTF_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_ID", Bs_DTF_a_1_1260_plus_piplus_0_ID, &b_Bs_DTF_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PE", Bs_DTF_a_1_1260_plus_piplus_0_PE, &b_Bs_DTF_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PX", Bs_DTF_a_1_1260_plus_piplus_0_PX, &b_Bs_DTF_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PY", Bs_DTF_a_1_1260_plus_piplus_0_PY, &b_Bs_DTF_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_0_PZ", Bs_DTF_a_1_1260_plus_piplus_0_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_ID", Bs_DTF_a_1_1260_plus_piplus_1_ID, &b_Bs_DTF_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PE", Bs_DTF_a_1_1260_plus_piplus_1_PE, &b_Bs_DTF_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PX", Bs_DTF_a_1_1260_plus_piplus_1_PX, &b_Bs_DTF_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PY", Bs_DTF_a_1_1260_plus_piplus_1_PY, &b_Bs_DTF_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_1_PZ", Bs_DTF_a_1_1260_plus_piplus_1_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_ID", Bs_DTF_a_1_1260_plus_piplus_ID, &b_Bs_DTF_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PE", Bs_DTF_a_1_1260_plus_piplus_PE, &b_Bs_DTF_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PX", Bs_DTF_a_1_1260_plus_piplus_PX, &b_Bs_DTF_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PY", Bs_DTF_a_1_1260_plus_piplus_PY, &b_Bs_DTF_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_DTF_a_1_1260_plus_piplus_PZ", Bs_DTF_a_1_1260_plus_piplus_PZ, &b_Bs_DTF_a_1_1260_plus_piplus_PZ);
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
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_ID", Bs_PV_Dplus_piplus_0_ID, &b_Bs_PV_Dplus_piplus_0_ID);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PE", Bs_PV_Dplus_piplus_0_PE, &b_Bs_PV_Dplus_piplus_0_PE);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PX", Bs_PV_Dplus_piplus_0_PX, &b_Bs_PV_Dplus_piplus_0_PX);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PY", Bs_PV_Dplus_piplus_0_PY, &b_Bs_PV_Dplus_piplus_0_PY);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_0_PZ", Bs_PV_Dplus_piplus_0_PZ, &b_Bs_PV_Dplus_piplus_0_PZ);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_ID", Bs_PV_Dplus_piplus_ID, &b_Bs_PV_Dplus_piplus_ID);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PE", Bs_PV_Dplus_piplus_PE, &b_Bs_PV_Dplus_piplus_PE);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PX", Bs_PV_Dplus_piplus_PX, &b_Bs_PV_Dplus_piplus_PX);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PY", Bs_PV_Dplus_piplus_PY, &b_Bs_PV_Dplus_piplus_PY);
	fChain->SetBranchAddress("Bs_PV_Dplus_piplus_PZ", Bs_PV_Dplus_piplus_PZ, &b_Bs_PV_Dplus_piplus_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_M", Bs_PV_a_1_1260_plus_M, &b_Bs_PV_a_1_1260_plus_M);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_MERR", Bs_PV_a_1_1260_plus_MERR, &b_Bs_PV_a_1_1260_plus_MERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_P", Bs_PV_a_1_1260_plus_P, &b_Bs_PV_a_1_1260_plus_P);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_PERR", Bs_PV_a_1_1260_plus_PERR, &b_Bs_PV_a_1_1260_plus_PERR);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctau", Bs_PV_a_1_1260_plus_ctau, &b_Bs_PV_a_1_1260_plus_ctau);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_ctauErr", Bs_PV_a_1_1260_plus_ctauErr, &b_Bs_PV_a_1_1260_plus_ctauErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLength", Bs_PV_a_1_1260_plus_decayLength, &b_Bs_PV_a_1_1260_plus_decayLength);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_decayLengthErr", Bs_PV_a_1_1260_plus_decayLengthErr, &b_Bs_PV_a_1_1260_plus_decayLengthErr);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_ID", Bs_PV_a_1_1260_plus_piplus_0_ID, &b_Bs_PV_a_1_1260_plus_piplus_0_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PE", Bs_PV_a_1_1260_plus_piplus_0_PE, &b_Bs_PV_a_1_1260_plus_piplus_0_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PX", Bs_PV_a_1_1260_plus_piplus_0_PX, &b_Bs_PV_a_1_1260_plus_piplus_0_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PY", Bs_PV_a_1_1260_plus_piplus_0_PY, &b_Bs_PV_a_1_1260_plus_piplus_0_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_0_PZ", Bs_PV_a_1_1260_plus_piplus_0_PZ, &b_Bs_PV_a_1_1260_plus_piplus_0_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_ID", Bs_PV_a_1_1260_plus_piplus_1_ID, &b_Bs_PV_a_1_1260_plus_piplus_1_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PE", Bs_PV_a_1_1260_plus_piplus_1_PE, &b_Bs_PV_a_1_1260_plus_piplus_1_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PX", Bs_PV_a_1_1260_plus_piplus_1_PX, &b_Bs_PV_a_1_1260_plus_piplus_1_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PY", Bs_PV_a_1_1260_plus_piplus_1_PY, &b_Bs_PV_a_1_1260_plus_piplus_1_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_1_PZ", Bs_PV_a_1_1260_plus_piplus_1_PZ, &b_Bs_PV_a_1_1260_plus_piplus_1_PZ);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_ID", Bs_PV_a_1_1260_plus_piplus_ID, &b_Bs_PV_a_1_1260_plus_piplus_ID);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PE", Bs_PV_a_1_1260_plus_piplus_PE, &b_Bs_PV_a_1_1260_plus_piplus_PE);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PX", Bs_PV_a_1_1260_plus_piplus_PX, &b_Bs_PV_a_1_1260_plus_piplus_PX);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PY", Bs_PV_a_1_1260_plus_piplus_PY, &b_Bs_PV_a_1_1260_plus_piplus_PY);
        fChain->SetBranchAddress("Bs_PV_a_1_1260_plus_piplus_PZ", Bs_PV_a_1_1260_plus_piplus_PZ, &b_Bs_PV_a_1_1260_plus_piplus_PZ);
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
	fChain->SetBranchAddress("Ds_DOCA1", &Ds_DOCA1, &b_Ds_DOCA1);
	fChain->SetBranchAddress("Ds_DOCA2", &Ds_DOCA2, &b_Ds_DOCA2);
	fChain->SetBranchAddress("Ds_DOCA3", &Ds_DOCA3, &b_Ds_DOCA3);
	fChain->SetBranchAddress("Ds_ETA", &Ds_ETA, &b_Ds_ETA);
/*	fChain->SetBranchAddress("Ds_TAU", &Ds_TAU, &b_Ds_TAU);
	fChain->SetBranchAddress("Ds_TAUERR", &Ds_TAUERR, &b_Ds_TAUERR);*/
// 	fChain->SetBranchAddress("Ds_CosTheta", &Ds_CosTheta, &b_Ds_CosTheta);
	fChain->SetBranchAddress("Ds_ENDVERTEX_X", &Ds_ENDVERTEX_X, &b_Ds_ENDVERTEX_X);
	fChain->SetBranchAddress("Ds_ENDVERTEX_Y", &Ds_ENDVERTEX_Y, &b_Ds_ENDVERTEX_Y);
	fChain->SetBranchAddress("Ds_ENDVERTEX_Z", &Ds_ENDVERTEX_Z, &b_Ds_ENDVERTEX_Z);
	fChain->SetBranchAddress("Ds_ENDVERTEX_XERR", &Ds_ENDVERTEX_XERR, &b_Ds_ENDVERTEX_XERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_YERR", &Ds_ENDVERTEX_YERR, &b_Ds_ENDVERTEX_YERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_ZERR", &Ds_ENDVERTEX_ZERR, &b_Ds_ENDVERTEX_ZERR);
	fChain->SetBranchAddress("Ds_ENDVERTEX_CHI2", &Ds_ENDVERTEX_CHI2, &b_Ds_ENDVERTEX_CHI2);
	fChain->SetBranchAddress("Ds_ENDVERTEX_NDOF", &Ds_ENDVERTEX_NDOF, &b_Ds_ENDVERTEX_NDOF);
// 	fChain->SetBranchAddress("Ds_ENDVERTEX_COV_", Ds_ENDVERTEX_COV_, &b_Ds_ENDVERTEX_COV_);
	fChain->SetBranchAddress("Ds_OWNPV_X", &Ds_OWNPV_X, &b_Ds_OWNPV_X);
	fChain->SetBranchAddress("Ds_OWNPV_Y", &Ds_OWNPV_Y, &b_Ds_OWNPV_Y);
	fChain->SetBranchAddress("Ds_OWNPV_Z", &Ds_OWNPV_Z, &b_Ds_OWNPV_Z);
	fChain->SetBranchAddress("Ds_OWNPV_XERR", &Ds_OWNPV_XERR, &b_Ds_OWNPV_XERR);
	fChain->SetBranchAddress("Ds_OWNPV_YERR", &Ds_OWNPV_YERR, &b_Ds_OWNPV_YERR);
	fChain->SetBranchAddress("Ds_OWNPV_ZERR", &Ds_OWNPV_ZERR, &b_Ds_OWNPV_ZERR);
	fChain->SetBranchAddress("Ds_OWNPV_CHI2", &Ds_OWNPV_CHI2, &b_Ds_OWNPV_CHI2);
	fChain->SetBranchAddress("Ds_OWNPV_NDOF", &Ds_OWNPV_NDOF, &b_Ds_OWNPV_NDOF);
// 	fChain->SetBranchAddress("Ds_OWNPV_COV_", Ds_OWNPV_COV_, &b_Ds_OWNPV_COV_);
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
// 	fChain->SetBranchAddress("Ds_ORIVX_COV_", Ds_ORIVX_COV_, &b_Ds_ORIVX_COV_);
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
	fChain->SetBranchAddress("Ds_ptasy_1.00", &Ds_ptasy_1_00, &b_Ds_ptasy_1_00);
	fChain->SetBranchAddress("pi_plus_fromDs_ETA", &pi_plus_fromDs_ETA, &b_pi_plus_fromDs_ETA);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNe", &pi_plus_fromDs_MC12TuneV2_ProbNNe, &b_pi_plus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNmu", &pi_plus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNpi", &pi_plus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNk", &pi_plus_fromDs_MC12TuneV2_ProbNNk, &b_pi_plus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNp", &pi_plus_fromDs_MC12TuneV2_ProbNNp, &b_pi_plus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV2_ProbNNghost", &pi_plus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNe", &pi_plus_fromDs_MC12TuneV3_ProbNNe, &b_pi_plus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNmu", &pi_plus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_plus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNpi", &pi_plus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_plus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNk", &pi_plus_fromDs_MC12TuneV3_ProbNNk, &b_pi_plus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNp", &pi_plus_fromDs_MC12TuneV3_ProbNNp, &b_pi_plus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_plus_fromDs_MC12TuneV3_ProbNNghost", &pi_plus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_plus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_fromDs_IP_OWNPV", &pi_plus_fromDs_IP_OWNPV, &b_pi_plus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("pi_plus_fromDs_IPCHI2_OWNPV", &pi_plus_fromDs_IPCHI2_OWNPV, &b_pi_plus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_plus_fromDs_P", &pi_plus_fromDs_P, &b_pi_plus_fromDs_P);
	fChain->SetBranchAddress("pi_plus_fromDs_PT", &pi_plus_fromDs_PT, &b_pi_plus_fromDs_PT);
	fChain->SetBranchAddress("pi_plus_fromDs_PE", &pi_plus_fromDs_PE, &b_pi_plus_fromDs_PE);
	fChain->SetBranchAddress("pi_plus_fromDs_PX", &pi_plus_fromDs_PX, &b_pi_plus_fromDs_PX);
	fChain->SetBranchAddress("pi_plus_fromDs_PY", &pi_plus_fromDs_PY, &b_pi_plus_fromDs_PY);
	fChain->SetBranchAddress("pi_plus_fromDs_PZ", &pi_plus_fromDs_PZ, &b_pi_plus_fromDs_PZ);
	fChain->SetBranchAddress("pi_plus_fromDs_ID", &pi_plus_fromDs_ID, &b_pi_plus_fromDs_ID);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDmu", &pi_plus_fromDs_PIDmu, &b_pi_plus_fromDs_PIDmu);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDK", &pi_plus_fromDs_PIDK, &b_pi_plus_fromDs_PIDK);
	fChain->SetBranchAddress("pi_plus_fromDs_PIDp", &pi_plus_fromDs_PIDp, &b_pi_plus_fromDs_PIDp);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNk", &pi_plus_fromDs_ProbNNk, &b_pi_plus_fromDs_ProbNNk);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNp", &pi_plus_fromDs_ProbNNp, &b_pi_plus_fromDs_ProbNNp);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNpi", &pi_plus_fromDs_ProbNNpi, &b_pi_plus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNmu", &pi_plus_fromDs_ProbNNmu, &b_pi_plus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("pi_plus_fromDs_ProbNNghost", &pi_plus_fromDs_ProbNNghost, &b_pi_plus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("pi_plus_fromDs_isMuon", &pi_plus_fromDs_isMuon, &b_pi_plus_fromDs_isMuon);
	fChain->SetBranchAddress("pi_plus_fromDs_TRACK_GhostProb", &pi_plus_fromDs_TRACK_GhostProb, &b_pi_plus_fromDs_TRACK_GhostProb);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA1", &pi_minus_fromDs_DOCA1, &b_pi_minus_fromDs_DOCA1);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA2", &pi_minus_fromDs_DOCA2, &b_pi_minus_fromDs_DOCA2);
// 	fChain->SetBranchAddress("pi_minus_fromDs_DOCA3", &pi_minus_fromDs_DOCA3, &b_pi_minus_fromDs_DOCA3);
	fChain->SetBranchAddress("pi_minus_fromDs_ETA", &pi_minus_fromDs_ETA, &b_pi_minus_fromDs_ETA);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNe", &pi_minus_fromDs_MC12TuneV2_ProbNNe, &b_pi_minus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNmu", &pi_minus_fromDs_MC12TuneV2_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNpi", &pi_minus_fromDs_MC12TuneV2_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNk", &pi_minus_fromDs_MC12TuneV2_ProbNNk, &b_pi_minus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNp", &pi_minus_fromDs_MC12TuneV2_ProbNNp, &b_pi_minus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV2_ProbNNghost", &pi_minus_fromDs_MC12TuneV2_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNe", &pi_minus_fromDs_MC12TuneV3_ProbNNe, &b_pi_minus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNmu", &pi_minus_fromDs_MC12TuneV3_ProbNNmu, &b_pi_minus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNpi", &pi_minus_fromDs_MC12TuneV3_ProbNNpi, &b_pi_minus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNk", &pi_minus_fromDs_MC12TuneV3_ProbNNk, &b_pi_minus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNp", &pi_minus_fromDs_MC12TuneV3_ProbNNp, &b_pi_minus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_fromDs_MC12TuneV3_ProbNNghost", &pi_minus_fromDs_MC12TuneV3_ProbNNghost, &b_pi_minus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_fromDs_IP_OWNPV", &pi_minus_fromDs_IP_OWNPV, &b_pi_minus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("pi_minus_fromDs_IPCHI2_OWNPV", &pi_minus_fromDs_IPCHI2_OWNPV, &b_pi_minus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_minus_fromDs_P", &pi_minus_fromDs_P, &b_pi_minus_fromDs_P);
	fChain->SetBranchAddress("pi_minus_fromDs_PT", &pi_minus_fromDs_PT, &b_pi_minus_fromDs_PT);
	fChain->SetBranchAddress("pi_minus_fromDs_PE", &pi_minus_fromDs_PE, &b_pi_minus_fromDs_PE);
	fChain->SetBranchAddress("pi_minus_fromDs_PX", &pi_minus_fromDs_PX, &b_pi_minus_fromDs_PX);
	fChain->SetBranchAddress("pi_minus_fromDs_PY", &pi_minus_fromDs_PY, &b_pi_minus_fromDs_PY);
	fChain->SetBranchAddress("pi_minus_fromDs_PZ", &pi_minus_fromDs_PZ, &b_pi_minus_fromDs_PZ);
	fChain->SetBranchAddress("pi_minus_fromDs_ID", &pi_minus_fromDs_ID, &b_pi_minus_fromDs_ID);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDmu", &pi_minus_fromDs_PIDmu, &b_pi_minus_fromDs_PIDmu);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDK", &pi_minus_fromDs_PIDK, &b_pi_minus_fromDs_PIDK);
	fChain->SetBranchAddress("pi_minus_fromDs_PIDp", &pi_minus_fromDs_PIDp, &b_pi_minus_fromDs_PIDp);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNk", &pi_minus_fromDs_ProbNNk, &b_pi_minus_fromDs_ProbNNk);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNp", &pi_minus_fromDs_ProbNNp, &b_pi_minus_fromDs_ProbNNp);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNpi", &pi_minus_fromDs_ProbNNpi, &b_pi_minus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNmu", &pi_minus_fromDs_ProbNNmu, &b_pi_minus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("pi_minus_fromDs_ProbNNghost", &pi_minus_fromDs_ProbNNghost, &b_pi_minus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_fromDs_isMuon", &pi_minus_fromDs_isMuon, &b_pi_minus_fromDs_isMuon);
	fChain->SetBranchAddress("pi_minus_fromDs_TRACK_GhostProb", &pi_minus_fromDs_TRACK_GhostProb, &b_pi_minus_fromDs_TRACK_GhostProb);
	fChain->SetBranchAddress("K_minus_fromDs_ETA", &K_minus_fromDs_ETA, &b_K_minus_fromDs_ETA);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNe", &K_minus_fromDs_MC12TuneV2_ProbNNe, &b_K_minus_fromDs_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNmu", &K_minus_fromDs_MC12TuneV2_ProbNNmu, &b_K_minus_fromDs_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNpi", &K_minus_fromDs_MC12TuneV2_ProbNNpi, &b_K_minus_fromDs_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNk", &K_minus_fromDs_MC12TuneV2_ProbNNk, &b_K_minus_fromDs_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNp", &K_minus_fromDs_MC12TuneV2_ProbNNp, &b_K_minus_fromDs_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV2_ProbNNghost", &K_minus_fromDs_MC12TuneV2_ProbNNghost, &b_K_minus_fromDs_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNe", &K_minus_fromDs_MC12TuneV3_ProbNNe, &b_K_minus_fromDs_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNmu", &K_minus_fromDs_MC12TuneV3_ProbNNmu, &b_K_minus_fromDs_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNpi", &K_minus_fromDs_MC12TuneV3_ProbNNpi, &b_K_minus_fromDs_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNk", &K_minus_fromDs_MC12TuneV3_ProbNNk, &b_K_minus_fromDs_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNp", &K_minus_fromDs_MC12TuneV3_ProbNNp, &b_K_minus_fromDs_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("K_minus_fromDs_MC12TuneV3_ProbNNghost", &K_minus_fromDs_MC12TuneV3_ProbNNghost, &b_K_minus_fromDs_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("K_minus_fromDs_IP_OWNPV", &K_minus_fromDs_IP_OWNPV, &b_K_minus_fromDs_IP_OWNPV);
	fChain->SetBranchAddress("K_minus_fromDs_IPCHI2_OWNPV", &K_minus_fromDs_IPCHI2_OWNPV, &b_K_minus_fromDs_IPCHI2_OWNPV);
	fChain->SetBranchAddress("K_minus_fromDs_P", &K_minus_fromDs_P, &b_K_minus_fromDs_P);
	fChain->SetBranchAddress("K_minus_fromDs_PT", &K_minus_fromDs_PT, &b_K_minus_fromDs_PT);
	fChain->SetBranchAddress("K_minus_fromDs_PE", &K_minus_fromDs_PE, &b_K_minus_fromDs_PE);
	fChain->SetBranchAddress("K_minus_fromDs_PX", &K_minus_fromDs_PX, &b_K_minus_fromDs_PX);
	fChain->SetBranchAddress("K_minus_fromDs_PY", &K_minus_fromDs_PY, &b_K_minus_fromDs_PY);
	fChain->SetBranchAddress("K_minus_fromDs_PZ", &K_minus_fromDs_PZ, &b_K_minus_fromDs_PZ);
	fChain->SetBranchAddress("K_minus_fromDs_ID", &K_minus_fromDs_ID, &b_K_minus_fromDs_ID);
	fChain->SetBranchAddress("K_minus_fromDs_PIDmu", &K_minus_fromDs_PIDmu, &b_K_minus_fromDs_PIDmu);
	fChain->SetBranchAddress("K_minus_fromDs_PIDK", &K_minus_fromDs_PIDK, &b_K_minus_fromDs_PIDK);
	fChain->SetBranchAddress("K_minus_fromDs_PIDp", &K_minus_fromDs_PIDp, &b_K_minus_fromDs_PIDp);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNk", &K_minus_fromDs_ProbNNk, &b_K_minus_fromDs_ProbNNk);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNp", &K_minus_fromDs_ProbNNp, &b_K_minus_fromDs_ProbNNp);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNpi", &K_minus_fromDs_ProbNNpi, &b_K_minus_fromDs_ProbNNpi);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNmu", &K_minus_fromDs_ProbNNmu, &b_K_minus_fromDs_ProbNNmu);
	fChain->SetBranchAddress("K_minus_fromDs_ProbNNghost", &K_minus_fromDs_ProbNNghost, &b_K_minus_fromDs_ProbNNghost);
	fChain->SetBranchAddress("K_minus_fromDs_isMuon", &K_minus_fromDs_isMuon, &b_K_minus_fromDs_isMuon);
// 	fChain->SetBranchAddress("K_minus_fromDs_hasRich", &K_minus_fromDs_hasRich, &b_K_minus_fromDs_hasRich);
// 	fChain->SetBranchAddress("K_minus_fromDs_hasCalo", &K_minus_fromDs_hasCalo, &b_K_minus_fromDs_hasCalo);
	fChain->SetBranchAddress("K_minus_fromDs_TRACK_GhostProb", &K_minus_fromDs_TRACK_GhostProb, &b_K_minus_fromDs_TRACK_GhostProb);
	fChain->SetBranchAddress("a_1_1260_plus_DOCA1", &a_1_1260_plus_DOCA1, &b_a_1_1260_plus_DOCA1);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA2", &a_1_1260_plus_DOCA2, &b_a_1_1260_plus_DOCA2);
        fChain->SetBranchAddress("a_1_1260_plus_DOCA3", &a_1_1260_plus_DOCA3, &b_a_1_1260_plus_DOCA3);
        fChain->SetBranchAddress("a_1_1260_plus_ETA", &a_1_1260_plus_ETA, &b_a_1_1260_plus_ETA);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_X", &a_1_1260_plus_ENDVERTEX_X, &b_a_1_1260_plus_ENDVERTEX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Y", &a_1_1260_plus_ENDVERTEX_Y, &b_a_1_1260_plus_ENDVERTEX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_Z", &a_1_1260_plus_ENDVERTEX_Z, &b_a_1_1260_plus_ENDVERTEX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_XERR", &a_1_1260_plus_ENDVERTEX_XERR, &b_a_1_1260_plus_ENDVERTEX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_YERR", &a_1_1260_plus_ENDVERTEX_YERR, &b_a_1_1260_plus_ENDVERTEX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_ZERR", &a_1_1260_plus_ENDVERTEX_ZERR, &b_a_1_1260_plus_ENDVERTEX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_CHI2", &a_1_1260_plus_ENDVERTEX_CHI2, &b_a_1_1260_plus_ENDVERTEX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ENDVERTEX_NDOF", &a_1_1260_plus_ENDVERTEX_NDOF, &b_a_1_1260_plus_ENDVERTEX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_X", &a_1_1260_plus_OWNPV_X, &b_a_1_1260_plus_OWNPV_X);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Y", &a_1_1260_plus_OWNPV_Y, &b_a_1_1260_plus_OWNPV_Y);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_Z", &a_1_1260_plus_OWNPV_Z, &b_a_1_1260_plus_OWNPV_Z);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_XERR", &a_1_1260_plus_OWNPV_XERR, &b_a_1_1260_plus_OWNPV_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_YERR", &a_1_1260_plus_OWNPV_YERR, &b_a_1_1260_plus_OWNPV_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_ZERR", &a_1_1260_plus_OWNPV_ZERR, &b_a_1_1260_plus_OWNPV_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_CHI2", &a_1_1260_plus_OWNPV_CHI2, &b_a_1_1260_plus_OWNPV_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_OWNPV_NDOF", &a_1_1260_plus_OWNPV_NDOF, &b_a_1_1260_plus_OWNPV_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_IP_OWNPV", &a_1_1260_plus_IP_OWNPV, &b_a_1_1260_plus_IP_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_IPCHI2_OWNPV", &a_1_1260_plus_IPCHI2_OWNPV, &b_a_1_1260_plus_IPCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FD_OWNPV", &a_1_1260_plus_FD_OWNPV, &b_a_1_1260_plus_FD_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_OWNPV", &a_1_1260_plus_FDCHI2_OWNPV, &b_a_1_1260_plus_FDCHI2_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_OWNPV", &a_1_1260_plus_DIRA_OWNPV, &b_a_1_1260_plus_DIRA_OWNPV);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_X", &a_1_1260_plus_ORIVX_X, &b_a_1_1260_plus_ORIVX_X);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Y", &a_1_1260_plus_ORIVX_Y, &b_a_1_1260_plus_ORIVX_Y);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_Z", &a_1_1260_plus_ORIVX_Z, &b_a_1_1260_plus_ORIVX_Z);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_XERR", &a_1_1260_plus_ORIVX_XERR, &b_a_1_1260_plus_ORIVX_XERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_YERR", &a_1_1260_plus_ORIVX_YERR, &b_a_1_1260_plus_ORIVX_YERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_ZERR", &a_1_1260_plus_ORIVX_ZERR, &b_a_1_1260_plus_ORIVX_ZERR);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_CHI2", &a_1_1260_plus_ORIVX_CHI2, &b_a_1_1260_plus_ORIVX_CHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ORIVX_NDOF", &a_1_1260_plus_ORIVX_NDOF, &b_a_1_1260_plus_ORIVX_NDOF);
        fChain->SetBranchAddress("a_1_1260_plus_FD_ORIVX", &a_1_1260_plus_FD_ORIVX, &b_a_1_1260_plus_FD_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_FDCHI2_ORIVX", &a_1_1260_plus_FDCHI2_ORIVX, &b_a_1_1260_plus_FDCHI2_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_DIRA_ORIVX", &a_1_1260_plus_DIRA_ORIVX, &b_a_1_1260_plus_DIRA_ORIVX);
        fChain->SetBranchAddress("a_1_1260_plus_P", &a_1_1260_plus_P, &b_a_1_1260_plus_P);
        fChain->SetBranchAddress("a_1_1260_plus_PT", &a_1_1260_plus_PT, &b_a_1_1260_plus_PT);
        fChain->SetBranchAddress("a_1_1260_plus_PE", &a_1_1260_plus_PE, &b_a_1_1260_plus_PE);
        fChain->SetBranchAddress("a_1_1260_plus_PX", &a_1_1260_plus_PX, &b_a_1_1260_plus_PX);
        fChain->SetBranchAddress("a_1_1260_plus_PY", &a_1_1260_plus_PY, &b_a_1_1260_plus_PY);
        fChain->SetBranchAddress("a_1_1260_plus_PZ", &a_1_1260_plus_PZ, &b_a_1_1260_plus_PZ);
        fChain->SetBranchAddress("a_1_1260_plus_MM", &a_1_1260_plus_MM, &b_a_1_1260_plus_MM);
        fChain->SetBranchAddress("a_1_1260_plus_MMERR", &a_1_1260_plus_MMERR, &b_a_1_1260_plus_MMERR);
        fChain->SetBranchAddress("a_1_1260_plus_ID", &a_1_1260_plus_ID, &b_a_1_1260_plus_ID);
/*        fChain->SetBranchAddress("a_1_1260_plus_TAU", &a_1_1260_plus_TAU, &b_a_1_1260_plus_TAU);
        fChain->SetBranchAddress("a_1_1260_plus_TAUERR", &a_1_1260_plus_TAUERR, &b_a_1_1260_plus_TAUERR);*/
        //fChain->SetBranchAddress("a_1_1260_plus_TAUCHI2", &a_1_1260_plus_TAUCHI2, &b_a_1_1260_plus_TAUCHI2);
        fChain->SetBranchAddress("a_1_1260_plus_ptasy_1.00", &a_1_1260_plus_ptasy_1_00, &b_a_1_1260_plus_ptasy_1_00);
        fChain->SetBranchAddress("pi_plus1_ETA", &pi_plus1_ETA, &b_pi_plus1_ETA);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNe", &K_plus_MC12TuneV2_ProbNNe, &b_K_plus_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNmu", &K_plus_MC12TuneV2_ProbNNmu, &b_K_plus_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNpi", &K_plus_MC12TuneV2_ProbNNpi, &b_K_plus_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNk", &K_plus_MC12TuneV2_ProbNNk, &b_K_plus_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNp", &K_plus_MC12TuneV2_ProbNNp, &b_K_plus_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV2_ProbNNghost", &K_plus_MC12TuneV2_ProbNNghost, &b_K_plus_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNe", &K_plus_MC12TuneV3_ProbNNe, &b_K_plus_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNmu", &K_plus_MC12TuneV3_ProbNNmu, &b_K_plus_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNpi", &K_plus_MC12TuneV3_ProbNNpi, &b_K_plus_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNk", &K_plus_MC12TuneV3_ProbNNk, &b_K_plus_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNp", &K_plus_MC12TuneV3_ProbNNp, &b_K_plus_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("K_plus_MC12TuneV3_ProbNNghost", &K_plus_MC12TuneV3_ProbNNghost, &b_K_plus_MC12TuneV3_ProbNNghost);
	    fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNmu", &pi_plus1_MC12TuneV2_ProbNNmu, &b_pi_plus1_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNpi", &pi_plus1_MC12TuneV2_ProbNNpi, &b_pi_plus1_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNk", &pi_plus1_MC12TuneV2_ProbNNk, &b_pi_plus1_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNp", &pi_plus1_MC12TuneV2_ProbNNp, &b_pi_plus1_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV2_ProbNNghost", &pi_plus1_MC12TuneV2_ProbNNghost, &b_pi_plus1_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNmu", &pi_plus1_MC12TuneV3_ProbNNmu, &b_pi_plus1_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNpi", &pi_plus1_MC12TuneV3_ProbNNpi, &b_pi_plus1_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNk", &pi_plus1_MC12TuneV3_ProbNNk, &b_pi_plus1_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNp", &pi_plus1_MC12TuneV3_ProbNNp, &b_pi_plus1_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_MC12TuneV3_ProbNNghost", &pi_plus1_MC12TuneV3_ProbNNghost, &b_pi_plus1_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_IP_OWNPV", &pi_plus1_IP_OWNPV, &b_pi_plus1_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus1_IPCHI2_OWNPV", &pi_plus1_IPCHI2_OWNPV, &b_pi_plus1_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus1_P", &pi_plus1_P, &b_pi_plus1_P);
        fChain->SetBranchAddress("pi_plus1_PT", &pi_plus1_PT, &b_pi_plus1_PT);
        fChain->SetBranchAddress("pi_plus1_PE", &pi_plus1_PE, &b_pi_plus1_PE);
        fChain->SetBranchAddress("pi_plus1_PX", &pi_plus1_PX, &b_pi_plus1_PX);
        fChain->SetBranchAddress("pi_plus1_PY", &pi_plus1_PY, &b_pi_plus1_PY);
        fChain->SetBranchAddress("pi_plus1_PZ", &pi_plus1_PZ, &b_pi_plus1_PZ);
        fChain->SetBranchAddress("pi_plus1_ID", &pi_plus1_ID, &b_pi_plus1_ID);
        fChain->SetBranchAddress("pi_plus1_PIDmu", &pi_plus1_PIDmu, &b_pi_plus1_PIDmu);
        fChain->SetBranchAddress("pi_plus1_PIDK", &pi_plus1_PIDK, &b_pi_plus1_PIDK);
        fChain->SetBranchAddress("pi_plus1_PIDp", &pi_plus1_PIDp, &b_pi_plus1_PIDp);
        fChain->SetBranchAddress("pi_plus1_ProbNNk", &pi_plus1_ProbNNk, &b_pi_plus1_ProbNNk);
        fChain->SetBranchAddress("pi_plus1_ProbNNp", &pi_plus1_ProbNNp, &b_pi_plus1_ProbNNp);
        fChain->SetBranchAddress("pi_plus1_ProbNNpi", &pi_plus1_ProbNNpi, &b_pi_plus1_ProbNNpi);
        fChain->SetBranchAddress("pi_plus1_ProbNNmu", &pi_plus1_ProbNNmu, &b_pi_plus1_ProbNNmu);
        fChain->SetBranchAddress("pi_plus1_ProbNNghost", &pi_plus1_ProbNNghost, &b_pi_plus1_ProbNNghost);
        fChain->SetBranchAddress("pi_plus1_isMuon", &pi_plus1_isMuon, &b_pi_plus1_isMuon);
        fChain->SetBranchAddress("pi_plus1_TRACK_CHI2NDOF", &pi_plus1_TRACK_CHI2NDOF, &b_pi_plus1_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus1_TRACK_GhostProb", &pi_plus1_TRACK_GhostProb, &b_pi_plus1_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus1_ptasy_1.00", &pi_plus1_ptasy_1_00, &b_pi_plus1_ptasy_1_00);
	 fChain->SetBranchAddress("pi_plus2_ETA", &pi_plus2_ETA, &b_pi_plus2_ETA);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNmu", &pi_plus2_MC12TuneV2_ProbNNmu, &b_pi_plus2_MC12TuneV2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNpi", &pi_plus2_MC12TuneV2_ProbNNpi, &b_pi_plus2_MC12TuneV2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNk", &pi_plus2_MC12TuneV2_ProbNNk, &b_pi_plus2_MC12TuneV2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNp", &pi_plus2_MC12TuneV2_ProbNNp, &b_pi_plus2_MC12TuneV2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV2_ProbNNghost", &pi_plus2_MC12TuneV2_ProbNNghost, &b_pi_plus2_MC12TuneV2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNmu", &pi_plus2_MC12TuneV3_ProbNNmu, &b_pi_plus2_MC12TuneV3_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNpi", &pi_plus2_MC12TuneV3_ProbNNpi, &b_pi_plus2_MC12TuneV3_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNk", &pi_plus2_MC12TuneV3_ProbNNk, &b_pi_plus2_MC12TuneV3_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNp", &pi_plus2_MC12TuneV3_ProbNNp, &b_pi_plus2_MC12TuneV3_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_MC12TuneV3_ProbNNghost", &pi_plus2_MC12TuneV3_ProbNNghost, &b_pi_plus2_MC12TuneV3_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_IP_OWNPV", &pi_plus2_IP_OWNPV, &b_pi_plus2_IP_OWNPV);
        fChain->SetBranchAddress("pi_plus2_IPCHI2_OWNPV", &pi_plus2_IPCHI2_OWNPV, &b_pi_plus2_IPCHI2_OWNPV);
        fChain->SetBranchAddress("pi_plus2_P", &pi_plus2_P, &b_pi_plus2_P);
        fChain->SetBranchAddress("pi_plus2_PT", &pi_plus2_PT, &b_pi_plus2_PT);
        fChain->SetBranchAddress("pi_plus2_PE", &pi_plus2_PE, &b_pi_plus2_PE);
        fChain->SetBranchAddress("pi_plus2_PX", &pi_plus2_PX, &b_pi_plus2_PX);
        fChain->SetBranchAddress("pi_plus2_PY", &pi_plus2_PY, &b_pi_plus2_PY);
        fChain->SetBranchAddress("pi_plus2_PZ", &pi_plus2_PZ, &b_pi_plus2_PZ);
        fChain->SetBranchAddress("pi_plus2_ID", &pi_plus2_ID, &b_pi_plus2_ID);
        fChain->SetBranchAddress("pi_plus2_PIDmu", &pi_plus2_PIDmu, &b_pi_plus2_PIDmu);
        fChain->SetBranchAddress("pi_plus2_PIDK", &pi_plus2_PIDK, &b_pi_plus2_PIDK);
        fChain->SetBranchAddress("pi_plus2_PIDp", &pi_plus2_PIDp, &b_pi_plus2_PIDp);
        fChain->SetBranchAddress("pi_plus2_ProbNNk", &pi_plus2_ProbNNk, &b_pi_plus2_ProbNNk);
        fChain->SetBranchAddress("pi_plus2_ProbNNp", &pi_plus2_ProbNNp, &b_pi_plus2_ProbNNp);
        fChain->SetBranchAddress("pi_plus2_ProbNNpi", &pi_plus2_ProbNNpi, &b_pi_plus2_ProbNNpi);
        fChain->SetBranchAddress("pi_plus2_ProbNNmu", &pi_plus2_ProbNNmu, &b_pi_plus2_ProbNNmu);
        fChain->SetBranchAddress("pi_plus2_ProbNNghost", &pi_plus2_ProbNNghost, &b_pi_plus2_ProbNNghost);
        fChain->SetBranchAddress("pi_plus2_isMuon", &pi_plus2_isMuon, &b_pi_plus2_isMuon);
        fChain->SetBranchAddress("pi_plus2_TRACK_CHI2NDOF", &pi_plus2_TRACK_CHI2NDOF, &b_pi_plus2_TRACK_CHI2NDOF);
        fChain->SetBranchAddress("pi_plus2_TRACK_GhostProb", &pi_plus2_TRACK_GhostProb, &b_pi_plus2_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_plus2_ptasy_1.00", &pi_plus2_ptasy_1_00, &b_pi_plus2_ptasy_1_00);
	fChain->SetBranchAddress("pi_minus_ETA", &pi_minus_ETA, &b_pi_minus_ETA);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNe", &pi_minus_MC12TuneV2_ProbNNe, &b_pi_minus_MC12TuneV2_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNmu", &pi_minus_MC12TuneV2_ProbNNmu, &b_pi_minus_MC12TuneV2_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNpi", &pi_minus_MC12TuneV2_ProbNNpi, &b_pi_minus_MC12TuneV2_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNk", &pi_minus_MC12TuneV2_ProbNNk, &b_pi_minus_MC12TuneV2_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNp", &pi_minus_MC12TuneV2_ProbNNp, &b_pi_minus_MC12TuneV2_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV2_ProbNNghost", &pi_minus_MC12TuneV2_ProbNNghost, &b_pi_minus_MC12TuneV2_ProbNNghost);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNe", &pi_minus_MC12TuneV3_ProbNNe, &b_pi_minus_MC12TuneV3_ProbNNe);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNmu", &pi_minus_MC12TuneV3_ProbNNmu, &b_pi_minus_MC12TuneV3_ProbNNmu);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNpi", &pi_minus_MC12TuneV3_ProbNNpi, &b_pi_minus_MC12TuneV3_ProbNNpi);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNk", &pi_minus_MC12TuneV3_ProbNNk, &b_pi_minus_MC12TuneV3_ProbNNk);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNp", &pi_minus_MC12TuneV3_ProbNNp, &b_pi_minus_MC12TuneV3_ProbNNp);
	//    fChain->SetBranchAddress("pi_minus_MC12TuneV3_ProbNNghost", &pi_minus_MC12TuneV3_ProbNNghost, &b_pi_minus_MC12TuneV3_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_IP_OWNPV", &pi_minus_IP_OWNPV, &b_pi_minus_IP_OWNPV);
	fChain->SetBranchAddress("pi_minus_IPCHI2_OWNPV", &pi_minus_IPCHI2_OWNPV, &b_pi_minus_IPCHI2_OWNPV);
	fChain->SetBranchAddress("pi_minus_P", &pi_minus_P, &b_pi_minus_P);
	fChain->SetBranchAddress("pi_minus_PT", &pi_minus_PT, &b_pi_minus_PT);
	fChain->SetBranchAddress("pi_minus_PE", &pi_minus_PE, &b_pi_minus_PE);
	fChain->SetBranchAddress("pi_minus_PX", &pi_minus_PX, &b_pi_minus_PX);
	fChain->SetBranchAddress("pi_minus_PY", &pi_minus_PY, &b_pi_minus_PY);
	fChain->SetBranchAddress("pi_minus_PZ", &pi_minus_PZ, &b_pi_minus_PZ);
	fChain->SetBranchAddress("pi_minus_ID", &pi_minus_ID, &b_pi_minus_ID);
	fChain->SetBranchAddress("pi_minus_PIDmu", &pi_minus_PIDmu, &b_pi_minus_PIDmu);
	fChain->SetBranchAddress("pi_minus_PIDK", &pi_minus_PIDK, &b_pi_minus_PIDK);
	fChain->SetBranchAddress("pi_minus_PIDp", &pi_minus_PIDp, &b_pi_minus_PIDp);
	fChain->SetBranchAddress("pi_minus_ProbNNk", &pi_minus_ProbNNk, &b_pi_minus_ProbNNk);
	fChain->SetBranchAddress("pi_minus_ProbNNp", &pi_minus_ProbNNp, &b_pi_minus_ProbNNp);
	fChain->SetBranchAddress("pi_minus_ProbNNpi", &pi_minus_ProbNNpi, &b_pi_minus_ProbNNpi);
	fChain->SetBranchAddress("pi_minus_ProbNNmu", &pi_minus_ProbNNmu, &b_pi_minus_ProbNNmu);
	fChain->SetBranchAddress("pi_minus_ProbNNghost", &pi_minus_ProbNNghost, &b_pi_minus_ProbNNghost);
	fChain->SetBranchAddress("pi_minus_isMuon", &pi_minus_isMuon, &b_pi_minus_isMuon);
	fChain->SetBranchAddress("pi_minus_TRACK_GhostProb", &pi_minus_TRACK_GhostProb, &b_pi_minus_TRACK_GhostProb);
        fChain->SetBranchAddress("pi_minus_ptasy_1.00", &pi_minus_ptasy_1_00, &b_pi_minus_ptasy_1_00);
	fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
	fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
	fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
	//    fChain->SetBranchAddress("nLong", &nLong, &b_nLong);
	fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
	fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
	fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
	fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
    }
	
	

    Notify();
}


#endif // #ifdef MiniDecayTree_cxx
