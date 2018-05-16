//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 27 11:31:06 2017 by ROOT version 5.34/36
// from TTree MCDecayTree/MCDecayTree
// found on file: DVntuple2.root
//////////////////////////////////////////////////////////

#ifndef MCDecayTreeTuple_h
#define MCDecayTreeTuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

static const double massKaon = 493.68;
static const double massPion = 139.57;

class MCDecayTreeTuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
    Double_t        B_s0_nPhotos;
    Double_t        B_s0_TRUEID;
    Int_t           B_s0_MC_MOTHER_ID;
    Int_t           B_s0_MC_MOTHER_KEY;
    Int_t           B_s0_MC_GD_MOTHER_ID;
    Int_t           B_s0_MC_GD_MOTHER_KEY;
    Int_t           B_s0_MC_GD_GD_MOTHER_ID;
    Int_t           B_s0_MC_GD_GD_MOTHER_KEY;
    Double_t        B_s0_TRUEP_E;
    Double_t        B_s0_TRUEP_X;
    Double_t        B_s0_TRUEP_Y;
    Double_t        B_s0_TRUEP_Z;
    Double_t        B_s0_TRUEPT;
    Double_t        B_s0_TRUEORIGINVERTEX_X;
    Double_t        B_s0_TRUEORIGINVERTEX_Y;
    Double_t        B_s0_TRUEORIGINVERTEX_Z;
    Double_t        B_s0_TRUEENDVERTEX_X;
    Double_t        B_s0_TRUEENDVERTEX_Y;
    Double_t        B_s0_TRUEENDVERTEX_Z;
    Bool_t          B_s0_TRUEISSTABLE;
    Double_t        B_s0_TRUETAU;
    Bool_t          B_s0_OSCIL;
    Double_t        Kplus_nPhotos;
    Double_t        Kplus_TRUEID;
    Int_t           Kplus_MC_MOTHER_ID;
    Int_t           Kplus_MC_MOTHER_KEY;
    Int_t           Kplus_MC_GD_MOTHER_ID;
    Int_t           Kplus_MC_GD_MOTHER_KEY;
    Int_t           Kplus_MC_GD_GD_MOTHER_ID;
    Int_t           Kplus_MC_GD_GD_MOTHER_KEY;
    Double_t        Kplus_TRUEP_E;
    Double_t        Kplus_TRUEP_X;
    Double_t        Kplus_TRUEP_Y;
    Double_t        Kplus_TRUEP_Z;
    Double_t        Kplus_TRUEPT;
    Double_t        Kplus_TRUEORIGINVERTEX_X;
    Double_t        Kplus_TRUEORIGINVERTEX_Y;
    Double_t        Kplus_TRUEORIGINVERTEX_Z;
    Double_t        Kplus_TRUEENDVERTEX_X;
    Double_t        Kplus_TRUEENDVERTEX_Y;
    Double_t        Kplus_TRUEENDVERTEX_Z;
    Bool_t          Kplus_TRUEISSTABLE;
    Double_t        Kplus_TRUETAU;
    Bool_t          Kplus_OSCIL;
    Double_t        piplus_nPhotos;
    Double_t        piplus_TRUEID;
    Int_t           piplus_MC_MOTHER_ID;
    Int_t           piplus_MC_MOTHER_KEY;
    Int_t           piplus_MC_GD_MOTHER_ID;
    Int_t           piplus_MC_GD_MOTHER_KEY;
    Int_t           piplus_MC_GD_GD_MOTHER_ID;
    Int_t           piplus_MC_GD_GD_MOTHER_KEY;
    Double_t        piplus_TRUEP_E;
    Double_t        piplus_TRUEP_X;
    Double_t        piplus_TRUEP_Y;
    Double_t        piplus_TRUEP_Z;
    Double_t        piplus_TRUEPT;
    Double_t        piplus_TRUEORIGINVERTEX_X;
    Double_t        piplus_TRUEORIGINVERTEX_Y;
    Double_t        piplus_TRUEORIGINVERTEX_Z;
    Double_t        piplus_TRUEENDVERTEX_X;
    Double_t        piplus_TRUEENDVERTEX_Y;
    Double_t        piplus_TRUEENDVERTEX_Z;
    Bool_t          piplus_TRUEISSTABLE;
    Double_t        piplus_TRUETAU;
    Bool_t          piplus_OSCIL;
    Double_t        piminus_nPhotos;
    Double_t        piminus_TRUEID;
    Int_t           piminus_MC_MOTHER_ID;
    Int_t           piminus_MC_MOTHER_KEY;
    Int_t           piminus_MC_GD_MOTHER_ID;
    Int_t           piminus_MC_GD_MOTHER_KEY;
    Int_t           piminus_MC_GD_GD_MOTHER_ID;
    Int_t           piminus_MC_GD_GD_MOTHER_KEY;
    Double_t        piminus_TRUEP_E;
    Double_t        piminus_TRUEP_X;
    Double_t        piminus_TRUEP_Y;
    Double_t        piminus_TRUEP_Z;
    Double_t        piminus_TRUEPT;
    Double_t        piminus_TRUEORIGINVERTEX_X;
    Double_t        piminus_TRUEORIGINVERTEX_Y;
    Double_t        piminus_TRUEORIGINVERTEX_Z;
    Double_t        piminus_TRUEENDVERTEX_X;
    Double_t        piminus_TRUEENDVERTEX_Y;
    Double_t        piminus_TRUEENDVERTEX_Z;
    Bool_t          piminus_TRUEISSTABLE;
    Double_t        piminus_TRUETAU;
    Bool_t          piminus_OSCIL;
    Double_t        D_sminus_nPhotos;
    Double_t        D_sminus_TRUEID;
    Int_t           D_sminus_MC_MOTHER_ID;
    Int_t           D_sminus_MC_MOTHER_KEY;
    Int_t           D_sminus_MC_GD_MOTHER_ID;
    Int_t           D_sminus_MC_GD_MOTHER_KEY;
    Int_t           D_sminus_MC_GD_GD_MOTHER_ID;
    Int_t           D_sminus_MC_GD_GD_MOTHER_KEY;
    Double_t        D_sminus_TRUEP_E;
    Double_t        D_sminus_TRUEP_X;
    Double_t        D_sminus_TRUEP_Y;
    Double_t        D_sminus_TRUEP_Z;
    Double_t        D_sminus_TRUEPT;
    Double_t        D_sminus_TRUEORIGINVERTEX_X;
    Double_t        D_sminus_TRUEORIGINVERTEX_Y;
    Double_t        D_sminus_TRUEORIGINVERTEX_Z;
    Double_t        D_sminus_TRUEENDVERTEX_X;
    Double_t        D_sminus_TRUEENDVERTEX_Y;
    Double_t        D_sminus_TRUEENDVERTEX_Z;
    Bool_t          D_sminus_TRUEISSTABLE;
    Double_t        D_sminus_TRUETAU;
    Bool_t          D_sminus_OSCIL;
    Double_t        Kplus0_nPhotos;
    Double_t        Kplus0_TRUEID;
    Int_t           Kplus0_MC_MOTHER_ID;
    Int_t           Kplus0_MC_MOTHER_KEY;
    Int_t           Kplus0_MC_GD_MOTHER_ID;
    Int_t           Kplus0_MC_GD_MOTHER_KEY;
    Int_t           Kplus0_MC_GD_GD_MOTHER_ID;
    Int_t           Kplus0_MC_GD_GD_MOTHER_KEY;
    Double_t        Kplus0_TRUEP_E;
    Double_t        Kplus0_TRUEP_X;
    Double_t        Kplus0_TRUEP_Y;
    Double_t        Kplus0_TRUEP_Z;
    Double_t        Kplus0_TRUEPT;
    Double_t        Kplus0_TRUEORIGINVERTEX_X;
    Double_t        Kplus0_TRUEORIGINVERTEX_Y;
    Double_t        Kplus0_TRUEORIGINVERTEX_Z;
    Double_t        Kplus0_TRUEENDVERTEX_X;
    Double_t        Kplus0_TRUEENDVERTEX_Y;
    Double_t        Kplus0_TRUEENDVERTEX_Z;
    Bool_t          Kplus0_TRUEISSTABLE;
    Double_t        Kplus0_TRUETAU;
    Bool_t          Kplus0_OSCIL;
    Double_t        Kminus_nPhotos;
    Double_t        Kminus_TRUEID;
    Int_t           Kminus_MC_MOTHER_ID;
    Int_t           Kminus_MC_MOTHER_KEY;
    Int_t           Kminus_MC_GD_MOTHER_ID;
    Int_t           Kminus_MC_GD_MOTHER_KEY;
    Int_t           Kminus_MC_GD_GD_MOTHER_ID;
    Int_t           Kminus_MC_GD_GD_MOTHER_KEY;
    Double_t        Kminus_TRUEP_E;
    Double_t        Kminus_TRUEP_X;
    Double_t        Kminus_TRUEP_Y;
    Double_t        Kminus_TRUEP_Z;
    Double_t        Kminus_TRUEPT;
    Double_t        Kminus_TRUEORIGINVERTEX_X;
    Double_t        Kminus_TRUEORIGINVERTEX_Y;
    Double_t        Kminus_TRUEORIGINVERTEX_Z;
    Double_t        Kminus_TRUEENDVERTEX_X;
    Double_t        Kminus_TRUEENDVERTEX_Y;
    Double_t        Kminus_TRUEENDVERTEX_Z;
    Bool_t          Kminus_TRUEISSTABLE;
    Double_t        Kminus_TRUETAU;
    Bool_t          Kminus_OSCIL;
    Double_t        piminus0_nPhotos;
    Double_t        piminus0_TRUEID;
    Int_t           piminus0_MC_MOTHER_ID;
    Int_t           piminus0_MC_MOTHER_KEY;
    Int_t           piminus0_MC_GD_MOTHER_ID;
    Int_t           piminus0_MC_GD_MOTHER_KEY;
    Int_t           piminus0_MC_GD_GD_MOTHER_ID;
    Int_t           piminus0_MC_GD_GD_MOTHER_KEY;
    Double_t        piminus0_TRUEP_E;
    Double_t        piminus0_TRUEP_X;
    Double_t        piminus0_TRUEP_Y;
    Double_t        piminus0_TRUEP_Z;
    Double_t        piminus0_TRUEPT;
    Double_t        piminus0_TRUEORIGINVERTEX_X;
    Double_t        piminus0_TRUEORIGINVERTEX_Y;
    Double_t        piminus0_TRUEORIGINVERTEX_Z;
    Double_t        piminus0_TRUEENDVERTEX_X;
    Double_t        piminus0_TRUEENDVERTEX_Y;
    Double_t        piminus0_TRUEENDVERTEX_Z;
    Bool_t          piminus0_TRUEISSTABLE;
    Double_t        piminus0_TRUETAU;
    Bool_t          piminus0_OSCIL;
    UInt_t          nCandidate;
    ULong64_t       totCandidates;
    ULong64_t       EventInSequence;
    
    // List of branches
    TBranch        *b_B_s0_nPhotos;   //!
    TBranch        *b_B_s0_TRUEID;   //!
    TBranch        *b_B_s0_MC_MOTHER_ID;   //!
    TBranch        *b_B_s0_MC_MOTHER_KEY;   //!
    TBranch        *b_B_s0_MC_GD_MOTHER_ID;   //!
    TBranch        *b_B_s0_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_B_s0_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_B_s0_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_B_s0_TRUEP_E;   //!
    TBranch        *b_B_s0_TRUEP_X;   //!
    TBranch        *b_B_s0_TRUEP_Y;   //!
    TBranch        *b_B_s0_TRUEP_Z;   //!
    TBranch        *b_B_s0_TRUEPT;   //!
    TBranch        *b_B_s0_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_B_s0_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_B_s0_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_B_s0_TRUEENDVERTEX_X;   //!
    TBranch        *b_B_s0_TRUEENDVERTEX_Y;   //!
    TBranch        *b_B_s0_TRUEENDVERTEX_Z;   //!
    TBranch        *b_B_s0_TRUEISSTABLE;   //!
    TBranch        *b_B_s0_TRUETAU;   //!
    TBranch        *b_B_s0_OSCIL;   //!
    TBranch        *b_Kplus_nPhotos;   //!
    TBranch        *b_Kplus_TRUEID;   //!
    TBranch        *b_Kplus_MC_MOTHER_ID;   //!
    TBranch        *b_Kplus_MC_MOTHER_KEY;   //!
    TBranch        *b_Kplus_MC_GD_MOTHER_ID;   //!
    TBranch        *b_Kplus_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_Kplus_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_Kplus_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_Kplus_TRUEP_E;   //!
    TBranch        *b_Kplus_TRUEP_X;   //!
    TBranch        *b_Kplus_TRUEP_Y;   //!
    TBranch        *b_Kplus_TRUEP_Z;   //!
    TBranch        *b_Kplus_TRUEPT;   //!
    TBranch        *b_Kplus_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_Kplus_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_Kplus_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_Kplus_TRUEENDVERTEX_X;   //!
    TBranch        *b_Kplus_TRUEENDVERTEX_Y;   //!
    TBranch        *b_Kplus_TRUEENDVERTEX_Z;   //!
    TBranch        *b_Kplus_TRUEISSTABLE;   //!
    TBranch        *b_Kplus_TRUETAU;   //!
    TBranch        *b_Kplus_OSCIL;   //!
    TBranch        *b_piplus_nPhotos;   //!
    TBranch        *b_piplus_TRUEID;   //!
    TBranch        *b_piplus_MC_MOTHER_ID;   //!
    TBranch        *b_piplus_MC_MOTHER_KEY;   //!
    TBranch        *b_piplus_MC_GD_MOTHER_ID;   //!
    TBranch        *b_piplus_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_piplus_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_piplus_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_piplus_TRUEP_E;   //!
    TBranch        *b_piplus_TRUEP_X;   //!
    TBranch        *b_piplus_TRUEP_Y;   //!
    TBranch        *b_piplus_TRUEP_Z;   //!
    TBranch        *b_piplus_TRUEPT;   //!
    TBranch        *b_piplus_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_piplus_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_piplus_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_piplus_TRUEENDVERTEX_X;   //!
    TBranch        *b_piplus_TRUEENDVERTEX_Y;   //!
    TBranch        *b_piplus_TRUEENDVERTEX_Z;   //!
    TBranch        *b_piplus_TRUEISSTABLE;   //!
    TBranch        *b_piplus_TRUETAU;   //!
    TBranch        *b_piplus_OSCIL;   //!
    TBranch        *b_piminus_nPhotos;   //!
    TBranch        *b_piminus_TRUEID;   //!
    TBranch        *b_piminus_MC_MOTHER_ID;   //!
    TBranch        *b_piminus_MC_MOTHER_KEY;   //!
    TBranch        *b_piminus_MC_GD_MOTHER_ID;   //!
    TBranch        *b_piminus_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_piminus_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_piminus_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_piminus_TRUEP_E;   //!
    TBranch        *b_piminus_TRUEP_X;   //!
    TBranch        *b_piminus_TRUEP_Y;   //!
    TBranch        *b_piminus_TRUEP_Z;   //!
    TBranch        *b_piminus_TRUEPT;   //!
    TBranch        *b_piminus_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_piminus_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_piminus_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_piminus_TRUEENDVERTEX_X;   //!
    TBranch        *b_piminus_TRUEENDVERTEX_Y;   //!
    TBranch        *b_piminus_TRUEENDVERTEX_Z;   //!
    TBranch        *b_piminus_TRUEISSTABLE;   //!
    TBranch        *b_piminus_TRUETAU;   //!
    TBranch        *b_piminus_OSCIL;   //!
    TBranch        *b_D_sminus_nPhotos;   //!
    TBranch        *b_D_sminus_TRUEID;   //!
    TBranch        *b_D_sminus_MC_MOTHER_ID;   //!
    TBranch        *b_D_sminus_MC_MOTHER_KEY;   //!
    TBranch        *b_D_sminus_MC_GD_MOTHER_ID;   //!
    TBranch        *b_D_sminus_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_D_sminus_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_D_sminus_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_D_sminus_TRUEP_E;   //!
    TBranch        *b_D_sminus_TRUEP_X;   //!
    TBranch        *b_D_sminus_TRUEP_Y;   //!
    TBranch        *b_D_sminus_TRUEP_Z;   //!
    TBranch        *b_D_sminus_TRUEPT;   //!
    TBranch        *b_D_sminus_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_D_sminus_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_D_sminus_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_D_sminus_TRUEENDVERTEX_X;   //!
    TBranch        *b_D_sminus_TRUEENDVERTEX_Y;   //!
    TBranch        *b_D_sminus_TRUEENDVERTEX_Z;   //!
    TBranch        *b_D_sminus_TRUEISSTABLE;   //!
    TBranch        *b_D_sminus_TRUETAU;   //!
    TBranch        *b_D_sminus_OSCIL;   //!
    TBranch        *b_Kplus0_nPhotos;   //!
    TBranch        *b_Kplus0_TRUEID;   //!
    TBranch        *b_Kplus0_MC_MOTHER_ID;   //!
    TBranch        *b_Kplus0_MC_MOTHER_KEY;   //!
    TBranch        *b_Kplus0_MC_GD_MOTHER_ID;   //!
    TBranch        *b_Kplus0_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_Kplus0_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_Kplus0_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_Kplus0_TRUEP_E;   //!
    TBranch        *b_Kplus0_TRUEP_X;   //!
    TBranch        *b_Kplus0_TRUEP_Y;   //!
    TBranch        *b_Kplus0_TRUEP_Z;   //!
    TBranch        *b_Kplus0_TRUEPT;   //!
    TBranch        *b_Kplus0_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_Kplus0_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_Kplus0_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_Kplus0_TRUEENDVERTEX_X;   //!
    TBranch        *b_Kplus0_TRUEENDVERTEX_Y;   //!
    TBranch        *b_Kplus0_TRUEENDVERTEX_Z;   //!
    TBranch        *b_Kplus0_TRUEISSTABLE;   //!
    TBranch        *b_Kplus0_TRUETAU;   //!
    TBranch        *b_Kplus0_OSCIL;   //!
    TBranch        *b_Kminus_nPhotos;   //!
    TBranch        *b_Kminus_TRUEID;   //!
    TBranch        *b_Kminus_MC_MOTHER_ID;   //!
    TBranch        *b_Kminus_MC_MOTHER_KEY;   //!
    TBranch        *b_Kminus_MC_GD_MOTHER_ID;   //!
    TBranch        *b_Kminus_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_Kminus_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_Kminus_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_Kminus_TRUEP_E;   //!
    TBranch        *b_Kminus_TRUEP_X;   //!
    TBranch        *b_Kminus_TRUEP_Y;   //!
    TBranch        *b_Kminus_TRUEP_Z;   //!
    TBranch        *b_Kminus_TRUEPT;   //!
    TBranch        *b_Kminus_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_Kminus_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_Kminus_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_Kminus_TRUEENDVERTEX_X;   //!
    TBranch        *b_Kminus_TRUEENDVERTEX_Y;   //!
    TBranch        *b_Kminus_TRUEENDVERTEX_Z;   //!
    TBranch        *b_Kminus_TRUEISSTABLE;   //!
    TBranch        *b_Kminus_TRUETAU;   //!
    TBranch        *b_Kminus_OSCIL;   //!
    TBranch        *b_piminus0_nPhotos;   //!
    TBranch        *b_piminus0_TRUEID;   //!
    TBranch        *b_piminus0_MC_MOTHER_ID;   //!
    TBranch        *b_piminus0_MC_MOTHER_KEY;   //!
    TBranch        *b_piminus0_MC_GD_MOTHER_ID;   //!
    TBranch        *b_piminus0_MC_GD_MOTHER_KEY;   //!
    TBranch        *b_piminus0_MC_GD_GD_MOTHER_ID;   //!
    TBranch        *b_piminus0_MC_GD_GD_MOTHER_KEY;   //!
    TBranch        *b_piminus0_TRUEP_E;   //!
    TBranch        *b_piminus0_TRUEP_X;   //!
    TBranch        *b_piminus0_TRUEP_Y;   //!
    TBranch        *b_piminus0_TRUEP_Z;   //!
    TBranch        *b_piminus0_TRUEPT;   //!
    TBranch        *b_piminus0_TRUEORIGINVERTEX_X;   //!
    TBranch        *b_piminus0_TRUEORIGINVERTEX_Y;   //!
    TBranch        *b_piminus0_TRUEORIGINVERTEX_Z;   //!
    TBranch        *b_piminus0_TRUEENDVERTEX_X;   //!
    TBranch        *b_piminus0_TRUEENDVERTEX_Y;   //!
    TBranch        *b_piminus0_TRUEENDVERTEX_Z;   //!
    TBranch        *b_piminus0_TRUEISSTABLE;   //!
    TBranch        *b_piminus0_TRUETAU;   //!
    TBranch        *b_piminus0_OSCIL;   //!
    TBranch        *b_nCandidate;   //!
    TBranch        *b_totCandidates;   //!
    TBranch        *b_EventInSequence;   //!

   MCDecayTreeTuple(TTree *tree=0);
   virtual ~MCDecayTreeTuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MCDecayTreeTuple_cxx
MCDecayTreeTuple::MCDecayTreeTuple(TTree *tree) : fChain(0) 
{   
   TChain* chain = new TChain("MCDecayTreeTuple/MCDecayTree");
   
   chain->Add("DVntuple4.root");
     
   Init((TTree*)chain);
}

MCDecayTreeTuple::~MCDecayTreeTuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MCDecayTreeTuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MCDecayTreeTuple::LoadTree(Long64_t entry)
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

void MCDecayTreeTuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

    fChain->SetBranchAddress("B_s0_nPhotos", &B_s0_nPhotos, &b_B_s0_nPhotos);
    fChain->SetBranchAddress("B_s0_TRUEID", &B_s0_TRUEID, &b_B_s0_TRUEID);
    fChain->SetBranchAddress("B_s0_MC_MOTHER_ID", &B_s0_MC_MOTHER_ID, &b_B_s0_MC_MOTHER_ID);
    fChain->SetBranchAddress("B_s0_MC_MOTHER_KEY", &B_s0_MC_MOTHER_KEY, &b_B_s0_MC_MOTHER_KEY);
    fChain->SetBranchAddress("B_s0_MC_GD_MOTHER_ID", &B_s0_MC_GD_MOTHER_ID, &b_B_s0_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("B_s0_MC_GD_MOTHER_KEY", &B_s0_MC_GD_MOTHER_KEY, &b_B_s0_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("B_s0_MC_GD_GD_MOTHER_ID", &B_s0_MC_GD_GD_MOTHER_ID, &b_B_s0_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("B_s0_MC_GD_GD_MOTHER_KEY", &B_s0_MC_GD_GD_MOTHER_KEY, &b_B_s0_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("B_s0_TRUEP_E", &B_s0_TRUEP_E, &b_B_s0_TRUEP_E);
    fChain->SetBranchAddress("B_s0_TRUEP_X", &B_s0_TRUEP_X, &b_B_s0_TRUEP_X);
    fChain->SetBranchAddress("B_s0_TRUEP_Y", &B_s0_TRUEP_Y, &b_B_s0_TRUEP_Y);
    fChain->SetBranchAddress("B_s0_TRUEP_Z", &B_s0_TRUEP_Z, &b_B_s0_TRUEP_Z);
    fChain->SetBranchAddress("B_s0_TRUEPT", &B_s0_TRUEPT, &b_B_s0_TRUEPT);
    fChain->SetBranchAddress("B_s0_TRUEORIGINVERTEX_X", &B_s0_TRUEORIGINVERTEX_X, &b_B_s0_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("B_s0_TRUEORIGINVERTEX_Y", &B_s0_TRUEORIGINVERTEX_Y, &b_B_s0_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("B_s0_TRUEORIGINVERTEX_Z", &B_s0_TRUEORIGINVERTEX_Z, &b_B_s0_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("B_s0_TRUEENDVERTEX_X", &B_s0_TRUEENDVERTEX_X, &b_B_s0_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("B_s0_TRUEENDVERTEX_Y", &B_s0_TRUEENDVERTEX_Y, &b_B_s0_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("B_s0_TRUEENDVERTEX_Z", &B_s0_TRUEENDVERTEX_Z, &b_B_s0_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("B_s0_TRUEISSTABLE", &B_s0_TRUEISSTABLE, &b_B_s0_TRUEISSTABLE);
    fChain->SetBranchAddress("B_s0_TRUETAU", &B_s0_TRUETAU, &b_B_s0_TRUETAU);
    fChain->SetBranchAddress("B_s0_OSCIL", &B_s0_OSCIL, &b_B_s0_OSCIL);
    fChain->SetBranchAddress("Kplus_nPhotos", &Kplus_nPhotos, &b_Kplus_nPhotos);
    fChain->SetBranchAddress("Kplus_TRUEID", &Kplus_TRUEID, &b_Kplus_TRUEID);
    fChain->SetBranchAddress("Kplus_MC_MOTHER_ID", &Kplus_MC_MOTHER_ID, &b_Kplus_MC_MOTHER_ID);
    fChain->SetBranchAddress("Kplus_MC_MOTHER_KEY", &Kplus_MC_MOTHER_KEY, &b_Kplus_MC_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus_MC_GD_MOTHER_ID", &Kplus_MC_GD_MOTHER_ID, &b_Kplus_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kplus_MC_GD_MOTHER_KEY", &Kplus_MC_GD_MOTHER_KEY, &b_Kplus_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus_MC_GD_GD_MOTHER_ID", &Kplus_MC_GD_GD_MOTHER_ID, &b_Kplus_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kplus_MC_GD_GD_MOTHER_KEY", &Kplus_MC_GD_GD_MOTHER_KEY, &b_Kplus_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus_TRUEP_E", &Kplus_TRUEP_E, &b_Kplus_TRUEP_E);
    fChain->SetBranchAddress("Kplus_TRUEP_X", &Kplus_TRUEP_X, &b_Kplus_TRUEP_X);
    fChain->SetBranchAddress("Kplus_TRUEP_Y", &Kplus_TRUEP_Y, &b_Kplus_TRUEP_Y);
    fChain->SetBranchAddress("Kplus_TRUEP_Z", &Kplus_TRUEP_Z, &b_Kplus_TRUEP_Z);
    fChain->SetBranchAddress("Kplus_TRUEPT", &Kplus_TRUEPT, &b_Kplus_TRUEPT);
    fChain->SetBranchAddress("Kplus_TRUEORIGINVERTEX_X", &Kplus_TRUEORIGINVERTEX_X, &b_Kplus_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("Kplus_TRUEORIGINVERTEX_Y", &Kplus_TRUEORIGINVERTEX_Y, &b_Kplus_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("Kplus_TRUEORIGINVERTEX_Z", &Kplus_TRUEORIGINVERTEX_Z, &b_Kplus_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("Kplus_TRUEENDVERTEX_X", &Kplus_TRUEENDVERTEX_X, &b_Kplus_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("Kplus_TRUEENDVERTEX_Y", &Kplus_TRUEENDVERTEX_Y, &b_Kplus_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("Kplus_TRUEENDVERTEX_Z", &Kplus_TRUEENDVERTEX_Z, &b_Kplus_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("Kplus_TRUEISSTABLE", &Kplus_TRUEISSTABLE, &b_Kplus_TRUEISSTABLE);
    fChain->SetBranchAddress("Kplus_TRUETAU", &Kplus_TRUETAU, &b_Kplus_TRUETAU);
    fChain->SetBranchAddress("Kplus_OSCIL", &Kplus_OSCIL, &b_Kplus_OSCIL);
    fChain->SetBranchAddress("piplus_nPhotos", &piplus_nPhotos, &b_piplus_nPhotos);
    fChain->SetBranchAddress("piplus_TRUEID", &piplus_TRUEID, &b_piplus_TRUEID);
    fChain->SetBranchAddress("piplus_MC_MOTHER_ID", &piplus_MC_MOTHER_ID, &b_piplus_MC_MOTHER_ID);
    fChain->SetBranchAddress("piplus_MC_MOTHER_KEY", &piplus_MC_MOTHER_KEY, &b_piplus_MC_MOTHER_KEY);
    fChain->SetBranchAddress("piplus_MC_GD_MOTHER_ID", &piplus_MC_GD_MOTHER_ID, &b_piplus_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("piplus_MC_GD_MOTHER_KEY", &piplus_MC_GD_MOTHER_KEY, &b_piplus_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piplus_MC_GD_GD_MOTHER_ID", &piplus_MC_GD_GD_MOTHER_ID, &b_piplus_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("piplus_MC_GD_GD_MOTHER_KEY", &piplus_MC_GD_GD_MOTHER_KEY, &b_piplus_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piplus_TRUEP_E", &piplus_TRUEP_E, &b_piplus_TRUEP_E);
    fChain->SetBranchAddress("piplus_TRUEP_X", &piplus_TRUEP_X, &b_piplus_TRUEP_X);
    fChain->SetBranchAddress("piplus_TRUEP_Y", &piplus_TRUEP_Y, &b_piplus_TRUEP_Y);
    fChain->SetBranchAddress("piplus_TRUEP_Z", &piplus_TRUEP_Z, &b_piplus_TRUEP_Z);
    fChain->SetBranchAddress("piplus_TRUEPT", &piplus_TRUEPT, &b_piplus_TRUEPT);
    fChain->SetBranchAddress("piplus_TRUEORIGINVERTEX_X", &piplus_TRUEORIGINVERTEX_X, &b_piplus_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("piplus_TRUEORIGINVERTEX_Y", &piplus_TRUEORIGINVERTEX_Y, &b_piplus_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("piplus_TRUEORIGINVERTEX_Z", &piplus_TRUEORIGINVERTEX_Z, &b_piplus_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("piplus_TRUEENDVERTEX_X", &piplus_TRUEENDVERTEX_X, &b_piplus_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("piplus_TRUEENDVERTEX_Y", &piplus_TRUEENDVERTEX_Y, &b_piplus_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("piplus_TRUEENDVERTEX_Z", &piplus_TRUEENDVERTEX_Z, &b_piplus_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("piplus_TRUEISSTABLE", &piplus_TRUEISSTABLE, &b_piplus_TRUEISSTABLE);
    fChain->SetBranchAddress("piplus_TRUETAU", &piplus_TRUETAU, &b_piplus_TRUETAU);
    fChain->SetBranchAddress("piplus_OSCIL", &piplus_OSCIL, &b_piplus_OSCIL);
    fChain->SetBranchAddress("piminus_nPhotos", &piminus_nPhotos, &b_piminus_nPhotos);
    fChain->SetBranchAddress("piminus_TRUEID", &piminus_TRUEID, &b_piminus_TRUEID);
    fChain->SetBranchAddress("piminus_MC_MOTHER_ID", &piminus_MC_MOTHER_ID, &b_piminus_MC_MOTHER_ID);
    fChain->SetBranchAddress("piminus_MC_MOTHER_KEY", &piminus_MC_MOTHER_KEY, &b_piminus_MC_MOTHER_KEY);
    fChain->SetBranchAddress("piminus_MC_GD_MOTHER_ID", &piminus_MC_GD_MOTHER_ID, &b_piminus_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("piminus_MC_GD_MOTHER_KEY", &piminus_MC_GD_MOTHER_KEY, &b_piminus_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piminus_MC_GD_GD_MOTHER_ID", &piminus_MC_GD_GD_MOTHER_ID, &b_piminus_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("piminus_MC_GD_GD_MOTHER_KEY", &piminus_MC_GD_GD_MOTHER_KEY, &b_piminus_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piminus_TRUEP_E", &piminus_TRUEP_E, &b_piminus_TRUEP_E);
    fChain->SetBranchAddress("piminus_TRUEP_X", &piminus_TRUEP_X, &b_piminus_TRUEP_X);
    fChain->SetBranchAddress("piminus_TRUEP_Y", &piminus_TRUEP_Y, &b_piminus_TRUEP_Y);
    fChain->SetBranchAddress("piminus_TRUEP_Z", &piminus_TRUEP_Z, &b_piminus_TRUEP_Z);
    fChain->SetBranchAddress("piminus_TRUEPT", &piminus_TRUEPT, &b_piminus_TRUEPT);
    fChain->SetBranchAddress("piminus_TRUEORIGINVERTEX_X", &piminus_TRUEORIGINVERTEX_X, &b_piminus_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("piminus_TRUEORIGINVERTEX_Y", &piminus_TRUEORIGINVERTEX_Y, &b_piminus_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("piminus_TRUEORIGINVERTEX_Z", &piminus_TRUEORIGINVERTEX_Z, &b_piminus_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("piminus_TRUEENDVERTEX_X", &piminus_TRUEENDVERTEX_X, &b_piminus_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("piminus_TRUEENDVERTEX_Y", &piminus_TRUEENDVERTEX_Y, &b_piminus_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("piminus_TRUEENDVERTEX_Z", &piminus_TRUEENDVERTEX_Z, &b_piminus_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("piminus_TRUEISSTABLE", &piminus_TRUEISSTABLE, &b_piminus_TRUEISSTABLE);
    fChain->SetBranchAddress("piminus_TRUETAU", &piminus_TRUETAU, &b_piminus_TRUETAU);
    fChain->SetBranchAddress("piminus_OSCIL", &piminus_OSCIL, &b_piminus_OSCIL);
    fChain->SetBranchAddress("D_sminus_nPhotos", &D_sminus_nPhotos, &b_D_sminus_nPhotos);
    fChain->SetBranchAddress("D_sminus_TRUEID", &D_sminus_TRUEID, &b_D_sminus_TRUEID);
    fChain->SetBranchAddress("D_sminus_MC_MOTHER_ID", &D_sminus_MC_MOTHER_ID, &b_D_sminus_MC_MOTHER_ID);
    fChain->SetBranchAddress("D_sminus_MC_MOTHER_KEY", &D_sminus_MC_MOTHER_KEY, &b_D_sminus_MC_MOTHER_KEY);
    fChain->SetBranchAddress("D_sminus_MC_GD_MOTHER_ID", &D_sminus_MC_GD_MOTHER_ID, &b_D_sminus_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("D_sminus_MC_GD_MOTHER_KEY", &D_sminus_MC_GD_MOTHER_KEY, &b_D_sminus_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("D_sminus_MC_GD_GD_MOTHER_ID", &D_sminus_MC_GD_GD_MOTHER_ID, &b_D_sminus_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("D_sminus_MC_GD_GD_MOTHER_KEY", &D_sminus_MC_GD_GD_MOTHER_KEY, &b_D_sminus_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("D_sminus_TRUEP_E", &D_sminus_TRUEP_E, &b_D_sminus_TRUEP_E);
    fChain->SetBranchAddress("D_sminus_TRUEP_X", &D_sminus_TRUEP_X, &b_D_sminus_TRUEP_X);
    fChain->SetBranchAddress("D_sminus_TRUEP_Y", &D_sminus_TRUEP_Y, &b_D_sminus_TRUEP_Y);
    fChain->SetBranchAddress("D_sminus_TRUEP_Z", &D_sminus_TRUEP_Z, &b_D_sminus_TRUEP_Z);
    fChain->SetBranchAddress("D_sminus_TRUEPT", &D_sminus_TRUEPT, &b_D_sminus_TRUEPT);
    fChain->SetBranchAddress("D_sminus_TRUEORIGINVERTEX_X", &D_sminus_TRUEORIGINVERTEX_X, &b_D_sminus_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("D_sminus_TRUEORIGINVERTEX_Y", &D_sminus_TRUEORIGINVERTEX_Y, &b_D_sminus_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("D_sminus_TRUEORIGINVERTEX_Z", &D_sminus_TRUEORIGINVERTEX_Z, &b_D_sminus_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("D_sminus_TRUEENDVERTEX_X", &D_sminus_TRUEENDVERTEX_X, &b_D_sminus_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("D_sminus_TRUEENDVERTEX_Y", &D_sminus_TRUEENDVERTEX_Y, &b_D_sminus_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("D_sminus_TRUEENDVERTEX_Z", &D_sminus_TRUEENDVERTEX_Z, &b_D_sminus_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("D_sminus_TRUEISSTABLE", &D_sminus_TRUEISSTABLE, &b_D_sminus_TRUEISSTABLE);
    fChain->SetBranchAddress("D_sminus_TRUETAU", &D_sminus_TRUETAU, &b_D_sminus_TRUETAU);
    fChain->SetBranchAddress("D_sminus_OSCIL", &D_sminus_OSCIL, &b_D_sminus_OSCIL);
    fChain->SetBranchAddress("Kplus0_nPhotos", &Kplus0_nPhotos, &b_Kplus0_nPhotos);
    fChain->SetBranchAddress("Kplus0_TRUEID", &Kplus0_TRUEID, &b_Kplus0_TRUEID);
    fChain->SetBranchAddress("Kplus0_MC_MOTHER_ID", &Kplus0_MC_MOTHER_ID, &b_Kplus0_MC_MOTHER_ID);
    fChain->SetBranchAddress("Kplus0_MC_MOTHER_KEY", &Kplus0_MC_MOTHER_KEY, &b_Kplus0_MC_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus0_MC_GD_MOTHER_ID", &Kplus0_MC_GD_MOTHER_ID, &b_Kplus0_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kplus0_MC_GD_MOTHER_KEY", &Kplus0_MC_GD_MOTHER_KEY, &b_Kplus0_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus0_MC_GD_GD_MOTHER_ID", &Kplus0_MC_GD_GD_MOTHER_ID, &b_Kplus0_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kplus0_MC_GD_GD_MOTHER_KEY", &Kplus0_MC_GD_GD_MOTHER_KEY, &b_Kplus0_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kplus0_TRUEP_E", &Kplus0_TRUEP_E, &b_Kplus0_TRUEP_E);
    fChain->SetBranchAddress("Kplus0_TRUEP_X", &Kplus0_TRUEP_X, &b_Kplus0_TRUEP_X);
    fChain->SetBranchAddress("Kplus0_TRUEP_Y", &Kplus0_TRUEP_Y, &b_Kplus0_TRUEP_Y);
    fChain->SetBranchAddress("Kplus0_TRUEP_Z", &Kplus0_TRUEP_Z, &b_Kplus0_TRUEP_Z);
    fChain->SetBranchAddress("Kplus0_TRUEPT", &Kplus0_TRUEPT, &b_Kplus0_TRUEPT);
    fChain->SetBranchAddress("Kplus0_TRUEORIGINVERTEX_X", &Kplus0_TRUEORIGINVERTEX_X, &b_Kplus0_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("Kplus0_TRUEORIGINVERTEX_Y", &Kplus0_TRUEORIGINVERTEX_Y, &b_Kplus0_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("Kplus0_TRUEORIGINVERTEX_Z", &Kplus0_TRUEORIGINVERTEX_Z, &b_Kplus0_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("Kplus0_TRUEENDVERTEX_X", &Kplus0_TRUEENDVERTEX_X, &b_Kplus0_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("Kplus0_TRUEENDVERTEX_Y", &Kplus0_TRUEENDVERTEX_Y, &b_Kplus0_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("Kplus0_TRUEENDVERTEX_Z", &Kplus0_TRUEENDVERTEX_Z, &b_Kplus0_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("Kplus0_TRUEISSTABLE", &Kplus0_TRUEISSTABLE, &b_Kplus0_TRUEISSTABLE);
    fChain->SetBranchAddress("Kplus0_TRUETAU", &Kplus0_TRUETAU, &b_Kplus0_TRUETAU);
    fChain->SetBranchAddress("Kplus0_OSCIL", &Kplus0_OSCIL, &b_Kplus0_OSCIL);
    fChain->SetBranchAddress("Kminus_nPhotos", &Kminus_nPhotos, &b_Kminus_nPhotos);
    fChain->SetBranchAddress("Kminus_TRUEID", &Kminus_TRUEID, &b_Kminus_TRUEID);
    fChain->SetBranchAddress("Kminus_MC_MOTHER_ID", &Kminus_MC_MOTHER_ID, &b_Kminus_MC_MOTHER_ID);
    fChain->SetBranchAddress("Kminus_MC_MOTHER_KEY", &Kminus_MC_MOTHER_KEY, &b_Kminus_MC_MOTHER_KEY);
    fChain->SetBranchAddress("Kminus_MC_GD_MOTHER_ID", &Kminus_MC_GD_MOTHER_ID, &b_Kminus_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kminus_MC_GD_MOTHER_KEY", &Kminus_MC_GD_MOTHER_KEY, &b_Kminus_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kminus_MC_GD_GD_MOTHER_ID", &Kminus_MC_GD_GD_MOTHER_ID, &b_Kminus_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("Kminus_MC_GD_GD_MOTHER_KEY", &Kminus_MC_GD_GD_MOTHER_KEY, &b_Kminus_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("Kminus_TRUEP_E", &Kminus_TRUEP_E, &b_Kminus_TRUEP_E);
    fChain->SetBranchAddress("Kminus_TRUEP_X", &Kminus_TRUEP_X, &b_Kminus_TRUEP_X);
    fChain->SetBranchAddress("Kminus_TRUEP_Y", &Kminus_TRUEP_Y, &b_Kminus_TRUEP_Y);
    fChain->SetBranchAddress("Kminus_TRUEP_Z", &Kminus_TRUEP_Z, &b_Kminus_TRUEP_Z);
    fChain->SetBranchAddress("Kminus_TRUEPT", &Kminus_TRUEPT, &b_Kminus_TRUEPT);
    fChain->SetBranchAddress("Kminus_TRUEORIGINVERTEX_X", &Kminus_TRUEORIGINVERTEX_X, &b_Kminus_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("Kminus_TRUEORIGINVERTEX_Y", &Kminus_TRUEORIGINVERTEX_Y, &b_Kminus_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("Kminus_TRUEORIGINVERTEX_Z", &Kminus_TRUEORIGINVERTEX_Z, &b_Kminus_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("Kminus_TRUEENDVERTEX_X", &Kminus_TRUEENDVERTEX_X, &b_Kminus_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("Kminus_TRUEENDVERTEX_Y", &Kminus_TRUEENDVERTEX_Y, &b_Kminus_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("Kminus_TRUEENDVERTEX_Z", &Kminus_TRUEENDVERTEX_Z, &b_Kminus_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("Kminus_TRUEISSTABLE", &Kminus_TRUEISSTABLE, &b_Kminus_TRUEISSTABLE);
    fChain->SetBranchAddress("Kminus_TRUETAU", &Kminus_TRUETAU, &b_Kminus_TRUETAU);
    fChain->SetBranchAddress("Kminus_OSCIL", &Kminus_OSCIL, &b_Kminus_OSCIL);
    fChain->SetBranchAddress("piminus0_nPhotos", &piminus0_nPhotos, &b_piminus0_nPhotos);
    fChain->SetBranchAddress("piminus0_TRUEID", &piminus0_TRUEID, &b_piminus0_TRUEID);
    fChain->SetBranchAddress("piminus0_MC_MOTHER_ID", &piminus0_MC_MOTHER_ID, &b_piminus0_MC_MOTHER_ID);
    fChain->SetBranchAddress("piminus0_MC_MOTHER_KEY", &piminus0_MC_MOTHER_KEY, &b_piminus0_MC_MOTHER_KEY);
    fChain->SetBranchAddress("piminus0_MC_GD_MOTHER_ID", &piminus0_MC_GD_MOTHER_ID, &b_piminus0_MC_GD_MOTHER_ID);
    fChain->SetBranchAddress("piminus0_MC_GD_MOTHER_KEY", &piminus0_MC_GD_MOTHER_KEY, &b_piminus0_MC_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piminus0_MC_GD_GD_MOTHER_ID", &piminus0_MC_GD_GD_MOTHER_ID, &b_piminus0_MC_GD_GD_MOTHER_ID);
    fChain->SetBranchAddress("piminus0_MC_GD_GD_MOTHER_KEY", &piminus0_MC_GD_GD_MOTHER_KEY, &b_piminus0_MC_GD_GD_MOTHER_KEY);
    fChain->SetBranchAddress("piminus0_TRUEP_E", &piminus0_TRUEP_E, &b_piminus0_TRUEP_E);
    fChain->SetBranchAddress("piminus0_TRUEP_X", &piminus0_TRUEP_X, &b_piminus0_TRUEP_X);
    fChain->SetBranchAddress("piminus0_TRUEP_Y", &piminus0_TRUEP_Y, &b_piminus0_TRUEP_Y);
    fChain->SetBranchAddress("piminus0_TRUEP_Z", &piminus0_TRUEP_Z, &b_piminus0_TRUEP_Z);
    fChain->SetBranchAddress("piminus0_TRUEPT", &piminus0_TRUEPT, &b_piminus0_TRUEPT);
    fChain->SetBranchAddress("piminus0_TRUEORIGINVERTEX_X", &piminus0_TRUEORIGINVERTEX_X, &b_piminus0_TRUEORIGINVERTEX_X);
    fChain->SetBranchAddress("piminus0_TRUEORIGINVERTEX_Y", &piminus0_TRUEORIGINVERTEX_Y, &b_piminus0_TRUEORIGINVERTEX_Y);
    fChain->SetBranchAddress("piminus0_TRUEORIGINVERTEX_Z", &piminus0_TRUEORIGINVERTEX_Z, &b_piminus0_TRUEORIGINVERTEX_Z);
    fChain->SetBranchAddress("piminus0_TRUEENDVERTEX_X", &piminus0_TRUEENDVERTEX_X, &b_piminus0_TRUEENDVERTEX_X);
    fChain->SetBranchAddress("piminus0_TRUEENDVERTEX_Y", &piminus0_TRUEENDVERTEX_Y, &b_piminus0_TRUEENDVERTEX_Y);
    fChain->SetBranchAddress("piminus0_TRUEENDVERTEX_Z", &piminus0_TRUEENDVERTEX_Z, &b_piminus0_TRUEENDVERTEX_Z);
    fChain->SetBranchAddress("piminus0_TRUEISSTABLE", &piminus0_TRUEISSTABLE, &b_piminus0_TRUEISSTABLE);
    fChain->SetBranchAddress("piminus0_TRUETAU", &piminus0_TRUETAU, &b_piminus0_TRUETAU);
    fChain->SetBranchAddress("piminus0_OSCIL", &piminus0_OSCIL, &b_piminus0_OSCIL);
    fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
    fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
    fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   Notify();
}

Bool_t MCDecayTreeTuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MCDecayTreeTuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MCDecayTreeTuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MCDecayTreeTuple_cxx
