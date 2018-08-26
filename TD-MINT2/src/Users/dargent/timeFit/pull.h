//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 24 11:58:28 2018 by ROOT version 5.34/10
// from TTree MinuitParameterSetNtp/MinuitParameterSetNtp
// found on file: signal_toy/pull_par0_1.root
//////////////////////////////////////////////////////////

#ifndef pull_h
#define pull_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

class pull {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   vector<TString> _paraNames; 
   TString _fileName;    
   vector<double*> _means;
   vector<double*> _inits;
   vector<double*> _errs;
   vector<double*> _pulls;
    
   pull(vector<TString> paraNames, TString fileName);
   virtual ~pull();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   
   TMatrixD getStatCov(TString label = "");
   TMatrixD getCov(TString label = "");
   TMatrixD getDeltaCov(TString refFileName,TString label = "");
   TMatrixD getDeltaCovChol(TString refFileName,TString label,int varPerParChol);
    
   TMatrixD 	combineCov_maxVal(vector<TMatrixD*> vec);

   TString latexName(TString s); 
    
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pull_cxx

pull::pull(vector<TString> paraNames, TString fileName):
    fChain(0), _paraNames(paraNames), _fileName(fileName)
{
    TChain* chain =  new TChain("MinuitParameterSetNtp");
    chain->Add(_fileName); 

    int N = chain->GetEntries();
    TChain* chain2 =  new TChain("MinuitParameterSetNtp");
    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = fileName;
	chain2->Add(file.ReplaceAll("*",index.str())); 
    }   
    else chain2->Add(_fileName); 

    Init((TTree*)chain2);
}	

pull::~pull()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pull::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pull::LoadTree(Long64_t entry)
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

void pull::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
    
    for (int i = 0 ; i < _paraNames.size(); i++) {
        double * mean = new double[1];
        _means.push_back(mean);
        fChain->SetBranchAddress(_paraNames[i]+"_mean", mean);

        double * init = new double[1];
        _inits.push_back(init);
        fChain->SetBranchAddress(_paraNames[i]+"_init", init);

        double * err = new double[1];
        _errs.push_back(err);
        fChain->SetBranchAddress(_paraNames[i]+"_err", err);

        double * pull = new double[1];
        _pulls.push_back(pull);
        fChain->SetBranchAddress(_paraNames[i]+"_pull", pull);
    } 
   Notify();
}

Bool_t pull::Notify()
{
   return kTRUE;
}

void pull::Show(Long64_t entry)
{
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pull::Cut(Long64_t entry)
{
   return 1;
}
#endif // #ifdef pull_cxx
