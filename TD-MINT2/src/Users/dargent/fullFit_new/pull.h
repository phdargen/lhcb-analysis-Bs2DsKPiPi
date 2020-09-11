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
   TString _treeName;    
   bool	_fraction;
   bool	_removeBar;
   int _nFiles;

   vector<double*> _means;
   vector<double*> _inits;
   vector<double*> _errs;
   vector<double*> _pulls;
   vector<int> _skip;    

   pull(vector<TString> paraNames, TString fileName, TString treeName = "MinuitParameterSetNtp", bool fraction = false, bool removeBar = false, int nFiles = -1);
   virtual ~pull();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);

   vector <double> getVals();   
   vector <double> getErrs();

   vector <double> getPullMean();
   vector <double> getPullSigma();

   TMatrixD getStatCov(TString label = "");
   TMatrixD getCov(TString label = "",double pull_max = -1, bool scaleCov = true);
   TMatrixD getDeltaCov(TString refFileName,TString label = "",double pull_max = 5, bool scaleCov = true);
   TMatrixD getDeltaCovChol(TString refFileName,TString label,int varPerParChol);
   TMatrixD getAbsDiff(TMatrixD cov1,TMatrixD cov2);

   TMatrixD 	combineCov_maxVal(vector<TMatrixD*> vec);

   vector<double> sampleMean();
   vector<double> sampleSigma();

   vector<double> sampleMean(vector< vector<double> > vec_vals);
   TMatrixD     sampleVariance(vector< vector<double> > vec_vals, int removeOutlier = 0);

   TString latexName(TString s); 
   TString latexNameMod(TString s); 
    
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pull_cxx

pull::pull(vector<TString> paraNames, TString fileName, TString treeName, bool fraction, bool removeBar, int nFiles):
    fChain(0), _paraNames(paraNames), _fileName(fileName), _treeName(treeName), _fraction(fraction), _removeBar(removeBar), _nFiles(nFiles)
{

    int N = _nFiles;
    if(_nFiles < 0){
	TChain* chain =  new TChain(treeName);
    	chain->Add(_fileName); 
    	N = chain->GetEntries();
    }
    if(removeBar)fileName.ReplaceAll("_Bar","");

    TChain* chain2 =  new TChain(treeName);
    TString lastFile;

    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = fileName;
	file.ReplaceAll("*",index.str());
	if(std::ifstream(((string)file).c_str()).good()){ 
		chain2->Add(file);
		lastFile = file;
	}
	else {
		chain2->Add(lastFile);
		_skip.push_back(i);
	} 
    }   
    else chain2->Add(fileName); 

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
	if(_fraction){
		fChain->SetBranchAddress(_paraNames[i], mean);
		continue;
	}
        int b_mean = fChain->SetBranchAddress(_paraNames[i]+"_mean", mean);
	if(b_mean) *mean = 0.;

        double * init = new double[1];
        _inits.push_back(init);
        int b_init = fChain->SetBranchAddress(_paraNames[i]+"_init", init);
	if(b_init) *init = 0.;

        double * err = new double[1];
        _errs.push_back(err);
        int b_err = fChain->SetBranchAddress(_paraNames[i]+"_err", err);
	if(b_err) *err = 0.;

        double * pull = new double[1];
        _pulls.push_back(pull);
        int b_pull = fChain->SetBranchAddress(_paraNames[i]+"_pull", pull);
	if(b_pull) *pull = 0.;
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
