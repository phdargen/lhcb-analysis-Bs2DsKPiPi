#define MiniDecayTree_cxx
#include "MiniDecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>

using namespace std; 


TTree* MiniDecayTree::GetInputTree(){
        
    cout << "Read from file " << _inFileName << endl;    
        
    TChain* chain= new TChain("DecayTree");
    chain->Add(_inFileName);
    //chain->Add("Signal/11D/*.root");
    //chain->Add("Signal/12U/*.root");
    //chain->Add("Signal/12D/*.root");
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

    for (Long64_t i=0; i<nentries;i++) {
        
        if(0ul == (i % 10000ul)) cout << "Read event " << i << "/" << nentries <<
        "  ( " << i/(double)nentries << " % )" << endl;

        fChain->GetEntry(i);   
        
        Long64_t j = LoadTree(i);
        if (j < 0) break;
        
        //if(!TriggerCuts(j)) continue;
        //if(!LooseCuts(j)) continue;
        
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
