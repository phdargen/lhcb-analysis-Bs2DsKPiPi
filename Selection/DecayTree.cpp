#define DecayTree_cxx
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>

using namespace std; 

TTree* DecayTree::GetInputTree(){
    
    TString tupleName;
    if(_decay==Decay::signal)tupleName="Bs2DsKpipi_";
    else tupleName = "Bs2Dspipipi_";
    if(_Ds_finalState==Ds_finalState::pipipi)tupleName+="Ds2pipipi_Tuple/DecayTree";
    if(_Ds_finalState==Ds_finalState::Kpipi)tupleName+="Ds2Kpipi_Tuple/DecayTree";
    else tupleName+="Ds2KKpi_Tuple/DecayTree";    

    TChain* chain = new TChain(tupleName);
    
    TString fileName = "/auto/data/dargent/BsDsKpipi/Stripped/";

    if(_decay==Decay::signal)fileName+="Signal/";
    else fileName += "Norm/";

    if(_data==DataType::data) {
        fileName+= "data/";
        fileName+= _year; 
    }
    else {
        fileName+= "mc/";
        fileName+= _year; 
        fileName+= "U/";
    }
    fileName+= "b2dhhh_*.root"; 

    cout << fileName << endl;
    //chain->Add("/auto/data/dargent/BsDsKpipi/Stripped/Norm/mc/11U/b2dhhh_11.root");
    //chain->Add("/auto/data/dargent/BsDsKpipi/Stripped/Norm/mc/11U/b2dhhh_10.root");
    chain->Add(fileName);
    if(_data!=DataType::data) chain->Add(fileName.ReplaceAll(TString("U/b2hhh"),TString("D/b2hhh")));
    if(chain==0){
        cout << "ERROR: No file found" << endl;
        throw "ERROR";
    }
    
    return (TTree*)chain;
}

DecayTree::DecayTree(Decay::Type decay, Year::Type year, Ds_finalState::Type finalState, DataType::Type dataType ) : 
fChain(0), _decay(decay), _year(year), _Ds_finalState(finalState), _data(dataType) 
{    
    cout << "Requested to process files with options: " << endl << endl;

    TString s1,s2,s3,s4;
    
    if(_data==DataType::data)s1="data";
    else s1 = "mc";
    if(_decay==Decay::signal)s2="signal";
    else s2 = "norm";
    if(_Ds_finalState==Ds_finalState::pipipi)s3="Ds2pipipi";
    if(_Ds_finalState==Ds_finalState::Kpipi)s3="Ds2Kpipi";
    else s3="Ds2KKpi";

    cout << "DataType: " << s1 << endl;
    cout << "Decay: " << s2 << endl;
    cout << "Ds finalState: " << s3 << endl;
    cout << "Year: " << _year << endl << endl;

    _outFileName = "/auto/data/dargent/BsDsKpipi/";    
    _outFileName += "Mini/";
    _outFileName += s1;  
    _outFileName += "/";  
    _outFileName += s2;   
    _outFileName += "_";  
    _outFileName += s3; 
    _outFileName += "_";     
    _outFileName += _year;
    _outFileName += ".root";    
}

Bool_t DecayTree::TriggerCuts(){

    if( (!Bs_L0Global_TIS) && (!Bs_L0HadronDecision_TOS)) return false;
    
    if(_year == 15 || _year == 16){
        if((!Bs_Hlt1TrackMVADecision_TOS) && (!Bs_Hlt1TwoTrackMVADecision_TOS) ) return false;
        if((!Bs_Hlt2Topo2BodyDecision_TOS) &&  (!Bs_Hlt2Topo3BodyDecision_TOS) && (!Bs_Hlt2Topo4BodyDecision_TOS)  
        && (!Bs_Hlt2PhiIncPhiDecision_TOS) ) return false;
    }
    
    else if(_year == 11 || _year == 12){
        if(!Bs_Hlt1TrackAllL0Decision_TOS) return false;
        if((!Bs_Hlt2Topo2BodyBBDTDecision_TOS) &&  (!Bs_Hlt2Topo3BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo4BodyBBDTDecision_TOS) 
        && (!Bs_Hlt2IncPhiDecision_TOS)) return false;
    }
    
    return true;
}

Bool_t DecayTree::LooseCuts(){

    if(Bs_DIRA_OWNPV<0.99994) return false;
    if(Bs_IPCHI2_OWNPV>20) return false;
    if(Bs_FDCHI2_OWNPV<100) return false;
    if((Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)> 8) return false;
    if(Bs_TAU < 0.0002) return false;

    if(Bs_MM < 4800. || Bs_MM > 6000.) return false;
    if(Bs_DTF_M[0] < 4800. || Bs_DTF_M[0] > 6000.) return false;
    if(Bs_BsDTF_M[0] < 4800. || Bs_BsDTF_M[0] > 6000.) return false;
    //if (Bs_PV_M[0] < 4800. || Bs_PV_M[0] > 6000.) return false;
    
    if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < -1) return false;
    if(fabs(Ds_MM-massDs) > 60 ) return false;
    //if(Ds_FDCHI2_ORIVX < 0) return false;

    //if(_data)if(K_plus_PIDK<5) return false;
    
    /*
    if(K_plus_PT < 250) return false;
    if(pi_plus_PT < 250) return false;
    if(pi_minus_PT < 250) return false;
    if(K_plus_fromDs_PT < 250) return false;
    if(K_minus_fromDs_PT < 250) return false;
    if(pi_minus_fromDs_PT < 250) return false;
*/

    return true;
}

void DecayTree::Loop()
{

   time_t startTime = time(0);
  
   Init();
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  // disable all branches
   
   // activate branchname
   fChain->SetBranchStatus("*ID",1);  

   fChain->SetBranchStatus("Bs_*_T*S",1);  
   fChain->SetBranchStatus("Bs_*Muon*_T*S",0);  
   fChain->SetBranchStatus("Bs_*Hlt*Phys*_T*S",0);  
   fChain->SetBranchStatus("Bs_*Hlt*Global*_T*S",0);  
    
   fChain->SetBranchStatus("Bs_TAG*",1);  
   fChain->SetBranchStatus("Bs_*DEC",0);  
   fChain->SetBranchStatus("Bs_*PROB",0);  

   fChain->SetBranchStatus("Bs_*DTF*",1);  
   fChain->SetBranchStatus("Bs_*PV*",1);  

   fChain->SetBranchStatus("Bs*ENDVERTEX*",1);  
   fChain->SetBranchStatus("Bs*OWNPV*",1);  
   fChain->SetBranchStatus("Bs*ENDVERTEX_COV*",0);  
   fChain->SetBranchStatus("Bs*OWNPV_COV*",0);  
       
   fChain->SetBranchStatus("Ds*ENDVERTEX*",1);  
   fChain->SetBranchStatus("Ds*OWNPV*",1);  
   fChain->SetBranchStatus("Ds*ENDVERTEX_COV*",0);  
   fChain->SetBranchStatus("Ds*OWNPV_COV*",0);  
   fChain->SetBranchStatus("Ds*ORIVX*",1);  
   fChain->SetBranchStatus("Ds*ORIVX_COV*",0);  
    
   fChain->SetBranchStatus("*_1_12*ENDVERTEX*",1);  
   fChain->SetBranchStatus("*_1_12*OWNPV*",1);  
   fChain->SetBranchStatus("*_1_12*ENDVERTEX_COV*",0);  
   fChain->SetBranchStatus("*_1_12*OWNPV_COV*",0);  
   fChain->SetBranchStatus("*_1_12*ORIVX*",1);  
   fChain->SetBranchStatus("*_1_12*ORIVX_COV*",0);  

   fChain->SetBranchStatus("*IP*",1);  
   fChain->SetBranchStatus("*IPCHI2*",1);  
   fChain->SetBranchStatus("*FD*",1);  
   fChain->SetBranchStatus("*FDCHI2*",1);  
   fChain->SetBranchStatus("*P",1);  
   fChain->SetBranchStatus("*PT",1);  
   fChain->SetBranchStatus("*PE",1);  
   fChain->SetBranchStatus("*PX",1);  
   fChain->SetBranchStatus("*PY",1);  
   fChain->SetBranchStatus("*PZ",1);  
   fChain->SetBranchStatus("*ETA",1);  
   fChain->SetBranchStatus("*MM*",1);  
   fChain->SetBranchStatus("*TAU*",1);  
   fChain->SetBranchStatus("*ptasy_1.00",1);  

   fChain->SetBranchStatus("*DIRA*",1);  
   fChain->SetBranchStatus("Ds_DOCA*",1);  
   fChain->SetBranchStatus("*_1_12*_DOCA*",1);     
    
   fChain->SetBranchStatus("*PID*",1);  
   fChain->SetBranchStatus("*PIDe*",0);  
   fChain->SetBranchStatus("*ProbNN*",1);  
   fChain->SetBranchStatus("*ProbNNe*",0);  
   fChain->SetBranchStatus("*TRACK_Ghost*",1);  
   fChain->SetBranchStatus("*TRACK_CHI2*",1);  
   fChain->SetBranchStatus("*isMuon*",1);  

   fChain->SetBranchStatus("nCandidate",1) ;
   fChain->SetBranchStatus("nTracks",1) ;
   fChain->SetBranchStatus("nPV",1) ;
   fChain->SetBranchStatus("eventNumber",1) ;
   fChain->SetBranchStatus("runNumber",1) ;
   fChain->SetBranchStatus("EventInSequence",1) ;
   fChain->SetBranchStatus("totCandidates",1) ;
   fChain->SetBranchStatus("Polarity",1) ;

   TFile* output = new TFile(_outFileName,"RECREATE");
   TTree* summary_tree = fChain->CloneTree(0);
    
   Long64_t nentries = fChain->GetEntries();
   cout << "Have " << nentries << " events" <<  endl << endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if(0ul == (jentry % 10000ul)) cout << "Read event " << jentry << "/" << nentries << endl;
      fChain->GetEntry(jentry);   

      //if(!TriggerCuts()) continue;
      //if(!LooseCuts()) continue;

      summary_tree->Fill();
      if(0ul == (jentry % 100000ul))summary_tree->AutoSave();
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
