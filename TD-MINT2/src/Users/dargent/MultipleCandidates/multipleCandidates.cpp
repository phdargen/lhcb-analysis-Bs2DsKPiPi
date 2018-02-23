#include "TEventList.h"
#include "TPaletteAxis.h"
#include "TProfile.h"
#include <vector>

using namespace std;


void writeMultipleDataCanidatesToFile(TString kind, TString year, double cutBDT, double cutPID){

  //write all the multiple canidated to a file with totCandidates>1 (after dm preselction and nShared cut for muons)                                                                    
  TString nameFout;

  nameFout = "candidateLists/"+year+"_"+kind+"Data.txt";
  std::ofstream fileOut(nameFout);
  TChain*  Chain = new TChain("BDT_Tree");
  //TChain*  preChain = new TChain("BDT_Tree");

  if(year=="2012" || year=="2011") Chain->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/"+year+"_"+kind+"_Run1BDT.root");
  else Chain->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/"+year+"_"+kind+"_Run2BDT.root");

  //TTree* Chain= (TTree*)preChain->CopyTree("");
  
  double m_D;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;


  int Dst_BKGCAT;
  ULong64_t totCandidates;

  double mu0_PX, mu0_PY, mu0_PZ;
  double mu1_PX, mu1_PY, mu1_PZ;
  double h0_PX, h0_PY, h0_PZ;
  double h1_PX, h1_PY, h1_PZ;
  double Slowpi_PX,Slowpi_PY,Slowpi_PZ;
  double deltaM;
  int mu1_MuonNShared,mu0_MuonNShared;
  double Slowpi_PT;
  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  bool mu1_Hlt1TrackMuonDecision_TOS, mu0_Hlt1TrackMuonDecision_TOS, D_Hlt1TrackAllL0Decision_TOS;
  bool Hlt2_TOS;
  int Slowpi_ID;

  double h0_ProbNNghost, h1_ProbNNghost,mu0_ProbNNghost, mu1_ProbNNghost, Slowpi_ProbNNghost;
  double h0_ProbNNk,h1_ProbNNk,h0_ProbNNpi,h1_ProbNNpi;

  Chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  Chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  Chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);
  Chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  Chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  Chain->SetBranchAddress("h0_ProbNNk",&h0_ProbNNk);
  Chain->SetBranchAddress("h1_ProbNNk",&h1_ProbNNk);
  Chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  Chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  Chain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  Chain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  Chain->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  Chain->SetBranchAddress("Slowpi_ID",&Slowpi_ID);

  Chain->SetBranchAddress("totCandidates",&totCandidates);
  Chain->SetBranchAddress("nCandidate",&nCandidate);
  Chain->SetBranchAddress("eventNumber",&eventNumber);
  Chain->SetBranchAddress("runNumber",&runNumber);
  Chain->SetBranchAddress("Dst_DTF_D0_M",&m_D);
  Chain->SetBranchAddress("deltaM",&deltaM);

  Chain->SetBranchAddress("BDT",&BDT);
  Chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  Chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);

  Chain->SetBranchAddress("Slowpi_PT",&Slowpi_PT);

  Chain->SetBranchAddress("mu0_PX",&mu0_PX);
  Chain->SetBranchAddress("mu0_PY",&mu0_PY);
  Chain->SetBranchAddress("mu0_PZ",&mu0_PZ);
  Chain->SetBranchAddress("mu1_PX",&mu1_PX);
  Chain->SetBranchAddress("mu1_PY",&mu1_PY);
  Chain->SetBranchAddress("mu1_PZ",&mu1_PZ);
  Chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  Chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);

  Chain->SetBranchAddress("h0_PX",&h0_PX);
  Chain->SetBranchAddress("h0_PY",&h0_PY);
  Chain->SetBranchAddress("h0_PZ",&h0_PZ);
  Chain->SetBranchAddress("h1_PX",&h1_PX);
  Chain->SetBranchAddress("h1_PY",&h1_PY);
  Chain->SetBranchAddress("h1_PZ",&h1_PZ);

                                                                                
  Long64_t nentries = Chain->GetEntries();
  int counter1=0;
  int counter2=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Chain->GetEntry(jentry);

    if(deltaM>146.5 || deltaM < 144.5) continue;
    if(BDT<cutBDT) continue;
    if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
    if(totCandidates==1) {counter1++; continue; }

    counter2++;

    fileOut << eventNumber << "\t" << runNumber << "\t" << m_D << "\t" <<Slowpi_PT << "\t" << nCandidate << "\t" << Slowpi_ID << "\t" << std::endl;

    
  }
  
  std::cout<<"single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl;
}

class Cand {
public :
  ULong64_t eventNumber;
  UInt_t    runNumber;
  Double_t  mass;
  UInt_t    nCandidate;
  double Slowpi_PT;
  int Slowpi_ID;
  Double_t  id() { return ((double) runNumber)/((double) eventNumber); }
};

void chose_multiple_Data_events(TString kind, TString year, double cutBDT, double cutPID, bool writeAllCandidatesToFile=false){

  vector<Cand> events;
  vector<Cand> chosen_events;

  Cand tmp;
  TString pathToFiles="candidateLists/";

  TString filename = pathToFiles+year+"_"+kind+"Data.txt";

  TString outputfilename;
  if(writeAllCandidatesToFile) outputfilename = pathToFiles+year+"_"+kind+"Data_AllMult.txt";
  else outputfilename = pathToFiles+year+"_"+kind+"Data_chosen.txt";
                                                                                    
  cout << "Readingfile "<<filename<<" ... to " <<outputfilename<< endl;
  ifstream rs(filename);
  if (!rs) {
    cout << "Unable to open " << filename << endl;
    return;
  }

  while (rs >> tmp.eventNumber >> tmp.runNumber  >>  tmp.mass >> tmp.Slowpi_PT >> tmp.nCandidate >> tmp.Slowpi_ID ) {
    events.push_back(tmp);
  }
  rs.close();
  cout << "Done: " << events.size() << " events found." << endl;



  ///___________________________________________________________________________________                                                                                                                                                                                                                                                                                                                                                                                                           
  TChain* chain = new TChain("BDT_Tree");
  if(year=="2012" || year=="2011") chain->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/"+year+"_"+kind+"_Run1BDT.root");
  else chain->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/"+year+"_"+kind+"_Run2BDT.root");

  //____________________________________________________________________________________                                                   

  double mD;
  ULong64_t loop_eventNumber;
  UInt_t loop_runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;
  int Slowpi_ID;
  double Slowpi_PT;
 
  Cand loopCand;
  vector<Cand> chosenCandidates;
  vector<Cand> multCandidates;

  double h0_ProbNNghost, h1_ProbNNghost,mu0_ProbNNghost, mu1_ProbNNghost, Slowpi_ProbNNghost;
  double h0_ProbNN,h1_ProbNN;
  double BDT, mu0_ProbNNmu, mu1_ProbNNmu;
  bool mu1_L0MuonDecision_TOS,mu0_L0MuonDecision_TOS;
  bool mu1_Hlt1TrackMuonDecision_TOS, mu0_Hlt1TrackMuonDecision_TOS, D_Hlt1TrackAllL0Decision_TOS;
  bool Hlt2_TOS;
  int mu1_MuonNShared,mu0_MuonNShared;
  double deltaM;

  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&loop_eventNumber);
  chain->SetBranchAddress("runNumber",&loop_runNumber);
  chain->SetBranchAddress("Dst_DTF_D0_M",&mD);
  chain->SetBranchAddress("Slowpi_ID",&Slowpi_ID);
  chain->SetBranchAddress("Slowpi_PT",&Slowpi_PT);

  chain->SetBranchAddress("totCandidates",&totCandidates);
  chain->SetBranchAddress("deltaM",&deltaM);
  chain->SetBranchAddress("mu0_MuonNShared",&mu0_MuonNShared);
  chain->SetBranchAddress("mu1_MuonNShared",&mu1_MuonNShared);
  chain->SetBranchAddress("BDT",&BDT);
  chain->SetBranchAddress("mu0_ProbNNmu",&mu0_ProbNNmu);
  chain->SetBranchAddress("mu1_ProbNNmu",&mu1_ProbNNmu);
  chain->SetBranchAddress("mu0_L0MuonDecision_TOS",&mu0_L0MuonDecision_TOS);
  chain->SetBranchAddress("mu1_L0MuonDecision_TOS",&mu1_L0MuonDecision_TOS);
  chain->SetBranchAddress("mu0_Hlt1TrackMuonDecision_TOS",&mu0_Hlt1TrackMuonDecision_TOS);
  chain->SetBranchAddress("mu1_Hlt1TrackMuonDecision_TOS",&mu1_Hlt1TrackMuonDecision_TOS);
  chain->SetBranchAddress("D_Hlt1TrackAllL0Decision_TOS",&D_Hlt1TrackAllL0Decision_TOS);
  chain->SetBranchAddress("mu0_ProbNNghost",&mu0_ProbNNghost);
  chain->SetBranchAddress("mu1_ProbNNghost",&mu1_ProbNNghost);
  chain->SetBranchAddress("h0_ProbNNghost",&h0_ProbNNghost);
  chain->SetBranchAddress("h1_ProbNNghost",&h1_ProbNNghost);
  chain->SetBranchAddress("Slowpi_ProbNNghost",&Slowpi_ProbNNghost);

  //Loop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  Long64_t nentries = chain->GetEntries();
  
  int counter1=0;int counter2=0;int counter3=0;


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    chain->GetEntry(jentry);
    Dst_BKGCAT=1;
    
    if(deltaM>146.5 || deltaM < 144.5) continue;
    if(BDT<cutBDT) continue;
    if(mu0_ProbNNmu<cutPID||mu1_ProbNNmu<cutPID) continue;
    if(totCandidates==1) continue;

    /// create temporary candidate to campare with candidates from file                                                                                                          
                                                                                                                                                                                                        
    loopCand.mass=mD;
    loopCand.eventNumber=loop_eventNumber;
    loopCand.runNumber=loop_runNumber;
    loopCand.nCandidate=nCandidate;
    loopCand.Slowpi_ID=Slowpi_ID;
    loopCand.Slowpi_PT=Slowpi_PT;

   ///every candidate is put into a vector which will be added by the candidates of the same event in the following                                                                                                                                                                                                                                                                                         
    multCandidates.push_back(loopCand);

    for (int i=1;i<events.size();++i) { ///compare to every candidate from file                                                                                                                                                                                                 
      if(events[i].eventNumber==0) continue; ///checks if candidate has already been matched (to avoid double counting)                 

      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate!=events[i].nCandidate ) {          ///find partner, but not the cand itself! (nCand)
      //if( (TMath::Abs((double)loop_eventNumber-(double)events[i].eventNumber)<1) && (TMath::Abs((double)loop_runNumber-(double)events[i].runNumber)<1) && (TMath::Abs((double)nCandidate-(double)events[i].nCandidate)>1) ) {  
	//cout<<jentry<<"  "<<i<<"  "<<nCandidate<<"  "<<events[i].nCandidate <<endl;               
	multCandidates.push_back(events[i]); /// if partner is found put into vector                                                    
	events[i].eventNumber=0; ///candidate has already found its partners, so set event number to 0                               
      }
      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate==events[i].nCandidate ) {  ///find itself and "delete" from list by setting eventnumber to 0     
	events[i].eventNumber=0;
      }

    }
    //cout<<"event "<<jentry<<endl;
    int nCand=multCandidates.size();
    //cout<<"nCand "<<nCand<<endl;
    counter1+=nCand;
    if(nCand>1) { ///multiple candidate found if in the vector multCand has more than one entry, chose one randomly                                                   
      //cout<<"in Loop"<<endl;
      counter2+=1;
      TRandom3 generator(loop_eventNumber); ///seed is event number                                                                                             
      double randomNr=generator.Rndm();
      int index=int(nCand*randomNr);
      //cout<<"chosen index "<<index<<" out of "<< nCand<<endl;                                                                                                               
      int littleCounter=0;
      for(int k=0;k<multCandidates.size();++k){
	//if one wants to have a closer look to all candidates, write them all to the file
	if(writeAllCandidatesToFile) chosenCandidates.push_back(multCandidates[k]);	
	else {
	 if(k!=index) {chosenCandidates.push_back(multCandidates[index]); ///chosen candidates contains the candidates to be removed from the events in the loop   
	  //std::cout<<"reject "<<k<<std::endl;
	  counter3+=1;
	  littleCounter+=1;}
	}
      }
      //cout<<"removed "<<littleCounter<<"  "<<littleCounter+1-nCand<<endl; 
    }

    multCandidates.clear(); ///clear vector for next event   

  }                                                                                                                                     
  ///write the chosen candidates to an output file                                                                              

  cout << "Writing output file..." << endl;
  ofstream fout(outputfilename);
  for (vector<Cand>::iterator it=chosenCandidates.begin(); it<chosenCandidates.end(); ++it) {
    fout << (*it).eventNumber << "\t" << (*it).runNumber << "\t" <<(*it).mass<< "\t" << (*it).Slowpi_PT<< "\t" <<(*it).nCandidate << "\t" << (*it).Slowpi_ID<<endl;
  }
  fout.close();

  cout << "Done. Removed candidates: " << chosenCandidates.size()<<endl;
  //cout << counter1 <<"  "<<counter2<<"  "<<counter3<< endl;
}

bool MatchToMultipleCandidate(vector<Cand> chosen_cand_list, ULong64_t eventNumber, UInt_t runNumber,UInt_t nCandidate){

  for (vector<Cand>::iterator it = chosen_cand_list.begin(); it < chosen_cand_list.end(); ++it) {
    if (eventNumber != (*it).eventNumber) continue;
    if (runNumber != (*it).runNumber) continue;
    if (nCandidate != (*it).nCandidate) continue;
    chosen_cand_list.erase(it);
    return true;
  }
  return false;
}


void addMultipleCandidateBranchData(TChain *chain, TString fOut,const char *f_chosenCand){

  //as written above, this one looks for the multiple candidates which are to be rejected and only saves the ones which are selected                                                                                                                       
  vector<Cand> matched_list;
  Cand event;
  std::ifstream file(f_chosenCand);

  if ( !file ) {
    std::cout << "Unable to open " << f_chosenCand << std::endl;
  } else
    std::cout << "Reading matched candidates from " << f_chosenCand  << "..." << std::endl;

  while ( (file >>  event.eventNumber >>  event.runNumber >> event.mass >>  event.Slowpi_PT  >>event.nCandidate  >> event.Slowpi_ID)) {
    matched_list.push_back(event);
  }
  file.close();
  std::cout << "Done." << std::endl;

  double mD;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  int Dst_BKGCAT;
  ULong64_t totCandidates;
  double Slowpi_PT;
  int Slowpi_ID;

  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&eventNumber);
  chain->SetBranchAddress("runNumber",&runNumber);
  chain->SetBranchAddress("Dst_DTF_D0_M",&mD);
  chain->SetBranchAddress("totCandidates",&totCandidates);
  chain->SetBranchAddress("Slowpi_PT",&Slowpi_PT);
  chain->SetBranchAddress("Slowpi_ID",&Slowpi_ID);

  TFile* targetFile = new TFile(fOut,"RECREATE");

  TTree* new_tree = chain->CopyTree("");
  bool isRejectedMultipleCandidate;
  TBranch* Bra = new_tree->Branch("isRejectedMultipleCandidate",&isRejectedMultipleCandidate);

  int numEvents = chain ->GetEntries();

  for(int i=0; i< numEvents; i++){

    chain->GetEntry(i);
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;

    if( MatchToMultipleCandidate(matched_list, eventNumber, runNumber, nCandidate) ) isRejectedMultipleCandidate=true;  // if is found on list, reject!                                                                                                    
    else isRejectedMultipleCandidate = false;
    Bra->Fill();
  }

  //here we only take the selected ones!!!!                                                                                                                                                                                                                
  TTree* new_tree2 = new_tree->CopyTree("!isRejectedMultipleCandidate");

  new_tree2->Write();
  targetFile->Write();
  targetFile->Close();
  //delete new_tree;                                                                                                                                                                                                                                       
  delete targetFile;

}


//DATA!!!                                                                                                                                                                                                                                                        
void createDataTuplesWithSelectedMultipleCand(){

  std::cout<<"creating tuples without multiple candidates..."<<std::endl;

  TChain *chain_pipi_12 = new TChain("BDT_Tree");
  TChain *chain_pipi_15 = new TChain("BDT_Tree");
  TChain *chain_pipi_16 = new TChain("BDT_Tree");

  chain_pipi_12->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2012_D2pipimumu_BDT.root");
  chain_pipi_15->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2015_D2pipimumu_Run2BDT.root");
  chain_pipi_16->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2016_D2pipimumu_Run2BDT.root");

  std::cout<<"...pipimumu 2012..."<<std::endl;
  addMultipleCandidateBranchData(chain_pipi_12,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2012_D2pipimumu_BDT_noMultCand.root","candidateLists/2012_D2pipimumuData_chosen.txt");
  std::cout<<"...pipimumu 2015..."<<std::endl;
  addMultipleCandidateBranchData(chain_pipi_15,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2015_D2pipimumu_Run2BDT_noMultCand.root","candidateLists/2015_D2pipimumuData_chosen.txt");
  std::cout<<"...pipimumu 2016..."<<std::endl;
  addMultipleCandidateBranchData(chain_pipi_16,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2016_D2pipimumu_Run2BDT_noMultCand.root","candidateLists/2016_D2pipimumuData_chosen.txt");



  TChain *chain_KK_12 = new TChain("BDT_Tree");
  TChain *chain_KK_15 = new TChain("BDT_Tree");
  TChain *chain_KK_16 = new TChain("BDT_Tree");

  chain_KK_12->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2012_D2KKmumu_BDT.root");
  chain_KK_15->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2015_D2KKmumu_Run2BDT.root");
  chain_KK_16->AddFile("/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2016_D2KKmumu_Run2BDT.root");

  std::cout<<"...KKmumu 2015..."<<std::endl;
  addMultipleCandidateBranchData(chain_KK_12,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2012_D2KKmumu_BDT_noMultCand.root","candidateLists/2012_D2KKmumuData_chosen.txt");
  std::cout<<"...KKmumu 2015..."<<std::endl;
  addMultipleCandidateBranchData(chain_KK_15,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2015_D2KKmumu_Run2BDT_noMultCand.root","candidateLists/2015_D2KKmumuData_chosen.txt");
  std::cout<<"...KKmumu 2016..."<<std::endl;
  addMultipleCandidateBranchData(chain_KK_16,"/auto/data/mitzel/D2hhmumu/new/AsymmetryMeasurements/2016_D2KKmumu_Run2BDT_noMultCand.root","candidateLists/2016_D2KKmumuData_chosen.txt");

}






