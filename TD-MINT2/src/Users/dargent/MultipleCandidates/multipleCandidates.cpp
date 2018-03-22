#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "TEventList.h"
#include "TPaletteAxis.h"
#include "TProfile.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TString.h"
#include <vector>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
using namespace std;
using namespace MINT;

TFile* output = 0;
TTree* out_tree = 0;

class Cand {
public :
  ULong64_t eventNumber;
  UInt_t    runNumber;
  Double_t  mass;
  UInt_t    nCandidate;
  double Ds_PT;
  int Ds_ID;
  Double_t  id() { return ((double) runNumber)/((double) eventNumber); }
};

void writeMultipleDataCanidatesToFile(TString listName = ""){
  /// Write all the multiple canidated to a file with totCandidates>1 
  TString nameFout = "candidateLists/"+listName+".txt";

  ofstream fileOut;
  fileOut.open(nameFout,ios_base::out );

  double m_Bs;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  ULong64_t totCandidates;
  int Ds_ID;
  double Ds_PT;

  out_tree->SetBranchAddress("Ds_ID",&Ds_ID);
  out_tree->SetBranchAddress("totCandidates",&totCandidates);
  out_tree->SetBranchAddress("nCandidate",&nCandidate);
  out_tree->SetBranchAddress("eventNumber",&eventNumber);
  out_tree->SetBranchAddress("runNumber",&runNumber);
  out_tree->SetBranchAddress("Bs_DTF_MM",&m_Bs);
  out_tree->SetBranchAddress("Ds_PT",&Ds_PT);

  out_tree->SetBranchStatus("*",0);
  out_tree->SetBranchStatus("Ds_ID",1);
  out_tree->SetBranchStatus("totCandidates",1);
  out_tree->SetBranchStatus("nCandidate",1);
  out_tree->SetBranchStatus("eventNumber",1);
  out_tree->SetBranchStatus("runNumber",1);
  out_tree->SetBranchStatus("Bs_DTF_MM",1);
  out_tree->SetBranchStatus("Ds_PT",1);

  Long64_t nentries = out_tree->GetEntries();
  int counter1=0;
  int counter2=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (0ul == (jentry % 10000ul)) std::cout << "Read event " << jentry << "/" << nentries << std::endl;

    out_tree->GetEntry(jentry);
    if(totCandidates==1) {counter1++; continue; }

    counter2++;
    fileOut << eventNumber << "\t" << runNumber << "\t" << m_Bs << "\t" <<Ds_PT << "\t" << nCandidate << "\t" << Ds_ID << "\t" << std::endl;
    
  }
  fileOut.close();
  cout<< endl << "single candidate events = "<<counter1<< endl << "events with potential multiple candidates = "<<counter2<< endl << "sum = "<< counter1+counter2<<endl;
}

void chose_multiple_Data_events(TString listName = "", bool writeAllCandidatesToFile=false){

  vector<Cand> events;
  vector<Cand> chosen_events;

  Cand tmp;
  TString filename = "candidateLists/"+listName+".txt";

  TString outputfilename;
  if(writeAllCandidatesToFile) outputfilename = "candidateLists/"+listName+"_AllMult.txt";
  else outputfilename = "candidateLists/"+listName+"_chosen.txt";

  cout << "Readingfile "<<filename<<" ... to " <<outputfilename<< endl;

  ifstream rs(filename);
  if (!rs) {
    cout << "Unable to open " << filename << endl;
    return;
  }

  while (rs >> tmp.eventNumber >> tmp.runNumber  >>  tmp.mass >> tmp.Ds_PT >> tmp.nCandidate >> tmp.Ds_ID ) {
    events.push_back(tmp);
  }
  rs.close();
  cout << "Done: " << events.size() << " events found." << endl;

  ///___________________________________________________________________________________
  Cand loopCand;
  vector<Cand> chosenCandidates;
  vector<Cand> multCandidates;

  double m_Bs;
  ULong64_t loop_eventNumber;
  UInt_t loop_runNumber;
  UInt_t nCandidate;
  ULong64_t totCandidates;
  int Ds_ID;
  double Ds_PT;

  out_tree->SetBranchAddress("Ds_ID",&Ds_ID);
  out_tree->SetBranchAddress("totCandidates",&totCandidates);
  out_tree->SetBranchAddress("nCandidate",&nCandidate);
  out_tree->SetBranchAddress("eventNumber",&loop_eventNumber);
  out_tree->SetBranchAddress("runNumber",&loop_runNumber);
  out_tree->SetBranchAddress("Bs_DTF_MM",&m_Bs);
  out_tree->SetBranchAddress("Ds_PT",&Ds_PT);

  out_tree->SetBranchStatus("*",0);
  out_tree->SetBranchStatus("Ds_ID",1);
  out_tree->SetBranchStatus("totCandidates",1);
  out_tree->SetBranchStatus("nCandidate",1);
  out_tree->SetBranchStatus("eventNumber",1);
  out_tree->SetBranchStatus("runNumber",1);
  out_tree->SetBranchStatus("Bs_DTF_MM",1);
  out_tree->SetBranchStatus("Ds_PT",1);

  Long64_t nentries = out_tree->GetEntries();
  int counter1=0;int counter2=0;int counter3=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if (0ul == (jentry % 10000ul)) std::cout << "Read event " << jentry << "/" << nentries << std::endl;

    out_tree->GetEntry(jentry);
    if(totCandidates==1) continue;

    /// create temporary candidate to compare with candidates from file
    loopCand.mass=m_Bs;
    loopCand.eventNumber=loop_eventNumber;
    loopCand.runNumber=loop_runNumber;
    loopCand.nCandidate=nCandidate;
    loopCand.Ds_ID=Ds_ID;
    loopCand.Ds_PT=Ds_PT;

    ///every candidate is put into a vector which will be added by the candidates of the same event in the following
    multCandidates.push_back(loopCand);

    for (int i=1;i<events.size();++i) { ///compare to every candidate from file

      if(events[i].eventNumber==0) continue; ///checks if candidate has already been matched (to avoid double counting)

      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate!=events[i].nCandidate ) {
	///find partner, but not the cand itself! (nCand)
      //if( (TMath::Abs((double)loop_eventNumber-(double)events[i].eventNumber)<1) && (TMath::Abs((double)loop_runNumber-(double)events[i].runNumber)<1) && (TMath::Abs((double)nCandidate-(double)events[i].nCandidate)>1) ) {  
	//cout<<jentry<<"  "<<i<<"  "<<nCandidate<<"  "<<events[i].nCandidate <<endl; 
	multCandidates.push_back(events[i]); /// if partner is found put into vector 
	events[i].eventNumber=0; ///candidate has already found its partners, so set event number to 0
      }
      if(loop_eventNumber==events[i].eventNumber && loop_runNumber==events[i].runNumber && nCandidate==events[i].nCandidate ) {  
	///find itself and "delete" from list by setting eventnumber to 0     
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
	 if(k!=index) {chosenCandidates.push_back(multCandidates[k]); ///chosen candidates contains the candidates to be removed from the events in the loop   
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
  ofstream fout;
  fout.open(outputfilename, ios_base::out);
  for (vector<Cand>::iterator it=chosenCandidates.begin(); it<chosenCandidates.end(); ++it) {
    fout << (*it).eventNumber << "\t" << (*it).runNumber << "\t" <<(*it).mass<< "\t" << (*it).Ds_PT<< "\t" <<(*it).nCandidate << "\t" << (*it).Ds_ID<<endl;
  }
  fout.close();

  cout << endl << "Removed " << counter3 << " events out of " << nentries << " ( " << counter3/(double)nentries * 100. << " %) " << endl; 
  cout << "Remaining events " << nentries-counter3 << " ( " << (nentries-counter3)/(double)nentries * 100. << " %) " << endl; 
  cout << "Fraction of events with multiple candidates = " << counter2/((double)nentries) * 100. << " %" << endl << endl;
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

void addMultipleCandidateBranchData(TString listName){

  /// As written above, this one looks for the multiple candidates which are to be rejected and only saves the ones which are selected
  vector<Cand> matched_list;
  Cand event;

  TString f_chosenCand = "candidateLists/"+listName+"_chosen.txt";
  std::ifstream file(f_chosenCand);

  if( !file ) std::cout << "Unable to open " << f_chosenCand << std::endl;
  else std::cout << "Reading matched candidates from " << f_chosenCand  << "..." << std::endl;

  while ( (file >>  event.eventNumber >>  event.runNumber >> event.mass >>  event.Ds_PT  >>event.nCandidate  >> event.Ds_ID)) {
    	matched_list.push_back(event);
  }
  file.close();
  std::cout << "Done." << std::endl;

  double m_Bs;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;
  ULong64_t totCandidates;
  int Ds_ID;
  double Ds_PT;

  out_tree->SetBranchAddress("Ds_ID",&Ds_ID);
  out_tree->SetBranchAddress("totCandidates",&totCandidates);
  out_tree->SetBranchAddress("nCandidate",&nCandidate);
  out_tree->SetBranchAddress("eventNumber",&eventNumber);
  out_tree->SetBranchAddress("runNumber",&runNumber);
  out_tree->SetBranchAddress("Bs_DTF_MM",&m_Bs);
  out_tree->SetBranchAddress("Ds_PT",&Ds_PT);

  bool isRejectedMultipleCandidate;
  TBranch* Bra = out_tree->Branch("isRejectedMultipleCandidate",&isRejectedMultipleCandidate);

  int counter = 0;
  int numEvents = out_tree->GetEntries();
  for(int i=0; i< numEvents; i++){
    out_tree->GetEntry(i);
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;

    if( MatchToMultipleCandidate(matched_list, eventNumber, runNumber, nCandidate) ) isRejectedMultipleCandidate=true;  // if is found on list, reject!
    else isRejectedMultipleCandidate = false;
    if(isRejectedMultipleCandidate == true)counter++;
    Bra->Fill();
  }

  out_tree->SetBranchStatus("*",1);
  out_tree->Write();
  output->Write();
  output->Close();
  delete output;
}

int main(int argc, char** argv){

  /// Options 
  NamedParameter<string> inFileName("inFileName", (string)"/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");
  NamedParameter<string> outFileName("outFileName", (string)"/auto/data/dargent/BsDsKpipi/BDT/Data/signal_MultCand.root");
  NamedParameter<string> cut("cut",(string)"");
  NamedParameter<int> makeNewList("makeNewList", 0);
  NamedParameter<string> listName("listName",(string)"");
  NamedParameter<int> matchList("matchList", 0);

  cout <<"Searching for multiple candidates in file "<< (string)inFileName <<  endl;
  /// Load file
  TFile *file = new TFile(((string)inFileName).c_str());
  TTree* tree = (TTree*) file->Get("DecayTree");	
  output = new TFile(((string)outFileName).c_str(),"RECREATE");
  out_tree = tree->CopyTree(((string)cut).c_str());

  if(makeNewList == 1) writeMultipleDataCanidatesToFile(TString((string)listName));
  
  if(matchList == 1) chose_multiple_Data_events(TString((string)listName));

  addMultipleCandidateBranchData(TString((string)listName));
   
  cout << endl << "Created new file "<< (string)outFileName <<  endl;

return 0;
}
