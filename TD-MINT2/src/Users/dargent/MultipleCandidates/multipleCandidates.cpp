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
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
#include "Mint/RooHILLdini.h"
#include "Mint/RooHORNSdini.h"

using namespace std;
using namespace MINT;


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


void writeMultipleDataCanidatesToFile(TString kind, int Year, double cutBDT){

  //write all the multiple canidated to a file with totCandidates>1 (after dm preselction and nShared cut for muons)                                                                    
  TString nameFout;
  TString StringYear;

  if(Year == 11) StringYear = "2011";
  if(Year == 12) StringYear = "2012";
  if(Year == 15) StringYear = "2015";
  if(Year == 16) StringYear = "2016";
  nameFout = "candidateLists/"+StringYear+"_"+kind+"Data.txt";
  ofstream fileOut;
  fileOut.open(nameFout,ios_base::out );
  TChain*  Chain = new TChain("DecayTree");

  if(kind == "norm") Chain->AddFile("/auto/data/dargent/BsDsKpipi/BDT/Data/norm.root");
  else Chain->AddFile("/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");

  double m_Bs;
  ULong64_t eventNumber;
  UInt_t runNumber;
  UInt_t nCandidate;

  ULong64_t totCandidates;


  int year;
  int Ds_ID;
  double Ds_PT;

  float BDTG_response;

  Chain->SetBranchAddress("Ds_ID",&Ds_ID);

  Chain->SetBranchAddress("totCandidates",&totCandidates);
  Chain->SetBranchAddress("nCandidate",&nCandidate);
  Chain->SetBranchAddress("eventNumber",&eventNumber);
  Chain->SetBranchAddress("runNumber",&runNumber);


  Chain->SetBranchAddress("Bs_DTF_MM",&m_Bs);

  Chain->SetBranchAddress("BDTG_response",&BDTG_response);
  Chain->SetBranchAddress("year",&year);

  Chain->SetBranchAddress("Ds_PT",&Ds_PT);

                                                                                
  Long64_t nentries = Chain->GetEntries();
  int counter1=0;
  int counter2=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (0ul == (jentry % 10000ul)) std::cout << "Read event " << jentry << "/" << nentries << std::endl;

    Chain->GetEntry(jentry);

    if(year != Year) continue;
    if(BDTG_response<cutBDT) continue;
    if(totCandidates==1) {counter1++; continue; }

    counter2++;

    fileOut << eventNumber << "\t" << runNumber << "\t" << m_Bs << "\t" <<Ds_PT << "\t" << nCandidate << "\t" << Ds_ID << "\t" << std::endl;

    
  }
  fileOut.close();
  std::cout<<"single candidate events "<<counter1<<" candidates in multiple cand. "<<counter2<<" sum "<< counter1+counter2<<std::endl;
}

void chose_multiple_Data_events(TString kind, int Year, double cutBDT, bool writeAllCandidatesToFile=false){

  vector<Cand> events;
  vector<Cand> chosen_events;

  Cand tmp;
  TString pathToFiles="candidateLists/";

  TString outputfilename;
  TString StringYear;

  if(Year == 11) StringYear = "2011";
  if(Year == 12) StringYear = "2012";
  if(Year == 15) StringYear = "2015";
  if(Year == 16) StringYear = "2016";

  TString filename = pathToFiles+StringYear+"_"+kind+"Data.txt";

  if(writeAllCandidatesToFile) outputfilename = pathToFiles+StringYear+"_"+kind+"Data_AllMult.txt";
  else outputfilename = pathToFiles+StringYear+"_"+kind+"Data_chosen.txt";
                                                                                    
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
  TChain* chain = new TChain("DecayTree");

  if(kind == "norm") chain->AddFile("/auto/data/dargent/BsDsKpipi/BDT/Data/norm.root");
  else chain->AddFile("/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");

  //____________________________________________________________________________________                                                   

  double m_Bs;
  ULong64_t loop_eventNumber;
  UInt_t loop_runNumber;
  UInt_t nCandidate;

  ULong64_t totCandidates;

  int year;
  int Ds_ID;
  double Ds_PT;

  float BDTG_response;
 
  Cand loopCand;
  vector<Cand> chosenCandidates;
  vector<Cand> multCandidates;

  chain->SetBranchAddress("Ds_ID",&Ds_ID);

  chain->SetBranchAddress("totCandidates",&totCandidates);
  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&loop_eventNumber);
  chain->SetBranchAddress("runNumber",&loop_runNumber);


  chain->SetBranchAddress("Bs_DTF_MM",&m_Bs);
  chain->SetBranchAddress("BDTG_response",&BDTG_response);
  chain->SetBranchAddress("year",&year);
  chain->SetBranchAddress("Ds_PT",&Ds_PT);

  //Loop                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
  Long64_t nentries = chain->GetEntries();
  
  int counter1=0;int counter2=0;int counter3=0;


  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if (0ul == (jentry % 10000ul)) std::cout << "Read event " << jentry << "/" << nentries << std::endl;

    chain->GetEntry(jentry);

    if(year != Year) continue;
    if(BDTG_response<cutBDT) continue;
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
  ofstream fout;
  fout.open(outputfilename, ios_base::out);
  for (vector<Cand>::iterator it=chosenCandidates.begin(); it<chosenCandidates.end(); ++it) {
    fout << (*it).eventNumber << "\t" << (*it).runNumber << "\t" <<(*it).mass<< "\t" << (*it).Ds_PT<< "\t" <<(*it).nCandidate << "\t" << (*it).Ds_ID<<endl;
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




void addMultipleCandidateBranchData(int Year, TChain *chain, TString fOut,const char *f_chosenCand){

  //as written above, this one looks for the multiple candidates which are to be rejected and only saves the ones which are selected                                                                                    

  vector<Cand> matched_list;
  Cand event;
  std::ifstream file(f_chosenCand);

  if ( !file ) {
    std::cout << "Unable to open " << f_chosenCand << std::endl;
  } else
    std::cout << "Reading matched candidates from " << f_chosenCand  << "..." << std::endl;

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
  int year;

  chain->SetBranchAddress("year",&year);
  chain->SetBranchAddress("Ds_ID",&Ds_ID);
  chain->SetBranchAddress("totCandidates",&totCandidates);
  chain->SetBranchAddress("nCandidate",&nCandidate);
  chain->SetBranchAddress("eventNumber",&eventNumber);
  chain->SetBranchAddress("runNumber",&runNumber);
  chain->SetBranchAddress("Bs_DTF_MM",&m_Bs);
  chain->SetBranchAddress("Ds_PT",&Ds_PT);


  TFile* targetFile = new TFile(fOut,"RECREATE");

  TString CutYear;

  if(Year == 11) CutYear = "year == 11";
  if(Year == 12) CutYear = "year == 12";
  if(Year == 15) CutYear = "year == 15";
  if(Year == 16) CutYear = "year == 16";

  TTree* new_tree = chain->CopyTree(CutYear);
  bool isRejectedMultipleCandidate;
  TBranch* Bra = new_tree->Branch("isRejectedMultipleCandidate",&isRejectedMultipleCandidate);

  int numEvents = chain ->GetEntries();

  for(int i=0; i< numEvents; i++){

    chain->GetEntry(i);
    if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
    if(year != Year) continue;

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

int main(int argc, char** argv){

  //set parameters
  NamedParameter<string> Channel("Channel", (std::string) "norm" , (char*) 0);
  NamedParameter<string> outFileLocation("outFileLocation", (std::string) " " , (char*) 0);
  NamedParameter<float> cutBDT("cutBDT",0.);
  NamedParameter<int> makeNewList("makeNewList", 0);
  NamedParameter<int> matchList("matchList", 0);

  cout <<"creating tuples without multiple candidates..."<< endl;

  TChain *chain = new TChain("DecayTree");
  TString path;
  string channel = Channel;
  string OutLocation = outFileLocation;

  string norm_out11 = "norm_11_noMultCand.root";
  string norm_out12 = "norm_12_noMultCand.root";
  string norm_out15 = "norm_15_noMultCand.root";
  string norm_out16 = "norm_16_noMultCand.root";

  string signal_out11 = "signal_11_noMultCand.root";
  string signal_out12 = "signal_12_noMultCand.root";
  string signal_out15 = "signal_15_noMultCand.root";
  string signal_out16 = "signal_16_noMultCand.root";

  if(channel == "norm") path = "/auto/data/dargent/BsDsKpipi/BDT/Data/norm.root";
  else path = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root";

  chain->AddFile(path);

  if(makeNewList == 1){
  writeMultipleDataCanidatesToFile(channel.c_str(), 11, cutBDT);
  writeMultipleDataCanidatesToFile(channel.c_str(), 12, cutBDT);
  writeMultipleDataCanidatesToFile(channel.c_str(), 15, cutBDT);
  writeMultipleDataCanidatesToFile(channel.c_str(), 16, cutBDT);
  }

  if(matchList == 1){
  chose_multiple_Data_events(channel.c_str(), 11, cutBDT);
  chose_multiple_Data_events(channel.c_str(), 12, cutBDT);
  chose_multiple_Data_events(channel.c_str(), 15, cutBDT);
  chose_multiple_Data_events(channel.c_str(), 16, cutBDT);
  }

  if(channel == "norm"){
  	addMultipleCandidateBranchData(11, chain,(OutLocation.append(norm_out11)).c_str(),"candidateLists/2011_normData_chosen.txt");
  	addMultipleCandidateBranchData(12, chain,(OutLocation.append(norm_out12)).c_str(),"candidateLists/2012_normData_chosen.txt");
  	addMultipleCandidateBranchData(15, chain,(OutLocation.append(norm_out15)).c_str(),"candidateLists/2015_normData_chosen.txt");
  	addMultipleCandidateBranchData(16, chain,(OutLocation.append(norm_out16)).c_str(),"candidateLists/2016_normData_chosen.txt");
  }

  else{
  	addMultipleCandidateBranchData(11, chain,(OutLocation.append(signal_out11)).c_str(),"candidateLists/2011_signalData_chosen.txt");
  	addMultipleCandidateBranchData(12, chain,(OutLocation.append(signal_out12)).c_str(),"candidateLists/2012_signalData_chosen.txt");
  	addMultipleCandidateBranchData(15, chain,(OutLocation.append(signal_out15)).c_str(),"candidateLists/2015_signalData_chosen.txt");
  	addMultipleCandidateBranchData(16, chain,(OutLocation.append(signal_out16)).c_str(),"candidateLists/2016_signalData_chosen.txt");
  }

return 0;
}
