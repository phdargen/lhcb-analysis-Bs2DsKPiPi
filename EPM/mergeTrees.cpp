#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;

void mergeTrees(string old_fileName , string tagging_fileName_Run1, string tagging_fileName_Run2){

	cout << "Adding OS combo branch to file" << old_fileName << ".root" << endl;

	string outputName = old_fileName;
	outputName.append("_tagged.root"); 

	string oldFileName = old_fileName.append(".root");
	TFile* old_file= new TFile(oldFileName.c_str());
	TTree* old_tree = (TTree*) old_file->Get("DecayTree");
	old_tree->SetBranchStatus("OS_Combination_*",0);
	
	string taggingFileName_Run1 = tagging_fileName_Run1.append(".root");
	TFile* tagging_file_Run1= new TFile(taggingFileName_Run1.c_str());
	TTree* tagging_tree_Run1 = (TTree*) tagging_file_Run1->Get("TaggingTree");

	string taggingFileName_Run2 = tagging_fileName_Run2.append(".root");
	TFile* tagging_file_Run2= new TFile(taggingFileName_Run2.c_str());
	TTree* tagging_tree_Run2 = (TTree*) tagging_file_Run2->Get("TaggingTree");

	Int_t OS_Combination_DEC;
	Double_t OS_Combination_ETA;
	//Int_t Combination_DEC;
	//Double_t Combination_ETA;
	Int_t OS_Muon_DEC;
	Double_t OS_Muon_ETA;

	tagging_tree_Run1->SetBranchAddress("OS_Combination_DEC" , &OS_Combination_DEC );
	tagging_tree_Run1->SetBranchAddress("OS_Combination_ETA" , &OS_Combination_ETA );
	tagging_tree_Run1->SetBranchAddress("OS_Muon_DEC" , &OS_Muon_DEC );
	tagging_tree_Run1->SetBranchAddress("OS_Muon_ETA" , &OS_Muon_ETA );
	//tagging_tree -> SetBranchAddress("Combination_DEC" , &Combination_DEC );
	//tagging_tree -> SetBranchAddress("Combination_ETA" , &Combination_ETA );
	tagging_tree_Run2->SetBranchAddress("OS_Combination_DEC" , &OS_Combination_DEC );
	tagging_tree_Run2->SetBranchAddress("OS_Combination_ETA" , &OS_Combination_ETA );
	tagging_tree_Run2->SetBranchAddress("OS_Muon_DEC" , &OS_Muon_DEC );
	tagging_tree_Run2->SetBranchAddress("OS_Muon_ETA" , &OS_Muon_ETA );

	TFile* output=new TFile(outputName.c_str(),"RECREATE");
	TTree* new_tree_Run1 = old_tree->CopyTree("run == 1");
	TTree* new_tree_Run2 = old_tree->CopyTree("run == 2");
	
	Int_t Bs_OS_Muon_DEC;
	Double_t Bs_OS_Muon_ETA;

	TBranch* b_OS_Muon_DEC_Run1; 
	TBranch* b_OS_Muon_PROB_Run1;
	TBranch* b_OS_Muon_DEC_Run2;
	TBranch* b_OS_Muon_PROB_Run2;

	new_tree_Run1->SetBranchAddress("OS_Muon_DEC" , &Bs_OS_Muon_DEC, &b_OS_Muon_DEC_Run1 );
	new_tree_Run1->SetBranchAddress("OS_Muon_PROB" , &Bs_OS_Muon_ETA, &b_OS_Muon_PROB_Run1);
	new_tree_Run2->SetBranchAddress("OS_Muon_DEC" , &Bs_OS_Muon_DEC, &b_OS_Muon_DEC_Run2);
	new_tree_Run2->SetBranchAddress("OS_Muon_PROB" , &Bs_OS_Muon_ETA, &b_OS_Muon_PROB_Run2);

	Int_t Bs_OS_Combination_DEC;
	Double_t Bs_OS_Combination_PROB;
	//Int_t Bs_Combination_DEC;
	//Double_t Bs_Combination_PROB;

	TBranch* OS_Comb_DEC_Branch_Run1 = new_tree_Run1->Branch("OS_Combination_DEC", &Bs_OS_Combination_DEC, "OS_Combination_DEC/I");
	TBranch* OS_Comb_PROB_Branch_Run1 = new_tree_Run1->Branch("OS_Combination_PROB", &Bs_OS_Combination_PROB, "OS_Combination_PROB/D");
	//TBranch* Comb_DEC_Branch = new_tree->Branch("Bs_Combination_DEC", &Bs_Combination_DEC, "Bs_Combination_DEC/I");
	//TBranch* Comb_PROB_Branch = new_tree->Branch("Bs_Combination_PROB", &Bs_Combination_PROB, "Bs_Combination_PROB/D");	
	TBranch* OS_Comb_DEC_Branch_Run2 = new_tree_Run2->Branch("OS_Combination_DEC", &Bs_OS_Combination_DEC, "OS_Combination_DEC/I");
	TBranch* OS_Comb_PROB_Branch_Run2 = new_tree_Run2->Branch("OS_Combination_PROB", &Bs_OS_Combination_PROB, "OS_Combination_PROB/D");

	int numEvents_Run1 = tagging_tree_Run1->GetEntries();
	if(new_tree_Run1->GetEntries() != numEvents_Run1){
		cout << "ERROR::Number of events in Run1 data tree and tagging tree are different" << endl <<  new_tree_Run1->GetEntries() << endl << numEvents_Run1 << endl;
		throw "ERROR";
	}

	for(int i=0; i< numEvents_Run1; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Run1 << endl;
		tagging_tree_Run1->GetEntry(i);
		
		b_OS_Muon_DEC_Run1->GetEntry(i);
		b_OS_Muon_PROB_Run1->GetEntry(i);
		if(Bs_OS_Muon_DEC != OS_Muon_DEC || Bs_OS_Muon_ETA != OS_Muon_ETA){
			cout << "ERROR::Inconsistent events" << endl;
			throw "ERROR";
		}

		Bs_OS_Combination_DEC = OS_Combination_DEC;
		Bs_OS_Combination_PROB = OS_Combination_ETA;
		//Bs_Combination_DEC = Combination_DEC;
		//Bs_Combination_PROB = Combination_ETA;
	
		OS_Comb_DEC_Branch_Run1->Fill();
		OS_Comb_PROB_Branch_Run1->Fill();
		//Comb_DEC_Branch->Fill();
		//Comb_PROB_Branch->Fill();
	}
	
	int numEvents_Run2 = tagging_tree_Run2->GetEntries();
	if(new_tree_Run2->GetEntries() != numEvents_Run2){
		 cout << "ERROR: Number of events in Run2 data tree and tagging tree are different" << endl <<  new_tree_Run2->GetEntries() << endl << numEvents_Run2 << endl;
		throw "ERROR";
	}

	for(int i=0; i< numEvents_Run2; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Run2 << endl;
		tagging_tree_Run2->GetEntry(i);
	
		b_OS_Muon_DEC_Run2->GetEntry(i);
		b_OS_Muon_PROB_Run2->GetEntry(i);
		if(Bs_OS_Muon_DEC != OS_Muon_DEC || Bs_OS_Muon_ETA != OS_Muon_ETA){
			cout << "ERROR::Inconsistent events" << endl;
			throw "ERROR";
		}

		Bs_OS_Combination_DEC = OS_Combination_DEC;
		Bs_OS_Combination_PROB = OS_Combination_ETA;
		//Bs_Combination_DEC = Combination_DEC;
		//Bs_Combination_PROB = Combination_ETA;
	
		OS_Comb_DEC_Branch_Run2->Fill();
		OS_Comb_PROB_Branch_Run2->Fill();
		//Comb_DEC_Branch->Fill();
		//Comb_PROB_Branch->Fill();
	}

  	TList *list = new TList; 
  	list->Add(new_tree_Run1); 
  	list->Add(new_tree_Run2); 
  
  	TTree *new_tree = TTree::MergeTrees(list); 
  	new_tree->Write();

	output->Close();	
	old_file->Close();
	tagging_file_Run1->Close();
	tagging_file_Run2->Close();

	cout << "Wrote output file " << outputName << endl;
}

int main(int argc, char** argv){
  	mergeTrees("/auto/data/dargent/BsDsKpipi/BDT/Data/signal_18_newBDT","OS_combo_Run1_signal","OS_combo_Run2_signal");
 	mergeTrees("/auto/data/dargent/BsDsKpipi/BDT/Data/norm_18_newBDT","OS_combo_Run1","OS_combo_Run2");
//  	mergeTrees("/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/TD-MINT2/src/Users/dargent/MassFits/test8","OS_combo_Run1_signal","OS_combo_Run2_signal");
	return 0;
}
