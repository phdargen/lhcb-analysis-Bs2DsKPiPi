void mergeTrees(string old_fileName , string tagging_fileName, string cut){

string outputName = old_fileName;

string oldFileName = old_fileName.append(".root");
string taggingFileName = tagging_fileName.append(".root");

TFile* old_file= new TFile(oldFileName.c_str());
TTree* old_tree = (TTree*) old_file->Get("DecayTree");

TFile* tagging_file= new TFile(taggingFileName.c_str());
TTree* tagging_tree = (TTree*) tagging_file->Get("TaggingTree");



Int_t OS_Combination_DEC;
Double_t OS_Combination_ETA;

Int_t Combination_DEC;
Double_t Combination_ETA;

tagging_tree -> SetBranchAddress("OS_Combination_DEC" , &OS_Combination_DEC );
tagging_tree -> SetBranchAddress("OS_Combination_ETA" , &OS_Combination_ETA );
tagging_tree -> SetBranchAddress("Combination_DEC" , &Combination_DEC );
tagging_tree -> SetBranchAddress("Combination_ETA" , &Combination_ETA );



outputName.append("_wTagging.root"); 

TFile* output=new TFile(outputName.c_str(),"RECREATE");
TTree* new_tree = old_tree->CopyTree(cut.c_str());

Int_t Bs_OS_Combination_DEC;
Double_t Bs_OS_Combination_PROB;

Int_t Bs_Combination_DEC;
Double_t Bs_Combination_PROB;


TBranch* OS_Comb_DEC_Branch = new_tree->Branch("Bs_OS_Combination_DEC", &Bs_OS_Combination_DEC, "Bs_OS_Combination_DEC/I");
TBranch* OS_Comb_PROB_Branch = new_tree->Branch("Bs_OS_Combination_PROB", &Bs_OS_Combination_PROB, "Bs_OS_Combination_PROB/D");

TBranch* Comb_DEC_Branch = new_tree->Branch("Bs_Combination_DEC", &Bs_Combination_DEC, "Bs_Combination_DEC/I");
TBranch* Comb_PROB_Branch = new_tree->Branch("Bs_Combination_PROB", &Bs_Combination_PROB, "Bs_Combination_PROB/D");


int numEvents = tagging_tree->GetEntries();
for(int i=0; i< numEvents; i++)
        {
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tagging_tree->GetEntry(i);


        Bs_OS_Combination_DEC = OS_Combination_DEC;
        Bs_OS_Combination_PROB = OS_Combination_ETA;
        Bs_Combination_DEC = Combination_DEC;
        Bs_Combination_PROB = Combination_ETA;


	OS_Comb_DEC_Branch->Fill();
	OS_Comb_PROB_Branch->Fill();
	Comb_DEC_Branch->Fill();
	Comb_PROB_Branch->Fill();

}

new_tree->Write();
output->Close();

old_file->Close();
tagging_file->Close();

}
