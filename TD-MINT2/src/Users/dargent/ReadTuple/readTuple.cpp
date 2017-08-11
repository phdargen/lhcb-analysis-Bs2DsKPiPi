// Fits Bs mass distribution and calculates sweights
// author: Philippe d'Argent, Matthieu Kecke
#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventList.h"
#include "TTree.h"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Mint/DalitzEventPattern.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "Mint/PlotSet.h"
#include "TChain.h"

using namespace std;
using namespace MINT;

int makeMINTtuple(){
    
    string outputDir = "/auto/data/dargent/BsDsKpipi/MINT/";
    
    bool dbThis=false;
    bool addSweight = false;
    bool bkg = false;
    
    int N=-1;
    if(dbThis) cout << "read ntuple" << endl;
	
    DalitzEventPattern pdg (531, -431, 321, 211, -211);
    if (dbThis) cout << "event pattern : " << pdg << endl;
	
    DalitzEventList eventList; 

    // Read the momenta from ntuple
    TChain* tree=new TChain("DecayTree");
    tree->Add("/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
   
    if (dbThis) cout << "Read the file" << endl;	
    // 4-momenta
    double K[4]; 
    double pip[4]; 
    double pim[4]; 
    double Ds_Kp[4],Ds_Km[4],Ds_pim[4];
    double sweight;
    double mB;
    
    tree->SetBranchAddress("Bs_DTF_MM",&mB);
    tree->SetBranchAddress("N_Bs_sw",&sweight);
    
    tree->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
    tree->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
    tree->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]); 
    tree->SetBranchAddress("BsDTF_Kplus_PE",&K[3]); 
	
    tree->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
    tree->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
    tree->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]); 
    tree->SetBranchAddress("BsDTF_piplus_PE",&pip[3]); 
	
    tree->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
    tree->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
    tree->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]); 
    tree->SetBranchAddress("BsDTF_piminus_PE",&pim[3]); 
	
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PX",&Ds_Kp[0]);
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PY",&Ds_Kp[1]);
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PZ",&Ds_Kp[2]); 
    tree->SetBranchAddress("BsDTF_Ds_Kplus_PE",&Ds_Kp[3]); 
    
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PX",&Ds_Km[0]);
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PY",&Ds_Km[1]);
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PZ",&Ds_Km[2]); 
    tree->SetBranchAddress("BsDTF_Ds_Kminus_PE",&Ds_Km[3]); 

    tree->SetBranchAddress("BsDTF_Ds_piminus_PX",&Ds_pim[0]);
    tree->SetBranchAddress("BsDTF_Ds_piminus_PY",&Ds_pim[1]);
    tree->SetBranchAddress("BsDTF_Ds_piminus_PZ",&Ds_pim[2]); 
    tree->SetBranchAddress("BsDTF_Ds_piminus_PE",&Ds_pim[3]); 
    

    int numEvents = tree->GetEntries();
    int numSelected =0;

    //loop over tree and fill eventList
    for(int i=0; i< numEvents; i++)
    {
	if(dbThis)cout << " getting " << i << " th entry" << endl;	
	tree->GetEntry(i);
        
        // Lorentz vectors: P=(Px,Py,Pz,E)
        TLorentzVector K_p(K[0],K[1],K[2],K[3]);
        TLorentzVector pip_p(pip[0],pip[1],pip[2],pip[3]);
	TLorentzVector pim_p(pim[0],pim[1],pim[2],pim[3]);
        TLorentzVector D_Kp_p(Ds_Kp[0],Ds_Kp[1],Ds_Kp[2],Ds_Kp[3]);
        TLorentzVector D_Km_p(Ds_Km[0],Ds_Km[1],Ds_Km[2],Ds_Km[3]);
        TLorentzVector D_pim_p(Ds_pim[0],Ds_pim[1],Ds_pim[2],Ds_pim[3]);
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 
	if(dbThis) cout <<	"the 4-momenta have been defined" << endl;
        //cuts
         if(!addSweight){
	    if(bkg){
	    	    if(mB< 5500)continue;
	    }
	    else {
	            if(mB<5360)continue;
        	    if(mB>5390)continue;		
	    }
        }
	// define the order of the vectors in the vectorOfvectors
        // include the 'MeV' to get the correct units, need to include CLHEPSystemOfUnits.h
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);

	if(dbThis) cout << "make event" << endl;
		
	DalitzEvent evt(pdg, vectorOfvectors);
        if(addSweight)evt.setWeight(sweight);
	//if(evt.phaseSpace()==0) cout << evt << endl;
	
	if(evt.s(1,2)<0 || evt.s(2,3)<0 || evt.s(3,4)<0 || evt.t(4,0)<0 || evt.t(1,0)< 0) 
	cout << "negative mass ?" << endl; 
        if(evt.s(1,2)< pdg.sijMin(1,2) || evt.s(1,2)> pdg.sijMax(1,2) ) continue;
        if(evt.s(2,3)< pdg.sijMin(2,3) || evt.s(2,3)> pdg.sijMax(2,3) )continue; 
        if(evt.s(3,4)< pdg.sijMin(3,4) || evt.s(3,4)> pdg.sijMax(3,4) ) continue;
        if(evt.t(0,4)< pdg.sijMin(1,2,3) || evt.t(4,0)> pdg.sijMax(1,2,3) ) continue;  
        if(evt.t(1,0)< pdg.sijMin(2,3,4) || evt.t(1,0)> pdg.sijMax(2,3,4) ) continue; 	
	if(dbThis) cout << "s12 =  " << ( evt.p(1) + evt.p(2) ).Mag2() << endl;
	if(dbThis) cout << "adding event " << evt << endl;
	
	eventList.Add(evt); // this fills the event list		
	if(dbThis) cout << " added event" << endl;
		
        numSelected++;
        if(numSelected==N)break;
    }
    
    string output = outputDir + "MINT_data";
    if(!addSweight){
	if(bkg) output+="_bkg";
	else output+="_3sigma";
    }
    if(N != -1)output += "_small";
    output+=".root";
    
    eventList.save(output.c_str());
    DalitzHistoSet datH = eventList.weightedHistoSet();
    datH.draw("data_","eps"); 
   
    cout << numSelected << " / " << numEvents << " events selected" << endl;
    cout << "Created File: " << output.c_str() << endl;    

    return 0;
}


int plot(){
   	 TH1::SetDefaultSumw2();

        DalitzEventPattern pat (531, -431, 321, 211, -211);
    	DalitzEventList eventListPhsp,eventList;
   	DalitzEventList eventList2,eventListPhsp2;
    	eventListPhsp.generatePhaseSpaceEvents(200000,pat);

        TFile *_InputFile =  TFile::Open("/auto/data/dargent/BsDsKpipi/MINT/MINT_data.root");
        TTree* in_tree, *in_tree2;
        in_tree=dynamic_cast<TTree*>(_InputFile->Get("DalitzEventList"));
	eventList.fromNtuple(in_tree,1);
        TFile *_InputFile2 =  TFile::Open("/auto/data/dargent/BsDsKpipi/MINT/MINT_data_bkg.root");
        in_tree2=dynamic_cast<TTree*>(_InputFile2->Get("DalitzEventList"));
	eventList2.fromNtuple(in_tree2,1);

	cout << "Have " << eventList.size() << " events " << endl;

	vector<int> s234;
        s234.push_back(2);
        s234.push_back(3);
        s234.push_back(4);

        TH1D* s_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.8,2);
        TH1D* s_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV);Events (norm.) ", 50,0.4,1.3);
        TH1D* s_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.,1.3);

	double sw = 0;
        double tw = 0;
        for (int i=0; i<eventList.size(); i++) {
	    if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 2 && sqrt(eventList[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventList[i].s(3,4)/(GeV*GeV) < 1.2)) sw += eventList[i].getWeight() ;     
     	    //if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) > 1.95) continue; 
     	    //if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 1.4) continue; 

	    s_Kpipi->Fill(sqrt(eventList[i].sij(s234)/(GeV*GeV)),eventList[i].getWeight());
            s_Kpi->Fill(sqrt(eventList[i].s(2,4)/(GeV*GeV)),eventList[i].getWeight());
            s_pipi->Fill(sqrt(eventList[i].s(3,4)/(GeV*GeV)),eventList[i].getWeight());
            tw += eventList[i].getWeight() ; 
	}    

	cout << "sw = " << sw << endl;
        cout << "tw = " << tw << endl;
        cout << "eff = " << sw/tw << endl;

        TH1D* s_Kpipi2 = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.8,2);
        TH1D* s_Kpi2 = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV);Events (norm.) ", 50,0.4,1.3);
        TH1D* s_pipi2 = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.,1.3);

        for (int i=0; i<eventList2.size(); i++) {
            s_Kpipi2->Fill(sqrt(eventList2[i].sij(s234)/(GeV*GeV)));
            s_Kpi2->Fill(sqrt(eventList2[i].s(2,4)/(GeV*GeV)));
            s_pipi2->Fill(sqrt(eventList2[i].s(3,4)/(GeV*GeV)));
        }    

        TH1D* s_Kpipi3 = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.8,2);
        TH1D* s_Kpi3 = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV);Events (norm.) ", 50,0.4,1.3);
        TH1D* s_pipi3 = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV);Events (norm.) ",50,0.,1.3);

        for (int i=0; i<eventListPhsp.size(); i++) {
            s_Kpipi3->Fill(sqrt(eventListPhsp[i].sij(s234)/(GeV*GeV)));
            s_Kpi3->Fill(sqrt(eventListPhsp[i].s(2,4)/(GeV*GeV)));
            s_pipi3->Fill(sqrt(eventListPhsp[i].s(3,4)/(GeV*GeV)));
	    if(sqrt(eventListPhsp[i].sij(s234)/(GeV*GeV)) < 2 && sqrt(eventListPhsp[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventListPhsp[i].s(3,4)/(GeV*GeV) < 1.2)) eventListPhsp2.Add(eventListPhsp[i]); 
        }    

	cout << s_Kpipi->Integral() << endl;

	cout << s_Kpipi2->Integral() << endl;


	TCanvas* c = new TCanvas();
	s_Kpipi->SetLineColor(kBlue);
	s_Kpipi->DrawNormalized("e1",1);
	s_Kpipi2->SetLineColor(kRed);
	s_Kpipi2->DrawNormalized("e1same",1);
	s_Kpipi3->SetLineColor(kBlack);
        s_Kpipi3->DrawNormalized("histcsame",1);
        c->Print("m_Kpipi.eps");
	//c->Print("/home/he/dargent/public_html/BsDsKpipi/m_Kpipi.pdf");

	s_Kpi->SetLineColor(kBlue);
	s_Kpi->DrawNormalized("e1");
	s_Kpi2->SetLineColor(kRed);
	s_Kpi2->DrawNormalized("e1same");
	s_Kpi3->SetLineColor(kBlack);
        s_Kpi3->DrawNormalized("histcsame");
        c->Print("m_Kpi.eps");
	//c->Print("/home/he/dargent/public_html/BsDsKpipi/m_Kpi.pdf");

	s_pipi->SetLineColor(kBlue);
	s_pipi->DrawNormalized("e1");
	s_pipi2->SetLineColor(kRed);
	s_pipi2->DrawNormalized("e1same");
	s_pipi3->SetLineColor(kBlack);
        s_pipi3->DrawNormalized("histcsame");
        c->Print("m_pipi.eps");
	//c->Print("/home/he/dargent/public_html/BsDsKpipi/m_pipi.pdf");

	cout << "phsps events with m_Kpipi < 2 = " << (double)eventListPhsp2.size()/eventListPhsp.size() ;

        //eventListPhsp.saveAsNtuple("phsp.root");
	//eventListPhsp2.saveAsNtuple("phsp_cut.root");
	

return 0;
}


int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    makeMINTtuple();
    plot();
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
