// Phaspace acceptance studies
// author: Philippe d'Argent, Matthieu Kecke
#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventList.h"
#include "Mint/DiskResidentEventList.h"
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
#include "Mint/NamedParameter.h"
#include <TROOT.h>
#include <sstream>
#include "TProfile.h"
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"

using namespace std;
using namespace MINT;

double cosThetaAngle(TLorentzVector p0,TLorentzVector p1,TLorentzVector p2,TLorentzVector p3){
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(TLorentzVector p0,TLorentzVector p1,TLorentzVector p2,TLorentzVector p3){
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}

double cosThetaAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}

// void prepareDataFilesForBDT(){
// 
//     NamedParameter<string> InputGenMC("InputGenMC", (std::string) "/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_13266008_PG.root");
// 
//     TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
//     tree_gen->Add(((string)InputGenMC).c_str());
//     double K_gen[5]; 
//     double pip_gen[5]; 
//     double pim_gen[5]; 
//     double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];
//     
//     tree_gen->SetBranchAddress("Kplus_TRUEP_X",&K_gen[0]);
//     tree_gen->SetBranchAddress("Kplus_TRUEP_Y",&K_gen[1]);
//     tree_gen->SetBranchAddress("Kplus_TRUEP_Z",&K_gen[2]); 
//     tree_gen->SetBranchAddress("Kplus_TRUEP_E",&K_gen[3]); 
//     tree_gen->SetBranchAddress("Kplus_TRUEPT",&K_gen[4]); 
// 	
//     tree_gen->SetBranchAddress("piplus_TRUEP_X",&pip_gen[0]);
//     tree_gen->SetBranchAddress("piplus_TRUEP_Y",&pip_gen[1]);
//     tree_gen->SetBranchAddress("piplus_TRUEP_Z",&pip_gen[2]); 
//     tree_gen->SetBranchAddress("piplus_TRUEP_E",&pip_gen[3]); 
//     tree_gen->SetBranchAddress("piplus_TRUEPT",&pip_gen[4]); 
// 
//     tree_gen->SetBranchAddress("piminus_TRUEP_X",&pim_gen[0]);
//     tree_gen->SetBranchAddress("piminus_TRUEP_Y",&pim_gen[1]);
//     tree_gen->SetBranchAddress("piminus_TRUEP_Z",&pim_gen[2]); 
//     tree_gen->SetBranchAddress("piminus_TRUEP_E",&pim_gen[3]); 
//     tree_gen->SetBranchAddress("piminus_TRUEPT",&pim_gen[4]); 
// 	
//     tree_gen->SetBranchAddress("Kplus0_TRUEP_X",&Ds_Kp_gen[0]);
//     tree_gen->SetBranchAddress("Kplus0_TRUEP_Y",&Ds_Kp_gen[1]);
//     tree_gen->SetBranchAddress("Kplus0_TRUEP_Z",&Ds_Kp_gen[2]); 
//     tree_gen->SetBranchAddress("Kplus0_TRUEP_E",&Ds_Kp_gen[3]); 
//     tree_gen->SetBranchAddress("Kplus0_TRUEPT",&Ds_Kp_gen[4]); 
//     
//     tree_gen->SetBranchAddress("Kminus_TRUEP_X",&Ds_Km_gen[0]);
//     tree_gen->SetBranchAddress("Kminus_TRUEP_Y",&Ds_Km_gen[1]);
//     tree_gen->SetBranchAddress("Kminus_TRUEP_Z",&Ds_Km_gen[2]); 
//     tree_gen->SetBranchAddress("Kminus_TRUEP_E",&Ds_Km_gen[3]); 
//     tree_gen->SetBranchAddress("Kminus_TRUEPT",&Ds_Km_gen[4]); 
// 
//     tree_gen->SetBranchAddress("piminus0_TRUEP_X",&Ds_pim_gen[0]);
//     tree_gen->SetBranchAddress("piminus0_TRUEP_Y",&Ds_pim_gen[1]);
//     tree_gen->SetBranchAddress("piminus0_TRUEP_Z",&Ds_pim_gen[2]); 
//     tree_gen->SetBranchAddress("piminus0_TRUEP_E",&Ds_pim_gen[3]); 
//     tree_gen->SetBranchAddress("piminus0_TRUEPT",&Ds_pim_gen[4]); 
// 
//     TFile* output_gen = new TFile("GenMC.root","RECREATE");
//     TTree* summary_tree_gen = tree_gen->CloneTree(0);
// 
//     double s_Kpipi,s_Kpi,s_pipi,s_Dspi,s_Dspipi,s_Kpip,s_Dspip,s_DsK,s_DsKpi;
//     double cos_theta_Kpi,cos_theta_Dspi,phi_Kpi_Dspi;
//     double cos_theta_pipi,cos_theta_DsK,phi_pipi_DsK;
//     double w;
// 
//     summary_tree_gen->Branch("s_Kpipi", &s_Kpipi, "s_Kpipi/D");
//     summary_tree_gen->Branch("s_Kpi", &s_Kpi, "s_Kpi/D");
//     summary_tree_gen->Branch("s_pipi", &s_pipi, "s_pipi/D");
//     summary_tree_gen->Branch("s_Dspi", &s_Dspi, "s_Dspi/D");
//     summary_tree_gen->Branch("s_Dspipi", &s_Dspipi, "s_Dspipi/D");
//     summary_tree_gen->Branch("s_Dspip", &s_Dspip, "s_Dspip/D");
//     summary_tree_gen->Branch("s_DsK", &s_DsK, "s_DsK/D");
//     summary_tree_gen->Branch("s_DsKpi", &s_DsKpi, "s_DsKpi/D");
//     summary_tree_gen->Branch("s_Kpip", &s_Kpip, "s_Kpip/D");
// 
//     summary_tree_gen->Branch("cos_theta_Kpi", &cos_theta_Kpi, "cos_theta_Kpi/D");
//     summary_tree_gen->Branch("cos_theta_Dspi", &cos_theta_Dspi, "cos_theta_Dspi/D");
//     summary_tree_gen->Branch("phi_Kpi_Dspi", &phi_Kpi_Dspi, "phi_Kpi_Dspi/D");
// 
//     summary_tree_gen->Branch("cos_theta_pipi", &cos_theta_pipi, "cos_theta_pipi/D");
//     summary_tree_gen->Branch("cos_theta_DsK", &cos_theta_DsK, "cos_theta_DsK/D");
//     summary_tree_gen->Branch("phi_pipi_DsK", &phi_pipi_DsK, "phi_pipi_DsK/D");
//     summary_tree_gen->Branch("weight", &w, "weight/D");
// 
//     for(int i=0; i< tree_gen->GetEntries(); i++)
//     {	
//         tree_gen->GetEntry(i);
// 
//         TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
//         TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
//         TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
//         TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
//         TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
//         TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
//         TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
//         TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
// 
// 	s_Kpipi = (K_p + pip_p + pim_p).M()*MeV;
// 	s_Kpi = (K_p + pim_p).M()*MeV;
// 	s_pipi = (pip_p + pim_p).M()*MeV;
// 	s_Dspi = (D_p + pip_p ).M()*MeV;
// 	s_Dspipi = (D_p + pip_p + pim_p).M()*MeV;
// 
// 	s_Dspip = (D_p + pip_p).M()*MeV;
// 	s_DsKpi = (D_p + K_p + pim_p).M()*MeV;
// 	s_DsK = (D_p + K_p).M()*MeV;
// 	s_Kpip = (K_p + pip_p).M()*MeV;
// 
// 	cos_theta_Kpi = cosThetaAngle(K_p,pim_p,D_p,pip_p);
// 	cos_theta_Dspi = cosThetaAngle(D_p,pip_p,K_p,pim_p);
// 	phi_Kpi_Dspi = acoplanarityAngle(K_p,pim_p,D_p,pip_p);
// 
// 	cos_theta_pipi = cosThetaAngle(pip_p,pim_p,D_p,K_p);
// 	cos_theta_DsK = cosThetaAngle(D_p,K_p,pip_p,pim_p);
// 	phi_pipi_DsK = acoplanarityAngle(pip_p,pim_p,D_p,K_p);
// 
// 	if(s_Kpipi > 1950. || s_Kpi > 1200. || s_pipi > 1200.) continue;
// 
// 	w = 1;
// 	summary_tree_gen->Fill();
//     }
// 
//     summary_tree_gen->Write();
//     output_gen->Close();
// 
//     NamedParameter<string> InputSelectedMC("InputSelectedMC", (std::string) "/auto/data/dargent/BsDsKpipi/Final/MC/signal_scaled.root");
//     TChain* tree=new TChain("DecayTree");
//     tree->Add(((string)InputSelectedMC).c_str());
// 
//     double K[5]; 
//     double pip[5]; 
//     double pim[5]; 
//     double Ds[5];
//     int cat;
//     tree->SetBranchAddress("Bs_BKGCAT",&cat);
//     tree->SetBranchAddress("K_plus_TRUEP_X",&K[0]);
//     tree->SetBranchAddress("K_plus_TRUEP_Y",&K[1]);
//     tree->SetBranchAddress("K_plus_TRUEP_Z",&K[2]); 
//     tree->SetBranchAddress("K_plus_TRUEP_E",&K[3]); 
//     tree->SetBranchAddress("K_plus_TRUEPT",&K[4]); 
// 	
//     tree->SetBranchAddress("pi_plus_TRUEP_X",&pip[0]);
//     tree->SetBranchAddress("pi_plus_TRUEP_Y",&pip[1]);
//     tree->SetBranchAddress("pi_plus_TRUEP_Z",&pip[2]); 
//     tree->SetBranchAddress("pi_plus_TRUEP_E",&pip[3]); 
//     tree->SetBranchAddress("pi_plus_TRUEPT",&pip[4]); 
// 
//     tree->SetBranchAddress("pi_minus_TRUEP_X",&pim[0]);
//     tree->SetBranchAddress("pi_minus_TRUEP_Y",&pim[1]);
//     tree->SetBranchAddress("pi_minus_TRUEP_Z",&pim[2]); 
//     tree->SetBranchAddress("pi_minus_TRUEP_E",&pim[3]); 
//     tree->SetBranchAddress("pi_minus_TRUEPT",&pim[4]); 
// 	
//     tree->SetBranchAddress("Ds_TRUEP_X",&Ds[0]);
//     tree->SetBranchAddress("Ds_TRUEP_Y",&Ds[1]);
//     tree->SetBranchAddress("Ds_TRUEP_Z",&Ds[2]); 
//     tree->SetBranchAddress("Ds_TRUEP_E",&Ds[3]); 
//     tree->SetBranchAddress("Ds_TRUEPT",&Ds[4]); 
// 
//     tree->SetBranchAddress("weight",&w); 
// 
//     TFile* output = new TFile("SelMC.root","RECREATE");
//     TTree* summary_tree = tree->CloneTree(0);
// 
//     summary_tree->Branch("s_Kpipi", &s_Kpipi, "s_Kpipi/D");
//     summary_tree->Branch("s_Kpi", &s_Kpi, "s_Kpi/D");
//     summary_tree->Branch("s_pipi", &s_pipi, "s_pipi/D");
//     summary_tree->Branch("s_Dspi", &s_Dspi, "s_Dspi/D");
//     summary_tree->Branch("s_Dspipi", &s_Dspipi, "s_Dspipi/D");
//     summary_tree->Branch("s_Dspip", &s_Dspip, "s_Dspip/D");
//     summary_tree->Branch("s_DsK", &s_DsK, "s_DsK/D");
//     summary_tree->Branch("s_DsKpi", &s_DsKpi, "s_DsKpi/D");
//     summary_tree->Branch("s_Kpip", &s_Kpip, "s_Kpip/D");
//     summary_tree->Branch("cos_theta_Kpi", &cos_theta_Kpi, "cos_theta_Kpi/D");
//     summary_tree->Branch("cos_theta_Dspi", &cos_theta_Dspi, "cos_theta_Dspi/D");
//     summary_tree->Branch("phi_Kpi_Dspi", &phi_Kpi_Dspi, "phi_Kpi_Dspi/D");
//     summary_tree->Branch("cos_theta_pipi", &cos_theta_pipi, "cos_theta_pipi/D");
//     summary_tree->Branch("cos_theta_DsK", &cos_theta_DsK, "cos_theta_DsK/D");
//     summary_tree->Branch("phi_pipi_DsK", &phi_pipi_DsK, "phi_pipi_DsK/D");
// 
//     summary_tree->Branch("weight", &w, "weight/D");
// 
//     int numEvents = tree->GetEntries();
//     for(int i=0; i< numEvents; i++)
//     {
//         tree->GetEntry(i);
//         
//         // Lorentz vectors: P=(Px,Py,Pz,E)
//         TLorentzVector K_p(K[0],K[1],K[2],K[3]);
//         TLorentzVector pip_p(pip[0],pip[1],pip[2],pip[3]);
//         TLorentzVector pim_p(pim[0],pim[1],pim[2],pim[3]);
//         TLorentzVector D_p(Ds[0],Ds[1],Ds[2],Ds[3]);
//         TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
// 	
// 	s_Kpipi = (K_p + pip_p + pim_p).M()*MeV;
// 	s_Kpi = (K_p + pim_p).M()*MeV;
// 	s_pipi = (pip_p + pim_p).M()*MeV;
// 	s_Dspi = (D_p + pip_p ).M()*MeV;
// 	s_Dspipi = (D_p + pip_p + pim_p).M()*MeV;
// 
// 	s_Dspip = (D_p + pip_p).M()*MeV;
// 	s_DsKpi = (D_p + K_p + pim_p).M()*MeV;
// 	s_DsK = (D_p + K_p).M()*MeV;
// 	s_Kpip = (K_p + pip_p).M()*MeV;
// 
// 	cos_theta_Kpi = cosThetaAngle(K_p,pim_p,D_p,pip_p);
// 	cos_theta_Dspi = cosThetaAngle(D_p,pip_p,K_p,pim_p);
// 	phi_Kpi_Dspi = acoplanarityAngle(K_p,pim_p,D_p,pip_p);
// 
// 	cos_theta_pipi = cosThetaAngle(pip_p,pim_p,D_p,K_p);
// 	cos_theta_DsK = cosThetaAngle(D_p,K_p,pip_p,pim_p);
// 	phi_pipi_DsK = acoplanarityAngle(pip_p,pim_p,D_p,K_p);
// 
// 	if(cat != 20) continue;
// 	if(s_Kpipi > 1950. || s_Kpi > 1200. || s_pipi > 1200.) continue;
// 
// 	summary_tree->Fill();
//     }
// 
//     summary_tree->Write();
//     output->Close();
// 
//     DalitzEventPattern pdg(531, -431, 321, 211, -211);
//     DalitzEventList eventList;
//     TFile *file_MINT =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
//     TTree* tree_MINT=dynamic_cast<TTree*>(file_MINT->Get("DalitzEventList"));
// 
//     TFile* output_MINT = new TFile("MintMC.root","RECREATE");
//     TTree* summary_tree_MINT = tree_MINT->CloneTree();
// 
//     eventList.fromNtuple(tree_MINT,1);
//     cout << " I've got " << eventList.size() << " events." << endl;
//     file_MINT->Close();
// 
//     TBranch* b_Kpipi = summary_tree_MINT->Branch("s_Kpipi", &s_Kpipi, "s_Kpipi/D");
//     TBranch* b_Kpi = summary_tree_MINT->Branch("s_Kpi", &s_Kpi, "s_Kpi/D");
//     TBranch* b_pipi =summary_tree_MINT->Branch("s_pipi", &s_pipi, "s_pipi/D");
//     TBranch* b_Dspi =summary_tree_MINT->Branch("s_Dspi", &s_Dspi, "s_Dspi/D");
//     TBranch* b_Dspipi =summary_tree_MINT->Branch("s_Dspipi", &s_Dspipi, "s_Dspipi/D");
//     TBranch* b_Dspip =summary_tree_MINT->Branch("s_Dspip", &s_Dspip, "s_Dspip/D");
//     TBranch* b_DsK =summary_tree_MINT->Branch("s_DsK", &s_DsK, "s_DsK/D");
//     TBranch* b_cos_theta_Kpi =summary_tree_MINT->Branch("cos_theta_Kpi", &cos_theta_Kpi, "cos_theta_Kpi/D");
//     TBranch* b_cos_theta_Dspi =summary_tree_MINT->Branch("cos_theta_Dspi", &cos_theta_Dspi, "cos_theta_Dspi/D");
//     TBranch* b_phi_Kpi_Dspi =summary_tree_MINT->Branch("phi_Kpi_Dspi", &phi_Kpi_Dspi, "phi_Kpi_Dspi/D");
// 
//     TBranch* b_w = summary_tree_MINT->Branch("weight", &w, "weight/D");
// 
// 
//     vector<int> s234;
//     s234.push_back(2);
//     s234.push_back(3);
//     s234.push_back(4);
//     vector<int> s134;
//     s134.push_back(1);
//     s134.push_back(3);
//     s134.push_back(4);
// 
//     for(int i=0; i<eventList.size();i++){
// 
// 	s_Kpipi = sqrt(eventList[i].sij(s234))*MeV;
// 	s_Kpi = sqrt(eventList[i].s(2,4))*MeV;
// 	s_pipi = sqrt(eventList[i].s(3,4))*MeV;
// 	s_Dspi = sqrt(eventList[i].s(1,3))*MeV;
// 	s_Dspipi = sqrt(eventList[i].sij(s134))*MeV;
// 	s_Dspip = sqrt(eventList[i].s(1,4))*MeV;
// 	s_DsK = sqrt(eventList[i].s(1,2))*MeV;
// 	cos_theta_Kpi = cosThetaAngle(eventList[i],2,4,1,3);
// 	cos_theta_Dspi = cosThetaAngle(eventList[i],1,3,2,4);
// 	phi_Kpi_Dspi = acoplanarityAngle(eventList[i],2,4,1,3);
// 	w = 1;
// 
// 	summary_tree_MINT->GetEntry(i);
// 
// 	b_w->Fill();
// 	b_Kpipi->Fill();
//         b_Kpi->Fill();
//         b_pipi->Fill();
//         b_Dspi->Fill();
//         b_Dspipi->Fill();
//         b_Dspip->Fill();
//         b_DsK->Fill();
//         b_cos_theta_Kpi->Fill();
//         b_cos_theta_Dspi->Fill();
//         b_phi_Kpi_Dspi->Fill();
//     }
//     summary_tree_MINT->Write();
//     output_MINT->Close();
// }
// 
//     NamedParameter<string> InputGenMC("InputGenMC", (std::string) "/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_13266008_PG.root");
// 
//     TChain* tree_gen=new TChain("DecayTree");
//     tree->Add(((string)InputFileName).c_str());
// 		tree->SetBranchStatus("*",0);
// 		tree->SetBranchStatus("N_Bs_sw",1);
// 		tree->SetBranchStatus("year",1);
// 		tree->SetBranchStatus("*DEC",1);
// 		tree->SetBranchStatus("*PROB",1);
// 		tree->SetBranchStatus("*OS*",1);
// 		tree->SetBranchStatus("*TAU*",1);
// 		tree->SetBranchStatus("*ID*",1);
// 		tree->SetBranchStatus("weight",1);
// 		tree->SetBranchStatus("Bs_DTF_MM",1);
//     		tree->SetBranchStatus("BsDTF_*P*",1);
//     double K[4];
//     double pip[4];
//     double pim[4];
//     double Ds_Kp[4],Ds_Km[4],Ds_pim[4],Ds[4];
//     
// 		tree->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
// 		tree->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
// 		tree->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]);
// 		tree->SetBranchAddress("BsDTF_Kplus_PE",&K[3]);
// 		
// 		tree->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
// 		tree->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
// 		tree->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]);
// 		tree->SetBranchAddress("BsDTF_piplus_PE",&pip[3]);
// 		
// 		tree->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
// 		tree->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
// 		tree->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]);
// 		tree->SetBranchAddress("BsDTF_piminus_PE",&pim[3]);
// 
// 
//     TFile* output_gen = new TFile("GenMC.root","RECREATE");
//     TTree* summary_tree_gen = tree_gen->CloneTree(0);
// 
//     double s_Kpipi,s_Kpi,s_pipi,s_Dspi,s_Dspipi,s_Kpip,s_Dspip,s_DsK,s_DsKpi;
//     double cos_theta_Kpi,cos_theta_Dspi,phi_Kpi_Dspi;
//     double cos_theta_pipi,cos_theta_DsK,phi_pipi_DsK;
//     double w;
// 
//     summary_tree_gen->Branch("s_Kpipi", &s_Kpipi, "s_Kpipi/D");
//     summary_tree_gen->Branch("s_Kpi", &s_Kpi, "s_Kpi/D");
//     summary_tree_gen->Branch("s_pipi", &s_pipi, "s_pipi/D");
//     summary_tree_gen->Branch("s_Dspi", &s_Dspi, "s_Dspi/D");
//     summary_tree_gen->Branch("s_Dspipi", &s_Dspipi, "s_Dspipi/D");
//     summary_tree_gen->Branch("s_Dspip", &s_Dspip, "s_Dspip/D");
//     summary_tree_gen->Branch("s_DsK", &s_DsK, "s_DsK/D");
//     summary_tree_gen->Branch("s_DsKpi", &s_DsKpi, "s_DsKpi/D");
//     summary_tree_gen->Branch("s_Kpip", &s_Kpip, "s_Kpip/D");
// 
//     summary_tree_gen->Branch("cos_theta_Kpi", &cos_theta_Kpi, "cos_theta_Kpi/D");
//     summary_tree_gen->Branch("cos_theta_Dspi", &cos_theta_Dspi, "cos_theta_Dspi/D");
//     summary_tree_gen->Branch("phi_Kpi_Dspi", &phi_Kpi_Dspi, "phi_Kpi_Dspi/D");
// 
//     summary_tree_gen->Branch("cos_theta_pipi", &cos_theta_pipi, "cos_theta_pipi/D");
//     summary_tree_gen->Branch("cos_theta_DsK", &cos_theta_DsK, "cos_theta_DsK/D");
//     summary_tree_gen->Branch("phi_pipi_DsK", &phi_pipi_DsK, "phi_pipi_DsK/D");
//     summary_tree_gen->Branch("weight", &w, "weight/D");
// 
//     for(int i=0; i< tree_gen->GetEntries(); i++)
//     {	
//         tree_gen->GetEntry(i);
// 
//         TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
//         TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
//         TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
//         TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
//         TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
//         TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
//         TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
//         TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
// 
// 	s_Kpipi = (K_p + pip_p + pim_p).M()*MeV;
// 	s_Kpi = (K_p + pim_p).M()*MeV;
// 	s_pipi = (pip_p + pim_p).M()*MeV;
// 	s_Dspi = (D_p + pip_p ).M()*MeV;
// 	s_Dspipi = (D_p + pip_p + pim_p).M()*MeV;
// 
// 	s_Dspip = (D_p + pip_p).M()*MeV;
// 	s_DsKpi = (D_p + K_p + pim_p).M()*MeV;
// 	s_DsK = (D_p + K_p).M()*MeV;
// 	s_Kpip = (K_p + pip_p).M()*MeV;
// 
// 	cos_theta_Kpi = cosThetaAngle(K_p,pim_p,D_p,pip_p);
// 	cos_theta_Dspi = cosThetaAngle(D_p,pip_p,K_p,pim_p);
// 	phi_Kpi_Dspi = acoplanarityAngle(K_p,pim_p,D_p,pip_p);
// 
// 	cos_theta_pipi = cosThetaAngle(pip_p,pim_p,D_p,K_p);
// 	cos_theta_DsK = cosThetaAngle(D_p,K_p,pip_p,pim_p);
// 	phi_pipi_DsK = acoplanarityAngle(pip_p,pim_p,D_p,K_p);
// 
// 	if(s_Kpipi > 1950. || s_Kpi > 1200. || s_pipi > 1200.) continue;
// 
// 	w = 1;
// 	summary_tree_gen->Fill();
//     }
// 
//     summary_tree_gen->Write();
//     output_gen->Close();
// }

void reweightGen(){
	
    double s_Kpipi,s_Kpi,s_pipi,s_Dspi,s_Dspipi,s_Kpip,s_Dspip,s_DsK,s_DsKpi;
    double cos_theta_Kpi,cos_theta_Dspi,phi_Kpi_Dspi;
    double BDTG,w;

    TChain* tree_gen=new TChain("MCDecayTree");
    tree_gen->Add("GenMC_BDT.root");

    tree_gen->SetBranchAddress("s_Kpipi",&s_Kpipi);
    tree_gen->SetBranchAddress("s_Kpi",&s_Kpi);
    tree_gen->SetBranchAddress("s_pipi",&s_pipi);
    tree_gen->SetBranchAddress("s_Dspipi",&s_Dspipi);
    tree_gen->SetBranchAddress("s_Dspi",&s_Dspi);
    tree_gen->SetBranchAddress( "cos_theta_Kpi", &cos_theta_Kpi );
    tree_gen->SetBranchAddress( "cos_theta_Dspi", &cos_theta_Dspi );
    tree_gen->SetBranchAddress( "phi_Kpi_Dspi", &phi_Kpi_Dspi );
    tree_gen->SetBranchAddress("BDTG",&BDTG);

    TChain* tree=new TChain("DecayTree");
    tree->Add("SelMC_BDT.root");

    tree->SetBranchAddress("s_Kpipi",&s_Kpipi);
    tree->SetBranchAddress("s_Kpi",&s_Kpi);
    tree->SetBranchAddress("s_pipi",&s_pipi);
    tree->SetBranchAddress("s_Dspipi",&s_Dspipi);
    tree->SetBranchAddress("s_Dspi",&s_Dspi);
    tree->SetBranchAddress( "cos_theta_Kpi", &cos_theta_Kpi );
    tree->SetBranchAddress( "cos_theta_Dspi", &cos_theta_Dspi );
    tree->SetBranchAddress( "phi_Kpi_Dspi", &phi_Kpi_Dspi );
    tree->SetBranchAddress("BDTG",&BDTG);
    tree->SetBranchAddress("weight",&w);

    TChain* tree_data=new TChain("DecayTree");
    tree_data->Add("signal_tagged_BDT.root");
    tree_data->SetBranchAddress("m_Kpipi",&s_Kpipi);
    tree_data->SetBranchAddress("m_Kpi",&s_Kpi);
    tree_data->SetBranchAddress("m_pipi",&s_pipi);
    tree_data->SetBranchAddress("m_Dspipi",&s_Dspipi);
    tree_data->SetBranchAddress("m_Dspi",&s_Dspi);
    tree_data->SetBranchAddress("BDTG_eff",&BDTG);
    tree_data->SetBranchAddress("N_Bs_sw",&w);


    TChain* tree_MINT=new TChain("DalitzEventList");
    tree_MINT->Add("MintMC_BDT.root");
    tree_MINT->SetBranchAddress("s_Kpipi",&s_Kpipi);
    tree_MINT->SetBranchAddress("s_Kpi",&s_Kpi);
    tree_MINT->SetBranchAddress("s_pipi",&s_pipi);
    tree_MINT->SetBranchAddress("s_Dspipi",&s_Dspipi);
    tree_MINT->SetBranchAddress("s_Dspi",&s_Dspi);
    tree_MINT->SetBranchAddress( "cos_theta_Kpi", &cos_theta_Kpi );
    tree_MINT->SetBranchAddress( "cos_theta_Dspi", &cos_theta_Dspi );
    tree_MINT->SetBranchAddress( "phi_Kpi_Dspi", &phi_Kpi_Dspi );
    tree_MINT->SetBranchAddress("BDTG",&BDTG);


    TH1D* h_Kpipi= new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (MeV);Events (norm.) ",40,1100,2000);
    TH1D* h_Kpipi_gen= new TH1D("","",40,1100,2000);
    TH1D* h_Kpipi_gen_rw= new TH1D("","",40,1100,2000);
    TH1D* h_Kpipi_MINT= new TH1D("","",40,1100,2000);
    TH1D* h_Kpipi_MINT_rw= new TH1D("","",40,1100,2000);

    TH1D* h_Kpi= new TH1D("",";#left[m(K^{+} #pi^{-})#right] (MeV);Events (norm.) ",40,600,1200);
    TH1D* h_Kpi_gen= new TH1D("","",40,600,1200);
    TH1D* h_Kpi_gen_rw= new TH1D("","",40,600,1200);
    TH1D* h_Kpi_MINT= new TH1D("","",40,600,1200);
    TH1D* h_Kpi_MINT_rw= new TH1D("","",40,600,1200);

    TH1D* h_pipi= new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (MeV);Events (norm.) ",40,200,1200);
    TH1D* h_pipi_gen= new TH1D("","",40,200,1200);
    TH1D* h_pipi_gen_rw= new TH1D("","",40,200,1200);
    TH1D* h_pipi_MINT= new TH1D("","",40,200,1200);
    TH1D* h_pipi_MINT_rw= new TH1D("","",40,200,1200);

    TH1D* h_Dspi= new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (MeV);Events (norm.) ",40,1900,5000);
    TH1D* h_Dspi_gen= new TH1D("","",40,1900,5000);
    TH1D* h_Dspi_gen_rw= new TH1D("","",40,1900,5000);
    TH1D* h_Dspi_MINT= new TH1D("","",40,1900,5000);
    TH1D* h_Dspi_MINT_rw= new TH1D("","",40,1900,5000);

    TH1D* h_Dspipi= new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (MeV);Events (norm.) ",40,2400,5100);
    TH1D* h_Dspipi_gen= new TH1D("","",40,2400,5100);
    TH1D* h_Dspipi_gen_rw= new TH1D("","",40,2400,5100);
    TH1D* h_Dspipi_MINT= new TH1D("","",40,2400,5100);
    TH1D* h_Dspipi_MINT_rw= new TH1D("","",40,2400,5100);

    TH1D* h_cosTheta_Kpi= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
    TH1D* h_cosTheta_Kpi_gen= new TH1D("","",40,-1,1);
    TH1D* h_cosTheta_Kpi_gen_rw= new TH1D("","",40,-1,1);
    TH1D* h_cosTheta_Kpi_MINT= new TH1D("","",40,-1,1);
    TH1D* h_cosTheta_Kpi_MINT_rw= new TH1D("","",40,-1,1);

    TH1D* h_cosTheta_Dspi= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    TH1D* h_cosTheta_Dspi_gen= new TH1D("","",40,0,1);
    TH1D* h_cosTheta_Dspi_gen_rw= new TH1D("","",40,0,1);
    TH1D* h_cosTheta_Dspi_MINT= new TH1D("","",40,0,1);
    TH1D* h_cosTheta_Dspi_MINT_rw= new TH1D("","",40,0,1);

    TH1D* h_phi_Kpi_Dspi= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);
    TH1D* h_phi_Kpi_Dspi_gen= new TH1D("","",40,-3.141,3.141);
    TH1D* h_phi_Kpi_Dspi_gen_rw= new TH1D("","",40,-3.141,3.141);
    TH1D* h_phi_Kpi_Dspi_MINT= new TH1D("","",40,-3.141,3.141);
    TH1D* h_phi_Kpi_Dspi_MINT_rw= new TH1D("","",40,-3.141,3.141);

    TH1D* bdt_gen= new TH1D("","; BTDG; Efficiency (norm.)",50,-0.65,0.25);
    TH1D* bdt_sel= new TH1D("","; BTDG; Events (norm.)",50,-0.65,0.25);
    TH1D* bdt_MINT= new TH1D("","",50,-0.65,0.25);

    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);
	bdt_gen->Fill(BDTG);
    }

    for(int i=0; i< tree->GetEntries(); i++)
    {	
        tree->GetEntry(i);
	bdt_sel->Fill(BDTG,w);
	h_Kpipi->Fill(s_Kpipi,w);
	h_Dspipi->Fill(s_Dspipi,w);
	h_Kpi->Fill(s_Kpi,w);
	h_pipi->Fill(s_pipi,w);
	h_Dspi->Fill(s_Dspi,w);
	h_cosTheta_Kpi->Fill(cos_theta_Kpi,w);
	h_cosTheta_Dspi->Fill(cos_theta_Dspi,w);
	h_phi_Kpi_Dspi->Fill(phi_Kpi_Dspi,w);
    }
    	
    TH1D *h_weight = (TH1D*)bdt_gen->Clone();
    h_weight->Divide(bdt_sel,bdt_gen);
    h_weight->Scale(1./h_weight->Integral());
    TCanvas* c = new TCanvas();
    h_weight->Draw();
    c->Print("eff.eps");
    c->Print("eff.pdf");

    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);
	h_Kpipi_gen->Fill(s_Kpipi);
	h_Dspipi_gen->Fill(s_Dspipi);
	h_Kpi_gen->Fill(s_Kpi);
	h_pipi_gen->Fill(s_pipi);
	h_Dspi_gen->Fill(s_Dspi);
	h_cosTheta_Kpi_gen->Fill(cos_theta_Kpi);
	h_cosTheta_Dspi_gen->Fill(cos_theta_Dspi);
	h_phi_Kpi_Dspi_gen->Fill(phi_Kpi_Dspi);

	double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
	h_Kpipi_gen_rw->Fill(s_Kpipi,weight);
	h_Dspipi_gen_rw->Fill(s_Dspipi,weight);
	h_Kpi_gen_rw->Fill(s_Kpi,weight);
	h_pipi_gen_rw->Fill(s_pipi,weight);
	h_Dspi_gen_rw->Fill(s_Dspi,weight);
	h_cosTheta_Kpi_gen_rw->Fill(cos_theta_Kpi,weight);
	h_cosTheta_Dspi_gen_rw->Fill(cos_theta_Dspi,weight);
	h_phi_Kpi_Dspi_gen_rw->Fill(phi_Kpi_Dspi,weight);
    }

   TFile *hFile = new TFile("signal_tagged_BDT_effCorr.root","RECREATE");
//    tree_data->SetBranchStatus("N_Bs_sw",0);
   TTree* new_tree = tree_data->CloneTree(0);
//    tree_data->SetBranchStatus("N_Bs_sw",1);

   double w_new;
   new_tree->Branch("N_Bs_sw_new",&w_new, "w_new/D");
 

   double sumw=0.;
   double sumw2=0.;
   for(int i=0; i< tree_data->GetEntries(); i++)
    {	
        tree_data->GetEntry(i);
	double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
	sumw = 1./weight;
	sumw2 = pow(1./weight,2);
    }

   for(int i=0; i< tree_data->GetEntries(); i++)
    {	
        tree_data->GetEntry(i);
	double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	

// 	w_new = (weight == 0.00000) ? w : w*1./weight;
	w_new = w * 1./weight * sumw/sumw2;
	new_tree->Fill();
    }


    DalitzEventPattern pdg(531, -431, 321, 211, -211);

    TFile *file_integ =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
    TTree* tree_integ=dynamic_cast<TTree*>(file_integ->Get("DalitzEventList"));

    DalitzEventList eventList,eventList_rw,eventList_rw_CP;
    eventList.fromNtuple(tree_integ,1);
    eventList[0].print();

    if(tree_MINT->GetEntries() != tree_integ->GetEntries()) throw "ERROR";

    for(int i=0; i< tree_MINT->GetEntries(); i++)
    {	
        tree_MINT->GetEntry(i);

	bdt_MINT->Fill(BDTG);

	h_Kpipi_MINT->Fill(s_Kpipi);
	h_Dspipi_MINT->Fill(s_Dspipi);
	h_Kpi_MINT->Fill(s_Kpi);
	h_pipi_MINT->Fill(s_pipi);
	h_Dspi_MINT->Fill(s_Dspi);
	h_cosTheta_Kpi_MINT->Fill(cos_theta_Kpi);
	h_cosTheta_Dspi_MINT->Fill(cos_theta_Dspi);
	h_phi_Kpi_Dspi_MINT->Fill(phi_Kpi_Dspi);

	double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
	h_Kpipi_MINT_rw->Fill(s_Kpipi,weight);
	h_Dspipi_MINT_rw->Fill(s_Dspipi,weight);
	h_Kpi_MINT_rw->Fill(s_Kpi,weight);
	h_pipi_MINT_rw->Fill(s_pipi,weight);
	h_Dspi_MINT_rw->Fill(s_Dspi,weight);
	h_cosTheta_Kpi_MINT_rw->Fill(cos_theta_Kpi,weight);
	h_cosTheta_Dspi_MINT_rw->Fill(cos_theta_Dspi,weight);
	h_phi_Kpi_Dspi_MINT_rw->Fill(phi_Kpi_Dspi,weight);

	DalitzEvent evt(eventList[i]);
	evt.setWeight(evt.getWeight()*weight);
	eventList_rw.Add(evt);
        evt.CP_conjugateYourself();
        eventList_rw_CP.Add(evt);
    }
    
    eventList_rw.saveAsNtuple("SignalIntegrationEvents_AccBDT_new.root");
    eventList_rw_CP.saveAsNtuple("SignalIntegrationEvents_AccBDT_CP_new.root");

    bdt_sel->SetLineColor(kBlue);
    bdt_sel->DrawNormalized("hist",1);
    bdt_gen->SetLineColor(kRed);
    bdt_gen->SetMarkerColor(kRed);
    bdt_gen->DrawNormalized("same",1);
    //bdt_MINT->DrawNormalized("histsame",1);
    c->Print("BDTG.eps");
    c->Print("BDTG.pdf");

    h_Kpipi->DrawNormalized("",1);
    h_Kpipi_gen->SetLineColor(kRed);
    h_Kpipi_gen->DrawNormalized("histsame",1);
    h_Kpipi_gen_rw->SetLineColor(kBlue);
    h_Kpipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_Kpipi.eps");
    c->Print("h_Kpipi.pdf");

    TH1D* h_Kpipi_eff = (TH1D*)h_Kpipi->Clone();
    h_Kpipi_eff->Divide(h_Kpipi,h_Kpipi_gen);
    h_Kpipi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Kpipi_eff->DrawNormalized("e",1);
    TH1D* h_Kpipi_MINT_eff = (TH1D*)h_Kpipi_MINT->Clone();
    h_Kpipi_MINT_eff->Divide(h_Kpipi_MINT_rw,h_Kpipi_MINT);
    h_Kpipi_MINT_eff->SetLineColor(kBlue);
    h_Kpipi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_Kpipi.eps");
    c->Print("eff_Kpipi.pdf");

    h_Kpi->DrawNormalized("",1);
    h_Kpi_gen->SetLineColor(kRed);
    h_Kpi_gen->DrawNormalized("histsame",1);
    h_Kpi_gen_rw->SetLineColor(kBlue);
    h_Kpi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_Kpi.eps");
    c->Print("h_Kpi.pdf");

    TH1D* h_Kpi_eff = (TH1D*)h_Kpi->Clone();
    h_Kpi_eff->Divide(h_Kpi,h_Kpi_gen);
    h_Kpi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Kpi_eff->DrawNormalized("e",1);
    TH1D* h_Kpi_MINT_eff = (TH1D*)h_Kpi_MINT->Clone();
    h_Kpi_MINT_eff->Divide(h_Kpi_MINT_rw,h_Kpi_MINT);
    h_Kpi_MINT_eff->SetLineColor(kBlue);
    h_Kpi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_Kpi.eps");
    c->Print("eff_Kpi.pdf");

    h_pipi->DrawNormalized("",1);
    h_pipi_gen->SetLineColor(kRed);
    h_pipi_gen->DrawNormalized("histsame",1);
    h_pipi_gen_rw->SetLineColor(kBlue);
    h_pipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_pipi.eps");
    c->Print("h_pipi.pdf");

    TH1D* h_pipi_eff = (TH1D*)h_pipi->Clone();
    h_pipi_eff->Divide(h_pipi,h_pipi_gen);
    h_pipi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_pipi_eff->DrawNormalized("e",1);
    TH1D* h_pipi_MINT_eff = (TH1D*)h_pipi_MINT->Clone();
    h_pipi_MINT_eff->Divide(h_pipi_MINT_rw,h_pipi_MINT);
    h_pipi_MINT_eff->SetLineColor(kBlue);
    h_pipi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_pipi.eps");
    c->Print("eff_pipi.pdf");

    h_Dspi->DrawNormalized("",1);
    h_Dspi_gen->SetLineColor(kRed);
    h_Dspi_gen->DrawNormalized("histsame",1);
    h_Dspi_gen_rw->SetLineColor(kBlue);
    h_Dspi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_Dspi.eps");
    c->Print("h_Dspi.pdf");

    TH1D* h_Dspi_eff = (TH1D*)h_Dspi->Clone();
    h_Dspi_eff->Divide(h_Dspi,h_Dspi_gen);
    h_Dspi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Dspi_eff->DrawNormalized("e",1);
    TH1D* h_Dspi_MINT_eff = (TH1D*)h_Dspi_MINT->Clone();
    h_Dspi_MINT_eff->Divide(h_Dspi_MINT_rw,h_Dspi_MINT);
    h_Dspi_MINT_eff->SetLineColor(kBlue);
    h_Dspi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_Dspi.eps");
    c->Print("eff_Dspi.pdf");

    h_Dspipi->DrawNormalized("",1);
    h_Dspipi_gen->SetLineColor(kRed);
    h_Dspipi_gen->DrawNormalized("histsame",1);
    h_Dspipi_gen_rw->SetLineColor(kBlue);
    h_Dspipi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_Dspipi.eps");
    c->Print("h_Dspipi.pdf");

    TH1D* h_Dspipi_eff = (TH1D*)h_Dspipi->Clone();
    h_Dspipi_eff->Divide(h_Dspipi,h_Dspipi_gen);
    h_Dspipi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Dspipi_eff->DrawNormalized("e",1);
    TH1D* h_Dspipi_MINT_eff = (TH1D*)h_Dspipi_MINT->Clone();
    h_Dspipi_MINT_eff->Divide(h_Dspipi_MINT_rw,h_Dspipi_MINT);
    h_Dspipi_MINT_eff->SetLineColor(kBlue);
    h_Dspipi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_Dspipi.eps");
    c->Print("eff_Dspipi.pdf");

    h_cosTheta_Kpi->DrawNormalized("",1);
    h_cosTheta_Kpi_gen->SetLineColor(kRed);
    h_cosTheta_Kpi_gen->DrawNormalized("histsame",1);
    h_cosTheta_Kpi_gen_rw->SetLineColor(kBlue);
    h_cosTheta_Kpi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_cosTheta_Kpi.eps");
    c->Print("h_cosTheta_Kpi.pdf");

    TH1D* h_cosTheta_Kpi_eff = (TH1D*)h_cosTheta_Kpi->Clone();
    h_cosTheta_Kpi_eff->Divide(h_cosTheta_Kpi,h_cosTheta_Kpi_gen);
    h_cosTheta_Kpi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_cosTheta_Kpi_eff->DrawNormalized("e",1);
    TH1D* h_cosTheta_Kpi_MINT_eff = (TH1D*)h_cosTheta_Kpi_MINT->Clone();
    h_cosTheta_Kpi_MINT_eff->Divide(h_cosTheta_Kpi_MINT_rw,h_cosTheta_Kpi_MINT);
    h_cosTheta_Kpi_MINT_eff->SetLineColor(kBlue);
    h_cosTheta_Kpi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_cosTheta_Kpi.eps");
    c->Print("eff_cosTheta_Kpi.pdf");

    h_cosTheta_Dspi->DrawNormalized("",1);
    h_cosTheta_Dspi_gen->SetLineColor(kRed);
    h_cosTheta_Dspi_gen->DrawNormalized("histsame",1);
    h_cosTheta_Dspi_gen_rw->SetLineColor(kBlue);
    h_cosTheta_Dspi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_cosTheta_Dspi.eps");
    c->Print("h_cosTheta_Dspi.pdf");

    TH1D* h_cosTheta_Dspi_eff = (TH1D*)h_cosTheta_Dspi->Clone();
    h_cosTheta_Dspi_eff->Divide(h_cosTheta_Dspi,h_cosTheta_Dspi_gen);
    h_cosTheta_Dspi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_cosTheta_Dspi_eff->DrawNormalized("e",1);
    TH1D* h_cosTheta_Dspi_MINT_eff = (TH1D*)h_cosTheta_Dspi_MINT->Clone();
    h_cosTheta_Dspi_MINT_eff->Divide(h_cosTheta_Dspi_MINT_rw,h_cosTheta_Dspi_MINT);
    h_cosTheta_Dspi_MINT_eff->SetLineColor(kBlue);
    h_cosTheta_Dspi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_cosTheta_Dspi.eps");
    c->Print("eff_cosTheta_Dspi.pdf");

    h_phi_Kpi_Dspi->DrawNormalized("",1);
    h_phi_Kpi_Dspi_gen->SetLineColor(kRed);
    h_phi_Kpi_Dspi_gen->DrawNormalized("histsame",1);
    h_phi_Kpi_Dspi_gen_rw->SetLineColor(kBlue);
    h_phi_Kpi_Dspi_gen_rw->DrawNormalized("histsame",1);
    c->Print("h_phi_Kpi_Dspi.eps");
    c->Print("h_phi_Kpi_Dspi.pdf");

    TH1D* h_phi_Kpi_Dspi_eff = (TH1D*)h_phi_Kpi_Dspi->Clone();
    h_phi_Kpi_Dspi_eff->Divide(h_phi_Kpi_Dspi,h_phi_Kpi_Dspi_gen);
    h_phi_Kpi_Dspi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_phi_Kpi_Dspi_eff->DrawNormalized("e",1);
    TH1D* h_phi_Kpi_Dspi_MINT_eff = (TH1D*)h_phi_Kpi_Dspi_MINT->Clone();
    h_phi_Kpi_Dspi_MINT_eff->Divide(h_phi_Kpi_Dspi_MINT_rw,h_phi_Kpi_Dspi_MINT);
    h_phi_Kpi_Dspi_MINT_eff->SetLineColor(kBlue);
    h_phi_Kpi_Dspi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_phi_Kpi_Dspi.eps");
    c->Print("eff_phi_Kpi_Dspi.pdf");
}

int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

//     prepareFilesForBDT();

    reweightGen();


    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}