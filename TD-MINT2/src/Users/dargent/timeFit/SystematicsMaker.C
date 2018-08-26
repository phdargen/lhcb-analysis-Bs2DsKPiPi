#include "pull.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>

using namespace std;

int main(int argc, char** argv){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    //gStyle->SetOptFit(111);
    //gStyle->UseCurrentStyle();

    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("C");
    paraNames.push_back("D");
    paraNames.push_back("D_bar");
    paraNames.push_back("S");
    paraNames.push_back("S_bar");

    vector<TMatrixD*> covs;

    /// Stat cov from toys
    pull p_stat(paraNames,"signal_toy6/pull__*.root");
    TMatrixD* cov_stat = new TMatrixD(p_stat.getStatCov());
    cov_stat->Print();

    /// Fit bias from toys
     pull p(paraNames,"signal_toy6/pull__*.root");
//     TMatrixD* cov = new TMatrixD(p.getCov());
//     cov->Print();
//     covs.push_back(cov);

    /// Systematics from data fits
    vector<TString> fileNames;
    /// Resolution systematics
    //fileNames.push_back("");
    /// ...
    //fileNames.push_back("");

    for (int i= 0; i<fileNames.size(); i++) {
        pull p(paraNames,fileNames[i]);
        TMatrixD* cov = new TMatrixD(p.getCov());
        cov->Print();
        covs.push_back(cov);
    }
    
    /// Acc systematics 
//     pull p_acc(paraNames,"signal_toy5/pullAcc_*.root");
//     TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("signal_toy5/pull_*.root","_acc"));
//     cov_acc->Print();
    //covs.push_back(cov_acc);

    /// Acc systematics (with cholesky)
//     pull p_acc_chol(paraNames,"signal_toy5/pullAccChol_*.root");
//     TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("signal_toy5/pull_*.root","_accChol",100));
//     cov_acc_chol->Print();
    //covs.push_back(cov_acc_chol);

    /// dms systematics 
    pull p_dm(paraNames,"signal_toy6/pull_dm_*.root");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("signal_toy6/pull__*.root","_dm"));
    cov_dm->Print();
    covs.push_back(cov_dm);

    /// asymmetry systematics
    pull p_production_asym_Run1(paraNames,"signal_toy6/pull_production_asym_Run1_*.root");
    TMatrixD* cov_production_asym_Run1 = new TMatrixD(p_production_asym_Run1.getDeltaCov("signal_toy6/pull__*.root","_production_asym_Run1"));
    cov_production_asym_Run1->Print();

    pull p_production_asym_Run2(paraNames,"signal_toy6/pull_production_asym_Run2_*.root");
    TMatrixD* cov_production_asym_Run2 = new TMatrixD(p_production_asym_Run2.getDeltaCov("signal_toy6/pull__*.root","_production_asym_Run2"));
    cov_production_asym_Run2->Print();

    pull p_detection_asym_Run1(paraNames,"signal_toy6/pull_detection_asym_Run1_*.root");
    TMatrixD* cov_detection_asym_Run1 = new TMatrixD(p_detection_asym_Run1.getDeltaCov("signal_toy6/pull__*.root","_detection_asym_Run1"));
    cov_detection_asym_Run1->Print();

    pull p_detection_asym_Run2(paraNames,"signal_toy6/pull_detection_asym_Run2_*.root");
    TMatrixD* cov_detection_asym_Run2 = new TMatrixD(p_detection_asym_Run2.getDeltaCov("signal_toy6/pull__*.root","_detection_asym_Run2"));
    cov_detection_asym_Run2->Print();

    TMatrixD cov_asym(*cov_production_asym_Run1);
    cov_asym +=  *cov_production_asym_Run2 ;
    cov_asym +=  *cov_detection_asym_Run1 ;
    cov_asym +=   *cov_detection_asym_Run2;

    covs.push_back(new TMatrixD(cov_asym));

    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"signal_sys_res_Run1_a/pull__1.root");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("signal/pull_1.root","_res_Run1_a"));
    cov_res_Run1_a->Print();

    pull p_res_Run1_b(paraNames,"signal_sys_res_Run1_b/pull__1.root");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("signal/pull_1.root","_res_Run1_b"));
    cov_res_Run1_b->Print();

    pull p_res_Run2_a(paraNames,"signal_sys_res_Run2_a/pull__1.root");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("signal/pull_1.root","_res_Run2_a"));
    cov_res_Run2_a->Print();

    pull p_res_Run2_b(paraNames,"signal_sys_res_Run2_b/pull__1.root");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("signal/pull_1.root","_res_Run2_b"));
    cov_res_Run2_b->Print();

    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    cov_res +=  p.combineCov_maxVal(covs_res_Run2) ;
    
    covs.push_back(new TMatrixD(cov_res));



    /// Total systematics table   
    vector<string> sysNames;
    sysNames.push_back("Fit bias");
    sysNames.push_back("Acceptance");
    sysNames.push_back("Resolution");
    sysNames.push_back("$\\Delta m_{s}$");
    sysNames.push_back("Asymmetries");
 
    ofstream SummaryFile;
    SummaryFile.open("pull_results/summary_table.tex",std::ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i <covs.size() ; i++) SummaryFile << " c " ;
    SummaryFile << " | c }" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Fit Parameter & " ;
    for(int i =0 ; i <covs.size() ; i++)  SummaryFile << sysNames[i] << " & " ;
    SummaryFile << " Total " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";
    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile << std::fixed << std::setprecision(2) << p.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile << sqrt((*covs[j])[i][i]) << " & ";  
        }
        SummaryFile << sqrt(tot) << " \\\\ " << "\n"; 
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";

    /// Total systematics table in terms of sigma_stat   
    ofstream SummaryFile2;
    SummaryFile2.open("pull_results/summary_table2.tex",std::ofstream::trunc);

    SummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i <covs.size() ; i++) SummaryFile2 << " c " ;
    SummaryFile2 << " | c }" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "Fit Parameter & " ;
    for(int i =0 ; i <covs.size() ; i++)  SummaryFile2 << sysNames[i] << " & " ;
    SummaryFile2 << " Total " << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";
    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile2 << std::fixed << std::setprecision(2) << p.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile2 << sqrt((*covs[j])[i][i])/sqrt((*cov_stat)[i][i]) << " & ";  
        }
        SummaryFile2 << sqrt(tot)/sqrt((*cov_stat)[i][i]) << " \\\\ " << "\n"; 
    }
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";

    return 0;
}
