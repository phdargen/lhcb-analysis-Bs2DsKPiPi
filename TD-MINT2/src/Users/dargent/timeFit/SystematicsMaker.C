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
    pull p_stat(paraNames,"signal_toy5/pull_*.root");
    TMatrixD* cov_stat = new TMatrixD(p_stat.getStatCov());
    cov_stat->Print();

    /// Fit bias from toys
    pull p(paraNames,"signal_toy5/pull_*.root");
    TMatrixD* cov = new TMatrixD(p.getCov());
    cov->Print();
    covs.push_back(cov);

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
    pull p_acc(paraNames,"signal_toy5/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("signal_toy5/pull_*.root","_acc"));
    cov_acc->Print();
    covs.push_back(cov_acc);

    /// Acc systematics (with cholesky)
    pull p_acc_chol(paraNames,"signal_toy5/pullAccChol_*.root");
    TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("signal_toy5/pull_*.root","_accChol",100));
    cov_acc_chol->Print();
    covs.push_back(cov_acc_chol);

    /// Total systematics table    
    ofstream SummaryFile;
    SummaryFile.open("pull_results/summary_table.tex",std::ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i <covs.size() ; i++) SummaryFile << " c " ;
    SummaryFile << " c }" << "\n";
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Fit Parameter & " ;
    for(int i =0 ; i <covs.size() ; i++)  SummaryFile << i+1 << " & " ;
    SummaryFile << " Total " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";
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
    SummaryFile << "\\end{tabular}" << "\n";

    /// Total systematics table in terms of sigma_stat   
    ofstream SummaryFile2;
    SummaryFile2.open("pull_results/summary_table2.tex",std::ofstream::trunc);

    SummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i <covs.size() ; i++) SummaryFile2 << " c " ;
    SummaryFile2 << " c }" << "\n";
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "Fit Parameter & " ;
    for(int i =0 ; i <covs.size() ; i++)  SummaryFile2 << i+1 << " & " ;
    SummaryFile2 << " Total " << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";
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
    SummaryFile2 << "\\end{tabular}" << "\n";

    
    return 0;
}
