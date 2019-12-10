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

void signal(){


    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("C");
    paraNames.push_back("D");
    paraNames.push_back("D_bar");
    paraNames.push_back("S");
    paraNames.push_back("S_bar");

    vector<TMatrixD*> covs;

    /// Data fit
    pull p_data(paraNames,"signal_test5/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();
    
    /// Stat cov from toys
    pull p_stat(paraNames,"signal_toy_ub/pull__*.root");
    TMatrixD* cov_stat = new TMatrixD(p_stat.getStatCov());
    cov_stat->Print();

    /// Fit bias from toys
    pull p(paraNames,"signal_toy_ub/pull__*.root");
    TMatrixD* cov = new TMatrixD(p.getCov());
//     cov->Print();
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov)[i][j] = (*cov)[i][j] * errs_stat[i] * errs_stat[j];
    }
    cov->Print();
    covs.push_back(cov);

    TMatrixD cor(*cov);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor[i][j] = (*cov)[i][j]/sqrt((*cov)[i][i])/sqrt((*cov)[j][j]);

    /// Acc systematics 
    pull p_acc(paraNames,"signal_toy_ub/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("signal_toy_ub/pull__*.root","_acc"));
    cov_acc->Print();

    TMatrixD cor_acc(*cov_acc);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_acc[i][j] = (*cov_acc)[i][j]/sqrt((*cov_acc)[i][i])/sqrt((*cov_acc)[j][j]);

    /// Acc systematics (with cholesky)
//     pull p_acc_chol(paraNames,"signal_toy/pullAccChol_*.root");
//     TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("signal_toy/pull__*.root","_accChol",50));
//     cov_acc_chol->Print();
//     covs.push_back(cov_acc_chol);

    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"signal_sys_res_Run1_a/pull__*.root");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("signal/pull__1.root","_res_Run1_a"));
//     cov_res_Run1_a->Print();

    pull p_res_Run1_b(paraNames,"signal_sys_res_Run1_b/pull__*.root");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("signal/pull__1.root","_res_Run1_b"));
//     cov_res_Run1_b->Print();

    pull p_res_Run2_a(paraNames,"signal_sys_res_Run2_a/pull__*.root");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("signal/pull__1.root","_res_Run2_a"));
//     cov_res_Run2_a->Print();

    pull p_res_Run2_b(paraNames,"signal_sys_res_Run2_b/pull__*.root");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("signal/pull__1.root","_res_Run2_b"));
//     cov_res_Run2_b->Print();

    pull p_res_Run2_c(paraNames,"signal_sys_res_Run2_c/pull__*.root");
    TMatrixD* cov_res_Run2_c = new TMatrixD(p_res_Run2_c.getDeltaCov("signal/pull__1.root","_res_Run2_c"));
//     cov_res_Run2_a->Print();

    pull p_res_Run2_d(paraNames,"signal_sys_res_Run2_d/pull__*.root");
    TMatrixD* cov_res_Run2_d = new TMatrixD(p_res_Run2_d.getDeltaCov("signal/pull__1.root","_res_Run2_d"));
//     cov_res_Run2_b->Print();

    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    covs_res_Run2.push_back(cov_res_Run2_c);
    covs_res_Run2.push_back(cov_res_Run2_d);

    cov_res +=  p.combineCov_maxVal(covs_res_Run2) ;

    TMatrixD cor_res(cov_res);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_res[i][j] = cov_res[i][j]/sqrt(cov_res[i][i])/sqrt(cov_res[j][j]);
    

    /// dms systematics 
    pull p_dm(paraNames,"signal_toy_ub/pull_dm_*.root");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("signal_toy_ub/pull__*.root","_dm"));
    cov_dm->Print();

    TMatrixD cor_dm(*cov_dm);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_dm[i][j] = (*cov_dm)[i][j]/sqrt((*cov_dm)[i][i])/sqrt((*cov_dm)[j][j]);
    
    /// asymmetry systematics
    pull p_production_asym_Run1(paraNames,"signal_toy_ub/pull_production_asym_Run1_*.root");
    TMatrixD* cov_production_asym_Run1 = new TMatrixD(p_production_asym_Run1.getDeltaCov("signal_toy_ub/pull__*.root","_production_asym_Run1"));
    cov_production_asym_Run1->Print();

    pull p_production_asym_Run2(paraNames,"signal_toy_ub/pull_production_asym_Run2_*.root");
    TMatrixD* cov_production_asym_Run2 = new TMatrixD(p_production_asym_Run2.getDeltaCov("signal_toy_ub/pull__*.root","_production_asym_Run2"));
    cov_production_asym_Run2->Print();

    pull p_detection_asym_Run1(paraNames,"signal_toy_ub/pull_detection_asym_Run1_*.root");
    TMatrixD* cov_detection_asym_Run1 = new TMatrixD(p_detection_asym_Run1.getDeltaCov("signal_toy_ub/pull__*.root","_detection_asym_Run1"));
    cov_detection_asym_Run1->Print();

    pull p_detection_asym_Run2(paraNames,"signal_toy_ub/pull_detection_asym_Run2_*.root");
    TMatrixD* cov_detection_asym_Run2 = new TMatrixD(p_detection_asym_Run2.getDeltaCov("signal_toy_ub/pull__*.root","_detection_asym_Run2"));
    cov_detection_asym_Run2->Print();

    pull p_tagging_asym_Run1(paraNames,"signal_toy_ub/pull_tagging_asym_Run1_*.root");
    TMatrixD* cov_tagging_asym_Run1 = new TMatrixD(p_tagging_asym_Run1.getDeltaCov("signal_toy_ub/pull__*.root","_tagging_asym_Run1"));
    cov_tagging_asym_Run1->Print();

    pull p_tagging_asym_Run2(paraNames,"signal_toy_ub/pull_tagging_asym_Run2_*.root");
    TMatrixD* cov_tagging_asym_Run2 = new TMatrixD(p_tagging_asym_Run1.getDeltaCov("signal_toy_ub/pull__*.root","_tagging_asym_Run2"));
    cov_tagging_asym_Run2->Print();


    TMatrixD cov_asym(*cov_production_asym_Run1);
    cov_asym +=  *cov_production_asym_Run2 ;
    cov_asym +=  *cov_detection_asym_Run1 ;
    cov_asym +=   *cov_detection_asym_Run2;
    cov_asym +=   *cov_tagging_asym_Run1;
    cov_asym +=   *cov_tagging_asym_Run2;


    TMatrixD cor_asym(cov_asym);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_asym[i][j] = cov_asym[i][j]/sqrt(cov_asym[i][i])/sqrt(cov_asym[j][j]);

    /// bkg systematics 
    pull p_bkg_1(paraNames,"signal_sys_bkg1/pull__*.root");
    TMatrixD* cov_bkg_1 = new TMatrixD(p_bkg_1.getDeltaCov("signal/pull__1.root","bkg_1"));
    vector<double> vals_bkg_1 = p_bkg_1.getVals();

    pull p_bkg_2(paraNames,"signal_sys_bkg2/pull__*.root");
    TMatrixD* cov_bkg_2 = new TMatrixD(p_bkg_2.getDeltaCov("signal/pull__1.root","bkg_2"));
    //cov_bkg_2->Print();
    vector<double> vals_bkg_2 = p_bkg_2.getVals();

    pull p_bkg_3(paraNames,"signal_sys_bkg3/pull__*.root");
    TMatrixD* cov_bkg_3 = new TMatrixD(p_bkg_3.getDeltaCov("signal/pull__1.root","bkg_3"));
    //cov_bkg_3->Print();
    vector<double> vals_bkg_3 = p_bkg_3.getVals();

    pull p_bkg_4(paraNames,"signal_sys_bkg4/pull__*.root");
    TMatrixD* cov_bkg_4 = new TMatrixD(p_bkg_4.getDeltaCov("signal/pull__1.root","bkg_4"));
    //cov_bkg_4->Print();
    vector<double> vals_bkg_4 = p_bkg_4.getVals();

    pull p_bkg_5(paraNames,"signal_sys_bkg5/pull__*.root");
    TMatrixD* cov_bkg_5 = new TMatrixD(p_bkg_5.getDeltaCov("signal/pull__1.root","bkg_5"));
    //cov_bkg_5->Print();
    vector<double> vals_bkg_5 = p_bkg_5.getVals();

    pull p_bkg_6(paraNames,"signal_sys_bkg6/pull__*.root");
    TMatrixD* cov_bkg_6 = new TMatrixD(p_bkg_6.getDeltaCov("signal/pull__1.root","bkg_6"));
    //cov_bkg_6->Print();
    vector<double> vals_bkg_6 = p_bkg_6.getVals();

    pull p_bkg_7(paraNames,"signal_sys_bkg7/pull__*.root");
    TMatrixD* cov_bkg_7 = new TMatrixD(p_bkg_7.getDeltaCov("signal/pull__1.root","bkg_7"));
    vector<double> vals_bkg_7 = p_bkg_7.getVals();

    pull p_bkg_8(paraNames,"signal_sys_bkg8/pull__*.root");
    TMatrixD* cov_bkg_8 = new TMatrixD(p_bkg_8.getDeltaCov("signal/pull__1.root","bkg_8"));
    vector<double> vals_bkg_8 = p_bkg_8.getVals();

//     pull p_bkg_9(paraNames,"signal_sys_bkg9/pull__*.root");
//     TMatrixD* cov_bkg_9 = new TMatrixD(p_bkg_9.getDeltaCov("signal/pull__1.root","bkg_9"));
//     vector<double> vals_bkg_9 = p_bkg_9.getVals();

    pull p_bkg_10(paraNames,"signal_sys_bkg10/pull__*.root");
    TMatrixD* cov_bkg_10 = new TMatrixD(p_bkg_10.getDeltaCov("signal/pull__1.root","bkg_10"));
    vector<double> vals_bkg_10 = p_bkg_10.getVals();


    // Take sample variance as systematic 
    vector< vector <double> > vec_vals_bkg;
    vec_vals_bkg.push_back(vals_bkg_1);
    vec_vals_bkg.push_back(vals_bkg_2);
    vec_vals_bkg.push_back(vals_bkg_3);
    vec_vals_bkg.push_back(vals_bkg_4);
    vec_vals_bkg.push_back(vals_bkg_5);
    vec_vals_bkg.push_back(vals_bkg_6);
    vec_vals_bkg.push_back(vals_bkg_7);
    vec_vals_bkg.push_back(vals_bkg_8);
//     vec_vals_bkg.push_back(vals_bkg_9);
    vec_vals_bkg.push_back(vals_bkg_10);

    TMatrixD cov_bkg(p.sampleVariance(vec_vals_bkg));
    cov_bkg.Print();

    TMatrixD cor_bkg(cov_bkg);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_bkg[i][j] = cov_bkg[i][j]/sqrt(cov_bkg[i][i])/sqrt(cov_bkg[j][j]);

    /// m,t correlations
//     pull p_corr(paraNames,"signal_toy_bkg7/pull__*.root");
//     TMatrixD* cov_corr = new TMatrixD(p_corr.getDeltaCov("signal_toy_bkg8/pull__*.root","_corr"));

    pull p_corr1(paraNames,"signal_toy_bkg7/pull__*.root");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"signal_toy_bkg8/pull__*.root");
    TMatrixD cov_corr2 = p_corr2.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
    }

    /// Total sys corr
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(new TMatrixD(cov_asym));
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(cov_dm);
    covs.push_back(cov_corr);

    TMatrixD cov_sys_tot(paraNames.size(),paraNames.size());
    TMatrixD cor_sys_tot(paraNames.size(),paraNames.size());

    for(int j =0 ; j <covs.size() ; j++) cov_sys_tot += *covs[j];

    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_sys_tot[i][j] = cov_sys_tot[i][j]/sqrt(cov_sys_tot[i][i])/sqrt(cov_sys_tot[j][j]);

    cout << "sys corr bias " << endl;
    cor.Print();

    cout << "sys corr asym " << endl;
    cor_asym.Print();

    cout << "sys corr bkg " << endl;
    cor_bkg.Print();

    cout << "sys corr dm " << endl;
    cor_dm.Print();

    cout << "sys corr res " << endl;
    cor_res.Print();

    cout << "sys corr acc " << endl;
    cor_acc.Print();

    cout << "sys corr tot " << endl;
    cor_sys_tot.Print();


    /// Total systematics table   
    vector<string> sysNames;
    sysNames.push_back("Fit bias");
    sysNames.push_back("Background");
    sysNames.push_back("Acceptance");
    sysNames.push_back("Resolution");
    sysNames.push_back("Asymmetries");
    sysNames.push_back("$\\Delta m_{s}$");
    sysNames.push_back("Correlations");
 
    ofstream SummaryFile;
    //SummaryFile.open("pull_results/sys_summary_table.tex",std::ofstream::trunc);
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/timeFit/signal/sys_summary_table.tex",std::ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i <covs.size() ; i++) SummaryFile << " c " ;
    SummaryFile << " | c }" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Fit Parameter & " ;
    for(int i =0 ; i <covs.size() ; i++)  SummaryFile << sysNames[i] << " & " ;
    SummaryFile << " Total " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";

    vector<double> errs_sys;    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile << std::fixed << std::setprecision(2) << p.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile << sqrt((*covs[j])[i][i]) << " & ";  
        }
        SummaryFile << sqrt(tot) << " \\\\ " << "\n";
	errs_sys.push_back(sqrt(tot)); 
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";

    /// Total systematics table in terms of sigma_stat   
    ofstream SummaryFile2;
    //SummaryFile2.open("pull_results/sys_summary_table2.tex",std::ofstream::trunc);
    SummaryFile2.open("../../../../../TD-AnaNote/latex/tables/timeFit/signal/sys_summary_table2.tex",std::ofstream::trunc);

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
            SummaryFile2 << sqrt((*covs[j])[i][i])/errs_stat[i] << " & ";  
        }
        SummaryFile2 << sqrt(tot)/errs_stat[i] << " \\\\ " << "\n"; 
    }
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";



    /// Result table   
    ofstream SummaryFile3;
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/timeFit/signal/result_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{c r } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "Fit Parameter & \\multicolumn{1}{c}{Value} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile3 << p.latexName(paraNames[i])  << " & " ;
        //SummaryFile3 << std::fixed << std::setprecision(2) << "x.xx" << " $\\pm$ " << errs_stat[i]  << " $\\pm$ " << errs_sys[i];
        SummaryFile3 << std::fixed << std::setprecision(2) << vals[i] << " $\\pm$ " << errs_stat[i]  << " $\\pm$ " << errs_sys[i];
        SummaryFile3  << " \\\\ " << "\n"; 
    }
    
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";
}


void norm() {

    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("p0_os_Run1");
    paraNames.push_back("p1_os_Run1");
    paraNames.push_back("delta_p0_os_Run1");
    paraNames.push_back("delta_p1_os_Run1");
    paraNames.push_back("tageff_os_Run1");
    paraNames.push_back("tageff_asym_os_Run1");

    paraNames.push_back("p0_ss_Run1");
    paraNames.push_back("p1_ss_Run1");
    paraNames.push_back("delta_p0_ss_Run1");
    paraNames.push_back("delta_p1_ss_Run1");
    paraNames.push_back("tageff_ss_Run1");
    paraNames.push_back("tageff_asym_ss_Run1");

    paraNames.push_back("p0_os_Run2");
    paraNames.push_back("p1_os_Run2");
    paraNames.push_back("delta_p0_os_Run2");
    paraNames.push_back("delta_p1_os_Run2");
    paraNames.push_back("tageff_os_Run2");
    paraNames.push_back("tageff_asym_os_Run2");

    paraNames.push_back("p0_ss_Run2");
    paraNames.push_back("p1_ss_Run2");
    paraNames.push_back("delta_p0_ss_Run2");
    paraNames.push_back("delta_p1_ss_Run2");
    paraNames.push_back("tageff_ss_Run2");
    paraNames.push_back("tageff_asym_ss_Run2");
    paraNames.push_back("production_asym_Run2");

    paraNames.push_back("dm");


    vector<TMatrixD*> covs;

    /// Data fit
    pull p_data(paraNames,"norm_taggingCalib/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();
    
    /// Fit bias from toys
    pull p(paraNames,"norm_toy/pull__*.root");
    TMatrixD* cov = new TMatrixD(p.getCov());
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov)[i][j] = (*cov)[i][j] * errs_stat[i] * errs_stat[j];
    }
    covs.push_back(cov);
   
    /// Acc systematics (with cholesky)
//     pull p_acc_chol(paraNames,"norm_toy/pullAccChol_*.root");
//     TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("norm_toy/pull__*.root","_accChol",50));
//     covs.push_back(cov_acc_chol);

    pull p_acc(paraNames,"norm_toy/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("norm_toy/pull__*.root","_acc"));

    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"norm_sys_res_Run1_a/pull__*.root");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run1_a"));

    pull p_res_Run1_b(paraNames,"norm_sys_res_Run1_b/pull__*.root");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run1_b"));

    pull p_res_Run2_a(paraNames,"norm_sys_res_Run2_a/pull__*.root");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_a"));

    pull p_res_Run2_b(paraNames,"norm_sys_res_Run2_b/pull__*.root");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_b"));

    pull p_res_Run2_c(paraNames,"norm_sys_res_Run2_c/pull__*.root");
    TMatrixD* cov_res_Run2_c = new TMatrixD(p_res_Run2_c.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_c"));

    pull p_res_Run2_d(paraNames,"norm_sys_res_Run2_d/pull__*.root");
    TMatrixD* cov_res_Run2_d = new TMatrixD(p_res_Run2_d.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_d"));

    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p_res_Run1_a.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    covs_res_Run2.push_back(cov_res_Run2_c);
    covs_res_Run2.push_back(cov_res_Run2_d);
    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2) ;
    
    /// asymmetry systematics
    pull p_production_asym_Run1(paraNames,"norm_toy/pull_production_asym_Run1_*.root");
    TMatrixD* cov_production_asym_Run1 = new TMatrixD(p_production_asym_Run1.getDeltaCov("norm_toy/pull__*.root","_production_asym_Run1"));

    pull p_detection_asym_Run1(paraNames,"norm_toy/pull_detection_asym_Run1_*.root");
    TMatrixD* cov_detection_asym_Run1 = new TMatrixD(p_detection_asym_Run1.getDeltaCov("norm_toy/pull__*.root","_detection_asym_Run1"));

    pull p_detection_asym_Run2(paraNames,"norm_toy/pull_detection_asym_Run2_*.root");
    TMatrixD* cov_detection_asym_Run2 = new TMatrixD(p_detection_asym_Run2.getDeltaCov("norm_toy/pull__*.root","_detection_asym_Run2"));

    TMatrixD cov_asym(*cov_production_asym_Run1);
    cov_asym +=  *cov_detection_asym_Run1 ;
    cov_asym +=   *cov_detection_asym_Run2;


    /// decay-time bias systematics
    pull p_scale_mean_dt_Run1(paraNames,"norm_toy/pull_scale_mean_dt_Run1_*.root");
    TMatrixD* cov_scale_mean_dt_Run1 = new TMatrixD(p_scale_mean_dt_Run1.getDeltaCov("norm_toy/pull__*.root","_scale_mean_dt_Run1"));

    pull p_scale_mean_dt_Run2(paraNames,"norm_toy/pull_scale_mean_dt_Run2a_*.root");
    TMatrixD* cov_scale_mean_dt_Run2 = new TMatrixD(p_scale_mean_dt_Run2.getDeltaCov("norm_toy/pull__*.root","_scale_mean_dt_Run2"));

    pull p_scale_mean_dt_Run2_17(paraNames,"norm_toy/pull_scale_mean_dt_Run2b_*.root");
    TMatrixD* cov_scale_mean_dt_Run2_17 = new TMatrixD(p_scale_mean_dt_Run2_17.getDeltaCov("norm_toy/pull__*.root","_scale_mean_dt_Run2_17"));

    TMatrixD cov_bias(*cov_scale_mean_dt_Run1);
    cov_bias +=   *cov_scale_mean_dt_Run2;
    cov_bias +=   *cov_scale_mean_dt_Run2_17;

    /// bkg systematics 
    pull p_bkg_1(paraNames,"norm_sys_bkg1/pull__*.root");
    TMatrixD* cov_bkg_1 = new TMatrixD(p_bkg_1.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_1"));
    vector<double> vals_bkg_1 = p_bkg_1.getVals();

    pull p_bkg_2(paraNames,"norm_sys_bkg2/pull__*.root");
    TMatrixD* cov_bkg_2 = new TMatrixD(p_bkg_2.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_2"));
    vector<double> vals_bkg_2 = p_bkg_2.getVals();

    pull p_bkg_3(paraNames,"norm_sys_bkg3/pull__*.root");
    TMatrixD* cov_bkg_3 = new TMatrixD(p_bkg_3.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_3"));
    vector<double> vals_bkg_3 = p_bkg_3.getVals();

    pull p_bkg_4(paraNames,"norm_sys_bkg4/pull__*.root");
    TMatrixD* cov_bkg_4 = new TMatrixD(p_bkg_4.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_4"));
    vector<double> vals_bkg_4 = p_bkg_4.getVals();

    pull p_bkg_5(paraNames,"norm_sys_bkg5/pull__*.root");
    TMatrixD* cov_bkg_5 = new TMatrixD(p_bkg_5.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_5"));
    vector<double> vals_bkg_5 = p_bkg_5.getVals();

    pull p_bkg_6(paraNames,"norm_sys_bkg6/pull__*.root");
    TMatrixD* cov_bkg_6 = new TMatrixD(p_bkg_6.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_6"));
    vector<double> vals_bkg_6 = p_bkg_6.getVals();

//     pull p_bkg_7(paraNames,"norm_sys_bkg_7/pull__*.root");
//     TMatrixD* cov_bkg_7 = new TMatrixD(p_bkg_7.getDeltaCov("norm_taggingCalib/pull__1.root","bkg_7"));
//     vector<double> vals_bkg_7 = p_bkg_7.getVals();

    // Take maximum as systematic
    vector<TMatrixD*> covs_bkg;
    covs_bkg.push_back(cov_bkg_1);
    covs_bkg.push_back(cov_bkg_2);
    covs_bkg.push_back(cov_bkg_3);
    covs_bkg.push_back(cov_bkg_4);
    covs_bkg.push_back(cov_bkg_5);
    covs_bkg.push_back(cov_bkg_6);
//     covs_bkg.push_back(cov_bkg_7);

    TMatrixD cov_bkg_max(p.combineCov_maxVal(covs_bkg));
    //cov_bkg_max.Print();
    //covs.push_back(new TMatrixD(cov_bkg_max));

    // Take sample variance as systematic 
    vector< vector <double> > vec_vals_bkg;
    vec_vals_bkg.push_back(vals_bkg_1);
    vec_vals_bkg.push_back(vals_bkg_2);
    vec_vals_bkg.push_back(vals_bkg_3);
    vec_vals_bkg.push_back(vals_bkg_4);
    vec_vals_bkg.push_back(vals_bkg_5);
    vec_vals_bkg.push_back(vals_bkg_6);
//     vec_vals_bkg.push_back(vals_bkg_7);

    TMatrixD cov_bkg(p.sampleVariance(vec_vals_bkg));
    //cov_bkg.Print();

    /// multiple candidates systematics 
//     pull p_mc(paraNames,"norm_sys_mc/pull__*.root");
//     TMatrixD* cov_mc = new TMatrixD(p_mc.getDeltaCov("norm_taggingCalib/pull__1.root"));
//     covs.push_back(cov_mc);

    /// m,t correlations
    pull p_corr1(paraNames,"norm_toy_bkg3/pull__*.root");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"norm_toy_bkg4/pull__*.root");
    TMatrixD cov_corr2 = p_corr2.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
    }

    /// Total systematics table   
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(new TMatrixD(cov_bias));
    covs.push_back(new TMatrixD(cov_asym));
    covs.push_back(cov_corr);

    vector<string> sysNames;
    sysNames.push_back("Fit-bias");
    sysNames.push_back("Background");
    sysNames.push_back("Acceptance");
    sysNames.push_back("Resolution");
    sysNames.push_back("Decay-time bias");
    sysNames.push_back("Asymmetries");
    sysNames.push_back("Correlations");
//     sysNames.push_back("Mult.-Cand.");
    sysNames.push_back("Mom./z-Scale");

    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/timeFit/norm_taggingCalib/sys_summary_table.tex",std::ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l " ;
    for(int i =0 ; i <=covs.size() ; i++) SummaryFile << " c " ;
    SummaryFile << " | c }" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Fit Parameter & " ;
    for(int i =0 ; i <=covs.size() ; i++)  SummaryFile << sysNames[i] << " & " ;
    SummaryFile << " Total " << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";

    vector<double> errs_sys;    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile << std::fixed << std::setprecision(4) << p_data.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile << sqrt((*covs[j])[i][i]) << " & ";  
        }
	if(i==25){  
                tot += pow(0.0056,2);
		SummaryFile << 0.0056 << " & " ;
	}
	else SummaryFile << " & ";
        SummaryFile << sqrt(tot) << " \\\\ " << "\n";
	errs_sys.push_back(sqrt(tot)); 
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";

    /// Total systematics table in terms of sigma_stat   
    ofstream SummaryFile2;
    SummaryFile2.open("../../../../../TD-AnaNote/latex/tables/timeFit/norm_taggingCalib/sys_summary_table2.tex",std::ofstream::trunc);

    SummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i <=covs.size() ; i++) SummaryFile2 << " c " ;
    SummaryFile2 << " | c }" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "Fit Parameter & " ;
    for(int i =0 ; i <=covs.size() ; i++)  SummaryFile2 << sysNames[i] << " & " ;
    SummaryFile2 << " Total " << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";
    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile2 << std::fixed << std::setprecision(2) << p_data.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile2 << sqrt((*covs[j])[i][i])/errs_stat[i] << " & ";  
        }
	if(i==25){  
                tot += pow(0.0056,2);
		SummaryFile2 << 0.0056/errs_stat[25] << " & " ;
	}
	else SummaryFile2 << " & ";
        SummaryFile2 << sqrt(tot)/errs_stat[i] << " \\\\ " << "\n"; 
    }
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";


    /// Result table   
    ofstream SummaryFile3;
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/timeFit/norm_taggingCalib/result_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{l r r } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Fit Parameter} & \\multicolumn{1}{c}{Run-I} & \\multicolumn{1}{c}{Run-II} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    
    SummaryFile3 << std::fixed << std::setprecision(3);
    for(int i =0 ; i < 12 ; i++){

	double scale = 1.;
	if(i == 4 || i == 5 || i == 10 || i == 11){
		scale = 100.;
	}

        SummaryFile3 << p_data.latexNameMod(paraNames[i])  << " & " ;
        SummaryFile3 << vals[i] * scale << " $\\pm$ " << errs_stat[i] * scale << " $\\pm$ " << errs_sys[i] * scale << " & ";
        SummaryFile3 << vals[i+12] * scale<< " $\\pm$ " << errs_stat[i+12] * scale << " $\\pm$ " << errs_sys[i+12] * scale;
        SummaryFile3  << " \\\\ " << "\n"; 

	if(i == 5)SummaryFile3  << " \\\\ " << "\n"; 
    }

    SummaryFile3  << " \\\\ " << "\n";
    SummaryFile3 << p_data.latexNameMod(paraNames[24])  << " & -0.045 (fixed) & " ;
    SummaryFile3 << vals[24] * 100 << " $\\pm$ " << errs_stat[24] * 100  << " $\\pm$ " << errs_sys[24] * 100;
    SummaryFile3  << " \\\\ " << "\n"; 
    SummaryFile3 << "\\hline" << "\n";

    SummaryFile3 << std::fixed << std::setprecision(4);
    SummaryFile3 << p_data.latexNameMod(paraNames[25])  << " & \\multicolumn{2}{c}{ " ;
    SummaryFile3 << "xx.xx" << " $\\pm$ " << errs_stat[25]  << " $\\pm$ " << errs_sys[25] << " } ";
    SummaryFile3  << " \\\\ " << "\n"; 

    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";
}

int main(int argc, char** argv){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    //gStyle->SetOptFit(111);
    //gStyle->UseCurrentStyle();

//    signal();
    norm();

    return 0;
}
