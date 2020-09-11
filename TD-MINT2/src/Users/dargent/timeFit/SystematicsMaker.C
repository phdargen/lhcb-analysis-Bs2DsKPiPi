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
//     pull p_data(paraNames,"signal_test5/pull__1.root");
    //pull p_data(paraNames,"signal_new/pull__2.root");
    pull p_data(paraNames,"signal_new2/pull__2.root");
    pull p_data_err(paraNames,"signal_new2/pull__1.root");

    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data_err.getErrs();
    
    /// Stat cov from toys
    pull p_stat(paraNames,"out_signal_new_toy/pull__*.root");
    TMatrixD* cov_stat = new TMatrixD(p_stat.getStatCov());
    cov_stat->Print();

    /// Fit bias from toys
    pull p(paraNames,"out_signal_new_toy/pull__*.root");
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
    pull p_acc(paraNames,"out_signal_new_toy/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("out_signal_new_toy/pull__*.root","_acc"));
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

    pull p_res_Run2_e(paraNames,"signal_sys_res_Run2_e/pull__*.root");
    TMatrixD* cov_res_Run2_e = new TMatrixD(p_res_Run2_e.getDeltaCov("signal/pull__1.root","_res_Run2_e"));

    pull p_res_Run2_f(paraNames,"signal_sys_res_Run2_f/pull__*.root");
    TMatrixD* cov_res_Run2_f = new TMatrixD(p_res_Run2_f.getDeltaCov("signal/pull__1.root","_res_Run2_f"));

    pull p_res_Run2_g(paraNames,"signal_sys_res_Run2_g/pull__*.root");
    TMatrixD* cov_res_Run2_g = new TMatrixD(p_res_Run2_g.getDeltaCov("signal/pull__1.root","_res_Run2_g"));

    pull p_res_Run2_h(paraNames,"signal_sys_res_Run2_h/pull__*.root");
    TMatrixD* cov_res_Run2_h = new TMatrixD(p_res_Run2_h.getDeltaCov("signal/pull__1.root","_res_Run2_h"));

    pull p_res_18_a(paraNames,"signal_new_sys_res_18_a/pull__*.root");
    TMatrixD* cov_res_18_a = new TMatrixD(p_res_18_a.getDeltaCov("signal_new/pull__1.root","_res_18_a"));

    pull p_res_18_b(paraNames,"signal_new_sys_res_18_b/pull__*.root");
    TMatrixD* cov_res_18_b = new TMatrixD(p_res_18_b.getDeltaCov("signal_new/pull__1.root","_res_18_b"));


    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2,covs_res_Run2_17,covs_res_Run2_18;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);

    covs_res_Run2_17.push_back(cov_res_Run2_c);
    covs_res_Run2_17.push_back(cov_res_Run2_d);

    covs_res_Run2_18.push_back(cov_res_18_a);
    covs_res_Run2_18.push_back(cov_res_18_b);

    cov_res +=  p.combineCov_maxVal(covs_res_Run2) ;
    cov_res +=  p.combineCov_maxVal(covs_res_Run2_17) ;
    cov_res +=  p.combineCov_maxVal(covs_res_Run2_18) ;

    TMatrixD cor_res(cov_res);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++)
		 cor_res[i][j] = cov_res[i][j]/sqrt(cov_res[i][i])/sqrt(cov_res[j][j]);
    

    /// decay time bias systematics 
    pull p_bias(paraNames,"signal_decayTimeBias/pull__2.root");
    TMatrixD cov_bias = p_bias.getDeltaCov("signal_decayTimeBias/pull__1.root","_bias");

    pull p_bias_18(paraNames,"out_signal_new_toy/pull_offset_mean_dt_Run218_*.root");
    TMatrixD cov_bias_18 = p_bias_18.getDeltaCov("out_signal_new_toy/pull__*.root","_bias18");

    cov_bias += cov_bias_18;

    /// dms systematics 
//     pull p_dm(paraNames,"signal_toy_ub2/pull_dm_*.root");
//     TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("signal_toy_ub/pull__*.root","_dm"));
    pull p_dm(paraNames,"out_signal_new_toy/pull_dm_*.root");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("out_signal_new_toy/pull__*.root","_dm"));
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

//     pull p_bkg_bs1(paraNames,"signal_toy_bkg_bs2/pull__*.root");
    pull p_bkg_bs1(paraNames,"out_signal_new_toy_bkg/sw_pull__*.root");
    vector<double> mean_bkg_bs1 = p_bkg_bs1.getPullMean();
    vector<double> sigma_bkg_bs1 = p_bkg_bs1.getPullSigma();

//     pull p_bkg_bs2(paraNames,"signal_toy_bkg_bs2/signal_pull__*.root");*/
    pull p_bkg_bs2(paraNames,"out_signal_new_toy_bkg/pull__*.root");   
    vector<double> mean_bkg_bs2 = p_bkg_bs2.getPullMean();
    vector<double> sigma_bkg_bs2 = p_bkg_bs2.getPullSigma();

    TMatrixD* cov_bkg_bs = new TMatrixD(cov_bkg);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++){
	if(i==j){
		 (*cov_bkg_bs)[i][j] = pow(mean_bkg_bs1[i]-mean_bkg_bs2[i],2)*errs_stat[i] * errs_stat[j];;
		 (*cov_bkg_bs)[i][j] += pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2)*errs_stat[i] * errs_stat[j];;
		cout << paraNames[i] << endl;
		cout << mean_bkg_bs1[i]-mean_bkg_bs2[i] << endl;
		cout << sigma_bkg_bs1[i]-sigma_bkg_bs2[i] << endl << endl;
	}
	else (*cov_bkg_bs)[i][j] = 0.;
     }

//      TMatrixD* cov_bkg_bs = new TMatrixD(p_bkg_bs2.getDeltaCov("out_signal_new_toy_bkg3/pull__*.root","_bkg_bs"));
     cov_bkg+= *cov_bkg_bs;


    /// m,t correlations
//      pull p_corr(paraNames,"out_signal_new_toy_bkg/sw_pull__*.root");
//      TMatrixD* cov_corr = new TMatrixD(p_corr.getDeltaCov("out_signal_new_toy_bkg2/sw_pull__*.root","_corr"));

//     pull p_corr1(paraNames,"signal_toy_bkg7/pull__*.root");
//     TMatrixD cov_corr1 = p_corr1.getCov();
// 
//     pull p_corr2(paraNames,"signal_toy_bkg8/pull__*.root");
//     TMatrixD cov_corr2 = p_corr2.getCov();
// 
//     pull p_corr3(paraNames,"signal_toy_bkg_bs0/pull__*.root");
//     TMatrixD cov_corr3 = p_corr3.getCov();
// 
//     pull p_corr4(paraNames,"signal_toy_bkg_bs1/pull__*.root");
//     TMatrixD cov_corr4 = p_corr4.getCov();
// 
//     TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
//     TMatrixD* cov_corr34 = new TMatrixD(p_corr3.getAbsDiff(cov_corr3,cov_corr4));
// 
//     for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
// 	if((*cov_corr)[i][j]>(*cov_corr34)[i][j]) (*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
// 	else (*cov_corr)[i][j] = (*cov_corr34)[i][j] * errs_stat[i] * errs_stat[j];
//     }

    pull p_corr1(paraNames,"out_signal_new_toy_bkg/sw_pull__*.root");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"out_signal_new_toy_bkg2/sw_pull__*.root");
    TMatrixD cov_corr2 = p_corr2.getCov();

    pull p_corr3(paraNames,"out_signal_new_toy_bkg3/sw_pull__*.root");
    TMatrixD cov_corr3 = p_corr3.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
     TMatrixD* cov_corr34 = new TMatrixD(p_corr3.getAbsDiff(cov_corr3,cov_corr1));

    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
// 	if((*cov_corr)[i][j]>(*cov_corr34)[i][j]) (*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
// 	else (*cov_corr)[i][j] = (*cov_corr34)[i][j] * errs_stat[i] * errs_stat[j];
    }

    /// Total sys corr
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(cov_corr);
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(new TMatrixD(cov_bias));
    covs.push_back(new TMatrixD(cov_asym));
    covs.push_back(cov_dm);

    TMatrixD cov_sys_tot(paraNames.size(),paraNames.size());
    TMatrixD cor_sys_tot(paraNames.size(),paraNames.size());

    vector<string> sysNames;
    sysNames.push_back("Fit bias");
    sysNames.push_back("Background");
    sysNames.push_back("Correlations");
    sysNames.push_back("Acceptance");
    sysNames.push_back("Resolution");
    sysNames.push_back("Decay-time bias");
    sysNames.push_back("Asymmetries");
    sysNames.push_back("$\\Delta m_{s}$");

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
        SummaryFile3 << std::fixed << std::setprecision(3) << vals[i] << " $\\pm$ " << errs_stat[i]  << " $\\pm$ " << errs_sys[i];
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
//     pull p_data(paraNames,"norm_taggingCalib/pull__1.root");
//    pull p_data(paraNames,"norm_decayTimeBias/pull__1.root");
     pull p_data(paraNames,"norm_new/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();
    
    /// Fit bias from toys
    pull p(paraNames,"out_norm_new_toy/pull__*.root");
    TMatrixD* cov = new TMatrixD(p.getCov("",5,false));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov)[i][j] = (*cov)[i][j] * errs_stat[i] * errs_stat[j];
    }
    covs.push_back(cov);
   
    /// Acc systematics (with cholesky)
//     pull p_acc_chol(paraNames,"norm_toy/pullAccChol_*.root");
//     TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("norm_toy/pull__*.root","_accChol",50));
//     covs.push_back(cov_acc_chol);

    pull p_acc(paraNames,"out_norm_new_toy/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("out_norm_new_toy/pull__*.root","_acc",0.25));

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

    pull p_res_Run2_e(paraNames,"norm_sys_res_Run2_e/pull__*.root");
    TMatrixD* cov_res_Run2_e = new TMatrixD(p_res_Run2_e.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_e"));

    pull p_res_Run2_f(paraNames,"norm_sys_res_Run2_f/pull__*.root");
    TMatrixD* cov_res_Run2_f = new TMatrixD(p_res_Run2_f.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_f"));

    pull p_res_Run2_g(paraNames,"norm_sys_res_Run2_g/pull__*.root");
    TMatrixD* cov_res_Run2_g = new TMatrixD(p_res_Run2_g.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_g"));

    pull p_res_Run2_h(paraNames,"norm_sys_res_Run2_h/pull__*.root");
    TMatrixD* cov_res_Run2_h = new TMatrixD(p_res_Run2_h.getDeltaCov("norm_taggingCalib/pull__1.root","_res_Run2_h"));

    pull p_res_18_a(paraNames,"norm_new_sys_res_18_a/pull__*.root");
    TMatrixD* cov_res_18_a = new TMatrixD(p_res_18_a.getDeltaCov("norm_new/pull__1.root","_res_18_a"));

    pull p_res_18_b(paraNames,"norm_new_sys_res_18_b/pull__*.root");
    TMatrixD* cov_res_18_b = new TMatrixD(p_res_18_b.getDeltaCov("norm_new/pull__1.root","_res_18_b"));


    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p_res_Run1_a.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2,covs_res_Run2_17,covs_res_Run2_18;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    //covs_res_Run2.push_back(cov_res_Run2_e);
    //covs_res_Run2.push_back(cov_res_Run2_f);

    covs_res_Run2.push_back(cov_res_Run2_c);
    covs_res_Run2.push_back(cov_res_Run2_d);
    //covs_res_Run2_17.push_back(cov_res_Run2_g);
    //covs_res_Run2_17.push_back(cov_res_Run2_h);

    covs_res_Run2.push_back(cov_res_18_a);
    covs_res_Run2.push_back(cov_res_18_b);

    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2) ;
    //cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2_17) ;
    //cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2_18) ;
    
    /// asymmetry systematics
    pull p_production_asym_Run1(paraNames,"norm_toy/pull_production_asym_Run1_*.root");
    TMatrixD* cov_production_asym_Run1 = new TMatrixD(p_production_asym_Run1.getDeltaCov("norm_toy/pull__*.root","_production_asym_Run1",0.25,true));

    pull p_detection_asym_Run1(paraNames,"norm_toy/pull_detection_asym_Run1_*.root");
    TMatrixD* cov_detection_asym_Run1 = new TMatrixD(p_detection_asym_Run1.getDeltaCov("norm_toy/pull__*.root","_detection_asym_Run1",0.25,true));

    pull p_detection_asym_Run2(paraNames,"norm_toy/pull_detection_asym_Run2_*.root");
    TMatrixD* cov_detection_asym_Run2 = new TMatrixD(p_detection_asym_Run2.getDeltaCov("norm_toy/pull__*.root","_detection_asym_Run2",0.25,true));

    TMatrixD cov_asym(*cov_production_asym_Run1);
    cov_asym +=  *cov_detection_asym_Run1 ;
    cov_asym +=   *cov_detection_asym_Run2;


    /// decay-time bias systematics
    pull p_scale_mean_dt_Run1(paraNames,"norm_toy/pull_scale_mean_dt_Run1_*.root");
    TMatrixD* cov_scale_mean_dt_Run1 = new TMatrixD(p_scale_mean_dt_Run1.getDeltaCov("norm_toy/pull__*.root","_scale_mean_dt_Run1"));

    pull p_scale_mean_dt_Run2(paraNames,"norm_toy/pull_scale_mean_dt_Run2a_*.root");
    TMatrixD* cov_scale_mean_dt_Run2 = new TMatrixD(p_scale_mean_dt_Run2.getDeltaCov("norm_toy/pull__*.root","_scale_mean_dt_Run2"));

    pull p_scale_mean_dt_Run2_17(paraNames,"norm_toy/pull_offset_mean_dt_Run217_*.root");
    TMatrixD* cov_scale_mean_dt_Run2_17 = new TMatrixD(p_scale_mean_dt_Run2_17.getDeltaCov("norm_toy/pull__*.root","_offset_mean_dt_Run217"));

    pull p_scale_mean_dt_Run2_18(paraNames,"out_norm_new_toy/pull_offset_mean_dt_Run218_*.root");
    TMatrixD* cov_scale_mean_dt_Run2_18 = new TMatrixD(p_scale_mean_dt_Run2_18.getDeltaCov("out_norm_new_toy/pull__*.root","_offset_mean_dt_Run218"));

    TMatrixD cov_bias(*cov_scale_mean_dt_Run1);
    cov_bias +=   *cov_scale_mean_dt_Run2;
    cov_bias +=   *cov_scale_mean_dt_Run2_17;
    cov_bias +=   *cov_scale_mean_dt_Run2_18;

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

//     pull p_bkg_bs1(paraNames,"norm_toy_bkg_bs2/pull__*.root");
    pull p_bkg_bs1(paraNames,"out_norm_new_toy_bkg/sw_pull__*.root");
    vector<double> mean_bkg_bs1 = p_bkg_bs1.getPullMean();
    vector<double> sigma_bkg_bs1 = p_bkg_bs1.getPullSigma();

//     pull p_bkg_bs2(paraNames,"norm_toy_bkg_bs2/signal_pull__*.root");
    pull p_bkg_bs2(paraNames,"out_norm_new_toy_bkg/pull__*.root");
    vector<double> mean_bkg_bs2 = p_bkg_bs2.getPullMean();
    vector<double> sigma_bkg_bs2 = p_bkg_bs2.getPullSigma();

    TMatrixD* cov_bkg_bs = new TMatrixD(cov_bkg);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++){
	if(i==j){
		 (*cov_bkg_bs)[i][j] = pow(mean_bkg_bs1[i]-mean_bkg_bs2[i],2)*errs_stat[i] * errs_stat[j];;
		 (*cov_bkg_bs)[i][j] += pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2)*errs_stat[i] * errs_stat[j];;
		cout << paraNames[i] << endl;
		cout << mean_bkg_bs1[i]-mean_bkg_bs2[i] << endl;
		cout << sigma_bkg_bs1[i]-sigma_bkg_bs2[i] << endl << endl;
	}
	else (*cov_bkg_bs)[i][j] = 0.;
     }
     cov_bkg+= *cov_bkg_bs;


    /// multiple candidates systematics 
//     pull p_mc(paraNames,"norm_sys_mc/pull__*.root");
//     TMatrixD* cov_mc = new TMatrixD(p_mc.getDeltaCov("norm_taggingCalib/pull__1.root"));
//     covs.push_back(cov_mc);

    /// m,t correlations
    //pull p_corr1(paraNames,"norm_toy_bkg3/pull__*.root");
    pull p_corr1(paraNames,"out_norm_new_toy_bkg/sw_pull__*.root");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"norm_toy_bkg4/pull__*.root");
    //pull p_corr2(paraNames,"out_norm_new_toy_bkg3/sw_pull__*.root");
    TMatrixD cov_corr2 = p_corr2.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
    }

    /// only dms
    TMatrixD cov_closure(paraNames.size(),paraNames.size());
    TMatrixD cov_velo(paraNames.size(),paraNames.size());
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	if(i==25 && j ==25){ 
		cov_closure[i][j] = pow(0.0027,2);
		cov_velo[i][j] = pow(0.003,2);
	}
	else { 
		cov_closure[i][j] = 0.;
		cov_velo[i][j] = 0.;		
	}
    }
    cov_bias +=   cov_velo;

    /// Total systematics table   
    covs.push_back(new TMatrixD(cov_closure));
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(cov_corr);
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(new TMatrixD(cov_bias));
    //covs.push_back(new TMatrixD(cov_velo));
    covs.push_back(new TMatrixD(cov_asym));

    vector<string> sysNames;
    sysNames.push_back("Fit-bias");
    sysNames.push_back("Closure");
    sysNames.push_back("Bkg.");
    sysNames.push_back("Correlations");
    sysNames.push_back("Acceptance");
    sysNames.push_back("Resolution");
    sysNames.push_back("Decay-time bias");
    //sysNames.push_back("VELO-Misalign.");
    sysNames.push_back("Asymmetries");
//     sysNames.push_back("Mult.-Cand.");
    sysNames.push_back("z-Scale");

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
            if((*covs[j])[i][i]>0)SummaryFile << sqrt((*covs[j])[i][i]) << " & ";  
	    else SummaryFile << " & ";
        }
	if(i==25){  
                tot += pow(vals[25]*0.0002,2);
		SummaryFile << vals[25]*0.0002 << " & " ;
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
            if((*covs[j])[i][i]>0)SummaryFile2 << sqrt((*covs[j])[i][i])/errs_stat[i] << " & ";  
	    else SummaryFile2 << " & ";
        }
	if(i==25){  
                tot += pow(vals[25]*0.0002,2);
		SummaryFile2 << vals[25]*0.0002/errs_stat[25] << " & " ;
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
    SummaryFile3 << vals[25] << " $\\pm$ " << errs_stat[25]  << " $\\pm$ " << errs_sys[25] << " } ";
    SummaryFile3  << " \\\\ " << "\n"; 

    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";


    for(int i =0 ; i < vals.size() ; i++){
        cout << p_data.latexNameMod(paraNames[i])  << " = " ;
        cout << vals[i] << " $\\pm$ " << sqrt(pow(errs_stat[i],2) + pow(errs_sys[i],2) ) << endl;
    }

}



void crossCheck_norm() {

    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("dm");

    pull p_data(paraNames,"norm_decayTimeBias/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();

    pull p_1(paraNames,"norm_crossCheck_trigger1/pull__1.root");
    vector<double> vals1 = p_1.getVals();
    vector<double> errs_stat1 = p_1.getErrs();

    pull p_2(paraNames,"norm_crossCheck_trigger2/pull__1.root");
    vector<double> vals2 = p_2.getVals();
    vector<double> errs_stat2 = p_2.getErrs();

    pull p_3(paraNames,"norm_crossCheck_run1/pull__1.root");
    vector<double> vals3 = p_3.getVals();
    vector<double> errs_stat3 = p_3.getErrs();

    pull p_4(paraNames,"norm_crossCheck_run2/pull__1.root");
    vector<double> vals4 = p_4.getVals();
    vector<double> errs_stat4 = p_4.getErrs();

    pull p_5(paraNames,"norm_crossCheck_Ds1/pull__1.root");
    vector<double> vals5 = p_5.getVals();
    vector<double> errs_stat5 = p_5.getErrs();

    pull p_6(paraNames,"norm_crossCheck_Ds2/pull__1.root");
    vector<double> vals6 = p_6.getVals();
    vector<double> errs_stat6 = p_6.getErrs();

    pull p_7(paraNames,"norm_crossCheck_OS/pull__1.root");
    vector<double> vals7 = p_7.getVals();
    vector<double> errs_stat7 = p_7.getErrs();

    pull p_8(paraNames,"norm_crossCheck_SS/pull__1.root");
    vector<double> vals8 = p_8.getVals();
    vector<double> errs_stat8 = p_8.getErrs();

    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/timeFit/norm_taggingCalib/crossCheck_table.tex",std::ofstream::trunc);

    for(int i =0 ; i <paraNames.size() ; i++){
        SummaryFile << p_1.latexName(paraNames[i])  << " & " ;
        SummaryFile << std::fixed << std::setprecision(3) << vals[i] << " $\\pm$ " << errs_stat[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals1[i] << " $\\pm$ " << errs_stat1[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals2[i] << " $\\pm$ " << errs_stat2[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals3[i] << " $\\pm$ " << errs_stat3[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals4[i] << " $\\pm$ " << errs_stat4[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals5[i] << " $\\pm$ " << errs_stat5[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals6[i] << " $\\pm$ " << errs_stat6[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals7[i] << " $\\pm$ " << errs_stat7[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals8[i] << " $\\pm$ " << errs_stat8[i];

        SummaryFile  << " \\\\ " << "\n"; 
    }
    
}

void crossCheck_signal() {

    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("C");
    paraNames.push_back("D");
    paraNames.push_back("D_bar");
    paraNames.push_back("S");
    paraNames.push_back("S_bar");

    pull p_data(paraNames,"signal_decayTimeBias/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();

    pull p_1(paraNames,"signal_crossCheck_trigger1/pull__1.root");
    vector<double> vals1 = p_1.getVals();
    vector<double> errs_stat1 = p_1.getErrs();

    pull p_2(paraNames,"signal_crossCheck_trigger2/pull__1.root");
    vector<double> vals2 = p_2.getVals();
    vector<double> errs_stat2 = p_2.getErrs();

    pull p_3(paraNames,"signal_crossCheck_run1/pull__1.root");
    vector<double> vals3 = p_3.getVals();
    vector<double> errs_stat3 = p_3.getErrs();

    pull p_4(paraNames,"signal_crossCheck_run2/pull__1.root");
    vector<double> vals4 = p_4.getVals();
    vector<double> errs_stat4 = p_4.getErrs();

    pull p_5(paraNames,"signal_crossCheck_Ds1/pull__1.root");
    vector<double> vals5 = p_5.getVals();
    vector<double> errs_stat5 = p_5.getErrs();

    pull p_6(paraNames,"signal_crossCheck_Ds2/pull__1.root");
    vector<double> vals6 = p_6.getVals();
    vector<double> errs_stat6 = p_6.getErrs();

    pull p_7(paraNames,"signal_crossCheck_OS/pull__1.root");
    vector<double> vals7 = p_7.getVals();
    vector<double> errs_stat7 = p_7.getErrs();

    pull p_8(paraNames,"signal_crossCheck_SS/pull__1.root");
    vector<double> vals8 = p_8.getVals();
    vector<double> errs_stat8 = p_8.getErrs();

    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/timeFit/signal/crossCheck_table.tex",std::ofstream::trunc);

    for(int i =0 ; i <paraNames.size() ; i++){
        SummaryFile << p_1.latexName(paraNames[i])  << " & " ;
        SummaryFile << std::fixed << std::setprecision(3) << vals[i] << " $\\pm$ " << errs_stat[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals1[i] << " $\\pm$ " << errs_stat1[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals2[i] << " $\\pm$ " << errs_stat2[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals3[i] << " $\\pm$ " << errs_stat3[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals4[i] << " $\\pm$ " << errs_stat4[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals5[i] << " $\\pm$ " << errs_stat5[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals6[i] << " $\\pm$ " << errs_stat6[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(3) << vals7[i] << " $\\pm$ " << errs_stat7[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(3) << vals8[i] << " $\\pm$ " << errs_stat8[i];

        SummaryFile  << " \\\\ " << "\n"; 
    }  
}


int main(int argc, char** argv){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    //gStyle->SetOptFit(111);
    //gStyle->UseCurrentStyle();

    signal();
   //norm();

//  crossCheck_norm();
// crossCheck_signal();

    return 0;
}