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

void fitParams(){

   /// Fit parameters    
    vector<TString> paraNames;

    paraNames.push_back("Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Amp");
    paraNames.push_back("Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Phase");
    paraNames.push_back("Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Amp");
    paraNames.push_back("Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Phase");

    paraNames.push_back("a_K1_1400_Amp");
    paraNames.push_back("a_K1_1400_Phase");
    paraNames.push_back("abar_K1_1400_Amp");
    paraNames.push_back("abar_K1_1400_Phase");

    paraNames.push_back("a_Ks_1410_Amp");
    paraNames.push_back("a_Ks_1410_Phase");
    paraNames.push_back("Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Amp");
    paraNames.push_back("Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Phase");

    paraNames.push_back("abar_K_1460_Amp");
    paraNames.push_back("abar_K_1460_Phase");

    paraNames.push_back("a_NS_Ks_Amp");
    paraNames.push_back("a_NS_Ks_Phase");
    paraNames.push_back("abar_NS_Ks_Amp");
    paraNames.push_back("abar_NS_Ks_Phase");

    paraNames.push_back("abar_NS_rho_Amp");
    paraNames.push_back("abar_NS_rho_Phase");

    paraNames.push_back("mass_K_1__1400_p");
    paraNames.push_back("width_K_1__1400_p");
    paraNames.push_back("mass_Ks_1410_p");
    paraNames.push_back("width_Ks_1410_p");

    paraNames.push_back("r");
    paraNames.push_back("delta");
    paraNames.push_back("gamma");
 

    vector<TMatrixD*> covs;
    /// Data fit
    pull p_data(paraNames,"signal/pull__1.root");
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data.getErrs();
    
    /// Stat cov from toys
//     pull p_stat(paraNames,"signal_toy/pull__*.root");
//     TMatrixD* cov_stat = new TMatrixD(p_stat.getStatCov());

    /// Fit bias from toys
//     pull p(paraNames,"signal_toy/pull__*.root");
//     TMatrixD* cov = new TMatrixD(p.getCov());
//     covs.push_back(cov);

    /// Acc systematics (with cholesky)
//     pull p_acc_chol(paraNames,"signal_toy/pullAccChol_*.root","MinuitParameterSetNtp",false, false, 900);
//     TMatrixD* cov_acc_chol = new TMatrixD(p_acc_chol.getDeltaCovChol("signal_toy/pull__*.root","_accChol",50));
//     covs.push_back(cov_acc_chol);

    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"signal_sys_res_Run1_a/pull__*.root");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("signal/pull__1.root","_res_Run1_a"));

    pull p_res_Run1_b(paraNames,"signal_sys_res_Run1_b/pull__*.root");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("signal/pull__1.root","_res_Run1_b"));

    pull p_res_Run2_a(paraNames,"signal_sys_res_Run2_a/pull__*.root");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("signal/pull__1.root","_res_Run2_a"));

    pull p_res_Run2_b(paraNames,"signal_sys_res_Run2_b/pull__*.root");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("signal/pull__1.root","_res_Run2_b"));

    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p_res_Run1_a.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2) ;
    
    covs.push_back(new TMatrixD(cov_res));


    /// dms systematics 
    pull p_dm(paraNames,"signal_toy/pull_dm_*.root");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("signal_toy/pull__*.root","_dm"));
    covs.push_back(cov_dm);

    /// asymmetry systematics
    pull p_production_asym_Run1(paraNames,"signal_toy/pull_production_asym_Run1_*.root");
    TMatrixD* cov_production_asym_Run1 = new TMatrixD(p_production_asym_Run1.getDeltaCov("signal_toy/pull__*.root","_production_asym_Run1"));

    pull p_production_asym_Run2(paraNames,"signal_toy/pull_production_asym_Run2_*.root");
    TMatrixD* cov_production_asym_Run2 = new TMatrixD(p_production_asym_Run2.getDeltaCov("signal_toy/pull__*.root","_production_asym_Run2"));

    pull p_detection_asym_Run1(paraNames,"signal_toy/pull_detection_asym_Run1_*.root");
    TMatrixD* cov_detection_asym_Run1 = new TMatrixD(p_detection_asym_Run1.getDeltaCov("signal_toy/pull__*.root","_detection_asym_Run1"));

    pull p_detection_asym_Run2(paraNames,"signal_toy/pull_detection_asym_Run2_*.root");
    TMatrixD* cov_detection_asym_Run2 = new TMatrixD(p_detection_asym_Run2.getDeltaCov("signal_toy/pull__*.root","_detection_asym_Run2"));

    TMatrixD cov_asym(*cov_production_asym_Run1);
    cov_asym +=  *cov_production_asym_Run2 ;
    cov_asym +=  *cov_detection_asym_Run1 ;
    cov_asym +=   *cov_detection_asym_Run2;

    covs.push_back(new TMatrixD(cov_asym));

    /// bkg systematics 
//     pull p_bkg_1(paraNames,"signal_sys_bkg_1/pull__*.root");
//     TMatrixD* cov_bkg_1 = new TMatrixD(p_bkg_1.getDeltaCov("signal/pull__1.root","bkg_1"));
//     vector<double> vals_bkg_1 = p_bkg_1.getVals();

    pull p_bkg_2(paraNames,"signal_sys_bkg_2/pull__*.root");
    TMatrixD* cov_bkg_2 = new TMatrixD(p_bkg_2.getDeltaCov("signal/pull__1.root","bkg_2"));
    vector<double> vals_bkg_2 = p_bkg_2.getVals();

    pull p_bkg_3(paraNames,"signal_sys_bkg_3/pull__*.root");
    TMatrixD* cov_bkg_3 = new TMatrixD(p_bkg_3.getDeltaCov("signal/pull__1.root","bkg_3"));
    vector<double> vals_bkg_3 = p_bkg_3.getVals();

    pull p_bkg_4(paraNames,"signal_sys_bkg_4/pull__*.root");
    TMatrixD* cov_bkg_4 = new TMatrixD(p_bkg_4.getDeltaCov("signal/pull__1.root","bkg_4"));
    vector<double> vals_bkg_4 = p_bkg_4.getVals();

    pull p_bkg_5(paraNames,"signal_sys_bkg_5/pull__*.root");
    TMatrixD* cov_bkg_5 = new TMatrixD(p_bkg_5.getDeltaCov("signal/pull__1.root","bkg_5"));
    vector<double> vals_bkg_5 = p_bkg_5.getVals();

    pull p_bkg_6(paraNames,"signal_sys_bkg_6/pull__*.root");
    TMatrixD* cov_bkg_6 = new TMatrixD(p_bkg_6.getDeltaCov("signal/pull__1.root","bkg_6"));
    vector<double> vals_bkg_6 = p_bkg_6.getVals();

    pull p_bkg_7(paraNames,"signal_sys_bkg_7/pull__*.root");
    TMatrixD* cov_bkg_7 = new TMatrixD(p_bkg_7.getDeltaCov("signal/pull__1.root","bkg_7"));
    vector<double> vals_bkg_7 = p_bkg_7.getVals();

    // Take maximum as systematic
    vector<TMatrixD*> covs_bkg;
    //covs_bkg.push_back(cov_bkg_1);
    covs_bkg.push_back(cov_bkg_2);
    covs_bkg.push_back(cov_bkg_3);
    covs_bkg.push_back(cov_bkg_4);
    covs_bkg.push_back(cov_bkg_5);
    covs_bkg.push_back(cov_bkg_6);
    covs_bkg.push_back(cov_bkg_7);

    TMatrixD cov_bkg_max(p_bkg_2.combineCov_maxVal(covs_bkg));
//     cov_bkg_max.Print();
    //covs.push_back(new TMatrixD(cov_bkg_max));

    // Take sample variance as systematic 
    vector< vector <double> > vec_vals_bkg;
    //vec_vals_bkg.push_back(vals_bkg_1);
    vec_vals_bkg.push_back(vals_bkg_2);
    vec_vals_bkg.push_back(vals_bkg_3);
    vec_vals_bkg.push_back(vals_bkg_4);
    vec_vals_bkg.push_back(vals_bkg_5);
    vec_vals_bkg.push_back(vals_bkg_6);
    vec_vals_bkg.push_back(vals_bkg_7);

    TMatrixD cov_bkg(p_bkg_2.sampleVariance(vec_vals_bkg));
//     cov_bkg.Print();
    covs.push_back(new TMatrixD(cov_bkg));


    /// Lineshape models systematics
    pull p_ls_1(paraNames,"signal_sys1/pull_*.root");
    TMatrixD* cov_ls_1 = new TMatrixD(p_ls_1.getDeltaCov("signal/pull__1.root","ls_1"));

    pull p_ls_2(paraNames,"signal_sys2/pull_*.root");
    TMatrixD* cov_ls_2 = new TMatrixD(p_ls_2.getDeltaCov("signal/pull__1.root","ls_2"));

    pull p_ls_3(paraNames,"signal_sys3/pull_*.root");
    TMatrixD* cov_ls_3 = new TMatrixD(p_ls_3.getDeltaCov("signal/pull__1.root","ls_3"));

    pull p_ls_4(paraNames,"signal_sys4/pull_*.root");
    TMatrixD* cov_ls_4 = new TMatrixD(p_ls_4.getDeltaCov("signal/pull__1.root","ls_4"));

    pull p_ls_5(paraNames,"signal_sys5/pull_*.root");
    TMatrixD* cov_ls_5 = new TMatrixD(p_ls_5.getDeltaCov("signal/pull__1.root","ls_5"));

    pull p_ls_6(paraNames,"signal_sys6/pull_*.root");
    TMatrixD* cov_ls_6 = new TMatrixD(p_ls_6.getDeltaCov("signal/pull__1.root","ls_6"));

    // Add
    TMatrixD cov_ls(*cov_ls_1);
    cov_ls +=  *cov_ls_2 ;
    cov_ls +=  *cov_ls_3 ;
    cov_ls +=  *cov_ls_4 ;
    cov_ls +=  *cov_ls_5 ;
    cov_ls +=  *cov_ls_6 ;
 
    covs.push_back(new TMatrixD(cov_ls));


    /// Resonance parameters systematics
    pull p_rp_1(paraNames,"signal_toy/pull_mass_K1_1270_*.root");
    TMatrixD* cov_rp_1 = new TMatrixD(p_rp_1.getDeltaCov("signal_toy/pull__*.root","_mass_K1_1270"));
    pull p_rp_2(paraNames,"signal_toy/pull_width_K1_1270_*.root");
    TMatrixD* cov_rp_2 = new TMatrixD(p_rp_2.getDeltaCov("signal_toy/pull__*.root","_width_K1_1270"));

    pull p_rp_3(paraNames,"signal_toy/pull_mass_K_1460_*.root");
    TMatrixD* cov_rp_3 = new TMatrixD(p_rp_3.getDeltaCov("signal_toy/pull__*.root","_mass_K_1460"));
    pull p_rp_4(paraNames,"signal_toy/pull_width_K_1460_*.root");
    TMatrixD* cov_rp_4 = new TMatrixD(p_rp_4.getDeltaCov("signal_toy/pull__*.root","_width_K_1460"));

    pull p_rp_5(paraNames,"signal_toy/pull_mass_Ks_*.root");
    TMatrixD* cov_rp_5 = new TMatrixD(p_rp_5.getDeltaCov("signal_toy/pull__*.root","_mass_Ks"));
    pull p_rp_6(paraNames,"signal_toy/pull_width_Ks_*.root");
    TMatrixD* cov_rp_6 = new TMatrixD(p_rp_6.getDeltaCov("signal_toy/pull__*.root","_width_Ks"));

    pull p_rp_7(paraNames,"signal_toy/pull_mass_rho_*.root");
    TMatrixD* cov_rp_7 = new TMatrixD(p_rp_7.getDeltaCov("signal_toy/pull__*.root","_mass_rho"));
    pull p_rp_8(paraNames,"signal_toy/pull_width_rho_*.root");
    TMatrixD* cov_rp_8 = new TMatrixD(p_rp_8.getDeltaCov("signal_toy/pull__*.root","_width_rho"));

    pull p_rp_9(paraNames,"signal_toy/pull_mass_K0s_*.root");
    TMatrixD* cov_rp_9 = new TMatrixD(p_rp_9.getDeltaCov("signal_toy/pull__*.root","_mass_K0s"));
    pull p_rp_10(paraNames,"signal_toy/pull_width_K0s_*.root");
    TMatrixD* cov_rp_10 = new TMatrixD(p_rp_10.getDeltaCov("signal_toy/pull__*.root","_width_K0s"));

    // Add
    TMatrixD cov_rp(*cov_rp_1);
    cov_rp +=  *cov_rp_2 ;
    cov_rp +=  *cov_rp_3 ;
//     cov_rp +=  *cov_rp_4 ;
    cov_rp +=  *cov_rp_5 ;
    cov_rp +=  *cov_rp_6 ;
//     cov_rp +=  *cov_rp_7 ;
    cov_rp +=  *cov_rp_8 ;
    cov_rp +=  *cov_rp_9 ;
    cov_rp +=  *cov_rp_10 ;
 
    covs.push_back(new TMatrixD(cov_rp));

    /// Form factor
    pull p_f_1(paraNames,"signal_toy/pull_BW_radius_*.root");
    TMatrixD* cov_f_1 = new TMatrixD(p_f_1.getDeltaCov("signal_toy/pull__*.root"));
    covs.push_back(cov_f_1);

    /// Phsp-Acc systematics
    pull p_phsp_acc_1(paraNames,"signal_sys8/pull_*.root");
    TMatrixD* cov_phsp_acc_1 = new TMatrixD(p_phsp_acc_1.getDeltaCov("signal/pull__1.root","phsp_acc_1"));
    covs.push_back(cov_phsp_acc_1);


    /// Alternative amp models 
    pull p_alt_1(paraNames,"signal_alt0/pull__*.root");
    TMatrixD* cov_alt_1 = new TMatrixD(p_alt_1.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_1 = p_alt_1.getVals();

    pull p_alt_2(paraNames,"signal_alt2/pull__*.root");
    TMatrixD* cov_alt_2 = new TMatrixD(p_alt_2.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_2 = p_alt_2.getVals();

    pull p_alt_3(paraNames,"signal_alt3/pull__*.root");
    TMatrixD* cov_alt_3 = new TMatrixD(p_alt_3.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_3 = p_alt_3.getVals();

    pull p_alt_4(paraNames,"signal_alt5/pull__*.root");
    TMatrixD* cov_alt_4 = new TMatrixD(p_alt_4.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_4 = p_alt_4.getVals();

    pull p_alt_5(paraNames,"signal_alt6/pull__*.root");
    TMatrixD* cov_alt_5 = new TMatrixD(p_alt_5.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_5 = p_alt_5.getVals();

    pull p_alt_6(paraNames,"signal_alt8/pull__*.root");
    TMatrixD* cov_alt_6 = new TMatrixD(p_alt_6.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_6 = p_alt_6.getVals();

    pull p_alt_7(paraNames,"signal_alt11/pull__*.root");
    TMatrixD* cov_alt_7 = new TMatrixD(p_alt_7.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_7 = p_alt_7.getVals();

    pull p_alt_8(paraNames,"signal_alt12/pull__*.root");
    TMatrixD* cov_alt_8 = new TMatrixD(p_alt_8.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_8 = p_alt_8.getVals();

    pull p_alt_9(paraNames,"signal_alt15/pull__*.root");
    TMatrixD* cov_alt_9 = new TMatrixD(p_alt_9.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_9 = p_alt_9.getVals();

    pull p_alt_10(paraNames,"signal_alt16/pull__*.root");
    TMatrixD* cov_alt_10 = new TMatrixD(p_alt_10.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_10 = p_alt_10.getVals();

    pull p_alt_11(paraNames,"signal_alt10/pull__*.root");
    TMatrixD* cov_alt_11 = new TMatrixD(p_alt_11.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_11 = p_alt_11.getVals();

    pull p_alt_12(paraNames,"signal_alt18/pull__*.root");
    TMatrixD* cov_alt_12 = new TMatrixD(p_alt_12.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_12 = p_alt_12.getVals();

    pull p_alt_13(paraNames,"signal_alt19/pull__*.root");
    TMatrixD* cov_alt_13 = new TMatrixD(p_alt_13.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_13 = p_alt_13.getVals();

    pull p_alt_14(paraNames,"signal_alt20/pull__*.root");
    TMatrixD* cov_alt_14 = new TMatrixD(p_alt_14.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_14 = p_alt_14.getVals();

    pull p_alt_15(paraNames,"signal_alt22/pull__*.root");
    TMatrixD* cov_alt_15 = new TMatrixD(p_alt_15.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_15 = p_alt_15.getVals();

    pull p_alt_16(paraNames,"signal_alt24/pull__*.root");
    TMatrixD* cov_alt_16 = new TMatrixD(p_alt_16.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_16 = p_alt_16.getVals();


    // Take maximum as systematic
    vector<TMatrixD*> covs_alt;
    covs_alt.push_back(cov_alt_1);
    covs_alt.push_back(cov_alt_2);
    covs_alt.push_back(cov_alt_3);
    covs_alt.push_back(cov_alt_4);
    covs_alt.push_back(cov_alt_5);
    covs_alt.push_back(cov_alt_6);
    covs_alt.push_back(cov_alt_7);
    covs_alt.push_back(cov_alt_8);
    covs_alt.push_back(cov_alt_9);
    covs_alt.push_back(cov_alt_10);
    covs_alt.push_back(cov_alt_11);
    covs_alt.push_back(cov_alt_12);
    covs_alt.push_back(cov_alt_13);
    covs_alt.push_back(cov_alt_14);
    covs_alt.push_back(cov_alt_15);
    covs_alt.push_back(cov_alt_16);
 
    cout << "Max cov " << endl;
    TMatrixD cov_alt_max(p_alt_1.combineCov_maxVal(covs_alt));
//     cov_alt_max.Print();
    //covs.push_back(new TMatrixD(cov_alt_max));

    // Take sample variance as systematic 
    vector< vector <double> > vec_vals_alt;
    vec_vals_alt.push_back(vals_alt_1);
    vec_vals_alt.push_back(vals_alt_2);
    vec_vals_alt.push_back(vals_alt_3);
    vec_vals_alt.push_back(vals_alt_4);
    vec_vals_alt.push_back(vals_alt_5);
    vec_vals_alt.push_back(vals_alt_6);
    vec_vals_alt.push_back(vals_alt_7);
    vec_vals_alt.push_back(vals_alt_8);
    vec_vals_alt.push_back(vals_alt_9);
    vec_vals_alt.push_back(vals_alt_10);
    vec_vals_alt.push_back(vals_alt_11);
    vec_vals_alt.push_back(vals_alt_12);
    vec_vals_alt.push_back(vals_alt_13);
    vec_vals_alt.push_back(vals_alt_14);
    vec_vals_alt.push_back(vals_alt_15);
    vec_vals_alt.push_back(vals_alt_16);

    cout << "Sample variance " << endl;
    TMatrixD cov_alt(p_alt_1.sampleVariance(vec_vals_alt));
//     cov_alt.Print();
//     covs.push_back(new TMatrixD(cov_alt));

    /// Total systematics table   
    vector<string> sysNames;
//     sysNames.push_back("Fit bias");
//     sysNames.push_back("Time-Acc.");
    sysNames.push_back("Resolution");
    sysNames.push_back("$\\Delta m_{s}$");
    sysNames.push_back("Asymmetries");
    sysNames.push_back("Background");

    sysNames.push_back("Lineshapes");
    sysNames.push_back("Resonances $m, \\Gamma$");
    sysNames.push_back("Form-Factors");
    sysNames.push_back("Phsp-Acc.");
    sysNames.push_back("Amp. Model");
 
    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal/sys_summary_table.tex",std::ofstream::trunc);

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
        SummaryFile << std::fixed << std::setprecision(2) << p_data.latexName(paraNames[i])  << " & " ;
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
    SummaryFile2.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal/sys_summary_table2.tex",std::ofstream::trunc);

    SummaryFile2 << "\\begin{tabular}{l " ;
    for(int i =0 ; i <= covs.size() ; i++) SummaryFile2 << " c " ;
    SummaryFile2 << " | c }" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "Fit Parameter & " ;
    for(int i =0 ; i <= covs.size() ; i++)  SummaryFile2 << sysNames[i] << " & " ;
    SummaryFile2 << " Total " << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";
    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        SummaryFile2 << std::fixed << std::setprecision(2) << p_data.latexName(paraNames[i])  << " & " ;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
            SummaryFile2 << sqrt((*covs[j])[i][i])/errs_stat[i] << " & ";  
        }
	if(i < 20) SummaryFile2 << " & ";
	else { 
	    	tot += cov_alt[i][i];
		SummaryFile2 << sqrt(cov_alt[i][i])/errs_stat[i] << " & " ;
	} 
        SummaryFile2 << sqrt(tot)/errs_stat[i] << " \\\\ " << "\n"; 
    }
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";


    /// Result table   
    ofstream SummaryFile3;
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal/result_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{l c c c c } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Decay Channel} & \\multicolumn{2}{c}{$A_{b \\to c}$} & \\multicolumn{2}{c}{$A_{b \\to u}$} " << " \\\\ " << "\n";
    SummaryFile3 << " & \\multicolumn{1}{c}{$\\vert a_i \\vert$}  & \\multicolumn{1}{c}{$arg(a_i) [\\degrees]$}  & \\multicolumn{1}{c}{$\\vert a_i \\vert$} & \\multicolumn{1}{c}{$arg(a_i) [\\degrees]$}"  << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";

    // K1(1270)
    SummaryFile3 << " $B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) ) $"  << " & ";
    SummaryFile3 << " 1.0 & 0.0 & 1.0 & 0.0 " << " \\\\ " << "\n";
    SummaryFile3 << p_data.latexNameMod(paraNames[0])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[0] << " $\\pm$ " << errs_stat[0]  << " $\\pm$ " << errs_sys[0] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[1] << " $\\pm$ " << errs_stat[1]  << " $\\pm$ " << errs_sys[1] << " & &  " << " \\\\ " << "\n";

    SummaryFile3 << p_data.latexNameMod(paraNames[2])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[2] << " $\\pm$ " << errs_stat[2]  << " $\\pm$ " << errs_sys[2] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[3] << " $\\pm$ " << errs_stat[3]  << " $\\pm$ " << errs_sys[3] << " & &  " << " \\\\ " << "\n";

    // K1(1400)
    SummaryFile3 << p_data.latexNameMod(paraNames[4])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[4] << " $\\pm$ " << errs_stat[4]  << " $\\pm$ " << errs_sys[4] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[5] << " $\\pm$ " << errs_stat[5]  << " $\\pm$ " << errs_sys[5] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[6] << " $\\pm$ " << errs_stat[6]  << " $\\pm$ " << errs_sys[6] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[7] << " $\\pm$ " << errs_stat[7]  << " $\\pm$ " << errs_sys[7] << " \\\\ " << "\n";

    // K1s
    SummaryFile3 << p_data.latexNameMod(paraNames[8])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[8] << " $\\pm$ " << errs_stat[8]  << " $\\pm$ " << errs_sys[8] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[9] << " $\\pm$ " << errs_stat[9]  << " $\\pm$ " << errs_sys[9] << " & ";
    SummaryFile3 << " & " << " \\\\ " << "\n";

    SummaryFile3 << p_data.latexNameMod(paraNames[10])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[10] << " $\\pm$ " << errs_stat[10]  << " $\\pm$ " << errs_sys[10] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[11] << " $\\pm$ " << errs_stat[11]  << " $\\pm$ " << errs_sys[11] << " & &  " << " \\\\ " << "\n";

    // K(1460)
    SummaryFile3 << p_data.latexNameMod(paraNames[12])  << " & & &";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[12] << " $\\pm$ " << errs_stat[12]  << " $\\pm$ " << errs_sys[12] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[13] << " $\\pm$ " << errs_stat[13]  << " $\\pm$ " << errs_sys[13] << " \\\\ " << "\n";

    // NV Ks
    SummaryFile3 << p_data.latexNameMod(paraNames[14])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[14] << " $\\pm$ " << errs_stat[14]  << " $\\pm$ " << errs_sys[14] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);  
    SummaryFile3 << vals[15] << " $\\pm$ " << errs_stat[15]  << " $\\pm$ " << errs_sys[15] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[16] << " $\\pm$ " << errs_stat[16]  << " $\\pm$ " << errs_sys[16] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[17] << " $\\pm$ " << errs_stat[17]  << " $\\pm$ " << errs_sys[17] << " \\\\ " << "\n";

    // NV rho
    SummaryFile3 << p_data.latexNameMod(paraNames[18])  << " & & &";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[18] << " $\\pm$ " << errs_stat[18]  << " $\\pm$ " << errs_sys[18] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[19] << " $\\pm$ " << errs_stat[19]  << " $\\pm$ " << errs_sys[19] << " \\\\ " << "\n";

    // masses,widths
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Fit parameter} & \\multicolumn{4}{c}{Value} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[20]) << "} & \\multicolumn{4}{c}{" << vals[20] << " $\\pm$ " << errs_stat[20] << " $\\pm$ " << errs_sys[20] << " $\\pm$ " << sqrt(cov_alt[20][20]) << "} \\\\ " << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[21]) << "} & \\multicolumn{4}{c}{" << vals[21] << " $\\pm$ " << errs_stat[21] << " $\\pm$ " << errs_sys[21] << " $\\pm$ " << sqrt(cov_alt[21][21]) << "} \\\\ " << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[22]) << "} & \\multicolumn{4}{c}{" << vals[22] << " $\\pm$ " << errs_stat[22] << " $\\pm$ " << errs_sys[22] << " $\\pm$ " << sqrt(cov_alt[22][22]) << "} \\\\ " << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[23]) << "} & \\multicolumn{4}{c}{" << vals[23] << " $\\pm$ " << errs_stat[23] << " $\\pm$ " << errs_sys[23] << " $\\pm$ " << sqrt(cov_alt[23][23]) << "} \\\\ " << "\n";
    SummaryFile3 << " \\\\ " << "\n";

    // r, delta, gamma
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[24]) << "} & \\multicolumn{4}{c}{" << "xx.xx" << " $\\pm$ " << errs_stat[24] << " $\\pm$ " << errs_sys[24] << " $\\pm$ " << sqrt(cov_alt[24][24]) << "} \\\\ " << "\n";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[25]) << "} & \\multicolumn{4}{c}{" << "xx.xx" << " $\\pm$ " << errs_stat[25] << " $\\pm$ " << errs_sys[25] << " $\\pm$ " << sqrt(cov_alt[25][25]) << "} \\\\ " << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[26]) << "} & \\multicolumn{4}{c}{" << "xx.xx" << " $\\pm$ " << errs_stat[26] << " $\\pm$ " << errs_sys[26] << " $\\pm$ " << sqrt(cov_alt[26][26]) << "} \\\\ " << "\n";
    
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";


}

void fractions(){

   /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Sum");
 
    vector<TString> paraNames_bar;
    paraNames_bar.push_back("Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames_bar.push_back("Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");
    paraNames_bar.push_back("Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");
    paraNames_bar.push_back("Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames_bar.push_back("Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames_bar.push_back("Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames_bar.push_back("Bs0_NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames_bar.push_back("Sum");

    pull p(paraNames,"signal_toy/fitFractions_*_Bar.root","fractions", true, true,100);
    vector<double> vals = p.sampleMean()  ;
    vector<double> errs_stat = p.sampleSigma()  ;

    pull p_bar(paraNames_bar,"signal_toy/fitFractions_*_Bar.root","fractions", true, false,100);
    vector<double> vals_bar = p_bar.sampleMean()  ;
    vector<double> errs_stat_bar = p_bar.sampleSigma()  ;


    /// Result table   
    ofstream SummaryFile3;
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal/fraction_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{l r r } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Decay Channel} & \\multicolumn{1}{c}{$F_{b \\to c} [\\%]$} & \\multicolumn{1}{c}{$F_{b \\to u} [\\%]$} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";

    SummaryFile3 << std::fixed << std::setprecision(1);
    // K1(1270)
    SummaryFile3 << p.latexNameMod(paraNames[0]) ; 
    SummaryFile3 << " & "  << vals[0] * 100. << " $\\pm$ " << errs_stat[0]* 100. ;
    SummaryFile3 << " & "  << vals_bar[0]* 100. << " $\\pm$ " << errs_stat_bar[0]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[1]) ; 
    SummaryFile3 << " & "  << vals[1] * 100. << " $\\pm$ " << errs_stat[1]* 100. ;
    SummaryFile3 << " & "  << vals_bar[1]* 100. << " $\\pm$ " << errs_stat_bar[1]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[2]) ; 
    SummaryFile3 << " & "  << vals[2] * 100. << " $\\pm$ " << errs_stat[2]* 100. ;
    SummaryFile3 << " & "  << vals_bar[2]* 100. << " $\\pm$ " << errs_stat_bar[2]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    // K1(1400)
    SummaryFile3 << p.latexNameMod(paraNames[3]) ; 
    SummaryFile3 << " & "  << vals[3] * 100. << " $\\pm$ " << errs_stat[3]* 100. ;
    SummaryFile3 << " & "  << vals_bar[3]* 100. << " $\\pm$ " << errs_stat_bar[3]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    // Ks(1410)
    SummaryFile3 << p.latexNameMod(paraNames[5]) ; 
    SummaryFile3 << " & "  << vals[5] * 100. << " $\\pm$ " << errs_stat[5]* 100. ;
    SummaryFile3 << " & "  ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[6]) ; 
    SummaryFile3 << " & "  << vals[6] * 100. << " $\\pm$ " << errs_stat[6]* 100. ;
    SummaryFile3 << " & "  ;
    SummaryFile3 << " \\\\ " << "\n";

    // K(1460)
    SummaryFile3 << p.latexNameMod(paraNames_bar[5]) ; 
    SummaryFile3 << " & "  ;
    SummaryFile3 << " & "  << vals_bar[5]* 100. << " $\\pm$ " << errs_stat_bar[5]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    // NS
    SummaryFile3 << p.latexNameMod(paraNames[4]) ; 
    SummaryFile3 << " & "  << vals[4] * 100. << " $\\pm$ " << errs_stat[4]* 100. ;
    SummaryFile3 << " & "  << vals_bar[4]* 100. << " $\\pm$ " << errs_stat_bar[4]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames_bar[6]) ; 
    SummaryFile3 << " & "   ;
    SummaryFile3 << " & "  << vals_bar[6]* 100. << " $\\pm$ " << errs_stat_bar[6]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << p.latexNameMod(paraNames[7]) ; 
    SummaryFile3 << " & "  << vals[7] * 100. << " $\\pm$ " << errs_stat[7]* 100. ;
    SummaryFile3 << " & "  << vals_bar[7]* 100. << " $\\pm$ " << errs_stat_bar[7]* 100. ;
    SummaryFile3 << " \\\\ " << "\n";
    
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

    fitParams();
//     fractions();

    return 0;
}
