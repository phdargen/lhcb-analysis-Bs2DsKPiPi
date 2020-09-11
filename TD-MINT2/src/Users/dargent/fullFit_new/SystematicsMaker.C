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

    paraNames.push_back("a_K1_1400_Amp");
    paraNames.push_back("a_K1_1400_Phase");
    paraNames.push_back("abar_K1_1400_Amp");
    paraNames.push_back("abar_K1_1400_Phase");
    paraNames.push_back("Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Amp");
    paraNames.push_back("Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Phase");
    paraNames.push_back("Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Amp");
    paraNames.push_back("Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Phase");

    paraNames.push_back("a_Ks_1410_Amp");
    paraNames.push_back("a_Ks_1410_Phase");
    paraNames.push_back("abar_Ks_1410_Amp");
    paraNames.push_back("abar_Ks_1410_Phase");
    paraNames.push_back("Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Amp");
    paraNames.push_back("Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Phase");

    paraNames.push_back("abar_K_1460_Amp");
    paraNames.push_back("abar_K_1460_Phase");

    paraNames.push_back("a_NS_Ks_Amp");
    paraNames.push_back("a_NS_Ks_Phase");
    paraNames.push_back("abar_NS_Ks_Amp");
    paraNames.push_back("abar_NS_Ks_Phase");

    paraNames.push_back("a_NS_rho_Amp");
    paraNames.push_back("a_NS_rho_Phase");

    paraNames.push_back("mass_K_1__1400_p");
    paraNames.push_back("width_K_1__1400_p");
    paraNames.push_back("mass_Ks_1410_p");
    paraNames.push_back("width_Ks_1410_p");

    paraNames.push_back("r");
    paraNames.push_back("delta");
    paraNames.push_back("gamma");
 

    vector<TMatrixD*> covs;
    /// Data fit
    //pull p_data(paraNames,"out/pull__1.root");  
    pull p_data(paraNames,"signal_new2/pull__2.root");   
    pull p_data_err(paraNames,"out/pull__2.root");   
 
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p_data_err.getErrs();
    errs_stat[paraNames.size()-1] = 10.4;
    errs_stat[paraNames.size()-2] = 9.9;   

    /// Fit bias from toys
    pull p(paraNames,"out_signal_new_toy/pull__*.root");
    TMatrixD* cov = new TMatrixD(p.getCov("",4,true));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov)[i][j] = (*cov)[i][j] * errs_stat[i] * errs_stat[j];
    }
    covs.push_back(cov);

    /// Acc systematics
    pull p_acc(paraNames,"out_signal_new_toy/pullAcc_*.root");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("out_signal_new_toy/pull__*.root","_acc"));

    
    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"out_sys_res/pull__1.root");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("out/pull__1.root","_res_Run1_a"));

    pull p_res_Run1_b(paraNames,"out_sys_res/pull__2.root");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("out/pull__1.root","_res_Run1_b"));

    pull p_res_Run2_a(paraNames,"out_sys_res/pull__3.root");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("out/pull__1.root","_res_Run2_a"));

    pull p_res_Run2_b(paraNames,"out_sys_res/pull__4.root");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("out/pull__1.root","_res_Run2_b"));

    pull p_res_Run2_c(paraNames,"out_sys_res/pull__5.root");
    TMatrixD* cov_res_Run2_c = new TMatrixD(p_res_Run2_c.getDeltaCov("out/pull__1.root","_res_Run2_c"));

    pull p_res_Run2_d(paraNames,"out_sys_res/pull__6.root");
    TMatrixD* cov_res_Run2_d = new TMatrixD(p_res_Run2_d.getDeltaCov("out/pull__1.root","_res_Run2_d"));

    pull p_res_Run2_e(paraNames,"out_sys_res/pull__7.root");
    TMatrixD* cov_res_Run2_e = new TMatrixD(p_res_Run2_e.getDeltaCov("out/pull__1.root","_res_Run2_e"));

    pull p_res_Run2_f(paraNames,"out_sys_res/pull__8.root");
    TMatrixD* cov_res_Run2_f = new TMatrixD(p_res_Run2_f.getDeltaCov("out/pull__1.root","_res_Run2_f"));

    vector<TMatrixD*> covs_res_Run1;
    covs_res_Run1.push_back(cov_res_Run1_a);
    covs_res_Run1.push_back(cov_res_Run1_b);
    TMatrixD cov_res(p_res_Run1_a.combineCov_maxVal(covs_res_Run1));

    vector<TMatrixD*> covs_res_Run2;
    covs_res_Run2.push_back(cov_res_Run2_a);
    covs_res_Run2.push_back(cov_res_Run2_b);
    covs_res_Run2.push_back(cov_res_Run2_c);
    covs_res_Run2.push_back(cov_res_Run2_d);
    //covs_res_Run2.push_back(cov_res_Run2_e);
    //covs_res_Run2.push_back(cov_res_Run2_f);
    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_Run2) ;

    vector<TMatrixD*> covs_res_18;
    covs_res_18.push_back(cov_res_Run2_e);
    covs_res_18.push_back(cov_res_Run2_f);
    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_18) ;

    /// decay time bias systematics 
    pull p_bias(paraNames,"out_sys_bias/pull__1.root");
    TMatrixD* cov_bias = new TMatrixD(p_bias.getDeltaCov("out/pull__1.root","_bias"));

    /// dms systematics 
    pull p_dm(paraNames,"out_signal_new_toy/pull_dm_*.root");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("out_signal_new_toy/pull__*.root","_dm",0.5));

    /// asymmetry systematics
    pull p_asym(paraNames,"out_signal_new_toy/pull_asym_*.root");
    TMatrixD* cov_asym = new TMatrixD(p_asym.getDeltaCov("out_signal_new_toy/pull__*.root","_asym",0.75));

    /// bkg systematics 
    pull p_bkg_1(paraNames,"out_sys_bkg/pull__1.root");
    TMatrixD* cov_bkg_1 = new TMatrixD(p_bkg_1.getDeltaCov("out_sys_bkg/pull__4.root","bkg_1"));
    vector<double> vals_bkg_1 = p_bkg_1.getVals();

    pull p_bkg_2(paraNames,"out_sys_bkg/pull__2.root");
    TMatrixD* cov_bkg_2 = new TMatrixD(p_bkg_2.getDeltaCov("out_sys_bkg/pull__4.root","bkg_2"));
    vector<double> vals_bkg_2 = p_bkg_2.getVals();

    pull p_bkg_3(paraNames,"out_sys_bkg/pull__3.root");
    TMatrixD* cov_bkg_3 = new TMatrixD(p_bkg_3.getDeltaCov("out_sys_bkg/pull__4.root","bkg_3"));
    vector<double> vals_bkg_3 = p_bkg_3.getVals();

    pull p_bkg_4(paraNames,"out_sys_bkg/pull__4.root");
    TMatrixD* cov_bkg_4 = new TMatrixD(p_bkg_4.getDeltaCov("out_sys_bkg/pull__4.root","bkg_4"));
    vector<double> vals_bkg_4 = p_bkg_4.getVals();

    pull p_bkg_5(paraNames,"out_sys_bkg/pull__5.root");
    TMatrixD* cov_bkg_5 = new TMatrixD(p_bkg_5.getDeltaCov("out_sys_bkg/pull__4.root","bkg_5"));
    vector<double> vals_bkg_5 = p_bkg_5.getVals();

    pull p_bkg_6(paraNames,"out_sys_bkg/pull__6.root");
    TMatrixD* cov_bkg_6 = new TMatrixD(p_bkg_6.getDeltaCov("out_sys_bkg/pull__4.root","bkg_6"));
    vector<double> vals_bkg_6 = p_bkg_6.getVals();

    pull p_bkg_7(paraNames,"out_sys_bkg/pull__7.root");
    TMatrixD* cov_bkg_7 = new TMatrixD(p_bkg_7.getDeltaCov("out_sys_bkg/pull__4.root","bkg_7"));
    vector<double> vals_bkg_7 = p_bkg_7.getVals();

    pull p_bkg_8(paraNames,"out_sys_bkg/pull__8.root");
    TMatrixD* cov_bkg_8 = new TMatrixD(p_bkg_8.getDeltaCov("out_sys_bkg/pull__4.root","bkg_8"));
    vector<double> vals_bkg_8 = p_bkg_8.getVals();

    pull p_bkg_9(paraNames,"out_sys_bkg/pull__9.root");
    TMatrixD* cov_bkg_9 = new TMatrixD(p_bkg_9.getDeltaCov("out_sys_bkg/pull__4.root","bkg_9"));
    vector<double> vals_bkg_9 = p_bkg_9.getVals();

    pull p_bkg_10(paraNames,"out_sys_bkg/pull__10.root");
    TMatrixD* cov_bkg_10 = new TMatrixD(p_bkg_10.getDeltaCov("out_sys_bkg/pull__4.root","bkg_10"));
    vector<double> vals_bkg_10 = p_bkg_10.getVals();

    // Take maximum as systematic
    vector<TMatrixD*> covs_bkg;
    covs_bkg.push_back(cov_bkg_1);
    covs_bkg.push_back(cov_bkg_2);
    covs_bkg.push_back(cov_bkg_3);
    covs_bkg.push_back(cov_bkg_4);
    covs_bkg.push_back(cov_bkg_5);
    covs_bkg.push_back(cov_bkg_6);
    covs_bkg.push_back(cov_bkg_7);
    covs_bkg.push_back(cov_bkg_8);
    covs_bkg.push_back(cov_bkg_9);
    covs_bkg.push_back(cov_bkg_10);

    TMatrixD cov_bkg_max(p_bkg_2.combineCov_maxVal(covs_bkg));
    //cov_bkg_max.Print();
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
    vec_vals_bkg.push_back(vals_bkg_8);
    //vec_vals_bkg.push_back(vals_bkg_9);
    //vec_vals_bkg.push_back(vals_bkg_10);

    TMatrixD cov_bkg(p_bkg_2.sampleVariance(vec_vals_bkg,1));
//     cov_bkg.Print();

    cout << endl << vals_bkg_1[3] << endl;
    cout << vals_bkg_2[3] << endl;
    cout << vals_bkg_3[3] << endl;
    cout << vals_bkg_4[3] << endl;
    cout << vals_bkg_5[3] << endl;
    cout << vals_bkg_6[3] << endl;
    cout << vals_bkg_7[3] << endl;
    cout << vals_bkg_8[3] << endl;
    cout << vals_bkg_9[3] << endl;
    cout << vals_bkg_10[3] << endl << endl;

    pull p_bkg_bs1(paraNames,"out_signal_new_toy/sw_pull__*.root");
    vector<double> mean_bkg_bs1 = p_bkg_bs1.getPullMean();
    vector<double> sigma_bkg_bs1 = p_bkg_bs1.getPullSigma();

    pull p_bkg_bs2(paraNames,"out_signal_new_toy/pull__*.root");
    vector<double> mean_bkg_bs2 = p_bkg_bs2.getPullMean();
    vector<double> sigma_bkg_bs2 = p_bkg_bs2.getPullSigma();

    TMatrixD* cov_bkg_bs = new TMatrixD(cov_bkg);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++){
	if(i==j){
 		//if(i<23)
		(*cov_bkg_bs)[i][j] = max(pow(min(mean_bkg_bs1[i]-mean_bkg_bs2[i],mean_bkg_bs1[i]),2) , pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2))*errs_stat[i] * errs_stat[j];
		//else (*cov_bkg_bs)[i][j] = max(sqrt((pow(min(mean_bkg_bs1[i]-mean_bkg_bs2[i],mean_bkg_bs1[i]),2))) + sqrt((pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2))) )*errs_stat[i] * errs_stat[j];

 		//if(i>22)(*cov_bkg_bs)[i][j] = max(pow(min(mean_bkg_bs1[i]-mean_bkg_bs2[i],mean_bkg_bs1[i]),2) , pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2))*errs_stat[i] * errs_stat[j];


/*
		cout << paraNames[i] << endl;
		cout << mean_bkg_bs1[i] << endl;
		cout << mean_bkg_bs2[i] << endl;
		cout << mean_bkg_bs1[i]-mean_bkg_bs2[i] << endl;
		cout << sigma_bkg_bs1[i]  << endl;
		cout << sigma_bkg_bs2[i]  << endl;
		cout << sigma_bkg_bs1[i]-sigma_bkg_bs2[i] << endl << endl;
		cout << (sqrt((pow(min(mean_bkg_bs1[i]-mean_bkg_bs2[i],mean_bkg_bs1[i]),2))) + sqrt((pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2))) )/2. << endl << endl;
*/
	}
	else (*cov_bkg_bs)[i][j] = 0.;
     }
     cov_bkg+= *cov_bkg_bs;


    /// Lineshape models systematics
    pull p_ls_1(paraNames,"out_sys_ls/pull__1.root");
    TMatrixD* cov_ls_1 = new TMatrixD(p_ls_1.getDeltaCov("out/pull__1.root","ls_1"));
    vector<double> vals_ls_1 = p_ls_1.getVals();
    cout << "ls 1" << endl;
    cout << vals_ls_1[3] << endl << endl;

    pull p_ls_2(paraNames,"out_sys_ls/pull__2.root");
    TMatrixD* cov_ls_2 = new TMatrixD(p_ls_2.getDeltaCov("out/pull__1.root","ls_2"));
    vector<double> vals_ls_2 = p_ls_2.getVals();
    cout << "ls 2" << endl;
    cout << vals_ls_2[3] << endl << endl;

    pull p_ls_3(paraNames,"out_sys_ls/pull__3.root");
    TMatrixD* cov_ls_3 = new TMatrixD(p_ls_3.getDeltaCov("out/pull__1.root","ls_3"));
    vector<double> vals_ls_3 = p_ls_3.getVals();
    cout << "ls 3" << endl;
    cout << vals_ls_3[3] << endl << endl;

    pull p_ls_4(paraNames,"out_sys_ls/pull__4.root");
    TMatrixD* cov_ls_4 = new TMatrixD(p_ls_4.getDeltaCov("out/pull__1.root","ls_4"));
    vector<double> vals_ls_4 = p_ls_4.getVals();
    cout << "ls 4" << endl;
    cout << vals_ls_4[3] << endl << endl;

    pull p_ls_5(paraNames,"out_sys_ls/pull__5.root");
    TMatrixD* cov_ls_5 = new TMatrixD(p_ls_5.getDeltaCov("out/pull__1.root","ls_5"));
    vector<double> vals_ls_5 = p_ls_5.getVals();
    cout << "ls 5" << endl;
    cout << vals_ls_5[3] << endl << endl;

    pull p_ls_6(paraNames,"out_sys_ls/pull__6.root");
    TMatrixD* cov_ls_6 = new TMatrixD(p_ls_6.getDeltaCov("out/pull__1.root","ls_6"));
    vector<double> vals_ls_6 = p_ls_6.getVals();
    cout << "ls 6" << endl;
    cout << vals_ls_6[3] << endl << endl;

    pull p_ls_7(paraNames,"out_sys_ls/pull__7.root");
    TMatrixD* cov_ls_7 = new TMatrixD(p_ls_7.getDeltaCov("out/pull__1.root","ls_7"));
    vector<double> vals_ls_7 = p_ls_7.getVals();
    cout << "ls 7" << endl;
    cout << vals_ls_7[3] << endl << endl;

    vector< vector <double> > vec_vals_ls;
    vec_vals_ls.push_back(vals_ls_1);
    vec_vals_ls.push_back(vals_ls_2);
    vec_vals_ls.push_back(vals_ls_3);
    vec_vals_ls.push_back(vals_ls_4);
    vec_vals_ls.push_back(vals_ls_5);
    vec_vals_ls.push_back(vals_ls_6);
    vec_vals_ls.push_back(vals_ls_7);
    TMatrixD cov_ls(p_ls_1.sampleVariance(vec_vals_ls,22));


    /// Resonance parameters systematics
    pull p_rp_1(paraNames,"out_signal_new_toy/pull_ls1_*.root");
    TMatrixD* cov_rp_1 = new TMatrixD(p_rp_1.getDeltaCov("out_signal_new_toy/pull__*.root","_ls1",0.5));

    pull p_rp_2(paraNames,"out_signal_new_toy/pull_ls2_*.root");
    TMatrixD* cov_rp_2 = new TMatrixD(p_rp_2.getDeltaCov("out_signal_new_toy/pull__*.root","_ls2",0.5));

    cout << "r 1 = " << sqrt((*cov_rp_1)[6][6]) << endl;
    cout << "r 2 = " << sqrt((*cov_rp_2)[6][6]) << endl;

    // Add
    TMatrixD cov_rp(*cov_rp_1);
    cov_rp +=  *cov_rp_2 ;
 
    /// Form factor
    pull p_f_1(paraNames,"out_signal_new_toy/pull_BW_radius_*.root");
    TMatrixD* cov_f_1 = new TMatrixD(p_f_1.getDeltaCov("out_signal_new_toy/pull__*.root",2));


    /// Phsp-Acc systematics
    pull p_phsp_acc_1(paraNames,"out_sys_acc/pull__1.root");
    TMatrixD* cov_phsp_acc_1 = new TMatrixD(p_phsp_acc_1.getDeltaCov("out/pull__1.root","phsp_acc_1"));
    vector<double> vals_ps_1 = p_phsp_acc_1.getVals();

    pull p_phsp_acc_2(paraNames,"out_sys_acc/pull__2.root");
    TMatrixD* cov_phsp_acc_2 = new TMatrixD(p_phsp_acc_2.getDeltaCov("out/pull__1.root","phsp_acc_2"));
    vector<double> vals_ps_2 = p_phsp_acc_2.getVals();

    pull p_phsp_acc_3(paraNames,"out_sys_acc/pull__3.root");
    TMatrixD* cov_phsp_acc_3 = new TMatrixD(p_phsp_acc_3.getDeltaCov("out/pull__1.root","phsp_acc_3"));
    vector<double> vals_ps_3 = p_phsp_acc_3.getVals();

    pull p_phsp_acc_4(paraNames,"out_sys_acc/pull__4.root");
    TMatrixD* cov_phsp_acc_4 = new TMatrixD(p_phsp_acc_4.getDeltaCov("out/pull__1.root","phsp_acc_4"));
    vector<double> vals_ps_4 = p_phsp_acc_4.getVals();

    pull p_phsp_acc_5(paraNames,"out_sys_acc/pull__5.root");
    TMatrixD* cov_phsp_acc_5 = new TMatrixD(p_phsp_acc_5.getDeltaCov("out/pull__1.root","phsp_acc_5"));
    vector<double> vals_ps_5 = p_phsp_acc_5.getVals();

    pull p_phsp_acc_6(paraNames,"out_sys_acc/pull__6.root");
    TMatrixD* cov_phsp_acc_6 = new TMatrixD(p_phsp_acc_6.getDeltaCov("out/pull__1.root","phsp_acc_6"));
    vector<double> vals_ps_6 = p_phsp_acc_6.getVals();

    pull p_phsp_acc_7(paraNames,"out_sys_acc/pull__7.root");
    TMatrixD* cov_phsp_acc_7 = new TMatrixD(p_phsp_acc_7.getDeltaCov("out/pull__1.root","phsp_acc_7"));
    vector<double> vals_ps_7 = p_phsp_acc_7.getVals();

    pull p_phsp_acc_8(paraNames,"out_sys_acc/pull__8.root");
    TMatrixD* cov_phsp_acc_8 = new TMatrixD(p_phsp_acc_8.getDeltaCov("out/pull__1.root","phsp_acc_8"));
    vector<double> vals_ps_8 = p_phsp_acc_8.getVals();

    pull p_phsp_acc_9(paraNames,"out_sys_acc/pull__9.root");
    TMatrixD* cov_phsp_acc_9 = new TMatrixD(p_phsp_acc_9.getDeltaCov("out/pull__1.root","phsp_acc_9"));
    vector<double> vals_ps_9 = p_phsp_acc_9.getVals();


    vector< vector <double> > vec_vals_ps;
    vec_vals_ps.push_back(vals_ps_1);
    vec_vals_ps.push_back(vals_ps_2);
    vec_vals_ps.push_back(vals_ps_3);
    vec_vals_ps.push_back(vals_ps_4);
    vec_vals_ps.push_back(vals_ps_5);
    vec_vals_ps.push_back(vals_ps_6);
    vec_vals_ps.push_back(vals_ps_7);
    vec_vals_ps.push_back(vals_ps_8);
    vec_vals_ps.push_back(vals_ps_9);

    TMatrixD cov_phsp(p_phsp_acc_1.sampleVariance(vec_vals_ps));

    /// Phsp acc factorization
    pull p_phsp_acc2_0(paraNames,"out_sys_acc2/pull__0.root");
    vector<double> vals_ps_acc0 = p_phsp_acc2_0.getVals();

    pull p_phsp_acc2_1(paraNames,"out_sys_acc2/pull__1.root");
    TMatrixD* cov_phsp_acc2_1 = new TMatrixD(p_phsp_acc2_1.getDeltaCov("out_sys_acc2/pull__0.root","phsp_acc2_1"));
    vector<double> vals_ps_acc1 = p_phsp_acc2_1.getVals();

    pull p_phsp_acc2_2(paraNames,"out_sys_acc2/pull__2.root");
    TMatrixD* cov_phsp_acc2_2 = new TMatrixD(p_phsp_acc2_2.getDeltaCov("out_sys_acc2/pull__0.root","phsp_acc2_2"));
    vector<double> vals_ps_acc2 = p_phsp_acc2_2.getVals();

    pull p_phsp_acc2_3(paraNames,"out_sys_acc2/pull__3.root");
    TMatrixD* cov_phsp_acc2_3 = new TMatrixD(p_phsp_acc2_3.getDeltaCov("out_sys_acc2/pull__0.root","phsp_acc2_3"));
    vector<double> vals_ps_acc3 = p_phsp_acc2_3.getVals();
	
    vector< vector <double> > vec_vals_ps2;
    vec_vals_ps2.push_back(vals_ps_acc0);
    vec_vals_ps2.push_back(vals_ps_acc1);
    vec_vals_ps2.push_back(vals_ps_acc2);
    vec_vals_ps2.push_back(vals_ps_acc3);

    TMatrixD cov_phsp2(p_phsp_acc2_1.sampleVariance(vec_vals_ps2,22));


    /// Alternative amp models 

    pull p_alt_1(paraNames,"out_alt/pull__1.root");
    TMatrixD* cov_alt_1 = new TMatrixD(p_alt_1.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_1 = p_alt_1.getVals();

    pull p_alt_2(paraNames,"out_alt/pull__9992.root");
    TMatrixD* cov_alt_2 = new TMatrixD(p_alt_2.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_2 = p_alt_2.getVals();

    pull p_alt_3(paraNames,"out_alt/pull__93.root");
    TMatrixD* cov_alt_3 = new TMatrixD(p_alt_3.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_3 = p_alt_3.getVals();

    pull p_alt_4(paraNames,"out_alt/pull__94.root");
    TMatrixD* cov_alt_4 = new TMatrixD(p_alt_4.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_4 = p_alt_4.getVals();

    pull p_alt_5(paraNames,"out_alt/pull__5.root");
    TMatrixD* cov_alt_5 = new TMatrixD(p_alt_5.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_5 = p_alt_5.getVals();

    pull p_alt_6(paraNames,"out_alt/pull__6.root");
    TMatrixD* cov_alt_6 = new TMatrixD(p_alt_6.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_6 = p_alt_6.getVals();

    pull p_alt_7(paraNames,"out_alt/pull__9997.root");
    TMatrixD* cov_alt_7 = new TMatrixD(p_alt_7.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_7 = p_alt_7.getVals();

    pull p_alt_8(paraNames,"out_alt/pull__9998.root");
    TMatrixD* cov_alt_8 = new TMatrixD(p_alt_8.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_8 = p_alt_8.getVals();

    pull p_alt_9(paraNames,"out_alt/pull__99.root");
    TMatrixD* cov_alt_9 = new TMatrixD(p_alt_9.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_9 = p_alt_9.getVals();

    pull p_alt_10(paraNames,"out_alt/pull__999997.root");
    TMatrixD* cov_alt_10 = new TMatrixD(p_alt_10.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_10 = p_alt_10.getVals();

    pull p_alt_11(paraNames,"out_alt/pull__14.root");
    TMatrixD* cov_alt_11 = new TMatrixD(p_alt_11.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_11 = p_alt_11.getVals();

    pull p_alt_12(paraNames,"out_alt/pull__912.root");
    TMatrixD* cov_alt_12 = new TMatrixD(p_alt_12.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_12 = p_alt_12.getVals();

    pull p_alt_13(paraNames,"out_alt13/pull__7.root");
    TMatrixD* cov_alt_13 = new TMatrixD(p_alt_13.getDeltaCov("out/pull__1.root"));
    vector<double> vals_alt_13 = p_alt_13.getVals();

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
    //vec_vals_alt.push_back(vals_alt_13);

    cout << "Sample variance " << endl;
    TMatrixD cov_alt(p_alt_1.sampleVariance(vec_vals_alt));
//     cov_alt.Print();
         //TMatrixD cov_alt(paraNames.size(),paraNames.size());

	cout << "alt models" << endl;

	cout << vals_alt_1[22] << endl;
	cout << vals_alt_2[22] << endl;
	cout << vals_alt_3[22] << endl;
	cout << vals_alt_4[22] << endl;
	cout << vals_alt_6[22] << endl;
	cout << vals_alt_6[22] << endl;
	cout << vals_alt_7[22] << endl;
	cout << vals_alt_8[22] << endl;
	cout << vals_alt_9[22] << endl;
	cout << vals_alt_10[22] << endl;
	cout << vals_alt_11[22] << endl;
	cout << vals_alt_12[22] << endl;
	cout << vals_alt_13[22] << endl;

	cout << sqrt(cov_alt[22][22]) << endl;



    /// m,t correlations
    pull p_corr1(paraNames,"out_signal_new_toy/sw_pull__*.root");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"out_signal_new_toy_bkg/sw_pull__*.root");
    TMatrixD cov_corr2 = p_corr2.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
    }

//      pull p_corr(paraNames,"signal_toy_bkg_bs0/pull__*.root");
//      TMatrixD* cov_corr = new TMatrixD(p_corr.getDeltaCov("signal_toy_bkg_bs1/pull__*.root","_corr"));


    /// Total systematics table   
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(cov_corr);
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(cov_bias);
    covs.push_back(cov_asym);
    covs.push_back(cov_dm);
    covs.push_back(new TMatrixD(cov_phsp));
    covs.push_back(new TMatrixD(cov_phsp2));
    covs.push_back(new TMatrixD(cov_ls));
    covs.push_back(new TMatrixD(cov_rp));
    covs.push_back(cov_f_1);

    vector<string> sysNames;
    sysNames.push_back("Fit bias");
    sysNames.push_back("Background");
    sysNames.push_back("Correlations");
    sysNames.push_back("Time-Acc.");
    sysNames.push_back("Resolution");
    sysNames.push_back("Decay-time bias");
    sysNames.push_back("Asymmetries");
    sysNames.push_back("$\\Delta m_{s}$");
    sysNames.push_back("Phsp-Acc.");
    sysNames.push_back("Acc. Factor.");
    sysNames.push_back("Lineshapes");
    sysNames.push_back("Resonances $m,\\Gamma$");
    sysNames.push_back("Form-Factors");
    sysNames.push_back("Amp. Model");
 
    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/sys_summary_table.tex",std::ofstream::trunc);

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
    SummaryFile2.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/sys_summary_table2.tex",std::ofstream::trunc);

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
	if(i < 22) SummaryFile2 << " & ";
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
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/result_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{l c c c c } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Decay Channel} & \\multicolumn{2}{c}{$A_{b \\to c}$} & \\multicolumn{2}{c}{$A_{b \\to u}$} " << " \\\\ " << "\n";
    SummaryFile3 << " & \\multicolumn{1}{c}{$\\vert a_i \\vert$}  & \\multicolumn{1}{c}{$arg(a_i) [\\degrees]$}  & \\multicolumn{1}{c}{$\\vert a_i \\vert$} & \\multicolumn{1}{c}{$arg(a_i) [\\degrees]$}"  << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";

	cout << "here 5" << endl;


    // K1(1270)
    SummaryFile3 << p_data.latexNameMod(paraNames[0])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[0] << " $\\pm$ " << errs_stat[0]  << " $\\pm$ " << errs_sys[0] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[1] << " $\\pm$ " << errs_stat[1]  << " $\\pm$ " << errs_sys[1] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[2] << " $\\pm$ " << errs_stat[2]  << " $\\pm$ " << errs_sys[2] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);
    SummaryFile3 << vals[3] << " $\\pm$ " << errs_stat[3]  << " $\\pm$ " << errs_sys[3] << " \\\\ " << "\n";

    SummaryFile3 << p_data.latexNameMod(paraNames[4])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[4] << " $\\pm$ " << errs_stat[4]  << " $\\pm$ " << errs_sys[4] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[5] << " $\\pm$ " << errs_stat[5]  << " $\\pm$ " << errs_sys[5] << " & ";
    SummaryFile3 << " & " << " \\\\ " << "\n";

    SummaryFile3 << p_data.latexNameMod(paraNames[6])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[6] << " $\\pm$ " << errs_stat[6]  << " $\\pm$ " << errs_sys[6] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[7] << " $\\pm$ " << errs_stat[7]  << " $\\pm$ " << errs_sys[7] << " & ";
    SummaryFile3 << " & " << " \\\\ " << "\n";

    // K1(1400)
    SummaryFile3 << " $B_s \\to D_s \\, ( K_1(1400) \\to K^{*}(892) \\, \\pi ) $"  << " & ";
    SummaryFile3 << " 1.0 & 0.0 & 1.0 & 0.0 " << " \\\\ " << "\n";


    // K1s
    SummaryFile3 << p_data.latexNameMod(paraNames[8])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[8] << " $\\pm$ " << errs_stat[8]  << " $\\pm$ " << errs_sys[8] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[9] << " $\\pm$ " << errs_stat[9]  << " $\\pm$ " << errs_sys[9] << " & ";
    SummaryFile3 << vals[10] << " $\\pm$ " << errs_stat[10]  << " $\\pm$ " << errs_sys[10] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[11] << " $\\pm$ " << errs_stat[11]  << " $\\pm$ " << errs_sys[11] << " \\\\ " << "\n";

    SummaryFile3 << p_data.latexNameMod(paraNames[12])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[12] << " $\\pm$ " << errs_stat[12]  << " $\\pm$ " << errs_sys[12] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[13] << " $\\pm$ " << errs_stat[13]  << " $\\pm$ " << errs_sys[13] << " & &  " << " \\\\ " << "\n";

    // K(1460)
    SummaryFile3 << p_data.latexNameMod(paraNames[14])  << " & & &";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[14] << " $\\pm$ " << errs_stat[14]  << " $\\pm$ " << errs_sys[14] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[15] << " $\\pm$ " << errs_stat[15]  << " $\\pm$ " << errs_sys[15] << " \\\\ " << "\n";

    // NV Ks
    SummaryFile3 << p_data.latexNameMod(paraNames[16])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[16] << " $\\pm$ " << errs_stat[16]  << " $\\pm$ " << errs_sys[16] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1);  
    SummaryFile3 << vals[17] << " $\\pm$ " << errs_stat[17]  << " $\\pm$ " << errs_sys[17] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[18] << " $\\pm$ " << errs_stat[18]  << " $\\pm$ " << errs_sys[18] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[19] << " $\\pm$ " << errs_stat[19]  << " $\\pm$ " << errs_sys[19] << " \\\\ " << "\n";

    // NV rho
    SummaryFile3 << p_data.latexNameMod(paraNames[20])  << " & ";
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << vals[20] << " $\\pm$ " << errs_stat[20]  << " $\\pm$ " << errs_sys[20] << " & ";
    SummaryFile3 << std::fixed << std::setprecision(1); 
    SummaryFile3 << vals[21] << " $\\pm$ " << errs_stat[21]  << " $\\pm$ " << errs_sys[21] << " & & \\\\ " << "\n";

    // masses,widths
    SummaryFile3 << std::fixed << std::setprecision(0);
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Fit parameter} & \\multicolumn{4}{c}{Value} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[22]) << "} & \\multicolumn{4}{c}{" << vals[22] << " $\\pm$ " << errs_stat[22] << " $\\pm$ " << errs_sys[22] << " $\\pm$ " << sqrt(cov_alt[22][22]) << "} \\\\ " << "\n";
    
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[23]) << "} & \\multicolumn{4}{c}{" << vals[23] << " $\\pm$ " << errs_stat[23] << " $\\pm$ " << errs_sys[23] << " $\\pm$ " << sqrt(cov_alt[23][23]) << "} \\\\ " << "\n";
    
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[24]) << "} & \\multicolumn{4}{c}{" << vals[24] << " $\\pm$ " << errs_stat[24] << " $\\pm$ " << errs_sys[24] << " $\\pm$ " << sqrt(cov_alt[24][24]) << "} \\\\ " << "\n";
    
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[25]) << "} & \\multicolumn{4}{c}{" << vals[25] << " $\\pm$ " << errs_stat[25] << " $\\pm$ " << errs_sys[25] << " $\\pm$ " << sqrt(cov_alt[25][25]) << "} \\\\ " << "\n";
    SummaryFile3 << " \\\\ " << "\n";

    // r, delta, gamma
    SummaryFile3 << std::fixed << std::setprecision(2);
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[26]) << "} & \\multicolumn{4}{c}{" << vals[26] << " $\\pm$ " << errs_stat[26] << " $\\pm$ " << errs_sys[26] << " $\\pm$ " << sqrt(cov_alt[26][26]) << "} \\\\ " << "\n";
    SummaryFile3 << std::fixed << std::setprecision(0);
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[27]) << "} & \\multicolumn{4}{c}{" << vals[27] << " $\\pm$ " << errs_stat[27] << " $\\pm$ " << errs_sys[27] << " $\\pm$ " << sqrt(cov_alt[27][27]) << "} \\\\ " << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{" << p_data.latexNameMod(paraNames[28]) << "} & \\multicolumn{4}{c}{" << vals[28] << " $\\pm$ " << errs_stat[28] << " $\\pm$ " << errs_sys[28] << " $\\pm$ " << sqrt(cov_alt[28][28]) << "} \\\\ " << "\n";
    
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";

    for(int i = 22 ; i <paraNames.size() ; i++){
        double tot = 0.;
        for(int j =8 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
	}
        tot += cov_alt[i][i];
        cout << "model sys = " << sqrt(tot) << endl;
     }

}

void fractions(){

   /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");

    paraNames.push_back("Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");

    paraNames.push_back("Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");

    paraNames.push_back("Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");

    paraNames.push_back("Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm");

    paraNames.push_back("Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm");

    paraNames.push_back("bar_Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm");

    paraNames.push_back("Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("bar_Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");

    paraNames.push_back("Bs0_NonResV0__DsmKp_rho_770_0__pippim_");

    paraNames.push_back("Sum");
    paraNames.push_back("bar_Sum");

    paraNames.push_back("k");


    vector<TMatrixD*> covs;

    pull p_data(paraNames,"signal_new2/pull__2.root","Coherence");
    pull p(paraNames,"out_signal_new_toy/pull__*.root","Coherence");
    //vector<double> vals = p.sampleMean()  ;
    vector<double> vals = p_data.getVals();
    vector<double> errs_stat = p.sampleSigma()  ;

    TMatrixD* cov = new TMatrixD(p.getCov());
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov)[i][j] = (*cov)[i][j] * errs_stat[i] * errs_stat[j];
    }
    covs.push_back(cov);

    /// Acc systematics (with cholesky)
    pull p_acc(paraNames,"out_signal_new_toy/pullAcc_*.root","Coherence");
    TMatrixD* cov_acc = new TMatrixD(p_acc.getDeltaCov("out_signal_new_toy/pull__*.root","_acc"));

    
    /// resolution systematics 
    pull p_res_Run1_a(paraNames,"out_sys_res/pull__1.root","Coherence");
    TMatrixD* cov_res_Run1_a = new TMatrixD(p_res_Run1_a.getDeltaCov("out/pull__1.root","_res_Run1_a"));

    pull p_res_Run1_b(paraNames,"out_sys_res/pull__2.root","Coherence");
    TMatrixD* cov_res_Run1_b = new TMatrixD(p_res_Run1_b.getDeltaCov("out/pull__1.root","_res_Run1_b"));

    pull p_res_Run2_a(paraNames,"out_sys_res/pull__3.root","Coherence");
    TMatrixD* cov_res_Run2_a = new TMatrixD(p_res_Run2_a.getDeltaCov("out/pull__1.root","_res_Run2_a"));

    pull p_res_Run2_b(paraNames,"out_sys_res/pull__4.root","Coherence");
    TMatrixD* cov_res_Run2_b = new TMatrixD(p_res_Run2_b.getDeltaCov("out/pull__1.root","_res_Run2_b"));

    pull p_res_Run2_c(paraNames,"out_sys_res/pull__5.root","Coherence");
    TMatrixD* cov_res_Run2_c = new TMatrixD(p_res_Run2_c.getDeltaCov("out/pull__1.root","_res_Run2_c"));

    pull p_res_Run2_d(paraNames,"out_sys_res/pull__6.root","Coherence");
    TMatrixD* cov_res_Run2_d = new TMatrixD(p_res_Run2_d.getDeltaCov("out/pull__1.root","_res_Run2_d"));

    pull p_res_Run2_e(paraNames,"out_sys_res/pull__7.root","Coherence");
    TMatrixD* cov_res_Run2_e = new TMatrixD(p_res_Run2_e.getDeltaCov("out/pull__1.root","_res_Run2_e"));

    pull p_res_Run2_f(paraNames,"out_sys_res/pull__8.root","Coherence");
    TMatrixD* cov_res_Run2_f = new TMatrixD(p_res_Run2_f.getDeltaCov("out/pull__1.root","_res_Run2_f"));


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

    vector<TMatrixD*> covs_res_18;
    covs_res_18.push_back(cov_res_Run2_e);
    covs_res_18.push_back(cov_res_Run2_f);
    cov_res +=  p_res_Run1_a.combineCov_maxVal(covs_res_18) ;

    /// decay time bias systematics 
    pull p_bias(paraNames,"out_sys_bias/pull__1.root","Coherence");
    TMatrixD* cov_bias = new TMatrixD(p_bias.getDeltaCov("out/pull__1.root","_bias"));

    /// dms systematics 
    pull p_dm(paraNames,"out_signal_new_toy/pull_dm_*.root","Coherence");
    TMatrixD* cov_dm = new TMatrixD(p_dm.getDeltaCov("out_signal_new_toy/pull__*.root","_dm"));

    /// asymmetry systematics
    pull p_asym(paraNames,"out_signal_new_toy/pull_asym_*.root","Coherence");
    TMatrixD* cov_asym = new TMatrixD(p_asym.getDeltaCov("out_signal_new_toy/pull__*.root","_asym"));

    /// bkg systematics 

    pull p_bkg_1(paraNames,"out_sys_bkg/pull__1.root","Coherence");
    vector<double> vals_bkg_1 = p_bkg_1.getVals();

    pull p_bkg_2(paraNames,"out_sys_bkg/pull__2.root","Coherence");
    vector<double> vals_bkg_2 = p_bkg_2.getVals();

    pull p_bkg_3(paraNames,"out_sys_bkg/pull__3.root","Coherence");
    vector<double> vals_bkg_3 = p_bkg_3.getVals();

    pull p_bkg_4(paraNames,"out_sys_bkg/pull__4.root","Coherence");
    vector<double> vals_bkg_4 = p_bkg_4.getVals();

    pull p_bkg_5(paraNames,"out_sys_bkg/pull__5.root","Coherence");
    vector<double> vals_bkg_5 = p_bkg_5.getVals();

    pull p_bkg_6(paraNames,"out_sys_bkg/pull__6.root","Coherence");
    vector<double> vals_bkg_6 = p_bkg_6.getVals();

    pull p_bkg_7(paraNames,"out_sys_bkg/pull__7.root","Coherence");
    vector<double> vals_bkg_7 = p_bkg_7.getVals();

    pull p_bkg_8(paraNames,"out_sys_bkg/pull__8.root","Coherence");
    vector<double> vals_bkg_8 = p_bkg_8.getVals();

    // Take sample variance as systematic 
    vector< vector <double> > vec_vals_bkg;
    //vec_vals_bkg.push_back(vals_bkg_1);
    vec_vals_bkg.push_back(vals_bkg_2);
    vec_vals_bkg.push_back(vals_bkg_3);
    vec_vals_bkg.push_back(vals_bkg_4);
    vec_vals_bkg.push_back(vals_bkg_5);
    vec_vals_bkg.push_back(vals_bkg_6);
    vec_vals_bkg.push_back(vals_bkg_7);
    vec_vals_bkg.push_back(vals_bkg_8);

    TMatrixD cov_bkg(p_bkg_2.sampleVariance(vec_vals_bkg,1));
//     cov_bkg.Print();

    pull p_bkg_bs1(paraNames,"out_signal_new_toy/sw_pull__*.root");
    vector<double> mean_bkg_bs1 = p_bkg_bs1.getPullMean();
    vector<double> sigma_bkg_bs1 = p_bkg_bs1.getPullSigma();

    pull p_bkg_bs2(paraNames,"out_signal_new_toy/pull__*.root");
    vector<double> mean_bkg_bs2 = p_bkg_bs2.getPullMean();
    vector<double> sigma_bkg_bs2 = p_bkg_bs2.getPullSigma();

    TMatrixD* cov_bkg_bs = new TMatrixD(cov_bkg);    
    for(int i =0 ; i <paraNames.size() ; i++)
	for(int j =0 ; j <paraNames.size() ; j++){
	if(i==j){
 		(*cov_bkg_bs)[i][j] = max(pow(min(mean_bkg_bs1[i]-mean_bkg_bs2[i],mean_bkg_bs1[i]),2) , pow(max(sigma_bkg_bs1[i]-sigma_bkg_bs2[i],0.),2))*errs_stat[i] * errs_stat[j];
	}
	else (*cov_bkg_bs)[i][j] = 0.;
     }
     cov_bkg+= *cov_bkg_bs;


    /// Lineshape models systematics

    pull p_ls_1(paraNames,"out_sys_ls/pull__1.root","Coherence");
    vector<double> vals_ls_1 = p_ls_1.getVals();
    cout << "ls 1" << endl;
    cout << vals_ls_1[4] << endl << endl;

    pull p_ls_2(paraNames,"out_sys_ls/pull__2.root","Coherence");
    vector<double> vals_ls_2 = p_ls_2.getVals();
    cout << "ls 2" << endl;
    cout << vals_ls_2[4] << endl << endl;

    pull p_ls_3(paraNames,"out_sys_ls/pull__3.root","Coherence");
    vector<double> vals_ls_3 = p_ls_3.getVals();
    cout << "ls 3" << endl;
    cout << vals_ls_3[4] << endl << endl;

    pull p_ls_4(paraNames,"out_sys_ls/pull__4.root","Coherence");
    vector<double> vals_ls_4 = p_ls_4.getVals();
    cout << vals_ls_4[4] << endl << endl;

    pull p_ls_5(paraNames,"out_sys_ls/pull__5.root","Coherence");
    vector<double> vals_ls_5 = p_ls_5.getVals();
    cout << vals_ls_5[4] << endl << endl;

    pull p_ls_6(paraNames,"out_sys_ls/pull__6.root","Coherence");
    vector<double> vals_ls_6 = p_ls_6.getVals();
    cout << vals_ls_6[4] << endl << endl;

    pull p_ls_7(paraNames,"out_sys_ls/pull__7.root","Coherence");
    vector<double> vals_ls_7 = p_ls_7.getVals();
    cout << vals_ls_7[4] << endl << endl;

    vector< vector <double> > vec_vals_ls;
    vec_vals_ls.push_back(vals_ls_1);
    vec_vals_ls.push_back(vals_ls_2);
    vec_vals_ls.push_back(vals_ls_3);
    vec_vals_ls.push_back(vals_ls_4);
    vec_vals_ls.push_back(vals_ls_5);
    vec_vals_ls.push_back(vals_ls_6);
    vec_vals_ls.push_back(vals_ls_7);
    TMatrixD cov_ls(p_ls_1.sampleVariance(vec_vals_ls));


    /// Resonance parameters systematics
    pull p_rp_1(paraNames,"out_signal_new_toy/pull_ls1_*.root","Coherence");
    TMatrixD* cov_rp_1 = new TMatrixD(p_rp_1.getDeltaCov("out_signal_new_toy/pull__*.root","_ls1",1));

    pull p_rp_2(paraNames,"out_signal_new_toy/pull_ls2_*.root","Coherence");
    TMatrixD* cov_rp_2 = new TMatrixD(p_rp_2.getDeltaCov("out_signal_new_toy/pull__*.root","_ls2",1));

    // Add
    TMatrixD cov_rp(*cov_rp_1);
    cov_rp +=  *cov_rp_2 ;
 
    /// Form factor
    pull p_f_1(paraNames,"out_signal_new_toy/pull_BW_radius_*.root","Coherence");
    TMatrixD* cov_f_1 = new TMatrixD(p_f_1.getDeltaCov("out_signal_new_toy/pull__*.root",5));

    /// Phsp-Acc systematics
    pull p_phsp_acc_1(paraNames,"out_sys_acc/pull__1.root","Coherence");
    vector<double> vals_ps_1 = p_phsp_acc_1.getVals();

    pull p_phsp_acc_2(paraNames,"out_sys_acc/pull__2.root","Coherence");
    vector<double> vals_ps_2 = p_phsp_acc_2.getVals();

    pull p_phsp_acc_3(paraNames,"out_sys_acc/pull__3.root","Coherence");
    vector<double> vals_ps_3 = p_phsp_acc_3.getVals();

    pull p_phsp_acc_4(paraNames,"out_sys_acc/pull__4.root","Coherence");
    vector<double> vals_ps_4 = p_phsp_acc_4.getVals();

    pull p_phsp_acc_5(paraNames,"out_sys_acc/pull__5.root","Coherence");
    vector<double> vals_ps_5 = p_phsp_acc_5.getVals();

    pull p_phsp_acc_6(paraNames,"out_sys_acc/pull__6.root","Coherence");
    vector<double> vals_ps_6 = p_phsp_acc_6.getVals();

    pull p_phsp_acc_7(paraNames,"out_sys_acc/pull__7.root","Coherence");
    vector<double> vals_ps_7 = p_phsp_acc_7.getVals();

    pull p_phsp_acc_8(paraNames,"out_sys_acc/pull__8.root","Coherence");
    vector<double> vals_ps_8 = p_phsp_acc_8.getVals();

    pull p_phsp_acc_9(paraNames,"out_sys_acc/pull__9.root","Coherence");
    vector<double> vals_ps_9 = p_phsp_acc_9.getVals();


    vector< vector <double> > vec_vals_ps;
    vec_vals_ps.push_back(vals_ps_1);
    vec_vals_ps.push_back(vals_ps_2);
    vec_vals_ps.push_back(vals_ps_3);
    vec_vals_ps.push_back(vals_ps_4);
    vec_vals_ps.push_back(vals_ps_5);
    vec_vals_ps.push_back(vals_ps_6);
    vec_vals_ps.push_back(vals_ps_7);
    vec_vals_ps.push_back(vals_ps_8);
    vec_vals_ps.push_back(vals_ps_9);

    TMatrixD cov_phsp(p_phsp_acc_1.sampleVariance(vec_vals_ps));

    /// Phsp acc factorization
    pull p_phsp_acc2_0(paraNames,"out_sys_acc2/pull__0.root","Coherence");
    vector<double> vals_ps_acc0 = p_phsp_acc2_0.getVals();

    pull p_phsp_acc2_1(paraNames,"out_sys_acc2/pull__1.root","Coherence");
    vector<double> vals_ps_acc1 = p_phsp_acc2_1.getVals();

    pull p_phsp_acc2_2(paraNames,"out_sys_acc2/pull__2.root","Coherence");
    vector<double> vals_ps_acc2 = p_phsp_acc2_2.getVals();

    pull p_phsp_acc2_3(paraNames,"out_sys_acc/pull__3.root","Coherence");
    vector<double> vals_ps_acc3 = p_phsp_acc2_3.getVals();
	
    vector< vector <double> > vec_vals_ps2;
    vec_vals_ps2.push_back(vals_ps_acc0);
    vec_vals_ps2.push_back(vals_ps_acc1);
    vec_vals_ps2.push_back(vals_ps_acc2);
    vec_vals_ps2.push_back(vals_ps_acc3);

    TMatrixD cov_phsp2(p_phsp_acc2_1.sampleVariance(vec_vals_ps2));


    /// Alternative amp models 
    pull p_alt_1(paraNames,"out_alt/pull__1.root","Coherence");
    vector<double> vals_alt_1 = p_alt_1.getVals();

    pull p_alt_2(paraNames,"out_alt/pull__9992.root","Coherence");
    vector<double> vals_alt_2 = p_alt_2.getVals();

    pull p_alt_3(paraNames,"out_alt/pull__93.root","Coherence");
    vector<double> vals_alt_3 = p_alt_3.getVals();

    pull p_alt_4(paraNames,"out_alt/pull__94.root","Coherence");
    vector<double> vals_alt_4 = p_alt_4.getVals();

    pull p_alt_5(paraNames,"out_alt/pull__5.root","Coherence");
    vector<double> vals_alt_5 = p_alt_5.getVals();

    pull p_alt_6(paraNames,"out_alt/pull__6.root","Coherence");
    vector<double> vals_alt_6 = p_alt_6.getVals();

    pull p_alt_7(paraNames,"out_alt/pull__9997.root","Coherence");
    vector<double> vals_alt_7 = p_alt_7.getVals();

    pull p_alt_8(paraNames,"out_alt/pull__9998.root","Coherence");
    vector<double> vals_alt_8 = p_alt_8.getVals();

    pull p_alt_9(paraNames,"out_alt/pull__99.root","Coherence");
    vector<double> vals_alt_9 = p_alt_9.getVals();

    pull p_alt_10(paraNames,"out_alt/pull__999997.root","Coherence");
    vector<double> vals_alt_10 = p_alt_10.getVals();

    pull p_alt_11(paraNames,"out_alt/pull__14.root","Coherence");
    vector<double> vals_alt_11 = p_alt_11.getVals();

    pull p_alt_12(paraNames,"out_alt/pull__912.root","Coherence");
    vector<double> vals_alt_12 = p_alt_12.getVals();

    pull p_alt_13(paraNames,"out_alt13/pull__7.root","Coherence");
    vector<double> vals_alt_13 = p_alt_13.getVals();

    //pull p_alt_14(paraNames,"out_alt/pull__14.root","Coherence");
    //vector<double> vals_alt_14 = p_alt_14.getVals();


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
    //vec_vals_alt.push_back(vals_alt_13);
    //vec_vals_alt.push_back(vals_alt_14);

    cout << "Sample variance " << endl;
    TMatrixD cov_alt(p_alt_1.sampleVariance(vec_vals_alt));

    /// m,t correlations
    pull p_corr1(paraNames,"out_signal_new_toy/sw_pull__*.root","Coherence");
    TMatrixD cov_corr1 = p_corr1.getCov();

    pull p_corr2(paraNames,"out_signal_new_toy_bkg/sw_pull__*.root","Coherence");
    TMatrixD cov_corr2 = p_corr2.getCov();

    TMatrixD* cov_corr = new TMatrixD(p_corr1.getAbsDiff(cov_corr1,cov_corr2));
    for(int i = 0 ; i < paraNames.size(); i++)for(int j = 0 ; j < paraNames.size(); j++){
	(*cov_corr)[i][j] = (*cov_corr)[i][j] * errs_stat[i] * errs_stat[j];
    }

    /// Total systematics table   
    covs.push_back(new TMatrixD(cov_bkg));
    covs.push_back(cov_corr);
    covs.push_back(cov_acc);
    covs.push_back(new TMatrixD(cov_res));
    covs.push_back(cov_bias);
    covs.push_back(cov_asym);
    covs.push_back(cov_dm);
    covs.push_back(new TMatrixD(cov_phsp));
    covs.push_back(new TMatrixD(cov_phsp2));
    covs.push_back(new TMatrixD(cov_ls));
    covs.push_back(new TMatrixD(cov_rp));
    covs.push_back(cov_f_1);

    vector<string> sysNames;
    sysNames.push_back("Fit bias");
    sysNames.push_back("Correlations");
    sysNames.push_back("Background");
     sysNames.push_back("Time-Acc.");
     sysNames.push_back("Resolution");
     sysNames.push_back("Decay-time bias");
     sysNames.push_back("Asymmetries");
     sysNames.push_back("$\\Delta m_{s}$");
    sysNames.push_back("Phsp-Acc.");
    sysNames.push_back("Acc. Factor.");
    sysNames.push_back("Lineshapes");
    sysNames.push_back("Resonances $m,\\Gamma$");
    sysNames.push_back("Form-Factors");
     sysNames.push_back("Amp. Model");

    vector<double> errs_sys;    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
        }
	errs_sys.push_back(sqrt(tot)); 
    }


    /// Result table   
    ofstream SummaryFile3;
    SummaryFile3.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/fraction_table.tex",std::ofstream::trunc);

    SummaryFile3 << "\\begin{tabular}{l r r } " << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\multicolumn{1}{c}{Decay Channel} & \\multicolumn{1}{c}{$F_{b \\to c} [\\%]$} & \\multicolumn{1}{c}{$F_{b \\to u} [\\%]$} " << " \\\\ " << "\n";
    SummaryFile3 << "\\hline" << "\n";

    SummaryFile3 << std::fixed << std::setprecision(1);
    // K1(1270)
    SummaryFile3 << p.latexNameMod(paraNames[0]) ; 
    SummaryFile3 << " & "  << vals[0] * 100. << " $\\pm$ " << errs_stat[0]* 100. << " $\\pm$ " << errs_sys[0]* 100. << " $\\pm$ " << sqrt(cov_alt[0][0]) * 100. ;
    SummaryFile3 << " & "  << vals[1]* 100. << " $\\pm$ " << errs_stat[1]* 100. << " $\\pm$ " << errs_sys[1]* 100. << " $\\pm$ " << sqrt(cov_alt[1][1]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[2]) ; 
    SummaryFile3 << " & "  << vals[2] * 100. << " $\\pm$ " << errs_stat[2]* 100. << " $\\pm$ " << errs_sys[2]* 100.  << " $\\pm$ " << sqrt(cov_alt[2][2]) * 100.;
    SummaryFile3 << " & "  << vals[3]* 100. << " $\\pm$ " << errs_stat[3]* 100. << " $\\pm$ " << errs_sys[3]* 100.  << " $\\pm$ " << sqrt(cov_alt[3][3]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[4]) ; 
    SummaryFile3 << " & "  << vals[4] * 100. << " $\\pm$ " << errs_stat[4]* 100. << " $\\pm$ " << errs_sys[4]* 100.  << " $\\pm$ " << sqrt(cov_alt[4][4]) * 100.;
    SummaryFile3 << " & "  << vals[5]* 100. << " $\\pm$ " << errs_stat[5]* 100. << " $\\pm$ " << errs_sys[5]* 100.  << " $\\pm$ " << sqrt(cov_alt[5][5]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    // K1(1400)
    SummaryFile3 << p.latexNameMod(paraNames[6]) ; 
    SummaryFile3 << " & "  << vals[6] * 100. << " $\\pm$ " << errs_stat[6]* 100. << " $\\pm$ " << errs_sys[6]* 100.  << " $\\pm$ " << sqrt(cov_alt[6][6]) * 100.;
    SummaryFile3 << " & "  << vals[7]* 100. << " $\\pm$ " << errs_stat[7]* 100. << " $\\pm$ " << errs_sys[7]* 100.  << " $\\pm$ " << sqrt(cov_alt[7][7]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    // Ks(1410)
    SummaryFile3 << p.latexNameMod(paraNames[8]) ; 
    SummaryFile3 << " & "  << vals[8] * 100. << " $\\pm$ " << errs_stat[8]* 100. << " $\\pm$ " << errs_sys[8]* 100.  << " $\\pm$ " << sqrt(cov_alt[8][8]) * 100.;
    SummaryFile3 << " & "  << vals[9] * 100. << " $\\pm$ " << errs_stat[9]* 100. << " $\\pm$ " << errs_sys[9]* 100.  << " $\\pm$ " << sqrt(cov_alt[9][9]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[10]) ; 
    SummaryFile3 << " & "  << vals[10] * 100. << " $\\pm$ " << errs_stat[10]* 100. << " $\\pm$ " << errs_sys[10]* 100.  << " $\\pm$ " << sqrt(cov_alt[10][10]) * 100.;
    SummaryFile3 << " & "  << vals[11] * 100. << " $\\pm$ " << errs_stat[11]* 100. << " $\\pm$ " << errs_sys[11]* 100.  << " $\\pm$ " << sqrt(cov_alt[11][11]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    // K(1460)
    SummaryFile3 << p.latexNameMod(paraNames[12]) ; 
    SummaryFile3 << " & "  ;
    SummaryFile3 << " & "  << vals[12]* 100. << " $\\pm$ " << errs_stat[12]* 100. << " $\\pm$ " << errs_sys[12]* 100.  << " $\\pm$ " << sqrt(cov_alt[12][12]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    // NS
    SummaryFile3 << p.latexNameMod(paraNames[13]) ; 
    SummaryFile3 << " & "  << vals[13] * 100. << " $\\pm$ " << errs_stat[13]* 100. << " $\\pm$ " << errs_sys[13]* 100.  << " $\\pm$ " << sqrt(cov_alt[13][13]) * 100.;
    SummaryFile3 << " & "  << vals[14]* 100. << " $\\pm$ " << errs_stat[14]* 100. << " $\\pm$ " << errs_sys[14]* 100.  << " $\\pm$ " << sqrt(cov_alt[14][14]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << p.latexNameMod(paraNames[15]) ; 
    SummaryFile3 << " & "  << vals[15]* 100. << " $\\pm$ " << errs_stat[15]* 100. << " $\\pm$ " << errs_sys[15]* 100.  << " $\\pm$ " << sqrt(cov_alt[15][15]) * 100.;
    SummaryFile3 << " & "   ;
    SummaryFile3 << " \\\\ " << "\n";

    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << p.latexNameMod(paraNames[16]) ; 
    SummaryFile3 << " & "  << vals[16] * 100. << " $\\pm$ " << errs_stat[16]* 100. << " $\\pm$ " << errs_sys[16]* 100.  << " $\\pm$ " << sqrt(cov_alt[16][16]) * 100.;
    SummaryFile3 << " & "  << vals[17]* 100. << " $\\pm$ " << errs_stat[17]* 100. << " $\\pm$ " << errs_sys[17]* 100.  << " $\\pm$ " << sqrt(cov_alt[17][17]) * 100.;
    SummaryFile3 << " \\\\ " << "\n";
    
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\hline" << "\n";
    SummaryFile3 << "\\end{tabular}" << "\n";

    cout << "kappa = " << vals[18] << " $\\pm$ " << errs_stat[18] << " $\\pm$ " << errs_sys[18]  << " $\\pm$ " << sqrt(cov_alt[18][18]) << endl;
    cout << "kappa = " << vals[18] << " $\\pm$ " << sqrt(pow(errs_stat[18],2) + pow(errs_sys[18],2) + cov_alt[18][18] ) << endl;  


    for(int i = paraNames.size()-1 ; i <paraNames.size() ; i++){
        double tot = 0.;
        for(int j =8 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
	}
        tot += cov_alt[i][i];
        cout << "model sys = " << sqrt(tot) << endl;
     }

}

void altModels(){

   /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p_D___Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p__rho_1450_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1__1400_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm");
//     paraNames.push_back("Bs0_K_2_s_1430_p__Ks_892_0__Kppim_pip_Dsm");
//     paraNames.push_back("Bs0_K_2_s_1430_p__rho_770_0__pippim_Kp_Dsm");    
    paraNames.push_back("Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_K_1460_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_K_1460_p__sigma10__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_Ks_1680_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("Bs0_Ks_1680_p__rho_770_0__pippim_Kp_Dsm");    
    paraNames.push_back("Bs0_K_2__1770_p__Ks_892_0__Kppim_pip_Dsm");
//     paraNames.push_back("Bs0_K_2__1770_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("Bs0_NonResS0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("Bs0_P__NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("Bs0_D__NonResV0__Dsmpip_Ks_892_0__Kppim_");
//     paraNames.push_back("Bs0_K_0_s_1430_0__Kppim_NonResS0__Dsmpip_");
    paraNames.push_back("Bs0_NonResS0__DsmKp_sigma10__pippim_");
    //paraNames.push_back("Bs0_NonResV0__DsmKp_sigma10__pippim_");
    paraNames.push_back("Bs0_NonResS0__DsmKp_f_0__980_0__pippim_");
    paraNames.push_back("Bs0_f_2__1270_0__pippim_NonResS0__DsmKp_");
    //paraNames.push_back("Bs0_f_2__1270_0__pippim_NonResV0__DsmKp_");
    //paraNames.push_back("Bs0_f_0__1370_0__pippim_NonResS0__DsmKp_");
    paraNames.push_back("Bs0_NonResS0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("Bs0_NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("Bs0_P__NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("Bs0_D__NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("Sum");
    int N_1 = paraNames.size();

    paraNames.push_back("bar_Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p_D___Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__rho_1450_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1__1400_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm");
//     paraNames.push_back("bar_Bs0_K_2_s_1430_p__Ks_892_0__Kppim_pip_Dsm");
//     paraNames.push_back("bar_Bs0_K_2_s_1430_p__rho_770_0__pippim_Kp_Dsm");    
    paraNames.push_back("bar_Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_K_1460_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_K_1460_p__sigma10__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1680_p__Ks_892_0__Kppim_pip_Dsm");
    paraNames.push_back("bar_Bs0_Ks_1680_p__rho_770_0__pippim_Kp_Dsm");    
    paraNames.push_back("bar_Bs0_K_2__1770_p__Ks_892_0__Kppim_pip_Dsm");
//     paraNames.push_back("bar_Bs0_K_2__1770_p__rho_770_0__pippim_Kp_Dsm");
    paraNames.push_back("bar_Bs0_NonResS0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("bar_Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("bar_Bs0_P__NonResV0__Dsmpip_Ks_892_0__Kppim_");
    paraNames.push_back("bar_Bs0_D__NonResV0__Dsmpip_Ks_892_0__Kppim_");
//     paraNames.push_back("bar_Bs0_K_0_s_1430_0__Kppim_NonResS0__Dsmpip_");
    paraNames.push_back("bar_Bs0_NonResS0__DsmKp_sigma10__pippim_");
    //paraNames.push_back("bar_Bs0_NonResV0__DsmKp_sigma10__pippim_");
    paraNames.push_back("bar_Bs0_NonResS0__DsmKp_f_0__980_0__pippim_");
    paraNames.push_back("bar_Bs0_f_2__1270_0__pippim_NonResS0__DsmKp_");
    //paraNames.push_back("bar_Bs0_f_2__1270_0__pippim_NonResV0__DsmKp_");
    //paraNames.push_back("bar_Bs0_f_0__1370_0__pippim_NonResS0__DsmKp_");
    paraNames.push_back("bar_Bs0_NonResS0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("bar_Bs0_NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("bar_Bs0_P__NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("bar_Bs0_D__NonResV0__DsmKp_rho_770_0__pippim_");
    paraNames.push_back("bar_Sum");
    int N_2 = paraNames.size();

    paraNames.push_back("mass_K_1__1400_p");
    paraNames.push_back("width_K_1__1400_p");
    paraNames.push_back("mass_Ks_1410_p");
    paraNames.push_back("width_Ks_1410_p");

    paraNames.push_back("r");
    paraNames.push_back("k");
    paraNames.push_back("delta");
    paraNames.push_back("gamma");
    paraNames.push_back("n2ll");

    pull p(paraNames,"signal_new2/pull__2.root","Coherence");
    pull p_err(paraNames,"out/pull__2.root","Coherence");
    vector<double> vals = p.getVals()  ;
    vector<double> errs_stat = p_err.getErrs()  ;

    pull p_1(paraNames,"out_alt/pull__1.root","Coherence");
    vector<double> vals_1 = p_1.getVals()  ;
    vector<double> errs_stat_1 = p_1.getErrs()  ;

    pull p_2(paraNames,"out_alt/pull__9992.root","Coherence");
    vector<double> vals_2 = p_2.getVals()  ;
    vector<double> errs_stat_2 = p_2.getErrs()  ;

    pull p_3(paraNames,"out_alt/pull__93.root","Coherence");
    vector<double> vals_3 = p_3.getVals()  ;
    vector<double> errs_stat_3 = p_3.getErrs()  ;

    pull p_4(paraNames,"out_alt/pull__94.root","Coherence");
    vector<double> vals_4 = p_4.getVals()  ;
    vector<double> errs_stat_4 = p_4.getErrs()  ;

    pull p_5(paraNames,"out_alt/pull__5.root","Coherence");
    vector<double> vals_5 = p_5.getVals()  ;
    vector<double> errs_stat_5 = p_5.getErrs()  ;

    pull p_6(paraNames,"out_alt/pull__6.root","Coherence");
    vector<double> vals_6 = p_6.getVals()  ;
    vector<double> errs_stat_6 = p_6.getErrs()  ;

    pull p_7(paraNames,"out_alt/pull__9997.root","Coherence");
    vector<double> vals_7 = p_7.getVals()  ;
    vector<double> errs_stat_7 = p_7.getErrs()  ;

    pull p_8(paraNames,"out_alt/pull__9998.root","Coherence");
    vector<double> vals_8 = p_8.getVals()  ;
    vector<double> errs_stat_8 = p_8.getErrs()  ;

    pull p_9(paraNames,"out_alt/pull__99.root","Coherence");
    vector<double> vals_9 = p_9.getVals()  ;
    vector<double> errs_stat_9 = p_9.getErrs()  ;

    pull p_10(paraNames,"out_alt/pull__999997.root","Coherence");
    vector<double> vals_10 = p_10.getVals()  ;
    vector<double> errs_stat_10 = p_10.getErrs()  ;

    pull p_11(paraNames,"out_alt/pull__11.root","Coherence");
    vector<double> vals_11 = p_11.getVals()  ;
    vector<double> errs_stat_11 = p_11.getErrs()  ;

    pull p_12(paraNames,"out_alt/pull__912.root","Coherence");
    vector<double> vals_12 = p_12.getVals()  ;
    vector<double> errs_stat_12 = p_12.getErrs()  ;

    pull p_13(paraNames,"out_alt13/pull__7.root","Coherence");
    vector<double> vals_13 = p_13.getVals()  ;
    vector<double> errs_stat_13 = p_13.getErrs()  ;

    pull p_14(paraNames,"out_alt/pull__14.root","Coherence");
    vector<double> vals_14 = p_14.getVals()  ;
    vector<double> errs_stat_14 = p_14.getErrs()  ;

  

    vector< vector<double> > vec_vals;
    vec_vals.push_back(vals);
    vec_vals.push_back(vals_1);
    vec_vals.push_back(vals_2);
    vec_vals.push_back(vals_3);
    vec_vals.push_back(vals_4);
    vec_vals.push_back(vals_5);
    vec_vals.push_back(vals_6);
    vec_vals.push_back(vals_7);
    vec_vals.push_back(vals_8);
    vec_vals.push_back(vals_9);
    vec_vals.push_back(vals_10);
    //vec_vals.push_back(vals_11);
    vec_vals.push_back(vals_12);
    vec_vals.push_back(vals_14);
    //vec_vals.push_back(vals_13);

    vector< vector<double> > vec_errs;
    vec_errs.push_back(errs_stat);
    vec_errs.push_back(errs_stat_1);
    vec_errs.push_back(errs_stat_2);
    vec_errs.push_back(errs_stat_3);
    vec_errs.push_back(errs_stat_4);
    vec_errs.push_back(errs_stat_5);
    vec_errs.push_back(errs_stat_6);
    vec_errs.push_back(errs_stat_7);
    vec_errs.push_back(errs_stat_8);
    vec_errs.push_back(errs_stat_9);
    vec_errs.push_back(errs_stat_10);
    //vec_errs.push_back(errs_stat_11);
    vec_errs.push_back(errs_stat_12);
    vec_errs.push_back(errs_stat_14);
    //vec_errs.push_back(errs_stat_13);

    int N_modelsPerTable = 7;
    /// Result table   
    ofstream SummaryFile;
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/altModel_table.tex",std::ofstream::trunc);

    SummaryFile << "\\begin{tabular}{l l " ;
    for(int j=0; j < N_modelsPerTable; j++) SummaryFile << " r " ;
    SummaryFile << " } " << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "& & \\multicolumn{1}{c}{Baseline} " ;
    for(int j=1; j < N_modelsPerTable; j++)SummaryFile << " & \\multicolumn{1}{c}{Alt." << j <<  "} ";
    SummaryFile << " \\\\ " << "\n";
    SummaryFile << "\\hline" << "\n";

    SummaryFile << std::fixed << std::setprecision(1);
 
    for(int i= 0; i < paraNames.size(); i++){

	    if(i==0) SummaryFile << "\\multirow{" << N_1 << "}{*}{$b \\to c$} ";
	    if(i==N_1-1) SummaryFile << "\\multirow{" << N_1 << "}{*}{$b \\to u$} ";
            SummaryFile << " & " ;
	    SummaryFile << p.latexNameMod(paraNames[i]) ;
	    double scale = (i>N_2-1) ? 1. : 100.; 
	    if(i == paraNames.size()-9)SummaryFile << std::fixed << std::setprecision(0);
	    if(i == paraNames.size()-5)SummaryFile << std::fixed << std::setprecision(2);
	    if(i == paraNames.size()-3)SummaryFile << std::fixed << std::setprecision(0);

	    for(int j=0; j < N_modelsPerTable; j++){
		    if(i < paraNames.size()-9 && vec_vals[j][i] < 0.001/100.) SummaryFile << " & - " ;
 		    //else if(i==paraNames.size()-1) SummaryFile << " & "  << vec_vals[j][i] - vec_vals[0][i];
		    else if(i<paraNames.size()-5)SummaryFile << " & "  << vec_vals[j][i] * scale;
		    else if(i==paraNames.size()-5)SummaryFile << " & "  << vec_vals[j][i];// - vec_vals[0][i]; // + 0.50;
		    else if(i==paraNames.size()-4)SummaryFile << " & "  << vec_vals[j][i];// - vec_vals[0][i]; // + 0.52;
		    else if(i==paraNames.size()-3)SummaryFile << " & "  << vec_vals[j][i];// - vec_vals[0][i]; // + 46.;
		    else if(i==paraNames.size()-2)SummaryFile << " & "  << vec_vals[j][i];// - vec_vals[0][i]; // + 61.;
	    }	    
	    SummaryFile << " \\\\ " << "\n";
	    if(i==N_1-1 || i==N_2-1) SummaryFile << "\\hline" << "\n";
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";


    /// Result table   
    ofstream SummaryFile2;
    SummaryFile2.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal_new/altModel_table2.tex",std::ofstream::trunc);

    SummaryFile2 << "\\begin{tabular}{l l " ;
    for(int j=N_modelsPerTable; j < vec_vals.size(); j++) SummaryFile2 << " r " ;
    SummaryFile2 << " } " << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "& " ;
    for(int j=N_modelsPerTable; j < vec_vals.size(); j++)SummaryFile2 << " & \\multicolumn{1}{c}{Alt." << j <<  "} ";
    SummaryFile2 << " \\\\ " << "\n";
    SummaryFile2 << "\\hline" << "\n";

    SummaryFile2 << std::fixed << std::setprecision(1);
 
    for(int i= 0; i < paraNames.size(); i++){

	    if(i==0) SummaryFile2 << "\\multirow{" << N_1 << "}{*}{$b \\to c$} ";
	    if(i==N_1-1) SummaryFile2 << "\\multirow{" << N_1 << "}{*}{$b \\to u$} ";
            SummaryFile2 << " & " ;
	    SummaryFile2 << p.latexNameMod(paraNames[i]) ;
	    double scale = (i>N_2-1) ? 1. : 100.; 
	    if(i == paraNames.size()-9)SummaryFile2 << std::fixed << std::setprecision(0);
	    if(i == paraNames.size()-5)SummaryFile2 << std::fixed << std::setprecision(2);
	    if(i == paraNames.size()-3)SummaryFile2 << std::fixed << std::setprecision(0);

	    for(int j=N_modelsPerTable; j < vec_vals.size(); j++){
		    if(i < paraNames.size()-9 && vec_vals[j][i] < 0.001/100.) SummaryFile2 << " & - " ;
 		    //else if(i==paraNames.size()-1) SummaryFile2 << " & "  << vec_vals[j][i] - vec_vals[0][i] ; 		 
		    else if(i<paraNames.size()-5)SummaryFile2 << " & "  << vec_vals[j][i] * scale;
		    else if(i==paraNames.size()-5)SummaryFile2 << " & "  << vec_vals[j][i];// - vec_vals[0][i];// + 0.50;
		    else if(i==paraNames.size()-4)SummaryFile2 << " & "  << vec_vals[j][i];// - vec_vals[0][i] ;//+ 0.52;
		    else if(i==paraNames.size()-3)SummaryFile2 << " & "  << vec_vals[j][i];// - vec_vals[0][i] ;//+ 46.;
		    else if(i==paraNames.size()-2)SummaryFile2 << " & "  << vec_vals[j][i];// - vec_vals[0][i] ;//+ 61.;
	    }	    
	    SummaryFile2 << " \\\\ " << "\n";
	    if(i==N_1-1 || i==N_2-1) SummaryFile2 << "\\hline" << "\n";
    }
    
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\hline" << "\n";
    SummaryFile2 << "\\end{tabular}" << "\n";


}

void corrCP(){

    vector<TString> paraNames;
    paraNames.push_back("r");
    paraNames.push_back("delta");
    paraNames.push_back("gamma");
    paraNames.push_back("k");

    pull p(paraNames,"signal_toy/pull__*.root","Coherence");
    vector<double> vals = p.sampleMean()  ;
    vector<double> sigma; // = p.sampleSigma()  ;

    TMatrixD cov_stat = p.getStatCov();
    cov_stat.Print();

    sigma.push_back(sqrt(cov_stat[0][0]));
    sigma.push_back(sqrt(cov_stat[1][1]));
    sigma.push_back(sqrt(cov_stat[2][2]));
    sigma.push_back(sqrt(cov_stat[3][3]));

    for(int i=0; i < cov_stat.GetNcols(); i++){
        for(int j=0; j < cov_stat.GetNcols(); j++){    
            cov_stat(i,j) = cov_stat(i,j) / sigma[i] / sigma[j];
        }
    }

    cov_stat.Print();


    vector<TMatrixD*> covs;

    /// Alternative amp models 
    pull p_alt_1(paraNames,"signal_alt0/pull__*.root","Coherence");
    TMatrixD* cov_alt_1 = new TMatrixD(p_alt_1.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_1 = p_alt_1.getVals();

    pull p_alt_2(paraNames,"signal_alt2/pull__*.root","Coherence");
    TMatrixD* cov_alt_2 = new TMatrixD(p_alt_2.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_2 = p_alt_2.getVals();

    pull p_alt_3(paraNames,"signal_alt3/pull__*.root","Coherence");
    TMatrixD* cov_alt_3 = new TMatrixD(p_alt_3.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_3 = p_alt_3.getVals();

    pull p_alt_4(paraNames,"signal_alt5/pull__*.root","Coherence");
    TMatrixD* cov_alt_4 = new TMatrixD(p_alt_4.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_4 = p_alt_4.getVals();

    pull p_alt_5(paraNames,"signal_alt6/pull__*.root","Coherence");
    TMatrixD* cov_alt_5 = new TMatrixD(p_alt_5.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_5 = p_alt_5.getVals();

    pull p_alt_6(paraNames,"signal_alt8/pull__*.root","Coherence");
    TMatrixD* cov_alt_6 = new TMatrixD(p_alt_6.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_6 = p_alt_6.getVals();

    pull p_alt_7(paraNames,"signal_alt11/pull__*.root","Coherence");
    TMatrixD* cov_alt_7 = new TMatrixD(p_alt_7.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_7 = p_alt_7.getVals();

    pull p_alt_8(paraNames,"signal_alt12/pull__*.root","Coherence");
    TMatrixD* cov_alt_8 = new TMatrixD(p_alt_8.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_8 = p_alt_8.getVals();

    pull p_alt_9(paraNames,"signal_alt15/pull__*.root","Coherence");
    TMatrixD* cov_alt_9 = new TMatrixD(p_alt_9.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_9 = p_alt_9.getVals();

    pull p_alt_10(paraNames,"signal_alt16/pull__*.root","Coherence");
    TMatrixD* cov_alt_10 = new TMatrixD(p_alt_10.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_10 = p_alt_10.getVals();

    pull p_alt_11(paraNames,"signal_alt10/pull__*.root","Coherence");
    TMatrixD* cov_alt_11 = new TMatrixD(p_alt_11.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_11 = p_alt_11.getVals();

    pull p_alt_12(paraNames,"signal_alt18/pull__*.root","Coherence");
    TMatrixD* cov_alt_12 = new TMatrixD(p_alt_12.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_12 = p_alt_12.getVals();

    pull p_alt_13(paraNames,"signal_alt19/pull__*.root","Coherence");
    TMatrixD* cov_alt_13 = new TMatrixD(p_alt_13.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_13 = p_alt_13.getVals();

    pull p_alt_14(paraNames,"signal_alt20/pull__*.root","Coherence");
    TMatrixD* cov_alt_14 = new TMatrixD(p_alt_14.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_14 = p_alt_14.getVals();

    pull p_alt_15(paraNames,"signal_alt22/pull__*.root","Coherence");
    TMatrixD* cov_alt_15 = new TMatrixD(p_alt_15.getDeltaCov("signal/pull__1.root"));
    vector<double> vals_alt_15 = p_alt_15.getVals();

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

    cout << "Sample variance " << endl;
    TMatrixD cov_alt(p_alt_1.sampleVariance(vec_vals_alt));


   /// Lineshape models systematics
    pull p_ls_1(paraNames,"signal_sys1/pull_*.root","Coherence");
    TMatrixD* cov_ls_1 = new TMatrixD(p_ls_1.getDeltaCov("signal/pull__1.root","ls_1"));
    vector<double> vals_ls_1 = p_ls_1.getVals();

    pull p_ls_2(paraNames,"signal_sys2/pull__2.root","Coherence");
    TMatrixD* cov_ls_2 = new TMatrixD(p_ls_2.getDeltaCov("signal/pull__1.root","ls_2"));
    vector<double> vals_ls_2 = p_ls_2.getVals();

    pull p_ls_3(paraNames,"signal_sys3/pull_*.root","Coherence");
    TMatrixD* cov_ls_3 = new TMatrixD(p_ls_3.getDeltaCov("signal/pull__1.root","ls_3"));
    vector<double> vals_ls_3 = p_ls_3.getVals();

    pull p_ls_4(paraNames,"signal_sys4/pull_*.root","Coherence");
    TMatrixD* cov_ls_4 = new TMatrixD(p_ls_4.getDeltaCov("signal/pull__1.root","ls_4"));
    vector<double> vals_ls_4 = p_ls_4.getVals();

    pull p_ls_5(paraNames,"signal_sys5/pull_*.root","Coherence");
    TMatrixD* cov_ls_5 = new TMatrixD(p_ls_5.getDeltaCov("signal/pull__1.root","ls_5"));
    vector<double> vals_ls_5 = p_ls_5.getVals();

    pull p_ls_6(paraNames,"signal_sys6/pull_*.root","Coherence");
    TMatrixD* cov_ls_6 = new TMatrixD(p_ls_6.getDeltaCov("signal/pull__1.root","ls_6"));
    vector<double> vals_ls_6 = p_ls_6.getVals();

    pull p_ls_7(paraNames,"signal_sys7/pull_*.root","Coherence");
    TMatrixD* cov_ls_7 = new TMatrixD(p_ls_7.getDeltaCov("signal/pull__1.root","ls_7"));
    vector<double> vals_ls_7 = p_ls_7.getVals();

    vector< vector <double> > vec_vals_ls;
    vec_vals_ls.push_back(vals_ls_1);
    vec_vals_ls.push_back(vals_ls_2);
    vec_vals_ls.push_back(vals_ls_3);
    vec_vals_ls.push_back(vals_ls_4);
    vec_vals_ls.push_back(vals_ls_5);
    vec_vals_ls.push_back(vals_ls_6);
    vec_vals_ls.push_back(vals_ls_7);
    TMatrixD cov_ls(p_ls_1.sampleVariance(vec_vals_ls));

    covs.push_back(new TMatrixD(cov_ls));


    /// Resonance parameters systematics
    pull p_rp_1(paraNames,"signal_toy/pull_mass_K1_1270_*.root","Coherence");
    TMatrixD* cov_rp_1 = new TMatrixD(p_rp_1.getDeltaCov("signal_toy/pull__*.root","_mass_K1_1270"));
    pull p_rp_2(paraNames,"signal_toy/pull_width_K1_1270_*.root","Coherence");
    TMatrixD* cov_rp_2 = new TMatrixD(p_rp_2.getDeltaCov("signal_toy/pull__*.root","_width_K1_1270"));

    pull p_rp_3(paraNames,"signal_toy/pull_mass_K_1460_*.root","Coherence");
    TMatrixD* cov_rp_3 = new TMatrixD(p_rp_3.getDeltaCov("signal_toy/pull__*.root","_mass_K_1460"));
    pull p_rp_4(paraNames,"signal_toy/pull_width_K1_1460_*.root","Coherence");
    TMatrixD* cov_rp_4 = new TMatrixD(p_rp_4.getDeltaCov("signal_toy/pull__*.root","_width_K_1460"));

    pull p_rp_5(paraNames,"signal_toy/pull_mass_Ks_*.root","Coherence");
    TMatrixD* cov_rp_5 = new TMatrixD(p_rp_5.getDeltaCov("signal_toy/pull__*.root","_mass_Ks"));
    pull p_rp_6(paraNames,"signal_toy/pull_width_Ks_*.root","Coherence");
    TMatrixD* cov_rp_6 = new TMatrixD(p_rp_6.getDeltaCov("signal_toy/pull__*.root","_width_Ks"));

    pull p_rp_7(paraNames,"signal_toy/pull_mass_rho_*.root","Coherence");
    TMatrixD* cov_rp_7 = new TMatrixD(p_rp_7.getDeltaCov("signal_toy/pull__*.root","_mass_rho"));
    pull p_rp_8(paraNames,"signal_toy/pull_width_rho_*.root","Coherence");
    TMatrixD* cov_rp_8 = new TMatrixD(p_rp_8.getDeltaCov("signal_toy/pull__*.root","_width_rho"));

    pull p_rp_9(paraNames,"signal_toy/pull_mass_K0s_*.root","Coherence");
    TMatrixD* cov_rp_9 = new TMatrixD(p_rp_9.getDeltaCov("signal_toy/pull__*.root","_mass_K0s"));
    pull p_rp_10(paraNames,"signal_toy/pull_width_K0s_*.root","Coherence");
    TMatrixD* cov_rp_10 = new TMatrixD(p_rp_10.getDeltaCov("signal_toy/pull__*.root","_width_K0s"));

    // Add
    TMatrixD cov_rp(*cov_rp_1);
    cov_rp +=  *cov_rp_2 ;
    cov_rp +=  *cov_rp_3 ;
    //cov_rp +=  *cov_rp_4 ;
    cov_rp +=  *cov_rp_5 ;
    cov_rp +=  *cov_rp_6 ;
    cov_rp +=  *cov_rp_7 ;
    cov_rp +=  *cov_rp_8 ;
    cov_rp +=  *cov_rp_9 ;
    cov_rp +=  *cov_rp_10 ;
 
    covs.push_back(new TMatrixD(cov_rp));

    /// Form factor
    pull p_f_1(paraNames,"signal_toy/pull_BW_radius_*.root","Coherence");
    TMatrixD* cov_f_1 = new TMatrixD(p_f_1.getDeltaCov("signal_toy/pull__*.root"));
    covs.push_back(cov_f_1);

    /// Phsp-Acc systematics
    pull p_phsp_acc_1(paraNames,"signal_sys_acc1/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_1 = new TMatrixD(p_phsp_acc_1.getDeltaCov("signal/pull__1.root","phsp_acc_1"));
    vector<double> vals_ps_1 = p_phsp_acc_1.getVals();

    pull p_phsp_acc_2(paraNames,"signal_sys_acc2/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_2 = new TMatrixD(p_phsp_acc_2.getDeltaCov("signal/pull__1.root","phsp_acc_2"));
    vector<double> vals_ps_2 = p_phsp_acc_2.getVals();

    pull p_phsp_acc_3(paraNames,"signal_sys_acc3/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_3 = new TMatrixD(p_phsp_acc_3.getDeltaCov("signal/pull__1.root","phsp_acc_3"));
    vector<double> vals_ps_3 = p_phsp_acc_3.getVals();

    pull p_phsp_acc_4(paraNames,"signal_sys_acc4/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_4 = new TMatrixD(p_phsp_acc_4.getDeltaCov("signal/pull__1.root","phsp_acc_4"));
    vector<double> vals_ps_4 = p_phsp_acc_4.getVals();

    pull p_phsp_acc_5(paraNames,"signal_sys_acc5/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_5 = new TMatrixD(p_phsp_acc_5.getDeltaCov("signal/pull__1.root","phsp_acc_5"));
    vector<double> vals_ps_5 = p_phsp_acc_5.getVals();

    pull p_phsp_acc_6(paraNames,"signal_sys_acc6/pull_*.root","Coherence");
    TMatrixD* cov_phsp_acc_6 = new TMatrixD(p_phsp_acc_6.getDeltaCov("signal/pull__1.root","phsp_acc_6"));
    vector<double> vals_ps_6 = p_phsp_acc_6.getVals();


    vector< vector <double> > vec_vals_ps;
    vec_vals_ps.push_back(vals_ps_1);
    vec_vals_ps.push_back(vals_ps_2);
    vec_vals_ps.push_back(vals_ps_3);
    vec_vals_ps.push_back(vals_ps_4);
    vec_vals_ps.push_back(vals_ps_5);
    vec_vals_ps.push_back(vals_ps_6);

    TMatrixD cov_phsp(p_phsp_acc_1.sampleVariance(vec_vals_ps));
    covs.push_back(new TMatrixD(cov_phsp));

    ///
    vector<double> errs_sys;    
    for(int i =0 ; i <paraNames.size() ; i++){
        double tot = 0.;
        for(int j =0 ; j <covs.size() ; j++){
            tot += (*covs[j])[i][i];
        }
        cout << "sigma_sys(" <<  paraNames[i] << ") = " <<  sqrt(tot) << endl;
    }

}

void crossCheck_signal() {

    /// Fit parameters    
    vector<TString> paraNames;
    paraNames.push_back("mass_K_1__1400_p");
    paraNames.push_back("width_K_1__1400_p");
    paraNames.push_back("mass_Ks_1410_p");
    paraNames.push_back("width_Ks_1410_p");

    paraNames.push_back("r");
    paraNames.push_back("delta");
    paraNames.push_back("gamma");

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
    SummaryFile.open("../../../../../TD-AnaNote/latex/tables/fullFit/signal/crossCheck_table.tex",std::ofstream::trunc);

    for(int i =0 ; i <paraNames.size() ; i++){
        SummaryFile << p_1.latexName(paraNames[i])  << " & " ;
        SummaryFile << std::fixed << std::setprecision(2) << vals[i] << " $\\pm$ " << errs_stat[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(2) << vals1[i] << " $\\pm$ " << errs_stat1[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(2) << vals2[i] << " $\\pm$ " << errs_stat2[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(2) << vals3[i] << " $\\pm$ " << errs_stat3[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(2) << vals4[i] << " $\\pm$ " << errs_stat4[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(2) << vals5[i] << " $\\pm$ " << errs_stat5[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(2) << vals6[i] << " $\\pm$ " << errs_stat6[i] << " & ";

        SummaryFile << std::fixed << std::setprecision(2) << vals7[i] << " $\\pm$ " << errs_stat7[i] << " & ";
        SummaryFile << std::fixed << std::setprecision(2) << vals8[i] << " $\\pm$ " << errs_stat8[i];

        SummaryFile  << " \\\\ " << "\n"; 
    }  
}


int main(int argc, char** argv){

    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    //gStyle->SetOptFit(111);
    //gStyle->UseCurrentStyle();

   //fitParams();
   fractions();

       //altModels();
    // corrCP();
//crossCheck_signal();

    return 0;
}
