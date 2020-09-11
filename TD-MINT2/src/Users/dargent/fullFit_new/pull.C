#define pull_cxx
#include "pull.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;


TString pull::latexName(TString s){
    if(s == "D_bar") return "$\\bar{D}$";
    if(s == "S_bar") return "$\\bar{S}$";
    if(s == "delta") return "$\\delta$";
    if(s == "gamma") return "$\\gamma - 2 \\beta_{s}$";

    if(s == "mass_K_1__1400_p") return "$m_{K_1(1400)} $";
    if(s == "width_K_1__1400_p") return "$\\Gamma_{K_1(1400)}$";
    if(s == "mass_Ks_1410_p") return "$m_{K^{*}(1410)}$";
    if(s == "width_Ks_1410_p") return "$\\Gamma_{K^{*}(1410)}$";

    if(s == "a_K1_1400_Amp") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) )  \\, \\text{Mag} (b \\to c)$";
    if(s == "a_K1_1400_Phase") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) ) \\, \\text{Phase} (b \\to c)$";
    if(s == "abar_K1_1400_Amp") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) ) \\, \\text{Mag} (b \\to u)$";
    if(s == "abar_K1_1400_Phase") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) ) \\, \\text{Phase} (b \\to u)$";
    if(s == "Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Amp") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}(892) \\, \\pi ) \\, \\text{Mag}$";
    if(s == "Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Phase") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}(892) \\, \\pi ) \\, \\text{Phase}$";
    if(s == "Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Amp") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi ) \\, \\text{Mag} $";
    if(s == "Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Phase") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi ) \\, \\text{Phase} $";

    if(s == "a_Ks_1410_Amp") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi ) \\, \\text{Mag} (b \\to c)$";
    if(s == "a_Ks_1410_Phase") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi ) \\, \\text{Phase} (b \\to c)$";
    if(s == "abar_Ks_1410_Amp") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi ) \\, \\text{Mag} (b \\to u)$";
    if(s == "abar_Ks_1410_Phase") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi ) \\, \\text{Phase} (b \\to u)$";
    if(s == "Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Amp") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K \\, \\rho(770) ) \\, \\text{Mag}$";
    if(s == "Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Phase") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K \\, \\rho(770) ) \\, \\text{Phase}$";

    if(s == "a_K_1460_Amp") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi ) \\, \\text{Mag} (b \\to c)$";
    if(s == "a_K_1460_Phase") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi ) \\, \\text{Phase} (b \\to c)$";
    if(s == "abar_K_1460_Amp") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi ) \\, \\text{Mag} (b \\to u)$";
    if(s == "abar_K_1460_Phase") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi ) \\, \\text{Phase} (b \\to u)$";

    if(s == "a_NS_Ks_Amp") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) \\, \\text{Mag} (b \\to c)$";
    if(s == "a_NS_Ks_Phase") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) \\, \\text{Phase} (b \\to c)$";
    if(s == "abar_NS_Ks_Amp") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) \\, \\text{Mag} (b \\to u)$";
    if(s == "abar_NS_Ks_Phase") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) \\, \\text{Phase} (b \\to u)$";

    if(s == "a_NS_rho_Amp") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) \\, \\text{Mag} (b \\to c)$";
    if(s == "a_NS_rho_Phase") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) \\, \\text{Phase} (b \\to c)$";
    if(s == "abar_NS_rho_Amp") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) \\, \\text{Mag} (b \\to u)$";
    if(s == "abar_NS_rho_Phase") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) \\, \\text{Phase} (b \\to u)$";

    if(s == "Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1270_p_D___Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270)[D] \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1__1270_p__rho_1450_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(1450) )$";
    if(s == "Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi )$";
    if(s == "Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1400_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_2_s_1430_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_2_s_1430_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1460_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1460_p__sigma10__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\sigma )$";
    if(s == "Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_Ks_1680_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_Ks_1680_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_2__1770_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_2__1770_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_NonResS0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{S} \\, \\, K^{*}(892)$";
    if(s == "Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_P__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[P] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_D__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[D] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_NonResS0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\rho(770)$";
    if(s == "Bs0_NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_P__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[P] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_D__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[D] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_NonResS0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\sigma$";
    if(s == "Bs0_NonResV0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\sigma$";
    if(s == "Bs0_NonResS0__DsmKp_f_0__980_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(980)$";
    if(s == "Bs0_f_2__1270_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_2(1270)$";
    if(s == "Bs0_f_2__1270_0__pippim_NonResV0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, f_2(1270)$";
    if(s == "Bs0_f_0__1370_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(1370)$";
    if(s == "Bs0_K_0_s_1430_0__Kppim_NonResS0__Dsmpip_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}_{0}(1430)$";

    if(s == "bar_Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_1__1270_p_D___Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270)[D] \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_K_1__1270_p__rho_1450_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(1450) )$";
    if(s == "bar_Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi )$";
    if(s == "bar_Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_1__1400_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_K_2_s_1430_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_2_s_1430_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_1460_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_K_1460_p__sigma10__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\sigma )$";
    if(s == "bar_Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_Ks_1680_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_Ks_1680_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_K_2__1770_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K^{*}(892) \\, \\pi )$";
    if(s == "bar_Bs0_K_2__1770_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K \\, \\rho(770) )$";
    if(s == "bar_Bs0_NonResS0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{S} \\, \\, K^{*}(892)$";
    if(s == "bar_Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "bar_Bs0_P__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[P] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "bar_Bs0_D__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[D] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "bar_Bs0_NonResS0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\rho(770)$";
    if(s == "bar_Bs0_NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "bar_Bs0_P__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[P] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "bar_Bs0_D__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[D] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "bar_Bs0_NonResS0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\sigma$";
    if(s == "bar_Bs0_NonResV0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\sigma$";
    if(s == "bar_Bs0_NonResS0__DsmKp_f_0__980_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(980)$";
    if(s == "bar_Bs0_f_2__1270_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_2(1270)$";
    if(s == "bar_Bs0_f_2__1270_0__pippim_NonResV0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, f_2(1270)$";
    if(s == "bar_Bs0_f_0__1370_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(1370)$";
    if(s == "bar_Bs0_K_0_s_1430_0__Kppim_NonResS0__Dsmpip_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}_{0}(1430)$";

    if(s == "Sum") return "$\\text{Sum}$";
    if(s == "bar_Sum") return "$\\text{Sum}$";

    if(s == "n2ll") return "$\\text{-2NLL}$";

    return "$" + s + "$";
}

TString pull::latexNameMod(TString s){

    if(s == "Bs0toK_1__1270_p_toKs_892_0_toKp_pim__pip__Dsm_Amp") return "$\\phantom{B_s \\to D_s \\, (} K_1(1270) \\to K^{*}(892) \\, \\pi \\phantom{)} $";
    if(s == "Bs0toK_1__1270_p_toK_0_s_1430_0_toKp_pim__pip__Dsm_Amp") return "$\\phantom{B_s \\to D_s \\, (} K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi \\phantom{)} $";
    if(s == "a_K1_1400_Amp") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) ) $";
    if(s == "a_Ks_1410_Amp") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi ) $";
    if(s == "Bs0toKs_1410_p_torho_770_0_topip_pim__Kp__Dsm_Amp") return "$\\phantom{B_s \\to D_s \\, (} K^{*}(1410) \\to K \\, \\rho(770) \\phantom{)} $";
    if(s == "abar_K_1460_Amp") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi ) $";
    if(s == "a_NS_Ks_Amp") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) $";
    if(s == "abar_NS_Ks_Amp") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892) $";
    if(s == "a_NS_rho_Amp") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) $";
    if(s == "abar_NS_rho_Amp") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770) $";

    if(s == "k") return "$\\kappa$";
    if(s == "delta") return "$\\delta \\, [\\degrees]$";
    if(s == "gamma") return "$\\gamma - 2 \\beta_{s} \\, [\\degrees]$";

    if(s == "mass_K_1__1400_p") return "$m_{K_1(1400)} \\, [\\text{MeV}]$";
    if(s == "width_K_1__1400_p") return "$\\Gamma_{K_1(1400)} \\, [\\text{MeV}]$";
    if(s == "mass_Ks_1410_p") return "$m_{K^{*}(1410)} \\, [\\text{MeV}]$";
    if(s == "width_Ks_1410_p") return "$\\Gamma_{K^{*}(1410)} \\, [\\text{MeV}]$";


    else return latexName(s);
}

TMatrixD pull::getAbsDiff(TMatrixD cov1,TMatrixD cov2){

	TMatrixD cov_diff(_paraNames.size(),_paraNames.size());

	for (int i = 0 ; i < _paraNames.size(); i++)
        	for (int j = 0 ; j < _paraNames.size(); j++){
			if(i==j)cov_diff[i][j] = abs(cov1[i][j]-cov2[i][j]);
			else cov_diff[i][j] = 0.;
		}

	return cov_diff;
}

TMatrixD pull::combineCov_maxVal(vector<TMatrixD*> vec){

	TMatrixD m(*vec[0]);

	for (int i = 0 ; i <  m.GetNcols(); i++) {
		for (int j = 0 ; j < m.GetNcols(); j++) {
			double max = 0.;
			for(int n = 0; n < vec.size(); n++){
				 if(abs((*vec[n])(i,j)) > abs(max)) max = (*vec[n])(i,j);
			}
			m(i,j) = max; 
		}
	}

	return m;
}

vector<double> pull::sampleMean(vector< vector<double> > vec_vals){

	vector<double> means;

	for(int i=0; i < vec_vals[0].size(); i++){
			double sum = 0.;
			for(int n = 0; n < vec_vals.size(); n++){
				sum += vec_vals[n][i];
			}
			means.push_back(sum/vec_vals.size());		
	}
	return means;
}


TMatrixD pull::sampleVariance(vector< vector<double> > vec_vals, int removeOutlier){

	TMatrixD m(vec_vals[0].size(),vec_vals[0].size());

	vector<double> vec_avg;
	vector<int> vec_remove;

	for(int i=0; i < vec_vals[0].size(); i++){
		double avg = 0;
		for(int n = 0; n < vec_vals.size(); n++){
				avg += vec_vals[n][i];
		}
		avg /= vec_vals.size();
		vec_avg.push_back(avg);
	}

	for(int i=0; i < vec_vals[0].size(); i++){
		int pos = 0;
		for(int n = 1; n < vec_vals.size(); n++){
				if(abs(vec_avg[i]-vec_vals[n][i])>abs(vec_avg[i]-vec_vals[pos][i]))pos = n;
		}
		vec_remove.push_back(pos);
	}

	if(removeOutlier == 1){
		for(int i=0; i < vec_vals[0].size(); i++){
			double sum = 0.;
			double sum_square = 0.;
			for(int n = 0; n < vec_vals.size(); n++){
				if(n==vec_remove[i])continue;
				sum += vec_vals[n][i];
				sum_square += pow(vec_vals[n][i],2);	
			}
			m(i,i) = (sum_square - pow(sum,2)/(vec_vals.size()-1.))/(vec_vals.size()-2.);		
		}
	}
	else if(removeOutlier > 1){
		for(int i=0; i < vec_vals[0].size(); i++){
			double sum = 0.;
			double sum_square = 0.;
			for(int n = 0; n < vec_vals.size(); n++){
				if(i < removeOutlier && n==vec_remove[i])continue;
				sum += vec_vals[n][i];
				sum_square += pow(vec_vals[n][i],2);	
			}
			if(i < removeOutlier)m(i,i) = (sum_square - pow(sum,2)/(vec_vals.size()-1.))/(vec_vals.size()-2.);	
			else m(i,i) = (sum_square - pow(sum,2)/vec_vals.size())/(vec_vals.size()-1.);	
		}
	}
	else for(int i=0; i < vec_vals[0].size(); i++){
			double sum = 0.;
			double sum_square = 0.;
			for(int n = 0; n < vec_vals.size(); n++){
				sum += vec_vals[n][i];
				sum_square += pow(vec_vals[n][i],2);	
			}
			m(i,i) = (sum_square - pow(sum,2)/vec_vals.size())/(vec_vals.size()-1.);		
	}

	return m;
}

vector<double> pull::sampleMean(){

   if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
   }
   int N = fChain->GetEntries();

   vector<double> sum(_paraNames.size(),0.);
   vector<double> sum_square(_paraNames.size(),0.);

   vector<double> var(_paraNames.size(),0.);
 
   for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
	for (int i = 0 ; i < _paraNames.size(); i++) {
			sum[i] += *_means[i];
			sum_square[i] += pow(*_means[i],2);	
	}
    }

    for (int i = 0 ; i < _paraNames.size(); i++) {
		var[i] = sum[i]/(double)N;
    }

    return var;
}

vector<double> pull::sampleSigma(){

    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    int N = fChain->GetEntries();

   vector<double> sum(_paraNames.size(),0.);
   vector<double> sum_square(_paraNames.size(),0.);

   vector<double> var(_paraNames.size(),0.);
 
   for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
	for (int i = 0 ; i < _paraNames.size(); i++) {
			sum[i] += *_means[i];
			sum_square[i] += pow(*_means[i],2);	
	}
    }

    for (int i = 0 ; i < _paraNames.size(); i++) {
		var[i] = sqrt((sum_square[i] - pow(sum[i],2)/(double)N)/((double)N-.1));
    }

    return var;
}


vector <double> pull::getPullMean(){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }

    TMatrixD cov(_paraNames.size(),_paraNames.size());
    int N = fChain->GetEntries();
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -3.,3.));
        
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        for (int i = 0 ; i < _paraNames.size(); i++){
            h_pulls[i]->Fill(*_pulls[i]);
        }
    }

    TF1 *gaussian = new TF1("gaussian","gaus",-5.,5.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
        
    vector<double> fit_means,fit_sigmas;
    for(int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(h_pulls[i]->GetMean());
        fit_sigmas.push_back(h_pulls[i]->GetStdDev());
//         fit_means.push_back(gaussian->GetParameter(1));
//         fit_sigmas.push_back(gaussian->GetParameter(2));
    }

    return fit_means;
}

vector <double> pull::getPullSigma(){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }

    TMatrixD cov(_paraNames.size(),_paraNames.size());
    int N = fChain->GetEntries();
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -3.,3.));
        
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        for (int i = 0 ; i < _paraNames.size(); i++){
            h_pulls[i]->Fill(*_pulls[i]);
        }
    }

    TF1 *gaussian = new TF1("gaussian","gaus",-5.,5.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
        
    vector<double> fit_means,fit_sigmas;
    for(int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(h_pulls[i]->GetMean());
        fit_sigmas.push_back(h_pulls[i]->GetStdDev());
//         fit_means.push_back(gaussian->GetParameter(1));
//         fit_sigmas.push_back(gaussian->GetParameter(2));
    }

    return fit_sigmas;
}



vector <double> pull::getVals(){
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    fChain->GetEntry(0);  
    vector<double> vals;
    for (int i = 0 ; i < _paraNames.size(); i++) 
	vals.push_back(*_means[i]);
    
    return vals;
}   

vector <double> pull::getErrs(){
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    fChain->GetEntry(0);  
    vector<double> vals;
    for (int i = 0 ; i < _paraNames.size(); i++) 
	vals.push_back(*_errs[i]);
    
    return vals;
}   


TMatrixD pull::getStatCov(TString label){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }
    int N = fChain->GetEntries();
          
   vector<double> inits(_paraNames.size(),0.); 
   for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
	for (int i = 0 ; i < _paraNames.size(); i++) {
			inits[i] += *_means[i]/(double)N;
	}
    }

    TMatrixD cov(_paraNames.size(),_paraNames.size());
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
	for (int i = 0 ; i < _paraNames.size(); i++) {
 		for (int j = 0 ; j < _paraNames.size(); j++) {
//  			cov[i][j] += (*_means[i] - *_inits[i]) * (*_means[j] - *_inits[j])/(N-1.);
 			cov[i][j] += (*_means[i] - inits[i]) * (*_means[j] - inits[j])/((double)N-1.);
		}
	}
    }

    return cov;
}

TMatrixD pull::getCov(TString label, double pull_max , bool scaleCov){
    
    if(fChain == 0){
            cout << "ERROR:: No file found" << endl;
            throw "ERROR";
    }

    TMatrixD cov(_paraNames.size(),_paraNames.size());
    int N = fChain->GetEntries();

    if(N == 1){
        fChain->GetEntry(0);  
        for (int i = 0 ; i < _paraNames.size(); i++) 
            for (int j = 0 ; j < _paraNames.size(); j++) {
                if(i==j)cov[i][j] = pow((*_means[i] - *_inits[i]),2);
                else cov[i][j] = 0.;
            }
        return cov;
    }
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -5.,5.));
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        for (int i = 0 ; i < _paraNames.size(); i++){
	    if(pull_max > 0 && abs(*_pulls[i]) > pull_max)continue;
            h_pulls[i]->Fill(*_pulls[i]);
            for (int j = 0 ; j < _paraNames.size(); j++) {
                cov[i][j] += (*_means[i] - *_inits[i]) * (*_means[j] - *_inits[j])/(N-1.);
            }
        }
    }
    
    TCanvas* c = new TCanvas();
    TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
        
    ofstream SummaryFile;
    SummaryFile.open("pull_results/pull_table"+label+".tex",std::ofstream::trunc);
    SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    
    vector<double> fit_means,fit_sigmas;
    for(int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(gaussian->GetParameter(1));
        fit_sigmas.push_back(gaussian->GetParameter(2));

        SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
        h_pulls[i]->Draw("");
        gaussian->Draw("SAME");
        c->Print("pull_results/pull_"+ _paraNames[i] + label + ".eps");
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";

    TMatrixD cov_prime(cov);
    for (int i = 0 ; i < _paraNames.size(); i++)
        for (int j = 0 ; j < _paraNames.size(); j++){
            //cov_prime[i][j] = cov[i][j]*abs(fit_means[i])*abs(fit_means[j]);
//             cov_prime[i][j] = fit_means[i]* fit_means[j];
            if(scaleCov)cov_prime[i][j] = cov[i][j]/sqrt(cov[i][i])/sqrt(cov[j][j])*abs(fit_means[i])*abs(fit_means[j]);
// 	    else cov_prime[i][j] = fit_means[i]* fit_means[j];
	    else cov[i][j]/sqrt(cov[i][i])/sqrt(cov[j][j]);
	}

    return cov_prime;
}

TMatrixD pull::getDeltaCov(TString refFileName,TString label, double pull_max , bool scaleCov){
    
    int N = fChain->GetEntries();

    TChain* chain =  new TChain(_treeName);
    if(_removeBar)refFileName.ReplaceAll("_Bar","");
    TString lastFile;

    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = refFileName;
	file.ReplaceAll("*",index.str()); 

	if(std::ifstream(((string)file).c_str()).good()){ 
		chain->Add(file);
		lastFile = file;
	}
	else {
		chain->Add(lastFile);
		_skip.push_back(i);
	} 
    }    
    else chain->Add(refFileName); 

    if(fChain == 0 || chain == 0){
        cout << "ERROR:: No file found" << endl;
        throw "ERROR";
    }

    if(N > chain->GetEntries()){
        cout << "ERROR:: Inconsistent number of entries" << endl;
        throw "ERROR";
    }
    
    vector<double*> means;
    vector<double*> inits;
    vector<double*> errs;
    vector<double*> pulls;
    
    for (int i = 0 ; i < _paraNames.size(); i++) {
        double * mean = new double[1];
        means.push_back(mean);
	if(_fraction){
		chain->SetBranchAddress(_paraNames[i], mean);
		continue;
	}
        chain->SetBranchAddress(_paraNames[i]+"_mean", mean);
        
        double * init = new double[1];
        inits.push_back(init);
        chain->SetBranchAddress(_paraNames[i]+"_init", init);
        
        double * err = new double[1];
        errs.push_back(err);
        chain->SetBranchAddress(_paraNames[i]+"_err", err);
        
        double * pull = new double[1];
        pulls.push_back(pull);
        chain->SetBranchAddress(_paraNames[i]+"_pull", pull);
    } 
    
    TMatrixD cov(_paraNames.size(),_paraNames.size());

    if(N == 1){
        fChain->GetEntry(0); 
        chain->GetEntry(0); 
        for (int i = 0 ; i < _paraNames.size(); i++) 
            for (int j = 0 ; j < _paraNames.size(); j++) {
                if(i==j)cov[i][j] = pow((*_means[i] - *means[i]),2);
                else cov[i][j] = 0.;
            }
        return cov;
    }

    if(_fraction){
	for (int n=0; n <N ;n++) {
		fChain->GetEntry(n);  
		chain->GetEntry(n);
		for (int i = 0 ; i < _paraNames.size(); i++)
		for (int j = 0 ; j < _paraNames.size(); j++) 
			cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
	}
	return cov;
    }

    double max = 0.;
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++){
	   double val = (*_means[i]-*means[i])/(*errs[i]);
           if(abs(val) > max) max = abs(val);
	}
    }
    //max *= 1.;
    max = max>5 ? 5: max;  
    
    vector<TH1D*> h_pulls;
    for (int i = 0 ; i < _paraNames.size(); i++) 
        h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 50, -pull_max,pull_max));
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++){
	    if(abs(*_pulls[i])>pull_max || abs(*pulls[i])>pull_max) continue;
            h_pulls[i]->Fill((*_means[i]-*means[i])/(*errs[i]));
	}
    }
    
    TCanvas* c = new TCanvas();
    TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
    gaussian->SetParameters(1.,0.,1.);
    gaussian->SetParLimits(1,-1., 1.);
    gaussian->SetParLimits(2, 0., 2.);
    gaussian->SetLineColor(kRed);
    
    ofstream SummaryFile;
    SummaryFile.open("pull_results/pull_table"+label+".tex",std::ofstream::trunc);
    SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\hline" << "\n";
    
    vector<double> fit_means,fit_sigmas;
    for (int i = 0 ; i < _paraNames.size(); i++) {
        h_pulls[i]->Fit(gaussian);
        fit_means.push_back(gaussian->GetParameter(1));
        fit_sigmas.push_back(gaussian->GetParameter(2));
        SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
        h_pulls[i]->Draw("");
        gaussian->Draw("SAME");
        c->Print("pull_results/pull_"+ _paraNames[i] + label + ".eps");
    }
    
    SummaryFile << "\\hline" << "\n";
    SummaryFile << "\\end{tabular}" << "\n";
    
    for (int n=0; n <N ;n++) {
        fChain->GetEntry(n);  
        chain->GetEntry(n);
        for (int i = 0 ; i < _paraNames.size(); i++)
            for (int j = 0 ; j < _paraNames.size(); j++) {
	     	if(abs(*_pulls[i])>pull_max || abs(*pulls[i])>pull_max) continue;
	     	if(abs(*_pulls[j])>pull_max || abs(*pulls[j])>pull_max) continue;
                cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
	    }
    }
    
    TMatrixD cov_prime(cov);
    for (int i = 0 ; i < _paraNames.size(); i++)
        for (int j = 0 ; j < _paraNames.size(); j++) 
            cov_prime[i][j] = cov[i][j]*sqrt(pow(fit_means[i],2)+pow(fit_sigmas[i],2))*sqrt(pow(fit_means[j],2)+pow(fit_sigmas[j],2))/fit_sigmas[i]/fit_sigmas[j];
    
     if(scaleCov)return cov_prime;
    else return cov;
}

TMatrixD pull::getDeltaCovChol(TString refFileName,TString label,int varPerParChol){
    
    int N = fChain->GetEntries();
    TChain* chain =  new TChain(_treeName);
    TString lastFile;

    if(N>1)for(int i = 1; i <= N; i++){
	stringstream index;
	index << i;
	TString file = refFileName;
	file.ReplaceAll("*",index.str()); 

	if(std::ifstream(((string)file).c_str()).good()){ 
		chain->Add(file);
		lastFile = file;
	}
	else {
		chain->Add(lastFile);
		_skip.push_back(i);
	} 
    }    
    else chain->Add(refFileName); 

    if(fChain == 0 || chain == 0){
        cout << "ERROR:: No file found" << endl;
        throw "ERROR";
    }

    if(N > chain->GetEntries()){
        cout << "ERROR:: Inconsistent number of entries" << endl;
	cout << N << " > " <<  chain->GetEntries() << endl ;
        throw "ERROR";
    }
    
    vector<double*> means;
    vector<double*> inits;
    vector<double*> errs;
    vector<double*> pulls;
    
    for (int i = 0 ; i < _paraNames.size(); i++) {
        double * mean = new double[1];
        means.push_back(mean);
        chain->SetBranchAddress(_paraNames[i]+"_mean", mean);
        
        double * init = new double[1];
        inits.push_back(init);
        chain->SetBranchAddress(_paraNames[i]+"_init", init);
        
        double * err = new double[1];
        errs.push_back(err);
        chain->SetBranchAddress(_paraNames[i]+"_err", err);
        
        double * pull = new double[1];
        pulls.push_back(pull);
        chain->SetBranchAddress(_paraNames[i]+"_pull", pull);
    } 
    
    int N_chol = N/varPerParChol;
    TMatrixD cov_tot(_paraNames.size(),_paraNames.size());

   if(_fraction){
	TMatrixD cov(_paraNames.size(),_paraNames.size());	
	for(int ic = 0 ; ic < N_chol; ic ++ ){
		for (int n= ic*varPerParChol; n < (ic+1)*varPerParChol ;n++) {
			vector<int>::iterator it = find(_skip.begin(),_skip.end(),n);
			if(it != _skip.end())continue;
		
			fChain->GetEntry(n);  
			chain->GetEntry(n);
			for (int i = 0 ; i < _paraNames.size(); i++)
				for (int j = 0 ; j < _paraNames.size(); j++) 
					cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
		}
	}
	return cov;
    }
   
    for(int ic = 0 ; ic < N_chol; ic ++ ){
    
        TMatrixD cov(_paraNames.size(),_paraNames.size());

        vector<TH1D*> h_pulls;
        for (int i = 0 ; i < _paraNames.size(); i++) 
            h_pulls.push_back(new TH1D("pull_"+_paraNames[i],"; Pull " + _paraNames[i] + "; Toy experiments", 40, -1,1));
        
        for (int n= ic*varPerParChol; n < (ic+1)*varPerParChol ;n++) {
	    vector<int>::iterator it = find(_skip.begin(),_skip.end(),n);
	    if(it != _skip.end())continue;

            fChain->GetEntry(n);  
            chain->GetEntry(n);

            for (int i = 0 ; i < _paraNames.size(); i++)
                h_pulls[i]->Fill((*_means[i]-*means[i])/(*errs[i]));
        }
        
        TCanvas* c = new TCanvas();
        TF1 *gaussian = new TF1("gaussian","gaus",-3.,3.);
        gaussian->SetParameters(1.,0.,1.);
        gaussian->SetParLimits(1,-1., 1.);
        gaussian->SetParLimits(2, 0., 2.);
        gaussian->SetLineColor(kRed);
        
        ofstream SummaryFile;
	stringstream number;
	number << ic;
        SummaryFile.open("pull_results/pull_table"+label+"_par_"+number.str()+".tex",std::ofstream::trunc);
        SummaryFile << "\\begin{tabular}{l  c  c}" << "\n";
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "Parameter & $\\mu$ of pull distribution & $\\sigma$ of pull distribution \\\\" << "\n";
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "\\hline" << "\n";
        
        vector<double> fit_means,fit_sigmas;
        
        for (int i = 0 ; i < _paraNames.size(); i++) {
            h_pulls[i]->Fit(gaussian);
            fit_means.push_back(gaussian->GetParameter(1));
            fit_sigmas.push_back(gaussian->GetParameter(2));
            SummaryFile << std::fixed << std::setprecision(2) << latexName(_paraNames[i]) << " & " << gaussian->GetParameter(1) << " $\\pm$ " << gaussian->GetParError(1) <<" & " << gaussian->GetParameter(2) << " $\\pm$ " <<  gaussian->GetParError(2) << " \\\\" << "\n";
            h_pulls[i]->Draw("");
            gaussian->Draw("SAME");
            c->Print("pull_results/pull_"+ _paraNames[i] + label + "_par_"+number.str()+ ".eps");
        }
        
        SummaryFile << "\\hline" << "\n";
        SummaryFile << "\\end{tabular}" << "\n";
        
        for (int n= ic*varPerParChol; n < (ic+1)*varPerParChol ;n++) {
	    vector<int>::iterator it = find(_skip.begin(),_skip.end(),n);
	    if(it != _skip.end())continue;

            fChain->GetEntry(n);  
            chain->GetEntry(n);
            for (int i = 0 ; i < _paraNames.size(); i++)
                for (int j = 0 ; j < _paraNames.size(); j++) 
                    cov[i][j] += (*_means[i] - *means[i]) * (*_means[j] - *means[j])/(N-1.);
        }
        
        TMatrixD cov_prime(cov);
        for (int i = 0 ; i < _paraNames.size(); i++)
            for (int j = 0 ; j < _paraNames.size(); j++) 
                cov_prime[i][j] = cov[i][j]*sqrt(pow(fit_means[i],2)+pow(fit_sigmas[i],2))*sqrt(pow(fit_means[j],2)+pow(fit_sigmas[j],2))/fit_sigmas[i]/fit_sigmas[j];
        
        cov_tot += cov_prime;
    }    
       
    return cov_tot;
}
