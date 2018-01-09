#ifndef TIMEPDFMASTER_HH
#define TIMEPDFMASTER_HH
// author: Philippe d'Argent

#include <complex>
#include <vector>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>

#include "Mint/IReturnRealForEvent.h"
#include "Mint/IReturnComplexForEvent.h"
#include "Mint/DalitzEvent.h"
#include "Mint/FitParRef.h"
#include "Mint/FitParDependent.h"
#include "Mint/IFitParRegister.h"
//#include "Mint/CachedByEvent.h"

#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooUniform.h"

#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/RooCubicSplineFun.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "Mint/TimePdfIntegrator.h"


enum basisType { 
      noBasis=0  ,  expBasis= 3
    , sinBasis=13,  cosBasis=23
    , sinhBasis=63, coshBasis=53 };

using namespace std;
using namespace RooFit ;
using namespace MINT;

class TimePdfMaster
{
 protected:
    
    // Fit parameters
    MINT::FitParameter _tau;
    MINT::FitParameter _dGamma;
    MINT::FitParameter _dm;
    
    MINT::FitParameter _offset_sigma_dt;    
    MINT::FitParameter _scale_mean_dt;
    MINT::FitParameter _scale_sigma_dt;
    
    MINT::FitParameter _c0;
    MINT::FitParameter _c1;
    MINT::FitParameter _c2;
    MINT::FitParameter _c3;
    MINT::FitParameter _c4;
    MINT::FitParameter _c5;
    MINT::FitParameter _c6;
    MINT::FitParameter _c7;
    MINT::FitParameter _c8;
    MINT::FitParameter _c9;
    
    MINT::FitParameter _p0_os;
    MINT::FitParameter _p1_os;
    MINT::FitParameter _delta_p0_os;
    MINT::FitParameter _delta_p1_os;
    MINT::FitParameter _avg_eta_os;
    MINT::FitParameter _tageff_os;
    MINT::FitParameter _tageff_asym_os;
    
    MINT::FitParameter _p0_ss;
    MINT::FitParameter _p1_ss;
    MINT::FitParameter _delta_p0_ss;
    MINT::FitParameter _delta_p1_ss;
    MINT::FitParameter _avg_eta_ss;
    MINT::FitParameter _tageff_ss;
    MINT::FitParameter _tageff_asym_ss;
    
    MINT::FitParameter _production_asym;
    MINT::FitParameter _detection_asym;

    /// Cast of MINT parameters to B2DX parameters
    RooRealVar* _r_t;
    RooRealVar* _r_dt;
    RooRealVar* _r_dt_scaled;
    RooCategory* _r_q_OS;
    RooRealVar* _r_eta_OS;
    RooCategory* _r_q_SS;
    RooRealVar* _r_eta_SS;
    RooCategory* _r_f;

    RooRealVar* _r_tau;
    RooRealVar* _r_dGamma;
    RooRealVar* _r_dm;
    
    RooRealVar* _r_scale_mean_dt;
    RooRealVar* _r_offset_sigma_dt;
    RooRealVar* _r_scale_sigma_dt;
    
    RooRealVar* _r_c0;
    RooRealVar* _r_c1;
    RooRealVar* _r_c2;
    RooRealVar* _r_c3;
    RooRealVar* _r_c4;
    RooRealVar* _r_c5;
    RooRealVar* _r_c6;
    RooRealVar* _r_c7;
    RooRealVar* _r_c8;
    RooRealVar* _r_c9;
    RooFormulaVar* _coeff_last;

    RooRealVar* _r_p0_os;
    RooRealVar* _r_p1_os;
    RooRealVar* _r_delta_p0_os;
    RooRealVar* _r_delta_p1_os;
    RooRealVar* _r_avg_eta_os;
    RooRealVar* _r_tageff_os;
    RooRealVar* _r_tageff_asym_os;
    RooRealVar* _r_p0_ss;
    RooRealVar* _r_p1_ss;
    RooRealVar* _r_delta_p0_ss;
    RooRealVar* _r_delta_p1_ss;
    RooRealVar* _r_avg_eta_ss;
    RooRealVar* _r_tageff_ss;
    RooRealVar* _r_tageff_asym_ss;
    
    RooRealVar* _r_production_asym;
    RooRealVar* _r_detection_asym;

    // Acceptance function
    RooCubicSplineFun* _spline;
    RooGaussEfficiencyModel* _efficiency;
    
    // CP coefficients
    RooRealVar* _r_C;
    RooRealVar* _r_C_bar;
    RooRealVar* _r_D;
    RooRealVar* _r_D_bar;
    RooRealVar* _r_S;
    RooRealVar* _r_S_bar;
    
    // Decay rate coefficients
    DecRateCoeff_Bd* _cos_coeff;
    DecRateCoeff_Bd* _cosh_coeff;
    DecRateCoeff_Bd* _sin_coeff;
    DecRateCoeff_Bd* _sinh_coeff;
    
    // Time pdf integrators
    TimePdfIntegrator* _cosh_term;
    TimePdfIntegrator* _sinh_term; 
    TimePdfIntegrator* _cos_term; 
    TimePdfIntegrator* _sin_term;
    
    // Marginal pdfs
    TH1D* _h_dt;
    RooDataHist* _r_h_dt;
    RooAbsPdf* _pdf_sigma_t;
    
    RooAbsPdf* _pdf_eta_OS;
    TH1D* _h_eta_OS;
    RooDataHist* _r_h_eta_OS;
    RooAbsPdf* _pdf_eta_OS_uniform;
    
    RooAbsPdf* _pdf_eta_SS;
    TH1D* _h_eta_SS;
    RooDataHist* _r_h_eta_SS;
    RooAbsPdf* _pdf_eta_SS_uniform;

    // Limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;

 public:
  TimePdfMaster(): 
        _tau("tau",2,1.509,0.1),
        _dGamma("dGamma",2,0.09,0.1),
        _dm("dm",2,17.757,0.1),
        _scale_mean_dt("scale_mean_dt",1,1,0.1),
        _offset_sigma_dt("offset_sigma_dt",1,0.,0.1),
        _scale_sigma_dt("scale_sigma_dt",1,1.2,0.1),
        _c0("c0",1,1,0.1),
        _c1("c1",1,1,0.1),
        _c2("c2",1,1,0.1),
        _c3("c3",1,1,0.1),
        _c4("c4",1,1,0.1),
        _c5("c5",1,1,0.1),
        _c6("c6",1,1,0.1),
        _c7("c7",1,1,0.1),
        _c8("c8",1,1,0.1),
        _c9("c9",1,1,0.1),
        _p0_os("p0_os",2,0.,0.),
        _p1_os("p1_os",2,0.,0.),
        _delta_p0_os("delta_p0_os",1,0.,0.),
        _delta_p1_os("delta_p1_os",1,0.,0.),
        _avg_eta_os("avg_eta_os",2,0.,0.),
        _tageff_os("tageff_os",2,0.,0.),
        _tageff_asym_os("tageff_asym_os",1,0.,0.),
        _p0_ss("p0_ss",2,0.,0.),
        _p1_ss("p1_ss",2,0.,0.),
        _delta_p0_ss("delta_p0_ss",1,0.,0.),
        _delta_p1_ss("delta_p1_ss",1,0.,0.),
        _avg_eta_ss("avg_eta_ss",2,0.,0.),
        _tageff_ss("tageff_ss",2,0.,0.),
        _tageff_asym_ss("tageff_asym_ss",1,0.,0.),
        _production_asym("production_asym",1,0.,0.),
        _detection_asym("detection_asym",1,0.,0.),
        _min_TAU("min_TAU", 0.4),
        _max_TAU("max_TAU", 10.)
    {
        /// Init B2DX parameters
        
        // Observables
        _r_t = new RooRealVar("t", "time", _min_TAU, _max_TAU);
        _r_dt = new RooRealVar("dt", "per-candidate time resolution estimate",0., 0.1);
        _r_f = new RooCategory("qf", "qf");
        _r_f->defineType("h+", +1);
        _r_f->defineType("h-", -1);
        _r_f->setRange("Range_p","h+");
        _r_f->setRange("Range_m","h-");
        _r_q_OS = new RooCategory("q_OS", "q_OS");
        _r_q_OS->defineType("B+", +1);
        _r_q_OS->defineType("B-", -1) ;
        _r_q_OS->defineType("untagged", 0);
        _r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5);
        _r_q_SS = new RooCategory("q_SS", "q_SS");
        _r_q_SS->defineType("B+", +1);
        _r_q_SS->defineType("B-", -1) ;
        _r_q_SS->defineType("untagged", 0);
        _r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5);
        
        // Fit parameters
        _r_tau = new RooRealVar("tau", "tau",_tau);
        _r_dGamma = new RooRealVar("dGamma", "dGamma",_dGamma);
        _r_dm = new RooRealVar("dm", "dm",_dm);
        
        _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", _scale_mean_dt);
        _r_offset_sigma_dt = new RooRealVar("offset_sigma_dt", "offset_sigma_dt", _offset_sigma_dt);
        _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", _scale_sigma_dt);
        
        _r_c0 = new RooRealVar("coeff_0", "coeff_0",_c0);
        _r_c1 = new RooRealVar("coeff_1", "coeff_1",_c1);
        _r_c2 = new RooRealVar("coeff_2", "coeff_2",_c2);
        _r_c3 = new RooRealVar("coeff_3", "coeff_3",_c3);
        _r_c4 = new RooRealVar("coeff_4", "coeff_4",_c4);
        _r_c5 = new RooRealVar("coeff_5", "coeff_5",_c5);
        _r_c6 = new RooRealVar("coeff_6", "coeff_6",_c6);
        _r_c7 = new RooRealVar("coeff_7", "coeff_7",_c7);
        _r_c8 = new RooRealVar("coeff_8", "coeff_8",_c8);
        _r_c9 = new RooRealVar("coeff_9", "coeff_9",_c9);
        
        _r_p0_os = new RooRealVar("p0_os", "p0_os",_p0_os);
        _r_p1_os = new RooRealVar("p1_os", "p1_os",_p1_os);
        _r_delta_p0_os = new RooRealVar("delta_p0_os", "delta_p0_os",_delta_p0_os);
        _r_delta_p1_os = new RooRealVar("delta_p1_os", "delta_p1_os",_delta_p1_os);
        _r_avg_eta_os = new RooRealVar("avg_eta_os", "avg_eta_os",_avg_eta_os);
        _r_tageff_os = new RooRealVar("tageff_os", "tageff_os",_tageff_os);
        _r_tageff_asym_os = new RooRealVar("tageff_asym_os", "tageff_asym_os",_tageff_asym_os);
        _r_p0_ss = new RooRealVar("p0_ss", "p0_ss",_p0_ss);
        _r_p1_ss = new RooRealVar("p1_ss", "p1_ss",_p1_ss);
        _r_delta_p0_ss = new RooRealVar("delta_p0_ss", "delta_p0_ss",_delta_p0_ss);
        _r_delta_p1_ss = new RooRealVar("delta_p1_ss", "delta_p1_ss",_delta_p1_ss);
        _r_avg_eta_ss = new RooRealVar("avg_eta_ss", "avg_eta_ss",_avg_eta_ss);
        _r_tageff_ss = new RooRealVar("tageff_ss", "tageff_ss",_tageff_ss);
        _r_tageff_asym_ss= new RooRealVar("tageff_asym_ss", "tageff_asym_ss",_tageff_asym_ss);
        
        _r_production_asym = new RooRealVar("production_asym", "production_asym",_production_asym);
        _r_detection_asym = new RooRealVar("detection_asym", "detection_asym",_detection_asym);
        
        /// Acceptance
        vector<RooRealVar*> v_coeff;
        v_coeff.push_back(_r_c0);
        v_coeff.push_back(_r_c1);
        v_coeff.push_back(_r_c2);
        v_coeff.push_back(_r_c3);
        v_coeff.push_back(_r_c4);
        v_coeff.push_back(_r_c5);
        v_coeff.push_back(_r_c6);
        v_coeff.push_back(_r_c7);
        v_coeff.push_back(_r_c8);
        v_coeff.push_back(_r_c9);
 
        //SPLINE KNOTS
        NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        vector<double> myBinning = knot_positions.getVector();
        
        RooArgList tacc_list;
        for(int i= 0; i<= myBinning.size(); i++){
            tacc_list.add(*v_coeff[i]);
        }
        
        _coeff_last = new RooFormulaVar(("coeff_"+anythingToString((int)myBinning.size()+1)).c_str(),("coeff_"+anythingToString((int)myBinning.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString((int)myBinning.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
        
        tacc_list.add(*_coeff_last);	
        
        // Cubic spline function
        _spline = new RooCubicSplineFun("splinePdf", "splinePdf", *_r_t, myBinning, tacc_list);        
        //_efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt, *_r_scale_mean_dt, *_r_scale_sigma_dt );
        
        _r_dt_scaled = (RooRealVar*) new RooFormulaVar( "r_dt_scaled","r_dt_scaled", "@0+@1*@2",RooArgList(*_r_offset_sigma_dt,*_r_scale_sigma_dt,*_r_dt));
        _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, RooRealConstant::value(0.), *_r_dt_scaled, *_r_scale_mean_dt, RooRealConstant::value(1.) );
        
        // CP coefficients
        _r_C = new RooRealVar("C", "C",0,-4,4);
        _r_C_bar = new RooRealVar("C", "C",0,-4,4);
        //_r_C_bar= (RooRealVar*) new RooFormulaVar("Cbar","-1. * @0",RooArgList(*_r_C));
        _r_D= new RooRealVar("D", "D",0,-4,4);
        _r_D_bar= new RooRealVar("Dbar", "Dbar",0,-4,4);
        _r_S= new RooRealVar("S", "S",0,-4,4);
        _r_S_bar= new RooRealVar("Sbar", "Sbar",0,-4,4);
        
        /// Decay rate coefficients
        _cosh_coeff = new DecRateCoeff_Bd("cosh_coeff",
                                          "cosh_coeff",
                                          DecRateCoeff_Bd::kCosh ,
                                          *_r_f,
                                          RooRealConstant::value(1.),
                                          RooRealConstant::value(1.),
                                          *_r_q_OS,
                                          *_r_eta_OS,
                                          *_r_p0_os,
                                          *_r_p1_os,
                                          *_r_delta_p0_os,
                                          *_r_delta_p1_os,
                                          *_r_avg_eta_os,
                                          *_r_tageff_os,
                                          *_r_tageff_asym_os,
                                          *_r_q_SS,
                                          *_r_eta_SS,
                                          *_r_p0_ss,
                                          *_r_p1_ss,
                                          *_r_delta_p0_ss,
                                          *_r_delta_p1_ss,
                                          *_r_avg_eta_ss,
                                          *_r_tageff_ss,
                                          *_r_tageff_asym_ss,
                                          *_r_production_asym,
                                          *_r_detection_asym );
        
        _cos_coeff = new DecRateCoeff_Bd("cos_coeff",
                                         "cos_coeff",
                                         DecRateCoeff_Bd::kCos ,
                                         *_r_f,
                                         *_r_C,
                                         *_r_C_bar,
                                         *_r_q_OS,
                                         *_r_eta_OS,
                                         *_r_p0_os,
                                         *_r_p1_os,
                                         *_r_delta_p0_os,
                                         *_r_delta_p1_os,
                                         *_r_avg_eta_os,
                                         *_r_tageff_os,
                                         *_r_tageff_asym_os,
                                         *_r_q_SS,
                                         *_r_eta_SS,
                                         *_r_p0_ss,
                                         *_r_p1_ss,
                                         *_r_delta_p0_ss,
                                         *_r_delta_p1_ss,
                                         *_r_avg_eta_ss,
                                         *_r_tageff_ss,
                                         *_r_tageff_asym_ss,
                                         *_r_production_asym,
                                         *_r_detection_asym );
        
        _sinh_coeff = new DecRateCoeff_Bd("sinh_coeff",
                                          "sinh_coeff",
                                          DecRateCoeff_Bd::kSinh ,
                                          *_r_f,
                                          *_r_D,
                                          *_r_D_bar,
                                          *_r_q_OS,
                                          *_r_eta_OS,
                                          *_r_p0_os,
                                          *_r_p1_os,
                                          *_r_delta_p0_os,
                                          *_r_delta_p1_os,
                                          *_r_avg_eta_os,
                                          *_r_tageff_os,
                                          *_r_tageff_asym_os,
                                          *_r_q_SS,
                                          *_r_eta_SS,
                                          *_r_p0_ss,
                                          *_r_p1_ss,
                                          *_r_delta_p0_ss,
                                          *_r_delta_p1_ss,
                                          *_r_avg_eta_ss,
                                          *_r_tageff_ss,
                                          *_r_tageff_asym_ss,
                                          *_r_production_asym,
                                          *_r_detection_asym );
        
        _sin_coeff = new DecRateCoeff_Bd("sin_coeff",
                                         "sin_coeff",
                                         DecRateCoeff_Bd::kSin,
                                         *_r_f,
                                         *_r_S,
                                         *_r_S_bar,
                                         *_r_q_OS,
                                         *_r_eta_OS,
                                         *_r_p0_os,
                                         *_r_p1_os,
                                         *_r_delta_p0_os,
                                         *_r_delta_p1_os,
                                         *_r_avg_eta_os,
                                         *_r_tageff_os,
                                         *_r_tageff_asym_os,
                                         *_r_q_SS,
                                         *_r_eta_SS,
                                         *_r_p0_ss,
                                         *_r_p1_ss,
                                         *_r_delta_p0_ss,
                                         *_r_delta_p1_ss,
                                         *_r_avg_eta_ss,
                                         *_r_tageff_ss,
                                         *_r_tageff_asym_ss,
                                         *_r_production_asym,
                                         *_r_detection_asym );
        
        /// Time pdf integrators
        _cosh_term = new TimePdfIntegrator(coshBasis,_tau,_dGamma,_dm,_efficiency);
        _cos_term = new TimePdfIntegrator(cosBasis,_tau,_dGamma,_dm,_efficiency);
        _sinh_term = new TimePdfIntegrator(sinhBasis,_tau,_dGamma,_dm,_efficiency);
        _sin_term = new TimePdfIntegrator(sinBasis,_tau,_dGamma,_dm,_efficiency);
        
        // Marginal pdfs
        TFile* f_pdfs = new TFile("Mistag_pdfs.root","OPEN");
        
        _h_dt = new TH1D( *((TH1D*) f_pdfs->Get("h_dt_norm")));
        _r_h_dt = new RooDataHist("r_h_dt","r_h_dt",*_r_dt,_h_dt);
        _pdf_sigma_t = (RooAbsPdf*) (new RooHistPdf("pdf_sigma_t","pdf_sigma_t",*_r_dt,*_r_h_dt));
        
        _h_eta_OS = new TH1D( *((TH1D*) f_pdfs->Get("h_w_OS_norm")));
        _r_h_eta_OS = new RooDataHist("r_eta_OS","r_eta_OS",*_r_eta_OS,_h_eta_OS);
        _pdf_eta_OS = (RooAbsPdf*) (new RooHistPdf("pdf_eta_OS","pdf_eta_OS",*_r_eta_OS,*_r_h_eta_OS));
        
        _h_eta_SS = new TH1D( *((TH1D*) f_pdfs->Get("h_w_SS_norm")));
        _r_h_eta_SS = new RooDataHist("r_eta_SS","r_eta_SS",*_r_eta_SS,_h_eta_SS);
        _pdf_eta_SS = (RooAbsPdf*) (new RooHistPdf("pdf_eta_SS","pdf_eta_SS",*_r_eta_SS,*_r_h_eta_SS));
        
        _pdf_eta_OS_uniform = new RooUniform("pdf_eta_OS_uniform","pdf_eta_OS_uniform",*_r_eta_OS);
        _pdf_eta_SS_uniform = new RooUniform("pdf_eta_SS_uniform","pdf_eta_SS_uniform",*_r_eta_SS);
        
        f_pdfs->Close();

    }
         
    double get_cosh_term_Val(IDalitzEvent& evt){
        return _cosh_coeff->evaluate() * _cosh_term->getVal(evt).real();
    }

    double get_cosh_term_Integral(IDalitzEvent& evt){
        return _cosh_coeff->analyticalIntegral(2) * _cosh_term->getVal(evt).imag();
    }

    double get_sinh_term_Val(IDalitzEvent& evt){
        return _sinh_coeff->evaluate() * _sinh_term->getVal(evt).real();
    }

    double get_sinh_term_Integral(IDalitzEvent& evt){
        return _sinh_coeff->analyticalIntegral(2) * _sinh_term->getVal(evt).imag();
    }
    
    double get_cos_term_Val(IDalitzEvent& evt){
        return _cos_coeff->evaluate() * _cos_term->getVal(evt).real();
    }
    
    double get_cos_term_Integral(IDalitzEvent& evt){
        return _cos_coeff->analyticalIntegral(2) * _cos_term->getVal(evt).imag();
    }
    
    double get_sin_term_Val(IDalitzEvent& evt){
        return _sin_coeff->evaluate() * _sin_term->getVal(evt).real();
    }
    
    double get_sin_term_Integral(IDalitzEvent& evt){
        return _sin_coeff->analyticalIntegral(2) * _sin_term->getVal(evt).imag();
    }
    
    double get_marginalPdfs_Val(IDalitzEvent& evt){
        const int q_OS = (int)evt.getValueFromVector(3);
        const int q_SS = (int)evt.getValueFromVector(5);
        
        return  _pdf_sigma_t->getVal()
                * ( abs(q_OS)/2. * _pdf_eta_OS->getVal() + ( 1. - abs(q_OS)) * _pdf_eta_OS_uniform->getVal() )
                * ( abs(q_SS)/2. * _pdf_eta_SS->getVal() + ( 1. - abs(q_SS)) * _pdf_eta_SS_uniform->getVal() );
    }
    
    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _cosh_coeff->calibrate(evt.getValueFromVector(4), _avg_eta_os, _p0_os, _p1_os, _delta_p0_os, _delta_p1_os);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _cosh_coeff->calibrate(evt.getValueFromVector(6), _avg_eta_ss, _p0_ss, _p1_ss, _delta_p0_ss, _delta_p1_ss);
    }
    
   void listFitParDependencies(){
       std::cout << "TimePDfMaster depends on these fitParameters:" << std::endl;
       _cosh_term->listFitParDependencies(std::cout);
       //std::cout << "TimePDfMaster depends on these fitParameters (sinh):" << std::endl;
       //_sinh_term->listFitParDependencies(std::cout);
       //std::cout << "" << std::endl;
   }
   
   void setAllFitParameters(){       
       _r_tau->setVal(_tau);
       _r_dGamma->setVal(_dGamma);
       _r_dm->setVal(_dm);
       
       _r_scale_mean_dt->setVal(_scale_mean_dt);
       _r_offset_sigma_dt->setVal(_offset_sigma_dt);
       _r_scale_sigma_dt->setVal(_scale_sigma_dt);
       _r_c0->setVal(_c0);
       _r_c1->setVal(_c1);
       _r_c2->setVal(_c2);
       _r_c3->setVal(_c3);
       _r_c4->setVal(_c4);
       _r_c5->setVal(_c5);
       _r_c6->setVal(_c6);
       _r_c7->setVal(_c7);
       _r_c8->setVal(_c8);
       _r_c9->setVal(_c9);
       
       _r_p0_os->setVal(_p0_os);
       _r_p1_os->setVal(_p1_os);
       _r_delta_p0_os->setVal(_delta_p0_os);
       _r_delta_p1_os->setVal(_delta_p1_os);
       _r_avg_eta_os->setVal(_avg_eta_os);
       _r_tageff_os->setVal(_tageff_os);
       _r_tageff_asym_os->setVal(_tageff_asym_os);
       _r_p0_ss->setVal(_p0_ss);
       _r_p1_ss->setVal(_p1_ss);
       _r_delta_p0_ss->setVal(_delta_p0_ss);
       _r_delta_p1_ss->setVal(_delta_p1_ss);
       _r_avg_eta_ss->setVal(_avg_eta_ss);
       _r_tageff_ss->setVal(_tageff_ss);
       _r_tageff_asym_ss->setVal(_tageff_asym_ss);
       
       _r_production_asym->setVal(_production_asym);
       _r_detection_asym->setVal(_detection_asym);
   }
   
   void setAllObservables(IDalitzEvent& evt){
       const double t = (double) evt.getValueFromVector(0);
       const double dt = (double) evt.getValueFromVector(1);
       const int f = (int)evt.getValueFromVector(2);
       const int q_OS = (int)evt.getValueFromVector(3);
       const double eta_OS = (double) evt.getValueFromVector(4);
       const int q_SS = (int)evt.getValueFromVector(5);
       const double eta_SS = (double) evt.getValueFromVector(6);
       
       _r_t->setVal(t);
       _r_dt->setVal(dt);
       _r_f->setIndex(f);
       _r_q_OS->setIndex(q_OS);
       _r_eta_OS->setVal(eta_OS);
       _r_q_SS->setIndex(q_SS);
       _r_eta_SS->setVal(eta_SS);
    }
    
    void setCP_coeff(double C,double C_bar,double D,double D_bar,double S,double S_bar ){
        _r_C->setVal(C);
        _r_C_bar->setVal(C_bar);
        _r_D->setVal(D);
        _r_D_bar->setVal(D_bar);
        _r_S->setVal(S);
        _r_S_bar->setVal(S_bar);
    }
    
    void setAllObservablesAndFitParameters(IDalitzEvent& evt){
        setAllObservables(evt);
        setAllFitParameters();
    }

   virtual ~TimePdfMaster(){
       // Plot acceptance
       TH1F *h_spline = new TH1F("", "", 100, _min_TAU, _max_TAU);
       for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
           _r_t->setVal(h_spline->GetXaxis()->GetBinCenter(i));
           h_spline->SetBinContent(i,_spline->getVal());
       }
       
       TCanvas* c = new TCanvas();
       h_spline->SetLineColor(kRed);
       h_spline->Draw("histc");
       c->Print("spline.eps");
       //c->Print("spline.pdf");
   }

};

#endif
//
