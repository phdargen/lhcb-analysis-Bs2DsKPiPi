#ifndef TIMEPDFINTEGRATOR_HH
#define TIMEPDFINTEGRATOR_HH
// author: Philippe d'Argent

#include <complex>
#include <vector>
#include <iostream>

#include "Mint/IReturnRealForEvent.h"
#include "Mint/IReturnComplexForEvent.h"
#include "Mint/DalitzEvent.h"

#include "Mint/FitParRef.h"
#include "Mint/FitParDependent.h"
#include "Mint/IFitParRegister.h"
#include "Mint/CachedByEvent.h"

#include "Mint/RooGaussEfficiencyModel.h"


class TimePdfIntegrator 
//: virtual public MINT::IReturnRealForEvent<IDalitzEvent>
: virtual public MINT::IReturnComplexForEvent<IDalitzEvent>
, public CachedByEvent<std::complex<double> >
, public MINT::FitParDependent
{
 protected:
    int _basisType;
    RooGaussEfficiencyModel* _efficiency;
    RooGaussEfficiencyModel* _efficiency2;
    RooGaussEfficiencyModel* _efficiency3;

    MINT::FitParRef _Gamma,_dGamma,_dm;

    MINT::FitParRef _offset_mean_dt;
    MINT::FitParRef _scale_mean_dt;
    MINT::FitParRef _scale_mean_2_dt;
    
    MINT::FitParRef _offset_sigma_dt;
    MINT::FitParRef _scale_sigma_dt;
    MINT::FitParRef _scale_sigma_2_dt;

    MINT::FitParRef _offset_sigma2_dt;
    MINT::FitParRef _scale_sigma2_dt;
    MINT::FitParRef _scale_sigma2_2_dt;

    MINT::FitParRef _offset_sigma3_dt;
    MINT::FitParRef _scale_sigma3_dt;
    MINT::FitParRef _scale_sigma3_2_dt;

    MINT::FitParRef _offset_f_dt;
    MINT::FitParRef _scale_f_dt;
    MINT::FitParRef _scale_f_2_dt;

    MINT::FitParRef _offset_f2_dt;
    MINT::FitParRef _scale_f2_dt;
    MINT::FitParRef _scale_f2_2_dt;


    MINT::FitParRef _c0;
    MINT::FitParRef _c1;
    MINT::FitParRef _c2;
    MINT::FitParRef _c3;
    MINT::FitParRef _c4;
    MINT::FitParRef _c5;
    MINT::FitParRef _c6;
    MINT::FitParRef _c7;
    MINT::FitParRef _c8;
    MINT::FitParRef _c9;

 public:
  TimePdfIntegrator( int basisType
                     ,RooGaussEfficiencyModel* efficiency
                     ,RooGaussEfficiencyModel* efficiency2
                     ,RooGaussEfficiencyModel* efficiency3
                     ,const MINT::FitParameter& Gamma, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
		     ,const MINT::FitParameter& offset_mean_dt,const MINT::FitParameter& scale_mean_dt,const MINT::FitParameter& scale_mean_2_dt
		     ,const MINT::FitParameter& offset_sigma_dt,const MINT::FitParameter& scale_sigma_dt,const MINT::FitParameter& scale_sigma_2_dt
		     ,const MINT::FitParameter& offset_sigma2_dt,const MINT::FitParameter& scale_sigma2_dt,const MINT::FitParameter& scale_sigma2_2_dt
		     ,const MINT::FitParameter& offset_sigma3_dt,const MINT::FitParameter& scale_sigma3_dt,const MINT::FitParameter& scale_sigma3_2_dt
		     ,const MINT::FitParameter& offset_f_dt,const MINT::FitParameter& scale_f_dt,const MINT::FitParameter& scale_f_2_dt
		     ,const MINT::FitParameter& offset_f2_dt,const MINT::FitParameter& scale_f2_dt,const MINT::FitParameter& scale_f2_2_dt
                     ,const MINT::FitParameter& c0, const MINT::FitParameter& c1, const MINT::FitParameter& c2
                     ,const MINT::FitParameter& c3, const MINT::FitParameter& c4, const MINT::FitParameter& c5
                     ,const MINT::FitParameter& c6, const MINT::FitParameter& c7, const MINT::FitParameter& c8
                     ,const MINT::FitParameter& c9, IFitParRegister* daddy=0):
    
                        FitParDependent(daddy)
                       ,_basisType(basisType)
                       ,_efficiency(efficiency)
                       ,_efficiency2(efficiency2)
                       ,_efficiency3(efficiency3)
                       ,_Gamma(Gamma,this),_dGamma(dGamma,this),_dm(dm,this)
		       ,_offset_mean_dt(offset_mean_dt,this),_scale_mean_dt(scale_mean_dt,this),_scale_mean_2_dt(scale_mean_2_dt,this)
                       ,_offset_sigma_dt(offset_sigma_dt,this),_scale_sigma_dt(scale_sigma_dt,this),_scale_sigma_2_dt(scale_sigma_2_dt,this)
                       ,_offset_sigma2_dt(offset_sigma2_dt,this),_scale_sigma2_dt(scale_sigma2_dt,this),_scale_sigma2_2_dt(scale_sigma2_2_dt,this)
                       ,_offset_sigma3_dt(offset_sigma3_dt,this),_scale_sigma3_dt(scale_sigma3_dt,this),_scale_sigma3_2_dt(scale_sigma3_2_dt,this)
                       ,_offset_f_dt(offset_f_dt,this),_scale_f_dt(scale_f_dt,this),_scale_f_2_dt(scale_f_2_dt,this)
                       ,_offset_f2_dt(offset_f2_dt,this),_scale_f2_dt(scale_f2_dt,this),_scale_f2_2_dt(scale_f2_2_dt,this)
                       ,_c0(c0,this),_c1(c1,this),_c2(c2,this)
                       ,_c3(c3,this),_c4(c4,this),_c5(c5,this)
                       ,_c6(c6,this),_c7(c7,this),_c8(c8,this)
                       ,_c9(c9,this)
                       {
                       }
         
  virtual std::complex<double> getVal(IDalitzEvent& evt){
    //return getNewVal(evt);// (for debugging)
    return getValWithCaching(evt);
  }

  virtual std::complex<double> getNewVal(IDalitzEvent& evt){
      //return std::complex<double>(_tau,_dGamma);

       std::complex<double> val = std::complex<double>(_efficiency->evaluate(_basisType,1./_Gamma,_dm,_dGamma),_efficiency->analyticalIntegral(_basisType,1./_Gamma,_dm,_dGamma));

       std::complex<double> val2 = std::complex<double>(_efficiency2->evaluate(_basisType,1./_Gamma,_dm,_dGamma),_efficiency2->analyticalIntegral(_basisType,1./_Gamma,_dm,_dGamma));

//        std::complex<double> val3 = std::complex<double>(_efficiency3->evaluate(_basisType,1./_Gamma,_dm,_dGamma),_efficiency3->analyticalIntegral(_basisType,1./_Gamma,_dm,_dGamma));

	double f = _offset_f_dt + _scale_f_dt * evt.getValueFromVector(1) + _scale_f_2_dt * evt.getValueFromVector(1) * evt.getValueFromVector(1);
	
	return f * val + (1.-f) * val2;


       return std::complex<double>(_efficiency->evaluate(_basisType,1./_Gamma,_dm,_dGamma),_efficiency->analyticalIntegral(_basisType,1./_Gamma,_dm,_dGamma));
  }

  virtual std::complex<double> ComplexVal(IDalitzEvent& evt){
      return getVal(evt);
  }


  virtual ~TimePdfIntegrator(){};

};

#endif
//
