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
    MINT::FitParRef _tau,_dGamma,_dm;

 public:
  TimePdfIntegrator( int basisType, 
                    const MINT::FitParameter& tau, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
                    , RooGaussEfficiencyModel* efficiency
                    , IFitParRegister* daddy=0):
                        FitParDependent(daddy),
                       _basisType(basisType), _tau(tau,this),_dGamma(dGamma,this),_dm(dm,this)
                       , _efficiency(efficiency) 
                       {
                       }
         
  virtual std::complex<double> getVal(IDalitzEvent& evt){
    //return getNewVal(evt);// (for debugging)
    return getValWithCaching(evt);
  }

  virtual std::complex<double> getNewVal(IDalitzEvent& evt){
      //return std::complex<double>(_tau,_dGamma);
       return std::complex<double>(_efficiency->evaluate(_basisType,_tau,_dm,_dGamma),_efficiency->analyticalIntegral(_basisType,_tau,_dm,_dGamma));
  }

  virtual std::complex<double> ComplexVal(IDalitzEvent& evt){
      return getVal(evt);
  }


  virtual ~TimePdfIntegrator(){};

};

#endif
//
