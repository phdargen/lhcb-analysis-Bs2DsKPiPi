// -- CLASS DESCRIPTION [PDF] --

#include "RooFit.h"
#include <iostream>
#include <math.h>

#include "Mint/RooHILLdini_misID.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"

RooHILLdini_misID::RooHILLdini_misID(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _a, RooAbsReal& _b, RooAbsReal& _csi, RooAbsReal& _m1, RooAbsReal& _s1, RooAbsReal& _m2, RooAbsReal& _s2, RooAbsReal& _m3, RooAbsReal& _s3, RooAbsReal& _m4, RooAbsReal& _s4 ,RooAbsReal& _f1, RooAbsReal& _f2, RooAbsReal& _f3 ) :

  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  a("a", "a", this, _a),
  b("b", "b", this, _b),
  csi("csi", "csi", this, _csi),
  m1("m1", "m1", this, _m1),
  s1("s1", "s1", this, _s1),
  m2("m2", "m2", this, _m2),
  s2("s2", "s2", this, _s2),
  m3("m3", "m3", this, _m3),
  s3("s3", "s3", this, _s3),
  m4("m4", "m4", this, _m4),
  s4("s4", "s4", this, _s4),
  f1("f1", "f1", this, _f1),
  f2("f2", "f2", this, _f2),
  f3("f3", "f3", this, _f3)
{
}


RooHILLdini_misID::RooHILLdini_misID(const RooHILLdini_misID& other, const char* name) :
RooAbsPdf(other, name), m("m", this, other.m), a("a", this, other.a), b("b", this, other.b), csi("csi", this, other.csi), m1("m1", this, other.m1), s1("s1", this, other.s1), m2("m2", this, other.m2), s2("s2", this, other.s2), m3("m3", this, other.m3), s3("s3", this, other.s3), m4("m4", this, other.m4), s4("s4", this, other.s4), f1("f1", this, other.f1), f2("f2", this, other.f2), f3("f3", this, other.f3)
{
}





Double_t RooHILLdini_misID::evaluate() const 
{

  ///////////////////////////////
  //double t = m;
  double a_new  = a;
  double b_new  = b;
  double B_NEW  = (a_new+b_new)/2;
  ///////////////////////////////

  ///////////////////////////////////////////////////////////////////////
  //mult           = ((1-csi)/(b_new-a_new)*m + (b_new*csi - a_new)/(b_new-a_new));
  //double CURVEG1 = fabs((1-csi)*secondG1    + (b_new*csi - a_new)*firstG1);
  ///////////////////////////////////////////////////////////////////////

  /////////////////// 
  double firstG1  = (2*exp(-((a_new-(m-m1))*(a_new-(m-m1))/(2*(s1*s1))))*s1*(b_new-(m-m1))+2*exp(-((b_new-(m-m1))*(b_new-(m-m1))/(2*(s1*s1))))*s1*(-a_new+(m-m1))-sqrt(2*TMath::Pi())*(a_new*b_new+(s1*s1)-(a_new+b_new)*(m-m1)+((m-m1)*(m-m1)))*TMath::Erf((-a_new+(m-m1))/(sqrt(2)*s1))+sqrt(2*TMath::Pi())*(a_new*b_new+(s1*s1)-(a_new+b_new)*(m-m1)+((m-m1)*(m-m1)))*TMath::Erf((-b_new+(m-m1))/(sqrt(2)*s1)))/(2*sqrt(2*TMath::Pi()));
  double  CURVEG1 =    fabs(   (1-csi)/(b_new-a_new)*(m-m1)  + (b_new*csi - a_new)/(b_new-a_new)  )*fabs(firstG1);
  ///////////////////  

  ///////////////////                                                                                                                                                                                                             
  double firstG2  = (2*exp(-((a_new-(m-m2))*(a_new-(m-m2))/(2*(s2*s2))))*s2*(b_new-(m-m2))+2*exp(-((b_new-(m-m2))*(b_new-(m-m2))/(2*(s2*s2))))*s2*(-a_new+(m-m2))-sqrt(2*TMath::Pi())*(a_new*b_new+(s2*s2)-(a_new+b_new)*(m-m2)+((m-m2)*(m-m2)))*TMath::Erf((-a_new+(m-m2))/(sqrt(2)*s2))+sqrt(2*TMath::Pi())*(a_new*b_new+(s2*s2)-(a_new+b_new)*(m-m2)+((m-m2)*(m-m2)))*TMath::Erf((-b_new+(m-m2))/(sqrt(2)*s2)))/(2*sqrt(2*TMath::Pi()));
  double  CURVEG2 =    fabs(   (1-csi)/(b_new-a_new)*(m-m2)  + (b_new*csi - a_new)/(b_new-a_new)  )*fabs(firstG2);
  ///////////////////                                                                                                    


  ///////////////////                                                                                           
  double firstG3  = (2*exp(-((a_new-(m-m3))*(a_new-(m-m3))/(2*(s3*s3))))*s3*(b_new-(m-m3))+2*exp(-((b_new-(m-m3))*(b_new-(m-m3))/(2*(s3*s3))))*s3*(-a_new+(m-m3))-sqrt(2*TMath::Pi())*(a_new*b_new+(s3*s3)-(a_new+b_new)*(m-m3)+((m-m3)*(m-m3)))*TMath::Erf((-a_new+(m-m3))/(sqrt(2)*s3))+sqrt(2*TMath::Pi())*(a_new*b_new+(s3*s3)-(a_new+b_new)*(m-m3)+((m-m3)*(m-m3)))*TMath::Erf((-b_new+(m-m3))/(sqrt(2)*s3)))/(2*sqrt(2*TMath::Pi()));
  double  CURVEG3 =    fabs(   (1-csi)/(b_new-a_new)*(m-m3)  + (b_new*csi - a_new)/(b_new-a_new)  )*fabs(firstG3);
  ///////////////////                                                                                                    

  ///////////////////
  double firstG4  = (2*exp(-((a_new-(m-m4))*(a_new-(m-m4))/(2*(s4*s4))))*s4*(b_new-(m-m4))+2*exp(-((b_new-(m-m4))*(b_new-(m-m4))/(2*(s4*s4))))*s4*(-a_new+(m-m4))-sqrt(2*TMath::Pi())*(a_new*b_new+(s4*s4)-(a_new+b_new)*(m-m4)+((m-m4)*(m-m4)))*TMath::Erf((-a_new+(m-m4))/(sqrt(2)*s4))+sqrt(2*TMath::Pi())*(a_new*b_new+(s4*s4)-(a_new+b_new)*(m-m4)+((m-m4)*(m-m4)))*TMath::Erf((-b_new+(m-m4))/(sqrt(2)*s4)))/(2*sqrt(2*TMath::Pi()));
  double  CURVEG4 =    fabs(   (1-csi)/(b_new-a_new)*(m-m4)  + (b_new*csi - a_new)/(b_new-a_new)  )*fabs(firstG4);
  ///////////////////                                                                                                    

  return fabs(f1*CURVEG1)+ fabs(f2*CURVEG2)+  fabs(f3*CURVEG3) + fabs((1-f1-f2-f3)*CURVEG4);
}
