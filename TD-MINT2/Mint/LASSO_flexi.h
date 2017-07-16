#ifndef LASSO_FLEXI_HH
#define LASSO_FLEXI_HH
// author: Philippe d'Argent (p.dargent@cern.ch)

#include "TMath.h"
#include "Mint/Minimisable.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/NamedParameter.h"
#include <vector>

namespace MINT{
    
    class LASSO_flexi: public Minimisable{
    protected:
        DalitzPdfBaseFlexiFastInteg* _pdf;
        
        double _lambda;
        
    public:
        LASSO_flexi(DalitzPdfBaseFlexiFastInteg* pdf, double lambda = 1.)
        : _pdf(pdf), _lambda(lambda) {
            _pdf->redoIntegrator();
        };
        
        virtual void beginFit(){
            _pdf->redoIntegrator();
        };
        virtual void parametersChanged(){
            _pdf->parametersChanged();
        };
        virtual void endFit(){};
        
        virtual double getVal();
        
        virtual double getNewVal(){ 
            parametersChanged();
            return getVal();
        }
        
        int numberOfFitFractionsLargerThanThreshold(double threshold);
        double absSumOfInterferenceFractions() {
            return _pdf->absSumOfInterferenceFractions();
        }
        double sumOfFitFractions(){
            return _pdf->sumOfFitFractions();
        }
        virtual ~LASSO_flexi(){}
        
    };
    
}// namespace MINT
#endif
//
