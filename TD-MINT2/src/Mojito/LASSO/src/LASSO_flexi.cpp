// author: Philippe d'Argent (p.dargent@cern.ch)
#include "Mint/Neg2LLSum.h"
#include "Mint/IMinimisable.h"
#include "Mint/LASSO_flexi.h"
#include "Mint/NamedParameter.h"


using namespace std;
using namespace MINT;


double LASSO_flexi::getVal(){
    return _lambda * _pdf->sumOfSqrtFitFractions();
}

int LASSO_flexi::numberOfFitFractionsLargerThanThreshold(double threshold){
    return _pdf->numberOfFitFractionsLargerThanThreshold(threshold);
}



//
