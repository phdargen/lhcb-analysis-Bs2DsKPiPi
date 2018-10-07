// author: Philippe d'Argent
#include "Mint/GLass.h"
#include <iostream>
#include <cmath>

using namespace std;
using namespace MINT;

std::complex<double> GLass::BreitWigner(){
    
    double q2              = pABSq();
    double q               = sqrt( q2 );
    double gamma           = GofM();
    double scatteringPhase = phiF() + atan( 2. * a() * q / ( 2. + a() * r() * q2 ) );
    double resonancePhase  = phiR() + atan( mumsMass() * gamma / ( mumsMass() * mumsMass() - mumsRecoMass2() ) );
    double rho             = sqrt( q2 / s );
    complex<double> returnValue     = ( F() * sin( scatteringPhase ) * std::polar(1.,scatteringPhase)  +
                                        R() * sin( resonancePhase ) * std::polar(1., resonancePhase + 2. * scatteringPhase ) ) / rho;
        
    std::vector<int> asi = _theDecay.getVal().asi();
    double min = getEvent()->eventPattern().sijMin(asi);
    double max = getEvent()->eventPattern().sijMax(asi);
    double x = 2.* (mumsRecoMass2() - min)/(max - min) - 1.;  

    returnValue *= exp(alpha1()*x) + exp(alpha2()*x*x) + exp(alpha3()*x*x*x);

    return returnValue;
}
