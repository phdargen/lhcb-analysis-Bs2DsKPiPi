----------------------------------------------------

Bs2DsKpipi.cc :
Various functions are defined in this macro:

void applyBDTcut(string "cutValue") - applies BDT cut of value "cutValue" to ntuple specified in the function 

void iterateBDT(double startValue, double stopValue, double stepSize) - iterates BDT cuts from startValue to stopValue in stepSize increments and does massfit to obtain S/sqrt(s + B) every iteration

double *fitBGShape(double fitValues[int num parameters]) - fits mass shape of Bs->Ds*Kpipi reconstructed as Bs->DsKpipi and fills fitValue array to use shape as input for other functions

double *fitBGShapeNorm(double fitValues[int numParameters]) - fits mass shape of Bs->Ds*pipipi reconstructed as Bs->Dspipipi and fills fitValue array to use shape as input for other functions 

double *fitBGShapeNormKpipi(double fitValues[int numParameters]) - fits mass shape of Bs->DsKpipi reconstructed as Bs->Dspipipi and fills fitValue array to use shape as input for other functions

double *fitBGShapeNormDstKpipi(double fitValues[int numParameters]) - fits mass shape of Bs->Ds*Kpipi reconstructed as Bs->Dspipipi and fills fitValue array to use shape as input for other functions

double *fitBGShapethreePi(double fitValues[int numParameters]) - fits mass shape of Bs->Dspipipi reconstructed as Bs->DsKpipi and fills fitValue array to use shape as input for other functions

double *fitBGShapethreePiDstar(double fitValues[int numParameters]) - fits mass shape of Bs->Ds*pipipi reconstructed as Bs->DsKpipi and fills fitValue array to use shape as input for other functions

void fitBDT() - fits fully selected signal mass distribution with all shapes and saves sWeights if option is set

double fitBDTNorm() - fits fully selected normalization channel mass shape with all shapes

void quickFit()	- implementation of quick gaussian+exp PDF to estimate signal component 

void makePlots() - makes plots

void MCStudies() - uses fake rates from PID calib to determin what percentage of faked events are within our massfit (4800-5800 MeV) region

 
