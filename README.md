All code used for the analysis should go in this repository along with the latest version of the Ana-Note and the latest paper draft.
Feel free to add usefull information to this README :)

Some comments on the code:

BDTSelection.cc :
Macro that is used for the complete (pre-) selection of Bs->DsKpipi signal candidates, which can then be further selected using TMVA. There are some bool variables at the top of the code where one can specify the Ds decay mode (Ds->3pi, Ds->KKpi , Ds->Kpipi), the data sample (7 TeV or 8 TeV) and whether we use MC or Data.

----------------------------------------------------

Bs2DsKpipi.cc :
Various functions are defined in this macro:

void preselect() - obsolete old function used to preselect signal candidates using cut strings, to prone to misstakes 

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

----------------------------------------------------

Bs2DsPiPiPi.cc :
Same as BDTSelection.cc , but for Bs->DsPiPiPi normalization channel and without nice bool variables and only for Ds->KKpi. Will be updated so that it can be used in same way as BDTSelection.cc .

----------------------------------------------------

MCvsData.cc :
Macro that compares distributions of observables in Data (s-weighted to get signal component) and MC (truth matched). A Kolmogorov test is included.

---------------------------------------------------

TMVAClassificationDsKpipi.C:
Does the BDT training for variables specified in the code. Is a bit messy with a lot of stuff commented out, but also its not clear whether variables are final

--------------------------------------------------

TMVAClassificationApplication_DsKpipi.C:
Applies the BDT weights to every event. Also a bit messy, see above :D

---------------------------------------------------
---------------------------------------------------

****Data samples:******
Are currently only available on the heidelberg /auto/data storage. We could upload them also to eos, but we will have to make a preselection to reduce the very big samples (~ 1TB for one polarity in 2012).
 
