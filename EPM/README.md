
This directory contains the options files in 'Options' used to perform the flavour tagging calibrations. Tutorial Slides can be found in the subdirectory 'Slides'.

After compiling the code, you can run from the EPM main folder using

./bin/SimpleEvaluator "path-to-optionsfile"


Options that you can set in the files are:
---------------------------------------------------
---------------------------------------------------

- General:
---------------------------------------------------
datafile = "path-to-rootfile"

TupleName = "TreeName"

CalibrationMode = "Bs"
  
DoCalibrations = 1  # set 0 for importing a calibration from CalibFile

CalibrationLink = "MISTAG"

CalibrationDegree = 1 # the usual linear calibration model, set 2 for quadratic model, 3 for cubic ...  

CalibrationModel = "POLY"

UseNewtonRaphson = 0 # use minuit (more stable), NewtonRaphson is good for big samples


Selection = "CutString" # e.g. year == 11 || year == 12 for Run1 only

Nmax = -1 


- Plotting (self-explaining):
---------------------------------------------------
PlotLabel = "LHCb"

PlotTitle = 0

PlotExtension = ".eps"

PlotStatBox = 0


- Evaluation Options:
---------------------------------------------------
 BranchID      = "Bs_ID" # Bs_TRUEID for MC

UseWeight      = 1 # using sWeights / MC weights

BranchWeight   = "N_Bs_sw" 


UseTau  = 1 # set 0 for MC, where TRUEID = ID at production 

TypeTau = "Double_t"

TauUnits = "ps"

BranchTau = "Bs_DTF_TAU"

UseTauErr = 1

TypeTauErr = "Double_t"

BranchTauErr = "Bs_DTF_TAUERR"


DrawOscillationPlots = 1 # draws  A = mixed-unmixed (normalized to 1) 

OscillationPlotsMaximum = 1 # default, draws one period  


- Resolution Scaling:
---------------------------------------------------
ResolutionGaussian1_A = 0.0103 # in ps (same as on our tuples)

ResolutionGaussian1_B = 1.28



This corresponds to scaling of the form: sigma_t_eff = A + B * sigma_t

In this example A&B are equal to Run1 Bs->DsK calibration

Could add more terms, e.g. C * sigma_t^2 
 

- Taggers to Evaluate:
---------------------------------------------------
 OS_Muon_Use = 1

OS_Muon_TypeDec          = "Short_t"

OS_Muon_BranchDec        = "Bs_OS_Muon_DEC"

OS_Muon_TypeProb        = "Float_t"

OS_Muon_BranchProb      = "Bs_OS_Muon_PROB"

...

This sets addresses and data types of taggers that should be evaluated 


- Offline Combinations:
---------------------------------------------------
PerformOfflineCombination_OS = 1 # set for combining following taggers to OS_Combination

OS_Muon_InOSComb = 1

OS_Electron_InOSComb = 1

OS_nnetKaon_InOSComb = 1

VtxCharge_InOSComb = 1



PerformOfflineCombination_OSplusSS = 1 # set for combining OS_Combination + one SS tagger

OS_Combination_InComb = 1

SS_nnetKaon_InComb = 1


-Save Calibration Output:
---------------------------------------------------
WriteCalibratedMistagBranches = 1

OS_Combination_Write = 1

Combination_Write = 1  

CalibratedOutputFile = <"outputfileName.root">


saves the combined tagging decision and probabilities in <outputfile>


-Import Calibration from file:
-------------------------------------------------
import "path-to-CalibFile.py"

In this case set DoCalibrations = 0 above.

A saved calibration is imported and used to calibrate the branches, no fits are performed.

-------------------------------------------------
-------------------------------------------------
-------------------------------------------------
-------------------------------------------------


After completing calibration you can merge the input tree and the tree with with the combined OS and combined OS+SS branches using mergeTrees.cpp:

open root shell

>.L mergeTrees.cpp

>mergeTrees("path-to-inputfile","path-to-taggingOutput","CutString")


Enter inputfile and taggingOutput without ".root" extension.

Make sure that CutString is identical to the cuts used in optionsfile (Selection = "CutString").


