#####################################
# FLAVOURTAGGINGTOOLS CONFIGURATION #
#####################################

########
# DATA #
########

#  This is the file/directory that you want to run (LAST ONE IS TAKEN):
#  if a directory is given all root files in it will be read:
  
datafile = "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root"
TupleName = "DecayTree"

###########
# GENERIC #
###########

CalibrationMode = "Bs"
DoCalibrations = 1
CalibrationLink = "MISTAG"
CalibrationDegree = 1
CalibrationModel = "POLY"
UseNewtonRaphson = 0

Selection = "run == 1"
#DBGLEVEL    =  5
Nmax        = -1  # Events to run, -1 means all

##################
# SET LHCB STYLE #
##################
PlotLabel = "LHCb"
PlotTitle = 0
PlotExtension = ".pdf"
PlotStatBox = 0

###################
# SIMPLEEVALUATOR #
###################

BranchID             = "Bs_ID"
UseWeight            = 1
BranchWeight         = "N_Bs_sw"
UseTau  = 1
TypeTau = "Double_t"
TauUnits = "ps"
BranchTau = "Bs_BsDTF_TAU"
UseTauErr = 1
TypeTauErr = "Double_t"
BranchTauErr = "Bs_BsDTF_TAUERR"

#Reso-scaling from Bs->DsK data:
ResolutionGaussian1_A = 0.0103
ResolutionGaussian1_B = 1.28

DrawOscillationPlots = 1
#OscillationPlotsMaximum = 1.1

###################
# SPECIFY TAGGERS #
###################

# TaggerName_NumBins = 20

### TAGGERS USED IN Bs->DsK
OS_Muon_Use = 1
OS_Muon_TypeDec          = "Int_t"
OS_Muon_BranchDec        = "OS_Muon_DEC"
OS_Muon_TypeProb        = "Double_t"
OS_Muon_BranchProb      = "OS_Muon_PROB"

OS_nnetKaon_Use = 1
OS_nnetKaon_TypeDec        = "Int_t"
OS_nnetKaon_BranchDec      = "OS_nnetKaon_DEC"
OS_nnetKaon_TypeProb      = "Double_t"
OS_nnetKaon_BranchProb    = "OS_nnetKaon_PROB"

OS_Electron_Use = 1
OS_Electron_TypeDec      = "Int_t"
OS_Electron_BranchDec    = "OS_Electron_DEC"
OS_Electron_TypeProb    = "Double_t"
OS_Electron_BranchProb  = "OS_Electron_PROB"

OS_Kaon_Use = 1
OS_Kaon_TypeDec        = "Int_t"
OS_Kaon_BranchDec      = "OS_Kaon_DEC"
OS_Kaon_TypeProb      = "Double_t"
OS_Kaon_BranchProb    = "OS_Kaon_PROB"

VtxCharge_Use = 1
VtxCharge_TypeDec     = "Int_t"
VtxCharge_BranchDec   = "OS_VtxCharge_DEC"
VtxCharge_TypeProb   = "Double_t"
VtxCharge_BranchProb = "OS_VtxCharge_PROB"

#SS_nnetKaon_Use = 1
#SS_nnetKaon_TypeDec      = "Int_t"
#SS_nnetKaon_BranchDec    = "SS_Kaon_DEC"
#SS_nnetKaon_TypeProb    = "Double_t"
#SS_nnetKaon_BranchProb  = "SS_Kaon_PROB"

#SS_Kaon_Use = 1
#SS_Kaon_TypeDec      = "Short_t"
#SS_Kaon_BranchDec    = "Bs_SS_Kaon_DEC"
#SS_Kaon_TypeProb    = "Float_t"
#SS_Kaon_BranchProb  = "Bs_SS_Kaon_PROB"

## BUGGY ?
#OS_Charm_Use = 1
#OS_Charm_TypeDec = "Short_t"
#OS_Charm_BranchDec = "Bs_OS_Charm_DEC"
#OS_Charm_TypeProb = "Float_t"
#OS_Charm_BranchProb = "Bs_OS_Charm_PROB"

#OS_Combination_Use  = 1
#OS_Combination_TypeDec	= "Int_t"
#OS_Combination_BranchDec  = "Bs_TAGDECISION_OS"
#OS_Combination_TypeProb	= "Double_t"
#OS_Combination_BranchProb = "Bs_TAGOMEGA_OS"


### OS AND OS+SS COMBINATION

#PerformOfflineCombination_OS = 1
#OS_Muon_InOSComb = 1
#OS_Electron_InOSComb = 1
#OS_nnetKaon_InOSComb = 1
#VtxCharge_InOSComb = 1

#PerformOfflineCombination_OSplusSS = 1
#OS_Combination_InComb = 1
#SS_nnetKaon_InComb = 1

############################
## SAVE CALIBRATION OUTPUT #
############################

#WriteCalibratedMistagBranches = 1
#OS_Combination_Write = 1
#OS_Muon_Write = 1
#OS_Electron_Write = 1
#OS_nnetKaon_Write = 1
#SS_nnetKaon_Write = 1
#Combination_Write = 1

#CalibratedOutputFile = "Bs2Dspipipi_TaggingCalibration_Output.root"


############################
# CALIBRATION INPUT VALUES #
############################

#import EspressoCalibrations.py
#SaveCalibrationsToXML = 1
