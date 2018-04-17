#####################################
# FLAVOURTAGGINGTOOLS CONFIGURATION #
#####################################

#  This is a python-like options file
#  # as FIRST char is used to comment out a specific line
#  Separators like # or "" are ignored
#  For names of the form A.B.C.Z, only Z is parsed
#  This allows copy-and-pasting of Gaudi python files, which might look like, e.g.
#  tag.BTaggingTool.CombineTaggersProbability.Eta_Cal_OS = 0.365
#  => Eta_Cal_OS = 0.365
#  Ideally, Gaudi will eventually use calibration parameters like mine, e.g.
#  <python object>.OS_Charm_Eta = X.YYYYY
#  which will directly modify the Tagger object via a Calibration object
#  I will advocate for this.

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


Selection = "year == 11 || year == 12"
#Selection = ""
#DBGLEVEL    =  5
Nmax        = -1  # Events to run, -1 means all


##################
# SET LHCB STYLE #
##################
#Dpesn't work for some reason ... 
PlotLabel = "LHCb"
PlotTitle = 0
PlotExtension = ".eps"
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
BranchTau = "Bs_DTF_TAU"
UseTauErr = 1
TypeTauErr = "Double_t"
BranchTauErr = "Bs_DTF_TAUERR"

#Reso-scaling from Bs->DsK data:
ResolutionGaussian1_A = 0.0103
ResolutionGaussian1_B = 1.28

#ResolutionGaussian1_A = 0
#ResolutionGaussian1_B = 1.25


DrawOscillationPlots = 1
#OscillationPlotsMaximum = 1.1
###################
# SPECIFY TAGGERS #
###################

# TaggerName_NumBins = 20

### TAGGERS USED IN Bs->DsK
OS_Muon_Use = 1
OS_Muon_TypeDec          = "Short_t"
OS_Muon_BranchDec        = "Bs_OS_Muon_DEC"
OS_Muon_TypeProb        = "Float_t"
OS_Muon_BranchProb      = "Bs_OS_Muon_PROB"

OS_nnetKaon_Use = 1
OS_nnetKaon_TypeDec        = "Short_t"
OS_nnetKaon_BranchDec      = "Bs_OS_nnetKaon_DEC"
OS_nnetKaon_TypeProb      = "Float_t"
OS_nnetKaon_BranchProb    = "Bs_OS_nnetKaon_PROB"

OS_Electron_Use = 1
OS_Electron_TypeDec      = "Short_t"
OS_Electron_BranchDec    = "Bs_OS_Electron_DEC"
OS_Electron_TypeProb    = "Float_t"
OS_Electron_BranchProb  = "Bs_OS_Electron_PROB"

#OS_Kaon_Use = 1
#OS_Kaon_TypeDec        = "Short_t"
#OS_Kaon_BranchDec      = "Bs_OS_Kaon_DEC"
#OS_Kaon_TypeProb      = "Float_t"
#OS_Kaon_BranchProb    = "Bs_OS_Kaon_PROB"

VtxCharge_Use = 1
VtxCharge_TypeDec     = "Short_t"
VtxCharge_BranchDec   = "Bs_VtxCharge_DEC"
VtxCharge_TypeProb   = "Float_t"
VtxCharge_BranchProb = "Bs_VtxCharge_PROB"



SS_nnetKaon_Use = 1
SS_nnetKaon_TypeDec      = "Short_t"
SS_nnetKaon_BranchDec    = "Bs_SS_nnetKaon_DEC"
SS_nnetKaon_TypeProb    = "Float_t"
SS_nnetKaon_BranchProb  = "Bs_SS_nnetKaon_PROB"

#SS_Kaon_Use = 1
#SS_Kaon_TypeDec      = "Short_t"
#SS_Kaon_BranchDec    = "Bs_SS_Kaon_DEC"
#SS_Kaon_TypeProb    = "Float_t"
#SS_Kaon_BranchProb  = "Bs_SS_Kaon_PROB"




### OTHER TAGGERS

## BUGGY ?
#OS_Charm_Use = 1
#OS_Charm_TypeDec = "Short_t"
#OS_Charm_BranchDec = "Bs_OS_Charm_DEC"
#OS_Charm_TypeProb = "Float_t"
#OS_Charm_BranchProb = "Bs_OS_Charm_PROB"

#VtxCharge_Use = 1
#VtxCharge_TypeDec     = "Short_t"
#VtxCharge_BranchDec   = "Bs_VtxCharge_DEC"
#VtxCharge_TypeProb   = "Float_t"
#VtxCharge_BranchProb = "Bs_VtxCharge_PROB"
#OS_Charm_Use = 0


#SS_Pion_Use = 1
#SS_Pion_TypeDec        = "Short_t"
#SS_Pion_BranchDec      = "Bs_SS_Pion_DEC"
#SS_Pion_TypeProb      = "Float_t"
#SS_Pion_BranchProb    = "Bs_SS_Pion_PROB"
#SS_PionBDT_Use         = 0
#SS_Proton_Use          = 0
#Combination_Use              = 0

#OS_Combination_Use  = 1
#OS_Combination_TypeDec	= "Int_t"
#OS_Combination_BranchDec  = "Bs_TAGDECISION_OS"
#OS_Combination_TypeProb	= "Double_t"
#OS_Combination_BranchProb = "Bs_TAGOMEGA_OS"


### OS AND OS+SS COMBINATION

PerformOfflineCombination_OS = 1
OS_Muon_InOSComb = 1
OS_Electron_InOSComb = 1
OS_nnetKaon_InOSComb = 1
VtxCharge_InOSComb = 1


PerformOfflineCombination_OSplusSS = 1
OS_Combination_InComb = 1
SS_nnetKaon_InComb = 1


###########################
# SAVE CALIBRATION OUTPUT #
###########################

WriteCalibratedMistagBranches = 1
OS_Combination_Write = 1
OS_Muon_Write = 1
OS_Electron_Write = 1
OS_nnetKaon_Write = 1
SS_nnetKaon_Write = 1
Combination_Write = 1

CalibratedOutputFile = "Bs2Dspipipi_TaggingCalibration_Output.root"


############################
# CALIBRATION INPUT VALUES #
############################

#import EspressoCalibrations.py

