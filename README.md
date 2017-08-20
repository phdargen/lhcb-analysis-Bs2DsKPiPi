All code used for the analysis should go in this repository along with the latest version of the Ana-Note and the latest paper draft.

Feel free to add usefull information to this README :)

For instructions on how to compile and run the code see README in the corresponding directories.


Workflow:
---------------------------------------------------
---------------------------------------------------

1) Minimize stripped data
---------------------------------------------------

in Selection/ :

source mini_data.sh

source mini_mc.sh

2) For MC: Add corrected PID vars
---------------------------------------------------

in PIDCalib/ :  

source runPIDCorr_signal.sh

source runPIDCorr_norm.sh


3) Preselection
---------------------------------------------------

in Selection/ :

source select_data.sh

source select_mc.sh

4) Fit normalization channel 
---------------------------------------------------

in TD-MINT2/src/Users/dargent/MassFits :

./massFit < fitPreselectedNorm.txt 

5) Reweight MC
---------------------------------------------------

in TD-MINT2/src/Users/dargent/ DataVsMC :

./dataVsMC < dataVsMC.txt 

6) Train and apply BDT
---------------------------------------------------

in Selection: 

Run in root:

TMVAClassification("BDTG")

TMVAClassificationApplication("Signal/Norm","Data/MC")

7) Optimize BDT cut
---------------------------------------------------

in TD-MINT2/src/Users/dargent/MassFits  :

./massFit < fitSignalForBDT.txt 

8) Fit final sample
---------------------------------------------------

in TD-MINT2/src/Users/dargent/MassFits  :

./massFit < massFit.txt 


Data samples:
---------------------------------------------------
---------------------------------------------------

Stripping output after loose preselection is applied and unnecessary branches are removed:
---------------------------------------------------
/auto/data/dargent/BsDsKpipi/Mini/Data(MC)/signal(norm)_Ds2KKpi(Ds2pipipi)_11(12/15/16).root 


After preselection with BDT variables added:
---------------------------------------------------
/auto/data/dargent/BsDsKpipi/Preselected/Data(MC)/signal(norm)_Ds2KKpi(Ds2pipipi)_11(12/15/16).root 


With BDT response:
---------------------------------------------------
/auto/data/dargent/BsDsKpipi/BDT/Data(MC)/signal(norm).root 


Final sample:
---------------------------------------------------
/auto/data/dargent/BsDsKpipi/Final/Data(MC)/signal(norm).root 


Final sample in MINT format:
---------------------------------------------------
/auto/data/dargent/BsDsKpipi/Mint/Data(MC)/signal(norm).root 



