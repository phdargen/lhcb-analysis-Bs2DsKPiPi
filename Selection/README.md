# compile code:
source make.sh

# Run over stripped data: 
removes not needed branches and applies loose preselection (chose one of the options)

./MiniMaker "Signal\Norm" "Ds2KKpi\Ds2pipipi\Ds2Kpipi" "Data\MC" 11\12\15\16

or submit all jobs to batch

source mini_mc.sh
source mini_data.sh

# Run over minimized sample: 
applies final preselection and adds variables for BDT

./SelectionMaker "Signal\Norm" "Ds2KKpi\Ds2pipipi\Ds2Kpipi" "Data\MC" 11\12\15\16

or submit all jobs to batch

source select_mc.sh
source select_data.sh

# train BDT   

(requires TMVA 4.2.0) 
source setup.sh < /work/dargent/TMVA-v4.2.0/

root -l
.L TMVAClassification.cpp
TMVAClassification("BDTG","Run1")



