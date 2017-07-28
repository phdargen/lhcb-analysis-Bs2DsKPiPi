----------------------------------------------------

Run code with ./massFit < massFit.txt

Options that can be tuned in the optionsfile:


Channel: specify "Normalization" or "Signal" depending on the channel you want to perform the massfit in

Year: takes integer 2011, 2012, 2015 or 2016. Only relevant in case you do not fit simultaneous

DsKKpi: choose 1 for Ds->KKpi final state or 0 for Ds->pipipi

fitSimultaneous: choose 1 for simultaneous fit or 0 for single fit (Year becomes relevant for 0 case)

sWeighting: choose 1 to create sWeights from massfits and save new Ntuple with weights, 0 otherwise

makePlots: choose 1 to save new plots, 0 otherwise

BDTScan: choose 1 if you want to scan BDT cut in Normalization channel, 0 otherwise 

fitCombined: choose 1 to combine years for Run1 and Run2, 0 otherwise 

 
