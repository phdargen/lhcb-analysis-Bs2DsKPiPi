# compile:
make

# options (set in TD_Bs.txt):
--------------------------------

# General options:

"Event Pattern": Gives PDG id of Bs and final state particles 

# Setup of MC integrator:

"IntegratorEventFile": File with MC integration events, if it doesn't exists it will be generated

"IntegratorEvents"	Number of generated integration events

"IntegPrecision": Required precision of MC integration


# Options of toy generator:

"InputFileName": Input file used for the fit, if none is given toys will be generated

"Nevents": Number of generated events

"saveEvents": Save generated events ?

"OutputRootFile": Name of output file with generated evenst

# Fit options:

"doTimeFit": Do phasespace integrated fit in addition to fit with full PDF ?

"do2DScan":	Do likelihood scan of gamma vs delta ?

"doPlots":  Produce plots?

"OutputDir":  Directory for plots

# Fit parameters:

Fit parameters are initiallized in the format:

name		Fix?	init	step	min	max

for example:

"gamma"		0 70 1 0.0 0.0

Fix = 0 == Fixed
Fix = 2 == Floating

# run it:
./ampFit 0 < TD_Bs.txt 

or to save std output in txt file:

./ampFit 0 < TD_Bs.txt > out.txt

