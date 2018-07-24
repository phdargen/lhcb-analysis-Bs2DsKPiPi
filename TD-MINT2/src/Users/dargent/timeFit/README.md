# Generate toys

./timeFit randomSeed "gen" < signal_toy.txt  

-> Outputfile: toys_randomSeed.root


# Fit toy sample

./timeFit randomSeed "fit" < signal_toy.txt

-> Fit parameters saved in pull_randomSeed.root


# Submit batch jobs (@ sigma0 -p26)

qsub qGen.sh

qsub qFit.sh


# Analyze pulls

root -l

.L pull.C

pull p
