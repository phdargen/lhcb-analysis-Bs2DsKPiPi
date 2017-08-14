import os

# List of input ROOT files with MC ntuples. Format: 
#   (inputfile, outputfile, dataset)
files = [
  ("/auto/data/dargent/BsDsKpipi/Mini/MC/norm_Ds2KKpi_11_PID_step1.root", "/auto/data/dargent/BsDsKpipi/Mini/MC/norm_Ds2KKpi_11_PID_step2.root", "MagUp_2012"), 
]

# Name of the input tree
# Could also include ROOT directory, e.g. "Dir/Ntuple"
input_tree = "DecayTree"

# Postfixes of the Pt, Eta and Ntracks variables (ntuple variable name w/o branch name)
# e.g. if the ntuple contains "pion_PT", it should be just "PT"
ptvar  = "PT"
etavar = "ETA"
##pvar   = "p"   # Could use P variable instead of eta
ntrvar = "nTracks" # This should correspond to the number of "Best tracks", not "Long tracks"!

seed = None   # No initial seed
##seed = 1    # Alternatively, could set initial random seed

# Dictionary of tracks with their PID variables, in the form {branch name}:{pidvars}
# For each track branch name, {pidvars} is a dictionary in the form {ntuple variable}:{pid config}, 
#   where 
#     {ntuple variable} is the name of the corresponding ntuple PID variable without branch name, 
#   and 
#     {pid_config} is the string describing the PID configuration. 
# Run PIDCorr.py without arguments to get the full list of PID configs
tracks = {
    'pi_minus'   : {
        "PIDK"       : "pi_CombDLLK", 
            "PIDp"       : "pi_CombDLLp", 
                },
'pi_plus1'   : {
    "PIDK"       : "pi_CombDLLK", 
    "PIDp"       : "pi_CombDLLp", 
        },
'pi_plus2'   : {
    "PIDK"       : "pi_CombDLLK", 
    "PIDp"       : "pi_CombDLLp", 
    },        
'K_plus_fromDs'   : {
    "PIDK"       : "K_CombDLLK", 
        "PIDp"       : "K_CombDLLp", 
            },                    
'K_minus_fromDs'   : {
    "PIDK"       : "K_CombDLLK", 
        "PIDp"       : "K_CombDLLp", 
            },            
'pi_minus_fromDs'   : {
    "PIDK"       : "pi_CombDLLK", 
        "PIDp"       : "pi_CombDLLp", 
            }                    
}

output_tree = input_tree.split("/")[-1]
treename = input_tree

for input_file, output_file, dataset in files : 
  tmpinfile = input_file
  tmpoutfile = "tmp1.root"
  for track, subst in tracks.iteritems() : 
    for var, config in subst.iteritems() : 
#      command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py"
      command = "python ./PIDGen.py"
      command += " -m %s_%s" % (track, ptvar)
      command += " -e %s_%s" % (track, etavar)
##      command += " -q %s_%s" % (track, pvar)   # Could also use P variable instead of eta
      command += " -n %s" % ntrvar
      command += " -t %s" % treename
      command += " -p %s_%s_gen_MagUp" % (track, var)
      command += " -c %s" % config
      command += " -d %s" % dataset
      command += " -i %s" % tmpinfile
      command += " -o %s" % tmpoutfile
      if seed : 
        command += " -s %d" % seed

      treename = output_tree
      tmpinfile = tmpoutfile
      if tmpoutfile == "tmp1.root" : 
        tmpoutfile = "tmp2.root"
      else : 
        tmpoutfile = "tmp1.root"

      print command
      os.system(command)

  os.system("rm %s" % tmpoutfile)
  os.system("mv %s %s" % (tmpinfile, output_file))
