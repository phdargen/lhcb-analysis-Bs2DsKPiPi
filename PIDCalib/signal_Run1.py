import os

# Dictionary of tracks with their PID variables, in the form {branch name}:{pidvars}
# For each track branch name, {pidvars} is a dictionary in the form {ntuple variable}:{pid config}, 
#   where 
#     {ntuple variable} is the name of the corresponding ntuple PID variable without branch name, 
#   and 
#     {pid_config} is the string describing the PID configuration. 
# Run PIDCorr.py without arguments to get the full list of PID configs
trackList = [ 
	# Bs->DsKpipi, Ds->KKpi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'K_plus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},

	# Bs->DsKpipi, Ds->pipipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus2_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
		# Bs->DsKpipi, Ds->Kpipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
		#
		#
		# Bs->Dspipipi, Ds->KKpi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'K_plus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},

		# Bs->Dspipipi, Ds->pipipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},          
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus2_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
		# Bs->Dspipipi, Ds->Kpipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},       
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
]


# List of input ROOT files with MC ntuples. Format: 
#   (inputfile, outputfile, dataset)
indir = "/auto/data/dargent/BsDsKpipi/Mini/MC/"
files = [
  (indir+"signal_Ds2KKpi_11.root",          indir+"/signal_Ds2KKpi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[0],"Gen"), 
  (indir+"signal_Ds2KKpi_11_PID_tmp.root",  indir+"/signal_Ds2KKpi_11_PID_tmp2.root", "MagUp_2011",   "gen_MagUp",  trackList[0],"Gen"), 
  (indir+"signal_Ds2KKpi_11_PID_tmp2.root", indir+"/signal_Ds2KKpi_11_PID_tmp3.root", "MagDown_2011", "corr_MagDown",trackList[0],"Corr"), 
  (indir+"signal_Ds2KKpi_11_PID_tmp3.root", indir+"/signal_Ds2KKpi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[0],"Corr"),
   
  (indir+"signal_Ds2KKpi_12.root",          indir+"/signal_Ds2KKpi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[0],"Gen"), 
  (indir+"signal_Ds2KKpi_12_PID_tmp.root",  indir+"/signal_Ds2KKpi_12_PID_tmp2.root", "MagUp_2012",   "gen_MagUp",  trackList[0],"Gen"), 
  (indir+"signal_Ds2KKpi_12_PID_tmp2.root", indir+"/signal_Ds2KKpi_12_PID_tmp3.root", "MagDown_2012", "corr_MagDown",trackList[0],"Corr"), 
  (indir+"signal_Ds2KKpi_12_PID_tmp3.root", indir+"/signal_Ds2KKpi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[0],"Corr"),
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

output_tree = input_tree.split("/")[-1]
treename = input_tree

simversion = "sim08"

for input_file, output_file, dataset, pidname, tracks, method in files : 
  tmpinfile = input_file
  tmpoutfile = "tmp3.root"
  for track, subst in tracks.iteritems() : 
    for var, config in subst.iteritems() : 
	if method == "Gen" :
		#      command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py"
		command = "python ./PIDGen.py"
		command += " -m %s_%s" % (track, ptvar)
		command += " -e %s_%s" % (track, etavar)
		##      command += " -q %s_%s" % (track, pvar)   # Could also use P variable instead of eta
		command += " -n %s" % ntrvar
		command += " -l %s" % config[1]
		command += " -z %s" % config[2]
		command += " -t %s" % treename
		command += " -p %s_%s_%s" % (track, var, pidname)
		command += " -c %s" % config[0]
		command += " -d %s" % dataset
		command += " -i %s" % tmpinfile
		command += " -o %s" % tmpoutfile
		if seed : 
			command += " -s %d" % seed

	else :
	    #command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDCorr.py"
            command = "python ./PIDCorr.py"
            command += " -m %s_%s" % (track, ptvar)
            command += " -e %s_%s" % (track, etavar)
            ## command += " -q %s_%s" % (track, pvar)   # Could also use P variable instead of eta
            command += " -n %s" % ntrvar
	    command += " -l %s" % config[1]
	    command += " -z %s" % config[2]
            command += " -t %s" % treename
	    command += " -p %s_%s_%s" % (track, var, pidname)
            command += " -s %s_%s" % (track, var)
            command += " -c %s" % config[0]
            command += " -d %s" % dataset
            command += " -i %s" % tmpinfile
            command += " -o %s" % tmpoutfile
            command += " -S %s" % simversion

      	treename = output_tree
      	tmpinfile = tmpoutfile
      	if tmpoutfile == "tmp3.root" : 
      	  tmpoutfile = "tmp4.root"
      	else : 
       	 tmpoutfile = "tmp3.root"

     	 print command
      	os.system(command)

  os.system("rm %s" % tmpoutfile)
  os.system("mv %s %s" % (tmpinfile, output_file))
