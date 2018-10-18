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
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'K_plus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},

	# Bs->DsKpipi, Ds->pipipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus2_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
		# Bs->DsKpipi, Ds->Kpipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'K_plus'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-2" , "100"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
		#
		#
		# Bs->Dspipipi, Ds->KKpi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},        
		'K_plus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},

		# Bs->Dspipipi, Ds->pipipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},          
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus2_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
		# Bs->Dspipipi, Ds->Kpipi
	{
		'pi_minus'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 	
			#         "PIDp"       : "pi_CombDLLp", 
		},
		'pi_plus1'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
			##    "PIDp"       : "K_CombDLLp", 
		},
		'pi_plus2'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "10"), 
				##    "PIDp"       : "pi_CombDLLp", 
		},       
		'pi_plus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20") , 
				##"PIDp"       : "K_CombDLLp", 
		},                    
		'K_minus_fromDs'   : {
				"PIDK"       : ("K_CombDLLK_Brunel" , "-10" , "100"), 
		##        "PIDp"       : "K_CombDLLp", 
		},            
		'pi_minus_fromDs'   : {
				"PIDK"       : ("pi_CombDLLK_Brunel" , "-100" , "20"), 
		##        "PIDp"       : "pi_CombDLLp", 
		}                   
	},
	
]


# List of input ROOT files with MC ntuples. Format: 
#   (inputfile, outputfile, dataset)
indir = "/eos/lhcb/user/p/phdargen/BsDsKpipi/" #/auto/data/dargent/BsDsKpipi/Mini/MC/"
files = [

  #(indir+"signal_Ds2KKpi_11.root",          indir+"/signal_Ds2KKpi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[0],"Gen"), 
  #(indir+"signal_Ds2KKpi_11_PID_tmp.root", indir+"/signal_Ds2KKpi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[0],"Corr"), 
  #(indir+"signal_Ds2KKpi_11_PID_tmp2.root", indir+"/signal_Ds2KKpi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[0],"Corr"),
   
  #(indir+"signal_Ds2KKpi_12.root",          indir+"/signal_Ds2KKpi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[0],"Gen"), 
  #(indir+"signal_Ds2KKpi_12_PID_tmp.root", indir+"/signal_Ds2KKpi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[0],"Corr"), 
  #(indir+"signal_Ds2KKpi_12_PID_tmp2.root", indir+"/signal_Ds2KKpi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[0],"Corr"),


  #(indir+"signal_Ds2pipipi_11.root",          indir+"/signal_Ds2pipipi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[1],"Gen"), 
  #(indir+"signal_Ds2pipipi_11_PID_tmp.root", indir+"/signal_Ds2pipipi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[1],"Corr"), 
  #(indir+"signal_Ds2pipipi_11_PID_tmp2.root", indir+"/signal_Ds2pipipi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[1],"Corr"),
   
  #(indir+"signal_Ds2pipipi_12.root",          indir+"/signal_Ds2pipipi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[1],"Gen"), 
  #(indir+"signal_Ds2pipipi_12_PID_tmp.root", indir+"/signal_Ds2pipipi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[1],"Corr"), 
  #(indir+"signal_Ds2pipipi_12_PID_tmp2.root", indir+"/signal_Ds2pipipi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[1],"Corr"),


  #(indir+"signal_Ds2Kpipi_11.root",          indir+"/signal_Ds2Kpipi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[2],"Gen"), 
  #(indir+"signal_Ds2Kpipi_11_PID_tmp.root", indir+"/signal_Ds2Kpipi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[2],"Corr"), 
  #(indir+"signal_Ds2Kpipi_11_PID_tmp2.root", indir+"/signal_Ds2Kpipi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[2],"Corr"),
   
  #(indir+"signal_Ds2Kpipi_12.root",          indir+"/signal_Ds2Kpipi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[2],"Gen"), 
  #(indir+"signal_Ds2Kpipi_12_PID_tmp.root", indir+"/signal_Ds2Kpipi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[2],"Corr"), 
  #(indir+"signal_Ds2Kpipi_12_PID_tmp2.root", indir+"/signal_Ds2Kpipi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[2],"Corr"),


  #(indir+"signal_Ds2KKpi_15.root",          indir+"/signal_Ds2KKpi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[0],"Gen"), 
  #(indir+"signal_Ds2KKpi_15_PID_tmp.root", indir+"/signal_Ds2KKpi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[0],"Corr"), 
  #(indir+"signal_Ds2KKpi_15_PID_tmp2.root", indir+"/signal_Ds2KKpi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[0],"Corr"),
   
  #(indir+"signal_Ds2KKpi_16.root",          indir+"/signal_Ds2KKpi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[0],"Gen"), 
  #(indir+"signal_Ds2KKpi_16_PID_tmp.root", indir+"/signal_Ds2KKpi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[0],"Corr"), 
  #(indir+"signal_Ds2KKpi_16_PID_tmp2.root", indir+"/signal_Ds2KKpi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[0],"Corr"),


  #(indir+"signal_Ds2pipipi_15.root",          indir+"/signal_Ds2pipipi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[1],"Gen"), 
  #(indir+"signal_Ds2pipipi_15_PID_tmp.root", indir+"/signal_Ds2pipipi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[1],"Corr"), 
  #(indir+"signal_Ds2pipipi_15_PID_tmp2.root", indir+"/signal_Ds2pipipi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[1],"Corr"),
   
  #(indir+"signal_Ds2pipipi_16.root",          indir+"/signal_Ds2pipipi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[1],"Gen"), 
  #(indir+"signal_Ds2pipipi_16_PID_tmp.root", indir+"/signal_Ds2pipipi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[1],"Corr"), 
  #(indir+"signal_Ds2pipipi_16_PID_tmp2.root", indir+"/signal_Ds2pipipi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[1],"Corr"),


  #(indir+"signal_Ds2Kpipi_15.root",          indir+"/signal_Ds2Kpipi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[2],"Gen"), 
  #(indir+"signal_Ds2Kpipi_15_PID_tmp.root", indir+"/signal_Ds2Kpipi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[2],"Corr"), 
  #(indir+"signal_Ds2Kpipi_15_PID_tmp2.root", indir+"/signal_Ds2Kpipi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[2],"Corr"),
   
  #(indir+"signal_Ds2Kpipi_16.root",          indir+"/signal_Ds2Kpipi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[2],"Gen"), 
  #(indir+"signal_Ds2Kpipi_16_PID_tmp.root", indir+"/signal_Ds2Kpipi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[2],"Corr"), 
  #(indir+"signal_Ds2Kpipi_16_PID_tmp2.root", indir+"/signal_Ds2Kpipi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[2],"Corr"),




  #(indir+"norm_Ds2KKpi_11.root",          indir+"/norm_Ds2KKpi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[3],"Gen"), 
  #(indir+"norm_Ds2KKpi_11_PID_tmp.root", indir+"/norm_Ds2KKpi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[3],"Corr"), 
  #(indir+"norm_Ds2KKpi_11_PID_tmp2.root", indir+"/norm_Ds2KKpi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[3],"Corr"),
   
  #(indir+"norm_Ds2KKpi_12.root",          indir+"/norm_Ds2KKpi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[3],"Gen"), 
  #(indir+"norm_Ds2KKpi_12_PID_tmp.root", indir+"/norm_Ds2KKpi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[3],"Corr"), 
  #(indir+"norm_Ds2KKpi_12_PID_tmp2.root", indir+"/norm_Ds2KKpi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[3],"Corr"),


  #(indir+"norm_Ds2pipipi_11.root",          indir+"/norm_Ds2pipipi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[4],"Gen"), 
  #(indir+"norm_Ds2pipipi_11_PID_tmp.root", indir+"/norm_Ds2pipipi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[4],"Corr"), 
  #(indir+"norm_Ds2pipipi_11_PID_tmp2.root", indir+"/norm_Ds2pipipi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[4],"Corr"),
   
  #(indir+"norm_Ds2pipipi_12.root",          indir+"/norm_Ds2pipipi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[4],"Gen"), 
  #(indir+"norm_Ds2pipipi_12_PID_tmp.root", indir+"/norm_Ds2pipipi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[4],"Corr"), 
  #(indir+"norm_Ds2pipipi_12_PID_tmp2.root", indir+"/norm_Ds2pipipi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[4],"Corr"),


  #(indir+"norm_Ds2Kpipi_11.root",          indir+"/norm_Ds2Kpipi_11_PID_tmp.root",  "MagDown_2011", "gen_MagDown",trackList[5],"Gen"), 
  #(indir+"norm_Ds2Kpipi_11_PID_tmp.root", indir+"/norm_Ds2Kpipi_11_PID_tmp2.root", "MagDown_2011", "corr_MagDown",trackList[5],"Corr"), 
  #(indir+"norm_Ds2Kpipi_11_PID_tmp2.root", indir+"/norm_Ds2Kpipi_11_PID.root",      "MagUp_2011",   "corr_MagUp",  trackList[5],"Corr"),
   
  #(indir+"norm_Ds2Kpipi_12.root",          indir+"/norm_Ds2Kpipi_12_PID_tmp.root",  "MagDown_2012", "gen_MagDown",trackList[5],"Gen"), 
  #(indir+"norm_Ds2Kpipi_12_PID_tmp.root", indir+"/norm_Ds2Kpipi_12_PID_tmp2.root", "MagDown_2012", "corr_MagDown",trackList[5],"Corr"), 
  #(indir+"norm_Ds2Kpipi_12_PID_tmp2.root", indir+"/norm_Ds2Kpipi_12_PID.root",      "MagUp_2012",   "corr_MagUp",  trackList[5],"Corr"),


  #(indir+"norm_Ds2KKpi_15.root",          indir+"/norm_Ds2KKpi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[3],"Gen"), 
  #(indir+"norm_Ds2KKpi_15_PID_tmp.root", indir+"/norm_Ds2KKpi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[3],"Corr"), 
  #(indir+"norm_Ds2KKpi_15_PID_tmp2.root", indir+"/norm_Ds2KKpi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[3],"Corr"),
   
  #(indir+"norm_Ds2KKpi_16.root",          indir+"/norm_Ds2KKpi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[3],"Gen"), 
  #(indir+"norm_Ds2KKpi_16_PID_tmp.root", indir+"/norm_Ds2KKpi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[3],"Corr"), 
  #(indir+"norm_Ds2KKpi_16_PID_tmp2.root", indir+"/norm_Ds2KKpi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[3],"Corr"),


  (indir+"norm_Ds2pipipi_15.root",          "norm_Ds2pipipi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[4],"Gen"), 
  ("norm_Ds2pipipi_15_PID_tmp.root", "norm_Ds2pipipi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[4],"Corr"), 
  ("norm_Ds2pipipi_15_PID_tmp2.root", "norm_Ds2pipipi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[4],"Corr"),
   
  #(indir+"norm_Ds2pipipi_16.root",          indir+"/norm_Ds2pipipi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[4],"Gen"), 
  #(indir+"norm_Ds2pipipi_16_PID_tmp.root", indir+"/norm_Ds2pipipi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[4],"Corr"), 
  #(indir+"norm_Ds2pipipi_16_PID_tmp2.root", indir+"/norm_Ds2pipipi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[4],"Corr"),


  #(indir+"norm_Ds2Kpipi_15.root",          indir+"/norm_Ds2Kpipi_15_PID_tmp.root",  "MagDown_2015", "gen_MagDown",trackList[5],"Gen"), 
  #(indir+"norm_Ds2Kpipi_15_PID_tmp.root", indir+"/norm_Ds2Kpipi_15_PID_tmp2.root", "MagDown_2015", "corr_MagDown",trackList[5],"Corr"), 
  #(indir+"norm_Ds2Kpipi_15_PID_tmp2.root", indir+"/norm_Ds2Kpipi_15_PID.root",      "MagUp_2015",   "corr_MagUp",  trackList[5],"Corr"),
   
  #(indir+"norm_Ds2Kpipi_16.root",          indir+"/norm_Ds2Kpipi_16_PID_tmp.root",  "MagDown_2016", "gen_MagDown",trackList[5],"Gen"), 
  #(indir+"norm_Ds2Kpipi_16_PID_tmp.root", indir+"/norm_Ds2Kpipi_16_PID_tmp2.root", "MagDown_2016", "corr_MagDown",trackList[5],"Corr"), 
  #(indir+"norm_Ds2Kpipi_16_PID_tmp2.root", indir+"/norm_Ds2Kpipi_16_PID.root",      "MagUp_2016",   "corr_MagUp",  trackList[5],"Corr"),


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

simversion = "run2"

# make sure we don't overwrite local files and prefix them with random strings
# IF ON LXPLUS: if /tmp exists and is accessible, use for faster processing
# IF NOT: use /tmp if you have enough RAM
# temp_folder = '/tmp'
# ELSE: use current folder
temp_folder = '.'
import string
import random
rand_string = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10))  # get 10 random chars for temp_file prefix
temp_file_prefix = temp_folder + '/' + rand_string  # prefix temp files with folder and unique ID


for input_file, output_file, dataset, pidname, tracks, method in files : 
  tmpinfile = input_file
  tmpoutfile = "%s_tmp1.root" % temp_file_prefix
  for track, subst in tracks.iteritems() : 
    for var, config in subst.iteritems() : 
	if method == "Gen" :
		#      command = "python $PIDPERFSCRIPTSROOT/scripts/python/PIDGenUser/PIDGen.py"
		command = "python ./PIDGen_new.py"
		command += " -m %s_%s" % (track, ptvar)
		command += " -e %s_%s" % (track, etavar)
		##      command += " -q %s_%s" % (track, pvar)   # Could also use P variable instead of eta
		command += " -n %s" % ntrvar
		#command += " -l %s" % config[1]
		#command += " -z %s" % config[2]
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
            command = "python ./PIDCorr_new.py"
            command += " -m %s_%s" % (track, ptvar)
            command += " -e %s_%s" % (track, etavar)
            ## command += " -q %s_%s" % (track, pvar)   # Could also use P variable instead of eta
            command += " -n %s" % ntrvar
            #command += " -l %s" % config[1]
            #command += " -z %s" % config[2]
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
      	if 'tmp1' in tmpoutfile:
        	tmpoutfile = tmpoutfile.replace('tmp1', 'tmp2')
      	else :
        	tmpoutfile = tmpoutfile.replace('tmp2', 'tmp1')
     	 
	print command
      	os.system(command)


  print("mv %s %s" % (tmpinfile, output_file))
  os.system("mv %s %s" % (tmpinfile, output_file))
  print("rm %s" % tmpoutfile)
  os.system("rm %s" % tmpoutfile)
  
  
  
for input_file, output_file, dataset, pidname, tracks, method in files :
  if 'tmp' in output_file:
    print("rm %s" % output_file)
    os.system("rm %s" % output_file)