#-- GAUDI jobOptions generated on Wed Jun  1 16:27:41 2016
#-- Contains event types : 
#--   13266007 - 28 files - 629199 events - 185.22 GBytes


#--  Extra information about the data processing phases:


#--  Processing Pass Step-127969 

#--  StepId : 127969 
#--  StepName : Reco14c for MC - 2012 
#--  ApplicationName : Brunel 
#--  ApplicationVersion : v43r2p11 
#--  OptionFiles : $APPCONFIGOPTS/Brunel/DataType-2012.py;$APPCONFIGOPTS/Brunel/MC-WithTruth.py;$APPCONFIGOPTS/Persistency/DST-multipleTCK-2012.py;$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py 
#--  DDDB : fromPreviousStep 
#--  CONDDB : fromPreviousStep 
#--  ExtraPackages : AppConfig.v3r218 
#--  Visible : Y 


#--  Processing Pass Step-127160 

#--  StepId : 127160 
#--  StepName : Stripping21-NoPrescalingFlagged for Sim08 - MU - Implicit merging. 
#--  ApplicationName : DaVinci 
#--  ApplicationVersion : v36r1p1 
#--  OptionFiles : $APPCONFIGOPTS/DaVinci/DV-Stripping21-Stripping-MC-NoPrescaling.py;$APPCONFIGOPTS/DaVinci/DV-RedoCaloPID-Stripping21.py;$APPCONFIGOPTS/DaVinci/DataType-2012.py;$APPCONFIGOPTS/DaVinci/InputType-DST.py;$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py 
#--  DDDB : fromPreviousStep 
#--  CONDDB : sim-20141210-1-vc-mu100 
#--  ExtraPackages : AppConfig.v3r205 
#--  Visible : Y 

from Gaudi.Configuration import * 
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000001_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000002_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000003_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000005_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000006_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000007_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000008_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000009_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000010_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000013_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000014_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000015_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000016_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000017_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000018_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000019_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000020_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000021_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000023_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000024_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000025_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000027_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000028_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000030_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000034_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000036_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000037_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051570/0000/00051570_00000038_2.AllStreams.dst'
], clear=True)
