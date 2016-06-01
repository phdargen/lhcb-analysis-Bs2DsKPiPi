#-- GAUDI jobOptions generated on Wed Jun  1 16:23:11 2016
#-- Contains event types : 
#--   13266006 - 27 files - 605946 events - 178.78 GBytes


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
IOHelper('ROOT').inputFiles(['LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000001_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000002_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000003_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000004_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000005_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000006_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000008_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000009_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000010_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000011_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000012_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000013_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000014_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000015_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000016_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000017_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000018_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000019_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000020_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000021_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000023_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000026_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000027_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000028_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000029_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000031_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051578/0000/00051578_00000032_2.AllStreams.dst'
], clear=True)
