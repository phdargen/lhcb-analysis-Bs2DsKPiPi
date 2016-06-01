#-- GAUDI jobOptions generated on Wed Jun  1 16:25:37 2016
#-- Contains event types : 
#--   13266007 - 28 files - 628045 events - 184.65 GBytes


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


#--  Processing Pass Step-127159 

#--  StepId : 127159 
#--  StepName : Stripping21-NoPrescalingFlagged for Sim08 - MD - Implicit merging. 
#--  ApplicationName : DaVinci 
#--  ApplicationVersion : v36r1p1 
#--  OptionFiles : $APPCONFIGOPTS/DaVinci/DV-Stripping21-Stripping-MC-NoPrescaling.py;$APPCONFIGOPTS/DaVinci/DV-RedoCaloPID-Stripping21.py;$APPCONFIGOPTS/DaVinci/DataType-2012.py;$APPCONFIGOPTS/DaVinci/InputType-DST.py;$APPCONFIGOPTS/Persistency/Compression-ZLIB-1.py 
#--  DDDB : fromPreviousStep 
#--  CONDDB : sim-20141210-1-vc-md100 
#--  ExtraPackages : AppConfig.v3r205 
#--  Visible : Y 

from Gaudi.Configuration import * 
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles(['LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000001_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000002_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000003_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000004_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000006_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000008_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000010_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000011_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000012_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000013_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000014_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000015_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000016_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000017_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000018_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000019_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000020_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000021_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000022_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000024_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000025_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000026_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000027_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000028_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000030_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000032_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000033_2.AllStreams.dst',
'LFN:/lhcb/MC/2012/ALLSTREAMS.DST/00051572/0000/00051572_00000034_2.AllStreams.dst'
], clear=True)
