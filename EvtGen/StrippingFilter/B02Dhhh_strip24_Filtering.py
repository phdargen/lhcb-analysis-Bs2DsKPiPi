"""
Stripping filtering file for B0 -> D0hhh MC
@author Philippe d'Argent
@date   2017-12-07
"""
#stripping version --> DaVinci v38r1p1
stripping='stripping24r1'

#use CommonParticlesArchive
from CommonParticlesArchive import CommonParticlesArchiveConf
CommonParticlesArchiveConf().redirect(stripping)

from Gaudi.Configuration import *
MessageSvc().Format = "% F%30W%S%7W%R%T %0W%M"

#Raw event juggler to split Other/RawEvent into Velo/RawEvent and Tracker/RawEvent
#
from Configurables import RawEventJuggler
juggler = RawEventJuggler( DataOnDemand=True, Input=2.0, Output=4.0 )
#
# Build the streams and stripping object
#
from StrippingConf.Configuration import StrippingConf, StrippingStream
from StrippingSettings.Utils import strippingConfiguration
from StrippingArchive.Utils import buildStreams, cloneLinesFromStream
from StrippingArchive import strippingArchive


#get the configuration dictionary from the database
config  = strippingConfiguration(stripping)
#get the line builders from the archive
archive = strippingArchive(stripping)

streams = buildStreams(stripping = config, archive = archive)
# Select my lines
AllStreams = StrippingStream("B02D0hhh.Strip")
linesToAdd = []
MyLines = [ 'StrippingB02DKPiPiD2HHHPIDBeauty2CharmLine	',
            'StrippingB02DPiPiPiD2HHHPIDBeauty2CharmLine'
          ]

for stream in streams:
    if 'Bhadron' in stream.name():
        for line in stream.lines:
            line._prescale = 1.0
            if line.name() in MyLines:
                linesToAdd.append(line)
AllStreams.appendLines(linesToAdd)



sc = StrippingConf( Streams = [ AllStreams ],
                    MaxCandidates = 2000,
                    TESPrefix = 'Strip'
                    )


AllStreams.sequence().IgnoreFilterPassed = False # so that we do not get all events written out

#

# Configuration of SelDSTWriter
#
enablePacking = True

from DSTWriters.microdstelements import *
from DSTWriters.Configuration import (SelDSTWriter,
                                      stripDSTStreamConf,
                                      stripDSTElements
                                      )

SelDSTWriterElements = {
    'default'               : stripDSTElements(pack=enablePacking)
    }

SelDSTWriterConf = {
    'default'               : stripDSTStreamConf(pack=enablePacking)
    }
#Items that might get lost when running the CALO+PROTO ReProcessing in DV
caloProtoReprocessLocs = [ "/Event/pRec/ProtoP#99", "/Event/pRec/Calo#99" ]

# Make sure they are present on full DST streams
SelDSTWriterConf['default'].extraItems += caloProtoReprocessLocs

dstWriter = SelDSTWriter( "MyDSTWriter",
                          StreamConf = SelDSTWriterConf,
                          MicroDSTElements = SelDSTWriterElements,
                          OutputFileSuffix ='Filtered',
                          SelectionSequences = sc.activeStreams()
                          )

# Add stripping TCK 
from Configurables import StrippingTCK
stck = StrippingTCK(HDRLocation = '/Event/Strip/Phys/DecReports', TCK=0x38162410)# TCK = 0xVVVVSSSS, where VVVV is the DaVinci version, and SSSS the Stripping number

#
# DaVinci Configuration
#
from Configurables import DaVinci
DaVinci().Simulation = True
DaVinci().EvtMax = -1                     # Number of events
DaVinci().HistogramFile = "DVHistos.root"
DaVinci().appendToMainSequence( [ sc.sequence() ] )
DaVinci().appendToMainSequence( [ stck ] )
DaVinci().appendToMainSequence( [ dstWriter.sequence() ] )
DaVinci().ProductionType = "Stripping"

# Change the column size of Timing table
from Configurables import TimingAuditor, SequencerTimerTool
TimingAuditor().addTool(SequencerTimerTool,name="TIMER")
TimingAuditor().TIMER.NameSize = 60
