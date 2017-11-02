#
# Options specific for a given job
# ie. setting of random number seed and name of output files
#

from Gauss.Configuration import *

#--Generator phase, set random numbers
GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1
GaussGen.RunNumber        = 1088

#--Number of events
nEvts = 150
LHCbApp().EvtMax = nEvts
Gauss().Production = 'PGUN'

#from Configurables import Gauss
#Gauss().Redecay['active'] = True
#Gauss().Redecay['N'] = 100


#Gauss().OutputType = 'NONE'
#Gauss().Histograms = 'NONE'
#--Set name of output files for given job (uncomment the lines)
#  Note that if you do not set it Gauss will make a name based on event type,
#  number of events and the date
idFile = 'Gen_5'
HistogramPersistencySvc().OutputFile = idFile+'-histos.root'
#
OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.xgen' TYP='POOL_ROOTTREE' OPT='RECREATE'"%idFile

#GenMonitor = GaudiSequencer( "GenMonitor" )
#SimMonitor = GaudiSequencer( "SimMonitor" )
#GenMonitor.Members += [ "GaussMonitor::CheckLifeTimeHepMC/HepMCLifeTime" ]
#SimMonitor.Members += [ "GaussMonitor::CheckLifeTimeMC/MCLifeTime" ]

