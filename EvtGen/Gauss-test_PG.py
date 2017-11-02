nEvts = 3000000
idFile = '/auto/data/dargent/BsDsKpipi/EvtGen/Gen_5'

from Gauss.Configuration import *

#--Generator phase, set random numbers
GaussGen = GenInit("GaussGen")
GaussGen.FirstEventNumber = 1
GaussGen.RunNumber        = 1090

LHCbApp().EvtMax = nEvts
Gauss().Production = 'PGUN'

#from Configurables import Gauss
#Gauss().Redecay['active'] = True
#Gauss().Redecay['N'] = 100
#Gauss().OutputType = 'NONE'
#Gauss().Histograms = 'NONE'
HistogramPersistencySvc().OutputFile = idFile+'-histos.root'
OutputStream("GaussTape").Output = "DATAFILE='PFN:%s.xgen' TYP='POOL_ROOTTREE' OPT='RECREATE'"%idFile

