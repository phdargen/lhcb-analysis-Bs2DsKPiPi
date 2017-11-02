# file /afs/cern.ch/user/p/phdargen/cmtuser/Gauss_v49r7/Gen/DecFiles/options/13266007.py generated: Mon, 24 Jul 2017 14:38:16
#
# Event Type: 13266007
#
# ASCII decay Descriptor: {[[B_s0]nos -> (D_s- => K+ K- pi-) K+ pi- pi+]cc, [[B_s0]os -> (D_s+ => K- K+ pi+) K- pi+ pi-]cc}
#
#from Gaudi.Configuration import *
#importOptions( "$DECFILESROOT/options/TracksInAccWithMinP.py" )
#from Configurables import Generation
#Generation().EventType = 13266007
#Generation().SampleGenerationTool = "SignalRepeatedHadronization"
#from Configurables import SignalRepeatedHadronization
#Generation().addTool( SignalRepeatedHadronization )
#Generation().SignalRepeatedHadronization.ProductionTool = "PythiaProduction"
#from Configurables import ToolSvc
#from Configurables import EvtGenDecay
#ToolSvc().addTool( EvtGenDecay )
#ToolSvc().EvtGenDecay.UserDecayFile = "test.dec"
#Generation().SignalRepeatedHadronization.CutTool = "DaughtersInLHCbAndWithMinP"
#Generation().SignalRepeatedHadronization.SignalPIDList = [ 531,-531 ]

# Ad-hoc particle gun code
#from Configurables import ParticleGun
#pgun = ParticleGun("ParticleGun")
#pgun.SignalPdgCode = 531
#pgun.DecayTool = "EvtGenDecay"
#pgun.GenCutTool = "DaughtersInLHCbAndWithMinP"
#pgun.addTool( Generation().SignalRepeatedHadronization.DaughtersInLHCbAndWithMinP.clone() )

#from Configurables import FlatNParticles
#pgun.NumberOfParticlesTool = "FlatNParticles"
#pgun.addTool( FlatNParticles , name = "FlatNParticles" )

#from Configurables import MomentumSpectrum
#pgun.ParticleGunTool = "MomentumSpectrum"
#pgun.addTool( MomentumSpectrum , name = "MomentumSpectrum" )
#pgun.MomentumSpectrum.PdgCodes = [ 531,-531 ]
#pgun.MomentumSpectrum.InputFile = "$PGUNSDATAROOT/data/Ebeam4000GeV/MomentumSpectrum_531.root"
#pgun.MomentumSpectrum.BinningVariables = "pteta"
#pgun.MomentumSpectrum.HistogramPath = "h_pteta"

#from Configurables import BeamSpotSmearVertex
#pgun.addTool(BeamSpotSmearVertex, name="BeamSpotSmearVertex")
#pgun.VertexSmearingTool = "BeamSpotSmearVertex"
#pgun.EventType = 13266007


from Configurables import ParticleGun, MomentumRange, FlatNParticles, ToolSvc, EvtGenDecay
from GaudiKernel import SystemOfUnits

pgun = ParticleGun()
pgun.ParticleGunTool = "MomentumRange"
pgun.addTool( MomentumRange , name = "MomentumRange" )

pgun.NumberOfParticlesTool = "FlatNParticles"
pgun.addTool( FlatNParticles , name = "FlatNParticles" )

pgun.MomentumRange.PdgCodes = [ 531 , -531 ]

tsvc = ToolSvc()
tsvc.addTool( EvtGenDecay , name = "EvtGenDecay" )
ParticleGun().EventType = 13266007
tsvc.EvtGenDecay.UserDecayFile = "Bs_DsKpipi-DDalitz=DecProdCut_pCut1600MeV_PG.dec"
pgun.DecayTool = "EvtGenDecay"

pgun.MomentumRange.MomentumMin = 20.0*SystemOfUnits.GeV
pgun.MomentumRange.MomentumMax = 140.0*SystemOfUnits.GeV
