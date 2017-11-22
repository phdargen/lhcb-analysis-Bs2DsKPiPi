datafile = "Gen.xgen"
tuplefile = "Gen.root"


from Configurables import (
    DaVinci,
    EventSelector,
    PrintMCTree,
    MCDecayTreeTuple, TupleToolDalitz
)
from DecayTreeTuple.Configuration import *

"""Configure the variables below with:                                                                                                                                                                 
decay: Decay you want to inspect, using 'newer' LoKi decay descriptor syntax,                                                                                                                          
decay_heads: Particles you'd like to see the decay tree of,                                                                                                                                            
datafile: Where the file created by the Gauss generation phase is, and                                                                                                                                 
year: What year the MC is simulating.                                                                                                                                                                  
"""

# https://twiki.cern.ch/twiki/bin/view/LHCb/FAQ/LoKiNewDecayFinders                                                                                                                                    
decay = "[[B_s0]cc ==> ^pi+ ^pi+ ^pi- ^(D_s- ==> ^K+ ^K- ^pi-) ]CC"
decay_heads = ["B_s0", "B_s~0"]
year = 2012

# For a quick and dirty check, you don't need to edit anything below here.                                                                                                                             
##########################################################################                                                                                                                             

# Create an MC DTT containing any candidates matching the decay descriptor                                                                                                                             
#mctuple = MCDecayTreeTuple("MCDecayTreeTuple")
#mctuple.Decay = decay
#mctuple.ToolList = [
#    "MCTupleToolHierarchy",
#    "LoKi::Hybrid::MCTupleTool/LoKi_Photos"
#]
# Add a 'number of photons' branch                                                                                                                                                                     
#mctuple.addTupleTool("MCTupleToolKinematic").Verbose = True
#mctuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_Photos").Variables = {
#    "nPhotos": "MCNINTREE(('gamma' == MCABSID))"
#}

# Print the decay tree for any particle in decay_heads                                                                                                                                                 
#printMC = PrintMCTree()
#printMC.ParticleNames = decay_heads

# Name of the .xgen file produced by Gauss                                                                                                                                                             
#EventSelector().Input = ["DATAFILE='{0}' TYP='POOL_ROOTTREE' Opt='READ'".format(datafile)]

# Create an MC DTT containing any candidates matching the decay descriptor                                                                                                                             
mctuple = MCDecayTreeTuple("MCDecayTreeTuple")
mctuple.Decay = decay
mctuple.ToolList = [
    "MCTupleToolHierarchy",
    "LoKi::Hybrid::MCTupleTool/LoKi_Photos"
]
# Add a 'number of photons' branch                                                                                                                                                                     
mctuple.addTupleTool("MCTupleToolKinematic").Verbose = True
mctuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_Photos").Variables = {
    "nPhotos": "MCNINTREE(('gamma' == MCABSID))"
}


mctuple.addTupleTool("LoKi::Hybrid::MCTupleTool/Loki_All")
mctuple.Loki_All.Variables = {'TRUEID' : 'MCID'}

# Print the decay tree for any particle in decay_heads                                                                                                                                                 
printMC = PrintMCTree()
printMC.ParticleNames = decay_heads

# Name of the .xgen file produced by Gauss                                                                                                                                                             
EventSelector().Input = ["DATAFILE='{0}' TYP='POOL_ROOTTREE' Opt='READ'".format(datafile)]

# Configure DaVinci                                                                                                                                                                                    
DaVinci().TupleFile = tuplefile
DaVinci().Simulation = True
DaVinci().Lumi = False
DaVinci().DataType = str(year)
DaVinci().UserAlgorithms = [mctuple]
#DaVinci().UserAlgorithms = [printMC, mctuple]
