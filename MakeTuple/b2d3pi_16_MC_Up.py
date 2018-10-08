from Gaudi.Configuration import *
from Configurables import DaVinci
from Configurables import CombineParticles

############# Global settings
year = "2016"
data = False
down = False
stream = "B02D0hhh.Strip"
if (data):
    stream = "BhadronCompleteEvent"

# Event filter
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence
from Configurables import FilterDesktop
line = 'B02DPiPiPiD2HHHPIDBeauty2CharmLine'
inputs = '/Event/'+stream+'/Phys/{0}/Particles'.format(line)
reqsel = AutomaticData(Location = inputs)
Bs2DsXSel = FilterDesktop("Bs2DsXSel", Code = "(INTREE((ABSID=='D+')&(M>1888))) & (INTREE((ABSID=='D+')&(M<2048))) & (INTREE((ABSID=='B0')&(M>5000))) & (INTREE((ABSID=='B0')&(M<6000))) & (INTREE((ABSID=='a_1(1260)+')&(M<3000)))")
MyFilterSel = Selection("MyFilterSel", Algorithm = Bs2DsXSel, RequiredSelections = [reqsel])

from Configurables import CheckPV
checkPVs = CheckPV("checkPVs")
checkPVs.MinPVs = 1
checkPVs.MaxPVs = -1

triggerlines_Run1 = [
		"L0HadronDecision", 
                "L0MuonDecision",
		"L0DiMuonDecision",
		"L0ElectronDecision",
		"L0PhotonDecision",
                "Hlt1TrackAllL0Decision", 
                'Hlt2IncPhiDecision', 
                'Hlt2Topo2BodyBBDTDecision', 
                'Hlt2Topo3BodyBBDTDecision', 
                'Hlt2Topo4BodyBBDTDecision' 
]

triggerlines_Run2 = [
		"L0HadronDecision", 
                "L0MuonDecision",
		"L0DiMuonDecision",
		"L0ElectronDecision",
		"L0PhotonDecision",
                "Hlt1TrackMVADecision",
                "Hlt1TwoTrackMVADecision",
		'Hlt2IncPhiDecision', 
		"Hlt2PhiIncPhiDecision",
		'Hlt2Topo2BodyDecision',
                'Hlt2Topo3BodyDecision',
                'Hlt2Topo4BodyDecision'
]

triggerlines = triggerlines_Run2

#from Configurables import (BTagging)   
#Btag = BTagging("BTagging")
#Btag.Inputs = [ "/Event/"+stream+"/Phys/"+line+"/Particles"]
#Btag.OutputLevel    = 6

#from FlavourTagging.Tunings import TuneTool
#def configureTaggingTools(BTuple, BsBd):
	
	    #from Configurables import BTagging, BTaggingTool
	    #toolname = BsBd+"TaggingTool"
	    #BTuple.TaggingToolName ="BTaggingTool/"+toolname
	    #BTuple.ExtraName = toolname
	    #if(data):
		    #TuneTool(BTuple,"Stripping21",toolname)
	    #else :
		    #TuneTool(BTuple,"Stripping21_MC",toolname)
	    #BTuple.Verbose = True
	    #BTuple.addTool(BTaggingTool,name=toolname)
	    #BTuple.__getattr__(toolname).ForceSignalID = BsBd


###############   Pre Filter, does not really do much except choose only candidates passing the Stripping line, maybe beneficial to performance
from Configurables import LoKi__HDRFilter as StripFilter
stripFilter = StripFilter( 'stripPassFilter',\
                           Code = "HLT_PASS('StrippingB02DPiPiPiD2HHHPIDBeauty2CharmLineDecision')",\
                           Location= "/Event/Strip/Phys/DecReports")

############# DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import TupleToolTISTOS
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence
from Configurables import PrintDecayTree, PrintDecayTreeTool
##subpid stuff
#from Configurables import SubPIDMMFilter
from Configurables import SubstitutePID,BTaggingTool
from Configurables import TupleToolDecayTreeFitter, TupleToolTrackIsolation, TupleToolTagging, TupleToolRecoStats, TupleToolKinematic, TupleToolGeometry,TupleToolVtxIsoln
from Configurables import LoKi__Hybrid__TupleTool


#B0 -> (D- -> K K pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipiTuple = DecayTreeTuple("Bs2Dspipipi_Ds2KKpi_Tuple")
b2dkpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^K+ ^K- ^pi-) ^(a_1(1260)+ -> ^pi+ ^pi+ ^pi-)]CC"
b2dkpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> K+ K- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC)" ,
"a_1_1260_plus" : "[[B0]cc -> (D- -> K+ K- pi-) ^(a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_plus1" : "[[B0]cc -> (D- -> K+ K- pi-) (a_1(1260)+ -> ^pi+ pi+ pi-)  ]CC",
"pi_plus2" : "[[B0]cc -> (D- -> K+ K- pi-) (a_1(1260)+ -> pi+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> K+ K- pi-) (a_1(1260)+ -> pi+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> K+ K- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"K_plus_fromDs" : "[[B0]cc -> (D- -> ^K+ K- pi-) (a_1(1260)+ -> pi+ pi+ pi-)  ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> K+ ^K- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> K+ K- ^pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC"
}
b2dkpipiTuple.ReFitPVs = True

#config tools
b2dkpipiTuple.ToolList +=  ["TupleToolGeometry", \
                            "TupleToolKinematic", \
                            "TupleToolPrimaries", \
                            "TupleToolEventInfo", \
                            "TupleToolTrackInfo", \
                            "TupleToolRecoStats", \
                            #"TupleToolAngles", \
                            "TupleToolPid", \
                            #"TupleToolPhotonInfo", \
                            "TupleToolTrackIsolation",
			    "TupleToolVtxIsoln" 
			    ]

from Configurables import MCMatchObjP2MCRelator

if (data==False):
    b2dkpipiTuple.ToolList +=  [
                           # "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo",
			     "TupleToolPhotonInfo"
				]
				
    # Add TupleToolMCTruth with fix for "Fatal error No valid data at '/Event/Hlt2/Long/Protos'"
    default_rel_locs = MCMatchObjP2MCRelator().getDefaultProperty('RelTableLocations')
    rel_locs = [loc for loc in default_rel_locs if 'Turbo' not in loc]

    mctruth = b2dkpipiTuple.addTupleTool('TupleToolMCTruth')
    mctruth.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
    ]
    mctruth.addTool(MCMatchObjP2MCRelator)
    mctruth.MCMatchObjP2MCRelator.RelTableLocations = rel_locs			
    #MCTruth = TupleToolMCTruth()
    #b2dkpipiTuple.addTool(MCTruth) 

b2dkpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipiTuple.TupleToolTrackIsolation.FillAsymmetry = True
b2dkpipiTuple.TupleToolTrackIsolation.FillDeltaAngles = False
b2dkpipiTuple.TupleToolTrackIsolation.MinConeAngle = 1.0

b2dkpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipiTuple.TupleToolGeometry.Verbose = True
b2dkpipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipiTuple.TupleToolGeometry.FillMultiPV = True

b2dkpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipiTuple.UseLabXSyntax = True                          
b2dkpipiTuple.RevertToPositiveID = False

b2dkpipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]        
b2dkpipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+'  }
b2dkpipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+", ]  
b2dkpipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipiTuple.Bs.DTF.Verbose = True

b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]
b2dkpipiTuple.Bs.B0DTF.constrainToOriginVertex = True
b2dkpipiTuple.Bs.B0DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
b2dkpipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0", "D_s+", "B~0" ]  
b2dkpipiTuple.Bs.B0DTF.UpdateDaughters = True
b2dkpipiTuple.Bs.B0DTF.Verbose = True

b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]         
b2dkpipiTuple.Bs.BsDTF.constrainToOriginVertex = True
b2dkpipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
b2dkpipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
b2dkpipiTuple.Bs.BsDTF.UpdateDaughters = True
b2dkpipiTuple.Bs.BsDTF.Verbose = True

b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]         
b2dkpipiTuple.Bs.PV.constrainToOriginVertex = True 
b2dkpipiTuple.Bs.PV.UpdateDaughters = True
b2dkpipiTuple.Bs.PV.Verbose = True

LoKiTool = b2dkpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { "ETA" : "ETA" };

LoKiToolBs = b2dkpipiTuple.Bs.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolBs")
LoKiToolBs.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                    	};

b2dkpipiTuple.addTool(TupleToolDecay, name="Ds")
LoKiToolDs = b2dkpipiTuple.Ds.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolDs")
LoKiToolDs.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"
                       };

b2dkpipiTuple.addTool(TupleToolDecay, name="a_1_1260_plus")
LoKiToola_1_1260_plus = b2dkpipiTuple.a_1_1260_plus.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToola_1_1260_plus")
LoKiToola_1_1260_plus.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"
                       };

#tagging config
from Configurables import BTaggingTool
tt_tagging = b2dkpipiTuple.addTupleTool("TupleToolTagging") 
tt_tagging.Verbose = True
tt_tagging.AddMVAFeatureInfo = False
tt_tagging.AddTagPartsInfo = False 

btagtool = tt_tagging.addTool(BTaggingTool , name = "MyBTaggingTool")
from FlavourTagging.Tunings import applyTuning as applyFTTuning # pick the right tuning here ...
applyFTTuning(btagtool , tuning_version="Summer2017Optimisation_v4_Run2")
tt_tagging.TaggingToolName = btagtool.getFullName ()

#trigger config
b2dkpipitt = b2dkpipiTuple.addTupleTool(TupleToolTISTOS)
b2dkpipitt.TriggerList = triggerlines
b2dkpipitt.FillL0 = True
b2dkpipitt.FillHlt1 = True
b2dkpipitt.FillHlt2 = True
b2dkpipitt.Verbose = True
b2dkpipitt.VerboseL0 = True
b2dkpipitt.VerboseHlt1 = True
b2dkpipitt.VerboseHlt2 = True

#b2dkpipiTuple.Bs.addTool(TupleToolTISTOS,name="TisTosB")
#b2dkpipiTuple.Bs.TisTosB.Verbose=True
#b2dkpipiTuple.Bs.TisTosB.TriggerList= triggers

#printer
#b2dkpipiprinter = PrintDecayTree("PrintB2Dkpipi")
#b2dkpipiprinter.addTool( PrintDecayTreeTool, name = "PrintDecay" )
#b2dkpipiprinter.PrintDecay.Information = "Name M P Px Py Pz Pt chi2"
#b2dkpipiprinter.Inputs = [ makebu2kpipimumuseq.outputLocation() ]

#main sequence
makeb2dkpipiseq = SelectionSequence("makeb2dkpipiseq", TopSelection = MyFilterSel)
b2dkpipiTuple.Inputs = [makeb2dkpipiseq.outputLocation()]
b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
#b2dkpipiseq.RootInTES = '/Event/{0}'.format(stream)
b2dkpipiseq.Members += [makeb2dkpipiseq.sequence(),b2dkpipiTuple]

#b2dkpipiTuple.Inputs = [ "/Event/"+stream+"/Phys/B02DPiPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
#b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
##bu2kpipimumuseq.Members += [makebu2kpipimumuseq.sequence(), bu2kpipimumuprinter, bu2kpipimumutuple]
#b2dkpipiseq.Members += [b2dkpipiTuple]


#
#B0 -> (D- -> pi pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2pipipiTuple = DecayTreeTuple("Bs2Dspipipi_Ds2pipipi_Tuple")
b2dkpipi_d2pipipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^pi-) ^(a_1(1260)+ -> ^pi+ ^pi+ ^pi-)]CC"
b2dkpipi_d2pipipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC)" ,
"a_1_1260_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-) ^(a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_plus1" : "[[B0]cc -> (D- -> pi+ pi- pi-) (a_1(1260)+ -> ^pi+ pi+ pi-)  ]CC",
"pi_plus2" : "[[B0]cc -> (D- -> pi+ pi- pi-) (a_1(1260)+ -> pi+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- pi-) (a_1(1260)+ -> pi+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- pi-) (a_1(1260)+ -> pi+ pi+ pi-)  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_minus2_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^pi-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC"
}
b2dkpipi_d2pipipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2pipipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              "TupleToolRecoStats", \
                              #"TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolTrackIsolation",
			       "TupleToolVtxIsoln"
			      ]

if (data==False):
    b2dkpipi_d2pipipiTuple.ToolList +=  [
                           # "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo",
			     "TupleToolPhotonInfo"
				]
				
    # Add TupleToolMCTruth with fix for "Fatal error No valid data at '/Event/Hlt2/Long/Protos'"
    default_rel_locs = MCMatchObjP2MCRelator().getDefaultProperty('RelTableLocations')
    rel_locs = [loc for loc in default_rel_locs if 'Turbo' not in loc]

    mctruth = b2dkpipi_d2pipipiTuple.addTupleTool('TupleToolMCTruth')
    mctruth.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
    ]
    mctruth.addTool(MCMatchObjP2MCRelator)
    mctruth.MCMatchObjP2MCRelator.RelTableLocations = rel_locs	

b2dkpipi_d2pipipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipi_d2pipipiTuple.TupleToolTrackIsolation.FillAsymmetry = True
b2dkpipi_d2pipipiTuple.TupleToolTrackIsolation.FillDeltaAngles = False
b2dkpipi_d2pipipiTuple.TupleToolTrackIsolation.MinConeAngle = 1.0

b2dkpipi_d2pipipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2pipipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2pipipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2pipipiTuple.TupleToolGeometry.Verbose = True
b2dkpipi_d2pipipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipi_d2pipipiTuple.TupleToolGeometry.FillMultiPV = True

b2dkpipi_d2pipipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipi_d2pipipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipi_d2pipipiTuple.UseLabXSyntax = True                          
b2dkpipi_d2pipipiTuple.RevertToPositiveID = False
            
	    
b2dkpipi_d2pipipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]      
b2dkpipi_d2pipipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
b2dkpipi_d2pipipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+" ]  
b2dkpipi_d2pipipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.DTF.Verbose = True

b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]        
b2dkpipi_d2pipipiTuple.Bs.B0DTF.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.B0DTF.Substitutions = {'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
b2dkpipi_d2pipipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0","D_s+", "B~0"]  
b2dkpipi_d2pipipiTuple.Bs.B0DTF.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.B0DTF.Verbose = True

b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]        
b2dkpipi_d2pipipiTuple.Bs.BsDTF.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
b2dkpipi_d2pipipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
b2dkpipi_d2pipipiTuple.Bs.BsDTF.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.BsDTF.Verbose = True

b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]        
b2dkpipi_d2pipipiTuple.Bs.PV.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.PV.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.PV.Verbose = True


LoKiTool = b2dkpipi_d2pipipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { "ETA" : "ETA" };

LoKiToolBs = b2dkpipi_d2pipipiTuple.Bs.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolBs")
LoKiToolBs.Variables = { "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()" };

b2dkpipi_d2pipipiTuple.addTool(TupleToolDecay, name="Ds")
LoKiToolDs = b2dkpipi_d2pipipiTuple.Ds.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolDs")
LoKiToolDs.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"};

b2dkpipi_d2pipipiTuple.addTool(TupleToolDecay, name="a_1_1260_plus")
LoKiToola_1_1260_plus = b2dkpipi_d2pipipiTuple.a_1_1260_plus.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToola_1_1260_plus")
LoKiToola_1_1260_plus.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)" };

#tagging config
from Configurables import BTaggingTool
tt_tagging = b2dkpipi_d2pipipiTuple.addTupleTool("TupleToolTagging") 
tt_tagging.Verbose = True
tt_tagging.AddMVAFeatureInfo = False
tt_tagging.AddTagPartsInfo = False 

btagtool = tt_tagging.addTool(BTaggingTool , name = "MyBTaggingTool")
from FlavourTagging.Tunings import applyTuning as applyFTTuning # pick the right tuning here ...
applyFTTuning(btagtool , tuning_version="Summer2017Optimisation_v4_Run2")
tt_tagging.TaggingToolName = btagtool.getFullName ()

#trigger config
b2dkpipi_d2pipipitt = b2dkpipi_d2pipipiTuple.addTupleTool(TupleToolTISTOS)
b2dkpipi_d2pipipitt.TriggerList = triggerlines
b2dkpipi_d2pipipitt.FillL0 = True
b2dkpipi_d2pipipitt.FillHlt1 = True
b2dkpipi_d2pipipitt.FillHlt2 = True
b2dkpipi_d2pipipitt.Verbose = True
b2dkpipi_d2pipipitt.VerboseL0 = True
b2dkpipi_d2pipipitt.VerboseHlt1 = True
b2dkpipi_d2pipipitt.VerboseHlt2 = True

#main sequence
makeb2dkpipi_d2pipipiseq = SelectionSequence("makeb2dkpipi_d2pipipiseq", TopSelection = MyFilterSel)
b2dkpipi_d2pipipiTuple.Inputs = [makeb2dkpipi_d2pipipiseq.outputLocation()]
b2dkpipi_d2pipipiseq = GaudiSequencer("B2dkpipi_d2pipipiSeq")
#b2dkpipiseq.RootInTES = '/Event/{0}'.format(stream)
b2dkpipi_d2pipipiseq.Members += [makeb2dkpipi_d2pipipiseq.sequence(),b2dkpipi_d2pipipiTuple]


#
#B0 -> (D- -> K pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2KpipiTuple = DecayTreeTuple("Bs2Dspipipi_Ds2Kpipi_Tuple")
b2dkpipi_d2KpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^K-) ^(a_1(1260)+ -> ^pi+ ^pi+ ^pi-)]CC"
b2dkpipi_d2KpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- K-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC)" ,
"a_1_1260_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) ^(a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_plus1" : "[[B0]cc -> (D- -> pi+ pi- K-) (a_1(1260)+ -> ^pi+ pi+ pi-)  ]CC",
"pi_plus2" : "[[B0]cc -> (D- -> pi+ pi- K-) (a_1(1260)+ -> pi+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- K-) (a_1(1260)+ -> pi+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- K-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- K-) (a_1(1260)+ -> pi+ pi+ pi-)  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- K-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^K-) (a_1(1260)+ -> pi+ pi+ pi-) ]CC"
}
b2dkpipi_d2KpipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2KpipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              "TupleToolRecoStats", \
                              #"TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolTrackIsolation",
	                      "TupleToolVtxIsoln"
                              #"TupleToolTagging" 
			      ]

if (data==False):
    b2dkpipi_d2KpipiTuple.ToolList +=  [
                           # "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo",
			     "TupleToolPhotonInfo"
				]
				
    # Add TupleToolMCTruth with fix for "Fatal error No valid data at '/Event/Hlt2/Long/Protos'"
    default_rel_locs = MCMatchObjP2MCRelator().getDefaultProperty('RelTableLocations')
    rel_locs = [loc for loc in default_rel_locs if 'Turbo' not in loc]

    mctruth = b2dkpipi_d2KpipiTuple.addTupleTool('TupleToolMCTruth')
    mctruth.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
    ]
    mctruth.addTool(MCMatchObjP2MCRelator)
    mctruth.MCMatchObjP2MCRelator.RelTableLocations = rel_locs	
			      
b2dkpipi_d2KpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipi_d2KpipiTuple.TupleToolTrackIsolation.FillAsymmetry = True
b2dkpipi_d2KpipiTuple.TupleToolTrackIsolation.FillDeltaAngles = False
b2dkpipi_d2KpipiTuple.TupleToolTrackIsolation.MinConeAngle = 1.0

b2dkpipi_d2KpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2KpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2KpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2KpipiTuple.TupleToolGeometry.Verbose = True
b2dkpipi_d2KpipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipi_d2KpipiTuple.TupleToolGeometry.FillMultiPV = True

b2dkpipi_d2KpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipi_d2KpipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipi_d2KpipiTuple.UseLabXSyntax = True                          
b2dkpipi_d2KpipiTuple.RevertToPositiveID = False


b2dkpipi_d2KpipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]         
b2dkpipi_d2KpipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
b2dkpipi_d2KpipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+" ]  
b2dkpipi_d2KpipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.DTF.Verbose = True

b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]        
b2dkpipi_d2KpipiTuple.Bs.B0DTF.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.B0DTF.Substitutions = {'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
b2dkpipi_d2KpipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0" , "D_s+", "B~0"]  
b2dkpipi_d2KpipiTuple.Bs.B0DTF.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.B0DTF.Verbose = True

b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]        
b2dkpipi_d2KpipiTuple.Bs.BsDTF.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
b2dkpipi_d2KpipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
b2dkpipi_d2KpipiTuple.Bs.BsDTF.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.BsDTF.Verbose = True

b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]        
b2dkpipi_d2KpipiTuple.Bs.PV.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.PV.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.PV.Verbose = True


LoKiTool = b2dkpipi_d2KpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { "ETA" : "ETA" };

LoKiToolBs = b2dkpipi_d2KpipiTuple.Bs.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolBs")
LoKiToolBs.Variables = { "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()" };

b2dkpipi_d2KpipiTuple.addTool(TupleToolDecay, name="Ds")
LoKiToolDs = b2dkpipi_d2KpipiTuple.Ds.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolDs")
LoKiToolDs.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"};

b2dkpipi_d2KpipiTuple.addTool(TupleToolDecay, name="a_1_1260_plus")
LoKiToola_1_1260_plus = b2dkpipi_d2KpipiTuple.a_1_1260_plus.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToola_1_1260_plus")
LoKiToola_1_1260_plus.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)" };

#tagging config
from Configurables import BTaggingTool
tt_tagging = b2dkpipi_d2KpipiTuple.addTupleTool("TupleToolTagging") 
tt_tagging.Verbose = True
tt_tagging.AddMVAFeatureInfo = False
tt_tagging.AddTagPartsInfo = False 

btagtool = tt_tagging.addTool(BTaggingTool , name = "MyBTaggingTool")
from FlavourTagging.Tunings import applyTuning as applyFTTuning # pick the right tuning here ...
applyFTTuning(btagtool , tuning_version="Summer2017Optimisation_v4_Run2")
tt_tagging.TaggingToolName = btagtool.getFullName ()

#trigger config
b2dkpipi_d2Kpipitt = b2dkpipi_d2KpipiTuple.addTupleTool(TupleToolTISTOS)
b2dkpipi_d2Kpipitt.TriggerList = triggerlines
b2dkpipi_d2Kpipitt.FillL0 = True
b2dkpipi_d2Kpipitt.FillHlt1 = True
b2dkpipi_d2Kpipitt.FillHlt2 = True
b2dkpipi_d2Kpipitt.Verbose = True
b2dkpipi_d2Kpipitt.VerboseL0 = True
b2dkpipi_d2Kpipitt.VerboseHlt1 = True
b2dkpipi_d2Kpipitt.VerboseHlt2 = True

#main sequence
makeb2dkpipi_d2Kpipiseq = SelectionSequence("makeb2dkpipi_d2Kpipiseq", TopSelection = MyFilterSel)
b2dkpipi_d2KpipiTuple.Inputs = [makeb2dkpipi_d2Kpipiseq.outputLocation()]
b2dkpipi_d2Kpipiseq = GaudiSequencer("B2dkpipi_d2KpipiSeq")
#b2dkpipiseq.RootInTES = '/Event/{0}'.format(stream)
b2dkpipi_d2Kpipiseq.Members += [makeb2dkpipi_d2Kpipiseq.sequence(),b2dkpipi_d2KpipiTuple]


#
#
#
DaVinci().EventPreFilters = [stripFilter,checkPVs]
DaVinci().UserAlgorithms += [ b2dkpipiseq,b2dkpipi_d2pipipiseq,b2dkpipi_d2Kpipiseq]

DaVinci().DataType = year
if (data):
    DaVinci().Simulation = False
else:
    DaVinci().Simulation = True
DaVinci().EvtMax = -1
#DaVinci().EvtMax = 5000
DaVinci().SkipEvents = 0
DaVinci().PrintFreq = 50000
DaVinci().TupleFile = "b2dhhh.root"
DaVinci().Lumi = True

from Configurables import CondDB, CondDBAccessSvc

if (data):
	CondDB().LatestGlobalTagByDataType = year

    
else:  
    #if(year == "2012"):
	    #DaVinci().DDDBtag = "dddb-20130929-1"
    	    #if (down):
            	#DaVinci().CondDBtag = "sim-20130522-1-vc-md100"
    	    #else:
        	#DaVinci().CondDBtag = "sim-20130522-1-vc-mu100" 

    #if(year == "2011"):
	    #DaVinci().DDDBtag = "dddb-20130929"
    	    #if (down):
            	#DaVinci().CondDBtag = "sim-20130522-vc-md100"
    	    #else:
        	#DaVinci().CondDBtag = "sim-20130522-vc-mu100" 

    if(year == "2012"):
	    DaVinci().DDDBtag = "dddb-20170721-2"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20160321-2-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20160321-2-vc-mu100" 

    if(year == "2011"):
	    DaVinci().DDDBtag = "dddb-20170721-1"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20160614-1-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20160614-1-vc-mu100" 

    if(year == "2015"):
	    DaVinci().DDDBtag = "dddb-20170721-3"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20161124-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20161124-vc-mu100" 

    if(year == "2016"):
	    DaVinci().DDDBtag = "dddb-20170721-3"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20170721-2-vc-mu100" 

	
		
## Use the local input data
#from GaudiConf import IOHelper
#IOHelper().inputFiles([
    #'00069595_00000014_1.bhadroncompleteevent.dst'
#], clear=True)
