from Gaudi.Configuration import *

from Configurables import DaVinci

from Configurables import FilterDesktop

from Configurables import CombineParticles

############# Global settings
year = "2011"
data = False
down = True
stream = "AllStreams"
if (data):
    stream = "BhadronCompleteEvent"

from Configurables import CheckPV
checkPVs = CheckPV("checkPVs")
checkPVs.MinPVs = 1
checkPVs.MaxPVs = -1

triggerlines = [ "L0HadronDecision", 
                "L0MuonDecision",
                "L0GlobalDecision",
                #"L0PhysicsDecision",
                "Hlt1TrackAllL0Decision", 
                #'Hlt1TrackAllL0TightDecision',
                #"Hlt1L0AnyDecision",
                #"Hlt1MBNoBiasDecision",
                #"Hlt1GlobalDecision",
                "Hlt1TrackMVADecision",
                "Hlt1TwoTrackMVADecision",
                "Hlt1TrackMVALooseDecision",
                "Hlt1TwoTrackMVALooseDecision",
                 #hlt2
                 #topo
                #'Hlt2CharmHadD2HHHHDecision', 
                #'Hlt2CharmHadD2HHHHWideMassDecision', 
                'Hlt2IncPhiDecision', 
                "Hlt2PhiIncPhiDecision",
                # 'Hlt2ExpressDs2PhiPiDecision', \
                # 'Hlt2Topo2BodySimpleDecision', \
                # 'Hlt2Topo3BodySimpleDecision', \
                # 'Hlt2Topo4BodySimpleDecision', \
                 'Hlt2Topo2BodyBBDTDecision', \
                 'Hlt2Topo3BodyBBDTDecision', \
                 'Hlt2Topo4BodyBBDTDecision', \
                # 'Hlt2TopoMu2BodyBBDTDecision', \
                # 'Hlt2TopoMu3BodyBBDTDecision', \
                # 'Hlt2TopoMu4BodyBBDTDecision', \
                # 'Hlt2TopoE2BodyBBDTDecision', \
                # 'Hlt2TopoE3BodyBBDTDecision', \
                # 'Hlt2TopoE4BodyBBDTDecision', \
                #'Hlt2Topo2BodyCombosDecision', \
                # 'Hlt2Topo3BodyCombosDecision', \
                # 'Hlt2Topo4BodyCombosDecision', 
                 'Hlt2Topo2BodyDecision', \
                 'Hlt2Topo3BodyDecision', \
                 'Hlt2Topo4BodyDecision'		 
                # 'Hlt2RadiativeTopoTrackTOSDecision', \
                # 'Hlt2RadiativeTopoPhotonL0Decision'
 ] 

from Configurables import (BTagging)   
Btag = BTagging("BTagging")
Btag.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHFULLDSTBeauty2CharmLine/Particles"]
Btag.OutputLevel    = 6

from FlavourTagging.Tunings import TuneTool
def configureTaggingTools(BTuple, BsBd):
	
	    from Configurables import BTagging, BTaggingTool
	    toolname = BsBd+"TaggingTool"
	    BTuple.TaggingToolName ="BTaggingTool/"+toolname
	    BTuple.ExtraName = toolname
	    if(data):
		    TuneTool(BTuple,"Stripping21",toolname)
	    else :
		    TuneTool(BTuple,"Stripping21_MC",toolname)
	    BTuple.Verbose = True
	    BTuple.addTool(BTaggingTool,name=toolname)
	    BTuple.__getattr__(toolname).ForceSignalID = BsBd


###############   Pre Filter, does not really do much except choose only candidates passing the Stripping line, maybe beneficial to performance
from Configurables import LoKi__HDRFilter as StripFilter
stripFilter = StripFilter( 'stripPassFilter',\
                           Code = "HLT_PASS('StrippingB02DKPiPiD2HHHPIDBeauty2CharmLineDecision')",\
                           Location= "/Event/Strip/Phys/DecReports")

############# DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import TupleToolTISTOS
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence
from Configurables import FilterDesktop
from Configurables import PrintDecayTree, PrintDecayTreeTool
##subpid stuff
#from Configurables import SubPIDMMFilter
from Configurables import SubstitutePID,BTaggingTool
from Configurables import TupleToolDecayTreeFitter, TupleToolTrackIsolation, TupleToolTagging, TupleToolRecoStats, TupleToolKinematic, TupleToolGeometry 
from Configurables import LoKi__Hybrid__TupleTool

reqsel = AutomaticData(Location = "/Event/"+stream+"/Phys/StrippingB02DKPiPiD2HHHPIDBeauty2CharmLine_Line/Particles")

from Configurables import FilterDesktop

#B0 -> (D- -> K K pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2KKpi_Tuple")
b2dkpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^K+ ^K- ^pi-) ^(K_1(1270)+ -> ^K+ ^pi+ ^pi-)]CC"
b2dkpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> K+ K- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC)" ,
"K_1_1270_plus" : "[[B0]cc -> (D- -> K+ K- pi-) ^(K_1(1270)+ -> K+ pi+ pi-) ]CC",
"K_plus" : "[[B0]cc -> (D- -> K+ K- pi-) (K_1(1270)+ -> ^K+ pi+ pi-)  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> K+ K- pi-) (K_1(1270)+ -> K+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> K+ K- pi-) (K_1(1270)+ -> K+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> K+ K- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"K_plus_fromDs" : "[[B0]cc -> (D- -> ^K+ K- pi-) (K_1(1270)+ -> K+ pi+ pi-)  ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> K+ ^K- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> K+ K- ^pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC"
}
b2dkpipiTuple.ReFitPVs = True

#config tools
b2dkpipiTuple.ToolList +=  ["TupleToolGeometry", \
                            "TupleToolKinematic", \
                            "TupleToolPrimaries", \
                            "TupleToolEventInfo", \
                            "TupleToolTrackInfo", \
                            "TupleToolRecoStats", \
                            "TupleToolAngles", \
                            "TupleToolPid", \
                            "TupleToolPhotonInfo", \
                            "TupleToolTrackIsolation",
                            "TupleToolTagging" ]

if (data==False):
    b2dkpipiTuple.ToolList +=  [
                            "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo"]
    MCTruth = TupleToolMCTruth()
    MCTruth.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
        ]
    b2dkpipiTuple.addTool(MCTruth) 

b2dkpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipiTuple.TupleToolGeometry.Verbose = False

b2dkpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipiTuple.UseLabXSyntax = True                          
b2dkpipiTuple.RevertToPositiveID = False

LoKiTool = b2dkpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                       };
                       
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

#tt_geometry = b2dkpipiTuple.addTupleTool('TupleToolGeometry')
#tt_geometry.Verbose = True
#tt_geometry.RefitPVs = True
#tt_geometry.FillMultiPV = True

#tagging config
b2dkpipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
b2dkpipiTuple.TupleToolTagging.Verbose = True
b2dkpipiTuple.TupleToolTagging.StoreTaggersInfo = True

tag=b2dkpipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
configureTaggingTools(tag, "Bs")

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
b2dkpipiTuple.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
#bu2kpipimumuseq.Members += [makebu2kpipimumuseq.sequence(), bu2kpipimumuprinter, bu2kpipimumutuple]
b2dkpipiseq.Members += [b2dkpipiTuple]


#
#B0 -> (D- -> pi pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2pipipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2pipipi_Tuple")
b2dkpipi_d2pipipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^pi-) ^(K_1(1270)+ -> ^K+ ^pi+ ^pi-)]CC"
b2dkpipi_d2pipipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC)" ,
"K_1_1270_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-) ^(K_1(1270)+ -> K+ pi+ pi-) ]CC",
"K_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-) (K_1(1270)+ -> ^K+ pi+ pi-)  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-) (K_1(1270)+ -> K+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- pi-) (K_1(1270)+ -> K+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- pi-) (K_1(1270)+ -> K+ pi+ pi-)  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"pi_minus2_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^pi-) (K_1(1270)+ -> K+ pi+ pi-) ]CC"
}
b2dkpipi_d2pipipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2pipipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              "TupleToolRecoStats", \
                              "TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolPhotonInfo", \
                              "TupleToolTrackIsolation",
                              "TupleToolTagging" ]

if(data==False):
    b2dkpipi_d2pipipiTuple.ToolList +=  [
                            "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo"]
    MCTruth_d2pipipi = TupleToolMCTruth()
    MCTruth_d2pipipi.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
        ]
    b2dkpipi_d2pipipiTuple.addTool(MCTruth_d2pipipi) 

b2dkpipi_d2pipipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipi_d2pipipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipi_d2pipipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2pipipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2pipipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2pipipiTuple.TupleToolGeometry.Verbose = False

b2dkpipi_d2pipipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipi_d2pipipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipi_d2pipipiTuple.UseLabXSyntax = True                          
b2dkpipi_d2pipipiTuple.RevertToPositiveID = False

LoKiTool_d2pipipi = b2dkpipi_d2pipipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool_d2pipipi.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
        };
                       
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

#tagging config
b2dkpipi_d2pipipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
b2dkpipi_d2pipipiTuple.TupleToolTagging.Verbose = True
b2dkpipi_d2pipipiTuple.TupleToolTagging.StoreTaggersInfo = True

tag_d2pipipi=b2dkpipi_d2pipipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
configureTaggingTools(tag_d2pipipi, "Bs")

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
b2dkpipi_d2pipipiTuple.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
b2dkpipi_d2pipipiseq = GaudiSequencer("B2dkpipi_d2pipipiSeq")
b2dkpipi_d2pipipiseq.Members += [b2dkpipi_d2pipipiTuple]


#
#B0 -> (D- -> K pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2KpipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2Kpipi_Tuple")
b2dkpipi_d2KpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^K-) ^(K_1(1270)+ -> ^K+ ^pi+ ^pi-)]CC"
b2dkpipi_d2KpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- K-) (K_1(1270)+ -> K+ pi+ pi-) ]CC)" ,
"K_1_1270_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) ^(K_1(1270)+ -> K+ pi+ pi-) ]CC",
"K_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) (K_1(1270)+ -> ^K+ pi+ pi-)  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) (K_1(1270)+ -> K+ ^pi+ pi-)  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- K-) (K_1(1270)+ -> K+ pi+ ^pi-) ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- K-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- K-) (K_1(1270)+ -> K+ pi+ pi-)  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- K-) (K_1(1270)+ -> K+ pi+ pi-) ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^K-) (K_1(1270)+ -> K+ pi+ pi-) ]CC"
}
b2dkpipi_d2KpipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2KpipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              "TupleToolRecoStats", \
                              "TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolPhotonInfo", \
                              "TupleToolTrackIsolation",
                              "TupleToolTagging" ]

if(data==False):
    b2dkpipi_d2KpipiTuple.ToolList +=  [
                            "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo"]
    MCTruth_d2Kpipi = TupleToolMCTruth()
    MCTruth_d2Kpipi.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
        ]
    b2dkpipi_d2KpipiTuple.addTool(MCTruth_d2Kpipi) 
			      
b2dkpipi_d2KpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipi_d2KpipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipi_d2KpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2KpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2KpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2KpipiTuple.TupleToolGeometry.Verbose = False

b2dkpipi_d2KpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipi_d2KpipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipi_d2KpipiTuple.UseLabXSyntax = True                          
b2dkpipi_d2KpipiTuple.RevertToPositiveID = False

LoKiTool_d2KpipiTuple = b2dkpipi_d2KpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool_d2KpipiTuple.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
        };

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

#tagging config
b2dkpipi_d2KpipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
b2dkpipi_d2KpipiTuple.TupleToolTagging.Verbose = True
b2dkpipi_d2KpipiTuple.TupleToolTagging.StoreTaggersInfo = True

tag_d2Kpipi=b2dkpipi_d2KpipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
configureTaggingTools(tag_d2Kpipi, "Bs")

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
b2dkpipi_d2KpipiTuple.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
b2dkpipi_d2Kpipiseq = GaudiSequencer("B2dkpipi_d2KpipiSeq")
b2dkpipi_d2Kpipiseq.Members += [b2dkpipi_d2KpipiTuple]


#
#
#
DaVinci().EventPreFilters = [stripFilter,checkPVs]
DaVinci().UserAlgorithms += [ b2dkpipiseq]

DaVinci().DataType = year
if (data):
    DaVinci().Simulation = False
else:
    DaVinci().Simulation = True
DaVinci().EvtMax = -1
#DaVinci().EvtMax = 1000
DaVinci().SkipEvents = 0
DaVinci().PrintFreq = 1000
DaVinci().TupleFile = "b2dhhh.root"
DaVinci().Lumi = True

if (data):
    if(year == "2016"):
	    DaVinci().CondDBtag = "cond-20161004"
    	    DaVinci().DDDBtag = "dddb-20150724"
    
    if(year == "2015"):
	    DaVinci().CondDBtag = "cond-20150828"
    	    DaVinci().DDDBtag = "dddb-20150724"
	    
    if(year == "2012"):
	    DaVinci().CondDBtag = "cond-20141107"
    	    DaVinci().DDDBtag = "dddb-20130929-1"

    if(year == "2011"):
	    DaVinci().CondDBtag = "cond-20141107"
    	    DaVinci().DDDBtag = "dddb-20130929"
    
else:  
    if(year == "2012"):
	    DaVinci().DDDBtag = "dddb-20130929-1"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20130522-1-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20130522-1-vc-mu100" 

    if(year == "2011"):
	    DaVinci().DDDBtag = "dddb-20130929"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20130522-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20130522-vc-mu100" 
