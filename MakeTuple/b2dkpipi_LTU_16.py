from Gaudi.Configuration import *
MessageSvc().Format = "% F%80W%S%7W%R%T %0W%M"

from Configurables import DaVinci
from Configurables import FilterDesktop
from Configurables import CombineParticles
from Configurables import EventNodeKiller
eventNodeKiller = EventNodeKiller('DAQkiller')
eventNodeKiller.Nodes = ['DAQ',
                         'pRec']

############# Global settings
year = "2016"
data = True
down = False
stream = "AllStreams"
if (data):
    stream = "Bhadron"


line = 'B02DsKPiPiLTUBD2HHHBeauty2CharmLine'

from Configurables import CheckPV
checkPVs = CheckPV("checkPVs")
checkPVs.MinPVs = 1
checkPVs.MaxPVs = -1

triggerlines = [ "L0HadronDecision", 
                "L0MuonDecision",
                #"L0GlobalDecision",
                #"L0PhysicsDecision",
                "Hlt1TrackAllL0Decision", 
                #'Hlt1TrackAllL0TightDecision',
                #"Hlt1L0AnyDecision",
                #"Hlt1MBNoBiasDecision",
                #"Hlt1GlobalDecision",
                "Hlt1TrackMVADecision",
                "Hlt1TwoTrackMVADecision",
                #"Hlt1TrackMVALooseDecision",
                #"Hlt1TwoTrackMVALooseDecision",
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

#from Configurables import (BTagging)   
#Btag = BTagging("BTagging")
#Btag.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHFULLDSTBeauty2CharmLine/Particles"]
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
                           Code = "HLT_PASS('StrippingB02DsKPiPiLTUBD2HHHBeauty2CharmLineDecision')",\
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
from Configurables import TupleToolDecayTreeFitter, TupleToolTrackIsolation, TupleToolRecoStats, TupleToolKinematic, TupleToolGeometry 
from Configurables import LoKi__Hybrid__TupleTool

#reqsel = AutomaticData(Location = "/Event/"+stream+"/Phys/StrippingB02DKPiPiD2HHHPIDBeauty2CharmLine_Line/Particles")

from Configurables import FilterDesktop

#B0 -> (D- -> K K pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2KKpi_Tuple")
b2dkpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^K+ ^K- ^pi-) ^K+ ^pi+ ^pi-]CC"
b2dkpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> K+ K- pi-) K+ pi+ pi- ]CC)" ,
"K_plus" : "[[B0]cc -> (D- -> K+ K- pi-) ^K+ pi+ pi-  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> K+ K- pi-)  K+ ^pi+ pi-  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> K+ K- pi-)  K+ pi+ ^pi- ]CC",
"Ds" : "[[B0]cc -> ^(D- -> K+ K- pi-)  K+ pi+ pi- ]CC",
"K_plus_fromDs" : "[[B0]cc -> (D- -> ^K+ K- pi-)  K+ pi+ pi-  ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> K+ ^K- pi-) K+ pi+ pi- ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> K+ K- ^pi-)  K+ pi+ pi- ]CC"
}
b2dkpipiTuple.ReFitPVs = True

#config tools
b2dkpipiTuple.ToolList +=  [ #"TupleToolGeometry", \
                            "TupleToolKinematic", \
                            "TupleToolPrimaries", \
                            "TupleToolEventInfo", \
                            "TupleToolTrackInfo", \
                            #"TupleToolRecoStats", \
                            "TupleToolAngles", \
                            "TupleToolPid", \
                            "TupleToolPhotonInfo", \
                            #"TupleToolTrackIsolation",
                            #"TupleToolTagging" 
			    ]

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

#b2dkpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
#b2dkpipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipiTuple.TupleToolGeometry.Verbose = True
b2dkpipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipiTuple.TupleToolGeometry.FillMultiPV = True

#b2dkpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
#b2dkpipiTuple.TupleToolRecoStats.Verbose = True
#b2dkpipiTuple.UseLabXSyntax = True                          
#b2dkpipiTuple.RevertToPositiveID = False

LoKiTool = b2dkpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                       };
                       
b2dkpipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]        
b2dkpipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ X+ X- ^(Charm & X-)' : 'D_s-', 'Beauty -> X- X- X+ ^(Charm & X+)' : 'D_s+'  }
b2dkpipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+", ]  
b2dkpipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipiTuple.Bs.DTF.Verbose = True

#b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
#b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]
#b2dkpipiTuple.Bs.B0DTF.constrainToOriginVertex = True
#b2dkpipiTuple.Bs.B0DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
#b2dkpipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0", "D_s+", "B~0" ]  
#b2dkpipiTuple.Bs.B0DTF.UpdateDaughters = True
#b2dkpipiTuple.Bs.B0DTF.Verbose = True

#b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
#b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]         
#b2dkpipiTuple.Bs.BsDTF.constrainToOriginVertex = True
#b2dkpipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
#b2dkpipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
#b2dkpipiTuple.Bs.BsDTF.UpdateDaughters = True
#b2dkpipiTuple.Bs.BsDTF.Verbose = True

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
#b2dkpipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
#b2dkpipiTuple.TupleToolTagging.Verbose = True
#b2dkpipiTuple.TupleToolTagging.StoreTaggersInfo = True

#tag=b2dkpipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
#configureTaggingTools(tag, "Bs")

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

from Configurables import LoKi__Hybrid__EvtTupleTool
LoKi_EvtTuple=LoKi__Hybrid__EvtTupleTool("LoKi_EvtTuple")
LoKi_EvtTuple.VOID_Variables = {
    # track information
    "nLong"       : "RECSUMMARY(LHCb.RecSummary.nLongTracks      , -999, '', False )"
    ,"nUpstream"   : "RECSUMMARY(LHCb.RecSummary.nUpstreamTracks  , -999, '', False )"
    ,"nDownstream" : "RECSUMMARY(LHCb.RecSummary.nDownstreamTracks, -999, '', False )"
    ,"nBackward"   : "RECSUMMARY(LHCb.RecSummary.nBackTracks      , -999, '', False )" 
    ,"nMuon"       : "RECSUMMARY(LHCb.RecSummary.nMuonTracks      , -999, '', False )"
    ,"nVELO"       : "RECSUMMARY(LHCb.RecSummary.nVeloTracks      , -999, '', False )"
    ,"nTracks"         : "RECSUMMARY( LHCb.RecSummary.nTracks,-1,'/Event/Rec/Summary',False )"
     # pileup
    ,"nPVs"        : "RECSUMMARY(LHCb.RecSummary.nPVs, -999, '', False )"
     
     # tracking multiplicities
    ,"nSpdDigits"  : "RECSUMMARY(LHCb.RecSummary.nSPDhits,    -999, '', False )"
    ,"nITClusters" : "RECSUMMARY(LHCb.RecSummary.nITClusters, -999, '', False )"
    ,"nTTClusters" : "RECSUMMARY(LHCb.RecSummary.nTTClusters, -999, '', False )"
    }

LoKi_EvtTuple.Preambulo +=['from LoKiTracks.decorators import *',
                           'from LoKiNumbers.decorators import *',
                           'from LoKiCore.functions  import *' ]
    
#b2dkpipiTuple.Bs.addTool(TupleToolTISTOS,name="TisTosB")
#b2dkpipiTuple.Bs.TisTosB.Verbose=True
#b2dkpipiTuple.Bs.TisTosB.TriggerList= triggers
#printer
#b2dkpipiprinter = PrintDecayTree("PrintB2Dkpipi")
#b2dkpipiprinter.addTool( PrintDecayTreeTool, name = "PrintDecay" )
#b2dkpipiprinter.PrintDecay.Information = "Name M P Px Py Pz Pt chi2"
#b2dkpipiprinter.Inputs = [ 'Phys/{0}/Particles'.format(line)]

# Put the GECs in the ntuples !!
b2dkpipiTuple.ToolList+=["LoKi::Hybrid::EvtTupleTool/LoKi_EvtTuple"]
b2dkpipiTuple.addTool(LoKi_EvtTuple)

b2dkpipiTuple.Bs.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipiTuple.Bs.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipiTuple.Bs.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}


b2dkpipiTuple.addTool(TupleToolDecay, name="Ds")
b2dkpipiTuple.Ds.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipiTuple.Ds.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipiTuple.Ds.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}


#main sequence
b2dkpipiTuple.Inputs = ['Phys/{0}/Particles'.format(line)]
#[ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
b2dkpipiseq.RootInTES = '/Event/{0}'.format(stream)
#bu2kpipimumuseq.Members += [makebu2kpipimumuseq.sequence(), bu2kpipimumuprinter, bu2kpipimumutuple]
b2dkpipiseq.Members += [b2dkpipiTuple]





#
#B0 -> (D- -> pi pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2pipipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2pipipi_Tuple")
b2dkpipi_d2pipipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^pi-) ^K+ ^pi+ ^pi-]CC"
b2dkpipi_d2pipipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- pi-) K+ pi+ pi- ]CC)" ,
"K_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-) ^K+ pi+ pi-  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> pi+ pi- pi-)  K+ ^pi+ pi-  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- pi-) K+ pi+ ^pi- ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- pi-) K+ pi+ pi- ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- pi-) K+ pi+ pi-  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- pi-) K+ pi+ pi- ]CC",
"pi_minus2_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^pi-) K+ pi+ pi- ]CC"
}
b2dkpipi_d2pipipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2pipipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              #"TupleToolRecoStats", \
                              "TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolPhotonInfo", \
                              #"TupleToolTrackIsolation",
                              #"TupleToolTagging" 
			      ]

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

#b2dkpipi_d2pipipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
#b2dkpipi_d2pipipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipi_d2pipipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2pipipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2pipipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2pipipiTuple.TupleToolGeometry.Verbose = True
b2dkpipi_d2pipipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipi_d2pipipiTuple.TupleToolGeometry.FillMultiPV = True

#b2dkpipi_d2pipipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
#b2dkpipi_d2pipipiTuple.TupleToolRecoStats.Verbose = True
#b2dkpipi_d2pipipiTuple.UseLabXSyntax = True                          
#b2dkpipi_d2pipipiTuple.RevertToPositiveID = False

LoKiTool_d2pipipi = b2dkpipi_d2pipipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool_d2pipipi.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
        };
                       
b2dkpipi_d2pipipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]      
b2dkpipi_d2pipipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ X+ X- ^(Charm & X-)' : 'D_s-', 'Beauty -> X- X- X+ ^(Charm & X+)' : 'D_s+'  }
b2dkpipi_d2pipipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+" ]  
b2dkpipi_d2pipipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.DTF.Verbose = True

#b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
#b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]        
#b2dkpipi_d2pipipiTuple.Bs.B0DTF.constrainToOriginVertex = True
#b2dkpipi_d2pipipiTuple.Bs.B0DTF.Substitutions = {'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
#b2dkpipi_d2pipipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0","D_s+", "B~0"]  
#b2dkpipi_d2pipipiTuple.Bs.B0DTF.UpdateDaughters = True
#b2dkpipi_d2pipipiTuple.Bs.B0DTF.Verbose = True

#b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
#b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]        
#b2dkpipi_d2pipipiTuple.Bs.BsDTF.constrainToOriginVertex = True
#b2dkpipi_d2pipipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
#b2dkpipi_d2pipipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
#b2dkpipi_d2pipipiTuple.Bs.BsDTF.UpdateDaughters = True
#b2dkpipi_d2pipipiTuple.Bs.BsDTF.Verbose = True

b2dkpipi_d2pipipiTuple.Bs.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]        
b2dkpipi_d2pipipiTuple.Bs.PV.constrainToOriginVertex = True
b2dkpipi_d2pipipiTuple.Bs.PV.UpdateDaughters = True
b2dkpipi_d2pipipiTuple.Bs.PV.Verbose = True

#tagging config
#b2dkpipi_d2pipipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
#b2dkpipi_d2pipipiTuple.TupleToolTagging.Verbose = True
#b2dkpipi_d2pipipiTuple.TupleToolTagging.StoreTaggersInfo = True

#tag_d2pipipi=b2dkpipi_d2pipipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
#configureTaggingTools(tag_d2pipipi, "Bs")

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

# Put the GECs in the ntuples !!
b2dkpipi_d2pipipiTuple.ToolList+=["LoKi::Hybrid::EvtTupleTool/LoKi_EvtTuple"]
b2dkpipi_d2pipipiTuple.addTool(LoKi_EvtTuple)

b2dkpipi_d2pipipiTuple.Bs.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipi_d2pipipiTuple.Bs.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipi_d2pipipiTuple.Bs.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}

b2dkpipi_d2pipipiTuple.addTool(TupleToolDecay, name="Ds")
b2dkpipi_d2pipipiTuple.Ds.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipi_d2pipipiTuple.Ds.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipi_d2pipipiTuple.Ds.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}


#main sequence
b2dkpipi_d2pipipiTuple.Inputs = ['Phys/{0}/Particles'.format(line)]
b2dkpipiseq_d2pipipi = GaudiSequencer("B2dkpipiSeq_d2pipipi")
b2dkpipiseq_d2pipipi.RootInTES = '/Event/{0}'.format(stream)
b2dkpipiseq_d2pipipi.Members += [b2dkpipi_d2pipipiTuple]

#
#
#B0 -> (D- -> K pi pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipi_d2KpipiTuple = DecayTreeTuple("Bs2DsKpipi_Ds2Kpipi_Tuple")
b2dkpipi_d2KpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^pi+ ^pi- ^K-) ^K+ ^pi+ ^pi-]CC"
b2dkpipi_d2KpipiTuple.Branches= {
"Bs" : "^([[B0]cc -> (D- -> pi+ pi- K-) K+ pi+ pi- ]CC)" ,
"K_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) ^K+ pi+ pi-  ]CC",
"pi_plus" : "[[B0]cc -> (D- -> pi+ pi- K-) K+ ^pi+ pi-  ]CC",
"pi_minus" : "[[B0]cc -> (D- -> pi+ pi- K-) K+ pi+ ^pi- ]CC",
"Ds" : "[[B0]cc -> ^(D- -> pi+ pi- K-) K+ pi+ pi- ]CC",
"pi_plus_fromDs" : "[[B0]cc -> (D- -> ^pi+ pi- K-)  K+ pi+ pi-  ]CC",
"pi_minus_fromDs" : "[[B0]cc -> (D- -> pi+ ^pi- K-)  K+ pi+ pi- ]CC",
"K_minus_fromDs" : "[[B0]cc -> (D- -> pi+ pi- ^K-)  K+ pi+ pi- ]CC"
}
b2dkpipi_d2KpipiTuple.ReFitPVs = True

#config tools
b2dkpipi_d2KpipiTuple.ToolList +=  ["TupleToolGeometry", \
                              "TupleToolKinematic", \
                              "TupleToolPrimaries", \
                              "TupleToolEventInfo", \
                              "TupleToolTrackInfo", \
                              #"TupleToolRecoStats", \
                              "TupleToolAngles", \
                              "TupleToolPid", \
                              "TupleToolPhotonInfo", \
                              #"TupleToolTrackIsolation",
                             # "TupleToolTagging" 
			     ]

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
			      
#b2dkpipi_d2KpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
#b2dkpipi_d2KpipiTuple.TupleToolTrackIsolation.Verbose = True

b2dkpipi_d2KpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipi_d2KpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipi_d2KpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipi_d2KpipiTuple.TupleToolGeometry.Verbose = True
b2dkpipi_d2KpipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipi_d2KpipiTuple.TupleToolGeometry.FillMultiPV = True

#b2dkpipi_d2KpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
#b2dkpipi_d2KpipiTuple.TupleToolRecoStats.Verbose = True
#b2dkpipi_d2KpipiTuple.UseLabXSyntax = True                          
#b2dkpipi_d2KpipiTuple.RevertToPositiveID = False

LoKiTool_d2KpipiTuple = b2dkpipi_d2KpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool_d2KpipiTuple.Variables = { #"inMuon" : "switch(INMUON, 1, 0)",
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
        };

b2dkpipi_d2KpipiTuple.addTool(TupleToolDecay, name="Bs")
b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]         
b2dkpipi_d2KpipiTuple.Bs.DTF.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ X+ X- ^(Charm & X-)' : 'D_s-', 'Beauty -> X- X- X+ ^(Charm & X+)' : 'D_s+'  }
b2dkpipi_d2KpipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+" ]  
b2dkpipi_d2KpipiTuple.Bs.DTF.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.DTF.Verbose = True

#b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
#b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]        
#b2dkpipi_d2KpipiTuple.Bs.B0DTF.constrainToOriginVertex = True
#b2dkpipi_d2KpipiTuple.Bs.B0DTF.Substitutions = {'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
#b2dkpipi_d2KpipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0" , "D_s+", "B~0"]  
#b2dkpipi_d2KpipiTuple.Bs.B0DTF.UpdateDaughters = True
#b2dkpipi_d2KpipiTuple.Bs.B0DTF.Verbose = True

#b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("BsDTF"))
#b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/BsDTF" ]        
#b2dkpipi_d2KpipiTuple.Bs.BsDTF.constrainToOriginVertex = True
#b2dkpipi_d2KpipiTuple.Bs.BsDTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' , 'B0 -> Hadron Charm' : 'B_s0' , 'B~0 -> Hadron Charm' : 'B_s~0'  }
#b2dkpipi_d2KpipiTuple.Bs.BsDTF.daughtersToConstrain = ["D_s-", "B_s0", "D_s+", "B_s~0" ]  
#b2dkpipi_d2KpipiTuple.Bs.BsDTF.UpdateDaughters = True
#b2dkpipi_d2KpipiTuple.Bs.BsDTF.Verbose = True

b2dkpipi_d2KpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]        
b2dkpipi_d2KpipiTuple.Bs.PV.constrainToOriginVertex = True
b2dkpipi_d2KpipiTuple.Bs.PV.UpdateDaughters = True
b2dkpipi_d2KpipiTuple.Bs.PV.Verbose = True

#tagging config
#b2dkpipi_d2KpipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
#b2dkpipi_d2KpipiTuple.TupleToolTagging.Verbose = True
#b2dkpipi_d2KpipiTuple.TupleToolTagging.StoreTaggersInfo = True

#tag_d2Kpipi=b2dkpipi_d2KpipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
#configureTaggingTools(tag_d2Kpipi, "Bs")

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

# Put the GECs in the ntuples !!
b2dkpipi_d2KpipiTuple.ToolList+=["LoKi::Hybrid::EvtTupleTool/LoKi_EvtTuple"]
b2dkpipi_d2KpipiTuple.addTool(LoKi_EvtTuple)

b2dkpipi_d2KpipiTuple.Bs.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipi_d2KpipiTuple.Bs.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipi_d2KpipiTuple.Bs.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}

b2dkpipi_d2KpipiTuple.addTool(TupleToolDecay, name="Ds")
b2dkpipi_d2KpipiTuple.Ds.addTool(LoKi__Hybrid__TupleTool('LoKi_Cone'))
b2dkpipi_d2KpipiTuple.Ds.ToolList +=  ["LoKi::Hybrid::TupleTool/LoKi_Cone"]
b2dkpipi_d2KpipiTuple.Ds.LoKi_Cone.Variables = {
  "CONEANGLE"      : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEANGLE',-1.)",
  "CONEMULT"       : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEMULT', -1.)",
  "ptasy_1.00"     : "RELINFO('/Event/Bhadron/Phys/B02DsKPiPiLTUBD2HHHBeauty2CharmLine/P2ConeVar3','CONEPTASYM',-1.)",
}

#main sequence
b2dkpipi_d2KpipiTuple.Inputs = ['Phys/{0}/Particles'.format(line)]
b2dkpipiseq_d2Kpipi = GaudiSequencer("B2dkpipiSeq_d2Kpipi")
b2dkpipiseq_d2Kpipi.RootInTES = '/Event/{0}'.format(stream)
b2dkpipiseq_d2Kpipi.Members += [b2dkpipi_d2KpipiTuple]



#
#
DaVinci().EventPreFilters = [stripFilter]
DaVinci().UserAlgorithms += [ eventNodeKiller, b2dkpipiseq,b2dkpipiseq_d2pipipi,b2dkpipiseq_d2Kpipi]


DaVinci().DataType = year
if (data):
    DaVinci().Simulation = False
else:
    DaVinci().Simulation = True
DaVinci().EvtMax = -1
#DaVinci().EvtMax = 1000
DaVinci().SkipEvents = 0
DaVinci().PrintFreq = 10000
DaVinci().TupleFile = "b2dhhh.root"
DaVinci().Lumi = True
DaVinci().InputType = "MDST"

from Configurables import CondDB, CondDBAccessSvc

if (data):
	CondDB().LatestGlobalTagByDataType = year
    
else:  
    if(year == "2012"):
	    DaVinci().DDDBtag = "dddb-20130929-1"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20141210-1-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20141210-1-vc-mu100" 

    if(year == "2011"):
	    DaVinci().DDDBtag = "dddb-20130929"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20141210-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20141210-vc-mu100" 


## Use the local input data
#from GaudiConf import IOHelper
#IOHelper().inputFiles([
    #'./00069527_00000006_1.bhadron.mdst'
#], clear=True)