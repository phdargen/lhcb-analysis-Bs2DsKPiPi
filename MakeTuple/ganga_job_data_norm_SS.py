j2=Job(name='SS12N')
j2.application=DaVinci(version="v36r1p1")
j2.application.platform = "x86_64-slc6-gcc48-opt"
j2.application.optsfile = [ 'b2dpipipi_SS_12.py']
datatmp3=BKQuery('/LHCb/Collision12/Beam4000GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()
datatmp4=BKQuery('/LHCb/Collision12/Beam4000GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()

for f in datatmp4.files:
	datatmp3.append(f)

j2.backend=Dirac()
j2.inputdata = datatmp3
j2.splitter = SplitByFiles( filesPerJob = 40 )
j2.splitter.ignoremissing= True
j2.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j2.do_auto_resubmit = False
j2.parallel_submit = True
j2.submit()


#j5 = Job(name = "SS16N")
#myApp = GaudiExec()
#myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
#j5.application = myApp
#j5.application.options = ["b2dpipipi_SS_16.py"]
#j5.backend=Dirac()
#j5.application.platform = "x86_64-slc6-gcc49-opt"

#datatmp9=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()
#datatmp10=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()

#for f in datatmp10.files:
	#datatmp9.append(f)

#j5.inputdata = datatmp9
#j5.splitter = SplitByFiles( filesPerJob = 50)
#j5.splitter.ignoremissing= True
#j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j5.do_auto_resubmit = False
#j5.parallel_submit = True
#j5.submit()


