j5 = Job(name = "16LTU")
myApp = GaudiExec()
myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
j5.application = myApp
j5.application.options = ["b2dkpipi_LTU_16.py"]
j5.backend=Dirac()
j5.application.platform = "x86_64-slc6-gcc49-opt"

datatmp9=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()
datatmp10=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()

for f in datatmp10.files:
	datatmp9.append(f)

j5.inputdata = datatmp9
j5.splitter = SplitByFiles( filesPerJob = 40)
j5.splitter.ignoremissing= True
j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j5.do_auto_resubmit = False
j5.parallel_submit = True
j5.submit()



j5 = Job(name = "17LTU")
#try:
	#myApp = prepareGaudiExec('DaVinci','v42r7p2', myPath='.')
#except:
	#myApp = GaudiExec()
	#myApp.directory = "./DaVinciDev_v42r7p2"
myApp = GaudiExec()
myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
j5.application = myApp
j5.application.options = ["b2dkpipi_LTU_17.py"]
j5.backend=Dirac()
j5.application.platform = "x86_64-slc6-gcc49-opt"

datatmp9=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()
datatmp10=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()

for f in datatmp10.files:
	datatmp9.append(f)

j5.inputdata = datatmp9
j5.splitter = SplitByFiles( filesPerJob = 50)
j5.splitter.ignoremissing= True
j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j5.do_auto_resubmit = False
j5.parallel_submit = True
j5.submit()