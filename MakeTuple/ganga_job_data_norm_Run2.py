#j3 = Job(name = "15N")
##try:
	##myApp = prepareGaudiExec('DaVinci','v38r1p4', myPath='.')
##except:
#myApp = GaudiExec()
#myApp.directory = "./FT_new3/DaVinciDev_v42r7p2"
#j3.application = myApp
#j3.application.options = ["b2d3pi_15.py"]
#j3.backend=Dirac()
#j5.application.platform = "x86_64-slc6-gcc49-opt"
#datatmp5=BKQuery('/LHCb/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a/Stripping24r1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
#datatmp6=BKQuery('/LHCb/Collision15/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco15a/Stripping24r1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

#for f in datatmp6.files:
	#datatmp5.append(f)

#j3.backend=Dirac()
#j3.inputdata = datatmp5
#j3.splitter = SplitByFiles( filesPerJob = 30 )
#j3.splitter.ignoremissing= True
#j3.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j3.do_auto_resubmit = False
#j3.parallel_submit = True
#j3.submit()

#j4 = Job(name = "16N")
##try:
	##myApp = prepareGaudiExec('DaVinci','v41r4p1', myPath='.')
##except:
#myApp = GaudiExec()
#myApp.directory = "./FT_new3/DaVinciDev_v42r7p2"
#j4.application = myApp
#j4.application.options = ["b2d3pi_16.py"]
#j4.backend=Dirac()
#j5.application.platform = "x86_64-slc6-gcc49-opt"
#datatmp7=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r1p1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
#datatmp8=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r1p1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

#for f in datatmp8.files:
	#datatmp7.append(f)

#j4.inputdata = datatmp7
#j4.splitter = SplitByFiles( filesPerJob = 50)
#j4.splitter.ignoremissing= True
#j4.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j4.do_auto_resubmit = False
#j4.parallel_submit = True
#j4.submit()


#j5 = Job(name = "17N")
#myApp = GaudiExec()
#myApp.directory = "./FT_new3/DaVinciDev_v42r7p2"
#j5.application = myApp
#j5.application.options = ["b2d3pi_17.py"]
#j5.backend=Dirac()
#j5.application.platform = "x86_64-slc6-gcc49-opt"
#datatmp9=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
#datatmp10=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

#for f in datatmp10.files:
	#datatmp9.append(f)

#j5.inputdata = datatmp9
#j5.splitter = SplitByFiles( filesPerJob = 50)
#j5.splitter.ignoremissing= True
#j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j5.do_auto_resubmit = False
#j5.parallel_submit = True
#j5.submit()














j5 = Job(name = "18N")
myApp = GaudiExec()
myApp.directory = "./FT_new3/DaVinciDev_v42r7p2"
j5.application = myApp
j5.application.options = ["b2dpipipi_18.py"]
j5.backend=Dirac()
j5.application.platform = "x86_64-slc6-gcc49-opt"
datatmp9=BKQuery('/LHCb/Collision18/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco18/Stripping34/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
datatmp10=BKQuery('/LHCb/Collision18/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco18/Stripping34/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

for f in datatmp10.files:
	datatmp9.append(f)

j5.inputdata = datatmp9
j5.splitter = SplitByFiles( filesPerJob = 100)
j5.splitter.ignoremissing= True
j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j5.do_auto_resubmit = False
j5.parallel_submit = True
j5.submit()
