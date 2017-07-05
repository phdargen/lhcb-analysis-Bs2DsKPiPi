j3 = Job(name = "15S")
try:
	myApp = prepareGaudiExec('DaVinci','v38r1p4', myPath='.')
except:
	myApp = GaudiExec()
	myApp.directory = "./DaVinciDev_v38r1p4"
j3.application = myApp
j3.application.options = ["b2dkpipi_15.py"]
j3.backend=Dirac()
j3.application.platform = "x86_64-slc6-gcc49-opt"
datatmp5=BKQuery('/LHCb/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a/Stripping24r0p1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
datatmp6=BKQuery('/LHCb/Collision15/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco15a/Stripping24r0p1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

for f in datatmp6.files:
	datatmp5.append(f)

j3.backend=Dirac()
j3.inputdata = datatmp5
j3.splitter = SplitByFiles( filesPerJob = 10 )
j3.splitter.ignoremissing= True
j3.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j3.do_auto_resubmit = True
j3.submit()


j4 = Job(name = "16S")
try:
	myApp = prepareGaudiExec('DaVinci','v41r4p1', myPath='.')
except:
	myApp = GaudiExec()
	myApp.directory = "./DaVinciDev_v41r4p1"
j4.application = myApp
j4.application.options = ["b2dkpipi_16.py"]
j4.backend=Dirac()
j4.application.platform = "x86_64-slc6-gcc49-opt"
datatmp7=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
datatmp8=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

for f in datatmp8.files:
	datatmp7.append(f)

j4.inputdata = datatmp7
j4.splitter = SplitByFiles( filesPerJob = 20)
j4.splitter.ignoremissing= True
j4.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j4.do_auto_resubmit = True
j4.submit()


