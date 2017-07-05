j=Job(name='11S')
#try:
	#myApp = prepareGaudiExec('DaVinci','v36r1', myPath='.')
#except:
	#myApp = GaudiExec()
	#myApp.directory = "./DaVinciDev_v36r1"
j.application=DaVinci(version="v36r1")
j.application.optsfile = ['b2dkpipi_11.py']
j.application.platform = "x86_64-slc6-gcc48-opt"
datatmp=BKQuery('/LHCb/Collision11/Beam3500GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21r1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/LHCb/Collision11/Beam3500GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21r1/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

j.backend=Dirac()
j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 10 )
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.do_auto_resubmit = True
j.submit()


j2=Job(name='12S')
#try:
	#myApp = prepareGaudiExec('DaVinci','v36r1', myPath='.')
#except:
	#myApp = GaudiExec()
	#myApp.directory = "./DaVinciDev_v36r1"
j.application=DaVinci(version="v36r1")
j.application.optsfile = ['b2dkpipi_11.py']
j2.application.platform = "x86_64-slc6-gcc48-opt"
j2.application.options = [ 'b2dkpipi_12.py']
datatmp3=BKQuery('/LHCb/Collision12/Beam4000GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()
datatmp4=BKQuery('/LHCb/Collision12/Beam4000GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21/90000000/BHADRONCOMPLETEEVENT.DST', dqflag=['OK']).getDataset()

for f in datatmp4.files:
	datatmp3.append(f)

j2.backend=Dirac()
j2.inputdata = datatmp3
j2.splitter = SplitByFiles( filesPerJob = 15 )
j2.splitter.ignoremissing= True
j2.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j2.do_auto_resubmit = True
j2.submit()



