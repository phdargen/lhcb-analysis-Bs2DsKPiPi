#j=Job(name='11UMC')
#j.application=DaVinci(version="v36r1p1")
#j.application.optsfile = ['b2d3pi_11_MC_Up.py']
#j.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp=BKQuery('', dqflag=['OK']).getDataset()
#datatmp2=BKQuery('', dqflag=['OK']).getDataset()
#datatmp3=BKQuery('', dqflag=['OK']).getDataset()

#for f in datatmp2.files:
	#datatmp.append(f)

#for f in datatmp3.files:
	#datatmp.append(f)

#j.backend=Dirac()
#j.inputdata = datatmp
#j.splitter = SplitByFiles( filesPerJob = 1)
#j.splitter.ignoremissing= True
#j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.do_auto_resubmit = False
#j.parallel_submit = True
#j.submit()

#j=Job(name='11DMC')
#j.application=DaVinci(version="v36r1p1")
#j.application.optsfile = ['b2d3pi_11_MC_Down.py']
#j.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp=BKQuery('', dqflag=['OK']).getDataset()
#datatmp2=BKQuery('', dqflag=['OK']).getDataset()
#datatmp3=BKQuery('', dqflag=['OK']).getDataset()

#for f in datatmp2.files:
	#datatmp.append(f)

#for f in datatmp3.files:
	#datatmp.append(f)

#j.backend=Dirac()
#j.inputdata = datatmp
#j.splitter = SplitByFiles( filesPerJob = 1)
#j.splitter.ignoremissing= True
#j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.do_auto_resubmit = False
#j.parallel_submit = True
#j.submit()


j=Job(name='12UMC')
j.application=DaVinci(version="v36r1p1")
j.application.optsfile = ['b2d3pi_12_MC_Up.py']
j.application.platform = "x86_64-slc6-gcc48-opt"
datatmp=BKQuery('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp3=BKQuery('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

for f in datatmp3.files:
	datatmp.append(f)

j.backend=Dirac()
j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 1)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()


j=Job(name='12DMC')
j.application=DaVinci(version="v36r1p1")
j.application.optsfile = ['b2d3pi_12_MC_Down.py']
j.application.platform = "x86_64-slc6-gcc48-opt"
datatmp=BKQuery('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp3=BKQuery('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim09c/Trig0x409f0045/Reco14c/Stripping21Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

for f in datatmp3.files:
	datatmp.append(f)

j.backend=Dirac()
j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 1)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()



j = Job(name = "15UMC")
myApp = GaudiExec()
myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
j.application = myApp
j.application.options = ["b2d3pi_15_MC_Up.py"]
j.backend=Dirac()
j.application.platform = "x86_64-slc6-gcc49-opt"
datatmp=BKQuery('/MC/2015/Beam6500GeV-2015-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/MC/2015/Beam6500GeV-2015-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp3=BKQuery('/MC/2015/Beam6500GeV-2015-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

for f in datatmp3.files:
	datatmp.append(f)

j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 1)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()


j = Job(name = "15DMC")
myApp = GaudiExec()
myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
j.application = myApp
j.application.options = ["b2d3pi_15_MC_Down.py"]
j.backend=Dirac()
j.application.platform = "x86_64-slc6-gcc49-opt"
datatmp=BKQuery('/MC/2015/Beam6500GeV-2015-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/MC/2015/Beam6500GeV-2015-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp3=BKQuery('/MC/2015/Beam6500GeV-2015-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

for f in datatmp3.files:
	datatmp.append(f)

j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 1)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()


#j = Job(name = "16UMC")
#myApp = GaudiExec()
#myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
#j.application = myApp
#j.application.options = ["b2d3pi_16_MC_Up.py"]
#j.backend=Dirac()
#j.application.platform = "x86_64-slc6-gcc49-opt"
#datatmp=BKQuery('/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
#datatmp2=BKQuery('/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
#datatmp3=BKQuery('/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

#for f in datatmp2.files:
	#datatmp.append(f)

#for f in datatmp3.files:
	#datatmp.append(f)

#j.inputdata = datatmp
#j.splitter = SplitByFiles( filesPerJob = 1)
#j.splitter.ignoremissing= True
#j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.do_auto_resubmit = False
#j.parallel_submit = True
#j.submit()


j = Job(name = "16DMC")
myApp = GaudiExec()
myApp.directory = "/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/FT/DaVinciDev_v42r7p2"
j.application = myApp
j.application.options = ["b2d3pi_16_MC_Down.py"]
j.backend=Dirac()
j.application.platform = "x86_64-slc6-gcc49-opt"
datatmp=BKQuery('/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266068/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp2=BKQuery('/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266078/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()
datatmp3=BKQuery('/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1p1Filtered/13266088/B02D0HHH.STRIP.DST', dqflag=['OK']).getDataset()

for f in datatmp2.files:
	datatmp.append(f)

for f in datatmp3.files:
	datatmp.append(f)

j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 1)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()







## old MC
#j=Job(name='11DMC',application=DaVinci(version="v36r1"), backend=Dirac() )
#j.application.optsfile = [ 'b2d3pi_11_MC_Down.py' ]
#j.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp=BKQuery('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia8/Sim08i/Digi13/Trig0x40760037/Reco14c/Stripping21r1NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
#j.inputdata = datatmp
#j.do_auto_resubmit = True
#j.splitter = SplitByFiles( filesPerJob = 1 )
#j.splitter.ignoremissing= True
#j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.submit()

#j2=Job(name='11UMC',application=DaVinci(version="v36r1"), backend=Dirac() )
#j2.application.optsfile = [ 'b2d3pi_11_MC_Up.py' ]
#j2.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp2=BKQuery('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia8/Sim08i/Digi13/Trig0x40760037/Reco14c/Stripping21r1NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
#j2.inputdata = datatmp2
#j2.splitter = SplitByFiles( filesPerJob = 1 )
#j2.splitter.ignoremissing= True
#j2.do_auto_resubmit = True
#j2.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j2.submit() 

#j3=Job(name='12DMC',application=DaVinci(version="v36r1"), backend=Dirac() )
#j3.application.optsfile = [ 'b2d3pi_12_MC_Down.py' ]
#j3.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp3=BKQuery('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim08i/Digi13/Trig0x409f0045/Reco14c/Stripping21NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
#j3.inputdata = datatmp3
#j3.do_auto_resubmit = True
#j3.splitter = SplitByFiles( filesPerJob = 1 )
#j3.splitter.ignoremissing= True
#j3.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j3.submit() 

#j4=Job(name='12UMC',application=DaVinci(version="v36r1"), backend=Dirac() )
#j4.application.optsfile = [ 'b2d3pi_12_MC_Up.py' ]
#j4.application.platform = "x86_64-slc6-gcc48-opt"
#datatmp4=BKQuery('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim08i/Digi13/Trig0x409f0045/Reco14c/Stripping21NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
#j4.inputdata = datatmp4
#j4.splitter = SplitByFiles( filesPerJob = 1 )
#j4.splitter.ignoremissing= True
#j4.do_auto_resubmit = True
#j4.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j4.submit() 
