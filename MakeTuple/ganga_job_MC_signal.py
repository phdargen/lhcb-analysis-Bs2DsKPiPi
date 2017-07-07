j=Job(name='11DMC',application=DaVinci(version="v36r1"), backend=Dirac() )
j.application.optsfile = [ 'b2dkpipi_11_MC_Down.py' ]
j.application.platform = "x86_64-slc6-gcc48-opt"
datatmp=BKQuery('/MC/2011/Beam3500GeV-2011-MagDown-Nu2-Pythia8/Sim08i/Digi13/Trig0x40760037/Reco14c/Stripping21r1NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
j.inputdata = datatmp
j.do_auto_resubmit = True
j.splitter = SplitByFiles( filesPerJob = 1 )
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.submit()

j2=Job(name='11UMC',application=DaVinci(version="v36r1"), backend=Dirac() )
j2.application.optsfile = [ 'b2dkpipi_11_MC_Up.py' ]
j2.application.platform = "x86_64-slc6-gcc48-opt"
datatmp2=BKQuery('/MC/2011/Beam3500GeV-2011-MagUp-Nu2-Pythia8/Sim08i/Digi13/Trig0x40760037/Reco14c/Stripping21r1NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
j2.inputdata = datatmp2
j2.splitter = SplitByFiles( filesPerJob = 1 )
j2.splitter.ignoremissing= True
j2.do_auto_resubmit = True
j2.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j2.submit() 

j3=Job(name='12DMC',application=DaVinci(version="v36r1"), backend=Dirac() )
j3.application.optsfile = [ 'b2dkpipi_12_MC_Down.py' ]
j3.application.platform = "x86_64-slc6-gcc48-opt"
datatmp3=BKQuery('/MC/2012/Beam4000GeV-2012-MagDown-Nu2.5-Pythia8/Sim08i/Digi13/Trig0x409f0045/Reco14c/Stripping21NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
j3.inputdata = datatmp3
j3.do_auto_resubmit = True
j3.splitter = SplitByFiles( filesPerJob = 1 )
j3.splitter.ignoremissing= True
j3.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j3.submit() 

j4=Job(name='12UMC',application=DaVinci(version="v36r1"), backend=Dirac() )
j4.application.optsfile = [ 'b2dkpipi_12_MC_Up.py' ]
j4.application.platform = "x86_64-slc6-gcc48-opt"
datatmp4=BKQuery('/MC/2012/Beam4000GeV-2012-MagUp-Nu2.5-Pythia8/Sim08i/Digi13/Trig0x409f0045/Reco14c/Stripping21NoPrescalingFlagged/13266007/ALLSTREAMS.DST', dqflag=['OK']).getDataset()
j4.inputdata = datatmp4
j4.splitter = SplitByFiles( filesPerJob = 1 )
j4.splitter.ignoremissing= True
j4.do_auto_resubmit = True
j4.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j4.submit() 
