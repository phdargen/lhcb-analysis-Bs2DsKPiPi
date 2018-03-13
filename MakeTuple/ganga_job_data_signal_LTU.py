j5 = Job(name = "16LTU")
try:
	myApp = prepareGaudiExec('DaVinci','v41r4p4', myPath='.')
except:
	myApp = GaudiExec()
	myApp.directory = "./DaVinciDev_v41r4p4"
j5.application = myApp
j5.application.options = ["b2dkpipi_LTU_16.py"]
j5.backend=Dirac()
j5.application.platform = "x86_64-slc6-gcc49-opt"

datatmp9=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()
datatmp10=BKQuery('/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST', dqflag=['OK']).getDataset()

for f in datatmp10.files:
	datatmp9.append(f)

j5.inputdata = datatmp9
#j5.inputdata = [
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000006_1.bhadron.mdst',
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000010_1.bhadron.mdst',
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000013_1.bhadron.mdst',
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000031_1.bhadron.mdst' ]
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000038_1.bhadron.mdst',
#'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00069527/0000/00069527_00000042_1.bhadron.mdst',
j5.splitter = SplitByFiles( filesPerJob = 40)
j5.splitter.ignoremissing= True
j5.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j5.do_auto_resubmit = False
j5.parallel_submit = True
j5.submit()


