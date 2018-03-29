# Collect list of LFNs (output stored in LFNs.txt):
# (on lxplus)

ganga -i collect-LFNs.py jobnumber

# faster way: in /gangadir/workspace/phdargen/LocalXML/jobid/
/bin/bash
for i in {0..1000}; do cat $i/output/__post* >> LFNs.txt ; done

#open in emacs
markieren control-leerzeichen : DiracFile:::*.root&&b2dhhh.root->
loeschen control-x rk

# Download files (modify outputLocation in collect-output.sh):

ssh -p28 sigma0
source collect-output.sh LFNs.txt

# ganga tricks 

start ganga without monitoring: 

ganga --no-mon
ganga --ganga-version 6.7.4 --no-mon

run Monitoring over selected jobs:

slice = j.subjobs.select(status='running')

runMonitoring(slice)

Run new job with inputdata from failed jobs:

ds = LHCbDataset()

for sj in j.subjobs.select(status='failed'):

             ds.extend(sj.inputdata)

j_new = j.copy()

j_new.inputdata = ds

j_new.splitter = SplitByFiles(filesPerJob = 10)

j_new.splitter.ignoremissing = True

if optionsfile has changed: 

j_new.unprepare()

j_new.submit()

paralles submit:

j.parallel_submit = True


# Get ddb tags used for MC production 

lhcb-proxy-init
lb-run LHCBDIRAC gaudirun.py get_bookkeeping_info.py 
lb-run LHCBDIRAC python get_bookkeeping_info.py

# Get generated lifetimes for given ddb tag

lb-run bender/latest python getGeneratedLifetimes.py


# Get dst file for testing

d = DiracFile('LFN:/lhcb/LHCb/Collision12/BHADRON.MDST/00041834/0000/00041834_00000001_1.bhadron.mdst')
d.localDir = '/work/dargent/Bs2DsKpipi/lhcb-analysis-Bs2DsKPiPi/MakeTuple/'
d.get()

lb-run DaVinci/v41r4p4 gaudirun.py b2dkpipi_16.py
