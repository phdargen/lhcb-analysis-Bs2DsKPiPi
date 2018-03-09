import argparse, sys, os, math

sys.path.append(os.environ["PIDPERFSCRIPTSROOT"] + "/scripts/python/PIDGenExpert/")
os.environ["ROOT_INCLUDE_PATH"] = os.pathsep + os.environ["MEERKATROOT"]

from ROOT import gROOT, TNtuple, TFile, TH1F, TH2F, TCanvas, TMath, TRandom3, gStyle, gSystem, RooRealVar

gSystem.Load("libMeerkatLib.so")

from ROOT import OneDimPhaseSpace, CombinedPhaseSpace, BinnedDensity, Logger

import Run1.Config as ConfigRun1
import Run2.Config as ConfigRun2

parser = argparse.ArgumentParser(description='PIDGen')
parser.add_argument('-i', '--input',   type=str, default = None, 
                    help='Input file name')
parser.add_argument('-t', '--tree',    type=str, default = "tree", 
                    help='Input tree name')
parser.add_argument('-o', '--output',  type=str, default = "output.root", 
                    help='Output file name')
parser.add_argument('-p', '--pidvar',  type=str, default = "PID_gen", 
                    help='PID variable')
parser.add_argument('-m', '--ptvar',   type=str, default = "Pt", 
                    help='Pt variable')
parser.add_argument('-q', '--pvar',    type=str, default = "P", 
                    help='P variable')
parser.add_argument('-e', '--etavar',  type=str, default = None, 
                    help='Eta variable (if None, calculated from P and Pt)')
parser.add_argument('-n', '--ntrvar',  type=str, default = "nTracks", 
                    help='Ntracks variable')
parser.add_argument('-l', '--lowerpid', type=str, default = None, 
                    help='Lower PID value to generate')
parser.add_argument('-z', '--higherpid', type=str, default = None, 
                    help='Higher PID value to generate')
parser.add_argument('-c', '--config',  type=str, default = "p_V3ProbNNp", 
                    help='PID response to sample')
parser.add_argument('-d', '--dataset', type=str, default = "MagDown_2011", 
                    help='Dataset (polarity_year)')
parser.add_argument('-v', '--var',     type=str, default = "default", 
                    help='Variation (default, syst_N, stat_N etc.)')
parser.add_argument('-s', '--seed',   type=str, default = None, 
                    help='Initial random seed')

parser.print_help()
args = parser.parse_args()

print args

infilename = args.input
intree = args.tree
outfilename = args.output
pidvar = args.pidvar
ptvar = args.ptvar
pvar  = args.pvar
etavar = args.etavar
ntrvar = args.ntrvar
minpid = args.lowerpid
maxpid = args.higherpid
config = args.config
dataset = args.dataset
variant = args.var
seed = args.seed

if not infilename : 
  print "Usage: PIDGen.py [options]"
#  print "  For the usage example, look at pid_resample.sh file"
  print "  Available PID configs are: "
  print "    For Run1 : "
  for i in sorted(ConfigRun1.configs.keys()) : 
    print "      ", i
  print "    For Run2 : "
  for i in sorted(ConfigRun2.configs.keys()) : 
    print "      ", i
  quit()

if variant == "default" : variant = "distrib"  # to do: change this name in CreatePIDPdf

year = None
run = None
try : 
  year = dataset.split("_")[1]
except : 
  print 'Dataset format "%s" not recognized. Should be {MagUp,MagDown}_[Year]' % dataset
  quit()
if year in ["2011", "2012"] : 
  run = 1
elif year in ["2015", "2016"] : 
  run = 2
else : 
  print 'Data taking year "%s" not recognized'
  quit()

print year, run, dataset

if run == 1 : 
  calibfilename = ConfigRun1.eosrootdir + "/" + config + "/" + "%s_%s.root" % (dataset, variant)
  transform_forward  = ConfigRun1.configs[config]['transform_forward']
  transform_backward = ConfigRun1.configs[config]['transform_backward']
  configs = ConfigRun1.configs
else : 
  calibfilename = ConfigRun2.eosrootdir + "/" + config + "/" + "%s_%s.root" % (dataset, variant)
  gamma = ConfigRun2.configs[config]['gamma']
  if gamma<0 : 
    transform_forward = "(1.-(1.-x)**%f)" % abs(gamma)
    transform_backward = "(1.-(1.-x)**%f)" % (1./abs(gamma))
  else : 
    transform_forward = "((x)**%f)" % abs(gamma)
    transform_backward = "((x)**%f)" % (1./abs(gamma))
  configs = ConfigRun2.configs

from math import sqrt

pidmin = 0.
pidmax = 1.
if 'limits' in configs[config] : 
  pidmin = configs[config]['limits'][0]
  pidmax = configs[config]['limits'][1]
if minpid == None : 
  minpid = pidmin
else :
  minpid = float(minpid)
  if minpid<pidmin : minpid = pidmin
if maxpid == None : 
  maxpid = pidmax
else :
  maxpid = float(maxpid)
  if maxpid>pidmax : maxpid = pidmax


# Calculate the minimum PID variable to generate (after transformation)
x = pidmin
pidmin = eval(transform_forward)
x = pidmax
pidmax = eval(transform_forward)
x = minpid
minpid = eval(transform_forward)

pid_phsp = OneDimPhaseSpace("PIDPhsp", minpid , maxpid )
mom_phsp = OneDimPhaseSpace("MomPhsp", 4.6, 10.5)
eta_phsp = OneDimPhaseSpace("EtaPhsp", 1.5, 5.5)
ntr_phsp = OneDimPhaseSpace("NtrPhsp", 3.0, 7.0)
pidmom_phsp = CombinedPhaseSpace("PIDMomPhsp", pid_phsp, mom_phsp) 
pidmometa_phsp = CombinedPhaseSpace("PIDMomEtaPhsp", pidmom_phsp, eta_phsp) 
phsp = CombinedPhaseSpace("FullPhsp", pidmometa_phsp, ntr_phsp) 

infile = TFile.Open(infilename)
tree = infile.Get(intree) 
if not tree :  
  print "Ntuple not found!"
  quit()

kde = BinnedDensity("KDEPDF", phsp, calibfilename)
nentries = tree.GetEntries()

gROOT.ProcessLine("struct MyStruct { \
  Double_t newpid; \
}")
from ROOT import MyStruct, AddressOf

s = MyStruct()

outfile = TFile.Open(outfilename, "RECREATE")
newtree = tree.CloneTree(0)
newtree.Branch(pidvar, AddressOf(s, "newpid"), pidvar + "/D")
infile.cd() 

h = TH1F("h", "h", 100, minpid, pidmax) 
rnd = TRandom3()
if seed != None : 
  rnd.SetSeed(int(seed))

from ROOT import std, Double
from math import log

#x_code = compile("i.%s" % pidvar, '<string>', 'eval')
pid_code = compile(transform_backward, '<string>', 'eval')
log_pt_code = compile("log(i.%s)" % ptvar, '<string>', 'eval')
if etavar == None : 
  p_code = compile("i.%s" % pvar, '<string>', 'eval')
  pt_code = compile("i.%s" % ptvar, '<string>', 'eval')
else : 
  eta_code = compile("i.%s" % etavar, '<string>', 'eval')
log_ntracks_code = compile("log(float(i.%s))" % ntrvar, '<string>', 'eval')

Logger.setLogLevel(1)

n = 0
for i in tree : 
  point = std.vector(Double)(4) 
  point[0] = (pidmin + pidmax)/2.
  point[1] = eval(log_pt_code)
  point[3] = eval(log_ntracks_code) 
  if etavar == None : 
    point[2] = -math.log( math.tan( math.asin( eval(pt_code)/eval(p_code) )/2. ) )
  else : 
    point[2] = eval(eta_code)

  h.Reset() 
  kde.slice(point, 0, h) 

  if h.Integral() > 0 : 
    x = h.GetRandom()
  else : 
    x = minpid + (pidmax-minpid)*rnd.Rndm()

  s.newpid = eval(pid_code)

  newtree.Fill() 

  if (n % 1000 == 0) : 
    print "Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, PIDGen=%f" % \
       (n, nentries, point[1], point[2], point[3], \
        s.newpid)

  n += 1

outfile.cd() 
newtree.Write() 
infile.Close() 
outfile.Close() 
