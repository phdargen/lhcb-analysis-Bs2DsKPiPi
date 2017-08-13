import argparse, sys, os, math

sys.path.append(os.environ["PIDPERFSCRIPTSROOT"] + "/scripts/python/PIDGenExpert/")
os.environ["ROOT_INCLUDE_PATH"] = os.pathsep + os.environ["MEERKATROOT"]

from ROOT import gROOT, TNtuple, TFile, TH1F, TH2F, TCanvas, TMath, gStyle, gSystem, RooRealVar

gSystem.Load("libMeerkatLib.so")

from ROOT import OneDimPhaseSpace, CombinedPhaseSpace, BinnedDensity, Logger

parser = argparse.ArgumentParser(description='PIDGen')
parser.add_argument('-i', '--input', type=str, default = None, 
                    help='Input file name')
parser.add_argument('-t', '--tree', type=str, default = "tree", 
                    help='Input tree name')
parser.add_argument('-o', '--output', type=str, default = "output.root", 
                    help='Output file name')
parser.add_argument('-p', '--pidvar', type=str, default = "PID_gen", 
                    help='Corrected PID variable')
parser.add_argument('-s', '--simpidvar', type=str, default = "PID", 
                    help='Simulated PID variable')
parser.add_argument('-m', '--ptvar', type=str, default = "Pt", 
                    help='Pt variable')
parser.add_argument('-q', '--pvar', type=str, default = "P", 
                    help='P variable')
parser.add_argument('-e', '--etavar', type=str, default = None, 
                    help='Eta variable')
parser.add_argument('-n', '--ntrvar', type=str, default = "nTracks", 
                    help='Ntracks variable')
parser.add_argument('-l', '--lowerpid', type=str, default = None, 
                    help='Lower PID value to generate')
parser.add_argument('-c', '--config', type=str, default = "p_V3ProbNNp", 
                    help='PID response to sample')
parser.add_argument('-d', '--dataset', type=str, default = "MagDown_2011", 
                    help='Dataset')
parser.add_argument('-v', '--var', type=str, default = "default", 
                    help='Variation (default, syst_N, stat_N etc.)')
parser.add_argument('-S', '--simversion', type=str, default = "sim08", 
                    help='Simulation version (sim08 or sim09 for run 1)')

parser.print_help()
args = parser.parse_args()

print args

infilename = args.input
intree = args.tree
outfilename = args.output
pidvar = args.pidvar
ptvar = args.ptvar
pvar = args.pvar
etavar = args.etavar
ntrvar = args.ntrvar
minpid = args.lowerpid
oldpidvar = args.simpidvar
conf = args.config
dataset = args.dataset
variant = args.var
simversion = args.simversion

import Run1.Config as Config
import Run1.ConfigMC as ConfigMCSim08
import Run1.ConfigMCSim09 as ConfigMCSim09

if not infilename : 
  print "Usage: PIDCorr.py [options]"
#  print "  For the usage example, look at pid_transform.sh file"
  print "  Available PID configs for sim08 are: "
  for i in sorted(Config.configs.keys()) : 
    if i in ConfigMCSim08.configs.keys() : 
      print "    ", i
  print "  Available PID configs for sim09 are: "
  for i in sorted(Config.configs.keys()) : 
    if i in ConfigMCSim09.configs.keys() : 
      print "    ", i
  quit()

if simversion == "sim08" : ConfigMC = ConfigMCSim08
elif simversion == "sim09" : ConfigMC = ConfigMCSim09
else : 
  print "Simulation version %s unknown" % simversion
  quit()

if variant == "default" : variant = "distrib"  # to do: change this name in CreatePIDPdf

datapdf = Config.eosrootdir + "/" + conf + "/" + dataset + "_" + variant + ".root"
simpdf = ConfigMC.eosrootdir + "/" + conf + "/" + dataset + "_" + variant + ".root"

transform_forward  = Config.configs[conf]["transform_forward"]
transform_backward = Config.configs[conf]["transform_backward"]

from math import sqrt

pidmin = 0.
pidmax = 1.
if 'limits' in Config.configs[conf] : 
  pidmin = Config.configs[conf]['limits'][0]
  pidmax = Config.configs[conf]['limits'][1]
if minpid == None : 
  minpid = pidmin
else :
  minpid = float(minpid)
  if minpid<pidmin : minpid = pidmin

# Calculate the minimum PID variable to generate (after transformation)
x = pidmin
pidmin = eval(transform_forward)
x = pidmax
pidmax = eval(transform_forward)
x = minpid
minpid = eval(transform_forward)

pid_phsp = OneDimPhaseSpace("PIDPhsp", pidmin , pidmax )
mom_phsp = OneDimPhaseSpace("MomPhsp", 5.5, 9.5)
eta_phsp = OneDimPhaseSpace("EtaPhsp", 1.5, 5.5)
ntr_phsp = OneDimPhaseSpace("NtrPhsp", 3.0, 6.5)
pidmom_phsp = CombinedPhaseSpace("PIDMomPhsp", pid_phsp, mom_phsp) 
pidmometa_phsp = CombinedPhaseSpace("PIDMomEtaPhsp", pidmom_phsp, eta_phsp) 
phsp = CombinedPhaseSpace("FullPhsp", pidmometa_phsp, ntr_phsp) 

infile = TFile.Open(infilename)
tree = infile.Get(intree) 
if not tree :  
  print "Ntuple not found!"
  quit()

nentries = tree.GetEntries()

datakde = BinnedDensity("KDEPDF", phsp, datapdf)
simkde = BinnedDensity("KDEPDF", phsp, simpdf)

gROOT.ProcessLine("struct MyStruct { \
  Double_t newpid; \
}")
from ROOT import MyStruct, AddressOf

s = MyStruct()

outfile = TFile.Open(outfilename, "RECREATE")
#if len(intree.split("/"))>1 : 
#  dirname = intree.split("/")[0]
#  outfile.mkdir(dirname)
#  outfile.cd(dirname)
newtree = tree.CloneTree(0)
newtree.Branch(pidvar, AddressOf(s, "newpid"), pidvar + "/D")
infile.cd() 

hdata = TH1F("hdata", "h", 100, minpid, pidmax) 
hsim = TH1F("hsim", "h", 100, minpid, pidmax) 

from ROOT import std, Double
from math import log

print transform_backward

var_code = compile("i.%s" % oldpidvar, '<string>', 'eval')
pid_code = compile(transform_backward, '<string>', 'eval')
oldpid_code = compile(transform_forward, '<string>', 'eval')
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

#  print point[0], point[1], point[2], point[3]

  hdata.Reset() 
  hsim.Reset() 
  datakde.slice(point, 0, hdata) 
  simkde.slice(point, 0, hsim) 

  x = eval(var_code)
  if transform_forward == "x" or x>=0 : 
    oldpid = eval(oldpid_code)
    if oldpid<pidmin or oldpid>pidmax : 
      x = oldpid
    else : 
      x = datakde.transform(hsim, hdata, oldpid)
    s.newpid = eval(pid_code)
  else : # The case for ProbNN<0, just leave as it is
    s.newpid = x

  newtree.Fill() 

  if (n % 1000 == 0) : 
    print "Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, OldPID=%f, PIDCorr=%f, X=%f" % \
       (n, nentries, point[1], point[2], point[3], \
        oldpid, s.newpid, x)

  n += 1

outfile.cd() 
newtree.Write() 
infile.Close() 
outfile.Close() 
