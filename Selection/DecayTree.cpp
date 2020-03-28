#define DecayTree_cxx
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>

using namespace std; 

TTree* DecayTree::GetInputTree(){
    
    TString tupleName;
    if(_decay==Decay::signal)tupleName="Bs2DsKpipi_";
    else tupleName = "Bs2Dspipipi_";
    if(_Ds_finalState==Ds_finalState::pipipi)tupleName+="Ds2pipipi_Tuple/DecayTree";
    else if(_Ds_finalState==Ds_finalState::Kpipi)tupleName+="Ds2Kpipi_Tuple/DecayTree";
    else tupleName+="Ds2KKpi_Tuple/DecayTree";    

    TChain* chain = new TChain(tupleName);

    if(_decay==Decay::signal && _data==DataType::data && _year == 11){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/11/";
        chain->Add(loc+"b2dhhh*.root");
    }
    
    else if(_decay==Decay::signal && _data==DataType::data && _year == 12 && _ss == false){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/12/";
        chain->Add(loc+"b2dhhh*.root");
	loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/12b/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 12 && _ss == true){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/12SS/";
        chain->Add(loc+"b2dhhh*.root");
	loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/12SSb/";
        chain->Add(loc+"b2dhhh*.root");
    }    

    else if(_decay==Decay::norm && _data==DataType::data && _year == 11){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Norm/Data/11/";
        chain->Add(loc+"b2dhhh*.root");
	loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Norm/Data/11b/";
        chain->Add(loc+"b2dhhh*.root");
    }
    
    else if(_decay==Decay::norm && _data==DataType::data && _year == 12){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Norm/Data/12/";
        chain->Add(loc+"b2dhhh*.root");
	loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Norm/Data/12b/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 15){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_10_18/Signal/Data/15/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 16 && _ltu == false && _ss == false){
        TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_10_18/Signal/Data/16/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 16 && _ltu == true){
	TString loc = "/auto/data/kecke/BsDsKpipi/LTU_16/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 16 && _ss == true){
	TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/16SS/";
        chain->Add(loc+"b2dhhh*.root");
	loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Signal/Data/16SSb/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 17 && _ltu == false){
        TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_10_18/Signal/Data/17/";
        chain->Add(loc+"b2dhhh*.root");
    }

    else if(_decay==Decay::signal && _data==DataType::data && _year == 17 && _ltu == true){
        TString loc = "/auto/data/kecke/BsDsKpipi/LTU_17/";
        chain->Add(loc+"b2dhhh*.root");
    }

    
    else if(_decay==Decay::signal && _data==DataType::data && _year == 18 && _ltu == false){
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154864/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154875/b2dhhh.root");
        chain->Add("root://hake7.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353154/353154879/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154886/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154890/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154897/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154910/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154914/b2dhhh.root");
        chain->Add("root://lobster4.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353154/353154916/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154921/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154925/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154928/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154932/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154935/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154939/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154943/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154945/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154948/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154953/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154956/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154962/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154965/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154968/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154974/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154977/b2dhhh.root");
        chain->Add("root://mouse6.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353154/353154981/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154987/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154990/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154993/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353154/353154997/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155005/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155008/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155011/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155015/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155018/b2dhhh.root");
        chain->Add("root://hake1.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155021/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155027/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155033/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155042/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155045/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155048/b2dhhh.root");
        chain->Add("root://hake2.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155055/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155062/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155066/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155076/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155081/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155083/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155087/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155090/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155093/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155099/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155102/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155105/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155108/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155111/b2dhhh.root");
        chain->Add("root://guppy9.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155114/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155117/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155120/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155123/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155125/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155133/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155134/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155137/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155139/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155141/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155143/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155146/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155149/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155152/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155156/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155159/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155162/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155164/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155170/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155172/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155175/b2dhhh.root");
        chain->Add("root://hake16.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155178/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155181/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155184/b2dhhh.root");
        chain->Add("root://shark11.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155190/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155193/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155198/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155201/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155204/b2dhhh.root");
        chain->Add("root://shark10.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155230/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155239/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155252/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155264/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155282/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155303/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155324/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155343/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155367/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155386/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155404/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155421/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155446/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155475/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155480/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155483/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155487/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155499/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155504/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155514/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155517/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155525/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155529/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155538/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155540/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155544/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155548/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155550/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155553/b2dhhh.root");
        chain->Add("root://guppy11.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155556/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155560/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155563/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155566/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155569/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155574/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155578/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155581/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155584/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155587/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155589/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155592/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155595/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155598/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155601/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155604/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155607/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155613/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155616/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155619/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155622/b2dhhh.root");
        chain->Add("root://hake12.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155625/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155628/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155631/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155633/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155638/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155643/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155646/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155650/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155652/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155656/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155662/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155665/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155669/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155676/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155683/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155686/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155688/b2dhhh.root");
        chain->Add("root://shark9.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155691/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155692/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155694/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155696/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155698/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155702/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155705/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155707/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155712/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155714/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155716/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155720/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155722/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155725/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155727/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155729/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155731/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155743/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155744/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155746/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155747/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155748/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155749/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155750/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155752/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155757/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155759/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155761/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155763/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155767/b2dhhh.root");
        chain->Add("root://lobster6.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155769/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155771/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155773/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155775/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155777/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155781/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155783/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155785/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155787/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155789/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155791/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155794/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155796/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155798/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155801/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155805/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155807/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155809/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155811/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155814/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155815/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155817/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155819/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155821/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155823/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155825/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155827/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155829/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155834/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155836/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155838/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155841/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155843/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155845/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155848/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155849/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155853/b2dhhh.root");
        chain->Add("root://guppy13.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353155/353155855/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155857/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155859/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155861/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155863/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155865/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155870/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155877/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155879/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155882/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155884/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155887/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155890/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155893/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155899/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155901/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155903/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155907/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155911/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155912/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155921/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155923/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155925/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155930/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155932/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155934/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353155/353155977/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156012/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156063/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156091/b2dhhh.root");
        chain->Add("root://mouse7.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353156/353156265/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156325/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156460/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156544/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353156/353156971/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157032/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157059/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157093/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157141/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157164/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157178/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157179/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157180/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157181/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157183/b2dhhh.root");
        chain->Add("root://shark1.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353157/353157185/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157186/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157187/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157188/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157189/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157190/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157191/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157192/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157193/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157194/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157195/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157196/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157197/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157198/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157199/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157200/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157201/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157207/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157215/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157228/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157236/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157252/b2dhhh.root");
        chain->Add("root://hake16.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353157/353157284/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157293/b2dhhh.root");
        chain->Add("root://hake16.grid.surfsara.nl:1094/pnfs/grid.sara.nl/data/lhcb/user/p/phdargen/2020_03/353157/353157325/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157339/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157352/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157366/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157375/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157384/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157389/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157400/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157408/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157415/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157430/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157432/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157436/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157440/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157442/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157443/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157444/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157445/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157446/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157447/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157448/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157449/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157450/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157451/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157453/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157454/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157455/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157457/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157460/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157463/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157464/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157465/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157466/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157591/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157620/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157669/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157706/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157755/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157792/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157804/b2dhhh.root");
        chain->Add("root://eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/p/phdargen/2020_03/353157/353157833/b2dhhh.root");
    }
    
    else if(_decay==Decay::signal && _data==DataType::data && _year == 18 && _ltu == true){
    }
    
    
    else if(_decay==Decay::norm && _data==DataType::data && _year == 15){
        TString loc = "/auto/data/kecke/B2DPiPiPi/15/";
        chain->Add(loc+"b2dhhh*.root");
   }

    else if(_decay==Decay::norm && _data==DataType::data && _year == 16){
        TString loc = "/auto/data/kecke/B2DPiPiPi/16/";
        chain->Add(loc+"b2dhhh*.root");
   }

    else if(_decay==Decay::norm && _data==DataType::data && _year == 17){
        TString loc = "/auto/data/dargent/BsDsKpipi/Stripped_04_18/Norm/Data/17/";
        chain->Add(loc+"b2dhhh*.root");
   }
    
   else if(_decay==Decay::norm && _data==DataType::data && _year == 18){
    }

  else if(_decay==Decay::norm && _data==DataType::mc && _year == 11){
        TString loc = "/auto/data/kecke/B2DPiPiPi/11DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DPiPiPi/11UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::norm && _data==DataType::mc && _year == 12){
        TString loc = "/auto/data/kecke/B2DPiPiPi/12DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DPiPiPi/12UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::norm && _data==DataType::mc && _year == 15){
        TString loc = "/auto/data/kecke/B2DPiPiPi/15DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DPiPiPi/15UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::norm && _data==DataType::mc && _year == 16){
        TString loc = "/auto/data/kecke/B2DPiPiPi/16DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DPiPiPi/16UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }

  else if(_decay==Decay::signal && _data==DataType::mc && _year == 11){
        TString loc = "/auto/data/kecke/B2DKPiPi/11DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DKPiPi/11UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::signal && _data==DataType::mc && _year == 12){
        TString loc = "/auto/data/kecke/B2DKPiPi/12DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DKPiPi/12UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::signal && _data==DataType::mc && _year == 15){
        TString loc = "/auto/data/kecke/B2DKPiPi/15DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DKPiPi/15UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }
  else if(_decay==Decay::signal && _data==DataType::mc && _year == 16){
        TString loc = "/auto/data/kecke/B2DKPiPi/16DMC/";
        chain->Add(loc+"b2dhhh*.root");
        loc = "/auto/data/kecke/B2DKPiPi/16UMC/";
        chain->Add(loc+"b2dhhh*.root");
   }




    else {
        TString fileName = _inFileLoc + "Stripped/";

        if(_decay==Decay::signal)fileName+="Signal/";
        else fileName += "Norm/";

        if(_data==DataType::data) {
            fileName+= "Data/";
            fileName+= _year; 
        }
        else {
            fileName+= "MC/";
            fileName+= _year; 
            fileName+= "U";
        }
        fileName+= "/b2dhhh_*.root"; 

	if(_bkg)fileName = "/auto/data/dargent/old_Bs2DsKpipi/MC/Norm/Bkg/Dst3pi-U/b2dhhh_*.root" ; 

        cout << "Using the files: " << endl;
        cout << fileName << endl << endl;
        chain->Add(fileName);
        if(_data!=DataType::data) chain->Add(fileName.ReplaceAll(TString("U/b2"),TString("D/b2")));
    }

    if(chain->GetEntries()==0){
        cout << "ERROR: No events found" << endl;
        throw "ERROR";
    }
    
    return (TTree*)chain;
}

DecayTree::DecayTree(Decay::Type decay, Year::Type year, Ds_finalState::Type finalState, DataType::Type dataType, TString polarity, TString inFileLoc, TString outFileLoc, Bool_t bkg, Bool_t ltu, Bool_t ss ) : 
fChain(0), _decay(decay), _year(year), _Ds_finalState(finalState), _data(dataType), _polarity(polarity), _inFileLoc(inFileLoc), _outFileLoc(outFileLoc), _bkg(bkg), _ltu(ltu), _ss(ss)
{    
    cout << "Requested to process files with options: " << endl << endl;

    TString s1,s2,s3,s4;
    
    if(_data==DataType::data)s1="Data";
    else s1 = "MC";
    if(_decay==Decay::signal)s2="signal";
    else s2 = "norm";
    if(_Ds_finalState==Ds_finalState::pipipi)s3="Ds2pipipi";
    else if(_Ds_finalState==Ds_finalState::Kpipi)s3="Ds2Kpipi";
    else s3="Ds2KKpi";
    stringstream ss_year;
    ss_year << _year;
    string str_year = ss_year.str();

    cout << "DataType: " << s1 << endl;
    cout << "Decay: " << s2 << endl;
    cout << "Ds finalState: " << s3 << endl;
    cout << "Year: " << _year << endl;
    cout << "Polarity: " << _polarity  << endl;
    cout << "Bkg: " << _bkg << endl;
    cout << "LTU: " << _ltu << endl;
    cout << "SS: " << _ss << endl << endl;

    _outFileName = _outFileLoc;    
    _outFileName += "Mini/";
    _outFileName += s1;  
    _outFileName += "/";  
    _outFileName += s2;   
    _outFileName += "_";  
    _outFileName += s3; 
    _outFileName += "_";     
    _outFileName += str_year;
    if(_polarity == "Up") _outFileName += "_up";
    if(_polarity == "Down") _outFileName += "_down";
    if(_bkg)_outFileName += "_Dstar_bkg";
    if(_ltu)_outFileName += "_LTU";
    if(_ss)_outFileName += "_SS";	
    _outFileName += ".root";    
}

inline Bool_t DecayTree::TriggerCuts(Long64_t i){

    b_Bs_L0Global_TIS->GetEntry(i);
    b_Bs_L0HadronDecision_TOS->GetEntry(i);
    if( (!Bs_L0Global_TIS) && (!Bs_L0HadronDecision_TOS)) return false;
    
    if(_year == 15 || _year == 16 || _year == 17 || _year == 18 ){
        b_Bs_Hlt1TrackMVADecision_TOS->GetEntry(i);
        b_Bs_Hlt1TwoTrackMVADecision_TOS->GetEntry(i);
        if((!Bs_Hlt1TrackMVADecision_TOS) && (!Bs_Hlt1TwoTrackMVADecision_TOS) ) return false;
        
        b_Bs_Hlt2Topo2BodyDecision_TOS->GetEntry(i);
        b_Bs_Hlt2Topo3BodyDecision_TOS->GetEntry(i);
        b_Bs_Hlt2Topo4BodyDecision_TOS->GetEntry(i);
	if(_year == 15){
		b_Bs_Hlt2IncPhiDecision_TOS->GetEntry(i);
        	if((!Bs_Hlt2Topo2BodyDecision_TOS) &&  (!Bs_Hlt2Topo3BodyDecision_TOS) && (!Bs_Hlt2Topo4BodyDecision_TOS) && (!Bs_Hlt2IncPhiDecision_TOS) ) return false;	
	}
	else {
		b_Bs_Hlt2PhiIncPhiDecision_TOS->GetEntry(i);
        	if((!Bs_Hlt2Topo2BodyDecision_TOS) &&  (!Bs_Hlt2Topo3BodyDecision_TOS) && (!Bs_Hlt2Topo4BodyDecision_TOS) && (!Bs_Hlt2PhiIncPhiDecision_TOS) ) return false;
	}
    }
    
    else if(_year == 11 || _year == 12){
        b_Bs_Hlt1TrackAllL0Decision_TOS->GetEntry(i);        
        if(!Bs_Hlt1TrackAllL0Decision_TOS) return false;
        
        b_Bs_Hlt2Topo2BodyBBDTDecision_TOS->GetEntry(i);
        b_Bs_Hlt2Topo3BodyBBDTDecision_TOS->GetEntry(i);
        b_Bs_Hlt2Topo4BodyBBDTDecision_TOS->GetEntry(i);
        b_Bs_Hlt2IncPhiDecision_TOS->GetEntry(i);
        if((!Bs_Hlt2Topo2BodyBBDTDecision_TOS) &&  (!Bs_Hlt2Topo3BodyBBDTDecision_TOS) && (!Bs_Hlt2Topo4BodyBBDTDecision_TOS) 
        && (!Bs_Hlt2IncPhiDecision_TOS)) return false;
    }
    
    return true;
}

inline Bool_t DecayTree::LooseCuts(Long64_t i){

    b_Bs_DIRA_OWNPV->GetEntry(i);
    if(Bs_DIRA_OWNPV<0.99994) return false;
    
    b_Bs_IPCHI2_OWNPV->GetEntry(i);
    if(Bs_IPCHI2_OWNPV>20) return false;
    
    b_Bs_FDCHI2_OWNPV->GetEntry(i);
    if(Bs_FDCHI2_OWNPV<100) return false;
    
    b_Bs_ENDVERTEX_CHI2->GetEntry(i);
    b_Bs_ENDVERTEX_NDOF->GetEntry(i);
    if((Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)> 8) return false;
    
//     b_Bs_TAU->GetEntry(i);
//     if(Bs_TAU < 0.0002) return false;

    b_Bs_MM->GetEntry(i);
    if(Bs_MM < 4800. || Bs_MM > 6000.) return false;
    
    b_Bs_DTF_M->GetEntry(i);
    if(Bs_DTF_M[0] < 4800. || Bs_DTF_M[0] > 6000.) return false;
    
    b_Bs_BsDTF_M->GetEntry(i);
    if(Bs_BsDTF_M[0] < 4800. || Bs_BsDTF_M[0] > 6000.) return false;
    
    if(!_bkg)b_Bs_PV_M->GetEntry(i);
    if(!_bkg)if (Bs_PV_M[0] < 4800. || Bs_PV_M[0] > 6000.) return false;
    
    b_Ds_ENDVERTEX_Z->GetEntry(i);
    b_Bs_ENDVERTEX_Z->GetEntry(i);
    if((Ds_ENDVERTEX_Z - Bs_ENDVERTEX_Z) < -1) return false;
    
    b_Ds_MM->GetEntry(i);
    if(fabs(Ds_MM-massDs) > 60 ) return false;
    
    //b_Ds_FDCHI2_ORIVX->GetEntry(i);
    //if(Ds_FDCHI2_ORIVX < 0) return false;

    if(_decay== Decay::signal){
	b_Bs_BsDTF_K_1_1270_plus_M->GetEntry(i);
	if(Bs_BsDTF_K_1_1270_plus_M[0] > 2000) return false;
        if(_data){
		b_K_plus_PIDK->GetEntry(i);
        	if(K_plus_PIDK<5) return false;
	}else {
		b_Bs_BKGCAT->GetEntry(i);
		if(Bs_BKGCAT != 20 && Bs_BKGCAT != 60) return false;
	}
    }

    else if(_decay == Decay::norm){
	b_Bs_BsDTF_a_1_1260_plus_M->GetEntry(i);
	if(Bs_BsDTF_a_1_1260_plus_M[0] > 2000) return false;

        if(_data){
		b_pi_plus1_PIDK->GetEntry(i);
        	b_pi_plus2_PIDK->GetEntry(i);
        	if(pi_plus1_PIDK > 5) return false;
        	if(pi_plus2_PIDK > 5) return false;
	}else {
		b_Bs_BKGCAT->GetEntry(i);
		if(Bs_BKGCAT != 20 && Bs_BKGCAT != 60) return false;
	}

    }
    
    return true;
}

inline Bool_t DecayTree::LooseCutsLTU(Long64_t i){

    //b_nPV->GetEntry(i);
    //if(nPV > 5) return false;

    b_Bs_DTF_ctau->GetEntry(i);
    if(Bs_DTF_ctau[0] * 3.33564095 > 0.05) return false;

    b_Bs_MM->GetEntry(i);
    if(Bs_MM < 5200. || Bs_MM > 5700.) return false;
        
    b_Bs_PV_M->GetEntry(i);
    if(Bs_PV_M[0] < 4800. || Bs_PV_M[0] > 6000.) return false;
    
    b_Bs_PV_Dplus_M->GetEntry(i);
    if(Bs_PV_Dplus_M[0] < 1500 || Bs_PV_Dplus_M[0] > 2500) return false; 

    b_Bs_ENDVERTEX_CHI2->GetEntry(i);
    b_Bs_ENDVERTEX_NDOF->GetEntry(i);
    if((Bs_ENDVERTEX_CHI2/Bs_ENDVERTEX_NDOF)> 8) return false;

    b_Bs_PV_chi2->GetEntry(i);
    b_Bs_PV_nDOF->GetEntry(i);
    if(Bs_PV_chi2[0]/Bs_PV_nDOF[0] > 15 )return false;

    b_Ds_IPCHI2_OWNPV->GetEntry(i);
    if(Ds_IPCHI2_OWNPV > 9)return false;

    b_Ds_FDCHI2_ORIVX->GetEntry(i);
    if(Ds_FDCHI2_ORIVX < 9) return false;

    b_Ds_PT->GetEntry(i);
    if(Ds_PT < 1800) return false;

    b_Ds_ENDVERTEX_CHI2->GetEntry(i);
    b_Ds_ENDVERTEX_NDOF->GetEntry(i);
    if((Ds_ENDVERTEX_CHI2/Ds_ENDVERTEX_NDOF)> 5) return false;

    if(_decay== Decay::signal){
        b_K_plus_PIDK->GetEntry(i);
        b_K_plus_isMuon->GetEntry(i);

        b_pi_plus_PIDK->GetEntry(i);
        b_pi_plus_isMuon->GetEntry(i);

        b_pi_minus_PIDK->GetEntry(i);
        b_pi_minus_isMuon->GetEntry(i);
        
	if(K_plus_PIDK < 10) return false;
        else if(pi_plus_PIDK > 5) return false;
        else if(pi_minus_PIDK > 5) return false;
	// remove events with no PID info
	if( fabs(K_plus_PIDK) > 200 ) return false;        
	if( fabs(pi_plus_PIDK) > 200 ) return false;        
	if( fabs(pi_minus_PIDK) > 200 ) return false;   

	if( K_plus_isMuon == 1 ) return false;
 	if( pi_plus_isMuon == 1 ) return false;
 	if( pi_minus_isMuon == 1 ) return false;
      
    }
    
    return true;
}


void DecayTree::Loop()
{

   time_t startTime = time(0);
  
   Init();
   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  // disable all branches
   
   // activate branchname
   fChain->SetBranchStatus("*ID",1);  

   fChain->SetBranchStatus("Bs_*_T*S",1);  
   fChain->SetBranchStatus("Bs_*Muon*_T*S",0);  
   fChain->SetBranchStatus("Bs_*Hlt*Phys*_T*S",0);  
   fChain->SetBranchStatus("Bs_*Hlt*Global*_T*S",0);  
   
   if(!_ltu && !_ss){ 
	fChain->SetBranchStatus("Bs_*TAG*",1);  
	fChain->SetBranchStatus("Bs_*DEC",1);  
	if(_year < 15)fChain->SetBranchStatus("Bs_*PROB",1);  
   }
   fChain->SetBranchStatus("Bs_*DTF*",1);  
   fChain->SetBranchStatus("Bs_*PV*",1);  

   fChain->SetBranchStatus("Bs*ENDVERTEX*",1);  
   fChain->SetBranchStatus("Bs*OWNPV*",1);  
   fChain->SetBranchStatus("Bs*ENDVERTEX_COV*",0);  
   fChain->SetBranchStatus("Bs*OWNPV_COV*",0);  

   fChain->SetBranchStatus("Ds_M",1);         
   fChain->SetBranchStatus("Ds*ENDVERTEX*",1);  
   fChain->SetBranchStatus("Ds*OWNPV*",1);  
   fChain->SetBranchStatus("Ds*ENDVERTEX_COV*",0);  
   fChain->SetBranchStatus("Ds*OWNPV_COV*",0);  
   fChain->SetBranchStatus("Ds*ORIVX*",1);  
   fChain->SetBranchStatus("Ds*ORIVX_COV*",0);  
    
   if(!_ltu){
	fChain->SetBranchStatus("*_1_12*ENDVERTEX*",1);  
	fChain->SetBranchStatus("*_1_12*OWNPV*",1);  
	fChain->SetBranchStatus("*_1_12*ENDVERTEX_COV*",0);  
	fChain->SetBranchStatus("*_1_12*OWNPV_COV*",0);  
	fChain->SetBranchStatus("*_1_12*ORIVX*",1);  
	fChain->SetBranchStatus("*_1_12*ORIVX_COV*",0);  
        fChain->SetBranchStatus("*_1_12*_DOCA*",1);  
	fChain->SetBranchStatus("*Chi2*",1);  
   }
   fChain->SetBranchStatus("*IP*",1);  
   fChain->SetBranchStatus("*IPCHI2*",1);  
   fChain->SetBranchStatus("*FD*",1);  
   fChain->SetBranchStatus("*FDCHI2*",1);  
   fChain->SetBranchStatus("*P",1);  
   fChain->SetBranchStatus("*PT",1);  
   fChain->SetBranchStatus("*PE",1);  
   fChain->SetBranchStatus("*PX",1);  
   fChain->SetBranchStatus("*PY",1);  
   fChain->SetBranchStatus("*PZ",1);  
   fChain->SetBranchStatus("*ETA",1);  
   fChain->SetBranchStatus("*MM*",1);  
   fChain->SetBranchStatus("*TAU*",1);  
   fChain->SetBranchStatus("*ptasy_1.00",1);  

   fChain->SetBranchStatus("*DIRA*",1);  
   fChain->SetBranchStatus("Ds_DOCA*",1);  
    
   fChain->SetBranchStatus("*PID*",1);  
   fChain->SetBranchStatus("*PIDe*",0);  
   fChain->SetBranchStatus("*ProbNN*",1);  
   fChain->SetBranchStatus("*ProbNNe*",0);  
   fChain->SetBranchStatus("*TRACK_Ghost*",1);  
   fChain->SetBranchStatus("*TRACK_CHI2*",1);  
   fChain->SetBranchStatus("*isMuon*",1);  
   fChain->SetBranchStatus("*hasRich",1);  

   fChain->SetBranchStatus("nCandidate",1) ;
   fChain->SetBranchStatus("nTracks",1) ;
   fChain->SetBranchStatus("nPV",1) ;
   fChain->SetBranchStatus("eventNumber",1) ;
   fChain->SetBranchStatus("runNumber",1) ;
   fChain->SetBranchStatus("EventInSequence",1) ;
   fChain->SetBranchStatus("totCandidates",1) ;
   fChain->SetBranchStatus("Polarity",1) ;

   if(!_data){
	   fChain->SetBranchStatus("*TRUE*",1) ;
	   fChain->SetBranchStatus("*BKG*",1) ;
   }

   TFile* output = new TFile(_outFileName,"RECREATE");
   TTree* summary_tree = fChain->CloneTree(0);
    
   Long64_t nentries = fChain->GetEntries();
   cout << "Have " << nentries << " events" <<  endl << endl;

   for (Long64_t i=0; i<nentries;i++) {

      if(0ul == (i % 100000ul)) cout << "Read event " << i << "/" << nentries <<
      "  ( " << i/(double)nentries * 100. << " % )" << endl;
      
      // Read from individual branches rather than whole tree,
      // messy and prone to errors but benefical to performance
      // fChain->GetEntry(i);   

      Long64_t j = LoadTree(i);
      if (j < 0) break;
       
      if(_ltu){
	if(!LooseCutsLTU(j)) continue;
      }
      else {
	if(!TriggerCuts(j)) continue;
	else if(!LooseCuts(j)) continue;
      }

      fChain->GetEntry(i);   
      summary_tree->Fill();
      if(0ul == (i % 100000ul))summary_tree->AutoSave();
   }

    cout << "Selected " << summary_tree->GetEntries() << " events" <<  endl;
    cout << "Efficiency = " << summary_tree->GetEntries()/(double)nentries * 100. << " %" <<  endl;
     
    summary_tree->Write();
    output->Close();
    
    cout << "Created new file " << _outFileName << endl;

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;

}
