// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk), Jack Benton (Jack.B.Benton@bristol.ac.uk)
// status:  Fri 28 Jun 2013 11:21:01 GMT
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/DalitzPdfBase.h"
#include "Mint/DalitzPdfBaseFastInteg.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/AmpRatios.h"
#include "Mint/IEventGenerator.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/SignalGenerator.h"
#include "Mint/FromFileGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/cexp.h"
#include "Mint/DalitzPdfNormChecker.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/DalitzPdfSaveInteg.h"
#include "Mint/Chi2Binning.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpList.h"
#include "Mint/DalitzPdfBaseMCInteg.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include <TStyle.h>
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TROOT.h>
#include "Mint/HyperHistogram.h"
//#include "Mint/GofTests.h"
//#include "Mint/PermutationTest.h"
#include "Mint/LASSO.h"
#include "Mint/LASSO_flexi.h"

using namespace std;
using namespace MINT;


class AmpsPdfFlexiFast
: public DalitzPdfBaseFlexiFastInteg
{
protected:
    TRandom* _localRnd;
    SignalGenerator* _sgGen;
    FromFileGenerator* _fileGen;
    IEventGenerator<IDalitzEvent>* _chosenGen;
    NamedParameter<std::string> _integratorSource;
    std::string _integratorFileName;
public:
    double un_normalised_noPs(IDalitzEvent& evt){
        double ampSq =  _amps->RealVal(evt);
        return ampSq; // * evt.phaseSpace();
    }
    
    std::complex<double> ComplexVal_un_normalised_noPs(IDalitzEvent& evt){
        return  _amps->ComplexVal(evt);
    }
    
    AmpsPdfFlexiFast(const DalitzEventPattern& pat
                     , IFastAmplitudeIntegrable* amps
                     , MinuitParameterSet* mps
                     , double precision=1.e-4
                     , std::string method="efficient"
                     , std::string fname =  "SignalIntegrationEvents.root", bool genMoreEvents = false
                     )
    : DalitzPdfBaseFlexiFastInteg(pat, 0, amps, precision, mps)
    , _localRnd(0)
    , _sgGen(0)
    , _fileGen(0)
    , _chosenGen(0)
    , _integratorSource("IntegratorSource", (std::string) "new", (char*) 0)
    , _integratorFileName(fname)
    {
        cout << " AmpsPdfFlexiFast with integ method " << method << endl;
        bool nonFlat = "efficient" == method;
        bool generateNew = ((string)_integratorSource == (string)"new");
        if(nonFlat){
            cout << "AmpsPdfFlexiFast uses nonFlat integration." << endl;
            if(generateNew){
                _sgGen =  new SignalGenerator(pat,amps);
                _sgGen->setWeighted();
                //_sgGen->setSaveEvents(_integratorFileName);
                _chosenGen = _sgGen;
            }else{
                // here, SignalGenerator is used by FromFileGenerator, to fill
                // up missing events in case more are needed than found in the
                // file.  Since we don't know with which random seed the
                // events in the file were generated, we supply a random
                // number generator with randomised seed.
                _localRnd = new TRandom3(time(0));
                _sgGen =  new SignalGenerator(pat, _localRnd);
                _sgGen->setWeighted();
                _sgGen->dontSaveEvents();// saving events is done by FromFileGenerator
                if(genMoreEvents) _fileGen   = new FromFileGenerator(_integratorFileName, _sgGen);
                else{
                    _fileGen = new FromFileGenerator(_integratorFileName, 0, "UPDATE");
                    cout << "not going to generate any more events" << endl;
                }
                _chosenGen = _fileGen;
            }
            this->setEventGenerator(_chosenGen);
        }else{
            cout << "AmpsPdfFlexiFast uses flat integration." << endl;
        }
    }
    
    IFastAmplitudeIntegrable* getAmpSum(){ return _amps;}
    
    ~AmpsPdfFlexiFast(){
        if(0 != _fileGen)  delete _fileGen;
        if(0 != _sgGen)    delete _sgGen;
        if(0 != _localRnd) delete _localRnd;
    }
};


double getChi2(DalitzEventList& data, DalitzEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 5;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
    //cout << " got event pattern: " << pdg << endl;
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 50);       
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4));
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4));
    HyperCuboid limits(min, max );
                      
    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI,  
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (0)                        
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),                      
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")                 
                         );

    //hist.save("histData.root");
    //HyperBinningHistogram binningHist("histData.root",5);    
    //HyperBinningHistogram dataHist( binningHist.getBinning() );
    //dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc[i];
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;

    return (double)chi2/(nBins-1.);
}

int ampFit(int step=0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    
    NamedParameter<string> InputFileName("InputFileName", (std::string) "");
    NamedParameter<string> InputTreeName("InputTreeName", (std::string) "DalitzEventList");
    std::string inputFile = InputFileName;
    std::string inputTreeName = InputTreeName;
    bool generateNew = (std::string) InputFileName == "";
    std::cout << "InputFileName: " << InputFileName << std::endl;
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile"
                                               , (std::string) "SignalIntegrationEvents.root"
                                               , (char*) 0);

    NamedParameter<string> OutputRootFile("OutputRootFile"
                                          , (std::string) "OutputRootFile.root"
                                          , (char*) 0);
    
    
    NamedParameter<int>  Nevents("Nevents", 1000);
    NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
    NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");
    NamedParameter<int> fitLineshapeParameters("FitLineshapeParameters", 0);
    
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    
    NamedParameter<int>  useLASSO("useLASSO", 1);
    NamedParameter<double>  lambda("lambda", 1.);
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
        
    DalitzEventList eventNorm;
    eventNorm.generatePhaseSpaceEvents(1,pat); 
    
    FitAmpIncoherentSum fas(pat);
    fas.getVal(eventNorm[0]);
    fas.print();
        
//     {
//         DalitzEventList eventNorm2;
//          eventNorm2.generatePhaseSpaceEvents(200000,pat); 
//          fas.normalizeAmps(eventNorm2);
//     }
        
    DalitzEventList eventList;
    
    if(generateNew){
        SignalGenerator sg(pat,&fas);
        cout << "Generating " << Nevents << " MC events." << endl;
        sg.FillEventList(eventList, Nevents);
        eventList.saveAsNtuple(OutputRootFile);
    }
    
    if(!generateNew){
        TFile *_InputFile =  TFile::Open(inputFile.c_str());
        TTree* in_tree;
        in_tree=dynamic_cast<TTree*>(_InputFile->Get(inputTreeName.c_str()));
        cout << "reading events from file " << inputFile << endl;
        eventList.fromNtuple(in_tree,1);
        cout << " I've got " << eventList.size() << " events." << endl;
        _InputFile->Close();
    }
    
    DalitzHistoSet datH = eventList.weightedHistoSet();
    
        AmpsPdfFlexiFast amps(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);        
        Neg2LL neg2ll(amps, eventList);
        
        double stepSize = 1;
        lambda = lambda + (step-1) * stepSize;
        LASSO_flexi lasso(&amps,lambda);
        Neg2LLSum fcn(&neg2ll,&lasso);
        
        Minimiser mini;
        if(useLASSO)mini.attachFunction(&fcn);
        else mini.attachFunction(&neg2ll);
        mini.doFit();
        
        mini.printResultVsInput();
          amps.doFinalStats(&mini);
                
            cout << "Now plotting:" << endl;
            
            DalitzHistoSet fitH = amps.histoSet();
            datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_l_"+anythingToString(step)+"_").c_str(),"eps");
            std::vector<DalitzHistoSet> EachAmpsHistos = amps.GetEachAmpsHistograms();
            datH.drawWithFitAndEachAmps(datH, fitH, EachAmpsHistos, ((string)OutputDir+(string)"WithAmps").c_str(), "eps");
            
            int nBins = 100;
            vector<int> s123;
            s123.push_back(1);
            s123.push_back(2);
            s123.push_back(3);
            
            vector<int> s234;
            s234.push_back(2);
            s234.push_back(3);
            s234.push_back(4);
            
            vector<int> s134;
            s134.push_back(1);
            s134.push_back(3);
            s134.push_back(4);
            
            vector<int> s124;
            s124.push_back(1);
            s124.push_back(2);
            s124.push_back(4);
            
            TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.4,12);
            TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,6);
            TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,10);
            TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
            TH1D* s_Dspi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            TH1D* s_Dspim = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            
	    TH1D* s_Kpipi_rw = (TH1D*)s_Kpipi->Clone();
	    TH1D* s_Kpi_rw = (TH1D*)s_Kpi->Clone();
	    TH1D* s_pipi_rw = (TH1D*)s_pipi->Clone();
	    TH1D* s_Dspipi_rw = (TH1D*)s_Dspipi->Clone();
	    TH1D* s_DsK_rw = (TH1D*)s_DsK->Clone();
	    TH1D* s_DsKpi_rw = (TH1D*)s_DsKpi->Clone();
	    TH1D* s_Dspi_rw = (TH1D*)s_Dspi->Clone();
	    TH1D* s_Dspim_rw = (TH1D*)s_Dspim->Clone();

	    TH1D* s_Kpipi_phsp = (TH1D*)s_Kpipi->Clone();
	    TH1D* s_Kpi_phsp = (TH1D*)s_Kpi->Clone();
	    TH1D* s_pipi_phsp = (TH1D*)s_pipi->Clone();
	    TH1D* s_Dspipi_phsp = (TH1D*)s_Dspipi->Clone();
	    TH1D* s_DsK_phsp = (TH1D*)s_DsK->Clone();
	    TH1D* s_DsKpi_phsp = (TH1D*)s_DsKpi->Clone();
	    TH1D* s_Dspi_phsp = (TH1D*)s_Dspi->Clone();
	    TH1D* s_Dspim_phsp = (TH1D*)s_Dspim->Clone();

            TH2D* s_Kpi_pipi = new TH2D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4}) ",60,0.2,1.6,60,0.,1.6);
            TH2D* s_DsKpi_Dspi = new TH2D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4}); #left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4}) ",80,5,30,80,0,25);
            TH2D* s_DsK_Dspi = new TH2D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4}); #left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4}) ",60,0,30,60,0,25);

	    int counter = 0;

            for (int i=0; i<eventList.size(); i++) {
		//if(eventList[i].phaseSpace()==0.000000) continue;
 		//else counter++;
		//eventList[i].setWeight(eventList[i].getWeight()/eventList[i].phaseSpace());
		s_Kpipi->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
                s_Kpi->Fill(eventList[i].s(2,4)/(GeV*GeV),eventList[i].getWeight());
                s_pipi->Fill(eventList[i].s(3,4)/(GeV*GeV),eventList[i].getWeight());
                s_Dspipi->Fill(eventList[i].sij(s134)/(GeV*GeV),eventList[i].getWeight());
                s_DsK->Fill(eventList[i].s(1,2)/(GeV*GeV),eventList[i].getWeight());
                s_DsKpi->Fill(eventList[i].sij(s124)/(GeV*GeV),eventList[i].getWeight());
                s_Dspi->Fill(eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());
                s_Dspim->Fill(eventList[i].s(1,4)/(GeV*GeV),eventList[i].getWeight());

                s_Kpi_pipi->Fill(eventList[i].s(2,4)/(GeV*GeV),eventList[i].s(3,4)/(GeV*GeV),eventList[i].getWeight());
                s_DsKpi_Dspi->Fill(eventList[i].sij(s124)/(GeV*GeV),eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());
                s_DsK_Dspi->Fill(eventList[i].s(1,2)/(GeV*GeV),eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());

// 		if(amps.RealVal(eventList[i])>0.)
		double weight;
		if(fas.getVal(eventList[i])>0.)weight = eventList[i].getWeight()/fas.getVal(eventList[i]);
		else weight = 0.;
// 		else { counter ++; eventList[i].setWeight(0);} 
		s_Kpipi_rw->Fill(eventList[i].sij(s234)/(GeV*GeV),weight);
                s_Kpi_rw->Fill(eventList[i].s(2,4)/(GeV*GeV),weight);
                s_pipi_rw->Fill(eventList[i].s(3,4)/(GeV*GeV),weight);
                s_Dspipi_rw->Fill(eventList[i].sij(s134)/(GeV*GeV),weight);
                s_DsK_rw->Fill(eventList[i].s(1,2)/(GeV*GeV),weight);
                s_DsKpi_rw->Fill(eventList[i].sij(s124)/(GeV*GeV),weight);
                s_Dspi_rw->Fill(eventList[i].s(1,3)/(GeV*GeV),weight);
                s_Dspim_rw->Fill(eventList[i].s(1,4)/(GeV*GeV),weight);

// 		eventList[i].setWeight(weight);
            }    

	    cout << "bad evt = " << counter << endl;            
            TH1D* s_Kpipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.4,12);
            TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,6);
            TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,10);
            TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsKpi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
            TH1D* s_Dspi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            TH1D* s_Dspim_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
                       
            //SignalGenerator sg(pat,&fas);
            //sg.setWeighted();
	    DalitzEventList eventListMC;
	    TFile *FileMC =  TFile::Open(((string) IntegratorEventFile).c_str());
	    TTree* treeMC = dynamic_cast<TTree*>(FileMC->Get("DalitzEventList"));
	    eventListMC.fromNtuple(treeMC,1);
	    FileMC->Close();
            
            for(int i = 0; i < eventListMC.size(); i++){
                                
                //counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
                //DalitzEvent evt(evtPtr.get());
                double weight = fas.getVal(eventListMC[i])*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                //double weight = eventListMC[i].phaseSpace()*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
		s_Kpipi_fit->Fill(eventListMC[i].sij(s234)/(GeV*GeV),weight);
                s_Kpi_fit->Fill(eventListMC[i].s(2,4)/(GeV*GeV),weight);
                s_pipi_fit->Fill(eventListMC[i].s(3,4)/(GeV*GeV),weight);
                s_Dspipi_fit->Fill(eventListMC[i].sij(s134)/(GeV*GeV),weight);
                s_DsK_fit->Fill(eventListMC[i].s(1,2)/(GeV*GeV),weight);
                s_DsKpi_fit->Fill(eventListMC[i].sij(s124)/(GeV*GeV),weight);
                s_Dspi_fit->Fill(eventListMC[i].s(1,3)/(GeV*GeV),weight);
                s_Dspim_fit->Fill(eventListMC[i].s(1,4)/(GeV*GeV),weight);

// 		weight = eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
// 		s_Kpipi_phsp->Fill(eventListMC[i].sij(s234)/(GeV*GeV),weight);
//                 s_Kpi_phsp->Fill(eventListMC[i].s(2,4)/(GeV*GeV),weight);
//                 s_pipi_phsp->Fill(eventListMC[i].s(3,4)/(GeV*GeV),weight);
//                 s_Dspipi_phsp->Fill(eventListMC[i].sij(s134)/(GeV*GeV),weight);
//                 s_DsK_phsp->Fill(eventListMC[i].s(1,2)/(GeV*GeV),weight);
//                 s_DsKpi_phsp->Fill(eventListMC[i].sij(s124)/(GeV*GeV),weight);
//                 s_Dspi_phsp->Fill(eventListMC[i].s(1,3)/(GeV*GeV),weight);
//                 s_Dspim_phsp->Fill(eventListMC[i].s(1,4)/(GeV*GeV),weight);
		eventListMC[i].setWeight(weight);
            }
             
	    DalitzEventList eventListPhsp;
	    TFile *FilePhsp =  TFile::Open("SignalIntegrationEvents_Phsp.root");
	    TTree* treePhsp = dynamic_cast<TTree*>(FilePhsp->Get("DalitzEventList"));
	    eventListPhsp.fromNtuple(treePhsp,1);
	    FilePhsp->Close();
            
            for(int i = 0; i < eventListPhsp.size(); i++){                                
                double weight = 1.;//eventListPhsp[i].getWeight()/eventListPhsp[i].getGeneratorPdfRelativeToPhaseSpace();
		s_Kpipi_phsp->Fill(eventListPhsp[i].sij(s234)/(GeV*GeV),weight);
                s_Kpi_phsp->Fill(eventListPhsp[i].s(2,4)/(GeV*GeV),weight);
                s_pipi_phsp->Fill(eventListPhsp[i].s(3,4)/(GeV*GeV),weight);
                s_Dspipi_phsp->Fill(eventListPhsp[i].sij(s134)/(GeV*GeV),weight);
                s_DsK_phsp->Fill(eventListPhsp[i].s(1,2)/(GeV*GeV),weight);
                s_DsKpi_phsp->Fill(eventListPhsp[i].sij(s124)/(GeV*GeV),weight);
                s_Dspi_phsp->Fill(eventListPhsp[i].s(1,3)/(GeV*GeV),weight);
                s_Dspim_phsp->Fill(eventListPhsp[i].s(1,4)/(GeV*GeV),weight);
            }


            TCanvas* c = new TCanvas();
                     
	    TF1* f = new TF1("f","1.0",0,30);
	    f->SetLineColor(kBlue);

            s_Kpi_pipi->SetMinimum(0.);
            s_Kpi_pipi->Draw("colz");
            c->Print(((string)OutputDir+"s_Kpi_pipi.eps").c_str());
            s_Kpi_pipi->Draw();
            c->Print(((string)OutputDir+"s_Kpi_pipi_scatter.eps").c_str());

            s_DsKpi_Dspi->SetMinimum(0);
            s_DsKpi_Dspi->Draw("colz");
            c->Print(((string)OutputDir+"s_DsKpi_Dspi.eps").c_str());

            s_DsK_Dspi->SetMinimum(0);
            s_DsK_Dspi->Draw("colz");
            c->Print(((string)OutputDir+"s_DsK_Dspi.eps").c_str());
            s_DsK_Dspi->Draw();
            c->Print(((string)OutputDir+"s_DsK_Dspi_scatter.eps").c_str());

  	    gPad->SetLogy(1);
            s_Kpipi->SetMinimum(0.1);
            s_Kpipi->SetLineColor(kBlack);
            s_Kpipi->DrawNormalized("e1",1);
            s_Kpipi_fit->SetLineColor(kBlue);
            s_Kpipi_fit->SetLineWidth(3);
            s_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi.eps").c_str());

            s_Kpipi_rw->SetMinimum(0.1);
            s_Kpipi_phsp->SetMinimum(0.1);
            s_Kpipi_rw->SetLineColor(kBlack);
            s_Kpipi_phsp->SetLineColor(kBlue);
            s_Kpipi_phsp->SetLineWidth(3);
            s_Kpipi_phsp->DrawNormalized("histc",1);
            s_Kpipi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_Kpipi.eps").c_str());

	    s_Kpipi_rw->Scale(1./s_Kpipi_rw->Integral());
    	    s_Kpipi_phsp->Scale(1./s_Kpipi_phsp->Integral());
   	    s_Kpipi_rw->Divide(s_Kpipi_rw,s_Kpipi_phsp);
    	    s_Kpipi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Kpipi.eps").c_str());


 	    s_Kpipi_fit->Scale(1./s_Kpipi_fit->Integral());
     	    s_Kpipi->Scale(1./s_Kpipi->Integral());
    	    s_Kpipi_fit->Divide(s_Kpipi_fit,s_Kpipi);
     	    s_Kpipi_fit->Draw("e");
     	    c->Print(((string)OutputDir+"eff2_Kpipi.eps").c_str());

            s_Kpipi->DrawNormalized("e1",1);
            c->Print(((string)OutputDir+"s_Kpipi_data.eps").c_str());
            
            s_Kpi->SetMinimum(0.1);
            s_Kpi->SetLineColor(kBlack);
            s_Kpi->DrawNormalized("e1",1);
            s_Kpi_fit->SetLineColor(kBlue);
            s_Kpi_fit->SetLineWidth(3);
            s_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi.eps").c_str());

            s_Kpi_rw->SetMinimum(0.1);
            s_Kpi_phsp->SetMinimum(0.1);
            s_Kpi_rw->SetLineColor(kBlack);
            s_Kpi_phsp->SetLineColor(kBlue);
            s_Kpi_phsp->SetLineWidth(3);
            s_Kpi_phsp->DrawNormalized("histc",1);
            s_Kpi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_Kpi.eps").c_str());

	    s_Kpi_rw->Scale(1./s_Kpi_rw->Integral());
    	    s_Kpi_phsp->Scale(1./s_Kpi_phsp->Integral());
   	    s_Kpi_rw->Divide(s_Kpi_rw,s_Kpi_phsp);
    	    s_Kpi_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_Kpi.eps").c_str());
            
	    s_pipi->SetMinimum(0.1);            
            s_pipi->SetLineColor(kBlack);
            s_pipi->DrawNormalized("e1",1);
            s_pipi_fit->SetLineColor(kBlue);
            s_pipi_fit->SetLineWidth(3);
            s_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi.eps").c_str());

            s_pipi_rw->SetMinimum(0.1);
            s_pipi_phsp->SetMinimum(0.1);
            s_pipi_rw->SetLineColor(kBlack);
            s_pipi_phsp->SetLineColor(kBlue);
            s_pipi_phsp->SetLineWidth(3);
            s_pipi_phsp->DrawNormalized("histc",1);
            s_pipi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_pipi.eps").c_str());

	    s_pipi_rw->Scale(1./s_pipi_rw->Integral());
    	    s_pipi_phsp->Scale(1./s_pipi_phsp->Integral());
   	    s_pipi_rw->Divide(s_pipi_rw,s_pipi_phsp);
    	    s_pipi_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_pipi.eps").c_str());
            
	    s_Dspipi->SetMinimum(0.1);            
            s_Dspipi->SetLineColor(kBlack);
            s_Dspipi->DrawNormalized("e1",1);
            s_Dspipi_fit->SetLineColor(kBlue);
            s_Dspipi_fit->SetLineWidth(3);
            s_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi.eps").c_str());

            s_Dspipi_rw->SetMinimum(0.1);
            s_Dspipi_phsp->SetMinimum(0.1);
            s_Dspipi_rw->SetLineColor(kBlack);
            s_Dspipi_phsp->SetLineColor(kBlue);
            s_Dspipi_phsp->SetLineWidth(3);
            s_Dspipi_phsp->DrawNormalized("histc",1);
            s_Dspipi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_Dspipi.eps").c_str());

	    s_Dspipi_rw->Scale(1./s_Dspipi_rw->Integral());
    	    s_Dspipi_phsp->Scale(1./s_Dspipi_phsp->Integral());
   	    s_Dspipi_rw->Divide(s_Dspipi_rw,s_Dspipi_phsp);
    	    s_Dspipi_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_Dspipi.eps").c_str());
           
	    s_DsK->SetMinimum(0.1);
	    s_DsK_fit->SetMinimum(0.1);
            s_DsK->SetLineColor(kBlack);
            s_DsK->DrawNormalized("e1",1);
            s_DsK_fit->SetLineColor(kBlue);
            s_DsK_fit->SetLineWidth(3);
            s_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK.eps").c_str());

            s_DsK_rw->SetMinimum(0.1);
            s_DsK_phsp->SetMinimum(0.1);
            s_DsK_rw->SetLineColor(kBlack);
            s_DsK_phsp->SetLineColor(kBlue);
            s_DsK_phsp->SetLineWidth(3);
            s_DsK_phsp->DrawNormalized("histc",1);
            s_DsK_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_DsK.eps").c_str());

	    s_DsK_rw->Scale(1./s_DsK_rw->Integral());
    	    s_DsK_phsp->Scale(1./s_DsK_phsp->Integral());
   	    s_DsK_rw->Divide(s_DsK_rw,s_DsK_phsp);
    	    s_DsK_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_DsK.eps").c_str());

 	    s_DsK_fit->Scale(1./s_DsK_fit->Integral());
     	    s_DsK->Scale(1./s_DsK->Integral());
    	    s_DsK_fit->Divide(s_DsK_fit,s_DsK);
     	    s_DsK_fit->Draw("e");
     	    c->Print(((string)OutputDir+"eff2_DsK.eps").c_str());

	    s_DsKpi->SetMinimum(0.1);            
            s_DsKpi->SetLineColor(kBlack);
            s_DsKpi->DrawNormalized("e1",1);
            s_DsKpi_fit->SetLineColor(kBlue);
            s_DsKpi_fit->SetLineWidth(3);
            s_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi.eps").c_str());

            s_DsKpi_rw->SetMinimum(0.1);
            s_DsKpi_phsp->SetMinimum(0.1);
            s_DsKpi_rw->SetLineColor(kBlack);
            s_DsKpi_phsp->SetLineColor(kBlue);
            s_DsKpi_phsp->SetLineWidth(3);
            s_DsKpi_phsp->DrawNormalized("histc",1);
            s_DsKpi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_DsKpi.eps").c_str());

	    s_DsKpi_rw->Scale(1./s_DsKpi_rw->Integral());
    	    s_DsKpi_phsp->Scale(1./s_DsKpi_phsp->Integral());
   	    s_DsKpi_rw->Divide(s_DsKpi_rw,s_DsKpi_phsp);
    	    s_DsKpi_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_DsKpi.eps").c_str());

	    s_Dspi->SetMinimum(0.1);
            s_Dspi->SetLineColor(kBlack);
            s_Dspi->DrawNormalized("e1",1);
            s_Dspi_fit->SetLineColor(kBlue);
            s_Dspi_fit->SetLineWidth(3);
            s_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi.eps").c_str());

            s_Dspi_rw->SetMinimum(0.1);
            s_Dspi_phsp->SetMinimum(0.1);
            s_Dspi_rw->SetLineColor(kBlack);
            s_Dspi_phsp->SetLineColor(kBlue);
            s_Dspi_phsp->SetLineWidth(3);
            s_Dspi_phsp->DrawNormalized("histc",1);
            s_Dspi_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_Dspi.eps").c_str());

	    s_Dspi_rw->Scale(1./s_Dspi_rw->Integral());
    	    s_Dspi_phsp->Scale(1./s_Dspi_phsp->Integral());
   	    s_Dspi_rw->Divide(s_Dspi_rw,s_Dspi_phsp);
    	    s_Dspi_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_Dspi.eps").c_str());

	    s_Dspim->SetMinimum(0.1);
            s_Dspim->SetLineColor(kBlack);
            s_Dspim->DrawNormalized("e1",1);
            s_Dspim_fit->SetLineColor(kBlue);
            s_Dspim_fit->SetLineWidth(3);
            s_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim.eps").c_str());

            s_Dspim_rw->SetMinimum(0.1);
            s_Dspim_phsp->SetMinimum(0.1);
            s_Dspim_rw->SetLineColor(kBlack);
            s_Dspim_phsp->SetLineColor(kBlue);
            s_Dspim_phsp->SetLineWidth(3);
            s_Dspim_phsp->DrawNormalized("histc",1);
            s_Dspim_rw->DrawNormalized("esame",1);
            c->Print(((string)OutputDir+"phsp_Dspim.eps").c_str());

	    s_Dspim_rw->Scale(1./s_Dspim_rw->Integral());
    	    s_Dspim_phsp->Scale(1./s_Dspim_phsp->Integral());
   	    s_Dspim_rw->Divide(s_Dspim_rw,s_Dspim_phsp);
    	    s_Dspim_rw->Draw("e");
    	    c->Print(((string)OutputDir+"eff_Dspim.eps").c_str());

 	    getChi2(eventList,eventListMC);        
        
    return 0;
}

void makeIntegratorFile(){
    
    FitAmplitude::AutogenerateFitFile();

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventListPhsp,eventList,eventList_cut;
    eventListPhsp.generatePhaseSpaceEvents(10,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.print();
    fas.getVal(eventListPhsp[0]);
    //fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);
    sg.FillEventList(eventList, IntegratorEvents);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    //for(int i = 0; i < eventList.size(); i++){
	//if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 1.95 && sqrt(eventList[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventList[i].s(3,4)/(GeV*GeV)) < 1.2)
	//eventList_cut.Add(eventList[i]);
    //}
    //cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    eventList.saveAsNtuple(IntegratorEventFile);
    return;
}


void makeIntegratorFilePhsp(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventList;

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    int i = 0;
    while(i < IntegratorEvents){
	DalitzEvent evt(pat);
	evt.generateThisToPhaseSpace();
	//if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2){ 
		eventList.Add(evt);
		i++;
	//}
    }

    cout << "Generated " << eventList.size() << " events inside selected phasespace region" << endl;
    
    eventList.saveAsNtuple(IntegratorEventFile);
    return;
}


int makeMINTtupleGen(){
    
    string outputDir = "/auto/data/dargent/BsDsKpipi/MINT/";

    bool dbThis=false;
    bool addSweight = true;
    bool bkg = false;
    
    int N=-1;
    if(dbThis) cout << "read ntuple" << endl;
	
    NamedParameter<int> EventPattern("Event Pattern", 521, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
    cout << " got event pattern: " << pdg << endl;
	
    DalitzEventList eventList; 

    // Read the momenta from ntuple
    TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
    tree_gen->Add("EvtGen.root");

   
    if (dbThis) cout << "Read the file" << endl;	

    double K_gen[5]; 
    double pip_gen[5]; 
    double pim_gen[5]; 
    double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];
    
    tree_gen->SetBranchAddress("Kplus_TRUEP_X",&K_gen[0]);
    tree_gen->SetBranchAddress("Kplus_TRUEP_Y",&K_gen[1]);
    tree_gen->SetBranchAddress("Kplus_TRUEP_Z",&K_gen[2]); 
    tree_gen->SetBranchAddress("Kplus_TRUEP_E",&K_gen[3]); 
    tree_gen->SetBranchAddress("Kplus_TRUEPT",&K_gen[4]); 
	
    tree_gen->SetBranchAddress("piplus_TRUEP_X",&pip_gen[0]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Y",&pip_gen[1]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Z",&pip_gen[2]); 
    tree_gen->SetBranchAddress("piplus_TRUEP_E",&pip_gen[3]); 
    tree_gen->SetBranchAddress("piplus_TRUEPT",&pip_gen[4]); 

    tree_gen->SetBranchAddress("piminus_TRUEP_X",&pim_gen[0]);
    tree_gen->SetBranchAddress("piminus_TRUEP_Y",&pim_gen[1]);
    tree_gen->SetBranchAddress("piminus_TRUEP_Z",&pim_gen[2]); 
    tree_gen->SetBranchAddress("piminus_TRUEP_E",&pim_gen[3]); 
    tree_gen->SetBranchAddress("piminus_TRUEPT",&pim_gen[4]); 
	
    tree_gen->SetBranchAddress("Kplus0_TRUEP_X",&Ds_Kp_gen[0]);
    tree_gen->SetBranchAddress("Kplus0_TRUEP_Y",&Ds_Kp_gen[1]);
    tree_gen->SetBranchAddress("Kplus0_TRUEP_Z",&Ds_Kp_gen[2]); 
    tree_gen->SetBranchAddress("Kplus0_TRUEP_E",&Ds_Kp_gen[3]); 
    tree_gen->SetBranchAddress("Kplus0_TRUEPT",&Ds_Kp_gen[4]); 
    
    tree_gen->SetBranchAddress("Kminus_TRUEP_X",&Ds_Km_gen[0]);
    tree_gen->SetBranchAddress("Kminus_TRUEP_Y",&Ds_Km_gen[1]);
    tree_gen->SetBranchAddress("Kminus_TRUEP_Z",&Ds_Km_gen[2]); 
    tree_gen->SetBranchAddress("Kminus_TRUEP_E",&Ds_Km_gen[3]); 
    tree_gen->SetBranchAddress("Kminus_TRUEPT",&Ds_Km_gen[4]); 

    tree_gen->SetBranchAddress("piminus0_TRUEP_X",&Ds_pim_gen[0]);
    tree_gen->SetBranchAddress("piminus0_TRUEP_Y",&Ds_pim_gen[1]);
    tree_gen->SetBranchAddress("piminus0_TRUEP_Z",&Ds_pim_gen[2]); 
    tree_gen->SetBranchAddress("piminus0_TRUEP_E",&Ds_pim_gen[3]); 
    tree_gen->SetBranchAddress("piminus0_TRUEPT",&Ds_pim_gen[4]); 
    
    int numEvents = tree_gen->GetEntries();
    int numSelected =0;

    //loop over tree and fill eventList
    for(int i=0; i< numEvents; i++)
    {
	if(dbThis)cout << " getting " << i << " th entry" << endl;	
	tree_gen->GetEntry(i);
        
        // Lorentz vectors: P=(Px,Py,Pz,E)
        TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
        TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
	TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
        TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
        TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
        TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 


	// define the order of the vectors in the vectorOfvectors
        // include the 'MeV' to get the correct units, need to include CLHEPSystemOfUnits.h
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
	vectorOfvectors.push_back(pip_p*MeV);
	vectorOfvectors.push_back(pim_p*MeV);

	if(dbThis) cout << "make event" << endl;
		
	DalitzEvent evt(pdg, vectorOfvectors);
	//if(evt.phaseSpace()==0) cout << evt << endl;
	/*
	if(evt.s(1,2)<0 || evt.s(2,3)<0 || evt.s(3,4)<0 || evt.t(4,0)<0 || evt.t(1,0)< 0) 
	cout << "negative mass ?" << endl; 
        if(evt.s(1,2)< pdg.sijMin(1,2) || evt.s(1,2)> pdg.sijMax(1,2) ) continue;
        if(evt.s(2,3)< pdg.sijMin(2,3) || evt.s(2,3)> pdg.sijMax(2,3) )continue; 
        if(evt.s(3,4)< pdg.sijMin(3,4) || evt.s(3,4)> pdg.sijMax(3,4) ) continue;
        if(evt.t(0,4)< pdg.sijMin(1,2,3) || evt.t(4,0)> pdg.sijMax(1,2,3) ) continue;  
        if(evt.t(1,0)< pdg.sijMin(2,3,4) || evt.t(1,0)> pdg.sijMax(2,3,4) ) continue; 	
	if(dbThis) cout << "s12 =  " << ( evt.p(1) + evt.p(2) ).Mag2() << endl;
	if(dbThis) cout << "adding event " << evt << endl;
	*/
 	//if(abs(D_p.M()-1968.2) < 0.2)
	eventList.Add(evt); // this fills the event list		
	if(dbThis) cout << " added event" << endl;

// 	if(D_p.M()-1968.2 > 0.2) cout <<  D_p.M() << endl;
		
        numSelected++;
        if(numSelected==N)break;
    }
    
    TString output = outputDir + "GenMC";
    if(!addSweight){
	if(bkg) output+="_bkg";
	else output+="_3sigma";
    }
    if(N != -1)output += "_small";
    output+=".root";
    
    eventList.save((string)output);
    DalitzHistoSet datH = eventList.weightedHistoSet();
    datH.draw("data_","eps"); 
   
    cout << numSelected << " / " << numEvents << " events selected" << endl;
    cout << "Created File: " << output << endl;    

    return 0;
}

void reweightEvtGen(){

  NamedParameter<int> EventPattern("Event Pattern", 421, -321, 211, 211, -211);
  DalitzEventPattern pat(EventPattern);

  FitAmpIncoherentSum fas(pat);

  DalitzEventList eventListMC;
  TFile *FileMC =  TFile::Open("/auto/data/dargent/BsDsKpipi/MINT/GenMC.root");
  TTree* treeMC = dynamic_cast<TTree*>(FileMC->Get("DalitzEventList"));
  eventListMC.fromNtuple(treeMC,1);
  FileMC->Close();
            
  for(int i = 0; i < eventListMC.size(); i++){                                
		eventListMC[i].setGeneratorPdfRelativeToPhaseSpace(fas.getVal(eventListMC[i]));
  }
  eventListMC.save("/auto/data/dargent/BsDsKpipi/MINT/GenMC_rw.root");
}

class FracLL : public Minimisable{
  FlexiFastAmplitudeIntegrator* _integ;
public:
  FracLL(FlexiFastAmplitudeIntegrator* integ) : _integ(integ){
    _integ->getVal();//makes sure we are initialised
  }
  double getVal(){
    return _integ->getFractionChi2();
  }
};

int fracFit(){

  NamedParameter<double> IntegPrecision("IntegPrecision", 1.e-3);
  NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);  
  NamedParameter<int> EventPattern("Event Pattern", 421, -321, 211, 211, -211);
  DalitzEventPattern pat(EventPattern);

  FitAmpIncoherentSum fas(pat);
  //SignalGenerator sg(pat,&fas);
  //sg.setWeighted();
  FromFileGenerator* fileGen = new FromFileGenerator(IntegratorEventFile, 0, "UPDATE");

  FlexiFastAmplitudeIntegrator integ(pat, &fas, fileGen, gRandom, IntegPrecision);

  cout << "integrator value: " << integ.getVal() << endl;
  cout << "now doing the fit" << endl;

  FracLL f(&integ);
  Minimiser mini(&f);
  mini.doFit();
  mini.printResultVsInput();
  integ.doFinalStats(&mini);

  return 0;
}


int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    if(! std::ifstream(((string)IntegratorEventFile).c_str()).good()) makeIntegratorFile();
  
//      makeMINTtupleGen();
//     ampFit(atoi(argv[1]));

    fracFit();
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
//
