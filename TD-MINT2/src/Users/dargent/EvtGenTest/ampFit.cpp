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

double getChi2(DiskResidentEventList& data, DiskResidentEventList& mc){
	
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
        DalitzEvent evt = data.getEvent(i);
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
	if(!(evt.phaseSpace() > 0.))continue;
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::MINT,  
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        ((pdg.sijMax(1,3)-pdg.sijMin(1,3))/100.),                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (4)                        
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),                      
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")                 
                         );

    dataHist.save("histData.root");

//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
	if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }


    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    mcHist.save("histMC.root");

    cout << mcHist.integral() << endl;

    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << dataHist.integral() << endl;
    cout << mcHist.integral() << endl;

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
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    
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
            
            TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,5);
            TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,4);
            TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
            TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
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

	    TH1D* s_Kpipi_pull = (TH1D*)s_Kpipi->Clone();
	    TH1D* s_Kpipi_pull2 = (TH1D*)s_Kpipi->Clone();

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

            TH1D* s_Kpipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,5);
            TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.,4);
            TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,4);
            TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_DsKpi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
            TH1D* s_Dspi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            TH1D* s_Dspim_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);

// {       
//   	DalitzEventList eventList;
//         TFile *_InputFile =  TFile::Open(inputFile.c_str());
//         TTree* in_tree;
//         in_tree=dynamic_cast<TTree*>(_InputFile->Get(inputTreeName.c_str()));
//         cout << "reading events from file " << inputFile << endl;
//         eventList.fromNtuple(in_tree,0.2);
//         cout << " I've got " << eventList.size() << " events." << endl;
//         _InputFile->Close();
// 
//          AmpsPdfFlexiFast amps(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);        
//          Neg2LL neg2ll(amps, eventList);
//          Minimiser mini(&neg2ll);    
//  
//          mini.doFit();        
//          mini.printResultVsInput();
//          amps.doFinalStats(&mini);
// }  
              
            //DalitzHistoSet fitH = amps.histoSet();
            //datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_l_"+anythingToString(step)+"_").c_str(),"eps");
            //std::vector<DalitzHistoSet> EachAmpsHistos = amps.GetEachAmpsHistograms();
            //datH.drawWithFitAndEachAmps(datH, fitH, EachAmpsHistos, ((string)OutputDir+(string)"WithAmps").c_str(), "eps");
	    int badEvents = 0;
	    double sumw = 0.;
	    double sumw2 = 0.;

            DiskResidentEventList eventList(inputFile.c_str(),"OPEN");
            for (int i=0; i<eventList.size(); i++) {
		DalitzEvent evt(eventList.getEvent(i));

 		if(evt.sij(s234)/(GeV*GeV)> 4) continue;
		if(!(evt.phaseSpace() > 0.)){
		 	//cout << "evt " << i << " 0 phsp " << endl << evt << endl;
			badEvents++;
			continue;
		}
		if(TMath::IsNaN(fas.getVal(evt))){
		 	//cout << "evt " << i << " isNaN " << endl << evt << endl;
			badEvents++;
			continue;
		}

		s_Kpipi->Fill(evt.sij(s234)/(GeV*GeV),evt.getWeight());
                s_Kpi->Fill(evt.s(2,4)/(GeV*GeV),evt.getWeight());
                s_pipi->Fill(evt.s(3,4)/(GeV*GeV),evt.getWeight());
                s_Dspipi->Fill(evt.sij(s134)/(GeV*GeV),evt.getWeight());
                s_DsK->Fill(evt.s(1,2)/(GeV*GeV),evt.getWeight());
                s_DsKpi->Fill(evt.sij(s124)/(GeV*GeV),evt.getWeight());
                s_Dspi->Fill(evt.s(1,3)/(GeV*GeV),evt.getWeight());
                s_Dspim->Fill(evt.s(1,4)/(GeV*GeV),evt.getWeight());

                s_Kpi_pipi->Fill(evt.s(2,4)/(GeV*GeV),evt.s(3,4)/(GeV*GeV),evt.getWeight());
                s_DsKpi_Dspi->Fill(evt.sij(s124)/(GeV*GeV),evt.s(1,3)/(GeV*GeV),evt.getWeight());
                s_DsK_Dspi->Fill(evt.s(1,2)/(GeV*GeV),evt.s(1,3)/(GeV*GeV),evt.getWeight());

		double weight;
		if(fas.getVal(evt)>0.)weight = evt.getWeight()/fas.getVal(evt);
		else weight = 0.;

		s_Kpipi_rw->Fill(evt.sij(s234)/(GeV*GeV),weight);
                s_Kpi_rw->Fill(evt.s(2,4)/(GeV*GeV),weight);
                s_pipi_rw->Fill(evt.s(3,4)/(GeV*GeV),weight);
                s_Dspipi_rw->Fill(evt.sij(s134)/(GeV*GeV),weight);
                s_DsK_rw->Fill(evt.s(1,2)/(GeV*GeV),weight);
                s_DsKpi_rw->Fill(evt.sij(s124)/(GeV*GeV),weight);
                s_Dspi_rw->Fill(evt.s(1,3)/(GeV*GeV),weight);
                s_Dspim_rw->Fill(evt.s(1,4)/(GeV*GeV),weight);
// 		evt.setWeight(weight);

		sumw += weight;
		sumw2 += weight*weight;
            }    
	    cout << endl << "bad EvtGen events " << badEvents << " ( " << badEvents/(double)eventList.size() << " %)" << endl << endl;

	    //cout << "sumw = " << sumw << endl;
	    //cout << "weff = " << sumw/sumw2 << endl;
	    //cout << "Neff = " << sumw*sumw/sumw2 << endl;
		
	    sumw =  0.;
	    sumw2 = 0.;
	    for(int i = 1; i <= s_Kpipi_rw->GetNbinsX(); i++){
		sumw += s_Kpipi_rw->GetBinContent(i);
		sumw2 += s_Kpipi_rw->GetBinError(i)*s_Kpipi_rw->GetBinError(i);
            
		//cout << "bin " << i << endl;
		//cout << " n = " <<  s_Kpipi_rw->GetBinContent(i) << endl;
		//cout << " e = " <<  s_Kpipi_rw->GetBinError(i) << endl;
		//cout << " neff = " <<  s_Kpipi_rw->GetBinContent(i)*s_Kpipi_rw->GetBinContent(i)/(s_Kpipi_rw->GetBinError(i)*s_Kpipi_rw->GetBinError(i)) << endl;
		}

	    //cout << "sumw = " << sumw << endl;
	    //cout << "weff = " << sumw/sumw2 << endl;
	    //cout << "Neff = " << sumw*sumw/sumw2 << endl;

            //SignalGenerator sg(pat,&fas);
            //sg.setWeighted();

	    // Need dummy file because of large event number
	    DiskResidentEventList eventListMC_rw(pat,("dummy_"+anythingToString(step)+".root").c_str(),"RECREATE");
/// Don't remove brackets, ensures memory is released
{
	    badEvents = 0;
	    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");

            for(int i = 0; i < eventListMC.size(); i++){
                //counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
                //DalitzEvent evt(evtPtr.get());
		DalitzEvent evt(eventListMC.getEvent(i));
		if(evt.sij(s234)/(GeV*GeV)> 4) continue;

		if(!(evt.phaseSpace() > 0.)){
		 	//cout << "evt " << i << " 0 phsp " << endl << evt << endl;
			badEvents++;
			continue;
		}
		if(TMath::IsNaN(fas.getVal(evt))){
		 	//cout << "evt " << i << " isNaN " << endl << evt << endl;
			badEvents++;
			continue;
		}

                double weight = fas.getVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
                //double weight = evt.phaseSpace()*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
		s_Kpipi_fit->Fill(evt.sij(s234)/(GeV*GeV),weight);
                s_Kpi_fit->Fill(evt.s(2,4)/(GeV*GeV),weight);
                s_pipi_fit->Fill(evt.s(3,4)/(GeV*GeV),weight);
                s_Dspipi_fit->Fill(evt.sij(s134)/(GeV*GeV),weight);
                s_DsK_fit->Fill(evt.s(1,2)/(GeV*GeV),weight);
                s_DsKpi_fit->Fill(evt.sij(s124)/(GeV*GeV),weight);
                s_Dspi_fit->Fill(evt.s(1,3)/(GeV*GeV),weight);
                s_Dspim_fit->Fill(evt.s(1,4)/(GeV*GeV),weight);

// 		weight = evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
// 		s_Kpipi_phsp->Fill(evt.sij(s234)/(GeV*GeV),weight);
//                 s_Kpi_phsp->Fill(evt.s(2,4)/(GeV*GeV),weight);
//                 s_pipi_phsp->Fill(evt.s(3,4)/(GeV*GeV),weight);
//                 s_Dspipi_phsp->Fill(evt.sij(s134)/(GeV*GeV),weight);
//                 s_DsK_phsp->Fill(evt.s(1,2)/(GeV*GeV),weight);
//                 s_DsKpi_phsp->Fill(evt.sij(s124)/(GeV*GeV),weight);
//                 s_Dspi_phsp->Fill(evt.s(1,3)/(GeV*GeV),weight);
//                 s_Dspim_phsp->Fill(evt.s(1,4)/(GeV*GeV),weight);
		evt.setWeight(weight);
		eventListMC_rw.Add(evt);
            }
	    cout << endl << "bad MINT events " << badEvents << " ( " << badEvents/(double)eventListMC.size() << " %)" << endl << endl;
}
/// Don't remove brackets, ensures memory is released
{  
	    badEvents = 0;
	    DiskResidentEventList eventListPhsp("/auto/data/dargent/BsDsKpipi/MINT/SignalIntegrationEvents_Phsp_15M_new.root","OPEN");
        
            for(int i = 0; i < eventListPhsp.size(); i++){                                
                double weight = 1.;//eventListPhsp[i].getWeight()/eventListPhsp[i].getGeneratorPdfRelativeToPhaseSpace();
		DalitzEvent evt(eventListPhsp.getEvent(i));
		if(evt.sij(s234)/(GeV*GeV)> 4) continue;

		if(!(evt.phaseSpace() > 0.)){
		 	//cout << "evt " << i << " 0 phsp " << endl << evt << endl;
			badEvents++;
			continue;
		}
		if(TMath::IsNaN(fas.getVal(evt))){
		 	//cout << "evt " << i << " isNaN " << endl << evt << endl;
			badEvents++;
			continue;
		}

		s_Kpipi_phsp->Fill(evt.sij(s234)/(GeV*GeV),weight);
                s_Kpi_phsp->Fill(evt.s(2,4)/(GeV*GeV),weight);
                s_pipi_phsp->Fill(evt.s(3,4)/(GeV*GeV),weight);
                s_Dspipi_phsp->Fill(evt.sij(s134)/(GeV*GeV),weight);
                s_DsK_phsp->Fill(evt.s(1,2)/(GeV*GeV),weight);
                s_DsKpi_phsp->Fill(evt.sij(s124)/(GeV*GeV),weight);
                s_Dspi_phsp->Fill(evt.s(1,3)/(GeV*GeV),weight);
                s_Dspim_phsp->Fill(evt.s(1,4)/(GeV*GeV),weight);
            }
	    cout << endl << "bad Phsp events " << badEvents << " ( " << badEvents/(double)eventListPhsp.size() << " %)" << endl << endl;
}

            TCanvas* c = new TCanvas();
                     
	    TF1* f = new TF1("f","1.0",0,30);
	    f->SetLineColor(kBlue);

  	    //gPad->SetLogy(1);
            s_Kpipi->SetMinimum(1);
            s_Kpipi->SetLineColor(kBlack);
            s_Kpipi->DrawNormalized("e",1);
            s_Kpipi_fit->SetMinimum(1);
            s_Kpipi_fit->SetLineColor(kBlue);
            s_Kpipi_fit->SetLineWidth(2);
            s_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi.eps").c_str());

  	    gPad->SetLogy(0);
            s_Kpipi_rw->SetMinimum(1);
            s_Kpipi_phsp->SetMinimum(1);
            s_Kpipi_rw->SetLineColor(kBlack);
            s_Kpipi_phsp->SetLineColor(kBlue);
            s_Kpipi_rw->DrawNormalized("e",1);
            s_Kpipi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_Kpipi.eps").c_str());

	    s_Kpipi_rw->Scale(1./s_Kpipi_rw->Integral());
    	    s_Kpipi_phsp->Scale(1./s_Kpipi_phsp->Integral());
   
	    for(int i = 1; i <= s_Kpipi_fit->GetNbinsX(); i++){
		double pull =  (s_Kpipi_rw->GetBinContent(i)-s_Kpipi_phsp->GetBinContent(i))/sqrt( pow(s_Kpipi_rw->GetBinError(i),2) + pow(s_Kpipi_phsp->GetBinError(i),2));
		s_Kpipi_pull->SetBinContent(i,pull);
		s_Kpipi_pull->SetBinError(i,1);
	    }
	    s_Kpipi_pull->Draw("e");
	    c->Print(((string)OutputDir+"pull_Kpipi.eps").c_str());
	    
	    s_Kpipi_rw->Divide(s_Kpipi_rw,s_Kpipi_phsp);
    	    s_Kpipi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_Kpipi_rw->SetMinimum(0.5);
	    s_Kpipi_rw->SetMaximum(1.5);
    	    s_Kpipi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Kpipi.eps").c_str());

 	    s_Kpipi_fit->Scale(1./s_Kpipi_fit->Integral());
     	    s_Kpipi->Scale(1./s_Kpipi->Integral());

	    for(int i = 1; i <= s_Kpipi_fit->GetNbinsX(); i++){
		double pull =  (s_Kpipi_fit->GetBinContent(i)-s_Kpipi->GetBinContent(i))/sqrt( pow(s_Kpipi_fit->GetBinError(i),2) + pow(s_Kpipi->GetBinError(i),2));
		s_Kpipi_pull2->SetBinContent(i,pull);
		s_Kpipi_pull2->SetBinError(i,1);
	    }
	    s_Kpipi_pull2->Draw("e");
	    c->Print(((string)OutputDir+"pull2_Kpipi.eps").c_str());

    	    s_Kpipi_fit->Divide(s_Kpipi_fit,s_Kpipi);
     	    s_Kpipi_fit->Draw("e");
     	    c->Print(((string)OutputDir+"eff2_Kpipi.eps").c_str());

// 	    gPad->SetLogy(1);
  	    s_Kpi->SetMinimum(1);
            s_Kpi->SetLineColor(kBlack);
            s_Kpi->DrawNormalized("e",1);
            s_Kpi_fit->SetMinimum(1);
            s_Kpi_fit->SetLineColor(kBlue);
            s_Kpi_fit->SetLineWidth(2);
            s_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi.eps").c_str());

  	    gPad->SetLogy(0);
            s_Kpi_rw->SetMinimum(1);
            s_Kpi_phsp->SetMinimum(1);
            s_Kpi_rw->SetLineColor(kBlack);
            s_Kpi_phsp->SetLineColor(kBlue);
            s_Kpi_rw->DrawNormalized("e",1);
            s_Kpi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_Kpi.eps").c_str());

	    s_Kpi_rw->Scale(1./s_Kpi_rw->Integral());
    	    s_Kpi_phsp->Scale(1./s_Kpi_phsp->Integral());
   	    s_Kpi_rw->Divide(s_Kpi_rw,s_Kpi_phsp);
    	    s_Kpi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_Kpi_rw->SetMinimum(0.5);
	    s_Kpi_rw->SetMaximum(1.5);
    	    s_Kpi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Kpi.eps").c_str());

//             gPad->SetLogy(1);
            s_pipi->SetMinimum(1);
            s_pipi->SetLineColor(kBlack);
            s_pipi->DrawNormalized("e",1);
            s_pipi_fit->SetMinimum(1);
            s_pipi_fit->SetLineColor(kBlue);
            s_pipi_fit->SetLineWidth(2);
            s_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi.eps").c_str());

  	    gPad->SetLogy(0);
            s_pipi_rw->SetMinimum(1);
            s_pipi_phsp->SetMinimum(1);
            s_pipi_rw->SetLineColor(kBlack);
            s_pipi_phsp->SetLineColor(kBlue);
            s_pipi_rw->DrawNormalized("e",1);
            s_pipi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_pipi.eps").c_str());

	    s_pipi_rw->Scale(1./s_pipi_rw->Integral());
    	    s_pipi_phsp->Scale(1./s_pipi_phsp->Integral());
   	    s_pipi_rw->Divide(s_pipi_rw,s_pipi_phsp);
    	    s_pipi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_pipi_rw->SetMinimum(0.5);
	    s_pipi_rw->SetMaximum(1.5);
    	    s_pipi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_pipi.eps").c_str());

//   	    gPad->SetLogy(1);
            s_Dspi->SetMinimum(1);
            s_Dspi->SetLineColor(kBlack);
            s_Dspi->DrawNormalized("e",1);
            s_Dspi_fit->SetMinimum(1);
            s_Dspi_fit->SetLineColor(kBlue);
            s_Dspi_fit->SetLineWidth(2);
            s_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi.eps").c_str());

  	    gPad->SetLogy(0);
            s_Dspi_rw->SetMinimum(1);
            s_Dspi_phsp->SetMinimum(1);
            s_Dspi_rw->SetLineColor(kBlack);
            s_Dspi_phsp->SetLineColor(kBlue);
            s_Dspi_rw->DrawNormalized("e",1);
            s_Dspi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_Dspi.eps").c_str());

	    s_Dspi_rw->Scale(1./s_Dspi_rw->Integral());
    	    s_Dspi_phsp->Scale(1./s_Dspi_phsp->Integral());
   	    s_Dspi_rw->Divide(s_Dspi_rw,s_Dspi_phsp);
    	    s_Dspi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_Dspi_rw->SetMinimum(0.5);
	    s_Dspi_rw->SetMaximum(1.5);
    	    s_Dspi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Dspi.eps").c_str());

//   	    gPad->SetLogy(1);
            s_Dspipi->SetMinimum(1);
            s_Dspipi->SetLineColor(kBlack);
            s_Dspipi->DrawNormalized("e",1);
            s_Dspipi_fit->SetMinimum(1);
            s_Dspipi_fit->SetLineColor(kBlue);
            s_Dspipi_fit->SetLineWidth(2);
            s_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi.eps").c_str());

  	    gPad->SetLogy(0);
            s_Dspipi_rw->SetMinimum(1);
            s_Dspipi_phsp->SetMinimum(1);
            s_Dspipi_rw->SetLineColor(kBlack);
            s_Dspipi_phsp->SetLineColor(kBlue);
            s_Dspipi_rw->DrawNormalized("e",1);
            s_Dspipi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_Dspipi.eps").c_str());

	    s_Dspipi_rw->Scale(1./s_Dspipi_rw->Integral());
    	    s_Dspipi_phsp->Scale(1./s_Dspipi_phsp->Integral());
   	    s_Dspipi_rw->Divide(s_Dspipi_rw,s_Dspipi_phsp);
    	    s_Dspipi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_Dspipi_rw->SetMinimum(0.5);
	    s_Dspipi_rw->SetMaximum(1.5);
    	    s_Dspipi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Dspipi.eps").c_str());

//   	    gPad->SetLogy(1);
            s_DsK->SetMinimum(1);
            s_DsK->SetLineColor(kBlack);
            s_DsK->DrawNormalized("e",1);
            s_DsK_fit->SetMinimum(1);
            s_DsK_fit->SetLineColor(kBlue);
            s_DsK_fit->SetLineWidth(2);
            s_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK.eps").c_str());

  	    gPad->SetLogy(0);
            s_DsK_rw->SetMinimum(1);
            s_DsK_phsp->SetMinimum(1);
            s_DsK_rw->SetLineColor(kBlack);
            s_DsK_phsp->SetLineColor(kBlue);
            s_DsK_rw->DrawNormalized("e",1);
            s_DsK_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_DsK.eps").c_str());

	    s_DsK_rw->Scale(1./s_DsK_rw->Integral());
    	    s_DsK_phsp->Scale(1./s_DsK_phsp->Integral());
   	    s_DsK_rw->Divide(s_DsK_rw,s_DsK_phsp);
    	    s_DsK_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_DsK_rw->SetMinimum(0.5);
	    s_DsK_rw->SetMaximum(1.5);
    	    s_DsK_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_DsK.eps").c_str());

//   	    gPad->SetLogy(1);
            s_DsKpi->SetMinimum(1);
            s_DsKpi->SetLineColor(kBlack);
            s_DsKpi->DrawNormalized("e",1);
            s_DsKpi_fit->SetMinimum(1);
            s_DsKpi_fit->SetLineColor(kBlue);
            s_DsKpi_fit->SetLineWidth(2);
            s_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi.eps").c_str());

  	    gPad->SetLogy(0);
            s_DsKpi_rw->SetMinimum(1);
            s_DsKpi_phsp->SetMinimum(1);
            s_DsKpi_rw->SetLineColor(kBlack);
            s_DsKpi_phsp->SetLineColor(kBlue);
            s_DsKpi_rw->DrawNormalized("e",1);
            s_DsKpi_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_DsKpi.eps").c_str());

	    s_DsKpi_rw->Scale(1./s_DsKpi_rw->Integral());
    	    s_DsKpi_phsp->Scale(1./s_DsKpi_phsp->Integral());
   	    s_DsKpi_rw->Divide(s_DsKpi_rw,s_DsKpi_phsp);
    	    s_DsKpi_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_DsKpi_rw->SetMinimum(0.5);
	    s_DsKpi_rw->SetMaximum(1.5);
    	    s_DsKpi_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_DsKpi.eps").c_str());

//   	    gPad->SetLogy(1);
            s_Dspim->SetMinimum(1);
            s_Dspim->SetLineColor(kBlack);
            s_Dspim->DrawNormalized("e",1);
            s_Dspim_fit->SetMinimum(1);
            s_Dspim_fit->SetLineColor(kBlue);
            s_Dspim_fit->SetLineWidth(2);
            s_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim.eps").c_str());

  	    gPad->SetLogy(0);
            s_Dspim_rw->SetMinimum(1);
            s_Dspim_phsp->SetMinimum(1);
            s_Dspim_rw->SetLineColor(kBlack);
            s_Dspim_phsp->SetLineColor(kBlue);
            s_Dspim_rw->DrawNormalized("e",1);
            s_Dspim_phsp->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"phsp_Dspim.eps").c_str());

	    s_Dspim_rw->Scale(1./s_Dspim_rw->Integral());
    	    s_Dspim_phsp->Scale(1./s_Dspim_phsp->Integral());
   	    s_Dspim_rw->Divide(s_Dspim_rw,s_Dspim_phsp);
    	    s_Dspim_rw->GetYaxis()->SetTitle("Reweighted EvtGen / PHSP");
	    s_Dspim_rw->SetMinimum(0.5);
	    s_Dspim_rw->SetMaximum(1.5);
    	    s_Dspim_rw->Draw("e");
	    f->Draw("same");
    	    c->Print(((string)OutputDir+"eff_Dspim.eps").c_str());

      	    getChi2(eventList,eventListMC_rw);        
        
    return 0;
}

void makeIntegratorFile(int step = 0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;

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
    

    TString outputName = (string)IntegratorEventFile;
    if(step>0) outputName.ReplaceAll(".root",("_" + anythingToString(step) + ".root").c_str());
    eventList.saveAsNtuple((string)outputName);
    return;
}

void makeIntegratorFilePhsp(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DiskResidentEventList eventList(pat,IntegratorEventFile,"RECREATE");

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
    
    eventList.save();
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
//     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_1.root");
//     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_2.root");
//    tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_3.root");
//     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_4.root");
//     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_5.root");
      tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_test.root");
//      tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_13266007.root");

    if (dbThis) cout << "Read the file" << endl;	
    double K_gen[5]; 
    double pip_gen[5]; 
    double pim_gen[5]; 
    double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];

    tree_gen->SetBranchStatus("*",0);
    tree_gen->SetBranchStatus("*TRUEP*",1);

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
 
/*      
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
*/
	TLorentzVector K_p;
	K_p.SetXYZM(K_gen[0],K_gen[1],K_gen[2],pdg[2].mass());
        TLorentzVector pip_p;
	pip_p.SetXYZM(pip_gen[0],pip_gen[1],pip_gen[2],pdg[3].mass());
	TLorentzVector pim_p;
	pim_p.SetXYZM(pim_gen[0],pim_gen[1],pim_gen[2],pdg[3].mass());
        TLorentzVector D_Kp_p;
	D_Kp_p.SetXYZM(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],pdg[2].mass());
        TLorentzVector D_Km_p;
	D_Km_p.SetXYZM(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],pdg[2].mass());
        TLorentzVector D_pim_p;
	D_pim_p.SetXYZM(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],pdg[3].mass());
	TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
	TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
        // array of vectors
	vector<TLorentzVector> vectorOfvectors; 

	//cout << std::fixed << std::setprecision(10) << K_p.M() << endl;
	//cout << std::fixed << std::setprecision(10) <<pip_p.M() << endl;
	//cout << std::fixed << std::setprecision(10) <<D_p.M() << endl;
	//cout << std::fixed << std::setprecision(10) <<B_p.M() << endl << endl;

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
    
    TString output = outputDir + "GenMC_EvtGen_test";
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
  TFile *FileMC =  TFile::Open("/auto/data/dargent/BsDsKpipi/MINT/GenMC_EvtGen.root");
  TTree* treeMC = dynamic_cast<TTree*>(FileMC->Get("DalitzEventList"));
  eventListMC.fromNtuple(treeMC,1);
  FileMC->Close();
            
  for(int i = 0; i < eventListMC.size(); i++){                                
		eventListMC[i].setGeneratorPdfRelativeToPhaseSpace(fas.getVal(eventListMC[i]));
  }
  eventListMC.save("/auto/data/dargent/BsDsKpipi/MINT/GenMC_EvtGen_rw.root");
}

int makeMINTtupleGenForToys(){
    
    string outputDir = "/auto/data/dargent/BsDsKpipi/MINT/";
	
    NamedParameter<int> EventPattern("Event Pattern", 521, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
    cout << " got event pattern: " << pdg << endl;

    // Read the momenta from ntuple
    TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_1.root");
     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_2.root");
    tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_3.root");
    tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_4.root");
    tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_5.root");
     tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_6.root");
//      tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/GenMC_13266007.root");

    double K_gen[5]; 
    double pip_gen[5]; 
    double pim_gen[5]; 
    double Ds_Kp_gen[5],Ds_Km_gen[5],Ds_pim_gen[5];

    tree_gen->SetBranchStatus("*",0);
    tree_gen->SetBranchStatus("*TRUEP*",1);

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
    int numEventsPerFile = 100000;
    int numFiles = numEvents/numEventsPerFile;
    cout << "I will produce " << numFiles << " files, each with " << numEventsPerFile << " events" << endl << endl;

    FitAmpIncoherentSum fas(pdg);
    fas.print();

    int counter = 0;
    for(int n= 0 ; n < numFiles; n++){

	DiskResidentEventList eventList(pdg,("/auto/data/dargent/BsDsKpipi/MINT/EvtGenToys/GenMC_" + anythingToString(n+1) + ".root").c_str(),"RECREATE");
	
	for(int i=0; i< numEventsPerFile; i++)
	{
		tree_gen->GetEntry(counter);
		
		TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
		TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
		TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
		TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
		TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
		TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
		TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
		TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
		vector<TLorentzVector> vectorOfvectors; 

		vectorOfvectors.push_back(B_p*MeV);      
		vectorOfvectors.push_back(D_p*MeV);
		vectorOfvectors.push_back(K_p*MeV); 
		vectorOfvectors.push_back(pip_p*MeV);
		vectorOfvectors.push_back(pim_p*MeV);
		
		DalitzEvent evt(pdg, vectorOfvectors);

		evt.setGeneratorPdfRelativeToPhaseSpace(fas.getVal(evt));

		eventList.Add(evt); 		
		if(counter < numEvents)counter++;
		else break;
	}    
	eventList.save();	
    }
        
    return 0;
}

void makeMINTtupleGenInAcc(){

    DalitzEventPattern pdg(531, -431, 321, 211, -211);
    vector<int> pdg_full_ids;
    pdg_full_ids.push_back(531); 
    pdg_full_ids.push_back(321); 
    pdg_full_ids.push_back(-321);
    pdg_full_ids.push_back(-211);
    pdg_full_ids.push_back(321); 
    pdg_full_ids.push_back(211);
    pdg_full_ids.push_back(-211);
    DalitzEventPattern pdg_full(pdg_full_ids);

    FitAmpIncoherentSum fas(pdg);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    TRandom3 ranLux;
    ranLux.SetSeed(0);
    gRandom = &ranLux;

    TFile* kinBs= new TFile("../AcceptancePhsp/LHCb8.root","Open");
    TH1D* h_Bs_pt = (TH1D*)kinBs->Get("pT");
    TH1D* h_Bs_eta = (TH1D*)kinBs->Get("eta");

    HyperHistogram hist_weights("../AcceptancePhsp/plots/weights.root");
    int dim = 1;

    TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
    tree_gen->Add("/auto/data/dargent/BsDsKpipi/EvtGen/Gen_6.root");
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

    DiskResidentEventList eventList(pdg,"SignalIntegrationEvents_EvtGen_inAcc.root","RECREATE");

    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);

        TLorentzVector K_p(K_gen[0],K_gen[1],K_gen[2],K_gen[3]);
        TLorentzVector pip_p(pip_gen[0],pip_gen[1],pip_gen[2],pip_gen[3]);
        TLorentzVector pim_p(pim_gen[0],pim_gen[1],pim_gen[2],pim_gen[3]);
        TLorentzVector D_Kp_p(Ds_Kp_gen[0],Ds_Kp_gen[1],Ds_Kp_gen[2],Ds_Kp_gen[3]);
        TLorentzVector D_Km_p(Ds_Km_gen[0],Ds_Km_gen[1],Ds_Km_gen[2],Ds_Km_gen[3]);
        TLorentzVector D_pim_p(Ds_pim_gen[0],Ds_pim_gen[1],Ds_pim_gen[2],Ds_pim_gen[3]);
        TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
        TLorentzVector B_p = K_p + pip_p + pim_p + D_p;

        // array of vectors
        vector<TLorentzVector> vectorOfvectors_full; 
        vectorOfvectors_full.push_back(B_p*MeV);      
        vectorOfvectors_full.push_back(D_Kp_p*MeV);
        vectorOfvectors_full.push_back(D_Km_p*MeV);
        vectorOfvectors_full.push_back(D_pim_p*MeV);
        vectorOfvectors_full.push_back(K_p*MeV); 
        vectorOfvectors_full.push_back(pip_p*MeV);
        vectorOfvectors_full.push_back(pim_p*MeV);
        DalitzEvent evt_full(pdg_full, vectorOfvectors_full);
         
        TLorentzVector mumsP4;
        mumsP4.SetPtEtaPhiM(h_Bs_pt->GetRandom()*1000.,h_Bs_eta->GetRandom(),ranLux.Uniform(0.,2*pi), pdg[0].mass());
        evt_full.setMothers3Momentum(mumsP4.Vect()); 
        
        double min_pt = evt_full.p(1).Pt(); 
        double max_pt = evt_full.p(1).Pt(); 
        for (int j = 1; j <= 6; j++) {
            if(evt_full.p(j).Pt()<min_pt)min_pt=evt_full.p(j).Pt();
            if(evt_full.p(j).Pt()>max_pt)max_pt=evt_full.p(j).Pt();
        }
         
        vector<TLorentzVector> vectorOfvectors; 
        vectorOfvectors.push_back(B_p*MeV);      
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(K_p*MeV); 
        vectorOfvectors.push_back(pip_p*MeV);
        vectorOfvectors.push_back(pim_p*MeV);
        DalitzEvent evt(pdg, vectorOfvectors);
      
	double w = 1;
        /// LHCb Acceptance
        for (int j = 1; j <= 6; j++) {            
            if (TMath::Abs(evt_full.p(j).Px()/evt_full.p(j).Pz()) > 0.3) w = 0;
            if (TMath::Abs(evt_full.p(j).Py()/evt_full.p(j).Pz()) > 0.25) w = 0;
            if (sqrt(pow(evt_full.p(j).Px()/evt_full.p(j).Pz(),2) + pow(evt_full.p(j).Py()/evt_full.p(j).Pz(),2)) <0.01) w = 0;
            if(evt_full.p(j).Pz() < 0.) w = 0;
        }
        
        /// Momentum cuts
        if( min_pt < 100  ) w = 0;
        if( max_pt < 1700  ) w = 0;
        
        int pt_counter = 0;
        if(evt_full.p(4).Pt() > 300)pt_counter++;
        if(evt_full.p(5).Pt() > 300)pt_counter++;
        if(evt_full.p(6).Pt() > 300)pt_counter++;
        if(pt_counter < 2) w = 0;
        
        if( min(evt_full.p(4).P(),min(evt_full.p(5).P(),evt_full.p(6).P())) < 2000 ) w = 0;
        if( evt_full.p(4).Pt()+evt_full.p(5).Pt()+evt_full.p(6).Pt() < 1250 ) w = 0;

        if( min(evt_full.p(1).P(),min(evt_full.p(2).P(),evt_full.p(3).P())) < 1000 ) w = 0;
        if( evt_full.p(1).Pt()+evt_full.p(2).Pt()+evt_full.p(3).Pt() < 1800 ) w = 0;

        /// PHSP cuts
        if(sqrt(evt.sij(s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) w = 0;
        //if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.) w = 0;            

	if(w==0)continue;

	HyperPoint point( dim );
        point.at(0)= 	log( min_pt  );
        /*
	if(K_gen[4] == min_pt) point.at(1)= 1/2. * log( (K_gen[3]+K_gen[2])/(K_gen[3]-K_gen[2]));
        else if(pip_gen[4] == min_pt) point.at(1)= 1/2. * log( (pip_gen[3]+pip_gen[2])/(pip_gen[3]-pip_gen[2]));
        else if(pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (pim_gen[3]+pim_gen[2])/(pim_gen[3]-pim_gen[2]));
        else if(Ds_Kp_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Kp_gen[3]+Ds_Kp_gen[2])/(Ds_Kp_gen[3]-Ds_Kp_gen[2]));
        else if(Ds_Km_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_Km_gen[3]+Ds_Km_gen[2])/(Ds_Km_gen[3]-Ds_Km_gen[2]));
        else if(Ds_pim_gen[4] == min_pt) point.at(1)= 1/2. * log( (Ds_pim_gen[3]+Ds_pim_gen[2])/(Ds_pim_gen[3]-Ds_pim_gen[2]));
	*/
        int bin = hist_weights.getBinning().getBinNum(point);
        if(hist_weights.checkBinNumber(bin)!= bin){
         	w = 0; //? should't happen
         	cout << "ERROR:: Event outside limits" << endl;	
        	cout << point.at(0) << endl;
        	cout << point.at(1) << endl << endl;
        }else w = hist_weights.getBinContent(bin);    
     
	evt.setWeight(w);
	evt.setGeneratorPdfRelativeToPhaseSpace(fas.getVal(evt));
	eventList.Add(evt);
    }

    eventList.save();
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
    gStyle->SetMarkerSize(0.8);
//      NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
//      if(! std::ifstream(((string)IntegratorEventFile).c_str()).good()) makeIntegratorFile();

// makeIntegratorFilePhsp();

//     HyperHistogram binningHist("histData.root",5);    
//  cout << endl << binningHist.integral();
// return 0;


   //    makeIntegratorFile(atoi(argv[1]));
 //    makeMINTtupleGen();
	makeMINTtupleGenInAcc();
  //     ampFit(atoi(argv[1]));
//      fracFit();
//     makeMINTtupleGenForToys();


    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
//
