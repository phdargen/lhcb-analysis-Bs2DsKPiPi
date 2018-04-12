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
#include "TFile.h"
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
        return ampSq;// * getEvent()->phaseSpace();
    }
    
    std::complex<double> ComplexVal_un_normalised_noPs(IDalitzEvent& evt){
        return  _amps->ComplexVal(evt);
    }
    
    AmpsPdfFlexiFast(const DalitzEventPattern& pat
                     , IFastAmplitudeIntegrable* amps
                     , MinuitParameterSet* mps
                     , double precision=1.e-4
                     , std::string method="efficient"
                     , std::string fname =  "SignalIntegrationEvents.root", bool generateNew = false, bool genMoreEvents = false
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
        //bool generateNew = ((string)_integratorSource == (string)"new");
        if(nonFlat){
            cout << "AmpsPdfFlexiFast uses nonFlat integration." << endl;
            if(generateNew){
                _sgGen =  new SignalGenerator(pat);
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
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 25);       
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
    
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
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
    
    FitAmpSum fas(pat);
    fas.getVal(eventNorm[0]);
    fas.print();
        
    FitAmpIncoherentSum fasBkg(pat);
    fasBkg.getVal(eventNorm[0]);

    {
        DalitzEventList eventNorm2;
        eventNorm2.generatePhaseSpaceEvents(200000,pat); 
        fas.normalizeAmps(eventNorm2);
        fasBkg.normalizeAmps(eventNorm2);
    }
    
    FitParameter sigfraction("SigFraction",2,1.,0.01);
    
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
        
        AmpsPdfFlexiFast ampsSig(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
        AmpsPdfFlexiFast ampsBkg(pat, &fasBkg, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
        DalitzSumPdf amps(sigfraction,ampsSig,ampsBkg);
        
        Neg2LL neg2ll(amps, eventList);
        
        double stepSize = 1;
        lambda = lambda + (step-1) * stepSize;
        LASSO_flexi lasso(&ampsSig,lambda);
        Neg2LLSum fcn(&neg2ll,&lasso);
        
        Minimiser mini;
        if(useLASSO)mini.attachFunction(&fcn);
        else mini.attachFunction(&neg2ll);
        mini.doFit();
        
        mini.printResultVsInput();
        ampsSig.doFinalStats(&mini);
                
        if(useLASSO){
            
            TFile* out_LASSO = new TFile(((string)OutputDir+"LASSO_"+anythingToString(step)+".root").c_str(),"RECREATE");
            Double_t x[1],y[1];
            vector<double> thresholds;
            thresholds.push_back(0.001);
            thresholds.push_back(0.002);
            thresholds.push_back(0.003);
            thresholds.push_back(0.004);
            thresholds.push_back(0.005);
            thresholds.push_back(0.006);
            thresholds.push_back(0.007);
            thresholds.push_back(0.008);
            thresholds.push_back(0.009);
            thresholds.push_back(0.01);
            thresholds.push_back(0.02);
            thresholds.push_back(0.05);
            
            x[0]=lambda;
            for(int i = 0; i < thresholds.size() ; i++){
                y[0]=neg2ll.getVal() + 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
                TGraph* aic = new TGraph(1,x,y);
                aic->SetName( ("AIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
                aic->SetTitle("");
                aic->GetXaxis()->SetTitle("#lambda");
                aic->GetXaxis()->SetTitleOffset(0.65);
                aic->GetYaxis()->SetTitle("AIC");
                aic->Draw("A*");
                aic->Write();
            }
            
            for(int i = 0; i < thresholds.size() ; i++){
                y[0]=neg2ll.getVal() + lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(eventList.size()*sigfraction);
                TGraph* bic = new TGraph(1,x,y);
                bic->SetName( ("BIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
                bic->SetTitle("");
                bic->GetXaxis()->SetTitle("#lambda");
                bic->GetXaxis()->SetTitleOffset(0.65);
                bic->GetYaxis()->SetTitle("BIC");
                bic->Draw("A*");
                bic->Write();
            }
            
            for(int i = 0; i < thresholds.size() ; i++){
                y[0]=lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
                TGraph* r = new TGraph(1,x,y);
                r->SetName( ("r_"+anythingToString((int) (thresholds[i]*1000))).c_str());
                r->SetTitle("");
                r->GetXaxis()->SetTitle("#lambda");
                r->GetXaxis()->SetTitleOffset(0.65);
                r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
                r->Draw("A*");
                r->Write();
            }
            
            y[0]=neg2ll.getVal() ;
            TGraph* nll = new TGraph(1,x,y);
            nll->SetName("Neg2LL");
            nll->SetTitle("");
            nll->GetXaxis()->SetTitle("#lambda");
            nll->GetXaxis()->SetTitleOffset(0.65);
            nll->GetYaxis()->SetTitle("-2 logL");
            nll->Draw("A*");
            nll->Write();
            
            y[0]= lasso.sumOfFitFractions() ;
            TGraph* ff = new TGraph(1,x,y);
            ff->SetName("SumOfFitFractions");
            ff->SetTitle("");
            ff->GetXaxis()->SetTitle("#lambda");
            ff->GetXaxis()->SetTitleOffset(0.65);
            ff->GetYaxis()->SetTitle("Total Fit Fraction");
            ff->Draw("A*");
            ff->Write();
            
            y[0]= lasso.absSumOfInterferenceFractions() ;
            TGraph* iff = new TGraph(1,x,y);
            iff->SetName("AbsSumOfInterferenceFractions");
            iff->SetTitle("");
            iff->GetXaxis()->SetTitle("#lambda");
            iff->GetXaxis()->SetTitleOffset(0.65);
            iff->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
            iff->Draw("A*");
            iff->Write();
            
            /*
            y[0]= getChi2( eventList, eventListMC) ;
            out_LASSO->cd();
            TGraph* chi2 = new TGraph(1,x,y);
            chi2->SetName("Chi2");
            chi2->SetTitle("");
            chi2->GetXaxis()->SetTitle("#lambda");
            chi2->GetXaxis()->SetTitleOffset(0.65);
            chi2->GetYaxis()->SetTitle("#Chi^{2}/Bin");
            chi2->Draw("A*");
            chi2->Write();
            
            for (int i = 0; i < eventListMC.size(); i++){
                //DalitzEvent evt = eventListMC[i];
                double ampVal=amps.getVal_noPs(eventListMC[i]);
                eventListMC[i].setWeight(ampVal*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace());
                //weightedMC.Add(evt);
            }
            //weightedMC.save(((string)OutputDir+"mcEvents.root").c_str());
            //DalitzHistoSet mcH = eventListMC.weightedHistoSet();
            //datH.drawWithFitNorm(mcH, ((string)OutputDir+"mcFit_l_"+anythingToString(step)+"_").c_str(),"eps");
            
            y[0]= getChi2_pionSwitched( eventList, eventListMC) ;
            out_LASSO->cd();
            TGraph* chi2_pionSwitched = new TGraph(1,x,y);
            chi2_pionSwitched->SetName("Chi2_pionSwitched");
            chi2_pionSwitched->SetTitle("");
            chi2_pionSwitched->GetXaxis()->SetTitle("#lambda");
            chi2_pionSwitched->GetXaxis()->SetTitleOffset(0.65);
            chi2_pionSwitched->GetYaxis()->SetTitle("#Chi^{2}/Bin (pion switched)");
            chi2_pionSwitched->Draw("A*");
            chi2_pionSwitched->Write();
            */
            
            out_LASSO->Write();
            out_LASSO->Close(); 
        }
        
        else{
        
            cout << "Now plotting:" << endl;
            
            DalitzHistoSet fitH = ampsSig.histoSet();
            datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_l_"+anythingToString(step)+"_").c_str(),"eps");
            std::vector<DalitzHistoSet> EachAmpsHistos = ampsSig.GetEachAmpsHistograms();
            datH.drawWithFitAndEachAmps(datH, fitH, EachAmpsHistos, ((string)OutputDir+(string)"WithAmps").c_str(), "eps");
            
            int nBins = 50;
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
            
	    TH1D* s_Kpipi = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,1,4);
	    TH1D* s_Kpi = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.3,1.6);
	    TH1D* s_pipi = new TH1D("",";#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,1.6);
	    TH1D* s_Dspipi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
	    TH1D* s_DsK = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30)  ;
	    TH1D* s_DsKpi = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
	    TH1D* s_Dspi = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
            TH1D* s_Dspim = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
		
 	    TH1D* m_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
            
            TH2D* s_Kpi_pipi = new TH2D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});#left[m^{2}(#pi^{+} #pi^{-})#right] (GeV^{2}/c^{4}) ",60,0.2,1.6,60,0.,1.6);
            TH2D* s_DsKpi_Dspi = new TH2D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4}); #left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4}) ",80,5,30,80,0,25);
            TH2D* s_DsK_Dspi = new TH2D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4}); #left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4}) ",60,0,30,60,0,25);

            for (int i=0; i<eventList.size(); i++) {
                s_Kpipi->Fill(eventList[i].sij(s234)/(GeV*GeV),eventList[i].getWeight());
                s_Kpi->Fill(eventList[i].s(2,4)/(GeV*GeV),eventList[i].getWeight());
                s_pipi->Fill(eventList[i].s(3,4)/(GeV*GeV),eventList[i].getWeight());
                s_Dspipi->Fill(eventList[i].sij(s134)/(GeV*GeV),eventList[i].getWeight());
                s_DsK->Fill(eventList[i].s(1,2)/(GeV*GeV),eventList[i].getWeight());
                s_DsKpi->Fill(eventList[i].sij(s124)/(GeV*GeV),eventList[i].getWeight());
                s_Dspi->Fill(eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());
                s_Dspim->Fill(eventList[i].s(1,4)/(GeV*GeV),eventList[i].getWeight());

                m_Kpipi->Fill(sqrt(eventList[i].sij(s234)/(GeV*GeV)),eventList[i].getWeight());
                m_Kpi->Fill(sqrt(eventList[i].s(2,4)/(GeV*GeV)),eventList[i].getWeight());
                m_pipi->Fill(sqrt(eventList[i].s(3,4)/(GeV*GeV)),eventList[i].getWeight());
                m_Dspipi->Fill(sqrt(eventList[i].sij(s134)/(GeV*GeV)),eventList[i].getWeight());
                m_DsK->Fill(sqrt(eventList[i].s(1,2)/(GeV*GeV)),eventList[i].getWeight());
                m_DsKpi->Fill(sqrt(eventList[i].sij(s124)/(GeV*GeV)),eventList[i].getWeight());
                m_Dspi->Fill(sqrt(eventList[i].s(1,3)/(GeV*GeV)),eventList[i].getWeight());
                m_Dspim->Fill(sqrt(eventList[i].s(1,4)/(GeV*GeV)),eventList[i].getWeight());

                s_Kpi_pipi->Fill(eventList[i].s(2,4)/(GeV*GeV),eventList[i].s(3,4)/(GeV*GeV),eventList[i].getWeight());
                s_DsKpi_Dspi->Fill(eventList[i].sij(s124)/(GeV*GeV),eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());
                s_DsK_Dspi->Fill(eventList[i].s(1,2)/(GeV*GeV),eventList[i].s(1,3)/(GeV*GeV),eventList[i].getWeight());

            }    
            
	    TH1D* s_Kpipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,1,4);
	    TH1D* s_Kpi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0.3,1.6);
	    TH1D* s_pipi_fit = new TH1D("",";#left[m^{2}(K^{+} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,1.6);
	    TH1D* s_Dspipi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
	    TH1D* s_DsK_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,30);
	    TH1D* s_DsKpi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} K^{+} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,5,30);
	    TH1D* s_Dspi_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{+})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
	    TH1D* s_Dspim_fit = new TH1D("",";#left[m^{2}(D_{s}^{-} #pi^{-})#right] (GeV^{2}/c^{4});Events (norm.) ",nBins,0,25);
	
	    TH1D* m_Kpipi_fit = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);		
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
                double weight = ampsSig.RealVal(eventListMC[i])*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                s_Kpipi_fit->Fill(eventListMC[i].sij(s234)/(GeV*GeV),weight);
                s_Kpi_fit->Fill(eventListMC[i].s(2,4)/(GeV*GeV),weight);
                s_pipi_fit->Fill(eventListMC[i].s(3,4)/(GeV*GeV),weight);
                s_Dspipi_fit->Fill(eventListMC[i].sij(s134)/(GeV*GeV),weight);
                s_DsK_fit->Fill(eventListMC[i].s(1,2)/(GeV*GeV),weight);
                s_DsKpi_fit->Fill(eventListMC[i].sij(s124)/(GeV*GeV),weight);
                s_Dspi_fit->Fill(eventListMC[i].s(1,3)/(GeV*GeV),weight);
                s_Dspim_fit->Fill(eventListMC[i].s(1,4)/(GeV*GeV),weight);

   	        m_Kpipi_fit->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight);
	        m_Kpi_fit->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight);
	        m_pipi_fit->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight);
	        m_Dspipi_fit->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight);
	        m_DsK_fit->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight);
	        m_DsKpi_fit->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight);
	        m_Dspi_fit->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight);
	        m_Dspim_fit->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight);
		
		eventListMC[i].setWeight(weight);
            }
            
            TCanvas* c = new TCanvas();
                     
            s_Kpi_pipi->SetMinimum(0);
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

            s_Kpipi->SetMinimum(0);
            s_Kpipi->SetLineColor(kBlack);
            s_Kpipi->DrawNormalized("e1",1);
            s_Kpipi_fit->SetLineColor(kBlue);
            s_Kpipi_fit->SetLineWidth(3);
            s_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_Kpipi.pdf").c_str());

            s_Kpi->SetMinimum(0);
            s_Kpi->SetLineColor(kBlack);
            s_Kpi->DrawNormalized("e1",1);
            s_Kpi_fit->SetLineColor(kBlue);
            s_Kpi_fit->SetLineWidth(3);
            s_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Kpi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_Kpi.pdf").c_str());

	    s_pipi->SetMinimum(0);            
            s_pipi->SetLineColor(kBlack);
            s_pipi->DrawNormalized("e1",1);
            s_pipi_fit->SetLineColor(kBlue);
            s_pipi_fit->SetLineWidth(3);
            s_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_pipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_pipi.pdf").c_str());

	    s_Dspipi->SetMinimum(0);            
            s_Dspipi->SetLineColor(kBlack);
            s_Dspipi->DrawNormalized("e1",1);
            s_Dspipi_fit->SetLineColor(kBlue);
            s_Dspipi_fit->SetLineWidth(3);
            s_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_Dspipi.pdf").c_str());

	    s_DsK->SetMinimum(0);
            s_DsK->SetLineColor(kBlack);
            s_DsK->DrawNormalized("e1",1);
            s_DsK_fit->SetLineColor(kBlue);
            s_DsK_fit->SetLineWidth(3);
            s_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsK.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_DsK.pdf").c_str());
	    
	    s_DsKpi->SetMinimum(0);            
            s_DsKpi->SetLineColor(kBlack);
            s_DsKpi->DrawNormalized("e1",1);
            s_DsKpi_fit->SetLineColor(kBlue);
            s_DsKpi_fit->SetLineWidth(3);
            s_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_DsKpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_DsKpi.pdf").c_str());
	    
	    s_Dspi->SetMinimum(0);
            s_Dspi->SetLineColor(kBlack);
            s_Dspi->DrawNormalized("e1",1);
            s_Dspi_fit->SetLineColor(kBlue);
            s_Dspi_fit->SetLineWidth(3);
            s_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_Dspi.pdf").c_str());
	    
	    s_Dspim->SetMinimum(0);
            s_Dspim->SetLineColor(kBlack);
            s_Dspim->DrawNormalized("e1",1);
            s_Dspim_fit->SetLineColor(kBlue);
            s_Dspim_fit->SetLineWidth(3);
            s_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"s_Dspim.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"s_Dspim.pdf").c_str());
	    

            m_Kpipi->SetMinimum(0);
            m_Kpipi->SetLineColor(kBlack);
            m_Kpipi->DrawNormalized("e1",1);
            m_Kpipi_fit->SetLineColor(kBlue);
            m_Kpipi_fit->SetLineWidth(3);
            m_Kpipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Kpipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpipi.pdf").c_str());

            m_Kpi->SetMinimum(0);
            m_Kpi->SetLineColor(kBlack);
            m_Kpi->DrawNormalized("e1",1);
            m_Kpi_fit->SetLineColor(kBlue);
            m_Kpi_fit->SetLineWidth(3);
            m_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Kpi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpi.pdf").c_str());

	    m_pipi->SetMinimum(0);            
            m_pipi->SetLineColor(kBlack);
            m_pipi->DrawNormalized("e1",1);
            m_pipi_fit->SetLineColor(kBlue);
            m_pipi_fit->SetLineWidth(3);
            m_pipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_pipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_pipi.pdf").c_str());

	    m_Dspipi->SetMinimum(0);            
            m_Dspipi->SetLineColor(kBlack);
            m_Dspipi->DrawNormalized("e1",1);
            m_Dspipi_fit->SetLineColor(kBlue);
            m_Dspipi_fit->SetLineWidth(3);
            m_Dspipi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspipi.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspipi.pdf").c_str());

	    m_DsK->SetMinimum(0);
            m_DsK->SetLineColor(kBlack);
            m_DsK->DrawNormalized("e1",1);
            m_DsK_fit->SetLineColor(kBlue);
            m_DsK_fit->SetLineWidth(3);
            m_DsK_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_DsK.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsK.pdf").c_str());
	    
	    m_DsKpi->SetMinimum(0);            
            m_DsKpi->SetLineColor(kBlack);
            m_DsKpi->DrawNormalized("e1",1);
            m_DsKpi_fit->SetLineColor(kBlue);
            m_DsKpi_fit->SetLineWidth(3);
            m_DsKpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_DsKpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsKpi.pdf").c_str());
	    
	    m_Dspi->SetMinimum(0);
            m_Dspi->SetLineColor(kBlack);
            m_Dspi->DrawNormalized("e1",1);
            m_Dspi_fit->SetLineColor(kBlue);
            m_Dspi_fit->SetLineWidth(3);
            m_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspi.pdf").c_str());
	    
	    m_Dspim->SetMinimum(0);
            m_Dspim->SetLineColor(kBlack);
            m_Dspim->DrawNormalized("e1",1);
            m_Dspim_fit->SetLineColor(kBlue);
            m_Dspim_fit->SetLineWidth(3);
            m_Dspim_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"m_Dspim.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspim.pdf").c_str());

	    getChi2(eventList,eventListMC);

        }
        
        
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
    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.print();
    fas.getVal(eventListPhsp[0]);
    fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);
    
    sg.FillEventList(eventList, IntegratorEvents);
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    for(int i = 0; i < eventList.size(); i++){
	if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 1.95 && sqrt(eventList[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventList[i].s(3,4)/(GeV*GeV)) < 1.2)eventList_cut.Add(eventList[i]);
    }

    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    eventList_cut.saveAsNtuple(IntegratorEventFile);
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
	if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2){ 
		if(evt.phaseSpace() > 0.){
			eventList.Add(evt);
			i++;
		}
	}
    }

    cout << "Generated " << eventList.size() << " events inside selected phasespace region" << endl;
    
    eventList.saveAsNtuple(IntegratorEventFile);
    return;
}

void makeIntegratorFileEvtGen(){
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventList_cut;
    DiskResidentEventList eventList(pat,"/auto/data/dargent/BsDsKpipi/MINT/SignalIntegrationEvents_EvtGen_new.root","OPEN");

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    for(int i = 0; i < IntegratorEvents; i++){

        DalitzEvent evt = eventList.getEvent(i);

	if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2)eventList_cut.Add(evt);
    }

    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    eventList_cut.saveAsNtuple("SignalIntegrationEvents_MINT_PhspCut.root");
    return;
}

int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

//makeIntegratorFileEvtGen(); return 0;

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    if(! std::ifstream(((string)IntegratorEventFile).c_str()).good()) makeIntegratorFilePhsp();
  
    ampFit(atoi(argv[1]));
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
//
