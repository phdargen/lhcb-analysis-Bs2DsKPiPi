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
#include <TChain.h>
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

void AddScaledAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name, counted_ptr<IReturnComplex>& scale){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.multiply(scale);
    fas.addAsList(fas_2,1.);
}

void AddAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.addAsList(fas_2,1.);
}

class AmpRatio : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    int _f;
public:
    AmpRatio(FitParameter& re, FitParameter& im,int f = 1)
    : _re(re), _im(im), _f(f) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result(_re,static_cast<double>(_f) * _im); 
        return result;
    }
};

double cosThetaAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}

int ampFit(int step=0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);

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
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);
    NamedParameter<int>  doPlots("doPlots", 0);
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
	
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
            
    FitAmpSum fas_tmp(pat);
    FitAmpIncoherentSum fasBkg(pat);
    {
	fas_tmp.print();    
    	DalitzEventList eventListNorm;
  	TFile *file =  TFile::Open("SignalIntegrationEvents_toys_phspCut.root");
  	TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
  	eventListNorm.fromNtuple(tree,0.5);
  	fas_tmp.normalizeAmps(eventListNorm);
        fasBkg.normalizeAmps(eventListNorm);
    }
    
    counted_ptr<FitAmpList> List_1 = fas_tmp.GetCloneOfSubsetSameFitParameters("K(1)(1400)+");
    FitAmpSum fas(*List_1);

    /// Fix relative decay modes from charm input
    FitParameter a_K1_1270_re("a_K1_1270_Re",2,1,0.01);
    FitParameter a_K1_1270_im("a_K1_1270_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_K1_1270 = new AmpRatio(a_K1_1270_re,a_K1_1270_im);

    FitParameter a_K_1460_re("a_K_1460_Re",2,1,0.01);
    FitParameter a_K_1460_im("a_K_1460_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_K_1460 = new AmpRatio(a_K_1460_re,a_K_1460_im);

    FitParameter a_Ks_1410_re("a_Ks_1410_Re",2,1,0.01);
    FitParameter a_Ks_1410_im("a_Ks_1410_Im",2,0,0.01); 
    counted_ptr<IReturnComplex> a_Ks_1410 = new AmpRatio(a_Ks_1410_re,a_Ks_1410_im);

    AddScaledAmpsToList(fas_tmp, fas, "K(1)(1270)+",a_K1_1270);
    AddScaledAmpsToList(fas_tmp, fas, "K(1460)",a_K_1460);
    AddScaledAmpsToList(fas_tmp, fas, "K*(1410)+",a_Ks_1410);
//     AddAmpsToList(fas_tmp, fas, "K*(1680)+");
/*    AddAmpsToList(fas_tmp, fas, "NonResS0(->Ds-,pi+)");*/
//     AddAmpsToList(fas_tmp, fas, "NonResV0(->Ds-,pi+)");
//     AddAmpsToList(fas_tmp, fas, "NonResS0(->Ds-,K+)");
//     AddAmpsToList(fas_tmp, fas, "NonResV0(->Ds-,K+)");
    AddAmpsToList(fas_tmp, fas, "NonRes");

    FitParameter sigfraction("SigFraction",2,1.,0.01);
    DalitzEventList eventList;
    
    if(generateNew){
        SignalGenerator sg(pat,&fas);
        cout << "Generating " << Nevents << " MC events." << endl;
        sg.FillEventList(eventList, Nevents);
        eventList.saveAsNtuple(OutputRootFile);
    }
    
    if(!generateNew){
        /// Load data
        double t,dt;
        int f;
        double Bs_ID,Ds_ID;
        int q_OS;
        Short_t q_SS;
        double eta_OS;
        Float_t eta_SS;
        double sw;
        int year,run,Ds_finalState,trigger;

        double K[4];
        double pip[4];
        double pim[4];
        double Ds_Kp[4],Ds_Km[4],Ds_pim[4];
        double mB;
        
        TChain* tree;
	tree=new TChain("DecayTree");
	tree->Add(((string)InputDir+"Data/"+(string)channel+".root").c_str());
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("N_Bs_sw",1);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*DEC",1);
	tree->SetBranchStatus("*PROB",1);
	tree->SetBranchStatus("*OS",1);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("*ID*",1);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	tree->SetBranchStatus("BsDTF_*P*",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);

	tree->SetBranchAddress("Bs_DTF_TAU",&t);
	tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
	tree->SetBranchAddress("Ds_ID",&f);
	tree->SetBranchAddress("N_Bs_sw",&sw);
	tree->SetBranchAddress("year",&year);
	tree->SetBranchAddress("run",&run);
	tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
	tree->SetBranchAddress("TriggerCat",&trigger);
	tree->SetBranchAddress("Bs_DTF_MM",&mB);
	
	tree->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
	tree->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
	tree->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]);
	tree->SetBranchAddress("BsDTF_Kplus_PE",&K[3]);
	
	tree->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
	tree->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
	tree->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]);
	tree->SetBranchAddress("BsDTF_piplus_PE",&pip[3]);
	
	tree->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
	tree->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
	tree->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]);
	tree->SetBranchAddress("BsDTF_piminus_PE",&pim[3]);
	
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PX",&Ds_Kp[0]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PY",&Ds_Kp[1]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PZ",&Ds_Kp[2]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PE",&Ds_Kp[3]);
	
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PX",&Ds_Km[0]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PY",&Ds_Km[1]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PZ",&Ds_Km[2]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PE",&Ds_Km[3]);

	tree->SetBranchAddress("BsDTF_Ds_piminus_PX",&Ds_pim[0]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PY",&Ds_pim[1]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PZ",&Ds_pim[2]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PE",&Ds_pim[3]);

	int N_sample = tree->GetEntries();
        vector<int> b_indices;
        while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(ranLux.Uniform(0,N_sample)));
        sort(b_indices.begin(), b_indices.end());
        if(doBootstrap)N_sample = b_indices.size();

	TRandom3 rndm;
	int badEvents = 0;
	for(int i=0; i< N_sample; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << N_sample << endl;
		if(doBootstrap) tree->GetEntry(b_indices[i]);
		else tree->GetEntry(i);

		double sign = 1.;
		if(f > 0) sign = -1.;
		TLorentzVector K_p(sign*K[0],sign*K[1],sign*K[2],K[3]);
		TLorentzVector pip_p(sign*pip[0],sign*pip[1],sign*pip[2],pip[3]);
		TLorentzVector pim_p(sign*pim[0],sign*pim[1],sign*pim[2],pim[3]);
		TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
		TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
		TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
		TLorentzVector D_p = D_Kp_p + D_Km_p + D_pim_p;
		TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
		// array of vectors
		vector<TLorentzVector> vectorOfvectors;
	
		if((string)channel=="norm"){
			TLorentzVector pip1_p, pip2_p;
			if(rndm.Rndm()<0.5) {
			pip1_p = K_p;
			pip2_p = pip_p;
			}
			else {
			pip1_p = pip_p;
			pip2_p = K_p;
			}
			K_p = pip1_p;
			pip_p = pip2_p;
		}
	
		vectorOfvectors.push_back(B_p*MeV);
		vectorOfvectors.push_back(D_p*MeV);
		vectorOfvectors.push_back(K_p*MeV);
		vectorOfvectors.push_back(pip_p*MeV);
		vectorOfvectors.push_back(pim_p*MeV);
	
		DalitzEvent evt = DalitzEvent(pat, vectorOfvectors);
	
		if(!(evt.phaseSpace() > 0.)){
			badEvents++;
			continue;
		}

// 		if(sqrt(evt.sij(s234)/(GeV*GeV))<1) continue;
		
		evt.setWeight(sw);
		eventList.Add(evt);	
	}
	cout << endl << "bad events " << badEvents << " ( " << badEvents/(double) N_sample * 100. << " %)" << endl << endl;
     }
    
        DalitzHistoSet datH = eventList.weightedHistoSet();
        
        AmpsPdfFlexiFast* ampsSig;
	if(useLASSO) ampsSig = new AmpsPdfFlexiFast(pat, &fas_tmp, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
	else ampsSig = new AmpsPdfFlexiFast(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
        AmpsPdfFlexiFast ampsBkg(pat, &fasBkg, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
        DalitzSumPdf amps(sigfraction,*ampsSig,ampsBkg);
        
        Neg2LL neg2ll(amps, eventList);
        
        double stepSize = 1;
        lambda = lambda + (step-1) * stepSize;
        LASSO_flexi lasso(ampsSig,lambda);
        Neg2LLSum fcn(&neg2ll,&lasso);
        fcn.addConstraints(); 
        
        Minimiser mini;
        if(useLASSO)mini.attachFunction(&fcn);
        else mini.attachFunction(&neg2ll);
        mini.doFit();
        
        mini.printResultVsInput();
        
	if(useLASSO)ampsSig->doFinalStats(&mini);
	else {
		string outTableName = (string)OutputDir+"FitAmpResults.tex";
        	if(updateAnaNote)outTableName = "../../../../../TD-AnaNote/latex/tables/lassoFit/"+(string)OutputDir+"fitFractions.tex";
        	ampsSig->doFinalStatsAndSave(&mini,outTableName,((string)OutputDir+"FitAmpResults.root").c_str());
        }
                
        if(useLASSO){
            
	    double N_sig = 0;
            for (int i=0; i<eventList.size(); i++) N_sig += eventList[i].getWeight();


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
                y[0]=neg2ll.getVal() +  lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
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
        
        if(doPlots){
        
            cout << "Now plotting:" << endl;
            
//             DalitzHistoSet fitH = ampsSig.histoSet();
//             datH.drawWithFitNorm(fitH, ((string)OutputDir+(string)"datFit_l_"+anythingToString(step)+"_").c_str(),"eps");
//             std::vector<DalitzHistoSet> EachAmpsHistos = ampsSig.GetEachAmpsHistograms();
//             datH.drawWithFitAndEachAmps(datH, fitH, EachAmpsHistos, ((string)OutputDir+(string)"WithAmps").c_str(), "eps");

            int nBins = 50;
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

	    TH1D* h_cosTheta_Kpi= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);
            
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

		h_cosTheta_Kpi->Fill(cosThetaAngle(eventList[i],2,4,1,3),eventList[i].getWeight());
		h_cosTheta_Dspi->Fill(cosThetaAngle(eventList[i],1,3,2,4),eventList[i].getWeight());
		h_phi_Kpi_Dspi->Fill(acoplanarityAngle(eventList[i],2,4,1,3),eventList[i].getWeight());
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
	    TH1D* h_cosTheta_Kpi_fit= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

	    TH1D* m_Kpipi_fit_1 = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit_1 = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit_1 = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit_1 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit_1 = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit_1 = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit_1 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit_1 = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);	
	    TH1D* h_cosTheta_Kpi_fit_1= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit_1= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit_1= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

	    TH1D* m_Kpipi_fit_2 = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit_2 = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit_2 = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit_2 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit_2 = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit_2 = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit_2 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit_2 = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);	
	    TH1D* h_cosTheta_Kpi_fit_2= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit_2= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit_2= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

	    TH1D* m_Kpipi_fit_3 = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit_3 = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit_3 = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit_3 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit_3 = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit_3 = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit_3 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit_3 = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);	
	    TH1D* h_cosTheta_Kpi_fit_3= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit_3= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit_3= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

	    TH1D* m_Kpipi_fit_4 = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit_4 = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit_4 = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit_4 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit_4 = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit_4 = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit_4 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit_4 = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);	
	    TH1D* h_cosTheta_Kpi_fit_4= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit_4= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit_4= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);

	    TH1D* m_Kpipi_fit_5 = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1,2);
	    TH1D* m_Kpi_fit_5 = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.6,1.3);
	    TH1D* m_pipi_fit_5 = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,0.2,1.3);
	    TH1D* m_Dspipi_fit_5 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_DsK_fit_5 = new TH1D("",";#left[m(D_{s}^{-} K^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5.)  ;
	    TH1D* m_DsKpi_fit_5 = new TH1D("",";#left[m(D_{s}^{-} K^{+} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,2,5.5);
	    TH1D* m_Dspi_fit_5 = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);
	    TH1D* m_Dspim_fit_5 = new TH1D("",";#left[m(D_{s}^{-} #pi^{-})#right] (GeV/c^{2});Events (norm.) ",nBins,1.5,5);	
	    TH1D* h_cosTheta_Kpi_fit_5= new TH1D("",";cos #theta_{K^{+}#pi^{-}}; Events (norm.) ",40,-1,1);
	    TH1D* h_cosTheta_Dspi_fit_5= new TH1D("",";cos #theta_{D_{s}#pi^{+}}; Events (norm.) ",40,0,1);
    	    TH1D* h_phi_Kpi_Dspi_fit_5= new TH1D("",";#phi_{K^{+}#pi^{-},D_{s}#pi^{+}}; Events (norm.)",40,-3.141,3.141);
	
            //SignalGenerator sg(pat,&fas);
            //sg.setWeighted();

	    DalitzEventList eventListMC;
	    TFile *FileMC =  TFile::Open(((string) IntegratorEventFile).c_str());
	    TTree* treeMC = dynamic_cast<TTree*>(FileMC->Get("DalitzEventList"));
	    eventListMC.fromNtuple(treeMC,1);
	    FileMC->Close();

	    vector<string> ampNames1;
	    ampNames1.push_back("K(1)(1270)+");

	    vector<string> ampNames2;
	    ampNames2.push_back("K(1)(1400)+");

	    vector<string> ampNames3;
	    ampNames3.push_back("K*(1410)+");

	    vector<string> ampNames4;
	    ampNames4.push_back("NonRes");

	    vector<string> ampNames5;
	    ampNames5.push_back("K(1460)+");
		           
            for(int i = 0; i < eventListMC.size(); i++){
                                
                //counted_ptr<IDalitzEvent> evtPtr(sg.newEvent());
                //DalitzEvent evt(evtPtr.get());
                double weight = ampsSig->un_normalised_noPs(eventListMC[i])*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
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
		h_cosTheta_Kpi_fit->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight);
		h_cosTheta_Dspi_fit->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight);
		h_phi_Kpi_Dspi_fit->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight);
 
		double weight1 = fas.getAmpSqr(eventListMC[i], ampNames1) 
			* eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
		double weight2 = fas.getAmpSqr(eventListMC[i], ampNames2) 
			* eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
 		double weight3 = fas.getAmpSqr(eventListMC[i], ampNames3) 
			* eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
 		double weight4 = fas.getAmpSqr(eventListMC[i], ampNames4) 
			* eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
 		double weight5 = fas.getAmpSqr(eventListMC[i], ampNames5) 
			* eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
  
		m_Kpipi_fit_1->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight1);
	        m_Kpi_fit_1->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight1);
	        m_pipi_fit_1->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight1);
	        m_Dspipi_fit_1->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight1);
	        m_DsK_fit_1->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight1);
	        m_DsKpi_fit_1->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight1);
	        m_Dspi_fit_1->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight1);
	        m_Dspim_fit_1->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight1);
		h_cosTheta_Kpi_fit_1->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight1);
		h_cosTheta_Dspi_fit_1->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight1);
		h_phi_Kpi_Dspi_fit_1->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight1);

   	        m_Kpipi_fit_2->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight2);
	        m_Kpi_fit_2->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight2);
	        m_pipi_fit_2->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight2);
	        m_Dspipi_fit_2->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight2);
	        m_DsK_fit_2->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight2);
	        m_DsKpi_fit_2->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight2);
	        m_Dspi_fit_2->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight2);
	        m_Dspim_fit_2->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight2);
		h_cosTheta_Kpi_fit_2->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight2);
		h_cosTheta_Dspi_fit_2->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight2);
		h_phi_Kpi_Dspi_fit_2->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight2);

   	        m_Kpipi_fit_3->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight3);
	        m_Kpi_fit_3->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight3);
	        m_pipi_fit_3->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight3);
	        m_Dspipi_fit_3->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight3);
	        m_DsK_fit_3->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight3);
	        m_DsKpi_fit_3->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight3);
	        m_Dspi_fit_3->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight3);
	        m_Dspim_fit_3->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight3);
		h_cosTheta_Kpi_fit_3->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight3);
		h_cosTheta_Dspi_fit_3->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight3);
		h_phi_Kpi_Dspi_fit_3->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight3);

   	        m_Kpipi_fit_4->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight4);
	        m_Kpi_fit_4->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight4);
	        m_pipi_fit_4->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight4);
	        m_Dspipi_fit_4->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight4);
	        m_DsK_fit_4->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight4);
	        m_DsKpi_fit_4->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight4);
	        m_Dspi_fit_4->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight4);
	        m_Dspim_fit_4->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight4);
		h_cosTheta_Kpi_fit_4->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight4);
		h_cosTheta_Dspi_fit_4->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight4);
		h_phi_Kpi_Dspi_fit_4->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight4);

   	        m_Kpipi_fit_5->Fill(sqrt(eventListMC[i].sij(s234)/(GeV*GeV)),weight5);
	        m_Kpi_fit_5->Fill(sqrt(eventListMC[i].s(2,4)/(GeV*GeV)),weight5);
	        m_pipi_fit_5->Fill(sqrt(eventListMC[i].s(3,4)/(GeV*GeV)),weight5);
	        m_Dspipi_fit_5->Fill(sqrt(eventListMC[i].sij(s134)/(GeV*GeV)),weight5);
	        m_DsK_fit_5->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight5);
	        m_DsKpi_fit_5->Fill(sqrt(eventListMC[i].sij(s124)/(GeV*GeV)),weight5);
	        m_Dspi_fit_5->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight5);
	        m_Dspim_fit_5->Fill(sqrt(eventListMC[i].s(1,4)/(GeV*GeV)),weight5);
		h_cosTheta_Kpi_fit_5->Fill(cosThetaAngle(eventListMC[i],2,4,1,3),weight5);
		h_cosTheta_Dspi_fit_5->Fill(cosThetaAngle(eventListMC[i],1,3,2,4),weight5);
		h_phi_Kpi_Dspi_fit_5->Fill(acoplanarityAngle(eventListMC[i],2,4,1,3),weight5);
		
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

	    h_cosTheta_Kpi->SetMinimum(0);
            h_cosTheta_Kpi->SetLineColor(kBlack);
            h_cosTheta_Kpi->DrawNormalized("e1",1);
            h_cosTheta_Kpi_fit->SetLineColor(kBlue);
            h_cosTheta_Kpi_fit->SetLineWidth(3);
            h_cosTheta_Kpi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_cosTheta_Kpi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Kpi.pdf").c_str());

	    h_cosTheta_Dspi->SetMinimum(0);
            h_cosTheta_Dspi->SetLineColor(kBlack);
            h_cosTheta_Dspi->DrawNormalized("e1",1);
            h_cosTheta_Dspi_fit->SetLineColor(kBlue);
            h_cosTheta_Dspi_fit->SetLineWidth(3);
            h_cosTheta_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_cosTheta_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Dspi.pdf").c_str());

	    h_phi_Kpi_Dspi->SetMinimum(0);
            h_phi_Kpi_Dspi->SetLineColor(kBlack);
            h_phi_Kpi_Dspi->DrawNormalized("e1",1);
            h_phi_Kpi_Dspi_fit->SetLineColor(kBlue);
            h_phi_Kpi_Dspi_fit->SetLineWidth(3);
            h_phi_Kpi_Dspi_fit->DrawNormalized("histcsame",1);
            c->Print(((string)OutputDir+"h_phi_Kpi_Dspi.eps").c_str());
	    if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_phi_Kpi_Dspi.pdf").c_str());


            m_Kpipi->SetMinimum(0.01);
            m_Kpipi->SetLineColor(kBlack);
            m_Kpipi->DrawNormalized("e1",1);
            m_Kpipi_fit->SetLineColor(kBlue);
            m_Kpipi_fit->SetLineWidth(3);
            m_Kpipi_fit->DrawNormalized("histcsame",1);
            m_Kpipi_fit_1->SetLineColor(kRed+1);
            m_Kpipi_fit_1->SetLineWidth(2);
            m_Kpipi_fit_1->SetFillColor(kRed+1);
            m_Kpipi_fit_1->SetFillStyle(3353);
            m_Kpipi_fit_1->DrawNormalized("histcsame",m_Kpipi_fit_1->Integral()/m_Kpipi_fit->Integral());
            m_Kpipi_fit_2->SetLineColor(kGreen+3);
            m_Kpipi_fit_2->SetLineWidth(2);
            m_Kpipi_fit_2->SetFillColor(kGreen+3);
            m_Kpipi_fit_2->SetFillStyle(3353);
            m_Kpipi_fit_2->DrawNormalized("histcsame",m_Kpipi_fit_2->Integral()/m_Kpipi_fit->Integral());
            m_Kpipi_fit_3->SetLineColor(kMagenta+3);
            m_Kpipi_fit_3->SetLineWidth(2);
            m_Kpipi_fit_3->SetFillColor(kMagenta+3);
            m_Kpipi_fit_3->SetFillStyle(3353);
            m_Kpipi_fit_3->DrawNormalized("histcsame",m_Kpipi_fit_3->Integral()/m_Kpipi_fit->Integral());
            m_Kpipi_fit_4->SetLineColor(kBlack);
            m_Kpipi_fit_4->SetLineWidth(3);
            m_Kpipi_fit_4->SetLineStyle(kDashed);
            m_Kpipi_fit_4->DrawNormalized("histcsame",m_Kpipi_fit_4->Integral()/m_Kpipi_fit->Integral());
            m_Kpipi_fit_5->SetLineColor(kGray+3);
            m_Kpipi_fit_5->SetLineWidth(2);
            m_Kpipi_fit_5->SetFillColor(kGray+3);
            m_Kpipi_fit_5->SetFillStyle(1001);
            m_Kpipi_fit_5->DrawNormalized("histcsame",m_Kpipi_fit_5->Integral()/m_Kpipi_fit->Integral());
            c->Print(((string)OutputDir+"m_Kpipi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpipi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Kpipi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpipi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_Kpi->SetMinimum(0.01);
            m_Kpi->SetLineColor(kBlack);
            m_Kpi->DrawNormalized("e1",1);
            m_Kpi_fit->SetLineColor(kBlue);
            m_Kpi_fit->SetLineWidth(3);
            m_Kpi_fit->DrawNormalized("histcsame",1);
            m_Kpi_fit_1->SetLineColor(kRed+1);
            m_Kpi_fit_1->SetLineWidth(2);
            m_Kpi_fit_1->SetFillColor(kRed+1);
            m_Kpi_fit_1->SetFillStyle(3353);
            m_Kpi_fit_1->DrawNormalized("histcsame",m_Kpi_fit_1->Integral()/m_Kpi_fit->Integral());
            m_Kpi_fit_2->SetLineColor(kGreen+3);
            m_Kpi_fit_2->SetLineWidth(2);
            m_Kpi_fit_2->SetFillColor(kGreen+3);
            m_Kpi_fit_2->SetFillStyle(3353);
            m_Kpi_fit_2->DrawNormalized("histcsame",m_Kpi_fit_2->Integral()/m_Kpi_fit->Integral());
            m_Kpi_fit_3->SetLineColor(kMagenta+3);
            m_Kpi_fit_3->SetLineWidth(2);
            m_Kpi_fit_3->SetFillColor(kMagenta+3);
            m_Kpi_fit_3->SetFillStyle(3353);
            m_Kpi_fit_3->DrawNormalized("histcsame",m_Kpi_fit_3->Integral()/m_Kpi_fit->Integral());
            m_Kpi_fit_4->SetLineColor(kBlack);
            m_Kpi_fit_4->SetLineWidth(3);
            m_Kpi_fit_4->SetLineStyle(kDashed);
            m_Kpi_fit_4->DrawNormalized("histcsame",m_Kpi_fit_4->Integral()/m_Kpi_fit->Integral());
            m_Kpi_fit_5->SetLineColor(kGray+3);
            m_Kpi_fit_5->SetLineWidth(2);
            m_Kpi_fit_5->SetFillColor(kGray+3);
            m_Kpi_fit_5->SetFillStyle(1001);
            m_Kpi_fit_5->DrawNormalized("histcsame",m_Kpi_fit_5->Integral()/m_Kpi_fit->Integral());
            c->Print(((string)OutputDir+"m_Kpi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Kpi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Kpi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_pipi->SetMinimum(0.01);
            m_pipi->SetLineColor(kBlack);
            m_pipi->DrawNormalized("e1",1);
            m_pipi_fit->SetLineColor(kBlue);
            m_pipi_fit->SetLineWidth(3);
            m_pipi_fit->DrawNormalized("histcsame",1);
            m_pipi_fit_1->SetLineColor(kRed+1);
            m_pipi_fit_1->SetLineWidth(2);
            m_pipi_fit_1->SetFillColor(kRed+1);
            m_pipi_fit_1->SetFillStyle(3353);
            m_pipi_fit_1->DrawNormalized("histcsame",m_pipi_fit_1->Integral()/m_pipi_fit->Integral());
            m_pipi_fit_2->SetLineColor(kGreen+3);
            m_pipi_fit_2->SetLineWidth(2);
            m_pipi_fit_2->SetFillColor(kGreen+3);
            m_pipi_fit_2->SetFillStyle(3353);
            m_pipi_fit_2->DrawNormalized("histcsame",m_pipi_fit_2->Integral()/m_pipi_fit->Integral());
            m_pipi_fit_3->SetLineColor(kMagenta+3);
            m_pipi_fit_3->SetLineWidth(2);
            m_pipi_fit_3->SetFillColor(kMagenta+3);
            m_pipi_fit_3->SetFillStyle(3353);
            m_pipi_fit_3->DrawNormalized("histcsame",m_pipi_fit_3->Integral()/m_pipi_fit->Integral());
            m_pipi_fit_4->SetLineColor(kBlack);
            m_pipi_fit_4->SetLineWidth(3);
            m_pipi_fit_4->SetLineStyle(kDashed);
            m_pipi_fit_4->DrawNormalized("histcsame",m_pipi_fit_4->Integral()/m_pipi_fit->Integral());
            m_pipi_fit_5->SetLineColor(kGray+3);
            m_pipi_fit_5->SetLineWidth(2);
            m_pipi_fit_5->SetFillColor(kGray+3);
            m_pipi_fit_5->SetFillStyle(1001);
            m_pipi_fit_5->DrawNormalized("histcsame",m_pipi_fit_5->Integral()/m_pipi_fit->Integral());
            c->Print(((string)OutputDir+"m_pipi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_pipi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_pipi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_pipi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_Dspipi->SetMinimum(0.01);
            m_Dspipi->SetLineColor(kBlack);
            m_Dspipi->DrawNormalized("e1",1);
            m_Dspipi_fit->SetLineColor(kBlue);
            m_Dspipi_fit->SetLineWidth(3);
            m_Dspipi_fit->DrawNormalized("histcsame",1);
            m_Dspipi_fit_1->SetLineColor(kRed+1);
            m_Dspipi_fit_1->SetLineWidth(2);
            m_Dspipi_fit_1->SetFillColor(kRed+1);
            m_Dspipi_fit_1->SetFillStyle(3353);
            m_Dspipi_fit_1->DrawNormalized("histcsame",m_Dspipi_fit_1->Integral()/m_Dspipi_fit->Integral());
            m_Dspipi_fit_2->SetLineColor(kGreen+3);
            m_Dspipi_fit_2->SetLineWidth(2);
            m_Dspipi_fit_2->SetFillColor(kGreen+3);
            m_Dspipi_fit_2->SetFillStyle(3353);
            m_Dspipi_fit_2->DrawNormalized("histcsame",m_Dspipi_fit_2->Integral()/m_Dspipi_fit->Integral());
            m_Dspipi_fit_3->SetLineColor(kMagenta+3);
            m_Dspipi_fit_3->SetLineWidth(2);
            m_Dspipi_fit_3->SetFillColor(kMagenta+3);
            m_Dspipi_fit_3->SetFillStyle(3353);
            m_Dspipi_fit_3->DrawNormalized("histcsame",m_Dspipi_fit_3->Integral()/m_Dspipi_fit->Integral());
            m_Dspipi_fit_4->SetLineColor(kBlack);
            m_Dspipi_fit_4->SetLineWidth(3);
            m_Dspipi_fit_4->SetLineStyle(kDashed);
            m_Dspipi_fit_4->DrawNormalized("histcsame",m_Dspipi_fit_4->Integral()/m_Dspipi_fit->Integral());
            m_Dspipi_fit_5->SetLineColor(kGray+3);
            m_Dspipi_fit_5->SetLineWidth(2);
            m_Dspipi_fit_5->SetFillColor(kGray+3);
            m_Dspipi_fit_5->SetFillStyle(1001);
            m_Dspipi_fit_5->DrawNormalized("histcsame",m_Dspipi_fit_5->Integral()/m_Dspipi_fit->Integral());
            c->Print(((string)OutputDir+"m_Dspipi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspipi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Dspipi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspipi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_Dspi->SetMinimum(0.01);
            m_Dspi->SetLineColor(kBlack);
            m_Dspi->DrawNormalized("e1",1);
            m_Dspi_fit->SetLineColor(kBlue);
            m_Dspi_fit->SetLineWidth(3);
            m_Dspi_fit->DrawNormalized("histcsame",1);
            m_Dspi_fit_1->SetLineColor(kRed+1);
            m_Dspi_fit_1->SetLineWidth(2);
            m_Dspi_fit_1->SetFillColor(kRed+1);
            m_Dspi_fit_1->SetFillStyle(3353);
            m_Dspi_fit_1->DrawNormalized("histcsame",m_Dspi_fit_1->Integral()/m_Dspi_fit->Integral());
            m_Dspi_fit_2->SetLineColor(kGreen+3);
            m_Dspi_fit_2->SetLineWidth(2);
            m_Dspi_fit_2->SetFillColor(kGreen+3);
            m_Dspi_fit_2->SetFillStyle(3353);
            m_Dspi_fit_2->DrawNormalized("histcsame",m_Dspi_fit_2->Integral()/m_Dspi_fit->Integral());
            m_Dspi_fit_3->SetLineColor(kMagenta+3);
            m_Dspi_fit_3->SetLineWidth(2);
            m_Dspi_fit_3->SetFillColor(kMagenta+3);
            m_Dspi_fit_3->SetFillStyle(3353);
            m_Dspi_fit_3->DrawNormalized("histcsame",m_Dspi_fit_3->Integral()/m_Dspi_fit->Integral());
            m_Dspi_fit_4->SetLineColor(kBlack);
            m_Dspi_fit_4->SetLineWidth(3);
            m_Dspi_fit_4->SetLineStyle(kDashed);
            m_Dspi_fit_4->DrawNormalized("histcsame",m_Dspi_fit_4->Integral()/m_Dspi_fit->Integral());
            m_Dspi_fit_5->SetLineColor(kGray+3);
            m_Dspi_fit_5->SetLineWidth(2);
            m_Dspi_fit_5->SetFillColor(kGray+3);
            m_Dspi_fit_5->SetFillStyle(1001);
            m_Dspi_fit_5->DrawNormalized("histcsame",m_Dspi_fit_5->Integral()/m_Dspi_fit->Integral());
            c->Print(((string)OutputDir+"m_Dspi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Dspi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_DsK->SetMinimum(0.01);
            m_DsK->SetLineColor(kBlack);
            m_DsK->DrawNormalized("e1",1);
            m_DsK_fit->SetLineColor(kBlue);
            m_DsK_fit->SetLineWidth(3);
            m_DsK_fit->DrawNormalized("histcsame",1);
            m_DsK_fit_1->SetLineColor(kRed+1);
            m_DsK_fit_1->SetLineWidth(2);
            m_DsK_fit_1->SetFillColor(kRed+1);
            m_DsK_fit_1->SetFillStyle(3353);
            m_DsK_fit_1->DrawNormalized("histcsame",m_DsK_fit_1->Integral()/m_DsK_fit->Integral());
            m_DsK_fit_2->SetLineColor(kGreen+3);
            m_DsK_fit_2->SetLineWidth(2);
            m_DsK_fit_2->SetFillColor(kGreen+3);
            m_DsK_fit_2->SetFillStyle(3353);
            m_DsK_fit_2->DrawNormalized("histcsame",m_DsK_fit_2->Integral()/m_DsK_fit->Integral());
            m_DsK_fit_3->SetLineColor(kMagenta+3);
            m_DsK_fit_3->SetLineWidth(2);
            m_DsK_fit_3->SetFillColor(kMagenta+3);
            m_DsK_fit_3->SetFillStyle(3353);
            m_DsK_fit_3->DrawNormalized("histcsame",m_DsK_fit_3->Integral()/m_DsK_fit->Integral());
            m_DsK_fit_4->SetLineColor(kBlack);
            m_DsK_fit_4->SetLineWidth(3);
            m_DsK_fit_4->SetLineStyle(kDashed);
            m_DsK_fit_4->DrawNormalized("histcsame",m_DsK_fit_4->Integral()/m_DsK_fit->Integral());
            m_DsK_fit_5->SetLineColor(kGray+3);
            m_DsK_fit_5->SetLineWidth(2);
            m_DsK_fit_5->SetFillColor(kGray+3);
            m_DsK_fit_5->SetFillStyle(1001);
            m_DsK_fit_5->DrawNormalized("histcsame",m_DsK_fit_5->Integral()/m_DsK_fit->Integral());
            c->Print(((string)OutputDir+"m_DsK_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsK_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_DsK_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsK_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_DsKpi->SetMinimum(0.01);
            m_DsKpi->SetLineColor(kBlack);
            m_DsKpi->DrawNormalized("e1",1);
            m_DsKpi_fit->SetLineColor(kBlue);
            m_DsKpi_fit->SetLineWidth(3);
            m_DsKpi_fit->DrawNormalized("histcsame",1);
            m_DsKpi_fit_1->SetLineColor(kRed+1);
            m_DsKpi_fit_1->SetLineWidth(2);
            m_DsKpi_fit_1->SetFillColor(kRed+1);
            m_DsKpi_fit_1->SetFillStyle(3353);
            m_DsKpi_fit_1->DrawNormalized("histcsame",m_DsKpi_fit_1->Integral()/m_DsKpi_fit->Integral());
            m_DsKpi_fit_2->SetLineColor(kGreen+3);
            m_DsKpi_fit_2->SetLineWidth(2);
            m_DsKpi_fit_2->SetFillColor(kGreen+3);
            m_DsKpi_fit_2->SetFillStyle(3353);
            m_DsKpi_fit_2->DrawNormalized("histcsame",m_DsKpi_fit_2->Integral()/m_DsKpi_fit->Integral());
            m_DsKpi_fit_3->SetLineColor(kMagenta+3);
            m_DsKpi_fit_3->SetLineWidth(2);
            m_DsKpi_fit_3->SetFillColor(kMagenta+3);
            m_DsKpi_fit_3->SetFillStyle(3353);
            m_DsKpi_fit_3->DrawNormalized("histcsame",m_DsKpi_fit_3->Integral()/m_DsKpi_fit->Integral());
            m_DsKpi_fit_4->SetLineColor(kBlack);
            m_DsKpi_fit_4->SetLineWidth(3);
            m_DsKpi_fit_4->SetLineStyle(kDashed);
            m_DsKpi_fit_4->DrawNormalized("histcsame",m_DsKpi_fit_4->Integral()/m_DsKpi_fit->Integral());
            m_DsKpi_fit_5->SetLineColor(kGray+3);
            m_DsKpi_fit_5->SetLineWidth(2);
            m_DsKpi_fit_5->SetFillColor(kGray+3);
            m_DsKpi_fit_5->SetFillStyle(1001);
            m_DsKpi_fit_5->DrawNormalized("histcsame",m_DsKpi_fit_5->Integral()/m_DsKpi_fit->Integral());
            c->Print(((string)OutputDir+"m_DsKpi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsKpi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_DsKpi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_DsKpi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            m_Dspim->SetMinimum(0.01);
            m_Dspim->SetLineColor(kBlack);
            m_Dspim->DrawNormalized("e1",1);
            m_Dspim_fit->SetLineColor(kBlue);
            m_Dspim_fit->SetLineWidth(3);
            m_Dspim_fit->DrawNormalized("histcsame",1);
            m_Dspim_fit_1->SetLineColor(kRed+1);
            m_Dspim_fit_1->SetLineWidth(2);
            m_Dspim_fit_1->SetFillColor(kRed+1);
            m_Dspim_fit_1->SetFillStyle(3353);
            m_Dspim_fit_1->DrawNormalized("histcsame",m_Dspim_fit_1->Integral()/m_Dspim_fit->Integral());
            m_Dspim_fit_2->SetLineColor(kGreen+3);
            m_Dspim_fit_2->SetLineWidth(2);
            m_Dspim_fit_2->SetFillColor(kGreen+3);
            m_Dspim_fit_2->SetFillStyle(3353);
            m_Dspim_fit_2->DrawNormalized("histcsame",m_Dspim_fit_2->Integral()/m_Dspim_fit->Integral());
            m_Dspim_fit_3->SetLineColor(kMagenta+3);
            m_Dspim_fit_3->SetLineWidth(2);
            m_Dspim_fit_3->SetFillColor(kMagenta+3);
            m_Dspim_fit_3->SetFillStyle(3353);
            m_Dspim_fit_3->DrawNormalized("histcsame",m_Dspim_fit_3->Integral()/m_Dspim_fit->Integral());
            m_Dspim_fit_4->SetLineColor(kBlack);
            m_Dspim_fit_4->SetLineWidth(3);
            m_Dspim_fit_4->SetLineStyle(kDashed);
            m_Dspim_fit_4->DrawNormalized("histcsame",m_Dspim_fit_4->Integral()/m_Dspim_fit->Integral());
            m_Dspim_fit_5->SetLineColor(kGray+3);
            m_Dspim_fit_5->SetLineWidth(2);
            m_Dspim_fit_5->SetFillColor(kGray+3);
            m_Dspim_fit_5->SetFillStyle(1001);
            m_Dspim_fit_5->DrawNormalized("histcsame",m_Dspim_fit_5->Integral()/m_Dspim_fit->Integral());
            c->Print(((string)OutputDir+"m_Dspim_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspim_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"m_Dspim_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"m_Dspim_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            h_cosTheta_Kpi->SetMinimum(0.01);
            h_cosTheta_Kpi->SetLineColor(kBlack);
            h_cosTheta_Kpi->DrawNormalized("e1",1);
            h_cosTheta_Kpi_fit->SetLineColor(kBlue);
            h_cosTheta_Kpi_fit->SetLineWidth(3);
            h_cosTheta_Kpi_fit->DrawNormalized("histcsame",1);
            h_cosTheta_Kpi_fit_1->SetLineColor(kRed+1);
            h_cosTheta_Kpi_fit_1->SetLineWidth(2);
            h_cosTheta_Kpi_fit_1->SetFillColor(kRed+1);
            h_cosTheta_Kpi_fit_1->SetFillStyle(3353);
            h_cosTheta_Kpi_fit_1->DrawNormalized("histcsame",h_cosTheta_Kpi_fit_1->Integral()/h_cosTheta_Kpi_fit->Integral());
            h_cosTheta_Kpi_fit_2->SetLineColor(kGreen+3);
            h_cosTheta_Kpi_fit_2->SetLineWidth(2);
            h_cosTheta_Kpi_fit_2->SetFillColor(kGreen+3);
            h_cosTheta_Kpi_fit_2->SetFillStyle(3353);
            h_cosTheta_Kpi_fit_2->DrawNormalized("histcsame",h_cosTheta_Kpi_fit_2->Integral()/h_cosTheta_Kpi_fit->Integral());
            h_cosTheta_Kpi_fit_3->SetLineColor(kMagenta+3);
            h_cosTheta_Kpi_fit_3->SetLineWidth(2);
            h_cosTheta_Kpi_fit_3->SetFillColor(kMagenta+3);
            h_cosTheta_Kpi_fit_3->SetFillStyle(3353);
            h_cosTheta_Kpi_fit_3->DrawNormalized("histcsame",h_cosTheta_Kpi_fit_3->Integral()/h_cosTheta_Kpi_fit->Integral());
            h_cosTheta_Kpi_fit_4->SetLineColor(kBlack);
            h_cosTheta_Kpi_fit_4->SetLineWidth(3);
            h_cosTheta_Kpi_fit_4->SetLineStyle(kDashed);
            h_cosTheta_Kpi_fit_4->DrawNormalized("histcsame",h_cosTheta_Kpi_fit_4->Integral()/h_cosTheta_Kpi_fit->Integral());
            h_cosTheta_Kpi_fit_5->SetLineColor(kGray+3);
            h_cosTheta_Kpi_fit_5->SetLineWidth(2);
            h_cosTheta_Kpi_fit_5->SetFillColor(kGray+3);
            h_cosTheta_Kpi_fit_5->SetFillStyle(1001);
            h_cosTheta_Kpi_fit_5->DrawNormalized("histcsame",h_cosTheta_Kpi_fit_5->Integral()/h_cosTheta_Kpi_fit->Integral());
            c->Print(((string)OutputDir+"h_cosTheta_Kpi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Kpi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"h_cosTheta_Kpi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Kpi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            h_cosTheta_Dspi->SetMinimum(0.01);
            h_cosTheta_Dspi->SetLineColor(kBlack);
            h_cosTheta_Dspi->DrawNormalized("e1",1);
            h_cosTheta_Dspi_fit->SetLineColor(kBlue);
            h_cosTheta_Dspi_fit->SetLineWidth(3);
            h_cosTheta_Dspi_fit->DrawNormalized("histcsame",1);
            h_cosTheta_Dspi_fit_1->SetLineColor(kRed+1);
            h_cosTheta_Dspi_fit_1->SetLineWidth(2);
            h_cosTheta_Dspi_fit_1->SetFillColor(kRed+1);
            h_cosTheta_Dspi_fit_1->SetFillStyle(3353);
            h_cosTheta_Dspi_fit_1->DrawNormalized("histcsame",h_cosTheta_Dspi_fit_1->Integral()/h_cosTheta_Dspi_fit->Integral());
            h_cosTheta_Dspi_fit_2->SetLineColor(kGreen+3);
            h_cosTheta_Dspi_fit_2->SetLineWidth(2);
            h_cosTheta_Dspi_fit_2->SetFillColor(kGreen+3);
            h_cosTheta_Dspi_fit_2->SetFillStyle(3353);
            h_cosTheta_Dspi_fit_2->DrawNormalized("histcsame",h_cosTheta_Dspi_fit_2->Integral()/h_cosTheta_Dspi_fit->Integral());
            h_cosTheta_Dspi_fit_3->SetLineColor(kMagenta+3);
            h_cosTheta_Dspi_fit_3->SetLineWidth(2);
            h_cosTheta_Dspi_fit_3->SetFillColor(kMagenta+3);
            h_cosTheta_Dspi_fit_3->SetFillStyle(3353);
            h_cosTheta_Dspi_fit_3->DrawNormalized("histcsame",h_cosTheta_Dspi_fit_3->Integral()/h_cosTheta_Dspi_fit->Integral());
            h_cosTheta_Dspi_fit_4->SetLineColor(kBlack);
            h_cosTheta_Dspi_fit_4->SetLineWidth(3);
            h_cosTheta_Dspi_fit_4->SetLineStyle(kDashed);
            h_cosTheta_Dspi_fit_4->DrawNormalized("histcsame",h_cosTheta_Dspi_fit_4->Integral()/h_cosTheta_Dspi_fit->Integral());
            h_cosTheta_Dspi_fit_5->SetLineColor(kGray+3);
            h_cosTheta_Dspi_fit_5->SetLineWidth(2);
            h_cosTheta_Dspi_fit_5->SetFillColor(kGray+3);
            h_cosTheta_Dspi_fit_5->SetFillStyle(1001);
            h_cosTheta_Dspi_fit_5->DrawNormalized("histcsame",h_cosTheta_Dspi_fit_5->Integral()/h_cosTheta_Dspi_fit->Integral());
            c->Print(((string)OutputDir+"h_cosTheta_Dspi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Dspi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"h_cosTheta_Dspi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_cosTheta_Dspi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

            h_phi_Kpi_Dspi->SetMinimum(0.01);
            h_phi_Kpi_Dspi->SetLineColor(kBlack);
            h_phi_Kpi_Dspi->DrawNormalized("e1",1);
            h_phi_Kpi_Dspi_fit->SetLineColor(kBlue);
            h_phi_Kpi_Dspi_fit->SetLineWidth(3);
            h_phi_Kpi_Dspi_fit->DrawNormalized("histcsame",1);
            h_phi_Kpi_Dspi_fit_1->SetLineColor(kRed+1);
            h_phi_Kpi_Dspi_fit_1->SetLineWidth(2);
            h_phi_Kpi_Dspi_fit_1->SetFillColor(kRed+1);
            h_phi_Kpi_Dspi_fit_1->SetFillStyle(3353);
            h_phi_Kpi_Dspi_fit_1->DrawNormalized("histcsame",h_phi_Kpi_Dspi_fit_1->Integral()/h_phi_Kpi_Dspi_fit->Integral());
            h_phi_Kpi_Dspi_fit_2->SetLineColor(kGreen+3);
            h_phi_Kpi_Dspi_fit_2->SetLineWidth(2);
            h_phi_Kpi_Dspi_fit_2->SetFillColor(kGreen+3);
            h_phi_Kpi_Dspi_fit_2->SetFillStyle(3353);
            h_phi_Kpi_Dspi_fit_2->DrawNormalized("histcsame",h_phi_Kpi_Dspi_fit_2->Integral()/h_phi_Kpi_Dspi_fit->Integral());
            h_phi_Kpi_Dspi_fit_3->SetLineColor(kMagenta+3);
            h_phi_Kpi_Dspi_fit_3->SetLineWidth(2);
            h_phi_Kpi_Dspi_fit_3->SetFillColor(kMagenta+3);
            h_phi_Kpi_Dspi_fit_3->SetFillStyle(3353);
            h_phi_Kpi_Dspi_fit_3->DrawNormalized("histcsame",h_phi_Kpi_Dspi_fit_3->Integral()/h_phi_Kpi_Dspi_fit->Integral());
            h_phi_Kpi_Dspi_fit_4->SetLineColor(kBlack);
            h_phi_Kpi_Dspi_fit_4->SetLineWidth(3);
            h_phi_Kpi_Dspi_fit_4->SetLineStyle(kDashed);
            h_phi_Kpi_Dspi_fit_4->DrawNormalized("histcsame",h_phi_Kpi_Dspi_fit_4->Integral()/h_phi_Kpi_Dspi_fit->Integral());
            h_phi_Kpi_Dspi_fit_5->SetLineColor(kGray+3);
            h_phi_Kpi_Dspi_fit_5->SetLineWidth(2);
            h_phi_Kpi_Dspi_fit_5->SetFillColor(kGray+3);
            h_phi_Kpi_Dspi_fit_5->SetFillStyle(1001);
            h_phi_Kpi_Dspi_fit_5->DrawNormalized("histcsame",h_phi_Kpi_Dspi_fit_5->Integral()/h_phi_Kpi_Dspi_fit->Integral());
            c->Print(((string)OutputDir+"h_phi_Kpi_Dspi_mod.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_phi_Kpi_Dspi_mod.pdf").c_str());
	    gPad->SetLogy(1);
            c->Print(((string)OutputDir+"h_phi_Kpi_Dspi_mod_log.eps").c_str());
            if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/lassoFit/"+(string)OutputDir +"h_phi_Kpi_Dspi_mod_log.pdf").c_str());
            gPad->SetLogy(0);

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
