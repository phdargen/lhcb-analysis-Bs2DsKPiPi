#ifndef MySPlot_HH
#define MySPlot_HH

class RooAbsReal;
class RooAbsPdf;
class RooFitResult;
class RooRealVar;
class RooSimultaneous;
// class TNamed;

#include "RooArgList.h"
#include "RooDataSet.h"


#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include <TNamed.h>

//using namespace RooStats;
using namespace std;
using namespace RooFit ;

//namespace RooStats{

class MySPlot: public TNamed {

  public:

    ~MySPlot();
    MySPlot();
    MySPlot(const MySPlot &other);
    MySPlot(const char* name, const char* title);
    MySPlot(const char* name, const char* title, const RooDataSet &data);
    MySPlot(const char* name, const char* title,RooDataSet& data, RooAbsPdf* pdf,
     const RooArgList &yieldsList,const RooArgSet &projDeps=RooArgSet(),
     bool useWeights=kTRUE, bool copyDataSet = kFALSE, const char* newName = "",
     const RooCmdArg& fitToarg5=RooCmdArg::none(),
     const RooCmdArg& fitToarg6=RooCmdArg::none(),
     const RooCmdArg& fitToarg7=RooCmdArg::none(),
     const RooCmdArg& fitToarg8=RooCmdArg::none());
     MySPlot(const char* name, const char* title,RooDataSet& data, RooAbsPdf* pdf,
     const RooArgList &allYieldsList, const RooArgList &fixedYields, const RooArgSet &projDeps=RooArgSet(),
     bool useWeights=kTRUE, bool copyDataSet = kFALSE, const char* newName = "",
     const RooCmdArg& fitToarg5=RooCmdArg::none(),
     const RooCmdArg& fitToarg6=RooCmdArg::none(),
     const RooCmdArg& fitToarg7=RooCmdArg::none(),
     const RooCmdArg& fitToarg8=RooCmdArg::none());

    RooDataSet* SetSData(RooDataSet* data);

    RooDataSet* GetSDataSet() const;

    RooArgList GetSWeightVars() const;

    Int_t GetNumSWeightVars() const;

    void AddSWeight(RooAbsPdf* pdf, const RooArgList &yieldsTmp,
          const RooArgSet &projDeps=RooArgSet(), bool includeWeights=kTRUE,
          const RooCmdArg& fitToarg5=RooCmdArg::none(),
          const RooCmdArg& fitToarg6=RooCmdArg::none(),
          const RooCmdArg& fitToarg7=RooCmdArg::none(),
          const RooCmdArg& fitToarg8=RooCmdArg::none());

    void AddSWeight(RooAbsPdf* pdf, const RooArgList &allYieldsList, const RooArgList &fixedYields,
          const RooArgSet &projDeps=RooArgSet(), bool includeWeights=kTRUE,
          const RooCmdArg& fitToarg5=RooCmdArg::none(),
          const RooCmdArg& fitToarg6=RooCmdArg::none(),
          const RooCmdArg& fitToarg7=RooCmdArg::none(),
          const RooCmdArg& fitToarg8=RooCmdArg::none());

    Double_t GetSumOfEventSWeight(Int_t numEvent) const;

    Double_t GetYieldFromSWeight(const char* sVariable) const;

    Double_t GetSWeight(Int_t numEvent, const char* sVariable) const;



  protected:

     enum {
        kOwnData = BIT(20)
     };

    RooArgList fSWeightVars;
    //  RooListProxy fSWeightVars;
    RooArgList fSWeightCoefs;


    RooDataSet* fSData;

    //ClassDef(MySPlot,1)   // Class used for making MySPlots


      };

//}
#endif