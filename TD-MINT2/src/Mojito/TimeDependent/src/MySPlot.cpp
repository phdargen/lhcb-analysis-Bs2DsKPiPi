#include <vector>
#include <map>

#include "Mint/MySPlot.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGlobalFunc.h"
#include "TTree.h"
#include "RooStats/RooStatsUtils.h"
#include "TCanvas.h"
#include "RooPlot.h"


#include "TMatrixD.h"


//ClassImp(RooStats::MySPlot); ;

using namespace RooStats;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

MySPlot::~MySPlot()
{
   if(TestBit(kOwnData) && fSData)
      delete fSData;

}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

MySPlot::MySPlot():
  TNamed()
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = NULL;

}

////////////////////////////////////////////////////////////////////////////////

MySPlot::MySPlot(const char* name, const char* title):
  TNamed(name, title)
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = NULL;

}

////////////////////////////////////////////////////////////////////////////////
///Constructor from a RooDataSet
///No sWeighted variables are present

MySPlot::MySPlot(const char* name, const char* title, const RooDataSet &data):
  TNamed(name, title)
{
  RooArgList Args;

  fSWeightVars = Args;

  fSData = (RooDataSet*) &data;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy Constructor from another MySPlot

MySPlot::MySPlot(const MySPlot &other):
  TNamed(other)
{
  RooArgList Args = (RooArgList) other.GetSWeightVars();

  fSWeightVars.addClone(Args);

  fSData = (RooDataSet*) other.GetSDataSet();

}

////////////////////////////////////////////////////////////////////////////////
///Construct a new MySPlot instance, calculate sWeights, and include them
///in the RooDataSet held by this instance.
///
/// The constructor automatically calls AddSWeight() to add s weights to the dataset.
/// These can be retrieved later using GetSWeight() or GetSDataSet().
///\param[in] name Name of the instance.
///\param[in] title Title of the instance.
///\param[in] data Dataset to fit to.
///\param[in] pdf PDF to compute s weights for.
///\param[in] yieldsList List of parameters in `pdf` that are yields.
///\param[in] projDeps Don't normalise over these parameters when calculating the sWeights. Will be passed on to AddSWeight().
///\param[in] useWeights Include weights of the input data in calculation of s weights.
///\param[in] cloneData Make a clone of the incoming data before adding weights.
///\param[in] newName New name for the data.
///\param[in] argX Additional arguments for the fitting step in AddSWeight().
MySPlot::MySPlot(const char* name, const char* title, RooDataSet& data, RooAbsPdf* pdf,
        const RooArgList &yieldsList, const RooArgSet &projDeps,
        bool useWeights, bool cloneData, const char* newName,
        const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8):
  TNamed(name, title)
{
   if(cloneData == 1) {
    fSData = (RooDataSet*) data.Clone(newName);
    SetBit(kOwnData);
   }
  else
    fSData = (RooDataSet*) &data;

  // Add check that yieldsList contains all RooRealVars
  TIterator* iter = yieldsList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dynamic_cast<RooRealVar*>(arg)) {
      coutE(InputArguments) << "SPlot::SPlot(" << GetName() << ") input argument " 
			    << arg->GetName() << " is not of type RooRealVar " << endl ;
      throw string(Form("SPlot::SPlot(%s) input argument %s is not of type RooRealVar",GetName(),arg->GetName())) ;
    }
  }
  delete iter ;

  //Construct a new MySPlot class,
  //calculate sWeights, and include them
  //in the RooDataSet of this class.

  this->AddSWeight(pdf, yieldsList, projDeps, useWeights, arg5, arg6, arg7, arg8);
}


MySPlot::MySPlot(const char* name, const char* title, RooDataSet& data, RooAbsPdf* pdf,
        const RooArgList &allYieldsList, const RooArgList &fixedYields, const RooArgSet &projDeps,
        bool useWeights, bool cloneData, const char* newName,
        const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8):
  TNamed(name, title)
{
   if(cloneData == 1) {
    fSData = (RooDataSet*) data.Clone(newName);
    SetBit(kOwnData);
   }
  else
    fSData = (RooDataSet*) &data;

  // Add check that yieldsList contains all RooRealVars
  TIterator* iter = allYieldsList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dynamic_cast<RooRealVar*>(arg)) {
      coutE(InputArguments) << "SPlot::SPlot(" << GetName() << ") input argument " 
			    << arg->GetName() << " is not of type RooRealVar " << endl ;
      throw string(Form("SPlot::SPlot(%s) input argument %s is not of type RooRealVar",GetName(),arg->GetName())) ;
    }
  }
  delete iter ;

  // add check that fixed yields are in allYieldsList and are constant
  iter = fixedYields.createIterator(); 
  while((arg=(RooAbsArg*)iter->Next())){
    RooRealVar* varTemp = ( dynamic_cast<RooRealVar*>(arg));
    if ( varTemp && varTemp->isConstant() == 0 ){
      varTemp->setConstant();
      coutE(InputArguments) << "SPlot::SPlot(" << GetName() << ") fixed yield input argument "
          << arg->GetName() << " was not constant, has been set constant. " << endl ;
    }
    if (!(allYieldsList.contains(*arg)))
      coutE(InputArguments)<<string(Form("SPlot::SPlot(%s) fixed yield %s is not in allYieldsList", GetName(), arg->GetName()))<<endl;
    }

  delete iter;

  //Construct a new MySPlot class,
  //calculate sWeights, and include them
  //in the RooDataSet of this class.

  this->AddSWeight(pdf, allYieldsList, fixedYields, projDeps, useWeights, arg5, arg6, arg7, arg8);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dataset (if not passed in constructor).
RooDataSet* MySPlot::SetSData(RooDataSet* data)
{
  if(data)    {
    fSData = (RooDataSet*) data;
    return fSData;
  }  else
    return NULL;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve s-weighted data.
/// It does **not** automatically call AddSWeight(). This needs to be done manually.
RooDataSet* MySPlot::GetSDataSet() const
{
  return fSData;
}

////////////////////////////////////////////////////////////////////////////////
/// Retrieve an s weight.
/// \param[in] numEvent Event number to retrieve s weight for.
/// \param[in] sVariable The yield parameter to retrieve the s weight for.
Double_t MySPlot::GetSWeight(Int_t numEvent, const char* sVariable) const
{
  if(numEvent > fSData->numEntries() )
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  if(numEvent < 0)
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  Double_t totalYield = 0;

  std::string varname(sVariable);
  varname += "_sw";


  if(fSWeightVars.find(sVariable) )
    {
      RooArgSet Row(*fSData->get(numEvent));
      totalYield += Row.getRealValue(sVariable);

      return totalYield;
    }

  if( fSWeightVars.find(varname.c_str())  )
    {

      RooArgSet Row(*fSData->get(numEvent));
      totalYield += Row.getRealValue(varname.c_str() );

      return totalYield;
    }

  else
    coutE(InputArguments) << "InputVariable not in list of sWeighted variables" << endl;

  return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Sum the SWeights for a particular event.
/// This sum should equal the total weight of that event.
/// This method is intended to be used as a check.

Double_t MySPlot::GetSumOfEventSWeight(Int_t numEvent) const
{
  if(numEvent > fSData->numEntries() )
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  if(numEvent < 0)
    {
      coutE(InputArguments)  << "Invalid Entry Number" << endl;
      return -1;
    }

  Int_t numSWeightVars = this->GetNumSWeightVars();

  Double_t eventSWeight = 0;

  RooArgSet Row(*fSData->get(numEvent));

  for (Int_t i = 0; i < numSWeightVars; i++)
    eventSWeight += Row.getRealValue(fSWeightVars.at(i)->GetName() );

  return  eventSWeight;
}

////////////////////////////////////////////////////////////////////////////////
/// Sum the SWeights for a particular species over all events.
/// This should equal the total (weighted) yield of that species.
/// This method is intended as a check.

Double_t MySPlot::GetYieldFromSWeight(const char* sVariable) const
{

  Double_t totalYield = 0;

  std::string varname(sVariable);
  varname += "_sw";


  if(fSWeightVars.find(sVariable) )
    {
      for(Int_t i=0; i < fSData->numEntries(); i++)
   {
     RooArgSet Row(*fSData->get(i));
     totalYield += Row.getRealValue(sVariable);
   }

      return totalYield;
    }

  if( fSWeightVars.find(varname.c_str())  )
    {
      for(Int_t i=0; i < fSData->numEntries(); i++)
   {
     RooArgSet Row(*fSData->get(i));
     totalYield += Row.getRealValue(varname.c_str() );
   }

      return totalYield;
    }

  else
    coutE(InputArguments) << "InputVariable not in list of sWeighted variables" << endl;

  return -1;
}


////////////////////////////////////////////////////////////////////////////////
/// Return a RooArgList containing all paramters that have s weights.

RooArgList MySPlot::GetSWeightVars() const
{

  RooArgList Args = fSWeightVars;

  return  Args;

}

////////////////////////////////////////////////////////////////////////////////
/// Return the number of SWeights
/// In other words, return the number of
/// species that we are trying to extract.

Int_t MySPlot::GetNumSWeightVars() const
{
  RooArgList Args = fSWeightVars;

  return Args.getSize();
}

////////////////////////////////////////////////////////////////////////////////
/// Method which adds the sWeights to the dataset.
///
/// The MySPlot will contain two new variables for each yield parameter:
/// - `L_<varname>` is the the likelihood for each event, *i.e.*, the pdf evaluated for the a given value of the variable "varname".
/// - `<varname>_sw` is the value of the sWeight for the variable "varname" for each event.
///
/// Find Parameters in the PDF to be considered fixed when calculating the SWeights
/// and be sure to NOT include the yields in that list.
///
/// After fixing non-yield parameters, this function will start a fit by calling
/// ```
/// pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1)).
/// ```
/// One can pass additional arguments to `fitTo`, such as `RooFit::Range("fitrange")`, as `arg5`, `arg6`, `arg7`, `arg8`.
///
/// \note A `RooFit::Range` may be necessary to get expected results if you initially fit in a range
/// and/or called `pdf->fixCoefRange("fitrange")` on `pdf`.
/// Pass `arg5`, `arg6`, `arg7`, `arg8` AT YOUR OWN RISK.
///
/// \param[in] pdf PDF to fit to data to compute s weights.
/// \param[in] yieldsTmp Yields to use to compute s weights.
/// \param[in] projDeps These will not be normalized over when calculating the sWeights,
/// and will be considered parameters, not observables.
/// \param[in] includeWeights Include weights of the input data in calculation of s weights.
/// \param[in] argX Optional additional arguments for the fitting step.
void MySPlot::AddSWeight( RooAbsPdf* pdf, const RooArgList &yieldsTmp,
         const RooArgSet &projDeps, bool includeWeights,
         const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8)
{

  RooFit::MsgLevel currentLevel =  RooMsgService::instance().globalKillBelow();

  // Find Parameters in the PDF to be considered fixed when calculating the SWeights
  // and be sure to NOT include the yields in that list
  RooArgList* constParameters = (RooArgList*)pdf->getParameters(fSData) ;
  constParameters->remove(yieldsTmp, kTRUE, kTRUE);


  // Set these parameters constant and store them so they can later
  // be set to not constant
  std::vector<RooRealVar*> constVarHolder;

  for(Int_t i = 0; i < constParameters->getSize(); i++)
    {
      RooRealVar* varTemp = ( dynamic_cast<RooRealVar*>( constParameters->at(i) ) );
      if(varTemp &&  varTemp->isConstant() == 0 )
   {
     varTemp->setConstant();
     constVarHolder.push_back(varTemp);
   }
    }

  // Fit yields to the data with all other variables held constant
  // This is necessary because MySPlot assumes the yields minimise -Log(likelihood)

  pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), arg5, arg6, arg7, arg8);
  TCanvas* c = new TCanvas();
  RooArgList* model_obs = (RooArgList*)pdf->getObservables(*fSData) ;
  model_obs->Print("v") ;
  RooPlot* frame= (dynamic_cast<RooRealVar*>(model_obs->at(0)))->frame();
  frame->SetTitle("");
  fSData->plotOn(frame,Binning(100));
  pdf->plotOn(frame,LineColor(kBlue));
  frame->Draw();
  c->Print("sWeightFit.eps");

  // Hold the value of the fitted yields
  std::vector<double> yieldsHolder;

  for(Int_t i = 0; i < yieldsTmp.getSize(); i++)
    yieldsHolder.push_back( ((RooRealVar*) yieldsTmp.at(i))->getVal());

  Int_t nspec = yieldsTmp.getSize();
  RooArgList yields = *(RooArgList*)yieldsTmp.snapshot(kFALSE);

  if(currentLevel <= RooFit::DEBUG)
    {
      coutI(InputArguments) << "Printing Yields" << endl;
      yields.Print();
    }

  // The list of variables to normalize over when calculating PDF values.

  RooArgSet vars(*fSData->get() );
  vars.remove(projDeps, kTRUE, kTRUE);

  // Attach data set

  // const_cast<RooAbsPdf*>(pdf)->attachDataSet(*fSData);

  pdf->attachDataSet(*fSData);

  // first calculate the pdf values for all species and all events
  std::vector<RooRealVar*> yieldvars ;
  RooArgSet* parameters = pdf->getParameters(fSData) ;

  std::vector<Double_t> yieldvalues ;
  for (Int_t k = 0; k < nspec; ++k)
    {
      RooRealVar* thisyield = dynamic_cast<RooRealVar*>(yields.at(k)) ;
      if (thisyield) {
         RooRealVar* yieldinpdf = dynamic_cast<RooRealVar*>(parameters->find(thisyield->GetName() )) ;

         if (yieldinpdf) {
            coutI(InputArguments)<< "yield in pdf: " << yieldinpdf->GetName() << " " << thisyield->getVal() << endl;

            yieldvars.push_back(yieldinpdf) ;
            yieldvalues.push_back(thisyield->getVal()) ;
         }
      }
    }

  Int_t numevents = fSData->numEntries() ;

  std::vector<std::vector<Double_t> > pdfvalues(numevents,std::vector<Double_t>(nspec,0)) ;


  // set all yield to zero
  for(Int_t m=0; m<nspec; ++m) yieldvars[m]->setVal(0) ;


  // For every event and for every specie,
  // calculate the value of the component pdf for that specie
  // by setting the yield of that specie to 1
  // and all others to 0.  Evaluate the pdf for each event
  // and store the values.

  RooArgSet * pdfvars = pdf->getVariables();

  for (Int_t ievt = 0; ievt <numevents; ievt++)
    {
      //   if (ievt % 100 == 0)
      //  coutP(Eval)  << ".";


      //FIX THIS PART, EVALUATION PROGRESS!!

      RooStats::SetParameters(fSData->get(ievt), pdfvars);

      //   RooArgSet row(*fSData->get(ievt));

      for(Int_t k = 0; k < nspec; ++k)
   {
     //Check that range of yields is at least (0,1), and fix otherwise
     if(yieldvars[k]->getMin() > 0)
       {
         coutW(InputArguments)  << "Minimum Range for " << yieldvars[k]->GetName() << " must be 0.  ";
         coutW(InputArguments)  << "Setting min range to 0" << std::endl;
         yieldvars[k]->setMin(0);
       }

     if(yieldvars[k]->getMax() < 1)
       {
         coutW(InputArguments)  << "Maximum Range for " << yieldvars[k]->GetName() << " must be 1.  ";
         coutW(InputArguments)  << "Setting max range to 1" << std::endl;
         yieldvars[k]->setMax(1);
       }

     // set this yield to 1
     yieldvars[k]->setVal( 1 ) ;
     // evaluate the pdf
     Double_t f_k = pdf->getVal(&vars) ;
     pdfvalues[ievt][k] = f_k ;
     if( !(f_k>1 || f_k<1) )
       coutW(InputArguments) << "Strange pdf value: " << ievt << " " << k << " " << f_k << std::endl ;
     yieldvars[k]->setVal( 0 ) ;
   }
    }
  delete pdfvars;

  // check that the likelihood normalization is fine
  std::vector<Double_t> norm(nspec,0) ;
  for (Int_t ievt = 0; ievt <numevents ; ievt++)
    {
      Double_t dnorm(0) ;
      for(Int_t k=0; k<nspec; ++k) dnorm += yieldvalues[k] * pdfvalues[ievt][k] ;
      for(Int_t j=0; j<nspec; ++j) norm[j] += pdfvalues[ievt][j]/dnorm ;
    }

  coutI(Contents) << "likelihood norms: "  ;

  for(Int_t k=0; k<nspec; ++k)  coutI(Contents) << norm[k] << " " ;
  coutI(Contents) << std::endl ;

  // Make a TMatrixD to hold the covariance matrix.
  TMatrixD covInv(nspec, nspec);
  for (Int_t i = 0; i < nspec; i++) for (Int_t j = 0; j < nspec; j++) covInv(i,j) = 0;

  coutI(Contents) << "Calculating covariance matrix";


  // Calculate the inverse covariance matrix, using weights
  for (Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // Calculate contribution to the inverse of the covariance
      // matrix. See BAD 509 V2 eqn. 15

      // Sum for the denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nspec; ++k)
   dsum += pdfvalues[ievt][k] * yieldvalues[k] ;

      for(Int_t n=0; n<nspec; ++n)
   for(Int_t j=0; j<nspec; ++j)
     {
       if(includeWeights)
         covInv(n,j) +=  fSData->weight()*pdfvalues[ievt][n]*pdfvalues[ievt][j]/(dsum*dsum) ;
       else
         covInv(n,j) +=                   pdfvalues[ievt][n]*pdfvalues[ievt][j]/(dsum*dsum) ;
     }

      //ADDED WEIGHT ABOVE

    }

  // Covariance inverse should now be computed!

  // Invert to get the covariance matrix
  if (covInv.Determinant() <=0)
    {
      coutE(Eval) << "MySPlot Error: covariance matrix is singular; I can't invert it!" << std::endl;
      covInv.Print();
      return;
    }

  TMatrixD covMatrix(TMatrixD::kInverted,covInv);

  //check cov normalization
  if(currentLevel <= RooFit::DEBUG)
    {
      coutI(Eval) << "Checking Likelihood normalization:  " << std::endl;
      coutI(Eval) << "Yield of specie  Sum of Row in Matrix   Norm" << std::endl;
      for(Int_t k=0; k<nspec; ++k)
   {
     Double_t covnorm(0) ;
     for(Int_t m=0; m<nspec; ++m) covnorm += covInv[k][m]*yieldvalues[m] ;
     Double_t sumrow(0) ;
     for(Int_t m = 0; m < nspec; ++m) sumrow += covMatrix[k][m] ;
     coutI(Eval)  << yieldvalues[k] << " " << sumrow << " " << covnorm << endl ;
   }
    }

  // calculate for each event the sWeight (BAD 509 V2 eq. 21)
  coutI(Eval) << "Calculating sWeight" << std::endl;
  std::vector<RooRealVar*> sweightvec ;
  std::vector<RooRealVar*> pdfvec ;
  RooArgSet sweightset ;

  // Create and label the variables
  // used to store the SWeights

  fSWeightVars.Clear();

  for(Int_t k=0; k<nspec; ++k)
    {
       std::string wname = std::string(yieldvars[k]->GetName()) + "_sw";
       RooRealVar* var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
       sweightvec.push_back( var) ;
       sweightset.add(*var) ;
       fSWeightVars.add(*var);

       wname = "L_" + std::string(yieldvars[k]->GetName());
       var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
       pdfvec.push_back( var) ;
       sweightset.add(*var) ;
    }

  // Create and fill a RooDataSet
  // with the SWeights

  RooDataSet* sWeightData = new RooDataSet("dataset", "dataset with sWeights", sweightset);

  for(Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // sum for denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nspec; ++k)   dsum +=  pdfvalues[ievt][k] * yieldvalues[k] ;
      // covariance weighted pdf for each specief
      for(Int_t n=0; n<nspec; ++n)
   {
     Double_t nsum(0) ;
     for(Int_t j=0; j<nspec; ++j) nsum += covMatrix(n,j) * pdfvalues[ievt][j] ;


     //Add the sWeights here!!
     //Include weights,
     //ie events weights are absorbed into sWeight


     if(includeWeights) sweightvec[n]->setVal(fSData->weight() * nsum/dsum) ;
     else  sweightvec[n]->setVal( nsum/dsum) ;

     pdfvec[n]->setVal( pdfvalues[ievt][n] ) ;

     if( !(fabs(nsum/dsum)>=0 ) )
       {
         coutE(Contents) << "error: " << nsum/dsum << endl ;
         return;
       }
   }

      sWeightData->add(sweightset) ;
    }


  // Add the SWeights to the original data set

  fSData->merge(sWeightData);

  delete sWeightData;

  //Restore yield values

  for(Int_t i = 0; i < yieldsTmp.getSize(); i++)
    ((RooRealVar*) yieldsTmp.at(i))->setVal(yieldsHolder.at(i));

  //Make any variables that were forced to constant no longer constant

  for(Int_t i=0; i < (Int_t) constVarHolder.size(); i++)
    constVarHolder.at(i)->setConstant(kFALSE);

  return;

}

void MySPlot::AddSWeight( RooAbsPdf* pdf, const RooArgList &allYieldsList,
			const RooArgList &fixedYields, const RooArgSet &projDeps, bool includeWeights,const RooCmdArg& arg5, const RooCmdArg& arg6, const RooCmdArg& arg7, const RooCmdArg& arg8)
{
  RooFit::MsgLevel currentLevel =  RooMsgService::instance().globalKillBelow();

  // Find Parameters in the PDF to be considered fixed when calculating the SWeights
  // and be sure to NOT include the yields in that list
  RooArgList* constParameters = (RooArgList*)pdf->getParameters(fSData) ;
  constParameters->remove(allYieldsList, kTRUE, kTRUE);


  // Set the parameters that are not constant constant 
  // and store them so they can later be set to not constant
  std::vector<RooRealVar*> constVarHolder;

  for(Int_t i = 0; i < constParameters->getSize(); i++)
    {
      RooRealVar* varTemp = ( dynamic_cast<RooRealVar*>( constParameters->at(i) ) );
      if(varTemp &&  varTemp->isConstant() == 0 )
   {
     varTemp->setConstant();
     constVarHolder.push_back(varTemp);
   }
    }

  TIterator *it = allYieldsList.createIterator();
  RooAbsArg* arg ;
  unsigned int iArg(0);
  std::vector<unsigned int> varIndexes ;
  std::vector<unsigned int> fixedIndexes ;
  while ((arg = (RooAbsArg*) it->Next()) != NULL){
    if (!(fixedYields.find(arg->GetName()))) 
      varIndexes.push_back(iArg);
    else fixedIndexes.push_back(iArg);
    iArg++;}
  // We now have two indexes over which we iterate
  Int_t nAllSpec = allYieldsList.getSize();
  Int_t nVarSpec = allYieldsList.getSize() - fixedYields.getSize();

  // Fit yields to the data with all other variables held constant
  // This is necessary because SPlot assumes the yields minimise -Log(likelihood)

  //pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::PrintLevel(-1), RooFit::PrintEvalErrors(-1), arg5, arg6, arg7, arg8);
  pdf->fitTo(*fSData, RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), arg5, arg6, arg7, arg8);


   TCanvas* c = new TCanvas();
    RooArgList* model_obs = (RooArgList*)pdf->getObservables(*fSData) ;
    model_obs->Print("v") ;

   	RooPlot* frame= (dynamic_cast<RooRealVar*>(model_obs->at(0)))->frame();
	frame->SetTitle("");
 	fSData->plotOn(frame,Binning(100));
	pdf->plotOn(frame,LineColor(kBlue));
	frame->Draw();
	c->Print("sWeightFit.eps");

  // Hold the value of the fitted yields
  std::vector<double> yieldsHolder;

  for(Int_t i = 0; i < nAllSpec; i++)
    yieldsHolder.push_back( ((RooRealVar*) allYieldsList.at(i))->getVal());

  // Int_t nspec = yieldsTmp.getSize();
  // no DeepCopy snapshot: kFalse
  RooArgList yields = *(RooArgList*)allYieldsList.snapshot(kFALSE);

  if(currentLevel <= RooFit::DEBUG)
    {
      coutI(InputArguments) << "Printing Yields" << endl;
      yields.Print();
    }

  // The list of variables to normalize over when calculating PDF values.

  RooArgSet vars(*fSData->get() );
  vars.remove(projDeps, kTRUE, kTRUE);

  // Attach data set

  // const_cast<RooAbsPdf*>(pdf)->attachDataSet(*fSData);

  pdf->attachDataSet(*fSData);

  // first calculate the pdf values for all species and all events
  std::vector<RooRealVar*> yieldvars ;
  RooArgSet* parameters = pdf->getParameters(fSData) ;

  std::vector<Double_t> yieldvalues ;
  for (Int_t k = 0; k < nAllSpec; ++k)
    {
      RooRealVar* thisyield = dynamic_cast<RooRealVar*>(yields.at(k)) ;
      if (thisyield) {
         RooRealVar* yieldinpdf = dynamic_cast<RooRealVar*>(parameters->find(thisyield->GetName() )) ;

         if (yieldinpdf) {
            coutI(InputArguments)<< "yield in pdf: " << yieldinpdf->GetName() << " " << thisyield->getVal() << endl;

            yieldvars.push_back(yieldinpdf) ;
            yieldvalues.push_back(thisyield->getVal()) ;
         }
      }
    }

  Int_t numevents = fSData->numEntries() ;

  std::vector<std::vector<Double_t> > pdfvalues(numevents,std::vector<Double_t>(nAllSpec,0)) ;


  // set all yield to zero
  for(Int_t m=0; m<nAllSpec; ++m) yieldvars[m]->setVal(0) ;


  // For every event and for every specie,
  // calculate the value of the component pdf for that specie
  // by setting the yield of that specie to 1
  // and all others to 0.  Evaluate the pdf for each event
  // and store the values.

  RooArgSet * pdfvars = pdf->getVariables();

  for (Int_t ievt = 0; ievt <numevents; ievt++)
    {
      //   if (ievt % 100 == 0)
      //  coutP(Eval)  << ".";


      //FIX THIS PART, EVALUATION PROGRESS!!

      RooStats::SetParameters(fSData->get(ievt), pdfvars);

      //   RooArgSet row(*fSData->get(ievt));

      for(Int_t k = 0; k < nAllSpec; ++k)
   {
     //Check that range of yields is at least (0,1), and fix otherwise
     if(yieldvars[k]->getMin() > 0)
       {
         coutW(InputArguments)  << "Minimum Range for " << yieldvars[k]->GetName() << " must be 0.  ";
         coutW(InputArguments)  << "Setting min range to 0" << std::endl;
         yieldvars[k]->setMin(0);
       }

     if(yieldvars[k]->getMax() < 1)
       {
         coutW(InputArguments)  << "Maximum Range for " << yieldvars[k]->GetName() << " must be 1.  ";
         coutW(InputArguments)  << "Setting max range to 1" << std::endl;
         yieldvars[k]->setMax(1);
       }

     // set this yield to 1
     yieldvars[k]->setVal( 1 ) ;
     // evaluate the pdf
     Double_t f_k = pdf->getVal(&vars) ;
     pdfvalues[ievt][k] = f_k ;
     if( !(f_k>1 || f_k<1) )
       coutW(InputArguments) << "Strange pdf value: " << ievt << " " << k << " " << f_k << " setting to 0 " << std::endl ;
     yieldvars[k]->setVal( 0 ) ;
   }
    }
  delete pdfvars;

  // check that the likelihood normalization is fine
  std::vector<Double_t> norm(nAllSpec,0) ;
  for (Int_t ievt = 0; ievt < numevents ; ievt++)
    {
      Double_t dnorm(0) ;
      for(Int_t k=0; k<nAllSpec; ++k) dnorm += yieldvalues[k] * pdfvalues[ievt][k] ;
      for(Int_t j=0; j<nAllSpec; ++j) norm[j] += pdfvalues[ievt][j]/dnorm ;
    }

  coutI(Contents) << "likelihood norms: "  ;

  for(Int_t k=0; k<nAllSpec; ++k)  coutI(Contents) << norm[k] << " " ;
  coutI(Contents) << std::endl ;

  // Make a TMatrixD to hold the inverse covariance matrix.
  TMatrixD covInv(nVarSpec, nVarSpec);
  for (Int_t i = 0; i < nVarSpec; i++) for (Int_t j = 0; j < nVarSpec; j++) covInv(i,j) = 0;

  coutI(Contents) << "Calculating covariance matrix";


  // Calculate the inverse covariance matrix, using weights
  for (Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // Sum for the denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nAllSpec; ++k)
	dsum += pdfvalues[ievt][k] * yieldvalues[k] ;

      for(Int_t n=0; n<nVarSpec; ++n)
	for(Int_t j=0; j<nVarSpec; ++j){
	  if(includeWeights == kTRUE)
	    covInv(n,j) +=  fSData->weight()*pdfvalues[ievt][varIndexes[n]]*pdfvalues[ievt][varIndexes[j]]/(dsum*dsum) ;
	  else
	    covInv(n,j) +=  pdfvalues[ievt][varIndexes[n]]*pdfvalues[ievt][varIndexes[j]]/(dsum*dsum) ;
	}
    }

  // Invert to get the covariance matrix
  if (covInv.Determinant() <=0)
    {
      coutE(Eval) << "SPlot Error: covariance matrix is singular; I can't invert it!" << std::endl;
      covInv.Print();
      return;
    }

  TMatrixD covMatrix(TMatrixD::kInverted,covInv);
coutI(Eval) << "Covariance matrix calculated" << endl;
covMatrix.Print();

  //check cov normalization
  if(currentLevel <= RooFit::DEBUG)
    {
      coutI(Eval) << "Checking Likelihood normalization:  " << std::endl;
      coutI(Eval) << "Yield of species  Sum of Row in Matrix   Norm" << std::endl;
      for(Int_t k=0; k<nVarSpec; ++k)
   {
     Double_t covnorm(0) ;
     for(Int_t m=0; m<nVarSpec; ++m) covnorm += covInv[k][m]*yieldvalues[varIndexes[m]] ;
     Double_t sumrow(0) ;
     for(Int_t m = 0; m < nVarSpec; ++m) sumrow += covMatrix[k][m] ;
     coutI(Eval)  << yieldvalues[varIndexes[k]] << " " << sumrow << " " << covnorm << endl ;
   }
    }

  // calculate for each event the sWeight (BAD 509 V2 eq. 21)
  coutI(Eval) << "Calculating sWeight" << std::endl;
  std::vector<RooRealVar*> sweightvec ;
  // vector for adapted sweights with fixed yields
  std::vector<RooRealVar*> sweightvec_fy ;
  std::vector<RooRealVar*> pdfvec ;
  RooArgSet sweightset ;
  // fixed yields weights
  RooArgSet sweightset_fy ;
  Double_t coeff_num_fy ;
  coeff_num_fy = numevents;

  // Create and label the variables
  // used to store the SWeights

  fSWeightVars.Clear();
  fSWeightCoefs.Clear();

  // can only define sweights and coefficients for the variable species
  for(Int_t k=0; k<nVarSpec; ++k)
    {
       std::string wname = std::string(yieldvars[varIndexes[k]]->GetName()) + "_sw_varyields";
       RooRealVar* var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
       sweightvec.push_back( var) ;
       sweightset.add(*var) ;
       fSWeightVars.add(*var);

       wname = std::string(yieldvars[varIndexes[k]]->GetName()) + "_c";
       var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
      double cVal = yieldvalues[varIndexes[k]];
      for(Int_t n=0; n<nVarSpec; ++n)
      {
        coeff_num_fy -= covMatrix[k][n];
	      cVal -= covMatrix[k][n];
      }
      var->setVal(cVal);
      coutI(Eval) << "Coefficient " << wname << " = " << cVal << endl;
      fSWeightCoefs.add(*var); // new attribute of the class
    }
  coutI(Eval) << "Denominator for sP0 N - sum(V_ij) = " << coeff_num_fy << endl;

  // can define extended sweights and pdf valuesfor fixed yields for all species
  for(Int_t k=0; k<nAllSpec; ++k)
    {
      std::string wname = std::string(yieldvars[k]->GetName()) + "_sw";
      RooRealVar* var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
      sweightvec_fy.push_back( var) ;
      sweightset_fy.add(*var) ;
      fSWeightVars.add(*var);

      wname = "L_" + std::string(yieldvars[k]->GetName());
      var = new RooRealVar(wname.c_str(),wname.c_str(),0) ;
      pdfvec.push_back( var) ;
      sweightset_fy.add(*var) ;
    }


  // Create and fill a RooDataSet
  // with the SWeights

  RooDataSet* sWeightData = new RooDataSet("dataset", "dataset with sWeights", sweightset_fy);

  for(Int_t ievt = 0; ievt < numevents; ++ievt)
    {

      fSData->get(ievt) ;

      // sum for denominator
      Double_t dsum(0);
      for(Int_t k = 0; k < nAllSpec; ++k)   
	      {
        dsum +=  pdfvalues[ievt][k] * yieldvalues[k] ;
        pdfvec[k]->setVal(pdfvalues[ievt][k]);
        }
      // covariance weighted pdf for each species
  // sweight-like attribute of fixed yields
  Double_t sP0(1);
  // calculate nominal sweights
  for(Int_t n=0; n<nVarSpec; ++n)
	{
	  Double_t nsum(0) ;
	  for(Int_t j=0; j<nVarSpec; ++j) 
	    nsum += covMatrix(n,j) * pdfvalues[ievt][varIndexes[j]] ;

	  //Add the sWeights here!!
	  //Include weights,
	  //ie events weights are absorbed into sWeight


	  if(includeWeights == kTRUE) 
	    sweightvec[n]->setVal(fSData->weight() * nsum/dsum) ;
	  else  
	    sweightvec[n]->setVal( nsum/dsum) ;

    sP0 -= nsum/dsum;

	  if( !(fabs(nsum/dsum)>=0 ) )
	    {
	      coutE(Contents) << "error: " << nsum/dsum << endl ;
	      return;
	    }

	}
  for(Int_t n=0; n<nVarSpec; ++n)
    sweightvec_fy[varIndexes[n]]->setVal(sweightvec[n]->getVal() + sP0*(dynamic_cast<RooRealVar*>(fSWeightCoefs.at(n))->getVal())/coeff_num_fy) ;

  // add sweights for the fixed species according to the ratio of their pdfs times yields for this event
  Double_t fixeddsum(0);
  
  for(Int_t h=0; h< nAllSpec-nVarSpec; ++h)
    fixeddsum +=  pdfvalues[ievt][fixedIndexes[h]] * yieldvalues[fixedIndexes[h]] ;
  for(Int_t h=0; h< nAllSpec-nVarSpec; ++h) {
    if(fixeddsum==0) {
        Double_t fixedweight(sP0/(nAllSpec-nVarSpec));
        // fixed pdfs are 0, setting fixed weights to sP0/nFixed
        sweightvec_fy[fixedIndexes[h]]->setVal(fixedweight);
        }
    else
        sweightvec_fy[fixedIndexes[h]]->setVal(sP0*(pdfvec[fixedIndexes[h]]->getVal()*yieldsHolder[fixedIndexes[h]]/fixeddsum));
    }
    sWeightData->add(sweightset_fy) ;
    }


  // Add the SWeights to the original data set


  fSData->merge(sWeightData);

  delete sWeightData;

  //Restore yield values

  for(Int_t i = 0; i < allYieldsList.getSize(); i++)
    ((RooRealVar*) allYieldsList.at(i))->setVal(yieldsHolder.at(i));

  //Make any variables that were forced to constant no longer constant

  for(Int_t i=0; i < (Int_t) constVarHolder.size(); i++)
    constVarHolder.at(i)->setConstant(kFALSE);

  return;

}
