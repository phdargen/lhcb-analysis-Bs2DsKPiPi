/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Uniform binning
 *
 **/

/** \class UniformBinning



*/

 
#ifndef UNIFORMBINNING_HH
#define UNIFORMBINNING_HH

// HyperPlot includes
#include "Mint/MessageService.h"
#include "Mint/HyperPoint.h"
#include "Mint/HyperPointSet.h"
#include "Mint/HyperCuboid.h"
#include "Mint/HyperVolume.h"
#include "Mint/RootPlotter1D.h"
#include "Mint/RootPlotter2D.h"
#include "Mint/HyperName.h"
#include "Mint/BinningBase.h"
#include "Mint/LoadingBar.h"
#include "Mint/CachedVar.h"


// Root includes
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

// std includes
#include <algorithm>
#include <sstream>

class UniformBinning : public BinningBase {

  private:
  
  HyperCuboid      _limits;
  std::vector<int> _nLocalBins; 

  public:
  
  UniformBinning(HyperCuboid limits, int nLocalBins);
  UniformBinning(HyperCuboid limits, std::vector<int> nLocalBins);
  
  virtual ~UniformBinning();

  int getNumLocalBins(int dimension) const;

  int getGlobalBinNumber( std::vector<int> localBinNumbers ) const;
  std::vector<int> getLocalBinNumbers( int globalBinNumber ) const;

  int getLocalBinNumber(int dim, double val) const;
  std::vector<int> getLocalBinNumbers(const HyperPoint& coords) const;

  double getLowBinEdgeLocal(int dim, int localBinNum) const;
  double getHighBinEdgeLocal(int dim, int localBinNum) const;

  HyperPoint getLowCorner(int globalBinNum) const;
  HyperPoint getHighCorner(int globalBinNum) const;

  /* Virtual functions that need to implemented from BinningBase */
  /*    These will be implemented in the derived classes */

  virtual void load(TString filename, TString option = "READ");
  virtual BinningBase* clone() const;

  virtual void save(TString filename) const;
  virtual void save() const; 

  virtual void mergeBinnings( const BinningBase& other );

  virtual int getNumBins() const;
  virtual int getBinNum(const HyperPoint& coords) const;

  virtual HyperVolume getBinHyperVolume(int binNumber) const;

  virtual HyperPoint  getAverageBinWidth() const;
  virtual HyperCuboid getLimits()          const;



};



#endif

