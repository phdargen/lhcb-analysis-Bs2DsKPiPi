/**
 * <B>HyperPlot</B>,
 * Author: Sam Harnew, sam.harnew@gmail.com ,
 * Date: Dec 2015
 *
 * Class plotting 1D HyperBinningHistograms 
 *
 **/
 
#ifndef HYPERBINNINGPAINTER1D_HH
#define HYPERBINNINGPAINTER1D_HH

// HyperPlot includes
#include "MessageService.h"
#include "RootPlotter1D.h"
#include "RootPlotter2D.h"
#include "HyperBinningPainter.h"


// Root includes
#include "TH1D.h"
#include "TH2D.h"

// std includes




class HyperBinningPainter1D : public HyperBinningPainter {
  
  protected:

  public:
  
  HyperBinningPainter1D(HyperHistogram* histogram);
  
  TH1D* getHistogram(TString histname);

  virtual void draw(TString path = "", TString option = "");

  virtual ~HyperBinningPainter1D();

};



#endif

