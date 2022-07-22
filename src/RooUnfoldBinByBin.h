//=====================================================================-*-C++-*-
//! \class  RooUnfoldBinByBinT
//! \brief  Unfolding class using the bin by bin method of conversion factors.
//! \author Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDBINBYBIN_H_
#define ROOUNFOLDBINBYBIN_H_

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "TVectorD.h"

#include "RooUnfoldFitHelpers.h"

template<class Hist, class Hist2D>
class RooUnfoldBinByBinT : public RooUnfoldT<Hist,Hist2D> {

public:
  RooUnfoldBinByBinT(); // default constructor
  RooUnfoldBinByBinT (const char*    name, const char*    title); // named constructor
  RooUnfoldBinByBinT (const TString& name, const TString& title); // named constructor
  RooUnfoldBinByBinT (const RooUnfoldBinByBinT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldBinByBinT(); // destructor
  RooUnfoldBinByBinT<Hist,Hist2D>& operator= (const RooUnfoldBinByBinT<Hist,Hist2D>& rhs); // assignment operator
  RooUnfoldBinByBinT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name=0, const char* title=0);

  TVectorD* Impl();
  virtual RooUnfolding::Algorithm GetAlgorithm() const override;
  
protected:
  virtual void Unfold() const override;
  virtual void GetCov() const override;

protected:
  // cache
  class Cache {
  public:
    TVectorD _factors;
  };
  mutable Cache _specialcache;  //!

public:
  ClassDefOverride (RooUnfoldBinByBinT, 1)  // Bin-by-bin unfolding
};

//! \class RooUnfoldBinByBin 
//! \brief specialization of RooUnfoldBinByBinT for TH1/TH2 objects
typedef RooUnfoldBinByBinT<TH1,TH2> RooUnfoldBinByBin;
#ifndef NOROOFIT
//! \class RooFitUnfoldBinByBin
//! \brief specialization of RooUnfoldBinByBinT for RooAbsReal objects
typedef RooUnfoldBinByBinT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldBinByBin;
#endif


#endif /*ROOUNFOLDBINBYBIN_H_*/
