//=====================================================================-*-C++-*-
//! \class RooUnfoldInvertT
//! \brief Unfolding class using inversion of the response matrix. This does not produce
//!      good results and is designed to illustrate the need for more sophisticated
//!      unfolding techniques
//! \author Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDINVERT_H_
#define ROOUNFOLDINVERT_H_

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

class TDecompSVD;

template<class Hist, class Hist2D>
class RooUnfoldInvertT : public RooUnfoldT<Hist,Hist2D> {

public:
  RooUnfoldInvertT(); // default constructor
  RooUnfoldInvertT (const char*    name, const char*    title); // named constructor
  RooUnfoldInvertT (const TString& name, const TString& title); // named constructor
  RooUnfoldInvertT (const RooUnfoldInvertT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldInvertT(); // destructor
  RooUnfoldInvertT& operator= (const RooUnfoldInvertT<Hist,Hist2D>& rhs); // assignment operator
  RooUnfoldInvertT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name=0, const char* title=0);

  virtual RooUnfolding::Algorithm GetAlgorithm() const override;  
  virtual void Reset() override;
  TDecompSVD* Impl();
  
protected:
  virtual void Unfold() const override;
  virtual void GetCov() const override;

private:
  void Init();
  Bool_t InvertResponse() const;

protected:
  // instance variables
  mutable TDecompSVD* _svd;
  mutable TMatrixD*   _resinv;

public:
  ClassDefOverride (RooUnfoldInvertT, 1)  // Unregularised unfolding
};

//! \class RooUnfoldInvert 
//! \brief specialization of RooUnfoldInvertT for TH1/TH2 objects
typedef RooUnfoldInvertT<TH1,TH2> RooUnfoldInvert;
#ifndef NOROOFIT
//! \class RooFitUnfoldInvert
//! \brief specialization of RooUnfoldInvertT for RooAbsReal objects
typedef RooUnfoldInvertT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldInvert;
#endif

#endif /*ROOUNFOLDINVERT_H_*/
