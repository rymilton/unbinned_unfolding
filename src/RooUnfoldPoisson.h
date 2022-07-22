//=====================================================================-*-C++-*-
//! \class RooUnfoldPoissonT
//==============================================================================

#ifndef ROOUNFOLDPOISSON_HH
#define ROOUNFOLDPOISSON_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"


#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TVectorD.h"
#include "TMatrixD.h"

class TH1;
class TH2;

template<class Hist, class Hist2D>
class RooUnfoldPoissonT : public RooUnfoldT<Hist,Hist2D> {

public:

  // Standard methods

  RooUnfoldPoissonT(); // default constructor
  RooUnfoldPoissonT (const char*    name, const char*    title); // named constructor
  RooUnfoldPoissonT (const TString& name, const TString& title); // named constructor
  RooUnfoldPoissonT (const RooUnfoldPoissonT<Hist,Hist2D>& rhs); // copy constructor
  RooUnfoldPoissonT& operator= (const RooUnfoldPoissonT<Hist,Hist2D>& rhs); // assignment operator

  // Special constructors

  RooUnfoldPoissonT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Double_t regparm= 1e-30, const char* name= 0, const char* title= 0);

  virtual void  SetRegParm (Double_t parm) override;
  virtual double GetRegParm() const override;
  virtual void SetMinimizerStart(const Hist* truth);
  virtual void SetMinimizerStart(TVectorD truth);
  virtual void SetPrintLevel(Int_t print = 0);
  virtual void Reset() override;
  // virtual void Print (Option_t* option= "") const override;
  virtual RooUnfolding::Algorithm GetAlgorithm() const override;  
  
protected:
  void Assign (const RooUnfoldPoissonT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void Unfold() const override ;
  virtual void GetCov() const override ;
  virtual void GetSettings() const override;

  void setup() const;

private:
  double* Rmu(const double* truth) const;
  Double_t NegativeLLH(const double* truth) const;
  Double_t TikhonovReg(const double* truth) const;
  Double_t RegLLH(const double* truth) const;

  void MinimizeRegLLH() const;

  void Init();
  void CopyData (const RooUnfoldPoissonT<Hist,Hist2D>& rhs);

protected:
  // instance variables
  mutable Int_t _min_status;
  mutable Int_t _min_print;
  mutable Double_t _RegLLH_factor;
  mutable double _regparm;
  mutable TMatrixD _response;
  mutable TVectorD _data;
  mutable TVectorD _truth_start;
  mutable TVectorD _unfolded;
  mutable TVectorD _truth_edges;
  mutable TMatrixD _covariance;

public:
  ClassDefOverride (RooUnfoldPoissonT, 1) 
};


//! \class RooUnfoldPoisson 
//! \brief specialization of RooUnfoldPoissonT for TH1/TH2 objects
typedef RooUnfoldPoissonT<TH1,TH2> RooUnfoldPoisson;
#ifndef NOROOFIT
//! \class RooFitUnfoldPoisson
//! \brief specialization of RooUnfoldPoissonT for RooAbsReal objects
typedef RooUnfoldPoissonT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldPoisson;
#endif

#endif
