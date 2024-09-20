//=====================================================================-*-C++-*-
//! \class RooUnfoldTUnfoldT
//! \brief Unfolding class using TUnfold from ROOT to do the actual unfolding.
//! \author Richard Claridge <richard.claridge@stfc.ac.uk> & Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDTUNFOLD_HH
#define ROOUNFOLDTUNFOLD_HH

#include "RooUnfold.h"
#include "RooUnfold/TUnfold/TUnfold.h"
#include "RooUnfoldResponse.h"

class TH1;
class TH1;
class TH2D;
class TGraph;
class TSpline;

template<class Hist, class Hist2D>
class RooUnfoldTUnfoldT : public RooUnfoldT<Hist,Hist2D> {

public:

  RooUnfoldTUnfoldT(); // default constructor
  RooUnfoldTUnfoldT (const char*    name, const char*    title); // named constructor
  RooUnfoldTUnfoldT (const TString& name, const TString& title); // named constructor
  RooUnfoldTUnfoldT (const RooUnfoldTUnfoldT& rhs); // copy constructor
  virtual ~RooUnfoldTUnfoldT(); // destructor
  RooUnfoldTUnfoldT& operator= (const RooUnfoldTUnfoldT& rhs); // assignment operator
  RooUnfoldTUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,TUnfold::ERegMode reg=TUnfold::kRegModeDerivative, Bool_t handleFakes= false, 
		     const char* name= 0, const char* title= 0);
  RooUnfoldTUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,Double_t tau,TUnfold::ERegMode reg=TUnfold::kRegModeDerivative, Bool_t handleFakes= false,
		     const char* name= 0, const char* title= 0);

  virtual void Reset() override;
  TUnfold* Impl();
  void FixTau(Double_t tau);
  void OptimiseTau();
  virtual void SetRegParm(Double_t parm) override;
  Double_t GetTau() const;
  const TGraph*  GetLCurve()  const;
  const TSpline* GetLogTauX() const;
  const TSpline* GetLogTauY() const;
  virtual Double_t GetRegParm() const override;
  void SetRegMethod (TUnfold::ERegMode regmethod);
  TUnfold::ERegMode GetRegMethod() const;

protected:
  void Init();
  void Destroy();
  virtual void Unfold() const override;
  virtual void GetCov() const override;
  virtual void GetSettings() const override;
  void Assign   (const RooUnfoldTUnfoldT& rhs); // implementation of assignment operator
  void CopyData (const RooUnfoldTUnfoldT& rhs);

private:
  TUnfold::ERegMode _reg_method; //Regularisation method
  mutable TUnfold* _unf; //! Implementation in TUnfold object (no streamer)
  mutable Bool_t tau_set;
  mutable Double_t _tau;
  mutable Bool_t _handleFakes;
  mutable   TSpline* _logTauX;
  mutable   TSpline* _logTauY;
  mutable   TGraph*  _lCurve;
  mutable   TGraph*  _logTauSURE;
  mutable   TGraph*  _df_chi2A;

public:

  ClassDefOverride (RooUnfoldTUnfoldT, 1)   // Interface to TUnfold
};

//! \class RooUnfoldTUnfold 
//! \brief specialization of RooUnfoldTUnfoldT for TH1/TH2 objects
typedef RooUnfoldTUnfoldT<TH1,TH2> RooUnfoldTUnfold;
#ifndef NOROOFIT
//! \class RooFitUnfoldTUnfold
//! \brief specialization of RooUnfoldTUnfoldT for RooAbsReal objects
typedef RooUnfoldTUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldTUnfold;
#endif

#endif
