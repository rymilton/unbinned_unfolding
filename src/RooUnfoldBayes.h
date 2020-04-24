//=====================================================================-*-C++-*-
//! \class RooUnfoldBayesT
//! \brief Bayesian unfolding. Just an interface to RooUnfoldBayesImpl.
//! \author Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDBAYES_HH
#define ROOUNFOLDBAYES_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

#include "TVectorD.h"
#include "TMatrixD.h"

class TH1;
class TH2;

template<class Hist, class Hist2D>
class RooUnfoldBayesT : public RooUnfoldT<Hist,Hist2D> {

public:

  // Standard methods

  RooUnfoldBayesT(); // default constructor
  RooUnfoldBayesT (const char*    name, const char*    title); // named constructor
  RooUnfoldBayesT (const TString& name, const TString& title); // named constructor
  RooUnfoldBayesT (const RooUnfoldBayesT<Hist,Hist2D>& rhs); // copy constructor
  RooUnfoldBayesT& operator= (const RooUnfoldBayesT<Hist,Hist2D>& rhs); // assignment operator

  // Special constructors

  RooUnfoldBayesT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t niter= 4, Bool_t smoothit= false,
		   Bool_t handleFakes= false, const char* name= 0, const char* title= 0);

  void SetIterations (Int_t niter= 4);
  void SetSmoothing  (Bool_t smoothit= false);
  void HandleFakes (Bool_t handleFakes);
  Int_t GetIterations() const;
  Int_t GetSmoothing()  const;
  const TMatrixD& UnfoldingMatrix() const;

  virtual void  SetRegParm (Double_t parm);
  virtual double GetRegParm() const;
  virtual void Reset();
  virtual void Print (Option_t* option= "") const;

protected:
  void Assign (const RooUnfoldBayesT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void Unfold() const override ;
  virtual void GetCov() const override ;
  virtual void GetSettings() const override;

  void setup() const;
  void unfold() const;
  void getCovariance() const;

  void smooth(TVectorD& PbarCi) const;
  double getChi2(const TVectorD& prob1,
                   const TVectorD& prob2,
                   double nevents) const;

private:
  void Init();
  void CopyData (const RooUnfoldBayesT<Hist,Hist2D>& rhs);

protected:
  // instance variables
  mutable int _niter;
  mutable int _smoothit;
  mutable int _handleFakes;

  mutable int _nc;              // number of causes  (same as _nt)
  mutable int _ne;              // number of effects (same as _nm)
  mutable double _N0C;          // number of events in prior
  mutable double _nbartrue;     // best estimate of number of true events

  mutable TVectorD _nEstj;        // Number of measured events from Effect E_j
  mutable TVectorD _nCi;          // Number of true events from cause C_i
  mutable TVectorD _nbarCi;       // Estimated number of true events from cause C_i
  mutable TVectorD _efficiencyCi; // efficiency for detecting cause C_i
  mutable TVectorD _P0C;          // prior before last iteration
  mutable TVectorD _UjInv;        // 1 / (folded prior) from last iteration

  mutable TMatrixD _Nji;          // mapping of causes to effects
  mutable TMatrixD _Mij;          // unfolding matrix
  mutable TMatrixD _Vij;          // covariance matrix
  mutable TMatrixD _VnEstij;      // covariance matrix of effects
  mutable TMatrixD _dnCidnEj;     // measurement error propagation matrix
  mutable TMatrixD _dnCidPjk;     // response error propagation matrix (stack j,k into each column)

public:
  ClassDefT (RooUnfoldBayesT, 1) // Bayesian Unfolding
};


//! \class RooUnfoldBayes 
//! \brief specialization of RooUnfoldBayesT for TH1/TH2 objects
typedef RooUnfoldBayesT<TH1,TH2> RooUnfoldBayes;
#ifndef NOROOFIT
//! \class RooFitUnfoldBayes
//! \brief specialization of RooUnfoldBayesT for RooAbsReal objects
typedef RooUnfoldBayesT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldBayes;
#endif

#endif
