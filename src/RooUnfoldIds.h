//=====================================================================-*-C++-*-
//! \class RooUnfoldIdsT
//! \brief Inspired by Tim Adye code for RooUnfoldSvd. For support, contact: chris.meyer@cern.ch
//! \author Bogdan Malaescu <bogdan.malaescu@cern.ch>, Christopher Meyer <chris.meyer@cern.ch>
//==============================================================================

#ifndef ROOUNFOLDIDS_HH
#define ROOUNFOLDIDS_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

#include "RooUnfoldFitHelpers.h"

#include "TRandom3.h"

template<class Hist, class Hist2D> 
  class RooUnfoldIdsT : public RooUnfoldT<Hist,Hist2D> {
  
 public:
  
  RooUnfoldIdsT(); // default constructor
  RooUnfoldIdsT(const char *name, const char *title); // named constructor
  RooUnfoldIdsT(const TString &name, const TString &title); // named constructor
  RooUnfoldIdsT(const RooUnfoldIdsT<Hist,Hist2D> &rhs); // copy constructor
  virtual ~RooUnfoldIdsT(); // destructor
  RooUnfoldIdsT<Hist,Hist2D> &operator=(const RooUnfoldIdsT<Hist,Hist2D> &rhs); // assignment operator
  //virtual RooUnfoldIds* Clone(const char *newname = NULL) const;
  
  // Special constructors

  RooUnfoldIdsT(const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist *meas, Int_t niter = 1, const char* name=0, const char* title=0);


  // Method-neutral method SetRegParm just calls SetNIter
  virtual void  SetRegParm (Double_t parm);
  virtual Double_t GetRegParm() const;
  
  void SetNIter(Int_t niter);
  Int_t GetNIter() const;
  
  void SetLambdaM(Double_t lambdaM);
  Double_t GetLambdaM() const;
  
  void SetLambdaU(Double_t lambdaU);
  Double_t GetLambdaU() const;
  
  void SetLambdaL(Double_t lambdaL);
  Double_t GetLambdaL() const;
  
  void SetLambdaS(Double_t lambdaS);
  Double_t GetLambdaS() const;
  
  virtual void Reset();

  TMatrixD* GetUnfoldCovMatrix(const TMatrixD *cov, Int_t ntoys, Int_t seed = 1) const;
  TMatrixD* GetAdetCovMatrix(Int_t ntoys, Int_t seed = 1) const;

 protected:
  void Assign(const RooUnfoldIdsT &rhs); // implementation of assignment operator
  virtual void Unfold() const override;
  virtual void GetCov() const override;

 private:
  void Init();
  void Destroy();
  void CopyData(const RooUnfoldIdsT &rhs);
  Hist* histNoOverflow(const Hist* hist, bool overflow);

  TVectorD* GetIDSUnfoldedSpectrum(TVectorD h_RecoMC, TVectorD h_TruthMC, TMatrixD h_2DSmear, TVectorD h_RecoData, TVectorD h_RecoDataError, Int_t iter) const;
  Double_t Probability(Double_t deviation, Double_t sigma, Double_t lambda) const;
  Double_t MCnormalizationCoeff(const TVectorD *vd, const TVectorD *errvd, const TVectorD *vRecmc, const Int_t dim, const Double_t estNknownd, const Double_t Nmc, const Double_t lambda, const TVectorD *soustr_ ) const;
  Double_t MCnormalizationCoeffIter(const TVectorD *vd, const TVectorD *errvd, const TVectorD *vRecmc, const Int_t dim, const Double_t estNknownd, const Double_t Nmc, const TVectorD *soustr_, Double_t lambdaN = 0., Int_t NiterMax = 5, Int_t messAct = 1) const;
  void IdsUnfold(const TVectorD &b, const TVectorD &errb, const TMatrixD &A, const Int_t dim, const Double_t lambda, TVectorD *soustr_, TVectorD *unf) const;
  void ComputeSoustrTrue(const TMatrixD *A, const TVectorD *unfres, const TVectorD *unfresErr, Int_t N, TVectorD *soustr_, Double_t lambdaS) const;
  void ModifyMatrix(TMatrixD *Am, const TMatrixD *A, const TVectorD *unfres, const TVectorD *unfresErr, Int_t N, const Double_t lambdaM_, TVectorD *soustr_, const Double_t lambdaS_) const;
  void PerformIterations(const TVectorD &data, const TVectorD &dataErr, const TMatrixD &A_, const Int_t &N_, Double_t lambdaL_, Int_t NstepsOptMin_, Double_t lambdaU_, Double_t lambdaM_, Double_t lambdaS_, TVectorD* unfres1IDS_, TVectorD* unfres2IDS_) const;
  TMatrixD* GetSqrtMatrix(const TMatrixD& covMat);
  void GenGaussRnd(TArrayD& v, const TMatrixD& sqrtMat, TRandom3& R) const;
  
 protected:
  mutable Int_t _niter;
  mutable Int_t _nb;
  
  mutable Double_t _lambdaL; // initial unfolding regularization (before folding matrix improvement)
  mutable Double_t _lambdaUmin; // regularize Unfolding
  mutable Double_t _lambdaMmin; // regularize Modification of folding matrix
  mutable Double_t _lambdaS; // regularize background Subtraction
  
 public:
  ClassDefT (RooUnfoldIdsT, 1)
};

//! \class RooUnfoldIds 
//! \brief specialization of RooUnfoldIdsT for TH1/TH2 objects
typedef RooUnfoldIdsT<TH1,TH2> RooUnfoldIds;
#ifndef NOROOFIT
//! \class RooFitUnfoldIds
//! \brief specialization of RooUnfoldIdsT for RooAbsReal objects
typedef RooUnfoldIdsT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldIds;
#endif

#endif
