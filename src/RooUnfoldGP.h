//=====================================================================-*-C++-*-
//! \class RooUnfoldGPT
//! \brief Unfolding class using Gaussian Processes. It uses the RooUnfoldInvertT class 
//!      to get an initial solution and uses a marginal likelihood minimization 
//!      to find the optimal regularization. 
//! \author Pim Verschuuren <pim.verschuuren@rhul.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDGP_H_
#define ROOUNFOLDGP_H_

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

class TDecompSVD;

template<class Hist, class Hist2D>
class RooUnfoldGPT : public RooUnfoldT<Hist,Hist2D> {

public:

  enum Kernel 
    {
      kRadial,
      kGibbs
    };

  RooUnfoldGPT(); // default constructor
  RooUnfoldGPT (const char*    name, const char*    title); // named constructor
  RooUnfoldGPT (const TString& name, const TString& title); // named constructor
  RooUnfoldGPT (const RooUnfoldGPT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldGPT(); // destructor
  RooUnfoldGPT& operator= (const RooUnfoldGPT<Hist,Hist2D>& rhs); // assignment operator
  RooUnfoldGPT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kernel = 1, const char* name=0, const char* title=0);
  virtual void SetRegParm(Double_t parm) override;
  virtual void Reset() override;
  TDecompSVD* Impl();
  virtual Double_t GetRegParm() const override;
  virtual RooUnfolding::Algorithm GetAlgorithm() const override;    

protected:
  virtual void Unfold() const override;
  virtual void GetCov() const override;
  void MLEstimator() const;
  void MLCovariance() const;
  void MAPEstimator() const;
  void MAPCovariance() const;
  Bool_t checkGP(const TVectorD& hist) const;
private:
  void Init();
  void SetBinCenters() const;
  void SetFitSettings() const;
  Bool_t InvertResponse() const;
  TVectorD MultiplyMatrixVector(const TMatrixD& matrix, const TVectorD& vector) const;
  TVectorD MultiplyVectorMatrix(const TVectorD& vector, const TMatrixD& matrix) const;
  

  Double_t MarginalLH(const double *params) const;
  void MinimizeMLH() const;
  void SetMinInit(std::vector<Double_t>& init) const;
  void SetMinStep(std::vector<Double_t>& step) const;

  void EvaluateK(const double *params) const;
  void EvaluateKstar(const double *params) const;
  void EvaluateKstarstar(const double *params) const;

  void EvaluateK(const std::vector<Double_t>&  params) const;
  void EvaluateKstar(const std::vector<Double_t>&  params) const;
  void EvaluateKstarstar(const std::vector<Double_t>&  params) const;

  void GibbsKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void RadialKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void GibbsKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;
  void RadialKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const;


  void printMatrix(const TMatrixD& matrix) const;
  void printVector(const TVectorD& vector) const;

protected:

  // cache
  class Cache {
  public:
  // instance variables
    TDecompSVD* _svd;
    TMatrixD*   _resinv;
    Bool_t      _haveMLEst;
    Bool_t      _haveMLCov;
    Bool_t      _MLHConverged;
    
    TVectorD    _MAPEst;
    TMatrixD    _MAPCov;
    
    TVectorD    _MLEst;
    TMatrixD    _MLCov;
    TMatrixD    _K;
    TMatrixD    _Kstar;
    TMatrixD    _Kstarstar;
    
    Double_t    _truMin;
    Double_t    _truMax;
    Double_t    _obsMin;
    Double_t    _obsMax;
    TVectorD    _truBinCenters;
    TVectorD    _obsBinCenters;
    
    std::vector<Double_t> _kernel_init;
    std::vector<Double_t> _kernel_step;
    std::vector<Double_t> _opt_params;
    ~Cache();
  };
  mutable Cache _specialcache;  //!
  Int_t _kernel;

public:
  ClassDefOverride (RooUnfoldGPT, 1)  // Unregularised unfolding
};

//! \class RooUnfoldGP 
//! \brief specialization of RooUnfoldGPT for TH1/TH2 objects
typedef RooUnfoldGPT<TH1,TH2> RooUnfoldGP;
#ifndef NOROOFIT
//! \class RooFitUnfoldGP
//! \brief specialization of RooUnfoldGPT for RooAbsReal objects
typedef RooUnfoldGPT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldGP;
#endif

#endif /*ROOUNFOLDGP_H_*/
