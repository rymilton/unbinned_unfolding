//=====================================================================-*-C++-*-
//! \class RooUnfoldSvdT
//! \brief SVD unfolding. Just an interface to SVDUnfold.
//! \author Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDSVD_HH
#define ROOUNFOLDSVD_HH

#include "RooUnfold.h"
#include "RooUnfoldResponse.h"

template<class Hist, class Hist2D>
class RooUnfoldSvdT : public RooUnfoldT<Hist,Hist2D> {
public:
  class SVDUnfold {

public:

   SVDUnfold( const Hist* bdat, const Hist* bini, const Hist* xini, const Hist2D* Adet );
   SVDUnfold( const TVectorD& bdat, const TMatrixD& Bcov, const TVectorD& bini, const TVectorD& xini, const TMatrixD& Mdet, const TMatrixD& MdetE );
   SVDUnfold( const Hist* bdat, const TMatrixD& Bcov, const Hist* bini, const Hist* xini, const Hist2D* Adet );
   SVDUnfold( const Hist* bdat, const TMatrixD& Bcov, const Hist* bini, const Hist* xini, const TMatrixD& Mdet, const TMatrixD& MdetE );    
   SVDUnfold( const SVDUnfold& other );

   // Destructor
   virtual ~SVDUnfold(); 

   // Set option to normalize unfolded spectrum to unit area
   // "normalize" - switch 
   void     SetNormalize ( Bool_t normalize ) { fNormalize = normalize; }

   // Do the unfolding
   // "kreg"   - number of singular values used (regularisation)
   //Hist*    Unfold       ( Int_t kreg );
   TVectorD    UnfoldV       ( Int_t kreg );    

   // Determine for given input error matrix covariance matrix of unfolded 
   // spectrum from toy simulation
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   TMatrixD    GetUnfoldCovMatrix( const TMatrixD& cov, Int_t ntoys, Int_t seed = 1 );

   // Determine covariance matrix of unfolded spectrum from finite statistics in 
   // response matrix
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // "uncmat" - matrix containing the uncertainty on the detector matrix elements if different from purely statistical without any weights
   TMatrixD    GetAdetCovMatrix( Int_t ntoys, Int_t seed=1, const TMatrixD* uncmat=0 );

   // Regularisation parameter
   Int_t    GetKReg() const { return fKReg; }

   // Obtain the distribution of |d| (for determining the regularization)
   Hist*    GetD() const;

   // Obtain the distribution of singular values
   TVectorD    GetSV() const;

   // Obtain the computed regularized covariance matrix
   const TMatrixD&    GetXtau() const;

   // Obtain the computed inverse of the covariance matrix
   const TMatrixD&    GetXinv() const;
   
   //Obtain the covariance matrix on the data
   const TMatrixD& GetBCov() const;

private: 

    // Helper functions for vector and matrix operations
   void            FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const;
   static Double_t GetCurvature       ( const TVectorD& vec, const TMatrixD& curv );

   // Helper functions
   static TMatrixD MatDivVec( const TMatrixD& mat, const TVectorD& vec, Int_t zero=0 );
   static TVectorD CompProd ( const TVectorD& vec1, const TVectorD& vec2 );

   static TVectorD VecDiv                 ( const TVectorD& vec1, const TVectorD& vec2, Int_t zero = 0 );
   static void     RegularisedSymMatInvert( TMatrixDSym& mat, Double_t eps = 1e-3 );
   
   // Class members
   Int_t       fMdim;        //! Reconstructed dimensions
   Int_t       fTdim;        //! Truth dimensions
   Int_t       fDdim;        //! Derivative for curvature matrix
   Bool_t      fNormalize;   //! Normalize unfolded spectrum to 1
   Int_t       fKReg;        //! Regularisation parameter
   Hist*       fDHist;       //! Distribution of d (for checking regularization)
   TVectorD       fSVHist;      //! Distribution of singular values
   TMatrixD       fXtau;        //! Computed regularized covariance matrix
   TMatrixD       fXinv;        //! Computed inverse of covariance matrix

   // Input histos

   const TVectorD fBdat;        // measured distribution (data)
   TMatrixD fBcov;        // covariance matrix of measured distribution (data)
   const TVectorD fBini;        // reconstructed distribution (MC)
   const TVectorD fXini;        // truth distribution (MC)
   TMatrixD fAdet;        // Detector response matrix
   TMatrixD fAdetE;      

   // Evaluation of covariance matrices
   TVectorD       fToyhisto;    //! Toy MC histogram
   TVectorD       fToyhistoE;    //! Toy MC histogram
   TMatrixD       fToymat;      //! Toy MC detector response matrix
   TMatrixD       fToymatE;      //! Toy MC detector response matrix
   Bool_t      fToyMode;     //! Internal switch for covariance matrix propagation
   Bool_t      fMatToyMode;  //! Internal switch for evaluation of statistical uncertainties from response matrix
};

public:

  // Standard methods

  RooUnfoldSvdT(); // default constructor
  RooUnfoldSvdT (const char*    name, const char*    title); // named constructor
  RooUnfoldSvdT (const TString& name, const TString& title); // named constructor
  RooUnfoldSvdT (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldSvdT(); // destructor
  RooUnfoldSvdT<Hist,Hist2D>& operator= (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // assignment operator

  // Special constructors

  RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg= 0,
                const char* name= 0, const char* title= 0);
  // compatibility constructor
  RooUnfoldSvdT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t kreg, Int_t ntoyssvd,
                const char* name= 0, const char* title= 0);

  void SetKterm (Int_t kreg);
  Int_t GetKterm() const;
  virtual void  SetRegParm (Double_t parm);
  virtual Double_t GetRegParm() const;
  virtual void Reset();
  SVDUnfold* Impl();

protected:
  void Assign (const RooUnfoldSvdT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void Unfold() const override;
  virtual void GetCov() const override;
  virtual void GetWgt() const override;
  virtual void GetSettings() const override;

private:
  void Init();
  void Destroy();
  void CopyData (const RooUnfoldSvdT<Hist,Hist2D>& rhs);
  void PrepareHistograms() const;
  
protected:
  // instance variables
  mutable SVDUnfold* _svd;  //! Implementation in TSVDUnfold object (no streamer)
  mutable Int_t _kreg;

  mutable const Hist *_meas1d, *_train1d, *_truth1d;
  mutable const Hist2D *_reshist, *_meascov;

public:
  ClassDefT (RooUnfoldSvdT, 1) // SVD Unfolding (interface to TSVDUnfold)
};

//! \class RooUnfoldSvd 
//! \brief specialization of RooUnfoldSvdT for TH1/TH2 objects
typedef RooUnfoldSvdT<TH1,TH2> RooUnfoldSvd;
#ifndef NOROOFIT
//! \class RooFitUnfoldSvd
//! \brief specialization of RooUnfoldSvdT for RooAbsReal objects
typedef RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldSvd;
#endif

#endif
