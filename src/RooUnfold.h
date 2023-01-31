//=====================================================================-*-C++-*-
//! \class RooUnfoldT
//! \brief Unfolding framework base class.
//! \author Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLD_HH
#define ROOUNFOLD_HH

#include "TNamed.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TMatrixD.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"

class TH1;

template<class Hist, class Hist2D>
class RooUnfoldT : public TNamed {
public:
  static RooUnfoldT<Hist,Hist2D>* New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Double_t regparm=-1e30,
                                       const char* name= 0, const char* title= 0);
  // typedefs to ensure compatibility with legacy code
  
  typedef RooUnfolding::Algorithm Algorithm;
  typedef RooUnfolding::ErrorTreatment ErrorTreatment;
  typedef RooUnfolding::BiasError BiasError;
  static const Algorithm kNone;
  static const Algorithm kBayes;
  static const Algorithm kSVD;
  static const Algorithm kBinByBin;
  static const Algorithm kTUnfold;
  static const Algorithm kInvert;
  //static const Algorithm kDagostini;
  static const Algorithm kIDS;
  static const Algorithm kGP;
  static const Algorithm kPoisson;
  static const ErrorTreatment kNoError;
  static const ErrorTreatment kErrors;
  static const ErrorTreatment kCovariance;
  static const ErrorTreatment kErrorsToys;
  static const ErrorTreatment kCovToys;
  static const ErrorTreatment kDefault;
  static const BiasError kBiasSD;
  static const BiasError kBiasSDM;
  static const BiasError kBiasRMS;
  
  // Standard methods

  RooUnfoldT(); // default constructor
  RooUnfoldT (const char*    name, const char*    title); // named constructor
  RooUnfoldT (const TString& name, const TString& title); // named constructor
  RooUnfoldT (const RooUnfoldT<Hist,Hist2D>& rhs); // copy constructor
  virtual ~RooUnfoldT(); // destructor
  RooUnfoldT<Hist,Hist2D>& operator= (const RooUnfoldT<Hist,Hist2D>& rhs); // assignment operator

  // Special constructors

  RooUnfoldT<Hist,Hist2D> (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name= 0, const char* title= 0);

  // Set up an existing object
  
  virtual RooUnfoldT<Hist,Hist2D>& Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas);
  virtual void SetMeasured (const Hist* meas);
  virtual void SetMeasured (const TVectorD& meas, const TMatrixD& cov);
  virtual void SetTruth (const Hist* truth);
  virtual void SetTruth (const TVectorD& truth, const TVectorD& err);
  virtual void SetBkg (const Hist* bkg);
  virtual void SetBkg (const TVectorD& bkg, const TVectorD& err);
  virtual void SetMeasured (const TVectorD& meas, const TVectorD& err);
  virtual void SetMeasuredCov (const TMatrixD& cov);
  virtual void SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership = false);
  virtual void Reset ();

  // Accessors
  virtual const RooUnfoldResponseT<Hist,Hist2D>* response() const;
  virtual RooUnfoldResponseT<Hist,Hist2D>* response();
  virtual const Hist* Hmeasured() const;
  virtual Hist* Hmeasured();
  virtual const Hist* Htruth() const;
  virtual Hist* Htruth(); 
  virtual const Hist* Hbkg() const;
  virtual Hist* Hbkg();
  virtual Hist* Hunfold (RooUnfolding::ErrorTreatment withError=RooUnfolding::kErrors);
  DEPRECATED("Hreco is considered ambigious and should no longer be used, please use Hunfold instead") virtual Hist* Hreco   () { return Hunfold(RooUnfolding::kErrors); }  

  const    TVectorD& Vmeasured() const;   // Measured distribution as a TVectorD
  const    TVectorD& Emeasured() const;   // Measured distribution errors as a TVectorD
  const    TVectorD& Vtruth() const;   // Truth distribution as a TVectorD
  const    TVectorD& Vbkg() const;   // Background distribution as a TVectorD
  const    TVectorD  Vbias() const;   // Bias distribution as a TVectorD
  //--------------------------------------------------------- !!!!!! attention Vpull() added by me D.Kavtaradze ----------------
  const    TVectorD  VAvgPull() const;   // Pull distribution as a TVectorD
  const    TMatrixD GetPull() const;    //Pull distribution as a matrix for each toy each bin
  const    TVectorD Ebias(RooUnfolding::BiasError E_type=RooUnfolding::kBiasSDM) const;   // Bias distribution errors as a TVectorD
  const    TMatrixD& GetMeasuredCov() const;   // Measured distribution covariance matrix

  virtual const TVectorD&  Vunfold() const;
  virtual TMatrixD   Eunfold  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance) const;
  virtual TVectorD   EunfoldV (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kErrors) const;
  virtual TMatrixD   Wunfold  (RooUnfolding::ErrorTreatment witherror=RooUnfolding::kCovariance) const;

  virtual TVectorD   CoverageProbV(Int_t sigma=1) const;
  virtual TVectorD   ScanCoverage(TVectorD& regparms, Int_t bin = -1, Int_t sigma=1) const;
  virtual TVectorD   ScanBias2Var(TVectorD& regparms, Int_t bin = -1) const;

  virtual Int_t      verbose() const;
  virtual void       SetVerbose (Int_t level);
  virtual void       SetOverflow(Int_t overflow);
  virtual void       IncludeSystematics (RooUnfolding::SystematicsTreatment dosys= RooUnfolding::kAll);
  inline void IncludeSystematics (int i) { return IncludeSystematics((RooUnfolding::SystematicsTreatment)i); }
  virtual Int_t      SystematicsIncluded() const;
  virtual Int_t      NToys() const;         // Number of toys
  virtual void       SetNToys (Int_t toys); // Set number of toys
  virtual Int_t      Overflow() const;
  virtual void       PrintTable (const Hist* hTrue= 0, RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault) const;  
  virtual void       PrintTable (std::ostream& o, const Hist* hTrue= 0, RooUnfolding::ErrorTreatment withError=RooUnfolding::kDefault) const;
  virtual void       SetRegParm (Double_t parm);
  virtual Double_t   GetRegParm() const; // Get Regularisation Parameter
  Double_t Chi2 (const Hist* hTrue,RooUnfolding::ErrorTreatment DoChi2=RooUnfolding::kCovariance) const;
  virtual void CalculateBias(Int_t ntoys = 50, const Hist* hTrue = 0) const; // Estimate bias by throwing toys.

  virtual RooUnfolding::Algorithm GetAlgorithm() const;
  Double_t GetMinParm() const;
  Double_t GetMaxParm() const;
  Double_t GetStepSizeParm() const;
  Double_t GetDefaultParm() const;
  double RunToy(TVectorD& x, TVectorD& xe) const;
  void RunRooFitToys(int ntoys, std::vector<TVectorD>& vx, std::vector<TVectorD>& vxe, std::vector<double>& chi2) const;
  void RunToys(int ntoys, std::vector<TVectorD>& vx, std::vector<TVectorD>& vxe, std::vector<double>& chi2) const;
  
  void Print(Option_t* opt="") const override;
  void Dump() const override;
  void ForceRecalculation() const ;

protected:
  void Assign (const RooUnfoldT<Hist,Hist2D>& rhs); // implementation of assignment operator
  virtual void SetNameTitleDefault();
  virtual void Unfold() const;
  virtual Bool_t UnfoldWithErrors (RooUnfolding::ErrorTreatment withError, bool getWeights=false) const;
  virtual void ResetUnfold();

  static TMatrixD CutZeros     (const TMatrixD& ereco);
  static Int_t    InvertMatrix (const TMatrixD& mat, TMatrixD& inv, const char* name="matrix", Int_t verbose=0);

private:
  void Init();
  void CopyData (const RooUnfoldT<Hist,Hist2D>& rhs);
  virtual void ClearCache() const;

protected:
  // cache 
  virtual void GetErrors() const; // Get 
  virtual void GetCov() const; // Get covariance matrix using errors on measured distribution
  void GetErrorsToys() const;
  void GetCovToys() const;

  void GetSampleVar(std::vector<TVectorD>& munfolded) const;
  void GetSampleCov(std::vector<TVectorD>& munfolded) const;

  virtual void GetSettings() const;
  virtual void GetErrMat() const; // Get covariance matrix using errors from residuals on reconstructed distribution
  virtual void GetWgt() const; // Get weight matrix using errors on measured distribution

  class Cache {
  public:
    Cache();
    ~Cache();
    typename RooUnfoldT<Hist,Hist2D>::Cache& operator=(const Cache&);
    Double_t _minparm;       // Minimum value to be used in RooUnfoldParms
    Double_t _maxparm;       // Maximum value to be used in RooUnfoldParms
    Double_t _stepsizeparm;  // StepSize value to be used in RooUnfoldParms
    Double_t _defaultparm;   // Recommended value for regularisation parameter
    Bool_t   _unfolded;      // unfolding done
    Bool_t   _fail;          // unfolding failed
    Bool_t   _haveCov;       // have _cov
    Bool_t   _haveWgt;       // have _wgt
    Bool_t   _have_err_mat;  // have _err_mat
    Bool_t   _haveErrors;    // have _variances
    Bool_t   _haveBias;      // have _bias
    //!!!!! attention next two lines added by me D.Kavtaradze
    Bool_t   _havePull;      // have _pull 
    Bool_t   _haveAvgPull;   //have _AvgPull
    // ---------------------------------------------end of addition
    TVectorD _bias;          // Estimated bias on each truth bin
    // !!!!! attention TVectorD _AvgPull; added by me D.Kavtaradze
    TVectorD _AvgPull;       // Estimated pull on each truth bin with average method
    TVectorD _sdbias;        // SD of the bias on each truth bin
    TVectorD _sdmbias;       // SD of the mean of the bias on each truth bin
    TVectorD _rmsbias;       // Root mean squared on each bin
    TVectorD _rec;           // Reconstructed distribution
    TMatrixD _cov;           // Reconstructed distribution covariance
    TMatrixD _wgt;           // Reconstructed distribution weights (inverse of _cov)
    TVectorD _variances;     // Error matrix diagonals
    TMatrixD _err_mat;       // Error matrix from toys
    TMatrixD _pull;          //Pull matrix for each toy each bin
    TVectorD* _vMes;         // Cached measured vector
    TVectorD* _eMes;         // Cached measured error
    TVectorD* _vTruth;       // Cached truth vector
    TVectorD* _vBkg;         // Cached bkg vector
    TMatrixD* _covL;         // Cached lower triangular matrix for which _covMes = _covL * _covL^T.
    TMatrixD* _covMes;       // Measurement covariance matrix    
  };
  mutable TRandom rnd; //!
  mutable Cache _cache; //!
  mutable RooUnfolding::ErrorTreatment _withError = RooUnfolding::kDefault; // type of error last calulcated
  
  TMatrixD* _covMes;                         // Measurement covariance matrix
  Int_t    _verbose;                         // Debug print level
  Int_t    _nm;                              // Total number of measured bins (including under/overflows if _overflow set)
  Int_t    _nt;                              // Total number of truth    bins (including under/overflows if _overflow set
  Int_t    _overflow;                        // Use histogram under/overflows if 1 (set from RooUnfoldResponse)
  Int_t    _NToys;                           // Number of toys to be used
  RooUnfolding::SystematicsTreatment _dosys; // include systematic errors from response matrix? use _dosys=2 to exclude measurement errors
  RooUnfoldResponseT<Hist,Hist2D>* _res;     // Response matrix (not owned)
  Hist*    _meas;                            // Measured distribution (not owned)
  Hist*    _bkg;                             // Estimated reconstructed background distribution (not owned)
  Hist*    _truth;                           // Estimated truth distribution. Used in some regularization schemes. (not owned)

public:

  ClassDefOverride (RooUnfoldT, 2) // Unfolding base class: implementations in RooUnfoldBayes, RooUnfoldSvd, RooUnfoldBinByBin, RooUnfoldTUnfold, RooUnfoldInvert, RooUnfoldIds
};

//! \class RooUnfoldBinByBin 
//! \brief specialization of RooUnfoldBinByBinT for TH1/TH2 objects
typedef RooUnfoldT<TH1,TH2> RooUnfold;
#endif
