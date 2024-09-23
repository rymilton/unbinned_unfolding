//=====================================================================-*-C++-*-
//! \class RooUnfoldResponseT
//! \brief Response Matrix
//! \author Tim Adye <T.J.Adye@rl.ac.uk>
//==============================================================================

#ifndef ROOUNFOLDRESPONSE_HH
#define ROOUNFOLDRESPONSE_HH

#include "TNamed.h"
#include "TMatrixD.h"
#include "RooUnfoldFitHelpers.h"
#include "RooRealVar.h"
#include "TRandom.h"
#include "RooAbsData.h"
#include "TH1.h"
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,0,0)
#include "TVectorDfwd.h"
#else
class TVectorD;
#endif
class TF1;
class TH2;
class TH2D;
class TCollection;

template<class Hist, class Hist2D>
class RooUnfoldResponseT : public TNamed {

public:
  // Standard methods

  RooUnfoldResponseT(); // default constructor
  virtual ~RooUnfoldResponseT(); // destructor
  RooUnfoldResponseT(const char* name, const char* title); // named constructor
  RooUnfoldResponseT(const TString& name, const TString& title); // named constructor
  RooUnfoldResponseT(const RooUnfoldResponseT<Hist,Hist2D>& rhs); // copy constructor
  RooUnfoldResponseT(const char* name, const char* title, Hist2D* response, Hist* truth, Hist* reco, bool overflow, bool density);
  // Accessors

  Int_t        GetDimensionMeasured() const;   // Dimensionality of the measured distribution
  Int_t        GetDimensionTruth()    const;   // Dimensionality of the truth distribution
  Int_t        GetNbinsMeasured()     const;   // Total number of bins in the measured distribution
  Int_t        GetNbinsTruth()        const;   // Total number of bins in the truth distribution

  const Hist*   Hmeasured()            const;   // Measured distribution, including fakes
  Hist*         Hmeasured();                    // Measured distribution, including fakes
  const Hist*   Hfakes()               const;   // Fakes distribution
  Hist*         Hfakes();                       // Fakes distribution
  const Hist*   Htruth()               const;   // Truth distribution, used for normalisation
  Hist*         Htruth();                       // Truth distribution, used for normalisation
  const Hist2D* Hresponse()            const;   // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  Hist2D*       Hresponse();                    // Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  Hist2D*       HresponseNoOverflow()  const;   // Response matrix with under/overflow bins moved into histogram body

  TVectorD        Vfolded(const TVectorD& truth) const;  // Fold a vector according to the response
  const TVectorD& Vpurity()           const;   // Distribution with the purities as a TVectorD
  const TVectorD& Vefficiency()       const;   // Distribution with the efficiencies as a TVectorD
  const TVectorD& Vmeasured()         const;   // Measured distribution as a TVectorD
  const TVectorD& Emeasured()         const;   // Measured distribution errors as a TVectorD
  const TVectorD& Vfakes()            const;   // Fakes distribution as a TVectorD
  const TVectorD& Vtruth()            const;   // Truth distribution as a TVectorD
  const TVectorD& Etruth()            const;   // Truth distribution errors as a TVectorD
  const TMatrixD& Mresponse(bool norm = true)         const;   // Response matrix as a TMa
  const TMatrixD& Mresponse(const Hist* truth)         const;   // Response matrix as a TMatrixD: (row,column)=(measured,truth)
  const TMatrixD& Eresponse(bool norm = true)         const;   // Response matrix errors as a TMatrixD: (row,column)=(measured,truth)

  Double_t operator() (Int_t r, Int_t t) const;// Response matrix element (measured,truth)

  void   UseOverflow (Bool_t set= kTRUE);      // Specify to use overflow bins
  Bool_t UseOverflowStatus() const;            // Get UseOverflow setting
  Bool_t UseDensityStatus() const;            // Get UseDensity setting
  bool HasFakes() const;                // Return number of bins with fakes
  virtual void Print (Option_t* option="") const override;

  Hist* ApplyToTruth (const Hist* truth= 0, const char* name= "AppliedResponse") const; // If argument is 0, applies itself to its own truth
  TF1* MakeFoldingFunction (TF1* func, Double_t eps=1e-12, Bool_t verbose=false) const;

  virtual void RunToy() const;
  virtual void ClearCache() const;

protected:

  virtual void setup();
  virtual bool Cached() const;
  virtual void SetNameTitleDefault (const char* defname= 0, const char* deftitle= 0);

  // instance variables

  Hist*  _mes = 0;      // Measured histogram
  Hist*  _fak = 0;      // Fakes    histogram
  Hist*  _tru = 0;      // Truth    histogra
  Hist2D*_res = 0;      // Response histogram
  Int_t _overflow = 0; // Use histogram under/overflows if 1
  bool  _density = false;

private:
  class Cache {
  public:
    Cache();
    ~Cache();
    typename RooUnfoldResponseT<Hist,Hist2D>::Cache& operator=(const Cache&);
    TVectorD* _vMes;   //! Cached measured vector
    TVectorD* _eMes;   //! Cached measured error
    TVectorD* _vFak;   //! Cached fakes    vector
    TVectorD* _vTru;   //! Cached truth    vector
    TVectorD* _eTru;   //! Cached truth    error
    TMatrixD* _mRes;   //! Cached response matrix
    TMatrixD* _eRes;   //! Cached response error
    TMatrixD* _mResNorm;   //! Cached normalized response matrix
    TMatrixD* _eResNorm;   //! Cached normalized response error
    TVectorD* _vPur;   //! Cached purity vector
    TVectorD* _vEff;   //! Cached efficiency vector        
  };
  mutable Cache _cache; //!
  mutable TRandom rnd; //!
  mutable Bool_t    _cached = false; //! We are using cached vectors/matrices

public:

  ClassDefOverride (RooUnfoldResponseT, 1) // Respose Matrix
};

//=====================================================================-*-C++-*-
//! \class RooUnfoldResponse
//! \brief specialization of RooUnfoldResponseT for TH1/TH2 objects
//==============================================================================
class RooUnfoldResponse : public RooUnfoldResponseT<TH1,TH2> {
public:

  RooUnfoldResponse(){}; // default constructor
  virtual ~RooUnfoldResponse(){}; // destructor
  RooUnfoldResponse(const char* name, const char* title) : RooUnfoldResponseT(name,title) {}; // named constructor
  RooUnfoldResponse(const TString& name, const TString& title) : RooUnfoldResponseT(name,title) {}; // named constructor

  // Special constructors
  
  RooUnfoldResponse(Int_t nb, Double_t xlo, Double_t xhi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case with same binning, measured vs truth
  RooUnfoldResponse(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi, const char* name= 0, const char* title= 0);  // constructor -  simple 1D case
  RooUnfoldResponse(const TH1* measured, const TH1* truth, const char* name= 0, const char* title= 0, bool overflow=false);  // constructor - measured and truth only used for shape
  RooUnfoldResponse(const TH1* measured, const TH1* truth, const TH2* response, const char* name= 0, const char* title= 0, bool overflow=false);  // create from already-filled histograms
  RooUnfoldResponse(const RooUnfoldResponse& rhs); // copy constructor
  virtual RooUnfoldResponse& operator= (const RooUnfoldResponse& rhs); // assignment operator

  // Set up an existing object

  virtual RooUnfoldResponse& Reset ();  // clear an existing object
  virtual RooUnfoldResponse& Setup (const RooUnfoldResponse& rhs);  // set up based on another instance
  virtual RooUnfoldResponse& Setup (Int_t nb, Double_t xlo, Double_t xhi);  // set up simple 1D case with same binning, measured vs truth
  virtual RooUnfoldResponse& Setup (Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi);  // set up simple 1D case
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth);  // set up - measured and truth only used for shape
  virtual RooUnfoldResponse& Setup (const TH1* measured, const TH1* truth, const TH2* response);  // set up from already-filled histograms

  // Fill with training data

  static Int_t   FindBin(const TH1*  h, Double_t x);  // return vector index for bin containing (x)
  static Int_t   FindBin(const TH1*  h, Double_t x, Double_t y);  // return vector index for bin containing (x,y)
  static Int_t   FindBin(const TH1*  h, Double_t x, Double_t y, Double_t z);  // return vector index for bin containing (x,y,z)

  virtual Int_t Fill (Double_t xr, Double_t xt, Double_t w= 1.0);  // Fill 1D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w= 1.0);  // Fill 2D Response Matrix
  virtual Int_t Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w= 1.0);  // Fill 3D Response Matrix

          Int_t Miss (Double_t xt);  // Fill missed event into 1D Response Matrix
          Int_t Miss (Double_t xt, Double_t w);  // Fill missed event into 1D (with weight) or 2D Response Matrix
          Int_t Miss (Double_t xt, Double_t yt, Double_t w);  // Fill missed event into 2D (with weight) or 3D Response Matrix
  virtual Int_t Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w);  // Fill missed event into 3D Response Matrix

          Int_t Fake (Double_t xr);  // Fill fake event into 1D Response Matrix
          Int_t Fake (Double_t xr, Double_t w);  // Fill fake event into 1D (with weight) or 2D Response Matrix
          Int_t Fake (Double_t xr, Double_t yr, Double_t w);  // Fill fake event into 2D (with weight) or 3D Response Matrix
  virtual Int_t Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w);  // Fill fake event into 3D Response Matrix

  virtual void Add (const RooUnfoldResponse& rhs);
  virtual Long64_t Merge (TCollection* others);
  virtual Long64_t FakeEntries() const;                // Return number of bins with fakes

private:
  virtual Int_t Miss1D (Double_t xt, Double_t w= 1.0);  // Fill missed event into 1D Response Matrix (with weight)
  virtual Int_t Miss2D (Double_t xt, Double_t yt, Double_t w= 1.0);  // Fill missed event into 2D Response Matrix (with weight)
  virtual Int_t Fake1D (Double_t xr, Double_t w= 1.0);  // Fill fake event into 1D Response Matrix (with weight)
  virtual Int_t Fake2D (Double_t xr, Double_t yr, Double_t w= 1.0);  // Fill fake event into 2D Response Matrix (with weight)

  ClassDef (RooUnfoldResponse, 1) // Respose Matrix
};


#ifndef NOROOFIT

class RooHistFunc;
class RooHistPdf;

//=====================================================================-*-C++-*-
//! \class RooFitUnfoldResponse
//! \brief specialization of RooUnfoldResponseT for RooAbsReal objects
//==============================================================================
class RooFitUnfoldResponse : public RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> {
public:

  RooFitUnfoldResponse(){}; // default constructor
  RooFitUnfoldResponse(const RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* res); // copy constructor
  virtual ~RooFitUnfoldResponse(){}; // destructor
  RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, const RooAbsCollection* observables, bool density = false);
  RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, RooRealVar* obs_truth, RooRealVar* obs_reco, bool density = false);  
  RooFitUnfoldResponse(const char* name, const char* title, RooUnfolding::RooFitHist* response, RooUnfolding::RooFitHist* truth, RooUnfolding::RooFitHist* reco, bool density = false);

  RooHistFunc* makeHistFunc(RooDataHist* dhist);  
  RooHistFunc* makeHistFuncTruth(const TH1* hist);
  RooHistFunc* makeHistFuncMeasured(const TH1* hist);
  RooAbsPdf* makeHistPdf(RooDataHist* dhist);  
  RooAbsPdf* makeHistPdfTruth(const TH1* hist);
  RooAbsPdf* makeHistPdfMeasured(const TH1* hist);
  RooUnfolding::RooFitHist* makeHistMeasured(const TH1* hist);
  RooUnfolding::RooFitHist* makeHistTruth(const TH1* hist);
  RooUnfolding::RooFitHist* makeHist(RooAbsReal* object);
  RooUnfolding::RooFitHist* makeHist(RooDataHist* object);
  RooUnfolding::RooFitHist* makeHistSum(RooAbsReal* a, RooAbsReal* b, double ca, double cb);  
  
  ClassDef (RooFitUnfoldResponse, 1) // Respose Matrix
  
};

#endif

#endif
