/*! \class RooUnfoldT
\brief A base class for several unfolding methods.
<p>The unfolding method can either use the constructors for individual unfolding algorithms or the New() method, specifiying the algorithm to be used.
<p>The resultant distribution can be displayed as a plot (Hunfold) or as a bin by bin breakdown of the true, measured and unfolded values (PrintTable)
<p>A covariance matrix can be returned using the Eunfold() method. A vector of its diagonals can be returned with the EunfoldV() method.
<p>A summary of the unfolding algorithms which inherit from this class is below:
<ul>
<li>RooUnfoldBayes: Uses the Bayes method of unfolding based on the method written by D'Agostini (<a href="http://www.slac.stanford.edu/spires/find/hep/www?j=NUIMA,A362,487">NIM A 362 (1995) 487</a>).
<ul>
<li>Works for 1 & 2 dimensional distributions
<li>Returned errors can be either as a diagonal matrix or as a full matrix of covariances
<li>Regularisation parameter sets the number of iterations used in the unfolding (default=4)
<li>Is able to account for bin migration and smearing
<li>Can unfold if test and measured distributions have different binning.
<li>Returns covariance matrices with conditions approximately that of the machine precision. This occasionally leads to very large chi squared values
</ul>
<li> RooUnfoldSVD: Uses the singular value decomposition method of Hocker and Kartvelishvili (<a href="http://arxiv.org/abs/hep-ph/9509307">NIM A 372 (1996) 469</a>)
<ul>
<li>Regularisation parameter defines the level at which values are deemed to be due to statistical fluctuations and are cut out. (Default= number of bins/2)
<li>Returns errors as a full matrix of covariances
<li>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
<li>Can account for both smearing and biasing
<li>Returns near singular covariance matrices, again leading to very large chi squared values
</ul>
<li> RooUnfoldIds: Uses the Bayes method of unfolding based on the method written by Malaescu (<a href="http://arxiv.org/abs/1106.3107">CERN-PH-EP-2011-111</a>)
<ul>
<li>Set the number of iterations used to improve the folding matrix
<li>Regularisation parameters define the level at which values are deemed to be due to statistical fluctuations. Used for modifying the folding matrix, as well as unfolding.
<li>Returns errors as a full matrix of covariances
<li>Error processing is much the same as with the kCovToy setting with 1000 toys. This is quite slow but can be switched off.
<li>Can handle 2 dimensional distributions
<li>True and measured distributions must have the same binning
<li>Can account for both smearing and biasing
</ul>
<li> RooUnfoldBinByBin: Unfolds using the method of correction factors.
<ul>
<li>Returns errors as a diagonal matrix.
<li>Is not able to handle bin migration caused by bias/smearing of the distribution
<li>Can only handle 1 dimensional distributions
<li>True and measured distributions must have the same binning
</ul>
<li> RooUnfoldTUnfold: Uses the unfolding method implemented in ROOT's <a href="http://root.cern.ch/root/html/TUnfold.html">TUnfold</a> class
<ul>
<li>Only included in ROOT versions 5.22 and higher
<li>Only able to unfold 1 dimensional distributions
<li>Can account for bin migration and smearing
<li>Errors come as a full covariance matrix.
<li>Will sometimes warn of "unlinked" bins. These are bins with 0 entries and do not effect the results of the unfolding
<li>Regularisation parameter can be either optimised internally by plotting log10(chi2 squared) against log10(tau). The 'kink' in this curve is deemed the optimum tau value. This value can also be set manually (FixTau)
<li>The latest version (TUnfold v15) requires that RooUnfoldResponse::SetOverflow=0. ROOT versions 5.26 or below use v13 and so should be safe to use overflows
</ul>
<li> RooUnfoldInvert: The simplest method of unfolding works by simply inverting the response matrix.
<ul>
<li>For small statistics, this method does not produce useful results.
<li>The inversion method is included largely to illustrate the necessity of a more effective method of unfolding</ul>
</ul>
 */

#include "RooUnfold.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>
#include <math.h>

#include "TClass.h"
#include "TMatrixD.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TRandom.h"
#include "TMath.h"
#include "Math/ProbFunc.h"
#include "RooRandom.h"
#include "RooMultiVarGaussian.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldErrors.h"
// Need subclasses just for RooUnfold::New()
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldBinByBin.h"
#include "RooUnfoldGP.h"
#include "RooUnfoldPoisson.h"
#ifndef NOTUNFOLD
#include "RooUnfoldTUnfold.h"
#endif
#ifdef HAVE_DAGOSTINI
#include "RooUnfoldDagostini.h"
#endif
#include "RooUnfoldIds.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::sqrt;
using std::fabs;

template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kNone = RooUnfolding::kNone;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kBayes = RooUnfolding::kBayes;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kSVD = RooUnfolding::kSVD;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kBinByBin = RooUnfolding::kBinByBin;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kTUnfold = RooUnfolding::kTUnfold;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kInvert = RooUnfolding::kInvert;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kDagostini = RooUnfolding::kDagostini;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kIDS = RooUnfolding::kIDS;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kGP = RooUnfolding::kGP; 
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::Algorithm RooUnfoldT<Hist,Hist2D>::kPoisson = RooUnfolding::kPoisson;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kNoError = RooUnfolding::kNoError;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kErrors = RooUnfolding::kErrors;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kCovariance = RooUnfolding::kCovariance;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kErrorsToys = RooUnfolding::kErrorsToys;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kCovToys = RooUnfolding::kCovToys;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::ErrorTreatment RooUnfoldT<Hist,Hist2D>::kDefault = RooUnfolding::kDefault;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::BiasError RooUnfoldT<Hist,Hist2D>::kBiasSD = RooUnfolding::kBiasSD;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::BiasError RooUnfoldT<Hist,Hist2D>::kBiasSDM = RooUnfolding::kBiasSDM;
template<class Hist,class Hist2D> const typename RooUnfoldT<Hist,Hist2D>::BiasError RooUnfoldT<Hist,Hist2D>::kBiasRMS = RooUnfolding::kBiasRMS;

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, const char* name, const char* title)
  : TNamed (name, title)
{

  //! Constructor with response matrix object and measured unfolding input histogram.
  //! Should not normally be used directly - instead, create an instance of one of RooUnfold's subclasses,
  //! or use the New() static constructor.
  Init();
  Setup (res, meas);
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>*
RooUnfoldT<Hist,Hist2D>::New (RooUnfolding::Algorithm alg, const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,Double_t regparm,
                           const char* name, const char* title)
{

    /*!Unfolds according to the value of the alg enum:
    0 = kNone:     dummy unfolding
    1 = kBayes:    Unfold via iterative application of Bayes theorem
    2 = kSVD:      Unfold using singlar value decomposition (SVD)
    3 = kBinByBin: Unfold bin by bin.
    4 = kTUnfold:  Unfold with TUnfold
    5 = kInvert:   Unfold using inversion of response matrix
    7 = kIDS:      Unfold using iterative dynamically stabilized (IDS) method
    8 = kGP:       Unfold using Gaussian Processes(GP)
    9 = kPoisson:  Unfold using Poisson-based likelihood and Tikhinov regularization
    */
  
  RooUnfoldT<Hist,Hist2D>* unfold(NULL);

  switch(alg) {

  case kNone:
    unfold= new RooUnfoldT<Hist,Hist2D>         (res, meas);
    break;

  case kBayes:
    unfold= new RooUnfoldBayesT<Hist,Hist2D>    (res, meas);
    break;

  case kSVD:
    unfold= new RooUnfoldSvdT<Hist,Hist2D>      (res, meas);
    break;

  case kBinByBin:
    unfold= new RooUnfoldBinByBinT<Hist,Hist2D> (res, meas);
    break;

  case kTUnfold:
#ifndef NOTUNFOLD
    unfold= new RooUnfoldTUnfoldT<Hist,Hist2D> (res,meas);
    break;
#else
    cerr << "TUnfold library is not available" << endl;
    return 0;
#endif
    
  case kInvert:
    unfold = new RooUnfoldInvertT<Hist,Hist2D>  (res,meas);
    break;
  case kGP:
    unfold = new RooUnfoldGPT<Hist,Hist2D> (res,meas);
    break;
  case kPoisson:
    unfold = new RooUnfoldPoissonT<Hist,Hist2D> (res,meas);
    break;    
  case kDagostini:
    cerr << "RooUnfoldDagostini is not available" << endl;
    return 0;
  
  case kIDS:
    unfold= new RooUnfoldIdsT<Hist,Hist2D>      (res, meas,4);
    break;

  default: 
    cerr << "Unknown RooUnfold method " << Int_t(alg) << endl;
    return 0;
  }

  if (name)  unfold->SetName  (name);
  if (title) unfold->SetTitle (title);
  unfold->SetAlgorithm(alg);
  if (regparm != -1e30){
    unfold->SetRegParm(regparm);
  }

  return unfold;
}


template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::Cache::Cache() :
  _minparm(0),
  _maxparm(0),
  _stepsizeparm(0),
  _defaultparm(0),
  _unfolded(false),
  _fail(false),
  _haveCov(false),
  _haveWgt(false),
  _have_err_mat(false),
  _haveBias(false),
  _haveErrors(false),
  _rec(1),
  _cov(1,1),
  _wgt(1,1),
  _variances(1),
  _err_mat(1,1),
  _bias(1),
  _sdbias(1),
  _sdmbias(1),
  _rmsbias(1),
  _vMes(0),
  _eMes(0),
  _vTruth(0),
  _vBkg(0),
  _covL(0),
  _covMes(0)
{
  //! default constructor
}

template<class Hist,class Hist2D> 
typename RooUnfoldT<Hist,Hist2D>::Cache& RooUnfoldT<Hist,Hist2D>::Cache::operator= ( const RooUnfoldT<Hist,Hist2D>::Cache & other ){
  //! assignment operator
  _minparm = other._minparm;
  _maxparm = other._maxparm;
  _stepsizeparm = other._stepsizeparm;
  _defaultparm = other._defaultparm;
  _unfolded = other._unfolded;
  _haveCov = other._haveCov;
  _fail = other._fail;
  _have_err_mat = other._have_err_mat;
  _haveBias =  other._haveBias;
  _haveErrors = other._haveErrors;
  _haveWgt = other._haveWgt;
  _rec.ResizeTo(other._rec);
  _cov.ResizeTo(other._cov);
  _wgt.ResizeTo(other._wgt);
  _variances.ResizeTo(other._variances);
  _err_mat.ResizeTo(other._err_mat);
  _rec = other._rec;
  _cov = other._cov;
  _wgt = other._wgt;
  _variances = other._variances;
  _err_mat = other._err_mat;
  _bias.ResizeTo(other._bias);
  _bias = other._bias;
  _sdbias.ResizeTo(other._sdbias);  
  _sdbias = other._sdbias;
  _sdmbias.ResizeTo(other._sdmbias);  
  _sdmbias = other._sdmbias;
  _rmsbias.ResizeTo(other._rmsbias);  
  _rmsbias = other._rmsbias;
  _vMes = other._vMes;
  _eMes = other._eMes;
  _vTruth = other._vTruth;
  _vBkg = other._vBkg;
  _covL = other._covL;
  _covMes = other._covMes;
  return *this;
}


template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::Cache::~Cache(){
  //! destructor
  delete this->_vMes;
  delete this->_eMes;
  delete this->_vTruth;
  delete this->_vBkg;
  delete this->_covMes;
  delete this->_covL;
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const RooUnfoldT<Hist,Hist2D>& rhs)
  : TNamed (rhs.GetName(), rhs.GetTitle())
{
  //! Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Assign (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! assign data from another unfolding object
  if (this == &rhs) return;
  Reset();
  SetNameTitle (rhs.GetName(), rhs.GetTitle());
  CopyData (rhs);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::CopyData (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! copy data from another unfolding object
  Setup (new RooUnfoldResponseT<Hist,Hist2D>(*(rhs.response())), clone(rhs.Hmeasured()));
  SetVerbose (rhs.verbose());
  SetNToys   (rhs.NToys());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetAlgorithm (RooUnfolding::Algorithm alg)
{
  //! set the unfolding algorithm to be used
  _alg = alg;
}

template<class Hist,class Hist2D> RooUnfolding::Algorithm
RooUnfoldT<Hist,Hist2D>::GetAlgorithm () const
{
  //! return the unfolding algorithm used
  return _alg;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Reset()
{
  //! clear and reinitialize
  ClearCache();
  Init();
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Init()
{
  //! initialize an object with zero
  _res= 0;
  _meas= 0;
  _bkg= 0;
  _truth= 0;
  _nm= _nt= 0;
  _verbose= 1;
  _overflow= 0;
  _dosys= kNoSystematics;
  _covMes= 0;
  _NToys=50;
  GetSettings();
}

template<class Hist,class Hist2D> RooUnfoldT<Hist,Hist2D>&
RooUnfoldT<Hist,Hist2D>::Setup (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas)
{
  //! setup object from a response
  Reset();
  SetResponse (res);
  SetMeasured (meas);
  return *this;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const Hist* meas)
{
  
  //! Set measured distribution and errors. RooUnfold does not own the histogram.
  _meas= clone(meas);
  _cache = Cache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TVectorD& err)
{
  //! Set measured distribution and errors. Should be called after setting response matrix.
  const Hist* orig = _res->Hmeasured();
  _meas = RooUnfolding::createHist<Hist>(meas,GetName(),GetTitle(),var(orig,X));
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetTruth (const Hist* truth)
{
  
  //! Set truth distribution and errors. RooUnfold does not own the histogram.
  _truth= clone(truth);
  _cache = Cache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetTruth (const TVectorD& truth, const TVectorD& err)
{
  //! Set truth distribution and errors. Should be called after setting response matrix.
  const Hist* orig = _res->Htruth();
  _truth = RooUnfolding::createHist<Hist>(truth,GetName(),GetTitle(),var(orig,X));
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetBkg (const Hist* bkg)
{
  //! Set background distribution and errors. RooUnfold does not own the histogram.
  _bkg= clone(bkg);
  _cache = Cache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetBkg (const TVectorD& bkg, const TVectorD& err)
{
  //! Set truth distribution and errors. Should be called after setting response matrix.
  const Hist* orig = _res->Htruth();
  _bkg = RooUnfolding::createHist<Hist>(bkg,GetName(),GetTitle(),var(orig,X));
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasured (const TVectorD& meas, const TMatrixD& cov)
{
  //! Set measured distribution and its covariance matrix. Should be called after setting response matrix.
  SetMeasuredCov (cov);
  SetMeasured (meas, Emeasured());
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetMeasuredCov (const TMatrixD& cov)
{
  //! Set covariance matrix on measured distribution.
  _cache = Cache();
  _covMes= new TMatrixD (cov);
}

template<class Hist,class Hist2D> const TMatrixD&
RooUnfoldT<Hist,Hist2D>::GetMeasuredCov() const
{
  //! Get covariance matrix on measured distribution.
  if (_covMes) return *_covMes;
  auto err(Emeasured());
  _cache._covMes= new TMatrixD (_nm,_nm);
  for (Int_t i= 0 ; i<_nm; i++) {
    Double_t e= err[i];
    (*_cache._covMes)(i,i)= e*e;
  }
  return *_cache._covMes;
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::ForceRecalculation () const {
  //! clear and rebuild the cache
  this->_cache = Cache();
  this->_res->ClearCache();
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetResponse (const RooUnfoldResponseT<Hist,Hist2D>* res, Bool_t takeOwnership){
  //! Set response matrix for unfolding, optionally taking ownership of the RooUnfoldResponseT<Hist,Hist2D> object
  if(!res) throw std::runtime_error("cannot set response to invalid value!");
  if(takeOwnership) _res= const_cast<RooUnfoldResponseT<Hist,Hist2D>*>(res);
  else _res = new RooUnfoldResponseT<Hist,Hist2D>(*res);
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  
  SetNameTitleDefault();
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Unfold() const
{
  //! Dummy unfolding - just copies input
  cout << "********************** " << ClassName() << ": dummy unfolding - just copy input **********************" << endl;

  _cache._rec.ResizeTo (_nt);
  Int_t nb= _nm < _nt ? _nm : _nt;
  TVectorD vmeas(Vmeasured());
  for (Int_t i= 0; i < nb; i++) {
    _cache._rec(i)= vmeas(i);
  }
  _cache._unfolded= true;

}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrors() const
{
  //!Creates vector of diagonals of covariance matrices.
  if(this->_withError != kErrors){
    throw std::runtime_error("unknown error propagation method!");
  }

  if (!_cache._haveCov) GetCov();
  if (!_cache._haveCov) return;
  _cache._variances.ResizeTo(_nt);
  for (Int_t i= 0; i < _nt; i++) {
    _cache._variances(i)= _cache._cov(i,i);
  }
  _cache._haveCov = false;
  _cache._haveErrors= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetCov() const
{
  //!Dummy routine to get covariance matrix. It should be overridden by derived classes.
  const TMatrixD& covmeas(GetMeasuredCov());
  Int_t nb= std::min(_nm,_nt);
  _cache._cov.ResizeTo (_nt, _nt);
  for (int i=0; i<nb; i++){
    for (int j=0; j<nb; j++){
      _cache._cov(i,j)= covmeas(i,j);
    }
  }
  
  _cache._haveCov= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetWgt() const
{
  //! Creates weight matrix
  //! This may be overridden if it can be computed directly without the need for inverting the matrix 
  if (!_cache._haveCov) GetCov();
  if (!_cache._haveCov) return;
  if (!InvertMatrix (_cache._cov, _cache._wgt, "covariance matrix", _verbose)) return;
  _cache._haveWgt= true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetErrMat() const
{
  //! Get covariance matrix from the variation of the results in toy MC tests
  if (_NToys<=1) return;
  _cache._err_mat.ResizeTo(_nt,_nt);
  TVectorD xisum (_nt);
  TMatrixD xijsum(_nt,_nt);
  for (Int_t k=0; k<_NToys; k++){
    TVectorD x(_nt),xe(_nt);
    this->RunToy(x,xe);
    for (Int_t i=0; i<_nt;i++){
      Double_t xi= x[i];
      xisum[i] += xi;
      for (Int_t j=0; j<_nt; j++) xijsum(i,j) += xi * x[j];
    }
  }
  for (Int_t i=0; i<_nt; i++){
    for (Int_t j=0; j<_nt; j++){
      _cache._err_mat(i,j)= (xijsum(i,j) - (xisum[i]*xisum[j])/_NToys) / (_NToys-1);
    }
  }
  _cache._have_err_mat=true;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::CalculateBias(Int_t ntoys, const Hist* hTrue) const
{
  //! There are two toy approaches of calculating the bias.

  //! Use the response matrix truth if not supplied..
  TVectorD vtruth(hTrue ? h2v(hTrue,false) : _res->Vtruth());

  TVectorD vreco2(this->response()->Vfolded(vtruth));
  TVectorD vrecoerr(vreco2);

  //! Create un unfolding instance with the reconstructed histogram
  //! set as the measured histogram.
  Hist* asimov = RooUnfolding::asimov1DClone(this->response()->Hmeasured(),this->response()->UseDensityStatus(),vreco2,vrecoerr);
  auto* toyFactory = this->New(this->GetAlgorithm(),this->response(),asimov,GetRegParm());
  toyFactory->SetVerbose(0);

  //! Resize the bias vectors.
  _cache._bias.ResizeTo(_nt);
  _cache._sdbias.ResizeTo(_nt);
  _cache._sdmbias.ResizeTo(_nt);
  _cache._rmsbias.ResizeTo(_nt);    

  //! An array that will contain all the unfolded toys.
  std::vector<TVectorD> munfolded, etoys;
  std::vector<double> chi2;

  this->RunToys(this->_NToys, munfolded, etoys, chi2);    

  //! Calculate the bias and its stat. error with 
  //! the unfolded toys.
  for (int i = 0; i < vtruth.GetNrows(); i++){

    Double_t av_unfolded = 0;

    //! Sum over all toys.
    for (int j = 0; j < ntoys; j++){
      av_unfolded += munfolded[j][i];
    }
    
    //! Divide to get the average of the unfolded histograms.
    av_unfolded = av_unfolded / ntoys;

    //! RMS
    Double_t rms = 0;

    //! Variance
    Double_t var = 0;
    
    //! Calculate the sample variance of the unfolded histograms.
    for (int j = 0; j < ntoys; j++){
      rms += (munfolded[j][i] - vtruth(i))*(munfolded[j][i] - vtruth(i));
      var += (munfolded[j][i] - av_unfolded) * (munfolded[j][i] - av_unfolded);
    }
      
    //! standard error on the mean.
    _cache._sdmbias(i) = sqrt(var /(ntoys*(ntoys - 1)));

    //! standard error.
    _cache._sdbias(i)= sqrt(var / (ntoys - 1));

    //! Estimate the bias with the average of the unfolded histograms.
    _cache._bias(i) = av_unfolded - vtruth(i);

    //! Get the rms.
    _cache._rmsbias(i) = sqrt(rms/ntoys);
  }


  delete asimov;
  delete toyFactory;
  
  this->_cache._haveBias=true;
}

template<class Hist,class Hist2D> Bool_t
RooUnfoldT<Hist,Hist2D>::UnfoldWithErrors (ErrorTreatment withError, bool getWeights) const
{
  //! This method initializes the unfolding with errors.
  if (!_cache._unfolded) {

    if (_cache._fail) return false;
    
    this->Unfold();

    if (!_cache._unfolded) {
      _cache._fail= true;
      return false;
    }
  }
  
  Bool_t ok;
  if(_withError != withError) _cache._haveErrors = false;
  _withError= withError;
  if (getWeights && (withError==kErrors || withError==kCovariance)) {
    if   (!_cache._haveWgt)    GetWgt();
    ok= _cache._haveWgt;
  } else {
    switch (withError) {
    case kErrors:
      if (!_cache._haveErrors)     GetErrors();
      ok=_cache._haveErrors;
      break;
    case kCovariance:
      if   (!_cache._haveCov)      GetCov();
      ok= _cache._haveCov;
      break;
    case kErrorsToys:
      if   (!_cache._haveErrors)   GetErrorsToys();
      ok= _cache._haveErrors;
      break;
    case kCovToys:
      if   (!_cache._haveCov)   GetCovToys();
      ok= _cache._haveCov;
      break;
    default:
      ok= true;
    }
  }

  if (!ok) _cache._fail= true;

  return ok;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::Chi2(const Hist* hTrue,ErrorTreatment DoChi2) const {
    /*!Calculates Chi squared. Method depends on value of DoChi2
    0: sum of (residuals/error)squared
    1: use errors propagated through the unfolding
    2: use covariance matrix returned from unfolding
    3: use covariance matrix from the variation of the results in toy MC tests
    Returns warnings for small determinants of covariance matrices and if the condition is very large.
    If a matrix has to be inverted also removes rows/cols with all their elements equal to 0*/
    if (!UnfoldWithErrors (DoChi2)) return -1.0;
    TVectorD res(subtract<Hist,TVectorD>(_cache._rec,hTrue,_overflow));
    Double_t chi2= 0.0;
    if (DoChi2==kCovariance || DoChi2==kCovToys) {
      TMatrixD wgt(Wunfold(DoChi2));
      if (_cache._fail) return -1.0;
      TMatrixD resmat(1,_nt), chi2mat(1,1);
      TMatrixDRow(resmat,0)= res;
      ABAT (resmat, wgt, chi2mat);
      chi2= chi2mat(0,0);
    } else {
      TVectorD eunfold(EunfoldV(DoChi2));
      if (_cache._fail) return -1.0;
      for (Int_t i = 0 ; i < _nt; i++) {
        Double_t e= eunfold[i];
        if (e<=0.0) continue;
        Double_t ypull = res[i] / e;
        chi2 += ypull*ypull;
      }
    }
    return chi2;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::PrintTable (const Hist* hTrue, RooUnfolding::ErrorTreatment withError) const {
  //! Prints entries from truth, measured, and unfolded data for each bin.
  this->PrintTable(std::cout);
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::PrintTable (std::ostream& o, const Hist* hTrue, ErrorTreatment withError) const
{
  //! Prints entries from truth, measured, and unfolded data for each bin.
  if (withError==kDefault) withError= _withError;
  if (withError==kDefault) withError= kErrors;
  if (!UnfoldWithErrors (withError)) withError= kNoError;
  if (!hTrue){
    hTrue = response()->Htruth();
  }

  const Hist* hTrainTrue = response()->Htruth();
  const Hist* hTrain = response()->Hmeasured();
  const Hist* hMeas = Hmeasured();

  int ntxb= nBins(_res->Htruth(),X)+2*this->_overflow;
  int ntyb= nBins(_res->Htruth(),Y)+2*this->_overflow;

  int d = dim(_res->Htruth());
  if (!_cache._unfolded) return;
  Double_t chi_squ= -999.0;
  if (hTrue && (withError==kCovariance || withError==kCovToys )) chi_squ = Chi2(hTrue,withError);

  printTable(o,d,
             ntxb,ntyb,
             h2v(hTrainTrue,this->_overflow, this->response()->UseDensityStatus()),
             h2v(hTrain,this->_overflow, this->response()->UseDensityStatus()),
             hTrue ? h2v(hTrue,this->_overflow, this->response()->UseDensityStatus()) : TVectorD(this->_nt) ,             
             h2v(hMeas,this->_overflow, this->response()->UseDensityStatus()),
             this->Vunfold(),
             withError,
             hTrue ? h2ve(hTrue,this->_overflow, this->response()->UseDensityStatus()) : TVectorD(this->_nt) ,
             this->EunfoldV(withError),
             chi_squ);
}


template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::SetNameTitleDefault()
{
  if (!_res) return;
  const char* s= GetName();
  if (s[0] == '\0') SetName (_res->GetName());
  s= GetTitle();
  if (s[0] == '\0') {
    TString title= "Unfold ";
    title += _res->GetTitle();
    SetTitle (title);
  }
}

template<class Hist,class Hist2D> Hist*
RooUnfoldT<Hist,Hist2D>::Hunfold (ErrorTreatment withError)
{
    /*!Creates unfolded distribution. Error calculation varies by withError:
    0: No errors
    1: Errors from the square root of the diagonals of the covariance matrix given by the unfolding
    2: Errors from the square root of of the covariance matrix given by the unfolding
    3: Errors from the square root of the covariance matrix from the variation of the results in toy MC tests
    */
  if (!UnfoldWithErrors (withError)) withError= kNoError;
  if (!_cache._unfolded){
    return clone(_res->Htruth());
  } else {
    TVectorD rec(this->Vunfold());
    TVectorD errors(this->EunfoldV());
    Hist* unfolded = RooUnfolding::createHist<Hist>(rec, errors, GetName(), GetTitle(), _res->Htruth(), this->_overflow);    
    return unfolded;
  }
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetSettings() const
{
    //!Gets maximum and minimum parameters and step size
    _cache._minparm=0;
    _cache._maxparm=0;
    _cache._stepsizeparm=0;
    _cache._defaultparm=0;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMinParm() const
{
    //!Get minimum regularisation parameter for unfolding method
    return _cache._minparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetMaxParm() const
{
    //!Get maximum regularisation parameter for unfolding method
    return _cache._maxparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetStepSizeParm() const
{
    //!Get suggested step size for unfolding distribution
    return _cache._stepsizeparm;
}

template<class Hist,class Hist2D> Double_t
RooUnfoldT<Hist,Hist2D>::GetDefaultParm() const
{
    //!Get suggested regularisation parameter.
    return _cache._defaultparm;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Print(Option_t* /*opt*/) const
{
  //! print a summary of the configuration
  cout << ClassName() << "::" << GetName() << " \"" << GetTitle()
       << "\"," << " regularisation parameter=" << GetRegParm() << ", ";
  if (_covMes) cout << "with measurement covariance, ";
  if (_dosys)      cout << "calculate systematic errors, ";
  if (dim(_meas)==1) cout << _nm;
  else {
    cout <<        nBins(_meas,X)
         << "x" << nBins(_meas,Y);
    if (dim(_meas)>=3)
      cout << "x" << nBins(_meas,Z);
    cout << " (" << _nm << ")";
  }
  cout << " bins measured, ";
  const Hist* rtrue= _res->Htruth();
  if (dim(rtrue)==1) cout << _nt;
  else {
    cout <<        nBins(rtrue,X)
         << "x" << nBins(rtrue,Y);
    if (dim(rtrue)>=3)
      cout << "x" << nBins(rtrue,Z);
    cout << " (" << _nt << ")";
  }
  cout << " bins truth";
  if (_overflow) cout << " including overflows";
  cout << endl;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Dump() const {
  //! dump the contents of the member variables
  std::cout << "covMes=" <<  _covMes << std::endl;
  std::cout << "verbose=" <<  _verbose << std::endl;
  std::cout << "nm=" <<  _nm << std::endl;
  std::cout << "nt=" <<  _nt << std::endl;
  std::cout << "overflow=" <<  _overflow << std::endl;
  std::cout << "NToys=" <<  _NToys << std::endl;
  std::cout << "dosys=" <<  _dosys << std::endl;
  std::cout << "res=" <<  _res << std::endl;
  std::cout << "meas=" <<  _meas << std::endl;
  std::cout << "bkg=" << _bkg << std::endl;
  std::cout << "truth=" << _truth << std::endl;
  _res->Print();
  _meas->Print();
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::CutZeros(const TMatrixD& ereco)
{
    //!Removes row & column if all their elements are 0.
    vector<int> diags;
        int missed=0;
        for (int i=0; i<ereco.GetNrows(); i++){
            double coltot=0;
            for (int j=0;j<ereco.GetNrows();j++){
                coltot+=ereco(i,j);
            }
            if (coltot==0){
                diags.push_back(i);
                missed++;
            }
        }
        int x=ereco.GetNrows()-missed;
        int y=ereco.GetNcols()-missed;
        TMatrixD ereco_cut(x,y);
        unsigned int v=0;
        for (int i=0;i<ereco.GetNrows();i++){
            if(v<diags.size() && diags[v]==i){
                v++;
            }
            else{
                for (int j=0; j<ereco_cut.GetNcols();j++){
                    ereco_cut(i-v,j)=ereco(i,j+v);
                    }
                }
        }
    return ereco_cut;
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::Eunfold(ErrorTreatment withError) const
{
    /*!Returns covariance matrices for error calculation of type withError
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */
  
    if (!UnfoldWithErrors (withError)) return TMatrixD(_nt,_nt);

    switch(withError){
    case kNoError: {
      TMatrixD Eunfold_m(_nt,_nt);
      for (int i=0; i<_nt; i++){
        Eunfold_m(i,i)=_cache._rec(i);
      }
      return Eunfold_m;
      break; }
    case kErrors:
    case kErrorsToys: {
      TMatrixD Eunfold_m(_nt,_nt);
      for (int i=0; i<_nt;i++){
        Eunfold_m(i,i)=_cache._variances(i);
      }
      return Eunfold_m;
      break; 
    }
    case kCovariance:
      return _cache._cov;
      break;
    case kCovToys: {
      return _cache._cov;
      break;
    }
    default:
      throw std::runtime_error(TString::Format("Error in RooUnfoldT::Wunfold, unrecognised error method '%d'",withError).Data());                      
    }
}

template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::EunfoldV(ErrorTreatment withError) const
{
    /*!Returns vector of unfolding errors computed according to the withError flag:
    0: Errors are the square root of the bin content
    1: Errors from the diagonals of the covariance matrix given by the unfolding
    2: Errors from the covariance matrix given by the unfolding
    3: Errors from the covariance matrix from the variation of the results in toy MC tests
    */

    TVectorD Eunfold_v(_nt);
    if (!UnfoldWithErrors (withError)) return Eunfold_v;
    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._rec(i)));
        }
        break;
      case kErrors:
      case kErrorsToys:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._variances(i)));
        }
        break;
      case kCovariance:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._cov(i,i)));
        }
        break;
      case kCovToys:
        for (int i=0; i<_nt; i++){
          Eunfold_v(i)=sqrt (fabs (_cache._err_mat(i,i)));
        }
        break;
      default:
        throw std::runtime_error(TString::Format("Error in RooUnfoldT::EunfoldV, unrecognised error method '%d'",withError).Data());        
    }

    return Eunfold_v;
}

template<class Hist,class Hist2D> TMatrixD
RooUnfoldT<Hist,Hist2D>::Wunfold(ErrorTreatment withError) const
{
  //! TODO: document
    TMatrixD Wunfold_m(_nt,_nt);
    if (!UnfoldWithErrors (withError, true)) return Wunfold_m;

    switch(withError){
      case kNoError:
        for (int i=0; i<_nt; i++){
          if (_cache._rec(i)!=0.0) Wunfold_m(i,i)=1.0/_cache._rec(i);
        }
        break;
      case kErrors:
      case kErrorsToys:
        for (int i=0; i<_nt;i++){
          Wunfold_m(i,i)=_cache._wgt(i,i);
        }
        break;
      case kCovariance:
        Wunfold_m=_cache._wgt;
        break;
      case kCovToys:
        InvertMatrix (_cache._err_mat, Wunfold_m, "covariance matrix from toys", _verbose);
        break;
      default:
        throw std::runtime_error(TString::Format("Error in RooUnfoldT::Wunfold, unrecognised error method '%d'",withError).Data());                
    }


    return Wunfold_m;
}

// !The coverage calculation is based on Mikael Kuusela's PhD thesis p. 84 
// !paragraph 6.4.2. The assummptions for this closed form probability solution
// !are that the estimator is a linear function of the observed data and that
// !the observed bin counts follow a Gaussian distribution.

// !The input argument defines the confidence level with 1 sigma indicating
// !a confidence level of 0.6827, 2 sigma 0.9545 and 3 sigma 0.9973.

// !The output is a vector which contains the coverage probability for each bin.
template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::CoverageProbV(Int_t sigma) const
{
  
  TVectorD bias(_cache._bias);

  TVectorD coverage(_cache._bias.GetNrows());

  // Calculate the bias if needed.
  if (!this->_cache._haveBias){
    this->CalculateBias(this->_NToys);
    return coverage;
  }

  TVectorD se(this->EunfoldV(RooUnfolding::kErrorsToys));

  if (sigma < 1){
    std::cout << "Pass a positive integer to define the confidence interval" << std::endl;
    return coverage;
  }

  for (int i = 0; i < coverage.GetNrows(); i++){

    if (se(i)){
      coverage(i) = ROOT::Math::normal_cdf(bias(i)/se(i) + sigma) - ROOT::Math::normal_cdf(bias(i)/se(i) - sigma);
    } else {
      coverage(i) = 0;
    }
  }
  
  return coverage;
}


//! Scan the coverage probability for a given set of regularisation parameter values.
//! One can either do so for a specified bin or averaging over all bins(bin=-1).
//! One can also specify the confidence level(sigma).
template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::ScanCoverage(TVectorD& regparms, Int_t bin, Int_t sigma) const
{
  TVectorD coverageprobs(regparms.GetNrows());
  
  for (int i = 0; i < regparms.GetNrows(); i++){
    
    auto* toy_unfold = this->New(this->GetAlgorithm(),this->response(),this->Hmeasured(),regparms(i));
    
    TVectorD cov(toy_unfold->CoverageProbV(sigma));
    
    if (bin > 0 && bin < cov.GetNrows()){
      coverageprobs(i) = cov(bin);
    } else {
      coverageprobs(i) = cov.Sum() / cov.GetNrows();
    }

    delete toy_unfold;
  }
  
  return coverageprobs;
}

//! Scan the coverage probability for a given set of regularisation parameter values.
//! One can either do so for a specified bin or averaging over all bins(bin=-1).
//! One can also specify the confidence level(sigma).
template<class Hist,class Hist2D> TVectorD
RooUnfoldT<Hist,Hist2D>::ScanBias2Var(TVectorD& regparms, Int_t bin) const
{
  TVectorD bias2var(regparms.GetNrows());
  
  for (int i = 0; i < regparms.GetNrows(); i++){
    
    auto* toy_unfold = this->New(this->GetAlgorithm(),this->response(),this->Hmeasured(),regparms(i));
    
    //! Calculate the bias.
    toy_unfold->CalculateBias(this->_NToys);
  
    //! Get the unfolded distribution.
    TVectorD unfold(toy_unfold->Vunfold());

    //! Get the bias.
    TVectorD bias(toy_unfold->Vbias());

    //! Get the error on the unfolded result.
    TVectorD se(toy_unfold->EunfoldV(RooUnfolding::kErrorsToys));
    
    if (bin > 0 && bin < bias.GetNrows()){
      bias2var(i) = bias(bin)*bias(bin) + (se(bin)/unfold(bin))*(se(bin)/unfold(bin));
    } else {
      
      Double_t bias2varsum = 0;

      for (int j = 0; j < bias.GetNrows(); j++){
	bias2varsum += bias(j)*bias(j) + (se(j)/unfold(j))*(se(j)/unfold(j));
      }

      bias2var(i) = bias2varsum / bias.GetNrows();
    }

    delete toy_unfold;
  }
  
  return bias2var;
}


template<class Hist,class Hist2D> Int_t
RooUnfoldT<Hist,Hist2D>::InvertMatrix(const TMatrixD& mat, TMatrixD& inv, const char* name, Int_t verbose)
{
  //! Invert a matrix using Single Value Decomposition: inv = mat^-1.
  //! Can use InvertMatrix(mat,mat) to invert in-place.
  Int_t ok= 1;
  TDecompSVD svd (mat);
  const Double_t cond_max= 1e17;
  Double_t cond= svd.Condition();
  if (verbose >= 1) {
    Double_t d1=0,d2=0;
    svd.Det(d1,d2);
    Double_t det= d1*TMath::Power(2.,d2);
    cout << name << " condition="<<cond<<", determinant="<<det;
    if (d2!=0.0) cout <<" ("<<d1<<"*2^"<<d2<<")";
    cout <<", tolerance="<<svd.GetTol()<<endl;
  }
  if        (cond<0.0){
    cerr <<"Warning: bad "<<name<<" condition ("<<cond<<")"<<endl;
    ok= 2;
  } else if (cond>cond_max) {
    cerr << "Warning: poorly conditioned "<<name<<" - inverse may be inaccurate (condition="<<cond<<")"<<endl;
    ok= 3;
  }
  inv.ResizeTo (mat.GetNcols(), mat.GetNrows());  // pseudo-inverse of A(r,c) is B(c,r)
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,4)  /* TDecompSVD::Invert() didn't have ok status before 5.13/04. */
  Bool_t okinv= false;
  inv= svd.Invert(okinv);
  if (!okinv) {
    cerr << name << " inversion failed" << endl;
    return 0;
  }
#else
  inv= svd.Invert();
#endif
  if (verbose>=1) {
    TMatrixD I (mat, TMatrixD::kMult, inv);
    if (verbose>=3) printMatrix(I,"V*V^-1");
    Double_t m= 0.0;
    for (Int_t i= 0; i<I.GetNrows(); i++) {
      Double_t d= fabs(I(i,i)-1.0);
      if (d>m) m= d;
      for (Int_t j= 0; j<i; j++) {
        d= fabs(I(i,j)); if (d>m) m= d;
        d= fabs(I(j,i)); if (d>m) m= d;
      }
    }
    cout << "Inverse " << name << " " << 100.0*m << "% maximum error" << endl;
  }
  return ok;
}

template<class Hist,class Hist2D> void
RooUnfoldT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  //! Stream an object of class RooUnfold.
  if (R__b.IsReading()) {
    RooUnfoldT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
  } else {
    RooUnfoldT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template<> void
RooUnfoldT<TH1,TH2>::Streamer (TBuffer &R__b)
{
  //! Stream an object of class RooUnfold.
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldT<TH1,TH2>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldT<TH1,TH2>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist,class Hist2D>
RooUnfoldT<Hist,Hist2D>::RooUnfoldT()
  : TNamed()
{
  //! Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const char*    name, const char*    title)
  : TNamed(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::RooUnfoldT (const TString& name, const TString& title)
  : TNamed(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>::~RooUnfoldT()
{
  // destructor
}

template<class Hist,class Hist2D> 
RooUnfoldT<Hist,Hist2D>& RooUnfoldT<Hist,Hist2D>::operator= (const RooUnfoldT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfold settings.
  Assign(rhs);
  return *this;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::verbose() const
{
  //! Get verbosity setting which controls amount of information to be printed
  return _verbose;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::NToys()     const
{
  //! Get number of toys used in kCovToy error calculation.
  return _NToys;
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::Overflow()  const
{
  //! Histogram under/overflow bins are used?
  return _overflow;
}

template<class Hist,class Hist2D> 
const RooUnfoldResponseT<Hist,Hist2D>* RooUnfoldT<Hist,Hist2D>::response()  const
{
   //! Response matrix object
  return _res;
}

template<class Hist,class Hist2D> 
RooUnfoldResponseT<Hist,Hist2D>* RooUnfoldT<Hist,Hist2D>::response()
{  
   //! Response matrix object
  return _res;
}

template<class Hist,class Hist2D> 
const Hist*               RooUnfoldT<Hist,Hist2D>::Hmeasured() const
{
  //! Measured Distribution as a histogram
  return _meas;
}

template<class Hist,class Hist2D> 
Hist*               RooUnfoldT<Hist,Hist2D>::Hmeasured()
{
  //! Measured Distribution as a histogram
  return _meas;
}

template<class Hist,class Hist2D> 
const Hist*               RooUnfoldT<Hist,Hist2D>::Htruth() const
{
  //! Measured Distribution as a histogram
  return _truth;
}

template<class Hist,class Hist2D> 
Hist*               RooUnfoldT<Hist,Hist2D>::Htruth()
{
  //! Measured Distribution as a histogram
  return _truth;
}

template<class Hist,class Hist2D> 
const Hist*               RooUnfoldT<Hist,Hist2D>::Hbkg() const
{
  //! Measured Distribution as a histogram
  return _bkg;
}

template<class Hist,class Hist2D> 
Hist*               RooUnfoldT<Hist,Hist2D>::Hbkg()
{
  //! Measured Distribution as a histogram
  return _bkg;
}

template<class Hist,class Hist2D> 
const TVectorD&                RooUnfoldT<Hist,Hist2D>::Vunfold() const
{
  //! Unfolded (reconstructed) distribution as a vector
  if (!_cache._unfolded) {
    if (!_cache._fail){
      
      this->Unfold();
    }
    if (!_cache._unfolded) {
      _cache._fail= true;
      if (_nt > 0 && _cache._rec.GetNrows() == 0) _cache._rec.ResizeTo(_nt);   // need something
    }
  }

  return _cache._rec;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Vmeasured() const
{
  //! Measured distribution as a vector.
  if (!_cache._vMes){
    _cache._vMes = new TVectorD(h2v (_meas, _overflow, this->response()->UseDensityStatus()));
  }
  return *_cache._vMes;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Vtruth() const
{
  //! Measured distribution as a vector.
  if (!_cache._vTruth){
    _cache._vTruth = new TVectorD(h2v (_truth, _overflow, this->response()->UseDensityStatus()));
  }
  return *_cache._vTruth;
}

template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Vbkg() const
{
  //! Measured distribution as a vector.
  if (!_cache._vBkg){
    _cache._vBkg = new TVectorD(h2v (_bkg, _overflow, this->response()->UseDensityStatus()));
  }
  return *_cache._vBkg;
}

template<class Hist,class Hist2D> 
const TVectorD          RooUnfoldT<Hist,Hist2D>::Vbias() const
{
  //! Bias distribution as a vector.
  if (!_cache._haveBias){
    throw std::runtime_error("calculate bias before attempting to retrieve it!");
  }

  return _cache._bias;
}

template<class Hist,class Hist2D> 
const TVectorD          RooUnfoldT<Hist,Hist2D>::Ebias(RooUnfolding::BiasError E_type) const
{
  //! Bias errors as a vector.
  if (!_cache._haveBias){
    throw std::runtime_error("calculate bias before attempting to retrieve it!");
  }

  switch (E_type) {
  case kBiasSD:
    return _cache._sdbias;
  case kBiasSDM:
    return _cache._sdmbias;
  case kBiasRMS:
    return _cache._rmsbias;
  }
}


template<class Hist,class Hist2D> 
const TVectorD&          RooUnfoldT<Hist,Hist2D>::Emeasured() const
{
  //! Measured errors as a vector.
  if (!_cache._eMes){
    if(_covMes){
      _cache._eMes= new TVectorD(_nm);
      for (Int_t i= 0; i<_nm; i++) {
        Double_t e= this->_cache._cov(i,i);
        if (e>0.0) (*_cache._eMes)[i]= sqrt(e);
      }
    } else {
      _cache._eMes = new TVectorD(h2ve (_meas, _overflow, this->response()->UseDensityStatus()));
    }
  }
  return *_cache._eMes;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetVerbose (Int_t level)
{
  //! Set verbosity level which controls amount of information to be printed
  _verbose= level;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetOverflow (Int_t overflow)
{
  //! set the usage of the overflow bin
  _overflow= overflow;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetNToys (Int_t toys)
{
  //! Set number of toys used in kCovToy error calculation.
  _NToys= toys;
}

template<class Hist,class Hist2D> 
void  RooUnfoldT<Hist,Hist2D>::SetRegParm (Double_t regparm)
{
  //! Set Regularisation parameter
}

template<class Hist,class Hist2D> 
Double_t RooUnfoldT<Hist,Hist2D>::GetRegParm() const
{
  //! Get regularisation parameter.
  return -1e30;
}

template<class Hist,class Hist2D> 
void RooUnfoldT<Hist,Hist2D>::ClearCache() const
{
  //! Clear the cache
  this->_cache = Cache();
}

template<class Hist,class Hist2D> 
void RooUnfoldT<Hist,Hist2D>::IncludeSystematics (RooUnfolding::SystematicsTreatment dosys)
{
  //! Include systematic errors from response matrix?
  //! Use dosys=2 to exclude measurement errors.
  if (dosys!=_dosys){
    this->ClearCache();
    _dosys= dosys;
  }
}

template<class Hist,class Hist2D> 
Int_t RooUnfoldT<Hist,Hist2D>::SystematicsIncluded() const
{
  //! return setting for whether to include systematic errors from response matrix
  return _dosys;
}

template class RooUnfoldT<TH1,TH2>;
ClassImp (RooUnfold)

#ifndef NOROOFIT
#include "RooFitResult.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
template<> void RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::SetResponse (const RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* res, Bool_t takeOwnership){
  //! Set response matrix for unfolding, optionally taking ownership of the RooUnfoldResponseT<Hist,Hist2D> object
  if(!res) throw std::runtime_error("cannot set response to invalid value!");
  _res = new RooFitUnfoldResponse(res);
  _overflow= _res->UseOverflowStatus() ? 1 : 0;
  _nm= _res->GetNbinsMeasured();
  _nt= _res->GetNbinsTruth();
  SetNameTitleDefault();
}
namespace {
  void getParameters(const RooUnfolding::RooFitHist* hist, RooArgSet& params, std::string type){
    if(hist){
      RooArgSet* args = hist->func()->getParameters((RooArgSet*)0);
      for(auto p:*args){
        if(params.find(*p)) continue;
        RooRealVar* rrv = dynamic_cast<RooRealVar*>(p);
        if(!rrv) continue;
        if(rrv->isConstant()) continue;
        if(rrv->getError() == 0.){
          throw std::runtime_error(TString::Format("unable to build covariance matrix for parameter '%s' with error 0 - is this an observable? please set constant",rrv->GetName()).Data());
        }
	std::string paramName(rrv->GetName());
	if (paramName.find(type) != 0 && type != "all") continue;
        params.add(*rrv);
      }
      delete args;
    }
  }
  class FitResultHack : public RooFitResult {
  public:
    void setCovariance(TMatrixDSym& m){
      this->setCovarianceMatrix(m);
    }
  };
}



template<> void
RooUnfoldT<TH1,TH2>::RunRooFitToys(int ntoys, std::vector<TVectorD>& vx, std::vector<TVectorD>& vxe, std::vector<double>& chi2) const {
  this->RunToys(ntoys, vx, vxe, chi2);
}

template<> void
RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::RunRooFitToys(int ntoys, std::vector<TVectorD>& vx, std::vector<TVectorD>& vxe, std::vector<double>& chi2) const {
    
  //! Create un unfolding instance with the reconstructed histogram
  //! set as the measured histogram.
  RooUnfolding::RooFitHist* asimov = RooUnfolding::asimovClone(this->response()->Hmeasured(),this->response()->UseDensityStatus());
  auto* toyFactory = this->New(this->GetAlgorithm(),this->response(),asimov,GetRegParm());
  TVectorD vbkg_nom(this->response()->Vmeasured().GetNrows());

  //! Get the background if its there.
  if (this->Hbkg()){
    toyFactory->SetBkg(this->Hbkg());
    vbkg_nom = this->Vbkg();
  }

  toyFactory->SetVerbose(0);

  //! run a number of toys, fill the values, errors and chi2 in the
  //! given vectors
  const auto* res = toyFactory->response();
  RooArgSet errorParams;

  //! Get all other possible systematic and statistical uncertainties.
  std::string paramType;

  switch(this->_dosys) {
  case kNoSystematics:
    this->RunToys(ntoys, vx, vxe, chi2);
    return;
  case kGammas:
    paramType = "gamma";
    break;
  case kAlphas:
    paramType = "alpha";
    break;
  case kAll:
  case kNoStatistics:
    paramType = "all";
    break;
  }

  getParameters(res->Hmeasured(),errorParams,paramType);
  getParameters(res->Hresponse(),errorParams,paramType);
  getParameters(res->Htruth(),errorParams,paramType);
  getParameters(res->Hfakes(),errorParams,paramType);
  getParameters(toyFactory->Hbkg(),errorParams,paramType);

  //! Save the parameter values.
  auto* snsh = errorParams.snapshot();

  RooArgList errorParamList(errorParams);
  RooArgList errorParamListMu(errorParams);

  RooFitResult * prefitResult = RooFitResult::prefitResult(errorParams);

  //! Evaluate this part with Carsten.

  // if(_cache._covMes && !this->_dosys==kNoMeasured){
  //   auto meas(this->Vmeasured());
  //   auto covMes = *(_cache._covMes);
  //   auto setCov(prefitResult->covarianceMatrix());
  //   auto gammas = this->Hmeasured()->nps();
  //   for(size_t i=0; i<covMes.GetNcols(); ++i){
  //     RooRealVar* p1 = gammas[i];
  //     int idx1 = errorParamList.index(p1);
  //     if(idx1<0) continue;
  //     for(size_t j=0; j<covMes.GetNrows(); ++j){
  //       RooRealVar* p2 = gammas[j];
  //       int idx2 = errorParamList.index(p2);
  //       if(idx2<0) continue;
  //       double val = covMes(i,j)/(meas[i]*meas[j]);
  //       setCov(idx1,idx2) = val;
  //     }
  //   }
    
  //   ((::FitResultHack*)prefitResult)->setCovariance(setCov);
  // }

  RooRandom::randomGenerator()->SetSeed();

  //! Create a multidimensional pdf of which each nuisance parameter
  //! represents one dimension.
  const TMatrixDSym& V = prefitResult->covarianceMatrix();
  RooAbsPdf* paramPdf = new RooMultiVarGaussian("mvg", "mvg", errorParamList, errorParamListMu, V);

  //! Sample new values for these nuisance parameters.
  RooDataSet* d = paramPdf->generate(errorParamList,ntoys);

  Int_t failed_toys = 0;
  auto errorType = _withError;
  _withError = kDefault;
  for(int i=0; i<ntoys; ++i){
    errorParams = (*d->get(i));
    toyFactory->ForceRecalculation();

    //! Get the varied reco distribution.
    TVectorD vtruth = res->Vtruth();
    TVectorD vreco = res->Vfolded(vtruth);

    //! Get the varied background reco distribution.
    TVectorD vbkg(toyFactory->Vbkg());

    //! Add the background.
    TVectorD vreco_tot = vreco + vbkg;

    //! Sample from Poisson p.d.f.s the toy data. This
    //! will include statistical uncertainty.
    if (this->_dosys == kAll){
      RooUnfolding::randomize(vreco_tot, this->rnd);
    }

    //! Set the parameters i.e. distributions used in the unfolding
    //! back to their nominal values for the unfolding.
    errorParams = *snsh;
    toyFactory->ForceRecalculation();

    //! Set the sampled toy data distribution after subtracting the
    //! nominal background.
    toyFactory->_cache._vMes = new TVectorD(vreco_tot-vbkg_nom);

    //! add this extra check in case a toy unfolding failed
    if (toyFactory->Vunfold().GetNrows() == 1){
      failed_toys++;
      continue;
    }

    //! Save the unfolded result.
    vx.push_back(toyFactory->Vunfold());
    if(errorType != kNoError){

      //! Save the errors and chi2. Note thtat the errors are
      //! calculated here with the method specific estimation.
      vxe.push_back(toyFactory->EunfoldV(RooUnfolding::kErrors));
      chi2.push_back(toyFactory->Chi2 (toyFactory->response()->Htruth(), RooUnfolding::kErrors));
    }
  }

  //! set an maximum amount of retries to avoid the retry of toys
  //! loop 
  Int_t max_retries = ntoys;

  //! run an extra loop for failed toys.
  while(failed_toys != 0 && max_retries != 0){

    RooDataSet* d_retry = paramPdf->generate(errorParams,1);

    errorParams = (*d_retry->get(0)) ;
    toyFactory->ForceRecalculation();

    //! Get the varied reco distribution.
    TVectorD vtruth = res->Vtruth();
    TVectorD vreco = res->Vfolded(vtruth);

    //! Get the varied background reco distribution.
    TVectorD vbkg(toyFactory->Vbkg());

    //! Add the background.
    TVectorD vreco_tot = vreco + vbkg;

    //! Sample from Poisson p.d.f.s the toy data. This
    //! will include statistical uncertainty.
    if (this->_dosys == kAll){
      RooUnfolding::randomize(vreco_tot, this->rnd);
    }

    //! Set the parameters i.e. distributions used in the unfolding
    //! back to their nominal values for the unfolding.
    errorParams = *snsh;
    toyFactory->ForceRecalculation();

    //! Set the sampled toy data distribution after subtracting the
    //! nominal background.
    toyFactory->_cache._vMes = new TVectorD(vreco_tot-vbkg_nom);
  
    //! add this extra check in case a toy unfolding failed
    if (toyFactory->Vunfold().GetNrows() == 1){
      max_retries--;
      delete d_retry;
      continue;
    } 

    vx.push_back(toyFactory->Vunfold());
    if(errorType != kNoError){
      vxe.push_back(toyFactory->EunfoldV(RooUnfolding::kErrors));
      chi2.push_back(toyFactory->Chi2 (toyFactory->response()->Htruth(), RooUnfolding::kErrors));
    }

    failed_toys--;
    delete d_retry;
  }

  _withError =  errorType;
  
  errorParams = *snsh;
  delete snsh;
  delete prefitResult;
  delete paramPdf;
  delete d;
  delete asimov;
  delete toyFactory;
}

template<class Hist, class Hist2D> void
RooUnfoldT<Hist, Hist2D>::RunToys(int ntoys, std::vector<TVectorD>& vx, std::vector<TVectorD>& vxe, std::vector<double>& chi2) const
{
  
  //! Get the reconstructed histogram from the response matrix.
  TVectorD vreco(h2v(this->response()->Hmeasured(),this->_overflow, this->response()->UseDensityStatus()));
    
  //! Create un unfolding instance with the reconstructed histogram
  //! set as the measured histogram.
  Hist* asimov = RooUnfolding::asimovClone(this->response()->Hmeasured(),this->response()->UseDensityStatus());
  auto* toyFactory = this->New(this->GetAlgorithm(),this->response(),asimov,GetRegParm());
  TVectorD vbkg(this->response()->Vmeasured().GetNrows());

  //! Get the background if its there.
  if (this->Hbkg()){
    vbkg = this->Vbkg();
  }

  toyFactory->SetVerbose(0);

  //! Throw toys around the reconstructed histogram.
  for (int i = 0; i < ntoys; i++){
    
    TVectorD vreco_tot = vreco + vbkg;
    toyFactory->_cache._unfolded = false;
    
    //! Get a new histogram by sampling from Poisson distributions.
    RooUnfolding::randomize(vreco_tot, this->rnd);

    toyFactory->_cache._vMes = new TVectorD(vreco_tot - vbkg);

    //! Unfold.
    TVectorD vunfolded(toyFactory->Vunfold());
    
    //! Save the unfolded result for the sample variance.
    vx.push_back(vunfolded);
    if(_withError != kNoError){
      vxe.push_back(this->EunfoldV(RooUnfolding::kErrors));
      chi2.push_back(this->Chi2 (this->response()->Htruth(), RooUnfolding::kErrors));
    }
  }


  delete asimov;
  delete toyFactory;
}


template<class Hist, class Hist2D> double
RooUnfoldT<Hist, Hist2D>::RunToy(TVectorD&x, TVectorD&xe) const
{
  //! run a single toy, fill the values and errors in the given vectors
  //! returns the chi2
  std::vector<TVectorD> vx, vxe;
  std::vector<double> chi2;
  this->RunToys(1,vx,vxe,chi2);
  x.ResizeTo(vx[0].GetNrows());
  xe.ResizeTo(vxe[0].GetNrows());
  x = vx[0];
  xe = vxe[0];
  return chi2[0];
}

template<class Hist, class Hist2D> void
RooUnfoldT<Hist, Hist2D>::GetErrorsToys() const
{
  
  //! A vector with the unfolded results.
  std::vector<TVectorD> munfolded, etoys;
  std::vector<double> chi2;
  Int_t ntoys = this->_NToys;

  this->RunRooFitToys(ntoys, munfolded, etoys, chi2);

  GetSampleVar(munfolded);

  _cache._haveErrors = true;
}

template<class Hist, class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetCovToys() const
{

  //! A vector with the unfolded results.
  std::vector<TVectorD> munfolded, etoys;
  std::vector<double> chi2;
  Int_t ntoys = this->_NToys;
  
  this->RunRooFitToys(ntoys, munfolded, etoys, chi2);

  GetSampleCov(munfolded);

  _cache._haveCov = true;
}

template<class Hist, class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetSampleVar(std::vector<TVectorD>& munfolded) const
{
  
  Int_t ntoys = munfolded.size();
  
  _cache._variances.ResizeTo(_nt);
  
  //! Loop over the unfolded results.
  for (int i=0 ; i<this->_nt ; ++i) {
    
    double sum = 0;
   
    for (int j=0 ; j<ntoys ; ++j) {
      sum += munfolded[j][i];
    }
    double mu = sum/ntoys;
    double sum2 = 0;
    for (int j=0 ; j<ntoys ; ++j) {
      sum2 += (munfolded[j][i] - mu)*(munfolded[j][i] - mu);
    }
    _cache._variances(i) = sum2/(ntoys-1);
  }
}

template<class Hist, class Hist2D> void
RooUnfoldT<Hist,Hist2D>::GetSampleCov(std::vector<TVectorD>& munfolded) const
{
  //! Get covariance matrix from the variation of the results in toy MC tests
  _cache._cov.ResizeTo(_nt,_nt);

  Int_t ntoys = munfolded.size();

  TVectorD mu_av(_nt);

  //! Loop over the unfolded results.
  for (int i=0 ; i<this->_nt ; ++i) {
    
    double sum = 0;
   
    for (int j=0 ; j<ntoys ; ++j) {
      sum += munfolded[j][i];
    }

    mu_av(i) = sum/ntoys;
  }

  for (int i = 0; i <this->_nt; i++){
    for (int j = 0; j <this->_nt; j++){
      double sum = 0;
      for (int k = 0; k < ntoys; k++){
	sum += (munfolded[k][i] - mu_av(i))*(munfolded[k][j] - mu_av(j));
      }
      _cache._cov(i,j) = sum/(ntoys - 1);
    }
  }
}


template class RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
typedef RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooUnfoldT_RooFitHist;
ClassImp (RooUnfoldT_RooFitHist)
#endif
