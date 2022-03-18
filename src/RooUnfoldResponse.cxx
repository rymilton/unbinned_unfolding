//=====================================================================-*-C++-*-
//
// Description:
//      Response Matrix
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

//____________________________________________________________
/*! \class RooUnfoldResponseT
    \brief Class to create response object as used in RooUnfold 
<p> Contains measured and truth distributions and the response matrix. Also contains methods for handling these data</p>
<p> Can handle 1,2 or 3 dimensional histograms and return vectors and matrices of their bin content and error.
 Conversely can also convert these vectors and matrices into TH1s and TH2Ds. </p>
<p> Can also take a variety of parameters as inputs. This includes maximum and minimum values, distributions and vectors/matrices of values. </p>
<p> This class does the numerical modifications needed to allow unfolding techniques to work in the unfolding routines used in RooUnfold. </p>
\class RooUnfoldResponse
\brief specialization of RooUnfoldResponseT for TH1/TH2 objects
 */

/////////////////////////////////////////////////////////////

#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#include "RooDataHist.h"
#endif

#include <iostream>
#include <assert.h>
#include <cmath>

#include "TClass.h"
#include "TNamed.h"
#include "TBuffer.h"
#include "TPRegexp.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TCollection.h"

#include "TH1.h"
#include "TH2.h"


using namespace RooUnfolding;

#if ROOT_VERSION_CODE >= ROOT_VERSION(5,18,0)
#define HAVE_RooUnfoldFoldingFunction
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::pow;
using std::sqrt;


#ifdef HAVE_RooUnfoldFoldingFunction
template<class Hist, class Hist2D>
class RooUnfoldFoldingFunction {
public:
  RooUnfoldFoldingFunction<Hist,Hist2D> (const RooUnfoldResponseT<Hist,Hist2D>* res, TF1* func, Double_t eps=1e-12, bool verbose=false)
    : _res(res), _func(func), _eps(eps), _verbose(verbose), _fvals(_res->GetNbinsMeasured()) {
    _ndim= dynamic_cast<TF3*>(_func) ? 3 :
      dynamic_cast<TF2*>(_func) ? 2 : 1;
    if (_ndim>=2 && eps==1e-12) eps= 0.000001;
    FVals();
  }

  double operator() (double* x, double* p) const {
    const Hist* mes= _res->Hmeasured();
    Int_t bin;
    if      (_ndim==1) bin= findBin (mes, x[0]);
    else if (_ndim==2) bin= findBin (mes, x[0], x[1]);
    else               bin= findBin (mes, x[0], x[1], x[2]);
    if (bin<0 || bin>=_res->GetNbinsMeasured()) return 0.0;
    for (Int_t i=0, n=_func->GetNpar(); i<n; i++) {
      if (p[i] == _func->GetParameter(i)) continue;
      _func->SetParameters(p);
      FVals();
      break;
    }
    Double_t fy= _fvals[bin];
    if (_verbose) cout << "x=" << x[0] << ", bin=" << bin << " -> " << fy << endl;
    return fy;
  }

private:
  void FVals() const {
    const Hist* tru= _res->Htruth();
    if (_verbose) {
      cout << "p=";
      for (int i=0, n=_func->GetNpar(); i<n; i++) cout <<_func->GetParameter(i)<<",";
      cout << " f=";
    }
    _fvals.Zero();
    for (Int_t i=0, n=_res->GetNbinsTruth(); i<n; i++) {
      Int_t jx, jy, jz;
      binXYZ(tru,i,jx,jy,jz);
      Double_t fv;
      if (_eps<=0.0) {
        if (_ndim>=2)
          fv= _func->Eval (binCenter(tru,jx,RooUnfolding::X),
                           binCenter(tru,jy,RooUnfolding::Y),
                           binCenter(tru,jz,RooUnfolding::Z));
        else
          fv= _func->Eval (binCenter(tru,jx,RooUnfolding::X));
      } else {
        if        (_ndim==1) {
          Double_t tw= binWidth(tru,jx,RooUnfolding::X), tlo= binLowEdge(tru,jx,RooUnfolding::X), thi= tlo+tw;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
          fv= _func->Integral (tlo, thi, _eps) / tw;
        } else {
          Double_t tlo[3]= { binLowEdge(tru,jx,RooUnfolding::X), binLowEdge(tru,jy,RooUnfolding::Y), binLowEdge(tru,jz,RooUnfolding::Z) };
          Double_t thi[3]= { binHighEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jz,RooUnfolding::Z) };
          Double_t relerr=0.0;
          fv= _func->IntegralMultiple (_ndim, tlo, thi, _eps, relerr);
          fv /= binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y);
          if (_ndim>=3) fv /= binWidth(tru,jz,RooUnfolding::Z);
#else
          fv= _func->Integral (tlo, thi, (Double_t*)0, _eps) / tw;
        } else if (_ndim==2) {
          fv= _func->Integral (binLowEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jx,RooUnfolding::X),
                               binLowEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jy,RooUnfolding::Y),
                               _eps);
          fv /=binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y);
        } else {
          fv= _func->Integral (binLowEdge(tru,jx,RooUnfolding::X), binHighEdge(tru,jx,RooUnfolding::X),
                               binLowEdge(tru,jy,RooUnfolding::Y), binHighEdge(tru,jy,RooUnfolding::Y),
                               binLowEdge(tru,jz,RooUnfolding::Z), binHighEdge(tru,jz,RooUnfolding::Z),
                               _eps);
          fv /=binWidth(tru,jx,RooUnfolding::X) * binWidth(tru,jy,RooUnfolding::Y) * binWidth(tru,jz,RooUnfolding::Z);          
#endif
        }
      }
      if (_verbose) cout << " " << fv;
      for (Int_t ibin=0, m=_res->GetNbinsMeasured(); ibin<m; ibin++) {
        _fvals[ibin] += fv * (*_res)(ibin,i);
      }
    }
    if (_verbose) cout << endl;
  }

  const RooUnfoldResponseT<Hist,Hist2D>* _res;
  TF1* _func;
  Double_t _eps;
  bool _verbose;
  mutable TVectorD _fvals;
  Int_t _ndim;
};
#endif  

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::setup()
{
  _tru= _mes= _fak= _tru_m = 0;
  _res= 0;
  this->ClearCache();
  SetNameTitleDefault ("response", "Response");
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::ClearCache() const
{
  _cache = Cache();
  _cached= false;
}

template <class Hist, class Hist2D> Hist2D*
RooUnfoldResponseT<Hist,Hist2D>::HresponseNoOverflow() const
{
  const Hist2D* res = Hresponse();
  TVectorD vals(h2v<Hist>(res,_overflow,_density));
  TVectorD errs(h2ve<Hist>(res,_overflow,_density));  
  return createHist<Hist2D>(vals,errs,name(res),title(res),vars(res),_overflow);
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Print (Option_t* /* option */) const
{
  std::cout << "RooUnfoldResponseT @ " << this << std::endl;
  printHistogram(this->Hmeasured());
  printHistogram(this->Hfakes());
  printHistogram(this->Htruth());
  printHistogram(this->Hresponse());
}

template <class Hist, class Hist2D> bool
RooUnfoldResponseT<Hist,Hist2D>::Cached () const
{
  return _cached;
}

template<class Hist,class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::Cache::Cache() :
  _vMes(0),   
  _eMes(0),   
  _vFak(0),   
  _vTru(0),   
  _eTru(0),   
  _mRes(0),   
  _eRes(0),   
  _mResNorm(0),
  _eResNorm(0),
  _vPur(0),
  _vEff(0)
{
  //! default constructor
}

template<class Hist,class Hist2D> 
typename RooUnfoldResponseT<Hist,Hist2D>::Cache& RooUnfoldResponseT<Hist,Hist2D>::Cache::operator= ( const RooUnfoldResponseT<Hist,Hist2D>::Cache & other ){
  //! assignment operator
  if(this->_vMes)     { delete this->_vMes;     }  this->_vMes     = other._vMes;    
  if(this->_eMes)     { delete this->_eMes;     }  this->_eMes     = other._eMes;    
  if(this->_vFak)     { delete this->_vFak;     }  this->_vFak     = other._vFak;    
  if(this->_vTru)     { delete this->_vTru;     }  this->_vTru     = other._vTru;    
  if(this->_eTru)     { delete this->_eTru;     }  this->_eTru     = other._eTru;    
  if(this->_mRes)     { delete this->_mRes;     }  this->_mRes     = other._mRes;    
  if(this->_eRes)     { delete this->_eRes;     }  this->_eRes     = other._eRes;    
  if(this->_mResNorm) { delete this->_mResNorm; }  this->_mResNorm = other._mResNorm;
  if(this->_eResNorm) { delete this->_eResNorm; }  this->_eResNorm = other._eResNorm;
  if(this->_vEff)     { delete this->_vEff;     }  this->_vEff     = other._vEff;    
  if(this->_vPur)     { delete this->_vPur;     }  this->_vPur     = other._vPur;      
  return *this;
}


template<class Hist,class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::Cache::~Cache(){
  //! destructor
  if(this->_vMes)     delete this->_vMes;    
  if(this->_eMes)     delete this->_eMes;    
  if(this->_vFak)     delete this->_vFak;    
  if(this->_vTru)     delete this->_vTru;    
  if(this->_eTru)     delete this->_eTru;    
  if(this->_mRes)     delete this->_mRes;    
  if(this->_eRes)     delete this->_eRes;    
  if(this->_mResNorm) delete this->_mResNorm;
  if(this->_eResNorm) delete this->_eResNorm;
  if(this->_vEff) delete this->_vEff;
  if(this->_vPur) delete this->_vPur;  
}

template <class Hist, class Hist2D> Hist*
RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth (const Hist* truth, const char* name) const
{
  //! Apply the response matrix to the truth
  //! Errors not set, since we assume original truth has no errors
  if (!Htruth()) return 0;  // Needed for checking binning if nothing else

  // If no truth histogram input, use training truth
  // If truth histogram input, make sure its binning is correct
  TVectorD resultvect;;
  if (truth) {
    if (nBins(truth,RooUnfolding::X) != nBins(_tru,RooUnfolding::X) ||
        nBins(truth,RooUnfolding::Y) != nBins(_tru,RooUnfolding::Y) ||
        nBins(truth,RooUnfolding::Z) != nBins(_tru,RooUnfolding::Z))
      cerr << "Warning: RooUnfoldResponseT<Hist,Hist2D>::ApplyToTruth truth histogram is a different size ("
           << (nBins(truth,RooUnfolding::X) * nBins(truth,RooUnfolding::Y) * nBins(truth,RooUnfolding::Z)) << " bins) or shape from response matrix truth ("
           << ( nBins(_tru,RooUnfolding::X) * nBins( _tru,RooUnfolding::Y) * nBins( _tru,RooUnfolding::Z)) << " bins)" << endl;
    resultvect= h2v (truth, _overflow,_density);
  } else {
    resultvect= Vtruth();
  }

  resultvect *= Mresponse();   // v= A*v

  // Turn results vector into properly binned histogram
  const Hist* t = Hmeasured();
  Hist* result= createHist<Hist>(resultvect, RooUnfolding::name(t),name, vars(t), _overflow);
  return result;
}


template <class Hist, class Hist2D> TF1*
RooUnfoldResponseT<Hist,Hist2D>::MakeFoldingFunction (TF1* func, Double_t eps, Bool_t verbose) const
{
  //! Creates a function object that applies the response matrix to a user parametric function.
  //! This can be fitted to the measured distribution as an alternative to unfolding.
  //! The returned object is owned by the caller. The function will be binned.
  //! Specify eps=0 to calculate function at bin centers; otherwise integrates over each bin (may be slow).
  //! Example:
  //!    TF1* func= new TF1 ("func", "gaus", 0, 10);
  //!    TF1* fold= respose->MakeFoldingFunction(func);
  //!    histMeasured->Fit(fold);
  //!    fold->Draw("h"); // draw function fitted to histMeasured
  //!    func->Draw();    // draw truth function
#ifdef HAVE_RooUnfoldFoldingFunction
  Int_t np= func->GetNpar();
  RooUnfoldFoldingFunction<Hist,Hist2D> ff (this, func, eps, verbose);
  TString name= func->GetName();
  name += "_folded";
  TF1* f;
  if        (TF3* func3= dynamic_cast<TF3*>(func))
    f= new TF3 (name, ROOT::Math::ParamFunctor(ff),
                func3->GetXmin(), func3->GetXmax(),
                func3->GetYmin(), func3->GetYmax(),
                func3->GetZmin(), func3->GetZmax(), np);
  else if (TF2* func2= dynamic_cast<TF2*>(func))
    f= new TF2 (name, ROOT::Math::ParamFunctor(ff),
                func2->GetXmin(), func2->GetXmax(),
                func2->GetYmin(), func2->GetYmax(), np);
  else
    f= new TF1 (name, ROOT::Math::ParamFunctor(ff),
                func ->GetXmin(), func ->GetXmax(), np);
  f->SetNpx (GetNbinsMeasured()<=2 ? 4 : GetNbinsMeasured()==3 ? 6 : GetNbinsMeasured());  // TF1 requires Npx>=4
  // Copy parameters in case we set them in func
  f->SetParameters (func->GetParameters());
  f->SetParErrors  (func->GetParErrors());
  for (Int_t i=0; i<np; i++) {
    Double_t plo=0.0, phi=0.0;
    func->GetParLimits (i, plo, phi);
    f   ->SetParLimits (i, plo, phi);
    f->SetParName (i, func->GetParName(i));
  }
  return f;
#else
  cerr << "RooUnfoldResponseT<Hist,Hist2D>::MakeFoldingFunction not supported in this version of ROOT" << endl;
  return 0;
#endif
}
template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::RunToy() const {
  this->ClearCache();
  Vmeasured()     ;
  Vtruth()        ;
  Mresponse(false) ;
  RooUnfolding::randomize(*_cache._vMes,*_cache._eMes,this->rnd);
  RooUnfolding::randomize(*_cache._vTru,*_cache._eTru,this->rnd);
  RooUnfolding::randomize(*_cache._mRes,*_cache._eRes,this->rnd);
  _cache._mResNorm = new TMatrixD(*_cache._mRes);
  mNorm(*_cache._mResNorm,*_cache._vTru);
}

template <class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::SetNameTitleDefault (const char* defname, const char* deftitle)
{
  //! Set object name and title
  const char* s= GetName();
  if (s[0] == '\0') {
    if (_res) s= name(_res);
    if (s[0] == '\0') {
      if (defname) SetName (defname);
      else if (_mes && _tru) {
        TString n= name(_mes);
        if (n.Length()) n.Append ("_");
        n.Append (name(_tru));
        if (!n.Length()) n= "response";
        SetName (n);
      }
    } else
      SetName (s);
  }
  s= GetTitle();
  if (s[0] == '\0') {
    if (_res) s= title(_res);
    if (s[0] == '\0') {
      if (deftitle) SetTitle (deftitle);
      else if (_mes && _tru) {
        TString n= title(_tru);
        if (n.Length()) n.Append (" #rightarrow ");
        n.Append (title(_mes));
        if (n.Length())
          n.Prepend ("Response ");
        else
          n= "Response";
        SetTitle (n);
      }
    } else
      SetTitle (s);
  }
}

template<class Hist, class Hist2D> void
RooUnfoldResponseT<Hist,Hist2D>::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    RooUnfoldResponseT<Hist,Hist2D>::Class()->ReadBuffer  (R__b, this);
  } else {
    RooUnfoldResponseT<Hist,Hist2D>::Class()->WriteBuffer (R__b, this);
  }
}

template <> void
RooUnfoldResponseT<TH1,TH2>::Streamer (TBuffer &R__b)
{
  if (R__b.IsReading()) {
    // Don't add our histograms to the currect directory.
    // We own them and we don't want them to disappear when the file is closed.
    Bool_t oldstat= TH1::AddDirectoryStatus();
    TH1::AddDirectory (kFALSE);
    RooUnfoldResponseT<TH1,TH2>::Class()->ReadBuffer  (R__b, this);
    TH1::AddDirectory (oldstat);
  } else {
    RooUnfoldResponseT<TH1,TH2>::Class()->WriteBuffer (R__b, this);
  }
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT()
  : TNamed()
{
  //!default constructor. Use setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const char*    name, const char*    title)
  : TNamed(name,title)
{
  //! default named constructor. Use setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const TString& name, const TString& title)
  : TNamed(name,title)
{
  //! default named constructor. Use setup() to set values.
  this->setup();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const RooUnfoldResponseT<Hist,Hist2D>& rhs) : 
  TNamed(rhs.GetName(),rhs.GetTitle()),
  _mes(clone(rhs._mes)),
  _mestru(clone(rhs._mestru)),
  _fak(clone(rhs._fak)),      
  _tru(clone(rhs._tru)),
  _tru_m(clone(rhs._tru_m)),
  _res(clone(rhs._res)),
  _overflow(rhs._overflow),
  _density(rhs._density)
{
  // copy constructor
  this->ClearCache();
}

template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::RooUnfoldResponseT(const char* name, const char* title, Hist2D* response, Hist* truth, Hist* reco, bool overflow,bool density) :
  TNamed(name,title),
  _mes(reco),
  _tru(truth),
  _res(response),
  _overflow(overflow),
  _density(density)
{
  // explicit constructor
  this->ClearCache();
}



template<class Hist, class Hist2D>
RooUnfoldResponseT<Hist,Hist2D>::~RooUnfoldResponseT()
{
  //! destructor
  ClearCache();
}


template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetDimensionMeasured() const
{
  //! Dimensionality of the measured distribution (1=1D, 2=2D, 3=3D)
  return ::dim(this->_mes);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetDimensionTruth() const
{
  //! Dimensionality of the truth distribution (1=1D, 2=2D, 3=3D)
  return ::dim(this->_tru);  
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetNbinsMeasured() const
{
  //! Total number of bins in the measured distribution
  return ::nBins(this->_mes,_overflow);
}

template<class Hist, class Hist2D>
Int_t RooUnfoldResponseT<Hist,Hist2D>::GetNbinsTruth() const
{
  //! Total number of bins in the truth distribution
  return ::nBins(this->_tru,_overflow);
}


template<class Hist, class Hist2D>
const Hist* RooUnfoldResponseT<Hist,Hist2D>::Hmeasured() const
{
  //! Measured distribution, including fakes
  return _mes;
}


template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Hmeasured()
{
  //! Measured distribution, including fakes
  return _mes;
}


template<class Hist, class Hist2D>
const Hist* RooUnfoldResponseT<Hist,Hist2D>::Hfakes() const
{
  //! Fakes distribution
  return _fak;
}


template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Hfakes()
{
  //! Fakes distribution
  return _fak;
}

template<class Hist, class Hist2D>
const Hist*   RooUnfoldResponseT<Hist,Hist2D>::Htruth() const
{
  //! Truth distribution, used for normalisation
  return _tru;
}

template<class Hist, class Hist2D>
Hist*         RooUnfoldResponseT<Hist,Hist2D>::Htruth()
{
  //! Truth distribution, used for normalisation
  return _tru;
}

template<class Hist, class Hist2D>
const Hist2D*   RooUnfoldResponseT<Hist,Hist2D>::Hresponse() const
{
  //! Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  return _res;
}

template<class Hist, class Hist2D>
Hist2D*         RooUnfoldResponseT<Hist,Hist2D>::Hresponse()
{
  //! Response matrix as a 2D-histogram: (x,y)=(measured,truth)
  return _res;
}

template<class Hist, class Hist2D>
TVectorD        RooUnfoldResponseT<Hist,Hist2D>::Vfolded(const TVectorD& truth) const {
  auto res = this->Mresponse(true);
  if(!truth.GetNrows() == res.GetNcols())  throw std::runtime_error("Error in RooUnfoldResponseT::Vfolded: invalid dimensionality in given truth vector!");
  auto prod = res*truth;
  return prod;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vmeasured() const
{
  //! Measured distribution as a TVectorD
  if (!_cache._vMes) _cached= (_cache._vMes= new TVectorD(h2v  (_mes, _overflow,_density)));
  if(_cache._vMes->GetNrows() != this->GetNbinsMeasured()) throw std::runtime_error("Error in RooUndfoldResponseT::Vmeasured: invalid dimensionality in measured vector!");
  return *_cache._vMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vefficiency() const
{
  if(!this->_cache._vEff){
    TMatrixD resp(h2m(_res, false, _density));
    TVectorD truth(h2v(_tru, false, _density));
    
    this->_cache._vEff = new TVectorD(resp.GetNcols());
    
    for (int i = 0; i < resp.GetNcols(); i++){
      
      double n_truth = truth[i];
      double n_truth_reco = 0;
      
      for (int i_r = 0; i_r < resp.GetNrows(); i_r++){
        
        n_truth_reco += resp[i_r][i];
      }
      
      (*this->_cache._vEff)[i] = (Double_t)n_truth_reco/n_truth;
    }
  }

  return *(this->_cache._vEff);
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vpurity() const
{
  if(!this->_cache._vPur){
    const TVectorD& reco(Vmeasured());
    const TMatrixD& resp(Mresponse(false));
    
    this->_cache._vPur = new TVectorD(resp.GetNcols());
    
    for (int i = 0; i < resp.GetNcols(); i++){
    Int_t n_recogen = resp[i][i];
    Int_t n_reco = reco[i];
    
    (*_cache._vPur)[i] = (Double_t)n_recogen/n_reco;
    }
  }
  
  return *(this->_cache._vPur);
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vfakes() const
{
  //! Fakes distribution as a TVectorD
  if(!_fak){
    _cached = (_cache._vFak = new TVectorD(this->GetNbinsMeasured()) );
  } else {
    if (!_cache._vFak) _cached= (_cache._vFak= new TVectorD(h2v  (_fak, _overflow,_density)));
  }
  if(_cache._vFak->GetNrows() != this->GetNbinsMeasured()) throw std::runtime_error("invalid dimensionality in fakes vector!");
  return *_cache._vFak;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Emeasured() const
{
  //! Measured distribution errors as a TVectorD
  if (!_cache._eMes) _cached= (_cache._eMes= new TVectorD(h2ve (_mes, _overflow,_density)));
  if(_cache._eMes->GetNrows() != this->GetNbinsMeasured()) throw std::runtime_error("invalid dimensionality in measured uncertainty vector!");
  return *_cache._eMes;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Vtruth() const
{
  //! Truth distribution as a TVectorD
  if (!_cache._vTru) _cached= (_cache._vTru= new TVectorD(h2v  (_tru, _overflow,_density))); 
  if(_cache._vTru->GetNrows() != this->GetNbinsTruth()) throw std::runtime_error("invalid dimensionality in truth vector!");
  return *_cache._vTru;
}

template<class Hist, class Hist2D>
const TVectorD& RooUnfoldResponseT<Hist,Hist2D>::Etruth() const
{
  //! Truth distribution errors as a TVectorD
  if (!_cache._eTru) _cached= (_cache._eTru= new TVectorD(h2ve (_tru, _overflow,_density))); 
  if(_cache._eTru->GetNrows() != this->GetNbinsTruth()) throw std::runtime_error("invalid dimensionality in truth uncertainty vector!");
  return *_cache._eTru;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Mresponse(bool norm) const
{
  //! Response matrix as a TMatrixD: (row,column)=(measured,truth)

  if(norm){
    if (!_cache._mResNorm) _cached= (_cache._mResNorm= new TMatrixD(h2mNorm  (_res, _tru, _overflow,_density)));
    return *_cache._mResNorm;
  } else {
    if (!_cache._mRes) _cached = (_cache._mRes = new TMatrixD(h2m (_res, _overflow, _density)));
    return *_cache._mRes;
  }
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Mresponse(const Hist* truth) const
{
  //! Response matrix as a TMatrixD: (row,column)=(measured,truth)

    _cached= (_cache._mResNorm= new TMatrixD(h2mNorm  (_res, truth, _overflow,_density)));
    return *_cache._mResNorm;
}

template<class Hist, class Hist2D>
const TMatrixD& RooUnfoldResponseT<Hist,Hist2D>::Eresponse(bool norm) const
{
  //! Response matrix errors as a TMatrixD: (row,column)=(measured,truth)
  if (norm){
    if(!_cache._eResNorm) _cached= (_cache._eResNorm= new TMatrixD(h2meNorm (_res, _tru, _overflow,_density)));
    return *_cache._eResNorm;
  } else {
    if(!_cache._eRes) _cached= (_cache._eRes= new TMatrixD(h2me (_res, _overflow,_density)));
    return *_cache._eRes;
  }
}


template<class Hist, class Hist2D>
Double_t RooUnfoldResponseT<Hist,Hist2D>::operator() (Int_t r, Int_t t) const
{
  //! Response matrix element (measured,truth)
  return Mresponse()(r,t);
}

template<class Hist, class Hist2D>
void RooUnfoldResponseT<Hist,Hist2D>::UseOverflow (Bool_t set)
{
  //! Specify to use overflow bins. Only supported for 1D truth and measured distributions.
  _overflow= (set ? 1 : 0);
}

template<class Hist, class Hist2D>
Bool_t RooUnfoldResponseT<Hist,Hist2D>::UseOverflowStatus() const
{
  //! Get UseOverflow setting
  return _overflow;
}

template<class Hist, class Hist2D>
Bool_t RooUnfoldResponseT<Hist,Hist2D>::UseDensityStatus() const
{
  //! Get UseDensity setting
  return _density;
}

template<class Hist, class Hist2D>
bool RooUnfoldResponseT<Hist,Hist2D>::HasFakes() const
{
  //! Return number of fake entries
  return _fak && !empty(_fak);
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
template class RooUnfoldResponseT<TH1,TH2>;

namespace {
  template<class Hist> int entries(const Hist* hist);

  template<class Hist> int entries(const Hist* hist){
    return hist->GetEntries();
  }

  template int entries<TH1>(TH1 const*);
  template int entries<TH2>(TH2 const*);

  template<class Hist, class Hist2D> void projectY(Hist2D* _res, Hist* _tru, bool overflow);
  template<class Hist, class Hist2D> void projectX(Hist2D* _res, Hist* _mes, bool overflow);
  template<class Hist, class Hist2D> void subtractProjectX(Hist2D* _res, Hist* _mes, Hist* _fak, bool overflow);

  template<> void projectY<TH1>(TH2* _res, TH1* _tru, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t j= 1-overflow; j<_res->GetNbinsY()+1+overflow; j++) {
      Double_t ntru= 0.0, wtru= 0.0;
      for (Int_t i= 0; i<_res->GetNbinsX()+2; i++) {
        ntru +=      _res->GetBinContent (i, j);
        if (s) wtru += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= bin<TH1>(_tru, j, overflow);
      _tru->SetBinContent (b,      ntru);
      if (s) _tru->SetBinError   (b, sqrt(wtru));
    }
  }
  template<> void projectX<TH1>(TH2* _res, TH1* _mes, bool overflow){
    Int_t s= _res->GetSumw2N();
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY(); j++) {
        nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= RooUnfolding::bin (_mes, i, overflow);
      _mes->SetBinContent (b,      nmes );
      if (s) _mes->SetBinError   (b, sqrt(wmes));
    }
  }
  template<> void subtractProjectX<TH1>(TH2* _res, TH1* _mes, TH1* _fak, bool overflow){
    Int_t s= _res->GetSumw2N();
    Int_t sm= _mes->GetSumw2N(), nfake=0;
    for (Int_t i= 1-overflow; i<_res->GetNbinsX()+1+overflow; i++) {
      Double_t nmes= 0.0, wmes= 0.0;
      for (Int_t j= 0; j<_res->GetNbinsY()+2; j++) {
        nmes +=      _res->GetBinContent (i, j);
        if (s) wmes += pow (_res->GetBinError   (i, j), 2);
      }
      Int_t b= RooUnfolding::bin (_mes, i, overflow);
      Double_t fake= _mes->GetBinContent (b) - nmes;
      if (fake!=0.0) nfake++;
      if (!s) wmes= nmes;
      _fak->SetBinContent (b, fake);
      _fak->SetBinError   (b, sqrt (wmes + (sm ? pow(_mes->GetBinError(b),2) : _mes->GetBinContent(b))));
    }
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,0)
    _fak->SetEntries (_fak->GetEffectiveEntries());  // 0 entries if 0 fakes
#else
    _fak->SetEntries (nfake);  // 0 entries if 0 fakes
#endif    
  }


  template<class Hist> int fill(Hist* hist, double x, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double w);
  template<class Hist> int fill(Hist* hist, double x, double y, double z, double w);
  template<> int fill<TH1>(TH1* hist, double x, double w){
    return hist->Fill (x, w);
  }
  template<> int fill<TH1>(TH1* hist, double x, double y, double w){
    return ((TH2*)hist)->Fill (x, y, w);
  }
  template<> int fill<TH1>(TH1* hist, double x, double y, double z, double w){
    return ((TH3*)hist)->Fill (x, y, z, w);
  }    
  template<> int fill<TH2>(TH2* hist, double x, double y, double w){
    return hist->Fill (x, y, w);
  } 


}




Int_t RooUnfoldResponse::Miss (Double_t xt)
{
  //! Fill missed event into 1D Response Matrix
  return Miss1D(xt);
}

Int_t RooUnfoldResponse::Miss (Double_t xt, Double_t w)
{
  //! Fill missed event into 1D (with weight) or 2D Response Matrix
  return GetDimensionTruth()==2 ? Miss2D(xt,w) : Miss1D(xt,w);
}

Int_t RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t w)
{
  //! Fill missed event into 2D (with weight) or 3D Response Matrix
  return GetDimensionTruth()==3 ? Miss(xt,yt,w,1.0) : Miss2D(xt,yt,w);
}


Int_t RooUnfoldResponse::Fake (Double_t xr)
{
  //! Fill fake event into 1D Response Matrix
  return Fake1D(xr);
}

Int_t RooUnfoldResponse::Fake (Double_t xr, Double_t w)
{
  //! Fill fake event into 1D (with weight) or 2D Response Matrix
  return GetDimensionMeasured()==2 ? Fake2D(xr,w) : Fake1D(xr,w);
}

Int_t RooUnfoldResponse::Fake (Double_t xr, Double_t yr, Double_t w)
{
  //! Fill fake event into 2D (with weight) or 3D Response Matrix
  return GetDimensionMeasured()==3 ? Fake(xr,yr,w,1.0) : Fake2D(xr,yr,w);
}




RooUnfoldResponse::RooUnfoldResponse(const RooUnfoldResponse& rhs)
  : RooUnfoldResponseT(rhs.GetName(), rhs.GetTitle())
{
  //! RooUnfoldResponseT<class Hist, class Hist2D> copy constructor
  Setup(rhs);
}

RooUnfoldResponse&
RooUnfoldResponse::operator=(const RooUnfoldResponse& rhs)
{
  //! RooUnfoldResponseT<class Hist, class Hist2D> assignment operator
  if (this == &rhs) return *this;
  Reset();
  SetNameTitle(rhs.GetName(), rhs.GetTitle());
  return Setup(rhs);
}



RooUnfoldResponse&
RooUnfoldResponse::Setup(const RooUnfoldResponse& rhs)
{
  //! Copy data from another RooUnfoldResponseT<class Hist, class Hist2D>
  _overflow= rhs._overflow;
  return Setup(clone(rhs.Hmeasured()), clone(rhs.Htruth()), clone(rhs.Hresponse()));
}


RooUnfoldResponse&
RooUnfoldResponse::Reset()
{
  //! Resets object to initial state.
  if(_mes) delete  _mes;
  if(_tru_m) delete _tru_m;
  if(_mestru) delete _mestru;
  if(_tru) delete  _tru;
  if(_res) delete  _res;
  if(_fak) delete _fak;
  setup();
  return *this;
}


RooUnfoldResponse::RooUnfoldResponse(Int_t nb, Double_t xlo, Double_t xhi,
                                                     const char* name, const char* title)
  : RooUnfoldResponseT(name, title)
{
  //! constructor - simple 1D case with same binning, measured vs truth
  Setup(nb, xlo, xhi);
}


RooUnfoldResponse::RooUnfoldResponse(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi,
                                                     const char* name, const char* title)
  : RooUnfoldResponseT(name, title)
{
  //! constructor - simple 1D case
  Setup(nm, mlo, mhi, nt, tlo, thi);
}


RooUnfoldResponse::RooUnfoldResponse(const TH1* measured, const TH1* truth, const TH2* response,
                                                     const char* name, const char* title, bool overflow)
  : RooUnfoldResponseT(name, title)
{
  //! RooUnfoldResponse constructor - create from already-filled histograms
  //! "response" gives the response matrix, measured X truth.
  //! "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  //! but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  //! in "truth" for unmeasured events (inefficiency).
  //! "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  //! to indicate, respectively, no fakes and/or no inefficiency.
  this->_overflow = overflow;
  Setup(measured, truth, response);
}


RooUnfoldResponse::RooUnfoldResponse(const TH1* measured, const TH1* truth,
                                                     const char* name, const char* title, bool overflow)
  : RooUnfoldResponseT(name, title)
{
  //! RooUnfoldResponse constructor - measured and truth only used for shape
  this->_overflow = overflow;
  Setup(measured, truth);
}


void RooUnfoldResponse::Add(const RooUnfoldResponse& rhs)
{
  //! Add another RooUnfoldResponse, accumulating contents
  if (_res == 0) {
    Setup(rhs);
    return;
  }
  assert (GetDimensionMeasured()==rhs.GetDimensionMeasured());
  assert (GetDimensionTruth()==rhs.GetDimensionTruth());
  assert (_mes != 0 && rhs._mes != 0);
  assert (_fak != 0 && rhs._fak != 0);
  assert (_tru != 0 && rhs._tru != 0);
  assert (_res != 0 && rhs._res != 0);
  if (Cached()) ClearCache();
  _mes->Add(rhs._mes);
  _tru_m->Add(rhs._tru_m);
  _mestru->Add(rhs._mestru);
  _fak->Add(rhs._fak);
  _tru->Add(rhs._tru);
  _res->Add(rhs._res);
}


 Long64_t
RooUnfoldResponse::Merge (TCollection* others)
{
  //! Add all RooUnfoldResponse<class TH1, class TH2> objects in the collection to this one.
  //! This allows merging with hadd and TFileMerger.
  for (TIter it= others; TObject* o= it();) {
    if (RooUnfoldResponse* other= dynamic_cast<RooUnfoldResponse*>(o))
      Add (*other);
  }
  return Long64_t(::entries(_res));
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(Int_t nb, Double_t xlo, Double_t xhi)
{
  //! constructor -  simple 1D case with same binning, measured vs truth
  return Setup(nb, xlo, xhi, nb, xlo, xhi);
}




 RooUnfoldResponse&
RooUnfoldResponse::Setup(Int_t nm, Double_t mlo, Double_t mhi, Int_t nt, Double_t tlo, Double_t thi)
{
  //! set up simple 1D case
  Reset();
  _mes= createHist<TH1>("measured", "Measured",   Variable<TH1>(nm, mlo, mhi,"xm"));
  _fak= createHist<TH1>("fakes",    "Fakes",      Variable<TH1>(nm, mlo, mhi,"xm"));
  _tru= createHist<TH1>("truth",    "Truth",      Variable<TH1>(nt, tlo, thi,"xt"));
  _tru_m= createHist<TH1>("truth_m","truth_m",Variable<TH1>(nm, tlo, thi,"xt"));
  _mestru= createHist<TH1>("meastru","Meastru",Variable<TH1>(nm, tlo, thi,"xt"));
  _res= createHist<TH2>("response", "Response", Variable<TH2>(nm, mlo, mhi, "xm"), Variable<TH2>(nt, tlo, thi, "xt"));
  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(const TH1* measured, const TH1* truth)
{
  //! set up - measured and truth only used for shape
  Reset();
  _mes= clone(measured);
  _fak= clone(measured);
  _tru= clone(truth);
  _tru_m= clone(measured);
  _mestru=clone(measured);
  _mes->Reset();
  _fak->Reset();
  _tru->Reset();
  _tru_m->Reset();
  _mestru->Reset();
  if (_overflow && (GetDimensionMeasured() > 1 || GetDimensionTruth() > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }
  

  SetNameTitleDefault();

  Double_t measbinedges[_mes->GetNbinsX()+1]; 
  Double_t truthbinedges[_tru->GetNbinsX()+1]; 
  
  for (int i = 0; i < _mes->GetNbinsX()+1; i++){
    measbinedges[i] = _mes->GetBinLowEdge(i+1);
  }
  for (int i = 0; i < _tru->GetNbinsX()+1; i++){
    truthbinedges[i] = _tru->GetBinLowEdge(i+1);
  }
  _res = new TH2D(GetName(), GetTitle(), _mes->GetNbinsX(), measbinedges, _tru->GetNbinsX(), truthbinedges);
  _res->Sumw2();

  return *this;
}

RooUnfoldResponse&
RooUnfoldResponse::Setup(const TH1* measured, const TH1* truth, const TH2* response)
{
  //! Set up from already-filled histograms.
  //! "response" gives the response matrix, measured X truth.
  //! "measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
  //! but with additional entries in "measured" for measurements with no corresponding truth (fakes/background) and
  //! in "truth" for unmeasured events (inefficiency).
  //! "measured" and/or "truth" can be specified as 0 (1D case only) or an empty histograms (no entries) as a shortcut
  //! to indicate, respectively, no fakes and/or no inefficiency.

  Reset();
  _res= clone(response);
  if (measured) {
    _mes= clone(measured);
  } else {
    _mes= createHist<TH1>("measured", "Measured", Variable<TH1>(nBins(response,RooUnfolding::X), 0.0, 1.0, "xm"));
  }

  if (truth) {
    _tru= clone(truth);
  } else {
    _tru= createHist<TH1>("truth",    "Truth",    Variable<TH1>(nBins(response,RooUnfolding::Y), 0.0, 1.0, "xt"));
  }
  if (_overflow && (GetDimensionMeasured() > 1 || GetDimensionTruth() > 1)) {
    cerr << "UseOverflow setting ignored for multi-dimensional distributions" << endl;
    _overflow= 0;
  }

  double mlow = measured->GetXaxis()->GetXmin();
  double mhigh = measured->GetXaxis()->GetXmax();
  double tlow = truth->GetXaxis()->GetXmin();
  double thigh = truth->GetXaxis()->GetXmax();

  _tru_m= createHist<TH1>("truth_m","truth_m",Variable<TH1>(GetNbinsMeasured(), tlow, thigh,"xt"));
  _mestru= createHist<TH1>("meastru","Meastru",Variable<TH1>(GetNbinsMeasured(), tlow, thigh,"xt"));
  _fak= createHist<TH1>("fakes","Fakes",Variable<TH1>(GetNbinsMeasured(), mlow, mhigh,"xm"));

  if (GetNbinsMeasured() != nBins(_res,RooUnfolding::X) || GetNbinsTruth() != nBins(_res,RooUnfolding::Y)) {
    cerr << "Warning: RooUnfoldResponse<class TH1, class TH2> measured X truth is " << GetNbinsMeasured() << " X " << GetNbinsTruth()
         << ", but matrix is " << nBins(_res,RooUnfolding::X)<< " X " << nBins(_res,RooUnfolding::Y) << endl;
  }

  if (!measured || ::entries(_mes) == 0.0) {
    // Similar to _res->ProjectionX() but without stupid reset of existing histograms
    // Always include under/overflows in sum of truth.
    projectX(_res,_mes,_overflow);
  } else {
    // Fill fakes from the difference of _mes - _res->ProjectionX()
    // Always include under/overflows in sum of truth.
    subtractProjectX(_res,_mes,_fak,_overflow);
  }

  if (!truth || ::entries(_tru) == 0.0) {
    // similar to _res->ProjectionY() but without stupid reset of existing histograms
    // Always include under/overflows in sum of measurements.
    projectY(_res,_tru,_overflow);
  }

  SetNameTitleDefault();
  return *this;
}


 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t xt, Double_t w)
{

  //! Fill 1D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (GetDimensionMeasured()==1 && GetDimensionTruth()==1);
  if (Cached()) ClearCache();
  Int_t bin_m = fill(_mes,xr,w);
  Int_t bin_t = fill(_tru,xt,w);
  if (GetNbinsMeasured() != GetNbinsTruth()){
    bin_t = fill(_tru_m,xt,w);
  }

  // Fill a histogram with events that are reconstructed
  // in the same truth bin.
  if (bin_m == bin_t){
    fill(_mestru,xr,w);
  }
  return fill(_res,xr,xt,w);
}

 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t xt, Double_t yt, Double_t w)
{
  //! Fill 2D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (GetDimensionMeasured()==2 && GetDimensionTruth()==2);
  if (Cached()) ClearCache();
  fill((TH2*)_mes,xr, yr, w);
  fill((TH2*)_tru,xt, yt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt)+1,RooUnfolding::Y), w);
}

 Int_t
RooUnfoldResponse::Fill (Double_t xr, Double_t yr, Double_t zr, Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  //! Fill 3D Response Matrix
  assert (_mes != 0 && _tru != 0);
  assert (GetDimensionMeasured()==3 && GetDimensionTruth()==3);
  if (Cached()) ClearCache();
  fill(_mes,xr, yr, zr, w);
  fill(_tru,xt, yt, zt, w);
  return fill(_res,binCenter(_res,findBin (_mes, xr, yr, zr)+1,RooUnfolding::X),binCenter(_res,findBin (_tru, xt, yt, zt)+1,RooUnfolding::Y), w);  
}


 Int_t
RooUnfoldResponse::FindBin(const TH1* h, Double_t x, Double_t y, Double_t z)
{
  //! Get vector index (0..nx*ny*nz-1) for bin containing (x,y,z) coordinates
  Int_t nx=   nBins(h,RooUnfolding::X);
  Int_t ny=   nBins(h,RooUnfolding::Y);
  Int_t nz=   nBins(h,RooUnfolding::Z);
  Int_t binx= findBin(h,x,RooUnfolding::X) - 1;
  if (binx <  0)  return -1;
  if (binx >= nx) return nx*ny*nz;
  Int_t biny= findBin(h,y,RooUnfolding::Y) - 1;
  if (biny <  0)  return -1;
  if (biny >= ny) return nx*ny*nz;
  Int_t binz= findBin(h,z,RooUnfolding::Z) - 1;
  if (binz <  0)  return -1;
  if (binz >= nz) return nx*ny*nz;
  return binx + nx*(biny + ny*binz);
}

Int_t
RooUnfoldResponse::Miss1D (Double_t xt, Double_t w)
{
  //! Fill missed event (not reconstructed due to detection inefficiencies) into 1D Response Matrix (with weight)
  assert (_tru != 0);
  assert (GetDimensionTruth()==1);
  if (Cached()) ClearCache();
  return fill(_tru, xt, w);
}

 Int_t
RooUnfoldResponse::Miss2D (Double_t xt, Double_t yt, Double_t w)
{
  //! Fill missed event (not reconstructed due to detection inefficiencies) into 2D Response Matrix (with weight)
  assert (_tru != 0);
  assert (GetDimensionTruth()==2);
  if (Cached()) ClearCache();
  return fill(_tru, xt, yt, w);
}

 Int_t
RooUnfoldResponse::Miss (Double_t xt, Double_t yt, Double_t zt, Double_t w)
{
  //! Fill missed event (not reconstructed due to detection inefficiencies) into 3D Response Matrix
  assert (_tru != 0);
  assert (GetDimensionTruth()==3);
  if (Cached()) ClearCache();
  return fill(_tru, xt, yt, zt, w);
}

 Int_t
RooUnfoldResponse::Fake1D (Double_t xr, Double_t w)
{
  //! Fill fake event (reconstructed event with no truth) into 1D Response Matrix (with weight)
  assert (_fak != 0 && _mes != 0);
  assert (GetDimensionMeasured()==1);
  if (Cached()) ClearCache();
  fill(_mes,xr, w);
  return fill(_fak, xr, w);
}

 Int_t
RooUnfoldResponse::Fake2D (Double_t xr, Double_t yr, Double_t w)
{
  //! Fill fake event (reconstructed event with no truth) into 2D Response Matrix (with weight)
  assert (_mes != 0);
  assert (GetDimensionMeasured()==2);
  if (Cached()) ClearCache();
  fill(_fak, xr, yr, w);
  return fill(_mes, xr, yr, w);
}

 Int_t
RooUnfoldResponse::Fake (Double_t xr, Double_t yr, Double_t zr, Double_t w)
{
  //! Fill fake event (reconstructed event with no truth) into 3D Response Matrix
  assert (_mes != 0);
  assert (GetDimensionMeasured()==3);
  if (Cached()) ClearCache();
  fill(_mes, xr, yr, zr, w);
  return fill(_fak, xr, yr, zr, w);
}


Long64_t RooUnfoldResponse::FakeEntries() const
{
  //! Return number of fake entries
  return _fak ? entries(_fak) : 0.0;
}

ClassImp (RooUnfoldResponse);


#ifndef NOROOFIT

#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
#include "RooAddition.h"
#include "RooConstVar.h"
#include "RooHistFunc.h"



namespace {
  template<class T> void makeList(const RooAbsCollection* c,std::vector<T*>& vars,TPRegexp* re = 0){
    RooFIter itr(c->fwdIterator());
    RooRealVar* obj = NULL;
    while((obj = (RooRealVar*)itr.next())){
      if(!obj->InheritsFrom(RooRealVar::Class())) continue;
      if(!re || re->Match(obj->GetName()))
        vars.push_back(obj);
    }
  }
  
  template<class T> std::vector<T*> makeList(const RooAbsCollection* c){
    std::vector<RooAbsArg*> vars;
    makeList<T>(c,vars);
    return vars;
  }
}

template <> RooUnfolding::RooFitHist*
RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::HresponseNoOverflow() const
{
  return const_cast<RooUnfolding::RooFitHist*>(Hresponse());
}



template class RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;

RooFitUnfoldResponse::RooFitUnfoldResponse(const RooUnfoldResponseT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* res) :
  RooUnfoldResponseT(res->GetName(),res->GetTitle(),static_cast<const RooFitUnfoldResponse*>(res)->_res->clone(),static_cast<const RooFitUnfoldResponse*>(res)->_tru->clone(),static_cast<const RooFitUnfoldResponse*>(res)->_mes->clone(),false,static_cast<const RooFitUnfoldResponse*>(res)->_density)  
{
}

RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooUnfolding::RooFitHist* response, RooUnfolding::RooFitHist* truth, RooUnfolding::RooFitHist* reco, bool density) : 
  RooUnfoldResponseT(name,title,response,truth,reco,false,density)
{
  //! inherited constructor
}


RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, RooRealVar* obs_truth, RooRealVar* obs_reco, bool density){
  if(!truth->dependsOn(*obs_truth)) throw std::runtime_error("truth histogram does not depend on truth observable!");
  if(!reco->dependsOn(*obs_reco)) throw std::runtime_error("reco histogram does not depend on reco observable!");
  if(!response->dependsOn(*obs_truth)) throw std::runtime_error("response histogram does not depend on truth observable!");
  if(!response->dependsOn(*obs_reco)) throw std::runtime_error("response histogram does not depend on reco observable!");
  
  TPRegexp gamma("gamma_stat_.*");

  std::vector<RooRealVar*> allvars;
  RooArgSet* c_reco_vars = reco->getVariables();
  makeList<RooRealVar>(c_reco_vars,allvars,&gamma);
  delete c_reco_vars;
  if(fakes){
    RooArgSet* c_fake_vars = fakes->getVariables();
    makeList<RooRealVar>(c_fake_vars,allvars,&gamma);
    delete c_fake_vars;
  }
  RooArgSet* c_truth_vars = truth->getVariables();
  makeList<RooRealVar>(c_truth_vars,allvars,&gamma);
  delete c_truth_vars;

  this->_mes = new RooFitHist(reco,obs_reco,allvars);
  this->_fak = fakes ? new RooFitHist(fakes,obs_reco,allvars) : 0;
  this->_tru = new RooFitHist(truth,obs_truth,allvars);
  this->_res = new RooFitHist(response,obs_truth,obs_reco,allvars);
  this->_overflow = 0;
  this->_density = density;
}

RooFitUnfoldResponse::RooFitUnfoldResponse(const char* name, const char* title, RooAbsReal* response, RooAbsReal* truth, RooAbsReal* reco, RooAbsReal* fakes, const RooAbsCollection* observables, bool density)
   : RooUnfoldResponseT(name,title) {
  //! standard constructor

  RooArgSet obs_in;
  obs_in.add(*observables);
  RooArgSet* obs = response->getObservables(&obs_in);
  RooArgSet* obsset_reco = reco->getObservables(&obs_in);
  RooArgSet* obsset_truth = truth->getObservables(&obs_in);    
 
  if(obs->getSize() != 2) throw std::runtime_error(TString::Format("unsupported dimensionality for response: %d",obs->getSize()).Data());
  if(obsset_reco->getSize() != 1) throw std::runtime_error(TString::Format("unsupported dimensionality for reco: %d",obsset_reco->getSize()).Data());
  if(obsset_truth->getSize() != 1) throw std::runtime_error(TString::Format("unsupported dimensionality for truth: %d",obsset_truth->getSize()).Data());

  RooRealVar* obs_reco  = (RooRealVar*)obsset_reco->first();
  RooRealVar* obs_truth = (RooRealVar*)obsset_truth->first();
  
  this->_mes = new RooFitHist(reco,obs_reco);
  this->_fak = fakes ? new RooFitHist(fakes,obs_reco) : 0;
  this->_tru = new RooFitHist(truth,obs_truth);
  this->_res = new RooFitHist(response,::makeList<RooAbsArg>(obs));
  this->_overflow = 0;
  this->_density = density;
}


RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistSum(RooAbsReal* a, RooAbsReal* b, double ca, double cb){
  TString name(TString::Format("%s_minus_%s",a->GetName(),b->GetName()));
  RooArgList funcs;
  funcs.add(*a);
  funcs.add(*b);
  RooConstVar* cA = new RooConstVar("cA","cA",ca);
  RooConstVar* cB = new RooConstVar("cB","cB",cb);  
  RooArgList coefs;
  coefs.add(*cA);
  coefs.add(*cB);
  RooAddition* add = new RooAddition(name.Data(),name.Data(),funcs,coefs);
  return this->makeHist(add);
}

RooHistFunc* RooFitUnfoldResponse::makeHistFunc(RooDataHist* dhist){
  if(!dhist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(dhist->get()->find(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(dhist->get()->find(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooDataHist '%s', does not seem to depend on any known observable!",dhist->GetName()).Data());
  }
  return RooUnfolding::makeHistFunc(dhist,v);
}

RooHistFunc* RooFitUnfoldResponse::makeHistFuncMeasured(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    v.push_back(this->_mes->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistFunc(dhist,v);
}
RooHistFunc* RooFitUnfoldResponse::makeHistFuncTruth(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_tru->dim(); ++i){
    v.push_back(this->_tru->obs(i));
  }

  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistFunc(dhist,v);
}

RooAbsPdf* RooFitUnfoldResponse::makeHistPdf(RooDataHist* dhist){
  if(!dhist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(dhist->get()->find(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(dhist->get()->find(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooDataHist '%s', does not seem to depend on any known observable!",dhist->GetName()).Data());
  }
  return RooUnfolding::makeHistPdf(dhist,v);
}

RooAbsPdf* RooFitUnfoldResponse::makeHistPdfMeasured(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    v.push_back(this->_mes->obs(i));
  }

  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistPdf(dhist,v);
}
RooAbsPdf* RooFitUnfoldResponse::makeHistPdfTruth(const TH1* hist){
  if(!hist) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_tru->dim(); ++i){
    v.push_back(this->_tru->obs(i));
  }
  RooDataHist* dhist = RooUnfolding::convertTH1(hist,v,this->_overflow,this->_density);
  return RooUnfolding::makeHistPdf(dhist,v);
}

RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistMeasured(const TH1* hist){
  return this->makeHist(this->makeHistFuncMeasured(hist));
}
RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHistTruth(const TH1* hist){
  return this->makeHist(this->makeHistFuncTruth(hist));  
}


RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHist(RooAbsReal* object){
  if(!object) return NULL;
  std::vector<RooAbsArg*> v;
  for(size_t i=0; i<this->_mes->dim(); ++i){
    if(object->dependsOn(*this->_mes->obs(i))) v.push_back(this->_mes->obs(i));
  }
  for(size_t i=0; i<this->_tru->dim(); ++i){
    if(object->dependsOn(*this->_tru->obs(i))) v.push_back(this->_tru->obs(i));
  }
  if(v.size() == 0){
    throw std::runtime_error(TString::Format("unable to construct histogram from RooAbsReal '%s', does not seem to depend on any known observable!",object->GetName()).Data());
  }
  return new RooUnfolding::RooFitHist(object,v);
}

RooUnfolding::RooFitHist* RooFitUnfoldResponse::makeHist(RooDataHist* object){
  if(!object) return NULL;
  std::vector<RooAbsArg*> v;
  TPRegexp obs("obs_.*");
  makeList<RooAbsArg>(object->get(),v,&obs);
  return new RooUnfolding::RooFitHist(object,v);
}

ClassImp (RooFitUnfoldResponse);

#endif
