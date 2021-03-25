/*! \class RooUnfoldTUnfold
<p>Uses the unfolding method implemented in ROOT's <a href="http://root.cern.ch/root/html/TUnfold.html">TUnfold</a> class
<p>Only included in ROOT versions 5.22 and higher
<p>Only able to reconstruct 1 dimensional distributions
<p>Can account for bin migration and smearing
<p>Errors come as a full covariance matrix.
<p>Will sometimes warn of "unlinked" bins. These are bins with  entries and do not effect the results of the unfolding
<p>Regularisation parameter can be either optimised internally by plotting log10(chi2 squared) against log10(tau). The 'kink' in this curve is deemed the optimum tau value. This value can also be set manually (FixTau)
<p>The latest version (TUnfold 15 in ROOT 2.27.04) will not handle plots with an additional underflow bin. As a result overflows must be turned off
if v15 of TUnfold is used. ROOT versions 5.26 or below use v13 and so should be safe to use overflows.</ul>
 */

#include "RooUnfoldTUnfold.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#endif

#include <iostream>
#include <iomanip>
#include <math.h>

#include "TNamed.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#ifndef NOTUNFOLDSYS
#include "TUnfoldSys.h"
#else
#include "TUnfold.h"
#endif
#include "TGraph.h"
#include "TSpline.h"

#include "RooUnfoldHelpers.h"
#include "RooUnfoldResponse.h"

using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT (const RooUnfoldTUnfoldT& rhs)
: RooUnfoldT<Hist,Hist2D> (rhs)
{
  // Copy constructor.
  Init();
  CopyData (rhs);
}

template<class Hist,class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, TUnfold::ERegMode reg, Bool_t handleFakes, const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title),_reg_method(reg), _handleFakes(handleFakes)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Double_t tau, TUnfold::ERegMode reg, Bool_t handleFakes, const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title),_reg_method(reg), _handleFakes(handleFakes)
{
  // Constructor with response matrix object and measured unfolding input histogram.
  // This one uses a fixed regularisation parameter, tau.
  Init();
  FixTau(tau);
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::Destroy()
{
  delete _unf;
  delete _lCurve;
  delete _logTauX;
  delete _logTauY;
  delete _logTauSURE;
  delete _df_chi2A;
}

// template<class Hist,class Hist2D>RooUnfoldTUnfoldT*
// RooUnfoldTUnfoldT<Hist,Hist2D>::Clone (const char* newname) const
// {
//   // Clones object
//   RooUnfoldTUnfoldT* unfold= new RooUnfoldTUnfoldT(*this);
//   if (newname && strlen(newname)) unfold->SetName(newname);
//   return unfold;
// }

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::Reset()
{
  // Resets all values
  Destroy();
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}


template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::Assign (const RooUnfoldTUnfoldT& rhs)
{
  RooUnfoldT<Hist,Hist2D>::Assign (rhs);
  CopyData (rhs);
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::CopyData (const RooUnfoldTUnfoldT& rhs)
{
  tau_set=rhs.tau_set;
  _tau=rhs._tau;
  _reg_method=rhs._reg_method;
  _lCurve  = (rhs._lCurve  ? dynamic_cast<TGraph*> (rhs._lCurve ->Clone()) : 0);
  _logTauX = (rhs._logTauX ? dynamic_cast<TSpline*>(rhs._logTauX->Clone()) : 0);
  _logTauY = (rhs._logTauY ? dynamic_cast<TSpline*>(rhs._logTauY->Clone()) : 0);
  _logTauSURE = (rhs._logTauSURE  ? dynamic_cast<TGraph*> (rhs._logTauSURE ->Clone()) : 0);
  _df_chi2A = (rhs._df_chi2A  ? dynamic_cast<TGraph*> (rhs._df_chi2A ->Clone()) : 0);
}


template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::Init()
{
  tau_set=false;
  _tau=-1e30;
  _unf=0;
  _lCurve = 0;
  _logTauX = 0;
  _logTauY = 0;
  _logTauSURE = 0;
  _df_chi2A = 0;
  GetSettings();
}

template<class Hist,class Hist2D>TUnfold* 
RooUnfoldTUnfoldT<Hist,Hist2D>::Impl()
{
  // Return TUnfold object used to implement the unfolding for RooUnfoldTUnfold.
  return _unf;
}

template<class Hist,class Hist2D>const TGraph* 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetLCurve() const
{
  // If an L curve scan has been done (tau not fixed by user),
  // return the L curve as graph
  return _lCurve;
}

template<class Hist,class Hist2D>const TSpline* 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetLogTauX() const
{
  // If an L curve scan has been done (tau not fixed by user),
  // return the spline of x-coordinates vs tau for the L curve
  return _logTauX;
}

template<class Hist,class Hist2D>const TSpline* 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetLogTauY() const
{
  // If an L curve scan has been done (tau not fixed by user),
  // return the spline of y-coordinates vs tau for the L curve
  return _logTauY;
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::Unfold() const
{
  // Does the unfolding. Uses the optimal value of the unfolding parameter unless a value has already been set using FixTau
  if (this->_nm<this->_nt)     std::cerr << "Warning: fewer measured bins than truth bins. TUnfold may not work correctly." << std::endl;
  if (this->_covMes) std::cerr << "Warning: TUnfold does not account for bin-bin correlations on measured input"    << std::endl;

  Bool_t oldstat= TH1::AddDirectoryStatus();
  TH1::AddDirectory (kFALSE);

  TVectorD vmeas = this->Vmeasured();
  TVectorD verr = this->Emeasured();
  const TVectorD& vtruth(this->_res->Vtruth());
  TMatrixD mres = this->_res->Mresponse(false);
  
  // The added empty bins are used to store inefficiencies.
  RooUnfolding::addEmptyBins(vmeas);
  RooUnfolding::addEmptyBins(verr);
  RooUnfolding::addEmptyBins(mres);  

  TH1::AddDirectory (oldstat);
    
  // Add inefficiencies to measured overflow bin
  for (Int_t j= 1; j<=this->_nt; j++) {
    Double_t ntru= 0.0;
    for (Int_t i= 1; i<=this->_nm; i++) {
      ntru += mres[i][j];
    }
    mres[this->_nm + 1][j] = vtruth[j - 1] - ntru;
  }
    
  // Subtract fakes from measured distribution
  if (this->_res->HasFakes() && _handleFakes) {
    TVectorD fakes= this->_res->Vfakes();
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) std::cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << std::endl;
    for (Int_t i = 1; i<=this->_nm; i++){
      vmeas[i] = vmeas[i] - (fac*fakes[i - 1]);
    }
  }

  Int_t ndim = dim(this->_meas);
  TUnfold::ERegMode reg= _reg_method;
  
  if (ndim == 2 || ndim == 3) reg= TUnfold::kRegModeNone;  // set explicitly

#ifndef NOTUNFOLDSYS
  if (this->_dosys){
    _unf= new TUnfoldSysV17(&mres,TUnfold::kHistMapOutputVert,reg);
  } else {
#endif
    _unf= new TUnfoldV17(&mres,TUnfold::kHistMapOutputVert,reg);
  }
  
  if        (ndim == 2) {
    Int_t nx= nBins(this->_meas,X), ny= nBins(this->_meas,Y);
    _unf->RegularizeBins2D (0, 1, nx, nx, ny, _reg_method);
  } else if (ndim == 3) {
    Int_t nx= nBins(this->_meas,X), ny= nBins(this->_meas,Y), nz= nBins(this->_meas,Z), nxy= nx*ny;
    for (Int_t i= 0; i<nx; i++) {
      _unf->RegularizeBins2D (    i, nx, ny, nxy, nz, _reg_method);
    }
    for (Int_t i= 0; i<ny; i++) {
      _unf->RegularizeBins2D ( nx*i,  1, nx, nxy, nz, _reg_method);
    }
    for (Int_t i= 0; i<nz; i++) {
      _unf->RegularizeBins2D (nxy*i,  1, nx,  nx, ny, _reg_method);
    }
  }

  Int_t nScan=30;
  
  // use automatic L-curve scan: start with taumin=taumax=0.0
  Double_t tauMin=0.0;
  Double_t tauMax=0.0;
 
  // this method scans the parameter tau and finds the kink in the L curve
  // finally, the unfolding is done for the best choice of tau
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,23,0)  /* TUnfold v6 (included in ROOT 5.22) didn't have setInput return value */

  Int_t stat= _unf->SetInput(&vmeas,&verr);
  if(stat>=10000) {
    std::cerr<<"Unfolding result may be wrong: " << stat/10000 << " unconstrained output bins\n";
  }
#else
  _unf->SetInput(&vmeas,&verr);
#endif
  _unf->SetConstraint(TUnfold::kEConstraintArea);

  if (!tau_set){
    delete _lCurve;  _lCurve  = 0;
    delete _logTauX; _logTauX = 0;
    delete _logTauY; _logTauY = 0;
    delete _logTauSURE; _logTauSURE = 0;
    delete _df_chi2A; _df_chi2A = 0;
    //Int_t bestPoint = _unf->ScanSURE(nScan,tauMin,tauMax,&_logTauSURE,&_df_chi2A,&_lCurve);
    Int_t bestPoint = _unf->ScanLcurve(nScan,tauMin,tauMax,&_lCurve,&_logTauX,&_logTauY);
    _tau=_unf->GetTau();  // save value, even if we don't use it unless tau_set
    std::cout <<"Lcurve scan chose tau= "<<_tau<<std::endl<<" at point "<<bestPoint<<std::endl;
    
    // Added the undersmoothing method developed by Junhyung Lyle Kim and Mikael Kuusela
    //_tau = _unf->UndersmoothTau(_tau, 0.01, 1000);
  }

  _unf->DoUnfold(_tau);

  TH1F reco("_cache._rec","reconstructed dist",this->_nt,0.0,this->_nt);
  _unf->GetOutput(&reco);
  this->_cache._rec.ResizeTo (this->_nt);
  for (int i=0;i<this->_nt;i++){
    this->_cache._rec(i)=(reco.GetBinContent(i + 1));
  }

  if (this->_verbose>=2) {
    const Hist* train1d = this->_res->Hmeasured();
    const Hist* truth1d = this->_res->Htruth();
    printTable (std::cout, h2v(this->_meas), h2v(train1d), h2v(truth1d), this->_cache._rec);
  }

  this->_cache._unfolded= true;
  this->_cache._haveCov=  false;
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::GetCov() const
{
  //! Get covariance matrix
  if (!_unf) return;
  TH2* ematrix= new TH2D ("ematrix","error matrix", this->_nt, 0.0, this->_nt, this->_nt, 0.0, this->_nt);
  if (this->_dosys!=2) _unf->GetEmatrix (ematrix);
  if (this->_dosys) {
#ifndef NOTUNFOLDSYS
    TUnfoldSys* unfsys= dynamic_cast<TUnfoldSys*>(_unf);
    if (unfsys)
      unfsys->GetEmatrixSysUncorr (ematrix,0,kFALSE);
    else
#endif
      std::cerr << "Did not use TUnfoldSys, so cannot calculate systematic errors" << std::endl;
  }
  this->_cache._cov.ResizeTo (this->_nt,this->_nt);
  for (Int_t i= 0; i<this->_nt; i++) {
    for (Int_t j= 0; j<this->_nt; j++) {
      this->_cache._cov(i,j)= ematrix->GetBinContent(i+1,j+1);
    }
  }
  delete ematrix;
  this->_cache._haveCov= true;
}


template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::FixTau(Double_t tau)
{
  // Fix regularisation parameter to a specified value
  _tau=tau;
  tau_set=true;
  this->ResetUnfold();
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::SetRegMethod(TUnfold::ERegMode regmethod)
{
  /*
    Specifies the regularisation method:

      regemthod setting             regularisation
      ===========================   ==============
      TUnfold::kRegModeNone         none
      TUnfold::kRegModeSize         minimize the size of (x-x0)
      TUnfold::kRegModeDerivative   minimize the 1st derivative of (x-x0)
      TUnfold::kRegModeCurvature    minimize the 2nd derivative of (x-x0)
   */
  _reg_method=regmethod;
  this->ResetUnfold();
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::OptimiseTau()
{
  // Choose optimal regularisation parameter
  tau_set=false;
  this->ResetUnfold();
}

template<class Hist,class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::GetSettings() const
{
    this->_cache._minparm=0;
    this->_cache._maxparm=1;
    this->_cache._stepsizeparm=1e-2;
    this->_cache._defaultparm=2;
}


template<class Hist, class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT()
  : RooUnfoldT<Hist,Hist2D>()
{
  // Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::RooUnfoldTUnfoldT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  // Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist, class Hist2D>RooUnfoldTUnfoldT<Hist,Hist2D>& 
RooUnfoldTUnfoldT<Hist,Hist2D>::operator= (const RooUnfoldTUnfoldT<Hist,Hist2D>& rhs)
{
  // Assignment operator for copying RooUnfoldTUnfold settings.
  Assign(rhs);
  return *this;
}

template<class Hist, class Hist2D>
RooUnfoldTUnfoldT<Hist,Hist2D>::~RooUnfoldTUnfoldT()
{
  Destroy();
}


template<class Hist, class Hist2D>void
RooUnfoldTUnfoldT<Hist,Hist2D>::SetRegParm(Double_t parm)
{
  // Set regularisation parameter (tau)
  FixTau(parm);
}

template<class Hist, class Hist2D>Double_t 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetTau() const
{
  // Return regularisation parameter (tau)
  return _tau;
}

template<class Hist, class Hist2D>Double_t 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetRegParm() const
{
  // Return regularisation parameter (tau)
  return _tau;
}

template<class Hist, class Hist2D>TUnfold::ERegMode 
RooUnfoldTUnfoldT<Hist,Hist2D>::GetRegMethod() const
{
  return _reg_method;
}


template class RooUnfoldTUnfoldT<TH1,TH2>;
ClassImp (RooUnfoldTUnfold)

#ifndef NOROOFIT
template class RooUnfoldTUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldTUnfold)
#endif
