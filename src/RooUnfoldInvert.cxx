/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2010–2025 CERN and the authors’ respective research institutions
 * Please refer to the CONTRIBUTORS file for details.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * END ROOUNFOLD COPYRIGHT
 */
/*===========================================================================*/

/* \class RooUnfoldInvertT
<p>The simplest method of unfolding works by simply inverting the response matrix.</p>
<p>This is not accurate for small matrices and produces inaccurate unfolded distributions.</p>
<p>The inversion method is included largely to illustrate the necessity of a more effective method of unfolding</p>
*/

#include "RooUnfoldInvert.h"

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"

using namespace RooUnfolding;

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::RooUnfoldInvertT(const RooUnfoldInvertT<Hist, Hist2D> &rhs)
   : RooUnfoldT<Hist, Hist2D>(rhs)
{
   // Copy constructor.
   Init();
}

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::RooUnfoldInvertT(const RooUnfoldResponseT<Hist, Hist2D> *res, const Hist *meas,
                                                 const char *name, const char *title)
   : RooUnfoldT<Hist, Hist2D>(res, meas, name, title)
{
   // Constructor with response matrix object and measured unfolding input histogram.
   Init();
}

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::~RooUnfoldInvertT()
{
   delete _svd;
   delete _resinv;
}

template <class Hist, class Hist2D>
void RooUnfoldInvertT<Hist, Hist2D>::Init()
{
   _svd = 0;
   _resinv = 0;
   this->GetSettings();
}

template <class Hist, class Hist2D>
void RooUnfoldInvertT<Hist, Hist2D>::Reset()
{
   delete _svd;
   delete _resinv;
   Init();
   RooUnfoldT<Hist, Hist2D>::Reset();
}

template <class Hist, class Hist2D>
TDecompSVD *RooUnfoldInvertT<Hist, Hist2D>::Impl()
{
   return _svd;
}

template <class Hist, class Hist2D>
RooUnfolding::Algorithm RooUnfoldInvertT<Hist, Hist2D>::GetAlgorithm() const
{
   //! return the unfolding algorithm used
   return kInvert;
}

template <class Hist, class Hist2D>
void RooUnfoldInvertT<Hist, Hist2D>::Unfold() const
{
   int nm = this->response()->Vmeasured().GetNrows();
   int nt = this->response()->Vtruth().GetNrows();

   this->_cache._rec.ResizeTo(nm);
   this->_cache._rec = this->Vmeasured();

   if (this->_res->HasFakes()) {
      TVectorD fakes = this->_res->Vfakes();
      Double_t fac = this->_res->Vmeasured().Sum();
      if (fac != 0.0)
         fac = this->Vmeasured().Sum() / fac;
      if (this->_verbose >= 1)
         std::cout << "Subtract " << fac * fakes.Sum() << " fakes from measured distribution" << std::endl;
      fakes *= fac;
      this->_cache._rec -= fakes;
   }

   Bool_t ok;
   ok = InvertResponse();
   if (nt > nm) {
      if (ok)
         this->_cache._rec *= *_resinv;
   } else {
      ok = _svd->Solve(this->_cache._rec);
   }

   this->_cache._rec.ResizeTo(nt);
   if (!ok) {
      std::cerr << "Response matrix Solve failed" << std::endl;
      return;
   }

   this->_cache._unfolded = true;
   this->_cache._haveCov = false;
}

template <class Hist, class Hist2D>
void RooUnfoldInvertT<Hist, Hist2D>::GetCov() const
{
   //! Get covariance matrix

   int nt = this->response()->Vtruth().GetNrows();
   if (!InvertResponse())
      return;
   this->_cache._cov.ResizeTo(nt, nt);
   ABAT(*_resinv, this->GetMeasuredCov(), this->_cache._cov);
   this->_cache._haveCov = true;
}

template <class Hist, class Hist2D>
const TMatrixD &RooUnfoldInvertT<Hist, Hist2D>::InverseResponse() const
{
   if (!this->InvertResponse()) {
      throw std::runtime_error("unable to invert response!");
   }
   return *_resinv;
}

template <class Hist, class Hist2D>
Bool_t RooUnfoldInvertT<Hist, Hist2D>::InvertResponse() const
{
   if (_resinv)
      return true;

   int nm = this->response()->Vmeasured().GetNrows();
   int nt = this->response()->Vtruth().GetNrows();

   TMatrixD res(this->_res->Mresponse(true));
   if (nt > nm) {
      TMatrixD resT(TMatrixD::kTransposed, res);
      _svd = new TDecompSVD(resT);
      delete _resinv;
      _resinv = 0;
   } else {
      _svd = new TDecompSVD(res);
   }
   double c = _svd->Condition();
   if (c < 0)
      std::cout << "WARNING: Response matrix is ill-conditioned. TDecompSVD condition number = " << c << std::endl;

   if (nt > nm)
      _resinv = new TMatrixD(nm, nt);
   else
      _resinv = new TMatrixD(nt, nm);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5, 13, 4) /* TDecompSVD::Invert() didn't have ok status before 5.13/04. */
   Bool_t ok;
   *_resinv = _svd->Invert(ok);
   if (!ok) {
      std::cerr << "response matrix inversion failed" << std::endl;
      return false;
   }
#else
   *_resinv = _svd->Invert();
#endif
   if (nt > nm)
      _resinv->T();
   return true;
}

// Inline method definitions

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::RooUnfoldInvertT() : RooUnfoldT<Hist, Hist2D>()
{
   // Default constructor. Use Setup() to prepare for unfolding.
   Init();
}

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::RooUnfoldInvertT(const char *name, const char *title)
   : RooUnfoldT<Hist, Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   Init();
}

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D>::RooUnfoldInvertT(const TString &name, const TString &title)
   : RooUnfoldT<Hist, Hist2D>(name, title)
{
   // Basic named constructor. Use Setup() to prepare for unfolding.
   Init();
}

template <class Hist, class Hist2D>
RooUnfoldInvertT<Hist, Hist2D> &RooUnfoldInvertT<Hist, Hist2D>::operator=(const RooUnfoldInvertT<Hist, Hist2D> &rhs)
{
   // Assignment operator for copying RooUnfoldInvertTsettings.
   this->Assign(rhs);
   return *this;
}

template class RooUnfoldInvertT<TH1, TH2>;
ClassImp(RooUnfoldInvert)

#ifndef NOROOFIT
   template class RooUnfoldInvertT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist>;
typedef RooUnfoldInvertT<RooUnfolding::RooFitHist, RooUnfolding::RooFitHist> RooFitUnfoldInvert;
ClassImp(RooFitUnfoldInvert)
#endif
