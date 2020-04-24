// Author: Kerstin Tackmann, Andreas Hoecker, Heiko Lacker

/**********************************************************************************
 *                                                                                *
 * Project: TSVDUnfold - data unfolding based on Singular Value Decomposition     *
 * Package: ROOT                                                                  *
 * Class  : TSVDUnfold                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Single class implementation of SVD data unfolding based on:               *
 *          A. Hoecker, V. Kartvelishvili,                                        *
 *          "SVD approach to data unfolding"                                      *
 *          NIM A372, 469 (1996) [hep-ph/9509307]                                 *
 *                                                                                *
 * Authors:                                                                       *
 *      Kerstin Tackmann <Kerstin.Tackmann@cern.ch>   - CERN, Switzerland         *
 *      Andreas Hoecker  <Andreas.Hoecker@cern.ch>    - CERN, Switzerland         *
 *      Heiko Lacker     <lacker@physik.hu-berlin.de> - Humboldt U, Germany       *
 *                                                                                *
 * Copyright (c) 2010:                                                            *
 *      CERN, Switzerland                                                         *
 *      Humboldt University, Germany                                              *
 *                                                                                *
 **********************************************************************************/

//_______________________________________________________________________
/* Begin_Html
<center><h2>SVD Approach to Data Unfolding</h2></center>
<p>
Reference: <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a>
<p>
TSVDUnfold implements the singular value decomposition based unfolding method (see reference). Currently, the unfolding of one-dimensional histograms is supported, with the same number of bins for the measured and the unfolded spectrum.
<p>
The unfolding procedure is based on singular value decomposition of the response matrix. The regularisation of the unfolding is implemented via a discrete minimum-curvature condition.
<p>
Monte Carlo inputs:
<ul>
<li><tt>xini</tt>: true underlying spectrum (TH1D, n bins)
<li><tt>bini</tt>: reconstructed spectrum (TH1D, n bins)
<li><tt>Adet</tt>: response matrix (TH2, nxn bins)
</ul>
Consider the unfolding of a measured spectrum <tt>bdat</tt> with covariance matrix <tt>Bcov</tt> (if not passed explicitly, a diagonal covariance will be built given the errors of <tt>bdat</tt>). The corresponding spectrum in the Monte Carlo is given by <tt>bini</tt>, with the true underlying spectrum given by <tt>xini</tt>. The detector response is described by <tt>Adet</tt>, with <tt>Adet</tt> filled with events (not probabilities) with the true observable on the y-axis and the reconstructed observable on the x-axis.
<p>
The measured distribution can be unfolded for any combination of resolution, efficiency and acceptance effects, provided an appropriate definition of <tt>xini</tt> and <tt>Adet</tt>.<br><br>
<p>
The unfolding can be performed by
<ul>
<pre>
TSVDUnfold *tsvdunf = new TSVDUnfold( bdat, Bcov, bini, xini, Adet );
TH1D* unfresult = tsvdunf->Unfold( kreg );
</pre>
</ul>
where <tt>kreg</tt> determines the regularisation of the unfolding. In general, overregularisation (too small <tt>kreg</tt>) will bias the unfolded spectrum towards the Monte Carlo input, while underregularisation (too large <tt>kreg</tt>) will lead to large fluctuations in the unfolded spectrum. The optimal regularisation can be determined following guidelines in <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a> using the distribution of the <tt>|d_i|<\tt> that can be obtained by <tt>tsvdunf->GetD()</tt> and/or using pseudo-experiments.
<p>
Covariance matrices on the measured spectrum (for either the total uncertainties or individual sources of uncertainties) can be propagated to covariance matrices using the <tt>GetUnfoldCovMatrix</tt> method, which uses pseudo experiments for the propagation. In addition, <tt>GetAdetCovMatrix</tt> allows for the propagation of the statistical uncertainties on the response matrix using pseudo experiments. The covariance matrix corresponding to <tt>Bcov</tt> is also computed as described in <a href="http://arXiv.org/abs/hep-ph/9509307">Nucl. Instrum. Meth. A372, 469 (1996) [hep-ph/9509307]</a> and can be obtained from <tt>tsvdunf->GetXtau()</tt> and its (regularisation independent) inverse from  <tt>tsvdunf->GetXinv()</tt>. The distribution of singular values can be retrieved using <tt>tsvdunf->GetSV()</tt>.
<p>
See also the tutorial for a toy example.
End_Html */
//_______________________________________________________________________


#include <iostream>


#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#ifndef NOROOFIT
#include "RooUnfoldFitHelpers.h"
#endif

#include "RooUnfoldSvd.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDecompSVD.h"
#include "TRandom3.h"
#include "TMath.h"

using namespace std;
using namespace RooUnfolding;

//_______________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::SVDUnfold( const Hist *bdat, const TMatrixD& Bcov, const Hist *bini, const Hist *xini, const TMatrixD& Mdet, const TMatrixD& MdetE)
  :
  fMdim       (nBins(bdat,X)),
  fTdim       (nBins(xini,X)),
  fDdim       (2),
  fNormalize  (kFALSE),
  fKReg       (-1),
  fDHist      (NULL),
  fBdat       (h2v(bdat)),
  fBcov       (Bcov), 
  fBini       (h2v(bini)),
  fXini       (h2v(xini)),
  fAdet       (Mdet),
  fAdetE      (MdetE),
  fToyMode    (kFALSE),
  fMatToyMode (kFALSE)
{

  if (fMdim != nBins(bini,X) || 
      fTdim != nBins(xini,X) ||
      fMdim != Bcov.GetNrows() ||
      fMdim != Bcov.GetNcols() ||
      fMdim != fAdet.GetNrows() ||
      fTdim != fAdet.GetNcols()) {
    TString msg = "All histograms must have equal dimension.\n";
    msg += Form( "  Found: dim(bdat)=%i\n",    nBins(bdat,X) );
    msg += Form( "  Found: dim(Bcov)=%i,%i\n", Bcov.GetNrows(), Bcov.GetNcols() );
    msg += Form( "  Found: dim(bini)=%i\n",    nBins(bini,X) );
    msg += Form( "  Found: dim(xini)=%i\n",    nBins(xini,X) );
    msg += Form( "  Found: dim(Adet)=%i,%i\n", fAdet.GetNrows(), fAdet.GetNcols() );
    msg += "Please start again!";

    throw std::runtime_error(msg.Data());
  }


  fSVHist.ResizeTo(fTdim);
  fXtau.ResizeTo(fTdim,fTdim);
  fXinv.ResizeTo(fTdim,fTdim);
  fBcov.ResizeTo(fMdim,fMdim);
   
  // Get the input histos
  fDdim = 2; // This is the derivative used to compute the curvature matrix
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::SVDUnfold( const TVectorD& bdat, const TMatrixD& bcov, const TVectorD& bini, const TVectorD& xini, const TMatrixD& Mdet, const TMatrixD& MdetE)
  :
  fMdim       (bdat.GetNrows()),
  fTdim       (xini.GetNrows()),
  fDdim       (2),
  fNormalize  (kFALSE),
  fKReg       (-1),
  fDHist      (NULL),
  fBdat       (bdat),
  fBcov       (bcov),
  fBini       (bini),
  fXini       (xini),
  fAdet       (Mdet),
  fAdetE      (MdetE),
  fToyMode    (kFALSE),
  fMatToyMode (kFALSE)
{


  if (fMdim != bini.GetNrows() || 
      fTdim != xini.GetNrows() ||
      fMdim != bcov.GetNrows() ||
      fMdim != bcov.GetNcols() ||
      fMdim != fAdet.GetNrows() ||
      fTdim != fAdet.GetNcols()) {
    TString msg = "All histograms must have equal dimension.\n";
    msg += Form( "  Found: dim(bdat)=%i\n",    bdat.GetNrows() );
    msg += Form( "  Found: dim(bcov)=%i,%i\n", bcov.GetNrows(), bcov.GetNcols() );
    msg += Form( "  Found: dim(bini)=%i\n",    bini.GetNrows() );
    msg += Form( "  Found: dim(xini)=%i\n",    xini.GetNrows() );
    msg += Form( "  Found: dim(Adet)=%i,%i\n", fAdet.GetNrows(), fAdet.GetNcols() );
    msg += "Please start again!";

    throw std::runtime_error(msg.Data());
  }


  fSVHist.ResizeTo(fTdim);
  fXtau.ResizeTo(fTdim,fTdim);
  fXinv.ResizeTo(fTdim,fTdim);
  fBcov.ResizeTo(fMdim,fMdim);
   
  // Get the input histos
  fDdim = 2; // This is the derivative used to compute the curvature matrix
}
  

//_______________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::SVDUnfold( const Hist *bdat, const TMatrixD& Bcov, const Hist *bini, const Hist *xini, const Hist2D *Adet ) :
  SVDUnfold(bdat,Bcov,bini,xini,h2m(Adet),h2me(Adet)) 
{
  if(!sumW2N(Adet)){
    fAdetE *= 0.;
  }
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::SVDUnfold( const SVDUnfold& other )
   :
     fMdim       (other.fMdim),
     fTdim       (other.fTdim),
     fDdim       (other.fDdim),
     fNormalize  (other.fNormalize),
     fKReg       (other.fKReg),
     fDHist      (other.fDHist),
     fSVHist     (other.fSVHist),
     fXtau       (other.fXtau),
     fXinv       (other.fXinv),
     fBdat       (other.fBdat),
     fBcov       (other.fBcov),
     fBini       (other.fBini),
     fXini       (other.fXini),
     fAdet       (other.fAdet),
     fToyhisto   (other.fToyhisto),
     fToyhistoE   (other.fToyhistoE),     
     fToymat     (other.fToymat),
     fToyMode    (other.fToyMode),
     fMatToyMode (other.fMatToyMode) 
{
   // Copy constructor
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::~SVDUnfold()
{
   if(fDHist){
      delete fDHist;
      fDHist = 0;
   }

}

namespace {
  void sanitizeNaN(TMatrixD& m,double eps){
    for(int i=0; i<m.GetNrows(); ++i){
      for(int j=0; j<m.GetNcols(); ++j){
        if(std::isnan(m(i,j))) m(i,j) = eps;
      }
    }
  }
  void sanitizeNaN(TVectorD& v,double eps){
    for(int j=0; j<v.GetNrows(); ++j){
      if(std::isnan(v(j))) v(j) = eps;
    }
  }
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TVectorD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::UnfoldV( Int_t kreg )
{
   // Perform the unfolding with regularisation parameter kreg
   fKReg = kreg;

   // Create vectors and matrices
   TVectorD vb(fMdim), vberr(fMdim);
   TMatrixD mB(fBcov), mA(fMdim, fTdim), mCurv(fTdim, fTdim), mC(fTdim, fTdim);
   Double_t eps = 1e-12;

   // Copy histogams entries into vector
   if (fToyMode) { vb = fToyhisto; vberr = fToyhistoE; }
   else          { 
     vb = fBdat; 
     for (int i = 0; i < fMdim; i++){
       vberr = fBcov[i][i];
     }
   }

   TVectorD vbini = fBini;
   TVectorD vxini = fXini;


   ::sanitizeNaN(mB,eps);
   ::sanitizeNaN(vberr,eps);

   if (fMatToyMode) mA = fToymat;
   else        mA=fAdet;   

   // Fill and invert the second derivative matrix
   FillCurvatureMatrix( mCurv, mC );

   // Inversion of mC by help of SVD
   TDecompSVD CSVD(mC);
   TMatrixD CUort(CSVD.GetU());
   TMatrixD CVort(CSVD.GetV());
   TVectorD CSV  (CSVD.GetSig());

   TMatrixD CSVM(fTdim, fTdim);
   for (Int_t i=0; i<fTdim; i++) CSVM(i,i) = 1/CSV(i);

   CUort.Transpose( CUort );
   TMatrixD mCinv((CVort*CSVM)*CUort);

   //Rescale using the data covariance matrix
   TDecompSVD BSVD( mB );   
   TMatrixD QT(BSVD.GetU());
   QT.Transpose(QT);
   TVectorD B2SV(BSVD.GetSig());
   TVectorD BSV(B2SV);

   for(int i=0; i<fMdim; i++){
     BSV(i) = TMath::Sqrt(B2SV(i));
   }

   TMatrixD mAtmp(fMdim,fTdim);
   TVectorD vbtmp(fMdim);
   mAtmp *= 0;
   vbtmp *= 0;
   for(int i=0; i<fMdim; i++){
     for(int j=0; j<fMdim; j++){
       if(BSV(i)){
         vbtmp(i) += QT(i,j)*vb(j)/BSV(i);
       }
     }

     for(int j=0; j<fTdim; j++){
       for(int m=0; m<fMdim; m++){
	 if(BSV(i)){
	   mAtmp(i,j) += QT(i,m)*mA(m,j)/BSV(i);
	 }
       }
     }
   }

   mA = mAtmp;
   vb = vbtmp;

   // Singular value decomposition and matrix operations
   TDecompSVD ASVD( mA*mCinv );
   TMatrixD Uort(ASVD.GetU());
   TMatrixD Vort(ASVD.GetV());
   TVectorD ASV (ASVD.GetSig());

   if (!fToyMode && !fMatToyMode) {
     fSVHist = ASV;
   }

   TMatrixD Vreg(mCinv*Vort);
   Uort.Transpose(Uort);
   TVectorD vd(Uort*vb);

   // if (!fToyMode && !fMatToyMode) {
   //   fDHist = createHist<Hist>(vd,name(fBdat),title(fBdat),vars(fBdat),false);
   // }

   // Damping coefficient
   Int_t k = GetKReg()-1; 

   // Damping factors
   TVectorD vdz(fMdim);
   TMatrixD Z(fTdim, fTdim);
   for (Int_t i=0; i<fTdim; i++) {
     Double_t sreg;
     if (ASV(i)<ASV(0)*eps) sreg = ASV(0)*eps;
     else                   sreg = ASV(i);
     vdz(i) = sreg/(sreg*sreg + ASV(k)*ASV(k));
     Z(i,i) = vdz(i)*vdz(i);
   }

   TVectorD vz(CompProd( vd, vdz ));
   
   TMatrixD VortT(Vort);
   VortT.Transpose(VortT);
   TMatrixD W(mCinv*Vort*Z*VortT*mCinv);

   TMatrixD Xtau(fTdim, fTdim);
   TMatrixD Xinv(fTdim, fTdim);
   Xtau *= 0;
   Xinv *= 0;
   for (Int_t i=0; i<fTdim; i++) {
     for (Int_t j=0; j<fTdim; j++) {
       Xtau(i,j) =  vxini(i) * vxini(j) * W(i,j);

       double a=0;
       for (Int_t m=0; m<fMdim; m++) {
         a += mA(m,i)*mA(m,j);
       }
       if(vxini(i) && vxini(j))
         Xinv(i,j) = a/vxini(i)/vxini(j);
     }
   }

   vz.ResizeTo(Vreg.GetNrows());
   
   // Compute the weights
   TVectorD vw(Vreg*vz);

   // Rescale by xini
   TVectorD vx(CompProd( vw, vxini ));
   
   if(fNormalize){ // Scale result to unit area
     Double_t scale = vx.Sum();
     if (scale > 0){
       vx *= 1.0/scale;
       Xtau *= 1./scale/scale;
       Xinv *= scale*scale;
     }
   }

   if (!fToyMode && !fMatToyMode) {
     fXtau = Xtau;
     fXinv = Xinv;
   }
   
   // Get Curvature and also chi2 in case of MC unfolding
//   if (!fToyMode && !fMatToyMode) {
//     std::cout << TString::Format( "Unfolding param: %i",k+1 ) << std::endl;
//     std::cout << TString::Format( "Curvature of weight distribution: %f", GetCurvature( vw, mCurv ) ) << std::endl;
//   }
   return vx;
}

//_______________________________________________________________________
// template<class Hist,class Hist2D>
// Hist* RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::Unfold( Int_t kreg )
// {

//   TVectorD vx(UnfoldV(kreg));
//   return createHist<Hist>(vx,"unfoldingresult",title(fBdat),vars(fBdat),false);
// }

//_______________________________________________________________________
template<class Hist,class Hist2D>
TMatrixD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetUnfoldCovMatrix( const TMatrixD& cov, Int_t ntoys, Int_t seed )
{
   // Determine for given input error matrix covariance matrix of unfolded 
   // spectrum from toy simulation given the passed covariance matrix on measured spectrum
   // "cov"    - covariance matrix on the measured spectrum, to be propagated
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // Note that this covariance matrix will contain effects of forced normalisation if spectrum is normalised to unit area. 

   fToyMode = true;
  
   // Code for generation of toys (taken from RooResult and modified)
   // Calculate the elements of the upper-triangular matrix L that
   // gives Lt*L = C, where Lt is the transpose of L (the "square-root method")  
   TMatrixD L(fMdim,fTdim); L *= 0;

   for (Int_t iPar= 0; iPar < fMdim; iPar++) {

      // Calculate the diagonal term first
     L(iPar,iPar) = cov(iPar,iPar);
     for (Int_t k=0; k<iPar; k++) L(iPar,iPar) -= TMath::Power( L(k,iPar), 2 );
     if (L(iPar,iPar) > 0.0) L(iPar,iPar) = TMath::Sqrt(L(iPar,iPar));
     else                    L(iPar,iPar) = 0.0;
     
     // ...then the off-diagonal terms
     for (Int_t jPar=iPar+1; jPar<fTdim; jPar++) {
       L(iPar,jPar) = cov(iPar,jPar);
       for (Int_t k=0; k<iPar; k++) L(iPar,jPar) -= L(k,iPar)*L(k,jPar);
       if (L(iPar,iPar)!=0.) L(iPar,jPar) /= L(iPar,iPar);
       else                  L(iPar,jPar) = 0;
     }
   }
   
   // Remember it
   TMatrixD Lt(TMatrixD::kTransposed,L);
   TRandom3 random(seed);

   fToyhisto = fBdat;
   fToyhistoE = fBdat;
   for (int i = 0; i < fToyhistoE.GetNrows(); i++){
     fToyhistoE[i] = fBcov[i][i];
   }

   TVectorD toymean(fMdim);

   // Get the mean of the toys first
   for (int i=1; i<=ntoys; i++) {

      // create a vector of unit Gaussian variables
      TVectorD g(fTdim);
      for (Int_t k= 0; k < fTdim; k++) g(k) = random.Gaus(0.,1.);

      // Multiply this vector by Lt to introduce the appropriate correlations
      g *= Lt;

      // Add the mean value offsets and store the results
      for (int j=0; j<fMdim; j++) {
	fToyhisto[j] = fBdat[j]+g(j-1);
        fToyhistoE[j] = fBcov[j][j];
      }

      TVectorD unfres(UnfoldV(GetKReg()));

      for (Int_t j=0; j<fMdim; j++) {
        toymean[j] = toymean[j] + unfres[j]/ntoys;
      }
   }

   // Reset the random seed
   random.SetSeed(seed);

   TMatrixD unfcov(fAdet.GetNrows(),fAdet.GetNcols());
   
   //Now the toys for the covariance matrix
   for (int i=1; i<=ntoys; i++) {

      // Create a vector of unit Gaussian variables
      TVectorD g(fMdim);
      for (Int_t k= 0; k < fMdim; k++) g(k) = random.Gaus(0.,1.);

      // Multiply this vector by Lt to introduce the appropriate correlations
      g *= Lt;

      // Add the mean value offsets and store the results
      for (int j=0; j<fMdim; j++) {
	fToyhisto[j] = fBdat[j]+g(j-1);
        fToyhistoE[j] = fBcov[j][j];
      }
      TVectorD unfres(UnfoldV(GetKReg()));
      
      for (Int_t j=0; j<fTdim; j++) {
        for (Int_t k=0; k<fTdim; k++) {
          unfcov(j,k) = unfcov(j,k) + ( (unfres[j] - toymean[j])* (unfres[k] - toymean[k])/(ntoys-1));
        }
      }
   }
   fToyMode = kFALSE;
   
   return unfcov;
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TMatrixD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetAdetCovMatrix( Int_t ntoys, Int_t seed, const TMatrixD* uncmat )
{
   // Determine covariance matrix of unfolded spectrum from finite statistics in 
   // response matrix using pseudo experiments
   // "ntoys"  - number of pseudo experiments used for the propagation
   // "seed"   - seed for pseudo experiments
   // "uncmat" - matrix to be interpreted as uncertainties on detector response matrix, to be propagated by toys, if no matrix passed, uncertainties on Adet will be used in toys with Gaussian smearing if Sumw2 is set for Adet, otherwise Poisson variations on Adet will be performed in toys
  if (uncmat && (uncmat->GetNrows() != fTdim || uncmat->GetNcols() != fTdim)) 
    {
      TString msg = "Uncertainty histogram must have the same dimension as all other histograms.\n";
      msg += Form( "  Found: dim(uncmat)=%i,%i\n", uncmat->GetNrows(), uncmat->GetNcols() );
      msg += Form( "  Found: dim(Adet)=%i,%i\n", fTdim, fTdim );
      msg += "Please start again!";
      throw std::runtime_error( msg.Data()) ;
    }
  
  fMatToyMode = true;
  
  //Now the toys for the detector response matrix
  TRandom3 random(seed);
  
  fToymat.ResizeTo(fAdet.GetNrows(),fAdet.GetNcols());
  TVectorD toymean(fXini.GetNrows());
    
  for (int i=0; i<ntoys; i++) {    
    for (Int_t k=0; k<fMdim; k++) {
      for (Int_t m=0; m<fTdim; m++) {
	if (fAdet(k,m)){
	  if(uncmat)
	    fToymat(k, m) =  fAdet(k,m)+random.Gaus(0.,(*uncmat)(k,m));
	  else if(fAdetE.GetNrows() == fAdet.GetNrows())
	    fToymat(k, m) = fAdet(k,m)+random.Gaus(0.,fAdetE(k,m));
	  else 
	    fToymat(k, m) = random.Poisson(fAdet(k,m));
	}
      }
    }

    TVectorD unfres(UnfoldV(GetKReg()));

    for (Int_t j=0; j<fTdim; j++) {
        toymean[j] = toymean[j] + unfres(j)/ntoys;
    }
  }

  // Reset the random seed
  random.SetSeed(seed);
  
  TMatrixD unfcov(fAdet.GetNrows(),fAdet.GetNcols());
  
  for (int i=1; i<=ntoys; i++) {
    for (Int_t k=0; k<fMdim; k++) {
      for (Int_t m=0; m<fTdim; m++) {
	if (fAdet(k,m)){
	  if(uncmat)
	    fToymat(k, m) = fAdet(k,m)+random.Gaus(0.,(*uncmat)(k,m));
	  else if(fAdetE.GetNrows() == fAdet.GetNrows())
	    fToymat(k, m) = fAdet(k,m)+random.Gaus(0.,fAdetE(k,m));
	  else 
	    fToymat(k, m) = random.Poisson(fAdet(k,m));
	}
      }
    }
    
    TVectorD unfres(UnfoldV(GetKReg()));
    
    for (Int_t j=0; j<fTdim; j++) {
      for (Int_t k=0; k<fTdim; k++) {
	unfcov(j,k) = unfcov(j,k) + ( (unfres[j] - toymean[j])*(unfres[k] - toymean[k])/(ntoys-1));
      }
    }
  }
  fMatToyMode = kFALSE;

  return unfcov;
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
Hist* RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetD() const 
{ 
  // Returns d vector (for choosing appropriate regularisation)
  return fDHist; 
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TVectorD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetSV() const 
{ 
   // Returns singular values vector
   return fSVHist; 
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
const TMatrixD& RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetXtau() const 
{ 
   // Returns the computed regularized covariance matrix corresponding to total uncertainties on measured spectrum as passed in the constructor.
  // Note that this covariance matrix will not contain the effects of forced normalization if spectrum is normalized to unit area.
  return fXtau; 
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
const TMatrixD& RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetXinv() const 
{ 
   // Returns the computed inverse of the covariance matrix
   return fXinv; 
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
const TMatrixD& RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetBCov() const 
{ 
   // Returns the covariance matrix
   return fBcov; 
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TVectorD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::VecDiv( const TVectorD& vec1, const TVectorD& vec2, Int_t zero )
{
   // Divide entries of two vectors
   TVectorD quot(vec1.GetNrows());
   for (Int_t i=0; i<vec1.GetNrows(); i++) {
      if (vec2(i) != 0) quot(i) = vec1(i) / vec2(i);
      else {
         if   (zero) quot(i) = 0;
         else        quot(i) = vec1(i);
      }
   }
   return quot;
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TMatrixD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::MatDivVec( const TMatrixD& mat, const TVectorD& vec, Int_t zero )
{
   // Divide matrix entries by vector
   TMatrixD quotmat(mat.GetNrows(), mat.GetNcols());
   for (Int_t i=0; i<mat.GetNrows(); i++) {
      for (Int_t j=0; j<mat.GetNcols(); j++) {
         if (vec(i) != 0) quotmat(i,j) = mat(i,j) / vec(i);
         else {
            if   (zero) quotmat(i,j) = 0;
            else        quotmat(i,j) = mat(i,j);
         }
      }
   }
   return quotmat;
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
TVectorD RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::CompProd( const TVectorD& vec1, const TVectorD& vec2 )
{
   // Multiply entries of two vectors
   TVectorD res(vec1.GetNrows());
   for (Int_t i=0; i<vec1.GetNrows(); i++) res(i) = vec1(i) * vec2(i);
   return res;
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
Double_t RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::GetCurvature(const TVectorD& vec, const TMatrixD& curv) 
{      
   // Compute curvature of vector
   return vec*(curv*vec);
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
void RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::FillCurvatureMatrix( TMatrixD& tCurv, TMatrixD& tC ) const
{
   Double_t eps = 0.00001;

   Int_t ndim = tCurv.GetNrows();

   // Init
   tCurv *= 0;
   tC    *= 0;

   if (fDdim == 0) for (Int_t i=0; i<ndim; i++) tC(i,i) = 1;
   else if (fDdim == 1) {
      for (Int_t i=0; i<ndim; i++) {
         if (i < ndim-1) tC(i,i+1) = 1.0;
         tC(i,i) = 1.0;
      }
   }
   else if (fDdim == 2) {
      for (Int_t i=0; i<ndim; i++) {
         if (i > 0)      tC(i,i-1) = 1.0;
         if (i < ndim-1) tC(i,i+1) = 1.0;
         tC(i,i) = -2.0;
      }
      tC(0,0) = -1.0;
      tC(ndim-1,ndim-1) = -1.0;
   }
   else if (fDdim == 3) {
      for (Int_t i=1; i<ndim-2; i++) {
         tC(i,i-1) =  1.0;
         tC(i,i)   = -3.0;
         tC(i,i+1) =  3.0;
         tC(i,i+2) = -1.0;
      }
   }
   else if (fDdim==4) {
      for (Int_t i=0; i<ndim; i++) {
         if (i > 0)      tC(i,i-1) = -4.0;
         if (i < ndim-1) tC(i,i+1) = -4.0;
         if (i > 1)      tC(i,i-2) =  1.0;
         if (i < ndim-2) tC(i,i+2) =  1.0;
         tC(i,i) = 6.0;
      }
      tC(0,0) = 2.0;
      tC(ndim-1,ndim-1) = 2.0;
      tC(0,1) = -3.0;
      tC(ndim-2,ndim-1) = -3.0;
      tC(1,0) = -3.0;
      tC(ndim-1,ndim-2) = -3.0;
      tC(1,1) =  6.0;
      tC(ndim-2,ndim-2) =  6.0;
   }
   else if (fDdim == 5) {
      for (Int_t i=2; i < ndim-3; i++) {
         tC(i,i-2) = 1.0;
         tC(i,i-1) = -5.0;
         tC(i,i)   = 10.0;
         tC(i,i+1) = -10.0;
         tC(i,i+2) = 5.0;
         tC(i,i+3) = -1.0;
      }
   }
   else if (fDdim == 6) {
      for (Int_t i = 3; i < ndim - 3; i++) {
         tC(i,i-3) = 1.0;
         tC(i,i-2) = -6.0;
         tC(i,i-1) = 15.0;
         tC(i,i)   = -20.0;
         tC(i,i+1) = 15.0;
         tC(i,i+2) = -6.0;
         tC(i,i+3) = 1.0;
      }
   }

   // Add epsilon to avoid singularities
   for (Int_t i=0; i<ndim; i++) tC(i,i) = tC(i,i) + eps;

   //Get curvature matrix
   for (Int_t i=0; i<ndim; i++) {
      for (Int_t j=0; j<ndim; j++) {
         for (Int_t k=0; k<ndim; k++) {
            tCurv(i,j) = tCurv(i,j) + tC(k,i)*tC(k,j);
         }
      }
   }
}

//_______________________________________________________________________
template<class Hist,class Hist2D>
void RooUnfoldSvdT<Hist,Hist2D>::SVDUnfold::RegularisedSymMatInvert( TMatrixDSym& mat, Double_t eps )
{
   // naive regularised inversion cuts off small elements

   // init reduced matrix
   const UInt_t n = mat.GetNrows();
   UInt_t nn = 0;   

   UInt_t *ipos = new UInt_t[n];
   //   UInt_t ipos[n];

   // find max diagonal entries
   Double_t ymax = 0;
   for (UInt_t i=0; i<n; i++) if (TMath::Abs(mat[i][i]) > ymax) ymax = TMath::Abs(mat[i][i]);

   for (UInt_t i=0; i<n; i++) {

         // save position of accepted entries
      if (TMath::Abs(mat[i][i])/ymax > eps) ipos[nn++] = i;
   }

   // effective matrix
   TMatrixDSym matwork( nn );
   for (UInt_t in=0; in<nn; in++) for (UInt_t jn=0; jn<nn; jn++) matwork(in,jn) = 0;

   // fill non-zero effective working matrix
   for (UInt_t in=0; in<nn; in++) {

      matwork[in][in] = mat[ipos[in]][ipos[in]];
      for (UInt_t jn=in+1; jn<nn; jn++) {
         matwork[in][jn] = mat[ipos[in]][ipos[jn]];
         matwork[jn][in] = matwork[in][jn];
      }
   }

   // invert
   matwork.Invert();

   // reinitialise old matrix
   for (UInt_t i=0; i<n; i++) for (UInt_t j=0; j<n; j++) mat[i][j] = 0;

   // refill inverted matrix in old one
   for (UInt_t in=0; in<nn; in++) {
      mat[ipos[in]][ipos[in]] = matwork[in][in];
      for (UInt_t jn=in+1; jn<nn; jn++) {
         mat[ipos[in]][ipos[jn]] = matwork[in][jn];
         mat[ipos[jn]][ipos[in]] = mat[ipos[in]][ipos[jn]];
      }
   }
   delete []  ipos;
}

template class RooUnfoldSvdT<TH1,TH2>::SVDUnfold;
ClassImp (RooUnfoldSvd::SVDUnfold)

#ifndef NOROOFIT
template class RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::SVDUnfold;
typedef RooUnfoldSvdT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldSvd;
ClassImp (RooFitUnfoldSvd::SVDUnfold)
#endif
