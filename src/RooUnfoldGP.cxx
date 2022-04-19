/*! \class RooUnfoldGPT
  Paper: Unfolding with Gaussian Processes by A. Bozson, G. Cowan and F. Spano
*/

#include "RooUnfoldGP.h"

#include <iostream>
#include <math.h>  
#include <stdlib.h>

#include "TH1.h"
#include "TH2.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldFitHelpers.h"

using namespace RooUnfolding;

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::RooUnfoldGPT()
  : RooUnfoldT<Hist,Hist2D>()
{
  //! Default constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::RooUnfoldGPT(const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::RooUnfoldGPT(const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D> RooUnfolding::Algorithm
RooUnfoldGPT<Hist,Hist2D>::GetAlgorithm () const
{
  //! return the unfolding algorithm used
  return kGP;
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>& RooUnfoldGPT<Hist,Hist2D>::operator= (const RooUnfoldGPT<Hist,Hist2D>& rhs)
{
  //! Assignment operator for copying RooUnfoldInvertTsettings.
  this->Assign(rhs);
  return *this;
}


template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::RooUnfoldGPT(const RooUnfoldGPT<Hist,Hist2D>& rhs)
  : RooUnfoldT<Hist,Hist2D> (rhs)
{
  //! Copy constructor.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::RooUnfoldGPT(const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas,
						    Int_t kernel,						    const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _kernel(kernel)
{
  //! Constructor with response matrix object and measured unfolding input histogram.
  Init();
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::~RooUnfoldGPT()
{
  //! destructor
}

template<class Hist, class Hist2D>
RooUnfoldGPT<Hist,Hist2D>::Cache::~Cache()
{
  //! destructor
  delete _svd;
  delete _resinv;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::Init()
{
  _specialcache._svd= 0;
  _specialcache._resinv= 0;
  _specialcache._haveMLEst = false;
  _specialcache._haveMLCov = false;
  _specialcache._MLHConverged = false;

  this->GetSettings();
}


template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::Reset()
{
  delete _specialcache._svd;
  delete _specialcache._resinv;
  Init();
  RooUnfoldT<Hist,Hist2D>::Reset();
}

template<class Hist, class Hist2D> TDecompSVD*
RooUnfoldGPT<Hist,Hist2D>::Impl()
{
  return _specialcache._svd;
}

template<class Hist, class Hist2D> Double_t
RooUnfoldGPT<Hist,Hist2D>::GetRegParm() const
{
  return (Double_t)_kernel;
}
template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::SetRegParm(Double_t regparm)
{
  _kernel = regparm;
}

template<class Hist, class Hist2D> Bool_t
RooUnfoldGPT<Hist,Hist2D>::checkGP(const TVectorD& hist) const
{

  for (int i = 0; i < hist.GetNrows(); i++){
    if (hist[i] < 12){
      return false;
    }
  }
  
  return true;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::SetBinCenters() const
{

  Double_t delta_truBins = 0.0;
  Double_t delta_obsBins = 0.0;

  if (_kernel == 1){
    delta_truBins = 1 / (Double_t)(this->_nt - 2);
    delta_obsBins = 1 / (Double_t)(this->_nm - 2);
  }
  if (_kernel == 2){
    Double_t trumin = ::min(this->_res->Htruth(),RooUnfolding::X);
    Double_t trumax = ::max(this->_res->Htruth(),RooUnfolding::X);
    Double_t obsmin = ::min(this->_res->Hmeasured(),RooUnfolding::X);
    Double_t obsmax = ::max(this->_res->Hmeasured(),RooUnfolding::X);

    delta_truBins = fabs(trumax - trumin) / (Double_t)(this->_nt);
    delta_obsBins = fabs(obsmax - obsmin) / (Double_t)(this->_nm);
  }
   
  Double_t truBin_0 = - delta_truBins/2.0;
  Double_t obsBin_0 = - delta_obsBins/2.0;

  _specialcache._truBinCenters.ResizeTo(this->_nt);
  _specialcache._obsBinCenters.ResizeTo(this->_nm);
  
  for (int i = 0; i < _specialcache._truBinCenters.GetNrows(); i++){
    _specialcache._truBinCenters[i] = truBin_0 + i * delta_truBins;
  }
  for (int i = 0; i < _specialcache._obsBinCenters.GetNrows(); i++){
    _specialcache._obsBinCenters[i] = obsBin_0 + i * delta_obsBins;
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::SetFitSettings() const
{
  
  if (_kernel == 1){
    if (_specialcache._kernel_init.size() == 2) return;

    _specialcache._kernel_init.push_back(13.0);
    _specialcache._kernel_init.push_back(0.01);
    _specialcache._kernel_step.push_back(0.001);
    _specialcache._kernel_step.push_back(0.0001);
  } else {
    if (_specialcache._kernel_init.size() == 3) return;
    _specialcache._kernel_init.push_back(10.0);
    _specialcache._kernel_init.push_back(1.0);
    _specialcache._kernel_init.push_back(1.0);
    _specialcache._kernel_step.push_back(0.1);
    _specialcache._kernel_step.push_back(0.1);
    _specialcache._kernel_step.push_back(0.1);
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::Unfold() const
{
  if (_kernel == 2){
    std::cout << "ERROR: The Gibbs kernel is not fully implemented yet. Please use the radial kernel for now." << std::endl;
    return;
  }

  SetBinCenters();
  
  // Get the maximum likelihood solution first through
  // matrix inversion. This is the same as the RooUnfoldInversion
  // class.
  MLEstimator();
  MLCovariance();
 
  SetFitSettings();

  // Minimize the marginalized likehihood to get the optimal
  // kernel parameters.
  MinimizeMLH();
  
  if (!_specialcache._MLHConverged){
    return;
  }

  // Evaluate the kernel matrices with the optimal parameters.
  EvaluateK(_specialcache._opt_params);
  // EvaluateKstar(_opt_params);
  // EvaluateKstarstar(_opt_params);  
  
  // Calculate the maximum a postiori estimator for the 
  // truth histogram and the corresponding covariance matrix.
  MAPEstimator();

  // Pass the solutions.
  this->_cache._rec.ResizeTo(this->_nt);
  this->_cache._rec = _specialcache._MAPEst;

  this->_cache._unfolded= true;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::GetCov() const
{
  MAPCovariance();
  this->_cache._cov.ResizeTo(this->_nt,this->_nt);
  this->_cache._cov = _specialcache._MAPCov;
  this->_cache._haveCov= true;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::MLEstimator() const
{

  TMatrixD res(this->_res->Mresponse(true));
  if (this->_nt>this->_nm) {
    TMatrixD resT (TMatrixD::kTransposed, res);
    _specialcache._svd= new TDecompSVD (resT);
    delete _specialcache._resinv; _specialcache._resinv= 0;
  } else
    _specialcache._svd= new TDecompSVD (res);

  double c = _specialcache._svd->Condition();
  if (c<0 && (this->_verbose>0)) std::cout << "WARNING: Response matrix is ill-conditioned. TDecompSVD condition number = " << c << std::endl;

  _specialcache._MLEst.ResizeTo(this->_nm);
  _specialcache._MLEst= this->Vmeasured();


  if (this->_res->HasFakes()) {
    TVectorD fakes= this->_res->Vfakes();
    Double_t fac= this->_res->Vmeasured().Sum();
    if (fac!=0.0) fac=  this->Vmeasured().Sum() / fac;
    if (this->_verbose>=1) std::cout << "Subtract " << fac*fakes.Sum() << " fakes from measured distribution" << std::endl;
    fakes *= fac;
    _specialcache._MLEst -= fakes;
  }

  if (!checkGP(_specialcache._MLEst) && this->_verbose) {
    std::cout << "WARNING! Some of the bin counts are very low. The Gaussian bin count assumption might not hold resulting in bad unfolding results" << std::endl;
  }

  Bool_t ok;
  if (this->_nt>this->_nm) {
    ok= InvertResponse();
    if (ok) _specialcache._MLEst *= *_specialcache._resinv;
  } else
    ok= _specialcache._svd->Solve (_specialcache._MLEst);

  _specialcache._MLEst.ResizeTo(this->_nt);
  if (!ok) {
    std::cerr << "Response matrix Solve failed" << std::endl;
    return;
  }
  _specialcache._haveMLEst = true;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::MLCovariance() const
{
  TMatrixD cov(this->_nm, this->_nm);  
  TVectorD nu = this->_res->Vmeasured();

  for (int i = 0; i < this->_nm; i++){
    cov(i,i) = nu(i);
  }

    if (!InvertResponse()) return;
    _specialcache._MLCov.ResizeTo(this->_nt,this->_nt);
    ABAT (*_specialcache._resinv, cov, _specialcache._MLCov);
    
    _specialcache._haveMLCov = true;
}


template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::MAPEstimator() const
{

  TMatrixD KUML_INV(_specialcache._K + _specialcache._MLCov);
  TMatrixD KT(_specialcache._K);
  KUML_INV.Invert();
  KT.T();

  _specialcache._MAPEst.ResizeTo(this->_nt);

  _specialcache._MAPEst = MultiplyMatrixVector( (KT * KUML_INV) , _specialcache._MLEst);
}


template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::MAPCovariance() const
{

  TMatrixD KUML_INV(_specialcache._K + _specialcache._MLCov);
  TMatrixD KT(_specialcache._K);
  KUML_INV.Invert();
  KT.T();


  _specialcache._MAPCov.ResizeTo(this->_nt,this->_nt);

  _specialcache._MAPCov = _specialcache._K - (KT * KUML_INV) * _specialcache._K;  
}

template<class Hist, class Hist2D> Bool_t
RooUnfoldGPT<Hist,Hist2D>::InvertResponse() const
{
    if (!_specialcache._svd)   return false;
    if (_specialcache._resinv) return true;
    if (this->_nt>this->_nm) _specialcache._resinv= new TMatrixD(this->_nm,this->_nt);
    else         _specialcache._resinv= new TMatrixD(this->_nt,this->_nm);
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,13,4)  /* TDecompSVD::Invert() didn't have ok status before 5.13/04. */
    Bool_t ok;
    *_specialcache._resinv=_specialcache._svd->Invert(ok);
    if (!ok) {
      std::cerr << "response matrix inversion failed" << std::endl;
      return false;
    }
#else
    *_specialcache._resinv=_specialcache._svd->Invert();
#endif
    if (this->_nt>this->_nm) _specialcache._resinv->T();
    return true;
}


template<class Hist, class Hist2D> TVectorD
RooUnfoldGPT<Hist,Hist2D>::MultiplyMatrixVector(const TMatrixD& matrix, const TVectorD& vector) const
{

  TVectorD prod(matrix.GetNrows());

  if (matrix.GetNcols() != vector.GetNrows()){
    std::cerr << "Vector and matrix have incompatible dimensions for multiplication.";
    return prod;
  }

  for (int i = 0; i < matrix.GetNrows(); i++){
    for (int j = 0; j < matrix.GetNcols(); j++){
      prod[i] += matrix[i][j] * vector[j];
    }
  }

  return prod;
}

template<class Hist, class Hist2D> TVectorD
RooUnfoldGPT<Hist,Hist2D>::MultiplyVectorMatrix(const TVectorD& vector, const TMatrixD& matrix) const
{
  
  TVectorD prod(matrix.GetNcols());

  if (matrix.GetNrows() != vector.GetNrows()){
    std::cerr << "Vector and matrix have incompatible dimensions for multiplication.";
    return prod;
  }

  for (int i = 0; i < matrix.GetNcols(); i++){
    for (int j = 0; j < matrix.GetNrows(); j++){
      prod[i] += matrix[j][i] * vector[j];
    }
  }

  return prod;
  
}


template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::printMatrix(const TMatrixD& matrix) const
{
  
  for (int i = 0; i < matrix.GetNrows(); i++){
    for (int j = 0; j < matrix.GetNcols(); j++){

      //if (i == j){
      std::cout << "ixj: " << i << "x" << j << " " << sqrt(matrix[i][j]) << std::endl;
	//}
    }
  }
  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::printVector(const TVectorD& vector) const
{
  
  for (int i = 0; i < vector.GetNrows(); i++){
    std::cout << "i: " << i << " " << vector[i] << std::endl;
  }
  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::SetMinInit(std::vector<Double_t>& init) const
{
  switch (_kernel) {

  case 1: { 
    if (init.size() != 2){
      std::cerr << init.size() << " parameter initial values have been passed for the Marginal LH minimization. This should be 2 for the radial kernel.";
      return;
    } 
  }

  case 2: { 
    if (init.size() != 3){
      std::cerr << init.size() << " parameter initial values have been passed for the Marginal LH minimization. This should be 3 for the Gibbs kernel.";
      return;
    }
  }
    
  }
  _specialcache._kernel_init = init;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::SetMinStep(std::vector<Double_t>& step) const
{
  switch (_kernel) {

  case 1: { 
    if (step.size() != 2){
      std::cerr << step.size() << " parameter step sizes have been passed. This should be 2 for the radial kernel.";
      return;
    } 
  }

  case 2: { 
    if (step.size() != 3){
      std::cerr << step.size() << " parameter step sizes have been passed. This should be 3 for the Gibbs kernel.";
      return;
    }
  }
  }
  _specialcache._kernel_step = step;
}


template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::MinimizeMLH() const
{

  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  min->SetMaxFunctionCalls(1000000);
  min->SetMaxIterations(100000);
  min->SetTolerance(0.01);
  
  //min->SetPrintLevel(1);
  
  ROOT::Math::Functor f(this, &RooUnfoldGPT<Hist,Hist2D>::MarginalLH,_specialcache._kernel_init.size()); 

  min->SetFunction(f);
 
  for (size_t param = 0; param < _specialcache._kernel_init.size(); param++){
    TString name;
    name.Form("%zu",param);
    min->SetVariable(param,name.Data(),_specialcache._kernel_init.at(param), _specialcache._kernel_step.at(param));
  }

  min->Minimize();

  //if (min->Status()){
  //std::cout << "Minimization of the marginal likelihood did not converge. Try again with new initial values and step sizes." << std::endl;
  //} else {
  const double* xs = min->X();
  
  if (this->_verbose >= 2) std::cout << "Marginal likelihood minimization converged." << std::endl;
    
  for (size_t i = 0; i < _specialcache._kernel_init.size(); i++){
    if (this->_verbose >= 2) std::cout << "Parameter " << i << ": " << xs[i] << std::endl;
    _specialcache._opt_params.push_back(xs[i]);
  }
  
  _specialcache._MLHConverged = true;
  //}

  delete min;
}


template<class Hist, class Hist2D> Double_t
RooUnfoldGPT<Hist,Hist2D>::MarginalLH(const double *params) const
{
  
  Double_t first_term;
  Double_t second_term;

  if (this->_nt > this->_nm){
    std::cerr << "To many truth bins w.r.t. observed bins.";
  }

  EvaluateK(params);

  TMatrixD KUML_INV(_specialcache._K + _specialcache._MLCov);
  Double_t KUML_DET = KUML_INV.Determinant();
  KUML_INV.Invert();

  first_term = 0.5 * _specialcache._MLEst * MultiplyMatrixVector(KUML_INV, _specialcache._MLEst);
  
  second_term = 0.5 * log(KUML_DET);

  return first_term + second_term;
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateK(const double *params) const 
{

  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
     break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateKstar(const double *params) const 
{
  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._truBinCenters, _specialcache._obsBinCenters, _specialcache._Kstar);
    break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateKstarstar(const double *params) const 
{
  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._obsBinCenters, _specialcache._obsBinCenters, _specialcache._Kstarstar);
    break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }
}
template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateK(const std::vector<Double_t>& params) const 
{

  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateKstar(const std::vector<Double_t>& params) const 
{
  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._truBinCenters, _specialcache._obsBinCenters, _specialcache._Kstar);
    break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }
  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::EvaluateKstarstar(const std::vector<Double_t>& params) const
{
  switch (_kernel) {
  case 1: RadialKernel(params, _specialcache._obsBinCenters, _specialcache._obsBinCenters, _specialcache._Kstarstar);
    break;
  case 2: GibbsKernel(params, _specialcache._truBinCenters, _specialcache._truBinCenters, _specialcache._K);
    break;
  default: std::cout << "Unknown kernel not available." << std::endl;
  }
  
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::RadialKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const
{
  
  matrix.ResizeTo(x.GetNrows(), xprime.GetNrows());

  Double_t A = params[0];
  Double_t l = params[1];

  for (int i = 0; i < x.GetNrows(); i++){
    for (int j = 0; j < xprime.GetNrows(); j++){
      
      Double_t sqdist = pow(fabs(x[i] - xprime[j]),2);

      matrix[i][j] = exp(A) * exp(- 0.5 * sqdist / pow(l,2));
    }
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::RadialKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const
{
  
  matrix.ResizeTo(x.GetNrows(), xprime.GetNrows());

  Double_t A = params.at(0);
  Double_t l = params.at(1);

  for (int i = 0; i < x.GetNrows(); i++){
    for (int j = 0; j < xprime.GetNrows(); j++){
      
      Double_t sqdist = pow(fabs(x[i] - xprime[j]),2);

      matrix[i][j] = exp(A) * exp(- 0.5 * sqdist / pow(l,2));
    }
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::GibbsKernel(const double *params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const
{

  matrix.ResizeTo(x.GetNrows(),xprime.GetNrows());

  Double_t A = params[0];
  Double_t b = params[1];
  Double_t c = params[2];
  
  for (int i = 0; i < x.GetNrows(); i++){
    for (int j = 0; j < xprime.GetNrows(); j++){
      
      Double_t lx = b * x[i] + c;
      Double_t lxprime = b * xprime[j] + c;

      Double_t sqdist = (x[i] - xprime[j]) * (x[i] - xprime[j]);

      Double_t sqrt_term = (2*lx*lxprime)/((lx*lx) + (lxprime*lxprime));
      Double_t exp_term = - sqdist / ((lx*lx) + (lxprime*lxprime));

      matrix[i][j] = exp(A) * sqrt(sqrt_term) * exp(exp_term);
    }
  }
}

template<class Hist, class Hist2D> void
RooUnfoldGPT<Hist,Hist2D>::GibbsKernel(const std::vector<Double_t>& params, const TVectorD& x, const TVectorD& xprime, TMatrixD& matrix) const
{

  matrix.ResizeTo(x.GetNrows(), xprime.GetNrows());
  
  Double_t A = params.at(0);
  Double_t b = params.at(1);
  Double_t c = params.at(2);

  for (int i = 0; i < x.GetNrows(); i++){
    for (int j = 0; j < xprime.GetNrows(); j++){
      
      Double_t lx = b * x[i] + c;
      Double_t lxprime = b * xprime[i] + c;

      Double_t sqdist = pow(fabs(x[i] - xprime[j]),2);

      Double_t sqrt_term = (2*lx*lxprime)/(pow(lx,2) + pow(lxprime,2));
      Double_t exp_term = - sqdist / (pow(lx,2) + pow(lxprime,2));

      matrix[i][j] = exp(A) * sqrt(sqrt_term) * exp(exp_term);
    }
  }
}



template class RooUnfoldGPT<TH1,TH2>;
ClassImp (RooUnfoldGP)

#ifndef NOROOFIT
template class RooUnfoldGPT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
typedef RooUnfoldGPT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldGP;
ClassImp (RooFitUnfoldGP)
#endif



