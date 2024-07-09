//=====================================================================-*-C++-*-
#ifndef ROOFitUNFOLD_HH
#define ROOFitUNFOLD_HH
#include "RooUnfold.h"

#ifndef NOROOFIT
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooHistFunc.h>

class RooProdPdf;
class THStack;

//! \class RooUnfoldSpec
//! \brief Specifications class to build models based on RooUnfold in a RooFit context
//! \author Carsten Burgard <cburgard@cern.ch>
class RooUnfoldSpec : public TNamed {
public:
  enum Contribution {
                     kBackground,
                     kData,
                     kResponse,
                     kTruth,
                     kMeasured
  };
  typedef std::vector<RooAbsReal*> ShapeSys;
  typedef std::vector<double> NormSys;

protected:
  bool _locked = false;
  void lockCheck();

  class HistContainer {
    friend RooUnfoldSpec;
    RooAbsReal* _nom = 0;
    RooAbsReal* _staterror = 0;
    RooArgSet _obs;
    std::vector<RooRealVar*> _gammas;    
    std::map<const std::string,ShapeSys> _shapes;
    std::map<const std::string,NormSys> _norms;
    ~HistContainer();
    void setNominal(RooAbsReal* nom, const RooArgList& obslist);
    void setNominal(const TH1* nom, const RooArgList& obslist, double errorThreshold = -1, bool includeUnderflowOverflow = false, bool useDensity = false);
    void setNominal(RooDataHist* data, const RooArgList& obslist);
    void addShape(const char* name, RooAbsReal* up, RooAbsReal* dn);
    void addNorm(const char* name, double up, double dn);
  };
  bool _includeUnderflowOverflow = false;
  bool _useDensity = false;
  double _errorThreshold = -1;
  RooArgList _obs_truth;
  RooArgList _obs_reco;    
  RooArgList _obs_all;
  RooArgList _alphas;
  RooArgList _gammas;

  HistContainer _bkg;  
  HistContainer _data;
  HistContainer _res;
  HistContainer _truth;
  HistContainer _reco;    

  class Cache {
    friend RooUnfoldSpec;
    RooUnfolding::RooFitHist* _bkg = 0;
    RooUnfolding::RooFitHist* _data = 0;
    RooUnfolding::RooFitHist* _res = 0;
    RooUnfolding::RooFitHist* _truth = 0;
    RooUnfolding::RooFitHist* _reco = 0;
    RooUnfolding::RooFitHist* _data_minus_bkg = 0;
    RooFitUnfoldResponse* _response = 0;
  };

  void makeBackground();
  void makeData();
  void makeResponse();
  void makeTruth();
  void makeReco();
  void makeDataMinusBackground();

  void checkConsistency(const HistContainer& cont, const TH1* hist);  

  Cache _cache;

  RooUnfolding::RooFitHist* makeHistogram(const HistContainer& source, double errorThreshold);

public:

  RooProdPdf* makeConstraints();

  RooAbsReal* getBackground();
  RooAbsReal* getData();
  RooAbsReal* getResponse();
  RooAbsReal* getTruth();
  RooAbsReal* getReco();
  RooAbsReal* getDataMinusBackground();

  void addGaussNP(RooRealVar* v);
  void addPoissonNP(RooRealVar* v);

  RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);  
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* bkg, const TH1* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg_th1, const TH1* data_th1, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, RooAbsReal* reco, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const RooArgList& reco_bins, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* measured, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& measured_bins, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);  
  RooUnfoldSpec(const char* name, const char* title, RooAbsReal* truth, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);  
  RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);  

  ~RooUnfoldSpec();
  RooHistFunc* makeHistFuncT(const TH1* hist);
  RooHistFunc* makeHistFuncM(const TH1* hist);
  void registerSystematic(Contribution c, const char* name, const TH1* up, const TH1* down);
  void registerSystematic(Contribution c, const char* name, double up, double dn);
  RooAbsPdf* makePdf(RooUnfolding::Algorithm alg, Double_t regparam=-1e30);
  RooAbsReal* makeFunc(RooUnfolding::Algorithm alg, Double_t regparam=-1e30);
  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfold(RooUnfolding::Algorithm alg, Double_t regparam = -1e30);
  RooUnfolding::RooFitHist* makeHistogram(const TH1* hist);
  RooHistFunc* makeHistFuncTruth(const TH1* hist);
  RooHistFunc* makeHistFuncMeasured(const TH1* hist);
  TMatrixD makeCovarianceMatrix() const;
  THStack* makeBreakdownHistogram() const;

  const RooArgList& getTruthObservables() const { return _obs_truth; }
  const RooArgList& getRecoObservables() const { return _obs_reco; }
  const RooArgList& getObservables() const { return _obs_all; }  
  
protected:
  void addToCovarianceMatrix(const HistContainer& histContainer, TMatrixD& covarianceMatrix) const;
  void addShapeToCovarianceMatrix(const HistContainer& cont, const ShapeSys& var, TMatrixD& covarianceMatrix) const;
  void addStatToCovarianceMatrix(const HistContainer& histContainer, TMatrixD& covarianceMatrix) const;
  void addNormToCovarianceMatrix(const HistContainer& cont, const NormSys& var, TMatrixD& covarianceMatrix) const;
  void setup(const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg_th1, const TH1* data_th1, bool includeUnderflowOverflow, double errorThreshold = -1, bool useDensity = false);
  ClassDef(RooUnfoldSpec,0)
};

//! \class RooUnfoldFunc
//! \brief RooAbsReal Wrapper object for RooUnfold
//! \author Carsten Burgard <cburgard@cern.ch>
class RooUnfoldFunc : public RooAbsReal {
protected:
  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* _unfolding;
  bool _useDensity = false;
  mutable const RooArgSet* _curNormSet ; //! 
  
public:
  
  const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfolding() const ;
  
  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const override;
  virtual Double_t evaluate() const override;
  virtual TObject* clone(const char* newname = 0) const override;
  virtual Double_t getValV(const RooArgSet* set=0) const override;
  
  virtual Bool_t checkObservables(const RooArgSet *nset) const override;
  virtual Bool_t forceAnalyticalInt(const RooAbsArg &arg) const override;
  virtual Int_t getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &numVars, const RooArgSet *normSet, const char *rangeName = 0) const override;
  virtual Double_t analyticalIntegralWN(Int_t code, const RooArgSet *normSet, const char *rangeName = 0) const override;
  virtual void printMetaArgs(std::ostream &os) const override;
  virtual RooAbsArg::CacheMode canNodeBeCached() const override;
  virtual void setCacheAndTrackHints(RooArgSet &) override;
  
  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) override;

  RooArgList* makeParameterList() const;
  bool isDensity() const;
  void setDensity(bool d);
  
  RooUnfoldFunc();    
  RooUnfoldFunc(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf,bool useDensity=false);
  RooUnfoldFunc(const RooUnfoldFunc& other);
  RooUnfoldFunc(const RooUnfoldFunc* other);    
  virtual ~RooUnfoldFunc();
  ClassDefOverride(RooUnfoldFunc,1)
};

#endif
#endif
