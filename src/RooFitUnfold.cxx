#ifndef NOROOFIT
#include <sstream>
#include "RooUnfold.h"
#include "RooFitUnfold.h"
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldFitHelpers.h"
#include "RooHistFunc.h"
#include "RooRealSumFunc.h"
#include "RooPrintable.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooPoisson.h"
#include "RooStats/HistFactory/PiecewiseInterpolation.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include "RooStats/HistFactory/ParamHistFunc.h"
#include "RooProduct.h"
#include "RooWrapperPdf.h"
#include "RooBinning.h"


using namespace RooUnfolding;

namespace { 
  void setBinning(RooRealVar* obs, const TAxis* ax, bool includeUnderflowOverflow){
    int n = ax->GetNbins()+(includeUnderflowOverflow?2:0);
    int offset = (includeUnderflowOverflow?0:1);
    if(ax->IsVariableBinSize()){
      std::vector<double> bounds;
      for(int i=0; i<n; ++i){
        bounds.push_back(ax->GetBinLowEdge(offset + i));
      }
      bounds.push_back(ax->GetBinUpEdge(n));
      RooBinning bins(n,&((bounds[0])));
      obs->setBinning(bins);
    } else {
      obs->setBins(n);
    }
  }
}

RooUnfoldFunc::RooUnfoldFunc(const char* name, const char* title, const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unf, bool useDensity) :
  RooAbsReal(name,title),
  _useDensity(useDensity)                                                                                                                                                               
{
  //! constructor
  auto unfolding = RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::New(unf->GetAlgorithm(),unf->response(),unf->Hmeasured(),unf->GetRegParm(),unf->GetName(),unf->GetTitle());
  if (unf->Hbkg()) unfolding->SetBkg(unf->Hbkg());
  unfolding->SetTruth(unf->Htruth());
  this->_unfolding = dynamic_cast<RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>*>(unfolding);
  this->_unfolding->SetVerbose(0);
  const RooUnfoldResponseT<RooFitHist,RooFitHist>* res = this->_unfolding->response();
  if(res){
    const RooFitHist* htruth = res->Htruth();
    if(htruth){
      this->addServer(*(htruth->func()));
      for(int i=0; i<dim(htruth); ++i){
        this->addServer(*htruth->obs(i));
      }
    }
    const RooFitHist* hfakes = res->Hfakes();
    if(hfakes){
      this->addServer(*(hfakes->func()));
      for(int i=0; i<dim(hfakes); ++i){
        this->addServer(*hfakes->obs(i));
      }
    }
    const RooFitHist* hresponse = res->Hresponse();
    if(hresponse){
      this->addServer(*(hresponse->func()));
      for(int i=0; i<dim(hresponse); ++i){
        this->addServer(*hresponse->obs(i));
      }
    }
    const RooFitHist* hmeasured = res->Hmeasured();
    if(hmeasured){
      this->addServer(*(hmeasured->func()));    
      for(int i=0; i<dim(hmeasured); ++i){
        this->addServer(*hmeasured->obs(i));
      }
    }
  }
  const RooFitHist* hmeasured = this->_unfolding->Hmeasured();

  if(hmeasured){
    this->addServer(*(hmeasured->func()));    
    for(int i=0; i<dim(hmeasured); ++i){
      this->addServer(*hmeasured->obs(i));
    }
  }
}
RooUnfoldFunc::RooUnfoldFunc() : _unfolding(NULL) {
  //! constructor
}
RooUnfoldFunc::~RooUnfoldFunc(){
  //! destructor
  delete _unfolding;
}


Bool_t RooUnfoldFunc::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive){
  //! redirect servers
  RooUnfoldResponseT<RooFitHist,RooFitHist>* res = this->_unfolding->response();
  if(res){
    RooFitHist* htruth = res->Htruth();
    if(htruth){
      htruth->replace(newServerList);
    }
    RooFitHist* hfakes = res->Hfakes();
    if(hfakes){
      hfakes->replace(newServerList);
    }
    RooFitHist* hresponse = res->Hresponse();
    if(hresponse){
      hresponse->replace(newServerList);
    }
    RooFitHist* hmeasured = res->Hmeasured();
    if(hmeasured){
      hmeasured->replace(newServerList);
    }
  }
  RooFitHist* hmeasured = this->_unfolding->Hmeasured();
  if(hmeasured){
    hmeasured->replace(newServerList);
  }
  return RooAbsReal::redirectServersHook(newServerList,mustReplaceAll,nameChange,isRecursive);
}
 




const RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* RooUnfoldFunc::unfolding() const { 
  //! retrieve the unfolding object
  return this->_unfolding;
}

std::list<Double_t>* RooUnfoldFunc::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  //! retrieve the list of bin boundaries
  return this->_unfolding->response()->Htruth()->func()->binBoundaries(obs,xlo,xhi);
}

std::list<Double_t>* RooUnfoldFunc::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  //! retrieve the sampling hint
  return this->_unfolding->response()->Htruth()->func()->plotSamplingHint(obs,xlo,xhi);
}

Double_t RooUnfoldFunc::getValV(const RooArgSet* set) const
{
  //! return the value
  this->_curNormSet = set ;
  return RooAbsReal::getValV(set) ;
}

RooArgList* RooUnfoldFunc::makeParameterList() const {
  //! return a list of all parameters in this function
  RooArgSet obs;
  for(int d=0; d<this->unfolding()->response()->Hresponse()->dim(); ++d){
    obs.add(*(this->unfolding()->response()->Hresponse()->obs(d)));
  }
  RooArgSet* pset = this->getParameters(&obs);

  RooArgList* list = new RooArgList(*pset);
  delete pset;
  return list;
}

bool RooUnfoldFunc::isDensity() const {
  //! return true if the return value is density-corrected, false otherwise
  return this->_useDensity;
}
void RooUnfoldFunc::setDensity(bool d){
  //! set if the return value should be density-corrected
  this->_useDensity = d;
}


Double_t RooUnfoldFunc::evaluate() const {
  //! call getVal on the internal function
  std::map<std::string,double> snapshot;
  this->_unfolding->response()->Hresponse()->saveSnapshot(snapshot);
  int bin = this->_unfolding->response()->Htruth()->bin();
  this->_unfolding->ForceRecalculation();
  this->_unfolding->response()->Htruth()->checkValidity();
  double v = this->_unfolding->Vunfold()[bin];
  if(this->_useDensity){
    v /= binVolume(this->_unfolding->response()->Htruth(),bin,false);
  }
  this->_unfolding->response()->Hresponse()->loadSnapshot(snapshot);
  return v;  
}

Bool_t  RooUnfoldFunc::isBinnedDistribution(const RooArgSet& obs) const {
  //! check if this PDF is a binned distribution in the given observable
  return this->_unfolding->response()->Hresponse()->func()->isBinnedDistribution(obs);
}


Bool_t RooUnfoldFunc::checkObservables(const RooArgSet *nset) const {
  //! call checkOvservables on the response
  return this->_unfolding->response()->Hresponse()->func()->checkObservables(nset);
}


Bool_t RooUnfoldFunc::forceAnalyticalInt(const RooAbsArg &arg) const {
  //! force the analytical integral
  return this->_unfolding->response()->Htruth()->func()->forceAnalyticalInt(arg);
}


Int_t RooUnfoldFunc::getAnalyticalIntegralWN(RooArgSet &allVars, RooArgSet &numVars, const RooArgSet *normSet, const char *rangeName) const {
  //! retrieve the analytical integral status
  return this->_unfolding->response()->Htruth()->func()->getAnalyticalIntegralWN(allVars,numVars,normSet,rangeName);
}


Double_t RooUnfoldFunc::analyticalIntegralWN(Int_t code, const RooArgSet *normSet, const char *rangeName) const {
  //! retrieve the analytical integral status
  double val = 0;
  auto vec = this->_unfolding->Vunfold();
  for(int i=0; i<vec.GetNrows(); ++i){
    // assuming that density correction has been applied already
    val += vec[i];
  }
  return val;
}


void RooUnfoldFunc::printMetaArgs(std::ostream &os) const {
  //! printing helper function
  return this->_unfolding->response()->Htruth()->func()->printMetaArgs(os);
}


RooAbsArg::CacheMode RooUnfoldFunc::canNodeBeCached() const {
  return this->_unfolding->response()->Htruth()->func()->canNodeBeCached();
}

void RooUnfoldFunc::setCacheAndTrackHints(RooArgSet& arg) {
  this->_unfolding->response()->Htruth()->func()->setCacheAndTrackHints(arg);
}
TObject* RooUnfoldFunc::clone(const char* newname) const {
  //! produce a clone (deep copy) of this object
  return new RooUnfoldFunc(newname ? newname : this->GetName(),this->GetTitle(),this->_unfolding,this->_useDensity);
}

namespace {
  bool readToken(TString& instr, std::vector<TString>& tokens){
    int pos = instr.First(",");
    if(pos == -1){
      tokens.push_back(instr);
      return false;
    } else {
      tokens.push_back(instr(0,pos));
      instr.Remove(0,pos+1);
      return true;
    }
  }
}

#include <RooAbsCollection.h>

#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
namespace {
  
  struct IteratorHelper {
    RooFIter itr;
    RooAbsArg* nextItem;
    IteratorHelper(const RooAbsCollection& c);
    IteratorHelper();
    RooAbsArg* operator++();
    bool operator!=(const IteratorHelper& other);
    bool operator!=(const RooAbsArg* other);  
    RooAbsArg* operator*();
  };
  
  IteratorHelper::IteratorHelper(const RooAbsCollection& c) : itr(c.fwdIterator()), nextItem(itr.next()) {}
  IteratorHelper::IteratorHelper() : itr(), nextItem(NULL) {}  
  RooAbsArg* IteratorHelper::operator++(){
    nextItem = itr.next();
    return nextItem;
  }
  bool IteratorHelper::operator!=(const IteratorHelper& other){
    return this->nextItem != other.nextItem;
  }
  bool IteratorHelper::operator!=(const RooAbsArg* other){
    return this->nextItem != other;
  }
  
  RooAbsArg* IteratorHelper::operator*(){
    return nextItem;
  }
}

::IteratorHelper begin(const RooAbsCollection& c){
  return ::IteratorHelper(c);
}

::IteratorHelper end(const RooAbsCollection&){
  return ::IteratorHelper();
}
#endif

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  RooUnfoldSpec(name,title,truth,obs_truth,reco,obs_reco,response,0,data,includeUnderflowOverflow,errorThreshold,useDensity)
{
  //! constructor forwarding
}


RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth, const char* obs_truth, const TH1* reco, const char* obs_reco, const TH2* response, const TH1* bkg, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  //! constructor
  int d =dim(truth);
  if(d!=dim(reco)){
    throw std::runtime_error("inconsistent dimensionality between truth and reco histograms!");
  }
  
  TString obs_truth_s(obs_truth);
  TString obs_reco_s(obs_reco);
  std::vector<TString> obs_truth_v;
  std::vector<TString> obs_reco_v;
  bool more_truth = false;
  bool more_reco = false;
  for(int i=0; i<d; ++i){
    more_truth = ::readToken(obs_truth_s,obs_truth_v);
    more_reco = ::readToken(obs_reco_s,obs_reco_v);
    if(!more_truth || !more_reco) break;
  }
  if(more_truth) throw std::runtime_error(TString::Format("encountered additional characters on truth observable list: '%s'",obs_truth_s.Data()).Data());
  if(more_reco) throw std::runtime_error(TString::Format("encountered additional characters on reco observable list: '%s'",obs_reco_s.Data()).Data());
  if((int)obs_truth_v.size() != d) throw std::runtime_error(TString::Format("truth observable list is too short for %d dimensions: '%s'",d,obs_truth).Data());
  if((int)obs_reco_v.size() != d) throw std::runtime_error(TString::Format("reco observable list is too short for %d dimensions: '%s'",d,obs_truth).Data());

  RooArgList truth_vars;
  for(int i=0; i<d; ++i){
    const TAxis* ax = RooUnfolding::getAxis(truth,(RooUnfolding::Dimension)i);
    double min = ax->GetBinLowEdge(!includeUnderflowOverflow);
    double max = ax->GetBinUpEdge(ax->GetNbins()+includeUnderflowOverflow);
    RooRealVar* obs = new RooRealVar(obs_truth_v[i],ax->GetTitle() ? ax->GetTitle() : obs_truth_v[i].Data(),min,min,max);
    setBinning(obs,ax,includeUnderflowOverflow);
    obs->setConstant(true);
    truth_vars.add(*obs);
  }

  RooArgList reco_vars;
  for(int i=0; i<d; ++i){
    const TAxis* ax = RooUnfolding::getAxis(reco,(RooUnfolding::Dimension)i);
    double min = ax->GetBinLowEdge(!includeUnderflowOverflow);
    double max = ax->GetBinUpEdge(ax->GetNbins()+includeUnderflowOverflow);
    RooRealVar* obs = new RooRealVar(obs_reco_v[i],ax->GetTitle() ? ax->GetTitle() : obs_reco_v[i].Data(),min,min,max);
    setBinning(obs,ax,includeUnderflowOverflow);
    obs->setConstant(true);
    reco_vars.add(*obs);
  }
  
  this->setup(truth,truth_vars,reco,reco_vars,response,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg, const TH1* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  //! constructor
  this->setup(truth_th1,obs_truth,reco_th1,obs_reco,response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  //! constructor
  this->_bkg.setNominal(bkg);
  this->_data.setNominal(data,obs_reco);
  this->setup(truth_th1,obs_truth,reco_th1,obs_reco,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, const RooArgList& obs_truth, RooAbsReal* reco, const RooArgList& obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  //! constructor
  this->_reco.setNominal(reco);
  this->_bkg.setNominal(bkg);
  this->_data.setNominal(data,obs_reco);
  this->setup(truth_th1,obs_truth,NULL,obs_reco,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : RooUnfoldSpec(name,title,truth_th1,RooArgList(*obs_truth),reco,RooArgList(*obs_reco),response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity) {
  //! constructor
}

namespace {
  RooAbsReal* makeParamHistFunc(const char* name, const char* title, const RooArgList& observables, const RooArgList& parameters, bool multiplyDensity){
    ParamHistFunc* phf = new ParamHistFunc(name,title,observables,parameters);
    if(!multiplyDensity) return phf;
    RooArgList inverseBinWidths;
    RooRealVar* obs = (RooRealVar*)(observables.first());
    for(int i=0; i<obs->getBinning().numBins(); ++i){
      double width = obs->getBinWidth(i);
      RooRealVar* bw = new RooRealVar(TString::Format("%s_bin_%d_widthCorr",obs->GetName(),i),TString::Format("density correction of %s in bin %d",obs->GetTitle(),i),1./width);
      bw->setConstant(true);
      inverseBinWidths.add(*bw);
    }
    ParamHistFunc* bwcorr = new ParamHistFunc(TString::Format("%s_widthCorr",obs->GetName()),TString::Format("density corrections for %s",obs->GetTitle()),observables,inverseBinWidths);
    RooArgList elems;
    elems.add(*phf);
    elems.add(*bwcorr);
    RooProduct* prod = new RooProduct(TString::Format("%s_x_%s",phf->GetName(),bwcorr->GetName()),title,elems);
    return prod;
  }
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const RooArgList& reco_bins, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_reco.setNominal(::makeParamHistFunc(TString::Format("signal_reco_%s_differential",obs_reco->GetName()).Data(),obs_reco->GetTitle(),obs_reco_list,reco_bins,useDensity));
  this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()).Data(),obs_reco->GetTitle(),obs_reco_list,bkg_bins,useDensity));
  this->_data.setNominal(data,obs_reco_list);
  this->setup(truth_th1,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* bkg, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : RooUnfoldSpec(name,title,truth_th1,RooArgList(*obs_truth),reco_th1,RooArgList(*obs_reco),response_th1,bkg,data,includeUnderflowOverflow,errorThreshold,useDensity) {}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& bkg_bins, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);  
  this->_bkg.setNominal(::makeParamHistFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),obs_reco_list,bkg_bins,useDensity));
  this->_data.setNominal(data,obs_reco_list);
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, RooAbsReal* measured, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_data.setNominal(measured);
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, const TH1* reco_th1, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgList& measured_bins, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) : 
  TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_data.setNominal(::makeParamHistFunc(TString::Format("measured_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),obs_reco_list,measured_bins,useDensity));
  this->setup(truth_th1,obs_truth_list,reco_th1,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

namespace {
  RooRealSumFunc* makeRooRealSumFunc(const char* name, const char* title, const RooArgSet& contributions, double densityCorrection){
    RooRealVar* binWidth = new RooRealVar(TString::Format("%s_binWidthCorrection",name),"bin width correction",densityCorrection);
    binWidth->setConstant(true);
    RooArgList functions;
    RooArgList coefs;
    RooFIter itr(contributions.fwdIterator());
    RooAbsArg* obj;
    while((obj = itr.next())){
      functions.add(*obj);
      coefs.add(*binWidth);
    }
    return new RooRealSumFunc(name,title,functions,coefs);
  }
}


RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, const TH1* truth_th1, RooAbsArg* obs_truth, RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_reco.setNominal(reco);
  this->_data.setNominal(data,obs_reco_list);
  double densityCorr = 1;
  if(useDensity && obs_reco->InheritsFrom(RooRealVar::Class())){
    densityCorr = 1./((RooRealVar*)(obs_reco))->getBinning().averageBinWidth();
  }
  this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),bkg_contributions,densityCorr));
  this->setup(truth_th1,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}

RooUnfoldSpec::RooUnfoldSpec(const char* name, const char* title, RooAbsReal* truth, RooAbsArg* obs_truth,  RooAbsReal* reco, RooAbsArg* obs_reco, const TH2* response_th1, const RooArgSet& bkg_contributions, RooDataHist* data, bool includeUnderflowOverflow, double errorThreshold, bool useDensity) :
  TNamed(name,title)
{
  //! constructor
  RooArgList obs_reco_list(*obs_reco);
  RooArgList obs_truth_list(*obs_truth);
  this->_truth.setNominal(truth);
  this->_reco.setNominal(reco);
  this->_data.setNominal(data,obs_reco_list);
  double densityCorr = 1;
  if(useDensity && obs_reco->InheritsFrom(RooRealVar::Class())){
    densityCorr = 1./((RooRealVar*)(obs_reco))->getBinning().averageBinWidth();
  }
  this->_bkg.setNominal(::makeRooRealSumFunc(TString::Format("bkg_reco_%s_differential",obs_reco->GetName()),obs_reco->GetTitle(),bkg_contributions,densityCorr));
  this->setup(NULL,obs_truth_list,NULL,obs_reco_list,response_th1,NULL,NULL,includeUnderflowOverflow,errorThreshold,useDensity);
}


void RooUnfoldSpec::setup(const TH1* truth_th1, const RooArgList& obs_truth, const TH1* reco_th1, const RooArgList& obs_reco, const TH2* response_th1, const TH1* bkg_th1, const TH1* data_th1, bool includeUnderflowOverflow, double errorThreshold, bool useDensity){
  //! setup helper function
  this->_includeUnderflowOverflow = includeUnderflowOverflow;
  this->_useDensity = useDensity;
  this->_errorThreshold = errorThreshold;
  if(truth_th1) this->_truth.setNominal(truth_th1,obs_truth,errorThreshold,includeUnderflowOverflow,this->_useDensity);
  if(reco_th1)  this->_reco.setNominal(reco_th1,obs_reco,errorThreshold,includeUnderflowOverflow,this->_useDensity);
  this->_obs_reco.add(obs_reco);
  this->_obs_all.add(obs_reco);  
  this->_obs_truth.add(obs_truth);  
  this->_obs_all.add(obs_truth);
  if(response_th1) this->_res.setNominal(response_th1,this->_obs_all,errorThreshold,includeUnderflowOverflow,this->_useDensity);
  if(bkg_th1)  this->_bkg.setNominal(bkg_th1,obs_reco,errorThreshold,includeUnderflowOverflow,this->_useDensity);

  if(data_th1) this->_data.setNominal(data_th1,obs_reco,0.,includeUnderflowOverflow,this->_useDensity);
}

RooUnfoldSpec::~RooUnfoldSpec(){
  //! destructor
}


RooUnfolding::RooFitHist* RooUnfoldSpec::makeHistogram(const TH1* hist){
  // convert a TH1 into a RooFitHist
  return new RooUnfolding::RooFitHist(hist, this->_obs_truth, this->_includeUnderflowOverflow, this->_errorThreshold, this->_useDensity);
}


RooUnfolding::RooFitHist* RooUnfoldSpec::makeHistogram(const HistContainer& source, double errorThreshold){
  // build a new RooFitHist based on the source HistContainer. relative bin errors above errorThreshold will be modelled as gamma parameters.
  RooAbsReal* hf = source._nom;
  RooAbsReal* func = hf;
  if(!hf) return 0;
  std::vector<RooAbsArg*> obs;
  RooArgList obslist;
  if(this->_obs_all.getSize() < 1){
    throw std::runtime_error("in RooUnfoldSpec::makeHistogram: no observables known!");
  }
  RooFIter itr(this->_obs_all.fwdIterator());
  RooAbsArg* arg = NULL;
  while((arg = itr.next())){
    if(!arg) continue;
    if(!hf->dependsOn(*arg)) continue;
    obs.push_back(arg);
    obslist.add(*arg);
  }
  if(obslist.getSize() < 1){
    std::stringstream ss;
    ss << "in RooUnfoldSpec::makeHistogram: function '";
    hf->printStream(ss,0,RooPrintable::kStandard,"");
    ss << "' does not depend on any of the known observables '";
    this->_obs_all.printStream(ss,0,RooPrintable::kStandard,"");
    ss << "'!";
    throw std::runtime_error(ss.str());
  }
  if(source._shapes.size() > 0){
    RooArgList up, dn;
    RooArgList params;
    for(auto var:source._shapes){
      TString sysname(var.first);
      if(var.second.size() != 2){
        throw std::runtime_error(TString::Format("unable to process systematics '%s' with size %d != 2",var.first.c_str(),(int)(var.second.size())).Data());
      }
      up.add(*var.second[0]);
      dn.add(*var.second[1]);
      TString name = TString::Format("alpha_%s",var.first.c_str());
      RooRealVar* p = (RooRealVar*)(this->_alphas.find(name));
      if(!p){
        p = new RooRealVar(name,name,0,-5,5);
        p->setError(1);
        this->addGaussNP(p);
      }
      params.add(*p);
    }
    TString name = TString::Format("%s_HistoSystematics",hf->GetName());
    hf->SetName(TString::Format("%s_Nominal",hf->GetName()));
    func = new PiecewiseInterpolation(name.Data(),name.Data(),*hf,up,dn,params);
  }
  RooArgList components; 
  if(source._norms.size() > 0){
    std::vector<double> up,dn;
    RooArgList params;
    for(auto var:source._norms){
      TString sysname(var.first);
      up.push_back(var.second.first);
      dn.push_back(var.second.second);
      TString name = TString::Format("alpha_%s",var.first.c_str());
      RooRealVar* p = (RooRealVar*)(this->_alphas.find(name));
      if(!p){
        p = new RooRealVar(name,name,0,-5,5);
        p->setError(1);
        this->addGaussNP(p);
      }
      params.add(*p);
    }
    TString name = TString::Format("%s_OverallSystematics",hf->GetName());
    components.add(*(new RooStats::HistFactory::FlexibleInterpVar(name.Data(),name.Data(),params,1.,up,dn)));
  }
  for(auto g:source._gammas){
    this->addPoissonNP(g);
  }
  if(source._staterror) components.add(*source._staterror);
  if(components.getSize() > 0){
    components.add(*func);
    TString name(hf->GetName());
    hf->SetName(TString::Format("%s_hist",hf->GetName()));
    func = new RooProduct(name.Data(),hf->GetTitle(),components);  
  }
  return new RooUnfolding::RooFitHist(func,obs,source._gammas);
}

RooHistFunc* RooUnfoldSpec::makeHistFuncTruth(const TH1* hist){
  //! create a new truth hist func
  if(!hist) return NULL;
  return RooUnfolding::makeHistFunc(hist,this->_obs_truth, this->_includeUnderflowOverflow, this->_useDensity);
}

RooHistFunc* RooUnfoldSpec::makeHistFuncMeasured(const TH1* hist){
  //! create a new measured hist func
  if(!hist) return NULL;
  return RooUnfolding::makeHistFunc(hist,this->_obs_reco, this->_includeUnderflowOverflow, this->_useDensity);
}

RooProdPdf* RooUnfoldSpec::makeConstraints(){
  //! create all the constraint terms
  RooArgList constraints;
  for(auto a:this->_alphas){
    RooRealVar* p = (RooRealVar*)(a);
    RooRealVar* mean = new RooRealVar(TString::Format("%s_nom",p->GetName()),TString::Format("Mean value for %s",p->GetName()),p->getVal());
    RooRealVar* sigma = new RooRealVar(TString::Format("%s_sigma",p->GetName()),TString::Format("Sigma for %s",p->GetName()),p->getError());
    mean->setConstant(true);
    sigma->setConstant(true);
    RooGaussian* gaus = new RooGaussian(TString::Format("%s_constraint",p->GetName()),TString::Format("Gaussian constraint term for %s",p->GetName()),*p,*mean,*sigma);
    constraints.add(*gaus);
  }
  for(auto g:this->_gammas){
    RooRealVar* p = (RooRealVar*)(g);
    double n = pow(p->getError(),-2);
    RooRealVar* tau = new RooRealVar(TString::Format("%s_tau",p->GetName()),TString::Format("Tau parameter value for %s",p->GetName()),n);
    tau->setConstant(true);
    RooArgList params(*p,*tau);
    RooProduct* prod = new RooProduct(TString::Format("%s_nEvents",p->GetName()),TString::Format("Number of events for %s",p->GetName()),params);
    RooPoisson* pois = new RooPoisson(TString::Format("%s_constraint",p->GetName()),TString::Format("Poisson constraint term for %s",p->GetName()),*prod,*tau);
    constraints.add(*pois);
  }
  return new RooProdPdf(TString::Format("%s_constraints",this->GetName()),"Unfolding constraint terms",constraints);
}

void RooUnfoldSpec::addGaussNP(RooRealVar* v){
  //! add a new gaussian NP
  if(v) this->_alphas.add(*v);
}
void RooUnfoldSpec::addPoissonNP(RooRealVar* v){ 
  //! add a new poisson NP
 if(v) this->_gammas.add(*v);
}

void RooUnfoldSpec::makeBackground(){
  //! create the background
  this->_locked = true;
  if(!this->_cache._bkg){
    this->_cache._bkg = this->makeHistogram(this->_bkg,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeData(){
  //! create the data
  this->_locked = true;
  if(!this->_cache._data){
    this->_cache._data = this->makeHistogram(this->_data,0);
  }
}
void RooUnfoldSpec::makeResponse(){
  //! create the response
  this->_locked = true;
  if(!this->_cache._res){
    this->makeReco();
    this->makeTruth();
    if(!this->_res._nom) throw std::runtime_error("no response input given!");        
    this->_cache._res = this->makeHistogram(this->_res,this->_errorThreshold);
    this->_cache._response = new RooFitUnfoldResponse(this->GetName(),this->GetTitle(),this->_cache._res,this->_cache._truth,this->_cache._reco,this->_useDensity);
  }
}
void RooUnfoldSpec::makeTruth(){
  //! create the truth
  this->_locked = true;
  if(!this->_cache._truth){
    if(!this->_truth._nom) throw std::runtime_error("no truth input given!");    
    this->_cache._truth = this->makeHistogram(this->_truth,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeReco(){
  //! create the reconstructed
  this->_locked = true;
  if(!this->_cache._reco){
    if(!this->_reco._nom) throw std::runtime_error("no measured input given!");
    this->_cache._reco = this->makeHistogram(this->_reco,this->_errorThreshold);
  }
}
void RooUnfoldSpec::makeDataMinusBackground(){
  //! create the data minus background
  this->_locked = true;
  this->makeData();
  if(!this->_cache._data_minus_bkg){
    if(this->_bkg._nom){
      this->makeResponse();
      this->makeBackground();
      this->_cache._data_minus_bkg = this->_cache._response->makeHistSum(this->_cache._data->func(),this->_cache._bkg->func(),1.,-1.);
    } else {
      this->_cache._data_minus_bkg = this->_cache._data;
    }
    TString name(TString::Format("%s_data_minus_bkg",this->GetName()));
    this->_cache._data_minus_bkg->func()->SetName(name);
    this->_cache._data_minus_bkg->func()->SetTitle(name);
  }
}

RooAbsReal* RooUnfoldSpec::getBackground(){
  //! retrieve the background
  this->makeBackground(); if(this->_cache._bkg) return this->_cache._bkg->func(); else return NULL;
}
RooAbsReal* RooUnfoldSpec::getData(){
  //! retrieve the background
  this->makeData(); return this->_cache._data->func(); 
}
RooAbsReal* RooUnfoldSpec::getResponse(){
  //! retrieve the response
  this->makeResponse(); return this->_cache._res->func(); 
}
RooAbsReal* RooUnfoldSpec::getTruth(){
  //! retrieve the truth
  this->makeTruth(); return this->_cache._truth->func(); 
}
RooAbsReal* RooUnfoldSpec::getReco(){
  //! retrieve the reconstructed
  this->makeReco(); return this->_cache._reco->func(); 
}
RooAbsReal* RooUnfoldSpec::getDataMinusBackground(){
  //! retrieve the data minus background
  this->makeDataMinusBackground(); return this->_cache._data_minus_bkg->func(); 
}

RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* RooUnfoldSpec::unfold(Algorithm alg, Double_t regparam){
  //! create the unfolding object
  this->makeResponse();
  this->makeDataMinusBackground();
  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfolding = RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>::New(alg,this->_cache._response,this->_cache._data_minus_bkg,regparam);
  if (this->_cache._bkg) {
    unfolding->SetBkg(this->_cache._bkg);
  }
  unfolding->SetTruth(this->_cache._truth);
  unfolding->SetOverflow(this->_includeUnderflowOverflow);

  return unfolding;
}


void RooUnfoldSpec::HistContainer::setNominal(RooAbsReal* nom){
  this->_nom = nom;
}

void RooUnfoldSpec::HistContainer::setNominal(RooDataHist* data, const RooArgList& obslist){
  this->_nom = RooUnfolding::makeHistFunc(data,obslist);
  this->_gammas = RooUnfolding::createGammas(data,obslist,0.);
  if(_gammas.size()>0){
    this->_staterror = RooUnfolding::makeParamHistFunc(TString::Format("%s_staterrors",_nom->GetName()),_nom->GetTitle(),obslist,_gammas);
  }
}

void RooUnfoldSpec::HistContainer::setNominal(const TH1* nom, const RooArgList& obslist,double errorThreshold, bool includeUnderflowOverflow, bool useDensity){
  this->_nom = RooUnfolding::makeHistFunc(nom,obslist,includeUnderflowOverflow,useDensity);  
  if(errorThreshold >= 0){
    this->_gammas = RooUnfolding::createGammas(nom,includeUnderflowOverflow,errorThreshold);
    if(_gammas.size()>0){
      this->_staterror = RooUnfolding::makeParamHistFunc(TString::Format("%s_staterrors",nom->GetName()),nom->GetTitle(),obslist,_gammas);
    }
  }
}


void RooUnfoldSpec::HistContainer::addShape(const char* name, RooAbsReal* up, RooAbsReal* dn){
  this->_shapes[name] = {up,dn};
}

void RooUnfoldSpec::HistContainer::addNorm(const char* name, double up, double dn){
  this->_norms[name] = {up,dn};
}

RooUnfoldSpec::HistContainer::~HistContainer(){}

void RooUnfoldSpec::lockCheck(){
  if(this->_locked){
    throw std::runtime_error("this instance of RooUnfoldSpec is locked - it has already been used to produce results and can no longer be modified. please create a new instance for modifications!");
  }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char* name, const TH1* up, const TH1* down){
  //! register a new shape systematic
  this->lockCheck();
  switch(c){
  case kTruth:
    this->_truth.addShape(name,
                          RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_up",up->GetName(),name),up,this->_obs_truth,this->_includeUnderflowOverflow,this->_useDensity),
                          RooUnfolding::makeHistFunc(TString::Format("truth_%s_%s_dn",down->GetName(),name),down,this->_obs_truth,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kMeasured:
    this->_reco.addShape(name,
                         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                         RooUnfolding::makeHistFunc(TString::Format("meas_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kData:
    this->_data.addShape(name,
                         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                         RooUnfolding::makeHistFunc(TString::Format("data_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kResponse:
    this->_res.addShape(name,
                        RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_up",up->GetName(),name),up,this->_obs_all,this->_includeUnderflowOverflow,this->_useDensity),
                        RooUnfolding::makeHistFunc(TString::Format("resp_%s_%s_dn",down->GetName(),name),down,this->_obs_all,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  case kBackground:
    this->_bkg.addShape(name,
                        RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_up",up->GetName(),name),up,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity),
                        RooUnfolding::makeHistFunc(TString::Format("bkg_%s_%s_dn",down->GetName(),name),down,this->_obs_reco,this->_includeUnderflowOverflow,this->_useDensity));
    break;
  }
}

void RooUnfoldSpec::registerSystematic(Contribution c, const char* name, double up, double down){
  //! register a new normalization systematic
  this->lockCheck();
  switch(c){
  case kTruth:
    this->_truth.addNorm(name,up,down);
    break;
  case kMeasured:
    this->_reco.addNorm(name,up,down);
    break;
  case kData:
    this->_data.addNorm(name,up,down);
    break;
  case kResponse:
    this->_res.addNorm(name,up,down);
    break;
  case kBackground:
    this->_bkg.addNorm(name,up,down);
    break;
  }
}




RooAbsPdf* RooUnfoldSpec::makePdf(Algorithm alg, Double_t regparam){
  //! create an unfolding pdf
  RooUnfoldFunc* func = static_cast<RooUnfoldFunc*>(this->makeFunc(alg,regparam));
  func->setDensity(true);
  RooWrapperPdf* pdf = new RooWrapperPdf(TString::Format("%s_pdf",func->GetName()),TString::Format("%s Pdf",func->GetTitle()),*func);
  RooAbsReal* integral = func->createIntegral(this->_obs_truth);
  RooExtendPdf* extpdf = new RooExtendPdf(TString::Format("%s_extpdf",func->GetName()),TString::Format("%s Extended Pdf",func->GetTitle()),*pdf,*integral);
  
  RooProdPdf* constraints = this->makeConstraints();
  RooArgList comps(*extpdf,*constraints);
  RooProdPdf* prod = new RooProdPdf(TString::Format("%s_x_constraints",this->GetName()),"Unfolding pdf, including constraints",comps);
  prod->setStringAttribute("source",func->GetName());
  return prod;
}


RooAbsReal* RooUnfoldSpec::makeFunc(Algorithm alg, Double_t regparam){
  //! create an unfolding function
  RooUnfoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>* unfold = this->unfold(alg, regparam);
  RooAbsReal* func = new RooUnfoldFunc(this->GetName(),this->GetTitle(),this->unfold(alg, regparam),false);
  func->setStringAttribute("source",func->GetName());
  delete unfold;
  return func;
}

ClassImp(RooUnfoldSpec)

#endif


