#ifndef NOROOFIT
#include "RooUnfoldHelpers.h"
#include "RooUnfoldFitHelpers.h"
#include "TH1.h"
#include "TAxis.h"
#include "RooUnfoldTH1Helpers.h"

#include <iostream>
#include <limits>

#include "TMath.h"
#include "RooAbsReal.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAbsDataStore.h"
#include "RooHistFunc.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooProduct.h"
#include "RooParamHistFunc.h"
#include "RooStats/HistFactory/ParamHistFunc.h"

#define nan std::numeric_limits<double>::quiet_NaN()

#include "RooPoisson.h"
#include "RooProduct.h"
#include "RooWorkspace.h"
#include "TPRegexp.h"

namespace {
  // hack to get access to RooPoisson
  template <typename RooPoissonTag>
  struct RooPoissonHackResult {
    typedef typename RooPoissonTag::type type;
    static type ptr;
  };
  
  template <typename RooPoissonTag>
  typename RooPoissonHackResult<RooPoissonTag>::type RooPoissonHackResult<RooPoissonTag>::ptr;
  
  template<typename RooPoissonTag, typename RooPoissonTag::type p>
  struct RooPoissonRob : RooPoissonHackResult<RooPoissonTag> {
    struct RooPoissonFiller {
      RooPoissonFiller() {RooPoissonHackResult<RooPoissonTag>::ptr = p;}
    };
    static RooPoissonFiller roopoissonfiller_obj;
  };
  
  template<typename RooPoissonTag, typename RooPoissonTag::type p>
  typename RooPoissonRob<RooPoissonTag, p>::RooPoissonFiller RooPoissonRob<RooPoissonTag, p>::roopoissonfiller_obj;
  
  //now expose some members of RooPoisson that we need to access
  struct RooPoissonMean { typedef RooRealProxy(RooPoisson::*type); };
  template struct RooPoissonRob<RooPoissonMean, &RooPoisson::mean>;
  //now expose some members of RooPoisson that we need to access
  struct RooPoissonX { typedef RooRealProxy(RooPoisson::*type); };
  template struct RooPoissonRob<RooPoissonX, &RooPoisson::x>;
  
  RooRealVar* getGammaParameter(RooPoisson* pois){
    const RooRealProxy& px = pois->*RooPoissonHackResult<RooPoissonMean>::ptr;
    RooProduct* prod = const_cast<RooProduct*>(dynamic_cast<const RooProduct*>(&(px.arg())));
    RooArgList compList(prod->components());
    RooFIter citr(compList.fwdIterator());
    RooAbsArg* arg;
    while((arg = citr.next())){
      RooRealVar* gamma = dynamic_cast<RooRealVar*>(arg);
      if(gamma) return gamma;
    }
    return NULL;
  }

  double getMean(RooPoisson* pois){
    const RooRealProxy& px = pois->*RooPoissonHackResult<RooPoissonMean>::ptr;
    return px.arg().getVal();
  }
}

namespace {
  // hack to get access to RooHistFunc
  template <typename RooHistFuncTag>
  struct RooHistFuncHackResult {
    typedef typename RooHistFuncTag::type type;
    static type ptr;
  };
  
  template <typename RooHistFuncTag>
  typename RooHistFuncHackResult<RooHistFuncTag>::type RooHistFuncHackResult<RooHistFuncTag>::ptr;
  
  template<typename RooHistFuncTag, typename RooHistFuncTag::type p>
  struct RooHistFuncRob : RooHistFuncHackResult<RooHistFuncTag> {
    struct RooHistFuncFiller {
      RooHistFuncFiller() {RooHistFuncHackResult<RooHistFuncTag>::ptr = p;}
    };
    static RooHistFuncFiller roohistfuncfiller_obj;
  };
  
  template<typename RooHistFuncTag, typename RooHistFuncTag::type p>
  typename RooHistFuncRob<RooHistFuncTag, p>::RooHistFuncFiller RooHistFuncRob<RooHistFuncTag, p>::roohistfuncfiller_obj;
  
  //now expose some members of RooHistFunc that we need to access
  struct RooHistFuncObs { typedef RooArgSet(RooHistFunc::*type); };
  template struct RooHistFuncRob<RooHistFuncObs, &RooHistFunc::_histObsList>;
  
}

namespace {
  template<class T> T* findLeafServer(RooAbsArg* rr, const char* name){
    if(!rr) return NULL;
    RooFIter sIter = rr->serverMIterator();
    RooAbsArg* server = NULL;
    while ((server=sIter.next())) {
      if(!server) continue;
      if(strcmp(server->GetName(),name) == 0 && (server->InheritsFrom(RooRealVar::Class()) || server->InheritsFrom(RooCategory::Class()))){
        return dynamic_cast<T*>(server);
      } else {
        T* rrv = findLeafServer<T>(server,name);
        if(rrv) return rrv;
      }
    }
    return NULL;
  }

  template<class T> RooArgList argList(const std::vector<T*>& vars){
    RooArgList retval;
    for(auto v:vars) retval.add(*v,true);
    return retval;
  }
}  


  


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace {
  int binNumber(RooAbsArg* arg, double x){
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getBinning().binNumber(x);
    } else {
      return x;
    }
  }
  int nBins(RooAbsArg* arg){
    if(!arg) return 0;
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getBinning().numBins();
    } else if(arg->InheritsFrom(RooCategory::Class())){
      return ((RooCategory*)(arg))->numTypes();
    }
    throw std::runtime_error("unknown argument type for nBins!");
  }
  double min(RooAbsArg* arg){
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getMin();
    } else {
      return 0.;
    }
  }
  double max(RooAbsArg* arg){
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getMax();
    } else if(arg->InheritsFrom(RooCategory::Class())){
      return ((RooCategory*)(arg))->numTypes()-1;      
    }
    throw std::runtime_error("unknown argument type for max!");    
  }
  double binCenter(RooAbsArg* arg, int i){
    if(arg->InheritsFrom(RooRealVar::Class())){    
      return ((RooRealVar*)(arg))->getBinning().binCenter(i);
    } else {
      return i;
    }
  }
  double binWidth(RooAbsArg* arg, int i){
    if(arg->InheritsFrom(RooRealVar::Class())){    
      return ((RooRealVar*)(arg))->getBinning().binWidth(i);
    } else {
      return i;
    }
  }
  double binHighEdge(RooAbsArg* arg, int i){
    if(arg->InheritsFrom(RooRealVar::Class())){    
      return ((RooRealVar*)(arg))->getBinning().binHigh(i);
    } else {
      return i;
    }
  }
  double binLowEdge(RooAbsArg* arg, int i){
    if(arg->InheritsFrom(RooRealVar::Class())){    
      return ((RooRealVar*)(arg))->getBinning().binLow(i);
    } else {
      return i;
    }
  }
  double getVal(RooAbsArg* arg){
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getVal();
    } else if(arg->InheritsFrom(RooCategory::Class())){
      return ((RooCategory*)(arg))->getIndex();
    }
    throw std::runtime_error("unknown argument type for getVal!");                
  }
  void setVal(RooAbsArg* arg, double val){
    if(arg->InheritsFrom(RooRealVar::Class())){
      ((RooRealVar*)(arg))->setVal(val);
    } else if(arg->InheritsFrom(RooCategory::Class())){
      ((RooCategory*)(arg))->setIndex(val);
    } else {
      throw std::runtime_error("unknown argument type for setVal");
    }
  }
  int getBin(RooAbsArg* arg){
    if(arg->InheritsFrom(RooRealVar::Class())){
      return ((RooRealVar*)(arg))->getBinning().binNumber(((RooRealVar*)(arg))->getVal());
    } else if(arg->InheritsFrom(RooCategory::Class())){
      return ((RooCategory*)(arg))->getIndex();
    }
    throw std::runtime_error("unknown argument type for getBin!");        
  }  
}


namespace RooUnfolding { // section 1: trivial helpers
  
  Variable<RooUnfolding::RooFitHist>::Variable(int nBins,double min,double max,const char* name) : _var(new RooRealVar(name,name,nBins,min,max)) {}
  Variable<RooUnfolding::RooFitHist>::Variable(RooAbsArg* var) : _var(var) {}

  template<> Variable<RooUnfolding::RooFitHist> var(const RooUnfolding::RooFitHist* h, Dimension d){
    return Variable<RooUnfolding::RooFitHist>(h->obs(d));
  }

  template<> const char* name<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h){
    return h->name();
  }
  template<> const char* title<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h){
    return h->title();
  }

  template<> int findBin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h, double x, RooUnfolding::Dimension d){
    if(d > h->dim()) throw std::runtime_error(TString::Format("unable to retrieve min for dimension %d for histogram %s with %d dimensions",d,name(h),(int)(h->dim())).Data());
    return ::binNumber(h->obs(d),x);
  }
  template<> int findBin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h, Double_t x){
    if(h->dim() != 1) throw std::runtime_error(TString::Format("inacceptable call: findBin with one coordinate for %d dimensional histogram",(int)(h->dim())).Data());
    return ::binNumber(h->obs(0),x);
  }
  template<> int findBin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h, Double_t x, Double_t y){
    if(h->dim() != 2) throw std::runtime_error(TString::Format("inacceptable call: findBin with two coordinates for %d dimensional histogram",(int)(h->dim())).Data());
    return ::binNumber(h->obs(0),x) + ::nBins(h->obs(0)) * ::binNumber(h->obs(1),y);
  }
  template<> int findBin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h, Double_t x, Double_t y, Double_t z){
    if(h->dim() != 3) throw std::runtime_error(TString::Format("inacceptable call: findBin with three coordinates for %d dimensional histogram",(int)(h->dim())).Data());
    return ::binNumber(h->obs(0),x) + ::nBins(h->obs(0)) * (::binNumber(h->obs(1),y) + ::nBins(h->obs(1)) * ::binNumber(h->obs(2),z));
  }

  template<> double min<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, RooUnfolding::Dimension d){
    if(d > hist->dim()) throw std::runtime_error(TString::Format("unable to retrieve min for dimension %d for histogram %s with %d dimensions",d,name(hist),(int)(hist->dim())).Data());
    return ::min(hist->obs(d));    
  }
  template<> double max<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, RooUnfolding::Dimension d){
    if(d > hist->dim()) throw std::runtime_error(TString::Format("unable to retrieve max for dimension %d for histogram %s with %d dimensions",d,name(hist),(int)(hist->dim())).Data());
    return ::max(hist->obs(d));
  }
  template<> int sumW2N<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist){
    return hist->weighted();
  }
  template<> bool empty<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* /*hist*/){
    return false;
  }
  template<> int dim<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist){
    // get number of observables
    return hist->dim();
  }
  template<> int nBins<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, bool /*overflow*/){
    // get number of bins for all observables
    int n = 1;
    for(size_t i=0; i<hist->dim(); ++i){
      n *= ::nBins(hist->obs(i));
    }
    return n;
  }
  template<> int nBins<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, RooUnfolding::Dimension d, bool /*overflow*/){
    // get number of bins for one observable
    if(hist->dim() <= d) return 1;
    return ::nBins(hist->obs(d));
  }
  template<> double binCenter<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist*hist, int i, RooUnfolding::Dimension d){
    if(hist->dim() <= d) return nan;
    return ::binCenter(hist->obs(d),i);
  }
  template<> double binWidth<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist*hist, int i, RooUnfolding::Dimension d){
    if(hist->dim() <= d) return nan;
    return ::binWidth(hist->obs(d),i);
  }
  template<> double binHighEdge<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist*hist, int i, RooUnfolding::Dimension d){
    if(hist->dim() <= d) return nan;
    return ::binHighEdge(hist->obs(d),i);
  }
  template<> double binLowEdge<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist*hist, int i, RooUnfolding::Dimension d){
    if(hist->dim() <= d) return nan;
    return ::binLowEdge(hist->obs(d),i);
  }
  template<> int bin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* /*hist*/, int i, Bool_t /*overflow*/){
    return i;
  }
  template<> int bin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, int i, int j, Bool_t /*overflow*/){
    int nx = nBins(hist,X);
    return j*nx+i;
  }
  template<> int bin<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, int i, int j, int k, Bool_t /*overflow*/){
    int nx = nBins(hist,X);
    int ny = nBins(hist,Y);    
    return k*nx*ny+j*nx+i;
  }    
  template<> void binXYZ<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, int binglobal, int& binx, int& biny, int& binz){
    // stole code from TH1::GetBinXYZ
    int d = dim(hist);
    int nx = nBins(hist,X);
    if (d == 1) {
       binx = binglobal%nx;
       biny = 0;
       binz = 0;
       return;
    } else if (d > 1) {
      int ny = nBins(hist,Y);
      if(d == 2){
        binx = binglobal%nx;
        biny = ((binglobal-binx)/nx)%ny;
        binz = 0;
        return;
      } else if(d == 3) {
        binx = binglobal%nx;
        biny = ((binglobal-binx)/nx)%ny;
        binz = ((binglobal-binx)/nx -biny)/ny;
      }
    }
  }

}
namespace { // interjunction: some additional helpers
  void setBin(RooAbsArg* arg, int ibin){
    if(arg->InheritsFrom(RooRealVar::Class())){
      ((RooRealVar*)(arg))->setVal(((RooRealVar*)(arg))->getBinning().binCenter(ibin));
    } else if(arg->InheritsFrom(RooCategory::Class())){
      ((RooCategory*)(arg))->setIndex(ibin);
    }
  }
  void setBin(const RooUnfolding::RooFitHist* h, int ibin){
    int d = RooUnfolding::dim(h);
    int ix,iy,iz;
    binXYZ(h,ibin,ix,iy,iz);
    setBin(h->obs(0),ix);
    if(d>1) setBin(h->obs(1),iy);
    if(d>2) setBin(h->obs(2),iz);
    if(d>3) throw std::runtime_error(TString::Format("unable to handle histograms %s: %d dimensions!",h->name(),d).Data());
    if(h->bin()!=ibin) throw std::runtime_error("internal binning error!");
  }
}

namespace RooUnfolding {
  template<> const char* varname<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* hist, Dimension d){
    if(d > hist->dim()) throw std::runtime_error(TString::Format("unable to retrieve name for dimension %d for histogram %s with %d dimensions",d,name(hist),(int)(hist->dim())).Data());
    return hist->obs(d)->GetName();        
  }
  template<> double binError<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h, Int_t i, Bool_t /*overflow*/)
  {
    setBin(h,i);
    return h->error();
  }  
  template<> double binContent<RooUnfolding::RooFitHist> (const RooUnfolding::RooFitHist* h, Int_t i, Bool_t /*overflow*/){
    setBin(h,i);
    return h->value();
  }
  template<> double binVolume<RooUnfolding::RooFitHist> (const RooUnfolding::RooFitHist* h, Int_t bin, Bool_t /*overflow*/){
    double v = 1;
    std::vector<int> ibins = {0,0,0};
    binXYZ(h,bin,ibins[0],ibins[1],ibins[2]);
    for(size_t d=0; d<h->dim(); ++d){
      v *= ::binWidth(h->obs(d),ibins[d]);
    }
    return v;
  }

  template<> void printHistogram<RooUnfolding::RooFitHist>(const RooUnfolding::RooFitHist* h){
    
    if(!h){
      std::cout << "RooFitHist (NULL)" << std::endl;
      return;
    }
    size_t n = nBins(h);
    
    std::cout << "RooFitHist " << name(h) << " @ " << h << " (" << h->func()->ClassName() << " @ " << h->func() << "): " << n << " bins in " << dim(h) << " dimensions (";
    
    for(int i=0; i<dim(h); ++i){
      if(i!=0) std::cout << ",";
      std::cout << h->obs(i)->GetName() << "@" << h->obs(i);
    }
    std::cout << ")";
    
    int nnps = h->nps().size();
    std::cout << " " << nnps << " bins affected by NPs" << std::endl;
    h->func()->Print("t");
    for(size_t i=0; i<n; ++i){
      setBin(h,i);
      std::cout << h->bin() << " " << h->value() << " " << h->error() << std::endl;
    }
  }
}

namespace { // interjunction: some additional helpers
  std::vector<RooRealVar*> findParameters(RooUnfolding::RooFitHist* hist, int ibin, const std::vector<RooRealVar*>& c){
    std::vector<RooRealVar*> vec;
    RooAbsReal* f = hist->func();
    for(auto p:c){
      if(!p || !f->dependsOn(*p)) continue;
      setBin(hist,ibin);
      double cval = p->getVal();
      double nomval = f->getVal();
      p->setVal(cval+fabs(p->getError()));
      setBin(hist,ibin);
      double upval = f->getVal();
      p->setVal(cval);
      if(!TMath::AreEqualRel(upval,nomval,1e-9)){
        vec.push_back(p);
      }
    }
    return vec;
  }
  
  void findDependantsPerBin(RooUnfolding::RooFitHist* hist, const std::vector<RooRealVar*>& allvars, std::vector<RooRealVar*>& out){
    if(!hist) throw std::runtime_error("unable to find dependants of invalid histogram!");
    size_t n = RooUnfolding::nBins<RooUnfolding::RooFitHist>(hist,false);
    bool found = false;
    for(size_t i=0; i<n; ++i){
      std::vector<RooRealVar*> params = findParameters(hist,i,allvars);
      if(params.size() == 0){
        out.push_back(NULL);
      } else if(params.size() == 1){
        found = true;
        out.push_back(params[0]);
      } else {
        throw std::runtime_error(TString::Format("found %d parameters influencing bin %d of histogram %s!",(int)(params.size()),(int)i,hist->name()).Data());
      }
    }
    if(!found) out.clear();
  }
  template<class T>void replace(T*& old, const RooAbsCollection& newlist){
    if(old){
      T* newObj = dynamic_cast<T*>(newlist.find(old->GetName()));
      if(newObj) old=newObj;
    }
  }
  inline double useIf(double val, bool selector){
    return selector ? val : 1.;
  }
}

namespace RooUnfolding { // section 2: non-trivial helpers

  const char* RooFitHist::name() const { return this->_func->GetName(); }
  const char* RooFitHist::title() const { return this->_func->GetTitle(); }
  RooAbsArg* RooFitHist::obs(size_t i) const { if(i>this->_obs.size()) throw std::runtime_error("attempt to access invalid observable!"); return this->_obs[i]; }
  size_t RooFitHist::dim() const { return this->_obs.size(); }
  
  const char* RooFitHist::GetName() const { return this->name(); }
  const char* RooFitHist::GetTitle() const { return this->title(); }
  void RooFitHist::Print(const char* /*opts*/) const { std::cout << this->GetName() << " " << this->GetTitle() << std::endl; this->printHist(); }
  
  void RooFitHist::saveSnapshot(std::map<std::string,double>& snsh) const {
    for(size_t i=0; i<this->_obs.size(); ++i){
      snsh[this->_obs[i]->GetName()] = ::getVal(this->_obs[i]);
    }
  }
  void RooFitHist::loadSnapshot(const std::map<std::string,double>& snsh){
    for(size_t i=0; i<this->_obs.size(); ++i){
      ::setVal(this->_obs[i],snsh.at(this->_obs[i]->GetName()));
    }
  }
    
  int RooFitHist::bin() const { 
    size_t offset = 1;
    size_t ibin = 0;
    for(size_t i=0; i<this->dim(); ++i){
      RooAbsArg* x = this->obs(i);
      size_t ix = ::getBin(x);
      ibin += offset*ix;
      offset *= ::nBins(x);
    }
    return ibin;
  }

  RooAbsReal* RooFitHist::func() const { return this->_func; }

  RooFitHist* RooFitHist::clone() const { return new RooFitHist(*this); }

  RooFitHist::RooFitHist() : _func(0) {
    // default constructor
  }

  RooFitHist::RooFitHist(const RooFitHist* h) : _func(h->_func), _obs(h->_obs) {
    // for(auto obs:_obs) if(!obs) obs->GetName();
  }
  RooFitHist::RooFitHist(RooAbsReal* f, const std::vector<RooAbsArg*>& obsvec) : _func(f), _obs(obsvec) {
    // for(auto obs:_obs) if(!obs) obs->GetName();
  }
  RooFitHist::RooFitHist(RooAbsReal* f, RooAbsArg* obsarg) : _func(f) {
    _obs.push_back(obsarg);
    // for(auto obs:_obs) if(!obs) obs->GetName();
  }
  RooFitHist::RooFitHist(RooAbsReal* f, RooAbsArg* obs1, RooAbsArg* obs2) : _func(f) { 
    _obs.push_back(obs1); _obs.push_back(obs2);
    // for(auto obs:_obs) if(!obs) obs->GetName();
  } 
  RooFitHist::RooFitHist(RooAbsReal* f, const std::vector<RooAbsArg*>& obs, const std::vector<RooRealVar*>& nps) : RooFitHist(f,obs) {
    findDependantsPerBin(this,nps,this->_gamma);
  }
  RooFitHist::RooFitHist(RooAbsReal* f, RooAbsArg* obs, const std::vector<RooRealVar*>& nps)  : RooFitHist(f,obs) {
    findDependantsPerBin(this,nps,this->_gamma);
  }
  RooFitHist::RooFitHist(RooAbsReal* f, RooAbsArg* obs1, RooAbsArg* obs2, const std::vector<RooRealVar*>& nps) :  RooFitHist(f,obs1,obs2) { 
    findDependantsPerBin(this,nps,this->_gamma);
  }

  double getIntegral(const TH1* histo, bool includeUnderflowOverflow, bool correctDensity){
    int offset = !includeUnderflowOverflow;    
    double integral = 0.;
    for (int ix=0 ; ix < histo->GetNbinsX() ; ix++) {
      if (histo->GetNbinsY() > 1) {
        for (int iy=0 ; iy < histo->GetNbinsY() ; iy++) {
          if (histo->GetNbinsZ() > 1) {          
            for (int iz=0 ; iz < histo->GetNbinsZ() ; iz++) {
              int bin = histo->GetBin(ix+offset,iy+offset,iz+offset);
              double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset)*histo->GetYaxis()->GetBinWidth(iy+offset)*histo->GetZaxis()->GetBinWidth(iz+offset),correctDensity);
              integral +=histo->GetBinContent(bin)/volume;
            }
          } else {
            int bin = histo->GetBin(ix+offset,iy+offset);
            double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset)*histo->GetYaxis()->GetBinWidth(iy+offset),correctDensity);
            integral +=histo->GetBinContent(bin)/volume;            
          }
        }
      } else {
        int bin = histo->GetBin(ix+offset);
        double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset),correctDensity);
        integral +=histo->GetBinContent(bin)/volume;                    
      }
    }
    return integral;
  }
  
  RooDataHist* convertTH1(const TH1* histo, const std::vector<RooAbsArg*>& vars, bool includeUnderflowOverflow, bool correctDensity, double scale){
    return convertTH1(histo,argList(vars),includeUnderflowOverflow,correctDensity,scale);
  }

  RooDataHist* convertTH1(const TH1* histo, const RooArgList& obs, bool includeUnderflowOverflow, bool correctDensity, double scale){
    TString name(histo->GetName());
    TString title(histo->GetTitle());
   
    // Define x,y,z as 1st, 2nd and 3rd observable
    RooAbsArg* xvar = obs.at(0);
    RooAbsArg* yvar = obs.at(1);
    RooAbsArg* zvar = obs.at(2);
 
    RooArgSet args(obs);
    RooDataHist* dh = new RooDataHist(TString::Format("%s_hist",name.Data()),title,args);

    int offset = !includeUnderflowOverflow;

    // Transfer contents
    for (int ix=0 ; ix < ::nBins(xvar) ; ix++) {
      ::setBin(xvar,ix) ;
      if (yvar) {
        for (int iy=0 ; iy < ::nBins(yvar) ; iy++) {
          ::setBin(yvar,iy) ;
          if (zvar) {
            for (int iz=0 ; iz < ::nBins(zvar) ; iz++) {
              ::setBin(zvar,iz) ;
              int bin = histo->GetBin(ix+offset,iy+offset,iz+offset);
              double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset)*histo->GetYaxis()->GetBinWidth(iy+offset)*histo->GetZaxis()->GetBinWidth(iz+offset),correctDensity);
              dh->add(obs,scale*histo->GetBinContent(bin)/volume,TMath::Power(scale*histo->GetBinError(bin)/volume,2)) ;
            }
          } else {
            int bin = histo->GetBin(ix+offset,iy+offset);
            double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset)*histo->GetYaxis()->GetBinWidth(iy+offset),correctDensity);
            dh->add(obs,scale*histo->GetBinContent(bin)/volume,TMath::Power(scale*histo->GetBinError(bin)/volume,2)) ;
          }
        }
      } else {
        int bin = histo->GetBin(ix+offset);
        double volume = useIf(histo->GetXaxis()->GetBinWidth(ix+offset),correctDensity);
        dh->add(obs,scale*histo->GetBinContent(bin)/volume,TMath::Power(scale*histo->GetBinError(bin)/volume,2)) ;     
      }
    }
    
    dh->removeSelfFromDir();

    return dh;
  }



  std::vector<RooRealVar*> createGammas(const TH1* histo, bool includeUnderflowOverflow, double uncThreshold){
    std::vector<RooRealVar*> gammas;
    int offset = !includeUnderflowOverflow;
    Int_t ix(0),iy(0),iz(0) ;
    for (ix=0 ; ix < histo->GetNbinsX()+2*includeUnderflowOverflow ; ix++) {
      if (dim(histo)>1) {
        for (iy=0 ; iy < histo->GetNbinsY()+2*includeUnderflowOverflow ; iy++) {
          if (dim(histo)>2) {
            for (iz=0 ; iz < histo->GetNbinsZ()+2*includeUnderflowOverflow ; iz++) {
              double val = histo->GetBinContent(ix+offset,iy+offset,iz+offset);
              double err = histo->GetBinError(ix+offset,iy+offset,iz+offset);
              if(val > 0 && err/val>uncThreshold){
                TString name = TString::Format("gamma_stat_%s_%d",histo->GetName(),(int)gammas.size());
                RooRealVar* g = new RooRealVar(name,name,1.);
                g->setError(err/val);
                gammas.push_back(g);                
              } else {
                gammas.push_back(0);
              }
            }
          } else {
            double val = histo->GetBinContent(ix+offset,iy+offset);
            double err = histo->GetBinError(ix+offset,iy+offset);            
            if(val > 0 && err/val>uncThreshold){
              TString name = TString::Format("gamma_stat_%s_%d",histo->GetName(),(int)gammas.size());
              RooRealVar* g = new RooRealVar(name,name,1.);
              g->setError(err/val);
              gammas.push_back(g);                
            } else {
              gammas.push_back(0);
            }
          }
        }
      } else {
        double val = histo->GetBinContent(ix+offset);
        double err = histo->GetBinError(ix+offset);              
        if(val > 0 && err/val>uncThreshold){
          TString name = TString::Format("gamma_stat_%s_%d",histo->GetName(),(int)gammas.size());
          RooRealVar* g = new RooRealVar(name,name,1.);
          g->setError(err/val);
          gammas.push_back(g);                
        } else {
          gammas.push_back(0);
        }
      }
    }
    bool allzero=true;
    for(auto g:gammas){
      if(g) allzero=false;
    }
    if(allzero) gammas.clear();
    return gammas;
  }

  std::vector<RooRealVar*> createGammas(const RooDataHist* dh, const RooArgList& obs, double uncThreshold){
    std::vector<RooRealVar*> gammas;

    // Define x,y,z as 1st, 2nd and 3rd observable
    RooAbsArg* xvar = obs.at(0);
    RooAbsArg* yvar = obs.at(1);
    RooAbsArg* zvar = obs.at(2);

    // Transfer contents
    Int_t ix(0),iy(0),iz(0) ;
    for (ix=0 ; ix < ::nBins(xvar) ; ix++) {
      ::setBin(xvar,ix) ;
      if (yvar) {
        for (iy=0 ; iy < ::nBins(yvar) ; iy++) {
          ::setBin(yvar,iy) ;
          if (zvar) {
            for (iz=0 ; iz < ::nBins(zvar) ; iz++) {
              ::setBin(zvar,iz) ;
              dh->get(obs);
              double val = dh->weight();
              double err = sqrt(dh->weightSquared());
              if(val > 0 && err/val>uncThreshold){
                TString name = TString::Format("gamma_stat_%s_%d",dh->GetName(),(int)gammas.size());
                RooRealVar* g = new RooRealVar(name,name,1.);
                g->setError(err/val);
                gammas.push_back(g);                
              } else {
                gammas.push_back(0);
              }
            }
          } else {
            dh->get(obs);
            double val = dh->weight();
            double err = sqrt(dh->weightSquared());
            if(val > 0 && err/val>uncThreshold){
              TString name = TString::Format("gamma_stat_%s_%d",dh->GetName(),(int)gammas.size());
              RooRealVar* g = new RooRealVar(name,name,1.);
              g->setError(err/val);
              gammas.push_back(g);                
            } else {
              gammas.push_back(0);
            }
          }
        }
      } else {
        dh->get(obs);
        double val = dh->weight();
        double err = sqrt(dh->weightSquared());
        if(val > 0 && err/val>uncThreshold){
          TString name = TString::Format("gamma_stat_%s_%d",dh->GetName(),(int)gammas.size());
          RooRealVar* g = new RooRealVar(name,name,1.);
          g->setError(err/val);
          gammas.push_back(g);                
        } else {
          gammas.push_back(0);
        }
      }
    }
    bool allzero=true;
    for(auto g:gammas){
      if(g) allzero=false;
    }
    if(allzero) gammas.clear();    
    return gammas;
  }  
  
  
  RooFitHist::RooFitHist(const TH1* hist, const std::vector<RooAbsArg*>& obs, bool includeUnderflowOverflow, double uncThreshold, bool correctDensity) : RooFitHist(hist,argList(obs),includeUnderflowOverflow,uncThreshold,correctDensity) {}
  RooFitHist::RooFitHist(const TH1* hist, const RooArgList& obslist, bool includeUnderflowOverflow, double uncThreshold, bool correctDensity) : RooFitHist(convertTH1(hist,obslist,includeUnderflowOverflow,correctDensity),obslist,uncThreshold){}


  RooFitHist::RooFitHist(RooHistFunc* hf, const RooArgList& obslist, double uncThreshold){
    RooFIter itr(obslist.fwdIterator());
    RooAbsArg* arg = NULL;
    while((arg = (RooAbsArg*)itr.next())){
      if(!arg) continue;
      if(!hf->dependsOn(*arg)) continue;
      this->_obs.push_back(arg);
    }
    if(uncThreshold >= 0){
      this->_func = this->setupErrors(hf,&(hf->dataHist()),uncThreshold);
    } else {
      this->_func = hf;
    }
  }
  
  RooFitHist::RooFitHist(RooDataHist* dh, const RooArgList& obslist, double uncThreshold){
    TString name(dh->GetName());
    TString title(dh->GetTitle());
    RooHistFunc* hf = new RooHistFunc(TString::Format("%s_Values",name.Data()),title,obslist,*dh);
    RooFIter itr(obslist.fwdIterator());
    RooAbsArg* arg = NULL;
    while((arg = (RooAbsArg*)itr.next())){
      if(!arg) continue;
      if(!hf->dependsOn(*arg)) continue;
      this->_obs.push_back(arg);
    }    
    if(uncThreshold >= 0){
      this->_func = this->setupErrors(hf,dh,uncThreshold);
    } else {
      this->_func = hf;
    }
  }

  RooAbsReal* makeParamHistFunc(const char* name, const char* title, const RooArgList& obslist, const std::vector<RooRealVar*>& gamma){
    RooArgList gammas;
    RooRealVar* const_g = NULL;
    for(auto g:gamma){
      if(g){
        g->setConstant(false);
        gammas.add(*g,true);
      } else {
        if(!const_g){
          const_g = new RooRealVar("gamma_stat_const","gamma_stat_const",1.);
          const_g->setConstant(true);
        }
        gammas.add(*const_g,true);
      }
    }
    ParamHistFunc* phf = new ParamHistFunc(TString::Format("%s_mcstat",name),title,obslist,gammas);
    phf->recursiveRedirectServers(obslist);

    return phf;
  }
  
  RooAbsReal* RooFitHist::setupErrors(const RooHistFunc* hf, const RooDataHist* dh, double uncThreshold){
    RooArgList obslist;
    for(auto obs:this->_obs) obslist.add(*obs,true);
    this->_gamma = createGammas(dh,obslist,uncThreshold);
    RooArgList components;
    components.add(*hf);
    if(this->_gamma.size()>0){
      RooAbsReal* phf =  makeParamHistFunc(hf->GetName(),hf->GetTitle(),obslist,this->_gamma);
      components.add(*phf);
    }
    RooProduct* prod = new RooProduct(TString::Format("%s_x_Uncertainties",hf->GetName()),hf->GetTitle(),components);
    return prod;
  }

  void RooFitHist::replace(const RooAbsCollection& newlist){
    ::replace(this->_func,newlist);
    for(size_t i=0; i<this->_obs.size(); ++i){
      ::replace(this->_obs[i],newlist);
    }
    for(size_t i=0; i<this->_gamma.size(); ++i){
      ::replace(this->_gamma[i],newlist);
    }
  }

  double RooFitHist::error() const {
    size_t i = this->bin();
    if(this->_gamma.size() <= (size_t)i) return sqrt(this->value());
    if(!this->_gamma[i]) return 0;
    double cval = this->_gamma[i]->getVal();
    double nomval = this->_func->getVal();
    this->_gamma[i]->setVal(cval+fabs(this->_gamma[i]->getError()));
    double eup = this->_func->getVal() - nomval;
    this->_gamma[i]->setVal(cval-fabs(this->_gamma[i]->getError()));
    double edn = this->_func->getVal() - nomval;
    double err = (fabs(eup)+fabs(edn))/2;
    this->_gamma[i]->setVal(cval);
    return err;
  }
  double RooFitHist::value() const {
    if(this->_func->InheritsFrom(RooAbsPdf::Class())){
      const RooAbsPdf* pdf = static_cast<const RooAbsPdf*>(this->_func);
      if(pdf->extendMode() == RooAbsPdf::CanBeExtended || pdf->extendMode() == RooAbsPdf::MustBeExtended){
        return pdf->getVal()*pdf->expectedEvents(0);
      }
    }
    return this->_func->getVal();
  }
  bool RooFitHist::weighted() const {
    return (this->_gamma.size() > 0);
  }
  void RooFitHist::printHist() const {
    printHistogram(this);
  }
  bool RooFitHist::checkValidity() const {
    for(auto obs:_obs){
      if(!this->_func->dependsOn(*obs)){
        throw std::runtime_error("observable mismatch!");
      }
    }
    for(auto np:_gamma){
      if(!np) continue;
      if(!this->_func->dependsOn(*np)){
        throw std::runtime_error("nuisance parameter mismatch!");
      }
    }
    return true;
  }

  RooHistFunc* makeHistFunc(const char* name, const TH1* histo, const RooArgList& obs, bool includeUnderflowOverflow, bool correctDensity){
    RooDataHist* dh = convertTH1(histo,obs,includeUnderflowOverflow,correctDensity);
    return new RooHistFunc(name,histo->GetTitle(),obs,*dh);
  }
  RooHistFunc* makeHistFunc(const TH1* histo, const RooArgList& obs, bool includeUnderflowOverflow, bool correctDensity){
    TString name(TString::Format("%s_histfunc",histo->GetName()));
    return makeHistFunc(name.Data(),histo,obs,includeUnderflowOverflow,correctDensity);
  }
  RooHistFunc* makeHistFunc(RooDataHist* dhist, const std::vector<RooAbsArg*>& obs){
    if(!dhist) return NULL;
    RooArgSet vars;
    for(size_t i=0; i<obs.size(); ++i){
      if(dhist->get()->find(*obs[i]))
        vars.add(*obs[i],true);
    }
    return new RooHistFunc(TString::Format("%s_func",dhist->GetName()),dhist->GetTitle(),vars,*dhist);
  }
  RooHistFunc* makeHistFunc(RooDataHist* dhist, const RooArgList& obs){
    if(!dhist) return NULL;
    return new RooHistFunc(TString::Format("%s_func",dhist->GetName()),dhist->GetTitle(),obs,*dhist);
  }
  RooAbsPdf* makeHistPdf(const char* name, const TH1* histo, const RooArgList& obs, bool includeUnderflowOverflow, bool correctDensity){
    double integral = getIntegral(histo,includeUnderflowOverflow,correctDensity);
    RooDataHist* dh = convertTH1(histo,obs,includeUnderflowOverflow,correctDensity,1./integral);
    RooHistPdf* hf = new RooHistPdf(TString::Format("%s_shape",name),histo->GetTitle(),obs,*dh);
    hf->setUnitNorm(true);
    RooRealVar* norm = new RooRealVar(TString::Format("%s_norm",dh->GetName()),dh->GetTitle(),integral);
    norm->setConstant(true);
    return new RooExtendPdf(name,histo->GetTitle(),*hf,*norm);
  }
  RooAbsPdf* makeHistPdf(const TH1* histo, const RooArgList& obs, bool includeUnderflowOverflow, bool correctDensity){
    return makeHistPdf(histo->GetName(),histo,obs,includeUnderflowOverflow,correctDensity);
  }
  RooAbsPdf* makeHistPdf(RooDataHist* dhist, const std::vector<RooAbsArg*>& obs){
    if(!dhist) return NULL;
    RooArgSet vars;
    for(size_t i=0; i<obs.size(); ++i){
      if(dhist->get()->find(*obs[i]))
        vars.add(*obs[i],true);
    }
    RooHistPdf* hf = new RooHistPdf(TString::Format("%s_shape",dhist->GetName()),dhist->GetTitle(),vars,*dhist);
    RooRealVar* norm = new RooRealVar(TString::Format("%s_norm",dhist->GetName()),dhist->GetTitle(),dhist->sumEntries());
    norm->setConstant(true);
    return new RooExtendPdf(TString::Format("%s_pdf",dhist->GetName()),dhist->GetTitle(),*hf,*norm);
  }

  RooFitHist::RooFitHist(RooDataHist* d, const std::vector<RooAbsArg*>& obs) : _func(0), _obs(obs) {
    this->_func = makeHistFunc(d,_obs);
  }
  RooFitHist::RooFitHist(RooDataHist* d, RooAbsArg* obs) : _func(0) {
    this->_obs.push_back(obs);
    this->_func = makeHistFunc(d,_obs);    
  }
  RooFitHist::RooFitHist(RooDataHist* d, RooAbsArg* obs1, RooAbsArg* obs2) : _func(0) { 
    this->_obs.push_back(obs1); this->_obs.push_back(obs2);
    this->_func = makeHistFunc(d,_obs);    
  } 


  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const char* name, const char* /*title*/, const Variable<RooUnfolding::RooFitHist>& x, const Variable<RooUnfolding::RooFitHist>& y){
    RooArgSet vars(*x._var,*y._var, name);
    std::runtime_error("createHist 1 called");
    return NULL;
  }
  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const char* /*name*/, const char* /*title*/, const std::vector<Variable<RooUnfolding::RooFitHist>>& /*x*/){
    std::runtime_error("createHist 2 called");
    return NULL;
  }
  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const TMatrixD& /*m*/, const char* /*name*/, const char* /*title*/, const Variable<RooUnfolding::RooFitHist>& /*x*/, const Variable<RooUnfolding::RooFitHist>& /*y*/){  
    // Sets the bin content of the histogram as that element of the input vector
    std::runtime_error("createHist 3 called");        
    return NULL;
  }
  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const TMatrixD& /*m*/, const TMatrixD& /*me*/, const char* /*name*/, const char* /*title*/, const Variable<RooUnfolding::RooFitHist>& /*x*/, const Variable<RooUnfolding::RooFitHist>& /*y*/){  
    // Sets the bin content of the histogram as that element of the input vector
    std::runtime_error("createHist 4 called");            
    return NULL;
  }

  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const TVectorD& /*v*/, const char* /*name*/, const char* /*title*/, const std::vector<Variable<RooUnfolding::RooFitHist>>& /*x*/, bool /*overflow*/){
    std::runtime_error("createHist 5 called");
    return 0;
  }
  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const TVectorD& /*v*/, const TVectorD& /*ve*/, const char* /*name*/, const char* /*title*/, const std::vector<Variable<RooUnfolding::RooFitHist>>& /*x*/, bool /*overflow*/){
    std::runtime_error("createHist 6 called");   
    return 0;
  }
  template<> RooUnfolding::RooFitHist* createHist<RooUnfolding::RooFitHist>(const TVectorD& /*v*/, const TVectorD& /*verr*/, const char* /*name*/, const char* /*title*/, const RooUnfolding::RooFitHist* /*origHist*/, bool /*overflow*/){
    std::runtime_error("createHist 7 called");
    return 0;
  }
  template<> TVectorD subtract<RooUnfolding::RooFitHist,TVectorD>(const TVectorD& orig, const RooUnfolding::RooFitHist* hist, bool overflow) {
    TVectorD res(orig);
    size_t n = nBins(hist);
    for(size_t i=0; i<n; ++i){
      res[i] -= binContent(hist,i,overflow);
    }
    return res;
  }

  template<> void h2v<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, TVectorD& v, bool overflow, bool correctDensity){
    size_t n = nBins(h);
    v.ResizeTo(n);
    for(size_t i=0; i<n; ++i){
      v[i] = binContent(h,i,overflow) * (correctDensity ? binVolume(h,i,overflow) : 1);
    }
  }
  template<> void h2ve<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, TVectorD& v, bool overflow, bool correctDensity){
    size_t n = nBins(h);
    v.ResizeTo(n);
    for(size_t i=0; i<n; ++i){
      v[i] = binError(h,i,overflow) * (correctDensity ? binVolume(h,i,overflow) : 1);      
    }    
  }    

  template<> TVectorD h2v<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, bool overflow, bool correctDensity){
    TVectorD v;
    h2v(h,v,overflow,correctDensity);
    return v;    
  }
  template<> TVectorD h2ve<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, bool overflow, bool correctDensity){
    TVectorD v;
    h2ve(h,v,overflow,correctDensity);
    return v;    
  }
  template<> void h2mNorm<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, TMatrixD& m, const RooUnfolding::RooFitHist* norm, bool overflow, bool correctDensity){
    // sets Matrix to values of bins in a 2D input histogram
    size_t nx = nBins(h,X);
    size_t ny = nBins(h,Y);
    m.ResizeTo(nx,ny);
    bool needSanitization = !norm;
    for(size_t j=0; j<ny; ++j){
      double fac = 1.;
      if (norm){
        int b = bin(norm,j,overflow);
        fac= binContent(norm,b,overflow) * (correctDensity ? binVolume(norm,b,overflow) : 1);
        if (fac != 0.0){
          fac= 1.0/fac;
        } else {
          needSanitization = true;
        }          
      }
      for(size_t i=0; i<nx; ++i){
        int b = bin(h,i,j,overflow);
        m(i,j) = binContent(h,b,overflow) * fac * (correctDensity ? binVolume(h,b,overflow) : 1);
      }
    }
    if(needSanitization) sanitize(m);
  }
  template<> void h2meNorm<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, TMatrixD& m, const RooUnfolding::RooFitHist* norm, bool overflow, bool correctDensity){
    // sets Matrix to errors of bins in a 2D input histogram
    size_t nx = nBins(h,X);
    size_t ny = nBins(h,Y);
    m.ResizeTo(nx,ny);
    for(size_t j=0; j<ny; ++j){
      double fac = 1.;
      if (norm){
        int b = bin(norm,j,overflow);
        fac= binContent(norm,b,overflow) * (correctDensity ? binVolume(norm,b,overflow) : 1);
        if (fac != 0.0) fac= 1.0/fac;
      }
      for(size_t i=0; i<nx; ++i){
        int b = bin(h,i,j,overflow);
        m(i,j) = binError(h,b,overflow) * fac * (correctDensity ? binVolume(h,b,overflow) : 1);
      }
    }
  }
  template<> TMatrixD h2mNorm<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, const RooUnfolding::RooFitHist* norm, bool overflow, bool correctDensity){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2mNorm(h,m,norm, overflow,correctDensity);
    return m;
  }
  template<> TMatrixD h2meNorm<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, const RooUnfolding::RooFitHist* norm, bool overflow, bool correctDensity){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2meNorm(h,m,norm, overflow,correctDensity);
    return m;
  }
  template<> void h2m  (const RooUnfolding::RooFitHist* h, TMatrixD& m, bool overflow, bool correctDensity) { h2mNorm (h,m,(const RooUnfolding::RooFitHist*)NULL,overflow,correctDensity);}
  template<> void h2me  (const RooUnfolding::RooFitHist* h, TMatrixD& m, bool overflow, bool correctDensity){ h2meNorm(h,m,(const RooUnfolding::RooFitHist*)NULL,overflow,correctDensity);}
  template<> TMatrixD h2m<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, bool overflow, bool correctDensity){
    // Returns Matrix of values of bins in a 2D input histogram
    TMatrixD m;
    h2m(h,m,overflow,correctDensity);
    return m;
  }
  template<> TMatrixD h2me<RooUnfolding::RooFitHist>  (const RooUnfolding::RooFitHist* h, bool overflow, bool correctDensity){
    // Returns Matrix of errors of bins in a 2D input histogram
    TMatrixD m;
    h2me(h,m,overflow,correctDensity);
    return m;
  }

  const std::vector<RooRealVar*>& RooUnfolding::RooFitHist::nps() const {
    return this->_gamma;
  }

  template<> RooUnfolding::RooFitHist* clone<RooUnfolding::RooFitHist>(RooUnfolding::RooFitHist const* h){
    if(!h) return NULL;
    RooAbsReal* rr = (RooAbsReal*)(h->func()->cloneTree());
    std::vector<RooAbsArg*> obs;
    for(int i=0; i<dim(h); ++i){
      RooAbsArg* rrv = findLeafServer<RooAbsArg>(rr,h->obs(i)->GetName());
      obs.push_back(rrv);
    }
    std::vector<RooRealVar*> nps;
    for(auto np:h->nps()){
      RooRealVar* new_np = NULL;
      if(np){
        new_np = findLeafServer<RooRealVar>(rr,np->GetName());
      } 
      nps.push_back(new_np);
    }    
    if(h->dim() != obs.size()) throw std::runtime_error("dimensionality mismatch!");
    if(nps.size() > 0 && nBins(h) != (int)nps.size()) throw std::runtime_error("bin number mismatch!");
    RooFitHist* newh = new RooFitHist(rr,obs,nps);
    newh->checkValidity();
    return newh;
  }
}  

template RooUnfolding::RooFitHist* RooUnfolding::createHist<RooUnfolding::RooFitHist>(TVectorT<double> const&, char const*, char const*, RooUnfolding::Variable<RooUnfolding::RooFitHist> const&, bool);
template std::vector<RooUnfolding::Variable<RooUnfolding::RooFitHist> > RooUnfolding::vars<RooUnfolding::RooFitHist>(RooUnfolding::RooFitHist const*); 
template void RooUnfolding::printTable<RooUnfolding::RooFitHist>(std::ostream&, RooUnfolding::RooFitHist const*, RooUnfolding::RooFitHist const*, RooUnfolding::RooFitHist const*, RooUnfolding::RooFitHist const*, RooUnfolding::RooFitHist const*, bool, RooUnfolding::ErrorTreatment, double);

#include "RooWorkspace.h"

void RooUnfolding::importToWorkspace(RooWorkspace* ws, RooAbsReal* object){
  // insert an object into a workspace (wrapper for RooWorkspace::import)
  if(!ws) return;
  if(!object) return;
  ws->import(*object,RooFit::RecycleConflictNodes());
}

void RooUnfolding::importToWorkspace(RooWorkspace* ws, RooAbsData* object){
  // insert an object into a workspace (wrapper for RooWorkspace::import)
  if(!ws) return;
  if(!object) return;
  ws->import(*object);
}

void RooUnfolding::setGammaUncertainties(RooWorkspace* ws){
  RooArgSet pdfs (ws->allPdfs());
  RooFIter itr(pdfs.fwdIterator());
  TPRegexp re("gamma_stat_.*");    
  RooAbsArg* obj = NULL;;
  while((obj = itr.next())){
    RooPoisson* pois = dynamic_cast<RooPoisson*>(obj);
    if(!pois) continue;
    if(!re.Match(pois->GetName())) continue;
    RooRealVar* gamma = getGammaParameter(pois);
    double val = getMean(pois);
    gamma->setError(sqrt(val)/val);
  }
}

const RooArgSet* RooUnfolding::getObservables(const RooHistFunc* f){
  return &(f->*RooHistFuncHackResult<RooHistFuncObs>::ptr);
}

void RooUnfolding::printClients(const RooAbsArg* obj){
  TIterator* itr = obj->clientIterator();
  TObject* x;
  std::cout << obj << " " << obj->ClassName() << " " << obj->GetName() << " has the following clients" << std::endl;
  while((x = itr->Next())){
    std::cout << "  " << x << " " << x->ClassName() << " " << x->GetName() << std::endl;
  }
}
void RooUnfolding::printServers(const RooAbsArg* obj){
  TIterator* itr = obj->serverIterator();
  TObject* x;
  std::cout << obj << " " << obj->ClassName() << " " << obj->GetName() << " has the following servers" << std::endl;
  while((x = itr->Next())){
    std::cout << "  " << x << " " << x->ClassName() << " " << x->GetName() << std::endl;
  }
}  

 template<class ObjT,class ListT>std::vector<ObjT*> RooUnfolding::matchingObjects(const ListT* c, const char* pattern){
  TPRegexp re(pattern);
  std::vector<ObjT*> retval;
  RooFIter itr(c->fwdIterator());
  ObjT* arg = NULL;
  while((arg = (ObjT*)itr.next())){
    if(!arg) continue;
    if(re.Match(arg->GetName())){
      retval.push_back(arg);
    }
  }
  return retval;
}
std::vector<RooAbsReal*> RooUnfolding::matchingObjects(const RooAbsCollection* c, const char* pattern){
  return matchingObjects<RooAbsReal,RooAbsCollection>(c,pattern);
}


RooArgSet RooUnfolding::allVars(RooWorkspace* ws, const char* pattern){
  TPRegexp re(pattern);
  RooArgSet allVars(ws->allVars());
  RooFIter itr(allVars.fwdIterator());
  RooAbsArg* arg = NULL;
  RooArgSet retval;
  while((arg = itr.next())){
    if(!arg) continue;
    if(re.Match(arg->GetName())){
      retval.add(*arg,true);
    }
  }
  return retval;
}

RooUnfolding::RooFitHist* RooUnfolding::RooFitHist::asimovClone(bool correctDensity) const {   
  // Define x,y,z as 1st, 2nd and 3rd observable
  RooAbsArg* xvar = _obs.at(0);
  RooAbsArg* yvar = _obs.size() > 1 ? _obs.at(1) : (RooAbsArg*)0;
  RooAbsArg* zvar = _obs.size() > 2 ? _obs.at(2) : (RooAbsArg*)0;
  
  RooArgSet args;
  RooArgList arglist;    
  for(auto v:this->_obs){
    args.add(*v);
    arglist.add(*v);      
  }
  TString name = TString::Format("%s_asimov",this->name());
  
  RooDataHist* dh = new RooDataHist(name,this->title(),args);
  
  // Transfer contents
  for (int ix=0 ; ix < ::nBins(xvar) ; ix++) {
    ::setBin(xvar,ix) ;
    if (yvar) {
      for (int iy=0 ; iy < ::nBins(yvar) ; iy++) {
        ::setBin(yvar,iy) ;
        if (zvar) {
          for (int iz=0 ; iz < ::nBins(zvar) ; iz++) {
            ::setBin(zvar,iz) ;
            double volume = ::useIf(dh->binVolume(),correctDensity);
            dh->add(args,this->value(),sqrt(volume*this->value())/volume) ;
          }
        } else {
          double volume = ::useIf(dh->binVolume(),correctDensity);
          dh->add(args,this->value(),sqrt(volume*this->value())/volume) ;
        }
      }
    } else {
      double volume = ::useIf(dh->binVolume(),correctDensity);
      dh->add(args,this->value(),sqrt(volume*this->value())/volume) ;
    }
  }
  
  dh->removeSelfFromDir();
  
  return new RooFitHist(dh,arglist,0);
}

RooUnfolding::RooFitHist* RooUnfolding::RooFitHist::asimov1DClone(bool /*correctDensity*/, TVectorD& val, TVectorD& err) const {   
  // Define x,y,z as 1st, 2nd and 3rd observable
  RooAbsArg* xvar = _obs.at(0);

  RooArgSet args;
  RooArgList arglist;    
  for(auto v:this->_obs){
    args.add(*v);
    arglist.add(*v);      
  }
  TString name = TString::Format("%s_asimov",this->name());
  
  RooDataHist* dh = new RooDataHist(name,this->title(),args);
  
  // Transfer contents
  for (int ix=0 ; ix < ::nBins(xvar) ; ix++) {
    ::setBin(xvar,ix) ;
    dh->add(args,val(ix),err(ix));
  }
  
  dh->removeSelfFromDir();
  
  return new RooFitHist(dh,arglist,0);
}
namespace RooUnfolding {
  template<> RooUnfolding::RooFitHist* asimovClone(const RooUnfolding::RooFitHist* hist, bool correctDensity){
    return hist->asimovClone(correctDensity);
  }
  template<> RooUnfolding::RooFitHist* asimov1DClone(const RooUnfolding::RooFitHist* hist, bool correctDensity, TVectorD& val, TVectorD& err){
    return hist->asimov1DClone(correctDensity, val, err);
  }
  std::vector<Variable<TH1>> convertTH1(const std::vector<Variable<RooUnfolding::RooFitHist> >& vars){
    std::vector<Variable<TH1> > outvars;
    for(auto var:vars){
      auto v = var._var;
      outvars.push_back(RooUnfolding::Variable<TH1>(::nBins(v),::min(v),::max(v),v->GetName()));
    }
    return outvars;
  }
  std::vector<Variable<TH2>> convertTH2(const std::vector<Variable<RooUnfolding::RooFitHist> >& vars){
    std::vector<Variable<TH2> > outvars;
    for(auto var:vars){
      auto v = var._var;
      outvars.push_back(RooUnfolding::Variable<TH2>(::nBins(v),::min(v),::max(v),v->GetName()));
    }
    return outvars;
  }
  TH1* convertTH1(const TVectorD& values, const TVectorD& errors, const RooUnfolding::RooFitHist* hist){
    return RooUnfolding::createHist<TH1>(values,errors,hist->GetName(),hist->GetTitle(),RooUnfolding::convertTH1(RooUnfolding::vars(hist)));
  }
  TH1* convertTH1(const TVectorD& values, const RooUnfolding::RooFitHist* hist){
    return RooUnfolding::createHist<TH1>(values,hist->GetName(),hist->GetTitle(),RooUnfolding::convertTH1(RooUnfolding::vars(hist)));
  }
  TH2* convertTH2(const TMatrixD& values, const TMatrixD& errors, const RooUnfolding::RooFitHist* hist){
    std::vector<Variable<TH2> > outvars = RooUnfolding::convertTH2(RooUnfolding::vars(hist));
    return RooUnfolding::createHist<TH2>(values,errors,hist->GetName(),hist->GetTitle(),outvars.at(0),outvars.at(1));
  }
}


ClassImp(RooUnfolding::RooFitHist)

#endif
