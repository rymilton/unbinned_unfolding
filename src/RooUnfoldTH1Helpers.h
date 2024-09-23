#ifndef ROOUNFOLDHELPERS_TH1_HH
#define ROOUNFOLDHELPERS_TH1_HH

#include "RooUnfoldHelpers.h"
class TAxis;
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace RooUnfolding {
  const TAxis* getAxis(const TH1* h, RooUnfolding::Dimension d);
  void print(int _nBins, double _min, double _max, const std::vector<double>& _bounds);
 
  template<> struct Variable<TH1> {
    int _nBins;
    double _min;
    double _max;
    std::vector<double> _bounds;
    bool irregular() const { return _bounds.size() > 0; }
    void print() { RooUnfolding::print(_nBins, _min, _max, _bounds); }
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
    Variable(int nBins, const double* bounds,const char*) : _nBins(nBins),_min(bounds[0]),_max(bounds[nBins]) { for(size_t i=0; i<=nBins; ++i){ _bounds.push_back(bounds[i]); } }
    Variable(const std::vector<double>& bounds,const char* s) : Variable(bounds.size()-1,&(bounds[0]),s) {};
  };
  template<> struct Variable<TH2> {
    int _nBins;
    double _min;
    double _max;
    std::vector<double> _bounds;    
    bool irregular() const { return _bounds.size() > 0; }
    void print() { RooUnfolding::print(_nBins, _min, _max, _bounds); }
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
    Variable(int nBins, const double* bounds,const char*) : _nBins(nBins),_min(bounds[0]),_max(bounds[nBins]) { for(size_t i=0; i<=nBins; ++i){ _bounds.push_back(bounds[i]); } }    
    Variable(const std::vector<double>& bounds,const char* s) : Variable(bounds.size()-1,&(bounds[0]),s) {};
  };
  template<> struct Variable<TH3> {
    int _nBins;
    double _min;
    double _max;
    std::vector<double> _bounds;
    bool irregular() const { return _bounds.size() > 0; }
    void print() { RooUnfolding::print(_nBins, _min, _max, _bounds); }
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
    Variable(int nBins, const double* bounds,const char*) : _nBins(nBins),_min(bounds[0]),_max(bounds[nBins]) { for(size_t i=0; i<=nBins; ++i){ _bounds.push_back(bounds[i]); } }
    Variable(const std::vector<double>& bounds,const char* s) : Variable(bounds.size()-1,&(bounds[0]),s) {};
  };

  Bool_t resizeAxis (TAxis* ax, Int_t nx);
  TH1* resize (TH1* h, Int_t nx, Int_t ny=1, Int_t nz=1);
  TH1* convertTH1(const TVectorD& values, const TVectorD& errors, const TH1* hist);
  TH1* convertTH1(const TVectorD& values, const TH1* hist);
  template<class Hist> bool irregular(const Hist* hist, RooUnfolding::Dimension d);
  template<class Hist> std::vector<double> binning(const Hist* hist, RooUnfolding::Dimension d);
}

#endif
