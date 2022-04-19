#ifndef ROOUNFOLDHELPERS_TH1_HH
#define ROOUNFOLDHELPERS_TH1_HH

#include "RooUnfoldHelpers.h"
class TAxis;
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

namespace RooUnfolding {
  const TAxis* getAxis(const TH1* h, RooUnfolding::Dimension d);
  
  template<> struct Variable<TH1> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
  };
  template<> struct Variable<TH2> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
  };
  template<> struct Variable<TH3> {
    int _nBins;
    double _min;
    double _max;
    Variable(int nBins,double min,double max,const char*) : _nBins(nBins),_min(min),_max(max){};
  };

  TH1* resize (TH1* h, Int_t nx, Int_t ny=1, Int_t nz=1);
  TH1* convertTH1(const TVectorD& values, const TVectorD& errors, const TH1* hist);
  TH1* convertTH1(const TVectorD& values, const TH1* hist);  
}

#endif
