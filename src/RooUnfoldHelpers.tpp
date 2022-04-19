// -*- mode: c++ -*-

#pragma once

namespace RooUnfolding {
  template<class Hist> Hist* createHist(const TVectorD& vec, const char* name, const char* title, const Variable<Hist>& x, bool overflow) { return createHist<Hist>(vec,name,title,std::vector<Variable<Hist>>{x},overflow); }
  template<class Hist> Hist* createHist(const TVectorD& vec, const TVectorD& errvec, const char* name, const char* title, const Variable<Hist>& x, bool overflow) { return createHist<Hist>(vec,errvec,name,title,std::vector<Variable<Hist>>{x},overflow); }
  template<class Hist> std::vector<RooUnfolding::Variable<Hist>> vars(const Hist* h){
    int d = dim(h);
    std::vector<Variable<Hist> > v;
    v.push_back(var(h,X));
    if(d>1) v.push_back(var(h,Y));
    if(d>2) v.push_back(var(h,Z));
    return v;
  }


  template<class Hist> Hist* histNoOverflow(const Hist* hist, bool overflow){
    return createHist<Hist>(h2v<Hist>(hist,overflow),h2ve<Hist>(hist,overflow),name(hist),title(hist),vars(hist),false);
  }
  
  template<class Hist> void printTable (std::ostream& o, const Hist* hTrainTrue, const Hist* hTrain,
                   const Hist* hTrue, const Hist* hMeas, const Hist* hReco,
                   Bool_t overflow,
                   ErrorTreatment withError, Double_t chi_squ)
  {
    // Prints entries from truth, measured, and reconstructed data for each bin.
    if(!hTrainTrue) throw std::runtime_error("unable to print: no trainTrue given!");
    if(!hTrain) throw std::runtime_error("unable to print: no train given!");
    if(!hMeas) throw std::runtime_error("unable to print: no measured given!");
    if(!hReco) throw std::runtime_error("unable to print: no reco given!");    
    if (withError==kDefault) withError= sumW2N(hReco) ? kErrors : kNoError;
    Int_t d= dim(hReco);
    int ntxb= nBins(hReco,X)+2*overflow, ntyb= nBins(hReco,Y)+2*overflow;
    if (dim(hMeas) != d || nBins(hMeas,X)+2*overflow != ntxb || nBins(hMeas,Y)+2*overflow != ntyb) d= 1;
    printTable(o,d,
               ntxb,ntyb,
               h2v(hTrainTrue,overflow),
               h2v(hTrain,overflow),
               h2v(hTrue,overflow),
               h2v(hMeas,overflow),
               h2v(hReco,overflow),
               withError,
               h2ve(hTrue,overflow),
               h2ve(hReco,overflow),               
               chi_squ,overflow);
  }

}


