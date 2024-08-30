//=====================================================================-*-C++-*-
//! \class Omnifold
//! \brief Unbinned & binned versions of Omnifold (ML unfolding) using decision trees
//! \author Ryan Milton <rmilt003@ucr.edu>
//==============================================================================

#ifndef OMNIFOLD_H_
#define OMNIFOLD_H_

#include "RooUnfoldResponse.h"
#include "TH1.h"

class Omnifold {
public:
    Omnifold();
    Omnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations);
    ~Omnifold();
    TH1D* BinnedOmnifold();
    void SetIterations(Int_t nIter) {_nIter = nIter;}
    Int_t GetIterations() const {return _nIter;}
    void SetMeasuredHist(TH1* measured_hist) {_measuredHist = measured_hist;}
    TH1* GetMeasuredHist() const {return _measuredHist;}
    void SetResponseMatrix(RooUnfoldResponse response) {_response = response;}
    RooUnfoldResponse GetResponseMatrix() const {return _response;}
    void EfficiencyCorrections(TH1* hist);

private:
    Int_t _nIter;
    RooUnfoldResponse _response;
    TH1* _measuredHist;
};
#endif
