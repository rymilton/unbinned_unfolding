//=====================================================================-*-C++-*-
//! \class Omnifold
//! \brief Unbinned & binned versions of Omnifold (ML unfolding) using decision trees
//! \author Ryan Milton <rmilt003@ucr.edu>
//==============================================================================

#ifndef OMNIFOLD_H_
#define OMNIFOLD_H_

#include "RooUnfoldResponse.h"
#include "TH1.h"
#include <tuple>
#include <TVector.h>
#include <TVectorD.h>

class Omnifold {
public:
    Omnifold();
    Omnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations);
    ~Omnifold();
    TH1D* BinnedOmnifold();
    std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> UnbinnedOmnifold(TVectorD MC_entries, TVectorD sim_entries, TVectorD measured_entries, TVector pass_reco, TVector pass_truth, Int_t num);
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
