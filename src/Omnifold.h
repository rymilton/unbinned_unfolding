//=====================================================================-*-C++-*-
//! \class Omnifold
//! \brief Unbinned & binned versions of Omnifold (ML unfolding) using decision trees
//! \author Ryan Milton <rmilt003@ucr.edu>
//==============================================================================

#ifndef OMNIFOLD_H_
#define OMNIFOLD_H_

#include <tuple>
#include "TH1.h"
#include "TVector.h"
#include "TVectorD.h"
#include "RooUnfoldResponse.h"
#include "TObjArray.h"

class Omnifold {
public:
    Omnifold();
    ~Omnifold();
    TH1D* BinnedOmnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations);
    std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> UnbinnedOmnifold(TVectorD MC_entries,
                                                                        TVectorD sim_entries,
                                                                        TVectorD measured_entries,
                                                                        TVector pass_reco,
                                                                        TVector pass_truth,
                                                                        Int_t num_iterations);
    std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> UnbinnedOmnifold(TObjArray MC_entries,
                                                                          TObjArray sim_entries,
                                                                          TObjArray measured_entries,
                                                                          TVector pass_reco,
                                                                          TVector pass_truth,
                                                                          Int_t num_iterations);

    std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> 
                                                         UnbinnedOmnifold(std::vector<std::vector<Double_t>> MC_entries,
                                                                          std::vector<std::vector<Double_t>> sim_entries,
                                                                          std::vector<std::vector<Double_t>> measured_entries,
                                                                          std::vector<Bool_t> pass_reco,
                                                                          std::vector<Bool_t> pass_truth,
                                                                          Int_t num_iterations);

    std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> 
                                                         UnbinnedOmnifold(std::vector<Double_t> MC_entries,
                                                                          std::vector<Double_t> sim_entries,
                                                                          std::vector<Double_t> measured_entries,
                                                                          std::vector<Bool_t> pass_reco,
                                                                          std::vector<Bool_t> pass_truth,
                                                                          Int_t num_iterations);
    // void SetIterations(Int_t nIter) {_nIter = nIter;}
    // Int_t GetIterations() const {return _nIter;}
    // void SetMeasuredHist(TH1* measured_hist) {_measuredHist = measured_hist;}
    // TH1* GetMeasuredHist() const {return _measuredHist;}
    // void SetResponseMatrix(RooUnfoldResponse response) {_response = response;}
    // RooUnfoldResponse GetResponseMatrix() const {return _response;}
    void EfficiencyCorrections(TH1* hist, RooUnfoldResponse response);
};
#endif
