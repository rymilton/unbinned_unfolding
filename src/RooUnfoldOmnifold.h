//=====================================================================-*-C++-*-
//! \class Omnifold
//! \brief Unbinned & binned versions of Omnifold (ML unfolding) using decision trees
//! \author Ryan Milton <rmilt003@ucr.edu>
//==============================================================================

#ifndef ROOUNFOLDOMNIFOLD_H_
#define ROOUNFOLDOMNIFOLD_H_

#include <tuple>
#include "TH1.h"
#include "TVector.h"
#include "TVectorD.h"
#include "RooUnfold.h"

#include "RooUnfoldResponse.h"
#include "TObjArray.h"

class TH1;
class TH2;
template<class Hist, class Hist2D>
class RooUnfoldOmnifoldT : public RooUnfoldT<Hist,Hist2D> {
public:
    RooUnfoldOmnifoldT(); // default constructor
    RooUnfoldOmnifoldT (const char*    name, const char*    title); // named constructor
    RooUnfoldOmnifoldT (const TString& name, const TString& title); // named constructor
    RooUnfoldOmnifoldT (const RooUnfoldOmnifoldT<Hist,Hist2D>& rhs); // copy constructor
    RooUnfoldOmnifoldT& operator= (const RooUnfoldOmnifoldT<Hist,Hist2D>& rhs); // assignment operator

    RooUnfoldOmnifoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t niter= 4,
                        const char* name= 0, const char* title= 0);


    // std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> UnbinnedOmnifold(TVectorD MC_entries,
    //                                                                     TVectorD sim_entries,
    //                                                                     TVectorD measured_entries,
    //                                                                     TVector pass_reco,
    //                                                                     TVector pass_truth,
    //                                                                     Int_t num_iterations);
    // std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> UnbinnedOmnifold(TObjArray MC_entries,
    //                                                                       TObjArray sim_entries,
    //                                                                       TObjArray measured_entries,
    //                                                                       TVector pass_reco,
    //                                                                       TVector pass_truth,
    //                                                                       Int_t num_iterations);

    // std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> 
    //                                                      UnbinnedOmnifold(std::vector<std::vector<Double_t>> MC_entries,
    //                                                                       std::vector<std::vector<Double_t>> sim_entries,
    //                                                                       std::vector<std::vector<Double_t>> measured_entries,
    //                                                                       std::vector<Bool_t> pass_reco,
    //                                                                       std::vector<Bool_t> pass_truth,
    //                                                                       Int_t num_iterations);

    // std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> 
    //                                                      UnbinnedOmnifold(std::vector<Double_t> MC_entries,
    //                                                                       std::vector<Double_t> sim_entries,
    //                                                                       std::vector<Double_t> measured_entries,
    //                                                                       std::vector<Bool_t> pass_reco,
    //                                                                       std::vector<Bool_t> pass_truth,
    //                                                                       Int_t num_iterations);
    // void EfficiencyCorrections(TH1* hist, RooUnfoldResponse response);

protected:
    virtual void Unfold() const override;
    // void EfficiencyCorrections(TH1* hist, RooUnfoldResponseT<Hist, Hist2D>* response);
    void BinnedOmnifold() const;

private:
    void Init();
    
protected:
    mutable int _niter;
public:
  ClassDefOverride (RooUnfoldOmnifoldT, 1)
};

typedef RooUnfoldOmnifoldT<TH1,TH2> RooUnfoldOmnifold;
#ifndef NOROOFIT
//! \class RooFitUnfoldBayes
//! \brief specialization of RooUnfoldBayesT for RooAbsReal objects
typedef RooUnfoldOmnifoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist> RooFitUnfoldOmnifold;
#endif
#endif
