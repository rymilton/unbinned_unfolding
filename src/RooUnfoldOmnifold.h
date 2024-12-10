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
#include <ROOT/RDataFrame.hxx>

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

    virtual RooUnfolding::Algorithm GetAlgorithm() const override;

    std::tuple<TVectorD, TVectorD> UnbinnedOmnifold();
    std::tuple<TVectorD, TVectorD> UnbinnedOmnifold(ROOT::RDataFrame MC_dataframe,
                          ROOT::RDataFrame sim_dataframe,
                          ROOT::RDataFrame measured_dataframe);
    
    void SetNumIterations(int num_iterations){this->_niter = num_iterations;}
    void SetMCDataFrame(ROOT::RDataFrame& MC){this->_MCDataFrame = MC;}
    void SetSimDataFrame(ROOT::RDataFrame& Sim){this->_SimDataFrame = Sim;};
    void SetMeasuredDataFrame(ROOT::RDataFrame& Measured){this->_MeasuredDataFrame = Measured;};
    void SetMCPassReco(TVector& MC_pass_reco){this->_MCPassReco.ResizeTo(MC_pass_reco);this->_MCPassReco = MC_pass_reco;}
    void SetMCPassTruth(TVector& MC_pass_truth){this->_MCPassTruth.ResizeTo(MC_pass_truth); this->_MCPassTruth = MC_pass_truth;};
    void SetMeasuredPassReco(TVector& measured_pass_reco){this->_MeasuredPassReco.ResizeTo(measured_pass_reco); this->_MeasuredPassReco = measured_pass_reco;};

    int GetNumIterations() {return this->_niter;}
    ROOT::RDataFrame GetMCDataFrame() { return this->_MCDataFrame;};
    ROOT::RDataFrame GetSimDataFrame()  { return this->_SimDataFrame;};
    ROOT::RDataFrame GetMeasuredDataFrame() { return this->_MeasuredDataFrame;};
    TVector GetMCPassReco() { return this->_MCPassReco;};
    TVector GetMCPassTruth()  { return this->_MCPassTruth;};
    TVector GetMeasuredPassReco() { return this->_MeasuredPassReco;};
    TVectorD GetUnbinnedStep1Weights() {return this->_unbinned_step1_weights;};
    TVectorD GetUnbinnedStep2Weights() {return this->_unbinned_step2_weights;};

protected:
    virtual void Unfold() const override;
    // void EfficiencyCorrections(TH1* hist, RooUnfoldResponseT<Hist, Hist2D>* response);
    void BinnedOmnifold() const;

private:
    void Init();
    ROOT::RDataFrame _MCDataFrame;
    ROOT::RDataFrame _SimDataFrame;
    ROOT::RDataFrame _MeasuredDataFrame;
    TVector _MCPassReco;
    TVector _MCPassTruth;
    TVector _MeasuredPassReco;
    TVectorD _unbinned_step1_weights;
    TVectorD _unbinned_step2_weights;

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
