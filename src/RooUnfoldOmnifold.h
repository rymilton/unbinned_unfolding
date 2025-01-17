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
#include <TMap.h>

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

    RooUnfoldOmnifoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t niter= 4, bool useDensity=false,
                        const char* name= 0, const char* title= 0 );

    virtual RooUnfolding::Algorithm GetAlgorithm() const override;

    std::tuple<TVectorD, TVectorD> UnbinnedOmnifold();
    std::tuple<TVectorD, TVectorD> UnbinnedOmnifold(ROOT::RDataFrame MCgen_dataframe,
                          ROOT::RDataFrame MCreco_dataframe,
                          ROOT::RDataFrame measured_dataframe);
    std::tuple<TVectorD, TVectorD> TestUnbinnedOmnifold();

    
    void SetNumIterations(int num_iterations){this->_niter = num_iterations;}
    void SetMCgenDataFrame(ROOT::RDataFrame& MCgen){this->_MCgenDataFrame = MCgen;}
    void SetMCrecoDataFrame(ROOT::RDataFrame& MCreco){this->_MCrecoDataFrame = MCreco;};
    void SetMeasuredDataFrame(ROOT::RDataFrame& Measured){this->_MeasuredDataFrame = Measured;};
    void SetMCPassReco(TVector& MC_pass_reco){this->_MCPassReco.ResizeTo(MC_pass_reco);this->_MCPassReco = MC_pass_reco;}
    void SetMCPassTruth(TVector& MC_pass_truth){this->_MCPassTruth.ResizeTo(MC_pass_truth); this->_MCPassTruth = MC_pass_truth;};
    void SetMeasuredPassReco(TVector& measured_pass_reco){this->_MeasuredPassReco.ResizeTo(measured_pass_reco); this->_MeasuredPassReco = measured_pass_reco;};
    void SetMCgenWeights(TVectorD& MCgen_weights){this->_MCgenWeights.ResizeTo(MCgen_weights);this->_MCgenWeights = MCgen_weights;}
    void SetMCrecoWeights(TVectorD& MCreco_weights){this->_MCrecoWeights.ResizeTo(MCreco_weights);this->_MCrecoWeights = MCreco_weights;}
    void SetMeasuredWeights(TVectorD& measured_weights){this->_MeasuredWeights.ResizeTo(measured_weights);this->_MeasuredWeights = measured_weights;}
    void SetModelSaving(bool option){this->_SaveUnbinnedModels = option;};
    void SetSaveDirectory(TString save_dir){this->_UnbinnedModelSaveDir = save_dir;};
    void SetModelName(TString model_name){this->_UnbinnedModelName = model_name;};
    void SetTestMCgenDataFrame(ROOT::RDataFrame& MCgen){this->_TestMCgenDataFrame = MCgen;}
    void SetTestMCrecoDataFrame(ROOT::RDataFrame& MCreco){this->_TestMCrecoDataFrame = MCreco;};
    void SetTestMCPassReco(TVector& MC_pass_reco){this->_TestMCPassReco.ResizeTo(MC_pass_reco);this->_TestMCPassReco = MC_pass_reco;}
    void SetStep1ClassifierParameters(TMap* parameters){this->_Step1ClassifierParameters=parameters;};
    void SetStep2ClassifierParameters(TMap* parameters){this->_Step2ClassifierParameters=parameters;};
    void SetStep1RegressorParameters(TMap* parameters){this->_Step1RegressorParameters=parameters;};

    int GetNumIterations() {return this->_niter;}
    ROOT::RDataFrame GetMCgenDataFrame() { return this->_MCgenDataFrame;};
    ROOT::RDataFrame GetMCrecoDataFrame()  { return this->_MCrecoDataFrame;};
    ROOT::RDataFrame GetMeasuredDataFrame() { return this->_MeasuredDataFrame;};
    TVector GetMCPassReco() { return this->_MCPassReco;};
    TVector GetMCPassTruth()  { return this->_MCPassTruth;};
    TVector GetMeasuredPassReco() { return this->_MeasuredPassReco;};
    TVectorD GetMCgenWeights() { return this->_MCgenWeights;};
    TVectorD GetMCrecoWeights()  { return this->_MCrecoWeights;};
    TVectorD GetMeasuredWeights() { return this->_MeasuredWeights;};

    TVectorD GetUnbinnedStep1Weights() {return this->_unbinned_step1_weights;};
    TVectorD GetUnbinnedStep2Weights() {return this->_unbinned_step2_weights;};

protected:
    virtual void Unfold() const override;
    void BinnedOmnifold() const;

private:
    void Init();
    bool _useDensity;
    ROOT::RDataFrame _MCgenDataFrame;
    ROOT::RDataFrame _MCrecoDataFrame;
    ROOT::RDataFrame _MeasuredDataFrame;
    TVector _MCPassReco;
    TVector _MCPassTruth;
    TVector _MeasuredPassReco;
    TVectorD _unbinned_step1_weights;
    TVectorD _unbinned_step2_weights;
    bool _SaveUnbinnedModels;
    TString _UnbinnedModelSaveDir;
    TString _UnbinnedModelName;
    ROOT::RDataFrame _TestMCgenDataFrame;
    ROOT::RDataFrame _TestMCrecoDataFrame;
    TVector _TestMCPassReco;
    TVectorD _unbinned_step1_test_weights;
    TVectorD _unbinned_step2_test_weights;
    TMap* _Step1ClassifierParameters;
    TMap* _Step2ClassifierParameters;
    TMap* _Step1RegressorParameters;
    TVectorD _MCgenWeights;
    TVectorD _MCrecoWeights;
    TVectorD _MeasuredWeights;


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
