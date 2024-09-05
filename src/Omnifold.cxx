#include <iostream>
#include "Omnifold.h"
#include "TParameter.h"
#include "TH1.h"
#include "TPython.h"

Omnifold::Omnifold()
{
}
Omnifold::Omnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations)
{
    this->SetMeasuredHist(measured_hist);
    this->SetResponseMatrix(response);
    this->SetIterations(num_iterations);
    TPython::LoadMacro( "../python/omnifold.py" );
}
Omnifold::~Omnifold()
{
}

// Performs binned unfolding using current response and measured histograms
// Returns an unfolded TH1* with efficiency corrections applied
TH1D* Omnifold::BinnedOmnifold()
{
    TH2 *response_hist = this->GetResponseMatrix().HresponseNoOverflow();
    // Sending histograms to Python
    TPython::Bind( response_hist, "response_hist" ); 
    TPython::Bind( this->GetMeasuredHist(), "measured_hist" );
    // Converting num_iterations to a TObject* to convert to Python
    TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", this->GetIterations());
    TPython::Bind(num_iterations_object, "num_iterations" );
    // Performing binned Omnifold
    TPython::Exec("unfolded_hist = binned_omnifold(response_hist, measured_hist, num_iterations.GetVal())");
    // Bringing histogram back to ROOT and correcting for efficiency
    TH1D *unfolded_hist = TPython::Eval("unfolded_hist");
    this->EfficiencyCorrections(unfolded_hist);
    return unfolded_hist;
}
// Divides histogram bin content by efficiency vector from current response
void Omnifold::EfficiencyCorrections(TH1* hist)
{
    auto efficiency = this->GetResponseMatrix().Vefficiency();
    for (Int_t i = 0; i < hist->GetNbinsX(); i++)
    {
        Double_t corrected_content = efficiency[i] > 0 ?
            hist->GetBinContent(i+1)/efficiency[i] : hist->GetBinContent(i+1);
        hist->SetBinContent(i+1, corrected_content);
    }
}
std::tuple<TVectorD, TVectorD> Omnifold::UnbinnedOmnifold(TVectorD MC_entries, TVectorD sim_entries, TVectorD measured_entries, TVector pass_reco, TVector pass_truth, Int_t num_iterations)
{
    TPython::Bind( &MC_entries, "MC_entries" );
    TPython::Bind( &MC_entries, "unbinned_MC_entries" );
    TPython::Bind( &sim_entries, "unbinned_sim_entries" );
    TPython::Bind( &measured_entries, "unbinned_measured_entries" );
    TPython::Bind( &pass_reco, "pass_reco_mask" );
    TPython::Exec( "pass_reco_mask = np.array(pass_reco_mask, dtype=bool)" );
    TPython::Bind( &pass_truth, "pass_truth_mask" );
    TPython::Exec( "pass_truth_mask = np.array(pass_truth_mask, dtype=bool)" );
    TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
    TPython::Bind(num_iterations_object, "num_iterations" );
    TPython::Exec("weights, MC_data, _, _ = unbinned_omnifold(unbinned_MC_entries, unbinned_sim_entries, unbinned_measured_entries, pass_reco_mask, pass_truth_mask, num_iterations.GetVal())");
    TPython::Exec("weights_MC_ROOT = convert_to_TVectorD(weights[-1][1])");
    TPython::Exec("unbinned_MC_ROOT = convert_to_TVectorD(MC_data.flatten()) ");
    TVectorD* weights_MC = TPython::Eval("weights_MC_ROOT");
    TVectorD* unbinned_MC = TPython::Eval("unbinned_MC_ROOT");
    TVectorD out_weights_MC = *weights_MC;
    TVectorD out_unbinned_MC = *unbinned_MC;
    std::tuple<TVectorD, TVectorD> out_tuple = std::make_tuple(out_weights_MC, out_unbinned_MC);
    return out_tuple;
    // std::cout<<"hi"<<std::endl;
}
                                                                       