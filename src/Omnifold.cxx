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
    TPython::LoadMacro( "../python/omnifold.py" );
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