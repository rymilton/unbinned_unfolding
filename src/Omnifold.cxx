#include <iostream>
#include "Omnifold.h"
#include "TParameter.h"
#include "TH1.h"
#include "TPython.h"
#include "TObjArray.h"

Omnifold::Omnifold()
{
    TPython::LoadMacro( "../python/omnifold.py" );
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
    delete num_iterations_object;
    
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
std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> Omnifold::UnbinnedOmnifold(TVectorD MC_entries,
                                                                              TVectorD sim_entries,
                                                                              TVectorD measured_entries,
                                                                              TVector pass_reco,
                                                                              TVector pass_truth,
                                                                              Int_t num_iterations)
{
    TPython::Bind( &MC_entries, "unbinned_MC_entries" );
    TPython::Exec("unbinned_MC_entries = np.array(unbinned_MC_entries)");
    TPython::Bind( &sim_entries, "unbinned_sim_entries" );
    TPython::Exec("unbinned_sim_entries = np.array(unbinned_sim_entries)");
    TPython::Bind( &measured_entries, "unbinned_measured_entries" );
    TPython::Exec("unbinned_measured_entries = np.array(unbinned_measured_entries)");
    TPython::Bind( &pass_reco, "pass_reco_mask" );
    TPython::Exec( "pass_reco_mask = np.array(pass_reco_mask, dtype=bool)" );
    TPython::Bind( &pass_truth, "pass_truth_mask" );
    TPython::Exec( "pass_truth_mask = np.array(pass_truth_mask, dtype=bool)" );
    TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
    TPython::Bind(num_iterations_object, "num_iterations" );
    TPython::Exec("weights_test, MC_test, sim_test, pass_reco_test = unbinned_omnifold(unbinned_MC_entries,\
                                                                                       unbinned_sim_entries,\
                                                                                       unbinned_measured_entries,\
                                                                                       pass_reco_mask,\
                                                                                       pass_truth_mask,\
                                                                                       num_iterations.GetVal())");
    TPython::Exec("weights_MC_TVectorD = convert_to_TVectorD(weights_test[-1, 1])");
    TPython::Exec("MC_TVectorD = convert_to_TVectorD(MC_test.flatten())");
    TPython::Exec("weights_sim_TVectorD = convert_to_TVectorD(weights_test[-1, 0][pass_reco_test])");
    TPython::Exec("sim_TVectorD = convert_to_TVectorD(sim_test.flatten())");
    
    TVectorD* weights_MC_test = TPython::Eval("weights_MC_TVectorD");
    TVectorD* MC_test = TPython::Eval("MC_TVectorD");
    TVectorD* weights_sim_test = TPython::Eval("weights_sim_TVectorD");
    TVectorD* sim_test = TPython::Eval("sim_TVectorD");

    TVectorD out_weights_MC_test = *weights_MC_test;
    TVectorD out_MC_test = *MC_test;
    TVectorD out_weights_sim_test = *weights_sim_test;
    TVectorD out_sim_test = *sim_test;

    delete num_iterations_object;

    std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> out_tuple = std::make_tuple(out_weights_MC_test,
                                                                                   out_MC_test,
                                                                                   out_weights_sim_test,
                                                                                   out_sim_test);
    return out_tuple;
}
std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> Omnifold::UnbinnedOmnifold(TObjArray MC_entries,
                                                                                     TObjArray sim_entries,
                                                                                     TObjArray measured_entries,
                                                                                     TVector pass_reco,
                                                                                     TVector pass_truth,
                                                                                     Int_t num_iterations)
{
    TPython::Bind( &MC_entries, "unbinned_MC_entries" );
    TPython::Exec("unbinned_MC_entries = np.array(unbinned_MC_entries)");
    TPython::Bind( &sim_entries, "unbinned_sim_entries" );
    TPython::Exec("unbinned_sim_entries = np.array(unbinned_sim_entries)");
    TPython::Bind( &measured_entries, "unbinned_measured_entries" );
    TPython::Exec("unbinned_measured_entries = np.array(unbinned_measured_entries)");
    TPython::Bind( &pass_reco, "pass_reco_mask" );
    TPython::Exec( "pass_reco_mask = np.array(pass_reco_mask, dtype=bool)" );
    TPython::Bind( &pass_truth, "pass_truth_mask" );
    TPython::Exec( "pass_truth_mask = np.array(pass_truth_mask, dtype=bool)" );
    TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
    TPython::Bind(num_iterations_object, "num_iterations" );

    TPython::Exec("weights_test, MC_test, sim_test, pass_reco_test = unbinned_omnifold(unbinned_MC_entries,\
                                                                                       unbinned_sim_entries,\
                                                                                       unbinned_measured_entries,\
                                                                                       pass_reco_mask,\
                                                                                       pass_truth_mask,\
                                                                                       num_iterations.GetVal())");

    TPython::Exec("weights_MC_TVectorD = convert_to_TVectorD(weights_test[-1, 1])");
    TPython::Exec("weights_sim_TVectorD = convert_to_TVectorD(weights_test[-1, 0][pass_reco_test])");

    TPython::Exec("MC_vectors = get_vectors(MC_test)");
    TPython::Exec("MC_TObjArray = ROOT.TObjArray(len(MC_vectors))\n"
                  "for i in range(len(MC_vectors)):\n"
                  "    MC_TObjArray[i] = MC_vectors[i]");
    TObjArray out_MC_test = *(TObjArray*) TPython::Eval("MC_TObjArray");

    TPython::Exec("sim_vectors = get_vectors(sim_test)");
    TPython::Exec("sim_TObjArray = ROOT.TObjArray(len(sim_vectors))\n"
                  "for i in range(len(sim_vectors)):\n"
                  "    sim_TObjArray[i] = sim_vectors[i]");
    TObjArray out_sim_test = *(TObjArray*) TPython::Eval("sim_TObjArray");
    
    TVectorD out_weights_MC_test = *(TVectorD*) TPython::Eval("weights_MC_TVectorD");
    TVectorD out_weights_sim_test = *(TVectorD*) TPython::Eval("weights_sim_TVectorD");

    std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> out_tuple = std::make_tuple(out_weights_MC_test,
                                                                                   out_MC_test,
                                                                                   out_weights_sim_test,
                                                                                   out_sim_test);
    return out_tuple;
}

// Does unbinned unfolding with std::vectors as input and output
// Converts the std::vectors to TObjects and then does the unfolding with TPython
std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>>
                                               Omnifold::UnbinnedOmnifold(std::vector<std::vector<Double_t>> MC_entries,
                                                                          std::vector<std::vector<Double_t>> sim_entries,
                                                                          std::vector<std::vector<Double_t>> measured_entries,
                                                                          std::vector<Bool_t> pass_reco,
                                                                          std::vector<Bool_t> pass_truth,
                                                                          Int_t num_iterations)
{
    // Converting the std::vectors to TObjArrays and TVectors
    TObjArray MC_TObjArray(MC_entries.size());
    TObjArray sim_TObjArray(sim_entries.size());
    TObjArray measured_TObjArray(measured_entries.size());

    for (int i = 0; i < MC_entries.size(); i++) {
        TVectorD *MC_event = new TVectorD(MC_entries[i].size());
        for(int j = 0; j < MC_entries[i].size(); j++)
        {
            (*MC_event)[j] = MC_entries[i][j];
        }
        MC_TObjArray[i] = MC_event;
    }

    for (int i = 0; i < sim_entries.size(); i++) {
        TVectorD *sim_event = new TVectorD(sim_entries[i].size());
        for(int j = 0; j < sim_entries[i].size(); j++)
        {
            (*sim_event)[j] = sim_entries[i][j];
        }
        sim_TObjArray[i] = sim_event;
    }

    for (int i = 0; i < measured_entries.size(); i++) {
        TVectorD *measured_event = new TVectorD(measured_entries[i].size());
        for(int j = 0; j < measured_entries[i].size(); j++)
        {
            (*measured_event)[j] = measured_entries[i][j];
        }
        measured_TObjArray[i] = measured_event;
    }

    TVector pass_reco_TVector(pass_reco.size());
    for (int i = 0; i < pass_reco.size(); i++) {
        pass_reco_TVector[i] = pass_reco[i];
    }

    TVector pass_truth_TVector(pass_truth.size());
    for (int i = 0; i < pass_truth.size(); i++) {
        pass_truth_TVector[i] = pass_truth[i];
    }

    // Running unbinned unfolding
    std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> out_tuple_TObject = this->UnbinnedOmnifold(MC_TObjArray,
                                                                                                    sim_TObjArray,
                                                                                                    measured_TObjArray,
                                                                                                    pass_reco_TVector,
                                                                                                    pass_truth_TVector,
                                                                                                    num_iterations);

    TVectorD weights_MC_TObject = std::get<0>(out_tuple_TObject);
    TObjArray MC_TObject = std::get<1>(out_tuple_TObject);
    TVectorD weights_sim_TObject = std::get<2>(out_tuple_TObject);
    TObjArray sim_TObject = std::get<3>(out_tuple_TObject); 

    // Converting TObjects back to std::vectors
    std::vector<std::vector<Double_t>> MC_vector(MC_TObject.GetEntries());
    std::vector<std::vector<Double_t>> sim_vector(sim_TObject.GetEntries());
    for(int i = 0; i < MC_TObject.GetEntries(); i++)
    {
        TVectorD event_TVectorD = *(TVectorD*) MC_TObject[i];
        std::vector<Double_t> event_vector(event_TVectorD.GetNoElements());
        for(int j = 0; j < event_TVectorD.GetNoElements(); j++)
        {
            event_vector[j] = event_TVectorD[j];
        }
        MC_vector[i] = event_vector;
    }

    for(int i = 0; i < sim_TObject.GetEntries(); i++)
    {
        TVectorD event_TVectorD = *(TVectorD*) sim_TObject[i];
        std::vector<Double_t> event_vector(event_TVectorD.GetNoElements());
        for(int j = 0; j < event_TVectorD.GetNoElements(); j++)
        {
            event_vector[j] = event_TVectorD[j];
        }
        sim_vector[i] = event_vector;
    }
    std::vector<Double_t> weights_MC_vector(weights_MC_TObject.GetNoElements());
    std::vector<Double_t> weights_sim_vector(weights_sim_TObject.GetNoElements());
    for(int i = 0; i < weights_MC_TObject.GetNoElements(); i++)
    {
        weights_MC_vector[i] = weights_MC_TObject[i];
    }

    for(int i = 0; i < weights_sim_TObject.GetNoElements(); i++)
    {
        weights_sim_vector[i] = weights_sim_TObject[i];
    }

    for (int i = 0; i < MC_TObjArray.GetEntries(); i++) {
        delete MC_TObjArray[i];
    }
    for (int i = 0; i < sim_TObjArray.GetEntries(); i++) {
        delete sim_TObjArray[i];
    }
    for (int i = 0; i < measured_TObjArray.GetEntries(); i++) {
        delete measured_TObjArray[i];
    }

    std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> out_tuple_vector = std::make_tuple(weights_MC_vector,
                                                                                                                                                                 MC_vector,
                                                                                                                                                                 weights_sim_vector,
                                                                                                                                                                 sim_vector);

    return out_tuple_vector;                                                                                                                                                                   
}
                                                                       