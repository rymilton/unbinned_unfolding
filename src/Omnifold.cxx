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

Omnifold::~Omnifold()
{
}

// Performs binned unfolding using current response and measured histograms
// Returns an unfolded TH1* with efficiency corrections applied
TH1D* Omnifold::BinnedOmnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations)
{
    TH2 *response_hist = response.HresponseNoOverflow();
    // Sending histograms to Python
    TPython::Bind( response_hist, "response_hist" ); 
    TPython::Bind( measured_hist, "measured_hist" );
    // Converting num_iterations to a TObject* to convert to Python
    TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
    TPython::Bind(num_iterations_object, "num_iterations" );
    // Performing binned Omnifold
    TPython::Exec("unfolded_hist = binned_omnifold(response_hist, measured_hist, num_iterations.GetVal())");
    // Bringing histogram back to ROOT and correcting for efficiency
    TH1D *unfolded_hist = TPython::Eval("unfolded_hist");
    this->EfficiencyCorrections(unfolded_hist, response);
    delete num_iterations_object;
    
    return unfolded_hist;
}
// Divides histogram bin content by efficiency vector from current response
void Omnifold::EfficiencyCorrections(TH1* hist, RooUnfoldResponse response)
{
    auto efficiency = response.Vefficiency();
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
    // std::vector<std::vector<Double_t>> to a TObjArray filled with TVectorDs
    auto nested_vector_to_TObjArray = [](std::vector<std::vector<Double_t>> vec)
    {
        TObjArray out_TObjArray(vec.size());
        for (size_t i = 0; i < vec.size(); i++)
        {
            TVectorD *subvec = new TVectorD(vec[i].size());
            for(int j = 0; j < vec[i].size(); j++)
            {
                (*subvec)[j] = vec[i][j];
            }
            out_TObjArray[i] = subvec;
        }
        return out_TObjArray;
    };

    auto TObjArray_to_nested_vector = [](TObjArray objarray)
    {
        std::vector<std::vector<Double_t>> nested_vector(objarray.GetEntries());
        for (size_t i = 0; i < objarray.GetEntries(); i++)
        {
            TVectorD subTVectorD = *(TVectorD*) objarray[i];
            std::vector<Double_t> subvec(subTVectorD.GetNoElements());
            for(int j = 0; j < subTVectorD.GetNoElements(); j++)
            {
                subvec[j] = subTVectorD[j];
            }
            nested_vector[i] = subvec;
        }
        return nested_vector;
    };

    auto MC_TObjArray = nested_vector_to_TObjArray(MC_entries);
    auto sim_TObjArray = nested_vector_to_TObjArray(sim_entries);
    auto measured_TObjArray = nested_vector_to_TObjArray(measured_entries);

    TVector pass_reco_TVector(pass_reco.size());
    for (size_t i = 0; i < pass_reco.size(); i++) {
        pass_reco_TVector[i] = pass_reco[i];
    }

    TVector pass_truth_TVector(pass_truth.size());
    for (size_t i = 0; i < pass_truth.size(); i++) {
        pass_truth_TVector[i] = pass_truth[i];
    }

    // Running unbinned unfolding
    std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> out_tuple_TObject = this->UnbinnedOmnifold(MC_TObjArray,
                                                                                                    sim_TObjArray,
                                                                                                    measured_TObjArray,
                                                                                                    pass_reco_TVector,
                                                                                                    pass_truth_TVector,
                                                                                                    num_iterations);

    TVectorD weights_MC_TVectorD = std::get<0>(out_tuple_TObject);
    TVectorD weights_sim_TVectorD = std::get<2>(out_tuple_TObject);

    // Converting TObjects back to std::vectors
    auto MC_vector = TObjArray_to_nested_vector(std::get<1>(out_tuple_TObject));
    auto sim_vector = TObjArray_to_nested_vector(std::get<3>(out_tuple_TObject));
    
    std::vector<Double_t> weights_MC_vector(weights_MC_TVectorD.GetNoElements());
    std::vector<Double_t> weights_sim_vector(weights_sim_TVectorD.GetNoElements());
    for(int i = 0; i < weights_MC_TVectorD.GetNoElements(); i++)
    {
        weights_MC_vector[i] = weights_MC_TVectorD[i];
    }

    for(int i = 0; i < weights_sim_TVectorD.GetNoElements(); i++)
    {
        weights_sim_vector[i] = weights_sim_TVectorD[i];
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

std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> 
                                                         Omnifold::UnbinnedOmnifold(std::vector<Double_t> MC_entries,
                                                                          std::vector<Double_t> sim_entries,
                                                                          std::vector<Double_t> measured_entries,
                                                                          std::vector<Bool_t> pass_reco,
                                                                          std::vector<Bool_t> pass_truth,
                                                                          Int_t num_iterations)
{

    auto nest_vector = [](std::vector<Double_t> vec)
    {
        std::vector<std::vector<Double_t>> nested_vector;
        for(const Double_t& entry: vec)
        {
            std::vector<Double_t> dummy_vector{entry};
            nested_vector.push_back(dummy_vector);
        }
        return nested_vector;
    };

    auto unnest_vector = [](std::vector<std::vector<Double_t>> vec)
    {
        std::vector<Double_t> unnested_vector;
        for(const std::vector<Double_t>& entry: vec)
        {
            unnested_vector.push_back(entry[0]);
        }
        return unnested_vector;
    };
    
    auto nested_MC_entries = nest_vector(MC_entries);
    auto nested_sim_entries = nest_vector(sim_entries);
    auto nested_measured_entries = nest_vector(measured_entries);

    std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> nested_out_tuple = this->UnbinnedOmnifold(nested_MC_entries,
                                                                                                                                                                               nested_sim_entries,
                                                                                                                                                                               nested_measured_entries,
                                                                                                                                                                               pass_reco,
                                                                                                                                                                               pass_truth,
                                                                                                                                                                               num_iterations);

    auto unnested_out_MC_entries = unnest_vector(std::get<1>(nested_out_tuple));
    auto unnested_out_sim_entries = unnest_vector(std::get<3>(nested_out_tuple));
    std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> out_tuple = std::make_tuple(std::get<0>(nested_out_tuple),
                                                                                                                                       unnested_out_MC_entries,
                                                                                                                                       std::get<2>(nested_out_tuple),
                                                                                                                                       unnested_out_sim_entries);

    return out_tuple;                                                                                                                                       
}
                                                                       