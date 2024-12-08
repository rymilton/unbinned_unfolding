/*
    Implementation of Omnifold using boosted decision trees (BDT)
    This follows the original Omnifold work (Andreassen et al. PRL 124, 182001 (2020))

    This class has BinnedOmnifold (binned unfolding) and UnbinnedOmnifold (unbinned unfolding).

    BinnedOmnifold is similar to RooUnfoldBayes, where it uses the response matrix and measured histogram to do unfolding 
    and it returns the unfolded TH1*.
    UnbinnedOmnifold instead takes a list of Monte Carlo data, reconstructed Monte Carlo data (called sim data), and
    measured data. It also takes masks that indicate an event passes generation level cuts and reconstruction.
    UnbinnedOmnifold is overloaded with 4 functions that adjust the data format. The four formats are:
        - 1D: TVectorD
        - 1D: std::vector<double>
        - Any dimension: TObjArray filled with TVectorD
        - Any dimension: std::vector<std::vector<double>>

    To do the unfolding, the input data is first moved to Python using TPython. This conversion needs TObjects, which
    is why TVectorD and TObjArray are used in UnbinnedOmnifold. The unfolding procedure is done in Python and the results
    are brought back to C++.

    BDTs are trained to output sets of weights that are applied to the data. 
    In BinnedOmnifold, these weights are applied within the function and an unfolded histogram is returned.
    In UnbinnedOmnifold, these weights are returned within an std::tuple of the form
    <Monte Carlo weights, Monte Carlo data, simulation/reconstructed weights, simulation/reconstructed data>.
    Within the std::tuple, there are either TObjects (TVectorD + TObjArray) or std::vectors, matching the input form.

    In both BinnedOmnifold and UnbinnedOmnifold, 50% of the Monte Carlo + simulation/reconstructed data is used for training
    the BDTs while the other 50% is used as test data to avoid overfitting.
    For UnbinnedOmnifold, the test data and its associated weights are returned in the std::tuple.
    
    To plot the UnbinnedOmnifold results in a histogram, the Monte Carlo data (second tuple entry) should be filled with
    the Monte Carlo weights (first tuple entry). The resulting histogram is the unfolded distribution. The measured/simulation
    data (fourth entry) should be filled with the simulation weights (third entry). This will compare to the measured data.
    The weights are event-by-event weights, so even if you have multi-dimensional data, there is only one weight for a given event.
    e.g. MC data = [[1, 2], [1.2, 2.5]], MC weights = [.5, .2]. The .5 (.2) is applied to both 1 and 2 (1.2 and 2.5) from [1, 2] ([1.2, 2.5]).
*/
#include <iostream>
#include "RooUnfoldOmnifold.h"
#include "TParameter.h"
#include "TH1.h"
#include "TPython.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <ROOT/RDataFrame.hxx>
using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT()
  : RooUnfoldT<Hist,Hist2D>(), _MCDataFrame(0), _SimDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0), _MCPassTruth(0), _MeasuredPassReco(0)
{

  //! Default constructor. Use Setup() to prepare for unfolding.]
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title), _MCDataFrame(0), _SimDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0), _MCPassTruth(0), _MeasuredPassReco(0)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title), _MCDataFrame(0), _SimDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0), _MCPassTruth(0), _MeasuredPassReco(0)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t niter,
                        const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _niter(niter), _MCDataFrame(0), _SimDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0), _MCPassTruth(0), _MeasuredPassReco(0)
{

  //! Constructor with response matrix object and measured unfolding input histogram.
  //! The regularisation parameter is niter (number of iterations).
  Init();
}

// Loading the omnifold.py script
template<class Hist,class Hist2D> void
RooUnfoldOmnifoldT<Hist,Hist2D>::Init()
{
  // A copy of ROOT/StringUtils.h since it's only available in ROOT 6.26 and later
  auto string_split = [](std::string str, std::string delims, bool skipEmpty /* =false */)
  {
    std::vector<std::string> out;

    std::size_t beg = 0;
    std::size_t end = 0;
    while ((end = str.find_first_of(delims, beg)) != std::string::npos) {
      if (!skipEmpty || end > beg)
          out.emplace_back(str.substr(beg, end - beg));
      beg = end + 1;
    }
    if (!skipEmpty || str.size() > beg)
      out.emplace_back(str.substr(beg, str.size() - beg));

   return out;
  };

  auto dynamic_paths = string_split(gSystem->GetDynamicPath(), ":", false);
  TString RooUnfold_install_path;
  for(auto &path : dynamic_paths)
  {
    TString RooUnfold_library_name = "libRooUnfold.so";
    if(gSystem->FindFile(path.c_str(), RooUnfold_library_name))
        RooUnfold_install_path = path;
  }
  TString omnifold_path;
  if (!RooUnfold_install_path.IsNull())
  {
    TString omnifold_script = "RooUnfold/omnifold.py";
    omnifold_path = gSystem->FindFile(RooUnfold_install_path, omnifold_script);
    if (omnifold_path.IsNull())
    {
      throw std::runtime_error("Cannot find omnifold.py. Please check your RooUnfold PATH environment variable directory,"
                                " and check whether /RooUnfold/omnifold.py is in that directory.");
    }
    else
    {
      TPython::Exec(Form("omnifold_path = '%s'", omnifold_path.Data()));
      TPython::Exec("exec(open(omnifold_path).read())");
    }
  }
  else
    throw std::runtime_error("Cannot find RooUnfold install directory. Did you source setup.sh?");
}

template<class Hist,class Hist2D> void
RooUnfoldOmnifoldT<Hist,Hist2D>::Unfold() const
{
  BinnedOmnifold();
}

// // RooUnfoldOmnifoldT::~RooUnfoldOmnifoldT()
// // {
// // }

template<class Hist,class Hist2D> RooUnfolding::Algorithm
RooUnfoldOmnifoldT<Hist,Hist2D>::GetAlgorithm () const
{
  //! return the unfolding algorithm used
  return kOmnifold;
}

// // Performs binned unfolding using current response and measured histograms
// // Returns an unfolded TH1* with efficiency corrections applied
template<class Hist,class Hist2D> void
RooUnfoldOmnifoldT<Hist,Hist2D>::BinnedOmnifold() const
{
  auto response = this->response();
  TH2D *response_hist = (TH2D*) response->HresponseNoOverflow();
  TH1D *measured_hist = (TH1D*) this->Hmeasured();
  // Sending histograms to Python
  TPython::Bind( response_hist, "response_hist" ); 
  TPython::Bind( measured_hist, "measured_hist" );
  // Converting num_iterations to a TObject* to convert to Python
  TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", _niter);
  TPython::Bind(num_iterations_object, "num_iterations" );
  // Performing binned Omnifold
  TPython::Exec("unfolded_hist = binned_omnifold(response_hist, measured_hist, num_iterations.GetVal())");
  // Bringing histogram back to ROOT and correcting for efficiency
  TH1D *unfolded_hist = TPython::Eval("unfolded_hist");
  // this->EfficiencyCorrections(unfolded_hist, response);
  auto efficiency = response->Vefficiency();
  for (Int_t i = 0; i < unfolded_hist->GetNbinsX(); i++)
  {
      Double_t corrected_content = efficiency[i] > 0 ?
          unfolded_hist->GetBinContent(i+1)/efficiency[i] : unfolded_hist->GetBinContent(i+1);
      unfolded_hist->SetBinContent(i+1, corrected_content);
  }
  delete num_iterations_object;
  
  const int num_bins = unfolded_hist->GetNbinsX();
  TVectorD unfolded_content(num_bins);
  for(int i = 0; i < num_bins; i++)
  {
      unfolded_content[i] = unfolded_hist->GetBinContent(i+1);
  }
  this->_cache._rec.ResizeTo(unfolded_content);
  this->_cache._rec = unfolded_content;
  this->_cache._unfolded= true;
}

template<class Hist,class Hist2D> std::tuple<TVectorD, TVectorD>
RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold()
{
  auto DataFrame_to_python = [](ROOT::RDataFrame df, TString data_type)
  {
    TPython::Exec(Form("data_type = '%s'", data_type.Data()));
    TPython::Exec("data_dict[data_type] = np.array([])");
    for(const auto& column_name : df.GetColumnNames())
    {
      std::string column_type = df.GetColumnType(column_name);
      std::cout<<column_type<<std::endl;
      TVectorD vector(*(df.Count()));
      if (column_type == "double")
      {
        auto column_values = *(df.Take<double>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("MC_vector = np.asarray(column_vector)");
        TPython::Exec("data_dict[data_type] = np.vstack([data_dict[data_type], column_vector] if data_dict[data_type].size else np.array([column_vector]))");
      }
      else if (column_type == "float")
      {
        auto column_values = *(df.Take<float>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("data_dict[data_type] = np.vstack([data_dict[data_type], column_vector] if data_dict[data_type].size else np.array([column_vector]))");
      }
      else if (column_type == "Long64_t")
      {
        auto column_values = *(df.Take<Long64_t>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("data_dict[data_type] = np.vstack([data_dict[data_type], column_vector] if data_dict[data_type].size else np.array([column_vector]))");
      }
      else if (column_type == "long")
      {
        auto column_values = *(df.Take<long>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("data_dict[data_type] = np.vstack([data_dict[data_type], column_vector] if data_dict[data_type].size else np.array([column_vector]))");
      }
      else
        throw std::runtime_error("Please make sure RDataFrame columns are type double, float, or int");
    }
  };
  TPython::Exec("data_dict = {}");
  DataFrame_to_python(this->_MCDataFrame, "MC");
  DataFrame_to_python(this->_SimDataFrame, "sim");
  DataFrame_to_python(this->_MeasuredDataFrame, "measured");

  TPython::Exec(Form("num_iterations = '%d'", this->_niter));

  DataFrame_to_python(this->_MCPassReco, "MC_pass_reco");
  TPython::Exec("data_dict['MC_pass_reco'] = data_dict['MC_pass_reco'].flatten()");
  DataFrame_to_python(this->_MCPassTruth, "sim_pass_truth");
  TPython::Exec("data_dict['sim_pass_truth'] = data_dict['sim_pass_truth'].flatten()");
  DataFrame_to_python(this->_MeasuredPassReco, "measured_pass_reco");
  TPython::Exec("data_dict['measured_pass_reco'] = data_dict['measured_pass_reco'].flatten()");
  TPython::Exec("print(data_dict['MC_pass_reco'])");
  TPython::Exec("weights, test_MC, test_sim, pass_reco_test = unbinned_omnifold(data_dict['MC'],\
                                                                                data_dict['sim'],\
                                                                                data_dict['measured'],\
                                                                                num_iterations,\
                                                                                data_dict['MC_pass_reco'],\
                                                                                data_dict['sim_pass_truth'],\
                                                                                data_dict['measured_pass_reco'])");
  TPython::Exec("weights_MC_TVectorD = convert_to_TVectorD(weights[-1, 1])");
  TPython::Exec("MC_TVectorD = convert_to_TVectorD(test_MC.flatten())");
  TPython::Exec("weights_sim_TVectorD = convert_to_TVectorD(weights[-1, 0][pass_reco_test])");
  TPython::Exec("sim_TVectorD = convert_to_TVectorD(test_sim.flatten())");
    
  TVectorD out_weights_MC_test = *(TVectorD*) TPython::Eval("weights_MC_TVectorD");
  TVectorD out_MC_test = *(TVectorD*) TPython::Eval("MC_TVectorD");
  TVectorD out_weights_sim_test = *(TVectorD*) TPython::Eval("weights_sim_TVectorD");
  TVectorD out_sim_test = *(TVectorD*) TPython::Eval("sim_TVectorD");

  std::tuple<TVectorD, TVectorD> out_tuple = std::make_tuple(out_weights_MC_test, out_MC_test);
  return out_tuple;
}

template<class Hist,class Hist2D> std::tuple<TVectorD, TVectorD>
RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(ROOT::RDataFrame MC_dataframe,
                                                  ROOT::RDataFrame sim_dataframe,
                                                  ROOT::RDataFrame measured_dataframe)
{
  this->SetMCDataFrame(MC_dataframe);
  this->SetSimDataFrame(sim_dataframe);
  this->SetMeasuredDataFrame(measured_dataframe);
  return this->UnbinnedOmnifold();
}


// // Divides histogram bin content by efficiency vector from current response
// template<class Hist,class Hist2D> void
// RooUnfoldOmnifoldT<Hist,Hist2D>::EfficiencyCorrections(TH1* hist, RooUnfoldResponseT<Hist, Hist2D>* response) const
// {
//     auto efficiency = response->Vefficiency();
//     for (Int_t i = 0; i < hist->GetNbinsX(); i++)
//     {
//         Double_t corrected_content = efficiency[i] > 0 ?
//             hist->GetBinContent(i+1)/efficiency[i] : hist->GetBinContent(i+1);
//         hist->SetBinContent(i+1, corrected_content);
//     }
// }

// // 1D unbinned unfolding using TVectors
// // The function will squeeze the vector from (N, ) to (N, 1) and then unsqueeze it
// // This will return a tuple of MC weights, MC data, sim weights, sim data, all as TVectorD 
// std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(TVectorD MC_entries,
//                                                                               TVectorD sim_entries,
//                                                                               TVectorD measured_entries,
//                                                                               TVector pass_reco,
//                                                                               TVector pass_truth,
//                                                                               Int_t num_iterations)
// {
//     TPython::Bind( &MC_entries, "unbinned_MC_entries" );
//     TPython::Exec("unbinned_MC_entries = np.array(unbinned_MC_entries)");
//     TPython::Bind( &sim_entries, "unbinned_sim_entries" );
//     TPython::Exec("unbinned_sim_entries = np.array(unbinned_sim_entries)");
//     TPython::Bind( &measured_entries, "unbinned_measured_entries" );
//     TPython::Exec("unbinned_measured_entries = np.array(unbinned_measured_entries)");
//     TPython::Bind( &pass_reco, "pass_reco_mask" );
//     TPython::Exec( "pass_reco_mask = np.array(pass_reco_mask, dtype=bool)" );
//     TPython::Bind( &pass_truth, "pass_truth_mask" );
//     TPython::Exec( "pass_truth_mask = np.array(pass_truth_mask, dtype=bool)" );
//     TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
//     TPython::Bind(num_iterations_object, "num_iterations" );
//     TPython::Exec("weights_test, MC_test, sim_test, pass_reco_test = unbinned_omnifold(unbinned_MC_entries,\
//                                                                                        unbinned_sim_entries,\
//                                                                                        unbinned_measured_entries,\
//                                                                                        pass_reco_mask,\
//                                                                                        pass_truth_mask,\
//                                                                                        num_iterations.GetVal())");
//     TPython::Exec("weights_MC_TVectorD = convert_to_TVectorD(weights_test[-1, 1])");
//     TPython::Exec("MC_TVectorD = convert_to_TVectorD(MC_test.flatten())");
//     TPython::Exec("weights_sim_TVectorD = convert_to_TVectorD(weights_test[-1, 0][pass_reco_test])");
//     TPython::Exec("sim_TVectorD = convert_to_TVectorD(sim_test.flatten())");
    
//     TVectorD out_weights_MC_test = *(TVectorD*) TPython::Eval("weights_MC_TVectorD");
//     TVectorD out_MC_test = *(TVectorD*) TPython::Eval("MC_TVectorD");
//     TVectorD out_weights_sim_test = *(TVectorD*) TPython::Eval("weights_sim_TVectorD");
//     TVectorD out_sim_test = *(TVectorD*) TPython::Eval("sim_TVectorD");

//     delete num_iterations_object;

//     std::tuple<TVectorD, TVectorD, TVectorD, TVectorD> out_tuple = std::make_tuple(out_weights_MC_test,
//                                                                                    out_MC_test,
//                                                                                    out_weights_sim_test,
//                                                                                    out_sim_test);
//     return out_tuple;
// }

// // Unbinned unfolding using TVectors for any dimension
// // To pass objects to Python via TPython, we must use TObjects, hence the use of TObjArray
// // Each TObjArray should be filled with a TVectorD of data
// // This will return a tuple of MC weights (TVectorD), MC data (TObjArray), sim weights (TVectorD), sim data (TObjArray)
// std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(TObjArray MC_entries,
//                                                                                 TObjArray sim_entries,
//                                                                                 TObjArray measured_entries,
//                                                                                 TVector pass_reco,
//                                                                                 TVector pass_truth,
//                                                                                 Int_t num_iterations)
// {
//     TPython::Bind( &MC_entries, "unbinned_MC_entries" );
//     TPython::Exec("unbinned_MC_entries = np.array(unbinned_MC_entries)");
//     TPython::Bind( &sim_entries, "unbinned_sim_entries" );
//     TPython::Exec("unbinned_sim_entries = np.array(unbinned_sim_entries)");
//     TPython::Bind( &measured_entries, "unbinned_measured_entries" );
//     TPython::Exec("unbinned_measured_entries = np.array(unbinned_measured_entries)");
//     TPython::Bind( &pass_reco, "pass_reco_mask" );
//     TPython::Exec( "pass_reco_mask = np.array(pass_reco_mask, dtype=bool)" );
//     TPython::Bind( &pass_truth, "pass_truth_mask" );
//     TPython::Exec( "pass_truth_mask = np.array(pass_truth_mask, dtype=bool)" );
//     TParameter<Int_t>* num_iterations_object = new TParameter<Int_t>("num_iterations", num_iterations);
//     TPython::Bind(num_iterations_object, "num_iterations" );

//     TPython::Exec("weights_test, MC_test, sim_test, pass_reco_test = unbinned_omnifold(unbinned_MC_entries,\
//                                                                                        unbinned_sim_entries,\
//                                                                                        unbinned_measured_entries,\
//                                                                                        pass_reco_mask,\
//                                                                                        pass_truth_mask,\
//                                                                                        num_iterations.GetVal())");

//     TPython::Exec("weights_MC_TVectorD = convert_to_TVectorD(weights_test[-1, 1])");
//     TPython::Exec("weights_sim_TVectorD = convert_to_TVectorD(weights_test[-1, 0][pass_reco_test])");

//     TPython::Exec("MC_vectors = get_vectors(MC_test)");
//     TPython::Exec("MC_TObjArray = ROOT.TObjArray(len(MC_vectors))\n"
//                   "for i in range(len(MC_vectors)):\n"
//                   "    MC_TObjArray[i] = MC_vectors[i]");
//     TObjArray out_MC_test = *(TObjArray*) TPython::Eval("MC_TObjArray");

//     TPython::Exec("sim_vectors = get_vectors(sim_test)");
//     TPython::Exec("sim_TObjArray = ROOT.TObjArray(len(sim_vectors))\n"
//                   "for i in range(len(sim_vectors)):\n"
//                   "    sim_TObjArray[i] = sim_vectors[i]");
//     TObjArray out_sim_test = *(TObjArray*) TPython::Eval("sim_TObjArray");
    
//     TVectorD out_weights_MC_test = *(TVectorD*) TPython::Eval("weights_MC_TVectorD");
//     TVectorD out_weights_sim_test = *(TVectorD*) TPython::Eval("weights_sim_TVectorD");

//     std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> out_tuple = std::make_tuple(out_weights_MC_test,
//                                                                                    out_MC_test,
//                                                                                    out_weights_sim_test,
//                                                                                    out_sim_test);
//     return out_tuple;
// }

// // Does unbinned unfolding with std::vectors as input and output
// // Converts the std::vectors to TObjects and then does unfolding with the TObjArray unfolding function
// std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>>
//                                                RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(std::vector<std::vector<Double_t>> MC_entries,
//                                                                           std::vector<std::vector<Double_t>> sim_entries,
//                                                                           std::vector<std::vector<Double_t>> measured_entries,
//                                                                           std::vector<Bool_t> pass_reco,
//                                                                           std::vector<Bool_t> pass_truth,
//                                                                           Int_t num_iterations)
// {
//     // Converting the std::vectors to TObjArrays and TVectors
//     // std::vector<std::vector<Double_t>> to a TObjArray filled with TVectorDs
//     auto nested_vector_to_TObjArray = [](std::vector<std::vector<Double_t>> vec)
//     {
//         TObjArray out_TObjArray(vec.size());
//         for (size_t i = 0; i < vec.size(); i++)
//         {
//             TVectorD *subvec = new TVectorD(vec[i].size());
//             for(size_t j = 0; j < vec[i].size(); j++)
//             {
//                 (*subvec)[j] = vec[i][j];
//             }
//             out_TObjArray[i] = subvec;
//         }
//         return out_TObjArray;
//     };

//     auto TObjArray_to_nested_vector = [](TObjArray objarray)
//     {
//         std::vector<std::vector<Double_t>> nested_vector(objarray.GetEntries());
//         for (int i = 0; i < objarray.GetEntries(); i++)
//         {
//             TVectorD subTVectorD = *(TVectorD*) objarray[i];
//             std::vector<Double_t> subvec(subTVectorD.GetNoElements());
//             for(int j = 0; j < subTVectorD.GetNoElements(); j++)
//             {
//                 subvec[j] = subTVectorD[j];
//             }
//             nested_vector[i] = subvec;
//         }
//         return nested_vector;
//     };

//     auto MC_TObjArray = nested_vector_to_TObjArray(MC_entries);
//     auto sim_TObjArray = nested_vector_to_TObjArray(sim_entries);
//     auto measured_TObjArray = nested_vector_to_TObjArray(measured_entries);

//     TVector pass_reco_TVector(pass_reco.size());
//     for (size_t i = 0; i < pass_reco.size(); i++) {
//         pass_reco_TVector[i] = pass_reco[i];
//     }

//     TVector pass_truth_TVector(pass_truth.size());
//     for (size_t i = 0; i < pass_truth.size(); i++) {
//         pass_truth_TVector[i] = pass_truth[i];
//     }

//     // Running unbinned unfolding
//     std::tuple<TVectorD, TObjArray, TVectorD, TObjArray> out_tuple_TObject = this->UnbinnedOmnifold(MC_TObjArray,
//                                                                                                     sim_TObjArray,
//                                                                                                     measured_TObjArray,
//                                                                                                     pass_reco_TVector,
//                                                                                                     pass_truth_TVector,
//                                                                                                     num_iterations);

//     TVectorD weights_MC_TVectorD = std::get<0>(out_tuple_TObject);
//     TVectorD weights_sim_TVectorD = std::get<2>(out_tuple_TObject);

//     // Converting TObjects back to std::vectors
//     auto MC_vector = TObjArray_to_nested_vector(std::get<1>(out_tuple_TObject));
//     auto sim_vector = TObjArray_to_nested_vector(std::get<3>(out_tuple_TObject));
    
//     std::vector<Double_t> weights_MC_vector(weights_MC_TVectorD.GetNoElements());
//     std::vector<Double_t> weights_sim_vector(weights_sim_TVectorD.GetNoElements());
//     for(int i = 0; i < weights_MC_TVectorD.GetNoElements(); i++)
//     {
//         weights_MC_vector[i] = weights_MC_TVectorD[i];
//     }

//     for(int i = 0; i < weights_sim_TVectorD.GetNoElements(); i++)
//     {
//         weights_sim_vector[i] = weights_sim_TVectorD[i];
//     }

//     for (int i = 0; i < MC_TObjArray.GetEntries(); i++) {
//         delete MC_TObjArray[i];
//     }
//     for (int i = 0; i < sim_TObjArray.GetEntries(); i++) {
//         delete sim_TObjArray[i];
//     }
//     for (int i = 0; i < measured_TObjArray.GetEntries(); i++) {
//         delete measured_TObjArray[i];
//     }

//     std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> out_tuple_vector = std::make_tuple(weights_MC_vector,
//                                                                                                                                                                  MC_vector,
//                                                                                                                                                                  weights_sim_vector,
//                                                                                                                                                                  sim_vector);

//     return out_tuple_vector;                                                                                                                                                                   
// }

// // 1D unbinned unfolding using std::vectors
// // Same as the TVectorD unbinned unfolding function with std::vectors instead of TVectorD
// std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> 
//                                                          RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(std::vector<Double_t> MC_entries,
//                                                                           std::vector<Double_t> sim_entries,
//                                                                           std::vector<Double_t> measured_entries,
//                                                                           std::vector<Bool_t> pass_reco,
//                                                                           std::vector<Bool_t> pass_truth,
//                                                                           Int_t num_iterations)
// {

//     auto nest_vector = [](std::vector< Double_t> vec)
//     {
//         std::vector<std::vector<Double_t>> nested_vector;
//         for(const Double_t& entry: vec)
//         {
//             std::vector<Double_t> dummy_vector{entry};
//             nested_vector.push_back(dummy_vector);
//         }
//         return nested_vector;
//     };

//     auto unnest_vector = [](std::vector<std::vector<Double_t>> vec)
//     {
//         std::vector<Double_t> unnested_vector;
//         for(const std::vector<Double_t>& entry: vec)
//         {
//             unnested_vector.push_back(entry[0]);
//         }
//         return unnested_vector;
//     };
    
//     auto nested_MC_entries = nest_vector(MC_entries);
//     auto nested_sim_entries = nest_vector(sim_entries);
//     auto nested_measured_entries = nest_vector(measured_entries);

//     std::tuple<std::vector<Double_t>, std::vector<std::vector<Double_t>>, std::vector<Double_t>, std::vector<std::vector<Double_t>>> nested_out_tuple = this->UnbinnedOmnifold(nested_MC_entries,
//                                                                                                                                                                                nested_sim_entries,
//                                                                                                                                                                                nested_measured_entries,
//                                                                                                                                                                                pass_reco,
//                                                                                                                                                                                pass_truth,
//                                                                                                                                                                                num_iterations);

//     auto unnested_out_MC_entries = unnest_vector(std::get<1>(nested_out_tuple));
//     auto unnested_out_sim_entries = unnest_vector(std::get<3>(nested_out_tuple));
//     std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> out_tuple = std::make_tuple(std::get<0>(nested_out_tuple),
//                                                                                                                                        unnested_out_MC_entries,
//                                                                                                                                        std::get<2>(nested_out_tuple),
//                                                                                                                                        unnested_out_sim_entries);

//     return out_tuple;                                                                                                                                       
// }
template class RooUnfoldOmnifoldT<TH1,TH2>;
ClassImp (RooUnfoldOmnifold)

#ifndef NOROOFIT
template class RooUnfoldOmnifoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldOmnifold)
#endif

                                                                       