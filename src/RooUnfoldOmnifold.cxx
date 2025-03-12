/*
    Implementation of Omnifold using boosted decision trees (BDT)
    This follows the original Omnifold work (Andreassen et al. PRL 124, 182001 (2020))

    This class has BinnedOmnifold (binned unfolding) and UnbinnedOmnifold (unbinned unfolding).

    BinnedOmnifold is similar to RooUnfoldBayes, where it uses the response matrix and measured histogram to do unfolding 
    and it returns the unfolded TH1*.
    UnbinnedOmnifold instead takes a RDataFrames of Monte Carlo data, reconstructed Monte Carlo data, and
    measured data. It also takes masks (in forms of TVectors) that indicate an event passes generation level cuts and reconstruction.
    
    To do the unfolding, the input data is first moved to Python using TPython. This conversion needs TObjects, which
    is why RDataFrames are used in UnbinnedOmnifold. The unfolding procedure is done in Python and the results
    are brought back to C++.

    BDTs are trained to output sets of weights that are applied to the data. 
    In BinnedOmnifold, these weights are applied within the function and an unfolded histogram is returned.
    In UnbinnedOmnifold, these weights are returned within an std::tuple of the form <step 1 weights, step 2 weights>,
    where each set of weights is a TVector.
    
    To plot the UnbinnedOmnifold results in a histogram, the Monte Carlo data at the gen level should be filled with
    the step 2 weights. The resulting histogram is the unfolded distribution. The Monte Carlo at reco level should be filled 
    with the step 1 weights. For the weights from the first iteration, this should closely match the measured distribution.
    The weights are event-by-event weights, so even if you have multi-dimensional data, there is only one weight for a given event.
    e.g. MC data = [[1, 2], [1.2, 2.5]], MC weights = [.5, .2]. The .5 (.2) is applied to both 1 and 2 (1.2 and 2.5) from [1, 2] ([1.2, 2.5]).

    The parameters of the BDTs can be adjusted using TMaps and the SetParameter functions (SetStep1ClassifierParameters,
    SetStep2ClassifierParameters, and SetStep1RegressorParameters).

    The trained models will be saved by default and they can be used to predict the weights of a new dataset using
    TestUnbinnedOmnifold.

    A notebook showing how to run this class is in the examples directory (examples/OmniFold_example.ipynb)
*/
#include <iostream>
#include "RooUnfoldHelpers.h"
#include "RooUnfoldTH1Helpers.h"
#include "RooUnfoldOmnifold.h"
#include "TParameter.h"
#include "TH1.h"
#include "TPython.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <ROOT/RDataFrame.hxx>
#include <TMap.h>
#include <TObjString.h>
using namespace RooUnfolding;

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT()
  : RooUnfoldT<Hist,Hist2D>(), _useDensity(0), _MCgenDataFrame(0), _MCrecoDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0),
    _MCPassTruth(0), _MeasuredPassReco(0),_unbinned_step1_weights(0), _unbinned_step2_weights(0), _SaveUnbinnedModels(true),
    _UnbinnedModelSaveDir("./"), _UnbinnedModelName("RooUnfoldOmnifold"), _TestMCgenDataFrame(0), _TestMCrecoDataFrame(0), _TestMCPassReco(0),
    _Step1ClassifierParameters(0), _Step2ClassifierParameters(0), _Step1RegressorParameters(0), _MCgenWeights(0), _MCrecoWeights(0), _MeasuredWeights(0)
{

  //! Default constructor. Use Setup() to prepare for unfolding.]
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D>(name,title), _useDensity(0), _MCgenDataFrame(0), _MCrecoDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0),
    _MCPassTruth(0), _MeasuredPassReco(0),_unbinned_step1_weights(0), _unbinned_step2_weights(0), _SaveUnbinnedModels(true),
    _UnbinnedModelSaveDir("./"), _UnbinnedModelName("RooUnfoldOmnifold"), _TestMCgenDataFrame(0), _TestMCrecoDataFrame(0), _TestMCPassReco(0),
    _Step1ClassifierParameters(0), _Step2ClassifierParameters(0), _Step1RegressorParameters(0), _MCgenWeights(0), _MCrecoWeights(0), _MeasuredWeights(0)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const TString& name, const TString& title)
  : RooUnfoldT<Hist,Hist2D>(name,title), _useDensity(0), _MCgenDataFrame(0), _MCrecoDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0),
    _MCPassTruth(0), _MeasuredPassReco(0),_unbinned_step1_weights(0), _unbinned_step2_weights(0), _SaveUnbinnedModels(true),
    _UnbinnedModelSaveDir("./"), _UnbinnedModelName("RooUnfoldOmnifold"), _TestMCgenDataFrame(0), _TestMCrecoDataFrame(0), _TestMCPassReco(0),
    _Step1ClassifierParameters(0), _Step2ClassifierParameters(0), _Step1RegressorParameters(0), _MCgenWeights(0), _MCrecoWeights(0), _MeasuredWeights(0)
{
  //! Basic named constructor. Use Setup() to prepare for unfolding.
  Init();
}

template<class Hist,class Hist2D>
RooUnfoldOmnifoldT<Hist,Hist2D>::RooUnfoldOmnifoldT (const RooUnfoldResponseT<Hist,Hist2D>* res, const Hist* meas, Int_t niter, bool useDensity,
                        const char* name, const char* title)
  : RooUnfoldT<Hist,Hist2D> (res, meas, name, title), _useDensity(useDensity), _niter(niter), _MCgenDataFrame(0), _MCrecoDataFrame(0), _MeasuredDataFrame(0), _MCPassReco(0),
    _MCPassTruth(0), _MeasuredPassReco(0),_unbinned_step1_weights(0), _unbinned_step2_weights(0), _SaveUnbinnedModels(true),
    _UnbinnedModelSaveDir("./"), _UnbinnedModelName("RooUnfoldOmnifold"), _TestMCgenDataFrame(0), _TestMCrecoDataFrame(0), _TestMCPassReco(0),
    _Step1ClassifierParameters(0), _Step2ClassifierParameters(0), _Step1RegressorParameters(0), _MCgenWeights(0), _MCrecoWeights(0), _MeasuredWeights(0)
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

  // Getting the directory of the RooUnfold installation so we can load the omnifold.py script
  auto dynamic_paths = string_split(gSystem->GetDynamicPath(), ":", false);
  TString RooUnfold_install_path;
  for(auto &path : dynamic_paths)
  {
    TString RooUnfold_library_name = "libRooUnfold.so";
    if(gSystem->FindFile(path.c_str(), RooUnfold_library_name))
    {
      if(!RooUnfold_install_path.IsNull())
      {
        std::cerr<<"WARNING: Multiple RooUnfold installations detected in PATH"<<"\n"
        "Previously found "<<RooUnfold_install_path<<". Will now use "<<path<<" for OmniFold."<<std::endl;
      }
      RooUnfold_install_path = path;
    }
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
      TPython::Exec(Form("exec(open('%s').read())", omnifold_path.Data()));
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
  
  const Hist2D* res = response->Hresponse();
  TVectorD vals(h2v<Hist>(res,false, false));
  TVectorD errs(h2ve<Hist>(res,false,false));
  TH2D *response_hist = (TH2D*) createHist<Hist2D>(vals,errs,name(res),title(res),vars(res),false);
  TH1D *measured_hist = (TH1D*) this->Hmeasured();
  // Sending histograms to Python
  TPython::Bind( response_hist, "response_hist" ); 
  TPython::Bind( measured_hist, "measured_hist" );
  // Converting num_iterations to a TObject* to convert to Python
  TPython::Exec(Form("num_iterations = %d", _niter));
  TPython::Exec(Form("binned_use_density = %d", this->_useDensity));
  TPython::Exec("binned_use_density = bool(binned_use_density)");
  // Performing binned Omnifold and bringing histogram back to ROOT and correcting for efficiency
  TH1D *unfolded_hist = TPython::Eval("OmniFold_helper_functions.binned_omnifold(response_hist, measured_hist, num_iterations, binned_use_density)");
  
  auto efficiency = response->Vefficiency();
  for (Int_t i = 0; i < unfolded_hist->GetNbinsX(); i++)
  {
      Double_t corrected_content = efficiency[i] > 0 ?
          unfolded_hist->GetBinContent(i+1)/efficiency[i] : unfolded_hist->GetBinContent(i+1);
      unfolded_hist->SetBinContent(i+1, corrected_content);
  }
  
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
      TVectorD vector(*(df.Count()));
      if (column_type == "double" || column_type == "Double_t")
      {
        auto column_values = *(df.Take<double>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("column_vector = np.asarray(column_vector)");
        TPython::Exec("data_dict[data_type] = column_vector.reshape(-1, 1) if data_dict[data_type].size == 0 else np.hstack([data_dict[data_type], column_vector.reshape(-1, 1)])");
      }
      else if (column_type == "float" || column_type == "Float_t")
      {
        auto column_values = *(df.Take<float>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("column_vector = np.asarray(column_vector)");
        TPython::Exec("data_dict[data_type] = column_vector.reshape(-1, 1) if data_dict[data_type].size == 0 else np.hstack([data_dict[data_type], column_vector.reshape(-1, 1)])");
      }
      else
        throw std::runtime_error("Please make sure RDataFrame columns are type double/Double_t or float/Float_t");
    }
  };
  TPython::Exec("data_dict = {}");
  DataFrame_to_python(this->_MCgenDataFrame, "MCgen");
  DataFrame_to_python(this->_MCrecoDataFrame, "MCreco");
  DataFrame_to_python(this->_MeasuredDataFrame, "measured");

  TPython::Exec(Form("num_iterations = %d", this->_niter));
  if(this->_MCPassReco.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MCPassReco), "MC_pass_reco");
    TPython::Exec("MC_pass_reco = np.array(MC_pass_reco, dtype=bool)");
  }
  else
    TPython::Exec("MC_pass_reco = None");
  
  if(this->_MCPassTruth.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MCPassTruth), "MC_pass_truth");
    TPython::Exec("MC_pass_truth = np.array(MC_pass_truth, dtype=bool)");
  }
  else
    TPython::Exec("MC_pass_truth = None");

  if(this->_MeasuredPassReco.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MeasuredPassReco), "measured_pass_reco");
    TPython::Exec("measured_pass_reco = np.array(measured_pass_reco, dtype=bool)");
  }
  else
    TPython::Exec("measured_pass_reco = None");

  if(this->_MCgenWeights.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MCgenWeights), "MCgen_weights");
    TPython::Exec("MCgen_weights = np.array(MCgen_weights, dtype=float)");
  }
  else
    TPython::Exec("MCgen_weights = None");
  if(this->_MCrecoWeights.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MCrecoWeights), "MCreco_weights");
    TPython::Exec("MCreco_weights = np.array(MCreco_weights, dtype=float)");
  }
  else
    TPython::Exec("MCreco_weights = None");
  if(this->_MeasuredWeights.GetNoElements() != 0)
  {
    TPython::Bind(&(this->_MeasuredWeights), "measured_weights");
    TPython::Exec("measured_weights = np.array(measured_weights, dtype=float)");
  }
  else
    TPython::Exec("measured_weights = None");

  // Moving save model info to Python
  TPython::Exec(Form("save_models = %d", this->_SaveUnbinnedModels? 1 : 0));
  TPython::Exec("save_models = True if save_models == 1 else False");
  TPython::Exec(Form("model_save_dir = '%s'", this->_UnbinnedModelSaveDir.Data()));
  TPython::Exec(Form("model_name = '%s'", this->_UnbinnedModelName.Data()));
  TPython::Exec("model_save_dict = {'save_models':save_models, 'save_dir':model_save_dir, 'model_name':model_name}");
  
  // Moving model parameter info to Python
  if(this->_Step1ClassifierParameters)
  {
    TPython::Bind(this->_Step1ClassifierParameters, "step1classifier_params_map");
    TString dict_script = Form(
      "step1classifier_params = {}\n"
      "iter = step1classifier_params_map.MakeIterator()\n"
      "num_keys = step1classifier_params_map.GetSize()\n"
      "counter = 0\n"
      "while counter < num_keys:\n"
      "\tkey = iter.Next()\n"
      "\tstep1classifier_params[str(key)]=str(step1classifier_params_map.GetValue(key))\n"
      "\tcounter += 1\n"
    );
    TPython::Exec(dict_script.Data());
  }
  else 
  {
    TPython::Exec("step1classifier_params = {}");
  }

  if(this->_Step2ClassifierParameters)
  {
    TPython::Bind(this->_Step2ClassifierParameters, "step2classifier_params_map");
    TString dict_script = Form(
      "step2classifier_params = {}\n"
      "iter = step2classifier_params_map.MakeIterator()\n"
      "num_keys = step2classifier_params_map.GetSize()\n"
      "counter = 0\n"
      "while counter < num_keys:\n"
      "\tkey = iter.Next()\n"
      "\tstep2classifier_params[str(key)]=str(step2classifier_params_map.GetValue(key))\n"
      "\tcounter += 1\n"
    );
    TPython::Exec(dict_script.Data());
  }
  else 
  {
    TPython::Exec("step2classifier_params = {}");
  }

  if(this->_Step1RegressorParameters)
  {
    TPython::Bind(this->_Step1RegressorParameters, "step1regressor_params_map");
    TString dict_script = Form(
      "step1regressor_params = {}\n"
      "iter = step1regressor_params_map.MakeIterator()\n"
      "num_keys = step1regressor_params_map.GetSize()\n"
      "counter = 0\n"
      "while counter < num_keys:\n"
      "\tkey = iter.Next()\n"
      "\tstep1regressor_params[str(key)]=str(step1regressor_params_map.GetValue(key))\n"
      "\tcounter += 1\n"
    );
    TPython::Exec(dict_script.Data());
  }
  else 
  {
    TPython::Exec("step1regressor_params = {}");
  }
	
    // Execute the Python script
  TPython::Exec("step1_weights, step2_weights = OmniFold_helper_functions.unbinned_omnifold(data_dict['MCgen'],\
                                                                   data_dict['MCreco'],\
                                                                   data_dict['measured'],\
                                                                   num_iterations,\
                                                                   MC_pass_reco,\
                                                                   MC_pass_truth,\
                                                                   measured_pass_reco,\
                                                                   MCgen_weights=MCgen_weights,\
                                                                   MCreco_weights=MCreco_weights,\
                                                                   measured_weights=measured_weights,\
                                                                   model_save_dict=model_save_dict,\
                                                                   classifier1_params=step1classifier_params,\
                                                                   classifier2_params=step2classifier_params,\
                                                                   regressor_params=step1regressor_params)");

    
  TVectorD step1_weights = *(std::unique_ptr<TVectorD>) TPython::Eval("OmniFold_helper_functions.convert_to_TVectorD(step1_weights)");
  TVectorD step2_weights = *(std::unique_ptr<TVectorD>) TPython::Eval("OmniFold_helper_functions.convert_to_TVectorD(step2_weights)");

  this->_unbinned_step1_weights.ResizeTo(step1_weights);
  this->_unbinned_step1_weights = step1_weights;
  this->_unbinned_step2_weights.ResizeTo(step2_weights);
  this->_unbinned_step2_weights = step2_weights;

  std::tuple<TVectorD, TVectorD> out_tuple = std::make_tuple(this->_unbinned_step1_weights, this->_unbinned_step2_weights);
  return out_tuple;
}

template<class Hist,class Hist2D> std::tuple<TVectorD, TVectorD>
RooUnfoldOmnifoldT<Hist,Hist2D>::UnbinnedOmnifold(ROOT::RDataFrame MCgen_dataframe,
                                                  ROOT::RDataFrame MCreco_dataframe,
                                                  ROOT::RDataFrame measured_dataframe)
{
  this->SetMCgenDataFrame(MCgen_dataframe);
  this->SetMCrecoDataFrame(MCreco_dataframe);
  this->SetMeasuredDataFrame(measured_dataframe);
  return this->UnbinnedOmnifold();
}

template<class Hist,class Hist2D> std::tuple<TVectorD, TVectorD>
RooUnfoldOmnifoldT<Hist,Hist2D>::TestUnbinnedOmnifold()
{
  auto DataFrame_to_python = [](ROOT::RDataFrame df, TString data_type)
  {
    TPython::Exec(Form("data_type = '%s'", data_type.Data()));
    TPython::Exec("data_dict[data_type] = np.array([])");
    for(const auto& column_name : df.GetColumnNames())
    {
      std::string column_type = df.GetColumnType(column_name);
      TVectorD vector(*(df.Count()));
      if (column_type == "double" || column_type == "Double_t")
      {
        auto column_values = *(df.Take<double>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("column_vector = np.asarray(column_vector)");
        TPython::Exec("data_dict[data_type] = column_vector.reshape(-1, 1) if data_dict[data_type].size == 0 else np.hstack([data_dict[data_type], column_vector.reshape(-1, 1)])");
      }
      else if (column_type == "float" || column_type == "Float_t")
      {
        auto column_values = *(df.Take<float>(column_name));
        for (size_t i = 0; i < column_values.size(); i++)
        {
          vector[i] = column_values[i];
        }
        TPython::Bind(&vector, "column_vector");
        TPython::Exec("column_vector = np.asarray(column_vector)");
        TPython::Exec("data_dict[data_type] = column_vector.reshape(-1, 1) if data_dict[data_type].size == 0 else np.hstack([data_dict[data_type], column_vector.reshape(-1, 1)])");
      }
      else
        throw std::runtime_error("Please make sure RDataFrame columns are type double/Double_t or float/Float_t");
    }
  };
  TPython::Exec("data_dict = {}");
  DataFrame_to_python(this->_TestMCgenDataFrame, "test_MCgen");
  DataFrame_to_python(this->_TestMCrecoDataFrame, "test_MCreco");

  TPython::Bind(&(this->_TestMCPassReco), "test_MC_pass_reco");
  TPython::Exec("test_MC_pass_reco = np.array(test_MC_pass_reco, dtype=bool)");

  TPython::Exec(Form("model_save_dir = '%s'", this->_UnbinnedModelSaveDir.Data()));
  TPython::Exec(Form("model_name = '%s'", this->_UnbinnedModelName.Data()));
  TPython::Exec("model_info_dict = {'save_dir':model_save_dir, 'model_name':model_name}");
  TPython::Exec("step1_test_weights = OmniFold_helper_functions.get_step1_predictions(data_dict['test_MCgen'],\
                                                            data_dict['test_MCreco'],\
                                                            model_info_dict,\
                                                            test_MC_pass_reco)");
  TPython::Exec("step2_test_weights = OmniFold_helper_functions.get_step2_predictions(data_dict['test_MCgen'], model_info_dict)");

  TVectorD step1_test_weights = *(std::unique_ptr<TVectorD>) TPython::Eval("OmniFold_helper_functions.convert_to_TVectorD(step1_test_weights)");
  TVectorD step2_test_weights = *(std::unique_ptr<TVectorD>) TPython::Eval("OmniFold_helper_functions.convert_to_TVectorD(step2_test_weights)");

  this->_unbinned_step1_test_weights.ResizeTo(step1_test_weights);
  this->_unbinned_step1_test_weights = step1_test_weights;
  this->_unbinned_step2_test_weights.ResizeTo(step2_test_weights);
  this->_unbinned_step2_test_weights = step2_test_weights;

  std::tuple<TVectorD, TVectorD> out_tuple = std::make_tuple(this->_unbinned_step1_test_weights, this->_unbinned_step2_test_weights);
  return out_tuple;
}

template class RooUnfoldOmnifoldT<TH1,TH2>;
ClassImp (RooUnfoldOmnifold)

#ifndef NOROOFIT
template class RooUnfoldOmnifoldT<RooUnfolding::RooFitHist,RooUnfolding::RooFitHist>;
ClassImp (RooFitUnfoldOmnifold)
#endif
