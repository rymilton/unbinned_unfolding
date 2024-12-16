import ROOT
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.model_selection import train_test_split
import pickle
import os

def TH1_to_numpy(hist):
    num_bins = hist.GetNbinsX()
    hist_counts = np.empty(num_bins)
    bin_centers = np.empty(num_bins)
    for i in range(num_bins):
        hist_counts[i] = hist.GetBinContent(i+1)
        bin_centers[i] = hist.GetBinCenter(i+1)
    return hist_counts, bin_centers

def TH2_to_numpy(hist):
    num_X_bins = hist.GetNbinsX()
    num_Y_bins = hist.GetNbinsY()
    hist_counts = np.empty(shape=(num_X_bins, num_Y_bins))
    bin_centers = np.empty(shape=(num_X_bins, num_Y_bins), dtype=object)
    for i in range(num_X_bins):
        for j in range(num_Y_bins):
            hist_counts[i, j] = hist.GetBinContent(i+1, j+1)
            bin_center_tuple = (hist.GetXaxis().GetBinCenter(i+1), hist.GetYaxis().GetBinCenter(j+1))
            bin_centers[i, j] = bin_center_tuple
    return hist_counts, bin_centers

def prepare_hist_data(counts, bin_centers):
    out_array = np.empty(shape=(int(np.sum(counts)), 1))
    entry_tracker = 0
    for (count, bin_center) in zip(counts, bin_centers):
        out_array[entry_tracker:int(entry_tracker+count)] = bin_center
        entry_tracker += int(count)
    return out_array

def prepare_response_data(counts, bin_centers):
    truth_array = np.empty(shape=(int(np.sum(counts)), 1))
    sim_array = np.empty(shape=(int(np.sum(counts)), 1))
    entry_tracker = 0
    for (count, bin_center) in zip(counts, bin_centers):
        sim_array[entry_tracker:int(entry_tracker+count)] = bin_center[0]
        truth_array[entry_tracker:int(entry_tracker+count)] = bin_center[1]
        entry_tracker += int(count)
    return truth_array, sim_array

def convert_to_TVectorD(array):
    vector = ROOT.TVectorD(len(array))
    for i, entry in enumerate(array):
        vector[i] = entry
    return vector

def get_vectors(array):
    vector_list = []
    for entry in array:
        vector = convert_to_TVectorD(entry)
        vector_list.append(vector)
    return vector_list

def reweight(events, classifier):
    class_probabilities = classifier.predict_proba(events)
    data_probability = class_probabilities[:,1]
    weights = data_probability / (1. - data_probability)
    return np.squeeze(np.nan_to_num(weights))

def omnifold(MC_entries, sim_entries, measured_entries, MC_pass_reco_mask, MC_pass_truth_mask, measured_pass_truth_mask, num_iterations, model_save_dict = None):
    # Removing events that don't pass generation level cuts
    sim_entries = sim_entries[MC_pass_truth_mask]
    MC_entries = MC_entries[MC_pass_truth_mask]
    MC_pass_reco_mask = MC_pass_reco_mask[MC_pass_truth_mask]
    
    MC_train, sim_train, pass_reco_train = MC_entries, sim_entries, MC_pass_reco_mask
    
    measured_labels = np.ones(len(measured_entries[measured_pass_truth_mask]))
    MC_labels = np.zeros(len(MC_train))

    weights_pull_train = np.ones(len(MC_train))
    weights_push_train = np.ones(len(MC_train))
    weights_train = np.empty(shape=(num_iterations, 2, len(MC_train)))

    step1_classifier = GradientBoostingClassifier()
    step2_classifier = GradientBoostingClassifier()
    use_regressor =  any(~pass_reco_train)
    if use_regressor:
        step1_regressor = GradientBoostingRegressor()
    for i in range(num_iterations):
        print(f"Starting iteration {i}") 
        step1_data = np.concatenate((sim_train[pass_reco_train], measured_entries[measured_pass_truth_mask]))
        step1_labels = np.concatenate((np.zeros(len(sim_train[pass_reco_train])), measured_labels))
        step1_weights = np.concatenate((weights_push_train[pass_reco_train], np.ones(len(measured_entries[measured_pass_truth_mask]))))
        
        # Training step 1 classifier and getting weights
        step1_classifier.fit(step1_data, step1_labels, sample_weight = step1_weights)
        new_weights_train = np.ones_like(weights_pull_train)
        new_weights_train[pass_reco_train] = reweight(sim_train[pass_reco_train], step1_classifier)
        
        # Training a regression model to predict the weights of the events that don't pass reconstruction
        if use_regressor:
            step1_regressor.fit(MC_train[pass_reco_train], new_weights_train[pass_reco_train])
            new_weights_train[~pass_reco_train] = step1_regressor.predict(MC_train[~pass_reco_train])
        weights_pull_train = np.multiply(weights_push_train, new_weights_train)
            
        # Training step 2 classifier
        step2_data = np.concatenate((MC_train, MC_train))
        step2_labels = np.concatenate((MC_labels, np.ones(len(MC_train))))
        step2_weights = np.concatenate((np.ones(len(MC_train)), weights_pull_train))
        step2_classifier.fit(step2_data, step2_labels, sample_weight = step2_weights)

        # Getting step 2 weights and storing iteration weights
        weights_push_train = reweight(MC_train, step2_classifier)
        weights_train[i, 0], weights_train[i, 1] = weights_pull_train, weights_push_train

    # Saving the models if the user wants to (saved by default)
    if model_save_dict is not None and model_save_dict['save_models']:
        base_path = model_save_dict['save_dir']
        os.makedirs(os.path.dirname(base_path ), exist_ok=True)
        models = {"step1_classifier": step1_classifier,
                  "step2_classifier": step2_classifier
                }
        if use_regressor:
            models['step1_regressor'] = step1_regressor
        with open(f"{model_save_dict['model_name']}_models.pkl", "wb") as outfile:
            pickle.dump(models, outfile)

    return weights_pull_train, weights_push_train


def binned_omnifold(response_hist, measured_hist, num_iterations):
    measured_counts, measured_bin_centers = TH1_to_numpy(measured_hist)
    response_counts, response_bin_centers = TH2_to_numpy(response_hist)
    MC_entries, sim_entries = prepare_response_data(response_counts.flatten(), response_bin_centers.flatten())
    measured_entries = prepare_hist_data(measured_counts, measured_bin_centers)
    MC_pass_reco_mask = np.full(MC_entries.shape[0], True)
    MC_pass_truth_mask = np.full(sim_entries.shape[0], True)
    MC_pass_measured_mask = np.full(measured_entries.shape[0], True)
    _, step2_weights = omnifold(MC_entries, sim_entries, measured_entries, MC_pass_reco_mask, MC_pass_truth_mask, MC_pass_measured_mask, num_iterations)
    unfolded_hist = ROOT.TH1D("unfolded_hist",
                              "unfolded_hist",
                              response_hist.GetNbinsY(),
                              response_hist.GetYaxis().GetBinLowEdge(1),
                              response_hist.GetYaxis().GetBinLowEdge(response_hist.GetNbinsY())
                              +response_hist.GetYaxis().GetBinWidth(response_hist.GetNbinsY())
                             )
    for (weight, MC) in zip(step2_weights, MC_entries.flatten()):
        unfolded_hist.Fill(MC, weight)
    return unfolded_hist
def unbinned_omnifold(MC_entries, sim_entries, measured_entries, num_iterations, MC_pass_reco_mask = None, MC_pass_truth_mask = None, measured_pass_reco_mask = None, model_save_dict = None):
    if MC_entries.ndim == 1:
        MC_entries = np.expand_dims(MC_entries, axis = 1)
    if sim_entries.ndim == 1:
        sim_entries = np.expand_dims(sim_entries, axis = 1)
    if measured_entries.ndim == 1:
        measured_entries = np.expand_dims(measured_entries, axis = 1)
    if MC_pass_reco_mask is None:
        MC_pass_reco_mask = np.full(MC_entries.shape[0], True, dtype=bool)
    if MC_pass_truth_mask is None:
        MC_pass_truth_mask = np.full(MC_entries.shape[0], True, dtype=bool)
    if measured_pass_reco_mask is None:
        measured_pass_reco_mask = np.full(measured_entries.shape[0], True, dtype=bool)
    return omnifold(MC_entries, sim_entries, measured_entries, MC_pass_reco_mask, MC_pass_truth_mask, measured_pass_reco_mask, num_iterations, model_save_dict)

def get_step1_predictions(MC_data, sim_data, model_info_dict, pass_reco = None):
    with open(f"{model_info_dict['model_name']}_models.pkl", "rb") as infile:
        loaded_models = pickle.load(infile)
    step1_test_weights = np.ones(len(sim_data))
    step1_test_weights[pass_reco] = reweight(sim_data[pass_reco], loaded_models['step1_classifier'])
    step1_test_weights[~pass_reco] = loaded_models['step1_regressor'].predict(MC_data[~pass_reco])
    return step1_test_weights

def get_step2_predictions(MC_data, model_info_dict):
    with open(f"{model_info_dict['model_name']}_models.pkl", "rb") as infile:
        loaded_models = pickle.load(infile)
    step2_test_weights = reweight(MC_data, loaded_models['step2_classifier'])
    return step2_test_weights