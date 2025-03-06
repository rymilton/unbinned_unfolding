from omnifold import DataLoader, MultiFold, MLP
import h5py as h5
import argparse 
import os 
import numpy as np
import horovod.tensorflow.keras as hvd
import time
import gzip
import pickle
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a PET model using Pythia and Herwig data.")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--model_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/final_weights", help="Folder where the model is saved")
    parser.add_argument("--output_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store predictions ")
    args = parser.parse_args()
    return args
def expit(x):
    return 1. / (1. + np.exp(-x))
def reweight(events,model,batch_size=None):
    f = expit(model.predict(events,batch_size=batch_size))
    weights = f / (1. - f)  # this is the crux of the reweight, approximates likelihood ratio
    weights = np.nan_to_num(weights[:,0],posinf=1)
    return weights

def main():
    flags = parse_arguments()
    data_dir = flags.data_dir
    gaussian_file = data_dir + "gaussian_data.npz"

    gaussian_data = np.load(gaussian_file)

    unbinned_MC_data = gaussian_data["MC"]
    unbinned_sim_data = gaussian_data["sim"]
    unbinned_truth_data = gaussian_data["truth"]
    measured_train = gaussian_data["measured"]
    sim_pass_reco = gaussian_data["sim_pass_reco"]
    measured_pass_reco = gaussian_data["measured_pass_reco"]
    
    if len(unbinned_MC_data) != len(unbinned_sim_data):
        raise ValueError("MC and sim data have different number of entries!")
    
    if len(unbinned_MC_data) != len(sim_pass_reco):
        raise ValueError("MC and sim pass reco have different number of entries!")
    
    if len(unbinned_truth_data) != len(measured_train):
        raise ValueError("Truth and measured data have different number of entries!")
    
    if len(unbinned_truth_data) != len(measured_pass_reco):
        raise ValueError("Truth and measured pass reco have different number of entries!")
    
    gen_model = MLP(1)
    gen_model.load_weights(flags.model_dir + "OmniFold_test_iter1_step2.weights.h5")

    MC_test = np.expand_dims(MC_test, axis=1)
    unfolded_weights = reweight(MC_test,gen_model)
    normalized_MC_test = (MC_test - np.mean(MC_test, axis=0))/np.std(MC_test, axis=0)
    unfolded_weights_normalized = reweight(normalized_MC_test,gen_model)

    # Combine the variables into a dictionary to save them in a single pickle file
    data_to_save = {
        "DNN_weights": unfolded_weights,
        "DNN_weights_normalized":unfolded_weights_normalized
    }

    # Saving the combined data as a compressed pickle.gz file
    output_file = f'{flags.output_dir}/saved_DNN_gaussian_weights_5iters.pickle.gz'

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(data_to_save, f)


if __name__ == '__main__':
    main()
