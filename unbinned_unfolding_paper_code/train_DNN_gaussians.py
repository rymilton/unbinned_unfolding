from omnifold import DataLoader, MultiFold, MLP
import h5py as h5
import argparse 
import numpy as np
import horovod.tensorflow.keras as hvd
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a DNN model using Gaussian data.")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--save_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store trained model weights")
    parser.add_argument("--num_train_data", type=int, default=1500000, help="Number of data to train with")
    parser.add_argument("--num_iterations", type=int, default=5, help="Number of iterations to use during training")
    args = parser.parse_args()
    return args

def main():
    flags = parse_arguments()
    num_train_data = flags.num_train_data
    itnum = flags.num_iterations
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
    
    
    MC_train = unbinned_MC_data[:num_train_data]
    sim_train = unbinned_sim_data[:num_train_data]
    pass_reco_train = sim_pass_reco[:num_train_data]
    measured_train = measured_train[:375000]
    measured_pass_reco_train = measured_pass_reco[:375000]
    MC_train = np.expand_dims(MC_train, axis=1)
    sim_train = np.expand_dims(sim_train, axis=1)
    measured_train = np.expand_dims(measured_train, axis=1)

    combined_step1_data = np.concatenate((sim_train[pass_reco_train], measured_train[measured_pass_reco_train]))
    normalized_MC_reco = (sim_train - np.mean(combined_step1_data, axis=0))/np.std(combined_step1_data, axis=0)
    normalized_data_reco = (measured_train - np.mean(combined_step1_data, axis=0))/np.std(combined_step1_data, axis=0)

    combined_step2_data = np.concatenate((MC_train, MC_train))
    normalized_MC_gen = (MC_train - np.mean(combined_step2_data, axis=0))/np.std(combined_step2_data, axis=0)

    MC_dataloader = DataLoader(reco = normalized_MC_reco,
                            gen = normalized_MC_gen,
                            pass_reco = pass_reco_train,
                            normalize = True,
                            rank=hvd.rank(),
                            size=hvd.size(),)
    data_dataloader = DataLoader(reco = normalized_data_reco,
                                pass_reco = measured_pass_reco_train,
                                normalize = True,
                                rank=hvd.rank(),
                                size=hvd.size(),)
    ndim = 1 #The number of features present in your dataset
    reco_model = MLP(ndim)
    model_gen = MLP(ndim)
    omnifold_dnn = MultiFold(
        "DNN_gaussians_02_2025",
        reco_model,
        model_gen,
        data_dataloader, # a dataloader instance containing the measured data
        MC_dataloader,
        niter = itnum,
        verbose=True, # a dataloader instance containing the simulation
        batch_size = 256,
        early_stop=3,
        rank=hvd.rank(),
        size=hvd.size(),
    )
    omnifold_dnn.Unfold()
if __name__ == '__main__':
    main()
