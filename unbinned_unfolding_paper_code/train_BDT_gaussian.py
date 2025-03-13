import h5py as h5
import numpy as np
import argparse 
import ROOT
import gc
import pickle
import gzip

sim_bins  = 40
MC_bins   = 40
sim_low   = -10
sim_high  = 10
MC_low    = -10
MC_high   = 10

def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a BDT using Gaussian data.")
    parser.add_argument("--data_dir", type=str, default="/home/ryan/unfolding_paper/unfolding/", help="Folder containing input files")
    parser.add_argument("--model_save_dir", type=str, default="/home/ryan/github_omnifold/BDT_models/", help="Folder to store trained model weights")
    parser.add_argument("--model_save_name", type=str, default="BDT_gaussian", help="Name of saved models")
    parser.add_argument("--unfolding_type", type=str, default="unbinned", help="unbinned/binned")
    parser.add_argument("--num_train_data", type=int, default=1_000_000, help="Number of data to train with")
    parser.add_argument("--num_iterations", type=int, default=5, help="Number of iterations to use during training")
    args = parser.parse_args()
    return args
def np_to_TVector(array):
    vector = ROOT.TVector(len(array))
    for i, entry in enumerate(array):
        vector[i] = entry
    return vector

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
    unbinned_measured_data = gaussian_data["measured"]
    sim_pass_reco = gaussian_data["sim_pass_reco"]
    measured_pass_reco = gaussian_data["measured_pass_reco"]
    
    if len(unbinned_MC_data) != len(unbinned_sim_data):
        raise ValueError("MC and sim data have different number of entries!")
    
    if len(unbinned_MC_data) != len(sim_pass_reco):
        raise ValueError("MC and sim pass reco have different number of entries!")
    
    if len(unbinned_truth_data) != len(unbinned_measured_data):
        raise ValueError("Truth and measured data have different number of entries!")
    
    if len(unbinned_truth_data) != len(measured_pass_reco):
        raise ValueError("Truth and measured pass reco have different number of entries!")
    
    response = ROOT.RooUnfoldResponse(sim_bins, sim_low, sim_high, MC_bins, MC_low, MC_high)
    measured_hist = ROOT.TH1D ("measured_hist", "measured_hist", sim_bins, sim_low, sim_high)

    for (MC, sim, pass_reco) in zip(unbinned_MC_data, unbinned_sim_data, sim_pass_reco):
        if pass_reco:
            response.Fill(sim, MC)
        else:
            response.Miss(MC)
    
    for (measured, pass_reco) in zip(unbinned_measured_data, measured_pass_reco):
        if pass_reco:
            measured_hist.Fill(measured)
    MC_train = unbinned_MC_data[:num_train_data]
    sim_train = unbinned_sim_data[:num_train_data]
    pass_reco_train = sim_pass_reco[:num_train_data]
    # Only use 75% of the measured data during training
    measured_train = unbinned_measured_data[:1000]
    measured_pass_reco_train = measured_pass_reco[:1000] 
    
    if flags.unfolding_type == "unbinned":
        df_MCgen = ROOT.RDF.FromNumpy({"MCgen": MC_train})
        df_MCreco = ROOT.RDF.FromNumpy({"MCreco": sim_train})
        df_measured = ROOT.RDF.FromNumpy({"measured": measured_train})

        MC_pass_reco_vector = np_to_TVector(pass_reco_train)
        measured_pass_reco_vector = np_to_TVector(measured_pass_reco_train)

        unbinned_BDT = ROOT.RooUnfoldOmnifold()
        unbinned_BDT.SetModelSaving(True) 
        unbinned_BDT.SetSaveDirectory(flags.model_save_dir) 
        unbinned_BDT.SetModelName(flags.model_save_name)
        unbinned_BDT.SetMCgenDataFrame(df_MCgen)
        unbinned_BDT.SetMCrecoDataFrame(df_MCreco)
        unbinned_BDT.SetMeasuredDataFrame(df_measured)
        unbinned_BDT.SetMCPassReco(MC_pass_reco_vector)
        unbinned_BDT.SetMeasuredPassReco(measured_pass_reco_vector)
        unbinned_BDT.SetNumIterations(itnum)
        _ = unbinned_BDT.UnbinnedOmnifold()
    elif flags.unfolding_type == "binned":
        BDT_binned = ROOT.RooUnfoldOmnifold(response, measured_hist, itnum)
        h_unfolded = BDT_binned.Hunfold()
        bin_centers = [h_unfolded.GetBinCenter(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
        counts = [h_unfolded.GetBinContent(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
        bin_edges = [h_unfolded.GetBinLowEdge(bin) for bin in range(1, h_unfolded.GetNbinsX() + 2)]
        bin_errors = [h_unfolded.GetBinError(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
        
        # Calculate the density-normalized histogram using np.histogram
        density_counts, _ = np.histogram(bin_centers, bins=bin_edges, weights=counts, density=True)
        bin_errors = bin_errors/(np.sum(counts)*(bin_edges[1]-bin_edges[0]))
        # Save counts, errors, and bin centers
        data_to_save = {"binned_BDT" : density_counts, "binned_BDT_unc" : bin_errors, "binned_BDT_bincenters" : bin_centers}
        output_file = f"./saved_binnedBDT_weights_gaussian.pickle.gz"
        with gzip.open(output_file, 'wb') as f:
            pickle.dump(data_to_save, f)
    else:
        raise ValueError("The only unfolding types available are binned and unbinned")
if __name__ == '__main__':
    main()