from omnifold import DataLoader, MultiFold, PET
import h5py as h5
import argparse 
import os 
import numpy as np
import horovod.tensorflow.keras as hvd
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a PET model using Pythia and Herwig data.")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--save_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store trained model weights")
    parser.add_argument("--num_data", type=int, default=-1, help="Number of data to train with")
    parser.add_argument("--num_iterations", type=int, default=10, help="Number of iterations to use during training")
    parser.add_argument("--start_N", type=int, default=0, help="Number of iteration to start with")
    args = parser.parse_args()
    return args

def main():
    flags = parse_arguments()
    num_data = flags.num_data
    itnum = flags.num_iterations
    data_dir = flags.data_dir
    synthetic_file_path = data_dir + "train_pythia.h5"
    nature_file_path = data_dir + "train_herwig.h5"

    synthetic  =  h5.File(synthetic_file_path, 'r')
    nature = h5.File(nature_file_path, 'r')
    synthetic_pass_reco = (synthetic['reco_jets'][:num_data,0]>150)
    nature_pass_reco = (nature['reco_jets'][:num_data,0]>150)

    synthetic_gen_parts = synthetic['gen'][:num_data]
    synthetic_reco_parts = synthetic['reco'][:num_data]
    nature_reco_parts = nature['reco'][:num_data]

    synthetic_parts_dataloader = DataLoader(reco = synthetic_reco_parts,
                                            gen = synthetic_gen_parts,
                                            pass_reco = synthetic_pass_reco,
                                            normalize = True,
                                            rank=hvd.rank(),
                                            size=hvd.size(),)
    nature_parts_dataloader = DataLoader(reco = nature_reco_parts,
                                        pass_reco = nature_pass_reco,
                                        normalize = True,
                                        rank=hvd.rank(),
                                        size=hvd.size(),)
    synthetic_PET_model = PET(synthetic_gen_parts.shape[2], num_part=synthetic_gen_parts.shape[1], num_heads = 4, num_transformer = 4, local = True, projection_dim = 128, K = 10)
    nature_PET_model = PET(synthetic_gen_parts.shape[2], num_part=synthetic_gen_parts.shape[1], num_heads = 4, num_transformer = 4, local = True, projection_dim = 128, K = 10)
    omnifold_PET = MultiFold(
        "PET_02_2025",
        model_reco = nature_PET_model,
        model_gen = synthetic_PET_model,
        data = nature_parts_dataloader,
        mc = synthetic_parts_dataloader,
        start = flags.start_N,
        niter = itnum,
        weights_folder = flags.save_dir,
        verbose=True,
        batch_size = 256,
        early_stop=3,
        rank=hvd.rank(),
        size=hvd.size(),
    )
    omnifold_PET.Unfold()

if __name__ == '__main__':
    main()