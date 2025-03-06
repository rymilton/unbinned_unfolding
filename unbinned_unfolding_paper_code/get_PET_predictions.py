from omnifold import PET
import h5py as h5
import argparse 
import numpy as np
import gzip
import pickle
import horovod.tensorflow.keras as hvd
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Make predictions using trained PET model")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--output_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/final_weights", help="Folder where trained model is saved")
    parser.add_argument("--model_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store trained model weights")
    parser.add_argument("--num_data", type=int, default=-1, help="Number of data to train with")
    parser.add_argument("--num_iterations", type=int, default=5, help="Number of iterations to use during training")
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
    num_data = flags.num_data
    data_dir = flags.data_dir

    synthetic_file_path = data_dir + "test_pythia.h5"
    synthetic  =  h5.File(synthetic_file_path, 'r')
    synthetic_gen_parts = synthetic['gen'][:num_data]
    model = PET(synthetic_gen_parts.shape[2], num_part=synthetic_gen_parts.shape[1], num_heads = 4, num_transformer = 4, local = True, projection_dim = 128, K = 10)
    model.load_weights(flags.model_dir + "OmniFold_PET_02_2025_iter4_step2.weights.h5")
    weights = reweight(synthetic_gen_parts, model)

    data_to_save = {
        "PET_weights": weights
    }

    output_file = f'{flags.output_dir}/saved_PET_weights_5iters.pickle.gz'

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(data_to_save, f)


if __name__ == '__main__':
    main()
