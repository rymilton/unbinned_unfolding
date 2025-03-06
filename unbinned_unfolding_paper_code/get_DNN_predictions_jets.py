from omnifold import MLP
import h5py as h5
import argparse 
import numpy as np
import horovod.tensorflow.keras as hvd
import gzip
import pickle
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a PET model using Pythia and Herwig data.")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--model_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/final_weights", help="Folder where the model is saved")
    parser.add_argument("--output_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store predictions ")
    parser.add_argument("--num_data", type=int, default=-1, help="Number of data to train with")
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
    obs_multifold = ['Mass', 'Mult', 'Width', 'Tau21', 'zg', 'SDMass']
    synthetic_obs = {'gen_mass':synthetic['gen_subs'][:num_data, 0], 'gen_width':synthetic['gen_subs'][:num_data, 1], 'gen_mult':synthetic['gen_subs'][:num_data, 2],
                    'gen_sdmass':synthetic['gen_subs'][:num_data, 3], 'gen_zg':synthetic['gen_subs'][:num_data, 4], 'gen_tau21':synthetic['gen_subs'][:num_data, 5],
                    'reco_mass':synthetic['reco_subs'][:num_data, 0], 'reco_width':synthetic['reco_subs'][:num_data, 1], 'reco_mult':synthetic['reco_subs'][:num_data, 2],
                    'reco_sdmass':synthetic['reco_subs'][:num_data, 3], 'reco_zg':synthetic['reco_subs'][:num_data, 4], 'reco_tau21':synthetic['reco_subs'][:num_data, 5]}
    # a dictionary to hold information about the observables
    obs = {}

    # the jet mass and histogram style information
    obs.setdefault('Mass', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_mass'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 75), 'ylim': (0, 0.065),
        'xlabel': r'Jet Mass $m$ [GeV]', 'symbol': r'$m$',
        'ylabel': r'Normalized Cross Section [GeV$^{-1}$]',
        'stamp_xy': (0.425, 0.65),
    })

    # the constituent multiplicity and histogram style information
    obs.setdefault('Mult', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_mult'],
        'nbins_det': 80, 'nbins_mc': 80,
        'xlim': (0, 80), 'ylim': (0, 0.065),
        'xlabel': 'Jet Constituent Multiplicity $M$', 'symbol': r'$M$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.42, 0.65),
    })

    # the jet width and histogram style information
    obs.setdefault('Width', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_width'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 0.6), 'ylim': (0, 10),
        'xlabel': r'Jet Width $w$', 'symbol': r'$w$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.425, 0.65),
    })

    # the N-subjettiness ratio and histogram style information
    obs.setdefault('Tau21', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_tau21'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 1.2), 'ylim': (0, 3),
        'xlabel': r'$N$-subjettiness Ratio $\tau_{21}^{(\beta=1)}$', 'symbol': r'$\tau_{21}^{(\beta=1)}$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.41, 0.92),
        'legend_loc': 'upper left', 'legend_ncol': 1,
    })

    # the groomed momentum fraction and histogram style information
    obs.setdefault('zg', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_zg'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 0.5), 'ylim': (0, 9),
        'xlabel': r'Groomed Jet Momentum Fraction $z_g$', 'symbol': r'$z_g$',
        'ylabel': 'Normalized Cross Section',
        'stamp_xy': (0.425, 0.65),
    })

    # the groomed jet mass and histogram style information
    obs.setdefault('SDMass', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_sdmass'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (-14, -2), 'ylim': (0, 0.3),
        'xlabel': r'Soft Drop Jet Mass $\ln\rho$', 'symbol': r'$\ln\rho$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.41, 0.92),
        'legend_loc': 'upper left', 'legend_ncol': 1,
    })

    # calculate quantities to be stored in obs
    for obkey,ob in obs.items():
        # calculate observable for GEN, SIM, DATA, and TRUE
        ob['genobs'], ob['simobs'] = ob['func'](synthetic_obs, 'gen'), ob['func'](synthetic_obs, 'reco')
        print('Done with', obkey)

    MC_gen = np.asarray([obs[obkey]['genobs'] for obkey in obs_multifold]).T

    normalized_MC_gen = (MC_gen - np.mean(MC_gen, axis=0))/np.std(MC_gen, axis=0)

    ndim = len(obs_multifold) #The number of features present in your dataset
    gen_model = MLP(ndim)
    gen_model.load_weights(flags.model_dir + "OmniFold_DNN_02_2025_iter4_step2.weights.h5")

    weights = reweight(normalized_MC_gen, gen_model)
    # Combine the variables into a dictionary to save them in a single pickle file
    data_to_save = {
        "DNN_weights": weights
    }

    # Saving the combined data as a compressed pickle.gz file
    output_file = f'{flags.output_dir}/saved_DNN_weights_5iters.pickle.gz'

    with gzip.open(output_file, 'wb') as f:
        pickle.dump(data_to_save, f)

if __name__ == '__main__':
    main()
