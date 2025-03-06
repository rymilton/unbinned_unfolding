from omnifold import DataLoader, MultiFold, MLP
import h5py as h5
import argparse 
import numpy as np
import horovod.tensorflow.keras as hvd
hvd.init()
def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a DNN model using Pythia and Herwig data.")
    parser.add_argument("--data_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/datasets/", help="Folder containing input files")
    parser.add_argument("--save_dir", type=str, default="/global/homes/r/rmilton/m3246/rmilton/omnifold_paper_plots/unfolding/weights/", help="Folder to store trained model weights")
    parser.add_argument("--num_data", type=int, default=-1, help="Number of data to train with")
    parser.add_argument("--num_iterations", type=int, default=5, help="Number of iterations to use during training")
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
    obs_multifold = ['Mass', 'Mult', 'Width', 'Tau21', 'zg', 'SDMass']
    synthetic_obs = {'gen_mass':synthetic['gen_subs'][:num_data, 0], 'gen_width':synthetic['gen_subs'][:num_data, 1], 'gen_mult':synthetic['gen_subs'][:num_data, 2],
                    'gen_sdmass':synthetic['gen_subs'][:num_data, 3], 'gen_zg':synthetic['gen_subs'][:num_data, 4], 'gen_tau21':synthetic['gen_subs'][:num_data, 5],
                    'reco_mass':synthetic['reco_subs'][:num_data, 0], 'reco_width':synthetic['reco_subs'][:num_data, 1], 'reco_mult':synthetic['reco_subs'][:num_data, 2],
                    'reco_sdmass':synthetic['reco_subs'][:num_data, 3], 'reco_zg':synthetic['reco_subs'][:num_data, 4], 'reco_tau21':synthetic['reco_subs'][:num_data, 5]}
    nature_obs    = {'gen_mass':nature['gen_subs'][:num_data, 0], 'gen_width':nature['gen_subs'][:num_data, 1], 'gen_mult':nature['gen_subs'][:num_data, 2],
                    'gen_sdmass':nature['gen_subs'][:num_data, 3], 'gen_zg':nature['gen_subs'][:num_data, 4], 'gen_tau21':nature['gen_subs'][:num_data, 5],
                    'reco_mass':nature['reco_subs'][:num_data, 0], 'reco_width':nature['reco_subs'][:num_data, 1], 'reco_mult':nature['reco_subs'][:num_data, 2],
                    'reco_sdmass':nature['reco_subs'][:num_data, 3], 'reco_zg':nature['reco_subs'][:num_data, 4], 'reco_tau21':nature['reco_subs'][:num_data, 5]}
    synthetic_pass_reco = (synthetic['reco_jets'][:num_data,0]>150)
    nature_pass_reco = (nature['reco_jets'][:num_data,0]>150)
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
        ob['truthobs'], ob['dataobs'] = ob['func'](nature_obs, 'gen'), ob['func'](nature_obs, 'reco')
        print('Done with', obkey)

    # Standardizing the data
    MC_gen = np.asarray([obs[obkey]['genobs'] for obkey in obs_multifold]).T
    MC_reco = np.asarray([obs[obkey]['simobs'] for obkey in obs_multifold]).T
    data_reco = np.asarray([obs[obkey]['dataobs'] for obkey in obs_multifold]).T
    combined_step1_data = np.concatenate((MC_reco[synthetic_pass_reco], data_reco[nature_pass_reco]), axis=0)
    normalized_MC_reco = (MC_reco - np.mean(combined_step1_data, axis=0))/np.std(combined_step1_data, axis=0)
    normalized_data_reco = (data_reco - np.mean(combined_step1_data, axis=0))/np.std(combined_step1_data, axis=0)

    combined_step2_data = np.concatenate((MC_gen, MC_gen), axis=0)
    normalized_MC_gen = (MC_gen - np.mean(combined_step2_data, axis=0))/np.std(combined_step2_data, axis=0)

    MC_dataloader = DataLoader(reco = normalized_MC_reco,
                            gen = normalized_MC_gen,
                            pass_reco = synthetic_pass_reco,
                            normalize = True,
                            rank=hvd.rank(),
                            size=hvd.size(),)
    data_dataloader = DataLoader(reco = normalized_data_reco,
                                pass_reco = nature_pass_reco,
                                normalize = True,
                                rank=hvd.rank(),
                                size=hvd.size(),)
    ndim = len(obs_multifold) #The number of features present in your dataset
    reco_model = MLP(ndim)
    gen_model = MLP(ndim)
    omnifold_dnn = MultiFold(
        "DNN_02_2025",
        reco_model,
        gen_model,
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
