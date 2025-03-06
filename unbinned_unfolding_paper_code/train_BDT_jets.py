import h5py as h5
import numpy as np
import argparse 
import ROOT
import gc
import pickle
import gzip

def parse_arguments():
    parser = argparse.ArgumentParser(description="Train a PET model using Pythia and Herwig data.")
    parser.add_argument("--data_dir", type=str, default="/media/miguel/Elements_2024/unfolding_data/Z+jets/", help="Folder containing input files")
    parser.add_argument("--model_save_dir", type=str, default="/home/ryan/unfolding_paper/unfolding/BDT_models/", help="Folder to store trained model weights")
    parser.add_argument("--model_save_name", type=str, default="BDT_5iterations_test", help="Name of saved models")
    parser.add_argument("--unfolding_type", type=str, default="multifold", help="unifold/multifold/binned")
    parser.add_argument("--num_train_data", type=int, default=-1, help="Number of data to train with")
    parser.add_argument("--num_test_data", type=int, default=-1, help="Number of data for model predictions. Only used for binned unfolding")
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

    # First, configure the data
    num_train_data = flags.num_train_data
    itnum = flags.num_iterations

    data_dir = flags.data_dir
    synthetic_file_path = data_dir + "train_pythia.h5"
    nature_file_path = data_dir + "train_herwig.h5"

    synthetic_test_file_path = data_dir + "test_pythia.h5"
    nature_test_file_path = data_dir + "test_herwig.h5"


    synthetic  =  h5.File(synthetic_file_path, 'r')
    nature = h5.File(nature_file_path, 'r')
    
    obs_multifold = ['Mass', 'Mult', 'Width', 'Tau21', 'zg', 'SDMass']

    synthetic_obs = {'gen_mass':synthetic['gen_subs'][:num_train_data, 0], 'gen_width':synthetic['gen_subs'][:num_train_data, 1], 'gen_mult':synthetic['gen_subs'][:num_train_data, 2],
                 'gen_sdmass':synthetic['gen_subs'][:num_train_data, 3], 'gen_zg':synthetic['gen_subs'][:num_train_data, 4], 'gen_tau21':synthetic['gen_subs'][:num_train_data, 5],
                 'reco_mass':synthetic['reco_subs'][:num_train_data, 0], 'reco_width':synthetic['reco_subs'][:num_train_data, 1], 'reco_mult':synthetic['reco_subs'][:num_train_data, 2],
                 'reco_sdmass':synthetic['reco_subs'][:num_train_data, 3], 'reco_zg':synthetic['reco_subs'][:num_train_data, 4], 'reco_tau21':synthetic['reco_subs'][:num_train_data, 5]}
    nature_obs    = {'gen_mass':nature['gen_subs'][:num_train_data, 0], 'gen_width':nature['gen_subs'][:num_train_data, 1], 'gen_mult':nature['gen_subs'][:num_train_data, 2],
                    'gen_sdmass':nature['gen_subs'][:num_train_data, 3], 'gen_zg':nature['gen_subs'][:num_train_data, 4], 'gen_tau21':nature['gen_subs'][:num_train_data, 5],
                    'reco_mass':nature['reco_subs'][:num_train_data, 0], 'reco_width':nature['reco_subs'][:num_train_data, 1], 'reco_mult':nature['reco_subs'][:num_train_data, 2],
                    'reco_sdmass':nature['reco_subs'][:num_train_data, 3], 'reco_zg':nature['reco_subs'][:num_train_data, 4], 'reco_tau21':nature['reco_subs'][:num_train_data, 5]}
    synthetic_pass_reco = (synthetic['reco_jets'][:num_train_data,0]>150)
    nature_pass_reco = (nature['reco_jets'][:num_train_data,0]>150)

    

    # a dictionary to hold information about the observables
    obs = {}

    # the jet mass and histogram style information
    obs.setdefault('Mass', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_mass'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 75), 'ylim': (0, 0.065),
        'xlabel': r'Jet Mass $m$ [GeV]', 'symbol': r'$m$',
        'ylabel': r'Normalized Cross Section [GeV$^{-1}$]',
        'stamp_xy': (0.425, 0.6),
    })

    # the constituent multiplicity and histogram style information
    obs.setdefault('Mult', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_mult'],
        'nbins_det': 80, 'nbins_mc': 80,
        'xlim': (0, 80), 'ylim': (0, 0.065),
        'xlabel': 'Jet Constituent Multiplicity $M$', 'symbol': r'$M$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.44, 0.6),
    })

    # the jet width and histogram style information
    obs.setdefault('Width', {}).update({
        'func': lambda dset, ptype: dset[ptype + '_width'],
        'nbins_det': 50, 'nbins_mc': 50,
        'xlim': (0, 0.6), 'ylim': (0, 10),
        'xlabel': r'Jet Width $w$', 'symbol': r'$w$',
        'ylabel': r'Normalized Cross Section',
        'stamp_xy': (0.425, 0.6),
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
        'stamp_xy': (0.425, 0.6),
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
    mc_entries = np.asarray([(obs[obkey]['genobs']) for obkey in obs_multifold]).T
    sim_entries = np.asarray([(obs[obkey]['simobs']) for obkey in obs_multifold]).T
    measured_entries = np.asarray([(obs[obkey]['dataobs']) for obkey in obs_multifold]).T
    pass_truth_mask = np.full(len(sim_entries), True)
    MC_pass_reco_vector_train = np_to_TVector(synthetic_pass_reco)
    MC_pass_gen_vector_train = np_to_TVector(pass_truth_mask)
    measured_pass_reco_vector = np_to_TVector(nature_pass_reco)
    
    # Use all of the data (including test) when training binned unfolding model
    if flags.unfolding_type == "binned":
        num_test_data = flags.num_test_data
        synthetic_test  =  h5.File(synthetic_test_file_path, 'r')
        nature_test = h5.File(nature_test_file_path, 'r')

        synthetic_test_obs = {'gen_mass':synthetic_test['gen_subs'][:num_test_data, 0], 'gen_width':synthetic_test['gen_subs'][:num_test_data, 1], 'gen_mult':synthetic_test['gen_subs'][:num_test_data, 2],
                    'gen_sdmass':synthetic_test['gen_subs'][:num_test_data, 3], 'gen_zg':synthetic_test['gen_subs'][:num_test_data, 4], 'gen_tau21':synthetic_test['gen_subs'][:num_test_data, 5],
                    'reco_mass':synthetic_test['reco_subs'][:num_test_data, 0], 'reco_width':synthetic_test['reco_subs'][:num_test_data, 1], 'reco_mult':synthetic_test['reco_subs'][:num_test_data, 2],
                    'reco_sdmass':synthetic_test['reco_subs'][:num_test_data, 3], 'reco_zg':synthetic_test['reco_subs'][:num_test_data, 4], 'reco_tau21':synthetic_test['reco_subs'][:num_test_data, 5]}
        nature_test_obs   = {'gen_mass':nature_test['gen_subs'][:num_test_data, 0], 'gen_width':nature_test['gen_subs'][:num_test_data, 1], 'gen_mult':nature_test['gen_subs'][:num_test_data, 2],
                        'gen_sdmass':nature_test['gen_subs'][:num_test_data, 3], 'gen_zg':nature_test['gen_subs'][:num_test_data, 4], 'gen_tau21':nature_test['gen_subs'][:num_test_data, 5],
                        'reco_mass':nature_test['reco_subs'][:num_test_data, 0], 'reco_width':nature_test['reco_subs'][:num_test_data, 1], 'reco_mult':nature_test['reco_subs'][:num_test_data, 2],
                        'reco_sdmass':nature_test['reco_subs'][:num_test_data, 3], 'reco_zg':nature_test['reco_subs'][:num_test_data, 4], 'reco_tau21':nature_test['reco_subs'][:num_test_data, 5]}
        synthetic_test_pass_reco = (synthetic_test['reco_jets'][:num_test_data,0]>150)
        nature_test_pass_reco = (nature_test['reco_jets'][:num_test_data,0]>150)
        obs_test = {}
        # the jet mass and histogram style information
        obs_test.setdefault('Mass', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_mass'],
            'nbins_det': 50, 'nbins_mc': 50,
            'xlim': (0, 75), 'ylim': (0, 0.065),
            'xlabel': r'Jet Mass $m$ [GeV]', 'symbol': r'$m$',
            'ylabel': r'Normalized Cross Section [GeV$^{-1}$]',
            'stamp_xy': (0.425, 0.6),
        })

        # the constituent multiplicity and histogram style information
        obs_test.setdefault('Mult', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_mult'],
            'nbins_det': 80, 'nbins_mc': 80,
            'xlim': (0, 80), 'ylim': (0, 0.065),
            'xlabel': 'Jet Constituent Multiplicity $M$', 'symbol': r'$M$',
            'ylabel': r'Normalized Cross Section',
            'stamp_xy': (0.44, 0.6),
        })

        # the jet width and histogram style information
        obs_test.setdefault('Width', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_width'],
            'nbins_det': 50, 'nbins_mc': 50,
            'xlim': (0, 0.6), 'ylim': (0, 10),
            'xlabel': r'Jet Width $w$', 'symbol': r'$w$',
            'ylabel': r'Normalized Cross Section',
            'stamp_xy': (0.425, 0.6),
        })

        # the N-subjettiness ratio and histogram style information
        obs_test.setdefault('Tau21', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_tau21'],
            'nbins_det': 50, 'nbins_mc': 50,
            'xlim': (0, 1.2), 'ylim': (0, 3),
            'xlabel': r'$N$-subjettiness Ratio $\tau_{21}^{(\beta=1)}$', 'symbol': r'$\tau_{21}^{(\beta=1)}$',
            'ylabel': r'Normalized Cross Section',
            'stamp_xy': (0.41, 0.92),
            'legend_loc': 'upper left', 'legend_ncol': 1,
        })

        # the groomed momentum fraction and histogram style information
        obs_test.setdefault('zg', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_zg'],
            'nbins_det': 50, 'nbins_mc': 50,
            'xlim': (0, 0.5), 'ylim': (0, 9),
            'xlabel': r'Groomed Jet Momentum Fraction $z_g$', 'symbol': r'$z_g$',
            'ylabel': 'Normalized Cross Section',
            'stamp_xy': (0.425, 0.6),
        })

        # the groomed jet mass and histogram style information
        obs_test.setdefault('SDMass', {}).update({
            'func': lambda dset, ptype: dset[ptype + '_sdmass'],
            'nbins_det': 50, 'nbins_mc': 50,
            'xlim': (-14, -2), 'ylim': (0, 0.3),
            'xlabel': r'Soft Drop Jet Mass $\ln\rho$', 'symbol': r'$\ln\rho$',
            'ylabel': r'Normalized Cross Section',
            'stamp_xy': (0.41, 0.92),
            'legend_loc': 'upper left', 'legend_ncol': 1,
        })

        # calculate quantities to be stored in obs test
        for obkey,ob in obs_test.items():
            # calculate observable for GEN, SIM, DATA, and TRUE
            ob['genobs'], ob['simobs'] = ob['func'](synthetic_test_obs, 'gen'), ob['func'](synthetic_test_obs, 'reco')
            ob['truthobs'], ob['dataobs'] = ob['func'](nature_test_obs, 'gen'), ob['func'](nature_test_obs, 'reco')
        
        mc_entries_test = np.asarray([(obs_test[obkey]['genobs']) for obkey in obs_multifold]).T
        sim_entries_test = np.asarray([(obs_test[obkey]['simobs']) for obkey in obs_multifold]).T
        measured_entries_test = np.asarray([(obs_test[obkey]['dataobs']) for obkey in obs_multifold]).T

        mc_entries_all = np.concatenate((mc_entries, mc_entries_test))
        sim_entries_all = np.concatenate((sim_entries, sim_entries_test))
        synthetic_pass_reco_all = np.concatenate((synthetic_pass_reco, synthetic_test_pass_reco))
        measured_entries_all = np.concatenate((measured_entries, measured_entries_test))
        measured_pass_reco_all = np.concatenate((nature_pass_reco, nature_test_pass_reco))

    # Training the BDTs
    if flags.unfolding_type == "multifold":
        df_MCgen_multi_train = ROOT.RDF.FromNumpy({name: mc_entries[:, i] for i, name in enumerate(obs_multifold)})
        df_MCreco_multi_train = ROOT.RDF.FromNumpy({name: sim_entries[:, i] for i, name in enumerate(obs_multifold)})
        df_measured_multi_train = ROOT.RDF.FromNumpy({name: measured_entries[:, i] for i, name in enumerate(obs_multifold)})
        multifold_bdt = ROOT.RooUnfoldOmnifold()
        multifold_bdt.SetModelSaving(True) 
        multifold_bdt.SetSaveDirectory(flags.model_save_dir) 
        multifold_bdt.SetModelName(flags.model_save_name)
        multifold_bdt.SetMCgenDataFrame(df_MCgen_multi_train)
        multifold_bdt.SetMCrecoDataFrame(df_MCreco_multi_train)
        multifold_bdt.SetMeasuredDataFrame(df_measured_multi_train)
        multifold_bdt.SetMCPassReco(MC_pass_reco_vector_train)
        multifold_bdt.SetMCPassTruth(MC_pass_gen_vector_train)
        multifold_bdt.SetMeasuredPassReco(measured_pass_reco_vector)
        multifold_bdt.SetNumIterations(itnum)
        _ = multifold_bdt.UnbinnedOmnifold()
        del multifold_bdt, df_MCgen_multi_train, df_MCreco_multi_train, df_measured_multi_train
        gc.collect()
    elif flags.unfolding_type == "unifold":
        unifold_bdt = ROOT.RooUnfoldOmnifold()

        for i, obs_name in enumerate(obs_multifold):
            print(f"Starting observable {obs_name}")
            df_MCgen_single = ROOT.RDF.FromNumpy({obs_multifold[i]: mc_entries[:, i]})
            df_MCreco_single = ROOT.RDF.FromNumpy({obs_multifold[i]: sim_entries[:, i]})
            df_measured_single = ROOT.RDF.FromNumpy({obs_multifold[i]: measured_entries[:, i]})
            unifold_bdt.SetMCgenDataFrame(df_MCgen_single)
            unifold_bdt.SetMCrecoDataFrame(df_MCreco_single)
            unifold_bdt.SetModelSaving(True) 
            unifold_bdt.SetSaveDirectory(flags.model_save_dir) 
            unifold_bdt.SetModelName(flags.model_save_name+f"_{obs_name}")
            unifold_bdt.SetMeasuredDataFrame(df_measured_single)
            unifold_bdt.SetMCPassReco(MC_pass_reco_vector_train)
            unifold_bdt.SetMCPassTruth(MC_pass_gen_vector_train)
            unifold_bdt.SetMeasuredPassReco(measured_pass_reco_vector)
            unifold_bdt.SetNumIterations(itnum)
            _ = unifold_bdt.UnbinnedOmnifold()
            del df_MCgen_single, df_MCreco_single, df_measured_single
            gc.collect()
    elif flags.unfolding_type == "binned":
        # Binning the data
        sim_data_mass, sim_data_width, sim_data_mult, sim_data_lnrho, sim_data_zgs, sim_data_tau21 = [], [], [], [], [], []
        mc_data_mass, mc_data_width, mc_data_mult, mc_data_lnrho, mc_data_zgs, mc_data_tau21 = [], [], [], [], [], []
        meas_data_mass, meas_data_width, meas_data_mult, meas_data_lnrho, meas_data_zgs, meas_data_tau21 = [], [], [], [], [], []

        for i in range(len(measured_entries_all)):
            meas_data_mass.append(measured_entries_all[i][0])
            meas_data_mult.append(measured_entries_all[i][1])
            meas_data_width.append(measured_entries_all[i][2])
            meas_data_tau21.append(measured_entries_all[i][3])
            meas_data_zgs.append(measured_entries_all[i][4])
            meas_data_lnrho.append(measured_entries_all[i][5])
            
        for i in range(len(mc_entries_all)):
            mc_data_mass.append(mc_entries_all[i][0])
            mc_data_mult.append(mc_entries_all[i][1])
            mc_data_width.append(mc_entries_all[i][2])
            mc_data_tau21.append(mc_entries_all[i][3])
            mc_data_zgs.append(mc_entries_all[i][4])
            mc_data_lnrho.append(mc_entries_all[i][5])
        for i in range(len(sim_entries_all)):
            sim_data_mass.append(sim_entries_all[i][0])
            sim_data_mult.append(sim_entries_all[i][1])
            sim_data_width.append(sim_entries_all[i][2])
            sim_data_tau21.append(sim_entries_all[i][3])
            sim_data_zgs.append(sim_entries_all[i][4])
            sim_data_lnrho.append(sim_entries_all[i][5])

        #Binned Data
        bins_mass, bin_low_mass, bin_high_mass = obs['Mass']['nbins_det'], obs['Mass']['xlim'][0], obs['Mass']['xlim'][1]
        bins_widths, bin_low_widths, bin_high_widths  = obs['Width']['nbins_det'], obs['Width']['xlim'][0], obs['Width']['xlim'][1]
        bins_mult, bin_low_mult, bin_high_mult  = obs['Mult']['nbins_det'], obs['Mult']['xlim'][0], obs['Mult']['xlim'][1]
        bins_lnrho, bin_low_lnrho, bin_high_lnrho  = obs['SDMass']['nbins_det'], obs['SDMass']['xlim'][0], obs['SDMass']['xlim'][1]
        bins_zgs, bin_low_zgs, bin_high_zgs  = obs['zg']['nbins_det'], obs['zg']['xlim'][0], obs['zg']['xlim'][1]
        bins_tau21, bin_low_tau21, bin_high_tau21  = obs['Tau21']['nbins_det'], obs['Tau21']['xlim'][0], obs['Tau21']['xlim'][1]

        meas_data_mass_binned = ROOT.TH1D("meas_mass_hist", "meas_mass_hist", bins_mass, bin_low_mass, bin_high_mass)
        meas_data_width_binned = ROOT.TH1D("meas_width_hist", "meas_width_hist", bins_widths, bin_low_widths, bin_high_widths)
        meas_data_mult_binned = ROOT.TH1D("meas_mult_hist", "meas_mult_hist", bins_mult, bin_low_mult, bin_high_mult)
        meas_data_lnrho_binned = ROOT.TH1D("meas_lnrho_hist", "meas_lnrho_hist", bins_lnrho, bin_low_lnrho, bin_high_lnrho)
        meas_data_zgs_binned = ROOT.TH1D("meas_zgs_hist", "meas_zgs_hist", bins_zgs, bin_low_zgs, bin_high_zgs)
        meas_data_tau21_binned = ROOT.TH1D("meas_tau21_hist", "meas_tau21_hist", bins_tau21, bin_low_tau21, bin_high_tau21)

        response_mass = ROOT.RooUnfoldResponse(bins_mass, bin_low_mass, bin_high_mass, bins_mass, bin_low_mass, bin_high_mass)
        response_width = ROOT.RooUnfoldResponse(bins_widths, bin_low_widths, bin_high_widths, bins_widths, bin_low_widths, bin_high_widths)
        response_mult = ROOT.RooUnfoldResponse(bins_mult, bin_low_mult, bin_high_mult, bins_mult, bin_low_mult, bin_high_mult)
        response_lnrho = ROOT.RooUnfoldResponse(bins_lnrho, bin_low_lnrho, bin_high_lnrho, bins_lnrho, bin_low_lnrho, bin_high_lnrho)
        response_zgs = ROOT.RooUnfoldResponse(bins_zgs, bin_low_zgs, bin_high_zgs, bins_zgs, bin_low_zgs, bin_high_zgs)
        response_tau21 = ROOT.RooUnfoldResponse(bins_tau21, bin_low_tau21, bin_high_tau21, bins_tau21, bin_low_tau21, bin_high_tau21)

        #Unbinned data and fill response matrices
        def build_data(
                response, 
                initial_mc_list, 
                initial_sim_list,
                synthetic_pass_reco,
                initial_measured_list,
                nature_pass_reco,
                measured_hist, 
            ):
            for (sim, MC, pass_reco) in zip(initial_sim_list, initial_mc_list, synthetic_pass_reco):
                if pass_reco:
                    response.Fill(sim, MC)
                else:
                    response.Miss(MC)
            for measured in np.asarray(initial_measured_list)[nature_pass_reco]:
                measured_hist.Fill(measured)
            return response, measured_hist

        response_mass, meas_data_mass_binned = build_data(
            response_mass,
            mc_data_mass,
            sim_data_mass,
            synthetic_pass_reco_all,
            meas_data_mass,
            measured_pass_reco_all,
            meas_data_mass_binned
        )

        response_mult, meas_data_mult_binned = build_data(
            response_mult,
            mc_data_mult,
            sim_data_mult,
            synthetic_pass_reco_all,
            meas_data_mult,
            measured_pass_reco_all,
            meas_data_mult_binned
        )

        response_width, meas_data_width_binned = build_data(
            response_width,
            mc_data_width,
            sim_data_width,
            synthetic_pass_reco_all,
            meas_data_width,
            measured_pass_reco_all,
            meas_data_width_binned
        )

        response_tau21, meas_data_tau21_binned = build_data(
            response_tau21,
            mc_data_tau21,
            sim_data_tau21,
            synthetic_pass_reco_all,
            meas_data_tau21,
            measured_pass_reco_all,
            meas_data_tau21_binned
        )

        response_zgs, meas_data_zgs_binned = build_data(
            response_zgs,
            mc_data_zgs,
            sim_data_zgs,
            synthetic_pass_reco_all,
            meas_data_zgs,
            measured_pass_reco_all,
            meas_data_zgs_binned
        )

        response_lnrho, meas_data_lnrho_binned = build_data(
            response_lnrho,
            mc_data_lnrho,
            sim_data_lnrho,
            synthetic_pass_reco_all,
            meas_data_lnrho,
            measured_pass_reco_all,
            meas_data_lnrho_binned
        )

        response_dict = {
            "Mass":response_mass,
            "Mult":response_mult,
            "Width":response_width,
            "Tau21":response_tau21,
            "zg":response_zgs,
            "SDMass":response_lnrho
        }

        measured_hist_dict = {
            "Mass":meas_data_mass_binned,
            "Mult":meas_data_mult_binned,
            "Width":meas_data_width_binned,
            "Tau21":meas_data_tau21_binned,
            "zg":meas_data_zgs_binned,
            "SDMass":meas_data_lnrho_binned
        }
        for i, (obkey,ob) in enumerate(obs.items()):
            print(f"Doing binned unfolding for {obkey}")
            BDT_binned = ROOT.RooUnfoldOmnifold(response_dict[obkey], measured_hist_dict[obkey], itnum)
            h_unfolded = BDT_binned.Hunfold()
            
            bin_centers = [h_unfolded.GetBinCenter(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
            counts = [h_unfolded.GetBinContent(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
            bin_edges = [h_unfolded.GetBinLowEdge(bin) for bin in range(1, h_unfolded.GetNbinsX() + 2)]
            bin_errors = [h_unfolded.GetBinError(bin) for bin in range(1, h_unfolded.GetNbinsX() + 1)]
            
            del h_unfolded, BDT_binned
            gc.collect()
            
            # Calculate the density-normalized histogram using np.histogram
            density_counts, _ = np.histogram(bin_centers, bins=bin_edges, weights=counts, density=True)
            bin_errors = bin_errors/(np.sum(counts)*(bin_edges[1]-bin_edges[0]))
            # Save counts, errors, and bin centers
            ob['binned_BDT'], ob['binned_BDT_unc'], ob['binned_BDT_bincenters'] = density_counts, bin_errors, bin_centers
            data_to_save = {"binned_BDT" : ob['binned_BDT'], "binned_BDT_unc" : ob['binned_BDT_unc'], "binned_BDT_bincenters" : ob['binned_BDT_bincenters']}
            output_file = f"./saved_binnedBDT_weights_{obkey}.pickle.gz"
            with gzip.open(output_file, 'wb') as f:
                pickle.dump(data_to_save, f)
    
    print("Done")
    # Deleting extra ROOT objects
    del MC_pass_reco_vector_train, MC_pass_gen_vector_train, measured_pass_reco_vector
    gc.collect()
if __name__ == '__main__':
    main()