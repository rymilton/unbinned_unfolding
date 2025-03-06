import numpy as np
import ROOT

def smear(xt):
  xeff = 0.3 + (1.0-0.3)/20*(xt+10.0)  #  efficiency                                                                                  
  x = ROOT.gRandom.Rndm()
  if x>xeff: return None
  xsmear = ROOT.gRandom.Gaus(-2.5,0.2)     #  bias and smear 
  return xt + xsmear

def CaloSmear(xt, xmin, xmax, a, b, c):
    # calculate efficiency
    diff = xmax - xmin
    center = xmin + diff/2
    xeff = 0.1 + 0.7 * (1 - ((center - xt)/diff)**2 )

    x = ROOT.gRandom.Rndm()
    if x>xeff: 
        return None
    
    sigma_x = xt * (a/ROOT.TMath.Sqrt(xt)) + b
    xsmear = ROOT.gRandom.Gaus(c, sigma_x)

    return xt + xsmear

def smearing_wrapper(xt, smearing_function):
    if smearing_function == "calo":
        return CaloSmear(xt, 0, 10, .15, .1, -1.25)
    else:
        return smear(xt)

num_MC_data_points  = 2000000
num_measured_data_points = 500000

smearing_function = "smear"

# Generating pseudo data
unbinned_MC_data = []
unbinned_sim_data = []
unbinned_true_data = []
unbinned_measured_data = []
sim_pass_reco = []
measured_pass_reco = []
for i in range(num_MC_data_points):
    # Generating MC data
    if smearing_function == "calo":
      xt = ROOT.gRandom.BreitWigner(5.0, 1.0)
    else:
      xt = ROOT.gRandom.BreitWigner(0.3, 2.5)
    unbinned_MC_data.append(xt)
    x = smearing_wrapper(xt, smearing_function)
    if x != None:
        unbinned_sim_data.append(x)
        sim_pass_reco.append(True)
    else:
        unbinned_sim_data.append(-9999) # Using -9999 as a filler for events that fail reco
        sim_pass_reco.append(False)

for i in range(int(num_measured_data_points)):
    if smearing_function == "calo":
        xt = ROOT.gRandom.Gaus (5.0, 1.0)
    else:
        xt = ROOT.gRandom.Gaus (0.0, 2.0)
    unbinned_true_data.append(xt)
    x = smearing_wrapper(xt, smearing_function)
    if x != None:
        unbinned_measured_data.append(x)
        measured_pass_reco.append(True)
    else:
        unbinned_measured_data.append(-9999)
        measured_pass_reco.append(False)

unbinned_MC_data = np.array(unbinned_MC_data, dtype=np.float32)
unbinned_sim_data = np.array(unbinned_sim_data, dtype=np.float32)
unbinned_true_data = np.array(unbinned_true_data, dtype=np.float32)
unbinned_measured_data = np.array(unbinned_measured_data, dtype=np.float32)
sim_pass_reco = np.array(sim_pass_reco, dtype=bool)
measured_pass_reco = np.array(measured_pass_reco, dtype=bool)

np.savez(
    "gaussian_data_test.npz",
    MC = unbinned_MC_data,
    sim = unbinned_sim_data,
    truth = unbinned_true_data, 
    measured = unbinned_measured_data,
    sim_pass_reco = sim_pass_reco,
    measured_pass_reco = measured_pass_reco
)