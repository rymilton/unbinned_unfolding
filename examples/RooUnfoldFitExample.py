#!/usr/bin/env python
# ==============================================================================
#  File and Version Information:
#       $Id$
#
#  Description:
#       This example shows a simple unfolding scenario. Truth data is generated
#       according to a bimodal p.d.f. and then simulates detector response on
#       these events in a naive way by applying a simple smearing and 
#       efficiency function. The resulting truth and reco level events are
#       then used to build a response matrix that can be used for unfolding.
#       The to-be unfolded data is simulated by Poisson smearing a disjoint
#       reco level dataset. This data is then unfolded with a unfolding method
#       of choice and, in case of regularization, a preset set of regularization
#       strengths.
#
#  Author: Pim Verschuuren <pim.verschuuren@rhul.ac.uk>
#
# ==============================================================================



# ====== Plotting ====== #


def plot_input(truth, reco, data, response, purity, eff):

  import ROOT

  c = ROOT.TCanvas("input","input", 900, 500)
  c.Divide(2,2)
  
  c.cd(1)

  ROOT.gStyle.SetOptStat(0)

  truth.SetLineColor(46)
  truth.GetYaxis().SetRangeUser(0,1.2*truth.GetMaximum())
  truth.GetYaxis().SetTitle("Events")
  truth.GetXaxis().SetTitle("#Delta#eta")
  truth.SetTitle("Input distributions")
  truth.Draw("HIST")

  reco.SetLineColor(38)
  reco.SetTitle("")
  reco.Draw("HISTSAME")

  data.SetMarkerColor(ROOT.kBlack)  
  data.SetLineColor(ROOT.kBlack)  
  data.SetMarkerStyle(20)
  data.SetMarkerSize(0.5)
  data.SetTitle("")
  data.Draw("SAMEPE")

  leg = ROOT.TLegend(0.65, 0.65, 0.85, 0.85)
  leg.SetLineWidth(0)
  leg.SetFillStyle(0)
  leg.AddEntry( truth, "Truth MC", "l" )
  leg.AddEntry( reco, "Reco MC", "l" )
  leg.AddEntry( data, "Measured Data", "pe" )
  leg.Draw()
  
  c.cd(2)

  ROOT.gStyle.SetOptStat(0)

  response.GetXaxis().SetTitle("#Delta#eta_{reco}")
  response.GetYaxis().SetTitle("#Delta#eta_{truth}")
  response.SetTitle("Response matrix")
  response.Draw("COLZ")
  c.cd(3)
  
  ROOT.gStyle.SetOptStat(0)
  
  purity.SetLineColor(30)
  purity.GetYaxis().SetRangeUser(0,1.2*purity.GetMaximum())
  purity.GetYaxis().SetTitle("Purity")
  purity.GetXaxis().SetTitle("#Delta#eta_{reco}")
  purity.SetTitle("Bin purity")
  purity.Draw()

  c.cd(4)

  ROOT.gStyle.SetOptStat(0)
  
  eff.SetLineColor(30)
  eff.GetYaxis().SetRangeUser(0,1.2*eff.GetMaximum())
  eff.GetYaxis().SetTitle("Efficiency")
  eff.GetXaxis().SetTitle("#Delta#eta_{truth}")
  eff.SetTitle("Reconstruction efficiency")
  eff.Draw()
  
  c.SaveAs("input.png")



def plot_result(theory, unfold_func, errors):

  import ROOT

  print('Plotting unfolding results...')

  truth_obs = unfold_func.unfolding().response().Htruth().obs(0)
  frame = truth_obs.frame()

  frame.GetXaxis().SetTitle("#Delta#eta")
  frame.GetYaxis().SetTitle("events")

  canvas_truth = ROOT.TCanvas("result","result")

  ROOT.gPad.SetLeftMargin(0.15)
 

  if errors:
    unfold_func.plotOn(frame,ROOT.RooFit.DrawOption("P"),ROOT.RooFit.Name("plot_unfolded"),ROOT.RooFit.MarkerStyle(20),ROOT.RooFit.VisualizeError(ROOT.RooFitResult.prefitResult(unfold_func.makeParameterList())))
  else:
    unfold_func.plotOn(frame,ROOT.RooFit.DrawOption("P"),ROOT.RooFit.Name("plot_unfolded"),ROOT.RooFit.MarkerStyle(20))

  frame.Draw()
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)  
  theory.SetLineColor(ROOT.kRed)
  theory.GetXaxis().SetTitle("#Delta#eta")
  theory.Draw("HISTSAME")
  leg_truth = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
  leg_truth.SetLineWidth(0)
  leg_truth.SetFillStyle(0)
  leg_truth.AddEntry( theory, "Truth MC", "l" )
  leg_truth.AddEntry( frame.findObject("plot_unfolded"), "Unfolded Data", "pe" )
  leg_truth.Draw()
  canvas_truth.SaveAs("result.png")




# ====== Efficiencies ====== #

# A quadratic efficiency.
def quad_eff(xt, xmin, xmax):

  diff = xmax - xmin

  center = xmin + diff/2

  xeff = 0.1 + 0.7 * (1 - ((center - xt)/diff)**2 )

  return xeff

# A linear efficiency.
def lin_eff(xt, xmin, xmax):

  xeff = 0.3 + (1.0-0.3)/(xmax - xmin)*(xt - xmin)

  return xeff


# ====== Smearing ====== #

def smear(xt, xmin, xmax, bias, sigma, eff):

  from ROOT import gRandom

  if eff == 'lin':
    xeff = lin_eff(xt, xmin, xmax)
  elif eff == 'quad':
    xeff = quad_eff(xt, xmin, xmax)
  else:
    print("Unknown passed efficiency type. Pass either lin or quad.")
    exit(0)

  x= gRandom.Rndm();

  if x>xeff: return None;

  xsmear= gRandom.Gaus(bias,sigma);     #  bias and smear
  
  return xt+xsmear



def prepare_bimodal(dict_data):

  print("Generating the data...")

  import ROOT

  # Get the parameters.
  xmin = dict_data['xmin']
  xmax = dict_data['xmax']
  recoBins = dict_data['rbins']
  truthBins = dict_data['tbins']
  bias = dict_data['bias']
  sigma = dict_data['sigma']
  frac = dict_data['frac']
  nevents = dict_data['nevents']
  bkg_events = dict_data['nbkg']
  mcLumiFactor = dict_data['mcLumiFactor']
  eff_type = dict_data['eff']

  # Define the response matrix as a 2D histogram.
  response= ROOT.RooUnfoldResponse (recoBins, xmin, xmax, truthBins, xmin, xmax);  

  # Define the truth histogram that will be used for testing.
  truth = ROOT.TH1D ("Truth", "Truth", truthBins, xmin, xmax);

  # Define the reco histogram to compare to the data.
  reco = ROOT.TH1D ("Reco", "Reco", recoBins, xmin, xmax);
  
  # Define the data histogram.
  data = ROOT.TH1D ("Data", "Data", recoBins, xmin, xmax);

  for i in range(1,truth.GetNbinsX() + 1):
    truth.SetBinContent(i, bkg_events)
  for i in range(reco.GetNbinsX() + 2):
    reco.SetBinContent(i, bkg_events)
    data.SetBinContent(i, bkg_events)

  # Generate events according too a simple double Gaussian
  # distribution.
  for i in range(nevents*mcLumiFactor):

    if i%2 == 0:
      peak_loc = 0.3
    else:
      peak_loc = 0.7

    x_truth = ROOT.gRandom.Gaus(xmin + (xmax-xmin)*peak_loc, (xmax-xmin)*0.1)

    # Use a specific smearing function that also models
    # detector inefficiencies. 
    x_reco = smear(x_truth, xmin, xmax, bias, sigma, eff_type)

    # The truth event is reconstructed.
    if x_reco!=None:
      response.Fill (x_reco, x_truth, 1./mcLumiFactor);

    # The truth event did not pass inefficiencies.
    else:
      response.Miss (x_truth, 1./mcLumiFactor);

    # Fill the testing histograms with seperate events.
    x_truth = ROOT.gRandom.Gaus(xmin + (xmax-xmin)*peak_loc, (xmax-xmin)*0.1)
    x_reco = smear(x_truth, xmin, xmax, bias, sigma, eff_type)
 
    if x_reco!=None: reco.Fill(x_reco, 1./mcLumiFactor)      

    truth.Fill(x_truth, 1./mcLumiFactor)

    # Introduce a shape difference for the data.
    if i > int(frac*nevents*mcLumiFactor):
      peak_loc = 0.3
    else:
      peak_loc = 0.7

    x_truth = ROOT.gRandom.Gaus(xmin + (xmax-xmin)*peak_loc, (xmax-xmin)*0.1)
    x_data = smear(x_truth, xmin, xmax, bias, sigma, eff_type)
    
    if x_data!=None: data.Fill(x_data, 1./mcLumiFactor)

     
  # Create Poisson distributed observed data events.
  for i in range(data.GetNbinsX() + 2):
    x_data = ROOT.gRandom.Poisson(data.GetBinContent(i))
    data.SetBinContent(i, x_data)

  eff = ROOT.RooUnfolding.convertTH1(response.Vefficiency(),response.Htruth())
  pur = ROOT.RooUnfolding.convertTH1(response.Vpurity(),response.Htruth())  
    
  # Save the histograms in a dictionary.
  histograms = {"truth_train":response.Htruth(),
                "reco_train":response.Hmeasured(),
                "truth_test":truth,
                "reco_test":reco,
                "data":data,
                "response":response.Hresponse(),
                "purity":pur,
                "eff":eff
              }
  
  
  for i in range(0,histograms["response"].GetNbinsX()+2):
    for j in range(0,histograms["response"].GetNbinsY()+2):
      histograms["response"].SetBinError(i,j,0)

  return histograms



# Unfolding methods.
def algorithm(method):

  import ROOT

  alg = None
  if method == "bayes":
    alg= ROOT.RooUnfolding.kBayes;
  elif method == "bbb":
    alg= ROOT.RooUnfolding.kBinByBin;
  elif method == "inv":
    alg= ROOT.RooUnfolding.kInvert;
  elif method == "svd":
    alg= ROOT.RooUnfolding.kSVD;
  elif method == "root":
    alg= ROOT.RooUnfolding.kTUnfold;
  elif method == "ids":
    alg= ROOT.RooUnfolding.kIDS;
  elif method == "gp":
    alg= ROOT.RooUnfolding.kGP;
  elif method == "ps":
    alg= ROOT.RooUnfolding.kPoisson;
  else:
    print("The passed unfolding method does not match any of the supported methods. Please pass one of the following methods:")
    print("bayes")
    print("bbb")
    print("inv")
    print("svd")
    print("root")
    print("ids")
    print("gp")
    print("ps")
    exit(0)

  return alg






# ======= Settings ======= #

dict_data = {}       

# Kinematic ranges. (delta eta)
dict_data['xmin'] = -4         
dict_data['xmax'] = 4          
                      
# Smearing parameters.
dict_data['bias'] = 0.2
dict_data['sigma'] = 0.5

# Efficiency type
dict_data['eff'] = 'lin'
                         
# Bimodal assymmetry.
dict_data['frac'] = 0.5
                                                  
# Overflow flag.         
overflow = True
                                                        
# Number of events.            
dict_data['nevents'] = 50000    
                                       
# Constant background bin count        
dict_data['nbkg'] = 0

# Truth/Reco bin count.  
dict_data['tbins'] = 20  
dict_data['rbins'] = 30  
                               
# How much more MC do we have than data
dict_data['mcLumiFactor'] = 1         


# ======================== #




def main(args):

  import ROOT

  # # Prepare the histograms and response matrix.
  histograms = prepare_bimodal(dict_data)

  # Plot the input distributions.
  #plot_input(histograms["truth_train"], histograms["reco_train"], histograms["data"], histograms["response"], histograms["purity"], histograms["eff"])  

  #return
  # Set the unfolding method and corresponding
  # regularization parameter values. See helpers.py
  # for the actual values.
  alg = algorithm(args.method)

  # Get the truth distribution for training.
  train_truth = histograms["truth_train"]
  
  # Get the reco distribution for training.
  train_reco = histograms["reco_train"]
  
  # Get the to-be unfolded data distribution.
  data = histograms["data"]

  # Get the response matrix.
  response = histograms["response"]

  # Set the regularization parameter.
  regparm = args.regparm
  
  print('Unfolding with the '+str(args.method)+r'-method and regularization strength of '+str(regparm))
  
  # Create a spectator object.  
  spec = ROOT.RooUnfoldSpec("unfold","unfold",train_truth,"obs_truth",train_reco,"obs_reco",response,0,data,overflow,0.0005)
  if not regparm == None:
    func = spec.makeFunc(alg, regparm)
  else:
    func = spec.makeFunc(alg)

  # Create a RooFitHist object as input for the unfolding
  # results printing.
  test_truth_RFH = spec.makeHistogram(histograms["truth_test"])

  unfold = func.unfolding()

  #unfold.SetTruth(test_truth_RFH)
  
  # Print the unfolding results and compare to a truth histogram.
  unfold.PrintTable(ROOT.cout, test_truth_RFH, ROOT.RooUnfolding.kRooFit)

  # Plot the results.
  #plot_result(histograms["truth_test"], func, True)

  


if __name__=="__main__":
  from argparse import ArgumentParser
  parser = ArgumentParser(description="RooUnfold testing script")
  parser.add_argument("method",default="inv",type=str)
  parser.add_argument("--regparm",type=float)
  main(parser.parse_args())
  
