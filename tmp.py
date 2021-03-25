import ROOT
import RooUnfold
f = ROOT.TFile.Open("tmp.root","READ")
unfold = f.Get("unfold")
h_unfolded = unfold.Hunfold()
