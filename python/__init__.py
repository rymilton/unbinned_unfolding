# BEGIN ROOUNFOLD COPYRIGHT
# RooUnfold — Unfolding library for particle-physics inverse problems
#
# Copyright © 2021–2025 CERN and the authors’ respective research institutions
# Please refer to the CONTRIBUTORS file for details.
#
# License: BSD-3-Clause
# SPDX-License-Identifier: BSD-3-Clause
#
# END ROOUNFOLD COPYRIGHT

import ROOT

ROOT.gSystem.Load("libRooUnfold")
RooUnfoldResponse = ROOT.RooUnfoldResponse
RooUnfoldBayes = ROOT.RooUnfoldBayes
RooUnfoldBinByBin = ROOT.RooUnfoldBayes
RooUnfoldInvert = ROOT.RooUnfoldInvert
RooUnfoldSvd = ROOT.RooUnfoldSvd
RooUnfoldTUnfold = ROOT.RooUnfoldTUnfold
RooUnfoldSpec = ROOT.RooUnfoldSpec
RooUnfoldFunc = ROOT.RooUnfoldFunc
RooFitUnfoldResponse = ROOT.RooFitUnfoldResponse
RooFitUnfoldBayes = ROOT.RooFitUnfoldBayes
RooFitUnfoldSvd = ROOT.RooFitUnfoldSvd
RooFitUnfoldBinByBin = ROOT.RooFitUnfoldBinByBin
RooFitHist = ROOT.RooUnfolding.RooFitHist
RooUnfolding = ROOT.RooUnfolding
