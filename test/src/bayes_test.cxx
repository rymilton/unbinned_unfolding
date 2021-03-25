//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unit tests for RooUnfoldBayes
//
// Authors: Vincent Croft <vincent.croft@cern.ch>
//
//==============================================================================

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TFile.h"
#include "TVector.h"
#include "unittests.h"

#include <iostream>
#include <fstream>

TVector BuildRooUnfoldBayes(int mode=0)
{
    char* filename;
    if (mode==0)
        filename = "response.root";
    else if (mode==1)
        filename = "response_var.root";
    TFile* f = new TFile(filename, "OPEN");
    TH1D* h_meas = (TH1D*)f->Get("meas");
    RooUnfoldResponse response = BuildRooUnfoldResponse(filename);
    RooUnfoldBayes unfold(&response, h_meas, 4);
    TH1D* h_unfolded = (TH1D*)unfold.Hunfold();
    TVector u(h_unfolded->GetNbinsX());
    for (int i=0; i<h_unfolded->GetNbinsX(); i++){
        u[i] = h_unfolded->GetBinContent(i+1);
    }
    return u;
}

void WriteRooUnfoldBayes(int mode=0)
{
    TVector u = BuildRooUnfoldBayes(mode);
    std::ofstream ref;
    char* filename;
    if (mode==0)
        filename = "../ref/bayes.ref";
    else if (mode==1)
        filename = "../ref/bayes_var.ref";
    ref.open(filename);
    for (int i=0;i<u.GetNrows();i++){
        ref << u[i]<<std::endl;
    }
    ref.close();
}
