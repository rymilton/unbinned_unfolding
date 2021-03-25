//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Unit tests for generating the RooUnfoldResponse
//
// Authors: Vincent Croft <vincent.croft@cern.ch>
//
//==============================================================================

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "RooUnfoldResponse.h"
#include "unittests.h"
#include <string>

RooUnfoldResponse BuildRooUnfoldResponse(std::string filename="response.root")
{
    TFile* f = new TFile(filename.c_str(),"OPEN");
    TH2D* h_response = (TH2D*)f->Get("res");
    TH1D* h_gen = (TH1D*)f->Get("gen");
    TH1D* h_sim = (TH1D*)f->Get("sim");
    RooUnfoldResponse response (h_sim,h_gen,h_response);
    f->Close();
    return response;
}