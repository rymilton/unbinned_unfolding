//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Generate base information for unit tests
//
// Authors: Vincent Croft <vincent.croft@cern.ch>
//
//==============================================================================

#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "RooUnfoldResponse.h"

#include "unittests.h"

const Double_t cutdummy= -99999.0;
TRandom randy(666.);

Double_t smear (Double_t xt)
{
  Double_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Double_t x= randy.Rndm();
  if (x>xeff) return cutdummy;
  Double_t xsmear= randy.Gaus(-2.5,0.2);     // bias and smear
  return xt+xsmear;
}

void RooUnfoldGenerate()
{
    TH1D* h_gen= new TH1D ("gen", "Train Truth",    40, -10.0, 10.0);
    TH1D* h_sim= new TH1D ("sim", "Train Measured", 40, -10.0, 10.0);
    TH2D* h_response= new TH2D ("res", "Response Hist", 40, -10.0, 10.0, 40, -10.0, 10.0);
    RooUnfoldResponse response (40, -10.0, 10.0);

    // Train with a Breit-Wigner, mean 0.3 and width 2.5.
    for (Int_t i= 0; i<100000; i++) {
        Double_t xt= randy.BreitWigner (0.3, 2.5);
        Double_t x= smear (xt);
        h_gen->Fill(xt);
        if (x!=cutdummy){
            response.Fill(x, xt);
            h_response->Fill(x, xt);
            h_sim->Fill(x);
        }
        else {
            response.Miss (xt);
        }
    }

    TH1D* h_true= new TH1D ("true", "Test Truth",    40, -10.0, 10.0);
    TH1D* h_meas= new TH1D ("meas", "Test Measured", 40, -10.0, 10.0);
    // Test with a Gaussian, mean 0 and width 2.
    for (Int_t i=0; i<10000; i++) {
        Double_t xt= randy.Gaus (0.0, 2.0), x= smear (xt);
        h_true->Fill(xt);
        if (x!=cutdummy) h_meas->Fill(x);
    }

    // Save objects to file
    TFile* output = new TFile("response.root","RECREATE");
    output->cd();
    h_response->Write();
    h_true->Write();
    h_meas->Write();
    h_gen->Write();
    h_sim->Write();
    output->Write();

}

void RooUnfoldGenerateVariable()
{
    const Int_t NBINS = 40;
    Double_t edges[NBINS + 1] = {-10.        ,  -9.87836329,  -9.74932883,  -9.61244672,
                                 -9.46723965,  -9.31320133,  -9.14979466,  -8.97644986,
                                 -8.79256251,  -8.59749144,  -8.39055647,  -8.17103605,
                                 -7.93816476,  -7.69113062,  -7.42907227,  -7.15107597,
                                 -6.85617238,  -6.54333325,  -6.21146775,  -5.85941874,
                                 -5.48595869,  -5.08978541,  -4.66951753,  -4.22368964,
                                 -3.75074723,  -3.24904125,  -2.71682233,  -2.15223474,
                                 -1.55330985,  -0.91795934,  -0.24396786,   0.47101468,
                                 1.22948128,   2.03407657,   2.88760603,   3.79304574,
                                 4.7535528 ,   5.77247633,   6.8533691 ,   8.        ,10.};
    TH1D* h_gen= new TH1D ("gen", "Train Truth",   NBINS, edges);
    TH1D* h_sim= new TH1D ("sim", "Train Measured", NBINS, edges);
    TH2D* h_response= new TH2D ("res", "Response Hist", NBINS, edges, NBINS, edges);
    RooUnfoldResponse response (h_gen, h_sim);

    // Train with a Breit-Wigner, mean 0.3 and width 2.5.
    for (Int_t i= 0; i<100000; i++) {
        Double_t xt= randy.BreitWigner (0.3, 2.5);
        Double_t x= smear (xt);
        h_gen->Fill(xt);
        if (x!=cutdummy){
            response.Fill(x, xt);
            h_response->Fill(x, xt);
            h_sim->Fill(x);
        }
        else {
            response.Miss (xt);
        }
    }

    TH1D* h_true= new TH1D ("true", "Test Truth",    NBINS, edges);
    TH1D* h_meas= new TH1D ("meas", "Test Measured", NBINS, edges);
    // Test with a Gaussian, mean 0 and width 2.
    for (Int_t i=0; i<10000; i++) {
        Double_t xt= randy.Gaus (0.0, 2.0), x= smear (xt);
        h_true->Fill(xt);
        if (x!=cutdummy) h_meas->Fill(x);
    }

    // Save objects to file
    TFile* output = new TFile("response_var.root","RECREATE");
    output->cd();
    h_response->Write();
    h_true->Write();
    h_meas->Write();
    h_gen->Write();
    h_sim->Write();
    output->Write();
    output->Close();
}