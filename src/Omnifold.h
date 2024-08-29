//=====================================================================-*-C++-*-
//! \class Omnifold
//! \brief Unbinned & binned versions of Omnifold (ML unfolding) using decision trees
//! \author Ryan Milton <rmilt003@ucr.edu>
//==============================================================================

#ifndef OMNIFOLD_H_
#define OMNIFOLD_H_

#include "RooUnfoldResponse.h"
#include "TH1.h"

class Omnifold {
public:
    void BinnedOmnifold(RooUnfoldResponse response, TH1* measured_hist, Int_t num_iterations);
};
#endif
