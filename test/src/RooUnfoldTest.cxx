/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2007–2025 CERN and the authors’ respective research institutions
 * Please refer to the CONTRIBUTORS file for details.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * END ROOUNFOLD COPYRIGHT
 */
/*===========================================================================*/

//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id$
//
// Description:
//      Tests RooUnfold package using toy MC generated according to PDFs
//      defined in RooUnfoldTestPdf.cxx or RooUnfoldTestPdfRooFit.cxx.
//      This is the main program. The actual tests are performed using the
//      RooUnfoldTestHarness class.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#include "RooUnfoldTestHarness.h"

RooUnfoldTestHarness *test = 0;

//==============================================================================
// Routine to run with parameters specified as a string
//==============================================================================

void RooUnfoldTest(const char *args = "")
{
   // If run interactively, remove canvas and all histograms that might have been
   // created with a previous invocation.
   delete test;
   test = 0;
   gDirectory->Clear();

   test = new RooUnfoldTestHarness("RooUnfoldTest", args);
   test->Run();
}

#ifndef __CINT__

//==============================================================================
// Main program when run stand-alone
//==============================================================================

int main(int argc, char **argv)
{
   RooUnfoldTestHarness maintest("RooUnfoldTest", argc, argv);
   return maintest.Run();
}

#endif
