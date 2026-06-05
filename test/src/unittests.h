/*===========================================================================*/
/*
 * BEGIN ROOUNFOLD COPYRIGHT
 * RooUnfold — Unfolding library for particle-physics inverse problems
 *
 * Copyright © 2021–2025 CERN and the authors’ respective research institutions
 * Please refer to the CONTRIBUTORS file for details.
 *
 * License: BSD-3-Clause
 * SPDX-License-Identifier: BSD-3-Clause
 *
 * END ROOUNFOLD COPYRIGHT
 */
/*===========================================================================*/

#ifndef __UNITTESTS_H__
#define __UNITTESTS_H__

#include "RooUnfoldResponse.h"
#include "TVector.h"
#include <string>

void RooUnfoldGenerate();
void RooUnfoldGenerateVariable();
RooUnfoldResponse BuildRooUnfoldResponse(std::string);
TVector BuildRooUnfoldBayes(int);
void WriteRooUnfoldBayes(int);
int TestBayes(int);
const char *test_bayes();
#endif
