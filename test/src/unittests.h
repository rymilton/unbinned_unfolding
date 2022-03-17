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
