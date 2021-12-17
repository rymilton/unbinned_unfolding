#include "src/unittests.h"
#include <iostream>
using std::cout;
using std::endl;

int main(){
    RooUnfoldGenerate();
    RooUnfoldGenerateVariable();
    WriteRooUnfoldBayes(0);
    WriteRooUnfoldBayes(1);
}
