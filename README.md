OmniFold Boosted Decision Tree implementation in RooUnfold
===
This repository contains the RooUnfoldOmnifold class that is an implementation of Omnifold using boosted decision trees (BDTs) within a fork of the RooUnfold framework. RooUnfold itself and associated details can be found in its [GitLab repository](https://gitlab.cern.ch/RooUnfold/RooUnfold). Thank you to Lydia Brenner and Vincent Croft of the RooUnfold team for providing guidance with making this implementation.

**If using the software in this repository, please cite R. Milton, et al., "Tools for Unbinned Unfolding", JINST 20 P05034 (2025).**

### Prerequisites
- [ROOT](https://root.cern/install/) (>= 6.10). This code has been tested with ROOT 6.32-6.36.
- [scikit-learn](https://scikit-learn.org/stable/install.html) (`pip install scikit-learn`)
- NumPy (`pip install numpy`)
  
Building the library
---

RooUnfold uses [ROOT](https://root.cern.ch/). The ROOT web site has [instructions](https://root.cern/install/)
for installing ROOT on different systems in various ways.
In particular, ROOT is already installed on [CERN lxplus](https://lxplusdoc.web.cern.ch/). Alternatively, if you have [CVMFS](https://cernvm.cern.ch/fs/), it is available using, eg. for CentOS7:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
```

Check out this repository

    git clone https://github.com/rymilton/unbinned_unfolding.git
    cd unbinned_unfolding

Build with `cmake`, using

    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..
    source build/setup.sh
    
Please make sure that you source the setup.sh script before using this class. The RooUnfoldOmnifold class looks for a RooUnfold installation in the PATH variables and will use choose the installation that is most recently sourced. 

Usage
---
Similar to the base RooUnfold, the RooUnfoldOmnifold class can be used in both ROOT in C++ and PyROOT. 
An example notebook for using the RooUnfoldOmnifold class in PyROOT can be found at `examples/OmniFold_example.ipynb`.

Scripts used for "Tools for Unbinned Unfolding"
---
This class is featured in the paper "Tools for Unbinned Unfolding". The Python scripts and Jupyter notebook used to generate the plots for that paper are in the `unbinned_unfolding_paper_code` directory. This directory contains the following:

Gaussian data scripts/notebooks:
- generate_gaussians.py -- Generates and saves the data used in the Gaussian results
- train_DNN_gaussians.py & get_DNN_predictions_gaussian.py -- Scripts used to train the DNN on the Gaussian data and make predictions on test Gaussian data with the trained model
- train_BDT_gaussian.py -- Used to train the binned and unbinned Gaussian BDT models
- final_gaussians.ipynb -- Generates the final Gaussian plots in the paper

Jet data scripts/notebooks:
- train_DNN_jets.py & get_DNN_predictions_jets.py -- Scripts used to train the DNN on the jet data and make predictions on test jet data with the trained model
- train_PET.py & get_PET_predictions.py -- Scripts used to train the PET on the jet data and make predictions on test jet data with the trained model
- train_BDT_jets.py -- Used to train the binned and unbinned jet BDT models
- final_jet_plots.ipynb -- Generates the final jet plots in the paper

The jet data can be obtained using the [preprocess_omnifold.py script from the OmniLearn respository](https://github.com/ViniciusMikuni/OmniLearn/blob/main/preprocessing/preprocess_omnifold.py). 

`modplot.py` contains plotting functions used in the notebooks.

The PET and DNN models can be installed from the [PyPI omnifold module](https://pypi.org/project/omnifold/) using `pip install omnifold`.

Common issues
---
1. `AttributeError: Failed to get attribute RooUnfoldOmnifold from ROOT`
   
   This likely means that the RooUnfold being used isn't the one installed from this repository. This typically happens if a different RooUnfold version's setup.sh script is used before using the RooUnfoldOmnifold class. To fix this, please source the setup.sh associated with this repository again.
   
Bugs/Issues
---

To report any bugs or issues found, please contact rmilt003@ucr.edu.
  
