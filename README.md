RooUnfold: ROOT Unfolding Framework
===

RooUnfold is a framework for unfolding (AKA "deconvolution" or
"unsmearing").  It can be used from the ROOT prompt, or compiled and
linked against the ROOT libraries.  It currently implements seven
methods:

  - iterative ("Bayesian");
  - singular value decomposition ([SVD](https://arxiv.org/abs/hep-ph/9509307), as proposed by Höcker and Kartvelishvili and implemented in [TSVDUnfold](https://root.cern.ch/doc/master/classTSVDUnfold.html));
  - bin-by-bin (simple correction factors);
  - an interface to the [TUnfold](https://root.cern.ch/doc/master/classTUnfold.html) method developed by Stefan Schmitt
  - simple inversion of the response matrix without regularisation
  - iterative dynamically stabilized (IDS) unfolding
  - usage of gaussian processes (GP) for regularizing a kernel, as developed by Adam Bosson
  - Poisson unfolding, a simple likelihood unfolding

RooUnfold was originally written by Tim Adye, Richard Claridge,
Kerstin Tackmann, and Fergus Wilson (2007). It has since received
contributions from many others.  It is currently maintained by Tim
Adye, Carsten Burgard, Lydia Brenner, and Vincent Croft. If you have
any additional methods, please contact us under [our support mailing
list](mailto:roounfold-support@cern.ch).

See this overview of RooUnfold or the references below for more
information. To cite the RooUnfold package in a publication, you can
refer to this web page and/or the paper:

[Comparison of unfolding methods using RooFitUnfold., L. Brenner et al., International Journal of Modern Physics A, Vol. 35, No. 24, 2050145 (2020)](https://doi.org/10.1142/S0217751X20501456), [ArXiV:1910.14654](https://arxiv.org/abs/1910.14654)

There is extensive documentation available online:
  - the [User Guide](https://gitlab.cern.ch/RooUnfold/documentation/-/blob/master/RooUnfold_user_guide.pdf).
  - the [RooUnfold tutorial](http://statisticalmethods.web.cern.ch/StatisticalMethods/unfolding/RooUnfold_01-Methods/) 
  - auto-generated [Doxygen class documentation](http://roounfold.web.cern.ch/hierarchy.html) for the RooUnfold package
  - RooUnfold package [README](README.md)
  - RooUnfold package [release notes](History.md)

For earlier documentation, check:

Tim Adye, in Proceedings of the PHYSTAT 2011 Workshop on
    Statistical Issues Related to Discovery Claims in Search
    Experiments and Unfolding, CERN, Geneva, Switzerland, 17–20
    January 2011, edited by H.B. Prosper and L. Lyons, CERN–2011–006,
    pp. 313–318. Proceedings [CDS](https://cdsweb.cern.ch/record/1306523), [refactored version](https://roounfold.web.cern.ch/phystat2011_adye.pdf); Slides [original](https://indico.cern.ch/event/107747/contributions/32673/), [refactored](https://roounfold.web.cern.ch/adye_tim.pdf).

Building the Library
---

RooUnfold uses [ROOT](https://root.cern.ch/). The ROOT web site has [instructions](https://root.cern/install/)
for installing ROOT on different systems in various ways.
In particular, ROOT is already installed on [CERN lxplus](https://lxplusdoc.web.cern.ch/). Alternatively, if you have [CVMFS](https://cernvm.cern.ch/fs/), it is available using, eg. for CentOS7:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc10-opt/setup.sh
```

Check out RooUnfold

    git clone ssh://git@gitlab.cern.ch:7999/RooUnfold/RooUnfold.git
    cd RooUnfold

Build with `cmake`, using

    mkdir build
    cd build
    cmake ..
    make -j4
    cd ..
    source build/setup.sh

For backwards compatibility, we also maintain a GNU `Makefile`.

Usage
---

Example usage instructions as well as a more detailed description and
recommendations can be found in the [User Guide](https://gitlab.cern.ch/RooUnfold/documentation/-/blob/master/RooUnfold_user_guide.pdf).

Versions
---

Many versions of RooUnfold are tagged, you can access them after with
`git checkout X.Y.Z` after cloning the repository.  The following list
contains some of the major milestone versions.

  - 3.0.0: Inclusion of RooFit, transition to RooFitUnfold. 
  - 2.0.0: General interface improvements
  - 1.1.1: Legacy port from SVN repository
  