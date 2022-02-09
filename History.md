Release notes
=========

Version [2.0.0](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/2.0.0)
---
16 September 2019

First new release in 8 years collects all changes since [version 1.1.1](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/1.1.1)
in the [SVN development version](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/svn-last)
and migration to [GitLab](https://gitlab.cern.ch/RooUnfold/RooUnfold).

##### New Features

  - Migrated repository to [GitLab](https://gitlab.cern.ch/RooUnfold/RooUnfold) from [SVN](https://svnsrv.desy.de/websvn/wsvn/General.unfolding/RooUnfold/trunk/), added continuous integration and unit tests.
  - Added compilation chain with CMake, both with standalone ROOT and in ATLAS environment
  - Added new [documentation webpage](http://roounfold.web.cern.ch/), auto-generated with doxygen and deployed with the CI
  - Include effect of response matrix errors in unfolded errors if
    `RooUnfold::IncludeSystematics()` is set.
    If `RooUnfold::IncludeSystematics(2)` is used, compute just the systematic
    error, leaving out the measurement errors.
    Use `RooUnfoldTest dosys=1` to test.
  - Add IDS (iterative dynamically stabilized) unfolding implemented
    in the `RooUnfoldIds` class by Chris Meyer (Indiana) using
    the algorithm from Bogdan Malaescu (LPNHE).
  - Add `RooUnfoldResponse::MakeFoldingFunction` which wraps a `TF1` to apply the
    response matrix to a user parametric function. This can be fitted to the
    measured distribution as an alternative to unfolding.\
    `RooUnfoldTest fit=1` tests this (test function definitions only implemented
    for the 1D case).

##### Improvements

  - Move documentation to GitLab [README.md](README.md) and [History.md](History.md) from [RAL web page](http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html). Old History file renamed to [History.svn](History.svn) and frozen.
  - Use a new local version of TSVDUnfold from Kerstin Tackmann which uses
    propagation of errors (for measurement and response) rather than toys.
    This version includes previous memory leak fixes.
    The local version is always used unless explicitly disabled with
    `make HAVE_TSVDUNFOLD=0`.
  - New `RooUnfold::Wreco()` method returns weight matrix of unfolded result.
    This is normally just the inverse of the covariance matrix, but can be
    calculated directly for RooUnfoldSvd using `TSVDUnfold::GetXinv`.
    The weights are used to calulcate the &chi;<sup>2</sup>.
  - Remove RooUnfoldSvd's `_ntoyssvd` parameter. Constructor argument and set/get
    methods are now depreciated (they just set/get RooUnfold's `_NToys`).
  - RooUnfoldSvd deals with fakes by subtracting from the measured input, rather
    than adding an extra effectless-causes bin, which didn't work very well if
    it was much different from its neighbouring bin.
    Giovanni Marchiori reported this problem.
  - Add RooUnfoldTUnfold methods for retrieving the L Curve and the tau values
    if an L curve scan is done.
    This improvement was kindly provided by Pedro Cipriano.
  - Add RooUnfoldTUnfold constructor that fixes regularisation parameter
    (equivalent to using `SetRegParm`).
  - Optimise calculation of the error propagation matrix in RooUnfoldBayes using
    TMatrixD multiplication. This can speed it up by up to a factor of 5.
    Includes another small speed-up if the measurement errors don't include bin-bin correlations.
  - RooUnfoldInverse now works even if there are more truth than measured bins.
    This is done by transposing the response matrix, inverting, and transposing
    back ((A<sup>-1</sup>)<sup>T</sup> = (A<sup>T</sup>)<sup>-1</sup> for pseudo-inverse).
  - Implement `RooUnfoldResponse::Merge` so `hadd` and `TFileMerger` can merge
    response matrix objects.
  - `RooUnfoldResponse::FindBin` is now public to help users fill their own
    response matrix histograms. For consistency, there's also now a trivial 1D
    version.
  - Add `.rootrc` file to define `ACLiC.IncludePaths`, so user doesn't have to
    say `.include src` to run code with ACLiC.

##### Diagnostics

  - Level 2 verbose `RooUnfold::PrintTable` output shows results table with algorithm binning.
  - `RooUnfoldResponse::Print` shows response matrix contents.
  - Matrix debug printout uses new `RooUnfoldResponse::PrintMatrix` method
    (borrowed from `RooUnfoldTestHarness`), which is neater than `TMatrix::Print`.
  - When inverting a matrix (eg. for &chi;<sup>2</sup>), show condition, determinant, and
    check inverse.
  - `RooUnfold::Print()` shows a nicer summary of the object.

##### Testing

  - RooUnfoldTest `seed` parameter allows random number seed to be set. `seed=0` uses a unique seed for each job.
  - `RooUnfoldTest doeff=0 addfakes=1` turn off/on all inefficiencies and fakes.
  - RooUnfoldTest2D+3D projections don't include under/overflow bins.
  - `RooUnfoldTest name=X` parameter allows one to specify a different name for
    the output files.
  - `RooUnfoldTest bincorr=999` tests with full measurement covariance matrix that
    is diagonal (can check it's identical to not using the covariance matrix).

##### Fixes

  - Various fixes for ROOT &le;5.12, 5.32, 5.34, 6.00, 6.08, 6.10, 6.12, 6.14; GCC 4.6; gmake &le;3.80; and Python 3.
  - RooUnfoldResponse::ReplaceAxis: keep axis name when copying from another
    axis. Unfolded histogram was getting the x-axis called 'yaxis' (from truth
    axis of response matrix), which breaks TAxis::SetRangeUser.
    Bug reported by Petr Balek.
  - 2+3D `RooUnfoldResponse::Fill` no longer assumes `_res` has the expected axis
    ranges, since this may be a user's histogram.
  - Protect against null RooUnfoldResponse fakes histogram pointer. This can
    arise if reading an object created with an old version of RooUnfold
    (prior to 1.1.0).
  - Don't use TUnfoldSys for systematic errors with TUnfold if not available
    (ROOT 5.22 - previously didn't have TUnfold, later versions had TUnfoldSys).
  - Fix 3D `RooUnfoldResponse::Fill`.
  - Change `RooUnfold::PrintTable` default to show errors as returned by the
    last unfolding (or else `kErrors`). If kNoError specified, don't show
    error at all. Previously &#8730;N was shown for the error unless specified
    otherwise in the withError parameter. [Bug reported by Stefan Kluth.]
  - Fix PrintTable "total pull". Bug reported by Giovanni Marchiori.
  - RooUnfoldTest3D: fix fake test, change test defaults, plot 2D projections.
  - Fix `RooUnfoldTest overflow=1` or `2` when `nmx`&ne;`ntx`.
  - Fix RooUnfoldErrors binning and &chi;<sup>2</sup> plot.
  - Preserve `Sumw2()` setting in histograms created by `RooUnfoldResponse::H2H1D`,
    `RooUnfoldResponse::HresponseNoOverflow`, and `RooUnfold::HistNoOverflow`.

Version [1.1.1](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/1.1.1)
---
10 October 2011

  - Allow `RooUnfoldResponse(measured,truth,response)` constructor and `Setup`
    method to take `0` or empty histogram for `measured` and/or `truth` in
    the case where there are no fakes and/or inefficiencies.
    With `0`, a 1D histogram is assumed.
  - Calculated histogram of fakes now includes truth under/overflows from user's
    response matrix histogram. Errors are also calculated, even though these are
    not used at present.
  - Fix uninitialsed pointer in `RooUnfoldResponse::Vfakes`. This fixes a crash
    reported by Katharina Bierwagen.
  - Add warning if the number of measured bins doesn't match the
    response matrix.

Version [1.1.0](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/1.1.0)
---
30 September 2011

  - All error calculations, except RooUnfoldTUnfold with error propagation,
    use the full measurement covariance matrix if specified.
    `RooUnfold::SetMeasured` can also be passed vectors for measured distribution
    and errors, or a matrix for the covariance matrix. `SetMeasuredCov` just
    sets the covariance matrix. `GetMeasuredCov` retrieves it.\
    RooUnfoldTest can include test measurement bin-bin correlations, eg.
    `bincorr=-0.5` gives an anti-correlation of 0.5 between neighbouring bins,
    a correlation of 0.25 between next-to-neighbours, etc. This input correlation matrix is also plotted
    for comparison with the unfolded output correlation matrix.
  - RooUnfoldResponse can now account for fakes (or background, ie. measured entries
    with no truth). Use `Fake(xmeas)` to fill, or else add fakes to measured input
    histogram.\
    This can be tested with, eg. `RooUnfoldTest fakexlo=0.1 fakexhi=0.2`
    (and similarly for `y` and `z` for 2D/3D) for linearly varying background level
    between 0.1 at xlo and 0.2 at xhi as a fraction of the number of measured events.
  - Rewrote RooUnfoldBayes to perform the unfolding itself instead of using old classes
    RooUnfoldBayesImpl and Array2D, which have now been removed from the package.
    New `RooUnfoldBayes::UnfoldingMatrix` accessor.
  - Improved RooUnfoldBayes error calculation to take account of dependence of
    the unfolding matrix (_M<sub>ij</sub>_) on the input measurements
    via previous iterations which update the prior (see [this note](bayes_errors.pdf) for details).
    The new errors now agree with the toy MC &mdash; previously they were significantly underestimated.
    In addition, the input measurement errors are now taken from those specified by the user
    (eg. with `SetMeasuredCov` or in the histogram errors)
    instead of assumed multinomial errors based on the number of measured events in each bin
    (as given by D'Agostini). If the user didn't specify any errors, then ROOT's default assumption of
    &#8730;N is used for each bin, which is similar to (and perhaps more appropriate than)
    a multinomial distribution if there are more than a few bins.
  - Allow RooUnfoldSVD to have different numbers of truth and measured bins.
    This is done by calling TSVDUnfold with a symmetric matrices with zeroes
    in the extra bins.
  - The Makefile creates `libRooUnfold.rootmap`, which can be used to load RooUnfold in PyROOT.
    Added `examples/RooUnfoldExample.py`, equivalent to `RooUnfoldExample.cxx`.
    If ROOT was built with PyROOT support and the correct version of Python is used,
    `examples/RooUnfoldExample.py` can be executed as-is.
  - Added `make help` target.
  - `RooUnfold::RunToy()` returns new RooUnfold object.
    `RunToy()` replaces old `Runtoy()` method, which returned a histogram.
    This also fixes a memory leak reported by Seth Zenz.
  - Add I/O streamers for RooUnfold and its subclasses.
  - Fix MacOSX shared library creation in GNUMakefile.
  - Fix problem with different versions of TSVDUnfold that caused a crash in ROOT 5.30.
  - Correct PDF normalisation due to change in RooFit with ROOT 5.29.02
    ([bug #83534](https://savannah.cern.ch/bugs/?83534)) for RooUnfoldTest.
  - Fix `RooUnfold::New` when name but no title is specified.
  - Fix `RooUnfold::PrintTable` for when truth histogram isn't specified.
  - RooUnfoldTest now writes histograms and unfold object to `RooUnfoldTest.root`
    (previously just response object was written there).
  - `RooUnfoldTest ploterrors=1` no longer calculates &chi;<sup>2</sup> distribution.
    This allows the comparison of the errors to be run alone, which can be very much faster.
    Use `ploterrors=2` to include the &chi;<sup>2</sup> plot.
  - RooUnfoldTest plots the unfolded correlation matrix. The full error matrix is printed if `verbose=2`.
  - RooUnfoldTest scales the measured histogram to match the truth if the number of bins is different.

Version [1.0.3](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/1.0.3)
---
14 January 2011

  - RooUnfoldSvd now uses [TSVDUnfold](http://root.cern.ch/root/html/TSVDUnfold.html),
    which was renamed from RooUnfHistoSvd.
    TSVDUnfold is included in ROOT 2.28 and later, but RooUnfold has a copy which
    is used with older versions of ROOT.
  - New RooUnfoldDagostini class. This is an interface to D'Agostini's `bayes.for`
    for comparison with our RooUnfoldBayes (they implement the same algorithm).
    RooUnfoldDagostini is not normally compiled, but will be if `bayes.for` and
    `bayes_c.for` ([download](http://www.roma1.infn.it/~dagos/bayes_distr.txt))
    are copied into the `src` directory.
  - New `RooUnfoldResponse::Add` method, suggested by Seth Zenz.
  - Fixes for building and running on different platforms: [MacOS X](http://www.apple.com/macosx/),
    non-standard RooFit installation, and ROOT versions 5.12 and before.
  - Input histograms can be any TH1, eg. TH1F and TH2F, not just TH1D, TH2D etc.
  - Make sure we don't continue/crash if the unfolding fails.
  - Fix RooUnfoldTUnfold under/overflow bin handling.
  - Reformat inline methods so they appear correctly in the ROOT class documentation.
    Add some missing method documentation comments.
  - Improve RooUnfoldTest generator and bias model.

Version [1.0.2](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/1.0.2)
---
13 September 2010

  - Provide different methods for computing the errors on the unfolded result.
    These are selected with a new `RooUnfold::ErrorTreatment` enum: no error treatment
    (`RooUnfold::kNoError`),
    use bin-by-bin errors (`RooUnfold::kErrors`) or full covariance matrix
    (`RooUnfold::kCovariance`) propagated through the unfolding, or
    covariance matrix from the variation of the results in toy MC tests (`RooUnfold::kCovToy`).
    This last method should be more accurate, especially for RooUnfoldBayes.
  - Added new unfolding algorithm classes: `RooUnfoldTUnfold` is a
    simple (though not yet fully-featured) interface to
    ROOT's [TUnfold](http://root.cern.ch/root/html/TUnfold.html) class (requires ROOT 5.22 and above).
    `RooUnfoldInvert` performs a simple inversion of the response matrix.
  - Simplify `RooUnfoldBinByBin` code (no longer uses `RooUnfoldBayesImpl`).
  - RooUnfold can now include the histogram underflow and overflow bins in the unfolding.
    This is enabled by calling `RooUnfoldResponse::UseOverflow()`.
    This currently only works for 1D histograms.
  - Named static constructor, `RooUnfold::New()`, creates unfolding object
    based on the `RooUnfold::Algorithm` enum (`RooUnfold::kNone`,
    `RooUnfold::kBayes`, `RooUnfold::kSVD`, `RooUnfold::kBinByBin`,
    `RooUnfold::kTUnfold`, or `RooUnfold::kInvert`).
  - New `Clone()` method for `RooUnfold` and its subclasses.
  - New RooUnfold public methods:
    `Chi2(hTrue)` calculates the &chi;<sup>2</sup> of the
    unfolded results with respect to a true distribution,
    `Vmeasured()` and `Emeasured()` return the measured distribution and its errors as vectors,
    `ErecoV()` returns unfolding errors as a vector,
    `SetResponse()` and `SetMeasured()` allow the unfolding inputs to be changed separately,
    `SetRegParm()` and `GetRegParm()` provide a common method to access the regularisation parameter, and
    `SetNToys()` and `NToys()` access the number of MC tests used in error calculation
    with the `RooUnfold::kCovToy` setting.
  - `Impl()` methods return the unfolding implementation object for some algorithms.
  - Added ROOT documentation comments to [RooUnfold](http://roounfold.web.cern.ch/classRooUnfold.html),
    [RooUnfoldErrors](http://roounfold.web.cern.ch/classRooUnfoldErrors.html), and [RooUnfoldParms](http://roounfold.web.cern.ch/classRooUnfoldParms.html).
  - `RooUnfold::PrintTable()` shows also the error, residual, and pull for each bin.
  - New RooUnfoldTest options: `wpaper` and `hpaper` (plot paper width and height),
    `draw=0` (disables histogram drawing), `ploterrors=1` (error analysis using `RooUnfoldErrors`),
    and `plotparms=1` (regularisation parameter analysis using `RooUnfoldParms`),
    `overflow` (`1`=unfolding uses overflows,
    `2`=show under/overflow bins on test histograms).
  - RooUnfoldTest bias model changed to reduce to zero at the edges of the histogram.
  - RooUnfoldTest2D shows projections of 2D pulls.

Version [0.2.2](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.2.2)
---
20 May 2010

  - Fix RooUnfoldBayes for when the measured and truth distributions have
    different numbers of bins. (Problem reported by Lluis Marti.)
  - Fix RooUnfoldSvd to fill `kterm` and `ntoys` correctly (broken in 0.2.1).
  - New accessor methods for unfolding parameters and implementation object.
  - Can control the verbosity of messages with `RooUnfold::SetVerbose(level)`:
    0=warnings, 1=verbose (default, as before), 2=debug, 3=detailed.

Version [0.2.1](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.2.1)
---
19 May 2010

  - Reorganise RooUnfold classes so as not to have hidden virtual methods
    (`Setup` and `Clear`) and to perform unfolding on demand, rather than on construction.
  - Various RooUnfoldTest changes to work round problems in different ROOT versions (especially with CINT).
  - RooUnfoldTest (and 2D, 3D) can now be run from the ROOT prompt with ACLiC (as an alternative to CINT),
    though this requires a bit more setup:

        root [0] gSystem->Load("libRooUnfold")
        root [1] .include src
        root [2] .L examples/RooUnfoldTest.cxx+
        root [3] RooUnfoldTest()

Version [0.2.0](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.2.0)
---
22nd January 2010

  - `RooUnfoldBayesImpl::train` should use data input rather than MC input for
    n(E<sub>j</sub>). This is `_nEstj[]` rather `_nEj[]`. The upshot of this bug was that only
    the final iteration of the Bayesian unfolding had any effect on the result (though a single iteration
    still goes some way!). Problem reported by Jan Kapitan.
  - That change entailed a small reorganisation of `RooUnfoldBayesImpl`.
    Client no longer uses `train()` or `trainBinByBin()` directly.
    This is now done in `unfold()` or
    `unfoldBinByBin()` (which specify the unfolding algorithm parameters,
    `iterations` and `smoothit`), since the training now requires the unfolding
    input. Of course users of `RooUnfoldBayes` and `RooUnfoldBinByBin` won't see any
    difference, since they wrap `RooUnfoldBayesImpl`.
  - Protect against division by zero error in `RooUnfoldBayesImpl::train`.
  - `RooUnfoldBayesImpl::train` should normalise `_nCi` to 1 for the initial P<sub>0</sub>(**C**)
    rather than N<sub>obs</sub>.
  - `RooUnfoldBayes`: fix if there are fewer measured than truth bins.
  - `RooUnfoldBayesImpl::getCovariance` option `doUnfoldSystematic` to enable
    systematic calculation. It remains disabled by default: I'm not sure it is
    correct, it is very slow, and the effect should be small with good MC
    statistics.
  - `RooUnfoldSvd` may not work very well for multi-dimensional distributions,
    so print warning if it's tried.
  - `RooUnfold::PrintTable` improvements: Don't show residuals and pulls for "empty" bins
    (both content=0 and error=0) and don't include them in the &chi;<sup>2</sup>.
    Fix bin numbering (no under/overflow, which aren't included in the unfolding) and show 2D/3D
    bin numbers. Also fixed for different number of measured and truth bins. Print test truth (`hTrue`), which
    is now optional.
  - `RooUnfoldBayesImpl::train`: remove redundant `TStopwatch` timer.
  - Tidy some default object names.
  - Many changes to the `examples/RooUnfoldTest`, `RooUnfoldTest2D`, and new
    `RooUnfoldTest3D`. They now use test harness classes, `RooUnfoldTestHarness`,
    `RooUnfoldTestHarness2D`, and `RooUnfoldTestHarness3D` respectively.
    Test parameters can be specified on the command-line (or ROOT prompt):
    use `RooUnfoldTest -h` or `RooUnfoldTest("-h")` for details.
    New PDFs, which now include a constant background by default. Improved plots.
  - Add `make html` target.

Version [0.1.9](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.1.9)
---
14th October 2009

  - `RooUnfold::PrintTable` prints a table of the results and &chi;<sup>2</sup>/DF.
    This is called from `RooUnfoldTest`.
  - Fix assignment of `RooUnfold` derived classes. Don't repeat `RooUnfold::Setup` when constructing derived classes.
  - Fix bin-by-bin, which was actually running the full bayes code.
  - SVD prints error if `kterm` is negative or greater than the number of bins.
  - Allow weight to be specified with `RooUnfoldResponse::Fill()` and `Miss()`.
    Defaults to 1.0, so no change if not specified. Support for variable binning in 1D case as
    suggested by Seth Zenz.
  - Add `RooUnfoldResponse::ApplyToTruth` based on idea and code from Seth Zenz.
  - Don't calculate covariance matrix unless required. This can be a slow
    operation for the Bayes algorithm. `RooUnfoldBayesImpl` does not call
    `getCovariance` after unfolding &mdash; now done automatically in `RooUnfoldBayes`.
  - Fix conversion of multi-dimensional distributions to 1D (for use in response
    matrix and in Bayes algorithm). Under/overflow bins are ignored.
    Uses new routines `RooUnfoldResponse::GetBin` and `FindBin`, which return global bin for multi-dimensional
    histogram corresponding to vector index or x value.
    Bug reported by Peter Waller.
  - (Fergus Wilson) `RooUnfoldBayesImpl`: added a `getChi2` method. Added 1D smoothing. Speeded
    up covarinace matrix.
  - Fixes for GCC 4.3.3, suggested by Ioana Maris. Add missing include file to `RooUnfoldTest`.
  - `RooUnfoldTest` can test different numbers of bins in truth and measured distributions.
  - `RooUnfoldTest2D` only calculates errors with the Bayes algorithm
    for 25 or fewer bins. It takes a long time for more bins (goes as the 4th power
    of the number of bins).
  - `RooUnfoldTest`: Added a pulls histogram. Draw a line at y=0 on residual plot. Added simple
    checking of some command line parameters.
  - (Fergus Wilson) Changed `RooUnfoldTest` smear method so it now works for different x-axis ranges.
    New test distribution of exponential decaying background and a resonance (i.e. Higgs-like).
  - More robust `GNUmakefile`. If `$ROOTSYS/test/Makefile.arch` doesn't exist, it
    gets settings from `root-config` as suggested by Peter Waller.
    Use `-lRooFitCore` if available (seems to be needed with ROOT 5.18 Cygwin).
  - Use `make ROOTBUILD=debug` for debug build. Allow `make VERBOSE=1` to display compilation commands.

Version [0.1.5](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.1.5)
---
23rd January 2008

  - Include missing headers for ROOT 5.16/00 (also verified to work in 5.18/00).
  - Fix `RooUnfoldTest2D.ps` file name.
  - Rename RooUnfoldSvd's `tau` parameter to `kterm` to match Hoecker and
    Kartvelishvili's usage (they use ***k*** for the last term used in the expansion).
  - Rebuild class documentation with ROOT 5.18/00.

Version [0.1.4](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.1.4)
---
2nd August 2007

##### changes to the SVD algorithm from Kerstin Tackmann

  - Replace use of copy constructors by use of `Clone`.
This should get rid of the segfaults that Jochen has been seeing.
  - Fixed bug in `TUnfHisto::GetCov` that could lead to small
numerical changes in the estimated uncertainties.
  - Tested on my spectrum that the changes, which were supposed to be only
clean-up between the previous two tags (Tim's and mine), do indeed give
identical results on the curvature of the weight distributions and the errors.

Version [0.1.3](https://gitlab.cern.ch/RooUnfold/RooUnfold/tree/0.1.3)
--- 
17th April 2007

First public version.
