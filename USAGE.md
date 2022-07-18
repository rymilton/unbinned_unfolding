Unfolding Method
---

We use unfolding to remove the known effects of measurement
resolutions, systematic biases, and detection efficiency to determine
the "true" distribution. We parametrise the measurement effects using
a response matrix that maps the (binned) true distribution onto the
measured one. For 1-dimensional true and measured distribution bins Tj
and Mi, the response matrix element Rij gives the fraction of events
from bin Tj that end up measured in bin Mi. The response matrix is
usually determined using Monte Carlo simulation (training), with the
true values coming from the generator output.

The unfolding procedure reconstructs the true Tj distribution from the
measured Mi distribution, taking into account the measurement
uncertainties due to statistical fluctuations in the finite measured
sample (without these uncertainties, the problem could be solved
uniquely by inverting the response matrix). RooUnfold provides several
algorithms for solving this problem.

The iterative and SVD unfolding algorithms require a regularisation
parameter to prevent the statistical fluctuations being interpreted as
structure in the true distribution. It is therefore necessary to
optimise this parameter for the number of bins and sample size, using
Monte Carlo samples of the same size as the data. These samples can
also be used to measure the effectiveness of the unfolding and hence
provide estimates of the systematic errors that result from the
procedure (testing).

Note that for this last step (in particular), it is important to use
Monte Carlo samples with truth distributions that are statistically
and systematically independent of the sample used in training (such
samples would anyway be used in a systematics analysis, eg. using a
different generator, or reweighting variations within the a-priori
uncertainties in the truth distribution). After all, if the Monte
Carlo were a perfect model of the data, we could use the Monte Carlo
truth information directly and dispense with unfolding altogether!

The bin-by-bin method assumes no migration of events between bins
(eg. resolution is much smaller than the bin size and no systematic
shifts). This is of course trivial to implement without resorting to
the RooUnfold machinery, but is included in the package to allow
simple comparison with the other methods.  

Using RooUnfold
---

To use RooUnfold, we must first supply the response matrix object
RooUnfoldResponse. It can be constructed like this:

    RooUnfoldResponse response (nbins, x_lo, x_hi);

or, if different truth and measured binning is required,

    RooUnfoldResponse response (nbins_measured, x_lo_measured, x_hi_measured,
                                nbins_true,     x_lo_true,     x_hi_true);

or, if different binning is required,

    RooUnfoldResponse response (hist_measured, hist_truth);

In that last case, `hist_measured` and `hist_truth` are used to specify
the dimensions of the distributions (the histogram contents are not
used here), eg. for 2D or 3D distributions or non-uniform binning.

This `RooUnfoldResponse` object is often most easily filled by looping
over the training sample and calling `response.Fill(x_measured,x_true)`
or, for events that were not measured due to detection inefficiency,
`response.Miss(x_true)`

    if (measurement_ok)
      response.Fill (x_measured, x_true);
    else
      response.Miss (x_true);

Alternatively, the response matrix can be constructed from a
pre-existing `TH2D` 2-dimensional histogram (with truth and measured
distribution `TH1D` histograms for normalisation).

This response object can be passed directly to the unfolding object,
or written to a ROOT file for lateruse.

To perform the unfolding (either to try different regularisation
parameters, for testing, or for real data), create a `RooUnfold` object
and pass it the test / measured distribution (as a histogram) and the
response object.

    RooUnfoldBayes    unfold (&response, hist_measured, iterations);

or

    RooUnfoldSvd      unfold (&response, hist_measured, kterm);

or

    RooUnfoldBinByBin unfold (&response, hist_measured);

Here, `hist_measured* is a pointer to a `TH1D` (or `TH2D` for the 2D case)
histogram of the measured distribution (it should have the same
binning as the response matrix). The classes `RooUnfoldBayes`,
`RooUnfoldSvd`, `RooUnfoldBinByBin`, etc. all inherit from `RooUnfold`
and implement the different algorithms. The integer `iterations` (for
`RooUnfoldBayes`) or `kterm` (`RooUnfoldSvd`) is the regularisation
parameter. (Note that `RooUnfoldSvd`'s `kterm` parameter is also known
as `tau` in the code. That usage is incompatible with the literature,
so we adopt `k` here.)

The reconstructed truth distribution (with errors) can be obtained
with the `Hunfold()` method.

    TH1D* hist_unfold= (TH1D*) unfold.Hunfold();

The result can also be obtained as as a `TVectorD` with full
`TMatrixD` covariance matrix.

Multi-dimensional distributions can also be unfolded, though this does
not work for the SVD method, and the interface is rather clumsy (we
hope to improve this).

See the [class documentation](https://roounfold.web.cern.ch) for
details of the `RooUnfold` and `RooUnfoldResponse public` methods.

A very simple example of RooUnfold's use is given in
`examples/RooUnfoldExample.cxx`. More complete tests, using different
toy MC distributions, are in `examples/RooUnfoldTest.cxx` and
`examples/RooUnfoldTest2D.cxx`.

Choice of Regularisation Parameter
---

The regularisation parameter determines the relative weight placed on
the data, compared to the training sample truth. Both RooUnfoldBayes
and RooUnfoldSvd take integer regularisation parameters, with small
values favouring the MC truth and larger values favouring the data.

`RooUnfoldBayes` Regularisation
---

For RooUnfoldBayes, the regularisation parameter specifies the number
of iterations, starting with the training sample truth
(`iterations=0`). You should choose a small integer greater than `0`
(we use `4` in the examples). Since only a few iterations are needed,
a reasonable performance can usually be obtained without fine-tuning
the parameter.

The optimal regularisation parameter can be selected by finding the
largest value up to which the errors remain reasonable (ie. do not
become much larger than previous values). This will give the smallest
systematic errors (reconstructed distribution least biased by the
training truth), without too-large statistical errors. Since the
statistical errors grow quite rapidly beyond this point, but the
systematic bias changes quite slowly below it, it can be prudent to
reduce the regularisation parameter a little below this optimal point.
RooUnfoldSvd Regularisation

For `RooUnfoldSvd`, the unfolding is something like a Fourier
expansion in "result to be obtained" vs "MC truth input". Low
frequencies are assumed to be systematic differences between the
training MC and the data, which should be retained in the output. High
frequencies are assumed to arise from statistical fluctuations in data
and unfortunately get numerically enhanced without proper
regularization. Choosing the regularization parameter, k (kterm),
effectively determines up to which frequencies the terms in the
expansion are kept. (Actually, this is not quite true, we don't use a
hard cut-off but a smooth one.)

The correct choice of k is of particular importance for the SVD
method. A too-small value will bias the unfolding result towards the
MC truth input, a too-large value will give a result that is dominated
by unphysically enhanced statistical fluctuations. This needs to be
tuned for any given distribution, number of bins, and approximate
sample size — with `k` between `2` and the number of bins. (Using
`k=1` means you get only the training truth input as result without
any corrections. You basically regularise away any differences, and
only keep the leading term which is, by construction, the MC truth
input.)

[Höcker and Kartvelishvili's
paper](https://arxiv.org/abs/hep-ph/9509307) (section 7) describes how
to choose the optimum value for `k`.

Running
--------

In an interactive ROOT shell, `RooUnfoldExample.cxx` makes a simple test of RooUnfold.

    % root
    root [0] .L examples/RooUnfoldExample.cxx
    root [1] RooUnfoldExample()

You can also use python

    % python examples/RooUnfoldExample.py
  
The example programs can also be run from the shell command line.

    ./RooUnfoldTest

and the output appears in `RooUnfoldTest.pdf`.

You can specify parameters for `RooUnfoldTest` (either as an argument to
the routine, or as parameters to the program), eg.

    % root
    RooUnfoldTest("method=2 ftestx=3")

or

    ./RooUnfoldTest method=2 ftestx=3

Use `RooUnfoldTest -h` or `RooUnfoldTest("-h")` to list all the parameters and their defaults.

`method` specifies the unfolding algorithm to use:

    0	  no unfolding (output copied from measured input)
    1	  Bayes
    2	  SVD
    3	  bin-by-bin
    4	  TUnfold
    5	  matrix inversion 

`ftrainx` and `ftestx` specify training and test PDFs:

    0	  flat distribution
    1	  top-hat distribution over mean ± 3 x width
    2	  Gaussian
    3	  Double-sided exponential decay
    4	  Breit-Wigner
    5	  Double Breit-Wigner, peaking at mean-width and mean+width
    6	  exponential
    7	  Gaussian resonance on an exponential background 

The centre and width scale of these signal PDFs can be specified with
the `mtrainx` and `wtrainx` (and `mtestx` and `wtestx`) parameters
respectively. A flat background of fraction `btrainx` (and `btestx`) is
added. Detector effects are modelled with a variable shift (`xbias` in
the centre, dropping to `0` near the edges), a smearing of `xsmear` bins,
as well as a variable efficiency (between `effxlo` at `xlo` and `effxhi` at
xhi).

For 2D and 3D examples look at RooUnfoldTest2D and RooUnfoldTest3D.

Limitations
---

   - RooUnfoldSvd does not work correctly for multi-dimensional distributions (gives a warning).