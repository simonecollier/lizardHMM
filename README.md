
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lizardHMM

<!-- badges: start -->
<!-- badges: end -->

The goal of lizardHMM is to fit lizard movement time series data,
composed of step-lengths per second, with hidden Markov models and
investigate the quality of fit that arises. This package can work with
other time series data including simulated data from the package itself.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simonecollier/lizardHMM")
```

## Key components

-   Simulate data from an HMM
    -   `norm_generate_sample()` produces data from an n-state HMM with
        the desired normal state dependent distributions and transition
        probabilities.
-   Fit data with an HMM
    -   `norm_fit_hmm()` and `gam_fit_hmm()` both work to fit an n-state
        HMM with  
        normal/gamma state dependent distributions to the data given.
    -   `gam0_fit_hmm()` fits the data with an n-state HMM with gamma
        state dependent distributions and includes a point mass on zero
        for the state with the smallest mean.
    -   `norm_viterbi()`, `gam_viterbi()`, and `gam0_viterbi()` use
        global decoding to find the most likely sequence of states that
        could have generated the data according the HMM that was fit.
-   Check fit
    -   `norm_ci()`, `gam_ci()`, and `gam0_ci()` produce confidence
        intervals for each of the fitted parameters.
    -   `norm_forecast_psr()`, `gam_forecast_psr()`, and
        `gam0_forecast_psr()` compute the normal forecast
        pseudo-residuals for the data fitted with the HMM.
-   Visualizations
    -   `timeseries_plot()` plots the time series data with colors
        corresponding to the states decoded by the viterbi algorithm.
    -   `norm_hist_ci()`, `gam_hist_ci()`, and `gam0_hist_ci()` plot the
        histogram of the data with the fitted state dependent
        distributions overlayed and their corresponding confidence
        intervals.
    -   The functions in `psr_plotting.R` produce visualizations of the
        normal forecast pseudo-residuals.

## Extra details

All of the distributions types are set up to handle multiple subjects,
variables, and covariates for the transition probabilities. The only
option is complete pooling of parameters when working with multiple
subjects, although there may be updates in the future that include more
options. The functions in `covariate_analysis.R` can be used to
investigate the effect of a single covariate on the transition
probabilities and stationary distribution.
