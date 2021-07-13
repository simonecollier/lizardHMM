
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lizardHMM

<!-- badges: start -->

<!-- badges: end -->

The goal of lizardHMM is to fit lizard movement data with Hidden Markov
Models.

## Installation

You can install the released version of lizardHMM from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("lizardHMM")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("simonecollier/lizardHMM")
```

## Example

Set up a basic HMM with data distributed according to normal state
dependent distributions. This HMM will have 3 states, 2 variables, 1
subject, and 2 covariates (included in the computation of the transition
probability matrix). We also define the design matrix which indicates
the values of the covariates at each point in time. Then we go about
generating data from this HMM.

``` r
library(lizardHMM)

hmm1 <- list(num_states = 3,
             num_variables = 2,
             num_subjects = 1,
             mu = list(matrix(c(12, 18, 22), ncol = 3, nrow = 1, byrow = TRUE), 
                       matrix(c(-12, -7, 0), ncol = 3, nrow = 1, byrow = TRUE)),
             sigma = list(matrix(c(3, 1, 1.5), ncol = 3, nrow = 1, byrow = TRUE), 
                          matrix(c(1, 3, 2), ncol = 3, nrow = 1, byrow = TRUE)),
             beta  = matrix(c(0.01, 0.02, 0.001,
                              0.01, 0.03, 0.004,
                              0.01, 0.01, 0.003,
                              0.01, 0.04, 0.002,
                              0.01, 0.01, 0.004,
                              0.01, 0.03, 0.001), ncol = 3, nrow = 6, byrow = TRUE),
             delta = list(c(0.3, 0.02, 0.5)))

num_sample <- 1000

design           <- list(matrix(0, nrow = 1000, ncol = 3))
design[[1]][, 1] <- 1 # First column is the intercept
design[[1]][, 2] <- sample(c(0, 1), size = 1000, 
                           prob = c(0.3, 0.7), replace = TRUE)
design[[1]][, 3] <- rnorm(1000, mean = 5, sd = 1)

sample <- norm_generate_sample(num_sample, hmm1, design)
x      <- sample$observ
```

We can plot the timeseries of the observations for each subect and
variable along with the corresponding states.

``` r
timeseries_plot(x, sample$state, hmm1$num_subjects, hmm1$num_variables)
#> [[1]]
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

    #> 
    #> [[2]]

<img src="man/figures/README-unnamed-chunk-2-2.png" width="100%" />

The histogram of the data for each subject and variable can be gnerated,
overlayed with the state dependent normal
distributions.

``` r
norm_hist(sample, hmm1$num_states,  hmm1$num_subjects, hmm1$num_variables, 
          hmm1, width = 1, x_step = 0.2)
#> [[1]]
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

    #> 
    #> [[2]]

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

Now we can try fitting the data we generated with a basic HMM with
reasonable guesses for the initial parameters.

``` r
num_states = 3
num_variables = 2
num_subjects = 1
num_covariates = 2
mu0 <- list(matrix(c(11, 19, 23), ncol = 3, nrow = 1, byrow = TRUE),
            matrix(c(-10, -5, 1), ncol = 3, nrow = 1, byrow = TRUE))
sigma0 <- list(matrix(c(3, 3, 3), ncol = 3, nrow = 1, byrow = TRUE),
               matrix(c(3, 3, 3), ncol = 3, nrow = 1, byrow = TRUE))
beta0 <- matrix(c(-2, 0, 0,
                  -2, 0, 0,
                  -2, 0, 0,
                  -2, 0, 0,
                  -2, 0, 0,
                  -2, 0, 0), ncol = 3, nrow = 6, byrow = TRUE)
delta0 <- list(c(1/3, 1/3, 1/3))

hmm_fit <- norm_fit_hmm(x, design, num_states, num_variables, num_subjects,
                        num_covariates, mu0, sigma0, beta0, delta0,
                        iterlim = 200, hessian = TRUE)
```

We can find the confidence intervals for our parameter estimates.

``` r
conf_intervals <- norm_ci(hmm_fit, state_dep_dist_pooled = FALSE, 
                          n = 100, level= 0.975) 
conf_intervals
#> $mu
#> $mu$estimate
#> $mu$estimate$`1`
#>         [,1]     [,2]     [,3]
#> [1,] 11.8435 17.92814 21.87435
#> 
#> $mu$estimate$`2`
#>           [,1]      [,2]       [,3]
#> [1,] -12.03413 -7.057236 -0.1166303
#> 
#> 
#> $mu$upper
#> $mu$upper$`1`
#>          [,1]     [,2]     [,3]
#> [1,] 12.21659 18.04655 22.05237
#> 
#> $mu$upper$`2`
#>           [,1]      [,2]       [,3]
#> [1,] -11.93159 -6.629203 0.04621965
#> 
#> 
#> $mu$lower
#> $mu$lower$`1`
#>          [,1]     [,2]     [,3]
#> [1,] 11.33377 17.79719 21.68645
#> 
#> $mu$lower$`2`
#>           [,1]      [,2]       [,3]
#> [1,] -12.12582 -7.452267 -0.3301668
#> 
#> 
#> 
#> $sigma
#> $sigma$estimate
#> $sigma$estimate$`1`
#>          [,1]     [,2]     [,3]
#> [1,] 2.805934 0.973412 1.470498
#> 
#> $sigma$estimate$`2`
#>         [,1]     [,2]     [,3]
#> [1,] 1.01955 2.750253 1.900116
#> 
#> 
#> $sigma$upper
#> $sigma$upper$`1`
#>          [,1]     [,2]     [,3]
#> [1,] 3.050136 1.054169 1.607198
#> 
#> $sigma$upper$`2`
#>          [,1]     [,2]     [,3]
#> [1,] 1.101603 3.009897 2.013362
#> 
#> 
#> $sigma$lower
#> $sigma$lower$`1`
#>          [,1]      [,2]    [,3]
#> [1,] 2.481494 0.9103495 1.34863
#> 
#> $sigma$lower$`2`
#>           [,1]     [,2]     [,3]
#> [1,] 0.9498778 2.505115 1.756415
#> 
#> 
#> 
#> $beta
#> $beta$estimate
#>             [,1]         [,2]        [,3]
#> [1,] -0.79685177  0.222763919  0.13535118
#> [2,] -1.97397817 -0.147895603  0.44100573
#> [3,]  0.09928201 -0.005353139  0.01007736
#> [4,]  0.17429201  0.013427091 -0.01592751
#> [5,] -0.18218420  0.183623444  0.02924870
#> [6,]  0.57678238 -0.151758725 -0.10271142
#> 
#> $beta$upper
#>            [,1]      [,2]      [,3]
#> [1,]  0.5959862 0.8908671 0.4630644
#> [2,] -0.4004136 0.3677746 0.7131303
#> [3,]  1.6727627 0.6199967 0.2422721
#> [4,]  1.7875848 0.5789854 0.2501439
#> [5,]  1.1355597 0.6278717 0.2747697
#> [6,]  2.0300227 0.2683707 0.1929391
#> 
#> $beta$lower
#>            [,1]       [,2]       [,3]
#> [1,] -2.2069678 -0.2916736 -0.1569097
#> [2,] -3.3313744 -0.8403146  0.1284087
#> [3,] -1.4692445 -0.6474915 -0.2887972
#> [4,] -1.6926678 -0.7082243 -0.3294184
#> [5,] -1.0688958 -0.4037150 -0.2339808
#> [6,] -0.7599779 -0.6564322 -0.3485810
```

Using the viterbi algorithm we can decode the hidden states according to
our fitted HMM and plot the resulting time series.

``` r
viterbi <- norm_viterbi(x, hmm_fit)
timeseries_plot(x, viterbi, num_subjects, num_variables)
#> [[1]]
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

    #> 
    #> [[2]]

<img src="man/figures/README-unnamed-chunk-6-2.png" width="100%" />

We can compute the confidence intervals of the state dependent normal
distributions and compare them to the histogram of
data.

``` r
norm_hist_ci(x, viterbi, num_states, num_subjects, num_variables, hmm_fit, 
             state_dep_dist_pooled = FALSE,
             width = 1, n = 100, level = 0.975, x_step = 0.2)
#> [[1]]
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

    #> 
    #> [[2]]

<img src="man/figures/README-unnamed-chunk-7-2.png" width="100%" />
