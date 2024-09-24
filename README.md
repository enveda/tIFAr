
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tIFA

The tIFA R package presents a set of functions that perform missing data
imputation using a truncated infinite factor analysis model designed for
high-dimensional metabolomics data.

## Installation

You can install the development version of tIFA like so:

``` r
# install remotes package if not already installed
# install.packages("remotes")
remotes::install_github("kfinucane/tIFA")

# load the library
library(tIFA)
```

## Example imputation

To illustrate the tIFA functionality, here a test dataset is generated
with two missing values.

``` r
# load the library
library(tIFA)

# generate example data
example_data <- matrix(abs(rnorm(100)), nrow = 5)

# add missingness to example data coded as 0.001
example_data[4, 2] <- 0.001
example_data[2, 18] <- 0.001
```

This example data can then be passed into the `tIFA_model()` function
which contains all package functionality. Here, `input_data` is the
dataset with missing values, in matrix format. The `coding` argument
refers to how the missing values are coded. If your missing values all
have a value of `NA`, for example, you should input `coding = NA`. Here,
for the example, missing entries are coded with the value `0.001`. Here,
the value `k.star = 3` refers to the practical non-infinite number of
latent factors used by the tIFA model; usually this defaults to
`k.star = 5` but for this small example dataset we will use a smaller
number.

The remaining parameters refer to the MCMC procedure; here
`n.iters = 100` runs the MCMC chain for 100 iterations, with `burn = 50`
controlling the number of draws discarded in a burn and `thin = 5`
stating that every fifth post-burn draw from the MCMC chain should be
retained.

``` r
# run tIFA model
# short chain for example
res <- tIFA_model(input_data = example_data, coding = 0.001, n.iters = 100, k.star = 3,
                  burn = 50, thin = 5)
```

The results of the tIFA imputation can be accessed as follows.

``` r
# checkout imputed dataset
res$imputed_dataset

# checkout further information on imputed entries
res$imputation_info
```

If one wishes to change the hyperparameters of the tIFA model from their
defaults, this can be done as follows. All parameters are as given in
the tIFA paper.

``` r
res <- tIFA_model(input_data = example_data, coding = 0.001, n.iters = 100, k.star = 3,
                  burn = 50, thin = 5, mu_varphi = 0.1, kappa_1 = 3L, kappa_2 = 2L,
                  a_sigma = 1L, b_sigma = 0.3, a_1 = 2.1, a_2 = 3.1)
```
