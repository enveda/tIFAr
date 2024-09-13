
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

This example data can then be passed into the `tIFA_model` function
which contains all package functionality. Here, `input_data` is the
dataset with missing values, in matrix format. The `coding` argument
refers to how the missing values are coded. If your missing values all
have a value of `NA`, for example, you should input `coding = NA`. Here,
for the example, missing entries are coded with the value `0.001`. Here,
the value `k = 3` refers to the number of effective factors assumed by
the tIFA model; usually this defaults to `k=5` but for this small
example dataset we will use a smaller number. The remaining parameters
refer to the MCMC procedure; this example runs a very short chain for
illustration purposes.

``` r
# run tIFA model
# short chain for example
tifa_res<- tIFA_model(input_data = example_data, coding = 0.001, k = 3, M = 100,
                      burn = 50, thin = 5)
```
