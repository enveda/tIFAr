
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tIFA

The tIFA R package presents a set of functions that perform missing data
imputation using a truncated infinite factor analysis model designed for
high-dimensional metabolomic data.

## Installation

You can install the development version of tIFA like so:

``` r
# install remotes package if not already installed
# install.packages("remotes")
remotes::install_github("kfinucane/tIFA")

# Load the library
library(tIFA)
```

## Example imputation

To illustrate the tIFA functionality, here a test dataset is generated
with two missing values.

``` r
library(tIFA)

# generate example data
example_data <- matrix(abs(rnorm(100)), nrow = 5)

# add missingness to example data coded as 0.001
example_data[4, 2] <- 0.001
example_data[2, 18] <- 0.001
```

The data are then processed using the FA_process function, which applies
relevant pre-processing to the data in preparation for tIFA imputation.

``` r
# apply data pre-processing function
input_data <- FA_process(example_data, coding = 0.001)
```

Zero-imputation is performed to pass as starting values to the tIFA
model.

``` r
# perform zero imputation to use as initial imputed data value in tIFA
# this is for example purposes and is not generally recommended
zero_imp <- input_data$data
zero_imp[is.na(zero_imp)] <- 0
```

The pre-processed data are then passed to the tIFA model function.

``` r
# run tIFA model
# here a very short chain length is used as an example
tifa_res <- tIFA_model(data = input_data$data, missingness = input_data$missingness, M = 300,
                      initial_dataset = zero_imp, k = 5)
```

The results are then post-processed to provide the final imputed dataset
with accompanying imputation information.

``` r
# process tIFA results
# here a very small burn is used as an example
processed_res <- tifa_post(tifa_res, burn = 60)
```
