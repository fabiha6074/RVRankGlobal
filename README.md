
# RVRankGlobal

<!-- badges: start -->
<!-- badges: end -->

The goal of RVRankGlobal is to implement prioritization of rare genetic variants related to the disease using penalized likelihood method.The method utillizes familial transmission patterns, function status of the variants and unusualness score defined by age-at-onset of disease.

## Installation

You can install the development version of RVRankGlobal like so:

``` r
install.packages("remotes")
remotes::install_github("fabiha6074/RVRankGlobal")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(RVRankGlobal)
## data <- make_data(...)
fit <- fit_model(data)
post <- variant_posteriors(fit, data)
results <- summarize_rv(fit, post, data)
```

