
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BACps

<!-- badges: start -->

<!-- badges: end -->
## Bayesian Adjustment for Confounding (BAC) in Bayesian Propensity Score Estimation

The goal of BACps is to estimate the average causal effect accounting for two sources of uncertainty:

<ul>
  <li> uncertainty regards the propensity score. The propensity score is a quantity that we do not observe and thus we have to estimate. So, the idea is to account for the fact that we are not using the true propensity score and thus we can be making a mistake  </li>
  <li> uncertainty regards the model. This is related to the uncertainty that we face when we decide the variables that we include in the model. Instead of fixing one model associated with one set of covariates/features, we consider all possible models. </li>

</ul>

## Installation

You can install the package version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("pablogonzalezginestet/BACps")
```

## Example

See the vignette for details: [online
vignette](https://pablogonzalezginestet.github.io/BACps/)
