
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

If every model is given the same weights, most of the time instrumental variables will be included in the propensity score since instrumental variable is associated to the exposure variable. However, the literature has shown that including these variables might increase the variance and amplify the bias of the estimate (Pearl 2010 and Brookhart et al. 2006)
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

## References

Brookhart MA, Schneeweiss S, Rothman KJ, Glynn RJ, Avorn J, Stürmer T. Variable selection for propensity score models. Am J Epidemiol. 2006 Jun 15;163(12):1149-56.

Pearl, J., P. Grunwald, and P. Spirtes. "Proceedings of the Twenty-Sixth Conference on Uncertainty in Artificial Intelligence (UAI 2010)." (2010): 417-24.
