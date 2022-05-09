
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgStandardized

<!-- badges: start -->
<!-- badges: end -->

`rgStandardized` implements standardized prevalence estimators that
correct for measurement error and selection bias. The estimators
leverage validation data on test sensitivity and specificity (false
positives and negatives) to correct for measurement error, following the
method of [Rogan and Gladen
(1978)](https://pubmed.ncbi.nlm.nih.gov/623091/). Assuming equal
probability of sample inclusion within strata, the estimators leverage
target population data to correct for selection bias in non-probability
samples. `rgStandardized` implements a nonparametric and logistic
regression estimator.

## Installation

The development version of `rgStandardized` can be installed with:

``` r
# install.packages("devtools")
library(devtools)
#> Loading required package: usethis
load_all()
#> ℹ Loading rgStandardized
#devtools::install_github("samrosin/rgStandardized")
```

## Citation

Rosin, S., Shook-Sa, B. E., Cole, S. R., and Hudgens, M. G. Estimating
SARS-CoV-2 Seroprevalence (2022). arXiV.

## Example

An example is a seroprevalence survey conducted in Belgium in 2020 to
estimate the prevalence of SARS-CoV-2 IgG antibodies. We estimate
seroprevalence in the first of seven collection rounds. The Rogan-Gladen
estimator corrects for measurement error (misclassification bias), but
not for selection bias. Note that *σ̂*<sub>*e*</sub> and
*σ̂*<sub>*p*</sub> are estimated sensitivity and specificity from sample
sizes of *n*<sub>1</sub> and *n*<sub>2</sub>, and that *ρ̂* is the sample
proportion of positives from the main survey of size *n*<sub>3</sub>.
The study used a stratified random sample, so not much selection bias is
expected.

``` r
library(rgStandardized)

serology <- belgium$serology
 belg_cr1 <- serology[serology$collection_round == 1,] # collection round 1
belg_cr1$x <- ifelse(belg_cr1$igg_cat=="positive", 1, 0) # data cleaning
 
rg_ests <- rg(rho_hat = sum(belg_cr1$x)/nrow(belg_cr1),
               154/181, 322/326, 181, 326, nrow(belg_cr1), TRUE)
 
paste(round(rg_ests$hat_pi_rg, 4),
      round(rg_ests$hat_var_pi_rg, 4))
#> [1] "0.0159 1e-04"
```

Not correcting for measurement error, the seroprevalence estimate is
0.016 with a variance estimate of 0.

The methods developed in the manuscript adjust for both measurement
error and selection bias. Estimates adjusting for age category and sex
are as follows:

``` r
library(rgStandardized)

srg_ests <- srg(belg_cr1, belgium$strataprops, 154/181, 322/326,
      181, 326, nrow(belg_cr1), vars_std = c("age_cat", "sex"), variance = TRUE)

paste(round(srg_ests$hat_pi_srg, 4),
      round(srg_ests$hat_var_pi_srg, 4))
#> [1] "0.0187 1e-04"

# logistic regression with interaction between age and sex
srgm_ests <- srgm(belg_cr1, belgium$strataprops, 154/181, 322/326,
      181, 326, nrow(belg_cr1), vars_std = c("age_cat", "sex"),
      mod_formula = formula("x ~ sex + age_cat + sex*age_cat"), variance = TRUE)

paste(round(srgm_ests$hat_pi_srgm, 4),
      round(srgm_ests$hat_var_pi_srgm, 4))
#> [1] "0.0187 1e-04"
```

After accounting for selection bias, the standardized seroprevalence
estimates are higher.
