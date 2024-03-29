% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srgm.R
\name{srgm}
\alias{srgm}
\title{Compute logistic regression standardized Rogan-Gladen estimates}
\usage{
srgm(
  sample,
  stratum_props,
  sigma_e_hat,
  sigma_p_hat,
  n_1,
  n_2,
  n_3,
  vars_std,
  mod_formula,
  variance = FALSE
)
}
\arguments{
\item{sample}{A dataframe with data from the sample, including the stratum_props (\eqn{\gamma_j} values) which will be used in standardization}

\item{stratum_props}{the stratum_props dataframe. note that the vars_std need to be the first columns from the left}

\item{sigma_e_hat}{Estimated sensitivity}

\item{sigma_p_hat}{Estimated specificity}

\item{n_1}{Sample size for sensitivity validation dataset}

\item{n_2}{Sample size for specificity validation dataset}

\item{n_3}{Sample size for main study from population}

\item{vars_std}{The names of variables to standardize over, corresponding to columns in sample}

\item{mod_formula}{An object of class "formula" that symbolically describes the model to be fitted}

\item{variance}{Should variance be estimated? (default TRUE)}
}
\value{
A list of estimated prevalence and variance; optionally,
if variance is FALSE, the truncated prevalence estimate
\itemize{
\item hat_pi_srgm - Estimated model-based standardized Rogan-Gladen prevalence
\item hat_var_pi_srgm - Estimate of the variance of the model-based standardized Rogan-Gladen estimator
}
}
\description{
Logistic regression outcome model for the mismeasured
stratum-specific prevalences \eqn{\rho_j}.
Corresponds to \eqn{\hat{\pi_{SRGM}}} and \eqn{\hat V_{\pi, SRGM}} in the manuscript.
}
\examples{

 serology <- belgium$serology
 belg_cr <- serology[serology$collection_round == 1,]
 belg_cr$x <- ifelse(belg_cr$igg_cat=="positive", 1, 0)
 belg_cr_srgm <- as.data.frame(belg_cr)

 srgm(belg_cr_srgm, belgium$strataprops, 154/181, 322/326,
      181, 326, nrow(belg_cr_srgm), vars_std = c("age_cat", "sex"),
      mod_formula = formula("x ~ sex + age_cat + sex*age_cat"), variance = TRUE)

}
