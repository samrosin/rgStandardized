#' Compute standardized Rogan-Gladen estimate
#'
#' @description
#' Computes nonparametric standardized Rogan-Gladen prevalence estimator.
#' Corresponds to \eqn{\hat{\pi_{SRG}}} and \eqn{\hat V_{\pi, SRG}} in the manuscript.
#'
#' @param sample A dataframe with data from the sample, including the stratum_props (\eqn{\gamma_j} values) which will be used in standardization
#' @param stratum_props The stratum_props dataframe
#' @param sigma_e_hat Estimated sensitivity
#' @param sigma_p_hat Estimated specificity
#' @param n_1 Sample size for sensitivity validation dataset
#' @param n_2 Sample size for specificity validation dataset
#' @param n_3 Sample size for main study from population
#' @param vars_std The names of variables to standardize over, corresponding to columns in sample
#' @param variance Should variance be estimated?
#'
#' @return A list of estimated prevalence and variance; optionally,
#'         if variance is FALSE, the truncated prevalence estimate
#' \itemize{
#'  \item hat_pi_srg - Estimated standardized Rogan-Gladen prevalence
#'  \item hat_var_pi_srg - Estimate of the variance of the standardized Rogan-Gladen estimator
#'  \item num_strata - Number of strata in the sample used for standardization
#' }
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom rlang .data
#'
#' @examples
#'
#'  serology <- belgium$serology
#'  belg_cr <- serology[serology$collection_round == 1,]
#'  belg_cr$x <- ifelse(belg_cr$igg_cat=="positive", 1, 0)
#'  belg_cr_srg <- as.data.frame(belg_cr)
#'
#'  srg(belg_cr_srg, belgium$strataprops, 154/181, 322/326,
#'       181, 326, nrow(belg_cr_srg), vars_std = c("age_cat", "sex"), variance = TRUE)
#'
#' @export
#'
srg <- function(sample, stratum_props, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, vars_std, variance = TRUE){

  # join stratum proportions
  sample <- left_join(sample, stratum_props, by = vars_std)

  #unique covariate strata in the sample data
  strata <- unique(sample[,c(vars_std,"stratum_prop")])

  #the stratum proportions (\gamma_js) have to first be redefined/recomputed to match the target dataset
  strata$stratum_prop <- strata$stratum_prop / sum(strata$stratum_prop)

  #strata, with number and proportion of positive tests
  strata_npos <- sample %>% dplyr::group_by_at(vars_std) %>%
    dplyr::summarise(n = n(), n_pos = sum(.data$x), .groups = "drop")

  #join the two datasets and make the standardization estimates, then sum them.
  #the formula is std_est = \hat \rho_j * \gamma_j
  strata <- dplyr::inner_join(strata, strata_npos, by = vars_std)

  strata$std_est <- (strata$n_pos / strata$n) * strata$stratum_prop
  pi_hat_st <- (sum(strata$std_est) - (1-sigma_p_hat))/(sigma_e_hat - (1-sigma_p_hat)) #point estimate
  if(variance == FALSE){
    return(c(truncate_01(pi_hat_st), nrow(strata))) # include truncated estimate and number of strata in return vector
  } else{
    # compute variance estimator by summing the three components. compare to formula in manuscript.
    a <- pi_hat_st^2 * sigma_e_hat * (1 - sigma_e_hat) / n_1
    b <- (1 - pi_hat_st)^2 * sigma_p_hat * (1 - sigma_p_hat) / n_2
    # c is the third component of the variance estimator. c = (gamma_j^2*\hat \rho_j * (1- \hat \rho_j))/(n_{z_j})
    strata <- strata %>% mutate(
      c = .data$stratum_prop^2 * (.data$n_pos / .data$n) * (1 - (.data$n_pos / .data$n) ) / .data$n
    )
    c <- sum(strata$c)
    var_pi_hat_st <- (a + b + c) * (sigma_e_hat + sigma_p_hat - 1)^(-2)
  }

  # return point and variance estimates
  # note point estimate is not truncated for use in CI construction
  return_list <- list(
    "hat_pi_srg" = pi_hat_st,
    "hat_var_pi_srg" = var_pi_hat_st,
    "num_strata" = nrow(strata)
  )
  return_list
}
