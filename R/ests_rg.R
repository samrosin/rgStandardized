########
# Implements the estimators proposed in the manuscript

#' Compute Rogan-Gladen estimate
#'
#' The estimator corresponds to the point and variance estimators
#' from Rogan and Gladen (1978)
#'
#' @param rho_hat Raw/apparent prevalence (positives divided by total samples)
#' @param sigma_e_hat Estimated sensitivity
#' @param sigma_p_hat Estimated specificity
#' @param n_1 Sample size for sensitivity validation dataset
#' @param n_2 Sample size for specificity validation dataset
#' @param n_3 Sample size for main study from population
#' @param variance Should variance be estimated? (default = TRUE)
#'
#' @return A list of estimated prevalence and variance; optionally,
#'         if variance is FALSE, the truncated prevalence estimate
#' \itemize{
#'  \item pi_hat - Estimated Rogan-Gladen prevalence
#'  \item var_pi_hat - Estimate of the variance of the Rogan-Gladen estimator
#' }
#' @examples
#' ests_rg(0.01, .99, .98, 30, 50, 1000, TRUE)
#'
#' @export

ests_rg <- function(rho_hat, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, variance = TRUE){
  pi_hat <- (rho_hat + sigma_p_hat - 1) / (sigma_e_hat + sigma_p_hat - 1) ###Rogan and Gladen formula
  if (variance == FALSE) {
    return(truncate_01(pi_hat))
  } else {
    # compute variance estimator
    a <- pi_hat^2 * sigma_e_hat * (1 - sigma_e_hat) / n_1
    b <- (1 - pi_hat)^2 * sigma_p_hat * (1 - sigma_p_hat) / n_2
    c <- rho_hat * (1 - rho_hat) / n_3
    var_pi_hat <- (a + b + c) * (sigma_e_hat + sigma_p_hat - 1)^(-2)
    return_list <- list(
      "pi_hat" = pi_hat,
      "var_pi_hat" = var_pi_hat
    )
    return_list
  }
}
