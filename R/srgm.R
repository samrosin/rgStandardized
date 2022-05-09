#' Compute logistic regression standardized Rogan-Gladen estimates
#'
#' @description
#' Logistic regression outcome model for the mismeasured
#' stratum-specific prevalences \eqn{\rho_j}.
#' Corresponds to \eqn{\hat{\pi_{SRGM}}} and \eqn{\hat V_{\pi, SRGM}} in the manuscript.
#'
#' @param sample A dataframe with data from the sample, including the stratum_props (\eqn{\gamma_j} values) which will be used in standardization
#' @param stratum_props the stratum_props dataframe. note that the vars_std need to be the first columns from the left
#' @param sigma_e_hat Estimated sensitivity
#' @param sigma_p_hat Estimated specificity
#' @param n_1 Sample size for sensitivity validation dataset
#' @param n_2 Sample size for specificity validation dataset
#' @param n_3 Sample size for main study from population
#' @param vars_std The names of variables to standardize over, corresponding to columns in sample
#' @param mod_formula An object of class "formula" that symbolically describes the model to be fitted
#' @param variance Should variance be estimated? (default TRUE)
#'
#' @return A list of estimated prevalence and variance; optionally,
#'         if variance is FALSE, the truncated prevalence estimate
#' \itemize{
#'  \item hat_pi_srgm - Estimated model-based standardized Rogan-Gladen prevalence
#'  \item hat_var_pi_srgm - Estimate of the variance of the model-based standardized Rogan-Gladen estimator
#' }
#'
#' @importFrom stats binomial glm model.frame model.matrix predict
#'
#' @examples
#'
#'  serology <- belgium$serology
#'  belg_cr <- serology[serology$collection_round == 1,]
#'  belg_cr$x <- ifelse(belg_cr$igg_cat=="positive", 1, 0)
#'  belg_cr_srgm <- as.data.frame(belg_cr)
#'
#'  srgm(belg_cr_srgm, belgium$strataprops, 154/181, 322/326,
#'       181, 326, nrow(belg_cr_srgm), vars_std = c("age_cat", "sex"),
#'       mod_formula = formula("x ~ sex + age_cat + sex*age_cat"), variance = TRUE)
#'
#' @export

srgm <- function(sample, stratum_props, sigma_e_hat,
                           sigma_p_hat, n_1, n_2, n_3, vars_std,
                           mod_formula, variance = FALSE){

  #fit the logistic regression
  x_model <- glm(formula = mod_formula, family = binomial, data = sample) #regression model

  # make predictions for all strata
  stratum_props$model_pred <- predict(x_model, newdata = stratum_props, type = "response")
  stratum_props$std_est_model <- stratum_props$model_pred * stratum_props$stratum_prop

  # point estimate
  pi_hat_mst <- (sum(stratum_props$std_est_model) - (1 - sigma_p_hat)) /
    (sigma_e_hat - (1 - sigma_p_hat))

  if(variance == FALSE){
    return(truncate_01(pi_hat_mst)) # truncated estimate
  } else{

    n <- n_1 + n_2 + n_3 #total sample size

    #####
    # estimate the different parts of the variance estimator
    coef <- as.matrix(coef(x_model)) # regression coefficients
    mat <- model.matrix(x_model) # design matrix
    a <- matrix(c(n_1 / n, 0, 0, n_2 / n), nrow = 2, ncol = 2) # estimate A

    # estimate B
    b <- matrix(NA, nrow = nrow(coef), ncol = nrow(coef)) # declare B-hat as an empty p * p matrix
    exp_coef_times_mat <- exp(crossprod(coef, t(mat))) # this is a 1 by n matrix with ith element exp(hat beta ^T %*% h(Z_i))
    for(j in 1:nrow(coef)){
      for(k in 1:nrow(coef)){
        unit_contrib_b <- rep(NA, nrow(mat)) # each unit's contribution to this jkth element of B
        for(i in 1:ncol(exp_coef_times_mat)){
          unit_contrib_b[i] <- mat[i,j] * mat[i,k] * exp_coef_times_mat[1,i] /
            (1 + exp_coef_times_mat[1,i])^2
        }
        b[j,k] <- (n_3 / n) * mean(unit_contrib_b)
      }
    }

    c <- matrix(c(0, pi_hat_mst, 0, -1 + pi_hat_mst), nrow = 2, ncol = 2) # estimate C

    # estimate D
    ff <- mod_formula[-2]
    m <- model.frame(ff, stratum_props)
    mat_strat <- model.matrix(ff, m) # model matrix for all strata, even those not in sample support

    d <- matrix(0, nrow = 2, ncol = nrow(coef)) # declare D as 2 * p matrix of 0s (convenient since the 2nd row is 0)
    for(l in 1:nrow(coef)){
      stratum_contrib <- rep(NA, nrow(mat_strat)) # contribution of the jth stratum to the sum
      exp_coef_times_mat_strat <- exp(crossprod(coef, t(mat_strat))) # this is a row vector with ith element exp(hat beta ^T %*% h(Z_i))
      for(j in 1:nrow(mat_strat)){
        stratum_contrib[j] <- as.matrix(mat_strat[j,l]) * exp_coef_times_mat_strat[1,j] /
          (1 + exp_coef_times_mat_strat[1,j])^2 * stratum_props$stratum_prop[j] # compare to formula in supporting info
      }
      d[1,l] <- -1 * sum(stratum_contrib)
    }

    e <- matrix(c(1, -1, 0, sigma_e_hat + sigma_p_hat -1), nrow = 2, ncol = 2) # estimate E
    f <- matrix(c(n_1 * sigma_e_hat * (1 - sigma_e_hat) / n, 0,
                  0, n_2 * sigma_p_hat * (1 - sigma_p_hat) / n ), nrow = 2, ncol = 2) # estimate F

    # estimate G
    g <- matrix(NA, nrow = nrow(coef), ncol = nrow(coef)) #declare G as an empty p * p matrix
    invlogit_coef_times_mat <- inv.logit(crossprod(coef, t(mat))) # 1 times n row vector
    for(j in 1:nrow(coef)){
      for(k in 1:nrow(coef)){
        unit_contrib_g <- rep(NA, nrow(mat)) # each unit's contribution to this jkth element of B
        for(i in 1:ncol(invlogit_coef_times_mat)){
          unit_contrib_g[i] <- mat[i,j] * mat[i,k] *
            (x_model$y[i] - invlogit_coef_times_mat[1,i])^2 # note that x_model$y is the model outcome, which is in our model denoted x
        }
        g[j,k] <- (n_3 / n) * mean(unit_contrib_g)
      }
    }

    ####
    # define a_mat, b_mat - the 'bread' and 'meat' matrices
    # of the sandwich variance estimator
    a_mat <- cbind(
      rbind(a, matrix(0, nrow = nrow(b), ncol = ncol(a)), c),
      rbind(matrix(0, nrow = 2, ncol = ncol(b)), b, d),
      rbind(matrix(0, nrow = 2, ncol = 2),
            matrix(0, nrow = nrow(b), ncol = 2), e)
    )

    b_mat <- cbind(
      rbind(f, matrix(0, nrow = nrow(g), ncol = 2),
            matrix(0, nrow = 2, ncol = 2)),
      rbind(matrix(0, nrow = 2, ncol = ncol(g)), g,
            matrix(0, nrow = 2, ncol = ncol(g))),
      rbind(matrix(0, nrow = 2, ncol = 2),
            matrix(0, nrow = nrow(g), ncol = 2),
            matrix(0, nrow = 2, ncol = 2))
    )

    var_hat_theta <- solve(a_mat) %*% b_mat %*% t(solve(a_mat)) # var-covar matrix
    var_hat_pi_mst <- var_hat_theta[nrow(var_hat_theta),ncol(var_hat_theta)] / n # the lower-right element divided by n is the variance estimator of interest

    # return point and variance estimates
    # note point estimate is not truncated for use in CI construction
    return_list <- list(
      "hat_pi_srgm" = pi_hat_mst,
      "hat_var_pi_srgm" = var_hat_pi_mst
    )
    return_list
  }
}
