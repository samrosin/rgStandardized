# analysis of belgian seroprevalence study

# setup and read data -----------------------------------------------------
library(tidyverse)
library(here)
library(binom)
source(here("estimation_fns.R"))

vars_std <- c("age_cat", "province", "sex")
methods <- c("raw", "herzog", "rg", "srg", "srgm")
srgm_formula <- formula("x ~ sex + age_cat + province + sex*age_cat")

# read data
strata_props <- read_csv(here("strataprops_clean.csv"))
belg <- read_csv(here("serology.csv"), guess_max = 10000) %>% 
  select(age_cat, sex, province, igg_cat, collection_round) %>% 
  mutate(igg_cat = fct_collapse(igg_cat, negative = c("negative", "borderline", "LoD")))  # per manuscript, set equivocal values as negative %>% 
 
belg$province <- iconv(belg$province, from = "UTF-8", to='ASCII//TRANSLIT') %>%
  str_replace_all("[[`']]", "") # remove special characters (accent marks)  

collection_rounds <- unique(belg$collection_round) 

cr_vec <- rep(NA, length(collection_rounds) * length(methods))
i <- 1
for(cr in collection_rounds){
  cr_vec[i:(i + length(methods) - 1)] <- rep(cr, length(methods))
  i <- i + length(methods)
}

nrow_res <- length(collection_rounds) * length(methods)
method <- rep(methods, length(collection_rounds))
mean <- rep(NA, nrow_res)
lower <- rep(NA, nrow_res)
upper <- rep(NA, nrow_res)
res <- data.frame(cr = cr_vec, method, mean, lower, upper)

# ests from herzog et al. manuscript using bayesian method, from Table 2
res[res$method == "herzog",]$mean <- c(.018, .053, .062, .047, .037, .033, .042) 
res[res$method == "herzog",]$lower <- c(.010, .042, .051, .037, .026, .023, .031)
res[res$method == "herzog",]$upper <- c(.026, .064, .073, .059, .047, .043, .052)

strats <- data.frame(
  cr = 1:7, strata_missing = NA
) # note 220 strata = 10 * 11 * 2 (age * province * sex)

for(cr in collection_rounds){
  belg_cr <- belg %>% filter(collection_round == cr)
  
  n_3 <- nrow(belg_cr)
  sum_x_delta3 <- nrow(belg_cr[belg_cr$igg_cat == "positive",])
  rho_hat <- sum_x_delta3 / n_3
  
  # validation data - using >= 1.1 EI cutoff for positive, see Table S1.1 
  n_1 <- 181
  sigma_e_hat <- 154/181
  n_2 <- 326
  sigma_p_hat <- 322/326
  
  # add to results df 
  res[(res$cr == cr) & (res$method == "raw"),]$mean <- rho_hat
  res[(res$cr == cr) & (res$method == "raw"),]$lower <-
    binom.confint(sum_x_delta3, n_3) %>%
    filter(method == "exact") %>%
    select(lower) %>%
    as.numeric()
  res[(res$cr == cr) & (res$method == "raw"),]$upper <-
      binom.confint(sum_x_delta3, n_3) %>%
      filter(method == "exact") %>%
      select(upper) %>%
      as.numeric()

  # rogan-gladen ------------------------------------------------------------
  
  rg_ests <- ests_rg(rho_hat, sigma_e_hat, sigma_p_hat, 
                     n_1, n_2, n_3, variance = TRUE)
  est_rg <- rg_ests[1]
  ASE_rg <- sqrt(rg_ests[2])
  ci_lower_rg <- max(0, est_rg - qnorm(.975) * ASE_rg)
  ci_upper_rg <- min(1, est_rg + qnorm(.975) * ASE_rg)
  
  #print(paste("Raw estimate: ", round(sum(belg_cr$igg_cat=="positive") / nrow(belg_cr), 4)))
  #print(paste("Rogan-Gladen estimate: ", round(est_rg, 4), " CI: ", round(ci_lower_rg, 4), round(ci_upper_rg, 4)))
  
  # add to results df
  res[(res$cr == cr) & (res$method == "rg"),]$mean <- est_rg
  res[(res$cr == cr) & (res$method == "rg"),]$lower <- ci_lower_rg
  res[(res$cr == cr) & (res$method == "rg"),]$upper <- ci_upper_rg
 
  # SRG ---------------------------------------------------------------------
  
  # use restriction to create new target population dataset
  strata_samp <- unique(belg_cr[, c(vars_std)]) # nrow = 213
  strata_tgt <- unique(strata_props[, c(vars_std)]) # nrow = 220
  # difference between the above two are the strata not in the sample; nrow = 3
  strata_missing_samp <- setdiff(strata_tgt, strata_samp) # nrow = 7
  
  # store strata data
  strats[strats$cr == cr,]$strata_missing <- nrow(strata_missing_samp)
  
  # create a new target dataset that only includes the restricted population, 
  # i.e. excludes the strata missing in the sample
  tgt_restricted <- anti_join(strata_props, strata_missing_samp, by = vars_std)
  tgt_restricted$stratum_prop <- tgt_restricted$count / sum(tgt_restricted$count)
  
  # include column for stratum proportion/gamma in the sample dataset
  belg_cr_srg <- left_join(belg_cr, tgt_restricted, by = vars_std) %>% 
    select(-c(count)) %>% 
    mutate(igg_cat = ifelse(igg_cat=="positive", 1, 0)) %>% 
    rename(x = igg_cat)
  
  std_ests <- ests_std(belg_cr_srg, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, 
                       vars_std, variance = TRUE)
  
  est_srg <- std_ests[1]
  ASE_srg <- sqrt(std_ests[2])
  ci_lower_srg <- max(0, est_srg - qnorm(.975) * ASE_srg)
  ci_upper_srg <- min(1, est_srg + qnorm(.975) * ASE_srg)
  
  #print(paste("Raw estimate: ", round(sum(belg_cr$igg_cat=="positive") / nrow(belg_cr), 4)))
  #print(paste("SRG estimate: ", round(est_srg, 4), " CI: ", round(ci_lower_srg, 4), round(ci_upper_srg, 4)))
  
  # add to results df
  res[(res$cr == cr) & (res$method == "srg"),]$mean <- est_srg
  res[(res$cr == cr) & (res$method == "srg"),]$lower <- ci_lower_srg
  res[(res$cr == cr) & (res$method == "srg"),]$upper <- ci_upper_srg
  
  # SRGM --------------------------------------------------------------------

  belg_cr_srgm <- belg_cr %>%
    mutate(igg_cat = ifelse(igg_cat=="positive", 1, 0)) %>%
    rename(x = igg_cat) %>%
    as.data.frame()

  model_ests <- ests_std_model(belg_cr_srgm, strata_props, sigma_e_hat,
                               sigma_p_hat, n_1, n_2, n_3, vars_std,
                               mod_formula = srgm_formula,
                               variance = TRUE)

  est_srgm <- model_ests[1]
  ASE_srgm <- sqrt(model_ests[2])
  ci_lower_srgm <- max(0, est_srgm - qnorm(.975) * ASE_srgm)
  ci_upper_srgm <- min(1, est_srgm + qnorm(.975) * ASE_srgm)

  #print(paste("Raw estimate: ", round(sum(belg_cr$igg_cat=="positive") / nrow(belg_cr), 4)))
  #print(paste("SRGM estimate: ", round(est_srgm, 4), " CI: ", round(ci_lower_srgm, 4), round(ci_upper_srgm, 4)))

  # add to results df
  res[(res$cr == cr) & (res$method == "srgm"),]$mean <- est_srgm
  res[(res$cr == cr) & (res$method == "srgm"),]$lower <- ci_lower_srgm
  res[(res$cr == cr) & (res$method == "srgm"),]$upper <- ci_upper_srgm

}

write_csv(res, here("ests.csv"))

# selection_bias_plot -----------------------------------------------------

# View(belg)
# View(strata_props)
# summary(strata_props)
# dfSummary(belg)
# dfSummary(strata_props)

belg_selectionbias_plot <- ggplot(data = sp3_05, aes(x=stratum_prop, y = sampling_prob)) + 
  geom_jitter(aes(size=prev_x), alpha = .6,
              color = "black", width = 0.001,
              show.legend = FALSE)  +
  geom_abline(colour = "grey50", size = 2) + 
  xlab("Stratum proportion in target pop.") + ylab("Sampling probability") + 
  #xlab(expression(gamma[j])) + 
  #ylab(expression(s[j])) + 
  theme_bw() + 
  theme(text = element_text(size=20),
        legend.position = c(.17,0.79),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16))# + 
scenario3_selectionbias_plot
