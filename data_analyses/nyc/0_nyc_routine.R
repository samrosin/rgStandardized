# analysis of nyc seroprevalence study

# setup and read data -----------------------------------------------------
library(tidyverse)
library(here)
library(binom)
library(summarytools)
source(here("estimation_fns.R"))

age_grp_cuts <- c(0, 20, 40, 60, 80, Inf) #how to bucket/discretize age 

methods <- c("naive", "rg", "srg", "srgm")
srgm_formula <- formula("x ~ sex + age + race + sex*age")


# load data
strata_props <- read_csv(here("strata_props_clean.csv")) 
strata_props$race <- recode(strata_props$race, "White or Caucasian" = "White") #recode race 
# strata_props$age <- recode(strata_props$age,
#                            "20-29" = "[20,30)", "30-39" = "[30,40)",
#                            "40-49" = "[40,50)", "50-59" = "[50,60)",
#                            "60-69" = "[60,70)", "70-79" = "[70,80)",
#                            "80-plus" = "[80,Inf)")

sero <- read_csv(here("2022-09-16 Serosurvey demographic data and titers until 07-05-20.csv")) %>%
  rename(x = spike_positive)



# format demographics
sero$age <- cut(sero$age, breaks = age_grp_cuts, right = FALSE) ### should match age groups 
sero <- sero %>% filter(sex %in% c("Male", "Female")) # filter out one individual with indeterminate sex
sero <- sero %>% filter(race != "Unknown") # remove unknown race

sero$race <- recode(sero$race, White="White", 
                    Ap = "Asian", 
                    B = "Black or African American",
       Ba = "Black or African American", "African-American" = "Black or African American", 
       B2 = "Black or African American", B3 = "Black or African American",
       B4 = "Black or African American", B7 = "Black or African American",
       B9 = "Black or African American", Ba = "Black or African American",
       Bb = "Black or African American", Be = "Black or African American",
       Bg = "Black or African American", Bj = "Black or African American",
       Bk = "Black or African American", Bm = "Black or African American",
       Bn = "Black or African American", Bo = "Black or African American",
       Bp = "Black or African American", Br = "Black or African American",
       Bs = "Black or African American", Bt = "Black or African American",
       Bu = "Black or African American", 
       "Asian Indian" = "Asian", Bangladeshi = "Asian", Bhutanese = "Asian", 
       Burmese = "Asian", Cambodian = "Asian", Chinese = "Asian",
       "Dominica Islander" = "Black or African American",
       Eritrean = "Black or African American", Filipino = "Asian", 
       Ghanaian = "Black or African American", Guinean = "Black or African American",
       Haitian = "Black or African American",
       Indonesian = "Asian", Japanese = "Asian", Jamaican = "Black or African American",
       Kenyan = "Black or African American", Korean = "Asian",
       Laotian = "Asian", Malaysian = "Asian", Malian = "Black or African American",
       Nepalese = "Asian", Nigerian = "Black or African American", 
       Okinawan = "Asian", Pakistani = "Asian", 
       "Other: North African" = "Black or African American",
       "Other: South African" = "Black or African American",
       "Other: West African" = "Black or African American",
       Senegalese = "Black or African American", "Sri lankan" = "Asian",
       "Sri Lankan" = "Asian", Taiwanese = "Asian", Thai = "Asian", 
       Trinidadian = "Black or African American", Ugandan = "Black or African American", 
       Vietnamese = "Asian", "West Indian" = "Black or African American",
       .default="Other")

# filter to routine care (vs. urgent care)
routine <- sero %>% filter(group_code == 0)
routine <- routine %>% select(-c(final_code, group_code, spike_titer))

# recode week into weeks grouping that roughly corresponds to month
routine$weeks <- recode(routine$week,
          "6" = "6-10", "7" = "6-10", "8" = "6-10", "9" = "6-10", "10" = "6-10", 
          "11" = "11-14", "12" = "11-14", "13" = "11-14", "14" = "11-14", 
          "15" = "15-18", "16" = "15-18", "17" = "15-18", "18" = "15-18", 
          "19" = "19-22", "20" = "19-22", "21" = "19-22", "22" = "19-22", 
          "23" = "23-27", "24" = "23-27", "25" = "23-27", "26"= "23-27",
          "27" = "23-27",
          .default = "NA")
routine$weeks <- factor(routine$weeks, 
                        levels = c("6-10", "11-14", 
                        "15-18", "19-22", "23-27"))

# setup results dataframe 1/2
month_vec <- rep(NA, length(unique(routine$weeks)) * length(methods))
i <- 1
for(month in sort(unique(routine$weeks))){
  month_vec[i:(i + length(methods) - 1)] <- rep(month, length(methods))
  i <- i + length(methods)
}

# setup results dataframe 2/2
nrow_res <- length(unique(routine$weeks)) * length(methods)
method <- rep(methods, length(unique(routine$weeks)))
mean <- rep(NA, nrow_res)
lower <- rep(NA, nrow_res)
upper <- rep(NA, nrow_res)
res <- data.frame(weeks = month_vec, method, mean, lower, upper)
strats <- data.frame(
  weeks = sort(unique(routine$weeks)), strata_missing = NA
)

# sensitivity data -- from Nature paper Extended Data Table 1
n_1 <- 40 # number of positives for SARS-CoV-2 antibodies
sum_x_delta1 <- 38 #number of positives that tested positive 
sigma_e_hat <- sum_x_delta1 / n_1 # sensitivity estimate

# sensitivity data -- from same source
n_2 <- 74 # number of people confirmed negative for antibodies 
sum_x_delta2 <- 0 # number of negatives that tested positive
sigma_p_hat <- 1 - (sum_x_delta2 / n_2) #sensitivity estimate

# main study data -- routine
# n_3 <- nrow(routine) #number of tests
# sum_x_delta3 <- sum(routine$x) # number positive
# rho_hat <- sum_x_delta3 / n_3
# binom.confint(sum_x_delta3, n_3) %>%
#   filter(method == "exact")
#i <- 1; n_3vec <- rep(NA, 5)
for(month in sort(unique(routine$weeks))){

  print(month)
  # filter data to correct weeks
  routine_weeks <- routine %>% filter(weeks == month)
  n_3 <- nrow(routine_weeks)
  #n_3vec[i] <- n_3
  i <- i + 1
  sum_x_delta3 <- sum(routine_weeks$x)
  rho_hat <- sum_x_delta3 / n_3
  
 # print(head(routine_weeks))
  # naive results
  res[(res$weeks == month) & (res$method == "naive"),]$mean <- rho_hat
  res[(res$weeks == month) & (res$method == "naive"),]$lower <-
    binom.confint(sum_x_delta3, n_3) %>%
    filter(method == "exact") %>%
    select(lower) %>%
    as.numeric()
  res[(res$weeks == month) & (res$method == "naive"),]$upper <-
    binom.confint(sum_x_delta3, n_3) %>%
    filter(method == "exact") %>%
    select(upper) %>%
    as.numeric()
  
  # pi_RG -------------------------------------------------------------------
  rg_ests <- ests_rg(rho_hat, sigma_e_hat, sigma_p_hat, 
                     n_1, n_2, n_3, variance = TRUE)
  est_rg <- rg_ests[1]
  ASE_rg <- sqrt(rg_ests[2])
  ci_lower_rg <- max(0, est_rg - qnorm(.975) * ASE_rg)
  ci_upper_rg <- min(1, est_rg + qnorm(.975) * ASE_rg)
  
  # add to results df
  res[(res$weeks == month) & (res$method == "rg"),]$mean <- est_rg
  res[(res$weeks == month) & (res$method == "rg"),]$lower <- ci_lower_rg
  res[(res$weeks == month) & (res$method == "rg"),]$upper <- ci_upper_rg
  
# # pi_SRG ------------------------------------------------------------------
  vars_std <- c("age", "race", "sex")
  
  # use restriction to create new target population dataset
  strata_samp <- unique(routine_weeks[, c(vars_std)]) # nrow = 213
  strata_tgt <- unique(strata_props[, c(vars_std)]) # nrow = 220
  # difference between the above two are the strata not in the sample; nrow = 3
  strata_missing_samp <- setdiff(strata_tgt, strata_samp) # nrow = 7
  
  # store strata data
  strats[strats$weeks == month,]$strata_missing <- nrow(strata_missing_samp)
  
  # create a new target dataset that only includes the restricted population, 
  # i.e. excludes the strata missing in the sample
  tgt_restricted <- anti_join(strata_props, strata_missing_samp, by = vars_std)
  tgt_restricted$stratum_prop <- tgt_restricted$pat_count / sum(tgt_restricted$pat_count)
  
  # include column for stratum proportion/gamma in the sample dataset
  routine_weeks_srg <- left_join(routine_weeks, tgt_restricted, by = vars_std) %>% 
    select(-c(pat_count)) 
  
  std_ests <- ests_std(routine_weeks_srg, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, 
                       vars_std, variance = TRUE)
  
  est_srg <- std_ests[1]
  ASE_srg <- sqrt(std_ests[2])
  ci_lower_srg <- max(0, est_srg - qnorm(.975) * ASE_srg)
  ci_upper_srg <- min(1, est_srg + qnorm(.975) * ASE_srg)
  
  # add to results df
  res[(res$weeks == month) & (res$method == "srg"),]$mean <- est_srg
  res[(res$weeks == month) & (res$method == "srg"),]$lower <- ci_lower_srg
  res[(res$weeks == month) & (res$method == "srg"),]$upper <- ci_upper_srg
  
  ests_srgm <- ests_std_model(routine_weeks, strata_props, sigma_e_hat,
                               sigma_p_hat, n_1, n_2, n_3, vars_std,
                               mod_formula = srgm_formula,
                               variance = TRUE)
  
  est_srgm <- ests_srgm[1]
  ASE_srgm <- sqrt(ests_srgm[2])
  ci_lower_srgm <- max(0, est_srgm - qnorm(.975) * ASE_srgm)
  ci_upper_srgm <- min(1, est_srgm + qnorm(.975) * ASE_srgm)

  # add to results df
  res[(res$weeks == month) & (res$method == "srgm"),]$mean <- est_srgm
  res[(res$weeks == month) & (res$method == "srgm"),]$lower <- ci_lower_srgm
  res[(res$weeks == month) & (res$method == "srgm"),]$upper <- ci_upper_srgm
  

}


# plot results ------------------------------------------------------------

# helpful: https://datascienceplus.com/lattice-like-forest-plot-using-ggplot2-in-r/

# # collection round dates
month <- factor(rep(c(1:5), each = 4),
             levels = c(1:5))
res$month <- month

res_routine <- res


