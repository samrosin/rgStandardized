# analysis of nyc seroprevalence study

# setup and read data -----------------------------------------------------
library(tidyverse)
library(here)
library(binom)
library(summarytools)
source(here("estimation_fns.R"))

#age_grp_cuts <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, Inf) #how to bucket/discretize age 
age_grp_cuts <- c(0, 20, 40, 60, 80, Inf) #how to bucket/discretize age 

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
                    Eritrean = "Black or African American", Filipino = "Asian", 
                    Ghanaian = "Black or African American", Guinean = "Black or African American",
                    Indonesian = "Asian", Japanese = "Asian",
                    Kenyan = "Black or African American", Korean = "Asian",
                    Laotian = "Asian", Malaysian = "Asian", Malian = "Black or African American",
                    Nepalese = "Asian", Nigerian = "Black or African American", 
                    Okinawan = "Asian", Pakistani = "Asian", 
                    Senegalese = "Black or African American", "Sri lankan" = "Asian",
                    "Sri Lankan" = "Asian", Taiwanese = "Asian", Thai = "Asian", 
                    Ugandan = "Black or African American", 
                    Vietnamese = "Asian",
                    .default="Other")

# filter to urgent care (vs. routine care)
urgent <- sero %>% filter(group_code == 1)
urgent <- urgent %>% select(-c(week, final_code, group_code, spike_titer))

# sensitivity data -- from Nature paper Extended Data Table 1
n_1 <- 40 # number of positives for SARS-CoV-2 antibodies
sum_x_delta1 <- 38 #number of positives that tested positive 
sigma_e_hat <- sum_x_delta1 / n_1 # sensitivity estimate

# sensitivity data -- from same source
n_2 <- 74 # number of people confirmed negative for antibodies 
sum_x_delta2 <- 0 # number of negatives that tested positive
sigma_p_hat <- 1 - (sum_x_delta2 / n_2) #sensitivity estimate

# main study data -- urgent
n_3 <- nrow(urgent) #number of tests
sum_x_delta3 <- sum(urgent$x) # number positive
rho_hat <- sum_x_delta3 / n_3

binom.confint(sum_x_delta3, n_3) %>%
  filter(method == "exact")

# pi_RG -------------------------------------------------------------------

# rg_ests_urgent is a vector containing the Rogan-Gladen point and variance estimates 
rg_ests <- ests_rg(rho_hat, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, 
                   variance = TRUE) 
est_rg <- rg_ests[1]
ASE_rg <- sqrt(rg_ests[2])
ci_lower_rg <- max(0, est_rg - qnorm(.975) * ASE_rg)
ci_upper_rg <- min(1, est_rg + qnorm(.975) * ASE_rg)


# pi_SRG ------------------------------------------------------------------
vars_std <- c("age", "race", "sex")

# FOR NOW, stratify to age 20+
#urgent <- urgent %>% filter(!is.na(age))
urgent_srg <- left_join(urgent, strata_props, 
                         by = vars_std) %>%
  select(-pat_count)

# use restriction *if necessary* to create new target population dataset
 strata_samp <- unique(urgent[, c(vars_std)]) # nrow = 213
strata_tgt <- unique(strata_props[, c(vars_std)]) # nrow = 220
# # difference between the above two are the strata not in the sample; nrow = 3
 strata_missing_samp <- setdiff(strata_tgt, strata_samp) # nrow = 7

srg_ests <- ests_std(urgent_srg, sigma_e_hat, sigma_p_hat, n_1, n_2, n_3, 
                     vars_std, variance = TRUE)
est_srg <- srg_ests[1]
ASE_srg <- sqrt(srg_ests[2])
ci_lower_srg <- max(0, est_srg - qnorm(.975) * ASE_srg)
ci_upper_srg <- min(1, est_srg + qnorm(.975) * ASE_srg)

srgm_formula <- formula("x ~ sex + age + race + age*race + sex*age + sex*race")
srgm_ests <- ests_std_model(urgent_srg, strata_props, sigma_e_hat,
                            sigma_p_hat, n_1, n_2, n_3, vars_std,
                            mod_formula = srgm_formula,
                            variance = TRUE)

est_srgm <- srgm_ests[1]
ASE_srgm <- sqrt(srgm_ests[2])
ci_lower_srgm <- max(0, est_srgm - qnorm(.975) * ASE_srgm)
ci_upper_srgm <- min(1, est_srgm + qnorm(.975) * ASE_srgm)

