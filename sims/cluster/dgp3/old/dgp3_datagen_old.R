##########
# Generates all sample data for simulation DGP 3.
# Relies on a SLURM-managed cluster.
# Data are stored in n_sim different files corresponding to the 
# n_sim simulations. 

# these two directories must be specified by the user 
libs <- "/nas/longleaf/home/srosin/RLibs/" # directory containing R packages
user_home_dir <- "/nas/longleaf/home/srosin/rgStandardized/" # top-level directory for the project

setwd(user_home_dir)
library(tidyverse)

# load helpful functions and parameter values from source files, 
# all of which are in the top-level user_home_dir
source("estimation_fns.R")
source("sim_fns.R")
source("param_values_var.R")

# sim parameter values
set.seed(2021)
n_sims <- 5 # number of simulations
n_strata <- 40 # number of strata for this dgp
vars_std <- c("z1", "z2", "z3")

# The known stratum proportions (the gamma_{zj}s) must be prespecified,
# and they are loaded here
gammas <- read_csv("dgp3_stratum_props.csv",
                   col_types = cols(
                     z1 = col_character(), 
                     z2 = col_character(), 
                     z3 = col_character(), 
                     stratum_prop = col_double(),
                     sampling_prob = col_double()
                   ))

# Create copies of the gamma (stratum_prop) dataframe,
# with stratum-specific prevalence created from a true logistic model. 
# The intercept of the logistic model varies to vary the marginal prevalence
prevs <- seq(.01, .20, by = .01)
stratum_props <- vector(mode = "list", length = length(prevs)) # create list of stratum proportion dataframes
for(p in 1:length(prevs)){
  s <- gammas %>% dplyr::mutate(
    prev = inv.logit(alpha_0[p]+alpha_1*(gammas$z1=="z11")+
                       alpha_2*(gammas$z2=="z20")+alpha_3*(gammas$z2=="z21")+
                       alpha_4*(gammas$z3=="z30")+alpha_5*(gammas$z3=="z31"))
  )
  stratum_props[[p]] <- s
}

# Uncomment to print prevalences, checking that they are, e.g., {.01, .05, .2}
for(s in 1:length(stratum_props)){
  print(sum(stratum_props[[s]]$stratum_prop * stratum_props[[s]]$prev))
}

# fully factorial combination of sample sizes and parameters, 
# where each row is a sub-scenario
sim_conditions <- tidyr::crossing(
  n_1, n_2, n_3, sigma_e, sigma_p, stratum_props) %>% 
  rowwise() %>% 
  mutate(prev = sum(stratum_props$stratum_prop * stratum_props$prev),
         sigma_e_hat = NA_real_,
         sigma_p_hat = NA_real_,
         rho_hat = NA_real_,
         hat_pi_RG = NA_real_,
         hat_pi_SRG = NA_real_, 
         hat_pi_SRGM = NA_real_
)

for(sim in 1:n_sims){
  print(paste("Simulation number ", sim, sep = ""))
  sample_list <- vector(mode = "list", length = nrow(sim_conditions)) # empty list to store each sample
  
  # simulate the data, iterating through the subscenarios
  for(i in 1:nrow(sim_conditions)){
    row <- sim_conditions[i,]
    dat <- gen_data_dgp3(row$n_1, row$sigma_e, row$n_2, row$sigma_p, 
                              row$n_3, as.data.frame(row$stratum_props), vars_std)
    sim_conditions[i,"sigma_e_hat"] <- dat$sigma_e_hat
    sim_conditions[i,"sigma_p_hat"] <- dat$sigma_p_hat
    sim_conditions[i,"rho_hat"] <- dat$rho_hat
    sample_list[[i]] <- dat$sample # store the entire sample 
  }
  # store this simulation's data in an .RData file. This data format is used because 
  # sim_conditions is a dataframe, itself containing a dataframe as a column (stratum_props)
  output_filename <- paste("dgp3_datasets/sim_",sim,".RData",sep="")
  save(sim_conditions,sample_list, file = output_filename)
}
