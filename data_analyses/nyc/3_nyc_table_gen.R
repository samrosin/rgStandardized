# analysis of nyc seroprevalence study

# setup and read data -----------------------------------------------------
library(tidyverse)
library(here)
library(binom)
source(here("estimation_fns.R"))

age_grp_cuts <- c(0, 20, 40, 60, 80, Inf) #how to bucket/discretize age 

# load data
strata_props <- read_csv(here("strata_props_clean.csv")) 
strata_props$race <- recode(strata_props$race, "White or Caucasian" = "White") #recode race 

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

urgent <- sero %>% filter(group_code == 1)
urgent <- urgent %>% select(-c(final_code, group_code, spike_titer))

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


# routine care data -------------------------------------------------------

routine_sex <- routine %>% 
  group_by(sex) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

routine_age <- routine %>% 
  group_by(age) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

routine_race <- routine %>% 
  group_by(race) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

# urgent care data --------------------------------------------------------

urgent_sex <- urgent %>% 
  group_by(sex) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

urgent_age <- urgent %>% 
  group_by(age) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

urgent_race <- urgent %>% 
  group_by(race) %>%
  summarise(pat_count = n()) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

# NYC population data -----------------------------------------------------

strata_props_sex <- strata_props %>% 
  group_by(sex) %>%
  summarise(pat_count = sum(pat_count)) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

strata_props_age <- strata_props %>% 
  group_by(age) %>%
  summarise(pat_count = sum(pat_count)) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

strata_props_race <- strata_props %>% 
  group_by(race) %>%
  summarise(pat_count = sum(pat_count)) %>%
  mutate(strata_prop = pat_count / sum(pat_count))

