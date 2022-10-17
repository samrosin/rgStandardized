library(tidyverse)
library(here)
library(readxl)

path <- here("strataprops_raw.xlsx")
age_cat_breaks <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, Inf)
  
# create male data --------------------------------------------------------

# extract male data from each Excel sheet
dat_m <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = path, range = cell_rows(120:231))
dat_m$Overview <- NULL # remove "Overview" sheet

# extract a single df for males with columns for province, age group, and count (of persons)
dat_m <- dat_m %>%
  map(select, ...1, `2020`) %>% 
  map(rename, age_cat = ...1, count = `2020`) %>% 
  map(dplyr::mutate, age_cat = as.numeric(str_replace_all(age_cat, "[^0-9.-]", ""))) %>% # remove characters ("ans") from age_cat column
  map_df(~ as.data.frame(.x), .id = "province") %>%
  mutate(age_cat = cut(age_cat, age_cat_breaks, right = FALSE)) %>%
  mutate(sex = "m") %>%
  group_by(province, age_cat, sex) %>% 
  summarise(count = sum(count)) %>%
  as_tibble()

# create female data ------------------------------------------------------
dat_f <- path %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = path, range = cell_rows(236:347))
dat_f$Overview <- NULL # remove "Overview" sheet

# extract a single df for males with columns for province, age_cat group, and count (of persons)
dat_f <- dat_f %>%
  map(select, ...1, `2020`) %>% 
  map(rename, age_cat = ...1, count = `2020`) %>% 
  map(dplyr::mutate, age_cat = as.numeric(str_replace_all(age_cat, "[^0-9.-]", ""))) %>% # remove characters ("ans") from age_cat column
  map_df(~ as.data.frame(.x), .id = "province") %>%
  mutate(age_cat = cut(age_cat, age_cat_breaks, right = FALSE)) %>%
  mutate(sex = "f") %>%
  group_by(province, age_cat, sex) %>% 
  summarise(count = sum(count)) %>%
  as_tibble()

# combine data, rename provinces, and write to csv -------------------------------------------

dat <- bind_rows(dat_m, dat_f) 
dat$province <- dat$province %>%
  iconv(from = "UTF-8", to='ASCII//TRANSLIT') %>%
  str_replace_all("[[`']]", "") %>% # remove special characters (accent marks)
  recode("Anvers" = "Antwerp",
         "Brabant flamand" = "Flemish Brabant",
         "Brabant wallon" = "Walloon Brabant",
         "Flandre occidentale" = "West Flanders",
         "Flandre orientale" = "East Flanders",
           "Limbourg" = "Limburg",
         "Region de Bruxelles-Capitale" = "Brussels")

# add stratum_prop column
dat$stratum_prop <- dat$count / sum(dat$count)

write_csv(dat, here("strataprops_clean.csv"))
