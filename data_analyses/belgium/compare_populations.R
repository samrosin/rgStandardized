library(tidyverse)
library(here)

# read data
strata_props <- read_csv(here("strataprops_clean.csv"))
belg <- read_csv(here("serology.csv"), guess_max = 10000) %>% 
  select(age_cat, sex, province, igg_cat, collection_round) %>% 
  mutate(igg_cat = fct_collapse(igg_cat, negative = c("negative", "borderline", "LoD")))  # per manuscript, set equivocal values as negative %>% 

belg$province <- iconv(belg$province, from = "UTF-8", to='ASCII//TRANSLIT') %>%
  str_replace_all("[[`']]", "") # remove special characters (accent marks)  

#dfSummary(belg) # look at data and check that it matches their Table 1 - from summaryTools library

belg_1 <- belg %>% filter(collection_round == 1)
belg_7 <- belg %>% filter(collection_round == 7)


# collection round 1 ------------------------------------------------------

nrow(belg_1) 

# sex
unique(belg_1$sex)
print(nrow(belg_1[belg_1$sex == "f",]))
nrow(belg_1[belg_1$sex == "f",]) / nrow(belg_1) # 54.0%
print(nrow(belg_1[belg_1$sex == "m",]))
nrow(belg_1[belg_1$sex == "m",]) / nrow(belg_1) # 46.0%

# age group
unique(belg_1$age_cat)
for(cat in unique(belg_1$age_cat)){
  print(cat)
  print(nrow(belg_1[belg_1$age_cat == cat,]))
  print(round(nrow(belg_1[belg_1$age_cat == cat,]) / nrow(belg_1), 4))
}

# province
unique(belg_1$province)
for(prov in sort(unique(belg_1$province))){
  print(prov)
  print(nrow(belg_1[belg_1$province == prov,]))
  print(round(nrow(belg_1[belg_1$province == prov,]) / nrow(belg_1), 4))
}


# collection round 7 ------------------------------------------------------


nrow(belg_7) 

# sex
unique(belg_7$sex)
print(nrow(belg_7[belg_7$sex == "f",]))
nrow(belg_7[belg_7$sex == "f",]) / nrow(belg_7) # 54.0%
print(nrow(belg_7[belg_7$sex == "m",]))
nrow(belg_7[belg_7$sex == "m",]) / nrow(belg_7) # 46.0%

# age group
unique(belg_7$age_cat)
for(cat in unique(belg_7$age_cat)){
  print(cat)
  print(nrow(belg_7[belg_7$age_cat == cat,]))
  print(round(nrow(belg_7[belg_7$age_cat == cat,]) / nrow(belg_7), 4))
}

# province
for(prov in sort(unique(belg_7$province))){
  print(prov)
  print(nrow(belg_7[belg_7$province == prov,]))
  print(round(nrow(belg_7[belg_7$province == prov,]) / nrow(belg_7), 4))
}

# stratum proportions -----------------------------------------------------
strata_props <- strata_props %>% select(-stratum_prop)

strata_long <- strata_props %>% 
  uncount(count)

nrow(strata_long) 

# sex
unique(strata_long$sex)
print(nrow(strata_long[strata_long$sex == "f",]))
nrow(strata_long[strata_long$sex == "f",]) / nrow(strata_long) # 54.0%
print(nrow(strata_long[strata_long$sex == "m",]))
nrow(strata_long[strata_long$sex == "m",]) / nrow(strata_long) # 46.0%

# age group
unique(strata_long$age_cat)
for(cat in unique(strata_long$age_cat)){
  print(cat)
  print(nrow(strata_long[strata_long$age_cat == cat,]))
  print(round(nrow(strata_long[strata_long$age_cat == cat,]) / nrow(strata_long), 4))
}

# province
for(prov in sort(unique(strata_long$province))){
  print(prov)
  print(nrow(strata_long[strata_long$province == prov,]))
  print(round(nrow(strata_long[strata_long$province == prov,]) / nrow(strata_long), 4))
}
