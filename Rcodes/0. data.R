library(dplyr)
library(ggh4x)
library(ggplot)
library(meta)
library(netmeta)
library(readxl)
library(tidyr)
library(tidyselect)

data_raw <- read_xlsx("NMA_data.xlsx")

data_long <- data_raw %>%
  pivot_longer(
    cols = -Name,   
    names_to = ".value",
    names_pattern = "(.+?)_?\\d+$"
  ) %>% 
  filter(is.na(arm) == F) 

outcome_names <- c("all_bleeding", "MajorBleeding", "DRT", "AllMortality")

study_names <- data_raw$Name; nstudies <- length(study_names)
treatments <- unique(data_long$arm); ntrts <- length(treatments)

data <- vector("list", length = length(outcome_names)) 
data[["Bleeding"]] <- data_long %>% 
                    select(Name, arm, all_bleeding, n_arm) %>% 
                    rename(events = all_bleeding)

data[["MajorBleeding"]] <- data_long %>% 
  select(Name, arm, MajorBleeding, n_arm) %>% 
  rename(events = MajorBleeding)

data[["DRT"]] <- data_long %>% 
  select(Name, arm, DRT, n_arm) %>% 
  rename(events = DRT)

data[["AllMortality"]] <- data_long %>% 
  select(Name, arm, AllMortality, n_arm) %>% 
  rename(events = AllMortality)

# Set data for Bayesian Analysis
data_recoded <- data
for(i in 1:length(data)){
  data_recoded[[i]] <- data_recoded[[i]] %>% 
    mutate(ID_study = match(Name, study_names),
           trt = match(arm, treatments))
}

# Calculate crude events rate
crude_events_rate <- data_long %>%  group_by(arm) %>%
  summarize(n = sum(n_arm), all_bleeding = sum(all_bleeding),
            Major_bleeding = sum(MajorBleeding)) %>% 
  mutate(bleeding_rate = all_bleeding/n,
         major_bleeding_rate = Major_bleeding/n)

# Data for meta-analysis
MA_data_long <- data_long %>% 
  mutate(arm_MA = ifelse(arm == "DAPT", "DAPT", "DOAC")) %>% 
  group_by(Name, arm_MA) %>% 
  summarize(n_arm = sum(n_arm),
            all_bleeding = sum(all_bleeding),
            MajorBleeding = sum(MajorBleeding),
            DRT = sum(DRT),
            AllMortality = sum(AllMortality))


# Bayesian analysis on the OR scale
data_MA <- list()
data_MA[["Bleeding"]] <- MA_data_long %>% 
  select(Name, arm_MA, all_bleeding, n_arm) %>% 
  rename(events = all_bleeding)

data_MA[["MajorBleeding"]] <- MA_data_long %>% 
  select(Name, arm_MA, MajorBleeding, n_arm) %>% 
  rename(events = MajorBleeding)

data_MA[["DRT"]] <- MA_data_long %>% 
  select(Name, arm_MA, DRT, n_arm) %>% 
  rename(events = DRT)

data_MA[["AllMortality"]] <- MA_data_long %>% 
  select(Name, arm_MA, AllMortality, n_arm) %>% 
  rename(events = AllMortality)

# Bayesian NMA
data_MA_recoded <- data_MA

for(i in 1:length(data)){
  data_MA_recoded[[i]] <- data_MA_recoded[[i]] %>% 
    mutate(ID_study = match(Name, study_names),
           trt = match(arm_MA, treatments))
}






  
