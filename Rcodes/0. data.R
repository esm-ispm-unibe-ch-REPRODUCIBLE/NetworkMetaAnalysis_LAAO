library(dplyr)
library(ggh4x)
library(ggplot)
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

outcome_names <- c("all_bleeding", "MajorBleeding", "DRT")

study_names <- data_raw$Name; nstudies <- length(study_names)
treatments <- unique(data_long$arm); ntrts <- length(treatments)

data_long %>% group_by(arm) %>% summarize(DRTs = sum(DRT),
                                          Bs = sum(all_bleeding),
                                          MBs = sum(MajorBleeding))

data <- list()
data[["Bleeding"]] <- data_long %>% 
                    select(Name, arm, all_bleeding, n_arm) %>% 
                    rename(events = all_bleeding)

data[["MajorBleeding"]] <- data_long %>% 
  select(Name, arm, MajorBleeding, n_arm) %>% 
  rename(events = MajorBleeding)

data[["DRT"]] <- data_long %>% 
  select(Name, arm, DRT, n_arm) %>% 
  rename(events = DRT)

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
            Major_bleeding = sum(MajorBleeding), DRT = sum(DRT)) %>% 
  mutate(bleeding_rate = all_bleeding/n,
         major_bleeding_rate = Major_bleeding/n,
         DRT_rate = DRT / n)
