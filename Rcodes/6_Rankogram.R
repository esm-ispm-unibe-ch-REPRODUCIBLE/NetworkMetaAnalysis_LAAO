# Rankograms
# All bleeding
bleeding_treatment_effectiveness <- NMA_bleeding_Bayesian %>% 
  select(starts_with("effectiveness")) %>% 
  colMeans() 

bleeding_treatment_rankings <- data.frame(
  outcome = "Any Bleeding",
  treatment = rep(1:3, times = 3),
  ranking = rep(1:3, each = 3),
  probability = bleeding_treatment_effectiveness
)


# Major bleeding
majorbleeding_treatment_effectiveness <- NMA_MajorBleeding_Bayesian %>% 
  select(starts_with("effectiveness")) %>% 
  colMeans() 

majorbleeding_treatment_rankings <- data.frame(
  outcome = "Major Bleeding",
  treatment = rep(1:3, times = 3),
  ranking = rep(1:3, each = 3),
  probability = majorbleeding_treatment_effectiveness
)


# DRT
DRT_treatment_effectiveness <- NMA_DRT_Bayesian %>% 
  select(starts_with("effectiveness")) %>% 
  colMeans() 

DRT_treatment_rankings <- data.frame(
  outcome = "Device-related Thrombosis",
  treatment = rep(1:3, times = 3),
  ranking = rep(1:3, each = 3),
  probability = DRT_treatment_effectiveness
)

# Mortality
mortality_treatment_effectiveness <- NMA_Mortality_Bayesian %>% 
  select(starts_with("effectiveness")) %>% 
  colMeans() 

mortality_treatment_rankings <- data.frame(
  outcome = "All mortality",
  treatment = rep(1:3, times = 3),
  ranking = rep(1:3, each = 3),
  probability = mortality_treatment_effectiveness
)

# Rankograms and plots
treatment_rankings <- bind_rows(bleeding_treatment_rankings, majorbleeding_treatment_rankings,
                                DRT_treatment_rankings, mortality_treatment_rankings)
treatment_rankings$treatment_factor <- factor(x = treatment_rankings$treatment,
                                              levels = unique(treatment_rankings$treatment),
                                              labels = c("DAPT", "Low-dose DOAC", "DOAC"))
treatment_rankings$treatment_string <- as.character(treatment_rankings$treatment_factor)

# library(ggplot2)
library(ggh4x)
rankogram <- ggplot(data = treatment_rankings, 
                    aes(x = as.factor(ranking), y = probability, col = outcome,
                        group = outcome, linetype = outcome)) +
  facet_wrap2(as.factor(treatment_string) ~ ., ncol = 1, strip.position = "top",
              axes = "all", remove_labels = "none") +
  geom_line(linewidth = 0.6) +
  geom_point() +
  scale_color_manual(values = c("salmon", "magenta", "dodgerblue", "indianred")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Ranking", y = "Probability", col = "", linetype = "") 

print(rankogram)




