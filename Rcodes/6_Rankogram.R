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

# Rankograms and plots
treatment_rankings <- bind_rows(bleeding_treatment_rankings, majorbleeding_treatment_rankings)
treatment_rankings$treatment_factor <- factor(x = treatment_rankings$treatment,
                                              levels = 1:3,
                                              labels = c("DAPT", "Low-dose DOAC", "DOAC"))

library(ggplot2)
library(ggh4x)
rankogram <- ggplot(data = treatment_rankings, 
                    aes(x = as.factor(ranking), y = probability, col = outcome,
                        group = outcome, linetype = outcome)) +
  facet_wrap2(treatment_factor ~ ., ncol = 1, strip.position = "top",
              axes = "all", remove_labels = "none") +
  geom_line(linewidth = 0.6) +
  geom_point() +
  scale_color_manual(values = c("salmon", "magenta")) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Ranking", y = "Probability", col = "", linetype = "") 

print(rankogram)




