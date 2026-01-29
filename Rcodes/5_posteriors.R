# Posteriors probability

# Probability that that active is better than control (p<0)


NMA_bleeding_Bayesian %>% 
  summarize(prob1 = mean(`RD[1]` > 0),
            prob2 = mean(`RD[2]` > 0),
            prob3 = mean(`RD[3]` > 0))

NMA_MajorBleeding_Bayesian %>% 
  summarize(prob1 = mean(`RD[1]` > 0),
            prob2 = mean(`RD[2]` > 0),
            prob3 = mean(`RD[3]` > 0))


NMA_DRT_Bayesian %>% 
  summarize(prob1 = mean(`RD[3,1]` > 0))

NMA_Mortality_Bayesian %>% 
  summarize(prob1 = mean(`RD[3,1]` > 0))


# PLOT
densities_df <- bind_rows(densities_all_bleeding, densities_Major_bleeding,
                          densities_drt, densities_mortality) %>%
  mutate(
    comparison = case_when(
      variable %in% c("RD[1]", "RD[3,1]") ~ "DOAC vs DAPT",
      variable %in% c("RD[2]", "RD[2,1]") ~ "LowDoseDOAC vs DAPT",
      TRUE                               ~ "DOAC vs LowDoseDOAC"
    ),
    outcome_order = case_when(
      outcome == "All mortality" ~ 1,
      outcome == "Any bleeding" ~ 2,
      outcome == "Major bleeding" ~ 3,
      outcome == "Device-related Thrombosis" ~ 4
    )
  )

densities_df$outcome_factor <- factor(densities_df$outcome_order,
                                      1:4,
                                      c("All mortality" ,"Any bleeding",
                                        "Major bleeding", "Device-related Thrombosis"))

densities_df$comparison_order <- ifelse(densities_df$comparison == "DOAC vs DAPT", 1,
                                        ifelse(densities_df$comparison == "LowDoseDOAC vs DAPT", 2, 3))

densities_df$comparison_factor <- factor(densities_df$comparison_order,
                                         1:3,
                                         c("DOAC vs DAPT", "Low-dose DOAC vs DAPT", "DOAC vs Low-dose DOAC"))

# library(dplyr)
densities_df <- densities_df %>% arrange(outcome_order, comparison_order)

posteriors <- ggplot(data = densities_df, aes(x = x, ymin = 0, ymax = y, fill = area,
                                              colour = area)) +
  geom_ribbon() +
  facet_wrap(~ comparison_factor + outcome_factor, scales = "free", axes = "all",
             ncol = 4) +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = 0, color = 'darkmagenta', linetype = "dashed") +
  labs(fill = "", y= "Posterior probability density", x = "Risk Difference",
       color = "") +
  scale_x_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5),
    legend.position = "bottom"
  ) 

print(posteriors)
