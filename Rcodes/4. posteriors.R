# Posteriors probability

# Probability that that active is better than control (p<0)

NMA_bleeding_Bayesian %>% 
  summarize(prob1 = mean(`RD[1]` < 0),
            prob2 = mean(`RD[2]` < 0),
            prob3 = mean(`RD[3]` < 0))

NMA_MajorBleeding_Bayesian %>% 
  summarize(prob1 = mean(`RD[1]` < 0),
            prob2 = mean(`RD[2]` < 0),
            prob3 = mean(`RD[3]` < 0))


NMA_DRT_Bayesian %>% 
  summarize(prob1 = mean(`RD[3,1]` < 0),
            prob2 = mean(`RD[2,1]` < 0),
            prob3 = mean(`RD[3,2]` < 0))


# Plot
densities_df <- bind_rows(densities_all_bleeding, densities_Major_bleeding,
                          densities_drt)%>%
  mutate(
    comparison = case_when(
      variable %in% c("RD[1]", "RD[3,1]") ~ "DOAC vs DAPT",
      variable %in% c("RD[2]", "RD[2,1]") ~ "LowDoseDOAC vs DAPT",
      TRUE                               ~ "DOAC vs LowDoseDOAC"
    )
  )

library(ggplot2)
ggplot(data = densities_df, aes(x = x, ymin = 0, ymax = y, fill = area,
                                colour = area)) +
  geom_ribbon() +
  facet_wrap(~ comparison + outcome, scales = "free", axes = "all",
             ncol = 3) +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = 0, color = 'darkmagenta', linetype = "dashed") +
  labs(fill = "", y= "Posterior probability density", x = "Risk Difference",
       color = "", title = "Network meta-analysis", subtitle = "Bayesian fixed effects model") +
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
