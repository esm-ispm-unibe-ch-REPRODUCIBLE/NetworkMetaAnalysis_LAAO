# Device-related Thrombosis
# Frequentist NMA - Fixed effect
DRT.data_pair <- pairwise(treat = arm, event = DRT,
                         n = n_arm, data = data_long,
                         studlab = Name, sm = "RD") 

DRT_NMA <- netmetabin(DRT.data_pair, sm = "RD", method = "Inverse",
                        common = T, reference.group = "DAPT")

summary(DRT_NMA)
netleague(DRT_NMA)
(r.DRT <- netsplit(DRT_NMA, method = "Back",
                  random = F)); plot(r.DRT, digits = 4)
decomp.design(DRT_NMA)

# Bayesian NMA - Fixed Effects
data_DRT_Bayesian <- build_nma_data(df = data_recoded[[3]])

DRT_Bayesian_object <- NMA_compute(data_DRT_Bayesian)  

NMA_DRT_Bayesian <- DRT_Bayesian_object[[2]]

apply(MARGIN = 2, 
      FUN = function(x){
        return(c(mean(x) ,quantile(x, probs = c(0.025, 0.975))))
      } ,
      X = NMA_DRT_Bayesian)

densities_drt <- Map(get_density, NMA_DRT_Bayesian, names(NMA_DRT_Bayesian)) %>% 
  bind_rows %>% filter(startsWith(variable, "RD")) %>% 
  mutate(area = ifelse(x >= 0, "Control better", "Active better"),
         outcome = "DRT")

DRT.mcmc_to_plot <- DRT_Bayesian_object[[1]]

# Convergence assessment
rd_order <- c("RD[3,1]", "RD[2,1]", "RD[3,2]")
rd_map <- c(
  "RD[3,1]" = "RD: DOAC vs DAPT",
  "RD[2,1]" = "RD: Low-dose DOAC vs DAPT",
  "RD[3,2]" = "RD: DOAC vs Low-dose DOAC"
)

new_names <- old_names <- colnames(DRT.mcmc_to_plot[[1]])
idx <- match(names(rd_map), old_names)

for (i in seq_along(DRT.mcmc_to_plot)) {
  colnames(DRT.mcmc_to_plot[[i]]) <- new_names
}

# Plot
params_to_plot <- c(
  "RD: DOAC vs DAPT",
  "RD: Low-dose DOAC vs DAPT",
  "RD: DOAC vs Low-dose DOAC"
); n_params <- length(params_to_plot)

# Set Layout
par(mfrow = c(n_params, 2), mar = c(4, 4, 2, 1))

for(i in params_to_plot) {
  # Trace plot
  traceplot(DRT.mcmc_to_plot[, i], main = paste("Trace of", i))
  
  # Density plot
  densplot(DRT.mcmc_to_plot[, i], show.obs = FALSE, main = paste("Density of", i),
           xlab = "Parameters estimate (%)")
}
par(mfrow = c(1,1))  # Reset layout

gelman.diag(DRT.mcmc_to_plot[, params_to_plot])

# Meta-analysis
DRT <- MA_data_long %>% select(Name, arm_MA, DRT, n_arm)

# Frequentist meta-analysis
drt.data_pair.MA <- pairwise(treat = arm_MA, event = DRT,
                            n = n_arm, data = DRT,
                            studlab = Name, sm = "OR")

DRT_MA_MH <- metabin(drt.data_pair.MA,
                    sm = "OR", method = "MH",
                    random = F, reference.group = "DAPT")

summary(DRT_MA_MH)

DRT_MA_RD <- update(DRT_MA_MH, sm = "RD")
summary(DRT_MA_RD)

# Bayesian meta-analysis
MA_data_DRT_Bayesian <- build_nma_data(df = data_MA_recoded[[3]])
set.seed(42)
DRT_MA_object <- MA_compute(data_object = MA_data_DRT_Bayesian,
                           DAPT_rate_event = crude_events_rate$DRT_rate[crude_events_rate$arm == "DAPT"])

MA_DRT_Bayesian <- DRT_MA_object[[2]]

apply(MARGIN = 2, 
      FUN = function(x){
        return(c(median(x), quantile(x, probs = c(0.025, 0.975))))
      } ,
      X = MA_DRT_Bayesian)

# Convergence assessment
plot(DRT_MA_object[[1]])
gelman.diag(DRT_MA_object[[1]])

mean(MA_DRT_Bayesian$RD < 0)
c <- get_density(MA_DRT_Bayesian$RD, "RD") %>% mutate(outcome = "Device-related Thrombus")


ggplot(data = c, aes(x = x, ymin = 0, ymax = y, fill = area, colour = area)) +
    geom_ribbon() +
    geom_line(aes(y = y)) +
    geom_vline(xintercept = 0, color = 'darkmagenta', linetype = "dashed") +
    labs(fill = "", y= "Posterior probability density", x = "Risk Difference",
         color = "") +
    scale_x_continuous(
      labels = function(x) paste0(x, "%")
    ) +
    scale_fill_manual(values = c("dodgerblue", "indianred2"))  + 
    scale_color_manual(values = c("dodgerblue",  "indianred2")) +
    theme(
      axis.title.x = element_text(color="black", size = 12),
      axis.title.y = element_text(color="black", size = 12),
      axis.text.x  = element_text(color = "black"),
      legend.text = element_text(size = 12),
      legend.position = "bottom",
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key        = element_rect(fill = "white", colour = NA),
      text = element_text(family = "Helvetica")
    )




