# All mortality
# Frequentist NMA - Fixed effect
# Overall mortality: Frequentist NMA
data_long %>% group_by(arm) %>% summarize(events_per_arm = sum(AllMortality))

AM.data_pair <- pairwise(treat = arm, event = AllMortality,
                         n = n_arm, data = data_long,
                         studlab = Name, sm = "RD") 

AM_NMA <- netmetabin(AM.data_pair, sm = "RD", method = "Inverse",
                     common = T, reference.group = "DAPT")

summary(AM_NMA)
netleague(AM_NMA)
(r.AM <- netsplit(AM_NMA, method = "Back",
                  random = F)); plot(r.AM, digits = 4)

decomp.design(AM_NMA)

# Bayesian NMA - Fixed effects
data_Mortality_Bayesian <- build_nma_data(df = data_recoded[[4]])

set.seed(42)
Mortality_Bayesian_object <- NMA_compute_risk_differences(data_Mortality_Bayesian)  

NMA_Mortality_Bayesian <- Mortality_Bayesian_object[[2]]

apply(MARGIN = 2, 
      FUN = function(x){
        return(c(mean(x) ,quantile(x, probs = c(0.025, 0.975))))
      } ,
      X = NMA_Mortality_Bayesian)

densities_mortality <- Map(get_density, NMA_Mortality_Bayesian, names(NMA_Mortality_Bayesian)) %>% 
  bind_rows %>% filter(startsWith(variable, "RD")) %>% 
  mutate(area = ifelse(x >= 0, "Control better", "Active better"),
         outcome = "All mortality")

mortality.mcmc_to_plot <- Mortality_Bayesian_object[[1]]


# Convergence assessment
rd_order <- c("RD[3,1]", "RD[2,1]", "RD[3,2]")
rd_map <- c(
  "RD[3,1]" = "RD: DOAC vs DAPT",
  "RD[2,1]" = "RD: Low-dose DOAC vs DAPT",
  "RD[3,2]" = "RD: DOAC vs Low-dose DOAC"
)

# Nomi originali (una catena basta)
old_names <- colnames(DRT.mcmc_to_plot[[1]])

# Nuovi nomi
new_names <- old_names
idx <- match(names(rd_map), old_names)
new_names[idx[!is.na(idx)]] <- rd_map[names(rd_map)[!is.na(idx)]]

# Applica a tutte le 5 catene
for (i in seq_along(DRT.mcmc_to_plot)) {
  colnames(mortality.mcmc_to_plot[[i]]) <- new_names
}

# Plot
params_to_plot <- c(
  "RD: DOAC vs DAPT",
  "RD: Low-dose DOAC vs DAPT",
  "RD: DOAC vs Low-dose DOAC"
)

n_params <- length(params_to_plot)

par(mfrow = c(n_params, 2), mar = c(4, 4, 2, 1))

for(i in params_to_plot) {
  # Trace plot
  traceplot(mortality.mcmc_to_plot[, i], main = paste("Trace of", i))
  
  # Density plot
  densplot(mortality.mcmc_to_plot[, i], show.obs = FALSE, main = paste("Density of", i),
           xlab = "Parameters estimate (%)")
}

par(mfrow = c(1,1))  



gelman.diag(mortality.mcmc_to_plot[, params_to_plot])