# DRT
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

DRT_Bayesian_object <- NMA_compute_risk_differences(data_DRT_Bayesian)  

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

old_names <- colnames(DRT.mcmc_to_plot[[1]])

# Nuovi nomi
new_names <- old_names
idx <- match(names(rd_map), old_names)

# Applica a tutte le 5 catene
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

