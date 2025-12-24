# All bleeding

## Frequentist Meta-Analysis

Bleeding <- data_long %>% select(Name, arm, all_bleeding,
                                     n_arm)
Bleeding.data_pair <- pairwise(treat = arm, event = all_bleeding,
                      n = n_arm, data = Bleeding,
                      studlab = Name, sm = "OR")

Bleeding_NMA_MH <- netmetabin(event1 = event1, n1 = n1,
                               event2 = event2, n2 = n2,
                               treat1 = treat1, treat2 = treat2,
                               studlab = Name, data = Bleeding.data_pair,
                               sm = "OR", method = "MH",
                               common = T, reference.group = "DAPT")
summary(Bleeding_NMA_MH);
netleague(Bleeding_NMA_MH)

Bleeding_NMA_LRP <- update(Bleeding_NMA_MH, method = "LRP")
summary(Bleeding_NMA_LRP); netleague(Bleeding_NMA_LRP)

Bleeding_NMA_IV <- update(Bleeding_NMA_MH, method = "Inverse") 
summary(Bleeding_NMA_IV); netleague(Bleeding_NMA_IV)
r <- netsplit(Bleeding_NMA_IV, method = "Back-calculation",
              random = F); s
decomp.design(Bleeding_NMA_IV)



# Bayesian NMA - fixed effects

# All bleedings
data_bleeding_Bayesian <- build_nma_data(df = data_recoded[[1]])

bleeding_Bayesian_object <- NMA_compute(data_bleeding_Bayesian,
                                        DAPT_rate_event = crude_events_rate$bleeding_rate[crude_events_rate$arm == "DAPT"])

NMA_bleeding_Bayesian <- bleeding_Bayesian_object[[2]] 

apply(MARGIN = 2, 
           FUN = function(x){
             return(c(mean(x) ,quantile(x, probs = c(0.025, 0.975))))
           } ,
           X = NMA_bleeding_Bayesian)

bleeding.mcmc_to_plot <- bleeding_Bayesian_object[[1]]



# Get densities
densities_all_bleeding <- Map(get_density, NMA_bleeding_Bayesian, 
                 names(NMA_bleeding_Bayesian)) %>% 
  bind_rows %>% filter(startsWith(variable, "RD")) %>% 
  mutate(area = ifelse(x >= 0, "Control better", "Active better"),
         outcome = "Bleedings")

# Convergence assessment
param_names <- colnames(NMA_bleeding_Bayesian)

old_names <- colnames(bleeding.mcmc_to_plot[[1]])

new_names <- old_names
new_names[old_names == "RD[1]"] <- "RD: DOAC vs DAPT"
new_names[old_names == "RD[2]"] <- "RD: Low-dose DOAC vs DAPT"
new_names[old_names == "RD[3]"] <- "RD: DOAC vs Low-dose DOAC"

for (i in seq_along(bleeding.mcmc_to_plot)) {
  colnames(bleeding.mcmc_to_plot[[i]]) <- new_names
}

params_to_plot <- grep("^(RD)", new_names, value = TRUE)
gelman.diag(bleeding.mcmc_to_plot[, params_to_plot])

# Set Layout
par(mfrow = c(length(params_to_plot), 2), mar = c(4, 4, 2, 1))

for(i in params_to_plot) {
  # Trace plot
  traceplot(bleeding.mcmc_to_plot[, i], main = paste("Trace of", i))
  
  # Density plot
  densplot(bleeding.mcmc_to_plot[, i], show.obs = FALSE, main = paste("Density of", i),
           xlab = "Parameters estimate (%)")
}

par(mfrow = c(1,1))  # Reset layout


















