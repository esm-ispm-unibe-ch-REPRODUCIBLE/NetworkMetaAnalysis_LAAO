# Major Bleeding

# Frequentist NMA - Fixed effects
Major.Bleeding <- data_long %>% select(Name, arm, MajorBleeding,
                                 n_arm)

MB.data_pair <- pairwise(treat = arm, event = MajorBleeding,
                      n = n_arm, data = Major.Bleeding,
                      studlab = Name, sm = "OR")

MB_NMA_MH <- netmetabin(event1 = event1, n1 = n1,
                               event2 = event2, n2 = n2,
                               treat1 = treat1, treat2 = treat2,
                               studlab = Name, data = MB.data_pair,
                               sm = "OR", method = "MH",
                               common = T, reference.group = "DAPT")

summary(MB_NMA_MH)
netleague(MB_NMA_MH)

MB_NMA_LRP <- update(MB_NMA_MH, method = "LRP")
summary(MB_NMA_LRP)
netleague(MB_NMA_LRP)

MB_NMA_IV <- update(MB_NMA_MH, method = "Inverse")
summary(MB_NMA_IV); netleague(MB_NMA_IV)
(r.MB <- netsplit(MB_NMA_IV, method = "Back",
                  random = F)); plot(r.MB)
decomp.design(MB_NMA_IV)


# Bayesian NMA - fixed effects
data_MajorBleeding_Bayesian <- build_nma_data(df = data_recoded[[2]])

MajorBleeding_Bayesian_object <- NMA_compute(data_MajorBleeding_Bayesian,
                                          DAPT_rate_event = crude_events_rate$major_bleeding_rate[crude_events_rate$arm == "DAPT"])

NMA_MajorBleeding_Bayesian <- MajorBleeding_Bayesian_object[[2]] 

apply(MARGIN = 2, 
      FUN = function(x){
        return(c(mean(x) ,quantile(x, probs = c(0.025, 0.975))))
      } ,
      X = NMA_MajorBleeding_Bayesian)

# Applica a ciascuna colonna, includendo il nome
densities_Major_bleeding <- Map(get_density, NMA_MajorBleeding_Bayesian, 
                                names(NMA_MajorBleeding_Bayesian)) %>% 
  bind_rows %>% filter(startsWith(variable, "RD")) %>% 
  mutate(area = ifelse(x >= 0, "Control better", "Active better"),
         outcome = "Major Bleedings")

MB.mcmc_to_plot <- MajorBleeding_Bayesian_object[[1]]


# Convergence assessment
param_names <- colnames(NMA_MajorBleeding_Bayesian)

old_names <- colnames(MB.mcmc_to_plot[[1]])

new_names <- old_names
new_names[old_names == "RD[1]"] <- "RD: DOAC vs DAPT"
new_names[old_names == "RD[2]"] <- "RD: Low-dose DOAC vs DAPT"
new_names[old_names == "RD[3]"] <- "RD: DOAC vs Low-dose DOAC"


for (i in seq_along(bleeding.mcmc_to_plot)) {
  colnames(MB.mcmc_to_plot[[i]]) <- new_names
}

params_to_plot <- grep("^(RD)", new_names, value = TRUE);
n_params <- length(params_to_plot)

gelman.diag(MB.mcmc_to_plot[, params_to_plot])


# Layout: n_params righe, 2 colonne
par(mfrow = c(n_params, 2), mar = c(4, 4, 2, 1))
for(i in params_to_plot) {
  # Trace plot
  traceplot(MB.mcmc_to_plot[, i], main = paste("Trace of", i))
  
  # Density plot
  densplot(MB.mcmc_to_plot[, i], show.obs = FALSE, main = paste("Density of", i),
           xlab = "Parameters estimate (%)")
}
par(mfrow = c(1,1))  # Reset layout


# Meta-analysis
MajorBleeding <- MA_data_long %>% select(Name, arm_MA, MajorBleeding,
                                    n_arm)

# Frequentist meta-analysis
MB.data_pair.MA <- pairwise(treat = arm_MA, event = MajorBleeding,
                                  n = n_arm, data = MajorBleeding,
                                  studlab = Name, sm = "OR")

MB_MA_MH <- metabin(MB.data_pair.MA,
                          sm = "OR", method = "MH",
                          random = F, reference.group = "DAPT")

summary(MB_MA_MH)

MB_MA_RD <- update(MB_MA_MH, sm = "RD")
summary(MB_MA_RD)

# Bayesian meta-analysis
MA_data_MB_Bayesian <- build_nma_data(df = data_MA_recoded[[2]])
set.seed(42)
MB_MA_object <- MA_compute(data_object = MA_data_MB_Bayesian,
                                    DAPT_rate_event = crude_events_rate$major_bleeding_rate[crude_events_rate$arm == "DAPT"])

MA_MB_Bayesian <- MB_MA_object[[2]]

apply(MARGIN = 2, 
      FUN = function(x){
        return(c(median(x), quantile(x, probs = c(0.025, 0.975))))
      } ,
      X = MA_MB_Bayesian)

# Convergence assessment
plot(MB_MA_object[[1]])
gelman.diag(MB_MA_object[[1]])

mean(MA_MB_Bayesian$RD < 0)


