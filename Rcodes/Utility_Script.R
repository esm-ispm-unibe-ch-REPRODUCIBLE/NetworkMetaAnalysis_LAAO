# Utility Script
build_nma_data <- function(df, ref_trt = 1){
  
  study_ids <- sort(unique(df$ID_study))
  trt_codes <- sort(unique(df$trt))
  
  nstudies <- length(study_ids)
  ntrts <- length(trt_codes)
  max_arms <- max(table(df$ID_study))
  
  t_mat <- matrix(NA_integer_, nstudies, max_arms)
  r_mat <- matrix(NA_integer_, nstudies, ntrts)
  n_mat <- matrix(NA_integer_, nstudies, ntrts)
  na <- integer(nstudies)
  
  for(i in seq_len(nrow(df))){
    
    s  <- match(df$ID_study[i], study_ids)
    tr <- match(df$trt[i], trt_codes)
    
    k <- na[s] + 1
    t_mat[s, k] <- tr
    
    r_mat[s, tr] <- df$events[i]
    n_mat[s, tr] <- df$n_arm[i]
    
    na[s] <- k
  }
  
  list(
    nstudies = nstudies,
    ntreatments = ntrts,
    narms = na,
    treatments = t_mat,
    events = r_mat,
    patients = n_mat,
    reference_treatment = ref_trt
  )
}


NMA_compute <- function(data_object, DAPT_rate_event){
  model.Bayesian_NMA <- "
  model{
    for(i in 1:ns){
      # Prior distribution for log-odds in baseline arm of study i
      u[i] ~ dnorm(0, 0.01)   
      
      # Binomial likelihood for number of events for each arm k for study i
      for(k in 1:na[i]){
        r[i, t[i, k]] ~ dbin(p[i, t[i, k]], n[i, t[i, k]])
      }
      
      # Parametrization of the true effect of each comparison
      # of arm k vs baseline arm (1) of study y   
      logit(p[i, t[i, 1]]) <- u[i] # baseline risk
      
      for(k in 2:na[i]){
        logit(p[i, t[i, k]]) <- u[i] + d[t[i, k]] - d[t[i, 1]]
      }
    }
  
    # prior distribution for basic parameters
    # except for the reference which is set to 0
    for(l in 2:nt){
      d[l] ~ dnorm(0, 0.01)
    }
    d[1] <- 0
    
    # OR for each comparison between treatment T(T-1)/2
    for(i in 1:(nt-1)){
      for(j in (i+1):nt){
        OR[j, i] <- exp(d[j] - d[i])
      }
    }
    
    odds_DAPT <- rate_DAPT / (1 - rate_DAPT)
    
    pLowDoseDOAC <- (OR[2, 1] * odds_DAPT) / (1 + OR[2, 1] * odds_DAPT)
    odds_LowDoseDOAC <- pLowDoseDOAC / (1 - pLowDoseDOAC)
    
    # Risk differences
    
    # DOAC vs DAPT
    RD[1] <- ((OR[3, 1] * odds_DAPT) / (1 + OR[3, 1] * odds_DAPT) - rate_DAPT) * 100
    
    # LowDoseDOAC vs DAPT
    RD[2] <- (pLowDoseDOAC - rate_DAPT) * 100
    
    # DOAC vs LowDoseDOAC
    RD[3] <- ((OR[3, 2] * odds_LowDoseDOAC) / (1 + OR[3, 2] * odds_LowDoseDOAC) - pLowDoseDOAC) * 100
    
  # Ranking of treatments
  order[1:nt] <- rank(d[1:nt]) # order by effect: ----> best treatment
  # smaller value means protective effect
  
  # treatment, rank
  for(i in 1:nt){    
    for(j in 1:nt){
      effectiveness[i, j] <- equals(order[i], j)
    }
  }
  
  }
"
  
  data_to_pass <- list(ns = data_object[["nstudies"]], 
                       r = data_object[["events"]], 
                       n = data_object[["patients"]],
                       nt = data_object[["ntreatments"]],
                       na = data_object[["narms"]],
                       t = data_object[["treatments"]],
                       rate_DAPT = DAPT_rate_event)
  
  model.NMA.spec <- textConnection(model.Bayesian_NMA)
  MA.jags.model <- jags.model(file = model.NMA.spec, data = data_to_pass,
                              n.chains = 5, quiet = T)
  
  params <- c("OR", "effectiveness", "RD")
  samps <- coda.samples(MA.jags.model, variable.names = params, n.iter = 5e4, 
                        thin = 5, progress.bar = "none")
  posterior_samples <- do.call(rbind, samps) %>% as.data.frame() 
  
  to_return <- list(samps, posterior_samples)
  
  closeAllConnections()
  
  return(to_return)
}


NMA_compute_risk_differences <- function(data_object){
  model.Bayesian_NMA <- "
  model{
    for(i in 1:ns){
      u[i] ~ dbeta(0.5, 0.5)
      
      # Binomial likelihood for number of events for each arm k for study i
      for(k in 1:na[i]){
        r[i, t[i, k]] ~ dbin(p[i, t[i, k]], n[i, t[i, k]])
      }
      
      # Parametrization of the true effect of each comparison
      for(k in 1:na[i]){
        p[i, t[i, k]] <- max(min(u[i] + d[t[i, k]], 0.999), 0.001)
      }
    }
  
    # prior distribution for basic parameters
    for(k in 2:nt){
      d[k] ~ dnorm(0, 100)
    }
    d[1] <- 0
    
    # RD for each comparison between treatment T(T-1)/2
    for(i in 1:(nt-1)){
      for(j in (i+1):nt){
        RD[j, i] <- (d[j] - d[i]) * 100
      }
    }
    
  # Ranking of treatments
  order[1:nt] <- rank(d[1:nt]) 
  # if negative we don't add anything because a smaller value has a protective
  # effect ----> best treatment
  
  for(i in 1:nt){
    for(j in 1:nt){
      effectiveness[i, j] <- equals(order[i], j)
    }
  }

  }
"
  
  input_data <- list(ns = data_object[["nstudies"]], 
                     r = data_object[["events"]], 
                     n = data_object[["patients"]],
                     nt = data_object[["ntreatments"]],
                     na = data_object[["narms"]],
                     t = data_object[["treatments"]]
  )
  
  model.NMA.spec <- textConnection(model.Bayesian_NMA)
  MA.jags.model <- jags.model(file = model.NMA.spec, data = input_data,
                              n.chains = 5, quiet = T)
  
  params <- c("RD", "effectiveness")
  samps <- coda.samples(MA.jags.model, variable.names = params, n.iter = 5e4, 
                        thin = 5, progress.bar = "none")
  posterior_samples <- do.call(rbind, samps) %>% as.data.frame()
  
  apply(MARGIN = 2, 
        FUN = function(x){
          return(c(mean(x) ,quantile(x, probs = c(0.025, 0.975))))
        } ,
        X = posterior_samples)
  
  
  closeAllConnections()
  to_return <- list(samps, posterior_samples)
  
  return(to_return)
}





