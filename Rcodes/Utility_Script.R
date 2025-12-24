# Utility Script
build_nma_data <- function(df, ref_trt = 1){
  ## studies and treatments
  study_ids <- sort(unique(df$ID_study)) # IDs of the studies
  trt_codes <- sort(unique(df$trt)) # IDs of the treatments
  
  nstudies <- length(study_ids) # Number of studies
  ntrts <- length(trt_codes) # number of treatments
  max_arms <- max(table(df$ID_study)) # Number of arms
  
  # matrix of the treatments
  matrix_of_treatments <- matrix(NA_integer_, nrow = nstudies, ncol = max_arms)
  # matrix of the sizes for each treatment in each study and of events
  number_of_patients <- number_of_events <- matrix(NA_integer_, nrow = nstudies, ncol = ntrts)
  
  # creates a vector of 0 with length equal to the number of studies
  # each element will be the number of arms for each study
  number_of_arms <- integer(nstudies)
  
  # fill the matrixes 
  for(i in 1:nrow(df)){
    # which study are we considering?
    s <- study_ids[df$ID_study[i]]
    # create k, an intermediate variable given by the unitary increase of the
    # current number of arms for that study in the vector na
    k <- number_of_arms[s] + 1 # next free arm slot
    # which treatment are we considering?
    tr <- trt_codes[df$trt[i]]
    
    # t-matrix
    # for the s-th study and the k-th arm I fill with the trt code 
    matrix_of_treatments[s, k] <- tr
    
    # events and sample matrixes
    # for the s-th study and k-th treatment I add the number of events per
    # treatment
    number_of_events[s, tr] <- df$events[i]
    # and also the number of patients for that treatment
    number_of_patients[s, tr] <- df$n_arm[i]
    
    number_of_arms[s] <- k # update the vector arm counter
  }
  
  to_return <- list(
    nstudies = nstudies,
    ntreatments = ntrts,
    narms = number_of_arms,
    treatments = matrix_of_treatments,
    patients = number_of_patients,
    events = number_of_events,
    reference_treatment = ref_trt # numeric reference code
  )
  
  return(to_return)
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
    #  logit(p[i, t[i, 1]]) <- u[i] # baseline risk
      
      for(k in 1:na[i]){
        logit(p[i, t[i, k]]) <- u[i] + d[t[i, k]] 
      }
    }
  
    # prior distribution for basic parameters
    # except for the reference which is set to 0
    for(k in c(2, 3)){
      d[k] ~ dnorm(0, 0.01)
    }
    d[ref] <- 0
    
    # OR for each comparison between treatment T(T-1)/2
    for(i in 1:(nt-1)){
      for(j in (i+1):nt){
        OR[j, i] <- exp(d[j] - d[i])
      }
    }
    
    odds_DAPT <- rate_DAPT / (1 - rate_DAPT)
    
    p_DOAC <- odds_DAPT / (OR[3, 1] + odds_DAPT)
    odds_DOAC <- p_DOAC / (1 - p_DOAC)
    
    pLowDoseDOAC <- odds_DAPT / (OR[3, 2] + odds_DAPT)
    odds_LowDoseDOAC <- pLowDoseDOAC / (1 - pLowDoseDOAC)
    
    # Risk differences
    
    # DOAC vs DAPT
    RD[1] <- ((OR[3, 1] * odds_DAPT) / (1 + OR[3, 1] * odds_DAPT) - rate_DAPT) * 100
    
    # LowDoseDOAC vs DAPT
    RD[2] <- ((OR[2, 1] * odds_DAPT) / (1 + OR[2, 1] * odds_DAPT) - rate_DAPT) * 100
    
    # DOAC vs LowDoseDOAC
    RD[3] <- ((OR[3, 2] * odds_DOAC) / (1 + OR[3, 2] * odds_DOAC) - p_DOAC) * 100
    
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
                       ref = data_object[["reference_treatment"]],
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
    # except for the reference which is set to 0
    for(k in 1:nt){
      d[k] ~ dnorm(0, 100)
    }
    # d[1] <- 0
    
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
        X = to_return)
  
  
  closeAllConnections()
  to_return <- list(samps, posterior_samples)
  
  return(to_return)
}

get_density <- function(x, varname){
  z <- density(x)
  data.frame(x = z$x, y = z$y, variable = varname)
}






