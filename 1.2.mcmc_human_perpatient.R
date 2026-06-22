## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(gridExtra)
library("readxl")
source("helper.R")
source("0.human_PBPK_model_4brain.R")

#--------------Maximum likelihood estimation (MLE) function for MCMC----

for (i in 1:200000) {
  p <- rnorm(length(theta.MCMC[-which_sig]), mean = theta.MCMC[-which_sig], sd = 0.1)
  names(p) = names(theta.MCMC)[-which_sig]
  out <- tryCatch({
    mod.human.brain %>% param(exp(p)) %>% mrgsim()
  }, error = function(e) e)
  if (inherits(out, "error")) print(paste("crash at sample", i))
}


mod.human.brain<- mcode ("brain", human.4B.code)

# load the human clinical data
Obs.df.human = read_excel("/Users/wuji/work/code/brain_PBPK/biodistribution_human_patient.xlsx")
Obs.df.human = data.frame(Obs.df.human)[,c(5,9,10)]
colnames(Obs.df.human) <- c("patient_id","Time", "C_Blood")

# use the  human physiological parameters from clinical data (181913)

tinterval = 1
TDoses = 1
DOSE_human = 5e8
V_Blood_human  = 5300   #:  ml
V_Tumor_human = 26.36 #ml, got from paper https://www.liebertpub.com/doi/10.1089/hum.2018.001#sec-3 33.5 ml 
sig_names <- c("sig2",paste0("sig_", names(params.init.human)))

#-------------- fitting --------
params2fit.human =params.init.human

tsamp_human=tgrid(0,max(Obs.df.human$Time),1)     ## Simulation time 24*7 hours (180 days)

# Get unique patient IDs
unique_patients <- unique(Obs.df.human$patient_id)
for (pid in unique_patients) {
#for (pid in c(201,202,204,205,209,213,217)){  
#for (pid in c(205)){
  if (pid == 216){next}
  patient_data <- subset(Obs.df.human, patient_id == pid)[,2:3]
  
  patient_data <- subset(patient_data, Time != 0)
  
  # ---- Process patient_data here ----
  print(paste("Processing:", pid))
  print(patient_data)
  # 205 Port: 1.23, nelder-mead:1.37, sann: 8.28
  # 213 nelder-mead: 1.076
  if (pid %in% c(207, 211)) {
    Fit.Result.A1<- modFit(f=MCcost.human, p=params2fit.human, obs=patient_data, 
                           method ="L-BFGS-B", control = nls.lm.control(nprint=1)) #"Nelder-Mead" L-BFGS-B for 207
  }else{
    Fit.Result.A1<- modFit(f=MCcost.human, p=params2fit.human, obs=patient_data, 
                         method ="Nelder-Mead", control = nls.lm.control(nprint=1)) #"Nelder-Mead" L-BFGS-B for 207
  }
  res.A1=MCcost.human(Fit.Result.A1$par, obs=patient_data)$residuals$res     ## Check the residual for each time points
  sum(res.A1^2) 
  
  #Fitted_output.A1 = pred.human.iv(par=Fit.Result.A1$par)
  
  Fitted_output.A1 = pred.human.iv(par=Fit.Result.A1$par,DOSE_human) 

  print(round(max(Fitted_output.A1$C_Blood/1e5)+1)*1e5)
  pdf(
    file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), paste0("initial_fitting_plot_", pid, ".pdf")),
    width = 6,height = 5
  )
  
  plot(patient_data$Time, patient_data$C_Blood, type = "p", col = "red", 
       xlab = "Time", ylab = "Concentration in Blood", 
       main = "Observed vs Fitted Blood Concentration", ylim = c(0, round(max(Fitted_output.A1$C_Blood/1e5)+1)*1e5))
  lines(Fitted_output.A1$Time, Fitted_output.A1$C_Blood, col = "green", type = "l")
  
  dev.off()

  
  print(round(max(Fitted_output.A1$C_Blood/1e5)+1)*1e5)
  
  pdf(file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), paste0("initial_fit_C_TE_", pid, ".pdf")),
    width = 5, height = 4)
  plot(Fitted_output.A1$Time, Fitted_output.A1$C_tumor_E, type = "l", col = "lightblue", 
       xlab = "Time", ylab = "Concentration in tumor extra vascular", 
       ylim = c(0, round(max(Fitted_output.A1$C_tumor_E/1e7, na.rm = TRUE)+1)*1e7))
  dev.off()
  

  pdf(file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), paste0("initial_fit_C_TV_", pid, ".pdf")),
      width = 5, height = 4)
  
  plot(Fitted_output.A1$Time, Fitted_output.A1$C_tumor_V, type = "l", col = "lightblue", 
       xlab = "Time", ylab = "Concentration in tumor vascular", 
        ylim = c(0, round(max(Fitted_output.A1$C_tumor_V/1e5, na.rm = TRUE)+1)*1e5))
  #lines(Fitted_output.A1$Time, Fitted_output.A1$C_tumor_E, col = "green", type = "l")
  dev.off()

  # ------------- mcmc
  params <- Fit.Result.A1$par

  #-------------sig value for parameters and errors between obs and pred
  
  # Define parameter distribution 
  sig_population = 0.3
  sig_values <- rep(sig_population, length(sig_names))
  names(sig_values) = sig_names
  
  # error with prediction and observation's deviation variance
  sig_pred = 0.1
  sig_values["sig2"] <- sig_pred
  
  # Create sig_list
  sig_list <- log(sig_values)
  theta.MCMC<- c(params, sig_list[1])
  which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
  
  
  ## Define the Prior distributions: either normal or log normal distribution
  ## normal distribution
  # sig_mean_population = 0.3
  
  # Set up 4 workers for parallel chains
  n_chains <- 4
  cl <- makeCluster(n_chains)
  registerDoParallel(cl)
  
  # MCMC options
  niter <- 200000
  burnin<- 100000
  #niter <- 2000
  #burnin<-  1000
  outputlen <- niter - burnin
  
  start_time <- Sys.time()
  
  system.time(
    MCMC_chains <- foreach(i = 1:n_chains, 
                           .packages = c("FME", "magrittr","mrgsolve","truncnorm",
                                         "EnvStats","invgamma","readxl"),
                           .export = c("patient_data"),
                           .errorhandling = "pass") %dopar% {
     mod.human.brain <- mod.human.brain
      modMCMC(
        f             = mcmc.fun,
        p             = theta.MCMC,
        niter         = niter,
        jump          = 0.001,# for 207 patient, the jump is 0.001, but for others, we set it as 0.01
        prior         = Prior,            ## prior function
        updatecov     = 100,
        ntrydr        = 1,                ## delayed Rejection
        burninlength  = burnin,
        outputlength  = outputlen,
        verbose=1500
      )
    }
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  stopCluster(cl)   
  
  #MCMC_chains
  saveRDS(MCMC_chains,
          file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
                           paste0("human.MCMC_", pid, ".rds")))
  
  # Combine into one matrix
  all_samples_matrix <- do.call(rbind, lapply(MCMC_chains, function(chain) {
    chain$pars
  }))
  
  
  parar_all_iters = exp(all_samples_matrix) # needs to exp form of parameters
  
  
  ## posterior distributions for parameters
  post_pars  <- parar_all_iters %>% apply(2,mean)
  
  
  which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
  
  Newtime.r   = pred.human.iv(theta.MCMC[-which_sig],DOSE_human)$Time  
  
  # this is the new time variable, now it has been changed to sample per day.
  nrwo.r = length (Newtime.r)
  outputlength = nrow(all_samples_matrix)
  
  # Define the matrix names
  matrix_names <- c("C_Blood","C_tumor_E")
  
  # Initialize lists to store matrices and plots
  MC_matrices <- list()
  MC_plots <- list()
  
  # Loop over matrix names
  for (name in matrix_names) {
    # Create the matrix
    MC_matrices[[name]] <- matrix(nrow = nrwo.r, ncol = outputlength / 1000)
    
    # Loop over columns
    for (i in 1:(outputlength / 1000)) {
      j <- i * 1000
      pars.mouse <- all_samples_matrix[j, ]
      MCdata <- pred.human.iv(pars.mouse,DOSE_human)
      MC_matrices[[name]][, i] <- MCdata[[name]]
    }
    
    # Calculate mean and standard deviation
    M <- apply(MC_matrices[[name]], 1, mean)
    SD <- apply(MC_matrices[[name]], 1, sd)
    
    # Create data for plot
    MC_plot <- cbind(
      Time = Newtime.r,
      as.data.frame(t(apply(MC_matrices[[name]], 1, function(y_est) c(
        median_est = median(y_est, na.rm = TRUE),
        ci_q1 = quantile(y_est, probs = 0.25, names = FALSE, na.rm = TRUE),
        ci_q3 = quantile(y_est, probs = 0.75, names = FALSE, na.rm = TRUE),
        ci_10 = quantile(y_est, probs = 0.1, names = FALSE, na.rm = TRUE),
        ci_90 = quantile(y_est, probs = 0.9, names = FALSE, na.rm = TRUE),
        ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE, na.rm = TRUE),
        ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE, na.rm = TRUE)
      ))))
    )
    
    # Store plot in list
    MC_plots[[name]] <- MC_plot
  }
  Obs.df <- subset(patient_data, Time != 0)
  # Create plots for each parameter
  p.r.L <- gen_MCMC_plot(Obs.df,MC_plots$C_Blood, "Blood", "C_Blood")
  p.r.L
  
  
  ggsave(file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
                   paste0("human.MCMC_", pid, ".pdf")), plot = p.r.L, width = 6, height = 4, units = "in")
}

# Initialize a list to store results
rhat_results <- list()
for (pid in unique_patients) {
  MCMC_chains = readRDS(file = file.path("mcmc_results", "per_patient", 
                                         paste0("patient_", pid), 
                                         paste0("human.MCMC_", pid, ".rds")))
  ## Performance four chains to check the convergences
  MC.mouse.1 = as.mcmc (MCMC_chains[[1]]$pars) # first  chain
  MC.mouse.2 = as.mcmc (MCMC_chains[[2]]$pars) # second chain
  MC.mouse.3 = as.mcmc (MCMC_chains[[3]]$pars) # third  chain
  MC.mouse.4 = as.mcmc (MCMC_chains[[4]]$pars) # fourth chain
  ## combine all chains
  combinedchains = mcmc.list(MC.mouse.1,MC.mouse.2,MC.mouse.3) 
  
  gd <- gelman.diag(combinedchains)

  
  # Store per-parameter point estimates and upper CIs in a data frame
  rhat_df <- data.frame(
    patient   = pid,
    parameter = rownames(gd$psrf),
    psrf = gd$psrf[, "Point est."],
    psrf_upper  = gd$psrf[, "Upper C.I."],
    mpsrf = gd$mpsrf,
    row.names = NULL
  )
  
  rhat_results[[as.character(pid)]] <- rhat_df
  print(pid)
  print(gelman.diag(combinedchains, multivariate = FALSE))
}  
# Combine into a single data frame
rhat_all <- do.call(rbind, rhat_results)

# Save to disk
saveRDS(rhat_all, file = "mcmc_results/rhat_results.rds")
write.csv(rhat_all, file = "mcmc_results/rhat_results.csv", row.names = FALSE)


for (pid in unique_patients) {
  MCMC_chains = readRDS(file = file.path("mcmc_results", "per_patient", 
                                         paste0("patient_", pid), 
                                         paste0("human.MCMC_", pid, ".rds")))
  ## Performance four chains to check the convergences
  MC.mouse.1 = as.mcmc (MCMC_chains[[1]]$pars) # first  chain
  MC.mouse.2 = as.mcmc (MCMC_chains[[2]]$pars) # second chain
  MC.mouse.3 = as.mcmc (MCMC_chains[[3]]$pars) # third  chain
  MC.mouse.4 = as.mcmc (MCMC_chains[[4]]$pars) # fourth chain
  ## combine all chains
  combinedchains = mcmc.list(MC.mouse.1,MC.mouse.2,MC.mouse.3,MC.mouse.4) 
  
  trace_plot <- mcmc_trace (
    combinedchains,
    pars = names(theta.MCMC),
    size = 0.5,
    facet_args = list(nrow = 4)) +
    ggplot2::scale_color_brewer() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))  # Adjust angle of x-axis labels
  
  ggsave(file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
                   paste0("human.MCMC_trace_", pid, ".pdf")), plot = trace_plot, width = 12, height = 8, units = "in")
}  

trace_plot
  plot_prob_chains <- mcmc_dens_overlay(
    combinedchains,
    pars = names(theta.MCMC),
    facet_args = list(nrow=5))  +
    #ggplot2::scale_color_brewer() +
    theme(axis.text.x = element_text(angle = 30))
