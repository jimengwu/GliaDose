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
source("0.human_PBPK_model_4brain.R")
tinterval = 1
TDoses = 1
DOSE_human = 5e8
V_Blood_human  = 5300   #:  ml
V_Tumor_human = 26.36 #ml, got from paper https://www.liebertpub.com/doi/10.1089/hum.2018.001#sec-3 33.5 ml 

params.init.human <- log(c(
  J_Lung =  32.49,
  J_Spleen  =  0.0016,
  J_Liver  =  0.38,
  J_Kidney  = 10.52,
  J_Tumor=2.957,
  J_Brain  =  2.385, 
  J_GI  =  0.000147, 
  J_Other  = 0.0000000484,  
  J_CSFBB = 1,
  J_BCSFB = 1,
  Kdep_Liver  =88.9, #The 1 st oder of depletion of CAR-T cells from the liver
  kon_pop = 1.08e-12,
  koff_pop = 6.48,
  Kg_Exp_Tumor_pop = 0.00291,
  Kkill_max_pop= 0.0445 ,
  IC50_pop =189.1,
  Kexp_CART_max_pop = 0.0604,
  EC50_pop = 0.000404,
  tau_kill_pop = 2.446,
  gamma_act_pop = 0.023,
  gamma_kill_pop = 0.99,
  psi_pop = 0.597,
  Kg_linear_pop = 101359
))



# ============================ for MCMC ==========================

mcmc.fun <- function (pars){
  
  
  out <- pred.human.iv(pars[-which_sig],DOSE_human)
  out = out[which(out$Time %in% patient_data$Time),]
  
  # making sure the predicted value and observed value matching with each other as
  # the time step in predicting function might not be small enough to cover the time step in the 
  # observed time step
  col_order <- names(patient_data)
  out = out[col_order]
  
  ## log-transformed prediction
  cols_to_log <- setdiff(names(out), "Time")
  out[cols_to_log] <- lapply(out[cols_to_log], function(x) ifelse(x != 0, log(x), NA))
  log.yhat <- unname(unlist(out[cols_to_log]))
  
  
  ## log-transformed experimental data
  patient_data[cols_to_log] <- lapply(patient_data[cols_to_log], function(x) ifelse(x != 0, log(x), NA))
  log.y <- unname(unlist(patient_data[cols_to_log]))
  
  non_nan_indices <- which(!is.na(log.y))
  log.yhat        <- log.yhat[non_nan_indices]
  log.y           <- log.y[non_nan_indices]
  
  sig2            <- as.numeric((exp(pars[which_sig][1])))
  
  log_likelihood  <- -2*sum((dnorm (log.y,
                                    mean = log.yhat,
                                    sd   = sqrt(sig2), 
                                    log=TRUE)))
  if (is.na(log_likelihood)){
    return(1e10)
  }else{
    return(log_likelihood)
  }
}
Prior <- function(pars) {
  
  ## Population level
  # The likelihood for population mean (parameters)
  # Calculate likelihoods of each parameters; P(u|mean,CV)
  pars.data = exp(pars[-which_sig]) # parameter value
  
  mean           = exp(theta.MCMC[-which_sig]) # parameter mean in the distribution
  
  CV             = 0.3  
  
  prior_pars     = dtruncnorm(pars.data, 
                              a = qnorm(0.025, mean = mean, sd = mean*CV), 
                              b = qnorm(0.975, mean = mean, sd = mean*CV), 
                              mean = mean, sd = mean*CV ) 
  

  # model residuals
  sig2 <- as.numeric (exp(pars[which_sig][1]))                    # error variances from model residual
  prior_sig2     = dunif (sig2, min = 0.01, max = 0.3)# error variances, Lower and upper boundary from Chiu et al., 2009; Chiu et al., 2014)   
  
  # log-transformed (log-likelihoods of each parameters)
  log.pri.pars   = log (prior_pars)
  log.pri.sig2   = log (prior_sig2)
  

  # maximum likelihood estimation (MLE): negative log-likelihood function, (-2 times sum of log-likelihoods)

  MLE = -2*sum(log.pri.pars,log.pri.sig2)
  return(MLE)
}

pred.human.iv <- function(pars, DOSE_human) {
  
  ## Get out of log domain
  pars %<>% lapply(exp) # important to have because we cannot have nagetive kinetic value
  
  ## Define the exposure scenario
  
  tinterval = 0.5
  TDoses = 1
  
  
  ex.iv_a.human<- ev(ID=1, amt= DOSE_human/V_Blood_human,                  ## Set up the exposure events
                     ii=tinterval, addl=TDoses-1, time = 0,
                     cmt="C_Blood", replicate = FALSE) 
  ex.iv_b.human<- ev(ID=1, amt= V_Tumor_human,                  ## Set up the exposure events
                     ii=tinterval, addl=TDoses-1,  time = 0,
                     cmt="V_Tumor", replicate = FALSE) 
  ex.iv_c.human<- ev(ID=1, amt= 1e-10,                  ## Set up the exposure events
                     ii=tinterval, addl=TDoses-1,  time = 0,
                     cmt="A_E_Tumor", replicate = FALSE) 
  
  ex.iv.human <- c(ex.iv_a.human, ex.iv_b.human,ex.iv_c.human)
  
  
  ## calculate the deposition volume
  results<- 
    mod.human.brain %>% 
    param(pars) %>%
    mrgsim_d(data = ex.iv.human, tgrid=tsamp_human)
  
  out_df <- data.frame(Time=results$time, 
                       C_Blood=results$OutputCART,
                       C_brain_V = results$OutputCART_brain_vascular,
                       V_tumor = results$OutputVolume,
                       C_tumor_tot = results$OutputCART_tumor,
                       C_tumor_E = results$OutputCART_tumor_ex,
                       C_tumor_V = results$OutputCART_tumor_vascular)
  return(out_df)
}




MCcost.human<-function (pars, obs){
  out<- pred.human.iv(pars,DOSE_human)
  cost<- modCost(model=out,obs=obs,weight='mean',x="Time")
  return(cost)
}


gen_MCMC_plot <- function(Obs.df,data, title, ylabel) {
  
  ggplot(data) + 
    geom_ribbon(aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
                fill = "lightblue", alpha = 0.3) +
    geom_ribbon(aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
                fill = "mediumblue", alpha = 0.3) +
    geom_line(aes(x = Time, y = median_est, color = "Median"), size = 1) +
    geom_point(data = Obs.df, aes(x = Time, y = !!sym(ylabel), color = "Observed"), 
               shape = 1, fill = "white", size = 2, stroke = 2) +
    scale_color_manual(values = c("95% CI" = alpha("lightblue", 0.3), 
                                  "50% CI" = alpha("mediumblue", 0.3), 
                                  "Median" = "blue", 
                                  "Observed" = "red"),
                       labels = c("95% CI", "50% CI", "Median", "Observed")) + 
    labs(y = NULL, x = "Time (h)") +
    theme_minimal() +
    theme(
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2,"cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
    xlim(0, round(max(data$Time),1)) +
    
    scale_y_continuous(labels = ~ sprintf(fmt = "%0.1e", .))+
    # Add the title inside the plot
    annotate("text", x = 0.9 * max(data$Time), y = 0.9 * max(data$ci_upper_est, na.rm = TRUE), 
             label = title, 
             size = 4, fontface = "bold", color = "black", hjust = 0.5)
}
