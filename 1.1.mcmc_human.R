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

# load the human clinical data without patient 216
Obs.df.human = read_excel("/Users/wuji/work/code/brain_PBPK/biodistribution_human.xlsx")
Obs.df.human = data.frame(Obs.df.human)
colnames(Obs.df.human) <- c("Time", "C_Blood")
Obs.df.human = Obs.df.human[2:8,] # exclude the first time point

# use the  human physiological parameters from clinical data (181913)

mod.human.brain<- mcode ("brain", human.4B.code)

tinterval = 1
TDoses = 1
DOSE_human = 5e8
V_Blood_human  = 5300   #:  ml
V_Tumor_human = 26.36 #ml, got from paper https://www.liebertpub.com/doi/10.1089/hum.2018.001#sec-3 33.5 ml 

sig_names <- c("sig2",paste0("sig_", names(params.init.human)))

#-------------- fitting --------
params2fit.human =params.init.human

tsamp_human=tgrid(0,700,1)     ## Simulation time 24*7 hours (180 days)



MCcost.human<-function (pars, obs){
  out<- pred.human.iv(pars,DOSE_human)
  cost<- modCost(model=out,obs=obs,weight='mean',x="Time")
  return(cost)
}




Fit.Result.A1<- modFit(f=MCcost.human, p=params2fit.human, obs=Obs.df.human, method ="L-BFGS-B", 
                       control = nls.lm.control(nprint=1)) #"Nelder-Mead"
res.A1=MCcost.human(Fit.Result.A1$par, obs=Obs.df.human)$residuals$res     ## Check the residual for each time points
sum(res.A1^2) 
saveRDS(Fit.Result.A1$par, file = "mcmc_results/avg/initial_fit_par.rds")

Fitted_output.A1 = pred.human.iv(par=Fit.Result.A1$par,DOSE_human)


#pdf("mcmc_results/patient_avg_plot.pdf", width = 6, height = 5)
#pdf("mcmc_results/patient_avg_plot_wo216.pdf", width = 6, height = 5)
pdf("mcmc_results/avg/inital_fitting_plot_wo216.pdf", width = 6, height = 5)


plot(Obs.df.human$Time, Obs.df.human$C_Blood, type = "p", col = "red", 
     xlab = "Time", ylab = "Concentration in Blood", 
     main = "Observed vs Fitted Blood Concentration",ylim=c(0,1.5e5))
lines(Fitted_output.A1$Time, Fitted_output.A1$C_Blood , col = "green", type = "l")
dev.off()



# ------------- mcmc
params <- Fit.Result.A1$par

#-------------sig value for parameters and errors between obs and pred

# Define parameter names

sig_population = 0.3

# Create a named vector with default values
sig_values <- rep(sig_population, length(sig_names))
names(sig_values) = sig_names

# error with prediction and observation's deviation variance
sig_pred = 0.1
sig_values["sig2"] <- sig_pred

# Create sig_list
sig_list <- log(sig_values)

theta.MCMC<- c(params, sig_list[1])
which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
# Create an environment to track iterations

#--------------Maximum likelihood estimation (MLE) function for MCMC----

## Define the Prior distributions: either normal or log normal distribution
## normal distribution
patient_data = Obs.df.human
# Set up 4 workers for parallel chains
n_chains <- 4
cl <- makeCluster(n_chains)
registerDoParallel(cl)

# MCMC options
niter <- 240000
burnin<-  120000
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
    jump          = 0.1,
    prior         = Prior,            ## prior function
    updatecov     = 100,
    ntrydr        = 1,                ## delayed Rejection
    burninlength  = burnin,
    outputlength  = outputlen,
    verbose=500
  )
}
)

end_time <- Sys.time()
print(end_time - start_time)

stopCluster(cl)   

#MCMC_chains
#saveRDS(MCMC_chains,file ='human.MCMC_wo216.rds')
saveRDS(MCMC_chains,file ='mcmc_results/avg/human.MCMC_wo216.rds')


MCMC_chains = readRDS(file ="mcmc_results/avg/human.MCMC_wo216.rds")
## Performance four chains to check the convergences
MC.mouse.1 = as.mcmc (MCMC_chains[[1]]$pars) # first  chain
MC.mouse.2 = as.mcmc (MCMC_chains[[2]]$pars) # second chain
MC.mouse.3 = as.mcmc (MCMC_chains[[3]]$pars) # third  chain
MC.mouse.4 = as.mcmc (MCMC_chains[[4]]$pars) # fourth chain
## combine all chains
combinedchains = mcmc.list(MC.mouse.1,MC.mouse.2,MC.mouse.3,MC.mouse.4) 

gd<-gelman.diag(combinedchains)

pdf("mcmc_diagnostics.pdf", width = 8, height = 10)
plot(combinedchains)
dev.off()

# Store per-parameter point estimates and upper CIs in a data frame
rhat_df <- data.frame(
  parameter = rownames(gd$psrf),
  point_est = gd$psrf[, "Point est."],
  upper_ci  = gd$psrf[, "Upper C.I."],
  mpsrf = gd$mpsrf,
  row.names = NULL
)



write.csv(rhat_df, file = "mcmc_results/rhat_results_avg.csv", row.names = FALSE)


trace_plot <- mcmc_trace (
  combinedchains,
  pars = names(theta.MCMC),
  size = 0.5,
  facet_args = list(nrow = 5)) +
  ggplot2::scale_color_brewer() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))  # Adjust angle of x-axis labels
trace_plot

plot_prob_chains <- mcmc_dens_overlay(
  combinedchains,
  pars = names(theta.MCMC),
  facet_args = list(nrow=5))  +
  #ggplot2::scale_color_brewer() +
  theme(axis.text.x = element_text(angle = 30))

ggsave(file.path("mcmc_results", "avg","human.MCMC_trace.pdf"), 
       plot = trace_plot, width = 12, height = 8, units = "in")
