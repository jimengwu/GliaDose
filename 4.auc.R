source("0.human_PBPK_model_4brain_v2.R")
source("helper.R")
library(ggplot2)

tsamp_human=tgrid(0,720,1)     ## Simulation time 24*7 hours (180 days)

mod.human.brain<- mcode ("brain", human.4B.code)

MCMC_chains =  readRDS(file ="mcmc_results/avg/human.MCMC_wo216.rds")

all_samples_matrix <- do.call(rbind, lapply(MCMC_chains, function(chain) {
  chain$pars
}))

parar_all_iters = exp(all_samples_matrix) # needs to exp form of parameters

sig_values <- data.frame(sig2 = 0.1)
# Create sig_list
sig_list <- log(sig_values)
theta.MCMC<- c(params.init.human, sig_list)

## posterior distributions for parameters
post_pars  <- parar_all_iters %>% apply(2,mean)
which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG
Newtime.r   = pred.human.iv(theta.MCMC[-which_sig],DOSE_human)$Time  

# this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)
# Define the matrix names
matrix_names <- c("C_tumor_E","C_tumor_tot","C_brain_V","C_Blood","V_tumor")



DOSE_human <- 5e8
# Initialize lists to store matrices and plots
MC_matrices <- list()
MC_plots <- list()

outputlength = nrow(all_samples_matrix)
# Loop over matrix names
for (name in matrix_names) {
  # Create the matrix
  MC_matrices[[name]] <- matrix(nrow = nrwo.r, ncol = outputlength / 1000)
  
  # Loop over columns
  for (i in 1:(outputlength / 1000)) {
    j <- i * 1000
    pars.mouse <- all_samples_matrix[j, 1:23]
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

# Subset 0–24h
df_auc <- subset(MC_plots$C_tumor_E, Time <= 24)
df_auc <- na.omit(df_auc)
# Trapezoidal AUC calculation
AUC_0_24 <- sum(diff(df_auc$Time) *
                  (head(df_auc$median_est, -1) + tail(df_auc$median_est, -1)) / 2)

AUC_0_24/24


# Subset 0–24h
df_auc <- subset(MC_plots$C_tumor_E, Time <= 240)
df_auc <- na.omit(df_auc)
# Trapezoidal AUC calculation
AUC_0_240 <- sum(diff(df_auc$Time) *
                  (head(df_auc$median_est, -1) + tail(df_auc$median_est, -1)) / 2)

AUC_0_240/240


# Subset 0–240h
df_auc <- subset(MC_plots$C_tumor_tot, Time <= 240)
df_auc <- na.omit(df_auc)
# Trapezoidal AUC calculation
AUC_0_240 <- sum(diff(df_auc$Time) *
                   (head(df_auc$median_est, -1) + tail(df_auc$median_est, -1)) / 2)

AUC_0_240/240


# ================== for every patient =========================

# load the human clinical data per patient 
Obs.df.human = read_excel("/Users/wuji/work/code/brain_PBPK/biodistribution_human_patient.xlsx")
Obs.df.human = data.frame(Obs.df.human)[,c(5,9,10)]
colnames(Obs.df.human) <- c("patient_id","Time", "C_Blood")



# initialize results storage
auc_results <- data.frame(
  pid = character(),
  whole_tumor = numeric(),
  extra_vascular = numeric(),
  stringsAsFactors = FALSE
)


tsamp_human=tgrid(0,max(Obs.df.human$Time),1)     ## Simulation time 24*7 hours (180 days)
unique_patients <- unique(Obs.df.human$patient_id)
for (pid in unique_patients) {
  #for (pid in c(209,211,213,216,217)){  
  MCMC_chains = readRDS(file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
                                         paste0("human.MCMC_", pid, ".rds")))
  #MC.mouse.1 = as.mcmc (MCMC_chains[[2]]$pars) # first  chain
  print("starting prediction results plotting...")
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
  matrix_names <- c("C_Blood","C_tumor_E","C_tumor_tot")
  
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
  
  # Subset 0–24h
  df_auc <- subset(MC_plots$C_tumor_tot, Time <= 240)
  df_auc <- na.omit(df_auc)
  # Trapezoidal AUC calculation
  AUC_tot <- sum(diff(df_auc$Time) *
                     (head(df_auc$median_est, -1) + tail(df_auc$median_est, -1)) / 2)
  


  # Subset 0–24h
  df_auc <- subset(MC_plots$C_tumor_E, Time <= 240)
  df_auc <- na.omit(df_auc)
  # Trapezoidal AUC calculation
  AUC_extra <- sum(diff(df_auc$Time) *
                     (head(df_auc$median_est, -1) + tail(df_auc$median_est, -1)) / 2)
  
  

  # Append to results
  auc_results <- rbind(auc_results, data.frame(
    pid = pid,
    whole_tumor = AUC_tot,
    extra_vascular = AUC_extra,
    stringsAsFactors = FALSE
  ))
  
}


