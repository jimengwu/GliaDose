source("0.human_PBPK_model_4brain.R")
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


# Define doses to test
dose_list <- c(1, 10, 100, 1000, 1e5, 1e6, 1e7,1e8,5e8,1e9)

# Initialize results storage
dose_response <- data.frame(DOSE_human = dose_list, Max_median_est = NA)

for (d in seq_along(dose_list)) {
  print(d)
  DOSE_human <- dose_list[d]
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
  # extract max median_est for C_tumor_E
  dose_response$Max_median_est[d] <- max(na.omit(MC_plots$C_tumor_E$median_est))

  p = ggplot(MC_plots$C_Blood, aes(x = Time)) +
    # 95% CI band
    geom_ribbon(aes(ymin = ci_lower_est, ymax = ci_upper_est),
                fill = "blue", alpha = 0.1) +
    # 25–75% IQR band
    geom_ribbon(aes(ymin = ci_q1, ymax = ci_q3),
                fill = "blue", alpha = 0.2) +
    # Median line
    geom_line(aes(y = median_est),
              color = "black", linewidth = 1, linetype = "dashed") +
    labs(
      x = "Time (h)",
      y = "Blood Concentration",
      title = "Blood Concentration"
    ) +
    ylim(0, 8e5) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  #ggsave(filename = paste0("plot/plot_dose_", d, "_",DOSE_human, "_blood.pdf"),
  #       plot = p, width = 7, height = 5)  # size in inches
  
  
  #plot(MC_plots$C_Blood$Time, MC_plots$V_tumor$median_est,
  #     col = "black", lwd = 2, lty = 2,
  #     xlab = "Time (h)", ylab = "Median Tumor Volume",
  #     main = paste("tumor volume -", name),
  #     ylim = c(0, 100000))
  
  p = ggplot(MC_plots$V_tumor, aes(x = Time)) +
    # 95% CI band
    geom_ribbon(aes(ymin = ci_lower_est, ymax = ci_upper_est),
                fill = "blue", alpha = 0.1) +
    # 25–75% IQR band
    geom_ribbon(aes(ymin = ci_q1, ymax = ci_q3),
                fill = "blue", alpha = 0.2) +
    # Median line
    geom_line(aes(y = median_est),
              color = "black", linewidth = 1, linetype = "dashed") +
    labs(
      x = "Time (h)",
      y = "Median Tumor Volume",
      title = "Tumor volume"
    ) +
    ylim(0, 8e4) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  
  #ggsave(filename = paste0("plot/plot_dose_", d, "_",DOSE_human, "_tumorv.pdf"), 
  #       plot = p, width = 7, height = 5)  # size in inches
  
  
  #plot(MC_plots$C_Blood$Time, MC_plots$C_tumor_E$median_est,
  #     col = "black", lwd = 2, lty = 2,
  #     xlab = "Time (h)", ylab = "Median concentration in tumor extravascular",
  #     main = paste("CTE -", name),
  #     ylim = c(0, 5e8))
  
  p = ggplot(MC_plots$C_tumor_E, aes(x = Time)) +
    # 95% CI band
    geom_ribbon(aes(ymin = ci_lower_est, ymax = ci_upper_est),
                fill = "blue", alpha = 0.1) +
    # 25–75% IQR band
    geom_ribbon(aes(ymin = ci_q1, ymax = ci_q3),
                fill = "blue", alpha = 0.2) +
    # Median line
    geom_line(aes(y = median_est),
              color = "black", linewidth = 1, linetype = "dashed") +
    labs(
      x = "Time (h)",
      y = "Median concentration in tumor extravascular",
      title = "CTE"
    ) +
    ylim(0, 1.8e9) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  #dev.off()
  #ggsave(file = paste0("plot/plot_dose_", d, "_",DOSE_human, "_cte.pdf"), 
  #       plot = p, width = 7, height = 5)  # size in inches
  
}


dose_response

avg_max_CTE = data.frame(
  max_CTE = max(na.omit(MC_plots$C_tumor_E$median_est)),
  max_CT = max(na.omit(MC_plots$C_tumor_tot$median_est)),
  max_CBV = max(na.omit(MC_plots$C_brain_V$median_est))
)





