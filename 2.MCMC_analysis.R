library(ggplot2)
library(grid)  # for unit()
source("0.human_PBPK_model_4brain_v2.R")
source("helper.R")

mod.human.brain<- mcode ("brain", human.4B.code)

# Create a named vector with default values
sig_population = 0.3
sig_names <- c("sig2",paste0("sig_", names(params.init.human)))
sig_values <- rep(sig_population, length(sig_names))
names(sig_values) = sig_names

# error with prediction and observation's deviation variance
sig_pred = 0.1
sig_values["sig2"] <- sig_pred

# Create sig_list
sig_list <- log(sig_values)

theta.MCMC<- c(params.init.human, sig_list)
which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG


plot_obs_pred <- function(df, r_squared_total, NRMSE) {
  ggplot(df, aes(x = obs, y = pred)) +
    geom_point(size = 3) +
    geom_abline(slope = 1, linetype = "dashed") +  # 1:1 line
    scale_shape_manual(values = c(1)) +  # if only one Variable
    labs(
      x = "Observed CART Concentration in Blood (#/ml)",
      y = "Predicted CART Concentration in Blood (#/ml)"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(fill = NA, color = "black", size = 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = c(0.8, 0.2),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_blank(),
      legend.spacing.x = unit(0, "cm")
    ) +
    scale_x_log10(limits = c(1e2, 1e6)) +
    scale_y_log10(limits = c(1e2, 1e6)) +
    annotate(
      "text", x = 1e2, y = 500000,
      label = paste0("Adj.R² = ", round(r_squared_total, 2),
                     "\nNRMSE = ", round(NRMSE, 2)),
      size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold"
    )
}


plot_ratio_vs_pred <- function(df) {
  ggplot(df, aes(x = pred, y = Ratio)) +
    geom_point(size = 3) +
    labs(
      x = "Predicted CART concentration in Blood (#/ml)",
      y = "Ratio of Prediction-to-Observation"
    ) +
    theme_minimal() +
    geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Reference lines
    scale_shape_manual(values = c(1, 2, 3, 4)) +  # Customize shape per organ
    scale_y_continuous(trans = "log10", limits = c(0.01, 100)) + 
    scale_x_continuous(trans = "log10") +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(fill = NA, color = "black", size = 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = c(0.8, 0.2),
      legend.background = element_rect(fill = "transparent"),
      legend.key = element_blank(),
      legend.spacing.x = unit(0, "cm")
    )
}



gen_MCMC_plot_wo_obs <- function(data, title, ylabel) {
  
  ggplot(data) + 
    geom_ribbon(aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est, color = "95% CI"), 
                fill = "lightblue", alpha = 0.3) +
    geom_ribbon(aes(x = Time, ymin = ci_q1, ymax = ci_q3, color = "50% CI"), 
                fill = "mediumblue", alpha = 0.3) +
    geom_line(aes(x = Time, y = median_est, color = "Median"), size = 1) +
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



# =========================== 1. for each patient ===============================
# load the human clinical data per patient 
Obs.df.human = read_excel("/Users/wuji/work/code/brain_PBPK/biodistribution_human_patient.xlsx")
Obs.df.human = data.frame(Obs.df.human)[,c(5,9,10)]
colnames(Obs.df.human) <- c("patient_id","Time", "C_Blood")


all_obs_pred <- data.frame(
  patient_id = character(),
  Time = numeric(),
  obs = numeric(),
  pred = numeric(),
  pred_CTE = numeric(),
  pred_CTE_lower = numeric(),
  pred_CTE_upper = numeric(),
  pred_CT = numeric(),
  pred_CT_lower = numeric(),
  pred_CT_upper = numeric(),
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
  
  
  patient_data <- subset(Obs.df.human, patient_id == pid)[,2:3]
  
  Obs.df <- subset(patient_data, Time != 0)

  
  
  # Create plots for each parameter
  #p.r.L <- gen_MCMC_plot(Obs.df, MC_plots$C_Blood, "Blood", "C_Blood")
  #p.r.L
  # Modify the first two plots to remove the x-axis
  #p.r.L_no_x <- p.r.L + theme(axis.title.x = element_blank())
  
  #ggsave(paste0("mcmc_results/human.MCMC_", pid, ".pdf"), plot = p.r.L, width = 6, height = 4, units = "in")
  
  # Create plots for each parameter
  p.r.CTE <- gen_MCMC_plot_wo_obs(MC_plots$C_tumor_E, "Tumor_E", "C_tumor_E")


  #ggsave(file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
  #                 paste0("human.MCMC_TE_", pid, ".pdf")),plot = p.r.CTE, width = 6, height = 4, units = "in")
  
  p.r.CT <- gen_MCMC_plot_wo_obs(MC_plots$C_tumor_tot, "Tumor", "C_tumor")

  
  #ggsave(file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
  #                 paste0("human.MCMC_T_", pid, ".pdf")),plot = p.r.CT, width = 6, height = 4, units = "in")
  
  
  

  # Interpolate predicted median estimates at obs time points
  interp_pred <- approx(
    x = MC_plot$Time,
    y = MC_plots$C_Blood$median_est,
    xout = Obs.df$Time,
    rule = 2  # extend with flat extrapolation if needed
  )
  

  # Build per-patient data frame
  patient_df <- data.frame(
    patient_id = pid,
    Time = Obs.df$Time,
    obs = Obs.df$C_Blood,
    pred = interp_pred$y,
    pred_CTE = max(na.omit(MC_plots$C_tumor_E$median_est)),
    pred_CTE_lower = max(na.omit(MC_plots$C_tumor_E$ci_lower_est)),
    pred_CTE_upper = max(na.omit(MC_plots$C_tumor_E$ci_upper_est)),
    pred_CT = max(na.omit(MC_plots$C_tumor_tot$median_est)),
    pred_CT_lower = max(na.omit(MC_plots$C_tumor_tot$ci_lower_est)),
    pred_CT_upper = max(na.omit(MC_plots$C_tumor_tot$ci_upper_est))
  )
  
  
  # Add to the full collection
  all_obs_pred <- rbind(all_obs_pred, patient_df)

  
  # -------------- plot the tumor extra and tumor region on the same plot ------
  
  library(ggplot2)
  MC_plots$C_tumor_E$Source <- "Tumor Extravscular"
  MC_plots$C_tumor_tot$Source <- "Tumor region"
  
  df_combined <- rbind(MC_plots$C_tumor_E, MC_plots$C_tumor_tot)
  
  comparison_tumor_plot = ggplot(df_combined) +
    geom_ribbon(aes(x = Time, ymin = ci_q1, ymax = ci_q3, fill = Source), alpha = 0.3) +
    geom_line(aes(x = Time, y = median_est, color = Source), size = 1) +
    labs(x = "Time (h)", y = "CART Concentration (#/ml)", color = "Compartment", fill = "Compartment") +
    theme_minimal() +
    theme(
      axis.ticks = element_line(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5)
    ) +
    xlim(0, round(max(df_combined$Time, na.rm = TRUE), 1)) +
    scale_y_continuous(labels = ~ sprintf(fmt = "%0.1e", .))
  #ggsave(file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
  #                 paste0("human.MCMC_T_vs_TE_", pid, ".pdf")),plot = comparison_tumor_plot, 
  #                 width = 6, height = 4, units = "in")
  
  
}
# ------ 1.1 for draw the accuracy of prediction plot ---
library(ggplot2)
library(patchwork)

# Remove NA rows if any
all_obs_pred <- na.omit(all_obs_pred)
all_obs_pred = all_obs_pred[all_obs_pred$obs != 0, ]


#without 216
all_obs_pred_o = all_obs_pred[all_obs_pred$patient_id != 216,]

# Left plot: per-patient boxplot
p1 <- ggplot(all_obs_pred_o, aes(x = as.factor(patient_id), y = pred_CTE)) +
  #geom_boxplot(fill = "lightgreen") +
  geom_point(size = 2, color = "#5B9279") +
  geom_errorbar(aes(ymin = pred_CTE_lower, ymax = pred_CTE_upper), 
                width = 0.3, color = "#8FCB9B", size = 1) +
  scale_y_log10() +
  labs(
    x = "Patient ID",
    y = "Predicted CART Concentration in Tumor extravascular (#/ml)",
    #title = "Predicted CE by Patient"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(size = 0.3, color = "gray80"),  # vertical lines
    panel.grid.major.y = element_line(size = 0.3, color = "gray80"),  # horizontal lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # box around plot
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

# Right plot: overall distribution, no Y-axis
p2 <- ggplot(all_obs_pred_o, aes(x = "", y = pred_CTE)) +
  geom_boxplot(fill = "lightblue", width = 0.5) +
  scale_y_log10() +
  labs(
    x = "overall",
    y = NULL,
    #title = "Overall"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_blank(),          # remove y-axis text
    axis.ticks.y = element_blank(),         # remove y-axis ticks
    axis.title.y = element_blank(),         # remove y-axis label
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(size = 0.3, color = "gray80"),  # vertical lines
    panel.grid.major.y = element_line(size = 0.3, color = "gray80"),  # horizontal lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # box around plot
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

# Combine both plots
combined_plot <- p1 + p2 + plot_layout(widths = c(3, 1))

# Show the combined plot
print(combined_plot)

#ggsave(paste0("mcmc_results/per_patient/", "max_CE_combined.pdf"), combined_plot, width = 9, height = 7)


# Left plot: per-patient boxplot
p1_CT <- ggplot(all_obs_pred_o, aes(x = as.factor(patient_id), y = pred_CT)) +
  #geom_boxplot(fill = "lightgreen") +
  geom_point(size = 2, color = "#5B9279") +
  geom_errorbar(aes(ymin = pred_CT_lower, ymax = pred_CT_upper), 
                width = 0.3, color = "#8FCB9B", size = 1) +
  scale_y_log10() +
  labs(
    x = "Patient ID",
    y = "Predicted CART Concentration in Tumor region (#/ml)",
    #title = "Predicted CE by Patient"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(size = 0.3, color = "gray80"),  # vertical lines
    panel.grid.major.y = element_line(size = 0.3, color = "gray80"),  # horizontal lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # box around plot
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

# Right plot: overall distribution, no Y-axis
p2_CT <- ggplot(all_obs_pred_o, aes(x = "", y = pred_CT)) +
  geom_boxplot(fill = "lightblue", width = 0.5) +
  scale_y_log10() +
  labs(
    x = "overall",
    y = NULL,
    #title = "Overall"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_blank(),          # remove y-axis text
    axis.ticks.y = element_blank(),         # remove y-axis ticks
    axis.title.y = element_blank(),         # remove y-axis label
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_line(size = 0.3, color = "gray80"),  # vertical lines
    panel.grid.major.y = element_line(size = 0.3, color = "gray80"),  # horizontal lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # box around plot
    axis.ticks = element_line(size = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

# Combine both plots
combined_CT <- p1_CT + p2_CT + plot_layout(widths = c(3, 1))

# Show the combined plot
print(combined_CT)

ggsave(paste0("mcmc_results/per_patient/", "max_CT_combined.pdf"), combined_CT, width = 8, height = 6)




library(dplyr)
# Group by patient_id and compute R² and NRMSE
patient_stats <- all_obs_pred %>%
  group_by(patient_id) %>%
  summarise(
    n = n(),
    r_squared = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
    NRMSE = sqrt(mean((obs - pred)^2)) / mean(obs),
    r_squared_log = 1 - sum((log(obs + 1) - log(pred + 1))^2) / sum((log(obs + 1) - mean(log(obs + 1)))^2),
    MAE = mean(abs(obs - pred)),
    log_RMSE = sqrt(mean((log(obs + 1) - log(pred + 1))^2))
  ) %>%
  ungroup()

# Show the result
print(patient_stats)


# Compute total R²
r_squared_total <- 1 - sum((all_obs_pred$obs - all_obs_pred$pred)^2) /
  sum((all_obs_pred$obs - mean(all_obs_pred$obs))^2)

cat("Total R² across all patients:", round(r_squared_total, 2), "\n")


r_squared_total <- 1 - sum((all_obs_pred_o$obs - all_obs_pred_o$pred)^2) /
  sum((all_obs_pred_o$obs - mean(all_obs_pred_o$obs))^2)

cat("Total R² across all patients:", round(r_squared_total, 2), "\n")


# 3. Compute NRMSE
NRMSE <- sqrt(mean((all_obs_pred$obs - all_obs_pred$pred)^2)) / mean(all_obs_pred$obs)
# 4. Create plot-friendly long-form dataframe
df_long <- all_obs_pred
df_long$Variable <- "C_Blood"  # Add dummy variable if only one type


# 5. Custom ggplot
obs_pred_plot <- ggplot(df_long, aes(x = obs, y = pred, shape = Variable)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, linetype = "dashed") +  # 1:1 line
  scale_shape_manual(values = c(1)) +  # if only one Variable
  labs(
    x = "Observed CART Concentration in Blood (#/ml)",
    y = "Predicted CART Concentration in Blood (#/ml)"
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(fill = NA, color = "black", size = 2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_blank(),
    legend.spacing.x = unit(0, "cm")
  ) +
  scale_x_log10(limits = c(1e2, 1e6)) +
  scale_y_log10(limits = c(1e2, 1e6)) +
  annotate("text", x = 1e2, y = 500000,
           label = paste0("Adj.R² = ", round(r_squared_total, 2),
                          "\nNRMSE = ", round(NRMSE, 2)),
           size = 4.5, color = "black", hjust = 0, vjust = 0, fontface = "bold")
ggsave(paste0("mcmc_results/per_patient/", "obs_pred_combined.pdf"), obs_pred_plot, width = 6, height = 5)

obs_pred_plot

df_long$Ratio <- df_long$pred / df_long$obs



# Plotting
p_combined_ratio_highlighted = 
  ggplot(df_long, aes(x = pred, y = Ratio, shape = Variable)) +
  geom_point(size = 3) +
  labs(
    x = "Predicted CART concentration in Blood (#/ml)",
    y = "Ratio of Prediction-to-Observation") +
  theme_minimal() +
  geom_hline(yintercept = c(0.5, 2), linetype = "dashed", color = "red") +  # Add dashed lines
  scale_shape_manual(values = c(1,2,3,4))  +  # Define different shapes for each organ
  scale_y_continuous(trans = "log10", limits = c(0.01, 100)) + 
  scale_x_continuous(trans = "log10") +
  theme(panel.background = element_rect(fill = "transparent"),  # Set background to transparent
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        axis.ticks = element_line(),
        axis.ticks.length = unit(0.2,"cm"),
        axis.text = element_text(size = 12),  # Increase text size for axis labels
        axis.title = element_text(size = 14, face = "bold"),  # Increase text size and bold axis titles
        legend.text = element_text(size = 10), # Increase text size for legend label
        legend.title = element_blank(),
        legend.position = c(0.8, 0.2),  # Adjust legend position (x, y)
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        legend.key = element_blank(),  # Remove legend key
        legend.spacing.x = unit(0, "cm")  # Adjust horizontal spacing of legend
  )  # Increase thickness of axis lines)  # Remove minor grid lines
p_combined_ratio_highlighted
ggsave(paste0("mcmc_results/per_patient/", "prediction_acc_combined.pdf"), p_combined_ratio_highlighted, width = 6, height = 5)




# =================== 2. FOR AVERAGE HUMAN CLINICAL DATA ========================



MCMC_chains = readRDS(MCMC_chains, file = "human.MCMC.rds")
MCMC_chains =  readRDS(file ="mcmc_results/avg/human.MCMC_wo216.rds")


print("starting prediction results plotting...")
# Combine into one matrix
all_samples_matrix <- do.call(rbind, lapply(MCMC_chains, function(chain) {
  chain$pars
}))

#all_samples_matrix <- MC.mouse.1
parar_all_iters = exp(all_samples_matrix) # needs to exp form of parameters

## posterior distributions for parameters
post_pars  <- parar_all_iters %>% apply(2,mean)


all_samples_matrix%>% apply(2,mean)
which_sig <- grep("sig", names(theta.MCMC)) # THE INDEX OF SIG

Newtime.r   = pred.human.iv(theta.MCMC[-which_sig],DOSE_human)$Time  

# this is the new time variable, now it has been changed to sample per day.
nrwo.r = length (Newtime.r)
outputlength = nrow(all_samples_matrix)
#outputlength = 100000

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

max_CTE_avg_wo216 <- max(na.omit(MC_plots$C_tumor_E$median_est))
max_CT_avg_wo216 <- max(na.omit(MC_plots$C_tumor_tot$median_est))
saveRDS(max_CTE_avg_wo216,file = file.path("mcmc_results", "avg", 
                                   paste0("max_CE_avg_wo216.rds"))) 
saveRDS(max_CTE_avg_wo216,file = file.path("mcmc_results", "avg", 
                                           paste0("max_CT_avg_wo216.rds"))) 

Obs.df.human = read_excel("/Users/wuji/work/code/brain_PBPK/biodistribution_human.xlsx")
Obs.df.human = data.frame(Obs.df.human)
colnames(Obs.df.human) <- c("Time", "C_Blood")
Obs.df.human = Obs.df.human[2:8,] # exclude the first time point

Obs.df = Obs.df.human


gen_MCMC_plot <- function(data, title, ylabel) {
  
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



# Create plots for each parameter
p.r.L <- gen_MCMC_plot(MC_plots$C_Blood, "Blood", "C_Blood")
p.r.L
ggsave(paste0("mcmc_results/avg/human.MCMC_avg_wo216.pdf"), plot = p.r.L, width = 6, height = 4, units = "in")
#ggsave(paste0("mcmc_results/human.MCMC_avg_wo216.pdf"), plot = p.r.L, width = 6, height = 4, units = "in")

p.r.CTE <- gen_MCMC_plot_wo_obs(MC_plots$C_tumor_E, "Tumor_E", "C_tumor_E")
#ggsave(paste0("mcmc_results/avg/human.MCMC_TE_avg_wo216.pdf"), plot = p.r.CTE, width = 6, height = 4, units = "in")


p.r.CT <- gen_MCMC_plot_wo_obs(MC_plots$C_tumor_tot, "Tumor", "C_tumor")
#ggsave(paste0("mcmc_results/avg/human.MCMC_T_avg_wo216.pdf"), plot = p.r.CT, width = 6, height = 4, units = "in")






# Interpolate predicted median estimates at obs time points
interp_pred <- approx(
  x = MC_plot$Time,
  y = MC_plots$C_Blood$median_est,
  xout = Obs.df$Time,
  rule = 2  # extend with flat extrapolation if needed
)$y

patient_df <- data.frame(
  Time = Obs.df$Time,
  obs = Obs.df$C_Blood,
  pred = interp_pred,
  pred_CE = max(na.omit(MC_plots$C_tumor_E$median_est))
)


# Compute total R²
r_squared_total <- 1 - sum((patient_df$obs - patient_df$pred)^2) /
  sum((patient_df$obs - mean(patient_df$obs))^2)

cat("Total R² across all patients:", round(r_squared_total, 2), "\n")

# Compute NRMSE
NRMSE <- sqrt(mean((patient_df$obs - patient_df$pred)^2)) / mean(patient_df$obs)
patient_df$Ratio <- patient_df$pred / patient_df$obs

obs_pred_plot_human_avg = plot_obs_pred(patient_df,r_squared_total,NRMSE)
p_ratio_human_avg = plot_ratio_vs_pred(patient_df)


#ggsave(paste0("mcmc_results/avg/", "obs_pred_human_avg.pdf"), obs_pred_plot_human_avg, width = 6, height = 5)

#ggsave(paste0("mcmc_results/avg/", "ratio_human_avg.pdf"), p_ratio_human_avg, width = 6, height = 5)

# ============= 3. for the extra vascular concentration calculate the CT & CTE ==============
max_CT_df <- data.frame(
  patient_id = character(),
  max_CTE = numeric(),
  max_CT = numeric(),
  max_CBV = numeric(),
  stringsAsFactors = FALSE
)

for (pid in unique_patients) {
  #for (pid in c(209,211,213,216,217)){  
  #MCMC_chains = readRDS(MCMC_chains, file = paste0("mcmc_results/human.MCMC_", pid, ".rds"))
  MCMC_chains = readRDS(file = file.path("mcmc_results", "per_patient", paste0("patient_", pid), 
                                         paste0("human.MCMC_", pid, ".rds")))
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
  matrix_names <- c("C_tumor_E","C_tumor_tot","C_brain_V")
  
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
  
  # Inside your patient loop (assumes `pid` and `MC_plots` are defined)
  
  # Append as one row
  max_CT_df <- rbind(max_CT_df, data.frame(
    patient_id = pid,
    max_CTE = max(na.omit(MC_plots$C_tumor_E$median_est)),
    max_CT = max(na.omit(MC_plots$C_tumor_tot$median_est)),
    max_CBV = max(na.omit(MC_plots$C_brain_V$median_est))
  ))
} 

saveRDS(max_CT_df,file = file.path("mcmc_results", "per_patient", 
                                   paste0("human.MCMC_max_CE.rds"))) 

readRDS(file = file.path("mcmc_results", "per_patient", 
                                   paste0("human.MCMC_max_CE.rds")))

# ============= 4. for average clinical data calculate the CT & CTE =============

MCMC_chains =  readRDS(file ="mcmc_results/avg/human.MCMC_wo216.rds")

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
matrix_names <- c("C_tumor_E","C_tumor_tot","C_brain_V","C_Blood")

# Initialize lists to store matrices and plots
MC_matrices <- list()
MC_plots <- list()
DOSE_human = 1
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

# Inside your patient loop (assumes `pid` and `MC_plots` are defined)

# Append as one row
avg_max_CTE = data.frame(
  max_CTE = max(na.omit(MC_plots$C_tumor_E$median_est)),
  max_CT = max(na.omit(MC_plots$C_tumor_tot$median_est)),
  max_CBV = max(na.omit(MC_plots$C_brain_V$median_est))
)

avg_max_CTE
# Base plot with log y-axis
plot(MC_plots$C_tumor_E$Time, MC_plots$C_tumor_E$median_est,
     type = "l", col = "black", lwd = 2,
     xlab = "Time (h)", ylab = "Median concentration",
     main = "Concentration–time profiles",
     log = "y")   # log scale for y-axis

# Add other curves
lines(MC_plots$C_brain_V$Time, MC_plots$C_brain_V$median_est,
      col = "blue", lwd = 2, lty = 2)
plot(MC_plots$C_Blood$Time, MC_plots$C_Blood$median_est,
      col = "green", lwd = 2, lty = 2)
lines(MC_plots$C_tumor_tot$Time, MC_plots$C_tumor_tot$median_est,
      col = "red", lwd = 2, lty = 3)

# Legend
legend("bottomright",
       legend = c("Tumor_E", "Brain_V", "Tumor_total"),
       col = c("black", "blue", "red"),
       lty = c(1, 2, 3), lwd = 2, bty = "n", inset = 0.02)
