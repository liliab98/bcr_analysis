# Survival Analysis, Liliane Bader, 15.Juni 2025

# Renal and Bladder Datasets ###################################################
source("config.R")
renal_data <- readRDS(renal_data_path)
bladder_data <- readRDS(bladder_data_path)

# Libraries ####################################################################
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(vegan)
library(survival)
library(survminer)
library(gridExtra)
################################################################################

# Kaplan Meier Curves ##########################################################
# compute Shannon Diversity
compute_total_shannon_per_sample <- function(data, sample_ids, chains = c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")) {
  diversity_scores <- sapply(seq_along(sample_ids), function(i) {
    combined_counts <- do.call(rbind, lapply(chains, function(chain) {
      df <- data[[chain]][[i]]
      if (!is.null(df) && "readCount" %in% colnames(df)) df["readCount"] else NULL
    }))
    if (!is.null(combined_counts) && nrow(combined_counts) > 0) {
      diversity(combined_counts$readCount)
    } else {
      NA
    }
  })
  
  return(data.frame(SampleID = sample_ids, Diversity = diversity_scores))
}
renal_shannon_combined <- compute_total_shannon_per_sample(renal_data, renal_data$SampleId)
bladder_shannon_combined <- compute_total_shannon_per_sample(bladder_data, bladder_data$sample)

# Merge survival time and death event into diversity dataframe
renal_shannon_combined$SurvivalTime <- renal_data$OS_time_year
renal_shannon_combined$Event <- renal_data$Death
bladder_data$Event <- ifelse(bladder_data$Death == "Yes", 1,
                             ifelse(bladder_data$Death == "No", 0, NA))
bladder_shannon_combined$SurvivalTime <- bladder_data$OS_months / 12
bladder_shannon_combined$Event <- bladder_data$Event

# Create a binary group based on median Shannon Index
create_diversity_group <- function(df) {
  median_val <- median(df$Diversity, na.rm = TRUE)
  df$DiversityGroup <- ifelse(df$Diversity > median_val, "High", "Low")
  return(df)
}
renal_shannon_combined <- create_diversity_group(renal_shannon_combined)
bladder_shannon_combined <- create_diversity_group(bladder_shannon_combined)

# Create Survival object
create_surv_object <- function(df) {
  Surv(time = df$SurvivalTime, event = df$Event)
}
surv_obj_renal <- create_surv_object(renal_shannon_combined)
surv_obj_bladder <- create_surv_object(bladder_shannon_combined)


# Fit KM curves
fit_renal <- survfit(surv_obj_renal ~ DiversityGroup, data = renal_shannon_combined)
fit_bladder <- survfit(surv_obj_bladder ~ DiversityGroup, data = bladder_shannon_combined)

# Plot
plot_survival <- function(fit, data, title = NULL) {
  ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    risk.table = TRUE,
    title = title,
    xlab = "Time (Years)",
    ylab = "Survival Probability",
    palette = c("blue", "orange")
  )
}
survival_shannon_renal_plot <- plot_survival(fit_renal, renal_shannon_combined, title = "Renal Dataset")
survival_shannon_bladder_plot <- plot_survival(fit_bladder, bladder_shannon_combined, title = "Bladder Dataset")

print(survival_shannon_renal_plot)
print(survival_shannon_bladder_plot)

# Save the plot as .pdf for the thesis
# Save just the KM survival plot (without risk table)
#ggsave("renal_survival_shannon.pdf", plot = survival_shannon_renal_plot$plot, width = 8, height = 6)
#ggsave("bladder_survival_shannon.pdf", plot = survival_shannon_bladder_plot$plot, width = 8, height = 6)
# Arrange plot and table vertically (plot on top, table below)
#pdf("renal_survival_shannon_full.pdf", width = 8, height = 6)
#grid.arrange(survival_shannon_renal_plot$plot, survival_shannon_renal_plot$table, 
#             ncol = 1, heights = c(2/3, 1/3))
#dev.off()
#pdf("bladder_survival_shannon_full.pdf", width = 8, height = 6)
#grid.arrange(survival_shannon_bladder_plot$plot, survival_shannon_bladder_plot$table, 
#             ncol = 1, heights = c(2/3, 1/3))
#dev.off()

# Save renal survival plot + table as PNG for the defense
# png("renal_survival_shannon_full.png", width = 8, height = 6, units = "in", res = 300)
# grid.arrange(
#   survival_shannon_renal_plot$plot,
#   survival_shannon_renal_plot$table,
#   ncol = 1,
#   heights = c(2/3, 1/3)
# )
# dev.off()
# png("bladder_survival_shannon_full.png", width = 8, height = 6, units = "in", res = 300)
# grid.arrange(
#   survival_shannon_bladder_plot$plot,
#   survival_shannon_bladder_plot$table,
#   ncol = 1,
#   heights = c(2/3, 1/3)
# )
# dev.off()

################################################################################

# Cox Model ####################################################################
# Prepare data
bladder_data$Event <- ifelse(bladder_data$Death == "Yes", 1,
                             ifelse(bladder_data$Death == "No", 0, NA))
bladder_shannon_combined$Age <- bladder_data$Age
bladder_shannon_combined$Stage <- factor(bladder_data$RC_Tstage)

renal_shannon_combined$Event <- ifelse(renal_data$Death == 1, 1, 0)
renal_shannon_combined$SurvivalTime <- renal_data$OS_time_year
renal_shannon_combined$Age <- renal_data$AgeSample
renal_shannon_combined$HealthStatus <- ifelse(renal_data$DonorStatus == "Healthy", 0,
                                              ifelse(renal_data$DonorStatus == "Cancer", 1, NA))
renal_shannon_combined$Charlson <- renal_data$CharlsonIndex  # keep numeric for now

# Define required columns
renal_vars <- c("SurvivalTime", "Event", "Diversity", "Age")
bladder_vars <- c("SurvivalTime", "Event", "Diversity", "Age")

# Clean data
prepare_cox_data <- function(df, required_vars) {
  df[complete.cases(df[, required_vars]), ]
}
renal_cox_data <- prepare_cox_data(renal_shannon_combined, renal_vars)
bladder_cox_data <- prepare_cox_data(bladder_shannon_combined, bladder_vars)

# Fit models
cox_renal <- coxph(Surv(SurvivalTime, Event) ~ Diversity + Age, data = renal_cox_data)
cox_bladder <- coxph(Surv(SurvivalTime, Event) ~ Diversity + Age, data = bladder_cox_data)


# Summaries
summary(cox_renal)
summary(cox_bladder)

# Plots
plot_forest <- function(cox_model, data, title = NULL, fontsize = 1.5) {
  ggforest(cox_model, data = data, main = title, fontsize = fontsize)
}
renal_hazard_ratio_plot <- plot_forest(cox_renal, renal_cox_data, title = "Renal Cox Model: Hazard Ratios")
bladder_hazard_ratio_plot <- plot_forest(cox_bladder, bladder_cox_data, title = "Bladder Cox Model: Hazard Ratios")

print(renal_hazard_ratio_plot)
print(bladder_hazard_ratio_plot)

# Save the plot as .pdf for the thesis
# ggsave("renal_hazard_ratio.pdf", plot = renal_hazard_ratio_plot, width = 10, height = 6)
# ggsave("bladder_hazard_ratio.pdf", plot = bladder_hazard_ratio_plot, width = 10, height = 6)

# Save the plot as .png for defense
ggsave("renal_hazard_ratio.png", plot = renal_hazard_ratio_plot, width = 9, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("bladder_hazard_ratio.png", plot = bladder_hazard_ratio_plot, width = 9, height = 5, units = "in", dpi = 300, bg = "white")

# Test assumptions (proportional hazards):
cox.zph(cox_renal)
cox.zph(cox_bladder)

plot(cox.zph(cox_renal))
plot(cox.zph(cox_bladder))


