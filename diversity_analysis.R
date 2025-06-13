# Diversity Analysis, Liliane Bader, 15.Juni 2025

# Renal and Bladder Datasets ###################################################
source("config.R")
renal_data <- readRDS(renal_data_path)
bladder_data <- readRDS(bladder_data_path)

# Libraries ####################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(vegan)
library(ineq)
library(patchwork)
library(ggpubr)
################################################################################

# Shannon Diversity ############################################################
compute_shannon_diversity <- function(data, sample_ids) {
  chains <- c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")
  
  diversity_scores <- map_dbl(seq_along(sample_ids), function(i) {
    combined_counts <- map_dfr(chains, function(chain) {
      df <- data[[chain]][[i]]
      if (!is.null(df) && "readCount" %in% names(df)) {
        df["readCount"]
      } else {
        NULL
      }
    })
    
    if (nrow(combined_counts) > 0) {
      diversity(combined_counts$readCount)
    } else {
      NA
    }
  })
  
  tibble(SampleID = sample_ids, Diversity = diversity_scores)
}

renal_shannon <- compute_shannon_diversity(renal_data, renal_data$SampleId)
bladder_shannon <- compute_shannon_diversity(bladder_data, bladder_data$sample)

# Plot
plot_shannon_histogram <- function(diversity_df, title = NULL) {
  ggplot(diversity_df, aes(x = Diversity)) +
    geom_histogram(binwidth = 0.1, fill = "blue", color = "white", alpha = 0.8) +
    theme_minimal() +
    labs(x = "Shannon Diversity", y = "Frequency", title = title) +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 16, face = "bold")
    )
}
renal_shannon_plot <- plot_shannon_histogram(renal_shannon, "Renal: Shannon Diversity")
bladder_shannon_plot <- plot_shannon_histogram(bladder_shannon, "Bladder: Shannon Diversity")

print(renal_shannon_plot)
print(bladder_shannon_plot)

# Save the plot as .pdf
#ggsave("renal_shannon_diversity.pdf", plot = renal_shannon_plot, width = 6, height = 5, units = "in")
#ggsave("bladder_shannon_diversity.pdf", plot = bladder_shannon_plot, width = 6, height = 5, units = "in")
################################################################################

# Normalized Shannon Diversity #################################################
compute_normalized_shannon <- function(read_counts) {
  if (length(read_counts) == 0) return(NA)
  H <- diversity(read_counts)
  S <- length(read_counts)
  
  if (S > 1) {
    H / log(S)
  } else {
    0
  }
}

compute_normalized_shannon_all_chains <- function(data, sample_ids) {
  chains <- c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")
  
  normalized_scores <- purrr::map_dbl(seq_along(sample_ids), function(i) {
    combined_counts <- purrr::map_dfr(chains, function(chain) {
      df <- data[[chain]][[i]]
      if (!is.null(df) && "readCount" %in% names(df)) {
        df["readCount"]
      } else {
        NULL
      }
    })
    compute_normalized_shannon(combined_counts$readCount)
  })
  
  tibble(SampleID = sample_ids, NormalizedShannon = normalized_scores)
}
renal_norm_shannon <- compute_normalized_shannon_all_chains(renal_data, renal_data$SampleId)
bladder_norm_shannon <- compute_normalized_shannon_all_chains(bladder_data, bladder_data$sample)

# Plot
plot_normalized_shannon <- function(data, title = NULL) {
  ggplot(data, aes(x = NormalizedShannon)) +
    geom_histogram(binwidth = 0.01, fill = "blue", color = "white", alpha = 0.8) +
    theme_minimal() +
    labs(x = "Normalized Shannon Diversity", 
         y = "Frequency", 
         #title = title
         ) +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 16, face = "bold")
    )
}
renal_plot <- plot_normalized_shannon(renal_norm_shannon, "Renal: Normalized Shannon")
bladder_plot <- plot_normalized_shannon(bladder_norm_shannon, "Bladder: Normalized Shannon")

print(renal_plot)
print(bladder_plot)

# Save the plot as .pdf
#ggsave("renal_normalized_shannon_diversity.pdf", plot = renal_plot, width = 6, height = 5, units = "in")
#ggsave("bladder_normalized_shannon_diversity.pdf", plot = bladder_plot, width = 6, height = 5, units = "in")
################################################################################

# Gini Index ###################################################################
compute_gini_per_sample <- function(data, sample_ids) {
  chains <- c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")
  
  gini_scores <- map_dbl(seq_along(sample_ids), function(i) {
    combined_counts <- map_dfr(chains, function(chain) {
      df <- data[[chain]][[i]]
      if (!is.null(df) && "readCount" %in% names(df)) {
        df["readCount"]
      } else {
        NULL
      }
    })
    
    if (nrow(combined_counts) > 0) {
      ineq::ineq(combined_counts$readCount, type = "Gini")
    } else {
      NA
    }
  })
  
  tibble(SampleID = sample_ids, Gini = gini_scores)
}

renal_gini <- compute_gini_per_sample(renal_data, renal_data$SampleId)
bladder_gini <- compute_gini_per_sample(bladder_data, bladder_data$sample)

plot_gini_histogram <- function(gini_df, title = NULL) {
  ggplot(gini_df, aes(x = Gini)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white", alpha = 0.8) +
    theme_minimal() +
    labs(x = "Gini Index", y = "Frequency", title = title) +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 16, face = "bold")
    )
}
renal_gini_plot <- plot_gini_histogram(renal_gini, "Renal: Gini Index")
bladder_gini_plot <- plot_gini_histogram(bladder_gini, "Bladder: Gini Index")

print(renal_gini_plot)
print(bladder_gini_plot)

# Save the plot as .pdf
#ggsave("renal_gini_index_diversity.pdf", plot = renal_gini_plot, width = 6, height = 5, units = "in")
#ggsave("bladder_gini_index_diversity.pdf", plot = bladder_gini_plot, width = 6, height = 5, units = "in")


################################################################################

# Clinical Features ############################################################
plot_continuous_clinical <- function(data, xvar, yvar, xlab, ylim_range, label_pos = c(Inf, Inf)) {
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.7, color = "blue") +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    stat_cor(method = "spearman", 
             label.x = label_pos[1], 
             label.y = label_pos[2], 
             size = 5) +
    coord_cartesian(ylim = ylim_range) +
    theme_minimal() +
    labs(x = xlab, y = "Normalized Shannon Diversity") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
}

plot_categorical_clinical <- function(data, xvar, yvar, xlab, test = "t.test", ylim_range) {
  ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_boxplot(alpha = 0.5, fill = "royalblue") +
    geom_jitter(width = 0.2, alpha = 0.7, color = "blue") +
    stat_compare_means(method = test, size = 5, label.y = max(ylim_range, na.rm = TRUE)) +
    coord_cartesian(ylim = ylim_range) +
    theme_minimal() +
    labs(x = xlab, y = "Normalized Shannon Diversity") +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 15))
}

Y <- "NormalizedShannon"

# Renal Data  ##################################################################
# Prepare Data 
renal_data_prepped <- renal_data %>%
  select(SampleId, 
         AgeSample, 
         Sex, 
         Pop1_DSkemaCReaktivtProtein, 
         DonorStatus,
         CharlsonIndex,
         Pop1_Fuhrgrad,
         Pop1_DSkemaHaemoglobin,
         Pop1_DSkemaSLactatdehydrogenase,
         Pop1_DSkemaSeNatrium,
         Leibscore,
         Pop1_DSkemaNeutrofileGranulocyt) %>%
  mutate(
    ShannonDiversiy = renal_shannon$Diversity,
    NormalizedShannon = renal_norm_shannon$NormalizedShannon,
    GiniIndex = renal_gini$Gini,
    Sex = as.factor(case_when(
      Sex == 2 ~ "female",   # Female
      Sex == 1 ~ "male"  # Male
    )),
    
    DonorStatus = factor(case_when(
      DonorStatus == "Healthy" ~ 0,  
      DonorStatus == "Cancer" ~ 1  
    ), levels = c(0, 1)),
    
    CharlsonIndex = as.factor(CharlsonIndex),
    
    Leibscore = as.factor(Leibscore),
    Pop1_DSkemaCReaktivtProtein = as.numeric(Pop1_DSkemaCReaktivtProtein),
    Pop1_DSkemaHaemoglobin = as.numeric(Pop1_DSkemaHaemoglobin),
    Pop1_DSkemaSLactatdehydrogenase = as.numeric(Pop1_DSkemaSLactatdehydrogenase),
    Pop1_DSkemaSeNatrium = as.numeric(Pop1_DSkemaSeNatrium),
    Pop1_DSkemaNeutrofileGranulocyt = as.numeric(Pop1_DSkemaNeutrofileGranulocyt)
  )

# PLot
ylim_range <- range(renal_data_prepped[[Y]], na.rm = TRUE)

age_plot <- plot_continuous_clinical(renal_data_prepped, "AgeSample", Y, "Age", ylim_range, c(30, 0.2))
sex_plot <- plot_categorical_clinical(renal_data_prepped, "Sex", Y, "Sex", "t.test", ylim_range)
charlson_plot <- plot_categorical_clinical(renal_data_prepped, "CharlsonIndex", Y, "Charlson Index", "kruskal.test", ylim_range)
leibscore_plot <- plot_categorical_clinical(renal_data_prepped, "Leibscore", Y, "Leibscore", "kruskal.test", ylim_range)
creaktivtprotein_plot <- plot_continuous_clinical(renal_data_prepped, "Pop1_DSkemaCReaktivtProtein", Y, "C-Reactive Protein", ylim_range, c(100, 0.2))
haemoglobin_plot <- plot_continuous_clinical(renal_data_prepped, "Pop1_DSkemaHaemoglobin", Y, "Haemoglobin", ylim_range, c(Inf, 0.2))
lactatdehydrogenase_plot <- plot_continuous_clinical(renal_data_prepped, "Pop1_DSkemaSLactatdehydrogenase", Y, "LDH", ylim_range, c(300, 0.2))
natrium_plot <- plot_continuous_clinical(renal_data_prepped, "Pop1_DSkemaSeNatrium", Y, "Se-Natrium", ylim_range, c(Inf, 0.2))
neutrofileGranulocyt_plot <- plot_continuous_clinical(renal_data_prepped, "Pop1_DSkemaNeutrofileGranulocyt", Y, "Neutrophil Granulocytes", ylim_range, c(10, 0.2))

renal_clin_vars_all <- (
  (age_plot + sex_plot + charlson_plot) /
    (leibscore_plot + creaktivtprotein_plot + haemoglobin_plot) /
    (lactatdehydrogenase_plot + natrium_plot + neutrofileGranulocyt_plot)
) + 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 45),
    axis.title = element_text(size = 45)
  )
#print(renal_clin_vars_all)
# ggsave("renal_clinical_shannon_grid.pdf", renal_clin_vars_all, width = 40, height = 30)


### Bladder Data ###############################################################
# Prepare Data
bladder_data_prepped <- bladder_data %>%
  select(sample, Age, Sex, BMI, Smoking_status, Death, Clinical_relapse,
         RFS_event, postCX_ctDNA_status, TURBT_Tstage) %>% 
  mutate(
    ShannonDiversiy = bladder_shannon$Diversity,
    NormalizedShannon = bladder_norm_shannon$NormalizedShannon,
    GiniIndex = bladder_gini$Gini,
    
    # Recode and factor variables
    Sex = factor(Sex, levels = c("female", "male")),
    Smoking_status = factor(Smoking_status, levels = c("never", "former", "current")),
    Death = factor(Death, levels = c("No", "Yes")),
    Clinical_relapse = factor(Clinical_relapse, levels = c("No", "Yes")),
    RFS_event = factor(RFS_event, levels = c(0, 1), labels = c("No event", "Event")),
    
    # Recode ctDNA status to shorter labels and set order
    postCX_ctDNA_status = factor(
      recode(postCX_ctDNA_status,
             "ctDNA negative" = "ctDNA neg.",
             "ctDNA positive" = "ctDNA pos."),
      levels = c("ctDNA neg.", "ctDNA pos.")
    ),
    
    # Ordered T stage
    TURBT_Tstage = factor(TURBT_Tstage, 
                          levels = c("T1b", "T2a", "T2b", "T3a", "T4a"), 
                          ordered = TRUE)
  )

# Plot
ylim_range <- range(bladder_data_prepped[[Y]], na.rm = TRUE)

age_plot <- plot_continuous_clinical(bladder_data_prepped, "Age", Y, "Age", ylim_range, c(30, 0.6))
bmi_plot <- plot_continuous_clinical(bladder_data_prepped, "BMI", Y, "BMI", ylim_range, c(25, 0.6))
sex_plot <- plot_categorical_clinical(bladder_data_prepped, "Sex", Y, "Sex", "t.test", ylim_range)
smoking_plot <- plot_categorical_clinical(bladder_data_prepped, "Smoking_status", Y, "Smoking Status", "kruskal.test", ylim_range)
tumor_stage_plot <- plot_categorical_clinical(bladder_data_prepped, "TURBT_Tstage", Y, "Tumor Stage", "kruskal.test", ylim_range)
bladder_data_prepped_clean <- bladder_data_prepped %>%
  filter(!is.na(Death))
death_plot <- plot_categorical_clinical(bladder_data_prepped_clean, "Death", Y, "Death", "t.test", ylim_range)
bladder_data_prepped_clean <- bladder_data_prepped %>%
  filter(!is.na(Clinical_relapse))
clinRelapse_plot <- plot_categorical_clinical(bladder_data_prepped_clean, "Clinical_relapse", Y, "Clinical Relapse", "t.test", ylim_range)
bladder_data_prepped_clean <- bladder_data_prepped %>%
  filter(!is.na(postCX_ctDNA_status))
ctDNAStatus_plot <- plot_categorical_clinical(bladder_data_prepped_clean, "postCX_ctDNA_status", Y, "Post-Chemo ctDNA", "t.test", ylim_range)

bladder_clin_vars_all <- (
  (age_plot + sex_plot + bmi_plot) /
    (death_plot + smoking_plot + clinRelapse_plot) /
    (ctDNAStatus_plot + tumor_stage_plot + plot_spacer())
) + 
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 45, face = "bold"),
    axis.text = element_text(size = 45),
    axis.title = element_text(size = 45)
  )
#print(bladder_clin_vars_all)
#ggsave("bladder_clinical_shannon_grid.pdf", bladder_clin_vars_all, width = 40, height = 30)


