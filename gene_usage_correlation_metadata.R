# Gene Usage Correlation with clinical Variables, Liliane Bader, 15.Juni 2025

# Renal and Bladder Datasets ###################################################
source("config.R")
renal_data <- readRDS(renal_data_path)
bladder_data <- readRDS(bladder_data_path)

# Libraries ####################################################################
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)
library(RColorBrewer)
################################################################################

# Add Healthy subset to Bladder dataset ########################################
bladder_data_healthy <- bladder_data
bladder_data_healthy$DonorStatus <- "Cancer"

healthy_donors <- renal_data %>%
  filter(DonorStatus == "Healthy") %>%
  mutate(
    Sex = case_when(
      Sex == 1 ~ "male",  
      Sex == 2 ~ "female"  
    ),
    Death = as.character(Death)
  ) %>%
  rename(
    sample = SampleId,
    Age = AgeSample
  ) %>%
  select(
    sample,
    Sex,
    Age,
    Death,
    DonorStatus,
    mixcr_df_IGH,
    mixcr_df_IGK,
    mixcr_df_IGL
  )

bladder_data_healthy <- bind_rows(bladder_data_healthy, healthy_donors)

summary(as.factor(bladder_data_healthy$DonorStatus))
#summary(subset(bladder_data_healthy, DonorStatus == "Healthy"))

clean_gene_names <- function(df, gene_columns) {
  df %>%
    mutate(across(all_of(gene_columns), ~ gsub("\\*00$", "", .x)))
}

bladder_data_healthy <- bladder_data_healthy %>%
  mutate(
    mixcr_df_IGH = map(mixcr_df_IGH, ~ clean_gene_names(.x, c("V.Name", "D.Name", "J.Name"))),
    mixcr_df_IGK = map(mixcr_df_IGK, ~ clean_gene_names(.x, c("V.Name", "J.Name"))),
    mixcr_df_IGL = map(mixcr_df_IGL, ~ clean_gene_names(.x, c("V.Name", "J.Name")))
  )
################################################################################

# Compute Gene Usage ###########################################################
extract_gene_usage <- function(df_list, sample_names, gene_column) {
  bind_rows(df_list, .id = "sample") %>%
    mutate(
      sample = sample_names[as.numeric(sample)],
    ) %>%
    group_by(sample, !!sym(gene_column)) %>%
    summarize(count = n(), .groups = "drop") %>%
    ungroup()
}

normalize_gene_usage <- function(gene_usage_df) {
  gene_usage_df %>%
    group_by(sample) %>%
    mutate(freq = count / sum(count)) %>%
    ungroup()
}

process_gene_usage <- function(data, sample_colname) {
  
  chain_data <- list(
    IGH = list(df = data$mixcr_df_IGH, genes = c("V.Name", "D.Name", "J.Name")),
    IGK = list(df = data$mixcr_df_IGK, genes = c("V.Name", "J.Name")),
    IGL = list(df = data$mixcr_df_IGL, genes = c("V.Name", "J.Name"))
  )
  
  gene_usage_results <- map(chain_data, function(chain) {
    map(chain$genes, ~ extract_gene_usage(chain$df, data[[sample_colname]], .x)) %>%
      set_names(chain$genes)
  })
  
  normalized_gene_usage <- map(gene_usage_results, function(chain) {
    map(chain, normalize_gene_usage)
  })
  
  return(normalized_gene_usage)
}

# Access flattened results
get_usage <- function(normalized_data, chain, gene) {
  normalized_data[[chain]][[gene]]
}

# Process bladder data
normalized_gene_usage_bladder <- process_gene_usage(bladder_data_healthy, "sample")

v_usage_IGH_bladder <- get_usage(normalized_gene_usage_bladder, "IGH", "V.Name")
d_usage_IGH_bladder <- get_usage(normalized_gene_usage_bladder, "IGH", "D.Name")
j_usage_IGH_bladder <- get_usage(normalized_gene_usage_bladder, "IGH", "J.Name")

v_usage_IGK_bladder <- get_usage(normalized_gene_usage_bladder, "IGK", "V.Name")
j_usage_IGK_bladder <- get_usage(normalized_gene_usage_bladder, "IGK", "J.Name")

v_usage_IGL_bladder <- get_usage(normalized_gene_usage_bladder, "IGL", "V.Name")
j_usage_IGL_bladder <- get_usage(normalized_gene_usage_bladder, "IGL", "J.Name")

# Process renal data
normalized_gene_usage_renal <- process_gene_usage(renal_data, "SampleId")

v_usage_IGH_renal <- get_usage(normalized_gene_usage_renal, "IGH", "V.Name")
d_usage_IGH_renal <- get_usage(normalized_gene_usage_renal, "IGH", "D.Name")
j_usage_IGH_renal <- get_usage(normalized_gene_usage_renal, "IGH", "J.Name")

v_usage_IGK_renal <- get_usage(normalized_gene_usage_renal, "IGK", "V.Name")
j_usage_IGK_renal <- get_usage(normalized_gene_usage_renal, "IGK", "J.Name")

v_usage_IGL_renal <- get_usage(normalized_gene_usage_renal, "IGL", "V.Name")
j_usage_IGL_renal <- get_usage(normalized_gene_usage_renal, "IGL", "J.Name")
################################################################################

# Prepare Clinical Variables ###################################################
# Bladder
metadata_bladder <- bladder_data_healthy %>% 
  select(sample, 
         Age, 
         Sex, 
         BMI,
         Smoking_status,
         Death,
         DonorStatus,
         Clinical_relapse,
         RFS_event,
         postCX_ctDNA_status,
         TURBT_Tstage
  ) %>% 
  mutate(
    # sex: female = 0, male = 1
    Sex = factor(ifelse(Sex == "female", 0, 1), levels = c(0, 1)),
    
    DonorStatus = factor(case_when(
      DonorStatus == "Healthy" ~ 0,  
      DonorStatus == "Cancer" ~ 1  
    ), levels = c(0, 1)),
    
    # smoking_status: never = 0, former = 1, current = 2
    Smoking_status = factor(case_when(
      Smoking_status == "never" ~ 0,
      Smoking_status == "former" ~ 1,
      Smoking_status == "current" ~ 2
    ), levels = c(0, 1, 2)),
    
    # death: No = 0 (better), Yes = 1 (worse)
    Death = factor(ifelse(Death == "No", 0, 1), levels = c(0, 1)),
    
    # clinical_relapse: No = 0, Yes = 1
    Clinical_relapse = factor(ifelse(Clinical_relapse == "No", 0, 1), levels = c(0, 1)),
    
    # RFS_event: 0 = no event (better), 1 = event (worse)
    RFS_event = factor(RFS_event, levels = c(0, 1)),
    
    # ctDNA_status: ctDNA negative = 0, ctDNA positive = 1
    postCX_ctDNA_status = factor(case_when(
      postCX_ctDNA_status == "ctDNA negative" ~ 0,
      postCX_ctDNA_status == "ctDNA positive" ~ 1
    ), levels = c(0, 1)),
    
    # TURBT_Tstage: 
    TURBT_Tstage = factor(TURBT_Tstage, 
                          levels = c("T1b", "T2a", "T2b", "T3a", "T4a"), 
                          ordered = TRUE),
    TURBT_Tstage_grouped = factor(case_when(
      TURBT_Tstage %in% c("T1b") ~ 0,
      TURBT_Tstage %in% c("T2a", "T2b") ~ 1,
      TURBT_Tstage %in% c("T3a", "T4a") ~ 2
    ), levels = c(0, 1, 2))
  )

# All have to be numeric, because Spearman correlation requires numeric inputs (ranks)
metadata_numeric_bladder <- metadata_bladder %>%
  select(sample, Age, BMI, Sex, DonorStatus, Smoking_status, Death, Clinical_relapse,
         RFS_event, postCX_ctDNA_status, TURBT_Tstage_grouped) %>%
  mutate(across(-sample, ~ as.numeric(as.character(.))))

# Renal
metadata_renal <- renal_data %>% 
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
    Sex = as.numeric(case_when(
      Sex == 2 ~ 0,   # Female
      Sex == 1 ~ 1  # Male
    )),
    
    DonorStatus = factor(case_when(
      DonorStatus == "Healthy" ~ 0,  
      DonorStatus == "Cancer" ~ 1  
    ), levels = c(0, 1)),
    
    CharlsonIndex = as.numeric(CharlsonIndex),
    
    Pop1_Fuhrgrad = na_if(Pop1_Fuhrgrad, "Unknown"),
    Pop1_Fuhrgrad = as.numeric(Pop1_Fuhrgrad),
    
    Leibscore = as.numeric(Leibscore)
  ) %>%
  rename(sample = SampleId,
         Age = AgeSample,
         Pop1_CReaktivtProtein = Pop1_DSkemaCReaktivtProtein,
         Pop1_Haemoglobin = Pop1_DSkemaHaemoglobin,
         Pop1_SLactatdehydrogenase = Pop1_DSkemaSLactatdehydrogenase,
         Pop1_SeNatrium = Pop1_DSkemaSeNatrium,
         Pop1_NeutrofileGranulocyt = Pop1_DSkemaNeutrofileGranulocyt)

# All have to be numeric, because Spearman correlation requires numeric inputs (ranks)
metadata_numeric_renal <- metadata_renal %>%
  select(sample, 
         Age, 
         Sex, 
         Pop1_CReaktivtProtein, 
         DonorStatus,
         CharlsonIndex,
         Pop1_Fuhrgrad,
         Pop1_Haemoglobin,
         Pop1_SLactatdehydrogenase,
         Pop1_SeNatrium,
         Leibscore,
         Pop1_NeutrofileGranulocyt) %>%
  mutate(across(-sample, ~ as.numeric(as.character(.))))

################################################################################

# Statistical tests (Spearmen) #################################################
# Spearman rank correlation coefficient. 
# This test assesses how well the relationship between two variables can be described,
# by a monotonic function (i.e., whether they tend to increase or decrease together, 
# but not necessarily in a linear fashion).
calc_spearman_stats <- function(gene_usage_df, clinical_df, gene_column, clinical_var) {
  gene_usage_df %>%
    left_join(clinical_df, by = "sample") %>%  # Use dynamic column for sample
    distinct(sample, !!sym(gene_column), .keep_all = TRUE) %>%
    group_by(!!sym(gene_column)) %>%
    filter(sum(!is.na(freq)) > 2, sum(!is.na(.data[[clinical_var]])) > 2) %>%
    group_modify(~{
      df <- .x
      x <- df$freq
      y <- df[[clinical_var]]
      
      if (length(x) > 2 && length(unique(y)) > 1) {
        test <- suppressWarnings(cor.test(x, y, method = "spearman"))
        tibble(estimate = test$estimate, p_value = test$p.value, stat = test$statistic)
      } else {
        message("Skipped gene ", unique(df[[gene_column]]), " for ", clinical_var, ": Not enough variation.")
        tibble(estimate = NA_real_, p_value = NA_real_, stat = NA_real_)
      }
    }) %>%
    ungroup()
}
# estimate tells you the strength and direction of the monotonic relationship.
# p_value tells you whether the correlation is statistically significant.
# statistic is a measure of the test statistic used to calculate the p-value.

# Bladder 
clinical_vars_bladder <- c("Age", "Sex", "DonorStatus", "BMI", "Death", "Smoking_status", "Clinical_relapse", "postCX_ctDNA_status", "TURBT_Tstage_grouped")

cor_results_bladder <- bind_rows(
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(v_usage_IGH_bladder, metadata_numeric_bladder, "V.Name", .x) %>% mutate(Gene_Type = "IGH-V", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(d_usage_IGH_bladder, metadata_numeric_bladder, "D.Name", .x) %>% mutate(Gene_Type = "IGH-D", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(j_usage_IGH_bladder, metadata_numeric_bladder, "J.Name", .x) %>% mutate(Gene_Type = "IGH-J", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(v_usage_IGK_bladder, metadata_numeric_bladder, "V.Name", .x) %>% mutate(Gene_Type = "IGK-V", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(j_usage_IGK_bladder, metadata_numeric_bladder, "J.Name", .x) %>% mutate(Gene_Type = "IGK-J", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(v_usage_IGL_bladder, metadata_numeric_bladder, "V.Name", .x) %>% mutate(Gene_Type = "IGL-V", Condition = .x)),
  map_df(clinical_vars_bladder, ~ calc_spearman_stats(j_usage_IGL_bladder, metadata_numeric_bladder, "J.Name", .x) %>% mutate(Gene_Type = "IGL-J", Condition = .x))
)

# Renal 
clinical_vars_renal <- c("Age", "Sex", "DonorStatus", "CharlsonIndex", "Leibscore", "Pop1_Fuhrgrad", "Pop1_CReaktivtProtein", "Pop1_Haemoglobin", "Pop1_SLactatdehydrogenase", "Pop1_SeNatrium", "Pop1_NeutrofileGranulocyt")

cor_results_renal <- bind_rows(
  map_df(clinical_vars_renal, ~ calc_spearman_stats(v_usage_IGH_renal, metadata_numeric_renal, "V.Name", .x) %>% mutate(Gene_Type = "IGH-V", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(d_usage_IGH_renal, metadata_numeric_renal, "D.Name", .x) %>% mutate(Gene_Type = "IGH-D", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(j_usage_IGH_renal, metadata_numeric_renal, "J.Name", .x) %>% mutate(Gene_Type = "IGH-J", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(v_usage_IGK_renal, metadata_numeric_renal, "V.Name", .x) %>% mutate(Gene_Type = "IGK-V", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(j_usage_IGK_renal, metadata_numeric_renal, "J.Name", .x) %>% mutate(Gene_Type = "IGK-J", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(v_usage_IGL_renal, metadata_numeric_renal, "V.Name", .x) %>% mutate(Gene_Type = "IGL-V", Condition = .x)),
  map_df(clinical_vars_renal, ~ calc_spearman_stats(j_usage_IGL_renal, metadata_numeric_renal, "J.Name", .x) %>% mutate(Gene_Type = "IGL-J", Condition = .x))
)
################################################################################

# Preprocessing: Renal #########################################################
# Multiple testing correction per clinical feature
# Benjamini-Hochberg (FDR)
# Testing many hypotheses (each gene–clinical variable pair is a separate test).
# The Benjamini-Hochberg (BH) method is less conservative than Bonferroni and 
# is more appropriate when you're willing to tolerate some false positives but 
# want to control the proportion of them.
cor_results_renal_corrected <- cor_results_renal %>%
  group_by(Condition) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Add direction for test statistic through the estimate
# Combine the three Gene Name columns into one
renal_cor <- cor_results_renal_corrected %>%
  mutate(Gene_Name = coalesce(V.Name, D.Name, J.Name), # Choose the first non-NA value
         signed_stat = ifelse(!is.na(estimate),
                              estimate, 
                              NA_real_                     
         ),
         p_val = p_adj) %>%  
  filter(Gene_Name != "") %>% 
  select(Gene_Name, Condition, signed_stat, p_val) 

# Winsorizing or winsorization is the transformation of statistics by 
# limiting extreme values in the statistical data to reduce the effect of 
# possibly spurious outliers. It is named after the engineer-turned-biostatistician 
# Charles P. Winsor (1895–1951). The effect is the same as clipping in signal processing.
# Compute 5th and 95th percentiles of signed_log_p
#p5_renal <- quantile(renal_cor$signed_stat, 0.025, na.rm = TRUE)
#p95_renal <- quantile(renal_cor$signed_stat, 0.975, na.rm = TRUE)

# Clamp values (winsorize)
#renal_cor <- renal_cor %>%
#  mutate(
#    signed_stat = pmin(pmax(signed_stat, p5_renal), p95_renal)
#  )
# Even more normalization to make values more comparable
#renal_cor <- renal_cor %>%
#  mutate(signed_stat = scale(signed_stat))  # Standardizing (z-scoring)

#renal_cor <- renal_cor %>%
#  group_by(Gene_Name) %>%
#  mutate(signed_stat = (signed_stat - mean(signed_stat, na.rm = TRUE)) / #sd(signed_stat, na.rm = TRUE))  # Gene-wise normalization

# Convert to wide format for heatmap
heatmap_data_renal <- renal_cor %>%
  select(-p_val) %>%
  pivot_wider(names_from = Gene_Name, values_from = signed_stat) %>%
  column_to_rownames(var = "Condition") 

# Matrix for each Gene (V, D, J)
heatmap_V_renal <- heatmap_data_renal %>% select(starts_with("IGHV")) %>% as.matrix()
heatmap_D_renal <- heatmap_data_renal %>% select(starts_with("IGHD")) %>% as.matrix()
heatmap_J_renal <- heatmap_data_renal %>% select(starts_with("IGHJ")) %>% as.matrix()

# Indicate significant p_values through stars 
signif_labels_renal <- renal_cor %>%
  mutate(
    signif = case_when(
      p_val <= 0.001 ~ "***",
      p_val <= 0.01  ~ "**",
      p_val <= 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(Condition, Gene_Name, signif) %>%
  pivot_wider(names_from = Gene_Name, values_from = signif, values_fill = "") %>%
  column_to_rownames("Condition")

signif_matrix_renal <- as.matrix(signif_labels_renal)

signif_V_renal <- signif_matrix_renal[rownames(heatmap_V_renal), colnames(heatmap_V_renal)]
signif_D_renal <- signif_matrix_renal[rownames(heatmap_D_renal), colnames(heatmap_D_renal)]
signif_J_renal <- signif_matrix_renal[rownames(heatmap_J_renal), colnames(heatmap_J_renal)]

# Order genes for age
# First, get gene order based on correlation with age
age_correlations_renal <- renal_cor %>%
  filter(Condition == "Age") %>%
  arrange(signed_stat) %>%   # Sort from most negative to most positive
  pull(Gene_Name)

# Then reorder the columns of your heatmap matrices
heatmap_V_renal <- heatmap_V_renal[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_V_renal)]]
heatmap_D_renal <- heatmap_D_renal[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_D_renal)]]
heatmap_J_renal<- heatmap_J_renal[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_J_renal)]]

# Same for the significance matrices
signif_V_renal <- signif_V_renal[, colnames(heatmap_V_renal)]
signif_D_renal <- signif_D_renal[, colnames(heatmap_D_renal)]
signif_J_renal <- signif_J_renal[, colnames(heatmap_J_renal)]
################################################################################

# Preprocessing: Bladder #######################################################
# Multiple testing correction per clinical feature
# Benjamini-Hochberg (FDR)
cor_results_bladder_corrected <- cor_results_bladder %>%
  group_by(Condition) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

# Add direction for test statistic through the estimate
# Combine the three Gene Name columns into one
bladder_cor <- cor_results_bladder_corrected %>%
  mutate(Gene_Name = coalesce(V.Name, D.Name, J.Name), # Choose the first non-NA value
         signed_stat = ifelse(!is.na(estimate),
                              estimate, 
                              NA_real_                     
         ),
         p_val = p_adj) %>%  
  filter(Gene_Name != "") %>% 
  select(Gene_Name, Condition, signed_stat, p_val) 

# Winsorizing 
# Compute 5th and 95th percentiles of signed_log_p
#p5_bladder <- quantile(bladder_cor$signed_stat, 0.025, na.rm = TRUE)
#p95_bladder <- quantile(bladder_cor$signed_stat, 0.975, na.rm = TRUE)

# Clamp values (winsorize)
#bladder_cor <- bladder_cor %>%
#  mutate(
#    signed_stat = pmin(pmax(signed_stat, p5_renal), p95_renal)
#  )
# Even more normalization to make values more comparable
#bladder_cor <- bladder_cor %>%
#  mutate(signed_stat = scale(signed_stat))  # Standardizing (z-scoring)

#bladder_cor <- bladder_cor %>%
#  group_by(Gene_Name) %>%
#  mutate(signed_stat = (signed_stat - mean(signed_stat, na.rm = TRUE)) / #sd(signed_stat, na.rm = TRUE))  # Gene-wise normalization

# Convert to wide format for heatmap
heatmap_data_bladder <- bladder_cor %>%
  select(-p_val) %>%
  pivot_wider(names_from = Gene_Name, values_from = signed_stat) %>%
  column_to_rownames(var = "Condition") 

# Matrix for each Gene (V, D, J)
heatmap_V_bladder <- heatmap_data_bladder %>% select(starts_with("IGHV")) %>% as.matrix()
heatmap_D_bladder <- heatmap_data_bladder %>% select(starts_with("IGHD")) %>% as.matrix()
heatmap_J_bladder <- heatmap_data_bladder %>% select(starts_with("IGHJ")) %>% as.matrix()

# Indicate significant p_values through stars 
signif_labels_bladder <- bladder_cor %>%
  mutate(
    signif = case_when(
      p_val <= 0.001 ~ "***",
      p_val <= 0.01  ~ "**",
      p_val <= 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(Condition, Gene_Name, signif) %>%
  pivot_wider(names_from = Gene_Name, values_from = signif, values_fill = "") %>%
  column_to_rownames("Condition")

signif_matrix_bladder <- as.matrix(signif_labels_bladder)

signif_V_bladder <- signif_matrix_bladder[rownames(heatmap_V_bladder), colnames(heatmap_V_bladder)]
signif_D_bladder <- signif_matrix_bladder[rownames(heatmap_D_bladder), colnames(heatmap_D_bladder)]
signif_J_bladder <- signif_matrix_bladder[rownames(heatmap_J_bladder), colnames(heatmap_J_bladder)]

# Order genes for age 
# First, get gene order based on correlation with age
age_correlations_bladder <- bladder_cor %>%
  filter(Condition == "Age") %>%
  arrange(signed_stat) %>%   # Sort from negative to positive
  pull(Gene_Name)

# Then reorder the columns of your heatmap matrices
heatmap_V_bladder <- heatmap_V_bladder[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_V_bladder)]]
heatmap_D_bladder <- heatmap_D_bladder[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_D_bladder)]]
heatmap_J_bladder <- heatmap_J_bladder[, age_correlations_renal[age_correlations_renal %in% colnames(heatmap_J_bladder)]]

# Same for the significance matrices
signif_V_bladder <- signif_V_bladder[, colnames(heatmap_V_bladder)]
signif_D_bladder <- signif_D_bladder[, colnames(heatmap_D_bladder)]
signif_J_bladder <- signif_J_bladder[, colnames(heatmap_J_bladder)]

# Normalize each row (z-score): subtract mean and divide by sd
#row_norm <- function(heatmap) {
#  t(apply(heatmap, 1, function(row) {
#  if (sd(row, na.rm = TRUE) == 0) {
#    return(rep(0, length(row)))  # Avoid division by zero
#  } else {
#    (row - mean(row, na.rm = TRUE)) / sd(row, na.rm = TRUE)
#  }
#}))
#}
################################################################################

# Heatmap: Bladder #############################################################
# range(bladder_cor$signed_stat, renal_cor$signed_stat, na.rm = TRUE)
my_breaks <- seq(-1, 1, by = 0.025)  # 81 breaks → 80 intervals → use 80 colors
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(my_breaks) - 1)

row_names_bladder <- c(
  "Age: 45 to 80", 
  "Sex: Female (0) / Male (1)", 
  "Health Status: Healthy (0) / Cancer (1)", 
  "BMI: 18 to 49", 
  "Death: No (0) / Yes (1)", 
  "Smoking Status:\nNever (0) / Former (1) / Current (2)", 
  "Clinical Relapse: No (0) / Yes (1)", 
  "Post-CX ctDNA Status: Negative (0) / Positive (1)", 
  "TURBT T-stage Grouped: T1 (0) / T2 (1) / T3+4 (2)"
)
rownames(heatmap_V_bladder) <- row_names_bladder
rownames(heatmap_D_bladder) <- row_names_bladder
rownames(heatmap_J_bladder) <- row_names_bladder

# Product of Spearman test statistic (S) and correlation coefficient (rho))
bladder_V_gene_usage_clinVars_plot <- pheatmap(heatmap_V_bladder,
                                               cluster_rows = FALSE,
                                               cluster_cols = FALSE,
                                               main = "Bladder IgH V Gene Usage",
                                               display_numbers = signif_V_bladder,
                                               fontsize_row = 10,
                                               fontsize_col = 10,
                                               color = my_palette,
                                               breaks = my_breaks,
                                               heatmap_legend_param = list(
                                                 title = "Spearman’s\np (rho)",
                                                 title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 11)
                                               )
)
# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))
print(bladder_V_gene_usage_clinVars_plot)

#pdf("bladder_V_gene_usage_clinVars.pdf", width = 12, height = 6)
#print(bladder_V_gene_usage_clinVars_plot)
#dev.off()

#rownames(heatmap_D_bladder) <- NULL

bladder_D_gene_usage_clinVars_plot <- pheatmap(heatmap_D_bladder,
                                               cluster_rows = FALSE,
                                               cluster_cols = FALSE,
                                               main = "Bladder IgH D Gene Usage",
                                               display_numbers = signif_D_bladder,
                                               fontsize_row = 10,
                                               fontsize_col = 15,
                                               cluster_row = FALSE,
                                               color = my_palette,
                                               breaks = my_breaks,
                                               legend = FALSE,
                                               heatmap_legend_param = list(
                                                 title = "Spearman’s\np (rho)",
                                                 title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 11)
                                                 )
)

# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))

print(bladder_D_gene_usage_clinVars_plot)

# pdf("bladder_D_gene_usage_clinVars.pdf", width = 8, height = 6)
# print(bladder_D_gene_usage_clinVars_plot)
# dev.off()

bladder_J_gene_usage_clinVars_plot <- pheatmap(heatmap_J_bladder,
                                               cluster_rows = FALSE,
                                               cluster_cols = FALSE,
                                               main = "Bladder IgH J Gene Usage",
                                               display_numbers = signif_J_bladder,
                                               fontsize_row = 10,
                                               fontsize_col = 15,
                                               color = my_palette,
                                               breaks = my_breaks,
                                               heatmap_legend_param = list(
                                                 title = "Spearman’s\np (rho)",
                                                 title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                 labels_gp = gpar(fontsize = 11)
                                               )
)
# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))

print(bladder_J_gene_usage_clinVars_plot)

#pdf("bladder_J_gene_usage_clinVars.pdf", width = 8, height = 6)
#print(bladder_J_gene_usage_clinVars_plot)
#dev.off()

# Save bladder
# png("bladder_V_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(bladder_V_gene_usage_clinVars_plot)
# dev.off()
# png("bladder_D_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(bladder_D_gene_usage_clinVars_plot)
# dev.off()
# png("bladder_J_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(bladder_J_gene_usage_clinVars_plot)
# dev.off()

################################################################################

# Heatmap: Renal ###############################################################
row_names_renal <- c(
  "Age: 32–86", 
  "Gender: Female (0) / Male (1)", 
  "Health Status: Healthy (0) / Cancer (1)", 
  "Charlson Index: 1–12", 
  "Leibscore: 0–11", 
  "Fuhrgrad: 1–4", 
  "C-Reactive Protein: 1–275", 
  "Hemoglobin: 4.6–11.4", 
  "Lactate Dehydrogenase: 99–544", 
  "Sodium: 130–145", 
  "Neutrophil Granulocytes: 1.6–21.9"
)
rownames(heatmap_V_renal) <- row_names_renal
rownames(heatmap_D_renal) <- row_names_renal
rownames(heatmap_J_renal) <- row_names_renal

# Product of Spearman test statistic (S) and correlation coefficient (rho))
renal_V_gene_usage_clinVars_plot <- pheatmap(heatmap_V_renal,
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             main = "Renal IgH V Gene Usage",
                                             display_numbers = signif_V_renal,
                                             fontsize_row =10,
                                             fontsize_col = 10,
                                             color = my_palette,
                                             breaks = my_breaks,
                                             heatmap_legend_param = list(
                                               title = "Spearman’s\np (rho)",
                                               title_gp = gpar(fontsize = 11, fontface = "bold"),
                                               labels_gp = gpar(fontsize = 11)
                                             )
)
# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))

# pdf("renal_V_gene_usage_clinVars.pdf", width = 12, height = 6)
# print(renal_V_gene_usage_clinVars_plot)
# dev.off()


renal_D_gene_usage_clinVars_plot <- pheatmap(heatmap_D_renal,
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             main = "Renal D Gene Usage",
                                             display_numbers = signif_D_renal,
                                             fontsize_row = 10,
                                             fontsize_col = 15,
                                             color = my_palette,
                                             breaks = my_breaks,
                                             legend = FALSE
)
#heatmap_legend_param = list(
# title = "Spearman’s\np (rho)",
#  title_gp = gpar(fontsize = 12, fontface = "bold"),
#  labels_gp = gpar(fontsize = 10)
#  )
#)
# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))

print(renal_D_gene_usage_clinVars_plot)

# pdf("renal_D_gene_usage_clinVars.pdf", width = 8, height = 6)
# print(renal_D_gene_usage_clinVars_plot)
# dev.off()

renal_J_gene_usage_clinVars_plot <- pheatmap(heatmap_J_renal,
                                             cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             main = "Renal IgH J Gene Usage",
                                             display_numbers = signif_J_renal,
                                             fontsize_row = 10,
                                             fontsize_col = 15,
                                             color = my_palette,
                                             breaks = my_breaks,
                                             heatmap_legend_param = list(
                                               title = "Spearman’s\np (rho)",
                                               title_gp = gpar(fontsize = 11, fontface = "bold"),
                                               labels_gp = gpar(fontsize = 11)
                                             )
)

# Add custom text annotations for significance symbols
#grid.text("Significance:\n* p-value <= 0.05\n** p-value <= 0.01\n*** p-value <= 0.001", x = 0.9, y = 0.1, gp = gpar(fontsize = 8))#, fontface = "bold"))

print(renal_J_gene_usage_clinVars_plot)

# pdf("renal_J_gene_usage_clinVars.pdf", width = 8, height = 6)
# print(renal_J_gene_usage_clinVars_plot)
# dev.off()

#Save renal
# png("renal_V_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(renal_V_gene_usage_clinVars_plot)
# dev.off()
# png("renal_D_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(renal_D_gene_usage_clinVars_plot)
# dev.off()
# png("renal_J_gene_usage_clinVars.png", width = 10, height = 5, units = "in", res = 300)
# print(renal_J_gene_usage_clinVars_plot)
# dev.off()

