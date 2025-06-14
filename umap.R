# Umap, Liliane Bader, 15.Juni 2025

# Renal and Bladder Datasets ###################################################
source("config.R")
renal_data <- readRDS(renal_data_path)
bladder_data <- readRDS(bladder_data_path)

# Libraries ####################################################################
library(umap)
library(ggplot2)
library(tidyverse)
################################################################################

# Gene Usage ###################################################################
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
normalized_gene_usage_bladder <- process_gene_usage(bladder_data, "sample")

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

# Umap: Individual chains ######################################################
prepare_umap_data <- function(gene_usage_df, gene_type) {
  # Filter out empty gene names and create a placeholder if needed
  gene_usage_df <- gene_usage_df %>%
    filter(!!sym(gene_type) != "") %>%  # Remove empty strings
    mutate(!!sym(gene_type) := ifelse(!!sym(gene_type) == "", 
                                      "Unknown", 
                                      !!sym(gene_type)))  # Or replace with "Unknown"
  
  # Proceed with pivoting
  gene_usage_df %>%
    select(sample, !!sym(gene_type), freq) %>%
    pivot_wider(names_from = !!sym(gene_type), 
                values_from = freq,
                values_fill = 0) %>%
    column_to_rownames("sample") %>%
    as.matrix()
}

# Modified plot_gene_umap to accept parameters
plot_gene_umap <- function(umap_data, title, n_neighbors = 15, min_dist = 0.1) {
  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- min_dist
  
  umap_results <- umap(umap_data, config = umap_config)
  umap_df <- as.data.frame(umap_results$layout)
  umap_df$sample <- rownames(umap_data)
  
  ggplot(umap_df, aes(x = V1, y = V2)) +
    geom_point() +
    ggtitle(title) +
    theme_minimal()
}

# Function to plot UMAP
plot_all_chains_genes_separately <- function(neighbors = 15, distance = 0.1, tissue = "both") {
  # Define all combinations we want to plot
  chain_gene_combinations <- list(
    list(data = v_usage_IGH_bladder, gene = "V.Name", chain = "IGH", type = "V", tissue_type = "Bladder"),
    list(data = d_usage_IGH_bladder, gene = "D.Name", chain = "IGH", type = "D", tissue_type = "Bladder"),
    list(data = j_usage_IGH_bladder, gene = "J.Name", chain = "IGH", type = "J", tissue_type = "Bladder"),
    list(data = v_usage_IGK_bladder, gene = "V.Name", chain = "IGK", type = "V", tissue_type = "Bladder"),
    list(data = j_usage_IGK_bladder, gene = "J.Name", chain = "IGK", type = "J", tissue_type = "Bladder"),
    list(data = v_usage_IGL_bladder, gene = "V.Name", chain = "IGL", type = "V", tissue_type = "Bladder"),
    list(data = j_usage_IGL_bladder, gene = "J.Name", chain = "IGL", type = "J", tissue_type = "Bladder"),
    list(data = v_usage_IGH_renal, gene = "V.Name", chain = "IGH", type = "V", tissue_type = "Renal"),
    list(data = d_usage_IGH_renal, gene = "D.Name", chain = "IGH", type = "D", tissue_type = "Renal"),
    list(data = j_usage_IGH_renal, gene = "J.Name", chain = "IGH", type = "J", tissue_type = "Renal"),
    list(data = v_usage_IGK_renal, gene = "V.Name", chain = "IGK", type = "V", tissue_type = "Renal"),
    list(data = j_usage_IGK_renal, gene = "J.Name", chain = "IGK", type = "J", tissue_type = "Renal"),
    list(data = v_usage_IGL_renal, gene = "V.Name", chain = "IGL", type = "V", tissue_type = "Renal"),
    list(data = j_usage_IGL_renal, gene = "J.Name", chain = "IGL", type = "J", tissue_type = "Renal")
  )
  
  # Filter based on tissue type if needed
  if (tissue %in% c("Bladder", "Renal")) {
    chain_gene_combinations <- Filter(function(x) x$tissue_type == tissue, chain_gene_combinations)
  } else if (tissue != "both") {
    stop("tissue must be 'Bladder', 'Renal', or 'both'")
  }
  
  # Process each combination
  for (combo in chain_gene_combinations) {
    # Prepare data
    matrix_data <- prepare_umap_data(combo$data, combo$gene)
    
    # Create title
    title <- sprintf("UMAP of %s %s Gene Usage %s (n=%d, dist=%.2f)",
                     combo$chain, combo$type, combo$tissue_type,
                     neighbors, distance)
    
    # Plot
    p <- plot_gene_umap(matrix_data, title, n_neighbors = neighbors, min_dist = distance)
    print(p)
  }
}
# n_neighbors: 5 to 50
# min_dist: 0.01 to 0.5
#plot_all_chains_genes_separately(neighbors = 50, distance = 0.1, tissue = "both")
################################################################################

# Umap: all chains combined ####################################################
# Combine all gene usage data
combined_usage_bladder <- bind_rows(
  v_usage_IGH_bladder %>% mutate(chain_gene = paste0("IGH_V_", V.Name)),
  d_usage_IGH_bladder %>% mutate(chain_gene = paste0("IGH_D_", D.Name)),
  j_usage_IGH_bladder %>% mutate(chain_gene = paste0("IGH_J_", J.Name)),
  v_usage_IGK_bladder %>% mutate(chain_gene = paste0("IGK_V_", V.Name)),
  j_usage_IGK_bladder %>% mutate(chain_gene = paste0("IGK_J_", J.Name)),
  v_usage_IGL_bladder %>% mutate(chain_gene = paste0("IGL_V_", V.Name)),
  j_usage_IGL_bladder %>% mutate(chain_gene = paste0("IGL_J_", J.Name))
) %>%
  select(sample, chain_gene, freq) %>%
  pivot_wider(names_from = chain_gene, values_from = freq, values_fill = 0)

combined_usage_renal <- bind_rows(
  v_usage_IGH_renal %>% mutate(chain_gene = paste0("IGH_V_", V.Name)),
  d_usage_IGH_renal %>% mutate(chain_gene = paste0("IGH_D_", D.Name)),
  j_usage_IGH_renal %>% mutate(chain_gene = paste0("IGH_J_", J.Name)),
  v_usage_IGK_renal %>% mutate(chain_gene = paste0("IGK_V_", V.Name)),
  j_usage_IGK_renal %>% mutate(chain_gene = paste0("IGK_J_", J.Name)),
  v_usage_IGL_renal %>% mutate(chain_gene = paste0("IGL_V_", V.Name)),
  j_usage_IGL_renal %>% mutate(chain_gene = paste0("IGL_J_", J.Name))
) %>%
  select(sample, chain_gene, freq) %>%
  pivot_wider(names_from = chain_gene, values_from = freq, values_fill = 0)

# Prepare matrix for UMAP
combined_matrix_bladder <- combined_usage_bladder %>%
  column_to_rownames("sample") %>%
  as.matrix()

combined_matrix_renal <- combined_usage_renal %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Run and plot UMAP
umap_results_bladder <- umap(combined_matrix_bladder)
umap_df_bladder <- as.data.frame(umap_results_bladder$layout)
umap_df_bladder$sample <- rownames(combined_matrix_bladder)

ggplot(umap_df_bladder, aes(x = V1, y = V2)) +
  geom_point() +
  ggtitle("UMAP of Combined BCR Gene Usage Bladder") +
  theme_minimal()

umap_results_renal <- umap(combined_matrix_renal)
umap_df_renal <- as.data.frame(umap_results_renal$layout)
umap_df_renal$sample <- rownames(combined_matrix_renal)

ggplot(umap_df_renal, aes(x = V1, y = V2)) +
  geom_point() +
  ggtitle("UMAP of Combined BCR Gene Usage Renal") +
  theme_minimal()
################################################################################

# Umap: IGHV colored by Sex and Age ############################################
prepare_umap_data <- function(usage_df, gene_column) {
  usage_df %>%
    select(sample, !!sym(gene_column), freq) %>%
    pivot_wider(names_from = !!sym(gene_column), values_from = freq, values_fill = 0) %>%
    column_to_rownames(var = "sample") %>%
    as.matrix()
}

run_umap_gene_usage <- function(usage_df, gene_column, config = umap.defaults) {
  usage_matrix <- prepare_umap_data(usage_df, gene_column)
  umap_result <- umap(usage_matrix, config = config)
  
  umap_df <- as.data.frame(umap_result$layout)
  umap_df$sample <- rownames(usage_matrix)
  
  return(umap_df)
}
plot_umap <- function(umap_df, meta_df, color_col, title, palette = "Set1", gradient = FALSE, size = 3) {
  merged_df <- left_join(umap_df, meta_df, by = "sample")
  
  p <- ggplot(merged_df, aes(x = V1, y = V2, color = .data[[color_col]])) +
    geom_point(size = size) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text = element_text(size = 13)
    )
  
  if (gradient) {
    p <- p + scale_color_gradient(low = "blue", high = "red")
  } else {
    p <- p + scale_color_brewer(palette = palette)
  }
  
  return(p)
}
# Config for UMAP
custom_config <- umap.defaults
custom_config$n_neighbors <- 10
custom_config$min_dist <- 0.2
custom_config$metric <- "cosine"

# Run UMAP on V gene
umap_bladder_IGHV <- run_umap_gene_usage(v_usage_IGH_bladder, "V.Name", config = custom_config)

# Metadata
meta_bladder <- bladder_data %>%
  distinct(sample, Sex, Age)

# Plot UMAPs
umap_bladder_sex_plot <- plot_umap(umap_bladder_IGHV, meta_bladder, "Sex", "IGH V Usage (Bladder) by Sex")
umap_bladder_age_plot <- plot_umap(umap_bladder_IGHV, meta_bladder, "Age", "IGH V Usage (Bladder) by Age", gradient = TRUE)

print(umap_bladder_sex_plot)
print(umap_bladder_age_plot)
umap_renal_IGHV <- run_umap_gene_usage(v_usage_IGH_renal, "V.Name", config = custom_config)

meta_renal <- renal_data %>%
  distinct(SampleId, DonorStatus, AgeSample) %>%
  rename(sample = SampleId, health_status = DonorStatus, Age = AgeSample)

umap_renal_health_plot <- plot_umap(umap_renal_IGHV, meta_renal, "health_status", "IGH V Usage (Renal) by Health")
umap_renal_age_plot <- plot_umap(umap_renal_IGHV, meta_renal, "Age", "IGH V Usage (Renal) by Age", gradient = TRUE)

print(umap_renal_health_plot)
print(umap_renal_age_plot)

#ggsave("umap_igh_v_bladder_age.pdf", plot = umap_bladder_age_plot, width = 8, height = 6, units = "in")
#ggsave("umap_igh_v_renal_health_status.pdf", plot = umap_renal_age_plot, width = 8, height = 6, units = "in")
#ggsave("umap_igh_v_renal_age.pdf", plot = umap_renal_health_plot, width = 8, height = 6, units = "in")
