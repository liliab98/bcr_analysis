# Cancer vs Healthy Immune Repertoire Analysis, Liliane Bader, 15.Juni 2025

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
library(scales)
################################################################################

# Combine Datasets ##############################################################
#Common columns?
# renal_data$SampleID             - bladder_data$sample
# renal_data$AgeSample            - bladder_data$Age
# renal_data$DonorStatus          - 
# renal_data$Sex                  - bladder_data$Sex
# renal_data$Death                - bladder_data$OS_event
# renal_data$rec_at_any_time_sum  - bladder_data$Clinical_relapse
# The three mixcr columns ofc

# Standardize column names
bladder_data$SampleId <- bladder_data$sample 
renal_data$Age <- renal_data$AgeSample
bladder_data$DonorStatus <- "Bladder" 
bladder_data$Death <- bladder_data$OS_event
renal_data$Clinical_relapse <- renal_data$rec_at_any_time_sum

levels(renal_data$DonorStatus)[levels(renal_data$DonorStatus) == "Cancer"] <- "Renal"

# Identify common columns
common_cols <- intersect(names(renal_data), names(bladder_data))

# Subset both datasets to only include common columns
renal_subset <- renal_data[, common_cols, drop = FALSE]
bladder_subset <- bladder_data[, common_cols, drop = FALSE]

# Merge datasets
combined_data <- rbind(renal_subset, bladder_subset)

# Check structure
table(combined_data$DonorStatus)  
names(combined_data)
################################################################################

my_comparisons <- list(
  c("Healthy", "Renal"),
  c("Healthy", "Bladder")
)

colors <- c("Renal" = "#D55E00", "Bladder" = "#0072B2", "Healthy" = "#009E73")

# Read Counts ##################################################################
# Function to sum read counts per sample for each gene
# Function to calculate total read counts
get_total_reads <- function(mixcr_list) {
  sapply(mixcr_list, function(df) {
    if (!is.null(df) && nrow(df) > 0 && "readCount" %in% names(df)) {
      sum(df$readCount, na.rm = TRUE)
    } else {
      0
    }
  })
}
# Apply function to each gene
total_reads_IGH <- get_total_reads(combined_data$mixcr_df_IGH)
total_reads_IGK <- get_total_reads(combined_data$mixcr_df_IGK)
total_reads_IGL <- get_total_reads(combined_data$mixcr_df_IGL)

# Add per-gene read counts to data
combined_data$reads_IGH <- total_reads_IGH
combined_data$reads_IGK <- total_reads_IGK
combined_data$reads_IGL <- total_reads_IGL

# Combine reads per sample
combined_data$reads_combined <- total_reads_IGH + total_reads_IGK + total_reads_IGL

# Function to create long-format data for plotting
prepare_long_reads <- function(df, id_col, group_col, colors) {
  df_long <- df %>%
    select(all_of(c(id_col, group_col, "reads_combined", "reads_IGH", "reads_IGK", "reads_IGL"))) %>%
    pivot_longer(cols = starts_with("reads_"),
                 names_to = "Gene", values_to = "ReadCount") %>%
    mutate(Gene = gsub("reads_", "", Gene),
           Group = .data[[group_col]]) %>%
    group_by(.data[[id_col]]) %>%
    mutate(ReadCount_norm = ReadCount / sum(ReadCount)) %>%
    ungroup()
  
  df_long[[group_col]] <- factor(df_long[[group_col]], levels = names(colors))
  df_long$Group <- factor(df_long$Group, levels = names(colors))
  
  return(df_long)
}
long_reads <- prepare_long_reads(combined_data, id_col = "SampleId", group_col = "DonorStatus", colors = colors)

# Plot
plot_read_counts <- function(df, comparisons, colors, gene_filter = NULL, title = NULL, size_base = 20) {
  if (!is.null(gene_filter)) df <- df %>% filter(Gene %in% gene_filter)
  
  p <- ggplot(df, aes(x = Group, y = ReadCount)) +
    geom_boxplot(aes(fill = DonorStatus), outlier.shape = NA, color = "black", alpha = 0.6) +
    geom_jitter(aes(color = DonorStatus), width = 0.2, size = 3, alpha = 0.6) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = comparisons, 
                       label = "p.signif",
                       tip.length = 0.02,
                       size = size_base * 0.3) +
    scale_y_continuous(labels = scales::comma) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    xlab(NULL) +
    ylab("Read Count") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = size_base * 1.25),
      legend.position = "none",
      axis.text.y = element_text(size = size_base),
      axis.text.x = element_text(size = size_base),
      axis.title = element_text(size = size_base * 1.2),
      plot.title = element_text(size = size_base * 1.2, hjust = 0.5)
    ) +
    ggtitle(title)
  
  # Only facet if there's more than one unique gene
  if (length(unique(df$Gene)) > 1) {
    p <- p + facet_wrap(~ Gene, scales = "free_y", ncol = 2)
  }
  
  return(p)
}
read_counts_plot_thesis <- plot_read_counts(long_reads, my_comparisons, colors, size_base = 10) #, title = "Read Counts by Group")
read_counts_all_chains_plot_thesis <- plot_read_counts(long_reads, my_comparisons, colors, gene_filter = "combined", title = "Total Read Count")

print(read_counts_plot_thesis)
print(read_counts_all_chains_plot_thesis)

# Save the plot as .pdf for the thesis
#ggsave("healthyVSsick_read_counts.pdf", plot = read_counts_plot_thesis, width = 25, height = 25, units = "in")
#ggsave("healthyVSsick_read_counts_all_chains.pdf", plot = read_counts_all_chains_plot_thesis, width = 20, height = 15, units = "in")

# Save the plot as .png for defense
#ggsave("healthyVSsick_read_counts.png", plot = read_counts_plot_thesis, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("healthyVSsick_read_counts_all_chains.pdf", plot = read_counts_all_chains_plot_thesis, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

################################################################################

# Shannon Diversity ############################################################
# Compute normalized nd raw Shannon Diversity Index 
compute_normalized_shannon <- function(read_counts) {
  if (length(read_counts) == 0) return(NA)
  H <- diversity(read_counts)
  S <- length(read_counts)
  if (S > 1) H / log(S) else 0
}

compute_shannon_and_normalized <- function(data, chains = c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")) {
  n <- nrow(data)
  out <- vector("list", length = n)
  
  for (i in seq_len(n)) {
    shannon_vals <- list()
    norm_vals <- list()
    combined_reads <- numeric(0)
    
    for (chain in chains) {
      df <- data[[chain]][[i]]
      rc <- if (!is.null(df) && "readCount" %in% names(df)) df$readCount else numeric(0)
      shannon_vals[[chain]] <- if (length(rc)) diversity(rc) else NA
      norm_vals[[chain]] <- compute_normalized_shannon(rc)
      combined_reads <- c(combined_reads, rc)
    }
    
    shannon_vals[["Combined"]] <- diversity(combined_reads)
    norm_vals[["Combined"]] <- compute_normalized_shannon(combined_reads)
    
    out[[i]] <- data.frame(
      SampleId = data$SampleId[i],
      DonorStatus = data$DonorStatus[i],
      CancerType = if ("CancerType" %in% names(data)) data$CancerType[i] else NA,
      Chain = c("IGH", "IGK", "IGL", "Combined"),
      Shannon_Diversity = unlist(shannon_vals),
      Normalized_Shannon = unlist(norm_vals)
    )
  }
  
  do.call(rbind, out)
}

diversity_data <- compute_shannon_diversity(combined_data)

diversity_data$DonorStatus <- factor(diversity_data$DonorStatus, levels = c("Renal", "Healthy", "Bladder"))

# Plot
plot_shannon_diversity <- function(div_df, comparisons, colors, gene_filter = NULL, title = NULL, size_base = 40) {
  if (!is.null(gene_filter)) div_df <- div_df %>% filter(Chain %in% gene_filter)
  
  ggplot(div_df, aes(x = DonorStatus, y = Shannon_Diversity)) +
    geom_boxplot(aes(fill = DonorStatus), outlier.shape = NA, color = "black", alpha = 0.6) +
    geom_jitter(aes(color = DonorStatus), width = 0.2, size = 3, alpha = 0.6) +
    facet_wrap(~ Chain, scales = "free_y", ncol = 2) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = comparisons, 
                       label = "p.signif",
                       tip.length = 0.02,
                       size = size_base * 0.35) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = scales::comma) +
    xlab(NULL) +
    ylab("Shannon Diversity") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = size_base * 1.25),
      legend.position = "none",
      axis.text.y = element_text(size = size_base),
      axis.text.x = element_text(size = size_base),
      axis.title = element_text(size = size_base * 1.2),
      plot.title = element_text(size = size_base * 1.2, hjust = 0.5)
    ) +
    ggtitle(title)
}
# Plot all chains
shannon_healthy_sick_plot_thesis <- plot_shannon_diversity(
  div_df = diversity_data,
  comparisons = my_comparisons,
  colors = colors,
  gene_filter = NULL,
  title = "Shannon Diversity by Chain and Donor Status"
)
# Plot combined only
shannon_all_chains_plot_thesis <- plot_shannon_diversity(
  div_df = diversity_data,
  comparisons = my_comparisons,
  colors = colors,
  gene_filter = "Combined",
  title = "Combined Shannon Diversity"
)
#print(shannon_all_chains_plot_thesis)
#print(shannon_healthy_sick_plot_thesis)

#ggsave("healthyVSsick_shannon.pdf", plot = shannon_healthy_sick_plot_thesis, width = 25, height = 25, units = "in")
#ggsave("healthyVSsick_shannon_all_chains.pdf", plot = shannon_all_chains_plot_thesis, width = 20, height = 15, units = "in")

# Normalized Shannon Diversity #################################################
# Compute full diversity dataset
diversity_data <- compute_shannon_and_normalized(combined_data)

# Set factor levels
diversity_data$DonorStatus <- factor(diversity_data$DonorStatus, levels = c("Renal", "Healthy", "Bladder"))

# Plot normalized diversity
plot_normalized_shannon <- function(df, comparisons, colors, gene_filter = NULL, title = NULL, size_base = 20) {
  if (!is.null(gene_filter)) df <- df %>% filter(Chain %in% gene_filter)
  
  p <- ggplot(df, aes(x = DonorStatus, y = Normalized_Shannon)) +
    geom_boxplot(aes(fill = DonorStatus), outlier.shape = NA, color = "black", alpha = 0.6) +
    geom_jitter(aes(color = DonorStatus), width = 0.2, size = 3, alpha = 0.6) +
    #facet_wrap(~ Chain, scales = "free_y", ncol = 2) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = comparisons, 
                       label = "p.signif",
                       tip.length = 0.02,
                       size = size_base * 0.3) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    xlab(NULL) +
    ylab("Normalized Shannon Diversity") +
    theme_minimal() +
    theme(
      strip.text = element_text(size = size_base * 1.25),
      legend.position = "none",
      axis.text.y = element_text(size = size_base),
      axis.text.x = element_text(size = size_base),
      axis.title = element_text(size = size_base * 1.2),
      plot.title = element_text(size = size_base * 1.2, hjust = 0.5)
    ) +
    ggtitle(title)
  # Only facet if there's more than one unique gene
  if (length(unique(df$Chain)) > 1) {
    p <- p + facet_wrap(~ Chain, scales = "free_y", ncol = 2)
  }
  
  return(p)
}
normalized_shannon_all <- plot_normalized_shannon(
  df = diversity_data,
  comparisons = my_comparisons,
  colors = colors,
  gene_filter = NULL,
  #title = "Normalized Shannon Diversity by Chain",
  size_base = 10
)

normalized_shannon_combined <- plot_normalized_shannon(
  df = diversity_data,
  comparisons = my_comparisons,
  colors = colors,
  gene_filter = "Combined",
  title = "Normalized Shannon Diversity"
)

print(normalized_shannon_all)
print(normalized_shannon_combined)

# Save the plot as .pdf for the thesis
# ggsave("normalized_shannon_all_chains.pdf", normalized_shannon_all, width = 25, height = 25, units = "in")
# ggsave("normalized_shannon_combined.pdf", normalized_shannon_combined, width = 20, height = 15, units = "in")

# Save the plot as .png for defense
#ggsave("normalized_shannon_all_chains.png", plot = normalized_shannon_all, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("normalized_shannon_combined.png", plot = normalized_shannon_combined, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

################################################################################

# Clonal Expansion #############################################################
# Top N Clones: Fraction of total reads from the top 10 clones.
# 50% of Reads Clones: Number of clones contributing to the first 50% of reads.
# Fraction of total reads from top N clones
compute_top_n_clones <- function(df, top_n = 10) {
  if (!"readCount" %in% names(df)) stop("Missing 'readCount' column in input.")
  sorted_reads <- sort(df$readCount, decreasing = TRUE)
  top_reads <- sum(sorted_reads[1:min(top_n, length(sorted_reads))])
  total_reads <- sum(df$readCount)
  return(top_reads / total_reads)
}
# Number of clones that contribute to 50% of reads
compute_50pct_clones <- function(df) {
  if (!"readCount" %in% names(df)) stop("Missing 'readCount' column in input.")
  sorted_counts <- sort(df$readCount, decreasing = TRUE)
  cumulative_fraction <- cumsum(sorted_counts) / sum(sorted_counts)
  num_clones <- which(cumulative_fraction >= 0.5)[1]
  return(num_clones)
}
compute_clonal_metrics <- function(data, chains = c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL"), top_n = 10) {
  top_n_list <- list()
  pct50_list <- list()
  
  for (i in seq_len(nrow(data))) {
    donor_status <- data$DonorStatus[i]
    
    df_list <- lapply(chains, function(chain) data[[chain]][[i]])
    names(df_list) <- gsub("mixcr_df_", "", chains)
    
    for (chain in names(df_list)) {
      df <- df_list[[chain]]
      if (!is.null(df) && nrow(df) > 0 && "readCount" %in% names(df)) {
        top_n_list[[length(top_n_list) + 1]] <- data.frame(
          ClonalExpansion = compute_top_n_clones(df, top_n),
          DonorStatus = donor_status,
          Chain = chain,
          Metric = "Top 10 Clones"
        )
        pct50_list[[length(pct50_list) + 1]] <- data.frame(
          ClonalExpansion = compute_50pct_clones(df),
          DonorStatus = donor_status,
          Chain = chain,
          Metric = "Clones for 50% Reads"
        )
      }
    }
    
    combined_reads <- unlist(lapply(df_list, function(df) df$readCount))
    if (length(combined_reads) > 0) {
      top_n_list[[length(top_n_list) + 1]] <- data.frame(
        ClonalExpansion = compute_top_n_clones(data.frame(readCount = combined_reads)),
        DonorStatus = donor_status,
        Chain = "Combined",
        Metric = "Top 10 Clones"
      )
      pct50_list[[length(pct50_list) + 1]] <- data.frame(
        ClonalExpansion = compute_50pct_clones(data.frame(readCount = combined_reads)),
        DonorStatus = donor_status,
        Chain = "Combined",
        Metric = "Clones for 50% Reads"
      )
    }
  }
  
  list(
    top_n_df = do.call(rbind, top_n_list),
    pct_50_df = do.call(rbind, pct50_list)
  )
}

# Compute metrics
clonal_results <- compute_clonal_metrics(combined_data)

top_n_df <- clonal_results$top_n_df %>%
  mutate(DonorStatus = factor(DonorStatus, levels = c("Renal", "Healthy", "Bladder")))

pct_50_df <- clonal_results$pct_50_df %>%
  mutate(DonorStatus = factor(DonorStatus, levels = c("Renal", "Healthy", "Bladder")))

# Plot
plot_clonal_expansion <- function(df, ylab_text, scale_log = FALSE, gene_filter = NULL, title = NULL, size_base = 20) {
  if (!is.null(gene_filter)) df <- df %>% filter(Chain %in% gene_filter)
  
  p <- ggplot(df, aes(x = DonorStatus, y = ClonalExpansion)) +
    geom_boxplot(aes(fill = DonorStatus), outlier.shape = NA, color = "black", alpha = 0.6) +
    geom_jitter(aes(color = DonorStatus), width = 0.2, size = 3, alpha = 0.6) +
    #facet_wrap(~ Chain, scales = "free_y", ncol = 2) +
    stat_compare_means(method = "wilcox.test", 
                       comparisons = my_comparisons, 
                       label = "p.signif",
                       tip.length = 0.02,
                       size = size_base * 0.3) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    xlab(NULL) +
    ylab(ylab_text) +
    ggtitle(title) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = size_base * 1.25),
      legend.position = "none",
      axis.text = element_text(size = size_base),
      axis.title = element_text(size = size_base * 1.2),
      plot.title = element_text(size = size_base * 1.2, hjust = 0.5)
    )
  
  if (scale_log) {
    p <- p + scale_y_log10(labels = comma)
  } else {
    p <- p + scale_y_continuous(labels = comma)
  }
  
  # Only facet if there's more than one unique gene
  if (length(unique(df$Chain)) > 1) {
    p <- p + facet_wrap(~ Chain, scales = "free_y", ncol = 2)
  }
  
  return(p)
}

# Plot Top 10 clones (all chains)
plot_top10 <- plot_clonal_expansion(top_n_df, 
                                    ylab_text = "Fraction of Reads (log10)",
                                    scale_log = TRUE, 
                                    #title = "Top 10 Clonal Expansion",
                                    size_base = 10)

# Plot Clones for 50% reads (all chains)
plot_50pct <- plot_clonal_expansion(pct_50_df, 
                                    ylab_text = "Number of Clones", 
                                    scale_log = FALSE, 
                                    #title = "Clones for 50% of Reads",
                                    size_base = 10)

# Plot combined only
plot_top10_combined <- plot_clonal_expansion(top_n_df, 
                                             gene_filter = "Combined",
                                             ylab_text = "Fraction of Reads (log10)",
                                             scale_log = TRUE,
                                             title = "Top 10 Clones")

plot_50pct_combined <- plot_clonal_expansion(pct_50_df, 
                                             gene_filter = "Combined",
                                             ylab_text = "Number of Clones",
                                             scale_log = FALSE,
                                             title = "Clones for 50% of Reads")

# Display
print(plot_top10)
print(plot_50pct)
print(plot_top10_combined)
print(plot_50pct_combined)

# Save the plot as .pdf for the thesis
#ggsave("top10_clonal_expansion_all.pdf", plot_top10, width = 25, height = 25, units = "in")
#ggsave("clones_50pct_reads_all.pdf", plot_50pct, width = 25, height = 25, units = "in")
#ggsave("top10_clonal_expansion_combined.pdf", plot_top10_combined, width = 20, height = 15, units = "in")
#ggsave("clones_50pct_reads_combined.pdf", plot_50pct_combined, width = 20, height = 15, units = "in")

# Save the plot as .png for defense
#ggsave("top10_clonal_expansion_all.png", plot = plot_top10, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("clones_50pct_reads_all.png", plot = plot_50pct, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("top10_clonal_expansion_combined.png", plot = plot_top10_combined, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("clones_50pct_reads_combined.png", plot = plot_50pct_combined, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

