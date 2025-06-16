# Clonotype Analysis, Liliane Bader, 15.Juni 2025

# Renal and Bladder Datasets ###################################################
source("config.R")
renal_data <- readRDS(renal_data_path)
bladder_data <- readRDS(bladder_data_path)

# Libraries ####################################################################
library(dplyr)
library(ggplot2)
library(tidyr)

################################################################################

# Unique clonotypes ############################################################
compute_unique_clonotypes <- function(data_name, data, sample_id){
  unique_clonotypes = data.frame(
    Sample = sample_id,
    Unique_IGH_Clonotypes = sapply(data$mixcr_df_IGH, nrow),
    Unique_IGK_Clonotypes = sapply(data$mixcr_df_IGK, nrow),
    Unique_IGL_Clonotypes = sapply(data$mixcr_df_IGL, nrow)
  )
  # Reshape data to long format for ggplot
  unique_clonotypes_long <- unique_clonotypes %>% 
    pivot_longer(cols = -Sample, names_to = "Chain", values_to = "Unique_Clonotypes")
  return(unique_clonotypes_long)
}
bladder_unique_clonotypes <- compute_unique_clonotypes("Bladder Data", bladder_data, bladder_data$sample)
renal_unique_clonotypes <- compute_unique_clonotypes("Renal Data", renal_data, renal_data$SampleId)

# Recode Chain values
recode_chains <- function(df_long) {
  df_long %>%
    mutate(Chain = recode(Chain,
                          "Unique_IGH_Clonotypes" = "IgH",
                          "Unique_IGK_Clonotypes" = "IgK",
                          "Unique_IGL_Clonotypes" = "IgL"))
}
bladder_clonotype_long <- recode_chains(bladder_unique_clonotypes)
renal_clonotype_long <- recode_chains(renal_unique_clonotypes)

# Boxplots
plot_clonotype_boxplot <- function(data_long, title = NULL) {
  ggplot(data_long, aes(x = Chain, y = Unique_Clonotypes, fill = Chain)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_minimal() +
    labs(title = title, x = NULL, y = "Unique Clonotypes") +
    #scale_x_discrete(position = "top") + 
    scale_fill_manual(values = c("blue", "orange", "yellow")) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      plot.title = element_text(size = 25)
    )
}
bladder_unique_clontype_plot <- plot_clonotype_boxplot(bladder_clonotype_long, title = "Bladder: Unique Clonotypes")
renal_unique_clontype_plot <- plot_clonotype_boxplot(renal_clonotype_long, title = "Renal: Unique Clonotypes")

print(bladder_unique_clontype_plot)
print(renal_unique_clontype_plot)

# Save the plot as .pdf for thesis
# ggsave("bladder_unique_clonotypes.pdf", plot = bladder_unique_clontype_plot, width = 6, height = 5, units = "in")
# ggsave("renal_unique_clonotypes.pdf", plot = renal_unique_clontype_plot, width = 6, height = 5, units = "in")

# Save the plot as .png for defense
#ggsave("bladder_unique_clonotypes.png", plot = bladder_unique_clontype_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("renal_unique_clonotypes.png", plot = renal_unique_clontype_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

################################################################################

# Read Depth ###################################################################
compute_read_depth <- function(data, sample_ids) {
  read_depth_summary = data.frame(
    Sample = sample_ids,
    IgH = sapply(data$mixcr_df_IGH, function(df) sum(df$readCount)),
    IgK = sapply(data$mixcr_df_IGK, function(df) sum(df$readCount)),
    IgL = sapply(data$mixcr_df_IGL, function(df) sum(df$readCount))
  )
  # Long format for ggplot
  read_depth_long <- read_depth_summary %>%
    pivot_longer(cols = -Sample, names_to = "Chain", values_to = "Total_Reads")
  return(read_depth_long)
}
read_depth_bladder <- compute_read_depth(bladder_data, bladder_data$sample)
read_depth_renal <- compute_read_depth(renal_data, renal_data$SampleId)

# Plot boxplot
plot_read_depth_boxplot <- function(data_long, title = NULL) {
  ggplot(data_long, aes(x = Chain, y = Total_Reads, fill = Chain)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_minimal() +
    labs(title = title, x = NULL, y = "Read Depth") +
    #scale_x_discrete(position = "top") + 
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values = c("blue", "orange", "yellow")) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 25),
      axis.text.y = element_text(size = 25),
      axis.text.x = element_text(size = 25),
      plot.title = element_text(size = 25)
    )
}
bladder_read_depth_plot <- plot_read_depth_boxplot(read_depth_bladder, title = "Bladder: Read Depth")
renal_read_depth_plot   <- plot_read_depth_boxplot(read_depth_renal, title = "Renal: Read Depth")

print(bladder_read_depth_plot)
print(renal_read_depth_plot)

# Save the plot as .pdf for the thesis
#ggsave("bladder_read_depth.pdf", plot = bladder_read_depth_plot, width = 6, height = 5, units = "in")
#ggsave("renal_read_depth.pdf", plot = renal_read_depth_plot, width = 6, height = 5, units = "in")

# Save the plot as .png for defense
#ggsave("bladder_read_depth.png", plot = bladder_read_depth_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("renal_read_depth.png", plot = renal_read_depth_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

################################################################################

# Expanded vs. Non-expanded ####################################################
compute_expanded_clonotypes_all_chains <- function(data, sample_ids) {
  chains <- c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")
  chain_labels <- c("IGH", "IGK", "IGL")
  
  bind_rows(lapply(seq_along(chains), function(chain_index) {
    chain_name <- chains[chain_index]
    chain_label <- chain_labels[chain_index]
    
    bind_rows(lapply(seq_along(data[[chain_name]]), function(i) {
      df <- data[[chain_name]][[i]]
      df %>%
        mutate(expanded = factor(ifelse(uniqueMoleculeCount > 1, "Expanded", "Non-Expanded"),
                                 levels = c("Expanded", "Non-Expanded")),
               SampleId = sample_ids[i],
               Chain = chain_label)
    }))
  }))
}
renal_expanded_clonotypes_all <- compute_expanded_clonotypes_all_chains(renal_data, renal_data$SampleId)
bladder_expanded_clonotypes_all <- compute_expanded_clonotypes_all_chains(bladder_data, bladder_data$sample)

# Bar Proportion plot
plot_expansion_proportion <- function(expanded_data, title = NULL) {
  ggplot(expanded_data, aes(x = "", fill = expanded)) +
    geom_bar(position = "fill", width = 0.4) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("Expanded" = "orangered", "Non-Expanded" = "royalblue")) +
    labs(title = title, x = NULL, y = NULL, fill = "Expansion Status") +
    facet_wrap(~ Chain, strip.position = "bottom") +
    theme_minimal() +
    theme(
      #legend.position = "none",
      strip.text = element_text(size = 23),
      axis.text.y = element_text(size = 23),
      plot.title = element_text(size = 23),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)
    )
}
bladder_expansion_proportion_plot <- plot_expansion_proportion(bladder_expanded_clonotypes_all, title = "Bladder: Proportion of Clonotypes")
renal_expansion_proportion_plot <- plot_expansion_proportion(renal_expanded_clonotypes_all, title = "Renal: Proportion of Clonotypes")

print(bladder_expansion_proportion_plot)
print(renal_expansion_proportion_plot)

# Save the plot as .pdf for the thesis
#ggsave("bladder_expansion_proportion.pdf", plot = bladder_expansion_proportion_plot, width = 6, height = 5, units = "in")
#ggsave("renal_expansion_proportion.pdf", plot = renal_expansion_proportion_plot, width = 6, height = 5, units = "in")

# Save the plot as .png for defense
#ggsave("bladder_expansion_proportion.png", plot = bladder_expansion_proportion_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("renal_expansion_proportion.png", plot = renal_expansion_proportion_plot, width = 6, height = 5, units = "in", dpi = 300, bg = "white")

################################################################################

# Expanded vs. Non-expanded: Frequency #########################################
compute_clonotype_frequencies <- function(data, sample_ids) {
  chains <- c("mixcr_df_IGH", "mixcr_df_IGK", "mixcr_df_IGL")
  chain_labels <- c("IGH", "IGK", "IGL")
  
  bind_rows(lapply(seq_along(chains), function(chain_index) {
    chain_name <- chains[chain_index]
    chain_label <- chain_labels[chain_index]
    
    bind_rows(lapply(seq_along(data[[chain_name]]), function(i) {
      df <- data[[chain_name]][[i]] %>%
        mutate(
          expanded = ifelse(uniqueMoleculeCount > 1, "Expanded", "Non-Expanded"),
          SampleId = sample_ids[i],
          Chain = chain_label
        )
      
      total_count <- sum(df$uniqueMoleculeCount, na.rm = TRUE)
      
      df %>%
        group_by(SampleId, expanded, Chain) %>%
        summarise(
          total_UMI = sum(uniqueMoleculeCount, na.rm = TRUE),
          freq = total_UMI / total_count,
          .groups = "drop"
        )
    }))
  }))
}
renal_freqs <- compute_clonotype_frequencies(renal_data, renal_data$SampleId)
bladder_freqs <- compute_clonotype_frequencies(bladder_data, bladder_data$sample)

compute_abundance_summary <- function(df_long) {
  df_long %>%
    group_by(Chain, expanded) %>%
    summarise(total_UMI = sum(uniqueMoleculeCount, na.rm = TRUE), .groups = "drop") %>%
    group_by(Chain) %>%
    mutate(freq = total_UMI / sum(total_UMI)) %>%
    ungroup()
}
bladder_abundance_summary <- compute_abundance_summary(bladder_expanded_clonotypes_all)
renal_abundance_summary   <- compute_abundance_summary(renal_expanded_clonotypes_all)

plot_expansion_frequency <- function(summary_df, title = NULL) {
  ggplot(summary_df, aes(x = Chain, y = freq, fill = expanded)) +
    geom_bar(stat = "identity", position = "fill", width = 0.4) +
    scale_y_continuous(labels = scales::percent_format()) +
    #scale_x_discrete(position = "top") + 
    scale_fill_manual(values = c("Expanded" = "orangered", "Non-Expanded" = "royalblue")) +
    labs(title = title, x = NULL, y = NULL, fill = "Expansion Status") +
    theme_minimal() +
    theme(
      #legend.position = "none",
      axis.text = element_text(size = 23),
      axis.text.y = element_text(size = 23), #element_blank()
      plot.title = element_text(size = 23),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
    )
}
bladder_expanded_freq_chain <- plot_expansion_frequency(bladder_abundance_summary, title = "Bladder: Proportion of Repertoire")
renal_expanded_freq_chain   <- plot_expansion_frequency(renal_abundance_summary, title = "Renal: Proportion of Repertoire")

print(bladder_expanded_freq_chain)
print(renal_expanded_freq_chain)

# Save the plot as .pdf for the thesis
#ggsave("bladder_expansion_abundance.pdf", plot = bladder_expanded_freq_chain, width = 6, height = 5, units = "in")
#ggsave("renal_expansion_abundance.pdf", plot = renal_expanded_freq_chain, width = 6, height = 5, units = "in")

# Save the plot as .png for defense
#ggsave("bladder_expansion_abundance.png", plot = bladder_expanded_freq_chain, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
#ggsave("renal_expansion_abundance.png", plot = renal_expanded_freq_chain, width = 6, height = 5, units = "in", dpi = 300, bg = "white")
