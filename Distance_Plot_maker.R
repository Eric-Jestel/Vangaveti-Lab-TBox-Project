# Load required libraries
library(tidyverse)
library(patchwork)
library(ggplot2)

# Define the base directory and groups
base_dir <- "."
group_dirs <- c("PKZ18_Base_Copy", "PKZ18_124_Copy")
group_labels <- c("base", "124")

# Initialize empty dataframes
system_df <- data.frame()


# Loop through each group
for (i in seq_along(group_dirs)) {
  group_dir <- file.path(base_dir, group_dirs[i])
  group_label <- group_labels[i]
  
  # Loop through each replicate
  for (rep in 1:4) {
    rep_dir <- file.path(group_dir, paste0("PKZ18_Rep", rep))
    
    # Process SYSTEM file
    system_file <- file.path(rep_dir, "SLdistDrug.xvg")
    if (file.exists(system_file)) {
      system_data <- read.table(system_file, skip = 18)
      colnames(system_data) <- c("time", "rmsd")
      system_data$replicate <- rep
      system_data$group <- group_label
      system_df <- rbind(system_df, system_data)
    } else {
      warning(paste("File not found:", system_file))
    }
    
  }
}

# Convert to wide format if needed
system_wide <- system_df %>%
  pivot_wider(
    id_cols = time,
    names_from = c(group, replicate),
    names_prefix = "rmsd_",
    values_from = rmsd
  )


# Print summary
cat("Data loaded successfully!\n")
cat("System data dimensions:", dim(system_df), "\n")


# Convert time to nanoseconds and RMSD to angstroms
system_df <- system_df %>%
  mutate(
    time_ns = time / 1000, # Convert to nanoseconds if needed
    rmsd_ang = rmsd * 10   # Convert to angstroms
  )


t=  theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          #legend.text = element_text(size = 20, face = "bold"),
          #   legend.title = element_blank(),
          plot.title=element_text(hjust=0.5,face="bold", size=14),
          axis.title=element_text(size=18, face="bold"),
          axis.text=element_text(size=20,face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))


# Create a color palette
base_colors <- c("#0D47A1", "#1976D2", "#42A5F5", "#90CAF9") # Blue shades
p124_colors <- c("#B71C1C", "#EF5350","#E53935", "#FFCDD2") # Red shades
all_colors <- c(base_colors, p124_colors)

# Create the plots
# 1. System RMSD for Base variants
p1 <- ggplot(system_df %>% filter(group == "base"), 
             aes(x = time_ns, y = rmsd_ang, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = base_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,40))+
  labs(title = "Distance - PKZ18-00",
       x = "Time (ns)",
       y = "Distance (Å)") +
  t

# 2. System RMSD for 124 variants
p2 <- ggplot(system_df %>% filter(group == "124"), 
             aes(x = time_ns, y = rmsd_ang, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = p124_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,40))+
  labs(title = "Distance - PKZ18-124h",
       x = "Time (ns)",
       y = "Distance (Å)") +
  t


# Combine all plots using patchwork
combined_plot <- (p1 + p2) + 
  plot_annotation(
    title = "RMSD Comparison Between Base and 124 Variants",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# Display the combined plot
combined_plot

# Save the plots
ggsave("system_base_rmsd.png", p1, width = 8, height = 6, dpi = 300)
ggsave("system_124_rmsd.png", p2, width = 8, height = 6, dpi = 300)
ggsave("specloop_base_rmsd.png", p3, width = 8, height = 6, dpi = 300)
ggsave("specloop_124_rmsd.png", p4, width = 8, height = 6, dpi = 300)
ggsave("combined_distance.png", combined_plot, width = 14, height = 10, dpi = 600)

# Calculate mean, standard deviation, and standard error for each group and replicate
stats_df <- system_df %>%
  group_by(group, replicate) %>%
  summarize(
    mean_rmsd = mean(rmsd_ang, na.rm = TRUE),
    sd_rmsd = sd(rmsd_ang, na.rm = TRUE),
    se_rmsd = sd_rmsd / sqrt(n()),
    .groups = 'drop'
  )

# Print statistics
print(stats_df)

# Create mean-subtracted dataframe
system_df_centered <- system_df %>%
  left_join(stats_df, by = c("group", "replicate")) %>%
  mutate(rmsd_ang_centered = rmsd_ang - mean_rmsd)

# Create the plots with centered data
# 1. System RMSD for Base variants (centered)
p1_centered <- ggplot(system_df_centered %>% filter(group == "base"), 
             aes(x = time_ns, y = rmsd_ang_centered, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = base_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-15, 40)) +
    labs(title = "Distance - PKZ18-00 (Mean-Centered)",
       x = "Time (ns)",
       y = "Distance (Å) - Mean Centered") +
  t

# 2. System RMSD for 124 variants (centered)
p2_centered <- ggplot(system_df_centered %>% filter(group == "124"), 
             aes(x = time_ns, y = rmsd_ang_centered, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = p124_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-15, 40)) +
    labs(title = "Distance - PKZ18-124h (Mean-Centered)",
       x = "Time (ns)",
       y = "Distance (Å) - Mean Centered") +
  t

# Combine centered plots
combined_plot_centered <- (p1_centered + p2_centered) + 
  plot_annotation(
    title = "Mean-Centered RMSD Comparison Between Base and 124 Variants",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# Display the combined centered plot
combined_plot_centered

# Save the centered plots
ggsave("system_base_rmsd_centered.png", p1_centered, width = 8, height = 6, dpi = 300)
ggsave("system_124_rmsd_centered.png", p2_centered, width = 8, height = 6, dpi = 300)
ggsave("combined_distance_centered.png", combined_plot_centered, width = 14, height = 10, dpi = 600)




