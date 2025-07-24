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
specloop_df <- data.frame()

# Loop through each group
for (i in seq_along(group_dirs)) {
  group_dir <- file.path(base_dir, group_dirs[i])
  group_label <- group_labels[i]
  
  # Loop through each replicate
  for (rep in 1:4) {
    rep_dir <- file.path(group_dir, paste0("PKZ18_Rep", rep))
    
    # Process SYSTEM file
    system_file <- file.path(rep_dir, "RMSD_SYSTEM.xvg")
    if (file.exists(system_file)) {
      system_data <- read.table(system_file, skip = 18)
      colnames(system_data) <- c("time", "rmsd")
      system_data$replicate <- rep
      system_data$group <- group_label
      system_df <- rbind(system_df, system_data)
    } else {
      warning(paste("File not found:", system_file))
    }
    
    # Process SPECLOOP file
    specloop_file <- file.path(rep_dir, "RMSD_SPECLOOP.xvg")
    if (file.exists(specloop_file)) {
      specloop_data <- read.table(specloop_file, skip = 18)
      colnames(specloop_data) <- c("time", "rmsd")
      specloop_data$replicate <- rep
      specloop_data$group <- group_label
      specloop_df <- rbind(specloop_df, specloop_data)
    } else {
      warning(paste("File not found:", specloop_file))
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

specloop_wide <- specloop_df %>%
  pivot_wider(
    id_cols = time,
    names_from = c(group, replicate),
    names_prefix = "rmsd_",
    values_from = rmsd
  )

# Print summary
cat("Data loaded successfully!\n")
cat("System data dimensions:", dim(system_df), "\n")
cat("Specloop data dimensions:", dim(specloop_df), "\n")

# Save the data
save(system_df, specloop_df, system_wide, specloop_wide, 
     file = "rmsd_data.RData")

# Convert time to nanoseconds and RMSD to angstroms
system_df <- system_df %>%
  mutate(
    time_ns = time / 1000, # Convert to nanoseconds if needed
    rmsd_ang = rmsd * 10   # Convert to angstroms
  )

specloop_df <- specloop_df %>%
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
  scale_y_continuous(limits = c(0,20))+
  labs(title = "System RMSD - Base Variant",
       x = "Time (ns)",
       y = "RMSD (Å)") +
  t

# 2. System RMSD for 124 variants
p2 <- ggplot(system_df %>% filter(group == "124"), 
             aes(x = time_ns, y = rmsd_ang, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = p124_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,20))+
  labs(title = "System RMSD - 124 Variant",
       x = "Time (ns)",
       y = "RMSD (Å)") +
  t

# 3. Specloop RMSD for Base variants
p3 <- ggplot(specloop_df %>% filter(group == "base"), 
             aes(x = time_ns, y = rmsd_ang, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = base_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,10))+
  labs(title = "Specifier Loop RMSD - Base Variant",
       x = "Time (ns)",
       y = "RMSD (Å)") +
  t

# 4. Specloop RMSD for 124 variants
p4 <- ggplot(specloop_df %>% filter(group == "124"), 
             aes(x = time_ns, y = rmsd_ang, color = as.factor(replicate))) +
  geom_line(linewidth = .5) +
  scale_color_manual(values = p124_colors, name = "Replicate") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0,10))+
  labs(title = "Specifier Loop RMSD - 124 Variant",
       x = "Time (ns)",
       y = "RMSD (Å)") +
  t

# Combine all plots using patchwork
combined_plot <- (p1 + p2) / (p3 + p4) + 
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
ggsave("combined_rmsd_plot.png", combined_plot, width = 14, height = 10, dpi = 600)




