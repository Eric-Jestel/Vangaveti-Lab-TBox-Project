# ------------------------------
# SETTINGS
# ------------------------------
base_dirs  <- c("PKZ18_124_Copy", "PKZ18_Base_Copy")
reps       <- paste0("PKZ18_Rep", 1:4)
skip_lines <- 17

# containers
distance_list <- list()
angle_list    <- list()
time_vec      <- NULL

# ------------------------------
# LOOP OVER DIRECTORIES
# ------------------------------
for (base in base_dirs) {
  for (rep in reps) {
    dir_path <- file.path(base, rep)
    xvg_files <- list.files(dir_path, pattern="*_*_*\\.xvg$", full.names = TRUE)
    
    for (f in xvg_files) {
      # Parse filename: e.g. "A17_Thiazole_Angle.xvg"
      fname <- basename(f)
      parts <- strsplit(tools::file_path_sans_ext(fname), "_")[[1]]
      
      print(dir_path)
      print(fname)
      
      # Safety check
      if (length(parts) != 3) {
        warning("Filename format unexpected: ", fname)
        next
      }
      
      groupA  <- parts[1]
      groupB  <- parts[2]
      measure <- parts[3]  # should be "Angle" or "Distance"
      
      # Create unique key
      sample_key <- paste(base, rep, groupA, groupB, measure, sep = "_")
      
      # Read and clean lines manually
      lines <- readLines(f)
      data_lines <- lines[(skip_lines + 1):length(lines)]
      data_lines <- grep("^[^@#]", data_lines, value = TRUE)  # exclude header/comments
      dat <- read.table(text = data_lines, header = FALSE)
      
      # Extract time and value
      t <- dat[[1]]
      v <- dat[[2]]
      
      # Check or initialize time vector
      if (is.null(time_vec)) {
        time_vec <- t
      }
      
      # Store into appropriate list
      if (measure == "Distance") {
        distance_list[[sample_key]] <- v
      } else if (measure == "Angle") {
        v <- ifelse(v > 90, 180 - v, v)
        angle_list[[sample_key]] <- v
      }
    }
  }
}

# ------------------------------
# BUILD WIDE TABLES
# ------------------------------
distance_df <- data.frame(Time = time_vec, do.call(cbind, distance_list))
angle_df    <- data.frame(Time = time_vec, do.call(cbind, angle_list))

# ------------------------------
# STACKING TABLE: angle <30 & distance <0.4
# ------------------------------
strip_suffix <- function(name) {
  sub("_(Distance|Angle)$", "", name)
}

distance_keys <- sapply(names(distance_list), strip_suffix)
angle_keys    <- sapply(names(angle_list), strip_suffix)

common_keys <- intersect(distance_keys, angle_keys)

# Now compute stacking only on matched keys
stack_list <- list()

for (key in common_keys) {
  dist_col <- paste0(key, "_Distance")
  angle_col <- paste0(key, "_Angle")
  
  if (dist_col %in% colnames(distance_df) && angle_col %in% colnames(angle_df)) {
    stack <- (angle_df[[angle_col]] < 30) & (distance_df[[dist_col]] < 0.4)
    stack_list[[key]] <- stack
  } else {
    warning("Missing angle or distance column for: ", key)
  }
}

# Combine into final stacking dataframe
stacking_df <- data.frame(Time = time_vec, do.call(cbind, stack_list))
# ------------------------------
# WRITE OUTPUTS
# ------------------------------
write.csv(distance_df, "distance_table.csv", row.names = FALSE)
write.csv(angle_df,    "angle_table.csv",    row.names = FALSE)
write.csv(stacking_df, "stacking_table.csv", row.names = FALSE)

message("Done! Files written: distance_table.csv, angle_table.csv, stacking_table.csv")

# ------------------------------
# ANGLE STATISTICS WHEN DISTANCE < 0.4
# ------------------------------

# ANGLE STATS WHEN DISTANCE < 0.4
angle_stats <- data.frame(
  Sample = character(),
  Mean   = numeric(),
  Median = numeric(),
  SD     = numeric(),
  Count  = integer(),
  stringsAsFactors = FALSE
)

# Extract base names by stripping _Distance and _Angle
distance_names <- names(distance_df)[-1]  # exclude Time
angle_names    <- names(angle_df)[-1]

get_base_name <- function(x) sub("_(Distance|Angle)$", "", x)

distance_bases <- setNames(distance_names, sapply(distance_names, get_base_name))
angle_bases    <- setNames(angle_names,    sapply(angle_names,    get_base_name))

# Find shared base names (i.e. A17_Thiazole_PKZ18_Rep1)
common_bases <- intersect(names(distance_bases), names(angle_bases))

for (base in common_bases) {
  dist_col <- distance_bases[[base]]
  angle_col <- angle_bases[[base]]
  
  dist_vals  <- as.numeric(distance_df[[dist_col]])
  angle_vals <- as.numeric(angle_df[[angle_col]])
  
  # Get valid indices where distance < 0.4 and both values are finite
  idx <- which(dist_vals < 0.4 & is.finite(dist_vals) & is.finite(angle_vals) & angle_vals>30)
  
  if (length(idx) > 0) {
    close_angles <- angle_vals[idx]
    angle_stats <- rbind(angle_stats, data.frame(
      Sample = base,
      Mean   = mean(close_angles),
      Median = median(close_angles),
      SD     = sd(close_angles),
      Count  = length(close_angles),
      stringsAsFactors = FALSE
    ))
  } else {
    angle_stats <- rbind(angle_stats, data.frame(
      Sample = base,
      Mean   = NA,
      Median = NA,
      SD     = NA,
      Count  = 0,
      stringsAsFactors = FALSE
    ))
  }
}

################################
#      DISTANCE  STATS         #
################################


# DISTANCE STATS WHEN ANGLE < 30
distance_stats <- data.frame(
  Sample = character(),
  Mean   = numeric(),
  Median = numeric(),
  SD     = numeric(),
  Count  = integer(),
  stringsAsFactors = FALSE
)

# Use same base name matching method
for (base in common_bases) {
  dist_col  <- distance_bases[[base]]
  angle_col <- angle_bases[[base]]
  
  dist_vals  <- as.numeric(distance_df[[dist_col]])
  angle_vals <- as.numeric(angle_df[[angle_col]])
  
  # Get valid indices where angle < 30 and both values are finite
  idx <- which(angle_vals < 30 & is.finite(angle_vals) & is.finite(dist_vals))
  
  if (length(idx) > 0) {
    close_dists <- dist_vals[idx]
    distance_stats <- rbind(distance_stats, data.frame(
      Sample = base,
      Mean   = mean(close_dists),
      Median = median(close_dists),
      SD     = sd(close_dists),
      Count  = length(close_dists),
      stringsAsFactors = FALSE
    ))
  } else {
    distance_stats <- rbind(distance_stats, data.frame(
      Sample = base,
      Mean   = NA,
      Median = NA,
      SD     = NA,
      Count  = 0,
      stringsAsFactors = FALSE
    ))
  }
}



write.csv(distance_stats, "distance_stats.csv", row.names = FALSE)

write.csv(angle_stats, "angle_stats.csv", row.names = FALSE)

library(ggplot2)

library(dplyr)
library(tidyr)
library(stringr)
library(scales)

# Convert wide to long for each dataframe
angle_long <- pivot_longer(angle_df, -Time, names_to = "Sample", values_to = "Angle")
angle_long <- angle_long %>%
  mutate(Sample = sub("_Angle$", "", Sample))

distance_long <- pivot_longer(distance_df, -Time, names_to = "Sample", values_to = "Distance")
distance_long <- distance_long %>%
  mutate(Sample = sub("_Distance$", "", Sample))

stacking_long <- pivot_longer(stacking_df, -Time, names_to = "Sample", values_to = "Stacked")

# Merge all
combined_df <- angle_long %>%
  inner_join(distance_long, by = c("Time", "Sample")) %>%
  inner_join(stacking_long, by = c("Time", "Sample"))

# Parse sample info
combined_df <- combined_df %>%
  mutate(
    Stacked = as.factor(Stacked),
    # Extract parts from sample name
    Base     = word(Sample, 2, sep = "_"),
    Replicate = str_extract(Sample, "Rep\\d"),
    GroupA   = word(Sample, 6, sep = "_"),
    GroupB   = word(Sample, 7, sep = "_"),
    SitePair = paste(GroupA, GroupB, sep = "_")
  )

df_long <- combined_df %>%
  filter(grepl("Base", Sample)) %>%
  mutate(
    AngleCloseness = pmin(1, 30 / pmax(Angle, 1)),         # closer to 1 is better
    DistanceCloseness = pmin(1, 0.4 / pmax(Distance, 1e-6)) # closer to 1 is better
  ) %>%
  mutate(
    AngleColor = scales::col_numeric(c("#351288", "#fdfc6b"), domain = c(0, 1))(AngleCloseness),
    DistanceColor = scales::col_numeric(c("#351288", "#fdfc6b"), domain = c(0, 1))(DistanceCloseness)
  )

df_long$Time = (df_long$Time -5000)/1000


stack_lines <- df_long %>%
  filter(as.logical(Stacked)) %>%
  distinct(Time, SitePair, Replicate)

p <- ggplot(df_long, aes(x = Time)) +
  # Angle line
  geom_line(aes(x = Time, y = 1, color = AngleColor), linewidth = 2) +
  # Distance line
  geom_line(aes(x = Time, y = 0, color = DistanceColor), linewidth = 2) +
  # Stacking vertical lines
  geom_vline(data = stack_lines, aes(xintercept = Time), color = "orange", linewidth = 0.02) +
  facet_grid(SitePair ~ Replicate) +
  scale_color_identity() +
  theme_minimal() +
  theme(panel.spacing = unit(1, "lines"))+
  scale_y_continuous(
    breaks = c(0, 1),
    labels = c("Distance", "Angle"),
    limits = c(-0.5, 1.5),
    name = NULL
  ) +
  labs(
    title = "Angle and Distance Over Time",
    subtitle = "Color from purple (far) to yellow (meets threshold); vertical lines = stacking",
    x = "Time"
  )

ggsave("Base_stacking.png", plot = p, width = 10, height = 8, units = "in", dpi=600)


