library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpattern)



# Define a function to process a single h-bond data file.
process_hbond_file <- function(filepath) {
  # ... (This function remains unchanged from the previous version)
  path_components <- str_split(filepath, "/")[[1]]
  condition_name <- path_components[length(path_components) - 2]
  replicate_name <- path_components[length(path_components) - 1]
  data <- read.table(filepath, skip = 1, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  processed_data <- data %>%
    mutate(occupancy = as.numeric(str_remove(occupancy, "%")) / 100) %>%
    mutate(condition = condition_name, replicate = replicate_name) %>%
    select(condition, replicate, donor, acceptor, occupancy)
  return(processed_data)
}

# Step 1: Discover all 'hbonds-details_system.dat' files
hbond_files <- list.files(path = ".", pattern = "^hbonds-details_system\\.dat$", recursive = TRUE, full.names = TRUE)

# Step 2: Process all files and combine them
all_hbond_data <- map_dfr(hbond_files, process_hbond_file)

# --- MODIFICATION: Replace U01 with UNK ---
all_hbond_data <- all_hbond_data %>%
  # Replace U01 with UNK
  mutate(across(c(donor, acceptor), ~ str_replace_all(., "U01", "UNK"))) %>%
  # Create canonical partners to handle A-B vs B-A duplicates
  rowwise() %>%
  mutate(
    partner1 = min(donor, acceptor),
    partner2 = max(donor, acceptor)
  ) %>%
  ungroup() %>%
  # Identify UNK partner and the source (nucleotide) partner
  mutate(
    unk_partner = case_when(
      str_detect(partner1, "UNK") ~ partner1,
      str_detect(partner2, "UNK") ~ partner2,
      TRUE ~ "N/A"
    ),
    source_partner = if_else(unk_partner == partner1, partner2, partner1)
  ) %>%
  # Filter out non-UNK interactions
  filter(unk_partner != "N/A") %>%
  # **NEW**: Generate the generalized atom type, combining O1P and O2P into OXP
  mutate(
    generalized_atom = case_when(
      str_detect(source_partner, "O[12]P") ~ "Phosphate", # Phosphate backbone oxygens
      str_detect(source_partner, "O[2345]'") ~ "Backbone", # Ribose sugar oxygens
      str_detect(source_partner, "-N[1-9]") ~ "Nucleobase",  # Base nitrogens
      str_detect(source_partner, "-O[1-9]") ~ "Nucleobase",  # Base oxygens
      TRUE ~ "Other"
    )
  )

all_hbond_data <- all_hbond_data[,c("condition", "replicate", "unk_partner", "source_partner", "generalized_atom", "occupancy")]

get_residue_numbers <- function(participant_string) {
  # Regex: Look for one of [A, C, G, U] followed by one or more digits (\d+).
  # The digits are "captured" for extraction.
  matches <- str_match_all(participant_string, "[ACGU]([1-9]\\d*)\\b")
  if (length(matches[[1]]) > 0) {
    # Return the captured group (the numbers) as a numeric vector
    return(as.numeric(matches[[1]][, 2]))
  } else {
    # Return an empty numeric vector if no standard nucleotide is found
    return(numeric(0))
  }
}


# --- MODIFICATION: Split into '124' and 'Base' dataframes ---
hbond_data_124 <- all_hbond_data %>% filter(str_detect(condition, "124"))
hbond_data_124 <- hbond_data_124[, !(names(hbond_data_124) %in% "condition")]

get_residue_numbers_vec <- Vectorize(get_residue_numbers)

hbond_data_124 <- hbond_data_124 %>%
  mutate(
    interactionSide = if_else(
      get_residue_numbers_vec(source_partner) > 50,
      "CS",
      "NCS"
    )
  )

hbond_data_base <- all_hbond_data %>% filter(str_detect(condition, "Base"))
hbond_data_base <- hbond_data_base[, !(names(hbond_data_base) %in% "condition")]

hbond_data_base <- hbond_data_base %>%
  mutate(
    interactionSide = if_else(
      get_residue_numbers_vec(source_partner) > 50,
      "CS",
      "NCS"
    )
  )

intra_replicate_summary_base <- hbond_data_base %>%
  # Group the data by the unique identifiers for each interaction
  # class within each simulation replicate.
  group_by(replicate, generalized_atom, interactionSide) %>%
  
  # Collapse each group into a single row, calculating the sum of
  # occupancies for that group.
  summarise(
    total_occupancy = sum(occupancy, na.rm = TRUE),
    # Explicitly drop all grouping structures after this summary.
    # This is a best practice for preventing unintended downstream effects.
    .groups = "drop"
  )


intra_replicate_summary_124 <- hbond_data_124 %>%
  # Group the data by the unique identifiers for each interaction
  # class within each simulation replicate.
  group_by(replicate, generalized_atom, interactionSide) %>%
  
  # Collapse each group into a single row, calculating the sum of
  # occupancies for that group.
  summarise(
    total_occupancy = sum(occupancy, na.rm = TRUE),
    # Explicitly drop all grouping structures after this summary.
    # This is a best practice for preventing unintended downstream effects.
    .groups = "drop"
  )


final_summary_base <- intra_replicate_summary_base %>%
  # Re-group the data, this time by the interaction class, to
  # aggregate across replicates.
  group_by(generalized_atom, interactionSide) %>%
  
  # Collapse the per-replicate data into a final summary row with
  # mean, standard deviation, and count statistics.
  summarise(
    # Calculate the mean occupancy across replicates.
    mean_occupancy = mean(total_occupancy, na.rm = TRUE),
    
    # Robustly calculate the standard deviation. If only one replicate
    # observed the interaction (n() == 1), sd() would return NA.
    # We handle this by assigning 0, as there is no variability.
    sd_occupancy = if_else(n() > 1, sd(total_occupancy, na.rm = TRUE), 0.0),
    
    # Count the number of replicates in which this interaction occurred.
    n_replicates = n(),
    
    # Again, explicitly drop grouping for a clean final output.
    .groups = "drop"
  )

final_summary_124 <- intra_replicate_summary_124 %>%
  # Re-group the data, this time by the interaction class, to
  # aggregate across replicates.
  group_by(generalized_atom, interactionSide) %>%
  
  # Collapse the per-replicate data into a final summary row with
  # mean, standard deviation, and count statistics.
  summarise(
    # Calculate the mean occupancy across replicates.
    mean_occupancy = mean(total_occupancy, na.rm = TRUE),
    
    # Robustly calculate the standard deviation. If only one replicate
    # observed the interaction (n() == 1), sd() would return NA.
    # We handle this by assigning 0, as there is no variability.
    sd_occupancy = if_else(n() > 1, sd(total_occupancy, na.rm = TRUE), 0.0),
    
    # Count the number of replicates in which this interaction occurred.
    n_replicates = n(),
    
    # Again, explicitly drop grouping for a clean final output.
    .groups = "drop"
  )

final_summary_124 = final_summary_124 %>% mutate(drugVariant="124")
final_summary_base = final_summary_base %>% mutate(drugVariant="base")

plot_data <- rbind(final_summary_124, final_summary_base)

plot_data_complete <- plot_data %>%
  complete(interactionSide, generalized_atom, drugVariant,
           fill = list(mean_occupancy = .001, sd_occupancy = 0, n_replicates = 4))

flip_adjacent_rows <- function(df) {
  n <- nrow(df)
  # Check if we have an odd number of rows
  if (n %% 2 != 0) {
    warning("Dataframe has odd number of rows. Last row will remain in place.")
  }
  
  # Create a new index vector that swaps adjacent rows
  new_idx <- numeric(n)
  for (i in seq(1, n-1, by=2)) {
    if (i+1 <= n) {
      new_idx[i] <- i+1
      new_idx[i+1] <- i
    } else {
      new_idx[i] <- i  # Handle the last odd row
    }
  }
  
  # Return reordered dataframe
  return(df[new_idx, ])
}

plot_data_complete <- flip_adjacent_rows(plot_data_complete)

plot_data_complete$drugVariant <- factor(plot_data_complete$drugVariant, levels = c("base", "124"))


p = ggplot(plot_data_complete, aes(x = generalized_atom, y = mean_occupancy, fill = n_replicates, group=drugVariant, pattern=drugVariant)) +
  
  # Create bars side-by-side. The completed data ensures no condition is skipped.
  # Create bars side-by-side with patterns
  geom_col_pattern(
    colour = "black",
    position = position_dodge(.7),
    width = 0.6,
    # Set pattern options
    aes(pattern_angle = drugVariant),  # Map angle to drugVariant
    pattern_density = 0.1,
    pattern_spacing = 0.025,
    pattern_colour = "black"
  ) +

  scale_pattern_manual(values = c("base" = "none", "124" = "pch"), 
                     name = "Drug Variant") +
  
  # Define angles for the stripes (135° for left, 45° for right)
  scale_pattern_angle_manual(values = c("base" = 135, "124" = 45),
                           guide = "none") +
  
  # Create separate plots for "CS" and "NCS"
  facet_wrap(~ interactionSide) +
  
  # Use a continuous, color-blind friendly gradient scale (Viridis is an excellent choice)
  # The "c" at the end of "_viridis_c" indicates a continuous scale.
  scale_fill_viridis_c(option="cividis", name = "Replicates") +
  
  # Set the y-axis range from 0 to 1 without removing any data
  coord_cartesian(ylim = c(0, 1)) +
  
  # Add titles and labels for clarity
  labs(
    title = "Hydrogen Bonding",
    x = "RNA Bond Target",
    y = "Mean HBond Occupancy"
  ) +
  
  # Use a clean theme and adjust text for readability
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(face="bold")
  )

p

ggsave("HBonding.png", plot = p, width = 10, height = 8, units = "in", dpi=600)




