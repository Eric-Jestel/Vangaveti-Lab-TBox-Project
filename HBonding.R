library(tidyverse)
library(stringr)
library(ggplot2)

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
      str_detect(source_partner, "O[12]P") ~ "OXP", # Phosphate backbone oxygens
      str_detect(source_partner, "O[2345]'") ~ "O'", # Ribose sugar oxygens
      str_detect(source_partner, "-N[1-9]") ~ "N",  # Base nitrogens
      str_detect(source_partner, "-O[1-9]") ~ "O",  # Base oxygens
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


# Function to summarize interactions on a per-atom basis (e.g., "A90-N6")
summarize_by_source_atom <- function(hbond_data) {
  # Step 1: Create a simplified ID for the source atom and sum within replicates
  summed_by_atom <- hbond_data %>%
    mutate(
      # Extract residue (e.g., "A90") and atom (e.g., "N6")
      residue_id = str_extract(source_partner, "[ACGU]\\d+"),
      atom_id = str_extract(source_partner, "[^-]+$"),
      # Combine them into a single identifier
      source_atom_id = paste(residue_id, atom_id, sep = "-")
    ) %>%
    # Filter out any rows where the regex didn't find a match
    filter(!is.na(source_atom_id)) %>%
    # Group by replicate and the new source atom ID
    group_by(replicate, source_atom_id) %>%
    # Sum occupancies for these groups
    summarise(total_occupancy = sum(occupancy), .groups = "drop")
  
  # Step 2: Average the summed occupancies across all replicates
  final_atom_summary <- summed_by_atom %>%
    group_by(source_atom_id) %>%
    summarise(
      avg_occupancy = mean(total_occupancy),
      sd_occupancy = sd(total_occupancy),
      .groups = "drop"
    ) %>%
    # Sort by occupancy for easier viewing
    arrange(desc(avg_occupancy))
  
  return(final_atom_summary)
}

# Function to summarize interactions on a per-residue basis (e.g., "A90")
summarize_by_source_residue <- function(hbond_data) {
  # Step 1: Create a simplified ID for the source residue and sum within replicates
  summed_by_residue <- hbond_data %>%
    mutate(
      # Extract just the residue ID (e.g., "A90")
      source_residue_id = str_extract(source_partner, "[ACGU]\\d+")
    ) %>%
    # Filter out any rows where the regex didn't find a match
    filter(!is.na(source_residue_id)) %>%
    # Group by replicate and the source residue ID
    group_by(replicate, source_residue_id) %>%
    # Sum occupancies for all atoms in that residue
    summarise(total_occupancy = sum(occupancy), .groups = "drop")
  
  # Step 2: Average the summed occupancies across all replicates
  final_residue_summary <- summed_by_residue %>%
    group_by(source_residue_id) %>%
    summarise(
      avg_occupancy = mean(total_occupancy),
      sd_occupancy = sd(total_occupancy),
      .groups = "drop"
    ) %>%
    # Sort by occupancy for easier viewing
    arrange(desc(avg_occupancy))
  
  return(final_residue_summary)
}

# --- Generate and Display Per-Atom Summaries ---
atom_summary_124 <- summarize_by_source_atom(hbond_data_124)
atom_summary_base <- summarize_by_source_atom(hbond_data_base)

cat("### Per-Atom Summary for '124' Condition ###\n")
print(atom_summary_124)
cat("\n### Per-Atom Summary for 'Base' Condition ###\n")
print(atom_summary_base)


# --- Generate and Display Per-Residue Summaries ---
residue_summary_124 <- summarize_by_source_residue(hbond_data_124)
residue_summary_base <- summarize_by_source_residue(hbond_data_base)

cat("\n\n### Per-Residue Summary for '124' Condition ###\n")
print(residue_summary_124)
cat("\n### Per-Residue Summary for 'Base' Condition ###\n")
print(residue_summary_base)
















# This function takes a dataframe, groups by the specific UNK atom,
# sums the occupancy within replicates, and then averages across them.
summarize_unk_occupancy <- function(hbond_data) {
  # Step 1: Sum occupancies for each UNK atom type within each replicate.
  summed_occupancy <- hbond_data %>%
    # Extract the specific atom from the UNK partner string (e.g., "N1", "O", "O3").
    mutate(unk_atom = str_extract(unk_partner, "[^-]+$")) %>%
    # Group by the replicate, the side of interaction, and the extracted UNK atom.
    group_by(replicate, unk_atom) %>%
    # Sum the occupancies for these groups.
    summarise(total_occupancy = sum(occupancy), .groups = "drop")

  # Step 2: Calculate the average and standard deviation of these summed
  # occupancies across all replicates.
  final_summary <- summed_occupancy %>%
    # Group by the interaction side and the specific UNK atom.
    group_by(unk_atom) %>%
    # Calculate the mean and standard deviation.
    summarise(
      avg_occupancy = mean(total_occupancy),
      sd_occupancy = sd(total_occupancy),
      .groups = "drop"
    )

  return(final_summary)
}

# Apply the function to both of your dataframes
final_summary_124 <- summarize_unk_occupancy(hbond_data_124)
final_summary_base <- summarize_unk_occupancy(hbond_data_base)

# --- Display the final results ---
cat("### Final Summary for '124' Condition ###\n")
print(final_summary_124)

cat("\n### Final Summary for 'Base' Condition ###\n")
print(final_summary_base)

print_alter_commands <- function(residue_summary_data) {
  # Remove rows where avg_occupancy is NA
  filtered_data <- residue_summary_data[!is.na(residue_summary_data$avg_occupancy), ]
  
  # Loop through each row and print formatted string
  for (i in seq_len(nrow(filtered_data))) {
    # Extract numeric part of the residue ID
    resid <- gsub("[^0-9]", "", filtered_data$source_residue_id[i])
    b_value <- filtered_data$avg_occupancy[i]
    cat(sprintf("alter resid %s, b=%.3f\n", resid, b_value))
  }
}

print_alter_commands(residue_summary_124)


print_atom_alter_commands <- function(atom_summary_data) {
  # Remove rows where avg_occupancy is NA
  filtered_data <- atom_summary_data[!is.na(atom_summary_data$avg_occupancy), ]
  
  # Loop through each row
  for (i in seq_len(nrow(filtered_data))) {
    # Split source_atom_id into residue ID and atom name
    parts <- unlist(strsplit(filtered_data$source_atom_id[i], "-"))
    resid <- gsub("[^0-9]", "", parts[1])   # Extract number from residue ID
    atom_name <- parts[2]
    b_value <- filtered_data$avg_occupancy[i]
    cat(sprintf("show sphere, resid %s and name %s\n", resid, atom_name))
    cat(sprintf("alter resid %s and name %s, b=%.3f\n", resid, atom_name, b_value))
  }
}

print_atom_alter_commands(atom_summary_base)
print_atom_alter_commands(atom_summary_124)

