# OLIGO DESIGN FOR CRISPR
# AUTOMATED TABLE

work_dir <- "C:/Users/thephillipslab/Documents/Projects/MGB2_screen/followup_screen/"
input_file_path <- file.path(work_dir, "oligo_design_crispr_hprt.csv")
output_file_path <- file.path(work_dir, "oligo_design_crispr_hprt_modified.csv")

# Set work directory
setwd(work_dir)

# Load necessary libraries
library(tidyverse)
library(stringi)

# Define helper functions to modify table
# Function to add G at 5' end
add_5prime_G <- function(seq) {
  # Vectorized operation with ifelse to handle multiple sequences
  ifelse(substr(seq, 1, 1) != "G", paste0("G", seq), seq)
}

# Function to create overhang (add 'CACC' before the sequence)
add_overhang <- function(seq) {
  paste0("CACC", seq)
}

# Function to get the reverse complement
reverse_complement <- function(seq) {
  comp_bases <- chartr("ACGT", "TGCA", seq)
  return(stri_reverse(comp_bases))
}


# Function to add overhang to anti-sense strand
add_5prime_overhang_antisense <- function(seq) {
  paste0("AAAC", seq)
}

# Load table with target gRNA sequences
df <- read.csv(input_file_path)

# Add new columns with transformations
df <- df %>%
  mutate(
    `Sense Add 5' G` = add_5prime_G(`Sequence`),
    `Sense Add 5' overhang` = add_overhang(`Sequence`),
    `Anti-sense` = reverse_complement(`Sequence`),
    `Anti-sense Add 5' overhang` = add_5prime_overhang_antisense(`Anti-sense`)
  )

# View the modified dataframe
print(df)


# Save the modified data as CSV if needed
write_csv(df, output_file_path)

