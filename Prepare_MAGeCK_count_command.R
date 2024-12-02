#==============================================================================
# PREPARING MAGeCK COUNT COMMAND FROM BARBODES FILE
#
# Giovani Luiz Genesi, MSc
# Phillips Lab @ UPenn
# 2023-06-09  
#
# This script will import an excel spreadsheet, a .csv file, or a .txt file,
# containing sample barcode information from a CRISPR screen and convert it 
# into a mageck count command line that extracts read counts from demultiplexed 
# samples from CRISPR screen.
#
# INSTRUCTIONS: The script assumes you are using the bsub_mageck_count script.
#               Change the variables below as needed. Run this R script.
#               Import the mageckCountCommand.txt into your HPC environment.
#               Copy the contents of the mageckCountCommand.txt into the 
#               bsub_mageck_count file under the "# Run code" comment.
#
#               Run inside your HPC terminal using: 
#                   
#                   perl bsub_mageck_count.pl
#
#==============================================================================

if(!require(readxl)){
  install.packages("readxl")
}

library(dplyr)

### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
path_to_barcodes = file.path("./CART_trial_sample_barcodes.xlsx")
project_name = "CART_screen"
dayZeroLabel = "Baseline"
###

# Load sample barcodes
if (grepl("\\.xlsx$", path_to_barcodes)) {
  library(readxl)
  barcodes <- read_excel(path_to_barcodes)
} else if (grepl("\\.csv$", path_to_barcodes)) {
  barcodes <- read.csv(path_to_barcodes)
} else if (grepl("\\.txt$", path_to_barcodes)) {
  barcodes <- read.table(path_to_barcodes)
} else {
  stop("Unsupported file format")
}

barcodes %>% head()

# Get vector of sample names
samples = paste0(gsub(" ", "_", 
                      barcodes$`sample names`))

# If a sample name starts with a number, add an "s" at the start of the name
samples = gsub("^([0-9])", "s\\1", samples)

# Remove "_A", "_B", or "_C" designating replicates
samples_prefix <- gsub("_[ABC]$", "", samples)

samples_prefix %>% head()

# Append "$fastq_dir/" prefix and ".fastq" extension to each sample
samples_files <- paste("$fastq_dir/", samples, ".fastq", sep = "")

# Create a list to store samples with the same prefix
prefix_list <- split(samples_files, samples_prefix)

# Create vector of sample labels
sample_labels = prefix_list %>% names()

# Separate replicates with commas
replicates <- sapply(prefix_list, function(x) paste(x, collapse = ","))

# Separate samples with spaces
fastq_files <- paste(replicates, collapse = " ")

# Print the final string
print(fastq_files)

# Save full mageck count command into a string
command_string <- paste0("mageck count -l $library -n ", 
                         project_name, 
                         " --control-sgrna $NTC --sample-label ",
                         paste(sample_labels, collapse = ","),
                         " --day0-label ",
                         dayZeroLabel,
                         " --fastq ",
                         fastq_files)

# Print command
command_string

# Save the command string to a text file
writeLines(command_string, "mageckCountCommand.txt")

# Copy the contents of the file into the bsub_mageck_count.pl script and run 
# it on the HPC