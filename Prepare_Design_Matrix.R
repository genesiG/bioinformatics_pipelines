#==============================================================================
# PREPARING DESIGN MATRIX FOR MAGeCK MLE
#
# Giovani Luiz Genesi, MSc
# Phillips Lab @ UPenn
# 2023-06-14  
#
# This script will import an excel spreadsheet, a .csv file, or a .txt file,
# containing sample barcode information from a CRISPR screen and convert it 
# into a design matrix to be input in the MAGeCK mle module
#
# INSTRUCTIONS: Use the exported design_matrix.txt as input in mageck mle
#
#==============================================================================

if(!require(readxl)){
  install.packages("readxl")
}

library(dplyr)

### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
original_fastq = "Phillips_James_20230509_S1_L001_R1_001.fastq.gz"
decompressed_fastq = "CART.fastq"
path_to_barcodes = file.path("./CART_trial_sample_barcodes.xlsx")
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
samples = paste0(gsub(" ", "_", barcodes$`sample names`))

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

sample_labels %>% head()

# Create data frame with sample labels
mat = data.frame(samples = sample_labels)

# Add baseline column
mat$baseline = 1

mat

# Add common column
# Add the "common" column
mat$common <- ifelse(mat$samples == "Baseline", 0, 1)

mat

# Add the "CD19" column
mat$CD19 <- ifelse(grepl("CD19", mat$samples), 1, 0)

# Add the "GD2" column
mat$GD2 <- ifelse(grepl("GD2", mat$samples), 1, 0)


# Define the late time points
late <- c("14_day", "21_day", "21d")

# Comparing late vs early time points
mat$CD19_latevsearly <- ifelse(grepl("CD19", mat$samples) & grepl(paste(late, collapse = "|"), mat$samples), 1,
                               ifelse(grepl("CD19", mat$samples), -1, 0))
mat$GD2_latevsearly <- ifelse(grepl("GD2", mat$samples) & grepl(paste(late, collapse = "|"), mat$samples), 1,
                               ifelse(grepl("GD2", mat$samples), -1, 0))

# Define the days elements
days <- c("14_day", "7_day", "21_day", "7d", "21d")

# Comparing days vs hours
mat$CD19_daysvshours <- ifelse(grepl("CD19", mat$samples) & grepl(paste(days, collapse = "|"), mat$samples), 1,
                               ifelse(grepl("CD19", mat$samples), -1, 0))
mat$GD2_daysvshours <- ifelse(grepl("GD2", mat$samples) & grepl(paste(days, collapse = "|"), mat$samples), 1,
                              ifelse(grepl("GD2", mat$samples), -1, 0))

# Compare CD19 at the later timepoint vs Baseline
mat$CD19_late <- ifelse(grepl("CD19_21d", mat$samples), 1, 0)

# Compare GD2 at the later timepoint vs Baseline
mat$GD2_late <- ifelse(grepl("GD2_21d", mat$samples), 1, 0)

# Compare CD19 to GD2 at all time points
mat$CD19vsGD2 <- ifelse(grepl("CD19", mat$samples), 1, ifelse(grepl("GD2", mat$samples), -1, 0))

# Compare CD19 to GD2 at each time point separately
mat$CD19vsGD2_21d <- ifelse(mat$samples %in% c("CD19_21d"), 1, ifelse(mat$samples %in% c("GD2_21d"), -1, 0))
mat$CD19vsGD2_48h <- ifelse(mat$samples %in% c("CD19_48h"), 1, ifelse(mat$samples %in% c("GD2_48h"), -1, 0))
mat$CD19vsGD2_72h <- ifelse(mat$samples %in% c("CD19_72h"), 1, ifelse(mat$samples %in% c("GD2_72h"), -1, 0))
mat$CD19vsGD2_7d <- ifelse(mat$samples %in% c("CD19_7d"), 1, ifelse(mat$samples %in% c("GD2_7d"), -1, 0))


mat


# Append to shell script
write.table(mat,
            file = "design_matrix.txt", 
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")

