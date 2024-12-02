#==============================================================================
# PREPARING SHELL SCRIPT FOR CRISPR DEMULTIPLEXING
#
# Giovani Luiz Genesi, MSc
# Phillips Lab @ UPenn
# 2023-06-08  
#
# This script will import an excel spreadsheet, a .csv file, or a .txt file,
# containing sample barcode information and convert it into a shell script that 
# can demultiplex a fastq file into multiple files, one fastq for each sample
#
# INSTRUCTIONS: The script assumes the barcodes are located 5' (upstream) of 
#               your sgRNA. The script decompresse a fastq.gz file and extracts 
#               sequences from the decompressed fastq.
#               It then deletes the decompressed file
#
#               Run inside your HPC terminal using: 
#                   
#                   sh demultiplex.sh
#
#==============================================================================

if(!require(readxl)){
  install.packages("readxl")
}

library(dplyr)

### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CUT&RUN_CBX2_2023/")
original_fastq = "Philips_Giovani_CUT_RUN_CRISPR_pool_S1_L001_R1_001.fastq"
decompressed_fastq = "EZH2i_screen.fastq"
path_to_barcodes = file.path("./sample_barcodes.csv")
# Specify column indexes (i.e., if your sample names are in the 2nd column, make samples_idx = 2)
barcodes_idx = 4
samples_idx = 1
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

# Start shell script
write.table("#!/bin/sh", 
            file = "demultiplex.sh", # must end with .sh
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.table(paste0("echo R script by Giovani Luiz Genesi, MSc. Phillips Lab, 2023\necho Based on extract_fastq.sh script by Qinglan Li, 2021"), 
            file = "demultiplex.sh",
            append = TRUE, 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "")

# Add command to decompress fastq.gz
write.table(paste0("gunzip -c ", 
                   original_fastq, " > ", decompressed_fastq, 
                   " | echo Decompressing ", 
                   original_fastq, " into ", decompressed_fastq),
            file = "demultiplex.sh",
            append = TRUE,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "")

# Create script tibble
script = barcodes

# Add grep command to extract sequences containing barcodes
script$grep = "grep -A 2 -B 1"

# Specify barcode between ' '
script$barcodes = paste0("'",script$Barcode,"'")

# Add name of decompressed fastq file.
script$fastq = paste(decompressed_fastq, 
                     # Delete the "--" separators introduced by grep
                     "| sed '/^--$/d' >")

# If sample names start with numbers, add an "s" before the number
script$samples = gsub("^([0-9])", "s\\1", script$`sample names`)

# Add file names, replacing spaces by underlines
script$samples = gsub(" ", "_", script$samples)

# Add file extension
script$samples = paste0(script$samples, ".fastq")

# Prevent trailing new line after file name
script$samples = paste0(script$samples, " | echo Demultiplexing ", script$samples, "...")

# Check tibble layout
script = script[c("grep","barcodes","fastq","samples")]
script

# Append to shell script
write.table(script,
            file = "demultiplex.sh", 
            append = TRUE, 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = " ")

# Add comand to delete decompressed fastq
write.table(paste0("rm ", decompressed_fastq, " && echo Removing ", decompressed_fastq), 
            file = "demultiplex.sh",
            append = TRUE, 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "")

write.table(paste0("echo Done."), 
            file = "demultiplex.sh",
            append = TRUE, 
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "")



