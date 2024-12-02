#==============================================================================
# PREPARE sgRNA LIBRARY FILE FOR MAGeCK
#
# Giovani Luiz Genesi, MSc
# Phillips Lab @ UPenn
# 2023-06-08  
#
# This script will import your sgRNA library from an excel spreadsheet, a .csv
# file, or a .txt file, and convert it into the required format for MAGeCK
# 
#==============================================================================

if(!require(readxl)){
  install.packages("readxl")
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}



### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
path_to_library = file.path("./sgRNA-Human_Epi_Library.xlsx")
ctrl = "Control" # any word or pattern within your Non-targeting Control sgRNAs
###

# Load your library file
if (grepl("\\.xlsx$", path_to_library)) {
  library(readxl)
  sgRNA_library <- read_excel(path_to_library)
} else if (grepl("\\.csv$", path_to_library)) {
  sgRNA_library <- read.csv(path_to_library)
} else if (grepl("\\.txt$", path_to_library)) {
  sgRNA_library <- read.table(path_to_library)
} else {
  stop("Unsupported file format")
}

# Each line is a sgRNA. You can see there are 4 sgRNAs per gene
sgRNA_library %>% head()

# Label each guide RNA according to its target gene
### Create a vector of target genes
symbols <- sgRNA_library$`Target Gene Symbol`

### Count the occurrences of each gene
counts <- table(symbols)  

### Create a vector of sgRNA ids
sgRNA_id <- vector()

### Loop through each gene
id <- 1
for (each_gene in seq_along(symbols)) {
  if (counts[symbols[each_gene]] > 0) {
    # Add label
    sgRNA_id[each_gene] <- paste(symbols[each_gene], id, sep = "_")
    counts[symbols[each_gene]] <- counts[symbols[each_gene]] - 1
    # Reset id to 1 for each new gene
    if (counts[symbols[each_gene]] == 0) {
      id <- 1  
    } else {
      id <- id + 1
    }
  }
}

# It will look like this
sgRNA_id[1:12]

# Add column of sgRNA ids
sgRNA_library$sgRNA_id = sgRNA_id

# Keep only 3 columns: sgRNA ids, Gene symbols, and target sequences, in this order
sgRNA_library = sgRNA_library[, c("sgRNA_id",
                                  "sgRNA Target Sequence",
                                  "Target Gene Symbol")]

### You can check that each NTC-gRNA was individually labelled
sgRNA_library[grepl(ctrl, sgRNA_library$sgRNA_id),]

# Save as tab-delimited file, excluding the column names and row names
write.table(sgRNA_library, 
            file = "sgRNA_library.txt", 
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)

# Create vector with only the names describing Non-targeting controls
NTCs = sgRNA_library$sgRNA_id[grepl(ctrl,
                                    sgRNA_library$sgRNA_id)]

# It will look like this
NTCs

# Save files describing Non-targeting controls
write.table(data.frame(id = NTCs), 
            file = "nonTargetingControls.txt", 
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)
