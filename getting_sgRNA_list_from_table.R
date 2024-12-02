#===
# GETTING LIST OF PURCHASED sgRNAs FROM ORDER TABLE
#===

#===
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
# Specify path to the table of ordered sgRNAs
path_to_library = file.path("./Updated Guide Orders.xlsx")
# Specify index of the column that contains sgRNA ids
id_col = 1
#===

library(dplyr)

# Load order table
if (grepl("\\.xlsx$", path_to_library)) {
  library(readxl)
  guides <- read_excel(path_to_library)
} else if (grepl("\\.csv$", path_to_library)) {
  guides <- read.csv(path_to_library)
} else if (grepl("\\.txt$", path_to_library)) {
  guides <- read.table(path_to_library)
} else {
  stop("Unsupported file format")
}

# Check if table was loaded properly
guides %>% head() 

# Get vector of sgRNA ids
guides = guides[id_col] %>% na.exclude()

# Export vector object to be used in other scripts
saveRDS(guides, file = "./selected_guides.rds") 