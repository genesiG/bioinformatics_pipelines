#===
# CHANGE FOR EACH PROJECT
#===
# Set working directory
setwd("C:/Users/thephillipslab/Documents/Projects/AAV/Bibawi_G5T_results/")
# Replace with the actual path to your CSV file
path_to_csv_file = file.path("./AAV_inserts.csv")
# Settings
U6_promoter = TRUE
seq_extension = ".fasta"
#===


# Load required libraries
library(stringr)

# Read the CSV file containing sample names and patterns
pattern_data <- read.csv(path_to_csv_file, stringsAsFactors = FALSE)


# Function to search for a pattern in a `seq_extension` file
search_pattern_in_seq <- function(file_path, pattern) {
  lines <- readLines(file_path)
  sequence <- paste(lines[-1], collapse = "")
  matches <- str_locate_all(sequence, pattern)[[1]]
  return(list(matches = matches, sequence = sequence))
}

# Function to classify the insert based on its position
classify_insert <- function(position, sequence) {
  if (position == -1) {
    return("Insert not found")
  } else if (grepl("CACC|CACCG", substr(sequence, position - 9, position - 0))) {
    return("Correct")
  } else {
    return("Check manually")
  }
}

# Initialize the results data frame
results <- data.frame(
  SampleName = character(),
  Pattern = character(),
  InsertFound = character(),
  Position = integer(),
  stringsAsFactors = FALSE
)

# Loop through the rows of the CSV file and search for patterns
for (i in 1:nrow(pattern_data)) {
  sample_name <- pattern_data$sample_name[i]
  pattern <- pattern_data$insert[i]
  seq_file_path <- paste0(sample_name, seq_extension)
  
  if (file.exists(seq_file_path)) {
    search_result <- search_pattern_in_seq(seq_file_path, pattern)
    matches <- search_result$matches
    sequence <- search_result$sequence
    
    position <- ifelse(length(matches) > 0, matches[1], -1)
    insert_found <- classify_insert(position, sequence)
    
    # Append the result to the data frame
    results <- rbind(results, data.frame(
      SampleName = sample_name,
      Pattern = pattern,
      InsertFound = insert_found,
      Position = position
    ))
  } else {
    print(paste("File not found for sample:", sample_name))
  }
}

# Print the results
print(results)
write.csv(results, 
          sub(".csv", "search_results.csv", path_to_csv_file), 
          quote = F )
