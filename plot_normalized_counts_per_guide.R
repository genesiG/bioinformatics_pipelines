#=============================#
# NORMALIZED COUNTS PER GUIDE #
#=============================#

#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen")
# sgRNAs of interest

  # Specify file path to list of purchased guides
  path_to_purchased_sgRNA_list = file.path("./selected_guides.rds")
  # Specify vector of sgRNA ids (if provided, will be used instead of purchased guides)
  guides_of_interest = c()

# Specify path to counts table
path_to_count_table = file.path("./edgeR/CART_counts_filtered_final.rds")
# Specify required variables for plotting
samples_to_compare = c("Baseline", "GD2_48h", "GD2_72h")
sample_colors = c("black", "gray", "white")
pdf_width = 4
pdf_height = 4
#===

# Load DGE List containing reads extrated with processAmplicon()
counts <- readRDS(path_to_count_table)

logCounts = cpm(counts$counts, log = TRUE)
logCounts %>% head()

# Load required packages
library(dplyr)
library(tidyverse)
library(ggplot2)

# Create data frame from matrix
df <- as.data.frame(logCounts)

df <- df %>% rownames_to_column(var = "Row")
df %>% head()

# Convert data frame into long format
df_long <- df %>% pivot_longer(cols = -Row, names_to = "Group", values_to = "Value")
df_long

# Extract time point and replicate
df_long <- df_long %>%
  mutate(Sample = sub("_[A-Z]$", "", Group))
df_long <- df_long %>% filter(Sample %in% samples_to_compare)

df_long

# Initialize an empty list to store the filtered data frames
filtered_list <- list()

# Define function to iterate over each row name, filter df_long, and store the filtered data frame in the list
getFilteredList <- function(row_names) {
  if (class(row_names) != "character") {
    message("Object guides_of_interest should be a character vector of sgRNA ids.", call. = FALSE)
  }
  
  if (!all(row_names %in% df_long$Row)) {
    stop("Oh oh! None of the sgRNAs match!\n  Did you provide the correct id for the sgRNA of interest?")
  }
  
  # Initialize an empty list to store the filtered data frames
  filtered_list <- list()
  for (row_to_plot in row_names) {
    df_filtered <- df_long %>% filter(Row == row_to_plot)
    filtered_list[[row_to_plot]] <- df_filtered
  }
  
  return(filtered_list)
}

# Define the vector of sgRNAs to iterate over
if (!is.null(guides_of_interest)) {
  if (class(guides_of_interest) != "character"){
    stop("Object guides_of_interest should be a character vector of sgRNA ids.", call. = F)
    
  } else {
    filtered_list <- getFilteredList(guides_of_interest)
    message("Filtering results table to include only selected sgRNAs")
  }
  
} else if (exists("path_to_purchased_sgRNA_list") && file.exists(path_to_purchased_sgRNA_list)) {
  row_names <- readRDS(path_to_purchased_sgRNA_list)
  filtered_list <- getFilteredList(row_names)
  message("Getting results table for purchased sgRNAs")
  
} else {
  
  stop("sgRNA ids needed but none provided. Please provide at least one.", call. = FALSE)
  
}

# Plot and save the filtered data frames as individual plots in a single PDF file
pdf("counts_per_guide.pdf", width = pdf_width, height = pdf_height) 

# Iterate over each filtered data frame, plot the counts per guide, and save the plot
for (row_to_plot in row_names) {
  df_filtered <- filtered_list[[row_to_plot]]
  
  plot <- ggplot(df_filtered, aes(x = Sample, y = Value, fill = Sample)) +
    geom_point(shape = 21, 
               size = 3, 
               alpha = 0.5,
               position = position_jitter(w=0.1,h=0)) +
    labs(x = "Group", y = "Log2(CPM)", title = row_to_plot) +
    theme_bw() +
    scale_fill_manual(
      values = sample_colors,
    ) +
    theme(text = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.title = element_text(face = "bold"),
          title = element_text(face = "bold")) +
    guides(fill = FALSE, alpha = FALSE)
  
  print(plot)
}

# Close the PDF device
dev.off()  


barplot(colSums(filter(counts$counts, 
                       colnames(counts$counts) == grepl("RP")
                       )
                ), las = 2, main = "Counts per index", col = cols, cex.names = 0.5, cex.axis = 0.8
        )