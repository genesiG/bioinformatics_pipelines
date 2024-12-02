library(MAGeCKFlute)

setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/mageck_test_edgeRcounts/")

gene_summary = ReadRRA("./GD2_48h_v_CD19_48h_edgeR.gene_summary.txt")
#gene_summary %>% head()
gene_summary$FDR %>% summary()
gene_summary %>% dplyr::arrange(FDR) %>% head()

sgrna_summary = ReadsgRRA("./CD19_48h_v_Baseline_edgeR.sgrna_summary.txt")
#sgrna_summary %>% head()
sgrna_summary$FDR %>% summary()
sgrna_summary %>% dplyr::arrange(FDR) %>% head()

setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/mageck_test_results/")

gene_summary = ReadRRA("./GD2_48h_v_Baseline.gene_summary.txt")
#gene_summary %>% head()
gene_summary$FDR %>% summary()
gene_summary %>% dplyr::arrange(FDR) %>% head()


sgrna_summary = ReadsgRRA("./GD2_48h_v_Baseline.sgrna_summary.txt")
#sgrna_summary %>% head()
sgrna_summary$FDR %>% summary()
sgrna_summary %>% dplyr::arrange(FDR) %>% head()

setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/edgeR/")

resultsedgeR = read.delim("./GD2_48h_vs_Baseline_topTags.txt")
#resultsedgeR %>% head()
resultsedgeR$FDR %>% summary()
resultsedgeR %>% dplyr::arrange(FDR) %>% head()





#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/edgeR")
path_to_count_table = file.path("./CART_counts_filtered_final.rds")
#===

#===
# LOAD COUNTS
#===

# Load DGE List containing reads extrated with processAmplicon()
counts <- readRDS(path_to_count_table)


df1 = as.data.frame(counts$counts[, c("GD2_48h_B")])
df2 = as.data.frame(counts$counts[, c("Baseline_A","Baseline_C")])

ggplot() +
  geom_bar(data = df1, aes(x = row.names(df1), y = rowSums(cpm(df1))),
           fill = "orange", alpha = 0.5, stat = "identity") +
  geom_bar(data = df2, aes(x = row.names(df2), y = rowSums(cpm(df2))),
           fill = "darkgray", alpha = 0.5, stat = "identity") +
  labs(x = "sgRNA", y = "Counts per guide", title = "Baseline A") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  theme(axis.text.x = element_blank())

ggplot() +
  geom_density(data = df1, aes(x = rowSums(cpm(df1))), fill = "orange", alpha = 0.5) +
  geom_density(data = df2, aes(x = rowSums(cpm(df2))), fill = "gray", alpha = 0.5) +
  labs(x = "Counts per guide", y = "Density", title = "Density Comparison") +
  theme(axis.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) + theme_minimal()








# Convert columns to data frames
df1 = as.data.frame(counts$counts[, "GD2_48h_B"])
df2 = as.data.frame(counts$counts[, "Baseline_A"])
df3 = as.data.frame(counts$counts[, "Baseline_C"])

df1 = cpm(df1)
df2 = cpm(df2)
df3 = cpm(df3)

# Create a data frame with all the columns except "Baseline_A" and "Baseline_C"
df_all = as.data.frame(counts$counts[, !(colnames(counts$counts) %in% c("Baseline_A", "Baseline_C"))])
df_all = cpm(df_all)
# Reshape the data frames and assign column names
df1_long <- reshape2::melt(df1)
colnames(df1_long) <- c("Var1", "Var2", "value")

df2_long <- reshape2::melt(df2)
colnames(df2_long) <- c("Var1", "Var2", "value")

df3_long <- reshape2::melt(df3)
colnames(df3_long) <- c("Var1", "Var2", "value")

df1 <- as.data.frame(cpm(counts$counts[, !(colnames(counts$counts) %in% c("Baseline_A", "Baseline_C"))]))
df2 <- as.data.frame(cpm(counts$counts[, "Baseline_A"]))
df3 <- as.data.frame(cpm(counts$counts[, "Baseline_C"]))

library(tidyr)
df <- as.data.frame(cpm(counts$counts))

# Convert df to long format and create a new column based on column names
df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "column") %>%
  mutate(value_type = case_when(
    grepl("^Baseline_A", column) ~ "Baseline_A",
    grepl("^Baseline_C", column) ~ "Baseline_C",
    TRUE ~ NA_character_
  ))


# Plot density plot with facets
ggplot(df_long, aes(x = value, fill = column)) +
  geom_density(data = df_long[df_long$column == "Baseline_C"], 
               alpha = 0.5) +
  geom_density(data = df_long[df_long$column != "Baseline_A"], 
               alpha = 0.5) +
  geom_density(data = df_long[df_long$column != "Baseline_C" & df_long$column != "Baseline_A"], 
               alpha = 0.5) +
  facet_wrap(~ column, scales = "free") +
  labs(x = "Value", y = "Density") +
  theme_bw()





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

#===
# PLOTTING LogFoldChange BETWEEN CONDITIONS/TIMEPOINTS
#===

#===
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
# Specify path to vector of sgRNA ids
path_to_purchased_sgRNA_list = file.path("./selected_guides.rds")

# Specify paths to results tables
path_to_results_table1 = file.path("./edgeR/GD2_48h_vs_Baseline.txt")
path_to_results_table2 = file.path("./edgeR/GD2_72h_vs_Baseline.txt")
#===

tab1 = read.table(path_to_results_table1)
tab2 = read.table(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

tab_merge = merge(tab1, tab2, by = "ID", suffixes = c(".48h",".72h"))
tab_merge %>% head()

tab_merge$enrichment = NA
tab_merge$enrichment[tab_merge$logFC.48h >= 0 & tab_merge$logFC.72h >= 0] = "enrich.all"
tab_merge$enrichment[tab_merge$logFC.48h >= 0 & tab_merge$logFC.72h < 0 ] = "enrich.48h"
tab_merge$enrichment[tab_merge$logFC.48h < 0 & tab_merge$logFC.72h >= 0] = "enrich.72h"
tab_merge$enrichment[tab_merge$logFC.48h < 0 & tab_merge$logFC.72h < 0] = "dep.all"
tab_merge$enrichment[tab_merge$logFC.48h < 0 & tab_merge$logFC.72h >= 0 ] = "dep.48h"
tab_merge$enrichment[tab_merge$logFC.48h >= 0 & tab_merge$logFC.48h < 0] = "dep.72h"

# Filter results table if list of purchased sgRNAs is defined
if (exists("path_to_purchased_sgRNA_list") && file.exists(path_to_purchased_sgRNA_list)) {
  message("Filtering results table to include only purchased sgRNAs")
  guides <- readRDS(path_to_purchased_sgRNA_list)
  tab_filt <- filter(tab_merge, ID %in% guides)
} else {
  warning("Couldn't detect list of purchased sgRNAs. Using all sgRNAs instead.", call. = FALSE)
  guides <- tab1$ID
  tab_filt <- filter(tab_merge, ID %in% guides)
}


# Create the scatter plot of LFC vs LFC
library(ggplot2)
library(ggrepel)

ggplot(tab_filt, aes(x = logFC.48h, 
                    y = logFC.72h, 
                    color = enrichment)
) +
  geom_point(aes(alpha = 0.75)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "48h treatment", y = "72h treatment") +
  scale_color_manual(
    values = c("gray", "darkgreen", "gray", "purple"),
  ) +
  ggtitle("Selected guides", 
          subtitle = "LogFC(GD2 vs Baseline)"
  ) + 
  geom_text_repel(aes(label = ID), 
                  max.overlaps = 15, size = 3.5, 
                  
  ) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, alpha = FALSE)

#===
# NORMALIZED COUNTS PER GUIDE
#===

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
    stop("Oh oh! None of the sgRNA ids provided match the ones in the results table!\n  Did you provide the correct id for the sgRNA of interest?")
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
  message("Filtering results table to include only purchased sgRNAs")
  
} else {
  
  stop("sgRNA ids needed but none provided. Please provide at least one.", call. = FALSE)
  
}

# Plot and save the filtered data frames as individual plots in a single PDF file
pdf("counts_per_guide.pdf", width = pdf_width, height = pdf_height) 

# Iterate over each filtered data frame, plot the counts per guide, and save the plot
for (row_to_plot in row_names) {
  df_filtered <- filtered_list[[row_to_plot]]
  
  plot <- ggplot(df_filtered, aes(x = Sample, y = Value, fill = Sample)) +
    geom_point(shape = 21, size = 3, alpha = 0.5) +
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
