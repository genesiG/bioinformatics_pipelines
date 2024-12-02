# PREPARE READ COUNT FILES FOR CaRpools


# Install required packages
packages = c("biomaRt",
             "tidyverse",
             "seqinr",
             "xlsx",
             "rJava",
             "xlsxjars",
             "stringi",
             "scatterplot3d",
             "MESS",
             "DESeq2",
             "rmarkdown",
             "knitr",
             "VennDiagram",
             "sm")

install = vector()
for (each in packages) {
  if (!require(each)){
    install = append(each, install)
  }
}

BiocManager::install(install)

# Prepare library file in the fasta format
### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
#dir.create("./CaRpools")
path_to_library = file.path("./CaRpools/humanEpiLibrary.xlsx")
path_to_readcounts = file.path("./CART_screen.count.txt")
# Specify any word or pattern within your Non-targeting Control sgRNAs
ctrl = "Control" 
# Specify column indexes (i.e., if your gene symbols are in the 2nd column, make gene_symbol = 2)
gene_id = 1
gene_symbol = 2
target_seq = 3
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

if (!require(caRpools)) {
  BiocManager::install("caRpools")
}

library(caRpools)

counts <- load.file(path_to_readcounts, header = TRUE, sep = "\t")

counts %>% head()

# Each line is a sgRNA. You can see there are 4 sgRNAs per gene
sgRNA_library %>% head()

# Label each guide RNA according to its target gene
### Create a vector of target genes
symbols <- sgRNA_library[gene_symbol]




# Remove missing symbols from data frame before merging
symbols <- symbols %>%
  filter(!external_gene_name %in% missing_symbols$external_gene_name)

# Add ENSEMBL ids to gene library
converted_symbols <- merge(symbols, mapping_table, by = colnames(symbols[1]), all.x = TRUE)

converted_symbols_seq <- merge(symbols, mapping_table, by = colnames(symbols[1]), all.x = TRUE)

# Perform a left join to match the sgRNA sequences with gene IDs
colname = colnames(sgRNA_library[gene_symbol])

# Perform an inner join to match the rows and retain the ensembl_gene_id column
fasta <- sgRNA_library %>%
  inner_join(converted_symbols, by = c("Target Gene Symbol" = "external_gene_name")) %>%
  group_by(`Target Gene Symbol`) %>%
  distinct(.keep_all = TRUE) %>%
  ungroup() %>%
  dplyr::select("sgRNA Target Sequence", "ensembl_gene_id")

# Print the matched data
fasta

library(tidyverse)

# Separate the gene symbol and number using the "." delimiter
counts <- counts %>%
  separate(sgRNA, into = c("gene_symbol", "number"), sep = "\\.", remove = FALSE)

# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert gene symbols to ENSEMBL IDs
mapping_table <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = counts$gene_symbol,
  mart = ensembl,
  uniqueRows = TRUE
)

# Check for missing gene symbols
missing_symbols <- setdiff(counts$gene_symbol, mapping_table["external_gene_name"])

# Remove missing symbols from data frame before merging
counts <- counts %>%
  filter(!gene_symbol %in% missing_symbols)

# Add ENSEMBL ids to gene library
colnames(counts)[2] <- "external_gene_name"
converted_symbols <- merge(counts, mapping_table, by = "external_gene_name", all.y = TRUE, all.x = FALSE)

converted_symbols$ids = paste(converted_symbols$ensembl_gene_id, converted_symbols$number, sep = "_")

converted_symbols %>% head()


###
library(tidyverse)
counts <- load.file(path_to_readcounts, header = TRUE, sep = "\t")
# Separate the gene symbol and number using the "." delimiter
counts <- counts %>%
  separate(sgRNA, into = c("gene_symbol", "number"), sep = "\\.", remove = FALSE)
converted_symbols = counts
converted_symbols$ids = paste(converted_symbols$gene_symbol, converted_symbols$number, sep = "_")

converted_symbols %>% head()

# Export day1 counts
write.table(converted_symbols[c("ids", "Baseline_C")], 
            file = "./Baseline_C_counts.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "\t")

# Export treat counts
write.table(converted_symbols[c("ids", "GD2_48h_B")], 
            file = "./GD2_48h_B_counts.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE, 
            sep = "\t")

# Perform the replacement by joining with the converted_symbols data frame
counts <- counts %>%
  inner_join(mapping_table, by = c("gene_symbol" = "external_gene_name")) %>%
  group_by("sgRNA") %>%
  distinct(.keep_all = TRUE) %>%
  ungroup() 

# Print the updated counts data frame
counts %>% head()



