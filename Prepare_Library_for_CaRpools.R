# PREPARE REFERENCE FASTA FILE FOR CaRpools

# Prepare library file in the fasta format
### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/")
dir.create("./CaRpools")
path_to_library = file.path("./CaRpools/humanEpiLibrary.xlsx")
# Specify any word or pattern within your Non-targeting Control sgRNAs
ctrl = "Control" 
# Specify column indexes (i.e., if your gene symbols are in the 2nd column, make gene_symbol = 2)
gene_id = 1
gene_symbol = 2
target_seq = 3
###


library(dplyr)

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
symbols <- sgRNA_library[gene_symbol]

# Load the biomaRt package
library(biomaRt)

# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert gene symbols to ENSEMBL IDs
mapping_table <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "external_gene_name",
  values = symbols[1],
  mart = ensembl,
  uniqueRows = TRUE
)

colnames(symbols)[1] <- "external_gene_name"

# Check for missing gene symbols
missing_symbols <- setdiff(symbols[1], mapping_table["external_gene_name"])

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

### Count the occurrences of each gene 
counts <- table(converted_symbols$ensembl_gene_id)  

### Create a vector of sgRNA ids
sgRNA_id <- vector()

### Loop through each gene
id <- 1
for (each_gene in seq_along(fasta$ensembl_gene_id)) {
  if (counts[fasta$ensembl_gene_id[each_gene]] > 0) {
    # Add label
    sgRNA_id[each_gene] <- paste(fasta$ensembl_gene_id[each_gene], id, sep = "_")
    counts[fasta$ensembl_gene_id[each_gene]] <- counts[fasta$ensembl_gene_id[each_gene]] - 1
    # Reset id to 1 for each new gene
    if (counts[fasta$ensembl_gene_id[each_gene]] == 0) {
      id <- 1  
    } else {
      id <- id + 1
    }
  }
}

# It will look like this
sgRNA_id[1:12]

# Add column of sgRNA ids
fasta$ids = sgRNA_id

# The sample names will be preceeded by a ">"
fasta$ids = paste0(">", fasta$ids)


# It will look like this
fasta = fasta[c("ids", "sgRNA Target Sequence")]
fasta

# Make list intercalating each sample with its barcode sequence
fasta_list = data.frame(C = c(rbind(fasta$ids, 
                                    fasta$`sgRNA Target Sequence`)))

fasta_list[1:20,1]

# Save as tab-delimited file, excluding the column names and row names, and removing quotation marks around characters
write.table(fasta_list, file = "./CaRpools/sgRNA_library.fa",
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)

