#==============================================================================
# CRISPR-SCREEN DOWNSTREAM ANALYSIS WITH MAGeCKFlute
#==============================================================================
if(!require("BiocManager")) {
  install.packages("BiocManager")
} 

if(!require("MAGeCKFlute")) {
  BiocManager::install("MAGeCKFlute")
} 
if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(MAGeCKFlute)
library(clusterProfiler)
library(fgsea)
library(ggplot2)
library(dplyr)

# Change for each projet
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/MAGeCK/")

# ## path to the gene summary file (required)
# gene_summary = file.path("C:/Users/thephillipslab/Documents/Projects/MAGeCK/demo/demo.gene_summary.txt")
# ## path to the sgRNA summary file (optional)
# sgRNA_summary = file.path("C:/Users/thephillipslab/Documents/Projects/MAGeCK/demo/demo.sgrna_summary.txt")

## path to the gene summary file (required)
gene_summary = file.path("./DMSO_samples_unfiltered/DMSO-day45_vs_day1-/DMSO-day45_vs_day1-.gene_summary.txt")
## path to the sgRNA summary file (optional)
# sgRNA_summary = file.path("./EPZ_screen/EPZ_total_d45vd1.sgrna_summary.txt")

# Run FluteRRA with only gene summary file
FluteRRA(gene_summary, 
         proj="DMSO_d45_v_baseline", # Prefix of output, can be your project name
         organism="hsa", # hsa = Homo sapiens, mmu = Mus musculus
         incorporateDepmap = TRUE, # Compare with Depmap data to potentially exclude lethal genes from analysis
         omitEssential = TRUE,
         top = 10,# exclude essential genes from analysis
         outdir = "./MAGeCKFlute")

# Run FluteRRA with both gene summary file and sgRNA summary file
FluteRRA(gene_summary,
         sgRNA_summary, 
         proj="EPZ_d45", 
         organism="hsa",  
         incorporateDepmap = TRUE,  
         omitEssential = TRUE,
         top = 10,
         outdir = "./EPZ_screen/")

#------
# STEP-BY-STEP
#------

countsummary = read.delim("./EPZ_screen.countsummary.txt", 
                          check.names = FALSE)
head(countsummary)

setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/MAGeCK/EPZ_screen/")

# QC results
pdf(file = "GiniIndex_sgRNA.pdf")
BarView(countsummary, x = "Label", y = "GiniIndex",
        ylab = "Gini index", main = "Evenness of sgRNA reads")
dev.off()

countsummary$Missed = log10(countsummary$Zerocounts)

pdf(file = "Missed_sgRNAs.pdf")
BarView(countsummary, x = "Label", y = "Missed", fill = "#394E80",
        ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")
dev.off()

# Read mapping
pdf(file = "MapRatesView.pdf")
MapRatesView(countsummary)
dev.off()

### Or
countsummary$Unmapped = countsummary$Reads - countsummary$Mapped
gg = reshape2::melt(countsummary[, c("Label", "Mapped", "Unmapped")], id.vars = "Label")
gg$variable = factor(gg$variable, levels = c("Unmapped", "Mapped"))
gg = gg[order(gg$Label, gg$variable), ]
p = BarView(gg, x = "Label", y = "value", fill = "variable", 
            position = "stack", xlab = NULL, ylab = "Reads", main = "Map ratio")

pdf(file = "MapRates_BarView.pdf")
p + scale_fill_manual(values = c("#9BC7E9", "#1C6DAB"))
dev.off()

#======
# Downstream analysis of MAGeCK RRA
#======

### Read the required data
## path to the gene summary file (required)
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/MAGeCK/")

## path to the gene summary file (required)
gene_summary = file.path("./EPZd45_median_total.gene_summary.txt")

gdata = ReadRRA(gene_summary)
head(gdata)
summary(gdata$FDR)

## path to the sgRNA summary file (optional)
sgRNA_summary = file.path("./EPZ_control_total_d45_v_d13.sgrna_summary.txt")
sdata = ReadsgRRA(sgRNA_summary)
head(sdata)

## Remove missing or duplicate human genes
idx = duplicated(gdata$id)#|is.na(gdata$id)
gdata = gdata[!idx, ]
depmap_similarity = ResembleDepmap(gdata, symbol = "id", score = "Score")
head(depmap_similarity)
gdata %>% head()

# Omit common essential genes from the data
gdata = OmitCommonEssential(gdata, symbol = "id")
sdata = OmitCommonEssential(sdata, symbol = "Gene")


## Visualization of negative selections and positive selections

gdata$LogFDR = -log10(gdata$FDR)
p1 = ScatterView(gdata, x = "Score", y = "LogFDR", label = "id", 
                 model = "volcano", top = 10)

#pdf("volcanoScatter.pdf")
p1
#dev.off()

#pdf("volcanoView.pdf")
VolcanoView(gdata, x = "Score", y = "LogFDR", label = "id"#, 
            #model = "volcano", top = 10
            )
#dev.off()

##Rank plot

###Rank all the genes based on their scores and label genes in the rank plot.
gdata$Rank = rank(gdata$Score)
p1 = ScatterView(gdata, x = "Rank", y = "Score", label = "id", 
                 top = 10, auto_cut_y = TRUE, ylab = "Log2FC", 
                 groups = c("top", "bottom")
                 )

pdf("rankplotScatter.pdf")
p1
dev.off()

### Label interested hits using parameter toplabels (in ScatterView) and genelist (in RankView).
pdf("rankplotLabeled.pdf")
ScatterView(gdata, x = "Rank", y = "Score", label = "id",
            top = 5, auto_cut_y = TRUE, ylab = "Log2FC", 
            groups = c("top", "bottom"), 
            toplabels = c("Pbrm1", "Arid2", "Brd7")
            )

### You can also use the function RankView to draw the figure.

geneList= gdata$Score
names(geneList) = gdata$id
p2 = RankView(geneList, top = 5, bottom = 10)

#or
# RankView(geneList, top = 0, bottom = 0, genelist = c("Pbrm1", "Arid2", "Brd7")) # label hits of interest

pdf("rankView.pdf")
p2
dev.off()

## Dot plot

###Visualize negative and positive selected genes separately.
gdata$RandomIndex = sample(1:nrow(gdata), nrow(gdata))
gdata = gdata[order(-gdata$Score), ]

gg = gdata[gdata$Score>0, ] # positive selection
p1 = ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "top", top = 5, ylab = "Log2FC")

pdf("dotPlot_LFC.pdf")
p1
dev.off()

gg = gdata[gdata$Score<0, ] # negative selection
p2 = ScatterView(gg, x = "RandomIndex", y = "Score", label = "id",
                 y_cut = CutoffCalling(gdata$Score,2), 
                 groups = "bottom", top = 5, ylab = "Log2FC")

pdf("dotPlot_negSelection.pdf")
p2
dev.off()

## sgRankView - visualize the rank of sgRNAs targeting top selected genes.

p2 = sgRankView(sdata, top = 4, bottom = 4)

pdf("sgRankView.pdf")
p2
dev.off()

## Enrichment analysis

### For more information about functional enrichment analysis in MAGeCKFlute, 
### please read the MAGeCKFlute_enrichment document 

geneList= gdata$Score
names(geneList) = gdata$id
enrich_pos = EnrichAnalyzer(geneList = geneList[geneList>0.5], 
                            method = "HGT", type = "KEGG")
enrich_neg = EnrichAnalyzer(geneList = geneList[geneList< -0.5], 
                            method = "HGT", type = "KEGG")

### Visualization of enrichment results
pdf("enrichView_pos.pdf")
EnrichedView(enrich_pos, mode = 1, top = 5, bottom = 0)
dev.off()

pdf("enrichView_neg.pdf")
EnrichedView(enrich_neg, mode = 1, top = 5, bottom = 0)
dev.off()


# PREPARE sgRNA LIBRARY FILE IN THE CORRECT FORMAT

# Load your library file
### Note that the first line is imported as column names
humanEPZlibrary = read.delim("C:/Users/thephillipslab/Desktop/humanEpiLibrary.txt")

# Each line is a sgRNA. You can see there are 4 sgRNAs per gene
humanEPZlibrary %>% head()

# Label each guide RNA according to its target gene
### Create a vector of target genes
symbols <- humanEPZlibrary$Target.Gene.Symbol

### Count the occurrences of each gene
counts <- table(symbols)  

### Create a vector of sgRNA ids
sgRNA_id <- vector()

### Loop through each gene
id <- 1
for (each_gene in seq_along(symbols)) {
  if (counts[symbols[each_gene]] > 0) {
    # Add label
    sgRNA_id[each_gene] <- paste(symbols[each_gene], id, sep = ".")
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
humanEPZlibrary$sgRNA_id = sgRNA_id

# Keep only 3 columns: sgRNA ids, Gene symbols, and target sequences, in this order
humanEPZlibrary = humanEPZlibrary[, c(5,3,2)]

### You can check that each NTC-gRNA was individually labelled
humanEPZlibrary[1188:1245,]

# Save as tab-delimited file, excluding the column names and row names
write.table(humanEPZlibrary, 
            file = "humanEPZlibrary.txt", 
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)

# Save files describing Non-targeting controls
NTCs = humanEPZlibrary$sgRNA_id[grepl("Control",
                                      humanEPZlibrary$sgRNA_id)]

# It will look like this
NTCs

write.table(data.frame(id = NTCs), 
            file = "nonTargetingControls.txt", 
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)









# Load your samples barcodes file
barcodes = read.delim("C:/Users/thephillipslab/Desktop/sample barcodes.txt")

# With the above code, the first line is imported as column names
barcodes = barcodes %>% na.exclude()
barcodes %>% head()

# The sample names will be preceeded by a ">"
barcodes$samples = paste0(">", barcodes$sample.names)
# Barcodes will be preceeded by a "^" to indicate they are at the 5' end
barcodes$sequences = paste0("^", barcodes$Barcode)

# It will look like this
barcodes[,c("samples","sequences")] %>% head()

# Make list intercalating each sample with its barcode sequence
barcodes_list = data.frame(C = c(rbind(barcodes$samples, 
                                       barcodes$sequences)))

barcodes_list %>% head()

# Save as tab-delimited file, excluding the column names and row names, and removing quotation marks around characters
write.table(barcodes_list, file = "samples_barcodes.fa",
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)








# Load your samples barcodes file
barcodes = read.delim("C:/Users/thephillipslab/Desktop/sample barcodes.txt")

# With the above code, the first line is imported as column names
barcodes = barcodes %>% na.exclude()
barcodes %>% head()

# Keep only two columns
barcodes_list = barcodes[,c("sample.names","Barcode")]
barcodes_list %>% head()

# Save as tab-delimited file, excluding the column names and row names, and removing quotation marks around characters
write.table(barcodes_list, file = "samples_barcodes.txt",
            sep = "\t", 
            col.names = F, 
            row.names = F, 
            quote = F)