if(!require(DESeq2)){
  BiocManager::install("DESeq2")
}
if(!require(showtext)){
install.packages("showtext", dependencies = T)
}
BiocManager::install("EnsDb")
BiocManager::install("EnsDb.Hsapiens.v75")

if(!require(tidytext)){
install.packages("tidytext")
}

library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(rtracklayer)
library(annotate)
library(org.Hs.eg.db)
library(showtext)
library(tidytext)
font_add_google(name = "Roboto", regular.wt = 400, bold.wt = 700)

# turn on showtext
showtext_auto()

#----------
# Prepare peak annotation data
#----------

# Set working directory
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN")

# Import broadPeak files as GRanges object
CBX2_DMSOpeaks <- import("./new_peaks/CBX2_DMSO_final.broadPeak",
                         format = "broadPeak")

CBX2_DMSOpeaks

CBX2_LostPeaks <- import("./new_peaks/CBX2_lost_DMSO_1uM.broadPeak",
                         format = "broadPeak")

CBX2_LostPeaks

CBX2_RetainedPeaks <- import("./new_peaks/CBX2_retained_DMSO_1uM.broadPeak",
                             format = "broadPeak")

CBX2_RetainedPeaks

# Import known genes from hg19 genome
hg19genes <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Annotate peaks
annotatedCBX2 <- annotatePeak(CBX2_DMSOpeaks, 
                              tssRegion = c(-3000,2000), 
                              TxDb = hg19genes, 
                              level = "gene", 
                              addFlankGeneInfo = TRUE, 
                              flankDistance = 0, 
                              overlap = "all")

annotatedLost <- annotatePeak(CBX2_LostPeaks, 
                              tssRegion = c(-3000,2000), 
                              TxDb = hg19genes, 
                              level = "gene", 
                              addFlankGeneInfo = TRUE, 
                              flankDistance = 0, 
                              overlap = "all")

annotatedRetained <- annotatePeak(CBX2_RetainedPeaks, 
                                  tssRegion = c(-3000,2000), 
                                  TxDb = hg19genes, 
                                  level = "gene", 
                                  addFlankGeneInfo = TRUE, 
                                  flankDistance = 0, 
                                  overlap = "all")


# Separate each annotated gene into one per row
annotatedDF <- annotatedCBX2 %>% 
  as.data.frame() %>% 
  separate_rows(19:20, sep = ";") %>% 
  as.data.frame()

annotatedDF %>% head()
annotatedDF %>% nrow()

annotatedLostDF <- annotatedLost %>% 
  as.data.frame() %>% 
  separate_rows(19:20, sep = ";") %>% 
  as.data.frame()

annotatedRetainedDF <- annotatedRetained %>% 
  as.data.frame() %>% 
  separate_rows(19:20, sep = ";") %>% 
  as.data.frame()

# Join all gene Ids into the geneId column.
CBX2_geneId <- merge(annotatedDF[c(6,17)],  
                     annotatedDF[c(6,19)], 
                     by.x = c("name","geneId"), 
                     by.y = c("name","flank_geneIds"), 
                     sort = F, 
                     all = T) %>% unique() %>% filter(geneId != "NA")

Lost_geneId <- merge(annotatedLostDF[c(6,17)],  
                     annotatedLostDF[c(6,19)], 
                     by.x = c("name","geneId"), 
                     by.y = c("name","flank_geneIds"), 
                     sort = F, 
                     all = T) %>% unique() %>% filter(geneId != "NA")

Retained_geneId <- merge(annotatedRetainedDF[c(6,17)],  
                         annotatedRetainedDF[c(6,19)], 
                         by.x = c("name","geneId"), 
                         by.y = c("name","flank_geneIds"), 
                         sort = F, 
                         all = T) %>% unique() %>% filter(geneId != "NA")

# Add gene symbols based on the gene IDs
CBX2_geneId <- CBX2_geneId %>% mutate("SYMBOL" = getSYMBOL(CBX2_geneId$geneId,
                                                           data='org.Hs.eg')) %>%
  filter(SYMBOL != "NA")

CBX2_geneId %>% nrow()

Lost_geneId <- Lost_geneId %>% mutate("SYMBOL" = getSYMBOL(Lost_geneId$geneId,
                                                           data='org.Hs.eg')) %>%
  filter(SYMBOL != "NA")

Retained_geneId <- Retained_geneId %>% mutate("SYMBOL" = getSYMBOL(Retained_geneId$geneId,
                                                                   data='org.Hs.eg')) %>%
  filter(SYMBOL != "NA")

# Get only gene names
CBX2_geneName <- CBX2_geneId["SYMBOL"] %>% unique()
CBX2_geneName %>% nrow()
CBX2_geneName %>% head()

Lost_geneName <- Lost_geneId["SYMBOL"] %>% unique()
Lost_geneName %>% nrow()
Lost_geneName %>% head()

Retained_geneName <- Retained_geneId["SYMBOL"] %>% unique()
Retained_geneName %>% nrow()
Retained_geneName %>% head()

# Filter out common genes among lost and retained peaks
Lost_uniqueGenes <- filter(Lost_geneName, !SYMBOL %in% Retained_geneName$SYMBOL) 
Lost_uniqueGenes %>% nrow()

Retained_uniqueGenes <- filter(Retained_geneName, !SYMBOL %in% Lost_geneName$SYMBOL) 
Retained_uniqueGenes %>% nrow()

#----------
# Prepare gene expression data
#----------
# Set the working directory

setwd("C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2")

# Load raw counts
# Raw counts are used to create DESeqDataSet
counts <- read.delim("cbx2vsLUCCounts.txt", header = T, 
                     row.names = 1)
counts %>% head() # check the loaded count matrix
colnames(counts) # make sure columns correctly correspond to samples' names

gene_length <- read.delim("gene_length.txt", header = F, 
                          row.names = 1, 
                          col.names = c("seqnames", 
                                        "strandness", 
                                        "basepairs"))
gene_length %>% head()

# Create the meta-data data frame with the comparative structure of the data.                     

metaData <- data.frame(row.names = colnames(counts), # samples names
                       treatment = factor(c("DMSO", # grouping
                                            "DMSO",
                                            "1uM",
                                            "1uM")),
                       sgRNA = factor(c("sgLUC","sgLUC","sgLUC","sgLUC",
                                        "sgCBX2.1","sgCBX2.1","sgCBX2.1","sgCBX2.1",
                                        "sgCBX2.3","sgCBX2.3","sgCBX2.3","sgCBX2.3")))

metaData$grouped <- factor(paste0(metaData$sgRNA, metaData$treatment))

metaData 



# Check the factor levels of the condition. The 1st level should be the control.
levels(metaData$sgRNA)
metaData$sgRNA <- relevel(metaData$sgRNA, "sgLUC")
levels(metaData$sgRNA)

levels(metaData$treatment)
metaData$treatment <- relevel(metaData$treatment, "DMSO")
levels(metaData$treatment)


levels(metaData$grouped)
metaData$grouped <- relevel(metaData$grouped, "sgLUCDMSO")
levels(metaData$grouped)

# Create DESeqDataSet

dataSet <- DESeqDataSetFromMatrix(countData = counts, # count matrix
                                  colData = metaData, # samples&grouping
                                  design = ~ grouped)

dataSet

# Add gene lengths as a metadata column named "basepairs"
mcols(dataSet) <- gene_length
mcols(dataSet)

# Pre-filtering genes with 1 read or less
nrow(dataSet)

keep <- rowSums(counts(dataSet)) > 1
dataSet <- dataSet[keep,]
nrow(dataSet)

#------------------------------------------------------------------------------#
#                             GENE EXPRESSION HEATMAP                          #
#------------------------------------------------------------------------------#

#----------
# Prepare data
#----------

# Normalize counts by FPKM
# Good to compare within sample expression, which is what we hope to see in the heatmaps
countsFPKM <- fpkm(dataSet, robust = TRUE) %>% 
  as.data.frame() 
countsFPKM %>% nrow()

# Filter to genes from DE analysis
countsFPKM <- countsFPKM %>% 
  dplyr::filter(row.names(countsFPKM) %in% ezh2iDE$annot.symbol)
countsFPKM %>% nrow()
countsFPKM %>% head()

#----------
# We ultimately want a heatmap where the different samples are shown along the 
# x-axis, the genes are shown along the y-axis, and the shading of the cell 
# reflects how much each gene is expressed within a subject. This latter value, 
# the measure of gene expression, is really just a third dimension. However, 
# instead of creating a 3-dimensional plot that can be difficult to visualize, we 
# instead use shading for our “z-axis”. To this end, we need our data formatted 
# so we have a column corresponding to each of these three dimensions:
#   
#   X: Sample ID
#   Y: Gene symbol
#   Z: Expression
#----------

# We start by transposing our data frame so that we have our samples as rows
# and genes as columns
countsFPKM <- countsFPKM %>%
  t() %>% 
  as.data.frame() 
countsFPKM %>% head()

# Calculate Z-Scores by gene
scaledFPKM <- countsFPKM %>% sapply(FUN = scale) %>% as.data.frame()
scaledFPKM %>% head()

row.names(scaledFPKM) <- row.names(countsFPKM)
scaledFPKM <- scaledFPKM %>% rownames_to_column()
scaledFPKM[,1:12] %>% print(max = 100) 

# # Calculate Log(FPKM)
# logFPKM <- countsFPKM %>% sapply(FUN = log) %>% as.data.frame()
# logFPKM %>% head()
# 
# row.names(logFPKM) <- row.names(countsFPKM)
# logFPKM <- logFPKM %>% rownames_to_column()
# logFPKM[,1:12] %>% print(max = 100) 

longerFPKM <- pivot_longer(data = scaledFPKM,
                           cols = -c(rowname),
                           names_to = "gene",
                           values_to = "Z-Score")

longerFPKM

# Add categories for Retained or Lost peaks
longerFPKM$peakCategory <- "Unnanotated genes"
longerFPKM$peakCategory[longerFPKM$gene %in% Lost_uniqueGenes$SYMBOL] <- "Lost CBX2 peaks" 
longerFPKM$peakCategory[longerFPKM$gene %in% Retained_uniqueGenes$SYMBOL] <- "Retained CBX2 peaks"

longerFPKM$Samples <- ""
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S01_VI_LUCA_DMSO"] <- "sgLUC\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S02_VI_LUCB_DMSO"] <- "sgLUC\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S03_VI_LUCA_1uM"] <- "sgLUC\n1uM EZH2i" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S04_VI_LUCB_1uM"] <- "sgLUC\n1uM EZH2i"
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S13_VI_CBX2.1A_DMSO"] <- "sgCBX2.1\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S14_VI_CBX2.1B_DMSO"] <- "sgCBX2.1\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S15_VI_CBX2.1A_1uM"] <- "sgCBX2.1\n1uM EZH2i" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S16_VI_CBX2.1B_1uM"] <- "sgCBX2.1\n1uM EZH2i" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S17_VI_CBX2.3A_DMSO"] <- "sgCBX2.3\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S18_VI_CBX2.3B_DMSO"] <- "sgCBX2.3\nDMSO" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S19_VVI_CBX2.3A_1uM"] <- "sgCBX2.3\n1uM EZH2i" 
longerFPKM$Samples[longerFPKM$rowname == "Sorted_S20_VI_CBX2.3B_1uM"] <- "sgCBX2.3\n1uM EZH2i" 
longerFPKM$Samples <- factor(longerFPKM$Samples, levels = c("sgLUC\nDMSO", 
                                                    "sgLUC\n1uM EZH2i", 
                                                    "sgCBX2.1\nDMSO", 
                                                    "sgCBX2.1\n1uM EZH2i", 
                                                    "sgCBX2.3\nDMSO", 
                                                    "sgCBX2.3\n1uM EZH2i"))
longerFPKM
levels(longerFPKM$Samples)

range(longerFPKM$`Z-Score`)
#longerFPKM <- longerFPKM %>% dplyr::filter(gene %in% ezh2iDE$SYMBOL)

# Start plot

ggplot(data = longerFPKM, mapping = aes(x = rowname,
                                        y = gene,
                                        #y = reorder_within(gene, by = expression,within = peakCategory),
                                        fill = `Z-Score`)) +
  geom_tile() + 
  theme_classic(base_family = "Roboto",
                base_line_size = 1)+
  scale_y_discrete(position = "right") +
  scale_x_discrete(labels = NULL, ) +
  facet_grid(peakCategory ~ Samples,
             # make the y labels be displayed on the right 
             switch = "y", 
             # let y labels be removed if they do not appear in one of the facets 
             scales = "free",
             # let facet size be proportional to amount of data 
             space = "free") +
  scale_fill_gradient2(low = "black", 
                      mid = "white", 
                      high = "purple", 
                      #midpoint = (max(longerFPKM$expression)+min(longerFPKM$expression))/2, 
                      guide = guide_colorbar(frame.colour = "black", 
                                             ticks.colour = "black")) +
  theme(rect = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(face = "bold.italic", 
                                  size = 11),
        strip.text.x = element_text(face = "bold", 
                                    size = 11),
        panel.spacing.x = unit(0, "cm"),
        panel.grid = element_blank(),
        axis.title = element_blank(), # Remove axes' titles
        axis.line.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.y = element_text(size = 8)) 
  




####################################################
#### TEST HEATMAP FUNCTIONS
####################################################

countsFPKM[1:5]
newTest <- countsFPKM[1:4,]
newTest[1:5]

newTest <- newTest %>% column_to_rownames("Sample") %>% t() %>% as.data.frame()
newTest %>% head()
  
  #mutate('Z-Score' = )

#heatmap(countsFPKM[,1:4], labCol = c("DMSO A", "DMSO B", "1uM A", "1uM B"))


testFPKM <- countsFPKM[1:4,] %>% column_to_rownames("Sample") %>% t() %>% as.data.frame() %>% rownames_to_column()  
testFPKM$peakCategory[testFPKM$rowname %in% Lost_uniqueGenes$SYMBOL] <- "Lost CBX2 peaks"
testFPKM$peakCategory[testFPKM$rowname %in% Retained_uniqueGenes$SYMBOL] <- "Retained CBX2 peaks"
testFPKM %>% head()
testFPKM <- group_by(testFPKM, peakCategory) %>% arrange(`1uM EZH2i A`, .by_group = T)
#testFPKM <- arrange(testFPKM, `1uM EZH2i A` )
testFPKM <- testFPKM %>% column_to_rownames("rowname")
testFPKM %>% head()

# heatmap(as.matrix(testFPKM), 
#         Rowv = NA, 
#         na.rm = T, col=rev(brewer.pal(9,"RdBu")), 
#         #reorderfun = reorder(row.names(testFPKM), testFPKM$`1uM EZH2i A`)
#         )


heatmap(as.matrix(
  dplyr::filter(
    testFPKM, 
    peakCategory != "NA" )[c(1,3)]
    ), 
        Rowv = NULL, 
  Colv = NA,
        na.rm = T, 
        col=rev(brewer.pal(9,"RdBu")))



library(gplots) 
########## TEST ###########
heatmap.2(as.matrix(testFPKM[1:4]), col=rev(brewer.pal(9,"RdBu")), scale="row", Rowv = NA, dendrogram = "none")

#########################################################
# Normalize counts by FPKM
# Good to compare within sample expression, which is what we hope to see in the heatmaps
testFPKM <- fpkm(dataSet, robust = TRUE) %>% 
  as.data.frame() %>% rownames_to_column() 
testFPKM %>% head()

testFPKM$peakCategory <- "Unnanotated genes"

testFPKM$peakCategory[testFPKM$rowname %in% Lost_uniqueGenes$SYMBOL] <- "Lost CBX2 peaks" 
testFPKM$peakCategory[testFPKM$rowname %in% Retained_uniqueGenes$SYMBOL] <- "Retained CBX2 peaks"

testFPKM <- testFPKM %>% dplyr::filter(testFPKM$rowname %in% ezh2iDE$SYMBOL)
testFPKM %>% nrow()

testFPKM %>% head()
testFPKM <- group_by(testFPKM, peakCategory) %>% 
  arrange(.by_group = T) %>% as.data.frame()
testFPKM[c(1,11:14)] %>% head()
row.names(testFPKM) = testFPKM$rowname
testFPKM[c(1,11:14)] %>% head()

if (!require("NMF")) {
  install.packages("NMF", dependencies = TRUE)
  library(NMF)
}

# Scaled within rows, not clustering columns
aheatmap(as.matrix(testFPKM[2:5]),
         color = rev(brewer.pal(9,"RdBu")), 
         scale="row", 
         annColors = "Set2", 
         annRow=testFPKM$peakCategory, 
         Rowv = NA, Colv = NA, labCol = NA )
