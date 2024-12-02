###############################################################################
#                       PEAK ANNOTATION WITH CHIPSEEKER                       #
###############################################################################

# Check if the package is installed by trying to load it with require()
# require() returns a logical (TRUE or FALSE) depending on if it's able to load 
# the package. If it fails, install the package.

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if(!require(annotatr)){
  BiocManager::install("annotatr")
  library(annotatr)
}

# Load required packages
library(annotatr)
library(ChIPseeker)
library(EnsDb.Hsapiens.v75)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(tidytext)
library(ggrepel)
library(ashr)





#----------
# Prepare gene expression data
#----------
# Set the working directory

setwd("C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2")

# Load raw counts
counts <- read.delim("cbx2vsLUCCounts.txt", header = T, row.names = 1)

counts %>% head() # check the loaded count matrix

#counts <- countsAll#[,-c(5:12)]
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
                                        "sgCBX2.3","sgCBX2.3","sgCBX2.3","sgCBX2.3")),
                       genotype = factor(c("WT","WT","WT","WT",
                                        "KO","KO","KO","KO",
                                        "KO","KO","KO","KO")))

metaData$grouped <- factor(paste0(metaData$genotype, metaData$treatment))

metaData 



# Check the factor levels of the condition. The 1st level should be the control.
levels(metaData$sgRNA)
metaData$sgRNA <- relevel(metaData$sgRNA, "sgLUC")
levels(metaData$sgRNA)

levels(metaData$treatment)
metaData$treatment <- relevel(metaData$treatment, "DMSO")
levels(metaData$treatment)

levels(metaData$genotype)
metaData$genotype <- relevel(metaData$genotype, "WT")
levels(metaData$genotype)

levels(metaData$grouped)
metaData$grouped <- relevel(metaData$grouped, "WTDMSO")
levels(metaData$grouped)

# Create DESeqDataSet
dataSet <- DESeqDataSetFromMatrix(countData = counts, # count matrix
                                  colData = metaData, # samples&grouping
                                  design = ~  grouped
                                  #design = ~  treatment + genotype + genotype:treatment
                                  )

dataSet

# Add gene lengths as a metadata column named "basepairs"
mcols(dataSet) <- gene_length
mcols(dataSet)

# Pre-filtering genes with 1 read or less
dataSet <- DESeq2::estimateSizeFactors(dataSet)
dataSet <- DESeq2::estimateDispersions(dataSet)

nrow(dataSet)

keep <- rowSums( counts(dataSet, normalized=TRUE) >= 5 ) >= 2
dataSet <- dataSet[keep,]
nrow(dataSet)

 # keep <- rowSums(counts(dataSet)) > 100
 # dataSet <- dataSet[keep,]
 # nrow(dataSet)


#------------------------------------------------------------------------------#
#                       DIFFERENTIAL EXPRESSION ANALYSIS                       #
#------------------------------------------------------------------------------#

DESeq <- DESeq(dataSet,
               #test="LRT", reduced=~sgRNA + treatment
               )

DESeq2::plotDispEsts(DESeq)
resultsNames(DESeq)

# How each condition affects read counts for CBX2?
plotCounts(DESeq, 
           gene = "CBX2", 
           intgroup = "treatment", 
           main = "CBX2")

plotCounts(DESeq, 
           gene = "CBX2", 
           intgroup = "grouped", 
           main = "CBX2")

plotCounts(DESeq, 
           gene = "CDKN2A", 
           intgroup = "treatment", 
           main = "CDKN2A")

plotCounts(DESeq, 
           gene = "CDKN2A", 
           intgroup = "grouped", 
           main = "CDKN2A")

saveRDS(DESeq,file = "./DESeq.rds")
DESeq <- readRDS("./DESeq.rds")
# Write results to comma-delimited file or tab-delimited text file
# write.csv(as.data.frame(results), 
#           file="results.csv")
# write.table(as.data.frame(results),
#             sep = "\t",
#             file="results.txt")

#----------
annotatr
#----------
library(annotatr)
# Change for each Project
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN")

# Import peaks
CBX2_diff <- GRangesList(CBX2_decreased = import("./csaw/DB/CBX2_decreased.merge10.q005.FC0.6.max10k.bed",
                                                format = "bed"),
                             CBX2_increased = import("./csaw/DB/CBX2_increased.merge10.q005.FC0.6.max10k.bed",
                                               format = "bed"))
CBX2_bdgdiff

CBX2_peakList <- GRangesList(CBX2_retained = import("./csaw/bed/CBX2_retained_1uM.bed",
                                                   format = "bed"), 
                            CBX2_decreased = import("./csaw/bed/CBX2_DB_decreased.filterControl2.0.merge1k.q001.FC.6.max50k.bed",
                                                    format = "bed"), 
                            CBX2_increased = import("./csaw/bed/CBX2_DB_increased.filterControl2.0.merge1k.q001.FC.6.max50k.bed",
                                                    format = "bed"))
CBX2_peakList

K27me3_peaks = import("./csaw/K27me3_DMSO.filterlocal1.0.merge600.q0001.FC4.max50k.bed", format = "bed")

# ChIP peaks coverage plot (from ChIPseeker)

#pdf('rplot.pdf')
#cplot <- covplot(CBX2_peakList, weightCol="score", chrs=c("chr2","chr7","chr18","chr22"), xlim = c(2e7, 2e8)                 ) 
#dev.off()

# # Import broadPeak files as GRanges object
# library(rtracklayer)
# CBX2_DMSOpeaks <- import("./new_peaks/CBX2_DMSO_final.broadPeak",
#                          format = "broadPeak")


# Select annotations for intersection with regions
annots = c( 
           'hg19_basicgenes', 'hg19_genes_intergenic',
           'hg19_cpg_islands',
           'hg19_enhancers_fantom'
           )

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg19', annotations = annots)
annotations

saveRDS(annotations, "./annotation/annotations_hg19.rds")

annotations = readRDS("./annotation/annotations_hg19.rds")

# Intersect the regions we read in with the annotations
increasedCBX2 = annotate_regions(regions = CBX2_peakList$CBX2_increased,
                                annotations = annotations,
                                ignore.strand = TRUE,
                                quiet = FALSE)

# A GRanges object is returned
print(increasedCBX2)

# Summarize annotations
summarize_annotations(
  annotated_regions = increasedCBX2,
  quiet = TRUE)

# Plot heatmap of pairs of annotations
plot_coannotations(annotated_regions = increasedCBX2,
                   axes_label = 'Annotations',
                   plot_title = 'Regions in Pairs of Annotations')

#----------
# The annotate_regions() function returns a GRanges, but it may be more 
# convenient to manipulate a coerced data.frame
#----------

# Coerce to a data.frame
df_annotatrCBX2 = data.frame(increasedCBX2)

# See the GRanges column of the annotation df expanded
print(head(df_annotatrCBX2))

# Repeat for the other peak sets
annotatrDecreased = annotate_regions(regions = CBX2_peakList$CBX2_decreased,
                                annotations = annotations,
                                ignore.strand = TRUE,
                                quiet = FALSE)

annotatrRetained = annotate_regions(regions = CBX2_peakList$CBX2_retained,
                                annotations = annotations,
                                ignore.strand = TRUE,
                                quiet = FALSE)

df_annotatrDecreased = data.frame(annotatrDecreased)
df_annotatrRetained = data.frame(annotatrRetained)


K27peaks = annotate_regions(regions = K27me3_peaks,
                 annotations = annotations,
                 ignore.strand = TRUE,
                 quiet = FALSE)

df_K27me3_peaks = data.frame(K27peaks)

PRC_geneName <- df_K27me3_peaks["annot.symbol"] %>% unique() %>% filter(!"annot.symbol" == "NA")
PRC_geneName %>% nrow() 
PRC_geneName %>% head()
# Get only gene names
Increased_geneName <- df_annotatrCBX2["annot.symbol"] %>% unique() %>% filter(!is.na("annot.symbol"))
Increased_geneName %>% nrow() 
Increased_geneName %>% head()

Lost_geneName <- df_annotatrDecreased["annot.symbol"] %>% unique() %>% filter(!"annot.symbol" == "NA")
Lost_geneName %>% nrow()
Lost_geneName %>% head()

Retained_geneName <- df_annotatrRetained["annot.symbol"] %>% unique() %>% filter(!"annot.symbol" == "NA")
Retained_geneName %>% nrow()
Retained_geneName %>% head()

# Filter out common genes among lost and retained peaks
Lost_uniqueGenes <- filter(Lost_geneName, !annot.symbol %in% Retained_geneName$annot.symbol & 
                             !annot.symbol %in% Increased_geneName$annot.symbol
                           ) 
Lost_uniqueGenes %>% nrow()

Retained_uniqueGenes <- filter(Retained_geneName, #!annot.symbol %in% Lost_geneName$annot.symbol & 
                                 !annot.symbol %in% Increased_geneName$annot.symbol) 
Retained_uniqueGenes %>% nrow()

Increased_uniqueGenes <- Increased_geneName
Increased_uniqueGenes %>% nrow()

CBX2_geneName = rbind(Increased_uniqueGenes, Retained_uniqueGenes, Lost_uniqueGenes)
CBX2_geneName %>% head()
CBX2_geneName %>% nrow()

### Volcano Plot
DESeqAnalysis <- function(resultsDF = resultsDF, 
                          alpha = 0.05,
                          comparison = character(),
                          FC.cutoff = 0.585,
                          colors = c("black", "gray", "red")
){
  
  FC.col = "magenta"
  if(FC.cutoff == 0) {
    FC.col = NA
  }
  
  
  
  # Merge annotated peaks and gene expression data by gene name
  peakExpression <- merge(CBX2_geneName,
                          resultsDF,
                          by.x = "annot.symbol",
                          by.y = "row", 
                          all.y = T,
                          sort = FALSE)
  
  #----------
  # Volcano Plot
  #----------
  # Add column for labels of differentially expressed genes
  peakExpression$Label <- NA
  peakExpression$Label[peakExpression$isDE != "Not DE"] <- peakExpression$annot.symbol[peakExpression$isDE != "Not DE"]
  
  labeledPeak <- peakExpression #%>% dplyr::filter(Peak != "")
  
  print(x = peakExpression[order(peakExpression$pvalue),-c(5,9)],
        max = 100,
        row.names = F)
  
  # Make volcano plot
  volcanoPlot <- ggplot(data=labeledPeak, 
                        aes(x=log2FoldChange, 
                            y=-log10(padj),
                            color = isDE, 
                            shape = Peak, 
                            label =  Label)) + 
    geom_point() + 
    theme_bw() +
    scale_shape_manual(values = c(1,16,4,17)) +
    ylab(expression(bold(bolditalic(-Log[10]) ~ "FDR"))) +
    xlab(expression(bold(bolditalic(Log[2]) ~ "FC"))) +
    ggtitle(comparison) + 
    scale_x_continuous(n.breaks = 7) +
    theme(legend.title = element_blank(), 
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"), 
          panel.grid.minor = element_blank()) 
  
  # We can organize the labels nicely using the "ggrepel" package and the 
  # geom_text_repel() function
  
  print(volcanoPlot + geom_text_repel() +
          # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
          geom_vline(xintercept=c(-FC.cutoff, FC.cutoff), col=FC.col, linetype = 2) +
          geom_hline(yintercept=-log10(alpha), col="black", linetype = 2) +
          # Change point color 
          scale_color_manual(values=c("Downregulated" = colors[1],
                                      "Not DE" = colors[2], 
                                      "Upregulated" = colors[3]))
  )
}

#------------------------------------------------------------------------------#
#                                1uM EZH2i VS DMSO                             #
#------------------------------------------------------------------------------#

# Conditions to compare
comparison = "1uM EZH2i vs DMSO"
# Fold-change cut-off
FC.cutoff = .5 
# q-value (FDR) cut-off
alpha = 0.05
# Contrast
contrast = c("grouped","WT1uM","WTDMSO")

# Building the results table 
results <- results(DESeq, 
                   contrast = contrast,
                   #lfcThreshold = .5, 
                   alpha = alpha) 
results # the comparison should be TREATED vs CONTROL

summary(results)
table(results$padj < alpha) 
table(results$padj < alpha & results$log2FoldChange > FC.cutoff) 

# Get DESeq results as data.frame with rownames as the first column
resultsNames(DESeq)

library(tidyverse)
resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")

resultsDF <- results(DESeq, 
                     #name = "treatment1uM.genotypeKO",
                    contrast = contrast,
                   #lfcThreshold = .5, 
                   alpha = alpha, 
                   tidy = T)
resultsDF %>% head()
resultsDF[order(resultsDF$log2FoldChange, decreasing = T),] %>% head()
resultsDF[order(resultsDF$padj),] %>% head()

# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

write.csv(na.omit(resultsDF), 
          file = "./resultsDF_1uM_vs_DMSO.csv", 
          row.names = F)

# Get ranked gene list for GSEA
resultsDF$rank = resultsDF$log2FoldChange * -1*log10(resultsDF$padj)

resultsRNK = resultsDF %>% filter(isDE != "Not DE")
resultsRNK = resultsRNK[c(1,9)]

resultsRNK = resultsDF[c(1,9)]
resultsRNK$row = factor(resultsRNK$row)
resultsRNK %>% head()

write.table(na.omit(resultsRNK), file = "./resultsRanked.rnk", sep = "\t", row.names = F)
write.csv(na.omit(resultsRNK), 
          file = "./resultsRanked.csv", 
          row.names = F)

resultsLogFC = resultsDF[c(1,3)]
resultsLogFC %>% head()

write.csv(na.omit(resultsLogFC), 
          file = "./resultsLogFC.csv", 
          row.names = F)

resultsDF$Condition = "EZH2i"
# Add lost/retained peaks category
resultsDF$Peak <- ""
resultsDF$Peak[resultsDF$row %in% Retained_uniqueGenes$annot.symbol] <- "Retained CBX2 peak"
resultsDF$Peak[resultsDF$row %in% Lost_uniqueGenes$annot.symbol] <- "Lost CBX2 peak"



# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "annot.symbol",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()

# Plot numerical
df_annotatrCBX2 %>% head()
df_annotatrCBX2 %>% nrow()
annotatrCBX2 <- merge(df_annotatrCBX2, resultsDF[c(1,3,7:9)], by.x = "annot.symbol", by.y = "row", sort = F) %>% 
  filter(isDE != "Not DE", Peak != "") %>%
  GRanges()
annotatrCBX2$padj <- -log10(annotatrCBX2$padj)
annotatrCBX2

resultsDF <- merge(resultsDF, df_annotatrCBX2[c(20:21)], by.y = "annot.symbol", by.x = "row", sort = F) %>% 
  filter(Peak != "") 

plot_numerical(
  annotated_regions = annotatrCBX2,
  x = 'log2FoldChange',
  y = 'padj',
  facet = c('annot.type','Peak'),
  facet_order = list(c(#'hg19_genes_1to5kb',
                       "hg19_genes_promoters"#,
                       #"hg19_genes_exons","hg19_genes_introns", "hg19_genes_5UTRs", "hg19_genes_3UTRs", "hg19_enhancers_fantom"
                       ), 
                  c("Retained CBX2 peak", "Lost CBX2 peak")),
  plot_title = 'Differentially expressed genes annotated to CBX2 peaks (EZH2i)',
  x_label = 'Log2FoldChange',
  y_label = 'Log10(FDR)'
  ) + scale_x_continuous(limits = c(-4,4))

plot_numerical(
  annotated_regions = annotatrCBX2,
  x = 'log2FoldChange',
  y = 'padj',
  facet = c('annot.type'),
  facet_order = "hg19_genes_promoters",
  plot_title = 'Genes annotated to CBX2 peaks (EZH2i)',
  x_label = 'Log2FoldChange',
  y_label = 'Log10(FDR)'
) + scale_x_continuous(limits = c(-4,4))
#----------
# Plot differential expression data with annotated peaks 
#----------
DESeqAnalysis <- function(resultsDF = resultsDF, 
                          alpha = 0.05,
                          comparison = character(),
                          FC.cutoff = 1,
                          colors = c("black", "gray", "red")
                          ){
  
  FC.col = "magenta"
  if(FC.cutoff == 0) {
    FC.col = NA
  }
  
  
  
  # Merge annotated peaks and gene expression data by gene name
  peakExpression <- merge(CBX2_geneName,
                          resultsDF,
                          by.x = "annot.symbol",
                          by.y = "row",
                          sort = FALSE)
  
  #----------
  # Volcano Plot
  #----------
  # Add column for labels of differentially expressed genes
  peakExpression$Label <- NA
  peakExpression$Label[peakExpression$isDE != "Not DE"] <- peakExpression$annot.symbol[peakExpression$isDE != "Not DE"]
  
  labeledPeak <- peakExpression #%>% dplyr::filter(Peak != "")
  
  print(x = peakExpression[order(peakExpression$pvalue),-c(5,9)],
        max = 100,
        row.names = F)
  
  # Make volcano plot
  volcanoPlot <- ggplot(data=labeledPeak, 
                        aes(x=log2FoldChange, 
                            y=-log10(padj),
                            color = isDE, 
                            shape = Peak, 
                            label =  Label)) + 
    geom_point() + 
    theme_bw() +
    scale_shape_manual(values = c(1,16,13)) +
    ylab(expression(bold(bolditalic(-Log[10]) ~ "FDR"))) +
    xlab(expression(bold(bolditalic(Log[2]) ~ "FC"))) +
    ggtitle(comparison) + 
    scale_x_continuous(n.breaks = 7) +
    theme(legend.title = element_blank(), 
          axis.title.x = element_text(face = "bold"), 
          axis.title.y = element_text(face = "bold"), 
          panel.grid.minor = element_blank()) 
  
  # We can organize the labels nicely using the "ggrepel" package and the 
  # geom_text_repel() function
  
  print(volcanoPlot + geom_text_repel() +
          # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
          geom_vline(xintercept=c(-FC.cutoff, FC.cutoff), col=FC.col, linetype = 2) +
          geom_hline(yintercept=-log10(alpha), col="black", linetype = 2) +
          # Change point color 
           scale_color_manual(values=c("Downregulated" = colors[1],
                                       "Not DE" = colors[2], 
                                       "Upregulated" = colors[3]))
        )
}


treatmentResults <- resultsDF
promoterResults <- resultsDF %>% filter(annot.type == "hg19_genes_promoters") %>% unique()
DESeqAnalysis(treatmentResults)
DESeqAnalysis(promoterResults)

#------------------------------------------------------------------------------#
#                         sgCBX2.1 (DMSO) VS sgLUC (DMSO)                        #
#------------------------------------------------------------------------------#

#----------
# Set parameters
#----------
# Conditions to compare
comparison = "sgCBX2.1 (DMSO) vs sgLUC (DMSO)"

# Contrast
contrast = c("grouped","sgCBX2.1DMSO","sgLUCDMSO")

# Get DESeq results as data.frame with rownames as the first column
results <- results(DESeq, 
                   contrast = contrast,
                   alpha = alpha, 
                  # independentFiltering = F
) 

results %>% summary()
resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")


# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

resultsDF %>% filter(row == "CDKN2A")
resultsDF$Condition = factor("CBX2 KO")



# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "SYMBOL",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()
ezh2iDE %>% nrow()

cbx2.1Results <- resultsDF
DESeqAnalysis(cbx2.1Results)

#------------------------------------------------------------------------------#
#                         sgCBX2.3 (DMSO) VS sgLUC (DMSO)                      #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.3 (DMSO) VS sgLUC (DMSO)"

# Contrast
contrast = c("grouped","sgCBX2.3DMSO","sgLUCDMSO")
contrast = c("grouped","KODMSO","WTDMSO")

# Building the results table #
results <- results(DESeq, 
                   contrast = contrast,
                   alpha = alpha, 
                   independentFiltering = F
) 
#results
print(summary(results))
table(results$padj < alpha) 
table(results$padj < alpha & results$log2FoldChange > FC.cutoff) 
# Prepare results data-frame
resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")

# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

resultsDF %>% filter(row == "CDKN2A")
resultsDF$Condition = "CBX2 KO"

# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "SYMBOL",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()
ezh2iDE %>% nrow()

cbx2.3Results <- resultsDF
DESeqAnalysis(cbx2.3Results)

#------------------------------------------------------------------------------#
#                          sgCBX2.1 (1uM) VS sgLUC (1uM)                       #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.1 (1uM EZH2i) VS sgLUC (1uM EZH2i)"
contrast = c("grouped","sgCBX2.11uM","sgLUC1uM")


# Building the results table #

#resultsNames(DESeq)
results <- results(DESeq, 
                   name = "sgRNAsgCBX2.1.treatment1uM",
                   #contrast = contrast,
                   alpha = alpha, 
                   #independentFiltering = F
                   ) 
#results
print(summary(results))
table(results$padj < alpha) 
table(results$padj < alpha & results$log2FoldChange > FC.cutoff) 
table(results$padj < alpha & results$log2FoldChange < -1*FC.cutoff)

resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")


# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

resultsDF %>% filter(row == "CDKN2A")
#results2.1 = resultsDF


resultsDF$Condition = "CBX2 KO + EZH2i / EZH2i"

# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "SYMBOL",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()

cbx2.1_1uM_Results <- resultsDF
DESeqAnalysis(cbx2.1_1uM_Results)

#------------------------------------------------------------------------------#
#                          sgCBX2.3 (1uM) VS sgLUC (1uM)                       #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.3 (1uM EZH2i) VS sgLUC (1uM EZH2i)"
comparison = "CBX2 KO (1uM EZH2i) vs WT (1uM EZH2i)"
# Fold-change cut-off
FC.cutoff = 0
# q-value (FDR) cut-off
alpha = 5e-2
# Contrast
contrast = c("grouped","KO1uM","WT1uM")


# Building the results table #

#resultsNames(DESeq)
results <- results(DESeq, 
                   #name = "sgRNAsgCBX2.3.treatment1uM",
                   contrast = contrast,
                   alpha = alpha,
                   #test = "Wald",
                   #lfcThreshold = 0.5, altHypothesis = "greaterAbs",
                   #independentFiltering = F,
                   ) 
#results
print(summary(results))
#table(results$padj < alpha)
table(results$padj < alpha & results$log2FoldChange > FC.cutoff) 
table(results$padj < alpha & results$log2FoldChange < -1* FC.cutoff) 

resultsDF <- lfcShrink(DESeq, 
          contrast = contrast,
          type = "ashr", 
          res = results
          ) %>% as.data.frame() %>% rownames_to_column("row")

contrast = c("grouped","WT1uM","WTDMSO")
contrast = c("grouped","KO1uM","WTDMSO")
resultsDF <- results(DESeq,
                     #name = "sgRNAsgCBX2.3.treatment1uM",
                     #test = "Wald",
                     contrast = contrast,
                     alpha = alpha,
                     tidy = TRUE)

# Add classes of differential expression
resultsDF$isDE <- "Not DE"

resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & 
                 resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & 
                 resultsDF$padj < alpha] <- "Downregulated"

# Add lost/retained peaks category  
resultsDF$Peak <- "No CBX2 occupancy"
resultsDF$Peak[resultsDF$row %in% CBX2_geneName$annot.symbol] <- "Increased CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Retained_uniqueGenes$annot.symbol] <- "Retained CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Lost_uniqueGenes$annot.symbol] <- "Decreased CBX2 binding"

resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & 
                 resultsDF$padj < alpha & resultsDF$Peak == "Increased CBX2 binding"] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & 
                 resultsDF$padj < alpha & resultsDF$Peak == "Retained CBX2 binding"] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & 
                 resultsDF$padj < alpha & resultsDF$Peak == "Increased CBX2 binding"] <- "Downregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & 
                 resultsDF$padj < alpha & resultsDF$Peak == "Retained CBX2 binding"] <- "Downregulated"



resultsDF %>% filter(row == "CDKN2A")
resultsDF$Condition = "CBX2 KO + EZH2i / DMSO"
resultsDF$Condition = "EZH2i / DMSO"
DESeqAnalysis(resultsDF)

# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "SYMBOL",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()

cbx2.3_1uM_Results <- resultsDF
DESeqAnalysis(cbx2.3_1uM_Results)


#------------------------------------------------------------------------------#
#                          sgCBX2.1 (1uM) VS sgLUC (DMSO)                      #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.1 (1uM EZH2i) VS sgLUC (DMSO)"

# Contrast
contrast = c("grouped","sgCBX2.11uM","sgLUCDMSO")


# Building the results table #

#resultsNames(DESeq)
results <- results(DESeq, 
                   contrast = contrast,
                   alpha = alpha, 
                   independentFiltering = F
) 
#results
print(summary(results))

resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")

# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$padj < alpha] <- "DE"
resultsDF$isDE[#resultsDF$log2FoldChange > 0 & 
                 resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[#resultsDF$log2FoldChange < -1*0 & 
                 resultsDF$padj < alpha] <- "Downregulated"

# Add lost/retained peaks category
resultsDF$Peak <- ""
resultsDF$Peak[resultsDF$row %in% CBX2_geneName$annot.symbol] <- "Increased CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Retained_uniqueGenes$annot.symbol] <- "Retained CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Lost_uniqueGenes$annot.symbol] <- "Decreased CBX2 binding"


results2.1_EZH2i = resultsDF %>% filter(isDE != "Not DE")

resultsDF$Condition = "CBX2 KO + EZH2i / DMSO"

# Get DE genes table
ezh2iDE <- merge(CBX2_geneName,
                 dplyr::filter(resultsDF, isDE == "Upregulated" | 
                                 isDE == "Downregulated"),
                 by.x = "SYMBOL",
                 by.y = "row",
                 sort = FALSE)

ezh2iDE %>% head()

treatmentCBX2.1vsDMSO_Results <- resultsDF
DESeqAnalysis(treatmentCBX2.1vsDMSO_Results)

#------------------------------------------------------------------------------#
#                          sgCBX2.3 (1uM) VS sgLUC (DMSO)                      #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.3 (1uM EZH2i) VS sgLUC (DMSO)"
# Contrast
contrast = c("grouped","sgCBX2.31uM","sgLUCDMSO")


# Building the results table #

#resultsNames(DESeq)
results <- results(DESeq, 
                   contrast = contrast,
                   alpha = alpha, 
                   independentFiltering = F
) 
#results
print(summary(results))

resultsDF <- lfcShrink(DESeq, 
                       contrast = contrast,
                       type = "ashr", 
                       res = results) %>% as.data.frame() %>% rownames_to_column("row")

# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$padj < alpha] <- "DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

resultsDF %>% nrow()

resultsKO_EZH2i = filter(resultsDF, isDE != "Not DE")
resultsKO_EZH2i = filter(resultsKO_EZH2i, row %in% results2.1_EZH2i$row)
resultsKO_EZH2i %>% nrow()

resultsDF$Condition = "CBX2 KO + EZH2i / DMSO"

resultsDF$Label <- NA
resultsDF$Label[resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$isDE != "Not DE"]

# Add lost/retained peaks category
resultsDF$Peak <- ""
#resultsDF$Peak[resultsDF$row %in% CBX2_geneName$annot.symbol] <- "CBX2 peaks at genes"
resultsDF$Peak[resultsDF$row %in% Retained_uniqueGenes$annot.symbol] <- "Increased CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Lost_uniqueGenes$annot.symbol] <- "Decreased CBX2 binding"

resultsDF %>% filter(isDE != "Not DE") %>% head()

# Get DE genes table
# ezh2iDE <- merge(CBX2_geneName, resultsDF,
#                 # dplyr::filter(resultsDF, isDE == "Upregulated" | 
#                 #                 isDE == "Downregulated"),
#                  by.x = "annot.symbol",
#                  by.y = "row",
#                  sort = FALSE)
# 
# ezh2iDE %>% head()

treatmentCBX2.3vsDMSO_Results <- resultsDF
DESeqAnalysis(treatmentCBX2.3vsDMSO_Results, 
              FC.cutoff = FC.cutoff, 
              comparison = comparison, 
              alpha = alpha)

#------------------------------------------------------------------------------#
#                          sgCBX2.3 (1uM) VS sgCBX2.3 (DMSO)                      #
#------------------------------------------------------------------------------#
# Conditions to compare
comparison = "sgCBX2.3 (1uM EZH2i) VS sgCBX2.3 (DMSO)"
# Fold-change cut-off
FC.cutoff = 0 
# q-value (FDR) cut-off
alpha = 0.05 
# Contrast
contrast = c("grouped","sgCBX2.11uM","sgCBX2.1DMSO")


# Building the results table #

resultsNames(DESeq)
results <- results(DESeq, 
                   contrast = contrast,
                   #lfcThreshold = 0,
                   alpha = 0.05, 
                   #pAdjustMethod = "bonferroni"
) 
results
print(summary(results))

resultsDF <- results(DESeq, 
                     contrast = contrast,
                     alpha = alpha, 
                     # lfcThreshold = 0,
                     # pAdjustMethod = "bonferroni",
                     tidy = TRUE)

# Add classes of differential expression
resultsDF$isDE <- "Not DE"
resultsDF$isDE[resultsDF$log2FoldChange > FC.cutoff & resultsDF$padj < alpha] <- "Upregulated"
resultsDF$isDE[resultsDF$log2FoldChange < -1*FC.cutoff & resultsDF$padj < alpha] <- "Downregulated"

resultsDF$Condition = "CBX2 KO + EZH2i / DMSO"

resultsDF$Label <- NA
resultsDF$Label[resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$isDE != "Not DE"]

# Add lost/retained peaks category
resultsDF$Peak <- ""
#resultsDF$Peak[resultsDF$row %in% CBX2_geneName$annot.symbol] <- "CBX2 peaks at genes"
resultsDF$Peak[resultsDF$row %in% Retained_uniqueGenes$annot.symbol] <- "Increased CBX2 binding"
resultsDF$Peak[resultsDF$row %in% Lost_uniqueGenes$annot.symbol] <- "Decreased CBX2 binding"

resultsDF %>% filter(isDE != "Not DE") %>% head()



################################################################################
################################################################################
# Plot numerical
df_annotatrCBX2 %>% head()
df_annotatrCBX2 %>% nrow()
annotatrCBX2 <- merge(df_annotatrCBX2, resultsDF[c(1,3,7:9)], by.x = "annot.symbol", by.y = "row", sort = F) %>% 
  filter(isDE != "Not DE", Peak != "") %>%
  GRanges()
annotatrCBX2$padj <- -log10(annotatrCBX2$padj)
annotatrCBX2

plot_numerical(
  annotated_regions = annotatrCBX2,
  x = 'log2FoldChange',
  y = 'padj',
  facet = c('annot.type', 'Peak'),
  facet_order = list(c('hg19_genes_1to5kb',"hg19_genes_promoters","hg19_genes_exons",
                       "hg19_genes_introns","hg19_genes_5UTRs","hg19_genes_3UTRs"), c("Retained CBX2 peak", "Lost CBX2 peak")),
  plot_title = 'Differentially expressed genes annotated to CBX2 peaks (sgCBX2.3/1uM vs sgLUC/DMSO)',
  x_label = 'Log2FoldChange',
  y_label = 'Log10(FDR)'
) + scale_x_continuous(limits = c(-5,5)) + scale_color_manual(values = "Log2FoldChange")


#------------------------------------------------------------------------------#
#                             GENE EXPRESSION HEATMAP                          #
#------------------------------------------------------------------------------#

#----------
# Prepare data
#----------

#########
if (!require("NMF")) {
  install.packages("NMF", dependencies = TRUE)
  }
library(NMF)

rlogData <- rlog(dataSet, blind = FALSE)
rlogData %>% assay() %>% head()

rlogCounts = assay(rlogData) %>% as.data.frame()
rlogCounts[1:6] %>% head()
rlogCounts %>% nrow()
rlogCounts = rlogCounts %>% mutate(label = rownames(rlogCounts))
rlogCounts[1:6] %>% head()

# annotatrResults <- merge(resultsDF, df_annotatrCBX2[c(20:21)], by.y = "annot.symbol", by.x = "row", sort = F) 
# annotatrResults %>% head()
# annotatrResults %>% nrow()
# 
# atPromoters <- annotatrResults %>% filter(annot.type == "hg19_genes_promoters")
# atPromoters %>% head()
# atPromoters %>% nrow()
# 
# rlogCounts$peakCategory <- "Without CBX2 occupancy"
# rlogCounts$peakCategory[rlogCounts$label %in% CBX2_geneName$annot.symbol] <- "CBX2 peaks"
# rlogCounts$peakCategory[rlogCounts$label %in% CBX2_geneName$annot.symbol &
#                        !rlogCounts$label %in% atPromoters$row] <- "CBX2 peaks NOT at promoters"
# rlogCounts$peakCategory[rlogCounts$label %in% Lost_uniqueGenes$annot.symbol &
#                        rlogCounts$label %in% atPromoters$row] <- "Lost CBX2 peaks at promoters"
# rlogCounts$peakCategory[rlogCounts$label %in% Retained_uniqueGenes$annot.symbol &
#                        rlogCounts$label %in% atPromoters$row] <- "Retained CBX2 peaks at promoters"

rlogCounts <- rlogCounts %>% 
#  group_by(peakCategory) %>% 
#  arrange(.by_group = TRUE) %>% 
  as.data.frame()
# rlogCounts <- rlogCounts[!(rlogCounts$Sorted_S01_VI_LUCA_DMSO == 0 | 
#                        rlogCounts $Sorted_S02_VI_LUCB_DMSO == 0 |
#                        rlogCounts$Sorted_S03_VI_LUCA_1uM == 0 |
#                        rlogCounts$Sorted_S04_VI_LUCB_1uM == 0 |
#                        rlogCounts$Sorted_S13_VI_CBX2.1A_DMSO == 0 |
#                        rlogCounts$Sorted_S14_VI_CBX2.1B_DMSO == 0 |
#                        rlogCounts$Sorted_S15_VI_CBX2.1A_1uM == 0 |
#                        rlogCounts$Sorted_S16_VI_CBX2.1B_1uM == 0 |
#                        rlogCounts$Sorted_S17_VI_CBX2.3A_DMSO == 0 |
#                        rlogCounts$Sorted_S18_VI_CBX2.3B_DMSO == 0 |
#                        rlogCounts$Sorted_S19_VVI_CBX2.3A_1uM == 0 |
#                        rlogCounts$Sorted_S20_VI_CBX2.3B_1uM == 0),] 
rlogCounts %>% head()
# rlogCounts$peakCategory %>% unique()
rlogCounts %>% nrow()

filter(resultsDF, isDE != "Not DE")$row %>% length()
#resultsDF = resultsDF %>% filter(row %in% results2.1$row)


# Subset to only DE expressed genes from DE analysis
rlogDE = rlogCounts %>% 
  filter(label %in% filter(resultsDF, isDE != "Not DE")$row) %>%   
  as.data.frame()

# Subset to only DE expressed genes from DE analysis
rlogCBX2 = rlogCounts %>% 
  filter(#peakCategory != "Without CBX2 occupancy",
    label %in% filter(resultsKO_EZH2i, isDE != "Not DE")$row
  ) %>%   
#  group_by(peakCategory) %>% 
#  arrange(.by_group = TRUE) %>% 
  as.data.frame()

rlogCBX2 %>% head()
rlogCBX2 %>% nrow()

#par(las=0)
aheatmap(as.matrix(rlogCBX2[1:12]),
         #color = c("black","darkgoldenrod1","cornsilk"),
         color = c("black","white","darkgoldenrod1"), 
         ##border_color = "black",
         scale="row", 
         #annColors = c("Set2"), 
         distfun = "pearson", 
         hclustfun = "average", 
         #main = "Heatmap of expression levels (rLog)", 
         main = paste0("Differentially expressed genes (CBX2 KO + EZH2i vs DMSO)\n N = ",nrow(rlogCBX2)), 
         #main = "Scaled rLog transformation", # to sum-up to one 
         #annRow=rlogCBX2$peakCategory, 
         #annCol = colnames(newFPKM),
         #Rowv = NA,  
         #Colv = NA,
         cexRow = .6,
         cexCol = .8,
         info = TRUE, 
         labRow = NA, gp = gpar(rot= 45),
         #labRow = rlogCBX2$label, 
         labCol = c("sgLUC + DMSO", "sgLUC + DMSO", "sgLUC + 1uM EZH2i", "sgLUC + 1uM EZH2i",
                    "sgCBX2.1 + DMSO", "sgCBX2.1 + DMSO", "sgCBX2.1 + 1uM EZH2i", "sgCBX2.1 + 1uM EZH2i",
                    "sgCBX2.3 + DMSO", "sgCBX2.3 + DMSO", "sgCBX2.3 + 1uM EZH2i", "sgCBX2.3 + 1uM EZH2i") )


## Try pheatmap function
library(pheatmap)
library(RColorBrewer)
pdf(file = "./heatmap_treat_names.pdf", height = 20, width = 8)

heatmap = pheatmap(as.matrix(rlogCounts[1:4]),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
         scale="row", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         cutree_rows = 2,
         show_rownames = F, 
         fontsize_row = 5,
         main = paste0("Differentially expressed genes (EZH2i vs DMSO) N = ",nrow(rlogCounts[1:4])),
         angle_col = "0", 
         #run_draw - T,
         labels_col = c("sgLUC + DMSO", "sgLUC + DMSO", "sgLUC + 1uM EZH2i", "sgLUC + 1uM EZH2i",
                        "sgCBX2.1 + DMSO", "sgCBX2.1 + DMSO", "sgCBX2.1 + 1uM EZH2i", "sgCBX2.1 + 1uM EZH2i",
                        "sgCBX2.3 + DMSO", "sgCBX2.3 + DMSO", "sgCBX2.3 + 1uM EZH2i", "sgCBX2.3 + 1uM EZH2i"))
dev.off()

geneOrder = data.frame("order" = heatmap$tree_row$order)
geneOrder$name = rlogCounts$label[geneOrder$order]

hclust = hclust(d = dist(as.matrix(rlogCounts[1:4]), method = "euclidean"), method = "ward.D")

geneOrder = geneOrder[1:2]
geneOrder %>% head()
geneOrder$cluster = NA
geneOrder$cluster[1:22] = "A"
geneOrder$cluster[23:329] = "B"

library("AnnotationDbi")
library("org.Hs.eg.db")
geneOrder$ensid = mapIds(org.Hs.eg.db,
                  keys=geneOrder$name, 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")

if (!require("EnsDb.Hsapiens.v75", quietly = TRUE)){
  BiocManager::install("EnsDb.Hsapiens.v75")
}
  


# Get hg19 gene coordinates and filter by those that are differentially expressed
library("EnsDb.Hsapiens.v75")

hg19genes = genes(EnsDb.Hsapiens.v75) %>% data.frame() %>% dplyr::filter(gene_id %in% geneOrder$ensid) %>% GRanges()
hg19genes

# Get TSS coordinates
geneTSS = promoters(hg19genes, upstream = 0, downstream = 1, use.names = T)

# Arrange in custom order (the order in the heatmap output)

geneTSS = geneTSS %>% data.frame() %>% arrange(factor(gene_id, levels = geneOrder$ensid)) %>% GRanges()

seqlevelsStyle(geneTSS) <- "UCSC"
geneTSS

geneTSS = setNames(geneTSS, geneTSS$gene_id)

library(EnrichedHeatmap)

# Import H3K27me3 normalized CUT&RUN signal over genomic coordinates
H3K27me3 = import.bedGraph("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/bedgraph/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.TMM.binned.bedgraph")
H3K27me3

H3K27me3_1uM = import.bedGraph("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/bedgraph/RP_054_H3K27me3_500K_1uMcpd_Rpcells_S54.TMM.binned.bedgraph")


mat1_trim = normalizeToMatrix(signal = H3K27me3, 
                              target = geneTSS, 
                              value_column = "score", 
                              extend = 5000, 
                              mean_mode = "w0", 
                              w = 100,
                              keep = c(0, 0.99),
                              smooth = TRUE 
                              )

mat2_trim = normalizeToMatrix(signal = H3K27me3_1uM, 
                              target = geneTSS, 
                              value_column = "score", 
                              extend = 5000, 
                              mean_mode = "w0", 
                              w = 100,
                              keep = c(0, 0.99),
                              smooth = TRUE)

pdf(file = "./mat1_trim.pdf", height = 6, width = 3)
print(EnrichedHeatmap(mat1_trim, 
                      col = c("white", "red"), 
                      name = "H3K27me3", 
                      row_split = geneOrder$cluster[1:269])
      )
dev.off()

#------
# GSEA
#------

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in data from deseq2
df = read.csv("./resultsRanked.csv", header=TRUE)
df %>% head()

df = read.csv("./resultsLogFC.csv", header=TRUE)

# we want the some ranked metric, in this case, it is a combination of the adjusted p-value and log2FC 
original_gene_list <- df$rank

original_gene_list <- df$log2FoldChange
# name the vector
names(original_gene_list) <- df$row
original_gene_list %>% head()

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# GENE SET ENRICHMENT
gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPermSimple = 100000, 
             eps = 0,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

# dotplot plot
require(DOSE)
dotplot(gse, showCategory=15, split=".sign") + facet_grid(.~.sign)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gse@result %>% head()



gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
