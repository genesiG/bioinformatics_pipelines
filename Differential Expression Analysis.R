# Set the working directory
directory <- "C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2"
setwd(directory)

if(!require(DESeq2)){
  BiocManager::install("DESeq2")
}

if(!require(AnnotationDbi)){
  BiocManager::install("AnnotationDbi")
}

if(!require(org.Hs.eg.db)){
  BiocManager::install("org.Hs.eg.db")
}

library(DESeq2)
library(RColorBrewer)
library(dplyr)
library(tidyverse)
library(ggplot2)

# Define cutoffs
FDR.cutOff=0.05			# 5% FDR cutoff
LFC.cutOff=0.5			# DE genes must have ? logFC in either direction

#------------------------------------------------------------------------------#
#                   CREATING DESeqDataSet FROM COUNTS MATRIX                   #
#------------------------------------------------------------------------------#


#----------
# We can extract the counts matrix from a DESeqDataSet                         
# For example, if a DESeqDataSet is sent by a collaborator, we can use:        
#                                                                              
# dds <- load("ddsFromCollaborator.RData") # load DESeqDataSet                 
# class(dds)                                                                     
# countData <- assay(dds) # extract counts matrix from DESeqDataSet            
#                                                                                
# colnames(countData) # check which columns correspond to CBX2 data              
#                                                                              
# cbx2Counts <- countData[,13:20] # extract columns for all CBX2 data          
# cbx2Counts %>% head()                                                        
#                                                                                
# * Alternatively:                                                             
#     cbx2Counts1 <- countData[,13:16] # extract columns for CBX2.1 data       
#     cbx2Counts3 <- countData[,17:20] # extract columns for CBX2.3 data       
#                                                                              
# colnames(cbx2Counts) # make sure which is the name of each column            
#                                                                              
#     colnames(cbx2Counts1) # make sure which is the name of each column       
#     colnames(cbx2Counts3) # make sure which is the name of each column       
#----------

# Alternatively we can load the count matrix generated from featureCounts
sgLUCCounts <- read.delim("lucCounts.txt", header = T, row.names = 1)

sgLUCCounts %>% head() # check the loaded count matrix
colnames(sgLUCCounts) # make sure columns correctly correspond to samples' names

#----------
# Create the meta-data data frame to help DESeq understand the comparative     
# structure of your data. Make sure the values in the condition argument are in 
# the same order as your samples appear in the colnames()                      
#----------

sgLUCMetaData <- data.frame(row.names = colnames(sgLUCCounts), # samples names
                           treatment = factor(c("DMSO", # grouping
                                                "DMSO",
                                                "1uM",
                                                "1uM")),
                           replicate = factor(c("A","B","A","B")))

sgLUCMetaData 

#----------
# This basically assigns the first two data columns (LUC.A_DMSO and 
# LUC.B_DMSO) to the same experimental group ("DMSO"); ditto with the other 
# data columns.
#----------

# Check the factor levels of the condition. The 1st level should be the control.
levels(sgLUCMetaData$treatment)

#----------
# DESeq will calculate the fold change of the LAST factor level over the first.
#----------

# If your control is NOT the 1st level, make it so before runnning the pipeline
sgLUCMetaData$treatment <- relevel(sgLUCMetaData$treatment, "DMSO")
levels(sgLUCMetaData$treatment)

#----------
# We now have all the ingredients to prepare our data object in a form that is 
# suitable for analysis, namely:                                               
#    countData: a table of the fragment counts and sample names as columns     
#    colData: a table of sample grouping (sample names as rows)                
#----------

#----------
# If the research aim is to determine for which genes the effect of treatment 
# is different across groups, then interaction terms can be included and tested 
# using a design such as ~ group + treatment + group:treatment
#----------

sgLUCDataSet <- DESeqDataSetFromMatrix(countData = sgLUCCounts, # count matrix
                                      colData = sgLUCMetaData, # samples&grouping
                                      design = ~ treatment)

#----------
# DESeqDataSetFromMatrix(countData, colData, design, tidy = F)                 
#    countData: the matrix with your summarized read counts                    
#    colData: a data.frame whose rows correspond to the columns in countData   
#    design: formula or matrix expressing how the counts for each gene relate  
#            # to the variables in colData                                     
#    tidy: whether the first column of countData is the rownames               
#----------

sgLUCDataSet

#------------------------------------------------------------------------------#
#                   EXPLORATORY ANALYSIS AND VISUALIZATION                     #
#------------------------------------------------------------------------------#

#---------------------------#
# Pre-filtering the dataset #
#---------------------------#

# Here we apply the most minimal filtering rule: removing rows of the 
# DESeqDataSet that have no counts, or only a single count across all samples.

nrow(sgLUCDataSet)

keep <- rowSums(counts(sgLUCDataSet)) > 1
sgLUCDataSet <- sgLUCDataSet[keep,]
nrow(sgLUCDataSet)

#------------------------------------------------------#
# The variance stabilizing transformation and the rlog #
#------------------------------------------------------#

?vst

?rlog

#----------
# Many common statistical methods for exploratory analysis of multidimensional 
# data, for example clustering and principal components analysis (PCA), work   
# best for data that generally has the same range of variance at different     
# ranges of the mean values.                                                   
#                                                                              
# When the expected amount of variance is approximately the same across        
# different mean values, the data is said to be HOMOSKEDASTIC. For RNA-seq     
# counts, however, the expected variance grows with the mean.                  
#                                                                              
# For example, if one performs PCA directly on a matrix of counts or normalized
# counts (e.g. correcting for differences in sequencing depth), the resulting    
# plot typically depends mostly on the genes with highest counts because they  
# show the largest absolute differences between samples.                        
#                                                                              
# A simple and often used strategy to avoid this is to take the logarithm of   
# the normalized count values plus a pseudocount of 1; however, depending on   
# the choice of pseudocount, now the genes with the very lowest counts will    
# contribute a great deal of noise to the resulting plot, because taking the   
# logarithm of small counts actually inflates their variance.                  
#                                                                              
# As a solution, DESeq2 offers two transformations for count data that         
# stabilize the variance across the mean: the variance stabilizing             
# transformation (VST) for negative binomial data with a dispersion-mean trend 
# (Anders and Huber 2010), implemented in the vst function, and the            
# regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014). 
#                                                                              
# For genes with high counts, both the VST and the rlog will give similar      
# result to the ordinary log2 transformation of normalized counts. For genes   
# with lower counts, however, the values are shrunken towards a middle value.  
# The VST or rlog-transformed data then become approximately homoskedastic     
# (more flat trend in the meanSdPlot), and can be used directly for computing  
# distances between samples, making PCA plots, or as input to downstream       
# methods which perform best with homoskedastic data.                          
#                                                                                 
# Which transformation to choose? The VST is much faster to compute and is less
# sensitive to high count outliers than the rlog. The rlog tends to work well  
# on small datasets (n < 30), potentially outperforming the VST when there is a
# wide range of sequencing depth across samples (an order of magnitude         
# difference). We therefore recommend the VST for medium-to-large datasets     
# (n > 30). You can perform both transformations and compare the meanSdPlot or 
# PCA plots generated, as described below.                                     
#                                                                              
# Note that the two transformations offered by DESeq2 are provided for         
# applications other than differential testing. For differential testing we    
# recommend the DESeq function applied to raw counts, as described later in    
# this workflow, which also takes into account the dependence of the variance  
# of counts on the mean value during the dispersion estimation step.           
#                                                                              
# Both vst and rlog return a DESeqTransform object which is based on the       
# SummarizedExperiment class. The transformed values are no longer counts, and 
# are stored in the assay slot. The colData that was attached to dds is still  
# accessible through the colData() function                                    
#----------

rlogDataSet <- rlog(sgLUCDataSet, blind = TRUE)

#----------
# We would specify blind = FALSE, if we expected that differences between the 
# variables in the design (say, cell lines and treatment) would not contribute 
# to the expected variance-mean trend of the experiment. The experimental design
# then is not used directly in the transformation, only in estimating the global
# amount of variability in the counts. 
#
# For a fully unsupervised transformation, we can set blind = TRUE (the default)                                 
#----------

#------------------#
# Sample distances #
#------------------#

#----------
# A useful first step in an RNA-seq analysis is often to assess overall        
# similarity between samples: Which samples are similar to each other, which    
# are different? Does this fit to the expectation from the experiment’s design?
#                                                                              
# Here we use the calculating sample distances is to use the Poisson Distance  
# (Witten 2011), implemented in the PoiClaClu package. This measure of         
# dissimilarity between counts also takes the inherent variance structure of   
# counts into consideration when calculating the distances between samples.    
#                                                                              
# If we wanted to calculate the Euclidean distances between the samples with   
# the R function dist(), we would need to use the vst- or rlog-transformed data
# instead. But the PoissonDistance function takes the ORIGINAL count matrix    
# (NOT normalized) instead.                                                    
#----------

if(!require(PoiClaClu)){
  install.packages("PoiClaClu")
}
library("PoiClaClu")

counts(sgLUCDataSet) %>% head() # see that in the count matrix the sample names
# are columns

sampleDistances <- PoissonDistance(x = t(counts(sgLUCDataSet)), # transpose matrix
                                   type = "deseq") 

#----------
# We need to transpose the matrix of values using t(), because the both dist   
# and the PoissonDistance functions expect the different samples to be rows of 
# its argument, and different dimensions (here, genes) to be columns.          
#----------

sampleDistances  
sampleDistances$dd # the n x n dissimilarity matrix

# We then visualize the distances using the the pheatmap package

sampleDistMatrix <- as.matrix(sampleDistances$dd)

rownames(sampleDistMatrix) <- paste(sgLUCDataSet$treatment, 
                                    sgLUCDataSet$experiment, 
                                    sep = " - " )
colnames(sampleDistMatrix) <- NULL

#----------
# In order to plot the sample distance matrix with the rows/columns arranged by
# the distances in our distance matrix, we manually provide sampleDists to the 
# clustering_distance argument of the pheatmap function. Otherwise the pheatmap
# function would assume that the matrix contains the data values themselves,   
# and would calculate distances between the rows/columns of the distance matrix
# which is not desired. We also manually specify a blue color palette using the
# colorRampPalette function from the RColorBrewer package.                     
#----------

library("pheatmap")
library("RColorBrewer")

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDistances$dd,
         clustering_distance_cols = sampleDistances$dd,
         col = colors)

#----------#
# PCA Plot #
#----------#

#----------
# Another way to visualize sample-to-sample distances is a principal components 
# analysis (PCA). In this ordination method, the data points (here, the        
# samples) are projected onto the 2D plane such that they spread out in the two
# directions that explain most of the differences.                             
#                                                                              
# The x-axis is the direction that separates the data points the most. The     
# values of the samples in this direction are written PC1.                     
#                                                                              
# The y-axis is an orthogonal direction to PC1 that separates the data the     
# second most. The values of the samples in this direction are written PC2.    
#                                                                              
# The percent of the total variance that is contained in the direction is      
# printed in the axis label. Note that these percentages do not add to 100%,   
# because there are more dimensions that contain the remaining variance        
# (although each of these remaining dimensions will explain less than the two  
# that we see).                                                                
#----------

plotPCA(rlogDataSet, intgroup = c("treatment", "experiment"))+
  aes(color = treatment, # color according to treatment condition
      shape = experiment) + # shape points according to experiment
  ggtitle(label = "Principal Component Analysis", 
          subtitle = "Top 500 genes with highest variance") +
  theme_minimal()

#----------
# Here, we have used the function plotPCA that comes with DESeq2. The two terms
# specified by intgroup are the interesting groups for labeling the samples;   
# they tell the function to use them to choose colors.                         
#----------

#--------------------------------#
# PCA Plot using Generalized PCA #
#--------------------------------#

#----------
# Another technique for performing dimension reduction on data that is NOT     
# Normally distributed (e.g. over-dispersed count data) is generalized         
# principal component analysis, or GLM-PCA, (Townes et al. 2019) as implemented
# in the CRAN package glmpca. This package takes as input the count matrix,    
# with features (genes) as rows, and observations (samples) as columns. It also
# requires the number of latent dimensions to fit (here, we specify 2).        
#                                                                              
# As stated by Townes et al. (2019):                                           
#                                                                              
#  "…we propose the use of GLM-PCA, a generalization of PCA to exponential     
#  family likelihoods. GLM-PCA operates on raw counts, avoiding the pitfalls   
#  of normalization. We also demonstrate that applying PCA to deviance or      
#  Pearson residuals provides a useful and fast approximation to GLM-PCA."     
#                                                                              
#----------

if(!require(glmpca)){
  install.packages("glmpca")
}
library("glmpca")

?glmpca # check what the arguments for the function mean

sgLUCGPCA <- glmpca(counts(sgLUCDataSet), # use the count matrix in DESeqDataSet 
                   L=2) # specify two dimensions for reduction

sgLUCGPCA # summary of the model
sgLUCGPCA$factors # matrix whose rows match the columns of the count matrix and
# whose columns are the different latent dimensions. Analogous
# to the principal components in PCA model
sgLUCGPCA$loadings # matrix with rows matching the rows (features/genes) of the
# count matrix. Columns are different latent dimensions

gpca.data <- sgLUCGPCA$factors # extract principal components from the GPCA model
gpca.data$treatment <- sgLUCDataSet$treatment
gpca.data$experiment <- sgLUCDataSet$experiment

ggplot(gpca.data, aes(x = dim1, 
                      y = dim2, 
                      color = treatment, 
                      shape = experiment)) +
  geom_point(size =3) + 
  coord_fixed() + 
  ggtitle("Generalized Principal Component Analysis") +
  theme_minimal()

#------------------------------------------------------------------------------#
#                       DIFFERENTIAL EXPRESSION ANALYSIS                       #
#------------------------------------------------------------------------------#

#----------------------------------------------#
# Running the differential expression pipeline #
#----------------------------------------------#

?DESeq

sgLUCDESeq <- DESeq(sgLUCDataSet)
sgLUCDESeq

#----------------------------#
# Building the results table #
#----------------------------#

sgLUCResults <- results(sgLUCDESeq)

#----------
# Calling results without any arguments will extract the estimated log2 fold   
# changes and p values for the last variable in the design formula. If there   
# are more than 2 levels for this variable, results will extract the results   
# table for a comparison of the last level over the first level. The comparison
# is printed at the top of the output.                                         
#----------

sgLUCResults # the comparison should be TREATED vs CONTROL

# As res is a DataFrame object, it carries metadata with information on the 
# meaning of the columns:

mcols(sgLUCResults, use.names = TRUE)

#----------
# The first column, baseMean, is just the average of the normalized count      
# values, divided by the size factors, taken over all samples in the           
# DESeqDataSet. The remaining four columns refer to a specific contrast, namely
# the comparison of the treatment level over the control level for the factor  
# variable treatment. We will find out below how to obtain other contrasts.    
#                                                                              
# The column log2FoldChange is the effect size estimate. It tells us how much  
# the gene’s expression seems to have changed due to treatment in comparison to
# untreated samples. This value is reported on a logarithmic scale to base 2:  
# for example, a log2 fold change of 1.5 means that the gene’s expression is   
# increased by a multiplicative factor of 21.5≈2.82.                           
#                                                                              
# Of course, this estimate has an uncertainty associated with it, which is     
# available in the column lfcSE, the standard error estimate for the log2 fold 
# change estimate.                                                             
#                                                                              
# We can also express the uncertainty of a particular effect size estimate as  
# the result of a statistical test. The purpose of a test for differential     
# expression is to test whether the data provides sufficient evidence to       
# conclude that this value is really different from zero. DESeq2 performs for  
# each gene a hypothesis test to see whether evidence is sufficient to decide  
# against the null hypothesis (that there is zero effect of the treatment on   
# the gene and that the observed difference between treatment and control was  
# merely caused by experimental variability, i.e., the type of variability that
# you can expect between different samples in the same treatment group).       
#                                                                                
# As usual in statistics, the result of this test is reported as a p value, and
# it is found in the column pvalue. Remember that a p value indicates the      
# probability that a fold change as strong as the observed one, or even        
# stronger, would be seen under the situation described by the null hypothesis.                                                                  #
#----------

summary(sgLUCResults)

#----------
# Here we see that many genes are differentially expressed in response to 
# treatment when using a log2foldChange threshold of 0 and a false discovery 
# rate (FDR) level of 10%
#----------

table(sgLUCResults$padj < 0.1)

#----------
# There are two ways to be more strict about which set of genes are considered 
# significant:
#   *lower the FDR threshold (the threshold on padj in the results table)
#   *raise the fold-change threshold using the lfcThreshold argument of results
#----------

# Reduce FDR (adjusted p-value) to 5%
sgLUCResults <- results(sgLUCDESeq,
                        alpha = 0.01) 

# How many genes are differentially expressed using this threshold?
table(sgLUCResults$padj < 0.01) 

# Raise the LFC threshold to 1 (only genes expressed 2x more or 50% less)
sgLUCResults <- results(sgLUCDESeq,
                        lfcThreshold = 0.5) 

# How many genes are differentially expressed using this LFC threshold?
table(sgLUCResults$padj < 0.1) 

# How many genes are differentially expressed if we combine both thresholds?
table(sgLUCResults$padj < 0.05) 

# Only genes with a 1.5x change in expression
# sgLUCResultsLFC.58 <- results(sgLUCDESeq, 
#                              lfcThreshold=0.58) 
# table(sgLUCResultsLFC.58$padj < 0.1)



#------------------------------------------#
# Plotting differential expression results #
#------------------------------------------#

sgLUCResults <- results(sgLUCDESeq,
                        alpha = 0.05, 
                        lfcThreshold = 0.6) 
table(sgLUCResults$padj < 0.05)

plotMA(sgLUCResults) # Create MA plot of differential expression

# Write results to comma-delimited file or tab-delimited text file
write.csv(as.data.frame(sgLUCResults), 
          file="sgLUCResults.csv")
write.table(as.data.frame(sgLUCResults),
            sep = "\t",
            file="sgLUCResults.txt")

#------------------------------------------------#
# Getting list of differentially expressed genes #
#------------------------------------------------#

# Get back the data used for plotting as a data frame with DE data
sgLUCDE <- plotMA(sgLUCResults, returnData = TRUE)
sgLUCDE %>% head() # note that there is a column showing which gene is DE

# Get DESeq results as data.frame with rownames as the first column
sgLUCResultsDF <- results(sgLUCDESeq,
                          alpha = 0.01,
                          #lfcThreshold = 0.5,
                          tidy = TRUE)
sgLUCResultsDF %>% head()

# Drop columns 4, 5, and 6
sgLUCResultsDF <- sgLUCResultsDF[,-c(4:6)]
sgLUCResultsDF %>% head()

# Check if both data frames have the same number of rows
nrow(sgLUCResultsDF) == nrow(sgLUCDE) # should return TRUE

# Combine by columns (cbind)
sgLUCResultsDF <- cbind(sgLUCResultsDF, sgLUCDE[,2:3])
sgLUCResultsDF %>% head()

# Check if columns for log2 fold-change match (they should)
(sgLUCResultsDF$log2FoldChange == sgLUCResultsDF$lfc)[FALSE] %>% 
  length() # should return zero
(sgLUCResultsDF$log2FoldChange == sgLUCResultsDF$lfc)[TRUE] %>% 
  length() # should be equal to number of rows of the data frame
(sgLUCResultsDF$log2FoldChange == sgLUCResultsDF$lfc)[TRUE] %>% 
  length() == nrow(sgLUCResultsDF) # should return TRUE

# Drop redundant lfc column
sgLUCResultsDF <- sgLUCResultsDF[,-5]
sgLUCResultsDF %>% head()
sgLUCResultsDF$row %>% unique() %>% length() # number of genes

# Make classes of differential expression
sgLUCResultsDF$isDE <- "Not DE"
# if log2Foldchange > 0.5 and FDR < 0.05, set as "UP" 
sgLUCResultsDF$isDE[sgLUCResultsDF$log2FoldChange > 0.5 & sgLUCResultsDF$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.5 and FDR < 0.05, set as "DOWN"
sgLUCResultsDF$isDE[sgLUCResultsDF$log2FoldChange < -0.5 & sgLUCResultsDF$padj < 0.05] <- "Downregulated"

sgLUCResultsDF %>% head()

#-------------------------------------------------------------------------------
# Annotate peaks to differential expression data
#-------------------------------------------------------------------------------

# Import peak annotation from Homer
CBX2_DMSO_peakAnno <- read.table("./Homer/CBX2_DMSO.txt", 
                                 sep = "\t", 
                                 header = T)


CBX2_lost_peakAnno <- read.table("./Homer/CBX2_lost_DE.txt", 
                                 sep = "\t", 
                                 header = T)
CBX2_lost_peakAnno %>% head()
CBX2_lost_peakAnno %>% nrow() # number of lost peaks

CBX2_retained_peakAnno <- read.table("./Homer/CBX2_retained.txt",
                                     sep = "\t", 
                                     header = T)
CBX2_retained_peakAnno %>% nrow() # number of retained peaks

CBX2_lost_peakAnno %>% nrow() + CBX2_retained_peakAnno %>% nrow() # number of peaks
CBX2_DMSO_peakAnno %>% nrow() # total number of peaks

# Get only the Gene Names column
CBX2_DMSO_peakAnno %>% colnames() 
CBX2_lost_peakAnno %>% colnames()
CBX2_retained_peakAnno %>% colnames()

CBX2_geneName <- CBX2_DMSO_peakAnno[16] %>% unique() # genes annotated to peaks
CBX2_geneName %>% head()
CBX2_geneName %>% nrow() # number of genes annotated to peaks

CBX2_lost_geneName <- CBX2_lost_peakAnno[16] 
CBX2_lost_geneName %>% head()
CBX2_lost_geneName %>% nrow()
CBX2_lost_uniqueGenes <- filter(CBX2_lost_geneName, !Gene.Name %in% CBX2_retained_geneName$Gene.Name) 
CBX2_lost_uniqueGenes %>% nrow()

CBX2_retained_geneName <- CBX2_retained_peakAnno[16] 
CBX2_retained_geneName %>% head()
CBX2_retained_geneName %>% nrow()
CBX2_retained_uniqueGenes <- filter(CBX2_retained_geneName, !Gene.Name %in% CBX2_lost_geneName$Gene.Name) 
CBX2_retained_uniqueGenes %>% nrow()

# Merge annotated peaks and gene expression data by gene name
sgLUCResultsDF %>% nrow() 
sgLUCResultsDF$row %>% unique() %>% length()# number of genes from DE analysis

sgLUC_DMSOPeakExpression <- merge(CBX2_geneName,
                                 sgLUCResultsDF,
                                 by.x = "Gene.Name",
                                 by.y = "row",
                                 sort = FALSE)
sgLUC_DMSOPeakExpression %>% nrow()
sgLUC_DMSOPeakExpression %>% head()

#-------------------------------------------------------------------------------
# Volcano Plot
#-------------------------------------------------------------------------------

if (!require("ggrepel")) {
  install.packages("ggrepel")
  library(ggrepel)
}
library(ggplot2)

# Data frame should have at least a column for log2fold change and one for FDR

# Add column for labels of differentially expressed genes
sgLUC_DMSOPeakExpression$Label <- NA

sgLUC_DMSOPeakExpression$Label[sgLUC_DMSOPeakExpression$isDE != "Not DE"] <- sgLUC_DMSOPeakExpression$Gene.Name[sgLUC_DMSOPeakExpression$isDE != "Not DE"]

sgLUC_DMSOPeakExpression %>% tail()

# Add lost/retained peaks category
sgLUC_DMSOPeakExpression$Peak <- ""
sgLUC_DMSOPeakExpression %>% head()
sgLUC_DMSOPeakExpression$Peak[sgLUC_DMSOPeakExpression$Gene.Name %in% CBX2_retained_uniqueGenes$Gene.Name] <- "Retained CBX2 peak"
sgLUC_DMSOPeakExpression$Peak[sgLUC_DMSOPeakExpression$Gene.Name %in% CBX2_lost_uniqueGenes$Gene.Name] <- "Lost CBX2 peak"

# Remove genes not annotated to either lost or retained peaks
sgLUC_DMSOPeakExpression <- filter(sgLUC_DMSOPeakExpression, Peak != "")
sgLUC_DMSOPeakExpression%>% tail()

# Make volcano plot

volcanoPlot <- ggplot(data=sgLUC_DMSOPeakExpression, 
                          aes(x=log2FoldChange, 
                              y=-log10(padj),
                              color = isDE, 
                              shape = Peak, 
                              label =  Label)) + 
  geom_point() + 
    theme_bw() +
  scale_shape_manual(values = c(1,16)) +
  ylab("log10(adjusted p-value)") +
  xlab("log2(Fold-Change)") +
  ggtitle("1 uM EZH2i vs DMSO") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

volcanoPlot + geom_text_repel() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  # Change point color 
  # By default, it is assigned to the categories in an alphabetical order):
  scale_color_manual(values=c("Downregulated" = "black","Not DE" = "gray", "Upregulated" = "orange")) 

#------------------------------------------------------------------------------#
#                         sgCBX2 (DMSO) VS sgLUC (DMSO)                        #
#------------------------------------------------------------------------------#

cbx2LUCCounts <- read.delim("cbx2vsLUCCounts.txt", header = T, row.names = 1)

cbx2LUCCounts %>% head() # check the loaded count matrix
cbx2LUCCounts %>% colnames()

DMSOCounts <- cbx2LUCCounts[,-c(3,4,5:8,11:12)]
colnames(DMSOCounts) # make sure columns correctly correspond to samples' names

# Create the meta-data data frame with the comparative structure of the data.                     

DMSOMetaData <- data.frame(row.names = colnames(DMSOCounts), # samples names
                            treatment = factor(c("DMSO", # grouping
                                                 "DMSO")),
                            sgRNA = factor(c("sgLUC", 
                                             "sgLUC", 
                                             "sgCBX2.3", 
                                             "sgCBX2.3")))

DMSOMetaData 

# Check the factor levels of the condition. The 1st level should be the control.
levels(DMSOMetaData$sgRNA)

DMSOMetaData$sgRNA <- relevel(DMSOMetaData$sgRNA, "sgLUC")
levels(DMSOMetaData$sgRNA)

DMSODataSet <- DESeqDataSetFromMatrix(countData = DMSOCounts, # count matrix
                                       colData = DMSOMetaData, # samples&grouping
                                       design = ~ sgRNA)

DMSODataSet

# Pre-filtering
nrow(DMSODataSet)

keep <- rowSums(counts(DMSODataSet)) > 1
DMSODataSet <- DMSODataSet[keep,]
nrow(DMSODataSet)

DMSODESeq <- DESeq(DMSODataSet)
DMSODESeq

# Building the results table #
DMSOResults <- results(DMSODESeq,
                        alpha = 0.05, 
                        lfcThreshold = 0.6) 
DMSOResults # the comparison should be TREATED vs CONTROL

summary(DMSOResults)
table(DMSOResults$padj < 0.05)

# Plotting differential expression results 

plotMA(DMSOResults) # Create MA plot of differential expression

# Write results to comma-delimited file or tab-delimited text file
# write.csv(as.data.frame(DMSOResults), 
#           file="DMSOResults.csv")
# write.table(as.data.frame(DMSOResults),
#             sep = "\t",
#             file="DMSOResults.txt")


#----------
# Getting list of differentially expressed genes 
#----------

# Get back the data used for plotting as a data frame with DE data
DMSODE <- plotMA(DMSOResults, returnData = TRUE)
DMSODE %>% head() # note that there is a column showing which gene is DE

# Get DESeq results as data.frame with rownames as the first column
DMSOResultsDF <- results(DMSODESeq, tidy = TRUE)
DMSOResultsDF %>% head()

# Drop columns 4, 5, and 6
DMSOResultsDF <- DMSOResultsDF[,-c(4:6)]
DMSOResultsDF %>% head()

# Check if both data frames have the same number of rows
nrow(DMSOResultsDF) == nrow(DMSODE) # should return TRUE

# Combine by columns (cbind)
DMSOResultsDF <- cbind(DMSOResultsDF, DMSODE[,2:3])
DMSOResultsDF %>% head()

# Drop redundant lfc column
DMSOResultsDF <- DMSOResultsDF[,-5]
DMSOResultsDF %>% head()
DMSOResultsDF$row %>% unique() %>% length() # number of genes

# Make classes of differential expression
DMSOResultsDF$isDE <- "Not DE"
# if log2Foldchange > 0.5 and FDR < 0.05, set as "UP" 
DMSOResultsDF$isDE[DMSOResultsDF$log2FoldChange > 0.6 & DMSOResultsDF$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.5 and FDR < 0.05, set as "DOWN"
DMSOResultsDF$isDE[DMSOResultsDF$log2FoldChange < -0.6 & DMSOResultsDF$padj < 0.05] <- "Downregulated"

DMSOResultsDF %>% head()

#----------
# Annotate peaks to differential expression data
#----------

# Merge annotated peaks and gene expression data by gene name
DMSOResultsDF %>% nrow() 
DMSOResultsDF$row %>% unique() %>% length()# number of genes from DE analysis

DMSOPeakExpression <- merge(CBX2_geneName,
                            DMSOResultsDF,
                            by.x = "Gene.Name",
                            by.y = "row",
                            sort = FALSE)
DMSOPeakExpression %>% nrow()
DMSOPeakExpression %>% head()

#----------
# Volcano Plot
#----------

# Add column for labels of differentially expressed genes
DMSOPeakExpression$Label <- NA

DMSOPeakExpression$Label[DMSOPeakExpression$isDE != "Not DE"] <- DMSOPeakExpression$Gene.Name[DMSOPeakExpression$isDE != "Not DE"]

DMSOPeakExpression %>% head()

# Add lost/retained peaks category
DMSOPeakExpression$Peak <- ""
DMSOPeakExpression %>% head()
DMSOPeakExpression$Peak[DMSOPeakExpression$Gene.Name %in% CBX2_retained_uniqueGenes$Gene.Name] <- "Retained CBX2 peak"
DMSOPeakExpression$Peak[DMSOPeakExpression$Gene.Name %in% CBX2_lost_uniqueGenes$Gene.Name] <- "Lost CBX2 peak"

# Remove genes not annotated to either lost or retained peaks
DMSOPeakExpression <- filter(DMSOPeakExpression, Peak != "")
DMSOPeakExpression%>% tail()

# Make volcano plot

volcanoPlotDMSO <- ggplot(data=DMSOPeakExpression, 
                          aes(x=log2FoldChange, 
                              y=-log10(padj),
                              color = isDE, 
                              shape = Peak, 
                              label =  Label)) + 
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values = c(1,16)) +
  ylab("log10(adjusted p-value)") +
  xlab("log2(Fold-Change)") +
  ggtitle("sgCBX2.3 (DMSO) vs sgLUC (DMSO)") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

volcanoPlotDMSO + geom_text_repel() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  # Change point color 
  # By default, it is assigned to the categories in an alphabetical order):
  scale_color_manual(values=c("black", "gray", "orange")) 

#------------------------------------------------------------------------------#
#                  sgCBX2 (1 uM EHZ2i) VS sgLUC (1 uM EZH2i)                   #
#------------------------------------------------------------------------------#

cbx2LUCCounts %>% head() # check the loaded count matrix
cbx2.1Counts <- cbx2LUCCounts[,-c(1,2,5:10)]
colnames(cbx2.1Counts) # make sure columns correctly correspond to samples' names

# Create the meta-data data frame with the comparative structure of the data.                     

cbx2.1MetaData <- data.frame(row.names = colnames(cbx2.1Counts), # samples names
                              treatment = factor(c("1uM")),
                              sgRNA = factor(c("sgLUC", 
                                               "sgLUC", 
                                               "sgCBX2.3",
                                               "sgCBX2.3")))

cbx2.1MetaData

# Check the factor levels of the condition. The 1st level should be the control.
levels(cbx2.1MetaData$sgRNA)

# If not, make the control the 1st factor level
cbx2.1MetaData$sgRNA <- relevel(cbx2.1MetaData$sgRNA, "sgLUC")
levels(cbx2.1MetaData$sgRNA)

cbx2.1DataSet <- DESeqDataSetFromMatrix(countData = cbx2.1Counts, # count matrix
                                      colData = cbx2.1MetaData, # samples&grouping
                                      design = ~ sgRNA)

cbx2.1DataSet

# Pre-filtering
nrow(cbx2.1DataSet)

keep <- rowSums(counts(cbx2.1DataSet)) > 1 # rows with a sum of counts > 1
cbx2.1DataSet <- cbx2.1DataSet[keep,] # only rows matching the criteria

nrow(cbx2.1DataSet)

# Create DESeq Data Set object
cbx2.1DESeq <- DESeq(cbx2.1DataSet)
cbx2.1DESeq

# Build the results table 
cbx2.1Results <- results(cbx2.1DESeq,
                       alpha = 0.05, 
                       lfcThreshold = 0.6) 
cbx2.1Results # the comparison should be TREATED vs CONTROL

summary(cbx2.1Results)
table(cbx2.1Results$padj < 0.05)

# Plotting differential expression results 

plotMA(cbx2.1Results) # Create MA plot of differential expression

# # Write results to comma-delimited file or tab-delimited text file
# write.csv(as.data.frame(cbx2.1Results), 
#           file="cbx2.1Results.csv")
# write.table(as.data.frame(cbx2.1Results),
#             sep = "\t",
#             file="cbx2.1Results.txt")


#----------
# Getting list of differentially expressed genes 
#----------

# Get back the data used for plotting as a data frame with DE data
cbx2.1DE <- plotMA(cbx2.1Results, returnData = TRUE)
cbx2.1DE %>% head() # note that there is a column showing which gene is DE

# Get DESeq results as data.frame with rownames as the first column
cbx2.1ResultsDF <- results(cbx2.1DESeq, tidy = TRUE)
cbx2.1ResultsDF %>% head()

# Drop columns 4, 5, and 6
cbx2.1ResultsDF <- cbx2.1ResultsDF[,-c(4:6)]
cbx2.1ResultsDF %>% head()

# Check if both data frames have the same number of rows
nrow(cbx2.1ResultsDF) == nrow(cbx2.1DE) # should return TRUE

# Combine by columns (cbind)
cbx2.1ResultsDF <- cbind(cbx2.1ResultsDF, cbx2.1DE[,2:3])
cbx2.1ResultsDF %>% head()

# Drop redundant lfc column
cbx2.1ResultsDF <- cbx2.1ResultsDF[,-5]
cbx2.1ResultsDF %>% head()
cbx2.1ResultsDF$row %>% unique() %>% length() # number of genes

# Make classes of differential expression
cbx2.1ResultsDF$isDE <- "Not DE"
# if log2Foldchange > 0.5 and FDR < 0.05, set as "UP" 
cbx2.1ResultsDF$isDE[cbx2.1ResultsDF$log2FoldChange > 0.6 & cbx2.1ResultsDF$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.5 and FDR < 0.05, set as "DOWN"
cbx2.1ResultsDF$isDE[cbx2.1ResultsDF$log2FoldChange < -0.6 & cbx2.1ResultsDF$padj < 0.05] <- "Downregulated"

cbx2.1ResultsDF %>% head()

#----------
# Annotate peaks to differential expression data
#----------

# Merge annotated peaks and gene expression data by gene name
cbx2.1ResultsDF %>% nrow() 
cbx2.1ResultsDF$row %>% unique() %>% length()# number of genes from DE analysis

cbx2.1PeakExpression <- merge(CBX2_geneName,
                            cbx2.1ResultsDF,
                            by.x = "Gene.Name",
                            by.y = "row",
                            sort = FALSE)
cbx2.1PeakExpression %>% nrow()
cbx2.1PeakExpression %>% head()

#----------
# Volcano Plot
#----------

# Add column for labels of differentially expressed genes
cbx2.1PeakExpression$Label <- NA

cbx2.1PeakExpression$Label[cbx2.1PeakExpression$isDE != "Not DE"] <- cbx2.1PeakExpression$Gene.Name[cbx2.1PeakExpression$isDE != "Not DE"]

cbx2.1PeakExpression %>% head()

# Add lost/retained peaks category
cbx2.1PeakExpression$Peak <- ""
cbx2.1PeakExpression %>% head()
cbx2.1PeakExpression$Peak[cbx2.1PeakExpression$Gene.Name %in% CBX2_retained_uniqueGenes$Gene.Name] <- "Retained CBX2 peak"
cbx2.1PeakExpression$Peak[cbx2.1PeakExpression$Gene.Name %in% CBX2_lost_uniqueGenes$Gene.Name] <- "Lost CBX2 peak"

# Remove genes not annotated to either lost or retained peaks
cbx2.1PeakExpression <- filter(cbx2.1PeakExpression, Peak != "")
cbx2.1PeakExpression%>% tail()

# Make volcano plot

volcanoPlotcbx2.1 <- ggplot(data=cbx2.1PeakExpression, 
                          aes(x=log2FoldChange, 
                              y=-log10(padj),
                              color = isDE, 
                              shape = Peak, 
                              label =  Label)) + 
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values = c(1,16)) +
  ylab("log10(adjusted p-value)") +
  xlab("log2(Fold-Change)") +
  ggtitle("sgCBX2.3 (1uM EZH2i) vs sgLUC (1uM EZH2i)") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

volcanoPlotcbx2.1 + geom_text_repel() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  # Change point color 
  # By default, it is assigned to the categories in an alphabetical order):
  scale_color_manual(values=c("black", "gray", "orange")) 

#------------------------------------------------------------------------------#
#            sgCBX2 (1 uM EHZ2i/DMSO) VS sgLUC (1 uM EZH2i/DMSO)               #
#------------------------------------------------------------------------------#

cbx2LUCCounts %>% head() # check the loaded count matrix
cbx2.1_vDMSOCounts <- cbx2LUCCounts[,-c(9:12)]
colnames(cbx2.1_vDMSOCounts) # make sure columns correctly correspond to samples' names

# Create the meta-data data frame with the comparative structure of the data.                     

cbx2.1_vDMSOMetaData <- data.frame(row.names = colnames(cbx2.1_vDMSOCounts), # samples names
                             treatment = factor(c("DMSO","DMSO","1uM","1uM")),
                             sgRNA = factor(c("sgLUC", 
                                              "sgLUC",
                                              "sgLUC",
                                              "sgLUC", 
                                              "sgCBX2.1",
                                              "sgCBX2.1",
                                              "sgCBX2.1",
                                              "sgCBX2.1")))

cbx2.1_vDMSOMetaData

# Check the factor levels of the condition. The 1st level should be the control.
levels(cbx2.1_vDMSOMetaData$sgRNA)
levels(cbx2.1_vDMSOMetaData$treatment)

# If not, make the control the 1st factor level
cbx2.1_vDMSOMetaData$sgRNA <- relevel(cbx2.1_vDMSOMetaData$sgRNA, "sgLUC")
levels(cbx2.1_vDMSOMetaData$sgRNA)

cbx2.1_vDMSOMetaData$treatment <- relevel(cbx2.1_vDMSOMetaData$treatment, "DMSO")
levels(cbx2.1_vDMSOMetaData$treatment)

# Create DESeq Data Set object
cbx2.1_vDMSODataSet <- DESeqDataSetFromMatrix(countData = cbx2.1_vDMSOCounts, # count matrix
                                        colData = cbx2.1_vDMSOMetaData, # samples&grouping
                                        design = ~ sgRNA:treatment + treatment * sgRNA)

cbx2.1_vDMSODataSet

#----------
# Exploratory analysis
#----------

# Pre-filtering
nrow(cbx2.1_vDMSODataSet)

keep <- rowSums(counts(cbx2.1_vDMSODataSet)) > 1 # rows with a sum of counts > 1

cbx2.1_vDMSODataSet <- cbx2.1_vDMSODataSet[keep,] # only rows matching the criteria
nrow(cbx2.1_vDMSODataSet)

# PCA plot

rlog_sgRNA_treatment <- rlog(cbx2.1_vDMSODataSet, blind = TRUE)

plotPCA(rlog_sgRNA_treatment, intgroup = c("treatment", "sgRNA"))+
  aes(color = treatment, # color according to treatment condition
      shape = sgRNA) + # shape points according to knockdown condition
  ggtitle(label = "PCA (~ sgRNA * treatment)", 
          subtitle = "Top 500 genes with highest variance") +
  theme_minimal()

#----------
# Differential expression analysis
#----------

cbx2.1_vDMSODESeq <- DESeq(cbx2.1_vDMSODataSet)
resultsNames(cbx2.1_vDMSODESeq)
cbx2.1_vDMSODESeq


# Treatment effect *plus* knockdown effect
cbx2.1_vDMSOResults_A <- results(cbx2.1_vDMSODESeq, 
                               contrast = list(c("sgRNAsgCBX2.1.treatment1uM",
                                                 "sgRNA_sgCBX2.1_vs_sgLUC")),
                               alpha = 0.05, 
                               lfcThreshold = 0.6) 
cbx2.1_vDMSOResults_A # the comparison should be TREATED vs CONTROL

summary(cbx2.1_vDMSOResults_A)
table(cbx2.1_vDMSOResults_A$padj < 0.05)

# Is the treatment effect *different* upon knockdown?
cbx2.1_vDMSOResults_B <- results(cbx2.1_vDMSODESeq, 
                               name = "sgRNAsgCBX2.1.treatment1uM",
                               alpha = 0.05, 
                               lfcThreshold = 0.6) 
cbx2.1_vDMSOResults_B # the comparison should be TREATED vs CONTROL

summary(cbx2.1_vDMSOResults_B)
table(cbx2.1_vDMSOResults_B$padj < 0.05)

cbx2.1_vDMSOResults_C <- results(cbx2.1_vDMSODESeq, 
                                 name = "sgRNAsgCBX2.1.treatment1uM",
                                 alpha = 0.05, 
                                 lfcThreshold = 0.6) 
cbx2.1_vDMSOResults_C # the comparison should be TREATED vs CONTROL

summary(cbx2.1_vDMSOResults_C)
table(cbx2.1_vDMSOResults_C$padj < 0.05)

# Plotting differential expression results 

plotMA(cbx2.1_vDMSOResults_A) # Create MA plot of differential expression
plotMA(cbx2.1_vDMSOResults_B) # Create MA plot of differential expression

# # Write results to comma-delimited file or tab-delimited text file
# write.csv(as.data.frame(cbx2.1_vDMSOResults), 
#           file="cbx2.1_vDMSOResults.csv")
# write.table(as.data.frame(cbx2.1_vDMSOResults),
#             sep = "\t",
#             file="cbx2.1_vDMSOResults.txt")


#----------
# Getting list of differentially expressed genes 
#----------

# Get back the data used for plotting as a data frame with DE data
cbx2.1_vDMSODE <- plotMA(cbx2.1_vDMSOResults_B, returnData = TRUE)
cbx2.1_vDMSODE %>% head() # note that there is a column showing which gene is DE

# Get DESeq results as data.frame with rownames as the first column
# Treatment effect *plus* the treatment effect in knockdown compared to control sgRNA
cbx2.1_vDMSOResultsDF <- results(cbx2.1_vDMSODESeq, 
                                 contrast = list(c("treatment_1uM_vs_DMSO", 
                                                   "sgRNAsgCBX2.1.treatment1uM")),
                                 alpha = 0.05, 
                                 lfcThreshold = 0.6, 
                                 tidy = T) 
cbx2.1_vDMSOResultsDF %>% head()


# Is the treatment effect *different* upon knockdown?
cbx2.1_vDMSOResultsDF <- results(cbx2.1_vDMSODESeq, 
                                 name = "sgRNAsgCBX2.1.treatment1uM",
                                 alpha = 0.05, 
                                 lfcThreshold = 0.6, 
                                 tidy = T)
cbx2.1_vDMSOResultsDF %>% head()

# Drop columns 4, 5, and 6
cbx2.1_vDMSOResultsDF <- cbx2.1_vDMSOResultsDF[,-c(4:6)]
cbx2.1_vDMSOResultsDF %>% head()

# Check if both data frames have the same number of rows
nrow(cbx2.1_vDMSOResultsDF) == nrow(cbx2.1_vDMSODE) # should return TRUE

# Combine by columns (cbind)
cbx2.1_vDMSOResultsDF <- cbind(cbx2.1_vDMSOResultsDF, cbx2.1_vDMSODE[,2:3])
cbx2.1_vDMSOResultsDF %>% head()

# Drop redundant lfc column
cbx2.1_vDMSOResultsDF <- cbx2.1_vDMSOResultsDF[,-5]
cbx2.1_vDMSOResultsDF %>% head()
cbx2.1_vDMSOResultsDF$row %>% unique() %>% length() # number of genes

# Make classes of differential expression
cbx2.1_vDMSOResultsDF$isDE <- "Not DE"
# if log2Foldchange > 0.5 and FDR < 0.05, set as "UP" 
cbx2.1_vDMSOResultsDF$isDE[cbx2.1_vDMSOResultsDF$log2FoldChange > 0.6 & cbx2.1_vDMSOResultsDF$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.5 and FDR < 0.05, set as "DOWN"
cbx2.1_vDMSOResultsDF$isDE[cbx2.1_vDMSOResultsDF$log2FoldChange < -0.6 & cbx2.1_vDMSOResultsDF$padj < 0.05] <- "Downregulated"

cbx2.1_vDMSOResultsDF %>% head()

#----------
# Annotate peaks to differential expression data
#----------

# Merge annotated peaks and gene expression data by gene name
cbx2.1_vDMSOResultsDF %>% nrow() 
cbx2.1_vDMSOResultsDF$row %>% unique() %>% length()# number of genes from DE analysis

cbx2.1_vDMSOPeakExpression <- merge(CBX2_geneName,
                              cbx2.1_vDMSOResultsDF,
                              by.x = "Gene.Name",
                              by.y = "row",
                              sort = FALSE)
cbx2.1_vDMSOPeakExpression %>% nrow()
cbx2.1_vDMSOPeakExpression %>% head()

#----------
# Volcano Plot
#----------

# Add column for labels of differentially expressed genes
cbx2.1_vDMSOPeakExpression$Label <- NA

cbx2.1_vDMSOPeakExpression$Label[cbx2.1_vDMSOPeakExpression$isDE != "Not DE"] <- cbx2.1_vDMSOPeakExpression$Gene.Name[cbx2.1_vDMSOPeakExpression$isDE != "Not DE"]

cbx2.1_vDMSOPeakExpression %>% head()

# Add lost/retained peaks category
cbx2.1_vDMSOPeakExpression$Peak <- ""
cbx2.1_vDMSOPeakExpression %>% head()
cbx2.1_vDMSOPeakExpression$Peak[cbx2.1_vDMSOPeakExpression$Gene.Name %in% CBX2_retained_uniqueGenes$Gene.Name] <- "Retained CBX2 peak"
cbx2.1_vDMSOPeakExpression$Peak[cbx2.1_vDMSOPeakExpression$Gene.Name %in% CBX2_lost_uniqueGenes$Gene.Name] <- "Lost CBX2 peak"

# Remove genes not annotated to either lost or retained peaks
cbx2.1_vDMSOPeakExpression <- filter(cbx2.1_vDMSOPeakExpression, Peak != "")
cbx2.1_vDMSOPeakExpression%>% tail()

# Make volcano plot

volcanoPlotcbx2.1_vDMSO <- ggplot(data=cbx2.1_vDMSOPeakExpression, 
                            aes(x=log2FoldChange, 
                                y=-log10(padj),
                                color = isDE, 
                                shape = Peak, 
                                label =  Label)) + 
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values = c(1,16)) +
  ylab("log10(adjusted p-value)") +
  xlab("log2(Fold-Change)") +
  ggtitle("sgCBX2.1 (1uM EZH2i) vs control sgRNA (DMSO)") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

volcanoPlotcbx2.1_vDMSO + geom_text_repel() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  # Change point color 
  # By default, it is assigned to the categories in an alphabetical order):
  scale_color_manual(values=c("black", "gray", "orange")) 

#------------------------------------------------------------------------------#
#                  sgCBX2 (1 uM EHZ2i) VS sgLUC (DMSO)                   #
#------------------------------------------------------------------------------#

cbx2LUCCounts %>% head() # check the loaded count matrix
cbx2.1Counts <- cbx2LUCCounts[,-c(3:6,9:12)]
colnames(cbx2.1Counts) # make sure columns correctly correspond to samples' names

# Create the meta-data data frame with the comparative structure of the data.                     

cbx2.1MetaData <- data.frame(row.names = colnames(cbx2.1Counts), # samples names
                             treatment = factor(c("DMSO","DMSO","1uM","1uM")),
                             sgRNA = factor(c("sgLUC", 
                                              "sgLUC", 
                                              "sgCBX2.1",
                                              "sgCBX2.1")))

cbx2.1MetaData

# Check the factor levels of the condition. The 1st level should be the control.
levels(cbx2.1MetaData$sgRNA)

# If not, make the control the 1st factor level
cbx2.1MetaData$sgRNA <- relevel(cbx2.1MetaData$sgRNA, "sgLUC")
levels(cbx2.1MetaData$sgRNA)
cbx2.1MetaData$treatment <- relevel(cbx2.1MetaData$treatment, "DMSO")
levels(cbx2.1MetaData$treatment)

cbx2.1DataSet <- DESeqDataSetFromMatrix(countData = cbx2.1Counts, # count matrix
                                        colData = cbx2.1MetaData, # samples&grouping
                                        design = ~ sgRNA)

cbx2.1DataSet

# Pre-filtering
nrow(cbx2.1DataSet)

keep <- rowSums(counts(cbx2.1DataSet)) > 1 # rows with a sum of counts > 1
cbx2.1DataSet <- cbx2.1DataSet[keep,] # only rows matching the criteria

nrow(cbx2.1DataSet)

# Create DESeq Data Set object
cbx2.1DESeq <- DESeq(cbx2.1DataSet)
cbx2.1DESeq
resultsNames(cbx2.1DESeq)
# Build the results table 
cbx2.1Results <- results(cbx2.1DESeq,
                         alpha = 0.05, 
                         lfcThreshold = 0.6) 
cbx2.1Results # the comparison should be TREATED vs CONTROL

summary(cbx2.1Results)
table(cbx2.1Results$padj < 0.05)

# Plotting differential expression results 

plotMA(cbx2.1Results) # Create MA plot of differential expression

# # Write results to comma-delimited file or tab-delimited text file
# write.csv(as.data.frame(cbx2.1Results), 
#           file="cbx2.1Results.csv")
# write.table(as.data.frame(cbx2.1Results),
#             sep = "\t",
#             file="cbx2.1Results.txt")


#----------
# Getting list of differentially expressed genes 
#----------

# Get back the data used for plotting as a data frame with DE data
cbx2.1DE <- plotMA(cbx2.1Results, returnData = TRUE)
cbx2.1DE %>% head() # note that there is a column showing which gene is DE

# Get DESeq results as data.frame with rownames as the first column
cbx2.1ResultsDF <- results(cbx2.1DESeq, tidy = TRUE)
cbx2.1ResultsDF %>% head()

# Drop columns 4, 5, and 6
cbx2.1ResultsDF <- cbx2.1ResultsDF[,-c(4:6)]
cbx2.1ResultsDF %>% head()

# Check if both data frames have the same number of rows
nrow(cbx2.1ResultsDF) == nrow(cbx2.1DE) # should return TRUE

# Combine by columns (cbind)
cbx2.1ResultsDF <- cbind(cbx2.1ResultsDF, cbx2.1DE[,2:3])
cbx2.1ResultsDF %>% head()

# Drop redundant lfc column
cbx2.1ResultsDF <- cbx2.1ResultsDF[,-5]
cbx2.1ResultsDF %>% head()
cbx2.1ResultsDF$row %>% unique() %>% length() # number of genes

# Make classes of differential expression
cbx2.1ResultsDF$isDE <- "Not DE"
# if log2Foldchange > 0.5 and FDR < 0.05, set as "UP" 
cbx2.1ResultsDF$isDE[cbx2.1ResultsDF$log2FoldChange > 0.6 & cbx2.1ResultsDF$padj < 0.05] <- "Upregulated"
# if log2Foldchange < -0.5 and FDR < 0.05, set as "DOWN"
cbx2.1ResultsDF$isDE[cbx2.1ResultsDF$log2FoldChange < -0.6 & cbx2.1ResultsDF$padj < 0.05] <- "Downregulated"

cbx2.1ResultsDF %>% head()

#----------
# Annotate peaks to differential expression data
#----------

# Merge annotated peaks and gene expression data by gene name
cbx2.1ResultsDF %>% nrow() 
cbx2.1ResultsDF$row %>% unique() %>% length()# number of genes from DE analysis

cbx2.1PeakExpression <- merge(CBX2_geneName,
                              cbx2.1ResultsDF,
                              by.x = "Gene.Name",
                              by.y = "row",
                              sort = FALSE)
cbx2.1PeakExpression %>% nrow()
cbx2.1PeakExpression %>% head()

#----------
# Volcano Plot
#----------

# Add column for labels of differentially expressed genes
cbx2.1PeakExpression$Label <- NA

cbx2.1PeakExpression$Label[cbx2.1PeakExpression$isDE != "Not DE"] <- cbx2.1PeakExpression$Gene.Name[cbx2.1PeakExpression$isDE != "Not DE"]

cbx2.1PeakExpression %>% head()

# Add lost/retained peaks category
cbx2.1PeakExpression$Peak <- ""
cbx2.1PeakExpression %>% head()
cbx2.1PeakExpression$Peak[cbx2.1PeakExpression$Gene.Name %in% CBX2_retained_uniqueGenes$Gene.Name] <- "Retained CBX2 peak"
cbx2.1PeakExpression$Peak[cbx2.1PeakExpression$Gene.Name %in% CBX2_lost_uniqueGenes$Gene.Name] <- "Lost CBX2 peak"

# Remove genes not annotated to either lost or retained peaks
cbx2.1PeakExpression <- filter(cbx2.1PeakExpression, Peak != "")
cbx2.1PeakExpression%>% tail()

# Make volcano plot

volcanoPlotcbx2.1 <- ggplot(data=cbx2.1PeakExpression, 
                            aes(x=log2FoldChange, 
                                y=-log10(padj),
                                color = isDE, 
                                shape = Peak, 
                                label =  Label)) + 
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values = c(1,16)) +
  ylab("log10(adjusted p-value)") +
  xlab("log2(Fold-Change)") +
  ggtitle("CBX2 knockdown (1uM EZH2i) vs control sgRNA (DMSO)") +
  theme(legend.title = element_blank(), 
        axis.title.x = element_text(face = "bold"), 
        axis.title.y = element_text(face = "bold"))

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

volcanoPlotcbx2.1 + geom_text_repel() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  # Change point color 
  # By default, it is assigned to the categories in an alphabetical order):
  scale_color_manual(values=c("black", "gray", "orange")) 
