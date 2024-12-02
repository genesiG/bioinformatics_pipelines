###############################################################################
#                       PEAK ANNOTATION WITH CHIPSEEKER                       #
###############################################################################

#####################
# PREPARE WORKSPACE #
#####################

# Check if the package is installed by trying to load it with require()
# require() returns a logical (TRUE or FALSE) depending on if it's able to load 
# the package. If it fails, install the package.

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  }

if(!require(ChIPseeker)){
  BiocManager::install("ChIPseeker")
  }

# Load required packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(rtracklayer)

# Set your working directory
# Change for each Project
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN")

#############
# LOAD DATA #
#############

# Provide the names of the BED files as a list
samplefiles <- list.files(path = "annotatePeaks", # name of folder
                          pattern = ".broadPeak", # file extension
                          full.names = TRUE)

samplefiles <- as.list(samplefiles)

################################################################################
#                                                                              #
# If the files are in a sub-folder within your working directory, then the     #
# "path" argument is just the name of folder, otherwise it's the full path     #
# name to the folder as a character vector.                                    #
#                                                                              #
# The "pattern" argument will fetch all file names within the path that have   #
# that file extension.                                                         #
#                                                                              #
################################################################################

# Name each element in the list
names(samplefiles) <- c("CBX2/1uM EZH2i", "CBX2/5uM EZH2i", "CBX2/DMSO")

# Check that each file is correctly named
print(samplefiles) 

# Load peaks as GRanges objects 
peakList <- lapply(samplefiles,          # apply function to this list
                   import,               # function to be applied
                   format = "broadPeak") # argument from the import() function


##############
# ANNOTATION #
##############

# Load database of known genes from the genome you used to align reads
hg19genes <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Use the annoatePeak function to annotate the list of peaks
annotatedPeaksList <- lapply(peakList,         # apply function to this list
                             annotatePeak,     # function to be applied
                             TxDb = hg19genes, # arguments from annotatePeak() 
                             tssRegion = c(-1000, 1000), 
                             verbose = FALSE)

# Plot pie chart of annotated peaks
plotAnnoPie(annotatedPeaksList$`CBX2/DMSO`, # use $ to specify which peak set
            main = "CBX2 (DMSO)",
            cex = 1.4) # legend size


###############################################################################
#                                VENN DIAGRAMS                                #
###############################################################################

###############
# IMPORT DATA #
###############

# Import broadPeak files as GRanges object
library(rtracklayer)
CBX2_DMSOpeaks <- import("./peaks/CBX2_DMSO.broadPeak",format = "broadPeak")

# Alternatively,
CBX2_1uMpeaks <- peakList$`CBX2/1uM EZH2i`
CBX2_5uMpeaks <- peakList$`CBX2/5uM EZH2i`

################################################################################
#                                                                              #  
# .broadPeak, .narrowPeak, and .gappedPeak files are types of BED files used by#
# the ENCODE project to provide called regions of signal enrichment based on   #
# pooled, normalized (interpreted) data. MACS2 outputs peaks as tab-separated  # 
# files.                                                                       #
#                                                                              #
# When importing peaks using read.delim(), the default for the col.names       #
# (column names) parameter is to use "V" followed by the column number.        #
#                                                                              #
# See below what each column means:                                            #
#                                                                              #
#                                                                              #
# V1  - chrom       - Name of the chromosome (or contig, scaffold, etc.).      #
# V2  - chromStart  - The starting position of the feature in the chromosome   # 
#                     or scaffold. The 1st base in a chromosome is numbered 0. #
# V3  - chromEnd    - The ending position of the feature in the chromosome or  #
#                     scaffold.                                                #
# V4  - name        - Name given to a region (preferably unique). Use "." if   # 
#                     no name is assigned.                                     #      
# V5  - score       - Indicates how dark the peak will be displayed in the     #
#                     genome browser or IGV (0-1000). Ideally the average      #
#                     signalValue per base spread is between 100-1000.         #
# V6  - strand      - +/- to denote strand or orientation (whenever applicable)#
#                     Use "." if no orientation is assigned.                   #
# V7  - signalValue - Measurement of overall (usually, average) enrichment for #
#                     the region.                                              #
# V8  - pValue      - Measurement of statistical significance (-log10). Use -1 #
#                     if no pValue is assigned.                                #
# V9  - qValue      - Measurement of statistical significance using false      #  
#                     discovery rate (-log10). Use -1 if no qValue is assigned.#
# V10 - peak        - (for narrowPeak files only) Point-source called for this #
#                     peak; 0-based offset from chromStart. Use -1 if no       # 
#                     point-source called.                                     #
#                                                                              #
################################################################################

#####################
# PLOT VENN DIAGRAM #
#####################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChIPpeakAnno")
library(ChIPpeakAnno)

makeVennDiagram(
  Peaks = peakList,
  NameOfPeaks = names(peakList),
  connectedPeaks = "min",
  plot = TRUE, 
  cat.fontfamily = "sans",
  cat.fontface = "bold",
  lty = "blank",
  fill = c("darkgreen", "yellow", "skyblue"),
  alpha = 0.55, 
  print.mode = c("raw", "percent")
)


################################################################################

library(dplyr)
library(rtracklayer)
CBX2_DMSOsharp <- import.bed("./peaks/CBX2_hclust2.mean.sharp.broadPeak")

CBX2_DMSOsharp %>% head()

CBX2_DMSObroad <- import.bed("./peaks/CBX2_hclust2.mean.broad.broadPeak")

CBX2_DMSObroad %>% head()

CBX2_DMSOpeaks <- import("./peaks/CBX2_min1k.b001.q1E5.final.broadPeak",
                         format = "broadPeak")


# Get peak widths
CBX2_DMSO_sharpWidth <- CBX2_DMSOsharp %>% ranges() %>% as.data.frame()
CBX2_DMSO_sharpWidth <- CBX2_DMSO_sharpWidth$width
CBX2_DMSO_sharpWidth %>% median() # median width of peaks
CBX2_DMSO_broadWidth <- CBX2_DMSObroad %>% ranges() %>% as.data.frame()
CBX2_DMSO_broadWidth <- CBX2_DMSO_broadWidth$width
CBX2_DMSO_broadWidth %>% median() # median width of peaks

CBX2_DMSOpeakWidth <- CBX2_DMSOpeaks %>% ranges() %>% as.data.frame()
CBX2_DMSOpeakWidth %>% head()
CBX2_DMSOpeakWidth <- CBX2_DMSOpeakWidth[3]
CBX2_DMSO_Widths <- t(CBX2_DMSOpeakWidth)
CBX2_DMSOpeakWidth %>% head()
CBX2_DMSOpeakWidth %>% median() # median width of peaks
CBX2_DMSOpeakWidth %>% mean() # mean width of peaks

# Calculate density
densityBroad <- density(CBX2_DMSO_broadWidth, 
                        n = length(CBX2_DMSO_broadWidth),
                        adjust = 1)
densitySharp <- density(CBX2_DMSO_sharpWidth, 
                        n = length(CBX2_DMSO_sharpWidth),
                        adjust = 1)

# Check data Normality
?qqnorm
qqnorm(CBX2_DMSO_sharpWidth, pch = 16, 
       col = rgb(0, 0, 0, alpha = 0.5)) #transparent grey
qqline(CBX2_DMSO_sharpWidth, 
       col = "red")

# Generalized PCA of peak widths
### for dimension reduction of non-normally distributed data.

library(glmpca)
?glmpca
peakWidthGPCA <- glmpca(CBX2_DMSO_Widths, # data frame of peak widths
                        L=2) # specify two dimensions for reduction

peakWidthGPCA # summary of the model
peakWidthGPCA$factors # matrix whose rows match the columns of the count matrix and
# whose columns are the different latent dimensions. Analogous
# to the principal components in PCA model
peakWidthGPCA$loadings # matrix with rows matching the rows (features/genes) of the
# count matrix. Columns are different latent dimensions

gpca.width <- peakWidthGPCA$factors # extract principal components from the GPCA model
gpca.width$width <- peakWidthGPCA$
gpca.width$experiment <- peakWidthGPCA$experiment

ggplot(gpca.width, aes(x = dim1, 
                      y = dim2, 
                      #color = treatment, 
                      #shape = experiment
                      )) +
  geom_point(size =3) + 
  coord_fixed() + 
  ggtitle("Generalized Principal Component Analysis") +
  theme_minimal()


################################################################################

# Detecting and removing outliers

# Plot histogram using the square root of the number of observations as the
# number of bins
hist(CBX2_DMSOpeakWidth$width,
     xlab = "Peak Width (bp)",
     main = "Histogram of Peak Widths ",
     breaks = sqrt(nrow(CBX2_DMSOpeakWidth))) # set number of bins

# Or detecting through a boxplot:
boxplot(CBX2_DMSOpeakWidth$width,
        xlab = "Peak Width (bp)")

# Observations considered as potential outliers by the interquartile range (IQR)
# criterion are displayed as points in the boxplot.

####
####


####

# This method of outliers detection is based on the percentiles. With the 
# percentiles method, all observations that lie outside the interval formed by 
# the 2.5 and 97.5 percentiles will be considered as potential outliers. Other 
# percentiles such as the 1 and 99, or the 5 and 95 percentiles can also be 
# considered to construct the interval.
 
# The values of the lower and upper percentiles (and thus the lower and upper 
# limits of the interval) can be computed with the quantile() function:
  
lower_bound <- quantile(CBX2_DMSOpeakWidth$width, 0.05)
lower_bound

upper_bound <- quantile(CBX2_DMSOpeakWidth$width, 0.95)
upper_bound

# According to this method, all observations below lower_boubd and above 
# upper_bound will be considered as potential outliers. The row numbers of the 
# observations outside of the interval can then be extracted using which()
#   

# Find which indexed (row numbers) correspond to outliers
outlier_ind <- which(CBX2_DMSOpeakWidth$width < lower_bound | CBX2_DMSOpeakWidth$width > upper_bound)

outlier_ind %>% head()

# Filter out outliers by matching row numbers
CBX2_peakWidthOutRm <- CBX2_DMSOpeakWidth %>% 
  mutate(rn = row_number()) %>% 
  filter(!rn %in% outlier_ind)

ggplot(CBX2_peakWidthOutRm, 
       aes(x= log(width))) + 
  geom_density(alpha = 0.25, 
               color = "red",
               fill = "red",
               size = 0.8,
               n = nrow(CBX2_peakWidthOutRm)) +
  #scale_x_continuous(limits = c(0,40000)) +
  theme_minimal() +
  labs(title = "CBX2 peaks (DMSO)") +
  xlab("Peak width (bp)") +
  ylab("Density") +
  theme(axis.title.y = element_text(face = "bold",),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", 
                             hjust = 0.5)) 

cut_clusters <- cutree(hclustOutRm, k = 2)

CBX2_DMSOpeakCluster <- mutate(CBX2_peakWidthOutRm, 
                               cluster = cut_clusters)

count(CBX2_DMSOpeakCluster, cluster) # see how many assigned to each cluster

CBX2_DMSOpeakCluster %>% filter(cluster == 2) %>% head()

CBX2_DMSOpeakCluster %>% filter(cluster == 1) %>% summary()
CBX2_DMSOpeakCluster %>% filter(cluster == 2) %>% summary()


ggplot(CBX2_DMSOpeakCluster, 
       aes(x= width,
           color = factor(cluster),
           fill = factor(cluster))) + 
  geom_density(alpha = 0.5, size = 1) +
  #scale_x_continuous(limits = c(0,40000)) +
  theme_minimal() +
  labs(title = "CBX2 peaks (DMSO)") +
  xlab("Peak width (bp)") +
  ylab("Density") +
  theme(axis.title.y = element_text(face = "bold",),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", 
                             hjust = 0.5)) 

################################################################################
################################################################################
# hclust of peak Widths
### Matrix of Widths

hclust <- hclust(dist(CBX2_DMSOpeakWidth$width), # euclidean distance matrix
                 method = "ward.D") # use ward linkage
hclustComplete <- hclust(dist(CBX2_DMSOpeakWidth$width), # euclidean distance matrix
                 method = "complete",
                 members = CBX2_DMSOpeakWidth$width) # use complete linkage
hclustCentroid <- hclust(dist(CBX2_DMSOpeakWidth$width), # euclidean distance matrix
                         method = "centroid",
                         members = CBX2_DMSOpeakWidth$width) # centroid linkage
hclustAvg <- hclust(dist(CBX2_DMSOpeakWidth$width), # euclidean distance matrix
                         method = "average",
                         members = CBX2_DMSOpeakWidth$width) # average linkage
hclustWardD2 <- hclust(dist(CBX2_DMSOpeakWidth$width), # euclidean distance matrix
                         method = "ward.D2",
                         members = CBX2_DMSOpeakWidth$width) # use Ward D2 linkage


hclustOutRm <- hclust(dist(CBX2_peakWidthOutRm$width), # euclidean distance matrix
                      method = "ward.D") # use ward linkage

plot(hclustOutRm, 
     main = "CBX2 (DMSO)", 
     xlab = "Individual peak widths (bp)",
     cex = 0.5,
     hang = -1,
     labels = F)

# If you visually want to see the clusters on the dendrogram you can use R's 
# rect.hclust() to superimpose rectangular compartments for each cluster on the 
# tree as shown in the following code:

plot(hclustWardD2,
     main = "CBX2 (DMSO)", 
     xlab = "Individual peak widths (bp)",
     cex = 0.5,
     hang = -1,
     labels = F)

rect.hclust(hclustOutRm, 
            k = 2, # cut the dendogram so that exactly k clusters are produced 
            border = 2:6) # vecotr with border colors for the rectangles

# You can also use the color_branches() function from the dendextend library to 
# visualize your tree with different colored branches.

if (!require("dendextend", quietly = TRUE)){
  install.packages("dendextend")}
library(dendextend)

dendogram <- as.dendrogram(hclust)
dendogramColored <- color_branches(dendogram, k = 2)
plot(dendogramColored,
     main = "CBX2 (DMSO)",
     xlab = "Individual peak widths (bp)",
     cex = 0.5)

# Use the cutree() function to cut the tree with hclust as one parameter and the
# other parameter as k = 2 to cut it into 2 cluster as you see in the dendogram

cut_clusters <- cutree(hclust, k = 2)

CBX2_DMSOpeakCluster <- mutate(CBX2_DMSOpeakWidth, 
                               cluster = cut_clusters)

count(CBX2_DMSOpeakCluster, cluster) # see how many assigned to each cluster

CBX2_DMSOpeakCluster %>% filter(cluster == 2) %>% head()

CBX2_DMSOpeakCluster %>% filter(cluster == 2) %>% summary()
CBX2_DMSOpeakCluster %>% filter(cluster == 1) %>% summary()

library(ggplot2)
ggplot(CBX2_DMSOpeakCluster, 
       aes(x= factor(cluster), 
           y = width)) + 
  geom_boxplot(varwidth = T) +
  scale_y_continuous(limits = c(0,50000)) +
  theme_minimal() +
  ggtitle("CBX2 peaks (DMSO)") +
  xlab("Cluster") +
  ylab("Peak width (bp)") +
  theme(axis.title.y = element_text(face = "bold",),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold")) 

ggplot(CBX2_DMSOpeakCluster, 
       aes(x= width,
           color = factor(cluster),
           fill = factor(cluster))) + 
  geom_density(alpha = 0.5, size = 1) +
  scale_x_continuous(limits = c(0,40000)) +
  theme_minimal() +
  labs(title = "CBX2 peaks (DMSO)") +
  xlab("Peak width (bp)") +
  ylab("Density") +
  theme(axis.title.y = element_text(face = "bold",),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", 
                             hjust = 0.5)) 
  

ggplot(CBX2_DMSOpeakWidth, 
       aes(x= log10(width))) + 
  geom_histogram(color = "black",
                 fill = "white",
                 binwidth = 100) #+
  geom_freqpoly(bins = sqrt(nrow(CBX2_DMSOpeakWidth)),
                size = 1,
                color = "red") +
  # stat_density(alpha = 0.25, 
  #              color = "red",
  #              fill = "red",
  #              size = 0.8,
  #              n = sqrt(nrow(CBX2_DMSOpeakWidth))) +
  scale_x_continuous(minor_breaks = NULL) +
  #scale_y_discrete() +
  theme_minimal() +
  labs(title = "CBX2 (DMSO)") +
  xlab("log10(Peak width)") +
  ylab("Frequency") +
  theme(axis.title.y = element_text(face = "bold",),
        axis.title.x = element_text(face = "bold"),
        title = element_text(face = "bold", 
                             hjust = 0.5)) 








# Plot densities
dx <- density(log10(CBX2_DMSOpeakWidth$width), 
              bw = "nrd0",
              adjust = 1)
plot(dx, main = "CBX2 (DMSO)")

plot(densityBroad, xlim = c(0,10000) )
lines(densitySharp, add = TRUE)

# Add density
lines(dx, lwd = 2, col = "red")


#############################################################################
if(!require(ChIPseeker)){
  BiocManager::install("ChIPseeker")
}

# Load required packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationDbi)
library(rtracklayer)

# Set your working directory
# Change for each Project
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN")

# Import broadPeak files as GRanges object
CBX2_DMSOpeaks <- import("./new_peaks/CBX2_DMSO_final.broadPeak",
                         format = "broadPeak")

CBX2_DMSOpeaks

