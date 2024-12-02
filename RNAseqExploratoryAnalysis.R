
########################################
# EXPLORATORY ANALYSIS OF RNA-seq DATA #
########################################

#----------
# Correlation Heatmap 
#----------

library("PoiClaClu") 

# See that in the count matrix the sample names are columns
counts(dataSet) %>% head() 
counts(dataSet) %>% nrow()

# We need to transpose the matrix of values using t(), because both the dist   
# and the PoissonDistance functions expect the different samples to be rows of 
# its argument, and different dimensions (here, genes) to be columns. 
sampleDistances <- PoissonDistance(x = t(counts(dataSet)), # transpose matrix
                                   type = "deseq") 

sampleDistances  
sampleDistances$dd # the n x n dissimilarity matrix

# We then visualize the distances using the the pheatmap package
library("pheatmap")
library("RColorBrewer")

# Get the dissimilarity matrix
sampleDistMatrix <- as.matrix(sampleDistances$dd)

# Tidy up samples' names
rownames(sampleDistMatrix) <- paste(dataSet$sgRNA, 
                                    dataSet$treatment, 
                                    sep = " / " )

# Omit gene names
colnames(sampleDistMatrix) <- NULL

# Define the colors for the heatmap
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Plot heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDistances$dd,
         clustering_distance_cols = sampleDistances$dd,
         col = colors)



#----------
# PCA Plot
#----------

rlogDataSet <- rlog(dataSet, 
                    blind = F # TRUE to compare samples in an unbiased way, as in Quality Assurance
                    )

plotPCA(rlogDataSet, intgroup = c("treatment", "sgRNA"))+
  aes(color = treatment, # color according to treatment condition
      shape = sgRNA) + # shape points according to experiment
  ggtitle(label = "Principal Component Analysis", 
          subtitle = "Top 500 genes with highest variance") +
  theme_bw()

#----------
# GLM PCA plot
#----------
library("glmpca")

?glmpca # check what the arguments for the function mean

GPCA <- glmpca(counts(dataSet), # use the count matrix in DESeqDataSet 
               L=2, # specify two dimensions for reduction
               fam = "nb", nb_theta = 100) # negative binomial as in DESEq2

#GPCA # summary of the model
#GPCA$factors # matrix whose rows match the columns of the count matrix and
# whose columns are the different latent dimensions. Analogous
# to the principal components in PCA model
#GPCA$loadings # matrix with rows matching the rows (features/genes) of the
# count matrix. Columns are different latent dimensions

gpca.data <- GPCA$factors # extract principal components from the GPCA model
gpca.data$treatment <- dataSet$treatment
gpca.data$genotype <- dataSet$genotype
gpca.data$sgRNA <- dataSet$sgRNA

ggplot(gpca.data, aes(x = dim1, 
                      y = dim2, 
                      color = treatment, 
                      #shape = genotype, 
                      shape = sgRNA)) +
  geom_point(size = 3) + 
  coord_fixed() + 
  ggtitle("Generalized Principal Component Analysis") +
  theme_bw()

#----------
# Scatterplot of transformed counts from two samples
#----------

vsd <- vst(dataSet, blind = FALSE)

a <- as_tibble(log2(counts(dataSet, normalized=TRUE)[, c(5,9)]+1)) %>%
  mutate(transformation = "log2(x + 1)")
b <- as_tibble(assay(vsd)[, c(5,9)]) %>% mutate(transformation = "vst")
c <- as_tibble(assay(rlogDataSet)[, c(5,9)]) %>% mutate(transformation = "rlog")

df <- bind_rows(a,b,c)
  
colnames(df)[1:2] <- c("sgCBX2.1/DMSO", "sgCBX2.3/DMSO")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = "sgCBX2.1/DMSO", y = "sgCBX2.3/DMSO")) + geom_hex(bins = 100) +
  coord_fixed() + facet_grid( . ~ transformation)
