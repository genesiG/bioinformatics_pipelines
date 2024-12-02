options(scipen = 999)

plotGenes <- function(counts){
  # Visualize data in a plot
  library(ggplot2)
  
  # Convert counts matrix to a data frame
  counts_df <- as.data.frame(counts)
  counts_df$gene <- rownames(counts_df)
  
  # Reshape data for ggplot
  counts_df_long <- tidyr::gather(counts_df, sample, count, -gene)
  
  # Merge with colData to get condition information
  counts_df_long <- merge(counts_df_long, colData, by.x = "sample", by.y = "row.names")
  
  # Create barplot
  ggplot(counts_df_long, aes(x = gene, y = count, fill = gene)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(y = "Count") +
    facet_wrap(~sample, dir = "v") +
    theme_classic() + 
    guides(fill = "none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())  # Rotate x-axis labels for better readability
}

naiveTest = function(counts){
  gene1_test = t.test(x = counts[1,3:4], 
                      y = counts[1,1:2], 
                      paired = T)
  lfc1 = log2(mean(counts[1,1:2])) - log2(mean(counts[1,3:4])) # LFC
  pval1 = gene1_test$p.value
  
  gene2_test = t.test(x = counts[2,3:4], 
                      y = counts[2,1:2], 
                      paired = T)
  lfc2 = log2(mean(counts[2,1:2])) - log2(mean(counts[2,3:4])) # LFC
  pval2 = gene2_test$p.value
  
  gene3_test = t.test(x = counts[3,3:4], 
                      y = counts[3,1:2], 
                      paired = T)
  lfc3 =log2(mean(counts[3,1:2])) - log2(mean(counts[3,3:4])) # LFC
  pval3 = gene3_test$p.value
  
  if(nrow(counts) == 3){
  df = data.frame(gene = row.names(counts),
             LFC = c(lfc1, lfc2, lfc3),
             p.value = c(pval1,pval2,pval3),
             FDR = c(p.adjust(pval1, method = "fdr", n = nrow(counts)),
                     p.adjust(pval2, method = "fdr", n = nrow(counts)),
                     p.adjust(pval3, method = "fdr", n = nrow(counts))))
  }
  
  if(nrow(counts) == 5){
    gene4_test = t.test(x = counts[4,3:4], 
                        y = counts[4,1:2], 
                        paired = T)
    lfc4 = log2(mean(counts[4,1:2])) - log2(mean(counts[4,3:4])) # LFC
    pval4 = gene4_test$p.value
    
    gene5_test = t.test(x = counts[5,3:4], 
                        y = counts[5,1:2], 
                        paired = T)
    lfc5 =log2(mean(counts[5,1:2])) - log2(mean(counts[5,3:4])) # LFC
    pval5 = gene5_test$p.value
    
    df = data.frame(gene = row.names(counts),
                    LFC = c(lfc1, lfc2, lfc3, lfc4, lfc5),
                    p.value = c(pval1,pval2,pval3,pval4,pval5),
                    FDR = c(p.adjust(pval1, method = "fdr", n = nrow(counts)),
                            p.adjust(pval2, method = "fdr", n = nrow(counts)),
                            p.adjust(pval3, method = "fdr", n = nrow(counts)),
                            p.adjust(pval4, method = "fdr", n = nrow(counts)),
                            p.adjust(pval5, method = "fdr", n = nrow(counts))))
    
  }
  
  print(df)
}

naiveRPKM = function(counts, gene.length = c(500,500,500)){
  
  rpkm_counts = rpkm(counts, 
                     gene.length = gene.length)
  rpkm_counts
  
  naiveTest(rpkm_counts)
}

naiveQuantile = function(counts){
  df_rank <- apply(as.data.frame(counts),2,rank,ties.method="min")
  df_sorted <- data.frame(apply(as.data.frame(counts), 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- row.names(df_rank)
  df_final
}


###########################
# EXAMPLE 1: Only gene 3 DE, different sequencing depths
###########################
# Create counts matrix
                  #30 #30 #18 #18
counts <- matrix(c(6, 6, 6, 6,  
                   5, 5, 5, 5,    
                   18, 18, 6, 6),
                 nrow = 3,
                 byrow = TRUE)
colnames(counts) <- c("Treated_1", "Treated_2", "Untreated_1", "Untreated_2")
rownames(counts) <- c("Gene1", "Gene2", "Gene3")

# Create colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Treated", 2), rep("Untreated", 2))
)

# Print counts matrix
print(counts)

# Visualize gene counts
plotGenes(counts)

# Total Counts norm
tc_counts = sweep(counts, 2, colSums(counts), "/")
tc_counts

#RPKM norm
# Same gene length
gene.length = c(500, 500, 500)
rpkm_counts = rpkm(counts, gene.length = gene.length)
rpkm_counts

# Quantile norm
naiveQuantile(counts)
naiveQuantile(tc_counts)
naiveQuantile(rpkm_counts)

# RPKM
# Unequal gene length
gene.length = c(500, 500, 1500)
rpkm_counts = rpkm(counts, gene.length = gene.length)

rpkm_counts
naiveQuantile(rpkm_counts)

###########################
# EXAMPLE 2: Only Gene 3 DE, same sequencing depth
###########################
# Create counts matrix
                  #18 #18 #18 #18
counts <- matrix(c(3, 3, 6, 6,  
                   3, 3, 6, 6,    
                   12, 12, 6, 6),
                 nrow = 3,
                 byrow = TRUE)
colnames(counts) <- c("Treated_1", "Treated_2", "Untreated_1", "Untreated_2")
rownames(counts) <- c("Gene1", "Gene2", "Gene3")

# Create colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Treated", 2), rep("Untreated", 2))
)

# Print counts matrix
print(counts)

# Visualize gene counts
plotGenes(counts)

# Total Counts norm
tc_counts = sweep(counts, 2, colSums(counts), "/")
tc_counts

#RPKM norm
# Same gene length
gene.length = c(500, 500, 500)
rpkm_counts = rpkm(counts, gene.length = gene.length)
rpkm_counts

# Quantile norm
naiveQuantile(counts)
naiveQuantile(tc_counts)
naiveQuantile(rpkm_counts)

# RPKM
# Unequal gene length
gene.length = c(500, 500, 1000)
rpkm_counts = rpkm(counts, gene.length = gene.length)

rpkm_counts
naiveQuantile(rpkm_counts)

###########################
# EXAMPLE 3
###########################
# Create counts matrix
                  #25 #25 #25 #25
counts <- matrix(c(1, 3, 8, 7,  
                   3, 1, 8, 8,    
                   21, 21, 9, 10),
                 nrow = 3,
                 byrow = TRUE)
colnames(counts) <- c("Treated_1", "Treated_2", "Untreated_1", "Untreated_2")
rownames(counts) <- c("Gene1", "Gene2", "Gene3")

# Create colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Treated", 2), rep("Untreated", 2))
)

# Print counts matrix
print(counts)

# Visualize gene counts
plotGenes(counts)

# Total Counts norm
tc_counts = sweep(counts, 2, colSums(counts), "/")
tc_counts

#RPKM norm
# Same gene length
gene.length = c(500, 500, 500)
rpkm_counts = rpkm(counts, gene.length = gene.length)
rpkm_counts

# Quantile norm
naiveQuantile(counts)
naiveQuantile(tc_counts)
naiveQuantile(rpkm_counts)

# Assuming same gene length, same lib size
naiveTest(counts)
naiveTest(tc_counts)
naiveTest(naiveQuantile(counts))
naiveTest(rpkm_counts)

###########################
# EXAMPLE 4
###########################
# Create counts matrix
#285  #285 #285 #570
counts <- matrix(c(39, 44, 80, 240,
                   25, 36, 60, 180,
                   47, 56, 100, 300,
                   80, 110, 180, 540,
                   264, 236, 50, 150),
                 nrow = 5,
                 byrow = TRUE)
colnames(counts) <- c("Treated_1", "Treated_2", "Untreated_1", "Untreated_2")
rownames(counts) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")

# Create colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Treated", 2), rep("Untreated", 2))
)

# Print counts matrix
print(counts)

# Visualize gene counts
plotGenes(counts)

# Total Counts norm
tc_counts = sweep(counts, 2, colSums(counts), "/")
tc_counts

#RPKM norm
# Same gene length
gene.length = c(500, 500, 500, 500, 500)
rpkm_counts = rpkm(counts, gene.length = gene.length)
rpkm_counts

# Quantile norm
naiveQuantile(counts)
naiveQuantile(tc_counts)
naiveQuantile(rpkm_counts)

# Convert to DESeq2 format
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "Untreated")

dds <- DESeq(dds, fitType = "mean")
res <- results(dds)
res

###########################
# EXAMPLE 5
###########################
# Create counts matrix
#285  #285 #285 #570
counts <- matrix(c(39, 0, 80, 240,
                   0, 36, 60, 180,
                   47, 56, 100, 300,
                   80, 110, 180, 540,
                   264, 236, 50, 150),
                 nrow = 5,
                 byrow = TRUE)
colnames(counts) <- c("Treated_1", "Treated_2", "Untreated_1", "Untreated_2")
rownames(counts) <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")

# Create colData
colData <- data.frame(
  row.names = colnames(counts),
  condition = c(rep("Treated", 2), rep("Untreated", 2))
)

# Print counts matrix
print(counts)

# Visualize gene counts
plotGenes(counts)

# Total Counts norm
tc_counts = sweep(counts, 2, colSums(counts), "/")
tc_counts

#RPKM norm
# Same gene length
gene.length = c(500, 500, 500, 500, 500)
rpkm_counts = rpkm(counts, gene.length = gene.length)
rpkm_counts

# Quantile norm
naiveQuantile(counts)
naiveQuantile(tc_counts)
naiveQuantile(rpkm_counts)

naiveTest(counts)
naiveTest(tc_counts)
naiveTest(naiveQuantile(counts))
naiveTest(rpkm_counts)
naiveTest(naiveQuantile(rpkm_counts))

# Convert to DESeq2 format
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)

dds$condition <- relevel(dds$condition, ref = "Untreated")

dds <- DESeq(dds)
res <- results(dds)
res
