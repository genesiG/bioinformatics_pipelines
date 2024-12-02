# Set the working directory
directory <- "C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2"
setwd(directory)

if(!require(DESeq2)){
  BiocManager::install("DESeq2")
}
if(!require(edgeR)){
  BiocManager::install("edgeR")
}

library("DESeq2")
library("edgeR")

library(edgeR)
library(DESeq2)
library(RColorBrewer)
library(EDASeq)
library(RColorBrewer)

# Define cutoffs
fdr.cutoff=0.05			# 5% FDR cutoff
logfc.cutoff=0			# DE genes must have ? logFC in either direction

### Important! Follow data format in datafile and samplefile
### Important! column names in samplefile should be same(first three) or similar (rest columns)
datafile = "./summary_count.txt"
samplefile = "./sample_list.txt"

count.table <- read.table(datafile, header=T, sep="\t", stringsAsFactors=F, row.name=1, )
#sample.list <- read.table(samplefile, header=T, sep="\t", stringsAsFactors=F, na.strings = "#NA", row.name=1)
sample.list <- read.table(samplefile, header=T, sep="\t", na.strings = "#NA", colClasses = "character")

# Subset the samples of interest
cbx2Counts <- count.table[13:20]
samplesCBX2 <- filter(sample.list, sample.list$sample.name %in% colnames(cbx2Counts))

# Check if the rows in the list of samples are in the same order as the samples in the columns of count.table
colnames(cbx2Counts)
samplesCBX2$sample.name

# If not, sort to make sure the rows in sample.list are in the same order as the columns in count.table
# sample.list[match(colnames(count.table), sample.list)]

attach(sample.list)
groups <- colnames(cbx2Counts) %>% as.factor()
groupA <- factor(sample.list$sample.name[1:12])
groupB <- factor(sample.list$sample.name[13:24])

if(length(levels(groups)) < 2){
  stop("The number of sample groups is too small, check your sample file, i.e. sample_list.txt")
}

colorset <- brewer.pal(max(length(levels(groups)),3), "Paired")

#define your contrasts here
#use values from levels of groups
refs = c()
cmps = c()
for(i in 2:length(groups)){
  t1=unique(groups[i])
  t1=na.omit(t1)
  t2=t1[t1 != ""]
  refs[i-1] <- t2[1]
  cmps[i-1] <- t2[2]
}

########################## DE analysis #################################
# result files : DESeq2.DE.result.xxx.csv & DESeq2.DE.result.xxx.pdf

# Use DESeqDataSetFromMatrix() if you used featureCounts to summarize read counts
# Use DESeqDataSetFromHTSeqCount() if you used HTSeqCounts to summarize read counts

# DESeqDataSetFromMatrix(countData, colData, design, tidy = F)
    # countData: the matrix with your summarized read counts
    # colData: a data.frame whose rows correspond to the columns in countData 
    # design: formula or matrix expressing how the counts for each gene relate to the variables in colData 
    # tidy: whether the first column of countData is the rownames for the count matrix

cbx2DataFrame = data.frame(row.names = colnames(cbx2Counts), 
                    condition = factor(c("DMSO_GuideA", "DMSO_GuideB",
                                  "1uM_GuideA", "1uM_GuideB",
                                  "DMSO_GuideA", "DMSO_GuideB",
                                  "1uM_GuideA", "1uM_GuideB")))

condition = cbx2DataFrame$condition

dds <- DESeqDataSetFromMatrix(countData = cbx2Counts, 
                              colData = cbx2DataFrame, 
                              design = ~condition)
dds <- DESeq(dds)

plotMA(dds, main = "CBX2 (DMSO vs 1uM EZH2i)")

cbx2DE_A <- results(dds, contrast = c("condition", "DMSO_GuideA", "1uM_GuideA"))
plotMA(cbx2DE_A, main = "CBX2 (DMSO vs 1uM EZH2i)")


#pdf("plotMDS_original.pdf")
#plotMDS(cbx2Counts, col=colorset[groups], main="original data")
#dev.off()

#pdf("plotMDS_normalized.pdf")
#plotMDS(counts(dds,normalized=TRUE), col=colorset[groups], main="DESeq normalized data")
#dev.off()

for(i in 1:length(refs)){		    
  ref=refs[i]
  cmp=cmps[i]
  res <- results(dds,contrast=c("groups",cmp,ref))
  
  pdf(paste("DESeq2.DE.result", paste(cmp, "vs", ref, sep="_"), "pdf",sep="."))
  plotMA(res, alpha=0.05, ylim=range(res$log2FoldChange[!is.na(res$log2FoldChange)]))
  hist(res[,"pvalue"], breaks=40, probability=T, plot=T, main="Distribution of P-Values", xlab="Gene P-values", ylab="Density", col="lightgrey")
  hist(res[,"padj"], breaks=40, probability=T, plot=T, main="Distribution of adjusted P-Values", xlab="Gene P-values", ylab="Density", col="lightgrey")
  dev.off()  
  
  res$DE.Call <- "Same"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] >= logfc.cutoff) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Up"
  res[(res[,"padj"] <= fdr.cutoff) & (res[,"log2FoldChange"] <= (-1 * logfc.cutoff)) & !is.na(res$padj) & !is.na(res$log2FoldChange),"DE.Call"] <- "Down"
  
  c <- DGEList(counts=count.table, group=groups)
  c$counts <- counts(dds, normalized=TRUE)
  c <- cpm(c,log=T)
  c <- c[,groups == cmp | groups == ref]
  c <- apply(c,1,mean)
  res$logCPM <- c
  
  DESeq2file = paste("DESeq2.DE.result", paste(cmp, "vs", ref, sep="_"), "csv",sep=".")	
  write.csv(res, file=DESeq2file)
}		
