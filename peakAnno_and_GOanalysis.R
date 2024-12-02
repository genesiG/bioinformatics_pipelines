###############
# IMPORT DATA #
###############

# Load required packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(rtracklayer)
library(tidyr)
library(annotate)

# Set your working directory
# Change for each Project
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN")

# Import broadPeak files as GRanges object
CBX2_DMSOpeaks <- import("./csaw/CBX2_DMSO.small.10s25.mid160s50.large500s400.local1k.filter1.0.merge600.q0001.FC2.max50k.bed",
                         format = "BED")

CBX2_DMSOpeaks

CBX2_1uMpeaks <- import("./csaw/CBX2_1uM.small.10s50.mid150s100.large500s400.local1k.filter1.0.merge600.q0001.FC2.max50k.bed",
                         format = "BED")


# Import known genes from hg19 genome
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

annotatedCBX2 <- annotatePeak(CBX2_DMSOpeaks, 
                              tssRegion=c(-3000, 2000),
                              TxDb=txdb, 
                              level = "gene",
                              annoDb="org.Hs.eg.db", 
                              overlap = "all",
                              addFlankGeneInfo = TRUE, 
                              flankDistance = 0) 
annotatedCBX2

annoCBX2_1uM <- annotatePeak(CBX2_1uMpeaks, 
                              tssRegion=c(-3000, 2000),
                              TxDb=txdb, 
                              level = "gene",
                              annoDb="org.Hs.eg.db", 
                              overlap = "all",
                              addFlankGeneInfo = TRUE, 
                              flankDistance = 0) 

# Annotate genes to DIPG H3K27me3 peaks
dipgAnno <- annotatePeak(K27M_vsIgG,
                         tssRegion=c(-3000, 2000),
                         TxDb=txdb, 
                         level = "gene",
                         annoDb="org.Hs.eg.db", 
                         overlap = "all",
                         addFlankGeneInfo = TRUE, 
                         flankDistance = 0) 

# Get data frame with peaks annotated to multiple genes with one flanking gene ID per row
cbx2Anno.df = annotatedCBX2 %>% as.data.frame()

cbx2Anno.genes = cbx2Anno.df[,c("SYMBOL", "score")] #%>% unique()

cbx2Anno.geneId = cbx2Anno.df[,c("geneId", "score")] #%>% unique()


cbx2_1uMAnno.df = annoCBX2_1uM %>% as.data.frame()

cbx2_1uMAnno.genes = cbx2_1uMAnno.df[,c("SYMBOL", "score")] #%>% unique()

cbx2_1uMAnno.geneId = cbx2_1uMAnno.df[,c("geneId", "score")] #%>% unique()


dipgAnno.df = dipgAnno %>% as.data.frame()

dipgAnno.genes = dipgAnno.df[,c("SYMBOL", "score")] #%>% unique()

dipgAnno.geneId = dipgAnno.df[,c("geneId", "score")] #%>% unique()

# Get list containing gene Ids for each condition
genes = list("CBX2 peaks\n(DMSO)" = cbx2Anno.geneId$geneId, 
             "CBX2 peaks\n(EPZ6438)" = cbx2_1uMAnno.geneId$geneId,
             "DIPG peaks\n(DMSO)" = dipgAnno.geneId$geneId,
             "DIPG peaks\n(EPZ6438)" = dipgEPZ.geneId$geneId
             )
names(genes) = sub("_", "\n", names(genes))

# Perform GO enrichment analysis
compGO <- compareCluster(geneCluster   = genes,
                           fun           = "enrichGO",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                         ont = "BP",
                         OrgDb = "org.Hs.eg.db")

# Save dotplot with results
pdf(file = "./images/enrichGO_CBX2vK27me3.pdf", height = 12, width = 10)
dotplot(compGO, 
        showCategory = 10, 
        title = "GO Enrichment Analysis",
        size = "geneRatio") # GeneRatio = genes of interest in the gene set / total genes of interest. 
dev.off()

pdf(file = "./images/enrichGO_DIPG_1uMonly.pdf", height = 8, width = 10)
enrichGO(dipgEPZ.geneId$geneId, pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         ont = "BP",
         OrgDb = "org.Hs.eg.db") %>% dotplot(showCategory = 10)
dev.off()
