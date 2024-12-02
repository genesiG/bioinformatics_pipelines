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
library(dplyr)
library(annotate)

# Set your working directory
# Change for each Project
setwd("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/")

# Import broadPeak files as GRanges object
K27M_vsIgG = rtracklayer::import("./peaks/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.min500.gap2k.broad001.q001.mapQ30.fe30.vsIgG_peaks.broadPeak",
                                 format = "broadPeak")

K27M_4fe = rtracklayer::import("./peaks/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.min500.gap2k.broad001.q001.mapQ30.fe40.vsIgG_peaks.broadPeak",
                                 format = "broadPeak")

K27M_10fe = rtracklayer::import("./peaks/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.min500.gap2k.broad001.q001.mapQ30.fe100.vsIgG_peaks.broadPeak",
                                   format = "broadPeak")

EPZ_10fe = rtracklayer::import("./peaks/RP_054_H3K27me3_500K_1uMcpd_Rpcells_S54.min500.gap2k.broad001.q001.mapQ30.fe100.vsIgG_peaks.broadPeak",
                                format = "broadPeak")

# Import known genes from hg19 genome
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

dipgAnno <- annotatePeak(K27M_vsIgG,
                         tssRegion=c(-3000, 2000),
                         TxDb=txdb,
                         level = "gene",
                         annoDb="org.Hs.eg.db",
                         overlap = "TSS",
                         addFlankGeneInfo = TRUE,
                         flankDistance = 0)

anno.fe4 <- annotatePeak(K27M_4fe, 
                              tssRegion=c(-3000, 2000),
                              TxDb=txdb, 
                              level = "gene",
                              annoDb="org.Hs.eg.db", 
                              overlap = "all",
                              addFlankGeneInfo = TRUE, 
                              flankDistance = 0) 


anno.fe10 <- annotatePeak(K27M_10fe, 
                             tssRegion=c(-3000, 2000),
                             TxDb=txdb, 
                             level = "gene",
                             annoDb="org.Hs.eg.db", 
                             overlap = "all",
                             addFlankGeneInfo = TRUE, 
                             flankDistance = 0) 

annoEPZ.fe10 <- annotatePeak(EPZ_10fe, 
                          tssRegion=c(-3000, 2000),
                          TxDb=txdb, 
                          level = "gene",
                          annoDb="org.Hs.eg.db", 
                          overlap = "all",
                          addFlankGeneInfo = TRUE, 
                          flankDistance = 0) 

# Get data frame with peaks annotated to multiple genes with one flanking gene ID per row
dipgAnno.df = dipgAnno %>% as.data.frame()
dipgAnno.promoters = filter(dipgAnno.df, 
                            annotation == "Promoter (<=1kb)" | 
                              annotation == "Promoter (1-2kb)" |
                              annotation == "Promoter (2-3kb)" 
                            ) 

dipgAnno.genes = dipgAnno.df[,c("SYMBOL", "score")] #%>% unique()
dipgAnno.geneId = dipgAnno.df[,c("geneId", "score")] #%>% unique()

dipgPromoter.genes = dipgAnno.promoters[,c("SYMBOL", "score")] #%>% unique()
dipgPromoter.geneId = dipgAnno.promoters[,c("geneId", "score")] #%>% unique()

fe4Anno.df = anno.fe4 %>% as.data.frame()
fe4Anno.genes = fe4Anno.df[,c("SYMBOL", "score")] #%>% unique()
fe4Anno.geneId = fe4Anno.df[,c("geneId", "score")] #%>% unique()

fe10Anno.df = anno.fe10 %>% as.data.frame()
fe10Anno.genes = fe10Anno.df[,c("SYMBOL", "score")] #%>% unique()
fe10Anno.geneId = fe10Anno.df[,c("geneId", "score")] #%>% unique()

fe10AnnoEPZ.df = annoEPZ.fe10 %>% as.data.frame()
fe10AnnoEPZ.genes = fe10AnnoEPZ.df[,c("SYMBOL", "score")] #%>% unique()
fe10AnnoEPZ.geneId = fe10AnnoEPZ.df[,c("geneId", "score")] #%>% unique()

# Get list containing gene Ids for each condition
genes = list("H3K27me3 (DIPG-VI)" = dipgAnno.geneId$geneId,
             "Mid-enrichment\nPolycomb-binding sites" = fe5Anno.geneId$geneId, 
             "High-enrichment\nPolycomb-binding sites" = fe10Anno.geneId$geneId,
             "EZH2i-resistant\nhigh-enrichment Polycomb sites" = fe10AnnoEPZ.geneId$geneId
             )
names(genes) = sub("_", "\n", names(genes))

# Perform GO enrichment analysis
compGO <- compareCluster(geneCluster   = genes[c(1,3:4)],
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05,
                         pAdjustMethod = "BH",
                         ont = "BP",
                         OrgDb = "org.Hs.eg.db")

# Save dotplot with results
pdf(file = "./images/enrichGO_K27me3_diffEnrich.pdf", height = 12, width = 11)
dotplot(compGO, 
        showCategory = 15, 
        title = "GO Enrichment Analysis",
        size = "geneRatio") # GeneRatio = genes of interest in the gene set / total genes of interest. 
dev.off()

pdf(file = "./images/enrichGO_K27me3_fe3.pdf", height = 8, width = 10)
enrichGO(genes[[1]], pvalueCutoff  = 0.05,
         pAdjustMethod = "BH",
         ont = "BP",
         OrgDb = "org.Hs.eg.db") %>% 
  dotplot(showCategory = 15,
          title = "DIPG-VI Polycomb-binding sites")
dev.off()

#------
# GSEA
#------

# Import H3K27me3 normalized CUT&RUN signal over genomic coordinates
H3K27me3 = import.bedGraph("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/bedgraph/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.TMM.binned.bedgraph")
H3K27me3

k27me3Anno <- annotatePeak(H3K27me3,
                         tssRegion=c(-3000, 2000),
                         TxDb=txdb,
                         level = "gene",
                         annoDb="org.Hs.eg.db",
                         overlap = "all",
                         addFlankGeneInfo = TRUE,
                         flankDistance = 0)

k27me3Anno.df = k27me3Anno %>% as.data.frame()
k27me3Anno.genes = k27me3Anno.df[,c("SYMBOL", "score")] #%>% unique()
k27me3Anno.geneId = k27me3Anno.df[,c("geneId", "score")] #%>% unique()

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)

# we want the some ranked metric, in this case, it is a combination of the adjusted p-value and log2FC 
original_gene_list <- dipgAnno.promoters$signalValue

# name the vector
names(original_gene_list) <- dipgAnno.promoters$SYMBOL 
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
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

# Use the `Gene Set` param for the index in the title, and as the value for geneSetId
gse@result %>% head()
gse$Description[1:20]


gseaplot2(gse, title = gse$Description[3], geneSetID = 3) #+ 
  # geom_vline(xintercept = gse@result$rank[3], 
  #            linetype = "dashed", 
  #            linewidth = 1, 
  #            color = "red") + 
  # geom_vline(xintercept = gse@result$enrichmentScore[3], 
  #            linetype = "dashed", 
  #            linewidth = 1, 
  #            color = "red")
gseaplot2(gse, title = gse$Description[17], geneSetID = 17)
gseaplot2(gse, title = gse$Description[27], geneSetID = 27)
gseaplot2(gse, title = gse$Description[116], geneSetID = 116)
gseaplot2(gse, title = gse$Description[80], geneSetID = 80)
gseaplot2(gse, title = gse$Description[69], geneSetID = 69)

gse@result %>% head()
#------
# HEATMAP
#------

# Load gene expression data set analyzed in DESeq2
resultsDF <- read.csv(file = "C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2/resultsDF.csv", 
                      row.names = 1)

resultsLFC = resultsDF[2] 
resultsLFC %>% head()

# Get matrix of log-transformed reads from RNA-seq data
rlogData <- rlog(dataSet, blind = FALSE)

rlogCounts = assay(rlogData) %>% as.data.frame()
rlogCounts = rlogCounts %>% mutate(label = rownames(rlogCounts))
rlogCounts[1:6] %>% head()

rlogCounts <- rlogCounts %>% 
  as.data.frame()

rlogCounts %>% head()
rlogCounts %>% nrow()

## Filter and scale (get Z-scores) gene expression matrix
rlogGenes <- filter(rlogCounts[1:4], rownames(rlogCounts) %in% dipgAnno@anno$SYMBOL) 
rlogGenes %>% nrow()
rlogGenes <- rlogGenes %>% t() %>% scale() %>% t()



# Import broadPeak files as GRanges object
setwd("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/")
K27M_vsIgG = rtracklayer::import("./peaks/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.min500.gap2k.broad001.q001.mapQ30.fe30.vsIgG_peaks.broadPeak",
                                 format = "broadPeak")

## resize peaks to 1 bp at the center
K27M_vsIgG@ranges = K27M_vsIgG@ranges %>% resize(width = 1, fix = "center")

# Import signal files
H3K27me3 = import.bedGraph("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/bedgraph/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.TMM.binned.bedgraph")
H3K27me3

H3K27me3_1uM = import.bedGraph("C:/Users/thephillipslab/Documents/Projects/H3K27me3_CD34_cells/bedgraph/RP_054_H3K27me3_500K_1uMcpd_Rpcells_S54.TMM.binned.bedgraph")

# Annotate peaks to hg19 genome
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

dipgAnno <- annotatePeak(K27M_vsIgG,
                         tssRegion=c(-3000, 2000),
                         TxDb=txdb,
                         level = "gene",
                         annoDb="org.Hs.eg.db",
                         overlap = "TSS",
                         #addFlankGeneInfo = TRUE,
                         flankDistance = 0
                         )

dipgAnno.df = dipgAnno %>% as.data.frame()
dipgAnno.atGenes = dipgAnno.df %>% filter(annotation != "Distal Intergenic")
dipgAnno.atGenes %>% nrow()
dipgAnno.atGenes$geneId %>% unique() %>% length()

## Sort Data Frame and filter out duplicated genes
dipgAnno.atGenes2 <- dipgAnno.atGenes[order(dipgAnno.atGenes$score, decreasing = T),]
dipgAnno.atGenes2 <- filter(dipgAnno.atGenes2, 
                            !duplicated(dipgAnno.atGenes2$SYMBOL), 
                            #!dipgAnno.atGenes2$SYMBOL %in% rownames(rlogGenes)
                            )[1:nrow(rlogGenes),]

dipgAnno.atGenes2 %>% nrow()


filter(dipgAnno.df,dipgAnno.df$SYMBOL %in% rownames(resultsLFC), !duplicated(dipgAnno.df$SYMBOL)) %>% nrow()
resultsLFC %>% filter(rownames(resultsDF) %in% dipgAnno.df$SYMBOL) %>% nrow()

dipgAnno.filt = filter(dipgAnno.df,dipgAnno.df$SYMBOL %in% rownames(resultsLFC), !duplicated(dipgAnno.df$SYMBOL))
resultsLFC.filt = resultsLFC %>% filter(rownames(resultsDF) %in% dipgAnno.df$SYMBOL)

dipgAnno.filt %>% head()
resultsLFC.filt %>% head()

rownames(resultsLFC.filt) == dipgAnno.filt$SYMBOL

dipgAnno.filt = arrange(dipgAnno.filt, 
                        factor(SYMBOL, levels = rownames(resultsLFC.filt)))


rownames(resultsLFC.filt) == dipgAnno.filt$SYMBOL


library(EnrichedHeatmap)

# Build matrix of signal over K27M peaks
mat1_trim = normalizeToMatrix(signal = H3K27me3, 
                              target = GRanges(dipgAnno.filt), 
                              value_column = "score", 
                              extend = 5000, 
                              include_target = F, 
                              mean_mode = "w0", 
                              w = 100,
                              keep = c(0.10, 0.99),
                              smooth = TRUE
                              )

mat2_trim = normalizeToMatrix(signal = H3K27me3_1uM, 
                              target = GRanges(dipgAnno.filt), 
                              value_column = "score", 
                              extend = 5000, # extend 5kbp upstream and downstream of target (peak)
                              include_target = F, # don't include whole peak region, extend around a single point in the target
                              mean_mode = "w0", 
                              w = 100,
                              keep = c(0.10, 0.99),
                              smooth = TRUE)


library(circlize)
col_fun = colorRamp2(quantile(mat1_trim, c(0, 1.0)), c("white", "navy"))

EnrichedHeatmap(mat1_trim, 
                col = col_fun, 
                name = "H3K27me3 (DMSO)", 
                axis_name = c("-5kb", "peak\ncenter", "+5kb"),
                km = 3,
                use_raster = T
                ) +
  EnrichedHeatmap(mat2_trim, 
                  col = col_fun, 
                  name = "H3K27me3 (1uM EZH2i)", 
                  axis_name = c("-5kb", "peak\ncenter", "+5kb"),
                  use_raster = T) +
  

Heatmap(as.matrix(rlogGenes), 
        col = c("white", "orange"), 
        name = "Log2 FC\nEZH2i vs DMSO",
        cluster_rows = T,
        cluster_row_slices = T,
        clustering_distance_rows = "euclidean",
        clustering_method_rows = "complete",
        show_row_names = FALSE, 
        show_column_names = F,
        cluster_columns = F,
        #row_order = rlogGenes[order(as.data.frame(rlogGenes)$Sorted_S04_VI_LUCB_1uM, decreasing = T),] %>% rownames(),
        #km = 2,
        width = unit(20, "mm"))

# rlogDE = rlogCounts %>% 
#   filter(rownames(rlogCounts) %in% filter(resultsDF, isDE != "Not DE")$row) %>%   
#   as.data.frame()
# 
# rlogDE = rlogDE[1:4] %>% t() %>% scale() %>% t()
# 
# Heatmap(rlogDE, 
#         col = c("blue", "white", "red"), 
#         name = "rLog(counts)",
#         cluster_rows = T,
#         cluster_row_slices = T,
#         clustering_distance_rows = "euclidean",
#         clustering_method_rows = "complete",
#         show_row_names = FALSE, 
#         show_column_dend = F,
#         cluster_columns = F, 
#         km = 2,
#         width = unit(100, "mm"))

clustering = hclust(mat1_trim, method = "complete")
partition = paste0("cluster", cutree(clustering, k = 2)$cluster)
###
mergedMat = cbind(mat1_trim, mat2_trim)
set.seed(123)
partition = paste0(kmeans(mat1_trim, centers = 2)$cluster)
lgd = Legend(at = c("Cluster 1", "Cluster 2"), title = "k-means\nclusters", 
             type = "lines", legend_gp = gpar(col = 2:4), direction = "horizontal")

ht_list = EnrichedHeatmap(mat1_trim, 
                  col = col_fun, 
                  name = "H3K27me3 (DMSO)", 
                  axis_name = c("-5kb", "peak\ncenter", "+5kb"),
                  use_raster = T,
                  #column_labels = "DMSO", 
                  #show_column_names = T, #column_names_rot = 0, column_names_side = "top", 
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),
                                                                           ylim = c(0,11)
                                                                           )
                                                     )
                  ) +
  EnrichedHeatmap(mat2_trim, 
                  col = col_fun, 
                  name = "H3K27me3 (1uM EZH2i)", 
                  axis_name = c("-5kb", "peak\ncenter", "+5kb"),
                  use_raster = T, #heatmap_legend_param = NA,
                  #column_labels = "RNA-seq", 
                  #show_column_names = T, #column_names_rot = 0, column_names_side = "top", 
                  top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),
                                                                           ylim = c(0,11)
                                                                           )
                                                     )
                  ) +
  Heatmap(as.matrix(scale(resultsLFC.filt)), 
          col = c("white", "orange"), 
          name = "Log2 FC\nEZH2i vs DMSO", 
          cluster_columns = F,
          show_row_names = FALSE, 
          column_labels = c("RNA-seq"), 
          #show_column_names = T, #column_names_rot = 0, #column_names_side = "top", 
          show_column_dend = F,
          width = unit(10, "mm"),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:4),
                                                                     outline = FALSE, 
                                                                     axis_param = list(side = "right")))
          )



pdf(file = "./images/H3K27me3_vs_RNAseq.pdf", height = 8, width = 7)
draw(ht_list, split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(8, 8, 8), "mm"))
dev.off()
