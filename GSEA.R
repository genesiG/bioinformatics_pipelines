#
# GENE SET ENRICHMENT ANALYSIS
#

BiocManager::install("clusterProfiler"#, 
                     #version = "3.8"
                     )
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

#### SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

### PREPARE INPUT

# reading in data from deseq2
#df = read.csv("drosphila_example_de.csv", header=TRUE)

resultsDF$metric = -log10(resultsDF$padj) * resultsDF$log2FoldChange/abs(resultsDF$log2FoldChange)
resultsDF %>% head()

# we want the log2 fold change 
original_gene_list <- resultsDF$metric

# name the vector
names(original_gene_list) <- resultsDF$row

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

### GENE SET ENRICHMENT

# ont one of “BP”, “MF”, “CC” or “ALL”
# nPerm the higher the number of permutations you set, the more accurate your 
# result will be, but the longer the analysis will take.
# minGSSize minimum number of genes in set (gene sets with lower than this many 
# genes in your dataset will be ignored).
# maxGSSize maximum number of genes in set (gene sets with greater than this 
# many genes in your dataset will be ignored).
# pvalueCutoff pvalue Cutoff.
# pAdjustMethod one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             nPerm = 100000, 
             minGSSize = 2, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "fdr")

### Dotplot

require(DOSE)

pdf("./heatmaps/GSEA_EZH2i.pdf", 
    #height = 12, width = 9
    )
dotplot(gse, showCategory=10, split=".sign", title = "GO Analysis (EZH2i)") + facet_grid(.~.sign)
dev.off()

#
# KEGG GENE SET ENRICHMENT ANALYSIS
#

#----------
# For KEGG pathway enrichment using the gseKEGG() function, we need to convert 
# id types. We can use the bitr function for this (included in clusterProfiler). 
# It is normal for this call to produce some messages / warnings.
# 
# In the bitr function, the param fromType should be the same as keyType from 
# the gseGO function above (the annotation source). This param is used again in 
# the next two steps: creating dedup_ids and df2.
# 
# toType in the bitr function has to be one of the available options from 
# keyTypes(org.Dm.eg.db) and must map to one of ‘kegg’, ‘ncbi-geneid’, ‘
# ncib-proteinid’ or ‘uniprot’ because gseKEGG() only accepts one of these 
# 4 options as it’s keytype parameter. In the case of org.Dm.eg.db, none of 
# those 4 types are available, but ‘ENTREZID’ are the same as ncbi-geneid for 
# org.Dm.eg.db so we use this for toType.
# 
# As our intial input, we use original_gene_list which we created above.
#----------

### PREPARE INPUT

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully 
# mapped using the bitr function above
df2 = resultsDF[resultsDF$row %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$metric

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


### CREATE gseKEGG OBJECT

# organism KEGG Organism Code: The full list is here: 
# https://www.genome.jp/kegg/catalog/org_list.html (need the 3 letter code). 
# I define this as kegg_organism first, because it is used again below when 
# making the pathview plots.
# nPerm the higher the number of permutations you set, the more accurate your 
# result will be, but the longer the analysis will take.
# minGSSize minimum number of genes in set (gene sets with lower than this many 
# genes in your dataset will be ignored).
# maxGSSize maximum number of genes in set (gene sets with greater than this 
# many genes in your dataset will be ignored).
# pvalueCutoff pvalue Cutoff.
# pAdjustMethod one of “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”.
# keyType one of ‘kegg’, ‘ncbi-geneid’, ‘ncib-proteinid’ or ‘uniprot’.

kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 100000,
               minGSSize    = 1,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "fdr",
               keyType       = "ncbi-geneid")
### Dotplot

pdf("./heatmaps/KEGG.pdf", height = 12, width = 9)
dotplot(kk2, showCategory = 10, title = "Enriched KEGG Pathways (EZH2i)" , 
        split=".sign") + facet_grid(.~.sign)
dev.off()
