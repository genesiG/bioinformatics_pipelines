# CRISPR SCREEN ANALYSIS WITH CaRpools

# Install required packages
packages = c("biomaRt",
             "tidyverse",
             "seqinr",
             "xlsx",
             "rJava",
             "xlsxjars",
             "stringi",
             "scatterplot3d",
             "MESS",
             "DESeq2",
             "rmarkdown",
             "knitr",
             "VennDiagram",
             "sm")

install = vector()
for (each in packages) {
  if (!require(each)){
    install = append(each, install)
  }
}

BiocManager::install(install)

# Prepare library file in the fasta format
### CHANGE FOR EACH PROJECT
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/")
dir.create("./CaRpools")
path_to_library = file.path("./CaRpools/sgRNA_library.fa")
path_to_readcounts = file.path("./CaRpools/EPZ_screen.count.txt")
###


CONTROL1 = load.file("./CaRpools/day1_A_counts.txt", header= TRUE, sep="\t")
CONTROL2 = load.file("./CaRpools/day1_B_counts.txt", header= TRUE, sep="\t")
CONTROL3 = load.file("./CaRpools/day1_C_counts.txt", header= TRUE, sep="\t")
TREAT1 = load.file("./CaRpools/EPZ_d45_A_counts.txt", header= TRUE, sep="\t")
TREAT2 = load.file("./CaRpools/EPZ_d45_B_counts.txt", header= TRUE, sep="\t")
TREAT3 = load.file("./CaRpools/EPZ_d45_C_counts.txt", header= TRUE, sep="\t")
###

CONTROL1 = load.file("./Baseline_A_counts.txt", header= TRUE, sep="\t")
CONTROL2 = load.file("./Baseline_B_counts.txt", header= TRUE, sep="\t")
CONTROL3 = load.file("./Baseline_C_counts.txt", header= TRUE, sep="\t")
TREAT1 = load.file("./GD2_48h_A_counts.txt", header= TRUE, sep="\t")
TREAT2 = load.file("./GD2_48h_B_counts.txt", header= TRUE, sep="\t")

###

# Don't forget the library reference
libFILE = load.file(path_to_library, header = FALSE, type="fastalib")

# Aggregating sgRNA read-count to Gene read-count
CONTROL1.g=aggregatetogenes(data.frame = CONTROL1, agg.function=sum,
                            extractpattern = expression("^(.+?)(_.+)"),type = "aggregate")

# Dataset statistics
U1.stats = stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
                      extractpattern=expression("^(.+?)_.+"), type="stats")
U2.stats = stats.data(dataset=CONTROL2, namecolumn = 1, fullmatchcolumn = 2,
                      extractpattern=expression("^(.+?)_.+"), type="stats")
T1.stats =stats.data(dataset=TREAT1, namecolumn = 1, fullmatchcolumn = 2,
                     extractpattern=expression("^(.+?)_.+"), type="stats")
T2.stats =stats.data(dataset=TREAT2, namecolumn = 1, fullmatchcolumn = 2,
                     extractpattern=expression("^(.+?)_.+"), type="stats")
# Combine Stats
combined.stats = cbind.data.frame(U1.stats[,1:2], U2.stats[,2], T1.stats[,2], T2.stats[,2])
colnames(combined.stats) = c("Readcount", d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2)

combined.stats %>% head()

# DESeq2
data.deseq = stat.DESeq(untreated.list = list(CONTROL1, CONTROL2,CONTROL3),
                        treated.list = list(TREAT1,TREAT2), 
                        namecolumn=1,
                        fullmatchcolumn=2, 
                        extractpattern=expression("^(.+?)(_.+)"),
                        sorting=FALSE, 
                        filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
                        fitType="parametric")

data.deseq
summary(data.deseq$genes$padj)

gene_deseq2 = data.deseq$genes
gene_deseq2_df = gene_deseq2 %>% as.data.frame()
gene_deseq2_df[c("log2FoldChange","padj")] %>% arrange(padj) %>% head(n = 25)

# MAGeCK
data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2), 
                          treated.list = list(TREAT1,TREAT2), 
                          namecolumn=1, 
                          fullmatchcolumn=2, 
                          norm.fun="median", 
                          extractpattern=expression("^(.+?)(_.+)"), 
                          mageckfolder="./",
                          sort.criteria="neg",
                          adjust.method="fdr",
                          filename = "GD2_48h_v_CD19_48h.tab", 
                          fdr.pval = 0.05)

gene_mageck
gene_mageck = data.mageck$genes
gene_mageck_df = gene_mageck %>% as.data.frame()
gene_mageck_df[c("log2FoldChange","padj")] %>% arrange(padj) %>% head(n = 25)