#==============================================================================
# CRISPR-SCREEN ANALYSIS WITH edgeR
#==============================================================================

if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if(!require("edgeR")) {
  BiocManager::install("edgeR")
} 
if(!require("pheatmap")) {
  BiocManager::install("pheatmap")
} 
#if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(edgeR)
#library(clusterProfiler)
#library(fgsea)
library(ggplot2)
library(dplyr)

#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_count_table = file.path("./EZH2i_screen_counts_filtered_final_1_in_10_cpm_cutoff.rds")
#===

#===
# LOAD COUNTS
#===

# Load DGE List containing reads extrated with processAmplicon()
counts <- readRDS(path_to_count_table)

counts %>% head()
counts <- counts[, grepl("RP", colnames(counts$counts))]
# Get all the column names containing the specified patterns
# selected_cols <- colnames(counts$counts)[grepl("conc1-day45", colnames(counts$counts)) |
#                                            grepl("DMSO-day45", colnames(counts$counts)) |
#                                            grepl("day1-", colnames(counts$counts))]

# Subset the 'counts' object with the selected columns
# counts <- counts[, selected_cols]

counts %>% head()

library(pheatmap)
corr_coeff = cor(cpm(counts$counts), method = "spearman")
as.dist(1-corr_coeff, upper = TRUE) %>%
  as.matrix %>%
  pheatmap::pheatmap(., main = "Spearman correlation")


# Resume exploratory analysis
subset_counts = counts

selc = "RP-conc1-day45-replicate_A"
subset_counts$counts <- subset_counts$counts[, !colnames(subset_counts$counts) == selc]
subset_counts$samples <- subset_counts$samples[subset_counts$samples["ID"] != selc, ]

subset_counts$samples$Grouped = paste(subset_counts$samples$Treatment, subset_counts$samples$Time, sep = "_")


# Set up colors by cell type
leg = subset_counts$samples$Grouped %>% factor() %>% unique()
cols = subset_counts$samples$Grouped %>% factor() %>% as.numeric() 

labels <- ifelse(grepl("_A", colnames(subset_counts$counts)), "A",
                 ifelse(grepl("_B", colnames(subset_counts$counts)), "B",
                        ifelse(grepl("_C", colnames(subset_counts$counts)), "C", "")))

# Setting up pdf file name to export
mds_name = paste(c(cell_type, "vs", cell_type2), collapse = "_")
mds_name = paste0(c(mds_name, "pdf"), collapse = ".")

# Export MDS plots to visualise relationships between replicate samples 
pdf(paste0("./plots/", mds_name), width = 5, height = 5)

plotMDS(subset_counts,
        labels = labels, 
        col = cols, # colors
        gene.selection = "common"
        )
legend("bottomright",fill = leg, 
       legend = leg, # unique() will only show unique legend labels 
       col = unique(cols), # one color per unique legend label
       )

dev.off()

# # Assuming 'counts' is your count matrix
# subset_counts <- counts
# 
# # Create the 'Grouped' column in 'samples'
# subset_counts$samples$Grouped <- paste(subset_counts$samples$Treatment, subset_counts$samples$Time, sep = "_")
# 
# # Set up colors by 'Grouped' column
# leg <- subset_counts$samples$Grouped %>% factor() %>% unique()
# cols <- subset_counts$samples$Grouped %>% factor() %>% as.numeric()
# 
# # Create 'labels' based on the column names in 'counts'
# labels <- ifelse(grepl(" A", colnames(subset_counts$counts)), "A",
#                  ifelse(grepl(" B", colnames(subset_counts$counts)), "B",
#                         ifelse(grepl(" C", colnames(subset_counts$counts)), "C", "")))
# 
# # Set up the pdf file name to export
# mds_name <- paste(c(cell_type, "vs", cell_type2), collapse = "_")
# mds_name <- paste0(c(mds_name, "pdf"), collapse = ".")
# 
# # Plot MDS
# plotMDS(subset_counts,
#         labels = labels, 
#         col = cols, # colors
#         gene.selection = "common"
# )

# Check sample labels to be used in the next steps
table(counts$samples$Treatment)


c = cpm(subset_counts$counts) %>% as.data.frame()
filter(c, grepl("SUV39H2", rownames(c)))
filter(c, grepl("CBX2", rownames(c)))
filter(c, grepl("CBX4", rownames(c)))

#===
# DEFINE WHICH COMPARISONS TO MAKE
#===
subset_counts = counts
selc = "RP-conc1-day45-replicate_A"
subset_counts$counts <- subset_counts$counts[, !colnames(subset_counts$counts) == selc]
subset_counts$samples <- subset_counts$samples[subset_counts$samples["ID"] != selc, ]

treat = "RP-conc1-day45"
control = "RP-day1"
subset_counts$counts <- subset_counts$counts[, c(grep(treat, colnames(subset_counts$counts)),
                            grep(control, colnames(subset_counts$counts)))]
subset_counts$samples <- subset_counts$samples[c(grep(treat, subset_counts$samples$ID),
                                          grep(control, subset_counts$samples$ID)),]

# Check the subsetted samples table
subset_counts$samples

table(subset_counts$samples$Treatment)

#===
# ACTUAL STATISTICAL TESTING
#===

# Begin differential representation analysis
# We will use GLMs to access additional functionaly within edgeR, such fold change testing with TREAT 

# Set up design matrix for GLM
treatment = as.factor(subset_counts$samples$Treatment)
treatment = relevel(treatment, ref = "Baseline")
treatment
#timepoints = as.factor(subset_counts$samples$Time)
des = model.matrix(~treatment)
des

# Estimate dispersions
xglm = estimateDisp(subset_counts, des)

# Plot BCVs versus abundance
sqrt(xglm$common.disp)
plotBCV(xglm, main = "Treat vs Control: BCV Plot")

# Fit negative bionomial GLM
fit = glmFit(xglm, des)

# Carry out Likelihood ratio test
lrt = glmLRT(fit, coef = 2)

# Show top ranked sgRNAs
topTags(lrt, n = 30)

topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "CBX2")
topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "CBX4")
topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "SUV39H2")
topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "Non-Targeting Control")
topTags(lrt, n = Inf) %>% as.data.frame() %>% nrow()

filter(c, grepl("SUV39H2", rownames(c)))
filter(c, grepl("CBX4", rownames(c)))

# Export results table
file_name = paste(c(treat, "vs", control, "1_in_10_cpm"), collapse = "_")
file_name = paste0("./edgeR/", file_name)
topTags(lrt, n = Inf) %>% write.table(file = paste(c(file_name, "txt"), 
                                                   collapse = ".", 
                                                   sep = "\t"), 
                                      row.names = T, 
                                      col.names = T)



# Rank-ordered plot
library(ggplot2)
library(ggrepel)

rankedTags = topTags(lrt, n = Inf) %>% 
  as.data.frame() %>% 
  arrange(logFC)

a = topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "Non-Targeting Control")
mean(a$logFC)
median(a$logFC)
rankedTags$Score = rankedTags$logFC - median(a$logFC)

#rankedTags %>% head()
filter(rankedTags, Gene == "CBX2")
filter(rankedTags, Gene == "CBX4")
filter(rankedTags, Gene == "SUV39H2")
# rankedTags$rank = 1:nrow(rankedTags)
# 
# # Get the indexes for the first 10 unique values in the Gene column
# unique_genes <- unique(rankedTags$Gene)
# depleted_indexes <- match(unique_genes[1:10], rankedTags$Gene)
# enriched_indexes <- match(unique_genes[(length(unique_genes)-10):length(unique_genes)], rankedTags$Gene)
# 
# # Label top 10 depleted
# rankedTags$label = NA
# rankedTags$label[depleted_indexes] = rankedTags$Gene[depleted_indexes]
# 
# # # Label top 10 enriched
# # rankedTags$label[enriched_indexes] = rankedTags$Gene[enriched_indexes]
# 
# # Color top 10 depleted
# rankedTags$category = NA
# rankedTags$category[depleted_indexes] = "depleted"
# 
# # # Color top 10 enriched
# # rankedTags$category[enriched_indexes] = "enriched"
# 
# rankedTags %>% filter(Gene == "CBX4")
# rankedTags %>% filter(Gene == "CBX2")
# rankedTags %>% filter(Gene == "SUV39H2" | Gene == "KMT1B")
# 
# # Get table of average logFC by gene
averages <- rankedTags %>%
  group_by(Gene) %>%
  slice(1:3) %>%
  summarize(Avg_Score = mean(Score)) %>% arrange(Avg_Score)

#averages <- aggregate(Score ~ Gene, data = rankedTags, FUN = mean) %>% arrange(Score)

averages$rank = 1:nrow(averages)

averages$label = NA
genes_of_interest = c("BMI1", "CBX2", "CDK9", "RPA3", "CBX4")
averages$label[averages$Gene %in% genes_of_interest] = averages$Gene[averages$Gene %in% genes_of_interest]

averages %>% head(n = 30)


filter(averages, Gene == "SUV39H2")
filter(averages, Gene == "CBX4")
filter(averages, Gene == "CBX2")

p = ggplot(data = averages, aes(x = rank, 
                              y = Avg_Score 
                              )) +
  geom_point(color = ifelse(is.na(averages$label), "black", "red"),
             size = ifelse(is.na(averages$label), 1.5, 2.0)) +
  xlab("Gene Rank") +
  ylab("Average Score") +
  geom_text_repel(aes(label = label),
                  color = "red", 
                  max.overlaps = 8, 
                  size = 4.5,
                  xlim = c(60,500), 
                  ylim = c(-8,-2),
                  force_pull = 0.5,
                  #nudge_x = 10, 
                  nudge_y = -1.5
  ) +
  ggtitle("EZH2i (day 45) vs Baseline (day 1)", subtitle = "edgeR") +
  theme_bw() + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, size = FALSE)

p


file_name = paste(c(treat, "vs", control, "rankedPlot"), collapse = "_")
plot_name = paste0("./plots/", file_name, ".pdf")




pdf(plot_name, width = 5, height = 5)
p
dev.off()



























#==============================================================================
# CRISPR-SCREEN ANALYSIS WITH edgeR
#==============================================================================

if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if(!require("edgeR")) {
  BiocManager::install("edgeR")
} 
if(!require("pheatmap")) {
  BiocManager::install("pheatmap")
} 
#if(!"clusterProfiler" %in% installed.packages()) BiocManager::install("clusterProfiler")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")

library(edgeR)
#library(clusterProfiler)
#library(fgsea)
library(ggplot2)
library(dplyr)

#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
# Count table when allowing mismatches
path_to_count_table = file.path("./counts/EZH2i_screen_counts_allowMismatches.rds")
#===

#===
# LOAD COUNTS
#===

# Load DGE List containing reads extrated with processAmplicon()
counts <- readRDS(path_to_count_table)

counts %>% head()
counts %>% nrow()
counts$samples %>% head()
#===
# DEFINE WHICH COMPARISONS TO MAKE
#===
subset_counts = counts

treat = "RP-DMSO-day45"
control = "RP-day1"
subset_counts$counts <- subset_counts$counts[, c(grep(treat, colnames(subset_counts$counts)),
                                                 grep(control, colnames(subset_counts$counts)))]
subset_counts$samples <- subset_counts$samples[c(grep(treat, subset_counts$samples$ID),
                                                 grep(control, subset_counts$samples$ID)),]

# Check the subsetted samples table
subset_counts$samples

table(subset_counts$samples$Treatment)

counts = subset_counts

#===
# FILTERING
#===
# Total read depth in the dataset BEFORE filtering
counts$samples$lib.size %>% sum()

##############################
# Filtering guides by low counts
keep <- filterByExpr(counts, group = counts$samples$Treatment)
#keep <- filterByExpr(counts, group = counts$samples$Treatment, min.count = 0, min.total.count = 3)
counts <- counts[keep, , keep.lib.sizes=FALSE] # update library sizes after filtering
# Total read depth in the dataset after filtering low counts
counts$samples$lib.size %>% sum()

# Define which samples to exclude
selc = c("DMSO-day45-replicate_C",
         "day1-replicate_A")
# Create a combined pattern using the | (OR) operator
combined_pattern <- paste(selc, collapse = "|")
# Use grepl with the combined pattern to exclude columns
counts$counts <- counts$counts[, !grepl(combined_pattern, colnames(counts$counts))]
counts$samples <- counts$samples[!grepl(combined_pattern, counts$samples$ID), ]
# Total read depth in the dataset AFTER filtering outlier samples
counts$samples$lib.size %>% sum()


# Formating for MAGeCK

# Export


##############################


# keep <- filterByExpr(counts, group = counts$samples$Treatment)
# counts <- counts[keep, , keep.lib.sizes=FALSE]
# counts %>% head()

# Checking correlation after filtering
library(pheatmap)
corr_coeff = cor(cpm(counts$counts), method = "pearson")
as.dist(1-corr_coeff, upper = TRUE) %>%
  as.matrix %>%
  pheatmap::pheatmap(., main = "Pearson correlation")

# Check how many guides were left for some genes of interest
# c = cpm(counts$counts) %>% as.data.frame()
# filter(c, grepl("SUV39H2", rownames(c)))
# filter(c, grepl("CBX2", rownames(c)))
# filter(c, grepl("CBX4", rownames(c)))
# #filter(c, grepl("Control", rownames(c)))
# 
# 
# #===
# # MDS PLOT
# #===
# subset_counts = counts
# 
# subset_counts$samples$Grouped = paste(subset_counts$samples$Treatment, subset_counts$samples$Time, sep = "_")
# 
# # Set up colors by cell type
# leg = subset_counts$samples$Grouped %>% factor() %>% unique()
# cols = subset_counts$samples$Grouped %>% factor() %>% as.numeric() 
# 
# labels <- ifelse(grepl("_A", colnames(subset_counts$counts)), "A",
#                  ifelse(grepl("_B", colnames(subset_counts$counts)), "B",
#                         ifelse(grepl("_C", colnames(subset_counts$counts)), "C", "")))
# 
# # Setting up pdf file name to export
# mds_name = paste(c(cell_type, "vs", cell_type2), collapse = "_")
# mds_name = paste0(c(mds_name, "pdf"), collapse = ".")
# 
# # Export MDS plots to visualise relationships between replicate samples 
# pdf(paste0("./plots/", mds_name), width = 5, height = 5)
# 
# plotMDS(subset_counts,
#         labels = labels, 
#         col = cols, # colors
#         gene.selection = "common"
# )
# legend("bottomright",fill = leg, 
#        legend = leg, # unique() will only show unique legend labels 
#        col = unique(cols), # one color per unique legend label
# )
# 
# dev.off()


#############################
subset_counts = counts
# Create column to store replicate information
# Replicate label should be the last character of the sample name (i.e. "_A") for this to work
subset_counts$samples$replicates <- substr(subset_counts$samples$ID,
                                           nchar(subset_counts$samples$ID),
                                           nchar(subset_counts$samples$ID)
)
# Create column to group samples by timepoint and replicate information
subset_counts$samples$grouping = paste(subset_counts$samples$Time,
                                       subset_counts$samples$replicates, sep = "_") %>% factor()
# Check table
subset_counts$samples %>% head()

# Set up colors by cell type
leg = subset_counts$samples$Treatment %>% factor() %>% unique()
cols = subset_counts$samples$Treatment %>% factor() %>% as.numeric()
# Set up labels for the plot
labels = subset_counts$samples$grouping %>% factor() #%>% as.numeric() + 12
leg2 = subset_counts$samples$Time %>% factor() %>% unique()
# Plot
plotMDS(subset_counts,
        labels = labels,
        col = cols, # colors
        gene.selection = "pairwise"
)
legend("topright",fill = leg,
       legend = leg, # unique() will only show unique legend labels
       col = unique(cols), # one color per unique legend label
)
#############################


#===
# ACTUAL STATISTICAL TESTING
#===

# Begin differential representation analysis
# We will use GLMs to access additional functionality within edgeR, such fold change testing with TREAT 

# Set up design matrix for GLM
treatment = as.factor(subset_counts$samples$Treatment)
treatment = relevel(treatment, ref = "Baseline")
treatment
#timepoints = as.factor(subset_counts$samples$Time)
des = model.matrix(~treatment)
des

# Estimate dispersions
xglm = estimateDisp(subset_counts, des)

# Plot BCVs versus abundance
sqrt(xglm$common.disp)
plotBCV(xglm, main = "Treat vs Control: BCV Plot")

# Fit negative bionomial GLM
fit = glmFit(xglm, des)

# Carry out Likelihood ratio test
lrt = glmLRT(fit, coef = 2)

# Show top 30 ranked sgRNAs
#topTags(lrt, n = 30)

topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "CBX2")
topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "CBX4")
topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "SUV39H2")
#topTags(lrt, n = Inf) %>% as.data.frame() %>% filter(Gene == "Non-Targeting Control") %>% head()
topTags(lrt, n = Inf) %>% as.data.frame() %>% nrow()

# Compute scores as Z-scored logFoldChanges
##########################
a = topTags(lrt, n = Inf) %>% as.data.frame() 
a$score = scale(a$logFC)
a %>% head()
a$score %>% summary()
a$logFC %>% summary()
filter(a, Gene == "CBX2")
filter(a, Gene == "CBX2")$logFC %>% mean() # average logFC
filter(a, Gene == "CBX2")$score %>% mean() # average Z-score

filter(a, Gene == "CBX4")
filter(a, Gene == "CBX4")$logFC %>% mean() # average logFC
filter(a, Gene == "CBX4")$score %>% mean() # average Z-score
##########################

### Compute CRISPR scores as how many SDs away from the the average LFC of NTC
##########################
# a = topTags(lrt, n = Inf) %>% as.data.frame() 
meanNTC = filter(a, grepl("Control", a$Gene))$logFC %>% mean()
meanNTC
sdNTC = filter(a, grepl("Control", a$Gene))$logFC %>% sd()
sdNTC

a$CRISPR_score = (a$logFC-meanNTC)/sdNTC
a %>% head()
a$CRISPR_score %>% summary()
a$logFC %>% summary()
filter(a, Gene == "CBX2")
filter(a, Gene == "CBX2")$CRISPR_score %>% mean() # average CRISPR-score

filter(a, Gene == "CBX4")
filter(a, Gene == "CBX4")$CRISPR_score %>% mean() # average CRISPR-score
#########################

### Compute CRISPR scores as LFC - median(LFC) of NTC
##########################
# a = topTags(lrt, n = Inf) %>% as.data.frame() 
medianNTC = filter(a, grepl("Control", a$Gene))$logFC %>% median()
medianNTC

a$score_medianNTC = (a$logFC-medianNTC)
a %>% head()
a$score_medianNTC %>% summary()
a$logFC %>% summary()
filter(a, Gene == "CBX2")
filter(a, Gene == "CBX2")$score_medianNTC %>% mean() # average CRISPR-score

filter(a, Gene == "CBX4")
filter(a, Gene == "CBX4")$score_medianNTC %>% mean() # average CRISPR-score
#########################



# Export results table
file_name = paste(c("filterByExpr_filterDMSODay45C_filterday1A", "mismatchesAllowed", treat, "vs", control), collapse = "_")
file_name = paste0("./edgeR/filter_vs_noFilter/", file_name)
a %>% write.table(file = paste(c(file_name, "txt"), 
                               collapse = ".", 
                               sep = "\t"), 
                  row.names = T, 
                  col.names = T)


#===
# EPIFACTORS CLASSIFICATION
#===
path_to_library = file.path("./humanEpiLibrary.xlsx")

# Load table of genes in the epifactor classification
epifactors = read.delim("./genes_epifactorClassiication.csv")
epifactors %>% head()
epifactors %>% nrow()

# Load file with sgRNA target sequences
if (grepl("\\.xlsx$", path_to_library)) {
  library(readxl)
  sgRNA_library <- read_excel(path_to_library)
} else if (grepl("\\.csv$", path_to_library)) {
  sgRNA_library <- read.csv(path_to_library)
} else if (grepl("\\.txt$", path_to_library)) {
  sgRNA_library <- read.table(path_to_library)
} else {
  stop("Unsupported file format")
} 

sgRNA_library$`Target Gene Symbol` %>% unique() %>% length()

screenGenes = filter(epifactors, HGNC.approved.symbol %in% sgRNA_library$`Target Gene Symbol`)
screenGenes %>% head()
screenGenes %>% nrow()

screenGenes %>% write.table(file = "./screenGenes_epifcator.txt",                   
                            row.names = T, 
                            col.names = T)

#===
# VOLCANO PLOT  
#===
results = read.table("./edgeR/RP-conc1-day45_vs_RP-day1_filterByExpr.txt")
results %>% head()

screenGenes = read.table("./screenGenes_epifcator.txt")
screenGenes %>% head()


results_merge = merge(results, screenGenes[,c("HGNC.approved.symbol", "Protein.complex")], 
                      by.x = "Gene", 
                      by.y = "HGNC.approved.symbol",
                      all.x = T)
results_merge %>% head()

results_merge$colors = NA
results_merge$colors[grepl("PRC1",results_merge$Protein.complex)] <- "PRC1"
results_merge$colors[grepl("PRC2",results_merge$Protein.complex)] <- "PRC2"
results_merge$colors[results_merge$FDR > 0.05] <- NA
results_merge$colors[results_merge$logFC > -1] <- NA

filter(results_merge, grepl("CBX", results_merge$Gene))

# Make volcano plot
comparison = "EZH2i vs baseline"

volcanoPlot <- ggplot(data=results_merge, 
                      aes(x=logFC, 
                          y=-log10(FDR),
                          color = colors, 
                          #shape = Peak, 
                          #label =  Label
                          )
                      ) + 
  geom_point() + 
  geom_point(data = filter(results_merge, !is.na(colors)), size = 2.5) + 
  theme_bw() +
  #scale_shape_manual(values = c(1,16,4,17)) +
  ylab(expression(bold(bolditalic(-Log[10]) ~ "FDR"))) +
  xlab(expression(bold(bolditalic(Log[2]) ~ "FC"))) +
  ggtitle(comparison) + 
  scale_x_continuous(n.breaks = 7, 
                     limits = c(min(results_merge$logFC),min(results_merge$logFC)*-1)
  ) +
  theme(legend.title = element_blank(),  
        #axis.title.x = element_text(face = "bold"), 
        #axis.title.y = element_text(face = "bold"), 
        panel.grid.minor = element_blank()) 

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

pdf_name = paste(c(treat, "vs", control, "volcanoPlot"), collapse = "_")
pdf_name = paste0("./plots/", pdf_name, ".pdf")
pdf(pdf_name, width = 5.75, height = 5)

volcanoPlot +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-1, 1), col="black", linetype = 2) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 2) +
  # Change point color
  scale_color_manual(values=c("PRC1" = "purple",
                              "PRC2" = "orange"))
dev.off()

#===
# LFC VS LFC PLOT
#===
results1 = read.table("./edgeR/RP-conc1-day45_vs_RP-day1_filterByExpr.txt")
results2 = read.table("./edgeR/RP-DMSO-day45_vs_RP-day1_filterByExpr.txt")

screenGenes = read.table("./screenGenes_epifcator.txt")
screenGenes %>% head()

library(dplyr)

merged <- merge(select(results1, -Sequences, -logCPM, -LR, -PValue),
                select(results2, -Sequences, -logCPM, -LR, -PValue),
                by = c("ID", "Gene"),
                suffixes = c(".EZH2i", ".DMSO"))
results_merge = merge(merged, screenGenes[,c("HGNC.approved.symbol", "Protein.complex")], 
                      by.x = "Gene", 
                      by.y = "HGNC.approved.symbol",
                      all.x = T)
results_merge %>% head()

results_merge$colors = "Other complexes"
results_merge$colors[grepl("PRC",results_merge$Protein.complex)] <- "PRC"

# Calculate the difference between the y-value and the reference line (slope = 1, intercept = 1)
# results_merge$above_reference_line <- ifelse(-results_merge$logFC.DMSO + -results_merge$logFC.EZH2i > 1, "Above", "Below")
results_merge$above_reference_line <- -results_merge$logFC.DMSO + (-results_merge$logFC.EZH2i)
# results_merge <- mutate(results_merge,
                        # above_reference_line = factor(above_reference_line))
# Get table of depleted guides
results_merge_dep = filter(results_merge, logFC.EZH2i < 0)

limits = c(min(range(results_merge_dep$logFC.EZH2i)), max(range(range(results_merge_dep$logFC.EZH2i)))) %>% rev()


p = ggplot(data = results_merge_dep, 
       aes(x = -logFC.DMSO, 
           y = -logFC.EZH2i,
           #color = colors
       )
) + 
  geom_point(data = results_merge_dep %>% filter(-logFC.DMSO - (-logFC.EZH2i) <= -1), 
             size = 2, 
             shape=21, 
             aes(fill = colors)) +
  geom_point(data = results_merge_dep %>% filter(-logFC.DMSO - (-logFC.EZH2i) <= -1 & colors == "PRC"), 
             size = 2.5, 
             shape=21, 
             # fill = "red", 
             aes(fill = colors, color = colors)) +
  geom_point(data = results_merge_dep %>% filter(-logFC.DMSO - (-logFC.EZH2i) > -1), 
             size = 2, 
             color = "gray",
             alpha = 0.6) +
  theme_bw() +
  ylab(expression(bold(bolditalic(-Log[2]) ~ "FC / EZH2i vs baseline"))) +
  xlab(expression(bold(bolditalic(-Log[2]) ~ "FC / DMSO vs baseline"))) +
  scale_x_continuous(limits = -limits) +
  scale_y_continuous(limits = -limits) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 13),
        panel.grid = element_blank()) +
  # Add reference lines
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "blue") +
  geom_abline(intercept = -1, slope = 1, linetype = "dashed", color = "blue") +
  # Change point color
  scale_color_manual(values = c("PRC" = "black", "Other genes" = "black")) +
  scale_fill_manual(values = c("PRC" = "purple", "Other complexes" = "black"),
                    guide = guide_legend(title = "Gene Group")) +
  guides(color = FALSE)


p


pdf("./plots/LFC_vs_LFC_EZH2i_DMSO.pdf", width = 6.65, height = 5)
p
dev.off()

#===
# RANK-ORDERED PLOT
#===

# Rank-ordered plot
library(ggplot2)
library(ggrepel)

results1 = read.table("./edgeR/RP-conc1-day45_vs_RP-day1_filterByExpr.txt")
results2 = read.table("./edgeR/RP-DMSO-day45_vs_RP-day1_filterByExpr.txt")

rankedTags = results2 %>% 
  arrange(logFC)

filter(rankedTags, Gene == "CBX2")
filter(rankedTags, Gene == "CBX4")
filter(rankedTags, Gene == "SUV39H2")
filter(rankedTags, grepl("Control", rankedTags$Gene))

averages <- rankedTags %>%
  group_by(Gene) %>%
  # slice(1:3) %>%
  summarize(Avg_Score = mean(score)) %>% arrange(Avg_Score)

#averages <- aggregate(Score ~ Gene, data = rankedTags, FUN = mean) %>% arrange(Score)
averages %>% head()
averages$rank = 1:nrow(averages)

averages$label = NA
genes_of_interest = c("BMI1", "CBX2", "CDK9", "RPA3", "CBX4")
averages$label[averages$Gene %in% genes_of_interest] = averages$Gene[averages$Gene %in% genes_of_interest]

averages %>% head(n = 30)


filter(averages, Gene == "SUV39H2")
filter(averages, Gene == "CBX4")
filter(averages, Gene == "CBX2")

p = ggplot(data = averages, aes(x = rank, 
                                y = Avg_Score 
)) +
  geom_point(color = ifelse(is.na(averages$label), "black", "red"),
             size = ifelse(is.na(averages$label), 1.5, 2.0)) +
  xlab("Gene Rank") +
  ylab("Average Score") +
  geom_text_repel(aes(label = label),
                  color = "red", 
                  max.overlaps = 8, 
                  size = 4.5,
                  xlim = c(60,500), 
                  ylim = c(-8,-2),
                  force_pull = 0.5,
                  #nudge_x = 10, 
                  nudge_y = -1.5
  ) +
  ggtitle("DMSO vs baseline") +
  theme_bw() + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, size = FALSE)

p


file_name = paste(c(treat, "vs", control, "rankedPlot"), collapse = "_")
plot_name = paste0("./plots/", file_name, ".pdf")




pdf(plot_name, width = 5, height = 5)
p
dev.off()

#===
# COMMON HITS
#===
results1 = read.table("./edgeR/RP-conc1-day45_vs_RP-day1_filterByExpr.txt")
results2 = read.table("./edgeR/RP-DMSO-day45_vs_RP-day1_filterByExpr.txt")

tab1_sig = filter(results1, FDR <= 0.05)
tab2_sig = filter(results2, FDR <= 0.05)

tab1_sig_dep = filter(tab1_sig, logFC <= -1)
tab2_sig_dep = filter(tab2_sig, logFC <= -1)

results1 %>% nrow()
tab1_sig %>% nrow()
tab1_sig_dep %>% nrow()
results2 %>% nrow()
tab2_sig %>% nrow()
tab2_sig_dep %>% nrow()


commonHits = filter(tab1_sig_dep, ID %in% tab2_sig_dep$ID)
commonHits %>% nrow()

dmso_only = filter(tab2_sig_dep, !ID %in% tab1_sig_dep$ID)
nrow(dmso_only) + nrow(commonHits) == nrow(tab2_sig_dep)

treat_only = filter(tab1_sig_dep, !ID %in% tab2_sig_dep$ID)
nrow(treat_only) + nrow(commonHits) == nrow(tab1_sig_dep)

nrow(dmso_only)
nrow(commonHits)
nrow(treat_only)

#===
# OVERREPRESENTATION ANALYSIS
#===
results = read.table("./edgeR/RP-conc1-day45_vs_RP-day1_filterByExpr.txt")
results %>% head()

screenGenes = read.table("./screenGenes_epifcator.txt")
screenGenes %>% head()


results_merge = merge(results, screenGenes[,c("HGNC.approved.symbol", "Protein.complex")], 
                      by.x = "Gene", 
                      by.y = "HGNC.approved.symbol",
                      all.x = T)
results_merge %>% head()

# Create the 2x2 contingency table
dep_PRC = filter(results_merge, logFC <= -1 & FDR <= 0.05 & grepl("PRC", results_merge$Protein.complex)) %>% nrow()
dep_notPRC = filter(results_merge, logFC <= -1 & FDR <= 0.05 ) %>% nrow() - dep_PRC
depleted_genes <- c(dep_PRC, dep_notPRC)
notDep_PRC = filter(results_merge, (logFC > -1 | FDR > 0.05) & grepl("PRC", results_merge$Protein.complex)) %>% nrow()
notDep_notPRC = filter(results_merge, (logFC > -1 | FDR > 0.05) ) %>% nrow() - notDep_PRC
non_depleted_genes <- c(notDep_PRC, notDep_notPRC)
contingency_table <- rbind(depleted_genes, non_depleted_genes)

contingency_table

# Perform Fisher's exact test
fisher_result <- fisher.test(contingency_table)

# Print the test result
print(fisher_result)

library(ggplot2)

# Create a data frame for the stacked bar plot
bar_data <- data.frame(
  Group = c("Depleted PRC Genes", "Depleted Non-PRC Genes", "Non-Depleted PRC Genes", "Non-Depleted Non-PRC Genes"),
  Count = c(dep_PRC, dep_notPRC, notDep_PRC, notDep_notPRC)
)

# Calculate the total number of depleted and non-depleted genes
total_depleted <- sum(bar_data$Count[1:2])
total_non_depleted <- sum(bar_data$Count[3:4])

# Calculate the proportion of PRC genes among depleted and non-depleted genes
bar_data$Proportion <- bar_data$Count / c(total_depleted, total_depleted, total_non_depleted, total_non_depleted)

# Create the bar plot
ggplot(bar_data, aes(x = Group, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Gene Group", y = "Proportion", fill = "Gene Group") +
  # ggtitle("Proportion of PRC Genes among Depleted and Non-Depleted Genes") +
  theme_minimal() + 
  theme(axis.text.x = element_blank(), 
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), 
        panel.border = element_rect(fill = NA, color = "black")) +
  scale_fill_manual(values = c("azure4", "mediumpurple4", "azure3", "mediumpurple1"))

#===
# ASSOCIATION PLOT
#===
if(!require(vcd)){
  install.packages("vcd")
  library(vcd)
}
# Convert the contingency table into a matrix
colnames(contingency_table) = c("PRC genes", "Non-PRC genes")
contingency_matrix <- as.matrix(contingency_table)
contingency_matrix
mosaic(contingency_table, shade=T, legend=T, row_vars = c(1,2), col_vars = c(1,2))
assoc(contingency_matrix, shade=T, legend=T, row_vars = c(1,2), col_vars = c(1,2))










# Create the contingency tables for Case 1 and Case 2
contingency_table_case1 <- matrix(c(5, 91, 18, 446), nrow = 2, byrow = TRUE)
contingency_table_case2 <- matrix(c(7, 60, 16, 477), nrow = 2, byrow = TRUE)

# Combine the two contingency tables into a single table
combined_contingency_table <- array(c(contingency_table_case1, contingency_table_case2), dim = c(2, 2, 2))

# Create row and column names for the combined table
rownames(combined_contingency_table) <- c("depleted_genes", "non_depleted_genes")
colnames(combined_contingency_table) <- c("Case1", "Case2")

# Load the required library
library(vcd)

# Create the mosaic plot for both cases
assoc(as.matrix(combined_contingency_table), 
       shade = TRUE, 
       legend = TRUE, 
       col = c("lightblue", "lightgreen"), 
       split_vertical = TRUE,
       main = "Mosaic Plot of Cases 1 and 2",
       sub = "depleted_genes vs. non_depleted_genes")







setwd("C:/Users/thephillipslab/Documents/Projects/Aurelie_CRISPR_screen/")
path_to_count_table = file.path("./AV_CRISPR_screen_counts.rds")

# Load DGE List containing reads extrated with processAmplicon()
counts <- readRDS(path_to_count_table)
counts %>% head()



# Dictionary mapping samples to their corresponding group of technical replicates
rep_list = list(("HGGx6C_1" = c("Aurelie_3", "Aurelie_4")), 
                ("HGGx6C_2" = c("Aurelie_5", "Aurelie_6"))
             )


rep_list <- list(
  setNames(nm = "HGGx6C_1", c("Aurelie_3", "Aurelie_4")),
  setNames(nm = "HGGx6C_2", c("Aurelie_5", "Aurelie_6"))
)

print(rep_list)
rep_list[1]
rep_list[[1]]
names(rep_list[[1]])[1]

names(rep_list[[1]])
rep_list %>% length()

for (sample in 1:length(rep_list)){
  
  sample_name = as.character(names(rep_list[[sample]])[1])
  replicates = as.character(rep_list[[sample]])
        
  print(sample_name)
  print(replicates)
  
  counts$counts = cbind(counts$counts, rowSums(counts$counts[, replicates]))
  colnames(counts$counts)[length(colnames(counts$counts))] = sample_name
  counts$counts <- counts$counts[, !(colnames(counts$counts) %in% replicates)]

}
head(counts$counts)

# Assuming your counts object is named 'counts'
cbind(counts$counts, "new_col" = rowMeans(counts$counts[, c("Aurelie_4", "Aurelie_5")]))

# Remove the original columns
counts$counts <- counts$counts[, !(colnames(counts$counts) %in% c("Aurelie_4", "Aurelie_5"))]

