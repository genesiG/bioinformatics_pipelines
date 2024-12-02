#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/CAR-T_screen/")
path_to_results_table1 = file.path("./edgeR/GD2_72h_vs_Baseline.txt")
path_to_results_table2 = file.path("./mageck_test_edgeRcounts/GD2_72h_v_Baseline_edgeR.median.sgrna_summary.txt")
#===

library(dplyr)

tab1 = read.table(path_to_results_table1)
tab2 = read.table(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

tab_merge = merge(tab1, tab2, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
tab_merge %>% head()

tab_merge$category = NA
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.all"
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck > 0.1] = "sig.edgeR"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.mageck"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck > 0.1] = "unsig"

tab_merge$enrichment = NA
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC >= 1] = "enrich.all"
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "enrich.edgeR"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC >= 1] = "enrich.mageck"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "unenrich"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC <= -1] = "dep.all"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC > -1 & tab_merge$LFC < 1] = "dep.edgeR"
tab_merge$enrichment[tab_merge$logFC > -1  & tab_merge$LFC < 1 & tab_merge$LFC <= -1] = "dep.mageck"

tab_merge$label = NA
tab_merge$label[tab_merge$category == "sig.all"] = tab_merge$ID[tab_merge$category == "sig.all"] 

tab_sig = filter(tab_merge, category == "sig.all")
tab_sig %>% head()
#write.table(tab_sig, "./GD2_72h_v_Baseline_commonHits.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

tab_pos = filter(tab_merge, logFC >=0 & LFC >= 0)
tab_neg = filter(tab_merge, logFC <0 & LFC < 0)

tab_sig$label = NA
tab_sig$label[tab_sig$enrichment == "enrich.all"] = tab_sig$ID[tab_sig$enrichment == "enrich.all"] 
tab_sig$label[tab_sig$enrichment == "dep.all"] = tab_sig$ID[tab_sig$enrichment == "dep.all"] 

library(ggplot2)
library(ggrepel)

# Create the scatter plot of FDR vs FDR
ggplot(tab_pos, aes(x = -log10(FDR.edgeR), 
                      y = -log10(FDR.mageck), 
                      color = category)) +
  geom_point(aes(alpha = 0.75)) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  geom_vline(xintercept = -log10(0.1), linetype = "dashed") +
  labs(x = "-Log10(FDR) / edgeR", y = "-Log10(FDR) / MAGeCK") +
  scale_x_continuous(limits = c(0, max(-log10(tab_merge$FDR.edgeR) + 0))) +
  scale_y_continuous(limits = c(0, max(-log10(tab_sig$FDR.mageck) + 5))) +
  scale_color_manual(
    values = c("red", "darkblue", "darkblue", "gray"),
    ) +
  ggtitle("Enriched sgRNAs", subtitle = "GD2 72h vs Baseline") + 
  geom_text_repel(aes(label = label)) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, alpha = FALSE)
  
# Create the scatter plot
ggplot(tab_sig, aes(x = logFC, 
                    y = LFC, 
                    color = enrichment)
       ) +
  geom_point(aes(alpha = 0.75)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(x = "-LogFC / edgeR", y = "-LogFC / MAGeCK") +
  # scale_x_continuous(limits = c(0, max(-log10(tab_merge$FDR.edgeR) + 0))) +
  # scale_y_continuous(limits = c(0, max(-log10(tab_sig$FDR.mageck) + 5))) +
  # scale_color_manual(
  #   values = c("purple", "darkgreen", "gray", "gray"),
  # ) +
  ggtitle("Common hits", 
          subtitle = "GD2 72h vs Baseline"
          ) + 
  geom_text_repel(aes(label = label), 
                  max.overlaps = 15, size = 3.5, 
                  #nudge_x = -0.2, nudge_y = -0.2
                  ) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, alpha = FALSE)















#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/edgeR/")
path_to_results_table1 = file.path("./results/RP-conc1-day45_vs_RP-day1.txt")
path_to_results_table2 = file.path("./mageck/DMSO_d45_vs_d1_5cpm_in_15.sgrna_summary.txt")
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)
tab2 = read.table(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

tab_merge = merge(tab1, tab2, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
tab_merge %>% head()

tab_merge$category = NA
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.all"
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck > 0.1] = "sig.edgeR"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.mageck"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck > 0.1] = "unsig"

tab_merge$enrichment = NA
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC >= 1] = "enrich.all"
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "enrich.edgeR"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC >= 1] = "enrich.mageck"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "unenrich"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC <= -1] = "dep.all"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC > -1 & tab_merge$LFC < 1] = "dep.edgeR"
tab_merge$enrichment[tab_merge$logFC > -1  & tab_merge$LFC < 1 & tab_merge$LFC <= -1] = "dep.mageck"

tab_merge$label = NA
tab_merge$label[tab_merge$category == "sig.all"] = tab_merge$ID[tab_merge$category == "sig.all"] 

tab_merge$AveLFC = (tab_merge$logFC + tab_merge$LFC) / 2

tab_sig = filter(tab_merge, category == "sig.all")
tab_sig %>% head()
write.table(tab_sig, "./EZH2i_day45_v_day1_commonHits.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

tab_pos = filter(tab_merge, logFC >=0 & LFC >= 0)
tab_neg = filter(tab_merge, logFC <0 & LFC < 0)

tab_sig$label = NA
tab_sig$label[tab_sig$enrichment == "enrich.all"] = tab_sig$ID[tab_sig$enrichment == "enrich.all"] 
tab_sig$label[tab_sig$enrichment == "dep.all"] = tab_sig$ID[tab_sig$enrichment == "dep.all"] 


# Rank-ordered plot
library(ggplot2)
library(ggrepel)

ranked_tab = tab_merge %>%
  arrange(AveLFC)

ranked_tab$rank = 1:nrow(ranked_tab)
ranked_tab %>% head()
# Get the indexes for the first 10 unique values in the Gene column
unique_genes <- unique(ranked_tab$Gene.edgeR)
depleted_indexes <- match(unique_genes[1:10], ranked_tab$Gene.edgeR)
enriched_indexes <- match(unique_genes[(length(unique_genes)-10):length(unique_genes)], ranked_tab$Gene.edgeR)

# Label top 10 depleted
ranked_tab$label = NA
ranked_tab$label[depleted_indexes] = ranked_tab$Gene.edgeR[depleted_indexes]

# # Label top 10 enriched
# rankedTags$label[enriched_indexes] = rankedTags$Gene[enriched_indexes]

# Color top 10 depleted
ranked_tab$category = NA
ranked_tab$category[depleted_indexes] = "depleted"

# # Color top 10 enriched
# rankedTags$category[enriched_indexes] = "enriched"

ranked_tab %>% filter(Gene.edgeR == "CBX4")
ranked_tab %>% filter(Gene.edgeR == "CBX2")
ranked_tab %>% filter(Gene.edgeR == "SUV39H2" | Gene.edgeR == "KMT1B")

# Get table of average logFC by gene
averages <- aggregate(logFC ~ Gene, data = my_table, FUN = mean)


p = ggplot(data = ranked_tab, aes(x = rank, 
                                  y = AveLFC 
)) +
  geom_point(color = ifelse(is.na(ranked_tab$label), "black", "red"),
             size = ifelse(is.na(ranked_tab$label), 1.5, 2.0)) +
  xlab("Guide Rank") +
  ylab("Average LogFC") +
  geom_text_repel(aes(label = label),
                  color = "red", 
                  max.overlaps = 8, 
                  size = 5,
                  nudge_x = 10, 
                  nudge_y = 0.5
  ) +
  ggtitle("EZH2i (day 45) vs Baseline (day 1)") +
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


pdf(plot_name, width = 6, height = 6)
p
dev.off()




#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
path_to_results_table2 = file.path("./MAGeCK/EZH2i_d45_vs_d1_1cpm_in_10.sgrna_summary.txt")
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)
tab2 = read.delim(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

tab_merge = merge(tab1, tab2, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
tab_merge %>% head()

tab_merge$category = NA
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.all"
tab_merge$category[tab_merge$FDR.edgeR <=0.1 & tab_merge$FDR.mageck > 0.1] = "sig.edgeR"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck <= 0.1] = "sig.mageck"
tab_merge$category[tab_merge$FDR.edgeR >0.1 & tab_merge$FDR.mageck > 0.1] = "unsig"

tab_merge$enrichment = NA
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC >= 1] = "enrich.all"
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "enrich.edgeR"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC >= 1] = "enrich.mageck"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "unenrich"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC <= -1] = "dep.all"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC > -1 & tab_merge$LFC < 1] = "dep.edgeR"
tab_merge$enrichment[tab_merge$logFC > -1  & tab_merge$LFC < 1 & tab_merge$LFC <= -1] = "dep.mageck"

tab_merge$label = NA
tab_merge$label[tab_merge$category == "sig.all"] = tab_merge$ID[tab_merge$category == "sig.all"] 

tab_merge$AveLFC = (tab_merge$logFC + tab_merge$LFC) / 2

tab_sig = filter(tab_merge, category == "sig.all")
tab_sig %>% head()
write.table(tab_sig, "./EZH2i_day45_v_day1_commonHits.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

tab_pos = filter(tab_merge, logFC >=0 & LFC >= 0)
tab_neg = filter(tab_merge, logFC <0 & LFC < 0)

tab_sig$label = NA
tab_sig$label[tab_sig$enrichment == "enrich.all"] = tab_sig$ID[tab_sig$enrichment == "enrich.all"] 
tab_sig$label[tab_sig$enrichment == "dep.all"] = tab_sig$ID[tab_sig$enrichment == "dep.all"] 

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
# Assuming 'tab_merge' contains the merged data
tab_merge$logFDR.edgeR = -log10(tab_merge$FDR.edgeR)
tab_merge$logFDR.mageck = -log10(tab_merge$FDR.mageck)

# Normalize logFDR values for each group (FDR source) to the range [0, 1]
tab_merge_norm <- tab_merge %>%
  mutate(
    logFDR.edgeR_norm = rescale(logFDR.edgeR),
    logFDR.mageck_norm = rescale(logFDR.mageck)
  )

# Reshape the data to long format
tab_merge_long <- tab_merge_norm %>%
  pivot_longer(
    cols = c(logFDR.edgeR_norm, logFDR.mageck_norm),
    names_to = "logFDR_source",
    values_to = "logFDR_value"
  ) %>%
  mutate(logFDR_source = gsub("logFDR.", "", logFDR_source)) 


# Create the dotplot using ggplot2
ggplot(tab_merge_long, aes(x = Gene.edgeR, y = logFDR_value, color = logFDR_source, group = ID)) +
  geom_point(alpha = 0.5, 
             position = position_dodge(width = 0)) +
  xlab("Gene") +
  ylab(expression("Normalized -log"[10]*" (FDR) / edgeR")) +
  ggtitle("Dotplot of -log10(FDR) for each ID within each Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    text = element_text(size = 12, face = "bold"),
    axis.ticks.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(), 
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_y_continuous(
    sec.axis = sec_axis(trans = ~ ., name = expression("Normalized -log"[10]*" (FDR) / MAGeCK"))
  ) +
  guides(color = guide_legend(title = "FDR Source", override.aes = list(shape = 16)))







ctrl = filter(tab_merge_norm, grepl("Control", tab_merge_norm$Gene.edgeR))
median_LFC_ctrl_edgeR = median(ctrl$logFC)
median_LFC_ctrl_mageck = median(ctrl$LFC)

tab_merge_norm$score.edgeR = tab_merge_norm$logFC - median_LFC_ctrl_edgeR
tab_merge_norm$score.mageck = tab_merge_norm$logFC - median_LFC_ctrl_mageck


# Reshape the data to long format with logFC and LFC columns
tab_merge_long <- tab_merge_norm %>%
  pivot_longer(
    cols = c(score.edgeR, score.mageck),
    names_to = "score_source",
    values_to = "score_value"
  ) %>%
  mutate(score_source = gsub("score.", "", score_source))

# Create the dotplot using ggplot2
min = min(tab_merge_long$score_value)
max = -min(tab_merge_long$score_value)

ggplot(tab_merge_long, aes(x = Gene.edgeR, y = score_value, color = score_source, group = ID)) +
  geom_point(alpha = 0.5, position = position_dodge(width = 0)) +
  xlab("Gene") +
  ylab(expression("CRISPR score / edgeR")) +
  ggtitle("Dotplot of logFC and LFC for each ID within each Gene") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    text = element_text(size = 12, face = "bold"),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(#limits = c(min, max),
                     sec.axis = sec_axis(trans = ~ ., name = expression("CRISPR score/ MAGeCK"))
  ) +
  guides(color = guide_legend(title = "", override.aes = list(shape = 16)))


file_name = paste(c(treat, "vs", control, "rankedPlot"), collapse = "_")
plot_name = paste0("./plots/", file_name, ".pdf")


pdf(plot_name, width = 6, height = 6)
p
dev.off()



#===
# RANK-ORDERED PLOT
#===
#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
path_to_results_table2 = file.path("./MAGeCK/EZH2i_d45_vs_d1_1cpm_in_10.sgrna_summary.txt")
treat = "EZH2i_day45"
control = "day1"
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)
tab2 = read.delim(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

# Get median LFC from control sgRNAs
ctrl_edgeR = filter(tab1, grepl("Control", tab1$Gene))
ctrl_mageck = filter(tab2, grepl("Control", tab2$Gene))

median_edgeR = median(ctrl_edgeR$logFC)
median_mageck = median(ctrl_mageck$LFC)

# Get CRISPR scores
tab1$score = tab1$logFC - median_edgeR
tab2$score = tab2$LFC - median_mageck

# Get table of average logFC by gene
tab1_avg <- aggregate(score ~ Gene, data = tab1, FUN = mean)
tab2_avg <- aggregate(score ~ Gene, data = tab2, FUN = mean)

tab1_avg %>% head()
tab2_avg %>% head()

# Merge tables
tab_merge = merge(tab1_avg, tab2_avg, by = "Gene", suffixes = c(".edgeR",".mageck"))
tab_merge$AveScore = (tab_merge$score.edgeR + tab_merge$score.mageck)/2
tab_merge %>% head()

# Rank genes by average score
ranked_tab = tab_merge %>% arrange(AveScore)
ranked_tab$rank = 1:nrow(ranked_tab)

# Add labels
ranked_tab$label = NA
genes_of_interest = c("BMI1", "CBX2", "CDK9", "RPA3", "CBX4")
ranked_tab$label[ranked_tab$Gene %in% genes_of_interest] = ranked_tab$Gene[ranked_tab$Gene %in% genes_of_interest]

# Rank-ordered plot
p = ggplot(data = ranked_tab, aes(x = rank, 
                                  y = AveScore 
)) +
  geom_point(color = "black",
    #color = ifelse(is.na(ranked_tab$label), "black", "red"),
             #fill = ifelse(ranked_tab$label == "CBX2", "red", "black"),
             size = 1.5) +
  geom_point(color = ifelse(ranked_tab$label == "CBX2", "red", NA),
             size = 2.0) +
  xlab("Gene Rank") +
  ylab("Average CRISPR score") +
  geom_text_repel(aes(label = label),
                  color = ifelse(ranked_tab$label == "CBX2", "red", "black"), 
                  max.overlaps = 8, 
                  size = 4.5,
                  xlim = c(60,500), 
                  ylim = c(-8,-2),
                  force_pull = 0.5,
                  #nudge_x = 10, 
                  nudge_y = -1.5
  ) +
  ggtitle("EZH2i vs baseline") +
  theme_bw() + 
  theme(text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, size = FALSE)


p


file_name = paste(c(treat, "vs", control, "1_in_10cpm", "rankedPlot"), collapse = "_")
plot_name = paste0("./plots/", file_name, ".pdf")


pdf(plot_name, width = 6, height = 6)
p
dev.off()

# FDR vs FDR plot to see common hits 


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

install.packages("ggplot2")
install.packages("plotly")

# Load required libraries
install.packages("ggplot2")
install.packages("plotly")
library(ggplot2)
library(plotly)

# Calculate percentages for each group in the Function column
function_counts <- table(screenGenes$Function)
function_percent <- round(100 * function_counts / sum(function_counts), 2)

# Create a data frame for the pie chart
pie_data <- data.frame(Function = names(function_counts), Percentage = function_percent)
pie_data %>% head()
# Create the pie chart using ggplot2 and plotly
pie_chart <- ggplot(pie_data, aes(x = "", y = Percentage.Freq, fill = factor(Function))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Percentage of Different Groups in Function",
       fill = "Function",
       x = NULL,
       y = NULL)

print(pie_chart)
# Convert ggplot2 chart to plotly
plotly_chart <- ggplotly(pie_chart)

# Display the plotly chart
print(plotly_chart)



screenGenes = filter(screenGenes, Function != "#")

# Calculate percentages for each group in the Function column
function_counts <- table(screenGenes$Function)
total_counts <- sum(function_counts)
function_percent <- round(100 * function_counts / total_counts, 2)

# Filter out data points with less than 1% counts
filtered_counts <- function_counts[function_percent >= 1]
filtered_percent <- function_percent[function_percent >= 1]

# Define custom colors for the chart
library("RColorBrewer")

custom_colors = colorRampPalette(brewer.pal(12,"Set3"))(length(names(filtered_counts))) %>% 
  rev()

color_annotation = data.frame(Function = names(filtered_counts), Colors = custom_colors)
# Create a data frame with names and percentages
pie_data <- data.frame(Function = names(filtered_counts), Percentage = filtered_percent)

# Arrange the data frame based on the percentages in descending order
pie_data <- pie_data[order(pie_data$Percentage.Freq, decreasing = T), ]

# Create the pie chart using plotly with the ordered data
plotly_chart <- plot_ly(labels = pie_data$Function, 
                        values = pie_data$Percentage.Freq, 
                        type = "pie",
                        marker = list(colors = custom_colors),
                        textinfo = "percent",
                        hoverinfo = "text",
                        textposition = "outside",
                        domain = list(x = c(0, 0.45), y = c(0, 1)))

plotly_chart #%>% hide_legend()


# Layout for the pie chart
layout <- list(
  showlegend = FALSE,
  margin = list(l = 10, r = 10, t = 10, b = 10),
  xaxis = list(showticklabels = FALSE, zeroline = FALSE),
  yaxis = list(showticklabels = FALSE, zeroline = FALSE)
)

# Create a separate figure for the legend
legend_chart <- plot_ly(showlegend = TRUE, 
                        type = "pie",
                        hole = 0.4,
                        marker = list(colors = custom_colors),
                        labels = pie_data$Function,
                        values = pie_data$Percentage.Freq,
                        domain = list(x = c(0.55, 1), y = c(0, 1)),
                        textinfo = "label",
                        hoverinfo = "text",
                        text = ~paste("Function: ", labels)
)

# Combine the pie chart and the legend chart using subplot
subplot_chart <- subplot(plotly_chart, legend_chart, layout = layout)

# Display the combined chart
print(subplot_chart)

#===
# FDR VS FDR PLOT
#===
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
path_to_results_table2 = file.path("./MAGeCK/EZH2i_d45_vs_d1_1cpm_in_10.sgrna_summary.txt")
treat = "EZH2i_day45"
control = "day1"
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)
tab2 = read.delim(path_to_results_table2, header = T)

tab1 %>% head()
tab2 %>% head()

# Get median LFC from control sgRNAs
ctrl_edgeR = filter(tab1, grepl("Control", tab1$Gene))
ctrl_mageck = filter(tab2, grepl("Control", tab2$Gene))

median_edgeR = median(ctrl_edgeR$logFC)
median_mageck = median(ctrl_mageck$LFC)

# Get CRISPR scores
tab1$score = tab1$logFC - median_edgeR
tab2$score = tab2$LFC - median_mageck

alpha1 = (0.05 / nrow(tab1))
alpha2 = (0.05 / nrow(tab2))
# tab1_sig = filter(tab1, FDR <= (0.05 / nrow(tab1)))
# tab2_sig = filter(tab2, FDR <= (0.05 / nrow(tab2)))
# 
# tab1 %>% nrow()
# tab1_sig %>% nrow()
# tab2 %>% nrow()
# tab2_sig %>% nrow()

# tab_sig_merge = merge(tab1_sig, tab2_sig, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
# tab_sig_merge %>% head()
# 
# tab_sig_merge$category = NA
# tab_sig_merge$category[tab_sig_merge$FDR.edgeR <=0.1 & tab_sig_merge$FDR.mageck <= 0.1] = "sig.all"
# tab_sig_merge$category[tab_sig_merge$FDR.edgeR <=0.1 & tab_sig_merge$FDR.mageck > 0.1] = "sig.edgeR"
# tab_sig_merge$category[tab_sig_merge$FDR.edgeR >0.1 & tab_sig_merge$FDR.mageck <= 0.1] = "sig.mageck"
# tab_sig_merge$category[tab_sig_merge$FDR.edgeR >0.1 & tab_sig_merge$FDR.mageck > 0.1] = "unsig"
# 
# tab_sig_merge$enrichment = NA
# tab_sig_merge$enrichment[tab_sig_merge$logFC >= 1 & tab_sig_merge$LFC >= 1] = "enrich.all"
# tab_sig_merge$enrichment[tab_sig_merge$logFC >= 1 & tab_sig_merge$LFC < 1 & tab_sig_merge$LFC > -1] = "enrich.edgeR"
# tab_sig_merge$enrichment[tab_sig_merge$logFC < 1 & tab_sig_merge$logFC > -1 & tab_sig_merge$LFC >= 1] = "enrich.mageck"
# tab_sig_merge$enrichment[tab_sig_merge$logFC < 1 & tab_sig_merge$logFC > -1 & tab_sig_merge$LFC < 1 & tab_sig_merge$LFC > -1] = "unenrich"
# tab_sig_merge$enrichment[tab_sig_merge$logFC <= -1 & tab_sig_merge$LFC <= -1] = "dep.all"
# tab_sig_merge$enrichment[tab_sig_merge$logFC <= -1 & tab_sig_merge$LFC > -1 & tab_sig_merge$LFC < 1] = "dep.edgeR"
# tab_sig_merge$enrichment[tab_sig_merge$logFC > -1  & tab_sig_merge$LFC < 1 & tab_sig_merge$LFC <= -1] = "dep.mageck"
# 
# tab_sig_merge$label = NA
# tab_sig_merge$label[tab_sig_merge$category == "sig.all"] = tab_sig_merge$ID[tab_sig_merge$category == "sig.all"] 
# 
# tab_sig = filter(tab_sig_merge, category == "sig.all")
# tab_sig %>% head()
# tab_sig$label = NA
# tab_sig$label[tab_sig$enrichment == "enrich.all"] = tab_sig$ID[tab_sig$enrichment == "enrich.all"] 
# tab_sig$label[tab_sig$enrichment == "dep.all"] = tab_sig$ID[tab_sig$enrichment == "dep.all"] 

tab_merge = merge(tab1, tab2, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck")) 
tab_neg = filter(tab_merge, score.edgeR < 0 & score.mageck < 0)

tab_neg$category = NA
tab_neg$category[tab_neg$FDR.edgeR <= alpha1 & tab_neg$FDR.mageck <= alpha2] = "sig.all"
tab_neg$category[tab_neg$FDR.edgeR <= alpha1 & tab_neg$FDR.mageck > alpha2] = "sig.edgeR"
tab_neg$category[tab_neg$FDR.edgeR >alpha1 & tab_neg$FDR.mageck <= alpha2] = "sig.mageck"
tab_neg$category[tab_neg$FDR.edgeR >alpha1 & tab_neg$FDR.mageck > alpha2] = "unsig"

tab_neg$enrichment = NA
tab_neg$enrichment[tab_neg$logFC >= 1 & tab_neg$LFC >= 1] = "enrich.all"
tab_neg$enrichment[tab_neg$logFC >= 1 & tab_neg$LFC < 1 & tab_neg$LFC > -1] = "enrich.edgeR"
tab_neg$enrichment[tab_neg$logFC < 1 & tab_neg$logFC > -1 & tab_neg$LFC >= 1] = "enrich.mageck"
tab_neg$enrichment[tab_neg$logFC < 1 & tab_neg$logFC > -1 & tab_neg$LFC < 1 & tab_neg$LFC > -1] = "unenrich"
tab_neg$enrichment[tab_neg$logFC <= -1 & tab_neg$LFC <= -1] = "dep.all"
tab_neg$enrichment[tab_neg$logFC <= -1 & tab_neg$LFC > -1 & tab_neg$LFC < 1] = "dep.edgeR"
tab_neg$enrichment[tab_neg$logFC > -1  & tab_neg$LFC < 1 & tab_neg$LFC <= -1] = "dep.mageck"

tab_neg$label = NA
tab_neg$label[tab_neg$category == "sig.all"] = tab_neg$ID[tab_neg$category == "sig.all"] 

tab_neg_function = merge(tab_neg, screenGenes[,c("HGNC.approved.symbol", "Protein.complex")], 
                         by.x = "Gene.edgeR", 
                         by.y = "HGNC.approved.symbol", all.x = T, all.y = F)
tab_neg %>% nrow()
tab_neg_function %>% nrow()

fun_counts = function_counts <- table(tab_neg_function$Protein.complex)
fun_counts = fun_counts[order(fun_counts, decreasing = T)]
# fun_counts = fun_counts[names(fun_counts) != "Chromatin remodeling"]
fun_counts %>% head()

# Create an empty column 'colors' in the tab_neg_function data frame
# # Loop through each row
# for (i in seq_len(nrow(tab_neg_function))) {
#   # Check if the value in the 'Function' column is in the first 4 names of fun_counts
#   if (tab_neg_function$Protein.complex[i] %in% names(fun_counts)[1:5]) {
#     # Set the value of 'colors' to the corresponding value in the 'Function' column
#     tab_neg_function$colors[i] <- as.character(tab_neg_function$Protein.complex[i])
#   } else {
#     tab_neg_function$colors[i] = NA
#   }
# }

tab_neg_function$colors <- NA
tab_neg_function$colors[grepl("PRC1",tab_neg_function$Protein.complex)] <- "PRC1"
tab_neg_function$colors[grepl("PRC2",tab_neg_function$Protein.complex)] <- "PRC2"

tab_neg_function$label = NA
#tab_neg_function$label[!is.na(tab_neg_function$colors)] <- tab_neg_function$ID[!is.na(tab_neg_function$colors)]
tab_neg_function$label[grepl("PRC",tab_neg_function$Protein.complex)] <- tab_neg_function$ID[grepl("PRC",tab_neg_function$Protein.complex)]
# tab_neg_function$label[grepl("PRC",tab_neg_function$Protein.complex)] <- tab_neg_function$ID[grepl("PRC",tab_neg_function$Protein.complex)]


filter(tab_neg_function, grepl("PRC", tab_neg_function$Protein.complex)) %>% head()

tab_neg_function_filt = filter(tab_neg_function, !is.na(tab_neg_function$colors))

p = ggplot(tab_neg_function, aes(x = -log10(FDR.edgeR), 
                    y = -log10(FDR.mageck), 
                    #color = colors
                    )
       ) +
  # geom_rect(data = NULL, 
  #           fill = alpha("purple", 0.2), 
  #           xmin = -Inf, xmax = -log10(alpha), ymin = -log10(alpha2), ymax = Inf
  # ) +  # Rectangle 1: x < alpha, y > alpha
  # geom_rect(data = NULL, 
  #           fill = alpha("blue", 0.2),
  #           xmin = -log10(alpha), xmax = Inf, ymin = -Inf, ymax = -log10(alpha)
  # ) +  # Rectangle 2: x > alpha, y < alpha
  # geom_rect(data = NULL, 
  #           fill = alpha("red", 0.2), 
  #           xmin = -Inf, xmax = -log10(alpha), ymin = -Inf, ymax = -log10(alpha)
  # ) +  # Rectangle 3: x < alpha, y < alpha
  geom_point(size = 1.5, color = "gray") +
  geom_point(data = tab_neg_function_filt, mapping = aes(color = colors),size = 2) +
  geom_hline(yintercept = -log10(0.05/644), linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05/644), linetype = "dashed") +
  labs(x = "-Log10(FDR) / edgeR", y = "-Log10(FDR) / MAGeCK") +
  scale_x_continuous(limits = c(0, max(-log10(tab_neg$FDR.edgeR) + 0))) +
  scale_y_continuous(limits = c(0, max(-log10(tab_neg$FDR.mageck) + 0))) +
  scale_color_manual(
    values = c("red", "purple", "orange", "green"), 
    na.value = "black"
  ) +
  scale_fill_manual(
    values = c("purple" = alpha("purple", 0.2), 
               "blue" = alpha("blue", 0.2), 
               "red" = alpha("red", 0.2))
  ) +
  geom_text_repel(aes(label = label)) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  guides(#color = FALSE, 
         alpha = FALSE, 
         fill = FALSE)


p

file_name = paste(c(treat, "vs", control, "1_in_10cpm", "FDR_vs_FDR"), collapse = "_")
plot_name = paste0("./plots/", file_name, ".pdf")


pdf(plot_name, width = 8.5, height = 6)
p
dev.off()


#===
# LFC vs LFC PLOT
#===
tab_merge = merge(tab1, tab2, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck")) 

tab_merge$category = NA
tab_merge$category[tab_merge$FDR.edgeR <= alpha1 & tab_merge$FDR.mageck <= alpha2] = "sig.all"
tab_merge$category[tab_merge$FDR.edgeR <= alpha1 & tab_merge$FDR.mageck > alpha2] = "sig.edgeR"
tab_merge$category[tab_merge$FDR.edgeR >alpha1 & tab_merge$FDR.mageck <= alpha2] = "sig.mageck"
tab_merge$category[tab_merge$FDR.edgeR >alpha1 & tab_merge$FDR.mageck > alpha2] = "unsig"

tab_merge$enrichment = NA
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC >= 1] = "enrich.all"
tab_merge$enrichment[tab_merge$logFC >= 1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "enrich.edgeR"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC >= 1] = "enrich.mageck"
tab_merge$enrichment[tab_merge$logFC < 1 & tab_merge$logFC > -1 & tab_merge$LFC < 1 & tab_merge$LFC > -1] = "unenrich"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC <= -1] = "dep.all"
tab_merge$enrichment[tab_merge$logFC <= -1 & tab_merge$LFC > -1 & tab_merge$LFC < 1] = "dep.edgeR"
tab_merge$enrichment[tab_merge$logFC > -1  & tab_merge$LFC < 1 & tab_merge$LFC <= -1] = "dep.mageck"

tab_merge$label = NA
tab_merge$label[tab_merge$category == "sig.all"] = tab_merge$ID[tab_merge$category == "sig.all"] 

tab_merge_function = merge(tab_merge, screenGenes[,c("HGNC.approved.symbol", "Protein.complex")], 
                         by.x = "Gene.edgeR", 
                         by.y = "HGNC.approved.symbol", all.x = T, all.y = F)
tab_merge %>% nrow()
tab_merge_function %>% nrow()

fun_counts = function_counts <- table(tab_merge_function$Protein.complex)
fun_counts = fun_counts[order(fun_counts, decreasing = T)]
# fun_counts = fun_counts[names(fun_counts) != "Chromatin remodeling"]
fun_counts %>% head()

tab_merge_function$colors <- NA
tab_merge_function$colors[grepl("PRC1",tab_merge_function$Protein.complex)] <- "PRC1"
tab_merge_function$colors[grepl("PRC2",tab_merge_function$Protein.complex)] <- "PRC2"
tab_merge_function$colors[tab_merge_function$logFC > -1] <- NA
tab_merge_function$colors[tab_merge_function$LFC > -1] <- NA

tab_merge_function$label = NA
#tab_merge_function$label[!is.na(tab_merge_function$colors)] <- tab_merge_function$ID[!is.na(tab_merge_function$colors)]
tab_merge_function$label[grepl("PRC",tab_merge_function$Protein.complex)] <- tab_merge_function$ID[grepl("PRC",tab_merge_function$Protein.complex)]
# tab_merge_function$label[grepl("PRC",tab_merge_function$Protein.complex)] <- tab_merge_function$ID[grepl("PRC",tab_merge_function$Protein.complex)]
tab_merge_function$label = NA
tab_merge_function$label[!is.na(tab_merge_function$colors)] <- tab_merge_function$ID[!is.na(tab_merge_function$colors)]

filter(tab_merge_function, grepl("PRC", tab_merge_function$Protein.complex)) %>% head()

tab_merge_function_filt = filter(tab_merge_function, !is.na(tab_merge_function$colors))


p = ggplot(tab_merge_function, aes(x = logFC, 
                                 y = LFC, 
                                 #color = colors
                                 )
           ) +
  # geom_rect(data = NULL, 
  #           fill = alpha("purple", 0.2), 
  #           xmin = -Inf, xmax = -log10(alpha), ymin = -log10(alpha2), ymax = Inf
  # ) +  # Rectangle 1: x < alpha, y > alpha
  # geom_rect(data = NULL, 
  #           fill = alpha("blue", 0.2),
  #           xmin = -log10(alpha), xmax = Inf, ymin = -Inf, ymax = -log10(alpha)
  # ) +  # Rectangle 2: x > alpha, y < alpha
  # geom_rect(data = NULL, 
  #           fill = alpha("red", 0.2), 
  #           xmin = -Inf, xmax = -log10(alpha), ymin = -Inf, ymax = -log10(alpha)
# ) +  # Rectangle 3: x < alpha, y < alpha
  geom_point(size = 1.5, color = "gray") +
  geom_point(data = tab_merge_function_filt, mapping = aes(color = colors),size = 2) +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "LogFC / edgeR", y = "LogFC / MAGeCK") +
  # scale_x_continuous(limits = c(0, max(-log10(tab_neg$FDR.edgeR) + 0))) +
  # scale_y_continuous(limits = c(0, max(-log10(tab_neg$FDR.mageck) + 0))) +
  scale_color_manual(
    values = c("red", "purple", "orange", "green"), 
    na.value = "black"
  ) +
  scale_fill_manual(
    values = c("purple" = alpha("purple", 0.2), 
               "blue" = alpha("blue", 0.2), 
               "red" = alpha("red", 0.2))
  ) +
  geom_text_repel(aes(label = label)) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"), 
        legend.title = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_rect(fill = NA, color = "black")) +
  guides(#color = FALSE, 
    alpha = FALSE, 
    fill = FALSE)


p




#===
# COMMON AND UNIQUE HITS DMSO VS EZH2i
#===
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
path_to_results_table2 = file.path("./MAGeCK/EZH2i_d45_vs_d1_1cpm_in_10.sgrna_summary.txt")
path_to_results_table3 = file.path("./edgeR/RP-DMSO-day45_vs_RP-day1.txt")
path_to_results_table4 = file.path("./MAGeCK/DMSO_d45_vs_d1_1cpm_in_10.sgrna_summary.txt")
# treat = "EZH2i_day45"
# control = "day1"
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)
tab2 = read.delim(path_to_results_table2, header = T)
tab3 = read.table(path_to_results_table3, header = T)
tab4 = read.delim(path_to_results_table4, header = T)

tab1 %>% head()
tab2 %>% head()

# Get median LFC from control sgRNAs
ctrl_edgeR_treat = filter(tab1, grepl("Control", tab1$Gene))
ctrl_mageck_treat = filter(tab2, grepl("Control", tab2$Gene))
ctrl_edgeR_dmso = filter(tab1, grepl("Control", tab1$Gene))
ctrl_mageck_dmso = filter(tab2, grepl("Control", tab2$Gene))

median_edgeR_treat = median(ctrl_edgeR_treat$logFC)
median_mageck_treat = median(ctrl_mageck_treat$LFC)
median_edgeR_dmso = median(ctrl_edgeR_dmso$logFC)
median_mageck_dmso = median(ctrl_mageck_dmso$LFC)

# Get CRISPR scores
tab1$score = tab1$logFC - median_edgeR_treat
tab2$score = tab2$LFC - median_mageck_treat
tab3$score = tab3$logFC - median_edgeR_dmso
tab4$score = tab4$LFC - median_mageck_dmso

alpha1 = (0.05 / nrow(tab1))
alpha2 = (0.05 / nrow(tab2))
alpha3 = (0.05 / nrow(tab3))
alpha4 = (0.05 / nrow(tab4))

tab1_sig = filter(tab1, FDR <= alpha1)
tab2_sig = filter(tab2, FDR <= alpha2)
tab3_sig = filter(tab3, FDR <= alpha3)
tab4_sig = filter(tab4, FDR <= alpha4)

tab1 %>% nrow()
tab1_sig %>% nrow()
tab2 %>% nrow()
tab2_sig %>% nrow()
tab3 %>% nrow()
tab3_sig %>% nrow()

tab_sig_treat = merge(tab1_sig, tab2_sig, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
tab_sig_treat %>% head()
tab_sig_dmso = merge(tab3_sig, tab4_sig, by.x = "ID", by.y = "sgrna", suffixes = c(".edgeR",".mageck"))
tab_sig_dmso %>% head()


# EZH2i
tab_sig_treat$category = NA
tab_sig_treat$category[tab_sig_treat$FDR.edgeR <=0.1 & tab_sig_treat$FDR.mageck <= 0.1] = "sig.all"
tab_sig_treat$category[tab_sig_treat$FDR.edgeR <=0.1 & tab_sig_treat$FDR.mageck > 0.1] = "sig.edgeR"
tab_sig_treat$category[tab_sig_treat$FDR.edgeR >0.1 & tab_sig_treat$FDR.mageck <= 0.1] = "sig.mageck"
tab_sig_treat$category[tab_sig_treat$FDR.edgeR >0.1 & tab_sig_treat$FDR.mageck > 0.1] = "unsig"

tab_sig_treat$enrichment = NA
tab_sig_treat$enrichment[tab_sig_treat$logFC >= 1 & tab_sig_treat$LFC >= 1] = "enrich.all"
tab_sig_treat$enrichment[tab_sig_treat$logFC >= 1 & tab_sig_treat$LFC < 1 & tab_sig_treat$LFC > -1] = "enrich.edgeR"
tab_sig_treat$enrichment[tab_sig_treat$logFC < 1 & tab_sig_treat$logFC > -1 & tab_sig_treat$LFC >= 1] = "enrich.mageck"
tab_sig_treat$enrichment[tab_sig_treat$logFC < 1 & tab_sig_treat$logFC > -1 & tab_sig_treat$LFC < 1 & tab_sig_treat$LFC > -1] = "unenrich"
tab_sig_treat$enrichment[tab_sig_treat$logFC <= -1 & tab_sig_treat$LFC <= -1] = "dep.all"
tab_sig_treat$enrichment[tab_sig_treat$logFC <= -1 & tab_sig_treat$LFC > -1 & tab_sig_treat$LFC < 1] = "dep.edgeR"
tab_sig_treat$enrichment[tab_sig_treat$logFC > -1  & tab_sig_treat$LFC < 1 & tab_sig_treat$LFC <= -1] = "dep.mageck"

tab_sig_treat$label = NA
tab_sig_treat$label[tab_sig_treat$category == "sig.all"] = tab_sig_treat$ID[tab_sig_treat$category == "sig.all"]

tab_sig_treat_final = filter(tab_sig_treat, category == "sig.all")
tab_sig_treat_final %>% nrow()

# DMSO
tab_sig_dmso$category = NA
tab_sig_dmso$category[tab_sig_dmso$FDR.edgeR <=0.1 & tab_sig_dmso$FDR.mageck <= 0.1] = "sig.all"
tab_sig_dmso$category[tab_sig_dmso$FDR.edgeR <=0.1 & tab_sig_dmso$FDR.mageck > 0.1] = "sig.edgeR"
tab_sig_dmso$category[tab_sig_dmso$FDR.edgeR >0.1 & tab_sig_dmso$FDR.mageck <= 0.1] = "sig.mageck"
tab_sig_dmso$category[tab_sig_dmso$FDR.edgeR >0.1 & tab_sig_dmso$FDR.mageck > 0.1] = "unsig"

tab_sig_dmso$enrichment = NA
tab_sig_dmso$enrichment[tab_sig_dmso$logFC >= 1 & tab_sig_dmso$LFC >= 1] = "enrich.all"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC >= 1 & tab_sig_dmso$LFC < 1 & tab_sig_dmso$LFC > -1] = "enrich.edgeR"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC < 1 & tab_sig_dmso$logFC > -1 & tab_sig_dmso$LFC >= 1] = "enrich.mageck"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC < 1 & tab_sig_dmso$logFC > -1 & tab_sig_dmso$LFC < 1 & tab_sig_dmso$LFC > -1] = "unenrich"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC <= -1 & tab_sig_dmso$LFC <= -1] = "dep.all"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC <= -1 & tab_sig_dmso$LFC > -1 & tab_sig_dmso$LFC < 1] = "dep.edgeR"
tab_sig_dmso$enrichment[tab_sig_dmso$logFC > -1  & tab_sig_dmso$LFC < 1 & tab_sig_dmso$LFC <= -1] = "dep.mageck"

tab_sig_dmso$label = NA
tab_sig_dmso$label[tab_sig_dmso$category == "sig.all"] = tab_sig_dmso$ID[tab_sig_dmso$category == "sig.all"]

tab_sig_dmso_final = filter(tab_sig_dmso, category == "sig.all")
tab_sig_dmso_final %>% nrow()


commonHits = filter(tab_sig_treat_final, ID %in% tab_sig_dmso_final$ID)
commonHits %>% nrow()

dmso_only = filter(tab_sig_dmso_final, !ID %in% tab_sig_treat_final$ID)
nrow(dmso_only) + nrow(commonHits) == nrow(tab_sig_dmso_final)

treat_only = filter(tab_sig_treat_final, !ID %in% tab_sig_dmso_final$ID)
nrow(treat_only) + nrow(commonHits) == nrow(tab_sig_treat_final)

nrow(dmso_only)
nrow(commonHits)
nrow(treat_only)

treat_only_function = merge(treat_only[,c("ID","Gene.edgeR")], screenGenes[,c("HGNC.approved.symbol", "Function")], 
                            by.x = "Gene.edgeR", by.y = "HGNC.approved.symbol",
                            all.x = T,
                            all.y = F)
treat_only_function %>% nrow()
treat_only_function %>% head()

# BAR PLOT OF EPIFACTOR ANNOTATION
library(dplyr)
library(ggplot2)

# Calculate total counts for each value in the 'Function' column
function_counts <- treat_only_function %>%
  count(Function, name = "Counts") %>%
  arrange(desc(Counts))

# Calculate percentage of each value in the 'Function' column
total_counts <- sum(function_counts$Counts)
function_counts <- function_counts %>%
  mutate(Percentage = (Counts / total_counts) * 100)

# Remove rows with NA in the 'Function' column
function_counts <- function_counts %>%
  filter(!is.na(Function))

# Merge function_counts with color_annotation based on the 'Function' column
function_counts <- merge(function_counts, color_annotation, by = "Function", all.x = TRUE)

# Create the bar plot using ggplot2 with custom colors
ggplot(function_counts, aes(x = reorder(Function, -Counts), y = Counts, fill = Colors)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(Counts, "%", sep = "")), position = position_stack(vjust = 1.05)) +
  xlab("") +
  ylab("Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_manual(values = function_counts$Colors)
# Create the bar plot using ggplot2
p = ggplot(function_counts, aes(x = reorder(Function, -Counts), y = Counts, fill = Function)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste(Counts, "%", sep = "")), position = position_stack(vjust = 1.05)) +
  xlab("") +
  ylab("Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    text = element_text(size = 12),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

p

pdf("./plots/EZH2i_only_epifactor.pdf", width = 13, height = 9)
p
dev.off()















# Create the scatter plot
ggplot(tab_sig, aes(x = logFC, 
                    y = LFC, 
                    color = enrichment)
) +
  geom_point(aes(alpha = 0.75)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(x = "-LogFC / edgeR", y = "-LogFC / MAGeCK") +
  # scale_x_continuous(limits = c(0, max(-log10(tab_sig_merge$FDR.edgeR) + 0))) +
  # scale_y_continuous(limits = c(0, max(-log10(tab_sig$FDR.mageck) + 5))) +
  # scale_color_manual(
  #   values = c("purple", "darkgreen", "gray", "gray"),
  # ) +
  ggtitle("Common hits", 
          subtitle = "GD2 72h vs Baseline"
  ) + 
  geom_text_repel(aes(label = label), 
                  max.overlaps = 15, size = 3.5, 
                  #nudge_x = -0.2, nudge_y = -0.2
  ) +
  theme_minimal() + 
  theme(axis.title = element_text(face = "bold"),
        title = element_text(face = "bold")) +
  guides(color = FALSE, alpha = FALSE)





#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
treat = "EZH2i_day45"
control = "day1"
#===

library(dplyr)

tab1 = read.table(path_to_results_table1, header = T)

tab1 %>% head()

# Get median LFC from control sgRNAs
ctrl = filter(tab1, grepl("Control", tab1$Gene))

medianLFC = median(ctrl_edgeR$logFC)

# Get CRISPR scores
tab1$score = tab1$logFC - medianLFC

# Get table of average logFC by gene
tab1_avg <- aggregate(score ~ Gene, data = tab1, FUN = mean)

tab1_avg %>% head()


# Rank genes by average score
ranked_tab = tab1_avg %>% arrange(score)
ranked_tab$rank = 1:nrow(ranked_tab)

# Add labels
ranked_tab$label = NA
genes_of_interest = c("BMI1", "CBX2", "CDK9", "RPA3", "CBX4")
ranked_tab$label[ranked_tab$Gene %in% genes_of_interest] = ranked_tab$Gene[ranked_tab$Gene %in% genes_of_interest]

# Rank-ordered plot
p = ggplot(data = ranked_tab, aes(x = rank, 
                                  y = score 
)) +
  geom_point(color = ifelse(is.na(ranked_tab$label), "black", "red"),
             size = ifelse(is.na(ranked_tab$label), 1.5, 2.0)) +
  xlab("Gene Rank") +
  ylab("Average CRISPR score") +
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
  ggtitle("EZH2i (day 45) vs Baseline (day 1)") +
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


pdf(plot_name, width = 6, height = 6)
p
dev.off()





#===
#===
# CHANGE FOR EACH PROJECT
#===
setwd("C:/Users/thephillipslab/Documents/Projects/screenResults/")
path_to_results_table1 = file.path("./edgeR/RP-conc1-day45_vs_RP-day1.txt")
treat = "EZH2i_day45"
control = "day1"
#===

tab1 = read.table(path_to_results_table1, header = T)

tab1 %>% head()

# Get median LFC from control sgRNAs
ctrl = filter(tab1, grepl("Control", tab1$Gene))

medianLFC = median(ctrl$logFC)

# Get CRISPR scores
tab1$score = tab1$logFC - medianLFC

tab1$category = NA
tab1$category[tab1$score > 1] = "enrich"
tab1$category[tab1$score < -1] = "dep"

# Set plot title for each comparison
comparison = "EZH2i vs baseline"
min_x = min(tab1$score)
max_x = -min_x

# Make volcano plot
volcanoPlot <- ggplot(data=tab1, 
                      aes(x=score, 
                          y=-log10(FDR),
                          color = category, 
                          #shape = Peak, 
                          #label =  Label
                          )
                      ) + 
  geom_point() + 
  theme_bw() +
  # scale_shape_manual(values = c(1,16,4,17)) +
  ylab(expression(bold(bolditalic(-Log[10]) ~ "FDR"))) +
  # xlab(expression(bold(bolditalic(Log[2]) ~ "FC"))) +
  xlab("sgRNA CRISPR score") +
  ggtitle(comparison) + 
  scale_x_continuous(n.breaks = 7, 
                     limits = c(min_x,max_x)
  ) +
  theme(legend.title = element_blank(),  
        #axis.title.x = element_text(face = "bold"), 
        #axis.title.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) 

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

print(volcanoPlot + # geom_text_repel(size = 3.5, max.overlaps = 10,
                      #              nudge_x = 0.3,
                       #             nudge_y = 0.3
# ) +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(-1, 1), col="black", linetype = 2) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = 2) #+
  # Change point color 
#   scale_color_manual(values=c("Downregulated" = "navyblue",
#                               "Not DE" = "gray", 
#                               "Upregulated" = "red"))
 )


tab_dep <- tab1 %>%
  # Get table of only depleted sgRNAs
  filter(score <= 0) %>%
  # Group table by Gene
  group_by(Gene) %>%
  # Create column that contain labels only for genes with two or more sgRNAs below the FDR threshold
  mutate(label = ifelse(sum(FDR < 5e-8) >= 2, "label_value", NA)) %>%
  ungroup()
# Make manhattan plot
ggplot(tab_dep, aes(x = Gene, y = -log10(FDR))) +
  geom_point(alpha = 0.5, position = position_dodge(width = 0),
             color = ifelse(is.na(tab_dep$label), "black", "red")) +
  xlab("Gene") +
  ylab(expression("-log10(FDR)")) +
  ggtitle("EZH2i vs Baseline") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    text = element_text(size = 12, face = "bold"),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  geom_hline(yintercept=-log10(5e-8), col="black", linetype = 2) +
  # scale_color_manual(values = c("black", "red")) +
  scale_y_continuous(#limits = c(min, max),
    # sec.axis = sec_axis(trans = ~ ., name = expression("CRISPR score/ MAGeCK"))
  ) +
  guides(color = guide_legend(title = "", override.aes = list(shape = 16)))
