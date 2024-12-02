#===
# CHANGE FOR EACH PROJECT
#===
workdir = "C:/Users/thephillipslab/Documents/Projects/CAR-T_DIPGVI_screen"
project_name = "Qplot_CART_DIPGVI_GD2_vs_CD19_invitro_twoMismatches_pooled"
path_to_results_table1 = file.path("./MAGeCK/CART_DIPGVI_invitro_twoMismatches_pooled/GD2_vs_baseline/GD2_vs_Baseline.sgrna_summary.txt")
path_to_results_table2 = file.path("./MAGeCK/CART_DIPGVI_invitro_twoMismatches_pooled/CD19_vs_baseline/CD19_vs_Baseline.sgrna_summary.txt")
path_to_select_guides = file.path("./Selected genes.xlsx")

# Plot settings
plot_title = "DIPGVI (in vitro)"
xlab = "LFC(GD2 vs Baseline)"
ylab = "LFC(CD19 vs Baseline)"
topguides = 50 # number of top depleted/enriched guides to label
#===

# Load the necessary packages if not already loaded
if(!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if(!"readxl" %in% installed.packages()) BiocManager::install("readxl")
if(!"ggplot2" %in% installed.packages()) BiocManager::install("ggplot2")
if(!"ggrepel" %in% installed.packages()) BiocManager::install("ggrepel")
if(!"tidyverse" %in% installed.packages()) BiocManager::install("tidyverse")
if(!"dplyr" %in% installed.packages()) BiocManager::install("dplyr")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readxl)

setwd(workdir)

output_dir = paste0(workdir, "/", "plots")
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

tab1 = read.delim(path_to_results_table1)
tab2 = read.delim(path_to_results_table2)

# tab1 %>% head()
tab1 %>% nrow()
# tab2 %>% head()
tab2 %>% nrow()

# Merge count tables
tab_merge = merge(tab1, tab2, by = "sgrna")

# tab_merge %>% head()
tab_merge %>% nrow()

# List of genes to check
# sel_genes = read_excel(path_to_select_guides)
# target_genes <- sel_genes$`Selected guides`
# 
# target_genes 
# 
# table(tab_merge$sgrna %in% target_genes)
# 
# gneg<-tab_merge
# # Add the label column and set its initial values to NA
# gneg$label <- NA
# 
# # Use %in% to check if each value in gneg$Gene.x is in the target_genes vector
# # and add an additional condition for LFC.x and LFC.y
# # gneg$label[gneg$sgrna %in% target_genes & 
# #              (gneg$LFC.x > 0.5 | gneg$LFC.x < -0.5 | gneg$LFC.y > 0.5 | gneg$LFC.y < -0.5) & 
# #              gneg$FDR.x < 0.1] <- gneg$sgrna[gneg$sgrna %in% target_genes & 
# #                                                (gneg$LFC.x > 0.5 | gneg$LFC.x < -0.5 | gneg$LFC.y > 0.5 | gneg$LFC.y < -0.5) & 
# #                                                gneg$FDR.x < 0.1]
# 
# gneg$label[gneg$sgrna %in% target_genes] <- gneg$sgrna[gneg$sgrna %in% target_genes]
# 
# gneg$label %>% unique()
# 
# label_data = subset(gneg, sgrna %in% target_genes)
# 
# # Create the scatter plot and label the data points for the specific genes
# p = ggplot(gneg) +
#         geom_point(aes(x = LFC.x, y = LFC.y, color = FDR.x),size=1.5) +
#         geom_vline(xintercept = 0,linetype="dashed") +
#         geom_hline(yintercept = 0,linetype="dashed") +
#         xlab(xlab) +
#         ylab(ylab) +
#         labs(title = plot_title) +
#         theme_minimal() +
#         geom_text_repel(
#                 data = label_data,
#                 aes(x = LFC.x, y = LFC.y, label = label, color = FDR.x),
#                 size = 6,  # Adjust the label size here
#                 nudge_y = -0.5,
#                 nudge_x = 2,
#                 max.overlaps = 8
#         ) +
#         scale_color_gradient(low = "darkblue", high = "red", 
#                              # limits = c(min(gneg$LFC.x), max(gneg$LFC.x)),
#                              limits = c(0, 0.1),
#                              breaks = c(0, 0.05, 0.1)) +
#         coord_equal(xlim = c(min(gneg$LFC.x), max(gneg$LFC.x)), 
#                     ylim = c(min(gneg$LFC.y), max(gneg$LFC.y))
#                     )+
#         theme(panel.grid.minor = element_blank(),  # Remove minor grid lines
#               plot.background = element_rect(color = "black", size = 1),
#               plot.title = element_text(hjust = 0.5))+
#         labs(color = "FDR")
# 
# 
# print(p)
# 
# 
# # Export as PDF
# pdf_name = paste0(output_dir, "/", project_name, ".pdf")
# pdf(file = pdf_name, 
#     width = 7, 
#     height = 7, 
#     title = project_name)
# 
# p
# 
# dev.off()







###
gneg<-tab_merge
label_data = gneg
# label_data = subset(gneg, LFC.y < 1 & LFC.y > -1)
label_data %>% nrow()
label_data = subset(label_data, LFC.x < -0.5)
label_data %>% nrow()
#label_data = subset(label_data, control_mean.x > 300)
label_data %>% nrow()

# Count occurrences of each Gene.x
gene_counts <- table(label_data$Gene.x)

# Identify genes occurring two times or more
genes_to_label <- names(gene_counts[gene_counts >= 2])

genes_to_label %>% length()


# label_data = subset(label_data, round(FDR.x, digits = 3) <= 0.1)
label_data = subset(label_data, !grepl("Control", label_data$Gene.x))
label_data %>% nrow()
label_data = subset(label_data, Gene.x %in% genes_to_label)
label_data %>% nrow()

# Order by lowest LFC
label_data <- label_data[order(label_data$LFC.x), ]
# Keep only top 20 depleted guides
# label_data <- label_data[1:topguides,]

# Add label
label_data$label <- label_data$sgrna
label_data = filter(label_data, !is.na(label_data$label))


# Number of genes represented in the top 20 most depleted guides
label_data$Gene.x %>% unique() %>% length()

# label_data = subset(label_data, Gene.x %in% genes_to_label)
maxlabels = min(topguides, nrow(label_data))
label_data2 = label_data[1:maxlabels,]
label_data2 = filter(label_data2, FDR.x <= 0.1)



# Create the scatter plot and label the data points for the specific genes
p = ggplot(gneg) +
  geom_point(aes(x = LFC.x, y = LFC.y), 
             color = "gray50",
             alpha = 0.5,
             size=1.5
             ) +
  geom_point(data = label_data2, aes(x = LFC.x, y = LFC.y, color = FDR.x), size=1.5) +
  geom_vline(xintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_hline(yintercept = 1,linetype="dashed", color = "blue") +
  geom_hline(yintercept = -1,linetype="dashed", color = "blue") +
  xlab(xlab) +
  ylab(ylab) +
  labs(title = plot_title) +
  theme_classic() +
  geom_text_repel(
    data = label_data2,
    aes(x = LFC.x, y = LFC.y, label = label, color = FDR.x),
    size = 3.5,  # Adjust the label size here
    nudge_y = ifelse(label_data2$LFC.y <= -2, -1.5, -1 * label_data2$LFC.y),
    # nudge_x = ifelse(label_data2$LFC.x <= -2, -0.5, 3),
    # nudge_y = 0.5 * label_data2$LFC.y,
    nudge_x = ifelse(label_data2$LFC.x/-4.7 <= 1, 0.25 * label_data2$LFC.x, -0.35 * label_data2$LFC.x),
    max.overlaps = 10
  ) +
  scale_color_gradient(low = "red", high = "darkblue", 
                       # limits = c(min(gneg$LFC.x), max(gneg$LFC.x)),
                       limits = c(0, 0.1),
                       breaks = c(0, 0.05, 0.1)) +
  coord_equal(xlim = c(min(gneg$LFC.x), max(gneg$LFC.x)), 
              ylim = c(min(gneg$LFC.y), max(gneg$LFC.y))
  )+
  theme(panel.grid.minor = element_blank(),  # Remove minor grid lines
        plot.title = element_text(hjust = 0.5))+
  labs(color = "FDR")


print(p)

head(label_data2)
head(label_data)

table_name = paste0(output_dir, "/label_data_controlMeanAbove300", project_name, ".txt")
output_table = label_data[,c("sgrna", "Gene.x", "FDR.x", "LFC.x", "FDR.y", "LFC.y")]
colnames(output_table)[2] = "Gene"
write.table(output_table, file = table_name, quote = FALSE,sep = "\t", row.names = FALSE)

# Export as PDF
pdf_name = paste0(output_dir, "/", project_name, ".pdf")
pdf(file = pdf_name,
    width = 7,
    height = 7,
    title = project_name)

p

dev.off()

# # Export as PNG
# png_name = paste0(output_dir, "/", project_name, ".png")
# png(file = png_name, 
#     width = 480, 
#     height = 480, 
#     type = "windows",
#     res = 72 #ppi
# )
# 
# p
# 
# dev.off()




# # Scatter plot for in vivo baseline vs in vitro baseline comparison
# p = ggplot(gneg) +
#   geom_point(aes(x = LFC.x, y = LFC.y), size=1.5, color = "gray50") +
#   geom_point(data = label_data, aes(x = LFC.x, y = LFC.y, color = score.x), size=1.5) +
#   geom_vline(xintercept = 0,linetype="dashed") +
#   geom_hline(yintercept = 0,linetype="dashed") +
#   geom_hline(yintercept = 1,linetype="dashed", color = "blue") +
#   geom_hline(yintercept = -1,linetype="dashed", color = "blue") +
#   geom_vline(xintercept = 1,linetype="dashed", color = "blue") +
#   geom_vline(xintercept = -1,linetype="dashed", color = "blue") +
#   xlab(xlab) +
#   ylab(ylab) +
#   labs(title = plot_title) +
#   theme_minimal() +
#   geom_text_repel(
#     data = label_data,
#     aes(x = LFC.x, y = LFC.y, label = label, color = score.x),
#     size = 4,  # Adjust the label size here
#     nudge_y = ifelse(gneg$LFC.y >= 0, 1, -2),
#     nudge_x = ifelse(gneg$LFC.x >= -2, -0.5, 1),
#     max.overlaps = 5
#   ) +
#   scale_color_gradient(high = "red", low = "darkblue", 
#                        # limits = c(min(gneg$LFC.x), max(gneg$LFC.x)),
#                        limits = c(min(gneg$score.x), max(gneg$score.x)),
#                        # limits = c(0, 0.1),
#                        #breaks = c(0, 0.05, 0.1)
#   ) +
#   coord_equal(xlim = c(min(gneg$LFC.x), max(gneg$LFC.x)), 
#               ylim = c(min(gneg$LFC.y), max(gneg$LFC.y))
#   )+
#   theme(panel.grid.minor = element_blank(),  # Remove minor grid lines
#         plot.background = element_rect(color = "black", size = 1),
#         plot.title = element_text(hjust = 0.5))+
#   labs(color = "CRISPR\nscore")
# 
# 
# 
# xmin = min(-gneg$LFC.y, -gneg$LFC.x)

# # Scatterplot LFC vs LFC
# p = ggplot(gneg) +
#   geom_point(aes(x = -LFC.x, y = -LFC.y), 
#              color = "gray50",
#              alpha = 0.5,
#              size=1.5
#   ) +
#   geom_point(data = label_data2, aes(x = LFC.x, y = LFC.y), 
#              color = "red",
#              size=1.5) +
#   # geom_vline(xintercept = 0,linetype="dashed") +
#   # geom_hline(yintercept = 0,linetype="dashed") +
#   # geom_hline(yintercept = 1,linetype="dashed", color = "blue") +
#   # geom_hline(yintercept = -1,linetype="dashed", color = "blue") +
#   geom_abline(slope = -gneg$LFC.y/(-gneg$LFC.x), intercept = min(-gneg$LFC.y), linetype = "dashed", color = "black") +
#   xlab(paste0("-",xlab)) +
#   ylab(paste0("-",ylab)) +
#   labs(title = plot_title) +
#   theme_minimal() +
#   geom_text_repel(
#     data = label_data2,
#     aes(x = LFC.x, y = LFC.y, label = label, color = FDR.x),
#     size = 3.5,  # Adjust the label size here
#     # nudge_y = ifelse(label_data2$LFC.y >= 0, 2, -3),
#     # nudge_x = ifelse(label_data2$LFC.x <= -2, -0.5, 3),
#     nudge_y = 1 * label_data2$LFC.y,
#     nudge_x = ifelse(label_data2$LFC.x <= -1, -0.5, 2),
#     max.overlaps = 10
#   ) +
#   scale_color_gradient(low = "red", high = "darkblue", 
#                        # limits = c(min(gneg$LFC.x), max(gneg$LFC.x)),
#                        limits = c(0, 0.1),
#                        breaks = c(0, 0.05, 0.1)) +
#   coord_equal(xlim = c(min(gneg$LFC.x), max(gneg$LFC.x)), 
#               ylim = c(min(gneg$LFC.y), max(gneg$LFC.y))
#   )+
#   theme(panel.grid.minor = element_blank(),  # Remove minor grid lines
#         #plot.background = element_rect(color = "black", size = 1),
#         plot.title = element_text(hjust = 0.5))+
#   labs(color = "FDR")





limits = c(min(range(label_data$LFC.x)), max(range(range(label_data$LFC.y)))) %>% rev()

p = ggplot(data = label_data, 
           aes(x = -LFC.y, 
               y = -LFC.x,
               #color = colors
           )
) + 
  geom_point(data = label_data %>% filter(-LFC.y - (-LFC.x) <= -1), 
             size = 2, 
             shape=21, fill = "red") +
  # geom_point(data = label_data %>% filter(-LFC.y - (-LFC.x) <= -1 ),
  #            size = 2.5,
  #            shape=21,
  #            fill = "black") +
  geom_point(data = gneg %>% filter(-LFC.y - (-LFC.x) > -1), 
             size = 2, 
             color = "gray",
             alpha = 0.6) +
  theme_bw() +
  geom_text_repel(
    data = label_data %>% filter(-LFC.y - (-LFC.x) <= -1),
    aes(x = -LFC.y, y = -LFC.x, label = label),
    color = "black",
    size = 3.5,  # Adjust the label size here
    nudge_y = ifelse(-1 * label_data2$LFC.x <= 3, -1, 0.5),
    nudge_x = ifelse(-1 * label_data2$LFC.y >= 2, 1.5, -1),
    # nudge_y = 0.5 * label_data2$LFC.y,
    # nudge_x = ifelse(-1 * label_data2$LFC.y <= 0, -2, 3),
    # nudge_y = -1,
    # nudge_x = 2,
    max.overlaps = 8
  ) +
  # ylab(expression(bold(bolditalic(-Log[2]) ~ "FC / EZH2i vs baseline"))) +
  # xlab(expression(bold(bolditalic(-Log[2]) ~ "FC / DMSO vs baseline"))) +
  xlab("Depletion at CD19 vs Baseline (-1 x LFC)") +
  ylab("Depletion at GD2 vs Baseline (-1 x LFC)") +
  labs(title = plot_title) +
  scale_x_continuous(limits = -limits) +
  scale_y_continuous(limits = -limits) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 13),
        panel.grid = element_blank()) +
  # Add reference lines
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "blue") 


p


# Export as PDF
pdf_name = paste0(output_dir, "/LFC_vs_LFC_", project_name, ".pdf")
pdf(file = pdf_name,
    width = 7,
    height = 7,
    title = project_name)

p

dev.off()
