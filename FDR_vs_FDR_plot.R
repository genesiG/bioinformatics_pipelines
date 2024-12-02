#=================#
# FDR vs FDR PLOT #
#=================#

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






