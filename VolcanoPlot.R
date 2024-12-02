setwd("C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2")
resultsDF$Label <- NA
#resultsDF$Label[resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$isDE != "Not DE"]
resultsDF$Label[resultsDF$Peak == "Retained CBX2 binding" &
                  resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$Peak == "Retained CBX2 binding" & resultsDF$isDE != "Not DE"]
resultsDF$Label[resultsDF$Peak == "Increased CBX2 binding" &
                  resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$Peak == "Increased CBX2 binding" & resultsDF$isDE != "Not DE"]



resultsDF %>% filter(row == "CDKN2A")
resultsDF$Label[resultsDF$row == "CDKN2A"] <- "CDKN2A"

# print(x = resultsDF[order(resultsDF$pvalue), ],
#       max = 100,
#       row.names = F)

# Set plot title for each comparison
comparison = "CBX2 KO + EZH2i vs EZH2i"

# Make volcano plot
volcanoPlot <- ggplot(data=resultsDF, 
                      aes(x=log2FoldChange, 
                          y=-log10(padj),
                          color = isDE, 
                          shape = Peak, 
                          label =  Label)) + 
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values = c(1,16,4,17)) +
  ylab(expression(bold(bolditalic(-Log[10]) ~ "FDR"))) +
  xlab(expression(bold(bolditalic(Log[2]) ~ "FC"))) +
  ggtitle(comparison) + 
  scale_x_continuous(n.breaks = 7, 
                     #limits = c(-8.2,8.2)
                     ) +
  theme(legend.title = element_blank(),  
        #axis.title.x = element_text(face = "bold"), 
        #axis.title.y = element_text(face = "bold"), 
        panel.grid.minor = element_blank()) 

# We can organize the labels nicely using the "ggrepel" package and the 
# geom_text_repel() function

print(volcanoPlot + geom_text_repel(size = 3.5, max.overlaps = 10,
                                    nudge_x = 0.3,
                                    nudge_y = 0.3
                                    ) +
        # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
        # geom_vline(xintercept=c(-FC.cutoff, FC.cutoff), col="magenta", linetype = 2) +
        # geom_hline(yintercept=-log10(alpha), col="black", linetype = 2) +
        # Change point color 
        scale_color_manual(values=c("Downregulated" = "navyblue",
                                    "Not DE" = "gray", 
                                    "Upregulated" = "red")) + guides(color = F)
)



#===
# LFC VS LFC PLOT
#===

results1 = resultsDF
results2 = resultsDF

results1 %>% head()
results2 %>% head()


library(dplyr)

merged <- merge(select(results1, -baseMean, -stat, -pvalue),
                select(results2, -baseMean, -stat, -pvalue),
                by = c("row", "Peak"),
                suffixes = c(".1", ".2"))
results_merge = merged
results_merge %>% head()

# results_merge$colors = "Other complexes"
# results_merge$colors[grepl("PRC",results_merge$Protein.complex)] <- "PRC"

# Calculate the difference between the y-value and the reference line (slope = 1, intercept = 1)
results_merge$above_reference_line <- -results_merge$log2FoldChange.2 + (-results_merge$log2FoldChange.1)

# Get table of depleted guides
results_merge_dep = filter(results_merge, log2FoldChange.2 > 0)
results_merge_dep = filter(results_merge, log2FoldChange.1 > 0)
results_merge_dep = filter(results_merge, isDE.1 != "Not DE" )

limits = c(0, max(range(range(results_merge_dep$log2FoldChange.1)))) 

results_merge_dep$Label <- NA
#resultsDF$Label[resultsDF$isDE != "Not DE"] <- resultsDF$row[resultsDF$isDE != "Not DE"]
# results_merge_dep$Label[results_merge_dep$Peak == "Retained CBX2 binding" &
#                   results_merge_dep$isDE != "Not DE"] <- results_merge_dep$row[results_merge_dep$Peak == "Retained CBX2 binding" & results_merge_dep$isDE != "Not DE"]
results_merge_dep$Label[results_merge_dep$Peak == "Increased CBX2 binding"] <- results_merge_dep$row[results_merge_dep$Peak == "Increased CBX2 binding"]
results_merge_dep$Label[results_merge_dep$log2FoldChange.2 - (results_merge_dep$log2FoldChange.1) > -0.5] <- NA

results_merge_dep$Peak = factor(x = results_merge_dep$Peak, 
                                levels = c("Increased CBX2 binding", 
                                           "Decreased CBX2 binding", 
                                           "Retained CBX2 binding", 
                                           "No CBX2 occupancy"))

results_merge_dep %>% filter(row == "CDKN2A")
# results_merge_dep$Label[results_merge_dep$row == "CDKN2A"] <- "CDKN2A"
# results_merge_dep %>% head()

p = ggplot(data = results_merge_dep, 
           aes(x = log2FoldChange.2, 
               y = log2FoldChange.1,
               )
) + 
  geom_point(data = results_merge_dep %>% filter(log2FoldChange.2 - (log2FoldChange.1) <= -0.5), 
             size = 2, 
             shape = 21,
             aes(fill = Peak, color = Peak)) +
  geom_point(data = results_merge_dep %>% filter(log2FoldChange.2 - (log2FoldChange.1) <= -0.5 & Peak == "Increased CBX2 binding"), 
             size = 3, 
             shape = 21,
             color = "black", fill = "red") +
  geom_point(data = results_merge_dep %>% filter(log2FoldChange.2 - (log2FoldChange.1) > -0.5),
             size = 2,
             color = "gray",
             alpha = 0.6) +
  theme_bw() +
  ylab(expression(bold(bolditalic(Log[2]) ~ "FC / EZH2i + CBX2 KO vs DMSO"))) +
  xlab(expression(bold(bolditalic(Log[2]) ~ "FC / EZH2i vs DMSO"))) +
  scale_x_continuous(limits = limits) +
  scale_y_continuous(limits = limits) +
  theme(axis.title = element_text(size = 13.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 13.5),
        panel.grid = element_blank()) +
  # Add reference lines
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_abline(intercept = 0.5, slope = 1, linetype = "dashed", color = "blue") +
  #geom_abline(intercept = -1, slope = 1, linetype = "dashed", color = "blue") +
  geom_text_repel(aes(label = Label), 
                  color = "red", 
                  xlim = c(2.5, 3.7),
                  nudge_y = -0.2, 
                  nudge_x = 0.5, 
                  max.overlaps = 1
                  ) +
  # Change point color
  scale_fill_manual(values = c("red", "white","black", "white")) +
  scale_color_manual(values = c("black", "red","black", "black")) 


p

setwd("C:/Users/thephillipslab/Documents/Projects/EZH2i & PRC1 RNAseq/DESeq2")

pdf("./plots/LFC_vs_LFC_RNAseq.pdf", width = 7.65, height = 5)
p
dev.off()












## BAR CHART OF DIFFERENTIALLY EXPRESSED GENES

# Define TOP20 DE genes in the dual treatment condition
dualTreatTOP20 <-resultsDF %>% filter(isDE != "Not DE")
dualTreatTOP20 <- dualTreatTOP20[order(dualTreatTOP20$padj), ]
dualTreatTOP20 = dualTreatTOP20[1:20,]
dualTreatTOP20

dualTreatTOP20 <- resultsKO_EZH2i[order(resultsKO_EZH2i$padj), ]
dualTreatTOP20 = dualTreatTOP20[1:20,]
dualTreatTOP20$Condition = "CBX2 KO + EZH2i"

# Create pooled DF with data from dual treatment
pooledResults = dualTreatTOP20
pooledResults %>% nrow()

# Reorder by log2FC
pooledResults = pooledResults[order(pooledResults$log2FoldChange),]

#
#
# Add other conditions to the DF
pooledResults = rbind(pooledResults, resultsDF)
pooledResults %>% nrow()

pooledResults$Condition <- factor(x=pooledResults$Condition, 
                                  levels = c("CBX2 KO",
                                             "EZH2i", 
                                             "CBX2 KO + EZH2i"))



# Filter to TOP20 DE genes
filteredResults = filter(pooledResults, row %in% dualTreatTOP20$row)
filteredResults %>% nrow()
filteredResults = filteredResults[order(filteredResults$padj, decreasing = T),]


# Add factor label for y axis
filteredResults$genes = factor(filteredResults$row, 
                             levels = unique(filteredResults$row))


# Plot bar graph
ggplot(filteredResults, aes(x = log2FoldChange, y = genes, fill = Condition)) +
  #geom_col(#stat = "identity",
           #position = "dodge") +
  geom_col(colour = "black", position = "dodge") +
  scale_fill_manual(values = c("black","white","gray")) +
  theme_classic() + ylab("") + 
  #theme(axis.text.y = element_text(face = "bold")) +
  xlab(expression(bolditalic(Log[2])~bold("(Fold-Change)"))) + 
  ggtitle("Top 20 differentially expressed genes", 
          subtitle = "(CBX2 KO + EZH2i vs DMSO)") #+ 
  geom_vline(xintercept = 0, 
             color = "red", 
             lwd = 1, lty = "dashed")





#######
print(x = results_Hist[order(results_Hist$pvalue), ],
      max = 100,
      row.names = F)
library(ggplot2)



