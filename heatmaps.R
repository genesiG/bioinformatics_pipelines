if(!require(profileplyr)){
  BiocManager::install("profileplyr")
}
if(!require(Cairo)){
  BiocManager::install("Cairo")
}
if(!require(magick)){
  BiocManager::install("magick")
}


library(profileplyr)
library(Cairo)
library(ggplot2)
library(dplyr)
library(tidytext)
library(ggrepel)
set.seed(123)

## Set working directory
setwd("C:/Users/thephillipslab/Documents/Projects/Enhancers_DIPG/")

## Load peaks as GRanges objects
dipg13_peaks = rtracklayer::import(con = "./macs3_hg19/dipg13_h3k27ac_SRR8551305_R1.q005.fe50.gap2k.min200_peaks.narrowPeak", 
                                   format = "narrowPeak")
iPSC_peaks = rtracklayer::import(con = "./macs3_hg19/iPSC-2242-1_h3k27ac_SRR9307340_R1.q005.fe50.gap2k.min200_peaks.narrowPeak", 
                                   format = "narrowPeak")
liver_peaks = rtracklayer::import(con = "./macs3_hg19/liver_NL1_h3k27ac_SRR6880492_R1.q005.fe50.gap2k.min200_peaks.narrowPeak", 
                                   format = "narrowPeak")
lung_peaks = rtracklayer::import(con = "./macs3_hg19/lung_HBEC_h3k27ac_SRR24049132_R1.q005.fe50.gap2k.min200_peaks.narrowPeak", 
                                   format = "narrowPeak")
neuron_peaks = rtracklayer::import(con = "./macs3_hg19/motor_neuron_h3k27ac_SRR22511762_R1.q005.fe50.gap2k.min200_peaks.narrowPeak", 
                                   format = "narrowPeak")

dipg_exclusive_peaks = rtracklayer::import(con = "./macs3_hg19_ENCODE/GB_exclusive_peaks.narrowPeak", 
                                                       format = "narrowPeak")

non_tumor_peaks = rtracklayer::import(con = "./macs3_hg19/non_tumor_peaks.narrowPeak", 
                                                     format = "narrowPeak")


## Get the union set of all non-dipg peaks
non_tumor_peaks = GenomicRanges::union(x = iPSC_peaks, 
                                       y = liver_peaks) %>% 
  GenomicRanges::union(y = lung_peaks) %>% 
  GenomicRanges::union(y = neuron_peaks)

## Take only the dipg peaks that DO NOT overlap the non-dipg peaks
dipg_exclusive_peaks = subsetByOverlaps(dipg13_peaks, 
                                        non_tumor_peaks, 
                                        invert = TRUE)


## Import signal quantification from deepTools as profileplyr object
dipg13_Monje <- import_deepToolsMat("./deepTools_hg19_ENCODE/matrix/GB_exclusive_vs_liver_lung_neuron.gz")

## Check params from computeMatrix call
dipg13_Monje@sampleData

## Check regions
dipg13_Monje@rowRanges

#===
# Unbiased clustering of enrichment signal
#===
kmeans_dipg13_Monje <- clusterRanges(dipg13_Monje, 
                                      fun = rowMeans,  
                                      # kmeans_k = 2
                                     )

mcols(dipg13_Monje)$dpGroup = "Glioblastoma-exclusive peaks (159)"
generateEnrichedHeatmap(dipg13_Monje, 
                        include_group_annotation = F, 
                        group_anno_color = "red", 
                        group_anno_width = 0,
                        # For sample names at the top of each heatmap
                        matrices_column_title_gp = grid::gpar(fontsize = 12, 
                                                              fontface = "bold", 
                                                              lwd = 2), 
                        raster_device = "tiff",
                        matrices_color = c("navyblue", "orange", "red"),
                        sample_names = c("GB01", "GB02", "GB04", "GB07", 
                                         "Liver", "Lung", "Brain"),
                        matrices_axis_name = c("-2kb","center","+2kb"))

#########

## Import signal quantification from deepTools as profileplyr object
a <- import_deepToolsMat("./hg19_aligned/deepTools/CBX2_DB_vsK27K4.gz")

# Check sample names
a@sampleData %>% rownames()

#===
# Clustering based on pre-defined regions
#===
mcols(a)$Group = "Retained CBX2"
mcols(a)$Group[grep("increased", mcols(a)$names)] = "Increased CBX2"
mcols(a)$Group[grep("decreased", mcols(a)$names)] = "Decreased CBX2"

## Heatmap is clustered based on the `dpGroup` column
mcols(a)$dpGroup = mcols(a)$Group
mcols(a)

a[,,c(1:3)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "purple"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)",  "CBX2\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))


a[,,c(4:6)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "goldenrod"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)", "H3K27me3\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(7:9)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "blue"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)", "H3K27me2\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(10:12)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "red"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K4me3\n(DMSO)","H3K4me3\n(1uM EZH2i)", "H3K4me3\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))



#===
# USING EPICYPHER'S CUT AND RUN DATA
#===

setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/")

## Import signal quantification from deepTools as profileplyr object
a <- import_deepToolsMat("./deepTools/CBX2_DB_vsK27K4.gz")

# Check sample names
a@sampleData %>% rownames()

#===
# Clustering based on pre-defined regions
#===
mcols(a)$Group = "Retained\nCBX2"
mcols(a)$Group[grep("increased", mcols(a)$names)] = "Increased\nCBX2"
mcols(a)$Group[grep("decreased", mcols(a)$names)] = "Decreased\nCBX2"

## Heatmap is clustered based on the `dpGroup` column
mcols(a)$dpGroup = mcols(a)$Group
mcols(a)

a[,,c(1:3)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "purple"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)",  "CBX2\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))


a[,,c(4:6)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "goldenrod"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)", "H3K27me3\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(7:9)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "blue"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)", "H3K27me2\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(10:12)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "red"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H3K4me3\n(DMSO)","H3K4me3\n(1uM EZH2i)", "H3K4me3\n(5uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))






























a[,,c(9:10)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "darkgreen"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("EZH2\n(DMSO)", "EZH2\n(1uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(11:12)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "pink"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("H2A119ub\n(DMSO)", "H2A119ub\n(1uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(13:14)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "darkgray"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("RFN2\n(DMSO)", "RFN2\n(1uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a[,,c(15:16)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "black"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("HP1beta\n(DMSO)", "HP1beta\n(1uM EZH2i)"),
                          matrices_axis_name = c("-5kb","center","+5kb"))



















