if(!require(profileplyr)){
  BiocManager::install("profileplyr")
}

library(profileplyr)
library(ggplot2)
library(dplyr)
library(tidytext)
library(ggrepel)

## Import signal quantification from deepTools as profileplyr object
setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/")

CBX2_5kb_mean <- import_deepToolsMat("./deepTools/CBX2_5kb_mean.gz")

# The profileplyr object is a form of the RangedSummarizedExperiment class, 
# which allows us to both store all of the relevant information that is imported
# from soGGi or deepTools, and have flexibility in how we manipulate this object
# in downstream analysis

library(SummarizedExperiment)

assays(CBX2_5kb_mean)

rowRanges(CBX2_5kb_mean)[1:3]

# The information associated with each sample is stored as a Dataframe within 
# with the parameters of the profileplyr object, and is accessed with sampleData().

sampleData(CBX2_5kb_mean)

# The ‘rowGroupsInUse’ element for the params slot indicates which column of the
# range metadata (mcols() or rowRanges()) that will be used for grouping in the 
# final output if that object were used for visualization. For a profileplyr 
# object created from a deepTools computeMatrix output, the inherited groups 
# correspond to the BED files over which the signal was quantified, and these 
# groups are contained in the ‘dpGroup’ column.

params(CBX2_5kb_mean)$rowGroupsInUse

# The ‘mcolToOrderBy’ slot of the profileplyr object parameters indicates which 
# column of the range metadata will be used for ordering the ranges as they are 
# exported to either deepTools or to EnrichedHeatmap (with the 
# generateEnrichedHeatmap() function within profileplyr). 

# This can be set using the orderBy() function, which requires a profileplyr 
# object and a character string matching a column name of the range metadata as 
# arguments. If groupBy is never used on an object, the default will be to order
# by the mean signal of each range (within each group). 

# In addition, until groupBy() has been used on a profileplyr object, the 
# ‘mcolToOrderBy’ slot of the parameters will be NULL and will not be seen 
# with params(proplyrObject).

params(CBX2_5kb_mean)$mcolToOrderBy # defaults to NULL

# To specify a column to reorder the ranges when they are exported:
CBX2_5kb_mean <- orderBy(CBX2_5kb_mean, "score")

params(CBX2_5kb_mean)$mcolToOrderBy # now will return "score"

## Subsetting the profileplyr object

# The profileplyr object can be subset either by sample, or by rows and columns
# of the matrix for each sample. This is done using the ‘[ ]’ brackets, with the
# first position being assay matrix rows, the second position being assay matrix
# columns, and the third position being the entire matrix for each sample 
# (in addition to the rest of the parameters and range data). 

CBX2_5kb_mean[1:10,1:10] # returns the first 10 rows and first 10 columns of the matrix

CBX2_5kb_mean[ , , 1:2] # returns the entire matrix for the first two samples

## Changing sample names

rownames(sampleData(CBX2_5kb_mean)) # returns sample names

rownames(sampleData(CBX2_5kb_mean)) <- c("CBX2/DMSO", "GST-CBX2/DMSO",
                                         "CBX2/1uM EZH2i","GST-CBX2/1uM EZH2i",
                                         "CBX2/5uM EZH2i", "GST-CBX2/5uM EZH2i") # reassigns samples' names

rownames(sampleData(CBX2_vsK27K4_5kb_median)) <- c("CBX2/DMSO", "GST-CBX2/DMSO",
                                                   "H3K27me3/DMSO", "H3K4me3/DMSO")

# Add width of reagions 

CBX2_ranges = rowRanges(CBX2_5kb_mean) %>% ranges() %>% as.data.frame()

mcols(CBX2_5kb_mean)$widths = CBX2_ranges$width
mcols(CBX2_5kb_mean)


#----------
# Export/Conversion of profileplyr object for heatmap visualization of ranges
#----------

# To generate a heatmap within R directly from the profileplyr object, the 
# generateEnrichedHeatmap() function should be used. This function takes a 
# profileplyr object and produces a heatmap using the EnrichedHeatmap package.
# It allows for easy export from the profileplyr object and simple inclusion of
# the metadata as heatmap annotations. 

heatmap_CBX2_5kb_mean <- generateEnrichedHeatmap(CBX2_5kb_mean)


# This function generates a multipanel heatmap with a variety of arguments that
# have been tailored to visualizing the profileplyr object. 
class(heatmap_CBX2_5kb_mean)

# The convertToEnrichedHeatmapMat() function within profileplyr takes a 
# profileplyr object as the only required input, and then converts the matrices 
# contained within the assays slot to a list of ‘normalizedMatrix’ class objects
# that can be used in the EnrichedHeatmap() function. See the EnrichedHeatmap 
# vignette for detailed examples for how heatmaps can be concatenated together 
# to visualize all of the data in one figure.

EH_mat <- convertToEnrichedHeatmapMat(CBX2_5kb_mean)
EH_mat[[1]]

#----------
# Summarize signal for ggplot or heatmap visualization
#----------

long_CBX2_5kb_mean <- profileplyr::summarize(CBX2_5kb_mean, 
                                             fun = rowMeans, 
                                             output = "long") 
long_CBX2_5kb_mean %>% head()


ggplot(long_CBX2_5kb_mean, 
       aes(x = Sample, 
           y = log(Signal))) + 
  geom_boxplot() + 
  theme_bw()

#----------
# Annotating of genomic ranges with clusters, genomic regions, and genes
#----------

# clusterRanges() takes a profileplyr object as its first argument as well as a 
# function to summarize the signal in each range (similar to the summarize() 
# function above). The pheatmap package is then used to cluster the ranges, and 
# the type of clustering depends on whether the user inputs a value for 
# ‘kmeans_k’ (for kmeans clustering) or ‘cutree_rows’ (for hierarchical 
# clustering using hclust). An integer value entered for either of these 
# arguments will specify the number of clusters used for each method. 

# If both ‘kmeans_k’ or ‘cutree_rows’ are left as NULL (default), then a heatmap
# will be printed with hierarchical clustering, but no distinct clusters 
# defined, and no profileplyr object will be returned. This might be helpful as 
# an initial and quick look at the ranges or as a means to determine the number
# of clusters to try.

set.seed(0)
clusterRanges(CBX2_5kb_mean, 
              fun = rowMeans, kmeans_k = 2, silent = F)

#----------
# Generate group-annotated heatmap in R directly with generateEnrichedHeatmap()
#----------

# After clustering the profileplyr object can be passed directly as an argument 
# into the generateEnrichedHeatmap() function, and by default the heatmap will 
# be grouped and annotated by these clusters, which were automatically set as 
# the ‘rowGroupsInUse’ in the clusterRanges() function. 

# Assuming that the ‘include_group_annotation’ argument of this function is set 
# to the default value of TRUE, whichever metadata column is set to the 
# ‘rowGroupsInUse’ slot will be used for the grouping and annotation seen below 
# to the left of the heatmap. If there are no range groups, then the user can 
# set this argument to be FALSE, and those color annotations will be absent.

# Further, the maximum value for the y-axis in the line plots at the top of each
# heatmap are also automatically set based on the highest mean range signal from
# all of the samples. This can be set manually with the ‘ylim’ argument, or if 
# ylim = NULL, the maximum will be inferred (i.e. be different) for each 
# individual heatmap.

kmeans_CBX2_5kb_mean <- clusterRanges(CBX2_5kb_mean, 
                                      fun = rowMeans,  
                                      kmeans_k = 2)

generateEnrichedHeatmap(kmeans_CBX2_5kb_mean, )

hclust_CBX2_5kb_mean <- clusterRanges(CBX2_5kb_mean, 
                                      fun = rowMeans,  
                                      cutree_rows = 2, silent = F)

generateEnrichedHeatmap(hclust_CBX2_5kb_mean)

CBX2_db_vK27K4 = import_deepToolsMat("./deepTools/CBX2_DB_vsK27K4_5kb_median.gz")

mcols(CBX2_db_vK27K4)$Group = "Retained"
mcols(CBX2_db_vK27K4)$Group[grep("increased", mcols(CBX2_db_vK27K4)$names)] = "Increased"
mcols(CBX2_db_vK27K4)$Group[grep("decreased", mcols(CBX2_db_vK27K4)$names)] = "Decreased"

mcols(CBX2_db_vK27K4)$dpGroup = mcols(CBX2_db_vK27K4)$Group
mcols(CBX2_db_vK27K4)


a<-clusterRanges(CBX2_db_vK27K4[,,c(1:3,7:9)], cutree_rows = 2 )
#pdf("./CBX2_DB_vsK27me3.pdf")
generateEnrichedHeatmap(a, matrices_color = c("white",
                                                            #"red",
                                                            "purple"),
                        ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                        all_color_scales_equal = F,
                        sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)","CBX2\n(5uM EZH2i)", 
                                         #"GST-CBX2\n(DMSO)","GST-CBX2\n(1uM EZH2i)","GST-CBX2\n(5uM EZH2i)"#,
                                         "H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)","H3K27me3\n(5uM EZH2i)"#,
                                         #"H3K4me3\n(DMSO)","H3K4me3\n(1uM EZH2i)","H3K4me3\n(5uM EZH2i)"#,
                                         #"H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                                         ),
                        matrices_axis_name = c("-5kb","center","+5kb")) %>% pdf(file = "./CBX2_DB_vsK27me3.pdf")
#dev.off()


CBX2_5kb_median = import_deepToolsMat("./deepTools/CBX2_vsGST_5kb_median.gz") #%>%

#########
a <- import_deepToolsMat("./deepTools/CBX2_DB_vPRC_5kb_mean.gz") #%>%

mcols(a)$Group = "Retained CBX2"
mcols(a)$Group[grep("increased", mcols(a)$names)] = "Increased CBX2"
mcols(a)$Group[grep("decreased", mcols(a)$names)] = "Decreased CBX2"
mcols(a)$dpGroup = mcols(a)$Group
mcols(a)

a[,,c(7:8)] %>% 
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "blue"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c(#"CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)"#, #"CBX2\n(5uM EZH2i)",
                                           #"CBX2-CD-GST\n(DMSO)", "CBX2-CD-GST\n(1uM EZH2i)", #"CBX2-CD-GST\n(5uM EZH2i)"
                                           #"H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)"#,
                                           #"EZH2\n(DMSO)", "EZH2\n(1uM EZH2i)", #"EZH2\n(5uM EZH2i)"#,
                                           #"H2A119ub\n(DMSO)", "H2A119ub\n(1uM EZH2i)",
                                           #"RFN2\n(DMSO)", "RFN2\n(1uM EZH2i)",
                                           #"HP1beta\n(DMSO)", "HP1beta\n(1uM EZH2i)"
                                           "H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                          ),
                          matrices_axis_name = c("-5kb","center","+5kb"))

a_EH <- convertToEnrichedHeatmapMat(a[,,c(1:2,5:8)], 
                                    sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)", #"CBX2\n(5uM EZH2i)",
                                                     #"CBX2-CD-GST\n(DMSO)", "CBX2-CD-GST\n(1uM EZH2i)", #"CBX2-CD-GST\n(5uM EZH2i)"
                                                     "H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)",
                                                     #"EZH2\n(DMSO)", "EZH2\n(1uM EZH2i)", #"EZH2\n(5uM EZH2i)"#,
                                                     #"H2A119ub\n(DMSO)", "H2A119ub\n(1uM EZH2i)",
                                                     #"RFN2\n(DMSO)", "RFN2\n(1uM EZH2i)",
                                                     #"HP1beta\n(DMSO)", "HP1beta\n(1uM EZH2i)"
                                                     "H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                                    ))
  
EnrichedHeatmap::EnrichedHeatmap(a_EH[], 
                                 col = c("white","purple"),
                                 axis_name = c("-5kb","center","+5kb")) 

EnrichedHeatmap::EnrichedHeatmap(a_EH[[3:4]], 
                                 col = c("white","goldenrod"),
                                 axis_name = c("-5kb","center","+5kb")) 

EnrichedHeatmap::EnrichedHeatmap(a_EH[[5:6]], 
                                 col = c("white","blue"),
                                 axis_name = c("-5kb","center","+5kb"))
 
  
a %>%
clusterRanges(#CBX2_vsK27K4_5kb_median[,,c(1,3,4)],
              fun = rowMeans,
              scaleRows = T,
              #cutree_rows = 2,
              silent = T,
              #cluster_method = "single",
              #clustering_distance_rows = "correlation",

              cluster_sample_subset = c(1:6),
              kmeans_k = 2,
              ) %>%
  generateEnrichedHeatmap(matrices_color = c("white",
                                             #"red",
                                             "purple"), 
                          #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                          sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)", #"CBX2\n(5uM EZH2i)",
                                           "CBX2-CD-GST\n(DMSO)", "CBX2-CD-GST\n(1uM EZH2i)", #"CBX2-CD-GST\n(5uM EZH2i)"
                                           "H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)",
                                           "EZH2\n(DMSO)", "EZH2\n(1uM EZH2i)", #"EZH2\n(5uM EZH2i)"#,
                                           "H2A119ub\n(DMSO)", "H2A119ub\n(1uM EZH2i)",
                                           "RFN2\n(DMSO)", "RFN2\n(1uM EZH2i)",
                                           "HP1beta\n(DMSO)", "HP1beta\n(1uM EZH2i)"
                                           #"H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                          ),
                          matrices_axis_name = c("-5kb","center","+5kb"))


rowRanges(CBX2_vsK27K4_5kb_median)
mcols(CBX2_vsK27K4_5kb_median)
rowData(CBX2_vsK27K4_5kb_median)

mcols(CBX2_5kb_median)$dpGroup = mcols(CBX2_vsK27K4_5kb_median)$cluster

CBX2_5kb_median %>% rowRanges()
profileplyr::groupBy(CBX2_5kb_median, group = "cluster")

generateEnrichedHeatmap(CBX2_5kb_median, matrices_color = c("white",
                                           #"red",
                                           "purple"), 
                        #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                        sample_names = c("CBX2/DMSO", "CBX2/1uM EZH2i","CBX2/5uM EZH2i", 
                                         "GST-CBX2/DMSO","GST-CBX2/1uM EZH2i","GST-CBX2/5uM EZH2i"),
                        matrices_axis_name = c("-5kb","center","+5kb"))


CBX2_DB_vK27me2 = import_deepToolsMat("./deepTools/CBX2_DB_vsPRC_5kb_mean.gz") 

mcols(CBX2_DB_vK27me2)$Group = NA
mcols(CBX2_DB_vK27me2)$Group[grep("increased", mcols(CBX2_DB_vK27me2)$names)] = "Increased CBX2"
mcols(CBX2_DB_vK27me2)$Group[grep("decreased", mcols(CBX2_DB_vK27me2)$names)] = "Decreased CBX2"

mcols(CBX2_DB_vK27me2)$dpGroup = mcols(CBX2_DB_vK27me2)$Group
mcols(CBX2_DB_vK27me2)$dpGroup %>% unique()
mcols(CBX2_DB_vK27me2)

generateEnrichedHeatmap(CBX2_DB_vK27me2[,,c(1:4,6,7,5,8,9:14)], 
                        matrices_color = c("white",
                                           "purple"), 
                        ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                        sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)", #"CBX2\n(5uM EZH2i)",
                                         "CBX2-CD-GST\n(DMSO)", "CBX2-CD-GST\n(1uM EZH2i)", #"CBX2-CD-GST\n(5uM EZH2i)"
                                         "H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)",
                                          "EZH2\n(DMSO)", "EZH2\n(1uM EZH2i)", #"EZH2\n(5uM EZH2i)"#,
                                         "H2A119ub\n(DMSO)", "H2A119ub\n(1uM EZH2i)",
                                         "RFN2\n(DMSO)", "RFN2\n(1uM EZH2i)",
                                         "HP1beta\n(DMSO)", "HP1beta\n(1uM EZH2i)"
                                         #"H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                        ),
                        #return_ht_list = T,
                        #all_color_scales_equal = F,
                        matrices_axis_name = c("-5kb","center","+5kb"))

CBX2_increased_vsK27me3me2 <- clusterRanges(CBX2_increased_vsK27me3me2, 
              fun = rowMeans, 
              #kmeans_k = 2, 
              cutree_rows = 3,
              silent = T)

CBX2_increased_vsK27_EH <- convertToEnrichedHeatmapMat(CBX2_increased_vsK27me3me2[,,1:2], 
                                                       sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)"#,#"CBX2\n(5uM EZH2i)",
                                                                        #"H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)",
                                                                        #"H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                                                       ))

generateEnrichedHeatmap(CBX2_increased_vsK27me3me2[], 
                        matrices_color = c("white",
                                           "purple"), 
                        #ylim = NULL, #list(c(0,3),c(0,3),c(0,3),c(0,6)),
                        sample_names = c("CBX2\n(DMSO)", "CBX2\n(1uM EZH2i)",#"CBX2\n(5uM EZH2i)",
                                         "H3K27me3\n(DMSO)","H3K27me3\n(1uM EZH2i)",
                                         "H3K27me2\n(DMSO)","H3K27me2\n(1uM EZH2i)"#,"H3K27me2\n(5uM EZH2i)"
                                         ),
                        #return_ht_list = T,
                        #all_color_scales_equal = F,
                        matrices_axis_name = c("-5kb","center","+5kb")) 

EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[1]], 
                                 col = c("white",
                                                      "purple"),
                axis_name = c("-5kb","center","+5kb")) 
  EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[2]], col = c("white",
                                                        "purple"),
                  axis_name = c("-5kb","center","+5kb")) +
  EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[3]], col = c("white",
                                                        "orange"),
                  axis_name = c("-5kb","center","+5kb")) +
  EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[4]], col = c("white",
                                                                         "orange"),
                                   axis_name = c("-5kb","center","+5kb")) +
  EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[5]], col = c("white",
                                                                         "green"),
                                   axis_name = c("-5kb","center","+5kb")) +
  EnrichedHeatmap::EnrichedHeatmap(CBX2_increased_vsK27_EH[[6]], col = c("white",
                                                                         "green"),
                                   axis_name = c("-5kb","center","+5kb")) 

  
  
CBX2_DB_vK27
CBX2_DB_vK27_long = profileplyr::summarize(CBX2_DB_vK27[,,4], fun = rowMeans, output = "long")
CBX2_DB_vK27_long %>% head()
CBX2_DB_vK27_long$Sample %>% unique()
mcols(CBX2_DB_vK27)$Group[mcols(CBX2_DB_vK27)$Group == "Increased"] = "Increased CBX2"

mcols(CBX2_DB_vK27)$dpGroup = mcols(CBX2_DB_vK27)$Group
mcols(CBX2_DB_vK27)
CBX2_DB_vK27
assays(CBX2_DB_vK27[,,4])[[1]] %>% head()


# Plot violin with data in the long format
ggplot(CBX2_DB_vK27_long, 
       aes(x = dpGroup, y = log(Signal))) + 
  geom_violin(draw_quantiles = TRUE, 
              scale = "count", 
              trim = F, 
              fill = "gray") +
  geom_boxplot(width=0.1) + 
  theme_classic() + 
  coord_flip() +
  ylab("H3K27me3 signal (DMSO)") + xlab("") +
  theme(axis.title = element_text(face = "bold"), 
        axis.text.y = element_text(face = "bold"), 
        panel.background = element_rect(fill = "white"))

# Run t-test in data in the long format, assuming unequal variances
t.test(log(Signal) ~ dpGroup, data = CBX2_DB_vK27_long)

# Scatter plot with correlation analysis 
signal <- profileplyr::summarize(CBX2_DB_vK27[,,4], fun = rowMeans, output = "object") %>% assay()
regions <- rowData(profileplyr::summarize(CBX2_DB_vK27[,,4], fun = rowMeans, output = "object"))$names
K27_vsDB = data.frame("name" = regions, "signal" = signal)
K27_vsDB = merge(y= K27_vsDB,
                 x=cbx2.diff.bed[c("name","rep.logFC")], 
                 by = "name",
                 )
K27_vsDB %>% head()
K27_vsDB %>% nrow()

## Plot
ggplot(K27_vsDB, aes(y = rep.logFC, 
                     x = log(RP_030_H3K27me3_500K_DMSO_Rpcells_S30.TMM.binned))) +
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  xlab(expression(bolditalic(Log)~bold("(H3K27me3 signal / DMSO)"))) + 
  ylab(expression(bold(bolditalic(Log[2]) ~ "FC"~"CBX2"~signal))) +
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold"),
        panel.grid = element_blank()) +
  geom_density_2d(#color = "red", 
                  size = .75
                  ) 
  #geom_smooth(method=lm)

cor(y = K27_vsDB$rep.logFC, x = (log(K27_vsDB$RP_030_H3K27me3_500K_DMSO_Rpcells_S30.TMM.binned)), method = "spearman" )




## With signal of both replicates

# Median signal per bin
DB_vK27me3 = import_deepToolsMat("./deepTools/CBX2_DB_vsK27me3_treat_5kb_mean.gz")
mcols(DB_vK27me3)$Group = NA
mcols(DB_vK27me3)$Group[grep("increased", mcols(DB_vK27me3)$names)] = "Increased CBX2"
mcols(DB_vK27me3)$Group[grep("decreased", mcols(DB_vK27me3)$names)] = "Decreased CBX2"
mcols(DB_vK27me3)$dpGroup = mcols(DB_vK27me3)$Group
mcols(DB_vK27me3)$dpGroup %>% unique()
mcols(DB_vK27me3)

# Mean signal per bin
DB_vK27me3 = import_deepToolsMat("./deepTools/CBX2_DB_vsK2me3_5kb_mean_keepZeros.gz")
mcols(DB_vK27me3)$Group = NA
mcols(DB_vK27me3)$Group[grep("increased", mcols(DB_vK27me3)$names)] = "Increased CBX2"
mcols(DB_vK27me3)$Group[grep("decreased", mcols(DB_vK27me3)$names)] = "Decreased CBX2"
mcols(DB_vK27me3)$dpGroup = mcols(DB_vK27me3)$Group
mcols(DB_vK27me3)$dpGroup %>% unique()
mcols(DB_vK27me3)

# Put data in the long format
DB_vK27me3_long = profileplyr::summarize(DB_vK27me3[,,4:5], 
                                         fun = rowMeans, # Average signal from all bins 
                                         output = "long")
DB_vK27me3_long %>% head()
DB_vK27me3_long$Sample %>% unique()
DB_vK27me3_long %>% na.exclude() %>% nrow()


# Run t-test in data in the long format, assuming unequal variances
t.test(log(Signal) ~ dpGroup, data = DB_vK27me3_long)

# Plot violin with data in the long format
ggplot(DB_vK27me3_long, 
       aes(x = dpGroup, y = log(Signal))) + 
  geom_violin(draw_quantiles = TRUE, 
              scale = "count", 
              trim = F, 
              fill = "gray") +
  geom_boxplot(width=0.1) + 
  theme_classic() + 
  coord_flip() +
  ylab(expression(bolditalic(Log)~bold("(H3K27me3 signal / 5uM EZH2i)"))) +
  xlab("") +
  # annotate(geom = "text", 
  #          x = 2.5,
  #          y = 3, 
  #          label = "p-value < 2.2e-16",
  #          fontface = "italic") +
  theme(axis.text.y = element_text(face = "bold", colour = "black"), 
        panel.background = element_rect(fill = "white"))


# Scatter plot with correlation analysis 

## Get the signal for each region in the deepTools matrix
signal <- profileplyr::summarize(DB_vK27me3[,,4:5], fun = rowMeans, output = "object") %>% assay()

## Get the regions
regions <- rowData(profileplyr::summarize(DB_vK27me3[,,4:5], fun = rowMeans, output = "object"))$names

## Create data table with signal for each region
K27_vsDB = data.frame("name" = regions, "signal" = signal)
K27_vsDB = merge(y= K27_vsDB,
                 x= cbx2.diff.bed[c("name","rep.logFC")], 
                 by = "name")

K27_vsDB$AvgSignal = rowMeans(K27_vsDB[3:4])

K27_vsDB %>% head()
K27_vsDB %>% nrow()

## Plot
ggplot(K27_vsDB, aes(y = rep.logFC, 
                     x = log(AvgSignal))) +
  geom_point(alpha = 0.5) + 
  theme_bw() + 
  xlab(expression(bolditalic(Log)~bold("(average H3K27me3 signal / 5uM EZH2i)"))) + 
  ylab(expression(bold(bolditalic(Log[2]) ~ "FC"~"CBX2 (1uM EZH2i / DMSO)"))) +
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold"),
        panel.grid = element_blank()) +
  #geom_density_2d(size = .75) +
  annotate(geom = "text", 
           x = 0.5,
           y = -8, 
           label = "Spearman rank correlation: 0.723", size = 3.5,
           fontface = "italic") + 
  geom_smooth(method=lm, color = "red")

cor(y = K27_vsDB$rep.logFC, x = (log(K27_vsDB$AvgSignal)), method = "spearman" )
# Mean signal per bin
