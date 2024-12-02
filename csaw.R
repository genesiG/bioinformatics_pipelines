BiocManager::install("csaw")

setwd("C:/Users/thephillipslab/Documents/Projects/CBX2 CUT&RUN/csaw")
library(dplyr)
library(csaw)
# Loading data and preparing contrast

dir.bam = "./bam"

CBX2.DMSO.1=file.path(dir.bam,"RP_031_AntiCBX2RP_500K_DMSO_Rpcells_S31.final.sorted.bam")
CBX2.DMSO.2=file.path(dir.bam,"RP_039_AntiCBX2RP_500K_DMSO_Rpcells_S39.final.sorted.bam")

K27me3.DMSO.1="./bam/RP_030_H3K27me3_500K_DMSO_Rpcells_S30.final.sorted.bam"
K27me3.DMSO.2="./bam/RP_038_H3K27me3_500K_DMSO_Rpcells_S38.final.sorted.bam"

IgG.1=file.path(dir.bam,"RP_025_RbIgG_500K_DMSO_Rpcells_S25.final.sorted.bam")
IgG.2=file.path(dir.bam,"RP_033_RbIgG_500K_DMSO_Rpcells_S33.final.sorted.bam")


CBX2.1uM.1="./bam/RP_055_AntiCBX2RP_500K_1uMcpd_Rpcells_S55.final.sorted.bam"
CBX2.1uM.2="./bam/RP_063_AntiCBX2RP_500K_1uMcpd_Rpcells_S63.final.sorted.bam"

IgG.1uM.1 = "./bam/RP_049_RbIgG_500K_1uMcpd_Rpcells_S49.final.sorted.bam"
IgG.1uM.2 = "./bam/RP_057_RbIgG_500K_1uMcpd_Rpcells_S57.final.sorted.bam"

GST.1="./bam/RP_045_GSTControl_500K_DMSO_Rpcells_S45.final.sorted.bam"
GST.2="./bam/RP_069_GSTControl_500K_1uMcpd_Rpcells_S69.final.sorted.bam"
GST.3="./bam/RP_093_GSTControl_500K_5uMcpd_Rpcells_S93.final.sorted.bam"



bam.files <- c(IgG.1, IgG.2, CBX2.DMSO.1, CBX2.DMSO.2)
bam.files <- c(IgG.1uM.1, IgG.1uM.2, CBX2.1uM.1, CBX2.1uM.2)
bam.files <- c(IgG.1, IgG.2, IgG.1uM.1, IgG.1uM.2, 
               CBX2.DMSO.1, CBX2.DMSO.2, CBX2.1uM.1, CBX2.1uM.2)
bam.files <- c(IgG.1, IgG.2, K27me3.DMSO.1, K27me3.DMSO.2)

# Get fragment sizes
library(csaw)

frag.len <- vector()
max.frag <- vector()
for (i in 1:length(bam.files)){
  PEsize <- getPESizes(bam.files[i])
  frag.len[i] = mean(PEsize$sizes)
  max.frag[i] = max(PEsize$sizes)
}

frag.len
max.frag

# Define readParam object
param <- readParam(minq=10, # minimum alingment quality score
                   pe = "both",  # paired-end data
                   max.frag = max(max.frag) + 1, # from getPEsizes analysis
                   dedup = F) # already removed

# Count reads into windows
data.small <- windowCounts(bam.files, 
                     ext = list(frag.len, NA), # ignored in PE data
                     width= 10, # size of window (bp)
                     spacing = 50, # space between each window (bp)
                     #filter = 10, # exclude windows with sum of reads in all libraries < filter
                     param = param)   

#data.small

data.mid <- windowCounts(bam.files, 
                           ext = list(frag.len, NA), # ignored in PE data
                           width= 160, # size of window (bp)
                           spacing = 50, # space between each window (bp)
                           #filter = 10, # exclude windows with sum of reads in all libraries < filter
                           param = param)  

data.large <- windowCounts(bam.files, 
                     ext = list(frag.len, NA), # ignored in PE data
                     width= 500, # size of window (bp)
                     spacing = 400, # space between each window (bp)
                     #filter = 12, # exclude windows with sum of reads in all libraries < filter
                     param = param)   

# data.small.F  = data.small
# data.mid.F = data.mid
# data.large.F = data.large
# data.small  = data.small[,5:8]
# data.mid = data.mid[,5:8]
# data.large = data.large[,5:8]

#----------
# FILTERING OUT LOW-ABUNDANCE WINDOWS
#----------

## Filtering with negative Controls
bins = windowCounts(bam.files,
                    bin=TRUE,
                    width=10000,
                    param=param)

chip.bins = bins[,5:8]

control.bins = bins[,1:4]

scale.info <- scaleControlFilter(chip.bins, control.bins)

filter.small <- filterWindowsControl(data.small[,5:8], 
                                    data.small[,1:4], 
                                    prior.count=5, 
                                    scale.info=scale.info)

summary(filter.small$filter)

filter.mid <- filterWindowsControl(data.mid[,5:8], 
                                    data.mid[,1:4], 
                                    prior.count=5, 
                                    scale.info=scale.info)

summary(filter.mid$filter)

filter.large <- filterWindowsControl(data.large[,5:8], 
                                    data.large[,1:4], 
                                    prior.count=5, 
                                    scale.info=scale.info)
summary(filter.large$filter)


#----------##----------##----------##
## Filtering by GLOBAL background enrichment

bins = windowCounts(bam.files,
                    bin=TRUE,
                    width=10000,
                    param=param)

chip.bins = bins[,5:8]

control.bins = bins[,1:4]

 filter.small <- filterWindowsGlobal(data.small[,5:8], chip.bins, prior.count = 2)
 summary(filter.small$filter)

 filter.mid <- filterWindowsGlobal(data.mid[,5:8], chip.bins, prior.count = 2)
 summary(filter.mid$filter)

 filter.large <- filterWindowsGlobal(data.large[,5:8], chip.bins, prior.count = 2)
 summary(filter.large$filter)

#----------##----------##----------##
## Filtering by LOCAL background enrichment
surrounds <- 1000 # in bp

neighbor <- suppressWarnings(resize(rowRanges(data.small), surrounds, fix="center"))
wider <- regionCounts(bam.files,
                      regions=neighbor,
                      #ext = list(frag.len, NA),
                      param=param)

filter.small <- filterWindowsLocal(data.small, wider, prior.count = 5)
summary(filter.small$filter)

neighbor <- suppressWarnings(resize(rowRanges(data.mid), surrounds, fix="center"))
wider <- regionCounts(bam.files,
                      regions=neighbor,
                      #ext = list(frag.len, NA),
                      param=param)

filter.mid <- filterWindowsLocal(data.mid, wider, prior.count = 5)
summary(filter.mid$filter)

neighbor <- suppressWarnings(resize(rowRanges(data.large), surrounds, fix="center"))
wider <- regionCounts(bam.files,
                      regions=neighbor,
                      #ext = list(frag.len, NA),
                      param=param)

filter.large <- filterWindowsLocal(data.large, wider, prior.count = 5)
summary(filter.large$filter)

#----------
# Check in the histogram if the chosen threshold is greater than the abundances 
# of most bins in the genome, presumably those corresponding to background 
# regions. This suggests that the filter will remove most windows lying within 
# background regions.
#----------
min.fc <- 1.5 # minimum fold-change threshold to filter

hist(filter.mid$filter, main="", breaks=50,
     xlab=expression("Background"~"abundance" ~ (log[2]~"CPM")))
abline(v=log2(min.fc), col="red", lwd = 2)

#----------
# We filter out the majority of windows in background regions upon applying a 
# modest fold-change threshold. This leaves a small set of relevant windows for 
# further analysis.
#----------
chip.small = data.small[,5:8]
chip.mid = data.mid[,5:8]
chip.large = data.large[,5:8]

keep.small <- filter.small$filter > log2(min.fc)
# To investigate the effectiveness of our filtering strategy:
summary(keep.small)

data.small.filt = data.small[keep.small,]
#data.small.filt = chip.small[keep.small,]

keep.mid <- filter.mid$filter > log2(min.fc)
# To investigate the effectiveness of our filtering strategy:
summary(keep.mid)

data.mid.filt = data.mid[keep.mid,]
#data.mid.filt = chip.mid[keep.mid,]

keep.large <- filter.large$filter > log2(min.fc)
# To investigate the effectiveness of our filtering strategy:
summary(keep.large)

data.large.filt = data.large[keep.large,]
#data.large.filt = chip.large[keep.large,]

# NORMALIZATION FOR COMPOSITION BIAS

#----------
# We normalize for composition biases resulting from imbalanced DB between 
# conditions. This is because we expect systematic DB in one direction as one
# is the ChIP sample and the other is the IgG control 
#----------

 bins = windowCounts(bam.files,
                     bin=TRUE,
                     width=10000,
                     param=param)

# Using the calculated bins 
data.small.filt <- normFactors(bins, se.out=data.small.filt)
data.mid.filt <- normFactors(bins, se.out=data.mid.filt)
data.large.filt <- normFactors(bins, se.out=data.large.filt)
(normfacs <- data.mid.filt$norm.factors)

# NORMALIZATION FOR IP EFFICIENCY BIAS 

# data.small.filt <- normFactors(data.small.filt)
# data.mid.filt <- normFactors(data.mid.filt)
# data.large.filt <- normFactors(data.large.filt)
# (normfacs <- data.large.filt$norm.factors)

 saveRDS(data.small.filt, "./data.K27me3.DMSO.filt_window10_spacing50_filter1.5_local1k.rds")
 saveRDS(data.mid.filt, "./data.K27me3.DMSO.filt_window160_spacing50_filter1.5_local1k.rds")
 saveRDS(data.large.filt, "./data.K27me3.DMSO.filt_window500_spacing400_filter1.5_local1k.rds")



#----------
# Let's see the effect of normalization on the relative enrichment between 
# pairs of samples. 
#----------

bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)

par(cex.lab=1.5, mfrow=c(1,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,2], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (1 vs 2)")
abline(h=log2(normfacs[1]/normfacs[2]), col="red")

smoothScatter(bin.ab, adjc[,1]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (1 vs 3)")
abline(h=log2(normfacs[1]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,1]-adjc[,4], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (1 vs 4)")
abline(h=log2(normfacs[1]/normfacs[4]), col="red")

#----------
# We see that log-ratio of normalization factors passes through the centre of 
# the cloud of background regions in each plot, indicating that the bias has 
# been successfully identified and removed.
# 
# Clouds at low A-values represent background and clouds at high A-values 
# represent bound regions
#----------

# STATISTICAL MODELING 

#----------
# We first convert our RangedSummarizedExperiment object into a DGEList for 
# modelling with edgeR.
#----------

library(edgeR)
y.small <- asDGEList(data.small.filt)
y.mid <- asDGEList(data.mid.filt)
y.large <- asDGEList(data.large.filt)

# We then construct a design matrix for our experimental design. Here, we use a 
# simple one-way layout with two groups of two replicates.

grouping <- factor(c('IgG','IgG',
                     #'IgG','IgG',
                     #'CBX2_DMSO','CBX2_DMSO',
                     #'CBX2_1uM','CBX2_1uM'
                     'K27me3_DMSO','K27me3_DMSO'
                     ))
design <- model.matrix(~0 + grouping)
colnames(design) <- levels(grouping)
design

#----------
# We estimate the negative binomial (NB) dispersions for each window.
#----------

y.small <- estimateDisp(y.small, design)
y.mid <- estimateDisp(y.mid, design)
y.large <- estimateDisp(y.large, design)

#----------
# We again observe an increasing trend in the NB dispersions with respect to 
# abundance (plotBCV).
#----------

summary(y.small$trended.dispersion)
summary(y.mid$trended.dispersion)
summary(y.large$trended.dispersion)

 plotBCV(y.small)
 plotBCV(y.mid)
 plotBCV(y.large)
# Abundance-dependent trend in the BCV for each window, represented 
# by the blue line. Common (red) and tagwise estimates (black) are 
# also shown.

plotMDS(cpm(y.small, log=TRUE), top=10000, labels=grouping,
        col=c("purple", "black")[as.integer(grouping)])
plotMDS(cpm(y.mid, log=TRUE), top=10000, labels=grouping,
        col=c("purple", "black")[as.integer(grouping)])
plotMDS(cpm(y.large, log=TRUE), top=10000, labels=grouping,
        col=c("purple", "black")[as.integer(grouping)])
#----------
# We estimate the quasi-likelihood (QL) dispersions for each window.
#----------

fit.small <- glmQLFit(y.small,
                      #dispersion=0.05,
                      design, 
                      robust=TRUE
                      )
fit.mid <- glmQLFit(y.mid, 
                    #dispersion=0.05,
                    design, 
                    robust=TRUE
                    )
fit.large <- glmQLFit(y.large, 
                      #dispersion=0.05,
                      design, 
                      robust=TRUE
                      )

#----------
# If the QL dispersions are strongly shrunk towards the trend, that indicates
# that there is little variability in the dispersions across windows.
#----------

summary(fit.small$df.prior)
summary(fit.mid$df.prior)
summary(fit.large$df.prior)

plotQLDisp(fit.small)
plotQLDisp(fit.mid)
plotQLDisp(fit.large)



#----------
# TESTING FOR DIFFERENTIAL BINDING
#----------

contrast <- makeContrasts(K27me3_DMSO - IgG, levels=design)

# We then test for DB between conditions in each window using the QL F-test.

res.small <- glmQLFTest(fit.small, contrast=contrast)
res.mid <- glmQLFTest(fit.mid, contrast=contrast)
res.large <- glmQLFTest(fit.large, contrast=contrast)

saveRDS(file="res.small-K27me3.DMSO-width10spacing50-filter1.5_local10k.rds", res.small)
saveRDS(file="res.mid-K27me3.DMSO-width160spacing50-filter1.5_local10.rds", res.mid)
saveRDS(file="res.large-K27me3.DMSO-width500spacing400-filter1.5_local10.rds", res.large)

  res.small = readRDS(file="res.small-K27me3.DMSO-width10spacing50-filter1.0_local10k.rds")
  res.mid = readRDS(file="res.mid-K27me3.DMSO-width160spacing50-filter1.0_local10.rds")
  res.large = readRDS(file="res.large-K27me3.DMSO-width500spacing400-filter1.0_local10.rds")

 data.small.filt <- readRDS("./data.K27me3.DMSO.filt_window10_spacing50_filter1.0_local1k.rds")
 data.mid.filt = readRDS("./data.K27me3.DMSO.filt_window160_spacing50_filter1.0_local1k.rds")
 data.large.filt = readRDS("./data.K27me3.DMSO.filt_window500_spacing400_filter1.0_local1k.rds")

#----------
# CONSOLIDATING RESULTS
#----------

merged <- mergeResultsList(list(data.small.filt, 
                                data.mid.filt, 
                                data.large.filt
                                ), 
                           tab.list=list(res.small$table, 
                                         res.mid$table, 
                                         res.large$table
                                         ),
                           equiweight=T, 
                           tol=600, 
                           merge.args=list(max.width=50000)
                           )
# merged$regions

tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.001
summary(is.sig)
summary((tabcom$rep.logFC >= 4)[is.sig]) # decreased

table(tabcom$direction[is.sig])

tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC >= 4)[is.sig] # increased
summary(is.sig.pos)

# Extract binding regions
cbx2.csaw.ranges <- merged$regions
mcols(cbx2.csaw.ranges) <- data.frame(tabcom, best.logFC=tabbest$rep.logFC)

#cbx2.csaw.ranges
#saveRDS(file="cbx2-small.w10s50-large.w150.s100-local10k-merge1k.rds", cbx2.csaw.ranges)

# cbx2.csaw.ranges <- readRDS("./cbx2-width100spacing25-filter12-q001-local10k-merge1k-logFC2.rds")

cbx2.csaw.bed =  cbx2.csaw.ranges %>% as.data.frame() %>% filter(#direction == "up",
                                                                 FDR <0.001, 
                                                                 rep.logFC >= 5,
                                                                 #rep.logFC <= -0.6,
                                                                 width >= 160
                                                                 #width >= mean(frag.len)
                                                                 )
cbx2.csaw.bed$score = -10*log10(1e-10+cbx2.csaw.bed$FDR)
cbx2.csaw.bed$name = paste("region",1:nrow(cbx2.csaw.bed))
#cbx2.csaw.bed$name = paste("decreased",1:nrow(cbx2.csaw.bed))

width = cbx2.csaw.bed[c("width", "score")] %>% arrange(desc(score))
summary(width)
plot(density(log10(width$width[])),
     lwd = 2,
     col = "red",
     #main = "Top 10000 binding regions", 
     #xlab = expression(bold("Region"~ "width")~ bolditalic((Log[10]))),
     ylab = expression(bold("Density")))
#abline(v = median(log10((width$width))), lwd = 1, lty = "dashed", col = "darkgreen")
#abline(v = log10(1000), lwd = 1, lty = "dashed", col = "darkgreen")

rtracklayer::export.bed(cbx2.csaw.bed, "./K27me3_DMSO.filterlocal1.0.merge600.q0001.FC5.max50k.bed")

#cbx2_regions <- import.bed("./cbx2-small.w10s50-large.w500.s100-local10k-merge1k.FC2.bed")
#cbx7.csaw.bed = cbx2.csaw.bed
cbx2.diff.bed = cbx2.csaw.bed
nrow(cbx2.diff.bed)
cbx2.diff.bed = rbind(cbx2.diff.bed, cbx2.csaw.bed)
head(cbx2.diff.bed)
tail(cbx2.diff.bed)

rtracklayer::export.bed(cbx2.diff.bed, "./DB/CBX2_diff.merge10.q005.FC0.6.max10k.bed")

cbx2.diff.bed[c("name","rep.logFC")]
lFC = cbx2.diff.bed$rep.logFC
