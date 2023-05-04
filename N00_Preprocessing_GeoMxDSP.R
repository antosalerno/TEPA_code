#### Pre-processing GeoMx DSP data ####
## author: Antonietta Salerno
## date: 06/01/2022
BiocManager::install("NanoStringNCTools")
BiocManager::install("GeomxTools")
BiocManager::install("GeoMxWorkflows")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)

# The ROIs 37,38,93, 94,95 have been removed for the quality of the samples: necrotic tumor
# We discard also the quality control samples Kidney and Liver from the analysis 

### 1 - Load data ####
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
datadir <- "TEPA_data"
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir), pattern = ".pkc$",
                                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <- file.path(datadir, "annotationFull.xlsx")

data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template",
                               phenoDataDccColName = "Sample_ID",
                               protocolDataColNames = c("area", "roi", "Tissue", "Group", 
                                                        "Core", "Infiltration", 
                                                        "AOISurfaceArea", "AOINucleiCount",
                                                        "ROICoordinateX", "ROICoordinateY"),
                               experimentDataColNames = c("Core", "aoi"))
  
### 2 - Study Design ####

library(knitr)
library(dplyr)
library(ggforce)

pkcs <- annotation(data)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

rep_str = c('T'='Infiltrated TME','F'='Non-Infiltrated TME')
pData(protocolData(data))$Infiltration <- str_replace_all(pData(protocolData(data))$Infiltration, rep_str)

count_mat <- dplyr::count(pData(protocolData(data)), Core, Condition, Infiltration) %>%
  mutate(Tissue = as.character(Core)) %>%
  mutate(Group = as.character(Condition)) %>% 
  mutate(Infiltration = as.character(Infiltration))

# gather the data and plot in order: 
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x)
levels(test_gr$x) = c("Core","Condition","Infiltration")

count_mat %>%
  group_by(Condition, Infiltration) %>%
  mutate(cond_inf = sum(n)) 

# tot needs to be 100 ROIs -> why C is missing?
# Percentage of infiltration: C = 29%, T = 83%

      
# plot Sankey
save = "N00_sankey"
p <- ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Infiltration), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 3.25, xend = 3.25,
           y = 0, yend = 105, lwd = 2) +
  annotate(geom = "text", x = 3.19, y = 63, angle = 90, size = 5,
           hjust = 0.5, label = "100 ROIs")
ggsave(p, file=paste0("TEPA_final_figures/", save, ".pdf"), width = 30, height = 30, units = "cm")

  
#### 3 - QC and preprocessing ####

# Shift counts to one
data <- shiftCountsOne(data, useDALogic = TRUE)

### 3.1 - Segment QC ###
# We first assess sequencing quality and adequate tissue sampling for every ROI/AOI segment.

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000) 
       percentTrimmed = 85,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%) ->75
       percentSaturation = 65, # Minimum sequencing saturation (50%)
       minNegativeCount = 1   # Minimum negative control counts (10) -> 1
       )         

data <- setSegmentQCFlags(data, qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

col_by <- "Group"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "ROIs, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans, limits=c(NA,1))
  }
  plt
}

plt <- QC_histogram(sData(data), "Trimmed (%)", col_by, 85)
ggsave(plt, file = "TEPA_plots/N00_histTrimmed.png",width = 30, height = 25, units = "cm" )

plt <- QC_histogram(sData(data), "Stitched (%)", col_by, 80)
ggsave(plt, file = "TEPA_plots/N00_histStitched.png",width = 30, height = 25, units = "cm" )

plt <- QC_histogram(sData(data), "Aligned (%)", col_by, 75)
ggsave(plt, file = "TEPA_plots/N00_histAligned.png",width = 30, height = 25, units = "cm" )

plt <- QC_histogram(sData(data), "Saturated (%)", col_by, 65) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
ggsave(plt, file = "TEPA_plots/N00_histSeqSat.png",width = 30, height = 25, units = "cm" )

# Calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 

negCols <- paste0("NegGeoMean_", modules)
pData(protocolData(data))[, negCols] <- negativeGeoMeans
for(ann in negCols) {
  plt <- QC_histogram(pData(protocolData(data)), ann, col_by, 2)
  ggsave(plt, file = paste0("TEPA_plots/N00_negCols.png"), width = 30, height = 25, units = "cm" )
}

# Detach neg_geomean columns ahead of aggregateCounts call
pData(data) <- pData(data)[, !colnames(pData(data)) %in% negCols]

data <- data[, QCResults$QCStatus == "PASS"]

# Subsetting our dataset has removed samples which did not pass QC
dim(data) # 108 >>> 101

### 3.2 - Probe QC ###

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
data <- setBioProbeQCFlags(data, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(data)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

data <- ProbeQCPassed 

### 3.3 - Create gene-level count data ###

# Check how many unique targets the object has
length(unique(featureData(data)[["TargetName"]])) # 19963

# collapse to targets
target_data <- aggregateCounts(data)
dim(target_data) # Features: 19963  ROIs: 101

### 3.4 - Limit of Quantification (Detection Rates) ###

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested 
LOQ <- data.frame(row.names = colnames(target_data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if(all(vars[1:2] %in% colnames(pData(target_data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_data)[, vars[1]] * 
             pData(target_data)[, vars[2]] ^ cutoff)
  }
}
pData(target_data)$LOQ <- LOQ

### 3.5 - Filtering out either segments and/or genes with abnormally low signal ###
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_data)$Module == module
  Mat_i <- t(esApply(target_data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_data)$TargetName, ]

### 3.6 - By segment gene detection rates ###

# Save detection rate information to pheno data
pData(target_data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_data)$GeneDetectionRate <-
  pData(target_data)$GenesDetected / nrow(target_data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_data)$DetectionThreshold <- 
  cut(pData(target_data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
plt <- ggplot(pData(target_data),
       aes(x = DetectionThreshold)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "ROIs, #")
ggsave(plt, file = paste0("TEPA_plots/N00_geneDetThr.png"), width = 30, height = 25, units = "cm" )

# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_data)$DetectionThreshold,
            pData(target_data)$class))

target_data <-
  target_data[, pData(target_data)$GeneDetectionRate >= .1]

dim(target_data)

# Plot Sankey with filtered dataset #

count_mat <- dplyr::count(pData(protocolData(data)), Tissue,  Group, Infiltration) %>%
  mutate(Tissue = as.character(Tissue)) %>%
  mutate(Group = as.character(Group)) %>% 
  mutate(Infiltration = as.character(Infiltration))

# gather the data and plot in order: 
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x)
levels(test_gr$x) = c("Tissue","Group", "Infiltration")

# plot Sankey
save = "N00_sankeyFilt"
p <- ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = Infiltration), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "#E3B9B1", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 3.25, xend = 3.25,
           y = 0, yend = 106, lwd = 2) +
  annotate(geom = "text", x = 3.19, y = 55, angle = 90, size = 5,
           hjust = 0.5, label = "101 ROIs")
ggsave(p, file=paste0("TEPA_plots/", save, ".png"), width = 30, height = 30, units = "cm")

### Gene detection rates ###
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_data)]
fData(target_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_data)$DetectionRate <-
  fData(target_data)$DetectedSegments / nrow(pData(target_data))

# Genes of interest detection table 

N1 <- c("S100a8", "S100a9", "Icam1", "Fas", "Tnf", "Isg15", "Isg20", 
        "Ccl3", "Ccl4",  "Cxcl2", "Cxcl3", "Cebpb", "Il1b", "Il1r2", "Il1rn", "Il6ra", "Il15",
        "Stat3", "Hif1a", "Ifitm1","Ifitm2",  "Ifitm3", "Ifitm6",  
        "Acod1", "Myd88", "Prkcd", "Mmp8", "Mmp9","Retnlg", "Arg2")

goi_df <- data.frame(
  Gene = N1,
  Number = fData(target_data)[N1, "DetectedSegments"],
  DetectionRate = percent(fData(target_data)[N1, "DetectionRate"]))

# Filter out low detection genes

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_data)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_data))
rownames(plot_detect) <- plot_detect$Freq

plt <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of ROIs",
       y = "Genes Detected, % of Panel > LOQ")

ggsave(plt, file=paste0("TEPA_plots/N00_N1_DetectionRates.png"), width = 30, height = 30, units = "cm")

# Focus on the genes detected in at least 10% of our segments: maybe too stringent

#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_data <- 
  target_data[fData(target_data)$DetectionRate >= 0.1 |
                    fData(target_data)$TargetName %in% neg_probes, ]
dim(target_data) # Genes: 14301  ROIs: 101 

# retain only detected genes of interest
goi <- N1[N1 %in% rownames(target_data)]


#### 4 - Normalisation ####

library(reshape2)  # for melt
library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "Group"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_data)),
             Segment = colnames(exprs(target_data)),
             Annotation = pData(protocolData(target_data))[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_data), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_data)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "ROIs, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plt <- plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
ggsave(plt, file=paste0("TEPA_plots/N00_Norm.png"), width = 50, height = 30, units = "cm")

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_data <- normalize(target_data ,
                             norm_method = "quant", 
                             desiredQuantile = .75,
                             toElt = "q_norm")

png("TEPA_plots/N00_NormFullQ3.png", h = 2000, w = 2500, res = 300)
boxplot(assayDataElement(target_data[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "ROI",
        ylab = "Counts, Q3 Normalized")
dev.off()

png("TEPA_plots/N00_NormFullNeg.png", h = 2000, w = 2500, res = 300)
boxplot(assayDataElement(target_data[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "ROI",
        ylab = "Counts, Neg. Normalized")
 dev.off()
 
 
