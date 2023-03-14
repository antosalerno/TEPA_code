#### Analysis GeoMx DSP data ####
## author: Antonietta Salerno
## date: 07/03/2022
BiocManager::install("NanoStringNCTools")
BiocManager::install("GeomxTools")
BiocManager::install("GeoMxWorkflows")

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(ggcharts)
library(MAST)
library(openxlsx)
library(patchwork)
library(clusterProfiler)
library("tibble")
library("Seurat")
library("SeuratDisk")
library("fgsea")

### 1 - Load data ####
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")

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
                               protocolDataColNames = c("area", "roi", "Condition", 
                                                        "Core", "Infiltration",
                                                        "AOISurfaceArea", "AOINucleiCount",
                                                        "ROICoordinateX", "ROICoordinateY"),
                               experimentDataColNames = c("Core", "aoi"))

# Clean the column with Core (C,C3,C7,T3,T7)
pData(protocolData(data))$Core <- 
  str_split(pData(protocolData(data))$Core, "\\-", simplify=T)[,1]

# Add a column with Group (Control, T3, T7)
group <- pData(protocolData(data))$Core
control <-  group %in% c("C","C3","C7")
group[control] <- "C"
pData(protocolData(data))$Group <- group

### 2 - Study Design ####

library(knitr)
library(dplyr)
library(ggforce)

count_mat <- dplyr::count(pData(protocolData(data)), Core, Group, Condition, Infiltration) %>%
  mutate(Condition = as.character(Condition)) %>% 
  mutate(Infiltration = as.character(Infiltration))

# gather the data and plot in order: 
test_gr <- gather_set_data(count_mat, 1:3)
test_gr$x <- factor(test_gr$x)
levels(test_gr$x) = c("Core","Group","Condition", "Infiltration")

# plot Sankey
save = "N00_sankeyCore"
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
           y = 0, yend = 105, lwd = 2) +
  annotate(geom = "text", x = 3.19, y = 63, angle = 90, size = 5,
           hjust = 0.5, label = "100 ROIs")
ggsave(p, file=paste0("TEPA_plots/", save, ".png"), width = 30, height = 30, units = "cm")

#### 3 - QC and preprocessing ####

### 3.1 - Shift counts to one###
data <- shiftCountsOne(data, useDALogic = TRUE)

### 3.2 - Flag low quality ROIs ###
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000) 
       percentTrimmed = 85,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%) ->75
       percentSaturation = 65, # Minimum sequencing saturation (50%)
       minNegativeCount = 1   # Minimum negative control counts (10) -> 1
  )   

data <- setSegmentQCFlags(data, qcCutoffs = QC_params)

### 3.3 - Flag low quality probes ###
data <- setBioProbeQCFlags(data, qcCutoffs = 
                             list(minProbeRatio = 0.1,percentFailGrubbs = 20), 
                           removeLocalOutliers = TRUE)

ProbeQCResults <- fData(data)[["QCFlags"]]
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df

### 3.4 - Remove low quality ROIs and probes ###

passedQC <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(passedQC)

data <- passedQC 

### 3.5 - Create gene-level count data ###
# Objects must be aggregated to Target level data before coercing. 
# This changes the row (gene) information to be the gene name rather than the probe ID.

target_data <- aggregateCounts(passedQC)

### 3.6 - Normalisation via 3rd quantile ###
norm_target_data <- normalize(target_data, norm_method="quant",
                                  desiredQuantile = .75, toElt = "q_norm")

#assayDataElementNames(norm_target_data)

# skipping removal of genes with low detection rates here

#### 4 - Coercion to Seurat object ####
seuset_nano <- as.Seurat(norm_target_data, normData = "q_norm", ident = "Group", 
                         coordinates = c("ROICoordinateX", "ROICoordinateY"))
Idents(seuset_nano) <- factor(x = Idents(seuset_nano), levels = c("C", "T3", "T7"))

SaveH5Seurat(seuset_nano, filename = "TEPA_results/N00_seusetNano.h5Seurat", overwrite = TRUE)

#head(seuset_nano@misc$QCMetrics$QCFlags) 

png("TEPA_plots/N00_countsROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "nCount_GeoMx", split.by = "Infiltration",
        pt.size = 0.1)
dev.off() # the number of genes is instead the same for all ROIs

png("TEPA_plots/N00_nucleiROIs.png", h = 3000, w = 2500, res = 300)
VlnPlot(seuset_nano, features = "AOINucleiCount", split.by = "Group",
        pt.size = 0.1)
dev.off() 

# ### Correct for batch effect ###
# 
# samples.list <- SplitObject(seuset_nano, split.by = "Group")
# 
# # Normalize and identify variable features for each dataset independently (Treatment vs Control)
# samples.list <- lapply(X = samples.list, FUN = function(x) {
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })
# 
# features <- SelectIntegrationFeatures(object.list = samples.list)
# nano.anchors <- FindIntegrationAnchors(object.list = samples.list,
#                                        dims = 1:4, 
#                                        k.anchor = 1, k.filter = 1,
#                                        anchor.features = features)
# nano.combined <- IntegrateData(anchorset = nano.anchors)


# Then go to a new file and first do the tests (N01) and 

#### 5 - Dimensionality reduction ####

seuset_nano <- FindVariableFeatures(seuset_nano)
seuset_nano <- ScaleData(seuset_nano)
seuset_nano <- RunPCA(seuset_nano, assay = "GeoMx", verbose = FALSE, approx=FALSE)

# Determine percent of variation associated with each PC
pct <- seuset_nano@reductions$pca@stdev / sum(seuset_nano@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=50) # Select 60 PCs to retain 73.8% of variability

seuset_nano <- FindNeighbors(seuset_nano, reduction = "pca", dims = 1:50)
seuset_nano <- FindClusters(seuset_nano, resolution = 0.8)
seuset_nano <- RunUMAP(seuset_nano, reduction = "pca", dims = 1:50)

png("TEPA_plots/N00_clustROIs.png", h = 3000, w = 2500, res = 300)
DimPlot(seuset_nano, reduction = "umap", pt.size = 5,
        label = F, group.by = "seurat_clusters")
dev.off()
# 3 groups, would they reflect Core?

png("TEPA_plots/N00_umapExplore.png", w = 6000, h = 4000, res = 300)
DimPlot(object = seuset_nano, pt.size = 5, reduction = 'umap', ncol = 2, 
        group.by = c("Group", "Infiltration","Condition","seurat_clusters"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_nano@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

Idents(seuset_nano) <- "seurat_clusters"
clusters = levels(Idents(seuset_nano))

# 5.1 - Find markers for every cluster 
nano.markers <- FindAllMarkers(seuset_nano, 
                                 only.pos = FALSE, 
                                 #min.pct = 0.5, 
                                 #min.diff.pct = 0.2,
                                 #logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="Core")
nano.markers$p_val_adj = p.adjust(nano.markers$p_val, method='BH')
write.csv(nano.markers, "TEPA_results/N00_DEA_clusterMarkers.csv")

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = nano.markers[nano.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/N00_DEA_clusterMarkers.xlsx", overwrite = TRUE)

png("TEPA_plots/N00_nanoMarkersHeatmap.png", h = 4000, w = 6000, res = 300)
top10 <- nano.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(object = subset(seuset_nano, downsample = 500), size = 6, 
           features = top10$gene) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

#### Differential expression analysis ####
# Check in parallel expression of DEA genes in the different cell types
seuset_immune <- LoadH5Seurat("TEPA_results/02_immuneAnn.h5Seurat")

# 5.2 Infiltration vs Non-Infiltration given Treatment

Idents(seuset_nano) <- "Infiltration"
seuset_nanoTEPA <- subset(seuset_nano, Condition == "Treatment")
res <- FindMarkers(seuset_nanoTEPA, ident.1 = "T", ident.2 = "F",
                   only.pos = FALSE, verbose = FALSE,
                   latent.vars="Core",
                   test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoInf_gCond_DEA.csv"))

log2FC = 0.25
save = "N00_nanoInf_gCond_"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA Infiltrated vs Non-Infiltrated TME (TEPA-treated ROIs)",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 3,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

ranked.genes<- res %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = ranked.genes,  split.by = "condition",
        scale=TRUE, dot.scale = 5,
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
Idents(seuset_nano) <- "Infiltration"
DoHeatmap(object = subset(seuset_nanoTEPA, downsample = 500), size = 6, 
          features = head(ranked.genes, 20)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# 5.4 Treatment vs Control given Infiltration

Idents(seuset_nano) <- "Condition"
seuset_nanoINF <- subset(seuset_nano, Infiltration == "T")
res <- FindMarkers(seuset_nanoINF, ident.1 = "Treatment", ident.2 = "Control",
                   only.pos = FALSE, verbose = FALSE,
                   latent.vars="Core",
                   test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoCond_gInf_DEA.csv"))

log2FC = 0.3
save = "N00_nanoCond_gInf_"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA Treatment vs Control given Infiltration",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 3,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

ranked.genes<- res %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = ranked.genes,  split.by = "condition",
        scale=TRUE, dot.scale = 5,
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
Idents(seuset_nanoINF) <- "Condition"
DoHeatmap(object = subset(seuset_nanoINF, downsample = 500), size = 6, 
          features = head(ranked.genes, 20)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# 5.5 Treatment vs Control

Idents(seuset_nano) <- "Condition"
res <- FindMarkers(seuset_nano, ident.1 = "Treatment", ident.2 = "Control",
                   only.pos = FALSE, verbose = FALSE,
                   latent.vars="Core",
                   test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/N00_nanoCond_DEA.csv"))

log2FC = 0.3
save = "N00_nanoCond_"

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA Treatment vs Control",
                     pCutoff = 0.05, 
                     FCcutoff = log2FC,
                     labFace = "bold",
                     labSize = 3,
                     col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                     colAlpha = 4/5,
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                     caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>log2FC&res$p_val_adj<=0.05,]), ' genes',
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

ranked.genes<- res %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  filter(avg_log2FC > abs(0.3))
ranked.genes <- rownames(ranked.genes)[!is.na(rownames(ranked.genes))]

png(paste0("TEPA_plots/",save,"DotPlot.png"), h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = ranked.genes,  split.by = "condition",
        scale=TRUE, dot.scale = 5,
        assay = "RNA", cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

png(paste0("TEPA_plots/",save,"Heatmap.png"), h = 4000, w = 6000, res = 300)
DoHeatmap(object = subset(seuset_nano, downsample = 500), size = 6, 
          features = head(ranked.genes, 20)) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()


SaveH5Seurat(seuset_nano, filename = "TEPA_results/N00_seusetNanoRed.h5Seurat", overwrite = TRUE)

#### 6 - Pathway enrichment Analysis ####
seuset_nano <- LoadH5Seurat("TEPA_results/N00_seusetNanoRed.h5Seurat")

### 6.1 - Cluster markers
### 6.2 - Infiltration vs Non-Infiltration given Treatment

sets1 <- read.gmt("TEPA_data/mh.all.v2022.1.Mm.symbols.gmt") # Mouse hallmark
#sets2 <- read.gmt("TEPA_data/m2.cp.reactome.v2022.1.Mm.symbols.gmt") # Reactome
#sets3 <- read.gmt("TEPA_data/m8.all.v2023.1.Mm.symbols.gmt") # Cell types

sets1$term <- as.character(sets1$term)
#sets2$term <- as.character(sets2$term)
#sets3$term <- as.character(sets3$term)

sets1 <- sets1 %>% split(x = .$gene, f = .$term)
#sets2 <- sets2 %>% split(x = .$gene, f = .$term)
#sets3 <- sets3 %>% split(x = .$gene, f = .$term)

#fgsea_sets <- append(sets1, sets2)
#fgsea_sets <- append(fgsea_sets, sets3)
fgsea_sets <- append(sets1, custom)

save = "N00_infTEPA_Enrichment_"
Idents(seuset_nano) <- "Infiltration"
seuset_nanoINF <- subset(seuset_nano, Infiltration == "T")
clusters = levels(Idents(seuset_nanoINF))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "nanoInf_gCond")

### 6.3 - Treatment vs Control given Infiltration
save = "N00_TEPAInf_Enrichment_"
Idents(seuset_nano) <- "Condition"
seuset_nanoTEPA <- subset(seuset_nano, Condition == "Treatment")
clusters = levels(Idents(seuset_nanoTEPA))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "nanoCond_gInf")

### 6.4 - Treatment vs Control given Infiltration
save = "N00_TEPA_Enrichment_"
Idents(seuset_nano) <- "Condition"
clusters = levels(Idents(seuset_nano))
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "nanoCond")

#### 7 - Spatial graphing ####

id <- "Core"
Idents(seuset_nano) <- id
plots <- list()
for(i in levels(Idents(seuset_nano))){
  p <- suppressMessages(SpatialFeaturePlot(seuset_nano[,seuset_nano[[id]] == i], 
                                      features = "nCount_GeoMx", pt.size.factor = 12) + 
                     labs(title = i) + 
                     theme(legend.position = "none") + 
                     scale_fill_continuous(type = "viridis",
                                           limits = c(min(seuset_nano$nCount_GeoMx), 
                                                      max(seuset_nano$nCount_GeoMx))))
  plots <- c(plots, list(p))
}
png("TEPA_plots/N00_spatialCore.png", h = 4000, w = 6000, res = 300)
wrap_plots(plots)
dev.off()

# Plot a specific feature
plots <- list()
gene = "Itgb7"
for(i in levels(Idents(seuset_nano))){
  p <- suppressMessages(SpatialFeaturePlot(seuset_nano[,seuset_nano[[id]] == i], 
                                           features = gene, pt.size.factor = 12) + 
                          labs(title = i) + 
                          theme(legend.position = "none") + 
                          scale_fill_continuous(type = "viridis",
                                                limits = c(min(seuset_nano@assays$GeoMx@counts[gene,]), 
                                                           max(seuset_nano@assays$GeoMx@counts[gene,]))))
  plots <- c(plots, list(p))
}
png("TEPA_plots/N00_spatialCore_Itgb7.png", h = 4000, w = 6000, res = 300)
wrap_plots(plots)
dev.off()

# no spatial meaningggg










