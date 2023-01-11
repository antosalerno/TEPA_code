###### Dimensionality reduction and clustering of immune cells ##
## author: Antonietta Salerno
## date: 22/12/2022

library("Seurat")
library("ggplot2")
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
# load("TEPA_results/00_immune.rda")
immune <- LoadH5Seurat("TEPA_results/00_immune.h5Seurat")

#### 1 - QC and filtering of immune cells ####

immune <- NormalizeData(immune) 
# Get gene names
genenames<-rownames(immune@assays$RNA@counts)
# Get mitochondrial gene names
mitogenes<-grep("^mt",genenames,value=TRUE)
# Add percent mito data in seurat object
immune[["percent.mt"]] <- PercentageFeatureSet(immune, pattern = "^mt.")

png("TEPA_plots/01_preFilterQC_scatter.png", w = 4000, h = 2000, res = 300)
plot1 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2
dev.off()

# Visualize the distribution of mitochondrial gene expression detected per cell
png("TEPA_plots/01_preFilterQC_mito.png", h = 3000, w = 4200, res = 300)
immune@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
dev.off()

immune <- subset(immune, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 &
                   percent.mt < 25 & nCount_RNA < 3000) # all above 50 genes per cell anyways

png("TEPA_plots/01_postFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/01_postFilterQC_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(immune, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 2)
dev.off()

png("TEPA_plots/01_postFilterQC_scatter.png", h = 3000, w = 4200, res = 300)
plot1 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(immune, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 # We can spot two separate data clouds because of the bimodal distribution of the data: tumour vs immune cells
dev.off()

# save(immune, file = "TEPA_results/01_immuneFilt.rda")
SaveH5Seurat(immune, filename = "TEPA_results/01_immuneFilt.h5Seurat", overwrite = TRUE)

#### 2 - Batch effect correction ####

# load("TEPA_results/01_immuneFilt.rda")
immune <- LoadH5Seurat("01_immuneFilt.h5Seurat")

samples.list <- SplitObject(immune, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = samples.list)
immune.anchors <- FindIntegrationAnchors(object.list = samples.list,anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

# save(immune.combined, file = "TEPA_results/01_immuneInt.rda")
SaveH5Seurat(immune.combined, filename = "TEPA_results/01_immuneInt.h5Seurat", overwrite = TRUE)

#### 3 - Dimensionality reduction ####

#load("TEPA_results/01_immuneInt.rda")
immune.combined <- LoadH5Seurat("01_immuneInt.h5Seurat")

# Scaling
seuset_immune <- ScaleData(immune.combined, verbose = FALSE)
seuset_immune <- RunPCA(seuset_immune, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_immune@reductions$pca@stdev / sum(seuset_immune@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 70.37% of variability

seuset_immune <- FindNeighbors(object = seuset_immune, dims = 1:60, reduction = 'pca')
seuset_immune <- FindClusters(object = seuset_immune, resolution = 0.28) # Louvain algorithm with multi-level refinement
seuset_immune <- RunUMAP(seuset_immune, dims = 1:60, reduction = "pca", verbose = FALSE)

png("TEPA_plots/01_umapExplore.png", w = 6000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("orig.ident", "condition", "seurat_clusters"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/01_umapClust.png", w = 4000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("seurat_clusters"), split.by= "condition",label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/01_umapCounts.png", w = 4000, h = 2000, res = 300)
FeaturePlot(seuset_immune, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

# save(seuset_immune, file = "TEPA_results/01_seusetImmune.rda")
SaveH5Seurat(seuset_immune, filename = "TEPA_results/01_seusetImmune.h5Seurat", overwrite = TRUE)






