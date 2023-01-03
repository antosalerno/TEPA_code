###### Subset tumour cells from immune cells ##
## author: Antonietta Salerno
## date: 21/12/2022

library("Seurat")
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")

#### 1 - Create Seurat object with all the samples ####

#d7c : day 7 - control
counts <- read.table("TEPA_data/Combined_D7_CONTROL_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
CF <- CreateSeuratObject(counts = t(counts), project="CF", min.cells = 3, min.features = 500)
#d7cm : day 7 - control myeloid 
counts <- read.table("TEPA_data/Combined_D7_CONTROLMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
CM <- CreateSeuratObject(counts = t(counts), project="CM", min.cells = 3, min.features = 500) 
#d7t : day 7 - tepa
counts <- read.table("TEPA_data/Combined_D7_TEPA_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
TF <- CreateSeuratObject(counts = t(counts), project="TF", min.cells = 3, min.features = 500)
#d7Tm : day 7 - tepa myeloid 
counts <- read.table("TEPA_data/Combined_D7_TEPAMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
TM <- CreateSeuratObject(counts = t(counts), project = "TM",min.cells = 3, min.features = 500)

# Create list object to be merged in a large Seurat object
seuset <- list(CM = CM, TF = TF, TM = TM)

combined <- merge(
  x = CF,
  y = seuset,
  add.cell.ids = c("CF",'CM','TF', 'TM')
)

seuset <- combined
seuset$condition <- ifelse(test = seuset$orig.ident %in% c("CF", "CM"), yes = "Control", no = "Treatment")

# table(seuset$orig.ident)

save(seuset, file = "TEPA_results/00_rawcounts.rda")

#### 2 - Data pre-processing ####

load("TEPA_results/00_rawcounts.rda")

png("TEPA_plots/00_preFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/00_preFilterQC_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(seuset, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005)
dev.off()

# Log-normalisation
seuset <- NormalizeData(seuset) 

# We only carry out filtering for immune cells after separating the tumour


#### 3 - Dimensionality reduction ####

seuset <- FindVariableFeatures(seuset, 
                               selection.method = "vst",
                               nfeatures = 2000)
seuset <- ScaleData(seuset, verbose = FALSE)
seuset <- RunPCA(seuset, npcs = 100, verbose = FALSE) 

# Determine percent of variation associated with each PC
pct <- seuset@reductions$pca@stdev / sum(seuset@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=40) # Select 40 PCs to retain 62.76% of variability

seuset <- FindNeighbors(object = seuset, dims = 1:40, reduction = 'pca')
seuset <- RunUMAP(seuset, dims = 1:40, reduction = "pca", verbose = FALSE)

png("TEPA_plots/00_umapExplore.png", w = 4000, h = 2000, res = 300)
DimPlot(object = seuset, pt.size = 0.0005, reduction = 'umap', ncol = 2,
        group.by = c("orig.ident", "condition"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/00_umapCounts.png", w = 4000, h = 2000, res = 300)
FeaturePlot(seuset, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

#### 4 - Create a gene set uniquely identifying Neuroblastoma tumour ####

NB = c("Mycn", "Myc", "Dcx", "Phox2b", "Ccnd1", "Alk",
           "Trp53", "Kif18a","Ezh2", "Nf1", "Sdhb",
           "Vangl1","Pde6g","Nras", "Brip1","Brca1","Brca2","Ptpn11","Apc")

# https://doi.org/10.1038/nrdp.2016.78

seuset <- AddModuleScore(seuset, features = list(NB), name="NB")
names(seuset@meta.data)[grep("NB", names(seuset@meta.data))] <- "NB"

png("TEPA_plots/00_umapNB.png", h = 2000, w = 6000, res = 300)
FeaturePlot(seuset, ncol = 3, pt.size = 0.0005,
            features = c("NB", "Mycn", "Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

save(seuset, file = "TEPA_results/00_seusetNB.rda")

#### 5 - Subset only immune cells by using the module score for NB ####

table(seuset$orig.ident) # pre-subsetting

immune <- subset(seuset, NB > 0, invert = TRUE) # we remove very little amount of cells from myeloid samples, it looks really specific to the tumour
immune <- subset(immune, Mycn > 2 & Ptprc < 1, invert = TRUE) # immune cells post-filtering: we double check if we missed some Mycn+/Cd45- cell

table(immune$orig.ident) # post-subsetting

png("TEPA_plots/00_postFilterImmune_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/00_postFilterImmune_Mycn_Cd45.png", h = 3000, w = 4200, res = 300)
FeatureScatter(immune, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 0.0005)
dev.off()

png("TEPA_plots/00_umapNB_Immune.png", h = 2000, w = 4000, res = 300)
FeaturePlot(immune, ncol = 2, pt.size = 0.0005,
            features = c("Mycn", "Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

save(immune, file = "TEPA_results/00_immune.rda")




