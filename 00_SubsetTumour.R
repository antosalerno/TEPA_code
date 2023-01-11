###### Subset tumour cells from immune cells ##
## author: Antonietta Salerno
## date: 21/12/2022

library("Seurat")
library("SeuratObject")
library("ggplot2")
library(dplyr)
library(SeuratDisk)
library(SeuratData)
library(RColorBrewer)

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
seuset$sampleType <- ifelse(test = seuset$orig.ident %in% c("TM", "CM"), yes = "Myeloid", no = "Full")


# table(seuset$orig.ident)

# save(seuset, file = "TEPA_results/00_rawcounts.rda")
SaveH5Seurat(seuset, filename = "TEPA_results/00_rawcounts.h5Seurat", overwrite = TRUE)

#### 2 - Data pre-processing ####

# load("TEPA_results/00_rawcounts.rda")
seuset <- LoadH5Seurat("TEPA_results/00_rawcounts.h5Seurat")

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

png("TEPA_plots/00_umapMycnPtprc.png", h = 2000, w = 4000, res = 300)
FeaturePlot(seuset, ncol = 2, pt.size = 0.0005,
            features = c("Mycn", "Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

SaveH5Seurat(seuset, filename = "TEPA_results/00_seusetRed.h5Seurat", overwrite = TRUE)

#### 4 - Identify tumour cells ####

seuset <- LoadH5Seurat("TEPA_results/00_seusetRed.h5Seurat")

table(seuset$orig.ident) # pre-subsetting: CF: 9100  CM: 7377  TF: 9412 TM: 11174 


hist_dens <- function(x, breaks = "Scott", 
                      main = "Mycn expression per number of cells",
                      xlab = "Mycn", ylab = "Cell counts") {
  dens <- density(x, na.rm = T)
  raw_hist <- hist(x, breaks = breaks, plot = F)
  scale <- max(raw_hist$counts)/max(raw_hist$density)
  hist(x, breaks = breaks, prob = F, main = main, xlab = xlab, ylab = ylab)
  lines(list(x = dens$x, y = scale * dens$y), col = "blue", lwd = 2)
}

png("TEPA_plots/00_MycnDistribution.png", h = 3000, w = 4200, res = 300)
X <- colSums(seuset["Mycn",])
hist_dens(X, breaks = 100)
dev.off()

tumorCells1 <- subset(seuset[,seuset$sampleType == "Full"], Mycn > 0 & Ptprc == 0)
tumorCells2 <- subset(tumorCells[,tumorCells$sampleType == "Full"], 
                     Cd3d > 0 | Cd3e > 0 | Cd3g > 0 | Cd8a > 0 | Cd8b1 > 0, invert = TRUE) # 13560 cells: 169 lymphocytic cells rescued
immuneCells <- setdiff(Cells(tumorCells1), Cells(tumorCells2))
seuset$class <- ifelse(test = colnames(seuset) %in% append(Cells(tumorCells1), Cells(tumorCells2)), "tumorCells", "noClass")

png("TEPA_plots/00_tumorCells.png", h = 2000, w = 4000, res = 300)
FeaturePlot(seuset[,seuset$class == "tumorCells"], ncol = 2, pt.size = 0.005,
            features = c("Mycn","Ptprc"), label = TRUE, repel = TRUE) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

immuneCells <- append(immuneCells, Cells(subset(seuset[,seuset$sampleType == "Full"], Ptprc > 0 & Mycn == 0)))
immuneCells <- append(immuneCells, Cells(subset(seuset[,seuset$sampleType == "Myeloid"])))
seuset$class <- ifelse(test = colnames(seuset) %in% immuneCells, "immuneCells", seuset$class)
seuset$class <- ifelse(test = colnames(seuset) %in% Cells(append(tumorCells,immuneCells)), "noClass", seuset$class)

pt <- table(Idents(seuset), seuset$class)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

png(paste0("TEPA_plots/00_classMemberMycnPtprc.png"), w=2500,h=2500, res=300)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Class") +
  ylab("Cell counts") +
  ggtitle("Class membership of cell counts") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

png("TEPA_plots/00_noClass_MycnPtprc.png", h = 3000, w = 4200, res = 300)
FeatureScatter(seuset[,seuset$class == "noClass"], feature1 = "Mycn", feature2 = "Ptprc", pt.size = 1)
dev.off()

SaveH5Seurat(seuset, filename = "TEPA_results/00_seusetClass.h5Seurat", overwrite = TRUE)

### 5 -  Select only immune cells ####

seuset <- LoadH5Seurat("TEPA_results/00_seusetClass.h5Seurat")

immune <- seuset[,seuset$class == "immuneCells"]

table(immune$orig.ident) 
# post-subsetting all samples: CF: 2147 CM: 7377  TF: 1484 TM: 11174

png("TEPA_plots/00_immuneCells_MycnPtprc.png", h = 3000, w = 4200, res = 300)
FeatureScatter(immune, feature1 = "Mycn", feature2 = "Ptprc", pt.size = 1)
dev.off()

png("TEPA_plots/00_postFilterImmune_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
dev.off()

# save(immune, file = "TEPA_results/00_immune.rda")
SaveH5Seurat(immune, filename = "TEPA_results/00_immune.h5Seurat", overwrite = TRUE)




