# Evaluate tumor heterogeneity #
## author: Antonietta Salerno
## date: 07/02/2023

library("Seurat")
library("SeuratDisk")
library("MAST")
library("writexl")
library(openxlsx)
library("EnhancedVolcano")
library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library(readr)
library(stringr)
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
library(ggplot2)
library(ggcharts)


setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")
tumor <- LoadH5Seurat("TEPA_results/S00_tumor.h5Seurat") #13560 cells and 21388 genes

#### 1 - QC and filtering of tumor cells ####

tumor <- NormalizeData(tumor) 
# Add percent mito data in seurat seuset_tumorect
tumor[["percent.mt"]] <- PercentageFeatureSet(tumor, pattern = "^mt.")
# Add percent ribo genes
tumor[["percent.ribo"]] <- PercentageFeatureSet(tumor, pattern = "^Rp.")

# Visualize the distribution of mitochondrial gene expression detected per cell
png("TEPA_plots/S05_tumorPreFilterQC_mito.png", h = 3000, w = 4200, res = 300)
tumor@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 25)
dev.off()

png("TEPA_plots/S05_tumorPreFilterQC_ribo.png", h = 3000, w = 4200, res = 300)
tumor@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.ribo, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 15)
dev.off()

tumor <- subset(tumor, subset = nFeature_RNA > 200 &
                  percent.mt < 25 & percent.ribo < 15) # 13535 cells

png("TEPA_plots/S05_tumorPostFilterQC_violin.png", h = 3000, w = 4200, res = 300)
VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "Mycn"), ncol = 3, pt.size = 0.000005)
# a little difference in the distribution -> batch effect 
dev.off()

png("TEPA_plots/S05_tumorPostFilterQC_scatter.png", h = 3000, w = 4200, res = 300)
plot1 <- FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.0005)
plot2 <- FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.0005)
plot1 + plot2 # We can spot two separate data clouds because of the bimodal distribution of the data: control vs treated cells
dev.off()


### Batch effect correction ###

samples.list <- Splitseuset_tumorect(tumor, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(seuset_tumorect.list = samples.list)
tumor.anchors <- FindIntegrationAnchors(seuset_tumorect.list = samples.list,anchor.features = features)
tumor.combined <- IntegrateData(anchorset = tumor.anchors)

#### 2 - Dimensionality reduction ####

# Scaling
seuset_tumor <- ScaleData(tumor.combined, verbose = FALSE)
seuset_tumor <- RunPCA(seuset_tumor, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_tumor@reductions$pca@stdev / sum(seuset_tumor@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=80) # Select 80 PCs to retain 86.09 % of variability

seuset_tumor <- FindNeighbors(seuset_tumorect = seuset_tumor, graph.name = "clust", dims = 1:80, reduction = 'pca')
seuset_tumor <- FindClusters(seuset_tumorect = seuset_tumor, graph.name = "clust", resolution = 0.3, algorithm = 1) # original plot is 0.7
seuset_tumor <- RunUMAP(seuset_tumor, dims = 1:80, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S05_tumorUmapClust.png", w = 4000, h = 4000, res = 350)
DimPlot(seuset_tumorect = seuset_tumor, pt.size = 0.5, reduction = 'umap', ncol = 1,
             group.by = c("seurat_clusters", "condition"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_tumor@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5)) 
dev.off()

SaveH5Seurat(seuset_tumor, filename = "TEPA_results/S05_seusetTumorClu.h5Seurat", overwrite = TRUE)

#### 3 - Clustering annotation ####

### 3.1 Inter-cluster DEA: get marker genes ###

DefaultAssay(seuset_tumor) <- "RNA"

save = "S05_tumorMarkers"
Idents(seuset_tumor) <- "seurat_clusters"

## A - Find markers for every cluster compared to all remaining cells
tumor.markers <- FindAllMarkers(seuset_tumor, 
                                only.pos = FALSE, 
                                min.pct = 0.5, 
                                min.diff.pct = 0.2,
                                logfc.threshold = 0.3, 
                                test.use="MAST",
                                latent.vars="orig.ident")

write.csv(tumor.markers, file=paste0("TEPA_results/", save,".csv"))

# Save results in different excel sheets 
clusters = levels(Idents(seuset_tumor))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = tumor.markers[tumor.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file=paste0("TEPA_results/",save,".xlsx"), overwrite = TRUE)

# B - Add module score to the Seurat seuset_tumorect for each cluster ###
seuset_tumor <- createSets(markers = tumor.markers, 
                           seuset_tumor=seuset_tumor, id = "seurat_clusters")

png("TEPA_plots/S05_tumorClustersAnn.png", h = 3000, w = 4500, res = 300)
clusters = unique(Idents(seuset_tumor))
patchwork::wrap_plots(FeaturePlot(seuset_tumor, ncol = 2, combine = TRUE, pt.size = 2, 
                                  features = as.character(clusters), label = TRUE, repel = TRUE)) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off() 

## C - Identify most important markers per cluster ###

markers <- getTopMarkers(tumor.markers, 5)

png("TEPA_plots/S05_tumorDotPlot.png", h = 2000, w = 2500, res = 300)
DotPlot(seuset_tumorect = seuset_tumor, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

## D - Plot Volcano of each cluster vs all the others: 
save = "S05_tumorMarkers"
plotVolcano(clusters, res = tumor.markers, type = "markers", input = "tumor", log2FC = 0.5, save = save)

#### 3.2 - Intra-cluster DEA with annotated dataset - Treatment vs Control ####
Idents(seuset_tumor) <- "seurat_clusters"
seuset_tumor$clusters.tepa <- paste(Idents(seuset_tumor), seuset_tumor$condition, sep = "_")
Idents(seuset_tumor) <- "clusters.tepa"
DefaultAssay(seuset_tumor) <- "RNA"

save = "S05_tumorCond_"
sheets <- list()
for (cluster in unique(seuset_tumor$seurat_clusters)){
  try({
    ident1 <- paste0(cluster,"_Treatment")
    ident2 <- paste0(cluster,"_Control")
    message(paste0("DEA control vs treatment for cluster: ", cluster))
    condition.diffgenes <- FindMarkers(seuset_tumor, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       #logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE,
                                       #latent.vars="orig.ident",
                                       #min.cells.feature = 1, min.cells.group = 1, 
                                       test.use="MAST")
    condition.diffgenes$p_val_adj = p.adjust(condition.diffgenes$p_val, method='BH')
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/", save,"DEA_",sub(" ", "_", cluster),".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, paste0("TEPA_results/", save,"DEA.xlsx"), rowNames=TRUE)

# Plot Volcano DEA by condition

Idents(seuset_tumor) <- "seurat_clusters"
clusters = unique(Idents(seuset_tumor))

save = "S05_tumor_"
plotVolcano(clusters, log2FC = 0.5, save = save, input = "tumor")

#### 3.3 DEA Treatment vs Control bulk dataset ####
counts <- GetAssayData(seuset_tumor, assay = "RNA")
counts <- counts[-(which(rownames(counts) == "Xist")),]
seuset_tumor <- subset(seuset_tumor, features = rownames(counts))

Idents(seuset_tumor) <- "condition"
res <- FindMarkers(seuset_tumor, 
                                ident.1 = "Treatment", ident.2 = "Control",
                                only.pos = FALSE, verbose = FALSE,
                                latent.vars="orig.ident",
                                test.use="MAST")
res$p_val_adj = p.adjust(res$p_val, method='BH')
write.csv(res, file=paste0("TEPA_results/S05_tumorBulkDEA.csv"))

log2FC = 0.5
save = "S05_tumorBulk_"

# volcano bulk doesn't work

p <- EnhancedVolcano(res, subtitle = "",
                     #selectLab = markers,
                     lab = rownames(res),
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     xlim = c(-2.5, 2.5),
                     title = "DEA bulk TEPA vs Control",
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
                                      "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -log2FC&res$p_val_adj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
ggsave(p, file=paste0("TEPA_plots/", save, "DEA.png"), width = 30, height = 25, units = "cm")

bulkDEAgenes <- head(rownames(res %>% arrange(avg_log2FC)), 30)
png(paste0("TEPA_plots/", save, "Heatmap.png"), h = 3000, w = 4500, res = 200)
DoHeatmap(seuset_tumorect = seuset_tumor, 
          features = bulkDEAgenes, label = TRUE)
dev.off()


#### 4 - Gene Set Enrichment Analysis ####

sets1 <- read.gmt("TEPA_data/mh.all.v2022.1.Mm.symbols.gmt") # Mouse hallmark
sets2 <- read.gmt("TEPA_data/REACTOME_NEUTROPHIL_DEGRANULATION.v2022.1.Mm.gmt") # The only Reactome pathway we're interested in

sets1$term <- as.character(sets1$term)
sets2$term <- as.character(sets2$term)

sets1 <- sets1 %>% split(x = .$gene, f = .$term)
sets2 <- sets2 %>% split(x = .$gene, f = .$term)

fgsea_sets <- append(sets1, sets2)
fgsea_sets <- append(fgsea_sets, custom)

Idents(seuset_tumor) <- "seurat_clusters"
clusters = unique(Idents(seuset_tumor))
save = "S05_tumorEnrichment_"
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "tumor")

#### 5 - Look at expression of specific genes ####

# Search all isoforms of gene of interest
grep(pattern = "Ifi2", 
     x = rownames(x = seuset_tumor@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

png("TEPA_plots/S05_tumorMt1.png", h = 2000, w = 2000, res = 200)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "Mt1", ncol = 1, pt.size = 0.000005) + 
  ylim(-1,5) +
  geom_signif(xmin = 1, xmax = 2, y_position = 4, annotations="***")
dev.off()

png("TEPA_plots/S05_tumorMycn.png", h = 2000, w = 2000, res = 200)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "Mycn", ncol = 1, pt.size = 0.000005) + 
  ylim(0,6) +
  geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***")
dev.off()

png("TEPA_plots/S05_tumorMycn.png", h = 2000, w = 2000, res = 200)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "Mycn", ncol = 1, pt.size = 0.000005) + 
  ylim(0,6) +
  geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***")
dev.off()

png("TEPA_plots/S05_tumorPdl1.png", h = 2000, w = 2000, res = 200)
Idents(seuset_tumor) <- "condition"
VlnPlot(seuset_tumor, features =  "Cd274", ncol = 1, pt.size = 0.000005) + 
  ylim(0,6) +
  geom_signif(xmin = 1, xmax = 2, y_position = 5.25, annotations="***")
dev.off()

#### 6 - Investigate whether the tumor is adrenergic or mesenchymal ####
seuset_tumor <- LoadH5Seurat("TEPA_results/S05_seusetTumorClu.h5Seurat")

MES <- c("Elk4", "Creg1", "Dcaf6", "Id1", "Smad3", "Six4", "Six1", "Maml2", "Notch2",
          "Cbfb", "Ifi203","Ifi204", "Ifi205", "Ifi206", "Ifi207", "Ifi208", "Ifi209",
         "Ifi211","Ifi214", "Zfp217", "Egr3", "Zfp36l1", "Wwtr1", "Prrx1", "Sox9",
          "Meox1", "Meox2", "Aebp1", "Col1a1", "Col1a2", "Fn1", "Vim", "Snai2", "Prom1",
          "Prrx1", "Yap1", "Hey1", "Egr3", "Fosl1", "Fosl2", "Tbx18", "Runx1", "Runx2",
         "Mef2d", "Irf1", "Irf2", "Irf3", "Fli1", "Glis3", "Maff", "Bhlhe41", "Nr3c1")
ADRN <- c("Zfp536", "Phox2a", "Hand1", "Ascl1", "Klf13", "Sox11", "Gata2", "Gata3", "Klf7",
          "Eya1", "Tfap2b", "Isl1", "Hey1", "Six3", "Dach1", "Phox2b", "Pbx3", "Satb1",
          "Chga", "Chgb", "Dbh", "Dlk1", "Tfap2b", "Six3")

# https://www.cell.com/cell-reports/pdf/S2211-1247(22)01296-7.pdf
# https://www.nature.com/articles/ng.3899
# https://www.nature.com/articles/ng.3899
# https://www.nature.com/articles/s41467-019-09470-w#:~:text=Mesenchymal%2Dtype%20(MES)%20neuroblastoma,a%20high%20transcriptional%20plasticity1.

# Run gsea
custom_tum <- list(MES,ADRN)
names(custom_tum) <- c("MESENCHYMAL TUMOR", "ADRENERGIC TUMOR")

Idents(seuset_tumor) <- "condition"
clusters = unique(Idents(seuset_tumor))
save = "S05_tumorAdreMese_"
gseaRES(clusters, fgsea_sets = custom_tum, save = save, 
        minSize = 5, input = "tumor")

# Add module score to the Seurat seuset_tumor for each signature ###

seuset_tumor <- AddModuleScore(seuset_tumor, assay = "RNA", features = list(ADRN), name=make.names("Adrenergic"))
names(seuset_tumor@meta.data)[grep(make.names("Adrenergic"), names(seuset_tumor@meta.data))] <- "Adrenergic"

seuset_tumor <- AddModuleScore(seuset_tumor, assay = "RNA", features = list(MES), name=make.names("Mesenchymal"))
names(seuset_tumor@meta.data)[grep(make.names("Mesenchymal"), names(seuset_tumor@meta.data))] <- "Mesenchymal"

library(RColorBrewer)
png("TEPA_plots/S05_tumorAdreMese_FeaturePlot.png", h = 3000, w = 4500, res = 200)
patchwork::wrap_plots(FeaturePlot(seuset_tumor, ncol = 2, combine = TRUE, pt.size = 5,
                                  features = c("Adrenergic", "Mesenchymal"), label = F,
                                  repel = TRUE)) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

Idents(seuset_tumor) <- "condition"
png("TEPA_plots/S05_tumorAdreMese_Dotplot.png", h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_tumor, features = unique(c(MES,ADRN)), 
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

Idents(seuset_tumor) <- "condition"
png("TEPA_plots/S05_tumorAdreMese_ScatterPlot.png", h = 3000, w = 4200, res = 300)
FeatureScatter(seuset_tumor, feature1 = "Adrenergic", feature2 = "Mesenchymal", pt.size = 0.0005)
dev.off()





