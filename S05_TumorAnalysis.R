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
# Add percent mito data in seurat object
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

samples.list <- SplitObject(tumor, split.by = "condition")

# Normalize and identify variable features for each dataset independently (Treatment vs Control)
samples.list <- lapply(X = samples.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = samples.list)
tumor.anchors <- FindIntegrationAnchors(object.list = samples.list,anchor.features = features)
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

seuset_tumor <- FindNeighbors(object = seuset_tumor, graph.name = "clust", dims = 1:80, reduction = 'pca')
seuset_tumor <- FindClusters(object = seuset_tumor, graph.name = "clust", resolution = 0.3, algorithm = 1) # original plot is 0.7
seuset_tumor <- RunUMAP(seuset_tumor, dims = 1:80, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S05_tumorUmapClust.png", w = 4000, h = 4000, res = 350)
DimPlot(object = seuset_tumor, pt.size = 0.5, reduction = 'umap', ncol = 1,
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

# B - Add module score to the Seurat object for each cluster ###
seuset_tumor <- createSets(markers = tumor.markers, 
                           obj=seuset_tumor, id = "seurat_clusters")

png("TEPA_plots/S05_tumorClustersAnn.png", h = 3000, w = 4500, res = 300)
clusters = unique(Idents(seuset_tumor))
patchwork::wrap_plots(FeaturePlot(seuset_tumor, ncol = 2, combine = TRUE, pt.size = 2, 
                                  features = as.character(clusters), label = TRUE, repel = TRUE)) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off() 

## C - Identify most important markers per cluster ###

markers <- getTopMarkers(tumor.markers, 5)

png("TEPA_plots/S05_tumorDotPlot.png", h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_tumor, features = unique(markers), # split.by = "condition",
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
DoHeatmap(object = seuset_tumor, 
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

N1 <- c("Icam1", "Cxcl10",  "Fas", "Tnf", "Ifng", "Ccl3", "Cxcl9", "Il12a", "Tnfaip2", "Cebpb",
        "Hif1a", "Tnfaip8", "Tnfaip1", "Tnfaip3", "Tnfaip6", "Tnfaip8l1", "Tnfaip8l2", 
        "Tnfaip8l3", "Ifit3", "Ifitm1", "Ifitm6", "Ifit1bl1", "Ifit2", "Ifit3", "Ifit3b", "Ifitm10", "Ifitm2", 
        "Ifitm3", "Ifitm5", "Ifitm6", "Ifitm7", "Ifit1bl2", "Ifi27l2a", "Il6ra", "Il15",
        "Il12b", "Csf2", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl4", "Cxcl11", "Cxcl1", "Ccl7", "Il18", "Il18bp",
        "Isg15", "Isg20", "Isg20l2")

N2 <- c("Il1rn","Tgfb1", "Tgfb1i1","Tgfb2", "Tgfb3", "Tgfbi", "Ccl2", "Mpo", "Cxcr4", 
        "Ccl5", "Slc27a2", "Ccl17", "Cxcl14", "Arg1", "Stat3", "Irf8", "Il17a", "Il17b","Il17c", "Il17d", "Il17f")

# https://www.sciencedirect.com/science/article/pii/S1535610809002153?via%3Dihub
# https://www.sciencedirect.com/science/article/pii/S1471490605002425
# https://ashpublications.org/blood/article/133/20/2159/273823/Neutrophil-plasticity-in-the-tumor

M1 <- c("Il12a", "Il12b", "Tnf", "Cxcl9", "Cxcl10", "Tlr2", "Tlr4", "Fcgr3", "Fcgr4", "Fcgr1",
        "Fcgr2b", "Cd80", "Cd86", "Il12a", "Il12b", "Il6", "Il1a", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
        "Cxcl8", "Cxcl9", "Cxcl10", "Cxcl11", "Ccr7")
M2 <- c("Il10", "Il1ra", "Ccl22", "Ccl24",  "Il4", "Il13", "Ccl16", 
        "Fcer2a", "Cd163", "Cxcr1", "Cxcr2", "Ccr2", "Arg1", "Il1rn")

# https://www.sciencedirect.com/science/article/pii/S1471490602023025

# Create gene sets list
custom <- list(N1, N2, M1, M2)
names(custom) <- c("N1 ANTI-TUMOR NEUTROPHILS", "N2 PRO-TUMOR NEUTROPHILS",
                   "M1 ANTI-TUMOR MACROPHAGES", "M2 PRO-TUMOR MACROPHAGES")

fgsea_sets <- append(fgsea_sets, custom)

Idents(seuset_tumor) <- "seurat_clusters"
clusters = unique(Idents(seuset_tumor))
save = "S05_tumorEnrichment_"
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save, input = "tumor")

#### 5 - Look at expression of specific genes ####

# Search all isoforms of gene of interest
grep(pattern = "B7", 
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



