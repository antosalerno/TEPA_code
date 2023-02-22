###### Dimensionality reduction and clustering of immune cells with higher resolution##
## author: Antonietta Salerno
## date: 13/01/2022

library("Seurat")
library("ggplot2")
library("writexl")
library(openxlsx)
library(HGNChelper)
library(dplyr)

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")

#### 1 - Dimensionality reduction ####

immune.combined <- LoadH5Seurat("TEPA_results/S01_immuneInt.h5Seurat")

# Scaling
seuset_immune <- ScaleData(immune.combined, verbose = FALSE)
seuset_immune <- RunPCA(seuset_immune, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_immune@reductions$pca@stdev / sum(seuset_immune@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 70.26% of variability

seuset_immune <- FindNeighbors(object = seuset_immune, graph.name = "clust", dims = 1:60, reduction = 'pca')

### 2 - Clustering ####

seuset_immune <- FindClusters(object = seuset_immune, graph.name = "clust", resolution = 0.3) 
seuset_immune <- RunUMAP(seuset_immune, dims = 1:60, reduction = "pca", verbose = FALSE)

png("TEPA_plots/S02_umapExplore.png", w = 6000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("orig.ident", "condition", "seurat_clusters"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("TEPA_plots/S02_umapCounts.png", w = 4000, h = 2000, res = 300)
FeaturePlot(seuset_immune, features = c("nCount_RNA", "nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
dev.off()

png("TEPA_plots/S02_umapClust.png", w = 4000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("seurat_clusters"), split.by= "condition",label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## Sub-clustering ##

Idents(seuset_immune) <- "clust_res.0.3"
seuset_immune <- FindSubCluster(seuset_immune, "3", "clust", subcluster.name = "seurat_clusters",  resolution = 0.2, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "3_0", "clust", subcluster.name = "seurat_clusters",  resolution = 0.3, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "2", "clust", subcluster.name = "seurat_clusters",  resolution = 0.35, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"

seuset_immune <- FindSubCluster(seuset_immune, "2_2", "clust", subcluster.name = "seurat_clusters",  resolution = 0.2, algorithm = 1)
DimPlot(seuset_immune, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_immune) <- "seurat_clusters"


### 4 - Annotate subclusters
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system"

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = seuset_immune[["integrated"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

####  4 - Get scType scores by cluster ####

cL_results = do.call("rbind", 
                     lapply(unique(seuset_immune@meta.data$seurat_clusters), 
                            function(cl){
                              es.max.cl = sort(rowSums(es.max[,rownames(seuset_immune@meta.data[seuset_immune@meta.data$seurat_clusters==cl, ])]),decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuset_immune@meta.data$seurat_clusters==cl)), 10)
                            }))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

# Manual curation of scType annotation
sctype_scores[3, "type"] <- "Cd4+ Naive T cells"
sctype_scores[4, "type"] <- "Cd4+ Memory T cells" #"Cd4+ Naive T cells"
sctype_scores[5, "type"] <- "Cd8+ Naive T cells" #"Cd4+ Memory T cells"
sctype_scores[6, "type"] <- "Cd8+ NkT-like cells"
sctype_scores[7, "type"] <- "DN Regulatory T cells" #"Cd8+ Naive T cells"
sctype_scores[8, "type"] <- "Cd4+ Naive T cells" #"Cd8+ NkT-like cells"
sctype_scores[10, "type"] <- "Macrophages" # or macrophages
sctype_scores[11, "type"] <- "Macrophages" # or macrophages
sctype_scores[12, "type"] <- "HLA-expressing cells" #"Dendritic cells" # or DCs
sctype_scores[13, "type"] <- "Dendritic cells" #"B cells"
sctype_scores[15, "type"] <- "B cells"

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(seuset_immune) <- seuset_immune@meta.data$scType

# We decide to remove the cluster HLA-expressing cells since it is probably ambient RNA -> mixture of B cells and macrophage markers
seuset_immune <- seuset_immune[,!seuset_immune$scType == "HLA-expressing cells"]
levels(Idents(seuset_immune)) # now 14 clusters rather than 15

options(ggrepel.max.overlaps = Inf)
png("TEPA_plots/S02_umapAnn.png", w = 3000, h = 3000, res = 300)
p <- DimPlot(object = seuset_immune, pt.size = 0.5, reduction = 'umap', ncol = 1,
             group.by = "scType", label = FALSE) +
  ggtitle("Cell types in NB Control and Treated samples") +
  theme(plot.title = element_text(hjust = 0.5)) 
LabelClusters(p, id = "scType", size = 5, repel = T,  box.padding = 1)
dev.off()

png("TEPA_plots/S02_umapCondAnn.png", w = 4000, h = 2000, res = 300)
p <- DimPlot(object = seuset_immune, pt.size = 0.05, reduction = 'umap', ncol = 2,
        group.by = "scType", split.by= "condition", label = FALSE) +
  ggtitle("Cell types in NB Control and Treated samples") +
  theme(plot.title = element_text(hjust = 0.5))
LabelClusters(p, id = "scType", size = 3, repel = T, box.padding = 1)
dev.off()

### Bar Plot with proportions
pt <- table(Idents(seuset_immune), seuset_immune$condition)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
#getPalette = colorRampPalette(brewer.pal(11, "PiYG"))
#colorCount = length(unique(pt$Var1))

png(paste0("TEPA_plots/S02_condAnnClusterFreq.png"), w=2500,h=2500, res=300)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Condition") +
  ylab("Cell type") +
  ggtitle("Cell types in TEPA vs Control") +
  #scale_fill_manual(values = getPalette(colorCount)) +
  theme(legend.title = element_blank())
dev.off()

SaveH5Seurat(seuset_immune, filename = "TEPA_results/S02_immuneAnn.h5Seurat", overwrite = TRUE)

