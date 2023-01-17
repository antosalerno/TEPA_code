###### Dimensionality reduction and clustering of immune cells with higher resolution##
## author: Antonietta Salerno
## date: 13/01/2022

library("Seurat")
library("ggplot2")
library("writexl")
library(openxlsx)
library(dplyr)

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
# load("TEPA_results/00_immune.rda")
immune <- LoadH5Seurat("TEPA_results/00_immune.h5Seurat")

#### 1 - Dimensionality reduction ####

#load("TEPA_results/01_immuneInt.rda")
immune.combined <- LoadH5Seurat("TEPA_results/01_immuneInt.h5Seurat")

# Scaling
seuset_immune <- ScaleData(immune.combined, verbose = FALSE)
seuset_immune <- RunPCA(seuset_immune, npcs = 100, verbose = FALSE, assay = "integrated") 

# Determine percent of variation associated with each PC
pct <- seuset_immune@reductions$pca@stdev / sum(seuset_immune@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
head(cum, n=60) # Select 60 PCs to retain 70.37% of variability

seuset_immune <- FindNeighbors(object = seuset_immune, graph.name = "clust", dims = 1:60, reduction = 'pca')

### 2 - Clustering ####

seuset_immune <- FindClusters(object = seuset_immune, graph.name = "clust", resolution = 0.3) 
seuset_immune <- RunUMAP(seuset_immune, dims = 1:60, reduction = "pca", verbose = FALSE)

png("TEPA_plots/01_umapClust2.png", w = 4000, h = 2000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 1,
        group.by = c("seurat_clusters"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### 3 - Find macrophages by sub-clustering Neutrophils and DCs ####
Idents(seuset_immune) <- "clust_res.0.3"
seuset_immune <- FindSubCluster(seuset_immune, "3", "clust", subcluster.name = "seurat_clusters",  resolution = 0.2, algorithm = 1)
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
sctype_scores[3, "type"] <- "Cd4+ Naive T-cells"
sctype_scores[4, "type"] <- "Cd4+ Naive T-cells"
sctype_scores[5, "type"] <- "Cd4+ Memory T-cells"
sctype_scores[6, "type"] <- "Cd4+ T-regulatory cells"
sctype_scores[7, "type"] <- "Cd8+ Naive T-cells"
sctype_scores[8, "type"] <- "Cd8+ NkT-like cells"
sctype_scores[10, "type"] <- "Macrophages" # or macrophages
sctype_scores[11, "type"] <- "Dendritic cells" # or DCs
sctype_scores[13, "type"] <- "B-cells"

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(seuset_immune) <- seuset_immune@meta.data$scType

png("TEPA_plots/02_umapClust2.png", w = 4000, h = 4000, res = 350)
p <- DimPlot(object = seuset_immune, pt.size = 0.5, reduction = 'umap', ncol = 1,
        group.by = c("scType"), label = FALSE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5)) 

LabelClusters(p, id = "scType", size = 5, repel = T,  box.padding = 1)
dev.off()


### 5 - Search marker genes per cluster ####
de.genes <- FindAllMarkers(seuset_immune, 
                           only.pos = FALSE, 
                           min.pct = 0.5, 
                           min.diff.pct = 0.2,
                           logfc.threshold = 0.5)

clusters = unique(Idents(seuset_immune))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = immune.markers[immune.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/02_DEA_clusterMarkers2.xlsx", overwrite = TRUE)


