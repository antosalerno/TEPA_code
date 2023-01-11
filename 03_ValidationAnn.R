#### Validate clustering annotation ####
## author: Antonietta Salerno
## date: 20/12/2022

library("Seurat")
library("ggplot2")
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
seuset_immune <- LoadH5Seurat("TEPA_results/02_immuneAnn.h5Seurat")
immune.markers <- read.csv("TEPA_results/02_DEA_clusterMarkers.csv")

DefaultAssay(seuset_immune) <- "RNA"

# Search all isoforms of gene of interest
grep(pattern = "Ly6", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)


# Add module score to the Seurat object for each cluster 
createSets <- function(markers = immune.markers,
                       obj = seuset_immune){
  Idents(seuset_immune) <- "scType"
  clusters = unique(Idents(seuset_immune))
  for(c in 1:length(clusters)){
    clusterDF = immune.markers[immune.markers$cluster == clusters[c],]
    clusterDF <- clusterDF[order(-clusterDF$avg_log2FC),]
    top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene)
    seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(top_genes), name=make.names(as.character(clusters[c])))
    names(seuset_immune@meta.data)[grep(make.names(as.character(clusters[c])), names(seuset_immune@meta.data))] <- as.character(clusters[c])
  }
  return(seuset_immune)
}

seuset_immune <- createSets()

# Plot results 

png("TEPA_plots/03_clustersAnnot.png", h = 3000, w = 4500, res = 200)
patchwork::wrap_plots(FeaturePlot(seuset_immune, ncol = 4, combine = TRUE,
            features = as.character(clusters), label = TRUE, repel = TRUE)) & theme_minimal() &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

# Get top marker genes by cluster
getTopMarkers <- function(df = immune.markers, geneNr){
  markers <- c()
  for(c in 1:length(clusters)){
    clusterDF = df[df$cluster == clusters[c],]
    clusterDF <- clusterDF[order(-clusterDF$avg_log2FC),]
    top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene, geneNr)
    markers <- append(markers, top_genes)
  }
  return(markers)
}
markers <- getTopMarkers(immune.markers, 5)

png("TEPA_plots/03_DotPlot.png", h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = unique(markers), # split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

SaveH5Seurat(seuset_immune, filename = "TEPA_results/03_seusetImmuneModule.h5Seurat", overwrite = TRUE)

