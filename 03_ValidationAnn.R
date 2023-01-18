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
grep(pattern = "Cx3cr1", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)


### 1 - Add module score to the Seurat object for each cluster ####
createSets <- function(markers = immune.markers,
                       obj = seuset_immune){
  Idents(seuset_immune) <- "scType"
  clusters = unique(Idents(seuset_immune))
  for(c in 1:length(clusters)){
    clusterDF <- immune.markers[immune.markers$cluster == clusters[c],]
    clusterDF <- clusterDF[clusterDF$p_val_adj < 0.05,]
    clusterDF <- clusterDF[order(-c(clusterDF$avg_log2FC)),]
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

#### 2 - Identify most important markers per cluster ####

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

### 3 - Plot density plot with key markers for important populations ####

png("TEPA_plots/03_densityMacrophagesJoint.png", h = 3000, w = 3500, res = 100)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = macro_markers_mouse, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMacrophagesClassic.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
classical <- c("Ccr2", "Ccl2")
plot_density(seuset_immune, features = classical, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMacrophagesNonClassic.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
non_classical <- c("Cx3cr1", "Spn")
plot_density(seuset_immune, features = non_classical, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityNeutrophils.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
neutro <- c("Cxcr1","Cxcr2", "Ccr1", "Ccrl2", "Ly6g", "Fcgr3")
plot_density(seuset_immune, features = neutro, joint = TRUE)
dev.off() # Cxcl8 is the most important markers for neutrophils (Il8 production) but no homolog in mouse -> analogs are Cxcl1, Cxcl2, Cxcl5

png("TEPA_plots/03_densityILCs.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
ilc <- c("Il7", "Il7r")
plot_density(seuset_immune, features = ilc, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMacro.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
macro <- c("Adgre1", "Fcgr1") # F4/80 (Adgre1) and Cd64 the most important markers for murine macrophages -> Adgre1 expressed also in eosinophils in mouse
plot_density(seuset_immune, features = macro, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMonocytes.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
mono <- c("Itgam")
plot_density(seuset_immune, features = mono, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityLeuko.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
leuko <- c("Cadm1", "Cadm2", "Cadm3", "Cadm4")
plot_density(seuset_immune, features = leuko, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityGranu.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
granu <- c("Ly6g", "Ly6c1", "Ly6c2")
plot_density(seuset_immune, features = granu, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityHST.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
hst <- c("Cd34")
plot_density(seuset_immune, features = hst, joint = TRUE)
dev.off()

png("TEPA_plots/03_densityPlatelets.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Itga2b", "Itgb3"), joint = TRUE)
dev.off()

png("TEPA_plots/03_densityErythro.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = "Gypa")
dev.off()

png("TEPA_plots/03_densityTreg.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il2ra","Foxp3"), joint = TRUE)
dev.off() # Treg higher in tumor and in inflammation: cd25 and foxp3

# Classical monocytes and Nks: Ccr2 important for neutrophils mobilitisation and recruitment to metastatic sites

png("TEPA_plots/03_densityEosino.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il5ra", "Ccr3", "Adgre1","Ccl6"), joint = TRUE)
dev.off()

png("TEPA_plots/03_densityMemory.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune)<-"RNA"
plot_density(seuset_immune, features = c("Il7r", "Ctla4"), joint = TRUE)
dev.off()

SaveH5Seurat(seuset_immune, filename = "TEPA_results/03_seusetImmuneModule.h5Seurat", overwrite = TRUE)



