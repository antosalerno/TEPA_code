#### Validate clustering annotation ####
## author: Antonietta Salerno
## date: 20/12/2022

library("Seurat")
library("SeuratDisk")
library("ggplot2")
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
library("Nebulosa")
library("ggsignif")
library("ggpubr")

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadH5Seurat("TEPA_results/02_immuneAnn.h5Seurat")
immune.markers <- read.csv("TEPA_results/02_DEA_clusterMarkers.csv")

DefaultAssay(seuset_immune) <- "RNA"

# Search all isoforms of gene of interest
grep(pattern = "Prnp", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

#  Plot density plot with key markers for important populations  

png("TEPA_plots/03_densityPlatelets.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <-"RNA"
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

png("TEPA_plots/03_Mt2_split.png", h = 2000, w = 3500, res = 200)
plot_density(seuset_immune, "Mt2") + 
  facet_grid(.~seuset_immune$condition)
#plot_density(seuset_immune, features = c("Mt1", "Mt2"), joint = TRUE)
dev.off()

png("TEPA_plots/03_Mt1Mt2_violin.png", h = 2000, w = 3500, res = 200)
Idents(seuset_immune) <- "condition"
VlnPlot(seuset_immune, features = c("Mt2", "Mt1"), ncol = 2,pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_Ifngr_violin.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = "Ifngr1", split.by="condition", ncol = 1,pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_Ceacam.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = c("Ceacam1","Ceacam10"), 
        split.by="condition", ncol = 2, pt.size = 0.000005) + geom_boxplot()
dev.off()

png("TEPA_plots/03_NeutrophilsActive.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Fut4", "Itgb2", "Itgam"), 
        split.by="condition", ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_NeutroMacro.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Ly6g", "Ly6g5b", "Ly6g6c"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S03_Prnp_ImmuneCond.png", h = 2000, w = 3500, res = 400)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
VlnPlot(seuset_immune, features = c("Prnp"), 
        split.by = "condition", ncol = 1, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S03_Mt1_ImmuneCond.png", h = 2000, w = 3500, res = 400)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
fig <- VlnPlot(seuset_immune, features = c("Mt1"), 
               split.by = "condition", ncol = 1, pt.size = 0.000005)
fig +  ylim(0, 5) +
  #geom_signif(xmin = 11.75, xmax = 12.25, y_position = 3.5, annotations="*") +
  geom_signif(xmin = 0.75, xmax = 1.25, y_position = 3.75, annotations="***") +
  geom_signif(xmin = 4.75, xmax = 5.25, y_position = 4.5, annotations="***")
dev.off()

png("TEPA_plots/S03_lympho.png", h = 4000, w = 6500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Cd4", "Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Cd3g"), ncol = 3, pt.size = 0.000005)
dev.off()

#### Create a signature of copper-related genes ####






