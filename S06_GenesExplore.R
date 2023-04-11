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
seuset_immune <- LoadH5Seurat("TEPA_results/S03_immuneDiff.h5Seurat")
immune.markers <- read.csv("TEPA_results/02_DEA_clusterMarkers.csv")

DefaultAssay(seuset_immune) <- "RNA"

# Search all isoforms of gene of interest
grep(pattern = "Il2", 
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
  geom_signif(xmin = 0.75, xmax = 1.25, y_position = 3.75, annotations="***") +
  geom_signif(xmin = 4.75, xmax = 5.25, y_position = 4.5, annotations="***")
dev.off()

png("TEPA_plots/S03_lympho.png", h = 4000, w = 6500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Cd4", "Cd8a", "Cd8b1", "Cd3d", "Cd3e", "Cd3g"), ncol = 3, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S06_LymphoMyelo.png", h = 2000, w = 3500, res = 200)
DefaultAssay(seuset_immune) <- "RNA"
Idents(seuset_immune) <- "scType"
FeaturePlot(seuset_immune, features = c("Cd3e", "Itgam"),
            ncol = 2, pt.size = 0.000005)
dev.off()


#### See difference in proportions of lymphoid and myeloid cells tumor vs immune cells ####
seuset_tumor <- LoadH5Seurat("TEPA_results/S05_seusetTumorClu.h5Seurat")

seuset_full <- merge(seuset_immune, y = seuset_tumor, 
                       add.cell.ids = c("immune", "tumor"), 
                       project = "singleCell")

seuset_full@assays$RNA@scale.data <- scale(seuset_full@assays$RNA@data, scale = TRUE)

SaveH5Seurat(seuset_full, filename = "TEPA_results/S06_seusetFull.h5Seurat", overwrite = TRUE)


png("TEPA_plots/S06_LymphoMyeloCompare.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "class"
VlnPlot(seuset_full, features = c("Cd3e", "Itgam"), split.by = "condition",
            ncol = 2, pt.size = 0.000005)
dev.off()


#### Create a signature of copper-related genes ####

copper_genes <- c("Slc31a1", "Atp7a", "Atp7b", "Sco1", "Cox11", "Steap3", "Commd1", "Mtf1", "Mtf2", "Sp1", "Sod1",
                  "Sod2", "Steap4", "Atox1", "Ccs", "Mt1", "Mt2", "Mt3", "Fdx1", "Lias", "Lipt1", "Dld", "Dlat",
                  "Pdha1", "Pdhb",  "Gls", "Cdkn2a")

copper_NBL <- c("Alb","Akt1","Ap1s1","Ang","Cdk1","Jun",
                "Map1lc3a","Mapt","Mtf2","F8","Tmprss6")

# Search all isoforms of gene of interest
grep(pattern = "Ang", 
     x = rownames(x = seuset_full@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)

Idents(seuset_immune) <- "scType"
png("TEPA_plots/S06_immuneGenesCopper.png", h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_immune, features = copper_genes, split.by = "condition",
        scale=TRUE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()


Idents(seuset_full) <- "class"
png("TEPA_plots/S06_fullGenesCopper.png", h = 2000, w = 2500, res = 300)
DotPlot(object = seuset_full, features = copper_NBL, split.by = "condition",
        scale=FALSE, col.min = -4, col.max = 4, 
        dot.min = 0, dot.scale = 5, cols = c("blue","red")) + RotatedAxis() + coord_flip() +
  theme(axis.text.x = element_text(size=7), axis.text.y = element_text(size=7))
dev.off()

neutroCells = subset(seuset_full, scType == "Neutrophils")
seuset_full_neutro <- merge(neutroCells, y = seuset_tumor, 
                            add.cell.ids = c("neutrophils", "tumor"), 
                            project = "singleCell")
seuset_full_neutro@meta.data$dataset[which(str_detect(seuset_full_neutro@meta.data$class, "immune"))] <- "neutrophils"
seuset_full_neutro@meta.data$dataset[which(str_detect(seuset_full_neutro@meta.data$class, "tumor"))] <- "tumor"

png("TEPA_plots/S06_Jun.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full_neutro) <- "dataset"
VlnPlot(seuset_full_neutro, features = c("Jun"), split.by = "condition",
        ncol = 1, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S06_Mt1_scType.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "scType"
VlnPlot(seuset_full, features = c("Mt1"), split.by = "condition",
        ncol = 1, pt.size = 0.000005)
dev.off()

png("TEPA_plots/S06_Sell_scType.png", h = 2000, w = 3500, res = 200)
Idents(seuset_full) <- "scType"
VlnPlot(seuset_full, features = c("Sell"), split.by = "condition",
        ncol = 1, pt.size = 0.000005)
dev.off()


png("TEPA_plots/S06_fullCopper.png", h = 4000, w = 6000, res = 300)
Idents(seuset_full) <- "condition"
#mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYBu"))(128)
DoHeatmap(object = seuset_full, size = 6,  group.by = "class", disp.min = -4, disp.max = 4,
          assay = "RNA", features = copper_genes, group.bar = T) +
  scale_fill_gradientn(colors = c("grey", "red"), breaks =  c(-4,-3, -2,-1, 0, 1, 2, 3, 4)) + 
  #scale_fill_continuous(breaks = c(-1, 0,1,2,3)) +
  theme(axis.text = element_text(size=15)) + 
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

png("TEPA_plots/S06_fullCopperImmune.png", h = 4000, w = 10000, res = 300)
Idents(seuset_immune) <- "condition"
DoHeatmap(object = seuset_immune, size = 6,  group.by = "scType",disp.min = -4, disp.max = 4,
          assay = "RNA", features = copper_genes, group.bar = T) +
  scale_fill_gradientn(colors = c("grey", "red"), breaks =  c(-4,-3, -2,-1, 0, 1, 2, 3, 4)) + 
  theme(axis.text = element_text(size=15)) + 
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()


