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
DefaultAssay(seuset_immune) <- "RNA"


grep(pattern = "Ly6", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)


#### Create module score ####
Macrophages <- c("Ccr2", "Ms4a6c", "Msr1", "Cd14")
Cd8_Tcells <- c("Top2a", "Pclaf", "Mki67")
Cd4_Tcells <- c("Cd3e", "Trac", "Trbc2", "Tcf7", "Cd4")
Ly6C_monocytes <- c("Fcgr4", "Pecam1")
Plasma_cells <- c("Igha", "Jchain", "Irf4", "Prdm1", "Sdc1", "Xbp1",
                  "Cd22", "Cd44", "Cxcr3", "Klf4", "Ly6a", "Ly6k") 
# https://www.biocompare.com/Editorial-Articles/585833-A-Guide-to-Plasma-Cell-Markers/
DCs <- c("Tmem176a", "Tmem176b", "Ciita", "Flt3")
Bcells <- c("Ighm", "Cd79a", "Ms4a1", "Cd19")
Th17_Tcells <- c("Tcrg.C1", "Rora", "Cxcr6", "Il2rb")
Neutrophils <- c("S100a8", "Cxcr2", "Csf3r")
NKs <- c("Ccl5", "Gzma", "Ncr1")
Basophils <- c("Il6", "Gata2", "Cyp11a1")
Eosinophils <- c("Alox15", "Itgam") #Itgam = Cd11b


#### Add module score to Seurat Object ####
seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Macrophages), name="Macrophages")
names(seuset_immune@meta.data)[grep("Macrophages", names(seuset_immune@meta.data))] <- "Macrophages"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Cd8_Tcells), name="Cd8_Tcells")
names(seuset_immune@meta.data)[grep("Cd8_Tcells", names(seuset_immune@meta.data))] <- "Cd8_Tcells"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Cd4_Tcells), name="Cd4_Tcells")
names(seuset_immune@meta.data)[grep("Cd4_Tcells", names(seuset_immune@meta.data))] <- "Cd4_Tcells"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Ly6C_monocytes), name="Ly6C_monocytes")
names(seuset_immune@meta.data)[grep("Ly6C_monocytes", names(seuset_immune@meta.data))] <- "Ly6C_monocytes"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Plasma_cells), name="Plasma_cells")
names(seuset_immune@meta.data)[grep("Plasma_cells", names(seuset_immune@meta.data))] <- "Plasma_cells"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(DCs), name="DCs")
names(seuset_immune@meta.data)[grep("DCs", names(seuset_immune@meta.data))] <- "DCs"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Bcells), name="Bcells")
names(seuset_immune@meta.data)[grep("Bcells", names(seuset_immune@meta.data))] <- "Bcells"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Th17_Tcells), name="Th17_Tcells")
names(seuset_immune@meta.data)[grep("Th17_Tcells", names(seuset_immune@meta.data))] <- "Th17_Tcells"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Neutrophils), name="Neutrophils")
names(seuset_immune@meta.data)[grep("Neutrophils", names(seuset_immune@meta.data))] <- "Neutrophils"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(NKs), name="NKs")
names(seuset_immune@meta.data)[grep("NKs", names(seuset_immune@meta.data))] <- "NKs"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Basophils), name="Basophils")
names(seuset_immune@meta.data)[grep("Basophils", names(seuset_immune@meta.data))] <- "Basophils"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Eosinophils), name="Eosinophils")
names(seuset_immune@meta.data)[grep("Eosinophils", names(seuset_immune@meta.data))] <- "Eosinophils"

seuset_immune <- AddModuleScore(seuset_immune, assay = "RNA", features = list(Platelets), name="Platelets")
names(seuset_immune@meta.data)[grep("Platelets", names(seuset_immune@meta.data))] <- "Platelets"


#### Plot scores ####
DefaultAssay(seuset_immune) <- "RNA"

png("TEPA_plots/clustersAnnot.png", h = 3000, w = 4200, res = 300)
FeaturePlot(seuset_immune, ncol = 4, 
            features = c("Macrophages", "Cd8_Tcells", "Cd4_Tcells", "Ly6C_monocytes", "Plasma_cells", "DCs", 
                         "Bcells", "Th17_Tcells", "Neutrophils", "NKs", "Basophils", "Eosinophils"), label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()

save(seuset_immune, file = "TEPA_results/seuset_ann_module.rda")
