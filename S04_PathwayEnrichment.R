#### Pathway Enrichment Analysis ####
## author: Antonietta Salerno
## date: 19/12/2022

library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library("Seurat")
library("SeuratDisk")
library(readr)
library(stringr)
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
library(ggplot2)
library(ggcharts)
library("EnhancedVolcano")
library(clusterProfiler)
library(pheatmap)

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadH5Seurat("TEPA_results/S03_immuneDiff.h5Seurat")
clusters = unique(seuset_immune@meta.data$scType)
immune.markers <- read.csv("TEPA_results/S03_DEA_clusterMarkers.csv")


#### 1 - Select the gene set collections of interest #### 

sets1 <- read.gmt("TEPA_data/mh.all.v2022.1.Mm.symbols.gmt") # Mouse hallmark
#sets2 <- read.gmt("TEPA_data/m2.cp.v2022.1.Mm.symbols.gmt") # Mouse curated canonical pathways M2-CP: BioCarrta, Reactome, WikiPathways
#sets3 <- read.gmt("TEPA_data/m5.all.v2022.1.Mm.symbols.gmt") # Mouse ontology gene sets M5: GO and tumor phenotype oncology
#sets4 <- read.gmt("TEPA_data/m2.cp.reactome.v2022.1.Mm.symbols.gmt") # Mouse Reactome subset of Canonical pathways
#sets5 <- msigdbr(species = "Mus musculus", category = "C7")
sets5 <- read.gmt("TEPA_data/REACTOME_NEUTROPHIL_DEGRANULATION.v2022.1.Mm.gmt") # The only Reactome pathway we're interested in

sets1$term <- as.character(sets1$term)
#sets4$term <- as.character(sets4$term)
sets5$term <- as.character(sets5$term)

sets1 <- sets1 %>% split(x = .$gene, f = .$term)
sets5 <- sets5 %>% split(x = .$gene, f = .$term)
#sets5 <- sets5 %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets <- append(sets1, sets5)
#fgsea_sets <- append(fgsea_sets, sets5)

#### 2 - Add custom gene set signatures ####
# Search all isoforms of gene of interest
grep(pattern = "Il1r", 
     x = rownames(x = seuset_immune@assays$RNA@data), 
     value = TRUE, ignore.case = TRUE)


N1 <- c("S100a8", "S100a9", "Icam1", "Fas", "Tnf", "Isg15", "Isg20", 
        "Ccl3", "Ccl4",  "Cxcl2", "Cxcl3", "Cebpb", "Il1a" ,"Il1b", "Il1r2", "Il1rn", "Il6ra", "Il15",
        "Stat3", "Hif1a", "Ifitm1","Ifitm2",  "Ifitm3", "Ifitm6",  
        "Acod1", "Myd88", "Prkcd", "Mmp8", "Mmp9","Retnlg", "Arg2")

# not in our dataset = "Cxcl10", "Il12a",  "Cxcl1",  "Ccl17", "Cxcl16", "Ccl7", "Il18", "Il18bp",

N2 <- c("Tgfb1", "Tgfb1i1","Tgfb2", "Tgfb3", "Ccl2",  "Ccl17","Cxcl14", 
        "Cxcl15",  "Il1r1", "Il2", "Il17a", "Mpo", "Slc27a2", "Arg1", 
         "Mrc1", "Chil3", "Elane", "Ctsg", "Retnla")

# Il1r1 = Il1 receptor antagonist
# Il8 = Cxcl15
# https://www.sciencedirect.com/science/article/pii/S1535610809002153?via%3Dihub
# https://www.sciencedirect.com/science/article/pii/S1471490605002425
# https://ashpublications.org/blood/article/133/20/2159/273823/Neutrophil-plasticity-in-the-tumor
# https://www.researchgate.net/figure/Transcriptomic-analysis-reveals-distinct-expression-signatures-during-neutrophil_fig3_364673856
# https://www.frontiersin.org/articles/10.3389/fimmu.2021.708770/full
# https://www.tandfonline.com/doi/pdf/10.1080/2162402X.2016.1232221
# Il1rn : https://www.spandidos-publications.com/10.3892/ijo.2015.2917 -> Il1a N2 (still expressed tho) and Il1b + Il1rn N1 ??

M1 <- c("Il12a", "Il12b", "Tnf", "Cxcl10", "Tlr2", "Tlr4", "Fcgr3", "Fcgr4", "Fcgr1",
         "Cd80", "Cd86", "Il12a",  "Il6", "Il1a", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
        "Cxcl9", "Cxcl10", "Ccr7", "Il1b", "Nos2", "Cxcl16", "Il23a")
M2 <- c("Il10", "Ccl22", "Il4ra", "Il13ra1", "Chil3", "Cd163", "Tgfb1",
        "Fcer2a", "Cd163", "Cxcr1", "Cxcr2", "Ccr2", "Arg1",  "Cd209a", "Mrc1", "Fcgr2b")

# https://www.sciencedirect.com/science/article/pii/S1471490602023025
# https://link.springer.com/protocol/10.1007/978-1-0716-0802-9_10/tables/1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3448669/

seuset_immune@assays$RNA@scale.data <- scale(seuset_immune@assays$RNA@data, scale = TRUE)
png("TEPA_plots/S04_neutrophilsPolarHeatmap.png", h = 4000, w = 6000, res = 300)
neutroCells = subset(seuset_immune, scType == "Neutrophils")
neutroCellsT = subset(seuset_immune, condition == "Treatment")
DoHeatmap(object = subset(neutroCellsT, downsample = 1000), size = 6, 
          assay = "RNA", features = c(N1,N2), group.bar = F) +
  scale_fill_gradientn(colors = c("blue", "black", "red")) + 
  theme(axis.text = element_text(size=15)) + NoLegend() +
  theme(plot.margin = margin(2,2,1.5,1.2, "cm"))
dev.off()

# Check if custom pathways match publicly available gene sets
pathways <- groupGO(M2, OrgDb = "org.Mm.eg.db", ont = "BP", keyType = "SYMBOL", level = 1, readable = T)
pathways@result$Description # just biological process...

# Create gene sets list
custom <- list(N1, N2, M1, M2)
names(custom) <- c("N1 ANTI-TUMOR NEUTROPHILS", "N2 PRO-TUMOR NEUTROPHILS",
                   "M1 ANTI-TUMOR MACROPHAGES", "M2 PRO-TUMOR MACROPHAGES")

fgsea_sets <- append(fgsea_sets, custom)

#### 3 - Run the custom gsea function ####

save = "S04_immuneEnrichment_"
gseaRES(clusters, fgsea_sets = fgsea_sets, save = save)

#gseaPlotRes(clusters)

#### 4 - Plot all clusters' results in a network ####

save = "S04_immuneJointNet"
gseaByCellType <- gseaJointNet(clusters, save = save)

# If you want to filter out some gene set manually!
df_gseaCT <- read.csv(paste0("TEPA_results/", save, "SHORT2.csv"), sep = ";")
gseaByCellType@compareClusterResult <- df_gseaCT

cnet <- cnetplot(gseaByCellType, showCategory=6, cex_label_category = 2, cex_label_gene = 1) 
file=paste0("TEPA_plots/",save,"2.png")
ggsave(cnet, file = file, width = 50, height = 40, units = "cm")


#### 5 - Clustered diverging bar plot all cell types

save = "S04_immuneJointBarplot"
gseaByType(clusters, save = save)

fgseaResByType = read.csv(paste0("TEPA_results/", save, "SHORT.csv"), sep = ";")
b <- barPlotGSEA(fgseaResByType, byType = TRUE)
ggsave(b, file=paste0("TEPA_plots/S04_barplotCellTypesEnrichedSHORT.png"),
       width = 40, height = 20, units = "cm", limitsize = F, dpi = 500)

### Only custom ###

# gsea
save = "S04_immuneEnrichmentCustom_"
gseaRES(clusters, fgsea_sets = custom, save = save, minSize = 7)

# barplot
save = "S04_immuneJointBarplotCustom"
gseaByType(clusters, save = save)

fgseaResByType = read.csv(paste0("TEPA_results/", save, ".csv"), sep = ",")
b <- barPlotGSEA(fgseaResByType, byType = TRUE)
ggsave(b, file=paste0("TEPA_plots/", save, ".png"),
       width = 30, height = 10, 
       units = "cm", limitsize = F, dpi = 500)

# gsea bulk
save = "S04_immuneEnrichmentBulkCustom_"
gseaRES("", fgsea_sets = custom, save = save, minSize = 7) # make it better



