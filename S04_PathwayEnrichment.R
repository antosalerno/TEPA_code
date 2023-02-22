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

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadH5Seurat("TEPA_results/S03_immuneDiff.h5Seurat")
clusters = unique(seuset_immune@meta.data$scType)

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


# Check if custom pathways match publicily available gene sets
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



