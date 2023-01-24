#### Pathway Enrichment Analysis ####
## author: Antonietta Salerno
## date: 19/12/2022


library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library("tibble")
library("Seurat")
library("SeuratData")
library("SeuratDisk")
library(readr)
library("DESeq2")
library(stringr)
library(organism, character.only = TRUE)
library(GOSemSim)
library("GOfuncR")
library(clusterProfiler)
library("pathview")
library(enrichplot)
library(ggplot2)
library("AnnotationDbi")
require(DOSE)
source("TEPA_code/supportFunctions.R")


setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
seuset_immune <- LoadH5Seurat("TEPA_results/03_seusetImmuneModule.h5Seurat")
clusters = unique(seuset_immune@meta.data$scType)

### DEA with less stringent method so we get more genes for GSEA ####
Idents(seuset_immune) <- "celltype.tepa"
DefaultAssay(seuset_immune) <- "RNA"

#sheets <- list()
for (cluster in unique(seuset_immune$celltype)){
  try({
    ident1 <- paste0(cluster,"_Control")
    ident2 <- paste0(cluster,"_Treatment")
    condition.diffgenes <- FindMarkers(seuset_immune, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE,
                                       #latent.vars="orig.ident",
                                       test.use="DESeq2")
    #sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/02_DEAclusterDESEQ2",cluster,".csv"))
  })
}


### fgsea ####

# 1 - GO database (best)

dbGO <- godata("org.Mm.eg.db", keytype = "SYMBOL", ont="BP")
dbGOAnn <- dbGO@geneAnno
dbGOAnn <- dbGOAnn[,-c(3:4)]
dbGOAnn <- as.data.frame(dbGOAnn)
dbGOAnn <- cbind(dbGOAnn,get_names(dbGOAnn$GO))
dbGOAnn <- dbGOAnn %>% distinct()
fgsea_sets<- dbGOAnn %>% split(x = .$SYMBOL, f = .$go_name)

gseaGO(clusters, fgsea_sets)

# 2 - MsigDBr database
m_df<- msigdbr(species = "Mus musculus", category = "C7")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

gseaGO(clusters, fgsea_sets)

# 3 - Reactome database

gseaReact(clusters)


### Alternative way with clusterProfiler::gseGO ####

# Produce plots
plot <- function(cluster, gse){
  r <- ridgeplot(gse) + labs(x = "enrichment distribution")
  ggsave(r, file=paste0("TEPA_plots/05_clProf_ridge_",cluster,".png"), width = 40, height = 20, units = "cm")
  
  gse_emap <- pairwise_termsim(gse, method="Wang", semData = d)
  emap <- emapplot(gse_emap, showCategory = 10)
  ggsave(emap, file=paste0("TEPA_plots/05_clProf_emap_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  d <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ggsave(d, file=paste0("TEPA_plots/05_clProf_dotplot_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  cnet <- cnetplot(gse, categorySize="p.adjust", color.params = list(ranked.genes))
  ggsave(cnet, file=paste0("TEPA_plots/05_clProf_cnet_",cluster,".png"), width = 20, height = 20, units = "cm")
  
}

# Run enrichment function

for (cluster in clusters){
  res <- read.csv(paste0("TEPA_results/02_DEAcluster",cluster,".csv"), sep=",")
  
  ranked.genes<- deframe(res %>%
                           dplyr::filter(p_val_adj <= 0.2) %>%
                           arrange(desc(avg_log2FC)) %>% 
                           dplyr::select(X, avg_log2FC))
  ranked.genes <- ranked.genes[!is.na(ranked.genes)]
  ranked.genes<- sort(ranked.genes, decreasing = TRUE)
  if (length(ranked.genes) > 3) {
    gse <- clusterProfiler::gseGO(geneList=ranked.genes, 
                 OrgDb = "org.Mm.eg.db",
               ont ="BP", 
               keyType = "SYMBOL", 
               eps = 0,
               minGSSize = 3, 
               #maxGSSize = 800, 
               pvalueCutoff = 1, 
               verbose = TRUE, 
               pAdjustMethod = "BH"
               )
  
  write.csv(gse, file = paste0("TEPA_results/05_clProf_",cluster,".csv"))
  if (nrow(gse@result) > 1){plot(cluster, gse)} # if there's an enrichment you run the plotting function
  }
}




