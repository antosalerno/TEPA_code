#### Pathway Enrichment Analysis ####
## author: Antonietta Salerno
## date: 19/12/2022


library("fgsea")
library(dplyr)
library(ggplot2)
library(devtools)
library(msigdbr)
library("tibble")

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
seuset_immune <- LoadH5Seurat("TEPA_results/03_seusetImmuneModule.h5Seurat")
clusters = unique(seuset_immune@meta.data$scType)

#### fgsea ####

# msigdbr_species()
m_df<- msigdbr(species = "Mus musculus", category = "C7")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


for (cluster in clusters){
  res <- read.csv(paste0("TEPA_results/02_DEAcluster",cluster,".csv"), sep=",")
  ranked.genes<- deframe(res %>%
                           dplyr::filter(cluster == cluster, p_val_adj < 0.05) %>%
                           arrange(desc(avg_log2FC)) %>% 
                           dplyr::select(X, avg_log2FC))
  
  fgseaRes<- fgsea(fgsea_sets, stats = ranked.genes)
  fgseaRes <- fgseaRes[fgseaRes$padj <=0.05] %>% 
    arrange(desc(NES))
  fgseaRes <- apply(fgseaRes,2,as.character)
  
  write.csv(fgseaRes, 
            file = paste0("TEPA_results/05_GSEAcluster",cluster,".csv"))
  
}

#### clusterProfiler ####

library(clusterProfiler)
library("pathview")
library(enrichplot)
library(ggplot2)
require(DOSE)

# Install organism annotation
organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)
library(GOSemSim)
d <- godata(organism, ont="BP")

# keytypes(org.Mm.eg.db)

# Produce plots
plot <- function(cluster, gse){
  r <- ridgeplot(gse) + labs(x = "enrichment distribution")
  ggsave(r, file=paste0("TEPA_plots/05_clProf_ridge_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  gse_emap <- pairwise_termsim(gse, method="Wang", semData = d)
  emap <- emapplot(gse_emap, showCategory = 10)
  ggsave(emap, file=paste0("TEPA_plots/05_clProf_emap_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  emap_cl <- emapplot_cluster(gse_emap)
  ggsave(emap_cl, file=paste0("TEPA_plots/05_clProf_emapCl_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  d <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  ggsave(d, file=paste0("TEPA_plots/05_clProf_dotplot_",cluster,".png"), width = 20, height = 20, units = "cm")
  
  cnet <- cnetplot(gse, categorySize="p.adjust", foldChange = ranked.genes)
  ggsave(cnet, file=paste0("TEPA_plots/05_clProf_cnet_",cluster,".png"), width = 20, height = 20, units = "cm")
  
}

# Run enrichment function

for (cluster in clusters){
  res <- read.csv(paste0("TEPA_results/02_DEAcluster",cluster,".csv"), sep=",")
  
  ranked.genes<- deframe(res %>%
                           dplyr::filter(cluster == cluster, p_val_adj <= 0.05) %>%
                           arrange(desc(avg_log2FC)) %>% 
                           dplyr::select(X, avg_log2FC))
  ranked.genes<- sort(ranked.genes, decreasing = TRUE)
  gse <- gseGO(geneList=ranked.genes, 
               ont ="BP", 
               keyType = "SYMBOL", 
               eps = 0,
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "BH")
  
  write.csv(gse, file = paste0("TEPA_results/05_clProf_",cluster,".csv"))
  if (nrow(gse@result) > 1){plot(cluster, gse)} # if there's an enrichment you run the plotting function

}






