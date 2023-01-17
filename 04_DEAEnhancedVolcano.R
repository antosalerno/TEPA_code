#### DEA with Enhanced Volcano ####
## author: Antonietta Salerno
## date: 16/12/2022

library("EnhancedVolcano")
setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
seuset_immune <- LoadH5Seurat("TEPA_results/03_seusetImmuneModule.h5Seurat")

Idents(seuset_immune) <- "scType"
clusters = unique(Idents(seuset_immune))

# Read the DEA tables for each cluster and then plot

for (cluster in clusters){
  res <- read.csv(paste0("TEPA_results/02_DEAcluster",cluster,".csv"), sep=",")
  rownames(res) <- res$X
  p <- EnhancedVolcano(res, subtitle = "",
                       #selectLab = markers,
                       lab = rownames(res),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       xlim = c(-2.5, 2.5),
                       title = paste0(cluster,', TEPA vs. CTRL '),
                       pCutoff = 0.05, #0.05 cutoff
                       FCcutoff = 0.5, # 2-fold change
                       labFace = "bold",
                       labSize = 3,
                       col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                       colAlpha = 4/5,
                       legendLabSize = 14,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 50,
                       caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC>0.75&res$p_val_adj<=0.05,]), ' genes',
                                        "\n",'Downregulated = ', nrow(res[res$avg_log2FC< -0.75&res$p_val_adj<=0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
  ggsave(p, file=paste0("TEPA_plots/04_DEAcluster_",cluster,".png"), width = 14, height = 20, units = "cm")
}
