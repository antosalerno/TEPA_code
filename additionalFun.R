# Additional function not used anymore #
# date: 20/02/2023
# author: Antonietta Salerno

# Run fgsea function from DEA results using Reactome database and produce networkPlot and barPlot 
gseaReact <- function(clusters){
  
  for (cluster in clusters){
    res <- read.csv(paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"), sep=",")
    ranked.genes<- deframe(res %>%
                             dplyr::filter(p_val_adj < 0.1) %>%
                             arrange(desc(avg_log2FC)) %>% 
                             dplyr::select(X, avg_log2FC))
    if (length(ranked.genes) > 0){
      # Convert into Entrez ID
      entrezID <- AnnotationDbi::select(org.Mm.eg.db, keys=names(ranked.genes), columns='ENTREZID', keytype='SYMBOL')
      ranked.genesENTREZ <- ranked.genes
      names(ranked.genesENTREZ) <- entrezID$ENTREZID
      ranked.genesENTREZ <- ranked.genesENTREZ[!is.na(names(ranked.genesENTREZ))]
      
      # Get Reactome pathways
      pathways <- reactomePathways(names(ranked.genesENTREZ))
      # Run GSEA
      fgseaRes<- fgsea(pathways, stats = ranked.genesENTREZ, minSize = 5)
      fgseaRes <- fgseaRes[fgseaRes$padj <=0.05] %>% arrange(desc(NES))
      fgseaRes <- as.data.frame(apply(fgseaRes,2,as.character))
      
      if (nrow(fgseaRes) > 0) {
        # Plot enrichment network after mapping genes back to symbols
        fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                           f = .$pathway)
        for (e in 1: length(fgseaResPlot)){
          leEd <- unlist(str_extract_all(fgseaRes$leadingEdge,"[[:digit:]]{4,}"))
          symbolID <- AnnotationDbi::select(org.Mm.eg.db, keys=leEd, columns='SYMBOL', keytype='ENTREZID')
          fgseaResPlot[e][[1]] <- symbolID$SYMBOL
        }
        
        ## Network
        cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
        ggsave(cnet, file=paste0("TEPA_plots/05_clProf_cnetReactome_",cluster,".png"), width = 20, height = 20, units = "cm")
        
        ## Barplot
        b <- barPlotGSEA(fgseaRes, cluster, name = cluster)
        ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotReactome_",cluster,".png"),width = 20, height = nrow(fgseaRes)/1.5, units = "cm")
        
        # Save dataframe
        df <- apply(fgseaRes,2,as.character)
        write.csv(df, 
                  file = paste0("TEPA_results/05_GSEAclusterReactome",cluster,".csv"))
      }
    }
  }
}


#### Subclustering ####

Idents(seuset_tumor) <- "clust_res.0.4"
seuset_tumor <- FindSubCluster(seuset_tumor, "3", "clust", subcluster.name = "seurat_clusters",  resolution = 0.5, algorithm = 3)
DimPlot(seuset_tumor, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3)
Idents(seuset_tumor) <- "seurat_clusters"

# Rename clusters
index = which(seuset_tumor@meta.data$seurat_clusters == "1")
seuset_tumor@meta.data$seurat_clusters <- replace(seuset_tumor@meta.data$seurat_clusters, index, "C1")
