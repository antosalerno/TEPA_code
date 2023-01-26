# Define a set of supporting function that will be used in the pipeline

# 5 - GSEA ####

# @ version of clusterProfiler:: cnetplot
networkPlotGSEA <- function(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 6){
  out <- tryCatch({
    cnet <- clusterProfiler::cnetplot(fgseaResPlot, showCategory = CNET_gene_cutoff, 
                                      categorySize="p.adjust", foldChange = ranked.genes) + 
      labs(title = paste("Significant Pathways in ", cluster ,sep = ""),
           subtitle = "by GSEA with GO database",
           caption = paste0("Based on ",
                            length(ranked.genes[ranked.genes > 0])," upregulated genes and ",
                            length(ranked.genes[ranked.genes < 0])," downregulated genes.",
                            " Showing ",CNET_gene_cutoff," genes per pathway."))+
      theme(plot.title = element_text(hjust = 0.5, size=20),
            plot.subtitle = element_text(hjust = 0.5, size=15),
            plot.caption = element_text(size=12))},
    finally = {
      return(cnet)
    })
  return(out)
}

barPlotGSEA <- function(fgseaRes, cluster){
  
  diverging_bar_chart(fgseaRes, pathway, as.numeric(NES), 
                      bar_color = c("#F9AF93", "#9BD6F0"),
                      text_color = "black") +
    ylab("Normalized enrichment score") +
    geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
    ggtitle(paste0("Enriched pathways in ", cluster," enriched by GSEA using GO database")) + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.title=element_text(size=20)) +
    theme_minimal() +
    theme(axis.text.y = element_blank(),  # Remove Y-axis texts
          axis.ticks.y = element_blank(), # Remove Y-axis ticks
          panel.grid.major.y = element_blank()) # Remove horizontal grid
  
}

# Run fgsea function from DEA results using GO database and produce networkPlot and barPlot 
gseaGO <- function(clusters, fgsea_sets){
  for (cluster in clusters){
    message(paste0("GSEA for cluster: ", cluster))
    res <- read.csv(paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"), sep=",")
    ranked.genes<- deframe(res %>%
                             dplyr::filter(p_val_adj < 0.1) %>%
                             arrange(desc(avg_log2FC)) %>% 
                             dplyr::select(X, avg_log2FC))
    ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
    
    # Run GSEA
    fgseaRes<- fgsea(fgsea_sets, stats = ranked.genes, minSize = 5)
    fgseaRes <- fgseaRes[fgseaRes$padj <=0.1] %>% arrange(desc(NES))
    fgseaRes <- as.data.frame(fgseaRes)
    
    
    # Generate enrichment plots 
    if (nrow(fgseaRes) > 0) {
      fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                         f = .$pathway)
      for (e in 1: length(fgseaResPlot)){
        leEd <- unlist(str_extract_all(fgseaRes$leadingEdge[e],"[[:alnum:]]{2,}"))
        fgseaResPlot[e][[1]] <- leEd
      }
      
      ## Network
      cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
      ggsave(cnet, file=paste0("TEPA_plots/05_clProf_cnetGO_",cluster,".png"), width = 30, height = 30, units = "cm")
      
      ## Barplot
      b <- barPlotGSEA(fgseaRes, cluster)
      ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotGO_",cluster,".png"),width = 20, height = nrow(fgseaRes)/1.5, units = "cm")
      
      # Save dataframe
      df <- apply(fgseaRes,2,as.character)
      write.csv(df, 
                file = paste0("TEPA_results/05_GSEAclusterGO",cluster,".csv"))
    }
  }
}

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
        b <- barPlotGSEA(fgseaRes, cluster)
        ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotReactome_",cluster,".png"),width = 20, height = nrow(fgseaRes)/1.5, units = "cm")
        
        # Save dataframe
        df <- apply(fgseaRes,2,as.character)
        write.csv(df, 
                  file = paste0("TEPA_results/05_GSEAclusterReactome",cluster,".csv"))
      }
    }
  }
}



