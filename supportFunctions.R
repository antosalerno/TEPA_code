# Define a set of supporting function that will be used in the pipeline #
## author: Antonietta Salerno
## date: 20/01/2023

# S00 - Density histogram ####

hist_dens <- function(x, breaks = "Scott", 
                      main = "Mycn expression per number of cells",
                      xlab = "Mycn", ylab = "Cell counts") {
  dens <- density(x, na.rm = T)
  raw_hist <- hist(x, breaks = breaks, plot = F)
  scale <- max(raw_hist$counts)/max(raw_hist$density)
  hist(x, breaks = breaks, prob = F, main = main, xlab = xlab, ylab = ylab)
  lines(list(x = dens$x, y = scale * dens$y), col = "blue", lwd = 2)
}

# S03 - Differential expression analysis ####

# Create signature of cell type with top markers
createSets <- function(markers = immune.markers,
                       obj = seuset_immune,
                       id = "scType"){
  Idents(obj) <- id
  clusters = unique(Idents(obj))
  for(c in 1:length(clusters)){
    clusterDF <- markers[markers$cluster == clusters[c],]
    if (nrow(clusterDF) > 0){
      clusterDF <- clusterDF[clusterDF$p_val_adj < 0.05,]
      clusterDF <- clusterDF[order(-c(clusterDF$avg_log2FC)),]
      top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene)
      obj <- AddModuleScore(obj, assay = "RNA", features = list(top_genes), name=make.names(as.character(clusters[c])))
      names(obj@meta.data)[grep(make.names(as.character(clusters[c])), names(obj@meta.data))] <- as.character(clusters[c])
    }
  }
  return(obj)
}

# Get top marker genes by cluster
getTopMarkers <- function(df = immune.markers, geneNr){
  markers <- c()
  for(c in 1:length(clusters)){
    clusterDF = df[df$cluster == clusters[c],]
    clusterDF <- clusterDF[order(-clusterDF$avg_log2FC),]
    top_genes <- head(clusterDF[clusterDF$cluster == clusters[c],]$gene, geneNr)
    markers <- append(markers, top_genes)
  }
  return(markers)
}

# Volcano plot

library("EnhancedVolcano")

plotVolcano <- function(clusters, res=NULL, type = "condition",
                        save = "", log2FC = 0.25, immune = TRUE, pval = 0.05){
  for (cluster in clusters){
    try({
      if (type == "condition"){
        if(immune){
          file=paste0("TEPA_results/S03_immuneCond_DEA_",sub(" ", "_", cluster),".csv")
        } else {
          file=paste0("TEPA_results/S05_tumorCond_DEA_",sub(" ", "_", cluster),".csv")
        } 
        res <- read.csv(file, sep=",")
        rownames(res) <- res$X
        title = paste0(cluster,', TEPA vs. CTRL ')
      } else if (type == "markers") {
        if (immune){res <- immune.markers[immune.markers$cluster == cluster,]}
        else {res <- tumor.markers[tumor.markers$cluster == cluster,]}
        title = paste0('DEA cluster ', cluster, " vs all")
      } else {message("Please specify the type of analysis: either 'markers' or 'condition'")}
      p <- EnhancedVolcano(res, subtitle = "",
                         #selectLab = markers,
                         lab = rownames(res),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         xlim = c(-2.5, 2.5),
                         title = title,
                         pCutoff = pval, 
                         FCcutoff = log2FC, # 2-fold change
                         labFace = "bold",
                         labSize = 3,
                         col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                         colAlpha = 4/5,
                         legendLabSize = 14,
                         legendIconSize = 4.0,
                         drawConnectors = TRUE,
                         widthConnectors = 0.3,colConnectors = 'gray51', maxoverlapsConnectors = 80,
                         caption = paste0('Upregulated = ', nrow(res[res$avg_log2FC > log2FC&res$p_val_adj<= 0.05,]), ' genes',
                                          "\n",'Downregulated = ', nrow(res[res$avg_log2FC < -log2FC &res$p_val_adj<= 0.05,]), ' genes'))+ theme(plot.title = element_text(hjust = 0.5)) + coord_flip()
      ggsave(p, file=paste0("TEPA_plots/",save, "DEA_",sub(" ", "_", cluster),".png"), width = 30, height = 25, units = "cm")
    })
  }
}

# S05 - GSEA ####

N1 <- c("S100a8", "S100a9", "Icam1", "Fas", "Tnf", "Isg15", "Isg20", 
        "Ccl3", "Ccl4",  "Cxcl2", "Cxcl3", "Cebpb", "Il1a" ,"Il1b", "Il1r2", "Il1rn", "Il6ra", "Il15",
        "Stat3", "Hif1a", "Ifitm1","Ifitm2",  "Ifitm3", "Ifitm6",  
        "Acod1", "Myd88", "Prkcd", "Mmp8", "Mmp9","Retnlg", "Arg2")
N2 <- c("Tgfb1", "Tgfb1i1","Tgfb2", "Tgfb3", "Ccl2",  "Ccl17","Cxcl14", 
        "Cxcl15",  "Il1r1", "Il2", "Il17a", "Mpo", "Slc27a2", "Arg1", 
        "Mrc1", "Chil3", "Elane", "Ctsg", "Retnla")
M1 <- c("Il12a", "Il12b", "Tnf", "Cxcl10", "Tlr2", "Tlr4", "Fcgr3", "Fcgr4", "Fcgr1",
        "Cd80", "Cd86", "Il12a",  "Il6", "Il1a", "Ccl2", "Ccl3", "Ccl4", "Ccl5",
        "Cxcl9", "Cxcl10", "Ccr7", "Il1b", "Nos2", "Cxcl16", "Il23a")
M2 <- c("Il10", "Ccl22", "Il4ra", "Il13ra1", "Chil3", "Cd163", "Tgfb1",
        "Fcer2a", "Cd163", "Cxcr1", "Cxcr2", "Ccr2", "Arg1",  "Cd209a", "Mrc1", "Fcgr2b")

custom <- list(N1, N2, M1, M2)
names(custom) <- c("N1 ANTI-TUMOR NEUTROPHILS", "N2 PRO-TUMOR NEUTROPHILS",
                   "M1 ANTI-TUMOR MACROPHAGES", "M2 PRO-TUMOR MACROPHAGES")


# @ version of clusterProfiler:: cnetplot
networkPlotGSEA <- function(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 6){
  out <- tryCatch({
    cnet <- clusterProfiler::cnetplot(fgseaResPlot, showCategory = CNET_gene_cutoff, 
                                      categorySize="p.adjust", color.params = list(foldChange = ranked.genes)) + 
      labs(title = paste("Significant Pathways in ", cluster ,sep = ""),
           subtitle = "enriched by GSEA",
           caption = paste0("Based on ",
                            length(ranked.genes[ranked.genes > 0])," upregulated genes and ",
                            length(ranked.genes[ranked.genes < 0])," downregulated genes.",
                            " Showing ",CNET_gene_cutoff," genes per pathway."))+
      theme(plot.title = element_text(hjust = 0.5, size=20),
            plot.subtitle = element_text(hjust = 0.5, size=15),
            plot.caption = element_text(size=12))
    },
    finally = {
      return(cnet)
    })
  return(out)
}

library(scales)

barPlotGSEA <- function(fgseaRes, cluster = NULL, name = "", byType = FALSE){
  
  # When you want to plot the enrichment of gene sets by cell type
  if (byType){
    
    clusters_colors = hue_pal()(length(clusters))
    order_cl <- sort(unlist(clusters))
    c <- setNames(clusters_colors, order_cl)
    
    fgseaRes %>% arrange(pathway, NES) %>%
      ggplot(aes(x = pathway, 
                 y = as.numeric(NES),
                 fill = reorder(cluster, as.numeric(NES)))) +
      scale_colour_manual(values = c) +
      geom_bar(position=position_dodge2(width = 1.5, preserve = "single", padding = 0),
               stat='identity', width = 1) +
      ylab("Normalized enrichment score") +
      geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
      scale_x_discrete(expand = c(.05, .05)) +
      coord_cartesian(xlim = c(-5, 5)) +
      coord_flip() +
      # scale_y_discrete(limits = factor(min(fgseaRes$NES), max(fgseaRes$NES)), breaks = seq(min(fgseaRes$NES), max(fgseaRes$NES),by = 1)) +
      labs(title = "Top pathways enriched by cell type",
           x = NULL,
           fill = NULL) +
      theme_minimal() +
      theme(panel.grid.major.y = element_blank(),
            legend.position = "right") +
      theme(axis.text.y = element_text(size = 10, hjust = 1),
            plot.margin = margin(rep(10, 8)))
  }
  
  # When you want to plot the enrichment of a specific cluster
  #if (is.null(cluster) == TRUE)
  else {
  diverging_bar_chart(fgseaRes, pathway, as.numeric(NES), 
                      bar_color = c("#F9AF93", "#9BD6F0"),
                      text_color = "black", text_size = 10) +
    ylab("Normalized enrichment score") +
    geom_hline(yintercept = 0, color = 1, lwd = 0.2) +
    ggtitle(paste0("Significant pathways in ", name," enriched by GSEA")) + 
    theme(plot.title = element_text(hjust = 1)) + 
    theme(axis.title=element_text(size=20)) +
    theme_minimal() +
    ylim(-5, 5) +
    theme(axis.text.y = element_blank(),  # Remove Y-axis texts
          axis.ticks.y = element_blank(), # Remove Y-axis ticks
          panel.grid.major.y = element_blank()) # Remove horizontal grid
  }
}


# Run fgsea function from DEA results using GO database and produce networkPlot and barPlot 
gseaRES <- function(clusters, markers = NULL, fgsea_sets, minSize = 12, adj = TRUE, 
                    type = "condition", input = "immune", save=""){
  for (cluster in clusters){
    try({
      if (type == "condition"){
        if(input == "immune"){
          
          file=paste0("TEPA_results/S03_immuneBulkDEA.csv")
          #file=paste0("TEPA_results/S03_immuneCond_DEA_",cluster,".csv")
        } else if (input == "tumor"){
            file=paste0("TEPA_results/S05_tumorCond_DEA_",cluster,".csv")
        } else if (input == "nanoCond_gInf"){
            file = paste0("TEPA_results/N00_nanoCond_gInf_DEA.csv")
        } else if (input == "nanoInf_gCond"){
            file = paste0("TEPA_results/N00_nanoInf_gCond_DEA.csv")
        } else if (input == "nanoCond"){
           file = paste0("TEPA_results/N00_nanoCond_DEA.csv")
        }
        res <- read.csv(file = file, sep=",")
      } else if (type == "markers") {
        res <- markers[markers$cluster == cluster,]
        res$X <- res$gene
      } else {message("Please specify the type of analysis: either 'markers' or 'condition'")}

      message(paste0("GSEA for cluster: ", cluster))
      ranked.genes<- deframe(res %>%
                             #dplyr::filter(p_val_adj < 0.1) %>%
                             arrange(desc(avg_log2FC)) %>% 
                             dplyr::select(X, avg_log2FC))
      ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
    
      # Run GSEA
      fgseaRes<- fgsea(fgsea_sets, stats = ranked.genes, minSize = minSize, maxSize = Inf)
      if (adj) {fgseaRes$padj = p.adjust(fgseaRes$pval, method='BH')}
      fgseaRes <- fgseaRes[fgseaRes$padj <= 0.2] %>% arrange(desc(NES))
      fgseaRes <- as.data.frame(fgseaRes)
    
      # Generate enrichment plots 
      if (nrow(fgseaRes) > 0) {
        fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                         f = .$pathway)
        for (e in 1: length(fgseaResPlot)){
          fgseaResPlot[[e]] <- unlist(fgseaResPlot[[e]])
        }
      
      ## Network
      cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
      file=paste0("TEPA_plots/",save, "cnet_",sub(" ", "_", cluster),".png")
      ggsave(cnet, file = file, width = 30, height = 30, units = "cm")
      
      ## Barplot
      b <- barPlotGSEA(fgseaRes, cluster, name = cluster)
      file=paste0("TEPA_plots/",save, "barplot_",sub(" ", "_", cluster),".png")
      ggsave(b, file = file,width = 20, height = ifelse(nrow(fgseaRes) > 4, nrow(fgseaRes)/1.5, nrow(fgseaRes) + 3), units = "cm")
      
      # Save dataframe
      df <- apply(fgseaRes,2,as.character)
      file=paste0("TEPA_results/",save, "GSEA_",sub(" ", "_", cluster),".csv")
      write.csv(df, file = file)
      }
    }
  )}
}

gseaPlotRes <- function(clusters, save=""){
  for (cluster in clusters){
    message(paste0("GSEA for cluster: ", cluster))
    
    # Load DEA analysis
    res <- read.csv(paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"), sep=",")
    ranked.genes<- deframe(res %>%
                             dplyr::filter(p_val_adj < 0.01) %>%
                             arrange(desc(avg_log2FC)) %>% 
                             dplyr::select(X, avg_log2FC))
    ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
    
    # Load GSEA analysis
    filename = paste0("TEPA_results/05_GSEAclusterMSIG",cluster,"SHORT.csv")
    
    if (file.exists(filename)){
      fgseaRes <- read.csv(filename, sep=";")
    
      fgseaResPlot <- fgseaRes %>% split(x = .$leadingEdge, 
                                         f = .$pathway)
      for (e in 1: length(fgseaResPlot)){
        leEd <- unlist(str_extract_all(fgseaRes$leadingEdge[e],"[[:alnum:]]{2,}"))
        fgseaResPlot[e][[1]] <- leEd
      }
      
      ## Network
      cnet <- suppressWarnings(networkPlotGSEA(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 10))
      ggsave(cnet, file=paste0("TEPA_plots/05_clProf_cnetMSIG_",paste(cluster, collapse = "_"),save,".png"), width = 30, height = 30, units = "cm")
      
      ## Barplot
      b <- barPlotGSEA(fgseaRes, cluster)
      ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotMSIG_",paste(cluster, collapse = "_"),save,".png"),width = 50, height = ifelse(nrow(fgseaRes) > 3, nrow(fgseaRes)/1.5, nrow(fgseaRes) + 1), units = "cm")
    }
  }
}
  

gseaByType <- function(clusters, save){
  
  fgseaResByType <- data.frame(
    pathway = character(0),
    NES = numeric(0),
    cluster = character(0)
  )
  for (cluster in clusters){
    index = which(clusters == cluster)
    filename = paste0("TEPA_results/S04_immuneEnrichment_GSEA_",sub(" ", "_", cluster),".csv")
    if (file.exists(filename)){
      clusterFile <- read.csv(filename, sep = ",")
      if (ncol(clusterFile) == 2){
        fgseaResByType[nrow(fgseaResByType)+1,"pathway"] <- clusterFile[1,2]
        fgseaResByType[nrow(fgseaResByType),"NES"] <- clusterFile[6,2]
        fgseaResByType[nrow(fgseaResByType),"cluster"] <- cluster
      }
      else {
        for (i in 1: nrow(clusterFile)){
          fgseaResByType[nrow(fgseaResByType)+1,"pathway"] <- clusterFile[i,"pathway"]
          fgseaResByType[nrow(fgseaResByType),"NES"] <- clusterFile[i,"NES"]
          fgseaResByType[nrow(fgseaResByType),"cluster"] <- cluster
        }
      }
    }
  }
  
  # Save dataframe
  df <- apply(fgseaResByType,2,as.character)
  write.csv(df, file = paste0("TEPA_results/", save, ".csv"))
  
}

# Lots of time -> try Reactome as well
gseaJointNet <- function (clusters, immune = TRUE, save=""){
  cell_types = list()
  for (cluster in clusters){
    if(immune){
      file=paste0("TEPA_results/S03_immuneCond_DEA_",cluster,".csv")
    } else {
      file=paste0("TEPA_results/S07_tumorCond_DEA_",cluster,".csv")
    } 
    res <- read.csv(file = file, sep=",")
    ranked.genes<- deframe(res %>%
                             dplyr::filter(p_val_adj < 0.05) %>%
                             dplyr::filter(avg_log2FC > 0.2) %>%
                             arrange(desc(avg_log2FC)) %>% 
                             dplyr::select(X, avg_log2FC))
    # We take only significant positively enriched genes since we are running ORA
    if (length(ranked.genes) > 0){
      cell_types[cluster] <- list(names(ranked.genes))
    
    }
  }
  
  gseaByCellType <- compareCluster(cell_types, fun="enrichGO", ont = "BP", readable = TRUE, keyType = "SYMBOL",
                                   OrgDb = "org.Mm.eg.db", pvalueCutoff=0.05, minGSSize = 15, maxGSSize = 500)
  
  df_gseaCT <- gseaByCellType@compareClusterResult
  df_gseaCT$geneNum <- as.numeric(str_extract_all(df_gseaCT$GeneRatio,"[[:digit:]]{1,}", simplify = T)[,1])
  df_gseaCT <- df_gseaCT %>%
    filter(p.adjust <=0.05) %>%
    filter(geneNum > 10) # remove gene sets matches having less than 10 genes 
  
  df_gseaCT <- as.data.frame(df_gseaCT)
  
  # Save dataframe
  df <- apply(df_gseaCT,2,as.character)
  file=paste0("TEPA_results/",save, ".csv")
  write.csv(df, file = file)
  
  return(gseaByCellType)
}

