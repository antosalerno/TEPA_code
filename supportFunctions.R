# Define a set of supporting function that will be used in the pipeline

# 5 - GSEA ####
# @ modified from genekitr package
plot_fgsea <- function (gsea_list, cluster, plot_type = c("volcano", "classic", "fgsea", "ridge", "bar"), 
                        stats_metric = c("p.adjust", "pvalue", "qvalue"), show_pathway = NULL, 
                        show_gene = NULL, colour = NULL, wrap_length = NULL, ...) {
  lst <- list(...)
  plot_type <- match.arg(plot_type)
  gsea_df <- gsea_list
  colnames(gsea_df) <- c("Description", "pvalue", "p.adjust", 
                         "log2err", "enrichmentScore", "NES", "setSize", "core_enrichment")
  ID <- gsea_df$Description
  gsea_df <- cbind(ID, gsea_df)
  gsea_df$p.adjust <- as.numeric(gsea_df$p.adjust)
  gsea_df$pvalue <- as.numeric(gsea_df$pvalue)
  gsea_df$NES <- as.numeric(gsea_df$NES)
  stats_metric_label <- ifelse(stats_metric == "pvalue", "pval", 
                               ifelse(stats_metric == "p.adjust", "padj", "FDR"))
  if (plot_type == "volcano") {
    if (!"main_text_size" %in% names(lst)) 
      lst$main_text_size <- 8
    if (is.null(show_pathway)) 
      show_pathway = 3
    gsea_df <- gsea_df[order(gsea_df$NES, decreasing = TRUE), 
    ]
    plot_df <- gsea_df
    if (is.numeric(show_pathway)) {
      plot_df$group <- c(rep("Up", show_pathway), rep("ignore", 
                                                      nrow(plot_df) - 2 * show_pathway), rep("Down", 
                                                                                             show_pathway))
    }
    else {
      nes <- plot_df[plot_df$ID %in% show_pathway, "NES"]
      plot_df$group <- "ignore"
      plot_df[plot_df$ID %in% show_pathway, "group"] <- sapply(nes, 
                                                               function(x) ifelse(x > 0, "Up", "Down"))
    }
    if (is.null(colour)) {
      message("\"colour\" is NULL, now using \"red\", \"grey\" and \"blue\"...")
      colour <- c("#E31A1C", "grey", "#1F78B4")
    }
    p <- ggplot(plot_df, aes(x = NES, y = -log10(eval(parse(text = stats_metric))), 
                             color = group)) + geom_point(alpha = 0.6, size = 3.5) + 
      xlab("NES") + ylab(paste0("-log10(", stats_metric_label, 
                                ")")) + ggrepel::geom_text_repel(data = plot_df[plot_df$group != 
                                                                                  "ignore", ], aes(label = ID), size = (lst$main_text_size/4), 
                                                                 color = "black", show.legend = F) + xlab("Normalized enrichment score") + 
      ylim(0, NA) + scale_color_manual(breaks = c("Up", 
                                                  "Down"), values = colour, name = "") + plot_theme(...)
  }
  if (plot_type == "classic") {
    gl1 = gl[, 2]
    names(gl1) = gl[, 1]
    genelist = gl1
    exponent <- as.numeric(gsea_list$exponent)
    org <- as.character(gsea_list$org)
    if (is.null(show_pathway)) 
      show_pathway = 3
    if (is.numeric(show_pathway)) {
      show_pathway <- gsea_df$ID[show_pathway]
    }
    else if (any(!show_pathway %in% gsea_df$ID)) {
      stop(paste0(show_pathway[!show_pathway %in% gsea_df$ID], 
                  " not in GSEA result!"))
    }
    if (is.null(colour)) {
      colour <- c("\\#5DA5DAFF", "\\#FAA43AFF", "\\#60BD68FF", 
                  "\\#F15854FF", "\\#B276B2FF", "\\#8D4B08FF", 
                  "\\#DECF3FFF", "\\#F17CB0FF", "\\#66E3D9FF", 
                  "\\#00FF7FFF", "\\#E31A1CFF", "\\#FFFF99FF")
      colour <- stringr::str_remove_all(colour, ".*#") %>% 
        paste0("#", .)
    }
    plot_df <- do.call(rbind, lapply(show_pathway, function(x) {
      calcScore(geneset, genelist, x, exponent, fortify = TRUE, 
                org)
    }))
    description_color <- table(plot_df$Description) %>% names()
    names(description_color) <- colour[seq_along(description_color)]
    p1 <- ggplot(plot_df, aes_(x = ~x)) + xlab(NULL) + geom_line(aes_(y = ~runningScore, 
                                                                      color = ~Description), size = 1) + scale_color_manual(values = names(description_color)) + 
      geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + 
      ylab("Enrichment\n Score") + plot_theme(remove_grid = T, 
                                              ...) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                           axis.line.x = element_blank(), legend.position = "top", 
                                                           legend.title = element_blank(), legend.background = element_rect(fill = "transparent"), 
                                                           plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                unit = "cm"))
    panel_size <- c(1.5, 0.5, 1.5)
    i <- 0
    for (term in unique(plot_df$Description)) {
      idx <- which(plot_df$ymin != 0 & plot_df$Description == 
                     term)
      plot_df[idx, "ymin"] <- i
      plot_df[idx, "ymax"] <- i + 1
      i <- i + 1
    }
    p2 <- ggplot(plot_df, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                              ymax = ~ymax, color = ~Description)) + xlab(NULL) + 
      ylab(NULL) + scale_color_manual(values = colour) + 
      plot_theme(remove_legend = T, remove_grid = T, remove_main_text = T, 
                 ...) + scale_y_continuous(expand = c(0, 0)) + 
      theme(plot.margin = margin(t = -0.1, b = 0, unit = "cm"), 
            axis.line.x = element_blank(), axis.ticks = element_blank())
    plot_df$y <- plot_df$genelist
    if (!is.null(show_gene)) {
      select_genes <- data.frame(gene = show_gene)
      select_genes <- merge(select_genes, plot_df, by = "gene")
      select_genes <- select_genes[select_genes$position == 
                                     1, ]
    }
    else {
      select_genes <- data.frame(gene = "no_select_genes")
      select_genes <- merge(select_genes, plot_df, by = "gene")
      select_genes <- select_genes[select_genes$position == 
                                     1, ]
    }
    gene_color <- names(description_color)[description_color %in% 
                                             select_genes$Description]
    if (!"main_text_size" %in% names(lst)) 
      lst$main_text_size <- 8
    p3 <- ggplot(plot_df) + geom_segment(aes_(x = ~x, xend = ~x, 
                                              y = ~y, yend = 0), color = "grey") + geom_bar(data = select_genes, 
                                                                                            aes(x = x, y = y, fill = Description, color = Description), 
                                                                                            position = "dodge", stat = "identity", width = 0.5) + 
      scale_fill_manual(values = gene_color, guide = "none") + 
      scale_color_manual(values = gene_color, guide = "none") + 
      geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + 
      ylab("Ranked list\n metric") + xlab("Rank in ordered dataset") + 
      plot_theme(remove_grid = T, ...) + ggrepel::geom_text_repel(data = select_genes, 
                                                                  aes(x = x, y = y, label = gene, color = Description), 
                                                                  show.legend = FALSE, direction = "x", ylim = c(2, 
                                                                                                                 NA), angle = 90, size = lst$main_text_size/3.5, 
                                                                  box.padding = unit(0.35, "lines"), point.padding = unit(0.3, 
                                                                                                                          "lines")) + theme(plot.margin = margin(t = -0.1, 
                                                                                                                                                                 r = 0.2, b = 0.2, l = 0.2, unit = "cm"))
    pl <- list(p1, p2, p3)
    n <- length(pl)
    pl[[n]] <- pl[[n]] + theme(axis.line.x = element_line(), 
                               axis.ticks.x = element_line())
    p <- cowplot::plot_grid(plotlist = pl, ncol = 1, align = "v", 
                            rel_heights = panel_size)
  }
  if (plot_type == "fgsea") {
    if (is.null(show_pathway)) 
      show_pathway = 3
    geneset_list <- geneset %>% split(.[[1]]) %>% lapply("[[", 
                                                         2)
    gl1 = gl[, 2]
    names(gl1) = gl[, 1]
    genelist = gl1
    fres <- suppressWarnings(fgsea::fgsea(pathways = geneset_list, 
                                          stats = genelist))
    if (is.numeric(show_pathway)) {
      up_path <- fres %>% dplyr::filter(ES > 0) %>% dplyr::arrange(., 
                                                                   pval) %>% dplyr::slice_head(n = show_pathway) %>% 
        dplyr::pull(pathway)
      down_path <- fres %>% dplyr::filter(ES < 0) %>% dplyr::arrange(., 
                                                                     pval) %>% dplyr::slice_head(n = show_pathway) %>% 
        dplyr::pull(pathway)
      all <- c(up_path, down_path)
    }
    else if (any(!show_pathway %in% fres$pathway)) {
      stop(paste0(show_pathway[!show_pathway %in% fres$pathway], 
                  " not in fgsea result!"))
    }
    else {
      all <- show_pathway
    }
    p <- fgsea::plotGseaTable(geneset_list[all], genelist, 
                              fres, render = F) %>% ggplotify::as.ggplot()
  }
  if (plot_type == "ridge") {
    if (is.null(show_pathway)) 
      show_pathway = 3
    if (is.numeric(show_pathway)) {
      if (length(show_pathway) > 1) {
        show_pathway <- gsea_df$ID[show_pathway]
      }
      else {
        show_pathway <- gsea_df$ID[1:show_pathway]
      }
    }
    check_gset_type <- suppressWarnings(ifelse(all(is.na(as.numeric(dplyr::slice_head(geneset, 
                                                                                      n = 2) %>% dplyr::pull(2)))), "symbol", "entrez"))
    if (check_gset_type == "symbol") {
      new_gsea_df <- gsea_df %>% dplyr::filter(ID %in% 
                                                 show_pathway) %>% dplyr::select(ID, dplyr::all_of(stats_metric), 
                                                                                 geneID_symbol) %>% tidyr::separate_rows(geneID_symbol, 
                                                                                                                         sep = "\\/") %>% dplyr::rename(geneID = geneID_symbol)
    }
    else {
      new_gsea_df <- gsea_df %>% dplyr::filter(ID %in% 
                                                 show_pathway) %>% dplyr::select(ID, dplyr::all_of(stats_metric), 
                                                                                 geneID) %>% tidyr::separate_rows(geneID, sep = "\\/")
    }
    logfc <- gsea_list$genelist
    plot_df <- merge(new_gsea_df, logfc, by.x = "geneID", 
                     by.y = "ID") %>% dplyr::select(-geneID)
    term_order <- plot_df %>% dplyr::group_by(ID) %>% dplyr::summarise(logfc2 = sum(logfc)) %>% 
      dplyr::arrange(desc(logfc2))
    plot_df$ID <- factor(plot_df$ID, levels = rev(term_order$ID))
    if (is.null(colour)) {
      colour = c("#E31A1C", "#1F78B4")
    }
    up_color = colour[1]
    down_color = colour[2]
    p <- ggplot(plot_df, aes_string(x = "logfc", y = "ID", 
                                    fill = stats_metric)) + ggridges::geom_density_ridges() + 
      scale_fill_continuous(low = up_color, high = down_color, 
                            name = stats_metric, guide = guide_colorbar(reverse = TRUE)) + 
      ylab(NULL) + xlab("log2(Fold Change)") + guides(fill = guide_colourbar(title = stats_metric_label, 
                                                                             reverse = T)) + plot_theme(remove_grid = T, ...)
  }
  if (plot_type == "bar") {
    gsea_df <- gsea_df %>% dplyr::select(ID, NES, "p.adjust") %>% 
      dplyr::mutate(padj.group = cut(.$p.adjust, breaks = c(-Inf, 
                                                            0.1, Inf), labels = c(1, 0)), nes.group = cut(.$NES, 
                                                                                                          breaks = c(-Inf, 0, Inf), labels = c(0, 1)), 
                    comb.group = paste0(padj.group, nes.group)) %>% 
      dplyr::mutate(group = dplyr::case_when(comb.group == 
                                               "10" ~ "DownStrong", comb.group == "11" ~ "UpStrong", 
                                             comb.group == "00" ~ "DownLight", comb.group == 
                                               "01" ~ "UpLight")) %>% dplyr::arrange(NES) %>% 
      dplyr::mutate(index = 1:dplyr::n())
    if (!"main_text_size" %in% names(lst)) 
      lst$main_text_size <- 8
    if (is.null(colour)) {
      colour <- c("#F04708", "#1476CD", "#F9AF93", "#9BD6F0")
      colour <- stringr::str_remove_all(colour, ".*#") %>% 
        paste0("#", .)
    }
    p <- ggplot(gsea_df, aes(x = index, y = NES, fill = group)) + 
      geom_bar(stat = "identity", width = 0.8) + scale_fill_manual(values = c(UpStrong = colour[1], 
                                                                              DownStrong = colour[2], UpLight = colour[3], DownLight = colour[4])) + 
      scale_x_discrete(expand = expansion(add = 0.5)) + 
      scale_y_continuous(breaks = seq(floor(min(gsea_df$NES)), 
                                      ceiling(max(gsea_df$NES)), ceiling((ceiling(max(gsea_df$NES)) - 
                                                                            floor(min(gsea_df$NES)))/6))) + coord_flip()
    pos_new <- sum(gsea_df$NES > 0)
    neg_nes <- sum(gsea_df$NES < 0)
    if (neg_nes == 0 & pos_new != 0) {
      p <- p + geom_text(data = subset(gsea_df, NES > 0), 
                         aes(x = index, y = 0, label = paste0("  ", ID), 
                             color = padj.group), size = lst$main_text_size/2, 
                         hjust = "outward")
    }
    else if (neg_nes != 0 & pos_new == 0) {
      p <- p + geom_text(data = subset(gsea_df, NES < 0), 
                         aes(x = index, y = 0, label = paste0("  ", ID), 
                             color = padj.group), size = lst$main_text_size/2, 
                         hjust = "inward")
    }
    else {
      if (pos_new <= neg_nes) {
        p <- p + geom_text(data = subset(gsea_df, NES > 
                                           0), aes(x = index, y = 0, label = paste0(ID, 
                                                                                    "  "), color = padj.group), size = lst$main_text_size/2, 
                           hjust = "inward") + geom_text(data = subset(gsea_df, 
                                                                       NES < 0), aes(x = index, y = 0, label = paste0("  ", 
                                                                                                                      ID), color = padj.group), size = lst$main_text_size/2, 
                                                         hjust = "outward")
      }
      else if (pos_new > neg_nes) {
        p <- p + geom_text(data = subset(gsea_df, NES > 
                                           0), aes(x = index, y = 0, label = paste0(ID, 
                                                                                    "  "), color = padj.group), size = lst$main_text_size/2, 
                           hjust = "outward") + geom_text(data = subset(gsea_df, 
                                                                        NES < 0), aes(x = index, y = 0, label = paste0("  ", 
                                                                                                                       ID), color = padj.group), size = lst$main_text_size/2, 
                                                          hjust = "inward")
      }
    }
    p <- p + scale_colour_manual(values = c("black", colour[3])) + 
      labs(x = "", y = "Normalized enrichment score") + 
      ggtitle(paste0("Top ", show_pathway, " pathways in ", cluster," enriched by GSEA using GO database")) + 
      plot_theme(remove_grid = T, remove_legend = T, ...) + 
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5)) + 
      theme(legend.position = "none") +
      theme(axis.title=element_text(size=10))

  }
  suppressMessages(print(p))
}

# @ version of clusterProfiler:: cnetplot
networkPlotGSEA <- function(fgseaResPlot, ranked.genes, cluster, CNET_gene_cutoff = 5){
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
            plot.caption = element_text(size=10))},
    finally = {
      return(cnet)
    })
  return(out)
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
    fgseaRes <- as.data.frame(apply(fgseaRes,2,as.character))
    
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
      b <- plot_fgsea(fgseaRes, cluster, plot_type = "bar", show_pathway = nrow(fgseaRes))
      ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotGO_",cluster,".png"), 
             width = 20, height = 2*show_pathway, units = "cm")
      
      # Save dataframe
      write.csv(fgseaRes, 
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
        b <- plot_fgsea(fgseaRes, cluster, plot_type = "bar", show_pathway = nrow(fgseaRes))
        ggsave(b, file=paste0("TEPA_plots/05_clProf_barplotReactome_",cluster,".png"), width = 20, height = 20, units = "cm")
        
        write.csv(fgseaRes, 
                  file = paste0("TEPA_results/05_GSEAclusterReactome",cluster,".csv"))
      }
    }
  }
}



