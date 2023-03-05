#### Pathway Enrichment Analysis ####
## author: Antonietta Salerno
## date: 19/12/2022


library("fgsea")
library(dplyr)
library(devtools)
library(msigdbr)
library("tibble")
library("Seurat")
library("SeuratDisk")
library(readr)
library(stringr)
library("org.Mm.eg.db", character.only = TRUE)
library(clusterProfiler)
# library("pathview")
# library(enrichplot)
library(ggplot2)
library(ggcharts)
library("ReactomePA")


setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")
seuset_immune <- LoadH5Seurat("TEPA_results/S03_immuneDiff.h5Seurat")
clusters = unique(seuset_immune@meta.data$scType)


###  Bar Plot with percentage
# Choose color palette (brewer.pal.info)
pt <- table(Idents(seuset_immune), seuset_immune$condition)
data_percentage<- apply(pt, 2, function(x){as.numeric(x)*100/sum(x,na.rm=T)})

# Make a stacked barplot--> it will be in %!
png(paste0("TEPA_plots/S02_condAnnClusterPerc.png"), w=2500,h=2500, res=300)
par(mar = c(5.1, 5.1, 4.1, 12))
barplot(data_percentage, col=getPalette(colorCount), border="white", 
        ylab="Percentage of cells per cell type",
        main = "Cell types in TEPA vs Control",
        legend = rownames(pt),
        args.legend = list(x = "topright",inset = c(-0.45, 0)))
dev.off()

### fgsea ####

# 1 - GO database 

dbGO <- godata("org.Mm.eg.db", keytype = "SYMBOL", ont="BP")
dbGOAnn <- dbGO@geneAnno
dbGOAnn <- dbGOAnn[,-c(3:4)]
dbGOAnn <- as.data.frame(dbGOAnn)
dbGOAnn <- cbind(dbGOAnn,get_names(dbGOAnn$GO))
dbGOAnn <- dbGOAnn %>% distinct()
fgsea_sets<- dbGOAnn %>% split(x = .$SYMBOL, f = .$go_name)

gseaRES(clusters, fgsea_sets)

# 2 - MsigDBr database
# Immunological signatures

sets1 <- read.gmt("TEPA_data/mh.all.v2022.1.Mm.symbols.gmt") # Mouse hallmark
#sets2 <- read.gmt("TEPA_data/m2.cp.v2022.1.Mm.symbols.gmt") # Mouse curated canonical pathways M2-CP: BioCarrta, Reactome, WikiPathways
#sets3 <- read.gmt("TEPA_data/m5.all.v2022.1.Mm.symbols.gmt") # Mouse ontology gene sets M5: GO and tumor phenotype oncology
sets4 <- read.gmt("TEPA_data/m2.cp.reactome.v2022.1.Mm.symbols.gmt") # Mouse Reactome subset of Canonical pathways
#sets5 <- msigdbr(species = "Mus musculus", category = "C7")

# tutte quelle che hai fatto più C7

sets1$term <- as.character(sets1$term)
sets4$term <- as.character(sets4$term)
#sets5$term <- as.character(sets5$term)

sets1 <- sets1 %>% split(x = .$gene, f = .$term)
sets4 <- sets4 %>% split(x = .$gene, f = .$term)
sets5 <- sets5 %>% split(x = .$gene_symbol, f = .$gs_name)

fgsea_sets <- append(sets1, sets4)
#fgsea_sets <- append(fgsea_sets, sets5)

gseaRES(clusters, fgsea_sets = fgsea_sets, save = "Full")

#gseaPlotRes(clusters)

# 3 - Reactome database
library("AnnotationDbi")

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

### Plot all clusters' results in a network ####

cell_types = list()
for (cluster in clusters){
  res <- read.csv(paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"), sep=",")
  ranked.genes<- deframe(res %>%
                           dplyr::filter(p_val_adj < 0.05) %>%
                           dplyr::filter(avg_log2FC > 0.2) %>%
                           arrange(desc(avg_log2FC)) %>% 
                           dplyr::select(X, avg_log2FC))
  if (length(ranked.genes) > 0){
    # Convert into Entrez ID
    entrezID <- AnnotationDbi::select(org.Mm.eg.db, keys=names(ranked.genes), columns='ENTREZID', keytype='SYMBOL')
    ranked.genesENTREZ <- ranked.genes
    names(ranked.genesENTREZ) <- entrezID$ENTREZID
    ranked.genesENTREZ <- ranked.genesENTREZ[!is.na(names(ranked.genesENTREZ))]
    
    # Append gene set to the cell types database
    cell_types[cluster] <- list(names(ranked.genesENTREZ))
    
  }
}

gseaByCellType <- compareCluster(cell_types, fun="enrichGO", ont = "BP",
                                 OrgDb = "org.Mm.eg.db", pvalueCutoff=0.05)

lEd <- gseaByCellType@compareClusterResult$geneID 
for (e in 1: length(lEd)){
  entrezID <- AnnotationDbi::select(org.Mm.eg.db, keys=unlist(str_extract_all(lEd[e],"[[:digit:]]{4,}")),
                                    columns='SYMBOL', keytype='ENTREZID')
  gseaByCellType@compareClusterResult$geneID[e] <- paste(entrezID$SYMBOL, collapse="/")
  
}

df_gseaCT <- gseaByCellType@compareClusterResult

cnet <- cnetplot(gseaByCellType, showCategory=5) # convert genes to gene names again-> remember is functional enrichment not gsea so take only positively enriched genes
ggsave(cnet, file=paste0("TEPA_plots/05_compareClusters.png"), width = 40, height = 40, units = "cm")

