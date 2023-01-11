#### Automatic clustering annotation with scType and DEA ##
## author: Antonietta Salerno
## date: 16/12/2022

library(openxlsx)
library(HGNChelper)
library("Seurat")
library("writexl")
library('limma')
library(dplyr)


setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
# load("TEPA_results/01_seusetImmune.rda")
seuset_immune <- LoadH5Seurat("TEPA_results/01_seusetImmune.h5Seurat")


source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = "Immune system"

gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = seuset_immune[["integrated"]]@scale.data,
                      scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

####  Get scType scores by cluster ####
cL_results = do.call("rbind", 
                     lapply(unique(seuset_immune@meta.data$seurat_clusters), 
                                     function(cl){
                                       es.max.cl = sort(rowSums(es.max[ ,rownames(seuset_immune@meta.data[seuset_immune@meta.data$seurat_clusters==cl, ])]),decreasing = !0)
                                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuset_immune@meta.data$seurat_clusters==cl)), 10)
                                       }))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

write.csv(cL_results, file=paste0("TEPA_results/02_scTypeAnn.csv"))

# Set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

# Manual curation of scType annotation
sctype_scores <- sctype_scores[-4,]
sctype_scores$type <- c("Dendritic cells", "Naive CD4+ T cells", "B cells", "γδ-T cells", "Memory CD4+ T cells", "CD8+ T cells", "Progenitor cells",
                        "Natural killer cells", "Neutrophils", "Eosinophils", "Basophils")

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(seuset_immune) <- seuset_immune@meta.data$scType

png("TEPA_plots/02_umapExploreAnn.png", w = 6000, h = 3000, res = 400)
DimPlot(object = seuset_immune, pt.size = 0.0000005, reduction = 'umap', ncol = 2,
        group.by = c("orig.ident", "scType"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells"))
dev.off()

png("TEPA_plots/02_umapClust.png", w = 6000, h = 3000, res = 300)
DimPlot(object = seuset_immune, pt.size = 0.5, reduction = 'umap', ncol = 2,
        group.by = c("scType"), split.by= "condition",label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pt <- table(Idents(seuset_immune), seuset_immune$condition)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

png(paste0("TEPA_plots/02_condAnnCluster.png"), w=2500,h=2500, res=300)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Condition") +
  ylab("Cell type") +
  ggtitle("Cell types in TEPA vs Control") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())
dev.off()

SaveH5Seurat(seuset_immune, filename = "TEPA_results/02_immuneAnn.h5Seurat", overwrite = TRUE)

#### Inter-cluster DEA: get marker genes ####

DefaultAssay(seuset_immune) <- "RNA"

# Find markers for every cluster compared to all remaining cells
immune.markers <- FindAllMarkers(seuset_immune, 
                                 only.pos = FALSE, 
                                 min.pct = 0.5, 
                                 min.diff.pct = 0.2,
                                 logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="orig.ident")

write.csv(immune.markers, "TEPA_results/02_DEA_clusterMarkers.csv")

# Save results in different excel sheets 
clusters = unique(Idents(seuset_immune))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = immune.markers[immune.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/02_DEA_clusterMarkers.xlsx", overwrite = TRUE)


#### Intra-cluster DEA with annotated dataset - Treatment vs Control ####

seuset_immune$celltype.tepa <- paste(Idents(seuset_immune), seuset_immune$condition, sep = "_")
seuset_immune$celltype <- Idents(seuset_immune)
Idents(seuset_immune) <- "celltype.tepa"
DefaultAssay(seuset_immune) <- "RNA"

sheets <- list()
for (cluster in unique(seuset_immune$celltype)){
  try({
    ident1 <- paste0(cluster,"_Control")
    ident2 <- paste0(cluster,"_Treatment")
    condition.diffgenes <- FindMarkers(seuset_immune, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       min.pct = 0.25, logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE,
                                       test.use="MAST", latent.vars="orig.ident")
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/02_DEAcluster",cluster,".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, "TEPA_results/02_DEA_TEPA.xlsx", rowNames=TRUE)

SaveH5Seurat(seuset_immune, filename = "TEPA_results/02_immuneAnn.h5Seurat", overwrite = TRUE)





