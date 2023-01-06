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
#sctype_scores$type <- c("Macrophages", "CD8+ T cells", "CD4+ T cells", "Platelets", "Plasma Cells", "DCs", "B cells", "γδ-T cells", "Neutrophils", "NKs","Basophils", "Eosinophils")

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

png("TEPA_plots/02_umapExploreAnn.png", w = 9000, h = 3500, res = 400)
DimPlot(object = seuset_immune, pt.size = 0.0005, reduction = 'umap', ncol = 3,
        group.by = c("orig.ident", "condition", "scType"), label = TRUE) +
  ggtitle(paste(as.character(nrow(seuset_immune@meta.data)), " cells"))
dev.off()

#### Inter-cluster DEA: get marker genes ####

immune.ann <- RenameIdents(seuset_immune,
                           `0` = "Intermediate monocytes",
                           `1` = "Eosinophils",
                           `2` = "ISG expressing immune cells",
                           `3` = "Naive CD4+ T cells",
                           `4` = "Ly6C monocytes",
                           `5` = "Naive CD8+ T cells",
                           `6` = "Memory CD4+ T cells",
                           `7` = "Dendritic cells",
                           `8` = "NKs",
                           `9` = "Progenitor cells",
                           `10` = "γδ-T cells",
                           `11` = "B-cells")

DefaultAssay(immune.ann) <- "RNA"

# Find markers for every cluster compared to all remaining cells
immune.markers <- FindAllMarkers(immune.ann, 
                                 only.pos = FALSE, 
                                 # min.pct = 0.5, 
                                 min.diff.pct = 0.2,
                                 logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="orig.ident")

write.csv(immune.markers, "TEPA_results/02_DEA_clusterMarkers.csv")

# Fast visualisation of markers by cluster
immune.markers %>%
  group_by(cluster) %>% 
  filter(p_val_adj<0.05)

# Save results in different excel sheets 
clusters = unique(Idents(immune.ann))

wb <- createWorkbook()
for(c in 1:length(clusters)){
  cluster = immune.markers[immune.markers$cluster == clusters[c],]
  addWorksheet(wb, as.character(clusters[c]))
  writeData(wb, as.character(clusters[c]), cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/02_DEA_clusterMarkers.xlsx", overwrite = TRUE)


#### Intra-cluster DEA with annotated dataset - Treatment vs Control ####

immune.ann$celltype.tepa <- paste(Idents(immune.ann), immune.ann$condition, sep = "_")
immune.ann$celltype <- Idents(immune.ann)
Idents(immune.ann) <- "celltype.tepa"
DefaultAssay(immune.ann) <- "RNA"

sheets <- list()
for (cluster in unique(immune.ann$celltype)){
  try({
    ident1 <- paste0(cluster,"_Control")
    ident2 <- paste0(cluster,"_Treatment")
    condition.diffgenes <- FindMarkers(immune.ann, 
                                       ident.1 = ident1, ident.2 = ident2,
                                       min.pct = 0.25, logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE)
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    # write.csv(condition.diffgenes, file=paste0("TEPA_results/02_DEAcluster",cluster,".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, "TEPA_results/02_DEA_TEPA.xlsx", rowNames=TRUE)

SaveH5Seurat(immune.ann, filename = "TEPA_results/02_immuneAnn.h5Seurat", overwrite = TRUE)





