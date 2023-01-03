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
load("TEPA_results/01_seusetImmune.rda")
#load("TEPA_results/cons_int_anchor_res0.1.rda")

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

# Set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
sctype_scores$type <- c("Macrophages", "CD8+ T cells", "CD4+ T cells", "Platelets", "Plasma Cells", "DCs", 
                        "B cells", "γδ-T cells", "Neutrophils", "NKs","Basophils", "Eosinophils")

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seuset_immune, reduction = "umap",
        label = TRUE, repel = TRUE,
        group.by = "scType")        

#### Intra-cluster DEA with annotated dataset ####

immune.ann <- RenameIdents(seuset_immune,
                           `0` = "Neutrophils",
                           `1` = "Macrophages",
                           `2` = "Eosinophils",
                           `3` = "Cd4+ T-cells",
                           `4` = "Ly6C+ monocytes",
                           `5` = "DCs",
                           `6` = "B-cells",
                           `7` = "Cd8+ T-cells",
                           `8` = "NKs",
                           `9` = "γδ-T cells",
                           `10` = "Basophils",
                           `11` = "Plasma cells")

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
    write.csv(condition.diffgenes, file=paste0("TEPA_results/DEA_",cluster,".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, "TEPA_results/DEA_singleCell.xlsx", rowNames=TRUE)

#### Inter-cluster DEA: get marker genes ####

DefaultAssay(seuset_immune) <- "RNA"

# Find markers for every cluster compared to all remaining cells
immune.markers <- FindAllMarkers(seuset_immune, 
                                 only.pos = FALSE, 
                                 min.pct = 0.5, 
                                 min.diff.pct = 0.5,
                                 logfc.threshold = 0.5, 
                                 test.use="MAST",
                                 latent.vars="orig.ident")

write.csv(immune.markers, "TEPA_results/immune_markers.csv")
#immune.markers <- read.csv("TEPA_results/immune_markers.csv")

immune.markers %>%
  group_by(cluster) %>% 
  filter(p_val_adj<0.05)

# Save results in different excel sheets 

clusters <- c("Neutrophils","Macrophages","Eosinophils","Cd4+ T-cells","Ly6C+ monocytes",
              "DCs","B-cells", "Cd8+ T-cells","NKs","γδ-T cells","Basophils","Plasma cells")

wb <- createWorkbook()
for(c in 1:12){
  cluster = immune.markers[immune.markers$cluster == c-1,]
  addWorksheet(wb, c)
  writeData(wb, c, cluster[,2:ncol(cluster)], colNames = TRUE)
}
saveWorkbook(wb, file="TEPA_results/markers_singleCell.xlsx", overwrite = TRUE)

