#### Automatic clustering annotation with scType and DEA ##
## author: Antonietta Salerno
## date: 16/12/2022

library(openxlsx)
library(HGNChelper)
library("Seurat")
library("writexl")
library('limma')
library(dplyr)
library("MAST")
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}


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
                              es.max.cl = sort(rowSums(es.max[,rownames(seuset_immune@meta.data[seuset_immune@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
                              head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seuset_immune@meta.data$seurat_clusters==cl)), 10)
                              }))

sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
sctype_scores = sctype_scores[order(sctype_scores$cluster),]

# Manual curation of scType annotation
sctype_scores[3, "type"] <- "Cd4+ Naive T cells"
sctype_scores[4, "type"] <- "Cd4+ Naive T cells"
sctype_scores[5, "type"] <- "Cd4+ Memory T cells"
sctype_scores[6, "type"] <- "Cd4+ Regulatory T cells"
sctype_scores[7, "type"] <- "Cd8+ Naive T cells"
sctype_scores[8, "type"] <- "Cd8+ NkT-like cells"
sctype_scores[10, "type"] <- "Macrophages" # or macrophages
sctype_scores[11, "type"] <- "Dendritic cells" # or DCs
sctype_scores[13, "type"] <- "B cells"
sctype_scores[15, "type"] <- "HLA-expressing cells"

seuset_immune@meta.data$scType = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seuset_immune@meta.data$scType[seuset_immune@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

Idents(seuset_immune) <- seuset_immune@meta.data$scType

# We decide to remove the cluster HLA-expressing cells since it is probably ambient RNA -> mixture of B cells and macrophage markers
seuset_immune <- seuset_immune[,!seuset_immune$scType == "HLA-expressing cells"]
levels(Idents(seuset_immune)) # now 14 clusters rather than 15

png("TEPA_plots/02_umapAnn.png", w = 3000, h = 3000, res = 300)
p <- DimPlot(object = seuset_immune, pt.size = 0.5, reduction = 'umap', ncol = 1,
             group.by = "scType", label = FALSE) +
  ggtitle("Cell types in NB Control and Treated samples") +
  theme(plot.title = element_text(hjust = 0.5)) 
LabelClusters(p, id = "scType", size = 5, repel = T,  box.padding = 1)
dev.off()

### Bar Plot with proportions
pt <- table(Idents(seuset_immune), seuset_immune$condition)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
getPalette = colorRampPalette(brewer.pal(11, "PiYG"))
colorCount = length(unique(pt$Var1))

png(paste0("TEPA_plots/02_condAnnClusterFreq.png"), w=2500,h=2500, res=300)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  geom_bar(stat = "identity")+
  xlab("Condition") +
  ylab("Cell type") +
  ggtitle("Cell types in TEPA vs Control") +
  scale_fill_manual(values = getPalette(colorCount)) +
  theme(legend.title = element_blank())
dev.off()

###  Bar Plot with percentage
# Choose color palette (brewer.pal.info)
pt <- table(Idents(seuset_immune), seuset_immune$condition)
data_percentage<- apply(pt, 2, function(x){as.numeric(x)*100/sum(x,na.rm=T)})

# Make a stacked barplot--> it will be in %!
png(paste0("TEPA_plots/02_condAnnClusterPerc.png"), w=2500,h=2500, res=300)
par(mar = c(5.1, 5.1, 4.1, 12))
barplot(data_percentage, col=getPalette(colorCount), border="white", 
        ylab="Percentage of cells per cell type",
        main = "Cell types in TEPA vs Control",
        legend = rownames(pt),
        args.legend = list(x = "topright",inset = c(-0.45, 0)))
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
clusters = levels(Idents(seuset_immune))

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
                                       logfc.threshold = 0.25, 
                                       only.pos = FALSE, verbose = FALSE,
                                       #latent.vars="orig.ident",
                                       min.cells.feature = 1, min.cells.group = 1, min.pct = 0,
                                       test.use="MAST")
    sheets[[cluster]] <- as.data.frame(condition.diffgenes)
    
    # Needed for plotting
    write.csv(condition.diffgenes, file=paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"))
  })
}
# Needed for manual curation
openxlsx::write.xlsx(sheets, "TEPA_results/02_DEA_TEPA_MAST.xlsx", rowNames=TRUE)




