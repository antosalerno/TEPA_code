
library(sumer)
file.copy(system.file("data", "sample.tgz", package="sumer"), ".")
untar("sample.tgz")

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
fgseaResByType = read.csv("TEPA_results/05_GSEAbyType.csv", sep = ",")

write.table(fgseaResByType[,2:3], 
            file ="05_GSEAbyType.txt", sep = "",
            row.names = F, col.names = F, quote = FALSE)

sumer("config.json", "TEPA_results")

library(DropletUtils)
write10xCounts(x = seuset_immune@assays$RNA@counts, path = "TEPA_data/expmatRhapsody.mtx")

expData <- Matrix(as.matrix(seuset_immune@assays$RNA@data), sparse = TRUE)
writeMM(obj = expData, file = "TEPA_data/expmatRhapsody.mtx")

####
devtools::install_github("aertslab/SCENIC@v1.1.2")
devtools::install_github("aertslab/AUCell")
devtools::install_github("aertslab/RcisTarget")
devtools::install_github("aertslab/GRNBoost")

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"), force = T)
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"), force = T)
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"), force = T)
# To export/visualize in http://scope.aertslab.org
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")

#dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
featherURL = "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/"
download.file(featherURL, destfile="TEPA_data/cisTarget_databases/mm9.refseq_r45.mc9nr.gene_based.feather") # saved in current dir

exprMat <- as.data.frame(as.matrix(GetAssayData(seuset_immune, assay="integrated", slot="data")))
cellInfo <- data.frame(CellType=Idents(seuset_immune))

cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)
cbind(table(cellInfo$CellType))
dir.create("TEPA_data/int")
saveRDS(cellInfo, file="TEPA_data/int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)

colVars <- list(CellType=c("Macrophages"="forestgreen", 
                           "Cd4+ Naive T cells"="darkorange", 
                           "B cells"="magenta4", 
                           "γδ-T cells" ="hotpink", 
                           "Cd4+ Regulatory T cells"="red3", 
                           "Cd8+ Naive T cells"="skyblue", 
                           "Cd4+ Memory T cells" = "red",
                           "Dendritic cells" ="magenta",
                           "Cd8+ NkT-like cells" ="orange",
                           "Progenitor cells"="darkblue",
                           "Natural killer  cells"="lightgreen",
                           "Neutrophils"="brown3",
                           "Eosinophils" ="blue1",
                           "Basophils" = "pink"
                           ))

colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="TEPA_data/int/colVars.Rds")
plot.new(); 
legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

library(SCENIC)
devtools::install_github("wesm/feather/R")
org <- "mgi" # or hgnc, or dmel
dbDir <- "TEPA_data/cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC TEPA analysis" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, 
                                  dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 


### Prepare data for Cytoscape ####
cluster = "Macrophages"
res <- read.csv(paste0("TEPA_results/02_DEAclusterMAST",cluster,".csv"), sep=",")
ranked.genes<- deframe(res %>%
                         dplyr::filter(p_val_adj < 0.1) %>%
                         arrange(desc(avg_log2FC)) %>% 
                         dplyr::select(X, avg_log2FC))
ranked.genes <- ranked.genes[!is.na(names(ranked.genes))]
write.table(ranked.genes, 
            file = paste0("TEPA_results/05_geneRanked",paste(cluster, collapse = "_"),".txt"),
            sep="\t", quote = F, row.names = T, col.names = F)

enriched.pos <- names(ranked.genes[ranked.genes < 0])
enriched.neg <- names(ranked.genes[ranked.genes < 0])
write.table(enriched.pos, 
          file = paste0("TEPA_results/05_enrichedPOS",paste(cluster, collapse = "_"),".txt"),
          sep="", quote = F, row.names = F, col.names = F)
write.table(enriched.neg, 
            file = paste0("TEPA_results/05_enrichedNEG",paste(cluster, collapse = "_"),".txt"),
            sep="", quote = F, row.names = F, col.names = F)












