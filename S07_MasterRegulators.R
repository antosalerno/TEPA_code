library(corto)
library("readr")
library(Seurat)
library(SeuratDisk)
#load(system.file("extdata","inmat.rda",package="corto"))

setwd("~/OneDrive - Childrens Cancer Institute Australia/OrazioLab")
source("TEPA_code/supportFunctions.R")

# Load raw data from neutrophils mice tumors and create seurat object
counts <- read.table("TEPA_data/neutrophils_tumors/GSM7021720_KP19_rep1_raw_counts.tsv", 
                     header = TRUE,  sep ="\t", row.names = "cell_id")
C1 <- CreateSeuratObject(counts = t(counts), 
                         project="Neutrophils", min.cells = 3, min.features = 500)

counts <- read.table("TEPA_data/neutrophils_tumors/GSM7021721_KP19_rep2_raw_counts.tsv",
                   header = TRUE,  sep ="\t", row.names = "cell_id")
C2 <- CreateSeuratObject(counts = t(counts), 
                         project="Neutrophils", min.cells = 3, min.features = 500)

counts <- read.table("TEPA_data/neutrophils_tumors/GSM7021724_aCD40_rep1_raw_counts.tsv",
                     header = TRUE,  sep ="\t", row.names = "cell_id")
T1 <- CreateSeuratObject(counts = t(counts),
                         project="Neutrophils", min.cells = 3, min.features = 500)

counts <- read.table("TEPA_data/neutrophils_tumors/GSM7021725_aCD40_rep2_raw_counts.tsv",
                     header = TRUE,  sep ="\t", row.names = "cell_id")
T2 <- CreateSeuratObject(counts = t(counts),
                         project="Neutrophils", min.cells = 3, min.features = 500)

# Create list object to be merged in a large Seurat object
seuset <- list(C2 = C2,
               T1 = T1, T2 = T2)

combined <- merge(
  x = C1,
  y = seuset,
  add.cell.ids = c("C1",'C2',"T1","T2")
)

inmat <- combined

inmat <- ScaleData(inmat)

SaveH5Seurat(inmat, filename = "TEPA_results/S07_scaledNeutrophilsInmat.h5Seurat", overwrite = TRUE)
#inmat <- LoadH5Seurat("TEPA_results/S07_scaledNeutrophilsInmat.h5Seurat")

inmat_counts <- inmat@assays$RNA@scale.data

# centroids are a list of genes of interest 
regulon<-corto(inmat_counts,centroids=copper_genes,nbootstraps=10,p=1e-30,nthreads=2)
save(regulon, file = "TEPA_data/regulon_NeutrophilsInmatCopper.rda")

# Get expression matrices of interest: our neutrophils tepa and control
seuset_immune <- LoadH5Seurat("TEPA_results/S04_immuneDiff.h5Seurat")
neutroCells <- subset(seuset_immune, scType == "Neutrophils")

assayNeutroTEPA <- as.matrix(subset(neutroCells, condition == "Treatment")@assays$RNA@scale.data)
assayNeutroCTRL <- as.matrix(subset(neutroCells, condition == "Control")@assays$RNA@scale.data)

# Predict Master regulators from gene co-expression network
predicted<-mra(assayNeutroTEPA,assayNeutroCTRL,regulon)
save(predicted, file = "TEPA_results/S07_predictedNeutrophils_MRA.rda")

png("TEPA_plots/MRA_neutrophils.png", h = 8000, w = 3000, res = 200)
mraplot(predicted, mrs = 18, pthr = 0.05, title = "Master TFs enriched in TEPA-treated Neutrophils")
dev.off()

#### MRA with tumor ####

load("TEPA_data/Regulon_NBL.rda") # find mouse
seuset_tumor <- LoadH5Seurat("TEPA_results/S05_seusetTumorClu.h5Seurat")

assayTumorTEPA <- as.matrix(subset(seuset_tumor, condition == "Treatment")@assays$integrated@scale.data)
assayTumorCTRL <- as.matrix(subset(seuset_tumor, condition == "Control")@assays$integrated@scale.data)

# Predict Master regulators from gene co-expression network
predicted<-mra(assayTumorTEPA,assayTumorCTRL,regulon)
save(predicted, file = "TEPA_results/S07_predictedTumor_MRA.rda")

png("TEPA_plots/S07_tumor_MRA.png", h = 8000, w = 3000, res = 200)
mraplot(predicted, mrs = 18, pthr = 0.05, title = "Master TFs enriched in TEPA-treated Tumor")
dev.off()


# plotta scores tutti insieme e expression plot dei geni per vedere se sono nei neutrofili o altrove

# differential expression cellule mieloidi tepa vs control -> top ones cerca copper-dependent

# fallo col tumore e poi con NK e Cd8+



