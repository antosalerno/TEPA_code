#### Unsupervised Analysis using GeoMxTools ####
## author: Antonietta Salerno
## date: 01/03/2022

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(umap)
library(Rtsne)

N1 <- c("S100a8", "S100a9", "Icam1", "Fas", "Tnf", "Isg15", "Isg20", 
        "Ccl3", "Ccl4",  "Cxcl2", "Cxcl3", "Cebpb", "Il1a" ,"Il1b", "Il1r2", "Il1rn", "Il6ra", "Il15",
        "Stat3", "Hif1a", "Ifitm1","Ifitm2",  "Ifitm3", "Ifitm6",  
        "Acod1", "Myd88", "Prkcd", "Mmp8", "Mmp9","Retnlg", "Arg2")

N2 <- c("Tgfb1", "Tgfb1i1","Tgfb2", "Tgfb3", "Ccl2",  "Ccl17","Cxcl14", 
        "Cxcl15",  "Il1r1", "Il2", "Il17a", "Mpo", "Slc27a2", "Arg1", 
        "Mrc1", "Chil3", "Elane", "Ctsg", "Retnla")


### Clustering high CV genes: most variable genes ###

library(pheatmap)  # for pheatmap
# create a log2 transform of the data for analysis
assayDataElement(object = target_data, elt = "log_q") <-
  assayDataApply(target_data, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_data,
                         elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]
#>  Kap     Aldob      Pck1       Alb    Cyp4b1 
# 0.5291805 0.4816749 0.4736630 0.4308604 0.4062199

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]

png("TEPA_plots/N01_heatmapGOI.png", h = 3000, w = 2500, res = 300)
pheatmap(assayDataElement(target_data[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(protocolData(target_data))[, c("Group", "Infiltration")])
dev.off()

#### Differential expression analysis ####
# 1- Within Slide Analysis: Glomeruli vs Tubules -> it will be infiltrated vs non-infiltrated
library(forcats)
# convert test variables to factors
pData(protocolData(target_data))$Infiltration <- 
  fct_rev(factor(pData(protocolData(target_data))$Infiltration))

pData(protocolData(target_data))$Group <- 
  fct_rev(factor(pData(protocolData(target_data))$Group))

assayDataElement(object = target_data, elt = "log_q") <-
  assayDataApply(target_data, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
tumor <- pData(protocolData(target_data))$Tissue == "Tumor"
mixedOutmc <-
  mixedModelDE(target_data[, tumor],
                 elt = "log_q",
                 modelFormula = ~ Infiltration + (1 + Infiltration | Group),
                 groupVar = "Infiltration",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
# format results as data.frame
results <- do.call(rbind, mixedOutmc["lsmeans", ])
tests <- rownames(results)
results <- as.data.frame(results)
results$Contrast <- tests

# use lapply in case you have multiple levels of your test factor to
# correctly associate gene name with it's row in the results table
results$Gene <- 
  unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
results$FDR <- p.adjust(results$`Pr(>|t|)`, method = "fdr")
results <- results[, c("Gene", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]

# Genes associated with infiltration given treatment that are in our N1 list (not really significant)
results[results$Gene %in% N1,] %>%
  group_by(Gene) %>% arrange(desc(Estimate)) %>%
  filter(`Pr(>|t|)` < 0.1) # Ifitm3 and then Isg20

# Genes associated with infiltration given treatment that are actually significant
resultsInf <- results %>% group_by(Gene)  %>% 
  #filter(`Pr(>|t|)` < 0.05)
  filter(FDR < 0.1) %>%
  arrange(desc(Estimate)) 

write.csv(resultsInf, "TEPA_results/N01_resultsInf.csv")

resultsInfFull <- results %>% group_by(Gene)  %>% 
  filter(`Pr(>|t|)` < 0.05)
  #filter(FDR < 0.1) %>%
  arrange(desc(Estimate)) 

kable(subset(results, Gene %in% N1), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)

# Heatmap of significant genes

# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `FDR` < 0.1)$Gene)
png("TEPA_plots/N01_heatmapGOI.png", h = 3000, w = 2500, res = 300)

pheatmap(log2(assayDataElement(target_data[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(protocolData(target_data))[, c("Infiltration", "Group")])
dev.off()



