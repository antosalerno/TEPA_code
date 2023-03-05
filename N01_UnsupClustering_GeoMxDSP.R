#### Unsupervised Analysis using GeoMxTools ####
## author: Antonietta Salerno
## date: 01/03/2022

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
library(umap)
library(Rtsne)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_data , elt = "q_norm"))),  
       config = custom_umap)
pData(protocolData(target_data))[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
ggplot(pData(protocolData(target_data)),
       aes(x = UMAP1, y = UMAP2, color = Group, shape = class)) +
  geom_point(size = 3) +
  theme_bw()


# run tSNE
set.seed(42) # set the seed for tSNE as well
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_data , elt = "q_norm"))),
        perplexity = ncol(target_data)*.15)
pData(target_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]
ggplot(pData(target_data),
       aes(x = tSNE1, y = tSNE2, color = region, shape = class)) +
  geom_point(size = 3) +
  theme_bw()

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
#>   CAMK2N1    AKR1C1      AQP2     GDF15       REN 
#> 0.5886006 0.5114973 0.4607206 0.4196469 0.4193216

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

# convert test variables to factors
pData(protocolData(target_data))$Infiltration <- 
  factor(pData(protocolData(target_data))$Infiltration)

pData(protocolData(target_data))$Group <- 
  factor(pData(protocolData(target_data))$Group)

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

# Try to run T vs C given Infiltration >> then try with seurat

### 2- Disease vs Healthy ###

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()

for(inf in c("T", "F")) {
  ind <- pData(protocolData(target_data))$Infiltration == T
  mixedOutmc <-
    mixedModelDE(target_data[, ind],
                 elt = "log_q",
                 modelFormula = ~ Group + (1 | Infiltration),
                 groupVar = "Group",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- inf
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}

kable(subset(results, Gene %in% goi & Subset == "T"), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)

# Graph results
png("TEPA_plots/N01_Group_DEA.png", h = 2000, w = 2500, res = 300)
ggplot(results,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in Control <- log2(FC) -> Enriched in Treatment",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") 
dev.off()

## Plotting genes of interest

kable(subset(results, Gene %in% c('PDHA1','ITGB1')), row.names = FALSE)

# show expression for a single target: PDHA1
ggplot(pData(target_data),
       aes(x = region, fill = region,
           y = assayDataElement(target_data["PDHA1", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "PDHA1 Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~class) +
  theme_bw()

glom <- pData(target_data)$region == "glomerulus"

# show expression of PDHA1 vs ITGB1
ggplot(pData(target_data),
       aes(x = assayDataElement(target_data["PDHA1", ],
                                elt = "q_norm"),
           y = assayDataElement(target_data["ITGB1", ],
                                elt = "q_norm"),
           color = region)) +
  geom_vline(xintercept =
               max(assayDataElement(target_data["PDHA1", glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_hline(yintercept =
               max(assayDataElement(target_data["ITGB1", !glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = "PDHA1 Expression", y = "ITGB1 Expression") +
  facet_wrap(~class)


# Heatmap of significant genes

# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `FDR` < 0.001)$Gene)
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
         annotation_col = pData(target_data)[, c("region", "class")])


