# Tumor vs Immune 

seuset_class <- LoadH5Seurat("TEPA_results/00_seusetClass.h5Seurat")
Idents(seuset_class) <- "class"

seuset_class <- seuset_class[,!seuset_class$class == "noClass"]

png("TEPA_plots/03_VEGF_violinTumorCondition.png", h = 2000, w = 3500, res = 200)
VlnPlot(seuset_class, features =c("Vegfa", "Vegfb", "Vegfc", "Vegfd"), ncol = 4, split.by="condition", pt.size = 0.000005)
dev.off()

png("TEPA_plots/03_VEGF_violinTumor.png", h = 2000, w = 3500, res = 200)
VlnPlot(seuset_class, features =c("Vegfa", "Vegfb", "Vegfc", "Vegfd"), ncol = 4, pt.size = 0.000005)
dev.off()

