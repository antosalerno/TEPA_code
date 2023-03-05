library(msigdbr)
library(scater)
library(Matrix)
library(data.table)
library(dplyr)
library(egg)
library(ggplot2)
library(openxlsx)
library(fgsea)
library(tidyverse)
source("~/Documents/DocumentsTEST/HCV_Data/GSEA/gsea_funcs.R")
data_path = "~/Documents/DocumentsTEST/HCV_Data/GSEA"
data_path2 = "~/Documents/DocumentsTEST/HCV_Data/DE_GENES"

pathways_of_interest = c(
  # "T.CD8.CYTOTOXIC",
  # "CD8_naive",
  # "Central_memory",
  # "Progenitor.exhausted",
  # "Terminally.exhausted",
  # "Tcell.exhaustion.(ZhangCell2017)",
  "PID_IL2_1PATHWAY",
  "PID_CXCR4_PATHWAY",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "PID_IL2_STAT5_PATHWAY",
  "REACTOME_RUNX3_REGULATES_BCL2L11_BIM_TRANSCRIPTION",
  "PID_NFAT_TFPATHWAY",
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_GLYCOLYSIS",
  "PID_AP1_PATHWAY",
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "KEGG_FATTY_ACID_METABOLISM",
  "REACTOME_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "REACTOME_PD_1_SIGNALING",
  "KEGG_PROTEASOME",
  "PID_TCR_PATHWAY" ,
  # "PID_TCR_CALCIUM_PATHWAY",
  # "PID_CD8_TCR_DOWNSTREAM_PATHWAY",
  "KEGG_JAK_STAT_SIGNALING_PATHWAY",
  # "KEGG_PROTEASOME",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",#
  "HALLMARK_G2M_CHECKPOINT",#
  "REACTOME_PD_1_SIGNALING",
  # metabolism 
  "HALLMARK_MTORC1_SIGNALING",
  # "KEGG_OXIDATIVE_PHOSPHORYLATION", #
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_FATTY_ACID_METABOLISM",
  # "HALLMARK_FATTY_ACID_METABOLISM",
  # "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MTORC1_SIGNALING", #
  "HALLMARK_MYC_TARGETS_V1",
  # "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
  # "REACTOME_NEGATIVE_REGULATION_OF_NOTCH4_SIGNALING", #
  # transcription factors 
  "PID_NFAT_TFPATHWAY", #(FOS JUN IFNG and TBX21) 
  "PID_AP1_PATHWAY",
  # "REACTOME_REGULATION_OF_RUNX2_EXPRESSION_AND_ACTIVITY", #
  # "REACTOME_RUNX3_REGULATES_BCL2L11_BIM_TRANSCRIPTION", 
  # functions 
  # "PID_IL12_2PATHWAY", #
  # "PID_IFNG_PATHWAY", #
  # "PID_IL2_STAT5_PATHWAY", 
  "HALLMARK_IL2_STAT5_SIGNALING", #
  "HALLMARK_TGF_BETA_SIGNALING",  #
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", 
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
  # "HALLMARK_INTERFERON_ALPHA_RESPONSE",  #
  "PID_TCR_PATHWAY",
  "CD8.cytotoxic",
  "NK.receptor",
  "NK.like.cytotoxicity",
  # "NK_like_Suppressive",
  # "nk_like_activating",
  # "NKT",
   "NK",
   "NK_COMPLETE",
   "MAIT_CELL",
   "Progenitor.exhausted"
)

# Get gene sets from databases and csv files
if (!exists("msig_data")) {
  msig_data = msigdbr(species = "Homo sapiens")
}
gene_sets = msig_data %>% filter(gs_name %in% pathways_of_interest) %>% split(x = .$gene_symbol, f = .$gs_name)

#custom_gene_sets = read.xlsx(file.path(data_path, "GeneSIgnatures_FINAL_merge.xlsx"), sheet=1, colNames=T)
custom_gene_sets = read.xlsx(file.path(data_path, "GENE_SIGNATURE_CURATED_HCV.xlsx"))
custom_gene_sets <- custom_gene_sets %>% 
  as.list() %>%
  lapply(function(x) x[!is.na(x)])
# custom_gene_sets = custom_gene_sets[names(custom_gene_sets) %in% pathways_of_interest]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "(etienne|stat3)"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "LONG"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "Dendrt"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "Rogue"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "TOO"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "NOO"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "TREG"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %like% "T CELL"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %like% "HCC"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %like% "NSCLC_TILs"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "Th1"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "CD4"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "Stress"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "acrophage"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) %ilike% "nocytes"]
custom_gene_sets = custom_gene_sets[! names(custom_gene_sets) == "ILC"]
#custom_gene_sets <- custom_gene_sets[sapply(custom_gene_sets, function(x) sum(x %in% dge$gene)) > 5] # Keep only gene sets with more than 10 genes im our data
all_gene_sets = c(gene_sets, custom_gene_sets)

all_gene_sets = all_gene_sets[! names(all_gene_sets) %ilike% "NOO"]
all_gene_sets = all_gene_sets[! names(all_gene_sets) %ilike% "Interferon.alpha_beta.signaling"]

# Read in DGE csv
#for (input_csv in c(
  # "DGE_MemLike_CH.csv",
  # "DGE_MemLike_MCRL.csv",
  # "DGE_mcrl_gpr_trajectories.csv",
  # "DGE_New_mcrl_gpr_PAGA_Clusters.csv",
  # "DGE_MemLike.csv",
  # "DGE_integrated_snn_res.1.csv",
 #  "DGE_ch_trajectories_split_EarlyLate.csv",
  # "DGE_ELISPOT_cat.csv",
  # "DGE_New_Ch_PAGA__NOINT_Clusters.csv",
  # "DGE_ch_tr_earlylate.csv",
#  "DGE_mcrl_tr_earlylate.csv",
  # "DGE_CL_tr_earlylate.csv"
#)) {
# input_csv = "DGE_clustOutcome_Stage.csv"
 # input_csv = "DGE_ch_tr_earlylate.csv"
input_csv = "DGE_clustATD_GPR_ALL_new.csv"

  dge = read.csv(file.path(data_path2, input_csv))
  clusters_ordered = as.character(sort(names(table(dge$cluster))))
  
  gsea_results = list()
  for (c in clusters_ordered) {
    dge_of_cluster = dge %>% filter(cluster == c)
    genelist = dge_of_cluster$avg_log2FC
    names(genelist) = dge_of_cluster$gene
    gsea_results[[c]] <- fgsea(pathways=all_gene_sets, stats=genelist, nperm=10000)
  }
  
  # Save GSEA results to CSV
  write.csv(
    apply(do.call(rbind,lapply(names(gsea_results), function(x) data.frame("condition" = x, gsea_results[[x]]))),2,as.character),
    file.path(data_path, gsub("[.]csv", "_GSEA_results.csv", input_csv)),
    row.names = F
  )
  
  # Take only GSEA results which have low adjusted pvalue for at least one condition
  gsea_plot_data = do.call(rbind,lapply(names(gsea_results), function(x) data.frame("condition" = x, gsea_results[[x]])))
  #gsea_plot_data = gsea_plot_data %>% filter(pathway %in% c(gsea_plot_data[gsea_plot_data$padj <= 0.05, "pathway"]))
  gsea_plot_data = gsea_plot_data %>% filter(pathway %in% c(gsea_plot_data[gsea_plot_data$pval <= 0.05, "pathway"]))
  # gsea_plot_data = getTopGeneSets(gsea_results, n = 10000, p_value_cutoff = 0.05)
  if (nrow(gsea_plot_data) > 0) {
    gsea_plot_data = gsea_plot_data %>% complete(condition,pathway) # add NAs for plot
    
    # Clean up pathway names:
    gsea_plot_data$pathway = gsub("HALLMARK_(.*)", "\\1 (HM)", gsea_plot_data$pathway)
    gsea_plot_data$pathway = gsub("[.][.]", " (", gsea_plot_data$pathway)
    gsea_plot_data$pathway = gsub("[.]$", ")", gsea_plot_data$pathway)
    gsea_plot_data$pathway = gsub("[._]+", " ", gsea_plot_data$pathway)
    gsea_plot_data$pathway = gsub("[^A-Za-z]nk ", "NK", gsea_plot_data$pathway)
    gsea_plot_data$pathway = gsub("NK like ", "NK-like", gsea_plot_data$pathway)
    
    # Order pathways by condition
    tmp = as.data.table(gsea_plot_data)
    tmp = tmp[,.SD[which.max(NES)],by=pathway]
    pathway_order = as.character(tmp[order(-tmp$condition, tmp$NES),]$pathway)
    gsea_plot_data$pathway = factor(gsea_plot_data$pathway, levels = pathway_order)
    gsea_plot_data$condition = factor(gsea_plot_data$condition, unique(gsea_plot_data$condition))
    
    
    # plot
    g = ggplot(gsea_plot_data, aes(y=pathway, x=condition)) +
      geom_tile(aes(fill = NES), colour = "white", size=0.25) + 
      theme_classic() +
      labs(x=NULL,y=NULL,title=NULL) +
      scale_fill_gradientn(colours=c("blue", "white", "red"), 
                           # colours=c("blue", "white", "red", "darkred"), 
                           na.value = 'darkgrey',
                           guide=guide_colourbar(),
                           name="NES",
                           limits=c(min(gsea_plot_data$NES, na.rm = TRUE), max(gsea_plot_data$NES, na.rm = TRUE))) + 
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0)) +
      coord_fixed() +
      theme(legend.position = "right",
            panel.border = element_rect(colour='black',size=0.7,fill=NA))
    if (max(nchar(as.character(gsea_plot_data$condition))) >= 3) {
      g = g + theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
    }
    pdf(file.path(data_path, gsub("[.]csv", "_GSEA_heatmap.pdf", input_csv)),
        height = (length(pathway_order)/9 + # for number of genes
          2)*1.3, # for axis labels and title
        width = ((length(unique(gsea_plot_data$condition)) * 0.28) + # for clusters
          0.5 + 0.45 + # For legend + whitespace
          max(nchar(pathway_order)) *0.05) * 1.3 # for gene names
    )
    print(g
          # +
          #  scale_x_discrete(limits = c("Ex", "PD-1int CD127low", "ML", "Mem", "PD-1low CD127low"),
          # labels = c(bquote(bold("T"["EX"])),
          #            bquote(bold("T"["PINT"])),
          #            bquote(bold("T"["ML"])),
          #            bquote(bold("T"["MEM"])),
          #            bquote(bold("T"["EFF"]))))
            )
    dev.off()
  }
#}

  