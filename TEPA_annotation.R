library("dplyr")
library("writexl")

load(file = "000_all.markers.immune.rda")
print(all.markers %>% group_by(cluster) 
      %>% top_n(5, avg_log2FC) %>% data.frame)

### how many markers per cluster
all.markers %>% group_by(cluster) %>% summarize(n_p05 = sum(p_val_adj<0.05),
                                                npe10 = sum(p_val_adj<1e-10),
                                                npe100 = sum(p_val_adj<1e-100),
                                                np0 = sum(p_val_adj==0))

### write file with markers having adj p < 0.05
top05 <- all.markers %>% group_by(cluster) %>% filter(p_val_adj<0.05)

### create a single file
marker.excel.pages <- list('clust0' = top05 %>% filter(cluster==0),
                           'clust1' = top05 %>% filter(cluster==1),
                           'clust2' = top05 %>% filter(cluster==2),
                           'clust3' = top05 %>% filter(cluster==3),
                           'clust4' = top05 %>% filter(cluster==4),
                           'clust5' = top05 %>% filter(cluster==5),
                           'clust6' = top05 %>% filter(cluster==6),
                           'clust7' = top05 %>% filter(cluster==7),
                           'clust8' = top05 %>% filter(cluster==8),
                           'clust9' = top05 %>% filter(cluster==9),
                           'clust10' = top05 %>% filter(cluster==10),
                           'clust11' = top05 %>% filter(cluster==11),
                           "clust12"= top05 %>% filter(cluster==12),
                           "clust13"= top05 %>% filter(cluster==13),
                           "clust14"= top05 %>% filter(cluster==14),
                           "clust15"= top05 %>% filter(cluster==15))

top05 <- top05[order(top05$avg_log2FC,decreasing = T),]
write_xlsx(top05, "000_topgenes_immune_cluster.xlsx")

cl0 <- top05[top05["cluster"] == 0,]
cl1 <- top05[top05["cluster"] == 1,]
cl2 <- top05[top05["cluster"] == 2,]
cl3 <- top05[top05["cluster"] == 3,]
cl4 <- top05[top05["cluster"] == 4,]
cl5 <- top05[top05["cluster"] == 5,]
cl6 <- top05[top05["cluster"] == 6,]
cl7 <- top05[top05["cluster"] == 7,]
cl8 <- top05[top05["cluster"] == 8,]
cl9 <- top05[top05["cluster"] == 9,]
cl10 <- top05[top05["cluster"] == 10,]
cl11 <- top05[top05["cluster"] == 11,]
cl12 <- top05[top05["cluster"] == 12,]
cl13 <- top05[top05["cluster"] == 13,]
cl14 <- top05[top05["cluster"] == 14,]
cl15 <- top05[top05["cluster"] == 15,]

#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
#install.packages("SeuratObject")

library("BiocManager")
library("multtest")
library("metap")
library("Seurat")






