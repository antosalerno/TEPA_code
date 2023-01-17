### 3D UMAP plot
# inspired by @Dragonmasterx87

# Install plot_ly
install.packages('plotly')

# Load plot_ly
library(plotly)
library("Seurat")
library("SeuratObject")
library(SeuratDisk)
library("ggplot2")

seuset_immune <- LoadH5Seurat("TEPA_results/03_seusetImmuneModule.h5Seurat")

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
seuset_immune <- RunUMAP(seuset_immune,
                            dims = 1:60,
                         reduction = "pca",
                            n.components = 3L)

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(seuset_immune, reduction = "umap")

### Plot clustering on 3D umap ####

# Prepare a dataframe for cell plotting
plot.data <- FetchData(seuset_immune, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "scType"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
               color = ~scType, 
               colors = c("lightseagreen",
                                         "red",
                                         "turquoise4",
                                         "royalblue1",
                                         "lightcyan3",
                                         "peachpuff3",
                                         "orange2",
                                         "yellow3",
                                         "darkorchid1",
                                         "plum2",
                                         "darkmagenta"),
                                         type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)

fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig
fig_cube

### Plot expression of a particular gene ####

DefaultAssay(seuset_immune) <- "integrated"
# create a dataframe
plot.data <- FetchData(seuset_immune, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "macro_markers_mouse"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
plot.data$changed <- ifelse(test = plot.data$macro_markers_mouse > 0, yes = plot.data$macro_markers_mouse, no = 0)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.data$label <- paste(rownames(plot.data)," - ", plot.data$macro_markers_mouse, sep="")

# Plot your data
dev.off()
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('purple', 'cyan'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 5, width=2), 
        text=~label,
        hoverinfo="text"
)

### Plot expression of a gene module ####

DefaultAssay(seuset_immune) <- "integrated"
goi <- "macro_markers_mouse"
yourseuratobject <- seuset_immune
plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
Cutoff <- quantile(plotting.data[,goi], probs = .95)
plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
plot_ly(data = plotting.data,
        # name = goi,
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
        opacity = .5,
        colors = c('purple', 'cyan'), 
        type = "scatter3d", 
        mode = "markers",
        marker = list(size = 1), 
        text=~label,
        hoverinfo="text"
) %>%layout(title=goi)
