# Load the libraries ----
library(SingleCellExperiment)
library(scDblFinder)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(BiocParallel)
library(SoupX)
#devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(HGNChelper)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(viridis)
library(ggthemes)
library(colorblindr)
library(Nebulosa)
library(qpdf)
set.seed("123")

# Load data
seurat_data <- Read10X(data.dir = "data/NATIVE/filtered_feature_bc_matrix")
NATIVE <- CreateSeuratObject(counts = seurat_data, 
                           min.features = 800,
                           min.cells = 1,
                           project = "NATIVE")

# Doublet Detection
sce <- scDblFinder(GetAssayData(NATIVE, layer = 'counts'))
NATIVE$scDblFinder.score <- sce$scDblFinder.score
NATIVE$scDblFinder.class <- sce$scDblFinder.class
rm(sce, seurat_data)

# Load mitochondrial genes
load("./data/mitoGenes.rda")

# Add number of genes per UMI for each cell to metadata
NATIVE$log10GenesPerUMI <- log10(NATIVE$nFeature_RNA) / log10(NATIVE$nCount_RNA)

# Add mitoRatio column for each cell
NATIVE@meta.data["mitoRatio"] = PercentageFeatureSet(object = NATIVE, features = mt_genes)

# Set memory limits
options(future.globals.maxSize = 4000 * 1024^2)

# Subset based on mitochondrial percentage
NATIVE <- subset(NATIVE, subset = (mitoRatio < 10 & nCount_RNA > 1000 & scDblFinder.class == "singlet"))
saveRDS(NATIVE, file = "results/NATIVE_prenormalization.RDS")

# Normalization and Scaling
NATIVE <- SCTransform(NATIVE, vars.to.regress = c("mitoRatio"), verbose = FALSE)

# Running PCA
NATIVE <- RunPCA(NATIVE, npcs = 100)
# ElbowPlot(NATIVE, ndims = 100, reduction = "pca")

# k nearest neighbor and clustering
NATIVE <- FindNeighbors(NATIVE, dims = 1:25, reduction = "pca")
NATIVE <- FindClusters(NATIVE, algorithm = 4, resolution = 0.2, cluster.name = "cluster_0.2")

# UMAP and tsne
NATIVE <- RunTSNE(NATIVE, dims = 1:75, reduction = "pca", seed.use = "123")
NATIVE <- RunUMAP(NATIVE, dims = 1:75, reduction = "pca", reduction.name = "umap", seed.use = "123")

# Visualizing UMAP and TSNE
DimPlot(NATIVE, reduction = "umap", group.by = "cluster_0.2", label = TRUE, seed = "123")
# TSNEPlot(NATIVE, reduction = "tsne", group.by = "cluster_0.2", label = TRUE, seed = "123")

# Setting default clustering identity
Idents(NATIVE) <- "cluster_0.2"

# Normalization for Heatmaps
DefaultAssay(NATIVE) <- "RNA"
NATIVE <- NormalizeData(NATIVE, verbose = FALSE)
NATIVE <- ScaleData(NATIVE, features = rownames(NATIVE))

# Labeling Clusters
NATIVE <- RenameIdents(NATIVE, "2" = "1", "5" = "1")
NATIVE <- RenameIdents(NATIVE, "1" = "Endothelial Cells", "6" = "Macrophages", "4" = "Fibroblasts", "8" = "B Lymphocytes", "3" = "T/NK Cells", "7" = "SMCs")
levels(NATIVE) <- c("Fibroblasts", "Endothelial Cells", "SMCs", "Macrophages", "B Lymphocytes", "T/NK Cells")
NATIVE@meta.data$seurat_clusters <- Idents(NATIVE)

# Saving the final Seurat File
saveRDS(NATIVE, file = "results/NATIVE_post_annotation.RDS")

# Differential Expression Testing for Cluster markers
DefaultAssay(NATIVE) <- "RNA"
diff_markers <- FindAllMarkers(object = NATIVE, only.pos = TRUE, densify = TRUE)
diff_markers %>%
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05, !str_detect(gene, "ENSOAR"), avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15

# Export the marker genes
write.csv(top15, "results/Top_15_Marker_genes.csv")
