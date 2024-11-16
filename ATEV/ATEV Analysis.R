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
seurat_data <- Read10X(data.dir = "data/ATEV/filtered_feature_bc_matrix")
ATEV <- CreateSeuratObject(counts = seurat_data, 
                           min.features = 600,
                           min.cells = 1,
                           project = "ATEV")

# Doublet Detection
sce <- scDblFinder(GetAssayData(ATEV, layer = 'counts'))
ATEV$scDblFinder.score <- sce$scDblFinder.score
ATEV$scDblFinder.class <- sce$scDblFinder.class
rm(sce, seurat_data)

# Load mitochondrial genes
load("./data/mitoGenes.rda")

# Add number of genes per UMI for each cell to metadata
ATEV$log10GenesPerUMI <- log10(ATEV$nFeature_RNA) / log10(ATEV$nCount_RNA)

# Add mitoRatio column for each cell
ATEV@meta.data["mitoRatio"] = PercentageFeatureSet(object = ATEV, features = mt_genes)

# Set memory limits
options(future.globals.maxSize = 4000 * 1024^2)

# Subset based on mitochondrial percentage
ATEV <- subset(ATEV, subset = (mitoRatio < 15 & nCount_RNA > 1200 & scDblFinder.class == "singlet"))
saveRDS(ATEV, file = "results/ATEV_prenormalization.RDS")

# Normalization and Scaling
ATEV <- SCTransform(ATEV, verbose = FALSE)

# Running PCA
ATEV <- RunPCA(ATEV, npcs = 100)
# ElbowPlot(ATEV, ndims = 100, reduction = "pca")

# k nearest neighbor and clustering
ATEV <- FindNeighbors(ATEV, dims = 1:75, reduction = "pca", k.param = 10)
ATEV <- FindClusters(ATEV, algorithm = 4, resolution = 1, cluster.name = "cluster_1")

# UMAP and tsne
ATEV <- RunTSNE(ATEV, dims = 1:75, reduction = "pca", seed.use = "123")
ATEV <- RunUMAP(ATEV, dims = 1:75, reduction = "pca", reduction.name = "umap", seed.use = "123")

# Visualizing UMAP and TSNE
DimPlot(ATEV, reduction = "umap", group.by = "cluster_1", label = TRUE, seed = "123")
# TSNEPlot(ATEV, reduction = "tsne", group.by = "cluster_1", label = TRUE, seed = "123")

# Setting default clustering identity
Idents(ATEV) <- "cluster_1"

# Normalization for Heatmaps
DefaultAssay(ATEV) <- "RNA"
ATEV <- NormalizeData(ATEV, verbose = FALSE)
ATEV <- ScaleData(ATEV, features = rownames(ATEV))

# Combining redundant clusters
ATEV <- RenameIdents(ATEV, "2" = "1", "6" = "1", "16" = "1", "17" = "1")
ATEV <- RenameIdents(ATEV, "5" = "3", "13" = "3", "8" = "7", "10" = "7", "12" = "7", "14" = "7", "15" = "7")

# Sub clustering
ATEV <- FindSubCluster(
  ATEV,
  cluster = 3,
  graph.name = "SCT_snn",
  subcluster.name = "sub.cluster_3",
  resolution = 0.5,
  algorithm = 4
)
Idents(ATEV) <- "sub.cluster_3"
ATEV <- RenameIdents(ATEV, "3_4" = "3_2", "3_3" = "3_1", "3_5" = "3_1", "3_6" = "3_1")

# Labeling Clusters
ATEV <- RenameIdents(ATEV, "7" = "Endothelial Cells", "1" = "Macrophages", "11" = "Fibroblasts", "3_2" = "Fibroblasts", "9" = "Neutrophils", "4" = "T/NK Cells", "3_1" = "SMCs")
levels(ATEV) <- c("Fibroblasts", "Endothelial Cells", "SMCs", "Macrophages", "Neutrophils", "T/NK Cells")
ATEV@meta.data$seurat_clusters <- Idents(ATEV)
saveRDS(ATEV, file = "results/ATEV_post_annotation.RDS")

# Differential Expression Testing for Cluster markers
DefaultAssay(ATEV) <- "RNA"
diff_markers <- FindAllMarkers(object = ATEV, densify = TRUE, only.pos = T)
diff_markers %>%
  group_by(cluster) %>% 
  dplyr::filter(p_val_adj < 0.05, !str_detect(gene, "ENSOAR"), avg_log2FC > 1) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15

# Export the marker genes
write.csv(top15, "results/Top_15_Marker_genes.csv")

