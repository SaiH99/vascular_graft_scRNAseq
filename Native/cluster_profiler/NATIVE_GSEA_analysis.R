# Load the libraries
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ggupset)
library(Seurat)
library(dplyr)

# Reading the genes in the sequencing genome
#NATIVE <- readRDS("results/NATIVE_prenormalization.RDS")
#genes_in_data <- Features(NATIVE)

# Loading the annotation and filtering as per the genome
#file <- "cluster_profiler/c2.cp.v2023.2.Hs.symbols.gmt"
#pwl2 <- read.gmt(file)
#pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]

# Saving the final Gene Annotation List
#filename <- "GSEA_all.RDS"
#saveRDS(pwl2, paste0("cluster_profiler/", filename))

# Preparing the datasets for GSEA
NATIVE <- readRDS("results/NATIVE_post_annotation.RDS")
DefaultAssay(NATIVE) <- "RNA"
diff_markers <- FindAllMarkers(object = NATIVE, densify = TRUE)
bg_genes <- readRDS("cluster_profiler/GSEA_all.RDS")
output_dir <- "cluster_profiler/GSEA_rough"

# The clusters to run GSEA on
cell_types <- c("Endothelial Cells", "Macrophages")

# Cluster specific function to carry out GSEA
GSEA_analysis <- function(cell_type) {
  
  # Filtering the differential expression matrix
  diff_filtered <- diff_markers %>%
    dplyr::filter(!str_detect(gene, "ENSOAR"), cluster == cell_type, pct.1 > 0.50)
  
  # Correcting for duplicates
  diff_filtered$gene <- sub("\\..*$", "", diff_filtered$gene)
  
  # Preparation of Ranked List
  genelist <- diff_filtered$avg_log2FC
  names(genelist) <- diff_filtered$gene
  genelist <- sort(genelist, decreasing = TRUE)
  
  # Running GSEA (fgsea)
  res <- GSEA(genelist, TERM2GENE = bg_genes, pvalueCutoff = 0.25, eps = 1e-300, minGSSize = 20, nPermSimple = 100000)
  
  # Export the result
  output_file <- paste0(output_dir, "/GSEA_complete_NATIVE_", cell_type, ".csv")
  write.csv(res@result, file = output_file, row.names = F)
  
  # Visualization of the pathway analysis
  
}

for (type in cell_types) {
  GSEA_analysis(type)
}
