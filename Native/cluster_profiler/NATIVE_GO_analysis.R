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
#file <- "cluster_profiler/c5.go.v2023.2.Hs.symbols.gmt"
#pwl2 <- read.gmt(file)
#pwl2 <- pwl2[pwl2$gene %in% genes_in_data,]

# Saving the final Gene Annotation List
#filename <- "GO.RDS"
#saveRDS(pwl2, paste0("cluster_profiler/", filename))

# Preparing the datasets for GO
NATIVE <- readRDS("results/NATIVE_post_annotation.RDS")
DefaultAssay(NATIVE) <- "RNA"
diff_markers <- FindAllMarkers(object = NATIVE, densify = TRUE)
bg_genes <- readRDS("cluster_profiler/GO.RDS")
output_dir <- "cluster_profiler/GO_rough"

# The clusters to run GO on
cell_types <- c("Endothelial Cells", "Macrophages")

# Cluster specific function to carry out GSEA
GO_analysis <- function(cell_type) {
  
  # Filtering the differential expression matrix
  diff_filtered <- diff_markers %>%
    dplyr::filter(!str_detect(gene, "ENSOAR"), cluster == cell_type, pct.1 > 0.50)
  
  # Correcting for duplicates
  diff_filtered$gene <- sub("\\..*$", "", diff_filtered$gene)
  
  # Marking the significantly UP and DOWN regulated genes
  diff_filtered <- diff_filtered %>% mutate(diffexpressed = case_when(
    avg_log2FC > 0 & p_val_adj < 0.05 ~ "UP",
    avg_log2FC < 0 & p_val_adj < 0.05 ~ "DOWN",
    p_val_adj > 0.05 ~ "NO"
  ))
  # Removing the non-significant genes
  diff_filtered <- diff_filtered[diff_filtered$diffexpressed != "NO",]
  
  # Splitting the dataset into UP and DOWN regulated genes
  diff_split <- split(diff_filtered, diff_filtered$diffexpressed)
  
  # Running enricher (GO Analysis)
  res <- lapply(names(diff_split), function(x) enricher(gene = diff_split[[x]]$gene, 
                                                        TERM2GENE = bg_genes, qvalueCutoff = 0.05))
  names(res) <- names(diff_split)
  
  # Splitting the result
  DOWN <- res$DOWN@result
  UP <- res$UP@result
  
  # Exporting the result
  write.csv(UP, file = paste0(output_dir, "/GO_UP_NATIVE_", cell_type, ".csv"), row.names = F)
  write.csv(DOWN, file = paste0(output_dir, "/GO_DOWN_NATIVE_", cell_type, ".csv"), row.names = F)
  
  # Visualization of the pathway analysis
  
}

for (type in cell_types) {
  GO_analysis(type)
}
