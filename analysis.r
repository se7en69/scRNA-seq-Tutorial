```r
#' Single-Cell RNA-Seq Analysis of B16 Melanoma Cells
#' 
#' A comprehensive analysis pipeline for processing and analyzing scRNA-seq data
#' from B16 melanoma cells using the Seurat package.

# Load required packages -------------------------------------------------------
required_packages <- c(
  "Seurat",          # Single-cell analysis
  "dplyr",           # Data manipulation
  "ggplot2",         # Data visualization
  "patchwork",       # Plot arrangement
  "SingleR",         # Cell type annotation
  "celldex",         # Reference datasets
  "RColorBrewer",    # Color palettes
  "Matrix"           # Sparse matrix support
)

# Install missing packages
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

# Set random seed for reproducibility
set.seed(1234)

# Configuration ---------------------------------------------------------------
config <- list(
  data_dir = "data/raw_feature_bc_matrix/",  # Path to 10X data
  project_name = "B16_Melanoma",
  min_cells = 3,        # Minimum cells expressing a gene
  min_features = 200,   # Minimum features per cell
  mt_pattern = "^mt-",  # Mitochondrial gene pattern (mouse)
  max_mt_percent = 15,  # Maximum mitochondrial percentage
  min_features = 200,   # Minimum features per cell
  max_features = 6000,  # Maximum features per cell
  nfeatures = 2000,     # Number of variable features
  npcs = 50,            # Number of PCs to compute
  resolution = 0.6,     # Clustering resolution
  output_dir = "results/"
)

# Create output directories if they don't exist
dir.create(config$output_dir, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(config$output_dir, "tables"), showWarnings = FALSE)

# Helper Functions ------------------------------------------------------------
save_plot <- function(plot, filename, width = 8, height = 6) {
  ggplot2::ggsave(
    file.path(config$output_dir, "figures", filename),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# Data Loading -----------------------------------------------------------------
cat("\nStep 1: Loading 10X Genomics data...\n")

# Read 10X data
raw_data <- Seurat::Read10X(data.dir = config$data_dir)

# Create Seurat object
seurat_obj <- Seurat::CreateSeuratObject(
  counts = raw_data,
  project = config$project_name,
  min.cells = config$min_cells,
  min.features = config$min_features
)

# Add metadata
seurat_obj$genotype <- "WT"  # Options: "WT" or "KO"
seurat_obj$sample <- "sample1"

# Quality Control -------------------------------------------------------------
cat("\nStep 2: Performing quality control...\n")

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj, 
  pattern = config$mt_pattern
)

# Visualize QC metrics
qc_plots <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.1,
  ncol = 3
)

save_plot(qc_plots, "qc_violin_plots.png")

# Filter cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > config$min_features &
    nFeature_RNA < config$max_features &
    percent.mt < config$max_mt_percent
)

# Normalization and Feature Selection ------------------------------------------
cat("\nStep 3: Normalizing data and selecting features...\n")

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = config$nfeatures
)

# Visualize variable features
var_plot <- VariableFeaturePlot(seurat_obj)
var_plot <- LabelPoints(
  plot = var_plot,
  points = head(VariableFeatures(seurat_obj), 20),
  repel = TRUE
)

save_plot(var_plot, "variable_features.png", width = 10, height = 6)

# Dimensionality Reduction and Clustering -------------------------------------
cat("\nStep 4: Performing dimensionality reduction and clustering...\n")

seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = config$npcs)

# Elbow plot to determine significant PCs
elbow_plot <- ElbowPlot(seurat_obj, ndims = config$npcs)
save_plot(elbow_plot, "elbow_plot.png")

# Cluster cells (using 15 PCs based on elbow plot)
pcs_used <- 15
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pcs_used)
seurat_obj <- FindClusters(seurat_obj, resolution = config$resolution)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pcs_used)

# Visualize clusters
umap_clusters <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP: Clustered Cells")

save_plot(umap_clusters, "umap_clusters.png")

# Cell Type Annotation --------------------------------------------------------
cat("\nStep 5: Annotating cell types...\n")

ref <- celldex::MouseRNAseqData()
annotations <- SingleR(
  test = GetAssayData(seurat_obj, slot = "data"),
  ref = ref,
  labels = ref$label.main
)

seurat_obj$celltype <- annotations$labels

# Visualize annotated cell types
umap_celltypes <- DimPlot(
  seurat_obj,
  group.by = "celltype",
  label = TRUE,
  repel = TRUE
) + ggtitle("UMAP: Annotated Cell Types")

save_plot(umap_celltypes, "umap_celltypes.png")

# Marker Gene Analysis ---------------------------------------------------------
cat("\nStep 6: Identifying marker genes...\n")

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save top markers
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(
  top_markers,
  file.path(config$output_dir, "tables", "cluster_markers.csv"),
  row.names = FALSE
)

# Differential Expression Analysis ---------------------------------------------
cat("\nPerforming differential expression analysis...\n")

de_genes <- FindMarkers(
  seurat_obj,
  ident.1 = 0,
  ident.2 = 1,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Add significance column
de_genes$gene <- rownames(de_genes)
de_genes$significant <- ifelse(
  de_genes$p_val_adj < 0.05 & abs(de_genes$avg_log2FC) > 0.5,
  "Yes",
  "No"
)

# Save DE results
write.csv(
  de_genes,
  file.path(config$output_dir, "tables", "cluster0_vs_cluster1_DEGs.csv"),
  row.names = FALSE
)

# Volcano plot
volcano_plot <- ggplot(de_genes, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Cluster 0 vs Cluster 1",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value"
  )

save_plot(volcano_plot, "volcano_plot.png")

# Save Seurat Object ----------------------------------------------------------
saveRDS(
  seurat_obj,
  file.path(config$output_dir, "seurat_object_processed.rds")
)

cat("\nAnalysis complete! Results saved to", config$output_dir, "\n")