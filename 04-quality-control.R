# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  04-quality-control.R
# Purpose: Quality control visualization and correlation analysis
# Input:   ST_LR_merged.rds
# Output:  results/figures/*.jpg (QC plots)
#          results/tables/QC_summary_13samples.csv
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Mitochondrial percentage
# ------------------------------------------------------------------------------

ST_merged[["percent.mt"]] <- PercentageFeatureSet(
    ST_merged,
    pattern = "Glyma.08G107000"
)

# ------------------------------------------------------------------------------
# QC violin plot
# ------------------------------------------------------------------------------

p_vln <- VlnPlot(
    ST_merged,
    features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
    group.by = "chip_id",
    pt.size = 0
) + theme_classic(base_size = 16)

ggsave(file.path(FIGURES_DIR, "qc_vlnplot_ST_merged.jpg"),
       plot = p_vln, height = 100, width = 200,
       device = "jpeg", dpi = 400, units = "mm")

# ------------------------------------------------------------------------------
# Spatial QC feature plots
# ------------------------------------------------------------------------------

p_ncount <- SpatialFeaturePlot(
    ST_merged,
    features = "nCount_Spatial",
    pt.size.factor = 16,
    ncol = 3,
    crop = FALSE
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "nCount_Spatial_ST_merged.jpg"),
       plot = p_ncount, height = 200, width = 200,
       device = "jpeg", dpi = 1200, units = "mm")

p_nfeature <- SpatialFeaturePlot(
    ST_merged,
    features = "nFeature_Spatial",
    pt.size.factor = 16,
    ncol = 3,
    crop = FALSE
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "nFeature_Spatial_ST_merged.jpg"),
       plot = p_nfeature, height = 200, width = 200,
       device = "jpeg", dpi = 1200, units = "mm")

p_pmt <- SpatialFeaturePlot(
    ST_merged,
    features = "percent.mt",
    pt.size.factor = 16,
    ncol = 3,
    crop = FALSE
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "percent_mt_ST_merged.jpg"),
       plot = p_pmt, height = 200, width = 200,
       device = "jpeg", dpi = 1200, units = "mm")

# ------------------------------------------------------------------------------
# Feature scatter
# ------------------------------------------------------------------------------

clusters <- levels(ST_merged$seurat_clusters)
cluster_colors <- cluster_colors_24[1:length(clusters)]
names(cluster_colors) <- clusters

p_scatter <- FeatureScatter(
    ST_merged,
    feature1 = "nCount_Spatial",
    feature2 = "nFeature_Spatial",
    cols = cluster_colors
) +
    geom_smooth(method = "lm", color = "red") +
    theme_classic(base_size = 16)

ggsave(file.path(FIGURES_DIR, "feature_scatter_ST_merged.jpg"),
       plot = p_scatter, height = 100, width = 130,
       device = "jpeg", dpi = 300, units = "mm")

# ------------------------------------------------------------------------------
# UMAP by chip ID
# ------------------------------------------------------------------------------

p_chip <- DimPlot(ST_merged, group.by = "chip_id")

ggsave(file.path(FIGURES_DIR, "umap_by_chip_id.jpg"),
       plot = p_chip, height = 100, width = 130,
       device = "jpeg", dpi = 300, units = "mm")

# ------------------------------------------------------------------------------
# Correlation heatmap
# ------------------------------------------------------------------------------

ST_merged <- JoinLayers(ST_merged, assay = "Spatial")

pseudo_bulk <- AggregateExpression(
    ST_merged,
    group.by = "chip_id",
    assays = "Spatial",
    fun = "sum"
)

expr_mat <- as.matrix(pseudo_bulk$Spatial)
cor_mat <- cor(log1p(expr_mat), method = "pearson")

jpeg(file.path(FIGURES_DIR, "sample_correlation_heatmap.jpg"),
     width = 8, height = 8, units = "in", res = 300)
pheatmap(cor_mat)
dev.off()

# ------------------------------------------------------------------------------
# QC summary table
# ------------------------------------------------------------------------------

meta_df <- ST_merged@meta.data

qc_summary <- meta_df %>%
    group_by(orig.ident) %>%
    summarise(
        n_spots = n(),
        mean_nCount = mean(nCount_Spatial, na.rm = TRUE),
        mean_nFeature = mean(nFeature_Spatial, na.rm = TRUE)
    )

write.csv(qc_summary,
          file = file.path(TABLES_DIR, "QC_summary_13samples.csv"),
          row.names = FALSE)

message("04-quality-control.R: Complete.")
