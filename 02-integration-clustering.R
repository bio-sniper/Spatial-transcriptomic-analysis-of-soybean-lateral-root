<<<<<<< HEAD
# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  02-integration-clustering.R
# Purpose: Data integration, Harmony, UMAP, clustering, marker identification
# Input:   ST_LR_merged.rds (Seurat object)
# Output:  ST_LR_merged.rds (annotated Seurat object with clusters)
#          results/figures/*.jpg (UMAP, spatial, marker plots)
#          results/tables/markers.csv
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

message("Samples in dataset:")
print(table(ST_merged$orig.ident))

# ------------------------------------------------------------------------------
# Integration and clustering
# ------------------------------------------------------------------------------

DefaultAssay(ST_merged) <- "SCT"

features <- SelectIntegrationFeatures(
    object.list = SplitObject(ST_merged, split.by = "orig.ident"),
    nfeatures = 3000
)

ST_merged <- RunPCA(ST_merged, npcs = 50, features = features, verbose = FALSE)

ST_merged <- RunHarmony(ST_merged, group.by.vars = "chip_id")

ST_merged <- RunUMAP(
    ST_merged,
    reduction = "harmony",
    dims = 1:30,
    seed.use = 42
)

ST_merged <- FindNeighbors(ST_merged, reduction = "harmony", dims = 1:30)

ST_merged <- FindClusters(ST_merged, resolution = 0.2, random.seed = 42)

# ------------------------------------------------------------------------------
# UMAP visualization
# ------------------------------------------------------------------------------

clusters <- levels(ST_merged$seurat_clusters)
cluster_colors <- cluster_colors_24[1:length(clusters)]
names(cluster_colors) <- clusters

p_umap <- DimPlot(
    ST_merged,
    reduction = "umap",
    pt.size = 0.2,
    group.by = "seurat_clusters",
    cols = cluster_colors
) +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

ggsave(file.path(FIGURES_DIR, "umap_clusters.jpg"),
       plot = p_umap, units = "mm", height = 200, width = 200,
       dpi = 300, device = "jpeg")

# ------------------------------------------------------------------------------
# Spatial visualization
# ------------------------------------------------------------------------------

p_spatial <- SpatialDimPlot(
    ST_merged,
    label.size = 3,
    pt.size.factor = 10,
    alpha = 0.8,
    stroke = 0.2,
    crop = FALSE,
    cols = cluster_colors,
    ncol = 3
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "spatial_cluster.jpg"),
       plot = p_spatial, units = "mm", height = 200, width = 200,
       dpi = 1200, device = "jpeg")

# ------------------------------------------------------------------------------
# Rename cluster identities
# ------------------------------------------------------------------------------

new.cluster.ids <- c(
    "Cortex-Endodermis", "Epidermis", "Pericycle", "Xylem",
    "Lateral root", "Phloem", "Lateral root junction"
)

names(new.cluster.ids) <- levels(ST_merged)
ST_merged <- RenameIdents(ST_merged, new.cluster.ids)

saveRDS(ST_merged, file = file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Marker gene identification
# ------------------------------------------------------------------------------

ST_merged <- PrepSCTFindMarkers(ST_merged)

markers <- FindAllMarkers(
    ST_merged,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
)

write.csv(markers, file = file.path(TABLES_DIR, "marker_genes.csv"),
          row.names = FALSE)

top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# ------------------------------------------------------------------------------
# Dot plot of top markers
# ------------------------------------------------------------------------------

p_dot <- DotPlot(
    ST_merged,
    features = unique(top_markers$gene),
    group.by = "seurat_clusters",
    cols = c("blue", "red"),
    dot.scale = 6,
    idents = NULL
) +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 16, angle = 60,
                                   hjust = 1, face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 40, unit = "pt")
    )

ggsave(file.path(FIGURES_DIR, "marker_gene_dotplot.jpg"),
       plot = p_dot, units = "mm", height = 150, width = 400,
       dpi = 300, device = "jpeg")

message("02-integration-clustering.R: Complete.")
=======
# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  02-integration-clustering.R
# Purpose: Data integration, Harmony, UMAP, clustering, marker identification
# Input:   ST_LR_merged.rds (Seurat object)
# Output:  ST_LR_merged.rds (annotated Seurat object with clusters)
#          results/figures/*.jpg (UMAP, spatial, marker plots)
#          results/tables/markers.csv
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

message("Samples in dataset:")
print(table(ST_merged$orig.ident))

# ------------------------------------------------------------------------------
# Integration and clustering
# ------------------------------------------------------------------------------

DefaultAssay(ST_merged) <- "SCT"

features <- SelectIntegrationFeatures(
    object.list = SplitObject(ST_merged, split.by = "orig.ident"),
    nfeatures = 3000
)

ST_merged <- RunPCA(ST_merged, npcs = 50, features = features, verbose = FALSE)

ST_merged <- RunHarmony(ST_merged, group.by.vars = "chip_id")

ST_merged <- RunUMAP(
    ST_merged,
    reduction = "harmony",
    dims = 1:30,
    seed.use = 42
)

ST_merged <- FindNeighbors(ST_merged, reduction = "harmony", dims = 1:30)

ST_merged <- FindClusters(ST_merged, resolution = 0.2, random.seed = 42)

# ------------------------------------------------------------------------------
# UMAP visualization
# ------------------------------------------------------------------------------

clusters <- levels(ST_merged$seurat_clusters)
cluster_colors <- cluster_colors_24[1:length(clusters)]
names(cluster_colors) <- clusters

p_umap <- DimPlot(
    ST_merged,
    reduction = "umap",
    pt.size = 0.2,
    group.by = "seurat_clusters",
    cols = cluster_colors
) +
    theme_classic() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    )

ggsave(file.path(FIGURES_DIR, "umap_clusters.jpg"),
       plot = p_umap, units = "mm", height = 200, width = 200,
       dpi = 300, device = "jpeg")

# ------------------------------------------------------------------------------
# Spatial visualization
# ------------------------------------------------------------------------------

p_spatial <- SpatialDimPlot(
    ST_merged,
    label.size = 3,
    pt.size.factor = 10,
    alpha = 0.8,
    stroke = 0.2,
    crop = FALSE,
    cols = cluster_colors,
    ncol = 3
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "spatial_cluster.jpg"),
       plot = p_spatial, units = "mm", height = 200, width = 200,
       dpi = 1200, device = "jpeg")

# ------------------------------------------------------------------------------
# Rename cluster identities
# ------------------------------------------------------------------------------

new.cluster.ids <- c(
    "Cortex-Endodermis", "Epidermis", "Pericycle", "Xylem",
    "Lateral root", "Phloem", "Lateral root junction"
)

names(new.cluster.ids) <- levels(ST_merged)
ST_merged <- RenameIdents(ST_merged, new.cluster.ids)

saveRDS(ST_merged, file = file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Marker gene identification
# ------------------------------------------------------------------------------

ST_merged <- PrepSCTFindMarkers(ST_merged)

markers <- FindAllMarkers(
    ST_merged,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
)

write.csv(markers, file = file.path(TABLES_DIR, "marker_genes.csv"),
          row.names = FALSE)

top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)

# ------------------------------------------------------------------------------
# Dot plot of top markers
# ------------------------------------------------------------------------------

p_dot <- DotPlot(
    ST_merged,
    features = unique(top_markers$gene),
    group.by = "seurat_clusters",
    cols = c("blue", "red"),
    dot.scale = 6,
    idents = NULL
) +
    theme(
        text = element_text(size = 12),
        axis.text.x = element_text(size = 16, angle = 60,
                                   hjust = 1, face = "italic"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 40, unit = "pt")
    )

ggsave(file.path(FIGURES_DIR, "marker_gene_dotplot.jpg"),
       plot = p_dot, units = "mm", height = 150, width = 400,
       dpi = 300, device = "jpeg")

message("02-integration-clustering.R: Complete.")
>>>>>>> 51aa9eb2afbce5797fc300fa3444dd5b09dd6cef
