# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  07-lr-subcluster-monocle.R
# Purpose: Subcluster lateral root cells and perform monocle3 trajectory analysis
# Input:   results/objects/ST_LR_merged.rds
#          LR_cluster.rds, LR_LRP.rds, LR_overlying_cortex.rds (external)
# Output:  results/figures/*.jpg (trajectory plots)
#          results/tables/*.csv (pseudotime DEGs, gene modules)
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")
source("scripts/utils.R")

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Subset lateral root clusters from pre-computed RDS files
# ------------------------------------------------------------------------------

LR_cluster       <- readRDS(file.path(OBJECTS_DIR, "LR_cluster.rds"))
LR_cluster_LRP   <- readRDS(file.path(OBJECTS_DIR, "LR_LRP.rds"))
LR_cluster_cortex <- readRDS(file.path(OBJECTS_DIR, "LR_overlying cortex.rds"))

# ==============================================================================
# LRP trajectory (lateral root primordium)
# ==============================================================================

cd_LRP <- as.cell_data_set(LR_cluster_LRP)
colData(cd_LRP)$seurat_clusters <- LR_cluster_LRP@meta.data$seurat_clusters

cd_LRP <- preprocess_cds(cd_LRP, num_dim = 10)
cd_LRP <- reduce_dimension(cd_LRP, reduction_method = "UMAP")
reducedDims(cd_LRP)$UMAP <- Embeddings(LR_cluster_LRP, reduction = "umap")

cd_LRP <- cluster_cells(cd_LRP, cluster_method = "leiden",
                        k = 40, resolution = 0.7)

cd_LRP <- learn_graph(cd_LRP,
                      learn_graph_control = list(
                          minimal_branch_len = 60,
                          euclidean_distance_ratio = 10
                      ))

root_cells_2 <- colnames(cd_LRP)[
    colData(cd_LRP)$seurat_clusters == "0"
]
cd_LRP <- order_cells(cd_LRP, root_cells = root_cells_2)

# ------------------------------------------------------------------------------
# Plot LRP trajectory
# ------------------------------------------------------------------------------

plot_cells(cd_LRP, color_cells_by = "seurat_clusters",
           label_groups_by_cluster = TRUE)

p_lrp_pseudotime <- plot_cells(
    cd_LRP,
    color_cells_by = "pseudotime",
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_groups_by_cluster = FALSE,
    cell_size = 1,
    label_roots = FALSE
) +
    scale_color_gradientn(
        colors = viridis::viridis(100, option = "viridis"),
        name = "Pseudotime"
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

ggsave(file.path(FIGURES_DIR, "LR_LRP_pseudotime_monocle3.jpg"),
       plot = p_lrp_pseudotime, height = 100, width = 130,
       units = "mm", device = "jpeg", dpi = 300)

# ------------------------------------------------------------------------------
# Pseudotime DEGs for LRP
# ------------------------------------------------------------------------------

deg_pseudotime_LRP <- graph_test(cd_LRP,
                                 neighbor_graph = "principal_graph",
                                 cores = 4)

write.csv(deg_pseudotime_LRP,
          file = file.path(TABLES_DIR, "pseudotime_DEG_LR_LRP.csv"))

# Gene modules for LRP
lrp_modules <- run_pseudotime_module_analysis(
    cds            = cd_LRP,
    deg_results    = deg_pseudotime_LRP,
    output_prefix  = file.path(TABLES_DIR, "LR_LRP")
)

# GO enrichment for LRP modules
run_module_go_loop(
    module_genes_list = split(lrp_modules$gene_modules$id,
                              lrp_modules$gene_modules$module),
    org_db        = org.Wm82.eg.db,
    num_modules   = 21,
    output_dir    = file.path(RESULTS_DIR, "GO_results_LRP"),
    plot_dir      = file.path(FIGURES_DIR, "GO_plots_LRP")
)

# ------------------------------------------------------------------------------
# Plot LRP gene expression modules
# ------------------------------------------------------------------------------

gene_module_df_LRP <- lrp_modules$gene_modules[, c("id", "module")]

p_lrp_modules <- tryCatch({
    plot_cells(
        cd_LRP,
        genes = gene_module_df_LRP,
        label_roots = FALSE,
        label_branch_points = FALSE,
        label_leaves = FALSE,
        label_groups_by_cluster = FALSE,
        label_principal_points = FALSE,
        labels_per_group = FALSE,
        graph_label_size = 16
    ) + theme(text = element_text(size = 16))
}, error = function(e) {
    message("LRP expression module plot skipped: ", e$message)
    NULL
})

if (!is.null(p_lrp_modules)) {
    ggsave(file.path(FIGURES_DIR, "LR_LRP_expression_modules.jpg"),
           plot = p_lrp_modules, width = 200, height = 200,
           units = "mm", device = "jpeg", dpi = 600)
}

# Save LRP heatmap
ggsave(file.path(FIGURES_DIR, "LR_LRP_module_heatmap.jpg"),
       plot = lrp_modules$heatmap, width = 100, height = 100,
       units = "mm", device = "jpeg", dpi = 400)

# ------------------------------------------------------------------------------
# Plot targeted genes in LRP trajectory
# ------------------------------------------------------------------------------

fData(cd_LRP)$gene_short_name <- rownames(cd_LRP)

targeted_genes_in_cd_LRP <- c(
    "Glyma.06G028000",  # GmWAT1
    "Glyma.03G063600",  # GmAUX1
    "Glyma.15G068500",  # GmRPL26A
    "Glyma.10G210600"   # GmARF16
)

p_lrp_targets <- tryCatch({
    plot_cells(
        cd_LRP,
        genes = targeted_genes_in_cd_LRP,
        cell_size = 2,
        label_roots = FALSE,
        label_branch_points = FALSE,
        label_leaves = FALSE,
        label_groups_by_cluster = FALSE,
        label_principal_points = FALSE,
        labels_per_group = FALSE,
        graph_label_size = 16
    ) + theme(text = element_text(size = 16))
}, error = function(e) {
    message("LRP targeted genes plot skipped: ", e$message)
    NULL
})

if (!is.null(p_lrp_targets)) {
    ggsave(file.path(FIGURES_DIR, "LR_LRP_targeted_genes.jpg"),
           plot = p_lrp_targets, width = 200, height = 200,
           units = "mm", device = "jpeg", dpi = 400)
}


# ==============================================================================
# Cortex trajectory (overlying cortex)
# ==============================================================================

cds_cortex <- as.cell_data_set(LR_cluster_cortex)
colData(cds_cortex)$seurat_clusters <-
    LR_cluster_cortex@meta.data$seurat_clusters

cds_cortex <- preprocess_cds(cds_cortex, num_dim = 10)
cds_cortex <- reduce_dimension(cds_cortex, reduction_method = "UMAP")
reducedDims(cds_cortex)$UMAP <- Embeddings(LR_cluster_cortex,
                                           reduction = "umap")

cds_cortex <- cluster_cells(cds_cortex, cluster_method = "leiden",
                            k = 40, resolution = 0.7)

cds_cortex <- learn_graph(cds_cortex,
                          learn_graph_control = list(
                              minimal_branch_len = 30,
                              euclidean_distance_ratio = 10
                          ))

root_cells_cortex <- colnames(cds_cortex)[
    colData(cds_cortex)$seurat_clusters == "4"
]
cds_cortex <- order_cells(cds_cortex, root_cells = root_cells_cortex)

# ------------------------------------------------------------------------------
# Plot cortex trajectory
# ------------------------------------------------------------------------------

plot_cells(cds_cortex, color_cells_by = "seurat_clusters",
           label_groups_by_cluster = TRUE)

p_cortex_pseudotime <- plot_cells(
    cds_cortex,
    color_cells_by = "pseudotime",
    label_branch_points = FALSE,
    label_leaves = FALSE,
    label_groups_by_cluster = FALSE,
    cell_size = 1,
    label_roots = FALSE
) +
    scale_color_gradientn(
        colors = viridis::viridis(100, option = "viridis"),
        name = "Pseudotime"
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )

ggsave(file.path(FIGURES_DIR, "LR_cortex_pseudotime_monocle3.jpg"),
       plot = p_cortex_pseudotime, height = 100, width = 130,
       units = "mm", device = "jpeg", dpi = 300)

# ------------------------------------------------------------------------------
# Pseudotime DEGs for cortex
# ------------------------------------------------------------------------------

deg_pseudotime_cortex <- graph_test(cds_cortex,
                                    neighbor_graph = "principal_graph",
                                    cores = 4)

write.csv(deg_pseudotime_cortex,
          file = file.path(TABLES_DIR,
                           "pseudotime_DEG_LR_overlying_cortex.csv"))

# Gene modules for cortex
cortex_modules <- run_pseudotime_module_analysis(
    cds            = cds_cortex,
    deg_results    = deg_pseudotime_cortex,
    output_prefix  = file.path(TABLES_DIR, "LR_cortex")
)

# GO enrichment for cortex modules
run_module_go_loop(
    module_genes_list = split(cortex_modules$gene_modules$id,
                              cortex_modules$gene_modules$module),
    org_db        = org.Wm82.eg.db,
    num_modules   = 20,
    qvalue_cutoff = 0.9,
    output_dir    = file.path(RESULTS_DIR, "GO_results_cortex"),
    plot_dir      = file.path(FIGURES_DIR, "GO_plots_cortex")
)

# ------------------------------------------------------------------------------
# Plot cortex gene expression modules
# ------------------------------------------------------------------------------

rowData(cds_cortex)$gene_short_name <- rownames(cds_cortex)

gene_module_df_cortex <- cortex_modules$gene_modules[, c("id", "module")]

p_cortex_modules <- tryCatch({
    plot_cells(
        cds_cortex,
        genes = gene_module_df_cortex,
        label_roots = FALSE,
        label_branch_points = FALSE,
        label_leaves = FALSE,
        label_groups_by_cluster = FALSE,
        label_principal_points = FALSE,
        labels_per_group = FALSE,
        graph_label_size = 16
    ) + theme(text = element_text(size = 16))
}, error = function(e) {
    message("Cortex expression module plot skipped: ", e$message)
    NULL
})

if (!is.null(p_cortex_modules)) {
    ggsave(file.path(FIGURES_DIR, "LR_cortex_expression_modules.jpg"),
           plot = p_cortex_modules, width = 200, height = 200,
           units = "mm", device = "jpeg", dpi = 600)
}

# Save cortex heatmap
ggsave(file.path(FIGURES_DIR, "LR_cortex_module_heatmap.jpg"),
       plot = cortex_modules$heatmap, width = 100, height = 100,
       units = "mm", device = "jpeg", dpi = 400)

message("07-lr-subcluster-monocle.R: Complete.")
