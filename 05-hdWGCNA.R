# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  05-hdWGCNA.R
# Purpose: hdWGCNA analysis integrating bulk RNA-seq module genes
# Input:   ST_LR_WGCNA_finished.rds, bulk_data (in workspace)
# Output:  results/figures/*.jpg, results/objects/ST_LR_WGCNA_finished.rds
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")

# ------------------------------------------------------------------------------
# Load pre-processed hdWGCNA object
# ------------------------------------------------------------------------------

hd_file <- file.path(OBJECTS_DIR, "ST_LR_WGCNA_finished.rds")
if (file.exists(hd_file)) {
    ST_merged_hd <- readRDS(hd_file)
} else {
    message("ST_LR_WGCNA_finished.rds not found, using ST_merged as base")
    ST_merged_hd <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))
}

# ------------------------------------------------------------------------------
# Import bulk WGCNA module genes
# ------------------------------------------------------------------------------

# NOTE: bulk_data should be loaded in the environment. If running fresh,
#       load from: results/objects/bulk_data.rds

cell_wall <- bulk_data[bulk_data$BULK_module == "turquoise", ]
cell_wall_gene <- cell_wall$Gene_ID

ribosome <- bulk_data[bulk_data$BULK_module == "brown", ]
ribosome_gene <- ribosome$Gene_ID

LR <- bulk_data[bulk_data$BULK_module == "pink", ]
LR_gene <- LR$Gene_ID

peptidyl <- bulk_data[bulk_data$BULK_module == "red", ]
peptidyl_gene <- peptidyl$Gene_ID

# ------------------------------------------------------------------------------
# Prepare for hdWGCNA
# ------------------------------------------------------------------------------

DimPlot(ST_merged_hd)

ST_merged_hd <- SeuratObject::UpdateSeuratObject(ST_merged_hd)
DefaultAssay(ST_merged_hd) <- "SCT"

my_candidate <- unique(c(cell_wall_gene, ribosome_gene, LR_gene, peptidyl_gene))
default_features <- VariableFeatures(ST_merged_hd, nfeatures = 3000)
existing_genes <- rownames(ST_merged_hd)
valid_candidate <- my_candidate[my_candidate %in% existing_genes]

message("Bulk module genes: ", length(my_candidate))
message("Matched spatial genes: ", length(valid_candidate))

combined_features <- unique(c(default_features, valid_candidate))
combined_features <- combined_features[combined_features %in% existing_genes]

ST_merged_hd <- SetupForWGCNA(
    ST_merged_hd,
    gene_select = "custom",
    features = combined_features,
    wgcna_name = "soybean_LR"
)

# ------------------------------------------------------------------------------
# Metacell construction
# ------------------------------------------------------------------------------

ST_merged_hd <- MetacellsByGroups(
    seurat_obj = ST_merged_hd,
    group.by = c("seurat_clusters", "orig.ident"),
    reduction = "harmony",
    k = 25,
    max_shared = 10,
    ident.group = "seurat_clusters",
    min_cells = 50,
    wgcna_name = "soybean_LR"
)

ST_merged_hd <- NormalizeMetacells(ST_merged_hd)

# Align ST with metacell
m_obj <- GetMetacellObject(ST_merged_hd)
m_obj$all_cells <- "all_groups"
DefaultAssay(m_obj) <- "SCT"
ST_merged_hd <- SetMetacellObject(ST_merged_hd, m_obj)

ST_merged_hd <- SetDatExpr(
    ST_merged_hd,
    group_name = "all_groups",
    group.by = "all_cells",
    assay = "SCT",
    layer = "data"
)

message("Genes in expression matrix: ",
        length(GetWGCNAGenes(ST_merged_hd)))

# ------------------------------------------------------------------------------
# Soft power selection
# ------------------------------------------------------------------------------

ST_merged_hd <- TestSoftPowers(ST_merged_hd, networkType = "signed")

plot_list <- PlotSoftPowers(ST_merged_hd)
wrap_plots(plot_list)

power_table <- GetPowerTable(ST_merged_hd)
head(power_table)

# ------------------------------------------------------------------------------
# Network construction
# ------------------------------------------------------------------------------

ST_merged_hd <- ConstructNetwork(
    ST_merged_hd,
    soft_power = 10,
    setDatExpr = FALSE,
    tom_name = "soybean_LR_network",
    tom_outdir = "TOM",
    overwrite_tom = TRUE,
    detectCutHeight = 0.995,
    minModuleSize = 30
)

tiff(file.path(FIGURES_DIR, "hdWGCNA_dendrogram.tif"),
     height = 4, width = 10, units = "in", res = 500, compression = "lzw")

PlotDendrogram(ST_merged_hd, main = "soybean_LR_network")
dev.off()

# ------------------------------------------------------------------------------
# Module eigengenes and connectivity
# ------------------------------------------------------------------------------

ST_merged_hd <- ScaleData(ST_merged_hd,
                          features = GetWGCNAGenes(ST_merged_hd))

ST_merged_hd <- ModuleEigengenes(
    ST_merged_hd,
    assay = "SCT",
    group.by.vars = NULL,
    verbose = TRUE
)

ST_merged_hd$all_cells <- "all_groups"

m_obj <- GetMetacellObject(ST_merged_hd)
m_obj$all_cells <- "all_groups"
ST_merged_hd <- SetMetacellObject(ST_merged_hd, m_obj)

ST_merged_hd <- ModuleConnectivity(
    ST_merged_hd,
    group.by = "all_cells",
    group_name = "all_groups"
)

saveRDS(ST_merged_hd,
        file = file.path(OBJECTS_DIR, "ST_LR_WGCNA_finished.rds"))

# ------------------------------------------------------------------------------
# Check candidate gene locations
# ------------------------------------------------------------------------------

modules <- GetModules(ST_merged_hd)
module_counts <- table(modules$color)
print(module_counts)

target_genes <- c(
    "Glyma.01G208200", "Glyma.03G185900", "Glyma.11G034000",
    "Glyma.13G325900", "Glyma.14G061500", "Glyma.06G250100"
)

target_modules <- modules %>%
    filter(gene_name %in% peptidyl_gene)
print(target_modules)

# ------------------------------------------------------------------------------
# Module expression scores and feature plots
# ------------------------------------------------------------------------------

ST_merged_hd <- ModuleExprScore(ST_merged_hd, n_genes = 25, method = "UCell")

plot_list_hmes <- ModuleFeaturePlot(ST_merged_hd, features = "hMEs", order = TRUE)
wrap_plots(plot_list_hmes, ncol = 4)

plot_list_scores <- ModuleFeaturePlot(
    ST_merged_hd,
    features = "scores",
    order = "shuffle",
    ucell = TRUE,
    point_size = 1
)

ggsave(file.path(FIGURES_DIR, "hdWGCNA_module_umap.jpg"),
       plot = last_plot(), height = 5, width = 10,
       dpi = 500, device = "jpeg")

# ------------------------------------------------------------------------------
# Bulk epigenetic score feature plot
# ------------------------------------------------------------------------------

ST_merged_hd <- AddModuleScore(
    ST_merged_hd,
    features = list(peptidyl_gene),
    name = "Bulk_Epigenetic_Score"
)

p_peptidyl <- FeaturePlot(
    ST_merged_hd,
    features = "Bulk_Epigenetic_Score1",
    order = TRUE,
    cols = c("lightgrey", "#CC79A7"),
    pt.size = 4
)

ggsave(file.path(FIGURES_DIR, "hdWGCNA_peptidyl_score.jpg"),
       plot = p_peptidyl, width = 15, height = 12,
       dpi = 500, device = "jpeg")

# ------------------------------------------------------------------------------
# Module radar plot
# ------------------------------------------------------------------------------

ST_merged_hd$celltype <- Idents(ST_merged_hd)

p_radar <- ModuleRadarPlot(
    ST_merged_hd,
    group.by = "celltype",
    axis.label.size = 5,
    grid.label.size = 4
)

ggsave(file.path(FIGURES_DIR, "hdWGCNA_radar_plot.jpg"),
       plot = p_radar, width = 19, height = 19,
       dpi = 500, device = "jpeg")

# ------------------------------------------------------------------------------
# Module UMAP embedding
# ------------------------------------------------------------------------------

ST_merged_hd <- RunModuleUMAP(
    ST_merged_hd,
    n_hubs = 10,
    n_neighbors = 15,
    min_dist = 0.6,
    spread = 1
)

umap_df <- GetModuleUMAP(ST_merged_hd)

p_module_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
    geom_point(
        color = umap_df$color,
        size = umap_df$kME * 2
    ) +
    umap_theme()

pdf(file.path(FIGURES_DIR, "hdWGCNA_network.pdf"),
    height = 10, width = 10)

ModuleUMAPPlot(
    ST_merged_hd,
    edge.alpha = 0.5,
    sample_edges = TRUE,
    edge_prop = 0.075,
    label_hubs = 1,
    keep_grey_edges = FALSE
)

dev.off()

# ------------------------------------------------------------------------------
# Hub genes export
# ------------------------------------------------------------------------------

modules_df <- GetModules(ST_merged_hd)

top25_hubs_brown <- GetHubGenes(ST_merged_hd, mods = "brown", n_hubs = 25)
write.csv(top25_hubs_brown,
          file = file.path(TABLES_DIR, "Top25_Hub_Genes_brown.csv"),
          row.names = FALSE)

top25_hubs_turquoise <- GetHubGenes(ST_merged_hd, mods = "turquoise", n_hubs = 25)
write.csv(top25_hubs_turquoise,
          file = file.path(TABLES_DIR, "Top25_Hub_Genes_turquoise.csv"),
          row.names = FALSE)

top25_hubs_blue <- GetHubGenes(ST_merged_hd, mods = "blue", n_hubs = 25)
write.csv(top25_hubs_blue,
          file = file.path(TABLES_DIR, "Top25_Hub_Genes_blue.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Module network plots
# ------------------------------------------------------------------------------

ModuleNetworkPlot(ST_merged_hd,
                  outdir = file.path(FIGURES_DIR, "ModuleNetworks"))

message("05-hdWGCNA.R: Complete.")
