# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  06-go-enrichment.R
# Purpose: GO enrichment analysis for all hdWGCNA modules and spatial projection
# Input:   results/objects/ST_LR_WGCNA_finished.rds
# Output:  results/tables/hdWGCNA_*_module_BP.csv
#          results/figures/go_*.jpg
#         results/figures/*_of_hdWGCNA.jpg
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")
source("scripts/utils.R")

ST_merged_hd <- readRDS(file.path(OBJECTS_DIR, "ST_LR_WGCNA_finished.rds"))
ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Define modules and their parameters
# ------------------------------------------------------------------------------

modules <- GetModules(ST_merged_hd)

module_config <- list(
    list(name = "green",     fill = "green",     text = "black", simplify = 0.7),
    list(name = "turquoise", fill = "turquoise",  text = "black", simplify = 0.7),
    list(name = "brown",    fill = "brown",     text = "black", simplify = 0.7),
    list(name = "blue",     fill = "blue",      text = "white", simplify = 0.7),
    list(name = "red",      fill = "red",       text = "black", simplify = 0.9),
    list(name = "yellow",   fill = "yellow",    text = "black", simplify = 0.7),
    list(name = "black",    fill = "black",     text = "white", simplify = 0.7)
)

# ------------------------------------------------------------------------------
# Run GO enrichment for each module
# ------------------------------------------------------------------------------

for (mod in module_config) {
    module_genes <- modules %>%
        filter(color == mod$name) %>%
        pull(gene_name)

    run_module_go_enrichment(
        module_genes  = module_genes,
        module_name   = mod$name,
        org_db        = org.Wm82.eg.db,
        bar_fill      = mod$fill,
        text_color    = mod$text,
        simplify_cutoff = mod$simplify,
        output_dir    = RESULTS_DIR
    )
}

# ------------------------------------------------------------------------------
# Spatial projection of module scores
# ------------------------------------------------------------------------------

# Reload ST_merged for score calculations
ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

spatial_modules <- c("turquoise", "blue", "brown", "black")

for (mod_name in spatial_modules) {
    mod_genes <- modules %>%
        filter(color == mod_name) %>%
        pull(gene_name)

    ST_merged <- AddModuleScore(
        ST_merged,
        features = list(mod_genes),
        name = paste0(mod_name, "_Score")
    )

    p <- SpatialFeaturePlot(
        ST_merged,
        features = paste0(mod_name, "_Score1"),
        crop = FALSE,
        pt.size.factor = 8,
        ncol = 3
    ) + coord_fixed()

    ggsave(file.path(FIGURES_DIR, paste0(mod_name, "_of_hdWGCNA.jpg")),
           plot = p, units = "mm", height = 200, width = 160,
           dpi = 1200, device = "jpeg")
}

# Peptidyl module score
ST_merged <- AddModuleScore(
    ST_merged,
    features = list(peptidyl_gene),
    name = "peptidyl_Score"
)

p_peptidyl <- SpatialFeaturePlot(
    ST_merged,
    features = "peptidyl_Score1",
    crop = FALSE,
    pt.size.factor = 8,
    ncol = 3
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "peptidyl_of_hdWGCNA.jpg"),
       plot = p_peptidyl, units = "mm", height = 200, width = 160,
       dpi = 1200, device = "jpeg")

# ------------------------------------------------------------------------------
# Plot KMEs
# ------------------------------------------------------------------------------

ST_merged_hd <- readRDS(file.path(OBJECTS_DIR, "ST_LR_WGCNA_finished.rds"))

p_kmes <- PlotKMEs(ST_merged_hd, ncol = 3)
print(p_kmes)

message("06-go-enrichment.R: Complete.")
