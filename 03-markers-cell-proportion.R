# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  03-markers-cell-proportion.R
# Purpose: Marker gene volcano plot (circular) and cell-type proportion analysis
# Input:   ST_LR_merged.rds
# Output:  results/figures/marker_gene_circle.tif
#          results/figures/cell_proportion.jpg
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Marker gene mapping and volcano plot
# ------------------------------------------------------------------------------

markers <- read.csv(file.path(TABLES_DIR, "marker_genes.csv"))

# Map Glyma IDs to common names
markers$gene <- ifelse(
    markers$gene %in% names(glyma_to_common),
    glyma_to_common[markers$gene],
    markers$gene
)

mygene <- c(
    "Glyma.17G139700", "Glyma.18G191800", "GmBGLU17", "GmGNS1",
    "Glyma.15G119100", "GmEXT3", "Glyma.17G030200", "GmCYP81D3",
    "GmPOT5", "Glyma.02G083702", "Glyma.07G139400", "GmRNS1",
    "Glyma.06G013200", "Glyma.18G190000", "Glyma.08G271600",
    "Glyma.07G096700", "GmCTL2", "GmXCP2", "GmSAHH1", "GmGER3",
    "GmPER54", "GmPER12", "Glyma.06G249700", "GmHDT2",
    "Glyma.13G236100", "GmMT3", "Glyma.05G216000", "Glyma.18G014000",
    "Glyma.13G094200", "GmAZI5", "GmGASA4", "Glyma.13G362800",
    "Glyma.17G008500", "GmFLA9"
)

p_volcano <- jjVolcano(
    diffData = markers,
    myMarkers = mygene,
    tile.col = c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00",
                 "#6a3d9a", "#b15928", "#a6cee3"),
    celltypeSize = 3,
    base_size = 8,
    col.type = "adjustP",
    fontface = "italic",
    polar = TRUE
) + ylim(-8, 10)

p_volcano$layers <- p_volcano$layers[
    !names(p_volcano$layers) %in% "geom_textpath"
]

ggsave(file.path(FIGURES_DIR, "marker_gene_circle.tif"),
       plot = p_volcano, width = 12, height = 12,
       dpi = 400, bg = "white", device = "tiff")

# ------------------------------------------------------------------------------
# Cell-type proportion by developmental stage
# ------------------------------------------------------------------------------

ST_merged$celltype <- Idents(ST_merged)
meta_df <- ST_merged@meta.data
meta_df$sample_group <- meta_df$orig.ident

# Map samples to developmental stages
stage_mapping <- list(
    "No LR"      = c("sample_1", "sample_12", "sample_13"),
    "Stage I"    = "sample_11",
    "Stage II"   = "sample_4",
    "Stage II&III" = "sample_14",
    "Stage III"  = c("sample_5", "sample_9"),
    "Stage III&IV" = "sample_2",
    "Stage V"    = c("sample_8", "sample_10"),
    "Stage VII"  = c("sample_6", "sample_7")
)

for (stage_name in names(stage_mapping)) {
    meta_df$sample_group[meta_df$orig.ident %in% stage_mapping[[stage_name]]] <-
        stage_name
}

# Compute proportions
prop_sample_df <- meta_df %>%
    group_by(sample_group, celltype) %>%
    summarise(cell_number = n(), .groups = "drop") %>%
    group_by(sample_group) %>%
    mutate(proportion = cell_number / sum(cell_number))

# Match cluster colors to cell types
cluster_colors <- setNames(
    celltype_colors,
    levels(ST_merged$celltype)
)

# Bar plot
p_proportion <- ggplot(
    prop_sample_df,
    aes(x = sample_group, y = proportion, fill = celltype)
) +
    geom_bar(stat = "identity", width = 0.95) +
    scale_fill_manual(values = cluster_colors) +
    scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        expand = c(0, 0)
    ) +
    theme_classic(base_size = 12) +
    labs(x = "Developmental stage", y = "Cell-type proportion") +
    theme(
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6)
    )

ggsave(file.path(FIGURES_DIR, "cell_proportion.jpg"),
       plot = p_proportion, width = 120, height = 100,
       dpi = 600, units = "mm", device = "jpeg")

# Save updated object
saveRDS(ST_merged, file = file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

message("03-markers-cell-proportion.R: Complete.")
