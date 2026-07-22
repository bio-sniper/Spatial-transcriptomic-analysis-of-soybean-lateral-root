# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  08-phytohormones.R
# Purpose: Phytohormone-related gene module scoring and visualization
# Input:   results/objects/ST_LR_merged.rds
# Output:  results/figures/*.jpg (violin, spatial plots)
#          results/tables/Hormone_ModuleScores.csv
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")
source("scripts/utils.R")

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Auxin-related gene sets
# ------------------------------------------------------------------------------

ARF_genes <- c(
    "Glyma.12G164100", "Glyma.16G000300", "Glyma.05G200800",
    "Glyma.13G174000", "Glyma.12G171000", "Glyma.14G217700",
    "Glyma.08G100100", "Glyma.17G047100", "Glyma.02G239600",
    "Glyma.14G208500", "Glyma.03G070500", "Glyma.11G145500",
    "Glyma.12G076200", "Glyma.10G210600", "Glyma.19G181900",
    "Glyma.07G054800", "Glyma.15G181000", "Glyma.07G130400"
)

IAA_genes <- c(
    "Glyma.02G142500", "Glyma.01G098000", "Glyma.10G180100",
    "Glyma.19G161100", "Glyma.20G210400", "Glyma.13G361100",
    "Glyma.19G221900", "Glyma.06G067700", "Glyma.01G019400",
    "Glyma.13G127000", "Glyma.20G225000", "Glyma.13G354100",
    "Glyma.08G036400"
)

auxin_other_signal_genes <- c(
    "Glyma.14G185600", "Glyma.17G247700", "Glyma.08G214600",
    "Glyma.13G367300", "Glyma.10G150000", "Glyma.17G112500"
)

auxin_biosynthesis_genes <- c(
    "Glyma.17G086500", "Glyma.10G128700", "Glyma.20G080000",
    "Glyma.03G169600", "Glyma.19G170800", "Glyma.19G206200",
    "Glyma.19G160800"
)

auxin_conjugation_genes <- c(
    "Glyma.08G197000", "Glyma.13G304000", "Glyma.06G141200"
)

auxin_transport_genes <- c(
    "Glyma.19G184300", "Glyma.03G183600", "Glyma.08G342600",
    "Glyma.19G021500", "Glyma.13G063700", "Glyma.10G290800",
    "Glyma.02G008000", "Glyma.18G017600", "Glyma.13G355000",
    "Glyma.08G201300", "Glyma.17G039200", "Glyma.18G198400",
    "Glyma.03G063600", "Glyma.20G094600", "Glyma.13G101900",
    "Glyma.U031115", "Glyma.05G070600", "Glyma.09G195600",
    "Glyma.09G116100", "Glyma.11G088600", "Glyma.09G271100",
    "Glyma.07G102500", "Glyma.08G054700", "Glyma.03G126000",
    "Glyma.14G216200", "Glyma.04G027800"
)

# Merge all auxin-related genes
auxin_gene_sets <- list(
    ARF               = ARF_genes,
    IAA               = IAA_genes,
    auxin_other_signal = auxin_other_signal_genes,
    auxin_conjugation  = auxin_conjugation_genes,
    auxin_biosynthesis = auxin_biosynthesis_genes,
    auxin_transport    = auxin_transport_genes
)

auxin_all_genes <- unique(unlist(auxin_gene_sets))

# ------------------------------------------------------------------------------
# Other hormone gene sets
# ------------------------------------------------------------------------------

cytokinin_genes <- c(
    "Glyma.05G033000", "Glyma.17G093900", "Glyma.04G062500",
    "Glyma.06G063500", "Glyma.09G098600", "Glyma.13G155400",
    "Glyma.14G110600", "Glyma.15G206200", "Glyma.17G076000",
    "Glyma.17G217100", "Glyma.07G253100", "Glyma.10G157800",
    "Glyma.13G212800", "Glyma.15G099800", "Glyma.20G230800",
    "Glyma.05G148100", "Glyma.08G105000", "Glyma.05G241600",
    "Glyma.08G049000", "Glyma.15G145200", "Glyma.07G079000",
    "Glyma.07G243300"
)

ABA_genes <- c(
    "Glyma.04G180400", "Glyma.08G036600", "Glyma.08G230600",
    "Glyma.12G217300", "Glyma.12G217400", "Glyma.14G140900",
    "Glyma.05G227100", "Glyma.08G033800", "Glyma.09G066500",
    "Glyma.13G106800", "Glyma.15G172500", "Glyma.02G131700",
    "Glyma.04G039300", "Glyma.06G040400", "Glyma.07G213100"
)

JA_genes <- c(
    "Glyma.01G018400", "Glyma.01G096600", "Glyma.08G271900",
    "Glyma.09G204500", "Glyma.17G209000", "Glyma.13G109700",
    "Glyma.17G050000"
)

SA_genes <- c(
    "Glyma.03G127600", "Glyma.03G128200", "Glyma.19G130200",
    "Glyma.04G223300", "Glyma.06G142000", "Glyma.09G274000",
    "Glyma.13G267400", "Glyma.13G267500", "Glyma.13G267700"
)

BR_genes <- c(
    "Glyma.11G067700", "Glyma.17G248900", "Glyma.11G103500",
    "Glyma.12G028300", "Glyma.04G218300", "Glyma.06G147600",
    "Glyma.01G178000", "Glyma.11G064300", "Glyma.04G063600",
    "Glyma.06G064800", "Glyma.16G180600", "Glyma.04G033800",
    "Glyma.06G034000", "Glyma.14G076900", "Glyma.04G225700",
    "Glyma.06G139100"
)

ga_genes <- c(
    "Glyma.02G151100", "Glyma.03G148300", "Glyma.10G022900",
    "Glyma.09G149200", "Glyma.16G200800", "Glyma.04G071000",
    "Glyma.13G361700", "Glyma.15G012100", "Glyma.17G205300",
    "Glyma.05G140400", "Glyma.08G095800", "Glyma.10G241100",
    "Glyma.20G230600"
)

ET_genes <- c(
    "Glyma.09G002600", "Glyma.12G241700", "Glyma.10G188500",
    "Glyma.20G087000", "Glyma.20G202200", "Glyma.03G181400",
    "Glyma.10G058300", "Glyma.13G145100", "Glyma.15G079100",
    "Glyma.08G018000", "Glyma.02G274600", "Glyma.03G245000",
    "Glyma.14G197100", "Glyma.18G059700", "Glyma.14G116800",
    "Glyma.17G211000", "Glyma.02G165800"
)

hormone_gene_sets <- list(
    Auxin     = auxin_all_genes,
    Cytokinin = cytokinin_genes,
    JA        = JA_genes,
    BR        = BR_genes,
    GA        = ga_genes,
    ET        = ET_genes,
    SA        = SA_genes,
    ABA       = ABA_genes
)

# ------------------------------------------------------------------------------
# Auxin sub-pathway scoring
# ------------------------------------------------------------------------------

ST_merged <- run_module_score_loop(
    seurat_obj    = ST_merged,
    gene_sets_list = auxin_gene_sets,
    color_vector  = celltype_colors,
    output_dir    = RESULTS_DIR
)

# ------------------------------------------------------------------------------
# Auxin sub-pathway spatial plots
# ------------------------------------------------------------------------------

run_spatial_score_loop(
    seurat_obj    = ST_merged,
    gene_sets_list = auxin_gene_sets,
    output_dir    = RESULTS_DIR,
    pt_size_factor = 16
)

# ------------------------------------------------------------------------------
# Global hormone scoring
# ------------------------------------------------------------------------------

ST_merged <- run_module_score_loop(
    seurat_obj    = ST_merged,
    gene_sets_list = hormone_gene_sets,
    color_vector  = celltype_colors,
    output_dir    = RESULTS_DIR
)

# ------------------------------------------------------------------------------
# Global hormone spatial plots
# ------------------------------------------------------------------------------

for (hormone in names(hormone_gene_sets)) {
    score_name <- paste0(hormone, "Score1")
    p <- SpatialFeaturePlot(ST_merged, features = score_name) +
        ggtitle(paste0(hormone, " pathway spatial distribution"))
    print(p)
}

# ------------------------------------------------------------------------------
# Export hormone scores
# ------------------------------------------------------------------------------

score_cols <- colnames(ST_merged@meta.data)[
    grepl("Score", colnames(ST_merged@meta.data))
]
hormone_scores <- FetchData(
    ST_merged,
    vars = c("orig.ident", score_cols)
)
write.csv(hormone_scores,
          file = file.path(TABLES_DIR, "Hormone_ModuleScores.csv"),
          row.names = TRUE)

# ------------------------------------------------------------------------------
# Additional score projections
# ------------------------------------------------------------------------------

ST_merged <- AddModuleScore(
    ST_merged,
    features = list(auxin_all_genes),
    name = "AuxinScore"
)
ST_merged <- AddModuleScore(
    ST_merged,
    features = list(cytokinin_genes),
    name = "CKScore"
)
ST_merged <- AddModuleScore(
    ST_merged,
    features = list(BR_genes),
    name = "BRScore"
)
ST_merged <- AddModuleScore(
    ST_merged,
    features = list(ABA_genes),
    name = "ABAScore"
)

# Spatial projections of combined scores
score_projections <- c("AuxinScore1", "ABAScore1", "BRScore1")
score_filenames  <- c("auxin_score_projection.jpg",
                      "ABA_score_projection.jpg",
                      "BR_score_projection.jpg")

for (i in seq_along(score_projections)) {
    p <- SpatialFeaturePlot(
        ST_merged,
        features = score_projections[i],
        crop = FALSE,
        pt.size.factor = 8,
        ncol = 3
    ) + coord_fixed()

    ggsave(file.path(FIGURES_DIR, score_filenames[i]),
           plot = p, units = "mm", height = 200, width = 160,
           dpi = 1200, device = "jpeg")
}

message("08-phytohormones.R: Complete.")
