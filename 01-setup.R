# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  01-setup.R
# Purpose: Load required packages, set global options, define paths
# Input:   None
# Output:  None (environment setup only)
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

# ------------------------------------------------------------------------------
# Global options
# ------------------------------------------------------------------------------

options(future.globals.maxSize = 10 * 1024^3)
set.seed(42)

# ------------------------------------------------------------------------------
# Path configuration
# ------------------------------------------------------------------------------

BASE_DIR <- getwd()
OUTPUT_BASE <- Sys.getenv("OUTPUT_DIR", unset = file.path(BASE_DIR, "results"))
RESULTS_DIR <- normalizePath(OUTPUT_BASE, mustWork = FALSE)
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
OBJECTS_DIR <- file.path(RESULTS_DIR, "objects")

dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OBJECTS_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# Package loading
# ------------------------------------------------------------------------------

# Core data science
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

# Integration
library(harmony)

# Trajectory
library(monocle3)
library(SeuratWrappers)
library(igraph)
library(matrixStats)

# Annotation & enrichment
library(clusterProfiler)
library(enrichplot)
library(org.Wm82.eg.db)

# Visualization
library(RColorBrewer)
library(pals)
library(ggsci)
library(scRNAtoolVis)
library(viridis)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(metR)

# hdWGCNA
library(hdWGCNA)
library(WGCNA)
library(UCell)

# Spatial interpolation
library(interp)
library(akima)
library(fields)
library(FNN)

# Modeling
library(pracma)
library(minpack.lm)
library(pbapply)
library(mgcv)

# ------------------------------------------------------------------------------
# Custom color palettes
# ------------------------------------------------------------------------------

cluster_colors_24 <- c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928",
    "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#ffff99",
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"
)

celltype_colors <- c(
    "Cortex-Endodermis"    = "#1f78b4",
    "Epidermis"            = "#33a02c",
    "Pericycle"            = "#e31a1c",
    "Xylem"                = "#ff7f00",
    "Lateral root"         = "#6a3d9a",
    "Phloem"               = "#b15928",
    "Lateral root junction" = "#a6cee3"
)

# ------------------------------------------------------------------------------
# Gene name mapping (Glyma → common name)
# ------------------------------------------------------------------------------

glyma_to_common <- c(
    "Glyma.03G132700" = "GmGNS1",
    "Glyma.08G150400" = "GmBGLU17",
    "Glyma.16G170100" = "GmEXT3",
    "Glyma.09G049100" = "GmCYP81D3",
    "Glyma.19G263100" = "GmPOT5",
    "Glyma.02G064100" = "GmRNS1",
    "Glyma.15G143600" = "GmCTL2",
    "Glyma.09G038500" = "GmCTL2",
    "Glyma.04G014800" = "GmXCP2",
    "Glyma.11G254700" = "GmSAHH1",
    "Glyma.07G038600" = "GmGER3",
    "Glyma.09G022500" = "GmPER54",
    "Glyma.16G164400" = "GmPER12",
    "Glyma.11G189500" = "GmHDT2",
    "Glyma.06G242900" = "GmMT3",
    "Glyma.05G057200" = "GmAZI5",
    "Glyma.09G238300" = "GmGASA4",
    "Glyma.12G207600" = "GmFLA9"
)

message("01-setup.R: Setup complete. Output directories configured.")
