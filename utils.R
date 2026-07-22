# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  utils.R
# Purpose: Shared utility functions for the project
# Input:   None (functions only)
# Output:  None (functions sourced by other scripts)
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

# ------------------------------------------------------------------------------
# GO Enrichment and Visualization (for hdWGCNA modules)
# ------------------------------------------------------------------------------

#' Run GO enrichment analysis for a module and generate bar plot
#'
#' @param module_genes Character vector of gene IDs
#' @param module_name  Module name/color (e.g., "green", "turquoise")
#' @param org_db       Organism database (e.g., org.Wm82.eg.db)
#' @param ont         Ontology (default "BP")
#' @param key_type    Key type (default "GID")
#' @param pvalue_cutoff P-value cutoff (default 0.05)
#' @param simplify_cutoff Similarity cutoff for GO simplification (default 0.7)
#' @param top_n       Number of top terms to plot (default 5)
#' @param bar_fill    Fill color for bar plot (default based on module_name)
#' @param text_color  Text color for bar plot (default "black")
#' @param output_dir  Output directory for CSV and figures
#' @param fig_width   Figure width in inches
#' @param fig_height  Figure height in inches
#' @param fig_dpi     Figure resolution
#'
#' @return Invisibly returns the enrichment result

run_module_go_enrichment <- function(
    module_genes,
    module_name,
    org_db,
    ont = "BP",
    key_type = "GID",
    pvalue_cutoff = 0.05,
    simplify_cutoff = 0.7,
    top_n = 5,
    bar_fill = NULL,
    text_color = "black",
    output_dir = "results",
    fig_width = 10,
    fig_height = 6,
    fig_dpi = 500
) {
    # Run enrichment
    ego <- enrichGO(
        gene = module_genes,
        OrgDb = org_db,
        keyType = key_type,
        ont = ont,
        readable = FALSE,
        pvalueCutoff = pvalue_cutoff
    )

    # Save full results
    csv_path <- file.path(output_dir, "tables",
                          paste0("hdWGCNA_", module_name, "_module_BP.csv"))
    dir.create(dirname(csv_path), showWarnings = FALSE, recursive = TRUE)
    write.csv(as.data.frame(ego), file = csv_path, row.names = FALSE)

    # Simplify
    ego_simplified <- simplify(ego, cutoff = simplify_cutoff,
                               by = "p.adjust", select_fun = min)

    # Prepare data frame for plotting
    df <- as.data.frame(ego_simplified) %>%
        top_n(top_n, wt = -p.adjust) %>%
        arrange(p.adjust)

    df$log_p <- -log10(df$p.adjust)
    df$Description <- factor(df$Description, levels = rev(df$Description))

    # Default fill color based on module name
    if (is.null(bar_fill)) {
        bar_fill <- module_name
    }

    # Generate bar plot
    p <- ggplot(df, aes(x = Description, y = log_p)) +
        geom_bar(stat = "identity", fill = bar_fill, width = 0.8) +
        geom_text(aes(label = Description, y = 0.05),
                  hjust = 0, size = 8, color = text_color) +
        coord_flip() +
        theme_classic() +
        theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.title.x = element_text(size = 18)
        ) +
        labs(y = expression(log[10](Enrichment))) +
        scale_y_continuous(expand = c(0, 0),
                           limits = c(0, max(df$log_p) * 1.1))

    # Save plot
    fig_path <- file.path(output_dir, "figures",
                          paste0("go_", module_name, ".jpg"))
    dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
    ggsave(fig_path, plot = p, width = fig_width, height = fig_height,
           dpi = fig_dpi, device = "jpeg")

    invisible(ego)
}


# ------------------------------------------------------------------------------
# Monocle3 Trajectory Analysis
# ------------------------------------------------------------------------------

#' Run full monocle3 trajectory analysis pipeline
#'
#' @param seurat_obj Seurat object
#' @param root_cluster Cluster ID to use as root
#' @param num_dim Number of dimensions for preprocessing
#' @param reduction Which reduction to use for UMAP
#' @param cluster_k kNN parameter for clustering
#' @param cluster_resolution Resolution for clustering
#' @param min_branch_len Minimal branch length
#' @param euclidean_ratio Euclidean distance ratio
#'
#' @return List containing cds object and pseudotime DEG results

run_monocle3_trajectory <- function(
    seurat_obj,
    root_cluster,
    num_dim = 10,
    reduction = "umap",
    cluster_k = 40,
    cluster_resolution = 0.7,
    min_branch_len = 60,
    euclidean_ratio = 10
) {
    cds <- as.cell_data_set(seurat_obj)
    colData(cds)$seurat_clusters <- seurat_obj@meta.data$seurat_clusters

    cds <- preprocess_cds(cds, num_dim = num_dim)
    cds <- reduce_dimension(cds, reduction_method = "UMAP")
    reducedDims(cds)$UMAP <- Embeddings(seurat_obj, reduction = reduction)

    cds <- cluster_cells(cds, cluster_method = "leiden",
                         k = cluster_k, resolution = cluster_resolution)

    cds <- learn_graph(cds,
                       learn_graph_control = list(
                           minimal_branch_len = min_branch_len,
                           euclidean_distance_ratio = euclidean_ratio
                       ))

    root_cells <- colnames(cds)[
        colData(cds)$seurat_clusters == as.character(root_cluster)
    ]
    cds <- order_cells(cds, root_cells = root_cells)

    deg <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

    list(cds = cds, deg = deg)
}


# ------------------------------------------------------------------------------
# Pseudotime Gene Module Analysis
# ------------------------------------------------------------------------------


run_pseudotime_module_analysis <- function(
    cds,
    deg_results,
    q_value_threshold = 0.01,
    resolution = c(10^seq(-6, -1)),
    output_prefix = "module"
) {
    res_ids <- row.names(subset(deg_results, q_value < q_value_threshold))
    gene_modules <- find_gene_modules(cds[res_ids, ], resolution = resolution)

    csv_path <- paste0(output_prefix, "_gene_expression_module_analysis.csv")
    write.csv(gene_modules, file = csv_path)

    cell_group_df <- tibble::tibble(
        cell = row.names(colData(cds)),
        cell_group = colData(cds)$seurat_clusters
    )

    agg_mat <- aggregate_gene_expression(cds, gene_modules, cell_group_df)
    row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))

    p <- pheatmap::pheatmap(
        agg_mat,
        scale = "column",
        clustering_method = "ward.D2"
    )

    list(
        gene_modules = gene_modules,
        agg_mat = agg_mat,
        heatmap = p,
        cell_group_df = cell_group_df
    )
}


# ------------------------------------------------------------------------------
# GO Enrichment Loop for Pseudotime Modules
# ------------------------------------------------------------------------------


run_module_go_loop <- function(
    module_genes_list,
    org_db,
    num_modules,
    key_type = "GID",
    ont = "BP",
    p_adjust_method = "BH",
    qvalue_cutoff = 0.05,
    output_dir = "GO_results",
    plot_dir = "GO_plots",
    plot_width = 8,
    plot_height = 6
) {
    dir.create(output_dir, showWarnings = FALSE)
    dir.create(plot_dir, showWarnings = FALSE)

    for (i in seq_len(num_modules)) {
        ego <- enrichGO(
            gene = module_genes_list[[i]],
            OrgDb = org_db,
            keyType = key_type,
            ont = ont,
            pAdjustMethod = p_adjust_method,
            qvalueCutoff = qvalue_cutoff
        )

        if (!is.null(ego) && nrow(ego) > 0) {
            csv_path <- file.path(output_dir, paste0("module_", i, "_GO.csv"))
            write.csv(as.data.frame(ego), file = csv_path, row.names = FALSE)

            p <- dotplot(ego, showCategory = 20) +
                ggtitle(paste0("Module ", i, " GO BP"))

            plot_path <- file.path(plot_dir,
                                   paste0("module_", i, "_dotplot.png"))
            ggsave(filename = plot_path, plot = p,
                   width = plot_width, height = plot_height)
        }

        cat("Module", i, "done.\n")
    }
}


# ------------------------------------------------------------------------------
# Module Score Loop (AddModuleScore + Violin Plot)
# ------------------------------------------------------------------------------


run_module_score_loop <- function(
    seurat_obj,
    gene_sets_list,
    color_vector,
    output_dir = "results",
    plot_width = 5,
    plot_height = 5,
    dpi = 300
) {
    for (set_name in names(gene_sets_list)) {
        score_name <- paste0(set_name, "Score1")
        seurat_obj <- AddModuleScore(
            seurat_obj,
            features = list(gene_sets_list[[set_name]]),
            name = paste0(set_name, "Score")
        )

        p <- VlnPlot(
            seurat_obj,
            features = score_name,
            pt.size = 0,
            cols = color_vector
        ) +
            ggtitle(paste0(set_name, " pathway activity")) +
            theme_classic() +
            theme(
                axis.title.x = element_blank(),
                plot.title = element_text(size = 16),
                axis.text.y = element_text(size = 14),
                axis.text.x = element_text(angle = 90, size = 14, hjust = 1),
                legend.position = "none"
            )

        fig_path <- file.path(output_dir, "figures",
                              paste0(set_name, "_violin_plot.jpg"))
        dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
        ggsave(fig_path, plot = p, width = plot_width,
               height = plot_height, dpi = dpi)
        print(p)
    }

    seurat_obj
}


# ------------------------------------------------------------------------------
# Spatial Feature Score Loop
# ------------------------------------------------------------------------------


run_spatial_score_loop <- function(
    seurat_obj,
    gene_sets_list,
    output_dir = "results",
    pt_size_factor = 16,
    plot_width = 30,
    plot_height = 5,
    dpi = 300
) {
    for (set_name in names(gene_sets_list)) {
        score_name <- paste0(set_name, "Score1")
        p <- SpatialFeaturePlot(
            seurat_obj,
            features = score_name,
            pt.size.factor = pt_size_factor
        ) +
            ggtitle(paste0(set_name, " pathway spatial distribution"))

        fig_path <- file.path(output_dir, "figures",
                              paste0(set_name, "_spatial_plot.jpg"))
        dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
        ggsave(fig_path, plot = p, width = plot_width,
               height = plot_height, dpi = dpi)
        print(p)
    }
}
