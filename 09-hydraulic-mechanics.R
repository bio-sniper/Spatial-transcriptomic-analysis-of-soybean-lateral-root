# ==============================================================================
# Project: Soybean Lateral Root Spatial Transcriptome
# Script:  09-hydraulic-mechanics.R
# Purpose: Hydraulic/mechanical force analysis and in silico perturbation
# Input:   results/objects/ST_LR_merged.rds
#          correlation Excel files (external data)
# Output:  results/figures/*.jpg (mechanical force plots)
#          results/tables/*.csv
# Author:  Lei Li
# Date:    2026-07-16
# ==============================================================================

source("scripts/01-setup.R")
library(readxl)

ST_merged <- readRDS(file.path(OBJECTS_DIR, "ST_LR_merged.rds"))

# ------------------------------------------------------------------------------
# Water-related gene module score
# ------------------------------------------------------------------------------

water_genes <- c(
    "Glyma.01G208200", "Glyma.03G185900", "Glyma.02G255000",
    "Glyma.13G325900", "Glyma.14G061500"
)

ST_merged <- AddModuleScore(
    ST_merged,
    features = list(water_genes),
    name = "water_Score"
)

# ------------------------------------------------------------------------------
# Asymmetry mechanics analysis function
# ------------------------------------------------------------------------------

analyze_asymmetry_mechanics <- function(chip_id, seurat_obj) {
    message(">>> Processing: ", chip_id)

    cells_to_keep <- rownames(seurat_obj@meta.data)[
        seurat_obj@meta.data$orig.ident == chip_id
    ]
    if (length(cells_to_keep) == 0) return(NULL)

    single_slice <- subset(seurat_obj, cells = cells_to_keep)
    coords <- GetTissueCoordinates(single_slice)

    water_data <- na.omit(data.frame(
        x = coords[, 1],
        y = coords[, 2],
        score = single_slice$water_Score1
    ))

    grid_res <- 150
    smoothed_surface <- interp(
        x = water_data$x, y = water_data$y, z = water_data$score,
        nx = grid_res, ny = grid_res, linear = FALSE
    )

    Z_matrix_smooth <- image.smooth(smoothed_surface$z, theta = 4)$z
    grad_field <- gradient(Z_matrix_smooth)

    vector_data <- data.frame(
        x = rep(smoothed_surface$x, times = grid_res),
        y = rep(smoothed_surface$y, each = grid_res),
        score = as.vector(Z_matrix_smooth),
        dx = as.vector(grad_field$X),
        dy = as.vector(grad_field$Y)
    )
    vector_data$force_magnitude <- sqrt(
        vector_data$dx^2 + vector_data$dy^2
    )

    real_coords <- water_data[, c("x", "y")]
    nn_dist <- get.knnx(
        data = real_coords,
        query = vector_data[, c("x", "y")],
        k = 1
    )$nn.dist

    tol_radius <- 0.03 * (max(real_coords$x) - min(real_coords$x))
    vector_data[nn_dist > tol_radius,
                c("score", "dx", "dy", "force_magnitude")] <- NA

    y_center <- median(vector_data$y[
        !is.na(vector_data$score)
    ])
    force_dorsal <- mean(
        vector_data$force_magnitude[vector_data$y > y_center],
        na.rm = TRUE
    )
    force_ventral <- mean(
        vector_data$force_magnitude[vector_data$y <= y_center],
        na.rm = TRUE
    )
    asymmetry_index <- force_dorsal / (force_ventral + 1e-6)

    nn_map <- get.knnx(
        data = vector_data[, c("x", "y")],
        query = real_coords,
        k = 1
    )
    single_slice$Mechanical_Force <- vector_data$force_magnitude[
        nn_map$nn.index
    ]

    threshold <- quantile(
        vector_data$force_magnitude, 0.25, na.rm = TRUE
    )

    p <- SpatialFeaturePlot(
        single_slice,
        features = "Mechanical_Force",
        pt.size.factor = 8,
        crop = FALSE
    ) +
        scale_fill_viridis_c(option = "magma", name = "Force") +
        metR::geom_streamline(
            data = dplyr::filter(
                vector_data,
                force_magnitude > threshold & !is.na(score)
            ),
            aes(x = x, y = y, dx = dx, dy = dy),
            inherit.aes = FALSE,
            L = 15, res = 2, size = 0.5, color = "white",
            arrow = arrow(length = unit(0.05, "inches"))
        ) +
        labs(title = paste(chip_id, "| AI:",
                           round(asymmetry_index, 2))) +
        theme(plot.title = element_text(size = 12, hjust = 0.5))

    list(
        obj          = single_slice,
        plot         = p,
        AI           = asymmetry_index,
        mean_force   = force_dorsal
    )
}


# ------------------------------------------------------------------------------
# Run asymmetry analysis across developmental stages
# ------------------------------------------------------------------------------

target_chips <- paste0("sample", 1:13)

time_mapping <- c(
    "sample1"  = 0, "sample2"  = 4, "sample3"  = 2,
    "sample4"  = 3, "sample5"  = 7, "sample6"  = 7,
    "sample7"  = 5, "sample8"  = 3, "sample9"  = 5,
    "sample10" = 1, "sample11" = 0, "sample12" = 0,
    "sample13" = 3
)

results_list <- list()
model2_df <- data.frame(
    Stage = character(),
    Time  = numeric(),
    AI    = numeric(),
    stringsAsFactors = FALSE
)

for (chip in target_chips) {
    res <- analyze_asymmetry_mechanics(chip, ST_merged)
    if (!is.null(res)) {
        results_list[[chip]] <- res
        model2_df <- rbind(model2_df, data.frame(
            Stage = chip,
            Time  = unname(time_mapping[chip]),
            AI    = res$AI
        ))
    }
}

# ------------------------------------------------------------------------------
# Mechanical force combined plot
# ------------------------------------------------------------------------------

p_force_combined <- wrap_plots(
    lapply(results_list, function(x) x$plot),
    ncol = 3
)

ggsave(file.path(FIGURES_DIR, "mechanical_force_all_stages.jpg"),
       plot = p_force_combined, width = 200, height = 240,
       units = "mm", dpi = 1200, device = "jpeg")

# ------------------------------------------------------------------------------
# Model 2: Asymmetric force burst threshold
# ------------------------------------------------------------------------------

model2_df <- model2_df[order(model2_df$Time), ]

smooth_fit <- smooth.spline(model2_df$Time, model2_df$AI, df = 3)

t_seq <- seq(min(model2_df$Time), max(model2_df$Time), length.out = 100)
pred_AI <- predict(smooth_fit, t_seq)$y
deriv_AI <- predict(smooth_fit, t_seq, deriv = 1)$y
rupture_point <- t_seq[which.max(deriv_AI)]

p_asymmetry <- ggplot() +
    geom_line(
        data = data.frame(t = t_seq, v = pred_AI),
        aes(x = t, y = v),
        color = "#d95f02", size = 1.2
    ) +
    geom_vline(
        xintercept = rupture_point,
        linetype = "dashed",
        color = "#1b9e77"
    ) +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
    ) +
    scale_x_continuous(breaks = seq(0, 7, by = 1)) +
    labs(
        title = "Model 2: Asymmetric Force Burst Threshold",
        subtitle = paste("Critical asymmetric drive at Stage",
                         round(rupture_point, 2)),
        x = "Developmental Progression",
        y = "Asymmetry Index (Dorsal/Ventral)"
    )

ggsave(file.path(FIGURES_DIR, "asymmetric_force_burst.jpg"),
       plot = p_asymmetry, width = 4, height = 5,
       dpi = 500, device = "jpeg")

# ------------------------------------------------------------------------------
# Mechanosensitive gene correlation (sample2 = Stage III&IV)
# ------------------------------------------------------------------------------

tryCatch({
target_chip_name <- "sample2"
chip_data <- results_list[[target_chip_name]]

if ("obj" %in% names(chip_data)) {
    target_obj <- chip_data$obj
} else if ("seurat_obj" %in% names(chip_data)) {
    target_obj <- chip_data$seurat_obj
} else {
    stop("Seurat object not found in results_list for sample2")
}

if (!"SCT" %in% Assays(target_obj)) {
    stop("SCT assay not found in target object")
}

expr_matrix <- GetAssayData(target_obj, assay = "SCT")

if (!"Mechanical_Force" %in% colnames(target_obj@meta.data)) {
    stop("Mechanical_Force column not found in metadata")
}
forces <- target_obj$Mechanical_Force

genes_to_test <- rownames(expr_matrix)[
    rowSums(expr_matrix > 0) > ncol(expr_matrix) * 0.05
]
filtered_expr <- as.matrix(expr_matrix[genes_to_test, ])
message(length(genes_to_test), " candidate genes for mechanical coupling")

cor_results <- pbapply::pblapply(rownames(filtered_expr), function(gene) {
    x <- filtered_expr[gene, ]
    test_res <- tryCatch(
        cor.test(x, forces, method = "spearman", exact = FALSE),
        error = function(e) list(estimate = NA_real_, p.value = NA_real_)
    )
    rho <- if (is.null(test_res$estimate)) NA_real_ else as.numeric(test_res$estimate)
    pv  <- if (is.null(test_res$p.value)) NA_real_ else test_res$p.value
    data.frame(Gene = gene, Spearman_Rho = rho, p_value = pv, stringsAsFactors = FALSE)
})

mechano_df <- do.call(rbind, cor_results)
mechano_df$p_adj <- p.adjust(mechano_df$p_value, method = "BH")
mechano_df <- mechano_df[!is.na(mechano_df$Spearman_Rho), ]
if (nrow(mechano_df) > 0) {
    significant_df <- mechano_df[mechano_df$p_adj < 0.05, ]
    significant_df <- significant_df[order(-significant_df$Spearman_Rho), ]
    print(head(significant_df, 15))
} else {
    message("No significant mechanosensitive genes found")
}
}, error = function(e) message("Mechanosensitive correlation skipped: ", e$message))

# ------------------------------------------------------------------------------
# Correlation heatmap across stages (skip if Excel files not found)
# ------------------------------------------------------------------------------

excel_dir <- "G:/open code project/project_4/data"
files <- list(
    "Stage_0"   = file.path(excel_dir, "correlation of sample 1 none LR.xlsx"),
    "Stage_IV"  = file.path(excel_dir, "correlation of sample 2_stage IV.xlsx"),
    "Stage_V"   = file.path(excel_dir, "correlation of sample 7_stage V.xlsx"),
    "Stage_VII" = file.path(excel_dir, "correlation of sample 5_stage VII.xlsx")
)

all_exist <- all(sapply(files, file.exists))
if (all_exist) {
    df_list <- list()
    for (stage in names(files)) {
        temp_df <- read_excel(files[[stage]])
        colnames(temp_df)[1:2] <- c("Gene", "Spearman_Rho")
        temp_df$Stage <- stage
        df_list[[stage]] <- temp_df[1:15, ]
    }

    merged_df <- bind_rows(df_list)

    matrix_df <- merged_df %>%
        select(Gene, Stage, Spearman_Rho) %>%
        pivot_wider(names_from = Stage,
                    values_from = Spearman_Rho,
                    values_fill = 0) %>%
        column_to_rownames("Gene")

    time_order <- c("Stage_0", "Stage_IV", "Stage_V", "Stage_VII")
    matrix_df <- matrix_df[, time_order]

    col_fun <- colorRamp2(c(0, max(matrix_df)), c("skyblue", "brown"))

    ht <- Heatmap(
        as.matrix(matrix_df),
        name = "Spearman Rho",
        col = col_fun,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 8),
        column_title = "Dynamic Shift of Mechanosensitive Genes",
        row_title = "Force-Coupled Genes"
    )

    draw(ht)

    # --------------------------------------------------------------------------
    # Bubble plot of mechanosensitive genes
    # --------------------------------------------------------------------------

    merged_df$label <- ifelse(
        is.na(merged_df$At_name) | merged_df$At_name == "" |
            grepl("^AT", merged_df$At_name),
        merged_df$Gene,
        paste0("Gm", merged_df$At_name)
    )

    p_bubble <- ggplot(
        merged_df,
        aes(x = factor(Stage, levels = time_order), y = Gene)
    ) +
        geom_point(aes(size = Spearman_Rho, color = Spearman_Rho)) +
        scale_color_viridis_c(option = "magma") +
        theme_bw() +
        labs(
            title = "Temporal Activation of Mechanotransduction Hubs",
            x = "Developmental Progression",
            y = ""
        ) +
        scale_y_discrete(labels = setNames(merged_df$label, merged_df$Gene)) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(hjust = 1, face = "italic", size = 12)
        )

    ggsave(file.path(FIGURES_DIR, "correlation_genes_mechanical.jpg"),
           plot = p_bubble, width = 5, height = 10,
           dpi = 500, device = "jpeg")
} else {
    message("Excel correlation files not found, skipping heatmap and bubble plot")
    message("Missing:", paste(names(files)[!sapply(files, file.exists)], collapse=", "))
}

# ------------------------------------------------------------------------------
# LOX1 module score
# ------------------------------------------------------------------------------

LOX1_genes <- c(
    "Glyma.03G237300", "Glyma.07G006900",
    "Glyma.07G007000", "Glyma.07G034800"
)

ST_merged <- AddModuleScore(
    ST_merged,
    features = list(LOX1_genes),
    name = "LOX1_Score"
)

p_lox1 <- SpatialFeaturePlot(
    ST_merged,
    features = "LOX1_Score1",
    crop = FALSE,
    pt.size.factor = 8,
    ncol = 3
) + coord_fixed()

ggsave(file.path(FIGURES_DIR, "LOX1_spatial_score.jpg"),
       plot = p_lox1, width = 160, height = 200,
       units = "mm", dpi = 1200, device = "jpeg")

# Sample2 stage-specific analysis
stage2_obj <- subset(ST_merged, orig.ident == "sample2")

p_stage2_vln <- VlnPlot(
    stage2_obj,
    features = c("Glyma.07G007000", "Glyma.08G230500"),
    cols = celltype_colors,
    pt.size = 0
)

ggsave(file.path(FIGURES_DIR, "stage2_LOX1_vlnplot.jpg"),
       plot = p_stage2_vln, width = 8, height = 6, dpi = 300)

# ==============================================================================
# In silico knockout simulation of LDL2
# ==============================================================================

local_fold_change <- 1.5
samples <- unique(ST_merged$orig.ident)

all_simulated_force <- data.frame(
    row.names = colnames(ST_merged),
    Force_Magnitude_KO = rep(NA, ncol(ST_merged))
)

plot_list <- list()

for (sam in samples) {
    message("Processing sample: ", sam)

    target_obj <- subset(ST_merged, orig.ident == sam)
    coords <- GetTissueCoordinates(target_obj)

    st_data <- data.frame(
        spot_id = rownames(coords),
        x = coords[, 1],
        y = coords[, 2],
        water_Score = target_obj$water_Score1,
        is_LRP = ifelse(Idents(target_obj) == "Lateral root", 1, 0)
    )

    # Apply fold change to LRP spots
    st_data_ko <- st_data
    st_data_ko$water_Score <- ifelse(
        st_data_ko$is_LRP == 1,
        st_data_ko$water_Score * local_fold_change,
        st_data_ko$water_Score
    )

    # GAM smoothing
    k_val <- min(300, floor(nrow(st_data_ko) * 0.8))
    gam_model_ko <- gam(
        water_Score ~ s(x, y, k = k_val),
        data = st_data_ko,
        method = "REML"
    )

    grid_res <- 150
    x_seq <- seq(min(st_data_ko$x), max(st_data_ko$x), length.out = grid_res)
    y_seq <- seq(min(st_data_ko$y), max(st_data_ko$y), length.out = grid_res)
    dense_grid <- expand.grid(x = x_seq, y = y_seq)
    dense_grid$score <- predict(gam_model_ko, newdata = dense_grid)

    Z_matrix <- matrix(dense_grid$score, nrow = grid_res, ncol = grid_res)

    # Compute gradient (force)
    grad_X <- matrix(0, nrow = grid_res, ncol = grid_res)
    grad_Y <- matrix(0, nrow = grid_res, ncol = grid_res)

    for (i in 2:(grid_res - 1)) {
        for (j in 2:(grid_res - 1)) {
            grad_X[i, j] <- (Z_matrix[i + 1, j] - Z_matrix[i - 1, j]) / 2
            grad_Y[i, j] <- (Z_matrix[i, j + 1] - Z_matrix[i, j - 1]) / 2
        }
    }

    force_matrix <- sqrt(grad_X^2 + grad_Y^2)
    grid_list <- list(x = x_seq, y = y_seq, z = force_matrix)
    spot_coords <- cbind(st_data$x, st_data$y)
    st_data$simulated_force_ko <- fields::interp.surface(grid_list, spot_coords)

    all_simulated_force[st_data$spot_id, "Force_Magnitude_KO"] <-
        st_data$simulated_force_ko

    target_obj <- AddMetaData(
        target_obj,
        metadata = st_data$simulated_force_ko,
        col.name = "Force_Magnitude_KO"
    )

    p <- SpatialFeaturePlot(
        target_obj,
        features = "Force_Magnitude_KO",
        pt.size.factor = 8,
        crop = FALSE
    ) +
        scale_fill_viridis_c(option = "magma", name = "Force (KO)") +
        ggtitle(sam) +
        theme(plot.title = element_text(size = 14,
                                        face = "bold", hjust = 0.5))

    plot_list[[sam]] <- p
}

ST_merged <- AddMetaData(
    ST_merged,
    metadata = all_simulated_force$Force_Magnitude_KO,
    col.name = "Force_Magnitude_KO"
)

# ------------------------------------------------------------------------------
# Combined KO simulation plot
# ------------------------------------------------------------------------------

p_ko_combined <- wrap_plots(plot_list, ncol = 4) +
    plot_annotation(
        title = paste(
            "In silico Knockout of LDL2 (LRP Fold Change =",
            local_fold_change, ")"
        ),
        subtitle = paste(
            "Pan-developmental Spatial Force Dynamics Across",
            length(samples), "Stages"
        ),
        theme = theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 16, hjust = 0.5)
        )
    )

ggsave(file.path(FIGURES_DIR, "LDL2_KO_mechanical_force.jpg"),
       plot = p_ko_combined, width = 500, height = 600,
       units = "mm", dpi = 500, device = "jpeg")

message("09-hydraulic-mechanics.R: Complete.")
