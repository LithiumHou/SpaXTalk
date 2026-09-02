#' Generate a categorical color palette
#'
#' Select colors from the built-in SpaXTalk palette and interpolate additional
#' colors when more categories are requested.
#'
#' @param n Integer number of colors to generate.
#'
#' @return A character vector of `n` color values.
scPalette <- function(n) {
    colorSpace <- c(
        '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#F29403', '#F781BF',
        '#BC9DCC', '#A65628', '#54B0E4', '#222F75', '#1B9E77', '#B2DF8A',
        '#E3BE00', '#FB9A99', '#E7298A', '#910241', '#00CDD1', '#A6CEE3',
        '#CE1261', '#5E4FA2', '#8CA77B', '#00441B', '#DEDC00', '#DCF0B9',
        '#8DD3C7', '#999999'
    )
    if (n <= length(colorSpace)) {
        colors <- colorSpace[1:n]
    } else {
        colors <- grDevices::colorRampPalette(colorSpace)(n)
    }
    return(colors)
}

#' Plot a spatial variable
#'
#' Display one continuous or categorical variable over spatial cell coordinates,
#' with cells lacking an active value shown as a gray background.
#'
#' @param data A named numeric, character, or factor vector whose names are cell
#'   identifiers.
#' @param meta A cell metadata data frame containing `x` and `y` columns and cell
#'   identifiers as row names.
#' @param colors Optional vector of gradient colors for numeric data or a named
#'   palette for categorical data.
#' @param pt_size Numeric size of active points.
#' @param legend_name Optional character scalar used as the color legend title.
#' @param plot_zero Logical; whether numeric zero values should be plotted as
#'   active observations.
#'
#' @return A `ggplot` object showing the variable in spatial coordinates.
PlotSpa_Scatter <- function(data, meta, colors = NULL, pt_size = 4,
    legend_name = "Specificity", plot_zero = FALSE) {

    df <- meta
    cells <- intersect(names(data), rownames(df))
    if (is.numeric(data)) {
        df$value <- 0
        df[cells, "value"] <- data[cells]
        active_cells <- if (plot_zero) cells else cells[df[cells, "value"] != 0]
        df <- df[order(df$value, decreasing = T), ]
    } else {
        df$value <- NA_character_
        df[cells, "value"] <- as.character(data[cells])
        active_cells <- cells[!is.na(df[cells, "value"])]
        df$value <- factor(df$value, levels = unique(as.character(data)))
    }
    background_cells <- setdiff(rownames(df), active_cells)

    p <- ggplot() +
        geom_point(
            data = df[background_cells, ],
            aes(x = x, y = y),
            color = "lightgray",
            alpha = 0.4,
            size = pt_size - 1
        ) +
        ggnewscale::new_scale("color") +
        geom_point(
            data = df[active_cells, ],
            aes(x = x, y = y, color = value),
            alpha = 0.9,
            size = pt_size
        )

    if (is.numeric(data)) {
        if (is.null(colors)) {
            colors <- c("yellow", "blue")
        }
        p <- p +
            scale_color_gradientn(colors = colors, name = legend_name) +
            guides(color = guide_colorbar(order = 1))
    } else {
        value_levels <- levels(df$value)
        if (is.null(colors)) {
            colors <- scPalette(length(value_levels))
        }
        if (is.null(names(colors))) {
            names(colors) <- value_levels
        }
        p <- p +
            scale_color_manual(
                values = colors,
                limits = names(colors),
                name = legend_name,
                drop = FALSE
            ) +
            guides(color = guide_legend(order = 1))
    }

    p +
        coord_fixed() +
        # Theme
        theme_bw() +
        theme(text = element_text(size = 24))

}

#' Plot spatial neighborhoods
#'
#' Highlight selected cells and all cells within a distance threshold while
#' retaining cell-type colors and labels for the selected cells.
#'
#' @param meta A cell metadata data frame containing `x`, `y`, and optionally
#'   `annotation`, with cell identifiers as row names.
#' @param cell_names Character vector of focal cell identifiers.
#' @param dist_threshold Numeric maximum distance defining a neighborhood.
#' @param point_size Numeric size of background cells.
#' @param neighbor_size Numeric size of neighboring cells.
#' @param sample_size Numeric size of focal-cell markers.
#' @param label_size Numeric size of focal-cell labels.
#' @param all_alpha Numeric alpha for all background cells.
#' @param neighbor_alpha Numeric alpha for highlighted neighbors.
#' @param ncol Optional facet-column setting retained for compatibility.
#' @param cell_palette Optional named vector mapping annotations to colors.
#'
#' @return A `ggplot` object showing focal cells and their spatial neighborhoods.
PlotSpa_neighbor <- function(
    meta,
    cell_names,
    dist_threshold = 50,
    point_size = 0.6,
    neighbor_size = 1.4,
    sample_size = 6,
    label_size = 4,
    all_alpha = 0.25,
    neighbor_alpha = 0.95,
    ncol = NULL,
    cell_palette = NULL
) {
    if (!all(c("x", "y") %in% colnames(meta))) {
        stop("meta must contain x and y columns.")
    }
    if (is.null(rownames(meta))) {
        stop("meta must have cell names as rownames.")
    }
    cell_names <- as.character(cell_names)
    if (!length(cell_names)) {
        stop("cell_names must contain at least one cell name.")
    }

    missing_cells <- setdiff(cell_names, rownames(meta))
    if (length(missing_cells) > 0) {
        stop(
            "These cell_names are missing from meta rownames: ",
            paste(missing_cells, collapse = ", ")
        )
    }

    df <- meta
    df$cell <- rownames(df)
    df$x <- as.numeric(df$x)
    df$y <- as.numeric(df$y)
    if ("annotation" %in% colnames(df)) {
        df$annotation <- as.character(df$annotation)
    } else {
        df$annotation <- "Cell"
    }
    df$annotation[is.na(df$annotation) | !nzchar(df$annotation)] <- "Unknown"
    df <- df[!is.na(df$x) & !is.na(df$y), ]

    if (length(setdiff(cell_names, df$cell)) > 0) {
        stop("Some input cells have missing or non-numeric x/y coordinates.")
    }

    sample_info <- data.frame(
        cell = cell_names,
        sample_label = paste0("Sample #", seq_along(cell_names)),
        stringsAsFactors = FALSE
    )

    sample_df <- df[match(sample_info$cell, df$cell), ]
    sample_df$sample_label <- sample_info$sample_label

    neighbor_list <- lapply(seq_len(nrow(sample_info)), function(i) {
            center <- df[df$cell == sample_info$cell[i], c("x", "y")]
            dist <- sqrt((df$x - center$x) ^ 2 + (df$y - center$y) ^ 2)
            neighbors <- df[dist <= dist_threshold & df$cell != sample_info$cell[i], ]
            neighbors$distance <- dist[match(neighbors$cell, df$cell)]
            neighbors$sample_cell <- rep(sample_info$cell[i], nrow(neighbors))
            neighbors$sample_label <- rep(sample_info$sample_label[i], nrow(neighbors))
            neighbors
    })
    neighbor_df <- do.call(rbind, neighbor_list)
    if (is.null(neighbor_df)) {
        neighbor_df <- df[0, ]
        neighbor_df$distance <- numeric(0)
        neighbor_df$sample_cell <- character(0)
        neighbor_df$sample_label <- character(0)
    }
    if (nrow(neighbor_df) > 0) {
        neighbor_df <- neighbor_df[!duplicated(neighbor_df$cell), ]
    }

    annotations <- unique(df$annotation)
    if (is.null(cell_palette)) {
        cell_palette <- setNames(scPalette(length(annotations)), annotations)
    } else if (is.null(names(cell_palette))) {
        names(cell_palette) <- annotations[seq_along(cell_palette)]
    }
    missing_annotations <- setdiff(annotations, names(cell_palette))
    if (length(missing_annotations) > 0) {
        cell_palette <- c(
            cell_palette,
            setNames(scPalette(length(missing_annotations)), missing_annotations)
        )
    }

    ggplot() +
        geom_point(
            data = df,
            aes(x = x, y = y, color = annotation),
            size = point_size,
            alpha = all_alpha,
            show.legend = TRUE
        ) +
        geom_point(
            data = neighbor_df,
            aes(x = x, y = y, color = annotation),
            size = neighbor_size,
            alpha = neighbor_alpha,
            show.legend = FALSE
        ) +
        geom_point(
            data = sample_df,
            aes(x = x, y = y),
            size = sample_size,
            shape = 21,
            fill = "#D73027",
            color = "black",
            stroke = 1,
            show.legend = FALSE
        ) +
        geom_label(
            data = sample_df,
            aes(x = x, y = y, label = sample_label),
            size = label_size,
            linewidth = 0.2,
            fill = "white",
            color = "black",
            show.legend = FALSE
        ) +
        scale_color_manual(values = cell_palette, name = "Cell type") +
        coord_fixed() +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            text = element_text(family = "Arial", size = 18),
            legend.text = element_text(size = 14, family = "Arial")
        ) +
        labs(x = "x", y = "y")
}

#' Plot spatial communication vectors
#'
#' Convert a sender-by-receiver communication matrix into a spatial vector field
#' showing the dominant communication direction for each sender.
#'
#' @param data A numeric sender-by-receiver communication matrix.
#' @param meta A cell metadata data frame containing `x` and `y` coordinates.
#' @param threshold Numeric cap applied to normalized communication strength.
#' @param local Logical; whether to omit inactive background cells.
#' @param colors Optional vector of colors for the strength gradient.
#'
#' @return A `ggplot` object containing communication vectors.
PlotSpa_SRVector <- function(data, meta, threshold = 0.02, local = F, colors = NULL) {

    # data: Sender to Receiver matrix
    senders <- rownames(data)
    non_senders <- setdiff(rownames(meta), senders)
    data_used <- rowSums(data)
    sender_id <- max.col(data, ties.method = "first")
    receiver_of_sender <- colnames(data)[sender_id]

    df <- meta
    df$xend <- df$x
    df$yend <- df$y

    df$value <- 0
    df[senders, "value"] <- data_used
    df$value <- df$value / max(df$value)
    df$value <- pmin(df$value, threshold)

    df[senders, "xend"] <- meta[receiver_of_sender, "x"]
    df[senders, "yend"] <- meta[receiver_of_sender, "y"]

    df$dx <- df$xend - df$x
    df$dy <- df$yend - df$y
    df$xend_scaled <- df$x + df$value * df$dx
    df$yend_scaled <- df$y + df$value * df$dy
    sender_real <- rownames(df[df$value > 1e-4, ])
    non_senders <- setdiff(rownames(meta), sender_real)
    df <- df[order(df$value, decreasing = TRUE), ]
    if (is.null(colors)) {
        colors <- RColorBrewer::brewer.pal(5, "Spectral")
    }
    p <- ggplot() +
        geom_segment(
            data = df[sender_real, ],
            aes(x = x, y = y, xend = xend_scaled, yend = yend_scaled, color = value),
            arrow = arrow(length = unit(0.15, "cm")),
            lineend = "round", size = 0.7) +
        scale_color_gradientn(colors = colors, name = "Strength") +
        coord_fixed() +
        theme_void() +
        theme(text = element_text(family = "Arial", size = 24))
    if (!local) {
        p <- ggplot() +
            geom_point(aes(x, y), data = df[non_senders, ], color = "lightgray",
                alpha = 0.2, size = 2, show.legend = FALSE) +
            geom_segment(
                data = df[sender_real, ],
                aes(x = x, y = y, xend = xend_scaled, yend = yend_scaled, color = value),
                arrow = arrow(length = unit(0.15, "cm")),
                lineend = "round", size = 0.7) +
            scale_color_gradientn(colors = colors, name = "Strength") +
            coord_fixed() +
            theme_void() +
            theme(text = element_text(family = "Arial", size = 24)) +
            theme(text = element_text(family = "Arial", size = 24)) +
            labs(x = "x", y = "y")
    }

    return(p)
}

#' Plot sender-to-receiver communication circles
#'
#' Aggregate selected ligand-receptor signals by sender annotation and draw
#' spatially anchored arrows or self-loops toward a target cell type.
#'
#' @param data A sender-by-ligand-receptor communication matrix.
#' @param meta A cell metadata data frame containing `x`, `y`, and `annotation`.
#' @param features Character vector of ligand-receptor features to aggregate.
#' @param target_type Character scalar naming the receiving annotation.
#'
#' @return A `ggplot` object showing aggregated sender-to-target communication.
PlotSpa_SenderLRCircle <- function(data, meta, features, target_type) {

    # data: Sender * LR

    features <- intersect(features, colnames(data))
    if (!length(features)) {
        stop("No features found!")
    }

    all_types <- meta$annotation %>% unique()
    num_cluster <- length(all_types)
    node_coords <- matrix(0, nrow = num_cluster, ncol = 2)
    median_x <- aggregate(meta$x, list(meta$annotation), median)
    median_y <- aggregate(meta$y, list(meta$annotation), median)
    colnames(median_x) <- colnames(median_y) <- c("Cluster", "value")
    node_coords <- data.frame(x = median_x$value, y = median_y$value, Sender = median_x$Cluster)
    rownames(node_coords) <- median_x$Cluster

    color_use = scPalette(length(all_types))
    names(color_use) <- all_types

    data_used <- rowSums(data[, features])
    df <- meta
    df$color <- color_use[df$annotation]
    df$value <- 0
    senders <- rownames(data)
    df[senders, "value"] <- data_used
    aggr_df <- aggregate(df$value, list(df$annotation), sum)
    colnames(aggr_df) <- c("Sender", "Value")
    aggr_df$Sender <- as.character(aggr_df$Sender)
    aggr_df$xstart <- node_coords[aggr_df$Sender, "x"]
    aggr_df$ystart <- node_coords[aggr_df$Sender, "y"]
    aggr_df$xend <- node_coords[target_type, "x"]
    aggr_df$yend <- node_coords[target_type, "y"]
    aggr_df$Sender <- factor(aggr_df$Sender, levels = all_types)
    aggr_df$Receiver <- factor(target_type, levels = all_types)
    aggr_df$Value <- 3 * aggr_df$Value / max(aggr_df$Value) + 0.4
    aggr_df$Value[aggr_df$Value == 0.4] <- 0

    self_loops <- aggr_df[aggr_df$xstart == aggr_df$xend & aggr_df$ystart == aggr_df$yend, ]
    normal_edges <- aggr_df[!(aggr_df$xstart == aggr_df$xend & aggr_df$ystart == aggr_df$yend), ]

    diff_min <- min(abs(c(
                normal_edges$xstart - normal_edges$xend,
                normal_edges$ystart - normal_edges$yend
    )))
    self_loops$xstart <- self_loops$xstart - diff_min / 2

    ggplot() +
        geom_point(aes(x, y, color = annotation), data = df,
            alpha = 0.1, size = .25, show.legend = FALSE) +
        scale_color_manual(
            values = color_use,
            name = "Cell Types"
        ) + theme(text = element_text(family = "Arial", size = 24)) +
        labs(x = "x", y = "y") +
        geom_curve(data = normal_edges,
            aes(x = xstart, y = ystart, xend = xend, yend = yend),
            size = normal_edges$Value,
            curvature = 0.25,
            arrow = arrow(length = unit(0.15, "cm")),
            colour = color_use[normal_edges$Sender]) +
        ggforce::geom_arc(data = self_loops,
            aes(x0 = xstart, y0 = ystart, r = diff_min / 2,
                start = pi / 2, end = 5 * pi / 2), color = color_use[self_loops$Sender],
            size = self_loops$Value,
            inherit.aes = FALSE,
            arrow = arrow(length = unit(0.15, "cm"))) +
        geom_point(aes(x, y, color = Sender), data = node_coords, size = 4, show.legend = TRUE) +
        guides(color = guide_legend(override.aes = list(size = 18))) +
        scale_color_manual(values = color_use, name = "Cell type") +
        coord_fixed() +
        theme_bw()

}

#' Plot a crosstalk alluvial diagram
#'
#' Draw the strongest regulator-mediator-target relationships as an alluvial
#' diagram.
#'
#' @param results A four-column data frame containing regulator, mediator, target,
#'   and weight values.
#' @param features Optional character vector of targets to include.
#' @param topk Integer maximum number of relationships to display.
#' @param vertical Logical orientation setting passed to the plotting workflow.
#'
#' @return The object returned by `alluvial::alluvial()`.
PlotXT_Alluvial <- function(results, features = NULL, topk = 20, vertical = T) {
    require(alluvial)
    colnames(results) <- c("Regulator", "Mediator", "Target", "Weight")
    if (is.null(features)) {
        features <- unique(results$Target)
    }
    CC_used <- results[results$Target %in% features, ]
    CC_used <- CC_used[order(CC_used$Weight, decreasing = T), ]
    topk <- min(topk, nrow(CC_used))
    XT_toplot <- CC_used[1:topk, ]
    cols <- scPalette(length(uf))
    col_map <- setNames(cols, uf)
    feature_colors <- col_map[XT_toplot$Ligand]
    p <- alluvial::alluvial(XT_toplot[, 1:3], gap.width = 0.02, xw = 0.2, cw = 0.2,
        freq = XT_toplot$Weight, border = NA, alpha = 0.75,
        col = feature_colors)
    return(p)
}

#' Compare two pathways in space
#'
#' Calculate the relative activity of two extracellular-to-transcriptional paths
#' in each cell and display the resulting ratio spatially.
#'
#' @param XT_results A ligand-receptor-TF result data frame.
#' @param path1 Character vector of paths in the first set.
#' @param path2 Character vector of paths in the second set.
#' @param meta A cell metadata data frame containing spatial coordinates.
#' @param path_flag Character scalar specifying `"LR"`, `"RT"`, or `"LRT"` paths.
#' @param value_flag Character scalar selecting communication (`"CCC"`), regulatory
#'   (`"GRN"`), or combined (`"XT"`) weights.
#' @param sep Character separator used to construct path identifiers.
#' @param colors_grad Optional vector of colors for the ratio gradient.
#'
#' @return A `ggplot` object showing the pathway activity ratio in space.
PlotSpa_Comparison <- function(XT_results, path1, path2, meta,
    path_flag = "RT", value_flag = "XT", sep = "=",
    colors_grad = NULL) {

    # path_flag: LR, RT or LRT
    if (path_flag == "LR") {
        XT_results$RT <- paste0(XT_results$Ligand, sep, XT_results$Receptor)
    } else if (path_flag == "RT") {
        XT_results$RT <- paste0(XT_results$Receptor, sep, XT_results$TF)
    } else {
        XT_results$RT <- paste0(XT_results$Ligand, sep, XT_results$Receptor, sep, XT_results$TF)
    }
    XT_used1 <- XT_results[XT_results$RT %in% path1, ]
    XT_used1$class <- "Path1"
    XT_used2 <- XT_results[XT_results$RT %in% path2, ]
    XT_used2$class <- "Path2"
    XT_used <- rbind(XT_used1, XT_used2)
    if (value_flag == "CCC") {
        XT_used <- aggregate(XT_used$Value, list(XT_used$Cell, XT_used$class), sum)
    } else if (value_flag == "GRN") {
        XT_used <- aggregate(XT_used$beta, list(XT_used$Cell, XT_used$class), sum)
    } else {
        XT_used <- aggregate(XT_used$Weight, list(XT_used$Cell, XT_used$class), sum)
    }

    colnames(XT_used) <- c("Cell", "Path", "Weight")

    XT_used <- filter(XT_used, Weight > 1e-6)
    XT_ext <- XT_used %>%
        tidyr::pivot_wider(names_from = Path, values_from = Weight)
    XT_ext[is.na(XT_ext)] <- 0
    XT_ext$ratio <- 2 * XT_ext[["Path1"]] / (XT_ext[["Path1"]] + XT_ext[["Path2"]]) - 1
    rownames(XT_ext) <- XT_ext$Cell

    ratio <- XT_ext$ratio
    names(ratio) <- rownames(XT_ext)
    mycolor <- if (is.null(colors_grad)) {
        c("red", "yellow", "green", "lightblue", "purple")
    } else {
        colors_grad
    }
    PlotSpa_Scatter(
        data = ratio,
        meta = meta,
        colors = mycolor,
        pt_size = 2,
        legend_name = "Ratio",
        plot_zero = TRUE
    )

}

#' Compare two communication pathways in space
#'
#' Compare two columns of a cell-by-pathway communication matrix using either a
#' continuous activity ratio or categorical co-activation states.
#'
#' @param CCC_mat A cell-by-pathway communication matrix.
#' @param path1 Character scalar naming the first pathway column.
#' @param path2 Character scalar naming the second pathway column.
#' @param meta A cell metadata data frame containing spatial coordinates.
#' @param mode Character scalar; `"ratio"` for a continuous comparison or another
#'   value for categorical activation states.
#' @param pt_size Numeric point size.
#'
#' @return A `ggplot` object showing the pathway comparison in space.
PlotSpa_Comparison2 <- function(CCC_mat, path1, path2, meta, mode = "ratio", pt_size = 3) {
    x <- CCC_mat[, path1]
    y <- CCC_mat[, path2]
    exped <- names(x[x>0 | y>0])
    if (mode == "ratio") {
        values <- (2 * x / (x + y) - 1)[exped]
        legname <- paste0(path1, " vs \n ", path2, " ratio")
        legname <- gsub("=", "-", legname)
        PlotSpa_Scatter(
            data = values,
            meta = meta,
            colors = c("blue", "pink", "red"),
            pt_size = pt_size + 0.5,
            legend_name = legname,
            plot_zero = TRUE
        )

    } else {
        cases <- c("No activation", paste0(path1, " only"), paste0(path2, " only"), "Coexistance")
        cases <- gsub("=", "_", cases)
        values <- cases[((x > 0) + 2 * (y > 0)) + 1]
        names(values) <- names(x)
        colors <- c("lightgray", "#4C72B0", "#DD8452", "#55A868")
        names(colors) <- cases
        PlotSpa_Scatter(
            data = values,
            meta = meta,
            colors = colors,
            pt_size = pt_size,
            legend_name = NULL,
            plot_zero = TRUE
        )
    }

}

#' Plot a crosstalk heatmap
#'
#' Aggregate weighted extracellular-to-transcriptional relationships and display
#' the strongest rows and columns in an annotated heatmap.
#'
#' @param results A four-column regulator-mediator-target-weight data frame.
#' @param features Optional character vector of feature values to retain.
#' @param feature_flag Character scalar naming the feature column to filter.
#' @param topk Integer maximum number of rows and columns to display.
#'
#' @return A `ComplexHeatmap::Heatmap` object.
PlotXT_Heatmap <- function(results, features = NULL, feature_flag = "Target", topk = 999) {

    require(ComplexHeatmap)
    colnames(results) <- c("Regulator", "Mediator", "Target", "Weight")
    other_cols <- setdiff(colnames(results), feature_flag)
    if (!is.null(features)) {
        results_TG <- results[results[[feature_flag]] %in% features, ]
    } else {
        results_TG <- results
    }

    results_aggr <- aggregate(
        results_TG$Weight,
        list(results_TG[[other_cols[1]]], results_TG[[other_cols[2]]]),
        sum
    )

    colnames(results_aggr) <- c("Var1", "Var2", "Value")
    temp_mat <- tidyfst::df_mat(results_aggr, row = Var1, col = Var2, value = Value)
    legend_name <- "Fidelity"
    temp_mat[is.na(temp_mat)] <- 0
    temp_mat <- temp_mat / sum(temp_mat)
    temp_mat2 <- temp_mat[
        order(apply(temp_mat / sum(temp_mat), 1, mean), decreasing = T),
        order(apply(temp_mat / sum(temp_mat), 2, mean), decreasing = T)
    ]
    topk_used <- min(topk, ncol(temp_mat2))
    temp_mat2 <- temp_mat2[, 1:topk_used]
    topk_used <- min(topk, nrow(temp_mat2))
    temp_mat2 <- temp_mat2[1:topk_used, ]
    # row (receptor) annotation
    ra = rowAnnotation(Fid = anno_boxplot((temp_mat2 / sum(temp_mat2)), height = unit(4, "cm")))
    # column (TF) annotation
    ha = HeatmapAnnotation(
        Fid = anno_boxplot(
            temp_mat2 / sum(temp_mat2), which = 'column', height = unit(4, "cm")
        )
    )

    p <- ComplexHeatmap::Heatmap(temp_mat2, right_annotation = ra, top_annotation = ha,
        cluster_rows = T, cluster_columns = T, name = legend_name,
        row_names_gp = gpar(fontsize = 18),
        column_names_gp = gpar(fontsize = 18))
    return(p)
}

#' Plot dominant ligand-receptor signals
#'
#' Aggregate a communication result matrix by ligand, receptor, or complete
#' ligand-receptor pair and display the strongest signals as bars.
#'
#' @param results A communication matrix whose columns identify ligand-receptor
#'   pairs.
#' @param flag Character scalar specifying `"Ligand"`, `"Receptor"`, or `"LR"`.
#' @param topk Integer maximum number of signals to display.
#'
#' @return A `ggplot` bar chart of dominant signals.
Plot_DominantLR <- function(results, flag = "LR", topk = 20) {

    # results could either be CCC_receiver or CCC_sender

    Ligands <- sub("=.*", "", colnames(results))
    Receptors <- sub(".*=", "", colnames(results))
    if (flag == "Ligand") {
        CCC_sum_LR <- rowsum(t(results), group = Ligands) %>% rowSums()
        df <- data.frame(Signal = Ligands, Strength = CCC_sum_LR)
    }
    if (flag == "Receptor") {
        CCC_sum_LR <- rowsum(t(results), group = Receptors) %>% rowSums()
        df <- data.frame(Signal = Receptors, Strength = CCC_sum_LR)
    }
    if (flag == "LR") {
        CCC_sum_LR <- colSums(results)
        df <- data.frame(Signal = names(CCC_sum_LR), Strength = CCC_sum_LR)
    }

    topk <- min(topk, nrow(df))
    df <- df[order(df$Strength, decreasing = T), ]
    df <- df[1:topk, ]
    df <- df %>% mutate(Signal = fct_reorder(Signal, Strength))

    p <- ggplot(
        data = df,
        aes(
            x = Signal, # First variable on the X-axis
            y = Strength, # Second variable on the X-axis
            fill = Signal
    )) +
        coord_flip() +
        geom_bar(stat = "identity", color = "black") +
        theme_bw() +
        labs(y = "Total strength", x = "CCC signal") +
        theme(text = element_text(family = "Arial", size = 24)) +
        theme(legend.position = "none")

    return(p)
}

#' Plot spatial group compositions
#'
#' Place signal or target composition pies at spatial group centers while showing
#' the cells assigned to each group.
#'
#' @param cluster_results A named list of group results containing cells, signals,
#'   and targets.
#' @param meta A cell metadata data frame containing `x` and `y` coordinates.
#' @param pie_disp Character scalar selecting `"signal"` or target compositions.
#' @param topk Integer number of dominant components retained per group.
#' @param pie_r Numeric scatter-pie radius.
#' @param size_target Numeric point size for grouped cells.
#' @param size_others Numeric point size for other cells.
#' @param alpha_target Numeric alpha for grouped cells.
#' @param alpha_others Numeric alpha for other cells.
#'
#' @return A `ggplot` object with spatially positioned composition pies.
PlotSpa_Scatterpie <- function(cluster_results, meta, pie_disp = "signal", topk = 5,
    pie_r = 5, size_target = 3, size_others = 2, alpha_target = 0.6, alpha_others = 0.3) {

    pie_disp <- tolower(pie_disp)
    # Each group contains cells plus signal and target summary data frames.
    group_names <- names(cluster_results)
    group_color <- setNames(scPalette(length(group_names)), group_names)

    Cellid <- lapply(group_names, function(x) {
            temp <- cluster_results[[x]][["Cells"]]
            data.frame(Cell = temp, Type = x)
    })
    Cellid <- do.call(rbind, Cellid)
    Cellid$x <- meta[Cellid$Cell, "x"]
    Cellid$y <- meta[Cellid$Cell, "y"]
    other_cells <- setdiff(rownames(meta), Cellid$Cell)

    median_x <- aggregate(Cellid$x, list(Cellid$Type), median)
    median_y <- aggregate(Cellid$y, list(Cellid$Type), median)
    colnames(median_x) <- colnames(median_y) <- c("Group", "value")
    node_coords <- data.frame(Group = median_x$Group, x = median_x$value, y = median_y$value)

    Signal_mat <- lapply(group_names, function(x) {
            if (pie_disp == "signal") {
                nam <- cluster_results[[x]][["Signals"]][["Signal"]]
                temp <- cluster_results[[x]][["Signals"]][["Value"]]
                data.frame(Signal = nam, Value = temp, Group = x)
            } else {
                nam <- cluster_results[[x]][["Targets"]][["Target"]]
                temp <- cluster_results[[x]][["Targets"]][["Value"]]
                data.frame(Signal = nam, Value = temp, Group = x)
            }

    })
    Signal_mat <- do.call(rbind, Signal_mat)
    Signal_mat <- tidyfst::df_mat(Signal_mat, row = Group, col = Signal, value = Value)
    Signal_mat[is.na(Signal_mat)] <- 0
    allsignals <- colnames(Signal_mat)

    # Filter the signal mat
    topk <- min(topk, ncol(Signal_mat))
    topk_list <- lapply(split(Signal_mat, row(Signal_mat)), function(row) {
            vals <- sort(row, decreasing = TRUE)[1:topk]
            names(vals) <- colnames(Signal_mat)[order(row, decreasing = TRUE)[1:topk]]
            vals
    })
    df_topk <- do.call(rbind, lapply(seq_along(topk_list), function(i) {
                data.frame(Group = names(topk_list)[i],
                    Signal = names(topk_list[[i]]),
                    Value = as.numeric(topk_list[[i]]))
    }))
    Signal_mat_filtered <- tidyfst::df_mat(df_topk, row = Group, col = Signal, value = Value)
    Signal_mat_filtered[is.na(Signal_mat_filtered)] <- 0

    df_target <- cbind(node_coords, Signal_mat_filtered)
    signal_used <- colnames(Signal_mat_filtered)
    palette_size <- length(group_names) + length(signal_used)
    signal_color <- scPalette(palette_size)[
        (length(group_names) + 1):palette_size
    ]
    names(signal_color) <- signal_used

    p <- ggplot() +
        # Layer 1: Other types
        geom_point(
            data = meta[other_cells, ],
            aes(x = x, y = y),
            alpha = alpha_others,
            size = size_others,
            color = "lightgray"
        ) +
        ggnewscale::new_scale("color") +
        geom_point(
            data = Cellid,
            aes(x = x, y = y, color = Type),
            alpha = alpha_target,
            size = size_target
        ) +
        scale_color_manual(
            values = group_color,
            name = "Group"
        ) +
        ggnewscale::new_scale("color") +
        geom_scatterpie(data = df_target, aes(x = x, y = y, r = pie_r),
            cols = signal_used) +
        scale_fill_manual(
            values = signal_color,
            name = ifelse(pie_disp == "signal", "Signal", "Target")
        ) +
        coord_fixed() +
        # Theme
        theme_bw() +
        theme(text = element_text(size = 24)) +
        theme(axis.text.x = element_blank()) +
        theme(axis.text.y = element_blank())
    return(p)
}

#' Plot grouped circular bars
#'
#' Display the strongest feature-specificity values for multiple groups in a
#' circular bar chart.
#'
#' @param df A three-column data frame containing individual, group, and
#'   specificity values.
#' @param KeyFactors Optional character vector of groups to include.
#' @param topk Integer number of individuals retained per group.
#' @param label_max Optional numeric radial maximum used to position labels.
#'
#' @return A `ggplot` circular bar chart.
PlotXT_MultiCircularBar <- function(df, KeyFactors = NULL, topk = 5, label_max = NULL) {

    colnames(df) <- c("individual", "group", "Specificity")
    if (is.null(KeyFactors)) {
        KeyFactors <- df$group %>% as.array() %>% unique()
    } else {
        KeyFactors <- intersect(KeyFactors, df$group)
    }
    df <- dplyr::filter(df, group %in% KeyFactors)

    data <- c()
    for (type in KeyFactors) {
        temp_df <- dplyr::filter(df, group == type)
        temp_df <- temp_df[order(temp_df$Specificity, decreasing = T), ]
        data <- rbind(data, temp_df[1:topk, ])
    }
    # Create dataset
    colnames(data) <- c("individual", "group", "value")
    data$group <- as.factor(data$group)
    rownames(data) <- NULL

    empty_bar <- 3
    to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group), ncol(data)))
    colnames(to_add) <- colnames(data)
    to_add$group <- rep(levels(data$group), each = empty_bar)
    data <- rbind(data, to_add)
    data <- data %>% arrange(group)
    data$id <- seq(1, nrow(data))

    label_data <- data
    number_of_bar <- nrow(label_data)
    # Offset by 0.5 so labels use the center angle of each bar.
    angle <- 90 - 360 * (label_data$id - 0.5) / number_of_bar
    label_data$hjust <- ifelse(angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle + 180, angle)

    base_data <- data %>%
        group_by(group) %>%
        summarize(start = min(id), end = max(id) - empty_bar) %>%
        rowwise() %>%
        mutate(title = mean(c(start, end)))

    grid_data <- base_data
    grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
    grid_data$start <- grid_data$start - 1
    grid_data <- grid_data[-1, ]
    # windowsFonts(A = windowsFont("Arial"), T = windowsFont("Times New Roman"))

    if (is.null(label_max)) {
        label_max <- max(data$value, na.rm = T) * 1.25
        label_max <- signif(label_max, 2)
    }

    # Treating id as a factor avoids a gap before the first bar.
    p <- ggplot(data, aes(x = as.factor(id), y = value, fill = group)) +

        geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity") +
        geom_segment(
            data = base_data,
            aes(x = start, y = label_max * 0.8, xend = end, yend = label_max * 0.8),
            colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
        ) +
        geom_segment(
            data = base_data,
            aes(x = start, y = label_max * 0.6, xend = end, yend = label_max * 0.6),
            colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
        ) +
        geom_segment(
            data = base_data,
            aes(x = start, y = label_max * 0.4, xend = end, yend = label_max * 0.4),
            colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
        ) +
        geom_segment(
            data = base_data,
            aes(x = start, y = label_max * 0.2, xend = end, yend = label_max * 0.2),
            colour = "grey", alpha = 1, linewidth = 0.3, inherit.aes = FALSE
        ) +
        # Add text showing the value of each 100/75/50/25 lines
        annotate("text",
            x = rep(max(data$id), 4),
            y = c(label_max * 0.2, label_max * 0.4, label_max * 0.6, label_max * 0.8),
            label = paste0("Spe=", label_max * c(0.2, 0.4, 0.6, 0.8)),
            color = "black", size = 7.5, angle = 0, fontface = "bold", hjust = 1
        ) +
        geom_bar(aes(x = as.factor(id), y = value, fill = group), stat = "identity") +
        ylim(-0.8 * label_max, 1.1 * label_max) +
        theme_void() +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        ) +
        coord_polar() +
        geom_text(
            data = label_data,
            aes(x = id, y = label_max, label = individual, hjust = hjust),
            color = "black", fontface = "bold", alpha = 1, size = 12,
            angle = label_data$angle, inherit.aes = FALSE
        ) +
        # Add base line information
        geom_segment(
            data = base_data,
            aes(x = start, y = -0.15 * label_max, xend = end, yend = -0.15 * label_max),
            colour = "black", alpha = 1, size = 0.6, inherit.aes = FALSE
        ) +
        geom_text(
            data = base_data,
            aes(x = title, y = -0.5 * label_max, label = group),
            hjust = 0.5, colour = "black", alpha = 1, size = 12,
            fontface = "bold", inherit.aes = FALSE
        )
    return(p)
}

#' Plot source-target chords
#'
#' Draw a chord diagram from a source-by-target matrix and label each sector.
#'
#' @param mat A numeric source-by-target matrix.
#' @param orders Optional character vector defining sector order.
#' @param edge_colors Optional vector or matrix of chord colors.
#'
#' @return Invisibly returns the result of the final `circlize` drawing call and
#'   draws the chord diagram on the active graphics device.
PlotXT_Chord <- function(mat, orders = NULL, edge_colors = NULL) {

    if (is.null(orders)) {
        orders <- c(rownames(mat), colnames(mat))
    }
    circos.clear()
    if (is.null(edge_colors)) {
        chordDiagram(t(mat),
            transparency = 0.25,
            order = orders,
            big.gap = 30,
            annotationTrack = "grid", # Show sector labels
            annotationTrackHeight = 0.05,
            preAllocateTracks = 1    # Reserve space for labels
        )
    } else {
        chordDiagram(t(mat),
            col = edge_colors,
            transparency = 0.25,
            order = orders,
            big.gap = 30,
            annotationTrack = "grid", # Show sector labels
            annotationTrackHeight = 0.05,
            preAllocateTracks = 1    # Reserve space for labels
        )
    }

    circos.track(track.index = 1, panel.fun = function(x, y) {
            xlim = get.cell.meta.data("xlim")
            xplot = get.cell.meta.data("xplot")
            ylim = get.cell.meta.data("ylim")
            sector.name = get.cell.meta.data("sector.index")

            if (abs(xplot[2] - xplot[1]) < 15) {
                circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise",
                    niceFacing = TRUE, adj = c(0, 0.5), col = "blue", cex = 2.5)
            } else {
                circos.text(mean(xlim), ylim[1], sector.name, facing = "inside",
                    niceFacing = TRUE, adj = c(0.5, 0), col = "blue", cex = 2.5)
            }
        }, bg.border = NA)
}

#' Plot a dominant-mediator bubble heatmap
#'
#' Summarize high-weight extracellular-to-transcriptional relationships and show
#' their total strength and dominant mediator in a bubble heatmap.
#'
#' @param results A four-column data frame containing two endpoint variables, a
#'   mediator, and a weight.
#' @param thresh_quantile Numeric quantile used to filter weak relationships.
#' @param topk_Rec Integer number of dominant mediators assigned distinct colors.
#'
#' @return A `ggplot` bubble heatmap.
Plot_BubbleHeatmap <- function(results, thresh_quantile = 0.95, topk_Rec = 5) {

    # result: XT_result aggregated by cells
    colnames(results) <- c("varx", "vary", "varz", "value")
    cutoff <- quantile(results$value, thresh_quantile, na.rm = T)
    results_used <- results[results$value >= cutoff, ]

    df_pair <- results_used %>%
        group_by(varx, varz, vary) %>%
        summarise(v = sum(value), .groups = "drop")

    df_bubble <- df_pair %>%
        group_by(varx, varz) %>%
        summarise(
            total_value = sum(v),
            winner_vary = vary[which.max(v)],
            .groups = "drop"
        )

    top_vary <- df_bubble %>%
        count(winner_vary, sort = TRUE) %>%
        slice_head(n = topk_Rec) %>%
        pull(winner_vary)

    varx_all <- sort(df_bubble$varx %>% unique())
    df_bubble2 <- df_bubble %>%
        mutate(
            winner_vary_collapsed = ifelse(winner_vary %in% top_vary, winner_vary, "Others"),
            winner_vary_collapsed = factor(winner_vary_collapsed, levels = c(top_vary, "Others")),
            varx = factor(varx, levels = rev(varx_all))
        )

    p <- ggplot(df_bubble2, aes(x = varz, y = varx)) +
        geom_point(aes(size = total_value, color = winner_vary_collapsed), alpha = 0.8) +
        scale_size_continuous(range = c(3, 10)) +
        theme_bw() +
        theme(text = element_text(size = 18)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(legend.text = element_text(size = 18, family = "Arial")) +
        guides(color = guide_legend(override.aes = list(size = 7))) +
        labs(
            x = "Target",
            y = "Signal",
            size = "Total Strength",
            color = "Dominant Receptor"
        )

    return(p)
}

#' Plot a spatial activity surface
#'
#' Aggregate selected result columns, optionally log-transform their values, and
#' display the resulting spatial activity as a raster surface.
#'
#' @param results A cell-by-feature numeric matrix or data frame.
#' @param meta A cell metadata data frame containing `x` and `y` coordinates.
#' @param target Optional character vector of result columns to aggregate.
#' @param logscale Logical; whether to apply `log1p()` to activity values.
#' @param nlines Integer contour-resolution setting retained for compatibility.
#'
#' @return A `ggplot` object showing the spatial activity surface.
PlotSpa_Contour <- function(results, meta, target = NULL, logscale = T, nlines = 10) {

    if (is.null(target)) {
        target <- colnames(results)
    } else {
        target <- intersect(target, colnames(results))
    }
    results_filtered <- results[, target] %>% as.data.frame()
    if (ncol(results_filtered) > 1) {
        results_filtered <- rowSums(results_filtered)
    }
    rownames(results_filtered) <- rownames(results)
    df <- meta
    df$Value <- 0
    df[rownames(results_filtered), "Value"] <- results_filtered[, 1]
    if (logscale) {
        df$Value <- log1p(df$Value)
    }
    df$Value[is.na(df$Value)] <- 0
    df$Value[is.infinite(df$Value)] <- 0

    colormap <- rev(RColorBrewer::brewer.pal(n = nlines + 1, name = "Spectral"))
    breaks_lines <- seq(min(df$Value), max(df$Value), by = (max(df$Value) - min(df$Value)) / nlines)
    p <- ggplot(data = df) +
        geom_raster(aes(x = x, y = y, fill = Value)) +
        scale_fill_viridis_c()
    # scale_color_gradientn(colours=colormap)+
    # geom_contour()+
    theme_bw()+
        theme(text = element_text(family = "Arial", size = 24))+
        theme(panel.grid.major.x = element_line(color = "grey85"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            axis.title = element_blank(),
            legend.position = "top",
            text = element_text(family = "Arial", size = 24))
    return(p)
}

#' Plot single-cell multilayer crosstalk
#'
#' Draw ligand, receptor, and TF nodes in separate layers with weighted arrows
#' connecting adjacent signaling layers.
#'
#' @param XT_filtered A ligand-receptor-TF data frame containing `Value` and
#'   `beta` weights.
#' @param x_scale Numeric horizontal spacing within each layer.
#' @param y_gap Numeric vertical spacing between layers.
#' @param rect_width Numeric node-box width.
#' @param rect_height Numeric node-box height.
#' @param start_trim Numeric distance trimmed from arrow starts.
#' @param end_trim Numeric distance trimmed from arrow ends.
#' @param arrow_len_cm Numeric arrowhead length in centimeters.
#' @param abs_value Logical; whether edge values should be treated as absolute.
#' @param layer_fill_colors Optional named vector of layer fill colors.
#' @param border_darken_factor Numeric factor used to darken node borders.
#' @param edge_width_range Numeric length-two range for edge widths.
#' @param box_linewidth Numeric node-box border width.
#' @param label_size Numeric node-label size.
#' @param label_face Character font face for node labels.
#' @param show_layer_legend Logical; whether to show the layer fill legend.
#'
#' @return A `ggplot` object showing multilayer signaling paths.
Plot_MultiLayer_Singlecell <- function(
    XT_filtered,
    x_scale = 1.0, # horizontal spacing within each layer
    y_gap = 1.0, # vertical spacing between layers
    rect_width = 0.9,
    rect_height = 0.35,
    start_trim = 0.1, # shorten arrow at source (in y units)
    end_trim = 0.1, # shorten arrow at target (in y units)
    arrow_len_cm = 0.15,
# whether use absolute value
    abs_value = False,
# styling
    layer_fill_colors = NULL, # named vector with names matching layer_levels; if NULL, auto
    border_darken_factor = 0.35, # 0..1, larger => darker border
    edge_width_range = c(0.4, 2.5),
    box_linewidth = 0.8,
    label_size = 6,
    label_face = "bold",
# legends
    show_layer_legend = TRUE      # if FALSE, hide fill legend too
) {

    layer_levels = c("Target", "Receptor", "Ligand")
    XT_filtered$beta <- abs(XT_filtered$beta)
    edge1 <- data.frame(from_name = XT_filtered$Ligand,
        from_layer = "Ligand",
        to_name = XT_filtered$Receptor,
        to_layer = "Receptor",
        weight_abs = XT_filtered$Value,
        weight = XT_filtered$Value / sum(XT_filtered$Value))
    edge1 <- edge1[!duplicated(edge1), ]
    edge2 <- data.frame(from_name = XT_filtered$Receptor,
        from_layer = "Receptor",
        to_name = XT_filtered$TF,
        to_layer = "Target",
        weight_abs = XT_filtered$beta,
        weight = XT_filtered$beta / sum(XT_filtered$beta))
    edge2 <- edge2[!duplicated(edge2), ]
    edges <- rbind(edge1, edge2)

    nodes <- bind_rows(
        edges %>%
            transmute(
                !!"node_name" := .data[["from_name"]],
                !!"layer" := .data[["from_layer"]]
            ),
        edges %>%
            transmute(
                !!"node_name" := .data[["to_name"]],
                !!"layer" := .data[["to_layer"]]
            )
    ) %>% distinct()
    # ---- 1) Standardize nodes, make layer-aware node_id, and compute positions ----
    nodes_std <- nodes %>%
        transmute(
            node_name = .data[["node_name"]],
            layer_raw = .data[["layer"]],
            value = 1
        ) %>%
        mutate(
            layer = factor(layer_raw, levels = layer_levels),
            node_id = paste(node_name, as.character(layer), sep = "__")
        ) %>%
        select(node_name, layer, node_id, value)

    nodes_pos <- nodes_std %>%
        group_by(layer) %>%
        arrange(node_name, .by_group = TRUE) %>%
        mutate(
            x = (seq_len(n()) - (n() + 1) / 2) * x_scale,
            y = as.numeric(layer) * y_gap
        ) %>%
        ungroup()

    # rectangle coords
    nodes_rect <- nodes_pos %>%
        mutate(
            xmin = x - rect_width / 2,
            xmax = x + rect_width / 2,
            ymin = y - rect_height / 2,
            ymax = y + rect_height / 2
        )

    # ---- 2) Standardize edges, drop non-adjacent edges, join coordinates ----
    edges_std <- edges %>%
        transmute(
            from_name = .data[["from_name"]],
            from_layer_raw = .data[["from_layer"]],
            to_name = .data[["to_name"]],
            to_layer_raw = .data[["to_layer"]],
            weight_abs = .data[["weight_abs"]],
            weight = as.numeric(.data[["weight"]])
        ) %>%
        mutate(
            from_layer = factor(from_layer_raw, levels = layer_levels),
            to_layer = factor(to_layer_raw, levels = layer_levels),
            from_i = as.numeric(from_layer),
            to_i = as.numeric(to_layer),
            from_id = paste(from_name, from_layer, sep = "__"),
            to_id = paste(to_name, to_layer, sep = "__")
        ) %>%
        # keep ONLY adjacent layers (top<->middle, middle<->bottom)
        filter(abs(from_i - to_i) == 1)

    edges_pos <- edges_std %>%
        left_join(nodes_pos %>% select(node_id, x_from = x, y_from = y),
            by = c("from_id" = "node_id")) %>%
        left_join(nodes_pos %>% select(node_id, x_to = x, y_to = y),
            by = c("to_id" = "node_id"))

    # sanity: drop edges that failed to join (missing nodes)
    edges_pos <- edges_pos %>%
        filter(!is.na(x_from), !is.na(y_from), !is.na(x_to), !is.na(y_to))

    # ---- 3) Build shortened equal-height arrow coordinates ----
    edges_arrow <- edges_pos %>%
        mutate(
            going_down = y_from > y_to,
            y_start = ifelse(going_down, y_from - start_trim, y_from + start_trim),
            y_end = ifelse(going_down, y_to + end_trim, y_to - end_trim),
            x_start = x_from,
            x_end = x_to
        )

    # ---- 4) Colors (fill + darker border) ----
    if (is.null(layer_fill_colors)) {
        # auto colors if not provided (stable order)
        auto_cols <- scales::hue_pal()(length(layer_levels))
        layer_fill_colors <- setNames(auto_cols, layer_levels)
    } else {
        # ensure named
        if (is.null(names(layer_fill_colors))) {
            stop(
                "layer_fill_colors must be a named vector, ",
                "e.g. c(top='#...', middle='#...', bottom='#...')."
            )
        }
    }
    layer_border_colors <- colorspace::darken(layer_fill_colors, amount = border_darken_factor)

    # ---- 5) Plot ----
    p <- ggplot() +
        # arrows: width encodes weight, color fixed
        geom_segment(
            data = edges_arrow,
            aes(x = x_start, y = y_start, xend = x_end, yend = y_end, linewidth = weight),
            arrow = arrow(length = unit(arrow_len_cm, "cm"), type = "closed"),
            lineend = "round"
        ) +
        # node rectangles: fill=layer, border=layer (darker palette)
        geom_rect(
            data = nodes_rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = layer, color = layer),
            linewidth = box_linewidth
        ) +
        # labels inside rectangles
        geom_text(
            data = nodes_pos,
            aes(x = x, y = y, label = node_name),
            size = label_size,
            fontface = label_face
        ) +
        scale_fill_manual(values = layer_fill_colors, name = "Layer") +
        scale_color_manual(values = layer_border_colors, guide = "none") +  # hide border legend
    scale_linewidth(range = edge_width_range, name = "Fidelity/Specificity") +
        scale_y_continuous(
            breaks = sort(unique(nodes_pos$y)),
            labels = layer_levels,
            expand = expansion(mult = 0.35)
        ) +
        labs(x = NULL, y = NULL) +
        coord_cartesian() +
        theme_minimal(base_size = 13) +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = 24),
            legend.text = element_text(size = 24)
        )

    if (!isTRUE(show_layer_legend)) {
        p <- p + guides(fill = "none")
    }

    p
}

#' Plot radial signaling around a receiver
#'
#' Draw sender cells and a receiver at their spatial coordinates with offset,
#' strength-scaled arrows for each active ligand-receptor signal.
#'
#' @param CCC_mat A ligand-receptor-by-sender communication matrix.
#' @param meta A cell metadata data frame containing `x`, `y`, and `annotation`.
#' @param receiver_name Optional receiver cell identifier; inferred when possible.
#' @param node_size Numeric node size.
#' @param offset_step Numeric relative spacing between parallel edges.
#' @param edge_trim Numeric relative amount trimmed from edge endpoints.
#' @param arrow_len_mm Numeric arrowhead length in millimeters.
#' @param width_range Numeric length-two range for edge widths.
#' @param node_palette Optional named vector mapping cell types to colors.
#' @param no_communication_alpha Numeric alpha for inactive sender cells.
#'
#' @return A `ggplot` object showing local signaling around the receiver.
Plot_Signaling_Radial <- function(
    CCC_mat,
    meta,
    receiver_name = NULL,
    node_size = 6,
    offset_step = 0.01,
    edge_trim = 0.1,
    arrow_len_mm = 3.5,
    width_range = c(0.8, 3.5),
    node_palette = NULL,
    no_communication_alpha = 0.45
) {
    if (!all(c("x", "y", "annotation") %in% colnames(meta))) {
        stop("meta must contain x, y, and annotation columns.")
    }

    CCC_mat <- as.matrix(CCC_mat)
    if (is.null(rownames(CCC_mat)) || is.null(colnames(CCC_mat))) {
        stop("CCC_mat must have ligand=receptor rownames and sender=receiver colnames.")
    }

    sender_receiver <- colnames(CCC_mat)
    has_receiver <- grepl("=", sender_receiver, fixed = TRUE)
    senders <- ifelse(has_receiver, sub("=.*", "", sender_receiver), sender_receiver)
    parsed_receivers <- ifelse(has_receiver, sub(".*=", "", sender_receiver), NA_character_)

    if (is.null(receiver_name)) {
        receiver_name <- unique(parsed_receivers[
                !is.na(parsed_receivers) & nzchar(parsed_receivers)
        ])
        if (length(receiver_name) != 1) {
            stop(
                "receiver_name must be supplied when CCC_mat colnames do not ",
                "contain one shared receiver."
            )
        }
    }
    receiver_name <- as.character(receiver_name[1])

    node_names <- unique(c(senders, receiver_name))
    missing_nodes <- setdiff(node_names, rownames(meta))
    if (length(missing_nodes) > 0) {
        stop(
            "These sender/receiver cells are missing from meta rownames: ",
            paste(missing_nodes, collapse = ", ")
        )
    }

    signals <- as.data.frame(as.table(CCC_mat), stringsAsFactors = FALSE)
    colnames(signals) <- c("signal_type", "sender_receiver", "strength")
    signals$sender <- senders[match(signals$sender_receiver, sender_receiver)]
    signals$receiver <- receiver_name

    coords <- meta[node_names, c("x", "y")]
    coords$node <- rownames(coords)
    # --- standardize coords ---
    coords_std <- coords %>%
        transmute(
            node = .data[["node"]],
            x = as.numeric(.data[["x"]]),
            y = as.numeric(.data[["y"]])
        ) %>%
        filter(!is.na(node), !is.na(x), !is.na(y))
    types_std <- data.frame(
        node = node_names,
        cell_type = as.character(meta[node_names, "annotation"])
    ) %>%
        transmute(
            node = .data[["node"]],
            cell_type = .data[["cell_type"]]
        )
    types_std$cell_type[is.na(types_std$cell_type) | !nzchar(types_std$cell_type)] <- "Unknown"

    edges <- signals %>%
        transmute(
            sender = .data[["sender"]],
            receiver = .data[["receiver"]],
            signal_type = .data[["signal_type"]],
            strength = as.numeric(.data[["strength"]])
        ) %>%
        filter(!is.na(strength), is.finite(strength), strength > 0) %>%
        left_join(coords_std %>% rename(sender = node, x_s = x, y_s = y), by = "sender") %>%
        left_join(coords_std %>% rename(receiver = node, x_r = x, y_r = y), by = "receiver") %>%
        filter(!is.na(x_s), !is.na(y_s), !is.na(x_r), !is.na(y_r))

    diag_len <- sqrt(diff(range(coords_std$x)) ^ 2 + diff(range(coords_std$y)) ^ 2)
    if (!is.finite(diag_len) || diag_len == 0) {
        diag_len <- 1
    }

    offset_step <- offset_step * diag_len
    trim_start <- edge_trim * diag_len
    trim_end <- edge_trim * diag_len

    if (nrow(edges) > 0) {
        edges2 <- edges %>%
            group_by(sender, receiver) %>%
            mutate(
                k = n(),
                i = row_number(),
                offset_index = i - (k + 1) / 2,
                offset_dist = offset_index * offset_step
            ) %>%
            ungroup() %>%
            rowwise() %>%
            mutate(
                dx = x_r - x_s,
                dy = y_r - y_s,
                len = sqrt(dx ^ 2 + dy ^ 2),
                nx = if (len > 0) - dy / len else 0,
                ny = if (len > 0)  dx / len else 0,
                x_s_off = x_s + nx * offset_dist,
                y_s_off = y_s + ny * offset_dist,
                x_r_off = x_r + nx * offset_dist,
                y_r_off = y_r + ny * offset_dist
            ) %>%
            ungroup()

        ## ---- trim arrows ----
        edges3 <- edges2 %>%
            rowwise() %>%
            mutate(
                dx2 = x_r_off - x_s_off,
                dy2 = y_r_off - y_s_off,
                L = sqrt(dx2 ^ 2 + dy2 ^ 2),
                ux = if (L > 0) dx2 / L else 0,
                uy = if (L > 0) dy2 / L else 0,
                trim_total = trim_start + trim_end,
                trim_scale = if (L > 0 && trim_total >= L) 0.9 * L / trim_total else 1,
                trim_start_i = trim_start * trim_scale,
                trim_end_i = trim_end * trim_scale,
                x_start = x_s_off + trim_start_i * ux,
                y_start = y_s_off + trim_start_i * uy,
                x_end = x_r_off - trim_end_i * ux,
                y_end = y_r_off - trim_end_i * uy
            ) %>%
            ungroup()
    } else {
        edges3 <- edges %>%
            mutate(
                x_start = numeric(0),
                y_start = numeric(0),
                x_end = numeric(0),
                y_end = numeric(0)
            )
    }
    edges3$signal_type <- gsub("=", " - ", edges3$signal_type)

    active_senders <- unique(edges$sender)
    nodes_plot <- coords_std %>%
        left_join(types_std, by = "node") %>%
        mutate(
            role = ifelse(.data[["node"]] == receiver_name, "receiver", "sender"),
            has_signal = .data[["node"]] == receiver_name | .data[["node"]] %in% active_senders,
            node_group = ifelse(
                .data[["role"]] == "sender" & !.data[["has_signal"]],
                "No communication",
                .data[["cell_type"]]
            ),
            node_alpha = ifelse(
                .data[["node_group"]] == "No communication",
                no_communication_alpha,
                1
            )
        )

    ## ---- color palette for node types ----
    if (is.null(node_palette)) {
        node_groups <- unique(nodes_plot$node_group)
        cell_groups <- node_groups[node_groups != "No communication"]
        node_palette <- character(0)
        if (length(cell_groups) > 0) {
            node_palette <- setNames(scPalette(length(cell_groups)), cell_groups)
        }
        if ("No communication" %in% node_groups) {
            node_palette <- c(node_palette, "No communication" = "grey70")
        }
    } else {
        if (is.null(names(node_palette))) {
            names(node_palette) <- unique(nodes_plot$cell_type)[seq_along(node_palette)]
        }
        missing_groups <- setdiff(unique(nodes_plot$node_group), names(node_palette))
        if ("No communication" %in% missing_groups) {
            node_palette <- c(node_palette, "No communication" = "grey70")
            missing_groups <- setdiff(missing_groups, "No communication")
        }
        if (length(missing_groups) > 0) {
            node_palette <- c(
                node_palette,
                setNames(scPalette(length(missing_groups)), missing_groups)
            )
        }
    }

    ggplot() +
        geom_point(
            data = nodes_plot,
            aes(x = x, y = y, fill = node_group, alpha = node_alpha),
            size = node_size,
            shape = 21,
            color = "black",
            stroke = 1
        ) +
        geom_segment(
            data = edges3,
            aes(
                x = x_start, y = y_start,
                xend = x_end, yend = y_end,
                color = signal_type,
                linewidth = strength
            ),
            arrow = arrow(length = unit(arrow_len_mm, "mm"), type = "closed"),
            lineend = "round",
            alpha = 0.9
        ) +
        scale_fill_manual(values = node_palette, name = "Sender cell type") +
        scale_alpha_identity() +
        scale_linewidth(range = width_range, name = "CCC Strength") +
        coord_equal() +
        theme_void(base_size = 14) +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(family = "Arial", size = 18),
            legend.text = element_text(size = 18, family = "Arial")
        ) +
        labs(
            color = "Ligand - Receptor",
            x = NULL,
            y = NULL
        )
}

#' Plot a two-sample GRN volcano comparison
#'
#' Compare gene-level regulatory specificity between two samples, calculate
#' per-gene significance, and display effect size against adjusted significance.
#'
#' @param df1 A data frame containing `Cell`, `Gene`, and `beta` for sample one.
#' @param df2 A data frame containing `Cell`, `Gene`, and `beta` for sample two.
#' @param min_cells_per_group Integer minimum observations required per sample.
#' @param pseudocount Numeric pseudocount used for log-fold changes.
#' @param test Character scalar selecting `"wilcox"` or a t-test fallback.
#' @param sig_lfc_cut Numeric absolute log2-fold-change significance threshold.
#' @param sig_p_cut Numeric unadjusted p-value significance threshold.
#' @param label_lfc_cut Numeric absolute log2-fold-change labeling threshold.
#' @param label_fdr_cut Numeric adjusted p-value labeling threshold.
#' @param cols Named color vector for change categories.
#' @param label_size Numeric text size for gene labels.
#' @param point_size Numeric point size.
#'
#' @return A `ggplot` volcano plot.
Plot_GRNVolcano <- function(
    df1, df2,
    min_cells_per_group = 10,
    pseudocount = 1e-6,
    test = "wilcox",
    sig_lfc_cut = 1,
    sig_p_cut = 0.05,
    label_lfc_cut = 5,
    label_fdr_cut = 0.01,
    cols = c(
        "Up-regulated" = "#D73027",
        "Down-regulated" = "#4575B4",
        "No change" = "grey70"
    ),
    label_size = 6,
    point_size = 1.5
) {
    df1 <- df1[, c("Cell", "Gene", "beta")]
    colnames(df1) <- c("Cell", "Gene", "Value")
    df2 <- df2[, c("Cell", "Gene", "beta")]
    colnames(df2) <- c("Cell", "Gene", "Value")
    df1$Value <- abs(df1$Value)
    df2$Value <- abs(df2$Value)
    df1 <- df1 %>%
        group_by(Cell) %>%
        mutate(Spe = Value / sum(Value)) %>%
        ungroup()
    df2 <- df2 %>%
        group_by(Cell) %>%
        mutate(Spe = Value / sum(Value)) %>%
        ungroup()
    x <- bind_rows(
        df1 %>% mutate(Sample = "S1"),
        df2 %>% mutate(Sample = "S2")
    ) %>%
        mutate(
            Cell = as.character(Cell),
            Gene = as.character(Gene),
            Value = as.numeric(Value)
        )
    res <- x %>%
        group_by(Gene) %>%
        group_modify(~{
                d <- .x
                v1 <- d$Value[d$Sample == "S1"]
                v2 <- d$Value[d$Sample == "S2"]
                n1 <- sum(!is.na(v1))
                n2 <- sum(!is.na(v2))

                mean1 <- mean(v1, na.rm = TRUE)
                mean2 <- mean(v2, na.rm = TRUE)

                if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
                    return(tibble(
                            n1 = n1,
                            n2 = n2,
                            mean1 = mean1,
                            mean2 = mean2,
                            log2FC = NA_real_,
                            p_value = NA_real_
                    ))
                }
                log2FC <- log2((mean2 + pseudocount) / (mean1 + pseudocount))
                p <- tryCatch({
                    if (test == "wilcox") {
                        wilcox.test(v2, v1)$p.value
                    } else {
                        t.test(v2, v1)$p.value
                    }
                }, error = function(e) NA_real_)
                tibble(
                    n1 = n1,
                    n2 = n2,
                    mean1 = mean1,
                    mean2 = mean2,
                    log2FC = log2FC,
                    p_value = p
                )
        }) %>%
        ungroup() %>%
        mutate(
            padj = p.adjust(p_value, method = "BH"),
            neglog10_p = -log10(p_value),
            neglog10_padj = -log10(padj)
        )
    res$Gene <- gsub("=", "-", res$Gene)
    data <- res %>%
        mutate(
            change = case_when(
                p_value < sig_p_cut & log2FC >  sig_lfc_cut ~ "Up-regulated",
                p_value < sig_p_cut & log2FC < -sig_lfc_cut ~ "Down-regulated",
                TRUE ~ "No change"
            ),
            change = factor(
                change,
                levels = c("Up-regulated", "Down-regulated", "No change")
            )
        )
    label_data <- data %>%
        filter(
            !is.na(log2FC),
            !is.na(padj),
            abs(log2FC) >= label_lfc_cut,
            padj <= label_fdr_cut
        )
    ggplot(data, aes(x = log2FC, y = -log10(padj))) +
        geom_hline(
            yintercept = -log10(sig_p_cut),
            linetype = "dashed",
            color = "#999999"
        ) +
        geom_vline(
            xintercept = c(-sig_lfc_cut, sig_lfc_cut),
            linetype = "dashed",
            color = "#999999"
        ) +
        geom_point(
            aes(color = change),
            size = point_size
        ) +
        geom_text_repel(
            data = label_data,
            aes(label = Gene, color = change),
            size = label_size,
            max.overlaps = getOption("ggrepel.max.overlaps", default = 20)
        ) +
        scale_color_manual(values = cols) +
        theme_bw(base_size = 12) +
        theme(
            panel.grid = element_blank(),
            legend.position = "top",
            text = element_text(family = "Arial", size = 18),
            legend.text = element_text(size = 18, family = "Arial")
        ) +
        xlab("Log2FC") +
        ylab("-Log10(FDR q-value)") +
        labs(color = "")
}

#' Plot bundled source-target edges
#'
#' Arrange source and target nodes on a circle and draw value-colored bundled
#' edges between them.
#'
#' @param df A three-column source-target-value data frame.
#' @param source_color Color used for source nodes.
#' @param target_color Color used for target nodes.
#' @param edge_colors Vector of colors for the edge-value gradient.
#' @param node_size Numeric node size.
#' @param text_size Numeric node-label size.
#' @param text_offset Numeric radial multiplier for label positions.
#' @param edge_width Numeric edge width.
#' @param edge_alpha Numeric edge alpha.
#' @param edge_strength Numeric arc-bundling strength.
#' @param legend_position Legend position passed to the plot theme.
#'
#' @return A `ggraph` plot object.
Plot_EdgeBundle <- function(df, source_color = "#0072B2",
    target_color = "#D55E00",
    edge_colors = c("lightgray", "yellow", "red"),
    node_size = 4,
    text_size = 3,
    text_offset = 1.12,
    edge_width = 0.7,
    edge_alpha = 0.7,
    edge_strength = 0.6,
    legend_position = "right") {

    colnames(df) <- c("From", "To", "Value")
    sources <- sort(unique(df[["From"]]))
    targets <- sort(setdiff(unique(df[["To"]]), sources))

    node_order <- c(sources, targets)

    nodes <- data.frame(name = node_order) %>%
        mutate(type = ifelse(name %in% sources, "source", "target"))

    edges <- df %>%
        transmute(
            from = .data[["From"]],
            to = .data[["To"]],
            Value = .data[["Value"]]
        )
    g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)
    layout_df <- create_layout(g, layout = "linear", circular = TRUE)

    ggraph(layout_df) +
        geom_edge_arc(
            aes(edge_colour = Value),
            alpha = edge_alpha,
            strength = edge_strength,
            width = edge_width
        ) +
        geom_node_point(aes(color = type), size = node_size) +
        geom_node_text(
            aes(label = name, x = text_offset * x, y = text_offset * y),
            size = text_size
        ) +
        scale_edge_colour_gradientn(
            colors = edge_colors,
            name = "Value"
        ) +
        scale_color_manual(
            values = c(source = source_color, target = target_color),
            name = "Node type"
        ) +
        coord_fixed() +
        theme_void() +
        theme(
            legend.position = legend_position
        )
}

#' Plot signaling fidelity or specificity
#'
#' For a selected ligand or TF, rank associated paths and display normalized
#' weights as a tiled heatmap.
#'
#' @param Aggr_XT_results An aggregated ligand-receptor-TF result data frame.
#' @param feature Character scalar naming the ligand or TF of interest.
#' @param mode Character scalar; `"Spe"` plots specificity and other values plot
#'   fidelity.
#' @param topk Integer maximum number of paths to display.
#'
#' @return A `ggplot` tile heatmap.
Plot_FidSpe <- function(Aggr_XT_results, feature, mode = "Spe", topk = 20) {

    if (mode == "Spe") {
        result <- Aggr_XT_results[Aggr_XT_results$Ligand == feature, ]
        result <- result[, c("Receptor", "TF", "Weight")] %>%
            mutate(Specificity = Weight / sum(Weight))
        result <- result[order(result$Specificity, decreasing = T), ] %>%
            slice_head(n = min(nrow(result), topk))
        ggplot(result, aes(x = Receptor, y = TF, fill = Specificity)) +
            geom_tile(color = "black") +
            scale_fill_gradientn(colors = c("yellow", "blue")) +
            coord_fixed() +
            theme(text = element_text(size = 24))
    } else {
        result <- Aggr_XT_results[Aggr_XT_results$TF == feature, ]
        result <- result[, c("Ligand", "Receptor", "Weight")] %>%
            mutate(Fidelity = Weight / sum(Weight))
        result <- result[order(result$Fidelity, decreasing = T), ] %>%
            slice_head(n = min(nrow(result), topk))
        ggplot(result, aes(x = Ligand, y = Receptor, fill = Fidelity)) +
            geom_tile(color = "black") +
            scale_fill_gradientn(colors = c("yellow", "blue")) +
            coord_fixed() +
            theme(text = element_text(size = 24))
    }

}
