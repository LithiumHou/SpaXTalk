#' Aggregate cell-cell communication results
#'
#' Combine per-receiver communication matrices and aggregate scores across
#' senders, receivers, or ligand-receptor pairs.
#'
#' @param CCC_origin A list of ligand-receptor communication matrices.
#' @param flag Character scalar specifying the aggregation dimension: `"Receiver"`,
#'   `"Sender"`, or `"LR"`.
#' @param LR_used Optional character vector of ligand-receptor pairs to retain.
#'
#' @return An aggregated numeric matrix whose dimensions depend on `flag`.
Aggr_CCC_results <- function(CCC_origin, flag = "Receiver", LR_used = NULL) {

    full_mat <- do.call(cbind, CCC_origin) %>% as.matrix()

    if (!is.null(LR_used)) {
        LR_used <- intersect(LR_used, rownames(full_mat))
        full_mat <- full_mat[LR_used, ] %>% as.matrix()
    }
    if (length(LR_used) == 1) {
        full_mat <- t(full_mat)
    }
    senders <- sub("=.*", "", colnames(full_mat))
    receivers <- sub(".*=", "", colnames(full_mat))

    # output: sender * LR, flatten the receiver layer
    if (flag == "Receiver") {
        agg_mat <- rowsum(t(full_mat), group = senders)
    }

    # output: receiver * LR, flatten the sender layer
    if (flag == "Sender") {
        agg_mat <- rowsum(t(full_mat), group = receivers)
    }

    # output: sender * receiver, flatten the LR layer
    if (flag == "LR") {
        df_long <- data.frame(
            Sender = senders,
            Receiver = receivers,
            value = colSums(full_mat)
        )
        row_levels <- unique(df_long$Sender)
        col_levels <- unique(df_long$Receiver)
        agg_mat <- Matrix::sparseMatrix(
            i = match(df_long$Sender, row_levels),
            j = match(df_long$Receiver, col_levels),
            x = df_long$value,
            dimnames = list(row_levels, col_levels)
        )
    }
    agg_mat[is.na(agg_mat)] <- 0
    return(agg_mat)
}

#' Process gene regulatory network results
#'
#' Convert a GRN coefficient matrix to long form and retain regulatory effects
#' whose absolute magnitude exceeds a relative threshold.
#'
#' @param GRN_original A list containing a coefficient matrix named `beta`.
#' @param beta_thres Numeric relative threshold applied to absolute beta values.
#'
#' @return A data frame containing `Gene`, `Cell`, `beta`, `Receptor`, and `TF`
#'   columns.
Process_GRN <- function(GRN_original, beta_thres = 0.005) {

    GRN_original$beta <- as.matrix(GRN_original$beta)
    GRN_I <- GRN_original$beta %>% tidyfst::mat_df()
    colnames(GRN_I) <- c("Gene", "Cell", "beta")

    GRN_filtered <- GRN_I[
        abs(GRN_I$beta) >= beta_thres * max(abs(GRN_I$beta)),
    ]
    GRN_filtered$Gene <- as.character(GRN_filtered$Gene)
    GRN_filtered$Receptor <- sub("=.*", "", GRN_filtered$Gene)
    GRN_filtered$TF <- sub(".*=", "", GRN_filtered$Gene)
    return(GRN_filtered)
}

#' Aggregate extracellular-to-transcriptional results
#'
#' Join aggregated ligand-receptor activity with filtered receptor-to-TF effects
#' to construct weighted ligand-receptor-TF relationships.
#'
#' @param agg_df A long-form communication data frame with cell, ligand-receptor,
#'   and score columns.
#' @param GRN_filtered A processed GRN data frame returned by [Process_GRN()].
#' @param threshold Numeric minimum communication score to retain.
#'
#' @return A data frame with cell, ligand, receptor, TF, communication score,
#'   regulatory coefficient, and combined weight columns.
Aggr_XT_results <- function(agg_df, GRN_filtered, threshold = 1e-4) {
    # agg_df: CCC df
    colnames(agg_df) <- c("Cell", "LR", "Value")
    agg_df <- filter(agg_df, Value > threshold)
    agg_df$Ligand <- sub("=.*", "", agg_df$LR)
    agg_df$Receptor <- sub(".*=", "", agg_df$LR)
    agg_df$comb <- paste0(agg_df$Receptor, "_", agg_df$Cell)
    GRN_filtered$comb <- paste0(GRN_filtered$Receptor, "_", GRN_filtered$Cell)
    merged_df <- merge(agg_df, GRN_filtered, by = "comb")
    merged_df <- merged_df[, c("Cell.x", "Ligand", "Receptor.x", "TF", "Value", "beta")]
    colnames(merged_df) <- c("Cell", "Ligand", "Receptor", "TF", "Value", "beta")
    merged_df <- merged_df[!duplicated(merged_df), ]
    merged_df$Weight <- abs(merged_df$Value * merged_df$beta)

    return(merged_df)
}

#' Aggregate extracellular-to-transcriptional results by sender type
#'
#' Summarize ligand-receptor communication by sender annotation and join it with
#' receptor-to-TF regulatory effects for each receiver cell.
#'
#' @param CCC_origin A list of ligand-receptor communication matrices, one per
#'   receiver cell.
#' @param GRN_filtered A processed GRN data frame returned by [Process_GRN()].
#' @param meta A cell metadata data frame containing an `annotation` column.
#' @param target Optional character vector of TFs to retain.
#'
#' @return A data frame containing ligand, receiver cell, sender type, TF, and
#'   combined weight.
Aggr_XT_senders <- function(CCC_origin, GRN_filtered, meta, target = NULL) {

    if (!is.null(target)) {
        GRN_target <- filter(GRN_filtered, TF %in% target)
    } else {
        GRN_target <- GRN_filtered
    }
    XT_senders <- list()
    for (j in 1:length(CCC_origin)) {
        x <- CCC_origin[[j]]
        temp <- colnames(x)
        if (!is.null(temp) && sum(x) > 0) {
            temp <- strsplit(temp, split = "=")
            temp <- purrr::map(temp, 1) %>% unlist()
            sample_type <- meta[temp, "annotation"]
            type_levels <- unique(sample_type)
            type_mat <- Matrix::Matrix(
                0, nrow = length(sample_type), ncol = length(type_levels),
                dimnames = list(NULL, type_levels), sparse = TRUE
            )

            for (i in seq_along(sample_type)) {
                type_mat[i, sample_type[i]] <- 1
            }
            agg_mat <- x %*% type_mat
            agg_df <- tidyfst::mat_df(as.matrix(agg_mat))
            agg_df <- filter(agg_df, value > 0.01 * max(agg_df$value))
            colnames(agg_df) <- c("LR", "Type", "value")
            agg_df$receiver <- names(CCC_origin)[j]
        } else {
            agg_df <- data.frame(
                LR = character(), Type = character(), value = numeric(),
                receiver = character()
            )
        }
        XT_senders[[j]] <- agg_df
    }
    XT_senders <- do.call(rbind, XT_senders)

    XT_sender_aggr <- aggregate(
        XT_senders$value,
        list(XT_senders$LR, XT_senders$Type, XT_senders$receiver),
        sum
    )
    colnames(XT_sender_aggr) <- c("LR", "Type", "Cell", "Value1")
    XT_sender_aggr$LR <- as.character(XT_sender_aggr$LR)
    XT_sender_aggr$Ligand <- sapply(XT_sender_aggr$LR, function(x) {
        temp <- strsplit(x, split = "=", fixed = T)
        unlist(temp)[[1]]
    })
    XT_sender_aggr$Receptor <- sapply(XT_sender_aggr$LR, function(x) {
        temp <- strsplit(x, split = "=", fixed = T)
        unlist(temp)[[2]]
    })

    GRN_target_use <- GRN_target[, c("Receptor", "Cell", "TF", "beta")]
    CC_results2 <- inner_join(
        XT_sender_aggr, GRN_target_use, by = c("Receptor", "Cell")
    )
    CC_results2 <- distinct(CC_results2)

    results <- data.frame(
        Ligand = CC_results2$Ligand,
        Cell = CC_results2$Cell,
        Type = CC_results2$Type,
        TF = CC_results2$TF,
        Weight = CC_results2$Value1 * CC_results2$beta
    )

    return(results)
}

#' Calculate pathway activity by sender
#'
#' Combine communication and regulatory effects for specified
#' ligand-receptor-TF paths and retain the originating sender cells.
#'
#' @param CCC_origin A list of ligand-receptor communication matrices.
#' @param GRN_filtered A processed GRN data frame returned by [Process_GRN()].
#' @param pathway_names Character vector of ligand-receptor-TF path identifiers.
#'
#' @return A data frame containing sender, receiver cell, ligand, receptor, TF,
#'   and path weight.
Pathway_Sender <- function(CCC_origin, GRN_filtered, pathway_names) {

    pathway_names <- intersect(pathway_names, GRN_filtered$Gene)

    results_all <- lapply(pathway_names, function(x) {

        LRT <- strsplit(x, split = "=") %>% unlist()
        LR <- paste0(LRT[1], "=", LRT[2])
        RT <- paste0(LRT[2], "=", LRT[3])
        results_list <- lapply(CCC_origin, function(x) {
            df <- data.frame(SR = colnames(x), Strength = x[LR, ] %>% as.vector())
        })
        results_df1 <- do.call(rbind, results_df)
        results_df1$Sender <- sub("=.*", "", results_df$SR)
        results_df1$Cell <- sub(".*=", "", results_df$SR)
        results_df1$Ligand <- LRT[1]
        results_df1$Receptor <- LRT[2]
        results_df2 <- filter(GRN_filtered, Gene == x)

        results_merged <- merge(results_df1, results_df2, by = c("Receptor", "Cell"))
        results_merged$Weight <- results_merged$beta * results_merged$SR
        results_merged <- results_merged[
            , c("Sender", "Cell", "Ligand", "Receptor", "TF", "Weight")
        ]
        results_merged
    })

    results_all <- do.call(rbind, results_all)
    return(results_all)
}

#' Group spatial signals and targets
#'
#' Partition receiver cells into spatial groups and summarize ligand and TF
#' weights within each group.
#'
#' @param results A ligand-receptor-TF result data frame containing `Cell`,
#'   `Ligand`, `TF`, and `Weight` columns.
#' @param meta A cell metadata data frame containing `x` and `y` coordinates.
#' @param features Optional character vector of TFs to include.
#' @param groups Optional named list of cell vectors defining precomputed groups.
#' @param cluster_k Integer number of spatial groups to create when `groups` is
#'   not supplied.
#'
#' @return A named list containing the cells and aggregated signal and target
#'   tables for each group.
Group_signals <- function(results, meta, features = NULL, groups = NULL, cluster_k = 3) {

    # Format of results: part of the XT_results
    # Format of group: a list of cells grouped by the cluster labels
    if (is.null(features)) {
        features <- results$TF %>% unique()
    }
    results_filtered <- filter(results, TF %in% features)
    cells <- unique(results_filtered$Cell) %>% as.character()

    # Perform the k means clustering
    if (is.null(groups) || (!is.list(groups))) {
        meta_used <- meta[cells, c("x", "y")]
        km_res <- kmeans(meta_used, centers = cluster_k)
        groups <- split(rownames(meta_used), km_res$cluster)
    }
    group_names <- names(groups) <- paste0("Group", 1:cluster_k)
    mylist <- lapply(group_names, function(x) {
        cells <- groups[[x]]
        temp <- results_filtered[results_filtered$Cell %in% cells, ]
        signal_s <- aggregate(temp$Weight, list(temp$Ligand), sum)
        colnames(signal_s) <- c("Signal", "Value")
        target_s <- aggregate(temp$Weight, list(temp$TF), sum)
        colnames(target_s) <- c("Target", "Value")

        list(Cells = cells, Signals = signal_s, Targets = target_s)
    })
    names(mylist) <- group_names
    return(mylist)
}

#' Convert pathways to receptor-TF pairs
#'
#' Retrieve the unique receptor-TF pairs associated with one or more pathways.
#'
#' @param Pathway_name Character vector of pathway names.
#' @param DB A pathway database containing `pathway`, `receptor`, and `tf`
#'   columns.
#'
#' @return A character vector of receptor-TF pairs separated by `"="`.
Path_To_Pair_GRN <- function(Pathway_name, DB) {

    DB_filtered <- filter(DB, pathway %in% Pathway_name)
    RTs <- paste(DB_filtered$receptor, DB_filtered$tf, sep = "=") %>% unique()

    return(RTs)
}

#' Convert pathways to ligand-receptor pairs
#'
#' Resolve pathway names or explicit ligand-receptor identifiers to the unique
#' ligand-receptor pairs represented in a communication database.
#'
#' @param Pathway_name Character vector of pathway names or ligand-receptor pairs.
#' @param DB A communication database containing ligand, receptor, and pathway
#'   information.
#'
#' @return A sorted character vector of unique ligand-receptor pairs.
Path_To_Pair_CCC <- function(Pathway_name, DB) {

    if (!"LR" %in% colnames(DB)) {
        DB$LR <- paste(DB$From, DB_filtered$To, sep = "=")
    }
    Pathway_name <- unique(Pathway_name)
    # Find the LR pairs that matches the pathway name
    LRP1 <- Pathway_name[Pathway_name %in% DB$LR]
    pathways <- Pathway_name[Pathway_name %in% DB$Pathway] %>% unique()

    if (length(pathways)) {
        LRP2 <- filter(DB, Pathway %in% pathways) %>% pull(LR)
    }
    LRs <- c(LRP1, LRP2) %>% unique() %>% sort()
    return(LRs)
}

#' Calculate pathway crosstalk coverage
#'
#' Compare two sets of ligand-receptor pairs across cells and TFs to quantify
#' their overlap and combined regulatory strength.
#'
#' @param XT_results A ligand-receptor-TF result data frame containing activity
#'   and regulatory weights.
#' @param pathway1 Character vector of ligand-receptor pairs in the first pathway.
#' @param pathway2 Character vector of ligand-receptor pairs in the second pathway.
#'
#' @return A data frame with TF target, cell coverage, and combined GRN strength,
#'   sorted by decreasing coverage.
Calculate_Coverage <- function(XT_results, pathway1, pathway2) {

    XT_results$LR <- paste(XT_results$Ligand, XT_results$Receptor, sep = "=")
    XT1 <- filter(XT_results, LR %in% pathway1)
    XT1 <- aggregate(cbind(Value, beta, Weight) ~ Cell + TF, data = XT1, FUN = sum)
    XT2 <- filter(XT_results, LR %in% pathway2)
    XT2 <- aggregate(cbind(Value, beta, Weight) ~ Cell + TF, data = XT2, FUN = sum)

    LR1_active_cells <- XT1$Cell %>% unique()
    LR2_active_cells <- XT2$Cell %>% unique()
    crosstalk_cells <- intersect(LR1_active_cells, LR2_active_cells)
    if (!length(crosstalk_cells)) {
        stop("There is no crosstalk between two pathways!")
    }

    TF_used <- intersect(XT1$TF, XT2$TF) %>% sort()
    cov_list <- lapply(TF_used, function(x) {
        temp1 <- filter(XT1, TF == x)
        temp2 <- filter(XT2, TF == x)
        colnames(temp2) <- c("Cell", "TF", "Value2", "beta2", "Weight2")
        temp1$Path <- "Pathway1"
        temp2$Path <- "Pathway2"
        # compare the overlap
        crosstalk_cells_used <- intersect(temp1$Cell, temp2$Cell) %>%
            intersect(crosstalk_cells)

        # ---------Need to fix this in the future----------
        strength_df <- temp1 %>%
            full_join(temp2, by = "Cell") %>%
            mutate(
                beta = coalesce(beta, 0),
                beta2 = coalesce(beta2, 0),
                combined_GRN = beta * beta2
            ) %>%
            arrange(Cell)
        data.frame(
            Target = x,
            Coverage = length(crosstalk_cells_used) / length(crosstalk_cells),
            Strength = sum(strength_df$combined_GRN)
        )
    })
    cov_list <- do.call(rbind, cov_list)
    cov_list <- cov_list %>% arrange(desc(Coverage))
    return(cov_list)
}

#' Calculate distance to an annotation border
#'
#' For each cell, calculate the distance to the nearest cell with a different
#' annotation.
#'
#' @param meta A cell metadata data frame containing `x`, `y`, and `annotation`
#'   columns, with cell identifiers as row names.
#'
#' @return A named numeric vector giving the nearest-border distance for each
#'   cell.
Dist_to_border <- function(meta) {
    dist_mat <- dist(meta[, c("x", "y")]) %>% as.matrix()
    result <- lapply(rownames(dist_mat), function(x) {
        hetero_cells <- meta[, "annotation"] != meta[x, "annotation"]
        min(dist_mat[x, hetero_cells])
    }) %>% unlist()
    names(result) <- rownames(dist_mat)
    return(result)
}

#' Calculate sender-type contributions
#'
#' Aggregate a cell-level communication matrix by sender annotation and
#' normalize the contribution of each annotation within each ligand-receptor row.
#'
#' @param CCC_LR A communication matrix with sender cells as rows.
#' @param meta A cell metadata data frame containing an `annotation` column.
#'
#' @return A row-normalized matrix of sender-annotation contributions.
Sender_contribution <- function(CCC_LR, meta) {
    annos <- meta[rownames(CCC_LR), "annotation"] %>% as.vector()
    result <- rowsum(as.matrix(CCC_LR), group = annos) %>% t()
    result <- result / rowSums(result)
    return(result)
}

#' Group cells by TF specificity
#'
#' Calculate cell-level TF specificity and cluster cells according to the
#' specificity of a selected TF.
#'
#' @param XT_results A ligand-receptor-TF result data frame containing `Cell`,
#'   `TF`, and `Weight` columns.
#' @param feature Character scalar naming the TF to group.
#' @param cluster_k Integer number of specificity groups to create.
#'
#' @return A named list containing the cells and minimum and maximum specificity
#'   values for each group.
Group_targets <- function(XT_results, feature, cluster_k = 4) {
    result <- aggregate(Weight ~ Cell + TF, data = XT_results, sum)
    result <- result %>%
        group_by(Cell) %>%
        mutate(Spe = Weight / sum(Weight)) %>%
        ungroup()

    if (is.null(feature) | length(feature) > 1) {
        stop("Should input one feature!")
    }

    results_filtered <- filter(result, TF == feature)
    df <- data.frame(Value = results_filtered$Spe, row.names = results_filtered$Cell)
    km_res <- kmeans(df, centers = cluster_k)
    groups <- split(rownames(df), km_res$cluster)
    group_names <- names(groups) <- paste0("Group", 1:cluster_k)
    mylist <- lapply(group_names, function(x) {
        cells <- groups[[x]]
        temp <- results_filtered[results_filtered$Cell %in% cells, ]
        list(
            Cells = cells,
            min_spe = min(temp$Spe),
            max_spe = max(temp$Spe)
        )
    })

    names(mylist) <- group_names
    return(mylist)
}

#' Fit group-specific TF regression models
#'
#' Select receptors associated with a TF and fit an elastic-net expression model
#' separately for each cell group.
#'
#' @param XT_results A ligand-receptor-TF result data frame containing receptor,
#'   TF, and specificity information.
#' @param group_results A named list of groups, each containing a `Cells` vector.
#' @param exp A gene-by-cell expression matrix.
#' @param feature Character scalar naming the TF response gene.
#' @param top_k Integer maximum number of receptors used as predictors.
#'
#' @return A named list of fitted `glmnet` models; entries are `NULL` when model
#'   fitting fails.
Regression_target <- function(XT_results, group_results, exp, feature, top_k = 10) {
    group_names <- names(group_results)
    XT_filtered <- XT_results[XT_results$TF == feature, ]
    temp <- aggregate(Spe ~ Receptor, data = XT_filtered, mean)
    top_k <- min(nrow(temp), top_k)
    receptors <- temp %>% arrange(Spe) %>% slice_tail(n = top_k) %>% pull(Receptor)
    y_all <- exp[feature, ]
    x_all <- exp[receptors, ]
    result <- lapply(group_results, function(x) {
        mymodel <- tryCatch({
            cell_used <- x$Cells
            y <- y_all[cell_used]
            x <- x_all[, cell_used]
            cvfit <- glmnet::cv.glmnet(
                t(x), y, alpha = 0.5, nfolds = 3, type.measure = "mse"
            )
            best_lambda <- cvfit$lambda.min
            model <- glmnet::glmnet(
                t(x), y, alpha = 0.5, lambda = best_lambda
            )
            model
        }, error = function(e) NULL)
        mymodel
    })
    names(result) <- group_names
    return(result)
}
