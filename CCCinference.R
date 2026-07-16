load_fit_spatial_local_effect <- function(dll_path = NULL) {
    dll_loaded <- any(vapply(getLoadedDLLs(), function(x) {
        identical(x[["name"]], "fit_spatial_local_effect")
    }, logical(1)))
    if (dll_loaded) return(invisible(TRUE))

    if (is.null(dll_path)) {
        candidates <- c(
            file.path(getwd(), "codes", "fit_spatial_local_effect.so"),
            file.path(getwd(), "fit_spatial_local_effect.so"),
            file.path(dirname(normalizePath("codes/CCCinference.R", mustWork = FALSE)), "fit_spatial_local_effect.so")
        )
        candidates <- unique(normalizePath(candidates, mustWork = FALSE))
        dll_path <- candidates[file.exists(candidates)][1]
    }

    if (is.na(dll_path) || !file.exists(dll_path)) {
        stop("Cannot find fit_spatial_local_effect.so. Compile it with: R CMD SHLIB codes/fit_spatial_local_effect.c")
    }
    dyn.load(dll_path)
    invisible(TRUE)
}

solve_spatial_ridge <- function(lhs, rhs, fallback_len) {
    solve_once <- function(mat){
        tryCatch(
            Matrix::solve(mat, rhs),
            error = function(e) NULL
        )
    }

    sol <- solve_once(lhs)
    if(!is.null(sol) && all(is.finite(sol))) return(sol)

    diag_scale <- suppressWarnings(max(abs(as.numeric(Matrix::diag(lhs))), na.rm = TRUE))
    if(!is.finite(diag_scale) || diag_scale <= 0) diag_scale <- 1
    ridge_values <- diag_scale * c(1e-8, 1e-6, 1e-4, 1e-2)
    I <- Matrix::Diagonal(nrow(lhs))
    for(ridge in ridge_values){
        sol <- solve_once(lhs + ridge * I)
        if(!is.null(sol) && all(is.finite(sol))) return(sol)
    }

    rep(NA_real_, fallback_len)
}

fit_spatial_local_effect_batch <- function(data_mat, pair_i, pair_j, neighbors, lambda = 1,
                                           progress = TRUE, update_every = 50) {
    load_fit_spatial_local_effect()
    .Call(
        "fit_spatial_local_effect_batch",
        data_mat, as.integer(pair_i), as.integer(pair_j),
        neighbors, as.numeric(lambda),
        as.logical(progress), as.integer(update_every)
    )
}

fit_multi_spatial_effect_batch <- function(data_mat, RTDB_used, neighbors, nsim = 99, lambda = 1,
                                           progress = TRUE, update_every = 50) {
    load_fit_spatial_local_effect()
    RTDB_used <- RTDB_used[, c("From", "To")]
    target_groups <- split(seq_len(nrow(RTDB_used)), RTDB_used$To)
    target_names <- names(target_groups)
    target_i <- match(target_names, rownames(data_mat))
    predictors_list <- lapply(target_groups, function(idx) {
        as.integer(match(RTDB_used$From[idx], rownames(data_mat)))
    })
    edge_group <- integer(nrow(RTDB_used))
    edge_pred_pos <- integer(nrow(RTDB_used))
    for(g in seq_along(target_groups)){
        idx <- target_groups[[g]]
        edge_group[idx] <- g
        edge_pred_pos[idx] <- seq_along(idx)
    }
    .Call(
        "fit_multi_spatial_effect_batch",
        data_mat, as.integer(target_i), predictors_list,
        as.integer(edge_group), as.integer(edge_pred_pos),
        neighbors, as.numeric(lambda), as.integer(nsim),
        as.logical(progress), as.integer(update_every)
    )
}

Prepare_data <- function(Object, receiver_cluster = NULL, spacex_flag = "x", spacey_flag = "y", cluster_flag = "annotation", assay_type = "SCT", layer_type = "data", filter_genes = FALSE){

    if(assay_type %in% names(Object@assays)){
        if(layer_type %in% slotNames(Object@assays[[assay_type]])){
            data_raw <- GetAssayData(Object, assay = assay_type, layer = layer_type)
        }else data_raw <- GetAssayData(Object, assay = assay_type, layer = "counts")
    }else{
        if(layer_type %in% slotNames(Object@assays[[assay_type]])){
            data_raw <- GetAssayData(Object, assay = "RNA", layer = layer_type)
        }else data_raw <- GetAssayData(Object, assay = "RNA", layer = "counts")
    }
    genes_raw <- rownames(data_raw)

    meta <- Object@meta.data
    meta <- meta[, c(spacex_flag,spacey_flag,cluster_flag)]
    colnames(meta) <- c("x","y","annotation")

    if(filter_genes){
        patterns <- c("^Hba-","^Hbb-","^Rps","^Rpl","^mt-")
        patterns <- paste(patterns, collapse = "|")
        genes_filtered <- genes_raw[!grepl(patterns, genes_raw)]
    }else genes_filtered <- genes_raw

    exp_data <<- data_raw[genes_filtered, ]
    meta <<- meta
    if(!is.null(receiver_cluster)){
        receiver <<- filter(meta, annotation %in% receiver_cluster) %>% rownames()
    }else receiver <<- rownames(meta)

}

Identify_neighbor_list <- function(locs, senders = NULL, receivers = NULL, dist_threshold = 5, self_neighbor = F,
    decay_type = NULL, decay_threshold = 1e-3, normalize = T,nclusters = 1){

    # decay type: exponential, reciprocal, inverse square, constant
    if(is.null(decay_type)){
        decay_type <- "constant"
    } else decay_type <- stringr::str_to_lower(decay_type)
    if(decay_type == "exp"){
        message("Using exponential decay of CCC strength coefficient ... \n")
        decay_type <- "exponential"
    } else if(decay_type == "rec"){
        message("Using reciprocal decay of CCC strength coefficient ... \n")
        decay_type <- "reciprocal"
    } else if(decay_type == "insq"){
        message("Using inverse-square decay of CCC strength coefficient ... \n")
        decay_type <- "inverse_square"
    } else if(decay_type == "rbf"){
        message("Using RBF decay of CCC strength coefficient ... \n")
        decay_type <- "RBF"
    } else {
        message("Using constant CCC strength coefficient ... \n")
        decay_type <- "constant"
    }

    allcells <- rownames(locs)
    if(is.null(senders)) senders <- allcells
    if(is.null(receivers)) receivers <- allcells

    senders <- intersect(senders, allcells)
    receivers <- intersect(receivers, allcells)

    dist_mat <- dist(locs) %>% as.matrix()
    if(!self_neighbor) diag(dist_mat) <- Inf # Will not include self
    
    message("Infer the receptive field of each cell ... \n")
    if(is.null(nclusters))
        nclusters <- ceiling(future::availableCores() / 2)
    if(nclusters > 1){
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(future::multisession, workers = nclusters)
    }

    infer_neighbors_one <- function(x){
        p()
        temp <- as.array(dist_mat[x, senders])
        names(temp) <- senders
        sender_dist <- sort(temp[temp <= dist_threshold])
        if(length(sender_dist) == 0) return(numeric())

        if(self_neighbor){
            decay_vector <- exp(-(sender_dist)^2)
        } else{
            if(decay_type == "exponential"){
                decay_vector <- exp(-(sender_dist)^2)
            }
            if(decay_type == "reciprocal"){
                decay_vector <- (sender_dist)^(-1)
            }
            if(decay_type == "inverse_square"){
                decay_vector <- (sender_dist)^(-2)
            }
            if(decay_type == "RBF"){
                id <- length(sender_dist)/3
                bandwidth <- sender_dist[ceiling(id)]
                decay_vector <- exp(bandwidth^2*(sender_dist)^(-2))
            }
            if(decay_type == "constant"){
                decay_vector <- rep(1, length(sender_dist))
            }
        }

        names(decay_vector) <- names(sender_dist)
        if(normalize){
            temp <- decay_vector/sum(decay_vector)
        }else{
            temp <- decay_vector
            temp <- temp[temp >= decay_threshold*sum(decay_vector)]
        }
        temp
    }

    with_progress({
        p <- progressor(steps = length(receivers))
        if(nclusters > 1){ # We recommend using parallel computing
            message("Using parallel computing ... ")
            neighbor_list <- future.apply::future_lapply(receivers, infer_neighbors_one, future.seed = F)
        }else{
            message("Using non-parallel ... ")
            neighbor_list <- lapply(receivers, infer_neighbors_one)
        }
    })
    names(neighbor_list) <- receivers
    return(neighbor_list)
    
}

InferCCC_Scriabin <- function(data, neighbor_list, LRDB, ligands = NULL, receptors = NULL, nclusters = 1){

    # Filter the dataset and database
    receivers <- names(neighbor_list)
    allgenes <- rownames(data)
    if(is.null(ligands)) ligands <- LRDB$From %>% unique()
    if(is.null(receptors)) receptors <- LRDB$To %>% unique()
    ligands <- intersect(ligands, allgenes); receptors <- intersect(receptors, allgenes)

    LRDB <- LRDB[!duplicated(LRDB), ]
    LRDB_filtered <- filter(LRDB, From %in% ligands) 
    LRDB_filtered <- filter(LRDB_filtered, To %in% receptors)

    if(is.null(nclusters))
        nclusters <- ceiling(future::availableCores() / 2)
    if(nclusters > 1){
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        future::plan(future::multisession, workers = nclusters)
    }

    message("Inferring the cell-cell communication patterns ... \n")
    infer_ccc_one <- function(x){
        p()
        decay_vector <- neighbor_list[[x]]
        if(length(decay_vector)){
            ligand_exp <- data[LRDB_filtered$From, names(decay_vector), drop = FALSE]
            ligand_exp[is.na(ligand_exp)] <- 0
            receptor_exp <- data[LRDB_filtered$To, x]
            receptor_exp[is.na(receptor_exp)] <- 0

            LR_mat <- t(t(ligand_exp) * as.vector(decay_vector)) * receptor_exp
            rownames(LR_mat) <- paste0(LRDB_filtered$From,"=",LRDB_filtered$To)
            colnames(LR_mat) <- paste0(names(decay_vector),"=",x)
            
        } else {
            LR_mat <- matrix(nrow = nrow(LRDB_filtered), ncol = 0)
            rownames(LR_mat) <- paste0(LRDB_filtered$From,"=",LRDB_filtered$To)
        }
        LR_mat
    }

    with_progress({
        p <- progressor(steps = length(receivers))
        if(nclusters > 1){
            message("Using parallel computing ... ")
            LR_list <- future.apply::future_lapply(receivers, infer_ccc_one, future.seed = F, future.packages = c("Matrix"))
        } else {
            message("Using non-parallel ... ")
            LR_list <- lapply(receivers, infer_ccc_one)
        }
    })

    names(LR_list) <- receivers
    return(LR_list)
}


InferGRN <- function(exps, neighbor_list, RTDB = NULL, Target_list = NULL, exp_thres = 0.01, PPR_top = 10, rec_only = F, nsim = 99, useC = T, lambda = 1, progress = TRUE, update_every = 3, nclusters = 1, use_delta_y = T){

    # Filter the neighbors
    allcells <- colnames(exps)
    receivers <- names(neighbor_list)

    message("Formulating the neighborhood list ... \n")
    nl_filtered <- lapply(neighbor_list, function(x){
        if(rec_only){
            xf <- x[names(x) %in% receivers]
        }else xf <- x     
        xf
    })

    # Filter the data with specific genes (expressed in >1% of the receivers)
    colnames(RTDB)[1:2] <- c("From","To")
    allgenes <- intersect(c(RTDB$From, RTDB$To), rownames(exps))
    if(!is.null(Target_list)) {
        allgenes <- intersect(allgenes, Target_list)
    }
    data_used <- exps[allgenes, allcells]
    data_used <- data_used[rowSums(data_used[ ,receivers] > 0) > exp_thres*length(receivers), ]
    # c(nrow(exps),nrow(data_used))
    RTDB_used <- filter(RTDB, From %in% rownames(data_used)) 
    RTDB_used <- filter(RTDB_used, To %in% rownames(data_used))
    RTDB_used <- RTDB_used %>% 
        group_by(To) %>%
        arrange(PPR) %>%
        slice_head(n = PPR_top) %>%
        ungroup()
    RTDB_used <- RTDB_used[, c("From", "To")]

    message("Inferring the pointwise GRN using multi-variable spatial regression ... \n")
    if(useC && use_delta_y){
        message("use_delta_y = TRUE requires the sparse Matrix solver; setting useC = FALSE. ")
        useC <- FALSE
    }
    if(useC){
        exps <- as.matrix(exps)
        storage.mode(exps) <- "double"
        results <- fit_multi_spatial_effect_batch(exps, RTDB_used, nl_filtered, nsim, lambda, progress, update_every)
    } else {
        if(!inherits(data_used, "sparseMatrix")){
            data_used <- as.matrix(data_used)
            storage.mode(data_used) <- "double"
        }
        results <- fit_multi_spatial_effect(data_used, RTDB_used, nl_filtered, nsim, lambda, progress, update_every, nclusters, use_delta_y)
    }

    alpha_df <- results$alpha %>% t()
    alpha_df[is.na(alpha_df)] <- 0
    beta_df <- results$beta %>% t()
    beta_df[is.na(beta_df)] <- 0
    pvalue_df <- results$pvalue %>% t()
    pvalue_df[is.na(pvalue_df)] <- 1
    rownames(alpha_df) <- rownames(beta_df) <- rownames(pvalue_df) <- paste0(RTDB_used$From, "=", RTDB_used$To)

    return(list(alpha = alpha_df, beta = beta_df, pvalue = pvalue_df))

}



Run_SCENIC <- function(SeuratObj, work_dir, tf_list = NULL, species = "human") {

    name_used <- names(SeuratObj@assays[["RNA"]]@layers)
    filen2 <- paste0(work_dir,"input_data.h5ad")
    scCustomize::as.anndata(SeuratObj, assay = "RNA", file_path = ".",
                            main_layer = name_used[1], other_layers = NULL,
                            file_name = filen2)
    if(!is.null(tf_list)){
        write.table(tf_list, file = paste0(work_dir,"tf_list.txt"), quote = F, row.names = F, col.names = F)
    }else if (species == "human"){
        data("human_TF_list")
        write.table(human_TF_list, file = paste0(work_dir,"tf_list.txt"), quote = F, row.names = F, col.names = F)
    }else if (species == "mouse"){
        data("mouse_TF_list")
        write.table(mouse_TF_list, file = paste0(work_dir,"tf_list.txt"), quote = F, row.names = F, col.names = F)
    }else{
        stop("Invalid species. Please choose one of: human, mouse, or provide your own TF list.\n")
    }

    # run scenic
    file_py <- paste0(work_dir, "run_scenic.py")
    reticulate::use_condaenv("SigXTalk_py")
    reticulate::source_python(file_py)
    
}

fit_multi_spatial_effect <- function(data_mat, RTDB_used, neighbors, nsim = 99, lambda = 1, progress = TRUE, update_every = 50, nclusters = 1, use_delta_y = FALSE) {
    # Multi-variable spatial regression for InferGRN.
    # For each target in RTDB_used$To, all matching RTDB_used$From genes are used as X.
    if(!inherits(data_mat, "sparseMatrix")){
        data_mat <- as.matrix(data_mat)
        storage.mode(data_mat) <- "double"
    }
    RTDB_used <- RTDB_used[, c("From", "To")]
    allcells <- colnames(data_mat)
    receivers <- names(neighbors)
    n <- length(allcells)
    edge_names <- paste0(RTDB_used$From, "=", RTDB_used$To)

    if(is.null(nclusters)){
        nclusters <- ceiling(future::availableCores() / 2)
    }
    nclusters <- as.integer(nclusters)
    if(is.na(nclusters) || nclusters < 1) nclusters <- 1

    w_i <- integer(0)
    w_j <- integer(0)
    w_x <- numeric(0)
    for(receiver in receivers){
        wt <- neighbors[[receiver]]
        receiver_idx <- match(receiver, allcells)
        if(length(wt) == 0 || is.na(receiver_idx)) next
        sender_idx <- match(names(wt), allcells)
        keep <- !is.na(sender_idx)
        if(!any(keep)) next
        w_i <- c(w_i, rep(receiver_idx, sum(keep)))
        w_j <- c(w_j, sender_idx[keep])
        w_x <- c(w_x, as.numeric(wt[keep]))
    }
    W <- Matrix::sparseMatrix(i = w_i, j = w_j, x = w_x,
                              dims = c(n, n), dimnames = list(allcells, allcells))
    L <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W

    alpha <- beta <- pvalue <- matrix(NA_real_, nrow = length(receivers), ncol = nrow(RTDB_used),
                                      dimnames = list(receivers, edge_names))
    target_groups <- split(seq_len(nrow(RTDB_used)), RTDB_used$To)
    target_names <- names(target_groups)

    fit_one_target <- function(target){
        edge_idx <- target_groups[[target]]
        predictors <- RTDB_used$From[edge_idx]
        y <- as.numeric(data_mat[target, allcells])
        if(use_delta_y){
            y <- y - mean(y, na.rm = TRUE)
        }
        p <- length(predictors)

        A <- Matrix::Matrix(1, nrow = n, ncol = 1, sparse = TRUE)
        for(k in seq_len(p)){
            xvals <- as.numeric(data_mat[predictors[k], allcells])
            A <- cbind(A, Matrix::Diagonal(x = xvals))
        }
        P <- do.call(Matrix::bdiag, c(list(Matrix::Matrix(0, 1, 1, sparse = TRUE)), replicate(p, lambda * L, simplify = FALSE)))
        lhs <- Matrix::crossprod(A) + P

        fit_once <- function(yvals){
            theta <- tryCatch(
                solve_spatial_ridge(lhs, Matrix::crossprod(A, as.numeric(yvals)), 1 + n * p),
                error = function(e) rep(NA_real_, 1 + n * p)
            )
            beta_mat <- matrix(as.numeric(theta[-1]), nrow = n, ncol = p)
            rownames(beta_mat) <- allcells
            colnames(beta_mat) <- predictors
            list(alpha = as.numeric(theta[1]), beta = beta_mat)
        }

        alpha_part <- beta_part <- pvalue_part <- matrix(NA_real_, nrow = length(receivers), ncol = length(edge_idx),
                                                         dimnames = list(receivers, edge_names[edge_idx]))

        obs <- fit_once(y)
        for(pos in seq_along(edge_idx)){
            alpha_part[, pos] <- obs$alpha
            beta_part[, pos] <- obs$beta[receivers, pos]
        }

        if(nsim > 0){
            perm_counts <- matrix(0, nrow = length(receivers), ncol = length(edge_idx),
                                  dimnames = list(receivers, edge_names[edge_idx]))
            valid_counts <- matrix(0, nrow = length(receivers), ncol = length(edge_idx),
                                   dimnames = list(receivers, edge_names[edge_idx]))
            for(i in seq_len(nsim)){
                perm <- fit_once(sample(y))$beta
                for(pos in seq_along(edge_idx)){
                    obs_beta <- beta_part[, pos]
                    perm_beta <- perm[receivers, pos]
                    valid <- !is.na(obs_beta) & !is.na(perm_beta)
                    perm_counts[valid, pos] <- perm_counts[valid, pos] + (abs(perm_beta[valid]) >= abs(obs_beta[valid]))
                    valid_counts[valid, pos] <- valid_counts[valid, pos] + 1
                }
            }
            for(pos in seq_along(edge_idx)){
                ok <- valid_counts[, pos] > 0
                pvalue_part[ok, pos] <- perm_counts[ok, pos] / valid_counts[ok, pos]
            }
        }

        list(edge_idx = edge_idx, alpha = alpha_part, beta = beta_part, pvalue = pvalue_part)
    }

    total_targets <- length(target_names)
    show_target_progress <- function(done){
        if(!isTRUE(progress) || total_targets == 0) return(invisible(NULL))
        width <- 40L
        filled <- floor(width * done / total_targets)
        bar <- paste0("[", strrep("=", filled), strrep(" ", width - filled), "]")
        cat(sprintf("\r%s %d/%d targets (y)", bar, done, total_targets))
        flush.console()
        if(done >= total_targets) cat("\n")
        invisible(NULL)
    }

    if(nclusters > 1 && total_targets > 1){
        if(!requireNamespace("future", quietly = TRUE)){
            stop("Parallel fit_multi_spatial_effect requires the future package.")
        }
        cat("Using parallel computing for multi-variable spatial regression ... \n")
        flush.console()
        show_target_progress(0L)
        workers <- min(nclusters, total_targets)
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        if(future::supportsMulticore()){
            future::plan(future::multicore, workers = workers)
        } else {
            future::plan(future::multisession, workers = workers)
        }
        target_chunks <- split(target_names, cut(seq_along(target_names), workers, labels = FALSE))
        futures <- lapply(target_chunks, function(target_chunk){
            future::future(lapply(target_chunk, fit_one_target), seed = TRUE, packages = c("Matrix"))
        })
        chunk_results <- vector("list", length(futures))
        finished <- rep(FALSE, length(futures))
        targets_done <- 0L
        while(!all(finished)){
            ready <- vapply(futures, future::resolved, logical(1))
            newly_finished <- which(ready & !finished)
            if(length(newly_finished) > 0){
                for(idx in newly_finished){
                    chunk_results[[idx]] <- future::value(futures[[idx]])
                    targets_done <- targets_done + length(chunk_results[[idx]])
                }
                finished[newly_finished] <- TRUE
                show_target_progress(targets_done)
            }
            if(!all(finished)) Sys.sleep(0.2)
        }
        result_list <- unlist(chunk_results, recursive = FALSE)
    } else {
        result_list <- vector("list", total_targets)
        names(result_list) <- target_names
        show_target_progress(0L)
        for(i in seq_along(target_names)){
            result_list[[i]] <- fit_one_target(target_names[[i]])
            show_target_progress(i)
        }
    }

    for(result in result_list){
        alpha[, result$edge_idx] <- result$alpha
        beta[, result$edge_idx] <- result$beta
        pvalue[, result$edge_idx] <- result$pvalue
    }

    list(alpha = alpha, beta = beta, pvalue = pvalue, lambda = lambda)
}
