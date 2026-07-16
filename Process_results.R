Aggr_CCC_results <- function(CCC_origin, flag = "Receiver", LR_used = NULL){

    full_mat <- do.call(cbind,CCC_origin) %>% as.matrix()

    if(!is.null(LR_used)){
        LR_used <- intersect(LR_used, rownames(full_mat))
        full_mat <- full_mat[LR_used, ] %>% as.matrix()
    }
    if(length(LR_used) == 1) full_mat <- t(full_mat)
    senders <- sub("=.*", "", colnames(full_mat))
    receivers <- sub(".*=", "", colnames(full_mat))

    # output: sender*LR, flatten the receiver layer
    if(flag == "Receiver"){
        agg_mat <- rowsum(t(full_mat), group = senders)
    }

    # output: receiver*LR, flatten the sender layer
    if(flag == "Sender"){
        agg_mat <- rowsum(t(full_mat), group = receivers)
    }

    # output: sender*receiver, flatten the LR layer
    if(flag == "LR"){

        df_long <- data.frame(
            Sender = senders,
            Receiver = receivers,
            value = colSums(full_mat)
        )
        row_levels <- unique(df_long$Sender)
        col_levels <- unique(df_long$Receiver)
        agg_mat <- Matrix::sparseMatrix(i = match(df_long$Sender, row_levels), 
                                        j = match(df_long$Receiver, col_levels), 
                                        x = df_long$value,
                                        dimnames = list(row_levels, col_levels))     
    }
    agg_mat[is.na(agg_mat)] <- 0
    return(agg_mat)
}

Process_GRN <- function(GRN_original, beta_thres = 0.005){

    GRN_original$beta <- as.matrix(GRN_original$beta)
    GRN_I <- GRN_original$beta %>% tidyfst::mat_df()
    colnames(GRN_I) <- c("Gene","Cell","beta")

    GRN_filtered <- GRN_I[abs(GRN_I$beta) >= beta_thres*max(abs(GRN_I$beta)),]
    GRN_filtered$Gene <- as.character(GRN_filtered$Gene)
    GRN_filtered$Receptor <- sub("=.*", "", GRN_filtered$Gene)
    GRN_filtered$SSC <- sub(".*=", "", GRN_filtered$Gene)
    return(GRN_filtered)
}

Aggr_XT_results <- function(agg_df, GRN_filtered, threshold = 1e-4){
    # agg_df: CCC df
    colnames(agg_df) <- c("Cell","LR","Value")
    agg_df <- filter(agg_df, Value > threshold)
    agg_df$Ligand <- sub("=.*", "", agg_df$LR)
    agg_df$Receptor <- sub(".*=", "", agg_df$LR)
    agg_df$comb <- paste0(agg_df$Receptor, "_", agg_df$Cell)
    GRN_filtered$comb <- paste0(GRN_filtered$Receptor, "_", GRN_filtered$Cell)
    merged_df <- merge(agg_df, GRN_filtered, by = "comb")
    merged_df <- merged_df[,c("Cell.x","Ligand","Receptor.x","SSC","Value","beta")]
    colnames(merged_df) <- c("Cell","Ligand","Receptor","SSC","Value","beta")
    merged_df <- merged_df[!duplicated(merged_df),]
    merged_df$Weight <- abs(merged_df$Value*merged_df$beta)

    return(merged_df)
}

Aggr_XT_senders <- function(CCC_origin, GRN_filtered, meta, target = NULL){

    if(!is.null(target)){
        GRN_target <- filter(GRN_filtered, SSC %in% target)
    }else GRN_target <- GRN_filtered
    XT_senders <- list()
    for(j in 1:length(CCC_origin)){
        x <- CCC_origin[[j]]
        temp <- colnames(x)
        if(!is.null(temp) && sum(x) > 0){
            temp <- strsplit(temp,split = "=")
            temp <- purrr::map(temp, 1) %>% unlist()
            sample_type <- meta[temp,"annotation"]
            type_levels <- unique(sample_type)
            type_mat <- Matrix::Matrix(0, nrow = length(sample_type), ncol = length(type_levels),
                            dimnames = list(NULL, type_levels), sparse = TRUE)

            for (i in seq_along(sample_type)) {
                type_mat[i, sample_type[i]] <- 1
            }
            agg_mat <- x %*% type_mat    
            agg_df <- tidyfst::mat_df(as.matrix(agg_mat))
            agg_df <- filter(agg_df, value > 0.01*max(agg_df$value))
            colnames(agg_df) <- c("LR","Type","value")
            agg_df$receiver <- names(CCC_origin)[j]
        }else{
            agg_df <- data.frame(LR = character(), Type = character(), value = numeric(),receiver = character())
        }
        XT_senders[[j]] <- agg_df
    }
    XT_senders <- do.call(rbind, XT_senders)

    XT_sender_aggr <- aggregate(XT_senders$value, list(XT_senders$LR, XT_senders$Type,XT_senders$receiver), sum)
    colnames(XT_sender_aggr) <- c("LR","Type","Cell","Value1")
    XT_sender_aggr$LR <- as.character(XT_sender_aggr$LR)
    XT_sender_aggr$Ligand <- sapply(XT_sender_aggr$LR, function(x){
        temp <- strsplit(x,split = "=",fixed = T)
        unlist(temp)[[1]]
    })
    XT_sender_aggr$Receptor <- sapply(XT_sender_aggr$LR, function(x){
        temp <- strsplit(x,split = "=",fixed = T)
        unlist(temp)[[2]]
    })

    GRN_target_use <- GRN_target[,c("Receptor","Cell","SSC","beta")]
    CC_results2 <- inner_join(XT_sender_aggr, GRN_target_use, by = c("Receptor","Cell"))
    CC_results2 <- distinct(CC_results2)

    results <- data.frame(Ligand = CC_results2$Ligand, Cell = CC_results2$Cell, Type = CC_results2$Type, TF = CC_results2$SSC, Weight = CC_results2$Value1*CC_results2$beta)

    return(results)
}

Pathway_Sender <- function(CCC_origin, GRN_filtered, pathway_names){

    pathway_names <- intersect(pathway_names, GRN_filtered$Gene)
    
    results_all <- lapply(pathway_names, function(x){

        LRT <- strsplit(x, split = "=") %>% unlist()
        LR <- paste0(LRT[1],"=",LRT[2])
        RT <- paste0(LRT[2],"=",LRT[3])
        results_list <- lapply(CCC_origin, function(x){
            df <- data.frame(SR = colnames(x), Strength = x[LR, ] %>% as.vector())
        })
        results_df1 <- do.call(rbind, results_df)
        results_df1$Sender <- sub("=.*", "", results_df$SR)
        results_df1$Cell <- sub(".*=", "", results_df$SR)
        results_df1$Ligand <- LRT[1]
        results_df1$Receptor <- LRT[2]
        results_df2 <- filter(GRN_filtered, Gene == x)

        results_merged <- merge(results_df1,results_df2, by = c("Receptor","Cell"))
        results_merged$Weight <- results_merged$beta * results_merged$SR
        results_merged <- results_merged[,c("Sender","Cell","Ligand","Receptor","SSC","Weight")]
        results_merged
    })

    results_all <- do.call(rbind, results_all)
    return(results_all)

}

Group_signals <- function(results, meta, features = NULL, groups = NULL, cluster_k = 3){

    # Format of results: part of the XT_results
    # Format of group: a list of cells grouped by the cluster labels
    if(is.null(features)) features <- results$SSC %>% unique()
    results_filtered <- filter(results, SSC %in% features)
    cells <- unique(results_filtered$Cell) %>% as.character()

    # Perform the k means clustering
    if(is.null(groups) || (!is.list(groups))){
        meta_used <- meta[cells,c("x","y")]
        km_res <- kmeans(meta_used, centers = cluster_k)
        groups <- split(rownames(meta_used), km_res$cluster)
    }
    group_names <- names(groups) <- paste0("Group", 1:cluster_k)
    mylist <- lapply(group_names, function(x){
        cells <- groups[[x]]
        temp <- results_filtered[results_filtered$Cell %in% cells, ]
        signal_s <- aggregate(temp$Weight, list(temp$Ligand), sum)
        colnames(signal_s) <- c("Signal","Value")
        target_s <- aggregate(temp$Weight, list(temp$SSC), sum)
        colnames(target_s) <- c("Target","Value")

        list(Cells = cells, Signals = signal_s, Targets = target_s)
    })
    names(mylist) <- group_names
    return(mylist)

}

Path_To_Pair_GRN <- function(Pathway_name, DB){

    DB_filtered <- filter(DB, pathway %in% Pathway_name)
    RTs <- paste(DB_filtered$receptor,DB_filtered$tf, sep = "=") %>% unique()

    return(RTs)
}

Path_To_Pair_CCC <- function(Pathway_name, DB){

    if(!"LR" %in% colnames(DB)){
        DB$LR <- paste(DB$From,DB_filtered$To, sep = "=")
    }
    Pathway_name <- unique(Pathway_name)
    # Find the LR pairs that matches the pathway name
    LRP1 <- Pathway_name[Pathway_name %in% DB$LR] 
    pathways <- Pathway_name[Pathway_name %in% DB$Pathway] %>% unique()

    if(length(pathways)){
        LRP2 <- filter(DB, Pathway %in% pathways) %>% pull(LR)
    }
    LRs <- c(LRP1,LRP2) %>% unique() %>% sort()
    return(LRs)
}

Calculate_Coverage <- function(XT_results, pathway1, pathway2){

    XT_results$LR <- paste(XT_results$Ligand, XT_results$Receptor, sep = "=")
    XT1 <- filter(XT_results, LR %in% pathway1)
    XT1 <- aggregate(cbind(Value, beta, Weight) ~ Cell + SSC, data = XT1, FUN = sum)
    XT2 <- filter(XT_results, LR %in% pathway2)
    XT2 <- aggregate(cbind(Value, beta, Weight) ~ Cell + SSC, data = XT2, FUN = sum)

    LR1_active_cells <- XT1$Cell %>% unique()
    LR2_active_cells <- XT2$Cell %>% unique()
    crosstalk_cells <- intersect(LR1_active_cells, LR2_active_cells)
    if(!length(crosstalk_cells)) stop("There is no crosstalk between two pathways!")

    TF_used <- intersect(XT1$SSC, XT2$SSC) %>% sort()
    cov_list <- lapply(TF_used, function(x){
        temp1 <- filter(XT1, SSC == x)
        temp2 <- filter(XT2, SSC == x)
        colnames(temp2) <- c("Cell","SSC","Value2","beta2","Weight2")
        temp1$Path <- "Pathway1"; temp2$Path <- "Pathway2"
        # compare the overlap 
        crosstalk_cells_used <- intersect(temp1$Cell, temp2$Cell) %>% intersect(crosstalk_cells)

        # ---------Need to fix this in the future----------
        strength_df <- temp1 %>%
            full_join(temp2, by = "Cell") %>%
            mutate(
                beta  = coalesce(beta, 0),
                beta2  = coalesce(beta2, 0),
                combined_GRN = beta * beta2
            ) %>%
            arrange(Cell)
        data.frame(
            Target = x,
            Coverage = length(crosstalk_cells_used)/length(crosstalk_cells),
            Strength = sum(strength_df$combined_GRN)          
        )
    })
    cov_list <- do.call(rbind, cov_list)
    cov_list <- cov_list %>% arrange(desc(Coverage))
    return(cov_list)
}

Dist_to_border <- function(meta){
    
    dist_mat <- dist(meta[,c("x","y")]) %>% as.matrix()
    result <- lapply(rownames(dist_mat),function(x){
        hetero_cells <- meta[,"annotation"] != meta[x,"annotation"]
        min(dist_mat[x, hetero_cells])
    }) %>% unlist()
    names(result) <- rownames(dist_mat)
    return(result)
    
}

Sender_contribution <- function(CCC_LR, meta){

    annos <- meta[rownames(CCC_LR),"annotation"] %>% as.vector()
    result <- rowsum(as.matrix(CCC_LR), group = annos) %>% t()
    result <- result/rowSums(result)
    return(result)
    
}

Group_targets <- function(XT_results, feature, cluster_k = 4){

    result <- aggregate(Weight ~ Cell + SSC, data = XT_results, sum) 
    result <- result %>%
        group_by(Cell) %>%
        mutate(Spe = Weight/sum(Weight)) %>%
        ungroup()

    if(is.null(feature) | length(feature) > 1) stop("Should input one feature!")
    
    results_filtered <- filter(result, SSC == feature)
    df <- data.frame(Value = results_filtered$Spe, row.names = results_filtered$Cell)
    km_res <- kmeans(df, centers = cluster_k)
    groups <- split(rownames(df), km_res$cluster)
    group_names <- names(groups) <- paste0("Group", 1:cluster_k)
    mylist <- lapply(group_names, function(x){
        cells <- groups[[x]]
        temp <- results_filtered[results_filtered$Cell %in% cells, ]
        list(Cells = cells, min_spe = min(temp$Spe),max_spe = max(temp$Spe))
    })

    names(mylist) <- group_names
    return(mylist)

}

Regression_target <- function(XT_results, group_results, exp, feature, top_k = 10){

    group_names <- names(group_results)
    XT_filtered <- XT_results[XT_results$SSC == feature, ]
    temp <- aggregate(Spe ~ Receptor, data = XT_filtered, mean)
    top_k <- min(nrow(temp), top_k)
    receptors <- temp %>% arrange(Spe) %>% slice_tail(n = top_k) %>% pull(Receptor)
    y_all <- exp[feature, ]
    x_all <- exp[receptors, ]
    result <- lapply(group_results, function(x){
        mymodel <- tryCatch({
            cell_used <- x$Cells
            y <- y_all[cell_used]
            x <- x_all[, cell_used]
            cvfit <- glmnet::cv.glmnet(t(x), y, alpha = 0.5,nfolds = 3,type.measure = "mse")
            best_lambda <- cvfit$lambda.min
            model <- glmnet::glmnet(t(x), y, alpha = 0.5,lambda = best_lambda)
            model
        }, error = function(e) NULL)
        mymodel
    })
    names(result) <- group_names
    return(result)
}
