# A quick start for SpaXTalk

## STEP 1: Load required packages
```r
library(SpaXTalk)
library(dplyr)
```

## STEP 2: Load the spatial transcriptomics data
An example dataset of mouse embryo Stereo-seq data is available [here](https://drive.google.com/file/d/1Xo40BXmw7u1EVjvJ2RlYh30Xrky8mSAF/view?usp=drive_link). 

This dataset is provided as a Seurat Object. For more information about the Seurat R package, please visit [here](https://satijalab.org/seurat/articles/get_started_v5_new).

You can now load the data using the following code:
```r
SeuratObj <- readRDS("path/to/SpaXTalk_SampleData.RDS")
```

## STEP 3: Prepare the input of SpaXTalk
The gene expression matrix (Gene-by-Cell), and cell meta data (x coordinates, y coordinates, and cell annotation) are required as the inputs for SpaXTalk analysis

By default, we assume the coordinates are stored in the "x" and "y" columns of the meta data, cell annotation is stored in the "annotation" column of the meta data. If not, you may use the column that stores them.
```r
exp_data <- SeuratObject::GetAssayData(
    SeuratObj,
    assay = "RNA",
    layer = "data"
)
meta_all <- SeuratObj@meta.data
meta <- data.frame(
    x = meta_all$x,
    y = meta_all$y,
    annotation = as.character(meta_all$annotation),
    row.names = rownames(meta_all),
    stringsAsFactors = FALSE
)

# Specify the receiver population of interest.
receiver <- rownames(meta)[meta$annotation %in% "Neural crest"]
```
If you use a Seurat Object to store the data, you may alternatively use the following command to generate all the input data you need.
```r
prepared <- Prepare_data(
    SeuratObj,
    receiver_cluster = "Neural crest",
    assay_type = "RNA",
    layer_type = "data"
)
exp_data <- prepared$exp_data
meta <- prepared$meta
receiver <- prepared$receiver
```
`Prepare_data()` returns the expression matrix, standardized metadata, and receiver names in a named list.

## STEP 4: Load the prior knowledge
SpaXTalk offers the built-in cell-cell communication database extracted from CellChat, and the intracellular gene-gene regulation database from KEGG.
The following code will import the two knowledge bases into the current workspace.
```r
data("CellChat_mouse", package = "SpaXTalk") # loads LRDB
data("KEGG_PPR_mouse", package = "SpaXTalk") # loads RTDB
```

## STEP 5: Perform the SpaXTalk analysis
```r
neighbor_list <- Identify_neighbor_list(
    locs = meta[, c("x", "y"), drop = FALSE],
    receivers = receiver,
    dist_threshold = 3,
    decay_type = "exp",
    normalize = TRUE,
    nclusters = 1
)
neighbor_list[[1]] # Check what the neighborhood looks like
CCC_origin <- InferCCC(
    data = exp_data,
    neighbor_list = neighbor_list,
    LRDB = LRDB,
    nclusters = 1
)
sum(CCC_origin[[1]])

GRNs <- InferGRN(
    exps = exp_data,
    neighbor_list = neighbor_list,
    RTDB = RTDB,
    PPR_top = 15,
    exp_thres = 0.05,
    nsim = 0,
    lambda = 1,
    nclusters = 1,
    use_delta_y = TRUE
)
betas <- GRNs$beta

# Save the results (optional)
allresults <- list(
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "SpaXTalk_results.rds")
```

## STEP 6: Process the results
```r
# Calculate the crosstalk
CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # receiver-by-LR
CCC_df <- tidyfst::mat_df(CCC_mat)
GRN_filtered <- Process_GRN(allresults$GRN, beta_thres = 1e-5)
XT_results <- Aggr_XT_results(CCC_df, GRN_filtered, threshold = 1e-8)
```

## STEP 7: Visualize the results
### 1. Plot the cell-cell communication
```r
# Plot the cluster-level communication pattern
Plot_DominantLR(results = CCC_mat)
```
We can also plot communication in space
```r
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = NULL)
PlotSpa_SRVector(meta = meta, data = data, threshold = 1)

# Plot a single LR pair
features <- c("Mif=Cd44", "Mdk=Ncl", "Mdk=Itgb1", "Igf2=Igf2r")
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[1])
PlotSpa_SRVector(meta = meta, data = data, threshold = 0.75)
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[2])
PlotSpa_SRVector(meta = meta, data = data, threshold = 0.75)
```

### 2. Plot the crosstalk
Crosstalk can be visualized at single-cell, or cell-group levels!
Plot crosstalk pattern of one individual cell
```r
color_use <- c(Ligand = "aquamarine", Receptor = "bisque", Target = "pink")
XT_filtered <- XT_results %>%
    filter(Cell == "185_151") %>%
    arrange(desc(Weight)) %>%
    slice_head(n = 6)
Plot_MultiLayer_Singlecell(
    XT_filtered,
    start_trim = 0.3,
    end_trim = 0.3,
    arrow_len_cm = 0.25,
    x_scale = 3.25,
    rect_width = 3,
    rect_height = 0.3,
    label_size = 6,
    layer_fill_colors = color_use
)
```
Plot crosstalk of the whole cell group
```r
aggr_XT <- aggregate(Weight ~ Ligand + Receptor + TF, data = XT_results, sum)
PlotXT_Alluvial(aggr_XT, features = NULL, topk = 10, vertical = FALSE)
```
Plot the spatial distribution of the specificity of a single gene
```r
tf <- XT_results$TF[which.max(XT_results$Weight)]
aggr_spe <- XT_results %>%
    group_by(Cell, TF) %>%
    summarise(Weight = sum(Weight), .groups = "drop_last") %>%
    mutate(Spe = Weight / sum(Weight)) %>%
    ungroup() %>%
    filter(TF == tf)

specificity <- setNames(aggr_spe$Spe, aggr_spe$Cell)
PlotSpa_Scatter(specificity, meta, legend_name = paste(tf, "specificity"))
```
