# A quick start for SpaXTalk

## STEP 1: Load required packages
```
library(SpaXTalkR)
library(Seurat)
library(dplyr)
```

## STEP 2: Load the spatial transcriptomics data
An example dataset of mouse embryo Stereo-seq data is available [here]. 

This dataset is provided as a Seurat Object. For more information about the Seurat R package, please visit [here](https://satijalab.org/seurat/articles/get_started_v5_new).

You can now load the data using the following code:
```
SeuratObj <- readRDS("E:/Scomics/spatialXT/E95.RDS")
```

## STEP 3: Prepare the input of SpaXTalk
The gene expression matrix (Gene-by-Cell), and cell meta data (x coordinates, y coordinates, and cell annotation) are required as the inputs for SpaXTalk analysis

By default, we assume the coordinates are stored in the "x" and "y" columns of the meta data, cell annotation is stored in the "annotation" column of the meta data. If not, you may use the column that stores them.
```
exp_data <- SeuratObj@assays$RNA$data
meta_all <- SeuratObj@meta.data
meta <- data.frame(x = meta_all$x,
                   y = meta_all$y,
                   annotation = meta_all$annotation %>% as.character(),
                   stringsAsFactors = FALSE)
receiver <- rownames(meta[meta$annotation %in% c("Neural crest"), ]) # You can specify the receivers of your interest.
```
If you use a Seurat Object to store the data, you may alternatively use the following command to generate all the input data you need.
```
receiver_cluster = "Neural crest"
Prepare_data(SeuratObj, receiver_cluster = receiver_cluster, assay_type = "RNA")
```
Please note that Prepare_data() does not have a return value. It will automatically load the required inputs as variables to the workspace.

## STEP 4: Load the prior knowledge
SpaXTalk offers the built-in cell-cell communication database extracted from CellChat, and the intracellular gene-gene regulation database from KEGG.
The following code will import the two knowledge bases into the current workspace.
```
data("CellChat_mouse")
data("KEGG_PPR_mouse")
```

## STEP 5: Perform the SpaXTalk analysis
```
neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 3,decay_type = "exp",normalize = T)
neighbor_list[[1]] # Check what the neighborhood looks like
CCC_origin <- InferCCC_Scriabin(exp_data, neighbor_list = neighbor_list, LRDB = LRDB, ligands = NULL, receptors = NULL)
(CCC_origin[[1]]) %>% sum()
# neighbor_list <- Identify_neighbor_list(locs = meta[,1:2],receivers = receiver,dist_threshold = 2,decay_type = "exp",normalize = T)
GRNs <- InferGRN(exps = exp_data, neighbor_list = neighbor_list, PPR_top = 15,exp_thres = 0.05,RTDB = RTDB, nsim = 0,lambda = 1, nclusters = 64, use_delta_y = T)
betas <- GRNs$beta

# Save the results (optional)
allresults <- list(
    CCC = CCC_origin,
    GRN = GRNs
)
saveRDS(allresults, "./results/Chen2022/E95_new.RDS")
```

## STEP 6: Process the results
```
# Calculate the crosstalk
CCC_mat <- Aggr_CCC_results(CCC_origin, flag = "Sender") # sender*LR
CCC_df <- tidyfst::mat_df(CCC_mat)
GRN_filtered <- Process_GRN(allresults$GRN, beta_thres = 1e-5)
XT_results <- Aggr_XT_results(CCC_df, GRN_filtered, threshold = 1e-8)
```

## STEP 7: Visualize the results
### 1. Plot the cell-cell communication
```
# Plot the cluster-level communication pattern
colnames(CCC_mat) <- gsub("="," - ", colnames(CCC_mat))
Plot_DominantLR(results = CCC_mat)

# If you have installed SigXTalk, you may also try the following code:
Ls <- sub(" - .*", "", colnames(CCC_mat))
Rs <- sub(".* - ", "", colnames(CCC_mat))
df <- data.frame(Ligand = Ls, Receptor = Rs, Weight = colSums(CCC_mat))
SigXTalkR::PlotCCI_CirclePlot(df, topk = 15)
```
We can also plot communication in space
```
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = NULL)
PlotSpa_SRVector(meta = meta, data = data,threshold = 1)

# Plot a single LR pair
features <- c("Mif=Cd44","Mdk=Ncl","Mdk=Itgb1","Igf2=Igf2r")
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[1])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
data <- Aggr_CCC_results(CCC_origin, flag = "LR", LR_used = features[2])
PlotSpa_SRVector(meta = meta, data = data,threshold = .75)
```

### 2. Plot the crosstalk
Crosstalk can be visualized at single-cell, or cell-group levels!
Plot crosstalk pattern of one individual cell
```
color_use <- c("aquamarine","bisque","pink")
names(color_use) <- c("Ligand","Receptor","Target")
XT_filtered <- XT_results[XT_results$Cell == "185_151", ] %>% arrange(desc(Weight)) %>% slice_head(n = 6)
Plot_MultiLayer_Singlecell(XT_filtered, start_trim = 0.3, end_trim = 0.3, arrow_len_cm = 0.25,
  x_scale = 3.25,rect_width = 3,rect_height = 0.3, label_size = 6,layer_fill_colors = color_use)
```
Plot crosstalk of the whole cell group
```
aggr_XT <- aggregate(Weight ~ Ligand+Receptor+SSC, data = XT_results, sum)
PlotXT_Alluvial(aggr_XT, features = NULL, topk = 10, vertical = F)
```
Plot the spatial distribution of the specificity of a single gene
```
XT_spe_features <- XT_results %>% 
    group_by(Cell, Ligand) %>% 
    mutate(Spe = Weight/sum(Weight)) %>%
    ungroup() %>% filter(TF %in% features)
aggr_spe <- aggregate(Spe ~ Cell + TF, data = XT_spe_features, mean)
PlotSpa_Scatter(aggr_spe, meta)
```
