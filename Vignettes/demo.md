# A quick start for SpaXTalk

A clean R script (no markdown instructions) of this tutorial is available through [here](./demo_pbmc.R).

## Load required packages
```
library(SpaXTalkR)
library(Seurat)
library(dplyr)
```

## Load the spatial transcriptomics data
An example dataset of mouse embryo Stereo-seq data is available [here]. 

This dataset is provided as a Seurat Object. For more information about the Seurat R package, please visit [here](https://satijalab.org/seurat/articles/get_started_v5_new).

You can now load the data using the following code:
```
SeuratObj <- readRDS("E:/Scomics/spatialXT/E95.RDS")
```

## Prepare the input of SpaXTalk
The gene expression matrix (Gene-by-Cell), and cell meta data (x coordinates, y coordinates, and cell annotation) are required as the inputs for SpaXTalk analysis

By default, we assume the coordinates are stored in the "x" and "y" columns of the meta data, cell annotation is stored in the "annotation" column of the meta data. If not, you may use the column that stores them.
```
exp_data <- SeuratObj@assays$RNA$data
meta_all <- SeuratObj@meta.data
meta <- data.frame(x = meta_all$x,
                   y = meta_all$y,
                   annotation = meta_all$annotation %>% as.character(),
                   stringsAsFactors = FALSE)
```

## Load the prior knowledge
SpaXTalk offers the built-in cell-cell communication database extracted from CellChat, and the intracellular gene-gene regulation database from KEGG.

```

```
