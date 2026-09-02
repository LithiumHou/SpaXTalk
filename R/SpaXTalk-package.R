#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' SpaXTalkR: Quantification of Signal Crosstalk using spatial data
#'
#' SpaXTalkR provides tools for quantifying downstream effects of
#' cell-cell communication signals and how they vary across space
#'
#' @import dplyr ggplot2
#' @importFrom circlize chordDiagram circos.clear circos.text circos.track
#'   get.cell.meta.data
#' @importFrom ComplexHeatmap anno_boxplot HeatmapAnnotation rowAnnotation
#' @importFrom forcats fct_reorder
#' @importFrom ggraph create_layout geom_edge_arc geom_node_point geom_node_text
#'   ggraph
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow gpar unit
#' @importFrom igraph graph_from_data_frame
#' @importFrom progressr progressor with_progress
#' @importFrom rlang .data
#' @importFrom scatterpie geom_scatterpie
#' @importFrom SeuratObject GetAssayData
#' @importFrom stats aggregate dist kmeans median p.adjust quantile setNames
#'   t.test wilcox.test
#' @importFrom tibble tibble
#' @keywords internal
"_PACKAGE"
## usethis namespace: end
NULL
