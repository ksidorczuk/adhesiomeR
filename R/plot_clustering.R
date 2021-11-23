#' Plot clustering with pathotypes
#' 
#' This function clusters analysed genomes with pathotyped genomes
#' from RefSeq database. Clustering is performed using UMAP algorithm.
#' @param presence_table data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function.
#' @param show_labels logical indicating if labels with filenames of 
#' the analysed genomes should be plotted. Default is \code{TRUE}.
#' @param pathotypes one of \code{"detailed"} or \code{"intestinal"}
#' indicating coloring scheme of the plot. By default \code{"detailed"}
#' is used, meaning that following pathotypes are indicated by colors:
#' aEPEC, tEPEC, EHEC, ETEC, EAEC, EIEC, DAEC, STEC, UPEC, NMEC, 
#' Nonpathogenic, NA (unknown). If \code{"intestinal"} is selected, 
#' colors indicate following groups: InPEC (intestinal pathogenic, 
#' includes aEPEC, tEPEC, EHEC, ETEC, EAEC, EIEC, DAEC, STEC), ExPEC
#' (extraintestinal pathogenic, includes UPEC and NMEC), Nonpathogenic
#' and NA (unknown).
#' @return A plot with the clustering results
#' @seealso get_presence_table
#' @importFrom stats predict setNames
#' @importFrom dplyr mutate bind_rows select
#' @importFrom ggplot2 aes geom_point theme_bw scale_color_manual scale_size_manual
#' @importFrom ggrepel geom_label_repel
#' @export
get_clustering_plot <- function(presence_table, show_labels = TRUE, pathotypes = "detailed") {
  if(!(pathotypes %in% c("detailed", "intestinal"))) {
    stop(paste0("Incorrect value passed as pathotypes argument. Allowed values are 'detailed' or 'intestinal'."))
  }
  plot_dat <- get_clustering_plot_data(presence_table)
  plot_clustering(plot_dat, show_labels = show_labels, pathotypes = pathotypes)
}
