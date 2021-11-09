#' Plot clustering with pathotypes
#' 
#' This function clusters analysed genomes with pathotyped genomes
#' from RefSeq database. Clustering is performed using UMAP algorithm.
#' @param presence_table data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function.
#' @param show_labels logical indicating if labels with filenames of 
#' the analysed genomes should be plotted. Default is \code{TRUE}.
#' @return A plot with the clustering results
#' @seealso get_presence_table
#' @importFrom stats predict setNames
#' @importFrom dplyr mutate bind_rows select
#' @importFrom ggplot2 aes geom_point theme_bw scale_color_manual scale_size_manual
#' @importFrom ggrepel geom_label_repel
#' @export
plot_clustering_with_pathotypes <- function(presence_table, show_labels = TRUE) {
  pred_res <- as.data.frame(predict(UMAP_data[["umap"]], select(presence_table, -"File")))
  colnames(pred_res) <- c("dim1", "dim2")
  plot_dat <- bind_rows(data.frame(dim1 = UMAP_data[["umap"]][["layout"]][,1],
                                   dim2 = UMAP_data[["umap"]][["layout"]][,2],
                                   pathotype = UMAP_data[["labels"]],
                                   label = NA),
                        mutate(pred_res,
                               pathotype = "new",
                               label = presence_table[["File"]]))
  p <- ggplot(plot_dat, aes(x = dim1, y = dim2, color = pathotype, label = label)) +
    geom_point(aes(size = pathotype == "new")) +
    theme_bw() +
    scale_color_manual("Pathotype", values = c(pathotype_colors, "new" = "black"),
                       breaks = sort(names(pathotype_colors))) +
    scale_size_manual(values = c(1, 3), guide = "none")
  if(show_labels == TRUE) {
    p +
      geom_label_repel(box.padding = 1.5,
                       show.legend = FALSE)
  } else{
    p
  }
}
