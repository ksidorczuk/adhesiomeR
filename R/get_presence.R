#' Get presence table
#' 
#' This function creates a presence/absence table of genes in analysed genomes
#' based on results of the BLAST search. Each row corresponds to one input file
#' and each column to one gene. Presence of a given gene is indicated as a 1, 
#' whereas absence as a 0. You can modify thresholds used to consider a gene
#' as present using \code{identity_threshold} and \code{evalue_threshold} 
#' arguments. By default, gene is considered to be present when it shares
#' over 70% of identity with a subject sequence and has E-value lower than
#' 1e-50. 
#' @param blast_res blast results obtained with \code{\link{get_blast_res}}
#' @param add_missing \code{logical} indicating if genes not found by BLAST 
#' should be added. By default \code{TRUE}, meaning that all genes are shown 
#' in the resulting table, even if they were not found in any genome. 
#' @param identity_threshold \code{numeric} indicating the percent of identity
#' used for labeling a gene as present or absent
#' @param evalue_threshold \code{numeric} indicating the E-value threshold
#' used for labeling a gene as present or absent
#' @return a data frame of gene presence/absence. The first column contains
#' the names of input files and the following correspond to analysed genes.
#' Presence of a gene is indicated by 1, whereas absence by 0. 
#' @importFrom dplyr group_by summarise mutate filter ungroup
#' @importFrom tidyr pivot_wider
#' @export
get_presence_table <- function(blast_res, add_missing = TRUE, identity_threshold = 75, evalue_threshold = 1e-100) {
  gene_groups <- adhesiomeR::gene_groups
  res <- mutate(
    summarise(
      group_by(blast_res, File, Subject),
      Presence = ifelse(any(`% identity` > identity_threshold & Evalue < evalue_threshold), 1, 0)),
    Gene = Subject)
  pivoted_res <- pivot_wider(res[, c("File", "Gene", "Presence")],
                             names_from = Gene, values_from = Presence, values_fill = 0)
  if("NA" %in% colnames(pivoted_res)) {
    pivoted_res <- pivoted_res[, which(colnames(pivoted_res) != "NA")]
  }
  # Check for genes that were not found
  pivoted_res <- ungroup(add_missing_genes(pivoted_res))
  # Group similar genes
  group_res <- do.call(cbind, lapply(names(gene_groups), function(ith_group) {
    x <- select(pivoted_res, gene_groups[[ith_group]])
    setNames(data.frame(group = ifelse(rowSums(x) > 0, 1, 0)), ith_group)
  }))
  updated_res <- cbind(pivoted_res[, colnames(pivoted_res)[which(!(colnames(pivoted_res) %in% unlist(unname(gene_groups))))]],
                       group_res)
  if(add_missing == FALSE) {
    updated_res <- updated_res[, c(TRUE, colSums(updated_res[, 2:ncol(updated_res)]) > 0)]
  }
  updated_res
}


#' @importFrom stats as.dendrogram hclust dist order.dendrogram
cluster_data <- function(df, data_to_cluster, var_name) {
  dendro_files <- as.dendrogram(hclust(d = dist(x = as.matrix(data_to_cluster[, 2:ncol(data_to_cluster)]))))
  files_order <- order.dendrogram(dendro_files)
  dendro_var <- as.dendrogram(hclust(d = dist(t(as.matrix(data_to_cluster[, 2:ncol(data_to_cluster)])))))
  var_order <- order.dendrogram(dendro_var)
  
  df[[var_name]] <- factor(df[[var_name]], 
                           levels = colnames(data_to_cluster[2:ncol(data_to_cluster)][var_order]),
                           ordered = TRUE)
  df[["File"]] <- factor(df[["File"]],
                         levels = data_to_cluster[["File"]][files_order],
                         ordered = TRUE) 
  df
}


#' @importFrom dplyr filter ungroup
#' @importFrom tidyr pivot_longer
#' @noRd
#' @export
get_presence_plot_data <- function(presence_table, systems = unique(adhesins_df[["System"]])) {
  
  selected_genes <- unique(filter(adhesins_df, System %in% systems)[["Gene"]])
  
  plot_dat <- filter(
    pivot_longer(presence_table, 2:ncol(presence_table), names_to = "Gene", values_to = "Presence"),
    Gene %in% selected_genes)
  
  if(nrow(presence_table) > 1) {
    plot_dat <- cluster_data(plot_dat, presence_table, "Gene")
  }
  
  plot_dat[["Presence"]] <- factor(ifelse(plot_dat[["Presence"]] == 1, "yes", "no"), levels = c("yes", "no"))
  ungroup(plot_dat)
}


#' Get plot with gene presence/absence
#' 
#' This function generates a heatmap showing the presence or absence of each 
#' gene in each of the analysed files. 
#' @param presence_table a data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function
#' @param systems a character vector with names of the systems that should
#' be shown on a heatmap. By default, all systems present in the database
#' @param presence_col color of the tiles representing present genes. Must be
#' specified as a hex color code
#' @param absence_col color of the tiles representing absent genes. Must be
#' specified as a hex color code
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual scale_x_discrete theme aes
#' @export
get_presence_plot <- function(presence_table, systems = unique(adhesins_df[["System"]]),
                              presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  plot_dat <- get_presence_plot_data(presence_table, systems)
  
  ggplot(plot_dat, aes(x = File, y = Gene, fill = Presence)) +
    geom_tile() +
    scale_fill_manual("Presence", values = c("yes" = presence_col, "no" = absence_col), drop = FALSE) +
    scale_x_discrete(position = "top") +
    plot_theme() +
    theme(legend.direction = "vertical")
}

