#' Get presence table
#' 
#' This function creates a presence/absence table of genes in analysed genomes
#' based on results of the BLAST search. Each row corresponds to one input file
#' and each column to one gene. Presence of a given gene is indicated as a 1, 
#' whereas absence as a 0. You can modify thresholds used to consider a gene
#' as present using \code{identity_threshold} and \code{evalue_threshold} 
#' arguments. 
#' @param blast_res blast results obtained with \code{\link{get_blast_res}}
#' @param add_missing \code{logical} indicating if genes not found by BLAST 
#' should be added. By default \code{TRUE}, meaning that all genes are shown 
#' in the resulting table, even if they were not found in any genome. 
#' @param count_copies \code{logical} indicating if occurences of gene 
#' copies should be counted. An occurence of gene is considered a separate
#' copy if its location do not overlap with other hit to the same gene.
#' @param coeff coefficient from the linear model used to set the bit score
#' threshold
#' @param interc intercept from the linear model used to set the bit score
#' threshold
#' @param n_threads number of threads for parallel processing. Default is one. 
#' The maximum number of threads is determined by \code{\link[parallel]{detectCores}}.
#' @return a data frame of gene presence/absence. The first column contains
#' the names of input files and the following correspond to analysed genes.
#' Presence of a gene is indicated by 1, whereas absence by 0. In case of
#' copy counts, numbers of gene copies is presented.
#' @importFrom dplyr group_by summarise mutate filter ungroup bind_rows select
#' @importFrom tidyr pivot_wider
#' @importFrom stats aggregate setNames
#' @importFrom future.apply future_apply
#' @importFrom future makeClusterPSOCK plan multisession sequential
#' @importFrom parallel stopCluster detectCores
#' @importFrom progressr with_progress progressor handlers handler_progress
#' @export
get_presence_table <- function(blast_res, add_missing = TRUE, count_copies = FALSE, coeff, interc, n_threads = 1) {
  max_nt <- detectCores(logical = FALSE)
  if(n_threads > max_nt) {
    stop(paste0("The number of threads you specified is too large. The maximum number of threads determined by parallel::detectCores function is: ", detectCores(logical = FALSE), ". 
  Please select a value between 1 and ", detectCores(logical = FALSE), "."))
  } else if (!(n_threads %in% 1L:detectCores(logical = FALSE)))  {
    stop("The number of threads is incorrect. Please make sure that you entered a valid number.")
  }
  
  problematic_genes <- adhesiomeR::problematic_genes
  nonproblematic_genes <- adhesiomeR::adhesins_df[["Gene"]][which(!(adhesiomeR::adhesins_df[["Gene"]] %in% unlist(problematic_genes)))]
  #len_groups <- adhesiomeR::len_groups
  #gene_groups <- adhesiomeR::gene_groups
  adhesins_lengths <- adhesiomeR::adhesins_lengths
  
  # short_res <- filter(blast_res, Subject %in% len_groups[["short"]], Evalue < 10^-35, `% identity` > identity)
  # medium_res <- filter(blast_res, Subject %in% len_groups[["medium"]], Evalue < 10^-70, `% identity` > identity)
  # long_res <- filter(blast_res, Subject %in% len_groups[["long"]], Evalue < 10^-100, `% identity` > identity)
  # all_blast_res <- bind_rows(bind_rows(short_res, medium_res), long_res)
  all_blast_res <- filter(left_join(blast_res, adhesins_lengths, by = c("Subject" = "Gene")), `Bit score` >= coeff*Length-interc)
  
  if(n_threads > 1) {
    plan(multisession, workers = n_threads, gc = TRUE)
    
    all_res <- bind_rows(
      future_lapply(problematic_genes, function(ith_set) {
        get_gene_presence_for_localizations(all_blast_res, "problematic", ith_set)
      }),
      future_lapply(nonproblematic_genes, function(ith_gene) {
        get_gene_presence_for_localizations(all_blast_res, "nonproblematic", ith_gene)
      })
    ) 
    plan(sequential)
  } else {
    all_res <- bind_rows(
      lapply(problematic_genes, function(ith_set) {
        get_gene_presence_for_localizations(all_blast_res, "problematic", ith_set)
      }),
      lapply(nonproblematic_genes, function(ith_gene) {
        get_gene_presence_for_localizations(all_blast_res, "nonproblematic", ith_gene)
      })
    ) 
  }
  
  # Group the most similar genes
  group_res <- do.call(cbind, lapply(names(gene_groups), function(ith_group) {
    x <- select(all_res, gene_groups[[ith_group]])
    setNames(data.frame(group = ifelse(rowSums(x) > 0, rowSums(x), 0)), ith_group)
  }))
  updated_res <- cbind(all_res[, colnames(all_res)[which(!(colnames(all_res) %in% unlist(unname(gene_groups))))]],
                       group_res)
  
  if(add_missing == FALSE) {
    updated_res <- updated_res[, c(TRUE, colSums(updated_res[, 2:ncol(updated_res)]) > 0)]
  }
  
  # Aggregate presence or copy counts
  updated_res[["File"]] <- as.factor(updated_res[["File"]])
  aggregated_res <- if(count_copies == FALSE) {
    aggregate(. ~ File, data = updated_res, FUN = function(i) ifelse(sum(i) > 0, 1, 0))
  } else {
    aggregate(. ~ File, data = updated_res, FUN = sum)
  }
  final_res <- mutate(aggregated_res, File = as.character(File))
  
  # Check for files without found genes
  files_to_add <- unique(blast_res[["File"]])[which(!(unique(blast_res[["File"]]) %in% unique(final_res[["File"]])))]
  if(length(files_to_add) > 0) {
    rbind(final_res, 
          setNames(
          cbind(data.frame("File" = files_to_add),
                data.frame(matrix(0, nrow = length(files_to_add), ncol = ncol(final_res)-1))),
          colnames(final_res)))
  } else {
    final_res
  }
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
get_presence_plot_data <- function(presence_table, systems = unique(adhesins_df_grouped[["System"]])) {
  
  selected_genes <- unique(filter(adhesins_df_grouped, System %in% systems)[["Gene"]])
  
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

