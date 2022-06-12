#' @export
#' @noRd
validate_input_file <- function(input_file) {
  if(!file.exists(input_file)) stop(paste0("File ", input_file, "does not exist.")) 
  x <- readLines(input_file)
  if(!(grepl("^>", x[1]))) stop("There is no '>' at the beginning of the file. Input file has to be in a FASTA format.")
  if(length(x) < 1) stop("Input file is empty.")
  if(length(x) < 2) stop("Input file has to be in a FASTA format.")
  if(sum(grepl("^>", x)) == 0) stop("Input file has to be in a FASTA format.")
}

#' @importFrom stats setNames
#' @export
#' @noRd
add_missing_genes <- function(results, type = "all") {
  if(type == "grouped") {
    missing <- unique(adhesins_df_grouped[["Gene"]])[which(!(unique(adhesins_df_grouped[["Gene"]]) %in% colnames(results)))]
  } else {
    missing <- unique(adhesins_df[["Gene"]])[which(!(unique(adhesins_df[["Gene"]]) %in% colnames(results)))]
  }
  if(length(missing) > 0) {
    cbind(results, setNames(lapply(missing, function(x) x = 0), missing))
  } else {
    results
  }
}

#' @importFrom ggplot2 theme_bw theme element_blank element_text
#' @export
#' @noRd
plot_theme <- function() {
  theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom")
}



#' @importFrom ggplot2 ggsave theme
#' @export
#' @noRd
generate_report_files <- function(presence_table, elements = c("summary_table", "summary_plot", 
                                                               "presence_table", "presence_plot"), 
                                  outdir = ".", hide_absent_genes = FALSE, hide_absent_systems = FALSE, 
                                  presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  if("summary_table" %in% elements) {
    summary_table <- get_summary_table(presence_table, hide_absent = hide_absent_systems)
    write.csv(summary_table, paste0(outdir, "/summary_table.csv"), row.names = FALSE)
  }
  
  if("presence_table" %in% elements) {
    if(hide_absent_genes == TRUE) {
      pres_table <- presence_table[, c(TRUE, colSums(presence_table[, 2:ncol(presence_table)]) > 0)]
    } else {
      pres_table <- presence_table
    }
    write.csv(pres_table, paste0(outdir, "/presence_table.csv"), row.names = FALSE)
  }
  
  if("summary_plot" %in% elements) {
    summary_plot <- get_summary_plot(presence_table, hide_absent = hide_absent_systems,
                                     presence_col = presence_col, absence_col = absence_col)
    t <- get_summary_table(presence_table, hide_absent = hide_absent_systems)
    ggsave(paste0(outdir, "/summary_plot.png"), summary_plot, 
           width = 80+5*ncol(t), height = 60 + 5*nrow(t), units = "mm", limitsize = FALSE)
  }
  
  if("presence_plot" %in% elements) {
    if(hide_absent_genes == TRUE) {
      plot_dat <- presence_table[, c(TRUE, colSums(presence_table[, 2:ncol(presence_table)]) > 0)]
    } else {
      plot_dat <- presence_table
    }
    presence_plot <- get_presence_plot(plot_dat, presence_col = presence_col, absence_col = absence_col) +
      theme(legend.position = "right")
    ggsave(paste0(outdir, "/presence_plot.png"), presence_plot,
           width = 50+5*nrow(plot_dat), height = 80 + 3*ncol(plot_dat), units = "mm", limitsize = FALSE)
  }
}

#' @export
#' @noRd
pathotype_colors <- c("aEPEC" = "#45e495", "STEC" = "#ca45e4", "NA" = "#949494", "EAEC" = "#e2ab35", "EHEC" = "#e44444", "NMEC" = "#e44496", 
                      "DAEC" = "#e46f44", "Nonpathogenic" = "#44b7e4", "ETEC" = "#8f44e4", "UPEC" = "#e1e444", "tEPEC" = "#9de444", "EIEC" = "#4471e4")


#' @importFrom dplyr select mutate bind_rows
#' @export
#' @noRd
get_clustering_plot_data <- function(presence_table) {
  pred_res <- as.data.frame(predict(UMAP_data[["umap"]], select(presence_table, -"File")))
  colnames(pred_res) <- c("dim1", "dim2")
  dat <- bind_rows(data.frame(dim1 = UMAP_data[["umap"]][["layout"]][,1],
                              dim2 = UMAP_data[["umap"]][["layout"]][,2],
                              pathotype = UMAP_data[["labels"]],
                              label = NA),
                   mutate(pred_res,
                          pathotype = "new",
                          label = presence_table[["File"]])) 
  mutate(dat,
         type = case_when(pathotype %in% c("UPEC", "NMEC") ~ "InPEC",
                          pathotype %in% c("aEPEC", "DAEC", "EAEC", "EHEC", "EIEC", "ETEC", "STEC", "tEPEC") ~ "ExPEC",
                          pathotype == "Nonpathogenic" ~ "Nonpathogenic",
                          pathotype == "NA" ~ "NA",
                          pathotype == "new" ~ "new"))
}

#' @importFrom ggplot2 ggplot aes geom_point theme_bw scale_size_manual scale_color_manual
#' @importFrom ggrepel geom_label_repel
#' @export
#' @noRd
plot_clustering <- function(clustering_plot_dat, show_labels = TRUE, pathotypes = "detailed") {
  
  color_by <- ifelse(pathotypes == "detailed", "pathotype", "type")
  
  p <- ggplot(clustering_plot_dat, aes(x = dim1, y = dim2, color = get(color_by), label = label)) +
    geom_point(aes(size = pathotype == "new")) +
    theme_bw() +
    scale_size_manual(values = c(1, 2), guide = "none")
  
  p <- if(pathotypes == "detailed") {
    p +
      scale_color_manual("Pathotype", values = c(pathotype_colors, "new" = "black"),
                         breaks = sort(names(pathotype_colors)))
  } else if(pathotypes == "intestinal") {
    p +
      scale_color_manual("Pathotype", values = c("InPEC" = "#e1e444", "NA" = "#949494", "ExPEC" = "#e44444", "Nonpathogenic" = "#44b7e4",
                                                 "new" = "black"), breaks = c("InPEC", "ExPEC", "Nonpathogenic", "NA"))
  }
  
  if(show_labels == TRUE) {
    p +
      geom_label_repel(box.padding = 1.5,
                       show.legend = FALSE)
  } else{
    p
  }
}

#' @importFrom dplyr between mutate
#' @export
#' @noRd
check_locations <- function(x) {
  max_row <- which(x[["Alignment length"]] == max(x[["Alignment length"]]))[1]
  start <- x[["Query start"]][max_row]-round(0.1*x[["Alignment length"]][max_row])
  end <- x[["Query end"]][max_row]+round(0.1*x[["Alignment length"]][max_row])
  mutate(x, 
         same_location = ifelse(between(`Query start`, start, end) | between(`Query end`, start, end), TRUE, FALSE))
}

#' Get presence from BLAST
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
get_presence_from_blast <- function(blast_res, add_missing = TRUE, identity_threshold = 75, evalue_threshold = 1e-100) {
  gene_groups <- adhesiomeR::gene_groups
  res <- mutate(
    summarise(
      group_by(blast_res, File, Subject),
      Presence = ifelse(any(`% identity` > identity_threshold & Evalue < evalue_threshold), 
                        `% identity`, 0)),
    Gene = Subject)
  pivoted_res <- pivot_wider(res[, c("File", "Gene", "Presence")],
                             names_from = Gene, values_from = Presence, values_fill = 0)
  if("NA" %in% colnames(pivoted_res)) {
    pivoted_res <- pivoted_res[, which(colnames(pivoted_res) != "NA")]
  }
  # Check for genes that were not found
  pivoted_res <- ungroup(add_missing_genes(pivoted_res))
  
  # Decide between similar genes
  y <- pivoted_res
  for (i in problematic_genes) {
    y <- data.frame(y[, which(!(colnames(y) %in% i))],
                    t(apply(y[, which(colnames(y) %in% i)], 1, 
                            function(x) replace(x, x != max(x, na.rm = TRUE), 0))),
                    check.names = FALSE)
  }
  y <- mutate(y, across(2:ncol(y), function(x) ifelse(x > 0, 1, 0)))
  
  # Group the most similar genes
  group_res <- do.call(cbind, lapply(names(gene_groups), function(ith_group) {
    x <- select(y, gene_groups[[ith_group]])
    setNames(data.frame(group = ifelse(rowSums(x) > 0, 1, 0)), ith_group)
  }))
  updated_res <- cbind(y[, colnames(y)[which(!(colnames(y) %in% unlist(unname(gene_groups))))]],
                       group_res)
  if(add_missing == FALSE) {
    updated_res <- updated_res[, c(TRUE, colSums(updated_res[, 2:ncol(updated_res)]) > 0)]
  }
  updated_res
}
