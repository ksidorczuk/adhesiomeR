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
add_missing_genes <- function(results) {
  missing <- unique(adhesins_df[["Gene"]])[which(!(unique(adhesins_df[["Gene"]]) %in% colnames(results)))]
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
