#' Validate input file
#' 
#' Checks input files for fasta format.
#' 
#' @param input_file path to input file
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

#' Add missing genes
#' 
#' Adds genes that were not found in any file to the presence table
#' 
#' @param type specifies type of added genes. All adds single genes,
#' grouped considers genes with the most similar ones grouped together.
#' @importFrom stats setNames
#' @export
#' @noRd
add_missing_genes <- function(results, type = "all") {
  if(!(type) %in% c("all", "grouped")) stop("Type must be 'all' or 'grouped'!")
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

#' Check gene locations
#' 
#' Checks if genes are located in the same region by comparing start
#' and end positions of the alignments. Genes are considered to be
#' located in the same region when their start or end position
#' overlap or are located within a close proximity (length equal to
#' 10% of the longest alignment for analysed gene). 
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



#' Get gene presence for localizations
#' 
#' This function translates BLAST search results into gene presence.
#' It considers gene localization and copies. For problematic genes it 
#' performs decision between them based on the identity percent. 
#' @param blast_res blast results obtained with \code{\link{get_blast_res}}
#' @param type if \code{problematic} (default) decides between genes based
#' on the identity percent, else no selection is made.
#' @param gene_set set of genes to consider
#' @return a data frame of gene presence/absence. The first column contains
#' the names of input files and the following correspond to analysed genes.
#' Presence of a gene is indicated by 1, whereas absence by 0. 
#' @importFrom dplyr mutate filter select bind_rows
#' @export
#' @noRd
get_gene_presence_for_localizations <- function(blast_res, type = "problematic", gene_set) {
  x <- mutate(
    filter(blast_res, Subject %in% gene_set), 
    same_location = FALSE)
  full_res <- data.frame()
  bind_rows(
    lapply(unique(x[["File"]]), function(ith_file) {
      while(any(x[["same_location"]] == FALSE)) {
        y <- filter(x, File == ith_file)
        locations <- check_locations(y)
        
        if(type == "problematic") {
          # Decide between similar genes
          res <- decide_between_problematic_genes(
            get_presence_from_blast(
              select(
                filter(locations, same_location == TRUE),
                -same_location), type = "identity"),
            gene_set)
          x <- filter(locations, same_location == FALSE)
          full_res <- bind_rows(full_res, res)
          
        } else {
          # Process the rest normally
          res <- get_presence_from_blast(
            select(
              filter(locations, same_location == TRUE),
              -same_location))
          x <- filter(locations, same_location == FALSE)
          full_res <- bind_rows(full_res, res)
        }
      }
      full_res
    })
  )
}



#' Decide between problematic genes
#' 
#' This function performs decision between a group of problematic genes
#' by selecting the hit with the highest identity percent to genes from the
#' group.
#' @param pivoted_res a data frame with identities from blast search,
#' obtained by \code{\link{get_presence_from_blast}} using \code{identity}
#' as type. 
#' @param problematic_genes set of problematic genes to consider.
#' @return a data frame of gene presence/absence modified to consider
#' only the best hit (with the highest identity) to one of the genes
#' from the considered group.
#' @importFrom dplyr mutate across
#' @export
#' @noRd
decide_between_problematic_genes <- function(pivoted_res, problematic_genes) {
  y <- data.frame(pivoted_res[, which(!(colnames(pivoted_res) %in% problematic_genes))],
                  t(apply(pivoted_res[, which(colnames(pivoted_res) %in% problematic_genes)], 1, 
                          function(x) replace(x, x != max(x, na.rm = TRUE), 0))),
                  check.names = FALSE)
  mutate(y, across(2:ncol(y), function(x) ifelse(x > 0, 1, 0)))
}


#' Get presence or identity from BLAST
#' 
#' This function summarises BLAST results and returns presence or
#' maximum Bit score for a given query.
#' @param blast_res blast results obtained with \code{\link{get_blast_res}}
#' @param type \code{"presence"} or \code{"identity"}
#' @return a data frame of gene presence/absence or the highest identity 
#' percent value for a given gene.
#' @importFrom dplyr mutate summarise group_by ungroup
#' @importFrom tidyr pivot_wider
#' @export
#' @noRd
get_presence_from_blast <- function(blast_res, type = "presence") {
  if(!(type %in% c("presence", "identity"))) stop("Type must be 'presence' or 'identity'!")
  res <- mutate(
    summarise(
      group_by(blast_res, File, Subject),
      Presence = ifelse(type == "presence", 1, max(`Bit score`))),
    Gene = Subject)
  
  pivoted_res <- pivot_wider(res[, c("File", "Gene", "Presence")],
                             names_from = Gene, values_from = Presence, values_fill = 0)
  
  # Check for genes that were not found
  ungroup(add_missing_genes(pivoted_res))
}

