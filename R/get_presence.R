#' @export
get_presence_table <- function(blast_res, add_missing = TRUE, identity_threshold = 70, evalue_treshold = 1e-50) {
  res <- blast_res %>% 
    group_by(File, Subject) %>% 
    summarise(Presence = ifelse(any(`% identity` > identity_threshold & Evalue < evalue_treshold), 1, 0)) %>% 
    mutate(Gene = sapply(Subject, function(x) strsplit(x, "~")[[1]][2]),
           System = sapply(Subject, function(x) strsplit(x, "~")[[1]][4]))
  pivoted_res <- res[, c("File", "Gene", "Presence")] %>% 
    pivot_wider(names_from = Gene, values_from = Presence, values_fill = 0)
  # Check for genes that were not found
  if(add_missing == TRUE) {
    if(length(unique(res[["Gene"]]) != length(unique(adhesins_df[["Gene"]])))) {
      pivoted_res <- add_missing_genes(pivoted_res)
    }
  }
  pivoted_res
}


#' @export
get_data_for_plots <- function(presence_table, systems = unique(adhesins_df[["System"]])) {
  all_genes <- unique(adhesins_df[["Gene"]])
  # if(show_missing == TRUE) {
  #   missing <- all_genes[which(!(all_genes %in% colnames(presence_table)))]
  #   presence_table <- cbind(presence_table, setNames(lapply(missing, function(x) x = 0), missing))
  # }
  selected_genes <- unique(filter(adhesins_df, System %in% systems)[["Gene"]])
  
  plot_dat <- presence_table %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    filter(Gene %in% selected_genes)
  
  if(nrow(presence_table) > 1) {
    dendro_files <- as.dendrogram(hclust(d = dist(x = as.matrix(presence_table[, 2:ncol(presence_table)]))))
    files_order <- order.dendrogram(dendro_files)
    dendro_genes <- as.dendrogram(hclust(d = dist(t(as.matrix(presence_table[, 2:ncol(presence_table)])))))
    genes_order <- order.dendrogram(dendro_genes)
    
    plot_dat[["Gene"]] <- factor(plot_dat[["Gene"]], 
                                 levels = colnames(presence_table[2:ncol(presence_table)][genes_order]),
                                 ordered = TRUE)
    plot_dat[["File"]] <- factor(plot_dat[["File"]],
                                 levels = presence_table[["File"]][files_order],
                                 ordered = TRUE) 
  }
  
  plot_dat[["Presence"]] <- as.factor(plot_dat[["Presence"]])
  plot_dat
}

#' @export
get_presence_plot <- function(plot_dat) {
  ggplot(plot_dat, aes(x = File, y = Gene, fill = Presence)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "bottom")
}

