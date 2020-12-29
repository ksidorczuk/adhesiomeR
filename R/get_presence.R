#' @export
get_presence_table <- function(blast_res, add_missing = TRUE, identity_threshold = 70, evalue_treshold = 1e-50) {
  res <- blast_res %>% 
    group_by(File, Subject) %>% 
    summarise(Presence = ifelse(any(`% identity` > identity_threshold & Evalue < evalue_treshold), 1, 0)) %>% 
    mutate(Gene = sapply(Subject, function(x) strsplit(x, "~")[[1]][2]),
           System = sapply(Subject, function(x) strsplit(x, "~")[[1]][4])) %>% 
    filter(Presence == 1)
  pivoted_res <- res[, c("File", "Gene", "Presence")] %>% 
    pivot_wider(names_from = Gene, values_from = Presence, values_fill = 0)
  # Check for genes that were not found
  if(add_missing == TRUE) {
    if(length(unique(res[["Gene"]]) != length(unique(adhesins_df[["Gene"]])))) {
      pivoted_res <- add_missing_genes(pivoted_res)
    }
  }
  ungroup(pivoted_res)
}


#' @export
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

#' @export
get_presence_plot_data <- function(presence_table, systems = unique(adhesins_df[["System"]])) {
  
  all_genes <- unique(adhesins_df[["Gene"]])
  selected_genes <- unique(filter(adhesins_df, System %in% systems)[["Gene"]])
  
  plot_dat <- presence_table %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    filter(Gene %in% selected_genes)
  
  if(nrow(presence_table) > 1) {
    plot_dat <- cluster_data(plot_dat, presence_table, "Gene")
  }
  
  plot_dat[["Presence"]] <- factor(ifelse(plot_dat[["Presence"]] == 1, "yes", "no"), levels = c("yes", "no"))
  ungroup(plot_dat)
}


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

