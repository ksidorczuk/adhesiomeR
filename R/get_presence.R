get_presence_table <- function(blast_res, add_missing = TRUE) {
  res <- blast_res %>% 
    group_by(File, Subject) %>% 
    summarise(Presence = ifelse(any(`% identity` > 70 & Evalue < 1e-50), 1, 0)) %>% 
    mutate(Gene = sapply(Subject, function(x) strsplit(x, "~")[[1]][2]),
           System = sapply(Subject, function(x) strsplit(x, "~")[[1]][4]))
  pivoted_res <- res[, c("File", "Gene", "Presence")] %>% 
    pivot_wider(names_from = Gene, values_from = Presence, values_fill = 0)
  # Check for genes that were not found
  if(add_missing == TRUE) {
    if(length(unique(res[["Gene"]]) != length(adhesins_df[["Gene"]]))) {
      missing <- adhesins_df[["Gene"]][which(!(adhesins_df[["Gene"]] %in% unique(res[["Gene"]])))]
      pivoted_res <- cbind(pivoted_res, setNames(lapply(missing, function(x) x = 0), missing))
    }
  }
  pivoted_res
}



get_data_for_plots <- function(presence_table, systems = unique(adhesins_df[["System"]])) {
  all_genes <- adhesins_df[["Gene"]]
  # if(show_missing == TRUE) {
  #   missing <- all_genes[which(!(all_genes %in% colnames(presence_table)))]
  #   presence_table <- cbind(presence_table, setNames(lapply(missing, function(x) x = 0), missing))
  # }
  selected_genes <- filter(adhesins_df, System %in% systems)[["Gene"]]
  
  plot_dat <- presence_table %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    filter(Gene %in% selected_genes)
  
  plot_dat[["Gene"]] <- factor(plot_dat[["Gene"]], levels = all_genes)
  plot_dat[["Presence"]] <- as.factor(plot_dat[["Presence"]]) 
  
  plot_dat
}


get_presence_plot <- function(plot_dat) {
  ggplot(plot_dat, aes(x = File, y = Gene, fill = Presence)) +
    geom_tile() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "bottom")
}

