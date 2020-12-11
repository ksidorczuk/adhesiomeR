#' @export
get_summary_table <- function(presence_table) {
  presence_table %>% 
    add_missing_genes() %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    left_join(adhesins_df) %>% 
    group_by(File, System) %>% 
    summarise(gene_percentage = round(sum(Presence)*100/n(), 2)) %>% 
    filter(gene_percentage > 0) %>% 
    pivot_wider(names_from = System, values_from = gene_percentage, values_fill = 0)
}


#' @export
get_summary_plot <- function(presence_table, presence_col = "#e42b24", absence_col = "#85c1ff") {
  presence_table %>% 
    get_summary_table() %>% 
    pivot_longer(., 2:ncol(.), names_to = "System", values_to = "Percentage of present genes") %>% 
    ggplot(aes(x = System, y = File, fill = as.numeric(`Percentage of present genes`))) +
    geom_tile() +
    scale_fill_gradient(low = absence_col, high = presence_col, name = "Percentage of present genes") +
    scale_x_discrete(position = "top") +
    plot_theme()

}

#' @export
get_count_table <- function(presence_table) {
  presence_table %>% 
    add_missing_genes() %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    left_join(adhesins_df) %>% 
    group_by(System) %>% 
    summarise(gene_count = sum(Presence))
}

#' @export
get_word_cloud <- function(count_table) {
  count_table[["System"]] <- gsub(" Adhesin| Adhesins| Fimbriae", "", count_table[["System"]])
  wordcloud(count_table[["System"]], count_table[["gene_count"]], 
            colors = c("blue", "green", "red", "orange", "purple"),
            scale = c(3, 0.5), rot.per = 0.3)
}
