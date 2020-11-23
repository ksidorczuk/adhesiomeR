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


get_count_table <- function(presence_table) {
  presence_table %>% 
    add_missing_genes() %>% 
    pivot_longer(., 2:ncol(.), names_to = "Gene", values_to = "Presence") %>% 
    left_join(adhesins_df) %>% 
    group_by(System) %>% 
    summarise(gene_count = sum(Presence))
}
#counts <- get_count_table(presence_tab1)

get_word_cloud <- function(count_table) {
  count_table[["System"]] <- gsub("_adhesin|_adhesins|_fimbriae|_Fimbriae", "", count_table[["System"]] )
  count_table[["System"]] <- gsub("_", " ", count_table[["System"]] )
  wordcloud(count_table[["System"]], count_table[["gene_count"]], 
            colors = c("blue", "green", "yellow", "red", "orange", "purple"),
            scale = c(4, 0.5), rot.per = 0.3)
}
