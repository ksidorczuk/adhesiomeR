#' Get table with summarized system presence
#' 
#' This functions creates a table with an overview of systems presence
#' in analysed genomes. Each system is summarized as a percent of genes
#' that were identified. 
#' @param presence_table a data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function
#' @param hide_absent \code{logical} indicating if columns representing 
#' systems that were not found in any file should be displayed. By default
#' \code{FALSE}
#' @return a data frame with systems presence indicated by a percentage of
#' found genes. First column contains names of the input files and the following
#' correspond to analysed systems. 
#' @importFrom dplyr left_join group_by summarise filter n case_when select mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
get_summary_table <- function(presence_table, hide_absent = FALSE) {
  presence_table <- select(presence_table, -c("faeA", "draP", "daaP"))
  res <- mutate(
    summarise(
      group_by(
        left_join(
          pivot_longer(add_missing_genes(presence_table, type = "grouped"), 
                       2:ncol(add_missing_genes(presence_table, type = "grouped")), 
                       names_to = "Gene", values_to = "Presence"), 
          adhesins_df_grouped, by = "Gene"),
        File, System),
      gene_percentage = round(sum(Presence)*100/n(), 2)),
    gene_percentage = case_when(gene_percentage == 100 ~ "Present",
                                gene_percentage > 0 & gene_percentage < 100 ~ "Partial",
                                gene_percentage == 0 ~ "Absent")
  )
  
  pivoted_res <- left_join(pivot_wider(res, File, names_from = "System", values_from = "gene_percentage", values_fill = "Absent"),
                           presence_table[, c("File", "eae")])
  
  updated_res <- mutate(
    pivoted_res,
    Dr_fimbriae = ifelse(`Afa-I_adhesins` == "Present" | `Afa-III_adhesins` == "Present" | F1845_fimbriae == "Present",
                         "Absent", Dr_fimbriae),
    `Afa-I_adhesins` = ifelse(`Afa-I_adhesins` == "Partial" & `Afa-III_adhesins` == "Partial" & F1845_fimbriae == "Partial" & Dr_fimbriae == "Present",
                              "Absent", `Afa-I_adhesins`),
    `Afa-III_adhesins` = ifelse(`Afa-I_adhesins` == "Partial" & `Afa-III_adhesins` == "Partial" & F1845_fimbriae == "Partial" & Dr_fimbriae == "Present",
                                "Absent", `Afa-III_adhesins`),
    F1845_fimbriae = ifelse(`Afa-I_adhesins` == "Partial" & `Afa-III_adhesins` == "Partial" & F1845_fimbriae == "Partial" & Dr_fimbriae == "Present",
                            "Absent", F1845_fimbriae),
    Intimin = ifelse(Intimin == "Partial" & eae == 0, "Absent", Intimin),
    CS1_fimbriae = ifelse(CS17_fimbriae == "Present" | CS19_fimbriae == "Present" | PCF071_fimbriae == "Present",
                          "Absent", CS1_fimbriae),
    CS17_fimbriae = ifelse(CS1_fimbriae == "Present" | CS19_fimbriae == "Present" | PCF071_fimbriae == "Present",
                           "Absent", CS17_fimbriae),
    CS19_fimbriae = ifelse(CS17_fimbriae == "Present" | CS1_fimbriae == "Present" | PCF071_fimbriae == "Present",
                           "Absent", CS19_fimbriae),
    PCF071_fimbriae = ifelse(CS17_fimbriae == "Present" | CS19_fimbriae == "Present" | CS1_fimbriae == "Present",
                             "Absent", PCF071_fimbriae),
    `CFA/I_fimbriae` = ifelse(CS14_fimbriae == "Present" | CS4_fimbriae == "Present",
                              "Absent", `CFA/I_fimbriae`),
    CS14_fimbriae = ifelse(`CFA/I_fimbriae` == "Present" | CS4_fimbriae == "Present",
                           "Absent", CS14_fimbriae),
    CS4_fimbriae = ifelse(`CFA/I_fimbriae` == "Present" | CS14_fimbriae == "Present",
                          "Absent", CS4_fimbriae),
    CS20_fimbriae = ifelse(CS28A_Fimbriae == "Present" | CS28B_Fimbriae == "Present",
                           "Absent", CS20_fimbriae),
    CS28A_Fimbriae = ifelse(CS20_fimbriae == "Present" | CS28B_Fimbriae == "Present",
                            "Absent", CS28A_Fimbriae),
    CS28B_Fimbriae = ifelse(CS28A_Fimbriae == "Present" | CS20_fimbriae == "Present",
                            "Absent", CS28B_Fimbriae),
    F17a_fimbriae = ifelse(F17b_fimbriae == "Present" | F17d_fimbriae == "Present",
                           "Absent", F17a_fimbriae),
    F17b_fimbriae = ifelse(F17a_fimbriae == "Present" | F17d_fimbriae == "Present",
                           "Absent", F17b_fimbriae),
    F17d_fimbriae = ifelse(F17b_fimbriae == "Present" | F17a_fimbriae == "Present",
                           "Absent", F17d_fimbriae),
    CS31A_Fimbriae = ifelse(F41_Fimbriae == "Present" | K88_Fimbriae == "Present" | CS23_adhesins == "Present" | CS13_fimbriae == "Present",
                            "Absent", CS31A_Fimbriae),
    F41_Fimbriae = ifelse(CS31A_Fimbriae == "Present" | K88_Fimbriae == "Present" | CS23_adhesins == "Present" | CS13_fimbriae == "Present",
                          "Absent", F41_Fimbriae),
    K88_Fimbriae = ifelse(CS31A_Fimbriae == "Present" | F41_Fimbriae == "Present" | CS23_adhesins == "Present" | CS13_fimbriae == "Present",
                          "Absent", K88_Fimbriae),
    CS23_adhesins = ifelse(CS31A_Fimbriae == "Present" | K88_Fimbriae == "Present" | F41_Fimbriae == "Present" | CS13_fimbriae == "Present",
                           "Absent", CS23_adhesins),
    CS13_fimbriae = ifelse(CS31A_Fimbriae == "Present" | K88_Fimbriae == "Present" | CS23_adhesins == "Present" | F41_Fimbriae == "Present",
                           "Absent", CS13_fimbriae),
    
    `CS27A/CS27B_fimbriae` = case_when((`CS27A_Fimbriae` == "Present" | `CS27B_Fimbriae` == "Present") | (`CS27A_Fimbriae` == "Present" & `CS27B_Fimbriae` == "Present" ) ~ "Present",
                                       (`CS27A_Fimbriae` == "Absent" & `CS27B_Fimbriae` == "Partial") | (`CS27A_Fimbriae` == "Partial" & `CS27B_Fimbriae` == "Absent") ~ "Partial",
                                       TRUE ~ "Absent"
    )
  )
  
  updated_res <- select(updated_res, -c("eae", "CS27A_Fimbriae", "CS27B_Fimbriae"))
  
  if(hide_absent == TRUE) {
    x <- colSums(updated_res == "Absent")
    updated_res <- updated_res[, which(x != nrow(updated_res))]
  }
  updated_res
}


#' Get plot with summarized system presence
#' 
#' This function generates a heatmap showing the percent of genes found
#' for each system in each analysed file. 
#' @param presence_table a data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function
#' @param hide_absent \code{logical} indicating if columns representing 
#' systems that were not found in any file should be displayed. By default
#' \code{FALSE}
#' @param presence_col color of the tiles representing present genes. Must be
#' specified as a hex color code
#' @param absence_col color of the tiles representing absent genes. Must be
#' specified as a hex color code
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_tile scale_fill_gradient scale_x_discrete aes
#' @export
get_summary_plot <- function(presence_table, hide_absent = FALSE, presence_col = "#e42b24", absence_col = "#85c1ff") {
  summary_dat <- get_summary_table(presence_table, hide_absent) 
  plot_dat <- pivot_longer(summary_dat, 2:ncol(summary_dat), names_to = "System", values_to = "Percentage of present genes") 
  
  if(nrow(presence_table) > 1) {
    plot_dat <- cluster_data(plot_dat, summary_dat, "System")
  }
  
  ggplot(plot_dat, aes(x = System, y = File, fill = as.numeric(`Percentage of present genes`))) +
    geom_tile() +
    scale_fill_gradient(low = absence_col, high = presence_col, name = "Percentage of present genes") +
    scale_x_discrete(position = "top") +
    plot_theme()
  
}

#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join group_by summarise
get_count_table <- function(presence_table) {
  summarise(
    group_by(
      left_join(
        pivot_longer(add_missing_genes(presence_table), 2:ncol(add_missing_genes(presence_table)), names_to = "Gene", values_to = "Presence"),
        adhesins_df), 
      System),
    gene_count = sum(Presence))
}


# get_word_cloud <- function(count_table) {
#   count_table[["System"]] <- gsub(" Adhesin| Adhesins| Fimbriae", "", count_table[["System"]])
#   wordcloud(count_table[["System"]], count_table[["gene_count"]], 
#             colors = c("blue", "green", "red", "orange", "purple"),
#             scale = c(3, 0.5), rot.per = 0.3)
# }
