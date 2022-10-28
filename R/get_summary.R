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
#' @importFrom dplyr left_join group_by summarise case_when select mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @export
get_summary_table <- function(presence_table, hide_absent = FALSE) {
  presence_table <- select(presence_table, -"faeA")
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
    Dr = ifelse(`Afa-I` == "Present" | `Afa-III` == "Present" | F1845 == "Present",
                "Absent", Dr),
    `Afa-I` = ifelse(`Afa-I` == "Partial" & `Afa-III` == "Partial" & F1845 == "Partial" & Dr == "Present",
                     "Absent", `Afa-I`),
    `Afa-III` = ifelse(`Afa-I` == "Partial" & `Afa-III` == "Partial" & F1845 == "Partial" & Dr == "Present",
                       "Absent", `Afa-III`),
    F1845 = ifelse(`Afa-I` == "Partial" & `Afa-III` == "Partial" & F1845 == "Partial" & Dr == "Present",
                   "Absent", F1845),
    Intimin = ifelse(Intimin == "Partial" & eae == 0, "Absent", Intimin),
    CS1 = ifelse(CS17 == "Present" | CS19 == "Present" | PCF071 == "Present",
                 "Absent", CS1),
    CS17 = ifelse(CS1 == "Present" | CS19 == "Present" | PCF071 == "Present",
                  "Absent", CS17),
    CS19 = ifelse(CS17 == "Present" | CS1 == "Present" | PCF071 == "Present",
                  "Absent", CS19),
    PCF071 = ifelse(CS17 == "Present" | CS19 == "Present" | CS1 == "Present",
                    "Absent", PCF071),
    `CFA/I` = ifelse(CS14 == "Present" | CS4 == "Present",
                     "Absent", `CFA/I`),
    CS14 = ifelse(`CFA/I` == "Present" | CS4 == "Present",
                  "Absent", CS14),
    CS4 = ifelse(`CFA/I` == "Present" | CS14 == "Present",
                 "Absent", CS4),
    CS20 = ifelse(CS28A == "Present" | CS28B == "Present",
                  "Absent", CS20),
    CS28A = ifelse(CS20 == "Present" | CS28B == "Present",
                   "Absent", CS28A),
    CS28B = ifelse(CS28A == "Present" | CS20 == "Present",
                   "Absent", CS28B),
    F17a = ifelse(F17b == "Present" | F17d == "Present",
                  "Absent", F17a),
    F17b = ifelse(F17a == "Present" | F17d == "Present",
                  "Absent", F17b),
    F17d = ifelse(F17b == "Present" | F17a == "Present",
                  "Absent", F17d),
    CS31A = ifelse(F41 == "Present" | `F4/K88` == "Present" | CS23 == "Present" | CS13 == "Present",
                   "Absent", CS31A),
    F41 = ifelse(CS31A == "Present" | `F4/K88` == "Present" | CS23 == "Present" | CS13 == "Present",
                 "Absent", F41),
    `F4/K88` = ifelse(CS31A == "Present" | F41 == "Present" | CS23 == "Present" | CS13 == "Present",
                      "Absent", `F4/K88`),
    CS23 = ifelse(CS31A == "Present" | `F4/K88` == "Present" | F41 == "Present" | CS13 == "Present",
                  "Absent", CS23),
    CS13 = ifelse(CS31A == "Present" | `F4/K88` == "Present" | CS23 == "Present" | F41 == "Present",
                  "Absent", CS13),
    
    `CS27A/CS27B` = case_when((`CS27A` == "Present" | `CS27B` == "Present") | (`CS27A` == "Present" & `CS27B` == "Present" ) ~ "Present",
                              (`CS27A` == "Absent" & `CS27B` == "Partial") | (`CS27A` == "Partial" & `CS27B` == "Absent") ~ "Partial",
                              TRUE ~ "Absent"
    )
  )
  
  updated_res <- select(updated_res, -c("eae", "CS27A", "CS27B"))
  
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
  plot_dat <- pivot_longer(summary_dat, 2:ncol(summary_dat), names_to = "System", values_to = "Presence") 
  plot_dat[["Presence"]] <- sapply(plot_dat[["Presence"]], function(i) case_when(i == "Absent" ~ 0,
                                                                                 i == "Partial" ~ 0.5,
                                                                                 i == "Present" ~ 1))
  if(nrow(presence_table) > 1) {
    plot_dat <- cluster_data(plot_dat, 
                             pivot_wider(plot_dat, names_from = "System", values_from = "Presence"), 
                             "System")
  }
  
  ggplot(plot_dat, aes(x = System, y = File, fill = Presence)) +
    geom_tile() +
    scale_fill_gradient(low = absence_col, high = presence_col, name = "Presence", breaks = c(0, 0.5, 1), labels = c("Absent", "Partial", "Present")) +
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
