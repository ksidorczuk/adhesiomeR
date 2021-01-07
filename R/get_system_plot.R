#' Get plot with gene presence/absence for one system
#' 
#' This function creates a heatmap indicating gene presence/absence for one
#' selected system. 
#' @param presence_table a data frame with gene presence/absence obtained 
#' using \code{\link{get_presence_table}} function
#' @param system \code{character} string indicating a name of selected system
#' @param presence_col color of the tiles representing present genes. Must be
#' specified as a hex color code
#' @param absence_col color of the tiles representing absent genes. Must be
#' specified as a hex color code
#' @importFrom dplyr left_join ungroup %>% filter
#' @importFrom ggplot2 ggplot geom_tile scale_fill_manual ggtitle theme
#' @export
get_system_plot <- function(presence_table, system, presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  plot_dat <- get_presence_plot_data(presence_table) %>% 
    left_join(., adhesins_df, by = "Gene") %>% 
    ungroup()
  
    filter(plot_dat, System == system) %>% 
    ggplot(aes(x = Gene, y = File, fill = Presence)) +
      geom_tile() +
      scale_fill_manual("Presence", values = c("yes" = presence_col, "no" = absence_col), drop = FALSE) +
      ggtitle(system) +
      plot_theme() +
      theme(legend.position = "right")
}
