#' @export
get_system_plot <- function(plot_dat, system, presence_col = "#e42b24", absence_col = "#85c1ff") {
    filter(plot_dat, System == system) %>% 
    ggplot(aes(x = Gene, y = File, fill = Presence)) +
      geom_tile() +
      scale_fill_manual("Presence", values = c("1" = presence_col, "0" = absence_col)) +
      ggtitle(system) +
      plot_theme()
}
