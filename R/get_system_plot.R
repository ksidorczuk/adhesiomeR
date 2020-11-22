get_system_plot <- function(plot_dat, system) {
    filter(plot_dat, System == system) %>% 
    ggplot(aes(x = Gene, y = File, fill = Presence)) +
      geom_tile() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "bottom") +
      ggtitle(system)
}