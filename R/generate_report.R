generate_report_files <- function(presence_table, elements = c("summary_table", "count_table", "summary_plot", 
                                                               "presence_table", "presence_plot"), 
                                  outdir = ".", hide_absent_genes = FALSE, hide_absent_systems = FALSE, 
                                  presence_col = "#e42b24", absence_col = "#85c1ff") {

  if("summary_table" %in% elements) {
    summary_table <- get_summary_table(presence_table, hide_absent = hide_absent_systems)
    write.csv(summary_table, paste0(outdir, "/summary_table.csv"), row.names = FALSE)
  }
  
  if("count_table" %in% elements) {
    count_table <- adhesiomeR:::get_count_table(presence_table)
    if(hide_absent_systems == TRUE) {
      count_table <- filter(count_table, gene_count > 0)
    }
    write.csv(count_table, paste0(outdir, "/count_table.csv"), row.names = FALSE)
  }
  
  if("presence_table" %in% elements) {
    if(hide_absent_genes == TRUE) {
      pres_table <- presence_table[, c(TRUE, col_sums(presence_table[, 2:ncol(presence_table)]) > 0)]
    } else {
      pres_table <- presence_table
    }
    write.csv(pres_table, paste0(outdir, "/presence_table.csv"), row.names = FALSE)
  }
  
  if("summary_plot" %in% elements) {
    summary_plot <- get_summary_plot(presence_table, hide_absent = hide_absent_systems,
                                     presence_col = presence_col, absence_col = absence_col)
    t <- get_summary_table(presence_table, hide_absent = hide_absent_systems)
    ggsave(paste0(outdir, "/summary_plot.png"), summary_plot, 
           width = 80+5*ncol(t), height = 60 + 5*nrow(t), units = "mm", limitsize = FALSE)
  }
  
  if("presence_plot" %in% elements) {
    if(hide_absent_genes == TRUE) {
      plot_dat <- presence_table[, c(TRUE, col_sums(presence_table[, 2:ncol(presence_table)]) > 0)]
    } else {
      plot_dat <- presence_table
    }
    presence_plot <- get_presence_plot(plot_dat, presence_col = presence_col, absence_col = absence_col) +
      theme(legend.position = "right")
    ggsave(paste0(outdir, "/presence_plot.png"), presence_plot,
           width = 50+5*nrow(plot_dat), height = 80 + 3*ncol(plot_dat), units = "mm", limitsize = FALSE)
  }
}


generate_report <- function(presence_table, elements = c("summary_table", "count_table", "summary_plot", 
                                                         "presence_table", "presence_plot"),
                            outdir = NULL, keep_intermediate_files = TRUE, hide_absent_genes = FALSE,
                            hide_absent_systems = FALSE, presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  # Create output directory for a report
  if(is.null(outdir)) {
    outdir <- paste0(getwd(), "/adhesiomeR-report-", Sys.Date())
    if(dir.exists(outdir)) {
      outdir <- paste0(outdir, "-", format(Sys.time(), format = "%H%M%S"))
    }
  }
  dir.create(outdir)

  # List of input files
  genome_files <- presence_table[["File"]]
  
  # Generate the files specified by the 'elements' argument
  generate_report_files(presence_table, elements = c("summary_table", "count_table", "summary_plot", 
                                                     "presence_table", "presence_plot"), 
                        outdir, hide_absent_genes, hide_absent_systems,
                        presence_col, absence_col)
  
  # Generate the report file
  report_template <- system.file("adhesiomeR/adhesiomeR-report.Rmd", package = "adhesiomeR")
  rmarkdown::render(report_template, output_format = "html_document", 
                    output_dir = paste0(outdir), quiet = FALSE, "adhesiomeR-results.html",
                    params = list(genome_files, outdir, elements))
  
}
