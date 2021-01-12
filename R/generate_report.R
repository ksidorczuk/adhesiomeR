#' Generate a HTML report with the results
#' 
#' @export
generate_report <- function(presence_table, elements = c("summary_table", "summary_plot", 
                                                         "presence_table", "presence_plot"),
                            outdir = NULL, remove_intermediate_files = FALSE, hide_absent_genes = FALSE,
                            hide_absent_systems = FALSE, presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  # Create output directory for a report
  if(is.null(outdir)) {
    outdir <- paste0(getwd(), "/adhesiomeR-report-", Sys.Date())
  }
  if(dir.exists(outdir)) {
    outdir <- paste0(outdir, "-", format(Sys.time(), format = "%H%M%S"))
  }

  dir.create(outdir)
  outdir <- normalizePath(outdir)

  # List of input files
  genome_files <- presence_table[["File"]]
  
  # Generate the files specified by the 'elements' argument
  generate_report_files(presence_table, elements = c("summary_table", "summary_plot", 
                                                     "presence_table", "presence_plot"), 
                        outdir = outdir, hide_absent_genes = hide_absent_genes, 
                        hide_absent_systems = hide_absent_systems,
                        presence_col = presence_col, absence_col = absence_col)
  
  # Generate the report file
  report_template <- system.file("adhesiomeR/adhesiomeR-report.Rmd", package = "adhesiomeR")
  rmarkdown::render(report_template, output_format = "html_document", 
                    output_dir = paste0(outdir), quiet = TRUE, "adhesiomeR-results.html",
                    params = list(genome_files, outdir, elements))
  
  if(remove_intermediate_files == TRUE) {
    fl <- list.files(outdir, full.names = TRUE)
    invisible(file.remove(fl[!grepl("adhesiomeR-results.html", fl)]))
  }
  
}
