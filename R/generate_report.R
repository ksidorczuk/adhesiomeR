#' Generate a HTML report with the results
#' 
#' This function generates a HTML report with the results of a performed
#' analysis, using an object obtained by \code{\link{get_presence_table}}
#' function.
#' 
#' @param presence_table an object obtained by \code{\link{get_presence_table}}
#' function.
#' @param elements a character vector of elements which should be included
#' in the report. The available values are: \code{summary_table}, 
#' \code{summary_plot}, \code{presence_table} and \code{presence_plot}.
#' By default, all elements are used.
#' @param outdir a string indicating the path and name of the directory,
#' in which the report will be generated.
#' @param remove_intermediate_files a logical indicating if the files
#' with tables and/or plots used for report generation should be removed
#' from the results directory. By default \code{FALSE}.
#' @param hide_absent_genes a logical indicating if the genes that were
#' not found in any of the files should be hidden from the results.
#' @param hide_absent_systems a logical indicating if the systems that were
#' not found in any of the files should be hidden from the results.
#' @param presence_col color of the tiles representing present genes. Must be
#' specified as a hex color code
#' @param absence_col color of the tiles representing absent genes. Must be
#' specified as a hex color code
#' @return No return value, called for its side effects.
#' @details This function uses object obtained by \code{get_presence_table}
#' function to create files with the results, i.e., plots and/or tables, 
#' depending on the selected elements. These files are placed in a newly 
#' created directory together with the report file, which location and name
#' may be specified using \code{outdir} argument. To get only the report file
#' in a results directory, \code{remove_intermediate_files} argument may be
#' used. 
#' 
#' The \code{elements} available to include in a report are:
#' \itemize{
#'   \item{summary_table}{ table with percents of genes found in each system}
#'   \item{summary_plot}{ plot with percents of genes found in each system}
#'   \item{presence_table}{ table with gene presence/absence}
#'   \item{presence_plot}{ plot with gene presence/absence.}}
#' @importFrom rmarkdown render
#' @importFrom utils write.csv
#' @export
generate_report <- function(presence_table, elements = c("summary_table", "summary_plot", 
                                                         "presence_table", "presence_plot"),
                            outdir = NULL, remove_intermediate_files = FALSE, hide_absent_genes = FALSE,
                            hide_absent_systems = FALSE, presence_col = "#e42b24", absence_col = "#85c1ff") {
  
  if(is.null(elements) | !all(elements %in% c("summary_table", "summary_plot", 
                                              "presence_table", "presence_plot"))) {
    stop("Incorrect elements selected. The 'elements' argument may contain only the following items: 'summary_table', 'summary_plot', 'presence_table', 'presence_plot'. Please make sure that you use at least one of these items.")
  }
  
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
  render(report_template, output_format = "html_document", 
         output_dir = paste0(outdir), quiet = TRUE, "adhesiomeR-results.html",
         params = list(genome_files, outdir, elements))
  
  if(remove_intermediate_files == TRUE) {
    fl <- list.files(outdir, full.names = TRUE)
    invisible(file.remove(fl[!grepl("adhesiomeR-results.html", fl)]))
  }
}
