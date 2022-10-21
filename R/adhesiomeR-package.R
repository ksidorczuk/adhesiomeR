#' Analysis of *Escherichia coli* adhesiome
#' 
#' For adhesiome analysis, adhesiomeR uses standalone BLAST and a manually curated 
#' database of adhesins. Our database contains **433 genes** comprising **74 systems** 
#' of fimbriae and other adhesins; therefore adhesiomeR provides the most 
#' comprehensive collection of *Escherichia coli* adhesins currently available. 
#' 
#' AdhesiomeR is available as an R package and shiny GUI.
#' 
#' @name adhesiomeR-package
#' @aliases adhesiomeR-package adhesiomeR
#' @docType package
#' @importFrom utils menu
#' @author
#' Maintainer: Katarzyna Sidorczuk <sidorczuk.katarzyna17@@gmail.com>
#' @keywords package
NULL

globalVariables(c("adhesins_df", "i", "Subject", "System", "Presence", "File", 
                  "Gene", "% identity", "gene_percentage", "Evalue", "Percentage of present genes"))

.onAttach <- function(libname = find.package("adhesiomeR"), pkgname = "adhesiomeR") {
  options(dplyr.summarise.inform = FALSE)
}
