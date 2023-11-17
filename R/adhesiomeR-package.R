#' Analysis of *Escherichia coli* adhesiome
#' 
#' For adhesiome analysis, adhesiomeR uses standalone BLAST and a manually curated 
#' database of adhesins. Our database contains **525 genes** comprising **102 systems** 
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
                  "Gene", "% identity", "gene_percentage", "Evalue", "Percentage of present genes",
                  "Query start", "Query end", "adhesins_df_grouped",
                  "Subject coverage", "Bit score", "Threshold", "Afa-I",
                  "Afa-III", "F1845", "Dr", "CS17", "CS19", "PCF071",
                  "CS1", "CS4", "CS14", "CFA/I", "CS28A", "CS28B", "CS20",
                  "F17b", "F17a", "F17d", "F41", "F4/K88", "CS23", "CS13",
                  "CS31A", "same_location", "Query"))

.onAttach <- function(libname = find.package("adhesiomeR"), pkgname = "adhesiomeR") {
  options(dplyr.summarise.inform = FALSE)
}
