#' @export
validate_input_file <- function(input_file) {
  if(!file.exists(input_file)) stop(paste0("File ", input_file, "does not exist.")) 
  x <- readLines(input_file)
  if(length(x) < 1) stop("Input file is empty.")
  if(length(x) < 2) stop("Input file has to be in a FASTA format.")
  if(sum(grepl("^>", x)) == 0) stop("Input file has to be in a FASTA format.")
}

#' @export
add_missing_genes <- function(results) {
  missing <- unique(adhesins_df[["Gene"]])[which(!(unique(adhesins_df[["Gene"]]) %in% colnames(results)))]
  if(length(missing) > 0) {
    cbind(results, setNames(lapply(missing, function(x) x = 0), missing))
  } else {
    results
  }
}
