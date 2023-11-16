#' Expand gene catalogue results to samples
#' 
#' This function extends blast results run on a gene catalogue to individual
#' samples that were used for gene catalogue creation. 
#' @param blast_res a data frame with BLAST search results obtained by using
#' \code{\link{get_blast_res}}. 
#' @param abundance_matrix A path to the abundance matrix file in
#' a tab-delimited format. The first column is expected to be named 'Gene' 
#' and contain gene names
#' @param threshold abundance threshold to consider a gene as present in a sample.
#' The default value is \code{NULL} and means that samples with any abundance
#' larger than 0 will be used. 
#' @return A data frame with BLAST results expanded to individual samples
#' @seealso get_blast_res
#' @importFrom dplyr bind_rows left_join select 
#' @importFrom utils read.delim
#' @export
gc_to_sample <- function(blast_res, abundance_matrix, threshold = NULL) {
  ab_mat <- read.delim(abundance_matrix)
  thresh <- if(is.null(threshold)) {
    0 
  } else {
    threshold
  }
  long_df <- bind_rows(
    lapply(1:nrow(ab_mat), function(i) {
      sel <- which(ab_mat[i, 2:ncol(ab_mat)] > thresh)
      if(length(sel) != 0) {
        data.frame(File = colnames(ab_mat)[sel+1],
                   Gene = ab_mat[["Gene"]][i])
      } else {
        data.frame(File = c(),
                   Gene = c())
      }
    })
  )
  blast_res <- filter(blast_res, Query %in% long_df[["Gene"]])
  left_join(select(blast_res, -File), long_df, by = c("Query" = "Gene"), relationship = "many-to-many")
}
