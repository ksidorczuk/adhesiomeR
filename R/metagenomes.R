#' Expand gene catalogue results to samples
#' 
#' This function extends blast results run on a gene catalogue to individual
#' samples that were used for gene catalogue creation. 
#' @param blast_res a data frame with BLAST search results obtained by using
#' \code{\link{get_blast_res}}. 
#' @param presence_absence_file A path to the abundance matrix file in
#' a tab-delimited format. The first column is expected to contain gene names.
#' @return A data frame with BLAST results expanded to individual samples
#' @seealso get_blast_res
#' @importFrom dplyr bind_rows left_join select 
#' @importFrom utils read.delim
#' @export
gc_to_sample <- function(blast_res, abundance_matrix, threshold = NULL) {
  ab_mat <- read.delim(abundance_matrix)
  first_col <- colnames(ab_mat)[1]
  tresh <- if(is.null(threshold)) {
    0 
  } else {
    threshold
  }
  long_df <- bind_rows(
    lapply(1:nrow(ab_mat), function(i) {
      print(i)
      sel <- which(ab_mat[i, 2:ncol(ab_mat)] > threshold)
      if(length(sel) != 0) {
        data.frame(File = colnames(ab_mat)[sel+1],
                   Gene = ab_mat[[first_col]][i])
      } else {
        data.frame(File = c(),
                   Gene = c())
      }
    })
  )
  
  left_join(select(blast_res, -File), long_df, by = c("Query" = first_col), relationship = "many-to-many")
}
