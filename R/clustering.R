#' Get cluster assignment for a single profile
#' 
#' Assigns adhesin cluster number based on gene presence/absence
#' 
#' @param model model with clustering result
#' @param new_data_single a single row of adhesiomeR results, without file column
#' @importFrom stats dist
#' @noRd
get_cluster <- function(model, new_data_single) {
  which.min(
    sapply(1:length(model[["i.med"]]), function(i) dist(rbind(data.frame(t(model[["medoids"]][i,]), check.names = FALSE),
                                                              new_data_single), method = "manhattan"))
  )
}

#' Assign genomes to adhesin profile clusters
#' 
#' Assigns analysed genomes to adhesin profile clusters using three subsets of genes:
#' all adhesins, fimbrial adhesins, nonfimbrial adhesins.
#' @param new_results a data frame of gene presence/absence generated using 
#' \code{\link{get_presence_table_strict}}
#' @param model_all clustering based on all adhesins
#' @param model_fimbrial clustering based on fimbrial adhesins
#' @param model_nonfimbrial clustering based on nonfimbrial adhesins
#' @return a data frame with four columns containing file name and assignments
#' to clusters for three types of clustering (all adhesins, fimbrial, nonfimbrial).
#' @export
get_adhesin_clusters <- function(new_results, model_all, model_fimbrial, model_nonfimbrial) {
  do.call(rbind, 
          lapply(1:nrow(new_results), function(i) {
            new_x <- new_results[which(colnames(new_results) != "File")][i,]
            data.frame(File = new_results[i,][["File"]],
                       Adhesins_cluster = get_cluster(model_all, new_x[which(colnames(new_x) %in% colnames(model_all[["data"]]))]),
                       Fimbrial_cluster = get_cluster(model_fimbrial, new_x[which(colnames(new_x) %in% colnames(model_fimbrial[["data"]]))]),
                       Nonfimbrial_cluster = get_cluster(model_nonfimbrial, new_x[which(colnames(new_x) %in% colnames(model_nonfimbrial[["data"]]))]))
          })) 
}

#' Assign genomes to adhesin profiles
#' 
#' Assigns analysed genomes to adhesin profiles using three subsets of genes:
#' all adhesins, fimbrial adhesins, nonfimbrial adhesins.
#' @param new_results a data frame of gene presence/absence generated using 
#' \code{\link{get_presence_table_strict}}
#' @param profiles a list of data frames with profiles for each subset of genes
#' @return a data frame with four columns containing file name and assignments
#' to profiles for three types of clustering (all adhesins, fimbrial, nonfimbrial).
#' @export
get_adhesin_profiles <- function(new_results, profiles = adhesiomeR::profiles) {
  data.frame(File = new_results[["File"]],
             A_profile = suppressMessages(left_join(new_results, profiles[["A"]])[["A_profile"]]),
             F_profile = suppressMessages(left_join(new_results, profiles[["F"]])[["F_profile"]]),
             N_profile = suppressMessages(left_join(new_results, profiles[["N"]])[["N_profile"]]))
}

