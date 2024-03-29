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
get_adhesin_clusters <- function(new_results, 
                                 model_all = adhesiomeR::clustering_all, 
                                 model_fimbrial = adhesiomeR::clustering_fimbrial, 
                                 model_nonfimbrial = adhesiomeR::clustering_nonfimbrial) {
  ver <- attributes(new_results)[["search_version"]]
  if(is.null(ver)) {
    stop(paste0("The input type is not supported for assignment of clustering. Please make sure to use results obtained with `get_presence_table_strict`."))
  } else if(ver == "relaxed") {
    stop("This function is not compatible with relaxed version of the search. To use cluster assignments, please generate presence table using `get_presence_table_strict` function.")
  } else if(ver == "strict") {
    do.call(rbind, 
            lapply(1:nrow(new_results), function(i) {
              new_x <- new_results[which(colnames(new_results) != "File")][i,]
              if(all(new_x == 0)) {
                data.frame(File = new_results[i,][["File"]],
                           Adhesins_cluster = NA,
                           Fimbrial_cluster = NA,
                           Nonfimbrial_cluster = NA)
              } else {
                data.frame(File = new_results[i,][["File"]],
                           Adhesins_cluster = c("A-G", "A-I", "A-E", "A-C", "A-A","A-D", "A-J", "A-B", "A-F", "A-H")
                           [get_cluster(model_all, new_x[which(colnames(new_x) %in% colnames(model_all[["data"]]))])],
                           Fimbrial_cluster = c("F-F", "F-G", "F-A", "F-D", "F-C", "F-E", "F-H", "F-B")
                           [get_cluster(model_fimbrial, new_x[which(colnames(new_x) %in% colnames(model_fimbrial[["data"]]))])],
                           Nonfimbrial_cluster = c("N-A", "N-D", "N-C", "N-B", "N-E")
                           [get_cluster(model_nonfimbrial, new_x[which(colnames(new_x) %in% colnames(model_nonfimbrial[["data"]]))])])
              }
            }))
  }
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
  ver <- attributes(new_results)[["search_version"]]
  if(is.null(ver)) {
    stop(paste0("The input type is not supported for assignment of clustering. Please make sure to use results obtained with `get_presence_table_strict`."))
  } else if(ver == "relaxed") {
    stop("This function is not compatible with relaxed version of the search. To use cluster assignments, please generate presence table using `get_presence_table_strict` function.")
  } else if(ver == "strict") {
    data.frame(File = new_results[["File"]],
               Adhesins_profile = suppressMessages(left_join(new_results, profiles[["A"]])[["Adhesins_profile"]]),
               Fimbrial_profile = suppressMessages(left_join(new_results, profiles[["F"]])[["Fimbrial_profile"]]),
               Nonfimbrial_profile = suppressMessages(left_join(new_results, profiles[["N"]])[["Nonfimbrial_profile"]]))
  }
}

