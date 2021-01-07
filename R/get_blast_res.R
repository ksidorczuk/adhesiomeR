#' Get BLAST results for many files
#' 
#' This function runs BLAST search on a provided list of files.
#' @param input_file_list A \code{list} of file names that will be used
#' as input files for BLAST search. File have to contain nucleotide sequences
#' in a FASTA format. 
#' @param nt number of threads used for running BLAST. 
#' @param blast_dir A path to the directory with BLAST executables. By default,
#' it tries to find the proper executable by running system command \code{which}
#' (except on Windows). If the \code{blastn} executable is not found by default,
#' you should provide a path to the directory in which it is located.
#' @return A data frame with BLAST results containing following information: 
#' query name, subject name, % of identity, alignment length, mismatches,
#' gap opens, query start, query end, subject start, subject end, E-value,
#' bit score (for detailed information about those outputs, please see BLAST
#' documentation), as well as name of a file that was used as an input.
#' @seealso do_blast_single
#' @importFrom pbapply pblapply
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach %dopar%
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
get_blast_res <- function(input_file_list, nt = 1, blast_dir = Sys.which("blastn")) {
  parallel_cluster <- makePSOCKcluster(nt)
  #registerDoParallel(parallel_cluster)
  registerDoSNOW(parallel_cluster)
  pb <- txtProgressBar(max = length(input_file_list), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  res <- foreach(i = 1:length(input_file_list), .combine = rbind, .packages = c('dplyr', 'adhesiomeR', 'biogram'), 
                 .options.snow=opts) %dopar% {
    cat(paste0("Processing ", input_file_list[[i]], ". (File ", i, "/", length(input_file_list), " files)"))
    do_blast_single(input_file_list[[i]], blast_dir)
  }
  stopCluster(parallel_cluster)
  res
    
}


#' Get BLAST results for a single file
#' 
#' This function runs BLAST search on a provided file.
#' @param input_file Name of a file that will be used as input files for BLAST search. 
#' File have to contain nucleotide sequences in a FASTA format. 
#' @param blast_dir A path to the directory with BLAST executables. By default,
#' it tries to find the proper executable by running system command \code{which}
#' (except on Windows). If the \code{blastn} executable is not found by default,
#' you should provide a path to the directory in which it is located.
#' @return A data frame with BLAST results containing following information: 
#' query name, subject name, % of identity, alignment length, mismatches,
#' gap opens, query start, query end, subject start, subject end, E-value,
#' bit score (for detailed information about those outputs, please see BLAST
#' documentation), as well as name of a file that was used as an input.
#' @seealso get_blast_res
#' @importFrom dplyr last
#' @importFrom utils read.delim
#' @export
do_blast_single <- function(input_file, blast_dir = Sys.which("blastn")) {
  validate_input_file(input_file)
  db_path <- paste0(normalizePath(system.file(package = "adhesiomeR")), "/db/adhesins")
  if(blast_dir == "") {
    stop("It seems that you do not have BLAST installed. To be able to use adhesiomeR,
    you should install standalone BLAST. If you have installed BLAST and still
    see this message, please use 'blast_dir' argument to provide a proper path to 
    the BLAST directory in which 'blastn' executable is located.")
  } else if(blast_dir == Sys.which("blastn")) {
    system(paste0(blast_dir, " -db ", db_path, " -query ", input_file, " -out ", input_file, ".blast -outfmt 6"))
  } else {
    if(grepl("/$", blast_dir)) blast_dir <- gsub("/$", "", blast_dir)
    tryCatch(system(paste0(blast_dir, "/blastn -db ", db_path, " -query ", input_file, " -out ", input_file, ".blast -outfmt 6")),
             warning = function(w) {
               msg <- conditionMessage(w)
               if(msg == "error in running command") message("Please check if the BLAST directory path you provided is correct.")
             })
  }
  name <- last(strsplit(input_file, "/")[[1]])
  res <- read.delim(paste0(input_file, ".blast"), header = FALSE)
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(paste0(input_file, ".blast"))
  mutate(res, File = name)
}

