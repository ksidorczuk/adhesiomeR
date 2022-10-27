#' Get BLAST results for many files
#' 
#' This function runs BLAST search on a provided list of files.
#' @param input_file_list A \code{list} of file names that will be used
#' as input files for BLAST search. File have to contain nucleotide sequences
#' in a FASTA format. 
#' @param nt number of threads used for running BLAST. Default is one. The 
#' maximum number of threads is determined by \code{\link[parallel]{detectCores}}.
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
#' @importFrom parallel makePSOCKcluster stopCluster detectCores
#' @importFrom future.apply future_lapply
#' @importFrom future makeClusterPSOCK plan
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom dplyr mutate
#' @importFrom progressr with_progress progressor handlers handler_progress
#' @export
get_blast_res <- function(input_file_list, nt = 1, blast_dir = Sys.which("blastn")) {
  max_nt <- detectCores(logical = FALSE)
  if(nt > max_nt) {
    stop(paste0("The number of threads you specified is too large. The maximum number of threads determined by parallel::detectCores function is: ", detectCores(logical = FALSE), ". 
  Please select a value between 1 and ", detectCores(logical = FALSE), "."))
  } else if (!(nt %in% 1L:detectCores(logical = FALSE)))  {
    stop("The number of threads is incorrect. Please make sure that you entered a valid number.")
  }
  
  if(nt > 1) {
    plan(multisession, workers = nt, gc = TRUE)
    
    handlers(handler_progress(format="[:bar] :percent :eta :message"))
    
    with_progress({ 
      p <- progressor(along = 1:length(input_file_list))
      res <- bind_rows(
        future_lapply(1:length(input_file_list), function(i) {
          p()
          do_blast_single(input_file_list[[i]], blast_dir)
        })
      )
    }) 
    plan(sequential)
  } else {
    with_progress({ 
      p <- progressor(along = 1:length(input_file_list))
      res <- bind_rows(
        lapply(1:length(input_file_list), function(i) {
          p()
          do_blast_single(input_file_list[[i]], blast_dir)
        })
      )
    }) 
  }
  
  mutate(res, Subject = sapply(Subject, function(i) strsplit(i, "~~~")[[1]][2]))
  
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
  input_file_name <- last(strsplit(input_file, "/")[[1]])
  if(blast_dir == "") {
    stop("It seems that you do not have BLAST installed. To be able to use adhesiomeR,
    you should install standalone BLAST. If you have installed BLAST and still
    see this message, please use 'blast_dir' argument to provide a proper path to 
    the BLAST directory in which 'blastn' executable is located.")
  } else if(blast_dir == Sys.which("blastn")) {
    system(paste0(blast_dir, " -db ", db_path, " -query ", gsub(" ", "\\ ", input_file, fixed = TRUE), " -out ", gsub(" ", "\\ ", input_file_name, fixed = TRUE), ".blast -outfmt 6"))
  } else {
    if(grepl("/$", blast_dir)) blast_dir <- gsub("/$", "", blast_dir)
    tryCatch(system(paste0(blast_dir, "/blastn -db ", db_path, " -query ", gsub(" ", "\\ ", input_file, fixed = TRUE), " -out ", gsub(" ", "\\ ", input_file_name, fixed = TRUE), ".blast -outfmt 6")),
             warning = function(w) {
               msg <- conditionMessage(w)
               if(msg == "error in running command") message("Please check if the BLAST directory path you provided is correct.")
             })
  }
  res <- tryCatch(
    read.delim(paste0(input_file_name, ".blast"), header = FALSE), 
    error = function(e) {
      msg <- conditionMessage(e)
      if(msg == "no lines available in input") as.data.frame(matrix(nrow = 1, ncol = 12))
    })
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(paste0(input_file_name, ".blast"))
  mutate(res, File = input_file_name)
}

