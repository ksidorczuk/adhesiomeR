options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")


# get_blast_res <- function(input_file, nt) {
#   validate_input_file(input_file)
#   input_seqs <- read_fasta(input_file)
#   input <- tempfile(tmpdir = getwd())
#   output <- tempfile(tmpdir = getwd())
#   db_path <- paste0(normalizePath(system.file(package = "adhesiomeR")), "/db/adhesins")
#   blast_command <- paste0("blastn -db ", db_path, " -query ", input, " -out ", output, " -outfmt 6")
#   write_fasta(input_seqs, input)
#   system(blast_command)
#   
#   res <- read.delim(output, header = FALSE)
#   colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
#                      "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
#   file.remove(input, output)
#   res
# }


run_blast <- function(input_files, nt, updateProgress = NULL) {
  parallel_cluster <- makePSOCKcluster(nt)
  registerDoParallel(parallel_cluster)
  
  res <- foreach(i = 1:length(input_files[, 1]), .combine = rbind, .packages = c('dplyr', 'adhesiomeR', 'biogram')) %dopar% {
    
    # if (is.function(updateProgress)) {
    #   text <- paste0("Processing file: ", input_files[[i, 1]])
    #   updateProgress(detail = text)
    # }
    # 
    adhesiomeR:::validate_input_file(input_files[[i, 4]])
    input_seqs <- read_fasta(input_files[[i, 4]])
    input <- tempfile(tmpdir = getwd())
    output <- tempfile(tmpdir = getwd())
    db_path <- paste0(normalizePath(system.file(package = "adhesiomeR")), "/db/adhesins")
    blast_command <- paste0("blastn -db ", db_path, " -query ", input, " -out ", output, " -outfmt 6")
    write_fasta(input_seqs, input)
    system(blast_command)
    
    res <- read.delim(output, header = FALSE)
    colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                       "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")

    file.remove(input, output)
   # updateProgress(value = 1)
    mutate(res, File = input_files[[i, 1]])
  } 
  stopCluster(parallel_cluster)
  res
}

