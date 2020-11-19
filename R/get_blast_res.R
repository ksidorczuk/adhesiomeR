get_blast_res <- function(input_file_list) {
  lapply(input_file_list, function(ith_file) {
    do_blast_single(ith_file)
  }) %>% bind_rows()
}

do_blast_single <- function(input_file) {
  system(paste0("blastn -db db/adhesins -query ", input_file, " -out ", input_file, ".blast -outfmt 6"))
  name <- last(strsplit(input_file, "/")[[1]])
  res <- read.delim(paste0(input_file, ".blast"), header = FALSE)
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(paste0(input_file, ".blast"))
  mutate(res, File = name)
}
