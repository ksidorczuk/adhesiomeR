options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")


get_blast_res <- function(input_file) {
  validate_input_file(input_file)
  input_seqs <- read_fasta(input_file)
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  blast_command <- paste0("blastn -db ../../db/adhesins -query ", input, " -out ", output, " -outfmt 6")
  write_fasta(input_seqs, input)
  system(blast_command)
  
  res <- read.delim(output, header = FALSE)
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(input, output)
  res
}


run_blast <- function(input, updateProgress = NULL) {
  
  res_df <- data.frame()
  for(i in 1:length(input[, 1])) {
    if (is.function(updateProgress)) {
      text <- paste0("Processing file: ", input[[i, 1]])
      updateProgress(detail = text)
    }
    temp <- get_blast_res(input[[i, 4]]) %>% 
      mutate(File = input[[i, 1]])
    updateProgress(value = 1)
    res_df <- rbind(res_df, temp)
  } 
  res_df
}


plot_theme <- function() {
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom")
}