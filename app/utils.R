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
  blast_command <- paste0("blastn -db ../db/adhesins -query ", input, " -out ", output, " -outfmt 6")
  write_fasta(input_seqs, input)
  system(blast_command)
  res <- read.delim(output, header = FALSE)
  colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                     "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
  file.remove(input, output)
  res
}
