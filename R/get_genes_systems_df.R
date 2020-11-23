#' Get gene-system data frame
#' 
#' Extract systems and corresponding genes from a sequence
#' file used for database creation.
#' @param seq_file File with sequences used for database creation 
#' in a FASTA format. See details section for FASTA header format requirements.
#' @return a data frame containing systems and corresponding genes.
#' It consists of three columns:
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{System}{System to which given gene is assigned}}
#' Each row contains only one gene and system, so in case of genes assigned to 
#' multiple systems, they will occupy several rows.
#' @details The sequence file should be in a FASTA format with headers 
#' containing database name, gene name, accession, system name and each element
#' should be separated by \code{~} without spaces, e.g. 
#' \code{>Adhesins~fimH~NC_000913.3:4548808-4549710~Type_1_Fimbriae}. If a gene
#' is assigned to multiple systems, system names should be separated using 
#' semicolons without spaces, e.g. \code{Adhesins~aggR~NC_019000.1:48472-49269~AA/I_Fimbriae;
#' AA/II_Fimbriae;AA/III_Fimbriae}.
#' 
#' @export
#' @examples 
#' get_genes_in_systems_db("data/Adhesins_sequences")
#' 
#' @importFrom dplyr mutate
get_genes_in_systems_db <- function(seq_file) {
  data_file <- readLines(seq_file)
  def_lines <- data_file[grepl(">", data_file)]
  unique_genes <- unique(sapply(def_lines, function(x) strsplit(x, "~")[[1]][2], USE.NAMES = FALSE)) 
  df <- data.frame(Gene = sapply(def_lines, function(x) strsplit(x, "~")[[1]][2], USE.NAMES = FALSE),
                   System = sapply(def_lines, function(x) strsplit(x, "~")[[1]][4], USE.NAMES = FALSE),
                   stringsAsFactors = FALSE) %>% 
    mutate(System = gsub("^ ", "", System))
  lapply(unique_genes, function(ith_gene) {
    data.frame(Gene = ith_gene,
               System = strsplit(df[["System"]][which(df[["Gene"]] == ith_gene)], ";")[[1]],
               stringsAsFactors = FALSE) %>% 
      mutate(System = gsub("_", " ", System))
  }) %>% bind_rows()
}
