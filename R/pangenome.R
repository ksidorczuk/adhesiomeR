#' Expand pangenome results to genomes
#' 
#' This function extends blast results run on a pangenome to individual
#' genomes that were used for pangenome creation. 
#' @param blast_res a data frame with BLAST search results obtained by using
#' \code{\link{get_blast_res}}. 
#' @param presence_absence_file A path to the gene presence absence file in
#' a csv format. 
#' @return A data frame with BLAST results expanded to individual genomes from
#' a pangenome.
#' @seealso get_blast_res
#' @importFrom dplyr bind_rows left_join select 
#' @importFrom utils read.csv
#' @export
pangenome_to_genome <- function(blast_res, presence_absence_file) {
  gene_presence_absence <- read.csv(presence_absence_file)
  genome_cols <- which(!(colnames(gene_presence_absence) %in% c("Gene", "Non.unique.Gene.name", "Annotation", "No..isolates", 
                                                                "No..sequences", "Avg.sequences.per.isolate", "Genome.Fragment", 
                                                                "Order.within.Fragment", "Accessory.Fragment", "Accessory.Order.with.Fragment",
                                                                "QC", "Min.group.size.nuc", "Max.group.size.nuc", "Avg.group.size.nuc")))
  gpa <- gene_presence_absence[, c(1, genome_cols)]
  long_df <- bind_rows(
    lapply(1:nrow(gpa), function(i) {
      sel <- which(gpa[i, 1:ncol(gpa)] != "")
      #sel_names <- colnames(gpa)[sel[2:length(sel)]]
      #multiple_genes <- sapply(gpa[i, sel[2:length(sel)]], function(i) grep(";", strsplit(i, "")[[1]]))
      #copies <- sapply(multiple_genes, function(i) )
      data.frame(File = colnames(gpa)[sel[2:length(sel)]],
                 Gene = gpa[["Gene"]][i])
    })
  )
  
  left_join(select(blast_res, -File), long_df, by = c("Query" = "Gene"), relationship = "many-to-many")
}
