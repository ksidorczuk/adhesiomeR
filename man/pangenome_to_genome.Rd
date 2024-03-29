% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pangenome.R
\name{pangenome_to_genome}
\alias{pangenome_to_genome}
\title{Expand pangenome results to genomes}
\usage{
pangenome_to_genome(blast_res, presence_absence_file)
}
\arguments{
\item{blast_res}{a data frame with BLAST search results obtained by using
\code{\link{get_blast_res}}.}

\item{presence_absence_file}{A path to the gene presence absence file in
a csv format.}
}
\value{
A data frame with BLAST results expanded to individual genomes from
a pangenome.
}
\description{
This function extends blast results run on a pangenome to individual
genomes that were used for pangenome creation.
}
\seealso{
get_blast_res
}
