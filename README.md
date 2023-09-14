# adhesiomeR

adhesiomeR is a tool for analysis of adhesiome in *Escherichia coli* genomes.

## Installation 

### Installing BLAST+ 

AdhesiomeR requires BLAST+ >= 2.10.0 

AdhesiomeR uses BLAST to search for hits to adhesins. Adhesin database included 
within the package is compatible with BLAST+ version 2.10.0 or newer. 

To install BLAST, follow instructions available at [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569861/).

### Installing adhesiomeR

You can install the latest development version of the package using:

``` r
if(!require("devtools")) install.packages("devtools")
devtools::install_github("https://github.com/ksidorczuk/adhesiomeR/")
```

## adhesiomeR 

Our tool is meant for analysis of *Escherichia coli* adhesin repertoire. It uses
genome files in a fasta format as an input and searches for adhesins using BLAST.
It considers a gene as present using curated gene-specific bit score thresholds 
or identity percent and coverage thresholds (see [Strict and relaxed versions] for
more details). 

### Adhesin database

AdhesiomeR includes the most comprehensive database of adhesins to allow characterization
of the full E. coli adhesiomes with high confidence. It includes **527 adhesin genes**
combined into **102 systems**, where system indicates a whole operon by which given
adhesin is encoded in case of fimbrial adhesins, or a single gene in case of autotransporters
and other nonfimbrial adhesins encoded by a single gene.

### Strict and relaxed versions

We offer two modes of search to allow users more flexibility when searching for adhesins.

The **strict** mode is the version we recommend for general use for identification
of known adhesins with high confidence. In this mode, genes are considered as present
if they pass a specific bit score threshold. These thresholds have been set manually for
each gene based on a set of reference adhesins (usually, genes forming a complete operon
on a single contig), a similar strategy has been used to determine curated cut-offs in 
[CARD Resistance Gene Identifier](https://github.com/arpcard/rgi).

The **relaxed** mode is a more flexible version. Here, genes are considered as present
using 75% identity and 75% coverage thresholds by default, but they may be modified by 
the user. This approach is more suitable for searching for putative adhesins and 
adhesin-like proteins.


## Usage

### Example analysis using ECOR genomes

Example below describes how to perform adhesiome analysis of 72 *Escherichia coli* reference
collection (ECOR) strains. 

#### Step 1. Download genomes 

This step is optional. If you already have your own genome fasta files you would like
to analyse with adhesiomeR, proceed to the next step. 

``` r
# rentrez package is required for downloading data
if (!require("rentrez")) install.packages("rentrez")
if (!require("R.utils")) install.packages("R.utils")
library(rentrez)
library(R.utils)
library(dplyr)

# get genome accessions from data file
ecor_data <- read.csv("./inst/ecor_data.csv", check.names = FALSE)
ecor_accessions <- ecor_data[["GenBank accession no."]]

# create directory for storing downloaded genomes
if(!dir.exists("./ecor_genomes/")) dir.create("./ecor_genomes/")

# download ECOR genomes in fasta format
# please be patient - it may take a few minutes
for(ith_acc in ecor_accessions) {
  search_id <- entrez_search(db = "nuccore", term = ith_acc)[["ids"]][1]
  link_to_assembly <- entrez_link(dbfrom = "nuccore", db = "assembly", id = search_id)[["links"]][["nuccore_assembly"]]
  ftp_link <- entrez_summary(db = "assembly", id = link_to_assembly)[["ftppath_refseq"]]
  assembly_accession <- entrez_summary(db = "assembly", id = link_to_assembly)[["assemblyaccession"]]
  full_link <- paste0(ftp_link, "/", last(strsplit(ftp_link, "/")[[1]]), "_genomic.fna.gz")
  file_name <- last(strsplit(full_link, "/")[[1]])
  download.file(url = full_link, destfile = paste0("./ecor_genomes/", file_name))
  gunzip(paste0("./ecor_genomes/", file_name))
  print(paste0("Downloaded ", assembly_accession))
}

```


#### Step 2. Run adhesiomeR analysis to search for known adhesins

You have downloaded *E. coli* genomes in fasta format. Now you can start analysing them.

```r
library(adhesiomeR)

# first, you need a vector of filenames you wish to analyse
# to obtain it for genomes downloaded in the previous step, run:
genomes <- list.files("./ecor_genomes/", full.names = TRUE)

# run blast
# you can specify number of threads to using nt argument
blast_results <- get_blast_res(genomes, n_threads = 4)

# get gene presence information from blast results
# here, 1 indicates gene presence and 0 its absence
presence_df <- get_presence_table_strict(blast_results, n_threads = 8)

# get gene copy number information from blast results
# here, if gene is present in multiple copies, its occurences are counted
presence_df_copies <- get_presence_table_strict(blast_results, n_threads = 8, count_copies = TRUE)

# get system information
system_df <- get_summary_table(presence_df)

```
#### Step 3. Plot the results

```r
# plot gene presence
get_presence_plot(presence_df)

# plot system presence
get_summary_plot(system_df, hide_absent = TRUE)

# by setting hide_absent parameter you may remove systems not found in any genome from the plot
get_summary_plot(system_df, hide_absent = TRUE)

```