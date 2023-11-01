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
of the full E. coli adhesiomes with high confidence. It includes **525 adhesin genes**
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

### Pangenome analysis

To analyse pangenomes with adhesiomeR, you will need protein sequences of the pangenome
and the gene presence absence table. You can generate this files with Roary or
panaroo. Note that adhesiomeR functions were developed to work on output files of
these tools! If you have generated your pangenome using different software, please
make sure that your files are in the same format as output files from Roary or panaroo. 

Please note that while using pangenome file, you will not be able to check multiple
occurrences of genes, using ```count_copies = TRUE``` on expanded results will not
be able to identify gene copies as the information about their original localization
is not retained in these files. 

You can download example files below:

- [pan_genome_reference.fa](https://www.dropbox.com/scl/fi/d7ydmolwgiqamfsq3q5y5/pan_genome_reference.fa?rlkey=v1svmtujkxy5gdu9xjebkeveq&dl=0) 
- [gene_presence_absence.csv](https://www.dropbox.com/scl/fi/6livd7r4comkp5dmkwev2/gene_presence_absence.csv?rlkey=q876glanp2arunf9sitovcqij&dl=0)

```r
# Run blast search on pangenome fasta file
blast_res <- get_blast_res("/path/to/pangenome/files/pan_genome_reference.fa", n_threads = 4)

# Expand results into individual genomes
blast_res_full <- pangenome_to_genome(blast_res, "path/to/pangenome/files/gene_presence_absence.csv")

# Get adhesin presence/absence table
presence_table <- get_presence_table(blast_res_full)

# Get adhesin system information
system_table <- get_summary_table(presence_table)
```

### Metagenomics - gene catalogue analysis

This example uses files published in [Hildebrand et al. 2021](https://doi.org/10.1016/j.chom.2021.05.008).
Please be aware that the following type of analysis is meant to be run on linux machine.
First, we run BLAST search on a gene catalogue. In the next step, we want to associate each hit with
samples that contain given gene. Since abundance matrices generated from metagenomic WGS analyses are
generally too big to load them into R, we first extract the subset of data from the original file by 
selecting  accessions/numbers of genes with hits to adhesins. This subset of the original abundance
matrix can be loaded into R and used for further processing. 
Please note that this approach will not be able to identify multiple copies of adhesin genes in samples.

#### Step 1. Run BLAST search

```r
# Run blast search on a gene catalogue 
blast_res <- get_blast_res("/path/to/gene/catalogue/compl.incompl.95.prot.faa", n_threads = 12)
# Save blast results to csv file for further processing
write.csv(blast_res, "gc_blast_results.csv", row.names = FALSE)
```

#### Step 2. Extract subset of data from abundance matrix

Notice that this step does not use R but bash. You can also use these commands
from within R using ```system``` function. 

```
# Extract gene names with hits to adhesins
cut -d$'\t' -f 1 gc_blast_results.csv | uniq | grep -v "Query" > genes.txt

# Extract gene abundance in samples
zcat Matrix.mat.gz | head -n 1 > gene_abundance.txt
zcat Matrix.mat.gz | sed -nf <(sed 's/$/p/' genes.txt) >> gene_abundance.txt
```

#### Step 3. Expand results to samples

```r
# Define path to the subset of abundance matrix
ab_mat <- "gene_abundance.txt"

# Extend blast results
extended_blast_res <- gc_to_sample(blast_res, ab_mat)
```

#### Step 4. Run adhesiomeR analysis

```r
# Get adhesin gene presence/absence table
presence_table <- get_presence_table_strict(extended_blast_res, n_threads = 12)

# Get system summary
system_table <- get_summary_table(presence_table)
```
