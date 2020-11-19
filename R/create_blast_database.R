create_blast_database <- function(sequence_file, dbname, blast_path = "/usr/bin") {
  command <- paste0(blast_path, "/makeblastdb -in ", sequence_file, " -dbtype nucl -out db/", dbname, " -title ", dbname)
  system(command)
}
