get_genes_in_systems_db <- function(seq_file) {
  data_file <- readLines(seq_file)
  def_lines <- data_file[grepl(">", data_file)]
  df <- data.frame(Gene = sapply(def_lines, function(x) strsplit(x, "~")[[1]][2], USE.NAMES = FALSE),
                   System = sapply(def_lines, function(x) strsplit(x, "~")[[1]][4], USE.NAMES = FALSE),
                   stringsAsFactors = FALSE) %>% 
    mutate(System = gsub("^ ", "", System))
  # unique_systems <- unique(df[["systems"]]) %>% 
  #   sapply(., function(x) strsplit(x, " ")[[1]], USE.NAMES = FALSE) %>% 
  #   unlist() %>% 
  #   unique()
  # lapply(unique_systems, function(ith_system) {
  #   genes <- filter(df, systems == ith_system)[["gene"]] %>% 
  #     paste0(collapse = ", ")
  #   data.frame(System = ith_system,
  #              Genes = genes,
  #              stringsAsFactors = FALSE)
  # }) %>% bind_rows()
  df
}
