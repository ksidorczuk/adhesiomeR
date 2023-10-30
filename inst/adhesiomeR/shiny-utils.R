options(DT.options = list(dom = "Brtip",
                          buttons = c("copy", "csv", "excel", "print"),
                          pageLength = 50, autoWidth = TRUE, bAutoWidth = FALSE
))

my_DT <- function(x, ...)
  datatable(x, ..., escape = FALSE, extensions = 'Buttons', filter = "top", rownames = FALSE,
            style = "bootstrap")




run_blast <- function(input_files, nt, updateProgress = NULL) {
  
    res <- lapply(1:nrow(input_files), function(i) {
      
      validate_input_file(input_files[[i, 4]])
      input_seqs <- read_fasta(input_files[[i, 4]])
      input <- tempfile(tmpdir = getwd())
      output <- tempfile(tmpdir = getwd())
      db_path <- paste0(normalizePath(system.file(package = "adhesiomeR")), "/db/adhesins")
      blast_command <- paste0("blastn -db ", db_path, " -query ", input, " -out ", output, " -outfmt 6")
      write_fasta(input_seqs, input)
      system(blast_command)
      
      res <- read.delim(output, header = FALSE)
      colnames(res) <- c("Query", "Subject", "% identity", "Alignment length", "Mismatches",
                         "Gap opens", "Query start", "Query end", "Subject start", "Subject end", "Evalue", "Bit score")
      
      file.remove(input, output)
      mutate(res, File = input_files[[i, 1]],
             Subject = sapply(Subject, function(i) strsplit(i, "~~~")[[1]][2]))
    }) %>% bind_rows() 
    
  res
  
}


run_analysis <- function(input_files, updateProgress = NULL) {
  blast_res <- run_blast(input_files)
 # print(paste0("blast done"))
  presence_table <- get_presence_table_strict(blast_res, n_threads = 1)
 # print(paste0("presence done"))
  return(list("blast_results" = blast_res, 
              "presence_table" = presence_table))
}


js <- c(
  "function(){", 
  "  this.api().columns().every(function(i){",
  "    var column = this;",
  "    var select = $('<select multiple=\"multiple\"><option value=\"\"></option></select>')",
  "      .appendTo( $('#th'+i).empty() )", 
  "      .on('change', function(){",
  "        var vals = $('option:selected', this).map(function(index,element){",
  "          return $.fn.dataTable.util.escapeRegex($(element).val());",
  "        }).toArray().join('|');",
  "        column.search(vals.length > 0 ? '^('+vals+')$' : '', true, false).draw();",
  "      });",
  "    var data = column.data();",
  "    if(i == 0){",
  "      data.each(function(d, j){",
  "        select.append('<option value=\"'+d+'\">'+d+'</option>');",
  "      });",
  "    }else{",
  "      data.unique().sort().each(function(d, j){",
  "        select.append('<option value=\"'+d+'\">'+d+'</option>');",
  "      });",
  "    }",
  "    select.select2({width: '100%', closeOnSelect: false});",
  "  });",
  "}")
