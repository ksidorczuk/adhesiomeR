library(shiny)
library(DT)
library(adhesiomeR)
library(biogram)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(tidyr)
library(wordcloud)
library(foreach)
library(doParallel)
library(pander)
library(rmarkdown)

data(adhesins_df)

source("shiny-utils.R")

options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output, session) {
  
  all_systems <- unique(adhesins_df[["System"]])
  
  filters <- reactiveValues(thresh = 50,
                            evalue = 1e-50,
                            systems = "Type 1 Fimbriae",
                            presence_col = "#e00e00",
                            absence_col = "#85c1ff",
                            old_files = NULL)
  
  hide_tabs <- eventReactive(input[["blast"]], {
    if(input[["blast"]] == 0) {
      return(0)
    } else {
      return(1)
    }
  }, ignoreNULL = FALSE)

  observe({
    if(hide_tabs() == 0) {
    hideTab("adhesiomer", "blast_res")
    hideTab("adhesiomer", "summary_plot")
    hideTab("adhesiomer", "all_genes")
    hideTab("adhesiomer", "systems")
    hideTab("adhesiomer", "summary_plot")
    hideTab("adhesiomer", "report")
    } else {
      showTab("adhesiomer", "blast_res")
      showTab("adhesiomer", "summary_plot")
      showTab("adhesiomer", "all_genes")
      showTab("adhesiomer", "systems")
      showTab("adhesiomer", "summary_plot")
      showTab("adhesiomer", "report")
    }
  })
  
  
  observeEvent(input[["filtering_button"]], {
    if(length(input[["systems"]]) < 1) {
      showModal(modalDialog(
        title = "Select at least one system!",
        "You have to select at least one system for the filtering to work.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    }
    validate(need(length(input[["systems"]]) > 0, ""))
    filters[["thresh"]] <- input[["identity"]]
    filters[["evalue"]] <- input[["evalue"]]
    filters[["systems"]] <- input[["systems"]]
    filters[["presence_col"]] <- input[["presence_col"]]
    filters[["absence_col"]] <- input[["absence_col"]]
  })
  
  
  output[["adhesins"]] <- renderDataTable({
    my_DT(adhesiomeR::adhesins_df)
  })
  
  output[["input_tab"]] <- renderTable({
    validate(need(input[["seq_file"]], "Please upload your files in a FASTA format."))
    data.frame(`Uploaded files` = input[["seq_file"]][["name"]], check.names = FALSE)
  })
  
  
  blast_results <- eventReactive(input[["blast"]], {
    req(input[["seq_file"]])
    if(length(input[["seq_file"]][["name"]]) > 100) {
      showModal(modalDialog(
        title = "Too many files!",
        "You can analyse up to 100 files at once using the GUI. For larger analyses please use
          the command-line interface of the adhesiomeR package.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    }
    if(all(input[["seq_file"]][["name"]] %in% filters[["old_input"]][["name"]]) &
       length(input[["seq_file"]][["name"]] == length(filters[["old_input"]][["name"]]))) {
      showModal(modalDialog(
        title = "You already have results for these files",
        "Uploaded files did not change since the last BLAST analysis. If you wish to analyse more 
          or different files, please be sure to upload them first.",
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
    } 
    req(!(all(input[["seq_file"]][["name"]] %in% filters[["old_input"]][["name"]])) &
          length(input[["seq_file"]][["name"]] != length(filters[["old_input"]][["name"]])) &
          length(input[["seq_file"]][["name"]]) <= 100)
    showModal(modalDialog(
      title = "Running BLAST...",
      "Please be patient - the calculations may take a few minutes.
        This window will disapear once calculations are completed.", 
      footer = NULL))

    res <- run_blast(input[["seq_file"]], input[["n_threads"]])
    removeModal()
    filters[["old_input"]] <- input[["seq_file"]]
    res
  })
  
  
  presence_tab <- reactive({
    #  validate(need(input[["blast"]], "Please run BLAST to see the results."))
    get_presence_table(blast_results(), 
                       identity_threshold = filters[["thresh"]], evalue_threshold = filters[["evalue"]])
  })
  
  
  presence_plot_dat <- reactive({
    validate(need(presence_tab, "Please run BLAST to see the results."))
    adhesiomeR:::get_presence_plot_data(presence_tab(), filters[["systems"]])
  })
  
  
  output[["blast_res"]] <- renderDataTable({
    #validate(need(nrow(blast_results() > 0), "Please run BLAST to see the results."))
    blast_results() %>% 
      mutate(Subject = sapply(Subject, function(i) strsplit(i, "~")[[1]][2])) %>% 
      my_DT()
  })
  
  
  # output[["presence_table"]] <- renderDataTable({
  #     my_DT(presence_tab())
  # })
  
  # output[["systems_summary_table"]] <- renderDataTable({
  #   get_summary_table(presence_tab()) %>% 
  #     my_DT()
  # })
  summary_table <- reactive({
    get_summary_table(presence_tab(), hide_absent = input[["systems_summary_hide_missing"]])
  })
  
  observe({
    output[["systems_summary_plot"]] <- renderPlot({
      # validate(need(is.null(blast_results), "Please run BLAST to see the results."))
      get_summary_plot(presence_tab(), hide_absent = input[["systems_summary_hide_missing"]],
                       presence_col = filters[["presence_col"]], 
                       absence_col = filters[["absence_col"]])
    }, height = 200+10*nrow(summary_table()), width = 300+10*ncol(summary_table()))
  })  
  
  absent_systems <- reactive({
    summary_tab <- get_summary_table(presence_tab())
    names(summary_tab[, 2:ncol(summary_tab)][which(colSums(summary_tab[, 2:ncol(summary_tab)]) == 0)])
  })

  
  
  all_genes_plot_dat <- reactive({
    get_presence_table(blast_results(), add_missing = !input[["all_genes_hide_missing"]],
                       identity_threshold = filters[["thresh"]], evalue_threshold = filters[["evalue"]])
  })
  
  scaling_dat <- reactive({
    req(all_genes_plot_dat)
    adhesiomeR:::get_presence_plot_data(all_genes_plot_dat(), systems = filters[["systems"]])
  })
  
  observe({
    req(all_genes_plot_dat)
    output[["presence_plot"]] <- renderPlot({
      #   validate(need(is.null(presence_plot_dat), "Please run BLAST to see the results."))
      get_presence_plot(all_genes_plot_dat(),
                        systems = filters[["systems"]],
                        presence_col = filters[["presence_col"]], 
                        absence_col = filters[["absence_col"]])
    }, height = 300+10*length(unique(scaling_dat()[["Gene"]])), width = 100+10*length(unique(scaling_dat()[["File"]])))
  })
  
  
  
  
  observe({
    plot_system_dat <- reactive({
      #req(presence_plot_dat)
      df <- left_join(presence_plot_dat(), adhesins_df, by = "Gene") 
      if(input[["systems_hide_missing"]] == TRUE) {
        df <- filter(df, !(System %in% absent_systems()))
      }
      ungroup(df)
    })
    
    output[["systems_plots"]] <- renderUI({
      nc <- reactive({ncol(plot_system_dat())})
      systems_plots_list <- lapply(1L:length(unique(plot_system_dat()[["System"]])), function(i) {
        list(plotOutput(paste0("systems_plot", i), width = 150+15*nc()))
      })
    })
    
    systems <- reactive({unique(plot_system_dat()[["System"]])})
    
    for(i in 1L:length(unique(plot_system_dat()[["System"]]))) {
      local({
        my_i <- i
        system_data <- reactive({
          filter(plot_system_dat(), System == systems()[my_i])
        })
        nr <- reactive({length(unique(system_data()[["File"]]))})
        nc <- reactive({length(unique(system_data()[["Gene"]]))})
        output[[paste0("systems_plot", my_i)]] <- renderPlot({
          get_system_plot(presence_tab(), systems()[my_i], 
                          presence_col = filters[["presence_col"]], 
                          absence_col = filters[["absence_col"]])
        }, width = 320+10*nc(), height = 60+15*nr())
      })
    }
  })
  
  
  output[["download"]] <- downloadHandler(
    filename = "adhesiomeR-results.html",
    content <- function(file) {
      
      src <- normalizePath("adhesiomeR-report.Rmd")
      
      input_files <- input[["seq_file"]][["name"]]
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      adhesiomeR:::generate_report_files(presence_table = presence_tab(), elements = input[["elements"]], outdir = owd,
                                         hide_absent_genes = input[["report_hide_genes"]], hide_absent_systems = input[["report_hide_systems"]],
                                         presence_col = filters[["presence_col"]], absence_col = filters[["absence_col"]])
      file.copy(src, "adhesiomeR-report.Rmd", overwrite = TRUE)
      genome_files <- input[["seq_file"]][["name"]]
      outdir <- owd
      elements <- input[["elements"]]
      out <- rmarkdown::render("adhesiomeR-report.Rmd", output_format = "html_document", 
                        file, quiet = TRUE,
                        params = list(genome_files, outdir, elements))
      
      fl <- list.files(outdir, full.names = TRUE)
      sapply(c("summary_table.csv", "summary_plot.png", "presence_table.csv", "presence_plot.png"), function(i) 
        invisible(file.remove(fl[grepl(i, fl)])))
      file.rename(out, file)
    }
  )
  
})
