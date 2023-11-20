library(shiny)
library(DT)
library(adhesiomeR)
library(biogram)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(tidyr)
library(pander)
library(rmarkdown)
library(future.apply)
library(gridExtra)
library(egg)
library(parallel)

pkgload::load_all()
source("shiny-utils.R")




options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output, session) {
  
  observe_helpers(withMathJax = TRUE, help_dir = ".")
  
  all_systems <- unique(adhesins_df[["System"]])
  
  filters <- reactiveValues(systems = unique(adhesins_df[["System"]]),
                            presence_col = "#e00e00",
                            absence_col = "#85c1ff",
                            old_files = NULL)
  
  
  hide_tabs <- reactiveVal(0)
  
  observe({
    if(hide_tabs() == 0) {
      hideTab("adhesiomer", "blast_res")
      hideTab("adhesiomer", "summary_plot")
      hideTab("adhesiomer", "all_genes")
      hideTab("adhesiomer", "systems")
      hideTab("adhesiomer", "summary_plot")
      hideTab("adhesiomer", "clusters")
      hideTab("adhesiomer", "profiles")
      hideTab("adhesiomer", "report")
    } else {
      showTab("adhesiomer", "blast_res")
      showTab("adhesiomer", "summary_plot")
      showTab("adhesiomer", "all_genes")
      showTab("adhesiomer", "systems")
      showTab("adhesiomer", "summary_plot")
      showTab("adhesiomer", "clusters")
      showTab("adhesiomer", "profiles")
      showTab("adhesiomer", "report")
    }
  })
  
  
  observeEvent(input[["filtering_button"]], {
    filters[["presence_col"]] <- input[["presence_col"]]
    filters[["absence_col"]] <- input[["absence_col"]]
  })
  
  
  output[["adhesins"]] <- renderDataTable({
    my_DT(df_genes)
  })
  
  output[["input_tab"]] <- renderDataTable({
    validate(need(input[["seq_file"]], "Please upload your files in a FASTA format."))
    my_DT(data.frame(`Uploaded files` = input[["seq_file"]][["name"]], check.names = FALSE))
  })
  
  
  analysis_results <- eventReactive(input[["blast"]], {
      req(input[["seq_file"]])
      if(length(input[["seq_file"]][["name"]]) > 10) {
        showModal(modalDialog(
          title = "Too many files!",
          "You can analyse up to 10 files at once using the GUI. For larger analyses please use
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
            length(input[["seq_file"]][["name"]]) <= 10)
      showModal(modalDialog(
        title = "Running analysis...",
        "Please be patient - the calculations may take a few minutes.
        This window will disappear once calculations are completed.", 
        footer = NULL))
      
      res <- run_analysis(input[["seq_file"]], input[["n_threads"]])
      hide_tabs(1)
      filters[["old_input"]] <- input[["seq_file"]]
      removeModal()
      res
    })
   

  
  output[["systems_df"]] <- renderDataTable({
    my_DT(df_systems, options = list(initComplete = JS(js)))
  })

  presence_tab <- reactive({
    #  validate(need(input[["blast"]], "Please run BLAST to see the results."))
    analysis_results()[[2]]
    # get_presence_table(blast_results())
  })
  
  
  presence_plot_dat <- reactive({
    validate(need(presence_tab, "Please run BLAST to see the results."))
    get_presence_plot_data(presence_tab(), filters[["systems"]])
  })
  
  clustering_plot_dat <- reactive({
    get_clustering_plot_data(presence_tab())
  })
  
  output[["blast_res"]] <- renderDataTable({
    #validate(need(nrow(blast_results() > 0), "Please run BLAST to see the results."))
    my_DT(analysis_results()[[1]])
    #my_DT(data.frame(blast_results()))
  })
  
  output[["profile_table"]] <- renderDataTable({
    my_DT(get_adhesin_profiles(presence_tab()), 
          options = list(columnDefs = list(list(
            targets = 1:2, 
            render = JS(
            "function(data, type, row, meta) {",
            "return data === null ? 'NA' : data;",
            "}")))))
  })
  
  output[["cluster_table"]] <- renderDataTable({
    my_DT(get_adhesin_clusters(presence_tab()))
  })
  
  output[["presence_table"]] <- renderDataTable({
    my_DT(presence_tab())
  })
  
  output[["systems_summary_table"]] <- renderDataTable({
    get_summary_table(presence_tab()) %>%
      my_DT()
  })
  
  
  summary_table <- reactive({
    get_summary_table(presence_tab(), hide_absent = input[["systems_summary_hide_missing"]])
  })
  
  observe({
    output[["systems_summary_plot"]] <- renderPlot({
      # validate(need(is.null(blast_results), "Please run BLAST to see the results."))
      get_summary_plot(
        presence_tab(), 
        hide_absent = input[["systems_summary_hide_missing"]],
        presence_col = filters[["presence_col"]],
        absence_col = filters[["absence_col"]])
    }, height = 150+10*nrow(summary_table()), width = 300+10*ncol(summary_table()))
  })
  
  absent_systems <- reactive({
    summary_tab <- get_summary_table(presence_tab())
    absent <- colnames(summary_tab)[sapply(colnames(summary_tab), function(i) ifelse(all(summary_tab[[i]] == "Absent"), TRUE, FALSE))]
    if("CS27A/CS27B" %in% absent) {
      absent <- c(absent, "CS27A", "CS27B")
    }
    absent
  })
  
  
  all_genes_plot_dat <- reactive({
    if(input[["all_genes_hide_missing"]]) {
      presence_tab()[c(1, which(colSums(presence_tab()[2:ncol(presence_tab())]) > 0) + 1)]
    } else {
      presence_tab()
    }
  })
  
  scaling_dat <- reactive({
    req(all_genes_plot_dat)
    get_presence_plot_data(all_genes_plot_dat(), systems = filters[["systems"]])
  })
  
  observe({
    req(all_genes_plot_dat)
    output[["presence_plot"]] <- renderPlot({
      #   validate(need(is.null(presence_plot_dat), "Please run BLAST to see the results."))
      get_presence_plot(all_genes_plot_dat(),
                        systems = filters[["systems"]],
                        presence_col = filters[["presence_col"]],
                        absence_col = filters[["absence_col"]])
    }, height = 300+10*length(unique(scaling_dat()[["Gene"]])), width = 200+10*length(unique(scaling_dat()[["File"]])))
  })
  
  
  
  
  observe({
    plot_system_dat <- reactive({
      df <- left_join(presence_plot_dat(), adhesins_df_grouped, by = "Gene", relationship = "many-to-many")
      if(input[["systems_hide_missing"]] == TRUE) {
        df <- filter(df, !(System %in% absent_systems()))
      }
      ungroup(df)
    })
    
    output[["systems_plots"]] <- renderUI({
      nc <- reactive({250 + 15*(length(unique(plot_system_dat()[["File"]])))})
      systems_plots_list <- lapply(1L:length(unique(plot_system_dat()[["System"]])), function(i) {
        list(withSpinner(plotOutput(paste0("systems_plot", i), height = nc())))
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
          grid.arrange(set_panel_size(
            get_system_plot(presence_tab(), systems()[my_i],
                            presence_col = filters[["presence_col"]],
                            absence_col = filters[["absence_col"]]),
            width = unit(20*nc(), "pt"), height = unit(20*nr(), "pt")))
        }, width = 350+10*nc(), height = 200+15*nr())
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
      generate_report_files(presence_table = presence_tab(), elements = input[["elements"]], outdir = owd,
                            hide_absent_genes = input[["report_hide_genes"]], hide_absent_systems = input[["report_hide_systems"]],
                            presence_col = filters[["presence_col"]], absence_col = filters[["absence_col"]])
      file.copy(src, "adhesiomeR-report.Rmd", overwrite = TRUE)
      genome_files <- input[["seq_file"]][["name"]]
      outdir <- owd
      elements <- input[["elements"]]
      out <- rmarkdown::render("adhesiomeR-report.Rmd", output_format = "html_document",
                               file, quiet = TRUE,
                               params = list(genome_files = genome_files, outdir = outdir, selected_elements = elements, presence_table = presence_tab()))
      
      fl <- list.files(outdir, full.names = TRUE)
      sapply(c("summary_table.csv", "summary_plot.png", "presence_table.csv", "presence_plot.png", "cluster_table.csv", "profile_table.csv"), function(i)
        invisible(file.remove(fl[grepl(i, fl)])))
      file.rename(out, file)
    }
  )
  
})
