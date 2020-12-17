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

data(adhesins_df)

source("utils.R")

options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output, session) {
    
    all_systems <- unique(adhesins_df[["System"]])
    
    plot_cols <- reactiveValues(presence_col = "#e00e00",
                                absence_col = "#85c1ff")
    
    observeEvent(input[["presence_col"]], {
      plot_cols[["presence_col"]] <- input[["presence_col"]]
    })
    
    observeEvent(input[["absence_col"]], {
      plot_cols[["absence_col"]] <- input[["absence_col"]]
    })
    
    output[["input_tab"]] <- renderTable({
      validate(need(input[["seq_file"]], "Please upload your files in a FASTA format."))
      data.frame(`Uploaded files` = input[["seq_file"]][["name"]], check.names = FALSE)
    })
    
    blast_results <- eventReactive(input[["blast"]], {
      req(input[["seq_file"]])
      run_blast(input[["seq_file"]], input[["n_threads"]], updateProgress)
      
    })
    

    presence_tab <- reactive({
    #  validate(need(input[["blast"]], "Please run BLAST to see the results."))
        get_presence_table(blast_results(), 
                           identity_threshold = input[["identity"]], evalue_treshold = input[["evalue"]])
    })
    
    selected_systems <- reactive({
        req(presence_tab)
        input[["systems"]]
    })
    
    presence_plot_dat <- reactive({
      #req(presence_tab)
      validate(need(presence_tab, "Please run BLAST to see the results."))
        get_data_for_plots(presence_tab(), input[["systems"]])
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
      get_summary_plot(presence_tab(), hide_absent = input[["systems_summary_hide_missing"]])
    }, height = 200+10*nrow(summary_table()), width = 300+10*ncol(summary_table()))
  })  

  absent_systems <- reactive({
    summary_tab <- get_summary_table(presence_tab())
    names(summary_tab[, 2:ncol(summary_tab)][which(colSums(summary_tab[, 2:ncol(summary_tab)]) == 0)])
  })
  
    output[["wordcloud"]] <- renderPlot({
      get_word_cloud(get_count_table(presence_tab()))
    })
    
    
    all_genes_plot_dat <- reactive({
      get_presence_table(blast_results(), add_missing = !input[["all_genes_hide_missing"]]) %>% 
        get_data_for_plots()
    })
    
  observe({
      output[["presence_plot"]] <- renderPlot({
        req(all_genes_plot_dat)
     #   validate(need(is.null(presence_plot_dat), "Please run BLAST to see the results."))
          get_presence_plot(all_genes_plot_dat(),
                            presence_col = input[["presence_col"]], 
                            absence_col = input[["absence_col"]])
      }, height = 100+10*length(unique(all_genes_plot_dat()[["Gene"]])), width = 100+10*length(unique(all_genes_plot_dat()[["File"]])))
  })

  
  plot_system_dat <- reactive({
    req(presence_plot_dat)
    df <- left_join(presence_plot_dat(), adhesins_df, by = "Gene") 
    if(input[["systems_hide_missing"]] == TRUE) {
      df <- filter(df, !(System %in% absent_systems()))
    }
    ungroup(df)
  })
  
  observe({
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
                get_system_plot(system_data(), systems()[my_i], 
                                presence_col = input[["presence_col"]], 
                                absence_col = input[["absence_col"]])
              }, width = 320+10*nc(), height = 60+15*nr())
          })
      }
  })

  
})
