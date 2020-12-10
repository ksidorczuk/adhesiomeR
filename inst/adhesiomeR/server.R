library(shiny)
library(DT)
library(adhesiomeR)
library(biogram)
library(dplyr)
library(ggplot2)
library(shinyWidgets)
library(tidyr)
library(wordcloud)

data(adhesins_df)
# 
# load("../data/adhesins_df.rda")
source("utils.R")
# source("../R/get_presence.R")
# source("../R/get_system_plot.R")
# source("../R/utils.R")

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
    
    output[["input_tab"]] <- renderText({
      validate(need(input[["seq_file"]], "Please upload your files in a FASTA format."))
      input[["seq_file"]][["name"]]
    })
    
      blast_results <- eventReactive(input[["blast"]], {
        validate(
          need(input[["seq_file"]], "Please provide a fasta file.")
        )
        progress <- shiny::Progress$new(min = 0, max = length(input[["seq_file"]][, 1]))
        progress$set(message = "Running BLAST", value = 0)
        on.exit(progress$close())
        
        updateProgress <- function(value = NULL, detail = NULL) {
          progress$set(value = value, detail = detail)
        }
        # if(length(input[["seq_file"]][, 1]) > 3) {
        #     stop("Too many files. You can analyze up to three genomes at once.")
        # }
        res <- run_blast(input[["seq_file"]], updateProgress)
        progress$set(value = progress[["getMax"]]())
        res

        # lapply(1:length(input[["seq_file"]][, 1]), function(i) {
        #   get_blast_res(input[["seq_file"]][[i, 4]]) %>% 
        #     mutate(File = input[["seq_file"]][[i, 1]])
        # }) %>% bind_rows()
      })
    


    presence_tab <- reactive({
        validate(need(blast_results, "Please run blast first."))
        get_presence_table(blast_results(), 
                           identity_threshold = input[["identity"]], evalue_treshold = input[["evalue"]])
    })
    
    selected_systems <- reactive({
        req(presence_tab)
        input[["systems"]]
    })
    
    presence_plot_dat <- reactive({
        req(presence_tab)
        get_data_for_plots(presence_tab(), input[["systems"]])
    })
    

    output[["blast_res"]] <- renderDataTable({
        validate(need(blast_results, "Please run blast first."))
         blast_results() %>% 
           mutate(Subject = sapply(Subject, function(i) strsplit(i, "~")[[1]][2])) %>% 
            my_DT()
    })
    
    
    output[["presence"]] <- renderDataTable({
        my_DT(presence_tab())
    })
    
    output[["systems_summary"]] <- renderDataTable({
      get_summary_table(presence_tab()) %>% 
        my_DT()
    })
    
    output[["wordcloud"]] <- renderPlot({
      get_word_cloud(get_count_table(presence_tab()))
    })
    
  observe({
      output[["presence_plot"]] <- renderPlot({
          get_presence_plot(presence_plot_dat(),
                            presence_col = input[["presence_col"]], 
                            absence_col = input[["absence_col"]])
      }, height = 300+10*nrow(presence_plot_dat()), width = 50+20*ncol(presence_plot_dat()))
  })

  
  plot_system_dat <- reactive({
      req(presence_plot_dat)
      left_join(presence_plot_dat(), adhesins_df, by = "Gene")
  })
  
  observe({
      output[["systems_plots"]] <- renderUI({
          nc <- reactive({ncol(plot_system_dat())})
          systems_plots_list <- lapply(1L:length(unique(plot_system_dat()[["System"]])), function(i) {
              list(plotOutput(paste0("systems_plot", i), height = 150+15*nc()))
          })
      })

      for(i in 1L:length(unique(plot_system_dat()[["System"]]))) {
          local({
              my_i <- i
              systems <- reactive({unique(plot_system_dat()[["System"]])})
              system_data <- reactive({
                filter(plot_system_dat(), System == systems()[[my_i]])
              })
              nr <- reactive({nrow(system_data())})
              nc <- reactive({ncol(system_data())})
              output[[paste0("systems_plot", my_i)]] <- renderPlot({
                get_system_plot(system_data(), systems()[[my_i]], 
                                presence_col = input[["presence_col"]], 
                                absence_col = input[["absence_col"]])
              }, width = 300+10*nr(), height = 100+15*nc())
          })
      }
  })

  
})
