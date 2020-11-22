library(shiny)
library(DT)
library(biogram)
library(dplyr)
library(ggplot2)
library(shinyWidgets)

load("../data/adhesins_df.rda")
source("utils.R")
source("../R/get_presence.R")
source("../R/get_system_plot.R")

options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output, session) {
    
    all_systems <- unique(adhesins_df[["System"]])


    
    blast_results <- reactive({
        req(input[["seq_file"]])
        if(length(input[["seq_file"]][, 1]) > 3) {
            stop("Too many files. You can analyze up to three genomes at once.")
        }
        lapply(1:length(input[["seq_file"]][, 1]), function(i) {
            get_blast_res(input[["seq_file"]][[i, 4]]) %>% 
                mutate(File = input[["seq_file"]][[i, 1]])
        }) %>% bind_rows()
    })
    


    presence_tab <- reactive({
        req(blast_results)
        get_presence_table(blast_results())
    })
    
    selected_systems <- reactive({
        req(presence_tab)
        input[["systems"]]
    })
    
    presence_plot_dat <- reactive({
        req(presence_tab)
        get_data_for_plots(presence_tab(), input[["systems"]])
    })
    
    output[["tab"]] <- renderTable({
        presence_plot_dat()
        # req(input[["seq_file"]])
        # input[["seq_file"]][["name"]]
    })
    
    output[["blast_res"]] <- renderDataTable({
        req(blast_results)
         blast_results() %>% 
            my_DT()
    })
    
    
    output[["presence"]] <- renderDataTable({
        my_DT(presence_tab())
    })
    
  observe({
      output[["presence_plot"]] <- renderPlot({
          get_presence_plot(presence_plot_dat())
      }, height = 300+10*nrow(presence_plot_dat()), width = 50+20*ncol(presence_plot_dat()))
  })

  
  plot_system_dat <- reactive({
      req(presence_plot_dat)
      left_join(presence_plot_dat(), adhesins_df, by = "Gene")
  })
  
  observe({
      output[["systems_plots"]] <- renderUI({
          systems_plots_list <- lapply(1L:length(unique(plot_system_dat()[["System"]])), function(i) {
              list(plotOutput(paste0("systems_plot", i)))
          })
      })

      for(i in 1L:length(unique(plot_system_dat()[["System"]]))) {
          local({
              my_i <- i
              systems <- reactive({unique(plot_system_dat()[["System"]])})
              output[[paste0("systems_plot", my_i)]] <- renderPlot({
                  filter(plot_system_dat(), System == systems()[[my_i]]) %>%
                      get_system_plot(., systems()[[my_i]])
              })
          })
      }
  })

  
})
