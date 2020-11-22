library(shiny)
library(DT)
library(biogram)
library(dplyr)
library(ggplot2)
library(shinyWidgets)

load("../data/adhesins_df.rda")
source("utils.R")
source("../R/get_presence.R")

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
        req(input[["seq_file"]])
        input[["seq_file"]][["name"]]
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
      }, height = 300+10*nrow(presence_plot_dat()), width = 50+10*ncol(presence_plot_dat()))
  })

})
