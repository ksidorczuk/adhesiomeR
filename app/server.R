library(shiny)
library(DT)
library(biogram)
library(dplyr)
library(ggplot2)

load("../data/adhesins_df.rda")
source("utils.R")

options(shiny.maxRequestSize=10*1024^2)

shinyServer(function(input, output) {
    
    output[["tab"]] <- renderTable({
        req(input[["seq_file"]])
        input[["seq_file"]][["name"]]
    })
    
    output[["blast_res"]] <- renderDataTable({
        req(input[["seq_file"]])
        if(length(input[["seq_file"]][, 1]) > 3) {
            stop("Too many files. You can analyze up to three genomes at once.")
        }
        lapply(1:length(input[["seq_file"]][, 1]), function(i) {
            get_blast_res(input[["seq_file"]][[i, 4]]) %>% 
                mutate(File = input[["seq_file"]][[i, 1]])
        }) %>% bind_rows() %>% 
            my_DT()
    })
})
