library(shiny)
library(DT)

shinyUI(navbarPage("AdhesiomeR",
                   tabPanel("Introduction",
                            "some text"),
                   tabPanel("Input & settings",
                            sidebarPanel(width = 5,
                                fileInput("seq_file", 
                                          "Upload your genome file(s) in a FASTA format. You can process up to 3 files at once.",
                                          multiple = TRUE
                                ),
                                pickerInput("systems",
                                            "Select systems to plot",
                                            choices = unique(adhesins_df[["System"]]),
                                            selected = "Type_1_fimbriae",
                                            multiple = TRUE,
                                            pickerOptions(actionsBox = TRUE,
                                                          selectedTextFormat = 'count',
                                                          countSelectedText = '{0} systems selected'))
                            ),
                            mainPanel(width = 6,
                                tableOutput("tab"))
                   ),
                   
                   tabPanel("Search results",
                            dataTableOutput("blast_res")),
                   tabPanel("Summary",
                            dataTableOutput("presence")),
                   tabPanel("All genes plot",
                            plotOutput("presence_plot"))
                   
                   
))
