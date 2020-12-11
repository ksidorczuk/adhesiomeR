library(shiny)
library(DT)
library(shinyWidgets)
library(adhesiomeR)

data(adhesins_df)

shinyUI(navbarPage("AdhesiomeR",
                   tabPanel("Introduction",
                            "some text about adhesins, our database, and why adhesiomeR is awesome"),
                   tabPanel("Input & settings",
                            sidebarPanel(width = 6,
                                         h5(tags$b("Step 1.")),
                                         fileInput("seq_file", 
                                                   "Upload your genome file(s) in a FASTA format. You can process up to x files at once.",
                                                   multiple = TRUE
                                         ),
                                         actionButton("blast", "Run BLAST!",
                                                      style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                         HTML("<hr>"),
                                         h5(tags$b("Step 2.")),
                                         pickerInput("systems",
                                                     "Select systems you want to analyse:",
                                                     choices = unique(adhesins_df[["System"]]),
                                                     selected = "Type 1 Fimbriae",
                                                     multiple = TRUE,
                                                     pickerOptions(actionsBox = TRUE,
                                                                   selectedTextFormat = 'count',
                                                                   countSelectedText = '{0} systems selected')),
                                         HTML("<hr>"),
                                         h5(tags$b("Step 3.")),
                                         h5(tags$b("Select thresholds for gene presence: ")),
                                         fluidRow(column(4,
                                                         numericInput("identity", "Identity [%]",
                                                                      value = 70,
                                                                      min = 1,
                                                                      max = 100,
                                                                      step = 1,
                                                                      width = "100px")),
                                                  column(4,
                                                         numericInput("evalue", "E-value",
                                                                      value = prettyNum(1e-50, scientific = TRUE),
                                                                      min = 0,
                                                                      max = 1,
                                                                      width = "100px"))),
                                         HTML("<hr>"),
                                         h5(tags$b("Step 4.")),
                                         h5(tags$b("Customise plot colors: ")),
                                         fluidRow(column(4,
                                                         textInput("presence_col", "Presence",
                                                                   value = "#e00e00",
                                                                   width = "100px")),
                                                  column(4,
                                                         textInput("absence_col", "Absence",
                                                                      value = "#85c1ff",
                                                                      width = "100px")))
                            ),
                            mainPanel(width = 6,
                                      textOutput("input_tab"))
                   ),
                   
                   tabPanel("Search results",
                            dataTableOutput("blast_res")),
                   tabPanel("Summary table",
                            dataTableOutput("presence_table"),
                            dataTableOutput("systems_summary_table")),
                   tabPanel("Summary plot",
                            plotOutput("systems_summary_plot"),
                            plotOutput("wordcloud")),
                   tabPanel("All genes plot",
                            plotOutput("presence_plot")),
                   tabPanel("Systems plots",
                            uiOutput("systems_plots"))
                   
                   
))
