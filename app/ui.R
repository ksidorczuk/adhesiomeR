library(shiny)
library(DT)

shinyUI(fluidPage(

    titlePanel("AdhesiomeR"),
    sidebarLayout(
        sidebarPanel(
            fileInput("seq_file", 
                      "Upload your genome file(s) in a FASTA format. You can process up to 3 files at once.",
                      multiple = TRUE
                      ),
            checkboxGroupInput("systems",
                               "Select systems to plot",
                               choices = unique(adhesins_df[["systems"]]),
                               selected = "Type_1_fimbriae"),
            tableOutput("tab")
        ),

        mainPanel(
            tabsetPanel(
                tabPanel("Full blast results",
                         dataTableOutput("blast_res")),
                tabPanel("Systems"
                         )
            )
        )
    )
))
