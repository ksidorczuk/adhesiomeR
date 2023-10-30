library(shiny)
library(DT)
library(shinyWidgets)
library(adhesiomeR)
library(shinycssloaders)
library(colourpicker)
library(shinythemes)

data(adhesins_df)

shinyUI(tagList(
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="https://cdn.datatables.net/1.12.1/css/jquery.dataTables.css")),
  navbarPage(theme = shinytheme("flatly"),
             title = "adhesiomeR",
             id = "adhesiomer",
             tabPanel("Adhesin database",
                      wellPanel(
                        includeMarkdown("intro.md")),
                      tabsetPanel(
                        tabPanel(h4("Genes"),
                                 wellPanel(
                                   includeMarkdown(
                                     "genes_table.md")),
                                 dataTableOutput("adhesins")),
                        tabPanel(h4("Systems"),
                                 wellPanel(
                                   includeMarkdown(
                                     "systems_table.md")),
                                 dataTableOutput("systems_df")))),
             tabPanel("Input & settings",
                      sidebarPanel(width = 4,
                                   fileInput("seq_file", 
                                             "Upload your genome file(s) in a FASTA format. You can process up to 10 files at once.",
                                             multiple = TRUE
                                   ),
                                   actionButton("blast", "Run analysis!"),
                                   HTML("<hr>"),
                                   h5(tags$b("Customise plot colors: ")),
                                   fluidRow(column(4,
                                                   colourInput("presence_col", "Presence",
                                                               value = "#e00e00")),
                                            column(4,
                                                   colourInput("absence_col", "Absence",
                                                               value = "#85c1ff"))),
                                   actionButton("filtering_button", "Apply changes!"),
                      ),
                      mainPanel(width = 8,
                                dataTableOutput("input_tab"),
                                textOutput("hidetabs"))
             ),
             
             tabPanel("Search results",
                      value = "blast_res",
                      dataTableOutput("blast_res")),
             tabPanel("Gene presence",
                      value = "all_genes",
                      tabsetPanel(
                        tabPanel("Plot",
                                 checkboxInput("all_genes_hide_missing",
                                               "Hide genes not found in any genome.",
                                               value = FALSE),
                                 plotOutput("presence_plot")),
                        tabPanel("Table",
                                 dataTableOutput("presence_table")))),
             tabPanel("System presence",
                      value = "summary_plot",
                      tabsetPanel(
                        tabPanel("Plot",
                                 checkboxInput("systems_summary_hide_missing",
                                               "Hide systems not found in any genome.",
                                               value = FALSE),
                                 plotOutput("systems_summary_plot")),
                        tabPanel("Table",
                                 dataTableOutput("systems_summary_table"))
                      )),
             tabPanel("Gene/system presence",
                      value = "systems",
                      checkboxInput("systems_hide_missing",
                                    "Hide systems not found in any genome.",
                                    value = TRUE),
                      uiOutput("systems_plots")),
             tabPanel("Profiles",
                      value = "profiles",
                      dataTableOutput("profile_table")),
             tabPanel("Clusters",
                      value = "clusters",
                      dataTableOutput("cluster_table")),
             tabPanel("Report",
                      value = "report",
                      checkboxGroupInput("elements",
                                         width = 900,
                                         "Select which elements you want to include in the report:",
                                         choices = c("Table with percents of genes found in each system" = "summary_table",
                                                     "Plot with percents of genes found in each system" = "summary_plot",
                                                     "Table with gene presence/absence" = "presence_table",
                                                     "Plot with gene presence/absence" = "presence_plot"),
                                         selected = c("summary_table", "summary_plot", "presence_table", "presence_plot")),
                      h5(tags$b("Options: ")),
                      checkboxInput("report_hide_genes",
                                    width = 900,
                                    "Hide genes that were not found in any genome.",
                                    value = FALSE),
                      checkboxInput("report_hide_systems",
                                    width = 900,
                                    "Hide systems that were not found in any genome.",
                                    value = FALSE),
                      downloadButton("download",
                                     "Download report"))
             
             
  ),
  tags$style(type='text/css',
        #      HTML(".dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing,.dataTables_wrapper .dataTables_paginate .paginate_button, .dataTables_wrapper .dataTables_paginate .paginate_button.disabled {
        #     color: #0000ff !important;
        # }"),
             HTML(".dataTables_wrapper .dataTables_paginate .paginate_button {color: #4A7443}"))
)
)
