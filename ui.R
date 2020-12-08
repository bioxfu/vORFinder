suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(DECIPHER))
suppressPackageStartupMessages(library(Gviz))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinycssloaders))

ui <- dashboardPage(
  dashboardHeader(title = "Virus ORF Alignment", disable = FALSE),
  dashboardSidebar(
    fileInput("file", "input virus genome files(.zip)"),
    textInput("minorflen", "Minimal ORF length (nt) for ORF finder:", "50"),
    textInput("evalue", "Expectation value (E) threshold for BLAST:", "1e-10"),
    #textInput("raw_name", "full name of virus genus", "Tomato yellow leaf curl"),
    #textInput("new_name", "short name of virus genus:", "TYLC"),
  
    div(style = "padding-left: 15px; padding-top: 70px;",
        h4("contact "),
        p(class="small"," "),
        p(class = "small","Xing Fu Ph.D. "),
        p(class = "small","Bioinfomatics Core Facility")
    )
  ),
  dashboardBody(
    fluidRow(
      tabBox(
        width = 12, height = NULL,
        tabPanel(
          "workflow",
          h3("Work Flow"),
          fluidRow(
            img(src="workflow", height="70%", width="60%"),
            verbatimTextOutput("preparation")
          )
        ),
        tabPanel(
          "orfdatabase", 
          h3("ORF database"),
          fluidRow(
            box(
              width = 5,
              title = "ORF database:",
              withSpinner(DT::dataTableOutput("orfdatabase"), type = 5)
            ),
            box(
              width = 7,
              title = "Phylogenetic Tree",
              solidHeader = TRUE,
              collapsible = TRUE, 
              withSpinner(plotOutput("plot1", height = "880px"), type = 5)  
            )
          )
        ),
        tabPanel(
          "blast",
          h3("BLAST"),
          fluidRow(
            box(
              width = 12,
              title = "ORF Viewer",
              solidHeader = TRUE,
              collapsible = TRUE, 
              withSpinner(plotOutput("plot2", height = "440px"), type = 5)
            ),
            box(
              width = 12,
              collapsible = TRUE,
              withSpinner(DT::dataTableOutput("blastp"), type = 5)
            )
          )
        ),
        tabPanel(
          "Hits Alignment & ORF Distribution",
          fluidRow(
            box(downloadButton("downloadHTML", "DownloadHTML")),
            box(
              width = 12,
              title = "ORF Hits Alignment", 
              collapsible = TRUE,
              solidHeader = TRUE,
              withSpinner(htmlOutput("msa"), type = 5)
            )
          ),
          fluidRow(
            h4("Phylo Tree"),
            box(downloadButton("downloadPDF", "DownloadPDF")),
            box(
              width = 12,
              title = "Phylogenetic Tree",
              solidHeader = TRUE,
              collapsible = TRUE,
              withSpinner(uiOutput("phylo"), type = 5)
            )
          )
        )
      )
    )
  )
)


