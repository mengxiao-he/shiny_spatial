#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
# Creeate a parameter tabs
# gene_name <- rownames(HD)
# Seurat_cluster <- HD@meta.data$seurat_cluster.projected
# Banksy_cluster <- HD@meta.data$banksy_cluster


# parameter_tabs <- tabsetPanel(
#   id = "params",
#   type = "hidden",
#   tabPanel("Transcript",
#            selectInput("gene_name", "Gene name", choices = rownames(HD)
#   ),
#   tabPanel("SeuratCluster", 
#            selectInput("Seurat_cluster", "Seurat")
#   ),
#   tabPanel("BanksyCluster",
#            numericInput("rate", "rate", value = 1, min = 0),
#   ),
#   tabPanel("Celltype",
#            numericInput("rate", "rate", value = 1, min = 0),
#   )
# )



# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Exploration of Visium HD data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          selectInput("HD_disply", "Disply", selected = "Label", choices = c("Transcript", "Label")),
          uiOutput("ui"),
          actionButton("loadData", "Load Data"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
)




