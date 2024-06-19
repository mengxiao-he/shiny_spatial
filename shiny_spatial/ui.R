#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

parameter_label <- tabsetPanel(
  id = "params_label",
  type = "hidden",

  tabPanel("Label",
           selectInput("label", "Label", choices = c("SeuratCluster", "BanksyCluster", "CellType")),
           ),
  tabPanel("Transcript",
           selectInput("gene", "Gene", choices = rownames(mb_hd))
           )
  )

parameter_cluster <- tabsetPanel(
  id = "params_cluster",
  type = "hidden",
  
  tabPanel("SeuratCluster",
           selectInput("SeuratCluster", "SeuratCluster", choices = levels(mb_hd@meta.data$seurat_cluster.projected)),
  ),
  tabPanel("BanksyCluster",
           selectInput("BanksyCluster", "BanksyCluster", choices = levels(mb_hd@meta.data$banksy_cluster))
  ),
  tabPanel("CellType",
           selectInput("CellType", "CellType", choices = levels(mb_hd@meta.data$full_first_type))
  ),
)




# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Exploration of Visium HD data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          selectInput("HD_disply", "Disply", selected = "Label", 
                      choices = c("Transcript", "Label")),
          
          # uiOutput("ui"),
          parameter_label,
          conditionalPanel(
            condition = "input.HD_disply == 'Label'",
            parameter_cluster
          ),
          #actionButton("loadData", "Load Data"),
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          
          imageOutput("umap",width="600px",height="500px"),
          imageOutput("spatial",width="600px",height="500px"),
        )
    )
)




