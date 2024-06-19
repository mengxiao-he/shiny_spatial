#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# Define server logic
function(input, output, session) {
  
output$ui <- renderUI({
  if(input$type=="Transcript") {
    selectInput("gene", "Gene", choices = rownames(mb_hd))
  } else {
    selectInput("label", "Label", choices = c("SeuratCluster", "BanksyCluster", "Celltype"))
  }
})

output$ui2 <- renderUI({
  if(input$label=="SeuratCluster") {
    selectInput("seurat_cluster", "Cluster", choices = c("All", levels(mb_hd$seurat_cluster.projected)), selected = "All")
  } else if(input$label=="BanksyCluster") {
    selectInput("banksy_cluster", "Cluster", choices = c("All", levels(mb_hd$banksy_cluster)), selected = "All")
  } else if(input$label=="Celltype") {
    selectInput("celltype", "Celltype", choices = levels(mb_hd$full_first_type))
  }
})
    
}
