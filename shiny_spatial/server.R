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
  if(input$HD_disply=="Transcript") {
    selectInput("gene", "Gene", choices = rownames(mb_hd))
  } else {
    selectInput("label", "Label", choices = c("SeuratCluster", "BanksyCluster", "Celltype"))
  }
  

  
})
  
}
