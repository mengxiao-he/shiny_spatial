#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)

# Define server logic
 function(input, output, session) {
   observeEvent(input$HD_disply, {
     updateTabsetPanel(inputId = "params_label", selected = input$HD_disply)
   })
   
   observeEvent(input$label, {
     updateTabsetPanel(inputId = "params_cluster", selected = input$label)
   })
   
#PLot
output$umap <- renderPlot({
  
if(input$HD_disply=="Label"){
  
  if(input$label=="SeuratCluster"){
    selected_cluster <- input$SeuratCluster
    Idents(mb_hd) <- "seurat_cluster.projected"
  } else if(input$label=="BanksyCluster"){
    selected_cluster <- input$BanksyCluster
    Idents(mb_hd) <- "banksy_cluster"
  } else if(input$label=="CellType"){
    selected_cluster <- input$CellType
    Idents(mb_hd) <- "full_first_type"
  }
  
  cells <- CellsByIdentities(mb_hd, idents = selected_cluster)
  
 DimPlot(mb_hd, reduction = "full.umap.sketch",cells.highlight = cells,
  cols.highlight = c("red", "azure1"), label = F) + 
  ggtitle("Projected clustering (full dataset)") + 
         theme(legend.position = "right")
 
 }
 
 else if(input$HD_disply=="Transcript"){
   selected_genes <- input$gene
   FeaturePlot(mb_hd, features = selected_genes, cols = c("azure1","red"))
   } 
  })


     ## Input spatial variable
output$spatial <- renderPlot({
  
if(input$HD_disply=="Label"){
  
  if(input$label=="SeuratCluster"){
    selected_cluster <- input$SeuratCluster
    Idents(mb_hd) <- "seurat_cluster.projected"
  } else if(input$label=="BanksyCluster"){
    selected_cluster <- input$BanksyCluster
    Idents(mb_hd) <- "banksy_cluster"
  } else if(input$label=="CellType"){
    selected_cluster <- input$CellType
    Idents(mb_hd) <- "full_first_type"
  }
  
  cells <- CellsByIdentities(mb_hd, idents = selected_cluster)
  SpatialDimPlot(mb_hd,
                 cells.highlight = cells,
                 cols.highlight = c("red", "azure1")
  ) + NoLegend()
}
  
  else if(input$HD_disply=="Transcript"){
    selected_genes <- input$gene
    SpatialFeaturePlot(mb_hd, features = selected_genes)+
      scale_fill_gradient2(low = "azure1",high = "red")
  }
})

}
