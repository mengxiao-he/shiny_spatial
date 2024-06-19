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
    selectInput("label", "Label", choices = c("SeuratCluster", "BanksyCluster", "Celltype"), selected = "SeuratCluster")
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



plotData <- eventReactive(input$loadPlot, {
  
  if (input$type=="Label") {
    
    if(input$label=="SeuratCluster") {
      group_by_param <- "seurat_cluster.projected"
    } else if(input$label=="BanksyCluster") {
      group_by_param <- "banksy_cluster"
    } else if(input$label=="Celltype") {
      group_by_param <- "full_first_type"
    }
    
    Idents(mb_hd) <- group_by_param
    
    if((input$label == "SeuratCluster" && input$seurat_cluster == "All") || 
       (input$label == "BanksyCluster" && input$banksy_cluster == "All")) {
      
      list(
        umapPlot = Seurat::DimPlot(mb_hd, 
                                   reduction = "full.umap.sketch", 
                                   group.by = group_by_param, 
                                   label = FALSE) + NoAxes(),
        spatialPlot = Seurat::SpatialDimPlot(mb_hd, 
                                             group.by = group_by_param, 
                                             label = FALSE)
      )
      
    } else {
      if(input$label=="SeuratCluster") {
        selected_cluster <- input$seurat_cluster
      } else if(input$label=="BanksyCluster") {
        selected_cluster <- input$banksy_cluster
      } else if(input$label=="Celltype") {
        selected_cluster <- input$celltype
      }
      
      cells_to_plot <- CellsByIdentities(mb_hd, idents = selected_cluster)
      
      list(
        umapPlot = Seurat::DimPlot(mb_hd, 
                                   reduction = "full.umap.sketch", 
                                   label = FALSE, 
                                   cells.highlight = cells_to_plot, 
                                   cols.highlight = c("#FFFF00", "grey50")) + NoLegend() + NoAxes(),
        spatialPlot = Seurat::SpatialDimPlot(mb_hd, 
                                             label = FALSE, 
                                             cells.highlight = cells_to_plot, 
                                             cols.highlight = c("#FFFF00", "grey50")) + NoLegend()
      )
    }
  } else if(input$type=="Transcript") {
    list(
      umapPlot = Seurat::FeaturePlot(mb_hd, 
                                     features = input$gene, 
                                     reduction = "full.umap.sketch", 
                                     cols = c("lightgray", "navyblue")) + NoAxes(),
      spatialPlot = Seurat::SpatialFeaturePlot(mb_hd, 
                                               features = input$gene) 
      & ggplot2::scale_fill_gradient2(high = "navyblue")
    )
  }
  
})

# Render the first plot in the UI
output$umapPlot <- renderPlot({
  plotData()$umapPlot
})

# Render the second plot in the UI
output$spatialPlot <- renderPlot({
  plotData()$spatialPlot
})


}

