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
   observeEvent(input$HD_disply, {
     updateTabsetPanel(inputId = "params_label", selected = input$HD_disply)
   })
   
   observeEvent(input$label, {
     updateTabsetPanel(inputId = "params_cluster", selected = input$label)
   })
   
     
     output$plot <- renderPlot({
       ggplot(iris, aes(Sepal.Length, Petal.Length, col = Species)) +
         geom_point() +
         facet_wrap(~Species) +
         labs(title = "Iris dataset", subtitle = "Scatterplots of Sepal and Petal lengths", caption = "The iris dataset by Edgar Anderson") +
         theme_grey(base_size = 16)
     })
     

#   
}
