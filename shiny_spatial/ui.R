#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)


# Define UI for application that draws a histogram
fluidPage(

    # Application title
    titlePanel("Exploration of Visium HD data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
          selectInput("type", "Type", selected = "Label", choices = c("Transcript", "Label")),
          uiOutput("ui"),
          conditionalPanel(
            condition = "input.type == 'Label'",
            uiOutput("ui2")
          ),
          actionButton("loadData", "Load Data"),
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("umapPlot"),
            plotOutput("spatialPlot")
        )
    )
)




