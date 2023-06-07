#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(reshape2)

#when adding heatmap functionality
# library(ComplexHeatmap)
# library(circlize)
# library(viridis)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title and other information
    titlePanel("Transcriptomic analyis of the spiral ganglion of hearing and deafened rats"),
    h3("Green Lab", "University of Iowa", "Department of Biology"),
    h4("Data collected by Erin Bailey"),
    h4("Data analyzed by Erin Bailey, Benjamin M. Gansemer, and Steven H. Green"),
    h5("App written by Benjamin M. Gansemer"),
    tags$a("Rahman, Bailey, Gansemer, et al. Anti-inflammatory Therapy Protects Spiral Ganglion Neurons After Aminoglycoside-Induced Hair Cell Loss. 2023. Neurotherapeutics.",
       href = "https://pubmed-ncbi-nlm-nih-gov.ezp3.lib.umn.edu/36697994/"),
    br(), br(),
    
    #generate panel for interacting with and viewing data
    fluidRow(
        column(4,
                textInput("geneInput", "Enter gene symbol"),
                selectInput("graphType", "Select graph type",
                    choices = c("column graph", "boxplot")),
                selectInput("groupInput", "select groupings",
                    choices = c("age and location", "age only")),
               ),
        column(8,
                plotOutput("expressionPlot"),
                br(),
                tableOutput("expressionTable1"),
                tableOutput("expressionTable2"),
               ),
    ),
    
    
    fluidRow(
        column(4,
               radioButtons("GOcat", "GO category:",
                            c("Biological Process" = "BP",
                              "Cellular Component" = "CC")),
               textInput("GOtermInput", "Enter GO term id"),
               ),
        column(8,
               plotOutput("heatmap"),
               ),
    )
    
    sidebarLayout(
        sidebarPanel(
            textInput("geneInput", "Enter gene symbol"),
            selectInput("graphType", "Select graph type",
                        choices = c("column graph", "boxplot")),
            selectInput("groupInput", "select groupings",
                        choices = c("age and location", "age only")),
            textInput("GOtermInput", "enter GO term ID"),
        ),
        mainPanel(
            plotOutput("expressionPlot"),
            br(), br(),
            tableOutput("expressionTable"),
        ),
)
