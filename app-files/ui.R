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
    h3("Green Lab"),
    h3("University of Iowa - Dept. of Biology"),
    h5("Data collected by Erin Bailey"),
    h5("Data analyzed by Erin Bailey, Benjamin M. Gansemer, and Steven H. Green"),
    h6("App written by Benjamin M. Gansemer"),
    tags$a("Rahman, Bailey, Gansemer, et al. Anti-inflammatory Therapy Protects Spiral Ganglion Neurons After Aminoglycoside-Induced Hair Cell Loss. 2023. Neurotherapeutics.",
       href = "https://pubmed-ncbi-nlm-nih-gov.ezp3.lib.umn.edu/36697994/"),
    br(), br(),
    
    #generate panel for interacting with and viewing data
    fluidRow(
      sidebarPanel(
                textInput("geneInput", "Enter gene symbol"),
                selectInput("graphType", "Select graph type",
                    choices = c("column graph", "boxplot")),
                selectInput("groupInput", "select groupings",
                    choices = c("age and location", "age only")),
      ),
        mainPanel(
                plotOutput("expressionPlot"),
                br(),
                tableOutput("expressionTable1"),
                tableOutput("expressionTable2"),
               ),
    ),

    fluidRow(
        sidebarPanel(
               radioButtons("gene_list", "Gene list type:",
                            c("GO Biological Process" = "BP",
                              "GO Cellular Component" = "CC",
                              "GO Molecular Function" = "MF")),
               textInput("GOtermInput", "Enter GO term id"),
               ),
        mainPanel(
               plotOutput("heatmap"),
               ),
    )
)
