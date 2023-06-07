#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#install packages that aren't installed
reqPack <- function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))
    need<-libs[req==FALSE]
    if(length(need)>0){ 
        install.packages(need)
        lapply(need,require,character.only=TRUE)
    }
}

reqPack("shiny", "ggplot2", "reshape2", "ComplexHeatmap", "circlize", "viridis")

#load packages
library(shiny)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

#read in necessary files
norm_expr <- read.delim("normalized_expr_avg.txt", header = T, row.names = 1)
gene_info <- read.delim("NM-symb-name.txt", stringsAsFactors = F)
designInfo <- read.delim("design-info.txt", header = T, row.names = 1)

GO_table <- as.data.frame(read.delim("GO_terms.txt", stringsAsFactors = F))
colnames(GO_table) <- gsub("\\.", ":", colnames(GO_table))
CC_table <- as.data.frame(read.delim("./CC_terms.txt", stringsAsFactors = F))
colnames(CC_table) <- gsub("\\.", ":", colnames(CC_table))
anno <- data.frame(condition = rep(c("P32D", "P32H", "P60D", "P60H"), each=6), row.names = rownames(designInfo))
HManno <- HeatmapAnnotation(apex_base = rep(c("apex", "base", "apex", "base", "apex", "base", "apex", "base"), each = 3),
                df=anno, name = "arrayanno", col = list(condition = c("P32D"="darkorange", "P32H"="darkorchid2",
                "P60D"="darkorange2", "P60H"="darkorchid4"), apex_base = c("apex" = "gray80", "base" = "gray50")),
                annotation_legend_param = list(nrow = 2))

norm_expr.copy <- norm_expr
newColnames <- c("E32DA1", "E32DA2", "E32DA3", "F32DB1", "F32DB2", "F32DB3", "A32HA1", "A32HA2", "A32HA3",
                 "B32HB1", "B32HB3", "B32HB3", "G60DA1", "G60DA2", "G60DA3", "H60DB1", "H60DB2", "H60DB3",
                 "C60HA1", "C60HA2", "C60HA3", "D60HB1", "D60HB2", "D60HB3")
colnames(norm_expr.copy) <- newColnames

GOterms <- read.delim("./GOids_terms.txt", stringsAsFactors = F)
CCterms <- read.delim("./CCids_terms.txt", stringsAsFactors = F)

#define necessary variables
groupingLoc <- c("P32A", "P32B", "P32A", "P32B", "P60A","P60B", "P60A", "P60B")
groupingNoLoc <- c("P32", "P32", "P60", "P60")
hearDeafLoc <- rep(c("Deaf", "Hearing", "Deaf", "Hearing"), each=2)
hearDeafNoLoc <- c("Deaf", "Hearing", "Deaf", "Hearing")


#define functions
SEM <- function(x) {
    sd(x)/sqrt(length(x))
}

getGeneInfo <- function(gene, data, design){
  
  #make sure gene is in dataset
  if (length(which(gene == rownames(data)) == 1)) {
    gene_data <- as.data.frame(t(data[which(gene == rownames(data)), ]))
    colnames(gene_data)[1] <- "value"
    gene_data[,2] <- design$Condition
    colnames(gene_data)[2] <- "Condition"
    gene_data[,3] <- design$age_location
    colnames(gene_data)[3] <- "age_location"
    gene_data[,4] <- design$age_HD
    colnames(gene_data)[4] <- "age_HD"
    gene_data[,5] <- design$age
    colnames(gene_data)[5] <- "age"
    gene_data[,6] <- design$Hear_Deaf
    colnames(gene_data)[6] <- "Hear_Deaf"
    
    #get averages by location
    gene_mean_loc <- aggregate(gene_data[,1], by = list(Condition=gene_data$Condition), FUN = mean)
    colnames(gene_mean_loc)[2] <- gene
    
    #get averages by hearing and deaf
    gene_mean_hd <- aggregate(gene_data[,1], 
                              by = list(Condition=gene_data$age_HD), FUN = mean)
    colnames(gene_mean_hd)[2] <- gene
    
    #get SEM by location
    gene_sem_loc <- aggregate(gene_data[,1], 
                              by = list(Condition=gene_data$Condition), FUN = SEM)
    colnames(gene_sem_loc)[2] <- gene
    
    #get SEM by hearing and deaf
    gene_sem_hd <- aggregate(gene_data[,1], 
                             by = list(Condition=gene_data$age_HD), FUN = SEM)
    colnames(gene_sem_hd)[2] <- gene
    
    #melt the dataframes
    mean_melt_loc <- melt(gene_mean_loc, id.vars = c("Condition"), 
                          variable.name = "gene", value.name = "value")
    mean_melt_loc[,4] <- groupingLoc
    colnames(mean_melt_loc)[4] <- "age_location"
    mean_melt_loc[,5] <- hearDeafLoc
    colnames(mean_melt_loc)[5] <- "Hear_Deaf"
    
    mean_melt_hd <- melt(gene_mean_hd, id.vars = c("Condition"), 
                         variable.name = "gene", value.name = "value")
    mean_melt_hd[,4] <- groupingNoLoc
    colnames(mean_melt_hd)[4] <- "age"
    mean_melt_hd[,5] <- hearDeafNoLoc
    colnames(mean_melt_hd)[5] <- "Hear_Deaf"
    
    sd.meltLoc <- melt(gene_sem_loc, id.vars = c("Condition"), 
                       variable.name = "gene", value.name = "value")
    sd.meltLoc[,4] <- groupingLoc
    colnames(sd.meltLoc)[4] <- "age_location"
    sd.meltHD <- melt(gene_sem_hd, id.vars = c("Condition"), 
                      variable.name = "gene", value.name = "value")
    sd.meltHD[,4] <- groupingNoLoc
    colnames(sd.meltHD)[4] <- "age"
    
    #set factor levels
    gene_data$Condition <- factor(gene_data$Condition, levels = c("P32HA", "P32DA", "P32HB", "P32DB", "P60HA", "P60DA", "P60HB", "P60DB"), ordered = T)
    gene_data$age_HD <- factor(gene_data$age_HD, levels = c("P32hearing", "P32deaf", "P60hearing", "P60deaf"), ordered = T)
    gene_data$Hear_Deaf <- factor(gene_data$Hear_Deaf, levels = c("Hearing", "Deaf"), ordered = T)
    mean_melt_loc$Condition <- factor(mean_melt_loc$Condition, levels = c("P32HA", "P32DA", "P32HB", "P32DB", "P60HA", "P60DA", "P60HB", "P60DB"), ordered = T)
    mean_melt_loc$Hear_Deaf <- factor(mean_melt_loc$Hear_Deaf, levels = c("Hearing", "Deaf"), ordered = T)
    mean_melt_hd$Condition <- factor(mean_melt_hd$Condition, levels = c("P32hearing", "P32deaf", "P60hearing", "P60deaf"), ordered = T)
    mean_melt_hd$Hear_Deaf <- factor(mean_melt_hd$Hear_Deaf, levels = c("Hearing", "Deaf"), ordered = T)
    
    return(list(gene_data, mean_melt_loc, mean_melt_hd, sd.meltLoc, sd.meltHD))
  } else {
    return("Gene not found")
  }
}

getGeneMatrix <- function(GOid, GOanno, exprData) {

    #get list of genes based on provided GO term ID
    geneList <- GOanno[, GOid][!is.na(GOanno[, GOid])]

    #get matrix of expr for genes in geneList
    geneMatrix <- as.matrix(exprData[ geneList, ])
    geneMatrix <- geneMatrix[complete.cases(geneMatrix), ]
    geneMatrix <- log2(geneMatrix)
    geneMatrix <- geneMatrix - rowMeans(geneMatrix)
    return(geneMatrix)
}

blwtrd = colorRamp2(c(-2, 0, 2), c('blue', 'white', 'red')) #color scheme for heatmaps

# Define server logic required to draw plots and generate tables
server <- function(input, output) {
    
    geneInfo <- reactive({
        validate(
            need(input$geneInput != "", "Please enter a gene")
        )
        getGeneInfo(input$geneInput, norm_expr, designInfo)
    })
    
    geneMatrix <- reactive({
        validate(
            need(input$GOtermInput != "", "Please enter a GO term ID")
        )
        if (input$GOcat == "BP"){
            getGeneMatrix(input$GOtermInput, GO_table, norm_expr.copy)
        } else if (input$GOcat == "CC") {
            getGeneMatrix(input$GOtermInput, CC_table, norm_expr.copy)
        }
    })

    term <- reactive({
        validate(
            need(input$GOtermInput != "", "Please enter a GO term ID")
        )
        if (input$GOcat == "BP") {
            GOterms$Term[which(input$GOtermInput == GOterms$Term.Accession)]
        } else if (input$GOcat == "CC") {
            CCterms$Term[which(input$GOtermInput == CCterms$Term.Accession)]
        }
    })
    
    #generate graph based on gene
    output$expressionPlot <- renderPlot({

      if (class(geneInfo()) == "character") {
        return(ggplot() + theme_void() + 
                 labs(title = "Gene not found"))
      } else {
        #get gene expression info
        geneData <- as.data.frame(geneInfo()[1])
        geneMeanLoc <- as.data.frame(geneInfo()[2])
        geneMeanHD <- as.data.frame(geneInfo()[3])
        geneSDLoc <- as.data.frame(geneInfo()[4])
        geneSDHD <- as.data.frame(geneInfo()[5])

        #specify error bars
        limitsLoc <- aes(ymax = geneMeanLoc[, "value"] + geneSDLoc[, "value"], ymin = geneMeanLoc[, "value"] - geneSDLoc[, "value"])
        limitsHD <- aes(ymax = geneMeanHD[, "value"] + geneSDHD[, "value"], ymin = geneMeanHD[, "value"] - geneSDHD[, "value"])
        
        #draw graphs
        if (input$graphType == "column graph"){
            if (input$groupInput == "age and location"){
            expr_plot <- ggplot(geneMeanLoc, aes(age_location, value, fill = Hear_Deaf)) + 
                      geom_bar(position="dodge", stat = "identity") + 
                      geom_errorbar(limitsLoc, position = "dodge") + 
                      labs(y = "Normalized Signal Intensity", title=paste("Gene:", input$geneInput, sep = " ")) + 
                      scale_fill_manual(values = c("darkorchid4", "darkorange2")) +
                      theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
            } else if (input$groupInput == "age only"){
                expr_plot <- ggplot(geneMeanHD, aes(age, value, fill = Hear_Deaf)) + 
                        geom_bar(position="dodge", stat = "identity") + 
                        geom_errorbar(limitsHD, position = "dodge") + 
                        labs(y = "Normalized Signal Intensity", title=paste("Gene:", input$geneInput, sep = " ")) + 
                        scale_fill_manual(values = c("darkorchid4", "darkorange2")) +
                        theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
            }
        } else if (input$graphType == "boxplot"){
            if (input$groupInput == "age and location"){
                expr_plot <- ggplot(geneData, aes(age_location, value, fill=Hear_Deaf)) + geom_boxplot() +
                        labs(y = "Normalized Signal Intensity", title=paste("Gene:", input$geneInput, sep = " ")) + 
                        scale_fill_manual(values = c("darkorchid4", "darkorange2")) +
                        theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
            } else if (input$groupInput == "age only"){
                expr_plot <- ggplot(geneData, aes(age_HD, value, fill=Hear_Deaf)) + geom_boxplot() +
                        labs(y = "Normalized Signal Intensity", title=paste("Gene:", input$geneInput, sep = " ")) + 
                        scale_fill_manual(values = c("darkorchid4", "darkorange2")) +
                        theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))
            }
        }
        return(expr_plot)
      }
    }) #closes output$expresssionPlot
    output$expressionTable1 <- renderTable({
        
        if (is.null(geneInfo())){
            return()
        }
        
        geneData <- as.data.frame(geneInfo()[1])
        
        #generate table
        gene.expr1 <- t(geneData[1])
        cols1 <- colnames(gene.expr1)
        #rownames(gene.expr) <- input$geneInput
        gene.expr1 <- t(gene.expr1[1:12])
        colnames(gene.expr1) <- cols1[1:12]
        gene.expr1
        
    })
    output$expressionTable2 <- renderTable({

        if (is.null(geneInfo())){
            return()
        }
        
        geneData <- as.data.frame(geneInfo()[1])
        
        #generate table
        gene.expr2 <- t(geneData[1])
        cols2 <- colnames(gene.expr2)
        gene.expr2 <- t(gene.expr2[13:24])
        colnames(gene.expr2) <- cols2[13:24]
        gene.expr2
    })
    
    
    output$heatmap <- renderPlot({

      if (class(geneMatrix() == "character")) {
        return(ggplot() + theme_void() +
                 labs(title = "Gene not found"))
      } else {
        HM <- Heatmap(geneMatrix(), column_title = paste(input$GOtermInput, term(), sep = " "), column_title_gp = gpar(fontsize = 20), top_annotation = HManno,
                      col = blwtrd, show_column_dend = FALSE, show_column_names = FALSE, cluster_columns = FALSE,
                      column_order = order(colnames(geneMatrix())),
                      heatmap_legend_param = list(legend_direction = "horizontal", title = "log2 norm expr", border = "gray30"))
        draw(HM, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = T)
      }
    })
    
}

