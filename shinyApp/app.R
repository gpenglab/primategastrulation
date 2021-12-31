##mfas HE and STGE-plot
#Guizhong Cui
#Bioland Laboratory （Guangzhou Regenerative Medicine and Health Guangdong Laboratory）
# Wed Dec 29 14:04:58 2021 ------------------------------
rm(list=ls())

library(shiny)
library(imager)
library(shinythemes)
library(shinyWidgets)
library(markdown)

source("stge.plot.R")
source("imgplot2.R")


ui <- navbarPage(
    'Geography and molecular anatomy of gastrulating primate',
    theme = shinytheme("flatly"),
    # titlePanel(title="Geography and molecular anatomy of gastrulating primate",windowTitle='Primate Gastrulation'),
    
    #Overview
    tabPanel(
        "Overview",
        includeMarkdown("Overview.md")
    ),
    
    ##serial H.E sections of gastrulating primate embryos
    tabPanel(
        "H.E",
        fluidRow(
            column(3, selectInput(
                "demo_dt", "Choose a Stage", choices = c(
                    "E17" = "data/E17/",
                    "E18" = "data/E18/",
                    "E19" = "data/E19/",
                    "E20" = "data/E20/",
                    "E21" = "data/E21/"
                )
            )),
            
        ),
        
        uiOutput("raster_panel")
    ),
    
    ## spatial gene expression in UMAP and STGE-plot 
    tabPanel(
      "STGE-plot",
      fluidRow(
        column(2, 
               selectizeInput("var", 
                              label = "Type a gene name and press enter to plot",selected=1,
                              choices=NULL)
               ),
        column(4,
               textOutput("selected_var"),
               plotOutput("umap_plot",width = "100%")
               ),
        column(2,
               imageOutput("Fig.umap")
      )
        ),
      fluidRow(
        column(4,
               plotOutput("stge_plot",width = "100%"),
               offset = 2),
        column(2,
               imageOutput("Fig.stge", width = "80%")
        )

      )
      ),
    
    # two gene expression
    tabPanel("Two-Gene STGE-plot", 
             fluidRow(
               column(2, 
                      selectizeInput("var1", 
                                     label = "Type a gene name and press enter to plot",selected=1,
                                     choices=NULL),
                      selectizeInput("var2",
                                     label = "Enter another gene to plot in the same STGE-plot (see Two-Gene expression tab)",selected=1,
                                     choices=NULL)
               ),
               column(4,
                      plotOutput("stge_2gene",width = "100%"),
                      plotOutput("stge_2genecolor",width = "100%", height = "200px")
                      ),
               column(2,
                      imageOutput("Fig.stge2")
               ),
             )
    ),
    
    ##Data download
    tabPanel("Data Download",br(),
             
             p("You can access the preprint ", a(href="https://www.biorxiv.org/", "here",target="https://www.biorxiv.org/"),style="text-align:center;color:white",
             style="text-align:justify;color:white;background-color:rgb(48,62,78);padding:15px;border-radius:15px"),
             br(),
             
             br(),
             p("You can access the Raw data ", a(href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182838", "here",target="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182838"),style="text-align:center;color:white",
             style="text-align:justify;color:white;background-color:rgb(48,62,78);padding:15px;border-radius:15px"),
             br(),
             
             downloadButton("downloadCell", label = "Download Cell counting of germ layers"),
             downloadButton("downloadData", label = "Download gene expression matrix(normalized)"),
             downloadButton("downloadInterspecies", label = "Download Table.S3-the interspecies differences")
             
             )
             
)

server <- function(input, output, session) {

  #H.E query
    output$raster_panel <- renderUI({
        callModule(raster3d_animation_Module, "E17", im = input$demo_dt)
        raster3d_animation_UI("E17")

    })
    
  # gene & signaling query
    updateSelectizeInput(session=session, 'var', choices =  c(Choose = '', gene.names), server = TRUE,selected="TBXT")
    updateSelectizeInput(session=session, 'var1', choices =  c(Choose = '', gene.co), server = TRUE,selected="TBXT")
    updateSelectizeInput(session=session, 'var2', choices =  c(Choose = '', gene.co), server = TRUE,selected="MESP1")
    
    output$umap_plot <- renderPlot({
      validate(
        need(input$var, "Please select another gene")
      )
      df_scatter<-data.table(x=express_vals$UMAP_1,
                             y=express_vals$UMAP_2,
                             z=unlist(express_vals[[input$var]]),
                             shp=express_vals$cluster)
      plot.std.col(df_scatter, xname="", yname="", title=paste0(input$var))
    })
    
    output$Fig.umap <- renderImage({
      filename <- 'image/UMAP.model2.png'
      list(src = filename)
      
    }, deleteFile = FALSE)
    
    output$stge_plot <- renderPlot({
      validate(
        need(input$var, "Please select target gene")
      )
      df_scatter<-data.table(x=express_vals$X,
                             y=express_vals$Y,
                             z=unlist(express_vals[[input$var]]),
                             shp=express_vals$cluster)
      plot.std.col(df_scatter, xname="", yname="", title=paste0(input$var))
    })
    
    output$Fig.stge <- renderImage({
      filename <- 'image/STGE.model2.png'
      list(src = filename)
    }, deleteFile = FALSE)
    
    output$Fig.stge2 <- renderImage({
      filename <- 'image/STGE.model2.png'
      list(src = filename)
    }, deleteFile = FALSE)
    
    output$stge_2gene<-renderPlot({
      plot_2gene(input$var1,input$var2)
    })
    
    output$stge_2genecolor<-renderPlot({
      plot_2genecolor(input$var1,input$var2)
    })
    
    ### Data download
    output$downloadCell <- downloadHandler(
      filename <- "Cell.counting.xlsx",
      
      content <- function(file) {
        file.copy(paste0(base,'Table.S1.xlsx'), file)
      }
    )
    
    output$downloadData <- downloadHandler(
      filename <- "Gastrulating_M.fas_TPM.txt",
      
      content <- function(file) {
        file.copy(paste0(base,'Mfas.TPM.matrix.txt'), file)
      }
    )
    
    output$downloadInterspecies <- downloadHandler(
      filename <- "Interspecies.differences.xlsx",
      
      content <- function(file) {
        file.copy(paste0(base,'Table.S3.xlsx'), file)
      }
    )
    
}

shinyApp(ui, server)

