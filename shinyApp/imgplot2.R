raster3d_animation_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, tagList(
        actionButton(ns("zoom_m"), "", icon = icon("minus-square-o"), width = "40px", 
                     style = "border-radius: 25px; padding: 0px;"),
        actionButton(ns("zoom_p"), "", icon = icon("plus-square-o"), width = "40px", 
                     style = "border-radius: 25px; padding: 0px;")
      )),
      column(4, offset = 4, uiOutput(ns("time_slider_UI")))
    ),
    fluidRow(
      lapply(c("x"), function(x) {
        column(6, offset=1, align = "center", uiOutput(ns(paste0(x, "_slider_UI"))))
      })
    ),
    fluidRow(
      lapply(c("x_plot"), function(x) {
        column(6, offset=1, align = "center", plotOutput(ns(x)))
      })
    )
  )
}

raster3d_animation_Module <- function(input, output, session, im) {
  fmax <- length(list.files(im))
  ns <- session$ns
  rv <- reactiveValues(zoom = 1)
  
  output$x_plot <- renderImage({
    req(input$x_slider)
    index_buff <- min(input$x_slider, fmax)
    if (index_buff > fmax) index_buff <- ceiling(fmax/2)
    path <- paste0(im,index_buff,".jpg")
    imgjpg <- load.image(path)
    
    list(
      src = path,
      width = width(imgjpg) * rv$zoom,
      height = height(imgjpg) * rv$zoom
    )
  }, deleteFile = F)
  
  output$x_slider_UI <- renderUI({
    req(im)
    sliderInput(
      ns("x_slider"), label = NULL, min = 1, max = fmax, step = 1,
      value = ceiling(fmax/ 2), 
      animate =T
    )
  })
  
  observeEvent(input$zoom_p, {
    if (rv$zoom < 1.5) rv$zoom <- rv$zoom + 0.05
  })
  observeEvent(input$zoom_m, {
    if (rv$zoom > 0.5) rv$zoom <- rv$zoom - 0.05
  })
 
 
}