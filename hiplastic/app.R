#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(Matrix)
library(RColorBrewer)
library(dcon)

library(GenomicRanges)
library(plotgardener)
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(shiny)
library(plotly)

library(shinydashboard)
library(shinyWidgets)

options(scipen=999)



load('shiny_data.RData')
b1 <- summary(b1_chr18) |> mutate(i=i*1e5, j=j*1e5)
l1 <- summary(l1_chr18) |> mutate(i=i*1e5, j=j*1e5)
tot <- summary(total_chr18) |> mutate(i=i*1e5, j=j*1e5)
b1_annot <- fread('/rafalab/lzou/rdsprite/B1_mm10.rmsk.annot.bed',
                  col.names = c('chrom','start','end','element')) |>
  makeGRangesFromDataFrame()
l1_annot <- fread('/rafalab/lzou/rdsprite/L1_mm10.rmsk.annot.bed',
                  col.names = c('chrom','start','end','element')) |>
  makeGRangesFromDataFrame()
windows <- fread('/rafalab/lzou/rdsprite/mm10')

my_heatmap <- function(mat, source='A') {
  (mat |>
     reshape2::melt() |>
     ggplot(aes(x = (Var2*1e5+input$r2[1])/1e6, y = (Var1*1e5+input$r1[1])/1e6)) +
     geom_raster(aes(fill = value,
                     text = paste("x:", (Var2*1e5+input$r2[1])/1e6, "<br>", 
                                  "y:", (Var1*1e5+input$r1[1])/1e6, "<br>", 
                                  "z:", value))) +
     scale_fill_distiller(name='prob', palette='YlOrRd', direction=1) +
     theme_classic() +
     scale_y_reverse(expand=c(0,0)) +
     scale_x_continuous(expand=c(0,0)) +
     xlab('position along chr (mb)') +
     ylab('position along chr (mb)')) |>
    ggplotly(tooltip = 'text', source=source)
}

# Define UI 
ui <- dashboardPage(
  dashboardHeader(title = "fake higlass"),
  dashboardSidebar(
    selectInput('chrom', label = 'Chr:', choices = 'chr18'),
    numericRangeInput('r1', label = 'Region 1', value = c(7e6,13e6)),
    numericRangeInput('r2', label = 'Region 2', value = c(22e6,28e6)),
    actionButton('Deconvolve', 'Deconvolve')
  ),
  dashboardBody(
    fluidRow(
      box(
        title = 'Interactive plot: total DNA-DNA contacts', status = 'primary',
        plotlyOutput('interactiveMat', width='500px', height='400px'),
        sliderInput('maxval', 'Max val:', min=0, max=1000, value=100)
      ),
      box(
        title = 'Blue loop region', status = 'primary',
        selectInput('type', label = 'Select contacts:', 
                    choices = c('Total', 'L1', 'L1/Total', 'B1/Total')),
        plotlyOutput('loopMat', width='500px', height='400px')
      )
    ),
    fluidRow(
      box(
        title = 'Deconvolved L1', status = 'primary',
        plotlyOutput('decon_full_L1', width='500px', height='400px'),
        verbatimTextOutput('selected_L1')
      ),
      box(
        title = 'Deconvolved total', status = 'primary',
        plotlyOutput('decon_full_total', width='500px', height='400px'),
        verbatimTextOutput('selected_total')
      )
    ),
    fluidRow(
      box(
        title = 'Deconvolved L1/total ratio', status = 'primary',
        plotlyOutput('decon_full_ratio', width='500px', height='400px')
      ),
      box(
      )
    )
    # fluidRow(
    #   box(
    #     title = 'Triangles: Both regions', status = 'primary',
    #     plotOutput("triangles", width='450px',height='600px')
    #     ),
    #   box(
    #     title = 'Zoom in on each region', status = 'primary',
    #     plotOutput('triangles_zoom', width='450px', height='600px')
    #   )
    # )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  big_m <- reactive({
    m <- as.matrix(total_chr18)
    m[m>input$maxval] <- input$maxval
    m
  })
  
  loop_coords <- reactive({
    r1 <- input$r1/1e5
    r2 <- input$r2/1e5
    if (r1[1] < r2[1]) {
      c(r1[1], r1[2], r2[1], r2[2])
    } else {
      c(r2[1], r2[2], r1[1], r1[2])
    }
  })
  
  little_m <- reactive({
    if (input$type == 'L1/Total') {
      m <- as.matrix(l1_chr18)/as.matrix(total_chr18)
    } else if (input$type == 'B1/Total') {
      m <- as.matrix(b1_chr18)/as.matrix(total_chr18)
    } else if (input$type=='L1') {
      m <- as.matrix(l1_chr18)
    }
    m[loop_coords()[1]:loop_coords()[2],
      loop_coords()[3]:loop_coords()[4]]
  })
  
  little_L1 <- eventReactive(input$Deconvolve, {
    as.matrix(l1_chr18)[loop_coords()[1]:loop_coords()[2],
                        loop_coords()[3]:loop_coords()[4]]
  })
  
  little_total <- eventReactive(input$Deconvolve, {
    as.matrix(total_chr18)[loop_coords()[1]:loop_coords()[2],
                           loop_coords()[3]:loop_coords()[4]]
  })
  
  decon_yax_L1 <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_L1())
    D <- as.matrix(total_chr18[loop_coords()[3]:loop_coords()[4],
                               loop_coords()[3]:loop_coords()[4]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, DF = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  decon_xax_L1 <- eventReactive(input$Deconvolve, {
    y <- colSums(little_L1())
    D <- as.matrix(total_chr18[loop_coords()[3]:loop_coords()[4],
                               loop_coords()[3]:loop_coords()[4]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, DF = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  output$decon_yax_L1 <- renderPlotly({
    d <- data.frame(x = 1:length(decon_yax_L1()), y = exp(decon_yax_L1()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1 y-axis')
  })
  
  output$decon_xax_L1 <- renderPlotly({
    d <- data.frame(x = 1:length(decon_xax_L1()), y = exp(decon_xax_L1()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1 x-axis')
  })
  
  output$decon_full_L1 <- renderPlotly({
    df1 <- length(decon_xax_L1())
    df2 <- length(decon_yax_L1())
    res <- exp(matrix(data = rep(decon_yax_L1(), df2), nrow=df2, ncol=df2) + 
                 t(matrix(data = rep(decon_xax_L1(), df1), nrow = df1, ncol = df1)))
    my_heatmap(res)
    
  })
  
  decon_yax_total <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_total())
    D <- as.matrix(total_chr18[loop_coords()[3]:loop_coords()[4],
                               loop_coords()[3]:loop_coords()[4]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, DF = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  decon_xax_total <- eventReactive(input$Deconvolve, {
    y <- colSums(little_total())
    D <- as.matrix(total_chr18[loop_coords()[3]:loop_coords()[4],
                               loop_coords()[3]:loop_coords()[4]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, DF = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  decon_yax_total_plot <- reactive({
    d <- data.frame(x = 1:length(decon_yax_total()), y = exp(decon_yax_total()))
    p <- plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved total y-axis')
    p
  })
  
  output$decon_yax_total <- renderPlotly({
    decon_yax_total_plot()
  })
  
  output$decon_xax_total <- renderPlotly({
    d <- data.frame(x = 1:length(decon_xax_total()), y = exp(decon_xax_total()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved total x-axis')
  })
  
  output$decon_full_total <- renderPlotly({
    df1 <- length(decon_xax_total())
    df2 <- length(decon_yax_total())
    res <- exp(matrix(data = rep(decon_yax_total(), df2), nrow=df2, ncol=df2) + 
                 t(matrix(data = rep(decon_xax_total(), df1), nrow = df1, ncol = df1)))
    my_heatmap(res)
    
  })
  
  output$decon_full_ratio <- renderPlotly({
    eps <- 1e6
    df1 <- length(decon_xax_total())
    df2 <- length(decon_yax_total())
    res <- exp(matrix(data = rep(decon_yax_total(), df2), nrow=df2, ncol=df2) + 
                 t(matrix(data = rep(decon_xax_total(), df1), nrow = df1, ncol = df1)))
    df1 <- length(decon_xax_L1())
    df2 <- length(decon_yax_L1())
    res2 <- exp(matrix(data = rep(decon_yax_L1(), df2), nrow=df2, ncol=df2) + 
                  t(matrix(data = rep(decon_xax_L1(), df1), nrow = df1, ncol = df1)))
    my_heatmap(res2-res)
    
  })
  
  output$interactiveMat <- renderPlotly({
    plot_ly(z = big_m(), 
            type = 'heatmap', colorscale = 'YlOrRd', reversescale=T,
            source = 'source') |>
      layout(yaxis = list(autorange='reversed'),
             shapes = list(type = 'rect',
                           fillcolor = NA, line = list(color='blue'),
                           y0 = loop_coords()[1], y1 = loop_coords()[2],
                           yref = 'y', x0 = loop_coords()[3], x1 = loop_coords()[4],
                           xref = 'x'))
  })
  
  output$loopMat <- renderPlotly({
    my_heatmap(little_m())
  })
  
  observeEvent(event_data('plotly_relayout', 'source'), {
    d <- event_data("plotly_relayout", "source")
    if (length(d)==4) {
      updateNumericRangeInput(session, 'r1', value=c(d[[3]]*1e5,d[[4]]*1e5))
      updateNumericRangeInput(session, 'r2', value=c(d[[1]]*1e5,d[[2]]*1e5))
    }
  })
  
  output$selected_L1 <- renderPrint(
    event_data("plotly_click", source = "L1_full")
  )
  
  output$selected_B1 <- renderPrint(
    event_data("plotly_click", source = "B1_full")
  )
  
  # observeEvent(event_data('plotly_click', 'L1_full'), {
  #   d <- event_data("plotly_click", "L1_full")
  #   if (length(d)==4) {
  #     updateNumericRangeInput(session, 'r1', value=c(d[[3]]*1e5,d[[4]]*1e5))
  #     updateNumericRangeInput(session, 'r2', value=c(d[[1]]*1e5,d[[2]]*1e5))
  #   }
  # })
  
  
  #updateNumericRangeInput()
  
  
  output$triangles <- renderPlot({
    minstart <- min(input$r1[1], input$r2[1])
    maxend <- max(input$r1[2], input$r2[2])
    pageCreate(width = 4.5, height = 8.5, default.units = "inches", bg = 'black')
    params <- pgParams(chrom = "chr18", chromstart = minstart, chromend = maxend, 
                       assembly = "mm10",
                       x = 0, width = 4.5, length = 4.5, default.units = "inches")
    plotHicTriangle(data = b1 |> filter(i >= minstart, j <= maxend), params = params,
                    resolution = 100000,
                    zrange = c(0,50),
                    y = 3.5, height = 3.5, just = c("left", "bottom"),
                    palette = colorRampPalette(brewer.pal(n=9, 'Reds')))
    plotHicTriangle(data = l1 |> filter(i >= minstart, j <= maxend), params = params,
                    resolution = 100000,
                    zrange = c(0,50),
                    y = 7, height = 3.5, just = c("left", "bottom"),
                    palette = colorRampPalette(brewer.pal(n=9, 'Blues')))
    plotGenes(
      params = params,
      y = 7, height= 1
    )
    
    plotGenomeLabel(
      params = params,
      y =8, scale = 'Mb'
    )
    
    # tot <- plot_ly(z = as.matrix(total_chr18), type = 'heatmap')
    # l1 <- plot_ly(z = as.matrix(l1_chr18), type = 'heatmap', colors='Blues')
    # b1 <- plot_ly(z = as.matrix(b1_chr18), type='heatmap', colors='Reds')
    # subplot(tot, l1, b1, nrows=3, shareX=T, shareY=T)
  })
  
  
  output$triangles_zoom <- renderPlot({
    pageCreate(width = 4.5, height = 8.5, default.units = "inches", bg = 'black')
    params <- pgParams(chrom = "chr18", chromstart = input$r1[1], chromend = input$r1[2], 
                       assembly = "mm10",
                       x = 0, width = 4.5, length = 4.5, default.units = "inches")
    plotHicTriangle(data = l1 |> filter(i >= input$r1[1], j <= input$r1[2]), params = params,
                    resolution = 100000,
                    zrange = c(0,50),
                    y = 3.5, height = 3.5, just = c("left", "bottom"),
                    palette = colorRampPalette(brewer.pal(n=9, 'Blues')))
    plotGenes(
      params = params,
      y = 3.5, height= 1
    )
    
    plotGenomeLabel(
      params = params,
      y =4.5, scale = 'Mb'
    )
    
    params <- pgParams(chrom = "chr18", chromstart = input$r2[1], chromend = input$r2[2], 
                       assembly = "mm10",
                       x = 0, width = 4.5, length = 4.5, default.units = "inches")
    plotHicTriangle(data = l1 |> filter(i >= input$r2[1], j <= input$r2[2]), params = params,
                    resolution = 100000,
                    zrange = c(0,50),
                    y = 7, height = 3.5, just = c("left", "bottom"),
                    palette = colorRampPalette(brewer.pal(n=9, 'Blues')))
    plotGenes(
      params = params,
      y = 7, height= 1
    )
    
    plotGenomeLabel(
      params = params,
      y =8, scale = 'Mb'
    )
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
