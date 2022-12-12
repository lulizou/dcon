#
# ADD IN TOTAL DNA-DNA CONTACTS REGARDLESS OF RNA
# for each off-diagonal bin, count total DNA contacts over expected O/E
# quantify total RNA in the bins
# total L1/total RNA in each of the bins


# O/E correlated with  line 4
# is line 5 predictive above and beyond line 4

# 3 ~ 4 + 5 regression - does L1 matter after accounting for RNA

# heavily transcribed regions tend to co-localize

# visualize with and without total RNA normalization

# can check known chromosome speckle locations 2018 cell paper to see if regions



library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(data.table)
library(Matrix)
library(RColorBrewer)
library(dcon)
library(GGally)

library(GenomicRanges)
library(plotgardener)
library("org.Mm.eg.db")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(shiny)
library(plotly)
library(DT)

library(shinydashboard)
library(shinyWidgets)

options(scipen=0)

load('shiny_data_v2.RData')
b1_annot <- fread('/rafalab/lzou/resources/mm10.100kb.windows.B1.coverage',
                  col.names = c('chrom','start','end','nfeatures', 'nbp',
                                'tot', 'frac')) |>
  filter(chrom=='chr18')
l1_annot <- fread('/rafalab/lzou/resources/mm10.100kb.windows.L1.coverage',
                  col.names = c('chrom','start','end','nfeatures', 'nbp',
                                'tot', 'frac')) |>
  filter(chrom=='chr18')
l1md_annot <- fread('/rafalab/lzou/resources/mm10.100kb.windows.L1Md.coverage',
                    col.names = c('chrom','start','end','nfeatures', 'nbp',
                                  'tot', 'frac'))|>
  filter(chrom=='chr18')
l1md_bed <- fread('/rafalab/lzou/resources/L1Md_mm10.rmsk.annot.bed',
                  col.names = c('chrom','start','end','element')) |>
  filter(chrom=='chr18') |>
  arrange(start)
windows <- fread('/rafalab/lzou/resources/mm10.100kb.windows', 
                 col.names = c('chrom','start','end')) |>
  filter(chrom=='chr18')

rna_offset <- fread('/rafalab/lzou/rdsprite/total_RNA_count_windows_100000_clusters_2-100.csv') |>
  filter(chrom=='chr18') |>
  select(-V1)

total_l1_rna <- fread('/rafalab/lzou/rdsprite/count_windows/L1,LINE|LINE,L1_RNA_count_windows_100000_clusters_2-100.csv') |>
  filter(chrom=='chr18')
total_l1md_rna <- fread('/rafalab/lzou/rdsprite/count_windows/L1Md_RNA_count_windows_100000_clusters_2-100.csv') |>
  filter(chrom=='chr18')
total_b1_rna <- fread('/rafalab/lzou/rdsprite/count_windows/Alu,SINE|SINE,Alu_RNA_count_windows_100000_clusters_2-100.csv') |>
  filter(chrom=='chr18')

my_heatmap <- function(mat, input=NULL, source='A', cutoff=0, legend='prob', middle = 0) {
  m <- mat |>
    reshape2::melt() |>
    filter(Var1>cutoff, Var2>cutoff, Var1<(nrow(mat)-cutoff), Var2<(ncol(mat)-cutoff)) |>
    ggplot(aes(x = ((Var2-1)*1e5+input$r2[1])/1e6, y = ((Var1-1)*1e5+input$r1[1])/1e6)) +
    geom_raster(aes(fill = value,
                    text = paste("x:", ((Var2-1)*1e5+input$r2[1])/1e6, "<br>", 
                                 "y:", ((Var1-1)*1e5+input$r1[1])/1e6, "<br>", 
                                 "z:", value))) +
    theme_classic() +
    scale_y_reverse(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    xlab('position along chr (mb)') +
    ylab('position along chr (mb)')
  if (middle==0) {
    m <- m + scale_fill_distiller(name=legend, palette='YlOrRd', direction=1)
  } else {
    m <- m + scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',
                                  midpoint = middle)
  }
  
  m |>
    ggplotly(tooltip = 'text', source=source)
}

# Define UI 
ui <- dashboardPage(
  dashboardHeader(title = "fake higlass"),
  dashboardSidebar(
    selectInput('resolution', label = 'Resolution:', choices = c('100kb', '10kb')),
    selectInput('chrom', label = 'Chr:', choices = 'chr18'),
    numericRangeInput('r1', label = 'Region 1', value = c(9.6e6,16.7e6)),
    numericRangeInput('r2', label = 'Region 2', value = c(53.1e6,59.4e6)),
    actionButton('Deconvolve', 'Deconvolve'),
    actionButton('Save', 'Save region'),
    actionButton('Delete', 'Delete region')
  ),
  dashboardBody(
    DTOutput('saved_regions'),
    fluidRow(
      box(
        title = 'Interactive plot', status = 'primary',
        selectInput('type', label = 'Select contacts:', 
                    choices = c('Total RD', 'Total DD', 'O/E DD', 'L1', 'B1', 'L1/Total',
                                'L1Md', 'L1Md/Total', 'B1/Total', 'Total RNA: L1',
                                'Total RNA: L1Md', 'Total RNA: B1', 'L1/B1'),
                    selected = 'L1Md/Total'),
        plotlyOutput('interactiveMat', width='500px', height='400px'),
        numericRangeInput('vizrange', 'Viz range:', value=c(0,1)),
        sliderInput('mintot', 'Minimum total contacts', min = 0, max = 1000, value=40)
      ),
      box(
        title = 'Blue region', status = 'primary',
        plotlyOutput('loopMat', width='500px', height='500px'),
        plotlyOutput('marginalRNA_yaxis', width='450px', height='120px')
      )
    ),
    # fluidRow(
    #   box(
    #     title = 'Deconvolved L1', status = 'primary',
    #     plotlyOutput('decon_full_L1', width='500px', height='400px'),
    #     verbatimTextOutput('selected_L1')
    #   ),
    #   box(
    #     title = 'Deconvolved total', status = 'primary',
    #     plotlyOutput('decon_full_total', width='500px', height='400px'),
    #     verbatimTextOutput('selected_total')
    #   )
    # ),
    fluidRow(
      box(
        title = 'Deconvolved L1Md with total offset', status = 'primary',
        plotlyOutput('decon_full_L1Md_offset', width='550px', height='500px'),
        plotlyOutput('l1mdfrac_r1', width='450px', height='120px')
      ),
      box(
        title = 'L1Md density', status = 'primary',
        plotlyOutput('l1md_matrix', width='550px', height='500px'),
        plotlyOutput('l1mdfrac_r1_2', width='450px', height='120px')
      )
    ),
    fluidRow(
      box(
        title = 'Deconvolved L1 with total offset', status = 'primary',
        plotlyOutput('decon_full_L1_offset', width='550px', height='500px'),
        plotlyOutput('l1frac_r1', width='450px', height='120px')
      ),
      box(
        title = 'L1 density', status = 'primary',
        plotlyOutput('l1_matrix', width='550px', height='500px'),
        plotlyOutput('l1frac_r1_2', width='450px', height='120px')
      )
    ),
    fluidRow(
      box(
        title = 'Deconvolved L1, no offset', status = 'primary',
        plotlyOutput('decon_full_L1_no_offset', width='550px', height='500px')
      ),
      box(
        title = 'Deconvolved L1 margins, with offset', status = 'primary',
        plotlyOutput('decon_xax_L1', width='500px', height='200px'),
        plotlyOutput('decon_yax_L1', width='500px', height='200px')
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
  
  #################### TABLE OF SAVED REGIONS ####################################
  
  observeEvent(input$Save, {
    data.frame(chrom=input$chrom, s1 = input$r1[1], e1 = input$r1[2],
               s2 = input$r2[1], e2 = input$r2[2]) |>
      write.table('saved_regions.tsv', append = T, quote=F, sep='\t', row.names=F,
                  col.names=F)
  })
  
  observeEvent(input$Delete, {
    fread('saved_regions.tsv') |>
      slice(-input$saved_regions_row_last_clicked) |>
      write.table('saved_regions.tsv', quote=F, sep='\t', row.names=F,
                  col.names=F)
  })
  
  regions <- eventReactive({
    input$Save
    input$Delete
  },{
    fread('saved_regions.tsv', col.names = c('chrom','s1','e1','s2','e2'))
  }, ignoreNULL=F)
  
  output$saved_regions <- renderDT(
    regions(), options = list(lengthChange = FALSE), selection = 'single'
  )
  
  observeEvent(input$saved_regions_row_last_clicked, {
    idx <- input$saved_regions_row_last_clicked
    updateNumericRangeInput(session, 'r1', 
                            value=c(regions()$s1[idx], regions()$e1[idx]))
    updateNumericRangeInput(session, 'r2', 
                            value=c(regions()$s2[idx], regions()$e2[idx]))
  })
  
  # regions <- eventReactive(input$Save, {
  #   fread('saved_regions.tsv')
  # })
  
  ###############################################################################
  
  main_mat <- reactive({
    if (input$type == 'L1/Total') {
      m <- as.matrix(l1_chr18)/as.matrix(total_chr18)
    } else if (input$type == 'B1/Total') {
      m <- as.matrix(b1_chr18)/as.matrix(total_chr18)
    } else if (input$type=='L1') {
      m <- as.matrix(l1_chr18)
    } else if (input$type == 'B1') {
      m <- as.matrix(b1_chr18)
    } else if (input$type=='Total RD') {
      m <- as.matrix(total_chr18)
    } else if (input$type =='Total DD') {
      m <- as.matrix(total_DNA_chr18)
    } else if (input$type == 'O/E DD') {
      m <- as.matrix(OE_mat)
    } else if (input$type == 'L1Md') {
      m <- as.matrix(l1md_chr18)
    } else if (input$type == 'L1Md/Total') {
      m <- as.matrix(l1md_chr18)/as.matrix(total_chr18)
    } else if (input$type == 'Total RNA: L1') {
      m <- matrix(rep(total_l1_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)) +
        t(matrix(rep(total_l1_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)))
    } else if (input$type == 'Total RNA: L1Md') {
      m <- matrix(rep(total_l1md_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)) +
        t(matrix(rep(total_l1md_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)))
    } else if (input$type == 'Total RNA: B1') {
      m <- matrix(rep(total_b1_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)) +
        t(matrix(rep(total_b1_rna$count, nrow(total_l1_rna)), nrow=nrow(total_l1_rna), ncol=nrow(total_l1_rna)))
    } else if (input$type == 'L1/B1') {
      m <- log2((as.matrix(l1_chr18)+1)/(as.matrix(b1_chr18)+1))
    }
    m
  })
  
  big_m <- reactive({
    m <- main_mat()
    m[m<input$vizrange[1]] <- 0
    m[as.matrix(total_chr18)<input$mintot] <- 0
    m[m>input$vizrange[2]] <- input$vizrange[2]
    m <- m + t(m) - diag(diag(m))
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
    m <- main_mat()
    m[loop_coords()[1]:loop_coords()[2],
      loop_coords()[3]:loop_coords()[4]]
  })
  
  
  D_r1 <- reactive({
    d <- as.matrix(total_DNA_chr18)[loop_coords()[1]:loop_coords()[2],
                                    loop_coords()[1]:loop_coords()[2]]
    return(d+t(d)-diag(diag(d)))
  })
  
  D_r2 <- reactive({
    d <- as.matrix(total_DNA_chr18)[loop_coords()[3]:loop_coords()[4],
                                    loop_coords()[3]:loop_coords()[4]]
    return(d+t(d)-diag(diag(d)))
  })
  
  offset_y <- reactive({
    #rowSums(D_r1())
    rna_offset$count[rna_offset$start >= input$r1[1] &
                       rna_offset$start <= input$r1[2]]
  })
  
  offset_x <- reactive({
    #rowSums(D_r2())
    rna_offset$count[rna_offset$start >= input$r2[1] &
                       rna_offset$start <= input$r2[2]]
  })
  
  ########################### L1 ##############################################
  
  little_L1 <- eventReactive(input$Deconvolve, {
    as.matrix(l1_chr18)[loop_coords()[1]:loop_coords()[2],
                        loop_coords()[3]:loop_coords()[4]]
  })
  
  
  decon_yax_L1 <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_L1())
    D <- D_r1()
    D <- dcon:::normalize_hic(D, gamma=0.8)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y), offset=offset_y())
    a <- fit$a
    #(B%*%a)-log(sum(exp(B%*%a)))
    B%*%a
  })
  
  decon_xax_L1 <- eventReactive(input$Deconvolve, {
    y <- colSums(little_L1())
    D <- D_r2()
    D <- dcon:::normalize_hic(D, gamma=0.8)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y), offset=offset_x())
    a <- fit$a
    #(B%*%a)-log(sum(exp(B%*%a)))
    B%*%a
  })
  
  decon_yax_L1_no_offset <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_L1())
    D <- D_r1()
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
    a <- fit$a
    #(B%*%a)-log(sum(exp(B%*%a)))
    B%*%a
  })
  
  decon_xax_L1_no_offset <- eventReactive(input$Deconvolve, {
    y <- colSums(little_L1())
    D <- D_r2()
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
    a <- fit$a
    #(B%*%a)-log(sum(exp(B%*%a)))
    B%*%a
  })
  
  output$decon_yax_L1 <- renderPlotly({
    d <- data.frame(x = 1:length(decon_yax_L1()), y = (decon_yax_L1()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1 y-axis')
  })
  
  output$decon_xax_L1 <- renderPlotly({
    d <- data.frame(x = 1:length(decon_xax_L1()), y = (decon_xax_L1()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1 x-axis')
  })
  
  output$decon_full_L1_offset <- renderPlotly({
    df1 <- length(decon_xax_L1())
    df2 <- length(decon_yax_L1())
    res <- exp(matrix(data = rep(decon_yax_L1(), df1), nrow=df2, ncol=df1) + 
                 t(matrix(data = rep(decon_xax_L1(), df2), nrow = df1, ncol = df2)))
    p1 <- my_heatmap(res, input=input, cutoff=2, legend='rate')
    p2 <- (l1_annot |>
             dplyr::filter(start>= input$r2[1], end <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = frac)) +
             geom_point() +
             geom_line() +
             xlab('position along chr18 (mb)') +
             ylab('fraction of 100kb window with L1') +
             theme_classic()) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
    
  })
  
  
  ########################### L1Md ##############################################
  little_L1Md <- eventReactive(input$Deconvolve, {
    as.matrix(l1md_chr18)[loop_coords()[1]:loop_coords()[2],
                          loop_coords()[3]:loop_coords()[4]]
  })
  
  little_total <- eventReactive(input$Deconvolve, {
    D()
  })
  
  decon_yax_L1Md <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_L1Md())
    D <- D_r1()
    D <- dcon:::normalize_hic(D, gamma=0.8)
    B <- dcon:::construct_basis(1:length(y), df = length(y)-20)
    fit <- dcon:::fit_decon(y, D, df = length(y)-20, offset=offset_y())
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
    #B%*%a
  })
  
  decon_xax_L1Md <- eventReactive(input$Deconvolve, {
    y <- colSums(little_L1Md())
    D <- D_r2()
    D <- dcon:::normalize_hic(D, gamma=0.8)
    B <- dcon:::construct_basis(1:length(y), df = length(y)-20)
    fit <- dcon:::fit_decon(y, D, df = length(y)-20, offset=offset_x())
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
    #B%*%a
  })
  
  decon_yax_L1Md_no_offset <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_L1Md())
    D <- D_r1()
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  decon_xax_L1Md_no_offset <- eventReactive(input$Deconvolve, {
    y <- colSums(little_L1Md())
    D <- D_r2()
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
    a <- fit$a
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  output$decon_yax_L1Md <- renderPlotly({
    d <- data.frame(x = 1:length(decon_yax_L1Md()), y = (decon_yax_L1Md()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1Md y-axis')
  })
  
  output$decon_xax_L1Md <- renderPlotly({
    d <- data.frame(x = 1:length(decon_xax_L1Md()), y = (decon_xax_L1Md()))
    plot_ly(d, x = ~x, y = ~y, type = 'scatter', mode = 'lines') |>
      layout(title = 'Deconvolved L1Md x-axis')
  })
  
  output$decon_full_L1Md_offset <- renderPlotly({
    df1 <- length(decon_xax_L1Md())
    df2 <- length(decon_yax_L1Md())
    res <- exp(matrix(data = rep(decon_yax_L1Md(), df1), nrow=df2, ncol=df1) + 
                 t(matrix(data = rep(decon_xax_L1Md(), df2), nrow = df1, ncol = df2)))
    p1 <- my_heatmap(res, input=input, cutoff=2, legend='rate')
    p2 <- (l1md_annot |>
             dplyr::filter(start>= input$r2[1], end <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = frac)) +
             geom_point() +
             geom_line() +
             xlab('position along chr18 (mb)') +
             ylab('fraction of 100kb window with L1Md') +
             theme_classic()) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
    
  })
  
  
  ##############################################################################
  
  output$decon_full_L1_no_offset <- renderPlotly({
    df1 <- length(decon_xax_L1_no_offset())
    df2 <- length(decon_yax_L1_no_offset())
    res <- exp(matrix(data = rep(decon_yax_L1_no_offset(), df1), nrow=df2, ncol=df1) + 
                 t(matrix(data = rep(decon_xax_L1_no_offset(), df2), nrow = df1, ncol = df2)))
    message(sum(res))
    p1 <- my_heatmap(res, input=input, cutoff=2, legend='prob')
    p2 <- (l1_annot |>
             dplyr::filter(start>= input$r2[1], end <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = frac)) +
             geom_point() +
             geom_line() +
             xlab('position along chr18 (mb)') +
             ylab('fraction of 100kb window with L1') +
             theme_classic()) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
    
  })
  
  decon_yax_total <- eventReactive(input$Deconvolve, {
    y <- rowSums(little_total())
    D <- as.matrix(total_chr18[loop_coords()[1]:loop_coords()[2],
                               loop_coords()[1]:loop_coords()[2]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
    a <- fit$a
    
    (B%*%a)-log(sum(exp(B%*%a)))
  })
  
  decon_xax_total <- eventReactive(input$Deconvolve, {
    y <- colSums(little_total())
    D <- as.matrix(total_chr18[loop_coords()[3]:loop_coords()[4],
                               loop_coords()[3]:loop_coords()[4]])
    D <- dcon:::normalize_hic(D, gamma=0.5)
    B <- dcon:::construct_basis(1:length(y), df = length(y))
    fit <- dcon:::fit_decon(y, D, df = length(y))
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
    res <- exp(matrix(data = rep(decon_yax_total(), df1), nrow=df2, ncol=df1) + 
                 t(matrix(data = rep(decon_xax_total(), df2), nrow = df1, ncol = df2)))
    my_heatmap(res, input=input)
    
  })
  
  
  output$interactiveMat <- renderPlotly({
    if (input$type == 'O/E DD') {
      plot_ly(z = big_m(),
              type = 'heatmap', colorscale = 'RdBu', zmid=1, 
              source = 'source') |>
        layout(yaxis = list(autorange='reversed'),
               shapes = list(type = 'rect',
                             fillcolor = NA, line = list(color='blue'),
                             y0 = loop_coords()[1], y1 = loop_coords()[2],
                             yref = 'y', x0 = loop_coords()[3], x1 = loop_coords()[4],
                             xref = 'x'))
    } else {
      plot_ly(z = big_m(), 
              type = 'heatmap', colorscale = 'YlOrRd', reversescale=T,
              source = 'source') |>
        layout(yaxis = list(autorange='reversed'),
               shapes = list(type = 'rect',
                             fillcolor = NA, line = list(color='blue'),
                             y0 = loop_coords()[1], y1 = loop_coords()[2],
                             yref = 'y', x0 = loop_coords()[3], x1 = loop_coords()[4],
                             xref = 'x'))
    }
  })
  
  output$loopMat <- renderPlotly({
    if (input$type=='O/E DD') {
      p1 <- my_heatmap(little_m(), input=input, middle=1)
    } else {
      p1 <- my_heatmap(little_m(), input=input)
    }
    p2 <- (rna_offset |>
             filter(start >= input$r2[1], start <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = count)) +
             geom_point() +
             geom_line() +
             theme_classic() +
             xlab('position along chr (mb)') +
             ylab('Total RNA')) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
  })
  
  output$marginalRNA_yaxis <- renderPlotly({
    (rna_offset |>
       filter(start >= input$r1[1], start <= input$r1[2]) |>
       ggplot(aes(x = start/1e6, y = count)) +
       geom_point() +
       geom_line() +
       theme_classic() +
       xlab('position along chr (mb)') +
       ylab('Total RNA')) |>
      ggplotly()
  })
  
  observeEvent(event_data('plotly_relayout', 'source'), {
    d <- event_data("plotly_relayout", "source")
    if (length(d)==4) {
      updateNumericRangeInput(session, 'r1', 
                              value=c(round(d[[3]])*1e5,round(d[[4]])*1e5))
      updateNumericRangeInput(session, 'r2', 
                              value=c(round(d[[1]])*1e5,round(d[[2]])*1e5))
    }
  })
  
  observeEvent(input$type, {
    if (input$type=='O/E DD') {
      updateNumericRangeInput(session, 'vizrange', value=c(0,3))
    } else if (input$type=='L1/Total') {
      updateNumericRangeInput(session, 'vizrange', value=c(0,1))
    }
  })
  
  output$selected_L1 <- renderPrint(
    event_data("plotly_click", source = "L1_full")
  )
  
  output$selected_B1 <- renderPrint(
    event_data("plotly_click", source = "B1_full")
  )
  
  output$l1frac_r1 <- renderPlotly({
    (l1_annot |>
       dplyr::filter(start>= input$r1[1], end <= input$r1[2]) |>
       ggplot(aes(x = start/1e6, y = frac)) +
       geom_point() +
       geom_line() +
       xlab('position along chr18 (mb)') +
       ylab('L1 frac, y-axis') +
       theme_classic()) |>
      ggplotly()
  })
  
  output$l1mdfrac_r1 <- renderPlotly({
    (l1md_annot |>
       dplyr::filter(start>= input$r1[1], end <= input$r1[2]) |>
       ggplot(aes(x = start/1e6, y = frac)) +
       geom_point() +
       geom_line() +
       xlab('position along chr18 (mb)') +
       ylab('L1Md frac, y-axis') +
       theme_classic()) |>
      ggplotly()
  })
  
  output$l1mdfrac_r1_2 <- renderPlotly({
    (l1md_annot |>
       dplyr::filter(start>= input$r1[1], end <= input$r1[2]) |>
       ggplot(aes(x = start/1e6, y = frac)) +
       geom_point() +
       geom_line() +
       xlab('position along chr18 (mb)') +
       ylab('L1Md frac, y-axis') +
       theme_classic()) |>
      ggplotly()
  })
  
  output$l1md_matrix <- renderPlotly({
    l1 <- l1md_annot |> 
      filter(start >= input$r1[1], start <= input$r1[2]) |>
      pull(frac)
    l2 <- l1md_annot |>
      filter(start >= input$r2[1], start <= input$r2[2]) |>
      pull(frac)
    m <- matrix(rep(l1, length(l2)), nrow=length(l1), ncol=length(l2)) +
      t(matrix(rep(l2, length(l1)), nrow = length(l2), ncol=length(l1)))
    p1 <- my_heatmap(m, input=input, legend='L1Md coverage')
    p2 <- (l1md_annot |>
             dplyr::filter(start>= input$r2[1], end <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = frac)) +
             geom_point() +
             geom_line() +
             xlab('position along chr18 (mb)') +
             ylab('fraction of 100kb window with L1Md') +
             theme_classic()) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
  })
  
  output$l1frac_r1_2 <- renderPlotly({
    (l1_annot |>
       dplyr::filter(start>= input$r1[1], end <= input$r1[2]) |>
       ggplot(aes(x = start/1e6, y = frac)) +
       geom_point() +
       geom_line() +
       xlab('position along chr18 (mb)') +
       ylab('L1 frac, y-axis') +
       theme_classic()) |>
      ggplotly()
  })
  
  output$l1_matrix <- renderPlotly({
    l1 <- l1_annot |> 
      filter(start >= input$r1[1], start <= input$r1[2]) |>
      pull(frac)
    l2 <- l1_annot |>
      filter(start >= input$r2[1], start <= input$r2[2]) |>
      pull(frac)
    m <- matrix(rep(l1, length(l2)), nrow=length(l1), ncol=length(l2)) +
      t(matrix(rep(l2, length(l1)), nrow = length(l2), ncol=length(l1)))
    p1 <- my_heatmap(m, input=input, legend='L1 coverage')
    p2 <- (l1_annot |>
             dplyr::filter(start>= input$r2[1], end <= input$r2[2]) |>
             ggplot(aes(x = start/1e6, y = frac)) +
             geom_point() +
             geom_line() +
             xlab('position along chr18 (mb)') +
             ylab('fraction of 100kb window with L1') +
             theme_classic()) |>
      ggplotly()
    subplot(p1, p2, nrows=2, shareX=T, heights=c(0.8,0.2))
  })
  
  
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
