library(shiny)
library(hdf5r)
library(deepG)
library(plotly)
library(reshape2)
library(rtracklayer)
library(DT)
library(matrixStats)
library(shinybusy)

# Make a list of genomes where we have metadata, CDS and states
cds_data <- list.files("files/tsv_renamed/")
cds_data_filename <-
  substr(cds_data, 1, nchar(cds_data) - 4) # get rid of file ending
annot <- readRDS("files/bacteria_viral_1_1_metadata.Rds")
annot$filename2 <-
  substr(annot$filename, 1, nchar(annot$filename) - 6) # get rid of file ending
state_data <- list.files("files/density_500_1")
state_data_filename <-
  substr(state_data, 1, nchar(state_data) - 9) # get rid of file ending
state_data_filename %in% cds_data_filename
sample_list <-
  state_data_filename[which(state_data_filename %in% cds_data_filename)]
gff_data <- list.files("before_2008_gff")
gff_data_filename <-
  substr(gff_data, 1, nchar(state_data) - 4) # get rid of file ending
sample_list <-
  sample_list[which(sample_list %in% gff_data_filename)]

# cosine distance function for loci to genome
distFun <- function(genome, sequence) {
  genome <- t(genome) / sqrt(colSums(genome ^ 2))
  sequence <- t(sequence) / sqrt(colSums(sequence ^ 2))
  slen <- nrow(sequence)
  glen <- nrow(genome)
  vapply(seq_len(glen - slen + 1), function(offset) {
    sum(genome[seq.int(offset, length.out = slen), , drop = FALSE] * sequence) / slen
  }, FUN.VALUE = numeric(1))
}


distFunL <- function(genome, sequence, lval = 1) {
  genome <- t(genome)
  sequence <- t(sequence)
  slen <- nrow(sequence)
  glen <- nrow(genome)
  vapply(seq_len(glen - slen + 1), function(offset) {
    mean(rowSums(abs(genome[seq.int(offset, length.out = slen), , drop = FALSE] - sequence) ^
                   lval) ^ (1 / lval))
  }, FUN.VALUE = numeric(1))
}


ui <- navbarPage("GenomeNet Search",
                 tabPanel("Precalculated DB",
                          sidebarLayout(
                            sidebarPanel(
                              width = 5,
                              helpText("First step: select a genome"),
                              selectInput("genome", "Choose a genome:", unique(sample_list)),
                              helpText(
                                "Second step: select a search query by clicking on a genetic element in the table"
                              ),
                              dataTableOutput('gfftable'),
                              helpText(
                                "Third step: select a method to quantify similarity between query and database"
                              ),
                              selectInput("method", "Method", c("Cosine similarity", "L1"))
                            ),
                            # Show a plot of the generated distribution
                            mainPanel(
                              width = 7,
                              plotlyOutput("variancePlot"),
                              #   plotlyOutput("statesPlot"),
                              uiOutput("searchResultsHelp"),
                              
                            
                              dataTableOutput("selectedGffEntryTable"),
                              #      dataTableOutput("selectedStatesSubsetTable"),
                              uiOutput("selectedStatesSubsetPlotHelp"),
                              
                              plotlyOutput("selectedStatesSubsetPlot"),
                              #    dataTableOutput("searchResultDataTable"),
                              uiOutput("searchResultDataAggTableHelp"),
                              dataTableOutput("searchResultDataAggTable"),
                              add_busy_bar(color = "#000000")
                            )
                          )
                          ),
    
                 tabPanel("Upload own data",
                          sidebarPanel(
                            width = 5,

                            helpText(
                              "Step 1: Please obtain states by follwoing the instructions here https://colab.research.google.com/drive/1DW01yDom_hsOHJT-IA3zwIeK3psxc_E5?usp=sharing" ),
                            helpText(
                              "Step 2: Please upload the states file you get from Step 1" ),
                            fileInput("file_states", "Choose hdf5 File",
                                      multiple = FALSE,
                                      accept = c(".h5", ".hdf5")),
                            helpText("Step 3: Please upload a .gff file e.g. you can obtain from running Prokka on the input FASTA sequence"),
                            fileInput("file_gff", "Choose gff File",
                                      multiple = FALSE,
                                      accept = c(".gff")),
                            dataTableOutput('gfftable2'),
                            selectInput("uploaded_method", "Method", c("Cosine similarity", "L1"))
                          ),   mainPanel(
                            width = 7,
                            plotlyOutput("variancePlotUploaded"),
#                            uiOutput("uploadedSearchResultsHelp"),
                            dataTableOutput("uploadedSelectedGffEntryTable"),
                            plotlyOutput("uploadedSelectedStatesSubsetPlot"),
                            dataTableOutput("uploadedSearchResultDataAggTable"),
                            add_busy_bar(color = "#000000")
                            
                          )
                 ),
                 tabPanel("Comparison",
                          sidebarPanel(
                            width = 5,
                            helpText("First step: select a genome"),
                            selectInput("genome", "Choose a genome:", unique(sample_list))),
                          mainPanel(
                            width = 7,
                            helpText("Results go here"),
                            plotlyOutput("densityPlot")
                          )
                          
                 ),
                 )

################################################################################
server <- function(input, output, session) {
  options(shiny.maxRequestSize=100*1024^2)
  
  # results of search
  res <- reactiveValues()
  
  gffData <- reactive({
    file <- paste0("before_2008_gff/", input$genome, ".gff")
    gff <- readGFF(file)
    gff
  })
  
  statesData <- reactive({
    file <- paste0("files/density_500_1/", input$genome, ".fasta.h5")
    states <-
      deepG::readRowsFromH5(
        h5_path = file,
        complete = TRUE,
        getTargetPositions = TRUE
      )
    df <- as.data.frame(states)
    df
  })
  
  annotData <- reactive({
    file <- paste0("files/tsv_renamed/", input$genome, ".tsv")
    annot <- read.csv2(file, sep = "\t", stringsAsFactors = F)
    annot
  })
  
  
  ################################################################################
  # Show the full GFF table
  output$gfftable <-
    DT::renderDataTable(
      selection = 'single',
      options = list(searchHighlight = TRUE),
      filter = 'top',
      {
        dat <- gffData()
        df <- as.data.frame(dat)
        df$seqid <- NULL
        df$strand <- NULL
        df$score <- NULL
        df$source <- NULL
        df$phase <- NULL
        df$db_xref <- NULL
        df$inference <- NULL
        df$ID <- NULL
        df$Name <- NULL
        df$locus_tag <- NULL
        df$eC_number <- NULL
        df$note <- NULL
        df$rpt_family <- NULL
        df$rpt_type <- NULL
        df$rpt_unit_seq <- NULL
        df$size <- df$end - df$start
        df
      }
    )
  
  # Get index of the subset selected by the gff table
  selectedGffEntry <- reactive({
    req(input$gfftable_rows_selected)
    gff <- as.data.frame(gffData())
    selected_entry <- gff[input$gfftable_rows_selected, ]
    selected_entry
  })
  
  output$selectedGffEntryTable <-
    DT::renderDataTable(selection = "none", {
      selectedGffEntry()
    })
  
  selectedStatesSubset <- reactive({
    req(selectedGffEntry())
    df <- statesData()
    df <- as.data.frame(df)
    targetPos <- df$targetPos
    positions <-
      which(targetPos > selectedGffEntry()$start &
              targetPos < selectedGffEntry()$end)
    df[positions, ]
  })
  
  output$selectedStatesSubsetTable <-
    DT::renderDataTable(selection = "none", {
      selectedStatesSubset()
    })
  
  output$selectedStatesSubsetPlot <- renderPlotly({
    mat <- selectedStatesSubset()
    mat$targetPos <- NULL
    plot_ly(z = data.matrix(mat),
            colors = "Greys",
            type = "heatmap")
  })
  
  searchResultData <- reactive({
    req(input$method)
    genome <- statesData()
    genome_targetPos <- genome$targetPos
    genome$targetPos <- NULL
    genome_mat <- data.matrix(t(genome))
    sequence <- selectedStatesSubset()
    sequence$targetPos <- NULL
    sequence_mat <- data.matrix(t(sequence))
    as.data.frame(sequence_mat)
    if (input$method == "L1") {
      distance <- distFunL(genome_mat, sequence_mat)
    } else {
      distance <- distFun(genome_mat, sequence_mat)
    }
    # as.data.frame(dat)
    
    length(distance) <- nrow(genome)
    
    #  full_distance <- tail(c(NA * 1:nrow(full_states), full_distance), nrow(full_states))
    
    df <-
      data.frame(
        pos = 1:length(distance),
        similarity = distance,
        genome_pos = genome_targetPos
      )
    if (input$method == "L1") {
      df <- df[order(df$similarity, decreasing = F), ]
    } else {
      df <- df[order(df$similarity, decreasing = T), ]
    }
    df
  })
  
  output$searchResultDataTable <-
    DT::renderDataTable(selection = "none", {
      dat <- searchResultData()
      as.data.frame(dat)
    })
  
  searchResultDataAgg <- reactive({
    gff_data <- as.data.frame(gffData())
    distances <- searchResultData()
    # annotate overlapping distances with gffs
    distances$start <- distances$genome_pos
    distances$end <- distances$genome_pos
    # annotate overlap with a ORF
    dt.gff <-
      data.table(
        start = gff_data$start,
        end = gff_data$end,
        feature = gff_data$product
      )
    distances <- as.data.table(distances)
    setkey(dt.gff, start, end)
    annotated <-
      foverlaps(distances, dt.gff, type = "within", mult = "first")
    annotated <- annotated[which(!is.na(annotated$feature)), ]
    #  annotated$feature <- ifelse(is.na(annotated$feature), "No annotation", annotated$feature)
    annotated$genome_pos <- NULL
    annotated$i.start <- NULL
    annotated$i.end <- NULL
    annotated$pos <- NULL
    annotated$similarity <- round(annotated$similarity, digits = 3)
    annotated
  })
  
  output$searchResultDataAggTable <-
    DT::renderDataTable(selection = "none", {
      dat <- searchResultDataAgg()
      
      as.data.frame(dat)
    })
  
  
  
  ################################################################################
  # states plot
  output$statesPlot <- renderPlotly({
    req(statesData())
    dm <- reshape2::melt(statesData(), id.vars = "targetPos")
    p <- ggplot(dm, aes(targetPos, value, group = variable)) #+
    p <- p + theme_classic()
    p <- p + geom_line(alpha = .1, size = .1)
    p <- p + ggtitle(input$genome)
    ggplotly(p, dynamicTicks = FALSE, tooltip = NULL)
  })
  
  # variance plot
  output$variancePlot <- renderPlotly({
    req(statesData())
    
    # prepare gff data
    gff <- gffData()
    gff <- as.data.frame(gff)
    states <- statesData()
    states_pos <- states$targetPos
    states$targetPos <- NULL
    df <-
      data.frame(variance = rowVars(data.matrix(states)), pos = states_pos)
    
    
    plot_ly(df,
            x = ~ pos,
            y = ~ variance,
            source = "variance") %>%
      add_lines() %>%
      layout(title = 'Variability in the genomic representations')  %>%
      rangeslider(start = as.numeric(gff[input$gfftable_rows_selected,]$start),
                  end = as.numeric(gff[input$gfftable_rows_selected,]$end))
  })
  
  
  
  ############# UI
  
  output$searchResultsHelp <- renderUI({
    req(searchResultDataAgg())
    helpText("A search is performed with the selected query shown below:")
  })
  
  output$selectedStatesSubsetPlotHelp <- renderUI({
    req(searchResultDataAgg())
    helpText("The representations for the selected loci visualized as a heatmap")
  })
  
  
  output$searchResultDataAggTableHelp <- renderUI({
    req(searchResultDataAgg())
    helpText(
      "Lcoi specified in the genome annoation are ranked based on the similarity to the query:"
    )
  })
 
  
  
  ##### Density
  
  # variance plot
  output$densityPlot <- renderPlotly({
    req(statesData())
    req(input$genome)
    # prepare gff data
    gff <- gffData()
    gff <- as.data.frame(gff)
    states <- statesData()
    states_pos <- states$targetPos
    states$targetPos <- NULL
    dm <- melt(states)
    p <- ggplot(dm, aes(x = value)) + geom_histogram()
 
    ggplotly(p)
  })
  
  
  ##### File uploaded
  
  
  uploadedStatesData <- reactive({
    # input$file will be NULL initially
    req(input$file_states)
    states <-
      deepG::readRowsFromH5(
        h5_path = input$file_states$datapath,
        complete = TRUE,
        getTargetPositions = TRUE)
    df <- as.data.frame(states)
      df
    })
    
 
  uploadedGffData <- reactive({
   req(input$file_gff)
   gff <- readGFF(input$file_gff$datapath)
   gff
 })
 
 
 # Show the full GFF table
  output$gfftable2 <-  DT::renderDataTable(
     selection = 'single',
     options = list(searchHighlight = TRUE),
     filter = 'top',
     {
       req(uploadedGffData())
       dat <- uploadedGffData()
       df <- as.data.frame(dat)
       df$seqid <- NULL
       df$strand <- NULL
       df$score <- NULL
       df$source <- NULL
       df$phase <- NULL
       df$db_xref <- NULL
       df$inference <- NULL
       df$ID <- NULL
       df$Name <- NULL
       df$locus_tag <- NULL
       df$eC_number <- NULL
       df$note <- NULL
       df$rpt_family <- NULL
       df$rpt_type <- NULL
       df$rpt_unit_seq <- NULL
       df$size <- df$end - df$start
       df
     }
   )
  
  # variance plot
  output$variancePlotUploaded <- renderPlotly({
    req(uploadedStatesData())
    req(uploadedGffData())
    # prepare gff data
    gff <- uploadedGffData()
    gff <- as.data.frame(gff)
    states <- uploadedStatesData()
    states_pos <- states$targetPos
    states$targetPos <- NULL
    df <-
      data.frame(variance = rowVars(data.matrix(states)), pos = states_pos)
    
    plot_ly(df,
            x = ~ pos,
            y = ~ variance,
            source = "uploadedVariance") %>%
      add_lines() %>%
      layout(title = 'Variability in the genomic representations')  %>%
      rangeslider(start = as.numeric(gff[input$gfftableuploaded_rows_selected,]$start),
                  end = as.numeric(gff[input$gfftableuploaded_rows_selected,]$end))
  })
  
  # Get index of the subset selected by the gff table
  uploadedSelectedGffEntry <- reactive({
    req(uploadedStatesData())
    req(uploadedGffData())
    req(input$gfftable2_rows_selected)
    gff <- as.data.frame(uploadedGffData())
    selected_entry <- gff[input$gfftable2_rows_selected, ]
    selected_entry
  })
  
  output$uploadedSelectedGffEntryTable <-
    DT::renderDataTable(selection = "none", {
      req(uploadedStatesData())
      req(uploadedGffData())
      uploadedSelectedGffEntry()
    })
  
  
  
  uploadedSelectedStatesSubset <- reactive({
    req(uploadedStatesData())
    req(uploadedGffData())
    req(uploadedSelectedGffEntry())
    df <- uploadedStatesData()
    df <- as.data.frame(df)
    targetPos <- df$targetPos
    positions <-
      which(targetPos > uploadedSelectedGffEntry()$start &
              targetPos < uploadedSelectedGffEntry()$end)
    df[positions, ]
  })
  
  output$uploadedSelectedStatesSubsetTable <-
    DT::renderDataTable(selection = "none", {
      req(uploadedStatesData())
      req(uploadedGffData())
      uploadedSelectedStatesSubset()
    })
  
  output$uploadedSelectedStatesSubsetPlot <- renderPlotly({
    req(uploadedStatesData())
    req(uploadedGffData())
    mat <- uploadedSelectedStatesSubset()
    mat$targetPos <- NULL
    plot_ly(z = data.matrix(mat),
            colors = "Greys",
            type = "heatmap")
  })
  
  uploadedSearchResultData <- reactive({
    req(uploadedStatesData())
    req(uploadedGffData())
    req(input$uploaded_method)
    genome <- uploadedStatesData()
    genome_targetPos <- genome$targetPos
    genome$targetPos <- NULL
    genome_mat <- data.matrix(t(genome))
    sequence <- uploadedSelectedStatesSubset()
    sequence$targetPos <- NULL
    sequence_mat <- data.matrix(t(sequence))
    as.data.frame(sequence_mat)
    if (input$uploaded_method == "L1") {
      distance <- distFunL(genome_mat, sequence_mat)
    } else {
      distance <- distFun(genome_mat, sequence_mat)
    }
    length(distance) <- nrow(genome)
    df <-
      data.frame(
        pos = 1:length(distance),
        similarity = distance,
        genome_pos = genome_targetPos
      )
    if (input$uploaded_method == "L1") {
      df <- df[order(df$similarity, decreasing = F), ]
    } else {
      df <- df[order(df$similarity, decreasing = T), ]
    }
    df
  })
  
  output$uploadedSearchResultDataTable <-
    DT::renderDataTable(selection = "none", {
      req(uploadedStatesData())
      req(uploadedGffData())
      dat <- uploadedSearchResultData()
      as.data.frame(dat)
    })
  
  uploadedSearchResultDataAgg <- reactive({
    req(uploadedStatesData())
    req(uploadedGffData())
    gff_data <- as.data.frame(uploadedGffData())
    distances <- uploadedSearchResultData()
    # annotate overlapping distances with gffs
    distances$start <- distances$genome_pos
    distances$end <- distances$genome_pos
    # annotate overlap with a ORF
    dt.gff <-
      data.table(
        start = gff_data$start,
        end = gff_data$end,
        feature = gff_data$product
      )
    distances <- as.data.table(distances)
    setkey(dt.gff, start, end)
    annotated <-
      foverlaps(distances, dt.gff, type = "within", mult = "first")
    annotated <- annotated[which(!is.na(annotated$feature)), ]
    #  annotated$feature <- ifelse(is.na(annotated$feature), "No annotation", annotated$feature)
    annotated$genome_pos <- NULL
    annotated$i.start <- NULL
    annotated$i.end <- NULL
    annotated$pos <- NULL
    annotated$similarity <- round(annotated$similarity, digits = 3)
    annotated
  })
  
  output$uploadedSearchResultDataAggTable <-
    DT::renderDataTable(selection = "none", {
      req(uploadedStatesData())
      req(uploadedGffData())
      dat <- uploadedSearchResultDataAgg()
      as.data.frame(dat)
    })
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)