library(shiny)
library(shinyalert)

# 1. Define UI
ui <- fluidPage(
  titlePanel("Multilayered approach to analyzing differential coexpression"),

  # Side bar layout with content
  sidebarLayout(
    sidebarPanel(
      # discussion of steps
      h4("1. Data & Coexpression"),

      radioButtons(
        "dataSource",
        "Expression data source:",
        choices = c("Example GTEx subset" = "example",
                    "Upload up to 3 raw matrices" = "upload"),
        selected = "example"
      ),

      actionButton(inputId = "datasets",
                   label = "Example dataset details"),

      conditionalPanel(
        "input.dataSource == 'upload'",
        h5("Layer 1"),
        textInput("layer1Name", "Layer 1 name:", value = "Brain"),
        fileInput("layer1File", "Layer 1 expression (.csv):",
                  accept = c(".csv")),
        tags$hr(),
        h5("Layer 2 (optional)"),
        textInput("layer2Name", "Layer 2 name:", value = "Heart"),
        fileInput("layer2File", "Layer 2 expression (.csv):",
                  accept = c(".csv")),
        tags$hr(),
        h5("Layer 3 (optional)"),
        textInput("layer3Name", "Layer 3 name:", value = "Liver"),
        fileInput("layer3File", "Layer 3 expression (.csv):",
                  accept = c(".csv"))
      ),

      selectInput(
        "corMethod",
        "Correlation method:",
        choices = c("spearman", "pearson", "biweight"),
        selected = "pearson"
      ),

      numericInput(
        "minExp",
        "Minimum expression (for preprocessing):",
        value = 10, min = 0, step = 1
      ),

      checkboxInput(
        "varianceFilter",
        "Filter low-variance genes",
        value = FALSE
      ),

      selectInput(
        "netType",
        "Network Type:",
        choices = c("signed", "unsigned", "signed hybrid"),
        selected = "signed"
      ),

      actionButton(inputId = "coexParameter",
                   label = "Information on coexpression parameters"),

      hr(),
      tags$p("This may take a while"),
      actionButton("runCoex", "Run coexpression network construction"),

      hr(),
      h4("2. Multilayer & Communities"),
      numericInput(
        "omega",
        "Inter-layer coupling (omega):",
        value = 0.5, min = 0, step = 0.1
      ),
      numericInput(
        "threshold",
        "Coexpression significance threshold filtering:",
        value = 0.05, min = 0, step = 0.01
      ),
      actionButton("runMultilayer", "Build multilayer & detect communities"),

      hr(),
      h4("Exploration options"),
      uiOutput("seedGeneUi"),
      numericInput("maxGenes", "Max genes around seed:", 10, min = 3, max = 50),
      numericInput("k", "Graph radius (k hops):", 1, min = 1, max = 3),

      actionButton(inputId = "multilayerParameter",
                   label = "Information on multilayer and exploration parameters"),
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("About",
                 h3("Shiny app for MultiCoex"),
                 tags$p("Here is a Shiny App for demonstrating visually and
                 interactively the components and functions of the R package
                        MultiCoex."),
                 br(),
                 tags$b("Description: MultiCoex is a tool used for analysis of
                        differential coexpression. It creates multilayered
                        networks by stacking the sparse coexpression networks,
                        and allows for visual inspection of gene neighbourhoods
                        and community based clustering though a Louvain
                        algorithm."),
                 br(),
                 br(),
                 tags$p("Instructions: Use the left hand bar to set the input
                        dataframes and change the parameters of the functions.
                        Run through the functions in order of their steps.
                        Look through the tabs to see explore the data and see
                        the results."),
                 br(),
                 tags$p("Step 1: Upload data or choose example data provided.
                         Select the parameters that
                        will be used for developing coexpression adjacency
                        matrices. Then click the button `Run coexpression
                        network construction`. This process may take a while.
                        Explore the `Coexpression` tab to see a summary of the
                        output when it is done running and explore the relevant
                        table."),
                 br(),
                 tags$p("Step 2: Enter parameters on the left for the multilayer
                        and communities section and click the
                        `Build multilayer & detect communities`
                        button. Then explore the summary and communities in
                        their respective tabs."),
                 br(),
                 tags$p("Step 3: Explore individual gene plots through the
                        exploration options section on the bottom left.
                        Input a gene of interest and the max genes around
                        the seed gene you wish to inspect."),
                 ),
        tabPanel("Data",
                 h3("Expression layers"),
                 tableOutput("exprSummary")),
        tabPanel("Coexpression",
                 h3("Coexpression network summary"),
                 tableOutput("coexSummary"),
                 selectInput(
                   "inspectLayer",
                   "Layer to inspect (this will only show the top 5 rows):",
                   choices = NULL   # we'll fill this from server
                 ),
                 DT::dataTableOutput("adjHead")),
        tabPanel("Multilayer & communities",
                 h3("Multilayer summary"),
                 verbatimTextOutput("mlSummary"),
                 h3("Communities"),
                 uiOutput("communityLayerUi"),
                 uiOutput("communityIdUi"),
                 tableOutput("communityTable")),
        tabPanel("Explore seed gene",
                 plotOutput("localPlot", height = "600px"))
      )
    )
  )
)

server <- function(input, output, session) {

  ## ---- 1. Expression data (upload or example) ------------------------)
  exprData <- reactive({
    if (input$dataSource == "example") {
      # Provided by the package MultiCoex
      # aggregated into a list
      data("GTExLiverTrimmed", package = "MultiCoex", envir = environment())
      data("GTExHeartTrimmed", package = "MultiCoex", envir = environment())
      data("GTExBrainTrimmed", package = "MultiCoex", envir = environment())
      exampleExprList <- list(Brain = GTExBrainTrimmed,
                                Heart = GTExHeartTrimmed,
                                Liver = GTExLiverTrimmed)
      exampleExprList
    } else {
      # read layers
      layers <- list()

      # helper to add one layer if both name + file are present
      addLayer <- function(fileInput, layerName) {
        if (is.null(fileInput) || layerName == "") return(NULL)

        # read CSV; assume first column = gene IDs (rownames)
        m <- read.csv(
          fileInput$datapath,
          row.names = 1,
          check.names = FALSE
        )

        m <- as.matrix(m)

        if (is.null(rownames(m))) {
          stop("Expression matrix for layer '", layerName,
               "' must have gene IDs in the first column (used as rownames).")
        }

        layers[[layerName]] <<- m
      }

      # fill from your UI inputs (adjust IDs to match your app)
      addLayer(input$layer1File, input$layer1Name)
      addLayer(input$layer2File, input$layer2Name)
      addLayer(input$layer3File, input$layer3Name)

      if (length(layers) == 0L) {
        validate("Please upload at least one CSV expression matrix.")
      }

      layers

    }
  })

  # Render summary table of the expression matrices that were used/uploaded
  output$exprSummary <- renderTable({
    ed <- exprData()
    req(ed)
    data.frame(
      layer = names(ed),
      Ngenes = vapply(ed, nrow, integer(1)),
      Nsamples = vapply(ed, ncol, integer(1))
    )
  })

  # Update the selections for later data inspection
  observe({
    ed <- exprData()
    req(ed)
    updateSelectInput(
      session,
      "inspectLayer",
      choices = names(ed),
      selected = names(ed)[1]
    )
  })

  ## ---- 2. Coexpression networks ------------------------)
  coexNets <- eventReactive(input$runCoex, {
    ed <- exprData()
    req(ed)

    # Run developCoexpressionNetwork for each layer
    nets <- lapply(names(ed), function(ly) {
      message("Building coexpression for layer: ", ly)
      MultiCoex::developCoexpressionNetwork(
        dataset = ed[[ly]],
        TPMNormalize = FALSE,
        logScale = FALSE,
        corMethod = input$corMethod,
        minExp = input$minExp,
        varianceFilter = input$varianceFilter,
        netType = input$netType
      )
    })
    names(nets) <- names(ed)

    # Extract adjacency matrices into a simple list
    adjList <- lapply(nets, function(x) x$adjacency)
    list(
      nets = nets,
      adjList = adjList
    )
  })

  # Render a small summary of the coexpression data for each layer
  output$coexSummary <- renderTable({
    cn <- coexNets()
    req(cn)

    adjList <- cn$adjList
    data.frame(
      layer = names(adjList),
      nGenes = vapply(adjList, nrow, integer(1))
    )
  })

  # Render the head of the adjacency table for the layer selected
  output$adjHead <- DT::renderDataTable({
    cn <- coexNets()
    req(cn, input$inspectLayer)
    adj <- cn$adjList[[input$inspectLayer]]
    DT::datatable(
      head(as.data.frame(adj)),
      options = list(pageLength = 5, scrollX = TRUE)
    )
  })

  ## ---- 3. Build multilayer & communities (reactive) -------------------)
  mlObj <- eventReactive(input$runMultilayer, {
    cn <- coexNets()
    req(cn)

    # Find common genes across layers
    adjList <- cn$adjList
    commonGenes <- Reduce(intersect, lapply(adjList, rownames))

    # Build multilayer object (your function)
    MultiCoex::buildMultilayerNetwork(adjList,
                                      genes = commonGenes,
                                      omega = input$omega,
                                      threshold = input$threshold,
                                      matchGenes = TRUE)
  })

  commRes <- eventReactive(input$runMultilayer, {
    ml <- mlObj()
    req(ml)
    MultiCoex::detectMultilayerCommunities(ml)
  })

  output$mlSummary <- renderPrint({
    ml <- mlObj()
    cr <- commRes()
    req(ml, cr)

    cat("Layers:", paste(ml$layerNames, collapse = ", "), "\n")
    cat("Genes:", length(ml$genes), "\n")
    cat("Communities:", length(unique(cr$stateMembership$community)), "\n")
  })

  ## ---- 4. Seed gene UI + local multilayer plot ------------------------)
  output$seedGeneUi <- renderUI({
    ml <- mlObj()
    req(ml)
    selectInput(
      "seedGene",
      "Seed gene:",
      choices = sort(ml$genes),
      selected = ml$genes[1]
    )
  })

  output$localPlot <- renderPlot({
    ml <- mlObj()
    req(ml, input$seedGene)
    MultiCoex::plotCoexpressionNetworks(
      ml,
      seedGene = input$seedGene,
      maxGenes = input$maxGenes,
      k = input$k
    )
  })


  ## ---- 5. Community-level table ---------------------------------------)
  output$communityLayerUi <- renderUI({
    ml <- mlObj()
    req(ml)
    selectInput(
      "communityLayer",
      "Layer:",
      choices = ml$layerNames,
      selected = ml$layerNames[1]
    )
  })

  output$communityIdUi <- renderUI({
    cr <- commRes()
    req(cr, input$communityLayer)
    sm <- cr$stateMembership
    ly <- subset(sm, layer == input$communityLayer)
    commIds <- sort(unique(ly$community))
    selectInput(
      "communityId",
      "Community ID:",
      choices = commIds,
      selected = commIds[1]
    )
  })

  output$communityTable <- renderTable({
    cr <- commRes()
    req(cr, input$communityLayer, input$communityId)
    sm <- cr$stateMembership
    subset(sm,
           layer == input$communityLayer &
             community == input$communityId)[, c("gene","layer","community")]
  })

  ## ---- 6. Catch alerts ---------------------------------------)
  observeEvent(input$datasets, {
    # Show a modal when the button is pressed
    shinyalert::shinyalert(title = "Example Datasets from GTEx",
               text = "Expression datasets from the brain, heart, and liver taken from GTEx, a comprehensive human RNA-seq dataset.
               This data has been trimmed for a smaller size.

               Citation: The GTEx Consortium. (2020). The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science, 369(6509), 1318–1330. https://doi.org/10.1126/science.aaz1776 ",
               type = "info")
  })
  observeEvent(input$coexParameter, {
    # Show a modal when the button is pressed
    shinyalert::shinyalert(title = "Run coexpression parameters",
               text = "Parameters used to create coexpression data information:

               Uploaded data must be TPM and log normalized data, with up to three layers each requiring a separate file. These files are gene expression files with genes as the rows and samples as the columns. These will be inputted as .csv.

               Correlation method describes the statistical test used to correlate gene expression. Use the default (pearson) for the fastest test on the example data.

               Minimum expression describes how much expression (in TPM) is required to retain that gene in the analysis.

               Filter low variance genes is a parameter that removes genes without much variance between samples. Keep this off for the example data because the trimmed data means there is not much variance retained.

               Network type is for how the gene-gene adjacency values are represented in the resulting output. Default is signed.
               ",
               type = "info")
  })
  observeEvent(input$multilayerParameter, {
    # Show a modal when the button is pressed
    shinyalert::shinyalert(title = "Run multilayer, communities, and exploration parameters",
               text = "Parameters used to create multilayer network, find gene communities, and explore graph:

               Inter-layer coupling (omega) is the strength of the edge between layers. This is usually calcualted based on the average edge weight of the network, but to get a naive sense of the algorithm you can use a value close to the size of the weights between genes as seen in the adjacency matrix.

               Coexpression signficiance threshold filtering is how many genes are filtered out to create the coexpression graph, a number from 0 - 1.

               Seed gene is the gene to start the plot analysis. This must be selected from available genes kept in analysis.

               Max genes around seed is a number showing the amount of genes that will be displayed along with the seed gene.

               Graph radius (k hops) is a number representing how far away the neighbor genes will be pulled to satisfy the max genes.
               ",
               type = "info")
  })



}

shinyApp(ui = ui, server = server)
# [END]
