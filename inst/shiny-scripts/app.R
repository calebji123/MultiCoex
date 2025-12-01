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
      numericInput("k", "Graph radius (k hops):", 1, min = 1, max = 3)
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
                 br(),
                 tags$p("Step 1: Upload data or choose example data provided.
                         Select the parameters that
                        will be used for developing coexpression adjacency
                        matrices")),
        tabPanel("Data",
                 h3("Expression layers"),
                 tableOutput("exprSummary")),
        tabPanel("Coexpression",
                 h3("Coexpression network summary"),
                 tableOutput("coexSummary")),
        tabPanel("Multilayer & Communities",
                 h3("Multilayer summary"),
                 verbatimTextOutput("mlSummary")),
        tabPanel("Explore seed gene",
                 plotOutput("localPlot", height = "600px"),
                 h4("Genes in selected seed neighborhood"),
                 tableOutput("neighTable")),
        tabPanel("Communities",
                 uiOutput("communityLayerUi"),
                 uiOutput("communityIdUi"),
                 tableOutput("communityTable"))
      )
    )
  )
)

server <- function(input, output, session) {

  # 1. Expression data (Upload or select example)
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
      add_layer <- function(file_input, layer_name) {
        if (is.null(file_input) || layer_name == "") return(NULL)

        # read CSV; assume first column = gene IDs (rownames)
        m <- read.csv(
          file_input$datapath,
          row.names = 1,
          check.names = FALSE
        )

        m <- as.matrix(m)

        if (is.null(rownames(m))) {
          stop("Expression matrix for layer '", layer_name,
               "' must have gene IDs in the first column (used as rownames).")
        }

        layers[[layer_name]] <<- m
      }

      # fill from your UI inputs (adjust IDs to match your app)
      add_layer(input$layer1File, input$layer1Name)
      add_layer(input$layer2File, input$layer2Name)
      add_layer(input$layer3File, input$layer3Name)

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
      n_genes = vapply(ed, nrow, integer(1)),
      n_samples = vapply(ed, ncol, integer(1))
    )
  })

  # 2. Coexpression networks
  coex_nets <- eventReactive(input$runCoex, {
    ed <- exprData()
    req(ed)

    # Run developCoexpressionNetwork for each layer
    nets <- lapply(names(ed), function(ly) {
      message("Building coexpression for layer: ", ly)
      MultiCoex::developCoexpressionNetwork(
        dataset = ed[[ly]],
        TPM_normalize = FALSE,
        log_scale = FALSE,
        cor_method = input$cor_method,
        min_exp = input$min_exp,
        variance_filter = input$variance_filter,
        net_type = input$net_type
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

  output$coexSummary <- renderTable({
    cn <- coexNets()
    req(cn)

    adjList <- cn$adjList
    data.frame(
      layer = names(adjList),
      nGenes = vapply(adjList, nrow, integer(1))
    )
  })

  ## ---- 3. Build multilayer & communities (reactive) --------------------
  ml_obj <- eventReactive(input$run_multilayer, {
    cn <- coex_nets()
    req(cn)

    # Find common genes across layers
    adj_list <- cn$adj_list
    common_genes <- Reduce(intersect, lapply(adj_list, rownames))

    # Build multilayer object (your function)
    MultiCoex::buildMultilayerNetwork(adj_list,
                                      genes = common_genes,
                                      omega = input$omega,
                                      threshold = input$threshold,
                                      match_genes = TRUE)
  })

  comm_res <- eventReactive(input$run_multilayer, {
    ml <- ml_obj()
    req(ml)
    MultiCoex::detectMultilayerCommunities(ml)
  })

  output$ml_summary <- renderPrint({
    ml <- ml_obj()
    cr <- comm_res()
    req(ml, cr)

    cat("Layers:", paste(ml$layer_names, collapse = ", "), "\n")
    cat("Genes:", length(ml$genes), "\n")
    cat("Communities:", length(unique(cr$state_membership$community)), "\n")
  })

  output$pmajor_hist <- renderPlot({
    cr <- comm_res()
    req(cr)
    gm <- cr$gene_membership
    hist(gm$p_major,
         breaks = 20,
         main = "Distribution of p_major",
         xlab = "p_major (community stability)",
         col = "steelblue")
  })

  ## ---- 4. Seed gene UI + local multilayer plot -------------------------
  output$seed_gene_ui <- renderUI({
    ml <- ml_obj()
    req(ml)
    selectInput(
      "seed_gene",
      "Seed gene:",
      choices = sort(ml$genes),
      selected = ml$genes[1]
    )
  })

  output$local_plot <- renderPlot({
    ml <- ml_obj()
    req(ml, input$seed_gene)
    MultiCoex::plotCoexpressionNetworks(
      ml,
      seedGene = input$seed_gene,
      maxGenes = input$max_genes,
      k = input$k_hops
    )
  })

  output$neigh_table <- renderTable({
    # You can adapt plot_local_multilayer to return the vertex set,
    # or reconstruct neighborhood here from ml and seed.
    ml <- ml_obj()
    req(ml, input$seed_gene)

    # Simplest: union graph, k-hop neighborhood
    genes_all <- ml$genes
    layer_adj <- ml$layer_adj
    # quick helper to compute union graph
    adj_to_edges <- function(A, nodes) {
      A <- as.matrix(A)
      rownames(A) <- colnames(A) <- nodes
      ut <- upper.tri(A, diag = FALSE)
      idx <- which(ut & (abs(A) > 0))
      if (!length(idx)) return(data.frame(from=character(), to=character()))
      rc <- arrayInd(idx, .dim = dim(A))
      data.frame(
        from = nodes[rc[,1]],
        to   = nodes[rc[,2]],
        stringsAsFactors = FALSE
      )
    }
    union_edges <- do.call(
      rbind,
      lapply(layer_adj, function(A) adj_to_edges(A, genes_all))
    )
    union_edges <- unique(union_edges)
    g_union <- igraph::graph_from_data_frame(union_edges, directed = FALSE,
                                             vertices = data.frame(name = genes_all))

    if (!input$seed_gene %in% igraph::V(g_union)$name) return(NULL)
    ego_nodes <- igraph::ego(g_union, order = input$k_hops,
                             nodes = input$seed_gene, mode = "all")[[1]]
    neigh_genes <- igraph::V(g_union)$name[ego_nodes]

    head(data.frame(gene = neigh_genes), n = input$max_genes)
  })

  ## ---- 5. Community-level table ----------------------------------------
  output$community_layer_ui <- renderUI({
    ml <- ml_obj()
    req(ml)
    selectInput(
      "community_layer",
      "Layer:",
      choices = ml$layer_names,
      selected = ml$layer_names[1]
    )
  })

  output$community_id_ui <- renderUI({
    cr <- comm_res()
    req(cr, input$community_layer)
    sm <- cr$state_membership
    ly <- subset(sm, layer == input$community_layer)
    comm_ids <- sort(unique(ly$community))
    selectInput(
      "community_id",
      "Community ID:",
      choices = comm_ids,
      selected = comm_ids[1]
    )
  })

  output$community_table <- renderTable({
    cr <- comm_res()
    req(cr, input$community_layer, input$community_id)
    sm <- cr$state_membership
    subset(sm,
           layer == input$community_layer &
             community == input$community_id)[, c("gene","layer","community")]
  })
}

shinyApp(ui = ui, server = server)
# [END]
