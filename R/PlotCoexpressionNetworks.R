#' Plot a multilayered graph in a neighborhood around a seed gene
#'
#' @description
#' This function will filter for `maxGenes` genes surrounding the `seedGene`in
#' the multilayered network. Then it will plot these genes into a graph,
#' with genes as the nodes and the coexpression relationships as the edges,
#' with relationships from different layers represented in different colors
#' and the width of the edges represent the strength of the coexpression.
#' The filtering occurs in two steps, first neighbours are aggregated based
#' on being less than `k` distance from the `seedGene`, then this list is
#' filtered for the top `maxGenes` by gene connectivity to retain maximum
#' information
#'
#' @param ml A `buildMultilayerNetworkResult` object
#'   from `buildMultilayerNetwork()`.
#' @param seedGene A character vector describing the gene ID to
#'   center the neighborhood on. Ensure this gene is in the network.
#' @param maxGenes An integer representing the maximum number of genes to
#'   include in the plot, excluding the seed gene. Default = 10.
#' @param k An integer representing maximum distance by which a gene is defined
#'   as the `seedGene`'s neighbour. Default = 1 (direct neighbors).
#'
#' @return Invisibly returns the igraph object used for plotting.
#'
#' @export
#' @import igraph
plotCoexpressionNetworks <- function(ml,
                                  seedGene,
                                  maxGenes = 10,
                                  k = 1,
                                  seed = 123) {
  set.seed(seed)
  # Basic checks
  if (!inherits(ml, "BuildMultilayerNetworkResult")) {
    stop("ml must be a `BuildMultilayerNetworkResult` object, see
         `BuildMultilayerNetwork()`")
  }
  if (!is.numeric(maxGenes)) {
    stop("`maxGenes` should be a numeric value that represents the maximum
         amount of genes to be plotted. Default = 10.")
  }
  if (maxGenes <= 0) {
    stop("`maxGenes` must be > 0 otherwise there is nothing to plot.")
  }
  if (!is.numeric(k)) {
    stop("`k` should be a numeric value representing the distance from
         the seed gene to define neighbours. Default = 1.")
  }
  if (k <= 0) {
    stop("`k` must be > 0, as there are no negative distances and k = 0
         would have nothing to plot.")
  }

  genesAll <- ml$genes
  layersAll <- ml$layer_names
  layerAdj <- ml$layer_adj

  if (!seedGene %in% genesAll) {
    stop("Seed gene '", seedGene, "' not found in ml$genes.")
  }

  # 1. Collect all the edges and build a union graph (forgetting about layers)

  # Union edges: ignore layer here, just need connectivity for neighborhood
  # Use the helper adjToEdges
  unionEdges <- do.call(
    rbind,
    lapply(layerAdj, function(A) adjToEdges(A, genesAll))
  )
  unionEdges <- unique(unionEdges[, c("from", "to")])

  gUnion <- igraph::graph_from_data_frame(
    d = unionEdges,
    directed = FALSE,
    vertices = data.frame(name = genesAll, stringsAsFactors = FALSE)
  )

  # 2. Get neighborhood around seed gene
  if (!seedGene %in% igraph::V(gUnion)$name) {
    stop("Seed gene '", seedGene, "' has no edges in the union of all layers")
  }
  # ego nodes are those in the neighbourhood
  egoNodes <- igraph::ego(gUnion, order = k,
                           nodes = seedGene,
                           mode = "all")[[1]]
  # Getting the gene names
  neighGenes <- igraph::V(gUnion)$name[egoNodes]

  # Limit to maxGenes (+ 1 to include the seed)
  if (length(neighGenes) > maxGenes + 1) {
    # keep seed and highest-degree neighbors
    degs <- igraph::degree(gUnion, v = neighGenes)
    # ensure seed is first
    neighGenes <- unique(c(seedGene,
                            setdiff(neighGenes[order(degs, decreasing = TRUE)],
                                    seedGene)))
    neighGenes <- neighGenes[seq_len(min(length(neighGenes), maxGenes))]
  }

  if (length(neighGenes) < 2L) {
    stop("No genes in neighbourhood except for seed. Nothing to plot,
         try another seed gene")
  }

  # 3. Collect edges from each layer restricted to these genes
  layerEdges <- lapply(layersAll, function(ly) {
    A <- as.matrix(layerAdj[[ly]])
    rownames(A) <- colnames(A) <- genesAll
    ASub <- A[neighGenes, neighGenes, drop = FALSE]

    ut <- upper.tri(ASub, diag = FALSE)
    idx <- which(ut & (abs(ASub) > 0))
    if (!length(idx)) {
      return(data.frame(from = character(), to = character(),
                        weight = numeric(), layer = character()))
    }
    rc <- arrayInd(idx, .dim = dim(ASub))
    data.frame(
      from = neighGenes[rc[, 1]],
      to = neighGenes[rc[, 2]],
      weight = ASub[idx],
      layer = ly,
      stringsAsFactors = FALSE
    )
  })

  edgesAll <- do.call(rbind, layerEdges)
  if (nrow(edgesAll) == 0L) {
    message("No nonzero edges among selected genes across layers.")
    return(invisible(NULL))
  }

  # 4. Build multilayer edge-colored graph
  g <- igraph::graph_from_data_frame(
    d = edgesAll,
    directed = FALSE,
    vertices = data.frame(name = neighGenes, stringsAsFactors = FALSE)
  )

  # Edge aesthetics: color by layer, width by |weight|
  layNames <- unique(edgesAll$layer)
  nL <- length(layNames)
  cols <- setNames(grDevices::rainbow(nL), layNames)

  igraph::E(g)$color <- cols[edgesAll$layer]
  w <- abs(edgesAll$weight)
  if (length(w)) {
    igraph::E(g)$width <- 0.5 + 2.5 * (w / max(w, na.rm = TRUE))
  } else {
    igraph::E(g)$width <- 1
  }

  # Highlight seed gene
  vcol <- ifelse(igraph::V(g)$name == seed_gene, "gold", "white")
  vframe <- ifelse(igraph::V(g)$name == seed_gene, "black", "grey40")

  layout <- igraph::layout_with_fr(g)

  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar), add = TRUE)
  par(mar = c(2, 2, 2, 2))

  plot(
    g,
    layout = layout,
    vertex.label = igraph::V(g)$name,
    vertex.size = 10,
    vertex.label.cex = 0.8,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.color = vcol,
    vertex.frame.color = vframe,
    edge.curved = 0.2,
    main = paste0("Local multilayer neighborhood of ", seedGene)
  )

  legend("topleft",
         legend = layNames,
         col = cols[layNames],
         lwd = 2,
         bty = "n",
         title = "Layer")

  invisible(g)
}

# Helper: adjacency -> edge list (upper triangle only)
adjToEdges <- function(A, nodes) {
  A <- as.matrix(A)
  rownames(A) <- colnames(A) <- nodes
  ut <- upper.tri(A, diag = FALSE)
  idx <- which(ut & (abs(A) > 0))
  if (!length(idx)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  rc <- arrayInd(idx, .dim = dim(A))
  data.frame(
    from = nodes[rc[, 1]],
    to   = nodes[rc[, 2]],
    weight = A[idx],
    stringsAsFactors = FALSE
  )
}


# [END]
