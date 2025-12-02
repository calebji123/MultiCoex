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
#' @param legendPosition A character vector describing the position of the
#'   legend in the figure. This is for custom positioning when the graph extends
#'   into the legend. One of "topleft", "topright", "bottomleft", "bottomright".
#'   Defailt = "topleft"
#' @param seed Random seed for reproducibility. Default = 123.
#'
#' @return Invisibly returns the igraph object used for plotting.
#'
#' @examples
#' \dontrun{
#' # Using example data from GTEx, trimmed to be lightweight and usable for
#' # examples
#' ?GTExBrainTrimmed
#' ?GTExHeartTrimmed
#' ?GTExLiverTrimmed
#' # Get the coexpression data for each dataset
#' brainNet <- developCoexpressionNetwork(
#'                SummarizedExperiment::assay(GTExBrainTrimmed, "tpm"),
#'                corMethod = "pearson")
#' heartNet <- developCoexpressionNetwork(
#'                SummarizedExperiment::assay(GTExHeartTrimmed, "tpm"),
#'                corMethod = "pearson")
#' liverNet <- developCoexpressionNetwork(
#'                SummarizedExperiment::assay(GTExLiverTrimmed, "tpm"),
#'                corMethod = "pearson")
#'
#' # Combine them together to prepare them for the multilayer network
#' layers <- list(Brain = brainNet$adjacency_matrix,
#'                Heart = heartNet$adjacency_matrix,
#'                Liver = liverNet$adjacency_matrix)
#' # Find the common genes
#' commonGenes <- Reduce(intersect, list(
#'                       rownames(brainNet$adjacency_matrix),
#'                       rownames(heartNet$adjacency_matrix),
#'                       rownames(liverNet$adjacency_matrix)
#'                      ))
#' #Building the multilayered network
#' multilayeredNetwork <- buildMultilayerNetwork(layers,
#'                                               genes = commonGenes,
#'                                               matchGenes = TRUE,
#'                                               omega = 0.5,
#'                                               threshold = 0.05)
#'
#' MultiCoex::plotCoexpressionNetworks(multilayeredNetwork,
#'                                     seedGene = "ENSG00000210082.2",
#'                                     maxGenes = 10,
#'                                     k = 1,
#'                                     legendPosition = "topright")
#'
#' }
#'
#' @references
#'
#' Silva, A. (2022; revised 2025) TestingPackage: An Example R Package For BCB410H.
#' Unpublished. URL https://github.com/anjalisilva/TestingPackage.
#'
#' Csárdi, G., & Nepusz, T. (2006). The igraph software package for complex
#' network research. InterJournal, Complex Systems, 1695. https://igraph.org
#'
#' Antonov, M., Csárdi, G., Horvát, S., Müller, K., Nepusz, T., Noom, D.,
#' Salmon, M., Traag, V., Welles, B. F., & Zanini, F. (2023). igraph enables
#' fast and robust network analysis across programming languages.
#' arXiv:2311.10260. https://doi.org/10.48550/arXiv.2311.10260
#'
#' Csárdi, G., Nepusz, T., Traag, V., Horvát, S., Zanini, F., Noom, D.,
#' Müller, K., Schoch, D., & Salmon, M. (2025). igraph: Network analysis and
#' visualization in R (Version 2.2.1). https://doi.org/10.5281/zenodo.7682609
#' https://CRAN.R-project.org/package=igraph
#'
#' @export
#' @import igraph
#' @import grDevices
#' @importFrom graphics legend par
#' @importFrom stats setNames
plotCoexpressionNetworks <- function(ml,
                                  seedGene,
                                  maxGenes = 10,
                                  k = 1,
                                  legendPosition = c("topleft", "topright",
                                                     "bottomleft",
                                                     "bottomright"),
                                  seed = 123) {
  set.seed(seed)
  legendPosition <- match.arg(legendPosition)

  # --- 1. Performs checks in inputs -----------------)
  if (!inherits(ml, "buildMultilayerNetworkResult")) {
    stop("ml must be a `buildMultilayerNetworkResult` object, see
         `buildMultilayerNetwork()`")
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
  layersAll <- ml$layerNames
  layerAdj <- ml$layerAdj

  if (!seedGene %in% genesAll) {
    stop("Seed gene '", seedGene, "' not found in ml$genes.")
  }

  # --- 2. Collect all edges to form a union graph -----------------)

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

  # --- 2. Get neighborhood around seed gene -----------------)
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

  # --- 3. Collect edges from each layer -----------------)
  # Restricted to the neighbourhood genes
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

  # --- 4. Build multilayer graph plot -----------------)
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


  # Change the curvature based on how many edges are between two nodes
  el <- as.data.frame(igraph::as_edgelist(g, names = TRUE))
  colnames(el) <- c("from", "to")
  pairId <- paste(pmin(el$from, el$to), pmax(el$from, el$to), sep = "||")

  # Initialize curvature vector of correct length
  curvVals <- numeric(igraph::ecount(g))

  # For each pair, if there are multiple edges, give them different curvature
  offset <- 0.2  # how "curvy" parallel edges are
  for (pid in unique(pairId)) {
    idx <- which(pairId == pid)
    n   <- length(idx)
    if (n == 1L) {
      curvVals[idx] <- 0.2          # normal edge
    } else {
      curvVals[idx] <- seq(-offset, offset, length.out = n)
    }
  }

  igraph::E(g)$curved <- curvVals

  # Highlight seed gene
  vcol <- ifelse(igraph::V(g)$name == seedGene, "gold", "white")
  vframe <- ifelse(igraph::V(g)$name == seedGene, "black", "grey40")

  layout <- igraph::layout_with_fr(g)

  oldPar <- par(no.readonly = TRUE)
  on.exit(par(oldPar), add = TRUE)
  par(mar = c(2, 2, 2, 2))

  # --- 5. Plot the graph -----------------)
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
    main = paste0("Local multilayer neighborhood of ", seedGene)
  )

  legend(legendPosition,
         legend = layNames,
         col = cols[layNames],
         lwd = 2,
         bty = "n",
         title = "Layer")

  invisible(g)
}

# Helper function: taking in an adjacency matrix `A`
# and a list `nodes`, return a list of edges and weights
# By parsing through the upper triangle
adjToEdges <- function(A, nodes) {
  A <- as.matrix(A)
  rownames(A) <- colnames(A) <- nodes
  ut <- upper.tri(A, diag = FALSE)
  idx <- which(ut & (abs(A) > 0))
  if (!length(idx)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  rc <- arrayInd(idx, .dim = dim(A))
  return(data.frame(
    from = nodes[rc[, 1]],
    to   = nodes[rc[, 2]],
    weight = A[idx],
    stringsAsFactors = FALSE
  ))
}


# [END]
