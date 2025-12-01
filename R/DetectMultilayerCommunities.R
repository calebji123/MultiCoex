#' Multilayer community detection using Louvain Algorithm
#'
#' @description
#' This function runs Louvain modularity optimization on a multilayer network.
#' It mirrors the approach of Russell et al. (2023), which used an iterated
#' GenLouvain function from a matlab github toolkit. tissue-specific
#' co-expression layers.
#'
#' @param ml A `buildMultilayerNetworkResult` object
#'    from `buildMultilayerNetwork()`.
#' @param omega Inter-layer coupling strength (if overriding ml$params$omega).
#' @param seed Random seed for reproducibility.
#'
#' @return An object of class `detectMultilayerCommunitiesResult`:
#' \itemize{
#'   \item \code{stateMembership} : data.frame of gene × layer to community
#'   \item \code{geneMembership}: consensus across layers
#'   \item \code{Q}: modularity score
#' }
#'
#' @examples
#' \dontrun{
#' # Using dummy datasets to keep the example code in the functions
#' # fast and consistent
#' genes <- paste0("Gene", 1:10)
#'
#' Brain  <- matrix(rnorm(10 * 5), nrow = 10,
#'                  dimnames = list(genes, paste0("B",1:5)))
#' Heart  <- matrix(rnorm(10 * 5), nrow = 10,
#'                  dimnames = list(genes, paste0("H",1:5)))
#' Liver  <- matrix(rnorm(10 * 5), nrow = 10,
#'                  dimnames = list(genes, paste0("L",1:5)))
#'
#' adjList <- list(Brain = Brain, Heart = Heart, Liver = Liver)
#'
#' ml <- buildMultilayerNetwork(adjList, threshold = 0.05, omega = 0.5)
#' comm <- detectMultilayerCommunities(ml)
#' }
#'
#' @export
#' @import Matrix
#' @import igraph
detectMultilayerCommunities <- function(
    ml,
    omega = NULL,
    seed = 123
){
  suppressPackageStartupMessages(requireNamespace("Matrix"))
  # --- 1. Performs checks on input -----------------)
  if (!inherits(ml, "buildMultilayerNetworkResult")) {
    stop("ml must be a `buildMultilayerNetworkResult` object, see `buildMultilayerNetwork()`")
  }

  set.seed(seed)

  # --- 2. Prepare multilayered adjacency matrix -----------------)
  A <- ml$supra
  Matrix::diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(as.matrix(A),
                                           mode = "undirected",
                                           weighted = TRUE)

  # --- 3. Run louvain on adjacency matrix -----------------)
  com <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  mem <- igraph::membership(com)
  Q   <- igraph::modularity(com)

  # --- 4. Map back to gene, layer -----------------)
  G <- length(ml$genes)
  L <- length(ml$layerNames)
  idxList <- lapply(seq_len(L), function(li){
    rng <- ml$blocks[[li]]$from:ml$blocks[[li]]$to
    data.frame(
      stateId = rng,
      gene = ml$genes,
      layer = ml$layerNames[li],
      stringsAsFactors = FALSE
    )
  })
  stateMap <- do.call(rbind, idxList)
  stateMembership <- transform(
    stateMap,
    community = mem[stateMap$stateId]
  )

  # --- 5. Aggregate and find consensus community -----------------)
  tabList <- table(stateMembership$community, stateMembership$gene)
  consensusComm <- apply(tabList, 2, FUN=function(tb) {
    as.integer(names(tb)[which.max(tb)])
  })
  pMajor <- apply(tabList, 2, function(tb) max(tb)/sum(tb))

  geneMembership <- data.frame(
    gene = names(consensusComm),
    consensusComm = consensusComm,
    pMajor = pMajor,
    stringsAsFactors = FALSE
  )

  # --- 6. Return -----------------)
  out <- list(
    stateMembership = stateMembership,
    geneMembership = geneMembership,
    Q = Q
  )
  class(out) <- "detectMultilayerCommunitiesResult"
  return(out)
}

# [END]
