#' Multilayer community detection using Louvain
#'
#' @description
#' Runs Louvain modularity optimization on a multilayer network.
#' This function mirrors the approach of Russell et al. (2023) for tissue-specific
#' co-expression layers.
#'
#' @param ml A `BuildMultilayerNetworkResult` object from `BuildMultilayerNetwork()`.
#' @param omega Inter-layer coupling strength (if overriding ml$params$omega).
#' @param seed Random seed for reproducibility.
#'
#' @return An object of class `DetectMultilayerCommunitiesResult`:
#' \itemize{
#'   \item \code{state_membership} : data.frame of gene × layer to community
#'   \item \code{gene_membership}: consensus across layers
#'   \item \code{Q}: modularity score
#' }
#'
#' @examples
#' \dontrun{
#' # With the multilayered network from BuildMultilayerNetwork()
#' comm <- DetectMultilayerCommunities(ml)
#' }
#'
#' @export
#' @import Matrix
#' @import igraph
DetectMultilayerCommunities <- function(
    ml,
    omega = NULL,
    seed = 123
){
  suppressPackageStartupMessages(requireNamespace("Matrix"))
  # Performing checks of user input
  if (!inherits(ml, "BuildMultilayerNetworkResult")) {
    stop("ml must be a `BuildMultilayerNetworkResult` object, see `BuildMultilayerNetwork()`")
  }

  set.seed(seed)

  # 1) Prepare multilayered adjacency matrix
  A <- ml$supra
  Matrix::diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(as.matrix(A),
                                           mode = "undirected",
                                           weighted = TRUE)

  # 2) Run Louvain on multilayered adjacency matrix
  com <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  mem <- igraph::membership(com)
  Q   <- igraph::modularity(com)

  # 3) Map back to (gene, layer)
  G <- length(ml$genes)
  L <- length(ml$layer_names)
  idx_list <- lapply(seq_len(L), function(li){
    rng <- ml$blocks[[li]]$from:ml$blocks[[li]]$to
    data.frame(
      state_id = rng,
      gene = ml$genes,
      layer = ml$layer_names[li],
      stringsAsFactors = FALSE
    )
  })
  state_map <- do.call(rbind, idx_list)
  state_membership <- transform(
    state_map,
    community = mem[state_map$state_id]
  )

  # 4) Aggregate: consensus community per gene
  tab_list <- table(state_membership$community, state_membership$gene)
  consensus_comm <- apply(tab_list, 2, FUN=function(tb) {
    as.integer(names(tb)[which.max(tb)])
  })
  p_major <- apply(tab_list, 2, function(tb) max(tb)/sum(tb))

  gene_membership <- data.frame(
    gene = names(consensus_comm),
    consensus_comm = consensus_comm,
    p_major = p_major,
    stringsAsFactors = FALSE
  )

  out <- list(
    state_membership = state_membership,
    gene_membership = gene_membership,
    Q = Q
  )
  class(out) <- "DetectMultilayerCommunitiesResult"
  return(out)
}

# [END]
