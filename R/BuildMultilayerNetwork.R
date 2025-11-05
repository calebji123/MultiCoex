#' Build a multilayer co-expression network (supra-adjacency)
#'
#' @description
#' Construct a multilayer network from per-layer adjacency matrices (or from
#' expression matrices by computing correlations first). Each layer is a context
#' (e.g., tissue/timepoint). The result includes a sparse supra-adjacency matrix
#' with omega-coupling between the *same gene* across layers.
#'
#' @param layers A **named** list. Each element is either:
#'   - a numeric gene x gene adjacency matrix (symmetric), or
#'   All layers must have the **same genes and order** (or provide `genes` and set `match_genes = TRUE`).
#' @param threshold How to sparsify layer adjacencies:
#'   - numeric in (0,1): keep top n% of edges per layer
#'   - "none": keep weights as-is (dense: not recommended for space).
#' @param omega Inter-layer coupling weight (numeric >= 0).
#' @param genes Optional character vector of gene IDs (length = nrow of one layer).
#' @param match_genes Logical; if TRUE and `genes` provided, reorder rows/cols in each layer to match.
#'
#' @returns An object of class `BuildMultilayerNetworkResult`:
#' \itemize{
#'   \item \code{supra} : dgCMatrix (sparse) supra-adjacency of size (G*L) x (G*L)
#'   \item \code{blocks}: list with block indices (layer offsets) for mapping
#'   \item \code{layer_adj}: list of per-layer sparse adjacencies (post-threshold)
#'   \item \code{genes}: character vector of genes (length G)
#'   \item \code{layer_names}: character vector (length L)
#' }
#'
#' @examples
#' \dontrun{
#' # Using datasets available in the package
#' adj_list <- list(
#'   liver = GTExLiverTrimmed,  # GxG symmetric
#'   brain = GTExBrainTrimmed,
#'   heart = GTExHeartTrimmed
#' )
#' ml <- build_multilayer(adj_list, threshold = "0.05", omega = 0.5)
#' }
#' @export
#' @import Matrix
BuildMultilayerNetwork <- function(
    layers,
    threshold = 0.05,
    omega = 1,
    genes = NULL,
    match_genes = FALSE
){
  suppressPackageStartupMessages(requireNamespace("Matrix"))
  # Performing checks of user input
  if (!is.list(layers)) {
    stop("`layers` must be a named list")
  }
  if (length(layers) < 2) {
    stop("`layers` must have two or more elements to perform multilayered analysis")
  }
  if (omega < 0) {
    stop("`omega` must be >= 0")
  }
  layer_names <- names(layers)
  if (is.null(layer_names) || any(layer_names == "")) {
    stop("`layers` must be a *named* list. Names are the IDs, i.e. tissue name, developmental stage")
  }
  if (is.na(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Threshold proportion must be in (0,1).")
  }

  # 1) Harmonize genes
  if (is.null(genes)) {
    # assume first layer defines gene order if none is provided
    genes <- rownames(layers[[1]])
    if (is.null(genes)) {
      stop("One must be provided: either `genes` must not be null, or the first
           element of layers must have gene row names")
    }
  }

  if (match_genes) {
    layers <- lapply(layers, function(m) {
      rn <- rownames(m)
      cn <- colnames(m)
      if (is.null(rn) || is.null(cn)) {
        stop("All matrices must have row/colnames to match genes.")
      }
      m <- m[genes, genes, drop = FALSE]
      m
    })
  } else {
    # quick sanity: all dimensions equal
    dims <- vapply(layers, nrow, integer(1))
    if (length(unique(dims)) != 1) {
      stop("All layers must have the same number of genes (rows).
           Otherwise, set match_genes to true")
    }
  }

  G <- length(genes)
  L <- length(layers)
  layer_names <- names(layers)

  # 2) Threshold adjacency matrices to create sparse networks
  layer_adj <- lapply(layers, function(A) .threshold_adj(A, threshold = threshold))

  # 3) Building multilayered network
  blocks <- vector("list", L)
  block_list <- vector("list", L)
  offset <- 0L
  for (i in seq_len(L)) {
    S <- Matrix::Matrix(layer_adj[[i]], sparse = TRUE)
    block_list[[i]] <- S
    blocks[[i]] <- list(layer = layer_names[i], from = offset + 1L, to = offset + G)
    offset <- offset + G
  }
  # supra --> mathematical representation of a multilayered network
  supra <- Matrix::bdiag(block_list)

  # Inter-layer diagonal couplings: connect same gene across layers with weight omega
  if (omega > 0) {
    # Add omega between (gene g in layer i) and (gene g in layer j) for all i != j
    for (i in seq_len(L - 1L)) {
      for (j in (i + 1L):L) {
        idx_i <- blocks[[i]]$from:blocks[[i]]$to
        idx_j <- blocks[[j]]$from:blocks[[j]]$to
        supra[idx_i, idx_j] <- supra[idx_i, idx_j] + omega * Matrix::Diagonal(G)
        supra[idx_j, idx_i] <- supra[idx_j, idx_i] + omega * Matrix::Diagonal(G)
      }
    }
  }

  # Return result
  out <- list(
    supra = supra,
    blocks = blocks,
    layer_adj = layer_adj,
    genes = genes,
    layer_names = layer_names
  )
  class(out) <- "BuildMultilayerNetworkResult"
  return(out)
}

# Helper function for thresholding an adjacency matrix
.threshold_adj <- function(A, threshold = 0.05) {
  diag(A) <- 0
  if (is.numeric(threshold)) {
    # keep top p of *upper triangle* absolute weights
    ut <- A[upper.tri(A, diag = FALSE)]
    k <- ceiling(length(ut) * threshold)
    cut <- sort(abs(ut), decreasing = TRUE, na.last = NA)[k]
    B <- (abs(A) >= cut) * A
  } else if (identical(threshold, "none")) {
    B <- A
  } else {
    stop("Unrecognized `threshold`. Use numeric in (0,1), or 'none'.")
  }
  # ensure symmetry
  B <- (B + t(B)) / 2
  return(B)
}

# [END]
