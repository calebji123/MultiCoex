#' Build a multilayer co-expression network (supra-adjacency)
#'
#' @description
#' Construct a multilayer network from per-layer adjacency matrices (or from
#' expression matrices by computing correlations first). Each layer is a context
#' (e.g., tissue/timepoint). The result includes a sparse supra-adjacency matrix
#' with omega-coupling between the *same gene* across layers.
#'
#' @param layers A **named** list. Each element is:
#'   a numeric gene x gene adjacency matrix (symmetric)
#'   All layers must have the **same genes and order** (or provide `genes` and
#'   set `matchGenes = TRUE`).
#' @param threshold How to sparsify layer adjacencies:
#'   - numeric in (0,1): keep top n% of edges per layer
#'   - "none": keep weights as-is (dense: not recommended for space).
#' @param omega Inter-layer coupling weight (numeric >= 0).
#' @param genes Optional character vector of gene IDs
#'   (length = nrow of one layer).
#' @param matchGenes Logical; if TRUE and `genes` provided, reorder rows/cols
#'   in each layer to match.
#'
#' @returns An object of class `buildMultilayerNetworkResult`:
#' \itemize{
#'   \item \code{supra} : dgCMatrix (sparse) supra-adjacency of size
#'   (G*L) x (G*L)
#'   \item \code{blocks}: list with block indices (layer offsets) for mapping
#'   \item \code{layerAdj}: list of per-layer sparse adjacencies
#'   (post-threshold)
#'   \item \code{genes}: character vector of genes (length G)
#'   \item \code{layerNames}: character vector (length L)
#' }
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
#' head(multilayered$supra)
#'
#' }
#'
#'
#' @references
#' Russell, M., Aqil, A., Saitou, M., Gokcumen, O., Masuda, N. (2023). Gene
#' communities in co-expression networks across different tissues. PLoS
#' Comput Biol, 19(11). https://doi.org/10.1371/journal.pcbi.1011616
#'
#' Bates, D., Maechler, M., & Jagan, M. (2025). Matrix: Sparse and dense matrix
#' classes and methods (Version 1.7-4).
#' https://doi.org/10.32614/CRAN.package.Matrix
#' https://CRAN.R-project.org/package=Matrix
#'
#' Silva, A. (2022; revised 2025) TestingPackage: An Example R Package For BCB410H.
#' Unpublished. URL https://github.com/anjalisilva/TestingPackage.
#'
#'
#' @export
#' @import Matrix
buildMultilayerNetwork <- function(
    layers,
    threshold = 0.05,
    omega = 0.5,
    genes = NULL,
    matchGenes = FALSE
){
  suppressPackageStartupMessages(requireNamespace("Matrix"))
  # --- 1. Performing checks of user input -----------------)
  if (!is.list(layers)) {
    stop("`layers` must be a named list")
  }
  if (length(layers) < 2) {
    stop("`layers` must have two or more elements to perform multilayered
         analysis")
  }
  if (omega < 0) {
    stop("`omega` must be >= 0")
  }
  layerNames <- names(layers)
  if (is.null(layerNames) || any(layerNames == "")) {
    stop("`layers` must be a *named* list. Names are the IDs, i.e. tissue name,
         developmental stage")
  }
  if (is.na(threshold) || threshold <= 0 || threshold >= 1) {
    stop("Threshold proportion must be in (0,1).")
  }

  # --- 2. Harmonize genes -----------------)
  if (is.null(genes)) {
    # assume first layer defines gene order if none is provided
    genes <- rownames(layers[[1]])
    if (is.null(genes)) {
      stop("One must be provided: either `genes` must not be null, or the first
           element of layers must have gene row names")
    }
  }

  if (matchGenes) {
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
  layerNames <- names(layers)

  # --- 3. Threshold to create a sparse adjacency matrix -----------------)
  layerAdj <- lapply(layers, function(A) .thresholdAdj(A,
                                                        threshold = threshold))

  # --- 3. Build multilayered network -----------------)
  blocks <- vector("list", L)
  blockList <- vector("list", L)
  offset <- 0L
  for (i in seq_len(L)) {
    S <- Matrix::Matrix(layerAdj[[i]], sparse = TRUE)
    blockList[[i]] <- S
    blocks[[i]] <- list(layer = layerNames[i], from = offset + 1L,
                        to = offset + G)
    offset <- offset + G
  }
  # supra --> mathematical representation of a multilayered network
  supra <- Matrix::bdiag(blockList)

  # Inter-layer diagonal couplings: connect same gene across
  # layers with weight omega
  if (omega > 0) {
    # Add omega between (gene g in layer i) and (gene g in layer j)
    # for all i != j
    for (i in seq_len(L - 1L)) {
      for (j in (i + 1L):L) {
        idxI <- blocks[[i]]$from:blocks[[i]]$to
        idxJ <- blocks[[j]]$from:blocks[[j]]$to
        supra[idxI, idxJ] <- supra[idxI, idxJ] + omega * Matrix::Diagonal(G)
        supra[idxJ, idxI] <- supra[idxJ, idxI] + omega * Matrix::Diagonal(G)
      }
    }
  }

  # --- 4. Return -----------------)
  out <- list(
    supra = supra,
    blocks = blocks,
    layerAdj = layerAdj,
    genes = genes,
    layerNames = layerNames
  )
  class(out) <- "buildMultilayerNetworkResult"
  return(out)
}

# Helper function
# Takes an adjacency matrix `A` and a `threshold`
# Thresholds the adjacency matrix so only those greater than or equal to that
# Threshold remains in the matrix
.thresholdAdj <- function(A, threshold = 0.05) {
  diag(A) <- 0
  if (is.numeric(threshold)) {
    # keep top p of *upper triangle* absolute weights
    ut <- A[upper.tri(A, diag = FALSE)]
    k <- ceiling(length(ut) * threshold)
    cut <- sort(abs(ut), decreasing = TRUE, na.last = NA)[k]
    B <- (abs(A) >= cut) * A
  } else if (threshold == "none") {
    B <- A
  } else {
    stop("Unrecognized `threshold`. Use numeric in (0,1), or 'none'.")
  }
  # ensure symmetry
  B <- (B + t(B)) / 2
  return(B)
}

# [END]
