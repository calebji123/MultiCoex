#' Develop a co-expression network from gene expression data
#'
#' @description
#' This function streamlines the process of converting gene expression data
#' into a processed co-expression network using the BioNERO framework.
#' It automatically checks and formats the input dataset,
#' performs optional TPM normalization and log-transformation, determines
#' an optimal soft-thresholding power to achieve scale-free topology,
#' and constructs a gene co-expression network.
#'
#'
#' @param dataset A numeric matrix or data frame of gene expression values,
#'   with genes in rows and samples in columns. Must have rownames (gene IDs)
#'   and ideally column names (sample IDs).
#' @param TPM_normalize Logical. If TRUE, performs TPM normalization assuming
#'   input data are raw read counts. Default = FALSE.
#' @param log_scale Logical. If TRUE, log2-transforms the data
#'   (after adding a pseudocount of 1). Recommended for TPM/FPKM data. Default = FALSE.
#' @param cor_method Correlation method used to compute co-expression similarity.
#'   One of "spearman", "pearson", or "biweight". Default = "spearman".
#' @param min_exp Minimum expression threshold (passed to BioNERO::exp_preprocess()).
#'   Genes with mean expression below this threshold are filtered out.
#' @param variance_filter Logical. Whether to filter low-variance genes during preprocessing.
#'   Default = FALSE.
#' @param net_type Network type for BioNERO::exp2gcn(). Options: "signed", "unsigned", or "signed hybrid".
#'   Default = "signed hybrid".
#' @param seed Random seed for reproducibility. Default = 123.
#'
#' @return A list (BioNERO network object) containing:
#' \itemize{
#'   \item \code{adjacency}: adjacency matrix of gene co-expression values.
#'   \item \code{modules}: module assignments (if available).
#'   \item \code{SFT_fit}: scale-free topology fit statistics.
#' }
#'
#' @examples
#' \dontrun{
#' # Using GTEx brain tissue dataset available within the package
#' dim(GTExBrainTrimmed)
#' net <- DevelopCoexpressionNetwork(exp_data,
#'                                   TPM_normalize = TRUE,
#'                                   log_scale = TRUE,
#'                                   cor_method = "spearman")
#' }
#'
#' @export
#' @import BioNERO
DevelopCoexpressionNetwork <- function(dataset,
                                       TPM_normalize = FALSE,
                                       log_scale = FALSE,
                                       cor_method = c("spearman", "pearson", "biweight"),
                                       min_exp = 10,
                                       variance_filter = FALSE,
                                       net_type = "signed",
                                       seed = 123) {
  set.seed(seed)
  cor_method <- match.arg(cor_method)

  # Performing checks of user input
  if (!is.matrix(dataset) && !is.data.frame(dataset)) {
    stop("Input `dataset` must be a numeric matrix or data.frame of gene expression values.")
  }
  dataset <- as.data.frame(dataset)

  if (is.null(rownames(dataset))) {
    stop("Input `dataset` must have rownames representing gene identifiers.")
  }

  if (!is.numeric(as.matrix(dataset))) {
    stop("Expression matrix must be numeric.")
  }


  # Optional processing to TPM normalize raw counts and log scale TPM
  if (TPM_normalize) {
    gene_lengths <- attr(dataset, "gene_length")
    if (is.null(gene_lengths)) {
      warning("No 'gene_length' attribute found; TPM normalization may be approximate.")
      gene_lengths <- rep(1, nrow(dataset))
    }
    rpk <- dataset / (gene_lengths / 1000)
    per_million <- colSums(rpk) / 1e6
    dataset <- sweep(rpk, 2, per_million, "/")
  }

  if (log_scale) {
    dataset <- log2(dataset + 1)
  }

  # Run BioNERO pipeline to get the coexpression network
  final_exp <- BioNERO::exp_preprocess(dataset,
                                       min_exp = min_exp,
                                       variance_filter = variance_filter)

  sft <- BioNERO::SFT_fit(final_exp,
                          net_type = net_type,
                          cor_method = cor_method)
  power <- sft$power

  net <- BioNERO::exp2gcn(final_exp,
                          net_type = net_type,
                          SFTpower = power,
                          cor_method = cor_method)
  return(net)
}

