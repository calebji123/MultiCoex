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
#' @param assay_name A string to refer to the name of the assay for the dataset
#'   if the dataset is of type SummarizedExperiment. If no processing is done,
#'   this name must be "tpm"
#' @param TPM_normalize Logical. If TRUE, performs TPM normalization assuming
#'   input data are raw read counts. Default = FALSE.
#' @param gene_lengths A vector of gene lengths in the order of the genes in the dataset.
#'   Must be provided if TPM normalize is set to TRUE. Default = NULL
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
#' net <- DevelopCoexpressionNetwork(GTExBrainTrimmed,
#'                                   cor_method = "pearson",
#'                                   seed = 123)
#' }
#'
#' @export
#' @import BioNERO
#' @import SummarizedExperiment
DevelopCoexpressionNetwork <- function(dataset,
                                       assay_name="tpm",
                                       TPM_normalize = FALSE,
                                       gene_lengths = NULL,
                                       log_scale = FALSE,
                                       cor_method = c("spearman", "pearson", "biweight"),
                                       min_exp = 10,
                                       variance_filter = FALSE,
                                       net_type = "signed",
                                       seed = 123) {
  set.seed(seed)
  cor_method <- match.arg(cor_method)

  # Performing checks of user input
  if (!is.matrix(dataset) && !is.data.frame(dataset) && !(class(dataset) == "SummarizedExperiment")) {
    stop("Input `dataset` must be a numeric matrix or data.frame of gene expression values.")
  }

  # If it is a summarized experiment, the assay must be in the summarized experiment
  # TODO

  if (is.null(rownames(dataset))) {
    stop("Input `dataset` must have rownames representing gene identifiers.")
  }

  if (TPM_normalize && is.null(gene_lengths)) {
    stop("Gene lengths must be provided if TPM_normalize is set to TRUE.")
  }

  # Optional processing to TPM normalize raw counts and log scale TPM
  if (TPM_normalize | log_scale) {
    if (class(dataset) == "SummarizedExperiment") {
      raw_data <- assay(dataset, assay_name)
    }
    else {
      raw_data <- as.data.frame(dataset)
    }
    # Check on the expression matrix
    if (!is.numeric(raw_data)) {
      stop("Expression matrix must be numeric.")
    }
    if (TPM_normalize) {
      rpk <- raw_data / (gene_lengths / 1000)
      per_million <- colSums(rpk) / 1e6
      raw_data <- sweep(rpk, 2, per_million, "/")
    }

    if (log_scale) {
      raw_data <- log2(raw_data + 1)
    }

    if (class(dataset) == "SummarizedExperiment") {
      assays(dataset, "tpm") <- raw_data
    }
    else {
      dataset <- raw_data
    }
  }

  # Run BioNERO pipeline to get the coexpression network
  # Source from the BioNERO vignette
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

# [END]
