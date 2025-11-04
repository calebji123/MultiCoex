#' Heart tissue TPM gene expression data from GTEx trimmed to a reduced size
#'
#' Filtered expression data from GTEx, a comprehensive RNA-seq dataset of
#' adult human tissue. A random 60 samples were chosen to reduce size.
#' Genes were filtered to the ones that code proteins, and only the top
#' 5000 variable genes were maintained to reduce size.
#'
#' @source GTEx Portal from the Broad Institute of MIT and Harvard.
#'
#' @format A SummaryExperiment object with:
#' \describe{
#'  \item{TPM assay}{A matrix of size 5000 x 60 for genes x samples}
#'  \item{colData}{The sample names}
#'  \item{rowData}{The gene names}
#' }
#'
#' @examples
#' \dontrun{
#'  GTExHeartTrimmed
#' }
"GTExHeartTrimmed"

#' Brain tissue TPM gene expression data from GTEx trimmed to a reduced size
#'
#' Filtered expression data from GTEx, a comprehensive RNA-seq dataset of
#' adult human tissue. A random 60 samples were chosen to reduce size.
#' Genes were filtered to the ones that code proteins, and only the top
#' 5000 variable genes were maintained to reduce size.
#'
#' @source GTEx Portal from the Broad Institute of MIT and Harvard.
#'
#' @format A SummaryExperiment object with:
#' \describe{
#'  \item{TPM assay}{A matrix of size 5000 x 60 for genes x samples}
#'  \item{colData}{The sample names}
#'  \item{rowData}{The gene names}
#' }
#'
#' @examples
#' \dontrun{
#'  GTExBrainTrimmed
#' }
"GTExBrainTrimmed"


#' Liver tissue TPM gene expression data from GTEx trimmed to a reduced size
#'
#' Filtered expression data from GTEx, a comprehensive RNA-seq dataset of
#' adult human tissue. A random 60 samples were chosen to reduce size.
#' Genes were filtered to the ones that code proteins, and only the top
#' 5000 variable genes were maintained to reduce size.
#'
#' @source GTEx Portal from the Broad Institute of MIT and Harvard.
#'
#' @format A SummaryExperiment object with:
#' \describe{
#'  \item{TPM assay}{A matrix of size 5000 x 60 for genes x samples}
#'  \item{colData}{The sample names}
#'  \item{rowData}{The gene names}
#' }
#'
#' @examples
#' \dontrun{
#'  GTExLiverTrimmed
#' }
"GTExLiverTrimmed"
# [END]
