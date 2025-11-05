
#' \dontrun {
#' library(SummarizedExperiment)
# library(edgeR)
# library(recount3)
# library(org.Hs.eg.db)
# library(dplyr)
#
# human_projects <- available_projects(organism = "human")
# gtex_projects <- human_projects %>% filter(file_source == "gtex")
#
# proj_info <- subset(gtex_projects, project == "HEART")
# proj_info <- subset(gtex_projects, project == "BRAIN")
# proj_info <- subset(gtex_projects, project == "LIVER")
#
# rse <- create_rse(
#   proj_info,
#   type = "gene",
# )
#
# set.seed(123)
# keep_samp <- sample(seq_len(ncol(rse)), size = min(60, ncol(rse)))
# rse <- rse[, keep_samp]
#
# # Extract Ensembl IDs (without version numbers)
# ens_ids <- sub("\\..*", "", rowData(rse)$gene_id)
#
# biotypes <- mapIds(
#   org.Hs.eg.db,
#   keys = ens_ids,
#   column = "GENETYPE",
#   keytype = "ENSEMBL",
#   multiVals = "first"
# )
#
# keep <- which(biotypes == "protein-coding")
# rse_coding <- rse[keep, ]
#
# gene_length_kb <- rowData(rse_coding)$bp_length / 1000
# raw_counts <- assay(rse, "raw_counts")
#
# rate <- sweep(raw_counts, 1, gene_length_kb, "/")
# tpm <- sweep(rate, 2, colSums(rate), "/") * 1e6
#
# # retain the top 10000 variable genes
# vars <- matrixStats::rowVars(tpm)
# keep_genes <- order(vars, decreasing = TRUE)[1:1000]
# tpm <- tpm[keep_genes,]
#
#
# # Create a simple summary experiment
# GTExLiverTrimmed <- SummarizedExperiment(
#   assays = list(tpm = tpm),
#   rowData = DataFrame(gene_id = rownames(tpm)),
#   colData = DataFrame(sample = colnames(tpm))
# )
#
# usethis::use_data(GTExLiverTrimmed, compress = "xz", overwrite = TRUE)

#' }




