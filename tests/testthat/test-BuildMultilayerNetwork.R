# Tests for buildMultilayerNetwork function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for BuildMultilayerNetwork")
test_that("buildMultilayerNetwork error on invalid user inputs", {
  # non named list input as layer
  testthat::expect_error(ml <- buildMultilayerNetwork(5))

  # omega is out of bounds
  testthat::expect_error(ml <- buildMultilayerNetwork(list(layer1 = 1:3, layer2 = 1:3),
                                                      omega = -0.5))

  # Threshold is out of bounds
  testthat::expect_error(ml <- buildMultilayerNetwork(list(layer1 = 1:3, layer2 = 1:3),
                                                      threshold = -0.5))
})


# Integration tests
context("Checking that function works with test data.")

test_that("buildMultilayerNetwork works on valid inputs",{
  brainNet <- developCoexpressionNetwork(
                 SummarizedExperiment::assay(GTExBrainTrimmed, "tpm"),
                 corMethod = "pearson")
  heartNet <- developCoexpressionNetwork(
                 SummarizedExperiment::assay(GTExHeartTrimmed, "tpm"),
                 corMethod = "pearson")
  liverNet <- developCoexpressionNetwork(
                 SummarizedExperiment::assay(GTExLiverTrimmed, "tpm"),
                 corMethod = "pearson")

  # Combine them together to prepare them for the multilayer network
  layers <- list(Brain = brainNet$adjacency_matrix,
                 Heart = heartNet$adjacency_matrix,
                 Liver = liverNet$adjacency_matrix)
  # Find the common genes
  commonGenes <- Reduce(intersect, list(
                        rownames(brainNet$adjacency_matrix),
                        rownames(heartNet$adjacency_matrix),
                        rownames(liverNet$adjacency_matrix)
                       ))
  # Building the multilayered network
  multilayeredNetwork <- buildMultilayerNetwork(layers,
                                                genes = commonGenes,
                                                matchGenes = TRUE,
                                                omega = 0.5,
                                                threshold = 0.05)
  testthat::expect_s3_class(multilayeredNetwork, "buildMultilayerNetworkResult")
  testthat::expect_length(multilayeredNetwork, 5)
  testthat::expect_length(multilayeredNetwork$layerNames, 3)
  testthat::expect_length(multilayeredNetwork$genes, 2279)
  testthat::expect_equal(dim(multilayeredNetwork$supra), c(6837, 6837))
  testthat::expect_length(multilayeredNetwork$blocks, 3)
})

# [END]
