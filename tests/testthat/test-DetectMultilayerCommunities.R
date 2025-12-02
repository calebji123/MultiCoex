# Tests for DetectMultilayerCommunities function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for detectMultilayerCommunities")
test_that("detectMultilayerCommunities error on invalid user inputs", {
  # not proper input for multilayer network
  testthat::expect_error(net <- detectMultilayerCommunities(5))

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
  communities <- detectMultilayerCommunities(multilayeredNetwork)
  testthat::expect_s3_class(communities, "detectMultilayerCommunitiesResult")
  testthat::expect_length(communities, 3)
  testthat::expect_equal(communities$Q, 0.569111283)
  testthat::expect_equal(dim(communities$stateMembership), c(6837, 4))
  testthat::expect_equal(dim(communities$geneMembership), c(2279, 3))
})


# [END]
