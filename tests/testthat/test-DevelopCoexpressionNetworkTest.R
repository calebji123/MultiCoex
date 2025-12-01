# Tests for developCoexpressionNetwork function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for DevelopCoexpressionNetwork")
test_that("developCoexpressionNetwork error on invalid user inputs", {
  # non matrix or data.fram inputted as dataset
  testthat::expect_error(net <- developCoexpressionNetwork(5))

  # data.frame with no row names inputted as dataset
  testthat::expect_error(net <- developCoexpressionNetwork(
    matrix(1:9, nrow = 3, ncol = 3)))

  # TPM normalize but no gene lengths
  testthat::expect_error(net <-
                           developCoexpressionNetwork(data.frame(first=1:3,
                                                                 second=1:3),
                                                      TPM_normalize = TRUE))
})


# Integration tests
context("Checking that function works with test data not used in example.")

test_that("developCoexpressionNetwork works on valid inputs",{
  net <- developCoexpressionNetwork(
     dataset = GTExBrainTrimmed,
     TPMNormalize = FALSE,
     logScale = FALSE,
     corMethod = "pearson",
     minExp = 10,
     varianceFilter = FALSE,
  )
  testthat::expect_type(net, "list")
  testthat::expect_length(net, 7)
  testthat::expect_equal(dim(net$adjacency_matrix), c(4754, 4754))
})
# [END]
