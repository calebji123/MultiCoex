# Tests for DevelopCoexpressionNetwork function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for DevelopCoexpressionNetwork")
test_that("DevelopCoexpressionNetwork error on invalid user inputs", {
  # non matrix or data.fram inputted as dataset
  testthat::expect_error(net <- DevelopCoexpressionNetwork(5))

  # data.frame with no row names inputted as dataset
  testthat::expect_error(net <- DevelopCoexpressionNetwork(matrix(1:9, nrow = 3, ncol = 3)))

  # TPM normalize but no gene lengths
  testthat::expect_error(net <- DevelopCoexpressionNetwork(data.frame(first=1:3,
                                                                      second=1:3),
                                                           TPM_normalize = TRUE))
})


# [END]
