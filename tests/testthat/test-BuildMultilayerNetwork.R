# Tests for BuildMultilayerNetwork function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for BuildMultilayerNetwork")
test_that("BuildMultilayerNetwork error on invalid user inputs", {
  # non named list input as layer
  testthat::expect_error(ml <- BuildMultilayerNetwork(5))

  # omega is out of bounds
  testthat::expect_error(ml <- BuildMultilayerNetwork(list(layer1 = 1:3, layer2 = 1:3),
                                                      omega = -0.5))

  # Threshold is out of bounds
  testthat::expect_error(ml <- BuildMultilayerNetwork(list(layer1 = 1:3, layer2 = 1:3),
                                                      threshold = -0.5))
})


# [END]
