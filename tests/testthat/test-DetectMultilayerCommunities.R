# Tests for DetectMultilayerCommunities function

library(MultiCoex)

# Unit tests

context("Checking for invalid user input for DetectMultilayerCommunities")
test_that("DetectMultilayerCommunities error on invalid user inputs", {
  # not proper input for multilayer network
  testthat::expect_error(net <- DetectMultilayerCommunities(5))

})


# [END]
