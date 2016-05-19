# setting R_TESTS to empty string because of https://github.com/hadley/testthat/issues/144
Sys.setenv('R_TESTS' = '')

library(testthat)
library(structurer)
test_check("structurer")
