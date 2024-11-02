
## test for the output columns
test_that("Column names are as expected", {
  testthat::expect_equal(colnames(buildDGP(
    export = TRUE,
    patterns = FALSE,
    slopes = FALSE
  )),
  c("pid", "changer", "timing", "t", "y_true", "y_obs"))

  testthat::expect_equal(colnames(buildDGP(
    export = FALSE,
    patterns = TRUE,
    slopes = FALSE
  )), c("patterns", "sim_counts"))

  testthat::expect_equal(colnames(buildDGP(
    export = FALSE,
    patterns = FALSE,
    slopes = TRUE
  )), c("pid", "estimate"))

})

## test for outcome variable
test_that("Outcome variable is binary", {
  testthat::expect_equal(names(table(buildDGP()$y_obs)), c("0", "1"))
})

## test for timing variable
test_that("Time variable is normalized", {
  testthat::expect_equal(max(buildDGP()$t), 1)
  testthat::expect_equal(min(buildDGP()$t), 0)
})

## no NA values on the `sim_counts`
test_that("Pattern generation works as expected", {
  testthat::expect_equal(sum(is.na(
    buildDGP(export = FALSE, patterns = TRUE)$sim_counts
  )), 0)
})
