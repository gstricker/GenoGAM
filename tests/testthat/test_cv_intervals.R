context("Testing generation of CV intervalls")

test_that("the CV intervals are produced correctly", {
    folds <- 10
    iv <- 20
    cv <- .leaveOutConsecutiveIntervals(folds, iv, 2000)
    expect_equal(length(cv), folds)
    expect_equal(length(cv[[1]]), folds * iv)
})
