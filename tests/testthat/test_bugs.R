context("Test for known bugs")

test_that("DelayedMatrix works as usual", {
    vec <- matrix(sample(1:3, 100, replace = TRUE), ncol = 1)
    h5vec <- HDF5Array::writeHDF5Array(vec)
    expect_error(pnorm(-h5vec))
    expect_true(all.equal(pnorm(-vec[,1]), pnorm(-h5vec[,1])))
})
