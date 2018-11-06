## context("Test the mergeRanges function")

## test_that("ranges are merged correctly", {
##           gr <- GRanges(c("chrI", "chrI", "chrI"),
##                         IRanges(c(1,10,27), c(20,25,40)))
##           mgr <- mergeRanges(gr)
##           mgr2 <- mergeRanges(gr, 5)
##           expect_equal(length(mgr), 2)
##           expect_equal(length(mgr2), 1)
## })
