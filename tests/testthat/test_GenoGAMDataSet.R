context("Testing GenoGAMDataSet-class")

test_that("The accessor functions work correctly", {
    ggd <- makeTestGenoGAMDataSet()
    expect_equal(ggd@design, design(ggd))
    expect_equal(ggd@sizeFactors, sizeFactors(ggd))
    sizeFactors(ggd) <- c(10, 10)
    design(ggd) <- ~s(x) + s(x, by = type)
    expect_equal(sizeFactors(ggd), c(10,10))
    expect_equal(design(ggd), ~s(x) + s(x, by = type))
})

test_that("Subset on GenoGAMDataSet class works correctly", {
    ggd <- makeTestGenoGAMDataSet()
    gr <- GRanges("chrI", IRanges(21, 40))
    sub <- subset(ggd, seqnames == "chrII" & pos <= 20)
    ov <- subsetByOverlaps(ggd, gr)
    expect_equal(length(rowRanges(sub)), 20)
    expect_equal(nrow(assay(sub)), 20)
    expect_equal(length(rowRanges(ov)), 37)
    expect_equal(nrow(assay(ov)), 37)
    expect_equal(getIndex(ov)$id, c(2,3,4))
})
