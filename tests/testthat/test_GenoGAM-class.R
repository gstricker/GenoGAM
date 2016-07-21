context("Testing GenoGAM-class")

test_that("The accessor functions work correctly", {
    gg <- makeTestGenoGAM()
    expect_equal(gg@positions, rowRanges(gg))
    expect_equal(gg@experimentDesign, design(gg))
    expect_equal(gg@fits, getFits(gg))
})

test_that("Subset on GenoGAM class works correctly", {
    gg <- makeTestGenoGAM()
    gr <- GRanges("chrI", IRanges(31, 50))
    sub <- subset(gg, pos >= 61)
    ov <- subsetByOverlaps(gg, gr)
    expect_equal(length(rowRanges(sub)), 40)
    expect_equal(nrow(getFits(sub)), 40)
    expect_equal(length(rowRanges(ov)), 20)
    expect_equal(nrow(getFits(ov)), 20)
})

test_that("Significance computation is working", {
    gg <- makeTestGenoGAM()
    gg <- computeSignificance(gg)
    expect_equal(ncol(getFits(gg)), 6)
    expect_is(getFits(gg)[,5], "numeric")
    expect_is(getFits(gg)[,6], "numeric")
})
    
    


    
