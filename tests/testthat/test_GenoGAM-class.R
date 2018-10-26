context("Testing GenoGAM-class")

test_that("GenoGAM builds without parameters", {
    expect_true(is(GenoGAM(), "GenoGAM"))
})

test_that("The accessor functions work correctly", {
    gg <- makeTestGenoGAM()
    
    expect_error(design())
    expect_true(is(design(gg), "formula"))
    expect_error(sizeFactors())
    expect_true(is(sizeFactors(gg), "numeric"))
    expect_error(getSettings())
    expect_true(is(getSettings(gg), "GenoGAMSettings"))
    expect_error(getFamily())
    expect_true(is(getFamily(gg), "GenoGAMFamily"))
    expect_error(colData())
    expect_true(is(colData(gg), "DataFrame"))
    expect_true(all(rownames(colData(gg)) != rownames(slot(gg, "colData"))))
    expect_error(getParams())
    expect_true(is(getParams(gg), "list"))
    expect_true(is(getCoefs(gg), "matrix") | is(getCoefs(gg), "HDF5Matrix"))
    expect_true(is(getKnots(gg), "numeric"))
    expect_error(rowRanges())
    expect_true(is(rowRanges(gg), "GPos"))
    expect_error(assay())
    expect_true(is(assay(gg), "DataFrame"))
    expect_true(all(colnames(gg) != rownames(colData(gg)), na.rm = TRUE))
    dnames <- c(names(gg), rownames(slot(gg, "colData")))
    expect_true(all(dimnames(gg) != dnames, na.rm = TRUE))
    expect_error(fits())
    expect_error(se())
    
    ggfits <- fits(gg)
    ggses <- se(gg)
    expect_true(is(ggfits, "DataFrame"))
    expect_true(is(ggses, "DataFrame"))
    expect_true(all(names(ggfits) == names(ggses)))

    assays(gg)[["fits"]] <- NULL
    assays(gg)[["se"]] <- NULL
    expect_true(length(assays(gg)) == 0)
    expect_true(is.null(fits(gg)))
    expect_true(is.null(se(gg)))
})

test_that("Subset on GenoGAM class works correctly", {
    gg <- makeTestGenoGAM()
    gr <- GRanges("chrXYZ", IRanges(101, 1000))
    subgg <- gg[gr]
    subrange <- range(rowRanges(subgg))
    expect_true(all(dim(subgg) == c(900, dim(gg)[2])))
    expect_true(start(subrange) == 101)
    expect_true(end(subrange) == 1000)
    
    subgg <- SummarizedExperiment::subset(subgg, pos >= 901)
    subrange <- range(rowRanges(subgg))
    expect_true(all(dim(subgg) == c(100, dim(subgg)[2])))
    expect_true(start(subrange) == 901)
    expect_true(end(subrange) == 1000)
})    
    


    
