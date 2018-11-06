context("Testing GenoGAMList-class")

test_that("GenoGAM builds without parameters", {
    expect_true(is(GenoGAMList(), "GenoGAMList"))
})

test_that("The accessor functions work correctly", {
    ggl <- makeTestGenoGAMList()
    
    expect_error(design())
    expect_true(is(design(ggl), "formula"))
    expect_error(sizeFactors())
    expect_true(is(sizeFactors(ggl), "numeric"))
    expect_error(getSettings())
    expect_true(is(getSettings(ggl), "GenoGAMSettings"))
    expect_error(getFamily())
    expect_true(is(getFamily(ggl), "GenoGAMFamily"))
    expect_error(colData())
    expect_true(is(colData(ggl), "DataFrame"))
    expect_true(all(rownames(colData(ggl)) != rownames(colData(ggl@data[[1]]))))
    expect_error(getParams())
    expect_true(is(getParams(ggl), "list"))
    expect_true(is(getCoefs(ggl), "matrix") | is(getCoefs(ggl), "HDF5Matrix"))
    expect_true(is(getKnots(ggl), "numeric"))
    expect_error(rowRanges())
    expect_true(is(rowRanges(ggl), "list"))
    expect_true(is(rowRanges(ggl)[[1]], "GPos"))
    expect_error(assay())
    expect_true(is(assay(ggl), "list"))
    expect_true(is(assay(ggl)[[1]], "DataFrame"))
    expect_true(all(colnames(ggl) == rownames(colData(ggl@data[[1]])), na.rm = TRUE))
    dnames <- c(names(ggl), rownames(colData(ggl@data[[1]])))
    expect_true(all(dimnames(ggl) != dnames, na.rm = TRUE))
    expect_error(fits())
    expect_error(se())
    
    ggfits <- fits(ggl)
    ggses <- se(ggl)
    expect_true(is(ggfits, "list"))
    expect_true(is(ggfits[[1]], "DataFrame"))
    expect_true(is(ggses, "list"))
    expect_true(is(ggses[[1]], "DataFrame"))
    expect_true(all(names(ggfits[[1]]) == names(ggses[[1]])))

    for(ii in 1:length(ggl@data)) {
        assays(ggl@data[[ii]])[["fits"]] <- NULL
        assays(ggl@data[[ii]])[["se"]] <- NULL
    }
    
    expect_true(all(sapply(assays(ggl), length) == 0))
    expect_true(all(sapply(fits(ggl), is.null)))
    expect_true(all(sapply(se(ggl), is.null)))
})

test_that("Subset on GenoGAM class works correctly", {
    ggl <- makeTestGenoGAMList()
    gr <- GRanges("chrX", IRanges(101, 1000))
    subggl <- ggl[gr]
    subrange <- range(rowRanges(subggl)[[1]])
    expect_true(all(dim(subggl) == c(900, dim(ggl)[2])))
    expect_true(start(subrange) == 101)
    expect_true(end(subrange) == 1000)
    
    ## subggl <- subset(subggl, pos >= 901)
    ## subrange <- range(rowRanges(subggl)[[1]])
    ## expect_true(all(dim(subggl) == c(100, dim(subggl)[2])))
    ## expect_true(start(subrange) == 901)
    ## expect_true(end(subrange) == 1000)
})    
    

