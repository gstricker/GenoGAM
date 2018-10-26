context("Test the correctness of some helper functions")

## Unittests
test_that("The parameter validation is correct", {
    l <- list(a = 1, b = 2, d = "a")
    params <- list(a = 5, b = 4, c = 10)
    true_params <- list(a = 1, b = 2, c = 10)
    validParams <- .fillParameters(l, params)
    valid <- sapply(1:length(validParams), function(ii) {
        true_params[[ii]] == validParams[[ii]]
    })
    expect_true(all(valid))
})

test_that("The subsetting coordinates function is correct", {
    ggdl <- makeTestGenoGAMDataSetList()
    coords <- .getCoordinates(ggdl)

    ## single chromosome subset at start
    id <- 1
    i <- start(coords[id,]):end(coords[id,])
    x <- assay(ggdl)
    res <- .subsetByCoords(x, i)
    expect_true(nrow(res) == length(i))
    expect_identical(assay(ggdl)[[1]][i,], res)

    ## single chromosome subset at end
    id <- length(coords)
    i <- start(coords[id,]):end(coords[id,])
    x <- assay(ggdl)
    res <- .subsetByCoords(x, i)
    
    expect_true(nrow(res) == length(i))
    a <- length(assay(ggdl))
    r <- ranges(getIndex(ggdl)[length(getIndex(ggdl))])
    expect_identical(assay(ggdl)[[a]][start(r):end(r),], res)

    ## multiple chromosome subset across entire coordinates
    i <- start(coords[1,]):end(coords[length(coords),])
    x <- assay(ggdl)
    res <- .subsetByCoords(x, i)
    expect_true(nrow(res) == length(i))
    test_res <- do.call("rbind", assay(ggdl))
    expect_identical(test_res, res)
})
