context("Test size factor computation")

test_that("Size factor computation is correct", {

    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    vec <- round(runif(dim(ggd)[1], 1, 10))
    assay(ggd)[,1] <- vec
    assay(ggd)[,2] <- 2*vec
    slot(ggd, "countMatrix") <- sum(ggd)
    
    ggd <- computeSizeFactors(ggd)
    sf <- sizeFactors(ggd)
    assay(ggd)[,1] <- assay(ggd)[,1]/exp(sf[1])
    assay(ggd)[,2] <- assay(ggd)[,2]/exp(sf[2])
    expect_true(all.equal(sum(assay(ggd)[,1]), sum(assay(ggd)[,2])))
})

test_that("Size factor function work with missing input", {
    ggd <- GenoGAMDataSet()
    new_ggd <- computeSizeFactors(ggd)

    expect_true(identical(ggd, new_ggd))

    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    colnames(ggd) <- NULL
    new_ggd <- computeSizeFactors(ggd)
    expect_true(identical(ggd, new_ggd))
})

test_that("Size factor function work in case of single input", {
    gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
    seqlengths(gr) <- 1e6
    df <- DataFrame(colA = 1:10000)
    se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
    ggd <- GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                      design = ~ s(x))

    expect_error(supressWarnings(computeSizeFactors(ggd)))
})


