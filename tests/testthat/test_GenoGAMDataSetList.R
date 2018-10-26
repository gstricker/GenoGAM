context("Testing the GenoGAMDataSetList-class")

## Unittests
wd <- system.file("extdata/Set1", package='GenoGAM')
dir <- file.path(wd, "bam")
config <- file.path(wd, "experimentDesign.txt")

region <- GRanges("chrI", IRanges(10000, 20000))
params <- Rsamtools::ScanBamParam(which = region)
settings <- GenoGAMSettings(bamParams = params)
## GenoGAMDataSetList
test_ggdl <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                       design = ~ s(x) + s(x, by = genotype), directory = dir,
                       settings = settings, split = TRUE)


test_that("The constructor works correctly in default mode", {
    df <- GenoGAMDataSetList()
    expect_identical(new("GenoGAMDataSetList"), df)
})

test_that("The constructor works correctly for config files and data.frames", {
    config <- system.file("extdata/Set1", "experimentDesign.txt", package = "GenoGAM")
    dir <- system.file("extdata/Set1/bam", package = "GenoGAM")

    region <- GRanges("chrI", IRanges(10000, 20000))
    params <- Rsamtools::ScanBamParam(which = region)
    settings <- GenoGAMSettings(bamParams = params)
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings, split = TRUE)
    expect_true(checkObject(ds))
    expect_equal(getTileNumber(ds), 10)
    expect_true(all(sapply(assay(ds)[[1]], sum) > 0))
    
    df <- read.table(config, header = TRUE, sep = '\t')
    ds <- GenoGAMDataSet(df, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings, split = TRUE)
    expect_true(checkObject(ds))
    expect_equal(getTileNumber(ds), 10)
    expect_true(all(sapply(assay(ds)[[1]], sum) > 0))

    region = GRanges("chrV", IRanges(105000, 108000))
    params <- Rsamtools::ScanBamParam(which = region)
    settings <- GenoGAMSettings(bamParams = params)
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings, split = TRUE)
    expect_false(checkObject(ds))
    expect_true(is.null(getTileNumber(ds)))
    expect_true(length(assays(ds)) == 0)

    gr <- GenomicRanges::GRanges("chrV", IRanges::IRanges(1, 100000))
    settings <- GenoGAMSettings(bamParams = Rsamtools::ScanBamParam(which = gr))
    
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings, split = TRUE)
    expect_false(checkObject(ds))
    expect_true(is.null(getTileNumber(ds)))
    expect_true(length(assays(ds)) == 0)

    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = myColumn), directory = dir,
                         split = TRUE)
    expect_false(checkObject(ds))
})


test_that("Settings checking functions work correctly", {
    ggdl <- makeTestGenoGAMDataSetList()

    expect_true(checkObject(ggdl))
    
    metadata(slot(ggdl, "index"))$check <- FALSE
    expect_false(checkObject(ggdl))

    metadata(slot(ggdl, "index"))$check <- NULL
    expect_false(checkObject(ggdl))

    metadata(slot(ggdl, "index"))$chunkSize <- 100
    expect_true(is(.checkChunkSize(ggdl), "character"))

    metadata(slot(ggdl, "index"))$tileSize <- 1000
    expect_true(is(.checkTileSize(ggdl), "character"))

    end(slot(ggdl, "index"))[1] <- 10000
    expect_true(is(.checkEqualityOfTiles(ggdl), "character"))

    gr2 <- GRanges("chr2", IRanges(1, 1000), id = 30)
    seqlengths(gr2) <- 1e6
    md <- metadata(slot(ggdl, "index"))
    suppressWarnings(slot(ggdl, "index") <- c(slot(test_ggdl, "index"), gr2))
    metadata(slot(ggdl, "index")) <- md
    expect_true(is(.checkChromosomes(ggdl), "character"))

    expect_true(is(.checkNumberOfTiles(ggdl), "character"))

    metadata(slot(ggdl, "index"))$chromosomes <- gr2
    expect_true(is(.checkTileRanges(ggdl), "character"))
})

test_that("Accessors return the right slots in line with GenoGAMDataSet", {

    config <- system.file("extdata/Set1", "experimentDesign.txt", package = "GenoGAM")
    dir <- system.file("extdata/Set1/bam", package = "GenoGAM")

    region <- GRanges("chrI", IRanges(10000, 20000))
    params <- Rsamtools::ScanBamParam(which = region)
    settings <- GenoGAMSettings(bamParams = params)

    ## GenoGAMDataSet
    ggd <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings)

    ggdl <- test_ggdl
    
    expect_identical(getIndex(ggd), getIndex(ggdl))
    expect_identical(tileSettings(ggd), tileSettings(ggdl))
    expect_identical(getChromosomes(ggd), getChromosomes(ggdl))
    expect_identical(getTileSize(ggd), getTileSize(ggdl))
    expect_identical(getChunkSize(ggd), getChunkSize(ggdl))
    expect_identical(getOverhangSize(ggd), getOverhangSize(ggdl))
    expect_identical(getTileNumber(ggd), getTileNumber(ggdl))
    expect_true(design(ggd) == design(ggdl))
    expect_true(all(sizeFactors(ggd) == sizeFactors(ggdl)))
    expect_true(length(ggd) == length(ggdl))
    expect_true(all(dim(ggd) == dim(ggdl)))
    expect_true(all(seqlengths(ggd) == seqlengths(ggdl)))
    expect_true(all(seqlevels(ggd) == seqlevels(ggdl)))
    expect_identical(colData(ggd), colData(ggdl))
    expect_identical(rowRanges(ggd), do.call("c", rowRanges(ggdl))[[1]])
    expect_identical(assay(ggd), assay(ggdl)[[1]])
    expect_identical(assays(ggd), assays(ggdl)[[1]])
    expect_true(all(colnames(ggd) == colnames(ggdl)))

    expect_error(design(ggdl) <- ~ s(x) + s(by = myColumn))
    getChunkSize(ggdl) <- 2000
    expect_equal(getTileNumber(ggdl), 5)
    getTileSize(ggdl) <- 10400
    expect_equal(getTileNumber(ggdl), 1)
    getOverhangSize(ggdl) <- 500
    expect_equal(median(width(getIndex(ggdl))), 10001)
    getTileNumber(ggdl) <- 10
    expect_equal(getChunkSize(ggdl), getChunkSize(ggd))
})

test_that("The subsetting methods work correct", {

    ## normal subsetting reduces the list of elements to those used only
    ggdl <- makeTestGenoGAMDataSetList()
    gr <- GRanges(c("chrX", "chrY"), IRanges(c(1,1), c(2000,2000)))
    res <- subsetByOverlaps(ggdl, gr)
    gr$id <- 1:2
    expect_true(all(seqlevelsInUse(res) == c("chrX", "chrY")))
    expect_true(all(gr == getIndex(res)))
    expect_true(length(rowRanges(res)) == 2)
    expect_true(length(assay(res)) == 2)

    ## res <- subset(ggdl, seqnames == "chrX")
    ## expect_true(seqlevelsInUse(res) == "chrX")
    ## expect_true(length(rowRanges(res)) == 1)
    ## expect_true(length(assay(res)) == 1)
    ## expect_identical(.extractGR(rowRanges(res)[[1]]), getChromosomes(res))

    ## further tests
    ggdl <- test_ggdl
    getOverhangSize(ggdl) <- 0
    gr <- GRanges("chrI", IRanges(11501,12500))
    seqlengths(gr) <- seqlengths(ggdl)

    res <- subsetByOverlaps(ggdl, gr)
    dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res)[[1]])),
                                 IRanges::ranges(rowRanges(res)[[1]])@pos_runs)
    seqinfo(dr) <- seqinfo(rowRanges(res)[[1]])
    
    expect_identical(gr, dr)
    expect_equal(getTileNumber(res), 1)

    ## test_gr <- GRanges("chrI", IRanges(10000,11999))
    ## seqlengths(test_gr) <- seqlengths(ggdl)
    ## res <- subset(ggdl, seqnames == "chrI" & pos < 12000)
    ## dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res)[[1]])),
    ##                              IRanges::ranges(rowRanges(res)[[1]])@pos_runs)
    ## seqinfo(dr) <- seqinfo(rowRanges(res)[[1]])
    
    ## expect_identical(test_gr, dr)
    ## expect_equal(getTileNumber(res), 2)

    ## res <- ggdl[gr]
    ## dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res)[[1]])),
    ##                              IRanges::ranges(rowRanges(res)[[1]])@pos_runs)
    ## seqinfo(dr) <- seqinfo(rowRanges(res)[[1]])
    
    ## expect_identical(gr, dr)
    ## expect_equal(getTileNumber(res), 1)

    ## res <- ggdl[[1]]
    ## dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res)[[1]])),
    ##                              IRanges::ranges(rowRanges(res)[[1]])@pos_runs)
    ## seqinfo(dr) <- seqinfo(rowRanges(res)[[1]])
    
    ## expect_identical(granges(getIndex(ggdl)[1]), dr)
    ## expect_equal(getTileNumber(res), 1)

    ## gr <- GRanges("chrI", IRanges(7501,8500))
    ## res <- subsetByOverlaps(ggdl, gr)
    ## expect_true(all(dim(res) == c(0,4)))
    ## expect_true(all(colnames(res) == colnames(ggdl)))

    ## res <- subset(ggdl, seqnames == "chrI" & pos <= 9000)
    ## expect_true(all(dim(res) == c(0,4)))
    ## expect_true(all(colnames(res) == colnames(ggdl)))
})

test_that("Metric computation works correct in normal case", {
    ## make ggd and select random tile
    ggdl <- makeTestGenoGAMDataSetList()
    tile <- sample(1:5, 1)
    l <- extractList(assay(ggdl)[[1]], ranges(getIndex(ggdl)[tile]))

    ## for sum
    res <- sum(ggdl)
    test_res <- sapply(l[[1]], sum)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for mean
    res <- mean(ggdl)
    test_res <- sapply(l[[1]], mean)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for var
    res <- var(ggdl)
    test_res <- sapply(l[[1]], var)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for sd
    res <- sd(ggdl)
    test_res <- sapply(l[[1]], sd)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for median
    res <- median(ggdl)
    test_res <- sapply(l[[1]], median)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for mad
    res <- mad(ggdl)
    test_res <- sapply(l[[1]], mad)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))

    ## for IQR
    res <- IQR(ggdl)
    test_res <- sapply(l[[1]], IQR)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res[tile,] == test_res))
    
})

test_that("Metric computation works correct in case of one tile", {

    ## TODO: Continue here, correct .makeTile to not exceed given regions
    
    ## make ggd and select random tile
    ggdl <- test_ggdl
    getTileNumber(ggdl) <- 1
    
    ## for sum
    res <- sum(ggdl)
    test_res <- sapply(assay(ggdl)[[1]], sum)
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res == test_res))

    ## for mean
    res <- as.vector(mean(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], mean))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all.equal(res, test_res))

    ## for var
    res <- as.vector(var(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], var))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(identical(res, test_res))

    ## for sd
    res <- as.vector(sd(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], sd))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(identical(res, test_res))

    ## for median
    res <- as.vector(median(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], median))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(all(res == test_res))

    ## for mad
    res <- as.vector(mad(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], mad))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(identical(res, test_res))

    ## for IQR
    res <- as.vector(IQR(ggdl))
    test_res <- as.vector(sapply(assay(ggdl)[[1]], IQR))
    
    expect_true(all(dim(res) == c(length(getIndex(ggdl)), ncol(assay(ggdl)[[1]]))))
    expect_true(identical(res, test_res))
})


test_that("Metric computation works correct in case of empty GenoGAMDataSet", {
    ## make ggd and select random tile
    ggdl <- GenoGAMDataSetList()
    
    ## for sum
    res <- sum(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for mean
    res <- mean(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for var
    res <- var(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for sd
    res <- sd(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for median
    res <- median(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for mad
    res <- mad(ggdl)
    expect_true(all(dim(res) == c(1,1)))

    ## for IQR
    res <- IQR(ggdl)
    expect_true(all(dim(res) == c(1,1)))
})

## test .getCoordinates (with and without gaps) and .getChunkCoords here
test_that("Coordinate and Chunk transformation work correctly", {
    ## empty input
    emptyGGDL <- GenoGAMDataSetList()
    emptyCoords <- .getCoordinates(emptyGGDL)
    emptyChunks <- .getChunkCoords(emptyCoords)
    expect_true(length(emptyCoords) == 0)
    expect_true(length(emptyChunks) == 0)

    ## normal input without overhang
    ggdl <- makeTestGenoGAMDataSetList()
    getOverhangSize(ggdl) <- 0
    index <- getIndex(ggdl)
    coords <- .getCoordinates(ggdl)
    chunks <- .getChunkCoords(coords)
    
    expect_true(length(coords) == length(index))
    widthCoords <- max(end(coords)) - min(start(coords)) + 1
    expect_true(widthCoords == sum(width(dataRange(ggdl))))
    expect_true(all(width(coords) == width(index)))

    expect_true(length(chunks) == length(index))
    widthChunks <- max(end(chunks)) - min(start(chunks)) + 1
    expect_true(widthChunks == sum(width(dataRange(ggdl))))
    expect_true(all(width(chunks) == width(index)))

    ## normal input with overhang
    ggdl <- makeTestGenoGAMDataSetList()
    coords <- .getCoordinates(ggdl)
    index <- getIndex(ggdl)

    expect_error(.getCoordinates())
    expect_true(length(coords) == length(index))
    widthCoords <- max(end(coords)) - min(start(coords)) + 1
    expect_true(widthCoords == sum(width(dataRange(ggdl))))
    expect_true(all(width(coords) == width(index)))

    chunks <- .getChunkCoords(coords)
    ## all indices must be of same length
    expect_true(all.equal(length(coords), length(index), length(chunks)))
    
    ## The first and the second chunks of a chromosome must be of same length
    ## as the last and the second-to-last. 
    expect_true(width(chunks[1,]) == width(chunks[length(chunks),]))
    expect_true(width(chunks[2,]) == width(chunks[length(chunks) - 1,]))

    ## The rest should be of size = chunkSize
    numEqualChunks <- length(which(width(chunks) == getChunkSize(ggdl)))
    ## 3*4 because we have 4 chunks (the first and last two beeing different)
    ## and that across 3 chromosomes, so it repeats three times
    expect_true(numEqualChunks == (length(chunks) - 3*4))

    ## But on average they should always be exactly chunkSize
    expect_true(all.equal(mean(width(chunks)), getChunkSize(ggdl)))
})
