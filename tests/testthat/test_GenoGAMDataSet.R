context("Testing the GenoGAMDataSet-class")

## Unittests
gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
seqlengths(gr) <- 1e6
df <- DataFrame(colA = 1:10000, colB = round(runif(10000)))
se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
ggd <- GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                      design = ~ s(x))

test_that("The constructor works generally correctly", {
    df <- GenoGAMDataSet()
    expect_identical(new("GenoGAMDataSet"), df)

    expect_error(GenoGAMDataSet(10))
    expect_error(GenoGAMDataSet(se))
})

test_that("The constructor works correctly for SummarizedExperiment", {
    test_ggd <- ggd
    expect_identical(assay(test_ggd), assay(se))
    expect_identical(rowRanges(test_ggd), rowRanges(se))
    expect_is(test_ggd, "GenoGAMDataSet")

    gr <- GPos(GRanges("chr1", IRanges(1, 10000)))
    df <- DataFrame(colA = 1:10000, colB = round(runif(10000)))
    se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
    expect_error(GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                                design = ~ s(x) + s(x, by = experiment)),
                 "Sequence lengths missing in the Seqinfo object of SummarizedExperiment")

    gr <- GPos(GRanges(c("chr1", "chr1"), IRanges(c(1,5001), c(10000, 10000))))
    seqlengths(gr) <- 1e6
    df <- DataFrame(colA = 1:15000, colB = round(runif(15000)))
    se <- SummarizedExperiment(rowRanges = gr, assays = list(df))
    expect_error(GenoGAMDataSet(se, chunkSize = 2000, overhangSize = 250, 
                                design = ~ s(x) + s(x, by = experiment)),
                 "Overlapping regions encountered. Please reduce ranges and data first.")
})

test_that("The constructor works correctly for config files and data.frames", {
    config <- system.file("extdata/Set1", "experimentDesign.txt", package = "GenoGAM")
    dir <- system.file("extdata/Set1/bam", package = "GenoGAM")

    region <- GRanges("chrI", IRanges(10000, 20000))
    params <- Rsamtools::ScanBamParam(which = region)
    settings <- GenoGAMSettings(bamParams = params)
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings)
    expect_true(checkObject(ds))
    expect_equal(getTileNumber(ds), 10)
    expect_true(all(sapply(assay(ds), sum) > 0))
    

    df <- read.table(config, header = TRUE, sep = '\t')
    ds <- GenoGAMDataSet(df, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings)
    expect_true(checkObject(ds))
    expect_equal(getTileNumber(ds), 10)
    expect_true(all(sapply(assay(ds), sum) > 0))

    region = GRanges("chrV", IRanges(105000, 108000))
    params <- Rsamtools::ScanBamParam(which = region)
    settings <- GenoGAMSettings(bamParams = params)
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings)
    expect_false(checkObject(ds))
    expect_true(is.null(getTileNumber(ds)))
    expect_true(length(assays(ds)) == 0)

    settings <- GenoGAMSettings(chromosomeList = c("chrV"))
    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = genotype), directory = dir,
                         settings = settings)
    expect_false(checkObject(ds))
    expect_true(is.null(getTileNumber(ds)))
    expect_true(length(assays(ds)) == 0)

    ds <- GenoGAMDataSet(config, chunkSize = 1000, overhangSize = 200,
                         design = ~ s(x) + s(x, by = myColumn), directory = dir)
    expect_false(checkObject(ds))
})

test_that("Tiling works correctly under full chromosome conditions", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 200, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, length(tiles))
    expect_equal(metadata(tiles)$numTiles, 20)
    expect_true(max(width(tiles)) == 1400 & min(width(tiles)) == 1400)
})

test_that("Tiling works correctly with no overhang", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 0, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, length(tiles))
    expect_equal(metadata(tiles)$numTiles, 20)
    expect_true(max(width(tiles)) == 1000 & min(width(tiles)) == 1000)
})

test_that("Tiling works correctly with incorrect input", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(10000, 10000)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = -100, chromosomes = chromosomes)
    expect_error(.makeTiles(l), "Overhang size must be equal or greater than 0")

    l <- list(chunkSize = 1000, overhangSize = 0, chromosomes = GRanges())
    expect_error(.makeTiles(l), "Chromosome list should contain at least one entry")
})

test_that("Tiling works correctly with tiles bigger than chromosomes", {
    chromosomes <- GRanges(c("chr1", "chr2"), IRanges(c(1,1), c(6000, 10000)))
    seqlengths(chromosomes) <- c(6000, 1e6)
    l <- list(chunkSize = 8000, overhangSize = 100, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    expect_true(length(tiles) == 3)
    expect_true(all(width(tiles) == c(6000, 6000, 6000)))
})

test_that("Tiling works correctly with a distinct set of regions", {
    chromosomes <- GRanges(c("chr1", "chr1", "chr2", "chr2"),
                           IRanges(c(1, 1200, 1, 10000), c(2000, 10000, 1101, 11100)))
    seqlengths(chromosomes) <- c(1e6, 1e6)
    l <- list(chunkSize = 1000, overhangSize = 200, chromosomes = chromosomes)
    tiles <- .makeTiles(l)
    
    expect_equal(metadata(tiles)$numTiles, 16)
    expect_true(max(width(tiles)) == 1101 & min(width(tiles)) == 1101)
})

test_that("Settings checking functions work correctly", {
    test_ggd <- ggd

    expect_true(checkObject(test_ggd))
    
    metadata(slot(test_ggd, "index"))$check <- FALSE
    expect_false(checkObject(test_ggd))

    metadata(slot(test_ggd, "index"))$check <- NULL
    expect_false(checkObject(test_ggd))

    metadata(slot(test_ggd, "index"))$chunkSize <- 100
    expect_true(is(.checkChunkSize(test_ggd), "character"))

    metadata(slot(test_ggd, "index"))$tileSize <- 1000
    expect_true(is(.checkTileSize(test_ggd), "character"))

    end(slot(test_ggd, "index"))[1] <- 10000
    expect_true(is(.checkEqualityOfTiles(test_ggd), "character"))

    gr2 <- GRanges("chr2", IRanges(1, 1000), id = 30)
    seqlengths(gr2) <- 1e6
    md <- metadata(slot(test_ggd, "index"))
    suppressWarnings(slot(test_ggd, "index") <- c(slot(test_ggd, "index"), gr2))
    metadata(slot(test_ggd, "index")) <- md
    expect_true(is(.checkChromosomes(test_ggd), "character"))

    expect_true(is(.checkNumberOfTiles(test_ggd), "character"))

    metadata(slot(test_ggd, "index"))$chromosomes <- gr2
    expect_true(is(.checkTileRanges(test_ggd), "character"))
})

test_that("Accessors return the right slots", {
    test_ggd <- ggd
    dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(test_ggd))),
                                 IRanges::ranges(rowRanges(test_ggd))@pos_runs)
    seqinfo(dr) <- seqinfo(rowRanges(test_ggd))

    expect_identical(getIndex(test_ggd), slot(test_ggd, "index"))
    expect_identical(tileSettings(test_ggd), metadata(slot(test_ggd, "index")))
    expect_identical(dataRange(test_ggd), dr)
    expect_identical(getChromosomes(test_ggd), metadata(slot(test_ggd, "index"))$chromosomes)
    expect_identical(getTileSize(test_ggd), metadata(slot(test_ggd, "index"))$tileSize)
    expect_identical(getChunkSize(test_ggd), metadata(slot(test_ggd, "index"))$chunkSize)
    expect_identical(getOverhangSize(test_ggd), metadata(slot(test_ggd, "index"))$overhangSize)
    expect_identical(getTileNumber(test_ggd), metadata(slot(test_ggd, "index"))$numTiles)
    expect_identical(design(test_ggd), slot(test_ggd, "design"))
    expect_true(all(sizeFactors(test_ggd) == slot(test_ggd, "sizeFactors")))

    expect_error(design(test_ggd) <- ~ s(x) + s(by = myColumn))
    getChunkSize(test_ggd) <- 2500
    expect_equal(getTileNumber(test_ggd), 4)
    getTileSize(test_ggd) <- 2500
    expect_equal(getTileNumber(test_ggd), 5)
    getOverhangSize(test_ggd) <- 500
    expect_equal(median(width(getIndex(test_ggd))), 3000)
    getTileNumber(test_ggd) <- 2
    expect_equal(getChunkSize(test_ggd), 5000)
})

test_that("The read-in functions work correctly", {
    config <- system.file("extdata/Set1", "experimentDesign.txt", package = "GenoGAM")
    config <- read.table(config, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    for(ii in 1:nrow(config)) {
        absPath <- system.file("extdata/Set1/bam", config$file[ii], package = "GenoGAM")
        config$file[ii] <- absPath
    }
    expect_true(is(readData(config), "DataFrame"))

    region <- GRanges("chrI", IRanges(10000, 20000))
    params <- Rsamtools::ScanBamParam(which = region)
    reads <- GenomicAlignments::readGAlignments(config$file[1], param = params)

    gr <- .processCountChunks(reads, center = TRUE)
    expect_true(length(gr) > 0)
    
    gr <- .centerFragments(reads, asMates = FALSE)
    expect_true(length(gr) > 0)
    expect_true(median(width(gr)) == 1)

    gr <- .countFragments(reads, asMates = FALSE)
    expect_true(length(gr) > 0)
    expect_true(median(width(gr)) > 1)
})

test_that("The subsetting methods work correct", {
    ggd <- makeTestGenoGAMDataSet()
    getOverhangSize(ggd) <- 0
    gr <- GRanges("chrI", IRanges(16501,17500))
    seqlengths(gr) <- seqlengths(ggd)

    res <- subsetByOverlaps(ggd, gr)
    dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res))),
                                 IRanges::ranges(rowRanges(res))@pos_runs)
    seqinfo(dr) <- seqinfo(rowRanges(res))
    
    expect_identical(gr, dr)
    expect_equal(getTileNumber(res), 1)

    test_gr <- GRanges("chrI", IRanges(10000,19999))
    seqlengths(test_gr) <- seqlengths(ggd)
    res <- subset(ggd, seqnames == "chrI" & pos < 20000)
    dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res))),
                                 IRanges::ranges(rowRanges(res))@pos_runs)
    seqinfo(dr) <- seqinfo(rowRanges(res))
    
    expect_identical(test_gr, dr)
    expect_equal(getTileNumber(res), 10)

    res <- ggd[gr]
    dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res))),
                                 IRanges::ranges(rowRanges(res))@pos_runs)
    seqinfo(dr) <- seqinfo(rowRanges(res))
    
    expect_identical(gr, dr)
    expect_equal(getTileNumber(res), 1)

    ## res <- ggd[[1]]
    ## dr <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(rowRanges(res))),
    ##                              IRanges::ranges(rowRanges(res))@pos_runs)
    ## seqinfo(dr) <- seqinfo(rowRanges(res))
    
    ## expect_identical(granges(getIndex(ggd)[1]), dr)
    ## expect_equal(getTileNumber(res), 1)

    gr <- GRanges("chrI", IRanges(7501,8500))
    res <- subsetByOverlaps(ggd, gr)
    expect_true(all(dim(res) == c(0,4)))
    expect_true(all(colnames(res) == colnames(ggd)))

    res <- subset(ggd, seqnames == "chrXIV" & pos <= 9000)
    expect_true(all(dim(res) == c(0,4)))
    expect_true(all(colnames(res) == colnames(ggd)))
})

test_that("Metric computation works correct in normal case", {
    ## make ggd and select random tile
    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    tile <- sample(length(getIndex(ggd)), 1)
    l <- extractList(assay(ggd), ranges(getIndex(ggd)[tile]))

    ## for sum
    res <- sum(ggd)
    test_res <- sapply(l[[1]], sum)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for mean
    res <- mean(ggd)
    test_res <- sapply(l[[1]], mean)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for var
    res <- var(ggd)
    test_res <- sapply(l[[1]], var)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for sd
    res <- sd(ggd)
    test_res <- sapply(l[[1]], sd)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for median
    res <- median(ggd)
    test_res <- sapply(l[[1]], median)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for mad
    res <- mad(ggd)
    test_res <- sapply(l[[1]], mad)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))

    ## for IQR
    res <- IQR(ggd)
    test_res <- sapply(l[[1]], IQR)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res[tile,] == test_res))
    
})

test_that("Metric computation works correct in case of one tile", {
    ## make ggd and select random tile
    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    getTileNumber(ggd) <- 1
    
    ## for sum
    res <- sum(ggd)
    test_res <- sapply(assay(ggd), sum)
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all(res == test_res))

    ## for mean
    res <- as.vector(mean(ggd))
    test_res <- as.vector(sapply(assay(ggd), mean))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(all.equal(res, test_res))

    ## for var
    res <- as.vector(var(ggd))
    test_res <- as.vector(sapply(assay(ggd), var))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(identical(res, test_res))

    ## for sd
    res <- as.vector(sd(ggd))
    test_res <- as.vector(sapply(assay(ggd), sd))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(identical(res, test_res))

    ## for median
    res <- as.vector(median(ggd))
    test_res <- as.vector(sapply(assay(ggd), median))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(identical(res, test_res))

    ## for mad
    res <- as.vector(mad(ggd))
    test_res <- as.vector(sapply(assay(ggd), mad))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(identical(res, test_res))

    ## for IQR
    res <- as.vector(IQR(ggd))
    test_res <- as.vector(sapply(assay(ggd), IQR))
    
    expect_true(all(dim(res) == c(length(getIndex(ggd)), ncol(assay(ggd)))))
    expect_true(identical(res, test_res))
})


test_that("Metric computation works correct in case of empty GenoGAMDataSet", {
    ## make ggd and select random tile
    ggd <- GenoGAMDataSet()
    
    ## for sum
    res <- sum(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for mean
    res <- mean(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for var
    res <- var(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for sd
    res <- sd(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for median
    res <- median(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for mad
    res <- mad(ggd)
    expect_true(all(dim(res) == c(1,1)))

    ## for IQR
    res <- IQR(ggd)
    expect_true(all(dim(res) == c(1,1)))
})

## test .getCoordinates (with and without gaps) and .getChunkCoords here
test_that("Coordinate and Chunk transformation work correctly", {
    ## empty input
    emptyGGD <- GenoGAMDataSet()
    emptyCoords <- .getCoordinates(emptyGGD)
    emptyChunks <- .getChunkCoords(emptyCoords)
    expect_true(length(emptyCoords) == 0)
    expect_true(length(emptyChunks) == 0)

    ## normal input without overhang
    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    getOverhangSize(ggd) <- 0
    index <- getIndex(ggd)
    coords <- .getCoordinates(ggd)
    chunks <- .getChunkCoords(coords)
    
    expect_true(length(coords) == length(index))
    widthCoords <- max(end(coords)) - min(start(coords)) + 1
    expect_true(widthCoords == length(rowRanges(ggd)))
    expect_true(all(width(coords) == width(index)))

    expect_true(length(chunks) == length(index))
    widthChunks <- max(end(chunks)) - min(start(chunks)) + 1
    expect_true(widthChunks == length(rowRanges(ggd)))
    expect_true(all(width(chunks) == width(index)))

    ## normal input with overhang
    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    coords <- .getCoordinates(ggd)
    index <- getIndex(ggd)

    expect_error(.getCoordinates())
    expect_true(length(coords) == length(index))
    widthCoords <- max(end(coords)) - min(start(coords)) + 1
    expect_true(widthCoords == length(rowRanges(ggd)))
    expect_true(all(width(coords) == width(index)))

    chunks <- .getChunkCoords(coords)
    ## all indices must be of same length
    expect_true(all.equal(length(coords), length(index), length(chunks)))
    
    ## The first and the second chunks of a chromosome must be of same length
    ## as the last and the second-to-last. 
    expect_true(width(chunks[1,]) == width(chunks[length(chunks),]))
    expect_true(width(chunks[2,]) == width(chunks[length(chunks) - 1,]))

    ## The rest should be of size = chunkSize
    numEqualChunks <- length(which(width(chunks) == getChunkSize(ggd)))
    expect_true(numEqualChunks == (length(chunks) - 4))

    ## But on average they should always be exactly chunkSize
    expect_true(all.equal(mean(width(chunks)), getChunkSize(ggd)))
})

