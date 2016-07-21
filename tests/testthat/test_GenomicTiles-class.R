context("Testing GenomicTiles-class")

test_that("The accessor functions work correctly", {
    gt <- makeTestGenomicTiles()
    expect_equal(getIndex(gt), gt@index)
    expect_true(all(width(getIndex(gt)) > width(getChunkIndex(gt))))
    expect_true(all(getChunkIndex(gt) == getIndex(untile(gt))))
    expect_is(tileSettings(gt), "list")
    expect_equal(length(tileSettings(gt)), 6)
    expect_equal(getCoordinates(gt), gt@coordinates)
    expect_equal(width(getIndexCoordinates(gt)), width(getIndex(gt)))
    expect_true(sum(getIndexCoordinates(gt) != getIndex(gt)) != length(getIndex(gt)))
    expect_true(all(width(getCoordinates(gt)), width(dataRange(gt))))
    expect_equal(getChromosomes(gt), tileSettings(gt)$chromosomes)
    expect_equal(getTileSize(gt), tileSettings(gt)$tileSize)
    expect_equal(getChunkSize(gt), tileSettings(gt)$chunkSize)
    expect_equal(getOverhangSize(gt), tileSettings(gt)$overhangSize)
    expect_equal(getTileNumber(gt), tileSettings(gt)$numTiles)
})

test_that("Setting changes work correctly", {
    gt <- makeTestGenomicTiles()
    test <- changeSettings(gt, "chunkSize", 10)
    expect_equal(round(mean(width(getChunkIndex(test)))), 10)
    test <- changeSettings(gt, "tileSize", 20)
    expect_true(all(width(getIndex(test)) == 20))
    test <- changeSettings(gt, "overhangSize", 0)
    expect_true(all(getIndex(test) == getChunkIndex(test)))
    test <- changeSettings(gt, "chromosomes", "chrI")
    expect_true(getChromosomes(test) == getChromosomes(gt)[1,])
    test <- changeSettings(gt, "numTiles", 4)
    expect_equal(length(getIndex(test)), 4)
})

test_that("Subsetting works correctly", {
    gt <- makeTestGenomicTiles()
    gr <- GRanges("chrI", IRanges(21, 40))
    sub <- subset(gt, pos <= 20)
    ov <- subsetByOverlaps(gt, gr)

    expect_equal(length(rowRanges(sub)), 40)
    expect_true(all(c("chrI", "chrII") %in% seqlevels(rowRanges(sub))))
    expect_equal(length(rowRanges(ov)), 37)
})

test_that("Coercion works correctly", {
    gt <- makeTestGenomicTiles()
    gr <- GRanges("chrI", IRanges(21, 40))
    sub <- subset(gt, seqnames == "chrI" & pos <= 40 & pos >= 21)
    test <- as.data.frame(gt)
    expect_is(test, "data.frame")
    test <- getTile(gt, 2)[[1]]
    gptest <- GPos(GRanges(test$seqnames[1], IRanges(min(test$pos), max(test$pos))))
    gpreal <- GPos(getIndex(gt)[2,])
    seqlevels(gpreal, force = TRUE) <- seqlevelsInUse(gpreal)
    expect_is(test, "DataFrame")
    expect_true(all(gptest == gpreal))
    expect_equal(test, gt[[2]][[1]])
    expect_true(all(getIndex(sub) == getIndex(gt[gr])))
})

test_that("Metrics computation works fine", {
    gt <- makeTestGenomicTiles()
    sgt <- sum(gt)[[1]]
    mgt <- mean(gt)[[1]]
    vgt <- var(gt)[[1]]
    sdgt <- sd(gt)[[1]]
    medgt <- median(gt)[[1]]
    madgt <- mad(gt)[[1]]
    expect_true(all.equal(dim(sgt), dim(mgt), dim(vgt), dim(sdgt), dim(medgt), dim(madgt)))
})


