### ===========================================================
### Genomic Tiles class
### ===========================================================

#' GenomicTiles class
#'
#' This class is designed to represent the entire genome (or a subset of
#' it) and any additional data associated with the samples or positions.
#' It extends the RangedSummarizedExperiment class and adds two additional
#' index slots to keep track of the data. The main change compared to
#' RangedSummarizedExperiment is the use of a GPos (basepair level) instead of
#' GRanges (ranges level) object as rowRanges and the use of two GRanges objects
#' as indices. The GPos object allows to store raw instead of summarized data
#' in the assays. Because of this the size of genomic data can increase
#' tremendously. Thus the GenomicTiles class automatically divides the data
#' in (overlapping) tiles, making any operation on this data easy executable
#' in parallel.
#' @slot index A GRanges object that stores the tiles ranges and their index
#' in the genome space. That is, ranges are the positions on the genome.
#' @slot coordinates A GRanges object that stores the tiles ranges and their index
#' in the DataFrame space. That is ranges are the row positions in the DataFrame.
#' @details For all other slots see \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @rdname GenomicTiles-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenomicTiles
setClass("GenomicTiles",
         contains = "RangedSummarizedExperiment",
         slots = list(index = "GRanges", coordinates = "GRanges"),
         prototype = list(index = GenomicRanges::GRanges(), coordinates = GenomicRanges::GRanges()))

## Validity
## ========

#' Validating the correct type
.validateIndexType <- function(object) {
    if(class(object@index) != "GRanges") {
        return("'index' must be a GRanges object")
    }
    NULL
}

.validateCoordinatesType <- function(object) {
    if(class(object@coordinates) != "GRanges") {
        return("'coordinates' must be a GRanges object")
    }
    NULL
}

.validateChromosomes <- function(object) {
    cindex <- GenomeInfoDb::seqlevels(object@index)
    ccoords <- GenomeInfoDb::seqlevels(object@coordinates)
    cobject <- GenomeInfoDb::seqlevels(rowRanges(object))
    if(!(all(cindex %in% ccoords) && all(cobject %in% ccoords))) {
        return("Different chromosomes for data and index objects.")
    }
    NULL
}

## general validate function
.validateGenomicTiles <- function(object) {
    c(.validateIndexType(object), .validateCoordinatesType(object),
      .validateChromosomes(object))
}

setValidity2("GenomicTiles", .validateGenomicTiles)

## Constructor
## ===========

#' GenomicTiles constructor.
#'
#' This is the constructor function for GenomicTiles. The easiest
#' construction is fropm SummarizedExperiment. However as the
#' class operates on basepair level, the rowRanges are restricted
#' to the GPos class.
#'
#' @param assays One of two things. Either directly an object of
#' type 'RangedSummarizedExperiment'. Or in case the object is
#' created from raw data, a 'list' or 'SimpleList' of matrix-like
#' elements, or a matrix-like object. All elements of the list must
#' have the same dimensions, and dimension names (if present) must
#' be consistent across elements and with the row names of 'rowRanges'
#' and 'colData'.
#' @param chunkSize An integer specifying the size of one chunk in bp.
#' @param overhangSize An integer specifying the size of the overhang in bp.
#' The overhang is regarded to be symmetric, such that only the overhang
#' of one side should be provided. 
#' @param ... Further parameters passed to the
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} constructor.
#' @return An object of class GenomicTiles.
#' @details Most, but not necessary all functionalities of
#' SummarizedExperiment are yet provided.
#' @examples
#' ## from raw data
#' gp <- GPos(GRanges(c("chrI", "chrII"), IRanges(c(1,1), c(5,5))))
#' assay <- matrix(1:10, 10, 1)
#' gt <- GenomicTiles(assay, chunkSize = 3, rowRanges = gp)
#'
#' ## from SummarizedExperiment
#' se <- SummarizedExperiment(assay, rowRanges = gp)
#' gt <- GenomicTiles(se, chunkSize = 3)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
GenomicTiles <- function(assays, chunkSize = 1e4, overhangSize = 0, ...) {

    if(missing(assays)) {
        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = GenomicRanges::GPos(GenomicRanges::GRanges()))
        return(new("GenomicTiles", se))
    }

    if(class(assays) == "RangedSummarizedExperiment") {
        gt <- .GenomicTilesFromSE(assays, chunkSize, overhangSize)
    }
    else {
        gt <- .GenomicTilesFromRawData(assays, chunkSize, overhangSize, ...)
    }
    
    return(gt)
}

#' Creating GenomicTiles from SummarizedEperiment
#'
#' @param assays A RangedSummarizedExperiment object
#' @param chunkSize An integer specifying the size of one chunk in bp.
#' Two out of tileSize, chunkSize and overhangSize must be specified.
#' @param overhangSize An integer specifying the size of the overhang in bp.
#' The overhang is regarded to be symmetric, such that only the overhang
#' of one side should be provided. Two out of tileSize, chunkSize and
#' overhang must be specified.
#' @return An object of class GenomicTiles.
.GenomicTilesFromSE <- function(assays, chunkSize, overhangSize) {
    
    chroms <- GenomeInfoDb::seqlevelsInUse(assays)
    gr <- rowRanges(assays)@pos_runs

    tileSize <- chunkSize + 2*overhangSize
    l <- list(chromosomes = gr,
              chunkSize = chunkSize,
              tileSize = tileSize,
              overhangSize = overhangSize,
              chunk = FALSE)
    tiles <- .makeTiles(l)
    coords <- .makeCoordinates(rowRanges(assays))
    gt <- new("GenomicTiles", index = tiles,
              coordinates = coords, assays)
    return(gt)
}

#' Creating GenomicTiles from SummarizedEperiment
#'
#' @param assays A 'list' or 'SimpleList' of matrix-like elements, or a
#' matrix-like object. All elements of the list must have the same
#' dimensions, and dimension names (if present) must be consistent
#' across elements and with the row names of 'rowRanges' and 'colData'.
#' @param chunkSize An integer specifying the size of one chunk in bp.
#' Two out of tileSize, chunkSize and overhangSize must be specified.
#' @param overhangSize An integer specifying the size of the overhang in bp.
#' The overhang is regarded to be symmetric, such that only the overhang
#' of one side should be provided. Two out of tileSize, chunkSize and
#' overhang must be specified.
#' @param ... Further parameters passed to the SummarizedExperiment
#' constructor.
#' @return An object of class GenomicTiles.
.GenomicTilesFromRawData <- function(assays, chunkSize, overhangSize, ...) {
    rowRanges <- list(...)$rowRanges
    gr <- rowRanges@pos_runs

    tileSize <- chunkSize + 2*overhangSize
    l <- list(chromosomes = gr,
              chunkSize = chunkSize,
              tileSize = tileSize,
              overhangSize = overhangSize,
              chunk = FALSE)
    tiles <- .makeTiles(l)
    coords <- .makeCoordinates(rowRanges)

    se <- SummarizedExperiment(assays = list(assays), ...)
    gt <- new("GenomicTiles", index = tiles,
              coordinates = coords, se)
    return(gt)
}


#' A function to produce a GRanges index from a list of settings.
#'
#' @param l A list of settings.
#' @return A Granges object of the tiles.
.makeTiles <- function(l) {

    if(length(l) == 0) return(GenomicRanges::GRanges())

    duplicates <- FALSE

    ## make chunks
    chroms <- width(l$chromosomes)
    names(chroms) <- GenomeInfoDb::seqnames(l$chromosomes)

    ## if regionsa re smaller than requested tileSize
    if(any(l$tileSize > chroms)) {
        shrink <- min(chroms)/l$tileSize
        l$tileSize <- round(l$tileSize * shrink)
        l$chunkSize <- round(l$chunkSize * shrink)
        l$overhangSize <- round(l$overhangSize * shrink)
    }

    if(any(duplicated(names(chroms)))) {
        backupChroms <- l$chromosomes
        names(chroms) <- paste(names(chroms), 1:length(chroms), sep = ".")
        l$chromosomes <- GenomicRanges::GRanges(names(chroms), ranges(l$chromosomes))
        duplicates <- TRUE
    }
    
    chunks <- unlist(GenomicRanges::tileGenome(chroms, tilewidth = l$chunkSize,
                        cut.last.tile.in.chrom = FALSE))
    chunks <- suppressWarnings(IRanges::shift(chunks, width(chunks)/2))

    ## add overhang and expand to tiles
    tiles <- suppressWarnings(trim(flank(chunks, round(l$tileSize/2),
                                         both = TRUE)))

    ## resize start and end tiles till correct tile size is reached
    startsToResize <- which(width(tiles) < l$tileSize & start(tiles) == 1)
    tiles[startsToResize] <- resize(tiles[startsToResize],
                                    width = l$tileSize)
    endsToResize <- which(width(tiles) < l$tileSize)
    tiles[endsToResize] <- resize(tiles[endsToResize],
                                  width = l$tileSize, fix = "end")

    tiles <- tiles[!duplicated(tiles)]

    ## if there was a shift such that seqlengths got deleted
    if(any(is.na(c(GenomeInfoDb::seqlengths(tiles), GenomeInfoDb::seqlengths(l$chromosomes)))) |
       any(GenomeInfoDb::seqlengths(tiles) != GenomeInfoDb::seqlengths(l$chromosomes))) {
        GenomeInfoDb::seqlengths(tiles) <- GenomeInfoDb::seqlengths(l$chromosomes)
        regions <- splitAsList(l$chromosomes, GenomeInfoDb::seqlevels(l$chromosomes))
        tiles <- do.call(c, lapply(GenomeInfoDb::seqlevels(tiles), function(y) {
            IRanges::shift(tiles[GenomeInfoDb::seqnames(tiles) == y], start(regions[[y]]) - 1)
                          }))
    }

    if(duplicates) {
        temp <- strsplit(as.character(GenomeInfoDb::seqnames(tiles)), split = "\\.")
        newSeqnames <- sapply(temp, function(y) y[1])
        tiles <- GenomicRanges::GRanges(newSeqnames, ranges(tiles))
        l$chromosomes <- backupChroms
    }
    GenomeInfoDb::seqlengths(tiles) <- GenomeInfoDb::seqlengths(l$chromosomes)

    ## add 'id' and 'dist' column and put settings in metadata
    mcols(tiles)$id <- 1:length(tiles)
    l$numTiles <- length(tiles)
    metadata(tiles) <- l
    return(tiles)
}

#' A function to produce the row coordinates from GPos object.
#' 
#' @param gp A GPos object.
#' @return A GRanges object of the row coordinates.
.makeCoordinates <- function(gp) {
    gpCoords <- gp@pos_runs
    gpDims <- length(gpCoords)
    if(length(gpCoords) > 0) {
        ir <- IRanges(start = cumsum(c(1, end(gpCoords)[-gpDims[1]])),
                      end = cumsum(width(gpCoords)))
    }
    else {
        ir <- IRanges()
    }
    coords <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(gpCoords), ir)
    return(coords)
}

#' Make an example /code{GenomicTile}
#'
#' @return A /code{GenomicTiles} object
#' @examples
#' test <- makeTestGenomicTiles()
#' @export
makeTestGenomicTiles <- function() {
    gp <- GenomicRanges::GPos(GenomicRanges::GRanges(c("chrI", "chrII"), IRanges::IRanges(c(1,1), c(50,50))))
    df <- DataFrame(a = 1:100, b = 101:200)
    gt2 <- GenomicTiles(assays = list(df), chunkSize = 15, rowRanges = gp)
    return(gt2)
}


## Accessors
## =========

#' @rdname getIndex
#' @export
setGeneric("getIndex", function(object, ...) standardGeneric("getIndex"))

#' Accessor to the 'index' slot
#'
#' The index holds the Granges object that splits the entire dataset in
#' tiles.
#'
#' @rdname getIndex
#' @param object A /code{GenomicTiles} object.
#' @param id A vector if tile ids. By default the complete index is returned.
#' @param ... Additional arguments
#' @return A /code{GRanges} object representing the index
#' @examples 
#' gt <- makeTestGenomicTiles()
#' getIndex(gt)
#' getIndex(gt, 1:3)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("getIndex", signature(object = "GenomicTiles"), function(object, id = NULL) {
    idx <- slot(object, "index")
    if(!is.null(id)) idx[mcols(idx)$id %in% id]
    return(idx)
})

#' @rdname getChunkIndex
#' @export
setGeneric("getChunkIndex", function(object, ...) standardGeneric("getChunkIndex"))

#' Compute the index for chunks instead tiles
#'
#' The chunk index holds the Granges object that splits the entire dataset in
#' chunk, that is non-overlapping intervals.
#'
#' @rdname getChunkIndex
#' @param object A /code{GenomicTiles} object.
#' @param id A vector if tile ids. By default the complete index is returned.
#' @param ... Additional arguments
#' @return A /code{GRanges} object representing the index
#' @examples 
#' gt <- makeTestGenomicTiles()
#' getChunkIndex(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("getChunkIndex", "GenomicTiles", function(object, id = NULL) {
    indx <- getIndex(object)
    if(length(indx) == 0) return(indx)
    
    splitIndx <- split(indx, GenomeInfoDb::seqnames(indx))
    
    chunkIndex <- lapply(splitIndx, function(y) {
        start <- c(start(y[1]), (end(y[-length(y)]) + start(y[-1]))/2)
        end <- c((end(y[-length(y)]) + start(y[-1]))/2 - 1, end(y[length(y)]))
        start(y) <- start
        end(y) <- end
        return(y)
    })

    res <- do.call(c, unname(chunkIndex))
    metadata(res) <- metadata(indx)

    if(!is.null(id)) res <- res[mcols(res)$id %in% id]
    return(res)
})

#' @rdname untile
#' @export
setGeneric("untile", function(object, ...) standardGeneric("untile"))

#' Set index to chunkIndex
#'
#' Replace the tile index with the chunk index in /code{GenomicTiles} object
#'
#' @rdname untile
#' @param object A /code{GenomicTiles} object.
#' @param id A vector if tile ids. By default the complete index is taken.
#' @param ... Additional arguments
#' @return A modified /code{GenomicTiles} object
#' @examples 
#' gt <- makeTestGenomicTiles()
#' newGT <- untile(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("untile", "GenomicTiles", function(object, id = NULL) {
    slot(object, "index") <- getChunkIndex(object, id)
    if(length(metadata(slot(object, "index"))) != 0) {
        metadata(slot(object, "index"))$chunk <- TRUE
    }
    return(object)
})

#' @rdname tileSettings
#' @export
setGeneric("tileSettings", function(object) standardGeneric("tileSettings"))

#' Return tile settings
#'
#' Returns a list settings used to generate the tile index
#'
#' @rdname tileSettings
#' @param object A /code{GenomicTiles} object.
#' @return A list of tile settings
#' @examples 
#' gt <- makeTestGenomicTiles()
#' tileSettings(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("tileSettings", "GenomicTiles", function(object) {
    metadata(getIndex(object))
})

#' @rdname getCoordinates
#' @export
setGeneric("getCoordinates", function(object) standardGeneric("getCoordinates")) 

#' Accessor to the /code{coordinates} slot
#'
#' The /code{coordinates} slot contains the row coordinates of each
#' chromosome in the data. Such that taken a genomic position from
#' a chromosome it's easy to detect the correct row in the assay
#'
#' @rdname getCoordinates
#' @param object A /code{GenomicTiles} object.
#' @return A /code{GRanges} object of row coordinates
#' @examples 
#' gt <- makeTestGenomicTiles()
#' getCoordinates(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("getCoordinates", "GenomicTiles", function(object) {
    slot(object, "coordinates")
})

#' @rdname getIndexCoordinates
#' @export
setGeneric("getIndexCoordinates", function(object, ...)
           standardGeneric("getIndexCoordinates")) 

#' Compute the row coordinates for a given index
#'
#' Given an index of genomic positions, this method computes
#' the corresponding row positions in the assay
#'
#' @rdname getIndexCoordinates
#' @param object A /code{GenomicTiles} object.
#' @param id A vector if tile ids. By default the complete index is returned.
#' @param index A /code{Granges} object representing an index of genomic positions.
#' @param ... Additional arguments
#' Usually the original index or the chunk index.
#' @return A /code{GRanges} object of row coordinates
#' @examples 
#' gt <- makeTestGenomicTiles()
#' getIndexCoordinates(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("getIndexCoordinates", "GenomicTiles", function(object, id = NULL, index = NULL) {
    if(is.null(index)) {
        index <- getIndex(object)
    }
    chromosomes <- GenomeInfoDb::seqlevels(index)
    rowCoords <- getCoordinates(object)
    GenomeInfoDb::seqlengths(index) <- rep(NA, length(GenomeInfoDb::seqlengths(index)))
    
    for(chr in chromosomes) {
        start <- min(start(index[GenomeInfoDb::seqnames(index) == chr]))
        index[GenomeInfoDb::seqnames(index) == chr] <- shift(index[GenomeInfoDb::seqnames(index) == chr],
                         start(rowCoords[GenomeInfoDb::seqnames(rowCoords) == chr]) - start)
    }

    if(!is.null(id)) index <- index[mcols(index)$id %in% id]
    return(index)
})

#' @rdname dataRange
#' @export
setGeneric("dataRange", function(object) standardGeneric("dataRange")) 

#' The /code{GRanges} of the underlying data
#'
#' Just like the /code{coordinates} slot but returns the genomic ranges
#' of the underlying data.
#'
#' @rdname dataRange
#' @param object A /code{GenomicTiles} object.
#' @return A /code{GRanges} object of genomic ranges of
#' the underlying data
#' @examples 
#' gt <- makeTestGenomicTiles()
#' dataRange(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("dataRange", "GenomicTiles", function(object) {
    rowRanges(object)@pos_runs
})

#' @rdname dataRange
setMethod("dataRange", "GPos", function(object) {
    ranges(object@pos_runs)
})

#' @rdname tileSettings-elements
#' @export
setGeneric("getChromosomes", function(object) standardGeneric("getChromosomes")) 

#' The single entries of the tile settings
#'
#' Returns the single elements of the tile settings
#'
#' @rdname tileSettings-elements
#' @param object A /code{GenomicTiles} object.
#' @param ... Additional arguments
#' @return An integer value, or in case of /code{getChromosomes}
#' a /code{GRanges} object
#' @examples 
#' gt <- makeTestGenomicTiles()
#' getChromosomes(gt)
#' getTileSize(gt)
#' getChunkSize(gt)
#' getOverhangSize(gt)
#' getTileNumber(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("getChromosomes", "GenomicTiles", function(object) {
    metadata(getIndex(object))$chromosomes
})

#' @rdname tileSettings-elements
#' @export
setGeneric("getTileSize", function(object) standardGeneric("getTileSize")) 

#' @rdname tileSettings-elements
#' @export
setMethod("getTileSize", "GenomicTiles", function(object) {
    metadata(getIndex(object))$tileSize
})

#' @rdname tileSettings-elements
#' @export
setGeneric("getChunkSize", function(object) standardGeneric("getChunkSize")) 

#' @rdname tileSettings-elements
#' @export
setMethod("getChunkSize", "GenomicTiles", function(object) {
    metadata(getIndex(object))$chunkSize
})

#' @rdname tileSettings-elements
#' @export
setGeneric("getOverhangSize", function(object) standardGeneric("getOverhangSize")) 

#' @rdname tileSettings-elements
#' @export
setMethod("getOverhangSize", "GenomicTiles", function(object) {
    metadata(getIndex(object))$overhangSize
})

#' @rdname tileSettings-elements
#' @export
setGeneric("getTileNumber", function(object, ...) standardGeneric("getTileNumber")) 

#' @rdname tileSettings-elements
#' @export
setMethod("getTileNumber", "GenomicTiles", function(object) {
    metadata(getIndex(object))$numTiles
})

## Index Modification
## ===============

## The underlying function to unlist indeces, see the genric unlistCoordinates function
.unlistCoordinates <- function(object, index = NULL, chromosomes = NULL) {
    if(is.null(index)) index <- getCoordinates(object)
    if(is.null(chromosomes)) chromosomes <- GenomeInfoDb::seqlevels(index)
    
    numAssays <- length(assays(object))
    dims <- dim(object)

    res <- lapply(chromosomes, function(y) {
        idx <- index[GenomeInfoDb::seqnames(index) == y]
        unlistedCoords <- do.call(c, lapply(1:numAssays, function(step) {
            if(step > 1) shift(ranges(idx), dims[1])
            else ranges(idx)
        }))
        mcols(unlistedCoords)$id <- Rle(factor(rep(idx$id, numAssays)))
        return(unlistedCoords)
    })
    names(res) <- chromosomes
    return(res)
}

setGeneric("unlistCoordinates", function(object, ...) standardGeneric("unlistCoordinates")) 

setMethod("unlistCoordinates", "GenomicTiles", function(object, index = NULL,
                                                        chromosomes = NULL) {
    .unlistCoordinates(object, index = index, chromosomes = chromosomes)
})

setGeneric("unlistIndexCoordinates",
           function(object, ...) standardGeneric("unlistIndexCoordinates")) 

setMethod("unlistIndexCoordinates", "GenomicTiles", function(object, chromosomes = NULL) {
    rindex <- getIndexCoordinates(object)
    .unlistCoordinates(object, chromosomes = chromosomes, index = rindex)
})

setGeneric("unlistIndex", function(object, ...) standardGeneric("unlistIndex")) 

setMethod("unlistIndex", "GenomicTiles", function(object, chromosomes = NULL) {
    tindex <- getIndex(object)
    .unlistCoordinates(object, chromosomes = chromosomes, index = tindex)
})

## check a specified setting
.checkSettings <- function(object, params = c("chunkSize", "tileSize",
                                       "overhangSize",
                                       "chromosomes", "numTiles")) {
    param <- match.arg(params)
    switch(param,
           chunkSize = .checkChunkSize(object),
           tileSize = .checkTileSize(object),
           overhangSize = .checkOverhangSize(object),
           chromosomes = .checkChromosomes(object),
           numTiles = .checkNumberOfTiles(object))
}

## check the chunk size and return a logical value
.checkChunkSize <- function(object) {
    widths <- width(getIndex(object))
    diffs <- (widths - 2*getOverhangSize(object)) - getChunkSize(object)
    res <- all.equal(widths, rep(0, length(diffs)), tolerance = 1)
    return(res)
}

## check the tile size and return a logical value
.checkTileSize <- function(object) {
    widths <- width(getIndex(object))
    diffs <- widths - getTileSize(object)
    res <- all.equal(diffs, rep(0, length(diffs)), tolerance = 1)
    return(res)
}

## check the overhang size and return a logical value
.checkOverhangSize <- function(object) {
    widths <- width(getIndex(object))
    diffs <- (widths - getChunkSize(object))/2 - getOverhangSize(object)
    res <- all.equal(diffs, rep(0, length(diffs)), tolerance = 1)
    return(res)
}

## check the chromosome list and return a logical value
.checkChromosomes <- function(object) {
    objChroms <- sort(GenomeInfoDb::seqlengths(object))
    indexChroms <- sort(GenomeInfoDb::seqlengths(getIndex(object)))
    validChroms <- sort(GenomeInfoDb::seqlengths(getChromosomes(object)))
    res1 <- all.equal(indexChroms, validChroms, objChroms)
    res2 <- all.equal(names(indexChroms), names(validChroms), names(objChroms))
    return(res1 & res2)
}

## check the number of tiles and return a logical value
.checkNumberOfTiles <- function(object) {
    tiles <- length(getIndex(object))
    validTiles <- getTileNumber(object)
    res <- tiles == validTiles
    return(res)
}

.checkGenomicTiles <- function(object) {
    futile.logger::flog.info("Check if tile settings match the data.")

    params = c("chunkSize", "tileSize",
        "overhangSize", "chromosomes", "numTiles")

    settings <- tileSettings(object)
    
    if(is.null(settings$check)) {
        futile.logger::flog.warn(paste("Checks dismissed due to empty object"))
        return(FALSE)
    }
    if(!settings$check) {
        futile.logger::flog.warn("Settings checking dismissed, probably due to different tile size. Computation on these tiles is discouraged.")
        return(FALSE)
    }
    res <- sapply(params, .checkSettings, object = object)
    if(all(res)) {
        futile.logger::flog.info("All checks passed.")
        return(TRUE)
    }
    else {
        futile.logger::flog.error(paste("Checks failed in",
                         paste(params[!res], collapse = ", ")))
        return(FALSE)
    }
}

#' @rdname checkSettings
#' @export
setGeneric("checkSettings", function(object) standardGeneric("checkSettings")) 

#' Check data compliance with tile settings
#'
#' Check if the indices were build correctly, according to the
#' specified parameters
#'
#' @rdname checkSettings
#' @param object A /code{GenomicTiles} object.
#' @return A logical value
#' @examples 
#' gt <- makeTestGenomicTiles()
#' checkSettings(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("checkSettings", "GenomicTiles", function(object) {
    .checkGenomicTiles(object)
})

.changeSettings <- function(object, params = c("chunkSize", "tileSize",
                                        "overhangSize", "chromosomes",
                                        "numTiles"), value) {

    ## writable settings, exclude config, regions and unit
    param <- match.arg(params)
    switch(param,
           chunkSize = .changeChunkSize(object, value),
           tileSize = .changeTileSize(object, value),
           overhangSize = .changeOverhangSize(object, value),
           chromosomes = .changeChromosomes(object, value),
           numTiles = .changeNumberOfTiles(object, value))
}

## change the chunk size and return the modified object
.changeChunkSize <- function(object, value) {
    settings <- tileSettings(object)
    settings$chunkSize <- value
    settings$tileSize <- value + 2*settings$overhangSize
    newIndex <- .makeTiles(settings)
    slot(object, "index") <- newIndex
    return(object)
}

## change the tile size and return the modified object
.changeTileSize <- function(object, value) {
    settings <- tileSettings(object)
    settings$tileSize <- value
    settings$chunkSize <- value - 2*settings$overhangSize
    newIndex <- .makeTiles(settings)
    slot(object, "index") <- newIndex
    return(object)
}

## change the overhang size and return the modified object
.changeOverhangSize <- function(object, value) {
    settings <- tileSettings(object)
    settings$overhangSize <- value
    settings$tileSize <- settings$chunkSize + 2*value
    newIndex <- .makeTiles(settings)
    slot(object, "index") <- newIndex
    return(object)
}

## change the chromosome list and return the modified object
.changeChromosomes <- function(object, value) {
    chroms <- GenomeInfoDb::seqlevels(object)
    matching <- value %in% chroms
    if(!all(matching)) {
        warning("Trying to add new chromosomes without data.
Only the ones matching will be used")
    }
    idx <- getIndex(object)
    GenomeInfoDb::seqlevels(idx, force = TRUE) <- value
    GenomeInfoDb::seqlevels(metadata(idx)$chromosomes, force = TRUE) <- value
    metadata(idx)$numTiles <- length(idx)
    slot(object, "index") <- idx
    return(object)
}

## change the number of tiles and return the modified object
.changeNumberOfTiles <- function(object, value) {
    settings <- tileSettings(object)
    size <- settings$chunkSize*settings$numTiles
    newChunkSize <- round(size/value)
    return(.changeChunkSize(object, newChunkSize))
}

#' @rdname changeSettings
#' @export
setGeneric("changeSettings", function(object, param, value) standardGeneric("changeSettings")) 

#' Check data compliance with tile settings
#'
#' Check if the indices were build correctly, according to the
#' specified parameters. This is the recommended way of changing
#' tile settings, as it triggers instant recomputation of the index.
#'
#' @rdname changeSettings
#' @param object A /code{GenomicTiles} object.
#' @param param The name of a tile settings parameter.
#' @param value An appropriate value. In most cases integer.
#' @return A /code{GenomicTiles} object
#' @examples 
#' gt <- makeTestGenomicTiles()
#' gt2 <- changeSettings(gt, "chunkSize", 20)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @export
setMethod("changeSettings", signature(object = "GenomicTiles",
                                      param = "character",
                                      value = "ANY"), function(object, param, value) {
    .changeSettings(object, param, value)
})

## Coercion
## ========

# converts the GenomicTiles object to a DataFrame
.GenomicTilesToDataFrame <- function(from) {
    df <- as.data.frame(rowRanges(from))
    if(all(dim(from) == c(0, 0))) {
        return(DataFrame(df))
    }
    
    numDataFrames <- length(assays(from))
    gtDim <- dim(from)
    
    res <- unlist(assays(from))
    expNames <- names(assays(from))
    if(is.null(expNames)) expNames <- paste0("V", 1:length(assays(from)))
    expNamesRepeats <- rep(gtDim[1], numDataFrames)
    res$assay <- Rle(factor(rep(expNames, expNamesRepeats)))

    coldt <- colData(from)
    sampleVars <- unique(names(coldt))
    meta <- lapply(sampleVars, function(y) {
        idx <- which(names(coldt) == y)
        return(coldt[,idx])
    })
    names(meta) <- sampleVars
    
    gtdf <- cbind(df, res)
    metadata(gtdf) <- meta  
    
    return(gtdf)
}

#' GenomicTiles to DataFrame
#' 
#' @name asDataFrame
#' @family GenomicTiles
setAs(from = "GenomicTiles", to = "DataFrame", def = function(from) {
    .GenomicTilesToDataFrame(from)
})

setMethod("as.data.frame", "GenomicTiles", function(x) {
    ## Bug in as(x, "data.frame"), replace with
    as.data.frame(as(x, "DataFrame"))
})

## Subsetting
## ==========

#' Method to subset the indeces of GenomicTiles.
#'
#' Subsetting the indeces of GenomicTiles based on a SummarizedExperiment.
#'
#' @param x A SummarizedExperiment object.
#' @param index The index of a GenomicTiles to be subsetted.
#' @return A list of two GRanges objects: the subsetted index and the
#' subsetted coordinates.
.subsetIndeces <- function(se, index) {
    gpCoords <- rowRanges(se)@pos_runs
    gpDims <- length(gpCoords)
    ir <- IRanges(start = cumsum(c(1, end(gpCoords)[-gpDims[1]])),
                  end = cumsum(width(gpCoords)))
    coords <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(gpCoords), ir)
    l <- metadata(index)
    l$chromosomes <- gpCoords
    
    indx <- .makeTiles(l)
    GenomeInfoDb::seqinfo(indx) <- GenomeInfoDb::seqinfo(gpCoords)
    GenomeInfoDb::seqlevels(coords) <- GenomeInfoDb::seqlevels(coords)[GenomeInfoDb::seqlevels(coords) %in%
                                           unique(GenomeInfoDb::seqnames(coords))]
    GenomeInfoDb::seqlevels(indx) <- GenomeInfoDb::seqlevels(indx)[GenomeInfoDb::seqlevels(indx) %in%
                                           unique(GenomeInfoDb::seqnames(indx))]
    return(list(index = indx, coordinates = coords))
}

#' Method to subset GenomicTiles.
#'
#' Subsetting the GenomicTiles based on the subset method of
#' SummarizedExperiment. Additionally the index and coordinates gets
#' subsetted.
#'
#' @param x A GenomicTiles object.
#' @param ... Any subset option accepted by SummarizedExperiment. Usually
#' in accordance to base::subset.
#' @return A subsetted GenomicTiles object.
.subsetGenomicTiles <- function(x,...){
    if(all(dim(x) == c(0, 0))) return(x)
    se <- subset(SummarizedExperiment(assays = assays(x),
                                      rowRanges = rowRanges(x),
                                      colData = colData(x)), ...)
    GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))

    indeces <- .subsetIndeces(se, getIndex(x))
    GenomeInfoDb::seqlevels(metadata(indeces$index)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
    return(new("GenomicTiles", index = indeces$index,
               coordinates = indeces$coordinates, se))
}


#' Subset method for /code{GenomciTiles}
#'
#' Subsetting the /code{GenomicTiles} by a logical statement
#'
#' @param x A /code{GenomicTiles} object.
#' @param ... Further arguments. Mostly a logical statement.
#' Note that the columnnames for chromosomes and positions
#' are: seqnames and pos.
#' @return A subsetted /code{GenomicTiles} object.
#' @examples
#' gt <- makeTestGenomicTiles()
#' res <- subset(gt, seqnames == "chrI" & pos <= 50)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subset", "GenomicTiles", function(x, ...) {
    .subsetGenomicTiles(x, ...)
})

#' Method to subset GenomicTiles by a GRanges object.
#'
#' Subsetting the GenomicTiles by a GRanges object based on the
#' subsetByOverlaps method applying to SummarizedExperiment.
#' Additionally the index and coordinates gets subsetted.
#'
#' @param query A GenomicTiles object.
#' @param subject A GRanges object.
#' @param ... Any other parameters applicable to subset a
#' SummarizedExperiment by overlaps.
#' @return A subsetted GenomicTiles object.
.subsetGenomicTilesByOverlaps <- function(query, subject, ...){
    se <- subsetByOverlaps(SummarizedExperiment(assays = assays(query),
                                      rowRanges = rowRanges(query)),
                           subject, ...)
    indeces <- .subsetIndeces(se, getIndex(query))
    return(new("GenomicTiles", index = indeces$index,
               coordinates = indeces$coordinates, se))
}

.subsetByOverlaps <- function(query, subject, maxgap, minoverlap,
                              type, ...) {
    rowRanges <- rowRanges(query)
    assay <- assays(query)
    index <- getIndex(query)
    colData <- colData(query)

    tiles <- subsetByOverlaps(index, subject)
    se <- subsetByOverlaps(SummarizedExperiment(assay, rowRanges = rowRanges, colData = colData),
                           tiles, ...)
    GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
    GenomeInfoDb::seqlevels(tiles, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
    if(length(metadata(tiles)) > 0) {
        metadata(tiles)$numTiles <- length(tiles)
        GenomeInfoDb::seqlevels(metadata(tiles)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
    }
    
    coords <- .makeCoordinates(rowRanges(se))
    
    res <- new("GenomicTiles", index = tiles,
               coordinates = coords, se)
    return(res)
}

.exactSubsetByOverlaps <- function(query, subject, ...) {
    rowRanges <- rowRanges(query)
    assay <- assays(query)
    settings <- tileSettings(query)
    colData <- colData(query)

    se <- subsetByOverlaps(SummarizedExperiment(assay, rowRanges = rowRanges, colData = colData),
                           subject, ...)
    tiles <- .makeTiles(settings)
    GenomeInfoDb::seqlevels(rowRanges(se), force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
    GenomeInfoDb::seqlevels(tiles, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(rowRanges(se))
    if(length(metadata(tiles)) > 0) {
        GenomeInfoDb::seqlevels(metadata(tiles)$chromosomes, force = TRUE) <- GenomeInfoDb::seqlevelsInUse(tiles)
    }
    
    coords <- .makeCoordinates(rowRanges(se))
    
    res <- new("GenomicTiles", index = tiles,
               coordinates = coords, se)
    return(res)
}

#' Subset by overlaps method for \code{GenomciTiles}
#'
#' Subsetting the \code{GenomicTiles} by a \code{GRanges} object
#'
#' @param query A \code{GenomicTiles} object.
#' @param subject A \code{GRanges} object
#' @param maxgap,minoverlap Intervals with a separation of \code{maxgap} or
#' less and a minimum of \code{minoverlap} overlapping positions, allowing for
#' \code{maxgap}, are considered to be overlapping.  \code{maxgap} should
#' be a scalar, non-negative, integer. \code{minoverlap} should be a scalar,
#' positive integer.
#' @param type By default, any overlap is accepted. By specifying the \code{type}
#' parameter, one can select for specific types of overlap. The types correspond
#' to operations in Allen's Interval Algebra (see references). If \code{type}
#' is \code{start} or \code{end}, the intervals are required to have matching
#' starts or ends, respectively. While this operation seems trivial, the naive
#' implementation using \code{outer} would be much less efficient. Specifying
#' \code{equal} as the type returns the intersection of the \code{start} and
#' \code{end} matches. If \code{type} is \code{within}, the query interval must
#' be wholly contained within the subject interval. Note that all matches must
#' additionally satisfy the \code{minoverlap} constraint described above.
#'
#' The \code{maxgap} parameter has special meaning with the special
#' overlap types. For \code{start}, \code{end}, and \code{equal}, it specifies
#' the maximum difference in the starts, ends or both, respectively. For
#' \code{within}, it is the maximum amount by which the query may be wider
#' than the subject.
#' @param ... Additional parameters
#' @return A subsetted \code{GenomicTiles} object.
#' @examples
#' gt <- makeTestGenomicTiles()
#' gr <- GRanges(c("chrI", "chrII"), IRanges(c(1, 120), c(40, 150)))
#' res <- subsetByOverlaps(gt, gr)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setMethod("subsetByOverlaps", c("GenomicTiles", "GRanges"),
          function(query, subject, maxgap=0L, minoverlap=1L,
                      type=c("any", "start", "end", "within", "equal"), ...) {
              .subsetByOverlaps(query, subject, maxgap = maxgap,
                                minoverlap = minoverlap,
                                type = type, ...)
          })

#' Method to extract all GenomicTiles at once.
#'
#' Extracting and converting all GenomicTiles at once using the complete
#' assays. This is more memory consuming, but faster for smaller datasets.
#'
#' @param gt A GenomicTiles object.
#' @param index A GRanges object.
#' @param chromosomes A character vector of chromosome names.
#' @return A SimpleList of DataFrames with as many elements
#' as there are tiles.
.extractFullGenomicTiles <- function(gt, index, chromosomes) {
    rowindx <- getIndexCoordinates(gt, index = index)
    coords <- unlistCoordinates(gt, index = rowindx,
                                chromosomes = chromosomes)
    fullIndex <- sort(do.call(c, unname(coords)))
    df <- as(gt, "DataFrame")
    meta <- metadata(df)
    vec <- rep(mcols(fullIndex)$id, width(fullIndex))

    applyVec <- IRanges(start=(0:(length(fullIndex)/length(rowindx) - 1))*length(rowindx) + 1,
                        end = (1:(length(fullIndex)/length(rowindx)))*length(rowindx))

    ## From in to out: For each assay make DataFrameList,
    ## combine and unlist the List to get a DataFrame of tiles that are consecutively placed.
    ## Because of overlaps this DataFrame is bigger than the original and the same size
    ## as "vec".
    dfList <- splitAsList(unlist(do.call(c, lapply(1:length(applyVec), function(y) {
        extractList(df, fullIndex[start(applyVec)[y] : end(applyVec)[y]])
    }))), vec)

    metadata(dfList) <- meta
    return(dfList)
}

#' Method to extract all GenomicTiles by blocks.
#'
#' Extracting and converting all GenomicTiles iteratively by blocks of
#' GenomicTiles. This is less memory consuming, but slower for smaller
#' datasets.
#'
#' @param gt A GenomicTiles object.
#' @param index A GRanges object.
#' @param blocks An integer value representing the number of blocks.
#' @return A SimpleList of DataFrames with as many elements
#' as there are tiles.
.extractGenomicTilesByBlocks <- function(gt, index, blocks) {
    blockSize <- length(index)/blocks
    blockList <- unname(split(mcols(index)$id,
                      ceiling(seq_along(mcols(index)$id)/blockSize)))
    
    dfList <- do.call(c, bplapply(blockList, function(y) {
        blockindx <- index[mcols(index)$id == y]
        subgt <- .subsetByOverlaps(gt, blockindx)
        df <- DataFrame(subgt)
        coords <- unlistCoordinates(subgt,
                                    getIndexCoordinates(subgt))
        fullIndex <- do.call(c, unname(coords))
        vec <- rep(mcols(fullIndex)$id, width(fullIndex))
        dfList <- splitAsList(df, vec)
        return(dfList)
    }))
    return(dfList)
}

.extractGenomicTilesByIndex <- function(gt, index, size = 3e7) {
    chromosomes <- GenomeInfoDb::seqlevels(index)
    gtDims <- dim(gt)
    if(all(gtDims == c(0, 0))) return(DataFrameList())
    space <- getChunkSize(gt)*length(index)
    
    if(space <= size) {
        dfList <- .extractFullGenomicTiles(gt, index, chromosomes)
    }
    else {
        dfList <- .extractGenomicTilesByBlocks(gt, index,
                                               ceiling(space/size))
    }
    return(dfList)
}

#' Get a tile.
#'
#' Extracting one tile from a GenomicTiles object.
#'
#' @param object A GenomicTiles Object.
#' @param id A tile id as an integer
#' @return A data.frame of the tile.
.getTile <- function(gt, id) {
    index <- getIndex(gt)
    range <- ranges(index[index$id == id,])
    chrom <- as.character(GenomeInfoDb::seqnames(index[index$id == id,]))
    gtSubset <- subset(gt, GenomeInfoDb::seqnames(gt) == chrom &
                       pos >= start(range) & pos <= end(range))

    df <- DataFrame(gtSubset)
    names(df)[1] <- "chromosome"
    ## df <- melt(df, measure.vars = rownames(colData(gtSubset)), variable.name = "sample")
    return(DataFrameList(df))
}

#' @rdname getTile
#' @export
setGeneric("getTile", function(object, id, ...) standardGeneric("getTile"))

#' Tile extraction as a DataFrame
#'
#' Extracting one or multiple tiles from a \code{GenomicTiles} object
#' and coercing them to a DataFrameList.
#'
#' @param object A \code{GenomicTiles} object
#' @param id A vector of tile ids
#' @param size The maximal number of rows that should be handled at once.
#' If the dataset is bigger it will be processed in chunks. This is to lower
#' memory consumption on big datasets, which in turn is slower.
#' @param ... Additional arguments
#' @return A \code{SimpleDataFrameList}
#' @examples
#' gt <- makeTestGenomicTiles()
#' getTile(gt, 1:3)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname getTile
#' @export
setMethod("getTile", "GenomicTiles", function(object, id, size = 3e7) {
    if(missing(id)) {
        indx <- getIndex(object)
        id <- mcols(indx)$id
        res <- .extractGenomicTilesByIndex(object,
                                           indx,
                                           size = size)
    }
    else {
        if(length(id) == 1) res <- .getTile(object, id)
        if(length(id) > 1) {
            indx <- getIndex(object, id)
            res <- .extractGenomicTilesByIndex(object, indx, size = size)
        }
    }
    meta <- tileSettings(object)
    meta$chunks <- getChunkIndex(object, id)
    metadata(res) <- meta
    return(res)
})

.subsetByDoubleBrackets <- function(x, i, j) {
    if(!missing(j)) {
        stop("Wrong number of dimensions")
    }

    if(length(i) > 1L) {
        stop("attempt to extract more than one element")
    }

    return(.getTile(x, i))
}

#' Providing pseudo-list functionality
#'
#' Getting a specific tile
#'
#' @param x A \code{GenomicTiles} object
#' @param i An integer (for '[[') or a \code{GRanges} object (for '[')
#' @return A \code{DataFrame} (for '[[') or a subsetted \code{GenomicTiles} object (for '[')
#' @rdname GenomicTiles-brackets
setMethod("[[", c("GenomicTiles", "numeric"),
          function(x, i) .subsetByDoubleBrackets(x, i))

#' @rdname GenomicTiles-brackets
setMethod("[", c("GenomicTiles", "GRanges"),
          function(x, i) .exactSubsetByOverlaps(x, i))

#' @rdname GenomicTiles-view
#' @export
setGeneric("view", function(object, ...) standardGeneric("view"))

#' View the dataset
#'
#' Cbinding the columns all together and coercing to data.frame
#'
#' @param object A \code{GenomicTiles} object
#' @param ranges A \code{GRanges} object. Makes it possible to
#' select regions by \code{GRanges}. Either ranges or seqnames, start and
#' end must be supplied
#' @param seqnames A chromosomes name. Either ranges or seqnames, start and
#' end must be supplied
#' @param start A start site. Either ranges or seqnames, start and
#' end must be supplied
#' @param end An end site. Either ranges or seqnames, start and
#' end must be supplied
#' @param ... Additional arguments
#' @return A data.frame of the selected data.
#' @examples
#' gt <- makeTestGenomicTiles()
#' gr <- GRanges(c("chrI", "chrII"), IRanges(c(1, 10), c(40, 30)))
#' head(view(gt, ranges = gr))
#' head(view(gt, seqnames = "chrI", start = 1, end = 20))
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenomicTiles-view
#' @export
setMethod("view", c("GenomicTiles"), function(object, ranges = NULL, seqnames = NULL,
                                              start = NULL, end = NULL) {
    if((is.null(seqnames) | is.null(start) | is.null(end)) & is.null(ranges)) {
        return(as.data.frame(object))
    }
    if(!is.null(ranges)) {
        res <- .exactSubsetByOverlaps(object, ranges)
    }
    else {
        res <- subset(object, seqnames == seqnames & pos >= start & pos <= end)
    }
    return(as.data.frame(res))
})

## Tile computation
## ================

#' compute metrics for each tile
.MetricsFun <- function(x, what, na.rm = FALSE) {

    indx <- getIndexCoordinates(x)
    
    res <- lapply(1:length(assays(x)), function(ii) {
        sapply(rownames(colData(x)), function(y) {
            rle <- extractList(assay(x, ii)[[y]], ranges(indx))
            eval(call(what, rle, na.rm = na.rm))
        })
    })

    names(res) <- names(assays(x))
    return(res)
}

#' Computing metrics
#'
#' Computing metrics on each tile of the \code{GenomicTiles} object.
#' So far all metrics from the Summary generics group, as well as
#' mean, var, sd, median, mad and IQR are supported.
#'
#' @param x A \code{GenomicTiles} object
#' @param ... Additional arguments
#' @param na.rm Should NAs be dropped. Otherwise the result is NA
#' @return A list of as many elements as there are assays.
#' Each element contains of a matrix with the specified
#' metric computed per tile per column of the assay data.
#' @examples
#' gt <- makeTestGenomicTiles()
#' sum(gt)
#' min(gt)
#' max(gt)
#' mean(gt)
#' var(gt)
#' sd(gt)
#' median(gt)
#' mad(gt)
#' IQR(gt)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenomicTiles-metrics
setMethod("Summary", "GenomicTiles", function(x, ..., na.rm = FALSE) {
    indx <- getIndexCoordinates(x)
    
    res <- lapply(1:length(assays(x)), function(ii) {
        sapply(rownames(colData(x)), function(y) {
            rle <- extractList(assay(x, ii)[[y]], ranges(indx))
            (getFunction(.Generic))(rle, na.rm = na.rm)
        })
    })

    names(res) <- names(assays(x))
    return(res)
})

#' @rdname GenomicTiles-metrics
setMethod("mean", "GenomicTiles", function(x) {
    .MetricsFun(x, "mean")
})

#' @rdname GenomicTiles-metrics
setMethod("var", "GenomicTiles", function(x) {
    .MetricsFun(x, "var")
})

#' @rdname GenomicTiles-metrics
setMethod("sd", "GenomicTiles", function(x) {
    .MetricsFun(x, "sd")
})

#' @rdname GenomicTiles-metrics
setMethod("median", "GenomicTiles", function(x) {
    .MetricsFun(x, "median")
})

#' @rdname GenomicTiles-metrics
setMethod("mad", "GenomicTiles", function(x) {
    .MetricsFun(x, "mad")
})

#' @rdname GenomicTiles-metrics
setMethod("IQR", "GenomicTiles", function(x) {
    .MetricsFun(x, "IQR")
})

## Cosmetics
## =========

.showGenomicTiles <- function(object) {
    cl <- class(object)
    dims <- dim(object)
    as <- assays(object)
    rd <- names(rowData(object))
    md <- unique(names(colData(object)))
    cnames <- colnames(object)
    cdata <- names(colData(object))

    if(length(tileSettings(object)) != 0) {
        tsize <- tileSettings(object)$tileSize
        tname <- "tiles"
        chunk <- tileSettings(object)$chunk
        if (chunk) {
            tsize <- tileSettings(object)$chunkSize
            tname <- "chunks"
        }
        unit <- ""
        if(!is.null(tsize)) {
            unit = "bp"
            if(tsize/1000 > 1) {
                unit <- "kbp"
                tsize <- tsize/1000
            }
            if(tsize/1e6 > 1) {
                unit <- "Mbp"
                tsize <- tsize/1e6
            }
        }
        chroms <- GenomeInfoDb::seqlevels(tileSettings(object)$chromosomes)
    }
    tnum <- length(getIndex(object))
    
    
    cat("class:", cl, "\n")
    cat("dimension:", dims, "\n")
    cat(paste0("assays(", length(as), "):"), names(as), "\n")
    cat(paste0("position variables(", length(rd), "):"), rd, "\n")
    cat(paste0("sample variables(", length(md), "):"), md, "\n")
    cat(paste0("samples(", length(cnames), "):"), cnames, "\n")
    if(length(tileSettings(object)) != 0) {
        cat(paste0(tname, " size: ", tsize, unit), "\n")
        cat(paste0("number of ", tname, ": ", tnum), "\n")
        cat("chromosomes:", chroms, "\n")
    }
}

## Show method for GenomicTiles.
setMethod("show", "GenomicTiles", function(object) {
    .showGenomicTiles(object)
})


