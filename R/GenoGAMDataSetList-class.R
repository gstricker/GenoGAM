## ====================
## GenoGAMDatasetList class
## ====================
#' @include GenoGAMSettings-class.R
#' @include GenoGAMDataSet-class.R
#' @include GenoGAM-class.R
NULL

setClassUnion("HDF5OrMatrix", c("matrix", "HDF5Matrix"))

#' GenoGAMDataSetList
#'
#' The GenoGAMDataSetList class contains the pre-processed raw data and
#' additional slots that define the input and framework for the model.
#' It extends upon the idea of the GenoGAMDataSet class to make it possible
#' to store genomes and data of size > 2^31 (maximum size of integers in R).
#' Thus the only difference to a GenoGAMDataSet is the arrangement as a
#' list of RangedSummarizedExperiments under the hood. On the surface is
#' does still behave like a GenoGAMDataSet. It is not intended to be used
#' by the user. For more information check the GenoGAMDataSet class documentation.
#' 
#' @slot settings The global and local settings that will be used to compute the
#' model.
#' @slot design The formula describing how to evaluate the data. See details.
#' @slot sizeFactors The normalized values for each sample. A named numeric vector.
#' @slot index A GRanges object representing an index of the ranges defined 
#' on the genome. Mostly used to store tiles.
#' @slot data A list of RangedSummarizedExperiment objects
#' @slot id A GRanges object keeping the identifiers assigning the regions to the
#' respective list elements
#' @slot hdf5 A logical slot indicating if the object should be stored as HDF5
#' @slot countMatrix Either a matrix or HDF5Matrix to store the sums of counts of
#' the regions (could also be seen as bins) for later use especially by DESeq2
#' @name GenoGAMDataSetList-class
#' @rdname GenoGAMDataSetList-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("GenoGAMDataSetList",
         slots = list(settings = "GenoGAMSettings",
                      design = "formula", sizeFactors = "numeric",
                      index = "GRanges", data = "list", id = "GRanges",
                      hdf5 = "logical", countMatrix = "HDF5OrMatrix"),
         prototype = list(settings = GenoGAMSettings(),
                          design = ~ s(x), sizeFactors = numeric(), 
                          index = GenomicRanges::GRanges(),
                          data = list(), id = GenomicRanges::GRanges(),
                          hdf5 = FALSE, countMatrix = matrix()))

## Validity
## ========

.validateDataType <- function(object) {
    if(!is(object@data, "list")) {
        return("'data' must be a list object and all elements must be of class RangedSummarizedExperiment")
    }
    NULL
}

.validateIDType <- function(object) {
    if(!is(object@id, "GRanges")) {
        return("'id' must be a GRanges object")
    }
    NULL
}

.validateGGDLChromosomes <- function(object) {
    cindex <- GenomeInfoDb::seqlevels(object@index)
    seqlev <- lapply(object@data, function(y) {
        GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(y))
    })
    cobject <- unique(unlist(seqlev))
    if(!all(cindex %in% cobject)) {
        return("Different chromosomes for data and index objects.")
    }
    NULL
}

## general validate function
.validateGenoGAMDataSetList <- function(object) {
    c(.validateSettingsType(object), .validateDesignType(object),
      .validateSFType(object), .validateIndexType(object),
      .validateGGDLChromosomes(object),
      .validateDataType(object),
      .validateIDType(object),
      .validateH5Type(object),
      .validateCountMatrixType(object))
}

S4Vectors::setValidity2("GenoGAMDataSetList", .validateGenoGAMDataSetList)

## Constructor
## ===========

#' GenoGAMDataSetList constructor.
#'
#' GenoGAMDataSetList is the constructor function for the GenoGAMDataSetList-class.
#'
#' @param ... The slots and their respective values
#' @param object,x For use of S4 methods. The GenoGAMDataSetList object.
#' @param value,i For use of S4 methods. The value to be assigned to the slot.
#' @param withDimnames For use of S4 methods.
#' 
#' @return An object of class GenoGAMDataSetList
#' @name GenoGAMDataSetList
#' @rdname GenoGAMDataSetList-class
GenoGAMDataSetList <- function(...) {
    return(new("GenoGAMDataSetList", ...))
}

#' Make an example /code{GenoGAMDataSet}
#'
#' @return A /code{GenoGAMDataSet} object
#' @examples
#' ggdl <- makeTestGenoGAMDataSetList()
#' @export
makeTestGenoGAMDataSetList <- function() {

    k <- 10000
    sinCurve <- sin(seq(-7, 5, length.out = k)) + 1
    ip <- rnbinom(k, size = 2, mu = sinCurve/max(sinCurve))
    sinCurve <- c(sin(seq(-7, -1, length.out = k/2)) + 1, runif(k/2, 0, 0.2))
    background <- rnbinom(k, size = 2, mu = sinCurve/max(sinCurve)/2)
    chroms <- c("chrX", "chrY", "chrZ")
    gr <- GenomicRanges::GPos(GenomicRanges::GRanges(chroms,
                                                     IRanges::IRanges(c(1,1,1), c(k, k, k))))
    GenomeInfoDb::seqlengths(gr) <- c(1e6, 2e6, 1.5e6)

    ## colData
    coldf <- S4Vectors::DataFrame(experiment = c(0, 1))
    rownames(coldf) <- c("input", "IP")
    
    selist <- lapply(chroms, function(chr) {
        add <- sample(0:3, 2)
        background <- background + add[1]
        ip <- ip + add[2]
        df <- S4Vectors::DataFrame(input = background, IP = ip)
        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = gr[GenomeInfoDb::seqnames(gr) == chr,],
                                                         assays = list(df), colData = coldf)
    })
    names(selist) <- chroms
    id <- .extractGR(gr)
    id$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(gr)))

    ## make tiles
    l <- list(chromosomes = .extractGR(gr),
              chunkSize = 2000,
              overhangSize = 250)
    tiles <- .makeTiles(l)

    ## size factors
    sf <- rep(0, length(colnames(selist[[1]])))
    
    ggdl <- GenoGAMDataSetList(data = selist, id = id, design = ~ s(x),
                              index = tiles, sizeFactors = sf)
    
    design(ggdl) <- ~ s(x) + s(x, by = experiment)
        
    return(ggdl)
}

## Check function
## ===============

#' @noRd
setMethod("checkObject", "GenoGAMDataSetList", function(object) {
    .checkGenoGAMDataSet(object)
})

## Accessors
## =========

#' @describeIn GenoGAMDataSetList Get the dimension of the object
setMethod("dim", "GenoGAMDataSetList", function(x) {
    dims <- sapply(x@data, dim)
    if(length(dims) == 0) {
        return(c(0, 0))
    }
    res <- c(sum(as.numeric(dims[1,])), max(as.numeric(dims[2,])))
    return(res)
})

#' @describeIn GenoGAMDataSetList The length of the object
setMethod("length", "GenoGAMDataSetList", function(x) {
    dim(x)[1]
})

#' @describeIn GenoGAMDataSetList The seqlengths of the object
setMethod("seqlengths", "GenoGAMDataSetList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(GenomeInfoDb::seqlengths(rowRanges(x)[[1]]))
    }

    integer()
})

#' @describeIn GenoGAMDataSetList The seqlevels of the object
setMethod("seqlevels", "GenoGAMDataSetList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(GenomeInfoDb::seqlevels(rowRanges(x)[[1]]))
    }

    character()
})

#' @describeIn GenoGAMDataSetList The seqlevelsInUse of the object
setMethod("seqlevelsInUse", "GenoGAMDataSetList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(names(rowRanges(x)))
    }

    character()
})

#' @describeIn GenoGAMDataSetList get colData from the first element of the
#' SummarizedExperiment list
setMethod("colData", "GenoGAMDataSetList", function(x, ...) {
    if(length(x@data) == 0) {
        return(DataFrame())
    }
    
    colData(x@data[[1]])
})

#' @describeIn GenoGAMDataSetList get a list of rowRanges from the
#' GenoGAMDataSetList object
setMethod("rowRanges", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@data, rowRanges)
})

#' @describeIn GenoGAMDataSetList get a list of assays from the
#' GenoGAMDataSetList object
setMethod("assay", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@data, assay)
})

#' @describeIn GenoGAMDataSetList get a list of list of assays from the
#' GenoGAMDataSetList object. Just for completeness, shouldn't be needed.
setMethod("assays", "GenoGAMDataSetList", function(x, ...) {
    lapply(x@data, assays)
})

#' @describeIn GenoGAMDataSetList get colnames from the first element of the
#' SummarizedExperiment list
setMethod("colnames", "GenoGAMDataSetList", function(x) {
    if(length(x@data) == 0) {
        return(NULL)
    }
        
    colnames(x@data[[1]])
})

#' @describeIn GenoGAMDataSetList accessor to the index slot
setMethod("getIndex", "GenoGAMDataSetList", function(object) {
    return(slot(object, "index"))
})

#' @describeIn GenoGAMDataSetList An accessor to the countMatrix slot
setMethod("getCountMatrix", signature(object = "GenoGAMDataSetList"), function(object) {
    res <- slot(object, "countMatrix")
    colnames(res) <- colnames(object)
    return(res)
})

#' @describeIn GenoGAMDataSetList The accessor to the list of settings, that
#' were used to generate the tiles.
setMethod("tileSettings", "GenoGAMDataSetList", function(object) {
    S4Vectors::metadata(getIndex(object))
})

.rowRangesFromList <- function(object) {
    res <- lapply(object, function(x) {
        .extractGR(rowRanges(x))
    })
    return(do.call("c", unname(res)))
}

#' @describeIn GenoGAMDataSetList The actual underlying GRanges showing the range of the data.
setMethod("dataRange", "GenoGAMDataSetList", function(object) {
    res <- .rowRangesFromList(object@data)
})

#' @describeIn GenoGAMDataSetList A GRanges object representing the chromosomes
#' or chromosome regions on which the model will be computed
setMethod("getChromosomes", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chromosomes
})

#' @describeIn GenoGAMDataSetList The size of the tiles
setMethod("getTileSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$tileSize
})

#' @describeIn GenoGAMDataSetList The size of the chunks
setMethod("getChunkSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$chunkSize
})

#' @describeIn GenoGAMDataSetList The size of the overhang (on one side)
setMethod("getOverhangSize", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$overhangSize
})

#' @describeIn GenoGAMDataSetList The total number of tiles
setMethod("getTileNumber", "GenoGAMDataSetList", function(object) {
    tileSettings(object)$numTiles
})

#' @describeIn GenoGAMDataSetList A boolean function that is true if object uses HDF5 backend
setMethod("is.HDF5", signature(object = "GenoGAMDataSetList"), function(object) {
    res <- slot(object, "hdf5")
    return(res)
})

#' @describeIn GenoGAMDataSetList Access to the design slot.
setMethod("design", "GenoGAMDataSetList", function(object) {
    slot(object, "design")
})

#' @describeIn GenoGAMDataSetList Replace method of the design slot.
setReplaceMethod("design", "GenoGAMDataSetList", function(object, value) {
    newCols <- as.vector(na.omit(.getVars(value)))
    if(!all(newCols %in% colnames(colData(object)))) {
        futile.logger::flog.error("'by' variables could not be found in colData")
        stop("'by' variables could not be found in colData")
    }
    slot(object, "design") <- value
    return(object)
})

#' @describeIn GenoGAMDataSetList Access to the sizeFactors slot
setMethod("sizeFactors", "GenoGAMDataSetList", function(object) {
    sf <- slot(object, "sizeFactors")
    names(sf) <- colnames(object)
    return(sf)
})

#' @describeIn GenoGAMDataSetList Replace method of the sizeFactors slot
setReplaceMethod("sizeFactors", "GenoGAMDataSetList", function(object, value) {
    slot(object, "sizeFactors") <- value
    return(object)
})


## change Settings
## ===============

#' @describeIn GenoGAMDataSetList Replace method of the chunkSize parameter,
#' that triggers a new computation of the tiles based on the new chunk size.
setReplaceMethod("getChunkSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$chunkSize <- value
                     settings$tileSize <- value + 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

#' @describeIn GenoGAMDataSetList Replace method of the tileSize parameter,
#' that triggers a new computation of the tiles based on the new tile size.
setReplaceMethod("getTileSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$tileSize <- value
                     settings$chunkSize <- value - 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

#' @describeIn GenoGAMDataSetList Replace method of the overhangSize parameter,
#' that triggers a new computation of the tiles based on the new overhang size.
setReplaceMethod("getOverhangSize", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     settings$overhangSize <- value
                     settings$tileSize <- settings$chunkSize + 2*value
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

#' @describeIn GenoGAMDataSetList Replace method of the tileNumber parameter,
#' that triggers a new computation of the tiles based on the new number of tiles.
setReplaceMethod("getTileNumber", signature = c("GenoGAMDataSetList", "numeric"),
                 function(object, value) {
                     settings <- tileSettings(object)
                     size <- min(width(settings$chromosomes))
                     if(size > sum(width(dataRange(object)))) {
                         warning("The settings indicated a longer total genome size than actually present. it was trimmed accordingly.")
                         size <- sum(width(dataRange(object)))
                     }
                     settings$chunkSize <- round(size/value)
                     settings$tileSize <- value + 2*settings$overhangSize
                     newIndex <- .makeTiles(settings)
                     slot(object, "index") <- newIndex
                     return(object)
                 })

## Subsetting
## ==========
#' Subset method for GenoGAMDataSetList
#'
#' @details
#' Those are various methods to subset the GenoGAMDataSetList object.
#' By logical statement or GRanges overlap. The '[' subsetter is
#' just a short version of 'subsetByOverlaps'.
#'
#' @param x A GenoGAMDataSetList object.
#' @param ranges,i A GRanges object. In case of subsetting by double brackets
#' 'i' is the index of the tile.
#' @param maxgap,minoverlap Intervals with a separation of 'maxgap' or
#' less and a minimum of 'minoverlap' overlapping positions, allowing for
#' 'maxgap', are considered to be overlapping. 'maxgap' should
#' be a scalar, non-negative, integer. 'minoverlap' should be a scalar,
#' positive integer.
#' @param type By default, any overlap is accepted. By specifying the 'type'
#' parameter, one can select for specific types of overlap. The types correspond
#' to operations in Allen's Interval Algebra (see references in). If \code{type}
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
#' @param invert If TRUE, keep only the query ranges that do _not_ overlap
#' the subject.
#' @param ... Further arguments. Mostly a logical statement
#' in case of the 'subset' function. Note that the columnnames
#' for chromosomes and positions are: 'seqnames' and 'pos'.
#' @references
#' Allen's Interval Algebra: James F. Allen: Maintaining knowledge
#' about temporal intervals. In: Communications of the ACM.
#' 26/11/1983. ACM Press. S. 832-843, ISSN 0001-0782
#' @return A subsetted GenoGAMDataSetList object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenoGAMDataSetList-subsetting
setMethod("subset", "GenoGAMDataSetList", function(x, ...) {
    if(all(dim(x) == c(0, 0))) return(x)
    
    settings <- slot(x, "settings")
    design <- design(x)
    sf <- sizeFactors(x)
    h5 <- slot(x, "hdf5")

    ## make data slot
    ## subset all SummarizedExperiments
    se <- lapply(x@data, function(se) {
        res <- subset(se, ...)
        return(res)
    })

    ## reduce the list to non-empty. However in order to mimic what happens
    ## if the complete subset is zero, we have to keep the result, as we might
    ## need it later. Otherwise it will be an empty list, which is not how
    ## SummarizedExperiment and GenoGAMDataSet deals with it. On the other hand
    ## we don't want to carry a bunch of empty SEs around if they are not needed anyway.
    reduced_se <- se[as.logical(vapply(se, length, integer(1)))]

    if(length(reduced_se) == 0) {
        reduced_se <- se[1]
        index <- GenomicRanges::GRanges()
    }
    else {
        reduced_se <- lapply(reduced_se, function(y) {
            ## set correct seqinfo
            GenomeInfoDb::seqlevels(rowRanges(y), pruning.mode="coarse") <- names(reduced_se)
            return(y)
        })
        slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(reduced_se[[1]])
        index <- .subsetIndex(reduced_se, getIndex(x))
    }

    ## make id slot
    splitid <- .rowRangesFromList(reduced_se)
    splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))

    ggdl <- new("GenoGAMDataSetList", settings = settings,
               design = design, sizeFactors = sf, index = index,
               data = reduced_se, id = splitid, hdf5 = h5)
    return(ggdl)
})

.subsetIndexGGDL <- function(se, index) {
    ## get the data range
    gpCoords <- .rowRangesFromList(se)
   
    l <- S4Vectors::metadata(index)
    l$chromosomes <- gpCoords
    minWidth <- min(width(gpCoords))
    if(minWidth < l$tileSize) {
        l$tileSize <- minWidth
        l$chunkSize <- minWidth - 2*l$overhangSize
    }
    indx <- .makeTiles(l)
    GenomeInfoDb::seqinfo(indx) <- GenomeInfoDb::seqinfo(gpCoords)
    
    return(indx)
}

## underlying function to subset by overlaps
.subsetByOverlapsGGDL <- function(query, subject, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {

    ## need to make even widths, otherwise subsetByOverlaps does it with a warning
    if(any((width(subject) %% 2) == 1)) {
        futile.logger::flog.info("Some subset ranges have odd widths. Rounding to the next even number.")
        idx <- which((width(subject) %% 2) == 1)
        width(subject)[idx] <- width(subject)[idx] + 1
    }          
    settings <- slot(query, "settings")
    design <- design(query)
    sf <- sizeFactors(query)
    h5 <- slot(query, "hdf5")

    ## iterate over all SummarizedExperiments and subset
    se <- lapply(query@data, function(se) {
        res <- subsetByOverlaps(se, subject, maxgap = maxgap,
                                minoverlap = minoverlap,
                                type=type, invert = invert, ...)
        return(res)
    })

    ## reduce the list to non-empty. However in order to mimic what happens
    ## if the complete subset is zero, we have to keep the result, as we might
    ## need it later. Otherwise it will be an empty list, which is not how
    ## SummarizedExperiment and GenoGAMDataSet deals with it. On the other hand
    ## we don't want to carry a bunch of empty SEs around if they are not needed anyway.
    reduced_se <- se[as.logical(vapply(se, length, integer(1)))]

    if(length(reduced_se) == 0) {
        reduced_se <- se[1]
        index <- GenomicRanges::GRanges()
    }
    else {
        reduced_se <- lapply(reduced_se, function(y) {
            ## set correct seqinfo
            GenomeInfoDb::seqlevels(rowRanges(y), pruning.mode="coarse") <- names(reduced_se)
            return(y)
        })
        slot(settings, "chromosomeList") <- GenomeInfoDb::seqlevels(reduced_se[[1]])
        index <- .subsetIndex(reduced_se, getIndex(query))
    }

    ## make id slot
    splitid <- .rowRangesFromList(reduced_se)
    splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))
    
    ggdl <- new("GenoGAMDataSetList", settings = settings,
               design = design, sizeFactors = sf, index = index,
               data = reduced_se, id = splitid, hdf5 = h5)
    
    return(ggdl)
}

#' @rdname GenoGAMDataSetList-subsetting
setMethod("subsetByOverlaps", signature(x = "GenoGAMDataSetList", ranges = "GRanges"),
          function(x, ranges, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {
              type <- match.arg(type)
              if(type == "any") {
                  maxgap <- -1L
                  minoverlap <- 0L
              }
              res <- .subsetByOverlapsGGDL(query = x, subject = ranges,
                                       maxgap = maxgap, minoverlap = minoverlap,
                                       type = type, invert = invert)
              return(res)
          })

#' @rdname GenoGAMDataSetList-subsetting
setMethod("[", c("GenoGAMDataSetList", "GRanges"), function(x, i) {
    ggdl <- subsetByOverlaps(x, i, type = "any", maxgap = 0L, minoverlap = 1L)
    return(ggdl)
})

## #' @rdname GenoGAMDataSetList-subsetting
## setMethod("[[", c("GenoGAMDataSetList", "numeric"), function(x, i) {
##     gr <- getIndex(x)[i]
##     ggd <- subsetByOverlaps(x,gr)
##     return(ggd)
## })

## ## Tile computation
## ## ================

## A function to enhance absoluteRanges function to values beyond 2^31
## by splitting the data and using own Coordinates class.
.absRanges <- function(x) {
    totalLen <- sum(as.numeric(GenomeInfoDb::seqlengths(x)))
    if(totalLen >= 2^31) {
        parts <- (totalLen %/% 2^31) + 1
        coords <- vector("list", parts)
        left <- GenomeInfoDb::seqlevels(x)
        for(ii in 1:parts) {
            if(length(left) == 0) {
                stop("No chromosomes left to compute coordinates")
            }
            idx <- getIndex(x)
            GenomeInfoDb::seqlevels(idx, pruning.mode = "coarse") <- left
            csum <- cumsum(as.numeric(GenomeInfoDb::seqlengths(idx)))
            id <- which(csum < 2^31)
            left <- GenomeInfoDb::seqlevels(idx)[-id]
            GenomeInfoDb::seqlevels(idx, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(idx)[id]
            coords[[ii]] <- as(GenomicRanges::absoluteRanges(idx), "Coordinates")
            if(ii > 1) {
                sh <- max(end(coords[[ii-1]]))
                end(coords[[ii]]) <- end(coords[[ii]]) + sh
                start(coords[[ii]]) <- start(coords[[ii]]) + sh
            }
        }
        coords <- do.call("rbind", coords)
    }
    else {
        coords <- GenomicRanges::absoluteRanges(getIndex(x))
        coords <- as(coords, "Coordinates")
    }
    return(coords)
}

#' Function to retrieve the row coordinates
#' @param x The GenoGAMDataSetList object
#' @return A Coordinates object specifying the row coordinates of
#' each tile
#' @noRd
.getCoordinatesGGDL <- function(x) {

    ## if genome is complete use the fast Bioconductor function
    totalLen <- sum(as.numeric(GenomeInfoDb::seqlengths(x)))
    if(totalLen == length(x)) {
        coords <- .absRanges(x)
    }
    ## otherwise the slower version 'by block'
    else {
        rr <- rowRanges(x)
        coords <- as(ranges(getIndex(x)), "Coordinates")
        current <- 1
        for(ii in 1:length(rr)) {
            ov <- IRanges::findOverlaps(rr[[ii]], getIndex(x))
            sh <- S4Vectors::subjectHits(ov)
            qh <- S4Vectors::queryHits(ov)
            l <- range(S4Vectors::splitAsList(qh, sh))
            l <- Coordinates(l[,1], l[,2])
            len <- length(l) + current - 1
            end(l) <- end(l) + max(end(coords)[current - 1], 0)
            start(l) <- start(l) + max(end(coords)[current - 1], 0)
            start(coords)[current:len] <- start(l)
            end(coords)[current:len] <- end(l)
            width(coords)[current:len] <- width(l)
            current <- len + 1
        }
    }
    return(coords)
}

#' Function to retrieve the row coordinates as a list
#' @param ind The GenoGAMDataSetList index
#' @param rr The GenoGAMDataSetList chromosome specific rowRanges
#' @return An integerList with the row numbers for each tile
#' @noRd
.getListedCoordinates <- function(ind, rr) {

    gr <- .extractGR(rr)
    ov <- IRanges::findOverlaps(rr, ind)
    sh <- Rle(S4Vectors::subjectHits(ov))
    qh <- queryHits(ov)
    l <- range(S4Vectors::splitAsList(qh, sh))
    l <- IRanges::IRanges(l[,1], l[,2])
    return(l)

}

#' compute metrics for each tile
#' @param x The GenoGAMDataSet object
#' @param what A character naming the metric
#' @param na.rm Should NAs be ignored
#' @return The metric value
#' @noRd
.MetricsFunGGDL <- function(x, what, na.rm = FALSE) {

    if(slot(x, "hdf5")) {
        res <- .MetricsFunHDF5(x, what, na.rm)
    }
    else {
        res <- .MetricsFunSimple(x, what, na.rm)
    }

    return(res)
}

#' compute metrics for each tile
#' @param x The GenoGAMDataSet object
#' @param what A character naming the metric
#' @param na.rm Should NAs be ignored
#' @return The metric value
#' @noRd
.MetricsFunSimple <- function(x, what, na.rm = FALSE) {

    ind <- getIndex(x)
    data <- assay(x)
    rr <- rowRanges(x)
                    
    res <- sapply(1:dim(x)[2], function(ii) {
        l <- unlist(lapply(names(data), function(d) {
            df <- data[[d]][[ii]]
            by <- ind[GenomeInfoDb::seqnames(ind) == d]
            by <- .getListedCoordinates(by, rr[[d]])
            rle <- IRanges::extractList(df, by)
            eval(call(what, rle, na.rm = na.rm))
        }))
        return(l)
    })
    
    if(!is(res, "matrix")) {
        if(is(res, "list") & is.null(dim(res))) {
            res <- matrix()
        }
        else {
            res <- t(matrix(res))
        }
    }

    colnames(res) <- colnames(x)
    rownames(res) <- getIndex(x)$id
    return(res)
}

#' compute metrics for each tile
#' @param x The GenoGAMDataSet object
#' @param what A character naming the metric
#' @param na.rm Should NAs be ignored
#' @return The metric value
#' @noRd
.MetricsFunHDF5 <- function(x, what, na.rm = FALSE) {

    ind <- getIndex(x)
    data <- assay(x)
    rr <- rowRanges(x)
    
    res <- sapply(1:dim(x)[2], function(ii) {
        l <- unlist(BiocParallel::bplapply(names(data), function(d) {
            df <- data[[d]][,ii]
            by <- ind[GenomeInfoDb::seqnames(ind) == d]
            by <- .getListedCoordinates(by, rr[[d]])
            rleDF <- as(df, "DataFrame")          
            rle <- IRanges::extractList(rleDF[,1], by)
            eval(call(what, rle, na.rm = na.rm))
        }))
        return(l)
    })
    
    if(!is(res, "matrix")) {
        if(is(res, "list") & is.null(dim(res))) {
            res <- matrix()
        }
        else {
            res <- t(matrix(res))
        }
    }

    colnames(res) <- colnames(x)
    rownames(res) <- getIndex(x)$id
    return(res)
}

#' Computing metrics
#'
#' Computing metrics on each tile of the GenoGAMDataSetList object.
#' All metrics from the Summary generics group, as well as
#' mean, var, sd, median, mad and IQR are supported.
#'
#' @param x A GenoGAMDataSetList object
#' @param ... Additional arguments
#' @param na.rm Should NAs be dropped. Otherwise the result is NA
#' @return A matrix with the specified metric computed per tile per column
#' of the assay data.
#' @examples
#' ggd <- makeTestGenoGAMDataSetList()
#' sum(ggd)
#' min(ggd)
#' max(ggd)
#' mean(ggd)
#' var(ggd)
#' sd(ggd)
#' median(ggd)
#' mad(ggd)
#' IQR(ggd)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenoGAMDataSetList-metrics
setMethod("Summary", "GenoGAMDataSetList", function(x, ..., na.rm = FALSE) {
    
    ind <- getIndex(x)
    data <- assay(x)
    rr <- rowRanges(x)
    
    if(slot(x, "hdf5")) {
        
        res <- sapply(1:dim(x)[2], function(ii) {
            l <- unlist(BiocParallel::bplapply(names(data), function(d) {
                df <- data[[d]][,ii]
                by <- ind[GenomeInfoDb::seqnames(ind) == d]
                by <- .getListedCoordinates(by, rr[[d]])
                rleDF <- as(df, "DataFrame")          
                rle <- IRanges::extractList(rleDF[,1], by)
                (getFunction(.Generic))(rle, na.rm = na.rm)
            }))
            return(l)
        })
    }
    else {

        res <- sapply(1:dim(x)[2], function(ii) {
            l <- unlist(lapply(names(data), function(d) {
                df <- data[[d]][[ii]]
                by <- ind[GenomeInfoDb::seqnames(ind) == d]
                by <- .getListedCoordinates(by, rr[[d]])
                rle <- IRanges::extractList(df, by)
                (getFunction(.Generic))(rle, na.rm = na.rm)
            }))
            return(l)
        })
    }

    
    if(!is(res, "matrix")) {
        if(is(res, "list") & is.null(dim(res))) {
            res <- matrix()
        }
        else {
            res <- t(matrix(res))
        }
    }

    colnames(res) <- colnames(x)
    rownames(res) <- getIndex(x)$id
    return(res)
})

#' @rdname GenoGAMDataSet-metrics
setMethod("mean", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "mean")
})

#' @rdname GenoGAMDataSet-metrics
setMethod("var", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "var")
})

#' @rdname GenoGAMDataSet-metrics
setMethod("sd", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "sd")
})

#' @rdname GenoGAMDataSet-metrics
setMethod("median", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "median")
})

#' @rdname GenoGAMDataSet-metrics
setMethod("mad", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "mad")
})

#' @rdname GenoGAMDataSet-metrics
setMethod("IQR", "GenoGAMDataSetList", function(x) {
    .MetricsFunGGDL(x, "IQR")
})

## Cosmetics
## =========

## Show method for GenomicTiles.
setMethod("show", "GenoGAMDataSetList", function(object) {
    .showGenoGAMDataSet(object)
})
