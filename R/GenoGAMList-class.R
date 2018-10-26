#############################
## GenoGAMList class 
############################

#' @include GenoGAMSettings-class.R
#' @include GenoGAM-class.R
NULL

#' GenoGAMList class
#'
#' This is the class the holds the complete model as well all hyperparameters
#' and settings that were used to fit it. It is basically identical to the
#' GenoGAM class, except for the data being inside a list of RangedSummarizedExperiment
#' objects.
#' 
#' @slot family The name of the distribution family used
#' @slot design The formula of the model
#' @slot sizeFactors The offset used in the model. 
#' @slot factorialDesign The factorial design used. The same as colData in the
#' GenoGAMDataSet
#' @slot params All hyperparameters used to fit the data. The parameters
#' estimated by cross validation can also be found here. But the parameters
#' used in cross validation are in the settings slot.
#' @slot settings A GenoGAMSettings object representing the global
#' settings that were used to compute the model.
#' @slot data A list of RangedSummarizedExperiment that holds the actual data
#' @slot id A GRanges object keeping the identifiers assigning the regions to the
#' respective list elements
#' @slot coefs The coefficients of the knots
#' @slot knots The relative knot positions
#' @name GenoGAMList-class
#' @rdname GenoGAMList-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("GenoGAMList",
         slots = list(family = "GenoGAMFamily",
             design = "formula",
             sizeFactors = "numeric",
             factorialDesign = "DataFrame",
             params = "list",
             settings = "GenoGAMSettings",
             data = "list", id = "GRanges",
             coefs = "HDF5OrMatrix",
             knots = "numeric",
             hdf5 = "logical"),
         prototype = prototype(family = GenoGAMFamily(),
             design = ~ s(x),
             sizeFactors = numeric(),
             factorialDesign = S4Vectors::DataFrame(),
             params = list(),
             settings = GenoGAMSettings(),
             data = list(), id = GenomicRanges::GRanges(),
             coefs = matrix(), knots = numeric(),
             hdf5 = FALSE))

## Validity
## ========

## general validate function
.validateGenoGAMList <- function(object) { 
    c(.validateFamilyType(object), ## all from GenoGAM-class.R
      .validateDesignType(object),
      .validateSFType(object),
      .validateFactorialDesign(object),
      .validateParamsType(object),
      .validateSettingsType(object),
      .validateCoefsType(object),
      .validateGGKnotsType(object),
      .validateDataType(object), ## from GenoGAMDataSetList-class.R
      .validateIDType(object), ## from GenoGAMDataSetList-class.R
      .validateH5Type(object)) ## from GenoGAMDataSet-class.R
}

S4Vectors::setValidity2("GenoGAMList", .validateGenoGAMList)


## Constructor
## ========

#' GenoGAMList constructor
#'
#' The GenoGAMList constructor, not designed to be actually used, by the user.
#' Rather to be a point of reference and documentation for slots and how
#' to access them.
#'
#' @param object,x For use of S4 methods. The GenoGAMList object.
#' @param withDimnames For use of S4 methods. The GenoGAMList object.
#' @param i A GRanges object (only for subsetting)
#' @param ... Slots of the GenoGAM class. See the slot description.
#' @param ggd The initial GenoGAMDataSet object. Only needed if fromHDF5 is TRUE.
#' @param fromHDF5 A convenience argument to create a GenoGAM object from the already
#' computed fits that are stored as HDF5 files
#' @return An object of the type GenoGAM.
#' @examples
#'
#' ## creating test GenoGAM object
#' gg <- makeTestGenoGAM()
#' gg
#'
#' ## using accessors
#' design(gg)
#' sizeFactors(gg)
#' getSettings(gg)
#' getFamily(gg)
#' colData(gg)
#' getParams(gg)
#' getCoefs(gg)
#' getKnots(gg)
#' rowRanges(gg)
#' assay(gg)
#' assays(gg)
#' fits(gg) 
#' se(gg)
#' @name GenoGAMList
#' @rdname GenoGAMList-class
GenoGAMList <- function(..., ggd = NULL, fromHDF5 = FALSE) {
    if(fromHDF5) {
        return(.GenoGAMListFromHDF5(ggd = ggd, ...))
    }
    return(new("GenoGAMList", ...))
}

## A function to create GenoGAMList from HDF5
.GenoGAMListFromHDF5 <- function(ggd, ...) {
    settings <- slot(ggd, "settings")
    
    ## get Identifier for the data
    path <- slot(settings, "hdf5Control")$dir
    ident <- .getIdentifier(path, fits = TRUE)

    ## finally build object
    rr <- rowRanges(ggd)
    splitid <- dataRange(ggd)
    splitid$id <- as.character(GenomeInfoDb::seqnames(splitid))
    
    selist <- lapply(splitid$id, function(id) {
        
        ## read HDF5 file
        h5file <- file.path(path, paste0("fits_", id, "_", ident))
        
        ## read HDF5 data and make SummarizedExperiment objects
        h5fits <- HDF5Array::HDF5Array(h5file, "fits")
        h5ses <- HDF5Array::HDF5Array(h5file, "ses")
        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rr[[id]],
                                                         assays = list(fits = h5fits,
                                                                       se = h5ses))
        return(se)
    })

    names(selist) <- splitid$id

    ## generate knots and get coefficients
    bpknots <- slot(settings, "dataControl")$bpknots
    knots <- .generateKnotPositions(ggd, bpknots = 20)

    h5coefs <- file.path(path, paste0("coefs_", ident))
    coefs <- HDF5Array::HDF5Array(h5coefs, "coefs")

    args <- list(...)
    if("params" %in% names(args)) {
        params <- c(args$params, chromosomes = granges(splitid))
        args$params <- NULL
    }
    else {
        params <- list(chromosomes = granges(splitid))
    }
    
    ggl <- do.call(new, c(list("GenoGAMList", design = design(ggd),
               sizeFactors = sizeFactors(ggd), factorialDesign = colData(ggd),
               settings = settings, data = selist, id = splitid,
               knots = knots[[1]], coefs = coefs, params = params,
               hdf5 = TRUE), args))

    return(ggl)
}


## Accessors
## =========

#' @describeIn GenoGAMList Get the dimension of the object
setMethod("dim", "GenoGAMList", function(x) {
    dims <- sapply(x@data, dim)
    if(length(dims) == 0) {
        return(c(0, 0))
    }
    res <- c(sum(as.numeric(dims[1,])), max(dims[2,]))
    return(res)
})

#' @describeIn GenoGAMList The length of the object
setMethod("length", "GenoGAMList", function(x) {
    dim(x)[1]
})

#' @describeIn GenoGAMList The seqlengths of the object
setMethod("seqlengths", "GenoGAMList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(GenomeInfoDb::seqlengths(rowRanges(x)[[1]]))
    }

    integer()
})

#' @describeIn GenoGAMList The seqlevels of the object
setMethod("seqlevels", "GenoGAMList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(GenomeInfoDb::seqlevels(rowRanges(x)[[1]]))
    }

    character()
})

#' @describeIn GenoGAMList The seqlevelsInUse of the object
setMethod("seqlevelsInUse", "GenoGAMList", function(x) {
    if(length(rowRanges(x)) > 0) {
        return(names(rowRanges(x)))
    }

    character()
})

#' @describeIn GenoGAMList An accessor to the design slot
setMethod("design", "GenoGAMList", function(object) {
    slot(object, "design")
})

#' @describeIn GenoGAMList An accessor to the sizeFactors slot
setMethod("sizeFactors", "GenoGAMList", function(object) {
    slot(object, "sizeFactors")
})

#' @describeIn GenoGAMList An accessor to the settings slot
setMethod("getSettings", "GenoGAMList", function(object) {
    slot(object, "settings")
})

#' @describeIn GenoGAMList An accessor to the family slot
setMethod("getFamily", "GenoGAMList", function(object) {
    slot(object, "family")
})

#' @describeIn GenoGAMList An accessor to the factorialDesign slot.
setMethod("colData", "GenoGAMList", function(x) {
    slot(x, "factorialDesign")
})

#' @describeIn GenoGAMList An accessor to the params slot
setMethod("getParams", "GenoGAMList", function(object) {
    slot(object, "params")
})

#' @describeIn GenoGAMList An accessor to the coefs slot
setMethod("getCoefs", "GenoGAMList", function(object) {
    slot(object, "coefs")
})

#' @describeIn GenoGAMList An accessor to the knots slot
setMethod("getKnots", "GenoGAMList", function(object) {
    slot(object, "knots")
})


#' @describeIn GenoGAMList The accessor to the fits and standard errors
setMethod("assay", c("GenoGAMList", "missing"), function(x, i) {
    lapply(x@data, assay)
})

#' @describeIn GenoGAMList get a list of list of assays from the
#' GenoGAMList object. Just for completeness, shouldn't be needed.
setMethod("assays", "GenoGAMList", function(x, ...) {
    lapply(x@data, assays)
})

#' @describeIn GenoGAMList get a list of rowRanges from the
#' GenoGAMList object
setMethod("rowRanges", "GenoGAMList", function(x, ...) {
    lapply(x@data, rowRanges)
})

#' @describeIn GenoGAMList An accessor to the fits
setMethod("fits", "GenoGAMList", function(object) {
    lapply(assays(object), function(y) y[["fits"]])
})

#' @describeIn GenoGAMList An accessor to the standard errors
setMethod("se", "GenoGAMList", function(object) {
    lapply(assays(object), function(y) y[["se"]])
})

#' @describeIn GenoGAMList An accessor to the pvalues
setMethod("pvalue", "GenoGAMList", function(object) {
    lapply(assays(object), function(y) y[["pvalue"]])
})

#' @describeIn GenoGAMList column names of GenoGAMList
setMethod("colnames", "GenoGAMList", function(x) {
    if(length(x@data) == 0) {
        return(NULL)
    }
    rownames(slot(x@data[[1]], "colData"))
})

#' @describeIn GenoGAMList The names of the dimensions of GenoGAMList
setMethod("dimnames", "GenoGAMList", function(x) {
    list(names(x@data[[1]]), colnames(x))
})

#' @describeIn GenoGAMList A boolean function that is true if object uses HDF5 backend
setMethod("is.HDF5", signature(object = "GenoGAMList"), function(object) {
    res <- slot(object, "hdf5")
    return(res)
})

## Test GenoGAM
## =============

#' Make an example /code{GenoGAMList}
#'
#' @return A /code{GenoGAMList} object
#' @examples
#' ggl <- makeTestGenoGAMList()
#' @export
makeTestGenoGAMList <- function() {

    k <- 10000
    ip <- sin(seq(-7, 5, length.out = k)) + 1
    reduce <- 100^(seq(1, 0.01, length.out = k))
    reduce <- reduce/max(reduce)
    background <- (sin(seq(-7, 5, length.out = k)) + 1) * reduce
    chroms <- c("chrX", "chrY", "chrZ")
    gr <- GenomicRanges::GPos(GenomicRanges::GRanges(chroms,
                                                     IRanges::IRanges(c(1,1,1), c(k, k, k))))
    GenomeInfoDb::seqlengths(gr) <- c(1e6, 2e6, 1.5e6)

    selist <- lapply(chroms, function(chr) {
        add <- sample(0:3, 2)
        ip <- ip + add[1]
        background <- background + add[2]
        sdx <- runif(k)
        sdexperiment <- runif(k)
        sdf <- S4Vectors::DataFrame(input = background, IP = ip)
        sedf <- S4Vectors::DataFrame(sdx = sdx, sde = sdexperiment)
        names(sdf) <- names(sedf) <- c("s(x)", "s(x):experiment")
        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = gr[GenomeInfoDb::seqnames(gr) == chr,],
                                                         assays = list(fits = sdf, se = sedf))
    })
    names(selist) <- chroms
    id <- .extractGR(gr)
    id$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(gr)))
    
    
    
    family <- GenoGAMFamily()
    design <- ~ s(x) + s(x, by = experiment)
    sizeFactors <- c(input = 0, IP = 0)
    params <- list(cv = FALSE)
    coldf <- S4Vectors::DataFrame(experiment = c(0, 1))
    rownames(coldf) <- c("input", "IP")
        
    ggl <- GenoGAMList(family = family, design = design,
                      sizeFactors = sizeFactors, params = params,
                      data = selist, id = id, factorialDesign = coldf)
        
    return(ggl)
}


## Subsetting 
## ==========

#' Subset method for GenoGAMList
#'
#' @details
#' Those are various methods to subset the GenoGAMList object.
#' By logical statement or GRanges overlap. The '[' subsetter is
#' just a short version of 'subsetByOverlaps'.
#'
#' @param x A GenoGAMList object.
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
#' @return A subsetted GenoGAMList object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname GenoGAMList-subsetting
setMethod("subset", "GenoGAMList", function(x, ...) {
    if(all(dim(x) == c(0, 0))) return(x)

    fam <- getFamily(x)
    settings <- getSettings(x)
    design <- design(x)
    sf <- sizeFactors(x)
    cd <- colData(x)
    params <- getParams(x)
    coefs <- getCoefs(x)
    knots <- getKnots(x)
    hdf5 <- is.HDF5(x)
    
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
    }

    ## make id slot
    splitid <- .rowRangesFromList(reduced_se)
    splitid$id <- as.character(S4Vectors::runValue(GenomeInfoDb::seqnames(splitid)))

    ggl <- new("GenoGAMList", settings = settings,
               design = design, sizeFactors = sf,
               family = fam, factorialDesign = cd,
               params = params, coefs = coefs,
               knots = knots, hdf5 = hdf5,
               data = reduced_se, id = splitid)

    return(ggl)
})

## underlying function to subset by overlaps
.subsetByOverlapsGGL <- function(query, subject, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {

    ## need to make even widths, otherwise subsetByOverlaps does it with a warning
    if(any((width(subject) %% 2) == 1)) {
        futile.logger::flog.debug("Some subset ranges have odd widths. Rounding to the next even number.")
        idx <- which((width(subject) %% 2) == 1)
        width(subject)[idx] <- width(subject)[idx] + 1
    }
    fam <- getFamily(query)
    settings <- getSettings(query)
    design <- design(query)
    sf <- sizeFactors(query)
    cd <- colData(query)
    params <- getParams(query)
    coefs <- getCoefs(query)
    knots <- getKnots(query)
    hdf5 <- is.HDF5(query)

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
    }

    ## make id slot
    splitid <- .rowRangesFromList(reduced_se)
    splitid$id <- as.character(GenomeInfoDb::seqnames(splitid))
    
    ggl <- new("GenoGAMList", settings = settings,
               design = design, sizeFactors = sf,
               family = fam, factorialDesign = cd,
               params = params, coefs = coefs,
               knots = knots, hdf5 = hdf5,
               data = reduced_se, id = splitid)
    
    return(ggl)
}

#' @rdname GenoGAMList-subsetting
setMethod("subsetByOverlaps", signature(x = "GenoGAMList", ranges = "GRanges"),
          function(x, ranges, maxgap = -1L, minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   invert = FALSE, ...) {
              type <- match.arg(type)
              if(type == "any") {
                  maxgap <- -1L
                  minoverlap <- 0L
              }
              res <- .subsetByOverlapsGGL(query = x, subject = ranges,
                                       maxgap = maxgap, minoverlap = minoverlap,
                                       type = type, invert = invert)
              return(res)
          })

##' @rdname GenoGAMList-subsetting
setMethod("[", c("GenoGAMList", "GRanges"), function(x, i) {
    ggl <- subsetByOverlaps(x, i)
    return(ggl)
})



## Cosmetics
## ==========

## returns NA if object is NULL and the object otherwise
.check <- function(x) {
    if(is.null(x)) {
        res <- NA
    }
    else {
        res <- x
    }
    return(res)
}

## The actual show function
.showGenoGAM <- function(gg) {
    params <- slot(gg, "params")
    
    if(is.null(params$cv)) {
        params$cv <- FALSE
    }
    
    form <- design(gg)
    cl <- class(gg)
    dims <- dim(gg)
    md <- unique(names(colData(gg)))
    tracks <- colnames(gg)
    samples <- rownames(colData(gg))
    sf <- sizeFactors(gg)
    fam <- getFamily(gg)
    
    cat("Family: ", fam@name, "\n")
    cat("Formula: ", paste(as.character(form), collapse = " "), "\n")

    tsize <- params$tileSize
    tname <- "Tiles"
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

    chroms <- getSettings(gg)@chromosomeList
    tnum <- params$numTiles
    
    cat("Class:", cl, "\n")
    cat("Dimension:", dims, "\n")
    cat(paste0("Samples(", length(samples), "):"), samples, "\n")
    cat(paste0("Design variables(", length(md), "):"), md, "\n")
    cat(paste0("Smooth functions(", length(tracks), "):"), tracks, "\n")
    cat("Chromosomes:", chroms, "\n")
    
    cat("\n")
    cat("Size factors:\n")
    show(sizeFactors(gg))

    if(params$cv) {
        performed <- "Performed"
    }
    else {
        performed <- "Not performed"
    }
    cat("\n")
    cat("Cross Validation: ", performed, "\n")
    if(params$cv) {
        cat("  Number of folds: ", .check(params$kfolds), "\n")
        cat("  Interval size: ", .check(params$intervalSize), "\n")
        cat("  Number of regions: ", .check(params$regions), "\n")
    }
    cat("\n")
    cat("Spline Parameters:\n")
    cat("  Knot spacing: ", .check(params$bpknots), "\n")
    cat("  B-spline order: ", .check(params$order), "\n")
    cat("  Penalization order: ", .check(params$penorder), "\n")

    cat("\n")
    cat("Tile settings:\n")
    cat("  Chunk size:", .check(params$chunkSize), "\n")
    cat("  Tile size:", .check(params$tileSize), "\n")
    cat("  Overhang size:", .check(params$overhangSize), "\n")
    cat("  Number of tiles:", .check(params$numTiles), "\n")
    cat("  Evaluated genome ranges:\n")
    show(.check(params$chromosomes))
}

setMethod("show", "GenoGAMList", function(object) {
    .showGenoGAM(object)
})
