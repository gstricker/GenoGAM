#############################
## GenoGAM class 
##################

#' @include GenoGAMSettings-class.R
#' @include GenoGAMFamily-class.R
NULL

setClassUnion("HDF5OrMatrix", c("matrix", "HDF5Matrix"))

#' GenoGAM class
#'
#' This is the class the holds the complete model as well all hyperparameters
#' and settings that were used to fit it. It extends the RangedSummarizedExperiment
#' class by adding a couple of more slots to hold hyperparameters and settings.
#' The 'assays' slot holds the basepair fit and standard deviation. Additionally
#' all knot positions and beta coefficients will be stored in the 'smooths' slot
#' in order to be able to make use of the piecewise function that produces the
#' fit. For information on the slots inherited from SummarizedExperiment
#' check the respective class.
#' 
#' @slot family A GenoGAMFamily object
#' @slot design The formula of the model
#' @slot sizeFactors The offset used in the model. 
#' @slot factorialDesign The factorial design used. The same as colData in the
#' GenoGAMDataSet
#' @slot params All hyperparameters used to fit the data. The parameters
#' estimated by cross validation can also be found here. But the parameters
#' used in cross validation are in the settings slot.
#' @slot settings A GenoGAMSettings object representing the global
#' settings that were used to compute the model.
#' @slot coefs The coefficients of the knots
#' @slot knots The relative knot positions
#' @slot hdf5 A logical slot indicating if the data is stored in a HDF5 format on hard drive
#' @name GenoGAM-class
#' @rdname GenoGAM-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAM
setClass("GenoGAM",
         contains = "RangedSummarizedExperiment",
         slots = list(family = "GenoGAMFamily",
             design = "formula",
             sizeFactors = "numeric",
             factorialDesign = "DataFrame",
             params = "list",
             settings = "GenoGAMSettings",
             coefs = "HDF5OrMatrix",
             knots = "numeric",
             hdf5 = "logical"),
         prototype = prototype(family = GenoGAMFamily(),
             design = ~ s(x),
             sizeFactors = numeric(),
             factorialDesign = S4Vectors::DataFrame(),
             params = list(),
             settings = GenoGAMSettings(),
             coefs = matrix(),
             knots = numeric(),
             hdf5 = FALSE))

## Validity
## ========

.validateFamilyType <- function(object) {
    if(!is(slot(object, "family"), "GenoGAMFamily")) {
        return("'family' must be a GenoGAMFamily object")
    }
    NULL
}

.validateDesignType <- function(object) {
    if(!is(slot(object, "design"), "formula")) {
        return("'design' must be a formula object")
    }
    NULL
}

.validateSFType <- function(object) {
    if(!is(slot(object, "sizeFactors"), "numeric")) {
        return("'sizeFactors' must be a numeric object")
    }
    NULL
}

.validateFactorialDesign <- function(object) {
    if(!is(slot(object, "factorialDesign"), "DataFrame")) {
        return("'factorialDesign' must be a DataFrame class")
    }
}

.validateParamsType <- function(object) {
    if(!is(slot(object, "params"), "list")) {
        return("'params' must be a list object")
    }
    NULL
}

.validateSettingsType <- function(object) {
    if(!is(slot(object, "settings"), "GenoGAMSettings")) {
        return("'settings' must be a GenoGAMSettings object")
    }
    NULL
}

.validateCoefsType <- function(object) {
    if(!is(slot(object, "coefs"), "HDF5Matrix") &
        !is(slot(object, "coefs"), "matrix")){
        return("'coefs' must be either a HDF5Matrix or matrix object")
    }
    NULL
}

.validateGGKnotsType <- function(object) {
    if(!is(slot(object, "knots"), "numeric")){
        return("'knots' must be a numeric vector")
    }
    NULL
}

## general validate function
.validateGenoGAM <- function(object) {
    c(.validateFamilyType(object),
      .validateDesignType(object),
      .validateSFType(object),
      .validateFactorialDesign(object),
      .validateParamsType(object),
      .validateSettingsType(object),
      .validateCoefsType(object),
      .validateGGKnotsType(object),
      .validateH5Type(object))
}

S4Vectors::setValidity2("GenoGAM", .validateGenoGAM)


## Constructor
## ========

#' GenoGAM constructor
#'
#' The GenoGAM constructor, not designed to be actually used, by the user.
#' Rather to be a point of reference and documentation for slots and how
#' to access them.
#'
#' @aliases getSettings getFamily colData getParams getKnots getCoefs fits se dimnames colnames
#' @param object,x For use of S4 methods. The GenoGAM object.
#' @param i A GRanges object (only for subsetting)
#' @param ggd The initial GenoGAMDataSet object. Only needed if fromHDF5 is TRUE.
#' @param fromHDF5 A convenience argument to create a GenoGAM object from the already
#' computed fits that are stored as HDF5 files
#' @param split A logical argument indicating if the model was fitted in a
#' per-chromosome fashion or not. Only needed if fromHDF5 is TRUE.
#' @param ... Slots of the GenoGAM class. See the slot description.
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
#' @name GenoGAM
#' @rdname GenoGAM-class
#' @export
GenoGAM <- function(..., ggd = NULL, fromHDF5 = FALSE, split = FALSE) {
    if(fromHDF5){
        return(.GenoGAMFromHDF5(ggd = ggd, split = split, ...))
    }
    return(new("GenoGAM", ...))
}

## A function to create GenoGAM from HDF5
.GenoGAMFromHDF5 <- function(ggd, split, ...) {
    if(split) {
        gg <- .GenoGAMListFromHDF5(ggd = ggd, ...)
    }
    
    else {
        settings <- slot(ggd, "settings")
    
        ## get Identifier for the data
        path <- slot(settings, "hdf5Control")$dir
        ident <- .getIdentifier(path, fits = TRUE)

        rr <- rowRanges(ggd)

        ## read HDF5 file
        h5file <- file.path(path, paste0("fits_dataset", "_", ident))

        ## read HDF5 data and make SummarizedExperiment objects
        fits <- HDF5Array::HDF5Array(h5file, "fits")
        ses <- HDF5Array::HDF5Array(h5file, "ses")

        se <- SummarizedExperiment::SummarizedExperiment(rowRanges = rr,
                                                         assays = list(fits = h5fits,
                                                                       se = h5ses))

        ## generate knots and get coefficients
        bpknots <- slot(settings, "dataControl")$bpknots
        knots <- .generateKnotPositions(ggd, bpknots = 20)

        h5coefs <- file.path(path, paste0("coefs_", ident))
        coefs <- HDF5Array::HDF5Array(h5coefs, "coefs")
    
        gg <- new("GenoGAM", se, design = design(ggd),
               sizeFactors = sizeFactors(ggd), factorialDesign = colData(ggd),
               settings = settings, knots = knots[[1]], coefs = coefs,
               hdf5 = TRUE, ...)
    }


    return(gg)
}


## Accessors
## =========

#' @describeIn GenoGAM An accessor to the design slot
setMethod("design", "GenoGAM", function(object) {
    slot(object, "design")
})

#' @describeIn GenoGAM An accessor to the sizeFactors slot
setMethod("sizeFactors", "GenoGAM", function(object) {
    slot(object, "sizeFactors")
})

#' @export
setGeneric("getSettings", function(object) standardGeneric("getSettings"))

#' @describeIn GenoGAM An accessor to the settings slot
setMethod("getSettings", "GenoGAM", function(object) {
    slot(object, "settings")
})

#' @export
setGeneric("getFamily", function(object) standardGeneric("getFamily"))

#' @describeIn GenoGAM An accessor to the family slot
setMethod("getFamily", "GenoGAM", function(object) {
    slot(object, "family")
})

#' @describeIn GenoGAM An accessor to the factorialDesign slot.
setMethod("colData", "GenoGAM", function(x) {
    slot(x, "factorialDesign")
})

#' @export
setGeneric("getParams", function(object) standardGeneric("getParams"))

#' @describeIn GenoGAM An accessor to the params slot
setMethod("getParams", "GenoGAM", function(object) {
    slot(object, "params")
})

#' @export
setGeneric("getCoefs", function(object) standardGeneric("getCoefs"))

#' @describeIn GenoGAM An accessor to the coefs slot
setMethod("getCoefs", "GenoGAM", function(object) {
    slot(object, "coefs")
})

#' @export
setGeneric("getKnots", function(object) standardGeneric("getKnots"))

#' @describeIn GenoGAM An accessor to the knots slot
setMethod("getKnots", "GenoGAM", function(object) {
    slot(object, "knots")
})

#' @describeIn GenoGAM The accessor to the fits and standard errors
setMethod("assay", c("GenoGAM", "missing"), function(x, i) {
    res <- assays(x)
    if(length(res) < 1) {
        return(NULL)
    }
    return(res[[1]])
})

#' @export
setGeneric("fits", function(object) standardGeneric("fits"))

#' @describeIn GenoGAM An accessor to the fits
setMethod("fits", "GenoGAM", function(object) {
    assays(object)[["fits"]]
})

#' @export
setGeneric("se", function(object) standardGeneric("se"))

#' @describeIn GenoGAM An accessor to the standard errors
setMethod("se", "GenoGAM", function(object) {
    assays(object)[["se"]]
})

#' @export
setGeneric("pvalue", function(object) standardGeneric("pvalue"))

#' @describeIn GenoGAM An accessor to the pvalues
setMethod("pvalue", "GenoGAM", function(object) {
    assays(object)[["pvalue"]]
})

#' @describeIn GenoGAM column names of GenoGAM
setMethod("colnames", "GenoGAM", function(x) {
    rownames(slot(x, "colData"))
})

#' @describeIn GenoGAM The names of the dimensions of GenoGAM
setMethod("dimnames", "GenoGAM", function(x) {
    list(names(x), rownames(slot(x, "colData")))
})

#' @export
setGeneric("is.HDF5", function(object) standardGeneric("is.HDF5"))

#' @describeIn GenoGAM A boolean function that is true if object uses HDF5 backend
setMethod("is.HDF5", signature(object = "GenoGAM"), function(object) {
    res <- slot(object, "hdf5")
    return(res)
})


## Test GenoGAM
## =============

#' Make an example /code{GenoGAM}
#'
#' @return A /code{GenoGAM} object
#' @examples
#' gg <- makeTestGenoGAM()
#' @export
makeTestGenoGAM <- function() {

    k <- 10000
    ip <- sin(seq(-7, 5, length.out = k)) + 1
    reduce <- 100^(seq(1, 0.01, length.out = k))
    reduce <- reduce/max(reduce)
    background <- (sin(seq(-7, 5, length.out = k)) + 1) * reduce
    sdx <- runif(k)
    sdexperiment <- runif(k)
    gr <- GenomicRanges::GPos(GenomicRanges::GRanges("chrXYZ", IRanges::IRanges(1, k)))
    GenomeInfoDb::seqlengths(gr) <- 1e6
    sdf <- S4Vectors::DataFrame(input = background, IP = ip)
    sedf <- S4Vectors::DataFrame(sdx = sdx, sde = sdexperiment)
    names(sdf) <- names(sedf) <- c("s(x)", "s(x):experiment")
    
    se <- SummarizedExperiment(rowRanges = gr, assays = list(fits = sdf, se = sedf))
    family <- GenoGAMFamily()
    design <- ~ s(x) + s(x, by = experiment)
    sizeFactors <- c(input = 0, IP = 0)
    params <- list(cv = FALSE)
        
    gg <- GenoGAM(se, family = family, design = design,
                   sizeFactors = sizeFactors, params = params)
    coldf <- S4Vectors::DataFrame(experiment = c(0, 1))
    rownames(coldf) <- c("input", "IP")
    slot(gg, "factorialDesign") <- coldf
    
    return(gg)
}


## Subsetting 
## ==========

#' @describeIn GenoGAM Additional subsetting by single brackets
setMethod("[", c("GenoGAM", "GRanges"), function(x, i) {
    gg <- subsetByOverlaps(x, i)
    slot(gg, "settings")@chromosomeList <- GenomeInfoDb::seqlevels(i)
    return(gg)
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

setMethod("show", "GenoGAM", function(object) {
    .showGenoGAM(object)
})

## View the dataset
## 
## Cbinding the columns all together and coercing to data.frame
## 
## @param object A \code{GenoGAM} object
## @param ranges A \code{GRanges} object. Makes it possible to
## select regions by \code{GRanges}. Either ranges or seqnames, start and
## end must be supplied
## @param seqnames A chromosomes name. Either ranges or seqnames, start and
## end must be supplied
## @param start A start site. Either ranges or seqnames, start and
## end must be supplied
## @param end An end site. Either ranges or seqnames, start and
## end must be supplied
## @return A data.frame of the selected data.
## @examples
## gg <- makeTestGenoGAM()
## gr <- GenomicRanges::GRanges("chrI", IRanges(1,40))
## head(view(gg, gr))
## @author Georg Stricker \email{georg.stricker@@in.tum.de}
## @rdname GenoGAM-view
## @export
## setMethod("view", "GenoGAM", function(object, ranges = NULL, seqnames = NULL,
##                                       start = NULL, end = NULL) {
##     ## keep "seqnames" for consistency with Bioc, but rename variable as subset in
##     ## .subsetByPosition does not work properly otherwise
##     chromosome <- seqnames 
##     if(is.null(seqnames) & is.null(start) & is.null(end) & is.null(ranges)) {
##         temp <- object
##     }
##     else {
##         if(!is.null(ranges)) {
##             temp <- .subsetByRanges(object, ranges)
##         }
##         else {
##             if(is.null(start)) {
##                 start <- 1
##             }
##             if(is.null(end)) {
##                 end <- Inf
##             }
##             if(is.null(seqnames)){
##                 temp <- .subsetByPosition(object, pos >= start & pos <= end)
##             }
##             else {
##                 temp <- .subsetByPosition(object, seqnames == chromosome & pos >= start & pos <= end)
##             }
##         }
##     }
##     res <- cbind(slot(temp, "positions"), slot(temp, "fits"))
##     return(res)
## })
