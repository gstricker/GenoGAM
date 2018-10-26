## ===============
## GenoGAMSettings
## ===============

#' @include helper.R
NULL

setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("logicalOrNULL", c("logical", "NULL"))
setClassUnion("functionOrNULL", c("function", "NULL"))

#' GenoGAMSettings
#'
#' This class is designed to store global settings for the computation of the
#' GenoGAM package
#'
#' @param ... Any parameters corresponding to the slots and their possible
#' values.
#' 
#' @slot center A logical or NULL value to specify if the raw data should be
#' centered, i.e. only the midpoint of the fragment will be used to
#' represent its coverage. See details.
#' @slot chromosomeList A character vector of chromosomes to be used.
#' NULL for all chromosomes.
#' @slot bamParams An object of class ScanBamParam.
#' See ?Rsamtools::ScanBamParam for possible settings. Usually used to set
#' specific ranges, to read in.
#' @slot processFunction A custom function on how to process raw data. Not
#' used if center is TRUE/FALSE. This is not intended for the user, but if
#' needed anyway, see details.
#' @slot optimMethod The optiomisation method to be used in cross validation.
#' See ?optim for a complete list.
#' @slot optimControl List of control settings for the optim function.
#' Almost all parameters are supported, with a couple of exceptions. See
#' details. For a complete list of parameters see ?optim.
#' @slot estimControl List of control settings for the parameter estimation algorithm.
#' @slot hdf5Control List of control settings for the HDF5 backend
#' @slot dataControl List of control settings for the processed data.
#' The size of the region to use for the computation of the
#' count matrix, that is later used by DESeq2. Also the size of the regions that
#' will be used for Cross Validation. And the spacing between knots.
#' @details Center can have three values: TRUE, FALSE, NULL. TRUE will
#' trigger the center function, FALSE will trigger the use of the entire
#' fragment. NULL should be used in case a custom process function is used.
#' In case a custom function is used, it has to satisfy the following:
#' It has to handle a GAlignments object as input and output a 
#' GRanges object of regions, e.g. fragments. This regions are in turn used
#' to compute the coverage via the IRanges::coverage function. Note, that there
#' is a difference between the GAlignments object in the single and paired end
#' case.
#'
#' optimControl has two maxit fields: 'maxit' refers to the maximal iterations
#' in the cross validation procedure. Convergence rarely exceeds 50. 'betaMaxit'
#' refers to the maximum iterations in estimation of the beta parameters.
#' Convergence will sometimes need up to 1000 iterations. Also the parameter
#' 'trace' is overwritten. Please use the threshold setting through
#' futile.logger::flog() to control trace information. All other parameters can
#' be used as specified in ?optim. Note that the method used to estimate the
#' beta vector is 'L-BFGS-B'. The method for cross validation can be changed.
#' Keep in mind however, that not all methods will work, as it is not gradient
#' based.
#' @name GenoGAMSettings-class
#' @rdname GenoGAMSettings-class
#' @examples
#' # Construct the class
#' GenoGAMSettings()
#' 
#' # Construct the class with custom parameters
#' ## specify chromosomes
#' center <- FALSE
#' chromosomeList <- c('chr1', 'chr2')
#' GenoGAMSettings(center = center, chromosomeList = chromosomeList)
#' 
#' ## Specify ranges
#' gr <- GenomicRanges::GRanges("chr1", IRanges(1, 10000))
#' bamParams <- Rsamtools::ScanBamParam(which = gr)
#' GenoGAMSettings(bamParams = bamParams, center = TRUE)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @exportClass GenoGAMSettings
setClass("GenoGAMSettings",
         slots = list(center = "logicalOrNULL", chromosomeList = "characterOrNULL",
             bamParams = "ScanBamParam", processFunction = "functionOrNULL", 
             optimMethod = "character", optimControl = "list", estimControl = "list",
             hdf5Control = "list", dataControl = "list"),
         prototype = list(center = TRUE, chromosomeList = NULL,
             bamParams = Rsamtools::ScanBamParam(what = c("pos", "qwidth")),
             processFunction = NULL, optimMethod = "Nelder-Mead",
             optimControl = list(maxit = 60L, fnscale = -1L, trace = 1L),
             estimControl = list(eps = 1e-6, maxiter = 1000L),
             hdf5Control = list(dir = NULL, level = NULL),
             dataControl = list(regionSize = 4000L, bpknots = 20L)))

## Validity
## ========

## Validating the correct type

.validateBAMParamsType <- function(object) {
    if(!is(object@bamParams, "ScanBamParam")) {
        return("'bamParams' must be a ScanBamParam object")
    }
    NULL
}

.validateOptimMethodType <- function(object) {
    if(!is(object@optimMethod, "character")) {
        return("'optimMethod' must be a character object")
    }
    NULL
}

.validateOptimControlType <- function(object) {
    if(!is(object@optimControl, "list")) {
        return("'optimControl' must be a list object")
    }
    NULL
}

.validateEstimControlType <- function(object) {
    if(!is(object@estimControl, "list")) {
        return("'estimControl' must be a list object")
    }
    NULL
}

.validateHDF5ControlType <- function(object) {
    if(!is(object@hdf5Control, "list")) {
        return("'hdf5Control' must be a list object")
    }
    NULL
}

.validateDataControlType <- function(object) {
    if(!is(object@dataControl, "list")) {
        return("'dataControl' must be an integer object")
    }
    NULL
}


## general validate function
.validateGenoGAMSettings <- function(object) {
    c(.validateBAMParamsType(object), .validateOptimMethodType(object),
      .validateOptimControlType(object) ,.validateEstimControlType(object),
      .validateHDF5ControlType(object), .validateDataControlType(object))
}

S4Vectors::setValidity2("GenoGAMSettings", .validateGenoGAMSettings)

## Constructor
## ==========

#' The Constructor function for GenoGAMSettings is merely a wrapper for
#' new("GenoGAMSettings", ...)
#' 
#' @return An object of class GenoGAMSettings
#' @name GenoGAMSettings
#' @rdname GenoGAMSettings-class
#' @export
GenoGAMSettings <- function(...) {
    ggs <- new("GenoGAMSettings", ...)

    ## check if all optim params are there
    params <- slot(ggs, "optimControl")
    optimControl = list(maxit = 50, fnscale = -1, trace = 1)
    params <- .fillParameters(l = params, optimControl)
    slot(ggs, "optimControl") <- params

    ## check if all estimation algo params are there
    params <- slot(ggs, "estimControl")
    estimControl = list(eps = 1e-6, maxiter = 1000L)
    params <- .fillParameters(l = params, estimControl)
    slot(ggs, "estimControl") <- params

    ## check if all data params are there
    params <- slot(ggs, "dataControl")
    dataControl = list(regionSize = 4000L, bpknots = 20L)
    params <- .fillParameters(l = params, dataControl)
    slot(ggs, "dataControl") <- params

    ## initialize HDF5 options globally
    slot(ggs, "hdf5Control") <- .initializeHDF5Params(slot(ggs, "hdf5Control"))
    return(ggs)
}

## initializing HDF5 control parameters
.initializeHDF5Params <- function(params) {
    currentDir <- HDF5Array::getHDF5DumpDir()
    currentLevel <- HDF5Array::getHDF5DumpCompressionLevel()

    ## set or get current dump directory
    if(is.null(params$dir)) {
        params$dir <- currentDir
    }
    else {
        if(currentDir != params$dir) {
            HDF5Array::setHDF5DumpDir(params$dir)
        }
    }

    ## set or get current level
    if(is.null(params$level)) {
        params$level <- currentLevel
    }
    else {
        if(currentLevel != params$level) {
            HDF5Array::setHDF5DumpCompressionLevel(params$level)
        }
    }

    return(params)
}

## Cosmetics
## =========

## show method for GenoGAMSettings
.showGenoGAMSettings <- function(ggs) {
    if(!is.null(ggs@center)) {
        processFun <- FALSE
    }
    else {
        processFun <- TRUE
    }
    
    cat("-------------------- Read-in parameters -----------------\n")
    cat("Center:", ggs@center, "\n")
    cat("Chromosomes:", paste(ggs@chromosomeList, collapse = ", "), "\n")
    cat("Custom process function:", ifelse(processFun, "enabled,", "disabled"), 
        "\n")
    cat("\n")
    cat("-------------------- BAM parameters ---------------------\n")
    show(ggs@bamParams)
    cat("\n")
    cat("-------------------- Parallel backend -------------------\n")
    show(BiocParallel::registered()[[1]])
    cat("\n")
    cat("-------------------- Logger settings --------------------\n")
    cat("Logger threshold:", futile.logger::flog.threshold(), "\n")
    cat("\n")
    cat("-------------------- HDF5 settings ----------------------\n")
    cat("HDF5 Directory:", ggs@hdf5Control$dir, "\n")
    cat("HDF5 Compression Level (0-9):", ggs@hdf5Control$level, "\n")
    cat("\n")
    cat("-------------------- Data settings ----------------------\n")
    cat("Region size:", ggs@dataControl$regionSize, "\n")
    cat("\n")
    cat("Basepairs per Knot:", ggs@dataControl$bpknots, "\n")
    cat("\n")
    cat("-------------------- Optimization parameters ------------\n")
    cat("Optimization method:", ggs@optimMethod, "\n")
    cat("Optimization control:\n")
    for(ii in 1:length(ggs@optimControl)) {
        cat(paste0("  ", names(ggs@optimControl)[ii], ": ",
                   ggs@optimControl[[ii]], "\n"))
    }
    cat("Parameter estimation control:\n")
    for(ii in 1:length(ggs@estimControl)) {
        cat(paste0("  ", names(ggs@estimControl)[ii], ": ",
                   ggs@estimControl[[ii]], "\n"))
    }
}

setMethod("show", "GenoGAMSettings", function(object) {
    .showGenoGAMSettings(object)
})
