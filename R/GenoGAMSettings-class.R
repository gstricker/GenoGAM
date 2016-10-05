#######################
## GenoGAMSettings
## ====================

setClassUnion("biocParallel", c("MulticoreParam", "SnowParam", "SerialParam", "BatchJobsParam"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("logicalOrNULL", c("logical", "NULL"))

#' GenoGAMSettings
#'
#' This class is designed to store settings for the computation of the
#' GenoGAM package
#'
#' @slot center A logical or NULL value to specify if the raw data should be
#' centered, i.e. only the midpoint of the fragment will be used to
#' represent its coverage. See details.
#' @slot chromosomeList A character vector of chromosomes to be used.
#' NULL for all chromosomes.
#' @slot bamParams An object of class ScanBamParam.
#' See ?Rsamtools::ScanBamParam.
#' @slot parallel A parallel backend of the respective class. See
#' BiocParalell for the options
#' @slot processFunction A custom function on how to process raw data. Not
#' used if center is TRUE/FALSE.
#' @slot optimMethod The optiomisation method to be used in cross validation.
#' @slot optimControl Settings for the optim() function.
#' @details Center can have three values: TRUE, FALSE, NULL. TRUE will
#' trigger the center function, FALSE will trigger the use of the entire
#' fragment. NULL should be used in case a custom process function is used.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("GenoGAMSettings",
         slots = list(center = "logicalOrNULL", chromosomeList = "characterOrNULL",
             bamParams = "ScanBamParam", parallel = "biocParallel",
             processFunction = "function", optimMethod = "character",
             optimControl = "list"),
         prototype = list(center = TRUE, chromosomeList = NULL,
             bamParams = Rsamtools::ScanBamParam(what = c("pos", "qwidth")),
             parallel = BiocParallel::registered()[[1]],
             processFunction = identity,
             optimMethod = "BFGS",
             optimControl = list(maxit=100, fnscale=-1, trace = 10)))

## Validity
## ========

#' Validating the correct type

.validateBAMParamsType <- function(object) {
    if(class(object@bamParams) != "ScanBamParam") {
        return("'bamParams' must be a ScanBamParam object")
    }
    NULL
}

.validateProcessType <- function(object) {
    if(class(object@processFunction) != "function") {
        return("'processFunction' must be a function object")
    }
    NULL
}

.validateOptimMethodType <- function(object) {
    if(class(object@optimMethod) != "character") {
        return("'optimMethod' must be a character object")
    }
    NULL
}

.validateOptimControlType <- function(object) {
    if(class(object@optimControl) != "list") {
        return("'optimControl' must be a list object")
    }
    NULL
}

## general validate function
.validateGenoGAMSettings <- function(object) {
    c(.validateBAMParamsType(object), .validateProcessType(object),
      .validateOptimMethodType(object), .validateOptimControlType(object))
}

setValidity2("GenoGAMSettings", .validateGenoGAMSettings)

## Constructor
## ==========

#' The constructor function for GenoGAMSettings
#' 
#' The constructor function for GenoGAMSettings
#'
#' @param ... Any parameters corresponding to the slots and their possible
#' values. See \linkS4class{GenoGAMSettings}
#' @return A GenoGAMSettings object.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
GenoGAMSettings <- function(...) {
    return(new("GenoGAMSettings", ...)) 
}

## Cosmetics
## =========

## show method for GenoGAMSettings
.showGenoGAMSettings <- function(gtd) {
    ident <- identity
    cat("-------------------- Read-in parameters -----------------\n")
    cat("center:", gtd@center, "\n")
    cat("chromosomes:", paste(gtd@chromosomeList, collapse = ", "), "\n")
    cat("process Function:", ifelse(identical(gtd@processFunction, ident),
                                    "identity", "custom"), "\n")
    cat("\n")
    cat("-------------------- BAM parameters ---------------------\n")
    show(gtd@bamParams)
    cat("\n")
    cat("-------------------- Parallel backend -------------------\n")
    show(gtd@parallel)
    cat("\n")
    cat("-------------------- Cross Validation parameters --------\n")
    cat("Optimization method:", gtd@optimMethod, "\n")
    cat("Optimization control:\n")
    for(ii in 1:length(gtd@optimControl)) {
        cat(paste0("  ", names(gtd@optimControl)[ii], ": ",
                   gtd@optimControl[[ii]], "\n"))
    }
}

setMethod("show", "GenoGAMSettings", function(object) {
    .showGenoGAMSettings(object)
})

## Getter and Setter
## =================

setGeneric("getDefaults", function(object, ...) standardGeneric("getDefaults"))

setMethod("getDefaults", "GenoGAMSettings",
          function(object, what = NULL) {
              if(is.null(what)) return(object)
              else return(slot(object, what))
          })

setGeneric("setDefaults", function(object, ...) {
    standardGeneric("setDefaults")
})

setMethod("setDefaults", "GenoGAMSettings",
          function(object, ...) {
              for(param in names(list(...))) {
                  slot(object, param) <- list(...)[[param]]
                  if(param == "parallel") {
                      BiocParallel::register(slot(object, "parallel"))
                  }
              }
              return(object)
          })
