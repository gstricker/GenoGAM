## ======================================================
## Family class for specific distributions for GenoGAM
## ======================================================

setClassUnion("FunctionOrNULL", c("function", "NULL"))

#' GenoGAMFamily class
#'
#' This class holds the distribution family with the log-likelihood functions
#' and the first two derivatives, the gradient vector and the Hessian matrix.
#' It is not meant to be used directly, but rather for development convenience.
#'
#' @slot ll The log-likelihood function. Gives a scalar
#' @slot gradient The first derivative of the log-likelihood functions. Gives a vector.
#' @slot hessian An integer indicating the family. Negative Binomial = 1, Quasi-Binomial = 2.
#' @name GenoGAMFamily-class
#' @rdname GenoGAMFamily-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("GenoGAMFamily",
         slots = list(ll = "FunctionOrNULL",
                      gradient = "FunctionOrNULL",
                      hessian = "integer",
                      name = "character"),
         prototype = list(ll = NULL,
                          gradient = NULL,
                          hessian = 0L,
                          name = NA_character_))

## ===========
## Validation
## ===========

.validateLogLik <- function(object) {
    if(!is(slot(object, "ll"), "function")) {
        return("'ll' must be a function object")
    }
    NULL
}

.validateGradient <- function(object) {
    if(!is(slot(object, "gradient"), "function")) {
        return("'gradient' must be a function object")
    }
    NULL
}

.validateHessian <- function(object) {
    if(!is(slot(object, "hessian"), "integer")) {
        return("'hessian' must be an integer object")
    }
    NULL
}

.validateName <- function(object) {
    if(!is(slot(object, "name"), "character")) {
        return("'name' must be a character object")
    }
    NULL
}

.validateGenoGAMFamily <- function(object){
    c(.validateLogLik(object),
      .validateGradient(object),
      .validateHessian(object),
      .validateName(object))
}

S4Vectors::setValidity2("GenoGAMFamily", .validateGenoGAMFamily)

## Constructor
GenoGAMFamily <- function(...) {
    return(new("GenoGAMFamily", ...))
}

