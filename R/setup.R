## TODO:
## - Do we need the covariance slot?

#################
## Setup class ##
#################

#' GenoGAMSEtup class
#'
#' A class to embody the setup for a GenoGAM fit
#' 
#' @slot params A list of all hyper-parameters which are either estimated
#' or fixed.
#' At the moment the smoothing parameter lambda and a second regularization
#' matrix H are provided.
#' @slot knots A list of knot positions on each chromosome.
#' @slot designMatrix The design matrix.
#' @slot beta The vector of coefficients to be estimated. Initialized.
#' @slot vcov The covariance matrix
#' @slot penaltyMatrix The penalty matrix S with penalization order r. 
#' By default r = 2.
#' @slot formula The formula of the model. Usually the same as the design of
#' the GenoGAMDataSet
#' @slot offset An offset of the samples
#' @slot family The distribution to be used. At the moment only "nb"
#' (Negative Binomial) is available.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
setClass("GenoGAMSetup",
         slots = list(params = "list", knots = "list",
                      designMatrix = "dgCMatrix", beta = "matrix",
                      vcov = "dgCMatrix", penaltyMatrix = "dgCMatrix", 
                      formula = "formula", offset = "numeric", family = "character"),
         prototype = list(params = list(lambda = 0, H = matrix()),
                          knots = list(), designMatrix = new("dgCMatrix"),
                          beta = matrix(), vcov = new("dgCMatrix"),
                          penaltyMatrix = new("dgCMatrix"), formula = ~1,
                          offset = numeric(), family = "nb"))

## Validity
## ========

#' Validating the correct type
.validateParamsType <- function(object) {
    if(class(slot(object, "params")) != "list") {
        return("'params' must be a list object")
    }
    NULL
}

.validateKnotsType <- function(object) {
    if(class(slot(object, "knots")) != "list") {
        return("'knots' must be a list object")
    }
    NULL
}

.validateDesignMatrixType <- function(object) {
    if(class(slot(object, "designMatrix")) != "dgCMatrix") {
        return("'designMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateBetaType <- function(object) {
    if(class(slot(object, "beta")) != "matrix") {
        return("'beta' must be a matrix object")
    }
    NULL
}

.validateCovarianceType <- function(object) {
    if(class(slot(object, "vcov")) != "dgCMatrix") {
        return("'vcov' must be a dgCMatrix object")
    }
    NULL
}

.validatePenaltyMatrixType <- function(object) {
    if(class(slot(object, "penaltyMatrix")) != "dgCMatrix") {
        return("'penaltyMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateFormulaType <- function(object) {
    if(class(slot(object, "formula")) != "formula") {
        return("'formula' must be a formula object")
    }
    NULL
}

.validateOffsetType <- function(object) {
    if(class(slot(object, "offset")) != "numeric") {
        return("'offset' must be a numeric object")
    }
    NULL
}

.validateFamilyType <- function(object) {
    if(class(slot(object, "family")) != "character") {
        return("'family' must be a character object")
    }
    NULL
}

## general validate function
.validateGenoGAMSetup <- function(object) {
    c(.validateParamsType(object), .validateKnotsType(object),
      .validateDesignMatrixType(object), .validateBetaType(object),
      .validateCovarianceType(object), .validatePenaltyMatrixType(object),
      .validateFormulaType(object), .validateOffsetType(object),
      .validateFamilyType(object))
}

setValidity2("GenoGAMSetup", .validateGenoGAMSetup)

#' Constructor
#' @noRd
GenoGAMSetup <- function(...) {
  return(new("GenoGAMSetup", ...))
}

#' Constructor function
#' @noRd
setupGenoGAM <- function(ggd, lambda = NULL, H = NULL, family = "nb", bpknots = 20, order = 2, penorder = 2) {

  ## knot placement
  chroms <- seqlevelsInUse(ggd)
  knots <- lapply(chroms, function(chr) {
    positions <- rowRanges(ggd)[seqnames(rowRanges(ggd)) == chr]
    nknots <- round(length(positions)/bpknots)
    x <- pos(positions)
    knots <- placeKnots(x = x, nknots = nknots)
    return(knots)
  })
  names(knots) <- chroms
  ## take the first tile and see how many knots are in there. then add 8 additional knots (4 to each side)
  nbetas <- length(which(knots[[1]] > start(getIndex(ggd)[1]) & knots[[1]] < end(getIndex(ggd)[1]))) + 2*order
  nknots <- nbetas + 4

  ## formula handling
  formula <- design(ggd)

  ## How many splines?
  vars <- .getVars(formula, "covar")
  nsplines <- length(vars)

  ## param initialization
  if(is.null(lambda)) lambda <- 0
  if(is.null(H)) {
    H <- matrix(0, nbetas*nsplines, nbetas*nsplines)
  }
  
  ## penaltyMatrix
  S <- buildSMatrix(nbetas, penorder)

  ggsetup <- GenoGAMSetup(params = list(lambda = lambda, H = H),
                          knots = knots, penaltyMatrix = S,
                          formula = formula, offset = sizeFactors(ggd),
                          family = family)
  
  return(ggsetup)
}

## Accessor function
getSlot <- function(object, slot = c("params", "knots", "X",
                                     "beta", "vcov", "S", "formula",
                                     "offset", "family")) {
  slot <- match.arg(slot)
  
  res <- switch(slot, 
                params = object@params,
                knots = object@knots,
                X = object@designMatrix,
                beta = object@beta,
                vcov = object@vcov,
                S = object@penaltyMatrix,
                formula = object@formula,
                offset = object@offset,
                family = object@family
                )
  return(res)
}

#' A function to place knots for P-Splines
#' Courtesy to Simon Wood (mgcv)
#' @noRd
placeKnots <- function(x, nknots, ord = 2) {
  m <- ord + 1
  nk <- nknots - ord
  xu <- max(x)
  xl <- min(x)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(nk - 1)
  k <- seq(xl - dx * m, xu + dx * m, length.out = nk + 2 * ord + 2)
  
  return(k)
}

#' A function to build the penalization matrix S
#' Courtesy to Simon Wood (mgcv)
#' @noRd
buildSMatrix <- function(p, order) {
  S = Matrix(diag(p), sparse = TRUE) ##initialize a diagonal identity matrix
  for (i in 1:order) S = diff(S) ## twice the difference
  S = t(S)%*%S ## square
}
