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
                      formula = "formula", offset = "numeric", 
                      family = "character", response = "numeric",
                      fits = "numeric"),
         prototype = list(params = list(lambda = 0, theta = 0, H = 0,
                                        order = 2, penorder = 2),
                          knots = list(), designMatrix = new("dgCMatrix"),
                          beta = matrix(), vcov = new("dgCMatrix"),
                          penaltyMatrix = new("dgCMatrix"), formula = ~1,
                          offset = numeric(), family = "nb", 
                          response = numeric(), fits = numeric()))

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

.validateResponseType <- function(object) {
    if(class(slot(object, "response")) != "numeric") {
        return("'response' must be a numeric object")
    }
    NULL
}

.validateFitsType <- function(object) {
    if(class(slot(object, "fits")) != "numeric") {
        return("'fits' must be a numeric object")
    }
    NULL
}

## general validate function
.validateGenoGAMSetup <- function(object) {
    c(.validateParamsType(object), .validateKnotsType(object),
      .validateDesignMatrixType(object), .validateBetaType(object),
      .validateCovarianceType(object), .validatePenaltyMatrixType(object),
      .validateFormulaType(object), .validateOffsetType(object),
      .validateFamilyType(object), .validateResponseType(object),
      .validateFitsType(object))
}

setValidity2("GenoGAMSetup", .validateGenoGAMSetup)

#' Constructor
#' @noRd
GenoGAMSetup <- function(...) {
  return(new("GenoGAMSetup", ...))
}

#' Constructor function
#' @noRd
setupGenoGAM <- function(ggd, lambda = NULL, theta = NULL, H = 0, family = "nb", bpknots = 20, order = 2, penorder = 2) {

  ## knot placement
  knots <- generateKnotPositions(ggd, bpknots)

  ## take the first tile and see how many knots are in there. then add 8 additional knots (4 to each side)
  nbetas <- length(which(knots[[1]] > start(getIndex(ggd)[1]) & knots[[1]] < end(getIndex(ggd)[1]))) + 2*order
  nknots <- nbetas + 4

  ## formula handling
  formula <- design(ggd)

    ggsetup <- GenoGAMSetup(params = list(lambda = lambda, theta = theta, H = H,
                                          order = order, penorder = penorder),
                            knots = knots, formula = formula, 
                            offset = sizeFactors(ggd), family = family)
  
  return(ggsetup)
}

#' A function to generate knots genome-wide for P-Splines,
#' ensuring same spacing everywhere independent of chromosome length.
#' @noRd
generateKnotPositions <- function(ggd, bpknots = 20){
  referenceChrom <- names(which(seqlengths(ggd) == max(seqlengths(ggd))))
  positions <- rowRanges(ggd)[seqnames(rowRanges(ggd)) == referenceChrom]
  nknots <- round(length(positions)/bpknots)
  x <- pos(positions)
  knots <- placeKnots(x = x, nknots = nknots)

  ## chroms <- seqlengths(ggd) 
  ## res <- lapply(names(chroms), function(chr) {
  ##   idx <- which(knots < chroms[chr])
  ##   ans <- knots[c(idx, c(1:4 + length(idx)))]
  ##   return(ans)
  ## })
  ## names(res) <- names(chroms)
  ## return(res)
  res <- list(knots)
  names(res) <- referenceChrom
  return(res)
}

#' A function to place knots for P-Splines
#' Courtesy to Simon Wood (mgcv). Slightly changed.
#' @noRd
placeKnots <- function(x, nknots, ord = 2) {
  m <- ord + 1
  nk <- nknots - ord
  xu <- max(x)
  xl <- min(x)
  xr <- xu - xl
  multFactor <- min(1/(10^floor(log10(xr))), 0.001)
  xl <- xl - xr * multFactor
  xu <- xu + xr * multFactor
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

#' A function to build the identity matrix I with multiple epsilon
#' Courtesy to Simon Wood (mgcv)
#' @noRd
buildIMatrix <- function(p, epsilon) {
  I = Matrix(diag(p), sparse = TRUE) ##initialize a diagonal identity matrix
  return(epsilon*I)
}


## Test
## ======

## library(GenoGAM)
## library(Matrix)
## FOLDER <- "/s/project/coreProm/Michi/thornton/align_STAR"
## CONFIG <- "~/workspace/analysis/diffBinding/config.txt"
## config <- data.table::fread(CONFIG)
## BiocParallel::register(BiocParallel::SnowParam(workers=4))

## BPPK <- 20
## CHUNKSIZE <- BPPK*1000
## OV <- BPPK*10

## ggd <- GenoGAMDataSet(CONFIG, chunkSize = CHUNKSIZE,
##                       overhangSize = OV, design = ~ s(x) + s(x, by = genotype),
##                       directory = FOLDER)
## ggd <- computeSizeFactors(ggd)
## ggs <- setupGenoGAM(ggd)
