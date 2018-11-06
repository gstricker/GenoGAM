#################
## Setup class ##
#################

#' @include helper.R
#' @include GenoGAMFamily-class.R
NULL

#' GenoGAMSEtup class
#'
#' A class to embody the setup for a GenoGAM fit
#' 
#' @slot params A list of all hyper-parameters which are either estimated
#' or fixed. At the moment the smoothing parameter lambda, the overdispersion
#' parameter theta the order and penalization order of the splines as well as
#' the second regularization parameter epsilon are provided.
#' @slot knots A list of knot positions on each chromosome.
#' @slot designMatrix The design matrix.
#' @slot beta The vector of coefficients to be estimated. Initialized.
#' @slot vcov The covariance matrix
#' @slot penaltyMatrix The penalty matrix S with penalization order r. 
#' By default r = 2.
#' @slot formula The formula of the model. Usually the same as the design of
#' the GenoGAMDataSet
#' @slot design The actual design used to model the data, obtained
#' from merging the formula into the colData
#' @slot offset An offset of the samples
#' @slot family The distribution to be used. At the moment only "nb"
#' (Negative Binomial) is available.
#' @slot response The response vector
#' @slot fits The vector of fits
#' @slot control A list of parameters to control the parameter estimation
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
setClass("GenoGAMSetup",
         slots = list(params = "list", knots = "list",
                      designMatrix = "dgCMatrix", beta = "matrix",
                      se = "list", penaltyMatrix = "dgCMatrix",
                      formula = "formula", design = "matrix",
                      offset = "numeric", family = "GenoGAMFamily",
                      response = "numeric", fits = "list",
                      control = "list"),
         prototype = list(params = list(lambda = 0, theta = 0, eps = 0,
                                        order = 2, penorder = 2),
                          knots = list(), designMatrix = new("dgCMatrix"),
                          beta = matrix(,0,0), se = list(),
                          penaltyMatrix = new("dgCMatrix"), formula = ~1,
                          design = matrix(,0,0),
                          offset = numeric(), family = GenoGAMFamily(), 
                          response = integer(), fits = list(),
                          control = list(eps = 1e-6, maxiter = 1000, alpha = 1,
                                         rho = 0.5, c = 1e-4, m = 6)))

## Validity
## ========

## Validating the correct type
.validateParamsType <- function(object) {
    if(!is(slot(object, "params"), "list")) {
        return("'params' must be a list object")
    }
    NULL
}

.validateParamsElements <- function(object) {
    if(!all(names(slot(object, "params")) %in% c("lambda", "theta", "eps", "order", "penorder"))) {
        return("'params' must contain the elements 'lambda', 'theta', 'eps', 'order' and 'penorder'")
    }
    NULL
}

.validateGGSKnotsType <- function(object) {
    if(!is(slot(object, "knots"), "list")) {
        return("'knots' must be a list object")
    }
    NULL
}

.validateDesignMatrixType <- function(object) {
    if(!is(slot(object, "designMatrix"), "dgCMatrix")) {
        return("'designMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateBetaType <- function(object) {
    if(!is(slot(object, "beta"), "matrix")) {
        return("'beta' must be a matrix object")
    }
    NULL
}

.validateSEType <- function(object) {
    if(!is(slot(object, "se"), "list")) {
        return("'se' must be a list object")
    }
    NULL
}

.validatePenaltyMatrixType <- function(object) {
    if(!is(slot(object, "penaltyMatrix"), "dgCMatrix")) {
        return("'penaltyMatrix' must be a dgCMatrix object")
    }
    NULL
}

.validateFormulaType <- function(object) {
    if(!is(slot(object, "formula"), "formula")) {
        return("'formula' must be a formula object")
    }
    NULL
}

.validateGGSDesignType <- function(object) {
    if(!is(slot(object, "design"), "matrix")) {
        return("'design' must be a matrix object")
    }
    NULL
}

.validateOffsetType <- function(object) {
    if(mode(slot(object, "offset")) != "numeric") {
        return("'offset' must be a numeric object")
    }
    NULL
}

.validateFamilyType <- function(object) {
    if(!is(slot(object, "family"), "GenoGAMFamily")) {
        return("'family' must be a GenoGAMFamily object")
    }
    NULL
}

.validateResponseType <- function(object) {
    if(mode(slot(object, "response")) != "numeric")  {
        return("'response' must be either a numeric or an Rle object")
    }
    NULL
}

.validateFitsType <- function(object) {
    if(mode(slot(object, "fits")) != "list") {
        return("'fits' must be a list object")
    }
    NULL
}

.validateControlType <- function(object) {
    if(!all(names(slot(object, "control")) %in% c("eps", "maxiter", "alpha", "rho", "c", "m"))) {
        return("'control' must contain the elements 'eps', 'maxiter', 'alpha', 'rho', 'c', 'm'")
    }
    NULL
}

## general validate function
.validateGenoGAMSetup <- function(object) {
    c(.validateParamsType(object), .validateGGSKnotsType(object),
      .validateDesignMatrixType(object), .validateBetaType(object),
      .validateSEType(object), .validatePenaltyMatrixType(object),
      .validateFormulaType(object), .validateGGSDesignType(object),
      .validateOffsetType(object), .validateFamilyType(object),
      .validateResponseType(object), .validateFitsType(object),
      .validateParamsElements(object), .validateControlType(object))
}

S4Vectors::setValidity2("GenoGAMSetup", .validateGenoGAMSetup)

## Constructor
GenoGAMSetup <- function(...) {
    ggs <- new("GenoGAMSetup", ...)
    params <- slot(ggs, "params")
    coreValues <- c(lambda = 0, theta = 0, eps = 0, order = 2,
                    penorder = 2)
    params <- .fillParameters(l = params, coreValues)
    slot(ggs, "params") <- params

    ## check if all estimation algo params are there
    params <- slot(ggs, "control")
    estimControl = list(eps = 1e-6, maxiter = 1000, alpha = 1, rho = 0.5,
                        c = 1e-4, m = 6)
    params <- .fillParameters(l = params, estimControl)
    slot(ggs, "control") <- params
    return(ggs)
}

## the dimension function
## @param x A GenoGAMSetup object
## @return The four dimensions of the object (designMatrix rows, designMatrix
## columns, experiment design rows, experiment design columns)
## @noRd
setMethod("dim", "GenoGAMSetup", function(x) {
    Xdim <- dim(slot(x, "designMatrix"))
    designDim <- dim(slot(x, "design"))
    blockDim <- c(Xdim[1]/max(1, designDim[1]), Xdim[2]/max(1, designDim[2]))
    return(c(blockDim, designDim))
})

## the length function
## @param x A GenoGAMSetup object
## @return The length of the object as the product of all dimensions
## @noRd
setMethod("length", "GenoGAMSetup", function(x) {
    return(prod(dim(x)))
})

## Get number of functions from GenoGAMSetup
.nfun <- function(ggs) {
    vars <- .getVars(slot(ggs, "formula"), type = "covar")
    return(length(vars))
}

## Get number of betas from GenoGAMSetup
.nbeta <- function(ggs) {
    betas <- slot(ggs, "beta")
    if(ncol(betas) > 1) {
        res <- nrow(betas)
    }
    else {
        funs <- .nfun(ggs)
        res <- nrow(betas)/funs
    }
    return(res)
}


## Constructor function
setupGenoGAM <- function(ggd, lambda = NULL, theta = NULL, eps = 0, family = "nb",
                         bpknots = 20, order = 2, penorder = 2, control = list()) {

    ## knot placement    
    positions <- ranges(getIndex(ggd))[1]
    x <- start(positions):end(positions)
    nknots <- round(length(x)/bpknots)
    knots <- .placeKnots(x = x, nknots = nknots)

    X <- .buildDesignMatrix(knots = knots, pos = x, order = order)
    des <- .getDesignFromFormula(design(ggd), colData(ggd))
    ## Number of betas = number of knots
    ## Number of functions = Count the functions in the formula
    nbetas <- nknots
    nfun <- length(.getVars(design(ggd), type = "covar"))
    S <- .buildSMatrix(nbetas, penorder)
    I <- .buildIMatrix(nbetas, eps)
    S <- S + I

    ## turn knots into list to comply with object requirements
    knots <- list(knots)
    names(knots) <- "1"

    offset <- rep(sizeFactors(ggd), each = getTileSize(ggd))

    if(family == "nb") {
        fam <- GenoGAMFamily(ll = ll_pen_nb,
                             gradient = gr_ll_pen_nb,
                             hessian = 1L,
                             name = "nb")
    }
    else {
        fam <- GenoGAMFamily()
    }

    ggsetup <- GenoGAMSetup(params = list(lambda = lambda, theta = theta, eps = eps,
                                          order = order, penorder = penorder),
                            knots = knots, formula = design(ggd),
                            design = des, offset = offset, family = fam,
                            designMatrix = X, penaltyMatrix = S, control = control)
  
    return(ggsetup)
}

## A function to generate knots for P-Splines from a GenoGAMDataSet object
.generateKnotPositions <- function(ggd, bpknots = 20){

    positions <- IRanges::ranges(getIndex(ggd))[1]
    x <- IRanges::start(positions):IRanges::end(positions)
    nknots <- round(length(x)/bpknots)
    knots <- .placeKnots(x = x, nknots = nknots)

    res <- list(knots)
    names(res) <- "1"
    return(res)
}

## A function to place knots for P-Splines
## Courtesy to Simon Wood (mgcv). Slightly changed.
.placeKnots <- function(x, nknots, ord = 2) {
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

## A function to build the penalization matrix S
## Courtesy to Simon Wood (mgcv). Slightly changed
.buildSMatrix <- function(p, order) {

    ##initialize a diagonal identity matrix
    S <- Matrix::bandSparse(p, k = 0, diag = c(list(rep(1, p))))
    ## S <- Matrix::Matrix(diag(p), sparse = TRUE) 
    for (i in 1:order) {
        S <- Matrix::diff(S) ## twice the difference   
    }
    S <- Matrix::t(S)%*%S ## square
    return(S)
}

## A function to build the identity matrix I with multiple epsilon
## Courtesy to Simon Wood (mgcv)
.buildIMatrix <- function(p, epsilon) {

    ##initialize a diagonal identity matrix
    I <- Matrix::bandSparse(p, k = 0, diag = c(list(rep(1, p)))) 
    return(epsilon*I)
}

## B spline basis
.bspline <- function(x, k, ord = 2, derivative = 0) {
  res <- splines::spline.des(k, x, ord + 2, rep(derivative,length(x)), sparse=TRUE)$design
  return(res)
}

## build a block matrix from a template submatrix and a design matrix
.blockMatrixFromDesignMatrix <- function(template, design) {
  ## ## create 4-dim array by 'inserting' the template into the desing matrix
  ## arr <- array(template, c(dim(template),dim(design)))
  ## dims <- dim(arr)
  ## multP <- c(3,4,1,2)
  ## reduceP <- c(3,1,4,2)
  ## ## permute array for correct multiplication
  ## multArr <- aperm(arr, multP)*as.vector(design)
  ## ## permute array for correct reduction
  ## reducedArr <- aperm(multArr, reduceP)
  ## ## reduce 4-dim array to 2-dim matrix
  ## dim(reducedArr) <- c(nrow(template)*nrow(design), ncol(template)*ncol(design))
  ## return(reducedArr)

    ## use kronecker product for matrices instead of own function
    ## very memory efficient
    res <- design %x% template
    return(res)
}

## Build design matrix from the data
.buildDesignMatrix <- function(knots, pos, order) {

    ## build matrix
    X <- as(.bspline(pos, knots, order),"dgCMatrix")
    ## design <- .getDesignFromFormula(design(ggd), colData(ggd))
    ## X <- as(.blockMatrixFromDesignMatrix(x, design), "dgCMatrix")
    return(X)
}

## get the design from formula
.getDesignFromFormula <- function(formula, design) {
    formulaCols <- .getVars(formula)
    designCols <- as.vector(na.omit(formulaCols))
    newDesign <- as.matrix(design[,designCols])
    colnames(newDesign) <- designCols

    if("s(x)" %in% names(formulaCols)) {
        control <- rep(1, nrow(newDesign))
        newDesign <- cbind(control, newDesign)
    }
    return(newDesign)
}
