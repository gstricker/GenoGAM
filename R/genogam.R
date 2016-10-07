## =========================
##  The main fitting method
## =========================

#' Function to initialize the Cross Validation
#'
#' @param formula A formula object.
#' @param data A DataFrame object.
#' @param family A distribution family object.
#' @param sp The penalization parameters, as many as there are splines.
#' @param fixedPars The parameters to be fixed.
#' @param fixedPars A logical vector indicating if a parameter is fixed or not.
#' @param sf A size factor vector
#' @return A list with the the parameters and their initial values.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
.initCV <- function(formula, data, family, sp, fixedPars, experimentDesign, sf) {

    idx <- sample(length(data), 1)
    initDataset <- .meltGTile(list(data[[idx]]), experimentDesign, sf, formula)

    ## The tryCatch not very elegant, but for now necessary to workaround
    ## the mgcv routine to estimate a first guess in low count regions.
    ## Otherwise the while loop in mgcv runs til divison by zero (see function mgcv::initial.sp)
    init <- tryCatch({
        mgcv::gam(formula, initDataset[[1]], family = family, sp = sp)
    }, error = function(e) {
        c(lambda = 100, theta = 5)
    })
    if(class(init)[1] == "numeric") {
        lambda <- init['lambda']
        theta <- init['theta']
    }
    else {
        theta <- init$family$getTheta()
        lambda <- init$full.sp
    }
    if(is.null(lambda)) lambda <- mean(init$sp)

    if(fixedPars[1]) {
        pars <- c(theta = theta)
        fpars <- list(lambda = lambda, theta = NULL)
    }
    if(fixedPars[2]) {
        pars <- c(lambda = log(lambda))
        fpars <- list(lambda = NULL, theta = exp(theta))
    }
    if(all(!fixedPars)) {
        pars <- c(lambda = log(lambda), theta = theta)
        fpars <- list(lambda = NULL, theta = NULL)
    }
    return(list(pars = pars, fpars = fpars))
}

#' genogam
#'
#' This is the fitting function for GenoGAMDataSet. It processes the data
#' in GenoGAMDataSet, estimates the overdispersion and the penalization parameter,
#' passes all information to mgcv::gam in parallel fashion and extracts
#' the results. So far the model is restricted to Negativ Binomial distribution (mgcv::nb()).
#' 
#' @param ggd A GenoGAMDataSet object to be fitted.
#' @param lambda The penalization parameter. Will be estimated if missing.
#' @param family A distribution family object. So far only mgcv::nb() is allowed.
#' @param bpknots Number of basepairs per one knot, that is,  how dense should the knots
#' be placed. The denser the knots, the more sensitive the fit. Note however, that computation
#' time increases approximately cubic with every additional knot.
#' @param kfolds An integer number giving the number of k-folds to be used in cross
#' validation, if parameters need to be estimated.
#' @param intervallSize The size of the intervalls to be used in cross validation.
#' Short intervalls are used instead of single points to be left out due to
#' spatial correlation. If replicates are present it is advised to make them bigger,
#' e.g. 2*fragment size. Otherwise, depending on the density of the data, they should
#' not exceed the size of a short read.
#' @param m The penalization order of the P-Splines.
#' @return A GenoGAM object containing the fits and parameters.
#' @examples
#' \dontrun{
#' ## simple example
#' config <- data.frame(ID = c("input", "IP"),
#'                      file = c("myInput.bam",
#'                               "myIP.bam"),,
#'                     paired = c(FALSE, FALSE),
#'                     type = c(0,1), stringsAsFactors = FALSE)
#' bpk <- 100 ## basepairs per one knot
#' chunkSize <- 5000
#' overhang <- round(7*chunkSize/bpk) ##overhang with 7 knots
#' knots <- chunkSize/bpk
#' ## build the GenoGAMDataSet
#' gtiles <- GenoGAMDataSet(config = config, chunkSize = chunkSize, overhangSize = overhang,
#'                          design = ~ s(x) + s(x, by = type))
#' gtiles <- computeSizeFactors(gtiles)
#' fits <- genogam(gtiles, bpknots = bpk)
#' }
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @rdname fitGenoGAM
#' @export
genogam <- function(ggd, lambda = NULL, family = mgcv::nb(),
                    bpknots = 20, kfolds = 10, intervallSize = 20, m = 2) {

  settings <- slot(ggd, "settings")
  chunkCoords <- getChunkIndex(ggd)
  tileIndex <- getIndex(ggd)
  tileSettings <- tileSettings(ggd)
  check <- checkSettings(ggd)
  if(!check) break
    
  knots <- round(tileSettings$tileSize/bpknots)
  if(knots > 100) {
    warning("Number of knots is above the recommended maximum of 100. Computation time will be high. Decrease tile size for faster computation.")
  }

  if(family$n.theta == 0) theta <- family$getTheta(trans = TRUE)
  else theta <- NULL
  FIXLAMBDA <- ifelse(is.null(lambda), FALSE, TRUE)
  FIXTHETA <- ifelse(is.null(theta), FALSE, TRUE)
  cv <- FALSE
    
  formula <- design(ggd)
  formula <- .updateFormula(formula, knots, m)
  
  ## How many splines?
  vars <- .getVars(formula)[-1]
  nsplines <- length(vars)

  ## convert data
  futile.logger::flog.info("Process data")
  data <- getTile(ggd)
  vec <- Rle(rep(tileIndex$id, width(tileIndex)))
  data@unlistData$id <- vec

  ncv <- min(20, length(data))

  if(is.null(lambda) | is.null(theta)) {
    futile.logger::flog.info("Estimating parameters")
    ## initialize parameters for CV
    sp <- rep(lambda, nsplines)
    init <- .initCV(formula, data, family, sp, c(FIXLAMBDA, FIXTHETA), colData(ggd), sizeFactors(ggd))
    if(!FIXLAMBDA) {
      init$pars['lambda'] <- min(init$pars['lambda'], log(knots))
    }

    ## get the tile ids for CV
    sumMatrix <- sum(ggd)
    if(ncv >= 20) {
      pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix[[1]], names(sizeFactors(ggd)))))
      ids <- order(pvals)[1:ncv]
    }
    else ids <- 1:length(data)
    
    ## perform CV
    params <- .doCrossValidation(init$pars, .loglik, data[ids], formula = formula,
                                 folds = kfolds, intervallSize = intervallSize,
                                 fixedpars = init$fpars,
                                 experimentDesign = colData(ggd), sf = sizeFactors(ggd),
                                 method = getDefaults(settings, "optimMethod"),
                                 control = getDefaults(settings, "optimControl"))
    names(params) <- c("lambda", "theta")
    lambda <- params[1]
    theta <- params[2]

    cv <- TRUE
    futile.logger::flog.info("Done")
  }
  if(!FIXTHETA) family <- mgcv::nb(theta = theta)
  cdata <- colData(ggd)
  sf <- sizeFactors(ggd)
  lambda_vec<- rep(lambda, nsplines)

  futile.logger::flog.info("Fitting model")

  
  lambdaFun <- function(data, colData, sf, formula, family, lambdas, 
                        chunkIndex) {
    suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
    id <- runValue(data$id)
    df <- .meltGTile(list(data), colData, sf, formula)
    mod <- mgcv::gam(formula, df[[1]], family = family, sp = lambdas)
    fits <- getFunctions(mod, df[[1]])
    rowindx <- match(c(start(chunkIndex[id]), end(chunkIndex[id])), data$pos)
    fits <- fits[seq(from = rowindx[1], to = rowindx[2]),]
    smooths <- extractSplines(mod)
    attr(smooths, "id") <- id
    vcov <- mgcv::vcov.gam(mod)
    upperTri_vcov <- upper.tri(vcov, diag = TRUE)
    vcov <- vcov[upperTri_vcov]
    attr(vcov, "smooths") <- na.omit(smooths)$smooth
    return(list(fits = fits, smooths = smooths, vcov = vcov))
  }

  res <- bplapply(data, lambdaFun, colData = cdata, sf = sf, formula = formula,
                  family = family, lambdas = lambda_vec, chunkIndex = chunkCoords)

  smooths <- list(chunkIndex = chunkCoords,
                  splines = list())
  smooths$splines <- lapply(res, function(y) y$smooths)
  names(smooths$splines) <- chunkCoords$id
  vcov <- lapply(res, function(y) y$vcov)
  names(vcov) <- chunkCoords$id
  fits <- lapply(res, function(y) y$fits)

  ##combine
  fits <- data.table::rbindlist(fits)

  ## create GenoGAM object
  fitparams <- c(lambda = unname(lambda), theta = unname(family$getTheta(trans = TRUE)),
                 CoV = sqrt(1/unname(family$getTheta(trans = TRUE))), penorder = m)
  cvparams <- c(kfolds = kfolds, ncv = ncv, size = intervallSize, cv = cv)
  tempFormula <- gsub("pos", "x", formula)
  saveFormula <- as.formula(paste(tempFormula[2], tempFormula[1], tempFormula[3]))
  genogamObject <- .GenoGAM(design = saveFormula, fits = fits,
                            positions = rowRanges(ggd), smooths = smooths, vcov = vcov,
                            experimentDesign = as.matrix(colData(ggd)[,na.omit(vars), drop = FALSE]),
                            fitparams = fitparams, family = family,
                            cvparams = cvparams, settings = settings,
                            tileSettings = tileSettings)

  futile.logger::flog.info("DONE")
  return(genogamObject)
}

#' Get all the splines and derivatives
#'
#' @param mod The mgcv model.
#' @param data The mgcv input data.
#' @param experimentMatrix A matrix object representing the experiment design.
#' @param m The penalization order of the P-splines.
#' @return A data.table of all the splines and derivatives.
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
getFunctions <- function(mod, data) {
    ## some constants
    samplenames <- unique(data$ID)
    vars <- .getVars(mod$formula)[-1]
    num_var <- length(mod$smooth)
    len <- length(unique(data$pos))

    ## get data
    pred <- mgcv::predict.gam(mod, type = "iterms", se.fit = TRUE)

    ## result data.table
    res <- data.table::data.table(matrix(0, len, num_var*2))
    nameslist <- c(colnames(pred$fit), paste("se", colnames(pred$se.fit), sep = "."))
    data.table::setnames(res, names(res), nameslist)    

    ## insert prediction values
    for(ii in 1:num_var) {
        if(ii == 1) intercept <- coefficients(mod)[1]
        else intercept <- 0
        label <- mod$smooth[[ii]]$label
        selabel <- paste("se", label, sep = ".")
        if(is.na(vars[ii])) {
            uniquepos <- !duplicated(data[, "pos"])
            res[[label]] <- pred$fit[uniquepos, label] + intercept
            res[[selabel]] <- pred$se.fit[uniquepos, label]
        }
        else {
            uniquepos <- !duplicated(data[as.logical(data[[vars[ii]]]), "pos"])
            res[[label]] <- pred$fit[as.logical(data[[vars[ii]]]),][uniquepos, label] + intercept
            res[[selabel]] <- pred$se.fit[as.logical(data[[vars[ii]]]),][uniquepos, label]
        }
    }
    newColNames <- gsub("pos", "x", names(res))
    data.table::setnames(res, names(res), newColNames)
    return(res)
}

#' The B-Spline function.
#'
#' A function to construct B-Spline bases.
#'
#' @param x The numeric vector of x values at which to evaluate the function
#' @param k A vector of knot positions
#' @param m The B-Spline basis order as order - 1, e.g. m = 2 is cubic
#' @param derivative The order of derivative
#' @return A matrix with dimensions length(x) * length(k)
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
bspline <- function(x, k, m = 2, derivative = 0) {
    res <- splines::spline.des(k, x, m + 2, rep(derivative,length(x)))$design
    return(res)
}    

#' Extract splines from model.
#'
#' A function to extract splines parameters from the gam model.
#'
#' @param mod A gam model object
#' @return A data.table of spline parameters
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
extractSplines <- function(mod){

    num_var <- length(mod$smooth)
    res <- NULL
    all_coef <- coefficients(mod)
    for (jj in 1:num_var) {
        smooth <- mod$smooth[[jj]]
        smoothType <- attr(smooth, "class")

        if(smoothType[1] == "pspline.smooth") {
            kn <- smooth$knots
            first <- smooth$first.para
            last <- smooth$last.para
            nConstraints <- attr(smooth, "nCons")
            coefs <- all_coef[first:last]
            name <- smooth$label
            name <- gsub("pos", "x", name)
        }

        if (nConstraints > 0) {
            qrc <- attr(smooth, "qrc")
            qrDim <- dim(qrc$qr)
            insertedZeros <- rep(0, nConstraints)
            y <- matrix(c(insertedZeros, coefs), qrDim)
            coefs <- qr.qy(qrc,y)
        }

        if(is.null(start)) start <- min(kn)
        if(is.null(end)) end <- max(kn)

        nas <- length(kn) - length(coefs)
        coefs <- c(coefs,rep(NA,nas))
        res <- rbind(res, data.frame(knots = kn, coefs = coefs, smooth = name))

    }
    attr(res, "intercept") <- all_coef["(Intercept)"]
    return(res)
}
