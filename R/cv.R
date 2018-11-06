## TODO: startup time of workers is approx 1/4 of the computation time. try to reduce
## e.g. package loading or passing setup object and initialize in the parallel

################################
## Cross Validation functions ##
################################

#' The cross validation function
#' 
#' @param ggd The GenoGAMDataSet object
#' @param setup The GenoGAMSetup object
#' @param coords The coordinates as a IRanges object
#' @param id The vector of ids
#' @param folds The number of folds
#' @param intervalSize The interval size of consecutive points to be left out
#' @param fn The log-likelihood function to be used
#' @param method The optim method
#' @param control Other control parameters according to optim
#' @param estimControl Control parameters according to parameter estimation method
#' @param ... Possibly other parameters
#' @return The optimized params (They are initially part of the setup object)
#' @noRd
.doCrossValidation <- function(ggd, setup, coords, id, folds, intervalSize,
                               fn, method = "Nelder-Mead",
                               control = list(maxit=100, fnscale=-1), ...) {
    
    setups <- vector("list", length(id))
    for (ii in 1:length(id)) {
        setups[[ii]] <- .initiate(ggd, setup, coords, id[ii])
    }
    names(setups) <- as.character(id)
    cvint <- .leaveOutConsecutiveIntervals(folds, intervalSize, 
                                           length(slot(setups[[1]],
                                                       "response")))
    ov <- getOverhangSize(ggd)
    par <- slot(setup, "params")
    initpars <- NULL
    if(is.null(par$lambda)) {
        initpars$lambda <- log(slot(setups[[1]], "params")[["lambda"]])
    }
    if(is.null(par$theta)) {
        initpars$theta <- log(slot(setups[[1]], "params")[["theta"]])
    }
    fixedpars <- list(lambda = par$lambda, theta = par$theta)
    estimControl <- slot(setup, "control")
    input <- paste0("Performing Cross-Validation with the following parameters:\n",
                    "  Tile indeces: ", paste(id, collapse = ","), "\n",
                    "  Number of folds: ", folds, "\n",
                    "  Interval size: ", intervalSize, "\n",
                    "  Optim method: ", method, "\n",
                    "  Maximal iterations: ", control$maxit, "\n",
                    "  Maximal L-BFGS iterations: ", estimControl$maxiter, "\n",
                    "  L-BFGS epsilon: ", estimControl$eps, "\n",
                    "  L-BFGS alpha: ", estimControl$alpha, "\n",
                    "  L-BFGS rho: ", estimControl$rho, "\n",
                    "  L-BFGS c: ", estimControl$c, "\n",
                    "  L-BFGS m: ", estimControl$m, "\n",
                    "  Parameters to optimize: ", paste(names(initpars), collapse = ","), "\n",
                    "  Value of fixed parameters: ", paste(fixedpars, collapse = ", "), "\n")
    futile.logger::flog.debug(input)

    if(flog.threshold() == "DEBUG" | flog.threshold() == "TRACE") {
        control$trace <- 10
    }

    ## reduce number of workers to reduce overhead if necessary
    ## Parameters
    currentBackend <- BiocParallel::registered()
    currentWorkers <- currentBackend[[1]]$workers
    currentFits <- folds * length(id)
    computationTime <- 0.2
    wakeUpTime <- 16

    ## Compute optimal number
    ## This is the solution to the equation argmin_m = (n * t)/m + u * m
    ## where
    ## n = currentFits = the total number of models to fit
    ## t = computationTime per Fit/Model
    ## m = number of workers (the parameter of interest)
    ## u = wakeUpTime = the time it takes a worker to start in SnowParam
    nworkers <- as.integer(sqrt((currentFits * computationTime)/wakeUpTime))
    backend <- BiocParallel::registered()[[1]]

    if(is(backend, "MulticoreParam") & currentWorkers > nworkers) {
        futile.logger::flog.debug(paste("Reducing number of workers during hyperparameter optimization to", nworkers))
        worker_backup <- currentBackend[[1]]$workers
        ## note, setting the workers in the variable, does change it in the
        ## overall setting because it is of class envRefClass
        currentBackend[[1]]$workers <- nworkers 
    }
    
    pars <- optim(initpars, fn, setup = setups, CV_intervals = cvint,
                  ov = ov, method = method, control = control, 
                  fixedpars = fixedpars, ...)
    params <- exp(pars$par)

    ## set the number of workers back to the specified number
    if(is(backend, "MulticoreParam") & currentWorkers > nworkers) {
        futile.logger::flog.debug(paste("Re-setting number of workers to", worker_backup))
        currentBackend[[1]]$workers <- worker_backup
    }

    futile.logger::flog.debug("Optimal parameter values:", params, capture = TRUE)
    
    if(length(params) == 1) {
        fixedpars[sapply(fixedpars, is.null)] <- params
        params <- unlist(fixedpars)
    }
    return(params)
}

#' The loglikelihood function
#' 
#' @param pars The parameters to be optimized. If .loglik is used within optim,
#' the pars argument will be coerced to vector even if list was given
#' @param setup The GenoGAMSetup object
#' @param CV_intervals A list of the indices of the data points split by
#' the folds
#' @param ov The size of the overlap. In order to be excluded from the
#' evaluation of the likelihood
#' @param fixedpars The parameters to be used in the model but kept
#' fixed, as they don't need optimization.
#' @param ... Other parameters
#' @return The mean log-likelihood over all models
#' @noRd
.loglik <- function(pars, setup, CV_intervals, ov, fixedpars, ...){

    if(is.null(fixedpars$lambda)) {
        fixedpars$lambda <- exp(pars[["lambda"]])
    }
    if(is.null(fixedpars$theta)) {
        fixedpars$theta <- exp(pars[["theta"]])
    }

    if(fixedpars$lambda/fixedpars$theta < 100) {
        fixedpars$lambda <- fixedpars$theta * 100
    }

    futile.logger::flog.debug(paste0("Using values: lambda = ", fixedpars$lambda, "and theta = ", fixedpars$theta))
    
    fullpred <- lapply(1:length(setup), function(y) {
        rep(NA, length(slot(setup[[1]], "response")))
    })
    names(fullpred) <- names(setup)
    ids <- expand.grid(folds = 1:length(CV_intervals), tiles = 1:length(setup))

    .local <- function(iter, ids, setup, CV_intervals) {
        suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
        id <- ids[iter,]
        
        trainsetup <- setup[[id$tiles]]
        trainX <- slot(trainsetup, "designMatrix")
        trainY <- slot(trainsetup, "response")
        testX <- trainX[CV_intervals[[id$folds]],]
        slot(trainsetup, "designMatrix") <- trainX[-CV_intervals[[id$folds]],]
        slot(trainsetup, "response") <- trainY[-CV_intervals[[id$folds]]]
        slot(trainsetup, "offset") <- slot(trainsetup, "offset")[-CV_intervals[[id$folds]]]
                    
        betas <- .estimateParams(trainsetup)
            
        pred <- exp(testX %*% betas$par)
        return(pred)
    }
    
    for(ii in 1:length(setup)) {
        slot(setup[[ii]], "params")$theta <- fixedpars$theta
        slot(setup[[ii]], "params")$lambda <- fixedpars$lambda
    }
   
    cvs <- BiocParallel::bplapply(1:nrow(ids), .local, ids = ids,
                                  setup = setup, CV_intervals = CV_intervals)

    for(ii in 1:length(cvs)) {
        id <- ids[ii,]
        fullpred[[id$tiles]][CV_intervals[[id$folds]]] <- cvs[[ii]]
    }

    res <- lapply(1:length(fullpred), function(y) {
        dens <- dnbinom(slot(setup[[y]], "response"), size = fixedpars$theta,
                        mu = fullpred[[y]], log = TRUE)
        return(dens)
    })

    #############################################################
    ## out-of-sample log-likelihood of the fit on the chunk part of a tile
    
    ## ll: sum CV log-lik by regions and parameter combination
    ## we take the average i.e we assume all regions have the same length.
    nfuns <- .nfun(setup[[1]])
    
    ll <- mean(sapply(res, function(y) {
        if(ov < 1) {
            res <- sum(y)
        }
        else {
            ymat <- matrix(y, ncol = nfuns)
            borders <- c(1:ov, (nrow(ymat) - ov + 1):nrow(ymat))
            if(length(borders) >= nrow(ymat)) {
                futile.logger::flog.error("The overhang size covers the entire tile. Change parameter to a lower meaningful value. See getOverhangSize().")
            }
            res <- sum(ymat[-borders,])
        }
        return(res)
    }))
    
    return(ll)
}


#' Create folds for Crossfold Validation bu consecutive intervals
#'
#' @param folds Number of folds.
#' @param intervalSize The size of the consecutive intervals.
#' @param tileSize The size of one tile.
#' @return A list of as many elements as there are folds. Each element
#' contains the rows that are part of this fold.
#' @noRd
.leaveOutConsecutiveIntervals <- function(folds, intervalSize, tileSize){
    if(intervalSize >= 50) {
        warning("CV using large interval may only make sense if all experiments have replicates")
    }
 
    ## number of intervals
    n = round(tileSize/intervalSize)
    ## break points
    bps = round(seq(1, tileSize + 1, length.out = n + 1))
    leftout_intervals  = split(sample(n),rep(1:folds, length = n))
    res = lapply(
        leftout_intervals, 
        function(lo){
            ## we return the indices in increasing order to preserve
            ## as much as possible the Genomic ranges for subsequent subsetting
            sort(unlist(lapply(lo, function(i) bps[i]:(bps[i+1]-1))))
        }
        )
    return(res)
}
