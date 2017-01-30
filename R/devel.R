library(GenoGAM)
load("/s/project/coreProm/peakCalling/sacCer2/rdata/tfiib_data.RData")
library(Matrix)
source("GenoGAMSetup-class.R")
source("sf.R")

ggs <- setupGenoGAM(ggd)
BiocParallel::register(BiocParallel::MulticoreParam(workers=4))

cvlog <- data.table::data.table(lambda = numeric(),
                                theta = numeric(),
                                ll = numeric())## to be maybe included in CV

#B spline basis
.bspline <- function(x, k, ord = 2, derivative = 0) {
  res <- splines::spline.des(k, x, ord + 2, rep(derivative,length(x)), sparse=TRUE)$design
  return(res)
}

#' build a block matrix from a template submatrix and a design matrix
#' @noRd
.blockMatrixFromDesignMatrix <- function(template, design) {
  ## create 4-dim array by 'inserting' the template into the desing matrix
  arr <- array(template, c(dim(template),dim(design)))
  dims <- dim(arr)
  multP <- c(3,4,1,2)
  reduceP <- c(3,1,4,2)
  ## permute array for correct multiplication
  multArr <- aperm(arr, multP)*as.vector(design)
  ## permute array for correct reduction
  reducedArr <- aperm(multArr, reduceP)
  ## reduce 4-dim array to 2-dim matrix
  dim(reducedArr) <- c(nrow(template)*nrow(design), ncol(template)*ncol(design))
  return(reducedArr)
}

genogam <- function(ggd, lambda = NULL, theta = NULL, family = "nb", H = 0,
                    bpknots = 20, kfolds = 10, intervalSize = 20, 
                    intervals = 20, order = 2, m = 2) {
    
    settings <- slot(ggd, "settings")
    tileSettings <- tileSettings(ggd)
    check <- checkSettings(ggd)
    coords <- getIndexCoordinates(ggd)

    if(!check) break

    ggs <- setupGenoGAM(ggd, lamda = lambda, theta = theta, family = family, 
                        H = H, bpknots = bpknots, order = order,
                        penorder = m)

    ## Cross Validation
    cv <- FALSE

    if(is.null(lambda) | is.null(theta)) {
        futile.logger::flog.info("Estimating parameters")
        
        ## get the tile ids for CV
        sumMatrix <- sum(ggd)
        ## ncv set to a number, such that CV does not compute more models than
        ## the actual genogam run. 
        ## 50 is the average expected number of iterations
        ncv <- min(intervals, ceiling(length(coords)/(kfolds*50)))
        if(ncv < length(coords)) {
            if(sum(sapply(colData(ggd), sum)) == nrow(colData(ggd))) {
                rsums <- rowSums(sumMatrix[[1]])
                ids <- order(rsums, decreasing = TRUE)[1:ncv]
            }
            else {
                pvals <- suppressMessages(suppressWarnings(.deseq(sumMatrix[[1]], colData(ggd))))
                ids <- order(pvals)[1:ncv]
            }
        }
        else ids <- 1:length(data)
    
        start <- proc.time()
        params <- .doCrossValidation(ggd, setup = ggs, coords = coords, 
                                     id = ids, folds = kfolds, 
                                     intervalSize = intervalSize,
                                     fn = .loglik, ov = getOverhangSize(ggd), 
                                     method = "Nelder-Mead", #GenoGAM:::getDefaults(settings, "optimMethod"),
                                     control = GenoGAM:::getDefaults(settings, "optimControl")) 
        proc.time() - start
        names(params) <- c("lambda", "theta")
        lambda <- params[1]
        theta <- params[2]
        slot(ggs, "params")$lambda <- lambda
        slot(ggs, "params")$theta <- theta
        
        cv <- TRUE
        futile.logger::flog.info("Done")
    }
    
    futile.logger::flog.info("Fitting model")

    lambdaFun <- function(id, data, init, coords) {
        ## suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))

        setup <- .initiate(data, init, coords, id)
        betas <- .estimateParams(setup)
        slot(setup, "beta") <- betas$par
        slot(setup, "fits") <- .getFits(setup)
        slot(setup, "designMatrix") <- new("dgCMatrix")
        slot(setup, "penaltyMatrix") <- new("dgCMatrix")
        slot(setup, "response") <- numeric()

        return(setup)
    }

    ids <- coords$id
    res <- BiocParallel::bplapply(ids, lambdaFun, 
                                  data = ggd, init = ggs, coords = coords)
}


.getFits <- function(setup) {
    nSplines <- length(.getVars(slot(setup, "formula"), "covar"))
    dims <- dim(slot(setup, "designMatrix"))
    rows <- split(1:dims[1], cut(1:dims[1], nSplines, labels = 1:nSplines))
    cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    fits <- as.vector(sapply(1:nSplines, function(y) {
        X <- slot(setup, "designMatrix")[rows[[y]], cols[[y]]]
        beta <- slot(setup, "beta")[cols[[y]], 1]
        res <- as.vector(X %*% beta)
        return(res)
    }))
    return(fits)
}


.buildDesignMatrix <- function(ggd, setup, pos) {
    ## get knots
    knots <- .getKnots(pos, setup)
    order <- slot(setup, "params")$order

    ## build matrix
    x <- as(.bspline(pos, knots, order),"dgCMatrix")
    design <- as.matrix(colData(ggd))
    control <- rep(1, nrow(design))
    design <- cbind(control, design)
    X <- as(.blockMatrixFromDesignMatrix(x, design), "dgCMatrix")
    return(X)
}

.getKnots <- function(pos, setup) {
    ## knots <- slot(setup, "knots")[[chrom]]
    knots <- slot(setup, "knots")[[1]]
    inner <- which(knots < max(pos) & knots > min(pos))
    knots <- knots[(inner[1] - 4):(inner[length(inner)] + 4)]
    return(knots)
}

.buildResponseVector <- function(ggd, gr) {
    y <- assay(ggd)[start(gr):end(gr),]
    
    Y <- unname(unlist(as.data.frame(y)))
    return(Y)
}

.estimateParams <- function(ggs) {
    ## turn the beta matrix into a 1-column matrix
    betas <- slot(ggs, "beta")

    distr <- slot(ggs, "family")
    X <- slot(ggs, "designMatrix")
    y <- slot(ggs, "response")
    offset <- slot(ggs, "offset")
    params <- slot(ggs, "params")
    S <- slot(ggs, "penaltyMatrix")

    if (distr == "nb") {    
        likelihood <- .likelihood_penalized_nb
        gradient <- .gradient_likelihood_penalized_nb
    }

    res <- optim(betas, likelihood, gradient, X = X, y = y, offset = offset,
                 theta = params$theta, lambda = params$lambda, S = S, 
                 method = "L-BFGS", control = list(fnscale=-1, maxit = 1000))
    return(res)
}

.initiate <- function(ggd, setup, index, id) {
    tile <- IRanges::ranges(index[index$id == id,])
    pos <- pos(SummarizedExperiment::rowRanges(ggd)[tile])
    ## chrom <- GenomeInfoDb::seqlevelsInUse(tile)

    ## initiate matrices and vectors
    slot(setup, "designMatrix") <- .buildDesignMatrix(ggd, setup, pos)
    slot(setup, "response") <- .buildResponseVector(ggd, tile)
    numBetas <- dim(slot(setup, "designMatrix"))[2]
    slot(setup, "beta") <- matrix(mean(slot(setup, "response"), na.rm = TRUE), numBetas, 1)
    slot(setup, "offset") <- unname(rep(sizeFactors(ggd)[colnames(ggd)], each = width(tile)))
    slot(setup, "knots")[[1]] <- .getKnots(pos, setup)
        
    ## iniate parameters
    params <- slot(setup, "params")
    if(is.null(params$lambda)) {
        params$lambda <- numBetas
    }
    if(is.null(params$theta)) {
        params$theta <- 1
    }
    slot(setup, "params") <- params

    ## penaltyMatrix
    S <- buildSMatrix(numBetas, params$penorder)

    ## add identity matrix to penalization
    if(params$H > 0) {
        I <- buildIMatrix(numBetas, params$H)
        S <- S + I
    }
    slot(setup, "penaltyMatrix") <- S

    return(setup)
}

#penalized likelihood to be maximized \ negbin
.likelihood_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  n <- dim(X)[1]
  eta <- offset + X%*%beta
  mu <- exp(eta)
  aux1 <- theta + y
  aux2 <- theta + mu
  ## pull the log inside to use gamma and factorial in log space due to 
  ## possibly very high numbers
  l <- sum(lgamma(aux1) - (lfactorial(y) + lgamma(theta))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
  pen <- t(beta) %*% S %*% beta
  return(l[1]-lambda*pen[1,1])
}  

#gradient of penalized likelihood \negbin
.gradient_likelihood_penalized_nb <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}

.doCrossValidation <- function(ggd, setup, coords, id, folds, intervalSize,
                               fn, method = "Nelder-Mead",
                               control = list(maxit=100, fnscale=-1), ...) {
    
    setups <- vector("list", length(id))
    for (ii in 1:length(id)) {
        setups[[ii]] <- .initiate(ggd, setup, coords, id[ii])
    }
    names(setups) <- as.character(id)
    cvint <- .leaveOutConsecutiveIntervals(folds, intervalSize, 
                                           length(slot(setups[[1]], "response")))

    par <- slot(setup, "params")
    initpars <- NULL
    if(is.null(par$lamda)) {
        initpars <- c(initpars, lambda = log(slot(setups[[1]], "params")[["lambda"]]))
    }
    if(is.null(par$theta)) {
        initpars <- c(initpars, theta = log(slot(setups[[1]], "params")[["theta"]]))
    }
    fixedpars <- c(par["lambda"], par["theta"])
    pars <- optim(initpars, fn, setup = setups, CV_intervals = cvint,
                  method = method, control = control, 
                  fixedpars = fixedpars, ...)
    params <- exp(pars$par)
    
    if(length(params) == 1) {
        fixedpars[sapply(fixedpars, is.null)] <- params
        params <- unlist(fixedpars)
    }
    return(params)
}

.loglik <- function(pars, setup, CV_intervals, ov, fixedpars, ...){

    if(is.null(fixedpars$lambda)) {
        fixedpars$lambda <- exp(pars[["lambda"]])
    }
    if(is.null(fixedpars$theta)) {
        fixedpars$theta <- exp(pars[["theta"]])
    }

    if(fixedpars$theta < 1e-3) {
        fixedpars$theta <- 1e-3
    }

    fullpred <- lapply(1:length(setup), function(y) {
        rep(NA, length(slot(setup[[1]], "response")))
    })
    names(fullpred) <- names(setup)
    ids <- expand.grid(folds = 1:length(CV_intervals), tiles = 1:length(setup))

    lambdaFun <- function(iter, ids, setup, CV_intervals) {
##        suppressPackageStartupMessages(require(GenoGAM, quietly = TRUE))
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
    
    cvs <- BiocParallel::bplapply(1:nrow(ids), lambdaFun, ids = ids, setup = setup,
                    CV_intervals = CV_intervals)

    for(ii in 1:length(cvs)) {
        id <- ids[ii,]
        fullpred[[id$tiles]][CV_intervals[[id$folds]]] <- cvs[[ii]]
    }

    res <- lapply(1:length(fullpred), function(y) {
        dens <- dnbinom(slot(setup[[y]], "response"), size = fixedpars$theta, mu = fullpred[[y]], log = TRUE)
        return(dens)
    })

    #############################################################
    ## out-of-sample log-likelihood of the fit on the chunk part of a tile
    
    ## ll: sum CV log-lik by regions and parameter combination
    ## we take the average i.e we assume all regions have the same length.
    ll <- mean(sapply(res, function(y) {
        borders <- c(1:ov, (length(y) - ov + 1):length(ov))
        sum(y[-borders])
    }))
    templog <- data.table::data.table(lambda = fixedpars$lambda,
                                      theta = fixedpars$theta,
                                      ll = ll)
    cvlog <<- rbind(cvlog, templog)
    return(ll)
}



############################ Tests

k <- placeKnots(1:5200, 260)
x1 <- rep(1:2660,2)
x2 <- rep(2541:5200,2)
x <- rep(1:5200,2)
y1 <- as.vector(assay(ggd)$input[1:5200])
y2 <- as.vector(assay(ggd)$IP[1:5200])
y <- c(y1, y2)
z1 <- rep(c(0,1), each = 2660)
z <- rep(c(0,1), each = 5200)
full <- mgcv::gam(y ~ s(x, bs = "ps", k = 260) + s(x, bs = "ps", k = 260, by = z), sp = rep(36.82, 2), knots = list(k, k), family = mgcv::nb(theta = 3.1))
left <- mgcv::gam(c(y1[1:2660], y2[1:2660]) ~ s(x1, bs = "ps", k = 135) + s(x1, bs = "ps", k = 135, by = z1), sp = rep(36.82*139/264, 2), knots = list(k[1:139], k[1:139]), family = mgcv::nb(theta = 3.1))
right <- mgcv::gam(c(y1[2541:5200], y2[2541:5200]) ~ s(x2, bs = "ps", k = 135) + s(x2, bs = "ps", k = 135, by = z1), sp = rep(36.82*139/264, 2), knots = list(k[126:264], k[126:264]), family = mgcv::nb(theta = 3.1))

pred.full <- predict(full, type = "iterms")
pred.left <- predict(left, type = "iterms")
pred.right <- predict(right, type = "iterms")

par(mfrow = c(2,1))
plot(pred.full[1:5200,1], type = 'l')
lines(x1[1:2660], pred.left[1:2660, 1], col = "red")
lines(x2[1:2660], pred.right[1:2660, 1], col = "green")

plot(pred.full[5201:10400,2], type = 'l')
lines(x1[1:2660], pred.left[2661:5320, 2], col = "red")
lines(x2[1:2660], pred.right[2661:5320, 2], col = "green")

full <- mgcv::gam(y ~ s(x, bs = "ps", k = 260) + s(x, bs = "ps", k = 260, by = z), knots = list(k, k), family = nb())
left <- mgcv::gam(c(y1[1:2660], y2[1:2660]) ~ s(x1, bs = "ps", k = 135) + s(x1, bs = "ps", k = 135, by = z1), knots = list(k[1:139], k[1:139]), family = nb())
right <- mgcv::gam(c(y1[2541:5200], y2[2541:5200]) ~ s(x2, bs = "ps", k = 135) + s(x2, bs = "ps", k = 135, by = z1), knots = list(k[126:264], k[126:264]), family = nb())



## lfbgs_sparse_and_compute_Hinv_opt_generic.R



## @author : Mathilde GALINIER   #


library(Matrix)

setClass(Class="result",
         representation(
           model_matrix = "dgCMatrix",
           beta = "matrix",
           Hinv.app = "dgCMatrix"
         )
)

setClass(Class="result2",
         representation(
           model_matrix = "dgCMatrix",
           beta = "matrix"
         )
)

## #create a banded matrix with 5 diagonals
## diag_5 <- function(upper1, upper2, main){
##   out <- Matrix(0,length(main),length(main),sparse=TRUE)
##   diag(out) <- main
##   indx <- seq.int(length(upper1))
##   out[cbind(indx+1,indx)] <- upper1
##   out[cbind(indx,indx+1)] <- upper1
##   indx <- seq.int(length(upper2))
##   out[cbind(indx+2,indx)] <- upper2
##   out[cbind(indx,indx+2)] <- upper2
##   return(out)
## }

## #Second differentiate of beta for the penalization term
## scnd_diff_beta <- function(n,nb_var=1){
##   res <- diag_5(rep(-4,n-1),rep(1,n-2),rep(6,n))
##   res[1,1]<-1
##   res[n,n]<-1
##   res[1,2]<--2
##   res[2,1]<--2
##   res[n,(n-1)]<--2
##   res[(n-1),n]<--2
##   res[2,2]<-5
##   res[(n-1),(n-1)]<-5
##   if (nb_var>1){
##     for (var in 1:(nb_var-1)){
##       coord <- dim(res)[2]/nb_var*var
##       res[(coord-2):(coord),(coord+1):(coord+2)] <- 0
##       res[(coord+1):(coord+2),(coord-2):(coord)] <- 0
##       res[coord-1,coord] <- -2 
##       res[coord,coord-1] <- -2
##       res[coord+2,coord+1] <- -2
##       res[coord+1,coord+2] <- -2
##       res[coord,coord] <- 1
##       res[coord+1,coord+1] <- 1
##       res[coord-1,coord-1] <- 5
##       res[coord+2,coord+2] <- 5}
##   }
##   return(res)
## }


#penalized likelihood to be maximized \ gaussian
likelihood_penalized_gaussian_log <- function(beta,X,y,lambda,S){
  n <- dim(X)[1]
  sigma2 <- var(y)
  l <- -n/2*log(2*pi*sigma2)-1/(2*sigma2)*sum((y-exp(X%*%beta))^2)
  pen <- t(beta) %*% S %*% beta
  return(l-lambda*pen[1,1])
}  

#penalized gradient of likelihood \ gaussian
gradient_likelihood_penalized_gaussian_log <- function(beta,X,y,lambda,S){
  sigma2 <- var(y)
  mu <- exp(X%*%beta)
  z <- 1/(sigma2[1])*mu*(y-mu)  
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}

#Compute the penalized Hessian in the gaussian case
compute_hessian_gaussian <- function(beta,X,y,lambda,S){
  sigma2 <- var(y)
  mu <- exp(X%*%beta)
  d <- mu*(y-2)/sigma2
  # S[t(X)%*%X==0]<-0 ####
  res <- t(X) %*% bandSparse(dim(X)[1], k = 0, diag = c(list(d[,1]))) %*% X  - 2*lambda*S
  return (res)
}

#penalized likelihood to be maximized \ negbin
likelihood_penalized_negbin_log <- function(beta,X,y,offset,theta,lambda,S){
  n <- dim(X)[1]
  eta <- offset + X%*%beta
  mu <- exp(eta)
  aux1 <- theta + y
  aux2 <- theta + mu
  l <- sum(log(gamma(aux1)/(factorial(y)*gamma(theta)))) + t(y) %*% eta + n*theta*log(theta) - t(aux1) %*% log(aux2)
  pen <- t(beta) %*% S %*% beta
  return(l[1]-lambda*pen[1,1])
}  

#gradient of penalized likelihood \negbin
gradient_likelihood_penalized_negbin_log <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  z <- (y-mu)/(1+mu/theta)
  res <- t(X)%*%z
  pen <- S %*% beta
  return (res[,1]-2*lambda*pen[,1])
}

#Compute the penalized Hessian in the negbin case
compute_hessian_negbin <- function(beta,X,y,offset,theta,lambda,S){
  eta <- offset + X%*%beta
  mu <- exp(eta)
  d <- - mu * (1 + y/theta) / ((1 + mu/theta)^2) 
  res <- t(X) %*% bandSparse(dim(X)[1], k = 0, diag = c(list(d[,1])))%*% X - 2*lambda*S 
  return (res)
}

#main function 
#struct_design_matrix is a matrix representing the dependencies between y and betas
lbfgs_sparse_and_compute_Hinv_opt_generic<-function(x_data, nb_var=1, nb_beta=1, nb_knots, ord=2, lambda, theta=1,distribution, 
                                                    list_y, list_offsets, struct_design_matrix, tol, hinv = "TRUE", method="solve", nb_spl = 10) 
{ 
  n <- length(x_data)
  
  # Generate x_knots 
  m <- ord +1 
  nk <- nb_knots - ord
  xu <- max(x_data)
  xl <- min(x_data)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl)/(nk - 1)
  x_knots <- seq(xl - dx * (m), xu + dx * (m), length = nk + 2 * m)
  nb_knots_upd <- length(x_knots)
  
  # Initialisation of variables
  if(nb_var > 1) {
    y <- NULL
    offset <- NULL
    for(var in 1:nb_var){
      y <- c(y, list_y[[var]])
      offset <- c(offset,list_offsets[[var]])
    }}
  beta <-matrix(log(mean(y)),nb_beta*nb_knots,1)
  
  # Create design matrices
  x <- as(bspline(x_data, x_knots, ord),"dgTMatrix")
  if (nb_var==1){
    X <- x
  }else{
    X <- Matrix(0,nb_var*n,nb_beta*nb_knots)
    X <- as(X,"dgTMatrix")
    for(var in 1:nb_var){
      for(b in 1:nb_beta){
        if(struct_design_matrix[var,b]==1){
          X@i <- c(X@i, x@i+as.integer((var-1)*n))
          X@j <- c(X@j, x@j+as.integer((b-1)*nb_knots))
          X@x <- c(X@x,x@x)
        }
      }
    }
  }
  X <- as(X,"CsparseMatrix")
  
  # Compute the second differences matrix of beta (for penalized term)
  S<- scnd_diff_beta(dim(X)[2],nb_beta)
  
  # Compute beta 
  if (distribution=="gaussian"){
    res.opt <- optim(beta,likelihood_penalized_gaussian_log,gradient_likelihood_penalized_gaussian_log,
                     X=X,y=y,lambda = lambda, S=S,method="L-BFGS",control=list(fnscale=-1))
  }else if (distribution=="negbin"){
    res.opt <- optim(beta,likelihood_penalized_negbin_log,gradient_likelihood_penalized_negbin_log,
                     X=X,y=y,offset = offset, theta = theta,lambda = lambda, S=S, method="L-BFGS",control=list(fnscale=-1))
    
  }
  beta <- res.opt$par
  
  if (hinv == "TRUE"){
    
    #   Compute the Hessian matrix
    if (distribution=="gaussian"){
      H <- compute_hessian_gaussian(beta,X,y,lambda,S)
    }else if (distribution=="negbin"){
      H <- compute_hessian_negbin(beta,X,y,offset=offset,theta,lambda,S)
    }
    
    if(method=="solve"){
      # ############ FIRST ATTEMPT ###################
      # Computation of the inverse of H
      Hinv <- solve(H)
      Hinv.app <- Hinv
      Hinv.app[abs(Hinv.app) < tol] <- 0
      # ##############################################  
    }else if(method=="cutoff") {
      ##############SECOND ATTEMPT ###################
  
      # C = chol(-H, pivot=TRUE)
      C <- Cholesky(-H)
      
      system.time(
        i_x <- lapply(
          1:ncol(H),
          function(j){
            ej = rep(0, nrow(H))
            ej[j] = 1
            s = solve(C,ej)
            keep = which(abs(s)>tol)
            list(i=keep, x=s[keep])
          }
        )
      )
      
      Hinv.app = sparseMatrix(
        i = unlist(lapply(i_x, function(l) l$i)),
        j = rep(1:nrow(H), times=sapply(i_x, function(l) length(l$i))),
        x = unlist(lapply(i_x, function(l) -l$x))
      )
      
      ##################################################
    }else if(method=="footprint") {
      ##############THIRD ATTEMPT ###################
      
      #### The code below is buggy 
      
    # #Create a sparsity pattern for Hinv
    #   ll <- length(which(X[,ceiling(ncol(H)/2)]!=0))
    #   aux <- seq(1:((nb_spl-1)*ll))
    #   keep <- lapply(1:ncol(H), function(i){
    #     ind <- which(X[,i]!=0)
    #     list_ind <- c(ind, aux+ind[length(ind)])
    #     ifelse(list_ind < nrow(X),list_ind,nrow(X))
    #     })
    #   
    #   big_X = sparseMatrix(
    #     i = unlist(keep),
    #     j = rep(1:ncol(H), times=sapply(keep, length)),
    #     x = 1
    #   )
    #   
    #   tX_X <- t(big_X)%*%big_X
    #   
    #   ind <- which(X[,50]!=0)
    #   keep <- lapply(1:ncol(H), function(i) which(tX_X[,i]!=0))
    # 
    #   Hinv.app = sparseMatrix(
    #     i = unlist(keep),
    #     j = rep(1:ncol(H), times=sapply(keep, length)),
    #     x = 1
    #   )
    #   
    # #Compute the approximation of the inverse of H
    #   x = Hinv.app@x
    #   for(j in 1:ncol(H)){
    #       ej = rep(0, nrow(H))
    #       ej[j] = 1
    #       s = solve(H,ej)
    #       x[1 + Hinv.app@p[j]:(Hinv.app@p[j+1]-1)] = s[keep[[j]]]
    #   }
    #   Hinv.app@x <- x
      
         ##################################################
    }
    return(new("result",model_matrix=X, beta = beta, Hinv.app = Hinv.app))
  }else{
    return(new("result2",model_matrix=X, beta = beta))
  }
}

######################################
### test_2_variables.R
#####################################

library(mgcv)
library(data.table)

n <- 1000
x_data <- seq(0, 1, length.out = n)
mu <- sin(seq(-6 , 6, length.out = n))+2
y1 <- rnbinom(length(x_data), mu=exp(mu), size=1)
y2 <- rnbinom(length(x_data), mu=exp(mu), size=1)

list_y <- list(y1,y2)
nb_var <- 2
nb_beta <- 2

nb_knots <- 50
ord<-2
tol <- 0.0001

lambda <- 0
theta <- 1

offset0 <- rep(0, length(y1))
list_offsets <- list(offset0,offset0)

struct_design_matrix <- matrix(c(1,1,0,1),2,2)

distribution <- "negbin"

res <- lbfgs_sparse_and_compute_Hinv_opt_generic(x_data, nb_var, nb_beta, nb_knots, ord=ord, lambda, theta, distribution, list_y, list_offsets, struct_design_matrix,tol,
                                                                                                      hinv = "TRUE", method="solve") 

##########################################################3333

##MGCV
dt <- data.table(y = c(y1,y2), x_data = c(x_data,x_data), X=rep(1,2*n), IP=c(rep(0,n),rep(1,n)))
mod <- gam(y ~ s(x_data , k = nb_knots, sp = 0, bs = "ps",by=X)+ s(x_data , k = nb_knots, sp = 0, bs = "ps", by=IP), family=nb(theta=1),scale=1,data=dt)
pred <- predict(mod, type = "iterms", se.fit = TRUE)
str(pred)
pred$fit

########################PLOT
plot.gam(mod,select=1)

# y1 : mcgv results :
X_m <- res@model_matrix
vcov <- vcov.gam(mod)
lines(x_data, pred$fit[1:n,1], type = "l", col = "blue")
lines(x_data, pred$fit[1:n,1] - 1.96*pred$se.fit[1:n,1], col = "blue", lty = "dashed")
lines(x_data, pred$fit[1:n,1] + 1.96*pred$se.fit[1:n,1], col = "blue", lty = "dashed")

# plot(se1[1:n], pred$se.fit[1:n,1])
# #y1 : my results :
Sigma <- -res@Hinv.app
beta <- res@beta
# #plot(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] ,col="white")

se1 <- diag(sqrt(abs(X_m[1:n,1:nb_knots]%*%Sigma[1:nb_knots,1:nb_knots]%*%t(X_m[1:n,1:nb_knots]))))
lines(x_data, pred$fit[1:n,1] - 1.96*se1[1:n], col = "red", lty = "dashed")
lines(x_data, pred$fit[1:n,1] + 1.96*se1[1:n], col = "red", lty = "dashed")
se2 <- diag(sqrt(abs(X_m%*%Sigma%*%t(X_m))))

# 
# plot(se2[1:n], pred$se.fit[1:n,1])
# abline(0,1,col='red')
# 
# plot(se2[(n+1):(2*n)], pred$se.fit[(n+1):(2*n),2])
# abline(0,1,col='red')

lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots]-attr(pred,"constant")  , type = "l", col = "green")
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] - 1.96*se2[1:n]-attr(pred,"constant"), col = "green", lty = "dashed")
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] + 1.96*se2[1:n] -attr(pred,"constant"), col = "green", lty = "dashed")

#############################
plot.gam(mod,select=2)

# y2 : mcgv results :
lines(x_data, pred$fit[(n+1):(2*n),2], type = "l", col = "blue")
lines(x_data, pred$fit[(n+1):(2*n),2] - 1.96*pred$se.fit[(n+1):(2*n),2], col = "blue", lty = "dashed")
lines(x_data, pred$fit[(n+1):(2*n),2] + 1.96*pred$se.fit[(n+1):(2*n),2], col = "blue", lty = "dashed")


#y2 : my results
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)], type = "l", col = "green")
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] -1.96*se2[(1+n):(2*n)], col = "green", lty = "dashed")
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] + 1.96*se2[(1+n):(2*n)], col = "green", lty = "dashed")

se3 <- diag(sqrt(abs(X_m[1:n,1:nb_knots]%*%Sigma[(nb_knots+1):(2*nb_knots),(nb_knots+1):(2*nb_knots)]%*%t(X_m[1:n,1:nb_knots]))))
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] -1.96*se3[1:n], col = "orange", lty = "dashed")
lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] + 1.96*se3[1:n], col = "orange", lty = "dashed")

###########################################################################################################################
# #FOR REPLICATES
###########################################################################################################################

# n <- 100
# x_data <- seq(0, 1, length.out = n)
# mu1 <- sin(seq(-6 , 6, length.out = n)) +2
# mu2 <- sin(seq(-5 , 5, length.out = n)) +2
# y1 <- rnbinom(length(x_data), mu=exp(mu1), size=1)
# y1 <- c(y1,y1)
# y2 <- rnbinom(length(x_data), mu=exp(mu2), size=1)
# y2 <- c(y2,y2)
# 
# list_y <- list(y1,y2)
# nb_var <- 2
# nb_beta <- 2
# 
# nb_knots <- 50
# ord<-2
# tol <- 0.0001
# 
# lambda <- 0
# 
# offset0 <- rep(0, length(y1))
# list_offsets <- list(offset0,offset0)
# 
# # struct_design_matrix <- matrix(c(1,1,0,1,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0),5,4)
# struct_design_matrix <- matrix(c(1,1,0,1),2,2)
# 
# distribution <- "negbin"
# 
# x<-c(x_data,x_data)
# res <- lbfgs_sparse_and_compute_Hinv_opt_generic(x, nb_var, nb_beta, nb_knots, ord=ord, lambda, distribution, list_y, list_offsets, struct_design_matrix,tol)
# 
# ##MGCV
# dt <- data.table(y = c(y1,y2), x_data = x, X=rep(1,4*n), IP=c(rep(0,2*n),rep(1,2*n)))
# # mod <- gam(y ~ s(X - 1 , k = nb_knots, sp = 0, bs = "ps") + s(X - 1 , k = nb_knots, sp = 0, bs = "ps",by="IP"), family=nb(theta=1),scale=1,data=dt)  #######DOES NOT WORK
# mod <- gam(y ~ s(x_data -1 , k = nb_knots, sp = 0, bs = "ps",by=X)+ s(x_data - 1 , k = nb_knots, sp = 0, bs = "ps", by=IP), family=nb(theta=1),scale=1,data=dt)
# 
# # plot(mod, pages = 1)
# pred <- predict(mod, type = "iterms", se.fit = TRUE)
# str(pred)
# pred$fit
# 
# ########################PLOT
# plot.gam(mod,select=1)
# # y1 : mcgv results :
# X_m <- res[[1]]
# vcov <- vcov.gam(mod)
# se1 <- diag(sqrt(abs(X_m%*%vcov%*%t(X_m))))
# lines(x_data, pred$fit[1:n,1], type = "l", col = "blue")
# lines(x_data, pred$fit[1:n,1] - 1.96*pred$se.fit[1:n,1], col = "blue", lty = "dashed")
# lines(x_data, pred$fit[1:n,1] + 1.96*pred$se.fit[1:n,1], col = "blue", lty = "dashed")
# 
# #y1 : my results :
# Sigma <- -res[[3]]
# beta <- res[[2]]
# se2 <- diag(sqrt(abs(X_m%*%Sigma%*%t(X_m))))
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] -attr(pred,"constant") , type = "l", col = "green")
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] - 1.96*se2[1:n] -attr(pred,"constant"), col = "green", lty = "dashed")
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[1:nb_knots] + 1.96*se2[1:n]-attr(pred,"constant"), col = "green", lty = "dashed")
# 
# #############################
# plot.gam(mod,select=2)
# 
# # y2 : mcgv results :
# lines(x_data, pred$fit[(n+1):(2*n),2], type = "l", col = "blue")
# lines(x_data, pred$fit[(n+1):(2*n),2] - 1.96*pred$se.fit[(n+1):(2*n),2], col = "blue", lty = "dashed")
# lines(x_data, pred$fit[(n+1):(2*n),2] + 1.96*pred$se.fit[(n+1):(2*n),2], col = "blue", lty = "dashed")
# 
# 
# #y2 : my results
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)], type = "l", col = "green")
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] -1.96*se2[(n+1):(2*n)], col = "green", lty = "dashed")
# lines(x_data, X_m[1:n,1:nb_knots]%*%beta[(nb_knots+1):(2*nb_knots)] + 1.96*se2[(n+1):(2*n)], col = "green", lty = "dashed")




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

##############################################
### test_generic_function.R
##############################################
#============================#
##### test for 1 track #######
#============================#

###############################
###variables initialization ###
###############################

n <- 1000
x_data <- seq(0, 1, length.out = n)
mu <- sin(seq(-6 , 6, length.out = n))
y <- rnbinom(length(x_data), mu=exp(mu), size=1)

nb_knots<- 500 
ord<-2
tol <- 1e-4

lambda <- 50
theta <- 1

nb_var <- 1 #number of tracks
nb_beta <- 1 #number of smooth functions
nb_spl <- 10 #how many splines a priori correlate with a given spline (used for hessian inverse with pre-defined sparse structure)

struct_design_matrix <- 1

distribution <- "negbin"
offset <- rep(0,length(y))

#############################
####implemented function##### 
#############################

n <- length(x_data)

# Generate x_knots 
m <- ord +1 
nk <- nb_knots - ord
xu <- max(x_data)
xl <- min(x_data)
xr <- xu - xl
xl <- xl - xr * 0.001
xu <- xu + xr * 0.001
dx <- (xu - xl)/(nk - 1)
x_knots <- seq(xl - dx * (m), xu + dx * (m), length = nk + 2 * m)
nb_knots_upd <- length(x_knots)

# Initialisation of variables
if(nb_var > 1) {
  y <- NULL
  offset <- NULL
  for(var in 1:nb_var){
    y <- c(y, list_y[[var]])
    offset <- c(offset,list_offsets[[var]])
  }}
beta <-matrix(log(mean(y)),nb_beta*nb_knots,1)

# Create design matrices
x <- as(bspline(x_data, x_knots, ord),"dgTMatrix")
if (nb_var==1){
  X <- x
}else{
  X <- Matrix(0,nb_var*n,nb_beta*nb_knots)
  X <- as(X,"dgTMatrix")
  for(var in 1:nb_var){
    for(b in 1:nb_beta){
      if(struct_design_matrix[var,b]==1){
        X@i <- c(X@i, x@i+as.integer((var-1)*n))
        X@j <- c(X@j, x@j+as.integer((b-1)*nb_knots))
        X@x <- c(X@x,x@x)
      }
    }
  }
}
X <- as(X,"CsparseMatrix")
image(X)

# Compute the second differences matrix of beta (for penalized term)
S<- scnd_diff_beta(dim(X)[2],nb_beta)
image(S)

# Compute beta 
system.time(
if (distribution=="gaussian"){
  res.opt <- optim(beta,likelihood_penalized_gaussian_log,gradient_likelihood_penalized_gaussian_log,
                   X=X,y=y,lambda = lambda, S=S,method="L-BFGS",control=list(fnscale=-1))
}else if (distribution=="negbin"){
  res.opt <- optim(beta,likelihood_penalized_negbin_log,gradient_likelihood_penalized_negbin_log,
                   X=X,y=y,offset = offset, theta = theta,lambda = lambda, S=S, method="L-BFGS",control=list(fnscale=-1))
  
})
beta <- res.opt$par

plot(x_data, log(y),ylim=c(-2, max(log(y))))
lines(x_data, mu, col='green', lwd = 2)
lines(x_data, X%*%beta, col='red',lty=4, lwd = 2)

system.time(
#   Compute the Hessian matrix
if (distribution=="gaussian"){
  H <- compute_hessian_gaussian(beta,X,y,lambda,S)
}else if (distribution=="negbin"){
  H <- compute_hessian_negbin(beta,X,y,offset=offset,theta,lambda,S)
})

image(H)


# ############ FIRST ATTEMPT : SOLVE ###################
# Computation of the inverse of H
system.time(
Hinv1 <- solve(H))
Hinv1.app <- Hinv1
Hinv1.app[abs(Hinv1.app) < tol] <- 0
image(Hinv1.app)
norm(Hinv1-Hinv1.app)/norm(Hinv1)
# ######################################################  

##############SECOND ATTEMPT : CUTOFF###################

# C = chol(-H, pivot=TRUE)
C <- Cholesky(-H)

i_x <- lapply(1:ncol(H),function(j){
  ej = rep(0, nrow(H))
  ej[j] = 1
  s = solve(C,ej)
  keep = which(abs(s)>tol)
  list(i=keep, x=s[keep])
}
)


Hinv2.app = sparseMatrix(
  i = unlist(lapply(i_x, function(l) l$i)),
  j = rep(1:nrow(H), times=sapply(i_x, function(l) length(l$i))),
  x = unlist(lapply(i_x, function(l) -l$x))
)

image(Hinv2.app)
norm(Hinv1-Hinv2.app)/norm(Hinv1)
#####################################################

##############THIRD ATTEMPT : FOOTPRINT #################

#### The code below is buggy 

#Create a sparsity pattern for Hinv

# ll <- length(which(X[,ceiling(ncol(H)/2)]!=0))
# aux <- seq(1:((nb_spl-1)*ll))
# keep <- lapply(1:ncol(H), function(i){
#   ind <- which(X[,i]!=0)
#   list_ind <- c(ind, aux+ind[length(ind)])
#   ifelse(list_ind < nrow(X),list_ind,nrow(X))
# })
# 
# big_X = sparseMatrix(
#   i = unlist(keep),
#   j = rep(1:ncol(H), times=sapply(keep, length)),
#   x = 1
# )
# 
# tX_X <- t(big_X)%*%big_X
# 
# ind <- which(X[,50]!=0)
# keep <- lapply(1:ncol(H), function(i) which(tX_X[,i]!=0))
# 
# Hinv3.app = sparseMatrix(
#   i = unlist(keep),
#   j = rep(1:ncol(H), times=sapply(keep, length)),
#   x = 1
# )
# 
# #Compute the approximation of the inverse of H
# x = Hinv3.app@x
# for(j in 1:ncol(H)){
#   ej = rep(0, nrow(H))
#   ej[j] = 1
#   s = solve(H,ej)
#   x[1 + Hinv3.app@p[j]:(Hinv3.app@p[j+1]-1)] = s[keep[[j]]]
# }
# Hinv3.app@x <- x
# 
# image(Hinv3.app)
# norm(Hinv1-Hinv3.app)/norm(Hinv1)
# ####################################################


#============================#
##### test for 2 tracks #######
#============================#

###############################
###variables initialization ###
###############################

n <- 1000
x_data <- seq(0, 1, length.out = n)
mu <- sin(seq(-6 , 6, length.out = n)) +1
y1 <- rnbinom(length(x_data), mu=exp(mu), size=1)
y2 <- rnbinom(length(x_data), mu=exp(mu), size=1)

list_y <- list(y1,y2)
nb_var <- 2
nb_beta <- 2

nb_knots <- 500
ord<-2
tol <- 0.0001

lambda <- 50
theta <- 1

offset0 <- rep(0, length(y1))
list_offsets <- list(offset0,offset0)

struct_design_matrix <- matrix(c(1,1,0,1),2,2)

distribution <- "negbin"

#############################
####implemented function##### 
#############################

n <- length(x_data)

# Generate x_knots 
m <- ord +1 
nk <- nb_knots - ord
xu <- max(x_data)
xl <- min(x_data)
xr <- xu - xl
xl <- xl - xr * 0.001
xu <- xu + xr * 0.001
dx <- (xu - xl)/(nk - 1)
x_knots <- seq(xl - dx * (m), xu + dx * (m), length = nk + 2 * m)
nb_knots_upd <- length(x_knots)

# Initialisation of variables
if(nb_var > 1) {
  y <- NULL
  offset <- NULL
  for(var in 1:nb_var){
    y <- c(y, list_y[[var]])
    offset <- c(offset,list_offsets[[var]])
  }}
beta <-matrix(mean(y),nb_beta*nb_knots,1)

# Create design matrices
x <- as(.bspline(x_data, x_knots, ord),"dgTMatrix")
if (nb_var==1){
  X <- x
}else{
  X <- Matrix(0,nb_var*n,nb_beta*nb_knots)
  X <- as(X,"dgTMatrix")
  for(var in 1:nb_var){
    for(b in 1:nb_beta){
      if(struct_design_matrix[var,b]==1){
        X@i <- c(X@i, x@i+as.integer((var-1)*n))
        X@j <- c(X@j, x@j+as.integer((b-1)*nb_knots))
        X@x <- c(X@x,x@x)
      }
    }
  }
}
X <- as(X,"CsparseMatrix")
image(X)

# Compute the second differences matrix of beta (for penalized term)
S<- scnd_diff_beta(dim(X)[2],nb_beta)
image(S)
image(S[(ceiling(dim(S)[2]/2)-50):(ceiling(dim(S)[2]/2)+50),(ceiling(dim(S)[2]/2)-50):(ceiling(dim(S)[2]/2)+50)])

# Compute beta 
system.time(
if (distribution=="gaussian"){
  res.opt <- optim(beta,likelihood_penalized_gaussian_log,gradient_likelihood_penalized_gaussian_log,
                   X=X,y=y,lambda = lambda, S=S,method="L-BFGS",control=list(fnscale=-1))
}else if (distribution=="negbin"){
  res.opt <- optim(beta,likelihood_penalized_negbin_log,gradient_likelihood_penalized_negbin_log,
                   X=X,y=y,offset = offset, theta = theta,lambda = lambda, S=S, method="L-BFGS",control=list(fnscale=-1))
  
})
beta <- res.opt$par

plot(x_data, log(y1),ylim=c(-2, max(log(y))))
lines(x_data, mu, col='green', lwd = 2)
lines(x_data, X[1:n,1:nb_knots]%*%beta[1:nb_knots] , col='red',lty=4, lwd = 2)

#   Compute the Hessian matrix
system.time(
if (distribution=="gaussian"){
  H <- compute_hessian_gaussian(beta,X,y,lambda,S)
}else if (distribution=="negbin"){
  H <- compute_hessian_negbin(beta,X,y,offset=offset,theta,lambda,S)
}
)
image(H)


# ############ FIRST ATTEMPT : SOLVE ###################
# Computation of the inverse of H
Hinv1 <- solve(H)
Hinv1.app <- Hinv1
Hinv1.app[abs(Hinv1.app) < tol] <- 0
image(Hinv1.app)

norm(Hinv1-Hinv1.app)/norm(Hinv1)
# ######################################################  

##############SECOND ATTEMPT : CUTOFF###################

# C = chol(-H, pivot=TRUE)
C <- Cholesky(-H)

system.time(
  i_x <- lapply(
    1:ncol(H),
    function(j){
      ej = rep(0, nrow(H))
      ej[j] = 1
      s = solve(C,ej)
      keep = which(abs(s)>tol)
      list(i=keep, x=s[keep])
    }
  )
)

Hinv2.app = sparseMatrix(
  i = unlist(lapply(i_x, function(l) l$i)),
  j = rep(1:nrow(H), times=sapply(i_x, function(l) length(l$i))),
  x = unlist(lapply(i_x, function(l) -l$x))
)
image(Hinv2.app)

norm(Hinv1-Hinv2.app)/norm(Hinv1)
#####################################################

##############THIRD ATTEMPT : FOOTPRINT #################
#Create a sparsity pattern for Hinv
ll <- length(which(X[,ceiling(ncol(H)/2)]!=0))
aux <- seq(1:((nb_spl-1)*ll))
keep <- lapply(1:ncol(H), function(i){
  ind <- which(X[,i]!=0)
  list_ind <- c(ind, aux+ind[length(ind)])
  ifelse(list_ind < nrow(X),list_ind,nrow(X))
})

big_X = sparseMatrix(
  i = unlist(keep),
  j = rep(1:ncol(H), times=sapply(keep, length)),
  x = 1
)

tX_X <- t(big_X)%*%big_X

ind <- which(X[,50]!=0)
keep <- lapply(1:ncol(H), function(i) which(tX_X[,i]!=0))

Hinv3.app = sparseMatrix(
  i = unlist(keep),
  j = rep(1:ncol(H), times=sapply(keep, length)),
  x = 1
)

Compute the approximation of the inverse of H
x = Hinv3.app@x
for(j in 1:ncol(H)){
  ej = re## p(0, nrow(H))
  ej[j] = 1
  s = solve(H,ej)
  x[1 + Hinv3.app@p[j]:(Hinv3.app@p[j+1]-1)] = s[keep[[j]]]
}
Hinv3.app@x <- x

image(Hinv3.app)

norm(Hinv1-Hinv3.app)/norm(Hinv1)
####################################################

###############
# Cholesky sparse factorization
# + trick from Wood, JASA, 2016
###############
n = 1e4
k = 2
u = c(1:k, k+1, k:1)

m = floor(n/length(u))
j = rep(1:m, each=length(u))
i = rep(1:length(u),m) + j - 1

X0 = sparseMatrix(i=i, j=j, x=rep(u,m))
Xleft = rbind(X0,X0)
Xright = rbind(sparseMatrix(i=NULL, j=NULL, x=0, dims=dim(X0)), X0)
X = cbind(Xleft,Xright)

D = Diagonal(x=exp(rnorm(nrow(X), mean=1, sd=0.5)))

A <- t(X)%*%D%*%X

lambda=1e-3
H = A + lambda*Diagonal(x=rep(1,nrow(A)))
H_chol = Cholesky(H, LDL = FALSE)
P = as(H_chol, "pMatrix")

C =  cbind(X0, sparseMatrix(i=NULL, j=NULL, x=0, dims=dim(X0)))

P_Ct = P%*%t(C)

## slow but one could handle more columns at the same time.
sigma = sapply(
  1:ncol(P_Ct),
  function(i){
    L_mt_P_Ct_i = solve(H_chol, P_Ct[,i], "L") # Wood uses the other convention L <-> L^T
    sum(L_mt_P_Ct_i**2)
  }
)

## compare against inversion then extracting the diagonal
## this can be faster for small n's but has a n**2 memory footprint
Sigma = C%*%solve(H_chol, t(C))
plot(sigma,diag(Sigma))
