## NOTE: All estimation and optimization is performed on the negative log-likelihood,
## as we are targeting a minimization problem (as is usually done in numerical optimization).
## The log-likelihood, the gradient and the hessian should usually be pre-multiplies by (-1)
## to get their normal form as derived on paper.

## estimate betas
.estimateParams <- function(ggs) {

    betas <- slot(ggs, "beta")

    distr <- slot(ggs, "family")
    X <- slot(ggs, "designMatrix")
    y <- slot(ggs, "response")
    offset <- slot(ggs, "offset")
    params <- slot(ggs, "params")
    S <- slot(ggs, "penaltyMatrix")
    control <- slot(ggs, "control")

    if(length(ggs) == 0) {
        res <- list(par = matrix(0,0,0), converged = FALSE, iterations = 0)
        return(res)
    }

    if(sum(slot(ggs, "response")) == 0){
        par <- rep(0, nrow(betas))
        res <- list(par = matrix(par, nrow = length(par), ncol = 1),
                    converged = TRUE, iterations = 0)
        matrix(res$par, nrow = length(res$par), ncol = 1)
    }
    else {
        H <- .compute_hessian(family = distr, beta = betas, X = X,
                              XT = Matrix::t(X), offset = unname(offset),
                              y = y, S = S, lambda = params$lambda,
                              theta = params$theta)

        res <- .newton(x0 = betas, H0 = H, fam = distr, X = X, y = y,
                       offset = offset, theta = params$theta,
                       lambda = params$lambda, S = S, fact = dim(X),
                       control = control)
    }
        
    return(res)
}

.compute_SE <- function(setup){
    if(length(setup) == 0) {
        return(numeric())
    }
    
    params <- slot(setup, "params")
    theta <- params$theta
    lambda <- params$lambda
    offset <- slot(setup, "offset")
    betas <- slot(setup, "beta")
    X <- slot(setup, "designMatrix")
    y <- slot(setup, "response")
    S <- slot(setup, "penaltyMatrix")
    distr <- slot(setup, "family")
    des <- slot(setup, "design")
  
    H <- .compute_hessian(family = distr, beta = betas, X = X,
                          XT = Matrix::t(X), offset = unname(offset),
                          y = y, S = S, lambda = params$lambda,
                          theta = params$theta)

    Hinv <- .invertHessian(H)

    ## get indeces of the design matrix and the Hessian
    ## by row and column according to the design
    
    nSplines <- ncol(des)
    nSamples <- nrow(des)
    dims <- dim(X)
    rows <- list(1:dims[1])
    cols <- list(1:dims[2])

    if(nSplines > 1) {
        cols <- split(1:dims[2], cut(1:dims[2], nSplines, labels = 1:nSplines))
    }
    if(nSamples > 1) {
        rows <- split(1:dims[1], cut(1:dims[1], nSamples, labels = 1:nSamples))
    }

    ## compute the standard error by column (splines) and row (samples)
    ## of the design matrix and the Hessian. All standard errors of the
    ## same spline are combined to a vector to have the same structure as
    ## the response. Each spline is one element of the list, that
    ## gets returned as the result.
    ses <- lapply(1:nSplines, function(j) {
        res <- sapply(1:nSamples, function(i) {
            Xsub <- X[rows[[i]], cols[[j]]]
            HinvSub <- Hinv[cols[[j]],cols[[j]]]
            se <- compute_stdError(Xsub, HinvSub)
            if(length(se@x) != 0) {
                return(se@x)
            }
            else {
                return(rep(0L, nrow(Xsub)))
            }
        })
        res <- as.vector(res)
        return(res)
    })
    varNames <- .makeNames(slot(setup, "formula"))
    names(ses) <- varNames

    return(ses)
}

## Compute the penalized Hessian
.compute_hessian <- function(family, ...){
    
    if(is.na(slot(family, "name"))) {
        res <- as(matrix(,0, 0), "dgCMatrix")
    }

    res <- compute_pen_hessian(hessid = slot(family, "hessian"), ...)
    
    return(res)
}

## Computation of the inverse of H
.invertHessian <- function(H) {
    if(length(H) == 0) {
        res <- as(matrix(,0, 0), "dgCMatrix")
    }
    else {
        res <- sparseinv::Takahashi_Davis(H)
    }
    return(res)
}

.Hupdate <- function(family, x, X, XT, ...) {
    
    H <- compute_pen_hessian(hessid = slot(family, "hessian"), beta = x, X = X, XT = XT, ...)
    return(H)
    
}

.newton <- function(x0, H0, fam, X, control = list(), fact = c(1,1), ...) {
    ## If list is empty then replace with the proper settings
    if(length(control) == 0) {
        control <- list(eps = 1e-6, maxiter = 1000)
    }

    ## initialize Hessian as identity matrix
    if(missing(H0)) {
        H0 <- Matrix::bandSparse(length(x0), k = 0)
    }

    ## set Family
    f <- slot(fam, "ll")
    gr <- slot(fam, "gradient")
    
    ## initialize variables
    k <- 1
    x <- x0
    H <- H0
    XT <- Matrix::t(X)
    epsilon <- control$eps
    converged <- TRUE

    lllast <- 0
    llk <- epsilon + 1
    lldiff <- abs(llk - lllast)
    normGrad <- epsilon + 1
    grad <- gr(beta = x, X = X, XT = XT, ...)
    
    ## Start L-BFGS loop
    while(lldiff > epsilon & normGrad > epsilon & k <= control$maxiter) {

        ## compute direction
        p <- (-1) * Matrix::solve(H, grad)@x
        
        ## update params
        xnext <- matrix(x + p, ncol = 1)

        ## update ll and gradient
        lllast <- llk
        llk <- f(xnext, X = X, ll_factor = fact[1], lambda_factor = fact[2], n = fact[1], ...)
        lldiff <- abs(llk - lllast)
        gradNext <- gr(beta = xnext, X = X, XT = XT, ...)
        normGrad <- sqrt(as.numeric(crossprod(gradNext))) 

        ## create new Hessian
        H <- .Hupdate(family = fam, x = xnext, X = X, XT = XT, ...)

        ## reset params variables
        x <- xnext
        grad <- gradNext

        ## printing
        futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))
        futile.logger::flog.debug(paste0("||gr|| = ", normGrad, "\n"))
        futile.logger::flog.debug(paste0("Likelihood = ", llk, "\n"))
        futile.logger::flog.debug("---------------------------\n")

        ## update iteration
        k <- k + 1
    }

    if(k == control$maxiter) {
        converged <- FALSE
    }
    return(list(par = x, converged = converged, iterations = k))
}
