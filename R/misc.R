## The IRLS algorithm
.irls_nb <- function(beta0, X, y, theta, lambda, S, offset, control = list(eps = 1e-6, maxiter = 1000)) {
    beta <- beta0
    beta_old <- beta - control$eps - 1
    dif <- max(abs(beta - beta_old))
    k <- 0
    converged <- TRUE
    
    while (dif > control$eps & k < control$maxiter){
        futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))

        eta <- offset + X %*% beta
        mu <- exp(eta)
        z <- 1/mu*(y - mu)+eta
  
        ## 1. Compute V 
        v <- as.numeric(mu^2/Matrix::diag(mu+(mu)^2/theta))
    
        ## 2. Update weight matrix
        W <- Matrix::bandSparse(dim(X)[1], k = 0, diag = c(list(v)))
        
        ## 3. Compute new value of beta
        beta_old <- beta
        H <- Matrix::t(X) %*% W %*% X + 2*lambda * S
        beta <- Matrix::solve(H, Matrix::t(X) %*% W %*% z)

        dif <- max(abs(beta - beta_old))

        futile.logger::flog.debug(paste0("Maximal Parameter Difference: ", dif, "\n"))
        futile.logger::flog.debug(paste0("---------------------------\n"))
    
        k <- k+1
    }

    if(k == control$maxiter) {
        converged <- FALSE
    }
    res <- list(par = as.numeric(beta), converged = converged, iterations = k + 1)
    return(res)
}


.computeDirection <- function(H, gradk, s, y, ro, len) {
    q <- gradk
    if(len == 0) {
        r <- Matrix::solve(H, q)
        return(r)
    }

    alpha <- numeric(len)
    
    for(ii in len:1) {
        alpha[ii] <- ro[ii] * crossprod(s[,ii], q)
        q <- q - alpha[ii] * y[,ii]
    }

    r <- Matrix::solve(H, q)@x
    
    for(jj in 1:len) {
        beta <- as.numeric(ro[jj] * crossprod(y[,jj], r))
        r <- r + s[,jj] * (alpha[jj] - beta)
    }
    return(r)
}

## the simple backtrack algorithm for linesearch
## to ensure Wolfe conditions. See Nocedal, Algorithm 3.1
.linesearch <- function(x, p, f, grad, alpha0 = 1, c = runif(1), rho = 0.9, ...) {
    alpha <- alpha0
    
    ## check if initial alpha is good enough
    betas <- as.numeric(x + alpha * p)
    fleft <- f(beta = betas, ...)
    fright <- as.numeric(f(beta = x, ...) + c * alpha * (Matrix::t(grad(beta = x, ...)) %*% p))

    ## decrease alpha by rho, to 'backtrack' it's optimal values
    while(fleft > fright) {
        alpha <- rho * alpha
        betas <- as.numeric(x + alpha * p)
        fleft <- f(beta = betas, ...)
        fright <- as.numeric(f(beta = x, ...) + c * alpha * (Matrix::t(grad(beta = x, ...)) %*% p))
    }
    return(alpha)
}

#' The modified lbfgs function for sparse optimization problems
#' @param x0 an initial vector
#' @param H0 an initial Hessian matrix. Can be missing, in which case
#' it will be the identity matrix
#' @param f the likelihood function
#' @param gr the likelihood gradient function
#' @param control parameters, if empty filled accordingly
#' @details Performs a modified L-BFGS. The modification stems largely from the fact,
#' that instead of an diagonal approximated inverse Hessian an approximated Hessian is
#' used and a linear system is solved via a direct solver is solved to obtain the
#' respective values. (see Nocedal, p. 178)
#' @noRd
.lbfgs <- function(x0, H0, f, gr, control = list(), fact = c(1,1), ...) {
    ## If list is empty then replace with the proper settings
    if(length(control) == 0) {
        control <- list(eps = 1e-6, maxiter = 1000, alpha = 1, rho = 0.5, c = 1e-4, m = 6)
    }

    ## initialize Hessian as identity matrix
    if(missing(H0)) {
        H0 <- Matrix::bandSparse(length(x0), k = 0)
    }
    
    ## initialize variables
    k <- 1
    idx <- 0
    x <- x0
    H <- H0
    epsilon <- control$eps
    m <- control$m
    converged <- TRUE

    ## the only vectors stored.
    s <- matrix(nrow = length(x0), ncol = m)
    y <- matrix(nrow = length(x0), ncol = m)
    ro <- numeric(m)

    lllast <- 0
    llk <- epsilon + 1
    lldiff <- abs(llk - lllast)
    normGrad <- epsilon + 1
    grad <- matrix(gr(beta = x, ...), ncol = 1)
    
    ## Start L-BFGS loop
    while(lldiff > epsilon & normGrad > epsilon & k <= control$maxiter) {

        ## compute direction and stepsize
        p <- (-1) * .computeDirection(H, grad, s, y, ro, idx)
        
        alphak <- .linesearch(x, p, f, gr, alpha0 = control$alpha, c = control$c, rho = control$rho, ...)

        ## update params
        idx <- min(k, m) ## are we within the first m iterations?
        xnext <- matrix(x + alphak*p, ncol = 1)

        ## update ll and gradient
        lllast <- llk
        llk <- f(xnext, fact = fact, ...)
        lldiff <- abs(llk - lllast)
        gradNext <- matrix(gr(beta = xnext, ...), ncol = 1)
        normGrad <- sqrt(as.numeric(crossprod(gradNext)))

        ## update s, y and rho vectors
        ## if 
        if(k > m) {
            s <- cbind(s[,-1], 0)
            y <- cbind(y[,-1], 0)
            ro <- ro[-1]
        }

        s[,idx] <- xnext - x
        y[,idx] <- gradNext - grad
        ro[idx] <- as.numeric(1/crossprod(y[,idx], s[,idx]))

        

        ## create new Hessian
        H <- .Hupdate(xnext, ...)

        ## reset params variables
        x <- xnext
        grad <- gradNext

        ## printing
        futile.logger::flog.debug(paste0("Iteration: ", k, "\n"))
        futile.logger::flog.debug(paste0("||gr|| = ", normGrad, "\n"))
        futile.logger::flog.debug(paste0("Alpha = ", alphak, "\n"))
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
