context("Testing the functionality of the main modelling functions")

ggd <- makeTestGenoGAMDataSet(sim = TRUE)
coords <- .getCoordinates(ggd)
emptyGGD <- GenoGAMDataSet()
emptyCoords <- .getCoordinates(emptyGGD)
ggs <- setupGenoGAM(ggd)
emptyGGS <- GenoGAMSetup()

settings <- GenoGAMSettings()
control <- slot(settings, "estimControl")

test_that("The response vector is build correctly", {
    ## all combinations of different empty inputs
    res1 <- .buildResponseVector(emptyGGD, coords, 1)
    res2 <- .buildResponseVector(emptyGGD, emptyCoords, 1)
    res3 <- .buildResponseVector(ggd, emptyCoords, 1)
    
    expect_true(all(c(length(res1), length(res2), length(res3)) == 0))

    id <- sample(length(coords),1)
    res <- .buildResponseVector(ggd, coords, id)
    correctLength <- width(coords[id,]) * dim(ggd)[2]
    expect_true(length(res) == correctLength)
})

test_that("The specific tile setup initializes correctly", {
    ## all combinations of different empty inputs
    res1 <- .initiate(emptyGGD, ggs, coords, 1)
    params <- slot(res1, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == dim(slot(res1, "designMatrix"))[2])
    expect_true(params$theta == 1)

    res2 <- .initiate(ggd, emptyGGS, coords, 1)
    params <- slot(res2, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)

    res3 <- .initiate(ggd, emptyGGS, emptyCoords, 1)
    params <- slot(res3, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)

    res4 <- .initiate(emptyGGD, ggs, emptyCoords, 1)
    params <- slot(res4, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == dim(slot(res1, "designMatrix"))[2])
    expect_true(params$theta == 1)
    
    res5 <- .initiate(ggd, emptyGGS, emptyCoords, 1)
    params <- slot(res5, "params")
    expect_true(length(slot(res1, "response")) == 0)
    expect_true(all(is.na(slot(res1, "beta"))))
    expect_true(params$lambda == 0)
    expect_true(params$theta == 0)

    ## from bug of empty regions. betas should be initialized as log of 1
    ## instead of mean of the response
    id <- sample(length(coords),1)
    test_ggd <- ggd
    assay(test_ggd)$input <- rep(0, length(test_ggd))
    assay(test_ggd)$IP <- rep(0, length(test_ggd))
    res6 <- .initiate(test_ggd, ggs, coords, id)
    trueBeta <- log(1)
    expect_true(all(slot(res6, "beta") == trueBeta))
})

test_that("Fits are correctly computed", {
    ## empty input
    res1 <- .getFits(emptyGGS)
    expect_true(length(res1) == 0)

    ## only empty betas
    res2 <- .initiate(emptyGGD, ggs, coords, 1)
    expect_true(length(res1) == 0)

    ## non-empty input with betas == 1, for easier checking
    ## fits should be just the row sums of the design matrix
    ## which should be 1.
    test_ggs <- ggs
    design <- .getDesignFromFormula(design(ggd), colData(ggd))
    X <- as(.blockMatrixFromDesignMatrix(slot(test_ggs, "designMatrix"), design), "dgCMatrix")
    slot(test_ggs, "designMatrix") <- X
    
    slot(test_ggs, "beta") <- matrix(1, dim(slot(test_ggs, "designMatrix"))[2], 1)
    res3 <- .getFits(test_ggs)
    X <- as.matrix(slot(test_ggs, "designMatrix"))
    fits1 <- rowSums(X[1:dim(X)[1], 1:(dim(X)[2]/2)])
    fits2 <- rowSums(X[1:dim(X)[1], (dim(X)[2]/2 + 1):dim(X)[2]])
    expect_true(length(res3) == 2)
    expect_true(all.equal(res3[[1]], fits1))
    expect_true(all.equal(res3[[2]], fits2))
    fitNames <- c("s(x)", paste("s(x)", colnames(colData(ggd)), sep = ":"))
    expect_true(all(names(res3) == fitNames))
})

test_that("Beta estimation work correct", {
    ## With empty input
    setup <- .initiate(emptyGGD, emptyGGS, coords, sample(5, 1))
    emptyBetas <- .estimateParams(setup)
    expect_true(length(emptyBetas$par) == 0)
    
    ## With normal input
    setup <- .initiate(ggd, ggs, coords, sample(5, 1))
    betas1 <- .estimateParams(setup)
    betas2 <- .estimateParams(setup)
    
    expect_true(all.equal(betas1$par, betas2$par))
    expect_true(all(c(betas1$convergence, betas2$convergence) == 0))
    expect_true(length(betas1$par) == length(slot(setup, "beta")))
})

test_that("Negative binomial log-likelihood and gradient give correct results", {
    setup <- .initiate(ggd, ggs, coords, sample(5, 1))
    X <- slot(setup, "designMatrix")
    S <- slot(setup, "penaltyMatrix")
    theta <- 1
    lambda <- 0
    y <- matrix(1, dim(X)[1], 1)
    betas <- rep(0, dim(X)[2])
    offset <- rep(0, dim(X)[1])

    ## true loglik with input from above
    loglik <- (2*dim(X)[1]*log(2))/dim(X)[1]
    ## computed loglik
    res <- ll_pen_nb(betas, X, y, offset, theta, lambda, S, dim(X)[1], dim(X)[2], dim(X)[1])
    expect_true(all.equal(loglik, res))
    
    ## return warnings due to inappropriate theta
    invalid <- ll_pen_nb(betas, X, y, offset, theta = 0, lambda, S,
                         dim(X)[1], dim(X)[2], dim(X)[1])
    
    expect_true(is.nan(invalid))

    newY <- matrix(3, dim(X)[1], 1)
    ## true gradient
    grad <- (-1) * colSums(as.matrix(X))
    grad <- matrix(grad, ncol = 1)
    ## computed gradients
    res <- gr_ll_pen_nb(betas, X, Matrix::t(X), newY, offset, theta, lambda, S)
    
    expect_true(all.equal(grad, res))
})

test_that("Hessian matrix computation is correct for empty spline", {
    ## with empty input
    setup <- .initiate(emptyGGD, emptyGGS, coords, sample(5, 1))

    params <- slot(setup, "params")
    theta <- params$theta
    lambda <- params$lambda
    offset <- slot(setup, "offset")
    betas <- slot(setup, "beta")
    X <- slot(setup, "designMatrix")
    y <- slot(setup, "response")
    S <- slot(setup, "penaltyMatrix")
    distr <- slot(setup, "family")

    hess <- compute_pen_hessian(betas, X, Matrix::t(X), offset, y, S, params$lambda, params$theta, 1)
    expect_true(length(hess) == 0)
    
    inv <- .invertHessian(hess)
    expect_true(all.equal(hess, inv))

    ## matrix inversion for a too high tolerance
    inv <- .invertHessian(hess)
    expect_true(length(inv) == 0)

    se <- .compute_SE(setup)
    expect_true(length(se) == 0)
})

## Doesn't make sense to test for one and more splines separately, if the computation
## is automatically done for all at once.

## test_that("Hessian matrix computation is correct for one spline", {
##     ## Check hessian
##     setup <- .initiate(ggd, ggs, coords, sample(5, 1))
##     X <- slot(setup, "designMatrix")
##     slot(setup, "fits") <- list("s(x)" = rep(0, dim(X)[1]))
##     slot(setup, "params")$lambda <- 0
##     slot(setup, "params")$theta <- 1
##     slot(setup, "response") <- rep(-5, dim(X)[1])

##     ## true Hessian with above inputs
##     Htrue <- (-1) * Matrix::t(X) %*% X
##     ## computed Hessian
##     params <- slot(setup, "params")
##     theta <- params$theta
##     lambda <- params$lambda
##     offset <- slot(setup, "offset")
##     betas <- slot(setup, "beta")
##     X <- slot(setup, "designMatrix")
##     y <- slot(setup, "response")
##     S <- slot(setup, "penaltyMatrix")
##     distr <- slot(setup, "family")

##     if(distr == "nb") {
##         f <- .hessian_nb
##         args <- list(y = y, theta = params$theta)
##     }
##     res <- do.call(.compute_hessian, c(list(betas, X, offset, S, params$lambda, f), args))
##     expect_true(all.equal(Htrue@x, res@x))

##     ## Check inversion of hessian
##     inv <- .invertHessian(res)
##     ## check that inverted matrix is covariance matrix
##     ## check for symmetry
##     l <- sort(inv[lower.tri(inv)])
##     u <- sort(inv[upper.tri(inv)])
##     expect_true(all.equal(u,l))
##     ## check for positive-definite
##     ## all eigenvalues are positive
##     e <- eigen(inv)
##     expect_true(all(e$values > 0))

##     ## Test computation of standard errors (SE)
##     slot(ggd, "design") <- ~ s(x)
##     ggs <- setupGenoGAM(ggd)
##     setup <- .initiate(ggd, ggs, coords, sample(5, 1))
    
##     slot(setup, "beta") <- .estimateParams(setup)$par
##     slot(setup, "fits") <- .getFits(setup)

##     se <- .computeSE(setup)

##     ## add all SEs at same position
##     seSums <- unique(se[[1]])
##     ## get the element indices for the first and last 1% of data
##     quants <- round(quantile(1:length(seSums), probs = c(0, 0.01,0.99, 1)))
##     data1 <- quants[1]:(quants[2] - 1)
##     data100 <- (quants[3] + 1):quants[4]

##     ## compute the differences of SE in the 1% quantile of data
##     diff1 <- diff(seSums[data1])
##     ## compute the differences of SE in the last 1% quantile of data
##     diff100 <- diff(seSums[data100])
##     ## compare the sum of the difference. The SE should grow at the borders
##     ## hence have a negative difference throughout the small intervals
##     ## close to the borders.
##     expect_true(sum(diff1) < 0 & sum(diff100) > 0)
## })

test_that("Hessian matrix computation is correct", {

    ## Check hessian
    setup <- .initiate(ggd, ggs, coords, sample(5, 1))
    X <- slot(setup, "designMatrix")
    slot(setup, "beta") <- matrix(rep(0, dim(X)[2]), ncol = 1)
    slot(setup, "params")$lambda <- 0
    slot(setup, "params")$theta <- 1
    slot(setup, "response") <- rep(-5, dim(X)[1])

    ## true Hessian with above inputs
    Htrue <- (-1) * Matrix::t(X) %*% X
    ## computed Hessian
    params <- slot(setup, "params")
    theta <- params$theta
    lambda <- params$lambda
    offset <- slot(setup, "offset")
    betas <- slot(setup, "beta")
    X <- slot(setup, "designMatrix")
    y <- slot(setup, "response")
    S <- slot(setup, "penaltyMatrix")
    distr <- slot(setup, "family")

    res <- compute_pen_hessian(betas, X, Matrix::t(X), offset, y, S,
                               params$lambda, params$theta, 1)    
    expect_true(all.equal(Htrue@x, res@x))

    ## Check inversion of hessian
    ## inv <- .invertHessian(res)

    ## ## check for symmetry
    ## l <- sort(inv[lower.tri(inv)])
    ## u <- sort(inv[upper.tri(inv)])
    ## expect_true(all.equal(u,l))
    ## ## check for positive-definite
    ## ## all eigenvalues are positive
    ## e <- eigen(inv)
    ## expect_true(all(e$values < 0))

    ## Test computation of standard errors (SE)
    slot(ggd, "design") <- ~ s(x) + s(x, by = experiment)
    ggs <- setupGenoGAM(ggd)
    setup <- .initiate(ggd, ggs, coords, sample(5, 1))
    slot(setup, "beta") <- .estimateParams(setup)$par
    slot(setup, "fits") <- .getFits(setup)

    se <- .compute_SE(setup)

    ## add all SEs at same position
    seSums <- unique(se[[1]]) + se[[2]][se[[2]] > 0]
    ## get the element indices for the first and last 1% of data
    quants <- round(quantile(1:length(seSums), probs = c(0, 0.01,0.99, 1)))
    data1 <- quants[1]:(quants[2] - 1)
    data100 <- (quants[3] + 1):quants[4]

    ## compute the differences of SE in the 1% quantile of data
    diff1 <- diff(seSums[data1])
    ## compute the differences of SE in the last 1% quantile of data
    diff100 <- diff(seSums[data100])
    ## compare the sum of the difference. The SE should grow at the borders
    ## hence have a negative difference throughout the small intervals
    ## close to the borders.
    expect_true(sum(diff1) < 0 & sum(diff100) > 0)
})

test_that("Data transformation works correct", {
    ## with empty input
    emptyChunks <- .getChunkCoords(emptyCoords)
    emptyFits <- .transformResults(list(emptyGGS), emptyChunks, what = "fits")
    expect_true(length(emptyFits) == 0)

    ## normal input
    chunks <- .getChunkCoords(coords)
    s <- start(chunks) %% start(coords) + 1
    e <- s + width(chunks) - 1
    relativeChunks <- IRanges::IRanges(s, e)

    id <- sample(5, 1)
    setup <- .initiate(ggd, ggs, coords, id)
    slot(setup, "beta") <- .estimateParams(setup)$par
    slot(setup, "fits") <- .getFits(setup)

    slot(setup, "se") <- .compute_SE(setup)
    slot(setup, "params")$id <- id

    fits <- .transformResults(list(setup), relativeChunks, what = "fits")
    se <- .transformResults(list(setup), relativeChunks, what = "se")

    expect_true(is(fits, "DataFrame") & is(se, "DataFrame"))
    expect_true(nrow(fits) == width(chunks[id,]) & nrow(se) == width(chunks[id,]))

    ## mixed input (should return same data.table as above due to first entry being empty)
    fits <- .transformResults(list(emptyGGS, setup), relativeChunks, what = "fits")
    se <- .transformResults(list(emptyGGS, setup), relativeChunks, what = "se")

    expect_true(is(fits, "DataFrame") & is(se, "DataFrame"))
    expect_true(nrow(fits) == width(chunks[id,]) & nrow(se) == width(chunks[id,]))
})
