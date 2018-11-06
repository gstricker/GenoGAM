context("Test the GenoGAMSetup functionality")

test_that("GenoGAMSetup constructor works correctly", {
    ggs <- GenoGAMSetup()

    expect_true(all.equal(slot(ggs, "params"),
                          list(lambda = 0, theta = 0, eps = 0,
                               order = 2, penorder = 2)))
    expect_true(length(slot(ggs, "knots")) == 0)
    expect_true(length(slot(ggs, "designMatrix")) == 0)
    expect_true(length(slot(ggs, "beta")) == 0)
    expect_true(length(slot(ggs, "se")) == 0)
    expect_true(length(slot(ggs, "penaltyMatrix")) == 0)
    expect_true(slot(ggs, "formula") == formula(~ 1))
    expect_true(length(slot(ggs, "offset")) == 0)
    expect_true(is(slot(ggs, "family"), "GenoGAMFamily"))
    expect_true(length(slot(ggs, "response")) == 0)
    expect_true(length(slot(ggs, "fits")) == 0)
    expect_true(length(slot(ggs, "control")) == 6)

    mat <- matrix(1:9, 3, 3)
    k <- 10
    control <- list(eps = 1e-6, maxiter = 1000, alpha = 1, rho = 0.5, c = 1e-4,
                    m = 6)
    ggs <- GenoGAMSetup(params = list(lambda = 5, order = 4, eps = 0.05),
                        knots = list(1:k), designMatrix = as(mat, "dgCMatrix"),
                        beta = mat, se = list(runif(k)),
                        penaltyMatrix = as(mat, "dgCMatrix"), formula = ~ s(x),
                        offset = 1:3, family = GenoGAMFamily(), response  = 1:k,
                        fits = list(runif(k)*3))

    expect_true(all.equal(slot(ggs, "params"),
                          list(lambda = 5, order = 4, eps = 0.05,
                               theta = 0, penorder = 2)))
    expect_true(all(slot(ggs, "knots")[[1]] == 1:10))
    expect_true(all(dim(slot(ggs, "designMatrix")) == c(3, 3)))
    expect_true(all(dim(slot(ggs, "beta")) == c(3, 3)))
    expect_true(length(slot(ggs, "se")[[1]]) == k)
    expect_true(all(dim(slot(ggs, "penaltyMatrix")) == c(3, 3)))
    expect_true(slot(ggs, "formula") == formula(~ s(x)))
    expect_true(all(slot(ggs, "offset") == 1:3))
    expect_true(is(slot(ggs, "family"), "GenoGAMFamily"))
    expect_true(all(slot(ggs, "response") == 1:10))
    expect_true(length(slot(ggs, "fits")[[1]]) == k)
    expect_true(length(slot(ggs, "control")) == length(control))
})

test_that("Knots are placed correctly", {
    x <- 1:1000
    nknots <- 100
    ord <- 2
    knots <- .placeKnots(x, nknots, ord)

    expect_true(length(knots[knots < 1]) == 4)
    expect_true(length(knots[knots > 1000]) == 4)
    expect_true(length(knots) == (nknots + 2*ord))
    ## check if all knots have same spacing up to a tolerance value
    lowerRange <- min(diff(knots))
    upperRange <- max(diff(knots))
    expect_true(all.equal(lowerRange, upperRange))
})

test_that("The penalty matrix is build correctly", {
    p <- 10
    order <- 2

    S <- .buildSMatrix(p, order)
    
    diagonal <- rep(c(1, 5, rep(6, 6), 5, 1), 2)

    expect_true(all(Matrix::diag(S) == diagonal))
    expect_true(all(Matrix::rowSums(S) == 0))
    expect_true(all(Matrix::colSums(S) == 0))

    I <- .buildIMatrix(p, 0.5)
    expect_true(all(Matrix::diag(I) == rep(0.5, p)))
})

test_that("The design matrix is build correctly", {
    design <- matrix(c(1,1,1,0), 2, 2)
    template <- matrix(1, 3, 3)
    X <- .blockMatrixFromDesignMatrix(template, design)

    expect_true(all(Matrix::diag(X) == rep(c(1, 0), each = 3)))
    expect_true(all(Matrix::colSums(X) == c(6, 6, 6, 3, 3, 3)))
    expect_true(all(Matrix::rowSums(X) == c(6, 6, 6, 3, 3, 3)))

    ggd <- makeTestGenoGAMDataSet(sim = TRUE)
    k <- 10
    order <- 2
    x <- pos(rowRanges(ggd))
    knots <- .placeKnots(x, k, order)
    X <- .buildDesignMatrix(knots, x, order)
    dims <- dim(ggd)
    
    expect_true(all.equal(Matrix::rowSums(X), rep(1, each = dims[1])))

    design(ggd) <- ~ s(x)
    X <- .buildDesignMatrix(knots, x, order)
    dimsX <- dim(X)
    expect_true(all(dimsX == c(dims[1], k)))

    design(ggd) <- ~ s(x, by = experiment)
    X <- .buildDesignMatrix(knots, x, order)
    dimsX <- dim(X)
    expect_true(all(dimsX == c(dims[1], k)))
    expect_true(sum(X[1:dims[1],]) == 10000)
})
