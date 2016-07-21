## n <- 6
## mat <- matrix(1:n, n, 1)
## x <- mat %*% t(mat)
## colnames(x) <- as.character(1:n)
## rownames(x) <- as.character(1:n)
## a <- upper.tri(x, diag = TRUE)
## trix <- x[a]

## Y <- diag(n)
## Y[upper.tri(Y, diag = TRUE)] <- trix
## Y <- Y + t(Y) - diag(diag(Y))
