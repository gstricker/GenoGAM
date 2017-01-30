## TODO:

############
## Functions to initiate and fill the GenoGAMSetup object for estimation

#' The B-Spline function.
#'
#' A function to construct B-Spline bases.
#'
#' @param x The numeric vector of x values at which to evaluate the function
#' @param k A vector of knot positions
#' @param ord The B-Spline basis order as order - 1, e.g. ord = 2 is cubic (default)
#' @param derivative The order of derivative
#' @return A dgCMatrix object with dimensions length(x) * length(k). dgCMatrix
#' is a sparse matrix object from the Matrix package
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
#' @noRd
.bspline <- function(x, k, ord = 2, derivative = 0) {
  res <- splines::spline.des(k, x, ord + 2, rep(derivative,length(x)), sparse=TRUE)$design
  return(res)
}
