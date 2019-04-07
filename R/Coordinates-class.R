## ====================
## Coordinates class
## ====================

#' Coordinates class
#'
#' This class mainly exists to overcome the lack of long integer use in R
#' and the IRanges class in particular. The object is basically a restricted
#' DataFrame object, with a few custom methods to mimic the behaviour of IRanges.
#' 
#' @name Coordinates-class
#' @rdname Coordinates-class
#' @author Georg Stricker \email{georg.stricker@@in.tum.de}
setClass("Coordinates",
         contains = "DataFrame",
         prototype = list(listData = list(start = NULL, end = NULL, width = NULL))
         )


.validateCoordinates <- function(object) {
    check <- all(c("start", "end", "width") %in% names(object))
    if(!check | ncol(object) != 3) {
        return("Invalid Coordinates object")
    }
    if(length(object) == 0) {
        if(any(end(object) < start(object))) {
            return("Start is greater than end")
        }
    }
    NULL
}

S4Vectors::setValidity2("Coordinates", .validateCoordinates)

## Constructor
## ===========

#' Coordinates constructor.
#'
#' Coordinates is the constructor function for the Coordinates-class.
#'
#' @param start The start vector, must be of same length as the end vector
#' @param end The end vector, must be of same length as the start vector
#' @param x The Coordinates object
#' @param value A numeric vector. For replace method in 'start' and 'end'.
#' 
#' @return An object of class Coordinates
#' @name Coordinates
#' @rdname Coordinates-class
Coordinates <- function(start = NULL, end = NULL) {
    if(any(is.null(start), is.null(end))) {
        return(new("Coordinates"))
    }
    df <- S4Vectors::DataFrame(start = as.numeric(start), end = as.numeric(end),
                               width = as.integer(end - start + 1))
    return(new("Coordinates", df))
}

## Methods
## ========

#' @describeIn Coordinates The length method
setMethod("length", "Coordinates", function(x) {
    length(x@listData$start)
})

#' @describeIn Coordinates The nrow method
setMethod("nrow", "Coordinates", function(x) {
    length(x@listData$start)
})

#' @describeIn Coordinates The ncol method
setMethod("ncol", "Coordinates", function(x) {
    length(x@listData)
})

#' @describeIn Coordinates The start accessor
setMethod("start", "Coordinates", function(x) {
    unname(x@listData$start)
})

#' @describeIn Coordinates The end accessor
setMethod("end", "Coordinates", function(x) {
    unname(x@listData$end)
})

#' @describeIn Coordinates The width accessor
setMethod("width", "Coordinates", function(x) {
    unname(x@listData$width)
})


#' @describeIn Coordinates Replacement method for end
setReplaceMethod("end", "Coordinates", function(x, value) {
    x@listData$end <- value
    return(x)
})

#' @describeIn Coordinates Replacement method for end
setReplaceMethod("start", "Coordinates", function(x, value) {
    x@listData$start <- value
    return(x)
})

#' @describeIn Coordinates Replacement method for width
setReplaceMethod("width", "Coordinates", function(x, value) {
    x@listData$width <- value
    return(x)
})

setMethod("rbind", "Coordinates", function(...) {
    args <- list(...)
    dflist <- lapply(args, DataFrame)
    df <- do.call("rbind", dflist)
    res <- new("Coordinates", df)
    return(res)
})

#' IRanges to Coordinates
#' 
#' @name asCoordinates
setAs(from = "IRanges", to = "Coordinates", def = function(from) {
    s <- start(from)
    e <- end(from)
    res <- Coordinates(s, e)
    return(res)
})

#' Coordinates to IRanges
#' 
#' @name asCoordinates
setAs(from = "Coordinates", to = "IRanges", def = function(from) {
    s <- start(from)
    e <- end(from)
    res <- IRanges::IRanges(s, e)
    return(res)
})
