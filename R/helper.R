#' get varibles from a formula object
#' @noRd
.getVars <- function(formula, type = c("by", "covar")) {
    type <- match.arg(type)
    chFormula <- as.character(formula)
    variables <- gsub(" ", "", strsplit(chFormula[length(chFormula)],"\\+")[[1]])
    r <- switch(type,
                by = "by=(.+)\\)",
                covar = "s\\((.+)\\)")
    if(type == "covar") variables <- gsub(",.+", ")", variables)
    temp <- sapply(sapply(variables, function(m) regmatches(m,regexec(r,m))), function(y) y[2])
    res <- gsub(",(.+)", "", temp)
    return(res)
}

#' makes column names from formula
#' @noRd
.makeNames <- function(formula) {
    vars <- .getVars(formula)
    vars[!is.na(vars)] <- paste("s(x)", vars[!is.na(vars)], sep = ":")
    vars[is.na(vars)] <- "s(x)"
    return(vars)
}

#' helper function to help fill a list of parameters
#' with the defaults if some are missing
#' @noRd
.fillParameters <- function(l, ...) {
    params <- c(...)
    allin <- names(params) %in% names(l)

    ## has to be done after the 'allin', since if all params in l
    ## are not valid, %in% would match an empty list and not fill
    ## it with the default values
    wrong_params <- names(l) %in% names(params)
    if(sum(!wrong_params) > 0) {
        futile.logger::flog.warn("Some supplied parameters aren't valid and won't be used")
        l <- l[wrong_params]
    }
    
    if(!all(allin)) {
        for(elem in names(params)) {
            if(is.null(l[[elem]])) {
                l[[elem]] <- params[[elem]]
            }
        }
    }
    return(l)
}

#' extract GenomicRanges from list of GPos objects
#' @noRd
.extractGR <- function(gp) {
    res <- GenomicRanges::GRanges(S4Vectors::runValue(GenomeInfoDb::seqnames(gp)), IRanges::ranges(gp)@pos_runs)
    GenomeInfoDb::seqinfo(res) <- GenomeInfoDb::seqinfo(gp)
    return(res)
}

#' helper function to subset by coordinates
#' @noRd
.subsetByCoords <- function(x, i) {
    
    if(length(i) > (2^31 - 1)) {
        stop("Subset to big for short integer vector. Please select a smaller region")
    }

    ## dealing with list objects (mostly from GenoGAMDataSetList)
    if(is(x, "list")) {
        ## build mapping function
        len <- cumsum(as.numeric(sapply(x, NROW)))
        .map <- stepfun(x = len, y = 1:(length(len) + 1), right = TRUE)

        ## map to correct list elements and rows
        ## make by name, as the list elements are supposed to be
        ## named in this application
        mapped_i <- .map(i)
        table_i <- table(mapped_i)

        ## initialize return object
        res <- vector("list", length(table_i))
            
        ## retrieve data
        for(jj in 1:length(table_i)) {
            idx <- as.integer(names(table_i)[jj])
            begin <- i[1] - max(len[idx - 1], 0)
            rows <- begin:(begin + table_i[jj] - 1)
            res[[jj]] <- x[[idx]][rows, ]
            i <- i[-(1:table_i[jj])]
        }

        ## check if everything went correct
        if(length(i) != 0) {
            warning("Subsetting by coordinates did not finish correctly, check the function or the coordinates")
        }

        res <- switch(class(res[[1]]),
                      GPos = do.call("c", res),
                      GRanges = do.call("c", res),
                      data.frame = do.call("rbind", res),
                      DataFrame = do.call("rbind", res),
                      DelayedMatrix = do.call("rbind", res))
    }

    ## for all non-list objects
    else {
        res <- x[i, ]
    }
    return(res)
}
