plot_ggplot2 <- function(df, scale = TRUE, pages = 1,
                      select = NULL) {
    
}

plot_base <- function(df, scale = TRUE, pages = 1,
                      select = NULL) {
    
    numPlots <- length(select)
    x <- df$pos
    y <- df[,select, drop = FALSE]
    upper <- y + 1.96*df[,paste("se", select, sep = "."), drop = FALSE]
    lower <- y - 1.96*df[,paste("se", select, sep = "."), drop = FALSE]
    xlab <- paste("Genomic poisition on", as.character(unique(df$seqnames)))

    if(scale) {
        ylim <- c(min(lower), max(upper))
    }
    
    if(pages <= 1) {
        par(mfrow = c(numPlots, 1))
    }
    for(track in select) {
        if(pages > 1) {
            X11()
        }
        if(!scale) {
            ylim = range(lower[[track]], upper[[track]])
        }
        plot(x, y[[track]], ylim = ylim, xlab = xlab, ylab = "log-fit", type = "l",
             main = track)
        lines(x, upper[[track]], lty = "dotted")
        lines(x, lower[[track]], lty = "dotted")
        abline(h = 0)
    }
}
