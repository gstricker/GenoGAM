gr <- GRanges(c("chrI", "chrI", "chrII", "chrIII", "chrII"),
              IRanges(start = c(1, 10000, 1000, 1000, 100000),
                      end = c(9234, 100001, 23456, 23456, 120123)))
subggd <- ggd[gr]
(if(all(width(getIndex(subggd)) == tileSettings(ggd)$tileSize)) {
     cat("SUCCESFULL\n")
 }
 else{
     cat("ERROR\n")
 })

