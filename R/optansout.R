optansout <- function(ansdf, filename) {
    ##### OPEN ISSUES: (any date order)
    
    ##### IMPLEMENTED: (reverse date order)
    
    #  A funtion to display and print to file (if present) the output of optimx
    if (!is.null(filename)) {
        sink(filename)
    }
    #  if (! exists(filename)) { sink(filename) } # may need to check paths
    tpar <- ansdf$par
    tdf <- ansdf
    ltvec <- length(tpar[[1]])
    for (i in 1:length(tpar)) {
        tvec <- tpar[[i]]
        ltvec <- length(tvec)
        if (ltvec > 5) {
            tvec <- tvec[1:5]
        }
        tpar[[i]] <- tvec
    }
    tdf$par <- tpar
    pardf <- tdf[1]
    if (ltvec > 5) {
        names(pardf)[1] <- "first.5.par"
    }
    else {
        names(pardf)[1] <- "par"
    }
    print(pardf)
    tdf$par <- NULL
    print(tdf)
    if (!is.null(filename)) {
        sink()
    }
    # if (! exists(filename)) { sink() } # to turn it off
    rm(tdf)
    rm(tpar)
    return(TRUE)
}

