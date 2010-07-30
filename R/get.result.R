##################################################################
get.result <- function(optimx.obj, method = NULL, attribute = NULL) {
    # optimx.obj = object returned by `optimx'
    # method = method for which all results are desired (e.g. `spg')
    # attribute = type of result desired for all methods (e.g. `fvalues')
    #
    if (!is.null(attribute) & is.null(method)) {
        allattr <- names(optimx.obj)[!(names(optimx.obj) %in% "method")]
        attrib <- unique(match.arg(attribute, allattr, several.ok = TRUE))
        
        sel.attr <- names(optimx.obj) %in% attrib
        ret <- cbind(optimx.obj["method"], optimx.obj[sel.attr])
        ord <- rev(order(unlist(ret[, 2])))
        ret <- ret[ord, ]
    }
    
    if (is.null(attribute) & !is.null(method)) {
        method <- unique(match.arg(method, unlist(optimx.obj$method), 
            several.ok = TRUE))
        sel.meth <- optimx.obj$method %in% method
        ret <- sapply(optimx.obj, function(x) x[sel.meth])
    }
    
    if (!is.null(attribute) & !is.null(method)) {
        warning("You cannot specify both `method' and `attribute' \n")
        ret <- NULL
    }

    if (is.null(attribute) & is.null(method)) {
        warning("You must specify one of `method' or `attribute' \n")
        ret <- NULL
    }
 
    return(ret)
}

##################################################################

