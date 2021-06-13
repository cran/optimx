## optimx-package.R 
## This file contains support routines (methods) for the optimx() function
##################################################################
summary.opm <- summary.optimx <- function(object, order = NULL, par.select = TRUE, ...) {
	
	# internally object is referred to as x and par.select as par
	x <- object
	par <- par.select

	# first npar columns of object are the parameters
        npar <- attr(x,"npar")

	if (is.character(par)) {
		idx <- which(par %in% names(x))
	} else if (is.logical(par)) { 
		idx <- which(rep(par, length = npar))
	} else if (is.numeric(par)) {
                idx <- intersect(par, seq_len(npar))
	} else stop("par.select must be character, logical or numeric")

        selidx <- union(idx, (npar+1):ncol(x))
        x <- x[, selidx]

	# xx same as x except:
	# - it has a rownames column
	# - it has a natural column reflecting the input ordering
	# - if objective maximized then value column negated
	xx <- cbind(rownames = rownames(x), x, natural=1:nrow(x))
        maximize <- attr(x,"maximize")
	if (maximize) xx$value <- - xx$value

	# try to evaluate order using standard evaluation
	order.try <- try(order, silent = TRUE)
	# did it work?
	if (is.null(order.try)) order.try <- "natural"
	e <- if (!inherits(order.try, "try-error") && is.character(order.try)) {


		# if all components are names then convert to a string
		if (all(order.try %in% names(xx))) {
			order.try <- paste0("list(", toString(order.try), ")")
		}

		# order.try is now the string representation of an R expression
		# so parse it
		e <- parse(text = order.try)
	} else substitute(order)

	# perform non-standard evaluation (as in transform and subset functions)
	order. <- eval(e, xx, parent.frame())

	# ensure order. is a list
	if (!is.list(order.)) order. <- list(order.)

	o <- do.call(base::order, order.)
	x <- x[o, ]

	# ensure details attribute corresponds to data
	attr(x, "details") <- attr(x, "details")[rownames(x), ]

	x
}

##################################################################
coef.opm <- coef.optimx <- function(object, ...) {
	npar <- attr(object, "npar")
	ix <- seq_len(npar)
        cc <- object[, ix]
        attr(cc,"details") <- NULL
        attr(cc,"maximize") <- NULL
        attr(cc,"npar") <- NULL
##         attr(cc,"follow.on") <- NULL # leave follow.on?
        cc<-as.matrix(cc) # coerce to matrix to accord with other uses 130406
        cc
}

"coef<-" <- function(x, value) UseMethod("coef<-")

"coef<-.opm" <- "coef<-.optimx" <- function(x, value) {
	npar <- attr(x, "npar")
	ix <- seq_len(npar)
	structure(cbind(value, x[, -ix, drop = FALSE]), 
		npar = NCOL(value),
		class = class(x))
}
"[.opm" <- "[.optimx" <- function(x, ...) {
	xx <- NextMethod()
	if (is.data.frame(xx)) {
		details <- attr(x, "details")
		# temporarily convert to data.frame so missing ..1 acts as
		#  in data frames rather than as in matrices
		if (is.matrix(details)) details <- 
			as.matrix(as.data.frame(details)[..1, , drop=FALSE])
		structure(xx,
			details = details,
			maximize = attr(x, "maximize"),
			npar = attr(x, "npar"),
			class = class(x))
	} else xx
}

##################################################################
as.data.frame.opm <- as.data.frame.optimx <- function(x, row.names = NULL, optional = FALSE, ...) {
	result <- do.call(data.frame, as.list(x))
	rownames(result) <- if (is.null(row.names)) rownames(x) else row.names
	result # NOTE: seems "details" are stripped away
}

##################################################################
