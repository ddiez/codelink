# normalize()
# wrapper function to apply normalization methods to Codelink objects.
normalize <- function(object, method = "quantiles", log.it = TRUE,
	preserve = FALSE, weights = NULL, verbose = FALSE) {

	if(!is(object,"Codelink")) stop("A Codelink object is needed.")
	if(is.null(object$Ri)) stop("Background corrected intensities needed.")
	if(log.it & object$method$log) stop("Intensities already log2.")
	
	method <- match.arg(method, c("loess", "quantiles", "median"))

	object$Ni <- object$Ri
	if(log.it) object$Ni <- log2(object$Ni)
	switch(method,
		loess = {
			object$Ni <- normalize.loess(object$Ni, log.it = FALSE,
				weights = weights, verbose = verbose)
			object$method$normalization <- "CyclicLoess"
		},
		quantiles = {
			object$Ni <- normalizeQuantiles(object$Ni)
            object$method$normalization <- "quantiles"
		},
		median = {
			# no weights used right now.
			object$Ni <- normalize.median(object$Ni)
            object$method$normalization <- "median"
		}
	)
	if(!preserve) object$Ri <- NULL
	if(log.it) object$method$log <- TRUE
	return(object)
}

# normalize.median()
# based on limma implementation.
normalize.median <- function(x, weights = NULL) {
	l <- dim(x)[2]
	
	if(is.null(weights)) {
		for(j in 1:l)
			x[, j] <- x[, j] - median(x[, j], na.rm = TRUE)
	} else {
		for(j in 1:l)
			x[, j] <- x[, j] - weighted.median(x[, j], weights[, j],
				na.rm = TRUE)
	}

	x
}

# normalize.loess()
# based on affy implementation.
# modified to allows missing values and weights.
normalize.loess <- function (mat, 
	subset = sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))),
	epsilon = 10^-2, maxit = 1, log.it = TRUE,
	verbose = FALSE, span = 2/3, family.loess = "symmetric", weights = NULL) 
{
    J <- dim(mat)[2]
    II <- dim(mat)[1]
    newData <- mat
    if (log.it) {
        mat <- log2(mat)
        newData <- log2(newData)
    }
    change <- epsilon + 1
    fs <- matrix(0, II, J)
    iter <- 0

	if(is.null(weights))
    	w <- c(0, rep(1, length(subset)), 0)
	else
		w <- weights
    
	while (iter < maxit) {
        iter <- iter + 1
        means <- matrix(0, II, J)
        for (j in 1:(J - 1)) {
            for (k in (j + 1):J) {
                y <- newData[, j] - newData[, k]
                x <- (newData[, j] + newData[, k])/2
				
				# select genes that are not set to NA
				sel <- which(!is.na(as.character(y)))
				y <- y[sel]
				x <- x[sel]
				ww <- w[sel]
                
				index <- c(order(x)[1], subset, order(-x)[1])
                xx <- x[index]
                yy <- y[index]
				# reorder weights.
				ww <- ww[index]

                aux <- loess(yy ~ xx, span = span, degree = 1, 
                  weights = ww, family = family.loess)
                aux <- predict(aux, data.frame(xx = x))/J
				# apply normalization to genes not NA.
                means[sel, j] <- means[sel, j] + aux
                means[sel, k] <- means[sel, k] - aux
                if (verbose) 
                  cat(" + Done with", j, "vs", k, " in iteration ", 
                    iter, "\n")
            }
        }
        fs <- fs + means
        newData <- mat - fs
        change <- max(colMeans((means[subset, ])^2))
        if (verbose) 
            cat(iter, change, "\n")
        oldfs <- fs
    }
    if ((change > epsilon) & (maxit > 1)) 
        warning(paste("!! No convergence after", maxit, "iterations.\n"))
    if (log.it) {
        return(2^newData)
    }
    else return(newData)
}

