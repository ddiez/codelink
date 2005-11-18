## normalize()
# Wrapper function to apply normalization methods to Codelink objects.
normalize <- function(object, method="quantiles", log.it=TRUE, type="ALL") {
	if(!is(object,"Codelink")) stop("A Codelink object is needed.")
	method <- match.arg(method,c("loess","quantiles"))

	type <- match.arg(type,c("ALL","DISCOVERY"))

	switch(type,
		DISCOVERY = {
			sel <- object$info[,"TYPE"] == "DISCOVERY"
		},
		ALL = {
			sel <- rep(TRUE,dim(object)[1])
		}
	)

	object$Normalized_intensity <- object$Raw_intensity
	if(log.it) object$Raw_intensity <- log2(object$Raw_intensity)
	switch(method,
		loess = {
			object$Normalized_intensity[sel,] <- normalize.loess(object$Raw_intensity[sel,], log.it=FALSE)
			object$Normalization_method <- "Cyclic Loess"
		},
		quantiles = {
			require(limma)
			object$Normalized_intensity[sel,] <- normalizeQuantiles(object$Raw_intensity[sel,])
                        object$Normalization_method <- "Quantiles"
		}
	)
	object$Raw_intensity <- NULL
	if(log.it) object$Log_transformed <- TRUE
	return(object)
}

## normalize.loess()
# modified from normalize.loess from affy package.
# Allow NA in input data.
normalize.loess <- function (mat, subset = sample(1:(dim(mat)[1]), min(c(5000, nrow(mat)))), 
    epsilon = 10^-2, maxit = 1, log.it = TRUE, verbose = TRUE, 
    span = 2/3, family.loess = "symmetric") 
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
    w <- c(0, rep(1, length(subset)), 0)
    while (iter < maxit) {
        iter <- iter + 1
        means <- matrix(0, II, J)
        for (j in 1:(J - 1)) {
            for (k in (j + 1):J) {
                y <- newData[, j] - newData[, k]
                x <- (newData[, j] + newData[, k])/2
		# Select genes that are not set to NA
		sel <- which(!is.na(as.character(y)))
		y <- y[sel]
		x <- x[sel]
		#
                index <- c(order(x)[1], subset, order(-x)[1])
                xx <- x[index]
                yy <- y[index]
                aux <- loess(yy ~ xx, span = span, degree = 1, 
                  weights = w, family = family.loess)
                aux <- predict(aux, data.frame(xx = x))/J
		# Apply normalization to genes not NA.
                means[sel, j] <- means[sel, j] + aux
                means[sel, k] <- means[sel, k] - aux
                if (verbose) 
                  cat("Done with", j, "vs", k, " in iteration ", 
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
        warning(paste("No convergence after", maxit, "iterations.\n"))
    if (log.it) {
        return(2^newData)
    }
    else return(newData)
}

