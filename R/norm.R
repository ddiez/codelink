# normalize()
# Wrapper function to apply normalization methods to Codelink objects.
normalize <- function(object, method = "quantiles", log.it = TRUE, preserve = FALSE) {
	if(!is(object,"Codelink")) stop("A Codelink object is needed.")
	if(is.null(object$Ri)) stop("Background corrected intensities needed.")
	if(log.it & object$method$log) stop("Intensities already log2.")
	
	method <- match.arg(method, c("loess", "quantiles", "median"))

	object$Ni <- object$Ri
	if(log.it) object$Ni <- log2(object$Ni)
	switch(method,
		loess = {
			object$Ni <- normalize.loess(object$Ni, log.it=FALSE)
			object$method$normalization <- "CyclicLoess"
		},
		quantiles = {
			object$Ni <- normalizeQuantiles(object$Ni)
            object$method$normalization <- "quantiles"
		},
		median = {
			# taken from limma.
			#if (is.null(weights))
            	for (j in 1:dim(object)[2]) object$Ni[, j] <- object$Ni[, 
                	j] - median(object$Ni[, j], na.rm = TRUE)
			#else
			#	for (j in 1:narrays) object$M[, j] <- object$M[, j] 
			#		- weighted.median(object$M[, j], weights[, j], na.rm = TRUE)
            object$method$normalization <- "median"
		}
	)
	if(!preserve) object$Ri <- NULL
	if(log.it) object$method$log <- TRUE
	return(object)
}

# normalize.loess()
# modified from normalize.loess() from affy package.
# allows NA in input data.
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

