# codCorrect-method.
setGeneric("codCorrect", function(x, method = "half")
	standardGeneric("codCorrect"))
setMethod("codCorrect", "CodelinkSet",
function(x, method = "half")
{
	if(getInfo(x, "background") != "NONE")
		warning("data already corrected.")
	
	method <- match.arg(method, c("none"< "subtract", "half", "normexp"))

	# take data.
	int <- assayDataElement(x, "intensity")
	bkgd <- assayDataElement(x, "background")

	# do correction stuff...
	switch(method,
		none = newInt <- int,
		subtract = newInt <- int - bkgd,
		half = newInt <- pmax(int - bkgd, 0.5),
		normexp = newInt <- bgcorrect.normexp(int, bkgd)
	)

	# reassign data.
	#assayDataElement(x, "background") <- NULL
	assayDataElement(x, "intensity") <- newInt
	
	# info.
	experimentData(x)@preprocessing[["background"]] <- method
	
	as(x, "CodelinkSet")
})

# getInfo-method.
setGeneric("getInfo", function(x, what)
	standardGeneric("getInfo"))
setMethod("getInfo", "CodelinkSet",
function(x, what)
{
	if(missing(what))
		experimentData(x)@preprocessing
	else
		experimentData(x)@preprocessing[[what]]
})

# codNormalize-method.
setGeneric("codNormalize", function(x, method = "quantile", log.it = TRUE)
	standardGeneric("codNormalize"))
setMethod("codNormalize", "CodelinkSet",
function(x, method = "quantile", log.it = TRUE)
{
	method <- match.arg(method, c("median", "quantile", "loess"))

	if(getInfo(x, "normalization") != "NONE")
		warning("data already normalized.")
	
	# take data.
	newInt <- assayDataElement(x, "intensity")

	# log.it?
	if(log.it & !getInfo(x, "log"))
		newInt <- log2(newInt)
	else
		if(log.it & getInfo(x, "log"))
			warning("already log2 data, skipping...")

	# do normalization stuff.
	switch(method,
		median = newInt <- normalize.median(newInt),
		quantile = newInt <- normalizeQuantiles(newInt),
		loess = newInt <- normalize.loess(newInt, log.it = FALSE,
			verbose = FALSE)
	)

	# reassign data.
	assayDataElement(x, "intensity") <- newInt
	
	# info.
	experimentData(x)@preprocessing[["normalization"]] <- method
	experimentData(x)@preprocessing[["log"]] <- log.it

	x
})

# codPreprocess-method.
setGeneric("codPreprocess", function(x, bg.method = "half", 
	norm.method = "quantile", log.it = TRUE) standardGeneric("codPreprocess"))
setMethod("codPreprocess", "CodelinkSet",
function(x, bg.method = "half", norm.method = "quantile", log.it = TRUE)
{
	x <- codCorrect(x, method = bg.method)
	x <- codNormalize(x, method = norm.method, log.it = log.it)
	x
})

# exprs method.
setMethod("exprs", "CodelinkSet",
function(object)
{
	getInt(object)
})

# getInt method.
setGeneric("getInt", function(object) standardGeneric("getInt"))
setMethod("getInt", "CodelinkSet",
function(object)
{
	assayDataElement(object, "intensity")
})

# getBkg method.
setGeneric("getBkg", function(object) standardGeneric("getBkg"))
setMethod("getBkg", "CodelinkSet",
function(object)
{
	assayDataElement(object, "background")
})

# getSNR method.
setGeneric("getSNR", function(object) standardGeneric("getSNR"))
setMethod("getSNR", "CodelinkSet",
function(object)
{
	assayDataElement(object, "snr")
})

# getFlag method.
setGeneric("getFlag", function(object) standardGeneric("getFlag"))
setMethod("getFlag", "CodelinkSet",
function(object)
{
	assayDataElement(object, "flag")
})

# probeNames-method.
setGeneric("probeNames", function(object) standardGeneric("probeNames"))
setMethod("probeNames", "CodelinkSet",
function(object)
{
	as.character(pData(featureData(object))[, "probeName"])
})

# probeTypes method to get feature types.
setGeneric("probeTypes", function(object) standardGeneric("probeTypes"))
setMethod("probeTypes", "CodelinkSet",
function(object)
{
	as.character(pData(featureData(object))[, "probeType"])
})

# probeTypes method to get feature types.
setGeneric("meanSNR", function(object) standardGeneric("meanSNR"))
setMethod("meanSNR", "CodelinkSet",
function(object)
{
	pData(featureData(object))[, "meanSNR"]
})

# experimental codPlot-method.
# this intends to be a general interface for all plotting utilities.
setGeneric("codPlot", function(x, array, what = "ma", ...)
	standardGeneric("codPlot"))
setMethod("codPlot", "CodelinkSet",
function(x, array, what = "ma", ...)
{
	switch(what,
		ma = { if(missing(array)) array = 1; codPlotMA(x, array, ...) },
		density = {
			if(missing(array)) array = NULL
			codPlotDensity(x, array, ...) },
		scatter = {
			if(missing(array)) array = 1;
			codPlotScatter(x, array, ...) },
		image = { if(missing(array)) array = 1; codPlotImage(x, array, ...) }
	)
})
# codPlotMA, CodelinkSet-method.
setGeneric("codPlotMA", function(x, array = 1, ...)
	standardGeneric("codPlotMA"))
setMethod("codPlotMA", "CodelinkSet",
function(x, array = 1, array2, label = "type", cutoff = c(-1, 1),
	snr.cutoff = 1,	legend.x, pch = ".", ...)
{
	label <- match.arg(label, c("type", "snr", "none"))
	type <- probeTypes(x)
	snr <- meanSNR(x)

	# compute MA values.
	X <- getInt(x)
	islog <- getInfo(x, "log")

	X1 <- X[, array]
	if(missing(array2)) {
		# compute mean in non-log scale and restore it if needed.
		X2 <- rowMeans(if(!islog) X else 2**X, na.rm = TRUE)
		if(islog) X2 <- log2(X2)
	} else
		X2 <- X[, array2]

	if(!islog) {
		X1 <- log2(X1)
		X2 <- log2(X2)
	}
	
	M <- X2 - X1
	A <- (X1 + X2) / 2

	plotma(A, M, label = label,	cutoff = cutoff, snr.cutoff = snr.cutoff,
		legend.x = legend.x, pch = pch, type = type, snr = snr, ...)
})
# codPlotMA, MArrayLM-method.
setMethod("codPlotMA", "MArrayLM",
function(x, array = 1, label = "type", cutoff = c(-1, 1),
	snr.cutofof = 1, legend.x, pch = ".", Amean, ...)
{
	label <- match.arg(label, c("type", "snr", "none"))

	M <- x$coef[, array]
	if(is.null(x$Amean))
		if(missing(Amean)) stop("no Amean information.")
		else A <- Amean
	else
		A <- x$Amean
	
	plotma(A, M, label = label,	cutoff = cutoff, snr.cutoff = snr.cutoff,
		legend.x = legend.x, pch = pch, type = type, snr = snr, ...)
})

# codPlotDensity-method.
setGeneric("codPlotDensity", function(x, array, ...)
	standardGeneric("codPlotDensity"))
setMethod("codPlotDensity", "CodelinkSet",
function(x, array = NULL, ...)
{
	X <- getInt(x)

	if(is.null(array))
		array <- 1:dim(X)[2]

	narray <- length(array)
	col <- rainbow(narray)

	d <- list()
	xlim <- c()
	ylim <- c()
	for(n in 1:narray) {
		d[[n]] <- density(X[, n], na.rm = TRUE)
		xlim <- c(xlim, d[[n]]$x)
		ylim <- c(ylim, d[[n]]$y)
	}

	xlim <- range(xlim, na.rm = TRUE)
	ylim <- range(ylim, na.rm = TRUE)

	plot(0, col = "white", xlim = xlim, ylim = ylim, xlab = "Intensity",
		ylab = "Density")
	
	for(n in 1:narray)
		lines(d[[n]], col = col[n], lwd = 2)
	
	legend("topright", legend = sampleNames(x), fill = col, bty = "n")
})
# codPlotImage-method.
setGeneric("codPlotImage", function(x, array = 1, ...)
	standardGeneric("codPlotImage"))
setMethod("codPlotImage", "CodelinkSet",
function(x, array = 1, ...)
{
	warning("FIXME: codPlotImage()  not implemented.")
})
