# codCorrect-method.
setGeneric("codCorrect", function(x, method = "half")
	standardGeneric("codCorrect"))
setMethod("codCorrect", "CodelinkSet",
function(x, method = "half")
{
	#if(getInfo(x, "background") != "NONE")
		#stop("data already corrected.")
	
	method <- match.arg(method, c("none"< "subtract", "half", "normexp"))

	# take data.
	newInt <- assayDataElement(x, "intensity")
	bkgd <- assayDataElement(x, "background")

	# do correction stuff...
	switch(method,
		none = ,
		subtract = newInt <- newInt - bkgd,
		half = { newInt <- newInt - bkgd; newInt[newInt < 0.5] <- 0.5 },
		normexp = warning("FIXME: not implemented.")
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
	experimentData(x)@preprocessing[[what]]
})

# codNormalize-method.
setGeneric("codNormalize", function(x, method = "quantile", log.it = TRUE)
	standardGeneric("codNormalize"))
setMethod("codNormalize", "CodelinkSet",
function(x, method = "quantile", log.it = TRUE)
{
	#if(getInfo(x, "normalization") != "NONE")
		#stop("data already normalized.")

	method <- match.arg(method, c("median", "quantile", "loess"))
	
	# take data.
	newInt <- assayDataElement(x, "intensity")
	
	# do correction stuff...
	switch(method,
		median = warning("FIXME: not implemented."),
		quantile = warning("FIXME: not implemented."),
		loess = warning("FIXME: not implemented.")
	)

	if(log.it)
		newInt <- log2(newInt)

	# reassign data.
	assayDataElement(x, "intensity") <- newInt
	
	# info.
	experimentData(x)@preprocessing[["normalization"]] <- method
	experimentData(x)@preprocessing[["log"]] <- log.it

	x
})

setGeneric("codImage", function(x, array = 1)
	standardGeneric("codImage"))
setMethod("codImage", "CodelinkSet",
function(x, array = 1)
{
	warning("FIXME: not implemented.")
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
	as.character(pData(featureData(foo2))[, "probeName"])
})

# probeTypes method to get feature types.
setGeneric("probeTypes", function(object) standardGeneric("probeTypes"))
setMethod("probeTypes", "CodelinkSet",
function(object)
{
	as.character(pData(featureData(foo2))[, "probeType"])
})
