# temporary wrapper.
#readCodelink2 <- function(..., phenodata = NULL, featuredata = NULL) {
#	#colnames = list(Signal = "Spot_mean", Background = "Bkgd_median")) {
#	tmp <- readCodelink(...)
#	#Codelink2eSet(tmp)
#	c2e(tmp, phenodata = phenodata, featuredata = featuredata)
#}

# readCodelinkSet <- function(targets, filename, columns = list(Signal = "Spot_mean", Background = "Bkgd_median"), phenoData, ...) {
# 	if(missing(filename) && missing(targets)) stop("argument 'targets' or 'filename' must be specified.")
#     if (missing(targets)) {
#         filename <- as.character(filename)
# 		if(missing(phenoData)) pdata <- NULL
# 		else pdata <- phenoData
#     } else {
#         pdata <- read.AnnotatedDataFrame(targets)
# 		if(! any(grep("FileName", varLabels(pdata)))) {
#             filename <- sampleNames(pdata)
#         } else {
#             filename <- pData(pdata)[, "FileName"]
#         }
#     }
#     if (is.null(filename)) stop("invalid filenames.")
# 	
# 	tmp <- readCodelink(files = filename, ...)
# 	
# 	c2e(tmp, phenodata = pdata, featuredata = NULL)
# }

readCodelinkSet <- function(filename, columns = list(Signal = "Spot_mean", Background = "Bkgd_median"), phenoData=NULL, ...) {
	if(missing(filename)) stop("argument 'filename' must be specified.")
	#tmp <- readCodelink(files = filename, ...)
	tmp = .readCodelinkRaw(files = filename, ...)
	Codelink2CodelinkSet(tmp, phenodata = phenoData)
}


Codelink2CodelinkSet <- function (object, annotation = NULL, phenodata = NULL, featuredata = NULL, intensity = "Smean") 
{
    if (class(object) != "Codelink") 
        stop("Codelink-object needed.")

    switch(intensity, 
    	"Smean" = int <- object$Smean,
		"Ri" = int <- object$Ri,
    	"Ni" = int <- object$Ni,
    )    
    bkg <- object$Bmedian

    if (is.null(phenodata)) {
	    phenodata <- data.frame(sample = object$sample)
    	phenodata.varMet <- data.frame(labelDescription = "sample names", row.names = "sample")
    	phenodata <- new("AnnotatedDataFrame", data = phenodata, varMetadata = phenodata.varMet)
    }
    
    if (is.null(featuredata)) {
    	featuredata <- data.frame(probeName = object$name, probeType = object$type, logicalRow = object$logical[, "row"], logicalCol = object$logical[, "col"], meanSNR = rowMeans(object$snr, na.rm = TRUE))
    	featuredata.feMet <- data.frame(labelDescription = c("probe names", "probe types", "probe row position", "probe column position", "mean snr"), row.names = c("probeName", "probeType", "logicalRow",      "logicalCol", "meanSNR"))
		featuredata <- new("AnnotatedDataFrame", data = featuredata, varMetadata = featuredata.feMet)

    }
    
    if (is.null(annotation))
    	chip <- annotation(object)
    
	tmp <- new("CodelinkSet", exprs = int, background = bkg, 
        flag = object$flag, weight=object$weight, snr = object$snr, annotation = chip)
    phenoData(tmp) <- phenodata
    featureData(tmp) <- featuredata
	experimentData(tmp)@preprocessing <- object$method
	experimentData(tmp)@other <- list("product" = object$product)
    tmp
}

# convert and old Codelink object to an CodelinkRawSet object.
Codelink2eSet <- function(object, annotation = NULL) {
	if(class(object) != "Codelink") stop("Codelink-object needed.")

	pD <- data.frame(sample=unique(object$sample))
	varMet <- data.frame(labelDescription="sample names", row.names="sample")
	pD <- new("AnnotatedDataFrame", data=pD, varMetadata=varMet)

	Rep <- data.frame(probeName = object$name, probeType = object$type,
		logicalRow = object$logical[, "row"],
		logicalCol = object$logical[, "col"], meanSNR = rowMeans(object$snr, na.rm = TRUE))
	feMet <- data.frame(labelDescription = c("probe names", "probe types", "probe row position", "probe column position", "mean snr"),
		row.names = c("probeName", "probeType", "logicalRow", "logicalCol", "meanSNR"))
	Rep <- new("AnnotatedDataFrame", data = Rep, varMetadata = feMet)

	if(is.null(annotation))
		chip <- annotation(object)

	int <- NULL
	bkg <- NULL
	if(!is.null(object$Smean)) int <- object$Smean
	if(!is.null(object$Ri)) int <- object$Ri
	if(!is.null(object$Ni)) int <- object$Ni
	if(!is.null(object$Bmedian)) bkg <- object$Bmedian
	else bkg <- int

	#tmp <- new("CodelinkRawSet", phenoData = pD, intensity = int,
	#tmp <- new("CodelinkSet", phenoData = pD, intensity = int,
	tmp <- new("CodelinkSet", phenoData = pD,  exprs = int,
		background = bkg, flag = object$flag, snr = object$snr,
		annotation = chip, featureData = Rep)
	
	featureNames(tmp) <- object$id
	sampleNames(tmp) <- object$sample
	
	experimentData(tmp)@preprocessing <- object$method
	experimentData(tmp)@other <- list("product" = object$product)

	tmp
}

# TODO: this function makes the opposite conversion.
eSet2Codelink <- function(object, annotation=NULL)
{
	tmp <- new("Codelink")

	warning("FIXME: not implemented.")

	tmp
}

# guess the annotation package associated.
setMethod("annotation", "Codelink",
function(object)
{
	base <- "cod"
	tmp <- object$product

	# organism.
	org <- ""
	if(any(grep("Rat", tmp))) org <- "r"
	if(any(grep("Mouse", tmp))) org <- "m"
	if(any(grep("Human", tmp))) org <- "h"

	# chip.
	chip <- ""
	if(any(grep("Whole Genome", tmp))) chip <- "wg"
	if(any(grep("UniSet.*I", tmp))) chip <- "10k"
	if(any(grep("UniSet.*II", tmp))) chip <- "20k"

	# return name.
	if(org != "" && chip != "")
		ann <- paste(org, chip, base, ".db", sep="")
	else ann <- ""
	
	ann
})
