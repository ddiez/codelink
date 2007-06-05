# temporary wrapper.
readCodelink2 <- function(...) {
	tmp <- readCodelink(...)
	Codelink2eSet(tmp)
}

# convert and old Codelink object to an CodelinkRawSet object.
Codelink2eSet <- function(object, annotation=NULL) {
	if(class(object) != "Codelink") stop("Codelink-object needed.")

	pD <- data.frame(sample=unique(object$sample))
	varMet <- data.frame(labelDescription="sample names", row.names="sample")
	pD <- new("AnnotatedDataFrame", data=pD, varMetadata=varMet)

	Rep <- data.frame(probeName = object$name, probeType = object$type,
		logicalRow = object$logical[, "row"],
		logicalCol = object$logical[, "col"])
	feMet <- data.frame(labelDescription = c("probe names", "probe types", "probe row position", "probe column position"),
		row.names = c("probeName", "probeType", "logicalRow", "logicalCol"))
	Rep <- new("AnnotatedDataFrame", data = Rep, varMetadata = feMet)

	if(is.null(annotation))
		tmp <- gessAnnotation(object$product)
		if(!is.null(tmp)) annotation <- tmp
		else annotation <- ""

	int <- NULL
	bkg <- NULL
	if(!is.null(object$Smean)) int <- object$Smean
	if(!is.null(object$Ri)) int <- object$Ri
	if(!is.null(object$Ni)) int <- object$Ni
	if(!is.null(object$Bmedian)) bkg <- object$Bmedian
	else bkg <- int

	#tmp <- new("CodelinkRawSet", phenoData = pD, intensity = int,
	tmp <- new("CodelinkSet", phenoData = pD, intensity = int,
		background = bkg, flag = object$flag, snr = object$snr,
		annotation = annotation, featureData = Rep)
	
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

# gess the annotation package associated.
gessAnnotation <- function(x)
{
	base <- "cod"

	# organism.
	org <- NULL
	if(any(grep("Rat", x))) org <- "r"
	if(any(grep("Mouse", x))) org <- "m"
	if(any(grep("Human", x))) org <- "h"

	# chip.
	chip <- NULL
	if(any(grep("Whole Genome", x))) chip <- "wg"
	if(any(grep("UniSet.*I", x))) chip <- "10k"
	if(any(grep("UniSet.*II", x))) chip <- "20k"

	# return name.
	if(!is.null(org) & !is.null(chip))
		return(paste(org, chip, base, sep=""))
	else
		return(NULL)
}
