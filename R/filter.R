# genefilter function for filterSNR().
filtersnr <- function(p = 1, A = 1, na.rm = TRUE) {
	function(x) {
		if(na.rm)
			x <- x[!is.na(x)]
		sum(x < A)/length(x) >= p
	}
}

# filter probes when all samples have SNR < cutoff [default: 1]
setGeneric("filterSNR", function(x, cutoff = 1) standardGeneric("filterSNR"))
setMethod("filterSNR", "Codelink",
function(x, cutoff = 1)
{
	filterSNR(x$snr, cutoff = cutoff)
})
setMethod("filterSNR", "CodelinkSet",
function(x, cutoff = 1)
{
	filterSNR(getSNR(x), cutoff = cutoff)
})
setMethod("filterSNR", "matrix",
function(x, cutoff = 1)
{
	f <- filtersnr(A = cutoff)
	!genefilter(x, filterfun(f))
})

# codUniverse()
# create a universe of EG terms, applying a non-specific filter.
setGeneric("codUniverse", function(x, chip) standardGeneric("codUniverse"))
setMethod("codUniverse", "Codelink",
function(x, chip)
{
	if(missing(chip))
		stop("chip type missing.")
	# filter not expressed probes.
	filtersnr <- filterSNR(x)
	x <- x[filtersnr,]

	probesname <- x$name
	mediansnr <- apply(x$snr, 1, function(z) median(z, na.rm = TRUE))
    mediansnr <- data.frame(probeName = probesname, medianSNR = mediansnr)

	createUniverse(chip, probesname, mediansnr)
})
setMethod("codUniverse", "CodelinkSet",
function(x, chip)
{
	if(missing(chip))
		chip <- annotation(x)
	if(nchar(chip) < 1)
		stop("chip type missing.")
	
	# filter not expressed probes.
	filtersnr <- filterSNR(x)
	x <- x[filtersnr,]

	probesname <- probeNames(x)
	# in CodelinkSet is currently stored meanSNR, may be good to have both.
	mediansnr <- apply(getSNR(x), 1, function(z) median(z, na.rm = TRUE))
    mediansnr <- data.frame(probeName = probesname, medianSNR = mediansnr)

	createUniverse(chip, probesname, mediansnr)
})

# more or less general function.
createUniverse <- function(chip, probesname, mediansnr) {
	if(!do.call(require, list(chip)))
		stop("package", chip, "not found.")

	# probes with EG.
	envirEG <- get(paste(chip, "ENTREZID", sep = ""))
	eg <- mget(probesname, envir = envirEG)
	
	haveEG <- sapply(eg, function(x) !is.na(x))
	haveEG <- unique(names(which(haveEG)))

	# probes with GO.
	envirGO <- get(paste(chip, "GO", sep = ""))
	go <- mget(probesname, envir = envirGO)
	haveGO <- sapply(go, function(x)
		ifelse(length(x) == 1 && is.na(x), FALSE, TRUE))
	haveGO <- unique(names(which(haveGO)))

	allProbes <- intersect(haveEG, haveGO)
	
	egSubset <- eg[allProbes]
	egSubsetUnique <- unique(unlist(egSubset))
	allProbesUnique <- sapply(egSubsetUnique, function(x) {
		selEG <- egSubset == x
		probesName <- names(egSubset)[selEG]
		if(length(probesName) == 1) probesName
		else {
			tmp <- mediansnr[mediansnr[, "probeName"] %in% probesName, ]
			if(all(is.na(tmp[, 2]))) probesName[1] 
			else {
				sel <- tmp[ ,"medianSNR"] == max(tmp[, "medianSNR"], na.rm = TRUE)
				sel <- na2false(sel)
				as.character(tmp[sel, "probeName"])
			}
		}
	})

	# now, get the ids from universe and check it's ok.
	egUniverse <- unlist(mget(allProbesUnique, envir = envirEG))
	if(any(duplicated(egUniverse)))
		stop("error in gene universe: can't have duplicated ids.")
	else
		egUniverse
}
