#
setGeneric("writeCluster", function(x, file, chip, ...)
	standardGeneric("writeCluster"))
setMethod("writeCluster", "Codelink",
function(x, file, chip)
{
	mat <- as.matrix(x)
	samples  <- x$sample
	probes <- x$name

	writecluster(file, mat, samples, probes, chip)
})

setMethod("writeCluster", "CodelinkSet",
function(x, file, chip)
{
	mat <- exprs(x)
	samples <- sampleNames(x)
	probes <- probeNames(x)

	if(missing(chip))
		chip <- annotation(x)
	if(chip == "")
		stop("invalid chip name ''")
	
	writecluster(file, mat, samples, probes, chip)
})

setMethod("writeCluster", "MArrayLM", 
function(x, file, chip, probes)
{
	mat <- x$coef
	samples <- colnames(x$coef)
	# try to get probe names from MArrayLM.
	if(missing(probes))
		probes <- x$genes[, "ID"]
	
	writecluster(file, mat, samples, probes, chip)
})

writecluster <- function(file, mat, samples, probes, chip) {
	if(!do.call("require", list(chip)))
		stop(chip, " package not found")

	symbol <- unlist(mget(probes, envir = get(paste(chip, "SYMBOL", sep = ""))))
	desc <- unlist(mget(probes, envir = get(paste(chip, "GENENAME", sep = ""))))
	info <- paste(probes, symbol, desc, sep = " | ")

	header <- c("ORF", "NAME", samples)

	mat <- cbind(probes, info, mat)
	mat <- rbind(header, mat)
	
	write.table(mat, file = file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}
