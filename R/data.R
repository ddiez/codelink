# mergeArray()
# returns a codelink data.frame mergin chips attending a function
# like median, mean, etc... default mean.
mergeArray <- function(object, class, names=NULL, method="mean",
			log.it=FALSE, merge.snr=TRUE) {
	if(!is(object,"Codelink")) stop("A Codelink object needed.")
	if(is.null(object$Ni) && is.null(object$Ri)) stop("Ni or Ri slots needed.")
	if(is.null(names)) stop("Group names needed.")
	if(log.it && object$method$log) stop("Data already log transformed.")
	if(is.null(object$snr) || !merge.snr) doSNR <- FALSE else doSNR <- TRUE

	l <- levels(as.factor(class))
	if(length(names)!=length(l))
		stop("Number of classes and class names disagree.")

	method = match.arg(method,c("mean"))
	dimx <- dim(object)[1]
	dimy <- dim(object)[2]
	if(!is.null(object$Ri)) data <- object$Ri else data <- object$Ni
	val <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	val.cv <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	if(doSNR) snr <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	cat("\n\n")
	switch(method,
		mean = {
			d <- 0
			for(n in l) {
				d <- d+1
				sel <- class==n
				cat("  Merging: ",object$sample[sel],"as", names[d],"\n")
				val[,as.numeric(n)] <- apply(data, 1, function(x) mean(x[sel], na.rm=TRUE))
				val.cv[,as.numeric(n)] <- apply(data, 1, function(x) sd(x[sel], na.rm=TRUE)/mean(x[sel],na.rm=TRUE))
				if(doSNR) snr[,as.numeric(n)] <- apply(object$snr, 1, function(x) mean(x[sel], na.rm=TRUE))
			}
			object$method$merge <- "mean"
		}
	)
	cat("\n")
	if(log.it) {
		if(!is.null(object$Ri)) object$Ri <- log2(val) else object$Ni <- log2(val)
		object$method$log <- TRUE
	} else {
		if(!is.null(object$Ri)) object$Ri <- val else object$Ni <- val
	}
	object$cv <- val.cv
	if(doSNR) object$snr <- snr
	object$sample <- names
	return(object)
}
## bkgdCorrect()
# Correct Spot intensity by background.
bkgdCorrect <- function(object, method = "half", preserve = FALSE,
                        verbose = FALSE, offset = 0) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	method <- match.arg(method, c("none", "subtract", "half", "normexp"))
	switch(method,
		none = {
			object$Ri <- object$Smean
			object$method$background <- "NONE"
		},
		subtract = {
			object$Ri <- object$Smean - object$Bmedian
			object$method$background <- "subtract"
		},
		half = {
			object$Ri <- pmax(object$Smean - object$Bmedian, 0.5)
			object$method$background <- "half"
		},
		normexp = {
			object$Ri <- object$Smean
			for (j in 1:ncol(object$Smean)) {
			    x <- object$Smean[, j] - object$Bmedian[, j]
			    out <- normexp.fit(x)
			    object$Ri[, j] <- normexp.signal(out$par, x)
			    if (verbose)
					cat(" + Corrected array", j, "\n")
			}
			object$method$background <- "normexp"
		}
	)
	if(!preserve) object$Smean <- NULL
	if(!preserve) object$Bmedian <- NULL
	if(offset) object$Ri <- object$Ri + offset
	return(object)
}

# bgcorrect.normexp()
bgcorrect.normexp <- function(x, y, verbose = FALSE) {
	#z <- matrix(NA, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
	for (j in 1:ncol(x)) {
		xc <- x[, j] - y[, j]
		out <- normexp.fit(xc)
		x[, j] <- normexp.signal(out$par, xc)
		if (verbose)
			cat(" + Corrected array", j, "\n")
	}
	x
}

## logCodelink()
# apply log2 to Codelink object.
logCodelink <- function(object) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(object$method$log) stop("Already log2 transformed.")

	if(!is.null(object$Smean)) object$Smean <- log2(object$Smean)
	if(!is.null(object$Ri)) object$Ri <- log2(object$Ri)
	if(!is.null(object$Ni)) object$Ni <- log2(object$Ni)

	object$method$log <- TRUE
	return(object)
}

## SNR()
# Compute Signal to Noise Ratio.
SNR <- function(Smean, Bmedian, Bstdev) {
	Smean / (Bmedian + 1.5 * Bstdev)
}
## fc2Cond()
# rewrite of fc to be more intuitive:
fc2Cond <- function(object, cond1=NULL, cond2=NULL, fc=1.0, verbose=FALSE) {
	if(!is(object, "Codelink")) stop("Codelink object needed.")
	if(!is(cond1, "numeric") && !is(cond2, "numeric")) {
		#cat("Cond1:", cond1,"\n")
		#cat("Cond2:", cond2,"\n")
		cond1.inx <- which(object$sample==cond1)
		cond2.inx <- which(object$sample==cond2)
		#cat("Index1:", cond1.inx, "\n")
		#cat("Index2:", cond2.inx, "\n")
		if(!any(cond1.inx)) stop("Unknown sample ",cond1)
		if(!any(cond2.inx)) stop("Unknown sample ",cond2)
		if(length(cond1.inx)>1) stop(cond1," matches more than one condition in sample slot.")
		if(length(cond2.inx)>1) stop(cond2," matches more than one condition in sample slot.")
	} else {
		cond1.inx <- cond1
		cond2.inx <- cond2
	}
	if(!object$method$log) {
		if(verbose) cat("Commputing log2...\n")
		object <- logCodelink(object)
	}
	if(verbose) cat("Computing FC for", cond1,"-",cond2,"...\n")
	fc.0 <- abs(object$Ni[,cond1.inx] - object$Ni[, cond2.inx])
	return(na2false(fc.0 >= fc))
}
## fc()
# Select denes based in fold change between two conditions.
#fc <- function(object, cond1=NULL, cond2=NULL, fc=1.0) {
#	if(!is(object, "Codelink")) stop("Codelink object needed.")
#	if(cond1 > dim(object)[2] || cond2 > dim(object)[2]) stop("Invalid conditions.")
#	if(cond1 == cond2) stop("Fold changes along same conditions are equal to 0.")
#	if(is.null(object$Ni)) stop("Normalize data first.")
#	if(object$method$log) {
#		fc.0 <- abs(object$Ni[,cond1] - object$Ni[,cond2])
#	} else {
#		stop("Only logged data at this moment supported.")
#	}
#	return(na2false(fc.0 >= fc))
#}

## cutCV()
# Calculate cutoff based on C.V.
cutCV <- function(object, subset=c(1:dim(object)[2])) {
	if(!is(object, "Codelink")) stop("Codelink object needed.")
	if(object$method$merge == "NONE") stop("Merged object needed.")
	cut.mean <- mean(apply(data$cv[,subset],2,function(x) median(x,na.rm=TRUE)), na.rm=TRUE)
	cut.sd <- sd(apply(data$cv[,subset],2,function(x) median(x,na.rm=TRUE)), na.rm=TRUE)
	cut <- cut.mean + 3*cut.sd
	return(cut)
}
## selCV()
# Select based on CV cutoff.
selCV <- function(object,cutoff) {
        cv.mean <- apply(object$cv,1,function(x) mean(x, na.rm=TRUE))
        return(na2false(cv.mean <= cutoff))
}
## na2false()
# Set all NA in a logical vector to FALSE.
na2false <- function(x) {
    x[which(is.na(x))] <- FALSE
    x	
}
## createWeights()
# Create weights based on Probe_type
# createWeights <- function(object, type.weights=NULL) {
# 	if(is.null(type))
# 		type.weights = c("DISCOVERY"=1,"FIDUCIAL"=0,"POSITIVE"=0,"NEGATIVE"=0, "OTHER"=0)
# 	
# 	tw = type.weights
# 	
# 	w <- array(1, dim(object))
# 	discovery <- object$type=="DISCOVERY"
# 	fiducial <- object$type=="FIDUCIAL"
# 	positive <- object$type=="POSITIVE"
# 	negative <- object$type=="NEGATIVE"
# 	other <- object$type=="OTHER"
# 	if(!is.null(type$DISCOVERY)) w[discovery,] <- tw$DISCOVERY
# 	if(!is.null(type$FIDUCIAL)) w[fiducial,] <- tw$FIDUCIAL
# 	if(!is.null(type$POSITIVE)) w[positive,] <- tw$POSITIVE
# 	if(!is.null(type$NEGATIVE)) w[negative,] <- tw$NEGATIVE
# 	if(!is.null(type$OTHER)) w[other,] <- tw$OTHER
# 	w
# }


setMethod("summary", "Codelink",
function(object, ...)
{
	summaryFlag(object)
})

setGeneric("summaryFlag", function(object) standardGeneric("summaryFlag"))
setMethod("summaryFlag", "Codelink",
function(object)
{
	flags <- c("G", "L", "S", "C", "I", "M", "X")
	res <- matrix(NA, length(flags), ncol(object),
		dimnames=list(flags, object$sample))
	for(flag in flags) {
		res[flag,] <- apply(object$flag, 2, function(z) length(grep(flag, z)))
	}
	res
})
