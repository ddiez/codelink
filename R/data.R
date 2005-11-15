# mergeArray.Codelink()
# returns a codelink data.frame mergin chips attending a function
# like median, mean, etc... default mean.
mergeArray.Codelink <- function(object, group, names=NULL, method="mean", log.it=FALSE) {
	if(!is(object,"Codelink")) stop("A Codelink object needed.")
	if(is.null(object$Normalized_intensity) && is.null(object$Raw_intensity)) stop("Normalized_intensity or Raw_intensity slots needed.")
	if(is.null(names)) stop("Group names needed.")
	if(log.it && object$Log_transformed) stop("Data already log transformed.")
	if(is.null(object$SNR)) SNR <- FALSE else SNR <- TRUE

	l <- levels(as.factor(group))
	if(length(names)!=length(l)) stop("Number of groups and group names disagree.")

	method = match.arg(method,c("mean"))
	dimx <- dim(object)[1]
	dimy <- dim(object)[2]
	if(!is.null(object$Raw_intensity)) data <- object$Raw_intensity else data <- object$Normalized_intensity
	val <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	val.cv <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	if(SNR) snr <- matrix(nrow=dimx, ncol=length(l), dimnames=list(rownames(data),l))
	cat("\n\n")
	switch(method,
		mean = {
			d <- 0
			for(n in l) {
				d <- d+1
				sel <- group==n
				cat("  Merging: ",object$Sample_name[sel],"as", names[d],"\n")
				val[,as.numeric(n)] <- apply(data, 1, function(x) mean(x[sel], na.rm=TRUE))
				val.cv[,as.numeric(n)] <- apply(data, 1, function(x) sd(x[sel], na.rm=TRUE)/mean(x[sel],na.rm=TRUE))
				if(SNR) snr[,as.numeric(n)] <- apply(object$SNR, 1, function(x) mean(x[sel], na.rm=TRUE))
			}
			object$Merge_method <- "Mean"
		}
	)
	cat("\n")
	if(log.it) {
		if(!is.null(object$Raw_intensity)) object$Raw_intensity <- log2(val) else object$Normalized_intensity <- log2(val)
		object$Log_transformed <- TRUE
	} else {
		if(!is.null(object$Raw_intensity)) object$Raw_intensity <- val else object$Normalized_intensity <- val
	}
	object$CV <- val.cv
	if(SNR) object$SNR <- snr
	object$Sample_name <- names
	return(object)
}
## bkgdCorrect.Codelink()
# Correct Spot intensity by background.
bkgdCorrect.Codelink <- function(object,method="half") {
	if(!is(object,"Codelink")) {
		stop("Codelink object needed.")
	}
	method <- match.arg(method, c("none","subtract","half","edwards","normexp"))
	switch(method,
		none = {
			object$Raw_intensity <- object$Spot_mean
			object$BkgdCorrect_method <- "None"
		},
		subtract = {
			object$Raw_intensity <- object$Spot_mean - object$Bkgd_median
			object$BkgdCorrection_method <- "Subtract"
		},
		half = {
			object$Raw_intensity <- pmax(object$Spot_mean - object$Bkgd_median, 0.5)
			object$BkgdCorrection_method <- "Half"
		}
	)
	object$Spot_mean <- NULL
	object$Bkgd_median <- NULL
	return(object)
}

## log2.Codelink()
# apply log2 to Codelink object.
log2.Codelink <- function(object) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(object$Log_transformed) stop("Already log2 transformed.")

	if(!is.null(object$Spot_mean)) object$Spot_mean <- log2(object$Spot_mean)
	if(!is.null(object$Raw_intensity)) object$Raw_intensity <- log2(object$Raw_intensity)
	if(!is.null(object$Normalized_intensity)) object$Normalized_intensity <- log2(object$Normalized_intensity)

	object$Log_transformed <- TRUE
	return(object)
}

## snr.Codelink()
# Calculate Codelink Flags: Signal-to-Noise Ratio = Spot_mean / (Bkgd_median + 1.5 * Bkgd_stdev)
# Spots with Signal-to-Noise Ratio < 1 are flagged as "L" (Not present)
snr.Codelink <- function(object) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(is.null(object$Spot_mean)) stop("Spot_mean missing.")
	if(is.null(object$Bkgd_median)) stop("Bkgd_median missing.")
	if(is.null(object$Bkgd_stdev)) stop("Bkgd_stdev missing.")

	snr <- object$Spot_mean / (object$Bkgd_median + 1.5 * object$Bkgd_stdev)
	object$SNR <- snr
	return(object)
}
# This apply inside read.Codelink()
SNR <- function(Spot_mean, Bkgd_median, Bkgd_stdev) {
	Spot_mean / (Bkgd_median + 1.5 * Bkgd_stdev)
}
## fc.Codelink()
# Select denes based in fold change between two conditions.
fc.Codelink <- function(object, cond1=NULL, cond2=NULL, fc=1.0) {
	if(!is(object, "Codelink")) stop("Codelink object needed.")
	if(cond1 > dim(object)[2] || cond2 > dim(object)[2]) stop("Invalid conditions.")
	if(cond1 == cond2) stop("Fold changes along same conditions are equal to 0.")
	if(is.null(object$Normalized_intensity)) stop("Normalize data first.")
	if(object$Log_transformed) {
		fc.0 <- abs(object$Normalized_intensity[,cond1] - object$Normalized_intensity[,cond2])
	} else {
		stop("Only logged data at this moment supported.")
	}
	return(na2false(fc.0 >= fc))
}

## cutCV.Codelink()
# Calculate cutoff based on C.V.
cutCV.Codelink <- function(object, subset=c(1:dim(object)[2])) {
	if(!is(object, "Codelink")) stop("Codelink object needed.")
	if(object$Merge_method=="NONE") stop("Merged object needed.")
	cut.mean <- mean(apply(data$CV[,subset],2,function(x) median(x,na.rm=TRUE)), na.rm=TRUE)
	cut.sd <- sd(apply(data$CV[,subset],2,function(x) median(x,na.rm=TRUE)), na.rm=TRUE)
	cut <- cut.mean + 3*cut.sd
	return(cut)
}
## selCV.Codelink()
# Select based on CV cutoff.
selCV.Codelink <- function(object,cutoff) {
        cv.mean <- apply(object$CV,1,function(x) mean(x, na.rm=TRUE))
        return(na2false(cv.mean <= cutoff))
}
## na2false()
# Set all NA in a logical vector to FALSE.
na2false <- function(x) {
    x[which(is.na(x))] <- FALSE
    x	
}
