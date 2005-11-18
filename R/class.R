#### To create Codelink Classes.
## I base it on the definitions found in limma v1.8.20
setClass("Codelink",representation("list"))
#setClass("LargeDataObject")
#setIs("Codelink","LargeDataObject")
printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  Gordon Smyth
#  May 2003.  Last modified 7 April 2004.
{
        if(is.atomic(x)) {
                d <- dim(x)
                if(length(d)<2) which <- "OneD"
                if(length(d)==2) which <- "TwoD"
                if(length(d)>2) which <- "Array"
        } else {
                if(inherits(x,"data.frame")) {
                        d <- dim(x)
                        which <- "TwoD"
                } else
                        which <- "Recursive"
        }
        switch(which,
        OneD={
                n <- length(x)
                if(n > 20) {
                        print(x[1:5])
                        cat(n-5,"more elements ...\n")
                } else
                        print(x)
        },
        TwoD={
                n <- d[1]
                if(n > 10) {
                        print(x[1:5,])
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Array={
                n <- d[1]
                if(n > 10) {
                        dn <- dimnames(x)
                        dim(x) <- c(d[1],prod(d[-1]))
                        x <- x[1:5,]
                        dim(x) <- c(5,d[-1])
                        if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
                        dimnames(x) <- dn
                        print(x)
                        cat(n-5,"more rows ...\n")
                } else
                        print(x)
        },
        Recursive=print(x)
        )
}

#setMethod("show","LargeDataObject",
setMethod("show","Codelink",
#  Print and show method large data objects
#  Based on same from:
#  Gordon Smyth
#  May 2003
function(object) {
        cat("An object of class \"",class(object),"\"\n",sep="")
        for (what in names(object)) {
                x <- object[[what]]
                cat("$",what,"\n",sep="")
                printHead(x)
                cat("\n")
        }
        for (what in setdiff(slotNames(object),".Data")) {
                x <- slot(object,what)
                if(length(x) > 0) {
                        cat("@",what,"\n",sep="")
                        printHead(x)
                        cat("\n")
                }
        }
})

dim.Codelink <- function(x) {
	if(is.null(x$Normalized_intensity) && is.null(x$Raw_intensity) && is.null(x$Spot_mean)) return(c(0,0))
	if(!is.null(x$Spot_mean)) return(dim(x$Spot_mean))
	if(!is.null(x$Raw_intensity)) return(dim(x$Raw_intensity))
	if(!is.null(x$Normalized_intensity)) return(dim(x$Normalized_intensity))
}

as.matrix.Codelink <- function(x) {
	if(is.null(x$Normalized_intensity) && is.null(x$Raw_intensity) && is.null(x$Spot_mean)) return(NULL)
	if(!is.null(x$Spot_mean)) return(x$Spot_mean)
	if(!is.null(x$Raw_intensity)) return(x$Raw_intensity)
	if(!is.null(x$Normalized_intensity)) return(x$Normalized_intensity)
}

setMethod("[","Codelink",
# Subsetting method.
function(x,i,j,drop=FALSE) {
        if(!missing(i)) {
                x$Spot_mean <- x$Spot_mean[i,]
                x$Bkgd_median <- x$Bkgd_median[i,]
                x$Bkgd_stdev <- x$Bkgd_stdev[i,]
                x$SNR <- x$SNR[i,]
                x$Raw_intensity <- x$Raw_intensity[i,]
                x$Normalized_intensity <- x$Normalized_intensity[i,]
                x$CV <- x$CV[i,]
                x$Quality_flag <- x$Quality_flag[i,]
                x$Probe_name <- x$Probe_name[i]
                x$Probe_type <- x$Probe_type[i]
        }
        if(!missing(j)) {
                x$Spot_mean <- x$Spot_mean[,j]
                x$Bkgd_median <- x$Bkgd_median[,j]
                x$Bkgd_stdev <- x$Bkgd_stdev[,j]
                x$SNR <- x$SNR[,j]
                x$Raw_intensity <- x$Raw_intensity[,j]
                x$Normalized_intensity <- x$Normalized_intensity[,j]
                x$CV <- x$CV[,j]
                x$Quality_flag <- x$Quality_flag[,j]
		x$File_name <- x$File_name[j]
		x$Sample_name <- x$Sample_name[j]
        }
        return(x)
})

