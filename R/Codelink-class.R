# Codelink-class definition.
# Based on the definitions found in limma v1.8.20 for RGList objects.
setClass("Codelink", representation("list"))

# Addapted? from limma:
#  Print leading 5 elements or rows of atomic object.
#  Gordon Smyth
#  May 2003.  Last modified 7 April 2004.
printHead <- function(x) {
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

# show method.
setMethod("show","Codelink",
#  Print and show method large data objects
#  Based on same from: limma
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

# [ method.
setMethod("[", "Codelink",
function(x, i, j, ..., drop=FALSE) {
	if(!missing(i)) {
		x$Smean <- x$Smean[i,]
		x$Bmedian <- x$Bmedian[i,]
		x$Bstdev <- x$Bstdev[i,]
		x$snr <- x$snr[i,]
		x$Ri <- x$Ri[i,]
		x$Ni <- x$Ni[i,]
		x$cv <- x$cv[i,]
		x$flag <- x$flag[i,]
		x$name <- x$name[i]
		x$type <- x$type[i]
	}
	if(!missing(j)) {
		x$Smean <- x$Smean[,j]
		x$Bmedian <- x$Bmedian[,j]
		x$Bstdev <- x$Bstdev[,j]
		x$snr <- x$snr[,j]
		x$Ri <- x$Ri[,j]
		x$Ni <- x$Ni[,j]
		x$cv <- x$cv[,j]
		x$flag <- x$flag[,j]
		x$file <- x$file[j]
		x$sample <- x$sample[j]
        }
        return(x)
})

# S3 methods.
# dim.
dim.Codelink <- function(x) {
	if(is.null(x$Ni) && is.null(x$Ri) && is.null(x$Smean)) return(c(0,0))
	if(!is.null(x$Smean)) return(dim(x$Smean))
	if(!is.null(x$Ri)) return(dim(x$Ri))
	if(!is.null(x$Ni)) return(dim(x$Ni))
}
# as.matrix.
as.matrix.Codelink <- function(x) {
	if(is.null(x$Ni) && is.null(x$Ri) && is.null(x$Smean)) return(NULL)
	if(!is.null(x$Smean)) return(x$Smean)
	if(!is.null(x$Ri)) return(x$Ri)
	if(!is.null(x$Ni)) return(x$Ni)
}

##
# onLoad stuff for S4 classes in NAMESPACE.
.onLoad <- function(lib, pkg) {
	require(methods, quietly=TRUE) || stop("Package methods unavailable!")
}
