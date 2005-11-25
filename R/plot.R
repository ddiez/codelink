## plotMA()
# MA plot of gene intensities.
plotMA <- function(object, array1=1, array2=2, cutoff=NULL, label="type", high.list=NULL, high.col="blue", high.pch="*", legend.x="bottomright", title=NULL, xlim=NULL, ylim=NULL) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
#	if(!is.null(high.list) && (!is(high.list,"logical") || !is(high.list,"vector"))) stop("logical vector needed.")
#	if(!is.null(high.list) && length(high.list) != dim(object)[1]) stop("high.list and number of genes differ.")
	label <- match.arg(label,c("type", "snr", "none"))

	if(!is.null(object$Smean)) {
                        val1 <- object$Smean[,array1]
                        val2 <- object$Smean[,array2]
			what <- "Smean"
	}
        if(!is.null(object$Ri)) {
                        val1 <- object$Ri[,array1]
                        val2 <- object$Ri[,array2]
			what <- "Ri"
	}
	if(!is.null(object$Ni)) {    
                        val1 <- object$Ni[,array1]
                        val2 <- object$Ni[,array2]
			what <- "Ni"
	}
	if(!object$method$log) {
		val1 <- log2(val1)
		val2 <- log2(val2)
	}
	# M, A computation.
	M <- val2 - val1
	A <- (val2 + val1)/2
        # Range computation.
	if(is.null(xlim)) xlim <- range(A, na.rm=TRUE)
	if(is.null(ylim)) ylim <- range(M, na.rm=TRUE)
	# Plotting.
	switch(label,
                type = {
                        negative <- object$type=="NEGATIVE"
                        positive <- object$type=="POSITIVE"
                        discovery <- object$type=="DISCOVERY"
                        fiducial <- object$type=="FIDUCIAL"
                        other <- object$type=="OTHER"
                        plot(A[discovery], M[discovery], xlim=xlim, ylim=ylim, xlab="A", ylab="M", pch=".")
                        points(A[negative], M[negative], col="red", pch=20)
                        points(A[positive], M[positive], col="blue", pch=20)
                        points(A[fiducial], M[fiducial], col="yellow", pch=20)
                        points(A[other], M[other], col="green", pch=20)
			legend.text <- c("DISCOVERY", "NEGATIVE", "POSITIVE", "FIDUCIAL", "OTHER")
			legend.fill <- c("black","red","blue","yellow","green")
		},
		snr = {
                        foo.mean <- apply(object$snr[,c(array1, array2)],1,mean)
                        foo.sel <- foo.mean >= 1
                        plot(A[foo.sel], M[foo.sel], xlim=xlim, ylim=ylim, col="black", xlab="A", ylab="M", pch=".")
                        points(A[!foo.sel], M[!foo.sel], col="orange", pch=".")
			legend.text <- c("SNR >= 1","SNR < 1")
			legend.fill <- c("black","orange")
		},
		none = {
			plot(A, M, xlab="A", ylab="M", pch=".");
			#legend.text <- "ALL PROBES"
			#legend.fill <- "black"
		}
	)

	# Highlighted genes.
        if(!is.null(high.list)) {
		#names <- object$name[high.list]
                #text(A[high.list], M[high.list], names, col="blue", cex=0.75)
		points(A[high.list],M[high.list],col=high.col,pch=high.pch)
        }

	# Lowess plot.
	sel <- which(!is.na(M))
	M <- M[sel]
	A <- A[sel]
	subset=sample(1:length(M),min(c(10000, length(M))))
	o <- order(A[subset])
	A <- A[subset][o]
	o <- which(!duplicated(A))
	lines(approx(lowess(A[o], M[o])), col = "red")
	
	# Misc.
        abline(h=0, col="blue")
	if(!is.null(cutoff)) {
	        abline(h=-cutoff, lty="dotted")
        	abline(h=cutoff, lty="dotted")
	}
        names <- object$sample[c(array1, array2)]
	#if(is.null(title)) title(paste(object$Experiment_name,"\n(",names[2],"-", names[1], ") MA Plot of ", what, sep=""))
	if(is.null(title)) title(paste("(",names[2],"-", names[1], ") MA Plot of ", what, sep=""))
	else title(title)
	if(label != "none") legend(x=legend.x, legend=legend.text, fill=legend.fill, inset=0.05)
}

## plotDensities()
# Densities plot of gene intensities.
plotDensities <- function(object, subset=c(1:dim(object)[2]), title=NULL, legend.cex=1) {
        if(!is(object,"Codelink")) stop("Codelink object needed.")
        if(!is.null(object$Smean)) {
		val <- object$Smean
		what <- "Smean"
	}
        if(!is.null(object$Ri)) {
		val <- object$Ri
		what <- "Ri"
	}
        if(!is.null(object$Ni)) {
		val <- object$Ni
		what <- "Ni"
	}
	if(!object$method$log) val <- log2(val)
        
	colors <- rainbow(length(subset))
	y.max <- c()
        x.max <- c()
        for(n in subset) {
                y.max[n] <- max(density(val[,n], na.rm=TRUE)$y)
                x.max[n] <- max(density(val[,n], na.rm=TRUE)$x)
        }
        y.pos <- order(y.max, decreasing=TRUE, na.last=NA)
        for(n in y.pos) {
                k <- which(y.pos==n)
                if(n==y.pos[1]) plot(density(val[,n], na.rm=TRUE), col=colors[k], main="")
                else lines(density(val[,n], na.rm=TRUE), col=colors[k])
        }
	if(is.null(title)) title(paste("Density Plot of",what))
	else title(title)
	legend(x="topright", legend=object$sample[subset], cex=legend.cex, fill=colors, inset=0.05)
}


## plotCorrelation()
# Scatter plot of intensities: One array in x compare to other in y axis.
plotCorrelation <- function(object, x=1, y=2, cutoff=FALSE, label="type", title=NULL, xlim=NULL, ylim=NULL) {
        if(!is(object,"Codelink")) stop("Codelink object needed.")       
        if(!is.null(object$Smean)) {
                        xval <- object$Smean[,x]
                        yval <- object$Smean[,y]
			what <- "Smean"
        }       
        if(!is.null(object$Ri)) {
                        xval <- object$Ri[,x]
                        yval <- object$Ri[,y]
			what <- "Ri"
        } 
        if(!is.null(object$Ni)) {    
                        xval <- object$Ni[,x]
                        yval <- object$Ni[,y]
			what <- "Ni"
        }
	if(!object$method$log) {
		xval <- log2(xval)
		yval <- log2(yval)
	}
        names <- object$sample[c(x,y)]
	 # Range computation.
        if(is.null(xlim)) xlim <- range(xval, na.rm=TRUE)
        if(is.null(ylim)) ylim <- range(yval, na.rm=TRUE)
        
	switch(label,
                type={
                        negative <- object$type=="NEGATIVE"
                        positive <- object$type=="POSITIVE"
                        discovery <- object$type=="DISCOVERY"
                        fiducial <- object$type=="FIDUCIAL"
                        other <- object$type=="OTHER"
                        plot(xval[discovery], yval[discovery], xlim=xlim, ylim=ylim, xlab=names[1], ylab=names[2], pch=".")
                        points(xval[negative],yval[negative],col="red",pch=20)
                        points(xval[positive],yval[positive],col="blue",pch=20)
                        points(xval[fiducial],yval[fiducial],col="yellow",pch=20)
                        points(xval[other],yval[other],col="green",pch=20)
			legend.text <- c("DISCOVERY","NEGATIVE","POSITIVE","FIDUCIAL","OTHER")
			legend.fill <- c("black","red","blue","yellow","green")
                },
		snr = {
                        foo.mean <- apply(object$snr[, c(x,y)], 1, mean)
                        foo.sel <- foo.mean >= 1
                        plot(xval[foo.sel], yval[foo.sel], xlim=xlim, ylim=ylim, col="black", xlab=names[1], ylab=names[2], pch=".")
                        points(xval[!foo.sel],yval[!foo.sel],col="orange",pch=".")
			legend.text <- c("SNR >= 1","SNR < 1")
			legend.fill <- c("black","orange")
                },
                none = {
			plot(xval, yval, pch=".", xlab=names[1], ylab=names[2])
			#legend.text <- "ALL PROBES"
			#legend.fill <- "black"
                }
        )

	abline(0,1,col="blue")
	if(!is.null(cutoff)) {
		abline(1,1, lty="dotted")
		abline(-1,1, lty="dotted")
	}
	#if(is.null(title)) title(paste(object$Experiment_name,"\n",names[1]," vs ", names[2], " - Scatter Plot of",what))
	if(is.null(title)) title(paste(names[1]," vs ", names[2], " - Scatter Plot of",what))
	else title(title)
	if(label != "none") legend(x="topleft", legend=legend.text, fill=legend.fill, inset=0.05)
}

## plotCV()
# density plot of C.V. of merged objects.
plotCV <- function(object, subset=c(1:dim(object)[2]), cutoff=NULL, title=NULL, legend.cex=1) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(object$method$merge == "NONE") stop("Merged object needed.")
	colors <- rainbow(length(subset))
	y.max <- c()
	x.max <- c()
	for(n in subset) {
        	y.max[n] <- max(density(object$cv[,n],na.rm=TRUE)$y)
	        x.max[n] <- max(density(object$cv[,n],na.rm=TRUE)$x)
	}
	y.pos <- order(y.max,decreasing=TRUE,na.last=NA)
	for(n in y.pos) {
        	k <- which(y.pos==n)
	        if(n==y.pos[1]) plot(density(object$cv[,n],na.rm=TRUE),col=colors[k], main="")
        	else lines(density(object$cv[,n],na.rm=TRUE),col=colors[k])
	}
	if(is.null(title)) title("CV of Merged Ni")
	else title(title)
	legend(x="topright", legend=object$sample[subset], cex=legend.cex, fill=colors, inset=0.05)
	if(!is.null(cutoff)) abline(v=cutoff, lty="dotted")
}

# Generate red-black-green colorscale
#redgreen <- colorRampPalette(c('red','black','green'))
# Generate green-black-red colorscale
#greenred <- colorRampPalette(c('green','black','red'))
# Generate blue white red colorscale
#bluered <- colorRampPalette(c('blue','white','red'))
