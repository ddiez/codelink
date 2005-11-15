## plotMA.Codelink()
# MA plot of gene intensities.
plotMA.Codelink <- function(object, array1=1, array2=2, cutoff=NULL, label="Probe_type", high.list=NULL, high.col="blue", high.pch="*", legend.x="bottomright", title=NULL, xlim=NULL, ylim=NULL) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
#	if(!is.null(high.list) && (!is(high.list,"logical") || !is(high.list,"vector"))) stop("logical vector needed.")
#	if(!is.null(high.list) && length(high.list) != dim(object)[1]) stop("high.list and number of genes differ.")
	label <- match.arg(label,c("Probe_type","SNR","none"))

	if(!is.null(object$Spot_mean)) {
                        val1 <- object$Spot_mean[,array1]
                        val2 <- object$Spot_mean[,array2]
			what <- "Spot_mean"
	}
        if(!is.null(object$Raw_intensity)) {
                        val1 <- object$Raw_intensity[,array1]
                        val2 <- object$Raw_intensity[,array2]
			what <- "Raw_intensity"
	}
if(!is.null(object$Normalized_intensity)) {    
                        val1 <- object$Normalized_intensity[,array1]
                        val2 <- object$Normalized_intensity[,array2]
			what <- "Normalized_intensity"
	}
	if(!object$Log_transformed) {
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
                Probe_type = {
                        negative <- object$Probe_type=="NEGATIVE"
                        positive <- object$Probe_type=="POSITIVE"
                        discovery <- object$Probe_type=="DISCOVERY"
                        fiducial <- object$Probe_type=="FIDUCIAL"
                        other <- object$Probe_type=="OTHER"
                        plot(A[discovery], M[discovery], xlim=xlim, ylim=ylim, xlab="A", ylab="M", pch=".")
                        points(A[negative], M[negative], col="red", pch=20)
                        points(A[positive], M[positive], col="blue", pch=20)
                        points(A[fiducial], M[fiducial], col="yellow", pch=20)
                        points(A[other], M[other], col="green", pch=20)
			legend.text <- c("DISCOVERY", "NEGATIVE", "POSITIVE", "FIDUCIAL", "OTHER")
			legend.fill <- c("black","red","blue","yellow","green")
		},
		SNR = {
                        foo.mean <- apply(object$SNR[,c(array1, array2)],1,mean)
                        foo.sel <- foo.mean >= 1
                        plot(A[foo.sel], M[foo.sel], xlim=xlim, ylim=ylim, col="black", xlab="A", ylab="M", pch=".")
                        points(A[!foo.sel], M[!foo.sel], col="orange", pch=".")
			legend.text <- c("SNR >= 1","SNR < 1")
			legend.fill <- c("black","orange")
		},
		none = {
			plot(A, M, xlab="A", ylab="M", pch=".");
			legend.text <- "ALL PROBES"
			legend.fill <- "black"
		}
	)

	# Highlighted genes.
        if(!is.null(high.list)) {
		#names <- object$Probe_name[high.list]
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
        names <- object$Sample_name[c(array1, array2)]
	#if(is.null(title)) title(paste(object$Experiment_name,"\n(",names[2],"-", names[1], ") MA Plot of ", what, sep=""))
	if(is.null(title)) title(paste("(",names[2],"-", names[1], ") MA Plot of ", what, sep=""))
	else title(title)
	legend(x=legend.x, legend=legend.text, fill=legend.fill, inset=0.05)
}

## plotDensities.Codelink()
# Densities plot of gene intensities.
plotDensities.Codelink <- function(object, subset=c(1:dim(object)[2]), title=NULL, legend.cex=1) {
        if(!is(object,"Codelink")) stop("Codelink object needed.")
        if(!is.null(object$Spot_mean)) {
		val <- object$Spot_mean
		what <- "Spot_mean"
	}
        if(!is.null(object$Raw_intensity)) {
		val <- object$Raw_intensity
		what <- "Raw_intensity"
	}
        if(!is.null(object$Normalized_intensity)) {
		val <- object$Normalized_intensity
		what <- "Normalized_intensity"
	}
	if(!object$Log_transformed) val <- log2(val)
        
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
	legend(x="topright", legend=object$Sample_name[subset], cex=legend.cex, fill=colors, inset=0.05)
}


## plotCorrelation.Codelink()
# Scatter plot of intensities: One array in x compare to other in y axis.
plotCorrelation.Codelink <- function(object, x=1, y=2, cutoff=FALSE, label="Probe_type", title=NULL, xlim=NULL, ylim=NULL) {
        if(!is(object,"Codelink")) stop("Codelink object needed.")       
        if(!is.null(object$Spot_mean)) {
                        xval <- object$Spot_mean[,x]
                        yval <- object$Spot_mean[,y]
			what <- "Spot_mean"
        }       
        if(!is.null(object$Raw_intensity)) {
                        xval <- object$Raw_intensity[,x]
                        yval <- object$Raw_intensity[,y]
			what <- "Raw_intensity"
        } 
        if(!is.null(object$Normalized_intensity)) {    
                        xval <- object$Normalized_intensity[,x]
                        yval <- object$Normalized_intensity[,y]
			what <- "Normalized_intensity"
        }
	if(!object$Log_transformed) {
		xval <- log2(xval)
		yval <- log2(yval)
	}
        names <- object$Sample_name[c(x,y)]
	 # Range computation.
        if(is.null(xlim)) xlim <- range(xval, na.rm=TRUE)
        if(is.null(ylim)) ylim <- range(yval, na.rm=TRUE)
        
	switch(label,
                Probe_type={
                        negative <- object$Probe_type=="NEGATIVE"
                        positive <- object$Probe_type=="POSITIVE"
                        discovery <- object$Probe_type=="DISCOVERY"
                        fiducial <- object$Probe_type=="FIDUCIAL"
                        other <- object$Probe_type=="OTHER"
                        plot(xval[discovery], yval[discovery], xlim=xlim, ylim=ylim, xlab=names[1], ylab=names[2], pch=".")
                        points(xval[negative],yval[negative],col="red",pch=20)
                        points(xval[positive],yval[positive],col="blue",pch=20)
                        points(xval[fiducial],yval[fiducial],col="yellow",pch=20)
                        points(xval[other],yval[other],col="green",pch=20)
			legend.text <- c("DISCOVERY","NEGATIVE","POSITIVE","FIDUCIAL","OTHER")
			legend.fill <- c("black","red","blue","yellow","green")
                },
		SNR = {
                        foo.mean <- apply(object$SNR[,c(x,y)],1,mean)
                        foo.sel <- foo.mean >= 1
                        plot(xval[foo.sel], yval[foo.sel], xlim=xlim, ylim=ylim, col="black", xlab=names[1], ylab=names[2], pch=".")
                        points(xval[!foo.sel],yval[!foo.sel],col="orange",pch=".")
			legend.text <- c("SNR >= 1","SNR < 1")
			legend.fill <- c("black","orange")
                },
                none = {
			plot(xval, yval, pch=".", xlab=names[1], ylab=names[2])
			legend.text <- "ALL PROBES"
			legend.fill <- "black"
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
	legend(x="topleft", legend=legend.text, fill=legend.fill, inset=0.05)
}

## plotCV.Codelink()
# density plot of C.V. of merged objects.
plotCV.Codelink <- function(object, subset=c(1:dim(object)[2]), cutoff=NULL, title=NULL, legend.cex=1) {
	if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(object$Merge_method=="NONE") stop("Merged object needed.")
	colors <- rainbow(length(subset))
	y.max <- c()
	x.max <- c()
	for(n in subset) {
        	y.max[n] <- max(density(object$CV[,n],na.rm=TRUE)$y)
	        x.max[n] <- max(density(object$CV[,n],na.rm=TRUE)$x)
	}
	y.pos <- order(y.max,decreasing=TRUE,na.last=NA)
	for(n in y.pos) {
        	k <- which(y.pos==n)
	        if(n==y.pos[1]) plot(density(object$CV[,n],na.rm=TRUE),col=colors[k], main="")
        	else lines(density(object$CV[,n],na.rm=TRUE),col=colors[k])
	}
	if(is.null(title)) title("CV of Merged Normalized_intensity")
	else title(title)
	legend(x="topright", legend=object$Sample_name[subset], cex=legend.cex, fill=colors, inset=0.05)
	if(!is.null(cutoff)) abline(v=cutoff, lty="dotted")
}

# Generate red-black-green colorscale
#redgreen <- colorRampPalette(c('red','black','green'))
# Generate green-black-red colorscale
#greenred <- colorRampPalette(c('green','black','red'))
# Generate blue white red colorscale
#bluered <- colorRampPalette(c('blue','white','red'))
