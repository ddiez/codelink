TYPE_COLOR = c(
	DISCOVERY = "black",
	POSITIVE = "gray",
	NEGATIVE = "gray",
	OTHER = "gray",
	FIDUCIAL = "gray"
)

TYPE_BG = c(
	DISCOVERY = "black",
	POSITIVE = "blue",
	NEGATIVE = "red",
	OTHER = "green",
	FIDUCIAL = "yellow"
)

TYPE_PCH = list(
	DISCOVERY = ".",
	POSITIVE = 21,
	NEGATIVE = 21,
	OTHER = 21,
	FIDUCIAL = 21
)
# plotMA()
# MA plot of gene intensities.
plotMA <- function(object, array1 = 1, array2 = NULL, cutoff = c(-1, 1),
	label = NULL, type = NULL, high.list = NULL, high.col = "gray",
	high.pch = 21, high.bg = "orange", snr = NULL, snr.cutoff = 1, 
	legend.x = NULL, pch = ".", subset = NULL, title = NULL, xlim = NULL,
	ylim = NULL)
{

	if(!is.null(type) && length(type) != dim(object)[1])
		warning("'type' has different length than object.")
	
	# get values.
	islog <- object$method$log
	switch(class(object),
		Codelink = {
			if(!is.null(object$Smean)) {
					val1 <- object$Smean[,array1]
					if(is.null(array2)) {
						val2 <- rowMeans(if(!islog) object$Smean
							else 2**object$Smean, na.rm = TRUE)
						if(islog) val2 <- log2(val2)
					} else
						val2 <- object$Smean[,array2]
			}
			if(!is.null(object$Ri)) {
					val1 <- object$Ri[,array1]
					if(is.null(array2)) {
						val2 <- rowMeans(if(!islog) object$Ri else
							2**object$Ri, na.rm = TRUE)
						if(islog) val2 <- log2(val2)
					} else
						val2 <- object$Ri[,array2]
			}
			if(!is.null(object$Ni)) {    
					val1 <- object$Ni[,array1]
					if(is.null(array2)) {
						val2 <- rowMeans(if(!islog) object$Ni else
							2**object$Ni, na.rm = TRUE)
						if(islog) val2 <- log2(val2)
					} else
						val2 <- object$Ni[,array2]
			}
			if(!islog) {
				val1 <- log2(val1)
				val2 <- log2(val2)
			}
			M <- val2 - val1
			A <- (val2 + val1)/2
			type <- object$type
			if(is.null(snr))
				snr <- rowMeans(object$snr, na.rm = TRUE)
			if(is.null(label) && !is.null(type))
				label = "type"
			if(is.null(array2))
				names <- paste("Mean Array vs.", object$sample[array1])
			else
				names <- paste(object$sample[array2], "vs.",
					object$sample[array1])
		},
		MArrayLM = {
			M <- as.matrix(object$coefficients)[, array1]
			A <- object$Amean
			# SNRmean can be added manually to the MArrayLM objects.
			if(is.null(snr))
				snr <- object$SNRmean
			names <- colnames(object$contrasts)[array1]
		},
		stop("invalid parameter object")
	)

	# check label.
	label <- match.arg(label,c("type", "snr", "none"))
	if(is.null(type) & label == "type") {
		warning("missing 'type' information for labelling spots.\n")
		label <- "none"
	}
	if(is.null(snr) & label == "snr") {
		warning("missing 'snr' information for labelling spots.\n")
		label <- "none"
	}

	# type based subset.
	if(!is.null(subset)) {
		if(is.null(type)) stop("need type information for subset.\n")
		else {
			subset <- match.arg(subset, levels(as.factor(type)), several.ok=TRUE)
			subset.sel <- type %in% as.character(subset)
			M <- M[subset.sel]
			A <- A[subset.sel]
			type <- type[subset.sel]
			snr <- snr[subset.sel]
			if(!is.null(high.list)) high.list <- high.list[subset.sel]
		}
	}

    # plot range.
	if(is.null(xlim)) xlim <- range(A, na.rm = TRUE)
	if(is.null(ylim)) ylim <- range(M, na.rm = TRUE)

	# basic plot.
	plot(0, col = "white", xlim = xlim, ylim = ylim, xlab="A", ylab="M")
	abline(h = 0, col = "steelblue", lwd = 2)
	if(!is.null(cutoff)) {
		abline(h=-cutoff, col = "gray", lty = "dotted")
		abline(h=cutoff, col = "gray", lty = "dotted")
	}

	# plot.
	switch(label,
		type = {
			for(level in names(TYPE_COLOR)) {
				sel <- type == level
				points(A[sel], M[sel], col = TYPE_COLOR[level],
					pch = TYPE_PCH[[level]], bg = TYPE_BG[level])
			}
			legend.text <- names(TYPE_COLOR)
			legend.fill <- TYPE_BG
		},
		snr = {
			g <- colorRampPalette(c("red", "orange"))
			s <- c(0, 0.85, snr.cutoff)
			col <- g(length(s) - 1)

			legend.text <- c()
			legend.fill <- c()

			for(n in 1:(length(s) - 1)) {
				sel <- snr >= s[n] & snr < s[n + 1]
				points(A[sel], M[sel], col = col[n], pch = pch);
				legend.text <- c(legend.text,
					paste(format(s[n], digits = 2, nsmall = 2), "<= SNR <",
					format(s[n + 1], digits = 2, nsmall = 2)))
				legend.fill <- c(legend.fill, col[n])
			}
			sel <- snr >= snr.cutoff
			points(A[sel], M[sel], col = "black", pch = pch);
			legend.text <- c(legend.text, paste("SNR >=", snr.cutoff))
			legend.fill <- c(legend.fill, "black")
		},
		none = {
			points(A, M, pch=pch);
		}
	)

	# lowess line block.
	# remove NA.
	sel <- which(!is.na(M))
	M.l <- M[sel]
	A.l <- A[sel]
	# take sample.
	subset = sample(1:length(M.l), min(c(5000, length(M.l))))
	A.l <- A.l[subset]
	M.l <- M.l[subset]
	# order it and remove duplicates.
	o <- order(A.l[subset])
	o <- which(!duplicated(A.l))
	# draw loess line.
	lines(approx(lowess(A.l[o], M.l[o])), col = "green", lwd=4)
	
	# highligh genes.
	if(!is.null(high.list))
		points(A[high.list], M[high.list], col=high.col, pch=high.pch, bg=high.bg)

	# title.
	if(is.null(title)) title(names)
	else title(title)

	# guess legend position.
	if(is.null(legend.x)) {
		if(abs(ylim[1]) > abs(ylim[2]))
			legend.x <- "bottomright"
		else
			legend.x <- "topright"
	}
		
	if(label != "none") legend(x=legend.x, legend=legend.text, fill=legend.fill, inset=0.05)
}

# basic plotma function.
plotma <- function(A, M, label = "type", cutoff = c(-1, 1), 
	snr.cutoff = 1, legend.x, pch = ".", xlim, ylim, type, snr, ...) {
	
	label <- match.arg(label, c("type", "snr", "none"))
	if(missing(type) || missing(snr) && label != "none") {
		label = "none"
		warning("missing type or snr info.")
	}
	
	# basic plot.
	if(missing(xlim)) xlim <- range(A, na.rm = TRUE)
	if(missing(ylim)) ylim <- range(M, na.rm = TRUE)

	plot(0, col = "white", xlim = xlim, ylim = ylim, xlab = "A", ylab = "M")
	abline(h = 0, col = "steelblue", lwd = 2)
	if(!is.null(cutoff))
		abline(h = cutoff, col = "gray", lty = "dotted")

	# plot MA values.
	switch(label,
		type = {
			#type <- probeTypes(x)
			for(level in names(TYPE_COLOR)) {
				sel <- type == level
				points(A[sel], M[sel], col = TYPE_COLOR[level],
					pch = TYPE_PCH[[level]], bg = TYPE_BG[level], ...)
			}
			legend.text = names(TYPE_COLOR)
			legend.fill = TYPE_BG
		},
		snr = {
			#snr <- meanSNR(x)
			s <- c(0, 0.85, snr.cutoff)
			col <- c("red", "orange")

			legend.text <- c()
			legend.fill <- c()

			for(n in 1:(length(s) - 1)) {
				sel <- snr >= s[n] & snr < s[n + 1]
				points(A[sel], M[sel], col = col[n], pch = pch, ...)
				legend.text <- c(legend.text,
					paste(format(s[n], digits = 2, nsmall = 2), "<= SNR <",
					format(s[n + 1], digits = 2, nsmall = 2)))
				legend.fill <- c(legend.fill, col[n])
			}

			sel <- snr >= snr.cutoff
			points(A[sel], M[sel], col = "black", pch = pch, ...)
			legend.text <- c(legend.text, paste("SNR >=", snr.cutoff))
			legend.fill <- c(legend.fill, "black")
		},
		none = points(A, M, pch = pch, ...)
	)

	# guess legend position.
	if(missing(legend.x)) {
		if(abs(ylim[1]) > abs(ylim[2]))
			legend.x <- "bottomright"
		else
			legend.x <- "topright"
	}

	if(label != "none")
		legend(x = legend.x, legend = legend.text, fill = legend.fill,
			bty = "n")
}

### plotMA()
## MA plot of gene intensities.
#plotMA <- function(object, array1=1, array2=2, cutoff=NULL, label="type",
	#type=NULL, high.list=NULL, high.col="gray", high.pch=21, high.bg="orange",
	#snr.cutoff=1, legend.x="bottomright", pch=".", subset=NULL, title=NULL, 
	#xlim=NULL, ylim=NULL)
#{
	##if(!is(object,"Codelink")) stop("Codelink object needed.")
##	if(!is.null(high.list) && (!is(high.list,"logical") || !is(high.list,"vector"))) stop("logical vector needed.")
##	if(!is.null(high.list) && length(high.list) != dim(object)[1]) stop("high.list and number of genes differ.")

	#switch(class(object),
		#Codelink={
			#if(!is.null(object$Smean)) {
					#val1 <- object$Smean[,array1]
					#val2 <- object$Smean[,array2]
					#what <- "Smean"
			#}
			#if(!is.null(object$Ri)) {
					#val1 <- object$Ri[,array1]
					#val2 <- object$Ri[,array2]
					#what <- "Ri"
			#}
			#if(!is.null(object$Ni)) {    
					#val1 <- object$Ni[,array1]
					#val2 <- object$Ni[,array2]
					#what <- "Ni"
			#}
			#if(!object$method$log) {
				#val1 <- log2(val1)
				#val2 <- log2(val2)
			#}
			## M, A computation.
			#M <- val2 - val1
			#A <- (val2 + val1)/2
			#type <- object$type
			#snr <- object$snr
		#},
		#MArrayLM={
			#M <- as.matrix(object$coefficients)[, array1]
			#A <- object$Amean
			#what <- "MArrayLM"
		#},
		#stop("invalid parameter object")
	#)
	
	#label <- match.arg(label,c("type", "snr", "none"))
	#if(is.null(type) & label!="none") {
		#cat("no type information, reverting to label 'none'\n")
		#label <- "none"
	#}

	#if(!is.null(subset)) {
		#if(is.null(type)) stop("need type information for subset.\n")
		#else {
			#subset <- match.arg(subset, levels(as.factor(type)), several.ok=TRUE)
			#subset.sel <- type %in% as.character(subset)
			#M <- M[subset.sel]
			#A <- A[subset.sel]
			#type <- type[subset.sel]
			#if(class(object)=="Codelink") snr <- snr[subset.sel,]
			#if(!is.null(high.list)) high.list <- high.list[subset.sel]
		#}
	#}


        ## Range computation.
	#if(is.null(xlim)) xlim <- range(A, na.rm=TRUE)
	#if(is.null(ylim)) ylim <- range(M, na.rm=TRUE)
	## Plotting.
	#switch(label,
                #type = {
                        #negative <- type=="NEGATIVE"
                        #positive <- type=="POSITIVE"
                        #discovery <- type=="DISCOVERY"
                        #fiducial <- type=="FIDUCIAL"
                        #other <- type=="OTHER"
                        #plot(A[discovery], M[discovery], xlim=xlim, ylim=ylim, xlab="A", ylab="M", pch=pch)
                        #points(A[negative], M[negative], col="red", pch=20)
                        #points(A[positive], M[positive], col="blue", pch=20)
                        #points(A[fiducial], M[fiducial], col="yellow", pch=20)
                        #points(A[other], M[other], col="green", pch=20)
			#legend.text <- c("DISCOVERY", "NEGATIVE", "POSITIVE", "FIDUCIAL", "OTHER")
			#legend.fill <- c("black","red","blue","yellow","green")
		#},
		#snr = {
			#sel.1 <- snr[, array1] >= snr.cutoff
			#sel.2 <- snr[, array2] >= snr.cutoff
			#plot(A[sel.1 & sel.2], M[sel.1 & sel.2], xlim=xlim, ylim=ylim, col="black", xlab="A", ylab="M", pch=pch)
                        #points(A[xor(sel.1, sel.2)], M[xor(sel.1, sel.2)], col="orange", pch=".")
                        #points(A[!sel.1 & !sel.2], M[!sel.1 & !sel.2], col="red", pch=".")
			#legend.text <- c("SNR >= 1 in all","SNR >= 1 in any","SNR < 1 in all")
			#legend.fill <- c("black", "orange", "red")
		#},
		#none = {
			#plot(A, M, xlab="A", ylab="M", pch=pch);
		#}
	#)

	## Misc.
    #abline(h=0, col="steelblue")
	#if(!is.null(cutoff)) {
	        #abline(h=-cutoff, col = "gray", lty="dotted")
        	#abline(h=cutoff, col = "gray", lty="dotted")
	#}
	## Highlighted genes.
	#if(!is.null(high.list)) {
		##names <- object$name[high.list]
		##text(A[high.list], M[high.list], names, col="blue", cex=0.75)
		#points(A[high.list], M[high.list], col=high.col, pch=high.pch,
			#bg=high.bg)
	#}

	### Lowess line.
	## Remove NA.
	#sel <- which(!is.na(M))
	#M <- M[sel]
	#A <- A[sel]
	## Take a sample.
	#subset=sample(1:length(M),min(c(10000, length(M))))
	#A <- A[subset]
	#M <- M[subset]
	## Order it and remove duplicates.
	#o <- order(A[subset])
	#o <- which(!duplicated(A))
	## draw the line.
	#lines(approx(lowess(A[o], M[o])), col = "green", lwd=4)
	
	#switch(class(object),
		#Codelink={
			#names <- paste(object$sample[array2],"-",object$sample[array1], sep="")
		#},
		#MArrayLM={
			#names <- colnames(object$contrasts)[array1]
		#}
	#)
	##if(is.null(title)) title(paste(names, " MA Plot of ", what, sep=""))
	#if(is.null(title)) title(names)
	#else title(title)
	#if(label != "none") legend(x=legend.x, legend=legend.text, fill=legend.fill, inset=0.05)
#}

## plotDensities()
# Densities plot of gene intensities.
plotDensities <- function(object, subset=1:dim(object)[2], title=NULL, legend.cex=1, what=NULL) {
        if(!is(object,"Codelink")) stop("Codelink object needed.")
	if(is.null(what)) {
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
	} else {
		what <- match.arg(what, c("bg","snr","smean","ri","ni"))
		switch(what,
			snr={
				val <- object$snr
				what <- "SNR"
			},
			bg={
				val <- object$Bmedian
				what <- "Bmedian"
			},
			smean={
				val <- object$Smean
				what <- "Smean"
			},
			ri={
				val <- object$Ri
				what <- "Ri"
			},
			ni={
				val <- object$Ni
				what <- "Ni"
			}
		)
	}
	if(what!="Ni") val <- log2(val)
	if(what=="Ni" & !object$method$log) val <- log2(val)
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

## imageCodelink()
# Function to plot images of arrays.
imageCodelink <- function (object, array = 1, what = "bg",
    low="black", high="white", mar = c(1, 1, 1, 1),
    gr=1, gc=1, log.it=FALSE, ...) {

	what <- match.arg(what, c("bg", "snr", "smean", "ri", "ni"))
	switch(what,
		bg = val <- object$Bmedian[, array],
		snr = val <- object$snr[, array],
		smean = val <- object$Smean[, array],
		ri = val <- object$Ri[,	array],
		ni = val <- object$Ni[, array]
	)
	if(!object$method$log & log.it)
		val <- log2(val)

	o.row <- range(object$logical[, "row"])
	o.col <- range(object$logical[, "col"])
	
	foo <- matrix(NA, nrow = o.row[2], ncol = o.col[2])

	for (n in 1:dim(object)[1]) {
		foo[object$logical[n, "row"],object$logical[n, "col"]] <- val[n]
	}

	col <- colorRampPalette(c(low, high))(123)
	old.par <- par(mar = mar)
	on.exit(par(old.par))

	sc <- o.col[2]/gc
	sr <- o.row[2]/gr

	image(0:(gr * sr), 0:(gc * sc), foo, col = col, xaxt="n", yaxt="n", ...)

	for (igrid in 0:gc) lines(c(0, gr * sr), rep(igrid * sc, 2))
		for (igrid in 0:gr) lines(rep(igrid * sr, 2), c(0, gc * sc))
	mtext(paste("Slide:", array," File:", object$file[array]), side=1, cex=0.8)
}
# arrayNew
# creates a suitable x11 device to see the chip with the correct dimensions.
arrayNew <- function(f=2, chip="rwgcod") {
	chip <- match.arg(chip, c("rwgcod", "mwgcod", "hwgcod", "h20kcod"))
	# This is hardcoded as it is not way to guess.
	switch(chip,
		rwgcod={
			gc <- 1
			gr <- 8
			sc <- 112
			sr <- 41
		},
		mwgcod={
			gc <- 1
			gr <- 10
			sc <- 112
			sr <- 41
		},
		hwgcod={
			gc <- 1
			gr <- 12
			sc <- 112
			sr <- 42
		},
		h20kcod={
			gc <- 1
			gr <- 1
			sc <- 71
			sr <- 332
		}
	)
	r <- sr / sc
	x11(width=f*gr*r, height=f*gc)
}
