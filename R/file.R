## readHeader()
# Read header information from codelink file.
readHeader <- function(file, dec=FALSE) {
	# Size of header:
	#nlines <- 0
	#repeat {
        #        foo <- scan(file, skip=nlines, nlines=1, sep="\t", what="", quiet=TRUE)
        #        nlines<-nlines+1
        #        #if(substr(foo[1],1,2)=="--") break # detect end of header.
##		if(any(grep("-{30}", foo))) {
#			break
#		} else {
#			if(nlines > 30) stop("30 lines without finding dashed header delimiter. Is this a Codelink exported file?")
#		}
#        }
	foo <- grep("-{80}", scan(file, nlines=30, flush=T, quiet=T, what=""))
	if(!any(foo)) stop("Not a Codelink exported file.")
	nlines <- foo
	# Return list:
	head <- list(header=NULL, nlines=NULL, product=NULL, sample=NULL, size=NULL, dec=NULL, columns=NULL)
        head$header <- scan(file, nlines=nlines, sep="\t", what="", quiet=TRUE)
	head$nlines <- nlines
	if(any(foo <- grep("PRODUCT", head$header))) head$product <- head$header[foo+1] else head$product <- "Unknown"
	if(any(foo <- grep("Sample Name", head$header))) head$sample <- head$header[foo+1] else head$sample <- file
	head$size <- arraySize(file, nlines)
	if(dec) head$dec <- decDetect(file, nlines)
	head$columns <- scan(file, skip=nlines, nlines=1, sep="\t", what="", quiet=TRUE)
	return(head)
}
# Read header information from codelink file.
# Read for text files exported from Excel.
#readHeaderExcel <- function(file, dec=FALSE) {
#	nlines <- 0
#	repeat {
#                foo <- scan(file, skip=nlines, nlines=1, sep="\t", what="", quiet=TRUE)
#                nlines<-nlines+1
#	}
#}
## decDetect()
# detect decimal point.
decDetect <- function(file, nlines) {
	foo <- read.table(file, skip=nlines, nrows=1, header=TRUE, sep="\t", na.strings="")
        val <- NULL
        if(!is.null(foo$Spot_mean)) val <- foo$Spot_mean
        if(!is.null(foo$Raw_intensity) && is.null(val)) val <- foo$Raw_intensity
        if(!is.null(foo$Normalized_intensity) && is.null(val)) val <- foo$NormalNormalized_intensity
        if(is(val,"numeric")) dec <- "." else dec <- ","
	return(dec)
}
## arraySize()
# calculate chip size.
arraySize <- function(file, nlines) {
        data <- scan(file, skip=nlines+1, sep="\t", what="integer", flush=TRUE, na.strings="", quiet=TRUE)
        #ngenes <- length(data)  # number of genes.
        ngenes <- length(data[!is.na(data)]) # Codelink exporter is a little buggy. 
}
## readCodelink()
# Dynamic detection of gene number.
readCodelink <- function(files=list.files(pattern="TXT"), sample.name=NULL, flag=list(M=NA,I=NA,C=NA), dec=NULL, type="Spot", preserve=FALSE, verbose=2) {
	if(length(files)==0) stop("Codelink files not found.")
	type <- match.arg(type,c("Spot", "Raw", "Norm"))
	nslides <- length(files)
	if(!is.null(sample.name) && (length(sample.name) != nslides)) stop("sample.name must have equal length as chips loaded.")
	
	# Read Header.
	#head <- readHeader(files[1])
	# Read arrays.
	for(n in 1:nslides) {
                if(verbose>0 && n>1) cat(paste("* Reading file", n, "of", nslides, ":", files[n], "\n"))
		if(is.null(dec)) head <- readHeader(files[n], dec=TRUE) else head <- readHeader(files[n])
		if(n==1) {
			product <- head$product
			ngenes <- head$size
			if(verbose>0) cat("* Product:", product, "\n")
			if(verbose>0) cat("* Chip size:", ngenes, "\n")
                	if(verbose>0) cat(paste("* Reading file", n, "of", nslides, ":", files[n], "\n"))

			# Define Codelink object.
			Y <- matrix(NA, nrow=ngenes, ncol=nslides, dimnames=list(1:ngenes,1:nslides))
			Z <- rep(NA, ngenes)
			X <- c(1:nslides)
			J <- matrix(NA, nrow=ngenes, ncol=2, dimnames=list(1:ngenes, c("row","col")))
			method.list <- list(background="NONE", normalization="NONE", merge="NONE", log=FALSE)
			head.list <- list(product="NONE", sample=X,  file=X,
					name=Z, type=Z, flag=Y, method=method.list, Bstdev=Y, snr=Y, logical=J)
			switch(type,
				Spot = data.list <- list(Smean=Y, Bmedian=Y),
				Raw = data.list <- list(Ri=Y),
				Norm = data.list <- list(Ni=Y)
			)
			codelink <- c(head.list,  data.list)
		}

		if(verbose>2) print(head)
        	if(verbose>1) cat(paste("  + Detected '", head$dec, "' as decimal symbol.\n",sep=""))

		if(head$product != product) stop("Different array type (", head$product, ")!")
		if(head$product == "Unknown") warning("Product type for ", files[n], " is unknown (missing PRODUCT field in header).")
		if(head$size != ngenes) stop("Mmm. Something is wrong. Different number of probes (", head$size, ")\n")

		if(is.null(sample.name)) codelink$sample[n] <- head$sample
		if(verbose>1) cat(paste("  + Sample Name:", codelink$sample[n],"\n"))

		# Read bulk data.
                data <- read.table(files[n], skip=head$nlines, sep="\t", header=TRUE, row.names=1, nrows=head$size, quote="", dec=head$dec)

		# Assign Flag values.
                codelink$flag[,n] <- as.character(data[,"Quality_flag"])
		# Flag information.
		#flag.m <- codelink$flag[,n]=="M"	# MSR masked spots.
		#flag.i <- codelink$flag[,n]=="I"	# Irregular spots.
		#flag.c <- codelink$flag[,n]=="C"	# Background contaminated
		#flag.s <- codelink$flag[,n]=="S"	# Saturated spots.
		#flag.g <- codelink$flag[,n]=="G"	# Good spots.
		#flag.l <- codelink$flag[,n]=="L"	# Limit spots.
		#flag.x <- codelink$flag[,n]=="X"	# User excluded spots.
		# Allow combination of different flags.
		flag.m <- grep("M", codelink$flag[,n])	# MSR masked spots.
		flag.i <- grep("I", codelink$flag[,n])	# Irregular spots.
		flag.c <- grep("C", codelink$flag[,n])	# Background contaminated
		flag.s <- grep("S", codelink$flag[,n])	# Saturated spots.
		flag.g <- grep("G", codelink$flag[,n])	# Good spots.
		flag.l <- grep("L", codelink$flag[,n])	# Limit spots.
		flag.x <- grep("X", codelink$flag[,n])	# User excluded spots.

		if(verbose>1) cat("  + Quality flags:\n")
		#if(verbose>1) cat(paste("      G:",length(which(flag.g)),"\t"))
		#if(verbose>1) cat(paste("      L:",length(which(flag.l)),"\t"))
		#if(verbose>1) cat(paste("      M:",length(which(flag.m)),"\t"))
		#if(verbose>1) cat(paste("      I:",length(which(flag.i)),"\n"))
		#if(verbose>1) cat(paste("      C:",length(which(flag.c)),"\t"))
		#if(verbose>1) cat(paste("      S:",length(which(flag.s)),"\t"))
		#if(verbose>1) cat(paste("      X:",length(which(flag.x)),"\n"))
		if(verbose>1) cat(paste("      G:",length(flag.g),"\t"))
		if(verbose>1) cat(paste("      L:",length(flag.l),"\t"))
		if(verbose>1) cat(paste("      M:",length(flag.m),"\t"))
		if(verbose>1) cat(paste("      I:",length(flag.i),"\n"))
		if(verbose>1) cat(paste("      C:",length(flag.c),"\t"))
		if(verbose>1) cat(paste("      S:",length(flag.s),"\t"))
		if(verbose>1) cat(paste("      X:",length(flag.x),"\n"))

		# Assignn Intensity values.
		switch(type,
			Spot = {
                		codelink$Smean[,n] <- data[,"Spot_mean"]
                		codelink$Bmedian[,n] <- data[,"Bkgd_median"]
				# If found Background sd get it to compute SNR.
				if(any(grep("Bkgd_stdev", head$columns))) codelink$Bstdev[,n] <- data[,"Bkgd_stdev"] else codelink$Bstdev[,n] <- NA
				# Set values based on Flags.
				if(!is.null(flag$M)) {
					#sel.flag <- codelink$Smean[flag.m, n] > flag$M
					#codelink$Smean[sel.flag, n] <- flag$M
					#codelink$Bmedian[sel.flag, n] <- flag$M
					codelink$Smean[flag.m, n] <- flag$M
					codelink$Bmedian[flag.m, n] <- flag$M
				}
				if(!is.null(flag$I)) {
					#sel.flag <- codelink$Smean[flag.i, n] > flag$I
					#codelink$Smean[sel.flag, n] <- flag$I
					#codelink$Bmedian[sel.flag, n] <- flag$I
					codelink$Smean[flag.i, n] <- flag$I
					codelink$Bmedian[flag.i, n] <- flag$I
				}
				if(!is.null(flag$C)) {
					#sel.flag <- codelink$Smean[flag.c, n] > flag$C
					#codelink$Smean[sel.flag, n] <- flag$C
					#codelink$Bmedian[sel.flag, n] <- flag$C
					codelink$Smean[flag.c, n] <- flag$C
					codelink$Bmedian[flag.c, n] <- flag$C
				}
				if(!is.null(flag$S)) {
      					#sel.flag <- codelink$Smean[flag.s, n] > flag$S
					#codelink$Smean[sel.flag, n] <- flag$S
					#codelink$Bmedian[sel.flag, n] <- flag$S
					codelink$Smean[flag.s, n] <- flag$S
                                        codelink$Bmedian[flag.s, n] <- flag$S
                                }
				if(!is.null(flag$G)) {
					#sel.flag <- codelink$Smean[flag.g, n] > flag$G
                                        #codelink$Smean[sel.flag, n] <- flag$G
                                        #codelink$Bmedian[sel.flag, n] <- flag$G
                                        codelink$Smean[flag.g, n] <- flag$G
                                        codelink$Bmedian[flag.g, n] <- flag$G
                                }
				if(!is.null(flag$L)) {
					#sel.flag <- codelink$Smean[flag.l, n] > flag$L
                                        #codelink$Smean[sel.flag, n] <- flag$L
                                        #codelink$Bmedian[sel.flag, n] <- flag$L
                                        codelink$Smean[flag.l, n] <- flag$L
                                        codelink$Bmedian[flag.l, n] <- flag$L
                                }
				if(!is.null(flag$X)) {
					#sel.flag <- codelink$Smean[flag.x, n] > flag$X
                                        #codelink$Smean[sel.flag, n] <- flag$X
                                        #codelink$Bmedian[sel.flag, n] <- flag$X
					codelink$Smean[flag.x, n] <- flag$X
					codelink$Bmedian[flag.x, n] <- flag$X
				}
			},
			Raw = {
                		codelink$Ri[,n] <- data[,"Raw_intensity"]
				codelink$method$backgrund <- "Codelink Subtract"
				# Set values based on Flags.
                                if(!is.null(flag$M)) {
					#sel.flag <- codelink$Ri[flag.m, n] > flag$M
					#codelink$Ri[sel.flag, n] <- flag$M
					codelink$Ri[flag.m, n] <- flag$M
				}
                                if(!is.null(flag$I)) {
					#sel.flag <- codelink$Ri[flag.i, n] > flag$I
                                        #codelink$Ri[sel.flag, n] <- flag$I
					codelink$Ri[flag.i, n] <- flag$I
				}
                                if(!is.null(flag$C)) {
					#sel.flag <- codelink$Ri[flag.c, n] > flag$C
                                        #codelink$Ri[sel.flag, n] <- flag$C
					codelink$Ri[flag.c, n] <- flag$C
				}
                                if(!is.null(flag$S)) {
                                        #sel.flag <- codelink$Ri[flag.s, n] > flag$S
                                        #codelink$Ri[sel.flag, n] <- flag$S
					codelink$Ri[flag.s, n] <- flag$S
				}
                                if(!is.null(flag$G)) {
                                        #sel.flag <- codelink$Ri[flag.g, n] > flag$G
                                        #codelink$Ri[sel.flag, n] <- flag$G
					codelink$Ri[flag.g, n] <- flag$G
				}
                                if(!is.null(flag$L)) {
                                        #sel.flag <- codelink$Ri[flag.l, n] > flag$L
                                        #codelink$Ri[sel.flag, n] <- flag$L
					codelink$Ri[flag.l, n] <- flag$L
				}
                                if(!is.null(flag$X)) {
                                        #sel.flag <- codelink$Ri[flag.x, n] > flag$X
                                        #codelink$Ri[sel.flag, n] <- flag$X
					codelink$Ri[flag.x, n] <- flag$X
				}
			},
			Norm = {
                		codelink$Ni[,n] <- data[,"Normalized_intensity"]
				codelink$method$background <- "Codelink Subtract"
				codelink$method$normalization <- "Codelink Median"
				# Set values based on Flags.
                                 if(!is.null(flag$M)) {
					#sel.flag <- codelink$Ni[flag.m, n] > flag$M
					#codelink$Ni[sel.flag, n] <- flag$M
					codelink$Ni[flag.m, n] <- flag$M
				}
                                if(!is.null(flag$I)) {
					#sel.flag <- codelink$Ni[flag.i, n] > flag$I
                                        #codelink$Ni[sel.flag, n] <- flag$I
					codelink$Ni[flag.i, n] <- flag$I
				}
                                if(!is.null(flag$C)) {
					#sel.flag <- codelink$Ni[flag.c, n] > flag$C
                                        #codelink$Ni[sel.flag, n] <- flag$C
					codelink$Ni[flag.c, n] <- flag$C
				}
                                if(!is.null(flag$S)) {
                                        #sel.flag <- codelink$Ni[flag.s, n] > flag$S
                                        #codelink$Ni[sel.flag, n] <- flag$S
					codelink$Ni[flag.s, n] <- flag$S
				}
                                if(!is.null(flag$G)) {
                                        #sel.flag <- codelink$Ni[flag.g, n] > flag$G
                                        #codelink$Ni[sel.flag, n] <- flag$G
					codelink$Ni[flag.g, n] <- flag$G
				}
                                if(!is.null(flag$L)) {
                                        #sel.flag <- codelink$Ni[flag.l, n] > flag$L
                                        #codelink$Ni[sel.flag, n] <- flag$L
					codelink$Ni[flag.l, n] <- flag$L
				}
                                if(!is.null(flag$X)) {
                                        #sel.flag <- codelink$Ni[flag.x, n] > flag$X
                                        #codelink$Ni[sel.flag, n] <- flag$X
					codelink$Ni[flag.x, n] <- flag$X
				}
 			}
		)
		if(n==1) {
                        codelink$name <- as.character(data[,"Probe_name"])
                        codelink$type <- as.character(data[,"Probe_type"])
			codelink$logical[,"row"] <- data[,"Logical_row"]
			codelink$logical[,"col"] <- data[,"Logical_col"]
                }
	}
	# Compute SNR.
        codelink$snr <- SNR(codelink$Smean, codelink$Bmedian, codelink$Bstdev)
        if(verbose>0) cat("* Computing SNR...\n")
	if(!preserve) codelink$Bstdev <- NULL
	
	if(!is.null(sample.name)) codelink$sample <- sample.name
	codelink$file <- files
	codelink$product <- product
        
	new("Codelink", codelink)
}
## writeCodelink()
#  write Codelink object to file
writeCodelink <- function(object, file=NULL, dec=".") {
	if(!is(object, "Codelink")) stop("A Codelink object is needed.")
	if(is.null(file)) stop("A file name is needed.")

	write(paste("BkgdCorrection_method", object$method$background, sep="\t"), file=file)
	write(paste("Normalization_method", object$method$normalization, sep="\t"), file=file, append=TRUE)
	write(paste("Log_transformed", object$method$log, sep="\t"), file=file, append=TRUE)

	val <- NULL
	if(!is.null(object$Smean)) val <- object$Smean
	if(!is.null(object$Ri)) val <- object$Ri
	if(!is.null(object$Ni)) val <- object$Ni
	tmp <- cbind(object$name, object$flag, val)
	head <- c("Probe_name", paste("Quality_flag", object$sample, sep="-"), paste("Intensity", object$sample, sep="-"))
	tmp2 <- rbind(Index=head, tmp)
	write.table(tmp2,file=file,append=TRUE,quote=FALSE,sep="\t",dec=dec,col.names=FALSE)
}

## reportCodelink()
# report output to HTML.
reportCodelink <- function(object, chip, filename=NULL, title="Main title", probe.type=FALSE, other=NULL, other.ord=NULL) {
	if(!is(object,"Codelink") && !is(object,"character")) stop("Codelink object or character vector needed.")
	if(probe.type && !is(object,"Codelink")) stop("Codelink object needed putting type")
	if(is.null(filename)) stop("Filename needed.")

	if(is(object,"Codelink")) genes <- object$name else genes <- object
	if(!is.null(other)) {
		if(!is.null(other.ord)) {
			ord <- order(other[[other.ord]])

			genes <- genes[ord]
			for(n in names(other)) {
				other[[n]] <- other[[n]][ord]
			}
		}
	}
	if(probe.type) genes.type <- object$type
	genes.acc <- lookUp(genes, chip, "ACCNUM")
#	genes.acc <- strsplit(unlist(genes.acc),"|",extended=FALSE)
	genes.ll <- lookUp(genes, chip, "LOCUSID")
	genes.ug <- lookUp(genes, chip, "UNIGENE")
	genes.sym <- lookUp(genes, chip, "SYMBOL")
	genes.name <- lookUp(genes, chip, "GENENAME")
#	genes.go <- lookUp(genes, chip, "GO")

	genes.list <- list(genes, genes.acc, as.character(genes.ll), as.character(genes.ug))
	head <- c("ID", "GENBANK", "LOCUSID", "UNIGENE", "SYMBOL", "NAME")
	other.list <- list(genes.sym, genes.name)
	
	if(probe.type) {
		head <- c(head, "TYPE")
		other.list <- c(other.list,list(genes.type))
	}
	if(!is.null(other)) {
		head <- c(head, names(other))
		other.list <- c(other.list, other)
	}
	htmlpage(genelist=genes.list, filename=filename, table.head=head, othernames = other.list, title=title, repository = list("gb","gb","ll","ug"))
}
