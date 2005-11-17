## readHeader.Codelink()
# Read header information from codelink file.
readHeader.Codelink <- function(file, dec=FALSE) {
	# Size of header:
	nlines <- 0
	repeat {
                foo <- scan(file, skip=nlines, nlines=1, sep="\t", what="", quiet=TRUE)
                nlines<-nlines+1
                if(substr(foo[1],1,2)=="--") break # detect end of header.
        }
	# Return list:
	head <- list(header=NULL, nlines=NULL, product=NULL, sample=NULL, size=NULL, dec=NULL, columns=NULL)
        head$header <- scan(file, nlines=nlines, sep="\t", what="", quiet=TRUE)
	head$nlines <- nlines
	head$product <- head$header[grep("PRODUCT", head$header)+1]
	if(any(foo <- grep("Sample Name", head$header))) head$sample <- head$header[foo+1] else head$sample <- file
	head$size <- arraySize.Codelink(file, nlines)
	if(dec) head$dec <- dec.Codelink(file, nlines)
	head$columns <- scan(file, skip=nlines, nlines=1, sep="\t", what="", quiet=TRUE)
	return(head)
}
## dec.Codelink()
# detect decimal point.
dec.Codelink <- function(file, nlines) {
	foo <- read.table(file, skip=nlines, nrows=1, header=TRUE, sep="\t", na.strings="")
        #        val <- NULL
        if(!is.null(foo$Spot_mean)) val <- foo$Raw_intensity
        if(!is.null(foo$Raw_intensity) && is.null(val)) val <- foo$Raw_intensity
        if(!is.null(foo$Normalized_intensity) && is.null(val)) val <- foo$Raw_intensity
        if(is(val,"numeric")) dec <- "." else dec <- ","
	return(dec)
}
## arraySize.Codelink()
# calculate chip size.
arraySize.Codelink <- function(file, nlines) {
        data <- scan(file, skip=nlines+1, sep="\t", what="integer", flush=TRUE, na.strings="", quiet=TRUE)
        #ngenes <- length(data)  # number of genes.
        ngenes <- length(data[!is.na(data)]) # Codelink exporter is a little buggy. 
}
## read.Codelink()
# Dynamic detection of gene number.
read.Codelink <- function(files, sample.name=NULL, flag=list(M=NA,I=NA,C=NA,X=NA), dec=NULL, type="Spot_mean", preserve=FALSE, verbose=2) {
	type <- match.arg(type,c("Spot_mean", "Raw_intensity", "Normalized_intensity"))
	nslides <- length(files)
	if(!is.null(sample.name) && (length(sample.name) != nslides)) stop("sample.name must have equal length as chips loaded.")
	
	# Read Header.
	head <- readHeader.Codelink(files[1])
	#print(head)
	product <- head$product
	ngenes <- head$size
	if(verbose>0) cat("* Product:", product, "\n")
	if(verbose>0) cat("* Chip size:", ngenes, "\n")

	# Define Codelink object.
        Y <- matrix(NA, nrow=ngenes, ncol=nslides, dimnames=list(1:ngenes,1:nslides))
        Z <- rep(NA, ngenes)
	X <- c(1:nslides)
	head.list <- list(Product_name="NONE", Sample_name=X,  File_name=X,
		BkgdCorrection_method="NONE", Normalization_method="NONE", Merge_method="NONE",
		Log_transformed=FALSE, Probe_name=Z, Probe_type=Z, Quality_flag=Y, Bkgd_stdev=Y,
		SNR=Y)
	switch(type,
		Spot_mean = {
			data.list <- list(Spot_mean=Y, Bkgd_median=Y)
		},
		Raw_intensity = {
			data.list <- list(Raw_intensity=Y)
		},
		Normalized_intensity = {
			data.list <- list(Normalized_intensity=Y)
		}
	)
	codelink <- c(head.list,  data.list)

	# Read arrays.
	for(n in 1:nslides) {
                if(verbose>0) cat(paste("* Reading file", n, "of", nslides, ":", files[n], "\n"))
		if(is.null(dec)) head <- readHeader.Codelink(files[n], dec=TRUE) else head <- readHeader.Codelink(files[n])
		if(verbose>2) print(head)
        	if(verbose>1) cat(paste("  + Detected '", head$dec, "' as decimal symbol.\n",sep=""))

		if(head$product != product) stop("Different array type (", head$product, ")!\n")
		if(head$size != ngenes) stop("Mmm. Something is wrong. Different number of probes (", head$size, ")\n")

		if(is.null(sample.name)) codelink$Sample_name[n] <- head$sample
		if(verbose>1) cat(paste("  + Sample Name:", codelink$Sample_name[n],"\n"))

		# Read bulk data.
                data <- read.table(files[n], skip=head$nlines, sep="\t", header=TRUE, row.names=1, nrows=head$size, quote="", dec=head$dec)

		# Assign Flag values.
                codelink$Quality_flag[,n] <- as.character(data[,"Quality_flag"])
		# Flag information.
		flag.m <- codelink$Quality_flag[,n]=="M"	# MSR masked spots.
		flag.i <- codelink$Quality_flag[,n]=="I"	# Irregular spots.
		flag.c <- codelink$Quality_flag[,n]=="C"	# Background contaminated
		flag.s <- codelink$Quality_flag[,n]=="S"	# Saturated spots.
		flag.g <- codelink$Quality_flag[,n]=="G"	# Good spots.
		flag.l <- codelink$Quality_flag[,n]=="L"	# Limit spots.
		flag.x <- codelink$Quality_flag[,n]=="X"	# User excluded spots.
		if(verbose>1) cat("  + Quality flags:\n")
		if(verbose>1) cat(paste("      G:",length(which(flag.g)),"\t"))
		if(verbose>1) cat(paste("      L:",length(which(flag.l)),"\t"))
		if(verbose>1) cat(paste("      M:",length(which(flag.m)),"\t"))
		if(verbose>1) cat(paste("      I:",length(which(flag.i)),"\n"))
		if(verbose>1) cat(paste("      C:",length(which(flag.c)),"\t"))
		if(verbose>1) cat(paste("      S:",length(which(flag.s)),"\t"))
		if(verbose>1) cat(paste("      X:",length(which(flag.x)),"\n"))

		# Assignn Intensity values.
		switch(type,
			Spot_mean = {
                		codelink$Spot_mean[,n] <- data[,"Spot_mean"]
                		codelink$Bkgd_median[,n] <- data[,"Bkgd_median"]
				if(any(grep("Bkgd_stdev", head$columns))) codelink$Bkgd_stdev[,n] <- data[,"Bkgd_stdev"] else codelink$Bkgd_stdev[,n] <- NA
		#		if(!is.null(snr.list)) codelink$Bkgd_stdev[,n] <- data[,"Bkgd_stdev"]
				# Set values based on Flags.
				if(!is.null(flag$M)) {
					codelink$Spot_mean[flag.m,n] <- flag$M		# Set M spots.
					codelink$Bkgd_median[flag.m,n] <- flag$M
					#if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.m,n] <- flag$M
				}
				if(!is.null(flag$I)) {
					codelink$Spot_mean[flag.i,n] <- flag$I		# Set I spots.
					codelink$Bkgd_median[flag.i,n] <- flag$I
					#if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.i,n] <- flag$I
				}
				if(!is.null(flag$C)) {
					codelink$Spot_mean[flag.c,n] <- flag$C		# Set C spots.
					codelink$Bkgd_median[flag.c,n] <- flag$C
					#if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.c,n] <- flag$C
				}
				if(!is.null(flag$S)) {
                                        codelink$Spot_mean[flag.s, n] <- flag$S          # Set S spots.
                                        codelink$Bkgd_median[flag.s, n] <- flag$S
                                        #if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.s,n] <- flag$S
                                }
				if(!is.null(flag$G)) {
                                        codelink$Spot_mean[flag.g, n] <- flag$G          # Set G spots.
                                        codelink$Bkgd_median[flag.g,n] <- flag$G
                                        #if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.g, n] <- flag$G
                                }
				if(!is.null(flag$L)) {
                                        codelink$Spot_mean[flag.l, n] <- flag$L          # Set L spots.
                                        codelink$Bkgd_median[flag.l, n] <- flag$L
                                        #if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.l, n] <- flag$L
                                }
				if(!is.null(flag$X)) {
					codelink$Spot_mean[flag.x,n] <- flag$X		# Set X spots.
					codelink$Bkgd_median[flag.x,n] <- flag$X
					#if(!is.null(snr.list)) codelink$Bkgd_stdev[flag.x,n] <- flag$X
				}
			},
			Raw_intensity = {
                		codelink$Raw_intensity[,n] <- data[,"Raw_intensity"]
				codelink$BkgdCorrection_method <- "Codelink Subtract"
				# Set values based on Flags.
                                if(!is.null(flag$M)) {
					codelink$Raw_intensity[flag.m, n] <- flag$M	# Set M spots.
                                }
                                if(!is.null(flag$I)) {
					codelink$Raw_intensity[flag.i, n] <- flag$I	# Set I spots.
                                }
                                if(!is.null(flag$C)) {
					codelink$Raw_intensity[flag.c, n] <- flag$C	# Set C spots.
                                }
                                if(!is.null(flag$S)) {
                                        codelink$Raw_intensity[flag.s, n] <- flag$S      # Set S spots.
                                }
                                if(!is.null(flag$G)) {
                                        codelink$Raw_intensity[flag.g, n] <- flag$G      # Set G spots.
                                }
                                if(!is.null(flag$L)) {
                                        codelink$Raw_intensity[flag.l, n] <- flag$L      # Set L spots.
                                }
                                if(!is.null(flag$X)) {
					codelink$Raw_intensity[flag.x, n] <- flag$X	# Set X spots.
                                }
			},
			Normalized_intensity = {
                		codelink$Normalized_intensity[,n] <- data[,"Normalized_intensity"]
				codelink$BkgdCorrection_method <- "Codelink Subtract"
				codelink$Normalization_method <- "Codelink Median"
				# Set values based on Flags.
                                if(!is.null(flag$M)) {
					codelink$Normalized_intensity[flag.m,n] <- flag$M	# Set M spots.
                                }
                                if(!is.null(flag$I)) {
					codelink$Normalized_intensity[flag.i,n] <- flag$I	# Set I spots.
                                }
                                if(!is.null(flag$C)) {
					codelink$Normalized_intensity[flag.c,n] <- flag$C	# Set C spots.
                                }
                                if(!is.null(flag$S)) {
                                        codelink$Normalized_intensity[flag.s, n] <- flag$S       # Set S spots.
                                }
                                if(!is.null(flag$G)) {
                                        codelink$Normalized_intensity[flag.g, n] <- flag$G       # Set G spots.
                                }
                                if(!is.null(flag$L)) {
                                        codelink$Normalized_intensity[flag.l, n] <- flag$L       # Set L spots.
                                }
                                if(!is.null(flag$X)) {
					codelink$Normalized_intensity[flag.x,n] <- flag$X	# Set X spots.
                                }
			}
		)
		if(n==1) {
                        codelink$Probe_name <- as.character(data[,"Probe_name"])
                        codelink$Probe_type <- as.character(data[,"Probe_type"])
                }
	}
	# Compute SNR.
        codelink$SNR <- SNR(codelink$Spot_mean, codelink$Bkgd_median, codelink$Bkgd_stdev)
        if(verbose>0) cat("* Computing SNR...\n")
	if(!preserve) codelink$Bkgd_stdev <- NULL
	
	if(!is.null(sample.name)) codelink$Sample_name <- sample.name
	codelink$File_name <- files
	codelink$Product_name <- product
        
	new("Codelink", codelink)
}
## write.Codelink()
# export Codelink object to file
write.Codelink <- function(object, file=NULL, dec=".") {
	if(!is(object, "Codelink")) stop("A Codelink object is needed.")
	if(is.null(file)) stop("A file name is needed.")

	write(paste("BkgdCorrection_method", object$BkgdCorrection_method, sep="\t"), file=file)
	write(paste("Normalization_method", object$Normalization_method, sep="\t"), file=file, append=TRUE)
	write(paste("Log_transformed", object$Log_transformed, sep="\t"), file=file, append=TRUE)

	val <- NULL
	if(!is.null(object$Spot_mean)) val <- object$Spot_mean
	if(!is.null(object$Raw_intensity)) val <- object$Raw_intensity
	if(!is.null(object$Normalized_intensity)) val <- object$Normalized_intensity
	tmp <- cbind(object$Probe_name,object$Quality_flag,val)
	head <- c("Probe_name", paste("Quality_flag",object$Sample_name,sep="-"),paste("Intensity",object$Sample_name,sep="-"))
	tmp2 <- rbind(Index=head, tmp)
	write.table(tmp2,file=file,append=TRUE,quote=FALSE,sep="\t",dec=dec,col.names=FALSE)
}

## report.Codelink()
# report output to HTML.
report.Codelink <- function(object, chip, filename=NULL, title="Main title", probe.type=FALSE, other=NULL, other.ord=NULL) {
	if(!is(object,"Codelink") && !is(object,"character")) stop("Codelink object or character vector needed.")
	if(probe.type && !is(object,"Codelink")) stop("Codelink object needed putting Probe_type")
	if(is.null(filename)) stop("Filename needed.")

	require(annotate)
	
	if(is(object,"Codelink")) genes <- object$Probe_name else genes <- object
	if(!is.null(other)) {
		if(!is.null(other.ord)) {
			ord <- order(other[[other.ord]])

			genes <- genes[ord]
			for(n in names(other)) {
				other[[n]] <- other[[n]][ord]
			}
		}
	}
	if(probe.type) genes.type <- object$Probe_type
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
