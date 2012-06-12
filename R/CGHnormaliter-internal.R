.readCghRaw <-
function (input) {
    if (class(input) == "character") {
        cat("Reading input file...\n")
        input <- read.table(input, header=TRUE, sep="\t", quote="\"", fill=TRUE)
    } else if (class(input) != "data.frame") {
        stop("Input should be either a data.frame or a character string ",
             "containing a filename.\nPlease read the documentation file ",
             "for the required format.\n")
    }
    
    # Check if sufficient columns
    if (ncol(input) < 6) {
        stop("The input file must have at least six columns.\n",
             "Please read the documentation for the required format.")
    }
    
    # Average duplicate clones
    IDfreqs <- table(input[, 1])
    if (any(IDfreqs > 1)) {
        IDfreqs <- IDfreqs[IDfreqs > 1]
        cat("Averaging", length(IDfreqs), "duplicated clones...\n")
	if (length(IDfreqs) > 100) {
	    cat("Please be patient...\n")
	}
        for (i in 1:length(IDfreqs)) {
            index <- which(input[, 1] == names(IDfreqs[i]))
            duplicates <- as.matrix(input[index, 5:ncol(input)])
            means <- apply(duplicates, 2, mean, na.rm=TRUE)
            input[index[1], 5:ncol(input)] <- means
            for (j in length(index):2) {
                input <- input[-index[j], ]
            }
        }
    }
    
    # Convert the X and Y chromosomes into numeric values
    chromosomes <- input[, 2]
    if (any(chromosomes == 'X') || any(chromosomes == 'Y')) {
        warning("Converting chromosomes `X' and `Y' into `23' and ",
                " `24' (assuming organism=human)",  immediate.=TRUE)
        chromosomes <- replace(chromosomes, chromosomes == 'X', 23)
        chromosomes <- replace(chromosomes, chromosomes == 'Y', 24)
        input[, 2] <- as.integer(as.vector(chromosomes))
    }
    
    # Sort by chromosome number and then by start position
    input <- input[order(input[, 2], input[, 3]), ]    
    
    make_cghRaw(input)
}


.calculateMA <-
function (raw.data) {
    intensities <- copynumber(raw.data)
    samples <- sampleNames(raw.data)
    features <- featureNames(raw.data)

    # Number of intensity columns must be even (ref and test per sample)
    if (ncol(intensities) %% 2 == 1) {
        stop("Two columns (ref and test) with intensities must be provided per sample")
    }    

    # Initialize matrix for log2 ratios (M) and average spot intensities (A)
    data.M <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features, NULL))
    data.A <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features, NULL))
    colnames.M <- c()
    colnames.A <- c()

    # Calculate M and A
    for (i in seq(1, ncol(intensities), 2)) {
        test <- intensities[, i]
        ref <- intensities[, i+1]
        data.M <- cbind(data.M, log2(test / ref))
        data.A <- cbind(data.A, log2(test * ref) / 2)
        colnames.M <- c(colnames.M, paste(samples[i], "_", samples[i+1], sep=""))
        colnames.A <- c(colnames.A, paste(samples[i], "_", samples[i+1], sep=""))
    }
    colnames(data.M) <- colnames.M
    colnames(data.A) <- colnames.A

    # Construct the new cghRaw object, including the average intensities (A)
    data.M <- new("cghRaw", copynumber=data.M, featureData=featureData(raw.data))
    
    list(M=data.M, A=data.A)
}


.runCGHcall <-
function (data.seg, extra.args) {
    cat("Start data calling ..\n")
    if (!"prior" %in% names(extra.args)) {
        extra.args$prior = "all";  # default in CGHnormaliter
    }
    if (!"robustsig" %in% names(extra.args)) {
        extra.args$robustsig = "no";  # default in CGHnormaliter
    }
    invisible(capture.output(data.call <- do.call("CGHcall", c(data.seg, extra.args))))
    invisible(capture.output(data.call <- ExpandCGHcall(data.call, data.seg)))
    
    data.call
}

    
.localLowess <-
function (data.call, data.A, max.losses) {
    calls <- calls(data.call)
    data.M <- copynumber(data.call)
    
    # Variable for normalization shift per sample
    shift.mean <- c()
    
    for (i in 1:ncol(data.M)) {
        # Extract normals
        nr.losses <- sum(calls[, i] == -1, na.rm=TRUE)
        nr.total <- nrow(data.call)
        frac.losses <- nr.losses / nr.total
        if (frac.losses <= max.losses[i]) {
            index.normals <- which(calls[, i] == 0)
        } else {  # Losses are considered normals
            index.normals <- which(calls[, i] == -1)
        }
        normals.M <- data.M[index.normals, i]
        normals.A <- data.A[index.normals, i]

        # Also include min and max values of A in the LOWESS
	# regression, otherwise predictions for M can yield NA's
        min.A <- min(data.A[, i])
        max.A <- max(data.A[, i])
        min.M <- data.M[which(data.A[, i] == min.A), i][1]
        max.M <- data.M[which(data.A[, i] == max.A), i][1]
        normals.A <- c(normals.A, min.A, max.A)
        normals.M <- c(normals.M, min.M, max.M)

        # Apply LOWESS regression based on the normals only
        regression <- loess(M ~ A, data.frame(A=normals.A, M=normals.M), span=0.2)
        
        # Compute the normalization values for the whole sample
        normalization.values <- predict(regression, data.frame(A=data.A[, i]))
        shift.mean <- c(shift.mean, mean(abs(normalization.values)))
        
        # Apply the normalization
        data.M[, i] <- data.M[, i] - normalization.values
    }
    
    # Store the renormalized data in a new cghCall object
    newAssayData <- new.env(parent=data.call@assayData)
    assign("copynumber", data.M, envir=newAssayData)
    data.normalized <- data.call
    data.normalized@assayData <- newAssayData

    list(data=data.normalized, shift=shift.mean)
}


.plotMA <-
function (data.A, M.before, cghCall.after) {
     # Find unique filename to store the plot in
    count <- 0
    basename <- "MAplot"
    file <- paste(basename, ".pdf", sep="")
    while (file.exists(file)) {
        file <- paste(basename, count<-count+1, ".pdf", sep="")
    }
    cat("Writing MA-plots to file:", file, "\n")
    pdf(file)

    # Generate plots
    M.after <- copynumber(cghCall.after)
    calls <- calls(cghCall.after)
    legend <- c("gain", "normal", "loss")
    palette(c("red", "black", "green"))
    for (i in 1:ncol(M.before)) {
        par(mfrow=c(2, 1))
        par(mar=c(3.5, 3.5, 2.5, 2))
        par(mgp=c(2, 0.6, 0))
        ylim <- range(-0.5, 0.5, range(M.before[,i]), range(M.after[,i]))
		    
        # MA-plot before normalization
        plot(data.A[, i], M.before[, i], ylim=ylim, xlab="A", ylab="M", pch='.', col="black")
        title(paste(sampleNames(cghCall.after)[i], "- Before normalization"), line=0.6)
        abline(h=0, col="orange", lty="dashed")
		    
        # MA-plot after normalization
        plot(data.A[, i], M.after[, i], ylim=ylim, xlab="A", ylab="M", pch='.', col=calls[,i]+2)
        title(paste(sampleNames(cghCall.after)[i], "- After normalization"), line=0.6)
        location <- ifelse(ylim[2] + ylim[1] > -0.4, "topright", "bottomright")
        legend(location, legend=legend, col=c(3,2,1), pch=20, cex=0.8, inset=0.02) 
        abline(h=0, col="orange", lty="dashed")
    }
    dev.off()
}
