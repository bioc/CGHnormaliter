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
    
    cghRaw(input)
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
    data.M <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features,NULL))
    data.A <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features,NULL))
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


.runCGHnormaliter <-
function (data.raw, cellularity, max.losses, stop.threshold, max.iterations, plot.MA) {
    # Perform an initial normalization, segmentation and calling
    cat("\nCGHnormaliter -- Running an initial segmentation and calling\n")
    invisible(capture.output(data.nor <-
        normalize(data.raw$M, cellularity=cellularity, smoothOutliers=FALSE)))
    data.seg <- segmentData(data.nor)
    rm(data.nor)
    data.seg <- postsegnormalize(data.seg)
    cghcall.robustsig <- "no"
    data.call <- .runCGHcall(data.seg, cghcall.robustsig)
    
    # Perform the iteration
    iteration <- 1
    repeat {
        cat("CGHnormaliter -- Iteration #", iteration, "\n")

        # Identify the normals and (re)normalize based to these normals
        data.normals <- .extractNormals(data.call, data.raw$A, max.losses)
        normalized <- .localLowess(data.call, data.raw$A, data.normals)
    
        # Print the mean normalization shift per sample
        cat("Mean normalization shift per sample:\n")
        samples <- sampleNames(data.raw$M)
        for (i in 1:length(normalized$shift)) {
            cat("  ", samples[i], ":", normalized$shift[i], "\n")
        }
        
        # Print message if an abortion criterion is reached
        convergence <- (sum(normalized$shift < stop.threshold) == length(normalized$shift))
        if (convergence) {
            cat("CGHnormaliter -- Reached convergence. ")
            cat("Running a final segmentation and calling...\n")
	    cghcall.robustsig <- "yes"
        } else if (iteration >= max.iterations) {
            cat("CGHnormaliter -- Max iterations (",max.iterations,") reached. ", sep="")
            cat("Running a final segmentation and calling...\n")
	    cghcall.robustsig <- "yes"
        }
        
        # Segment new data again and repeat the calling procedure
        data.seg <- segmentData(normalized$data)
        data.call <- .runCGHcall(data.seg, cghcall.robustsig)
	
        # If abortion criterion reached, draw MA-plots and leave iteration
        if (convergence || iteration >= max.iterations) {
	    rm(data.seg, data.normals, normalized)
            if (plot.MA) {
                .plotMA(data.raw$A, copynumber(data.raw$M), data.call)
            }
            cat("CGHnormaliter -- FINISHED\n")
            break
        }
        iteration <- iteration + 1
    }
    data.call
}


.runCGHcall <-
function (data.seg, robustsig) {
    cat("Start data calling ..\n")
    if (compareVersion(package.version("CGHcall"), "2.9.2") >= 0) {
        invisible(capture.output(data.call <- CGHcall(data.seg, robustsig=robustsig)))
        invisible(capture.output(data.call <- ExpandCGHcall(data.call, data.seg)))
    } else if (compareVersion(package.version("CGHcall"), "2.6.0") >= 0) {
        invisible(capture.output(data.call <- CGHcall(data.seg, robustsig=robustsig)))
    } else {
        invisible(capture.output(data.call <- CGHcall(data.seg)))
    }
    data.call
}

    
.extractNormals <-
function (data.call, data.A, max.losses) {
    calls <- calls(data.call)
    data.M <- copynumber(data.call)

    # A list to store the extracted copy numbers
    normals.M <- list()
    normals.A <- list()

    # Extract normals in each sample
    for (i in 1:ncol(data.M)) {
        nr.losses <- sum(calls[, i] == -1, na.rm=TRUE)
        nr.total <- nrow(data.call)
        frac.losses <- nr.losses / nr.total
        if (frac.losses <= max.losses) {
            index.normals <- which(calls[, i] == 0)
        } else {  # Losses are considered normals
            index.normals <- which(calls[, i] == -1)
        }
        normals.M[[i]] <- data.M[, i][index.normals]
        normals.A[[i]] <- data.A[, i][index.normals]
    }
    list(M=normals.M, A=normals.A)
}


.localLowess <-
function (data.call, data.A, normals) {
    data.M <- copynumber(data.call)
    
    # Variable for normalization shift per sample
    shift.mean <- c()

    # Apply LOWESS in each sample
    for (i in 1:length(normals$M)) {
        # Make sure that the min and max value for A are also included in the
        # LOWESS regression, otherwise predictions for M might lead to NA's
        min.A <- range(data.A[, i])[1]
        max.A <- range(data.A[, i])[2]
        min.M <- data.M[which(data.A == min.A)][1]
        max.M <- data.M[which(data.A == max.A)][1]
        normals$A[[i]] <- c(normals$A[[i]], min.A, max.A)
        normals$M[[i]] <- c(normals$M[[i]], min.M, max.M)

        # Apply LOWESS regression based on the normals only
        regression <- loess(M ~ A, data.frame(A=normals$A[[i]], M=normals$M[[i]]), span=0.2)
        
        # Compute the normalization values for the whole sample
        normalization.values <- predict(regression, data.frame(A=data.A[, i]))
        shift.mean <- c(shift.mean, mean(abs(normalization.values)))
        
        # Apply the normalization
        data.M[, i] <- data.M[, i] - normalization.values
    }
    
    # Create a new environment to store the renormalized data in
    newAssayData <- new.env(parent=data.call@assayData)
    assign("copynumber", data.M, envir=newAssayData)

    # Store the new environment in the return object
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
