.cghRaw_read <-
function (input) {
    if (class(input) == "character") {
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
        cat("\nAveraging duplicated clones...\n")
        IDfreqs <- IDfreqs[IDfreqs > 1]
        IDlabels <- names(IDfreqs)
        for (i in 1:length(IDfreqs)) {
            index <- which(input[, 1] == IDlabels[i])
            duplicates <- as.matrix(input[index, 5:ncol(input)])
            means <- apply(duplicates, 2, mean, na.rm=TRUE)
            input[index[1], 5:ncol(input)] <- means
            for (j in 2:length(index)) {
                   input <- input[-index[j],]
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
        input[,2] <- as.integer(as.vector(chromosomes))
    }

    cghRaw(input)
}


.cghRaw_ma <-
function (cghRaw_obj) {
    intensities <- copynumber(cghRaw_obj)
    samples <- sampleNames(cghRaw_obj)
    features <- featureNames(cghRaw_obj)

    # Number of intensity columns must be even (ref and test per sample)
    if (ncol(intensities) %% 2 == 1) {
        stop("Two columns (ref and test) with intensities must be provided per sample")
    }    

    # Initialize matrix for log2 ratios (M) and average spot intensities (A)
    M <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features,NULL))
    A <- matrix(0, nrow=nrow(intensities), ncol=0, dimnames=list(features,NULL))
    colnames_M <- c()
    colnames_A <- c()

    # Calculate M and A
    for (i in seq(1, ncol(intensities), 2)) {
        test <- intensities[,i]
        ref <- intensities[,i+1]
        M <- cbind(M, log2(test / ref))
        A <- cbind(A, log2(test * ref) / 2)
        colnames_M <- c(colnames_M, paste(samples[i],"_",samples[i+1], sep=""))
        colnames_A <- c(colnames_A, paste(samples[i],"_",samples[i+1], sep=""))
    }
    colnames(M) <- colnames_M
    colnames(A) <- colnames_A

    # Construct the new cghRaw object, including the average intensities (A)
    cghRaw_m <- new("cghRaw", copynumber=M, featureData=featureData(cghRaw_obj))
    
    list(M=cghRaw_m, A=A)
}

.iterate_normalize_call <-
function (data.raw, stop_threshold, max_iterations, plotMA) {
    # Perform an initial normalization, segmentation and calling
    cat("\nCGHnormaliter -- Running an initial segmentation and calling\n")
    invisible(capture.output(data.nor <- normalize(data.raw$M, smoothOutliers=FALSE)))
    data.seg <- segmentData(data.nor)
    rm(data.nor)
    cat("Start data calling ...\n")
    if (compareVersion(package.version("CGHcall"), "2.6.0") >= 0) {
        invisible(capture.output(data.call <- CGHcall(data.seg, robustsig="no")))
    } else {
        invisible(capture.output(data.call <- CGHcall(data.seg)))
    }
        
    # Perform the iteration
    iteration <- 1
    repeat {
        cat("CGHnormaliter -- Iteration #", iteration, "\n")

        # Identify the normals and (re)normalize based to these normals
        data.normals <- .extract_normals(data.call, data.raw$A)
        normalized <- .local_lowess(data.call, data.raw$A, data.normals)
        
        # Print the mean normalization shift per sample
        cat("Mean normalization shift per sample:\n")
        samples <- sampleNames(data.raw$M)
        for (i in 1:length(normalized$shift)) {
            cat("  ", samples[i], ":", normalized$shift[i], "\n")
        }
        
        # Print message if a abortion criterion is reached
        convergence <- (sum(normalized$shift < stop_threshold) == length(normalized$shift))
        if (convergence) {
            cat("CGHnormaliter -- Reached convergence. ")
            cat("Running a final segmentation and calling...\n")
        }
        else if (iteration >= max_iterations) {
            cat("CGHnormaliter -- Max iterations (",max_iterations,") reached. ", sep="")
            cat("Running a final segmentation and calling...\n")
        }
        
        # Segment new data again and repeat the calling procedure
        data.seg <- segmentData(normalized$data)
        cat("Start data calling ...\n")
        if (compareVersion(package.version("CGHcall"), "2.6.0") >= 0) {
            invisible(capture.output(data.call <- CGHcall(data.seg, robustsig="no")))
        } else {
            invisible(capture.output(data.call <- CGHcall(data.seg)))
        }	
	
        # If abortion criterion reached, draw MA-plots and leave iteration
        if (convergence || iteration >= max_iterations) {
	    rm(data.seg, data.normals, normalized)
            if (plotMA) {
                .ma_plots(data.raw$A, copynumber(data.raw$M), data.call)
            }
            cat("CGHnormaliter -- FINISHED\n")
            break
        }
        iteration <- iteration + 1
    }

    # Return the normalized data, including segmentation and calling results
    data.call
}

.extract_normals <-
function (cghCall_obj, avg_intensities) {
    # Extract the calls and the copynumbers from the cghCall object
    calls <- calls(cghCall_obj)
    copynumbers <- copynumber(cghCall_obj)

    # A list to store the extracted copy numbers
    M.normals <- list()
    A.normals <- list()

    # A variable to store the processed calls
    calls.processed <- list()

    # Check if 50% of the calls are normals
    for (i in 1:ncol(copynumbers)) {
        # Get the ith call
        ith_call <- calls[,i]

        # Get the number of losses, normals and gains and calculate
        # their fractions
        nr.losses <- length(which(ith_call == -1))
        nr.normals <- length(which(ith_call == 0))
        nr.gains <- length(which(ith_call == 1))
        nr.total <- nrow(copynumbers)
        
        frac.losses <- nr.losses / nr.total
        frac.normals <- nr.normals / nr.total
        frac.gains <- nr.gains / nr.total

        # Check if the normals are the majority
        if (frac.normals < frac.losses || frac.normals < frac.gains) {
            if (frac.losses > frac.gains) {
                ith_call <- ith_call + 1  # Make the losses normals
            }
            else {
                ith_call <- ith_call - 1  # Make the gains normals
            }
        }

        # Set the new call
        calls.processed[[i]] <- ith_call
    }

    # For every sample extract the normals
    for (i in 1:ncol(copynumbers)) {
        ith_call <- calls.processed[[i]]
        M.normals[[i]] <- copynumbers[,i][which(ith_call == 0)]
        A.normals[[i]] <- avg_intensities[,i][which(ith_call == 0)]
    }

    # Return the extracted copynumbers that are identified as normals
    list(M=M.normals, A=A.normals)
}

.local_lowess <-
function (cghCall_obj, avg_intensities, normals) {
    M.all <- copynumber(cghCall_obj)
    A.all <- avg_intensities
    M.normals <- normals$M
    A.normals <- normals$A
    
    # Variable for normalization shift per sample
    shift.mean <- c()

    # Normalize each sample
    for (i in 1:length(M.normals)) {
        # Make sure that the min and max value for A are also included in the
        # LOWESS regression, otherwise predictions for M might lead to NA's
        A.min <- range(A.all[,i])[1]
        A.max <- range(A.all[,i])[2]
        M.min <- M.all[which(A.all == A.min)][1]
        M.max <- M.all[which(A.all == A.max)][1]
        A.normals[[i]] <- c(A.normals[[i]], A.min, A.max)
        M.normals[[i]] <- c(M.normals[[i]], M.min, M.max)

        # Apply LOWESS regression based on the normals only
        regression <- loess(M ~ A, data.frame(A=A.normals[[i]], M=M.normals[[i]]), span=0.2)
        
        # Compute the normalization values for the whole sample
        normalization.values <- predict(regression, data.frame(A=A.all[,i]))
        shift.mean <- c(shift.mean, mean(abs(normalization.values)))
        
        # Apply the normalization
        M.all[,i] <- M.all[,i] - normalization.values
    }
    
    # Create a new environment to store the renormalized data in
    newAssayData <- new.env(parent=cghCall_obj@assayData)
    assign("copynumber", M.all, envir=newAssayData)

    # Store the new environment in the return object
    data.normalized <- cghCall_obj
    data.normalized@assayData <- newAssayData

    # Return the return object
    list(data=data.normalized, shift=shift.mean)
}

.ma_plots <-
function (A, M.before, cghCall.after) {
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
    legend <- c("gain","normal","loss")
    palette(c("red","black","green"))
    for (i in 1:ncol(M.before)) {
        par(mfrow=c(2,1))
        par(mar=c(3.5,3.5,2.5,2))
        par(mgp=c(2,0.6,0))
        ylim <- range(-0.5, 0.5, range(M.before[,i]), range(M.after[,i]))
		    
        # MA-plot before normalization
        plot(A[,i], M.before[,i], ylim=ylim, xlab="A", ylab="M", pch='.', col="black")
        title(paste(sampleNames(cghCall.after)[i], "- Before normalization"), line=0.6)
        abline(h=0, col="orange", lty="dashed")
		    
        # MA-plot after normalization
        plot(A[,i], M.after[,i], ylim=ylim, xlab="A", ylab="M", pch='.', col=calls[,i]+2)
        title(paste(sampleNames(cghCall.after)[i], "- After normalization"), line=0.6)
        location <- ifelse(ylim[2] + ylim[1] > -0.4, "topright", "bottomright")
        legend(location, legend=legend, col=c(3,2,1), pch=20, cex=0.8, inset=0.02) 
        abline(h=0, col="orange", lty="dashed")
    }
    dev.off()
}
