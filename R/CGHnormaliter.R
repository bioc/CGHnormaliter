CGHnormaliter <-
function (data, nchrom = 24, cellularity = 1, max.losses = 0.3, plot.MA=TRUE, ...) {
    max.iterations <- 5

    # Read the raw intensity data and preprocess
    data.raw <- .readCghRaw(data)
    invisible(capture.output(data.prep <- preprocess(data.raw, nchrom=nchrom)))
    rm(data.raw)
    
    # Convert intensities into log2 ratios (M) and average intensities (A)
    data.ma <- .calculateMA(data.prep)
    rm(data.prep)
    
    # Initial normalization, segmentation and calling
    formals()
    cat("\nCGHnormaliter -- Running an initial segmentation and calling\n")
    invisible(capture.output(data.nor <-
        normalize(data.ma$M, cellularity=cellularity, smoothOutliers=FALSE)))
    formals(segmentData) <- c(formals(segmentData), alist(... = ))
    data.seg <- segmentData(data.nor, ...)
    rm(data.nor)
    data.seg <- postsegnormalize(data.seg)
    data.call <- .runCGHcall(data.seg, ...)
    
    # Perform the iteration
    for (iteration in 1:max.iterations) {
        cat("CGHnormaliter -- Iteration #", iteration, "\n")

        # Identify the normals and (re)normalize based on these normals
        data.normals <- .extractNormals(data.call, data.ma$A, max.losses)
        normalized <- .localLowess(data.call, data.ma$A, data.normals)
    
        # Print the mean normalization shift per sample
        cat("Mean normalization shift per sample:\n")
        samples <- sampleNames(data.ma$M)
        for (i in 1:length(normalized$shift)) {
            cat("  ", samples[i], ":", normalized$shift[i], "\n")
        }
        
        # Print message if an abortion criterion is reached
        convergence <- (sum(normalized$shift < 0.01) == length(normalized$shift))
        if (convergence) {
            cat("CGHnormaliter -- Reached convergence. ")
            cat("Running a final segmentation and calling...\n")
        } else if (iteration >= max.iterations) {
            cat("CGHnormaliter -- Max iterations (",max.iterations,") reached. ", sep="")
            cat("Running a final segmentation and calling...\n")
        }
        
        # Repeat the segmentation and calling procedure
        data.seg <- segmentData(normalized$data, ...)
        data.call <- .runCGHcall(data.seg, ...)
	
	# Stop iteration if convergence reached
	if (convergence) {
	    break;
	}
    }
    
    # Draw MA-plots, if asked for
    if (plot.MA) {
        rm(data.seg, data.normals, normalized)
        .plotMA(data.ma$A, copynumber(data.ma$M), data.call)
    }
    
    cat("CGHnormaliter -- FINISHED\n")
    data.call
}



CGHnormaliter.write.table <-
function (input, data.type = c("normalized", "segmented", "called"),
          file = paste(data.type,".txt", sep="")) {
    # First do a type check
    if (class(input) != "cghCall") {
        stop("Input should be an object of type cghCall.")
    }
    
    # Extract proper data from input
    data.type <- match.arg(data.type)
    if (data.type == "normalized") {
        data <- copynumber(input)
        data.string <- "normalized log2 ratios"
    } else if (data.type == "segmented") {
        data <- segmented(input)
        data.string <- "segmented log2 ratios"
    } else if (data.type == "called") {
        data <- calls(input)
        data.string <- "calls"
    }
    
    # Rebuild data frame and write data to file
    rownames(data) <- NULL
    fd <- featureData(input)
    data <- data.frame(probeID=featureNames(fd), Chromosome=fd$Chromosome,
                       Start=fd$Start, End=fd$End, data)
    cat("Saving", data.string, "to file:", file, "\n")
    write.table(file=file, data, sep="\t", quote=F, row.names=F)
}

# Calculate intensities back from log2 ratios
# intensity.ratio = 2 ^ copynumber(result)
# int1 = sqrt(2 ^ data.prep$A) / sqrt(intensity.ratio)
# int2 = int1 * intensity.ratio
