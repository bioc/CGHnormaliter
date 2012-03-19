CGHnormaliter <-
function (data, nchrom = 24, cellularity = 1, max.losses = 0.3, plot.MA = TRUE, ...) {
    max.iterations <- 5
    
    # Extract proper arguments for segment (DNAcopy) and/or CGHcall
    args.extra <- list(...)
    args.segment <- args.extra[!is.na(pmatch(names(args.extra), names(formals(segment))))]
    args.CGHcall <- args.extra[!is.na(pmatch(names(args.extra), names(formals(CGHcall))))]
    args.extra[names(args.segment)] <- NULL
    args.extra[names(args.CGHcall)] <- NULL
    if (length(args.extra) > 0) {
        names.unused <- paste(unlist(names(args.extra)), collapse=", ")
        warning("unused argument(s): ", names.unused, immediate.=TRUE)
    }
    
    # Read the raw intensity data and preprocess
    data.raw <- .readCghRaw(data)
    invisible(capture.output(data.prep <- preprocess(data.raw, nchrom=nchrom)))
    rm(data.raw)
    
    # Convert intensities into log2 ratios (M) and average intensities (A)
    data.ma <- .calculateMA(data.prep)
    rm(data.prep)
    
    # Expand or reduce size of 'max.losses' to number of samples
    if (sum(max.losses < 0 || max.losses > 1) > 0) {
        warning("Some values of max.losses are outside range [0,1]", immediate.=TRUE)
    }
    if (length(max.losses) < ncol(data.ma$M)) max.losses <- rep(max.losses, ncol(data.ma$M));
    if (length(max.losses) > ncol(data.ma$M)) max.losses <- max.losses[1:ncol(data.ma$M)];    
    
    # Initial normalization, segmentation and calling
    cat("\nCGHnormaliter -- Running an initial segmentation and calling\n")
    invisible(capture.output(data.nor <-
        normalize(data.ma$M, cellularity=cellularity, smoothOutliers=FALSE)))
    data.seg <- do.call("segmentData", c(data.nor, args.segment))
    rm(data.nor)
    data.seg <- postsegnormalize(data.seg)
    data.call <- .runCGHcall(data.seg, args.CGHcall)
    
    # Perform the iteration
    for (iteration in 1:max.iterations) {
        cat("CGHnormaliter -- Iteration #", iteration, "\n")
	
	# LOWESS based on normals only
        normalized <- .localLowess(data.call, data.ma$A, max.losses)
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
        
	# Redo segmentation and calling
        data.seg <- do.call("segmentData", c(normalized$data, args.segment))
        data.call <- .runCGHcall(data.seg, args.CGHcall)
	
	if (convergence) {
	    break;
	}
    }
    
    # Draw MA-plots, if asked for
    if (plot.MA) {
        rm(data.seg, normalized)
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

### How to calculate intensities back from log2 ratios ###
# intensity.ratio = 2 ^ copynumber(result)
# int1 = sqrt(2 ^ data.prep$A) / sqrt(intensity.ratio)
# int2 = int1 * intensity.ratio
