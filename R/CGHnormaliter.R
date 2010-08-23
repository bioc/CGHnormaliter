CGHnormaliter <-
function (data, nchrom=24, stop_threshold=0.01, max_iterations=5) {
    # Read the raw intensity data and preprocess
    data.raw <- .cghRaw_read(data)
    invisible(capture.output(data.prep <- preprocess(data.raw, nchrom=nchrom)))

    # Convert into log2 ratios (M) and average intensities (A)
    data.ma <- .cghRaw_ma(data.prep)

    # Perform the iteration
    result <- .iterate_normalize_call(data.ma, stop_threshold, max_iterations)
    
    # Return a cghCall object with the normalized log2 ratios, as
    # well as the inferred segmentations and calls
    return (result)
}

CGHnormaliter.write.table <-
function (input, data.type=c("normalized","segmented","called"),
                           file=paste(data.type,".txt", sep="")) {
    # First do a type check
    if (class(input) != "cghCall") {
        stop("Input should be an object of type cghCall.")
    }
    
    # Obtain wanted data from input
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
