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
function (result, file="normalized.txt") {
    fd <- featureData(result)
    copynumbers <- copynumber(result)
    rownames(copynumbers) <- NULL
    normalized <- data.frame(probeID=featureNames(fd), Chromosome=fd$Chromosome,
                             Start=fd$Start, End=fd$End, copynumbers)
    cat("Saving normalized log2 ratios to file:", file, "\n")
    write.table(file=file, normalized, sep="\t", quote=F, row.names=F)
}
