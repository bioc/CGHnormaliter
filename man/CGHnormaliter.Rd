\name{CGHnormaliter}
\alias{CGHnormaliter}

\title{
  Iterative normalization of aCGH data
}

\description{
  Implementation of an iterative algorithm for normalization of aCGH data displaying imbalanced aberrations.
}

\usage{
  CGHnormaliter(data, nchrom = 24, cellularity = 1, max.losses = 0.3, stop.threshold = 0.01, max.iterations = 5, plot.MA = TRUE)
}

\arguments{
  \item{data}{either a dataframe or character string containing a filename. See section Details for the format.}
  \item{nchrom}{number of chromosomes.}
  \item{cellularity}{vector of cellularities ranging from 0 to 1 to define the contamination of your sample with healthy cells (1 = no contamination). See details for more information.}
  \item{max.losses}{maximum fraction of losses in your samples. If this fraction is exceeded, losses are considerd normals.}
  \item{stop.threshold}{threshold value for the mean difference between the LOWESS regression curves from two consecutive iterations. The iteration is terminated if this difference is below the \code{stop.threshold} for all samples.}
  \item{max.iterations}{maximum number of iterations.}
   \item{plot.MA}{logical. If \sQuote{plot.MA = TRUE}, MA-plots before and after normalization are generated of each sample. The plots are saved into a PDF file.}
}

\details{
  The input should be either a \code{data.frame} or the file name of a tabseparated text file (text files must contain a header). The first four columns should contain the name, chromosome and the start and end position in bp for each array target respectively. The position columns must contain numbers only. Following these are two columns with the raw test and reference intensities for each of your samples. These intensities must be numeric as well. If the input type is a text file, missing values should be represented as 'NA' or an empty field. There is a \code{CGHnormaliter.write.table} method that prints the results in a tabular format.
  The cellularity and baselevel parameters should both be a vector of length \code{n} where \code{n} is the number of samples in your dataset. Each vector is recycled if there are not enough values in it, or truncated if there are too many.
}

\value{
  This function returns a matrix of objects of class \code{cghCall} with dimension (number of clones) * (number of samples). Each object contains the following components (See section Examples on how to access them):
    \item{normalized data}{A matrix with the normalized log2 intensity ratios for each profile.}
    \item{segments}{A matrix with the segments for each profile.}
    \item{calls}{A \code{data.frame} with the calls for each profile. Values are \code{-1} (loss), \code{0} (normal) or \code{1} (gain).}
    \item{probabilities}{A \code{data.frame} with 3 columns of probe information (name, chromosome and position), followed by 3 columns with aberration probabilities for each sample.}
}

\author{
  Thomas W. Binsl, Bart P.P. van Houte, Hannes Hettling
}

\references{
  Bart P.P. van Houte, Thomas W. Binsl, Hannes Hettling, Walter Pirovano and Jaap Heringa. CGHnormaliter: an iterative strategy to enhance normalization of array CGH data with imbalanced aberrations. BMC Genomics, 10:401, 2009.
}

\examples{
  data(Leukemia)
  ## Normalize the raw intensity values of the first 3 chromosomes.
  result <- CGHnormaliter(Leukemia, nchrom=3)
  ## Get the normalized log2 intensity ratios, segments and calls
  normalized.data <- copynumber(result)
  segmented.data <- segmented(result)
  called.data <- calls(result)
  ## Plot the normalization result of sample 2
  plot(result[,2])
  ## Write the normalized log2 intensity ratios to file
  CGHnormaliter.write.table(result)
}

\keyword{methods}
\keyword{iteration}

