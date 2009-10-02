\name{CGHnormaliter}
\alias{CGHnormaliter}

\title{
  Iterative normalization of aCGH data
}

\description{
  Implementation of an iterative algorithm for normalization of aCGH data displaying imbalanced aberrations.
}

\usage{
  CGHnormaliter(data, nchrom = 24, stop_threshold = 0.01, max_iterations = 5)
}

\arguments{
  \item{data}{Either a dataframe or character string containing a filename. See section Details for the format.}
  \item{nchrom}{Number of chromosomes.}
  \item{stop_threshold}{Threshold value for the mean difference between the LOWESS regression curves from two consecutive iterations. The iteration is terminated if this difference is below the \code{stop_threshold} for all samples.}
  \item{max_iterations}{Maximum number of iterations.}
}

\details{
  The input should be either a \code{data.frame} or the file name of a tabseparated text file (text files must contain a header). The first four columns should contain the name, chromosome and the start and end position in bp for each array target respectively. The position columns must contain numbers only. Following these are two columns with the raw test and reference intensities for each of your samples. These intensities must be numeric as well. If the input type is a text file, missing values should be represented as 'NA' or an empty field. There is a \code{CGHnormaliter.write.table} method that prints the results in a tabular format.
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

