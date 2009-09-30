\name{CGHnormaliter.write.table}
\alias{CGHnormaliter.write.table}

\title{Print normalized data to a file}

\description{
  This function stores the results of a CGHnormaliter normalization into a plain text file.
}

\usage{
  CGHnormaliter.write.table(result, file="normalized.txt")
}

\arguments{
  \item{result}{Result object of a CGHnormaliter normalization.}
  \item{file}{Filename to store the data in.}
}

\details{
  This function prints the normalized log2 intensity ratios to a tabseparated file with the specified file name. 
}

\author{
  Thomas W. Binsl, Bart P.P. van Houte, Hannes Hettling
}

\references{
  Bart P.P. van Houte, Thomas W. Binsl, Hannes Hettling, Walter Pirovano and Jaap Heringa. CGHnormaliter: an iterative strategy to enhance normalization of array CGH data with imbalanced aberrations. BMC Genomics, 10:401, 2009.
}

\seealso{
  \code{\link[CGHnormaliter:CGHnormaliter-package]{CGHnormaliter}}
}

\examples{
  data(Leukemia)
  ## Normalize the intensity ratios
  \dontrun{result <- CGHnormaliter(Leukemia)}
  ## Write the result to a file
  \dontrun{CGHnormaliter.write.table(result, "normalization_results.txt")}
}
  
\keyword{IO}
\keyword{file}

