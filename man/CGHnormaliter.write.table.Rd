\name{CGHnormaliter.write.table}
\alias{CGHnormaliter.write.table}

\title{Print normalized aCGH data to a file}

\description{
  This function stores the results of a CGHnormaliter normalization into a plain text file.
}

\usage{
  CGHnormaliter.write.table(input, data.type=c("normalized","segmented","called"), file=paste(data.type,".txt", sep=""))
}

\arguments{
  \item{input}{Result object of a CGHnormaliter normalization.}
  \item{data.type}{Type of aCGH data to store.}
  \item{file}{File name to store the data in.}
}

\details{
  Available aCGH data types are:
  \itemize{
    \item normalized: The normalized log2 intensity ratios
    \item segmented: The segmented log2 intensity ratios (after normalization)
    \item called: The calls (after normalization)
  }
  The selected data are stored into a tabseparated file with the specified file name.
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
  \dontrun{CGHnormaliter.write.table(result)}
  \dontrun{CGHnormaliter.write.table(result, type="segmented")}
  \dontrun{CGHnormaliter.write.table(result, type="called", file="mycalls.txt")}
}
  
\keyword{IO}
\keyword{file}

