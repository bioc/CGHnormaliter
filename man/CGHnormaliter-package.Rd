\name{CGHnormaliter-package}
\alias{CGHnormaliter-package}
\docType{package}

\title{
  This package implements the CGHnormaliter algorithm which is a strategy for improved normalization of array Comparative Genomic Hybridization (aCGH) data
}

\description{
  CGHnormaliter is an iterative algorithm for normalization of aCGH data displaying imbalanced aberrations. First, provisory balanced copy numbers are identified and subsequently used for normalization based on LOWESS. These two steps are then iterated to refine the normalization. The assumption here is that the temporary exclusion of aberrations allows for a more appropriate calculation of the LOWESS regression curve. As a result, after normalization, the log2 intensity ratios of the normals will generally be closer to zero and better reflect the biological reality.
}

\details{
  \tabular{ll}{
    Package: \tab CGHnormaliter\cr
    Type: \tab Package\cr
    Version: \tab 1.3.6\cr
    Date: \tab 2010-09-08\cr
    License: \tab GPL (>= 3)\cr
  }
  The package contains two public functions. The function \code{CGHnormaliter} performs the iterative normalization of aCGH data, while the function \code{CGHnormaliter.write.table} prints normalized aCGH data to a file. See function documentation for details.
}

\author{
  Thomas W. Binsl, Bart P.P. van Houte, Hannes Hettling
}

\references{
  Bart P.P. van Houte, Thomas W. Binsl, Hannes Hettling, Walter Pirovano and Jaap Heringa. CGHnormaliter: an iterative strategy to enhance normalization of array CGH data with imbalanced aberrations. BMC Genomics, 10:401, 2009.
}

\seealso{
  \code{CGHcall},
  \code{DNAcopy}
}

\keyword{package}
