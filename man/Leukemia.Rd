\name{Leukemia}
\docType{data}
\alias{Leukemia}

\title{
  Array CGH experiment data on childhood acute lymphoblastic leukemia (ALL) in humans
}

\description{
  A dataframe containing 30180 rows and 10 columns, representing the array CGH data of 3 ALL samples.
}

\usage{Leukemia}

\format{
  A dataframe containing the following 10 columns:
  \describe{
    \item{CloneID}{The unique identifiers of array elements.}
    \item{Chromosome}{Chromosome number of each array element.}
    \item{Start}{Chromosomal start position in bp of each array element.}
    \item{End}{Chromosomal end position in bp of each array element.}
    \item{Case1.test}{Background corrected test intensity values for sample 1.}
    \item{Case1.ref}{Background corrected reference intensity values for sample 1.}
    \item{Case2.test}{Background corrected test intensity values for sample 2.}
    \item{Case2.ref}{Background corrected reference intensity values for sample 2.}
    \item{Case3.test}{Background corrected test intensity values for sample 3.}
    \item{Case3.ref}{Background corrected reference intensity values for sample 3.}
  }
}

\source{Provided by the authors (see references).}

\references{
  Paulsson K, Heidenblad M, Morse H, Borg A, Fioretos T, Johansson B: Identication of cryptic aberrations and characterization of translocation breakpoints using array CGH in high hyperdiploid childhood acute lymphoblastic leukemia. Leukemia 2006, 20:2002-2007.
}

\keyword{datasets}
