\name{splice}
\alias{splice}
\docType{data}
\title{Dataset of human donor splice sites}
\description{
Dataset of 400 human donor splice sites with a sequence length of 7 base
pairs.
}
\usage{data(splice)}
\format{
  \describe{
    \item{y}{binary response. True (1) or false (0) splice site.}
    \item{Pos.x}{DNA letter (A, C, G, T) at position x, where x ranges
      from 1 to 7.} 

  }
}
\details{The dataset is a random subset of the MEMset Donor dataset
  used in Gene et al. (2004). 
}
\references{Gene, Y. and Burge, C. (2004) \emph{Maximum Entropy Modeling
    of Short Sequence Motifs with Applications to RNA Splicing Signals},
  Journal of Computational Biology, 11, 475 - 494.
}

\examples{
  data(splice)
}
\keyword{datasets}
