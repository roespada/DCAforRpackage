\name{msa.removed.gaps}
\alias{msa.removed.gaps}
\title{Remove positions with excess of gaps}
\description{
  Removes positons with excess of gaps. 
}
\usage{msa.removed.gaps(msa,cutoff.gaps=0.7,aa2rm="-")
}
\arguments{
  \item{msa}{Matrix of characters (aminoacids). It is intended to be an MSA as returned by read.fasta(file)$ali.}
  \item{cutoff.gaps}{Maximum number of sequences allowed with gap in a position.}
  \item{aa2rm}{Symbol representing the gaps. Usually "-" or 1.}
}
\value{
  Matrix equivalent to the input msa, with the same or less columns.
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("-","Y","V","G","S"),c("-","F","D","G","F"),c("-","Y","D","G","F"),c("G","Y","Q","G","G"))
  #random sequences
  msa.removed.gaps(msa.letters)
  

}
\keyword{msa}
\keyword{gap}
