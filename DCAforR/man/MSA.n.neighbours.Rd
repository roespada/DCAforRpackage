\name{MSA.n.neighbours}
\alias{MSA.n.neighbours}
\title{MSA N-neigbours repeats}
\description{
  Creates a MSA of repeats n-neighbouring from a MSA from PFAM.
}
\usage{
  MSA.n.neighbours(msa,n=2,outfile=paste("MSA_neighb_",n,".fasta",sep=""))
}
\arguments{
  \item{msa}{A multiple sequence alignment of proteins in numeric code, as obtined from \code{msa2num}. Generally, an integer matrix.}
  \item{n}{how distant neighbours are the repeats to paste}
  \item{outfile}{Name for the file with MSA obtained. }
}
\value{
0
}
\references{
Espada et al
}
\examples{
}
\keyword{repeat}\keyword{proteins}\keyword{neighbours}

