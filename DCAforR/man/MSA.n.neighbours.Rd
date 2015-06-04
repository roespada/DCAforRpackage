\name{MSA.n.neighbours}
\alias{MSA.n.neighbours}
\title{MSA N-neigbours repeats}
\description{
  Creates a MSA of repeats n-neighbouring from a MSA from PFAM.
}
\usage{
  MSA.n.neighbours(msa,n=2,outfile=paste("MSA_neighb_",n,".fasta",sep=""),maxdist=ncol(msa)*(n*1.5-1),mindist=ncol(msa)*(n-1))
}
\arguments{
  \item{msa}{A multiple sequence alignment of proteins in numeric code, as obtined from \code{msa2num}. Generally, an integer matrix.}
  \item{n}{how distant neighbours are the repeats to paste}
  \item{outfile}{Name for the file with MSA obtained. }
  \item{maxdist}{Maximum number of residues between consecutive repeats.}
  \item{mindist}{Minimum number of residues between consecutive repeats.}
}
\value{
0
}
\references{
Espada et al
}
\examples{
file<-system.file("PF00023_ncbi.txt", package = "DCAforR")
msa<-readAlignment(file)
MSA.n.neighbours(msa = msa,n = 1,outfile = "PF00023_2reps.txt")
#outfile will be save in working directory
}
\keyword{repeat}\keyword{proteins}\keyword{neighbours}

