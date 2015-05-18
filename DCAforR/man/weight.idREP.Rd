\name{weight.idREP}
\alias{weight.idREP}
\title{Correction for weights of repeat proteins.}
\description{
  Assigns Henikoff weights to sequences on a MSA.
}
\usage{
  weight.idREP(msa.real)
}
\arguments{
  \item{msa.real}{MSA of pairs of repeats of a protein. Generally, an integer matrix.}
}
\value{
  numeric
}
\references{
Espada et al
}
\examples{
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
  w<- weight.idREP(msa)
}
\keyword{weights}\keyword{msa}

