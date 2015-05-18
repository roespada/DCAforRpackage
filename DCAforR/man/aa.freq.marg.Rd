\name{aa.freq.marg}
\alias{aa.freq.marg}
\title{Marginal frequency in weighted MSA}
\description{
Calculates marginal frequencies of aminoacid in each position (column) of a multiple sequence alignment (MSA)
}
\usage{
aa.freq.marg(M,weight=rep(1,nrow(M)),npos=dim(M)[2],naa=21,colnames.aa=T)
}
\arguments{
  \item{M}{ integer matrix, of MSA on numeric-code as obtained with msa2num.}
  \item{weight}{ numeric vector assigning a weight to each row (sequence).  Default is all sequences with weight=1. }
  \item{npos}{ length of sequences on the MSA, default number of columns of the MSA. }
  \item{naa}{ number of aminoacids used. Default 21.}
  \item{colnames.aa}{ logical, should the results have the letter-coded aminoacids names.}
}
\value{
 Matrix with positions at rows, aminoacids at columns and frequencies at values.
}
\references{
}
\seealso{
  \code{\link{Henikoff.w}}
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
  
  ##Calculate marginal frequencies, all sequences weighted the same:
  aa.freq.marg(msa)
  
  ##Calculate marginal frequencies, sequences weighted according Henikoff:
  aa.freq.marg(msa,weight=Henikoff.w(msa)$w)
}
\keyword{frequency}
