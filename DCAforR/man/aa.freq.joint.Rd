\name{aa.freq.joint}
\alias{aa.freq.joint}
\title{Joint frequency in weighted MSA}
\description{
Calculates joint frequencies of aminoacid in each position (column) of a multiple sequence alignment (MSA)
}
\usage{
aa.freq.joint(M,weight=rep(1,nrow(M)),naa=21,npos=dim(M)[2],colnames.aa=T)
}
\arguments{
  \item{M }{integer matrix. For protein alignments, a matrix as obtained with msa2num.}
  \item{weight}{numeric vector assigning a weight to each row (sequence). Default is all weights equal 1. }
  \item{naa}{number of aminoacids used. Default 21.}
  \item{npos}{length of sequences on the MSA, default number of columns of the MSA. }
  \item{colnames.aa}{ logical, should the results have the letter-coded aminoacids names.}
}
\value{
  Matrix with joint frequencies between pairs of aminoacids-at-positions. Columns/rows goes over positions, and each position over aminoacids. Ie, they go pos1.aa1, pos1.aa1, ... , pos1.aa21, pos2.aa1, pos2.aa2, ...  See example below.
}
\seealso{
  \code{\link{aa.freq.marg}}
}
\examples{
 ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
  
  ##Calculate marginal frequencies, all sequences weighted the same:
  fj<-aa.freq.joint(msa)
  fj[1:10,1:10]
  dim(fj) # 5x21 = 105
  
  ##Calculate marginal frequencies, sequences weighted according Henikoff:
  fj<-aa.freq.joint(msa,weight=Henikoff.w(msa)$w)
  fj[1:10,1:10]
  dim(fj) # 5x21 = 105
}
\keyword{}
