\name{load}
\alias{msa2num}
\title{AminoAcids matrix to integer matrix}
\description{
  Changes the letter-code of a amino-acids multiple alignment to a number-coded.
}
\usage{
msa2num(msa)
}
\arguments{
  \item{msa}{A multiple sequence alignment of proteins as obtined from \code{readAlignment}.}
}
\value{
  Integer matrix with the following code: 
  1=-; 2=A ,3=C, 4=D, 5=E, 6=F, 7=G, 8=H, 9=I, 10=K,
  11=L, 12=M, 13=N, 14=P, 15=Q, 16=R, 17=S, 18=T, 19=V, 20=W, 21=Y
  
}

\examples{
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
}
\keyword{msa}