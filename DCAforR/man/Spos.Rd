\name{Spos}
\alias{Spos}
\title{Entropy per position for a weighted multiple sequence alignment}
\description{
}
\usage{Spos(freq.marg,npos=nrow(frec.marg))
}
\arguments{
  \item{freq.marg}{Matrix with frequency normalized of residues on each position of the MSA, as the one returned by freq.marg.correction.}
  \item{npos}{length of sequences on the MSA, default number of rows of freq.marg of the MSA. }
}
\value{
Numeric with entropy value for each column of the matrix used for calculate freq.marg (ie for each position of the MSA). Calculation is \eqn{latex}{-\sum x log_2(x)}.
}
\seealso{
  \code{\link{aa.freq.marg}}
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
  
  ##Calculate marginal frequencies, all sequences weighted the same:
  fm<-aa.freq.marg(msa)
  ## Correct frequencies
  fmc<-freq.marg.correction(fm,Mef=nrow(msa))
  ##Entropy per position
  Spos(fmc)
}
\keyword{}
