\name{freq.marg.correction}
\alias{freq.marg.correction}
\title{Correction for marginal frequency in weighted MSA}
\description{
  Corrects marginal frequencies of aminoacids in each position (column) of a multiple sequence alignment (MSA).
}
\usage{
  freq.marg.correction(f,Mef,lambda=Mef,naa=21)
}
\arguments{
  \item{f}{Matrix with Frequency of residues on each position of the MSA, as the one returned by aa.freq.marg.}
  \item{Mef}{Number of effective sequences on the original MSA. Recommended values: sum the weights given to each sequence, or value returned by Henikoff.w$w. }
  \item{lambda}{Pseudocounter used to correct frequencies. Default is equivalent to Mef. }
  \item{naa}{number of aminoacids used. Default 21.}
}
\value{
  Matrix with positions at rows, aminoacids at columns and corrected frequencies at values. The correction is made according to \deqn{latex}{f.corrected=((\lambda/naa)+f)/(\lambda+Mef)}
}
\seealso{
  \code{\link{aa.freq.marg},\link{Henikoff.w}}
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
    
  ##Calculate marginal frequencies, sequences weighted according Henikoff:
  fm<-aa.freq.marg(msa,weight=Henikoff.w(msa)$w)
  
  ## Corrected marginal frequencies
  freq.marg.correction(fm,Mef=Henikoff.w(msa)$mef)
}
\keyword{frequency}
