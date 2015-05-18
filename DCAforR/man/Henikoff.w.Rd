\name{Henikoff.w}
\alias{Henikoff.w}
\title{Weights proteins sequences according to Henikoff method.}
\description{
  Assigns Henikoff weights to sequences on a MSA.
}
\usage{
  Henikoff.w(msa,fmarg,normalize=FALSE)
}
\arguments{
  \item{msa}{ An integer matrix representing a multiple sequence alignment of proteins in numeric code, as obtined from \code{msa2ali}.}
  \item{fmarg}{Marginal frequencies, such as returned by aa.freq.marg}
  \item{normalize}{Logical, should the weights be normalized such that  \eqn{latex}{\sum _i w_i=1} }
}
\value{
  A list with the following components:
  \item{w}{ vector of weights.}
  \item{mef}{ Effective number of sequences. If normalzed=T, calculated as Henikoff proposal 
  \eqn{latex}{exp(\sum(-w log(w)))}. If normalized=F, the sum of $w.}
}
\references{
  Henikoff, Steven, and Jorja G. Henikoff. "Amino acid substitution matrices from protein blocks." Proceedings of the National Academy of Sciences 89.22 (1992): 10915-10919.
}
\examples{
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
  
  w<-Henikoff.w(msa,fmarg=aa.freq.marg(msa))
  w$w
}
\keyword{weights}\keyword{msa}


