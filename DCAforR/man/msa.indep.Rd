\name{msa.indep}
\alias{msa.indep}
\title{Generate MSA of random sequences}
\description{
Generates a MSA of random sequences keeping the marginal frequencies. This function only scrambles the residues on a position between sequences. 
}
\usage{msa.indep(M)
}
\arguments{
  \item{M}{Matrix. It is intended to be an MSA as returned by read.fasta(file)$ali or msa2num.}
}
\value{
Matrix with different sequences (rows) but keeping equal the marginal frequencies at each column.
}
\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  #random sequences
  msa.indep(msa.letters)
  
  #The same on number-code
  msa<-msa2num(msa.letters)
  msa.indep(M=msa)
}
\keyword{msa}
\keyword{scramble}
