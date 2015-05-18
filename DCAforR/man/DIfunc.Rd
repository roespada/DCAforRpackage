\name{DIfunc}
\alias{DIfunc}
\title{Calculates Direct Information}
\description{Calculates Direct Information from the marginal frequencies and coupling fields.
}
\usage{DIfunc(fmarg,K,naa=21,npos=nrow(fmarg))
}
\arguments{
  \item{fmarg}{Matrix with frequency normalized of residues on each position of the MSA, as the one returned by freq.marg.correction.}
  \item{K}{Matrix with coupling energies, as the one returned by coupling.fields()$eij.}
  \item{naa}{Number of aminoacids present on the MSA, default 21.}
  \item{npos}{Length of sequences at the MSA (ie, number of columns at the MSA).}
}
\value{
Direct Information Matrix.
}
\references{
Morcos, Faruck, et al. "Direct-coupling analysis of residue coevolution captures native contacts across many protein families." Proceedings of the National Academy of Sciences 108.49 (2011): E1293-E1301.

Espada et al.
}

\examples{
  ##Example of MSA
  msa.letters<-rbind(c("G","Y","V","G","S"),c("G","F","D","G","F"),c("G","Y","D","G","F"),c("G","Y","Q","G","G"))
  msa<-msa2num(msa.letters)
    
  ##Calculate marginal frequencies, sequences weighted according Henikoff:
  fm<-aa.freq.marg(msa,weight=Henikoff.w(msa)$w)
  ## Corrected marginal frequencies
  fmc<-freq.marg.correction(fm,Mef=Henikoff.w(msa)$mef)
  
  ##Calculate joint frequencies, sequences weighted according Henikoff:
  fj<-aa.freq.joint(msa,weight=Henikoff.w(msa)$w)
  ## Corrected joint frequencies
  fjc<-freq.joint.correction(fj,Mef=Henikoff.w(msa)$mef)
  
  ## Calculate connected-correlation matrix
  Cmatrix<-Cmat(fmc,fjc)

  ## Calculate coupling fields
  cf<-coupling.fields(Cmatrix)
  
  DI<-DIfunc(fmarg=fmc,K=cf$eij,naa=21)
  image(x=1:5,y=1:5,z=DI$DI,xlab="Position",ylab="Position")
}
\keyword{di}\keyword{direct}\keyword{information}
